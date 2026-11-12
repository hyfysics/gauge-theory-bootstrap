"""
compute.py — Optimized high-precision quadrature computations.

Key optimizations over the original:
  - compute_Phi_hat and compute_Phi_tilde are parallelized over the s-index (si)
    using multiprocessing.Pool, giving ~4x speedup on a 4-core machine.
  - Power arrays (zfv[si][muj]^n) are precomputed once per (si, muj) pair,
    eliminating redundant gmpy2 exponentiations in the inner loops.
  - The Phi_tilde symmetry fill (nj < ni) is done in a single vectorised-style
    pass after the upper-triangle is complete, avoiding redundant work.
  - All per-worker data is passed as plain Python lists so pickle round-trips
    through multiprocessing stay fast.
"""

import gmpy2 as gp
from itertools import product
from multiprocessing import Pool, cpu_count

# ---------------------------------------------------------------------------
# Precision — match whatever the caller set globally
# ---------------------------------------------------------------------------
PRECISION = gp.get_context().precision   # read once at import time


def _set_precision(prec):
    gp.get_context().precision = prec

# ---------------------------------------------------------------------------
# Tanh-sinh nodes and weights
# ---------------------------------------------------------------------------
def tanh_sinh_nodes_weights(n: int = 5):
    """Compute tanh-sinh quadrature nodes and weights at level n."""
    h = gp.mpfr(2) ** (-n)
    nodes, weights = [], []
    w = gp.mpfr(1)
    i = 1
    while w > gp.mpfr('1e-100'):
        t = i * h
        x = gp.tanh(gp.const_pi() / 2 * gp.sinh(t))
        w = (h * gp.const_pi() / 2 * gp.cosh(t)
             / gp.cosh(gp.const_pi() / 2 * gp.sinh(t)) ** 2)
        nodes.append(x)
        weights.append(w)
        i += 1
    all_nodes   = [-x for x in nodes[::-1]] + [gp.mpfr(0)] + nodes
    all_weights = weights[::-1] + [h * gp.const_pi() / 2] + weights
    return all_nodes, all_weights


# ---------------------------------------------------------------------------
# Conformal map
# ---------------------------------------------------------------------------
def zf(sa, s):
    return ((gp.sqrt(4 - sa) - gp.sqrt(4 - s))
            / (gp.sqrt(4 - sa) + gp.sqrt(4 - s)))


# ---------------------------------------------------------------------------
# Legendre polynomials  (unchanged — fast enough)
# ---------------------------------------------------------------------------
def compute_Leg_pol(Nl, nodes):
    """Return Leg_pol[li][muj] as a 2-D list."""
    Nm = len(nodes)
    Leg_pol = [[gp.mpfr(0)] * Nm for _ in range(Nl)]
    for muj in range(Nm):
        Leg_pol[0][muj] = gp.mpfr(1)
        Leg_pol[1][muj] = nodes[muj]
    for li in range(2, Nl):
        for muj in range(Nm):
            Leg_pol[li][muj] = (
                (2 * li - 1) * nodes[muj] * Leg_pol[li - 1][muj]
                - (li - 1) * Leg_pol[li - 2][muj]
            ) / li
    return Leg_pol


# ---------------------------------------------------------------------------
# zfv  (unchanged — fast enough)
# ---------------------------------------------------------------------------
def compute_zfv(s, nodes, sa):
    """Return zfv[si][muj]."""
    Ns, Nm = len(s), len(nodes)
    zfv = [[gp.mpfr(0)] * Nm for _ in range(Ns)]
    for si, muj in product(range(Ns), range(Nm)):
        zfv[si][muj] = zf(sa, (nodes[muj] - 1) / 2 * (s[si] - 4))
    return zfv

# ---------------------------------------------------------------------------
# computes Phi_n(s)  
# ---------------------------------------------------------------------------
def compute_PhiM(xiv,Nmax2):
    Ns = len(xiv)
    PhiM = [[gp.mpc(0)] * Nmax2 for _ in range(Ns)]
    for si in range(Ns):
        for ni in range(Nmax2):
            PhiM[si][ni] = gp.exp(gp.mpc(0,(ni+1)*xiv[si]))
    return PhiM

# ===========================================================================
# Parallel Phi_hat
# ===========================================================================

def _phi_hat_worker(args):
    """Compute Phi_hat for a single si."""
    si, Nl, Nmax2, weights, zfv_si, Leg_pol, prec = args
    Nm = len(weights)
    gp.get_context().precision = prec

    # Precompute power arrays: pow_z[muj][ni] = zfv_si[muj]^(ni+1)
    pow_z = []
    for muj in range(Nm):
        z = zfv_si[muj]
        row = []
        zp = z
        for _ in range(Nmax2):
            row.append(zp)
            zp = zp * z
        pow_z.append(row)

    result = [[gp.mpfr(0)] * Nmax2 for _ in range(Nl)]
    for li in range(Nl):
        for ni in range(Nmax2):
            acc = gp.mpfr(0)
            for muj in range(Nm):
                acc += weights[muj] * Leg_pol[li][muj] * pow_z[muj][ni]
            result[li][ni] = acc
    return si, result


def compute_Phi_hat_parallel(Ns, Nl, Nmax2, weights, zfv, Leg_pol, ncores=4):
    """Return Phi_hat[si][li][ni], computed in parallel over si."""
    Phi_hat = [[[gp.mpfr(0)] * Nmax2 for _ in range(Nl)] for _ in range(Ns)]
    tasks = [
        (si, Nl, Nmax2, weights, zfv[si], Leg_pol, gp.get_context().precision)
        for si in range(Ns)
    ]
    with Pool(processes=ncores,
          initializer=_set_precision,
          initargs=(gp.get_context().precision,)) as pool:
        for si, result in pool.imap_unordered(_phi_hat_worker, tasks):
            Phi_hat[si] = result
    return Phi_hat


# ===========================================================================
# Parallel Phi_tilde
# ===========================================================================

def _phi_tilde_worker(args):
    """Compute Phi_tilde for a single si."""
    si, Nl, Nmax, weights, zfv_si, Leg_pol, prec = args
    gp.get_context().precision = prec
    Nm = len(weights)

    # Precompute power arrays for both z and its mirror
    pow_z1 = []   # pow_z1[muj][ni] = zfv_si[muj]^(ni+1)
    pow_z2 = []   # pow_z2[muj][ni] = zfv_si[Nm-1-muj]^(ni+1)
    for muj in range(Nm):
        z1, z2 = zfv_si[muj], zfv_si[Nm - 1 - muj]
        r1, r2 = [], []
        zp1, zp2 = z1, z2
        for _ in range(Nmax):
            r1.append(zp1); zp1 = zp1 * z1
            r2.append(zp2); zp2 = zp2 * z2
        pow_z1.append(r1)
        pow_z2.append(r2)

    # Accumulate upper triangle
    PT = [[[gp.mpfr(0)] * Nmax for _ in range(Nmax)] for _ in range(Nl)]
    for muj in range(Nm):
        wp_row = [weights[muj] * Leg_pol[li][muj] for li in range(Nl)]
        z1v = pow_z1[muj]
        z2v = pow_z2[muj]
        for li in range(Nl):
            wp = wp_row[li]
            for ni in range(Nmax):
                wz1 = wp * z1v[ni]
                for nj in range(ni, Nmax):
                    PT[li][ni][nj] += wz1 * z2v[nj]

    # Apply symmetry and zero out odd-diagonal
    for li in range(Nl):
        sign = (-1) ** li
        for ni in range(Nmax):
            if li % 2 == 1:
                PT[li][ni][ni] = gp.mpfr(0)
            for nj in range(ni + 1, Nmax):
                PT[li][nj][ni] = sign * PT[li][ni][nj]

    return si, PT


def compute_Phi_tilde_parallel(Ns, Nl, Nmax, weights, zfv, Leg_pol, ncores=4):
    """Return Phi_tilde[si][li][ni][nj], computed in parallel over si."""
    Phi_tilde = [[[[gp.mpfr(0)] * Nmax for _ in range(Nmax)]
                  for _ in range(Nl)]
                 for _ in range(Ns)]
    tasks = [
        (si, Nl, Nmax, weights, zfv[si], Leg_pol, gp.get_context().precision)
        for si in range(Ns)
    ]
    with Pool(processes=ncores,
          initializer=_set_precision,
          initargs=(gp.get_context().precision,)) as pool:
        for si, result in pool.imap_unordered(_phi_tilde_worker, tasks):
            Phi_tilde[si] = result
    return Phi_tilde


# ===========================================================================
# File I/O
# ===========================================================================

# def save_Phi_data(filename, Ns, Nl, Nmax, Nmax2, PhiM, Phi_hat, Phi_tilde):
#     print("Writing Phi matrices to text file …")
#     jwrite = 0
#     with open(filename, "w") as fh:
#         for si in range(Ns):
#             for m1 in range(Nmax2):
#                 val = PhiM[si][m1]
#                 fh.write(f"{float(val.real)}\t{float(val.imag)}\n") 
#                 jwrite +=1
#         for si in range(Ns):
#             for li in range(Nl):
#                 for ni in range(Nmax2):
#                     val = Phi_hat[si][li][ni]
#                     fh.write(f"{float(val.real)}\t{float(val.imag)}\n")  
#                     jwrite += 1
#         for si in range(Ns):
#             for li in range(Nl):
#                 for ni in range(Nmax):
#                     for mi in range(Nmax):
#                         val = Phi_tilde[si][li][ni][mi]
#                         fh.write(f"{float(val.real)}\t{float(val.imag)}\n")
#                         jwrite += 1
#     print(f"Written {jwrite} lines — done.")

def save_Phi_data(filename, Ns, Nl, Nmax, Nmax2, PhiM, Phi_hat, Phi_tilde):
    print("Writing Phi matrices to text file …")
    jwrite = 0
    with open(filename, "w") as fh:
        for si in range(Ns):
            for m1 in range(Nmax2):
                fh.write(str(PhiM[si][m1].real)+'\n')  # Write line with the precision data
                fh.write(str(PhiM[si][m1].imag)+'\n')  # Write line with the precision data
                jwrite +=1
        for si in range(Ns):
            for li in range(Nl):
                for ni in range(Nmax2):
                    fh.write(str(Phi_hat[si][li][ni]) + "\n")
                    jwrite += 1
        for si in range(Ns):
            for li in range(Nl):
                for ni in range(Nmax):
                    for mi in range(Nmax):
                        fh.write(str(Phi_tilde[si][li][ni][mi]) + "\n")
                        jwrite += 1
    print(f"Written {jwrite} lines — done.")

# ===========================================================================
# Main entry point
# ===========================================================================

def compute_PhiHat_PhiTilde_r(
    Ns: int = 50,
    lmax: int = 21,
    Nmax: int = 50,
    Nmax2: int = 50,
    sa: gp.mpfr = gp.mpfr(-20),
    basefilename: str = "C5Dgmpy2_",
    ncores: int = 4,
):
    """
    Compute PhiHat[si][li][ni] and PhiTilde[si][li][ni][nj] and save to disk.

    Parameters
    ----------
    Ns          : number of s-points
    lmax        : maximum angular momentum (Nl = lmax+1)
    Nmax        : matrix size for Phi_tilde
    Nmax2       : vector size for Phi_hat
    sa          : subtraction point
    basefilename: prefix for the output file
    mq_factor   : mass factor used in the filename
    ncores      : parallel workers (set to your CPU core count)
    """
    Nl  = lmax + 1
    Nl2 = Nl // 2

    filename = (
        f"{basefilename}sa_{int(-sa)}"
        f"_Phisn_PhiHat_sln_PhiTilde_slnm"
        f"_{Ns}_{Nl2}_{Nmax}_{Nmax2}.dat"
    )
    print("Output file:", filename)

    nodes, weights = tanh_sinh_nodes_weights(n=5)
    print(f"Nodes: {len(nodes)}")

    gpi = gp.const_pi()
    xiv = [gpi / 2 / Ns + i * gpi / Ns for i in range(Ns)]
    s   = [4 + (4 - sa) * gp.tan(xiv[i] / 2) ** 2 for i in range(Ns)]

    print("Computing Legendre polynomials …")
    Leg_pol = compute_Leg_pol(Nl, nodes)

    print("Computing zfv …")
    zfv = compute_zfv(s, nodes, sa)

    print("Computing PhiM …")
    PhiM = compute_PhiM(xiv,Nmax2)
    print("PhiM done") 

    print(f"Computing Phi_hat in parallel (ncores={ncores}) …")
    Phi_hat = compute_Phi_hat_parallel(Ns, Nl, Nmax2, weights, zfv, Leg_pol, ncores)
    print("Phi_hat done.")

    print(f"Computing Phi_tilde in parallel (ncores={ncores}) …")
    Phi_tilde = compute_Phi_tilde_parallel(Ns, Nl, Nmax, weights, zfv, Leg_pol, ncores)
    print("Phi_tilde done.")

    save_Phi_data(filename, Ns, Nl, Nmax, Nmax2, PhiM, Phi_hat, Phi_tilde)
    print("All done.")
    return PhiM, Phi_hat, Phi_tilde
