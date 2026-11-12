# Example of usage
# requires gmpy2
# That can be installed with 
# pip install gmpy2

import json
import time
import gmpy2 as gp
from computePhi import computePhi_unphy

# 1. Global Setup
gp.get_context().precision = 334
gpi = gp.const_pi()

def load_inputs(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    # Recombine lists into high-precision complex numbers
    s0v = [gp.mpc(gp.mpfr(str(re), precision=334)) for re in data["sunphys"]]
    
    return {
        "sunphys": s0v,
        "lmax": int(data["lmax"]),
        "Nmax": int(data["Nmax"]),
        "Nmax2": int(data["Nmax2"]),
        "sa": gp.mpfr(data["sa"]),
        "ncores": int(data["ncores"])
    }

# 2. Main Execution Block
if __name__ == "__main__":
    # --- A. Load and Initialize ---
    cfg = load_inputs("python_inputs_u.json")
    
    s0v = cfg["sunphys"]
    lmax = cfg["lmax"]
    Nmax = cfg["Nmax"]
    Nmax2 = cfg["Nmax2"]
    sa = cfg["sa"]
    ncores = cfg["ncores"] # Dynamic from Mathematica

    N0v = len(s0v)
    Nl = lmax + 1
    Nl2 = Nl // 2

    print(f"Parameters Initialized: N0v={N0v}, Nl={Nl}, Nl2={Nl2}")

    # --- B. Run Computation ---
    start = time.perf_counter()
    
    # Pass variables directly (using ncores from JSON)
    PhiM, Phi_hat, Phi_tilde = computePhi_unphy.compute_PhiHat_PhiTilde_c(
        s=s0v, 
        lmax=lmax, 
        Nmax=Nmax, 
        Nmax2=Nmax2, 
        sa=sa, 
        basefilename="Phiunphys_sa_", 
        ncores=ncores 
    )
    
    end = time.perf_counter()
    print(f"Function took {end - start:.6f} seconds.")

    print("All done")
