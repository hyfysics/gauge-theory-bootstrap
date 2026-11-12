# Example of usage
# requires gmpy2
# That can be installed with 
# pip install gmpy2

import json
import time
import gmpy2 as gp
from computePhi import computePhi

# 1. Global Setup
PRECISION = 334
gp.get_context().precision = PRECISION
gpi = gp.const_pi()

def load_inputs(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)

    return {
        "M": int(data["M"]),
        "lmax": int(data["lmax"]),
        "Nmax": int(data["Nmax"]),
        "Nmax2": int(data["Nmax2"]),
        "sa": gp.mpfr(data["sa"]),
        "ncores": int(data["ncores"])
    }

# 2. Main Execution Block
if __name__ == "__main__":
    # --- A. Load and Initialize ---
    cfg = load_inputs("python_inputs.json")
    
    M = cfg["M"]
    lmax = cfg["lmax"]
    Nmax = cfg["Nmax"]
    Nmax2 = cfg["Nmax2"]
    sa = cfg["sa"]
    ncores = cfg["ncores"] # Dynamic from Mathematica

    Ns = M
    Nl = lmax + 1
    Nl2 = Nl // 2

    print(f"Parameters Initialized: M={Ns}, Nl={Nl}, Nl2={Nl2}")

    # --- B. Run Computation ---
    start = time.perf_counter()
    
    # Pass variables directly (using ncores from JSON)
    PhiM, Phi_hat, Phi_tilde = computePhi.compute_PhiHat_PhiTilde_r(
        Ns=Ns, 
        lmax=lmax, 
        Nmax=Nmax, 
        Nmax2=Nmax2, 
        sa=sa, 
        basefilename="Phiphys_", 
        ncores=ncores
    )
    
    end = time.perf_counter()
    print(f"Function took {end - start:.6f} seconds.")

    print("All done")