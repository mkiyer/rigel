"""A/B the factor-1 bedrock + the two toy arms across HEAD (log r)^2 vs the working-tree M5+DL."""
import importlib.util, sys, os
import numpy as np

HEAD = "/Users/mkiyer/proj/rigel/scratchpad/audit_dledit_head_bp_solver.py"

def load_head():
    spec = importlib.util.spec_from_file_location("rigel.calibration.bp_solver", HEAD)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["rigel.calibration.bp_solver"] = mod
    spec.loader.exec_module(mod)

if os.environ.get("ARM") == "head":
    import rigel.calibration  # ensure the package exists first
    load_head()

sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests/calibration")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests")
import test_bp_solver as T

l, r = T._factor1_uniform_rho()
print(f"ARM={os.environ.get('ARM','worktree')}")
print(f"  factor1 intergenic rho_l={l[[1,5]]}  (truth 0.5)")
print(f"  factor1 AMBIG      rho_l={l[3]:.4f} rho_r={r[3]:.4f}  (truth 0.5, tol 0.05)")
fin, cap = T._sweep(T._mature_exon_chain(spliced=True))
print(f"  mature toy: exon f_g={fin.f_g[3]:.4f} (truth~0.32)  introns f_g={fin.f_g[[1,5]]} (truth 1.0)")
ok, _ = T._sweep(T._mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=1.0))
lo, _ = T._sweep(T._mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=0.25))
print(f"  disagreement-silenced |dfg|={abs(float(lo.f_g[3])-float(ok.f_g[3])):.4f} (strict-xfail bound 0.05)")
