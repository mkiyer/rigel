"""Quick test of FL model behavior under different data conditions."""
import numpy as np
from rigel.frag_length_model import FragmentLengthModel

# Case 1: Empty model, finalized
m = FragmentLengthModel(max_size=1000)
m.finalize()
print(f"Empty finalized: total_weight={m.total_weight}")
print(f"  log_prob[200]={m._log_prob[200]:.4f}")
print(f"  uniform: -log(1001)={-np.log(1001):.4f}")

# Case 2: Never finalized → _log_prob is None
m2 = FragmentLengthModel(max_size=1000)
print(f"\nUnfinalized: _log_prob is None = {m2._log_prob is None}")

# Case 3: Sparse data (11 obs)
m3 = FragmentLengthModel(max_size=1000)
for fl in [180, 190, 195, 200, 205, 210, 220, 230, 185, 215, 225]:
    m3.observe(fl)
m3.finalize()
print(f"\nSparse (11 obs): total_weight={m3.total_weight}")
print(f"  log_prob[200]={m3._log_prob[200]:.4f}")
print(f"  log_prob[500]={m3._log_prob[500]:.4f}")
print(f"  log_prob[0]={m3._log_prob[0]:.4f}")

# Case 4: Rich data (1000 obs)
np.random.seed(42)
m4 = FragmentLengthModel(max_size=1000)
flens = np.random.normal(200, 30, 1000).astype(int)
for fl in flens:
    if 0 <= fl <= 1000:
        m4.observe(fl)
m4.finalize()
print(f"\nRich (1000 obs): total_weight={m4.total_weight}")
print(f"  log_prob[200]={m4._log_prob[200]:.4f}")
print(f"  log_prob[500]={m4._log_prob[500]:.4f}")
print(f"  log_prob[0]={m4._log_prob[0]:.4f}")

# Key comparison: what does C++ get when _log_prob is None?
print(f"\n--- C++ behavior ---")
print(f"When _log_prob is None (passed as None to constructor):")
print(f"  frag_len_log_lik() returns 0.0 (has_fl_lut_ = false)")
print(f"When _log_prob is a flat array of -6.91:")
print(f"  frag_len_log_lik() returns -6.91 (has_fl_lut_ = true)")
print(f"  ASYMMETRY: 6.91 log-unit advantage for the missing model!")
