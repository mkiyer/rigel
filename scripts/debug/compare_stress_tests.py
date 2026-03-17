"""Compare two stress test results (v3 vs v5)."""
import pandas as pd

v3 = pd.read_csv("results/stress_test_v3/stress_test_results.csv")
v5 = pd.read_csv("results/stress_test_v5/stress_test_results.csv")
m = v3.merge(v5, on="label", suffixes=("_v3", "_v5"))

# Overall
print("=== OVERALL ===")
for col in ["auroc", "gdna_f1", "rna_f1", "pi_error"]:
    c3 = col + "_v3"
    c4 = col + "_v5"
    print(f"  {col:12s}: v3={m[c3].mean():.4f}  v4={m[c4].mean():.4f}  delta={m[c4].mean()-m[c3].mean():+.4f}")

# By SS
print("\n=== BY STRAND SPECIFICITY ===")
for ss in sorted(m["strand_specificity_v3"].unique()):
    sub = m[m["strand_specificity_v3"] == ss]
    n = len(sub)
    g3 = sub["gdna_f1_v3"].mean()
    g4 = sub["gdna_f1_v5"].mean()
    r3 = sub["rna_f1_v3"].mean()
    r4 = sub["rna_f1_v5"].mean()
    print(f"  SS={ss:.2f} (n={n:2d}): gF1 {g3:.3f}->{g4:.3f} ({g4-g3:+.3f})  rF1 {r3:.3f}->{r4:.3f} ({r4-r3:+.3f})")

# By FL overlap
print("\n=== BY FL OVERLAP (gDNA FL mean) ===")
for fl in sorted(m["gdna_fl_mean_v3"].unique()):
    sub = m[m["gdna_fl_mean_v3"] == fl]
    n = len(sub)
    g3 = sub["gdna_f1_v3"].mean()
    g4 = sub["gdna_f1_v5"].mean()
    r3 = sub["rna_f1_v3"].mean()
    r4 = sub["rna_f1_v5"].mean()
    print(f"  FL={fl:.0f} (n={n:2d}): gF1 {g3:.3f}->{g4:.3f} ({g4-g3:+.3f})  rF1 {r3:.3f}->{r4:.3f} ({r4-r3:+.3f})")

# Regressions and improvements
m["gdna_delta"] = m["gdna_f1_v5"] - m["gdna_f1_v3"]

print("\n=== NOTABLE REGRESSIONS (v4 gF1 dropped >0.05) ===")
reg = m[m["gdna_delta"] < -0.05].sort_values("gdna_delta")
for _, r in reg.head(15).iterrows():
    sc = r["label"]
    print(f"  {sc:55s}: gF1 {r['gdna_f1_v3']:.3f}->{r['gdna_f1_v5']:.3f} ({r['gdna_delta']:+.3f})")

print("\n=== NOTABLE IMPROVEMENTS (v4 gF1 improved >0.05) ===")
imp = m[m["gdna_delta"] > 0.05].sort_values("gdna_delta", ascending=False)
for _, r in imp.head(15).iterrows():
    sc = r["label"]
    print(f"  {sc:55s}: gF1 {r['gdna_f1_v3']:.3f}->{r['gdna_f1_v5']:.3f} ({r['gdna_delta']:+.3f})")

# Count wins/ties/losses
wins = (m["gdna_delta"] > 0.01).sum()
losses = (m["gdna_delta"] < -0.01).sum()
ties = len(m) - wins - losses
print(f"\nWin/Tie/Loss (|delta|>0.01): {wins}/{ties}/{losses}")

# FL-identical subset
print("\n=== FL-IDENTICAL SUBSET (hardest cases) ===")
fl_id = m[m["gdna_fl_mean_v3"] == 250]
print(f"  n={len(fl_id)}")
g3 = fl_id["gdna_f1_v3"].mean()
g4 = fl_id["gdna_f1_v5"].mean()
r3 = fl_id["rna_f1_v3"].mean()
r4 = fl_id["rna_f1_v5"].mean()
print(f"  gF1: {g3:.3f} -> {g4:.3f} ({g4-g3:+.3f})")
print(f"  rF1: {r3:.3f} -> {r4:.3f} ({r4-r3:+.3f})")
for ss in sorted(fl_id["strand_specificity_v3"].unique()):
    sub = fl_id[fl_id["strand_specificity_v3"] == ss]
    g3s = sub["gdna_f1_v3"].mean()
    g4s = sub["gdna_f1_v5"].mean()
    print(f"    SS={ss:.2f} (n={len(sub):2d}): gF1 {g3s:.3f}->{g4s:.3f} ({g4s-g3s:+.3f})")
