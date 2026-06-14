"""Is the calibration gDNA prior accurate, or does the EM overrun an accurate prior?

Per locus, compare gdna_prior_count (calibration's gDNA split, before EM) against gdna_expected
(ground truth) and gdna_observed (after EM). If prior ~ truth but observed << truth, the leak is
EM-side; if prior << truth, the leak is calibration-side (count-mean bias)."""
import numpy as np
import pandas as pd

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
loc = pd.read_csv(f"{S}/net_flow_per_locus.tsv", sep="\t")

for cond in ["gdna_gdna300_ss_0.99_nrna_none_capture_on",
             "gdna_gdna300_ss_0.50_nrna_none_capture_on",
             "gdna_gdna300_ss_0.99_nrna_none_capture_off"]:
    d = loc[(loc.condition == cond) & (loc.locus_id >= 0)].copy()
    d = d[d.gdna_expected > 100]  # loci with meaningful gDNA
    pri = d.gdna_prior_count
    exp = d.gdna_expected
    obs = d.gdna_observed
    print(f"\n=== {cond} ===")
    print(f"  loci with gdna_expected>100: {len(d)}")
    print(f"  SUM   prior={pri.sum():,.0f}  truth={exp.sum():,.0f}  observed={obs.sum():,.0f}")
    print(f"  prior/truth    = {pri.sum()/exp.sum():.3f}   (1.0 = calibration prior matches truth)")
    print(f"  observed/truth = {obs.sum()/exp.sum():.3f}   (1.0 = EM recovers truth)")
    print(f"  observed/prior = {obs.sum()/pri.sum():.3f}   (<1 = EM leaks AWAY from an accurate prior)")
    # per-locus ratio distribution
    r_pt = (pri / exp).replace([np.inf, -np.inf], np.nan).dropna()
    r_op = (obs / pri.clip(lower=1)).replace([np.inf, -np.inf], np.nan).dropna()
    print(f"  per-locus prior/truth   median={r_pt.median():.3f}  IQR=[{r_pt.quantile(.25):.3f},{r_pt.quantile(.75):.3f}]")
    print(f"  per-locus observed/prior median={r_op.median():.3f}  IQR=[{r_op.quantile(.25):.3f},{r_op.quantile(.75):.3f}]")
