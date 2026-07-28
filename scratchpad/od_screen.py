"""A PRINCIPLED SEED SCREEN for the strand overdispersion — reject on the seed's OWN premise.

`od_outliers.py` showed the top gDNA seeds on real data are things like `N = 1,523, sense = 5`
(sense fraction **0.003**) carrying `gdna_weight = 1.000`.  gDNA cannot be 0.3 % sense: it is genomic DNA,
symmetric around ½ by construction.  Such a seed is transcription the annotation does not record.

So the screen should NOT be a percentile of excess variance (that is a knob, and `od_outliers.py` measured
that trimming the top 1–5 % fails to release two of four samples).  It should be the seed's **own premise**:

> a seed entering a **mean-½** gDNA fit claims its gDNA fragments are Binomial(n_g, ½).
> Test that claim.  A seed that contradicts it is not pure gDNA and must not set the dispersion.

The test statistic is the seed's own standardized deviation from its claimed mean,
`z = (sense − N·mean) / sqrt(N·mean·(1−mean))`, which is what the estimator already forms — `excess_var`
is literally `N·mean(1−mean)·(z² − 1)`.  Rejecting on |z| is therefore not a new quantity, it is a
DIRECTION on the one the estimator already computes: the estimator currently *sums* z²−1 and lets the
outliers win; a robust version asks whether each seed is consistent before pooling it.

⚠ A cutoff on |z| IS a constant and would need owner sign-off. This script does not propose one — it
measures the whole curve so the shape of the decision is visible, and reports what the classical
robust choice (the MEDIAN of z²−1, which needs no cutoff) would give.

Run: OMP_NUM_THREADS=1 python scratchpad/od_screen.py
"""
from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from od_outliers import CF, SUITE, capture, terms  # noqa: E402

from rigel.calibration import gdna_strand as GS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

_EPS = 1e-12
_CEIL = GS._MAX_OVERDISPERSION


def od_from(ex, sc, keep, prior_od, prior_w):
    num, den = ex[keep].sum(), sc[keep].sum()
    mom = num / den if den > 0 else 0.0
    n = int(keep.sum())
    od = (n * mom + prior_w * prior_od) / (n + prior_w) if (n + prior_w) > 0 else prior_od
    return float(np.clip(od, 0.0, _CEIL))


def analyse(tag, sense, total, weight, kappa, prior_od, prior_w):
    node_mean = 0.5 * weight + kappa * (1.0 - weight)
    sd = np.sqrt(np.maximum(total * node_mean * (1.0 - node_mean), _EPS))
    z = (sense - total * node_mean) / sd
    ex, sc, ok = terms(sense, total, weight, kappa, 0.5)
    ex, sc, z, w_, t_ = ex[ok], sc[ok], z[ok], weight[ok], total[ok]
    n = ex.size
    base = od_from(ex, sc, np.ones(n, bool), prior_od, prior_w)
    row = [f"{tag:26s} {n:7d} {base:8.4f}"]
    for zc in (10.0, 5.0, 3.0):
        row.append(f"{od_from(ex, sc, np.abs(z) <= zc, prior_od, prior_w):8.4f}")
    # the cutoff-free classical alternative: the MEDIAN of the per-seed moment ratio
    good = sc > 0
    med = float(np.median((ex[good] / sc[good]))) if good.any() else 0.0
    nn = int(good.sum())
    od_med = float(np.clip((nn * med + prior_w * prior_od) / (nn + prior_w), 0.0, _CEIL))
    row.append(f"{od_med:9.4f}")
    # how much is rejected, and what it looks like
    rej = np.abs(z) > 5.0
    frac_sense = np.where(t_ > 0, sense[ok] / np.maximum(t_, _EPS), np.nan)
    row.append(f"{rej.mean():8.2%} {np.median(frac_sense[rej]) if rej.any() else np.nan:9.3f}"
               f" {np.median(w_[rej]) if rej.any() else np.nan:8.2f}")
    print("".join(row))


def main():
    cfg = PipelineConfig()
    pri = GS.overdispersion_for_beta(cfg.calibration.gdna_strand_prior_alpha_beta)
    pw = cfg.calibration.gdna_strand_prior_weight
    print("=== gDNA overdispersion under a SEED-CONSISTENCY screen (reject |z| > zc), and the")
    print("    cutoff-free MEDIAN-of-moments alternative.  Ceiling = 0.2000.\n")
    print(f"{'sample':26s} {'n seed':>7s} {'od now':>8s} {'|z|<=10':>8s} {'|z|<=5':>8s} {'|z|<=3':>8s} "
          f"{'MEDIAN':>9s} | {'%rej':>8s} {'their':>9s} {'their':>8s}")
    print(f"{'':26s} {'':>7s} {'':>8s} {'':>8s} {'':>8s} {'':>8s} {'(no cutoff)':>9s} | "
          f"{'(z>5)':>8s} {'sensefr':>9s} {'gdna wt':>8s}")
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, w, k, _ = g["gdna"]
        analyse(f"REAL  {name}", s, t, w, k, pri, pw)

    from selfsolve_diag import _scan_and_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    for cond in ("gdna_gdna100_ss_0.99_nrna_present_capture_on",
                 "gdna_gdna100_ss_0.50_nrna_none_capture_off"):
        inp = _scan_and_truth(SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        g = capture(lambda inp=inp: calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, w, k, _ = g["gdna"]
        analyse(f"SYNTH {cond[5:24]}", s, t, w, k, pri, pw)

    print("\n  'their sensefr' = median sense fraction of the REJECTED seeds. A pure-gDNA seed must sit")
    print("  at ~0.5; anything far from it is transcription, and 'their gdna wt' shows the count module")
    print("  had labelled it ~100 % gDNA. The screen is testing the seed's OWN premise, not a threshold")
    print("  on the answer.")


if __name__ == "__main__":
    main()
