"""WHERE DOES THE STRAND OVERDISPERSION COME FROM?  The pre-registered outlier test.

OWNER'S HYPOTHESIS (2026-07-28): gDNA is symmetric around 50/50 **by construction** — it is genomic DNA,
there is no strand preference to disperse.  So excess strand variance on a seed we believe is pure gDNA is
not evidence about gDNA; it is evidence that **the seed is contaminated** by transcription the annotated
reference has no record of (antisense, unannotated genes, readthrough).  If so, a **small population of
outliers** should be driving the fit, because the estimator is a pooled ratio of sums with no trimming:

        od_mom = Σ_s excess_var_s / Σ_s var_scale_s          (`gdna_strand._fit_overdispersion`)

⚠ This matters because the strand likelihood is the ONLY intrinsic gDNA/RNA information source
(`CALIBRATION_ARCHITECTURE.md`).  Inflating its overdispersion to absorb annotation error DILUTES our
strongest evidence.  The ceiling (`Beta(2,2)`, od = 0.2) is currently the only thing preventing that.

THE TEST — concentration.  What share of `Σ excess_var` is carried by the top 0.1 % / 1 % / 5 % of seeds?
And what does `od` become if those are trimmed?

THE CONTROL, and it is the decisive half: run the identical analysis on the SYNTHETIC suite, where the
simulator emits only annotated transcripts, so unannotated transcription **cannot** exist.  If the real
data is outlier-driven and the synthetic is not, the diagnosis holds.

FALSIFIER: if `Σ excess_var` is spread evenly across seeds on real data, the fit is measuring something
real and broad, the outlier story is wrong, and trimming would be unjustified.

Run: OMP_NUM_THREADS=1 python scratchpad/od_outliers.py
"""
from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

from rigel.calibration import gdna_strand as GS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-12
_CEIL = GS._MAX_OVERDISPERSION


def capture(run):
    """Run `calibrate` with the two strand fitters wrapped, returning their exact production inputs."""
    grabbed: dict = {}
    g_orig, r_orig = GS.fit_gdna_strand_overdispersion, GS.fit_rna_strand_overdispersion

    def g_spy(sense, total, gdna_weight, rna_sense_frac, **kw):
        grabbed["gdna"] = (np.asarray(sense, float).copy(), np.asarray(total, float).copy(),
                           np.clip(np.asarray(gdna_weight, float), 0, 1).copy(),
                           float(rna_sense_frac), dict(kw))
        return g_orig(sense, total, gdna_weight, rna_sense_frac, **kw)

    def r_spy(sense, total, rna_sense_frac, **kw):
        grabbed["rna"] = (np.asarray(sense, float).copy(), np.asarray(total, float).copy(),
                          None, float(rna_sense_frac), dict(kw))
        return r_orig(sense, total, rna_sense_frac, **kw)

    GS.fit_gdna_strand_overdispersion = g_spy
    GS.fit_rna_strand_overdispersion = r_spy
    try:
        run()
    finally:
        GS.fit_gdna_strand_overdispersion = g_orig
        GS.fit_rna_strand_overdispersion = r_orig
    return grabbed


def terms(sense, total, weight, kappa, component_mean):
    """The estimator's own per-seed numerator and denominator (`_fit_overdispersion`, verbatim)."""
    if weight is None:                      # RNA: pure-spliced seeds, component_frac = 1
        node_mean = np.full(total.shape, kappa)
        frac = np.ones_like(total)
    else:
        node_mean = 0.5 * weight + kappa * (1.0 - weight)
        frac = weight
    binom = total * node_mean * (1.0 - node_mean)
    excess = (sense - total * node_mean) ** 2 - binom
    nc = frac * total
    scale = np.maximum(nc * (nc - 1.0), 0.0) * (component_mean * (1.0 - component_mean))
    ok = total > 0.0
    return excess, scale, ok


def report(tag, sense, total, weight, kappa, component_mean, prior_od, prior_w):
    ex, sc, ok = terms(sense, total, weight, kappa, component_mean)
    ex, sc = ex[ok], sc[ok]
    n = ex.size
    if n == 0 or sc.sum() <= 0:
        print(f"{tag:34s}  (no seeds)")
        return
    order = np.argsort(-ex)                  # most positive excess first — what inflates od
    cum = np.cumsum(ex[order]) / max(ex.sum(), _EPS)

    def share_at(p):
        k = max(int(np.ceil(p * n)), 1)
        return cum[k - 1]

    def od_trim(p):
        """od with the top p of seeds by excess_var removed — the estimator, honestly robustified."""
        k = int(np.ceil(p * n))
        keep = np.ones(n, bool)
        keep[order[:k]] = False
        num, den = ex[keep].sum(), sc[keep].sum()
        mom = num / den if den > 0 else 0.0
        nn = int(keep.sum())
        od = (nn * mom + prior_w * prior_od) / (nn + prior_w) if (nn + prior_w) > 0 else prior_od
        return float(np.clip(od, 0.0, _CEIL))

    print(f"{tag:34s} {n:7d} {od_trim(0.0):8.4f} | {share_at(0.001):7.1%} {share_at(0.01):7.1%} "
          f"{share_at(0.05):7.1%} | {od_trim(0.001):8.4f} {od_trim(0.01):8.4f} {od_trim(0.05):8.4f}")


def main():
    cfg = PipelineConfig()
    pri_g = GS.overdispersion_for_beta(cfg.calibration.gdna_strand_prior_alpha_beta)
    pri_r = GS.overdispersion_for_beta(cfg.calibration.rna_strand_prior_alpha_beta)
    hdr = (f"{'sample':34s} {'n seed':>7s} {'od':>8s} | {'top0.1%':>7s} {'top1%':>7s} {'top5%':>7s} | "
           f"{'od-0.1%':>8s} {'od-1%':>8s} {'od-5%':>8s}")

    print("=== gDNA STRAND OVERDISPERSION — share of the MoM numerator carried by the top seeds ===")
    print("    'od-x%' = the fit with the top x% of seeds by excess variance removed.")
    print(f"    ceiling = {_CEIL:.4f} (Beta(2,2)).  gDNA's true strand mean is 1/2 BY CONSTRUCTION.\n")
    print(hdr)
    keep = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        keep[name] = g
        s, t, w, k, _ = g["gdna"]
        report(f"REAL  {name}", s, t, w, k, 0.5, pri_g, cfg.calibration.gdna_strand_prior_weight)

    # ── the control: the synthetic suite CANNOT contain unannotated transcription ──────────────────
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
        keep[cond] = g
        s, t, w, k, _ = g["gdna"]
        report(f"SYNTH {cond[5:29]}", s, t, w, k, 0.5, pri_g, cfg.calibration.gdna_strand_prior_weight)

    print("\n\n=== RNA STRAND OVERDISPERSION (pure-spliced seeds, component mean kappa) ===\n")
    print(hdr)
    for name, g in keep.items():
        if "rna" not in g:
            continue
        s, t, _, k, _ = g["rna"]
        lab = ("REAL  " + name) if (CF / f"{name}.pkl").exists() else ("SYNTH " + name[5:29])
        report(lab, s, t, None, k, k, pri_r, cfg.calibration.rna_strand_prior_weight)

    print("\n\n=== WHAT ARE THE TOP gDNA SEEDS?  (real data, the worst sample) ===\n")
    worst = max((n for n in keep if (CF / f"{n}.pkl").exists()),
                key=lambda n: terms(*keep[n]["gdna"][:3], keep[n]["gdna"][3], 0.5)[0].max())
    s, t, w, k, _ = keep[worst]["gdna"]
    ex, sc, ok = terms(s, t, w, k, 0.5)
    idx = np.argsort(-np.where(ok, ex, -np.inf))[:12]
    print(f"  sample = {worst}   (total seeds {int(ok.sum()):,})")
    print(f"  {'excess_var':>14s} {'share of num':>13s} {'N':>9s} {'sense':>9s} {'sense frac':>11s} "
          f"{'gdna wt':>8s}")
    tot = ex[ok].sum()
    for i in idx:
        i = int(i)
        print(f"  {ex[i]:14,.0f} {ex[i] / max(tot, _EPS):12.2%} {t[i]:9,.0f} {s[i]:9,.0f} "
              f"{s[i] / max(t[i], _EPS):11.3f} {w[i]:8.3f}")
    print("\n  A pure-gDNA seed should sit at sense frac ~ 0.5. A seed far from 0.5 with large N is")
    print("  transcription the annotation does not record — exactly what must NOT set the dispersion.")


if __name__ == "__main__":
    main()
