"""DESIGNING THE ROBUST STRAND-OVERDISPERSION ESTIMATOR — the candidate bake-off.

⚠ TWO EARLIER CONCLUSIONS WERE ARTIFACTS OF SEED SIZE, and this script exists because of them:

  * "60-87 % of seeds are DISPLACED (|sense frac − ½| > 0.3)" — the median such seed has **N = 1-2
    fragments**. With N = 2 the only outcomes are sf ∈ {0, ½, 1}, so ¾ of *fair coins* land there. That is
    counting, not contamination.
  * "the MEDIAN of per-seed moments still saturates, so >50 % of seeds are overdispersed" — with n_c = 2
    the per-seed ratio `excess/scale` is exactly ±1, so the unweighted median of a ±1 discrete population
    is ±1 or 0 and carries no information. The median was measuring seed size, not dispersion.

WHAT IS ACTUALLY TRUE (`od_shape.py`): median z = 0.000 on every real sample and IQR(z)/1.349 = 1.3-1.5,
which small-N discreteness alone explains. The seed population is NOT broadly overdispersed. The pooled
estimate is decided by a HANDFUL OF LARGE-N SEEDS, because the estimator's denominator is `n_c(n_c−1)` —
fragment-SQUARED weighting with no bound on any single seed's influence. In transcribed territory the
large-N seeds are transcripts.

So the design problem is bounded influence, and the right currency is INFORMATION, not seed count.

CANDIDATES (all estimate the same quantity; none introduces a tuned constant except where marked):
  pooled       the shipped estimator, Σexcess / Σscale
  wmedian      the INFORMATION-WEIGHTED median of the per-seed ratio (weights = scale_s)
  wtrim        information-weighted mean of the central 50 % by ratio (a weighted IQR mean)
  huber        bounded-influence: cap each seed's weight share of Σscale at 1/sqrt(n_informative)
  zscreen+X    any of the above after dropping seeds whose own |z| contradicts Binomial(½)

Truth is known on the synthetic suite (the simulator is multinomial ⇒ od = 0 by construction), so the
synthetic rows are the accuracy check: anything materially above 0 there is the estimator inventing
dispersion.

Run: OMP_NUM_THREADS=1 python scratchpad/od_design.py
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


def wquantile(v, w, q):
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    if c[-1] <= 0:
        return np.nan
    return float(np.interp(q * c[-1], c, v))


def estimators(sense, total, weight, kappa):
    ex, sc, ok = terms(sense, total, weight, kappa, 0.5)
    node_mean = 0.5 * weight + kappa * (1.0 - weight)
    sd = np.sqrt(np.maximum(total * node_mean * (1.0 - node_mean), _EPS))
    z = (sense - total * node_mean) / sd
    ex, sc, z = ex[ok], sc[ok], z[ok]
    inf = sc > 0                       # a seed with n_c < 2 carries NO dispersion information
    ex, sc, z = ex[inf], sc[inf], z[inf]
    if ex.size == 0 or sc.sum() <= 0:
        return {}
    r = ex / np.maximum(sc, _EPS)
    out = {}
    out["pooled"] = ex.sum() / sc.sum()
    out["wmedian"] = wquantile(r, sc, 0.5)
    lo, hi = wquantile(r, sc, 0.25), wquantile(r, sc, 0.75)
    m = (r >= lo) & (r <= hi)
    out["wtrim"] = (ex[m].sum() / sc[m].sum()) if sc[m].sum() > 0 else np.nan
    cap = sc.sum() / max(np.sqrt(sc.size), 1.0)          # no single seed may exceed 1/sqrt(n) of the info
    scc = np.minimum(sc, cap)
    out["huber"] = float(np.sum(r * scc) / max(scc.sum(), _EPS))
    keep = np.abs(z) <= 3.0
    out["zscr+pooled"] = (ex[keep].sum() / sc[keep].sum()) if sc[keep].sum() > 0 else np.nan
    out["zscr+wmed"] = wquantile(r[keep], sc[keep], 0.5) if keep.any() else np.nan
    out["_n_inf"] = float(ex.size)
    out["_n_all"] = float(ok.sum())
    return out


NAMES = ["pooled", "wmedian", "wtrim", "huber", "zscr+pooled", "zscr+wmed"]


def main():
    cfg = PipelineConfig()
    rows = {}
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        rows[f"REAL  {name}"] = estimators(*g["gdna"][:3], g["gdna"][3])
    from selfsolve_diag import _scan_and_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    for cond in ("gdna_gdna100_ss_0.99_nrna_present_capture_on",
                 "gdna_gdna100_ss_0.50_nrna_none_capture_off",
                 "gdna_gdna300_ss_0.99_nrna_none_capture_off"):
        inp = _scan_and_truth(SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        g = capture(lambda inp=inp: calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        rows[f"SYNTH {cond[5:27]}"] = estimators(*g["gdna"][:3], g["gdna"][3])

    print("=== HOW MANY SEEDS CARRY ANY DISPERSION INFORMATION AT ALL (n_c >= 2)? ===\n")
    print(f"{'sample':30s} {'seeds':>9s} {'informative':>12s} {'share':>8s}")
    for tag, r in rows.items():
        if not r:
            continue
        print(f"{tag:30s} {r['_n_all']:9,.0f} {r['_n_inf']:12,.0f} "
              f"{r['_n_inf'] / max(r['_n_all'], 1):8.1%}")

    print("\n\n=== THE BAKE-OFF (raw estimates, BEFORE the prior shrinkage and the ceiling clip) ===")
    print("    SYNTH truth = 0 by construction (multinomial simulator). Ceiling = 0.2000.\n")
    print(f"{'sample':30s}" + "".join(f"{n:>13s}" for n in NAMES))
    for tag, r in rows.items():
        if not r:
            continue
        print(f"{tag:30s}" + "".join(f"{r.get(n, np.nan):13.4f}" for n in NAMES))

    print("\n  pooled      = shipped. wmedian/wtrim/huber = information-weighted robust variants.")
    print("  zscr+X      = X after dropping seeds whose own |z| > 3 contradicts the Binomial(1/2)")
    print("                premise they enter a mean-1/2 fit on.")


if __name__ == "__main__":
    main()
