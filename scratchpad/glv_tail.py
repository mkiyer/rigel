"""Is the tail of E[d^2] driven by SPARSE (noisy) pairs or by WELL-COUNTED ones?

This decides two things at once:
  * whether inverse-variance weighting (proper DerSimonian-Laird) would help the fit or throw away signal;
  * whether a per-node omega could ever be better than the pooled one (it can only be if the large |d|
    pairs are well-counted, i.e. their disagreement is REAL and not a counting artifact).

    OMP_NUM_THREADS=1 python scratchpad/glv_tail.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    f = dbg["capture"]["_glv"][0]  # + strand, relay fit
    ok = f["ok"]
    d2, v = f["d"][ok] ** 2, f["noise"][ok]
    n = d2.shape[0]
    order = np.argsort(-d2)
    top = order[: max(n // 10, 1)]
    print(f"\n=== {cond[5:]}   n={n}  omega={f['omega']:.4f}  E[d2]={d2.mean():.4f} "
          f"E[noise]={v.mean():.4f}")
    print(f"  top decile of d^2 carries {d2[top].sum() / d2.sum():>6.1%} of sum(d^2); "
          f"its mean NOISE is {v[top].mean():.4f} vs {v.mean():.4f} overall "
          f"({v[top].mean() / max(v.mean(), 1e-9):.2f}x)")
    print(f"  corr(d^2, noise) = {np.corrcoef(d2, v)[0, 1]:+.3f}    "
          f"Spearman = {np.corrcoef(np.argsort(np.argsort(d2)), np.argsort(np.argsort(v)))[0, 1]:+.3f}")
    # excess per noise-quartile: is the PREMISE error itself count-dependent?
    q = np.quantile(v, [0.25, 0.5, 0.75])
    print(f"  {'noise quartile':<18}{'n':>5}{'E[d2]':>9}{'E[noise]':>10}{'-> omega':>10}")
    lo = -np.inf
    for k, hi in enumerate(list(q) + [np.inf]):
        m = (v > lo) & (v <= hi)
        if m.sum():
            print(f"  Q{k + 1:<17}{int(m.sum()):>5}{d2[m].mean():>9.4f}{v[m].mean():>10.4f}"
                  f"{max(0.0, d2[m].mean() - v[m].mean()) / 2:>10.4f}")
        lo = hi
    # what a proper inverse-variance (DerSimonian-Laird) fit would give instead of the plain difference
    w = 1.0 / np.maximum(v, 1e-6)
    Q = float((w * d2).sum())
    c = float(w.sum() - (w * w).sum() / w.sum())
    dl = max(0.0, (Q - (n - 1)) / max(c, 1e-9)) / 2
    print(f"  plain moment omega = {f['omega']:.4f}   DerSimonian-Laird omega = {dl:.4f}")
