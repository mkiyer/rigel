"""Plot the pass-0 gDNA-rate NPMLE prior P(ρ) for every condition of a sim suite (a grid, for comparison).

Scans + calibrates each condition (capturing the fitted `DensityNPMLE` via the calibrate `_debug` hook) and
plots its P(log ρ) — so we can see the pass-0 prior's shape (depleted / middle / enriched modes) across gDNA
level, strand specificity, capture, and nascent. Curves are coloured by the condition's true gDNA level.

    OMP_NUM_THREADS=1 python scripts/debug/plot_pass0_priors.py --suite DIR --out FIG.png
"""

from __future__ import annotations

import argparse
import os
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

_LN10 = np.log(10.0)
_GDNA_COLOR = {"none": "tab:green", "gdna100": "tab:orange", "gdna300": "tab:red"}


def _fit_prior(bam, index, cfg):
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _s, sm, flm, _b, payload = scan_and_buffer(bam, index, sc)
    sm.finalize()
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    dbg: dict = {}
    calibrate(
        payload=payload, region_arrays=ra, strand_model=sm,
        gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration, _debug=dbg,
    )
    return dbg["gdna_prior"]


def _gdna_level(cond):
    for k in ("gdna300", "gdna100", "none"):
        if k in cond:
            return k
    return "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    suite = Path(args.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    conds = sorted(
        p.name for p in suite.iterdir() if p.is_dir() and (p / "sim_oracle.bam").exists()
    )
    n = len(conds)
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow), squeeze=False)
    for i, cond in enumerate(conds):
        ax = axes[i // ncol][i % ncol]
        try:
            prior = _fit_prior(str(suite / cond / "sim_oracle.bam"), index, cfg)
            x = prior.log_rho / _LN10
            p = np.exp(prior.logP - prior.logP.max())
            ax.fill_between(x, p, color=_GDNA_COLOR[_gdna_level(cond)], alpha=0.55)
            ax.plot(x, p, color="black", lw=0.7)
            ax.set_yticks([])
            ax.set_title(cond.replace("gdna_", "").replace("_capture", "_cap"), fontsize=6.5)
            ax.tick_params(axis="x", labelsize=6)
            print(f"{cond:48s} cells={prior.n_cells}", flush=True)
        except Exception as e:  # noqa: BLE001 — a diagnostic; keep plotting the rest
            ax.set_title(f"{cond}\nFAIL: {e}", fontsize=6)
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"pass-0 gDNA-rate NPMLE prior P(log ρ) — {suite.name}", fontsize=11)
    fig.supxlabel("log10 gDNA rate ρ", fontsize=9)
    fig.tight_layout(rect=(0, 0.01, 1, 0.99))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"-> {args.out}", flush=True)


if __name__ == "__main__":
    main()
