"""Plot the COMBINED pass-0 gDNA prior — NPMLE + pinned background component + the one-sided floor anchor —
against the oracle gDNA density distribution, per condition.

Per panel (log10 ρ axis):
  * light fill  — the TOTAL unspliced density M/E (what the pass-0 NPMLE is fit on; f_g=1);
  * dashed blue — the plain NPMLE prior (no background);
  * solid blue  — the COMBINED prior density (NPMLE + pinned component at ρ_bg);
  * red         — the ORACLE gDNA density G/E distribution (the truth the deconvolution must reach);
  * green line  — log ρ_bg, the floor/anchor location; the shaded band left of it is the one-sided floor's
    penalty zone (ρ_g < ρ_bg pushed up).

The gap between the pass-0 prior (total-based) and the oracle gDNA is exactly what the deconvolution + refit
must close; the floor/pinned anchor pins the low end at the measured background.

    OMP_NUM_THREADS=1 python background_prior_plot.py [--out DIR]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
import refit_loop_study as R  # noqa: E402
from rigel.calibration.background_reference import measure_background  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_LN10 = np.log(10.0)
_EPS = 1e-12

CONDS = [
    ("gdna_gdna300_ss_0.50_nrna_present_capture_on", "FLAGSHIP: gdna300 unstr capON"),
    ("gdna_gdna300_ss_0.50_nrna_present_capture_off", "gdna300 unstr capOFF"),
    ("gdna_gdna300_ss_0.99_nrna_present_capture_on", "stranded gdna300 capON"),
    ("gdna_gdna5_ss_0.50_nrna_present_capture_on", "low-DNA: gdna5 unstr capON"),
    ("gdna_gdna1_ss_0.50_nrna_present_capture_verystrong", "gdna1 unstr VERYSTRONG"),
    ("gdna_none_ss_0.50_nrna_present_capture_on", "ZERO-gDNA: none unstr capON"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    axes = axes.ravel()
    for ax, (c, title) in zip(axes, CONDS):
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        s = R._setup(inp, index, cfg)
        M, E, G, Rn = s["mass_g"], s["eff_g"], s["G"], s["R"]
        live = (E > 1e-9 * 1.001) & (M > _EPS)
        ra = s["ra"]
        sub = CalibrationSubstrate.from_payload(inp["payload"], ra)
        reg_el = region_eff_length(ra.region_size_bp, np.asarray(inp["gdna_fl_pmf"]))
        bg = measure_background(sub, ra, reg_el)
        pr_bg = DensityNPMLE.fit(M, E, background=bg, bandwidth=cfg.calibration.npmle_bandwidth)
        pr_0 = DensityNPMLE.fit(M, E, bandwidth=cfg.calibration.npmle_bandwidth)

        tot = (M[live] / np.maximum(E[live], _EPS))
        ltot = np.log10(tot[tot > 0])
        gdl = G[live] / np.maximum(E[live], _EPS)
        lg = np.log10(gdl[gdl > 0])
        lo = min(ltot.min(), (lg.min() if lg.size else ltot.min())) - 0.5
        hi = ltot.max() + 0.5
        bins = np.linspace(lo, hi, 60)
        ax.hist(ltot, bins=bins, density=True, color="0.85", label="total M/E")
        if lg.size:
            ax.hist(lg, bins=bins, density=True, histtype="step", color="tab:red", lw=1.8,
                    label="oracle gDNA G/E")
        # prior densities (normalise the peak to the hist scale for visual comparison)
        def _curve(pr, **kw):
            x = pr.log_rho / _LN10
            P = np.exp(pr.logP - pr.logP.max())
            m = (x >= lo) & (x <= hi)
            ax.plot(x[m], P[m] * ax.get_ylim()[1] * 0.9, **kw)
        _curve(pr_0, color="tab:blue", ls="--", lw=1.3, label="NPMLE (no bg)")
        _curve(pr_bg, color="tab:blue", ls="-", lw=2.0, label="COMBINED (NPMLE+pinned)")
        if np.isfinite(bg.log_rho_bg):
            xb = bg.log_rho_bg / _LN10
            ax.axvline(xb, color="tab:green", lw=1.6, label=f"ρ_bg anchor ({xb:.2f})")
            ax.axvspan(lo, xb, color="tab:green", alpha=0.06)  # the floor's penalty zone
        else:
            ax.text(0.03, 0.9, "ρ_bg = 0 (floor DORMANT)", transform=ax.transAxes, color="tab:green",
                    fontsize=9)
        fgt = G[live].sum() / max(G[live].sum() + Rn[live].sum(), 1.0)
        ax.set_title(f"{title}\n[true gDNAfrac={fgt:.2f}, Σg_bg={bg.n_counts:.0f}]", fontsize=9)
        ax.set_xlabel("log10 ρ", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper right")
        ax.tick_params(labelsize=7)
        print(f"{c:52s} plotted (ρ_bg log10={bg.log_rho_bg/_LN10 if np.isfinite(bg.log_rho_bg) else -99:.2f})",
              flush=True)
    fig.suptitle(
        "Combined pass-0 gDNA prior (NPMLE + pinned bg + floor anchor) vs oracle gDNA density. "
        "Gap total→oracle-gDNA = what deconvolution/refit closes; green = ρ_bg anchor.",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = outdir / "background_prior_plot.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
