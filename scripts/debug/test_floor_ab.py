"""FULL-SOLVE A/B of Phase-2: the one-sided background FLOOR in `DensityNPMLE.logprior`.

Per condition: measure the background (`measure_background`, intergenic-only), then run the pass-0 solve with
the prior built floor-OFF (`background=None`, today's behavior) vs floor-ON (`background=bg`), and compare the
per-region mass-weighted |f_g − oracle| (mwae). The floor is default-off in production, so this A/B is the gate
before wiring `calibrate` to pass the background.

Gates: (1) gdna_none stays ≈ current at every capture (dormant floor, ρ_bg=0); (2) low-DNA/strong-capture rows
do NOT manufacture gDNA (mwae must not jump); (3) off-target/depleted crush relieved where ρ_bg is substantial;
(4) enriched flagship: partial improvement only (the honest limit — the floor is small for enriched nodes).

    OMP_NUM_THREADS=1 python test_floor_ab.py
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

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


def _pass0(s, bg):
    prior = DensityNPMLE.fit(s["mass_g"], s["eff_g"], background=bg, bandwidth=0.15)
    belief = R._sweep(s, s["b0"], prior, transfer_variance=True)
    return R._measure(s, belief, prior, None)


def _group(c):
    return ("gNONE" if "gdna_none" in c else "g1" if "gdna_gdna1_" in c else "g5" if "gdna_gdna5_" in c
            else "g100" if "gdna100" in c else "g300")


def main():
    suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = sorted(p.name for p in suite.iterdir() if p.is_dir() and p.name.startswith("gdna_"))
    print(f"{'condition':46s} {'ρ_bg':>9} {'off':>8} {'FLOOR':>8} {'Δ':>8}  note")
    rows = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        s = R._setup(inp, index, cfg)
        ra = s["ra"]
        sub = CalibrationSubstrate.from_payload(inp["payload"], ra)
        reg_el = region_eff_length(ra.region_size_bp, np.asarray(inp["gdna_fl_pmf"]))
        bg = measure_background(sub, ra, reg_el)
        rho = 0.0 if not np.isfinite(bg.log_rho_bg) else float(np.exp(bg.log_rho_bg))
        cur = _pass0(s, None)
        flr = _pass0(s, bg)
        d = flr["mwae"] - cur["mwae"]
        g = _group(c)
        note = "ZERO-gDNA FP-check" if g == "gNONE" else (
            "FLAGSHIP" if ("_ss_0.50_" in c and "capture_on" in c) else
            "low-DNA/strong-cap" if (g in ("g1", "g5") and "verystrong" in c) else "")
        rows.append((c, g, cur["mwae"], flr["mwae"], d))
        print(f"{c.replace('gdna_',''):46s} {rho:>9.2e} {cur['mwae']:>8.4f} {flr['mwae']:>8.4f} "
              f"{d:>+8.4f}  {note}", flush=True)

    print("\n=== SUMMARY by group (mean mwae, off → FLOOR) ===")
    for g in ("gNONE", "g1", "g5", "g100", "g300"):
        gr = [r for r in rows if r[1] == g]
        if not gr:
            continue
        mc = float(np.mean([r[2] for r in gr]))
        mf = float(np.mean([r[3] for r in gr]))
        worst = max(gr, key=lambda r: r[4])
        print(f"  {g:6s} n={len(gr):2d}  off={mc:.4f}  FLOOR={mf:.4f}  Δ={mf - mc:+.4f}   "
              f"worst: {worst[0].replace('gdna_', '')} {worst[4]:+.4f}")


if __name__ == "__main__":
    main()
