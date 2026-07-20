"""FULL-SOLVE A/B of the V4 fix: keep the NPMLE gDNA arm but ADD BACK the ½log(f_g) Jeffreys barrier (so the
reference measure is symmetric: ½log f_g + ½log(1−f_g)), instead of the NPMLE REPLACING the gDNA arm.

Monkeypatches `simplex_logodds._gdna_arm` (no production edit — safe vs the concurrent editor) and runs all 7
cached conditions at pass-0, comparing to CURRENT (unpatched). The fix should LIFT the flagship (unstranded
capON) without REGRESSING the other six (the barrier is dormant except near f_g→0).

    OMP_NUM_THREADS=1 python test_v4_symmetric.py
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
import rigel.calibration.simplex_logodds as slo  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_orig_gdna_arm = slo._gdna_arm


def _v4_gdna_arm(lam, global_logprior):
    """V4: NPMLE population shape + the ½log(f_g) Jeffreys boundary barrier (restore the symmetric measure)."""
    base = _orig_gdna_arm(lam, global_logprior)
    if global_logprior is not None:
        base = np.asarray(base, np.float64) + slo._JEFFREYS_REF * slo._log_fg(lam)[None, :]
    return base


def _pass0(s, bw, patched):
    slo._gdna_arm = _v4_gdna_arm if patched else _orig_gdna_arm
    try:
        prior = DensityNPMLE.fit(s["mass_g"], s["eff_g"], bandwidth=bw)
        belief = R._sweep(s, s["b0"], prior, transfer_variance=True)
        return R._measure(s, belief, prior, None)
    finally:
        slo._gdna_arm = _orig_gdna_arm


def main():
    suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    # ALL 24 scenarios (gdna_none = the zero-gDNA false-positive control; gdna100/gdna300 = moderate/heavy).
    conds = sorted(p.name for p in suite.iterdir() if p.is_dir() and p.name.startswith("gdna_"))
    print(f"{'condition':44s} {'current':>8} {'V4':>8} {'Δ':>8}   note")
    rows = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        s = R._setup(inp, index, cfg)
        cur = _pass0(s, 0.15, patched=False)
        v4 = _pass0(s, 0.15, patched=True)
        d = v4["mwae"] - cur["mwae"]
        gd = "gNONE" if "gdna_none" in c else ("g100" if "gdna100" in c else "g300")
        # for gdna_none, oracle f_g = 0 ⇒ mwae IS the false-positive gDNA (want ~0, must NOT rise under V4)
        note = "ZERO-gDNA FP-check" if gd == "gNONE" else (
            "FLAGSHIP" if ("0.50" in c and "capture_on" in c) else "")
        rows.append((c, gd, cur["mwae"], v4["mwae"], d))
        print(f"{c.replace('gdna_',''):44s} {cur['mwae']:>8.4f} {v4['mwae']:>8.4f} {d:>+8.4f}   {note}",
              flush=True)

    print("\n=== SUMMARY by group (mean mwae) ===")
    for g in ("gNONE", "g100", "g300"):
        gr = [r for r in rows if r[1] == g]
        mc = float(np.mean([r[2] for r in gr]))
        mv = float(np.mean([r[3] for r in gr]))
        worst = max(gr, key=lambda r: r[4])
        print(f"  {g:6s} n={len(gr):2d}  current={mc:.4f}  V4={mv:.4f}  Δ={mv-mc:+.4f}   "
              f"worst-regression: {worst[0].replace('gdna_','')} {worst[4]:+.4f}")
    allr = rows
    print(f"  ALL    n={len(allr):2d}  current={np.mean([r[2] for r in allr]):.4f}  "
          f"V4={np.mean([r[3] for r in allr]):.4f}  "
          f"Δ={np.mean([r[3] for r in allr])-np.mean([r[2] for r in allr]):+.4f}")


if __name__ == "__main__":
    main()
