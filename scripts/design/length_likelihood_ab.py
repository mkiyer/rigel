"""P2 A/B — does the fragment-length likelihood reach the mass the strand channel cannot?

    Plan: `docs/SOLVER_OBSERVABLES_PLAN.md` §6.4   ·   Before-picture: `LEDGER.md`, the P0 entry

⭐ **ONE ARGUMENT DIFFERS BETWEEN THE ARMS.** ``CalibrationConfig.length_likelihood`` False/True; every
other line of code is shared, and the False arm is byte-identical to the P1 path.

⛔ **DO NOT SCORE ON THE ZERO-gDNA ARM ALONE** (`CARRY_FORWARD.md` §3 trap 19). All four `none`
conditions are saturated at 100 % blind, and on a library with no gDNA *any* change that lowers the gDNA
fraction scores better — that one-sidedness has reversed a published verdict in this project once. The
clean comparison is the **gdna100** arm; `gdna100 ss0.50 capture_on` (98.2 % blind, ``f_gdna`` 0.3754
against ~0.50, the row S5.f recorded as "unexplained") is the sharpest single target in the suite.

⚠ **The low-count stratum is scored separately, and it is not a formality.** The likelihood is a
Gaussian on a sum of ``N`` i.i.d. draws, so it is asymptotic in ``N``; measured, the heteroscedastic
``−½ log det`` term displaces the peak by **0.32 at N = 1** and 0.004 at N = 50
(`tests/calibration/test_length_likelihood.py`). The median node's entire evidence is ~0.6 contained
fragments. If this channel does harm, that is where it will show.

    python scripts/design/length_likelihood_ab.py --index IDX --cache-root CACHE_DIR
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

_EPS = 1.0e-9


def _one_arm(inputs, on: bool):
    cfg = dataclasses.replace(CalibrationConfig(), length_likelihood=on)
    dbg: dict = {}
    res = calibrate(**inputs, config=cfg, _debug=dbg)
    cap, chain = dbg["capture"], dbg["chain"]
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    cnt = np.asarray(cap["count"], np.float64).sum(axis=1)
    lock = (~np.asarray(cap["solvable"], bool)) & (np.asarray(chain.kind) == NODE)
    g = float(np.asarray(res.mass_gdna_node).sum() + np.asarray(res.mass_gdna_edge).sum())
    r = float(np.asarray(res.mass_rna_node).sum() + np.asarray(res.mass_rna_edge).sum())
    return {
        "f_gdna": g / (g + r) if (g + r) > 0 else 0.0,
        "no_ev": float(cnt[(tau <= _EPS) & ~lock].sum()) / max(float(cnt.sum()), _EPS),
        "fg_slot": np.asarray(cap["f_g"], np.float64),
        "count": cnt,
    }


def ab_one(index: TranscriptIndex, cache_dir: Path) -> dict:
    cache = read_scan_cache(cache_dir, index)
    inputs = calibration_inputs(cache, index)
    off, on = _one_arm(inputs, False), _one_arm(inputs, True)
    d = np.abs(on["fg_slot"] - off["fg_slot"])
    c = off["count"]
    thin, thick = (c > 0) & (c <= 3), c > 3

    def mass_wt(mask):
        w = c[mask]
        return float((d[mask] * w).sum() / w.sum()) if w.sum() > 0 else 0.0

    return {
        "condition": cache_dir.name,
        "off_fgdna": off["f_gdna"],
        "on_fgdna": on["f_gdna"],
        "off_noev": off["no_ev"],
        "on_noev": on["no_ev"],
        "d_thin": mass_wt(thin),
        "d_thick": mass_wt(thick),
        "thin_mass": float(c[thin].sum() / max(c.sum(), _EPS)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--cache-root", type=Path, required=True)
    ap.add_argument("--conditions", nargs="*", default=None)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    dirs = sorted(d for d in args.cache_root.iterdir() if (d / "payload.npz").exists())
    if args.conditions:
        dirs = [d for d in dirs if d.name in set(args.conditions)]
    rows = [ab_one(index, d) for d in dirs]

    print(
        f"\n{'condition':46s} {'truth':>6s} {'OFF':>8s} {'ON':>8s} {'delta':>8s}   "
        f"{'no-ev OFF':>9s} {'no-ev ON':>9s}"
    )
    print("-" * 46 + " " + "-" * 60)
    for r in rows:
        truth = 0.0 if "gdna_none" in r["condition"] else float("nan")
        t = f"{truth:6.3f}" if truth == truth else "  ~.52"
        print(
            f"{r['condition']:46s} {t} {r['off_fgdna']:8.4f} {r['on_fgdna']:8.4f} "
            f"{r['on_fgdna'] - r['off_fgdna']:+8.4f}   {r['off_noev']:9.1%} {r['on_noev']:9.1%}"
        )

    print(f"\n{'condition':46s} {'thin mass':>10s} {'|d f_g| N<=3':>13s} {'|d f_g| N>3':>12s}")
    print("-" * 46 + " " + "-" * 38)
    for r in rows:
        print(
            f"{r['condition']:46s} {r['thin_mass']:10.1%} {r['d_thin']:13.4f} {r['d_thick']:12.4f}"
        )

    zero = [r for r in rows if "gdna_none" in r["condition"]]
    if zero:
        a = float(np.mean([abs(r["off_fgdna"]) for r in zero]))
        b = float(np.mean([abs(r["on_fgdna"]) for r in zero]))
        print(f"\nzero-gDNA arm mean |f_gdna| (truth 0):  OFF {a:.4f}  ->  ON {b:.4f}")
        print("⛔ one-sided (trap 19) — read it WITH the gdna100 arm, never alone.")


if __name__ == "__main__":
    main()
