"""THE CALIBRATION WALK — every stage of the solve scored against the CERTIFIED oracle, in order.

⭐⭐⭐ **WHICH STAGE INTRODUCES THE ERROR?** Owner, 2026-08-18: calibration shows catastrophic error,
so it is to be walked function by function from the accumulator state forward, against a proven ground
truth, until the stage that breaks is NAMED. This instrument is that walk. It runs `calibrate` under a
2x2 of the two big switches — landscape refits x message propagation — and reads the per-slot belief at
every rung the solver already exposes, so the ladder is:

    A  fg_init      the initialisation belief (before any solve)
    B  fg_strand    the strand likelihood alone
    C  fg_loc       the message-free LOCAL solve (strand + density + reference), refits 0
    D  f_g          refits 0,      messages ON     -> C→D is what the RELAY does at pass-0
    E  f_g          shipped refits, messages OFF   -> C→E is what the LANDSCAPE refit does alone
    F  f_g          shipped refits, messages ON    -> the SHIPPED tool

⛔ **TRUTH COMES ONLY FROM THE CERTIFIED TABLE** (`calibration_oracle.py`'s ``slot_truth.npz``) — this
file recomputes nothing about truth and REFUSES to run on a condition whose table is missing, because a
walk against an uncertified truth debugs the wrong thing. The table's ``field_certified`` stamp is
printed with every run: ``true_f_g`` is composition-certified either way; densities are only comparable
when the field gate passed.

⛔ Every arm asserts what ran: ``_uni`` is written only under ``HeadPolicy``, and a muted arm must
reproduce ``f_g == fg_loc`` bit for bit (TRAPS: an-ablation-that-never-ran).

Errors are ``Sum |f_g - true_f_g| * mass`` in FRAGMENTS, total and per stratum, never pooled across
strata in the verdict line. The per-stage DELTA column is the point of the file: the stage whose delta
is large is where to open the code next.

Usage::

    python scripts/design/calibration_walk.py --condition NAME
    python scripts/design/calibration_walk.py --condition NAME --strata "R intron" "R exon"
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


OC = _sibling("object_composition.py")

from rigel.calibration import calibrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX


def run_arm(payload, kw, *, refits: int | None, messages: bool) -> dict:
    """One `calibrate` under one (refits, messages) cell, with the ran/did-not-run assertions."""
    cfg = CalibrationConfig()
    if refits is not None:
        cfg = dataclasses.replace(cfg, calib_refit_iters=refits)
    cfg = dataclasses.replace(cfg, message_propagation=messages)
    debug: dict = {}
    calibrate(payload=payload, config=cfg, _debug=debug,
              **{k: v for k, v in kw.items() if k != "payload"})
    cap = debug["capture"]
    if ("_uni" in cap) != messages:
        raise AssertionError(f"messages={messages} but the relay "
                             f"{'ran' if '_uni' in cap else 'did not run'} — inert or leaking arm")
    if not messages and not np.array_equal(np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])):
        raise AssertionError("muted arm's final belief differs from its local solve — not muted")
    return cap


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--strata", nargs="*", default=None,
                    help="restrict the per-stratum table (default: all seven)")
    args = ap.parse_args()

    truth_path = args.suite / "oracle_cache" / args.condition / "slot_truth.npz"
    if not truth_path.is_file():
        print(f"⛔ no certified truth at {truth_path} — run calibration_oracle.py first. "
              "This instrument does not derive truth, on purpose.")
        return 2
    t = np.load(truth_path, allow_pickle=True)
    true_fg = np.asarray(t["true_f_g"], np.float64)
    mass = (np.asarray(t["n_gdna"], np.float64) + np.asarray(t["n_nrna"], np.float64)
            + np.asarray(t["n_mrna"], np.float64))
    live = np.asarray(t["live"], bool)
    stratum = np.asarray(t["stratum"]).astype(str)
    field_ok = bool(t["field_certified"]) if "field_certified" in t else False

    index = TranscriptIndex.load(args.index)
    cache = read_scan_cache(args.suite / "scan_cache" / args.condition, index)
    kw = calibration_inputs(cache, index)

    shipped_refits = int(CalibrationConfig().calib_refit_iters)
    cell_cd = run_arm(cache.payload, kw, refits=0, messages=False)
    cell_d = run_arm(cache.payload, kw, refits=0, messages=True)
    cell_e = run_arm(cache.payload, kw, refits=None, messages=False)
    cell_f = run_arm(cache.payload, kw, refits=None, messages=True)

    stages = [
        ("A init", np.asarray(cell_cd["fg_init"], np.float64)),
        ("B strand", np.asarray(cell_cd["fg_strand"], np.float64)),
        ("C local (refit0)", np.asarray(cell_cd["fg_loc"], np.float64)),
        ("D C+messages", np.asarray(cell_d["f_g"], np.float64)),
        (f"E C+refits({shipped_refits})", np.asarray(cell_e["f_g"], np.float64)),
        ("F SHIPPED (refits+msgs)", np.asarray(cell_f["f_g"], np.float64)),
    ]

    def err(fg, sel):
        return float((np.abs(fg[sel] - true_fg[sel]) * mass[sel]).sum())

    names = args.strata or sorted(set(stratum[live].tolist()))
    print(f"\n== {args.condition}   truth: {truth_path.name} "
          f"[{'COMPOSITION+FIELD' if field_ok else 'COMPOSITION only'}]   "
          f"live mass {mass[live].sum():,.0f} fragments")
    print("\n   Sum|f_g − true|·mass, FRAGMENTS.  Δ = versus the stage ABOVE it on the ladder "
          "(D and E are both deltas from C; F from E).")
    print(f"   {'stage':<26} {'ALL live':>14} {'Δ':>12}", end="")
    for s in names:
        print(f" {s:>15}", end="")
    print()
    base_for = {"A init": None, "B strand": "A init", "C local (refit0)": "B strand",
                "D C+messages": "C local (refit0)",
                f"E C+refits({shipped_refits})": "C local (refit0)",
                "F SHIPPED (refits+msgs)": f"E C+refits({shipped_refits})"}
    total = {}
    for name, fg in stages:
        total[name] = err(fg, live)
        prev = base_for[name]
        d = "" if prev is None else f"{total[name] - total[prev]:>+12,.0f}"
        print(f"   {name:<26} {total[name]:>14,.0f} {d:>12}", end="")
        for s in names:
            print(f" {err(fg, live & (stratum == s)):>15,.0f}", end="")
        print()
    # ⛔ A and B are INGREDIENT views, not stages C is built on — the strand-only rung cannot see the
    #   structural lock, so B "loses" intergenic and C "recovers" it, and scoring those deltas as
    #   introductions blamed the wrong rung on the first real run. The sequential chain is C→D/E→F.
    seq = [n for n in total if base_for[n] and not n.startswith(("A ", "B ", "C "))]
    worst = max(((n, total[n] - total[base_for[n]]) for n in seq), key=lambda x: x[1])
    print(f"\n   ⭐ largest error INTRODUCED by a SEQUENTIAL stage (C onward): {worst[0]}  "
          f"({worst[1]:+,.0f} fragments)")
    print(f"   ⭐ error already present at C, the message-free local solve: {total['C local (refit0)']:,.0f} "
          "— everything above C inherits what C could not remove")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
