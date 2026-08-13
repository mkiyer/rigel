#!/usr/bin/env python
"""⭐⭐ **STEP 3 OF THE DEBUG LOOP** — dissect one condition's error down to individual objects.

The loop is: run the panel → measure the full error table → take the worst conditions → **dissect
them to individual regions/boundaries** → find the mechanism → fix → start again.
``pass0_vs_oracle.py`` does the table and ranks CLASSES; this does the last step, which no instrument
covered. A class share says *where* the error lives; it cannot say *why*, and "the relay" is a name
for a set of objects, not a mechanism.

⭐ **RANKED BY ERROR MASS, NEVER BY MEAN ERROR.** A 1 bp region holding two fragments can carry a
``|Δf_g|`` of 1.0 and be worth two fragments of error; an exon holding 40,000 at ``|Δf_g| = 0.05``
carries a thousand times more. Ranking by rate puts the first at the top and sends the reader after
noise. The ranking key here is ``|gDNA_pred − gDNA_true|`` **in fragments**, which is also exactly the
quantity that sums to the ``Σ|err|`` the table reports — so the top of this list is literally the top
of that number.

⭐ **THE FIRST THING TO READ IS THE CONCENTRATION CURVE**, before any individual row. If the top
hundred objects hold most of the error, there is a mechanism to find and a handful of objects that
demonstrate it. If the error is spread evenly over tens of thousands, there is a systematic bias and
chasing individual objects is wasted effort. Those two findings call for completely different work,
and the curve distinguishes them in one boundary.

⚠ **EVERY COLUMN IS AN INPUT THE SOLVER ACTUALLY SAW, OR AN OUTPUT IT ACTUALLY PRODUCED.** Nothing
here is re-derived: the counts come from the substrate, the divisors from the geometry the solver
divided by, the classes from ``pass0_vs_oracle`` (one definition, not a second copy), and the belief
from ``_debug["capture"]``. A dissection tool that recomputes its own version of the solver's inputs
is debugging a different program.

⚠ **THE NEIGHBOUR COLUMNS ARE THE POINT ON RELAY-ONLY OBJECTS.** An object with no own evidence takes
its answer from its neighbours, so its error is only interpretable beside theirs. For a REGION the
neighbours are the two flanking BOUNDARY slots and vice versa — the chain is ``N E N E … N`` per
reference, so "neighbour" is unambiguous and needs no graph traversal.

⛔ Undrained, like the instrument it builds on, and for the same forced reason (see that module).

Usage::

    python scripts/design/worst_objects.py --condition NAME [--suite DIR] [--top 40]
    python scripts/design/worst_objects.py --condition NAME --arm final --axis region
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    """Load a sibling instrument by path. ⚠ Registered in ``sys.modules`` before execution, or its
    dataclasses fail to resolve their own module."""
    key = name[:-3]
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


P0 = _sibling("pass0_vs_oracle.py")

from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.types import Strand  # noqa: E402

TYPE_NAME = {0: "intergenic", 1: "intron", 2: "exon"}
#: ⚠ ``region_arrays.strand_class`` is the region's TRANSCRIPT-strand class over the whole ``Strand``
#: enum, AMBIGUOUS included — and AMBIGUOUS is the interesting one here, because it is where the
#: strand λ-term is gated off by the Schur argument. Naming all four rather than two keeps that
#: visible in the table instead of collapsing it into "other".
STRAND_NAME = {
    int(Strand.NONE): "none",
    int(Strand.POS): "+",
    int(Strand.NEG): "-",
    int(Strand.AMBIGUOUS): "amb",
}

#: Where the concentration curve is sampled. ⚠ Not thresholds — nothing branches on these. They are
#: read points on a cumulative curve, chosen to span three orders of magnitude so the SHAPE is
#: visible rather than a single number that could be read either way.
CUTS = (10, 100, 1_000, 10_000)


def concentration(err: np.ndarray) -> list[tuple[int, float, float]]:
    """``(k, share of Σ|err| in the top k objects, share of objects that is)`` at each cut."""
    order = np.argsort(-np.abs(err))
    ranked = np.abs(err)[order]
    total = ranked.sum()
    n = ranked.shape[0]
    out = []
    for k in CUTS:
        if k > n:
            break
        out.append((k, float(ranked[:k].sum() / total) if total > 0 else 0.0, k / n))
    return out


def _slot_lookup(chain, n_regions: int, n_boundaries: int):
    """``(region_slot, boundary_slot)`` — the chain slot each object occupies. The chain is ``N E N E … N``
    per reference, so this is a bijection per axis and there is nothing to pool."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    region_slot = np.full(n_regions, -1, np.int64)
    boundary_slot = np.full(n_boundaries, -1, np.int64)
    region_slot[obj[kind == REGION]] = np.flatnonzero(kind == REGION)
    boundary_slot[obj[kind == BOUNDARY]] = np.flatnonzero(kind == BOUNDARY)
    return region_slot, boundary_slot


def dissect(m, *, axis: str, arm: str, top: int, index) -> dict:
    """Rank one axis of one arm by error MASS and assemble the per-object diagnostic table."""
    truth, res = m.truth, m.arms[arm]
    pred_g = np.asarray(getattr(res, f"mass_gdna_{axis}"), np.float64)
    true_g = np.asarray(getattr(truth, f"mass_gdna_{axis}"), np.float64)
    true_r = np.asarray(getattr(truth, f"mass_rna_{axis}"), np.float64)
    err = pred_g - true_g
    total = true_g + true_r
    live = total > 0

    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    ra = m.debug_pass0["region_arrays"]
    n_regions, n_boundaries = int(m.payload.n_regions), int(m.payload.n_boundaries)
    region_slot, boundary_slot = _slot_lookup(chain, n_regions, n_boundaries)
    slot_of = region_slot if axis == "region" else boundary_slot

    solver, info = m.solver_masks[axis], m.info_masks[axis]
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    var_g = np.asarray(cap["var_g"], np.float64)
    fg_loc = np.asarray(cap["fg_loc"], np.float64)
    counts = np.asarray(cap["count"], np.float64)
    eff_g = np.asarray(cap["eff_gdna"], np.float64)
    eff_r = np.asarray(cap["eff_rna"], np.float64)

    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    strand_class = np.asarray(ra.strand_class, np.int64)
    ref_id, start, end = (np.asarray(ra.ref_id), np.asarray(ra.start), np.asarray(ra.end))
    ref_name = {v: k for k, v in index.ref_name_to_id.items()}

    # per-slot error, so a neighbour's error can be read off directly
    slot_err = np.zeros(chain.kind.shape[0], np.float64)
    for ax, sl in (("region", region_slot), ("boundary", boundary_slot)):
        tg = np.asarray(getattr(truth, f"mass_gdna_{ax}"), np.float64)
        pg = np.asarray(getattr(res, f"mass_gdna_{ax}"), np.float64)
        ok = sl >= 0
        slot_err[sl[ok]] = (pg - tg)[ok]

    order = np.argsort(-np.abs(np.where(live, err, 0.0)))[:top]
    rows = []
    for obj in order:
        if not live[obj]:
            continue
        s = int(slot_of[obj])
        cls = next((c for c in P0.SOLVER_CLASSES if solver[c][obj]), "?")
        icls = next((c for c in P0.INFO_CLASSES if info[c][obj]), "?")
        if axis == "region":
            where = f"{ref_name.get(int(ref_id[obj]), '?')}:{int(start[obj])}-{int(end[obj])}"
            kind_txt = TYPE_NAME.get(int(rtype[obj]), "?")
            sc = STRAND_NAME.get(int(strand_class[obj]), "?")
            span = int(end[obj] - start[obj])
        else:  # a contiguous boundary is a 0-bp boundary; name it by the slot's neighbours
            where, kind_txt, sc, span = f"boundary#{obj}", "boundary", "-", 0
        rows.append(
            {
                "obj": int(obj), "where": where, "type": kind_txt, "strand": sc, "bp": span,
                "n_pos": float(counts[s, 0]), "n_neg": float(counts[s, 1]),
                "eff_g": float(eff_g[s]), "eff_r": float(eff_r[s]),
                "true_g": float(true_g[obj]), "true_r": float(true_r[obj]),
                "pred_fg": float(pred_g[obj] / total[obj]), "true_fg": float(true_g[obj] / total[obj]),
                "err": float(err[obj]), "cls": cls, "info": icls,
                "tau": float(tau[s]), "var_g": float(var_g[s]), "fg_loc": float(fg_loc[s]),
                # ⭐ the region's own density in the gDNA frame -- the quantity the density model
                # compares against the intergenic background. N / E_g, nothing re-derived.
                "rho": float((counts[s, 0] + counts[s, 1]) / max(eff_g[s], 1e-12)),
                "nb_err": [float(slot_err[s - 1]) if s > 0 else float("nan"),
                           float(slot_err[s + 1]) if s + 1 < slot_err.shape[0] else float("nan")],
            }
        )
    return {
        "rows": rows,
        "concentration": concentration(np.where(live, err, 0.0)),
        "n_live": int(live.sum()),
        "abs_err": float(np.abs(err[live]).sum()),
        "solver": solver,
        "info": info,
        "err": err,
        "live": live,
        "total": total,
    }


def _profile(d: dict, masks: dict, names) -> list[tuple[str, float, float]]:
    """``(class, share of the TOP objects' error, share of ALL error)`` — what the worst objects have
    in common, against the background rate. ⚠ Both shares, always: a class that is 60 % of the worst
    objects and 60 % of everything is not a finding."""
    err, live = np.abs(d["err"]) * d["live"], d["live"]
    top = np.zeros_like(live)
    top[[r["obj"] for r in d["rows"]]] = True
    t_tot, a_tot = err[top].sum(), err.sum()
    out = []
    for name in names:
        m = masks[name] & live
        out.append(
            (
                name,
                float(err[m & top].sum() / t_tot) if t_tot > 0 else 0.0,
                float(err[m].sum() / a_tot) if a_tot > 0 else 0.0,
            )
        )
    return out


def report(m, d: dict, axis: str, arm: str, top: int) -> None:
    print()
    print("=" * 120)
    print(f"⭐⭐ WORST OBJECTS — {m.condition}   axis={axis}  arm={arm}")
    print("=" * 120)
    print(f"   {d['n_live']:,} objects carry mass; Σ|err| = {d['abs_err']:,.0f} fragments")
    print()
    print("   ⭐ CONCENTRATION — read this FIRST. Concentrated ⇒ there is a mechanism and a few")
    print("      objects demonstrate it. Diffuse ⇒ a systematic bias, and individual rows are noise.")
    for k, share, obj_share in d["concentration"]:
        print(f"      top {k:>6,} objects ({obj_share:>6.2%} of them) hold {share:>6.2%} of Σ|err|")

    print()
    print("   ⭐ WHAT THE WORST OBJECTS HAVE IN COMMON (share of the top's error vs of ALL error —")
    print("      a class that matches the background is not a finding)")
    for label, masks, names in (
        ("solver class", d["solver"], P0.SOLVER_CLASSES),
        ("C_info class", d["info"], P0.INFO_CLASSES),
    ):
        print(f"      {label}")
        for name, t, a in _profile(d, masks, names):
            flag = " ⭐" if a > 0 and t / max(a, 1e-12) > 1.5 else ""
            print(f"        {name:<22} top {t:>7.1%}   all {a:>7.1%}   ratio "
                  f"{t / a if a > 0 else float('nan'):>5.1f}x{flag}")

    print()
    print(f"   ⭐ THE TOP {len(d['rows'])} BY ERROR MASS")
    print(f"   {'#':>3} {'where':<26} {'type':<10} {'st':<3} {'n+':>8} {'n-':>8} "
          f"{'eff_g':>8} {'rho':>7} {'true_fg':>8} {'fg_loc':>7} {'pred_fg':>8} {'err':>10} "
          f"{'class':<13} {'C_info':<20} {'nb err':>18}")
    print("   " + "-" * 165)
    for i, r in enumerate(d["rows"], 1):
        nb = "/".join("—" if np.isnan(x) else f"{x:+,.0f}" for x in r["nb_err"])
        print(
            f"   {i:>3} {r['where']:<26} {r['type']:<10} {r['strand']:<3} {r['n_pos']:>8,.0f} "
            f"{r['n_neg']:>8,.0f} {r['eff_g']:>8.1f} {r['rho']:>7.2f} {r['true_fg']:>8.3f} "
            f"{r['fg_loc']:>7.3f} {r['pred_fg']:>8.3f} {r['err']:>+10,.0f} {r['cls']:<13} "
            f"{r['info']:<20} {nb:>18}"
        )
    print()
    print("   fg_loc = the MESSAGE-FREE local self-solve; pred_fg = after the sweep. ⭐ If the two")
    print("   agree, the messages are innocent and the local solve is the defect — and vice versa.")
    print("   rho = N / E_g, the region's own density in the gDNA frame (vs the intergenic background).")
    print("   nb err = the two ADJACENT chain slots' gDNA mass errors. ⭐ On a relay-only object this")
    print("   is the whole explanation: it has no evidence of its own, so its answer came from these.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--arm", default="pass0", choices=("pass0", "final"))
    ap.add_argument("--axis", default="region", choices=("region", "boundary", "both"))
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")))
    ap.add_argument("--oracle-cache", type=Path, default=None)
    args = ap.parse_args()

    cond_dir = args.suite / args.condition
    bam = str(cond_dir / "sim_oracle.bam")
    if not Path(bam).is_file():
        print(f"no BAM at {bam}", file=sys.stderr)
        return 2
    index = TranscriptIndex.load(str(args.index))
    m = P0.measure_condition(
        bam=bam,
        index=index,
        pipeline_config=PipelineConfig(),
        calibration_config=CalibrationConfig(),
        work_dir=args.work_dir / "rigel_pass0_oracle",
        tag=args.condition,
        truth_pmfs=lambda size, d=cond_dir: (
            P0.truth_length_pmf(d, "gdna", size),
            P0.truth_length_pmf(d, "rna", size),
        ),
        oracle_cache=args.oracle_cache,
    )
    for axis in (("region", "boundary") if args.axis == "both" else (args.axis,)):
        report(m, dissect(m, axis=axis, arm=args.arm, top=args.top, index=index), axis, args.arm,
               args.top)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
