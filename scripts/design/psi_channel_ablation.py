#!/usr/bin/env python
"""⭐⭐⭐ WHICH ψ CHANNEL IS DOING THE WORK? — the ablation, at the ψ boundary, on the SHIPPED solver.

Step 3 of the debug loop, one level below `worst_objects.py`. That instrument says WHICH OBJECTS carry
the error; this says WHICH CHANNEL put it there. It records every argument of the final ψ combine, then
re-solves with one imputed channel nulled at a time — so an attribution is a re-solve of the real call
rather than a second implementation of it.

⛔ **TRAPS: byte-identity-gate — the ``as-is`` arm must reproduce the run BIT-IDENTICALLY, and reproducing it means reproducing
the WRITE-BACK too.** `node_sweep` keeps the incoming belief wherever ``solvable`` is False, so comparing
raw ψ output against the stored belief differs by up to 1.0 at every unsolvable slot. The first version of
this script reported ``max |Δ| = 1.0`` and its numbers were on a different basis than the panel's; the
`wb()` helper is that fix and the identity line is the gate.
⛔ **TRAPS: could-the-arm-have-fired — every ablation prints how many LIVE slots it COULD have moved**, so "no effect" is separable
from "never fired".

⭐ **Read the per-slot table, not only the totals.** TRAPS: all-small-singly-large-jointly: when every single ablation is small
and the joint one is large, the channels are all built from one upstream quantity and ablating consumers
tells you nothing. That is exactly what this prints at `g01 ss0.50 capture_on`: −2.3 / −6.8 / −4.9 / −0.1 %
singly against −60.7 % jointly, while the per-slot rows show the certified-RNA message alone taking eight
of the twelve worst slots from ~0.0086 to ~0.32 against truths of ~0.002.

Usage::

    python scripts/design/psi_channel_ablation.py --condition NAME [--arm pass0|final]
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_DESIGN = Path("/Users/mkiyer/proj/rigel/scripts/design")


def _sibling(name: str):
    key = name[:-3]
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, _DESIGN / name)
    m = importlib.util.module_from_spec(spec)
    sys.modules[key] = m
    spec.loader.exec_module(m)
    return m


P0 = _sibling("pass0_vs_oracle.py")

import rigel.calibration.node_init as NI  # noqa: E402
import rigel.calibration.simplex_logodds as SL  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
CACHE = SUITE / "oracle_cache"
TYPE_NAME = {0: "intergenic", 1: "intron", 2: "exon"}

CALLS: list[dict] = []
_orig = SL._solve_nodes_logodds_all


def _rec(*a, **kw):
    CALLS.append({"args": a, "kw": dict(kw)})
    return _orig(*a, **kw)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="gdna_g01_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--suite", type=Path, default=SUITE)
    ap.add_argument("--index", type=Path, default=INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=CACHE)
    ap.add_argument("--arm", default="final", choices=("pass0", "final"))
    args = ap.parse_args()

    import rigel.calibration.sweep as BP
    patched = []
    for mod in (SL, NI, BP):
        if hasattr(mod, "_solve_nodes_logodds_all"):
            setattr(mod, "_solve_nodes_logodds_all", _rec)
            patched.append(mod.__name__)
    print(f"   [TRAPS: an-ablation-that-never-ran] recording ψ in {patched}")

    index = TranscriptIndex.load(str(args.index))
    m = P0.measure_condition(
        bam=str(args.suite / args.condition / "sim_oracle.bam"), index=index,
        pipeline_config=PipelineConfig(), calibration_config=CalibrationConfig(),
        work_dir=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "eta_dissect",
        tag=args.condition, truth_pmfs=None, oracle_cache=args.oracle_cache,
    )
    print(f"\n⭐⭐ HEAD channel ablation — {args.condition}  arm={args.arm}")
    print(f"   library f_gdna  truth {m.library_f_gdna['T']:.4f}   HEAD {m.library_f_gdna[args.arm]:.4f}")

    dbg = m.debug_final if args.arm == "final" else m.debug_pass0
    chain, cap = dbg["chain"], dbg["capture"]
    ra = dbg["region_arrays"]
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n_slots = kind.shape[0]
    n_nodes, n_edges = int(m.payload.n_nodes), int(m.payload.n_edges)
    node_slot = np.full(n_nodes, -1, np.int64)
    edge_slot = np.full(n_edges, -1, np.int64)
    node_slot[obj[kind == NODE]] = np.flatnonzero(kind == NODE)
    edge_slot[obj[kind != NODE]] = np.flatnonzero(kind != NODE)

    def field(res, which):
        out = np.zeros(n_slots)
        for ax, sl in (("node", node_slot), ("edge", edge_slot)):
            ok = sl >= 0
            out[sl[ok]] = np.asarray(getattr(res, f"mass_{which}_{ax}"), np.float64)[ok]
        return out

    tg, tr = field(m.truth, "gdna"), field(m.truth, "rna")
    pg = field(m.arms[args.arm], "gdna")
    total = tg + tr
    live = total > 0
    err = np.where(live, pg - tg, 0.0)

    combine = [c for c in CALLS if c["kw"].get("lam_imp_prec") is not None]
    C = combine[-1]

    def solve(drop=()):
        kk = dict(C["kw"])
        for n in drop:
            kk[f"{n}_imp_mode"] = None
            kk[f"{n}_imp_prec"] = None
        return np.asarray(_orig(*C["args"], **kk).gdna_frac, np.float64)

    print("   capture keys:", sorted(k for k in cap if isinstance(cap[k], np.ndarray)))
    fg = np.asarray(cap["f_g"], np.float64)
    # ⛔ the WRITE-BACK, reproduced: `node_sweep` keeps the incoming belief wherever `solvable` is False,
    # so comparing raw psi output against the stored belief differs by up to 1.0 at every unsolvable slot.
    cnt = np.asarray(cap["count"], np.float64)
    solvable = (np.asarray(cap["free_pos"], bool) | np.asarray(cap["free_neg"], bool)) & (
        cnt.sum(axis=1) > 0.0)
    def wb(x):
        return np.where(solvable, np.clip(x, 0.0, 1.0), fg)
    base_re = wb(solve())
    d = np.abs(base_re - fg)
    print(f"   ⛔ TRAPS: byte-identity-gate — the re-solve reproduces HEAD's own answer: {np.array_equal(base_re, fg)}"
          f"   (max |Δ| = {d.max():.3e} at slot {int(np.argmax(d))})")

    def sigma(x):
        return float(np.abs(x[live] * total[live] - tg[live]).sum())

    def could(drop):
        n = 0
        for nm in drop:
            p = C["kw"].get(f"{nm}_imp_prec")
            if p is None:
                continue
            for q in (p if isinstance(p, tuple) else (p,)):
                n = max(n, int(((np.asarray(q, np.float64) > 0) & live).sum()))
        return n

    print("\n   ⭐⭐ WHICH CHANNEL IS DOING THE WORK IN HEAD")
    print(f"      {'arm':<32} {'Σ|err| frag':>14} {'vs as-is':>10} {'could move':>11}")
    ref = sigma(wb(solve()))
    for name, drop in (("as-is (HEAD)", ()), ("− gdna_imp (the LEVEL)", ("gdna",)),
                       ("− rna_imp (certified RNA)", ("rna",)), ("− lam_imp (composition)", ("lam",)),
                       ("− theta_imp (tilt)", ("theta",)),
                       ("− ALL messages", ("gdna", "rna", "lam", "theta"))):
        s = sigma(wb(solve(drop)))
        print(f"      {name:<32} {s:>14,.0f} {(s - ref) / max(ref, 1) * 100:>+9.1f}% {could(drop):>11,}")

    # per-channel reach, by slot type
    rt = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    st = np.where(kind == NODE, rt[np.clip(obj, 0, rt.shape[0] - 1)], -1)
    print("\n   ⭐ WHO RECEIVES EACH CHANNEL IN HEAD (nonzero precision), by slot type")
    print(f"      {'channel':<26} {'EDGE':>9} {'intergenic':>11} {'intron':>9} {'exon':>9}")
    chans = {
        "gdna_imp (LEVEL)": np.asarray(C["kw"]["gdna_imp_prec"], np.float64),
        "rna_imp (certified)": np.asarray(C["kw"]["rna_imp_prec"][0], np.float64)
        + np.asarray(C["kw"]["rna_imp_prec"][1], np.float64),
        "lam_imp (composition)": np.asarray(C["kw"]["lam_imp_prec"], np.float64),
        "theta_imp (tilt)": np.asarray(C["kw"]["theta_imp_prec"], np.float64),
    }
    for nm, p in chans.items():
        row = [int(((st == c) & (p > 0)).sum()) for c in (-1, 0, 1, 2)]
        print(f"      {nm:<26} {row[0]:>9,} {row[1]:>11,} {row[2]:>9,} {row[3]:>9,}")
    tot = [int((st == c).sum()) for c in (-1, 0, 1, 2)]
    print(f"      {'(slots of each type)':<26} {tot[0]:>9,} {tot[1]:>11,} {tot[2]:>9,} {tot[3]:>9,}")

    # the top slots, with each of HEAD's channels visible
    order = np.argsort(-np.abs(err))[:12]
    x = {n: wb(solve((n,))) for n in ("gdna", "rna", "lam", "theta")}
    fgl = np.asarray(cap["fg_loc"], np.float64)
    print("\n   ⭐ TOP 12 SLOTS — HEAD's own channels, and what removing each does")
    print(f"   {'#':>3} {'slot':>7} {'type':<11} {'M':>10} {'trueFg':>7} {'own':>7} {'HEAD':>7} "
          f"{'pLevel':>9} {'pRna':>9} {'pLam':>9} {'pThe':>9} "
          f"{'noLvl':>7} {'noRna':>7} {'noLam':>7}")
    print("   " + "-" * 140)
    for r, s in enumerate(order, 1):
        s = int(s)
        print(f"   {r:>3} {s:>7,} {TYPE_NAME.get(int(st[s]), 'EDGE'):<11} {total[s]:>10,.0f} "
              f"{tg[s] / total[s]:>7.4f} {fgl[s]:>7.4f} {fg[s]:>7.4f} "
              f"{chans['gdna_imp (LEVEL)'][s]:>9.3g} {chans['rna_imp (certified)'][s]:>9.3g} "
              f"{chans['lam_imp (composition)'][s]:>9.3g} {chans['theta_imp (tilt)'][s]:>9.3g} "
              f"{x['gdna'][s]:>7.4f} {x['rna'][s]:>7.4f} {x['lam'][s]:>7.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
