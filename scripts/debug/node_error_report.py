"""Per-node calibration error vs the oracle — BEFORE vs AFTER the mature-crossing gate.

A bug-finding lens, not an accuracy metric. For each cached condition (with oracle truth) it runs the real
prior-free ``node_sweep`` TWICE — once with the gate live (`mrna_active_s` mask) and once with the mask forced
all-True (``send_s`` always True ⇒ the gate is a structural no-op ⇒ the un-gated baseline) — and compares each
node's solved ``f_g`` to the oracle's true ``f_g``. No production change: the A/B is a `dataclasses.replace` on
the frozen `NodeStatics`.

It reports (per condition and pooled):
  * the mass-weighted mean |Δf_g| per node class, before vs after (did the gate help / where);
  * the TOP residual nodes after the gate (the worst remaining error — the next bugs to chase);
  * the TOP gate-WORSENED nodes (where the gate made a node worse — a gate-induced regression or a latent bug
    the gate exposed, e.g. strand-blind spliced routing at AMBIG boundaries — D4).

Each row carries the diagnostic context needed to triage: node kind, region class, strand class (AMBIG?),
`mrna_active_{pos,neg}`, `free_{pos,neg}`, the facing RNA eff-length (the short-region fabrication bug), mass.

    OMP_NUM_THREADS=1 python scripts/debug/node_error_report.py [--suite DIR] [--conditions a,b] [--top N]
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from selfsolve_diag import _scan_and_truth, _true_fg

from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.bp_solver import node_sweep
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def _solve(inp, index, cfg, *, gate: bool):
    """Run the prior-free sweep once. ``gate=False`` forces the `mrna_active_*` mask all-True (un-gated)."""
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    st = build_node_statics(chain, sub, bsub, ra)
    if not gate:  # UN-GATED baseline: send_s = (True or not True) = True everywhere
        n = st.n_nodes
        st = dataclasses.replace(
            st,
            mrna_active_pos=np.ones(n, bool),
            mrna_active_neg=np.ones(n, bool),
        )
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    cc = cfg.calibration
    belief = init_beliefs(
        chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window, statics=st,
    )
    final = node_sweep(
        chain, st, geom, belief, ra, rna_sense_frac=kappa, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt, gdna_prior=None,
    )
    return np.asarray(final.f_g, float), (chain, ra, st, geom)


def _node_frame(inp, index, cfg):
    """Per-node truth + diagnostic context (independent of the gate)."""
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    st = build_node_statics(chain, sub, bsub, ra)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    is_r = kind == REGION
    tf_r, m_r = _true_fg(inp["region_pools"])
    tf_b, m_b = _true_fg(inp["boundary_pools"])
    ri = np.clip(idx, 0, tf_r.shape[0] - 1)
    bi = np.clip(idx, 0, tf_b.shape[0] - 1)
    true_fg = np.where(is_r, tf_r[ri], tf_b[bi])
    mass = np.where(is_r, m_r[ri], m_b[bi])
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    sc = np.asarray(ra.strand_class)
    cls = np.where(is_r, np.array([_RTYPE.get(int(t), "?") for t in rtype])[ri], "boundary")
    ambig = np.where(is_r, (sc[ri] == TS_AMBIG), False)  # region AMBIG; boundaries flagged via flanks below
    # a boundary is "AMBIG-adjacent" if either flank region is AMBIG (where D4 strand-routing bites)
    lr = np.asarray(chain.left)
    rr = np.asarray(chain.right)

    def flank_ambig(nbr):
        out = np.zeros(chain.n_nodes, bool)
        ok = nbr >= 0
        nidx = idx[np.clip(nbr, 0, None)]
        reg_nbr = ok & (kind[np.clip(nbr, 0, None)] == REGION)
        out[reg_nbr] = sc[np.clip(nidx[reg_nbr], 0, len(sc) - 1)] == TS_AMBIG
        return out

    b_ambig = (~is_r) & (flank_ambig(lr) | flank_ambig(rr))
    ambig = ambig | b_ambig
    # facing RNA eff-length (min of the two faces — the short-region fabrication lever)
    eff_rna = np.minimum(np.asarray(geom.eff_rna_left), np.asarray(geom.eff_rna_right))
    return dict(
        kind=np.where(is_r, "region", "boundary"),
        cls=cls,
        ambig=ambig,
        true_fg=true_fg,
        mass=mass,
        mrp=np.asarray(st.mrna_active_pos, bool),
        mrn=np.asarray(st.mrna_active_neg, bool),
        fp=np.asarray(st.free_pos, bool),
        fn=np.asarray(st.free_neg, bool),
        eff_rna=eff_rna,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    frames = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, Path("/tmp/rigel_selfsolve"), cache)
        fg_on, _ = _solve(inp, index, cfg, gate=True)
        fg_off, _ = _solve(inp, index, cfg, gate=False)
        d = _node_frame(inp, index, cfg)
        ok = np.isfinite(d["true_fg"]) & (d["mass"] > _EPS)
        df = pd.DataFrame(
            dict(
                cond=c,
                kind=d["kind"],
                cls=d["cls"],
                ambig=d["ambig"],
                mass=d["mass"],
                true_fg=d["true_fg"],
                fg_before=fg_off,
                fg_after=fg_on,
                mrp=d["mrp"],
                mrn=d["mrn"],
                fp=d["fp"],
                fn=d["fn"],
                eff_rna=d["eff_rna"],
            )
        )[ok]
        df["err_before"] = (df["fg_before"] - df["true_fg"]).abs()
        df["err_after"] = (df["fg_after"] - df["true_fg"]).abs()
        df["delta"] = df["err_after"] - df["err_before"]  # <0 = gate improved
        frames.append(df)
    allrows = pd.concat(frames, ignore_index=True)

    def mwae(g):
        w = g["mass"]
        ws = max(w.sum(), _EPS)
        return pd.Series(
            dict(
                n=len(g),
                mass=w.sum(),
                before=(w * g["err_before"]).sum() / ws,
                after=(w * g["err_after"]).sum() / ws,
            )
        )

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    print("\n=== mass-weighted |Δf_g| vs oracle, BEFORE vs AFTER the gate, by class (pooled 7 cond) ===")
    agg = allrows.groupby("cls").apply(mwae, include_groups=False)
    agg["improve"] = agg["before"] - agg["after"]
    print(agg.to_string(float_format=lambda x: f"{x:,.4f}"))
    print("\n--- ALL nodes pooled ---")
    print(mwae(allrows).to_string(float_format=lambda x: f"{x:,.4f}"))

    show = ["cond", "kind", "cls", "ambig", "mass", "true_fg", "fg_before", "fg_after",
            "err_after", "delta", "mrp", "mrn", "eff_rna"]
    print(f"\n=== TOP {a.top} RESIDUAL nodes after the gate (worst remaining |err|, mass-weighted rank) ===")
    allrows["rank"] = allrows["err_after"] * allrows["mass"]
    print(allrows.nlargest(a.top, "rank")[show].to_string(index=False, float_format=lambda x: f"{x:,.3f}"))

    print(f"\n=== TOP {a.top} gate-WORSENED nodes (delta>0, mass-weighted — the bug-hunt: why did the gate hurt?) ===")
    w = allrows[allrows["delta"] > 1e-6].copy()
    w["wrank"] = w["delta"] * w["mass"]
    print(w.nlargest(a.top, "wrank")[show].to_string(index=False, float_format=lambda x: f"{x:,.3f}") if len(w)
          else "(none — the gate worsened no node above threshold)")

    allrows.to_csv("/tmp/node_error_report.tsv", sep="\t", index=False)
    print("\n(full table -> /tmp/node_error_report.tsv)")


if __name__ == "__main__":
    main()
