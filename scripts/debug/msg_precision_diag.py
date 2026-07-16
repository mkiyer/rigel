"""Is σ²_imp (the message precision) inflated by RNA? — the unstranded-peeling diagnostic.

`bp_solver.adjacent_disagreement_variance` fits σ²_imp on the **total** density (`mass/eff_gdna`) of adjacent
node pairs — "the naive total-density frame, RNA included" (its own docstring). But σ²_imp is then used as the
reliability of the **gDNA** (and per-strand RNA) density messages: `pr = n_src/(n_src·σ²_imp + 1)`, which
saturates at `1/σ²_imp` no matter how confident the source is. Adjacent TOTAL densities disagree wildly (an
expressed exon beside an intron differs by orders of magnitude) whereas true gDNA density is genomically
smooth — so σ²_imp is over-estimated, the precision cap collapses, and messages become impotent. That is the
suspected root of the unstranded gDNA over-call (RNA can't be peeled without a strand OR a strong message).

This measures, per condition:
  * σ²_imp as PRODUCTION computes it (all adjacent pairs, total density), and its precision cap 1/σ²_imp;
  * σ²_imp on gDNA-CLEAN structural pairs (both endpoints intergenic/intron, where total ≈ gDNA) — the honest
    gDNA imputation spread — and its cap.
A large ratio ⇒ production needlessly throttles every message.

    OMP_NUM_THREADS=1 python scripts/debug/msg_precision_diag.py --suite DIR --cache-dir DIR
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rigel.calibration.bp_solver import (
    _adjacent_edges,
    _poisson_moment_var,
    adjacent_disagreement_variance,
    build_node_geometry,
)
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1.0e-9


def _gdna_clean_nodes(chain, region_arrays, bsub):
    """Per chain node: is its mass gDNA-clean BY STRUCTURE? A REGION iff intergenic/intron (rtype 0/1); a
    BOUNDARY iff BOTH flanking regions are intergenic/intron (no exon RNA on either side)."""
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    ref = np.asarray(chain.ref_idx, np.int64)
    node_rtype, rtype = _node_region_type(chain, region_arrays)
    clean_reg = is_reg & ((node_rtype == 0) | (node_rtype == 1))

    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    B, R = lr.shape[0], rtype.shape[0]
    bi = np.clip(ref, 0, B - 1)
    lt = np.where(lr[bi] >= 0, rtype[np.clip(lr[bi], 0, R - 1)], -1)
    rt = np.where(rr[bi] >= 0, rtype[np.clip(rr[bi], 0, R - 1)], -1)
    clean_bnd = (~is_reg) & np.isin(lt, [0, 1]) & np.isin(rt, [0, 1])
    return clean_reg | clean_bnd


def _sigma2_on_subset(chain, geometry, keep_node):
    """σ²_imp restricted to adjacent pairs whose BOTH endpoints are in ``keep_node`` (same estimator)."""
    ML = np.asarray(geometry.mass_left)
    MR = np.asarray(geometry.mass_right)
    EGL = np.asarray(geometry.eff_gdna_left)
    EGR = np.asarray(geometry.eff_gdna_right)
    src, dst, s_bnd = _adjacent_edges(chain)
    both = keep_node[src] & keep_node[dst]
    src, dst, s_bnd = src[both], dst[both], s_bnd[both]
    n_i, e_i, n_j, e_j = MR[src], EGR[src], ML[dst], EGL[dst]
    ok = (n_i > _EPS) & (n_j > _EPS) & (e_i > _EPS) & (e_j > _EPS)
    n_i, e_i, n_j, e_j, s_bnd = n_i[ok], e_i[ok], n_j[ok], e_j[ok], s_bnd[ok]
    if n_i.size < 2:
        return float("nan"), 0
    lr_i, lr_j = np.log(n_i / e_i), np.log(n_j / e_j)
    resid = np.where(s_bnd, lr_i - lr_j, lr_j - lr_i)
    return _poisson_moment_var(resid, n_i, n_j), int(n_i.size)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--cache-dir", required=True, help="calib_pool_benchmark --cache-dir (has payloads)")
    ap.add_argument("--conditions", default=None)
    args = ap.parse_args()
    suite = Path(args.suite)
    cache_dir = Path(args.cache_dir)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cfg = PipelineConfig()  # noqa: F841 (kept for parity with the benchmark wiring)

    conds = (
        args.conditions.split(",")
        if args.conditions
        else sorted(p.stem for p in cache_dir.glob("*.pkl"))
    )
    print(f"{'condition':46s} {'σ²_prod':>9s} {'cap':>7s} {'σ²_clean':>9s} {'cap':>8s} {'ratio':>7s} {'n_pairs':>8s}")
    for cond in conds:
        p = cache_dir / f"{cond}.pkl"
        if not p.exists():
            continue
        with open(p, "rb") as fh:
            inp = pickle.load(fh)
        payload = inp["payload"]
        sub = CalibrationSubstrate.from_payload(payload, ra)
        bsub = BoundarySubstrate.from_payload(payload)
        chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
        geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])

        s2_prod = adjacent_disagreement_variance(chain, geom)
        clean = _gdna_clean_nodes(chain, ra, bsub)
        s2_clean, npairs = _sigma2_on_subset(chain, geom, clean)
        cap_p = 1.0 / max(s2_prod, _EPS)
        cap_c = 1.0 / max(s2_clean, _EPS) if np.isfinite(s2_clean) else float("nan")
        ratio = s2_prod / s2_clean if (np.isfinite(s2_clean) and s2_clean > 0) else float("nan")
        print(
            f"{cond.replace('gdna_','').replace('_capture','_cap'):46s} {s2_prod:9.3f} {cap_p:7.2f} "
            f"{s2_clean:9.3f} {cap_c:8.2f} {ratio:7.2f} {npairs:8d}"
        )


if __name__ == "__main__":
    main()
