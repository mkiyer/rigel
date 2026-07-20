"""PROTOTYPE — is a BELIEF-FREE total-density regime a good gDNA-regime stratifier for sigma^2_transfer?

The derivation question (source-only vs pair): the oracle table (transfer_variance_diag.py) proves the gDNA
transfer variance is STRATIFIED under capture (dep-dep 0.33, enr-enr ~1.6 bio, MIXED 10-25) but FLAT off
capture (~0.04 everywhere). A pair-stratified model needs a regime label per node. The only belief-free label
available at APPLY time is the TOTAL density (mass/eff, f_g=1) — but total = gDNA + RNA, so a high-RNA exon
with depleted gDNA is mis-labelled "enriched".

This tool proves the mis-labelling is BENIGN exactly where it happens:
  * under CAPTURE gDNA is enriched at the same probes that make total enriched, so total-regime == gDNA-regime
    (high agreement) — the stratifier is accurate exactly where sigma^2 is stratified and large.
  * OFF capture total-regime and gDNA-regime DIVERGE (RNA contaminates total) BUT sigma^2_transfer,g is flat
    (~0.04) across all regimes, so the mis-label costs nothing.

It reports, per condition:
  (1) agreement between belief-free total-regime and oracle gDNA-regime, on ALL nodes and on exon nodes only;
  (2) the gDNA transfer variance (oracle log-rate, Poisson-removed) on ALL edges, stratified by
      belief-free-total-regime vs oracle-gDNA-regime — do they land in the same place?
  (3) the SOURCE-ONLY blend: sigma^2(regime_src) averaged over that source's edges — showing it conflates
      dep-dep with the crossing.

    OMP_NUM_THREADS=1 python sigma2_transfer_proto.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402

from _disagreement_variance import _adjacent_edges, _poisson_moment_var
from rigel.calibration.bp_solver import node_global_geometry
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import _node_region_type, build_node_geometry
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _pools(p):
    g = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
    r = sum(np.asarray(p[k], float) for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg"))
    return g, r


def _tv(lo, src, dst, mask, n):
    """Poisson-removed transfer variance of log-rate lo over a pair mask."""
    if int(mask.sum()) < 2:
        return np.nan, int(mask.sum())
    return float(_poisson_moment_var(lo[dst][mask] - lo[src][mask], n[src][mask], n[dst][mask])), int(
        mask.sum()
    )


def analyse(inp, index):
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    mass, eff = node_global_geometry(chain, geom)
    e = np.maximum(eff, _EPS)

    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    gR, rR = _pools(inp["region_pools"])
    gB, rB = _pools(inp["boundary_pools"])
    ri, bi = np.clip(idx, 0, gR.shape[0] - 1), np.clip(idx, 0, gB.shape[0] - 1)
    g = np.where(isr, gR[ri], gB[bi])
    r = np.where(isr, rR[ri], rB[bi])

    rtype_node, _ = _node_region_type(chain, ra)
    rtype = np.where(isr, rtype_node, -1)  # 0 interg / 1 intron / 2 exon / -1 boundary
    is_exon = rtype == 2

    rho_g = np.maximum(g, 1.0) / e   # oracle gDNA rate (opportunity floor: 1 fragment)
    rho_t = np.maximum(mass, 1.0) / e  # belief-free TOTAL rate (f_g=1)
    Lg = np.where(g > 0, np.log(rho_g), np.log(1.0 / e))
    Lt = np.where(mass > 0, np.log(rho_t), np.log(1.0 / e))

    src, dst, _ = _adjacent_edges(chain)
    obs = np.asarray(eff, float) > 1e-9 * 1.001
    live = obs[src] & obs[dst]
    src, dst = src[live], dst[live]

    # ---- regime labels (median split; the NPMLE modes are the production version of this) ----
    node_ids = np.unique(np.concatenate([src, dst]))
    reg_or = (Lg > np.median(Lg[node_ids])).astype(int)  # ORACLE gDNA regime (the ceiling)
    reg_bf = (Lt > np.median(Lt[node_ids])).astype(int)  # BELIEF-FREE total regime (implementable)

    # (1) agreement of the belief-free regime with the oracle gDNA regime
    def agree(m):
        return float((reg_or[m] == reg_bf[m]).mean()) if m.any() else np.nan

    all_nodes = np.zeros(len(kind), bool)
    all_nodes[node_ids] = True
    agr_all = agree(all_nodes)
    agr_ex = agree(all_nodes & is_exon)

    # (2) gDNA transfer variance on ALL edges, stratified two ways
    def strata(reg):
        s, d = reg[src], reg[dst]
        return {
            "dep-dep": (s == 0) & (d == 0),
            "enr-enr": (s == 1) & (d == 1),
            "MIXED": s != d,
        }

    rows = []
    for label, reg in (("oracle-gDNA-regime", reg_or), ("belief-free-total-regime", reg_bf)):
        st = strata(reg)
        vals = {k: _tv(Lg, src, dst, m, g)[0] for k, m in st.items()}
        ns = {k: int(m.sum()) for k, m in st.items()}
        rows.append(dict(stratifier=label, **{f"s2g_{k}": vals[k] for k in st}, **{f"n_{k}": ns[k] for k in st}))

    # (3) source-only blend (belief-free regime): variance over all edges from a source in each regime
    so = {}
    for rlab, rv in (("dep-src", 0), ("enr-src", 1)):
        m = reg_bf[src] == rv
        so[rlab] = _tv(Lg, src, dst, m, g)[0]

    return pd.DataFrame(rows), agr_all, agr_ex, so, float(is_exon[node_ids].mean())


def main():
    suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = sorted(p.stem for p in cache.glob("*.pkl"))
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        df, agr_all, agr_ex, so, exfrac = analyse(inp, index)
        cap = "capON " if "capture_on" in c else "capOFF"
        print(f"\n{'=' * 100}\n{c}   [{cap}]")
        print(
            f"  belief-free regime agreement w/ oracle gDNA regime:  ALL nodes={agr_all:.3f}   "
            f"EXON nodes={agr_ex:.3f}   (exon frac of nodes={exfrac:.2f})"
        )
        print(
            df.to_string(
                index=False,
                float_format=lambda x: f"{x:,.3g}",
                columns=[
                    "stratifier", "s2g_dep-dep", "s2g_enr-enr", "s2g_MIXED",
                    "n_dep-dep", "n_enr-enr", "n_MIXED",
                ],
            )
        )
        print(
            f"  SOURCE-ONLY blend (belief-free regime):  dep-src={so['dep-src']:.3g}   "
            f"enr-src={so['enr-src']:.3g}   "
            f"(compare to the pair strata above — source-only blends dep-dep with the crossing)"
        )


if __name__ == "__main__":
    main()
