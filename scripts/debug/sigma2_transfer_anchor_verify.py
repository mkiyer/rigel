"""Checkpoint 1 + Landmine 1 — does the BELIEF-FREE STRUCTURAL anchor reproduce the strata, and are the
strata stable across coverage depth (or just low-count sampling bias)?

The proof so far used the ORACLE RNA-free mask (r<=1e-9) for the measurement SET. Production must use a
belief-free structural RNA-free proxy. This tool measures sigma^2_transfer,g (oracle log-rate, Poisson-removed)
three ways per regime stratum, and cross-checks coverage-depth stability:

  A. ORACLE anchor      : both endpoints r<=1e-9 (the ceiling; what the proof used)
  B. STRUCTURAL anchor  : both endpoints belief-free RNA-free = intergenic region OR RNA-free single-strand SJ
                          boundary (spliced=0). NO oracle, NO belief.
  C. TOTAL-density (all): every live edge, using TOTAL density as the gDNA substrate (the fallback — valid
                          only where gDNA dominates total, i.e. enriched/crossing under capture).

Landmine-1 check: within each regime stratum, split pairs by coverage depth (mass below/above the median) and
report sigma^2 for each half. If the [0,1.6,25] strata are structural physics, they should be stable across
depth; if they are low-count sampling bias, the low-depth half will be inflated.

    OMP_NUM_THREADS=1 python anchor_verify.py
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
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, nrna_active_strands
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _sj_boundary_mask(bsub, ra):
    """Per-boundary single-strand intron↔exon splice-junction mask (belief-free; inlined from npmle_fusion)."""
    sig = np.asarray(ra.signature).astype(np.int64)
    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    ex_p_l, ex_p_r = (sig_l & BIT_EXON_POS) != 0, (sig_r & BIT_EXON_POS) != 0
    ex_n_l, ex_n_r = (sig_l & BIT_EXON_NEG) != 0, (sig_r & BIT_EXON_NEG) != 0
    sj_pos = (nrp_l & nrp_r) & (ex_p_l ^ ex_p_r)
    sj_neg = (nrn_l & nrn_r) & (ex_n_l ^ ex_n_r)
    return sj_pos ^ sj_neg


def _pools(p):
    g = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
    r = sum(np.asarray(p[k], float) for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg"))
    return g, r


def _tv(lo, src, dst, mask, n):
    if int(mask.sum()) < 2:
        return np.nan, int(mask.sum())
    return float(_poisson_moment_var(lo[dst][mask] - lo[src][mask], n[src][mask], n[dst][mask])), int(mask.sum())


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
    rtype = np.where(isr, rtype_node, -1)

    # belief-free structural RNA-free anchor: intergenic regions + RNA-free single-strand SJ boundaries
    B = np.asarray(bsub.left_region).shape[0]
    bic = np.clip(idx, 0, B - 1)
    spl_b = np.asarray(bsub.left.mass_spliced, float) + np.asarray(bsub.right.mass_spliced, float)
    sj0_b = _sj_boundary_mask(bsub, ra) & (spl_b <= 1e-9)
    struct_rnafree = ((rtype == 0)) | ((~isr) & sj0_b[bic])

    rho_g = np.maximum(g, 1.0) / e
    rho_t = np.maximum(mass, 1.0) / e
    Lg = np.where(g > 0, np.log(rho_g), np.log(1.0 / e))
    Lt = np.where(mass > 0, np.log(rho_t), np.log(1.0 / e))

    src, dst, _ = _adjacent_edges(chain)
    obs = np.asarray(eff, float) > 1e-9 * 1.001
    live = obs[src] & obs[dst]
    src, dst = src[live], dst[live]

    node_ids = np.unique(np.concatenate([src, dst]))
    reg = (Lt > np.median(Lt[node_ids])).astype(int)  # belief-free regime
    s, d = reg[src], reg[dst]
    strat = {"dep-dep": (s == 0) & (d == 0), "enr-enr": (s == 1) & (d == 1), "MIXED": s != d}

    or_free = (r[src] <= 1e-9) & (r[dst] <= 1e-9)
    st_free = struct_rnafree[src] & struct_rnafree[dst]

    # coverage depth: below/above median mass on the edge (min of endpoints)
    mm = np.minimum(mass[src], mass[dst])
    lowcov = mm <= np.median(mm)

    rows = []
    for name, base in strat.items():
        row = {"stratum": name}
        # A oracle anchor (gDNA substrate), B structural anchor (gDNA), C total-density all-edges
        row["A_oracle"], row["nA"] = _tv(Lg, src, dst, base & or_free, g)
        row["B_struct"], row["nB"] = _tv(Lg, src, dst, base & st_free, g)
        row["C_total"], row["nC"] = _tv(Lt, src, dst, base, mass)
        # Landmine-1: coverage-depth stability of the ORACLE-anchor gDNA sigma^2
        row["depth_lo"], _ = _tv(Lg, src, dst, base & or_free & lowcov, g)
        row["depth_hi"], _ = _tv(Lg, src, dst, base & or_free & ~lowcov, g)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    for c in sorted(p.stem for p in cache.glob("*.pkl")):
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        df = analyse(inp, index)
        cap = "capON " if "capture_on" in c else "capOFF"
        nas = "nas" if "nrna_present" in c else "   "
        print(f"\n{'=' * 108}\n{c}   [{cap} {nas}]")
        print(
            df.to_string(
                index=False,
                float_format=lambda x: f"{x:,.3g}",
                columns=["stratum", "A_oracle", "B_struct", "C_total", "nA", "nB", "nC", "depth_lo", "depth_hi"],
            )
        )
    print(
        "\nA_oracle = oracle-RNA-free anchor (ceiling). B_struct = BELIEF-FREE structural anchor. "
        "C_total = total-density all-edges fallback."
    )
    print("depth_lo/hi = oracle-anchor gDNA sigma^2 on the low/high-coverage half (Landmine-1 stability).")


if __name__ == "__main__":
    main()
