"""Channel-ablation replay — the message-layer A/B tripwire for the calibration sweep.

Runs the real prior-free ``node_sweep`` ONCE per condition (capturing its phase-D inputs via the ``_capture``
hook), then re-solves phase D with the gDNA / RNA imputation channels individually silenced (their message
precision forced to 0). Reports the mass-weighted ``|f_g − true_f_g|`` per (node class × condition) for four
arms: ``msgfree`` / ``gdna_only`` / ``rna_only`` / ``both``.

Two properties make it the trustworthy A/B basis (both verified by independent reviewers):
  * **Fidelity.** ``both`` reuses the captured ``(mode, prec, global_lp, solvable)`` and re-runs only the
    final solve, so ``max|both − shipped f_g| = 0.000e+00`` on every condition — printed each run as the
    tripwire. A non-zero value means the snapshot drifted (a concurrent edit) — discard the run.
  * **Basis.** ``msgfree`` here is phase-D with all message precisions zeroed — NOT ``selfsolve_diag --stage
    self`` (which runs ``init_beliefs`` alone). The two differ; this one (``0.0098`` intron on
    ``gdna300_ss0.99_none_capOFF``) is the message-free floor the gate is judged against, reported as
    ``msgfree``. See ``docs/calibration/archive/mature_crossing_gate.md`` §4.4.

Truth basis: the oracle region/boundary pools in ``<suite>/_selfsolve_cache`` (7 conditions, gdna300, with
truth). Reproduces the item-1 doc's §1 4-row table.

    OMP_NUM_THREADS=1 python scripts/debug/ablate_replay.py [--suite DIR] [--conditions a,b]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from selfsolve_diag import _scan_and_truth, _true_fg

from rigel.calibration.bp_solver import node_sweep
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def build(inp, index, cfg):
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    st = build_node_statics(chain, sub, bsub, ra)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    cc = cfg.calibration
    belief = init_beliefs(
        chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window, statics=st,
    )
    cap: dict = {}
    node_sweep(
        chain, st, geom, belief, ra, rna_sense_frac=kappa, n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt, gdna_prior=None, _capture=cap,
    )
    return ra, sub, bsub, chain, st, geom, kappa, belief, cap, cc


def replay(st, kappa, cc, belief, cap, use_g: bool, use_r: bool):
    """Re-run phase-D final solve with channels ablated (prec forced 0). Everything else identical."""
    z = np.zeros_like(np.asarray(cap["prec_g"], float))
    gm = np.asarray(cap["mode_g"], float) if use_g else z
    gp = np.asarray(cap["prec_g"], float) if use_g else z
    rm = (np.asarray(cap["mode_p"], float), np.asarray(cap["mode_n"], float)) if use_r else (z, z)
    rp = (np.asarray(cap["prec_p"], float), np.asarray(cap["prec_n"], float)) if use_r else (z, z)
    dc = _solve_nodes_logodds_all(
        st.u_pos, st.u_neg, st.free_pos, st.free_neg, st.mass_unspliced, st.mass_spliced,
        kappa=kappa, od_g=0.0, od_r=0.0, n_grid=int(cc.sweep_n_grid),
        L=float(cc.sweep_logodds_window), n_tilt=cc.sweep_n_tilt,
        n_grid_ss=cc.sweep_n_grid_single_strand, global_logprior=cap["global_lp"],
        gdna_imp_mode=gm, gdna_imp_prec=gp, rna_imp_mode=rm, rna_imp_prec=rp,
        fg_ref=np.asarray(belief.f_g, float), fpos_ref=np.asarray(belief.f_pos, float),
        fneg_ref=np.asarray(belief.f_neg, float),
    )
    solvable = np.asarray(cap["solvable"], bool)
    return np.where(solvable, np.clip(dc.gdna_frac, 0.0, 1.0), np.asarray(belief.f_g, float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    rows = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, Path("/tmp/rigel_selfsolve"), cache)
        ra, sub, bsub, chain, st, geom, kappa, belief, cap, cc = build(inp, index, cfg)
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
        cls = np.where(is_r, np.array([_RTYPE.get(int(t), "?") for t in rtype])[ri], "boundary")

        variants = {
            "msgfree": replay(st, kappa, cc, belief, cap, False, False),
            "gdna_only": replay(st, kappa, cc, belief, cap, True, False),
            "rna_only": replay(st, kappa, cc, belief, cap, False, True),
            "both": replay(st, kappa, cc, belief, cap, True, True),
        }
        shipped = np.asarray(cap["f_g"], float)
        ok = np.isfinite(true_fg) & (mass > _EPS)
        print(f"\n[{c}] kappa={kappa:.4f}  max|both - shipped| = "
              f"{np.max(np.abs(variants['both'] - shipped)):.3e}  (replay fidelity check)")
        for cl in ["intron", "exon", "intergenic", "boundary", "ALL"]:
            m = ok if cl == "ALL" else (ok & (cls == cl))
            if m.sum() == 0:
                continue
            w = mass[m]
            ws = max(w.sum(), _EPS)
            r = dict(condition=c, cls=cl, n=int(m.sum()), mass=float(w.sum()),
                     true_fg=float((w * true_fg[m]).sum() / ws))
            for k, v in variants.items():
                r[k] = float((w * np.abs(v[m] - true_fg[m])).sum() / ws)
            rows.append(r)
    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print("\n" + "=" * 130)
    print(df.to_string(index=False, float_format=lambda x: f"{x:,.4f}"))
    df.to_csv("/tmp/ablate_matrix.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
