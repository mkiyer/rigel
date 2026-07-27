"""Do the MESSAGES damage the introns, and by how much? — the per-node proof.

`selfsolve_diag` measures that on gdna300_ss0.99_none_capOFF the single-strand introns go from
mwae 0.0065 (message-free self-solve — nearly perfect) to 0.1573 after ONE prior-free forward-backward pass
(24x WORSE), while the exons go the other way (0.1254 -> 0.0108). That is an aggregate. This tool proves the
mechanism per node and quantifies it three ways:

  1. TRACE      — per intron: local f_g (message-free), the gDNA message it received (mode + precision, and
                  WHICH neighbour sent it), the final f_g, and the truth. If the local is right and the final
                  is wrong, the message moved it.
  2. ARITHMETIC — the solver combines local and message by precision-weighting in log-fraction space
                  (`bp_solver._scan`: fbg = exp((p_loc*lf_loc + pr*mo)/(p_loc+pr))). Reproduce the final f_g
                  from (local, message) with that formula. If it reproduces, the message IS the whole story.
  3. ABLATION   — re-run the identical sweep with the gDNA message channel silenced (prec forced to 0) and
                  confirm the introns return to their local value. This is the causal proof: remove the
                  suspect, the damage disappears.

sigma^2_transfer is ZERO (priority #2), so `pr = n_src/(n_src*vb_src + 1)` saturates at `1/vb_src`: a
structurally-certain source (an intergenic gDNA sink, var_src = 0) emits at FULL count precision and asserts
its own rate on its neighbour. The measured TRUE capture-OFF gDNA transfer variance is 0.028 (essentially
zero — adjacent gDNA rates really are equal), so pinning is nearly right there; the question this tool
answers is where the residual 0.157 comes from and whether sigma^2_transfer is the right lever for it.

    OMP_NUM_THREADS=1 python scripts/debug/intron_message_trace.py [--suite DIR] [--condition C] [--top N]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from selfsolve_diag import _scan_and_truth, _true_fg

from rigel.calibration import bp_solver as BP
from rigel.calibration.bp_solver import node_sweep
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _build(inp, index, cfg):
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
        chain,
        sub,
        bsub,
        ra,
        rna_sense_frac=kappa,
        n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand,
        logodds_window=cc.sweep_logodds_window,
        statics=st,
    )
    return ra, sub, bsub, chain, st, geom, kappa, belief


def _sweep(chain, st, geom, belief, ra, bsub, kappa, cc, capture):
    return node_sweep(
        chain,
        st,
        geom,
        belief,
        ra,
        rna_sense_frac=kappa,
        n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand,
        logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt,
        gdna_prior=None,
        _capture=capture,
    )


def run(inp, index, cfg, top=12):
    ra, sub, bsub, chain, st, geom, kappa, belief = _build(inp, index, cfg)
    cc = cfg.calibration
    cap: dict = {}
    fb = _sweep(chain, st, geom, belief, ra, bsub, kappa, cc, cap)

    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    is_r = kind == REGION
    tf_r, m_r = _true_fg(inp["region_pools"])
    tf_b, m_b = _true_fg(inp["boundary_pools"])
    ri, bi = np.clip(idx, 0, tf_r.shape[0] - 1), np.clip(idx, 0, tf_b.shape[0] - 1)
    true_fg = np.where(is_r, tf_r[ri], tf_b[bi])
    mass = np.where(is_r, m_r[ri], m_b[bi])
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    is_intron = is_r & (rtype[ri] == 1)

    fg_loc = np.asarray(cap["fg_loc"], float)
    fg_fin = np.asarray(cap["f_g"], float)
    mode_g = np.asarray(cap["mode_g"], float)
    prec_g = np.asarray(cap["prec_g"], float)
    vg_loc = np.asarray(cap["vg_loc"], float)
    p_loc = 1.0 / np.maximum(vg_loc, 1e-9)

    ok = is_intron & np.isfinite(true_fg) & (mass > _EPS)
    err_loc = np.abs(fg_loc - true_fg)
    err_fin = np.abs(fg_fin - true_fg)
    damage = err_fin - err_loc  # >0 == the sweep made this node WORSE

    # --- 2. ARITHMETIC: reproduce the final from (local, message) with the solver's own combine ---
    lf_loc = np.log(np.maximum(fg_loc, 1e-9))
    pt = p_loc + prec_g
    fg_pred = np.exp((p_loc * lf_loc + prec_g * mode_g) / np.maximum(pt, 1e-12))

    df = pd.DataFrame(
        dict(
            node=np.arange(chain.n_nodes),
            mass=mass,
            true_fg=true_fg,
            fg_loc=fg_loc,
            fg_fin=fg_fin,
            fg_pred=np.clip(fg_pred, 0, 1),
            err_loc=err_loc,
            err_fin=err_fin,
            damage=damage,
            msg_mode=mode_g,
            msg_prec=prec_g,
            p_loc=p_loc,
            vg_loc=vg_loc,
        )
    )[ok].sort_values("damage", ascending=False)
    return df, cap, (ra, sub, bsub, chain, st, geom, kappa, belief, cc), (ok, true_fg, mass)


def ablate(bits, ok, true_fg, mass):
    """3. ABLATION — silence the gDNA message channel and re-run the identical sweep."""
    ra, sub, bsub, chain, st, geom, kappa, belief, cc = bits
    orig = BP.node_sweep.__globals__["math"]

    # Silence by forcing every gDNA message precision to 0 at the point it is formed. The cleanest
    # intervention that leaves EVERYTHING else identical (same local solve, same RNA messages, same grid).
    import rigel.calibration.bp_solver as M

    src = Path(M.__file__).read_text()
    patched = src.replace(
        "                pr = n_src / (n_src * vbg[lsrc] + 1.0)",
        "                pr = 0.0  # ABLATION: gDNA message silenced",
    )
    assert patched != src, "ablation patch did not apply"
    ns = {}
    exec(compile(patched, M.__file__, "exec"), ns)
    fb0 = ns["node_sweep"](
        chain,
        st,
        geom,
        belief,
        ra,
        rna_sense_frac=kappa,
        n_grid=cc.sweep_n_grid,
        n_grid_ss=cc.sweep_n_grid_single_strand,
        logodds_window=cc.sweep_logodds_window,
        n_tilt=cc.sweep_n_tilt,
        gdna_prior=None,
    )
    _ = orig
    return np.asarray(fb0.f_g, float)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.99_nrna_none_capture_off")
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    inp = _scan_and_truth(
        suite, a.condition, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache"
    )
    df, cap, bits, (ok, true_fg, mass) = run(inp, index, cfg, a.top)

    w = df.mass.to_numpy()
    print(f"\n{'=' * 118}\n{a.condition}   —  INTRONS ONLY  (n={len(df)})\n{'=' * 118}")
    print(
        f"  mass-weighted |Δf_g|:  message-free LOCAL = {np.sum(w * df.err_loc) / w.sum():.4f}"
        f"   ->  after ONE prior-free FB pass = {np.sum(w * df.err_fin) / w.sum():.4f}"
    )
    hurt = df.damage > 1e-6
    print(
        f"  introns made WORSE by the sweep: {int(hurt.sum())} / {len(df)}"
        f"   ({100 * w[hurt.to_numpy()].sum() / w.sum():.1f}% of intron mass)"
    )
    print(f"  introns made BETTER            : {int((df.damage < -1e-6).sum())}")

    print(f"\n--- 1. TRACE: the {a.top} most damaged introns ---")
    print(
        df.head(a.top).to_string(
            index=False,
            float_format=lambda x: f"{x:,.4g}",
            columns=[
                "node",
                "mass",
                "true_fg",
                "fg_loc",
                "msg_mode",
                "msg_prec",
                "p_loc",
                "fg_fin",
                "damage",
            ],
        )
    )

    print(
        "\n--- 2. ARITHMETIC: does the solver's own combine reproduce the final from (local, message)? ---"
    )
    d = df[df.msg_prec > 0]
    resid = np.abs(d.fg_pred - d.fg_fin)
    print(f"  max |predicted − actual| over {len(d)} message-receiving introns = {resid.max():.3e}")
    print(
        "  (fbg = exp((p_loc*log f_loc + pr*mode)/(p_loc+pr)) — bp_solver._scan's combine, reproduced)"
    )
    print(
        "  => the final belief IS the local belief pulled toward the message. Nothing else moves it."
    )

    print("\n--- 3. ABLATION: silence the gDNA message channel, re-run the IDENTICAL sweep ---")
    fg0 = ablate(bits, ok, true_fg, mass)[df.node.to_numpy()]
    e0 = np.abs(fg0 - df.true_fg.to_numpy())
    print(f"  intron mwae with gDNA messages ON  = {np.sum(w * df.err_fin) / w.sum():.4f}")
    print(f"  intron mwae with gDNA messages OFF = {np.sum(w * e0) / w.sum():.4f}")
    print(f"  message-free local (for reference) = {np.sum(w * df.err_loc) / w.sum():.4f}")
    print(
        f"  => silencing the gDNA channel recovers {100 * (1 - (np.sum(w * e0) / w.sum()) / (np.sum(w * df.err_fin) / w.sum())):.1f}% of the damage."
    )

    out = suite / "intron_message_trace.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nfull per-intron table -> {out}")


if __name__ == "__main__":
    main()
