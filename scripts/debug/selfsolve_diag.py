"""Steps 3 & 5 — the MESSAGE-FREE SELF-SOLVE and the PRIOR-FREE FORWARD-BACKWARD SOLVE, vs the oracle.

The question this answers (calibration roadmap step 2/3): **before any message propagation, which nodes can
solve themselves, and with what precision?**

    "Boundaries should solve and have precision. Strand-specific data should solve."

The message-free self-solve IS ``node_geometry.init_beliefs`` — the bare per-node ψ (strand likelihood + the
derived reference; NO prior, NO messages). This tool runs it on cached real payloads and reports, per
(node kind × region class), the solved ``f_g`` against the validated oracle truth, plus the two quantities
that decide whether a node can TEACH its neighbours:

  * ``var_g``  = ``Var(log f_g)`` — the node's own belief width. The no-information value is the reference's
    own ``Var(log f_g)`` (Beta(½,½) truncated to the λ window, ≈2.8 at L=10): a node at ~2.8 learned NOTHING.
  * ``pr``     = ``n/(n·var_g + 1)`` = ``1/(var_g + 1/n)`` — the precision its message would carry
    (``bp_solver._scan``). Saturates at ``1/var_g``, so the reference-only ceiling is ``1/2.8 ≈ 0.357``
    pseudo-observations. **pr >> 0.357 means the node genuinely knows something.**

Truth basis: :meth:`oracle.OracleTruth.region_pools` / :meth:`~oracle.OracleTruth.boundary_pools` — the
production accumulator split by true fragment origin (sum-to-full validated). Region nodes use contained
mass; boundary nodes use ``left + right`` summed, matching ``_boundary_strand_stats``' pie base.

`--stage self` (step 3) runs `init_beliefs` alone: the bare per-node psi (strand + the derived reference; NO
prior, NO messages). `--stage fb` (step 5) additionally runs the single forward-backward `node_sweep` with
`gdna_prior=None` — prior-free but NOT reference-free — so the two columns isolate exactly what MESSAGES add.

sigma^2_transfer is currently ZERO (priority #2), so a structurally-certain source (var_src=0, e.g. an
intergenic gDNA sink) emits at FULL count precision and PINS its neighbour rather than informing it. The
`fb` stage is where that shows up, and measuring its cost is the point: the sigma^2_transfer choice is meant
to come from this diagnostic rather than be derived a priori.

    OMP_NUM_THREADS=1 python scripts/debug/selfsolve_diag.py [--suite DIR] [--conditions a,b]
                                                             [--stage self|fb|both]

The cache (``<suite>/_selfsolve_cache``) holds the scan + BOTH truth pools, so re-running after a solver
change costs only init_beliefs + metrics. It is separate from ``_calib_cache`` because that one predates
``boundary_pools`` and stores region truth only.
"""

from __future__ import annotations

import argparse
import os
import pickle
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from oracle import OracleTruth

from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.bp_solver import node_sweep
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import (
    build_node_geometry,
    build_node_statics,
    init_beliefs,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

_EPS = 1e-12
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def _scan_and_truth(suite: Path, cond: str, index, cfg, work_dir: Path, cache_dir: Path) -> dict:
    """The expensive, solver-INDEPENDENT half: one production scan + the oracle's region AND boundary truth."""
    cache = cache_dir / f"{cond}.pkl"
    if cache.exists():
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    bam = str(suite / cond / "sim_oracle.bam")
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, sm, flm, _buf, payload = scan_and_buffer(bam, index, sc)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    orc = OracleTruth.from_bam(
        bam, index, cfg, work_dir, cond, full_payload=payload, boundary_mass_tol=0.5
    )
    out = dict(
        payload=payload,
        strand_model=sm,
        gdna_fl_pmf=np.asarray(fl.gdna_pmf),
        rna_fl_pmf=np.asarray(fl.rna_pmf),
        region_pools={k: np.asarray(v) for k, v in orc.region_pools().items()},
        boundary_pools={k: np.asarray(v) for k, v in orc.boundary_pools().items()},
    )
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return out


def _true_fg(pools):
    """TRUE gDNA fraction of the UNSPLICED mass (the deconvolution's competition basis) + that mass."""
    g = pools["gdna_pos"] + pools["gdna_neg"]
    r = pools["mat_uns_pos"] + pools["mat_uns_neg"] + pools["nas_uns_pos"] + pools["nas_uns_neg"]
    tot = g + r
    return np.where(tot > _EPS, g / np.maximum(tot, _EPS), np.nan), tot


def self_solve(inp, index, cfg, stage="both"):
    """Run the message-free self-solve (and optionally the prior-free FB sweep), aligned to oracle truth."""
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = inp["payload"]
    sub = CalibrationSubstrate.from_payload(payload, ra)
    bsub = BoundarySubstrate.from_payload(payload)
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    st = build_node_statics(chain, sub, bsub, ra)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    cc = cfg.calibration
    kappa = float(fit_strand_balance(inp["strand_model"]).rna_sense_frac)
    # THE MESSAGE-FREE SELF-SOLVE: bare psi = strand + reference. No prior, no messages.
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
    fb = None
    if stage in ("fb", "both"):
        # STEP 5: the PRIOR-FREE forward-backward pass. gdna_prior=None => psi carries the derived
        # reference alone on both arms (prior-free is not reference-free). sigma^2_transfer is 0.
        fb = node_sweep(
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
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    tf_r, m_r = _true_fg(inp["region_pools"])
    tf_b, m_b = _true_fg(inp["boundary_pools"])
    is_r = kind == REGION
    R, B = tf_r.shape[0], tf_b.shape[0]
    ri, bi = np.clip(idx, 0, R - 1), np.clip(idx, 0, B - 1)
    true_fg = np.where(is_r, tf_r[ri], tf_b[bi])
    true_mass = np.where(is_r, m_r[ri], m_b[bi])
    # node class: regions by coarse type; boundaries as one class (their flanks vary)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    cls = np.where(is_r, np.array([_RTYPE.get(int(t), "?") for t in rtype])[ri], "boundary")
    n = np.asarray(st.u_pos, float) + np.asarray(st.u_neg, float)
    var_g = np.asarray(belief.var_gdna, float)
    return pd.DataFrame(
        dict(
            kind=np.where(is_r, "region", "boundary"),
            cls=cls,
            fg=np.asarray(belief.f_g, float),
            fg_fb=(np.asarray(fb.f_g, float) if fb is not None else np.nan),
            var_g_fb=(np.asarray(fb.var_gdna, float) if fb is not None else np.nan),
            true_fg=true_fg,
            mass=true_mass,
            n=n,
            var_g=var_g,
            # the precision this node's gDNA message WOULD carry (bp_solver._scan). var_g=inf marks a
            # node with NO information (an AMBIG init holds {0,0,1} at MAX variance) -> it sends nothing;
            # n=0 likewise ("zero density is not a measurement"). Both -> pr = 0 (n*inf would be NaN).
            pr=np.where(
                (n > 0.0) & np.isfinite(var_g),
                n / (n * np.where(np.isfinite(var_g), var_g, 0.0) + 1.0),
                0.0,
            ),
            solvable=(np.asarray(st.free_pos, bool) | np.asarray(st.free_neg, bool))
            & (np.asarray(st.mass_unspliced, float) > 0.0),
        )
    ), kappa


def report(df, cond, kappa, var_ref):
    ok = np.isfinite(df.true_fg) & (df.mass > _EPS)
    d = df[ok]
    # kappa near 0 vs near 1 is the same STRAND INFORMATION (an antisense library convention); what the
    # solve can use is |1/2 - kappa|, and I(f_g) ~ n(1/2 - kappa)^2 vanishes at kappa = 1/2 exactly.
    print(
        f"\n{'=' * 100}\n{cond}\n   kappa = {kappa:.3f}   |1/2 - kappa| = {abs(0.5 - kappa):.3f}"
        f"   (0 = NO strand information, whatever the count)\n{'=' * 100}"
    )
    print(
        f"  reference-only ceiling: var_g = {var_ref:.3f}  ->  a node that learned NOTHING sends at most "
        f"pr = {1.0 / var_ref:.3f}"
    )
    rows = []
    for (k, c), g in d.groupby(["kind", "cls"]):
        w = g.mass.to_numpy()
        wsum = max(w.sum(), _EPS)
        rows.append(
            dict(
                kind=k,
                cls=c,
                n_nodes=len(g),
                mass_frac=wsum / max(d.mass.sum(), _EPS),
                mwae_fg=float(np.sum(w * np.abs(g.fg - g.true_fg)) / wsum),
                mwae_fb=float(np.sum(w * np.abs(g.fg_fb - g.true_fg)) / wsum),
                solved_fb=float(np.sum(w * g.fg_fb) / wsum),
                true_fg=float(np.sum(w * g.true_fg) / wsum),
                solved_fg=float(np.sum(w * g.fg) / wsum),
                med_var_g=float(np.median(g.var_g)),
                med_pr=float(np.median(g.pr)),
                # "knows something" = its message beats the reference-only ceiling by 2x
                frac_informed=float(np.mean(g.pr > 2.0 / var_ref)),
            )
        )
    r = pd.DataFrame(rows).sort_values(["kind", "cls"])
    print(
        r.to_string(
            index=False,
            float_format=lambda x: f"{x:,.4g}",
            columns=[
                "kind",
                "cls",
                "n_nodes",
                "mass_frac",
                "true_fg",
                "solved_fg",
                "mwae_fg",
                "solved_fb",
                "mwae_fb",
                "med_var_g",
                "med_pr",
            ],
        )
    )
    return r.assign(condition=cond, kappa=kappa)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None, help="comma-separated subset")
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--stage",
        default="both",
        choices=["self", "fb", "both"],
        help="self = step 3 (message-free); fb = + the prior-free forward-backward pass (step 5)",
    )
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = (
        a.conditions.split(",")
        if a.conditions
        else sorted(
            d.name
            for d in suite.iterdir()
            if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_")
        )
    )
    # the reference-only Var(log f_g): Beta(1/2,1/2) truncated to the solver's lambda window
    from rigel.calibration.simplex_logodds import _log_fg, _logodds_grid

    lam, _ = _logodds_grid(4001, cfg.calibration.sweep_logodds_window)
    lp = 0.5 * _log_fg(lam) + 0.5 * _log_fg(-lam)
    p = np.exp(lp - lp.max())
    p /= p.sum()
    L = _log_fg(lam)
    var_ref = float(p @ (L * L) - (p @ L) ** 2)

    out = []
    for i, c in enumerate(conds, 1):
        print(f"[{i}/{len(conds)}] {c} ...", flush=True)
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        df, kappa = self_solve(inp, index, cfg, stage=a.stage)
        out.append(report(df, c, kappa, var_ref))
    res = pd.concat(out, ignore_index=True)
    dest = Path(a.out) if a.out else suite / "selfsolve_diag.tsv"
    res.to_csv(dest, sep="\t", index=False)
    print(f"\nfull table -> {dest}")


if __name__ == "__main__":
    main()
