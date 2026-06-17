"""Illustrate the increment-4 simplex propagation on a flagship condition (4 panels → PNG).

Shows how the count ``var~mean`` fit and the gDNA-density propagation look on real data:

  A. the var~mean FIT — per-node triplet disagreement (μ, raw_var) in log-log, coloured by enrichment
     class (on-target exon vs off-target), with the LOESS σ²_density(μ) curve overlaid. **Use this to
     check the point-vs-interval scale caveat:** if the count-observable *contained* points (triplet) and
     the boundary-only points form two separated clouds, the eff-length normalisations disagree.
  B. the derived per-hop process variance Q vs local density (what the propagation decay uses).
  C. propagation accuracy — smoothed ρ_g vs ORACLE gDNA density per region (identity line).
  D. f_g calibration — predicted region f_g vs oracle gDNA fraction (identity line), by strand class.

Usage:  python scripts/debug/plot_simplex_var_mean.py [condition] [out.png]
"""
import dataclasses
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

from rigel.calibration.density_model import (
    count_observable_masks, density_variance_curve, node_gdna_density,
)
from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.run_fill import same_ref_left_right
from rigel.calibration.signature import RegionType, TS_AMBIG, TS_NEG, TS_NONE, TS_POS, coarse_type_array
from rigel.calibration.simplex_propagate import _coupling_process_var, propagate_simplex
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
SC = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}
SC_COLOR = {TS_POS: "tab:blue", TS_NEG: "tab:orange", TS_NONE: "0.6", TS_AMBIG: "tab:red"}


def oracle_contained(bam, starts, ids, R):
    g, m, n = np.zeros(R), np.zeros(R), np.zeros(R)
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i])
            k = parse_origin(r.query_name).kind
            (g if k == "gdna" else n if k == "nrna" else m)[rid] += 1.0
    return g, m, n


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    out = sys.argv[2] if len(sys.argv) > 2 else f"{SUITE}/{cond}/simplex_var_mean.png"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    ccfg = cfg.calibration
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    reg_el = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bnd_el = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion

    regions, _left, _right = propagate_simplex(
        sub, ra, gdna_region_eff_len=reg_el, gdna_boundary_side_eff_len=bnd_el, gdna_fl_mean=fl_mean,
        rna_sense_frac=kappa, gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
        n_grid=ccfg.n_grid,
    )

    # Recompute the triplet disagreement points + the fitted curve (for panel A/B).
    sig = np.asarray(ra.signature)
    ref_id = np.asarray(ra.ref_id)
    R = sig.shape[0]
    L = np.maximum(reg_el, 1e-9)
    inv_fl = 1.0 / fl_mean

    def total(view):
        return view.n_unspliced_pos.astype(float) + view.n_unspliced_neg.astype(float)

    region_obs, boundary_obs = count_observable_masks(sig, ref_id)
    lsame, rsame = same_ref_left_right(ref_id)
    la = np.zeros(R, bool)
    rb = np.zeros(R, bool)
    if R > 1:
        la[1:] = boundary_obs[:-1] & lsame[1:]
        rb[:-1] = boundary_obs[:-1] & rsame[:-1]
    d_left = np.where(la, total(sub.left) * inv_fl, np.nan)
    d_right = np.where(rb, total(sub.right) * inv_fl, np.nan)
    c_ok = region_obs & (L > 1e-9)
    contained = np.where(c_ok, total(sub.contained) / L, np.nan)
    # per-node mean/raw_var over available clean observations (same statistic as density_variance_curve)
    obs = np.stack([d_left, d_right, contained])
    msk = np.stack([la, rb, c_ok])
    k = msk.sum(0).astype(float)
    mu = np.where(msk, np.where(np.isnan(obs), 0.0, obs), 0.0).sum(0) / np.maximum(k, 1.0)
    dev2 = np.where(msk, (np.nan_to_num(obs) - mu) ** 2, 0.0).sum(0)
    raw_var = np.where(k > 1.0, dev2 / np.maximum(k - 1.0, 1e-9), np.nan) / np.maximum(k, 1.0)
    is_triplet = k >= 3.0
    fit_pt = (k >= 2.0) & np.isfinite(mu) & (mu > 1e-9) & np.isfinite(raw_var) & (raw_var > 1e-9)

    q = _coupling_process_var(sub, ra, reg_el, fl_mean)
    on_target = coarse_type_array(sig) == int(RegionType.EXON)
    sigma_curve = density_variance_curve(  # σ²_density(μ) evaluated at each node's density (the curve)
        np.asarray(nd.density), d_left=d_left, d_right=d_right, left_ok=la, right_ok=rb,
        contained=contained, contained_ok=c_ok,
    )

    # oracle
    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    og, om, on = oracle_contained(bam, starts, ids, R)
    o_tot = og + om + on
    o_gfrac = np.where(o_tot > 0, og / np.maximum(o_tot, 1e-9), np.nan)
    o_gdens = np.where(reg_el > 0, og / np.maximum(reg_el, 1e-9), np.nan)
    rho_post = np.asarray(regions.gdna_mass) / np.maximum(reg_el, 1e-9)  # gdna density (mass = count)
    fg = np.asarray(regions.gdna_frac)
    ts = np.asarray(ra.strand_class)

    fig, ax = plt.subplots(2, 2, figsize=(13, 11))
    fig.suptitle(f"simplex propagation — {cond}\nκ={kappa:.3f}  od_g={od_g:.3f}  od_r={od_r:.3f}",
                 fontsize=11)

    # A — var~mean fit (log-log), coloured by triplet vs pair and on/off-target
    a = ax[0, 0]
    pair = fit_pt & ~is_triplet
    trip = fit_pt & is_triplet
    a.scatter(mu[pair], raw_var[pair], s=12, c="0.5", alpha=0.5, label="boundary pair (k=2)")
    a.scatter(mu[trip & on_target], raw_var[trip & on_target], s=14, c="tab:green", alpha=0.6,
              label="triplet, on-target")
    a.scatter(mu[trip & ~on_target], raw_var[trip & ~on_target], s=14, c="tab:purple", alpha=0.6,
              label="triplet, off-target")
    order = np.argsort(nd.density)
    good = np.isfinite(sigma_curve[order]) & (nd.density[order] > 1e-9)
    a.plot(nd.density[order][good], sigma_curve[order][good], "k-", lw=2, label="LOESS σ²_density(μ)")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("density μ (gDNA frags / eff-bp)")
    a.set_ylabel("raw_var = s²/k  (variance of the mean)")
    a.set_title("A. var~mean fit (triplet disagreement + LOESS)")
    a.legend(fontsize=7, loc="best")

    # B — derived Q vs density
    b = ax[0, 1]
    finite_q = np.isfinite(q) & (q > 0)
    b.scatter(nd.density[finite_q & on_target], q[finite_q & on_target], s=10, c="tab:green",
              alpha=0.5, label="on-target")
    b.scatter(nd.density[finite_q & ~on_target], q[finite_q & ~on_target], s=10, c="tab:purple",
              alpha=0.5, label="off-target")
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel("local density μ")
    b.set_ylabel("per-hop process variance Q = σ²_density(μ)")
    b.set_title(f"B. derived coupling Q  ({int(finite_q.sum())}/{R} from curve)")
    b.legend(fontsize=7)

    # C — propagated ρ_g vs oracle density
    c = ax[1, 0]
    for cv in (TS_NONE, TS_POS, TS_NEG, TS_AMBIG):
        m = (ts == cv) & np.isfinite(o_gdens) & (o_tot > 0)
        c.scatter(o_gdens[m], rho_post[m], s=10, alpha=0.5, c=SC_COLOR[cv], label=SC[cv])
    lim = np.nanpercentile(np.concatenate([o_gdens[o_tot > 0], rho_post[o_tot > 0]]), 99)
    c.plot([0, lim], [0, lim], "k--", lw=1)
    c.set_xlim(0, lim)
    c.set_ylim(0, lim)
    c.set_xlabel("oracle gDNA density (og / eff-len)")
    c.set_ylabel("propagated ρ_g")
    c.set_title("C. propagation accuracy (ρ_g vs oracle)")
    c.legend(fontsize=7)

    # D — f_g vs oracle gDNA fraction
    d = ax[1, 1]
    for cv in (TS_NONE, TS_POS, TS_NEG, TS_AMBIG):
        m = (ts == cv) & np.isfinite(o_gfrac)
        d.scatter(o_gfrac[m], fg[m], s=10, alpha=0.5, c=SC_COLOR[cv], label=SC[cv])
    d.plot([0, 1], [0, 1], "k--", lw=1)
    d.set_xlim(-0.02, 1.02)
    d.set_ylim(-0.02, 1.02)
    d.set_xlabel("oracle gDNA fraction")
    d.set_ylabel("predicted f_g")
    d.set_title("D. f_g calibration")
    d.legend(fontsize=7)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")
    # quick numeric summary
    amb = (ts == TS_AMBIG) & (o_tot > 0)
    print(f"AMBIG nodes={int(amb.sum())}  mean f_g={np.nanmean(fg[amb]):.3f}  "
          f"mean oracle gfrac={np.nanmean(o_gfrac[amb]):.3f}")
    print(f"Q from curve: {int((np.isfinite(q) & (q > 0)).sum())}/{R} nodes; "
          f"triplet fit points={int(trip.sum())}, pair fit points={int(pair.sum())}")


if __name__ == "__main__":
    main()
