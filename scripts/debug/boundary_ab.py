"""A/B the BOUNDARY self-solve at initialization: current production init vs the §6 boundary self-solve,
both vs oracle boundary-crossing truth. Shows per-boundary VALUES (f_g) and PRECISIONS.

A boundary's unspliced-crossing spans the intron => nascent + gDNA ONLY (mature crossing is spliced).
So the truth is f_g = gdna_unspliced_crossing / (gdna + nascent unspliced_crossing); nrna_none => 1.
The §6 self-solve: VALUE = nascent~0 default (f_g->1), strand peels nascent (confidence-weighted, no
blowup at k=1/2); PRECISION = composition of the two independent count pools (unspliced + spliced).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
from dataclasses import replace as dc
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.effective_length import region_eff_length, boundary_side_eff_length, boundary_eff_length
from rigel.calibration.node_chain import build_node_chain, BOUNDARY
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta)
from rigel.calibration.density_model import node_gdna_density
from rigel.splice import SpliceType
from oracle import OracleTruth

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
EPS = 1e-9


def estimate_phi(chain, geometry):
    """Overdispersion φ from adjacent-pair log-gDNA-density disagreement, Poisson-subtracted, MEDIAN
    (robust → the WITHIN-regime typical disagreement, not the across-regime enrichment jumps that the
    pooled mean is dominated by). Var_within ≈ 2φ (two adjacent nodes' overdispersions) ⇒ φ = median/2.
    (An upper bound: within-regime pairs still carry some real local density variation.)"""
    from rigel.calibration.bp_solver import _adjacent_log_density_residuals
    resid, n_i, n_j = _adjacent_log_density_residuals(chain, geometry)
    dc = resid - np.median(resid)
    pois = 1.0 / n_i + 1.0 / n_j
    excess = np.maximum(dc * dc - pois, 0.0)
    return float(np.median(excess) / 2.0)


def boundary_self_solve(u_pos, u_neg, free_pos, free_neg, mass_uns, mass_spl, kappa, phi):
    """§6 boundary self-solve with the HONEST count precision.
    VALUE: nascent~0 default f_g=1; strand peels nascent ONLY on single-strand nodes (AMBIG magnitude is
      strand-INVISIBLE — p depends only on the tilt — so AMBIG holds f_g=1, needs the sweep). Peel weighted
      by strand confidence conf=(2k-1)²N/(1+(2k-1)²N) → 0 at k=½ (no blowup; gDNA symmetric ⇒ f_g=1).
    PRECISION (per-component COUNT precision, capped by φ): gDNA from the unspliced-crossing (all gDNA under
      nascent~0) prec=N_u/(1+N_u·φ); mature from the spliced prec=N_s/(1+N_s·φ) (separate RNA channel)."""
    N = u_pos + u_neg
    single = free_pos ^ free_neg
    tk = 2.0 * kappa - 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        mom = np.clip(np.abs(u_pos - u_neg) / np.maximum(abs(tk) * N, EPS), 0.0, 1.0)
    fisher = (tk * tk) * N
    conf = fisher / (1.0 + fisher)
    nascent_frac = np.where(single & (N > 0), conf * mom, 0.0)     # AMBIG → 0 (strand-invisible)
    f_g = 1.0 - nascent_frac
    prec_gdna = np.where(N > EPS, N / (1.0 + N * phi), 0.0)        # honest gDNA count precision
    var_g = np.where(prec_gdna > EPS, 1.0 / prec_gdna, np.inf)
    prec_mat = np.where(mass_spl > EPS, mass_spl / (1.0 + mass_spl * phi), 0.0)  # mature channel precision
    return f_g, var_g, prec_gdna, prec_mat


def run(cond):
    bam = str(SUITE / cond / "sim_oracle.bam")
    idx = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cfg = CalibrationConfig()
    _st, sm, flm, _buf, pl = scan_and_buffer(
        bam, idx, dc(BamScanConfig(), sj_strand_tag=_native_detect_sj_tag(bam)))
    sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    reg_eff = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bnd_eff = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    flmean = boundary_eff_length(fl.gdna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    ndens = node_gdna_density(sub, ra, reg_eff, flmean)
    gs = fit_gdna_strand_from_substrate(sub, ra, ndens, bnd_eff, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta),
        prior_weight=cfg.gdna_strand_prior_weight)
    rs = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(cfg.rna_strand_prior_alpha_beta),
        prior_weight=cfg.rna_strand_prior_weight)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    statics = build_node_statics(chain, sub, bsub, ra)
    geometry = build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
    phi = estimate_phi(chain, geometry)

    prod = init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa,
        gdna_strand_overdispersion=gs.gdna_strand_overdispersion,
        rna_strand_overdispersion=rs.rna_strand_overdispersion,
        n_grid=cfg.sweep_n_grid, n_grid_ss=cfg.sweep_n_grid_single_strand,
        logodds_window=cfg.sweep_logodds_window, statics=statics)

    # §6 boundary self-solve over ALL chain nodes (only boundary nodes examined below)
    s6_fg, s6_var, s6_pc, s6_pm = boundary_self_solve(
        np.asarray(statics.u_pos), np.asarray(statics.u_neg),
        np.asarray(statics.free_pos), np.asarray(statics.free_neg),
        np.asarray(statics.mass_unspliced), np.asarray(statics.mass_spliced), kappa, phi)

    # oracle boundary-crossing truth: f_g = gdna_uns / (gdna_uns + nascent_uns), per boundary
    orc = OracleTruth.from_bam(bam, idx, PipelineConfig(), Path("/tmp/bnd_split"), cond,
                               full_payload=pl, boundary_mass_tol=0.5)
    bg = BoundarySubstrate.from_payload(orc.parts["gdna"])
    bn = BoundarySubstrate.from_payload(orc.parts["nrna"])
    gd_cross = np.asarray(bg.left.mass_unspliced, float) + np.asarray(bg.right.mass_unspliced, float)
    na_cross = np.asarray(bn.left.mass_unspliced, float) + np.asarray(bn.right.mass_unspliced, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        fg_true_b = np.where((gd_cross + na_cross) > 0, gd_cross / (gd_cross + na_cross), np.nan)

    # chain boundary nodes -> boundary index
    kind = np.asarray(chain.kind); bnd = kind == BOUNDARY
    bidx = np.asarray(chain.ref_idx, np.int64)[bnd]
    prod_fg = np.asarray(prod.f_g, float)[bnd]; prod_var = np.asarray(prod.var_gdna, float)[bnd]
    s6fg = s6_fg[bnd]; s6var = s6_var[bnd]; s6pc = s6_pc[bnd]
    mass_spl_b = np.asarray(statics.mass_spliced, float)[bnd]
    mass_uns_b = np.asarray(statics.mass_unspliced, float)[bnd]
    ftrue = fg_true_b[bidx]
    ts = np.asarray(ra.strand_class)
    # boundary continuity class from statics
    fp = np.asarray(statics.free_pos, bool)[bnd]; fn = np.asarray(statics.free_neg, bool)[bnd]
    ambig = fp & fn; single = fp ^ fn

    print(f"\n=== {cond}   κ={kappa:.3f}   φ(overdispersion)={phi:.4f}  → precision cap 1/φ={1/max(phi,1e-9):.1f} ===")
    has = np.isfinite(ftrue) & (mass_uns_b > 0)
    def block(name, m):
        m = m & has
        if not m.any():
            print(f"    {name:20} (none)"); return
        w = mass_uns_b[m]; ws = max(w.sum(), EPS)
        mae_p = float(np.sum(w * np.abs(prod_fg[m] - ftrue[m])) / ws)
        mae_6 = float(np.sum(w * np.abs(s6fg[m] - ftrue[m])) / ws)
        print(f"    {name:20} n={int(m.sum()):4d}  massU={w.sum():9.0f}   MAE prod={mae_p:.3f} §6={mae_6:.3f}   "
              f"mean f_g true/prod/§6 = {np.sum(w*ftrue[m])/ws:.3f}/{np.sum(w*prod_fg[m])/ws:.3f}/"
              f"{np.sum(w*s6fg[m])/ws:.3f}")
    for nm, msk in (("ALL boundaries", np.ones_like(has)), ("single-strand", single), ("AMBIG", ambig)):
        block(nm, msk)
    # precision: §6 gДNA precision = honest COUNT precision N_u/(1+N_u·φ), finite for any unspliced.
    dat = (mass_uns_b + mass_spl_b) > 0
    pc = s6pc[dat & (mass_uns_b > 0)]
    print(f"    §6 gDNA precision = count N_u/(1+N_u·φ)  [median {np.median(pc):.1f}, cap 1/φ={1/max(phi,1e-9):.1f}]"
          f"   finite: §6 {int((s6pc[dat] > 0).sum())}/{int(dat.sum())}  prod {int(np.isfinite(prod_var[dat]).sum())}")
    # example boundaries (single-strand, with data) — values + honest precisions side by side
    ex = np.where(has & single & (mass_uns_b > 5))[0][:4]
    for i in ex:
        pp = 1 / prod_var[i] if np.isfinite(prod_var[i]) and prod_var[i] > 0 else 0.0
        print(f"      ex b{bidx[i]:5d}: N_uns={mass_uns_b[i]:6.0f} N_spl={mass_spl_b[i]:5.0f}  "
              f"f_g true={ftrue[i]:.3f} prod={prod_fg[i]:.3f} §6={s6fg[i]:.3f}   "
              f"prec: prod(grid)={pp:8.2f}  §6(count)={s6pc[i]:7.2f}")


if __name__ == "__main__":
    conds = sys.argv[1:] or [
        "gdna_gdna300_ss_0.50_nrna_none_capture_on",       # flagship: unstranded, no nascent -> true f_g=1
        "gdna_gdna300_ss_0.99_nrna_none_capture_on",       # stranded, no nascent -> true f_g=1
        "gdna_gdna300_ss_0.99_nrna_present_capture_on",    # stranded + nascent -> strand should peel
        "gdna_gdna300_ss_0.50_nrna_present_capture_on",    # unstranded + nascent -> §6 limit (overshoot)
    ]
    for c in conds:
        run(c)
