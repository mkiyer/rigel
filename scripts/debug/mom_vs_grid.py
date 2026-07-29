"""A/B the INIT self-solve: the proposed closed-form MoM construction vs the current full GRID solver,
both vs the oracle, at the message-free initialization stage. Value (f_g) is the primary comparison;
precision is reported for structure (finite=solved / inf=unsolved). Runs on real ambig-suite BAMs.

Node classes: G1 (locked sink), single-strand (G2), AMBIG (G3); region vs boundary.
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
from rigel.calibration.node_geometry import build_node_statics, init_beliefs
from rigel.calibration.bp_solver import chain_region_deconv
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_sj_table)
from rigel.calibration.density_model import node_gdna_density
from rigel.splice import SpliceType
from oracle import OracleTruth

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
EPS = 1e-9


def mom_init(st, kappa):
    """The proposed closed-form init: nascent~0 default value + strand peel; composition (φ=0) + strand-Fisher
    precision. Returns (f_g, f_pos, f_neg, var_gdna)."""
    up = np.asarray(st.u_pos, float); un = np.asarray(st.u_neg, float)
    N = up + un
    Ns = np.asarray(st.mass_spliced, float)          # spliced count proxy (independent RNA pool)
    fpF = np.asarray(st.free_pos, bool); fnF = np.asarray(st.free_neg, bool)
    ss_pos = fpF & ~fnF; ss_neg = fnF & ~fpF; ambig = fpF & fnF; g1 = ~fpF & ~fnF
    # --- VALUE: nascent~0 default (f_g=1), strand peels nascent for single-strand ---
    f_g = np.ones_like(N); f_pos = np.zeros_like(N); f_neg = np.zeros_like(N)
    tk = 2.0 * kappa - 1.0
    if abs(tk) > 1e-6:
        with np.errstate(invalid="ignore", divide="ignore"):
            nasc_p = np.clip(np.where(N > 0, (up - un) / tk, 0.0), 0.0, N)   # + node: excess sense→nascent+
            nasc_n = np.clip(np.where(N > 0, (un - up) / tk, 0.0), 0.0, N)
            fr_p = np.where(N > 0, nasc_p / np.maximum(N, EPS), 0.0)
            fr_n = np.where(N > 0, nasc_n / np.maximum(N, EPS), 0.0)
        f_g = np.where(ss_pos, 1.0 - fr_p, f_g); f_pos = np.where(ss_pos, fr_p, f_pos)
        f_g = np.where(ss_neg, 1.0 - fr_n, f_g); f_neg = np.where(ss_neg, fr_n, f_neg)
    # AMBIG keeps f_g=1 (strand can't see the magnitude — nullspace); G1 keeps f_g=1 (locked below)
    # --- PRECISION (φ=0 pure-Poisson composition + strand Fisher; structure only) ---
    pi_comp = np.where((N > EPS) & (Ns > EPS), N * Ns / (N + Ns), 0.0)      # harmonic; 0 if no spliced
    p = 0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    single = ss_pos | ss_neg
    pi_strand = np.where(single, f_g**2 * N * (kappa - 0.5)**2 / np.maximum(p * (1.0 - p), EPS), 0.0)
    prec = pi_comp + pi_strand
    var_g = np.where(prec > EPS, 1.0 / prec, np.inf)
    var_g = np.where(g1, 0.0, var_g)                                        # locked gDNA sink
    return f_g, f_pos, f_neg, var_g


def region_fg(belief, chain, substrate):
    """Per-region f_g (gDNA fraction of contained unspliced) from a chain belief, via the production projection."""
    d = chain_region_deconv(chain, belief, substrate)
    return np.asarray(d.gdna_frac, float)


def run(cond):
    bam = str(SUITE / cond / "sim_oracle.bam")
    idx = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_index(idx)
    cfg = CalibrationConfig()
    st_m, sm, flm, buf, pl = scan_and_buffer(
        bam, idx, dc(BamScanConfig(), sj_strand_tag=_native_detect_sj_tag(bam)))
    sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    # replicate calibrate.py fitting up to init
    reg_eff = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bnd_eff = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    flmean = boundary_eff_length(fl.gdna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    ndens = node_gdna_density(sub, ra, reg_eff, flmean)
    gs = fit_gdna_strand_from_substrate(sub, ra, ndens, bnd_eff, rna_sense_frac=kappa)
    rs = fit_rna_strand_from_sj_table(sm.sj_table, rna_sense_frac=kappa)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    statics = build_node_statics(chain, sub, bsub, ra)

    grid = init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa,
        gdna_strand_overdispersion=gs.gdna_strand_overdispersion,
        rna_strand_overdispersion=rs.rna_strand_overdispersion,
        n_grid=cfg.sweep_n_grid, n_grid_ss=cfg.sweep_n_grid_single_strand,
        logodds_window=cfg.sweep_logodds_window, statics=statics)
    from rigel.calibration.node_geometry import NodeBelief
    mfg, mfp, mfn, mvg = mom_init(statics, kappa)
    mom = NodeBelief(f_pos=mfp, f_neg=mfn, f_g=mfg,
                     var_pos=np.zeros_like(mfg), var_neg=np.zeros_like(mfg), var_gdna=mvg)

    fg_grid_r = region_fg(grid, chain, sub)
    fg_mom_r = region_fg(mom, chain, sub)

    orc = OracleTruth.from_bam(bam, idx, PipelineConfig(), Path("/tmp/mvg_split"), cond,
                               full_payload=pl, boundary_mass_tol=0.5)
    p = orc.region_pools()
    gd = p["gdna_pos"] + p["gdna_neg"]
    rna_u = p["mat_uns_pos"] + p["mat_uns_neg"] + p["nas_uns_pos"] + p["nas_uns_neg"]
    with np.errstate(invalid="ignore", divide="ignore"):
        fg_true = np.where((gd + rna_u) > 0, gd / (gd + rna_u), np.nan)

    # per-region strand class + mass, to break down the comparison
    ts = np.asarray(ra.strand_class)
    massU = np.asarray(sub.contained.mass_unspliced, float)
    has = (gd + rna_u) > 0
    def cls_mask(name):
        from rigel.calibration.signature import TS_POS, TS_NEG, TS_AMBIG
        if name == "AMBIG": return has & (ts == TS_AMBIG)
        if name == "single": return has & ((ts == TS_POS) | (ts == TS_NEG))
        return has & np.ones_like(has)
    print(f"\n=== {cond}   κ={kappa:.3f}  od_g={gs.gdna_strand_overdispersion:.3f} ===")
    print(f"    {'class':10} {'n_reg':>6} {'massU':>10}   {'MAE(grid)':>10} {'MAE(MoM)':>10}   "
          f"{'mean f_g: true / grid / MoM':>32}")
    for name in ("ALL", "single", "AMBIG"):
        m = cls_mask(name) & np.isfinite(fg_true)
        if not m.any():
            continue
        w = massU[m]; wsum = max(w.sum(), EPS)
        mae_g = float(np.sum(w * np.abs(fg_grid_r[m] - fg_true[m])) / wsum)
        mae_m = float(np.sum(w * np.abs(fg_mom_r[m] - fg_true[m])) / wsum)
        mt = float(np.sum(w * fg_true[m]) / wsum)
        mg = float(np.sum(w * fg_grid_r[m]) / wsum)
        mm = float(np.sum(w * fg_mom_r[m]) / wsum)
        print(f"    {name:10} {int(m.sum()):6d} {w.sum():10.0f}   {mae_g:10.3f} {mae_m:10.3f}   "
              f"{mt:8.3f} /{mg:8.3f} /{mm:8.3f}")
    # boundary precision structure: how many boundaries become 'solved' (finite var) under MoM
    kind = np.asarray(chain.kind); bnd = kind == BOUNDARY
    gvar = np.asarray(grid.var_gdna, float)[bnd]; mvar = np.asarray(mom.var_gdna, float)[bnd]
    sig = (np.asarray(statics.u_pos) + np.asarray(statics.u_neg) + np.asarray(statics.mass_unspliced)
           + np.asarray(statics.mass_spliced))[bnd] > 0
    print(f"    boundaries with data: {int(sig.sum())}   finite-var (solved): "
          f"grid {int(np.isfinite(gvar[sig]).sum())}  MoM {int(np.isfinite(mvar[sig]).sum())}")


if __name__ == "__main__":
    conds = sys.argv[1:] or [
        "gdna_gdna300_ss_0.50_nrna_none_capture_on",     # flagship: unstranded + capture
        "gdna_gdna300_ss_0.99_nrna_none_capture_on",     # stranded control
    ]
    for c in conds:
        run(c)
