"""Pass-by-pass calibration error trace against the sim oracle.

Faithfully replicates calibrate()'s internal two-pass flow (init -> pass1 sweep -> KDE train ->
pass2 sweep) but captures the per-node gDNA fraction f_g at EVERY stage via node_sweep's _capture
hook, joins to the oracle per-region true gDNA, ranks regions by FINAL |gDNA mass error|, and
decomposes each top culprit's error across stages:

  init  -> p1_strand (strand-likelihood only) -> p1_loc (+global floor prior) -> p1_final (+messages)
        -> p2_loc (+KDE mixture prior)         -> p2_final (+messages)   [== production output]

    OMP_NUM_THREADS=1 python pass_trace.py <condition_dir> [topN] [deep_region_ids_csv]
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np, pandas as pd, pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import coarse_type_from_signature
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration.bp_solver import (
    build_node_geometry, build_node_statics, init_beliefs, node_sweep, node_densities,
)
from rigel.calibration.gdna_density_prior import build_training_substrate, GdnaDensityPrior
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

cond = Path(sys.argv[1]); suite = cond.parent
topN = int(sys.argv[2]) if len(sys.argv) > 2 else 15
deep_ids = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else None

idx = TranscriptIndex.load(suite / "rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
roff = np.asarray(ra.ref_offsets, np.int64); refid = np.asarray(ra.ref_id)
n = ra.n_regions

# ---- oracle TRUE per-region contained gDNA / RNA mass ----
true_g = np.zeros(n); true_r = np.zeros(n)
with pysam.AlignmentFile(str(cond / "sim_oracle.bam"), "rb") as f:
    dref = f.references[0]; seen: set[str] = set()
    for r in f:
        q = r.query_name
        if q in seen: continue
        seen.add(q); o = parse_origin(q)
        if o.start is None: continue
        rid = idx.ref_name_to_id.get(str(o.ref if o.ref is not None else dref))
        if rid is None: continue
        lo0, hi0 = int(roff[rid]), int(roff[rid + 1]); a, b = int(o.start), int(o.end)
        rr = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
        if rr < hi0 and starts[rr] <= a and b <= ends[rr]:
            (true_g if o.kind == "gdna" else true_r)[rr] += 1.0

# ---- scan + FL models ----
# CalibrationConfig; RIGEL_MIX_BRIDGE env overrides the gDNA-prior mixture-bridge ε for A/B (else the
# production config default). The bridge is now a config param (not a solver env var), so it is threaded
# into GdnaDensityPrior.fit below.
_eps = os.environ.get("RIGEL_MIX_BRIDGE")
_ng = os.environ.get("RIGEL_N_GRID")  # override sweep_n_grid for the grid-resolution study
_ckw = {}
if _eps is not None:
    _ckw["gdna_prior_mixture_bridge"] = float(_eps)
if _ng is not None:
    _ckw["sweep_n_grid"] = int(_ng)
cfg = CalibrationConfig(**_ckw)
st, sm, flm, buf, pl = scan_and_buffer(str(cond / "sim_oracle.bam"), idx, BamScanConfig())
sub = CalibrationSubstrate.from_payload(pl, ra)
bsub = BoundarySubstrate.from_payload(pl)
raw = np.asarray(sub.contained.mass_unspliced, float)          # contained unspliced (the deconv split base)
npos = np.asarray(sub.contained.n_unspliced_pos, float)
nneg = np.asarray(sub.contained.n_unspliced_neg, float)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
gdna_fl_pmf, rna_fl_pmf = fl.gdna_pmf, fl.rna_pmf

# ---- replicate calibrate() internals ----
region_eff_len = region_eff_length(ra.region_size_bp, gdna_fl_pmf)
boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, ra.region_size_bp)
fl_mean = boundary_eff_length(gdna_fl_pmf)
balance = fit_strand_balance(sm); kappa = float(balance.rna_sense_frac)
node_density_raw = node_gdna_density(sub, ra, region_eff_len, fl_mean)
gdna_strand = fit_gdna_strand_from_substrate(
    sub, ra, node_density_raw, boundary_eff_len, rna_sense_frac=kappa,
    prior_overdispersion=overdispersion_for_beta(cfg.gdna_strand_prior_alpha_beta),
    prior_weight=cfg.gdna_strand_prior_weight)
od_g = gdna_strand.gdna_strand_overdispersion
rna_strand = fit_rna_strand_from_substrate(
    sub, rna_sense_frac=kappa,
    prior_overdispersion=overdispersion_for_beta(cfg.rna_strand_prior_alpha_beta),
    prior_weight=cfg.rna_strand_prior_weight)
od_r = rna_strand.rna_strand_overdispersion

chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geometry = build_node_geometry(chain, sub, bsub, ra, gdna_fl_pmf, rna_fl_pmf)
statics = build_node_statics(chain, sub, bsub, ra)
belief0 = init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa,
                       gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
                       n_grid=cfg.sweep_n_grid, n_grid_ss=cfg.sweep_n_grid_single_strand,
                       logodds_window=cfg.sweep_logodds_window, statics=statics)

from rigel.calibration.bp_solver import adjacent_disagreement_variance
# Single production path: the belief-free total-density σ²_imp message precision (σ²_msg = σ²_imp + 1/n_src).
_sig_total = adjacent_disagreement_variance(chain, geometry)
print(f"total-density sigma2_imp={_sig_total:.4f}")

def _sweep(prior, belief, cap, dis2):
    return node_sweep(chain, statics, geometry, belief, ra, bsub, rna_sense_frac=kappa,
                      gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
                      n_grid=cfg.sweep_n_grid, max_passes=cfg.sweep_max_passes,
                      convergence_delta=cfg.sweep_convergence_delta,
                      logodds_window=cfg.sweep_logodds_window, n_tilt=cfg.sweep_n_tilt,
                      n_grid_ss=cfg.sweep_n_grid_single_strand,
                      gdna_prior=prior, disagreement_sigma2=dis2, _capture=cap)

cap1 = {}
belief1 = _sweep(None, belief0, cap1, _sig_total)
_dis2 = _sig_total  # single total-density basis for both passes
# KDE training on pass-1 belief
train_sub = build_training_substrate(chain, belief1, geometry, statics, ra, bsub, min_eff_length=fl_mean)
gdna_prior = GdnaDensityPrior.fit(train_sub, bandwidth=cfg.gdna_prior_bandwidth,
                                  mixture_bridge=cfg.gdna_prior_mixture_bridge)
cap2 = {}
belief2 = _sweep(gdna_prior, belief1, cap2, _dis2)

# ---- chain-node -> region mapping ----
kind = np.asarray(chain.kind); cref = np.asarray(chain.ref_idx, np.int64)
reg_mask = kind == REGION
reg_node = np.full(n, -1, np.int64)
reg_node[cref[reg_mask]] = np.where(reg_mask)[0]

def per_region(node_arr):
    out = np.full(n, np.nan)
    ok = reg_node >= 0
    out[ok] = np.asarray(node_arr)[reg_node[ok]]
    return out

# per-region f_g at each stage
init_fg   = per_region(belief0.f_g)
p1_strand = per_region(cap1["fg_strand"])
p1_loc    = per_region(cap1["fg_loc"])
p1_final  = per_region(cap1["f_g"])
p2_strand = per_region(cap2["fg_strand"])
p2_loc    = per_region(cap2["fg_loc"])
p2_final  = per_region(cap2["f_g"])   # == production mass_gdna_contained / raw

true_g = np.minimum(true_g, raw)
true_gfrac = true_g / np.maximum(raw, 1e-9)
obs_final = p2_final * raw
err = obs_final - true_g

sig = idx.region_df["signature"].to_numpy()
cls = np.array(["ig", "intron", "exon"])[[coarse_type_from_signature(int(s)) for s in sig]]
S = np.maximum(region_eff_len, 1e-9)
sense = npos / np.maximum(npos + nneg, 1e-9)

df = pd.DataFrame(dict(
    region=np.arange(n), cls=cls, ref=refid, start=starts, end=ends, S=S, raw=raw,
    true_g=true_g, true_gf=true_gfrac, sense=sense,
    init=init_fg, p1_str=p1_strand, p1_loc=p1_loc, p1_fin=p1_final,
    p2_loc=p2_loc, p2_fin=p2_final, obs_g=obs_final, err=err,
))
# error attributable to each stage transition, in gfrac units (x raw = mass)
df["d_strand"] = (df.p1_str - df.true_gf)         # strand solve vs truth
df["d_glob1"]  = (df.p1_loc - df.p1_str)          # pass1 global floor prior effect
df["d_msg1"]   = (df.p1_fin - df.p1_loc)          # pass1 message effect
df["d_kde"]    = (df.p2_loc - df.p1_fin)          # KDE effect on local (pass2 loc vs pass1 final)
df["d_msg2"]   = (df.p2_fin - df.p2_loc)          # pass2 message effect
df["gf_err"]   = (df.p2_fin - df.true_gf)

pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
print(f"=== {cond.name} ===")
print(f"totals: obs_g={obs_final.sum():,.0f} true_g={true_g.sum():,.0f}  Σ|err|={np.abs(err).sum():,.0f}")
print(f"net err by class:  " + "  ".join(f"{c}:{df[df.cls==c].err.sum():+,.0f}" for c in ["exon","intron","ig"]))
print(f"strand overdisp od_g={od_g:.4g} od_r={od_r:.4g}  kappa={kappa:.4f}  rho_global(p2)={cap2['rho_global']:.4g} rho_floor(p2)={cap2['rho_floor']:.4g}")
print(f"KDE: bandwidth={gdna_prior.bandwidth:.4f}  n_train={train_sub.n}  n_eff={gdna_prior.n_eff:.0f}  modes(logrho,logP)={[(round(m,3),round(p,2)) for m,p in gdna_prior.modes[:5]]}")

top = df.reindex(df.err.abs().sort_values(ascending=False).index).head(topN)
cols = ["region","cls","raw","true_g","true_gf","sense","init","p1_str","p1_loc","p1_fin","p2_loc","p2_fin","err"]
print(f"\nTOP {topN} regions by |final gDNA mass error| (f_g fractions per stage):")
print(top[cols].to_string(index=False, float_format=lambda x: f"{x:,.2f}"))

# stage decomposition, mass-weighted over exon nodes with real gDNA
ex = df[(df.cls == "exon") & (df.raw > 20)]
print("\nStage error decomposition (mass = Δgfrac x raw, summed over exon nodes raw>20):")
for lab, c in [("strand vs truth","d_strand"),("+global(p1)","d_glob1"),("+messages(p1)","d_msg1"),
               ("+KDE(p2 loc)","d_kde"),("+messages(p2)","d_msg2"),("FINAL vs truth","gf_err")]:
    net = float(np.sum(ex[c] * ex.raw)); absm = float(np.sum(np.abs(ex[c]) * ex.raw))
    print(f"  {lab:16s}  net={net:+11,.0f}  |mass|={absm:11,.0f}")

# ---- deep dive on requested regions (default: top culprit) ----
deep = deep_ids if deep_ids is not None else [int(top.iloc[0].region)]
left = np.asarray(chain.left); right = np.asarray(chain.right)
dens1 = node_densities(belief1, geometry)
dens2 = node_densities(belief2, geometry)
eff_global = cap2["eff_global"]; mass_global = cap2["mass_global"]
free_pos = cap2["free_pos"]; free_neg = cap2["free_neg"]
def msg_fg(cap, nd):
    """recover the message-implied f_g/f_pos/f_neg (mode is log-fraction in dst frame)."""
    return np.exp(cap["mode_g"][nd]), np.exp(cap["mode_p"][nd]), np.exp(cap["mode_n"][nd])
for w in deep:
    nd = reg_node[w]
    print(f"\n=== DEEP: region {w} ({cls[w]}) [{starts[w]}-{ends[w]}) chain_node={nd} ===")
    print(f"  raw={raw[w]:,.0f} true_g={true_g[w]:,.0f} true_gf={true_gfrac[w]:.3f}  "
          f"sense={sense[w]:.3f} (npos={npos[w]:,.0f} nneg={nneg[w]:,.0f})  free_pos={bool(free_pos[nd])} free_neg={bool(free_neg[nd])}")
    print(f"  f_g by stage:  init={init_fg[w]:.3f} -> p1_str={p1_strand[w]:.3f} -> p1_loc={p1_loc[w]:.3f} "
          f"-> p1_fin={p1_final[w]:.3f} -> p2_loc={p2_loc[w]:.3f} -> p2_fin={p2_final[w]:.3f}  (true {true_gfrac[w]:.3f})")
    md, eg = mass_global[nd], eff_global[nd]
    print(f"  global geom: mass_global={md:,.0f} eff_global={eg:.1f}  M/E(total dens)={md/max(eg,1e-9):.2f}  "
          f"implied f_g@rho_global={min(cap2['rho_global']*eg/max(md,1e-9),1):.3f} @rho_floor={min(cap2['rho_floor']*eg/max(md,1e-9),1):.3f}")
    d_tot = md / max(eg, 1e-9)
    for tag, fgv in [("p1_fin", p1_final[w]), ("true", true_gfrac[w])]:
        lr = np.log(max(fgv * d_tot, 1e-9))
        lp = float(gdna_prior.logpdf_kernel(np.array([lr]))[0])
        print(f"    KDE @ f_g={fgv:.3f} -> log_rho_g={lr:.3f}  logP_kernel={lp:.2f}")
    # component beliefs + all 3 messages, BOTH passes
    for pn, cap, bel in [("PASS1", cap1, belief1), ("PASS2", cap2, belief2)]:
        fgm, fpm, fnm = msg_fg(cap, nd)
        print(f"  [{pn}] local pie: f_g={cap['f_g'][nd] if 'f_g' in cap else float('nan'):.3f}  "
              f"fg_loc={cap['fg_loc'][nd]:.3f} fp_loc={cap['fp_loc'][nd]:.3f} fn_loc={cap['fn_loc'][nd]:.3f}")
        print(f"         msgs -> gDNA: f_g={fgm:.3f} prec={cap['prec_g'][nd]:.3f} | "
              f"RNA+: f={fpm:.3f} prec={cap['prec_p'][nd]:.3f} | RNA-: f={fnm:.3f} prec={cap['prec_n'][nd]:.3f}  "
              f"(local prec_g~{1.0/max(cap['vg_loc'][nd],1e-9):.1f})")
    print("  NEIGHBOURS (chain adjacency), f_g/f_pos/f_neg + facing gDNA density (pass1 -> pass2):")
    for tag, nb in [("L", left[nd]), ("R", right[nd])]:
        if nb < 0: continue
        rr = cref[nb] if kind[nb] == REGION else -1
        knd = "REG" if kind[nb] == REGION else "BND"
        rg1 = float(dens1.rho_g_right[nb]) if tag == "L" else float(dens1.rho_g_left[nb])  # face toward w
        rg2 = float(dens2.rho_g_right[nb]) if tag == "L" else float(dens2.rho_g_left[nb])
        rp2 = float(dens2.rho_pos_right[nb]) if tag == "L" else float(dens2.rho_pos_left[nb])
        rn2 = float(dens2.rho_neg_right[nb]) if tag == "L" else float(dens2.rho_neg_left[nb])
        info = f" reg={rr} true_gf={true_gfrac[rr]:.2f} sense={sense[rr]:.2f} raw={raw[rr]:,.0f}" if rr >= 0 else ""
        print(f"    {tag}: node {nb} [{knd}] p1(fg={belief1.f_g[nb]:.3f}) p2(fg={belief2.f_g[nb]:.3f} fp={belief2.f_pos[nb]:.3f} fn={belief2.f_neg[nb]:.3f})  "
              f"rho_g_face {rg1:.2f}->{rg2:.2f}  rho_p/n_face(p2) {rp2:.2f}/{rn2:.2f}{info}")

# dump full table
out = cond / "pass_trace.tsv"; df.to_csv(out, sep="\t", index=False)
print(f"\nfull per-region trace -> {out}")

# ---- self-contained NPZ for downstream verification (no BAM re-scan needed) ----
npz = Path(os.environ.get("PASS_TRACE_NPZ", str(cond / "pass_trace.npz")))
np.savez_compressed(
    npz,
    # per-region (length n)
    reg_true_g=true_g, reg_true_r=true_r, reg_raw=raw, reg_cls=cls.astype("U8"),
    reg_sense=sense, reg_S=S, reg_start=starts, reg_end=ends, reg_ref=refid,
    reg_init=init_fg, reg_p1_str=p1_strand, reg_p1_loc=p1_loc, reg_p1_fin=p1_final,
    reg_p2_loc=p2_loc, reg_p2_fin=p2_final, reg_node=reg_node,
    reg_npos=npos, reg_nneg=nneg,
    # per chain-node (length n_nodes)
    chain_kind=kind, chain_ref_idx=cref, chain_left=left, chain_right=right,
    b1_fg=belief1.f_g, b1_fp=belief1.f_pos, b1_fn=belief1.f_neg,
    b2_fg=belief2.f_g, b2_fp=belief2.f_pos, b2_fn=belief2.f_neg,
    b1_vg=belief1.var_gdna, b2_vg=belief2.var_gdna,
    eff_global=eff_global, mass_global=mass_global,
    free_pos=free_pos.astype(bool), free_neg=free_neg.astype(bool),
    # pass-1 capture
    p1_fg_strand=cap1["fg_strand"], p1_fg_loc=cap1["fg_loc"], p1_fp_loc=cap1["fp_loc"], p1_fn_loc=cap1["fn_loc"],
    p1_mode_g=cap1["mode_g"], p1_prec_g=cap1["prec_g"], p1_mode_p=cap1["mode_p"], p1_prec_p=cap1["prec_p"],
    p1_mode_n=cap1["mode_n"], p1_prec_n=cap1["prec_n"], p1_vg_loc=cap1["vg_loc"], p1_vp_loc=cap1["vp_loc"], p1_vn_loc=cap1["vn_loc"],
    # pass-2 capture
    p2_fg_strand=cap2["fg_strand"], p2_fg_loc=cap2["fg_loc"], p2_fp_loc=cap2["fp_loc"], p2_fn_loc=cap2["fn_loc"],
    p2_mode_g=cap2["mode_g"], p2_prec_g=cap2["prec_g"], p2_mode_p=cap2["mode_p"], p2_prec_p=cap2["prec_p"],
    p2_mode_n=cap2["mode_n"], p2_prec_n=cap2["prec_n"], p2_vg_loc=cap2["vg_loc"], p2_vp_loc=cap2["vp_loc"], p2_vn_loc=cap2["vn_loc"],
    # node face densities (belief2 + belief1)
    d2_rho_g_l=dens2.rho_g_left, d2_rho_g_r=dens2.rho_g_right,
    d2_rho_p_l=dens2.rho_pos_left, d2_rho_p_r=dens2.rho_pos_right,
    d2_rho_n_l=dens2.rho_neg_left, d2_rho_n_r=dens2.rho_neg_right,
    d1_rho_g_l=dens1.rho_g_left, d1_rho_g_r=dens1.rho_g_right,
    # KDE
    kde_x=gdna_prior.x_grid, kde_logP=gdna_prior.logP_grid, kde_bw=gdna_prior.bandwidth,
    kde_train_x=train_sub.log_rho, kde_train_kind=train_sub.node_kind,
    kde_modes=np.array(gdna_prior.modes, float) if gdna_prior.modes else np.zeros((0, 2)),
    # scalars
    kappa=kappa, od_g=od_g, od_r=od_r, rho_global=cap2["rho_global"], rho_floor=cap2["rho_floor"],
    fl_mean=fl_mean,
)
print(f"self-contained NPZ -> {npz}  ({npz.stat().st_size/1e6:.1f} MB)")
# dump KDE training distribution summary
tk = train_sub.node_kind; tx = train_sub.log_rho
print("\nKDE training-node log_rho by kind:")
for kc, nm in [(0,"intergenic"),(1,"intron"),(2,"exon"),(3,"boundary")]:
    m = tk == kc
    if m.sum():
        xr = tx[m]
        print(f"  {nm:11s} n={m.sum():5d}  log_rho: p10={np.percentile(xr,10):.2f} med={np.median(xr):.2f} p90={np.percentile(xr,90):.2f}  (rho med={np.exp(np.median(xr)):.2f})")
