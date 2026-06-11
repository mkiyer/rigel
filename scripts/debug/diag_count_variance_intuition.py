"""Phase-2 research: understand the count posterior variance + the strand/count blend.

Reconstructs the calibration internals per region for several SS regimes on the complex quick suite and
answers three questions, with a figure (/tmp/count_variance_intuition.png):

  (A) What does the count variance look like? Poisson floor + imputation (variance~mean) LOESS, by node
      type (observable own-count vs imputed-from-anchors). Boundary Poisson floor shown for scale.
  (B) The blend: current w=(2κ−1)² (a per-library scalar) vs the Fisher/inverse-variance blend
      w_fisher = σ²_count/(σ²_strand+σ²_count) (per node). Does it hand strand→count off as SS drops?
  (C) THE decision: is count CONFIDENCE (small σ_count) a proxy for count ACCURACY (|μ_count−truth|)?
      If yes, the Fisher blend is trustworthy; if confident=confidently-biased, keep (2κ−1)².

Truth: oracle gDNA fraction of the unspliced mass per region (gdna* vs GENE* read names).
"""
import numpy as np, pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import replace as dcr
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length, boundary_side_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import deconv_regions
from rigel.calibration.gdna_strand import fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
EXON = BIT_EXON_POS | BIT_EXON_NEG
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
df = idx.region_df
sig = df["signature"].to_numpy()
starts, ids = {}, {}
for ref, g in df.groupby("ref_name"):
    starts[ref] = g["start"].to_numpy(); ids[ref] = g["region_id"].to_numpy()
R = len(df)


def oracle_unspliced_gdna_frac(cond):
    ug = np.zeros(R); um = np.zeros(R)
    for r in pysam.AlignmentFile(f"{SUITE}/{cond}/sim_oracle.bam", "rb").fetch(until_eof=True):
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
        (ug if r.query_name.startswith("gdna") else um)[rid] += 1.0
    U = ug + um
    return np.where(U > 0, ug / np.maximum(U, 1e-9), np.nan), U


def reconstruct(cond):
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    cfg = PipelineConfig(); scan = dcr(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    rel = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bel = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    nd = node_gdna_density(sub, ra, rel, boundary_eff_length(fl.gdna_pmf))
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    cfgc = cfg.calibration
    gd = fit_gdna_strand_from_substrate(sub, ra, nd, bel, rna_sense_frac=kappa,
            prior_overdispersion=overdispersion_for_beta(cfgc.gdna_strand_prior_alpha_beta),
            prior_weight=cfgc.gdna_strand_prior_weight)
    rs = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
            prior_overdispersion=overdispersion_for_beta(cfgc.rna_strand_prior_alpha_beta),
            prior_weight=cfgc.rna_strand_prior_weight)
    reg = deconv_regions(sub, ra, nd, rna_sense_frac=kappa,
            gdna_strand_overdispersion=gd.gdna_strand_overdispersion,
            rna_strand_overdispersion=rs.rna_strand_overdispersion, n_grid=cfgc.n_grid)
    truth, U = oracle_unspliced_gdna_frac(cond)
    return dict(kappa=kappa, mu_c=nd.count_gdna_frac, var_c=nd.count_gdna_frac_var,
                var_s=np.asarray(reg.gdna_frac_var), contained=np.asarray(sub.contained.n_unspliced)+0.0,
                truth=truth, U=U, nonobs=(sig & EXON) != 0,
                cross=(np.asarray(sub.left.n_unspliced)+np.asarray(sub.right.n_unspliced))+0.0)


CONDS = [("gdna_gdna300_ss_0.99_nrna_none_capture_on", "ss0.99 cap_on"),
         ("gdna_gdna300_ss_0.50_nrna_none_capture_on", "ss0.50 cap_on"),
         ("gdna_gdna300_ss_0.99_nrna_none_capture_off", "ss0.99 cap_off")]
D = {lab: reconstruct(c) for c, lab in CONDS}

print("=== (A) count variance + (B) blend weights, per condition ===")
for lab, d in D.items():
    k = d["kappa"]; wcur = (2 * k - 1) ** 2
    var_s = d["var_s"]; var_c = d["var_c"]
    so = var_s > 0  # strand-observable (posterior ran)
    wf = np.where(so, var_c / np.maximum(var_s + var_c, 1e-12), 0.0)
    print(f"\n{lab}: kappa={k:.3f}  w_current=(2k-1)^2={wcur:.3f}")
    print(f"  count sigma (sqrt var_c): median={np.median(np.sqrt(var_c)):.3f}  "
          f"observable-region Poisson floor median 1/sqrt(N)~{np.median(1/np.sqrt(np.maximum(d['contained'][~d['nonobs']],1))):.3f}")
    print(f"  boundary Poisson floor 1/sqrt(crossing): median={np.median(1/np.sqrt(np.maximum(d['cross'][d['cross']>0],1))):.3f}")
    if so.any():
        print(f"  w_fisher (strand-observable nodes): median={np.median(wf[so]):.3f}  "
              f"IQR=[{np.percentile(wf[so],25):.3f},{np.percentile(wf[so],75):.3f}]   vs w_current={wcur:.3f}")

# (C) accuracy vs confidence — does small sigma_count predict small |mu_count - truth|?
print("\n=== (C) is count CONFIDENCE a proxy for count ACCURACY? (imputed exon nodes, U>30) ===")
from scipy.stats import spearmanr
for lab, d in D.items():
    sel = d["nonobs"] & (d["U"] > 30) & np.isfinite(d["truth"]) & (d["var_c"] > 0)
    if sel.sum() < 10:
        print(f"  {lab}: too few"); continue
    sig_c = np.sqrt(d["var_c"][sel]); err = np.abs(d["mu_c"][sel] - d["truth"][sel])
    rho, _ = spearmanr(sig_c, err)
    print(f"  {lab}: n={sel.sum()}  Spearman(sigma_count, |error|)={rho:+.3f}  "
          f"(want >0 for Fisher-blend to be trustworthy; <0 ⇒ confident-but-biased ⇒ keep (2k-1)^2)")

# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
d = D["ss0.50 cap_on"]
nb = d["nonobs"]; ob = ~nb
ax[0].scatter(np.maximum(d["contained"][ob], 0.5), np.sqrt(d["var_c"][ob]), s=8, alpha=.4, label="observable region", color="#4477aa")
ax[0].scatter(np.maximum(d["contained"][nb], 0.5), np.sqrt(d["var_c"][nb]), s=8, alpha=.4, label="imputed exon", color="#cc3311")
ax[0].set_xscale("log"); ax[0].set_xlabel("region contained count"); ax[0].set_ylabel("sigma_count (sqrt var_c)")
ax[0].set_title("(A) count posterior std by node type\n(ss0.50 cap_on)"); ax[0].legend(fontsize=8)
for lab, dd in D.items():
    k = dd["kappa"]; so = dd["var_s"] > 0
    wf = dd["var_c"][so] / np.maximum(dd["var_s"][so] + dd["var_c"][so], 1e-12)
    if so.any():
        ax[1].hist(wf, bins=30, alpha=.5, label=f"{lab} (w_cur={(2*k-1)**2:.2f})", density=True)
ax[1].set_xlabel("w_fisher = var_c/(var_s+var_c)  [strand-observable nodes]")
ax[1].set_title("(B) Fisher blend weight vs current (2k-1)^2\n(vertical lines = current w)"); ax[1].legend(fontsize=7)
for lab, dd in D.items():
    k = dd["kappa"]; ax[1].axvline((2*k-1)**2, ls="--", alpha=.6)
d = D["ss0.50 cap_on"]
sel = d["nonobs"] & (d["U"] > 30) & np.isfinite(d["truth"]) & (d["var_c"] > 0)
ax[2].scatter(np.sqrt(d["var_c"][sel]), np.abs(d["mu_c"][sel] - d["truth"][sel]), s=10, alpha=.4, color="#228833")
ax[2].set_xlabel("sigma_count (confidence: small=confident)"); ax[2].set_ylabel("|mu_count - truth| (error)")
ax[2].set_title("(C) does confidence predict accuracy?\n(ss0.50 cap_on imputed nodes)")
fig.tight_layout(); fig.savefig("/tmp/count_variance_intuition.png", dpi=110)
print("\nsaved /tmp/count_variance_intuition.png")
