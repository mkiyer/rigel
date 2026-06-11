"""Phase 1 acceptance test: does the count posterior variance σ_g CALIBRATE?

For each non-observable region we have the count module's gDNA-fraction estimate μ_g
(count_gdna_frac), its predicted std σ_g = sqrt(count_gdna_frac_var), and the ORACLE gDNA fraction of
the unspliced mass (gdna* vs GENE* read names). A well-behaved variance predicts error: binning nodes
by σ_g, the realized |μ_g − truth| should INCREASE with σ_g (Spearman > 0), and σ_g should be on the
same scale as the realized RMS error (rough coverage), not orders off.
"""
import numpy as np, pysam
from dataclasses import replace as dcr
from scipy.stats import spearmanr
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
COND = "gdna_gdna1000_ss_0.50_nrna_none_capture_on"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
df = idx.region_df
sig = df["signature"].to_numpy()
EXON = BIT_EXON_POS | BIT_EXON_NEG
starts, ids = {}, {}
for ref, g in df.groupby("ref_name"):
    starts[ref] = g["start"].to_numpy(); ids[ref] = g["region_id"].to_numpy()

R = len(df)
u_gdna = np.zeros(R); u_mat = np.zeros(R)
for r in pysam.AlignmentFile(f"{SUITE}/{COND}/sim_oracle.bam", "rb").fetch(until_eof=True):
    if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
        continue
    ref = r.reference_name
    if ref not in starts:
        continue
    if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
        continue  # spliced = mature junction, not part of the unspliced fraction
    i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
    if i < 0:
        continue
    rid = int(ids[ref][i])
    if r.query_name.startswith("gdna"):
        u_gdna[rid] += 1.0
    else:
        u_mat[rid] += 1.0

bam = f"{SUITE}/{COND}/sim_oracle.bam"
cfg = PipelineConfig(); scan = dcr(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats, sm, flm, buffer, payload = scan_and_buffer(bam, idx, scan)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
sub = CalibrationSubstrate.from_payload(payload, ra)
rel = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
nd = node_gdna_density(sub, ra, rel, boundary_eff_length(fl.gdna_pmf))

U = u_gdna + u_mat
truth = np.where(U > 0, u_gdna / np.maximum(U, 1e-9), np.nan)
mu = nd.count_gdna_frac
sigma = np.sqrt(np.maximum(nd.count_gdna_frac_var, 0.0))
nonobs = (sig & EXON) != 0
sel = nonobs & (U > 30) & np.isfinite(truth)
err = np.abs(mu[sel] - truth[sel])
sg = sigma[sel]

print(f"=== {COND} ===")
print(f"fitted α_impute = {nd.count_rel_var_alpha:.4f}")
print(f"non-observable nodes with truth (U>30): {sel.sum()}")
print(f"realized mean |μ_g − truth| = {err.mean():.3f}  (RMS = {np.sqrt((err**2).mean()):.3f})")
print(f"predicted mean σ_g = {sg.mean():.3f}")
rho, p = spearmanr(sg, err)
print(f"Spearman(σ_g, |error|) = {rho:+.3f}  (p={p:.1e})  — should be POSITIVE (σ predicts error)")
order = np.argsort(sg)
q = np.array_split(order, 4)
print("\nquartile of σ_g:   mean σ_g   mean|error|   RMS error   n")
for k, idxs in enumerate(q):
    print(f"  Q{k+1}: {sg[idxs].mean():.3f}      {err[idxs].mean():.3f}        "
          f"{np.sqrt((err[idxs]**2).mean()):.3f}      {len(idxs)}")
