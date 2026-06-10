"""Phase 4-mean make-or-break: can the SPLICED channel debias the exon gDNA estimate?

Hypothesis: an exon's unspliced mass = gDNA + mature-exon-body RNA, and the mature-exon-body reads
are predictable from the (clean, unconflated) SPLICED junction reads via a structural ratio
r = U_mat / S.  Then gDNA = U_total − r·S, using the exon's OWN (capture-enriched) counts.

CONSTRAINT (must): r is only defined where splice junctions exist — i.e. exon regions that abut a
true splice junction (spliced reads present, S > 0).  Single-exon / junction-free exons are excluded
(no denominator).

Oracle read origins (nrna_none): 'gdna*' = gDNA, 'GENE*' = mature RNA.  Spliced (N in CIGAR) = mature
junction reads; unspliced = gDNA + mature-exon-body.  We check (1) is r predictable (vs region
length), and (2) does U_total − r̂·S recover the oracle exon gDNA fraction (biased ≈0.41 → truth ≈0.91)?
"""
import numpy as np, pysam
from dataclasses import replace as _dc_replace
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.density_model import count_observable_masks
from scipy.stats import spearmanr

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
COND = "gdna_gdna1000_ss_0.99_nrna_none_capture_on"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
df = idx.region_df
EXON_BITS = 0b0011
sig = df["signature"].to_numpy(); observable = (sig & EXON_BITS) == 0
length = (df["end"] - df["start"]).to_numpy().astype(float)
starts = {}; ids = {}
for ref, g in df.groupby("ref_name"):
    starts[ref] = g["start"].to_numpy(); ids[ref] = g["region_id"].to_numpy()
R = len(df)
u_gdna = np.zeros(R); u_mat = np.zeros(R); spliced = np.zeros(R)

bam = f"{SUITE}/{COND}/sim_oracle.bam"
for r in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True):
    if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
        continue
    ref = r.reference_name
    if ref not in starts:
        continue
    i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
    if i < 0:
        continue
    rid = int(ids[ref][i])
    is_spliced = bool(r.cigartuples) and any(op == 3 for op, _ in r.cigartuples)
    is_gdna = r.query_name.startswith("gdna")
    if is_spliced:
        spliced[rid] += 1.0          # mature junction read (gDNA/nascent never splice)
    elif is_gdna:
        u_gdna[rid] += 1.0           # unspliced gDNA
    else:
        u_mat[rid] += 1.0            # unspliced mature exon-body

U = u_gdna + u_mat                   # total unspliced
# Junction-bearing exon regions only (the MUST constraint): non-observable, with spliced reads.
sel = (~observable) & (spliced > 0) & (u_mat > 0) & (U > 30)
r_emp = u_mat[sel] / spliced[sel]    # empirical mature within-exon / junction ratio
L = length[sel]
true_frac = u_gdna[sel] / U[sel]
w = U[sel]                           # mass-weight by unspliced count

print(f"=== {COND} ===")
print(f"junction-bearing exon regions: {sel.sum()} (of {(~observable).sum()} non-observable)")
print(f"empirical r = U_mat/S:  median={np.median(r_emp):.3f}  IQR=[{np.percentile(r_emp,25):.3f},"
      f"{np.percentile(r_emp,75):.3f}]  CV={r_emp.std()/r_emp.mean():.2f}")
rho_L, _ = spearmanr(r_emp, L)
print(f"  Spearman(r, region_length) = {rho_L:+.3f}   (is r predictable from length?)")
print(f"\ntrue exon gDNA fraction (oracle):     mass-wtd = {np.average(true_frac, weights=w):.3f}")
# Current (biased) count estimate proxy: gDNA frac if we used boundary-imputed density — for context
print(f"  (the boundary-imputed count gave ~0.41 here; truth ~0.91 — diag_imputation_truth.py)")

def recovered(rhat):
    gdna_est = np.maximum(U[sel] - rhat * spliced[sel], 0.0)
    frac = gdna_est / U[sel]
    return np.average(frac, weights=w)

# (a) ORACLE r (per-exon true ratio): algebra sanity — must recover the truth exactly.
print(f"\nrecovery U_total − r̂·S, mass-wtd gDNA fraction:")
print(f"  (a) r̂ = per-exon TRUE r (sanity):       {recovered(r_emp):.3f}  (must ≈ true)")
# (b) CONSTANT r̂ = median empirical r (no structure): how far off?
print(f"  (b) r̂ = constant median(r):             {recovered(np.median(r_emp)):.3f}")
# (c) STRUCTURAL r̂ from a length fit r = a·L + b (predicted, leave-one-out-ish via global fit):
A = np.vstack([L, np.ones_like(L)]).T
coef, *_ = np.linalg.lstsq(A, r_emp, rcond=None)
rhat_L = A @ coef
print(f"  (c) r̂ = length fit (a·L+b), a={coef[0]:.2e} b={coef[1]:.2f}: {recovered(rhat_L):.3f}")
print(f"      Spearman(true_frac, recovered_with_length_fit) per-exon = "
      f"{spearmanr((U[sel]-rhat_L*spliced[sel])/U[sel], true_frac)[0]:+.3f}")
