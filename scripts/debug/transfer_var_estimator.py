"""Break the circularity: estimate the de-conflated gDNA transfer variance from OBSERVABLE anchor pairs
(no solve), then VALIDATE against the oracle that it generalizes to the unobservable RNA-present pairs
(the homogeneity assumption). Also show the total-density σ²_imp inflation on real data.

Estimator uses ONLY observables: total density, per-node flanking SPLICED density (RNA proxy), node type.
Oracle is used ONLY to validate (compute the TRUE gDNA transfer variance on all pairs, incl. RNA-present).
Condition: gdna_gdna300_ss_0.50_nrna_none_capture_on (the flagship unstranded+capture).
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
from rigel.calibration.effective_length import region_eff_length
from rigel.calibration.node_chain import build_node_chain, REGION, BOUNDARY
from rigel.calibration.bp_solver import build_node_geometry
from rigel.splice import SpliceType
from oracle import OracleTruth

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = "gdna_gdna300_ss_0.50_nrna_none_capture_on"
bam = str(SUITE / COND / "sim_oracle.bam"); idx = TranscriptIndex.load(str(SUITE / "rigel_index"))
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id); cfg = CalibrationConfig()
st, sm, flm, buf, pl = scan_and_buffer(bam, idx, dc(BamScanConfig(), sj_strand_tag=_native_detect_sj_tag(bam)))
sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
fl = build_fl_models(global_counts=flm.global_model.counts, rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts, gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geom = build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
orc = OracleTruth.from_bam(bam, idx, PipelineConfig(), Path("/tmp/tv_split"), COND, full_payload=pl, boundary_mass_tol=0.5)
p = orc.region_pools()
n = ra.n_regions
rel = np.maximum(region_eff_length(ra.region_size_bp, fl.gdna_pmf), 1e-9)
M = np.asarray(sub.contained.mass_unspliced, float)
total_dens = M / rel                                  # OBSERVABLE total density
true_g = p["gdna_pos"] + p["gdna_neg"]
true_gdna_dens = true_g / rel                         # ORACLE gDNA density (validation only)
true_rna = (p["mat_uns_pos"] + p["mat_uns_neg"] + p["nas_uns_pos"] + p["nas_uns_neg"])
true_rna_frac = np.where(M > 0, true_rna / np.maximum(M, 1e-9), 0.0)

# per-region flanking SPLICED density (OBSERVABLE RNA proxy) — spliced mass at its two boundary neighbours
kind = np.asarray(chain.kind); cref = np.asarray(chain.ref_idx, np.int64)
reg_node = np.full(n, -1, np.int64); rmask = kind == REGION; reg_node[cref[rmask]] = np.where(rmask)[0]
left = np.asarray(chain.left); right = np.asarray(chain.right)
spl = (np.asarray(geom.spliced_pos_left) + np.asarray(geom.spliced_neg_left)
       + np.asarray(geom.spliced_pos_right) + np.asarray(geom.spliced_neg_right))
esp = np.maximum(np.asarray(geom.eff_spl_left) + np.asarray(geom.eff_spl_right), 1e-9)
spl_dens_node = spl / esp
flank_spl = np.zeros(n)
for r in range(n):
    ndc = reg_node[r]
    if ndc < 0:
        continue
    s = 0.0
    for nb in (left[ndc], right[ndc]):
        if nb >= 0 and kind[nb] == BOUNDARY:
            s = max(s, spl_dens_node[nb])
    flank_spl[r] = s

ref = np.asarray(ra.ref_id) if hasattr(ra, "ref_id") else idx.region_df["ref_name"].to_numpy()

# adjacent same-ref region pairs with real mass on both sides
FLOOR = 1e-3
pairs = []
for r in range(n - 1):
    if reg_node[r] < 0 or reg_node[r + 1] < 0:
        continue
    if ref[r] != ref[r + 1] or M[r] <= 0 or M[r + 1] <= 0:
        continue
    pairs.append((r, r + 1))
pairs = np.array(pairs)
a, b = pairs[:, 0], pairs[:, 1]

def v(x):  # variance of Δlog across a boolean-selected pair set
    return float(np.var(x)) if len(x) else float("nan")

dlog_total = np.log(np.maximum(total_dens[a], FLOOR)) - np.log(np.maximum(total_dens[b], FLOOR))
dlog_gtrue = np.log(np.maximum(true_gdna_dens[a], FLOOR)) - np.log(np.maximum(true_gdna_dens[b], FLOOR))

# OBSERVABLE RNA-free gate (spliced density below a low threshold on BOTH sides)
thr = np.percentile(flank_spl[flank_spl > 0], 25) if (flank_spl > 0).any() else 0.1
obs_rnafree = (flank_spl < max(thr, 0.05))
# TRUE RNA-free (oracle) — to check the gate
true_rnafree = true_rna_frac < 0.05

print(f"=== transfer-variance estimator — {COND} ===")
print(f"pairs (adjacent same-ref, mass>0): {len(pairs)};  observable RNA-free thr(spliced dens)={max(thr,0.05):.3f}")
print("\ngate check: observable-RNA-free vs true-RNA-free regions —")
print(f"  observable RNA-free: {int(obs_rnafree.sum())}  true RNA-free: {int(true_rnafree.sum())}  "
      f"agreement(obs⊆true): {float(np.mean(true_rnafree[obs_rnafree])):.2%}")

both_obs = obs_rnafree[a] & obs_rnafree[b]     # ANCHOR pairs: both sides observably RNA-free
both_true = true_rnafree[a] & true_rnafree[b]
print("\n=== VARIANCE of Δlog density across adjacent region pairs ===")
print(f"  TOTAL density, ALL pairs          [current σ²_imp basis]  = {v(dlog_total):.3f}   (n={len(a)})")
print(f"  TRUE gDNA density, ALL pairs       [the target variance]   = {v(dlog_gtrue):.3f}")
print(f"  ANCHOR: TOTAL density on OBS-RNA-free pairs [ESTIMATOR]     = {v(dlog_total[both_obs]):.3f}   (n={int(both_obs.sum())})")
print(f"  (cross-check) TRUE gDNA on TRUE-RNA-free pairs             = {v(dlog_gtrue[both_true]):.3f}   (n={int(both_true.sum())})")
print(f"\n  ⇒ HOMOGENEITY: does the observable anchor estimate ({v(dlog_total[both_obs]):.3f}) match the true gDNA"
      f"\n    transfer variance over ALL pairs ({v(dlog_gtrue):.3f})?  ratio = {v(dlog_total[both_obs])/max(v(dlog_gtrue),1e-9):.2f}×")
print(f"  ⇒ INFLATION: total-density σ²_imp ({v(dlog_total):.3f}) / true gDNA σ² ({v(dlog_gtrue):.3f}) "
      f"= {v(dlog_total)/max(v(dlog_gtrue),1e-9):.1f}×  (the composition conflation)")

# stratify the TRUE gDNA transfer variance by enrichment regime (both-enriched / mixed / both-depleted)
en = total_dens > 1.0
strat = {"both-enriched": en[a] & en[b], "mixed": en[a] ^ en[b], "both-depleted": (~en[a]) & (~en[b])}
print("\n=== TRUE gDNA transfer variance by enrichment regime (observable stratification by total density) ===")
for name, msk in strat.items():
    print(f"  {name:14}: n={int(msk.sum()):5d}  Var(Δlog gDNA)={v(dlog_gtrue[msk]):.3f}  "
          f"Var(Δlog total)={v(dlog_total[msk]):.3f}")
