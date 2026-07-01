"""Gating measurement for the AMBIG anchor fix (defect c): does the ê(z) predictor SATURATE at high
enrichment? Compare candidate predictors against the oracle TRUE contained gDNA density on exon regions:
  z_flux      = the CURRENT ê predictor (flanking clean-boundary crossing density; cap['z_enrich'])
  rho_contain = the node's own contained density M/E_gdna  (the §2.4 'contained interior' candidate)
  ehat_pred   = what ê(z) actually predicts
vs true_rho_g = oracle contained gDNA / E_gdna. If z_flux saturates (flat at high true density) but
rho_contain stays monotone, the PREDICTOR is the bottleneck (swap it); else the teacher-weighting is."""
from pathlib import Path
import collections
import numpy as np
import pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration import bp_solver as bp
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, TS_AMBIG, TS_POS, TS_NEG
from rigel.splice import SpliceType

ROOT = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
F = ROOT / "gdna_gdna300_ss_0.99_nrna_none_capture_on"
idx = TranscriptIndex.load(str(ROOT / "rigel_index"))
st_, sm, flm, buf, pl = scan_and_buffer(str(F / "sim_oracle.bam"), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
statics = bp.build_node_statics(chain, sub, bsub, ra)
kappa = float(fit_strand_balance(sm).rna_sense_frac)
b0 = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
cap = {}
bp.node_sweep(chain, statics, geom, b0, ra, bsub, rna_sense_frac=kappa, n_grid=60, _capture=cap)

R = ra.signature.shape[0]
sig = np.asarray(ra.signature).astype(np.int64); scl = np.asarray(ra.strand_class)
EXON = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
kind = np.asarray(chain.kind); ridx = np.asarray(chain.ref_idx); reg_node = kind == REGION
reg_of = ridx[reg_node]
EG = np.asarray(geom.eff_gdna_left)        # region contained gDNA eff-len at region node
mass_u = np.asarray(sub.contained.mass_unspliced)
# map per-node quantities to region index
z_reg = np.full(R, np.nan); eg_reg = np.zeros(R); ehat_reg = np.full(R, np.nan)
z_node = np.asarray(cap["z_enrich"]); ehat = cap["ehat"]
z_reg[reg_of] = z_node[reg_node]
eg_reg[reg_of] = EG[reg_node]
fin = np.isfinite(z_reg)
ehat_reg[fin] = ehat.predict(z_reg[fin])

# oracle TRUE contained gDNA per region (fragment fully inside region, via mate+TLEN)
starts = np.asarray(ra.start); ends = np.asarray(ra.end); rid = np.asarray(ra.ref_id)
truth_g = np.zeros(R)
by_ref = collections.defaultdict(list)
for i in range(R):
    by_ref[int(rid[i])].append(i)
for k in list(by_ref):
    arr = sorted(by_ref[k], key=lambda i: starts[i])
    by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]), np.array(arr))
name2ref = {n: i for i, n in enumerate(idx.ref_names)}
bam = pysam.AlignmentFile(str(F / "sim_oracle.bam"), "rb")
for rd in bam.fetch(until_eof=True):
    if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1 or rd.mate_is_unmapped:
        continue
    if "gdna" not in rd.query_name.lower():
        continue
    rf = name2ref.get(bam.get_reference_name(rd.reference_id))
    if rf is None or rf not in by_ref:
        continue
    tl = abs(rd.template_length)
    if tl == 0:
        continue
    lo = min(rd.reference_start, rd.next_reference_start); hi = lo + tl
    s, e, ii = by_ref[rf]
    j = np.searchsorted(s, lo, side="right") - 1
    if 0 <= j < len(ii) and s[j] <= lo and hi <= e[j]:
        truth_g[ii[j]] += 1.0
bam.close()

eg_safe = np.maximum(eg_reg, 1e-9)
true_rho = truth_g / eg_safe                 # ORACLE contained gDNA density
contain_rho = mass_u / eg_safe               # candidate predictor: total contained density M/E

def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    if m.sum() < 5:
        return float("nan"), int(m.sum())
    return float(np.corrcoef(np.log(a[m]), np.log(b[m]))[0, 1]), int(m.sum())

for label, cls in [("SINGLE-STRAND exons (the teacher set)", (scl == TS_POS) | (scl == TS_NEG)),
                   ("AMBIG exons (the leak set)", scl == TS_AMBIG)]:
    m = EXON & cls & (mass_u > 0) & np.isfinite(z_reg)
    print(f"\n===== {label}:  n={int(m.sum())} =====")
    cz, nz = corr(z_reg, true_rho); cc, nc = corr(contain_rho, true_rho)
    print(f"  corr(log z_flux,    log true_rho_g) = {cz:.3f}   <- the CURRENT predictor")
    print(f"  corr(log rho_contain, log true_rho_g) = {cc:.3f}   <- candidate (contained density)")
    # saturation: bin by TRUE density quartiles, show mean predictor + ê prediction per bin
    tr = true_rho[m]
    qs = np.quantile(tr[tr > 0], [0, .25, .5, .75, 1.0]) if (tr > 0).any() else None
    if qs is not None:
        print(f"  {'true_rho bin':<22}{'n':>5}{'mean_true':>11}{'mean_z_flux':>12}{'mean_contain':>13}{'mean_ehat':>11}")
        zz, cc2, ee, tt = z_reg[m], contain_rho[m], ehat_reg[m], true_rho[m]
        for a, b in zip(qs[:-1], qs[1:]):
            bm = (tt >= a) & (tt <= b) if b == qs[-1] else (tt >= a) & (tt < b)
            if bm.sum() == 0:
                continue
            print(f"  [{a:8.2f},{b:8.2f}]{'':<3}{int(bm.sum()):>5}{np.mean(tt[bm]):>11.2f}"
                  f"{np.nanmean(zz[bm]):>12.2f}{np.nanmean(cc2[bm]):>13.2f}{np.nanmean(ee[bm]):>11.2f}")
print("\nSATURATION TEST: if mean_z_flux flattens across rising true_rho bins but mean_contain keeps rising,")
print("the PREDICTOR z is the bottleneck (swap to contained density); ê can't lift no matter the teacher.")
