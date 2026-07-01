"""Flagship gDNA leak — node-level root cause. Re-run the flagship calibration, capture the per-node
belief + class, and aggregate gDNA-assigned by node class. Also pull per-region TRUTH gDNA from the
oracle BAM read names (parse_origin) to compute the per-class under-call directly."""
from pathlib import Path
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
from rigel.calibration.signature import (
    BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG, TS_AMBIG, TS_POS, TS_NEG, TS_NONE,
)
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
belief = bp.node_sweep(chain, statics, geom, b0, ra, bsub, rna_sense_frac=kappa, n_grid=60, _capture=cap)

# --- per-REGION calibration state ---
R = ra.signature.shape[0]
sig = np.asarray(ra.signature).astype(np.int64)
sc = np.asarray(ra.strand_class)
EXON = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
INTRON = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
rtype = np.where(EXON, "exon", np.where(INTRON, "intron", "intergenic"))
# region nodes on the chain → region index
kind = np.asarray(chain.kind); ridx = np.asarray(chain.ref_idx)
reg_node = kind == REGION
reg_of = ridx[reg_node]
fg = np.zeros(R); fgloc = np.zeros(R); fgstr = np.zeros(R)
fg[reg_of] = np.asarray(belief.f_g)[reg_node]
fgloc[reg_of] = np.asarray(cap["fg_loc"])[reg_node]
fgstr[reg_of] = np.asarray(cap["fg_strand"])[reg_node]
mass_u = np.asarray(sub.contained.mass_unspliced)
gdna_assigned = fg * mass_u   # calibration's gDNA mass per region (contained)

def cls(i):
    if rtype[i] == "exon":
        return "exon_AMBIG" if sc[i] == TS_AMBIG else ("exon_SS" if sc[i] in (TS_POS, TS_NEG) else "exon_NONE")
    return rtype[i]

classes = np.array([cls(i) for i in range(R)])

# --- per-REGION TRUTH gDNA from the oracle BAM: bin each gDNA FRAGMENT (full span via mate+TLEN) as
#     CONTAINED in region j iff [leftmost, leftmost+|tlen|) ⊆ [start_j, end_j) — matches the accumulator's
#     contained deposit (so it is comparable to mass_unspl, the contained mass). Crossing fragments → boundary. ---
import collections
starts = np.asarray(ra.start); ends = np.asarray(ra.end); rid = np.asarray(ra.ref_id)
truth_g = np.zeros(R)          # contained gDNA fragments per region
by_ref = collections.defaultdict(list)
for i in range(R):
    by_ref[int(rid[i])].append(i)
for k in list(by_ref):
    arr = sorted(by_ref[k], key=lambda i: starts[i])
    by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]), np.array(arr))
name2ref = {n: i for i, n in enumerate(idx.ref_names)}
bcross = {}  # (region_left, region_right) -> crossing gDNA fragment count
bam = pysam.AlignmentFile(str(F / "sim_oracle.bam"), "rb")
n_g = n_contained = 0
for rd in bam.fetch(until_eof=True):
    if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1 or rd.mate_is_unmapped:
        continue
    if "gdna" not in rd.query_name.lower():
        continue
    n_g += 1
    rf = name2ref.get(bam.get_reference_name(rd.reference_id))
    if rf is None or rf not in by_ref:
        continue
    tl = abs(rd.template_length)
    if tl == 0:
        continue
    lo = min(rd.reference_start, rd.next_reference_start)
    hi = lo + tl
    s, e, ii = by_ref[rf]
    j = np.searchsorted(s, lo, side="right") - 1
    if 0 <= j < len(ii) and s[j] <= lo and hi <= e[j]:   # fully CONTAINED
        truth_g[ii[j]] += 1.0
        n_contained += 1
    elif 0 <= j < len(ii) - 1 and s[j] <= lo < e[j] and s[j + 1] <= hi <= e[j + 1]:
        # crosses exactly the j↔j+1 boundary → record on both flank regions (by region pair)
        bcross[(int(ii[j]), int(ii[j + 1]))] = bcross.get((int(ii[j]), int(ii[j + 1])), 0) + 1
bam.close()
print(f"oracle gDNA fragments: {n_g:,};  contained-in-a-region: {n_contained:,} "
      f"(rest cross boundaries / land outside the region partition)")

# CONTAINED-mass comparison: gdna_assigned = f_g·mass_u; truth f_g = contained gDNA / contained mass.
# fg_loc = message-free anchor (strand+ê); f_g = final (anchor⊗messages). So:
#   anchor error   = truth_fg − fg_loc   (the ê / global-hyperprior gap)
#   message error  = fg_loc − f_g        (how far the messages then drag it, + = drag toward RNA)
print(f"\n{'class':<11}{'n_reg':>6}{'mass_u':>11}{'g_assign':>10}{'g_TRUTH':>9}{'under':>9}"
      f"{'truth_fg':>9}{'fg_loc':>8}{'f_g_fin':>8}{'anchorErr':>10}{'msgDrag':>8}")
for c in ["exon_AMBIG", "exon_SS", "exon_NONE", "intron", "intergenic"]:
    m = classes == c
    if not m.any():
        continue
    w = mass_u[m] + 1e-9
    ga, gt = gdna_assigned[m].sum(), truth_g[m].sum()
    tfg = gt / max(mass_u[m].sum(), 1e-9)           # truth gDNA fraction of contained mass
    afg = np.average(fg[m], weights=w); aloc = np.average(fgloc[m], weights=w)
    print(f"{c:<11}{int(m.sum()):>6}{mass_u[m].sum():>11.0f}{ga:>10.0f}{gt:>9.0f}{ga-gt:>9.0f}"
          f"{tfg:>9.3f}{aloc:>8.3f}{afg:>8.3f}{tfg-aloc:>10.3f}{aloc-afg:>8.3f}")
tot_assigned = gdna_assigned.sum(); tot_truth = truth_g.sum()
print(f"{'TOTAL':<11}{'':>6}{mass_u.sum():>11.0f}{tot_assigned:>10.0f}{tot_truth:>9.0f}{tot_assigned-tot_truth:>9.0f}")
print("  anchorErr>0 ⇒ ê/global under-calls gDNA (anchor too low);  msgDrag>0 ⇒ messages pull f_g toward RNA")
print(f"\nenrich_w (ê applied strength) = {cap['enrich_w']:.3f}   (0 ⇒ ê collapses to flat ρ_global)")
print(f"rho_global = {cap['rho_global']:.4f}")

# --- BOUNDARY nodes: same AMBIG signature? f_g + assigned vs crossing-gDNA truth, by flank class ---
from rigel.calibration.node_chain import BOUNDARY
bnode = kind == BOUNDARY
b_of = ridx[bnode]                       # boundary index per boundary chain-node
fg_node = np.asarray(belief.f_g)
ml = np.asarray(geom.mass_left); mr = np.asarray(geom.mass_right)
blr = np.asarray(bsub.left_region); brr = np.asarray(bsub.right_region)
rows = []
chain_b_nodes = np.where(bnode)[0]
for nidx in chain_b_nodes:
    b = int(ridx[nidx]); lr = int(blr[b]); rr = int(brr[b])
    mass = ml[nidx] + mr[nidx]
    if mass <= 0:
        continue
    lc = classes[lr] if lr >= 0 else "edge"; rc = classes[rr] if rr >= 0 else "edge"
    if "exon_AMBIG" in (lc, rc):
        cl = "AMBIG-adj"
    elif "exon_SS" in (lc, rc):
        cl = "SSexon-adj"
    elif "exon_NONE" in (lc, rc):
        cl = "exonNONE-adj"
    else:
        cl = "intron/intergenic"
    truth = bcross.get((lr, rr), 0) + bcross.get((rr, lr), 0)
    rows.append((cl, mass, fg_node[nidx] * mass, truth, fg_node[nidx]))
import pandas as _pd
B = _pd.DataFrame(rows, columns=["cl", "mass", "assigned", "truth", "fg"])
print(f"\n=== BOUNDARY nodes by flank class (crossing gDNA) ===")
print(f"{'flank class':<18}{'n':>6}{'mass':>11}{'g_assign':>10}{'g_TRUTH':>9}{'under':>9}{'mean_fg':>9}")
for cl, g in B.groupby("cl"):
    print(f"{cl:<18}{len(g):>6}{g['mass'].sum():>11.0f}{g['assigned'].sum():>10.0f}{g['truth'].sum():>9.0f}"
          f"{g['assigned'].sum()-g['truth'].sum():>9.0f}{np.average(g['fg'],weights=g['mass']+1e-9):>9.3f}")
print(f"{'TOTAL(boundary)':<18}{len(B):>6}{B['mass'].sum():>11.0f}{B['assigned'].sum():>10.0f}{B['truth'].sum():>9.0f}"
      f"{B['assigned'].sum()-B['truth'].sum():>9.0f}")
