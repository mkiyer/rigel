"""Verify the user's point: corr(contained_density, true_gDNA_density)=0.995 is INFLATED by the flagship's
3:1 gDNA:RNA ratio. contained = gDNA + RNA, so as the gDNA share shrinks, contained tracks RNA, not gDNA.
Simulate lower contamination by rescaling per-region gDNA by s and recomputing the correlation."""
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
R = ra.signature.shape[0]
sig = np.asarray(ra.signature).astype(np.int64); scl = np.asarray(ra.strand_class)
EXON = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
kind = np.asarray(chain.kind); ridx = np.asarray(chain.ref_idx); reg_node = kind == REGION
EG = np.full(R, 1.0); EG[ridx[reg_node]] = np.asarray(geom.eff_gdna_left)[reg_node]
EG = np.maximum(EG, 1e-9)

# oracle contained counts per region, split gDNA vs RNA (mrna), via mate+TLEN containment
starts = np.asarray(ra.start); ends = np.asarray(ra.end); rid = np.asarray(ra.ref_id)
g_cnt = np.zeros(R); r_cnt = np.zeros(R)
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
    is_g = "gdna" in rd.query_name.lower()
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
        (g_cnt if is_g else r_cnt)[ii[j]] += 1.0
bam.close()

def corr(a, b):
    m = (a > 0) & (b > 0) & np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(np.log(a[m]), np.log(b[m]))[0, 1]) if m.sum() > 5 else float("nan")

for label, cls in [("ALL exons", EXON), ("AMBIG exons", EXON & (scl == TS_AMBIG)),
                   ("single-strand exons", EXON & ((scl == TS_POS) | (scl == TS_NEG)))]:
    m = cls & ((g_cnt + r_cnt) > 0)
    g = g_cnt[m] / EG[m]; r = r_cnt[m] / EG[m]
    print(f"\n===== {label}: n={int(m.sum())},  gDNA:RNA contained = {g_cnt[m].sum():.0f}:{r_cnt[m].sum():.0f}"
          f" ({g_cnt[m].sum()/max(r_cnt[m].sum(),1):.2f}:1) =====")
    print(f"  corr(log contained_density, log RNA_density)  = {corr(g + r, r):.3f}  (contained tracks RNA too)")
    print(f"  {'gDNA scale s':>12}{'gdna:rna':>10}{'corr(contained_s, gDNA)':>26}")
    for s in [3.0, 1.0, 0.3, 0.1, 0.03]:
        contained_s = s * g + r
        gdr = s * g_cnt[m].sum() / max(r_cnt[m].sum(), 1)
        print(f"  {s:>12.2f}{gdr:>9.2f}:1{corr(contained_s, g):>26.3f}")
print("\nIf corr(contained_s, gDNA) falls toward corr(contained,RNA) as s drops, the contained-density")
print("predictor is RATIO-INFLATED and fails at low contamination — confirming the depleted-residual is required.")
