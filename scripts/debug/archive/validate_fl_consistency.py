"""Validate the FL-consistency fix on REAL data (calibration level, undiluted by the blend).

On the gdna_shortfl condition (gDNA FL≈100 vs RNA FL≈250 — a large gap) it reconstructs, per eligible
splice-junction exon region, the BARE boundary density fraction f_b and the CORRECTED region count
fraction (density_frac_to_count_frac with the region eff-lengths), and compares BOTH to the ORACLE
gDNA fraction of the region's contained unspliced mass (from read-origin names). The fix should bring
the estimate closer to the oracle, most at short exons. (This is the calibration-level effect; the
pipeline output additionally down-weights it by the blend w=(2κ−1)² on stranded nodes.)
"""
import numpy as np
import pysam
from dataclasses import replace as dcr

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.splice_junction import (
    splice_junction_eligibility, boundary_gdna_fraction, density_frac_to_count_frac,
)
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_shortfl_5mb"
COND = "gdna_gdna400_ss_0.99_nrna_none_capture_on"
EXON = BIT_EXON_POS | BIT_EXON_NEG
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
df = idx.region_df
starts = {ref: g["start"].to_numpy() for ref, g in df.groupby("ref_name")}
ids = {ref: g["region_id"].to_numpy() for ref, g in df.groupby("ref_name")}
R = len(df)


def oracle_region_gdna_frac():
    """Per region: gDNA fraction of CONTAINED UNSPLICED fragments (gdna* vs others), from origin names."""
    ug, um = np.zeros(R), np.zeros(R)
    bam = pysam.AlignmentFile(f"{SUITE}/{COND}/sim_oracle.bam", "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
            continue
        if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
            continue  # spliced — not contained-unspliced
        ref = r.reference_name
        if ref not in starts:
            continue
        i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
        if i < 0:
            continue
        rid = int(ids[ref][i])
        (ug if r.query_name.startswith("gdna") else um)[rid] += 1.0
    tot = ug + um
    return np.where(tot > 0, ug / np.maximum(tot, 1e-9), np.nan), tot


bam = f"{SUITE}/{COND}/sim_oracle.bam"
cfg = PipelineConfig()
scan = dcr(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
sub = CalibrationSubstrate.from_payload(pl, ra)
fl_mean_g = boundary_eff_length(fl.gdna_pmf)
fl_mean_r = boundary_eff_length(fl.rna_pmf)
eg_region = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
er_region = region_eff_length(ra.region_size_bp, fl.rna_pmf)
print(f"gDNA fl_mean={fl_mean_g:.1f}  RNA fl_mean={fl_mean_r:.1f}  (gap={fl_mean_g-fl_mean_r:+.0f})")

sig = np.asarray(ra.signature)
ref_id = np.asarray(ra.ref_id)
L = np.asarray(ra.region_size_bp)
oracle, n_contained = oracle_region_gdna_frac()
left, right = sub.left, sub.right
lu = (left.n_unspliced_pos + left.n_unspliced_neg).astype(float)
ls = (left.n_spliced_sense + left.n_spliced_antisense).astype(float)
ru = (right.n_unspliced_pos + right.n_unspliced_neg).astype(float)
rs = (right.n_spliced_sense + right.n_spliced_antisense).astype(float)


def anchor_fb(unspl, spl):
    if spl <= 0.0:
        return None
    f = boundary_gdna_fraction(unspliced_gdna=unspl, unspliced_rna=0.0, spliced=spl,
                               eff_gdna=fl_mean_g, eff_rna=fl_mean_r)
    return None if np.isnan(f) else f


rows = []  # (L, oracle, bare, corrected)
for i in range(R):
    anchors = []
    if i > 0 and ref_id[i] == ref_id[i - 1]:
        e = splice_junction_eligibility(int(sig[i - 1]), int(sig[i]))
        if e is not None and e.exon_side == "R":
            a = anchor_fb(lu[i], ls[i])
            if a is not None:
                anchors.append(a)
    if i < R - 1 and ref_id[i] == ref_id[i + 1]:
        e = splice_junction_eligibility(int(sig[i]), int(sig[i + 1]))
        if e is not None and e.exon_side == "L":
            a = anchor_fb(ru[i], rs[i])
            if a is not None:
                anchors.append(a)
    if anchors and np.isfinite(oracle[i]) and n_contained[i] >= 10:
        fb = float(np.mean(anchors))
        corr = density_frac_to_count_frac(fb, eg_region[i], er_region[i])
        rows.append((L[i], oracle[i], fb, corr))

rows = np.array(rows)
print(f"\neligible exon regions with oracle & >=10 contained: n={len(rows)}")
print(f"{'exon L band':>14} {'n':>4} {'oracle gf':>10} {'bare f_b':>9} {'corrected':>10} "
      f"{'|bare-or|':>10} {'|corr-or|':>10}")
for lo, hi in [(0, 200), (200, 400), (400, 800), (800, 1e9)]:
    m = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
    if m.sum() == 0:
        continue
    sub_ = rows[m]
    err_bare = np.abs(sub_[:, 2] - sub_[:, 1]).mean()
    err_corr = np.abs(sub_[:, 3] - sub_[:, 1]).mean()
    band = f"{lo}-{hi if hi < 1e9 else '+'}"
    print(f"{band:>14} {int(m.sum()):>4} {sub_[:,1].mean():>10.3f} {sub_[:,2].mean():>9.3f} "
          f"{sub_[:,3].mean():>10.3f} {err_bare:>10.3f} {err_corr:>10.3f}")
eb = np.abs(rows[:, 2] - rows[:, 1]).mean()
ec = np.abs(rows[:, 3] - rows[:, 1]).mean()
print(f"\nALL eligible: mean |bare − oracle| = {eb:.3f}   mean |corrected − oracle| = {ec:.3f}   "
      f"({'IMPROVED' if ec < eb else 'WORSE'} by {eb - ec:+.3f})")
