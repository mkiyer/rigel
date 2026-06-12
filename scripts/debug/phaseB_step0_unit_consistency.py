"""Phase B / Step 0: is predicted mature M (from spliced flux) in the SAME accumulator units as U?

The subtraction U − M is only valid if M and the contained count U share units. Both are accumulator
fragment counts, so they should — but the locus-21 per-region check (vs a read1 oracle) showed scatter,
possibly a read1-vs-fragment artifact OR a geometry bug at gene-edge / single-junction exons.

This resolves it rigorously: build a MATURE-ONLY BAM (filter sim_oracle.bam to mrna-origin fragments) and
scan it — the accumulator's contained-unspliced of the mature-only BAM IS the substrate's contained mature,
in exact accumulator units. Compare to M = mature_pos + mature_neg from the FULL scan's spliced flux
(spliced is mature-only anyway). M / substrate_mature ≈ 1 ⇒ units consistent, subtraction valid. A
systematic per-region bias (esp. single-junction / gene-edge exons) ⇒ a geometry fix is needed first.

Usage:  python scripts/debug/phaseB_step0_unit_consistency.py [condition]
"""
import dataclasses
import pathlib
import sys
import tempfile

import numpy as np
import pysam

from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.mature_density import mature_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
_EXON = BIT_EXON_POS | BIT_EXON_NEG


def _mature_only_bam(src, dst):
    """Write a BAM with only mrna-origin reads (origin is per-fragment, same for both mates)."""
    with pysam.AlignmentFile(src, "rb") as bam_in:
        with pysam.AlignmentFile(dst, "wb", template=bam_in) as bam_out:
            for r in bam_in.fetch(until_eof=True):
                if parse_origin(r.query_name).kind == "mrna":
                    bam_out.write(r)


def _scan(bam, idx, cfg):
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, _sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    return flm, pl


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)

    # FULL scan → predicted M from spliced flux.
    flm, pl = _scan(bam, idx, cfg)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    md = mature_density(sub, ra, e_mu, boundary_eff_length(fl.rna_pmf))
    M = md.mature_pos + md.mature_neg

    # MATURE-ONLY scan → substrate contained mature (exact accumulator units).
    tmp = pathlib.Path(tempfile.mkdtemp()) / "mature_only.bam"
    _mature_only_bam(bam, str(tmp))
    _flm2, pl2 = _scan(str(tmp), idx, cfg)
    sub_m = CalibrationSubstrate.from_payload(pl2, ra)
    U_mat = (sub_m.contained.n_unspliced_pos + sub_m.contained.n_unspliced_neg).astype(float)

    sig = np.asarray(ra.signature)
    df = idx.region_df
    rid = df["region_id"].to_numpy()
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()

    # aggregate accuracy on substantial exon regions
    exon = ((sig & _EXON) != 0) & (U_mat > 30.0)
    ratio = M[exon] / np.maximum(U_mat[exon], 1e-9)
    print(f"=== {cond}: M (predicted) vs SUBSTRATE contained mature (mature-only scan) ===")
    print(f"exon regions with substrate mature > 30: n={int(exon.sum())}")
    print(f"  M/U_mature  median={np.median(ratio):.3f}  mean={np.mean(ratio):.3f}  "
          f"p10={np.percentile(ratio,10):.3f}  p90={np.percentile(ratio,90):.3f}")
    print(f"  total M={M[exon].sum():.0f}  total substrate mature={U_mat[exon].sum():.0f}  "
          f"ratio={M[exon].sum()/max(U_mat[exon].sum(),1):.3f}")

    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165) & ((sig & _EXON) != 0))
    print("\n  locus 21 exon regions (rid start  U_mature_substrate  M_pred  ratio  n_junc-hint):")
    for i in sel:
        rr = M[i] / U_mat[i] if U_mat[i] > 0 else float("nan")
        print(f"   {rid[i]:>4}{s[i]:>8}{U_mat[i]:>10.0f}{M[i]:>9.0f}{rr:>7.2f}")
    tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
