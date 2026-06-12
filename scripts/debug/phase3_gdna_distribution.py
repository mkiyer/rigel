"""Phase-3 gDNA-distribution forensic: is the per-region gDNA mass concentrated on the on-target exons?

The flagship leak is dominated (78%) by the EM under-recovering gDNA, because the gDNA effective length
(the IPR `(Σg)²/Σ(g²/L)`) is ≈ the full locus span instead of contracting to the on-target exons. The
IPR squares the per-region gDNA mass `g`, so a small UNDER-concentration (gDNA spread to introns instead
of concentrated on exons) is super-linearly amplified into a too-long eff_len → the gDNA competes at
too-low a density → leak. This tool asks the root question: after the boundary-flux transport, is the
per-region gDNA mass concentrated on the exon-signature regions (matching where capture puts it), or
spread to the introns? It reports the gDNA by region class at each stage: contained, +sides,
+transport, vs the oracle.

Usage:  python scripts/debug/phase3_gdna_distribution.py [condition]
"""
import dataclasses
import sys

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.priors import _transport_boundary_flux
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS,
)
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

from phase3_flagship_debug import SUITE, oracle_contained_unspliced  # noqa: E402

EXON = BIT_EXON_POS | BIT_EXON_NEG
INTRON = BIT_INTRON_POS | BIT_INTRON_NEG


def region_class(sig: int) -> str:
    if sig & EXON:
        return "EXON"
    if sig & INTRON:
        return "INTRON"
    return "INTERGENIC"


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    r = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)

    # gDNA mass at each stage, per region
    contained = np.asarray(r.mass_gdna_contained)
    sides = np.asarray(r.mass_gdna_left) + np.asarray(r.mass_gdna_right)
    transported = _transport_boundary_flux(
        r.mass_gdna_contained, r.mass_gdna_left, r.mass_gdna_right,
        np.asarray(ra.region_size_bp, dtype=np.float64), r.gdna_boundary_len, np.asarray(ra.ref_id),
    )

    # oracle gDNA per region (contained-unspliced, by start)
    df = idx.region_df
    starts = {ref: gg["start"].to_numpy() for ref, gg in df.groupby("ref_name")}
    ids = {ref: gg["region_id"].to_numpy() for ref, gg in df.groupby("ref_name")}
    og, om, on = oracle_contained_unspliced(bam, starts, ids, len(df))

    sig = np.asarray(ra.signature)
    cls = np.array([region_class(int(s)) for s in sig])
    print(f"\n===== {cond}: gDNA mass by region class (where does the gDNA land?) =====")
    print(f"{'class':>11}{'n':>6}{'oracle_gdna':>12}{'contained':>12}{'+sides':>12}{'transported':>13}")
    for c in ("EXON", "INTRON", "INTERGENIC"):
        m = cls == c
        print(f"{c:>11}{int(m.sum()):>6}{og[m].sum():>12.0f}{contained[m].sum():>12.0f}"
              f"{(contained+sides)[m].sum():>12.0f}{transported[m].sum():>13.0f}")
    tot = transported.sum()
    print(f"{'TOTAL':>11}{len(sig):>6}{og.sum():>12.0f}{contained.sum():>12.0f}"
          f"{(contained+sides).sum():>12.0f}{transported.sum():>13.0f}")
    # concentration: fraction of transported gDNA on EXON regions vs oracle fraction on EXON
    ex = cls == "EXON"
    print(f"\n  gDNA fraction on EXON regions:  oracle={og[ex].sum()/max(og.sum(),1):.2f}  "
          f"transported={transported[ex].sum()/max(tot,1):.2f}")
    print(f"  (if transported << oracle on exons, the gDNA is under-concentrated → eff_len won't contract)")


if __name__ == "__main__":
    main()
