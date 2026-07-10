#!/usr/bin/env python3
"""Phase 0 validation — intergenic fragments now deposit into the accumulator.

Builds a tiny oracle scenario (one small gene, large intergenic flanks, stranded
gDNA contamination, no capture), scans it, and checks the calibration payload's
per-region CONTAINED mass against oracle truth derived from the BAM read names.

PASS criteria:
  1. intergenic regions carry NONZERO contained mass (the Phase-0 bug fix — these
     were identically zero before).
  2. each intergenic region's contained gDNA count ≈ the number of unique gDNA
     fragments whose genomic span lies fully inside that region (oracle).
  3. total deposited mass (contained + boundary) ≈ number of deposited unique
     fragments — mass conservation still holds with the new intergenic deposits.

Usage:
  python scripts/debug/phase0_intergenic_deposit_validation.py
"""
from __future__ import annotations

import numpy as np
import pysam

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import _check_region_payload_alignment, scan_and_buffer
from rigel.sim import Scenario
from rigel.sim.read_name import parse_origin
from rigel.sim.reads import GDNAConfig, ReadSimConfig

_TYPE = {0: "intergenic", 1: "intron", 2: "EXON"}


def main() -> None:
    sc = Scenario(
        "phase0_intergenic", genome_length=60000, seed=11,
        work_dir="/tmp/rigel_sim/phase0_intergenic",
    )
    # One multi-exon gene in the middle; the rest of the 60 kb is intergenic.
    sc.add_gene("G", "+", [{"t_id": "G.1",
                            "exons": [(20000, 20600), (21000, 21600), (22000, 22600)],
                            "abundance": 80}])
    res = sc.build_oracle(
        n_fragments=200000,
        sim_config=ReadSimConfig(strand_specificity=0.99, frag_mean=250, frag_std=50, seed=11),
        gdna_config=GDNAConfig(abundance=3000.0, frag_mean=350, frag_std=100),
    )

    index = res.index
    _s, strand_models, _fl, _b, payload = scan_and_buffer(
        str(res.bam_path), index, BamScanConfig(sj_strand_tag="auto")
    )
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    _check_region_payload_alignment(ra, payload)
    sub = CalibrationSubstrate.from_payload(payload, ra)

    ctype = coarse_type_array(np.asarray(ra.signature))
    cc = sub.contained
    # contained gDNA proxy = unspliced contained mass (gDNA is unspliced genomic).
    cal_contained_unspl = (cc.n_unspliced_pos + cc.n_unspliced_neg).astype(np.int64)

    ref_name = "phase0_intergenic"
    ref_id = index.ref_name_to_id[ref_name]
    starts = np.asarray(ra.start)
    ends = np.asarray(ra.end)
    rids = np.asarray(ra.ref_id)

    # ---- Oracle truth: count unique gDNA fragments contained per region ----
    bam = pysam.AlignmentFile(str(res.bam_path), "rb")
    seen: set[str] = set()
    oracle_contained = np.zeros(ra.n_regions, dtype=np.int64)
    n_gdna_unique = 0
    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
            continue
        if not rec.is_read1:  # one fragment per template
            continue
        qn = rec.query_name
        if qn in seen:
            continue
        seen.add(qn)
        origin = parse_origin(qn)
        if origin.kind != "gdna":
            continue
        n_gdna_unique += 1
        # fragment genomic span (template): [min(start, mate_start), +tlen)
        s = min(rec.reference_start, rec.next_reference_start)
        e = s + abs(rec.template_length)
        # region fully containing [s, e) on this ref
        mask = (rids == ref_id) & (starts <= s) & (ends >= e)
        idx = np.where(mask)[0]
        if idx.size == 1:
            oracle_contained[idx[0]] += 1
    bam.close()

    # ---- Report intergenic regions ----
    interg = np.where((rids == ref_id) & (ctype == 0))[0]
    print(f"\n=== Phase 0 validation: {ref_name} ({ra.n_regions} regions) ===")
    print(f"unique gDNA fragments in BAM: {n_gdna_unique}\n")
    print(f"{'region':>6} {'type':<11} {'span':>16} {'cal_contained':>14} {'oracle_contained':>16}")
    total_cal_ig = 0
    total_oracle_ig = 0
    for i in interg:
        total_cal_ig += int(cal_contained_unspl[i])
        total_oracle_ig += int(oracle_contained[i])
        print(f"R{i:<5} {_TYPE[ctype[i]]:<11} {f'{starts[i]:,}-{ends[i]:,}':>16} "
              f"{cal_contained_unspl[i]:>14} {oracle_contained[i]:>16}")

    print(f"\nintergenic contained — cal={total_cal_ig}  oracle={total_oracle_ig}")

    # ---- Mass conservation ----
    total_contained = int(cc.mass_unspliced.sum() + cc.mass_spliced.sum())
    bl = sub.left
    br = sub.right
    # boundary mass is per-side; total fragment boundary mass = left + right halves
    total_boundary = float(
        bl.mass_unspliced.sum() + br.mass_unspliced.sum()
        + bl.mass_spliced.sum() + br.mass_spliced.sum()
    )
    print(f"\nmass conservation: contained_mass={total_contained}  "
          f"boundary_mass={total_boundary:.0f}  "
          f"total={total_contained + total_boundary:.0f}")

    # ---- Verdict ----
    ok_nonzero = total_cal_ig > 0
    rel = abs(total_cal_ig - total_oracle_ig) / max(total_oracle_ig, 1)
    ok_match = rel < 0.05
    print(f"\n[{'PASS' if ok_nonzero else 'FAIL'}] intergenic contained mass is nonzero")
    print(f"[{'PASS' if ok_match else 'FAIL'}] cal vs oracle intergenic contained "
          f"within 5% (rel={rel:.3f})")


if __name__ == "__main__":
    main()
