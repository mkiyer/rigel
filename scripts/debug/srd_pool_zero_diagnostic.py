"""Debug: dump exon_bp/intron_bp per-candidate for unique unspliced fragments.

Stops after the FIRST chunk (don't scan whole 18 GB BAM). Targets a
gDNA-rich library where unique-mapper unspliced intronic/intergenic
reads should be abundant.
"""
from __future__ import annotations

import sys

import numpy as np

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType
from rigel.calibration._categorize import (
    FragmentCategory,
    N_CATEGORIES,
    build_t_to_ref_id,
    categorize_chunk,
)


def main():
    bam = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna80m/rigel/annotated.bam"
    )
    index_dir = (
        sys.argv[2]
        if len(sys.argv) > 2
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
    )
    print(f"BAM:   {bam}")
    print(f"INDEX: {index_dir}")

    index = TranscriptIndex.load(index_dir)
    scan_cfg = BamScanConfig(sj_strand_tag="auto")

    print("Scanning ...")
    _stats, _strand_models, _frag_models, buffer = scan_and_buffer(
        bam, index, scan_cfg,
    )

    t_to_strand = np.asarray(index.t_to_strand_arr, dtype=np.int8)
    t_to_ref_id = build_t_to_ref_id(index.t_df)
    ref_codes = index.t_df["ref"].cat.categories.tolist()

    chunk = next(buffer.iter_chunks())

    print(f"\nChunk #1: size={chunk.size}")
    cc = categorize_chunk(chunk, t_to_strand, t_to_ref_id)
    cnts = cc.n_per_category[:N_CATEGORIES]
    print("\nCategories on chunk #1 (no unique-map filter yet):")
    for c in FragmentCategory:
        print(f"  {c.name:32s} = {cnts[c]}")

    unique_mask = chunk.num_hits == 1
    print(f"\n  num_hits==1: {unique_mask.sum()} / {chunk.size}")
    cnts_unique = np.bincount(cc.category[unique_mask], minlength=N_CATEGORIES)
    print("\nCategories on chunk #1 AFTER unique-map filter:")
    for c in FragmentCategory:
        print(f"  {c.name:32s} = {cnts_unique[c]}")

    splice_type = chunk.splice_type
    spliced = (splice_type == SpliceType.SPLICED_ANNOT) | (splice_type == SpliceType.SPLICED_UNANNOT)
    t_offsets = chunk.t_offsets
    n_cands = np.diff(t_offsets).astype(np.intp)
    has_cands = n_cands > 0
    target = unique_mask & ~spliced & has_cands
    print(f"\nUNIQUE & UNSPLICED & has_cands: {target.sum()}")

    if not target.any():
        print("Nothing to inspect.")
        return

    safe_starts = np.where(has_cands, t_offsets[:-1], 0).astype(np.intp)
    min_intron = np.minimum.reduceat(chunk.intron_bp, safe_starts)
    max_exon = np.maximum.reduceat(chunk.exon_bp, safe_starts)

    t_min_intron = min_intron[target]
    t_max_exon = max_exon[target]
    t_n_cands = n_cands[target]

    print(f"\n  unique-unspliced n_cands percentiles 10/50/90/99: {np.percentile(t_n_cands, [10,50,90,99])}")
    print(f"  unique-unspliced min(intron_bp) percentiles 10/50/75/90/95/99/99.9/100: "
          f"{np.percentile(t_min_intron, [10,50,75,90,95,99,99.9,100])}")
    print(f"  unique-unspliced max(exon_bp)   percentiles 10/50/75/90/95/99: "
          f"{np.percentile(t_max_exon, [10,50,75,90,95,99])}")
    print(f"  -> n with min_intron > 5: {(t_min_intron > 5).sum()}")
    print(f"  -> n with min_intron > 50: {(t_min_intron > 50).sum()}")
    print(f"  -> n with max_exon == 0:  {(t_max_exon == 0).sum()}")

    # Per-candidate dissection.
    print("\n=== 10 unique-unspliced frags with n_cands >= 5 ===")
    high_cand_idx = np.where(target & (n_cands >= 5))[0]
    rng = np.random.default_rng(0)
    sample = rng.choice(high_cand_idx, size=min(10, len(high_cand_idx)), replace=False) if len(high_cand_idx) else []
    for i in sorted(sample):
        s = int(t_offsets[i]); e = int(t_offsets[i+1])
        ibp = chunk.intron_bp[s:e]
        ebp = chunk.exon_bp[s:e]
        ti = chunk.t_indices[s:e]
        try:
            ref_name = ref_codes[t_to_ref_id[ti[0]]]
        except Exception:
            ref_name = "?"
        rl = int(chunk.read_length[i])
        gs = int(chunk.genomic_start[i])
        gf = int(chunk.genomic_footprint[i])
        all_intronic_5 = int((ibp > 5).all())
        any_exonic_fit = int((ibp <= 5).any())
        print(f"\n  i={i} ref={ref_name} rl={rl} g0={gs} gf={gf} n_cands={e-s}")
        print(f"    intron_bp: {ibp.tolist()}")
        print(f"    exon_bp:   {ebp.tolist()}")
        print(f"    every cand intronic (>5)?={all_intronic_5}  any cand exon_fit (<=5)?={any_exonic_fit}")

    print("\n=== 10 unique-unspliced frags with n_cands == 1 ===")
    low_cand_idx = np.where(target & (n_cands == 1))[0]
    if len(low_cand_idx):
        sample2 = rng.choice(low_cand_idx, size=min(10, len(low_cand_idx)), replace=False)
        for i in sorted(sample2):
            s = int(t_offsets[i]); e = int(t_offsets[i+1])
            ibp = chunk.intron_bp[s:e]
            ebp = chunk.exon_bp[s:e]
            ti = chunk.t_indices[s:e]
            try:
                ref_name = ref_codes[t_to_ref_id[ti[0]]]
            except Exception:
                ref_name = "?"
            rl = int(chunk.read_length[i])
            gs = int(chunk.genomic_start[i])
            gf = int(chunk.genomic_footprint[i])
            print(f"  i={i} ref={ref_name} rl={rl} g0={gs} gf={gf} t={int(ti[0])} intron={int(ibp[0])} exon={int(ebp[0])}")


if __name__ == "__main__":
    main()
