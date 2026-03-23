#!/usr/bin/env python3
"""Diagnose fragment length bugs by tracing FL computation on benchmark data.

Traces the FL pipeline for oracle and minimap2 BAMs from the pristine benchmark.
Identifies:
  B1: nRNA gets FL=genomic_footprint (no gap SJ correction) → inflates get_unique_frag_length
  B2: Overlap fix incorrectly credits partial intron matches
  B3: Multiple gap SJs create ambiguous FLs that should be excluded from training

Usage:
    conda activate rigel
    python scripts/debug/diagnose_fl_bugs.py
"""
import sys
import logging
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, "src")

from rigel.index import TranscriptIndex

# --- Paths ---
BENCH = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3")
INDEX = BENCH / "rigel_index"
COND = BENCH / "gdna_none_ss_0.95_nrna_none"
ORACLE_BAM = COND / "align_oracle" / "reads_namesort.bam"
MM2_BAM = COND / "align_minimap2" / "reads_namesort.bam"

# --- Load index ---
logger.info("Loading index...")
idx = TranscriptIndex.load(str(INDEX))
ctx = idx.resolver

# Build some lookup tables
tx_df = idx.t_df
is_nrna = dict(zip(tx_df["t_index"], tx_df["is_synthetic_nrna"]))
t_id_map = dict(zip(tx_df["t_index"], tx_df["t_id"]))

logger.info(f"Index: {len(tx_df)} transcripts, "
            f"{sum(1 for v in is_nrna.values() if v)} nRNA synthetics")

# --- Scan BAM and trace FL computation ---
def trace_fl_from_bam(bam_path, label, max_frags=500000):
    """Scan BAM and collect detailed FL diagnostics."""
    logger.info(f"Scanning {label}: {bam_path}")

    # Use the native scanner to get resolved fragments
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer

    cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True)

    # We can't easily hook into the middle of the scan, so let's manually
    # trace by using the resolver on individual fragments
    # Instead, let's collect the FL observations that the scanner produces
    import pysam

    stats = {
        "total_pairs": 0,
        "unique_mappers": 0,
        "fl_collected": 0,
        "fl_ambiguous": 0,  # get_unique_frag_length returns -1
        "nrna_in_tinds": 0,  # Fragments where nRNA is in t_inds
        "nrna_only_tinds": 0,  # Fragments where ONLY nRNA is in t_inds
        "nrna_different_fl": 0,  # nRNA has different FL than mRNA
        "nrna_inflated_fl": 0,  # nRNA FL > 2x mRNA FL
        "spliced_annot_with_nrna": 0,
        "fl_vals": [],
        "fl_splice_types": [],
        "nrna_fl_vals": [],
        "mrna_fl_vals": [],
        "gap_count_distribution": Counter(),  # How many gaps per fragment
        "multi_gap_sj_count": 0,  # Gaps with multiple possible SJ corrections
    }

    # Manual resolution approach: process fragments through resolver
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    # Group reads by name for paired-end
    current_name = None
    current_reads = []

    def process_pair(reads):
        stats["total_pairs"] += 1
        if stats["total_pairs"] > max_frags:
            return

        # Check if unique mapper
        r1 = reads[0]
        nh = r1.get_tag("NH") if r1.has_tag("NH") else 1
        is_unique = (nh == 1)

        if not is_unique:
            return

        stats["unique_mappers"] += 1

        # Build exon/intron blocks from the reads
        exon_ref_ids = []
        exon_starts = []
        exon_ends = []
        exon_strands = []
        intron_ref_ids = []
        intron_starts = []
        intron_ends = []
        intron_strands = []

        for read in reads:
            if read.is_unmapped:
                continue
            ref_name = read.reference_name
            ref_id = ctx.ref_id(ref_name) if hasattr(ctx, 'ref_id') else -1
            if ref_id < 0:
                continue

            strand = 1 if read.is_reverse else 0  # STRAND_POS=0, STRAND_NEG=1
            pos = read.reference_start

            # Parse CIGAR into exon and intron blocks
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # M, =, X
                    exon_ref_ids.append(ref_id)
                    exon_starts.append(pos)
                    exon_ends.append(pos + length)
                    exon_strands.append(strand)
                    pos += length
                elif op == 3:  # N (intron)
                    intron_ref_ids.append(ref_id)
                    intron_starts.append(pos)
                    intron_ends.append(pos + length)
                    intron_strands.append(strand)
                    pos += length
                elif op == 2:  # D
                    pos += length
                elif op == 1:  # I
                    pass
                elif op == 4 or op == 5:  # S, H
                    pass
                else:
                    pos += length

        if not exon_ref_ids:
            return

        # Compute genomic footprint
        all_positions = list(zip(exon_starts, exon_ends))
        gfp = max(e for _, e in all_positions) - min(s for s, _ in all_positions)

        # Resolve through the native resolver
        result = ctx.resolve(
            exon_ref_ids, exon_starts, exon_ends, exon_strands,
            intron_ref_ids, intron_starts, intron_ends, intron_strands,
            gfp
        )

        if result is None:
            return

        # result is a ResolvedFragment or tuple
        # Access t_inds and frag_lengths
        try:
            t_inds = list(result.t_inds)
            frag_lengths = list(result.frag_lengths)
            splice_type = int(result.splice_type)
        except AttributeError:
            return

        if not t_inds:
            return

        # Classify: does this fragment have nRNA in t_inds?
        has_nrna = any(is_nrna.get(t, False) for t in t_inds)
        mrna_inds = [t for t in t_inds if not is_nrna.get(t, False)]
        nrna_inds = [t for t in t_inds if is_nrna.get(t, False)]

        if has_nrna:
            stats["nrna_in_tinds"] += 1
            if not mrna_inds:
                stats["nrna_only_tinds"] += 1

        # Get FL per transcript
        tid_to_fl = {}
        for i, t in enumerate(t_inds):
            if i < len(frag_lengths):
                tid_to_fl[t] = frag_lengths[i]

        mrna_fls = [tid_to_fl.get(t, -1) for t in mrna_inds if tid_to_fl.get(t, -1) > 0]
        nrna_fls = [tid_to_fl.get(t, -1) for t in nrna_inds if tid_to_fl.get(t, -1) > 0]

        if has_nrna and mrna_fls and nrna_fls:
            # Compare nRNA FL vs mRNA FL
            mrna_fl = mrna_fls[0]
            nrna_fl = nrna_fls[0]
            if nrna_fl != mrna_fl:
                stats["nrna_different_fl"] += 1
            if nrna_fl > 2 * mrna_fl:
                stats["nrna_inflated_fl"] += 1

        # get_unique_frag_length equivalent
        ufl = -1
        for fl in frag_lengths:
            if fl > 0:
                if ufl == -1:
                    ufl = fl
                elif fl != ufl:
                    ufl = -1
                    break

        if ufl > 0:
            stats["fl_collected"] += 1
            stats["fl_vals"].append(ufl)
            stats["fl_splice_types"].append(splice_type)

            if has_nrna and nrna_fls:
                stats["nrna_fl_vals"].append(nrna_fls[0])
            if mrna_fls:
                stats["mrna_fl_vals"].append(mrna_fls[0])

            if splice_type == 0 and has_nrna:  # SPLICE_SPLICED_ANNOT = 0
                stats["spliced_annot_with_nrna"] += 1
        else:
            stats["fl_ambiguous"] += 1

        if stats["total_pairs"] % 100000 == 0:
            logger.info(f"  {label}: {stats['total_pairs']} pairs processed, "
                        f"{stats['fl_collected']} FL collected")

    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        name = read.query_name
        if current_name is None:
            current_name = name
            current_reads = [read]
        elif name == current_name:
            current_reads.append(read)
        else:
            if len(current_reads) >= 2:
                process_pair(current_reads)
            current_name = name
            current_reads = [read]
            if stats["total_pairs"] > max_frags:
                break

    if len(current_reads) >= 2 and stats["total_pairs"] <= max_frags:
        process_pair(current_reads)

    bam.close()
    return stats

# --- Run traces ---
logger.info("=" * 70)
logger.info("ORACLE BAM TRACE")
logger.info("=" * 70)
orc_stats = trace_fl_from_bam(ORACLE_BAM, "oracle", max_frags=200000)

logger.info("\n" + "=" * 70)
logger.info("MINIMAP2 BAM TRACE")
logger.info("=" * 70)
mm2_stats = trace_fl_from_bam(MM2_BAM, "minimap2", max_frags=200000)

# --- Analysis ---
print("\n" + "=" * 80)
print("FRAGMENT LENGTH DIAGNOSTIC RESULTS")
print("=" * 80)

for label, stats in [("Oracle", orc_stats), ("Minimap2", mm2_stats)]:
    print(f"\n--- {label} ---")
    print(f"  Total pairs processed:         {stats['total_pairs']:>10,}")
    print(f"  Unique mappers:                {stats['unique_mappers']:>10,}")
    print(f"  FL collected (unambiguous):     {stats['fl_collected']:>10,}")
    print(f"  FL ambiguous (disagree):        {stats['fl_ambiguous']:>10,}")
    print(f"  Fragments with nRNA in t_inds: {stats['nrna_in_tinds']:>10,}")
    print(f"  Fragments with nRNA ONLY:      {stats['nrna_only_tinds']:>10,}")
    print(f"  nRNA has different FL than mRNA:{stats['nrna_different_fl']:>10,}")
    print(f"  nRNA FL > 2x mRNA FL:          {stats['nrna_inflated_fl']:>10,}")
    print(f"  SPLICED_ANNOT with nRNA:       {stats['spliced_annot_with_nrna']:>10,}")

    if stats["fl_vals"]:
        fl_arr = np.array(stats["fl_vals"])
        sp_arr = np.array(stats["fl_splice_types"])

        # Overall FL distribution
        print(f"\n  Overall FL: mean={fl_arr.mean():.1f}, median={np.median(fl_arr):.1f}, "
              f"mode≈{Counter(fl_arr).most_common(1)[0][0]}")

        # Per splice type
        for st, st_name in [(0, "SPLICED_ANNOT"), (1, "SPLICED_UNANNOT"), (2, "UNSPLICED")]:
            mask = sp_arr == st
            if mask.sum() > 0:
                sub = fl_arr[mask]
                mode_val = Counter(sub).most_common(1)[0][0]
                print(f"  {st_name}: n={mask.sum()}, mean={sub.mean():.1f}, "
                      f"median={np.median(sub):.1f}, mode={mode_val}")

        # FL > 1000 (likely corrupted)
        n_big = (fl_arr > 1000).sum()
        print(f"\n  FL > 1000 (potential corruption): {n_big} ({n_big/len(fl_arr)*100:.2f}%)")
        if n_big > 0:
            big_sp = sp_arr[fl_arr > 1000]
            for st, st_name in [(0, "SPLICED_ANNOT"), (1, "SPLICED_UNANNOT"), (2, "UNSPLICED")]:
                print(f"    {st_name}: {(big_sp == st).sum()}")

    if stats["nrna_fl_vals"]:
        nfl = np.array(stats["nrna_fl_vals"])
        print(f"\n  nRNA FL values: n={len(nfl)}, mean={nfl.mean():.1f}, "
              f"median={np.median(nfl):.1f}")
        print(f"  nRNA FL > 1000: {(nfl > 1000).sum()}")

    if stats["mrna_fl_vals"]:
        mfl = np.array(stats["mrna_fl_vals"])
        print(f"  mRNA FL values: n={len(mfl)}, mean={mfl.mean():.1f}, "
              f"median={np.median(mfl):.1f}")

# --- BUG-SPECIFIC DIAGNOSIS ---
print("\n" + "=" * 80)
print("BUG DIAGNOSIS")
print("=" * 80)

print("""
BUG 1: nRNA FL Contamination
  nRNA transcripts have a single giant "exon" covering the entire locus.
  In compute_frag_lengths():
    - The nRNA has NO introns in the SJ gap index
    - For paired-end fragments with a gap, the gap SJ correction finds
      ZERO overlapping SJs for the nRNA transcript
    - Result: nRNA FL = footprint (no correction) = genomic distance between mates
    - This is typically ~5000+ bp for fragments spanning introns
  
  In get_unique_frag_length():
    - If nRNA FL != mRNA FL, returns -1 (ambiguous) → NOT collected, so far OK
    - BUT: If nRNA is the ONLY transcript in t_inds (pruned by overhang gate),
      then FL = nRNA's inflated value → COLLECTED AS TRAINING DATA
    - This is the catastrophic path: inflated FL enters SPLICED_ANNOT training
""")

print(f"Oracle: nRNA-only fragments: {orc_stats['nrna_only_tinds']}")
print(f"MM2:    nRNA-only fragments: {mm2_stats['nrna_only_tinds']}")
print(f"Oracle: nRNA different FL:   {orc_stats['nrna_different_fl']}")
print(f"MM2:    nRNA different FL:   {mm2_stats['nrna_different_fl']}")

print(f"""
BUG 2: Overlap Fix Over-Credits Partial Introns
  The overlap fix changed:
    BEFORE: if (hs >= gs && he <= ge) t_gap_size[ti] += (he - hs)
    AFTER:  overlap = min(he,ge) - max(hs,gs); if (overlap > 0) t_gap_size[ti] += overlap
  
  Problem: When a gap only partially overlaps an intron, the overlap amount
  is credited. But the gap itself is NOT the intron — the gap is the
  unsequenced region between mates. If the intron is larger than the gap,
  crediting the overlap creates an UNDER-corrected FL (FL too large).
  
  For large introns that fully contain the gap, both methods give gap_size.
  For introns partially overlapping the gap, the overlap fix credit is wrong.
  The original strict containment was more conservative but could also miss
  legitimate corrections.

BUG 3: Multiple Gap SJs Create Ambiguous FLs
  When a gap overlaps multiple SJs for the same transcript, ALL their
  overlaps are summed in t_gap_size. This can over-correct or create
  incorrect FL values. Fragments with multiple possible gap SJ corrections
  should probably be excluded from FL training.
""")

print("=== DIAGNOSIS COMPLETE ===")
