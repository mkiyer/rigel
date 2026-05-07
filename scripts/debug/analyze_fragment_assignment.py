#!/usr/bin/env python3
"""Correct fragment-level assignment analysis using ZF bitfield."""

import pysam
import sys
from collections import Counter

# ZF bitfield constants from rigel annotate.py
IS_RESOLVED = 0x01
IS_MRNA = 0x02
IS_GDNA = 0x04
IS_NRNA = 0x08
IS_SYNTHETIC = 0x10
IS_INTERGENIC = 0x20
IS_CHIMERIC = 0x40
IS_MULTIMAP_DROP = 0x80


def analyze_condition(bam_path, label):
    """Analyze fragment assignment accuracy for one condition."""
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Confusion: (true_pool, assigned_pool) -> count
    confusion = Counter()
    # ZF breakdown for true-gDNA reads
    zf_for_gdna = Counter()
    # Track gDNA reads that go to mRNA: which transcripts absorb them?
    gdna_to_tx = Counter()
    # Track read category (ZC) for gDNA reads assigned to mRNA
    gdna_as_mrna_zc = Counter()

    n = 0
    for read in bam:
        if read.is_read2 or read.is_secondary or read.is_supplementary:
            continue
        n += 1

        qname = read.query_name
        true_is_gdna = qname.startswith("gdna")

        try:
            zf = read.get_tag("ZF")
        except KeyError:
            zf = 0

        # Determine assigned pool from ZF
        if zf & IS_GDNA:
            assigned = "gDNA"
        elif zf & IS_NRNA:
            assigned = "nRNA"
        elif zf & IS_MRNA:
            assigned = "mRNA"
        elif zf & IS_CHIMERIC:
            assigned = "chimeric"
        elif zf & IS_MULTIMAP_DROP:
            assigned = "mm_drop"
        else:
            assigned = "unresolved"

        true_pool = "gDNA" if true_is_gdna else "mRNA"
        confusion[(true_pool, assigned)] += 1

        if true_is_gdna:
            zf_for_gdna[zf] += 1
            if assigned == "mRNA":
                try:
                    tx = read.get_tag("ZT")
                except KeyError:
                    tx = "unknown"
                gdna_to_tx[tx] += 1
                try:
                    zc = read.get_tag("ZC")
                except KeyError:
                    zc = "unknown"
                gdna_as_mrna_zc[zc] += 1

    bam.close()

    print(f"\n{'='*80}")
    print(f"  {label}")
    print(f"  BAM: {bam_path}")
    print(f"{'='*80}")
    print(f"\n  Total R1 fragments: {n:,}")

    # Confusion matrix
    print(f"\n  CONFUSION MATRIX:")
    print(f"  {'True':<8} {'Assigned':<12} {'Count':>12} {'%':>8}")
    print(f"  {'-'*45}")
    for k in sorted(confusion.keys()):
        pct = confusion[k] / n * 100
        print(f"  {k[0]:<8} {k[1]:<12} {confusion[k]:>12,}  {pct:>7.2f}%")

    # gDNA specifics
    total_gdna = sum(v for (k, _), v in confusion.items() if k == "gDNA")
    if total_gdna > 0:
        gdna_as_mrna = confusion.get(("gDNA", "mRNA"), 0)
        gdna_as_gdna = confusion.get(("gDNA", "gDNA"), 0)
        gdna_as_nrna = confusion.get(("gDNA", "nRNA"), 0)

        print(f"\n  gDNA FRAGMENT FATE ({total_gdna:,} total):")
        print(f"    Correctly → gDNA:    {gdna_as_gdna:>10,} ({gdna_as_gdna/total_gdna*100:.2f}%)")
        print(f"    Leaked → mRNA:       {gdna_as_mrna:>10,} ({gdna_as_mrna/total_gdna*100:.2f}%)")
        print(f"    Leaked → nRNA:       {gdna_as_nrna:>10,} ({gdna_as_nrna/total_gdna*100:.2f}%)")

        # ZF breakdown
        print(f"\n  ZF BITFIELD BREAKDOWN for true-gDNA reads:")
        for zf, cnt in zf_for_gdna.most_common(10):
            bits = []
            if zf & IS_RESOLVED:
                bits.append("resolved")
            if zf & IS_MRNA:
                bits.append("mRNA")
            if zf & IS_GDNA:
                bits.append("gDNA")
            if zf & IS_NRNA:
                bits.append("nRNA")
            if zf & IS_INTERGENIC:
                bits.append("intergenic")
            if zf & IS_CHIMERIC:
                bits.append("chimeric")
            print(f"    ZF=0x{zf:02x} ({' | '.join(bits):<35}) : {cnt:>10,} ({cnt/total_gdna*100:.1f}%)")

        # Pre-EM category (ZC) of gDNA reads that ended up as mRNA
        if gdna_as_mrna > 0:
            print(f"\n  PRE-EM CATEGORY (ZC) of gDNA reads misassigned to mRNA:")
            for zc, cnt in gdna_as_mrna_zc.most_common(10):
                print(f"    {zc:<25} : {cnt:>10,} ({cnt/gdna_as_mrna*100:.1f}%)")

        # Top transcript absorbers of gDNA
        if gdna_to_tx:
            print(f"\n  TOP TRANSCRIPTS absorbing gDNA fragments:")
            for tx, cnt in gdna_to_tx.most_common(15):
                print(f"    {tx:<25} : {cnt:>8,} ({cnt/gdna_as_mrna*100:.2f}%)")

    # mRNA specifics
    total_mrna = sum(v for (k, _), v in confusion.items() if k == "mRNA")
    if total_mrna > 0:
        mrna_as_gdna = confusion.get(("mRNA", "gDNA"), 0)
        mrna_as_mrna = confusion.get(("mRNA", "mRNA"), 0)
        mrna_as_nrna = confusion.get(("mRNA", "nRNA"), 0)
        print(f"\n  mRNA FRAGMENT FATE ({total_mrna:,} total):")
        print(f"    Correctly → mRNA:    {mrna_as_mrna:>10,} ({mrna_as_mrna/total_mrna*100:.2f}%)")
        print(f"    Misclass → gDNA:     {mrna_as_gdna:>10,} ({mrna_as_gdna/total_mrna*100:.2f}%)")
        print(f"    Misclass → nRNA:     {mrna_as_nrna:>10,} ({mrna_as_nrna/total_mrna*100:.2f}%)")


if __name__ == "__main__":
    base = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic"
    conditions = [
        ("gdna_none_ss_0.99_nrna_none", "CLEAN (no gDNA, SS=0.99)"),
        ("gdna_low_ss_0.99_nrna_none", "LOW gDNA (10%, SS=0.99)"),
        ("gdna_med_ss_0.99_nrna_none", "MED gDNA (50%, SS=0.99)"),
        ("gdna_high_ss_0.99_nrna_none", "HIGH gDNA (200%, SS=0.99)"),
        ("gdna_high_ss_0.50_nrna_none", "HIGH gDNA (200%, SS=0.50)"),
    ]

    for cond_name, label in conditions:
        bam_path = f"{base}/{cond_name}/annotated.bam"
        analyze_condition(bam_path, label)
