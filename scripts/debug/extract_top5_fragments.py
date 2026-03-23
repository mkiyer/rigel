#!/usr/bin/env python3
"""Extract and analyze fragments for top error-gap transcripts from annotated BAM.

For each of the top 5 transcripts, extract ALL reads whose true source
(from read name) is that transcript, and tabulate:
- Where minimap2 aligned them (chromosome, position)
- What rigel annotated them as (ZT=transcript, ZP=pool, ZW=posterior, ZS=splice)
- The number of candidates (ZN), locus (ZL)
- The alignment mapping quality, NM, and CIGAR

This gives us a per-fragment view of WHY rigel-minimap2 is failing.
"""

import pysam
import sys
from collections import Counter, defaultdict
from pathlib import Path

BAM = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"

# Top 5 transcripts by error gap
TOP5 = [
    "ENST00000229239.10",  # GAPDH
    "ENST00000573283.7",   # ACTG1
    "ENST00000631211.1",   # ENSG00000280800
    "ENST00000646664.1",   # ACTB
    "ENST00000436459.2",   # EEF1A1P5
]

OUTDIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/fragment_analysis")
OUTDIR.mkdir(exist_ok=True)


def get_tag(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def extract_transcript_fragments(bam_path, target_tx_prefix):
    """Extract all reads whose name starts with target_tx_prefix."""
    fragments = defaultdict(dict)  # qname -> {R1: read, R2: read, secondaries: [...]}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qname = read.query_name
            # Read names: ENST00000229239.10:1-229:f:57
            if not qname.startswith(target_tx_prefix):
                continue

            key = qname
            if key not in fragments:
                fragments[key] = {"reads": [], "seen": set()}

            # Deduplicate: same read can appear multiple times (primary + supplementary)
            read_id = (read.reference_name, read.reference_start, read.is_read1,
                       read.is_secondary, read.is_supplementary)
            if read_id not in fragments[key]["seen"]:
                fragments[key]["reads"].append(read)
                fragments[key]["seen"].add(read_id)

    return fragments


def analyze_fragments(fragments, target_tx, outpath):
    """Analyze extracted fragments and write detailed report."""
    # Parse read names to get truth fragment info
    # Format: TXID:frag_start-frag_end:strand:frag_id
    n_frags = len(fragments)
    print(f"\n  Total fragment groups: {n_frags}")

    # Classify by rigel annotation
    assigned_to = Counter()  # ZT value
    pool_counts = Counter()  # ZP value
    class_counts = Counter()  # ZC value
    splice_counts = Counter()  # ZS value
    locus_counts = Counter()  # ZL value
    ncand_hist = Counter()  # ZN value
    mapq_hist = Counter()
    chrom_counts = Counter()

    # Per-fragment detail
    details = []

    for qname, frag in fragments.items():
        reads = frag["reads"]
        # Get primary reads (not secondary/supplementary)
        primaries = [r for r in reads if not r.is_secondary and not r.is_supplementary]
        secondaries = [r for r in reads if r.is_secondary]

        if not primaries:
            continue

        # Use first primary for Z-tags (they're set per-fragment, same for both mates)
        p = primaries[0]
        zt = get_tag(p, "ZT", ".")
        zg = get_tag(p, "ZG", ".")
        zp = get_tag(p, "ZP", "intergenic")
        zw = get_tag(p, "ZW", 0.0)
        zc = get_tag(p, "ZC", "unknown")
        zh = get_tag(p, "ZH", 0)
        zn = get_tag(p, "ZN", 0)
        zs = get_tag(p, "ZS", "unknown")
        zl = get_tag(p, "ZL", -1)
        mapq = p.mapping_quality
        chrom = p.reference_name or "unmapped"
        nm = get_tag(p, "NM", 0)
        cigar = p.cigarstring or "unmapped"

        # Count secondary alignments (multimapper loci)
        n_secondary = len(secondaries)
        secondary_chroms = [r.reference_name for r in secondaries]

        assigned_to[zt] += 1
        pool_counts[zp] += 1
        class_counts[zc] += 1
        splice_counts[zs] += 1
        locus_counts[zl] += 1
        ncand_hist[zn] += 1
        mapq_hist[mapq] += 1
        chrom_counts[chrom] += 1

        details.append({
            "qname": qname,
            "chrom": chrom,
            "pos": p.reference_start,
            "mapq": mapq,
            "cigar": cigar,
            "nm": nm,
            "zt": zt,
            "zg": zg,
            "zp": zp,
            "zw": zw,
            "zc": zc,
            "zn": zn,
            "zs": zs,
            "zl": zl,
            "n_secondary": n_secondary,
            "sec_chroms": ",".join(secondary_chroms[:5]),
        })

    print(f"  Fragments with primaries: {len(details)}")

    # Summary report
    print(f"\n  === POOL ASSIGNMENT ===")
    for pool, cnt in pool_counts.most_common():
        pct = 100 * cnt / len(details)
        print(f"    {pool:>20s}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === FRAGMENT CLASS ===")
    for cls, cnt in class_counts.most_common():
        pct = 100 * cnt / len(details)
        print(f"    {cls:>20s}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === SPLICE TYPE ===")
    for sp, cnt in splice_counts.most_common():
        pct = 100 * cnt / len(details)
        print(f"    {sp:>20s}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === TOP ASSIGNED TRANSCRIPTS (ZT) ===")
    for zt, cnt in assigned_to.most_common(20):
        pct = 100 * cnt / len(details)
        print(f"    {zt:>35s}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === N_CANDIDATES (ZN) HISTOGRAM ===")
    for zn, cnt in sorted(ncand_hist.items()):
        pct = 100 * cnt / len(details)
        print(f"    ZN={zn:>4d}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === CHROMOSOME DISTRIBUTION ===")
    for chrom, cnt in chrom_counts.most_common(15):
        pct = 100 * cnt / len(details)
        print(f"    {chrom:>10s}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === MAPPING QUALITY ===")
    for mq, cnt in sorted(mapq_hist.items()):
        pct = 100 * cnt / len(details)
        print(f"    MQ={mq:>3d}: {cnt:>6d} ({pct:>5.1f}%)")

    print(f"\n  === N_SECONDARY ALIGNMENTS ===")
    sec_hist = Counter(d["n_secondary"] for d in details)
    for ns, cnt in sorted(sec_hist.items())[:15]:
        pct = 100 * cnt / len(details)
        print(f"    sec={ns:>3d}: {cnt:>6d} ({pct:>5.1f}%)")

    # Look at fragments NOT assigned to the correct transcript
    correct_tx = target_tx.split(".")[0]  # strip version for partial match
    misassigned = [d for d in details if not d["zt"].startswith(correct_tx)]
    n_correct = len(details) - len(misassigned)
    print(f"\n  === CORRECT vs MISASSIGNED ===")
    print(f"    Correctly assigned to {target_tx}: {n_correct} ({100*n_correct/len(details):.1f}%)")
    print(f"    Misassigned: {len(misassigned)} ({100*len(misassigned)/len(details):.1f}%)")

    if misassigned:
        # Where did misassigned fragments go?
        mis_zt = Counter(d["zt"] for d in misassigned)
        print(f"\n    Top misassignment targets:")
        for zt, cnt in mis_zt.most_common(20):
            pct = 100 * cnt / len(misassigned)
            print(f"      {zt:>35s}: {cnt:>6d} ({pct:>5.1f}%)")

        # Splice type of misassigned fragments
        mis_sp = Counter(d["zs"] for d in misassigned)
        print(f"\n    Splice type of misassigned fragments:")
        for sp, cnt in mis_sp.most_common():
            pct = 100 * cnt / len(misassigned)
            print(f"      {sp:>20s}: {cnt:>6d} ({pct:>5.1f}%)")

        # Pool of misassigned fragments
        mis_pool = Counter(d["zp"] for d in misassigned)
        print(f"\n    Pool of misassigned fragments:")
        for pool, cnt in mis_pool.most_common():
            pct = 100 * cnt / len(misassigned)
            print(f"      {pool:>20s}: {cnt:>6d} ({pct:>5.1f}%)")

        # Chromosome distribution of misassigned
        mis_chrom = Counter(d["chrom"] for d in misassigned)
        print(f"\n    Chromosome distribution of misassigned:")
        for chrom, cnt in mis_chrom.most_common(15):
            pct = 100 * cnt / len(misassigned)
            print(f"      {chrom:>10s}: {cnt:>6d} ({pct:>5.1f}%)")

    # Write detailed per-fragment TSV
    import csv
    tsv_path = outpath / f"{target_tx}_fragments.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=details[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(details)
    print(f"\n  Wrote {len(details)} fragments to {tsv_path}")


def main():
    print("Loading BAM (name-sorted, streaming through entire file)...")
    print("This will take a few minutes for a 710MB BAM...")

    # Collect fragments for ALL top 5 transcripts in a single pass
    target_prefixes = {tx: tx.split(".")[0] for tx in TOP5}
    all_fragments = {tx: defaultdict(lambda: {"reads": [], "seen": set()}) for tx in TOP5}

    count = 0
    with pysam.AlignmentFile(BAM, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            count += 1
            if count % 5_000_000 == 0:
                print(f"  Processed {count:,} reads...")

            qname = read.query_name
            for tx, prefix in target_prefixes.items():
                if qname.startswith(prefix):
                    key = qname
                    frag = all_fragments[tx][key]
                    read_id = (read.reference_name, read.reference_start,
                               read.is_read1, read.is_secondary, read.is_supplementary)
                    if read_id not in frag["seen"]:
                        frag["reads"].append(read)
                        frag["seen"].add(read_id)
                    break

    print(f"\n  Total reads processed: {count:,}")

    for tx in TOP5:
        gene = {"ENST00000229239.10": "GAPDH",
                "ENST00000573283.7": "ACTG1",
                "ENST00000631211.1": "ENSG00000280800",
                "ENST00000646664.1": "ACTB",
                "ENST00000436459.2": "EEF1A1P5"}[tx]
        print(f"\n\n{'='*100}")
        print(f"  TRANSCRIPT: {tx} ({gene})")
        print(f"{'='*100}")
        analyze_fragments(all_fragments[tx], tx, OUTDIR)


if __name__ == "__main__":
    main()
