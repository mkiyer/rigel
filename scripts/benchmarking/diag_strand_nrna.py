#!/usr/bin/env python3
"""Diagnose strand behavior for nRNA-leaked LINC02912 reads."""
import pysam

BAM = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/align_oracle/annotated_default.sorted.bam"
)

nrna_reads = []
mrna_reads = []

with pysam.AlignmentFile(BAM, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_read2:
            continue
        tx = read.get_tag("ZT")
        pool = read.get_tag("ZP")
        gene = read.get_tag("ZG")

        if tx == "ENST00000624314.1" or (
            pool == "nRNA" and "ENSG00000249859" in gene
        ):
            flag = read.flag
            is_reverse = read.is_reverse
            is_read1 = read.is_read1
            pos = read.reference_start
            cigar = read.cigarstring
            posterior = read.get_tag("ZW")
            n_cand = read.get_tag("ZN")
            fclass = read.get_tag("ZC")
            splice = read.get_tag("ZS")

            entry = {
                "qname": read.query_name,
                "flag": flag,
                "is_reverse": is_reverse,
                "is_read1": is_read1,
                "pos": pos,
                "pool": pool,
                "tx": tx,
                "gene": gene,
                "posterior": posterior,
                "n_cand": n_cand,
                "class": fclass,
                "splice": splice,
                "cigar": cigar,
            }

            if pool == "nRNA":
                if len(nrna_reads) < 10:
                    nrna_reads.append(entry)
            else:
                if len(mrna_reads) < 5:
                    mrna_reads.append(entry)

        if len(nrna_reads) >= 10 and len(mrna_reads) >= 5:
            break

print("=== nRNA-assigned reads (R1 only) ===")
for r in nrna_reads:
    print(
        f"  {r['qname'][:50]:50s}  flag={r['flag']:4d}  "
        f"is_reverse={r['is_reverse']}  is_read1={r['is_read1']}  "
        f"pos={r['pos']:8d}  tx={r['tx']:25s}  gene={r['gene']:25s}  "
        f"pool={r['pool']}  post={r['posterior']:.4f}  n_cand={r['n_cand']}  "
        f"splice={r['splice']}  class={r['class']}"
    )

print()
print("=== mRNA-assigned ENST00000624314.1 reads (R1 only) ===")
for r in mrna_reads:
    print(
        f"  {r['qname'][:50]:50s}  flag={r['flag']:4d}  "
        f"is_reverse={r['is_reverse']}  is_read1={r['is_read1']}  "
        f"pos={r['pos']:8d}  tx={r['tx']:25s}  gene={r['gene']:25s}  "
        f"pool={r['pool']}  post={r['posterior']:.4f}  n_cand={r['n_cand']}  "
        f"splice={r['splice']}  class={r['class']}"
    )

# Now check: what strand is the exon_strand for these reads?
# exon_strand in the pipeline is derived from the R1 strand of the mapped read.
# For RF library, R1 maps antisense. LINC02912 is on - strand.
# So R1 of LINC02912 should map to + strand (is_reverse=False for R1).
print()
print("=== Strand analysis ===")
print("LINC02912 gene strand: - (negative)")
print("PVT1 gene strand:      + (positive)")
print()
for r in nrna_reads[:3]:
    r1_maps_to = "+" if not r["is_reverse"] else "-"
    print(
        f"  Read {r['qname'][:40]}  R1 maps to {r1_maps_to} strand  "
        f"(is_reverse={r['is_reverse']})"
    )
    print(f"    → exon_strand would be: {'POS(1)' if not r['is_reverse'] else 'NEG(2)'}")
    print(f"    → vs PVT1 gene strand (+): {'SAME' if r1_maps_to == '+' else 'DIFF'}")
    print(f"    → vs LINC02912 gene strand (-): {'SAME' if r1_maps_to == '-' else 'DIFF'}")
