#!/usr/bin/env python3
"""Audit oracle BAM SEQ orientation against the reference FASTA."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from rigel.sim.genome import reverse_complement

MATCH_OPS = {0, 7, 8}
REF_SKIP_OPS = {2, 3}
QUERY_SKIP_OPS = {1, 4, 5}


def reference_sequence_for_record(record: pysam.AlignedSegment, fasta: pysam.FastaFile) -> str:
    """Return reference bases consumed by query-aligned CIGAR operations."""
    ref_pos = record.reference_start
    parts: list[str] = []
    for op, length in record.cigartuples or []:
        if op in MATCH_OPS:
            parts.append(fasta.fetch(record.reference_name, ref_pos, ref_pos + length).upper())
            ref_pos += length
        elif op in REF_SKIP_OPS:
            ref_pos += length
        elif op in QUERY_SKIP_OPS:
            continue
        else:
            raise ValueError(f"unsupported CIGAR op {op} in {record.query_name}")
    return "".join(parts)


def category(record: pysam.AlignedSegment) -> str:
    if record.query_name.startswith("gdna:"):
        return "gdna"
    if record.query_name.startswith("nrna_"):
        return "nrna"
    return "mrna"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--max-records", type=int, default=10000, help="0 means all records")
    parser.add_argument("--examples", type=int, default=8)
    args = parser.parse_args()

    counts: Counter[tuple[str, str, str]] = Counter()
    examples: dict[str, list[str]] = defaultdict(list)

    with pysam.FastaFile(str(args.fasta)) as fasta, pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for record_index, record in enumerate(bam):
            if args.max_records > 0 and record_index >= args.max_records:
                break
            ref_seq = reference_sequence_for_record(record, fasta)
            query_seq = (record.query_sequence or "").upper()
            if query_seq == ref_seq:
                status = "matches_reference"
            elif reverse_complement(query_seq) == ref_seq:
                status = "matches_reference_after_rc"
            else:
                status = "mismatch_other"
            read = "read1" if record.is_read1 else "read2"
            strand = "reverse" if record.is_reverse else "forward"
            counts[(category(record), read, strand, status)] += 1
            if status != "matches_reference" and len(examples[status]) < args.examples:
                examples[status].append(
                    f"{record.query_name}\tflag={record.flag}\t{record.reference_name}:"
                    f"{record.reference_start + 1}\t{record.cigarstring}\t{status}"
                )

    print("category\tread\tstrand\tstatus\tcount")
    for key, count in sorted(counts.items()):
        print("\t".join((*key, str(count))))
    for status, rows in examples.items():
        print(f"\nexamples: {status}")
        for row in rows:
            print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
