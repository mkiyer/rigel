"""Sample flowcell IDs from a BAM via pysam (no shell pipes)."""

from __future__ import annotations

import argparse
import collections

import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--n", type=int, default=50000)
    args = parser.parse_args()

    counts: collections.Counter[str] = collections.Counter()
    bam = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    for i, rec in enumerate(bam):
        if i >= args.n:
            break
        parts = rec.query_name.split(":")
        fc = parts[2] if len(parts) >= 3 else rec.query_name
        counts[fc] += 1
    for fc, n in counts.most_common():
        print(f"{n}\t{fc}")


if __name__ == "__main__":
    main()
