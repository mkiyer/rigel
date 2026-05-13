#!/usr/bin/env python3
"""Drop duplicate-exon transcripts from a GTF.

Two transcripts are considered duplicates if they share (chrom, strand) and
have identical sorted exon coordinates. The first transcript encountered per
group is kept; all others are removed (along with their exon/CDS/UTR lines).
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict


def _open(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def parse_attrs(attrs: str) -> dict:
    out = {}
    for field in attrs.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " in field:
            k, v = field.split(" ", 1)
            out[k] = v.strip().strip('"')
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Pass 1: gather exon coords per transcript_id
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    tx_meta: dict[str, tuple[str, str]] = {}  # tx_id -> (chrom, strand)
    with _open(args.inp, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = parse_attrs(f[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            exons[tid].append((int(f[3]), int(f[4])))
            tx_meta[tid] = (f[0], f[6])

    # Group identical transcripts; keep first per group
    groups: dict[tuple, list[str]] = defaultdict(list)
    for tid, ex in exons.items():
        chrom, strand = tx_meta[tid]
        key = (chrom, strand, tuple(sorted(ex)))
        groups[key].append(tid)

    drop: set[str] = set()
    n_groups_dup = 0
    for tids in groups.values():
        if len(tids) > 1:
            n_groups_dup += 1
            for t in tids[1:]:
                drop.add(t)

    print(f"duplicate groups: {n_groups_dup}", file=sys.stderr)
    print(f"transcripts dropped: {len(drop)}", file=sys.stderr)
    print(f"total transcripts: {len(exons)}", file=sys.stderr)

    # Pass 2: rewrite GTF
    kept_lines = 0
    dropped_lines = 0
    with _open(args.inp, "r") as fin, _open(args.out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.split("\t")
            if len(f) < 9:
                fout.write(line)
                continue
            attrs = parse_attrs(f[8])
            tid = attrs.get("transcript_id")
            if tid and tid in drop:
                dropped_lines += 1
                continue
            fout.write(line)
            kept_lines += 1
    print(f"lines kept: {kept_lines}, dropped: {dropped_lines}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
