#!/usr/bin/env python
"""Decode the FATE of gDNA fragments (em_strand investigation).

For each gDNA-truth fragment, classify:
  - location: EXONIC (overlaps an exon), INTRONIC (overlaps a transcript span but no
    exon), or INTERGENIC (no span overlap)
  - strand: co-stranded / antisense (genomic orientation vs overlapping transcript strand)
  - fate (ZF bitfield): intergenic-fastpath (0x25), em-gDNA (0x05), mRNA (0x03),
    nRNA (0x09/0x19), chimeric (0x40), unresolved/other.

Answers: are exon-overlapping co-stranded gDNA fragments entering the EM? Where do
antisense exon gDNA go, and is it strand-dependent (ss99 vs ss50)?

Usage: python scripts/debug/em_strand_fate.py [--bam-name annotated.baseline_cur.bam]
"""

from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from rigel.sim.truth import parse_origin

DEFAULT_SUITE = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"
CONDITIONS = [
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_on",
]


def zf_fate(zf: int) -> str:
    if zf == 0x25:
        return "intergenic(0x25)"
    if zf == 0x05:
        return "em-gDNA(0x05)"
    if zf == 0x03:
        return "mRNA(0x03)"
    if zf in (0x09, 0x19):
        return "nRNA"
    if zf == 0x40:
        return "chimeric"
    if zf == 0x00:
        return "unresolved(0x00)"
    return f"other(0x{zf:02x})"


def load_intervals(gtf_path: Path, feature: str):
    by_chrom = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != feature:
                continue
            by_chrom[f[0]].append((int(f[3]) - 1, int(f[4]), f[6]))
    out = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        out[chrom] = ([iv[0] for iv in ivs], ivs)
    return out


def overlaps(iv_map, chrom, start, end):
    """Return set of strands of intervals overlapping [start,end)."""
    if chrom not in iv_map:
        return set()
    starts, ivs = iv_map[chrom]
    hi = bisect.bisect_right(starts, end)
    out = set()
    for i in range(hi - 1, -1, -1):
        s, e, strand = ivs[i]
        if e <= start:
            if s < start - 2_000_000:
                break
            continue
        if s < end and e > start:
            out.add(strand)
    return out


def analyze(cond_dir: Path, spans, exons, bam_name: str):
    bam_path = cond_dir / bam_name
    if not bam_path.exists():
        bam_path = cond_dir / "annotated.bam"
    # table[(location, strand)] -> Counter(fate)
    table = defaultdict(Counter)
    seen = set()
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        if q in seen:
            continue
        seen.add(q)
        org = parse_origin(q)
        if org is None or org.kind != "gdna" or not r.has_tag("ZF"):
            continue
        zf = int(r.get_tag("ZF"))
        span_str = overlaps(spans, org.ref, org.start, org.end)
        if not span_str:
            location = "INTERGENIC"
        else:
            ex_str = overlaps(exons, org.ref, org.start, org.end)
            location = "EXONIC" if ex_str else "INTRONIC"
        # strand vs overlapping transcript (use span strands; skip ambiguous)
        ref_str = span_str if span_str else set()
        if len(ref_str) != 1:
            strand = "ambig/none"
        else:
            gd = "+" if org.strand == "f" else "-"
            strand = "co" if gd in ref_str else "anti"
        table[(location, strand)][zf_fate(zf)] += 1
    bam.close()
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=DEFAULT_SUITE)
    ap.add_argument("--bam-name", default="annotated.baseline_cur.bam")
    args = ap.parse_args()
    suite = Path(args.suite)
    gtf = suite / "reference" / "genes.gtf"
    spans = load_intervals(gtf, "transcript")
    exons = load_intervals(gtf, "exon")

    for cond in CONDITIONS:
        cdir = suite / cond
        if not cdir.exists():
            continue
        table = analyze(cdir, spans, exons, args.bam_name)
        print(f"\n{'='*78}\n=== {cond} ===\n{'='*78}")
        for loc in ("EXONIC", "INTRONIC", "INTERGENIC"):
            for strand in ("co", "anti", "ambig/none"):
                c = table.get((loc, strand))
                if not c:
                    continue
                tot = sum(c.values())
                parts = ", ".join(f"{k} {v} ({100*v/tot:.0f}%)"
                                  for k, v in c.most_common())
                print(f"  {loc:>10} {strand:>10}  n={tot:7d}  | {parts}")


if __name__ == "__main__":
    main()
