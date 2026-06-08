#!/usr/bin/env python
"""Strand-stratified gDNA->RNA leak diagnostic (em_strand line of work).

Tests the falsifiable prediction of the rho=0 spurious-penalty theory
(docs/em_strand/01_theory_and_fix.md):

  At a stranded library (ss99), in-locus gDNA that LEAKS into the RNA pool should be
  *co-stranded-enriched* -- a gDNA fragment whose genomic orientation matches the
  overlapping transcript's strand gets the RNA per-fragment reward log(p_sense)>>log(0.5)
  and is stolen, while antisense gDNA gets log(p_antisense)<<log(0.5) and is caught.
  At an unstranded library (ss50) the per-fragment strand factor cancels, so the leak
  should be ~strand-symmetric.

It reads the benchmark's annotated BAM (truth in qname, ZF bit 0x04 = called gDNA,
ZL>=0 = in a locus) and the reference GTF (gene strands), and reports, for each
condition, the gDNA->RNA leak rate split by co-stranded vs antisense.

Usage:
  python scripts/debug/em_strand_leak_stratify.py [--suite DIR]
"""

from __future__ import annotations

import argparse
import bisect
from collections import defaultdict
from pathlib import Path

import pysam

from rigel.sim.truth import parse_origin

AF_GDNA_BIT = 0x04
DEFAULT_SUITE = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"
CONDITIONS = [
    "gdna_high_ss_0.99_nrna_none_capture_off",
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_off",
    "gdna_high_ss_0.50_nrna_none_capture_on",
]


def load_gene_intervals(gtf_path: Path):
    """Return {chrom: (sorted_starts, [(start,end,strand), ...])} from gene lines."""
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] not in ("gene", "transcript"):
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            by_chrom[chrom].append((start, end, strand))
    out = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        out[chrom] = ([iv[0] for iv in ivs], ivs)
    return out


def overlapping_strands(genes, chrom: str, start: int, end: int) -> set[str]:
    """Set of gene strands overlapping [start,end) on chrom."""
    if chrom not in genes:
        return set()
    starts, ivs = genes[chrom]
    # genes overlapping must have start < end_query; scan a window left of end.
    hi = bisect.bisect_right(starts, end)
    strands = set()
    # walk back while gene.end > start (genes can be long; small synthetic genome)
    for i in range(hi - 1, -1, -1):
        g_start, g_end, g_strand = ivs[i]
        if g_end <= start:
            # since sorted by start, earlier ones can still overlap (long genes);
            # but for this synthetic genome genes are short -- bound the walk.
            if g_start < start - 1_000_000:
                break
            continue
        if g_start < end and g_end > start:
            strands.add(g_strand)
    return strands


def analyze_condition(cond_dir: Path, genes, bam_name: str = "annotated.bam"):
    bam_path = cond_dir / bam_name
    if not bam_path.exists():
        return None

    # buckets: [co-stranded, antisense] x [leaked(called RNA), caught(called gDNA)]
    stats = {
        "co": [0, 0],   # [leaked, caught]
        "anti": [0, 0],
        "ambig": [0, 0],  # overlaps both strands
        "nostrand": [0, 0],  # in-locus but no gene strand resolved
    }
    seen: set[str] = set()
    n_gdna_inlocus = 0
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        if q in seen:
            continue
        seen.add(q)
        org = parse_origin(q)
        if org is None or org.kind != "gdna":
            continue
        try:
            zf = int(r.get_tag("ZF")) if r.has_tag("ZF") else 0
        except KeyError:
            continue
        # "in-locus" = genomically overlaps a transcript span (truth-based), NOT the
        # pipeline's ZL: at ss99 the pipeline kicks antisense gDNA out to intergenic
        # (ZL=-1), which would bias the denominator. We want every gDNA fragment that
        # physically sits on a gene, regardless of how the pipeline routed it.
        gd_strand = "+" if org.strand == "f" else "-"
        strands = overlapping_strands(genes, org.ref, org.start, org.end)
        if not strands:
            continue  # intergenic by genomic position; not the co-location case
        n_gdna_inlocus += 1
        called_gdna = bool(zf & AF_GDNA_BIT)
        leaked_idx = 0 if not called_gdna else 1  # 0=leaked(RNA), 1=caught(gDNA)

        if len(strands) > 1:
            stats["ambig"][leaked_idx] += 1
        elif gd_strand in strands:
            stats["co"][leaked_idx] += 1
        else:
            stats["anti"][leaked_idx] += 1
    bam.close()
    return stats, n_gdna_inlocus


def fmt_bucket(name, b):
    leaked, caught = b
    tot = leaked + caught
    rate = 100.0 * leaked / tot if tot else 0.0
    return f"  {name:>9}: n={tot:7d}  leak->RNA={leaked:7d} ({rate:5.1f}%)  caught={caught:7d}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=DEFAULT_SUITE)
    ap.add_argument("--bam-name", default="annotated.bam",
                    help="BAM filename within each condition dir (e.g. annotated.baseline_cur.bam)")
    args = ap.parse_args()
    suite = Path(args.suite)
    gtf = suite / "reference" / "genes.gtf"
    if not gtf.exists():
        # try per-condition reference
        gtf = suite / CONDITIONS[0] / "reference" / "genes.gtf"
    genes = load_gene_intervals(gtf)
    n_genes = sum(len(v[1]) for v in genes.values())
    print(f"Loaded {n_genes} genes from {gtf}\n")

    for cond in CONDITIONS:
        cdir = suite / cond
        if not cdir.exists():
            print(f"[skip] {cond} (missing)")
            continue
        res = analyze_condition(cdir, genes, args.bam_name)
        if res is None:
            print(f"[skip] {cond} (no annotated bam)")
            continue
        stats, n_inlocus = res
        print(f"=== {cond} ===")
        print(f"  in-locus gDNA fragments: {n_inlocus}")
        for k in ("co", "anti", "ambig", "nostrand"):
            print(fmt_bucket(k, stats[k]))
        co_rate = 100.0 * stats["co"][0] / max(sum(stats["co"]), 1)
        anti_rate = 100.0 * stats["anti"][0] / max(sum(stats["anti"]), 1)
        print(f"  --> co-stranded leak {co_rate:.1f}%  vs  antisense leak {anti_rate:.1f}%"
              f"   (asymmetry = {co_rate - anti_rate:+.1f} pts)\n")


if __name__ == "__main__":
    main()
