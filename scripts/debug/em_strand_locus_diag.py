#!/usr/bin/env python
"""Per-locus strand/budget diagnostic for the em_strand M-step design.

Answers the two questions that decide the M-step strand-balance design
(docs/em_strand/01_theory_and_fix.md §5-6):

  Q1 (anchor location): where does ANTISENSE gDNA live -- inside the locus (ZL>=0) or
     routed to intergenic (ZL=-1)? A locus-level strand-balance term needs the antisense
     gDNA *in the locus* to detect the 50/50 balance. If it's intergenic, the locus sees
     an all-sense pile and the term cannot work as designed.

  Q2 (budget vs assignment): at the loci where co-stranded gDNA leaks, is the calibration
     gDNA budget (loci.tsv gdna_prior_count) sized to the gDNA that is actually IN the
     locus, and how much does the EM assign (loci.tsv gdna)?

Inputs per condition dir: annotated BAM (truth qname + ZF bit 0x04 + ZL locus tag),
loci.tsv, and the reference GTF (transcript strands).

Usage:
  python scripts/debug/em_strand_locus_diag.py [--suite DIR] [--bam-name NAME]
"""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import defaultdict
from pathlib import Path

import pysam

from rigel.sim.truth import parse_origin

AF_GDNA_BIT = 0x04
DEFAULT_SUITE = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"
CONDITIONS = [
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_on",
]


def load_gene_intervals(gtf_path: Path):
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


def overlapping_strands(genes, chrom, start, end) -> set[str]:
    if chrom not in genes:
        return set()
    starts, ivs = genes[chrom]
    hi = bisect.bisect_right(starts, end)
    strands = set()
    for i in range(hi - 1, -1, -1):
        g_start, g_end, g_strand = ivs[i]
        if g_end <= start:
            if g_start < start - 1_000_000:
                break
            continue
        if g_start < end and g_end > start:
            strands.add(g_strand)
    return strands


def load_loci(loci_tsv: Path) -> dict[int, dict]:
    out = {}
    with open(loci_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[int(row["locus_id"])] = row
    return out


def fnum(d, k):
    try:
        return float(d.get(k, 0.0))
    except (TypeError, ValueError):
        return 0.0


def analyze(cond_dir: Path, genes, bam_name: str):
    bam_path = cond_dir / bam_name
    if not bam_path.exists():
        bam_path = cond_dir / "annotated.bam"
    loci = load_loci(cond_dir / "rigel_out" / "loci.tsv")

    # Q1 cross-tab: gDNA-truth fragments overlapping a gene, by (ZL in-locus?) x (co/anti)
    xtab = {"inlocus": {"co": [0, 0], "anti": [0, 0]},  # [leaked, caught]
            "interg": {"co": [0, 0], "anti": [0, 0]}}
    # per-locus accumulation (only fragments with ZL>=0)
    per_locus = defaultdict(lambda: {"co": 0, "anti": 0, "co_leak": 0, "anti_leak": 0})
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
        if org is None or org.kind != "gdna":
            continue
        if not r.has_tag("ZF"):
            continue
        zf = int(r.get_tag("ZF"))
        zl = int(r.get_tag("ZL")) if r.has_tag("ZL") else -1
        strands = overlapping_strands(genes, org.ref, org.start, org.end)
        if not strands or len(strands) > 1:
            continue  # need an unambiguous single-strand gene context
        gd = "+" if org.strand == "f" else "-"
        bucket = "co" if gd in strands else "anti"
        called_gdna = bool(zf & AF_GDNA_BIT)
        idx = 1 if called_gdna else 0  # 1=caught, 0=leaked
        loc = "inlocus" if zl >= 0 else "interg"
        xtab[loc][bucket][idx] += 1
        if zl >= 0:
            pl = per_locus[zl]
            pl[bucket] += 1
            if not called_gdna:
                pl[bucket + "_leak"] += 1
    bam.close()
    return xtab, per_locus, loci


def pct(leaked, caught):
    t = leaked + caught
    return 100.0 * leaked / t if t else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=DEFAULT_SUITE)
    ap.add_argument("--bam-name", default="annotated.baseline_cur.bam")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    suite = Path(args.suite)
    gtf = suite / "reference" / "genes.gtf"
    genes = load_gene_intervals(gtf)

    for cond in CONDITIONS:
        cdir = suite / cond
        if not cdir.exists():
            continue
        xtab, per_locus, loci = analyze(cdir, genes, args.bam_name)
        print(f"\n{'='*90}\n=== {cond} ===\n{'='*90}")
        print("Q1 — where gDNA-truth (gene-overlapping) fragments live & their leak:")
        print(f"  {'':>10} {'co-stranded (n / leak%)':>28} {'antisense (n / leak%)':>28}")
        for loc in ("inlocus", "interg"):
            co, anti = xtab[loc]["co"], xtab[loc]["anti"]
            print(f"  {loc:>10} {sum(co):>10d} / {pct(*co):5.1f}%        "
                  f"{sum(anti):>10d} / {pct(*anti):5.1f}%")
        # anchor question: of antisense gDNA, what frac is intergenic (ZL=-1)?
        anti_in = sum(xtab["inlocus"]["anti"])
        anti_ig = sum(xtab["interg"]["anti"])
        co_in = sum(xtab["inlocus"]["co"])
        if anti_in + anti_ig:
            print(f"  --> antisense gDNA routed to intergenic: "
                  f"{100.0*anti_ig/(anti_in+anti_ig):.1f}%  "
                  f"(anchor {'ABSENT' if anti_ig > anti_in else 'present'} in loci)")

        # Q2 — top leaky loci: budget vs truth vs assignment
        rows = []
        for zl, pl in per_locus.items():
            true_gdna_inlocus = pl["co"] + pl["anti"]
            leak = pl["co_leak"] + pl["anti_leak"]
            L = loci.get(zl, {})
            rows.append((leak, zl, pl, true_gdna_inlocus, L))
        rows.sort(reverse=True)
        print(f"\nQ2 — top {args.top} loci by co-stranded gDNA leak "
              f"(truth counts are IN-LOCUS only, ZL-matched):")
        print(f"  {'locus':>6} {'mrna_em':>9} {'gdna_em':>8} {'gdna_prior':>10} "
              f"{'true_gd_in':>10} {'co':>6} {'anti':>6} {'co_leak':>7} "
              f"{'budget/true':>11} {'em/true':>8}")
        for leak, zl, pl, true_in, L in rows[:args.top]:
            mrna = fnum(L, "mrna")
            gdna_em = fnum(L, "gdna")
            prior = fnum(L, "gdna_prior_count")
            br = prior / true_in if true_in else 0.0
            er = gdna_em / true_in if true_in else 0.0
            print(f"  {zl:>6d} {mrna:>9.0f} {gdna_em:>8.0f} {prior:>10.1f} "
                  f"{true_in:>10d} {pl['co']:>6d} {pl['anti']:>6d} {pl['co_leak']:>7d} "
                  f"{br:>11.2f} {er:>8.2f}")


if __name__ == "__main__":
    main()
