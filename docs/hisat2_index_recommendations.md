# HISAT2 Index & Alignment Recommendations for hulkrna

## TL;DR

**Build a plain HISAT2 index** (no `--ss` / `--exon` at build time).
Pass `--known-splicesite-infile` **at alignment time** only.
No other parameter changes are needed beyond the defaults used by the
benchmark pipeline (`-k 50 --secondary --no-unal`).

## Background

HISAT2 supports two index modes:

| Mode | Build flags | Description |
|------|-------------|-------------|
| **HGFM** (graph) | `--ss splice_sites.txt --exon exons.txt` | Embeds known splice junctions into the FM-index graph |
| **Plain** (BWT) | *(none)* | Standard BWT index; junctions discovered at alignment time |

The HGFM mode is the recommended default in the HISAT2 documentation and
generally performs well.  However, we discovered a severe failure mode in
regions with **paralogous genes sharing exonic sequence** (e.g., the HBB
locus on chr11).

## The Problem: HGFM Graph Bias

When the HGFM index embeds splice sites from closely related genes, the
graph search preferentially traverses **shorter-intron** junctions first.
If two genes share exonic sequence but have introns of very different
lengths, reads that should map across a long intron are instead routed
through the shorter junction — and the correct alignment is never explored.

### Case Study: ENST00000642908 (HBB Locus)

- **ENST00000642908** (gene ENSG00000284931): 3 exons, 5.8 kb intron
- **HBG1 / HBG2**: paralogous genes with shared exon sequence and
  sub-1 kb introns

With HGFM, only **1.0%** (152/15,105) of true ENST00000642908 reads
produce an alignment containing the correct 5.8 kb intron junction —
in *any* alignment (primary or secondary).  The remaining 99% are
captured exclusively by HBG1/HBG2 junctions.

## The Solution: Plain Index

Switching to a plain BWT index (no `--ss`/`--exon` at build time)
removes the graph bias entirely.  HISAT2 discovers junctions *de novo*
during alignment and evaluates them on equal footing regardless of intron
length.

With the plain index, **66.2%** (10,014/15,129) of true ENST00000642908
reads produce an alignment with the correct junction — a **66× improvement**.

The remaining 33.8% are reads that land entirely within shared exonic
regions and have no splice junction to disambiguate them.  These are
correctly handled by hulkrna's multimap EM algorithm.

## Parameter Sweep Results

We tested 20 configurations combining the plain index with various
alignment-time parameters:

| Parameter | Options tested | Effect on recovery |
|-----------|---------------|-------------------|
| `--rna-strandness RF` | on / off | No change |
| `--pen-canintronlen` | `G,0,0` / `G,-4,1` / `G,-8,0.5` / default | No change |
| `--no-templatelen-adjustment` | on / off | No change |
| `--no-temp-splicesite` | on / off | No change |
| `--pen-cansplice` | 0 / 3 | No change |
| `--min-intronlen` | 20 (default) / 50 | No change |
| `-k` | 10 / 50 / 100 | No change |
| `--known-splicesite-infile` | on / off | Marginal (+36 reads, +0.2%) |

**Conclusion:** The only parameter that matters is the index type.
Alignment-time parameters have negligible effect once the plain index
is used.  Providing `--known-splicesite-infile` at alignment time gives
a small benefit and is recommended.

## Global Alignment Metrics

| Metric | HGFM | Plain |
|--------|------|-------|
| Overall alignment rate | 99.98% | 100.00% |
| Primary alignments | 99,976 | 100,000 |
| Spliced primary alignments | 52,664 | 53,044 |
| Multimapping primary reads | 10,444 | 33,746 |

The plain index produces more multimapping reads (expected since splice
junctions aren't baked into the graph).  hulkrna's multimap-aware EM
algorithm is designed to resolve these correctly.

## Recommended HISAT2 Invocation

### Index building

```bash
hisat2-build -p <threads> <region.fa> <index_prefix>
```

Do **not** pass `--ss` or `--exon`.

### Alignment

```bash
hisat2 \
  -x <index_prefix> \
  -1 reads_R1.fq -2 reads_R2.fq \
  -p <threads> \
  -k 50 \
  --secondary \
  --no-unal \
  --known-splicesite-infile <splice_sites.txt>
```

## Diagnostic Scripts

The analysis supporting these recommendations is in the following scripts:

- `scripts/_diag_hisat2_options.py` — Round 1: alignment-time params with HGFM index (all 1.0%)
- `scripts/_diag_hisat2_options2.py` — Round 2: plain vs HGFM index (66.2% vs 1.0%)
- `scripts/_diag_hisat2_sweep.py` — Round 3: 20-configuration plain-index sweep (all ~66.2%)

---

*Last updated: 2025-06-24*
