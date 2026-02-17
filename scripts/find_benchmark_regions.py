#!/usr/bin/env python3
"""Find diverse gene-dense regions for benchmarking."""
import gzip
from collections import defaultdict

GTF = "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz"

# Parse genes from GTF
genes = defaultdict(list)
transcripts_per_gene = defaultdict(int)

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        if not chrom.startswith("chr"):
            continue
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]

        gene_name = gene_type = gene_id = tid = ""
        for a in attrs.split(";"):
            a = a.strip()
            if a.startswith('gene_name "'):
                gene_name = a.split('"')[1]
            elif a.startswith('gene_type "'):
                gene_type = a.split('"')[1]
            elif a.startswith('gene_id "'):
                gene_id = a.split('"')[1]
            elif a.startswith('transcript_id "'):
                tid = a.split('"')[1]

        if feature == "gene" and gene_type == "protein_coding":
            genes[chrom].append((start, end, gene_name, strand))
        elif feature == "transcript" and gene_type == "protein_coding":
            transcripts_per_gene[gene_id] += 1

# Curated well-known regions for benchmarking diversity
curated = [
    # FGFR2: complex alternative splicing, many isoforms
    ("chr10", 121476640, 121601584, "FGFR2"),
    # BRCA1: long gene, many exons
    ("chr17", 43044295, 43170245, "BRCA1"),
    # HBB cluster: multiple overlapping globin genes
    ("chr11", 5225000, 5310000, "HBB_cluster"),
    # TP53: tumor suppressor, moderate complexity
    ("chr17", 7668402, 7687550, "TP53"),
    # EGFR: receptor tyrosine kinase, many isoforms
    ("chr7", 55019017, 55211628, "EGFR"),
    # CD44: extreme alternative splicing
    ("chr11", 35138870, 35232402, "CD44"),
    # HOXA cluster: tandem homeobox genes, overlapping
    ("chr7", 27090070, 27254935, "HOXA_cluster"),
    # MYC: oncogene, small gene + neighbors
    ("chr8", 127735434, 127842951, "MYC_locus"),
    # ACTB: housekeeping, high expression baseline
    ("chr7", 5527148, 5563902, "ACTB"),
    # GAPDH: housekeeping
    ("chr12", 6534512, 6538371, "GAPDH"),
]

# Check which have transcripts in the GTF
print("# Curated regions and contained gene counts:")
print("# chrom\tstart\tend\tlabel\tcontained_genes")
for chrom, start, end, label in curated:
    if chrom not in genes:
        print(f"# SKIP {label}: {chrom} not in GTF")
        continue
    contained = [
        (s, e, n, st) for s, e, n, st in genes[chrom]
        if s >= start and e <= end
    ]
    if contained:
        gene_names = ", ".join(n for _, _, n, _ in contained[:6])
        strands = set(st for _, _, _, st in contained)
        print(f"{chrom}\t{start}\t{end}\t{label}\t# {len(contained)} genes ({gene_names}), strands={strands}")
    else:
        # Try expanding
        nearby = [
            (s, e, n, st) for s, e, n, st in genes[chrom]
            if abs(s - start) < 200000 or abs(e - end) < 200000
        ]
        if nearby:
            min_s = min(s for s, _, _, _ in nearby)
            max_e = max(e for _, e, _, _ in nearby)
            print(f"# {label}: no genes fully contained, nearest at {min_s}-{max_e}")
        else:
            print(f"# {label}: no nearby genes found")

# Also find high-isoform-count loci automatically
print("\n# Auto-discovered dense regions (>= 5 protein-coding genes in 150kb):")
for chrom in sorted(genes.keys()):
    if "_" in chrom:
        continue
    g_list = sorted(genes[chrom])
    if len(g_list) < 5:
        continue
    for i in range(len(g_list)):
        s0 = g_list[i][0]
        window_end = s0 + 150000
        contained = [(s, e, n, st) for s, e, n, st in g_list if s >= s0 and e <= window_end]
        if len(contained) >= 8:
            max_e = max(e for _, e, _, _ in contained)
            gene_names = ", ".join(n for _, _, n, _ in contained[:8])
            strands = set(st for _, _, _, st in contained)
            print(f"{chrom}\t{s0}\t{max_e + 5000}\t{chrom}_dense_{s0}\t# {len(contained)} genes ({gene_names}), strands={strands}")
            break
