#!/usr/bin/env python3
"""Find benchmark loci from a GTF file.

Identifies regions around target genes, expanding to capture nearby/overlapping
genes for realistic benchmarking. Outputs region specs suitable for the
benchmark YAML config.
"""
import gzip
import sys
from collections import defaultdict

GTF_PATH = "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz"

# Target anchor genes for benchmark regions
ANCHOR_GENES = [
    "EGFR", "FGFR1", "FGFR2", "FGFR3",
    "HOXA1", "HOXB1", "HOXC4", "HOXD1",
    "HBB", "HBA1", "HBA2",
    "MALAT1", "NEAT1",
    "MYC", "PVT1",
    "TP53", "BRCA1", "BRCA2",
    "KRAS", "BRAF",
    "CDK4", "CDK6", "CCND1",
    "BCL2", "BCL6",
    "RB1", "PTEN",
    "CDKN2A", "CDKN2B",
    "TERT",
    "ALK", "ROS1",
    "ERBB2", "ESR1",
    "XIST", "TSIX",
    "H19", "IGF2",
    "SOX2", "NANOG",
    "GAPDH", "ACTB",
    "FMR1",
    "DMD",
    "SNHG1", "GAS5",
    "HOTAIR",
    "APOB", "TTN",
]


def main():
    # Step 1: Parse all genes from GTF
    print("Parsing GTF...", file=sys.stderr)
    all_genes = {}  # gene_name -> (chrom, start, end, strand, gene_id, gene_type)
    gene_by_chrom = defaultdict(list)  # chrom -> [(start, end, gene_name, strand)]

    with gzip.open(GTF_PATH, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            attrs = parts[8]
            gene_name = None
            gene_id = None
            gene_type = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                elif attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1]
                elif attr.startswith("gene_type") or attr.startswith("gene_biotype"):
                    gene_type = attr.split('"')[1]

            if gene_name:
                all_genes[gene_name] = (chrom, start, end, strand, gene_id, gene_type)
                gene_by_chrom[chrom].append((start, end, gene_name, strand))

    # Step 2: Report target genes found
    print("\n=== Target genes found ===", file=sys.stderr)
    target_set = set(ANCHOR_GENES)
    found = {name: all_genes[name] for name in ANCHOR_GENES if name in all_genes}
    missing = [name for name in ANCHOR_GENES if name not in all_genes]
    if missing:
        print(f"Missing: {missing}", file=sys.stderr)

    for name in sorted(found.keys()):
        chrom, start, end, strand, gene_id, gene_type = found[name]
        print(f"  {name}\t{chrom}:{start}-{end}\t{strand}\t{gene_type}\t{end-start}bp", file=sys.stderr)

    # Step 3: Design benchmark regions
    # For each target gene, expand to capture overlapping/nearby genes
    # with padding, then merge nearby regions sharing a label
    print("\n=== Designing benchmark regions ===", file=sys.stderr)

    # Sort genes on each chromosome by position
    for chrom in gene_by_chrom:
        gene_by_chrom[chrom].sort()

    def find_genes_in_region(chrom, start, end):
        """Find all genes overlapping a region."""
        genes = []
        for gs, ge, gn, gstr in gene_by_chrom.get(chrom, []):
            if ge > start and gs < end:
                genes.append((gs, ge, gn, gstr))
        return genes

    def expand_region(chrom, start, end, padding=50000):
        """Expand region to capture overlapping genes + padding."""
        expanded_start = max(1, start - padding)
        expanded_end = end + padding

        # Find all genes overlapping the expanded region
        genes = find_genes_in_region(chrom, expanded_start, expanded_end)
        if genes:
            region_start = max(1, min(g[0] for g in genes) - 10000)
            region_end = max(g[1] for g in genes) + 10000
        else:
            region_start = expanded_start
            region_end = expanded_end

        return region_start, region_end, genes

    # Define the 20 benchmark regions with labels
    regions = []

    # 1. PVT1/MYC (chr8) - existing benchmark
    if "MYC" in found and "PVT1" in found:
        chrom = found["MYC"][0]
        s = min(found["MYC"][1], found["PVT1"][1])
        e = max(found["MYC"][2], found["PVT1"][2])
        rs, re, genes = expand_region(chrom, s, e, padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("PVT1_MYC", chrom, rs, re, gene_names))

    # 2. EGFR (chr7) - large gene, many isoforms
    if "EGFR" in found:
        chrom, s, e = found["EGFR"][0], found["EGFR"][1], found["EGFR"][2]
        rs, re, genes = expand_region(chrom, s, e, padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("EGFR", chrom, rs, re, gene_names))

    # 3. MALAT1/NEAT1 (chr11) - long non-coding RNAs, overlapping
    if "MALAT1" in found and "NEAT1" in found:
        chrom = found["MALAT1"][0]
        s = min(found["MALAT1"][1], found["NEAT1"][1])
        e = max(found["MALAT1"][2], found["NEAT1"][2])
        rs, re, genes = expand_region(chrom, s, e, padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("MALAT1_NEAT1", chrom, rs, re, gene_names))

    # 4. HBB cluster (chr11 - beta-globin locus)
    if "HBB" in found:
        chrom = found["HBB"][0]
        rs, re, genes = expand_region(chrom, found["HBB"][1], found["HBB"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("HBB_cluster", chrom, rs, re, gene_names))

    # 5. HOXA cluster (chr7) - tightly packed, same-strand homeobox genes
    if "HOXA1" in found:
        chrom = found["HOXA1"][0]
        # Expand to capture the whole HOX cluster
        rs, re, genes = expand_region(chrom, found["HOXA1"][1], found["HOXA1"][2], padding=200000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("HOXA_cluster", chrom, rs, re, gene_names))

    # 6. HOXD cluster (chr2)
    if "HOXD1" in found:
        chrom = found["HOXD1"][0]
        rs, re, genes = expand_region(chrom, found["HOXD1"][1], found["HOXD1"][2], padding=200000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("HOXD_cluster", chrom, rs, re, gene_names))

    # 7. TP53 (chr17) - tumor suppressor, many isoforms
    if "TP53" in found:
        chrom = found["TP53"][0]
        rs, re, genes = expand_region(chrom, found["TP53"][1], found["TP53"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("TP53", chrom, rs, re, gene_names))

    # 8. BRCA1 (chr17) - large gene with antisense overlaps
    if "BRCA1" in found:
        chrom = found["BRCA1"][0]
        rs, re, genes = expand_region(chrom, found["BRCA1"][1], found["BRCA1"][2], padding=50000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("BRCA1", chrom, rs, re, gene_names))

    # 9. KRAS (chr12)
    if "KRAS" in found:
        chrom = found["KRAS"][0]
        rs, re, genes = expand_region(chrom, found["KRAS"][1], found["KRAS"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("KRAS", chrom, rs, re, gene_names))

    # 10. CDKN2A/CDKN2B (chr9) - overlapping antisense + readthrough
    if "CDKN2A" in found:
        chrom = found["CDKN2A"][0]
        s = found["CDKN2A"][1]
        e = found["CDKN2A"][2]
        if "CDKN2B" in found and found["CDKN2B"][0] == chrom:
            s = min(s, found["CDKN2B"][1])
            e = max(e, found["CDKN2B"][2])
        rs, re, genes = expand_region(chrom, s, e, padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("CDKN2A_2B", chrom, rs, re, gene_names))

    # 11. ERBB2/HER2 (chr17) - amplified in breast cancer
    if "ERBB2" in found:
        chrom = found["ERBB2"][0]
        rs, re, genes = expand_region(chrom, found["ERBB2"][1], found["ERBB2"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("ERBB2", chrom, rs, re, gene_names))

    # 12. FGFR1 (chr8) - receptor tyrosine kinase, fusions
    if "FGFR1" in found:
        chrom = found["FGFR1"][0]
        rs, re, genes = expand_region(chrom, found["FGFR1"][1], found["FGFR1"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("FGFR1", chrom, rs, re, gene_names))

    # 13. BCL2 (chr18) - apoptosis regulator
    if "BCL2" in found:
        chrom = found["BCL2"][0]
        rs, re, genes = expand_region(chrom, found["BCL2"][1], found["BCL2"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("BCL2", chrom, rs, re, gene_names))

    # 14. PTEN (chr10) - tumor suppressor
    if "PTEN" in found:
        chrom = found["PTEN"][0]
        rs, re, genes = expand_region(chrom, found["PTEN"][1], found["PTEN"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("PTEN", chrom, rs, re, gene_names))

    # 15. H19/IGF2 (chr11) - imprinted locus, antisense
    if "H19" in found and "IGF2" in found:
        chrom = found["H19"][0]
        s = min(found["H19"][1], found["IGF2"][1])
        e = max(found["H19"][2], found["IGF2"][2])
        rs, re, genes = expand_region(chrom, s, e, padding=50000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("H19_IGF2", chrom, rs, re, gene_names))

    # 16. XIST/TSIX (chrX) - X-inactivation locus, antisense pair
    if "XIST" in found:
        chrom = found["XIST"][0]
        s = found["XIST"][1]
        e = found["XIST"][2]
        if "TSIX" in found and found["TSIX"][0] == chrom:
            s = min(s, found["TSIX"][1])
            e = max(e, found["TSIX"][2])
        rs, re, genes = expand_region(chrom, s, e, padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("XIST_TSIX", chrom, rs, re, gene_names))

    # 17. BRAF (chr7) - kinase, many isoforms
    if "BRAF" in found:
        chrom = found["BRAF"][0]
        rs, re, genes = expand_region(chrom, found["BRAF"][1], found["BRAF"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("BRAF", chrom, rs, re, gene_names))

    # 18. TERT (chr5) - telomerase, promoter mutations
    if "TERT" in found:
        chrom = found["TERT"][0]
        rs, re, genes = expand_region(chrom, found["TERT"][1], found["TERT"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("TERT", chrom, rs, re, gene_names))

    # 19. TTN (chr2) - largest human gene
    if "TTN" in found:
        chrom = found["TTN"][0]
        rs, re, genes = expand_region(chrom, found["TTN"][1], found["TTN"][2], padding=50000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("TTN", chrom, rs, re, gene_names))

    # 20. GAPDH/ACTB - housekeeping controls
    if "GAPDH" in found:
        chrom = found["GAPDH"][0]
        rs, re, genes = expand_region(chrom, found["GAPDH"][1], found["GAPDH"][2], padding=100000)
        gene_names = sorted(set(g[2] for g in genes))
        regions.append(("GAPDH", chrom, rs, re, gene_names))

    # Print results
    print("\n=== Benchmark Regions ===")
    for i, (label, chrom, start, end, gene_names) in enumerate(regions, 1):
        span = end - start
        n_genes = len(gene_names)
        top_genes = ", ".join(gene_names[:8])
        if len(gene_names) > 8:
            top_genes += f", ... (+{len(gene_names)-8})"
        print(f"{i:2d}. {label:20s}  {chrom}:{start}-{end}  ({span/1000:.0f}kb, {n_genes} genes)  [{top_genes}]")

    # Print YAML region block
    print("\n=== YAML region block ===")
    print("region:")
    for label, chrom, start, end, _ in regions:
        print(f"  - {label}: {chrom}:{start}-{end}")


if __name__ == "__main__":
    main()
