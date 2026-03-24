#!/usr/bin/env python3
"""Root cause analysis for top 10 minimap2 regression genes.

For each gene, reports:
1. Transcript coordinates and exon structure from rigel index
2. Alignment fate: where do reads from target transcripts actually land?
3. Classification: pseudogene absorption, isoform switch, micro-exon, PAR, etc.
"""
import sys
import subprocess
import pyarrow.feather as pf
import pandas as pd
import numpy as np
from pathlib import Path

COND = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"
IDX  = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/rigel_index"
MM2_BAM = f"{COND}/align_minimap2/reads_namesort.bam"
ORA_BAM = f"{COND}/align_oracle/reads_namesort.bam"

# Top 10 target transcripts from find_minimap2_regressions.py
TARGETS = [
    ("ENST00000631211.1",  "ENSG00000280800", "ENSG00000280800"),
    ("ENST00000292807.9",  "AP2M1",           "ENSG00000161203.14"),
    ("ENST00000674977.2",  "POLR2A",          "ENSG00000181222.19"),
    ("ENST00000535788.1",  "UBB",             "ENSG00000170315.14"),
    ("ENST00000381401.11", "SLC25A6",         "ENSG00000169100.14"),
    ("ENST00000350320.11", "SEPTIN7",         "ENSG00000122545.22"),
    ("ENST00000546591.6",  "RPL41",           "ENSG00000229117.9"),
    ("ENST00000286448.12", "VAMP7",           "ENSG00000124333.16"),
    ("ENST00000602676.6",  "OAZ1",            "ENSG00000104904.12"),
    ("ENST00000446868.8",  "CYTH1",           "ENSG00000108669.18"),
]

def load_index():
    tx_df = pf.read_feather(f"{IDX}/transcripts.feather").to_pandas()
    return tx_df

def get_tx_info(tx_df, t_id):
    rows = tx_df[tx_df["t_id"] == t_id]
    if len(rows) == 0:
        return None
    return rows.iloc[0]

def format_exons(row):
    """Format exon structure from transcript row."""
    chrom = row.get("chrom", row.get("chr", "?"))
    start = row.get("tx_start", row.get("start", 0))
    end   = row.get("tx_end", row.get("end", 0))
    strand = row.get("strand", "?")
    exon_starts = row.get("exon_starts", None)
    exon_ends   = row.get("exon_ends", None)
    if exon_starts is not None and exon_ends is not None:
        starts = [int(x) for x in str(exon_starts).split(",") if x]
        ends   = [int(x) for x in str(exon_ends).split(",") if x]
        exon_sizes = [e - s for s, e in zip(starts, ends)]
        n_exons = len(exon_sizes)
        return (chrom, start, end, strand, n_exons, exon_sizes, starts, ends)
    return (chrom, start, end, strand, "?", [], [], [])

def cigar_stats(cigar):
    """Parse CIGAR to get N operations (intron spans) and exon bases."""
    import re
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    n_operations = [(int(n), op) for n, op in ops if op == 'N']
    exon_bases = sum(int(n) for n, op in ops if op in ('M', '=', 'X'))
    return n_operations, exon_bases

def run_samtools_view(bam, region=None, flags=None):
    """Run samtools view and return lines."""
    cmd = ["samtools", "view"]
    if flags:
        cmd += ["-F", str(flags)]
    cmd.append(bam)
    if region:
        cmd.append(region)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return result.stdout.strip().split("\n") if result.stdout.strip() else []
    except Exception as e:
        return []

def analyze_alignment_fate(t_id, chrom, start, end, truth_count):
    """For a target transcript, analyze where its reads land in mm2 BAM."""
    region = f"{chrom}:{start}-{end}"
    # Note: bam is name-sorted so we can't do region queries easily
    # Instead use samtools with -F 256 (not secondary) but BAM is name-sorted...
    # We can still get alignments in the region if we accept it may be slow for name-sorted
    lines = run_samtools_view(MM2_BAM, region=region)
    if not lines or lines == ['']:
        return {"n_reads": 0, "region": region}

    n_total = len(lines)
    primary_mapq = []
    ref_positions = []
    cigars = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        flag  = int(parts[1])
        rname = parts[2]
        pos   = int(parts[3])
        mapq  = int(parts[4])
        cigar = parts[5]
        if not (flag & 256):  # primary
            primary_mapq.append(mapq)
            ref_positions.append(pos)
            cigars.append(cigar)

    n_spliced = sum(1 for c in cigars if 'N' in c)
    return {
        "region": region,
        "n_alignments": n_total,
        "n_primary": len(primary_mapq),
        "n_spliced": n_spliced,
        "mapq_mean": np.mean(primary_mapq) if primary_mapq else 0,
        "mapq_0_frac": sum(1 for m in primary_mapq if m == 0) / max(len(primary_mapq), 1),
    }

def main():
    print("Loading rigel index...")
    tx_df = load_index()
    print(f"  Loaded {len(tx_df)} transcripts")
    print(f"  Columns: {list(tx_df.columns[:15])}")

    # Load counts
    oracle_df = pd.read_csv(f"{COND}/per_transcript_counts_oracle.csv")
    mm2_df    = pd.read_csv(f"{COND}/per_transcript_counts_minimap2.csv")
    counts = oracle_df.merge(mm2_df[["transcript_id", "rigel_minimap2"]], on="transcript_id", how="left")

    print("\n" + "=" * 110)
    print("ROOT CAUSE ANALYSIS: TOP 10 REGRESSION GENES")
    print("=" * 110)

    for t_id, gene_name, gene_id in TARGETS:
        print(f"\n{'─'*110}")

        # Get transcript coordinates
        row = get_tx_info(tx_df, t_id)
        if row is not None:
            info = format_exons(row)
            chrom, start, end, strand, n_exons, exon_sizes, ex_starts, ex_ends = info
            # small exons
            tiny_exons = [(i+1, s) for i, s in enumerate(exon_sizes) if s < 60]
            # intron sizes
            introns = []
            for i in range(len(ex_starts)-1):
                introns.append(ex_ends[i+1-1+1] - ex_starts[i+1] if i+1 < len(ex_starts) else 0)
                # Actually:
            introns = [ex_starts[i+1] - ex_ends[i] for i in range(len(ex_starts)-1)]
            tiny_introns = [(i+1, s) for i, s in enumerate(introns) if s < 300]
            coord_str = f"{chrom}:{start}-{end}({strand})"
        else:
            coord_str = "(not found in index)"
            tiny_exons = []
            tiny_introns = []
            chrom = start = end = None
            n_exons = "?"
            exon_sizes = []
            introns = []

        # Get counts for this transcript
        cnt_row = counts[counts["transcript_id"] == t_id]
        if len(cnt_row):
            c = cnt_row.iloc[0]
            truth = c.mrna_truth
            oracle = c.rigel_oracle
            salmon = c.salmon
            kalli  = c.kallisto
            mm2    = c.rigel_minimap2
        else:
            truth = oracle = salmon = kalli = mm2 = 0

        # Get all transcripts in this gene for context
        gene_counts = counts[counts["gene_id"].astype(str) == gene_id].sort_values("mrna_truth", ascending=False)

        print(f"  GENE: {gene_name} ({gene_id})")
        print(f"  TARGET: {t_id}  |  {coord_str}  |  n_exons={n_exons}")
        if exon_sizes:
            print(f"  Exon sizes: {exon_sizes}")
        if introns:
            print(f"  Intron sizes: {introns}")
        if tiny_exons:
            print(f"  *** TINY EXONS (< 60bp): {tiny_exons}")
        if tiny_introns:
            print(f"  *** TINY INTRONS (< 300bp): {tiny_introns}")

        print(f"\n  Counts:  truth={truth:.0f}  oracle={oracle:.0f}  salmon={salmon:.0f}  kallisto={kalli:.0f}  mm2={mm2:.0f}")
        if truth > 0:
            print(f"  Ratios:  oracle={oracle/truth:.3f}  salmon={salmon/truth:.3f}  kallisto={kalli/truth:.3f}  mm2={mm2/truth:.3f}")

        # All transcripts in gene
        if len(gene_counts) > 1:
            print(f"\n  All transcripts in {gene_name}:")
            for _, gr in gene_counts.iterrows():
                marker = " <-- TARGET" if gr.transcript_id == t_id else ""
                print(f"    {gr.transcript_id:42s}  truth={gr.mrna_truth:7.0f}  oracle={gr.rigel_oracle:7.0f}  "
                      f"mm2={gr.rigel_minimap2:7.0f}  ratio={gr.rigel_minimap2/max(gr.mrna_truth,1):.3f}{marker}")

        # Check if there's a companion pseudogene in the BED
        if chrom is not None:
            ann_bed = f"/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/minimap2_annotation.bed"
            try:
                result = subprocess.run(
                    ["awk", f'$1=="{chrom}" && $2>={max(0,start-5000)} && $3<={end+5000}', ann_bed],
                    capture_output=True, text=True, timeout=30
                )
                nearby_genes = result.stdout.strip().split("\n") if result.stdout.strip() else []
                if nearby_genes:
                    uniq_names = set()
                    for ent in nearby_genes:
                        parts = ent.split("\t")
                        if len(parts) >= 4:
                            uniq_names.add(parts[3].split("|")[-1] if "|" in parts[3] else parts[3])
                    print(f"\n  Nearby transcripts in paftools.js BED: {sorted(uniq_names)[:10]}")
            except Exception:
                pass

        # Analyze alignment fate from mm2 BAM (coordinate query works even in name-sorted)
        if chrom is not None:
            fate = analyze_alignment_fate(t_id, chrom, start, end, truth)
            print(f"\n  mm2 BAM alignments at {fate['region']}:")
            print(f"    total={fate.get('n_alignments',0)}  primary={fate.get('n_primary',0)}  "
                  f"spliced={fate.get('n_spliced',0)}  mapq_mean={fate.get('mapq_mean',0):.1f}  "
                  f"mapq0_frac={fate.get('mapq_0_frac',0):.1%}")

        # Root cause diagnosis
        print(f"\n  *** ROOT CAUSE HYPOTHESIS:")
        mm2_ratio = mm2 / max(truth, 1) if truth > 0 else float('nan')
        if mm2_ratio < 0.1:
            if tiny_exons:
                print(f"    MICRO-EXON ANCHOR FAILURE: transcript has exons {[e[1] for e in tiny_exons]}bp — "
                      f"reads spanning these junctions misaligned, MANE absent from candidates")
            else:
                print(f"    NEAR-COMPLETE SIPHONING: mm2_ratio={mm2_ratio:.3f} — "
                      f"reads mapped to pseudogene or other locus; NM indistinguishable")
        elif mm2_ratio < 0.55:
            # Check if there's a para-locus (gene duplication)
            print(f"    PARTIAL SIPHONING (ratio={mm2_ratio:.3f}): likely 2-copy locus (PAR, duplication); "
                  f"50% of reads mapped to paralog/pseudogene; NM penalty insufficient")
        elif mm2_ratio > 1.5:
            print(f"    FRAGMENT GAIN (ratio={mm2_ratio:.3f}): reads from another transcript mismapped here; "
                  f"EM seeding artifact from wrong-mapping fragments")
        else:
            print(f"    PARTIAL LOSS (ratio={mm2_ratio:.3f}): unclear; inspect manually")

        # Check over-counted transcripts in same gene
        overcounted = gene_counts[(gene_counts["mrna_truth"] < 50) & (gene_counts["rigel_minimap2"] > 100)]
        if not overcounted.empty:
            print(f"\n  *** OVER-COUNTED (truth<50, mm2>100) transcripts in same gene:")
            for _, oc in overcounted.iterrows():
                print(f"    {oc.transcript_id}  truth={oc.mrna_truth:.0f}  mm2={oc.rigel_minimap2:.0f}")

if __name__ == "__main__":
    main()
