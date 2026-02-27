#!/usr/bin/env python3
"""Generate abundance TSV for chr8:126,445,441-128,594,175 (PVT1/MYC locus).

For each gene with G_n isoforms:
  - ~50% of isoforms are set to abundance 0 (negative controls)
  - Remaining isoforms get abundance = 2^U(0, 12)  (log-linear, 1 to 4096)

Output: abundance.tsv with columns transcript_id, abundance, gene_id, gene_name, is_negative_control
"""

import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REGION_CHR = "chr8"
REGION_START = 126_445_441   # 1-based
REGION_END   = 128_594_175

GTF_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz")
OUTPUT = Path("/Users/mkiyer/proj/hulkrna/scripts/benchmarking/chr8_pvt1_myc_abundance.tsv")

SEED = 42
NEG_FRAC = 0.5          # fraction of isoforms per gene to zero out
EXP_MIN = 0.0           # min exponent
EXP_MAX = 12.0          # max exponent → abundance ∈ [1, 4096]


def main():
    # ── Parse GTF for transcripts in region ──────────────────────────
    # We use 0-based half-open internally, GTF is 1-based closed
    start0 = REGION_START - 1
    end0 = REGION_END

    tx_gene: dict[str, str] = {}            # transcript_id → gene_id
    tx_gname: dict[str, str] = {}           # transcript_id → gene_name
    tx_exons: dict[str, list[tuple]] = defaultdict(list)

    with gzip.open(GTF_PATH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start_s, end_s, _, strand, _, attrs_str = parts
            if chrom != REGION_CHR or feature != "exon":
                continue
            s0 = int(start_s) - 1
            e0 = int(end_s)
            if s0 < start0 or e0 > end0:
                continue

            # Parse attributes
            tid = gid = gname = ""
            for a in attrs_str.strip(";").split(";"):
                a = a.strip()
                m = re.match(r'(\w+)\s+"([^"]*)"', a)
                if m:
                    if m.group(1) == "transcript_id":
                        tid = m.group(2)
                    elif m.group(1) == "gene_id":
                        gid = m.group(2)
                    elif m.group(1) == "gene_name":
                        gname = m.group(2)
            if tid and gid:
                tx_gene[tid] = gid
                tx_gname[tid] = gname or gid
                tx_exons[tid].append((s0, e0))

    # Filter to transcripts fully within region
    valid_tids = sorted(
        tid for tid, exons in tx_exons.items()
        if all(s >= start0 and e <= end0 for s, e in exons)
    )

    # ── Group transcripts by gene ────────────────────────────────────
    gene_to_txs: dict[str, list[str]] = defaultdict(list)
    for tid in valid_tids:
        gene_to_txs[tx_gene[tid]].append(tid)

    # ── Assign abundances ────────────────────────────────────────────
    rng = np.random.default_rng(SEED)
    rows = []

    for gid in sorted(gene_to_txs):
        txs = gene_to_txs[gid]
        n = len(txs)

        # For genes with only 1 isoform, always keep it expressed
        if n == 1:
            exp = rng.uniform(EXP_MIN, EXP_MAX)
            abundance = 2.0 ** exp
            rows.append((txs[0], abundance, gid, tx_gname[txs[0]], False))
            continue

        # Randomly select ~50% to be negative controls (zero)
        n_neg = max(1, round(n * NEG_FRAC))  # at least 1 negative control
        n_pos = n - n_neg
        if n_pos == 0:
            n_pos = 1
            n_neg = n - 1

        # Shuffle and partition
        shuffled = list(txs)
        rng.shuffle(shuffled)
        neg_set = set(shuffled[:n_neg])

        for tid in txs:
            if tid in neg_set:
                rows.append((tid, 0.0, gid, tx_gname[tid], True))
            else:
                exp = rng.uniform(EXP_MIN, EXP_MAX)
                abundance = 2.0 ** exp
                rows.append((tid, abundance, gid, tx_gname[tid], False))

    # ── Write output ─────────────────────────────────────────────────
    with open(OUTPUT, "w") as f:
        f.write("transcript_id\tabundance\tgene_id\tgene_name\tis_negative_control\n")
        for tid, ab, gid, gname, is_neg in rows:
            f.write(f"{tid}\t{ab:.4f}\t{gid}\t{gname}\t{is_neg}\n")

    # ── Summary ──────────────────────────────────────────────────────
    n_total = len(rows)
    n_neg = sum(1 for _, ab, _, _, _ in rows if ab == 0.0)
    n_pos = n_total - n_neg
    n_genes = len(gene_to_txs)
    pos_abs = [ab for _, ab, _, _, _ in rows if ab > 0]
    print(f"Region: {REGION_CHR}:{REGION_START}-{REGION_END}")
    print(f"Genes: {n_genes}")
    print(f"Transcripts: {n_total} ({n_pos} expressed, {n_neg} negative controls)")
    print(f"Abundance range (expressed): {min(pos_abs):.1f} – {max(pos_abs):.1f}")
    print(f"Abundance median (expressed): {np.median(pos_abs):.1f}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
