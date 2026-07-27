#!/usr/bin/env python3
"""Build the exact-geometry toy reference for the gDNA effective-length study.

A 1 Mb single-chromosome genome with **10 genes**, each a single transcript of
**two 1 kb exons separated by a 20 kb intron** (genomic span 22 kb), spaced 100 kb
apart on the + strand. This is the deliberately trivial scenario for deriving and
making contractual the IPR / exposure-weighted gDNA effective length under a
hybrid-capture energy sweep.

Why a custom builder (not ``simulate_suite.py``'s generator): the random generator
caps introns at 10 kb and cannot emit an exact 2-exon / 1 kb / 20 kb structure. We
reuse the suite's own ``GeneDef``/``TranscriptDef``/``ExonDef`` + writers so the
emitted ``reference/genome.fa`` + ``reference/genes.gtf`` are byte-identical in
format to a normal suite reference; ``simulate_suite.py --skip-existing`` then reuses
them (designs capture probes on these exact exons, runs the sweep).

Usage::

    python scripts/sim/build_toy_2exon_reference.py --outdir <suite_dir>
    python scripts/sim/simulate_suite.py --config scripts/sim/configs/toy_1mb_2exon_capture_sweep.yaml --skip-existing
    python scripts/sim/bench_calibration.py --sim-base <suite_dir> --run --force

The geometry constants below are the contractual ground truth (exon footprint
2 kb/gene, full span 22 kb/gene → ~11x eff-len contraction at perfect capture).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from rigel.sim.synthetic_genome import (
    REF_NAME,
    ExonDef,
    GeneDef,
    TranscriptDef,
    generate_random_genome,
    inject_splice_sites,
    write_fasta,
    write_gtf,
)

GENOME_LENGTH = 1_000_000
N_GENES = 10
GENE_PITCH = 100_000  # one gene per 100 kb slot
GENE_OFFSET = 30_000  # exon1 start within each slot
EXON_LEN = 1_000
INTRON_LEN = 20_000
SEED = 20260605


def build_genes() -> list[GeneDef]:
    """10 genes, each two 1 kb exons + a 20 kb intron, + strand, single isoform."""
    genes: list[GeneDef] = []
    for i in range(N_GENES):
        base = i * GENE_PITCH + GENE_OFFSET
        e1 = ExonDef(start=base, end=base + EXON_LEN)
        e2 = ExonDef(start=base + EXON_LEN + INTRON_LEN, end=base + 2 * EXON_LEN + INTRON_LEN)
        gene_id = f"GENE{i + 1:04d}"
        tx = TranscriptDef(
            t_id=f"{gene_id}.1",
            gene_id=gene_id,
            gene_name=f"Gene{i + 1}",
            strand="+",
            exons=[e1, e2],
        )
        genes.append(
            GeneDef(
                gene_id=gene_id,
                gene_name=f"Gene{i + 1}",
                strand="+",
                start=e1.start,
                end=e2.end,
                transcripts=[tx],
            )
        )
    return genes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True, help="suite dir (reference/ written under it)")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    ref_dir = args.outdir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    genes = build_genes()
    all_tx = [tx for g in genes for tx in g.transcripts]

    seq_array = np.frombuffer(generate_random_genome(GENOME_LENGTH, args.seed).encode("ascii"),
                              dtype=np.uint8).copy()
    inject_splice_sites(seq_array, all_tx)
    genome_seq = seq_array.tobytes().decode("ascii")

    fasta_path = write_fasta(genome_seq, REF_NAME, ref_dir)
    gtf_path = write_gtf(genes, REF_NAME, ref_dir)

    print(f"  Genome:      {fasta_path} ({GENOME_LENGTH:,} bp, ref={REF_NAME})")
    print(f"  Annotation:  {gtf_path}")
    print(f"  Genes:       {len(genes)} x (2 exons {EXON_LEN} bp + {INTRON_LEN} bp intron, "
          f"span {2 * EXON_LEN + INTRON_LEN:,} bp; exonic {2 * EXON_LEN:,} bp)")
    print(f"  Contract:    full span {2 * EXON_LEN + INTRON_LEN:,} bp/gene; "
          f"exon footprint {2 * EXON_LEN:,} bp/gene "
          f"(~{(2 * EXON_LEN + INTRON_LEN) / (2 * EXON_LEN):.0f}x eff-len contraction at perfect capture)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
