"""Build the mini-genome reference bundle.

Steps
-----
1. Build a controlled mini genome (FASTA + GTF + blocks.tsv) via
   :mod:`build_mini`.
2. Build a STAR index for the mini genome.
3. Run ``alignable compute --aligner star`` (in the ``alignable``
   conda env) on the mini genome to produce a per-position
   mappability zarr.
4. Build a rigel index that incorporates the mappability zarr.

The output directory ends up with::

    <out>/genome.fa              # genome FASTA (+ .fai)
    <out>/genes.gtf              # mini annotation
    <out>/blocks.tsv             # truth block table (for analysis)
    <out>/alignable/…            # alignable zarr store
    <out>/rigel_index/…          # rigel TranscriptIndex (with mappability)

Usage
-----
::

    python -m scripts.calibration.stress_mini.build_reference \\
        --out /tmp/mini_stress --seed 42
"""
from __future__ import annotations

import argparse
import logging
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

from rigel.index import TranscriptIndex

from .build_mini import build_mini_genome

logger = logging.getLogger("build_reference")


# alignable/STAR read-length bins — keep in sync with the simulator
READ_LENGTHS = "100"       # simulator uses read_length=100
SJDB_OVERHANG = 99         # read_length − 1
FRAG_LEN_MEAN = 200
FRAG_LEN_SD = 40
FRAG_LEN_MIN = 80
FRAG_LEN_MAX = 500

# STAR parameter file — match the reference benchmark config so
# mappability reflects the same aligner behaviour we run at scale.
STAR_PARAMS_SRC = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/alignable/star_arriba_parameters.txt"
)

# Conda env containing STAR + alignable
ALIGNABLE_ENV = "alignable"


def write_genome_and_gtf(out: Path, seed: int,
                         n_transcripts: int,
                         n_unique_intergenic: int) -> None:
    """Step 1: build + write genome.fa, genes.gtf, blocks.tsv."""
    genome, builder, _blocks, _specs, block_df = build_mini_genome(
        n_transcripts=n_transcripts,
        n_unique_intergenic=n_unique_intergenic,
        seed=seed,
    )

    out.mkdir(parents=True, exist_ok=True)

    fasta_path = genome.write_fasta(out)
    if fasta_path.name != "genome.fa":
        # MutableGenome names it ``{genome.name}.fa`` — rename.
        target = out / "genome.fa"
        fasta_path.rename(target)
        # Re-faidx
        subprocess.run(["samtools", "faidx", str(target)], check=True)
        fasta_path = target

    gtf_path = builder.write_gtf(out)
    if gtf_path.name != "genes.gtf":
        gtf_path.rename(out / "genes.gtf")
        gtf_path = out / "genes.gtf"

    block_df.to_csv(out / "blocks.tsv", sep="\t", index=False)
    logger.info("Wrote genome (%d bp), GTF (%d lines), blocks.tsv (%d rows)",
                len(genome), sum(1 for _ in open(gtf_path)), len(block_df))


def build_star_index(out: Path, threads: int) -> Path:
    """Step 2: build a STAR genome index for the mini genome.

    Uses a small-genome-safe ``--genomeSAindexNbases`` per STAR's
    recommendation: ``min(14, log2(L)/2 - 1)``.
    """
    star_dir = out / "star_index"
    if star_dir.exists():
        shutil.rmtree(star_dir)
    star_dir.mkdir(parents=True)

    genome_len = (out / "genome.fa.fai").read_text().split("\t")[1]
    genome_len = int(genome_len)
    sa_index = min(14, int(math.log2(genome_len) / 2 - 1))

    cmd = [
        "conda", "run", "-n", ALIGNABLE_ENV, "--no-capture-output",
        "STAR",
        "--runMode", "genomeGenerate",
        "--runThreadN", str(threads),
        "--genomeDir", str(star_dir),
        "--genomeFastaFiles", str(out / "genome.fa"),
        "--sjdbGTFfile", str(out / "genes.gtf"),
        "--sjdbOverhang", str(SJDB_OVERHANG),
        "--genomeSAindexNbases", str(sa_index),
    ]
    logger.info("Building STAR index (L=%d, SAindex=%d): %s",
                genome_len, sa_index, " ".join(cmd))
    subprocess.run(cmd, check=True)
    return star_dir


def run_alignable(out: Path, star_dir: Path,
                  threads: int, parallel: int) -> Path:
    """Step 3: run ``alignable compute --aligner star``."""
    aln_dir = out / "alignable"
    if aln_dir.exists():
        shutil.rmtree(aln_dir)

    if not STAR_PARAMS_SRC.exists():
        raise FileNotFoundError(
            f"STAR parameter file not found: {STAR_PARAMS_SRC}"
        )

    cmd = [
        "conda", "run", "-n", ALIGNABLE_ENV, "--no-capture-output",
        "alignable", "compute",
        "--genome", str(out / "genome.fa"),
        "--aligner", "star",
        "--aligner-index", str(star_dir),
        "--aligner-config", str(STAR_PARAMS_SRC),
        "--read-lengths", READ_LENGTHS,
        "--frag-len-mean", str(FRAG_LEN_MEAN),
        "--frag-len-sd", str(FRAG_LEN_SD),
        "--frag-len-min", str(FRAG_LEN_MIN),
        "--frag-len-max", str(FRAG_LEN_MAX),
        "--error-rate", "0.0",     # noise-free — isolates mappability from seq errors
        "--coverage", "5",
        "--tolerance", "3",
        "--threads", str(threads),
        "--parallel", str(parallel),
        "--output-dir", str(aln_dir),
    ]
    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # ``alignable compute`` writes the alignable output *directory* at
    # <aln_dir>, with a ``mappability.zarr`` subdir inside it.  Rigel's
    # index builder expects the parent directory (the alignable store).
    if not (aln_dir / "mappability.zarr").exists():
        raise FileNotFoundError(
            f"Expected alignable store layout at {aln_dir} "
            "(missing mappability.zarr subdirectory)"
        )
    logger.info("alignable store at %s", aln_dir)
    return aln_dir


def build_rigel_index(out: Path, alignable_zarr: Path, read_length: int) -> Path:
    """Step 3: build rigel TranscriptIndex with mappability."""
    index_dir = out / "rigel_index"
    if index_dir.exists():
        shutil.rmtree(index_dir)
    logger.info("Building rigel index with mappability (read_length=%d)…", read_length)
    TranscriptIndex.build(
        fasta_file=out / "genome.fa",
        gtf_file=out / "genes.gtf",
        output_dir=index_dir,
        alignable_zarr_path=alignable_zarr,
        mappability_read_length=read_length,
        write_tsv=False,
    )
    return index_dir


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, required=True,
                   help="Output root directory for the reference bundle.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-transcripts", type=int, default=60)
    p.add_argument("--n-unique-intergenic", type=int, default=20)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--parallel", type=int, default=4)
    p.add_argument("--read-length", type=int, default=100,
                   help="Read length for rigel index mappability query.")
    p.add_argument("--skip-alignable", action="store_true",
                   help="(debug) build genome + index without mappability.")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    out = args.out.resolve()
    write_genome_and_gtf(
        out, seed=args.seed,
        n_transcripts=args.n_transcripts,
        n_unique_intergenic=args.n_unique_intergenic,
    )
    if args.skip_alignable:
        logger.warning("Skipping alignable — index will have no mappability.")
        TranscriptIndex.build(
            fasta_file=out / "genome.fa",
            gtf_file=out / "genes.gtf",
            output_dir=out / "rigel_index",
            write_tsv=False,
        )
        return 0

    star_dir = build_star_index(out, threads=args.threads)
    zarr_path = run_alignable(out, star_dir, threads=args.threads, parallel=args.parallel)
    build_rigel_index(out, zarr_path, read_length=args.read_length)
    logger.info("Reference bundle ready in %s", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
