#!/usr/bin/env python3
"""Generate whole-genome simulated reads from an existing FASTA/GTF config."""

from __future__ import annotations

import argparse
import logging
import sys
import time

from rigel.sim.whole_genome import parse_yaml_config, run_simulation


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Whole-genome RNA-seq read simulator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", required=True, help="YAML configuration file")
    p.add_argument("--genome", help="Genome FASTA (overrides YAML)")
    p.add_argument("--gtf", help="Gene annotation GTF (overrides YAML)")
    p.add_argument("--outdir", help="Output directory (overrides YAML)")
    p.add_argument("--n-rna", type=int, help="Number of RNA fragments (overrides YAML)")
    p.add_argument(
        "-j",
        "--n-workers",
        type=int,
        default=None,
        help="Worker processes for parallel read generation (overrides YAML)",
    )
    p.add_argument("--no-oracle", action="store_true", help="Skip oracle BAM generation")
    p.add_argument("--verbose", action="store_true", default=None, help="Verbose logging")
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    cfg = parse_yaml_config(args.config)

    # CLI overrides
    if args.genome:
        cfg.genome = args.genome
    if args.gtf:
        cfg.gtf = args.gtf
    if args.outdir:
        cfg.outdir = args.outdir
    if args.n_rna is not None:
        cfg.simulation.n_rna_fragments = args.n_rna
    if args.n_workers is not None:
        cfg.simulation.n_workers = args.n_workers
    if args.no_oracle:
        cfg.oracle_bam = False
    if args.verbose is not None:
        cfg.verbose = args.verbose

    level = logging.DEBUG if cfg.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    if not cfg.genome:
        print("Error: genome FASTA not specified", file=sys.stderr)
        return 1
    if not cfg.gtf:
        print("Error: GTF not specified", file=sys.stderr)
        return 1

    print("Whole-genome RNA-seq simulator", flush=True)
    print(f"  Genome:           {cfg.genome}", flush=True)
    print(f"  GTF:              {cfg.gtf}", flush=True)
    print(f"  Output:           {cfg.outdir}", flush=True)
    print(f"  RNA fragments:    {cfg.simulation.n_rna_fragments:,}", flush=True)
    print(f"  Workers:          {cfg.simulation.n_workers}", flush=True)
    print(f"  gDNA rates:       {cfg.gdna.rates}", flush=True)
    print(f"  Strand specs:     {cfg.strand_specificities}", flush=True)
    # Print the configured nRNA mode, whichever it is. Whether the nascent sweep is actually run
    # cannot be reported here: it is decided later, by whether the loaded abundance file carried an
    # `nrna_abundance` column, which is why the file case prints a conditional line of its own.
    if cfg.nrna.mode == "additive_ratio":
        print(f"  nRNA ratios:      {cfg.nrna.ratios}", flush=True)
    elif cfg.nrna.mode == "sparse":
        print(
            f"  nRNA ranges:      {cfg.nrna.abundance_ranges} (LOG-uniform, absolute)", flush=True
        )
        print(f"  nRNA on_fraction: {cfg.nrna.on_fraction} of gene spans", flush=True)
    if cfg.abundance.mode == "file":
        print(
            "  nRNA:             the sweep above is SKIPPED if the abundance file supplies "
            "nrna_abundance",
            flush=True,
        )
    print(f"  Oracle BAM:       {cfg.oracle_bam}", flush=True)

    t0 = time.monotonic()
    conditions = run_simulation(cfg)
    elapsed = time.monotonic() - t0

    print(f"\nSimulation complete in {elapsed:.1f}s", flush=True)
    print(f"  {len(conditions)} conditions generated → {cfg.outdir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())