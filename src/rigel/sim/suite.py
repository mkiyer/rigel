#!/usr/bin/env python3
"""Generate synthetic mini-genome simulation suites.

Steps:
  1. Generate synthetic genome (10 Mb, ~50 genes, ~250 transcripts)
    2. Simulate reads across strand x gDNA x nRNA-ratio conditions
  3. Each condition produces oracle BAM + paired-end FASTQ

Usage:
    python scripts/sim/simulate_suite.py \\
    --outdir /Users/mkiyer/Downloads/rigel_runs/sim_synthetic
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from rigel.sim.capture import CaptureConfig
from rigel.sim.synthetic_genome import (
    ANTISENSE_OVERLAP_FRAC,
    GENOME_LENGTH,
    MAX_ISOFORMS,
    MIN_ISOFORMS,
    N_GENES,
    REF_NAME,
    SEED,
    TARGET_TRANSCRIPTS,
    generate_genes,
    generate_random_genome,
    inject_splice_sites,
    write_fasta,
    write_gtf,
)
from rigel.sim.manifest import condition_dir_name
from rigel.sim.whole_genome import (
    AbundanceConfig,
    GDNASimConfig,
    NRNAConfig,
    WholeGenomeSimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    apply_nrna_ratio,
    assign_random_abundances,
    load_transcripts,
    write_truth_abundances,
    write_manifest,
)
from rigel.transcript import Transcript

logger = logging.getLogger(__name__)

GDNA_RATES = [0.0, 0.1, 0.5, 1.0, 2.0]
GDNA_LABELS = ["none", "low", "med", "equal", "high"]
STRAND_SPECIFICITIES = [0.99, 0.50]
PROFILE_NRNA = {
    "smoke": ([0.0], ["none"]),
    "full": ([0.0, 0.1, 0.5, 1.0], ["none", "low", "med", "high"]),
}


def stable_seed(base_seed: int, *parts: object) -> int:
    """Derive a reproducible 32-bit seed from a base seed and string parts."""
    text = "\0".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") & 0xFFFF_FFFF


def _parse_csv_floats(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def _parse_csv_labels(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


def _as_float_list(value: object) -> list[float]:
    if isinstance(value, str):
        return _parse_csv_floats(value)
    return [float(item) for item in value]  # type: ignore[union-attr]


def _as_label_list(value: object) -> list[str]:
    if isinstance(value, str):
        return _parse_csv_labels(value)
    return [str(item) for item in value]  # type: ignore[union-attr]


def _arg_defaults(parser: argparse.ArgumentParser) -> dict[str, object]:
    return {action.dest: action.default for action in parser._actions}


def _apply_config_defaults(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    values: dict[str, object],
) -> None:
    defaults = _arg_defaults(parser)
    for key, value in values.items():
        if getattr(args, key, None) == defaults.get(key):
            setattr(args, key, value)


def _load_suite_config(path: Path) -> dict[str, object]:
    """Load a simulate_suite YAML config into argparse destination names."""
    with open(path) as handle:
        raw = yaml.safe_load(handle) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"suite config must be a YAML mapping: {path}")

    values: dict[str, object] = {}

    def put(section: str | None, source_key: str, dest: str | None = None) -> None:
        dest_key = dest or source_key
        if section is not None and isinstance(raw.get(section), dict):
            section_map = raw[section]
            if source_key in section_map:
                values[dest_key] = section_map[source_key]
                return
        if source_key in raw:
            values[dest_key] = raw[source_key]

    for key in (
        "outdir",
        "genome_length",
        "n_genes",
        "min_isoforms",
        "max_isoforms",
        "target_transcripts",
        "antisense_fraction",
        "seed",
        "profile",
        "conditions",
        "reference_only",
        "skip_existing",
    ):
        put("suite", key)
    put("suite", "n_rna")
    put("simulation", "n_rna_fragments", "n_rna")
    put("simulation", "n_workers")
    for key in ("frag_mean", "frag_std", "frag_min", "frag_max", "read_length", "error_rate"):
        put("simulation", key)

    put("abundance", "seed", "abundance_seed")
    put("abundance", "min", "abundance_min")
    put("abundance", "max", "abundance_max")
    put("abundance", "frac_expressed")

    put("gdna", "rates", "gdna_rates")
    put("gdna", "rate_labels", "gdna_labels")
    for key in ("gdna_frag_mean", "gdna_frag_std", "gdna_frag_min", "gdna_frag_max"):
        put("suite", key)
    put("gdna", "frag_mean", "gdna_frag_mean")
    put("gdna", "frag_std", "gdna_frag_std")
    put("gdna", "frag_min", "gdna_frag_min")
    put("gdna", "frag_max", "gdna_frag_max")

    put(None, "strand_specificities")
    put("nrna", "ratios", "nrna_ratios")
    put("nrna", "ratio_labels", "nrna_labels")

    put("capture", "fraction", "capture_fraction")
    put("capture", "capture_fraction")
    put("capture", "probe_length")
    put("capture", "probe_density")
    put("capture", "off_target_weight", "capture_off_target_weight")
    put("capture", "binding_per_base", "capture_binding_per_base")
    put("capture", "gdna_split_penalty", "capture_gdna_split_penalty")
    put("capture", "min_overlap", "capture_min_overlap")

    return values


@dataclass(frozen=True)
class CaptureProbeDesignResult:
    """Summary of generated random-transcriptome capture probes."""

    path: Path
    n_transcripts: int
    n_eligible: int
    n_captured: int
    n_probes: int


def design_capture_probe_intervals(
    transcript_length: int,
    probe_length: int = 120,
    probe_density: float = 1.0,
) -> list[tuple[int, int]]:
    """Return non-overlapping, centered probe intervals on transcript coordinates."""
    if probe_length <= 0:
        raise ValueError("probe_length must be > 0")
    if not 0.0 <= probe_density <= 1.0:
        raise ValueError("probe_density must be between 0 and 1")
    if transcript_length < probe_length or probe_density <= 0:
        return []

    max_probes = int(transcript_length) // int(probe_length)
    n_probes = int(np.floor(max_probes * probe_density))
    n_probes = max(1, n_probes)
    total_probe_bases = n_probes * probe_length
    total_gap = int(transcript_length) - total_probe_bases
    gap = total_gap / (n_probes + 1)

    intervals: list[tuple[int, int]] = []
    previous_end = 0
    for i in range(n_probes):
        start = int(round(gap * (i + 1) + probe_length * i))
        start = max(start, previous_end)
        start = min(start, int(transcript_length) - probe_length)
        end = start + probe_length
        if start < previous_end:
            raise RuntimeError("probe tiling produced overlapping probes")
        intervals.append((start, end))
        previous_end = end
    return intervals


def write_random_capture_probes(
    transcripts: list[Transcript],
    path: Path,
    *,
    capture_fraction: float,
    probe_length: int = 120,
    probe_density: float = 1.0,
    seed: int = 42,
) -> CaptureProbeDesignResult:
    """Select captured transcripts and write transcript-coordinate probe TSV."""
    if not 0.0 <= capture_fraction <= 1.0:
        raise ValueError("capture_fraction must be between 0 and 1")
    if probe_length <= 0:
        raise ValueError("probe_length must be > 0")
    if not 0.0 <= probe_density <= 1.0:
        raise ValueError("probe_density must be between 0 and 1")

    eligible: list[tuple[int, Transcript, int]] = []
    for idx, transcript in enumerate(transcripts):
        length = int(transcript.length or transcript.compute_length())
        if length >= probe_length and transcript.t_id is not None:
            eligible.append((idx, transcript, length))

    if capture_fraction <= 0.0 or probe_density <= 0.0 or not eligible:
        path.write_text("transcript_id\tstart\tend\n")
        return CaptureProbeDesignResult(path, len(transcripts), len(eligible), 0, 0)

    n_captured = int(np.floor(len(eligible) * capture_fraction))
    n_captured = min(len(eligible), max(1, n_captured))
    rng = np.random.default_rng(seed)
    selected_positions = set(rng.choice(len(eligible), size=n_captured, replace=False))

    n_probes = 0
    with open(path, "w") as handle:
        handle.write("transcript_id\tstart\tend\n")
        for pos, (_, transcript, length) in enumerate(eligible):
            if pos not in selected_positions:
                continue
            for start, end in design_capture_probe_intervals(
                length, probe_length=probe_length, probe_density=probe_density,
            ):
                handle.write(f"{transcript.t_id}\t{start}\t{end}\n")
                n_probes += 1

    return CaptureProbeDesignResult(path, len(transcripts), len(eligible), n_captured, n_probes)


def main():
    parser = argparse.ArgumentParser(
        description="Full synthetic mini-genome simulation pipeline"
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="YAML suite config; explicit CLI flags override matching config values"
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic"),
        help="Base output directory"
    )
    parser.add_argument(
        "--genome-length", type=int, default=GENOME_LENGTH,
        help=f"Genome length in bp (default: {GENOME_LENGTH:,})"
    )
    parser.add_argument(
        "--n-genes", type=int, default=N_GENES,
        help=f"Number of genes (default: {N_GENES})"
    )
    parser.add_argument(
        "--min-isoforms", type=int, default=MIN_ISOFORMS,
        help=f"Minimum transcript isoforms per primary gene (default: {MIN_ISOFORMS})"
    )
    parser.add_argument(
        "--max-isoforms", type=int, default=MAX_ISOFORMS,
        help=f"Maximum transcript isoforms per primary gene (default: {MAX_ISOFORMS})"
    )
    parser.add_argument(
        "--target-transcripts", type=int, default=TARGET_TRANSCRIPTS,
        help=(
            f"Approximate target transcript count for isoform sampling "
            f"(default: {TARGET_TRANSCRIPTS}; use 0 for uniform min/max sampling)"
        )
    )
    parser.add_argument(
        "--antisense-fraction", type=float, default=ANTISENSE_OVERLAP_FRAC,
        help=(
            "Fraction of primary genes with an antisense overlapping partner "
            f"(default: {ANTISENSE_OVERLAP_FRAC})"
        )
    )
    parser.add_argument(
        "--seed", type=int, default=SEED,
        help=f"Random seed (default: {SEED})"
    )
    parser.add_argument(
        "--profile", choices=("smoke", "full"), default="smoke",
        help="Condition profile: smoke = 10 no-nRNA baseline conditions, full = 40 conditions"
    )
    parser.add_argument(
        "--n-rna", type=int, default=1_000_000,
        help="Number of mature RNA fragments per condition (default: 1,000,000)"
    )
    parser.add_argument("--frag-mean", type=float, default=250.0)
    parser.add_argument("--frag-std", type=float, default=50.0)
    parser.add_argument("--frag-min", type=int, default=80)
    parser.add_argument("--frag-max", type=int, default=800)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--error-rate", type=float, default=0.0)
    parser.add_argument("--abundance-seed", type=int, default=None)
    parser.add_argument("--abundance-min", type=float, default=1.0)
    parser.add_argument("--abundance-max", type=float, default=10000.0)
    parser.add_argument("--frac-expressed", type=float, default=0.5)
    parser.add_argument("--gdna-rates", type=str, default=None)
    parser.add_argument("--gdna-labels", type=str, default=None)
    parser.add_argument("--gdna-frag-mean", type=float, default=350.0)
    parser.add_argument("--gdna-frag-std", type=float, default=100.0)
    parser.add_argument("--gdna-frag-min", type=int, default=100)
    parser.add_argument("--gdna-frag-max", type=int, default=1000)
    parser.add_argument("--strand-specificities", type=str, default=None)
    parser.add_argument(
        "--capture-fraction", type=float, default=0.0,
        help="Fraction of eligible transcripts to target with generated probes (default: 0.0)"
    )
    parser.add_argument(
        "--probe-length", type=int, default=120,
        help="Generated capture probe length in transcript coordinates (default: 120)"
    )
    parser.add_argument(
        "--probe-density", type=float, default=1.0,
        help="Fraction of each captured transcript's non-overlapping probe tiling to emit (default: 1.0)"
    )
    parser.add_argument(
        "--capture-off-target-weight", type=float, default=1.0,
        help="Hybrid-capture baseline weight for every legal fragment start"
    )
    parser.add_argument(
        "--capture-binding-per-base", type=float, default=10.0,
        help="Hybrid-capture extra weight per overlapping probe base"
    )
    parser.add_argument(
        "--capture-gdna-split-penalty", type=float, default=0.2,
        help="Hybrid-capture multiplier for probes split by introns in gDNA/pre-mRNA"
    )
    parser.add_argument(
        "--capture-min-overlap", type=int, default=1,
        help="Minimum fragment/probe overlap required to add capture weight"
    )
    parser.add_argument(
        "--nrna-ratios", type=str, default=None,
        help="Comma-separated additive nRNA:mRNA fragment ratios (overrides --profile)"
    )
    parser.add_argument(
        "--nrna-labels", type=str, default=None,
        help="Comma-separated labels for --nrna-ratios (required when overriding ratios)"
    )
    parser.add_argument(
        "--conditions", nargs="*", default=None,
        help="Optional condition names to generate from the selected profile/grid"
    )
    parser.add_argument(
        "-j", "--n-workers", type=int, default=1,
        help="Worker processes for parallel read generation"
    )
    parser.add_argument(
        "--reference-only", action="store_true",
        help="Generate only the synthetic reference FASTA/GTF and exit"
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip conditions that already have output files"
    )
    args = parser.parse_args()

    if args.config is not None:
        _apply_config_defaults(args, parser, _load_suite_config(args.config))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    ref_dir = outdir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()
    if args.nrna_ratios is None:
        nrna_ratios, nrna_labels = PROFILE_NRNA[args.profile]
    else:
        if args.nrna_labels is None:
            raise ValueError("--nrna-labels is required when --nrna-ratios is set")
        nrna_ratios = _as_float_list(args.nrna_ratios)
        nrna_labels = _as_label_list(args.nrna_labels)
    if len(nrna_labels) != len(nrna_ratios):
        raise ValueError("--nrna-labels must have the same length as --nrna-ratios")
    if args.min_isoforms <= 0:
        raise ValueError("--min-isoforms must be > 0")
    if args.max_isoforms < args.min_isoforms:
        raise ValueError("--max-isoforms must be >= --min-isoforms")
    if not 0.0 <= args.antisense_fraction <= 1.0:
        raise ValueError("--antisense-fraction must be between 0 and 1")
    target_transcripts = args.target_transcripts
    if target_transcripts is not None and target_transcripts <= 0:
        target_transcripts = None
    if not 0.0 <= args.capture_fraction <= 1.0:
        raise ValueError("--capture-fraction must be between 0 and 1")
    if args.probe_length <= 0:
        raise ValueError("--probe-length must be > 0")
    if not 0.0 <= args.probe_density <= 1.0:
        raise ValueError("--probe-density must be between 0 and 1")
    if not 0.0 <= args.frac_expressed <= 1.0:
        raise ValueError("--frac-expressed must be between 0 and 1")
    gdna_rates = _as_float_list(args.gdna_rates) if args.gdna_rates is not None else GDNA_RATES
    gdna_labels = _as_label_list(args.gdna_labels) if args.gdna_labels is not None else GDNA_LABELS
    strand_specs = (
        _as_float_list(args.strand_specificities)
        if args.strand_specificities is not None else STRAND_SPECIFICITIES
    )
    if len(gdna_labels) != len(gdna_rates):
        raise ValueError("gDNA labels must have the same length as gDNA rates")
    selected_conditions = set(args.conditions or [])

    # ══════════════════════════════════════════════════════════════════════
    # Step 1: Generate synthetic genome
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 70)
    print("STEP 1: Generating synthetic genome")
    print("=" * 70)

    fasta_path = ref_dir / "genome.fa"
    gtf_path = ref_dir / "genes.gtf"

    if fasta_path.exists() and gtf_path.exists() and args.skip_existing:
        print(f"  Genome already exists at {ref_dir}, skipping generation")
    else:
        genes = generate_genes(
            args.genome_length,
            args.n_genes,
            args.seed,
            min_isoforms=args.min_isoforms,
            max_isoforms=args.max_isoforms,
            target_transcripts=target_transcripts,
            antisense_overlap_frac=args.antisense_fraction,
        )
        n_transcripts = sum(len(g.transcripts) for g in genes)
        n_antisense = sum(1 for g in genes if g.gene_name.endswith("as"))

        # Generate genome sequence
        genome_seq = generate_random_genome(args.genome_length, args.seed)

        # Inject splice sites
        all_transcripts = [tx for g in genes for tx in g.transcripts]
        seq_array = np.frombuffer(genome_seq.encode("ascii"), dtype=np.uint8).copy()
        inject_splice_sites(seq_array, all_transcripts)
        genome_seq = seq_array.tobytes().decode("ascii")

        # Write files
        fasta_path = write_fasta(genome_seq, REF_NAME, ref_dir)
        gtf_path = write_gtf(genes, REF_NAME, ref_dir)

        iso_counts = [len(g.transcripts) for g in genes]
        n_multi_exon = sum(1 for tx in all_transcripts if len(tx.exons) > 1)

        print(f"  Genome:       {fasta_path} ({args.genome_length:,} bp)")
        print(f"  Annotation:   {gtf_path}")
        print(f"  Genes:        {len(genes)} ({n_antisense} antisense overlapping)")
        print(f"  Transcripts:  {n_transcripts} ({n_multi_exon} multi-exon)")
        if iso_counts:
            print(
                f"  Isoforms/gene: {min(iso_counts)}-{max(iso_counts)} "
                f"(mean {np.mean(iso_counts):.1f})"
            )
        else:
            print("  Isoforms/gene: n/a")

    transcripts: list[Transcript] | None = None
    capture_result: CaptureProbeDesignResult | None = None
    capture_config = CaptureConfig()
    if args.capture_fraction > 0.0 and args.probe_density > 0.0:
        print(f"\n  Loading transcripts from {gtf_path} for capture probe design...")
        transcripts = load_transcripts(gtf_path)
        probe_path = ref_dir / "capture_probes.tsv"
        capture_result = write_random_capture_probes(
            transcripts,
            probe_path,
            capture_fraction=args.capture_fraction,
            probe_length=args.probe_length,
            probe_density=args.probe_density,
            seed=stable_seed(args.seed, "capture", args.capture_fraction, args.probe_length,
                             args.probe_density),
        )
        if capture_result.n_probes > 0:
            capture_config = CaptureConfig(
                probes=str(probe_path),
                probe_format="transcript",
                off_target_weight=args.capture_off_target_weight,
                binding_per_base=args.capture_binding_per_base,
                gdna_split_penalty=args.capture_gdna_split_penalty,
                min_overlap=args.capture_min_overlap,
            )
            print(
                f"  Capture probes: {probe_path} "
                f"({capture_result.n_probes:,} probes across "
                f"{capture_result.n_captured:,}/{capture_result.n_eligible:,} eligible transcripts)"
            )
        else:
            print("  Capture requested but no eligible probes were generated; capture disabled")

    if args.reference_only:
        print(f"\nReference generated at {ref_dir}")
        return 0

    # ══════════════════════════════════════════════════════════════════════
    # Step 2: Run simulation across conditions
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("STEP 2: Simulating reads across conditions")
    print("=" * 70)

    # Configuration
    n_rna = args.n_rna

    # Build config object for manifest
    cfg = WholeGenomeSimConfig(
        genome=str(fasta_path),
        gtf=str(gtf_path),
        outdir=str(outdir),
        transcript_filter="all",
        simulation=SimulationParams(
            n_rna_fragments=n_rna,
            sim_seed=args.seed,
            frag_mean=args.frag_mean,
            frag_std=args.frag_std,
            frag_min=args.frag_min,
            frag_max=args.frag_max,
            read_length=args.read_length,
            error_rate=args.error_rate,
            n_workers=args.n_workers,
        ),
        abundance=AbundanceConfig(
            mode="random",
            seed=(
                args.abundance_seed if args.abundance_seed is not None
                else stable_seed(args.seed, "abundance")
            ),
            min=args.abundance_min,
            max=args.abundance_max,
            frac_expressed=args.frac_expressed,
        ),
        gdna=GDNASimConfig(
            rates=gdna_rates,
            rate_labels=gdna_labels,
            frag_mean=args.gdna_frag_mean,
            frag_std=args.gdna_frag_std,
            frag_min=args.gdna_frag_min,
            frag_max=args.gdna_frag_max,
        ),
        nrna=NRNAConfig(
            mode="additive_ratio",
            ratios=nrna_ratios,
            ratio_labels=nrna_labels,
        ),
        capture=capture_config,
        strand_specificities=strand_specs,
        oracle_bam=True,
    )

    # Load transcripts and assign abundances
    if transcripts is None:
        print(f"\n  Loading transcripts from {gtf_path}...")
        transcripts = load_transcripts(gtf_path)
    assign_random_abundances(transcripts, cfg.abundance)

    n_expressed = sum(1 for t in transcripts if (t.abundance or 0) > 0)
    print(f"  Expressed: {n_expressed}/{len(transcripts)} transcripts")
    if cfg.capture.probes:
        print(f"  Hybrid capture: {cfg.capture.probes}")

    # Run conditions
    all_condition_names = [
        condition_dir_name(gdna_label, strand_spec, nrna_label)
        for nrna_label in nrna_labels
        for gdna_label in gdna_labels
        for strand_spec in strand_specs
    ]
    unknown_conditions = selected_conditions.difference(all_condition_names)
    if unknown_conditions:
        unknown = ", ".join(sorted(unknown_conditions))
        raise ValueError(f"Unknown condition(s) for selected profile/grid: {unknown}")
    total_conditions = (
        len(selected_conditions) if selected_conditions else len(all_condition_names)
    )
    conditions: list[dict] = []
    cond_num = 0

    for nrna_ratio, nrna_label in zip(nrna_ratios, nrna_labels):
        cond_transcripts = copy.deepcopy(transcripts)
        apply_nrna_ratio(cond_transcripts, nrna_ratio)
        truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        truth_path = outdir / truth_name
        write_truth_abundances(cond_transcripts, truth_path)
        print(f"  Truth abundances ({nrna_label}): {truth_path}")

        for i, gdna_rate in enumerate(gdna_rates):
            gdna_label = gdna_labels[i]
            for strand_spec in strand_specs:
                n_mrna = n_rna
                n_nrna = round(n_mrna * nrna_ratio)
                n_rna_total = n_mrna + n_nrna
                n_gdna = round(gdna_rate * n_mrna)

                cond_name = condition_dir_name(gdna_label, strand_spec, nrna_label)
                if selected_conditions and cond_name not in selected_conditions:
                    continue
                cond_num += 1
                condition_seed = stable_seed(args.seed, cond_name)
                cond_dir = outdir / cond_name

                print(
                    f"\n  [{cond_num}/{total_conditions}] {cond_name}: "
                    f"mRNA={n_mrna:,} nRNA={n_nrna:,} gDNA={n_gdna:,} "
                    f"SS={strand_spec:.2f}",
                    flush=True,
                )

                cond_entry = {
                    "name": cond_name,
                    "gdna_label": gdna_label,
                    "gdna_rate": gdna_rate,
                    "strand_specificity": strand_spec,
                    "nrna_label": nrna_label,
                    "nrna_mode": "additive_ratio",
                    "nrna_ratio": nrna_ratio,
                    "n_mrna": n_mrna,
                    "n_nrna": n_nrna,
                    "n_rna": n_rna_total,
                    "n_gdna": n_gdna,
                    "n_total": n_rna_total + n_gdna,
                    "seed": condition_seed,
                    "truth_abundances": truth_name,
                    "fastq_r1": f"{cond_name}/sim_R1.fq.gz",
                    "fastq_r2": f"{cond_name}/sim_R2.fq.gz",
                }

                # Check if already done
                if (cond_dir / "sim_R1.fq.gz").exists() and args.skip_existing:
                    print("    Output exists, skipping", flush=True)
                    cond_entry["oracle_bam"] = f"{cond_name}/sim_oracle.bam"
                else:
                    print("    Simulating...", end="", flush=True)
                    t_cond = time.monotonic()
                    cond_sim = copy.deepcopy(cfg.simulation)
                    cond_sim.sim_seed = condition_seed
                    simulator = WholeGenomeSimulator(
                        fasta_path, cond_transcripts, cond_sim, cfg.gdna,
                        strand_specificity=strand_spec,
                        capture_config=cfg.capture,
                    )
                    _, _, bam_path = simulator.simulate_and_write(
                        cond_dir, n_rna_total, n_gdna,
                        oracle_bam=True, prefix="sim",
                        n_mrna=n_mrna,
                        n_nrna=n_nrna,
                        n_workers=cond_sim.n_workers,
                    )
                    simulator.close()
                    cond_entry["oracle_bam"] = (
                        f"{cond_name}/sim_oracle.bam" if bam_path else None
                    )
                    elapsed_cond = time.monotonic() - t_cond
                    print(f" done ({elapsed_cond:.1f}s)", flush=True)

                conditions.append(cond_entry)

    # Write manifest
    write_manifest(outdir, cfg, conditions)

    elapsed_total = time.monotonic() - t0
    print(f"\n{'=' * 70}")
    print("SIMULATION COMPLETE")
    print(f"  Output:     {outdir}")
    print(f"  Conditions: {total_conditions}")
    print(f"  Elapsed:    {elapsed_total:.1f}s")
    print(f"{'=' * 70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
