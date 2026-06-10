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
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from rigel.sim.capture import CaptureConfig, CaptureScenario
from rigel.sim.capture.design import CaptureProbeDesignResult, write_random_capture_probes
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
from rigel.sim.orchestrator import run_condition_grid, stable_seed
from rigel.sim.orchestrator import (  # noqa: F401  re-exported for tests/back-compat
    capture_paired_condition_seed,
)
from rigel.sim.whole_genome import (
    AbundanceConfig,
    GDNASimConfig,
    NRNAConfig,
    WholeGenomeSimConfig,
    SimulationParams,
    assign_random_abundances,
    load_transcripts,
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
                value = section_map[source_key]
                values[dest_key] = Path(value) if dest_key == "outdir" else value
                return
        if source_key in raw:
            value = raw[source_key]
            values[dest_key] = Path(value) if dest_key == "outdir" else value

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
    put("gdna", "strand_overdispersions", "gdna_strand_overdispersions")
    put("gdna", "strand_overdispersion_labels", "gdna_strand_overdispersion_labels")

    put(None, "strand_specificities")
    put("nrna", "ratios", "nrna_ratios")
    put("nrna", "ratio_labels", "nrna_labels")
    put("nrna", "mode", "nrna_mode")
    put("nrna", "ratio_ranges", "nrna_ratio_ranges")
    put("nrna", "eligible_fraction", "nrna_eligible_fraction")
    put("nrna", "seed", "nrna_seed")

    put("capture", "fraction", "capture_fraction")
    put("capture", "capture_fraction")
    put("capture", "probes", "capture_probes")
    put("capture", "probe_format", "capture_probe_format")
    put("capture", "format", "capture_probe_format")
    put("capture", "probe_length")
    put("capture", "probe_density")
    put("capture", "off_target_weight", "capture_off_target_weight")
    put("capture", "binding_per_base", "capture_binding_per_base")
    put("capture", "gdna_split_penalty", "capture_gdna_split_penalty")
    put("capture", "min_overlap", "capture_min_overlap")
    put("capture", "configs", "capture_configs")
    put("capture", "scenarios", "capture_configs")

    return values


@dataclass(frozen=True)
class SuiteCaptureSpec:
    """Capture setting for either provided panels or generated probes."""

    label: str
    fraction: float
    probe_length: int
    probe_density: float
    off_target_weight: float
    binding_per_base: float
    gdna_split_penalty: float
    min_overlap: int
    probes: str | None = None
    probe_format: str = "auto"

    @property
    def enabled(self) -> bool:
        return self.uses_provided_probes or self.generates_probes

    @property
    def uses_provided_probes(self) -> bool:
        return bool(self.probes)

    @property
    def generates_probes(self) -> bool:
        return not self.uses_provided_probes and self.fraction > 0.0 and self.probe_density > 0.0


CaptureProbeGroupKey = tuple[float, int, float]


def _probe_float_key(value: float) -> float:
    return round(float(value), 12)


def capture_probe_group_key(spec: SuiteCaptureSpec) -> CaptureProbeGroupKey:
    """Return the probe-geometry key shared by compatible capture scenarios."""
    return (
        _probe_float_key(spec.fraction),
        int(spec.probe_length),
        _probe_float_key(spec.probe_density),
    )


def _float_slug(value: float) -> str:
    return f"{value:.12g}".replace("-", "m").replace(".", "p")


def _capture_probe_group_slug(key: CaptureProbeGroupKey) -> str:
    fraction, probe_length, probe_density = key
    return f"f{_float_slug(fraction)}_len{probe_length}_den{_float_slug(probe_density)}"


def _capture_probe_group_path(
    ref_dir: Path,
    key: CaptureProbeGroupKey,
    *,
    include_capture_in_names: bool,
) -> Path:
    if not include_capture_in_names:
        return ref_dir / "capture_probes.tsv"
    return ref_dir / f"capture_probes_{_capture_probe_group_slug(key)}.tsv"


def _provided_probe_kind(probes: str, probe_format: str) -> str:
    """Best-effort manifest classification for a provided capture panel."""
    fmt = probe_format.lower().replace("-", "_")
    if fmt in {"bed", "bed12"}:
        return "bed"
    if fmt in {"transcript", "transcript_tsv", "tsv"}:
        return "transcript"
    suffix = Path(probes).suffix.lower()
    if suffix in {".bed", ".bed12"}:
        return "bed"
    if suffix in {".tsv", ".txt"}:
        return "transcript"
    return "unknown"


def _capture_label(raw: dict, index: int) -> str:
    value = raw.get("label", raw.get("name", f"c{index}"))
    label = "on" if value is True else "off" if value is False else str(value)
    if not label or "/" in label:
        raise ValueError("capture config labels must be non-empty path-safe strings")
    return label


def _suite_capture_specs(args: argparse.Namespace) -> tuple[list[SuiteCaptureSpec], bool]:
    """Return generated-capture specs and whether names should include labels."""
    raw_configs = getattr(args, "capture_configs", None)

    def from_mapping(raw: dict, index: int) -> SuiteCaptureSpec:
        if not isinstance(raw, dict):
            raise ValueError("each capture.configs entry must be a mapping")
        label = _capture_label(raw, index)
        enabled = bool(raw.get("enabled", True))
        fraction = float(raw.get("fraction", raw.get("capture_fraction", args.capture_fraction)))
        probes = raw.get("probes", args.capture_probes)
        probe_format = str(
            raw.get("probe_format", raw.get("format", args.capture_probe_format))
        )
        if not enabled:
            fraction = 0.0
            probes = None
        spec = SuiteCaptureSpec(
            label=label,
            fraction=fraction,
            probe_length=int(raw.get("probe_length", args.probe_length)),
            probe_density=float(raw.get("probe_density", args.probe_density)),
            off_target_weight=float(
                raw.get("off_target_weight", args.capture_off_target_weight)
            ),
            binding_per_base=float(
                raw.get("binding_per_base", args.capture_binding_per_base)
            ),
            gdna_split_penalty=float(
                raw.get("gdna_split_penalty", args.capture_gdna_split_penalty)
            ),
            min_overlap=int(raw.get("min_overlap", args.capture_min_overlap)),
            probes=str(probes) if probes else None,
            probe_format=probe_format,
        )
        if not 0.0 <= spec.fraction <= 1.0:
            raise ValueError(f"capture config '{label}' fraction must be between 0 and 1")
        if spec.probe_length <= 0:
            raise ValueError(f"capture config '{label}' probe_length must be > 0")
        if not 0.0 <= spec.probe_density <= 1.0:
            raise ValueError(f"capture config '{label}' probe_density must be between 0 and 1")
        return spec

    if raw_configs is None:
        label = (
            "on" if args.capture_probes or (args.capture_fraction > 0.0 and args.probe_density > 0.0)
            else "off"
        )
        return [from_mapping({"label": label}, 1)], False
    if not isinstance(raw_configs, list) or not raw_configs:
        raise ValueError("capture.configs must be a non-empty list")

    specs = [from_mapping(raw, index) for index, raw in enumerate(raw_configs, start=1)]
    labels = [spec.label for spec in specs]
    if len(set(labels)) != len(labels):
        raise ValueError("capture config labels must be unique")
    return specs, True


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
    parser.add_argument("--gdna-strand-overdispersions", type=str, default=None,
                        help="Comma-separated gDNA strand overdispersion sweep values in [0, 1).")
    parser.add_argument("--gdna-strand-overdispersion-labels", type=str, default=None,
                        help="Comma-separated labels for --gdna-strand-overdispersions.")
    parser.add_argument("--strand-specificities", type=str, default=None)
    parser.add_argument(
        "--capture-fraction", type=float, default=0.0,
        help="Fraction of eligible genes to target with generated probes (default: 0.0)"
    )
    parser.add_argument(
        "--capture-probes", type=str, default=None,
        help="Provided capture probe panel path; skips generated probe design when set"
    )
    parser.add_argument(
        "--capture-probe-format", type=str, default="auto",
        help="Provided capture probe format: auto, transcript, or bed12 (default: auto)"
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
        help="Comma-separated labels for --nrna-ratios / --nrna-ratio-ranges"
    )
    parser.add_argument(
        "--nrna-mode", type=str, default="additive_ratio",
        choices=("additive_ratio", "random_fraction"),
        help="Nascent sweep mode: additive_ratio (one global ratio for all expressed multi-exon "
             "transcripts) or random_fraction (a random per-transcript ratio on a fraction of them)"
    )
    parser.add_argument(
        "--nrna-ratio-ranges", type=str, default=None,
        help="random_fraction: semicolon-separated lo,hi ranges, e.g. '0,0;0.1,1.0' "
             "(per-transcript ratio drawn ~Uniform(lo,hi))"
    )
    parser.add_argument(
        "--nrna-eligible-fraction", type=float, default=1.0,
        help="random_fraction: fraction of expressed multi-exon transcripts that carry nascent RNA"
    )
    parser.add_argument(
        "--nrna-seed", type=int, default=42,
        help="random_fraction: RNG seed for the per-transcript nascent assignment"
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
    nrna_mode = getattr(args, "nrna_mode", None) or "additive_ratio"
    if nrna_mode == "random_fraction":
        if args.nrna_ratio_ranges is None or args.nrna_labels is None:
            raise ValueError(
                "nrna.ratio_ranges and nrna.ratio_labels are required for nrna.mode=random_fraction"
            )
        raw = args.nrna_ratio_ranges
        if isinstance(raw, str):
            raw = [pair.split(",") for pair in raw.split(";")]
        nrna_ratio_ranges = [(float(lo), float(hi)) for lo, hi in raw]
        nrna_labels = _as_label_list(args.nrna_labels)
        nrna_ratios = None
        if len(nrna_labels) != len(nrna_ratio_ranges):
            raise ValueError("nrna.ratio_labels must match nrna.ratio_ranges in length")
    else:
        nrna_ratio_ranges = None
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

    # gDNA strand-overdispersion sweep axis (one condition per value). Default: a single value
    # of 0 (no overdispersion) — keeps condition names unchanged for suites that do not sweep it.
    gdna_ods = (
        _as_float_list(args.gdna_strand_overdispersions)
        if args.gdna_strand_overdispersions is not None
        else [0.0]
    )
    gdna_od_labels = (
        _as_label_list(args.gdna_strand_overdispersion_labels)
        if args.gdna_strand_overdispersion_labels is not None
        else None
    )
    gdna_od_pairs: list[tuple[str, float]] = []
    for i, od in enumerate(gdna_ods):
        if not (0.0 <= od < 1.0):
            raise ValueError(f"gDNA strand_overdispersion must be in [0, 1); got {od}")
        od_label = (
            gdna_od_labels[i]
            if gdna_od_labels is not None and i < len(gdna_od_labels)
            else ("od00" if od == 0.0 else f"od{od:g}".replace(".", "p"))
        )
        gdna_od_pairs.append((od_label, od))

    selected_conditions = set(args.conditions or [])
    capture_specs, include_capture_in_names = _suite_capture_specs(args)
    abundance_config = AbundanceConfig(
        mode="random",
        seed=(
            args.abundance_seed if args.abundance_seed is not None
            else stable_seed(args.seed, "abundance")
        ),
        min=args.abundance_min,
        max=args.abundance_max,
        frac_expressed=args.frac_expressed,
    )

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
    capture_scenarios: list[CaptureScenario] = []
    abundances_assigned = False
    if any(spec.generates_probes for spec in capture_specs):
        print(f"\n  Loading transcripts from {gtf_path} for capture probe design...")
        transcripts = load_transcripts(gtf_path)
        assign_random_abundances(transcripts, abundance_config)
        abundances_assigned = True

    capture_probe_results: dict[CaptureProbeGroupKey, CaptureProbeDesignResult] = {}
    capture_probe_groups: dict[CaptureProbeGroupKey, list[SuiteCaptureSpec]] = {}
    for spec in capture_specs:
        if spec.generates_probes:
            capture_probe_groups.setdefault(capture_probe_group_key(spec), []).append(spec)

    for group_key, group_specs in capture_probe_groups.items():
        if transcripts is None:
            continue
        probe_path = _capture_probe_group_path(
            ref_dir,
            group_key,
            include_capture_in_names=include_capture_in_names,
        )
        bed_path = probe_path.with_suffix(".bed")
        capture_probe_results[group_key] = write_random_capture_probes(
            transcripts,
            probe_path,
            capture_fraction=group_key[0],
            probe_length=group_key[1],
            probe_density=group_key[2],
            seed=stable_seed(args.seed, "capture", *group_key),
            bed_path=bed_path,
        )
        labels = ", ".join(spec.label for spec in group_specs)
        result = capture_probe_results[group_key]
        if result.n_probes > 0:
            print(
                f"  Capture probe set {_capture_probe_group_slug(group_key)}: "
                f"{result.path} + {result.bed_path} "
                f"({result.n_probes:,} probes across "
                f"{result.n_captured_genes:,}/{result.n_eligible_genes:,} eligible genes, "
                f"{result.n_captured:,}/{result.n_eligible:,} eligible transcripts; "
                f"scenarios: {labels})"
            )

    capture_probe_tsv_by_label: dict[str, str] = {}
    capture_probe_bed_by_label: dict[str, str] = {}
    capture_probe_panel_by_label: dict[str, str] = {}
    capture_probe_source_by_label: dict[str, str] = {}
    for spec in capture_specs:
        capture_config = CaptureConfig()
        if spec.uses_provided_probes and spec.probes is not None:
            capture_config = CaptureConfig(
                probes=spec.probes,
                probe_format=spec.probe_format,
                off_target_weight=spec.off_target_weight,
                binding_per_base=spec.binding_per_base,
                gdna_split_penalty=spec.gdna_split_penalty,
                min_overlap=spec.min_overlap,
            )
            capture_probe_panel_by_label[spec.label] = spec.probes
            capture_probe_source_by_label[spec.label] = "provided"
            probe_kind = _provided_probe_kind(spec.probes, spec.probe_format)
            if probe_kind == "bed":
                capture_probe_bed_by_label[spec.label] = spec.probes
            elif probe_kind == "transcript":
                capture_probe_tsv_by_label[spec.label] = spec.probes
            print(
                f"  Capture {spec.label}: using provided probe panel "
                f"{spec.probes} ({spec.probe_format})"
            )
        elif spec.generates_probes and transcripts is not None:
            capture_result = capture_probe_results.get(capture_probe_group_key(spec))
            if capture_result is not None and capture_result.n_probes > 0:
                capture_probe_tsv_by_label[spec.label] = str(capture_result.path)
                if capture_result.bed_path is not None:
                    capture_probe_bed_by_label[spec.label] = str(capture_result.bed_path)
                    capture_probe_panel_by_label[spec.label] = str(capture_result.bed_path)
                else:
                    capture_probe_panel_by_label[spec.label] = str(capture_result.path)
                capture_probe_source_by_label[spec.label] = "generated"
                capture_config = CaptureConfig(
                    probes=str(capture_result.bed_path),
                    probe_format="bed12",
                    off_target_weight=spec.off_target_weight,
                    binding_per_base=spec.binding_per_base,
                    gdna_split_penalty=spec.gdna_split_penalty,
                    min_overlap=spec.min_overlap,
                )
                print(
                    f"  Capture {spec.label}: using {capture_result.bed_path}"
                )
            else:
                print(
                    f"  Capture {spec.label} requested but no eligible probes were "
                    "generated; disabled"
                )
        else:
            print(
                f"  Capture {spec.label}: disabled"
                if include_capture_in_names else "  Capture: disabled"
            )
        capture_scenarios.append(CaptureScenario(label=spec.label, config=capture_config))

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
        abundance=abundance_config,
        gdna=GDNASimConfig(
            rates=gdna_rates,
            rate_labels=gdna_labels,
            frag_mean=args.gdna_frag_mean,
            frag_std=args.gdna_frag_std,
            frag_min=args.gdna_frag_min,
            frag_max=args.gdna_frag_max,
        ),
        nrna=NRNAConfig(
            mode=nrna_mode,
            ratios=nrna_ratios if nrna_ratios is not None else [0.0],
            ratio_ranges=nrna_ratio_ranges,
            ratio_labels=nrna_labels,
            eligible_fraction=float(getattr(args, "nrna_eligible_fraction", 1.0) or 1.0),
            seed=int(getattr(args, "nrna_seed", 42) or 42),
        ),
        capture=(
            capture_scenarios[0].config
            if capture_scenarios else CaptureConfig()
        ),
        capture_configs=(
            capture_scenarios if include_capture_in_names else []
        ),
        strand_specificities=strand_specs,
        oracle_bam=True,
    )

    # Load transcripts and assign abundances
    if transcripts is None:
        print(f"\n  Loading transcripts from {gtf_path}...")
        transcripts = load_transcripts(gtf_path)
    if not abundances_assigned:
        assign_random_abundances(transcripts, cfg.abundance)
        abundances_assigned = True

    n_expressed = sum(1 for t in transcripts if (t.abundance or 0) > 0)
    print(f"  Expressed: {n_expressed}/{len(transcripts)} transcripts")
    if include_capture_in_names:
        labels = ", ".join(
            f"{scenario.label}={'on' if scenario.config.probes else 'off'}"
            for scenario in capture_scenarios
        )
        print(f"  Capture configs: {labels}")
    elif cfg.capture.probes:
        print(f"  Hybrid capture: {cfg.capture.probes}")

    # Run conditions
    all_condition_names = [
        condition_dir_name(
            gdna_label,
            strand_spec,
            nrna_label,
            capture_scenario.label if include_capture_in_names else None,
            gdna_strand_overdispersion=gdna_od,
        )
        for nrna_label in nrna_labels
        for gdna_label in gdna_labels
        for _od_label, gdna_od in gdna_od_pairs
        for strand_spec in strand_specs
        for capture_scenario in capture_scenarios
    ]
    unknown_conditions = selected_conditions.difference(all_condition_names)
    if unknown_conditions:
        unknown = ", ".join(sorted(unknown_conditions))
        raise ValueError(f"Unknown condition(s) for selected profile/grid: {unknown}")
    total_conditions = (
        len(selected_conditions) if selected_conditions else len(all_condition_names)
    )
    capture_meta_by_label = {
        scenario.label: {
            "source": capture_probe_source_by_label.get(scenario.label),
            "panel": capture_probe_panel_by_label.get(scenario.label),
            "tsv": capture_probe_tsv_by_label.get(scenario.label),
            "bed": capture_probe_bed_by_label.get(scenario.label),
        }
        for scenario in capture_scenarios
    }
    base_abundances = [(t.abundance or 0.0, t.nrna_abundance) for t in transcripts]
    if nrna_mode == "random_fraction":
        nrna_pairs = [
            (label, "random_fraction", ratio_range, idx)
            for idx, (ratio_range, label) in enumerate(zip(nrna_ratio_ranges, nrna_labels))
        ]
    else:
        nrna_pairs = [
            (label, "additive_ratio", ratio, idx)
            for idx, (ratio, label) in enumerate(zip(nrna_ratios, nrna_labels))
        ]
    conditions = run_condition_grid(
        outdir=outdir,
        genome_path=fasta_path,
        transcripts=transcripts,
        base_abundances=base_abundances,
        sim=cfg.simulation,
        gdna=cfg.gdna,
        nrna=cfg.nrna,
        gdna_pairs=list(zip(gdna_labels, gdna_rates)),
        gdna_od_pairs=gdna_od_pairs,
        strand_specificities=strand_specs,
        nrna_pairs=nrna_pairs,
        capture_scenarios=capture_scenarios,
        include_capture_in_names=include_capture_in_names,
        base_seed=args.seed,
        oracle_bam=True,
        skip_existing=args.skip_existing,
        selected_conditions=selected_conditions or None,
        capture_meta_by_label=capture_meta_by_label,
    )

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
