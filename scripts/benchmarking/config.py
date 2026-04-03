"""Benchmark configuration: YAML loading and dataclass definitions."""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

DEFAULT_BENCHMARK_DIR = (
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks"
)

# Default BAM location relative to each condition directory
DEFAULT_BAM_RELPATH = "rigel_star/annotated.bam"

# rigel quant CLI flags that named configs can set.
# Maps config YAML key → rigel CLI flag.
RIGEL_CLI_FLAGS: dict[str, str] = {
    "seed": "--seed",
    "prior_pseudocount": "--prior-pseudocount",
    "em_iterations": "--em-iterations",
    "em_mode": "--em-mode",
    "em_convergence_delta": "--em-convergence-delta",
    "assignment_mode": "--assignment-mode",
    "assignment_min_posterior": "--assignment-min-posterior",
    "overhang_alpha": "--overhang-alpha",
    "mismatch_alpha": "--mismatch-alpha",
    "gdna_splice_penalty_unannot": "--gdna-splice-penalty-unannot",
    "pruning_min_posterior": "--pruning-min-posterior",
    "threads": "--threads",
    "buffer_size": "--buffer-size",
    "include_multimap": "--include-multimap",
    "no_include_multimap": "--no-include-multimap",
    "keep_duplicates": "--keep-duplicates",
    "no_keep_duplicates": "--no-keep-duplicates",
    "sj_strand_tag": "--sj-strand-tag",
    "annotated_bam": "--annotated-bam",
    "tsv": "--tsv",
    "tmpdir": "--tmpdir",
}


@dataclass
class RigelNamedConfig:
    """A named rigel configuration mapping to ``rigel quant`` CLI params."""

    name: str
    params: dict[str, object] = field(default_factory=dict)
    bam_relpath: str | None = None  # override per-config BAM path

    def to_cli_args(self) -> list[str]:
        """Convert params to ``rigel quant`` CLI argument list."""
        args: list[str] = []
        for key, value in self.params.items():
            flag = RIGEL_CLI_FLAGS.get(key)
            if flag is None:
                logger.warning("Unknown rigel param %r — skipping", key)
                continue
            if isinstance(value, bool):
                if value:
                    args.append(flag)
                # For bool False with --no- prefix, skip
                continue
            if isinstance(value, (list, tuple)):
                args.append(flag)
                args.extend(str(v) for v in value)
            else:
                args.extend([flag, str(value)])
        return args


@dataclass
class AnalysisConfig:
    """Options for the analysis step."""

    log2_pseudocount: float = 1.0
    parse_fastq_truth: bool = False
    output_dir: str = "results/benchmark_report"


@dataclass
class BenchmarkConfig:
    """Top-level benchmark configuration."""

    benchmark_dir: Path
    rigel_index: Path
    bam_relpath: str = DEFAULT_BAM_RELPATH
    rigel_configs: dict[str, RigelNamedConfig] = field(default_factory=dict)
    conditions: list[str] | None = None  # None = discover all
    threads: int = 4
    seed: int = 42
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)

    # External tool indexes (optional)
    salmon_index: str | None = None
    kallisto_index: str | None = None

    @property
    def runs_dir(self) -> Path:
        return self.benchmark_dir / "runs"

    def condition_dir(self, condition: str) -> Path:
        return self.runs_dir / condition

    def bam_path(self, condition: str, config: RigelNamedConfig | None = None) -> Path:
        """Resolve BAM input path for a condition + config."""
        relpath = DEFAULT_BAM_RELPATH
        if config and config.bam_relpath:
            relpath = config.bam_relpath
        elif self.bam_relpath:
            relpath = self.bam_relpath
        return self.condition_dir(condition) / relpath

    def output_dir(self, condition: str, config_name: str) -> Path:
        """Output directory for a rigel run: runs/<cond>/rigel/<config>/."""
        return self.condition_dir(condition) / "rigel" / config_name

    def discover_conditions(self) -> list[str]:
        """Discover available conditions from the runs directory."""
        if not self.runs_dir.exists():
            raise FileNotFoundError(f"Runs directory not found: {self.runs_dir}")
        conditions = sorted(
            d.name
            for d in self.runs_dir.iterdir()
            if d.is_dir() and d.name.startswith("gdna_")
        )
        return conditions

    def get_conditions(self) -> list[str]:
        """Return explicit condition list or discover all."""
        if self.conditions:
            return self.conditions
        return self.discover_conditions()

    def load_manifest(self) -> dict:
        """Load manifest.json from the benchmark directory."""
        path = self.benchmark_dir / "manifest.json"
        if not path.exists():
            raise FileNotFoundError(f"Manifest not found: {path}")
        with open(path) as f:
            return json.load(f)


def load_config(path: str | Path) -> BenchmarkConfig:
    """Load a BenchmarkConfig from a YAML file."""
    path = Path(path)
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    # Benchmark directory
    benchmark_dir = Path(raw.get("benchmark_dir", DEFAULT_BENCHMARK_DIR))
    rigel_index = Path(raw.get("rigel_index", str(benchmark_dir / "rigel_index")))
    bam_relpath = raw.get("bam_relpath", DEFAULT_BAM_RELPATH)
    threads = raw.get("threads", 4)
    seed = raw.get("seed", 42)
    conditions = raw.get("conditions")

    # Parse named rigel configs
    rigel_configs: dict[str, RigelNamedConfig] = {}
    for name, cfg_dict in raw.get("rigel_configs", {}).items():
        if not isinstance(cfg_dict, dict):
            cfg_dict = {}
        bam_override = cfg_dict.pop("bam_relpath", None)
        rigel_configs[name] = RigelNamedConfig(
            name=name,
            params=cfg_dict,
            bam_relpath=bam_override,
        )

    # Analysis config
    analysis_raw = raw.get("analysis", {})
    analysis = AnalysisConfig(
        log2_pseudocount=analysis_raw.get("log2_pseudocount", 1.0),
        parse_fastq_truth=analysis_raw.get("parse_fastq_truth", False),
        output_dir=analysis_raw.get("output_dir", "results/benchmark_report"),
    )

    # External tool indexes
    salmon_index = raw.get("salmon_index")
    kallisto_index = raw.get("kallisto_index")

    return BenchmarkConfig(
        benchmark_dir=benchmark_dir,
        rigel_index=rigel_index,
        bam_relpath=bam_relpath,
        rigel_configs=rigel_configs,
        conditions=conditions,
        threads=threads,
        seed=seed,
        analysis=analysis,
        salmon_index=salmon_index,
        kallisto_index=kallisto_index,
    )
