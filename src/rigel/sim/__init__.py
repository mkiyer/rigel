"""
rigel.sim — Simulation framework for synthetic RNA-seq test scenarios.

Modules
-------
genome
    ``MutableGenome`` — random DNA generation with positional editing.
annotation
    ``GeneBuilder`` — gene/transcript annotation with splice-motif injection.
reads
    ``ReadSimulator``, ``SimConfig`` — paired-end FASTQ generation.
scenario
    ``Scenario``, ``ScenarioResult`` — end-to-end orchestration
    (genome → GTF → FASTQ → BAM → TranscriptIndex).
benchmark
    ``BenchmarkResult``, ``TranscriptAccuracy``, ``run_benchmark`` —
    accuracy benchmarking comparing observed vs expected counts.

Quick Start
-----------
>>> from rigel.sim import Scenario
>>> with Scenario("test1", genome_length=5000, seed=42) as sc:
...     sc.add_gene("g1", "+", [
...         {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
...     ])
...     result = sc.build(n_fragments=500)
"""

from .annotation import GeneBuilder
from .benchmark import BenchmarkResult, TranscriptAccuracy, run_benchmark
from .oracle_bam import OracleBamSimulator
from .genome import MutableGenome, reverse_complement
from .manifest import (
    condition_dir_name,
    condition_manifest_map,
    load_manifest,
    write_manifest,
)
from .reads import GDNAConfig, ReadSimulator, SimConfig
from .scenario import Scenario, ScenarioResult
from .truth import Origin, parse_origin

__all__ = [
    "BenchmarkResult",
    "OracleBamSimulator",
    "GDNAConfig",
    "GeneBuilder",
    "MutableGenome",
    "Origin",
    "ReadSimulator",
    "ScenarioResult",
    "Scenario",
    "SimConfig",
    "TranscriptAccuracy",
    "condition_dir_name",
    "condition_manifest_map",
    "load_manifest",
    "parse_origin",
    "reverse_complement",
    "run_benchmark",
    "write_manifest",
]
