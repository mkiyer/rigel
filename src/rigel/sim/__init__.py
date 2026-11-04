"""Simulation of synthetic RNA-seq libraries with per-fragment truth.

``Scenario`` builds a small random genome, its annotation, the index and the reads end to end; ``whole_genome`` runs
a YAML-configured condition grid on a real reference. Both drive one engine, ``wgs_engine.WholeGenomeSimulator``,
and every read name carries its fragment's origin (``read_name``).

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
from .capture import CaptureConfig, CaptureSampler, ProbeInterval
from .genome import MutableGenome, reverse_complement
from .manifest import condition_dir_name, write_manifest
from .reads import GDNAConfig, ReadSimConfig
from .scenario import Scenario, ScenarioResult
from .read_name import Origin, parse_origin

__all__ = [
    "BenchmarkResult",
    "CaptureConfig",
    "CaptureSampler",
    "GDNAConfig",
    "GeneBuilder",
    "MutableGenome",
    "Origin",
    "ScenarioResult",
    "Scenario",
    "ReadSimConfig",
    "TranscriptAccuracy",
    "ProbeInterval",
    "condition_dir_name",
    "parse_origin",
    "reverse_complement",
    "run_benchmark",
    "write_manifest",
]
