"""Simulation of synthetic RNA-seq libraries with per-fragment truth.

``Scenario`` builds a small random genome, its annotation, the index and the reads end to end; it is what
``rigel sim`` runs. Its engine, ``wgs_engine.WholeGenomeSimulator``, also drives the whole-genome condition grids the
lab's panels are simulated with (``scripts/sim/``), and every read name carries its fragment's origin
(``read_name``).

Quick Start
-----------
>>> from rigel.sim import Scenario
>>> with Scenario("test1", genome_length=5000, seed=42) as sc:
...     sc.add_gene("g1", "+", [
...         {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
...     ])
...     result = sc.build(n_fragments=500)
"""

from .reads import GDNAConfig, ReadSimConfig
from .scenario import Scenario, ScenarioResult

__all__ = ["GDNAConfig", "ReadSimConfig", "Scenario", "ScenarioResult"]
