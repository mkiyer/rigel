"""The instruments' window onto a calibration solve: two records the solve fills when handed one.

Production never passes a record, and nothing the solve computes reads one. ``calibrate(diagnostics=...)`` fills a
:class:`CalibrationDiagnostics` with the solve's intermediates; ``sweep.solve_chain(diagnostics=...)`` fills a
:class:`SweepDiagnostics` with what one sweep's kernel saw, and ``calibrate`` hands every sweep the same record, so it
ends holding the last.

Gate: ``tests/calibration/test_sweep_backbone.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .gdna_density import BackgroundLaw
    from .intron_departure import DepartureFit
    from .landscape import DensityLandscape
    from .lock_departure import LockFit
    from .region_arrays import RegionArrays
    from .region_chain import RegionChain
    from .region_geometry import RegionBelief, RegionGeometry, RegionStatics
    from .simplex_logodds import CubeRows
    from .splice_graph import SpliceJunctionGeometry
    from .substrate import CalibrationSubstrate
    from .sweep import AssertionCounts

__all__ = ["CalibrationDiagnostics", "SweepDiagnostics"]


@dataclass(slots=True)
class SweepDiagnostics:
    """What one sweep's kernel saw, filled by ``sweep.solve_chain``.

    Per slot (``(n,)`` unless stated): the message-free self-solve's ``fg_loc`` and ``tau_lam`` (the strand term's
    precision at ``fg_loc``); the final solve's ``f_g``; ``solvable``; the observations ``count`` ``(n, 2)``,
    ``free_pos``, ``free_neg``; and the gDNA support's ``mass_global``. The message layer's delivery: ``lam_rows``
    ``(n, K)`` (zero rows where a block delivered nothing, ``None`` when no block did), ``cube_rows`` (a ``CubeRows``
    table keyed to the chain), and the received tables ``from_left`` / ``from_right``, each a dict of the kernel's
    arrays over the chain (``has_neighbour``, ``has_composition``, ``composition``, and per level lane ``present``,
    ``profile``, ``count``, ``opportunity``, ``has_witness``, ``rna_count``, ``rna_count_var``). Then the backbone's
    assertion counts and the name of the policy that ran."""

    fg_loc: np.ndarray | None = None
    f_g: np.ndarray | None = None
    tau_lam: np.ndarray | None = None
    solvable: np.ndarray | None = None
    count: np.ndarray | None = None
    free_pos: np.ndarray | None = None
    free_neg: np.ndarray | None = None
    mass_global: np.ndarray | None = None
    lam_rows: np.ndarray | None = None
    cube_rows: CubeRows | None = None
    from_left: dict | None = None
    from_right: dict | None = None
    backbone_assertions: AssertionCounts | None = None
    policy_name: str | None = None


@dataclass(slots=True)
class CalibrationDiagnostics:
    """The solve's intermediates, filled by ``calibrate``: the chain, its geometry and statics, the substrate, the sj
    axis and region frame it ran on, the final belief (after the refits), the last refit's landscape (``None`` when no
    refit fitted one), the library's RNA sense fraction, the intergenic background's law as the sweep read it (its
    departure prior the introns'), the fit of that departure prior, the locks' departure on the slot axis, and the last
    sweep's :class:`SweepDiagnostics`, whose counts are the sweep's own: a locked slot's at its gDNA share."""

    chain: RegionChain | None = None
    geometry: RegionGeometry | None = None
    statics: RegionStatics | None = None
    substrate: CalibrationSubstrate | None = None
    sj: SpliceJunctionGeometry | None = None
    region_arrays: RegionArrays | None = None
    belief: RegionBelief | None = None
    landscape: DensityLandscape | None = None
    rna_sense_frac: float | None = None
    background: BackgroundLaw | None = None
    departure: DepartureFit | None = None
    lock: LockFit | None = None
    sweep: SweepDiagnostics = field(default_factory=SweepDiagnostics)
