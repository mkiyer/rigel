"""The chain's arrays as a policy may read them, and the diagnostics capture of a sweep.

:func:`view_fields` is every per-slot array a policy may read (:class:`~.messages.ChainView`), for the whole chain;
:class:`SweepCapture` is the instruments' view of one sweep, filled by the backbone (``sweep.solve_chain``) from the
kernel's capture arrays. Neither knows what a solve or a message is.

Gate: ``tests/calibration/test_sweep_backbone.py`` (the block tests).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

from .region_chain import REGION
from .simplex_logodds import CubeRows

__all__ = ["SweepCapture", "view_fields"]


def view_fields(chain, statics, geometry, structure) -> dict:
    """Every per-slot array a policy may read, under the two belief-free headings of
    :class:`~.messages.ChainView`, for the whole chain; the kernel reads it a block at a time."""
    return dict(
        # observations
        eff_gdna=np.asarray(geometry.eff_gdna, np.float64),
        eff_rna=np.asarray(geometry.eff_rna, np.float64),
        sj_count=np.asarray(geometry.sj_count, np.float64),  # [n, 2] by transcript strand
        sj_count_lo=np.asarray(geometry.sj_count_lo, np.float64),
        sj_count_hi=np.asarray(geometry.sj_count_hi, np.float64),
        route_rate_lo=np.asarray(geometry.route_rate_lo, np.float64),
        route_rate_hi=np.asarray(geometry.route_rate_hi, np.float64),
        unspliced_count=np.asarray(geometry.unspliced_count, np.float64),  # [n, 2] by genome strand
        spliced_count=np.asarray(geometry.spliced_count, np.float64),
        # geometry / structure
        left=np.asarray(chain.left, np.int64),
        right=np.asarray(chain.right, np.int64),
        is_boundary=np.asarray(chain.kind) != REGION,
        is_exon_region=np.asarray(structure.is_exon_region, bool),
        free_pos=np.asarray(statics.free_pos, bool),
        free_neg=np.asarray(statics.free_neg, bool),
        exon_pos=np.asarray(structure.exon_pos, bool),
        exon_neg=np.asarray(structure.exon_neg, bool),
        boundary_flags=statics.boundary_flags,
    )


@dataclass(slots=True)
class SweepCapture:
    """The diagnostic capture of one sweep, for instruments and tests only: production passes no capture, so the
    sweep never asks the kernel for one. A caller hands ``sweep.solve_chain`` an empty record and the sweep fills it.

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
    backbone_assertions: object = None
    policy_name: str | None = None

    @staticmethod
    def concat_received(parts: list) -> dict:
        """The blocks' received tables as the chain's: each array concatenated over the blocks' owned slots in block
        order (a dict of arrays, its lanes dicts of arrays)."""
        out = {}
        for name, first in parts[0].items():
            if isinstance(first, dict):
                out[name] = SweepCapture.concat_received([p[name] for p in parts])
            else:
                out[name] = np.concatenate([np.asarray(p[name]) for p in parts], axis=0)
        return out

    def fill(self, other: "SweepCapture") -> None:
        """Copy every field of ``other`` onto this record."""
        for f in dataclasses.fields(self):
            setattr(self, f.name, getattr(other, f.name))
