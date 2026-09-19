"""The diagnostics capture of the locus solve — the instruments' view of a sweep — and the chain's arrays
as a policy may read them.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (the block tests)

The chain is solved a LOCUS BLOCK at a time in one native call (`sweep.solve_chain`,
`region_chain.locus_blocks`). Two things live here and nothing about what a solve or a message is:
`view_fields` is every per-slot array a policy may read (:class:`~.messages.ChainView`), for the whole
chain; `SweepCapture` is the diagnostic capture of a sweep — the instruments' view of the solve, filled by
the backbone from the kernel's capture arrays.
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
    :class:`~.messages.ChainView` — the whole chain, which is what the kernel reads a block at a time."""
    return dict(
        # observations
        eff_gdna=np.asarray(geometry.eff_gdna, np.float64),
        eff_rna=np.asarray(geometry.eff_rna, np.float64),
        sj_count=np.asarray(geometry.sj_count, np.float64),  # [n, 2] by TRANSCRIPT strand
        sj_count_lo=np.asarray(geometry.sj_count_lo, np.float64),
        sj_count_hi=np.asarray(geometry.sj_count_hi, np.float64),
        route_rate_lo=np.asarray(geometry.route_rate_lo, np.float64),
        route_rate_hi=np.asarray(geometry.route_rate_hi, np.float64),
        unspliced_count=np.asarray(geometry.unspliced_count, np.float64),  # [n, 2] by GENOME strand
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
    """The diagnostic capture of one sweep — the instruments' view of the solve, never read by
    production (`pipeline` passes no debug dict, so the sweep never asks the kernel for it). A caller
    hands `sweep.solve_chain` an empty record; the sweep fills it from the kernel's capture mode. Every
    field has a reader in an instrument or a test; a field nothing reads is dead surface and does not
    belong here.

    Per slot (``(n,)`` unless stated): the message-free SELF-SOLVE's ``fg_loc`` and its ``tau_lam`` (the
    slot's own composition evidence: the strand term's precision at ``fg_loc`` plus the factory row's
    curvature, ``tau_fac``, which is published apart so an instrument can split the two sources); the
    STRAND-ONLY solve's ``fg_strand`` (no prior, no messages, to split the local error into the strand
    likelihood against the prior's contribution); the FINAL solve's ``f_g`` and ``var_g``; the incoming
    belief's ``fg_init``; ``solvable``; the observations ``count`` ``(n, 2)``, ``spliced``, ``mature``,
    ``free_pos``, ``free_neg``, ``eff_gdna``, ``eff_rna``, and the gDNA support ``mass_global``,
    ``eff_global``. The message layer's delivery: ``lam_rows`` ``(n, K)`` (zero rows where a block
    delivered nothing, ``None`` when no block did), ``cube_rows`` (a `CubeRows` table keyed to the chain),
    and the two received tables ``from_left`` / ``from_right`` — each a dict of the kernel's arrays over
    the chain (``has_neighbour``, ``has_composition``, ``composition``, and per level lane ``present``,
    ``profile``, ``count``, ``opportunity``, ``has_witness``, ``rna_count``, ``rna_count_var``). The chain:
    its adjacency ``left`` / ``right``, the backbone's assertion counts, the name of the policy that RAN
    (the witness an instrument's "the arm ran" assertion needs), and the solve grid ``f_g = σ(λ)``."""

    fg_loc: np.ndarray | None = None
    fg_strand: np.ndarray | None = None
    f_g: np.ndarray | None = None
    var_g: np.ndarray | None = None
    tau_lam: np.ndarray | None = None
    tau_fac: np.ndarray | None = None
    fg_init: np.ndarray | None = None
    solvable: np.ndarray | None = None
    count: np.ndarray | None = None
    spliced: np.ndarray | None = None
    mature: np.ndarray | None = None
    free_pos: np.ndarray | None = None
    free_neg: np.ndarray | None = None
    eff_gdna: np.ndarray | None = None
    eff_rna: np.ndarray | None = None
    mass_global: np.ndarray | None = None
    eff_global: np.ndarray | None = None
    lam_rows: np.ndarray | None = None
    cube_rows: CubeRows | None = None
    from_left: dict | None = None
    from_right: dict | None = None
    left: np.ndarray | None = None
    right: np.ndarray | None = None
    backbone_assertions: object = None
    policy_name: str | None = None
    solve_grid: np.ndarray | None = None

    @staticmethod
    def concat_received(parts: list) -> dict:
        """The blocks' received tables as the chain's: each array concatenated over the blocks' owned
        slots, in block order (a dict of arrays, its lanes dicts of arrays)."""
        out = {}
        for name, first in parts[0].items():
            if isinstance(first, dict):
                out[name] = SweepCapture.concat_received([p[name] for p in parts])
            else:
                out[name] = np.concatenate([np.asarray(p[name]) for p in parts], axis=0)
        return out

    def fill(self, other: "SweepCapture") -> None:
        """Every field of ``other`` copied onto this record — how the backbone fills the record a
        caller handed it."""
        for f in dataclasses.fields(self):
            setattr(self, f.name, getattr(other, f.name))
