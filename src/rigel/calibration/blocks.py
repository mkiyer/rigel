"""The block plumbing of the locus solve — how one block of the chain is cut out of it and how its
results are put back.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (the block tests)

The chain is solved a LOCUS BLOCK at a time (`sweep.solve_chain`, `region_chain.locus_blocks`). Two
helpers and one record, and nothing about what a solve or a message is: `view_fields` is every per-slot
array a policy may read, for the whole chain or one block; `block_slice` restricts a chain, its statics,
its geometry or a belief to one block's slots; `SweepCapture` is the diagnostic capture of a sweep — the
instruments' view of the solve, a block's or, through `SweepCapture.gather`, the chain's.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

from .messages import Received
from .region_chain import REGION, RegionChain

__all__ = ["SweepCapture", "block_slice", "view_fields"]


def view_fields(chain, statics, geometry, structure) -> dict:
    """Every per-slot array a policy may read, under the two belief-free headings of
    :class:`~.messages.ChainView` — for the whole chain (the library) or one block (the context)."""
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


def block_slice(obj, sl: slice):
    """``obj`` — a chain, its statics, its geometry or a belief — restricted to the slots ``sl``: every
    per-slot array sliced (a view, no copy), ``n_slots`` updated, and a chain's links re-based so that
    a neighbour outside the block is no neighbour (``-1``). The block's first slot is a terminal or a
    reference start and its last read slot a terminal, so a link that leaves the block always leaves
    it at a node that receives nothing, and what that node holds is not read by any consumer."""
    fields = {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
    arrays = {k: v for k, v in fields.items() if isinstance(v, np.ndarray) and v.ndim >= 1}
    n = next(iter(arrays.values())).shape[0]
    new = {k: v[sl] for k, v in arrays.items() if v.shape[0] == n}
    if "n_slots" in fields:
        new["n_slots"] = len(range(*sl.indices(n)))
    if isinstance(obj, RegionChain):
        m = new["kind"].shape[0]
        for k in ("left", "right"):
            a = new[k] - sl.start
            new[k] = np.where((a >= 0) & (a < m), a, -1)
    return dataclasses.replace(obj, **new)


@dataclass(slots=True)
class SweepCapture:
    """The diagnostic capture of one sweep — the instruments' view of the solve, never read by
    production (`pipeline` passes no debug dict, so the sweep never builds one). A caller hands
    `sweep.solve_chain` an empty record; the sweep fills it, one block at a time and then gathered to the
    chain. Every field has a reader in an instrument or a test; a field nothing reads is dead surface
    and does not belong here.

    Per slot (``(n,)`` unless stated): the message-free SELF-SOLVE's ``fg_loc`` and its ``tau_lam`` (the
    slot's own composition evidence, `RegionInit.tau_lam`); the STRAND-ONLY solve's ``fg_strand`` (no
    prior, no messages, to split the local error into the strand likelihood against the prior's
    contribution); the FINAL solve's ``f_g`` and ``var_g``; the incoming belief's ``fg_init``;
    ``solvable``; the observations ``count`` ``(n, 2)``, ``spliced``, ``mature``, ``free_pos``,
    ``free_neg``, ``eff_gdna``, ``eff_rna``, and the gDNA support ``mass_global``, ``eff_global``.
    The message layer's delivery: ``lam_rows`` ``(n, K)`` (zero rows where a block delivered nothing,
    ``None`` when no block did), ``cube_rows`` ``{slot: CubeRow}`` keyed to the chain, and the two
    :class:`~.messages.Received` tables ``from_left`` / ``from_right``. The chain: its adjacency
    ``left`` / ``right``, the backbone's assertion counts, the name of the policy that RAN (the witness
    an instrument's "the arm ran" assertion needs — `TRAPS: an-ablation-that-never-ran`), the solve
    grid ``f_g = σ(λ)``, and the intron factory's rows."""

    fg_loc: np.ndarray | None = None
    fg_strand: np.ndarray | None = None
    f_g: np.ndarray | None = None
    var_g: np.ndarray | None = None
    tau_lam: np.ndarray | None = None
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
    cube_rows: dict | None = None
    from_left: Received | None = None
    from_right: Received | None = None
    left: np.ndarray | None = None
    right: np.ndarray | None = None
    backbone_assertions: object = None
    policy_name: str | None = None
    solve_grid: np.ndarray | None = None
    intron_prior: np.ndarray | None = None

    #: the per-slot arrays, concatenated over each block's OWNED slots when the chain's is gathered
    PER_SLOT = (
        "fg_loc",
        "fg_strand",
        "f_g",
        "var_g",
        "tau_lam",
        "fg_init",
        "solvable",
        "count",
        "spliced",
        "mature",
        "free_pos",
        "free_neg",
        "eff_gdna",
        "eff_rna",
        "mass_global",
        "eff_global",
    )

    @classmethod
    def gather(cls, blocks: list) -> "SweepCapture":
        """The blocks' captures as the chain's: ``blocks`` is ``[(locus_block, capture), ...]`` in
        chain order. Per-slot arrays are concatenated over each block's owned slots, the delivered rows
        zero-filled where a block delivered nothing (``None`` only when no block delivered), the cube
        rows re-keyed to the chain, the two tables concatenated. The chain-level fields are the
        backbone's to set."""
        out = cls()
        if not blocks:
            return out
        owned = [(b.stop - b.start, c) for b, c in blocks]
        for name in cls.PER_SLOT:
            setattr(out, name, np.concatenate([np.asarray(getattr(c, name))[:k] for k, c in owned]))
        out.from_left = Received.concat([c.from_left.take(k) for k, c in owned])
        out.from_right = Received.concat([c.from_right.take(k) for k, c in owned])
        if any(c.lam_rows is not None for _k, c in owned):
            K = next(np.asarray(c.lam_rows).shape[1] for _k, c in owned if c.lam_rows is not None)
            out.lam_rows = np.concatenate(
                [
                    np.zeros((k, K)) if c.lam_rows is None else np.asarray(c.lam_rows)[:k]
                    for k, c in owned
                ]
            )
        cube: dict = {}
        for b, c in blocks:
            for slot, row in (c.cube_rows or {}).items():
                if int(slot) < b.stop - b.start:
                    cube[int(slot) + b.start] = row
        out.cube_rows = cube or None
        return out

    def fill(self, other: "SweepCapture") -> None:
        """Every field of ``other`` copied onto this record — how the backbone fills the record a
        caller handed it."""
        for f in dataclasses.fields(self):
            setattr(self, f.name, getattr(other, f.name))
