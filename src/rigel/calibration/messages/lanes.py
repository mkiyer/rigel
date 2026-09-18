"""The LEVEL LANES — the gDNA lane and one RNA lane per strand — the class they are built on.

       Gate: ``tests/calibration/test_transfer_rna_lanes.py``, ``test_transfer_policy.py``

A level is a population's density as an ABSOLUTE profile over ``u = log(rho / rho_ref)`` on the solve
grid, so it needs no map and no recipient and crosses the faces composition cannot (`LevelLane`). The
gDNA lane serves every directed face left without a composition rule (`faces.Faces`); each RNA lane
serves its strand's faces read off the flag bits. Both are built by `native.transfer_prepare`
(`native/transfer_kernel.cpp`'s ``gdna_lane`` and ``rna_lane``) into the tables `TransferPolicy.prepare`
allocates on a `LevelLane` — its faces, own levels, junction flux levels and flux witnesses — and carried
across each face by the pass kernel (``lane_emit`` and ``lane_receive``, the same file); its coordinates,
witnesses and emptiness are the chain's arrays.
"""

from __future__ import annotations

import numpy as np

from . import Received
from .faces import RowTable, side_of

__all__ = ["LevelLane"]


_FIELD = {"gdna": "level_gdna", "pos": "level_rna_pos", "neg": "level_rna_neg"}


class LevelLane:
    """ONE population's LEVEL LANE — the gDNA lane, or one strand's RNA lane — for one sweep, as the
    tables the native pass reads. A level travels as an ABSOLUTE profile over ``u = log(rho / rho_ref)``
    on the solve grid, so it needs no map and no recipient: an EMPTY node forwards it unchanged, a full
    node emits the INTERSECTION of its own level and what it holds — ⛔ bounds INTERSECT, they do not
    multiply: a product of one-sided claims ratchets along a chain into a hard bound at the noisiest
    node's mode — and a full recipient prices the hop and takes the level as a LOWER bound; the hop's
    price is both witness counts' counting plus the disagreement beyond it between the two nodes'
    estimates of the population's abundance, per hop and nothing pooled (``lane_emit`` and
    ``lane_receive`` in `native/transfer_kernel.cpp`).

    The two kinds are one class because they differ by PARAMETERS and nothing else: the WITNESS the
    hop's price reads — every node's ``count`` on its opportunity ``a`` (the totals on the gDNA lane;
    on an RNA lane the strand's read column, and where the library's strand channel is live the
    column split's asymmetry against ``other``, the other column — the column holds gDNA's half, which
    jumps with every probe edge whether or not this strand's RNA is there, while the asymmetry is this
    strand's RNA alone) — and ``two_sided``, the faces the WHOLE profile crosses (none on the gDNA lane:
    every gDNA hop is a lower bound; an intron and its own boundary on an RNA lane: one shared
    unspliced population). ``total`` is every node's total, through which the solve reads a held level
    back as the node's composition; ``own_level`` each node's own level (a `RowTable`); ``face`` the
    directed faces the lane serves; ``flux``, per ``(exon, side)`` — the junction's priced estimate of
    the exon's RNA, kept per FACE (a `RowTable` over the node's two sides) so the solve can tell which
    face's composition already carries it; ``flux_witness``, per EMPTY node whose own level is a flux
    level, the witness that level travels with — the pooled spliced count on the pooled route
    opportunity of the node's lit junctions (an empty node has no count of its own to stamp a level
    with), a `RowTable` of ``(count, opportunity)`` pairs. Every array is in the layout the native pass
    reads (`tables`)."""

    __slots__ = (
        "population",
        "field",
        "u",
        "lam",
        "rho_ref",
        "count",
        "total",
        "a",
        "empty",
        "own_level",
        "face",
        "two_sided",
        "flux",
        "other",
        "flux_witness",
    )

    def __init__(
        self,
        population: str,
        u,
        lam,
        rho_ref,
        count,
        a,
        empty,
        own_level,
        face,
        *,
        two_sided=None,
        total=None,
        flux=None,
        other=None,
        flux_witness=None,
    ):
        n = len(own_level)
        self.population, self.field = population, _FIELD[population]
        self.u, self.lam, self.rho_ref = u, lam, float(rho_ref)
        self.count, self.a, self.empty = count, a, empty
        self.total = count if total is None else total
        self.own_level = own_level
        # the faces the lane serves and the two-sided ones, per (destination, side)
        self.face = np.asarray(face, bool).reshape(n, 2)
        self.two_sided = (
            np.zeros((n, 2), bool) if two_sided is None else np.asarray(two_sided, bool)
        )
        # the junction flux levels, per (exon, side of its junction)
        self.flux = RowTable((n, 2), self.own_level.rows.shape[1]) if flux is None else flux
        self.other = other
        self.flux_witness = RowTable(n, 2) if flux_witness is None else flux_witness

    def tables(self) -> tuple:
        """The lane as the native pass reads it (`native.transfer_pass`'s ``lanes`` entry): which
        `Received` table it writes, its faces, its two-sided faces, its emptiness, the own levels with
        their mask, the witness counts and opportunities, the other column or ``None``, and the flux
        witnesses with their mask — the arrays themselves, nothing copied."""
        return (
            Received.LANES.index(self.field),
            self.face,
            self.two_sided,
            self.empty,
            self.own_level.rows,
            self.own_level.mask,
            self.count,
            self.a,
            self.other,
            self.flux_witness.rows,
            self.flux_witness.mask,
        )

    def serves(self, s: int, x: int) -> bool:
        """Does the lane carry its level across the face into ``x`` from ``s``?"""
        return bool(self.face[int(x), side_of(int(s), int(x))])

    def flux_at(self, x: int, side: int):
        """The junction flux level ``x`` holds at its ``side`` (0: from its left junction, 1: its right),
        or ``None``."""
        return self.flux[int(x), int(side)]
