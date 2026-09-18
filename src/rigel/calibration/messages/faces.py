"""The composition rules as TYPED TABLES — what a recipient does with a composition that arrives at
each of its two sides — and the three small functions every reader of a face needs.

       Gate: ``tests/calibration/test_transfer_faces.py``

A directed face is ``(destination, side)``: a node hears from its LEFT neighbour (side 0, the forward
pass) or its RIGHT (side 1, the backward pass), so a rule needs no pair and no lookup (`side_of`). A
rule is a KIND and its parameters; `Faces.apply` is the one home of the rule arithmetic, on the row
constructors of `transfer_rows`. The builders that write the rules live in `transfer`; the level lanes
that serve every face left without one live in `lanes`. `RowTable` is the one shape of an OPTIONAL row
per node — a claim, a level, a witness — as the native pass reads it: a matrix and a presence mask,
which the builders write directly.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .transfer_rows import EPS, blur_row, level_row, splice_out_row, transport_row

__all__ = [
    "EDGE",
    "FORWARD",
    "LEVEL",
    "NONE",
    "RULE_NAMES",
    "SPLICE_OUT",
    "TRANSPORT",
    "FaceRule",
    "Faces",
    "RowTable",
    "fuse",
    "norm",
    "side_of",
]


def norm(row):
    row = np.asarray(row, np.float64)
    return row - row.max()


def fuse(parts):
    """Independent witnesses about one slot: log-profiles add, then re-normalise."""
    out = None
    for p in parts:
        out = p if out is None else out + p
    return None if out is None else norm(out)


def side_of(s: int, i: int) -> int:
    """The side of ``i`` a hop from ``s`` arrives on: ``0`` from its LEFT neighbour (the forward pass),
    ``1`` from its RIGHT (the backward pass). Slot ids are genomic order, so the left neighbour is the
    smaller id — which is why a directed face needs no pair and no lookup: it IS ``(destination,
    side)``, and every table below is indexed that way."""
    return 0 if s < i else 1


class RowTable:
    """OPTIONAL ROWS over the nodes, in the layout the native pass reads: an ``(n, K)`` matrix and an
    ``(n,)`` presence mask, written by the builders directly. ``t[i]`` is node ``i``'s row (a view into
    the matrix) or ``None``; ``t[i] = row`` writes it and marks it present; ``t[i] = None`` clears it.
    THE MASK, never the matrix, says whether a node has a row: the matrix is allocated UNFILLED
    (``np.empty`` — a sweep allocates these tables per block, and filling them was a cost with no
    reader), so an absent row's cells are unspecified and every reader — the kernels, the solve, the
    gates — reads the mask before the row. ``shape`` may be ``(n, 2)`` for a row per directed face —
    ``t[x, side]`` — as the lanes' junction flux levels are kept."""

    __slots__ = ("rows", "mask")

    def __init__(self, shape, K: int):
        shape = (int(shape),) if np.ndim(shape) == 0 else tuple(int(s) for s in shape)
        self.rows = np.empty((*shape, int(K)))
        self.mask = np.zeros(shape, bool)

    def __len__(self) -> int:
        return self.mask.shape[0]

    def __getitem__(self, i):
        return self.rows[i] if self.mask[i] else None

    def __setitem__(self, i, row) -> None:
        if row is None:
            self.mask[i] = False
        else:
            self.rows[i] = row
            self.mask[i] = True


#: the five KINDS of composition rule a directed face can carry (``NONE``: the face has no rule and the
#: level lane serves it). Each is one arithmetic of `transfer_rows` with the parameters `Faces` holds.
NONE, FORWARD, TRANSPORT, SPLICE_OUT, EDGE, LEVEL = range(6)
RULE_NAMES = ("none", "forward", "transport", "splice_out", "edge", "level")


class FaceRule(NamedTuple):
    """One face's rule read out of the tables (`Faces.at`): its kind and parameters, the rows resolved."""

    kind: int
    n_u: float
    n_s: float
    a_b: float
    a_x: float
    width: float
    var: float
    row: "np.ndarray | None"
    row2: "np.ndarray | None"


class Faces:
    """THE COMPOSITION RULES AS TYPED TABLES — the port's data layout for what a recipient does with
    a composition that arrives at each of its two sides.

    Every array is ``(n, 2)`` over ``(destination, side)``: ``kind`` (one of the five above), the
    scalar parameters ``n_u`` / ``n_s`` / ``a_b`` / ``a_x`` (the face's unspliced and spliced counts,
    the boundary's and the far region's gDNA opportunity), ``width`` (a blur variance beyond counting)
    and ``var`` (the level rule's width), and ``row`` / ``row2``, indices into ``rows`` — the ``(K,)``
    maps a rule needs (a face map's λ image, a level map, the edge's level, the level rule's bound),
    a matrix of which the first ``n_rows`` are written.
    ``nbr[i, side]`` names the neighbour each face hears from, so a builder cannot write a rule at a
    face that does not exist; and a face carries ONE rule — a second write is refused — because the
    builders' faces are disjoint by construction (the splice faces serve intron|exon pairs, the edge
    rule gene edges, the terminus rules unlicensed faces, the alternative splice site junctions with
    no terminus), and a precedence that nothing exercises is a hidden assumption, not a rule.

    :meth:`apply` is the one place a rule's arithmetic lives:

    ==========  ===============================================================================
    FORWARD     the identity: the sender's claim and what it holds, fused
    TRANSPORT   boundary → region through the face map ``rows[row]`` (`transport_row`), blurred by
                ``width`` where a pair's discrepancy adds one
    SPLICE_OUT  region → boundary, the face map read backwards (`splice_out_row`), likewise
    EDGE        the intergenic|exon edge's one-sided level ``rows[row]``, a constant; a flat one
                is no claim
    LEVEL       the terminus's level rule: the sender's OWN claim through the level map
                ``rows[row]`` at width ``var`` (`level_row`), or the crossing total's bound
                ``rows[row2]`` when it has none — never what it holds
    ==========  ===============================================================================
    """

    __slots__ = (
        "lam",
        "nbr",
        "kind",
        "row",
        "row2",
        "n_u",
        "n_s",
        "a_b",
        "a_x",
        "width",
        "var",
        "rows",
        "n_rows",
    )

    def __init__(self, lam, left, right):
        n = int(np.asarray(left).shape[0])
        self.lam = np.asarray(lam, np.float64)
        self.nbr = np.stack((np.asarray(left, np.int64), np.asarray(right, np.int64)), axis=1)
        self.kind = np.zeros((n, 2), np.int8)
        self.row = np.full((n, 2), -1, np.int32)
        self.row2 = np.full((n, 2), -1, np.int32)
        self.n_u = np.zeros((n, 2))
        self.n_s = np.zeros((n, 2))
        self.a_b = np.zeros((n, 2))
        self.a_x = np.zeros((n, 2))
        self.width = np.zeros((n, 2))
        self.var = np.zeros((n, 2))
        # the row store, as the native pass reads it: a rule keeps at most two rows (a map and a
        # bound) and a face carries one rule, so two rows per face is its capacity; the unwritten
        # rows are never read
        self.rows = np.empty((2 * int((self.nbr >= 0).sum()), self.lam.shape[0]))
        self.n_rows = 0

    def set(
        self,
        s,
        i,
        kind,
        *,
        row=None,
        row2=None,
        n_u=0.0,
        n_s=0.0,
        a_b=0.0,
        a_x=0.0,
        width=0.0,
        var=0.0,
    ):
        """The rule at the face into ``i`` from ``s`` — which must be ``i``'s neighbour on that side."""
        s, i = int(s), int(i)
        side = side_of(s, i)
        if self.nbr[i, side] != s:
            raise ValueError(
                f"no face into {i} from {s}: its neighbour on that side is {self.nbr[i, side]}"
            )
        if self.kind[i, side] != NONE:
            raise ValueError(
                f"the face into {i} from {s} already carries a {RULE_NAMES[self.kind[i, side]]} rule; "
                f"a face carries one rule, and the builders' faces are disjoint"
            )
        self.kind[i, side] = kind
        self.row[i, side] = self._keep(row)
        self.row2[i, side] = self._keep(row2)
        self.n_u[i, side], self.n_s[i, side] = float(n_u), float(n_s)
        self.a_b[i, side], self.a_x[i, side] = float(a_b), float(a_x)
        self.width[i, side], self.var[i, side] = float(width), float(var)

    def _keep(self, row) -> int:
        if row is None:
            return -1
        k = self.n_rows
        self.rows[k] = row
        self.n_rows = k + 1
        return k

    def any(self) -> bool:
        return bool((self.kind != NONE).any())

    def kind_at(self, s, i) -> int:
        return int(self.kind[int(i), side_of(int(s), int(i))])

    def has(self, s, i) -> bool:
        return self.kind_at(s, i) != NONE

    def at(self, s, i) -> "FaceRule":
        """The rule at one face as a record — its kind, parameters and resolved rows — for a reader or a
        gate; the passes read the tables directly."""
        i, side = int(i), side_of(int(s), int(i))
        r, r2 = int(self.row[i, side]), int(self.row2[i, side])
        return FaceRule(
            int(self.kind[i, side]),
            float(self.n_u[i, side]),
            float(self.n_s[i, side]),
            float(self.a_b[i, side]),
            float(self.a_x[i, side]),
            float(self.width[i, side]),
            float(self.var[i, side]),
            None if r < 0 else self.rows[r],
            None if r2 < 0 else self.rows[r2],
        )

    def pairs(self) -> list:
        """Every directed face ``(s, i)`` that carries a rule, in table order."""
        i_idx, sides = np.nonzero(self.kind != NONE)
        return [(int(self.nbr[i, side]), int(i)) for i, side in zip(i_idx.tolist(), sides.tolist())]

    def apply(self, s, i, own, held):
        """The rule at the face into ``i`` from ``s``, applied to what ``s`` sends: its OWN claim and what
        it HOLDS from its far side — each ``None`` where absent. Returns the row for ``i``, or ``None``
        for no claim."""
        i, side = int(i), side_of(int(s), int(i))
        k = self.kind[i, side]
        if k == FORWARD:
            return fuse([r for r in (own, held) if r is not None])
        if k == EDGE:
            r = self.rows[self.row[i, side]]
            return r if np.ptp(r) > EPS else None
        if k == LEVEL:
            if own is None:
                return self.rows[self.row2[i, side]]
            return level_row(own, self.lam, self.rows[self.row[i, side]], float(self.var[i, side]))
        sending = fuse([r for r in (own, held) if r is not None])
        if sending is None:
            return None
        n_u, n_s = float(self.n_u[i, side]), float(self.n_s[i, side])
        if k == TRANSPORT:
            out = transport_row(sending, self.lam, self.rows[self.row[i, side]], n_u, n_s)
        elif k == SPLICE_OUT:
            out = splice_out_row(
                sending, self.lam, n_u, n_s, float(self.a_b[i, side]), float(self.a_x[i, side])
            )
        else:
            return None
        w = float(self.width[i, side])
        return blur_row(out, self.lam, w) if w > 0.0 else out
