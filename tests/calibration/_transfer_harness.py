"""The shared harness of the transfer policy's gate files (`test_transfer_*.py`): ONE real
`solve_chain` call captured from a calibrate run on the toy — the backbone-parity pattern, so every
gate re-runs the sweep with a different policy on byte-identical inputs — plus the containers the gates
read the kernel's tables through and the drivers they share: the builders on one block's context, whole
passes and single hops of the native pass, the solve. The kernel is the ONE implementation
(`native/transfer_kernel.h`, run per block by `native/solve_kernel.cpp`); the bindings `transfer_prepare`,
`transfer_pass` and `transfer_solve` run it on tables the call allocates, and the containers here —
`RowTable`, `Faces`, `LevelLane`, `Received`, `Levels` — hold and index those arrays and compute nothing.
Not a test module: `capture_sweep_inputs` is what each gate file's ``sweep_inputs`` fixture returns, and
nothing here asserts."""

from __future__ import annotations

import dataclasses
import sys
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

import rigel.calibration.sweep as SW
from rigel.calibration.blocks import SweepCapture, view_fields
from rigel.calibration.messages import ChainView
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION
from rigel.calibration.region_init import strand_discriminability
from rigel.calibration.simplex_logodds import CubeRows
from rigel.native import transfer_pass, transfer_prepare, transfer_rows, transfer_solve

#: the five KINDS of composition rule a directed face can carry (``NONE``: the face has no rule and the
#: level lane serves it) — the kernel's `transfer_rows.h` constants
NONE, FORWARD, TRANSPORT, SPLICE_OUT, EDGE, LEVEL = range(6)
RULE_NAMES = ("none", "forward", "transport", "splice_out", "edge", "level")

#: the received tables' three level lanes, by field name, in the kernel's order
LANES = ("level_gdna", "level_rna_pos", "level_rna_neg")
_FIELD = {"gdna": "level_gdna", "pos": "level_rna_pos", "neg": "level_rna_neg"}


def side_of(s: int, i: int) -> int:
    """The side of ``i`` a hop from ``s`` arrives on: ``0`` from its LEFT neighbour (the forward pass),
    ``1`` from its RIGHT (the backward pass)."""
    return 0 if s < i else 1


@dataclass(frozen=True, slots=True, kw_only=True)
class BlockContext(ChainView):
    """One block of the chain as the kernel's builders read it: the :class:`ChainView` plus the two
    source-side inputs — every node's own-evidence bit and the incoming belief — and the factory's rows."""

    has_own_composition: np.ndarray
    belief_fg: np.ndarray
    factory_rows: np.ndarray | None = None


# ── the containers: the kernel's tables, held and indexed ───────────────────────────────────────────────


class RowTable:
    """OPTIONAL ROWS over the nodes as the kernel keeps them: an ``(n, K)`` matrix (or ``(n, 2, K)`` for a
    row per directed face) and a presence mask; ``t[i]`` is the row or ``None``, ``t[i] = row`` writes it
    and marks it present, ``t[i] = None`` clears it. The matrix is unfilled where the mask is False."""

    __slots__ = ("rows", "mask")

    def __init__(self, shape, K: int, rows=None, mask=None):
        shape = (int(shape),) if np.ndim(shape) == 0 else tuple(int(s) for s in shape)
        self.rows = np.empty((*shape, int(K))) if rows is None else rows
        self.mask = np.zeros(shape, bool) if mask is None else mask

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


class FaceRule(NamedTuple):
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
    """The composition rules as TYPED TABLES over ``(destination, side)`` — the kernel's ``faces`` dict:
    ``kind``, the scalar parameters, ``row`` / ``row2`` as indices into the row store ``rows``."""

    def __init__(self, lam, left, right, tables: dict | None = None):
        n = int(np.asarray(left).shape[0])
        self.lam = np.ascontiguousarray(lam, np.float64)
        self.nbr = np.stack((np.asarray(left, np.int64), np.asarray(right, np.int64)), axis=1)
        if tables is None:  # a table with no rule anywhere
            tables = dict(
                kind=np.zeros((n, 2), np.int8),
                row=np.full((n, 2), -1, np.int32),
                row2=np.full((n, 2), -1, np.int32),
                rows=np.zeros((0, self.lam.shape[0])),
                **{k: np.zeros((n, 2)) for k in ("n_u", "n_s", "a_b", "a_x", "width", "var")},
            )
        self.tables = tables
        for k in ("kind", "row", "row2", "n_u", "n_s", "a_b", "a_x", "width", "var", "rows"):
            setattr(self, k, tables[k])

    @property
    def n_rows(self) -> int:
        return int(self.rows.shape[0])

    def any(self) -> bool:
        return bool((self.kind != NONE).any())

    def kind_at(self, s, i) -> int:
        return int(self.kind[int(i), side_of(int(s), int(i))])

    def has(self, s, i) -> bool:
        return self.kind_at(s, i) != NONE

    def at(self, s, i) -> FaceRule:
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


class LevelLane:
    """ONE population's LEVEL LANE as the kernel keeps it — its faces and two-sided faces per
    ``(destination, side)``, every node's emptiness, own level, witness count and opportunity, the other
    column where the strand channel is live, the junction flux levels per ``(exon, side)`` in a row store
    with an index, and the flux witnesses. Built from the kernel's ``lanes`` dict (`from_kernel`) or by
    hand from arrays, and handed back to the pass and the solve as that dict (`as_kernel`)."""

    def __init__(
        self,
        population: str,
        u,
        lam,
        rho_ref,
        count,
        a,
        empty,
        own_level: RowTable,
        face,
        *,
        two_sided=None,
        total=None,
        flux: RowTable | None = None,
        other=None,
        flux_witness: RowTable | None = None,
        flux_rows=None,
        flux_index=None,
    ):
        n = len(own_level)
        self.population, self.field = population, _FIELD[population]
        self.u, self.lam, self.rho_ref = u, lam, float(rho_ref)
        self.count, self.a, self.empty = count, a, np.asarray(empty, bool)
        self.total = count if total is None else total
        self.own_level = own_level
        self.face = np.asarray(face, bool).reshape(n, 2)
        self.two_sided = (
            np.zeros((n, 2), bool) if two_sided is None else np.asarray(two_sided, bool)
        )
        self.other = other
        self.flux_witness = RowTable(n, 2) if flux_witness is None else flux_witness
        K = own_level.rows.shape[-1]
        if flux is not None:  # a hand-built (n, 2, K) table into the kernel's store and index
            idx = np.full((n, 2), -1, np.int32)
            store = []
            for x, side in zip(*np.nonzero(flux.mask)):
                idx[x, side] = len(store)
                store.append(np.asarray(flux.rows[x, side], np.float64))
            self.flux_rows = np.array(store).reshape(len(store), K)
            self.flux_index = idx
        else:
            self.flux_rows = np.zeros((0, K)) if flux_rows is None else flux_rows
            self.flux_index = np.full((n, 2), -1, np.int32) if flux_index is None else flux_index

    @classmethod
    def from_kernel(cls, population: str, d: dict, lam) -> "LevelLane":
        n = int(d["own_mask"].shape[0])
        return cls(
            population,
            lam,
            lam,
            float(d["rho_ref"]),
            d["count"],
            d["a"],
            d["empty"],
            RowTable(n, lam.shape[0], rows=d["own_rows"], mask=d["own_mask"]),
            d["face"],
            two_sided=d["two_sided"],
            total=d["total"],
            other=d["other"],
            flux_witness=RowTable(n, 2, rows=d["witness"], mask=d["witness_mask"]),
            flux_rows=d["flux_rows"],
            flux_index=d["flux_index"],
        )

    def as_kernel(self) -> dict:
        """The lane as the pass and the solve read it."""
        return dict(
            field=LANES.index(self.field),
            face=np.ascontiguousarray(self.face, bool),
            two_sided=np.ascontiguousarray(self.two_sided, bool),
            empty=np.ascontiguousarray(self.empty, bool),
            own_rows=np.ascontiguousarray(self.own_level.rows, np.float64),
            own_mask=np.ascontiguousarray(self.own_level.mask, bool),
            count=np.ascontiguousarray(self.count, np.float64),
            a=np.ascontiguousarray(self.a, np.float64),
            other=None if self.other is None else np.ascontiguousarray(self.other, np.float64),
            total=np.ascontiguousarray(self.total, np.float64),
            rho_ref=float(self.rho_ref),
            flux_rows=np.ascontiguousarray(self.flux_rows, np.float64),
            flux_index=np.ascontiguousarray(self.flux_index, np.int32),
            witness=np.ascontiguousarray(self.flux_witness.rows, np.float64),
            witness_mask=np.ascontiguousarray(self.flux_witness.mask, bool),
        )

    def serves(self, s: int, x: int) -> bool:
        """Does the lane carry its level across the face into ``x`` from ``s``?"""
        return bool(self.face[int(x), side_of(int(s), int(x))])

    def flux_at(self, x: int, side: int):
        """The junction flux level ``x`` holds at its ``side`` (0: from its left junction, 1: its right),
        or ``None``."""
        r = int(self.flux_index[int(x), int(side)])
        return None if r < 0 else self.flux_rows[r]


@dataclass(slots=True)
class Levels:
    """One population's LEVEL lane as RECEIVED: row ``i`` is the level node ``i`` holds on this lane from
    one side after a pass, or nothing (``present[i]`` False); ``count`` and ``opportunity`` the witness of
    the last full node the claim passed through; the RNA witness where ``has_witness``."""

    present: np.ndarray
    profile: np.ndarray
    count: np.ndarray
    opportunity: np.ndarray
    has_witness: np.ndarray
    rna_count: np.ndarray
    rna_count_var: np.ndarray

    @classmethod
    def empty(cls, n: int, K: int) -> "Levels":
        return cls(
            np.zeros(n, bool),
            np.empty((n, K)),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n, bool),
            np.full(n, np.nan),
            np.full(n, np.nan),
        )

    def write(self, i, profile, count, opportunity, rna_count=None, rna_count_var=None) -> None:
        i = int(i)
        self.present[i] = True
        self.profile[i] = profile
        self.count[i], self.opportunity[i] = float(count), float(opportunity)
        self.has_witness[i] = rna_count is not None
        self.rna_count[i] = np.nan if rna_count is None else float(rna_count)
        self.rna_count_var[i] = np.nan if rna_count_var is None else float(rna_count_var)

    def as_kernel(self) -> dict:
        return {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}


@dataclass(slots=True)
class Received:
    """What every node RECEIVED from one side after one pass, as a table: the backbone's neighbour bit,
    the composition where ``has_composition``, and the three level lanes. The two states where a node
    heard nothing: SILENCE (a neighbour, nothing present) and NO NEIGHBOUR (an open side)."""

    has_neighbour: np.ndarray
    has_composition: np.ndarray
    composition: np.ndarray
    level_gdna: Levels
    level_rna_pos: Levels
    level_rna_neg: Levels

    LANES = LANES

    @classmethod
    def empty(cls, n: int, K: int) -> "Received":
        return cls(
            np.zeros(n, bool),
            np.zeros(n, bool),
            np.empty((n, K)),
            Levels.empty(n, K),
            Levels.empty(n, K),
            Levels.empty(n, K),
        )

    @classmethod
    def from_kernel(cls, d: dict) -> "Received":
        """The kernel's received table (a sweep capture's ``from_left`` / ``from_right``) as a table."""
        return cls(
            np.asarray(d["has_neighbour"], bool),
            np.asarray(d["has_composition"], bool),
            np.asarray(d["composition"]),
            *(Levels(**{k: np.asarray(v) for k, v in d[lane].items()}) for lane in LANES),
        )

    def as_kernel(self) -> dict:
        return dict(
            has_neighbour=self.has_neighbour,
            has_composition=self.has_composition,
            composition=self.composition,
            **{lane: getattr(self, lane).as_kernel() for lane in LANES},
        )

    @property
    def has_level(self) -> np.ndarray:
        return self.level_gdna.present | self.level_rna_pos.present | self.level_rna_neg.present

    @property
    def heard(self) -> np.ndarray:
        return self.has_composition | self.has_level

    @property
    def silence(self) -> np.ndarray:
        return self.has_neighbour & ~self.heard

    @property
    def no_neighbour(self) -> np.ndarray:
        return ~self.has_neighbour


class Delivered(NamedTuple):
    """What the solve hands ψ: the fused λ rows (``None`` when nothing fused) and the cube table."""

    lam_rows: "np.ndarray | None"
    cube_rows: "CubeRows | None"


class Prepared:
    """One block's tables as the kernel built them — every node's own claim (`RowTable`), the rules per
    directed face (`Faces`), the lanes by population — with the pass and the solve run on them through the
    bindings. Built by `_prepared` from a context, or by hand (`hand_built`)."""

    def __init__(self, tables: dict, left, right, free_pos, free_neg):
        self.tables = tables
        lam = np.asarray(tables["lam"], np.float64)
        n = int(tables["own"]["mask"].shape[0])
        self.own = RowTable(n, lam.shape[0], rows=tables["own"]["rows"], mask=tables["own"]["mask"])
        self.faces = Faces(lam, left, right, tables["faces"])
        self.lanes = {
            name: LevelLane.from_kernel(name, d, lam) for name, d in tables["lanes"].items()
        }
        self.free_pos = np.ascontiguousarray(free_pos, bool)
        self.free_neg = np.ascontiguousarray(free_neg, bool)
        self.ambig = self.free_pos & self.free_neg

    @classmethod
    def hand_built(cls, own: RowTable, faces: Faces, lanes: dict, free_pos, free_neg) -> "Prepared":
        n = len(own)
        tables = dict(
            lam=faces.lam,
            n_u=np.zeros(n),
            own=dict(rows=own.rows, mask=own.mask),
            faces=faces.tables,
            lanes={name: ln.as_kernel() for name, ln in lanes.items()},
        )
        p = cls(tables, faces.nbr[:, 0], faces.nbr[:, 1], free_pos, free_neg)
        p.lanes = dict(lanes)  # the hand-built lanes themselves, so a gate reads what it wrote
        return p

    def _tables(self) -> dict:
        return dict(self.tables, lanes={name: ln.as_kernel() for name, ln in self.lanes.items()})

    def run_pass(self, received: Received, seq, nbr, terminal, *, backward: bool) -> None:
        """PHASE 1 on the table in one native call: for every destination in ``seq`` (chain order) whose
        neighbour ``nbr[i] >= 0`` and which is not a terminal, the recipient receives what its neighbour
        sends. ``backward`` names the pass; the tables say the rest."""
        transfer_pass(
            self._tables(),
            received.as_kernel(),
            np.ascontiguousarray(seq, np.int64),
            np.ascontiguousarray(nbr, np.int64),
            np.ascontiguousarray(terminal, bool),
        )

    def solve(self, from_left: Received, from_right: Received) -> Delivered:
        """PHASE 2, the policy's half: the two held tables into ψ's two channels."""
        live, rows, cube = transfer_solve(
            self._tables(),
            from_left.as_kernel(),
            from_right.as_kernel(),
            self.free_pos,
            self.free_neg,
        )
        return Delivered(rows if live else None, None if cube is None else CubeRows(**cube))


# ── the captured toy sweep ─────────────────────────────────────────────────────────────────────────────


def capture_sweep_inputs(tmp_path_factory):
    """ONE real `solve_chain` call captured from a calibrate run on the toy — the backbone-parity
    pattern: every gate re-runs the sweep with a different policy on byte-identical inputs. Each gate
    file wraps this in its own module-scoped ``sweep_inputs`` fixture."""
    from _prior_toy import build_toy

    toy = build_toy(tmp_path_factory)

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )
    from rigel.config import CalibrationConfig, PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

    index = toy.index
    scan_cfg = dataclasses.replace(
        PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(str(toy.bam_path))
    )
    _stats, strand_model, _buf, payload = scan_and_buffer(str(toy.bam_path), index, scan_cfg)
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    grabbed: list = []
    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    orig = SW.solve_chain

    def spy(chain, statics, geometry, belief, region_arrays, **kw):
        if not grabbed:
            grabbed.append((chain, statics, geometry, belief, region_arrays, dict(kw)))
        return orig(chain, statics, geometry, belief, region_arrays, **kw)

    calibrate_mod.solve_chain = spy
    try:
        calibrate_mod.calibrate(
            payload=payload,
            config=CalibrationConfig(),
            region_arrays=ra,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf,
            rna_fl_pmf=fl.rna_pmf,
            sj=build_sj_geometry_arrays(index),
            boundary_flags=build_boundary_flags_array(index),
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert grabbed, "the spy never fired"
    chain, statics, geometry, belief, region_arrays, kw = grabbed[0]
    kw = {k: v for k, v in kw.items() if k not in ("policy", "_capture")}
    return dict(
        args=(chain, statics, geometry, belief, region_arrays),
        kw=kw,
        payload=payload,
        calibrate_kw=dict(
            region_arrays=ra,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf,
            rna_fl_pmf=fl.rna_pmf,
            sj=build_sj_geometry_arrays(index),
            boundary_flags=build_boundary_flags_array(index),
        ),
    )


def _run(si, policy):
    out = SW.solve_chain(*si["args"], **si["kw"], policy=policy)
    return {f: np.asarray(getattr(out, f)) for f in ("f_g", "f_pos", "f_neg", "var_gdna")}


def _expected_pairs(si):
    """The intron|exon pairs derived INDEPENDENTLY of the policy (its falsification power):
    a BOUNDARY whose one flank is an exon REGION and whose other flank is an intron REGION
    that admits RNA."""
    from rigel.calibration.signature import coarse_type_array

    chain, statics, _geometry, _belief, region_arrays = si["args"]
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_exon = is_reg & (rtype[np.clip(obj, 0, rtype.shape[0] - 1)] == 2)
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    is_intron = is_reg & ~is_exon & (fp | fn)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    pairs = []
    for i in np.flatnonzero(~is_reg):
        lo, hi = left[i], right[i]
        if lo < 0 or hi < 0:
            continue
        if is_exon[lo] and is_intron[hi]:
            pairs.append((int(i), int(hi)))
        elif is_exon[hi] and is_intron[lo]:
            pairs.append((int(i), int(lo)))
    return pairs, int(chain.n_slots)


def _live_rows(si, n_grid, window):
    """Synthetic factory rows: a distinct non-flat row at every intron REGION slot, on the sweep's
    grid — what the live toy's context carries as ``factory_rows``."""
    from rigel.calibration.simplex_logodds import _logodds_grid

    pairs, n_slots = _expected_pairs(si)
    lam, _ = _logodds_grid(n_grid, window)
    rows = np.zeros((n_slots, lam.shape[0]))
    for _b, j in pairs:
        rows[j] = -0.05 * (lam - (0.1 * (j % 7) - 0.3)) ** 2  # non-flat, slot-distinct
    return rows


def _bits(n, pairs):
    """A lane's face table from directed ``(source, destination)`` pairs: ``(n, 2)`` bits over
    (destination, side) — the form `LevelLane` holds its faces in."""
    out = np.zeros((int(n), 2), bool)
    for s, i in pairs:
        out[int(i), side_of(int(s), int(i))] = True
    return out


def _pairs(bits, ctx):
    """The directed ``(source, destination)`` pairs a ``(n, 2)`` face table holds, read back through the
    chain's two neighbour arrays."""
    nbr = np.stack((np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)), axis=1)
    return {(int(nbr[i, sd]), int(i)) for i, sd in zip(*np.nonzero(np.asarray(bits, bool)))}


def _two_sided(lane, s, x) -> bool:
    return bool(lane.two_sided[int(x), side_of(int(s), int(x))])


def _prepared(pol, ctx, library=None) -> Prepared:
    """The kernel's builders on ``ctx`` — the whole chain as one block, in every gate here — with the
    library reduced over that same context unless one is given, exactly as the backbone pairs the two."""
    lib = pol.library(ctx) if library is None else library
    strand = pol.strand
    tables = transfer_prepare(
        lam=np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), int(ctx.n_grid)),
        is_boundary=np.ascontiguousarray(ctx.is_boundary, bool),
        is_exon=np.ascontiguousarray(ctx.is_exon_region, bool),
        free_pos=np.ascontiguousarray(ctx.free_pos, bool),
        free_neg=np.ascontiguousarray(ctx.free_neg, bool),
        exon_pos=np.ascontiguousarray(ctx.exon_pos, bool),
        exon_neg=np.ascontiguousarray(ctx.exon_neg, bool),
        left=np.ascontiguousarray(ctx.left, np.int64),
        right=np.ascontiguousarray(ctx.right, np.int64),
        flags=np.ascontiguousarray(ctx.boundary_flags, np.uint16),
        cnt=np.ascontiguousarray(ctx.unspliced_count, np.float64),
        spliced=np.ascontiguousarray(ctx.spliced_count, np.float64),
        sj_count=np.ascontiguousarray(ctx.sj_count, np.float64),
        sj_count_lo=np.ascontiguousarray(ctx.sj_count_lo, np.float64),
        sj_count_hi=np.ascontiguousarray(ctx.sj_count_hi, np.float64),
        route_rate_lo=np.ascontiguousarray(ctx.route_rate_lo, np.float64),
        route_rate_hi=np.ascontiguousarray(ctx.route_rate_hi, np.float64),
        eff_gdna=np.ascontiguousarray(ctx.eff_gdna, np.float64),
        eff_rna=np.ascontiguousarray(ctx.eff_rna, np.float64),
        belief_fg=np.ascontiguousarray(ctx.belief_fg, np.float64),
        has_own_composition=np.ascontiguousarray(ctx.has_own_composition, bool),
        factory_rows=None
        if ctx.factory_rows is None
        else np.ascontiguousarray(ctx.factory_rows, np.float64),
        has_strand=strand is not None,
        kappa=0.0 if strand is None else float(strand[0]),
        od_g=0.0 if strand is None else float(strand[1]),
        od_r=0.0 if strand is None else float(strand[2]),
        rho_gdna=0.0 if lib is None else float(lib.rho_gdna),
        rho_rna=0.0 if lib is None else float(lib.rho_rna),
        split_live=False if lib is None else bool(lib.split_live),
    )
    return Prepared(tables, ctx.left, ctx.right, ctx.free_pos, ctx.free_neg)


def _ctx_of(si) -> BlockContext:
    """The context exactly as the backbone builds it for the kernel — the chain's view, the incoming belief,
    every node's own-evidence bit read off a captured sweep — with the live toy's synthetic factory rows
    attached, so every gate's builders read the same rows the independent recomputes read
    (``ctx.factory_rows``)."""
    chain, statics, geometry, belief, ra = si["args"]
    kw = si["kw"]
    cap = SweepCapture()
    SW.solve_chain(*si["args"], **kw, policy=SilentPolicy(), _capture=cap)
    disc = strand_discriminability(float(kw["rna_sense_frac"]), float(kw.get("n_rna_obs", 0.0)))
    structure = SW._structure(chain, statics, ra)
    rows = _live_rows(si, int(kw["n_grid"]), float(kw["logodds_window"]))
    return BlockContext(
        **view_fields(chain, statics, geometry, structure),
        n_grid=int(kw["n_grid"]),
        logodds_window=float(kw["logodds_window"]),
        strand_live=disc > 0.0,
        has_own_composition=np.asarray(cap.tau_lam, np.float64) > 0.0,
        belief_fg=np.asarray(belief.f_g, np.float64),
        factory_rows=rows,
    )


def _strand_of(si):
    kw = si["kw"]
    return (
        float(kw["rna_sense_frac"]),
        float(kw.get("gdna_strand_overdispersion", 0.0)),
        float(kw.get("rna_strand_overdispersion", 0.0)),
    )


def _intron_mask(ctx):
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    fp, fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
    return ~is_bnd & ~is_exon & (fp | fn)


def _dead_boundaries(ctx):
    """The context with every BOUNDARY's strand channel declared dead and every region's intact."""
    live = np.asarray(ctx.has_own_composition, bool).copy()
    live[np.asarray(ctx.is_boundary, bool)] = False
    return dataclasses.replace(ctx, has_own_composition=live)


def _hop(prepared: Prepared, s, i, received=None) -> Received:
    """ONE hop of the native pass — the recipient ``i`` receives from ``s`` (its neighbour on that side):
    the face's composition rule on what ``s`` sends, and every lane face into ``i`` — on ``received``, a
    fresh table when none is given, so what ``s`` holds is whatever the caller wrote into row ``s``.
    Returns the table."""
    n, K = len(prepared.own), prepared.own.rows.shape[1]
    if received is None:
        received = Received.empty(n, K)
    nbr = np.full(n, -1, np.int64)
    nbr[int(i)] = int(s)
    received.has_neighbour[int(i)] = True
    prepared.run_pass(
        received, np.array([int(i)], np.int64), nbr, np.zeros(n, bool), backward=int(s) > int(i)
    )
    return received


_KEEP = object()


def _rule(prepared: Prepared, s, i, *, own=_KEEP, held=None):
    """The composition rule at the face into ``i`` from ``s``, applied to what ``s`` sends: its own claim
    — or ``own``, a substituted row, or ``None`` for no claim — composed with ``held``, what it holds from
    its far side (a row, or nothing). One hop of the native pass on a fresh table; the row ``i``
    receives, max-normalised, or ``None`` for no claim."""
    table = prepared.own
    kept_mask, kept_row = bool(table.mask[s]), table.rows[s].copy()
    if own is not _KEEP:
        table[s] = own
    received = Received.empty(len(table), table.rows.shape[1])
    if held is not None:
        received.composition[s] = held
        received.has_composition[s] = True
    try:
        _hop(prepared, s, i, received)
    finally:
        table.rows[s], table.mask[s] = kept_row, kept_mask
    return received.composition[int(i)].copy() if received.has_composition[int(i)] else None


def _lane_prepared(lane: LevelLane, left, right) -> Prepared:
    """A hand-built lane as the pass reads it: tables with no claims and no composition rules, that one
    lane, and the chain's two neighbour arrays."""
    n, K = len(lane.own_level), lane.own_level.rows.shape[1]
    return Prepared.hand_built(
        RowTable(n, K),
        Faces(lane.lam, left, right),
        {lane.population: lane},
        np.ones(n, bool),
        np.zeros(n, bool),
    )


def norm(row):
    """A row max-normalised (the oracle)."""
    row = np.asarray(row, np.float64)
    return row - row.max()


def fuse(parts):
    """Independent witnesses about one slot: log-profiles add, then re-normalise (the oracle)."""
    return norm(sum(np.asarray(q, np.float64) for q in parts))


def intersect(bounds):
    """Two or more bounds on ONE density combine by INTERSECTION: the pointwise minimum of their
    log-profiles, max-normalised — bounds intersect, they do not multiply (the oracle)."""
    out = np.asarray(bounds[0], np.float64)
    for b in bounds[1:]:
        out = np.minimum(out, np.asarray(b, np.float64))
    return norm(out)


def _native_passes(prepared: Prepared, ctx):
    """The backbone's two directional passes through the kernel on fresh tables; no terminal. Returns
    the two tables, ``(from_left, from_right)``."""
    order = np.arange(int(ctx.n_slots), dtype=np.int64)
    tables = []
    for nbr, seq, backward in ((ctx.left, order, False), (ctx.right, order[::-1], True)):
        nbr = np.asarray(nbr, np.int64)
        received = Received.empty(order.size, int(ctx.n_grid))
        received.has_neighbour[seq] = nbr[seq] >= 0
        prepared.run_pass(received, seq, nbr, np.zeros(order.size, bool), backward=backward)
        tables.append(received)
    return tuple(tables)


#: a row matrix of a `Received` table and the presence bits that say which of its rows are rows at all
_ROWS_OF = {"composition": "has_composition", "profile": "present"}


def _leaves(a, b, prefix=""):
    """Every array of two `Received` tables side by side, as ``(name, x, y)``, a row matrix restricted
    to the rows ``a``'s bits say are present — the matrices are allocated unfilled, so an absent row is
    not compared; the bits themselves are leaves and are compared whole."""
    for f in dataclasses.fields(a):
        x, y = getattr(a, f.name), getattr(b, f.name)
        if dataclasses.is_dataclass(x):
            yield from _leaves(x, y, prefix + f.name + ".")
            continue
        x, y = np.asarray(x), np.asarray(y)
        if f.name in _ROWS_OF:
            keep = np.asarray(getattr(a, _ROWS_OF[f.name]), bool)
            x, y = x[keep], y[keep]
        yield prefix + f.name, x, y


def _drive(prepared: Prepared, ctx):
    """The backbone's own contract, reproduced: the two passes and the solve, which receives the two
    tables. Returns ``(rows, from_left, from_right)`` — the delivered rows (zeros when the policy is
    silent) and the tables."""
    from_left, from_right = _native_passes(prepared, ctx)
    msg = prepared.solve(from_left, from_right)
    rows = (
        np.zeros((int(ctx.n_slots), int(ctx.n_grid)))
        if msg.lam_rows is None
        else np.asarray(msg.lam_rows)
    )
    return rows, from_left, from_right


def _drive_the_backbone(prepared: Prepared, ctx):
    """`_drive`'s rows alone."""
    return _drive(prepared, ctx)[0]


def _strand_row_of(ctx, strand, lam, x):
    """A slot's own strand profile, recomputed independently (the frozen-variance count form)."""
    kappa, od_g, od_r = strand
    fp = np.asarray(ctx.free_pos, bool)
    cnt = np.asarray(ctx.unspliced_count, np.float64)
    belief = np.asarray(ctx.belief_fg, np.float64)
    fg = 1.0 / (1.0 + np.exp(-lam))
    n = cnt[x].sum()
    ks = kappa if fp[x] else 1.0 - kappa
    f_ref = float(np.clip(belief[x], 1e-9, 1 - 1e-9))
    p = 0.5 * fg + ks * (1 - fg)
    p_ref = 0.5 * f_ref + ks * (1 - f_ref)
    var = max(
        n * p_ref * (1 - p_ref)
        + (n * f_ref) ** 2 * 0.25 * od_g
        + (n * (1 - f_ref)) ** 2 * ks * (1 - ks) * od_r,
        1e-9,
    )
    row = -0.5 * (cnt[x, 0] - n * p) ** 2 / var
    return row - row.max()


def _full_policy(sweep_inputs):
    """The shipped policy with the toy's strand model, the live rows it will find on `_ctx_of`'s
    context, and the grid: ``(policy, rows, n_grid, window)``."""
    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    rows = _live_rows(sweep_inputs, n_grid, window)
    return TransferPolicy(strand=_strand_of(sweep_inputs)), rows, n_grid, window


def _rna_lanes_of(sweep_inputs, ctx=None):
    pol = _full_policy(sweep_inputs)[0]
    ctx = _ctx_of(sweep_inputs) if ctx is None else ctx
    prepared = _prepared(pol, ctx)
    assert set(prepared.lanes) == {"gdna", "pos", "neg"}
    return ctx, prepared


def _rna(prepared) -> dict:
    """The two RNA lanes by strand."""
    return {name: prepared.lanes[name] for name in ("pos", "neg")}


def _strand_intron(ctx, name):
    """The PER-STRAND intron test the lane uses: a region that admits ``s`` and carries no exon of ``s``."""
    is_bnd = np.asarray(ctx.is_boundary, bool)
    free = np.asarray(ctx.free_pos if name == "pos" else ctx.free_neg, bool)
    exon_s = np.asarray(ctx.exon_pos if name == "pos" else ctx.exon_neg, bool)
    return ~is_bnd & free & ~exon_s


def _held_rna(held, x, field):
    """Row ``x`` of the lane ``field`` of the table ``held``, or ``None`` where nothing is present."""
    lv = getattr(held, field)
    return lv.profile[int(x)].copy() if lv.present[int(x)] else None


def _with_populated_inside(ctx):
    """The toy's inside pieces are EMPTY (shorter than a fragment: no contained fragment, a contained
    opportunity below one base), so the level rule serves nothing there and a gate on the bare toy proves
    nothing. A context is data: populate the two inside slots with consistent counts (both strand
    channels live, so the step's spread is fitted from two pairs) and a plausible contained
    opportunity, for the policy and the independent recompute alike."""
    R = transfer_rows
    is_bnd = np.asarray(ctx.is_boundary, bool)
    is_exon = np.asarray(ctx.is_exon_region, bool)
    left, right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
    flags = np.asarray(ctx.boundary_flags, np.uint16)
    cnt = np.asarray(ctx.unspliced_count, np.float64).copy()
    a_g = np.asarray(ctx.eff_gdna, np.float64).copy()
    a_r = np.asarray(ctx.eff_rna, np.float64).copy()
    live = np.asarray(ctx.has_own_composition, bool).copy()
    fills = iter(
        [(20.0, 180.0), (80.0, 120.0), (25.0, 175.0), (70.0, 130.0)]
    )  # two pairs that DISAGREE
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]):
            continue
        _o, i = R.outside_flank(int(flags[b]), int(lo), int(hi))
        if i is None:
            continue
        cnt[i] = next(fills)
        a_g[i] = a_r[i] = 150.0
        live[i] = True
    return dataclasses.replace(
        ctx,
        unspliced_count=cnt,
        eff_gdna=a_g,
        eff_rna=a_r,
        has_own_composition=live,
    )


def _with_alt_splice_sites(ctx):
    """The toy carries no alternative splice site: turn its two exon|exon terminus boundaries into a
    DONOR (intron to the right) and an ACCEPTOR (intron to the left) with a route flux each, and
    populate the pieces beyond them — a context is data; the policy and the recompute read the same."""
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_NEG, FLAG_DONOR_POS, FLAG_TERMINUS

    base = _with_populated_inside(ctx)
    is_bnd = np.asarray(base.is_boundary, bool)
    is_exon = np.asarray(base.is_exon_region, bool)
    left, right = np.asarray(base.left, np.int64), np.asarray(base.right, np.int64)
    flags = np.asarray(base.boundary_flags, np.uint16).copy()
    sjc = np.asarray(base.sj_count, np.float64).copy()
    kinds = iter([(FLAG_DONOR_POS, (40.0, 0.0)), (FLAG_ACCEPTOR_NEG, (0.0, 50.0))])
    for b in np.flatnonzero(is_bnd):
        lo, hi = left[b], right[b]
        if lo < 0 or hi < 0 or not (is_exon[lo] and is_exon[hi]) or not (flags[b] & FLAG_TERMINUS):
            continue
        kind, fl = next(kinds)
        flags[b] = kind
        sjc[b] = fl
    return dataclasses.replace(base, boundary_flags=flags, sj_count=sjc)
