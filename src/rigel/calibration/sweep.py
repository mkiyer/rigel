"""The backbone of the calibration solve: the self-solve, two directional passes, one ψ solve and one write-back.

Each slot's unspliced fragment mass is deconvolved into ``(f_pos, f_neg, f_g)`` (sense RNA, antisense RNA, gDNA) over
the ``N E N E … N`` chain (``region_chain``). Every node states its own claim; a forward then a backward pass lets
each recipient receive what its neighbour sends; and each node is solved once from its own evidence, the two held
messages and the prior. The chain is a forest of linear paths, so this is exact belief propagation, not an iteration.

The sweep is one native call (``native.solve_blocks``, ``native/solve_kernel.cpp``). This module cuts the chain into
locus blocks, reduces the policy's library over the whole chain, hands the kernel the chain's arrays, the blocks and
the landscape's curve, and reads back the belief, ``has_composition`` and the assertions' counts. The kernel solves
the blocks on a thread pool and reduces only integer counts across them, so the answer is bit-identical at every
thread count.

What a message says belongs to a :mod:`~.messages` policy, built in the kernel; this module must not name what a
message is about. It owns the shape of the solve and two invariants no policy may break:

- every delivered row is one finite row per slot, which catches a NaN reaching ψ (counted per block by the kernel and
  judged in :class:`AssertionCounts`);
- the kernel sees only the two neighbour states: each destination's row is built from the source's claim and what the
  source holds, never from the destination's own belief (enforced by construction).

A terminal, a region that admits no RNA strand and so is structurally pure gDNA, receives nothing. Nothing crosses a
terminal, so the chain splits into independent loci there (``region_chain.locus_blocks``).

Both gates on a slot read the region signature, never the counts. The solve gate (``solvable``): a slot is solved iff
it admits at least one RNA strand and has unspliced mass. A slot with no admissible RNA strand is a locked all-gDNA
object that keeps its signature-binary init. The emission gate, which messages a slot sends, is the policy's.

Gate: ``tests/calibration/test_sweep_backbone.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..native import POLICY_SILENT, POLICY_TRANSFER, solve_blocks
from .diagnostics import SweepDiagnostics
from .messages import ChainView, view_fields
from .region_geometry import (
    RegionBelief,
    RegionGeometry,
    RegionStatics,
    region_gdna_geometry,
    slot_signature,
)
from .strand_balance import strand_discriminability
from .signature import BIT_EXON_NEG, BIT_EXON_POS, RegionType, coarse_type_array
from .simplex_logodds import TILT_NODES, CubeRows, logodds_grid
from .region_chain import (
    BOUNDARY,
    REGION,
    RegionChain,
    RegionDeconv,
    locus_blocks,
    on_object_axis,
)

__all__ = [
    "AssertionCounts",
    "chain_boundary_deconv",
    "chain_region_deconv",
    "solve_chain",
]


#: the kernel's code for each policy name (`native/solve_kernel.cpp`): the silent policy runs no layer
POLICY_KERNEL = {"silent": POLICY_SILENT, "transfer": POLICY_TRANSFER}

#: the two assertions that run only where a channel was delivered — the kernel says per block whether it was
_DELIVERED = {"lam_rows_finite": "rows_delivered", "cube_rows_finite": "cube_delivered"}


class AssertionCounts(dict):
    """How many slots violated each backbone assertion, published into the diagnostics capture.

    Each count carries the number of eligible slots beside it, because zero violations where the predicate can never
    fire is not evidence.
    """

    def note(self, name: str, violations: int, eligible: int) -> None:
        n_v, n_e = int(violations), int(eligible)
        self._add(name, n_v, n_e)
        if n_v:
            raise AssertionError(
                f"backbone assertion {name!r} violated at {n_v:,} of {n_e:,} eligible slots. "
                f"The policy may not do this — see rigel.calibration.sweep's module docstring for what "
                f"each assertion catches; never widen the predicate."
            )

    def _add(self, name: str, n_v: int, n_e: int) -> None:
        cur = self.get(name, {"violations": 0, "eligible": 0})
        self[name] = {"violations": cur["violations"] + n_v, "eligible": cur["eligible"] + n_e}

    @classmethod
    def of_blocks(cls, names, counts, delivered: dict, n_owned) -> "AssertionCounts":
        """The kernel's per-block ``(violations, eligible)`` summed over the chain, each assertion noted once (so a
        violation raises).

        An assertion runs only where its channel was delivered and is absent when no block delivered it. For the λ-row
        check, a block that delivered no rows counts its owned slots as eligible (zero rows are finite rows), so the
        published count does not depend on how the chain was cut.
        """
        out = cls()
        counts = np.asarray(counts, np.int64)
        n_owned = np.asarray(n_owned, np.int64)
        for a, name in enumerate(names):
            flag = _DELIVERED[name]
            if not delivered[flag].any():
                continue
            v, e = int(counts[:, a, 0].sum()), int(counts[:, a, 1].sum())
            if name == "lam_rows_finite":
                e += int(n_owned[~delivered[flag]].sum())
            out.note(name, v, e)
        return out


def solve_chain(
    chain: RegionChain,
    statics: RegionStatics,
    geometry: RegionGeometry,
    belief: RegionBelief,
    region_arrays,
    *,
    rna_sense_frac: float,
    n_rna_obs: float = 0.0,
    n_grid: int,
    logodds_window: float,
    gdna_prior=None,
    policy,
    block_slots: int | None = None,
    n_threads: int = 1,
    diagnostics: SweepDiagnostics | None = None,
) -> RegionBelief:
    """One forward-backward sweep over the chain; returns the resolved :class:`RegionBelief`.

    The chain breaks into loci at every terminal, and ``region_chain.locus_blocks`` merges consecutive loci into blocks
    of up to ``block_slots`` slots (``None``: the whole chain as one block). Each block is solved on its own slice of
    every input, reading one slot past its range where that slot is the terminal its last node receives from. Only the
    policy's library (``policy.library``, reduced once over the whole chain) crosses a block boundary, and ψ's read-out
    is chunk-exact (``simplex_logodds``), so the answer is the same for every ``block_slots``; the block size trades
    working set against per-block overhead.

    ``policy`` is the message policy (:mod:`~.messages`): a name the kernel switches on, a strand model and a library
    reduction. It has no default; ``calibrate`` installs the one the config names,
    :class:`~.messages.transfer.TransferPolicy` or :class:`~.messages.silent.SilentPolicy`.

    ``n_threads`` is the thread budget (``CalibrationConfig.n_threads``; 0 is every core).

    ``gdna_prior`` is the fitted landscape (``landscape.DensityLandscape``), or ``None`` for a prior-free solve in which
    ψ carries the Jeffreys reference measure alone on both arms. The prior-free pass only trains the population gDNA
    hyperprior, so it need not answer objects it cannot solve. The kernel reads the curve at every cell on each slot's
    own gDNA support (``region_geometry.region_gdna_geometry``).

    ``diagnostics``, a record instruments and tests pass, asks the kernel for its diagnostic arrays and is filled with
    them (:class:`~.diagnostics.SweepDiagnostics`); production passes none.
    """
    kappa = float(rna_sense_frac)
    disc = strand_discriminability(kappa, n_rna_obs)
    # the library: the policy's reductions over the whole chain, from observations and geometry alone; the only
    # information a message may carry across a locus boundary
    view = ChainView(
        **view_fields(chain, statics, geometry, slot_structure(chain, region_arrays)),
        strand_live=disc > 0.0,
    )
    library = policy.library(view)
    lam, _ = logodds_grid(int(n_grid), float(logodds_window))
    terminal = view.terminal
    blocks = locus_blocks(chain, terminal, block_slots)
    out = {
        k: np.ascontiguousarray(getattr(belief, k), dtype=np.float64).copy()
        for k in ("f_pos", "f_neg", "f_g", "var_gdna")
    }
    has_composition = np.zeros(int(chain.n_slots), dtype=bool)
    buffers = None if diagnostics is None else _diagnostic_buffers(chain.n_slots, lam.shape[0])
    res = solve_blocks(
        **_chain_arrays(view, terminal),
        belief_fpos=np.ascontiguousarray(belief.f_pos, np.float64),
        belief_fneg=np.ascontiguousarray(belief.f_neg, np.float64),
        belief_fg=np.ascontiguousarray(belief.f_g, np.float64),
        blocks=np.array([(b.start, b.stop, b.end) for b in blocks], np.int64).reshape(-1, 3),
        kappa=kappa,
        disc=float(disc),
        lam=lam,
        n_tilt=int(TILT_NODES),
        policy=POLICY_KERNEL[policy.name],
        **_library_args(library),
        gdna=_landscape_curve(gdna_prior),
        out_fpos=out["f_pos"],
        out_fneg=out["f_neg"],
        out_fg=out["f_g"],
        out_var=out["var_gdna"],
        out_has_composition=has_composition,
        diagnostics=buffers,
        n_threads=int(n_threads),
    )
    n_owned = np.array([b.stop - b.start for b in blocks], np.int64)
    delivered = {k: np.asarray(res[k], bool) for k in ("rows_delivered", "cube_delivered")}
    counts = AssertionCounts.of_blocks(res["assertions"], res["counts"], delivered, n_owned)
    if diagnostics is not None:
        _record_sweep(diagnostics, view, out, buffers, res, blocks, geometry, lam)
        diagnostics.backbone_assertions = counts
        diagnostics.policy_name = str(policy.name)
    return RegionBelief(**out, has_composition=has_composition)


def _library_args(library) -> dict:
    """The policy's library as the kernel takes it; the silent policy has none."""
    if library is None:
        return dict(rho_gdna=0.0, rho_rna=0.0, split_live=False)
    return dict(
        rho_gdna=float(library.rho_gdna),
        rho_rna=float(library.rho_rna),
        split_live=bool(library.split_live),
    )


def _landscape_curve(gdna_prior):
    """The landscape's curve ``(log_rho, logP)`` as the kernel takes it, or ``None`` for the prior-free solve."""
    if gdna_prior is None:
        return None
    return (
        np.ascontiguousarray(gdna_prior.log_rho, np.float64),
        np.ascontiguousarray(gdna_prior.logP, np.float64),
    )


def _chain_arrays(view: ChainView, terminal) -> dict:
    """The chain's arrays under the kernel's argument names, each C-contiguous in the kernel's dtype: the view's
    observations and geometry, and the terminal bits."""

    def f64(a):
        return np.ascontiguousarray(a, np.float64)

    def bits(a):
        return np.ascontiguousarray(a, bool)

    return dict(
        is_boundary=bits(view.is_boundary),
        is_exon=bits(view.is_exon_region),
        free_pos=bits(view.free_pos),
        free_neg=bits(view.free_neg),
        exon_pos=bits(view.exon_pos),
        exon_neg=bits(view.exon_neg),
        terminal=bits(terminal),
        left=np.ascontiguousarray(view.left, np.int64),
        right=np.ascontiguousarray(view.right, np.int64),
        flags=np.ascontiguousarray(view.boundary_flags, np.uint16),
        cnt=f64(view.unspliced_count),
        spliced=f64(view.spliced_count),
        sj_count=f64(view.sj_count),
        sj_count_lo=f64(view.sj_count_lo),
        sj_count_hi=f64(view.sj_count_hi),
        route_rate_lo=f64(view.route_rate_lo),
        route_rate_hi=f64(view.route_rate_hi),
        eff_gdna=f64(view.eff_gdna),
        eff_rna=f64(view.eff_rna),
    )


@dataclass(frozen=True, slots=True)
class SlotStructure:
    """The per-slot structural classes, from the signature: the exon region (the splice-in message's destination
    class, a policy input rather than a gate) and the signature's two exon bits, which tell a strand's intron from its
    exon inside an overlapping locus. The terminal is the view's own predicate (``ChainView.terminal``)."""

    is_exon_region: np.ndarray
    exon_pos: np.ndarray
    exon_neg: np.ndarray


def slot_structure(chain, region_arrays) -> SlotStructure:
    """The structural classes from each slot's signature; a boundary slot's signature is 0, so it is none of them."""
    sig = slot_signature(chain, region_arrays.signature)
    return SlotStructure(
        is_exon_region=coarse_type_array(sig) == RegionType.EXON,
        exon_pos=(sig & BIT_EXON_POS) > 0,
        exon_neg=(sig & BIT_EXON_NEG) > 0,
    )


def _diagnostic_buffers(n: int, n_grid: int) -> dict:
    """The per-slot arrays the kernel's diagnostics mode fills."""
    return {
        "fg_loc": np.zeros(n),
        "tau_lam": np.zeros(n),
        "lam_rows": np.zeros((n, n_grid)),
    }


def _record_sweep(diagnostics: SweepDiagnostics, view, out, buffers, res, blocks, geometry, lam):
    """Fill the sweep's diagnostics record from the kernel's per-slot arrays and its per-block cube rows and received
    tables."""
    d = diagnostics
    d.fg_loc = buffers["fg_loc"]
    d.tau_lam = buffers["tau_lam"]
    d.f_g = out["f_g"].copy()
    d.free_pos = np.asarray(view.free_pos, bool)
    d.free_neg = np.asarray(view.free_neg, bool)
    d.solvable = (d.free_pos | d.free_neg) & (view.n_slot > 0.0)
    d.count = np.asarray(view.unspliced_count, np.float64)
    d.mass_global, _ = region_gdna_geometry(geometry)
    d.lam_rows = buffers["lam_rows"] if np.asarray(res["rows_delivered"], bool).any() else None
    cubes = [
        _owned_cube_rows(cube, b, lam) for b, cube in zip(blocks, res["cubes"]) if cube is not None
    ]
    d.cube_rows = CubeRows.concat(cubes) if cubes else None
    received = [r for r in res["received"] if r is not None]
    d.from_left = _concat_received([r[0] for r in received]) if received else None
    d.from_right = _concat_received([r[1] for r in received]) if received else None


def _owned_cube_rows(cube, block, lam) -> CubeRows:
    """One block's cube rows at its own slots (not its read-ahead terminal), keyed to the chain."""
    slot, pos, has_pos, neg, has_neg, total, opportunity, rho = cube
    keep = np.asarray(slot) < (block.stop - block.start)
    return CubeRows(
        slot=np.asarray(slot)[keep] + block.start,
        profile_pos=np.asarray(pos)[keep],
        has_pos=np.asarray(has_pos)[keep],
        profile_neg=np.asarray(neg)[keep],
        has_neg=np.asarray(has_neg)[keep],
        u=lam,
        total=np.asarray(total)[keep],
        opportunity=np.asarray(opportunity)[keep],
        rho_ref=np.asarray(rho)[keep],
    )


def _concat_received(parts: list) -> dict:
    """The blocks' received tables as the chain's: each array concatenated over the blocks' owned slots in block order
    (a dict of arrays, its lanes dicts of arrays)."""
    out = {}
    for name, first in parts[0].items():
        if isinstance(first, dict):
            out[name] = _concat_received([p[name] for p in parts])
        else:
            out[name] = np.concatenate([np.asarray(p[name]) for p in parts], axis=0)
    return out


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# The output contract: the chain belief projected back onto the two payload axes for ``CalibrationResult``,
# ``priors`` and ``derive``. Not part of the solve.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def chain_region_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's region slots onto the region axis as a :class:`RegionDeconv`.

    The RNA mass has no spliced term: the accumulator credits ``region_contained`` only when the fragment used no sj,
    so a contained fragment is unspliced by construction.
    """
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    f_g, f_pos, f_neg = (
        on_object_axis(chain, f, REGION, count.shape[0])
        for f in (belief.f_g, belief.f_pos, belief.f_neg)
    )
    return RegionDeconv(
        gdna_mass=f_g * count,
        rna_mass=(1.0 - f_g) * count,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_boundary_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's boundary slots onto the contiguous-boundary axis as a :class:`RegionDeconv`.

    One result per boundary, not a ``(left, right)`` pair of per-region halves: splitting a boundary's flux onto its
    flanks only for ``priors`` to pool them again is a no-op that invites a lost factor of 2.

    The RNA mass includes the boundary's spliced crossings (``boundary_spliced``), which are RNA whatever the unspliced
    mixture resolves to, since gDNA cannot be spliced. ψ solves the full simplex at every slot (the population set
    depends on the ``free_*`` bits, never on the slot's kind), so a boundary carries a per-strand RNA split like a
    region.
    """
    unspliced = np.asarray(substrate.boundary_unspliced.count, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(axis=1)
    f_g, f_pos, f_neg = (
        on_object_axis(chain, f, BOUNDARY, unspliced.shape[0])
        for f in (belief.f_g, belief.f_pos, belief.f_neg)
    )
    return RegionDeconv(
        gdna_mass=f_g * unspliced,
        rna_mass=(1.0 - f_g) * unspliced + spliced,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )
