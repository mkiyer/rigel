"""THE BACKBONE — the self-solve, two directional passes, one ψ solve, one write-back, four assertions.

       Gate: ``tests/calibration/test_sweep_backbone.py``

Each slot's unspliced fragment mass is deconvolved into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the ``N E N E … N`` chain (`region_chain`), on the TWO-PHASE shape: every
node's own claim, a forward pass then a backward pass in which each RECIPIENT receives what its
neighbour sends, and ONE solve per node from its own evidence, the two held messages and the prior. The
chain is a forest of linear paths, so that is exact belief propagation, not an iteration.

THE SWEEP IS ONE NATIVE CALL (`native.solve_blocks`, `native/solve_kernel.cpp`). This file cuts the chain
into locus blocks, reduces the policy's library over the whole chain, asks the message cache which blocks
it can serve, hands the kernel the chain's arrays, the blocks, the priors' INPUTS (the landscape's curve;
the intron factory's background, mask, counts and opportunities) and the served deliveries, and reads
back the belief, ``has_composition``, the assertions' counts and — for the cache — each block's delivery.
Inside the call a pool of threads pulls the blocks one at a time; each block runs the whole pipeline on
its own thread and arena — the prior rows, the self-solve ψ, the own-evidence precision, the policy's
builders, the two passes, the solve, the final ψ, the write-back, the counts — and nothing is reduced
across blocks but integer counts, so the answer is BIT-IDENTICAL at every thread count.

This file knows nothing about capture, splice in, levels, lanes or enrichment — those words do not appear
in it, and that is the design rather than tidiness. Everything about *what a message says* is a
:mod:`~.messages` policy, built in the kernel. What is left here is the shape of the solve and the four
invariants no policy may break, counted in the kernel on every block's owned slots and judged here:

===================================================  ====================================================
the backbone asserts                                 it would have caught
===================================================  ====================================================
every delivered row is one row per slot, finite      a NaN reaching ψ (a row off the solve grid cannot be
                                                     built: the kernel's rows are the grid's)
``|T| <= 3``                                         AXIOM 0, made executable
the write-back touches only ``solvable`` slots       a replay that read the untouched mask as a difference
the kernel sees only the two NEIGHBOUR states        a message built from the destination's own belief —
                                                     enforced by construction: the pass builds each
                                                     destination's row from the source's claim and what
                                                     the source holds, and the only belief the layer reads
                                                     is the incoming one at a node's OWN claim
===================================================  ====================================================

The assertions' semantics live HERE (`AssertionCounts`), not in the policy, and that is the entire point:
a future policy can be as wrong as it likes and still cannot commit any of them, and each has shipped at
least once.

And ONE structural rule, which is what lets the chain be solved a locus at a time: a TERMINAL — a REGION
that admits no RNA strand, so it is structurally pure gDNA, solved and fixed before any message exists —
RECEIVES NOTHING. The pass never asks the kernel for a hop into one; the terminal holds silence from a side
it has a neighbour on. Nothing can therefore cross a terminal, and the chain breaks into independent loci
at every one (`region_chain.locus_blocks`). Measured on the human chain before it was made structural: of
1,206,202 composition faces and 4,621,302 lane faces the shipped policy built, none delivered into a
terminal, so the rule moved no number; it is written so that no future policy can move the boundary
condition a block solve depends on.

Two gates on this slot in the pipeline
--------------------------------------
Both come from the region SIGNATURE and never from the counts:

* the SOLVE gate (``solvable``) — a slot deconvolves its own split iff it admits >= 1 RNA strand and has
  unspliced mass. A slot with no admissible RNA strand (an intergenic region, or a gene boundary) is a
  LOCKED all-gDNA object; it is not solved and keeps its signature-binary init, because RNA cannot cross a
  gene boundary so its unspliced mass is purely gDNA.
* the EMISSION gate — which MESSAGES a slot sends. That is a policy question and lives there.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np

from ..native import solve_blocks
from .blocks import SweepCapture, view_fields
from .message_cache import MessageCache
from .messages import ChainView
from .messages.silent import SilentPolicy
from .region_geometry import (
    RegionBelief,
    RegionGeometry,
    RegionStatics,
    g1_locked,
    region_gdna_geometry,
)
from .region_init import strand_discriminability
from .signature import BIT_EXON_NEG, BIT_EXON_POS, coarse_type_array
from .simplex_logodds import _TILT_NODES, CubeRows, _logodds_grid
from .region_chain import BOUNDARY, REGION, RegionChain, RegionDeconv, locus_blocks

__all__ = [
    "AssertionCounts",
    "chain_boundary_deconv",
    "chain_region_deconv",
    "solve_chain",
]


# ⛔ ASSERTIONS A SHIPPED POLICY IS KNOWN TO VIOLATE, each entered with the measurement that proved it.
# An entry is COUNTED and PUBLISHED rather than raised, because widening an assertion to fit a defect is
# how a gate becomes vacuous; an entry also
# carries a STRICT xfail in the gate file, this project's convention for a PROVEN defect whose fix is
# panel-negative on its own. The dict is EMPTY, so anything violated raises.
#: ``name -> why it is not fatal yet``.
_KNOWN_VIOLATIONS: dict[str, str] = {}

#: the kernel's name for each policy (`native/solve_kernel.cpp`): the silent policy runs no layer
_POLICY_KERNEL = {"silent": 0, "transfer": 1}

#: the two assertions that run only where a channel was delivered — the kernel says per block whether it was
_DELIVERED = {"lam_rows_finite": "rows_delivered", "cube_rows_finite": "cube_delivered"}


class AssertionCounts(dict):
    """How many slots violated each backbone assertion, published into the diagnostics capture.

    A count rather than a bool: before believing "the arm changed nothing", check it COULD have changed
    something. An assertion reporting 0 violations on a substrate where the predicate can never fire is
    not evidence, so the report also carries how many slots were ELIGIBLE for each check.
    """

    def note(self, name: str, violations: int, eligible: int) -> None:
        n_v, n_e = int(violations), int(eligible)
        self._add(name, n_v, n_e)
        if n_v and name not in _KNOWN_VIOLATIONS:
            raise AssertionError(
                f"backbone assertion {name!r} violated at {n_v:,} of {n_e:,} eligible slots. "
                f"The policy may not do this — see rigel.calibration.sweep's module docstring for what "
                f"each assertion catches. If this is a NEW and PROVEN defect whose fix is panel-negative "
                f"alone, add it to _KNOWN_VIOLATIONS with the measurement, never widen the predicate."
            )

    def _add(self, name: str, n_v: int, n_e: int) -> None:
        cur = self.get(name, {"violations": 0, "eligible": 0})
        self[name] = {"violations": cur["violations"] + n_v, "eligible": cur["eligible"] + n_e}

    @classmethod
    def of_blocks(cls, names, counts, delivered: dict, n_owned) -> "AssertionCounts":
        """The kernel's per-block ``(violations, eligible)`` summed to the chain's, each assertion noted
        once (so a violation raises). An assertion that runs only where a channel was delivered is
        absent from the report when no block delivered it; the λ-row check's ELIGIBLE set reads as the
        chain's — where any block delivered rows, a block that delivered none holds zero rows, which are
        finite rows that were checked — so the published count does not depend on how the chain was
        cut."""
        out = cls()
        counts = np.asarray(counts, np.int64)
        n_owned = np.asarray(n_owned, np.int64)
        for a, name in enumerate(names):
            flag = _DELIVERED.get(name)
            if flag is not None and not delivered[flag].any():
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
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_rna_obs: float = 0.0,
    n_grid: int,
    logodds_window: float = 10.0,
    gdna_prior=None,
    intron_prior=None,
    policy=None,
    block_slots: int | None = None,
    message_cache: "MessageCache | None" = None,
    n_threads: int = 1,
    _capture: SweepCapture | None = None,
) -> RegionBelief:
    """One forward-backward sweep over the chain. Returns the resolved :class:`RegionBelief`.

    ``message_cache`` shares the message layer's output between sweeps whose inputs to it are identical
    (:class:`MessageCache`): a block whose digest is held is served its delivery and pays only its two ψ
    solves. ``None`` runs the layer for every block, as does a diagnostic capture.

    THE SWEEP IS SOLVED A LOCUS BLOCK AT A TIME, IN ONE NATIVE CALL. The chain breaks at every TERMINAL —
    a region that admits no RNA strand, structurally pure gDNA, which receives nothing (the module
    docstring) — into loci, and `region_chain.locus_blocks` merges consecutive loci into blocks of up to
    ``block_slots`` slots (``None``: the whole chain as one block). Each block is solved on its own slice
    of every input — its own claims, the policy's claims and rules, the two passes, the solve and the
    write-back — reading one slot beyond its own range where that slot is the terminal its last node
    receives from. The only information that crosses a block boundary is the policy's LIBRARY
    (:meth:`~.messages.Policy.library`), reduced once over the whole chain from observations and
    geometry, and the answer is the same for every ``block_slots`` because ψ's read-out is
    chunk-exact (`simplex_logodds`); the block size sets only the working set — a block's
    ``(slots, K)`` arrays instead of the chain's — against the per-block overhead.

    ``policy`` is the message-composition policy (:mod:`~.messages`): a NAME the kernel switches on, a
    strand model and a library reduction. It defaults to :class:`~.messages.silent.SilentPolicy`, which
    sends nothing. The shipped answer is :class:`~.messages.transfer.TransferPolicy`, which ``calibrate``
    passes explicitly.

    ``n_threads`` is the thread budget (`CalibrationConfig.n_threads`; 0 is every core): the blocks are
    pulled one at a time by a pool, bit-identically at every count.

    ``gdna_prior`` is the fitted landscape (`landscape.DensityLandscape`) or ``None`` — a first-class
    PRIOR-FREE solve: ψ then carries the Jeffreys reference measure alone on both arms. Prior-free is not
    reference-free. That pass's only job is to be a training substrate for the population gDNA hyperprior
    — it is not the deliverable, and it does not have to answer objects it cannot solve. The kernel reads
    the curve at every cell itself, on each slot's own gDNA support (`region_geometry.region_gdna_geometry`).

    ``intron_prior`` is the intron factory (`calibrate.FactoryRows`: its inputs, the rows built per block
    inside the kernel), one ``(n, K)`` array of rows (a gate's), or ``None``. ⛔ ψ carries NO reference
    location: the reference is the symmetric Jeffreys measure and asserts nothing; background information
    enters as the factory's λ-factor, a likelihood whose precision scales with counts.
    """
    policy = policy if policy is not None else SilentPolicy()
    n = int(chain.n_slots)
    structure = _structure(chain, statics, region_arrays)
    kappa = float(rna_sense_frac)
    disc = strand_discriminability(kappa, n_rna_obs)
    # THE LIBRARY — the policy's reductions over the WHOLE chain, once, from observations and geometry
    # alone: the only information a message may carry across a locus boundary
    view = ChainView(
        **view_fields(chain, statics, geometry, structure),
        n_grid=int(n_grid),
        logodds_window=float(logodds_window),
        strand_live=disc > 0.0,
    )
    library = policy.library(view)
    lam, _ = _logodds_grid(int(n_grid), float(logodds_window))
    blocks = locus_blocks(chain, structure.terminal, block_slots)
    factory = _factory_of(intron_prior)
    cache = None if _capture is not None else message_cache
    keys = served = None
    if cache is not None:
        keys = [
            cache.key(
                view, b, belief.f_g, None if factory is None else factory.digest(b), library, policy
            )
            for b in blocks
        ]
        served = [cache.served(k) for k in keys]

    out = {
        k: np.ascontiguousarray(getattr(belief, k), dtype=np.float64).copy()
        for k in ("f_pos", "f_neg", "f_g", "var_gdna")
    }
    has_composition = np.zeros(n, dtype=bool)
    diag = None
    if _capture is not None:  # the kernel's diagnostics mode fills these
        diag = {k: np.zeros(n) for k in ("fg_loc", "fg_strand", "tau_lam", "tau_fac")}
        diag["lam_rows"] = np.zeros((n, lam.shape[0]))
    strand = policy.strand
    res = solve_blocks(
        **_chain_arrays(view, structure.terminal),
        belief_fpos=np.ascontiguousarray(belief.f_pos, np.float64),
        belief_fneg=np.ascontiguousarray(belief.f_neg, np.float64),
        belief_fg=np.ascontiguousarray(belief.f_g, np.float64),
        blocks=np.array([(b.start, b.stop, b.end) for b in blocks], np.int64).reshape(-1, 3),
        kappa=kappa,
        od_g=float(gdna_strand_overdispersion),
        od_r=float(rna_strand_overdispersion),
        disc=float(disc),
        lam=lam,
        n_tilt=int(_TILT_NODES),
        policy=_POLICY_KERNEL[policy.name],
        has_strand=strand is not None,
        policy_kappa=float(strand[0]) if strand is not None else 0.0,
        policy_od_g=float(strand[1]) if strand is not None else 0.0,
        policy_od_r=float(strand[2]) if strand is not None else 0.0,
        rho_gdna=float(library.rho_gdna) if library is not None else 0.0,
        rho_rna=float(library.rho_rna) if library is not None else 0.0,
        split_live=bool(library.split_live) if library is not None else False,
        gdna=None
        if gdna_prior is None
        else (
            np.ascontiguousarray(gdna_prior.log_rho, np.float64),
            np.ascontiguousarray(gdna_prior.logP, np.float64),
        ),
        factory=None if factory is None else factory.kernel(),
        served=[None] * len(blocks) if served is None else served,
        out_fpos=out["f_pos"],
        out_fneg=out["f_neg"],
        out_fg=out["f_g"],
        out_var=out["var_gdna"],
        out_has_composition=has_composition,
        deliveries=cache is not None,
        diagnostics=diag,
        n_threads=int(n_threads),
    )
    n_owned = np.array([b.stop - b.start for b in blocks], np.int64)
    delivered = {k: np.asarray(res[k], bool) for k in ("rows_delivered", "cube_delivered")}
    counts = AssertionCounts.of_blocks(res["assertions"], res["counts"], delivered, n_owned)
    if cache is not None:
        for key, delivery in zip(keys, res["deliveries"]):
            if delivery is not None:
                cache.put(key, delivery)
    if _capture is not None:  # inert diagnostic hook
        _capture.fill(
            _gather_diagnostics(chain, view, belief, out, diag, res, blocks, geometry, lam)
        )
        _capture.backbone_assertions = counts
        # which policy ran, read off the artifact — the witness an instrument's "the arm ran"
        # assertion needs, never a config flag it did not thread
        _capture.policy_name = str(policy.name)
        _capture.solve_grid = _logodds_grid(int(n_grid), float(logodds_window))[1]
    return RegionBelief(**out, has_composition=has_composition)


def _chain_arrays(view: ChainView, terminal) -> dict:
    """The chain's arrays by the kernel's argument names, each C-contiguous in the kernel's dtype — the
    observations and the geometry of the :class:`~.messages.ChainView`, and the terminal bits."""

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
class _Structure:
    """The structural classes per slot, from the signature and once for the chain: the EXON region (the
    SPLICE IN's destination class, a policy input rather than a gate); the signature's two EXON bits per
    slot — a strand's opportunity geometry, which the RNA level lanes read to tell a strand's intron
    from its exon inside an overlapping locus; and the TERMINAL — a region with no admissible RNA
    strand, the same predicate the SOLVE gate locks — which receives nothing and where the chain
    breaks."""

    is_exon_region: np.ndarray
    exon_pos: np.ndarray
    exon_neg: np.ndarray
    terminal: np.ndarray


def _structure(chain, statics, region_arrays) -> _Structure:
    is_region = np.asarray(chain.kind) == REGION
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, rtype.shape[0] - 1)
    sig = np.asarray(region_arrays.signature).astype(np.int64)[ri]
    fp, fn = np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool)
    return _Structure(
        is_exon_region=is_region & (rtype[ri] == 2),
        exon_pos=is_region & ((sig & BIT_EXON_POS) > 0),
        exon_neg=is_region & ((sig & BIT_EXON_NEG) > 0),
        terminal=is_region & g1_locked(fp, fn),
    )


class _RowsOfArray:
    """The λ-factor rows given as ONE ``(n, K)`` array (a gate's synthetic rows): handed to the kernel whole
    and digested per block by content — the two calls `calibrate.FactoryRows` answers from its inputs."""

    __slots__ = ("_rows",)

    def __init__(self, rows):
        self._rows = np.ascontiguousarray(rows, np.float64)

    def kernel(self) -> tuple:
        return ("rows", self._rows)

    def digest(self, block) -> bytes:
        a = np.ascontiguousarray(self._rows[block.start : block.end])
        h = hashlib.blake2b(digest_size=16)
        h.update(f"{a.dtype.str}{a.shape}".encode())
        h.update(a)
        return h.digest()


def _factory_of(intron_prior):
    """The λ-factor as the sweep hands it to the kernel — ``kernel()`` and ``digest(block)`` — from the
    factory (`calibrate.FactoryRows`: the rows built per block inside the kernel from their inputs, the
    inputs digested) or from one array of rows; ``None`` is no factory."""
    if intron_prior is None:
        return None
    return intron_prior if hasattr(intron_prior, "digest") else _RowsOfArray(intron_prior)


def _gather_diagnostics(chain, view, belief, out, diag, res, blocks, geometry, lam) -> SweepCapture:
    """The diagnostic capture of one sweep — the instruments' view (:class:`~.blocks.SweepCapture`),
    assembled from the kernel's per-slot arrays and its per-block deliveries and received tables. One extra
    solve lives in the kernel's capture mode and nowhere in production: the strand-ONLY belief (no prior, no
    messages), to split the local error into the strand likelihood against the prior's contribution."""
    fp, fn = np.asarray(view.free_pos, bool), np.asarray(view.free_neg, bool)
    mass_global, eff_global = region_gdna_geometry(geometry)
    cubes = []
    for b, delivery in zip(blocks, res["deliveries"]):
        if delivery is None or delivery[3] is None:
            continue
        slot, pos, has_pos, neg, has_neg, total, opportunity, rho = delivery[3]
        keep = np.asarray(slot) < (
            b.stop - b.start
        )  # the block's own slots, not its read-ahead terminal
        cubes.append(
            CubeRows(
                slot=np.asarray(slot)[keep] + b.start,
                profile_pos=np.asarray(pos)[keep],
                has_pos=np.asarray(has_pos)[keep],
                profile_neg=np.asarray(neg)[keep],
                has_neg=np.asarray(has_neg)[keep],
                u=lam,
                total=np.asarray(total)[keep],
                opportunity=np.asarray(opportunity)[keep],
                rho_ref=np.asarray(rho)[keep],
            )
        )
    received = [r for r in res["received"] if r is not None]
    delivered_rows = bool(np.asarray(res["rows_delivered"], bool).any())
    return SweepCapture(
        fg_loc=diag["fg_loc"],
        fg_strand=diag["fg_strand"],
        f_g=out["f_g"].copy(),
        var_g=out["var_gdna"].copy(),
        tau_lam=diag["tau_lam"],
        tau_fac=diag["tau_fac"],
        fg_init=np.asarray(belief.f_g, np.float64),
        solvable=(fp | fn) & (view.n_slot > 0.0),
        count=np.asarray(view.unspliced_count, np.float64),
        spliced=view.spliced_slot,
        mature=np.asarray(view.sj_count, np.float64).sum(axis=1),
        free_pos=fp,
        free_neg=fn,
        eff_gdna=np.asarray(view.eff_gdna, np.float64),
        eff_rna=np.asarray(view.eff_rna, np.float64),
        mass_global=mass_global,
        eff_global=eff_global,
        lam_rows=diag["lam_rows"] if delivered_rows else None,
        cube_rows=CubeRows.concat(cubes) if cubes else None,
        from_left=SweepCapture.concat_received([r[0] for r in received]) if received else None,
        from_right=SweepCapture.concat_received([r[1] for r in received]) if received else None,
        left=np.asarray(chain.left, np.int64),
        right=np.asarray(chain.right, np.int64),
    )


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# THE OUTPUT CONTRACT — projecting the chain belief back onto the two payload axes.
# Not part of the solve: this is what ``CalibrationResult`` / ``priors`` / ``derive`` consume.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def chain_region_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's REGION slots back onto the REGION axis as a :class:`RegionDeconv` — what
    ``CalibrationResult`` / ``priors`` / ``derive`` consume.

    A region's contained population carries no spliced term, and that is structural: the accumulator
    credits ``region_contained`` only when the fragment used no sj, so a contained fragment is unspliced
    by construction. ⛔ Do not add a ``+ mass_spliced`` term here — the quantity is identically zero on
    the region axis, so adding it would be adding a channel that cannot exist.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    reg = kind == REGION
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    n = count.shape[0]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return RegionDeconv(
        gdna_mass=f_g * count,
        rna_mass=(1.0 - f_g) * count,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_boundary_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's BOUNDARY slots onto the CONTIGUOUS-BOUNDARY axis — the crossing flux that
    ``priors`` and ``derive`` consume.

    ONE per-boundary result, not a ``(left, right)`` pair of per-region ones: splitting each boundary's
    flux onto its two flanking regions only for ``priors`` to pool the halves back together is a no-op,
    and that sum-then-halve pattern is what hides a factor of 2. ``CalibrationResult`` carries
    per-boundary arrays for this reason.

    The RNA mass is spliced-inclusive: a boundary's certified-RNA crossings (``boundary_spliced``) are RNA
    whatever the unspliced mixture resolves to, since gDNA cannot be spliced.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    boundary = kind == BOUNDARY
    unspliced = np.asarray(substrate.boundary_unspliced.count, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(axis=1)
    n = unspliced.shape[0]
    ei = idx[boundary]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g[ei] = np.asarray(belief.f_g, dtype=np.float64)[boundary]
    # THE PER-STRAND RNA SPLIT, PROJECTED ON THIS AXIS TOO. ψ solves the simplex
    # ``(f_g, f_pos, f_neg)`` at EVERY slot — AXIOM 0's `T(slot)`, which is a function of the two
    # `free_*` bits and never of the slot's kind — so a BOUNDARY slot has the same three-way composition
    # a REGION slot does. ⛔ Emitting zeros for both RNA strands here would publish a composition on the
    # crossing axis that sums to ``f_g`` alone, which a per-transcript prior reading composition per
    # object then believes.
    f_pos[ei] = np.asarray(belief.f_pos, dtype=np.float64)[boundary]
    f_neg[ei] = np.asarray(belief.f_neg, dtype=np.float64)[boundary]
    return RegionDeconv(
        gdna_mass=f_g * unspliced,
        rna_mass=(1.0 - f_g) * unspliced + spliced,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )
