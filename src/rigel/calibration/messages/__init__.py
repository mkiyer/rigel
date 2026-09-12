"""rigel.calibration.messages — the message-composition POLICY, behind one interface.

       Gate: ``tests/calibration/test_sweep_backbone.py``

The backbone (:mod:`rigel.calibration.sweep`) owns the SHAPE of the solve — the self-solve, two
directional passes over the ``N E N E … N`` chain, one ψ solve, one write-back, and four assertions.
Everything about *what a message says* is a policy, and it lives here.

Two policies exist; `CalibrationConfig.message_policy` selects which one `message_propagation = True`
installs, and an unknown name raises. The default is `"transfer"`.

* :class:`~.transfer.TransferPolicy` — the shipped policy: every message is either a composition
  profile carried across one face by a derived map, or a population's LEVEL carried where composition
  cannot cross, each hop priced by the two nodes' counting and their own disagreement
  (`transfer_rows` holds the pure row constructors).
* :class:`~.silent.SilentPolicy` — sends nothing; the OFF state and the measured floor. Five
  boundaries long: a reader who holds ``sweep.py`` plus ``silent.py`` in their head holds the entire
  working system.

THE TWO PHASES. Phase 1, PROPAGATE: a forward pass then a backward pass; at each hop the RECIPIENT
receives what its neighbour sends — the sender's own claim composed with what the sender holds from
its far side — and decides to STOP, FORWARD or MODIFY it; beliefs do not change; when both passes end
every node holds one message from each neighbour it has. Phase 2, SOLVE: every node once, from its own
evidence, the two held messages and the gDNA hyperprior. The backbone owns the passes and the solve's
shape; a policy owns what a message says and what a recipient does with it.

The interface
-------------
::

    library  = policy.library(view)               # once per sweep, over the WHOLE chain: the only
                                                  #   cross-block reductions a message may use
    prepared = policy.prepare(ctx, library)       # per block: every node's OWN claim
    receive  = prepared.propagate(received, backward=False)  # phase 1: the recipient's kernel writing
                                                  #   rows of the pass's table, or None ⇒ all silence
    receive(source, destination)                  # ... the BACKBONE runs the pass, in chain order
    evidence = prepared.solve(from_left, from_right)   # phase 2, the policy's half -> PsiMessage

The chain is solved a LOCUS BLOCK at a time (`sweep.solve_chain`, `region_chain.locus_blocks`), and the
two calls above are the two scopes a policy sees: ``library`` reads a :class:`ChainView` of the whole
chain — observations and geometry, NO beliefs, which is what makes a cross-block reduction over beliefs
unwritable — and returns whatever library-wide facts its messages need (the transfer policy's: three
reference densities and whether the strand split is live). ``prepare`` reads a :class:`BlockContext` of
one block, beliefs included, plus that library. Everything else a message reads is per slot or per face.

⛔ THE CONTRACT, and it is TRAPS: a-message-from-the-destinations-belief, a lesson that has recurred
nine times in nine costumes:

    A message may use the destination's **CONSTANTS** (geometry, effective lengths) and its
    **OBSERVATIONS** (counts, mass). It may **NEVER** use the destination's **BELIEFS**.

The rule is about SENDING, not about RECEIVING, which resolves what otherwise reads as a contradiction
with the shipped policy:

    A SENDER publishes its claim UNCHANGED: it does not tailor, scale or hedge it for whoever is
    listening. Deciding how much of an arriving claim to BELIEVE is the RECIPIENT's job, and a
    recipient necessarily reads its own belief to do it.

So a destination that receives a claim wildly at odds with what its own data says may DISCOUNT it, and
that is reception rather than a message built from the destination — the transfer policy's hop price
does exactly this: the recipient's own counts price the arriving claim's width. The line the trap
draws is that a claim's VALUE may never be built from the destination's belief, because that
manufactures agreement out of nothing. A reception step is safe when it can only ever WIDEN a claim
and never move its mode — it can discard information, never invent it.

:class:`BlockContext` splits its fields under exactly those three headings, and the heading is what
turns the trap from a discipline into something a reader — and the backbone — can check. The backbone
enforces the half that is enforceable: the kernel is called with two INDICES and builds the message
into the destination's row from the SOURCE's claim and what the source holds (its own row of the same
table, written one step earlier); the backbone owns the table and the order, and the policy never
reaches past its hop, so a message built from the destination's belief has nowhere to come from.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "BlockContext",
    "ChainView",
    "Levels",
    "Policy",
    "Prepared",
    "PsiMessage",
    "Received",
]


@dataclass(frozen=True, slots=True)
class PsiMessage:
    """What the two neighbours jointly tell ψ about this slot, and **nothing else** — two row channels,
    each a max-normalised log-profile on the solve grid, each ``None`` when the policy has no claim:

    * ``lam_rows`` — ``(n_slots, K)`` over ψ's log-odds grid ``λ``, added into the FINAL solve only
      (never phase-A, never the own-evidence precision);
    * ``cube_rows`` — ``{slot: (K, K_t) row}`` over the ``(λ, θ)`` cube, for AMBIG slots only, added
      inside the AMBIG solve the same way.

    A profile on ψ's own grid cannot be delivered off-grid and cannot claim an over-unit share, so the
    lesson of TRAPS: off-grid-message-mode is structural on this channel rather than asserted. A
    fully-``None`` message is :meth:`silent` — the floor the whole message layer is priced against.
    """

    #: THE λ ROWS — an ``(n_slots, K)`` λ-factor row array in ψ's general evidence currency
    #: (θ-independent, finite, an all-zero row is inert), or ``None`` for no claim: what the two held
    #: messages say about a slot's composition, delivered as a row over the solve grid (the transfer
    #: policy's `solve` sums the held profiles per slot into it). The backbone adds the rows into the
    #: FINAL solve only: never phase-A, never the own-evidence precision — that citizenship is the
    #: entire difference from the intron factory's factor.
    lam_rows: np.ndarray | None = None
    #: THE CUBE CHANNEL: ``{slot: (K, K_t) row}`` for AMBIG slots only — a max-normalised log-profile
    #: over ψ's ``(λ, θ)`` cube, the delivery of the RNA LEVEL lanes at a node where both strands are
    #: live (each held strand level read at the density every cell implies: ``f_s = (1 − σ)(1 ± τ)/2``,
    #: ``ρ_s = f_s n / a_r``). The backbone adds it to ψ inside the AMBIG solve, FINAL solve only, like
    #: ``lam_rows``; ``None`` or an absent slot leaves the solve exactly as it is without the channel.
    #: A single-strand slot has no cube and may not appear here.
    cube_rows: dict | None = None

    @classmethod
    def silent(cls) -> PsiMessage:
        """No claim on either channel — what :class:`~.silent.SilentPolicy` delivers, the measured floor."""
        return cls()

    @property
    def is_silent(self) -> bool:
        return self.lam_rows is None and self.cube_rows is None


@dataclass(slots=True)
class Levels:
    """One population's LEVEL lane as RECEIVED: row ``i`` is the level node ``i`` holds on this lane from
    one side after a pass, or nothing (``present[i]`` False). A level is a max-normalised log-profile
    over ``u = log(rho / rho_ref)`` — the population's density in counts per base of its opportunity,
    relative to the library's structurally pure gDNA density ``rho_ref`` — on the solve grid (``K``
    points, the ``lam`` window: a coordinate choice, no constant). ``count`` and ``opportunity`` are the
    total and the opportunity of the last node WITH a total the claim passed through: the next recipient
    prices its hop from them — both totals' counting, and the abundance discrepancy beyond it, per hop
    and nothing pooled. An EMPTY node (no total) forwards a level unchanged and leaves them as they
    were: a few bases of the same gDNA density — unless it is itself a flux SOURCE on an RNA lane, whose
    level travels with the flux's witness (the spliced count on the route rate's opportunity).

    ``rna_count`` / ``rna_count_var`` are an RNA lane's witness of ITS strand's abundance at that same
    last full node — the strand's RNA count read from the node's column split (its asymmetry over the
    protocol's strand contrast) and that estimate's Poisson variance — the pair the next recipient's
    price compares with its own split. ``has_witness`` is False on the gDNA lane and where the library's
    strand channel is dead (the derived deadband), where the column count is the witness."""

    present: np.ndarray  # (n,) bool
    profile: np.ndarray  # (n, K) f64
    count: np.ndarray  # (n,) f64
    opportunity: np.ndarray  # (n,) f64
    has_witness: np.ndarray  # (n,) bool
    rna_count: np.ndarray  # (n,) f64, NaN without a witness
    rna_count_var: np.ndarray  # (n,) f64, NaN without a witness

    @classmethod
    def empty(cls, n: int, K: int) -> Levels:
        return cls(
            np.zeros(n, bool),
            np.zeros((n, K)),
            np.zeros(n),
            np.zeros(n),
            np.zeros(n, bool),
            np.full(n, np.nan),
            np.full(n, np.nan),
        )

    def write(self, i, profile, count, opportunity, rna_count=None, rna_count_var=None) -> None:
        """Row ``i`` holds this level (a hop wrote it, or re-priced what it had written)."""
        i = int(i)
        self.present[i] = True
        self.profile[i] = profile
        self.count[i], self.opportunity[i] = float(count), float(opportunity)
        self.has_witness[i] = rna_count is not None
        self.rna_count[i] = np.nan if rna_count is None else float(rna_count)
        self.rna_count_var[i] = np.nan if rna_count_var is None else float(rna_count_var)

    def forward(self, s, i) -> None:
        """Row ``i`` holds exactly what row ``s`` holds — an EMPTY node forwards a level unchanged."""
        s, i = int(s), int(i)
        for f in dataclasses.fields(self):
            getattr(self, f.name)[i] = getattr(self, f.name)[s]

    def take(self, n: int) -> Levels:
        return Levels(*(getattr(self, f.name)[:n] for f in dataclasses.fields(self)))

    @classmethod
    def concat(cls, parts: list) -> Levels:
        return cls(
            *(
                np.concatenate([getattr(q, f.name) for q in parts], axis=0)
                for f in dataclasses.fields(cls)
            )
        )


@dataclass(slots=True)
class Received:
    """What every node of a block RECEIVED from one side after one pass — the transfer policy's message,
    as a table: row ``i`` is what node ``i`` holds from its neighbour on that side. ``from_left`` and
    ``from_right`` are two of these per block.

    THE LANES. A node's unknown is its COMPOSITION on the simplex ``(f_g, f_+, f_-)`` — two degrees of
    freedom where both strands are live, one where a single strand is — and, where composition cannot
    cross a face, the LEVELS of the three populations. So a row carries up to four lanes, every one
    optional:

    * ``composition`` (where ``has_composition``) — the gDNA-versus-RNA PROFILE: a max-normalised
      log-likelihood over the solve grid of the destination's gDNA share (``lam = log f_g/(1-f_g)``,
      ``K = n_grid`` points). Scale-free, so it crosses a face by a derived map and never carries a
      level across a capture cliff.
    * ``level_gdna``, ``level_rna_pos``, ``level_rna_neg`` — a LEVEL claim per population
      (:class:`Levels`: a PROFILE over the log density relative to the library's structurally pure gDNA
      density, on the same grid as ``lam``), for faces composition cannot cross: gDNA is genomically
      continuous across ANY face; a strand's RNA continues across a face where that strand's
      population is unchanged (an AMBIG region's two degrees of freedom are imputed by exactly these).
      A level is ABSOLUTE — it needs no map and no knowledge of its recipient, which is what lets it
      cross a node that has no total at all — and it is a profile, not a Gaussian pair, because the
      claims that travel on it are LOWER-ONLY (a level says "at least this much gDNA") and a Gaussian
      summary of a one-sided profile invents a value.

    THE RNA LANES: a single-strand node's own claim read as its live strand's RNA level, and the
    certified flux at an exon's junctions as that strand's level at the exon; per-strand faces from the
    flag bits; two-sided only between an intron and its own boundary; delivered at AMBIG nodes on ψ's
    cube (`PsiMessage.cube_rows`). The tilt — the RNA+ versus RNA− degree of freedom at a both-stranded
    node — has NO lane: the two RNA levels constrain it inside the cube, and a tilt profile from the
    same witnesses would count them twice.

    THE TWO STATES the table expresses where a node ``heard`` nothing. SILENCE (:attr:`silence`): the
    node HAS a neighbour on this side (``has_neighbour``) and nothing is present — the neighbour spoke
    and had nothing to say, or the node is a terminal, which receives nothing. NO NEIGHBOUR
    (:attr:`no_neighbour`): the side is open — a reference start or end, or a block's edge — and there
    was no hop at all. The backbone writes ``has_neighbour``; a policy's kernel writes the lanes and
    nothing else; the two states need no word of their own.
    """

    has_neighbour: np.ndarray  # (n,) bool — the backbone's: the side exists
    has_composition: np.ndarray  # (n,) bool
    composition: np.ndarray  # (n, K) f64, where has_composition
    level_gdna: Levels
    level_rna_pos: Levels
    level_rna_neg: Levels

    #: the three level lanes, by field name
    LANES = ("level_gdna", "level_rna_pos", "level_rna_neg")

    @classmethod
    def empty(cls, n: int, K: int) -> Received:
        return cls(
            np.zeros(n, bool),
            np.zeros(n, bool),
            np.zeros((n, K)),
            Levels.empty(n, K),
            Levels.empty(n, K),
            Levels.empty(n, K),
        )

    @property
    def has_level(self) -> np.ndarray:
        """A level is present on some lane — a BOUND arrived (never a composition)."""
        return self.level_gdna.present | self.level_rna_pos.present | self.level_rna_neg.present

    @property
    def heard(self) -> np.ndarray:
        """Something arrived on this side: a composition or a level."""
        return self.has_composition | self.has_level

    @property
    def silence(self) -> np.ndarray:
        """A neighbour on this side, and nothing heard from it."""
        return self.has_neighbour & ~self.heard

    @property
    def no_neighbour(self) -> np.ndarray:
        """No side to hear from: an open end of the chain, or of the block."""
        return ~self.has_neighbour

    def take(self, n: int) -> Received:
        """The first ``n`` rows — a block's OWNED slots, without the terminal it read beyond them."""
        return Received(
            self.has_neighbour[:n],
            self.has_composition[:n],
            self.composition[:n],
            *(getattr(self, lane).take(n) for lane in self.LANES),
        )

    @classmethod
    def concat(cls, parts: list) -> Received:
        """The blocks' tables as the chain's, in block order."""
        return cls(
            np.concatenate([q.has_neighbour for q in parts]),
            np.concatenate([q.has_composition for q in parts]),
            np.concatenate([q.composition for q in parts], axis=0),
            *(Levels.concat([getattr(q, lane) for q in parts]) for lane in cls.LANES),
        )


@dataclass(frozen=True, slots=True)
class ChainView:
    """A stretch of the chain as a policy may read it WITHOUT beliefs: the observations and the
    geometry, under the two headings that make TRAPS: a-message-from-the-destinations-belief legible,
    plus the solve's own scalars. `Policy.library` receives the WHOLE chain in this form, so the only
    cross-block information a policy can build is a reduction over observations and geometry — a
    reduction over beliefs has no field to read. :class:`BlockContext` adds the beliefs for one block.

    ⛔ The headings are load-bearing. ``observations`` and ``geometry`` may be indexed at either end of
    a hop; ``beliefs`` may be indexed at the SOURCE only. A policy that reads a ``beliefs`` field at the
    destination is committing TRAPS: a-message-from-the-destinations-belief, and the field's heading is
    what makes that visible in review. The shipped policy reads ``belief_fg`` once, at ``prepare``, for
    the variance freeze of each node's OWN strand profile — a source-side read by construction, since
    the profile is the node's claim before any hop.

    Every field here has a reader in the transfer policy or the backbone.
    """

    # ── OBSERVATIONS — readable at either end of a hop ────────────────────────────────────────────────
    eff_gdna: np.ndarray  # per-slot gDNA opportunity (the capture-blind effective length)
    eff_rna: np.ndarray  # per-slot RNA opportunity
    sj_count: np.ndarray  # [n, 2] sj fragment count by TRANSCRIPT strand (both faces)
    #: the same flux split by which genomic END of its sj this boundary is — the count that
    #: matches `route_rate_lo`/`route_rate_hi`, so a per-face rate is priced on its own count
    sj_count_lo: np.ndarray
    sj_count_hi: np.ndarray
    #: [n, 2] the ROUTE-SUMMED certified rate per face by transcript strand (Σ flux_J/A_J over
    #: the face's disjoint routes) — the pooled sj_count/eff_sj ratio under-reads k-route faces
    #: ~k×, so consumers of a face's RATE read these, never the ratio
    route_rate_lo: np.ndarray
    route_rate_hi: np.ndarray
    unspliced_count: (
        np.ndarray
    )  # [n, 2] unspliced count by GENOME strand — the density numerator AND n
    spliced_count: np.ndarray  # [n, 2] spliced count by strand

    # ── GEOMETRY / STRUCTURE — readable at either end, and belief-free by construction ────────────────
    left: np.ndarray  # adjacent slot of the other kind, -1 at a reference start
    right: np.ndarray
    is_boundary: np.ndarray  # ~is_region; the chain strictly alternates N E N E … N
    is_exon_region: (
        np.ndarray
    )  # a REGION whose region signature is EXON — the SPLICE IN's destination
    free_pos: np.ndarray  # does the annotation admit +RNA here?  one of AXIOM 0's TWO BITS
    free_neg: np.ndarray  # …and -RNA?                            the other
    #: the region signature's two EXON bits per slot (False at a boundary) — a strand's OPPORTUNITY
    #: geometry, not a population: a REGION that admits strand ``s`` and carries no exon of ``s`` is
    #: strand ``s``'s INTRON whatever the other strand does there (the h-intron ∩ a-exon piece of an
    #: overlapping locus is the host strand's intron and the antisense's exon at once). The coarse
    #: ``is_exon_region`` cannot say this, and the RNA level lanes need it.
    exon_pos: np.ndarray
    exon_neg: np.ndarray
    #: the terminus and junction bits per BOUNDARY slot (0 at a region): which faces composition may
    #: cross, the outside flank of a terminus, the junction's exon side (`transfer_rows`)
    boundary_flags: np.ndarray

    # ── the solve's own scalars (neither observation nor belief) ──────────────────────────────────────
    n_grid: int
    logodds_window: float
    #: the AMBIG cube's tilt-grid size ``K_t`` (``None`` ⇒ ``n_grid``, as the solver reads it) — what a
    #: policy needs to lay a ``cube_rows`` row on the grid ψ will evaluate it on
    n_tilt: int | None = None
    #: the intron factory's per-slot λ-factor rows on THIS grid, ``(n_slots, K)`` or ``None`` — the
    #: same array ψ adds as its own λ-factor (`sweep.solve_chain`'s ``intron_prior``), so an intron's
    #: own claim and the solver's factor cannot drift apart; ``None`` is no factory
    factory_rows: np.ndarray | None = None
    #: the derived strand deadband's verdict for the LIBRARY: does the strand split carry composition
    #: information at all (`region_init.strand_discriminability` > 0)? A κ within its noise of ½ makes
    #: every single-strand exon's strand precision exactly zero, and a policy reading the split as an
    #: RNA witness must know that before it prices a single hop
    strand_live: bool = False

    @property
    def n_slots(self) -> int:
        return int(self.unspliced_count.shape[0])

    @property
    def n_slot(self) -> np.ndarray:
        """The per-slot unspliced count over both strands — the density numerator AND the Poisson n:
        one number, not a fractional mass plus a separate integer flux."""
        return self.unspliced_count.sum(axis=1)

    @property
    def spliced_slot(self) -> np.ndarray:
        """The per-slot spliced count over both strands."""
        return self.spliced_count.sum(axis=1)

    def population_size(self) -> np.ndarray:
        """``|T(slot)|`` — AXIOM 0 made arithmetic.

        ``T(slot) = {gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}``, so the size is
        ``1 + free_pos + free_neg`` and is ≤ 3 for every slot, always, because it is a function of TWO
        BITS. That is what makes the three-population rule structural rather than something to remember.
        ⛔ There is no fourth population: "mature" and "nascent" are not species, and RNA inside an intron
        is RNA that has not spliced *at that position*.
        """
        return (
            1
            + np.asarray(self.free_pos, bool).astype(np.int64)
            + np.asarray(self.free_neg, bool).astype(np.int64)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BlockContext(ChainView):
    """One block of the chain as `Policy.prepare` reads it: the :class:`ChainView` plus the BELIEFS —
    SOURCE-SIDE ONLY (TRAPS: a-message-from-the-destinations-belief)."""

    #: does this node have OWN composition evidence — `RegionInit.tau_lam > 0`, the one bit of the
    #: message-free self-solve a policy may know: the strand term (the node has counts and the library's
    #: deadband is open) or the factory's row, so it does not depend on the prior a sweep carries. A
    #: policy is not handed the self-solve's fractions or precisions, and that is what lets the message
    #: layer be shared across the refit sweeps (`message_cache.MessageCache`): every input it reads is on this
    #: context and can be digested.
    has_own_composition: np.ndarray
    belief_fg: np.ndarray  # the INCOMING belief: the variance freeze of a node's own strand profile


@runtime_checkable
class Prepared(Protocol):
    """A policy's per-sweep working object: every node's own claim, the propagate kernel, the solve."""

    def propagate(self, received: Received, *, backward: bool):
        """PHASE 1. Return ``receive(source, destination)`` for one direction — the kernel that writes
        row ``destination`` of ``received``, the pass's table — or ``None`` when this policy sends
        nothing (every node with a neighbour then holds silence from that side).

        The BACKBONE owns the table and runs the pass: it marks ``has_neighbour`` and, in chain order,
        for every destination with a neighbour on that side that is not a terminal, calls
        ``receive(source, destination)``. Inside ``receive`` the policy composes what the source sends
        — its own claim with what the source holds from ITS far side, row ``source`` of the same
        table, written by this same pass one step earlier — and applies the recipient's decision for
        the face: STOP (write nothing: the row stays silent), FORWARD, or MODIFY.

        ONE pass per direction: the forward pass reads each node's LOW neighbour and the backward pass
        its HIGH one; on a chain that IS forward-backward, and nothing here iterates.
        """

    def solve(self, from_left: Received, from_right: Received) -> PsiMessage:
        """PHASE 2, the policy's half: the ψ channels at every slot from the two tables — row ``i`` of
        ``from_left`` is what slot ``i`` holds from its LOW neighbour (no neighbour at a reference
        start), of ``from_right`` from its HIGH one. Never the destination's belief
        (TRAPS: a-message-from-the-destinations-belief)."""


@runtime_checkable
class Policy(Protocol):
    """A message-composition policy. ``name`` appears in diagnostics and in arm output."""

    name: str

    def library(self, view: ChainView):
        """Once per sweep, over the WHOLE chain: whatever library-wide facts this policy's messages
        need, reduced from observations and geometry alone (the view carries no belief), or ``None``.
        This is the ONLY place a policy may look across the chain; ``prepare`` sees one block."""

    def prepare(self, ctx: BlockContext, library) -> Prepared:
        """Derive whatever this policy needs for the block ``ctx`` covers, given its own ``library``:
        every node's OWN claim, the rules per face, the lanes."""
