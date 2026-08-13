"""rigel.calibration.messages — the message-composition POLICY, behind one interface.

       Gate: ``tests/calibration/test_sweep_backbone.py``

The backbone (:mod:`rigel.calibration.sweep`) owns the SHAPE of the solve — two directional scans over the
``N E N E … N`` chain, one combine, one ψ solve, one write-back, and five assertions. Everything about
*what a message says* is a policy, and it lives here.

Three policies ship:

* :class:`~.silent.SilentPolicy` — sends nothing. ⭐ **THE DEFAULT**, and five lines long. A reader who
  holds ``sweep.py`` plus ``silent.py`` in their head holds the entire working system.
* :class:`~.head.HeadPolicy` — every operator the evolved solver carried, each behind a NAMED switch, so
  ``ladder_arm_ab.py`` can price them ONE AT A TIME instead of as a block.
* :mod:`~.variance` — the shared variance arithmetic both of the above draw on. Not a policy; a toolbox.

⭐⭐ **WHY THE SPLIT IS SHAPED THIS WAY, and it is a measurement rather than a taste.** The message layer
at the prior-free pass is worth **+0.2 %** of the shipped answer while moving that pass's own error by
**77.5 %**; muted everywhere it is a net *harm* on three of the four strata and its entire value sits in
one. So the operators in ``head.py`` are not load-bearing as a group — they have to be priced individually,
and a switch per operator is the only way to do that.

The interface
-------------
::

    relay = policy.prepare(ctx)                  # one working object per sweep
    step, publish = relay.scan(backward=False)   # the per-hop kernel; None ⇒ nothing to relay
    ...                                          # the BACKBONE runs the loop
    msg = relay.deliver(left_state, right_state) # -> PsiMessage

⛔⛔ **THE CONTRACT, AND IT IS TRAPS: a-message-from-the-destinations-belief — a lesson that has recurred NINE times in nine costumes:**

    A message may use the destination's **CONSTANTS** (geometry, effective lengths) and its
    **OBSERVATIONS** (counts, mass). It may **NEVER** use the destination's **BELIEFS**.

:class:`StepContext` splits its fields under exactly those three headings, and the heading is what turns
TRAPS: a-message-from-the-destinations-belief from a discipline into something a reader — and the backbone — can check. The backbone enforces the
half that is enforceable: it hands :meth:`deliver` the relayed belief state **already indexed at the
SOURCE**, so no policy can read a neighbour's relayed belief at the destination however it is written.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "NeighbourState",
    "Policy",
    "PsiMessage",
    "Relay",
    "StepContext",
]


@dataclass(frozen=True, slots=True)
class PsiMessage:
    """What the two neighbours jointly tell ψ about this slot, and **nothing else**.

    Four channels, each a ``(mode, precision)`` pair, each in ITS OWN COORDINATE — which is the whole
    reason the modes are stated separately rather than read back off a fused density:

    ==================  =====================================  ==============================
    channel             coordinate of the mode                 grid domain
    ==================  =====================================  ==============================
    ``gdna_*``          ``log`` of the gDNA SHARE of the mass  ``(-inf, 0]`` for a real share
    ``rna_*``           the same, per strand — a 2-tuple       ``(-inf, 0]`` per strand
    ``lam_*``           ``lambda = log(f_g / f_R)``            ``[-L, +L]``, L = logodds_window
    ``theta_*``         the tilt ANGLE ``arcsin(tau)``         ``[-pi/2, +pi/2]`` exactly
    ==================  =====================================  ==============================

    ⛔⛔ **THE COORDINATE IS NOT A DETAIL — it is TRAPS: off-grid-message-mode, which cost 74 % of a zero control's
    error.** A mode delivered outside its grid's domain is not a weak claim, it is a PIN AT THE BOUNDARY:
    the penalty ``-1/2 p (x - m)^2`` with ``m`` off-grid is monotone across every grid point, so it has no
    interior minimum and precision buys a corner rather than a location. The backbone asserts the domain
    for exactly this reason.

    ``None`` in any field means "this channel says nothing", and a fully-``None`` message is
    :meth:`silent` — the floor the whole message layer is priced against.
    """

    gdna_mode: np.ndarray | None = None
    gdna_prec: np.ndarray | None = None
    rna_mode: tuple[np.ndarray, np.ndarray] | None = None
    rna_prec: tuple[np.ndarray, np.ndarray] | None = None
    lam_mode: np.ndarray | None = None
    lam_prec: np.ndarray | None = None
    theta_mode: np.ndarray | None = None
    theta_prec: np.ndarray | None = None

    @classmethod
    def silent(cls) -> PsiMessage:
        """No claim on any channel. ⭐ Byte-identical to muting ψ's four imputed arguments, which is what
        makes :class:`~.silent.SilentPolicy` the measured floor rather than a new code path."""
        return cls()

    @property
    def is_silent(self) -> bool:
        return all(
            getattr(self, f) is None
            for f in (
                "gdna_mode",
                "gdna_prec",
                "rna_mode",
                "rna_prec",
                "lam_mode",
                "lam_prec",
                "theta_mode",
                "theta_prec",
            )
        )


@dataclass(frozen=True, slots=True)
class NeighbourState:
    """One neighbour's relayed belief, **already indexed at the SOURCE**, plus the validity mask.

    ⭐⭐ **THIS TYPE IS THE TRAPS: a-message-from-the-destinations-belief ASSERTION.** The backbone gathers each relayed array at the source slot
    before handing it over, so a policy holding a :class:`NeighbourState` is holding values for the source
    and has no way to ask the same array about the destination. TRAPS: a-message-from-the-destinations-belief's nine costumes were all a message built
    from the destination's own relayed/fused belief; none of them is expressible through this type.

    ``valid`` is ``False`` at a reference terminal, where ``chain.left`` / ``chain.right`` is ``-1`` and the
    gathered value is whatever the clipped index happened to hit — so it must be masked, never read.
    """

    state: tuple[np.ndarray, ...]
    valid: np.ndarray
    src: np.ndarray


@dataclass(frozen=True, slots=True)
class StepContext:
    """Everything a policy may read, under the three headings that make TRAPS: a-message-from-the-destinations-belief legible.

    ⛔ **The headings are load-bearing.** ``observations`` and ``geometry`` may be indexed at either end of
    a hop; ``beliefs`` may be indexed at the SOURCE only. A policy that reads a ``beliefs`` field at the
    destination is committing TRAPS: a-message-from-the-destinations-belief, and the field's heading is what makes that visible in review.

    ⚠ **One field in ``beliefs`` is read at the destination by the shipped policy and it is a KNOWN,
    MEASURED DEBT, not an oversight**: ``belief_fg`` reaches the reframe's frame pair, so the frame at a
    hop is a function of the destination's belief. The operator ledger prices it — slots where a *solved*
    belief rather than the ``{0,0,1}`` default sets the frame carry **57–77 % of library mass** — and it is
    named here so the next reader finds it recorded rather than discovers it again.
    """

    # ── OBSERVATIONS — readable at either end of a hop ────────────────────────────────────────────────
    mass: np.ndarray  # per-slot gDNA-support mass (the pin's and the share's denominator)
    eff_gdna_global: np.ndarray  # the matching gDNA opportunity
    eff_rna: np.ndarray  # per-slot RNA effective length
    eff_gdna: np.ndarray  # per-slot gDNA effective length (per-face geometry, diagnostics)
    eff_junction: np.ndarray  # [n, 2] junction opportunity by TRANSCRIPT strand
    junction_count: np.ndarray  # [n, 2] junction fragment count by TRANSCRIPT strand
    unspliced_count: (
        np.ndarray
    )  # [n, 2] unspliced count by GENOME strand — the density numerator AND n
    n_slot: np.ndarray  # unspliced_count.sum(axis=1)
    spliced_slot: np.ndarray  # per-slot spliced count, summed over strands

    # ── GEOMETRY / STRUCTURE — readable at either end, and belief-free by construction ────────────────
    left: np.ndarray  # adjacent slot of the other kind, -1 at a reference start
    right: np.ndarray
    is_boundary: np.ndarray  # ~is_region; the chain strictly alternates N E N E … N
    is_exon_region: np.ndarray  # a REGION whose region signature is EXON — the graft's destination
    free_pos: np.ndarray  # does the annotation admit +RNA here?  ⭐ one of AXIOM 0's TWO BITS
    free_neg: np.ndarray  # …and -RNA?                            ⭐ the other
    g1_locked: np.ndarray  # structurally pure-gDNA, BOTH axes (an intergenic|exon BOUNDARY included)
    boundary_flags: np.ndarray  # for terminus_flank_gain — does a flank's RNA population grow?
    geometry: object  # RegionGeometry, for the frame pair (a policy-owned derivation)
    order: list  # the genomic visiting order — slot ids ARE it, so this is range(n)
    left_list: list  # ``left`` as a Python list: the scan reads it one element at a time
    right_list: list

    # ── BELIEFS — SOURCE-SIDE ONLY (TRAPS: a-message-from-the-destinations-belief) ──────────────────────────────────────────────────────────────
    own: object  # RegionInit: the message-free self-solve — rho_*, prec_*, tau_lam, struct_lock, f_*
    belief_fg: (
        np.ndarray
    )  # the INCOMING belief. ⚠ the frame pair reads it at BOTH ends — the debt above

    # ── the solve's own scalars and fitted library constants (neither observation nor belief) ─────────
    n_grid: int
    logodds_window: float
    solve_grid: np.ndarray
    capture: dict | None = None  # the diagnostics hook; inert in production

    @property
    def n_slots(self) -> int:
        return int(self.n_slot.shape[0])

    def population_size(self) -> np.ndarray:
        """``|T(slot)|`` — AXIOM 0 made arithmetic.

        ``T(slot) = {gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}``, so the size is
        ``1 + free_pos + free_neg`` and is **≤ 3 for every slot, always**, because it is a function of TWO
        BITS. That is what makes the three-population rule structural rather than something to remember.
        ⛔ There is no fourth population: "mature" and "nascent" are not species, and RNA inside an intron
        is RNA that has not spliced *at that position*.
        """
        return (
            1
            + np.asarray(self.free_pos, bool).astype(np.int64)
            + np.asarray(self.free_neg, bool).astype(np.int64)
        )


@runtime_checkable
class Relay(Protocol):
    """A policy's per-sweep working object: prepared arrays, the scan kernel, and the combine."""

    def scan(self, *, backward: bool):
        """Return ``(step, publish)`` for one direction, or ``None`` to relay nothing.

        ``step(s, i)`` performs ONE hop from source slot ``s`` into destination slot ``i``, accumulating in
        place — which is the forward half of forward-backward on a chain. ``publish()`` returns the relayed
        state as a tuple of arrays.

        ⛔ **TRAPS: a-comment-quoted-as-a-finding: "in place" is not "iterative".** Both directions are ONE pass. A source comment
        calling the scan "Gauss-Seidel" meant *un-vectorisable*, that word crossed into a design doc as if
        it were a structural finding, and a reviewer then correctly derived a build plan for a defect that
        does not exist.
        """

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        """The four ψ channels at every slot, from the two NEIGHBOUR states only (TRAPS: a-message-from-the-destinations-belief)."""


@runtime_checkable
class Policy(Protocol):
    """A message-composition policy. ``name`` appears in diagnostics and in arm output."""

    name: str

    def prepare(self, ctx: StepContext) -> Relay:
        """Derive whatever this policy needs from ``ctx``, once per sweep."""
