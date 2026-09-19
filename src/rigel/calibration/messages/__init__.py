"""rigel.calibration.messages — the message-composition POLICY, behind one interface.

       Gate: ``tests/calibration/test_sweep_backbone.py``

The backbone (:mod:`rigel.calibration.sweep`) owns the SHAPE of the solve — the self-solve, two
directional passes over the ``N E N E … N`` chain, one ψ solve, one write-back, and four assertions —
and runs it in ONE native call per sweep (`native.solve_blocks`, `native/solve_kernel.cpp`). Everything
about *what a message says* is a policy, and the policy's arithmetic — its builders, its passes, its
solve — is the kernel's (`native/transfer_kernel.h`). What lives here is what the kernel is TOLD: which
policy runs, the strand model its own claims read, and its LIBRARY.

Two policies exist; `CalibrationConfig.message_policy` selects which one `calibrate` installs, and an
unknown name raises. The default is `"transfer"`.

* :class:`~.transfer.TransferPolicy` — the shipped policy: every message is either a composition
  profile carried across one face by a derived map, or a population's LEVEL carried where composition
  cannot cross, each hop priced by the two nodes' counting and their own disagreement
  (its row constructors are `native/transfer_rows.h`).
* :class:`~.silent.SilentPolicy` — sends nothing; the OFF state and the measured floor: the kernel runs
  no layer for it, and ψ solves every slot on its own evidence and the prior alone. Five lines long: a
  reader who holds ``sweep.py`` plus ``silent.py`` in their head holds the entire working system.

THE TWO PHASES. Phase 1, PROPAGATE: a forward pass then a backward pass; at each hop the RECIPIENT
receives what its neighbour sends — the sender's own claim composed with what the sender holds from
its far side — and decides to STOP, FORWARD or MODIFY it; beliefs do not change; when both passes end
every node holds one message from each neighbour it has. Phase 2, SOLVE: every node once, from its own
evidence, the two held messages and the gDNA hyperprior. The backbone owns the passes and the solve's
shape; a policy owns what a message says and what a recipient does with it.

The interface
-------------
::

    library = policy.library(view)   # once per sweep, over the WHOLE chain: the only cross-block
                                     #   reductions a message may use
    policy.name                      # which layer the kernel runs for the chain's blocks
    policy.strand                    # the strand model its own claims read — (κ, od_gdna, od_rna) or None

The chain is solved a LOCUS BLOCK at a time (`sweep.solve_chain`, `region_chain.locus_blocks`), and
``library`` is the one scope a policy sees in Python: a :class:`ChainView` of the whole chain —
observations and geometry, NO beliefs, which is what makes a cross-block reduction over beliefs
unwritable — from which it returns whatever library-wide facts its messages need (the transfer
policy's: three reference densities and whether the strand split is live). The kernel then sees one
block at a time: the same arrays, the incoming belief, and the library.

⛔ THE CONTRACT:

    A message may use the destination's **CONSTANTS** (geometry, effective lengths) and its
    **OBSERVATIONS** (counts, mass). It may **NEVER** use the destination's **BELIEFS**.

The rule is about SENDING, not about RECEIVING, which resolves what otherwise reads as a contradiction
with the shipped policy:

    A SENDER publishes its claim UNCHANGED: it does not tailor, scale or hedge it for whoever is
    listening. Deciding how much of an arriving claim to BELIEVE is the RECIPIENT's job, and a
    recipient necessarily reads its own belief to do it.

So a destination that receives a claim wildly at odds with what its own data says may DISCOUNT it, and
that is reception rather than a message built from the destination — the transfer policy's hop price
does exactly this: the recipient's own counts price the arriving claim's width. The line the contract
draws is that a claim's VALUE may never be built from the destination's belief, because that
manufactures agreement out of nothing. A reception step is safe when it can only ever WIDEN a claim
and never move its mode — it can discard information, never invent it.

The kernel enforces the half that is enforceable, BY CONSTRUCTION: the pass builds each destination's
row from the SOURCE's claim and what the source holds (its own row of the same received table, written
by this pass one hop earlier), in the backbone's chain order on the backbone's table; and the only
belief the layer reads is the incoming belief at a node's OWN claim — the variance freeze of its own
strand profile, a source-side read, since the profile is the node's claim before any hop. A message
built from the destination's belief has no input to come from.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "ChainView",
    "Policy",
]


@dataclass(frozen=True, slots=True)
class ChainView:
    """The chain as a policy may read it WITHOUT beliefs: the observations and the geometry, under the two
    headings that make that contract legible, plus the solve's own scalars. `Policy.library` receives the
    WHOLE chain in this form, so the only cross-block information a policy can build is a reduction over
    observations and geometry — a reduction over beliefs has no field to read. These are also the arrays
    the backbone hands the kernel (`sweep.solve_chain`), which reads them a block at a time beside the
    incoming belief; and the arrays the message cache digests per block (`message_cache.MessageCache`).

    ⛔ The headings are load-bearing. ``observations`` and ``geometry`` may be indexed at either end of
    a hop; a belief may be read at the SOURCE only, and the kernel reads the incoming belief once, at a
    node's OWN strand claim, before any hop.

    Every field here has a reader in the transfer policy, the kernel or the backbone.
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
    #: cross, the outside flank of a terminus, the junction's exon side (the builders,
    #: `native/transfer_kernel.h`)
    boundary_flags: np.ndarray

    # ── the solve's own scalars (neither observation nor belief) ──────────────────────────────────────
    n_grid: int
    logodds_window: float
    #: the strand protocol decision for the LIBRARY: does the spliced 2×2 read the protocol as
    #: strand-preserving (`region_init.strand_discriminability` > 0)? An unstranded verdict makes every
    #: single-strand exon's strand precision exactly zero, and a policy reading the split as an RNA
    #: witness must know that before it prices a single hop
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


@runtime_checkable
class Policy(Protocol):
    """A message-composition policy: what the kernel is told. ``name`` selects the layer the kernel runs
    (``"silent"``: none; ``"transfer"``: the composition transfer) and appears in diagnostics and in arm
    output; ``strand`` is the strand model the policy's own claims read — ``(κ, od_gdna, od_rna)`` — or
    ``None`` for no strand claim; ``library`` is its one look across the chain."""

    name: str
    strand: tuple | None

    def library(self, view: ChainView):
        """Once per sweep, over the WHOLE chain: whatever library-wide facts this policy's messages
        need, reduced from observations and geometry alone (the view carries no belief), or ``None``.
        This is the ONLY place a policy may look across the chain; the kernel sees one block."""
