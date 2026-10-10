"""Message policies supplied to the native two-pass solve.

The backbone (`sweep.solve_chain`) owns traversal and count inference. A policy
selects the message layer, its RNA sense fraction and its library coordinates:

    library = policy.library(view)  # once per sweep, from observations and geometry
    policy.name                     # "transfer" or "silent"
    policy.kappa                    # RNA sense fraction, or None for no strand claims

Transfer carries composition profiles or component levels across the existing
chain. Each forward/backward pass combines a sender's own observations with its
far-side message; the solve then reads both received messages and the density
prior. Silent sends nothing and supplies the comparison floor.

A message builder may read counts and geometry at either end of a hop, but no
posterior beliefs. Own strand claims are conditional-binomial likelihoods of the
observed columns. Counting noise and observed disagreement can widen a received
claim; they cannot manufacture agreement by moving its mode toward a belief.

`ChainView` provides the whole chain for library reductions. The native
`transfer_kernel::Chain` provides one block to message builders. Neither has a
belief field. The count solver has its own separate belief inputs.

Gate: `tests/calibration/test_sweep_backbone.py`.
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
    """Observations and geometry available to a policy, without posterior beliefs.

    `Policy.library` sees the whole chain; native message builders read the same
    fields a block at a time. Counts and geometry may be read at either end of a hop.
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
    #: the face's disjoint routes)
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

    # ── the library's strand verdict (neither observation nor belief) ─────────────────────────────────
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


@runtime_checkable
class Policy(Protocol):
    """The message layer, observed-strand protocol and library reduction."""

    name: str
    kappa: float | None

    def library(self, view: ChainView):
        """Once per sweep, over the WHOLE chain: whatever library-wide facts this policy's messages
        need, reduced from observations and geometry alone (the view carries no belief), or ``None``.
        This is the ONLY place a policy may look across the chain; the kernel sees one block."""
