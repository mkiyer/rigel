"""The message-policy interface for the native two-pass solve, and the chain view a policy reads.

The backbone (``sweep.solve_chain``) owns traversal and the solve. A policy supplies a name the kernel switches on
(``"transfer"`` or ``"silent"``) and a library, ``policy.library(view)``, reduced once per sweep from observations and
geometry. A message builder may read counts and geometry at either end of a hop but never a posterior belief: own
strand claims are conditional-binomial likelihoods of the observed columns, and counting noise or observed
disagreement may widen a received claim but never move its mode toward a belief. :class:`ChainView` holds the whole
chain for the library (built from :func:`view_fields`); the native ``transfer_kernel::Chain`` holds one block for the
builders. Neither has a belief field.

Gate: ``tests/calibration/test_sweep_backbone.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..region_chain import REGION
from ..region_geometry import g1_locked

__all__ = [
    "ChainView",
    "view_fields",
]


@dataclass(frozen=True, slots=True)
class ChainView:
    """Observations and geometry available to a policy, without posterior beliefs.

    A policy's ``library`` sees the whole chain; native message builders read the same
    fields a block at a time. Counts and geometry may be read at either end of a hop.
    """

    # ── observations: readable at either end of a hop ──────────────────────────────────────────────────────────
    eff_gdna: np.ndarray  # per-slot gDNA opportunity (the capture-blind effective length)
    eff_rna: np.ndarray  # per-slot RNA opportunity
    sj_count: np.ndarray  # [n, 2] sj fragment count by transcript strand (both faces)
    #: the same flux split by which genomic end of its sj this boundary is: the count matching
    #: ``route_rate_lo`` / ``route_rate_hi``, so a per-face rate is priced on its own count
    sj_count_lo: np.ndarray
    sj_count_hi: np.ndarray
    #: [n, 2] the route-summed certified rate per face by transcript strand (Σ flux_J/A_J over the face's disjoint
    #: routes)
    route_rate_lo: np.ndarray
    route_rate_hi: np.ndarray
    unspliced_count: (
        np.ndarray
    )  # [n, 2] unspliced count by genome strand: the density numerator and the Poisson n
    spliced_count: np.ndarray  # [n, 2] spliced count by strand

    # ── geometry and structure: readable at either end, belief-free by construction ───────────────────────────
    left: np.ndarray  # adjacent slot of the other kind, -1 at a reference start
    right: np.ndarray
    is_boundary: np.ndarray  # ~is_region; the chain strictly alternates N E N E … N
    is_exon_region: (
        np.ndarray
    )  # a region whose signature is exon: the splice-in message's destination
    free_pos: np.ndarray  # does the annotation admit RNA+ here? (one of the two population bits)
    free_neg: np.ndarray  # …and RNA−?
    #: the region signature's two exon bits per slot (False at a boundary): a strand's opportunity geometry. A region
    #: that admits strand ``s`` and carries no exon of ``s`` is ``s``'s intron whatever the other strand does there,
    #: which the coarse ``is_exon_region`` cannot say and the RNA level lanes need
    exon_pos: np.ndarray
    exon_neg: np.ndarray
    #: the terminus and junction bits per boundary slot (0 at a region): which faces composition may cross, the outside
    #: flank of a terminus, the junction's exon side (read by the builders in ``native/transfer_kernel.h``)
    boundary_flags: np.ndarray

    # ── the library's strand verdict (neither observation nor belief) ─────────────────────────────────────────────
    #: the library's strand protocol decision (``strand_balance.strand_discriminability`` > 0). When it is False every
    #: single-strand exon's strand precision is exactly zero, which a policy reading the split as an RNA witness must
    #: know before it prices a hop
    strand_live: bool = False

    @property
    def terminal(self) -> np.ndarray:
        """A region that admits no RNA strand: intergenic, structurally pure gDNA (``region_geometry.g1_locked``'s
        class). The one definition: the sweep cuts the chain into loci at these and never delivers into one, the
        library reads its gDNA coordinate off them, and the kernel's builders take this array as their intergenic
        class."""
        return ~np.asarray(self.is_boundary, bool) & g1_locked(self.free_pos, self.free_neg)

    @property
    def n_slots(self) -> int:
        return int(self.unspliced_count.shape[0])

    @property
    def n_slot(self) -> np.ndarray:
        """The per-slot unspliced count over both strands: the density numerator and the Poisson n."""
        return self.unspliced_count.sum(axis=1)


def view_fields(chain, statics, geometry, structure) -> dict:
    """Every per-slot array a policy may read, under the two belief-free headings of :class:`ChainView`, for the whole
    chain; the kernel reads it a block at a time. ``structure`` is the sweep's per-slot structural classes."""
    return dict(
        # observations
        eff_gdna=np.asarray(geometry.eff_gdna, np.float64),
        eff_rna=np.asarray(geometry.eff_rna, np.float64),
        sj_count=np.asarray(geometry.sj_count, np.float64),  # [n, 2] by transcript strand
        sj_count_lo=np.asarray(geometry.sj_count_lo, np.float64),
        sj_count_hi=np.asarray(geometry.sj_count_hi, np.float64),
        route_rate_lo=np.asarray(geometry.route_rate_lo, np.float64),
        route_rate_hi=np.asarray(geometry.route_rate_hi, np.float64),
        unspliced_count=np.asarray(geometry.unspliced_count, np.float64),  # [n, 2] by genome strand
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
