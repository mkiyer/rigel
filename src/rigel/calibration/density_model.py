"""Per-region gDNA density from observed counts — the seed for the gDNA strand-overdispersion fit.

The gDNA density is read **directly** from count-observable regions (where fragments are gDNA by
construction) and **imputed locally** for the rest. Its output (``count_gdna_frac`` + the
count-observability masks) is consumed only by the gDNA strand-overdispersion fit (``gdna_strand``) as
its seed selector and by the global prior as the signature partition — NOT by the per-region
deconvolution itself (the belief-propagation sweep owns that). It never consults a global
``gdna_density_global * L`` product, so there is no density->deconv->density feedback loop.

Count-observability is a property of the region **signature** (4-bit exon/intron ± flags):

* **region** is observable ⇔ it has **no exon bit** (intergenic or intron-only). Its
  contained unspliced mass is gDNA (+ nascent RNA — an upper bias the strand clue removes);
  an exonic region's contained mass is contaminated by mature RNA.
* **boundary** is observable ⇔ **no exon bit is shared** across its two sides → no single
  exon-strand continues across it → no *unspliced mature RNA* crosses (a single-exon
  transcript spanning the seam would put unspliced mature RNA there). Its crossing
  unspliced mass is then gDNA(+nascent).

Raw counts (no strand cleaning): the strand channel is owned by the belief-propagation sweep; this
module works on the raw unspliced count purely to seed the strand-overdispersion fit. The local
imputation (no global sweep):

* **observable region** (intron/intergenic): its own contained ``count / region_eff_len``.
* **non-observable region** (exon/AMBIG): the gDNA density of each *observable boundary side*
  (``crossing count / fl_mean`` — the accumulator deposits the molecule's true span, so the
  one-side crossing flux is an unbiased density estimator), averaged over the available sides.
* **run interiors** (consecutive non-observable regions, no observable side): the nearest anchored
  neighbour, carried inward from the run's observable edges (forward + reverse, averaged).
* **no anchor in the whole reference**: the global gDNA density (the count-weighted mean of the
  count-observable regions' densities) — a rare fallback.

Counts → density via the gDNA FL effective lengths: contained ``count ÷ region_eff_len``, crossing
``count ÷ fl_mean``. For uniform gDNA at a given density both recover that density.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .region_chain import REGION
from .region_arrays import edge_region_indices
from .run_fill import runfill_bidirectional
from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class RegionDensity:
    """Per-region gDNA density (count clue) after local imputation."""

    density: np.ndarray  # float64[N] — local gDNA density (fragments per effective bp)
    count_gdna_frac: (
        np.ndarray
    )  # float64[N] — count module's gDNA fraction g_count = clip(density·eff_gdna /
    #   contained count): the gDNA fraction of the contained unspliced mass from the (raw) local gDNA
    #   density. KEPT as the gDNA strand-overdispersion fit SEED selector (gdna_strand.py) — NOT a
    #   gDNA-fraction vote in the solve (the BP sweep owns the gDNA/RNA call).
    region_count_observable: np.ndarray  # bool[N] — count-observable region (non-exonic)
    edge_count_observable: np.ndarray  # bool[E] — count-observable contiguous edge


def count_observable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-observability, on the two axes it describes.

    Returns ``(region_count_observable[N], edge_count_observable[E])``.

    ⭐ **The edge mask is on the EDGE axis now.** It used to be a ``bool[R]`` in which entry ``r``
    described the seam ``(r, r+1)`` and the last entry of every reference was a padding ``False`` — a
    region-shaped array standing in for an edge-shaped one, with the ``same``-reference test carried
    inside it. A contiguous edge is a first-class object with its own axis, so the mask is
    ``bool[E]`` and the ``same`` test is :func:`~rigel.calibration.region_arrays.edge_region_indices`'s
    job rather than this function's.
    """
    sig = np.asarray(signature).astype(np.int64)
    region_count_observable = (sig & _EXON_BITS) == 0
    lo, hi = edge_region_indices(ref_id)
    # observable ⇔ no exon bit is SHARED across the line ⇒ no single exon-strand continues across it
    # ⇒ no unspliced mature RNA crosses, so the crossing count is gDNA (+ nascent).
    edge_count_observable = (sig[lo] & sig[hi] & _EXON_BITS) == 0
    return region_count_observable, edge_count_observable


def region_gdna_density(chain, geometry, region_arrays) -> RegionDensity:
    """Per-region gDNA density from the count clue via LOCAL edge-anchored imputation.

    Density is read directly from count-observable regions and imputed locally for the rest; a region with
    no local anchor anywhere takes the count-weighted-mean observable density as a global fallback. See
    the module docstring.

    ⭐ **Two arguments dissolved with the faces (S5.e).** ``region_eff_len`` and ``fl_mean`` were passed
    in separately; both are now :class:`~rigel.calibration.region_geometry.RegionGeometry`'s own
    ``eff_gdna`` — the contained placements at a REGION slot, the crossing placements at an EDGE slot.
    One divisor rule, one source, so a caller can no longer hand this a length model that disagrees
    with the one the solver uses. The flanking edges are found through the chain's own adjacency, so
    the ``left``/``right`` per-region side views are gone too.

    Output ``count_gdna_frac`` (the density→fraction MEAN) is consumed only as the gDNA
    strand-overdispersion fit SEED selector (``gdna_strand.py``) and the signature-partition masks.
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    region_count_observable, edge_count_observable = count_observable_masks(sig, ref_id)

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    slot_of_region = np.zeros(r, dtype=np.int64)
    slot_of_region[obj[is_region]] = np.flatnonzero(is_region)

    count = np.asarray(geometry.unspliced_count, dtype=np.float64).sum(axis=1)
    eff = np.asarray(geometry.eff_gdna, dtype=np.float64)

    contained_gdna = count[slot_of_region]
    region_eff_len = eff[slot_of_region]

    def flank(neighbour: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The gDNA density at a region's flanking EDGE, and whether that edge can anchor it.

        ⭐ **The edge is reached through the chain's own adjacency**, and its observability is read on
        the EDGE axis at that slot's ``obj_idx`` — so "which seam is to my left" is answered by the
        chain rather than by region-index arithmetic that has to know about reference terminals.

        ⚠ An object with no opportunity contributes nothing rather than a floored division
        so a zero divisor makes the flank unusable, not infinite.
        """
        slots = neighbour[slot_of_region]
        exists = slots >= 0
        safe = np.clip(slots, 0, max(count.shape[0] - 1, 0))
        e = np.where(exists, eff[safe], 0.0)
        usable = exists & (e > _EPS)
        edge_idx = np.clip(obj[safe], 0, max(edge_count_observable.shape[0] - 1, 0))
        observable = (
            usable & edge_count_observable[edge_idx]
            if edge_count_observable.size
            else np.zeros(r, dtype=bool)
        )
        return np.where(usable, count[safe] / np.where(usable, e, 1.0), 0.0), observable

    left_gdna, left_anchor = flank(np.asarray(chain.left, dtype=np.int64))
    right_gdna, right_anchor = flank(np.asarray(chain.right, dtype=np.int64))

    density = np.full(r, np.nan, dtype=np.float64)
    # Observable region with a usable contained length → its own contained density. (Exons are NOT
    # count-observable and are imputed from the edges below; the strand for an exon enters the
    # deconvolution as ``g_strand`` in the combine, not via the count density — so this stays the
    # signature count-observable set, and ``g_count`` carries count magnitude only, no double-count.)
    own = region_count_observable & (region_eff_len > _EPS)
    density[own] = contained_gdna[own] / region_eff_len[own]
    # Everything else: anchor from the available observable edges — the AVERAGE of whichever sides are
    # anchored (both → mean, one → that side, none → stays NaN for the run-fill below).
    todo = ~own
    la, ra = todo & left_anchor, todo & right_anchor
    n_sides = la.astype(np.float64) + ra.astype(np.float64)
    side_sum = np.where(la, left_gdna, 0.0) + np.where(ra, right_gdna, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        density = np.where(n_sides > 0.0, side_sum / np.maximum(n_sides, 1.0), density)

    # Run-fill: a region still unset (a run interior with no observable side) inherits the nearest
    # anchored neighbour, carried inward from both directions within its reference and averaged.
    density = runfill_bidirectional(density, ref_id)

    # A region still unset has NO local anchor anywhere in its reference. It takes the GLOBAL gDNA
    # density (the count-weighted mean of the count-observable regions' own densities — a sensible
    # baseline, not 0). A zero-gDNA library has ~no observable mass ⇒ baseline ≈ 0 ⇒ no manufactured
    # gDNA. (A rare fallback — local imputation handles the rest.)
    anchored = ~np.isnan(density)
    own_len = float(np.sum(region_eff_len[own]))
    global_density = float(np.sum(contained_gdna[own]) / own_len) if own_len > 0.0 else 0.0
    density = np.where(anchored, density, global_density)

    # Count-prior MEAN g_count = clip(density·eff_len / contained count): the gDNA fraction of the
    # contained unspliced population implied by the local gDNA density. Consumed as the gDNA strand-fit
    # seed selector (gdna_strand.py) — NOT a gDNA-fraction vote in the solve (the count prior was removed).
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(
            np.where(contained_gdna > 0.0, density * region_eff_len / contained_gdna, 0.0), 0.0, 1.0
        )

    return RegionDensity(
        density=density,
        count_gdna_frac=count_gdna_frac,
        region_count_observable=region_count_observable,
        edge_count_observable=edge_count_observable,
    )


__all__ = [
    "RegionDensity",
    "count_observable_masks",
    "region_gdna_density",
]
