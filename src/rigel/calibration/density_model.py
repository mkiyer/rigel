"""Phase 1 — the density model ("count clue"): per-region gDNA density from OBSERVED counts.

Acyclic by construction: the gDNA density is read **directly** from count-observable nodes
(where fragments are gDNA by construction) and **imputed locally** for the rest. It never consults
a global ``gdna_density_global * L`` product, so there is no density->deconv->density feedback loop.

Count-observability is a property of the region **signature** (4-bit exon/intron ± flags):

* **region** is observable ⇔ it has **no exon bit** (intergenic or intron-only). Its
  contained unspliced mass is gDNA (+ nascent RNA — an upper bias the strand clue removes);
  an exonic region's contained mass is contaminated by mature RNA.
* **boundary** is observable ⇔ **no exon bit is shared** across its two sides → no single
  exon-strand continues across it → no *unspliced mature RNA* crosses (a single-exon
  transcript spanning the seam would put unspliced mature RNA there). Its crossing
  unspliced mass is then gDNA(+nascent).

Raw counts (no strand cleaning): in the decoupled architecture the strand module owns the strand
channel; this count module is the fallback for strand-unobservable / unstranded nodes and works on
the raw unspliced count. The local imputation (no global sweep):

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

from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class NodeDensity:
    """Per-region gDNA density (count clue) after local imputation."""

    density: np.ndarray  # float64[R] — local gDNA density (fragments per effective bp)
    count_gdna_frac: (
        np.ndarray
    )  # float64[R] — count module's gDNA fraction g_count = clip(density·region_eff_len /
    #   contained_mass): the gDNA fraction of the contained unspliced mass from the (raw) local gDNA
    #   density. Consumed by the count module (strand-unobservable / unstranded nodes) and the gDNA
    #   strand-fit seed weight.
    region_count_observable: np.ndarray  # bool[R] — count-observable region (non-exonic)
    boundary_count_observable: np.ndarray  # bool[R] — count-observable boundary right of region r
    n_region_count_observable: int
    n_boundary_count_observable: int


def count_observable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-observability for regions and (right-) boundaries.

    Returns ``(region_count_observable, boundary_count_observable)``, both ``bool[R]``. ``boundary_count_observable[r]``
    describes the internal boundary between region ``r`` and ``r+1`` (defined iff same ref).
    """
    sig = np.asarray(signature).astype(np.int64)
    ref = np.asarray(ref_id)
    r = sig.shape[0]
    region_count_observable = (sig & _EXON_BITS) == 0
    boundary_count_observable = np.zeros(r, dtype=bool)
    if r > 1:
        same = ref[:-1] == ref[1:]
        shared_exon = (sig[:-1] & sig[1:] & _EXON_BITS) != 0
        boundary_count_observable[:-1] = same & ~shared_exon
    return region_count_observable, boundary_count_observable


def node_gdna_density(
    substrate,
    region_arrays,
    region_eff_len: np.ndarray,
    fl_mean: float,
) -> NodeDensity:
    """Per-region gDNA density from the count clue via LOCAL boundary-anchored imputation.

    The count module estimates gDNA density from **raw** unspliced counts — no strand cleaning (in
    the decoupled architecture the strand module owns the strand channel; the count module is the
    fallback for strand-unobservable / unstranded nodes). The density is read directly from
    count-observable regions and imputed locally for the rest; a region with no local anchor anywhere
    takes the count-weighted-mean observable density as a global fallback. See the module docstring.
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    region_count_observable, boundary_count_observable = count_observable_masks(sig, ref_id)

    # Raw unspliced COUNT per node (no strand cleaning — the count module works on raw counts).
    def total_count(view) -> np.ndarray:
        pos = np.asarray(view.n_unspliced_pos, dtype=np.float64)
        neg = np.asarray(view.n_unspliced_neg, dtype=np.float64)
        return pos + neg

    contained_gdna = total_count(substrate.contained)
    left_gdna = total_count(substrate.left)  # right side of region r's LEFT boundary
    right_gdna = total_count(substrate.right)  # left side of region r's RIGHT boundary

    # Per-side boundary observability for region r: its LEFT side uses boundary (r−1, r); its RIGHT
    # side uses boundary (r, r+1). boundary_count_observable[k] describes boundary (k, k+1).
    left_anchor = np.zeros(r, dtype=bool)
    right_anchor = np.zeros(r, dtype=bool)
    if r > 1:
        left_anchor[1:] = boundary_count_observable[:-1] & (ref_id[1:] == ref_id[:-1])
        right_anchor[:-1] = boundary_count_observable[:-1] & (ref_id[:-1] == ref_id[1:])

    density = np.full(r, np.nan, dtype=np.float64)
    # Observable region with a usable contained length → its own contained density.
    own = region_count_observable & (region_eff_len > _EPS)
    density[own] = contained_gdna[own] / region_eff_len[own]
    # Everything else: anchor from the available observable boundary sides (crossing count / fl_mean).
    inv_fl = 1.0 / fl_mean if fl_mean > 0.0 else 0.0
    for i in np.where(~own)[0]:
        est = []
        if left_anchor[i]:
            est.append(left_gdna[i] * inv_fl)
        if right_anchor[i]:
            est.append(right_gdna[i] * inv_fl)
        if est:
            density[i] = float(np.mean(est))

    # Run-fill: a region still unset (a run interior with no observable side) inherits the nearest
    # anchored neighbour, carried inward from both directions within its reference and averaged.
    fwd = density.copy()
    for i in range(1, r):
        if np.isnan(fwd[i]) and ref_id[i] == ref_id[i - 1]:
            fwd[i] = fwd[i - 1]
    rev = density.copy()
    for i in range(r - 2, -1, -1):
        if np.isnan(rev[i]) and ref_id[i] == ref_id[i + 1]:
            rev[i] = rev[i + 1]
    stack = np.vstack([fwd, rev])
    valid = ~np.isnan(stack)
    n_valid = valid.sum(axis=0)
    carried = np.where(n_valid > 0, np.nansum(stack, axis=0) / np.maximum(n_valid, 1), np.nan)
    density = np.where(np.isnan(density), carried, density)

    # A region still unset has NO local anchor anywhere in its reference. It takes the GLOBAL gDNA
    # density (the count-weighted mean of the count-observable regions' own densities — a sensible
    # baseline, not 0). A zero-gDNA library has ~no observable mass ⇒ baseline ≈ 0 ⇒ no manufactured
    # gDNA. (A rare fallback — local imputation handles the rest.)
    anchored = ~np.isnan(density)
    own_len = float(np.sum(region_eff_len[own]))
    global_density = float(np.sum(contained_gdna[own]) / own_len) if own_len > 0.0 else 0.0
    density = np.where(anchored, density, global_density)

    # Count-prior MEAN g_count = clip(density·eff_len / contained_mass): the gDNA fraction of the
    # contained unspliced mass implied by the local gDNA density. Consumed by the count module
    # (strand-unobservable / unstranded nodes) and by the gDNA strand-fit seed weight.
    contained_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(
            np.where(contained_mass > 0.0, density * region_eff_len / contained_mass, 0.0), 0.0, 1.0
        )

    return NodeDensity(
        density=density,
        count_gdna_frac=count_gdna_frac,
        region_count_observable=region_count_observable,
        boundary_count_observable=boundary_count_observable,
        n_region_count_observable=int(region_count_observable.sum()),
        n_boundary_count_observable=int(boundary_count_observable.sum()),
    )


__all__ = [
    "NodeDensity",
    "count_observable_masks",
    "node_gdna_density",
]
