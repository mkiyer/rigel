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

The local imputation (no global sweep):

* **observable region** (intron/intergenic): its own contained gDNA ``count / region_eff_len``,
  strand-cleaned.
* **non-observable region** (exon/AMBIG): the gDNA density of each *observable boundary side*
  (``crossing gDNA count / fl_mean`` — the accumulator deposits the molecule's true span, so the
  one-side crossing flux is an unbiased density estimator), averaged over the available sides.
* **run interiors** (consecutive non-observable regions, no observable side): the nearest anchored
  neighbour, carried inward from the run's observable edges (forward + reverse, averaged).
* **no anchor in the whole reference**: the global gDNA density (the count-weighted mean of the
  observable regions' densities) for the MEAN, but with **no count support** ⇒ the count prior
  collapses to Jeffreys Beta(½,½) and the strand clue governs.

Counts → density via the gDNA FL effective lengths: contained ``count ÷ region_eff_len``, crossing
``count ÷ fl_mean``. For uniform gDNA at a given density both recover that density.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_AMBIG, TS_NEG, TS_NONE

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_EPS = 1.0e-9
#: |½ − κ| below which the strand clean is degenerate (κ ≈ ½, unstranded) and cannot separate gDNA
#: from RNA — the linear unmix denominator vanishes, so the count clue keeps the raw count (frac 1)
#: and defers identifiability to whatever strand signal exists. Distinct from _EPS (a divide floor):
#: this is a "κ is effectively ½" detector on the sense-rate scale.
_UNSTRANDED_TOL = 1.0e-6


@dataclass(frozen=True, slots=True)
class NodeDensity:
    """Per-region gDNA density (count clue) after local imputation."""

    density: np.ndarray  # float64[R] — local gDNA density (fragments per effective bp)
    count_gdna_frac: (
        np.ndarray
    )  # float64[R] — count-prior MEAN: clip(density·region_eff_len / contained_mass). The
    #   strand-cleaned gDNA fraction of the contained mass; consumed by the joint deconv (the
    #   prior mean) AND the gDNA strand-fit seed weight. Separated from the concentration so the
    #   latter can carry the overdispersion-honest precision (see count_overdispersion plan).
    count_support: (
        np.ndarray
    )  # float64[R] — the HONEST gDNA count behind the density estimate: the contained gDNA count
    #   ``density·region_eff_len`` for a count-observable region, else the anchoring crossing count
    #   ``density·fl_mean`` for an imputed/carried region. ``0`` for a region with no anchor anywhere
    #   (global-baseline density, no local evidence). ``calibrate`` turns it into the count-prior
    #   concentration ``count_evidence = N/(1+α·N)`` once the count overdispersion α is fit.
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


def strand_clean_gdna_frac(
    sense: np.ndarray, total: np.ndarray, rna_sense_frac: float
) -> np.ndarray:
    """Closed-form strand-cleaned gDNA fraction of unspliced mass (the count clue's cleaner).

    A two-component linear unmix of the observed sense fraction ``s = sense/total`` between gDNA
    (symmetric, ``s = ½``) and RNA (``s = κ``)::

        gdna_frac = clip( (s − κ) / (½ − κ),  0, 1 )

    Clamped to ``[0, 1]`` — a fraction cannot be negative or exceed 1, so a noise-skewed node reads
    as 0 % gDNA (all RNA) or 100 % gDNA rather than an unphysical value. Unstranded (``κ ≈ ½``) is
    degenerate (gDNA and RNA both symmetric) → returns ``1.0`` (cannot strand-clean; keep the raw
    count). Empty nodes (``total = 0``) take ``s = ½``. The caller overrides ``TS_NONE`` (intergenic
    ⇒ pure gDNA, frac 1) and ``TS_AMBIG`` (no defined sense ⇒ defer, frac 1).
    """
    sense = np.asarray(sense, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)
    sense_frac = np.where(total > 0.0, sense / np.maximum(total, _EPS), 0.5)
    denom = 0.5 - rna_sense_frac
    if abs(denom) < _UNSTRANDED_TOL:  # unstranded — strand cannot clean
        return np.ones_like(sense_frac)
    return np.clip((sense_frac - rna_sense_frac) / denom, 0.0, 1.0)


def node_gdna_density(
    substrate,
    region_arrays,
    region_eff_len: np.ndarray,
    fl_mean: float,
    rna_sense_frac: float,
) -> NodeDensity:
    """Per-region gDNA density from the count clue via LOCAL boundary-anchored imputation.

    ``rna_sense_frac`` (κ) is the global RNA sense fraction; it strand-cleans the nascent-RNA upper
    bias from every node's unspliced count before the density is read, so the density is clean gDNA.
    See the module docstring for the estimator.
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    ts = np.asarray(region_arrays.strand_class)
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    region_count_observable, boundary_count_observable = count_observable_masks(sig, ref_id)

    # Strand-cleaned gDNA COUNT per node, oriented by region r's transcript strand. Intergenic
    # (NONE, no transcript) is pure gDNA; AMBIG (both strands, no defined sense) defers — both
    # keep the full unspliced count (frac 1), as does unstranded data (handled in the helper).
    def clean_count(view) -> np.ndarray:
        pos = np.asarray(view.n_unspliced_pos, dtype=np.float64)
        neg = np.asarray(view.n_unspliced_neg, dtype=np.float64)
        total = pos + neg
        sense = np.where(ts == TS_NEG, neg, pos)
        gf = strand_clean_gdna_frac(sense, total, rna_sense_frac)
        gf = np.where((ts == TS_NONE) | (ts == TS_AMBIG), 1.0, gf)
        return gf * total

    contained_gdna = clean_count(substrate.contained)
    left_gdna = clean_count(substrate.left)  # right side of region r's LEFT boundary
    right_gdna = clean_count(substrate.right)  # left side of region r's RIGHT boundary

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
    # baseline, not 0) for its MEAN, but carries NO count evidence (count_support 0 below) so its
    # count prior collapses to Jeffreys Beta(½,½) and the strand clue governs where present.
    anchored = ~np.isnan(density)
    own_count = contained_gdna[own]
    own_len = region_eff_len[own]
    global_density = float(np.sum(own_count) / np.sum(own_len)) if np.sum(own_len) > 0.0 else 0.0
    density = np.where(anchored, density, global_density)

    # Count-prior MEAN: the strand-cleaned gDNA fraction of the contained mass,
    # clip(density·eff_len / contained_mass). Computed here (consumed by the joint deconv as the
    # prior mean AND by the gDNA strand-fit seed weight) so the concentration is free to carry the
    # overdispersion-honest precision (count_evidence = N/(1+α·N)) downstream.
    contained_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(
            np.where(contained_mass > 0.0, density * region_eff_len / contained_mass, 0.0), 0.0, 1.0
        )

    # Count-prior SUPPORT: the honest gDNA count behind the estimate. A count-observable region uses
    # its directly-observed contained gDNA count ``density·region_eff_len``; an imputed/carried
    # region uses the anchoring crossing count ``density·fl_mean`` (its eff length is the gDNA-FL
    # crossing mean). A no-anchor region carries 0 (global baseline, no local evidence). calibrate
    # turns this into the overdispersion-limited concentration N/(1+α·N) — replacing the old
    # categorical defer_to_strand zeroing with one principled precision.
    eff_type = np.where(region_count_observable, region_eff_len, fl_mean)
    count_support = np.where(anchored, density * eff_type, 0.0)

    return NodeDensity(
        density=density,
        count_gdna_frac=count_gdna_frac,
        count_support=count_support,
        region_count_observable=region_count_observable,
        boundary_count_observable=boundary_count_observable,
        n_region_count_observable=int(region_count_observable.sum()),
        n_boundary_count_observable=int(boundary_count_observable.sum()),
    )


__all__ = [
    "NodeDensity",
    "count_observable_masks",
    "strand_clean_gdna_frac",
    "node_gdna_density",
]
