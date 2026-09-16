"""rigel.calibration.abundance_landscape — the pre-pass-0 total-density field, and its mode census.

The one question this module answers: what does the library's total fragment density look like over the
genome, before anything is solved, and which regions sit on which mode? Under hybrid capture the field
is bimodal by construction — a depleted off-target level and an enriched on-target one, a couple of
decades apart — and off capture it is unimodal. The census reads that structure off a fitted density:
`rho_0`, the depleted mode, which the pooled intergenic anchors also measure independently; the span
`R`, the ratio between the two modes, read off the field and never from in-gene anchors, which
under-read the enriched level; and per region a responsibility `w_i` for the enriched basin.

The estimator is `landscape.fit_landscape`, reused as it stands. It is deliberately component-agnostic
and every decision in it transfers: zero-native Poisson kernels (a wall-exact region that sequenced
nothing says "below the resolution wall", not "at 1/E"), knn population resolution (what suppresses
combing with no tuning), the one-pseudo-region Laplace floor, and the grid derived from the data's own
support. What is particular to this module is the inputs and the census:

* the inputs are the measured totals — `total_abundance.region_counts_and_exposure`'s side-selected
  START/END counts over the region's own length, REGIONs only (boundaries cross rather than contain and
  have zero genomic measure, the same geometry rule the gDNA hyperprior's substrate follows), restricted
  to the wall-exact (`model_free`) population;
* `var = 0` — a direct measurement has no deconvolution ambiguity, so the reliability weights are
  honestly flat. The hyperprior's `Var(log f_g)` weight exists because its training data came out of a
  solve; this data did not.

⛔ No significance threshold exists in the census, anywhere. Every interior local maximum is a mode, the
grid is partitioned into basins at the minima between them, and masses carry every verdict continuously.
A phantom wiggle above the bulk owns near-zero basin mass and so yields `w ≈ 0`, which is harmless, and
the capture-OFF unimodality gate measures that rather than a constant asserting it. The depleted mode is
picked by an independent measurement — the basin containing the pooled intergenic anchor rate — and the
anchor-consistency verdict's tolerance is the depleted mode's own fitted width, the density's statement
of its own resolution, never a chosen number.

The field conflates enrichment with expression, and that is stated up front: a hot unprobed exon and a
probed cold one can land in the same basin, so the failure direction of any consumer is permissive (an
over-wide enriched basin) and never a hard exclusion.

`_KNN_SCALE` / `_S0` are inherited from `landscape` and have only ever been validated on gDNA-shaped
data — that module's own warning — so every result reported off this landscape carries the caveat until
it is priced.

Nothing here decides anything in the solve. Its consumers are the QC surface
(`CalibrationDiagnostics.from_abundance_landscape`) and the injection substrate. Its own falsification
is ``tests/calibration/test_abundance_landscape.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .landscape import (
    _KNN_SCALE,
    _LOCATED_VAR,
    DensityLandscape,
    _poisson_kernels,
    fit_landscape,
    knn_widths,
)
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS
from .total_abundance import RegionWallMask, region_counts_and_exposure

_EPS = 1e-12
_GENE_BITS = BIT_INTRON_POS | BIT_INTRON_NEG | BIT_EXON_POS | BIT_EXON_NEG


__all__ = [
    "AbundanceLandscape",
    "AbundanceMode",
    "LocatedMode",
    "fit_abundance_landscape",
    "located_enriched_mode",
    "split_basins",
]


@dataclass(frozen=True, slots=True)
class AbundanceMode:
    """One mode of the fitted total-density field: a local maximum and its basin.

    ``log_rho`` is the peak's location (natural-log rate); ``lo``/``hi`` the basin bounds (the
    interior minima flanking it — adjacent basins share a bound, so the modes partition the grid);
    ``basin_mass`` the share of the fitted density inside the basin; ``width`` the mass-weighted
    standard deviation of ``log_rho`` within the basin — the fit's own statement of how precisely this
    mode is located, and the only tolerance any consumer is given.
    """

    log_rho: float
    basin_mass: float
    width: float
    lo: float
    hi: float


@dataclass(frozen=True, slots=True)
class AbundanceLandscape:
    """The fitted total-density field plus its census. See the module docstring for every rule.

    ``w_slot`` is per-REGION (the accumulator's region axis, not the chain): the region's own Poisson
    kernel times the fitted density, normalised, integrated over the ENRICHED basin — an honest
    posterior responsibility. ``0`` everywhere when the field is unimodal; ``NaN`` where the region is
    not model-free (a double-walled region has no trustworthy total and therefore no reading).

    ``anchor_log_rho`` is the independent depleted-level estimator — the pooled rate over intergenic
    model-free regions (the same composition-free pool ``fit_intron_background`` uses) — and
    ``anchor_consistent`` says whether it falls within the depleted mode's own width. ``NaN`` (and a
    ``False``-free fallback to the largest basin) when no intergenic region exists, which is a toy.

    ``train_log_rho`` / ``train_class`` are the population this landscape was fitted on: one entry per
    selected region, the kernel centre ``log(max(count,1)) − log(exposure)`` in natural log, and its
    coarse class (``0`` intergenic / ``1`` intron / ``2`` exon — the report's own rug codes, where ``3``
    is a boundary and this substrate has none, being REGIONs only). They are published because a
    consumer that plots the fit needs to plot what it was fitted on, and re-deriving the centres
    elsewhere would be a second copy of the selection rule.
    """

    landscape: DensityLandscape
    modes: tuple[AbundanceMode, ...]
    depleted: AbundanceMode
    enriched: AbundanceMode | None
    rho_0: float
    span_R: float
    anchor_log_rho: float
    anchor_gap_nats: float
    anchor_consistent: bool
    w_slot: np.ndarray
    n_train: int
    train_log_rho: np.ndarray
    train_class: np.ndarray


def _census(landscape: DensityLandscape) -> tuple[AbundanceMode, ...]:
    """Every interior local maximum with its basin. Basins split at the interior minima between
    adjacent maxima, so they partition the grid and their masses sum to one — asserted by a gate,
    never renormalised here."""
    x = np.asarray(landscape.log_rho, dtype=np.float64)
    p = np.exp(np.asarray(landscape.logP, dtype=np.float64))
    p = p / max(float(p.sum()), _EPS)
    interior = np.where((p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:]))[0] + 1
    if interior.size == 0:
        interior = np.array([int(np.argmax(p))])
    peaks = np.sort(interior)
    # basin bounds: the minimum of the density between each adjacent pair of peaks
    cuts = [0]
    for a, b in zip(peaks[:-1], peaks[1:], strict=False):
        cuts.append(int(a + np.argmin(p[a : b + 1])))
    cuts.append(p.size - 1)
    modes = []
    for i, pk in enumerate(peaks):
        s, e = cuts[i], cuts[i + 1]
        # half-open segments so the shared cut bin is counted once — the last basin takes the final
        # grid point. Without this the basin masses sum to 1 + (the mass at each cut), which is what the
        # partition gate catches.
        hi_idx = e + 1 if i == len(peaks) - 1 else e
        seg_p, seg_x = p[s:hi_idx], x[s:hi_idx]
        m = float(seg_p.sum())
        mu = float((seg_p * seg_x).sum() / max(m, _EPS))
        width = float(np.sqrt((seg_p * (seg_x - mu) ** 2).sum() / max(m, _EPS)))
        modes.append(
            AbundanceMode(
                log_rho=float(x[pk]),
                basin_mass=m,
                width=width,
                lo=float(x[s]),
                hi=float(x[e]),
            )
        )
    # shared bounds: mode i's hi is mode i+1's lo, both cut at one index. The construction guarantees
    # it; the partition gate asserts it.
    return tuple(modes)


def split_basins(
    modes: tuple[AbundanceMode, ...], anchor_log_rho: float
) -> tuple[AbundanceMode, AbundanceMode | None]:
    """The depleted/enriched selection rule, kept on its own so a scorer can apply it to another
    estimator's census: depleted is the basin containing the anchor rate (the nearest mode if the anchor
    falls between basins, the largest-mass basin when the anchor is NaN, which happens only on a toy);
    enriched is the largest-mass basin strictly above the depleted one, ``None`` meaning unimodal."""
    if np.isfinite(anchor_log_rho):
        inside = [m for m in modes if m.lo <= anchor_log_rho <= m.hi]
        depleted = (
            inside[0] if inside else min(modes, key=lambda m: abs(m.log_rho - anchor_log_rho))
        )
    else:
        depleted = max(modes, key=lambda m: m.basin_mass)
    above = [m for m in modes if m.lo >= depleted.hi - _EPS and m is not depleted]
    enriched = max(above, key=lambda m: m.basin_mass) if above else None
    return depleted, enriched


@dataclass(frozen=True, slots=True)
class LocatedMode:
    """A mode of the located population: the basin, and the number of located kernels behind it — the
    regime a consumer publishes beside the reference it reads off ``mode.log_rho``."""

    mode: AbundanceMode
    n_members: int


def located_enriched_mode(landscape: DensityLandscape) -> LocatedMode | None:
    """The mode a gDNA-density consumer may take as the fully-captured level, or ``None``.

    The census names the depleted basin as the largest by rendered mass (for gDNA the unprobed
    regions always outnumber the probed ones, and the zero-count anchors are that population's own
    statement). Above it the candidate is the basin holding the most LOCATED kernels — a kernel with a
    location, ``landscape.located``, is one that counted at least a fragment; an anchor's or a
    sub-fragment kernel's centre is its resolution wall ``1/E`` and is no member of anything. A human
    index trains a quarter of a million anchors whose walls span every decade, and a basin above the
    bulk can be packed with them around ten measured kernels (`ISSUES: the-ruler-reference-on-sparse-real-libraries`).

    The candidate is a MODE only if its members resolve it at the located population's own
    resolution: with ``k = √n_located`` (:func:`~.landscape.knn_widths`' population ``k``), each
    member's width is half the distance to its ``k``-th nearest MEMBER, and the median of those
    satisfies ``width² ≤`` :data:`~.landscape._LOCATED_VAR`, the one-fragment floor in the population's
    own variable. A basin with ``k`` members or fewer has no ``k``-th neighbour inside itself — the
    cluster smaller than ``√n`` that reaches outside itself — and is no mode, however narrow the
    rendered density's cut made it; the within-basin spread is NOT the statement, since a basin cut by
    the grid's edge is narrow whatever its kernels. ``None`` is "no enriched mode" — the capture-OFF
    field, the gDNA-free field, and a library whose gDNA is too sparse to locate its probed level —
    and the consumer then contracts nothing (`capture_eff_length`, `priors`). The verdict and the peak
    are stable across an 8× range of the render resolution
    (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`).
    """
    modes = _census(landscape)
    depleted, _enriched = split_basins(modes, float("nan"))
    above = [m for m in modes if m.lo >= depleted.hi - _EPS and m is not depleted]
    if not above:
        return None
    centre = np.asarray(landscape.centre, dtype=np.float64)
    located = np.asarray(landscape.located, dtype=bool)
    n_located = int(located.sum())
    k = max(int(round(np.sqrt(n_located))), 2)

    def members(m: AbundanceMode) -> np.ndarray:
        return located & (centre >= m.lo) & (centre <= m.hi)

    enriched = max(above, key=lambda m: int(members(m).sum()))
    mem = members(enriched)
    n_members = int(mem.sum())
    if n_members <= k:
        return None
    step = float(landscape.log_rho[1] - landscape.log_rho[0])
    width = knn_widths(centre[mem], step, k=k)
    if float(np.median(width)) ** 2 > _LOCATED_VAR:
        return None
    return LocatedMode(mode=enriched, n_members=n_members)


def fit_abundance_landscape(
    substrate, region_arrays, wall_mask: RegionWallMask, *, knn_scale: float = _KNN_SCALE
) -> AbundanceLandscape | None:
    """Fit the total-density field on the wall-exact measured totals, and read the census off it.

    Returns ``None`` when fewer than two model-free regions exist (a fit over one point is not a
    population), exactly as the underlying estimator does.
    """
    counts, exposure, model_free = region_counts_and_exposure(substrate, region_arrays, wall_mask)
    sel = np.asarray(model_free, dtype=bool) & (exposure > 0.0)
    n_regions = counts.shape[0]
    if int(sel.sum()) < 2:
        return None

    c, e = counts[sel], exposure[sel]
    landscape = fit_landscape(
        c,
        c,  # a TOTAL's ceiling is the observed density itself: mass ≡ count fixes the grid top
        e,
        np.zeros_like(c),  # a direct measurement has no deconvolution ambiguity
        anchor=(c == 0.0),
        knn_scale=knn_scale,
    )
    if landscape is None:
        return None

    modes = _census(landscape)

    # ── the independent depleted-level estimator: the pooled intergenic rate, composition-free
    sig = np.asarray(region_arrays.signature, dtype=np.int64)
    anchors = sel & ((sig & _GENE_BITS) == 0)
    if anchors.any() and float(exposure[anchors].sum()) > 0.0 and counts[anchors].sum() > 0.0:
        anchor_log_rho = float(np.log(counts[anchors].sum()) - np.log(exposure[anchors].sum()))
    else:
        anchor_log_rho = float("nan")

    # ── depleted/enriched: one rule with one home — `split_basins`
    depleted, enriched = split_basins(modes, anchor_log_rho)

    span_R = float(np.exp(enriched.log_rho - depleted.log_rho)) if enriched is not None else 1.0
    gap = abs(depleted.log_rho - anchor_log_rho) if np.isfinite(anchor_log_rho) else float("nan")
    # the tolerance is the depleted mode's own fitted width — the density's statement of its
    # resolution — floored at one grid step, below which nothing is representable at all.
    step = float(landscape.log_rho[1] - landscape.log_rho[0])
    consistent = bool(np.isfinite(gap) and gap <= max(depleted.width, step))

    # ── per-region enriched-basin responsibility, on the REGION axis
    w_slot = np.full(n_regions, np.nan, dtype=np.float64)
    if enriched is not None:
        grid10 = landscape.log_rho / np.log(10.0)  # _poisson_kernels takes a log10 grid
        kern = _poisson_kernels(c, e, grid10)
        post = kern * np.exp(landscape.logP)[None, :]
        post /= np.maximum(post.sum(axis=1, keepdims=True), _EPS)
        in_basin = (landscape.log_rho >= enriched.lo) & (landscape.log_rho <= enriched.hi)
        w_slot[sel] = post[:, in_basin].sum(axis=1)
    else:
        w_slot[sel] = 0.0

    # ── the training population, published: the kernel centres this fit was built from, in natural
    # log, with each region's coarse class. The centre expression is `fit_landscape`'s own
    # (`log10(max(count,1)) − log10(eff)`, then to nats) — the same floor, because a zero-count region
    # sits at its resolution wall rather than at −inf.
    train_log_rho = np.log(np.maximum(c, 1.0)) - np.log(e)
    exon = (sig[sel] & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    intron = (sig[sel] & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    # exon wins over intron, matching `signature.coarse_type_array`'s rule (imported semantics, not a
    # second table): 0 intergenic, 1 intron, 2 exon.
    train_class = np.where(exon, 2, np.where(intron, 1, 0)).astype(np.int64)

    return AbundanceLandscape(
        landscape=landscape,
        modes=modes,
        depleted=depleted,
        enriched=enriched,
        rho_0=float(np.exp(depleted.log_rho)),
        span_R=span_R,
        anchor_log_rho=anchor_log_rho,
        anchor_gap_nats=float(gap) if np.isfinite(gap) else float("nan"),
        anchor_consistent=consistent,
        w_slot=w_slot,
        n_train=int(landscape.n_train),
        train_log_rho=train_log_rho,
        train_class=train_class,
    )
