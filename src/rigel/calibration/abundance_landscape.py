"""rigel.calibration.abundance_landscape — the pre-pass-0 TOTAL-density field, and its mode census.

⭐⭐ **The one question this module answers: what does the library's TOTAL fragment density look like
over the genome, BEFORE anything is solved — and which regions sit on which mode?** Under hybrid
capture the field is bimodal by construction (a depleted off-target level and an enriched on-target
one, ~2–3 decades apart); off capture it is unimodal. The census reads that structure off a fitted
density: `rho_0` (the depleted mode — the level the pooled intergenic anchors also measure,
independently), the span `R` (the mode ratio — never the 2.6–3.6×-under-reading in-gene anchors), and
per region a responsibility `w_i` for the enriched basin. Those three are what the pass-0 measured
reference consumes.

⭐ **The estimator is `landscape.fit_landscape`, reused verbatim.** It is deliberately
component-agnostic, and every decision in it transfers: zero-native Poisson kernels (a wall-exact
region that sequenced nothing says "below the resolution wall", not "at 1/E"), knn population
resolution (what suppresses combing with no tuning), the one-pseudo-region Laplace floor, and the grid
derived from the data's own support. What is NEW here is only the inputs and the census:

* the inputs are the MEASURED TOTALS — `total_abundance.region_counts_and_exposure`'s side-selected
  START/END counts over the region's own length, REGIONs only (boundaries cross rather than contain
  and have zero genomic measure — the same geometry ruling the gDNA hyperprior's substrate carries),
  restricted to the wall-exact (`model_free`) population;
* `var = 0` — ⭐ a DIRECT measurement has no deconvolution ambiguity, so the reliability weights are
  honestly flat. The hyperprior's `Var(log f_g)` weight exists because its training data came out of a
  solve; this data did not.

⛔ **NO SIGNIFICANCE THRESHOLD EXISTS IN THE CENSUS, ANYWHERE.** Every interior local maximum is a
mode; the grid is partitioned into basins at the minima between them; masses carry every verdict
continuously. A phantom wiggle above the bulk owns ~zero basin mass, so it yields `w ≈ 0` — harmless —
and the capture-OFF unimodality gate MEASURES that instead of a constant asserting it. The depleted
mode is picked by an independent measurement (the basin containing the pooled intergenic anchor rate),
and the anchor-consistency verdict's tolerance is the depleted mode's own fitted width — the density's
statement of its resolution, never a chosen number.

⚠ **T CONFLATES ENRICHMENT WITH EXPRESSION, AND THAT IS RECORDED UP FRONT**: a hot unprobed exon and a
probed cold one can land in the same basin, so the failure direction of any consumer is PERMISSIVE (an
over-wide enriched basin), never a hard exclusion. Pricing that confound is the plan's rung-5 read.

⚠ `_KNN_SCALE` / `_S0` are inherited from `landscape` and have only ever been validated on gDNA-shaped
data (that module's own warning) — every result reported off this landscape carries that caveat until
it is priced.

⛔ **Nothing here decides anything in the solve.** Its consumers are the census, the QC/injection
surface, and (behind its own flag, separately gated) the pass-0 measured reference. Its own
falsification is ``tests/calibration/test_abundance_landscape.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .landscape import DensityLandscape, _poisson_kernels, fit_landscape
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS
from .total_abundance import RegionWallMask, region_counts_and_exposure

_EPS = 1e-12
_GENE_BITS = BIT_INTRON_POS | BIT_INTRON_NEG | BIT_EXON_POS | BIT_EXON_NEG


__all__ = [
    "AbundanceLandscape",
    "AbundanceMode",
    "fit_abundance_landscape",
]


@dataclass(frozen=True, slots=True)
class AbundanceMode:
    """One mode of the fitted total-density field: a local maximum and its basin.

    ``log_rho`` is the peak's location (natural-log rate); ``lo``/``hi`` the basin bounds (the
    interior minima flanking it — adjacent basins share a bound, so the modes PARTITION the grid);
    ``basin_mass`` the share of the fitted density inside the basin; ``width`` the mass-weighted
    standard deviation of ``log_rho`` WITHIN the basin — the fit's own statement of how precisely this
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

    ``anchor_log_rho`` is the INDEPENDENT depleted-level estimator — the pooled rate over intergenic
    model-free regions (the same composition-free pool ``fit_intron_background`` uses) — and
    ``anchor_consistent`` says whether it falls within the depleted mode's own width. ``NaN`` (and a
    ``False``-free fallback to the largest basin) when no intergenic region exists, which is a toy.
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
        # half-open segments so the shared cut bin is counted ONCE — the last basin takes the final
        # grid point. Without this the basin masses sum to 1 + (mass at each cut), which the
        # partition gate caught at 1.0000047.
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
    # ⚠ shared bounds: mode i's hi IS mode i+1's lo (they were cut at one index). The construction
    # guarantees it; the partition gate asserts it.
    return tuple(modes)


def fit_abundance_landscape(
    substrate, region_arrays, wall_mask: RegionWallMask
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
        strength=1.0,
    )
    if landscape is None:
        return None

    modes = _census(landscape)

    # ── the INDEPENDENT depleted-level estimator: the pooled intergenic rate, composition-free
    sig = np.asarray(region_arrays.signature, dtype=np.int64)
    anchors = sel & ((sig & _GENE_BITS) == 0)
    if anchors.any() and float(exposure[anchors].sum()) > 0.0 and counts[anchors].sum() > 0.0:
        anchor_log_rho = float(np.log(counts[anchors].sum()) - np.log(exposure[anchors].sum()))
    else:
        anchor_log_rho = float("nan")

    # ── depleted = the basin CONTAINING the anchor rate; fallback (a toy): the largest basin
    if np.isfinite(anchor_log_rho):
        inside = [m for m in modes if m.lo <= anchor_log_rho <= m.hi]
        depleted = (
            inside[0] if inside else min(modes, key=lambda m: abs(m.log_rho - anchor_log_rho))
        )
    else:
        depleted = max(modes, key=lambda m: m.basin_mass)

    # ── enriched = the largest-mass basin strictly ABOVE the depleted one; None ⇒ unimodal
    above = [m for m in modes if m.lo >= depleted.hi - _EPS and m is not depleted]
    enriched = max(above, key=lambda m: m.basin_mass) if above else None

    span_R = float(np.exp(enriched.log_rho - depleted.log_rho)) if enriched is not None else 1.0
    gap = abs(depleted.log_rho - anchor_log_rho) if np.isfinite(anchor_log_rho) else float("nan")
    # the tolerance is the depleted mode's OWN fitted width — the density's statement of its
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
    )
