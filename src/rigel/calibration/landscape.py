"""The landscape: the population hyperprior on one component's density, ``P(log ρ_c)``.

It is fitted from a previous solve's deconvolved mass for one component and fed to that component's ψ arm on the
re-solve. The census (:func:`census`, :func:`split_basins`) splits the fitted density into basins; the capture reader
(``capture_efficiency``) takes the gDNA landscape's largest basin as the bulk and every basin above it as the
enriched spectrum. Layer 5.

The estimator is component-agnostic arithmetic on ``(count, mass, eff)``; which objects train it and what an anchor
means are the caller's. A component's density is typically a sharp mass plus a broad bump decades away (for gDNA, the
depleted level with capture lifting covered regions above it), so the estimator must resolve a spike and a wide bump
on one axis. Two rules follow: above the location floor (``CONSTANTS.landscape.located_var``) precision is a
continuous weight, never an admission threshold (:func:`_reliability`); and resolution is a population quantity, not
a measurement one (:func:`knn_widths`). Its constants are ``CONSTANTS.landscape``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from ..config import CONSTANTS
from ..native import transfer_rows as _rows
from .simplex_logodds import block_rows

#: the arm's clips, shared with the kernel arm that evaluates this prior (``psi_kernel.h``)
_EPS = float(_rows.ARM_EPS)
_LN10 = np.log(10.0)


@dataclass(frozen=True)
class DensityLandscape:
    """A fitted population gDNA-density hyperprior: ``logP`` over a natural-log rate grid, entering ψ as exact Bayes.

    ψ reads the curve at the density each cell of each slot implies, ``log ρ_c = log f_g + log M − log E`` on the
    slot's gDNA support, by linear interpolation with the ends held constant (``native/psi_kernel.h``, ``Arm``).
    ``logP`` is bare: a density in log-rate needs no Jacobian, since its conversion to linear rate cancels ψ's
    ``log σ'(λ)`` change of variable exactly, and ψ's arm adds the reference measure itself. The grid spans the
    data's support (:func:`log_rate_grid`), so the held ends act only at ψ's extreme fractions, where the prior has no
    more information."""

    log_rho: np.ndarray
    logP: np.ndarray
    n_train: int
    #: the population's own rendered density on the same grid, without the pseudo-region of ignorance
    #: (:func:`fit_landscape`); ``-inf`` where no kernel puts mass. ψ reads ``logP``; the capture reader reads this,
    #: because the pseudo-region's uniform mass is not a level any object was measured at.
    logP_free: np.ndarray

    def required_logodds_window(self, mass, eff) -> float:
        """The λ bracket this prior's own support demands, derived with no constant.

        ψ reads the curve at ``log ρ_i(f) = log f + log M_i − log E_i`` for ``f ∈ [σ(−L), σ(L)]``, so the prior is
        expressible at slot ``i`` only if ψ can place the slot on the landscape's floor ``ρ_floor = exp(log_rho[0])``::

            σ(−L)·M_i/E_i ≤ ρ_floor   ⇔   σ(−L) ≤ f_i := ρ_floor·E_i/M_i   ⇔   L ≥ log((1 − f_i)/f_i)

        The bracket is the maximum over the slots with mass and opportunity; for ``f_i ≪ 1`` it is the landscape's own
        log-dynamic range ``max_i log(M_i/E_i) − log ρ_floor``. A slot with ``f_i ≥ 1`` imposes nothing (``L_i ≤ 0``).

        The caller must take ``max(configured L, this)``: the configured bracket keeps the Beta(½,½) reference proper,
        and this is the fitted prior's additional demand. It is not an accuracy knob: it makes the answer independent
        of ``L``, so widening past it must change nothing.
        """
        m = np.asarray(mass, dtype=np.float64)
        e = np.asarray(eff, dtype=np.float64)
        live = (m > 0.0) & (e > 0.0)
        if not live.any():
            return 0.0
        rho_floor = float(np.exp(self.log_rho[0]))
        # f_i is a fraction, so clip into (0, 1): at f_i ≥ 1 the slot imposes nothing (L_i ≤ 0), and the lower clip
        # is a float guard, not a floor on the answer.
        f = np.clip(rho_floor * e[live] / m[live], _EPS, 1.0 - _EPS)
        return float(max(0.0, np.max(np.log1p(-f) - np.log(f))))


def log_rate_grid(mass: np.ndarray, eff: np.ndarray) -> np.ndarray:
    """The log10 rate axis: exactly the domain ψ's fitted arm can ask the curve about.

    ψ evaluates the prior at ``ρ_g = f_g·M/E`` for ``f_g ∈ (0, 1]``. The top is ``max_i log10(M_i/E_i)``: no region
    can sit above its own total density. The bottom is the resolution wall ``min_i log10(1/E_i)``: a region of
    effective length ``E`` with no gDNA says only ``ρ ≲ 1/E``, and since a zero-count kernel ``e^{−ρE}`` decreases
    monotonically, the floor is where the depleted anchor deposits its mass. Both bounds dominate every kernel centre
    ``log10(max(count, 1)/E)``, since ``count ≤ M``.

    The grid must not be padded by a kernel width: the pad would be set by the widest kernel, so one isolated region
    would stretch the span, coarsen the fixed-count grid and over-smooth the whole landscape.
    """
    lo = float(np.min(-np.log10(np.maximum(eff, _EPS))))
    hi = float(np.max(np.log10(np.maximum(mass, _EPS)) - np.log10(np.maximum(eff, _EPS))))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < _EPS:
        lo, hi = lo - 0.5, lo + 0.5
    return np.linspace(lo, hi, CONSTANTS.landscape.grid_points)


def poisson_kernels(count: np.ndarray, eff: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Per-region ``P(count | ρ·E)`` on the grid, each row normalised to unit mass; ``(n, K)``.

    At ``count = 0`` the kernel is ``e^{−ρE}``, which says "ρ is anywhere below the resolution wall"; a point-estimate
    KDE would invent a location at ``1/E`` instead.
    """
    lam = np.exp(grid * _LN10)[None, :] * eff[:, None]
    ll = count[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(count[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    return pn / np.maximum(pn.sum(1, keepdims=True), _EPS)


def knn_widths(centres: np.ndarray, grid_step: float) -> np.ndarray:
    """The population resolution: ``h_i = scale · dist(a_i, k-th nearest neighbour)``, ``k = √n`` and ``scale`` the
    ``CONSTANTS.landscape.knn_scale``.

    The per-region Poisson likelihood is a measurement width, ``1/(√g·ln10)`` decades on the log axis, which shrinks
    as ρ^(−1/2) and falls below one grid step where the enriched mode lives; measurement widths alone render a smooth
    bulk and a comb over the enriched band, reporting one broad mode as many. A population's resolution is how finely
    the sample resolves it, which is the nearest-neighbour spacing, and it self-corrects: fewer regions give farther
    neighbours and wider kernels. The ``grid_step`` floor is forced by the axis: a kernel narrower than one cell is a
    delta at the wrong height.

    The distance must be the true k-th-nearest-neighbour distance. The far edge of a 2k window,
    ``max(a_i − a_{i−k}, a_{i+k} − a_i)``, is systematically larger and hands the widest kernel to the most isolated
    region, smearing one observation across decades.
    """
    n = centres.size
    k = max(int(round(np.sqrt(n))), 2)
    srt = np.sort(centres)
    if n <= k:
        return np.full(n, max(float(srt[-1] - srt[0]), grid_step))
    pos = np.searchsorted(srt, centres)
    # The k+1 nearest points are a contiguous window [j, j+k] containing `pos`, so the distance is
    #     min over j in [pos−k, pos] of  max(a_pos − a_j, a_{j+k} − a_pos).
    # The first term falls and the second rises with j, so the minimum sits where `a_j + a_{j+k} − 2·a_pos` crosses
    # zero; bisection finds it exactly in O(n log k) rather than O(n·k).
    a = srt[np.minimum(pos, n - 1)]
    lo = np.maximum(pos - k, 0)
    hi = np.minimum(pos, n - 1 - k)
    j0, j1 = lo.copy(), np.maximum(hi, lo)
    while np.any(j0 < j1):
        mid = (j0 + j1) // 2
        rise = srt[mid] + srt[mid + k] - 2.0 * a >= 0.0
        j1 = np.where(rise, mid, j1)
        j0 = np.where(rise, j0, mid + 1)

    def _radius(j):
        j = np.clip(j, lo, np.maximum(hi, lo))
        return np.maximum(a - srt[j], srt[j + k] - a)

    d_k = np.minimum(_radius(j0), _radius(j0 - 1))
    return np.maximum(CONSTANTS.landscape.knn_scale * d_k, grid_step)


def _reliability(count: np.ndarray, var: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Per-region weight ``w = ref/(v + ref)``: the irreducible share of the log-rate variance against the
    deconvolution ambiguity ``v = Var(log f_g)``.

    ``ref`` is the region's Poisson counting floor ``1/max(count, 1)`` plus the reference variance
    (``CONSTANTS.landscape.reliability_sd_decades``). A confident region keeps its mass, so a real enriched mode
    survives; a region the deconvolution could not resolve (``v ≫ ref``) collapses toward zero. The zero-count
    structural anchor carries ``w = 1``: its density is 0 for every ``f_g``, so its composition ambiguity is
    irrelevant.

    Above the location floor (``CONSTANTS.landscape.located_var``) precision enters only as this continuous weight,
    never as an admission threshold or a kernel width. At the floor it is admission: a slot whose log-count is
    uncertain by more than one nat² has no location, and the caller does not train on it. On unstranded data the
    weight separates region classes (exons against introns and intergenic) rather than enriched from depleted within
    one, so it is informative but coarse, as coarse as the reference variance sets.
    """
    v = np.maximum(np.nan_to_num(var, nan=np.inf, posinf=np.inf), 0.0)
    reference_var = (CONSTANTS.landscape.reliability_sd_decades * _LN10) ** 2
    ref = 1.0 / np.maximum(count, 1.0) + reference_var
    return np.where(anchor, 1.0, ref / (v + ref))


def estep_kernels(kernels: np.ndarray, count: np.ndarray, grid: np.ndarray, prev) -> np.ndarray:
    """The E-step on the kernels that have no location.

    A region trained at less than one fragment is centred at its resolution wall, and its Poisson kernel is flat
    below it; normalised, it spreads its mass uniformly under the wall, which is the flat-prior posterior rather than a
    population statement, and short empty regions would deposit spurious mass at high density. With the previous fit
    at hand, such a kernel is placed where the population is: ``kernel × P_prev(ρ)``, renormalised. A counted kernel
    keeps its own location, so an enriched minority cannot be competed away by the bulk. ``prev`` is ``None`` at the
    first fit, where nothing changes. No constant enters: the one-fragment floor is the estimator's own.
    """
    if prev is None:
        return kernels
    p_prev = np.exp(
        np.interp(grid * _LN10, prev.log_rho, prev.logP, left=prev.logP[0], right=prev.logP[-1])
    )
    m = count < 1.0
    if not m.any():
        return kernels
    k = kernels[m] * p_prev[None, :]
    kernels[m] = k / np.maximum(k.sum(1, keepdims=True), _EPS)
    return kernels


def _render(count, eff, grid, prev, weights, widths) -> np.ndarray:
    """Sum the weighted kernels, each widened to the population resolution, into an unnormalised density.

    Convolution is linear, so the kernels are grouped by width (``CONSTANTS.landscape.width_bins`` groups) and each
    group's sum is convolved once. The kernels are built and summed a row tile at a time
    (``simplex_logodds.block_rows``), so the training set never exists as an ``(n, K)`` matrix. Each row joins its
    group's sum in row order whatever the tile height, so the tile height sets the working set and changes no bit.
    """
    n_bins = CONSTANTS.landscape.width_bins
    step = float(grid[1] - grid[0])
    boundaries = np.quantile(widths, np.linspace(0.0, 1.0, n_bins + 1))
    # the bin of each width: the last bin includes its upper boundary, every other bin excludes it
    bin_of = np.clip(np.searchsorted(boundaries, widths, side="right") - 1, 0, n_bins - 1)
    sums = np.zeros((n_bins, grid.size))
    rows = block_rows(grid.size, sums.itemsize)
    for r0 in range(0, count.size, rows):
        sl = slice(r0, r0 + rows)
        kernels = estep_kernels(poisson_kernels(count[sl], eff[sl], grid), count[sl], grid, prev)
        contrib = weights[sl, None] * kernels
        for b in np.unique(bin_of[sl]):
            part = contrib[bin_of[sl] == b]
            # the group's running sum enters as the part's first row and the sum runs down the rows in order, so no
            # tile boundary regroups an addition (adding a tile's own partial sum into the total would)
            part[0] += sums[b]
            sums[b] = part.sum(0)
    out = np.zeros_like(grid)
    for b in range(n_bins):
        m = bin_of == b
        if not m.any():
            continue
        d = sums[b]
        h = float(np.mean(widths[m]))
        if h > step:
            k = np.exp(-0.5 * ((grid[:, None] - grid[None, :]) / h) ** 2)
            d = (k / np.maximum(k.sum(0, keepdims=True), _EPS)) @ d
        out += d
    return out


def fit_landscape(
    count,
    mass,
    eff,
    var,
    *,
    anchor,
    domain: tuple | None = None,
    prev: "DensityLandscape | None" = None,
) -> "DensityLandscape | None":
    """Fit the landscape from a previous solve's per-region deconvolved gDNA; ``None`` if it cannot be fit.

    One entry per training region: ``count`` the deconvolved gDNA mass ``f_g·M``, ``mass`` the total unspliced mass
    ``M`` (which bounds the achievable density and so the grid top), ``eff`` the effective length, ``var`` the
    belief's ``Var(log f_g)``, and ``anchor`` the zero-mass structural regions (:func:`_reliability`). The caller
    selects the substrate, since that needs the chain.

    ``domain`` is the ``(mass, eff)`` of the population the prior will be read at, when it is wider than the training
    population (a slot with no composition does not train). The grid then spans the consumers while the kernels stay
    the training set's; without it, a consumer above the training set's top density would read a flat prior.

    ``prev`` is the previous refit's landscape, or ``None`` at the first fit (:func:`estep_kernels`).

    The estimator is a weighted sum of zero-native per-region kernels at the population resolution, with no EM, no
    competition between components and no iteration: it is deterministic, every region's contribution is traceable,
    and a capture-enriched minority cannot be competed away by the bulk.
    """
    count = np.asarray(count, dtype=np.float64)
    mass = np.asarray(mass, dtype=np.float64)
    eff = np.asarray(eff, dtype=np.float64)
    anchor = np.asarray(anchor, dtype=bool)
    live = np.isfinite(count) & np.isfinite(eff) & np.isfinite(mass) & (eff > _EPS)
    if int(live.sum()) < 2:
        return None
    count, mass, eff, anchor = np.maximum(count[live], 0.0), mass[live], eff[live], anchor[live]
    var = np.asarray(var, dtype=np.float64)[live]

    if domain is None:
        grid = log_rate_grid(mass, eff)
    else:
        d_mass = np.asarray(domain[0], dtype=np.float64)
        d_eff = np.asarray(domain[1], dtype=np.float64)
        d_live = np.isfinite(d_mass) & np.isfinite(d_eff) & (d_eff > _EPS)
        grid = (
            log_rate_grid(d_mass[d_live], d_eff[d_live])
            if d_live.any()
            else log_rate_grid(mass, eff)
        )
    centres = np.clip(np.log10(np.maximum(count, 1.0)) - np.log10(eff), grid[0], grid[-1])
    widths = knn_widths(centres, float(grid[1] - grid[0]))
    weights = _reliability(count, var, anchor)
    density = _render(count, eff, grid, prev, weights, widths)
    total = float(density.sum())
    if not (total > 0.0 and np.isfinite(total)):
        return None

    # one pseudo-region of complete ignorance, spread uniformly: Laplace smoothing where the "one" is one region, so
    # no constant. It bounds how hard the prior can push: a population of `W` weighted regions cannot resolve a cell
    # rarer than ~1/W, so the log-range is bounded by log W and the prior stays weak and correctable.
    free = density / total
    density = (density + total / (grid.size * max(float(weights.sum()), 1.0))) / total
    with np.errstate(divide="ignore"):
        log_free = np.log(free / free.sum())
    return DensityLandscape(
        log_rho=grid * _LN10,
        logP=np.log(density / density.sum()),
        n_train=int(live.sum()),
        logP_free=log_free,
    )


# ── The census: the bulk and the enriched basins ─────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class LandscapeMode:
    """One mode of a fitted density: the peak ``log_rho`` (natural log) and its basin ``[lo, hi]``, the
    minima flanking it, with ``basin_mass`` the share of the density inside. Adjacent basins share a bound,
    so the modes partition the grid."""

    log_rho: float
    basin_mass: float
    lo: float
    hi: float


def census(landscape: DensityLandscape) -> tuple[LandscapeMode, ...]:
    """Every local maximum of ``logP`` with its basin, in grid order. Basins split at the minimum between
    adjacent maxima, so they partition the grid and their masses sum to one.

    A maximum is an interior cell above its left neighbour and at least its right one (a plateau peaks at
    its left end), or an edge cell strictly above its one neighbour: a density that is highest at the edge
    of the grid peaks there, whatever else it holds. A flat edge is no peak. A density with no maximum at
    all peaks at its arg-max."""
    x = np.asarray(landscape.log_rho, dtype=np.float64)
    p = np.exp(np.asarray(landscape.logP, dtype=np.float64))
    p = p / max(float(p.sum()), _EPS)
    peaks = list(np.where((p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:]))[0] + 1)
    if p.size > 1 and p[0] > p[1]:
        peaks.insert(0, 0)
    if p.size > 1 and p[-1] > p[-2]:
        peaks.append(p.size - 1)
    if not peaks:
        peaks = [int(np.argmax(p))]
    cuts = [0]
    for a, b in zip(peaks[:-1], peaks[1:], strict=False):
        cuts.append(int(a + np.argmin(p[a : b + 1])))
    cuts.append(p.size - 1)
    modes = []
    for i, pk in enumerate(peaks):
        s, e = cuts[i], cuts[i + 1]
        # half-open segments, so a shared cut bin is counted once; the last basin takes the final point
        hi_idx = e + 1 if i == len(peaks) - 1 else e
        modes.append(
            LandscapeMode(
                log_rho=float(x[pk]),
                basin_mass=float(p[s:hi_idx].sum()),
                lo=float(x[s]),
                hi=float(x[e]),
            )
        )
    return tuple(modes)


def split_basins(
    modes: tuple[LandscapeMode, ...],
) -> tuple[LandscapeMode, tuple[LandscapeMode, ...]]:
    """The bulk basin, the largest by rendered mass, and the basins wholly above it, in grid order (empty
    when nothing lies above the bulk)."""
    bulk = max(modes, key=lambda m: m.basin_mass)
    above = tuple(m for m in modes if m.lo >= bulk.hi - _EPS and m is not bulk)
    return bulk, above
