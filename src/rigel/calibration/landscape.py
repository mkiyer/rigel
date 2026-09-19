"""The population component-density hyperprior — the landscape.

Fit `P(log ρ_c)` over the population from the previous solve's deconvolved mass for one component, then
feed it to that component's ψ composition arm on the re-solve.

The estimator is component-agnostic, deliberately so. Every step below is arithmetic on
``(count, mass, eff)`` for whichever component the caller selected: :func:`_grid` is ``mass/eff``,
:func:`_poisson_kernels` is `P(count | ρ·E)`, :func:`knn_widths` is nearest-neighbour spacing on a 1-D
axis, :func:`_render` is convolution. Nothing here knows which component it is fitting. The
component-specific reasoning — which objects train it, and what an anchor means — belongs to the caller
and lives beside the caller's substrate selector, where the chain is: the REGIONs-only / AMBIG-excluded /
anchor-bearing selection, and what a zero-count object is evidence of.

The shape of the truth it has to represent: a component's density is typically one near point mass plus
one broad bump some decades away — for gDNA, a uniform depleted level with hybrid capture lifting the
covered regions a couple of decades above it; for RNA, a silent majority against an expressed spread.
Two components, one sharp and one broad, several decades apart, so the estimator must resolve a spike
and a wide bump on the same axis. That requirement fixes the two rules the design rests on:

1. Above the location floor (:data:`_LOCATED_VAR`) precision is a continuous weight, never a tuned
   admission threshold: a tuned cutoff scores worse than ignoring precision. See :func:`_reliability`.
2. Resolution is a population quantity, not a measurement one. See :func:`knn_widths`, worth reading
   before touching the kernel.

⚠ The two modelling constants below (``_KNN_SCALE``, ``_S0``) were selected against gDNA-shaped data.
They are not component-specific by construction, but they have only ever been validated on one
component; a second caller inherits them and must say so in its results rather than discover it later.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from .simplex_logodds import _block_rows

_EPS = 1e-12
_LN10 = np.log(10.0)

# ── Computational budgets (discretization, not modelling — they trade cost for exactness) ────────────────
#: Points on the log-rate grid. Same role as the solver's λ lattice: finer is strictly more
#: faithful and strictly slower.
_N_GRID = 260
#: Kernels are grouped into this many equal-count width bins and each bin convolved once, instead of one
#: convolution per region. Pure speed: the cost goes from O(n·K²) to O(bins·K²), which is what makes this
#: affordable at genome scale.
_WIDTH_BINS = 12

# ── Modelling constants ──────────────────────────────────────────────────────────────────────────────────
#: Population-resolution scale for :func:`knn_widths`, selected by shape against a reference that is
#: itself validated against ground truth: at 0.5 the fit renders the enriched mode at the width the truth
#: has; below it the landscape combs, above it the two modes merge and the enriched mass collapses.
#: EMD does not discriminate here — it is monotone in smoothing at every reference — so do not
#: re-select on it.
_KNN_SCALE = 0.5
#: Reference variance scale in the reliability weight, as a log-rate variance (0.15 decades).
#: ⚠ It is a tuning constant in disguise. It reads as "the kernel resolution floor", but the actual
#: rendering resolution is the grid step (~0.025 dec), and substituting that makes the weight more
#: aggressive and the census worse. What it really does is cap how far a confident region can be
#: down-weighted. Changing it is its own measured experiment; it must not ride along with anything else.
_S0 = (0.15 * _LN10) ** 2

#: THE LOCATION FLOOR, in the variable every solve reports. The estimator's resolution wall is one fragment
#: (``max(count, 1)`` centres a kernel; ``count < 1`` is location-free and E-step placed), and a Poisson count
#: ``c`` has ``Var(log c) = 1/c``, so "below one fragment" is "the log-count is uncertain by more than one nat²".
#: `RegionBelief.var_gdna` is ``Var(log f_g)`` — at fixed mass, ``Var(log count)`` — whatever produced the
#: solve (the strand term, the factory, a delivered row), so a slot wider than this has no location by the
#: same floor the count rule applies, and does not train (`calibrate._fit_gdna_hyperprior`). Not a tuned
#: constant: the identity's value at the wall. Gate: `test_landscape_training_population.py`.
_LOCATED_VAR = 1.0


@dataclass(frozen=True)
class DensityLandscape:
    """A fitted population gDNA-density hyperprior: ``logP`` over a natural-log rate grid, entering ψ
    as exact Bayes (a temperature on the term was a tunable nothing ever moved; retired 2026-09-13).

    ψ reads the curve ITSELF, at the density every cell of every slot implies — ``log ρ_c = log f_g + log M −
    log E`` on the slot's gDNA support — with numpy's interpolation and the ends held constant off the grid
    (`native/psi_kernel.h`, ``Arm``; `sweep._gdna_arm` hands the kernel the curve and the per-slot support).
    Bare: no reference prior, no measure term, no Jacobian — ``logP`` is a density in log-rate, so its
    conversion to a linear-rate density cancels the ``log σ'(λ)`` change of variable exactly, per component,
    and ψ's arm adds the reference itself. The grid already spans the data's own support (:func:`_grid`), so
    the held ends fire only on ψ's extreme fractions, where the honest statement is "no more information out
    here", not a linear extension of the last slope."""

    log_rho: np.ndarray
    logP: np.ndarray
    n_train: int
    #: The kernels this density was rendered from, one per training region: the centre
    #: ``log(max(count, 1) / E)`` in nats, and whether it is a LOCATION — a count of at least one
    #: fragment, the same wall :data:`_LOCATED_VAR` is read at. A zero-count anchor or a sub-fragment
    #: kernel is centred at its resolution wall ``1/E``, which says where the kernel could not see, not
    #: where a density is; a consumer reading a mode off the density (`abundance_landscape.located_enriched_mode`)
    #: counts members among the located kernels only, and never re-derives either array.
    centre: np.ndarray
    located: np.ndarray

    def required_logodds_window(self, mass, eff) -> float:
        """The λ bracket this prior's own support demands — derived, with no constant chosen.

        ψ reads the curve at ``log ρ_i(f) = log f + log M_i − log E_i`` and can only offer
        ``f ∈ [σ(−L), σ(L)]``. So this prior is expressible at slot ``i`` only if ψ can place that slot on
        the landscape's own floor ``ρ_floor = exp(log_rho[0])`` — :func:`_grid`'s resolution wall::

            σ(−L)·M_i/E_i ≤ ρ_floor   ⇔   σ(−L) ≤ f_i := ρ_floor·E_i/M_i   ⇔   L ≥ log((1 − f_i)/f_i)

        and the bracket the chain needs is the maximum over the slots ψ evaluates. For the typical
        ``f_i ≪ 1`` that collapses to ``max_i log(M_i/E_i) − log ρ_floor``, this landscape's own
        log-dynamic range, which is why nothing is chosen here.

        A slot with ``f_i ≥ 1`` imposes no constraint: its own total density already sits at or below the
        wall, so ψ can express it at any bracket. Those return a non-positive ``L_i`` and are dominated by
        the maximum, which is why no mask is needed beyond "has mass and has opportunity".

        The caller must take ``max(configured L, this)``, never this alone. The configured bracket is
        the state-space range the Beta(½,½) reference needs to stay proper; this is the additional demand
        the fitted prior makes, and narrowing below the configured value trades one defect for another.

        This is not an accuracy knob and must not become one. It restores the property
        :mod:`.simplex_logodds` already claims as its acceptance test — that the answer does not depend on
        ``L`` — by making the bracket wide enough that the prior is representable. Widening past this must
        therefore change nothing, and that saturation is the gate.
        """
        m = np.asarray(mass, dtype=np.float64)
        e = np.asarray(eff, dtype=np.float64)
        live = (m > 0.0) & (e > 0.0)
        if not live.any():
            return 0.0
        rho_floor = float(np.exp(self.log_rho[0]))
        # f_i is a FRACTION, so clip into (0, 1): at f_i ≥ 1 the slot imposes nothing (L_i ≤ 0), and the
        # lower clip is the float guard, not a floor on the answer.
        f = np.clip(rho_floor * e[live] / m[live], _EPS, 1.0 - _EPS)
        return float(max(0.0, np.max(np.log1p(-f) - np.log(f))))


def _grid(mass: np.ndarray, eff: np.ndarray) -> np.ndarray:
    """The log10 rate axis: exactly the domain ψ's fitted arm can ask the curve about. No
    asserted range, and nothing left over to choose.

    ψ evaluates the prior at ``ρ_g = f_g·M/E`` for ``f_g ∈ (0, 1]``, so:

    * the top is a hard bound, ``max_i log10(M_i/E_i)``. No region can ever be placed above its own total
      density, so there is nothing above it to represent.
    * the bottom is the resolution wall, ``min_i log10(1/E_i)``. A region of effective length ``E`` that
      sequenced no gDNA says only "ρ ≲ 1/E"; below the deepest such wall nothing is distinguishable from
      anything else. This end matters more than it looks — a zero-count region's kernel ``e^{−ρE}`` is
      monotone decreasing, so where the floor sits is where the depleted anchor deposits its mass.

    Both bounds also dominate every kernel centre ``log10(max(count,1)/E)``, since ``count ≤ M``, so no
    centre is ever truncated.

    Do not pad this by a kernel width. The pad would be set by the single widest kernel, one isolated
    region would stretch the span by a couple of decades, and since the point count is fixed the grid step
    coarsens and the whole landscape over-smooths. A lone region's broad kernel is 1/n of the mass spread
    over decades; spending resolution to render its tails costs resolution everywhere that matters.
    """
    lo = float(np.min(-np.log10(np.maximum(eff, _EPS))))
    hi = float(np.max(np.log10(np.maximum(mass, _EPS)) - np.log10(np.maximum(eff, _EPS))))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < _EPS:
        lo, hi = lo - 0.5, lo + 0.5
    return np.linspace(lo, hi, _N_GRID)


def _poisson_kernels(count: np.ndarray, eff: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Per-region ``P(count | ρ·E)`` on the grid, row-normalised to unit mass → ``(n, K)``.

    Zero-native, and that is the point: at ``count = 0`` this is ``e^{−ρE}``, a monotone decay that says
    "ρ is anything below the resolution wall" — the honest statement, and the one a point-estimate KDE
    cannot make, since it would invent a location at ``1/E``.
    """
    lam = np.exp(grid * _LN10)[None, :] * eff[:, None]
    ll = count[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(count[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    return pn / np.maximum(pn.sum(1, keepdims=True), _EPS)


def knn_widths(
    centres: np.ndarray, grid_step: float, scale: float = _KNN_SCALE, k: int | None = None
) -> np.ndarray:
    """The population resolution: ``h_i = scale · dist(a_i, k-th nearest neighbour)``, ``k = √n``.

    ``k`` is the population's own by default; a caller asking about a SUBSET of the population at the
    population's resolution passes the population's ``k`` (`abundance_landscape.located_enriched_mode`
    reads a basin's members at the located population's ``√n``, so a cluster smaller than ``k`` is not
    narrow but unresolved).

    ⛔ Read this before changing the kernel. The per-region Poisson likelihood is a measurement width, and
    on the log axis it is ``1/(√g·ln10)`` decades, so it shrinks as ρ^(−1/2) — by well over an order of
    magnitude across a library, dropping below one grid step exactly where the enriched mode lives, so
    that the regions carrying that band's mass become deltas in a single cell. A landscape built from
    measurement widths alone is therefore smooth through the depleted bulk and a comb over the enriched
    one, and it reports one broad mode as dozens.

    A population's resolution is set by how finely the sample can resolve it, never by how precisely one
    member happens to be measured. Nearest-neighbour spacing is exactly that quantity — it is what decides
    whether two kernels merge into one mode or stand apart as two spurious ones — and it self-corrects
    with no tuning: fewer regions ⇒ farther neighbours ⇒ wider kernels.

    The ``grid_step`` floor is forced by the axis, not chosen: nothing narrower than one cell is
    representable, so a kernel below it is not a narrow density but a delta at the wrong height.

    Two alternatives are measured and rejected: a single global width, which flattens the genuinely sharp
    depleted bulk in order to serve the sparse tail, one width being unable to serve regions whose
    required resolution differs by an order of magnitude; and Abramson sample-point adaptivity
    ``h ∝ f(a)^(−1/2)``, which is about neutral, because a spike is high pilot density at its own location
    and so keeps itself narrow, reinforcing exactly what needs merging.

    ⛔ This must be the true k-th-nearest-neighbour distance, and getting it wrong is not cosmetic. Taking
    ``max(a_i − a_{i−k}, a_{i+k} − a_i)`` instead — the far boundary of a 2k window — is systematically
    larger and, for a region with no near neighbours on one side, reaches all the way back into the bulk.
    That hands the widest kernel in the fit to the most isolated region, which is precisely backwards: an
    isolated observation is one observation, and smearing it across decades asserts population mass
    everywhere it touches. On a zero-gDNA library, where every such region is a false positive, the
    handful of over-called regions is then rendered as a permissive plateau across the whole enriched
    range.
    """
    n = centres.size
    k = max(int(round(np.sqrt(n))), 2) if k is None else int(k)
    srt = np.sort(centres)
    if n <= k:
        return np.full(n, max(float(srt[-1] - srt[0]), grid_step))
    pos = np.searchsorted(srt, centres)
    # The k+1 nearest points are a contiguous window [j, j+k] containing `pos`, so the distance is
    #     min over j in [pos−k, pos] of  max(a_pos − a_j, a_{j+k} − a_pos).
    # The first term falls and the second rises with j, so the objective is V-shaped and the minimum sits at
    # the crossing of `a_j + a_{j+k} − 2·a_pos`. Bisect for it — exact, and O(n log k) rather than O(n·k).
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
    return np.maximum(scale * d_k, grid_step)


def _reliability(count: np.ndarray, var: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Per-region mass ``w = ref/(v + ref)`` — the irreducible share of the log-rate variance against the
    deconvolution ambiguity ``v = Var(log f_g)``.

    ``ref`` sums the region's own Poisson counting floor ``1/max(count,1)`` and the reference scale
    :data:`_S0`. A confident region keeps mass so a real enriched mode survives; a give-up region
    (``v ≫ ref``) collapses toward zero. The zero-count structural anchor is the trusted "no gDNA here"
    statement and carries ``w = 1``: its density is ``0`` for every ``f_g``, so its composition ambiguity
    is irrelevant.

    Precision enters twice, and the two roles are distinct. ABOVE the location floor (:data:`_LOCATED_VAR`)
    it is this continuous weight and nothing else: a tuned admission threshold on it was measured worse than
    ignoring precision, and moving it into the kernel width is refuted by its own control (a single constant
    width performs identically, so that form is a global bandwidth under another name, and it inflates false
    enrichment on zero-gDNA libraries). AT the floor it is admission — a slot whose log-count is uncertain by
    more than one nat² has no location, as a count below one fragment has none, and it is not a training
    slot (2026-09-14: the four ladder zero controls 500 → 282, 211 → 194, 550 → 265, 231 → 172 with every
    contaminated stratum unchanged or better).

    On unstranded data the weight does not separate enriched from depleted within a region class; it
    separates classes (exons being down-weighted against introns and intergenic), and every enriched
    region is an exon. It is therefore informative but coarse, and how coarse is set by :data:`_S0`.
    """
    v = np.maximum(np.nan_to_num(var, nan=np.inf, posinf=np.inf), 0.0)
    ref = 1.0 / np.maximum(count, 1.0) + _S0
    return np.where(anchor, 1.0, ref / (v + ref))


def _estep_kernels(kernels: np.ndarray, count: np.ndarray, grid: np.ndarray, prev) -> np.ndarray:
    """The E-step on the kernels that have no location.

    A region trained at less than one fragment has no location of its own — the estimator already centres
    it at its resolution wall, ``max(count, 1)`` — and its Poisson kernel is flat below that wall.
    Normalised to unit mass, such a kernel spreads that mass uniformly under the wall, which is the
    flat-prior posterior and not a population statement: short empty regions, whose wall sits at exon
    densities, then deposit a tail of landscape mass at high density on a library that has none. That
    tail is the reason a blind exon is not given a posterior median.

    The refit loop already holds the previous fit, and the deconvolution's E-step places a location-free
    kernel where the population is: ``kernel × P_prev(ρ)``, renormalised. A counted kernel keeps its own
    location and is not touched, so an enriched minority cannot be competed away by the bulk — which is
    what happens when the E-step is applied to every kernel. ``prev`` is ``None`` at the first fit, where
    nothing changes. What acts is where the previous fit puts the mass: a mirrored previous landscape, as
    a control, makes the capture-ON zero rows far worse. No constant enters, the one-fragment floor being
    the estimator's own."""
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
    """Sum the weighted kernels, each widened to the population resolution → an unnormalised density.

    Convolution is linear, so widening every kernel and then summing equals summing and then convolving
    — which is why kernels are grouped by width and each group convolved once (:data:`_WIDTH_BINS`). The
    kernels are built and summed a row tile at a time (`simplex_logodds._block_rows`, the same
    working-set rule ψ tiles by): a training set of a million regions never exists as an ``(n, K)``
    matrix, only its ``_WIDTH_BINS`` weighted sums do.
    """
    step = float(grid[1] - grid[0])
    boundaries = np.quantile(widths, np.linspace(0.0, 1.0, _WIDTH_BINS + 1))
    # the bin of each width: the last bin includes its upper boundary, every other bin excludes it
    bin_of = np.clip(np.searchsorted(boundaries, widths, side="right") - 1, 0, _WIDTH_BINS - 1)
    sums = np.zeros((_WIDTH_BINS, grid.size))
    rows = _block_rows(grid.size, 8)
    for r0 in range(0, count.size, rows):
        sl = slice(r0, r0 + rows)
        kernels = _estep_kernels(_poisson_kernels(count[sl], eff[sl], grid), count[sl], grid, prev)
        for b in np.unique(bin_of[sl]):
            m = bin_of[sl] == b
            sums[b] += (weights[sl][m][:, None] * kernels[m]).sum(0)
    out = np.zeros_like(grid)
    for b in range(_WIDTH_BINS):
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
    knn_scale: float = _KNN_SCALE,
    domain: tuple | None = None,
    prev: "DensityLandscape | None" = None,
) -> "DensityLandscape | None":
    """Fit the landscape from pass-0's per-region deconvolved gDNA. Returns ``None`` if it cannot be fit.

    Parameters mirror one training region each: ``count`` the deconvolved gDNA mass ``f_g·M``, ``mass`` the
    region's total unspliced mass ``M`` (which bounds the achievable density and so fixes the grid top),
    ``eff`` the effective length, ``var`` the belief's ``Var(log f_g)``, and ``anchor`` the zero-mass
    structural regions (see :func:`_reliability`). Substrate selection is the caller's job, since it needs
    the chain.

    ``domain`` is the ``(mass, eff)`` of the population the prior will be read at, when that is wider than
    the population it is fitted on. :func:`_grid` spans the data's own support and :meth:`logprior` clamps
    flat beyond it, which is right only while the two populations coincide; where the training population
    is a subset — a slot with no composition does not train — a consumer above the training set's top
    density would read a flat prior and fall back to ψ's reference. On a gDNA-free library whose exons are
    all blind, that collapses the grid to a decade at the floor and invents gDNA. With ``domain`` the grid
    is the consumers' and the kernels are the training set's.

    ``prev`` is the previous refit's landscape, or ``None`` at the first fit: the E-step on the
    location-free kernels (:func:`_estep_kernels`).

    The estimator is a weighted sum of zero-native per-region kernels at the population resolution. There
    is no EM, no competition between components and no iteration, so it is deterministic and every
    region's contribution is traceable — and a capture-enriched minority cannot be competed away by the
    bulk, which is how a mixture fitted by EM loses one.
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
        grid = _grid(mass, eff)
    else:
        d_mass = np.asarray(domain[0], dtype=np.float64)
        d_eff = np.asarray(domain[1], dtype=np.float64)
        d_live = np.isfinite(d_mass) & np.isfinite(d_eff) & (d_eff > _EPS)
        grid = _grid(d_mass[d_live], d_eff[d_live]) if d_live.any() else _grid(mass, eff)
    centres = np.clip(np.log10(np.maximum(count, 1.0)) - np.log10(eff), grid[0], grid[-1])
    widths = knn_widths(centres, float(grid[1] - grid[0]), knn_scale)
    weights = _reliability(count, var, anchor)
    density = _render(count, eff, grid, prev, weights, widths)
    total = float(density.sum())
    if not (total > 0.0 and np.isfinite(total)):
        return None

    # One pseudo-region of complete ignorance, spread uniformly. Ordinary Laplace smoothing, and the "one"
    # is one region — no constant. It bounds how hard this prior can push: without it the empty cells take
    # whatever floor the arithmetic happens to underflow to, which is an assertion the sample cannot
    # support. A population of `W` weighted regions cannot resolve a cell rarer than ~1/W, so the log-range
    # is bounded by log W and the prior stays weak and correctable — the governing principle for pass-0
    # output, which is exactly what this is fitted from.
    density = (density + total / (grid.size * max(float(weights.sum()), 1.0))) / total
    return DensityLandscape(
        log_rho=grid * _LN10,
        logP=np.log(density / density.sum()),
        n_train=int(live.sum()),
        centre=centres * _LN10,
        located=count >= 1.0,
    )
