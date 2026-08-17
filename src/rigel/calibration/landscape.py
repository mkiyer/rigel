"""The population component-density hyperprior — the **landscape**.

Fit `P(log ρ_c)` over the population from the previous solve's deconvolved mass for ONE component, then
feed it to that component's ψ composition arm on the re-solve. It replaces the δ-pin `DensityNPMLE` in
that role (which is retired, not deleted — see :mod:`.npmle`).

⭐⭐ **THE ESTIMATOR IS COMPONENT-AGNOSTIC, AND DELIBERATELY SO.** Every step below is arithmetic on
``(count, mass, eff)`` for whichever component the caller selected: :func:`_grid` is ``mass/eff``,
:func:`_poisson_kernels` is `P(count | ρ·E)`, :func:`knn_widths` is nearest-neighbour spacing on a
1-D axis, :func:`_render` is convolution. **Nothing here knows which component it is fitting.**
⛔ **The component-specific reasoning — which objects train it, and what an anchor MEANS — is the
caller's**, and it lives beside each caller's substrate selector where the chain is. Two decisions that
used to be stated here are stated there instead: the REGIONs-only / AMBIG-excluded / anchor-bearing
selection, and what a zero-count object is evidence *of*.

**The shape of the truth it has to represent.** A component's density is typically one near point mass
plus one broad bump some decades away — for gDNA, a uniform depleted level with hybrid capture lifting
the covered regions ~2.7 decades above it; for RNA, a silent majority against an expressed spread. Two
components, one sharp and one broad, several decades apart — so the estimator must resolve a spike *and*
a wide bump on the same axis, which is what fixes the design below.

**Two decisions, each measured**:

1. **PRECISION IS A CONTINUOUS WEIGHT, NEVER AN ADMISSION RULE.** A hard precision cutoff scores *worse
   than ignoring precision altogether* (+0.175 vs +0.001 ambig; +0.359 vs +0.098 quick). See
   :func:`_reliability`.

2. **RESOLUTION IS A POPULATION QUANTITY, NOT A MEASUREMENT ONE.** See :func:`knn_widths` — this is the
   one that had a real bug in it, and it is worth reading before touching the kernel.

⚠ **The two constants below were SELECTED AGAINST gDNA-SHAPED DATA** (`_KNN_SCALE`, `_S0`). They are
not component-specific by construction, but they have only ever been validated on one component; a
second caller inherits them and must say so in its results rather than discover it later.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

_EPS = 1e-12
_LN10 = np.log(10.0)

# ── Computational budgets (discretization, not modelling — they trade cost for exactness) ────────────────
#: Points on the log-rate grid. Same role as the solver's ``n_grid`` / ``n_tilt``: finer is strictly more
#: faithful and strictly slower. 260 is what every measurement in the production plan was taken at.
_N_GRID = 260
#: Kernels are grouped into this many equal-count width bins and each bin convolved once, instead of one
#: convolution per region. Pure speed: the cost goes from O(n·K²) to O(bins·K²), which is what makes this
#: affordable at genome scale (1.5 M regions).
_WIDTH_BINS = 12

# ── Modelling constants (both tracked in the production plan's constants ledger) ─────────────────────────
#: Population-resolution scale for :func:`knn_widths`. **Selected by SHAPE against a reference that is
#: itself validated against ground truth** (N5): at 0.5 both the fit and the reference render the enriched
#: mode at the width the truth actually has (0.356 / 0.328 vs a true 0.307); below it the landscape combs,
#: above it the two modes start merging and the enriched mass collapses. ⚠ EMD does NOT discriminate here —
#: it is monotone in smoothing at every reference — so do not re-select on it.
_KNN_SCALE = 0.5
#: Reference variance scale in the reliability weight, as a log-rate variance (0.15 decades).
#: ⚠ It is a tuning constant in disguise: it is presented as "the kernel resolution floor", but the actual
#: rendering resolution is the grid step (~0.025 dec), and substituting that makes the weight *more*
#: aggressive and the census *worse*. What it really does is cap how far a confident region can be
#: down-weighted. Changing it is its own measured experiment; it must not ride along with anything else.
#:
#: ⛔ **IT IS NO LONGER THE BINDING CONSTRAINT ON THE ENRICHED CENSUS, AND THIS NOTE USED TO SAY IT WAS.**
#: The old reading — enriched recovery 0.69 at capture-ON and 0.29 at VSTRONG against 0.84 / 1.24 with flat
#: weights — was taken on the retired 32-condition suite and predates the 2026-08-06 pass-0 change, which
#: shrank ``Var(log f_g)`` across the board and so drove ``ref/(v+ref)`` toward 1. Re-measured on the
#: 36-condition ladder through the production refit: the mean weight is **0.90–0.93** (not the 0.46 this
#: file records for exons) and shipped-vs-flat moves the enriched prior mass by **0.97× / 1.07× / 1.08×**
#: at ``g75`` capture-ON / capture-OFF / ``g00`` — i.e. nothing. ⚠ VSTRONG is not a ladder condition, so the
#: old worst case is *unreproduced here*, not disproven.
_S0 = (0.15 * _LN10) ** 2


@dataclass(frozen=True)
class DensityLandscape:
    """A fitted population gDNA-density hyperprior: ``logP`` over a natural-log rate grid.

    ``strength`` is a temperature on the whole term. Default 1 is exact Bayes; below 1 tempers a prior that
    was, after all, fitted from *biased* pass-0 output, which is robustness rather than a fudge — it is what
    lets real data overcome a wrong prior. It lives here rather than at the call site so the object is
    self-describing and the sweep needs no knowledge of it.
    """

    log_rho: np.ndarray
    logP: np.ndarray
    n_train: int
    strength: float = 1.0

    def logprior(self, frac_grid, mass, eff) -> np.ndarray:
        """Project onto the ψ solve grid → ``(n_slots, K)`` additive term ``= log P(log ρ_c)`` evaluated at
        ``ρ_c = f_c·M/E_c``. **Bare** — no reference prior, no measure term, no Jacobian; ψ's arm adds the
        reference itself.

        ⭐⭐ **``frac_grid`` is the COMPONENT's fraction, which is what makes this serve either component.**
        The caller passes ``f_g`` for gDNA and ``1 − f_g`` for RNA, with that component's own opportunity
        as ``eff``. Nothing else differs between the two: the arithmetic below is
        ``log ρ_c = log f_c + log M − log E_c`` for whichever component is asked about.

        ⚠ **``(n_slots, K)``, not ``(n_regions, K)``** — the first axis is the unified region+boundary
        CHAIN, the same axis as ``u_pos``/``u_neg``. This docstring said "regions" until 2026-08-15 and a
        reader sizing a new array off it would have built the wrong shape.

        The latents are rates, conditioning on ``M`` contributes no ``f_c``-dependent Jacobian, and because
        ``logP`` is a density in **log**-rate its conversion to a linear-rate density cancels the
        ``log σ'(λ)`` change-of-variable exactly — so neither is written, here or in the caller. ⭐ That
        cancellation is *per component*, which is why it keeps holding when both arms carry a fitted prior.

        Off-grid values clamp to the end values rather than extrapolating: the grid already spans the data's
        own support (:func:`_grid`), so this only fires on the ψ grid's extreme fractions, where the honest
        statement is "no more information out here", not a linear extension of the last slope.
        """
        eff = np.maximum(np.asarray(eff, dtype=np.float64), _EPS)
        mass = np.maximum(np.asarray(mass, dtype=np.float64), _EPS)
        frac = np.clip(np.asarray(frac_grid, dtype=np.float64), _EPS, 1.0 - _EPS)
        log_rho_c = np.log(frac)[None, :] + (np.log(mass) - np.log(eff))[:, None]
        lp = np.interp(
            log_rho_c.ravel(), self.log_rho, self.logP, left=self.logP[0], right=self.logP[-1]
        ).reshape(log_rho_c.shape)
        return lp if self.strength == 1.0 else float(self.strength) * lp


def _grid(mass: np.ndarray, eff: np.ndarray) -> np.ndarray:
    """The log10 rate axis = **exactly the domain :meth:`DensityLandscape.logprior` can be asked about.** No
    asserted range, and nothing left over to choose.

    ψ evaluates the prior at ``ρ_g = f_g·M/E`` for ``f_g ∈ (0, 1]``, so:

    * **the top is a hard bound**, ``max_i log10(M_i/E_i)``. No region can ever be placed above its own total
      density, so there is nothing above it to represent.
    * **the bottom is the resolution wall**, ``min_i log10(1/E_i)``. A region of effective length ``E`` that
      sequenced no gDNA says only "ρ ≲ 1/E"; below the deepest such wall nothing is distinguishable from
      anything else. This end matters more than it looks — a zero-count region's kernel ``e^{−ρE}`` is
      monotone decreasing, so where the floor sits is where the depleted anchor deposits its mass.

    Both bounds also dominate every kernel centre ``log10(max(count,1)/E)``, since ``count ≤ M``, so no
    centre is ever truncated.

    ⚠ **Do not pad this by a kernel width.** That was tried: the pad is set by the single widest kernel, one
    isolated region then stretches the span by 2–3 decades, and since the point count is fixed the grid step
    coarsens and the whole landscape over-smooths (measured +0.007 / +0.056 EMD, enriched width 0.42 → 0.54).
    A lone region's broad kernel is 1/n of the mass spread over decades; spending resolution to render its
    tails costs resolution everywhere that matters.
    """
    lo = float(np.min(-np.log10(np.maximum(eff, _EPS))))
    hi = float(np.max(np.log10(np.maximum(mass, _EPS)) - np.log10(np.maximum(eff, _EPS))))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < _EPS:
        lo, hi = lo - 0.5, lo + 0.5
    return np.linspace(lo, hi, _N_GRID)


def _poisson_kernels(count: np.ndarray, eff: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Per-region ``P(count | ρ·E)`` on the grid, row-normalised to unit mass → ``(n, K)``.

    **Zero-native, and that is the point**: at ``count = 0`` this is ``e^{−ρE}``, a monotone decay that says
    "ρ is anything below the resolution wall" — the honest statement, and the one a point-estimate KDE
    cannot make (it would invent a location at ``1/E``).
    """
    lam = np.exp(grid * _LN10)[None, :] * eff[:, None]
    ll = count[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(count[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    return pn / np.maximum(pn.sum(1, keepdims=True), _EPS)


def knn_widths(centres: np.ndarray, grid_step: float, scale: float = _KNN_SCALE) -> np.ndarray:
    """THE population resolution: ``h_i = scale · dist(a_i, k-th nearest neighbour)``, ``k = √n``.

    ⚠⚠ **READ THIS BEFORE CHANGING THE KERNEL.** The per-region Poisson likelihood is a *measurement* width,
    and on the log axis it is ``1/(√g·ln10)`` decades — so it **shrinks as ρ^(−1/2)**. Across the benchmark
    that is a **41× collapse** (0.317 → 0.0078 dec) which drops **below one grid step** exactly where the
    enriched mode lives: above +1 decade, 88 % of regions carrying 100 % of that band's mass become deltas in
    a single cell. A landscape built from measurement widths alone is therefore smooth through the depleted
    bulk and a **comb** over the enriched one (roughness 46.9 vs the 2–4 of a smooth bump), and it reports
    that one broad mode as dozens.

    A population's resolution is set by how finely the *sample* can resolve it, never by how precisely one
    member happens to be measured. Nearest-neighbour spacing is exactly that quantity — it is what decides
    whether two kernels merge into one mode or stand apart as two spurious ones — and it self-corrects with
    no tuning: fewer regions ⇒ farther neighbours ⇒ wider kernels (mode-count-vs-sample-size slope +5.6 → −1.0).

    The `grid_step` floor is **forced by the axis, not chosen**: nothing narrower than one cell is
    representable, so a kernel below it is not a narrow density but a delta at the wrong height.

    ⚠ Two alternatives, both measured and both rejected: a single GLOBAL width (flattens the genuinely-sharp
    depleted bulk 7–9× to fix the sparse tail — one width cannot serve regions whose required resolution
    differs 10×), and ABRAMSON sample-point adaptivity ``h ∝ f(a)^(−1/2)`` (~neutral, because a spike *is*
    high pilot density at its own location and so keeps itself narrow — self-reinforcing exactly what needs
    merging).

    ⚠⚠ **This must be the true k-th-nearest-neighbour distance, and getting it wrong is not cosmetic.** The
    first implementation used ``max(a_i − a_{i−k}, a_{i+k} − a_i)`` — the FAR boundary of a 2k window, which is
    systematically larger and, for a region with no near neighbours on one side, reaches all the way back into
    the bulk. That hands **the widest kernel in the fit to the most isolated region**, which is precisely
    backwards: an isolated observation is one observation, and smearing it across decades asserts population
    mass everywhere it touches. Measured on zero-gDNA libraries, where every such region is a false positive:
    p99 width **1.04 decades → 0.21–0.46**, and the handful of over-called regions (6–126 of ~1100, carrying
    0.03–1.5 % of the weight) were being rendered as a permissive plateau across the whole enriched range.
    """
    n = centres.size
    k = max(int(round(np.sqrt(n))), 2)
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
    deconvolution AMBIGUITY ``v = Var(log f_g)``.

    ``ref`` sums the region's own Poisson counting floor ``1/max(count,1)`` and the reference scale
    :data:`_S0`. A confident region keeps mass so a real enriched mode survives; a give-up region (``v ≫ ref``)
    collapses toward zero. The zero-count structural anchor is the trusted "no gDNA here" statement and
    carries ``w = 1``: its density is ``0`` for *every* ``f_g``, so its composition ambiguity is irrelevant.

    ⚠ **Precision belongs here — as a continuous weight — and NOWHERE ELSE.** Expressing it as an admission
    threshold is measurably worse than ignoring precision entirely. And moving it into the kernel *width*
    instead was derived, implemented and **refuted by its own control**: a single constant width performs
    identically, so it is a global bandwidth under another name, and it inflates false enrichment on
    zero-gDNA libraries from 6.2× to 35×.

    ⚠ On unstranded data this weight does not separate enriched from depleted *within* a region class (exon
    0.466 vs 0.459) — it separates *classes* (exon 0.46 / intron 0.95 / intergenic 1.00), and every enriched
    region is an exon. So it is an *informative* weight, and the standing deficit in :data:`_S0` is its size.
    """
    v = np.maximum(np.nan_to_num(var, nan=np.inf, posinf=np.inf), 0.0)
    ref = 1.0 / np.maximum(count, 1.0) + _S0
    return np.where(anchor, 1.0, ref / (v + ref))


def _render(
    kernels: np.ndarray, weights: np.ndarray, widths: np.ndarray, grid: np.ndarray
) -> np.ndarray:
    """Sum the weighted kernels, each widened to the population resolution → an unnormalised density.

    Convolution is linear, so widening every kernel and then summing equals summing and then convolving —
    which is why kernels can be grouped by width and each group convolved once (:data:`_WIDTH_BINS`).
    """
    step = float(grid[1] - grid[0])
    boundaries = np.quantile(widths, np.linspace(0.0, 1.0, _WIDTH_BINS + 1))
    out = np.zeros_like(grid)
    for b in range(_WIDTH_BINS):
        upper = widths <= boundaries[b + 1] if b == _WIDTH_BINS - 1 else widths < boundaries[b + 1]
        m = (widths >= boundaries[b]) & upper
        if not m.any():
            continue
        d = (weights[m][:, None] * kernels[m]).sum(0)
        h = float(np.mean(widths[m]))
        if h > step:
            k = np.exp(-0.5 * ((grid[:, None] - grid[None, :]) / h) ** 2)
            d = (k / np.maximum(k.sum(0, keepdims=True), _EPS)) @ d
        out += d
    return out


def fit_landscape(
    count, mass, eff, var, *, anchor, strength: float = 1.0, knn_scale: float = _KNN_SCALE
) -> "DensityLandscape | None":
    """Fit the landscape from pass-0's per-region deconvolved gDNA. Returns ``None`` if it cannot be fit.

    Parameters mirror one training region each: ``count`` the deconvolved gDNA mass ``f_g·M``, ``mass`` the
    region's total unspliced mass ``M`` (which bounds the achievable density and so fixes the grid top),
    ``eff`` the effective length, ``var`` the belief's ``Var(log f_g)``, and ``anchor`` the zero-mass
    structural regions (see :func:`_reliability`). Substrate selection is the caller's job — it needs the chain.

    The estimator is a weighted sum of zero-native per-region kernels at the population resolution; there is
    no EM, no competition between components and no iteration, so it is deterministic and every region's
    contribution is traceable. A capture-enriched minority therefore cannot be competed away, which is the
    failure the δ-pin predecessor had.
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

    grid = _grid(mass, eff)
    centres = np.clip(np.log10(np.maximum(count, 1.0)) - np.log10(eff), grid[0], grid[-1])
    widths = knn_widths(centres, float(grid[1] - grid[0]), knn_scale)
    weights = _reliability(count, var, anchor)
    density = _render(_poisson_kernels(count, eff, grid), weights, widths, grid)
    total = float(density.sum())
    if not (total > 0.0 and np.isfinite(total)):
        return None

    # ONE PSEUDO-REGION OF COMPLETE IGNORANCE, spread uniformly. Ordinary Laplace smoothing, and the "one" is
    # one region — no constant. It is what bounds how hard this prior can push: without it the empty cells
    # take whatever floor the arithmetic happens to underflow to, which is an assertion the sample cannot
    # support. A population of `W` weighted regions cannot resolve a cell rarer than ~1/W, so the log-range is
    # bounded by log W and the prior stays WEAK AND CORRECTABLE — the governing principle for pass-0 output,
    # which is exactly what this is fitted from.
    density = (density + total / (grid.size * max(float(weights.sum()), 1.0))) / total
    return DensityLandscape(
        log_rho=grid * _LN10,
        logP=np.log(density / density.sum()),
        n_train=int(live.sum()),
        strength=float(strength),
    )
