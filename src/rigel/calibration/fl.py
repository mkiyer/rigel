"""The gDNA and RNA fragment-length laws, and each law's sampling variance.

Builds the library-wide length laws from the accumulator's five length pools (the axis ``rigel.scan_payload`` owns):
two contained gDNA pools (intergenic, intronic; dominant off capture), two crossing gDNA pools ({intron, exon},
{intergenic, exon}; dominant under capture) and the certified RNA pool (an annotated splice). Calibration's effective
lengths read the laws (gDNA opportunities on the gDNA law, RNA opportunities on the RNA law), and so does the EM's
scorer; this is not a per-fragment length likelihood. RNA's law is read from its own spliced counts
(:func:`own_count_law`). gDNA's is the uniform-frame law, what a gDNA fragment looks like before capture selects it,
recovered by a two-pool contrast (:func:`deconvolved_gdna_counts`) and empirical-Bayes-shrunk toward the all-fragment
law, so a thin pool falls back to that anchor without a threshold. No gDNA pool is pure, and nothing here may assume
one is. The EM's two length tables are decided later (:mod:`rigel.calibration.shared_length_law`), because that
decision needs calibration's capture efficiencies. Layer 2, opportunity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..config import CONSTANTS
from .effective_length import observable_contained_eff_length
from .gdna_density import one_sided_rate
from .region_arrays import observable_length
from .signature import RegionType, coarse_type_array
from .sj_opportunity import detilt_pool
from ..scan_payload import (
    N_FRAGMENT_POOLS,
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)

if TYPE_CHECKING:
    from ..frag_length_model import FragmentLengthModel
    from ..scan_payload import AccumulatorPayload
    from .gdna_opportunity import GdnaOpportunity
    from .region_arrays import RegionArrays

__all__ = [
    "FLModels",
    "build_fl_models",
    "gdna_fl_mass",
    "rna_fl_mass",
]

#: gDNA contained in exactly one intergenic or intronic region. Dominant off capture.
_GDNA_CONTAINED_POOLS = (POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC)

#: gDNA crossing exactly one boundary whose flanks are {intron, exon} or {intergenic, exon}. Dominant under capture,
#: where a fragment beside a probe reaches the exon boundary. gDNA-dominated, not pure: RNA crosses such a boundary
#: where it has not spliced or runs past its annotation.
_GDNA_CROSSING_POOLS = (POOL_DNA_INTRON_EXON, POOL_DNA_INTERGENIC_EXON)

#: All four, in ``rigel.scan_payload`` pool order so they pair 1:1 with ``GdnaOpportunity.pools``.
_GDNA_POOLS = _GDNA_CONTAINED_POOLS + _GDNA_CROSSING_POOLS

#: The certified RNA pool: a splice across an annotated sj on the fragment's one surviving path.
_RNA_POOLS = (POOL_RNA_SPLICED,)


@dataclass(frozen=True, slots=True)
class FLModels:
    """Library-wide FL distributions (float64[max_size + 1]) and their pool totals.

    * ``*_pmf`` — the laws for scoring and calibration: RNA's from its own counts (:func:`own_count_law`), gDNA's
      EB-smoothed toward ``global_pmf`` (:func:`build_fl_models`).
    * ``*_counts`` — the raw histograms they were built from (RNA's de-tilted), surfaced as QC through
      :meth:`rna_model`, :meth:`gdna_model` and :meth:`global_model`. ``gdna_counts`` is the uniform-frame histogram:
      the two-pool contrast, or the de-tilted four-pool sum where the contrast declines.

    ``rna_pmf_var`` and ``gdna_pmf_var`` are each pmf's per-length sampling variance, which
    :mod:`rigel.calibration.shared_length_law` reads to tell a difference from noise. RNA's is the own-count law's
    posterior variance at its fitted smoothing (zero when ``n_rna`` is zero, which is no precision); gDNA's is its raw
    law's plug-in variance times the square of the empirical-Bayes weight ``n / (n + prior_ess)``, and 0 at a length
    with no count.
    """

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # the spliced pool's own counts, read as counts (own_count_law)
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    global_counts: np.ndarray  # the accumulator's deposited_lengths — every deposited fragment at L
    rna_counts: np.ndarray  # raw spliced-annotated histogram
    gdna_counts: np.ndarray  # the uniform-frame gDNA histogram (see above)
    n_rna: float
    n_gdna: float
    max_size: int
    #: float64[N_FRAGMENT_POOLS, max_size + 1] — the five pools unaggregated, straight off the payload, so the report
    #: can show each pool's length law on its own.
    pool_counts: np.ndarray
    rna_pmf_var: np.ndarray
    gdna_pmf_var: np.ndarray

    def _empirical(self, counts: np.ndarray) -> "FragmentLengthModel":
        """Wrap a raw count vector as an unfinalized ``FragmentLengthModel``, so its summary statistics read the
        raw histogram with no scoring-side smoothing, like the scanner's QC models."""
        from ..frag_length_model import FragmentLengthModel

        return FragmentLengthModel(max_size=self.max_size, counts=counts.copy())

    def global_model(self) -> "FragmentLengthModel":
        """Empirical global FL distribution as a ``FragmentLengthModel`` (QC)."""
        return self._empirical(self.global_counts)

    def rna_model(self) -> "FragmentLengthModel":
        """Empirical RNA (spliced-annotated) FL distribution (QC)."""
        return self._empirical(self.rna_counts)

    def gdna_model(self) -> "FragmentLengthModel":
        """Empirical gDNA FL distribution (QC): ``gdna_counts``, the histogram the gDNA pmf is shrunk
        from."""
        return self._empirical(self.gdna_counts)


def _pool_sum(payload: "AccumulatorPayload", pools) -> np.ndarray:
    """Sum the named rows of ``payload.pool_lengths`` into one float64 histogram over ``L``.

    Binned at ``L``, the molecule length, never the covered length, which would collapse gDNA to a spike at twice the
    read length. The accumulator's ``L`` includes the mate gap and excludes spliced-out introns, so it is the molecule
    length for both origins.
    """
    return np.asarray(payload.pool_lengths, dtype=np.float64)[list(pools)].sum(axis=0)


def gdna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The gDNA length histogram: all four structural gDNA pools, summed.

    Meaningful only with :meth:`GdnaOpportunity.combined_probability` as its divisor: the pools tilt in opposite
    directions, so the raw sum is biased long. The two-pool contrast's fallback (:func:`build_fl_models`).
    """
    return _pool_sum(payload, _GDNA_POOLS)


def rna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The certified RNA length histogram: fragments that used an annotated sj on their one surviving
    path."""
    return _pool_sum(payload, _RNA_POOLS)


def resolution_weight(signal_sq: float, noise_sq: float) -> float:
    """``S / (S + N)``: how much of an observed split is signal rather than its own sampling noise.

    It replaces a ``|split| > standard error`` threshold with the same comparison made continuously: 0 for pure noise,
    1 when the noise vanishes, ½ where the threshold would sit. Both arguments are variances the data supply, so no
    constant enters.
    """
    s = max(float(signal_sq), 0.0)
    n = max(float(noise_sq), 0.0)
    if s <= 0.0:
        return 0.0
    if n <= 0.0:
        return 1.0
    return s / (s + n)


def _region_counts(payload, regions: "RegionArrays") -> np.ndarray:
    """``float64[R]`` — each region's contained count, genome-strand columns summed: gDNA is strand-symmetric and these
    estimators read density, so they also work on an unstranded library. Refused unless *regions* is the partition the
    payload was deposited into."""
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)
    counts = counts.sum(axis=1) if counts.ndim > 1 else counts
    if counts.shape[0] != regions.n_regions:
        raise ValueError(
            f"the frame has {regions.n_regions:,} regions and the payload {counts.shape[0]:,}: they address "
            "different partitions. Pass RegionArrays.from_index of the index the scan used."
        )
    return counts


def contained_purities(payload: "AccumulatorPayload", regions: "RegionArrays"):
    """The one-sided gDNA rate over the two contained pools' observable support, and each pool's gDNA share at it.

    Supports are summed over each region's observable pieces (``region_arrays.observable_pieces``); a region masked
    whole carries no exposure, since over its full length it would look like the cleanest object of all and collapse
    the purities. Returns ``(rate, (a_0, a_1), (n_0, n_1))`` for the intergenic and intronic pools, a pool of count
    ``n_p`` and support ``E_p`` having share ``a_p = min(rate·E_p / n_p, 1)`` (0 when empty), or ``None`` where the
    rate found no support, as in a zero-gDNA library."""
    types = coarse_type_array(regions.signature)
    counts = _region_counts(payload, regions)
    observable = observable_length(regions)

    is_ig = (types == int(RegionType.INTERGENIC)) & (observable > 0.0)
    is_in = (types == int(RegionType.INTRON)) & (observable > 0.0)
    # the contained pools' own law prices the regions they come from
    support = observable_contained_eff_length(regions, _pool_sum(payload, _GDNA_CONTAINED_POOLS))
    fit = one_sided_rate(counts[is_ig | is_in], support[is_ig | is_in])
    if not fit.informative:
        return None

    weights, totals = [], []
    for mask in (is_ig, is_in):
        n_p = float(counts[mask].sum())
        e_p = float(support[mask].sum())
        weights.append(min(fit.rate * e_p / n_p, 1.0) if n_p > 0.0 else 0.0)
        totals.append(n_p)
    return fit.rate, (weights[0], weights[1]), (totals[0], totals[1])


def deconvolved_gdna_counts(
    payload: "AccumulatorPayload",
    gdna_opportunity: "GdnaOpportunity",
    regions: "RegionArrays",
) -> "tuple[np.ndarray | None, np.ndarray | None]":
    """The two-pool contrast: the gDNA length histogram with the contaminant divided out.

    No contained pool is pure (the intronic pool carries RNA inside introns, the intergenic pool unannotated
    transcription, biasing a pool's mean by ``RNA_share × (len_RNA − len_gDNA)``), but the two share one opportunity
    geometry, since RNA there is genomically contiguous as gDNA is. After de-tilting each is the same two shapes at a
    different gDNA share ``a_p = rho·E_p/n_p`` (:func:`contained_purities`), so the contaminant cancels::

        f_0 = a_0*g + (1-a_0)*r        f_1 = a_1*g + (1-a_1)*r
        =>  g = [ (1-a_1)*f_0 - (1-a_0)*f_1 ] / (a_0 - a_1)

    The divisor is the separation of the purities, not a purity, so it stays finite as a pool gets dirty. The
    inversion is blended toward the pools' own mixture by the resolution weight ``sep² / (sep² + SE²)``, with
    ``SE² ~ a_0²/n_0 + a_1²/n_1`` by the delta method (``Var(n_p) ~ n_p``; the shared ``rho`` is common-mode and
    cancels): the estimator fades out at a near-pure library, and identical purities return the mixture itself, so the
    law is continuous through a purity tie.

    Under hybrid capture the shared-contaminant premise is false, and the contrast is safe only because capture
    depletes the intergenic pool rather than making it impure: ``a_0`` exceeds 1 and clips, and at ``a_0 = 1`` the
    formula is ``g = f_0``. The intronic pool's coefficient is then zero, so the term that needed the premise is gone.
    A probe panel that put RNA back into intergenic space would break this silently.

    Returns ``(counts, variance)``, both ``None`` where it declines and the caller falls back to the four-pool sum:
    with no gDNA density (a zero-gDNA library has no gDNA law to estimate), an empty pool, or nothing left after the
    projection back onto the cone. ``variance`` is the normalised law's per-length plug-in variance through the
    contrast's own coefficients: each de-tilted pool contributes ``f_p² / c_p`` (de-tilting rescales a length and adds
    no evidence), combined coefficient by coefficient so their covariance is included.
    """
    purities = contained_purities(payload, regions)
    if purities is None:  # no gDNA density
        return None, None
    _rate, (a0, a1), totals = purities
    if totals[0] <= 0.0 or totals[1] <= 0.0:  # an empty pool
        return None, None

    sep = a0 - a1
    se = float(np.sqrt(a0 * a0 / totals[0] + a1 * a1 / totals[1]))
    # blend the inversion toward the mixture by how well the separation is resolved: lam = 1 accepts it whole, lam = 0
    # declines it, continuously between
    lam_sep = resolution_weight(sep * sep, se * se)

    # each contained pool takes its own opportunity, never the combined divisor of the four-pool sum
    total = np.asarray(gdna_opportunity.total, dtype=np.float64)
    raw = np.asarray(payload.pool_lengths, dtype=np.float64)
    f = []
    for pool in _GDNA_CONTAINED_POOLS:
        opp = np.asarray(gdna_opportunity.pools[pool], dtype=np.float64)
        prob = np.zeros_like(opp)
        np.divide(opp, total, out=prob, where=total > 0.0)
        f.append(_normalized(detilt_pool(raw[pool], prob)))
    mixture = _normalized(totals[0] * f[0] + totals[1] * f[1])
    # the law's per-length coefficients on the two pools, for its variance: the mixture's weight everywhere, the
    # inverse's only where the clip leaves it open, through the renormalisation
    w0, w1 = totals[0] / (totals[0] + totals[1]), totals[1] / (totals[0] + totals[1])
    c0, c1 = np.full(f[0].size, w0), np.full(f[0].size, w1)
    if sep == 0.0:
        # the fade's own limit: at zero separation the resolution weight is zero
        g = mixture
    else:
        inverse = ((1.0 - a1) * f[0] - (1.0 - a0) * f[1]) / sep
        # negative excursions are sampling noise on a density; clipping is the cheapest projection back onto the cone
        g = np.clip(inverse, 0.0, None)
        s = g.sum()
        if s > 0.0:
            open_ = (inverse > 0.0) * (lam_sep / (s * sep))
            c0 = (1.0 - lam_sep) * w0 + open_ * (1.0 - a1)
            c1 = (1.0 - lam_sep) * w1 - open_ * (1.0 - a0)
        g = mixture + lam_sep * (_normalized(g) - mixture) if s > 0.0 else mixture
    g = np.clip(g, 0.0, None)
    if not g.sum() > 0.0:  # empty after the contrast
        return None, None
    variance = c0 * c0 * _plugin_variance(
        f[0], raw[_GDNA_CONTAINED_POOLS[0]]
    ) + c1 * c1 * _plugin_variance(f[1], raw[_GDNA_CONTAINED_POOLS[1]])
    # rescaled to the pool mass the four-pool path would carry, because the EB shrinkage reads the total as the
    # evidence behind the shape
    return _normalized(g) * float(raw[list(_GDNA_POOLS)].sum()), variance


def _aligned(counts: np.ndarray, max_size: int) -> np.ndarray:
    """Align a raw FL count vector to ``max_size + 1`` bins (overflow → last bin)."""
    out = np.zeros(max_size + 1, dtype=np.float64)
    c = np.asarray(counts, dtype=np.float64)
    n = min(out.size, c.size)
    out[:n] = c[:n]
    if c.size > out.size:
        out[max_size] += float(np.sum(c[out.size :], dtype=np.float64))
    return out


def _normalized(v: np.ndarray) -> np.ndarray:
    s = float(v.sum())
    return v / s if s > 0.0 else np.full(v.size, 1.0 / v.size, dtype=np.float64)


def _plugin_variance(law: np.ndarray, raw: np.ndarray) -> np.ndarray:
    """A normalised law's per-length plug-in variance ``law² / raw``: a length's relative variance is one over the raw
    count behind it, which de-tilting rescales and does not add to. A length with no raw count carries none."""
    law = np.asarray(law, dtype=np.float64)
    c = _fitted(raw, law.size)
    return np.divide(law * law, c, out=np.zeros(law.size), where=c > 0.0)


def _fitted(raw: np.ndarray, size: int) -> np.ndarray:
    """``raw`` as a float vector of exactly ``size`` entries, zero-padded or cut."""
    out = np.zeros(size, dtype=np.float64)
    v = np.asarray(raw, dtype=np.float64)
    n = min(size, v.size)
    out[:n] = v[:n]
    return out


def _walk_hessian(mu: np.ndarray, step: float) -> np.ndarray:
    """``diag(mu) + K / step`` in upper banded form, ``K`` the first-difference precision of a path over the lengths."""
    n = mu.size
    degree = np.full(n, 2.0)
    degree[0] = degree[-1] = 1.0
    if n == 1:
        degree[0] = 0.0
    ab = np.zeros((2, n), dtype=np.float64)
    ab[1] = mu + degree / step
    ab[0, 1:] = -1.0 / step
    return ab


def own_count_law(
    raw: np.ndarray, selection: np.ndarray | None, anchor: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """An origin's length law from its own counts, read as counts: ``(pmf, per-length variance, step)``.

    The pool's count at length ``l`` is Poisson with mean ``N · s(l) · law(l)``, ``s`` the pool's selection: for RNA's
    spliced pool the chance a length-``l`` fragment crosses an annotated junction, the divisor
    :func:`~rigel.calibration.sj_opportunity.detilt_pool` de-tilts by (``None`` is no selection). The law is the
    all-fragment law times a smooth deviation, ``law = m · e^d / Z``, ``d`` a first-order random walk over the lengths
    ``m`` holds, its step variance the Laplace-evidence maximum within ``CONSTANTS.fragment_length.law_step_min`` and
    ``law_step_max``. So:

    * a length with counts reads the pool's own law, smoothed at the scale the data support;
    * a length the pool could reach and did not (``s > 0``, no count) is evidence that the law is low there;
    * a length it cannot reach (``s = 0``) carries none, and ``d`` holds its neighbour's value there: ``m``'s shape;
    * an empty pool leaves ``d`` at zero, so the law is ``m`` exactly, with no pseudo-count.

    ``m`` enters only where the pool says nothing and is never weighed as evidence about the pool's origin, as a
    pseudo-count toward it would be (a minority origin's law would take the majority's shape). The variance is the
    normalised law's posterior variance at the fitted step, by the delta method through ``J = diag(law) − law lawᵀ``:
    ``Var law(l) = law(l)² (e_l − law)ᵀ H⁻¹ (e_l − law)``, ``H`` the posterior Hessian in ``d``. The deviation's level
    cancels in it, so a one-length law carries none, and each term is a squared norm through ``H``'s Cholesky factor,
    so none is negative. An empty pool returns ``m`` with zero variance and step 0: a fallback no count measured, which
    the length decision reads as no RNA measurement, never as a precise law.

    Gate: ``tests/calibration/test_fl.py``.
    """
    m = np.asarray(anchor, dtype=np.float64)
    size = m.size
    counts = _fitted(raw, size)
    if np.any((counts > 0.0) & (m <= 0.0)):
        raise ValueError(
            "a pool count at a length the anchor holds no fragment at: the anchor is every deposited fragment and the "
            "pool is drawn from it, so the two describe different payloads"
        )
    reach = np.ones(size) if selection is None else _fitted(selection, size)
    held = np.flatnonzero(m > 0.0)
    total = float(counts.sum())
    if held.size == 0 or total <= 0.0:
        return m.copy(), np.zeros(size, dtype=np.float64), 0.0
    lo, hi = int(held[0]), int(held[-1]) + 1
    live = (reach[lo:hi] > 0.0) & (m[lo:hi] > 0.0)
    offset = np.zeros(hi - lo)
    offset[live] = np.log(total * reach[lo:hi][live] * m[lo:hi][live])
    step, d, factor = _fit_step(_Deviation(counts[lo:hi], live, offset))
    law = np.zeros(size, dtype=np.float64)
    var = np.zeros(size, dtype=np.float64)
    law[lo:hi], var[lo:hi] = _law_and_variance(m[lo:hi], d, factor)
    return law, var, step


#: The log of the largest double: the clip that keeps a Poisson mean ``e^eta`` finite.
_LOG_MAX = float(np.log(np.finfo(np.float64).max))


@dataclass(frozen=True)
class _Deviation:
    """The deviation ``d`` of a pool's law from the anchor over the anchor's held lengths: the pool's Poisson
    likelihood (``c`` its counts, ``offset`` the log mean at ``d = 0``, ``live`` where the pool can reach a length)
    under a first-order random walk over the lengths, of step variance ``step``."""

    c: np.ndarray
    live: np.ndarray
    offset: np.ndarray

    def objective(self, d: np.ndarray, step: float) -> tuple[float, np.ndarray]:
        """The negative log posterior at ``d``, up to a constant, and the Poisson means."""
        eta = np.clip(self.offset + d, -_LOG_MAX, _LOG_MAX)
        mu = np.where(self.live, np.exp(eta), 0.0)
        penalty = float(np.sum(np.diff(d) ** 2)) / (2.0 * step)
        return float(np.sum(mu - np.where(self.live, self.c * eta, 0.0))) + penalty, mu

    def mode(self, step: float, d: np.ndarray):
        """The deviation's mode at one step variance: Newton with step halving on a convex objective."""
        from scipy.linalg import cho_solve_banded, cholesky_banded

        k = CONSTANTS.fragment_length
        f, mu = self.objective(d, step)
        for _ in range(k.law_newton_steps):
            walk = np.zeros_like(d)
            diff = np.diff(d)
            walk[:-1] -= diff
            walk[1:] += diff
            grad = (mu - np.where(self.live, self.c, 0.0)) + walk / step
            direction = cho_solve_banded((cholesky_banded(_walk_hessian(mu, step)), False), -grad)
            scale = 1.0
            while True:
                trial = d + scale * direction
                f_trial, mu_trial = self.objective(trial, step)
                # halve until the objective falls, or until the step no longer moves the iterate at all
                if f_trial <= f or scale * np.abs(direction).max() <= np.finfo(np.float64).eps * (
                    1.0 + np.abs(d).max()
                ):
                    break
                scale *= 0.5
            if f_trial > f:
                break  # no step lowers the objective: the mode, to rounding
            change = f - f_trial
            d, f, mu = trial, f_trial, mu_trial
            if change <= k.law_newton_tolerance * max(1.0, abs(f)):
                break
        return d, f, mu

    def evidence(self, log_step: float, d0: np.ndarray):
        """``(log evidence, mode, Hessian factor)`` at step variance ``exp(log_step)``, Newton started from ``d0``."""
        from scipy.linalg import cholesky_banded

        step = float(np.exp(log_step))
        d, f, mu = self.mode(step, d0)
        factor = cholesky_banded(_walk_hessian(mu, step))
        width = d.size
        log_det = 2.0 * float(np.log(factor[1]).sum())
        # the path's first-difference precision has pseudo-determinant `width`; its constant mode is the likelihood's
        log_prior = 0.5 * (np.log(width) - (width - 1) * np.log(step))
        return -f + log_prior - 0.5 * log_det, d, factor


def _fit_step(walk: _Deviation) -> tuple[float, np.ndarray, np.ndarray]:
    """``(step, mode, Hessian factor)`` at the step variance of maximal Laplace evidence, searched between
    ``CONSTANTS.fragment_length.law_step_min`` and ``law_step_max``. Each evaluation starts Newton from the last one's
    mode."""
    from scipy.optimize import minimize_scalar

    k = CONSTANTS.fragment_length
    state = {"d": np.zeros(walk.c.size)}

    def evidence(log_step: float):
        value, d, factor = walk.evidence(log_step, state["d"])
        state["d"] = d
        return value, d, factor

    search = minimize_scalar(
        lambda x: -evidence(x)[0],
        bounds=(np.log(k.law_step_min), np.log(k.law_step_max)),
        method="bounded",
        options={"xatol": k.law_step_tolerance},
    )
    step = float(np.exp(search.x))
    _, d, factor = evidence(float(search.x))
    return step, d, factor


def _law_and_variance(m: np.ndarray, d: np.ndarray, factor: np.ndarray):
    """The normalised law ``m · e^d / Z`` and its per-length delta-method variance through the Hessian's factor."""
    from scipy.linalg import solve_banded

    w = m * np.exp(np.clip(d - d.max(), -_LOG_MAX, 0.0))
    f = w / w.sum()
    width = f.size
    if width == 1:
        # one length holds probability one, whatever the factor: no variance (the solve below can read it as 0/0)
        return f, np.zeros(1)
    # H = UᵀU with U the upper banded factor, so (e_l − f)ᵀ H⁻¹ (e_l − f) = |U⁻ᵀ (e_l − f)|², a block of lengths per
    # solve: the right-hand side is dense, and the whole support at once would hold width² doubles. The expansion
    # diag(H⁻¹) − 2 H⁻¹f + fᵀH⁻¹f is linear in the width but cancels where the law is near a point mass; the squared
    # norm never goes negative.
    lower = np.zeros((2, width), dtype=np.float64)
    lower[0] = factor[1]
    lower[1, :-1] = factor[0, 1:]
    quadratic = np.empty(width, dtype=np.float64)
    block = CONSTANTS.fragment_length.law_variance_block
    for first in range(0, width, block):
        last = min(width, first + block)
        rhs = np.repeat(-f[:, None], last - first, axis=1)
        rhs[np.arange(first, last), np.arange(last - first)] += 1.0
        whitened = solve_banded((1, 0), lower, rhs)
        quadratic[first:last] = np.einsum("ij,ij->j", whitened, whitened)
    return f, f * f * quadratic


def _smooth_eb(aligned: np.ndarray, global_pmf: np.ndarray, prior_ess: float):
    """Smooth EB pmf: ``(counts + prior_ess·global_pmf) / (total + prior_ess)``.

    ``aligned`` is an already ``max_size``-aligned count vector. Continuous in the
    pool total — no quality threshold. Returns ``(pmf, pool_total)``.
    """
    total = float(aligned.sum())
    denom = total + prior_ess
    pmf = (aligned + prior_ess * global_pmf) / denom if denom > 0.0 else global_pmf.copy()
    return _normalized(pmf), total


def build_fl_models(
    payload: "AccumulatorPayload",
    *,
    sj_opportunity: np.ndarray,
    gdna_opportunity: "GdnaOpportunity",
    regions: "RegionArrays",
) -> FLModels:
    """Build the global, RNA and gDNA FL pmfs from one payload, in one frame.

    All three histograms come off the same payload, so they agree on what a fragment length is. The anchor is
    ``payload.deposited_lengths``, every deposited fragment at its own ``L``, and the component pools
    (:func:`rna_fl_mass`, :func:`gdna_fl_mass`) are drawn from exactly that population. It is unconditional given
    deposit: it excludes what the accumulator rejects (too long, ambiguous path, strand-undefined, empty), each counted
    in ``payload.qc``. No outside anchor can be passed: a scanner-frame anchor measures length by other rules, and
    shrinking toward it would read a ruler mismatch as composition.

    RNA's law is its own counts' at their sj opportunity (:func:`own_count_law`); gDNA's is EB-shrunk toward the anchor
    with one Dirichlet ``CONSTANTS.fragment_length.pool_prior_ess``. Each pool is divided by its own opportunity,
    because the two selections differ:

    * ``sj_opportunity`` — ``pi(w)``, the chance a uniformly placed length-``w`` fragment crosses an annotated sj
      (:mod:`rigel.calibration.sj_opportunity`); longer fragments do so more often.
    * ``gdna_opportunity`` — the four gDNA pools' opportunities and the reference total
      (:mod:`rigel.calibration.gdna_opportunity`). A contained pool's opportunity falls with length and a crossing
      pool's rises, so one divisor over the pooled sum reads the gDNA mean too long.

    gDNA's histogram is the two-pool contrast (:func:`deconvolved_gdna_counts`) over ``regions``, the index's frame
    (:meth:`~rigel.calibration.region_arrays.RegionArrays.from_index`): the partition the scanner deposits into,
    carrying the unobservable mask every contained support is read over. Where the contrast declines it is the four
    pools de-tilted by their combined opportunity. :func:`fl_models_from_histograms` is the kernel over three free
    histograms, for tests.
    """
    gdna_counts, gdna_var = deconvolved_gdna_counts(payload, gdna_opportunity, regions)
    if gdna_counts is None:
        gdna_raw = gdna_fl_mass(payload)
        gdna_counts = detilt_pool(gdna_raw, gdna_opportunity.combined_probability())
        gdna_var = _plugin_variance(_normalized(gdna_counts), gdna_raw)
    return fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=rna_fl_mass(payload),
        rna_selection=sj_opportunity,
        gdna_counts=gdna_counts,
        max_size=int(payload.max_length),
        pool_counts=payload.pool_lengths,
        gdna_var=gdna_var,
    )


def fl_models_from_histograms(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    pool_counts: np.ndarray | None = None,
    rna_selection: np.ndarray | None = None,
    gdna_var: np.ndarray | None = None,
) -> FLModels:
    """The kernel: three histograms in, three pmfs out, with each component pmf's sampling variance.

    Production reaches it only through :func:`build_fl_models`; it exists so tests can exercise the laws over anchors
    and pools no real payload produces. ``rna_counts`` are raw counts at their selection ``rna_selection`` (``None``:
    no selection), read by :func:`own_count_law`; ``FLModels.rna_counts`` keeps them de-tilted for the report.
    ``gdna_counts`` is gDNA's histogram, shrunk toward the anchor, and ``gdna_var`` its normalised law's per-length
    variance (:func:`build_fl_models` passes the contrast's own); with none, it is read as raw counts.
    """
    global_aligned = _aligned(global_counts, max_size)
    rna_raw = _aligned(rna_counts, max_size)
    # the report's histogram, de-tilted with its pool total preserved
    rna_aligned = rna_raw if rna_selection is None else detilt_pool(rna_raw, rna_selection)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, rna_var, _step = own_count_law(rna_raw, rna_selection, global_pmf)
    n_rna = float(rna_raw.sum())
    prior_ess = CONSTANTS.fragment_length.pool_prior_ess
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)

    # gDNA's law's variance: its raw law's, at the square of the empirical-Bayes weight its counts carry
    raw_var = (
        _plugin_variance(_normalized(gdna_aligned), gdna_aligned)
        if gdna_var is None
        else _aligned(gdna_var, max_size)
    )
    weight = n_gdna / (n_gdna + prior_ess) if n_gdna + prior_ess > 0.0 else 0.0

    return FLModels(
        global_pmf=global_pmf,
        rna_pmf=rna_pmf,
        gdna_pmf=gdna_pmf,
        global_counts=global_aligned,
        rna_counts=rna_aligned,
        gdna_counts=gdna_aligned,
        pool_counts=(
            np.zeros((N_FRAGMENT_POOLS, max_size + 1), dtype=np.float64)
            if pool_counts is None
            else np.asarray(pool_counts, dtype=np.float64)
        ),
        n_rna=n_rna,
        n_gdna=n_gdna,
        max_size=int(max_size),
        rna_pmf_var=rna_var,
        gdna_pmf_var=weight * weight * raw_var,
    )
