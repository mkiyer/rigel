"""Strand Beta-Binomial overdispersion fits for both components, gDNA and RNA.

Both strand splits are modelled as Beta-Binomial, ``K⁺ | N ~ BetaBinom(N, mean, overdispersion)``, with
the intra-class correlation in ``[0, 1)`` (``0`` being Binomial). The mean differs by component — gDNA is
unstranded, biologically fixed at ½, RNA is skewed at ``rna_sense_frac`` = κ from the spliced channel —
but the overdispersion machinery is shared: a per-region or per-boundary sense split is more spread out
than Binomial (local capture bias, PCR, pair-anchoring noise), and that excess variance is the
overdispersion.

Both components carry a fitted overdispersion, and both must. Forcing RNA Binomial while gDNA is
Beta-Binomial makes the strand likelihood spuriously informative on unstranded data, because the
``−½·log var`` term prefers the lower-variance component and pulls balanced regions toward RNA.

The two are fit by different estimators, because the RNA side has a certified pure population to fit
from and the gDNA side has none: RNA by the pooled moment over certified spliced seeds
(:func:`fit_rna_strand_from_sj_table`, the per-sj SJ strand table whose marginal is κ itself), gDNA by
the raw away-half moment below.

Breaking the circularity — the away-half estimator. Fitting the gDNA overdispersion needs seeds whose
RNA content does not masquerade as gDNA strand spread, and a fit that picks a structural class and
asserts it pure cannot supply them: purity is a property of the annotation and the sample, not of the
genome. Pervasive transcription is real, the
intergenic space is whatever the user's GTF leaves over, and most genes are off in any one sample but
nobody knows which; on a contig fed unannotated transcripts every purity-based fit moves and this one
does not.

So the estimator trusts no class. The lemma: on a strand-specific library, RNA of a seed's own gene
pulls its sense fraction only toward κ. Orient ``d = (K − N/2)·sign(½ − κ)`` so that RNA pulls ``d``
negative. Under pure gDNA ``d`` is symmetric about 0 and the moment excess ``d² − N/4`` is an even
function of ``d``, so the pooled method of moments restricted to the away half (``d > 0``, with a tie
``d = 0`` at weight ½, which keeps even ``N`` unbiased) has exactly the same expectation as the full
moment under the null — while a contaminated seed reaches the away side only by noise, with a small
``d``, biasing the estimate down and never up. That is unbiased for ρ_g under any distribution of RNA
content across the seeds, with no weight, no purity class and no reference (:func:`away_half_moment`).

The seed set is every count- and strand-observable genic object in whatever GTF is supplied: intron
regions of single-strand genes, exon|intron boundaries and gene-edge boundaries alike, the lemma
protecting them all. Intergenic regions are out — they have no gene strand to orient ``d`` by, so the
lemma cannot protect them, and they are exactly where unannotated transcription lives — and so are
``TS_AMBIG`` objects, which have no defined sense. The known limit, recorded rather than hidden:
unannotated antisense transcription pushes toward the away side and inflates ρ̂. Half the pairs are
used, so the estimator's null information is half the pair count.

Estimator. For seed ``s`` with oriented ``d_s`` and weight ``a_s = 1[d_s > 0] + ½·1[d_s = 0]``::

    od = clip( Σ a_s·(d_s² − N_s/4)  /  Σ a_s·N_s(N_s − 1)/4 ,  0, _MAX_OVERDISPERSION )

The gDNA fit is that raw moment clipped to ``[0, _MAX_OVERDISPERSION]`` — the ``Beta(2,2)`` ceiling,
od = 0.2 — with no shrinkage target: a target picked rather than measured over-shrinks exactly where
gDNA information is scarce. With no pair in the away half the fit returns the ceiling at zero
information, and :func:`reconcile_overdispersions` hands the component the other one's measured value.

The currency in which the two fits are weighed against each other is information, not seed count, and
that is a units fix rather than a tuning choice: overdispersion is a correlation between fragments in a
region, so a seed with one fragment carries none of it, having nothing to be correlated with. Weighting
by seed count lets a mass of singletons outvote everything on a library that has almost no pairs. See
:func:`_null_information` for the null form and :func:`influence_weights` for the fitted one.

The method of moments is closed-form, ``O(n_seed_regions)``, and uses the same variance decomposition
the deconvolution applies (ψ's ``simplex_logodds._mixture_strand_loglik``, whose two-component reference
is :mod:`strand_likelihood`), so fit and application stay consistent. The ceiling is the only asserted
constant left in the pair of fits, and ``tests/calibration/test_sweep.py`` pins that.

Each component pairs a pure estimator with a thin seed-extraction wrapper, so the estimator itself is
trivially testable: :func:`fit_gdna_strand_overdispersion` / :func:`fit_gdna_strand_from_substrate` for
gDNA (the wrapper needs the region/boundary geometry), and :func:`fit_rna_strand_overdispersion` /
:func:`fit_rna_strand_from_sj_table` for RNA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .region_arrays import boundary_region_indices
from .signature import TS_NEG, TS_NONE, TS_POS

#: Hard floor on overdispersion — the Binomial limit (no negative overdispersion).
_OVERDISPERSION_FLOOR: float = 0.0
#: Symmetric ``Beta(a, a)`` shape at the most overdispersion allowed (the per-region sense rate is at
#: most this spread out). ``a = 2`` ⇒ ``od = 1/(2·2+1) = 0.2``, the overdispersion ceiling. Every fitted
#: or reconciled overdispersion is clamped to ``[0, _MAX_OVERDISPERSION]``.
_CEIL_ALPHA_BETA: float = 2.0


def overdispersion_for_beta(alpha_beta: float) -> float:
    """Beta-Binomial intra-class correlation (overdispersion) for a symmetric ``Beta(a, a)``.

    ``od = 1 / (2·a + 1)``: ``a = 1`` (uniform) → 1/3; ``a = 2`` → 1/5 = 0.2 (max overdispersion);
    ``a = 3`` → 1/7 ≈ 0.143; ``a → ∞`` → 0 (Binomial). Inverse of ``beta_concentration``.
    """
    return 1.0 / (2.0 * float(alpha_beta) + 1.0)


#: Overdispersion ceiling (most overdispersion allowed) = ``Beta(2, 2)`` ⇒ ``od = 0.2``.
_MAX_OVERDISPERSION: float = overdispersion_for_beta(_CEIL_ALPHA_BETA)


@dataclass(frozen=True, slots=True)
class GdnaStrandModel:
    """Fitted gDNA strand model: the global Beta-Binomial overdispersion + fit provenance."""

    gdna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
    #: The precision of the estimate actually made — the currency in which this fit is weighed against
    #: the RNA fit's (:func:`reconcile_overdispersions`). 0 when the seeds carry no pair.
    information: float
    n_seed_regions: int  # seed regions that carried gDNA-eligible fragments
    n_seed_fragments: int  # total gDNA-eligible fragments across seed regions
    fallback_used: (
        bool  # True ⇒ no gDNA strand signal ⇒ returned the caller's fallback, never a constant
    )
    #: The moment before the clip — what the data actually said. ``nan`` on the fallback path.
    raw_overdispersion: float = float("nan")
    #: True when the raw moment was above the ceiling, so the value carried is a clamp and not a
    #: measurement. A clamp is a legitimate answer to a misspecified model — a genuine intra-class
    #: correlation is depth-invariant, and on real libraries the moment rises with seed depth — but it
    #: has to announce itself rather than be read as a fit.
    clamped_at_ceiling: bool = False
    #: The effective number of seeds behind the estimate (:func:`seed_participation`).
    effective_seeds: float = float("nan")

    def beta_concentration(self) -> float:
        """Symmetric Beta concentration ``a`` of ``Beta(a, a)`` for the gDNA sense rate.

        ``a = ½·(1 − od)/od``; larger ``a`` ⇒ tighter about ½ (less overdispersion). Returns
        ``+inf`` at ``od = 0`` (the Binomial limit — a point mass at ½).
        """
        od = self.gdna_strand_overdispersion
        return float("inf") if od <= 0.0 else 0.5 * (1.0 - od) / od


@dataclass(frozen=True, slots=True)
class RnaStrandModel:
    """Fitted RNA strand model: the global Beta-Binomial overdispersion of the spliced (pure-RNA)
    strand split + fit provenance. The RNA sense *mean* (``rna_sense_frac``) lives in
    :class:`strand_balance.StrandBalance`; this carries only the between-SJ overdispersion."""

    rna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
    #: The moment BEFORE the clip — what the spliced seeds actually said. ``nan`` with no pair.
    raw_overdispersion: float
    #: The estimator's null information, comparable with :attr:`GdnaStrandModel.information`.
    information: float
    n_seed_regions: int  # splice junctions that carried strand-qualified fragments
    n_seed_fragments: (
        int  # total strand-qualified fragments (exactly one per fragment, no double-count)
    )
    fallback_used: bool  # True ⇒ no spliced strand signal ⇒ returned the prior overdispersion

    def beta_concentration(self) -> float:
        """Symmetric Beta concentration ``a`` of ``Beta(a, a)`` for the RNA sense rate (about κ)."""
        od = self.rna_strand_overdispersion
        return float("inf") if od <= 0.0 else 0.5 * (1.0 - od) / od


def _null_information(comp_frags: np.ndarray, comp_var: np.ndarray | float) -> float:
    """The estimator's null information ``I = 1/Var(od_mom)|₀`` — exact, not asymptotic.

    ``od_mom = Σe_s / Σc_s`` with ``c_s = n_c(n_c−1)·pq`` deterministic, so ``Var(od_mom) = Var(Σe_s)/(Σc_s)²``.
    Under the null (``od = 0``) the sense count is Binomial, and with ``X = K − n·μ``,
    ``E[X²] = n·pq`` and ``E[X⁴] = 3(n·pq)² + n·pq(1 − 6pq)``, so ``Var(e_s) = 2(n·pq)² + n·pq(1 − 6pq)``:

        I = (Σ n(n−1)·pq)² / Σ [ 2n²pq² + n·pq − 6n·pq² ]

    At μ = ½ this reduces exactly to the pair count ``Σ n(n−1)/2``, the intuition being that
    overdispersion is a correlation between fragments, so the unit of evidence is a pair: a seed of one
    fragment contributes nothing and a seed of ``n`` contributes ``n(n−1)/2``. Verified algebraically and
    by Monte Carlo.

    ⛔ Away from μ = ½ the pair count is wrong and must not be substituted — at the RNA fit's κ it
    overstates the information by an order of magnitude. That is why this is the general form rather than
    a pair count, and why the two components share one code path. Gate:
    ``tests/calibration/test_sweep.py``'s pair-count-at-symmetric-mean case.
    """
    n = np.asarray(comp_frags, dtype=np.float64)
    pq = np.asarray(comp_var, dtype=np.float64)
    num = float(np.sum(np.maximum(n * (n - 1.0), 0.0) * pq)) ** 2
    den = float(np.sum(2.0 * n * n * pq * pq + n * pq - 6.0 * n * pq * pq))
    return num / den if den > 0.0 else 0.0


def _away_half_parts(
    sense: np.ndarray, total: np.ndarray, rna_sense_frac: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(a_s, e_s, c_s, n_s)`` — the away-half weight, moment excess, pair scale and seed depth.

    ``d = (K − N/2)·sign(½ − κ)`` so RNA of the seed's own gene pulls ``d`` negative; ``a_s`` is 1 on the
    away half, ½ on a tie, 0 otherwise; ``e_s = d² − N/4``; ``c_s = N(N−1)/4``. At ``κ = ½`` the orientation
    is degenerate and every seed enters at weight 1 — :func:`away_half_moment` says why.
    """
    k = np.asarray(sense, dtype=np.float64)
    n = np.asarray(total, dtype=np.float64)
    valid = n > 0.0
    orient = float(np.sign(0.5 - float(rna_sense_frac)))
    d = (k - 0.5 * n) * (orient if orient != 0.0 else 1.0)
    if orient == 0.0:  # κ = ½ — RNA is symmetric; the full two-sided moment
        a = valid.astype(np.float64)
    else:
        a = np.where(d > 0.0, 1.0, np.where(d == 0.0, 0.5, 0.0)) * valid
    return a, d * d - 0.25 * n, np.maximum(n * (n - 1.0), 0.0) * 0.25, n


def between_seed_variance(overdispersion: float, mean: float = 0.5) -> float:
    """``V∞(ρ, μ)`` — the variance one seed's ``od̂_s`` retains at UNLIMITED depth.

        V∞(ρ, μ) = 3ρ²·[2ρ + μ(1−μ)(1 − 7ρ)] / [μ(1−μ)(1+ρ)(1+2ρ)]  −  ρ²

    The mean matters enormously and the two components do not share it. At μ = ½ this reduces
    algebraically to ``2ρ²(1−ρ)/(1+2ρ)``, the gDNA case, but the RNA component sits at κ, where the same ρ
    gives a variance nearly two orders of magnitude larger: a seed at an extreme mean carries far less
    information about ρ than its pair count suggests. Comparing the two components on pair counts alone
    therefore over-credits whichever sits nearer ½, which is why :func:`reconcile_overdispersions` is fed
    the estimators' own precisions and not their null informations.

    Derived, and it introduces no constant. Given the seed's latent sense rate ``p`` (``u = p − ½``),
    ``E[od̂_s | p] = 4u²`` exactly, because ``E[d² − N/4 | p] = N(N−1)u²`` and ``od̂_s`` divides by
    ``N(N−1)/4``. So the between-seed term is ``Var(4u²) = 16(E[u⁴] − (E[u²])²)``, and for the symmetric
    ``Beta`` the strand model already asserts — ``E[u²] = ρ/4``, ``E[u⁴] = 3ρ²/(16(1+2ρ))`` — that collapses
    to ``2ρ²(1−ρ)/(1+2ρ)``. Both forms are verified against Monte Carlo.

    This is why a deep seed is not worth its pair count. Sampling noise vanishes as ``2/(N(N−1))``, but
    this term does not depend on ``N`` at all: past ``c_s·V∞ ≫ ½`` a seed's information about ρ saturates,
    so weighting by pairs (∝ N²) lets a handful of deep seeds decide the answer.
    """
    r = float(overdispersion)
    q = float(mean) * (1.0 - float(mean))
    if q <= 0.0:
        return float("inf")
    return 3.0 * r * r * (2.0 * r + q * (1.0 - 7.0 * r)) / (q * (1.0 + r) * (1.0 + 2.0 * r)) - r * r


def influence_weights(
    pair_scale: np.ndarray,
    overdispersion: float,
    *,
    binom_scale: np.ndarray | float = 0.5,
    mean: float = 0.5,
) -> np.ndarray:
    """``w_s = 1/(½ + c_s·V∞(ρ))`` — minimum-variance pooling weights, the only freedom in the fit.

    The pooled estimate is a weighted mean of the per-seed unbiased estimates ``od̂_s``, carried at weight
    ``w_s·c_s``. Gauss–Markov fixes those weights at ``1/Var(od̂_s | ρ)``: not a preference but the unique
    minimum-variance choice, with any other weighting a strictly worse estimator of the same quantity. By the
    law of total variance ``Var(od̂_s | ρ) = V∞(ρ) + E_p[Var(od̂_s | p)]``, whose second term is the binomial
    sampling noise ``2/(N(N−1))`` at ``p = ½``; substituting and clearing ``c_s`` gives this expression.

    At ρ = 0 the weight is a constant (``w ≡ 2``) that cancels from the ratio, so the pair-count estimator
    is this one with ρ pinned at 0: the weighting does not add an assumption, it removes one.

    The one approximation — the sampling term at ``p = ½`` rather than integrated over ``p`` — cannot bias
    the fit, because the weights depend only on ``n_s`` and ρ and never on a seed's own data, so the ratio
    has expectation ρ under any weight function. A wrong weight costs efficiency, never correctness, and
    this one is measurably more accurate than pair-count weighting when deep seeds are present.
    """
    return 1.0 / (
        np.asarray(binom_scale, dtype=np.float64)
        + np.asarray(pair_scale, dtype=np.float64) * between_seed_variance(overdispersion, mean)
    )


def away_half_moment(
    sense: np.ndarray,
    total: np.ndarray,
    rna_sense_frac: float,
    *,
    overdispersion: float = 0.0,
) -> tuple[float, float]:
    """The away-half pooled moment ``(od_mom, information)`` at an assumed ``overdispersion`` — unclipped.

    ``od_mom = Σ w_s a_s e_s / Σ w_s a_s c_s`` with the influence weights of :func:`influence_weights`.
    The default ``overdispersion = 0`` makes those weights constant, so the default is exactly the
    pair-count moment, this estimator's ρ = 0 special case; :func:`fit_gdna_strand_overdispersion` solves
    for the self-consistent value.

    ``d = (K − N/2)·sign(½ − κ)`` so that RNA of the seed's own gene pulls ``d`` negative; the away half is
    ``d > 0`` at weight 1, a tie ``d = 0`` at weight ½ (the lemma is in the module docstring). Returns
    ``nan`` when nothing carries a pair.

    The degenerate orientation is not a failure, it is the unstranded case. ``κ = ½`` exactly is
    reachable — ``rna_sense_frac`` is the posterior mean ``(n_same + 1)/(n_obs + 2)``, exactly ½ whenever
    ``2·n_same = n_obs``, the modal outcome on an unstranded library — and there ``sign(½ − κ) = 0``. At
    that κ the RNA is unstranded too, so its contamination is symmetric about ½: there is no skew to
    orient by, and the full two-sided moment is used instead. That keeps the one-sided guarantee, since a
    same-mean RNA component contributes only ``n_g(n_g−1)/[N(N−1)] ≤ 1`` of the excess and so biases od
    down, never up. Without this branch every residual collapses to 0, every seed enters as a tie with
    excess ``−N/4``, and the fit returns a hard ``od = 0`` — a claim of perfect Binomiality, the most
    confident strand likelihood assertable — with ``fallback_used`` still False, so nothing signals it.

    The information is the null information of all the seeds (:func:`_null_information` at μ = ½, hence
    the pair count), halved on the away-half branch because only half of each seed's symmetric noise is
    used: with ``Var(e_s)|₀ = n(n−1)/8`` and ``E[a_s] = ½``, ``Var(num) = P/8`` and ``E[den] = P/4``, so
    ``Var(od_mom) = 2/P`` and ``I = P/2`` for the total pair count ``P``. It is not the pair count of the
    seeds that happened to land on the away side: halving that halves twice and overstates the standard
    error by √2. ``1/sqrt(information)`` is the null standard error of ``od_mom``; on the full-moment
    branch no halving applies.
    """
    a, e, c, n = _away_half_parts(sense, total, rna_sense_frac)
    w = influence_weights(c, overdispersion)
    num = float(np.sum(w * a * e))
    den = float(np.sum(w * a * c))
    valid = n > 0.0
    degenerate = float(np.sign(0.5 - float(rna_sense_frac))) == 0.0
    info = _null_information(n[valid], 0.25) * (1.0 if degenerate else 0.5)
    return (num / den if den > 0.0 else float("nan")), info


def seed_participation(
    sense: np.ndarray,
    total: np.ndarray,
    rna_sense_frac: float,
    *,
    overdispersion: float = 0.0,
) -> float:
    """The effective number of seeds behind the pooled moment — its participation ratio.

    ``n_eff = (Σ|x_s|)² / Σ x_s²`` over each seed's signed contribution ``x_s = a_s·(d_s² − n_s/4)`` to the
    numerator, at the fit's own ``overdispersion`` so that it describes the estimate actually made: the
    inverse Simpson index of the contribution masses. ``n_eff = S`` when all ``S`` seeds contribute alike,
    and ``n_eff → 1`` when one seed is the estimate.

    Threshold-free, and that is the point: a "top-k share" would have to pick ``k``, and this picks
    nothing, so the diagnostic introduces no constant.

    It is reported, never acted on. An estimate resting on a handful of seeds is a fact the caller must
    see, not a licence to drop them: trimming the upper tail of the distribution whose mean is the
    estimand would bias the fit down. The concentration is a real-data phenomenon — a simulated panel at a
    true overdispersion of 0 with tens of thousands of seeds cannot produce it.

    Returns ``nan`` when nothing contributes.
    """
    a, e, c, _n = _away_half_parts(sense, total, rna_sense_frac)
    x = influence_weights(c, overdispersion) * a * e
    s1 = float(np.sum(np.abs(x)))
    s2 = float(np.sum(x * x))
    return (s1 * s1 / s2) if s2 > 0.0 else float("nan")


def fit_gdna_strand_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    rna_sense_frac: float,
) -> GdnaStrandModel:
    """Pooled away-half method-of-moments fit of the global gDNA strand overdispersion: the raw moment
    over the RNA-free half of each seed's noise, clipped to the physical support, with no location prior.

    Per seed ``s`` with sense count ``K_s`` of ``N_s`` unspliced fragments, oriented so RNA pulls down::

        d_s      = (K_s − N_s/2)·sign(½ − κ)
        a_s      = 1[d_s > 0] + ½·1[d_s = 0]
        od       = clip( Σ a_s·(d_s² − N_s/4) / Σ a_s·N_s(N_s − 1)/4 ,  0, _MAX_OVERDISPERSION )

    No weight and no purity class: the lemma in the module docstring makes the away half unbiased for
    ρ_g under any RNA content of the seeds, which is what a fit that must hold whatever reference
    transcriptome it is handed requires.

    No shrinkage target either. Where there is no gDNA strand information at all — no seed carrying a
    pair on the away side — the fit returns the ceiling at zero information, the widest strand likelihood
    the model admits, so the channel says nothing; ``calibrate`` then hands this component the RNA fit's
    measured value through :func:`reconcile_overdispersions`, never a conjured one.

    Gate: ``tests/calibration/test_gdna_strand_fit.py``, where the two-sided moment is the perturbation
    each property is watched failing against.
    """
    sense = np.asarray(sense, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)

    def clipped_at(rho: float) -> float:
        m = away_half_moment(sense, total, rna_sense_frac, overdispersion=rho)[0]
        return (
            float(np.clip(m, _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION)) if np.isfinite(m) else m
        )

    raw = away_half_moment(sense, total, rna_sense_frac)[0]  # ρ = 0: the pair-count moment
    fallback = not np.isfinite(raw)
    if fallback:
        # No conjured fallback. With no pair on the away side this component has measured nothing, and
        # the least-committal answer is the ceiling — the widest strand likelihood the model admits, so
        # the channel says nothing rather than something confident. `calibrate` then reconciles it against
        # the RNA fit (:func:`reconcile_overdispersions`), which usually has measured something.
        od = _MAX_OVERDISPERSION
    else:
        # Bisection on a bracket that exists by construction, so no iteration limit is asserted:
        # g(ρ) = clip(moment(ρ)) − ρ has g(0) = clip(raw) ≥ 0 and g(ceiling) ≤ 0, because the clip caps the
        # moment at the ceiling. Halving until the interval is one ULP wide terminates in ~60 steps by
        # floating-point exhaustion. On real libraries g has a single sign change on [0, ceiling] and the
        # root agrees with plain fixed-point iteration.
        lo, hi = _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION
        while hi - lo > np.spacing(hi):
            mid = 0.5 * (lo + hi)
            if clipped_at(mid) - mid > 0.0:
                lo = mid
            else:
                hi = mid
        od = 0.5 * (lo + hi)
    # The precision of the estimate that was made, not the null pair count: for inverse-variance pooling
    # ``Var(pooled) = 1/Σ(1/V_s)`` and ``1/V_s = w_s·c_s``, so the information is the weighted denominator.
    # This is the currency `reconcile_overdispersions` weighs the two components in, and it must be the
    # fitted-ρ one — a null information over-credits whichever component sits nearer mean ½.
    a, _e, c, _n = _away_half_parts(sense, total, rna_sense_frac)
    information = float(np.sum(influence_weights(c, od) * a * c))
    seeded = total > 0.0
    return GdnaStrandModel(
        gdna_strand_overdispersion=od,
        information=float(information),
        n_seed_regions=int(np.sum(seeded)),
        n_seed_fragments=int(np.sum(total[seeded])),
        fallback_used=bool(fallback),
        raw_overdispersion=float(raw),
        clamped_at_ceiling=bool(not fallback and od >= _MAX_OVERDISPERSION),
        effective_seeds=seed_participation(sense, total, rna_sense_frac, overdispersion=od),
    )


def _region_seeds(substrate, region_arrays, region_count_observable):
    """``(sense, total)`` from the count-observable GENIC contained regions — the single-strand intron
    regions (``TS_POS`` / ``TS_NEG``).

    Intergenic (``TS_NONE``) regions are not seeds: with no gene strand there is nothing to orient the
    away half by, so the lemma cannot protect them, and they are where unannotated transcription lives.
    ``TS_AMBIG`` regions have no defined sense and are not seeds either. The count columns are genome
    strand, so orienting them to transcript sense is this function's job.
    """
    ts = np.asarray(region_arrays.strand_class)
    count = np.asarray(substrate.region_contained.count, dtype=np.float64)
    pos, neg = count[:, 0], count[:, 1]
    total = pos + neg
    sense = np.where(ts == TS_NEG, neg, pos)
    seed = (
        np.asarray(region_count_observable, bool)
        & ((ts == TS_POS) | (ts == TS_NEG))
        & (total > 0.0)
    )
    return sense[seed], total[seed]


def fit_gdna_strand_from_substrate(
    substrate,
    region_arrays,
    *,
    region_count_observable: np.ndarray,
    boundary_count_observable: np.ndarray,
    rna_sense_frac: float,
) -> GdnaStrandModel:
    """Fit the global gDNA strand overdispersion from the calibration substrate.

    Pools every count- and strand-observable genic object — intron regions (:func:`_region_seeds`) and
    contiguous boundaries with a defined sense, exon|intron and gene-edge alike
    (:func:`boundary_seeds`) — into the away-half estimator. No weight is computed and no
    class is asserted pure. With no pair on the away side the fit returns the ceiling at zero information,
    which is how ``calibrate`` knows to hand this component the RNA fit's measured value instead
    (:func:`reconcile_overdispersions`).

    The two count-observability masks are the whole input beyond the counts, and they come straight from
    the signature (:func:`density_model.count_observable_masks`). No density, no imputation and no
    opportunity model enter: the away-half moment needs none.
    """
    n_sense, n_total = _region_seeds(substrate, region_arrays, region_count_observable)
    e_sense, e_total = boundary_seeds(substrate, region_arrays, boundary_count_observable)
    return fit_gdna_strand_overdispersion(
        np.concatenate([n_sense, e_sense]),
        np.concatenate([n_total, e_total]),
        rna_sense_frac=rna_sense_frac,
    )


def reconcile_overdispersions(
    rna_moment: float,
    rna_information: float,
    gdna_moment: float,
    gdna_information: float,
) -> tuple[float, float]:
    """``(od_rna, od_gdna)`` — each component's dispersion, with the weaker one shrunk toward the better
    measured one. The reference is a measurement of the same library and the weight is that measurement's
    own information, so no target constant and no asserted prior weight enters.

    Why the other component is the right reference. Physically, both components ride the same library
    preparation — the same local capture bias, the same pair-anchoring noise, the same residual clumping —
    so each is genuine evidence about the other. Structurally, and this is the stronger reason, the
    composition channel reads the ``−½·log var`` difference between the two dispersions: a technical
    artefact shared by both cancels, whereas an asserted difference does not. Under no information the two
    coinciding is therefore exactly the required behaviour, and this rule produces it rather than ruling
    it.

    The rule. The better-informed component keeps its own clipped moment; the weaker one borrows only its
    deficit of information, ``od_w' = (I_w·od_w + (I_s − I_w)·od_s) / I_s``.

    The deficit, not the whole information, and the difference matters at the top end: under
    ``(I_w·od_w + I_s·od_s)/(I_w + I_s)`` two equally well-measured components would still see the weaker
    one dragged to their midpoint while the stronger one did not move — asymmetric, and a claim neither
    measurement supports. Here the borrow weight is ``(I_s − I_w)/I_s``, which is 0 when the two are
    equally informed (neither moves) and 1 when the weak one has measured nothing (it takes the other's
    value outright, subsuming the no-evidence case).

    Not a symmetric blend toward a pooled mean. That is complete pooling and it erases a real
    difference, which the two components do show on real libraries. The weak one moves and the strong one
    does not, which is the limit a hierarchical model reaches when the between-component variance cannot
    be estimated — and with only two groups it cannot.

    A component with no pair at all (a ``nan`` moment, or zero information) takes the other's value
    outright. With neither measured, both return the ceiling: any common value leaves the strand channel
    uninformative, and the ceiling is the one already asserted, so no constant is introduced.
    """
    r_ok = bool(np.isfinite(rna_moment)) and float(rna_information) > 0.0
    g_ok = bool(np.isfinite(gdna_moment)) and float(gdna_information) > 0.0
    clip = lambda v: float(np.clip(v, _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION))  # noqa: E731
    if not r_ok and not g_ok:
        return _MAX_OVERDISPERSION, _MAX_OVERDISPERSION
    if not r_ok:
        return clip(gdna_moment), clip(gdna_moment)
    if not g_ok:
        return clip(rna_moment), clip(rna_moment)
    r, g = clip(rna_moment), clip(gdna_moment)
    i_r, i_g = float(rna_information), float(gdna_information)
    if i_r >= i_g:  # RNA is the better-measured component; gDNA borrows its deficit from it
        return r, clip((i_g * g + (i_r - i_g) * r) / i_r)
    return clip((i_r * r + (i_g - i_r) * g) / i_g), g


def fit_rna_strand_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    rna_sense_frac: float,
) -> RnaStrandModel:
    """Pooled method-of-moments fit of the global RNA strand overdispersion, with prior shrinkage.

    The twin of :func:`fit_gdna_strand_overdispersion`, on pure-RNA spliced seeds: an annotated splice
    junction proves RNA origin, so each seed is all-RNA — the seed mean is ``rna_sense_frac`` (κ) and the
    component fraction is 1, the whole seed being the RNA component. The excess variance of the sense
    split beyond ``Binomial(N, κ)`` identifies the overdispersion. This is the raw pooled moment clamped
    to ``[0, _MAX_OVERDISPERSION]``, exactly as on the gDNA side and with no location prior of its own;
    the shrinkage between the two components happens afterwards, in
    :func:`reconcile_overdispersions`. With no seed carrying a pair the fit returns the ceiling and sets
    ``fallback_used``.

    The information reported here is not the pair count. ``I = Σ n(n−1)/2`` holds only at ``μ = ½``,
    and at the RNA fit's κ the pair count overstates the information by an order of magnitude, so
    :func:`_null_information`'s general form is required (see its docstring).

    Parameters
    ----------
    sense, total : np.ndarray
        Per-seed motif-relative sense count and total qualified count ``N``. The production caller
        is :func:`fit_rna_strand_from_sj_table`, one seed per splice junction.
    rna_sense_frac : float
        Library RNA sense fraction ``κ`` (the spliced-channel ``StrandModel`` mean) — the region mean.
    """
    sense = np.asarray(sense, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)
    kappa = float(rna_sense_frac)
    pq = kappa * (1.0 - kappa)  # the component's own μ(1−μ): the seed is certified pure RNA
    valid = total > 0.0
    excess = (sense - total * kappa) ** 2 - total * pq
    scale = np.maximum(total * (total - 1.0), 0.0) * pq
    num = float(np.sum(excess[valid]))
    den = float(np.sum(scale[valid]))
    raw = num / den if den > 0.0 else float("nan")
    fallback = not np.isfinite(raw)
    od = (
        _MAX_OVERDISPERSION
        if fallback
        else float(np.clip(raw, _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION))
    )
    return RnaStrandModel(
        rna_strand_overdispersion=od,
        raw_overdispersion=float(raw),
        information=_null_information(total[valid], pq),
        n_seed_regions=int(np.sum(valid)),
        n_seed_fragments=int(np.sum(total[valid])),
        fallback_used=bool(fallback),
    )


def fit_rna_strand_from_sj_table(
    sj_table,
    *,
    rna_sense_frac: float,
) -> RnaStrandModel:
    """Fit the global RNA strand overdispersion from the per-sj SJ strand table.

    One population, one source of truth: ``rna_sense_frac`` (κ, the mean) is the marginal of this same
    table, so both halves of the Beta-Binomial are estimated from the same strand-qualified population
    (``SPLICE_SPLICED_ANNOT``, unique-mapper, unambiguous exon and SJ strand, non-chimeric), one
    observation per fragment. Each sj ``j`` is one seed ``(sense_j, n_j)``, and the dispersion is the
    spread of those splits across sj at mean κ.

    ⛔ The accumulator's boundary spliced channels are not an alternative source. Their population also
    pools ``SPLICED_UNANNOT`` and — the damaging one — ``SPLICED_IMPLICIT``: an implicit splice's mate gap
    spans an annotated intron, so it lands on exactly the annotated boundaries where the genuine sj live,
    and its motif was never sequenced (the annotation supplies the strand), so its sense bit is arbitrary
    and sits near ½. Fitting a dispersion about a κ near 0 from seeds sitting at ½ returns a raw moment
    far above 1, impossible for an intra-class correlation, and clips at the ceiling on every real
    library: a mean misfit reported as dispersion, carrying no information at all.

    The accumulator's spliced mass channel is deliberately untouched — it legitimately wants implicit and
    novel RNA for the SPLICE OUT, the SPLICE IN and the mature-RNA floor. This belongs in the strand
    model, not in the deposit predicate.

    A sj with one fragment contributes no pair to correlate and therefore no information, which
    :func:`_null_information` handles with no depth threshold.
    """
    n_sense = np.asarray(sj_table.n_sense, dtype=np.float64)
    n_total = n_sense + np.asarray(sj_table.n_antisense, dtype=np.float64)
    return fit_rna_strand_overdispersion(n_sense, n_total, rna_sense_frac)


def boundary_strand_orientation(
    ts_lo: np.ndarray, ts_hi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """``(orient_neg, strand_observable)`` per contiguous boundary, from its two flanks' strand classes.

    A boundary is strand-observable iff its two flanks define a single consistent transcript sense:
    ``{POS,POS}`` / ``{NEG,NEG}``, or a gene-boundary ``{POS,NONE}`` / ``{NEG,NONE}``. An intergenic
    (``TS_NONE``) flank is a strand wildcard — it carries no transcript, so it is compatible with
    either strand and the boundary is oriented by its gene flank. Only a genuine conflict — opposite
    strands ``{POS,NEG}``, or a ``TS_AMBIG`` flank (overlapping ± transcripts) — leaves the sense
    undefined, and such a boundary cannot seed the fit at all.

    The rule is symmetric in the two flanks, so a boundary's orientation does not depend on which side
    it is read from.

    ``orient_neg`` selects the NEG genome-strand column as "sense"; where it is ``False`` (and the
    boundary is observable) the POS column is sense.

    There is no explicit AMBIG guard and there must not be one: ``TS_AMBIG`` is a fourth distinct value,
    so ``(ts == TS_POS) | (ts == TS_NONE)`` is already ``False`` on it. An ``~either_ambig`` term would be
    dead code reading as load-bearing, which is worse than absent because it invites the belief that the
    rule lives there. ``test_an_AMBIG_flank_cannot_seed`` pins the rule instead.
    """
    lo_pos_or_none = (ts_lo == TS_POS) | (ts_lo == TS_NONE)
    hi_pos_or_none = (ts_hi == TS_POS) | (ts_hi == TS_NONE)
    lo_neg_or_none = (ts_lo == TS_NEG) | (ts_lo == TS_NONE)
    hi_neg_or_none = (ts_hi == TS_NEG) | (ts_hi == TS_NONE)
    cons_pos = lo_pos_or_none & hi_pos_or_none & ((ts_lo == TS_POS) | (ts_hi == TS_POS))
    cons_neg = lo_neg_or_none & hi_neg_or_none & ((ts_lo == TS_NEG) | (ts_hi == TS_NEG))
    return cons_neg, cons_pos | cons_neg


def boundary_seeds(
    substrate, region_arrays, boundary_count_observable
) -> tuple[np.ndarray, np.ndarray]:
    """``(sense, total)`` — one seed per count- and strand-observable contiguous boundary.

    The exon|intron and exon|intergenic boundary seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the intron-region seeds (needed under hybrid capture, which
    depletes off-target intronic gDNA while a boundary abutting a probed exon stays counted). Sense is
    the boundary's gene strand (:func:`boundary_strand_orientation`); a boundary with no defined sense
    is not a seed. Gate: ``tests/calibration/test_gdna_strand_fit.py``.
    """
    ts = np.asarray(region_arrays.strand_class)
    lo, hi = boundary_region_indices(np.asarray(region_arrays.ref_id))
    count = np.asarray(substrate.boundary_unspliced.count, dtype=np.float64)
    pos, neg = count[:, 0], count[:, 1]
    total = pos + neg

    orient_neg, strand_observable = boundary_strand_orientation(ts[lo], ts[hi])
    count_observable = np.asarray(boundary_count_observable, dtype=bool)
    seed = count_observable & strand_observable & (total > 0.0)
    return np.where(orient_neg, neg, pos)[seed], total[seed]


__all__ = [
    "GdnaStrandModel",
    "boundary_seeds",
    "boundary_strand_orientation",
    "RnaStrandModel",
    "away_half_moment",
    "reconcile_overdispersions",
    "between_seed_variance",
    "influence_weights",
    "seed_participation",
    "fit_gdna_strand_overdispersion",
    "fit_gdna_strand_from_substrate",
    "fit_rna_strand_overdispersion",
    "fit_rna_strand_from_sj_table",
    "overdispersion_for_beta",
]
