"""Strand Beta-Binomial overdispersion fits for both components (gDNA and RNA).

Both gDNA and RNA strand splits are modelled as **Beta-Binomial**:
``K⁺ | N ~ BetaBinom(N, mean, overdispersion)`` with the intra-class correlation in ``[0, 1)``
(``0`` ⇒ Binomial). The **mean** differs by component — gDNA is unstranded (biologically fixed at
½), RNA is skewed (``rna_sense_frac`` = κ from the spliced channel) — but the **overdispersion**
machinery is shared: a per-region/per-boundary sense split is more spread out than Binomial
(local capture bias, PCR, pair-anchoring noise), and that excess variance is the overdispersion.

**Symmetric by design.** Earlier the RNA side was forced Binomial while gDNA was Beta-Binomial;
that asymmetry made the strand likelihood *spuriously informative on unstranded data* (the
``−½·log var`` term prefers the lower-variance component, pulling balanced regions toward RNA). Both
components now carry a fitted overdispersion with the **same default prior**, so under sparse data
the two collapse to the same distribution and an unstranded region is uninformative — as it must be.
The gDNA mean ½ and RNA mean κ are unchanged; only RNA gains the overdispersion term.

⭐ **The RNA overdispersion is fit from the PER-JUNCTION SJ strand table — the same population κ is the
marginal of** (:func:`fit_rna_strand_from_sj_table`). It used to be scavenged from the accumulator's
boundary spliced channels, which pool unannotated and implicit splices; an implicit splice has no
sequenced motif, so its sense bit is arbitrary, and fitting a dispersion about κ ≈ 0.002 from seeds
sitting at ½ produced a raw ``od_mom`` of 10.7–79.9 — impossible for an intra-class correlation — clipped
to the ceiling on 4/4 real libraries. Two halves of one Beta-Binomial, two different populations.


**Breaking the circularity** (gDNA only). Fitting the overdispersion needs to know which fragments are gDNA,
which is what the deconvolution determines — circular. We break it with the count⊥strand
conditional independence the engine already relies on: the **count module** supplies a per-region
gDNA fraction ``gdna_weight`` (``count_gdna_frac``, a raw count/density ratio) that is *independent*
of the strand overdispersion (it uses no strand information at all). Given those weights and the RNA
sense rate, the overdispersion is identified from the **excess variance of the sense split beyond
Binomial**, attributable to the gDNA fragments.

**Estimator — pooled method of moments.** For each seed region ``s`` (a count-observable region or
boundary side — intergenic, intronic, exon–intron / exon–intergenic boundary) with ``sense_s`` of
``n_s`` gDNA-eligible unspliced fragments and count-derived gDNA weight ``w_s``:

    mean_s        = ½·w_s + rna_sense_frac·(1 − w_s)        # mixture sense rate
    excess_var_s  = (sense_s − n_s·mean_s)²  −  n_s·mean_s·(1 − mean_s)   # observed − Binomial
    gdna_var_s    = (w_s·n_s)·(w_s·n_s − 1)·¼                # BetaBinom excess-variance scale

    od_mom = Σ_s excess_var_s / Σ_s gdna_var_s        # pooled point estimate

The point estimate is then **shrunk toward a prior overdispersion** ``od₀`` = ``Beta(14,14)`` ⇒ 0.0345,
**weighted by INFORMATION** (:func:`_null_information`), and clamped to ``[0, _MAX_OVERDISPERSION]`` (the
``Beta(2,2)`` ceiling, od = 0.2, the most overdispersion allowed).

⭐ **The shrinkage currency is INFORMATION, not seed count, and that is a units fix, not a tuning choice.**
Overdispersion is a correlation *between* fragments in a region, so a seed with ONE fragment carries none of
it — it has nothing to be correlated with. Weighting by ``n_seed_regions`` let 160 k singleton seeds outvote
the prior by four orders of magnitude on a library that has almost no pairs. With ``I = 1/Var(od_mom)|₀``
the expression becomes an exact Normal–Normal posterior mean and the prior's weight ``W`` is **DERIVED**
from the two asserted constants rather than asserted itself (:func:`_prior_information`).

⚠ **Measured: this changes nothing on real data** (identical to 4 d.p.), because real libraries carry
0.7 M–101 M information units against a prior worth ~909. It binds on the zero-seed fallback and on unit
fixtures — which is exactly where it should. The saturation seen on real cfRNA is **bias from a
contaminated seed channel**, not a shrinkage problem.

The MoM is closed-form, ``O(n_seed_regions)``, and uses the **same variance decomposition the deconv
applies** (ψ's ``simplex_logodds._mixture_strand_loglik``; its two-component reference is
:mod:`strand_likelihood`), so fit and application are consistent. The two constants live in this
module, next to the estimator they parameterise.

Each component pairs a pure estimator with a thin seed-extraction wrapper, so the estimator itself
is trivially testable: :func:`fit_gdna_strand_overdispersion` / :func:`fit_gdna_strand_from_substrate`
for gDNA (the wrapper needs the region/boundary geometry), and
:func:`fit_rna_strand_overdispersion` / :func:`fit_rna_strand_from_sj_table` for RNA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .strand_deconv import boundary_seeds
from .signature import TS_AMBIG, TS_NEG

#: Hard floor on overdispersion — the Binomial limit (no negative overdispersion).
_OVERDISPERSION_FLOOR: float = 0.0
#: Symmetric ``Beta(a, a)`` shape at the **most overdispersion we allow** (the per-region sense
#: rate is at most this spread out). ``a = 2`` ⇒ ``od = 1/(2·2+1) = 0.2`` — the overdispersion
#: ceiling. The fitted/shrunk overdispersion is clamped to ``[0, _MAX_OVERDISPERSION]``.
_CEIL_ALPHA_BETA: float = 2.0


def overdispersion_for_beta(alpha_beta: float) -> float:
    """Beta-Binomial intra-class correlation (overdispersion) for a symmetric ``Beta(a, a)``.

    ``od = 1 / (2·a + 1)``: ``a = 1`` (uniform) → 1/3; ``a = 2`` → 1/5 = 0.2 (max overdispersion);
    ``a = 3`` → 1/7 ≈ 0.143; ``a → ∞`` → 0 (Binomial). Inverse of ``beta_concentration``.
    """
    return 1.0 / (2.0 * float(alpha_beta) + 1.0)


#: Overdispersion ceiling (most overdispersion allowed) = ``Beta(2, 2)`` ⇒ ``od = 0.2``.
_MAX_OVERDISPERSION: float = overdispersion_for_beta(_CEIL_ALPHA_BETA)

#: **K2 — the shrinkage TARGET**, the symmetric ``Beta(a, a)`` the fit falls back to when the data carry
#: no information. ``a = 14`` ⇒ ``od₀ = 1/29 = 0.0345``. ⭐ **No longer arbitrary**: two assumption-light
#: measurements now bracket the truth — gDNA from the exact ``od(n=2) = 2·P(both same strand) − 1`` readout
#: (no Beta assumption, no estimated mean, since μ = ½ is asserted biology) gives a plateau of
#: **0.007–0.028**; RNA from junctions deep enough to see the minority strand gives **0.0011–0.0158**, and
#: the synthetic suite (true od = 0 by construction) gives 0.0008–0.0017. ``od₀`` sits **1.2–30× ABOVE** the
#: top of every honest measurement, i.e. at the conservative end of measured reality, which is what a
#: fallback should be. It remains ASSERTED — the measurements bracket it, they do not derive it.
#: **Shared by the gDNA and RNA fits, and that is required, not a convenience**: at κ = ½ with no data the
#: two components must coincide, or ψ's ``−½·log var`` term hands an unstranded region a spurious gDNA/RNA
# preference (see:mod:`.strand_likelihood`). Full record:.
_PRIOR_ALPHA_BETA: float = 14.0
_PRIOR_OVERDISPERSION: float = overdispersion_for_beta(_PRIOR_ALPHA_BETA)


def _prior_information() -> float:
    """**DERIVED** prior precision ``W``, in the data's own information units — not asserted.

    The shrinkage ``(I·od_mom + W·od₀)/(I + W)`` is a Normal–Normal posterior mean once ``I`` is the
    estimator's null information (:func:`_null_information`), because ``Var(od_mom)|₀ = 1/I`` exactly.
    So ``W = 1/τ²`` with ``τ²`` the PRIOR variance of ``od`` — and the prior is already fully pinned by the
    two asserted constants: it lives on ``[0, _MAX_OVERDISPERSION]`` (K1) with mean ``od₀`` (K2). The
    least-committal distribution under exactly those two constraints is the maximum-entropy one (a
    truncated exponential), whose variance closes the derivation with **no third constant**.

    ⚠ **Name the approximation honestly:** reducing a bounded prior to two moments is *Bühlmann linear
    credibility*, not exact Bayes (exact Bayes under a bounded prior at these information levels is simply
    ``clip(od_mom, 0, od_max)``, i.e. no shrinkage at all). The linear rule is adopted because it degrades
    gracefully in the low-information corner, where the clip alone would return a hard ``0`` — a claim of
    perfect Binomiality, the *most* confident strand likelihood we could possibly assert.

    ⚠ **And it is measurably INERT on real data**: the fitted overdispersion is identical to four decimals
    across ``W ∈ {30, 100, 300, 588, 2658}`` and with the shrinkage deleted entirely, because a real library
    carries 0.7 M–101 M information units against this ~909. It binds only on the zero-seed fallback and on
    unit fixtures. That inertness is why the *shape* assertion above costs nothing.
    """
    b, m = _MAX_OVERDISPERSION, _PRIOR_OVERDISPERSION

    # max-entropy density on [0, b] with mean m is ∝ exp(−λx); solve mean(λ) = m by bisection.
    def _mean(lam: float) -> float:
        z = lam * b
        if abs(z) < 1e-9:  # λ → 0 is the uniform limit
            return b / 2.0
        if z > 700.0:  # e^z overflows; the b/(e^z − 1) term is then exactly negligible
            return 1.0 / lam
        return 1.0 / lam - b / np.expm1(z)

    lo, hi = 1e-6, 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _mean(mid) > m:
            lo = mid
        else:
            hi = mid
    lam = 0.5 * (lo + hi)
    z = lam * b
    var = 1.0 / lam**2 - b * b * np.exp(z) / np.expm1(z) ** 2 if z < 700 else 1.0 / lam**2
    return float(1.0 / max(var, 1e-12))


#: ``W`` — the prior's weight in information units. ≈ 909; computed, never asserted.
_PRIOR_INFORMATION: float = _prior_information()


@dataclass(frozen=True, slots=True)
class GdnaStrandModel:
    """Fitted gDNA strand model: the global Beta-Binomial overdispersion + fit provenance."""

    gdna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
    n_seed_regions: int  # seed regions that carried gDNA-eligible fragments
    n_seed_fragments: int  # total gDNA-eligible fragments across seed regions
    fallback_used: bool  # True ⇒ no gDNA strand signal ⇒ returned the prior overdispersion

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
    :class:`strand_balance.StrandBalance`; this carries only the between-JUNCTION overdispersion."""

    rna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
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
    """The estimator's **null information** ``I = 1/Var(od_mom)|₀`` — exact, not asymptotic.

    ``od_mom = Σe_s / Σc_s`` with ``c_s = n_c(n_c−1)·pq`` deterministic, so ``Var(od_mom) = Var(Σe_s)/(Σc_s)²``.
    Under the null (``od = 0``) the sense count is Binomial, and with ``X = K − n·μ``,
    ``E[X²] = n·pq`` and ``E[X⁴] = 3(n·pq)² + n·pq(1 − 6pq)``, so ``Var(e_s) = 2(n·pq)² + n·pq(1 − 6pq)``:

        I = (Σ n(n−1)·pq)² / Σ [ 2n²pq² + n·pq − 6n·pq² ]

    ⭐ **At μ = ½ this reduces EXACTLY to the pair count** ``Σ n(n−1)/2`` — the intuition being that
    overdispersion is a correlation *between* fragments, so the unit of evidence is a PAIR (a seed of one
    fragment contributes nothing; a seed of ``n`` contributes ``n(n−1)/2``). Verified algebraically and by
    Monte Carlo (ratio 0.994–1.003 over μ ∈ [0.002, 0.5] and mixed seed-size populations).

    ⚠ **Away from μ = ½ the pair count is WRONG and must not be substituted**: at the RNA fit's κ the
    measured ``I/pairs`` is **0.05–0.14**, i.e. the raw pair count overstates the information 7–20×. That is
    why this is a general form rather than a pair count, and why the two components share one code path.
    """
    n = np.asarray(comp_frags, dtype=np.float64)
    pq = np.asarray(comp_var, dtype=np.float64)
    num = float(np.sum(np.maximum(n * (n - 1.0), 0.0) * pq)) ** 2
    den = float(np.sum(2.0 * n * n * pq * pq + n * pq - 6.0 * n * pq * pq))
    return num / den if den > 0.0 else 0.0


def _fit_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    region_mean: np.ndarray,
    component_frac: np.ndarray,
    component_mean: np.ndarray,
    *,
    prior_overdispersion: float = _PRIOR_OVERDISPERSION,
    prior_weight: float = _PRIOR_INFORMATION,
) -> tuple[float, int, int, bool]:
    """Shared pooled-MoM + prior-shrinkage core for one component's strand overdispersion.

    Component-agnostic, with three per-region inputs:

    * ``region_mean`` — the region's mixture sense rate, used to subtract the Binomial variance:
      ½·w + κ·(1−w) for a gDNA seed (a gDNA/RNA mix), κ for a pure-RNA spliced seed.
    * ``component_frac`` — the fraction of the region's fragments in the component being fit (the
      gDNA weight ``w`` for gDNA, ``1`` for pure-RNA spliced); sets ``n_c = component_frac·N``.
    * ``component_mean`` — the component's *own* sense mean ``μ_c`` (½ for gDNA, κ for RNA); the
      BetaBinom excess variance of ``n_c`` correlated fragments scales as
      ``n_c·(n_c − 1)·μ_c·(1 − μ_c)`` — **not** ¼ unless ``μ_c = ½``.

    Returns ``(overdispersion, n_seed_regions, n_seed_fragments, fallback_used)``; clamped to
    ``[0, _MAX_OVERDISPERSION]``; fallback (pooled denominator ≤ 0) ⇒ ``prior_overdispersion``.
    """
    sense = np.asarray(sense, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)
    region_mean = np.asarray(region_mean, dtype=np.float64)
    component_frac = np.asarray(component_frac, dtype=np.float64)
    component_mean = np.asarray(component_mean, dtype=np.float64)

    valid = total > 0.0
    binom_var = total * region_mean * (1.0 - region_mean)
    excess_var = (sense - total * region_mean) ** 2 - binom_var
    # ⚠ INTEGER component counts. ``component_frac`` is a float that is 1.0 only to within rounding
    # (measured min − 1 = −2.2e−16), so ``w·N`` on a 2-fragment seed is 1.9999999999999996 and every
    # ``n_c ≥ 2`` test silently under-counts. Round: a fragment count is an integer.
    comp_frags = np.rint(component_frac * total)
    comp_var = component_mean * (1.0 - component_mean)  # μ_c(1−μ_c): ¼ for gDNA, κ(1−κ) for RNA
    # BetaBinom excess-variance scale n_c(n_c − 1)·μ_c(1−μ_c) (clipped ≥ 0 for tiny component mass).
    var_scale = np.maximum(comp_frags * (comp_frags - 1.0), 0.0) * comp_var

    num = float(np.sum(excess_var[valid]))
    denom = float(np.sum(var_scale[valid]))
    n_seed_regions = int(np.sum(valid & (comp_frags > 0.0)))
    n_seed_fragments = int(np.sum(total[valid & (comp_frags > 0.0)]))
    prior_overdispersion = float(prior_overdispersion)
    prior_weight = max(float(prior_weight), 0.0)

    fallback = denom <= 0.0 or not np.isfinite(num)
    if fallback:
        # No component strand information in the seeds → fall back entirely to the prior.
        od = prior_overdispersion
    else:
        od_mom = num / denom
        # ── SHRINKAGE IN THE DATA'S OWN CURRENCY (the units fix) ────────────────────────────────────
        # This used to weight the data by ``n_seed_regions``, which is the wrong measure of evidence
        # about a SECOND moment: overdispersion is a correlation BETWEEN fragments, so a seed with one
        # fragment carries none of it — it has nothing to be correlated with. Counting seeds let
        # 160k singleton seeds outvote the prior by four orders of magnitude on a library that has
        # almost no pairs. ``_null_information`` is the honest measure and makes the expression an
        # exact Normal–Normal posterior mean (see :func:`_prior_information`).
        info = _null_information(comp_frags[valid], np.broadcast_to(comp_var, total.shape)[valid])
        total_weight = info + prior_weight
        od = (
            (info * od_mom + prior_weight * prior_overdispersion) / total_weight
            if total_weight > 0.0
            else od_mom
        )
    od = float(np.clip(od, _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION))
    return od, n_seed_regions, n_seed_fragments, fallback


def fit_gdna_strand_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    gdna_weight: np.ndarray,
    rna_sense_frac: float,
    *,
    prior_overdispersion: float = _PRIOR_OVERDISPERSION,
    prior_weight: float = _PRIOR_INFORMATION,
) -> GdnaStrandModel:
    """Pooled method-of-moments fit of the global gDNA strand overdispersion, with prior shrinkage.

    The pooled MoM point estimate is shrunk toward ``prior_overdispersion`` weighted by
    INFORMATION (:func:`_null_information`), not by seed-region count:
    ``od = (I·od_mom + prior_weight·prior_overdispersion) / (I + prior_weight)``.
    Seed sets carrying little information lean on the prior; abundant ones on the fit. This
    replaces the earlier hard min-region / significance gates (no thresholds — graceful with data).
    The result is clamped to ``[0, _MAX_OVERDISPERSION]`` (the ``Beta(2, 2)`` ceiling).

    Parameters
    ----------
    sense, total : np.ndarray
        Per-seed-region sense count ``K⁺`` and total gDNA-eligible unspliced count ``N``.
    gdna_weight : np.ndarray
        Per-seed-region count-clue gDNA fraction ``∈ [0, 1]`` (independent of the overdispersion).
    rna_sense_frac : float
        Library RNA sense fraction ``κ`` (the spliced-channel ``StrandModel`` mean).
    prior_overdispersion : float
        Prior overdispersion to shrink toward (the ``Beta(a, a)`` "floor"; see
        :func:`overdispersion_for_beta`). ``0`` with ``prior_weight = 0`` ⇒ pure MoM.
    prior_weight : float
        Prior strength in the estimator's own INFORMATION units (:func:`_prior_information`
        derives the production value). ``0`` ⇒ no shrinkage.

    Returns
    -------
    GdnaStrandModel
        ``gdna_strand_overdispersion`` in ``[0, _MAX_OVERDISPERSION]``; ``fallback_used`` when the
        seeds carry no gDNA strand signal (pooled denominator ≤ 0) ⇒ the prior is returned.
    """
    weight = np.clip(np.asarray(gdna_weight, dtype=np.float64), 0.0, 1.0)
    # gDNA component: the region's sense mean is the mixture ½·w + κ·(1−w); the component whose
    # shared rate inflates the variance is the gDNA fraction w, with its own mean ½ (⇒ scale ¼).
    region_mean = 0.5 * weight + float(rna_sense_frac) * (1.0 - weight)
    component_mean = np.full(region_mean.shape, 0.5, dtype=np.float64)
    od, n_seed_regions, n_seed_fragments, fallback = _fit_overdispersion(
        sense,
        total,
        region_mean,
        weight,
        component_mean,
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )
    return GdnaStrandModel(
        gdna_strand_overdispersion=od,
        n_seed_regions=n_seed_regions,
        n_seed_fragments=n_seed_fragments,
        fallback_used=fallback,
    )


def _region_seeds(substrate, region_arrays, region_density):
    """``(sense, total, gdna_weight)`` from the count-observable CONTAINED regions.

    Intergenic (``TS_NONE``) and intron-only (``TS_POS``/``TS_NEG``) regions — i.e.
    ``region_density.region_count_observable`` — excluding ``TS_AMBIG`` (both strands, no defined
    sense). The weight is the count-clue gDNA fraction ``region_density.count_gdna_frac`` (=
    ``clip(density·eff_gdna / count)``, density cleaned by the strand *mean* ½, not the dispersion).
    It reads the explicit count-prior MEAN (``count_gdna_frac``) directly — decoupled from the
    count-prior concentration, which carries the overdispersion-honest precision.

    ⚠ The columns are GENOME strand, so orienting to transcript sense is this function's job
    (arbitrary but consistent for ``TS_NONE`` — an intergenic region has no transcript to be sense to,
    and gDNA's strand mean is ½ either way).
    """
    ts = np.asarray(region_arrays.strand_class)
    count = np.asarray(substrate.region_contained.count, dtype=np.float64)
    pos, neg = count[:, 0], count[:, 1]
    total = pos + neg
    sense = np.where(ts == TS_NEG, neg, pos)
    weight = np.clip(np.asarray(region_density.count_gdna_frac, dtype=np.float64), 0.0, 1.0)
    seed = np.asarray(region_density.region_count_observable) & (ts != TS_AMBIG)
    return sense[seed], total[seed], weight[seed]


def fit_gdna_strand_from_substrate(
    substrate,
    region_arrays,
    region_density,
    *,
    rna_sense_frac: float,
    prior_overdispersion: float = _PRIOR_OVERDISPERSION,
    prior_weight: float = _PRIOR_INFORMATION,
) -> GdnaStrandModel:
    """Fit the global gDNA strand overdispersion from the calibration substrate.

    Pools two kinds of count-observable seed (the same seeds the density estimator trusts):

    * **contained regions** — intergenic + intron-only (:func:`_region_seeds`);
    * **contiguous boundaries** — exon–intron / exon–intergenic lines
      (:func:`strand_deconv.boundary_seeds`), needed under hybrid capture, which depletes off-target
      intergenic/intronic gDNA.

    Both contribute ``(sense, total, gdna_weight)`` on the same footing, and the pooled estimator fits
    one global overdispersion, shrunk toward ``prior_overdispersion`` (strength ``prior_weight``).

    ⚠ **The boundary arm contributes one seed per line, not two per boundary** (S5.f). The predecessor
    counted each physical crossing twice — once from each face — which inflated the pooled sample size
    2× and paired every observation with a perfectly correlated twin. A dispersion estimator reading
    its own duplication is the failure mode this removes; see :mod:`strand_deconv`.
    """
    n_sense, n_total, n_weight = _region_seeds(substrate, region_arrays, region_density)
    e_sense, e_total, e_weight = boundary_seeds(substrate, region_arrays, region_density)
    sense = np.concatenate([n_sense, e_sense])
    total = np.concatenate([n_total, e_total])
    weight = np.concatenate([n_weight, e_weight])
    return fit_gdna_strand_overdispersion(
        sense,
        total,
        weight,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )


def fit_rna_strand_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    rna_sense_frac: float,
    *,
    prior_overdispersion: float = _PRIOR_OVERDISPERSION,
    prior_weight: float = _PRIOR_INFORMATION,
) -> RnaStrandModel:
    """Pooled method-of-moments fit of the global RNA strand overdispersion, with prior shrinkage.

    The twin of :func:`fit_gdna_strand_overdispersion`, on **pure-RNA spliced** seeds: an annotated
    splice junction proves RNA origin, so each seed is all-RNA — ``region_mean = rna_sense_frac`` (κ)
    and the component fraction is ``1`` (the whole seed is the RNA component). The excess variance of
    the sense split beyond ``Binomial(N, κ)`` identifies the overdispersion. Shrinks toward
    ``prior_overdispersion`` by INFORMATION and clamps to ``[0, _MAX_OVERDISPERSION]`` exactly
    like the gDNA fit, so under sparse data the two collapse to the same prior.

    ⚠ The information here is NOT the pair count. ``I = Σ n(n−1)/2`` holds only at ``μ = ½``; at the
    RNA fit's κ the measured ``I/pairs`` is 0.05–0.14, so :func:`_null_information`'s general form is
    required (see its docstring).

    Parameters
    ----------
    sense, total : np.ndarray
        Per-seed motif-relative sense count and total qualified count ``N``. The production caller
        is :func:`fit_rna_strand_from_sj_table`, one seed per splice junction.
    rna_sense_frac : float
        Library RNA sense fraction ``κ`` (the spliced-channel ``StrandModel`` mean) — the region mean.
    """
    total = np.asarray(total, dtype=np.float64)
    # Pure-RNA spliced region: mixture mean AND component mean are both κ; component fraction is 1.
    region_mean = np.full(total.shape, float(rna_sense_frac), dtype=np.float64)
    component_frac = np.ones(total.shape, dtype=np.float64)  # pure RNA
    od, n_seed_regions, n_seed_fragments, fallback = _fit_overdispersion(
        sense,
        total,
        region_mean,
        component_frac,
        region_mean,  # component_mean = κ ⇒ scale κ(1−κ)
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )
    return RnaStrandModel(
        rna_strand_overdispersion=od,
        n_seed_regions=n_seed_regions,
        n_seed_fragments=n_seed_fragments,
        fallback_used=fallback,
    )


def fit_rna_strand_from_sj_table(
    sj_table,
    *,
    rna_sense_frac: float,
    prior_overdispersion: float = _PRIOR_OVERDISPERSION,
    prior_weight: float = _PRIOR_INFORMATION,
) -> RnaStrandModel:
    """Fit the global RNA strand overdispersion from the **per-junction SJ strand table**.

    ⭐ **One population, one source of truth.** ``rna_sense_frac`` (κ, the mean) is the marginal of
    this same table — both halves of the Beta-Binomial are now estimated from the same
    strand-qualified population (``SPLICE_SPLICED_ANNOT``, unique-mapper, unambiguous exon and SJ
    strand, non-chimeric), one observation per fragment. Each junction ``j`` is one seed:
    ``(sense_j, n_j)``, and the dispersion is the spread of those splits ACROSS junctions at mean κ.

    ⛔ **What this replaces, and why.** ``od_r`` used to be scavenged from the accumulator's
    boundary spliced channels, whose population also pools ``SPLICED_UNANNOT`` and — the damaging
    one — ``SPLICED_IMPLICIT``. An implicit splice's mate gap spans an annotated intron, so it lands
    on exactly the annotated boundaries where the genuine junctions live, and its motif was never
    sequenced (the annotation supplies the strand), so its sense bit is arbitrary: measured ~0.49 at
    depth ≥ 10 on real cfRNA. Fitting a dispersion about κ ≈ 0.002 from seeds sitting at ½ gave a
    raw pooled ``od_mom`` of **10.7–79.9** — impossible for an intra-class correlation, which is
    bounded by 1 — clipped to the 0.2 ceiling on 4/4 real libraries. That was a mean misfit being
    reported as dispersion, and it left ``od_r`` carrying no information about any real library.


    ⚠ The accumulator's spliced MASS channel is deliberately untouched: it legitimately wants
    implicit and novel RNA for the peel, the graft and the mature-RNA floor. This fix belongs in the
    strand model, not in the deposit predicate.

    A junction with one fragment contributes no pair to correlate and so no information; the
    information-currency shrinkage (:func:`_null_information`) handles that with no depth threshold.
    """
    n_sense = np.asarray(sj_table.n_sense, dtype=np.float64)
    n_total = n_sense + np.asarray(sj_table.n_antisense, dtype=np.float64)
    return fit_rna_strand_overdispersion(
        n_sense,
        n_total,
        rna_sense_frac,
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )


__all__ = [
    "GdnaStrandModel",
    "RnaStrandModel",
    "fit_gdna_strand_overdispersion",
    "fit_gdna_strand_from_substrate",
    "fit_rna_strand_overdispersion",
    "fit_rna_strand_from_sj_table",
    "overdispersion_for_beta",
]
