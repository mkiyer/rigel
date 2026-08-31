"""rigel.calibration.fl — gDNA / RNA fragment-length distributions (PR 4c).

Produces the library-wide FL distributions the calibrator's **effective lengths**
need (boundary ``fl_mean``, region ``E_f[max(0, L−ℓ)]``): the **gDNA FL** (from
gDNA-dominated regions/boundaries) and the **RNA FL** (spliced fragments), each
**smoothly empirical-Bayes-shrunk** toward the global FL.

This is **not** a per-fragment FL likelihood — that channel is deliberately
excluded from calibration. BOTH FLs drive per-region effective lengths in the sweep
(``region_geometry.build_region_geometry``): gDNA eff-lengths use the gDNA FL, RNA (nascent
unspliced + spliced) eff-lengths use the RNA FL.

⛔⛔ **NO POOL IS PURE, AND THIS MODULE USED TO ASSERT ONE WAS.** The claim that stood here — *"every pool
is PURE BY CONSTRUCTION, and that purity is what removes the circularity"* — is false and was measured
false: against an origin-split oracle the intronic pool runs **95 % nascent RNA** and the intergenic pool
**53 % mature**, because "intergenic" is whatever the annotation leaves over and nascent RNA sits inside
introns by definition (`TRAPS: purity-is-a-property-of-the-annotation`). The resulting bias is
``RNA_share x (len_RNA - len_gDNA)``, which is why it went unseen: the panels give both components EQUAL
fragment lengths on purpose, so the second factor is zero there and a 95 %-contaminated pool reads under a
bp of error.

⭐⭐ **WHAT REPLACES IT: a two-pool CONTRAST, which needs no pure pool and no template.** The two CONTAINED
pools share one opportunity geometry — nascent RNA in an intron and unannotated transcription in
intergenic space are genomically contiguous exactly as gDNA is — so de-tilted, each is the SAME two shapes
at a DIFFERENT mixing weight, and the contaminant cancels::

    f_0 = a_0*g + (1-a_0)*r        f_1 = a_1*g + (1-a_1)*r
    =>  g = [ (1-a_1)*f_0 - (1-a_0)*f_1 ] / (a_0 - a_1)

⭐ The divisor is the SEPARATION of the two purities, not a purity, so it does not blow up as a pool gets
dirty; and with ``f_0 = f_1`` it returns them unchanged, so a pool pair with nothing to say changes
nothing. The weights come from :mod:`rigel.calibration.gdna_density` — ``a_p = rho_g*E_p/n_p`` with
``rho_g`` read off the low side of the per-object density, where a contaminant that only ADDS cannot
reach. ⚠ Under hybrid capture the premise weakens — the probes reshape the RNA lengths, so the two
contaminants stop resembling each other — and the contrast is applied anyway; measured against the
four-pool sum it is a large WIN there too, but for a reason it does not model, so read
:func:`_deconvolved_gdna_counts` before relying on it.

* **gDNA** = **all FOUR** gDNA pools — intergenic contained, intronic contained, intron-exon crossing,
  intergenic-exon crossing — each divided by **its own** opportunity and then combined
  (:mod:`rigel.calibration.gdna_opportunity`). The contained pair dominates off capture; the crossing
  pair dominates under it, because a fragment beside a probe *reaches* the exon boundary and stops being
  contained. Fitting from the contained pair alone measures the short half of one population.
  ⛔ **Pooling the four RAW is a different operation and it is wrong**: the contained opportunity falls
  with length and the crossing opportunity rises, and one divisor over the sum read a gDNA mean of
  **146.05** where the contained pool said **88.0**. Divide each, then combine.
* **RNA** = the annotated-sj pool, splice OBSERVED. gDNA cannot be spliced.
  ⚠ ``sj_implicit`` fragments are already excluded by the accumulator — a splice that was never
  sequenced is a product of the very model this pool is used to fit.

The pool axis itself lives in :mod:`rigel.scan_payload`, with the schema, because it is the
accumulator's own enum and a private copy here is how three files come to disagree about which row is
which.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .gdna_density import contained_opportunity, one_sided_rate
from .signature import RegionType
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

__all__ = [
    "FLModels",
    "GdnaContrast",
    "POOL_EB_PRIOR_ESS",
    "build_fl_models",
    "gdna_contained_fl_mass",
    "gdna_fl_mass",
    "rna_fl_mass",
    "splash_fl_mass",
]

#: gDNA contained in exactly one intergenic or intronic region. Dominant OFF capture.
_GDNA_CONTAINED_POOLS = (POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC)

#: gDNA crossing exactly one boundary whose flanks are {intron, exon} or {intergenic, exon}. ⭐ Dominant
#: UNDER capture: a fragment beside a probe reaches the exon boundary, so it leaves the contained pools
#: and arrives here. Mature RNA never crosses an exon<->intron boundary, so these are gDNA by construction.
_GDNA_CROSSING_POOLS = (POOL_DNA_INTRON_EXON, POOL_DNA_INTERGENIC_EXON)

#: ⭐ All four, in ``rigel.scan_payload`` pool order so they pair 1:1 with ``GdnaOpportunity.pools``.
_GDNA_POOLS = _GDNA_CONTAINED_POOLS + _GDNA_CROSSING_POOLS

#: The pure RNA pool: an OBSERVED splice across an annotated sj.
_RNA_POOLS = (POOL_RNA_SPLICED,)

#: Kept as a name because the report shows on-target gDNA separately; it is now also FITTED, via
#: ``_GDNA_CROSSING_POOLS``.
_SPLASH_POOLS = _GDNA_CROSSING_POOLS

#: Dirichlet pseudo-count for the smooth EB shrink toward the global FL. Not a
#: cliff: ``pool_total ≫ prior_ess`` → empirical; ``≪`` → the global anchor; ``= 0``
#: → global. Revisit the default on real-data pool sizes (PR04c decision 5).
POOL_EB_PRIOR_ESS: float = 1000.0


@dataclass(frozen=True, slots=True)
class GdnaContrast:
    """What the two-pool contrast did, or why it declined — QC, never an input to anything.

    ⭐ ``rate_over_pooled`` is the number to read: how much contamination the density fit found. A value
    near 1 says the pools were already clean, which is a measurement and not an inaction.
    """

    applied: bool
    declined_because: str
    gdna_density: float
    rate_over_pooled: float
    intergenic_gdna_share: float
    intronic_gdna_share: float
    separation: float


@dataclass(frozen=True, slots=True)
class FLModels:
    """Library-wide FL distributions (float64[max_size + 1]) + their pool totals.

    Two views of each distribution are carried:

    * ``*_pmf`` — the EB-smoothed pmf used for **scoring / calibration** (RNA and
      gDNA are shrunk toward ``global_pmf``; see :func:`build_fl_models`).
    * ``*_counts`` — the **raw, unsmoothed** histograms (aligned to
      ``max_size + 1`` bins) the pmfs were built from. These are the honest
      empirical distributions surfaced as QC (:meth:`rna_model` / :meth:`gdna_model`
      / :meth:`global_model`); the EB smoothing is a <1% perturbation at real
      library scale and stays internal to scoring.

    ``gdna_counts`` is the pure contained-pool histogram (:func:`gdna_fl_mass`) — intergenic + intronic,
    so gDNA *plus* whatever nascent RNA sits in an intron, not a deconvolved pure-gDNA distribution.
    ``rna_counts`` is the annotated-sj histogram (:func:`rna_fl_mass`).
    """

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # spliced, EB-shrunk toward global
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    global_counts: np.ndarray  # the accumulator's deposited_lengths — every deposited fragment at L
    rna_counts: np.ndarray  # raw spliced-annotated histogram
    gdna_counts: np.ndarray  # raw structural pool (intergenic + intronic)
    n_global: float
    n_rna: float
    n_gdna: float
    max_size: int
    #: float64[N_FRAGMENT_POOLS, max_size + 1] — the five pure pools UNAGGREGATED, straight off the
    #: payload. ``rna_counts`` and ``gdna_counts`` above are sums of subsets of these rows; the two
    #: crossing pools (:func:`splash_fl_mass`, on-target gDNA) are in **neither** sum and appear
    #: nowhere else. ⭐ They are carried here so the report can show each pool separately — that
    #: function's docstring asks for exactly this ("makes that comparison an output instead of an
    #: assumption") and nothing was doing it.
    pool_counts: np.ndarray = None
    #: What the two-pool contrast did. ``None`` when it was never offered the per-region inputs (the
    #: second pass and most tests), which is a different state from having declined.
    gdna_contrast: "GdnaContrast | None" = None

    def _empirical(self, counts: np.ndarray) -> "FragmentLengthModel":
        """Wrap a raw count vector as an unfinalized ``FragmentLengthModel``.

        Unfinalized on purpose: summary stats (mean/std/median/mode) and
        ``to_dict()`` then read the *raw* histogram, matching how the scanner's
        global + per-category models are reported (no scoring-side pseudo-count
        smoothing enters the QC numbers).
        """
        from ..frag_length_model import FragmentLengthModel

        return FragmentLengthModel(max_size=self.max_size, counts=counts.copy())

    def global_model(self) -> "FragmentLengthModel":
        """Empirical global FL distribution as a ``FragmentLengthModel`` (QC)."""
        return self._empirical(self.global_counts)

    def rna_model(self) -> "FragmentLengthModel":
        """Empirical RNA (spliced-annotated) FL distribution (QC)."""
        return self._empirical(self.rna_counts)

    def gdna_model(self) -> "FragmentLengthModel":
        """Empirical gDNA structural-pool FL distribution (QC).

        The pool is the intergenic + intronic CONTAINED mass (:func:`gdna_fl_mass`): a gDNA-dominated
        proxy that also includes intronic/nascent RNA, not a deconvolved pure-gDNA distribution.
        """
        return self._empirical(self.gdna_counts)


def _pool_sum(payload: "AccumulatorPayload", pools) -> np.ndarray:
    """Sum the named rows of ``payload.pool_lengths`` into one float64 histogram over ``L``.

    ⭐ Binned at ``L``, the molecule length — not at the covered length. Binning at covered length
    collapses the gDNA histogram to a spike at twice the read length, so every long gDNA fragment scores
    as RNA. The accumulator's ``L`` already includes the mate gap and
    excludes region_bound introns, so it is the molecule length for both components under one rule.
    """
    pool_lengths = payload.pool_lengths
    if pool_lengths is None:
        return np.zeros(0, dtype=np.float64)
    return np.asarray(pool_lengths, dtype=np.float64)[list(pools)].sum(axis=0)


def gdna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The gDNA length histogram: **all four** pure gDNA pools, summed.

    ⛔ **This sum is only meaningful paired with :meth:`GdnaOpportunity.combined_probability`.** The four
    pools tilt in opposite directions, so the raw sum is biased long — that is the 146.05-against-88.0
    defect. :func:`build_fl_models` sums the counts and divides by the summed opportunity, which is a
    different operation; it never uses this histogram on its own.
    """
    return _pool_sum(payload, _GDNA_POOLS)


def gdna_contained_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """Only the two CONTAINED gDNA pools — the fallback when no annotation divisor is available.

    ⚠ Correct to within ~0.5 % off capture and ~15 % short under it, because capture moves the long half
    of the population into the crossing pools. Use the four-pool form whenever an index is at hand.
    """
    return _pool_sum(payload, _GDNA_CONTAINED_POOLS)


def rna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The pure RNA length histogram: fragments that used an annotated sj, splice OBSERVED."""
    return _pool_sum(payload, _RNA_POOLS)


def splash_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The ON-TARGET gDNA length histogram — the two crossing pools, for QC.

    ⚠ Never fold this into :func:`gdna_fl_mass`. Under capture the intergenic pool is *depleted, not
    impure*, so composition stays clean while coverage does not, and on-target gDNA fragments run ~42 bp
    shorter. A model fitted off-target is therefore mis-centred for exactly the fragments that leak —
    and having this as a named pool makes that comparison an output instead of an assumption.
    """
    return _pool_sum(payload, _SPLASH_POOLS)


def _deconvolved_gdna_counts(
    payload: "AccumulatorPayload",
    gdna_opportunity: "GdnaOpportunity",
    region_lengths: np.ndarray,
    region_types: np.ndarray,
) -> "tuple[np.ndarray | None, GdnaContrast]":
    """The two-pool contrast: the gDNA length histogram with the contaminant divided out.

    Returns ``(counts_or_None, diagnostics)``. ⭐ **``None`` means DECLINED and the caller must fall back
    to the four-pool sum** — declining is a real answer here, not a failure, and every decline is named in
    the returned :class:`GdnaContrast` so a run that corrected nothing cannot be mistaken for one that
    corrected everything.

    **The three declines, each derived rather than chosen:**

    * *no density* — the one-sided rate found no support, which is what a zero-gDNA library looks like.
      There is no gDNA length distribution to estimate and inventing one is the failure mode to avoid.
    * *purities not separated* — the contrast divides by ``a_0 - a_1``, so it needs the two pools to have
      measurably different compositions. ⛔ The comparison is against that difference's OWN sampling
      error, never a chosen floor: with ``a_p = rho*E_p/n_p`` and ``Var(n_p) ~ n_p``, the delta method
      gives ``Var(a_0 - a_1) ~ a_0^2/n_0 + a_1^2/n_1`` — the shared ``rho`` term is common-mode and
      cancels out of the difference — and the fit stands down when the separation does not exceed its own
      standard error. Same shape as the strand channel's derived noise-floor deadband. ⭐ This is what
      makes the estimator stand down by itself at a near-pure library, where nothing needs correcting.
    * *degenerate* — a pool with no fragments at all.

    ⭐⭐ **UNDER HYBRID CAPTURE THE CONTRAST DEGENERATES TO THE INTERGENIC POOL ALONE, AND THAT IS WHY IT
    IS SAFE THERE.** The premise is that both pools' contaminants share a length distribution; off capture
    they agree to a total variation of 0.06-0.14, under capture only 0.95-0.97. That would be alarming
    except that under capture the intergenic pool is **depleted, not impure** (its measured gDNA purity
    goes 0.48 -> 0.955 at ``g05``), so ``a_0 = rho*E_0/n_0`` exceeds 1 and CLIPS. Put ``a_0 = 1`` in the
    formula above and it collapses exactly::

        g = [(1-a_1) f_0 - 0] / (1 - a_1) = f_0

    ⭐ **The intronic pool's coefficient becomes zero, so its contaminant never enters the answer** — which
    is precisely the term the shared-``r`` premise was needed for. Verified: under capture the returned
    histogram equals the de-tilted intergenic pool to **3e-17**, machine epsilon. ⚠ So the estimator does
    not "get away with" a false premise; the clip removes the term that depended on it, and what is left is
    a one-pool estimate that capture happens to make nearly pure. ⛔ **The safety therefore rests on
    ``a_0`` clipping**, i.e. on the intergenic pool really being near-pure under capture; a probe panel
    that put RNA back into intergenic space would silently break it, and nothing here would notice.

    ⛔ **Two candidate detectors were built and BOTH refuted by measurement**, so do not re-propose them:
    projecting onto the non-negative cone (it reproduces the contrast wherever the contrast is already
    feasible, and degrades capture further where it is not), and testing the recovered components'
    negative mass against its own Poisson noise floor (the floor scales with the violation, so the ratio
    sits flat at 0.2-0.7 on and off capture and never separates).
    """
    types = np.asarray(region_types).ravel()
    lengths = np.asarray(region_lengths, dtype=np.float64).ravel()
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)
    # ⚠ genome-strand columns summed: gDNA is strand-symmetric and this estimator is about DENSITY, which
    # is why it also works on an unstranded library, where the strand axis carries nothing.
    counts = counts.sum(axis=1) if counts.ndim > 1 else counts
    n = min(types.size, lengths.size, counts.size)
    types, lengths, counts = types[:n], lengths[:n], counts[:n]

    is_ig = (types == int(RegionType.INTERGENIC)) & (lengths > 0.0)
    is_in = (types == int(RegionType.INTRON)) & (lengths > 0.0)
    pooled_pmf = _normalized(_pool_sum(payload, _GDNA_CONTAINED_POOLS))
    fit = one_sided_rate(
        counts[is_ig | is_in], contained_opportunity(pooled_pmf, lengths[is_ig | is_in])
    )
    if not fit.informative:
        return None, GdnaContrast(
            False, "no gDNA density", fit.rate, fit.rate_over_pooled, 0.0, 0.0, 0.0
        )

    weights, totals = [], []
    for mask in (is_ig, is_in):
        n_p = float(counts[mask].sum())
        e_p = float(contained_opportunity(pooled_pmf, lengths[mask]).sum())
        weights.append(min(fit.rate * e_p / n_p, 1.0) if n_p > 0.0 else 0.0)
        totals.append(n_p)
    a0, a1 = weights
    if totals[0] <= 0.0 or totals[1] <= 0.0:
        return None, GdnaContrast(
            False, "an empty pool", fit.rate, fit.rate_over_pooled, a0, a1, 0.0
        )

    sep = a0 - a1
    se = float(np.sqrt(a0 * a0 / totals[0] + a1 * a1 / totals[1]))
    if not abs(sep) > se:
        return None, GdnaContrast(
            False, "purities not separated", fit.rate, fit.rate_over_pooled, a0, a1, sep
        )

    pi = gdna_opportunity.combined_probability()
    del pi  # the combined divisor is for the four-pool sum; each pool here takes its OWN
    total = np.asarray(gdna_opportunity.total, dtype=np.float64)
    raw = np.asarray(payload.pool_lengths, dtype=np.float64)
    f = []
    for pool in _GDNA_CONTAINED_POOLS:
        opp = np.asarray(gdna_opportunity.pools[pool], dtype=np.float64)
        prob = np.zeros_like(opp)
        np.divide(opp, total, out=prob, where=total > 0.0)
        f.append(_normalized(detilt_pool(raw[pool], prob)))
    g = ((1.0 - a1) * f[0] - (1.0 - a0) * f[1]) / sep
    # ⚠ The negative excursions are sampling noise on a quantity that is a density; clipping is the
    # cheapest projection back onto the cone and was measured equal-or-better than a least-squares one.
    g = np.clip(g, 0.0, None)
    if not g.sum() > 0.0:
        return None, GdnaContrast(
            False, "empty after the contrast", fit.rate, fit.rate_over_pooled, a0, a1, sep
        )
    # ⭐ Rescaled to the pool mass the four-pool path would have carried, because the EB shrinkage reads
    # the TOTAL as "how much evidence stands behind this shape".
    return _normalized(g) * float(raw[list(_GDNA_POOLS)].sum()), GdnaContrast(
        True, "", fit.rate, fit.rate_over_pooled, a0, a1, sep
    )


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


def _smooth_eb(aligned: np.ndarray, global_pmf: np.ndarray, prior_ess: float):
    """Smooth EB pmf: ``(counts + prior_ess·global_pmf) / (total + prior_ess)``.

    ``aligned`` is an already ``max_size``-aligned count vector. Continuous in the
    pool total — no quality threshold (PR04c decision 5). Returns ``(pmf, pool_total)``.
    """
    total = float(aligned.sum())
    denom = total + prior_ess
    pmf = (aligned + prior_ess * global_pmf) / denom if denom > 0.0 else global_pmf.copy()
    return _normalized(pmf), total


def build_fl_models(
    payload: "AccumulatorPayload",
    *,
    sj_opportunity: np.ndarray | None = None,
    gdna_opportunity: "GdnaOpportunity | None" = None,
    region_lengths: np.ndarray | None = None,
    region_types: np.ndarray | None = None,
    prior_ess: float = POOL_EB_PRIOR_ESS,
) -> FLModels:
    """Build the global / RNA / gDNA FL pmfs from ONE payload, in ONE frame.

     ⭐ **All three histograms come off the same object, so they cannot disagree about what a
     fragment length IS.** The anchor is ``payload.deposited_lengths`` — every deposited fragment
     binned at its own ``L`` with no purity condition (TRAPS: a-purity-filter-is-a-length-filter); the two component pools are
     :func:`rna_fl_mass` and :func:`gdna_fl_mass`, drawn from exactly that same population. RNA and
     gDNA are EB-shrunk toward the anchor with a single Dirichlet ``prior_ess``.

     ⚠ **The payload is the only argument on purpose**. The
     anchor used to be passed in separately and was taken from the **scanner's** histogram, which
     measures fragment length by two other rules — ``frag.genomic_footprint()`` for one subset and a
     transcript-space length for a disjoint one — over a population that was never stated. That is
     accumulator-frame pools shrunk toward a scanner-frame anchor.
     in shipped code, and it is what made the length likelihood read a ruler mismatch as composition
    Removing the parameter is what makes the mixed-frame call unrepresentable
     rather than merely discouraged.

     ⚠ The anchor is unconditional **given deposit**, not unconditional: it excludes what the
     accumulator rejects (too long, ambiguous path, strand-undefined, empty), each counted in
     ``payload.qc``. That is precisely the population the pools are drawn from, which is what makes
     it the right anchor rather than a merely convenient one.

     ⭐ **Each component pool is divided by ITS OWN opportunity, and the two divisors are different
     objects because the two selections are different.**

     * ``sj_opportunity`` — ``pi(w)``, the chance a uniformly placed length-``w`` fragment crosses
       an annotated sj at all (:mod:`rigel.calibration.sj_opportunity`). The RNA pool is
       selected on *"used an annotated sj"*, which longer fragments do more often.
     * ``gdna_opportunity`` — the four gDNA pools' opportunities and the reference total
       (:mod:`rigel.calibration.gdna_opportunity`). Two of those pools are *contained in one region*, whose
       opportunity **falls** with length; two are *crossing exactly one boundary*, whose opportunity
       **rises**. ⛔ Folding one divisor into the other, or one divisor over the pooled sum, is a
       category error — it is the defect that read a gDNA mean of 146.05 where the contained pool said
       88.0.

     ⭐ **``region_lengths`` and ``region_types`` are what enable the two-pool CONTRAST** (the module
     docstring derives it; :func:`_deconvolved_gdna_counts` implements it and owns every reason it
     declines). Both come straight from
     :func:`~rigel.calibration.splice_graph.build_region_partition_arrays`, the same partition the
     scanner deposits into, so they cannot disagree with the banks they index. ⛔ **Omitting them is a
     supported state, not a degraded one** — the second pass and most tests do — and it falls back to the
     four-pool sum, which is what shipped before.

     ⚠ **``None`` means no annotation was offered, and the fallback is the honest one, not the
     convenient one**: the RNA pool stays tilted and the gDNA pool falls back to the CONTAINED pair
     alone (:func:`gdna_contained_fl_mass`). ⛔ It does **not** fall back to the four pools pooled raw,
     because that is measurably worse than either.

     For the EB kernel over three free histograms — the shape a unit test needs and production never
     has — see :func:`_fl_models_from_histograms`.
    """
    rna_counts = rna_fl_mass(payload)
    # ⚠ One de-tilt implementation, shared: it preserves the pool TOTAL (the EB shrinkage reads that as
    # "how much evidence stands behind this shape") and drops bins the opportunity says are impossible.
    from .sj_opportunity import detilt_pool

    if sj_opportunity is not None:
        rna_counts = detilt_pool(rna_counts, sj_opportunity)

    contrast = None
    if gdna_opportunity is None:
        gdna_counts = gdna_contained_fl_mass(payload)
    else:
        gdna_counts = detilt_pool(gdna_fl_mass(payload), gdna_opportunity.combined_probability())
        if region_lengths is not None and region_types is not None:
            deconvolved, contrast = _deconvolved_gdna_counts(
                payload, gdna_opportunity, region_lengths, region_types
            )
            if deconvolved is not None:
                gdna_counts = deconvolved

    return _fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=int(payload.max_length),
        prior_ess=prior_ess,
        pool_counts=payload.pool_lengths,
        gdna_contrast=contrast,
    )


def _fl_models_from_histograms(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    prior_ess: float = POOL_EB_PRIOR_ESS,
    pool_counts: np.ndarray | None = None,
    gdna_contrast: "GdnaContrast | None" = None,
) -> FLModels:
    """The smooth-EB kernel: three histograms in, three pmfs out.

    ⛔ **Not a production entry point.** Production has exactly one source for all three histograms
    and reaches it through :func:`build_fl_models`; this exists so the shrinkage policy itself can be
    exercised over anchors and pools that no real payload would produce.
    """
    global_aligned = _aligned(global_counts, max_size)
    rna_aligned = _aligned(rna_counts, max_size)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, n_rna = _smooth_eb(rna_aligned, global_pmf, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)
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
        n_global=float(global_aligned.sum()),
        n_rna=n_rna,
        n_gdna=n_gdna,
        max_size=int(max_size),
        gdna_contrast=gdna_contrast,
    )
