"""The gDNA and RNA fragment-length laws, and each law's sampling variance.

This module produces the library-wide length laws the calibrator's effective lengths and the EM's
per-fragment scorer read, from the accumulator's five length pools. It is not a per-fragment length
likelihood: that composition channel is deliberately outside calibration. Both laws drive per-region
effective lengths in the sweep (``region_geometry.build_region_geometry``) — gDNA opportunities use
the gDNA law, RNA opportunities the RNA law — and each is smoothly empirical-Bayes-shrunk toward the
global length distribution, so a thin pool falls back to the anchor without a threshold.

The five pools, as ``rigel.scan_payload`` orders them:

* two CONTAINED gDNA pools, intergenic and intronic — gDNA wholly inside one region. Dominant off
  capture.
* two CROSSING gDNA pools, {intron, exon} and {intergenic, exon} flanks — gDNA across exactly one
  boundary. Dominant under capture, because a fragment beside a probe reaches the exon boundary and
  so stops being contained. Mature RNA never crosses an exon/intron boundary, so these are gDNA by
  structure.
* one RNA pool: a splice across an annotated sj on the fragment's one surviving path, sequenced or
  implied. gDNA cannot splice, so this pool is certified RNA.

No pool is pure, and nothing here may assume one is. Against an origin-split oracle the intronic
pool is dominated by RNA inside introns and the intergenic pool carries RNA too, because
"intergenic" is whatever the annotation leaves over. The resulting bias is
``RNA_share x (len_RNA - len_gDNA)``, so it is invisible on any panel that gives the two components
equal fragment lengths.

What replaces purity is a two-pool CONTRAST, which needs no pure pool and no template. The two
contained pools share one opportunity geometry — RNA inside an intron and unannotated transcription
in intergenic space are genomically contiguous exactly as gDNA is — so after de-tilting each is the
same two shapes at a different mixing weight, and the contaminant cancels::

    f_0 = a_0*g + (1-a_0)*r        f_1 = a_1*g + (1-a_1)*r
    =>  g = [ (1-a_1)*f_0 - (1-a_0)*f_1 ] / (a_0 - a_1)

The divisor is the SEPARATION of the two purities, not a purity, so it does not blow up as a pool
gets dirty; and with ``f_0 == f_1`` it returns them unchanged, so a pool pair with nothing to say
changes nothing. The weights come from :mod:`rigel.calibration.gdna_density`: ``a_p = rho_g*E_p/n_p``
with ``rho_g`` read off the low side of the per-object density, where a contaminant that only ADDS
cannot reach. Under hybrid capture the premise weakens, because the probes reshape the RNA lengths
and the two contaminants stop resembling each other; the contrast is applied anyway, and
:func:`_deconvolved_gdna_counts` states what that costs.

Never pool the four gDNA histograms raw. The contained opportunity falls with length and the
crossing opportunity rises, so one divisor over the sum reads the gDNA mean far too long. Divide
each pool by its own opportunity (:mod:`rigel.calibration.gdna_opportunity`), then combine.

Every region is read over its OBSERVABLE pieces (`region_arrays.observable_pieces`): no fragment is observed
on a base of the index's unobservable mask, so a region's count is over its observable starts alone, and every
contained support here — the one-sided rate's exposure and the purities — is summed over those pieces. A region
masked whole carries no exposure and is no object: read over its whole length it would be the cleanest object of
all, a zero count against a large exposure, and the rate it pulls down collapses the purities and hands the
contrast's noise to the law.

``gdna_pmf`` is the UNIFORM-FRAME law, what a gDNA fragment looks like before capture selects it — the
contrast's alone, nothing read beside the probes mixed in: the one the opportunity and prior mathematics assumes.
The EM's two length tables are not built here. They are decided after calibration, from this law, ``rna_pmf`` and
their sampling variances (:mod:`rigel.calibration.shared_length_law`), because whether the two origins' laws differ
by more than capture can make needs calibration's capture efficiencies.

The pool axis itself lives in :mod:`rigel.scan_payload`, with the schema, because it is the
accumulator's own enum and a private copy here is how three files come to disagree about which row
is which.
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
    "gdna_contained_fl_mass",
    "gdna_fl_mass",
    "rna_fl_mass",
]

#: gDNA contained in exactly one intergenic or intronic region. Dominant OFF capture.
_GDNA_CONTAINED_POOLS = (POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC)

#: gDNA crossing exactly one boundary whose flanks are {intron, exon} or {intergenic, exon}. Dominant
#: UNDER capture: a fragment beside a probe reaches the exon boundary, so it leaves the contained pools
#: and arrives here. Mature RNA never crosses an exon/intron boundary, so these are gDNA by structure.
_GDNA_CROSSING_POOLS = (POOL_DNA_INTRON_EXON, POOL_DNA_INTERGENIC_EXON)

#: All four, in ``rigel.scan_payload`` pool order so they pair 1:1 with ``GdnaOpportunity.pools``.
_GDNA_POOLS = _GDNA_CONTAINED_POOLS + _GDNA_CROSSING_POOLS

#: The certified RNA pool: a splice across an annotated sj on the fragment's one surviving path.
_RNA_POOLS = (POOL_RNA_SPLICED,)


@dataclass(frozen=True, slots=True)
class GdnaContrast:
    """What the two-pool contrast did, or why it declined — test surface: nothing downstream reads it.

    ``rate_over_pooled`` is the number to read: how much contamination the density fit found. A
    value near 1 says the pools were already clean, which is a measurement and not an inaction.
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

    * ``*_pmf`` — the EB-smoothed pmf used for scoring and calibration; RNA and gDNA are shrunk
      toward ``global_pmf``, see :func:`build_fl_models`.
    * ``*_counts`` — the raw, unsmoothed histograms, aligned to ``max_size + 1`` bins, that the pmfs
      were built from. These are the honest empirical distributions surfaced as QC
      (:meth:`rna_model` / :meth:`gdna_model` / :meth:`global_model`); the EB smoothing is a small
      perturbation at real library scale and stays internal to scoring.

    ``gdna_counts`` is the uniform-frame gDNA histogram the pmf is shrunk from, as
    :func:`build_fl_models` builds it: the two-pool contrast when the region inputs are given, else the
    de-tilted four-pool sum, or the contained pair with no opportunity at all. ``rna_counts`` is the
    annotated-sj histogram (:func:`rna_fl_mass`).

    ``rna_pmf_var`` and ``gdna_pmf_var`` are each pmf's per-length sampling variance: its raw law's plug-in
    variance, scaled by the square of the empirical-Bayes weight ``n / (n + prior_ess)`` its counts carry
    against the anchor, and 0 at a length the raw law has no count at. The EM's two length tables read
    them (:mod:`rigel.calibration.shared_length_law`) to tell a difference from sampling noise.
    """

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # spliced, EB-shrunk toward global
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    global_counts: np.ndarray  # the accumulator's deposited_lengths — every deposited fragment at L
    rna_counts: np.ndarray  # raw spliced-annotated histogram
    gdna_counts: np.ndarray  # the uniform-frame gDNA histogram (see above)
    n_global: float
    n_rna: float
    n_gdna: float
    max_size: int
    #: float64[N_FRAGMENT_POOLS, max_size + 1] — the five pools UNAGGREGATED, straight off the
    #: payload. ``rna_counts`` and ``gdna_counts`` above are sums of subsets of these rows. They are
    #: carried separately so the report can show each pool on its own, which is what makes the
    #: off-target / on-target length comparison an output rather than an assumption.
    pool_counts: np.ndarray = None
    rna_pmf_var: np.ndarray = None
    gdna_pmf_var: np.ndarray = None

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
        """Empirical gDNA FL distribution (QC): ``gdna_counts``, the histogram the gDNA pmf is shrunk
        from."""
        return self._empirical(self.gdna_counts)


def _pool_sum(payload: "AccumulatorPayload", pools) -> np.ndarray:
    """Sum the named rows of ``payload.pool_lengths`` into one float64 histogram over ``L``.

    Binned at ``L``, the molecule length, never at the covered length. Binning at covered length
    collapses the gDNA histogram to a spike at twice the read length, so every long gDNA fragment
    scores as RNA. The accumulator's ``L`` already includes the mate gap and excludes spliced-out
    introns, so it is the molecule length for both components under one rule.
    """
    return np.asarray(payload.pool_lengths, dtype=np.float64)[list(pools)].sum(axis=0)


def gdna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The gDNA length histogram: all four structural gDNA pools, summed.

    This sum is only meaningful paired with :meth:`GdnaOpportunity.combined_probability`. The four
    pools tilt in opposite directions, so the raw sum is biased long. :func:`build_fl_models` sums
    the counts and divides by the summed opportunity, which is a different operation; it never uses
    this histogram on its own.
    """
    return _pool_sum(payload, _GDNA_POOLS)


def gdna_contained_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """Only the two CONTAINED gDNA pools — the fallback when no annotation divisor is available.

    Accurate off capture, but short under it, because capture moves the long half of the population
    into the crossing pools. Use the four-pool form whenever an index is at hand.
    """
    return _pool_sum(payload, _GDNA_CONTAINED_POOLS)


def rna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The certified RNA length histogram: fragments that used an annotated sj on their one surviving
    path."""
    return _pool_sum(payload, _RNA_POOLS)


def _resolution_weight(signal_sq: float, noise_sq: float) -> float:
    """``S / (S + N)`` — how much of an observed split is signal rather than its own sampling noise.

    This is what replaces every "is it big enough?" threshold in this module. A test of the form
    ``|split| > standard error`` is a cliff: the estimator's behaviour changes discontinuously as
    data accumulates, and where the cliff sits is a property of the sample rather than of the tool.
    This is the same comparison made continuously — the signal-to-noise weight, 0 when the split is
    pure noise, 1 when the noise vanishes, 1/2 exactly where such a threshold would sit, and
    monotone between. No constant is introduced: both arguments are variances the data supplies.
    """
    s = max(float(signal_sq), 0.0)
    n = max(float(noise_sq), 0.0)
    if s <= 0.0:
        return 0.0
    if n <= 0.0:
        return 1.0
    return s / (s + n)


def _region_counts(payload, regions: "RegionArrays") -> np.ndarray:
    """``float64[R]`` — each region's contained count, its genome-strand columns summed: gDNA is strand-symmetric and
    these estimators are about DENSITY, which is why they also work on an unstranded library, where the strand axis
    carries nothing. Refused unless *regions* is the partition the payload was deposited into."""
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)
    counts = counts.sum(axis=1) if counts.ndim > 1 else counts
    if counts.shape[0] != regions.n_regions:
        raise ValueError(
            f"the frame has {regions.n_regions:,} regions and the payload {counts.shape[0]:,}: they address "
            "different partitions. Pass RegionArrays.from_index of the index the scan used."
        )
    return counts


def _deconvolved_gdna_counts(
    payload: "AccumulatorPayload",
    gdna_opportunity: "GdnaOpportunity",
    regions: "RegionArrays",
) -> "tuple[np.ndarray | None, np.ndarray | None, GdnaContrast]":
    """The two-pool contrast: the gDNA length histogram with the contaminant divided out.

    Returns ``(counts_or_None, variance_or_None, diagnostics)``. ``None`` means DECLINED and the caller
    must fall back to the four-pool sum; declining is a real answer here, not a failure, and every
    decline is named in the returned :class:`GdnaContrast`, so a run that corrected nothing cannot be
    mistaken for one that corrected everything. ``variance`` is the normalised law's per-length plug-in
    variance through the contrast's own coefficients: each de-tilted pool's ``f_p² / c_p`` (de-tilting
    rescales a length, it adds no evidence; a length with no count carries none), the inversion and the
    mixture combined coefficient by coefficient, so their covariance is in it.

    The declines, each derived rather than chosen:

    * *no density* — the one-sided rate found no support, which is what a zero-gDNA library looks like.
      There is no gDNA length distribution to estimate and inventing one is the failure mode to avoid.
    * *an empty pool* — a pool with no fragments at all.
    * *empty after the contrast* — nothing survives the projection back onto the cone.

    Identical purities are NOT a decline. The contrast divides by ``a_0 - a_1``, and there is no
    threshold short of zero separation: the inversion is blended toward the pools' own mixture by the
    resolution weight ``sep^2 / (sep^2 + SE^2)``, where with ``a_p = rho*E_p/n_p`` and ``Var(n_p) ~ n_p``
    the delta method gives ``SE^2 ~ a_0^2/n_0 + a_1^2/n_1`` — the shared ``rho`` term is common-mode and
    cancels out of the difference. So the estimator fades out by itself at a near-pure library, where
    nothing needs correcting, and at a separation of exactly zero the weight is exactly zero and the
    answer is the mixture itself: the limit the fade approaches, and what the boundary pair returns at
    its own tie. Both purities clip at 1 in a near-pure captured library (``g98 ss.99 ON``), and a
    decline there handed the uniform-frame law to the capture-selected four-pool census, 245 bp against
    a true 217 (``tests/calibration/test_fl.py::test_the_uniform_law_is_continuous_through_a_purity_tie``).

    Under hybrid capture the contrast degenerates to the intergenic pool alone, and that is why it
    is safe there. The premise is that both pools' contaminants share a length distribution, and
    under capture that premise is false. It survives because under capture the intergenic pool is
    depleted rather than made impure, so ``a_0 = rho*E_0/n_0`` exceeds 1 and CLIPS. Put ``a_0 = 1``
    in the formula above and it collapses exactly::

        g = [(1-a_1) f_0 - 0] / (1 - a_1) = f_0

    The intronic pool's coefficient becomes zero, so its contaminant never enters the answer, and
    that is precisely the term the shared-contaminant premise was needed for; under capture the
    returned histogram equals the de-tilted intergenic pool to machine epsilon. So the estimator
    does not get away with a false premise — the clip removes the term that depended on it, and what
    is left is a one-pool estimate that capture happens to make nearly pure. The safety therefore
    rests on ``a_0`` clipping, i.e. on the intergenic pool really being near-pure under capture; a
    probe panel that put RNA back into intergenic space would break it silently, and nothing here
    would notice.

    Two candidate detectors for that failure were built and both refuted by measurement, so do not
    re-propose them: projecting onto the non-negative cone, and testing the recovered components'
    negative mass against its own Poisson noise floor.
    """
    types = coarse_type_array(regions.signature)
    counts = _region_counts(payload, regions)
    observable = observable_length(regions)

    is_ig = (types == int(RegionType.INTERGENIC)) & (observable > 0.0)
    is_in = (types == int(RegionType.INTRON)) & (observable > 0.0)
    # the contained pools' own law prices the regions they come from
    support = observable_contained_eff_length(regions, _pool_sum(payload, _GDNA_CONTAINED_POOLS))
    fit = one_sided_rate(counts[is_ig | is_in], support[is_ig | is_in])
    if not fit.informative:
        return (
            None,
            None,
            GdnaContrast(False, "no gDNA density", fit.rate, fit.rate_over_pooled, 0.0, 0.0, 0.0),
        )

    weights, totals = [], []
    for mask in (is_ig, is_in):
        n_p = float(counts[mask].sum())
        e_p = float(support[mask].sum())
        weights.append(min(fit.rate * e_p / n_p, 1.0) if n_p > 0.0 else 0.0)
        totals.append(n_p)
    a0, a1 = weights
    if totals[0] <= 0.0 or totals[1] <= 0.0:
        return (
            None,
            None,
            GdnaContrast(False, "an empty pool", fit.rate, fit.rate_over_pooled, a0, a1, 0.0),
        )

    sep = a0 - a1
    se = float(np.sqrt(a0 * a0 / totals[0] + a1 * a1 / totals[1]))
    # No threshold. The inversion is blended toward the pools' own mixture by how well the purity
    # separation is resolved against its own standard error, so a pair that says nothing changes
    # nothing and a pair that says a little changes a little. `lam = 1` accepts the inversion whole
    # and `lam = 0` declines it, joined continuously rather than switched between.
    lam_sep = _resolution_weight(sep * sep, se * se)

    # each contained pool takes its OWN opportunity, never the combined divisor of the four-pool sum
    total = np.asarray(gdna_opportunity.total, dtype=np.float64)
    raw = np.asarray(payload.pool_lengths, dtype=np.float64)
    f = []
    for pool in _GDNA_CONTAINED_POOLS:
        opp = np.asarray(gdna_opportunity.pools[pool], dtype=np.float64)
        prob = np.zeros_like(opp)
        np.divide(opp, total, out=prob, where=total > 0.0)
        f.append(_normalized(detilt_pool(raw[pool], prob)))
    mixture = _normalized(totals[0] * f[0] + totals[1] * f[1])
    # The law's per-length coefficients on the two pools, for its variance: the estimator differentiated as it is,
    # the mixture's weight everywhere and the inverse's only where the clip leaves it open, through the renormalisation.
    w0, w1 = totals[0] / (totals[0] + totals[1]), totals[1] / (totals[0] + totals[1])
    c0, c1 = np.full(f[0].size, w0), np.full(f[0].size, w1)
    if sep == 0.0:
        # the fade's own limit: at zero separation the resolution weight is zero
        g = mixture
    else:
        inverse = ((1.0 - a1) * f[0] - (1.0 - a0) * f[1]) / sep
        # The negative excursions are sampling noise on a quantity that is a density; clipping is the
        # cheapest projection back onto the cone, and measures no worse than a least-squares one.
        g = np.clip(inverse, 0.0, None)
        s = g.sum()
        if s > 0.0:
            open_ = (inverse > 0.0) * (lam_sep / (s * sep))
            c0 = (1.0 - lam_sep) * w0 + open_ * (1.0 - a1)
            c1 = (1.0 - lam_sep) * w1 - open_ * (1.0 - a0)
        g = mixture + lam_sep * (_normalized(g) - mixture) if s > 0.0 else mixture
    g = np.clip(g, 0.0, None)
    if not g.sum() > 0.0:
        return (
            None,
            None,
            GdnaContrast(
                False, "empty after the contrast", fit.rate, fit.rate_over_pooled, a0, a1, sep
            ),
        )
    variance = c0 * c0 * _plugin_variance(
        f[0], raw[_GDNA_CONTAINED_POOLS[0]]
    ) + c1 * c1 * _plugin_variance(f[1], raw[_GDNA_CONTAINED_POOLS[1]])
    # Rescaled to the pool mass the four-pool path would have carried, because the EB shrinkage reads
    # the TOTAL as "how much evidence stands behind this shape".
    return (
        _normalized(g) * float(raw[list(_GDNA_POOLS)].sum()),
        variance,
        GdnaContrast(True, "", fit.rate, fit.rate_over_pooled, a0, a1, sep),
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


def _plugin_variance(law: np.ndarray, raw: np.ndarray) -> np.ndarray:
    """A normalised law's per-length plug-in variance ``law² / raw``: a length's relative variance is one over the raw
    count behind it, which de-tilting rescales and does not add to. A length with no raw count carries none."""
    law = np.asarray(law, dtype=np.float64)
    n = min(law.size, np.asarray(raw).size)
    c = np.zeros(law.size)
    c[:n] = np.asarray(raw, dtype=np.float64)[:n]
    return np.divide(law * law, c, out=np.zeros(law.size), where=c > 0.0)


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
    sj_opportunity: np.ndarray | None = None,
    gdna_opportunity: "GdnaOpportunity | None" = None,
    regions: "RegionArrays | None" = None,
    prior_ess: float = CONSTANTS.fragment_length.pool_prior_ess,
) -> FLModels:
    """Build the global / RNA / gDNA FL pmfs from ONE payload, in ONE frame.

    All three histograms come off the same object, so they cannot disagree about what a fragment
    length IS. The anchor is ``payload.deposited_lengths`` — every deposited fragment binned at its
    own ``L`` with no purity condition — and the two
    component pools, :func:`rna_fl_mass` and :func:`gdna_fl_mass`, are drawn from exactly that
    population. RNA and gDNA are EB-shrunk toward the anchor with a single Dirichlet ``prior_ess``.

    The payload is the only source of histograms on purpose. An anchor passed in separately would
    come from the scanner, which measures fragment length by other rules over a population nobody
    states, and shrinking accumulator-frame pools toward a scanner-frame anchor makes a ruler
    mismatch read as composition. Having no such parameter is what makes that call unrepresentable
    rather than merely discouraged.

    The anchor is unconditional GIVEN DEPOSIT, not unconditional: it excludes what the accumulator
    rejects (too long, ambiguous path, strand-undefined, empty), each counted in ``payload.qc``.
    That is precisely the population the pools are drawn from, which is what makes it the right
    anchor rather than a merely convenient one.

    Each component pool is divided by its OWN opportunity, and the two divisors are different
    objects because the two selections are different:

    * ``sj_opportunity`` — ``pi(w)``, the chance a uniformly placed length-``w`` fragment crosses an
      annotated sj at all (:mod:`rigel.calibration.sj_opportunity`). The RNA pool is selected on
      "used an annotated sj", which longer fragments do more often.
    * ``gdna_opportunity`` — the four gDNA pools' opportunities and the reference total
      (:mod:`rigel.calibration.gdna_opportunity`). Two of those pools are contained in one region,
      whose opportunity falls with length; two are crossing exactly one boundary, whose opportunity
      rises. Folding one divisor into the other, or applying one divisor over the pooled sum, is a
      category error that reads the gDNA mean far too long.

    ``regions`` is what enables the two-pool CONTRAST: the module docstring derives it, and
    :func:`_deconvolved_gdna_counts` implements it and owns every reason it declines. It is the index's
    frame (:meth:`~rigel.calibration.region_arrays.RegionArrays.from_index`) — the partition the scanner
    deposits into, so it cannot disagree with the banks it indexes, carrying the unobservable mask every
    contained support here is read over.

    For the EB kernel over three free histograms — the shape a unit test needs and production never
    has — see :func:`_fl_models_from_histograms`.
    """
    rna_raw = rna_fl_mass(payload)
    rna_counts = rna_raw
    # One de-tilt implementation, shared: it preserves the pool TOTAL (the EB shrinkage reads that as
    # "how much evidence stands behind this shape") and drops bins the opportunity says are impossible.
    if sj_opportunity is not None:
        rna_counts = detilt_pool(rna_counts, sj_opportunity)

    if gdna_opportunity is None:
        gdna_counts = gdna_contained_fl_mass(payload)
        gdna_var = _plugin_variance(_normalized(gdna_counts), gdna_counts)
    else:
        gdna_raw = gdna_fl_mass(payload)
        gdna_counts = detilt_pool(gdna_raw, gdna_opportunity.combined_probability())
        gdna_var = _plugin_variance(_normalized(gdna_counts), gdna_raw)
        if regions is not None:
            deconvolved, deconvolved_var, _contrast = _deconvolved_gdna_counts(
                payload, gdna_opportunity, regions
            )
            if deconvolved is not None:
                gdna_counts, gdna_var = deconvolved, deconvolved_var

    return _fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=int(payload.max_length),
        prior_ess=prior_ess,
        pool_counts=payload.pool_lengths,
        rna_var=_plugin_variance(_normalized(rna_counts), rna_raw),
        gdna_var=gdna_var,
    )


def _fl_models_from_histograms(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    prior_ess: float = CONSTANTS.fragment_length.pool_prior_ess,
    pool_counts: np.ndarray | None = None,
    rna_var: np.ndarray | None = None,
    gdna_var: np.ndarray | None = None,
) -> FLModels:
    """The smooth-EB kernel: three histograms in, three pmfs out, with each component pmf's sampling variance.

    Not a production entry point. Production has exactly one source for all three histograms and
    reaches it through :func:`build_fl_models`; this exists so the shrinkage policy itself can be
    exercised over anchors and pools that no real payload would produce. ``rna_var`` / ``gdna_var`` are
    each normalised raw law's per-length variance (:func:`build_fl_models` passes the estimators' own);
    given none, a histogram is read as raw counts.
    """
    global_aligned = _aligned(global_counts, max_size)
    rna_aligned = _aligned(rna_counts, max_size)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, n_rna = _smooth_eb(rna_aligned, global_pmf, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)

    def pmf_var(var, aligned, n):
        raw = (
            _plugin_variance(_normalized(aligned), aligned)
            if var is None
            else _aligned(var, max_size)
        )
        weight = n / (n + prior_ess) if n + prior_ess > 0.0 else 0.0
        return weight * weight * raw

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
        rna_pmf_var=pmf_var(rna_var, rna_aligned, n_rna),
        gdna_pmf_var=pmf_var(gdna_var, gdna_aligned, n_gdna),
    )
