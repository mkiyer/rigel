"""The gDNA and RNA fragment-length distributions, and the two gDNA estimands.

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
* one RNA pool: an OBSERVED splice across an annotated sj. gDNA cannot splice. Implicit splices are
  excluded by the accumulator, since a splice that was never sequenced is a product of the model
  this pool is used to fit.

No pool is pure, and nothing here may assume one is. Against an origin-split oracle the intronic
pool is dominated by RNA inside introns and the intergenic pool carries mature RNA too, because
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

TWO gDNA ESTIMANDS, and routing them apart is the point:

* ``gdna_pmf`` is the UNIFORM-FRAME law, what a gDNA fragment looks like before capture selects it.
  This is the one the opportunity and prior mathematics assumes, and the only one geometry may read.
* ``gdna_realized_pmf`` is the LIBRARY-CENSUS law, what a SEQUENCED gDNA fragment looks like,
  capture selection included. This is the one the EM's per-fragment scorer conditions on.

Off capture the two coincide exactly, and the realized field is never ``None``, so the scorer reads
it unconditionally. Feeding the realized law to geometry is a large regression; one name over the
two quantities is how that happens.

The pool axis itself lives in :mod:`rigel.scan_payload`, with the schema, because it is the
accumulator's own enum and a private copy here is how three files come to disagree about which row
is which.
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
    "GdnaRealized",
    "POOL_EB_PRIOR_ESS",
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

#: The pure RNA pool: an OBSERVED splice across an annotated sj.
_RNA_POOLS = (POOL_RNA_SPLICED,)

#: Dirichlet pseudo-count for the smooth EB shrink toward the global length law. Not a cliff: a pool
#: total far above it gives the empirical law, far below it the global anchor, and 0 the anchor
#: exactly.
POOL_EB_PRIOR_ESS: float = 1000.0


@dataclass(frozen=True, slots=True)
class GdnaContrast:
    """What the two-pool contrast did, or why it declined — QC, never an input to anything.

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
class GdnaRealized:
    """What the LIBRARY-CENSUS (realized) gDNA law's estimator did, or why it declined — QC only.

    ``ontarget_share`` is the number to read: the fraction of the realized law's mass carried by the
    EXCESS-enrichment exon classes — the capture-only part, identically 0 when the boundaries carry no
    enrichment excess. ``boundary_share`` is the sampled crossing pair's mass, which is legitimately
    nonzero at every capture level. Near-zero ``ontarget_share`` says "this library shows no capture
    excess and the two estimands coincide", which is a measurement, not an inaction.
    """

    applied: bool
    declined_because: str
    boundary_share: float
    ontarget_share: float
    intron_exon_share: float
    intergenic_exon_share: float


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
    #: What the two-pool contrast did. ``None`` when it was never offered the per-region inputs (the
    #: second pass and most tests), which is a different state from having declined.
    gdna_contrast: "GdnaContrast | None" = None
    #: The second estimand — the LIBRARY-CENSUS law: what a sequenced gDNA fragment looks like,
    #: capture selection included. ``gdna_pmf`` above is the UNIFORM-FRAME law the opportunity and
    #: prior mathematics assumes; this one is what the EM's per-fragment scorer conditions on. Off
    #: capture the two coincide and this field equals ``gdna_pmf`` exactly; it is never ``None``, so
    #: the scorer reads it unconditionally. Feeding this law to GEOMETRY is a large regression, and
    #: one name over the two quantities is how that happens.
    gdna_realized_pmf: np.ndarray = None
    gdna_realized: "GdnaRealized | None" = None

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
    """The pure RNA length histogram: fragments that used an annotated sj, splice OBSERVED."""
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


def _couple_estimands(g_uniform, mass_uniform: float, g_boundary, mass_boundary: float):
    """Couple the two gDNA length estimands so they CONVERGE when their split is not measurable.

    The problem this solves: capture is a spectrum, and at either end one of the two strata has no
    data — at zero capture the boundary pools are nearly empty, and at very strong capture the
    contained pools are. A hard "use the realized law, else fall back to the uniform law" switch is
    then both a cliff and a lie, because the estimates do not merely become uncertain, they become
    the SAME estimate: nothing in the data distinguishes them any more.

    So the split is carried explicitly and weighted by how well it is resolved. With ``M`` the total
    gDNA mass and ``lam`` the resolution weight of ``Δ = g_boundary − g_uniform``::

        g_common = (m_u·g_uniform + m_b·g_boundary) / M          the precision-weighted consensus
        uniform  = g_uniform + (1 − lam)·(g_common − g_uniform)   shrinks toward consensus as lam → 0
        realized = g_common                                       the census, always the mixture

    At ``lam = 1`` this is exactly the uncoupled behaviour: ``uniform`` is the contained-derived law
    and ``realized`` the mass-weighted census. At ``lam = 0`` both are ``g_common`` — identical
    arrays, whichever stratum is the starved one, because ``g_common`` is precision-weighted and the
    starved stratum contributes nothing to it. The noise term ``1/m_u + 1/m_b`` diverges when either
    mass collapses, which is what makes the convergence automatic from both ends.

    Returns ``(uniform, realized, lam)``; the caller adds any capture-only correction to ``realized``,
    itself scaled by ``lam`` so that it too vanishes when the split is unresolved.
    """
    g_u = np.asarray(g_uniform, dtype=np.float64)
    g_b = np.asarray(g_boundary, dtype=np.float64)
    n = min(g_u.size, g_b.size)
    g_u, g_b = g_u[:n], g_b[:n]
    m_u, m_b = max(float(mass_uniform), 0.0), max(float(mass_boundary), 0.0)
    total = m_u + m_b
    if total <= 0.0:
        return g_u.copy(), g_u.copy(), 0.0
    g_common = (m_u * g_u + m_b * g_b) / total
    delta = g_b - g_u
    # the split's own sampling variance: a multinomial pmf from `m` counts has Var ~ p/m per bin, so
    # summing over bins gives 1/m for a normalised law. It DIVERGES when either mass collapses.
    noise = (1.0 / m_u if m_u > 0.0 else np.inf) + (1.0 / m_b if m_b > 0.0 else np.inf)
    lam = _resolution_weight(float((delta * delta).sum()), noise)
    uniform = g_u + (1.0 - lam) * (g_common - g_u)
    return uniform, g_common, lam


def _deconvolved_gdna_counts(
    payload: "AccumulatorPayload",
    gdna_opportunity: "GdnaOpportunity",
    region_lengths: np.ndarray,
    region_types: np.ndarray,
) -> "tuple[np.ndarray | None, GdnaContrast]":
    """The two-pool contrast: the gDNA length histogram with the contaminant divided out.

    Returns ``(counts_or_None, diagnostics)``. ``None`` means DECLINED and the caller must fall back
    to the four-pool sum; declining is a real answer here, not a failure, and every decline is named
    in the returned :class:`GdnaContrast`, so a run that corrected nothing cannot be mistaken for
    one that corrected everything.

    The declines, each derived rather than chosen:

    * *no density* — the one-sided rate found no support, which is what a zero-gDNA library looks like.
      There is no gDNA length distribution to estimate and inventing one is the failure mode to avoid.
    * *an empty pool* — a pool with no fragments at all.
    * *purities identical* — the contrast divides by ``a_0 - a_1``. Short of zero separation there is
      no threshold: the inversion is blended toward the pools' own mixture by the resolution weight
      ``sep^2 / (sep^2 + SE^2)``, where with ``a_p = rho*E_p/n_p`` and ``Var(n_p) ~ n_p`` the delta
      method gives ``SE^2 ~ a_0^2/n_0 + a_1^2/n_1`` — the shared ``rho`` term is common-mode and
      cancels out of the difference. So the estimator fades out by itself at a near-pure library,
      where nothing needs correcting.
    * *empty after the contrast* — nothing survives the projection back onto the cone.

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
    types = np.asarray(region_types).ravel()
    lengths = np.asarray(region_lengths, dtype=np.float64).ravel()
    counts = np.asarray(payload.region_contained_count, dtype=np.float64)
    # genome-strand columns summed: gDNA is strand-symmetric and this estimator is about DENSITY, which
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
    # No threshold. The inversion is blended toward the pools' own mixture by how well the purity
    # separation is resolved against its own standard error, so a pair that says nothing changes
    # nothing and a pair that says a little changes a little. `lam = 1` accepts the inversion whole
    # and `lam = 0` declines it, joined continuously rather than switched between.
    lam_sep = _resolution_weight(sep * sep, se * se)
    if sep == 0.0:
        return None, GdnaContrast(
            False, "purities identical", fit.rate, fit.rate_over_pooled, a0, a1, sep
        )

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
    g = ((1.0 - a1) * f[0] - (1.0 - a0) * f[1]) / sep
    # The negative excursions are sampling noise on a quantity that is a density; clipping is the
    # cheapest projection back onto the cone, and measures no worse than a least-squares one.
    g = np.clip(g, 0.0, None)
    s = g.sum()
    g = mixture + lam_sep * (_normalized(g) - mixture) if s > 0.0 else mixture
    g = np.clip(g, 0.0, None)
    if not g.sum() > 0.0:
        return None, GdnaContrast(
            False, "empty after the contrast", fit.rate, fit.rate_over_pooled, a0, a1, sep
        )
    # Rescaled to the pool mass the four-pool path would have carried, because the EB shrinkage reads
    # the TOTAL as "how much evidence stands behind this shape".
    return _normalized(g) * float(raw[list(_GDNA_POOLS)].sum()), GdnaContrast(
        True, "", fit.rate, fit.rate_over_pooled, a0, a1, sep
    )


def _adjacent_pairs(region_offsets, boundary_offsets) -> tuple[np.ndarray, np.ndarray]:
    """Every adjacent region pair of every reference, as ``(left region, the boundary between them)``.

    A reference contributes its regions in order and one boundary between each consecutive pair, so the
    pair table is a property of the PARTITION and not of any fit — built once, read by both passes of
    the census. A reference with fewer than two regions contributes nothing, and one whose boundary
    count disagrees with its region count is skipped rather than guessed at: the two axes are built by
    the same partition, so a disagreement is a corrupt payload, not a case to interpret.
    """
    roff = np.asarray(region_offsets, dtype=np.int64)
    boff = np.asarray(boundary_offsets, dtype=np.int64)
    n_regions = np.diff(roff)
    n_bounds = np.diff(boff)
    usable = (n_regions >= 2) & (n_bounds == n_regions - 1)
    if not usable.any():
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    counts = (n_regions - 1)[usable]
    starts_r = roff[:-1][usable]
    starts_b = boff[:-1][usable]
    # one contiguous run per reference, concatenated in reference order — the order the loop walked
    within = np.arange(int(counts.sum()), dtype=np.int64) - np.repeat(
        np.concatenate(([0], np.cumsum(counts)[:-1])), counts
    )
    return np.repeat(starts_r, counts) + within, np.repeat(starts_b, counts) + within


def _realized_gdna_counts(
    payload,
    gdna_opportunity: "GdnaOpportunity",
    region_lengths: np.ndarray,
    region_types: np.ndarray,
    rna_pmf: np.ndarray,
    uniform_counts: np.ndarray,
) -> "tuple[np.ndarray | None, np.ndarray | None, GdnaRealized]":
    """The LIBRARY-CENSUS gDNA law: the uniform-frame estimate plus everything capture SELECTED.

    Three ingredients, blended by realized gDNA mass so the capture spectrum needs no switch:

    * the OFF-TARGET stratum — the uniform-frame law ``phi`` (``uniform_counts``, stage 1's output),
      weighted by the contained pools' own gDNA mass;
    * the BOUNDARY stratum — the two crossing pools, deconvolved by the same two-shape contrast, with
      per-boundary composition CALIBRATED BY THE REGIONS: probes bind nucleic acid indiscriminately, so
      at one boundary gDNA and nascent share the enrichment and it CANCELS from the composition —
      ``a_b = 1/(1 + R_b)`` with ``R_b`` the RNA:gDNA odds transported from the adjacent region's own
      one-sided RNA excess. The pair is solved COMPLEMENT-FIRST (the contaminant clipped to the cone,
      then subtracted), which damps the ``1/sep`` noise by ``(1 - a_mix)`` — small exactly when the
      separation is small, so the conditioning fix is the algebra's own;
    * the ON-TARGET EXCESS — the exon classes no pool samples (contained-in-exon, spanning-the-exon),
      priced by each exon's own boundaries' measured enrichment and entering at ``(eps - 1)+`` ONLY:
      the sampled union already carries every class's uniform part, so at ``eps = 1`` (no capture) the
      correction vanishes IDENTICALLY and the realized law collapses to the sampled blend.

    Returns ``(realized_counts, uniform_counts, diagnostics)`` — BOTH estimands, because they are
    coupled: :func:`_couple_estimands` shrinks them toward one another by how well their split is
    resolved, so the uniform law returned here may differ from the ``uniform_counts`` handed in when the
    contained stratum is thin. That coupling is what removes the last cliff: there is no data at
    which behaviour switches, only a weight that fades.

    Declining is still a real answer at literally zero gDNA — there is no census to take — and the
    caller then keeps the uniform-frame law for both estimands.
    """
    obs_pools = np.asarray(payload.pool_lengths, dtype=np.float64)
    ty = np.asarray(region_types).ravel()
    ell = np.asarray(region_lengths, dtype=np.float64).ravel()
    cnt = np.asarray(payload.region_contained_count, dtype=np.float64)
    cnt = cnt.sum(axis=1) if cnt.ndim > 1 else cnt
    n = min(ty.size, ell.size, cnt.size)
    ty, ell, cnt = ty[:n], ell[:n], cnt[:n]
    is_ig = (ty == int(RegionType.INTERGENIC)) & (ell > 0.0)
    is_in = (ty == int(RegionType.INTRON)) & (ell > 0.0)

    g_C = _normalized(np.asarray(uniform_counts, dtype=np.float64))
    L_axis = np.arange(g_C.size, dtype=np.float64)
    e_g = contained_opportunity(g_C, ell)
    fit = one_sided_rate(cnt[is_ig | is_in], e_g[is_ig | is_in])
    if not fit.informative:
        return None, None, GdnaRealized(False, "no gDNA density", 0.0, 0.0, 0.0, 0.0)
    rho_off = fit.rate

    m_C = 0.0
    for mask in (is_ig, is_in):
        n_p = float(cnt[mask].sum())
        if n_p > 0.0:
            m_C += min(rho_off * float(e_g[mask].sum()) / n_p, 1.0) * n_p

    # ── the boundary stratum: per-boundary composition, regions calibrating boundaries
    rna = np.asarray(rna_pmf, dtype=np.float64)
    mu_r = float((rna * np.arange(rna.size)).sum() / max(rna.sum(), 1e-30))
    e_r = contained_opportunity(rna[: g_C.size], ell)
    excess = np.clip(cnt - rho_off * e_g, 0.0, None)
    rho_r = np.zeros_like(excess)
    np.divide(excess, e_r, out=rho_r, where=e_r > 0.0)

    bnd = np.asarray(payload.boundary_unspliced_count, dtype=np.float64)
    bnd = bnd.sum(axis=1) if bnd.ndim > 1 else bnd
    roff = np.asarray(payload.ref_region_offsets, dtype=np.int64)
    boff = np.asarray(payload.ref_boundary_offsets, dtype=np.int64)
    exon = int(RegionType.EXON)

    # THE EXON-FLANKING PAIRS, as arrays. Walking every adjacent region pair of every reference in
    # Python costs a second a fit on a human index and says nothing a mask cannot: the pair table is a
    # property of the partition, so it is built once, and the arithmetic over it is the same sums in the
    # same order (`np.bincount` accumulates in input order, which is the order the loop visited).
    left, boundary = _adjacent_pairs(roff, boff)
    right = left + 1
    tl, tr = ty[left], ty[right]
    l_exon, r_exon = tl == exon, tr == exon
    keep = (l_exon & ~r_exon) | (r_exon & ~l_exon)
    nb_all = bnd[boundary]
    keep &= nb_all > 0.0
    left, right, boundary = left[keep], right[keep], boundary[keep]
    l_exon = l_exon[keep]
    adj = np.where(l_exon, right, left)
    other = np.where(l_exon, ty[right], ty[left])
    exon_region = np.where(l_exon, left, right)
    nb_pair = bnd[boundary]
    # class 2 is an exon against an INTRON, class 3 an exon against anything else
    cls = np.where(other == int(RegionType.INTRON), 0, 1)
    rho_adj = rho_r[adj]
    n_exons = ty.shape[0]

    mu_g = float((g_C * L_axis).sum())
    g_B, a2 = None, float("nan")
    a3 = float("nan")
    eps_sum = np.zeros(n_exons)
    eps_count = np.zeros(n_exons)
    weight_sum = np.zeros(n_exons)
    for _ in range(2):  # one refresh of mu_g from the boundary law; measured stable
        r_b = (rho_adj / max(rho_off, 1e-30)) * max(mu_r - 1.0, 1e-9) / max(mu_g - 1.0, 1e-9)
        a_b = 1.0 / (1.0 + r_b)
        weighted = a_b * nb_pair
        num = np.bincount(cls, weights=weighted, minlength=2)
        den = np.bincount(cls, weights=nb_pair, minlength=2)
        # SIGNED enrichment ratio (1 = uniform); the clip lives at the exon mean so noise
        # cancels instead of accumulating one-sidedly. An exon has at most two flanking boundaries, so
        # the per-exon mean below is a sum of at most two terms — the same arithmetic as the list it
        # replaces, in the same order.
        eps_pair = weighted / max(rho_off * max(mu_g - 1.0, 1e-9), 1e-30)
        eps_sum = np.bincount(exon_region, weights=eps_pair, minlength=n_exons)
        eps_count = np.bincount(exon_region, minlength=n_exons).astype(np.float64)
        weight_sum = np.bincount(exon_region, weights=weighted, minlength=n_exons)
        if den[0] <= 0.0 or den[1] <= 0.0:
            break
        a2, a3 = float(num[0]) / float(den[0]), float(num[1]) / float(den[1])
        n2 = float(obs_pools[_GDNA_CROSSING_POOLS[0]].sum())
        n3 = float(obs_pools[_GDNA_CROSSING_POOLS[1]].sum())
        total_opp = np.asarray(gdna_opportunity.total, dtype=np.float64)
        f = []
        for pool in _GDNA_CROSSING_POOLS:
            a_of_l = np.asarray(gdna_opportunity.pools[pool], dtype=np.float64)
            prob = np.zeros_like(a_of_l)
            np.divide(a_of_l, total_opp, out=prob, where=total_opp > 0.0)
            f.append(_normalized(detilt_pool(obs_pools[pool], prob)))
        f2, f3 = f[0][: g_C.size], f[1][: g_C.size]
        sep = a2 - a3
        f_mix = _normalized(n2 * f2 + n3 * f3)
        if sep == 0.0:
            g_B = f_mix
        else:
            # the same fade as the contained pair: resolve the separation against its own error and
            # blend the complement-first inversion toward the pools' mixture by that weight.
            se_b = float(np.sqrt(a2 * a2 / max(n2, 1.0) + a3 * a3 / max(n3, 1.0)))
            lam_b = _resolution_weight(sep * sep, se_b * se_b)
            a_mix = (a2 * n2 + a3 * n3) / max(n2 + n3, 1.0)
            r_hat = _normalized(np.clip((a3 * f2 - a2 * f3) / (a3 - a2), 0.0, None))
            inverted = _normalized(np.clip(f_mix - (1.0 - a_mix) * r_hat, 0.0, None))
            g_B = _normalized(np.clip(f_mix + lam_b * (inverted - f_mix), 0.0, None))
        mu_next = float((g_B * L_axis[: g_B.size]).sum())
        if not np.isfinite(mu_next) or abs(mu_next - mu_g) < 0.25:
            break
        mu_g = mu_next
    m_B = 0.0 if g_B is None else a2 * n2 + a3 * n3

    # ── the on-target excess: exon classes no pool samples, at (eps - 1)+ only
    is_ex = (ty == exon) & (ell > 0.0)
    h_E = np.zeros_like(g_C)
    m_E = 0.0
    # ⛔ An exon no boundary pair witnessed contributes NOTHING, and that is arithmetic rather than a
    # choice: its weight is `_resolution_weight(signal, inf) = 0`, so its excess is 0 and it is skipped.
    witnessed = is_ex & (eps_count > 0.0)
    if witnessed.any():
        for e_idx in np.flatnonzero(witnessed):
            eps_e = float(eps_sum[e_idx] / eps_count[e_idx])
            n_e = float(weight_sum[e_idx])
            # The excess has its OWN resolution weight, and it is not `lam`: `lam` asks whether the
            # two strata's LAWS differ, while this asks whether this exon's ENRICHMENT differs from
            # 1, and gating one on the other suppresses a real correction. Same helper, its own
            # signal and its own noise: a ratio estimated from `n_e` gDNA crossings has relative
            # variance ~ 1/n_e, so `Var(eps) ~ eps^2/n_e`.
            w_e = _resolution_weight(
                (eps_e - 1.0) ** 2, (eps_e * eps_e / n_e) if n_e > 0.0 else np.inf
            )
            excess_e = w_e * max(eps_e - 1.0, 0.0)
            if excess_e <= 0.0:
                continue
            le = ell[e_idx]
            shape = g_C * (
                np.clip(le - L_axis + 1.0, 0.0, None) + np.clip(L_axis - le - 1.0, 0.0, None)
            )
            h_E += rho_off * excess_e * shape
            m_E += rho_off * excess_e * float(shape.sum())

    if m_C + m_B <= 0.0:
        return None, None, GdnaRealized(False, "no gDNA mass anywhere", 0.0, 0.0, a2, a3)

    # Couple the two estimands. `lam` is how well the capture-induced split between the strata is
    # resolved; at `lam = 0` the two laws are the SAME ARRAY, which is the honest answer whenever one
    # stratum is starved — at zero capture the boundaries are empty, at very strong capture the
    # contained pools are, and in neither case does the data distinguish a chemistry law from a census.
    if g_B is None:
        uniform, realized, lam = g_C.copy(), g_C.copy(), 0.0
    else:
        uniform, realized, lam = _couple_estimands(g_C, m_C, g_B, m_B)

    m0 = min(uniform.size, realized.size, h_E.size)
    # the excess rides its OWN resolution weight (applied per exon above), not `lam`
    realized = realized[:m0] + h_E[:m0] / max(m_C + m_B, 1e-30)
    if not realized.sum() > 0.0 or not uniform[:m0].sum() > 0.0:
        return None, None, GdnaRealized(False, "empty census", 0.0, 0.0, a2, a3)
    scale = float(obs_pools[list(_GDNA_POOLS)].sum())
    total_mass = m_C + m_B + m_E
    return (
        _normalized(realized) * scale,
        _normalized(uniform[:m0]) * scale,
        GdnaRealized(True, "", m_B / max(m_C + m_B, 1e-30), m_E / total_mass, a2, a3),
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
    region_lengths: np.ndarray | None = None,
    region_types: np.ndarray | None = None,
    prior_ess: float = POOL_EB_PRIOR_ESS,
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

    ``region_lengths`` and ``region_types`` are what enable the two-pool CONTRAST: the module
    docstring derives it, and :func:`_deconvolved_gdna_counts` implements it and owns every reason
    it declines. Both come straight from
    :func:`~rigel.calibration.splice_graph.build_region_partition_arrays`, the same partition the
    scanner deposits into, so they cannot disagree with the banks they index.

    Omitting them is a supported state, not a degraded one — the second pass and most tests do —
    and ``None`` then means no annotation was offered. The fallback is the honest one rather than
    the convenient one: the RNA pool stays tilted and the gDNA pool falls back to the CONTAINED pair
    alone (:func:`gdna_contained_fl_mass`). It does not fall back to the four pools pooled raw,
    which is worse than either.

    For the EB kernel over three free histograms — the shape a unit test needs and production never
    has — see :func:`_fl_models_from_histograms`.
    """
    rna_counts = rna_fl_mass(payload)
    # One de-tilt implementation, shared: it preserves the pool TOTAL (the EB shrinkage reads that as
    # "how much evidence stands behind this shape") and drops bins the opportunity says are impossible.
    if sj_opportunity is not None:
        rna_counts = detilt_pool(rna_counts, sj_opportunity)

    contrast = None
    realized_counts, realized = None, None
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
            # the SECOND estimand: the library-census law for the scorer. It reads the uniform-frame
            # result and the same banks; on decline the two estimands coincide, which is the honest
            # off-capture answer rather than a degraded one.
            realized_counts, coupled_uniform, realized = _realized_gdna_counts(
                payload,
                gdna_opportunity,
                region_lengths,
                region_types,
                rna_counts,
                gdna_counts,
            )
            # the coupling can move the UNIFORM law too, and that is the point: when the contained
            # stratum is starved the chemistry law is not estimable and must borrow the one that is.
            if coupled_uniform is not None:
                gdna_counts = coupled_uniform

    return _fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=rna_counts,
        gdna_counts=gdna_counts,
        max_size=int(payload.max_length),
        prior_ess=prior_ess,
        pool_counts=payload.pool_lengths,
        gdna_contrast=contrast,
        gdna_realized_counts=realized_counts,
        gdna_realized=realized,
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
    gdna_realized_counts: np.ndarray | None = None,
    gdna_realized: "GdnaRealized | None" = None,
) -> FLModels:
    """The smooth-EB kernel: three histograms in, three pmfs out.

    Not a production entry point. Production has exactly one source for all three histograms and
    reaches it through :func:`build_fl_models`; this exists so the shrinkage policy itself can be
    exercised over anchors and pools that no real payload would produce.
    """
    global_aligned = _aligned(global_counts, max_size)
    rna_aligned = _aligned(rna_counts, max_size)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, n_rna = _smooth_eb(rna_aligned, global_pmf, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)
    # the realized law is shrunk exactly like its sibling; with no realized estimate the two
    # estimands COINCIDE — same array values, so an off-capture or input-starved build behaves
    # byte-identically to a build with a single law.
    if gdna_realized_counts is not None:
        gdna_realized_pmf, _ = _smooth_eb(
            _aligned(gdna_realized_counts, max_size), global_pmf, prior_ess
        )
    else:
        gdna_realized_pmf = gdna_pmf.copy()
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
        gdna_realized_pmf=gdna_realized_pmf,
        gdna_realized=gdna_realized,
    )
