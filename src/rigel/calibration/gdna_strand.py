"""Strand Beta-Binomial overdispersion fits for both components (gDNA and RNA).

Both gDNA and RNA strand splits are modelled as **Beta-Binomial**:
``K⁺ | N ~ BetaBinom(N, mean, overdispersion)`` with the intra-class correlation in ``[0, 1)``
(``0`` ⇒ Binomial). The **mean** differs by component — gDNA is unstranded (biologically fixed at
½), RNA is skewed (``rna_sense_frac`` = κ from the spliced channel) — but the **overdispersion**
machinery is shared: a per-region/per-boundary sense split is more spread out than Binomial
(local capture bias, PCR, pair-anchoring noise), and that excess variance is the overdispersion.

**Symmetric by design.** Earlier the RNA side was forced Binomial while gDNA was Beta-Binomial;
that asymmetry made the strand likelihood *spuriously informative on unstranded data* (the
``−½·log var`` term prefers the lower-variance component, pulling balanced nodes toward RNA). Both
components now carry a fitted overdispersion with the **same default prior**, so under sparse data
the two collapse to the same distribution and an unstranded node is uninformative — as it must be.
The gDNA mean ½ and RNA mean κ are unchanged; only RNA gains the overdispersion term. The RNA
overdispersion is fit from **boundary-side spliced counts** (spliced ⇒ pure RNA); see
:func:`fit_rna_strand_overdispersion` and the substrate wrapper in :mod:`calibrate`.

**Breaking the circularity** (gDNA only). Fitting the overdispersion needs to know which fragments are gDNA,
which is what the deconvolution determines — circular. We break it with the count⊥strand
conditional independence the engine already relies on: the **count module** supplies a per-node
gDNA fraction ``gdna_weight`` (``count_gdna_frac``, a raw count/density ratio) that is *independent*
of the strand overdispersion (it uses no strand information at all). Given those weights and the RNA
sense rate, the overdispersion is identified from the **excess variance of the sense split beyond
Binomial**, attributable to the gDNA fragments.

**Estimator — pooled method of moments.** For each seed node ``s`` (a count-observable region or
boundary side — intergenic, intronic, exon–intron / exon–intergenic seam) with ``sense_s`` of
``n_s`` gDNA-eligible unspliced fragments and count-derived gDNA weight ``w_s``:

    mean_s        = ½·w_s + rna_sense_frac·(1 − w_s)        # mixture sense rate
    excess_var_s  = (sense_s − n_s·mean_s)²  −  n_s·mean_s·(1 − mean_s)   # observed − Binomial
    gdna_var_s    = (w_s·n_s)·(w_s·n_s − 1)·¼                # BetaBinom excess-variance scale

    od_mom = Σ_s excess_var_s / Σ_s gdna_var_s        # pooled point estimate

The point estimate is then **shrunk toward a prior overdispersion** ``od₀`` (a symmetric
``Beta(a, a)`` "floor"; default ``a = 3`` ⇒ ``od₀ = 1/7``) by seed-node count — sparse/low-signal
libraries lean on the prior, abundant ones on the fit — and clamped to ``[0, _MAX_OVERDISPERSION]``
(the ``Beta(2, 2)`` ceiling, ``od = 0.2``, the most overdispersion allowed). This continuous
shrinkage replaces the earlier hard min-seed-node / significance gates. The MoM is closed-form,
``O(n_seed_nodes)``, and uses the **same variance decomposition the decode applies**
(:mod:`strand_likelihood`), so fit and application are consistent. The prior parameters live on
``CalibrationConfig`` (``gdna_strand_prior_alpha_beta``, ``gdna_strand_prior_weight``).

The substrate→seed-node extraction wrapper lives in :mod:`calibrate` (it needs the region/
boundary geometry); this module is the pure estimator + model so it is trivially testable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .strand_deconv import boundary_side_seeds
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


@dataclass(frozen=True, slots=True)
class GdnaStrandModel:
    """Fitted gDNA strand model: the global Beta-Binomial overdispersion + fit provenance."""

    gdna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
    n_seed_nodes: int  # seed nodes that carried gDNA-eligible fragments
    n_seed_fragments: int  # total gDNA-eligible fragments across seed nodes
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
    :class:`strand_balance.StrandBalance`; this carries only the between-boundary overdispersion."""

    rna_strand_overdispersion: float  # intra-class correlation in [0, _MAX_OVERDISPERSION]
    n_seed_nodes: int  # boundary sides that carried spliced fragments
    n_seed_fragments: int  # total spliced fragments across seed sides (with junction double-count)
    fallback_used: bool  # True ⇒ no spliced strand signal ⇒ returned the prior overdispersion

    def beta_concentration(self) -> float:
        """Symmetric Beta concentration ``a`` of ``Beta(a, a)`` for the RNA sense rate (about κ)."""
        od = self.rna_strand_overdispersion
        return float("inf") if od <= 0.0 else 0.5 * (1.0 - od) / od


def _fit_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    node_mean: np.ndarray,
    component_frac: np.ndarray,
    component_mean: np.ndarray,
    *,
    prior_overdispersion: float,
    prior_weight: float,
) -> tuple[float, int, int, bool]:
    """Shared pooled-MoM + prior-shrinkage core for one component's strand overdispersion.

    Component-agnostic, with three per-node inputs:

    * ``node_mean`` — the node's mixture sense rate, used to subtract the Binomial variance:
      ½·w + κ·(1−w) for a gDNA seed (a gDNA/RNA mix), κ for a pure-RNA spliced seed.
    * ``component_frac`` — the fraction of the node's fragments in the component being fit (the
      gDNA weight ``w`` for gDNA, ``1`` for pure-RNA spliced); sets ``n_c = component_frac·N``.
    * ``component_mean`` — the component's *own* sense mean ``μ_c`` (½ for gDNA, κ for RNA); the
      BetaBinom excess variance of ``n_c`` correlated fragments scales as
      ``n_c·(n_c − 1)·μ_c·(1 − μ_c)`` — **not** ¼ unless ``μ_c = ½``.

    Returns ``(overdispersion, n_seed_nodes, n_seed_fragments, fallback_used)``; clamped to
    ``[0, _MAX_OVERDISPERSION]``; fallback (pooled denominator ≤ 0) ⇒ ``prior_overdispersion``.
    """
    sense = np.asarray(sense, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)
    node_mean = np.asarray(node_mean, dtype=np.float64)
    component_frac = np.asarray(component_frac, dtype=np.float64)
    component_mean = np.asarray(component_mean, dtype=np.float64)

    valid = total > 0.0
    binom_var = total * node_mean * (1.0 - node_mean)
    excess_var = (sense - total * node_mean) ** 2 - binom_var
    comp_frags = component_frac * total
    comp_var = component_mean * (1.0 - component_mean)  # μ_c(1−μ_c): ¼ for gDNA, κ(1−κ) for RNA
    # BetaBinom excess-variance scale n_c(n_c − 1)·μ_c(1−μ_c) (clipped ≥ 0 for tiny component mass).
    var_scale = np.maximum(comp_frags * (comp_frags - 1.0), 0.0) * comp_var

    num = float(np.sum(excess_var[valid]))
    denom = float(np.sum(var_scale[valid]))
    n_seed_nodes = int(np.sum(valid & (comp_frags > 0.0)))
    n_seed_fragments = int(np.sum(total[valid & (comp_frags > 0.0)]))
    prior_overdispersion = float(prior_overdispersion)
    prior_weight = max(float(prior_weight), 0.0)

    fallback = denom <= 0.0 or not np.isfinite(num)
    if fallback:
        # No component strand information in the seeds → fall back entirely to the prior.
        od = prior_overdispersion
    else:
        od_mom = num / denom
        # Precision-weighted shrinkage toward the prior overdispersion, by seed-node count.
        total_weight = n_seed_nodes + prior_weight
        od = (
            (n_seed_nodes * od_mom + prior_weight * prior_overdispersion) / total_weight
            if total_weight > 0.0
            else od_mom
        )
    od = float(np.clip(od, _OVERDISPERSION_FLOOR, _MAX_OVERDISPERSION))
    return od, n_seed_nodes, n_seed_fragments, fallback


def fit_gdna_strand_overdispersion(
    sense: np.ndarray,
    total: np.ndarray,
    gdna_weight: np.ndarray,
    rna_sense_frac: float,
    *,
    prior_overdispersion: float = 0.0,
    prior_weight: float = 0.0,
) -> GdnaStrandModel:
    """Pooled method-of-moments fit of the global gDNA strand overdispersion, with prior shrinkage.

    The pooled MoM point estimate is shrunk toward ``prior_overdispersion`` by seed-node count:
    ``od = (n_seed_nodes·od_mom + prior_weight·prior_overdispersion) / (n_seed_nodes + prior_weight)``.
    Sparse seed sets (few informative nodes) lean on the prior; abundant sets on the fit. This
    replaces the earlier hard min-node / significance gates (no thresholds — graceful with data).
    The result is clamped to ``[0, _MAX_OVERDISPERSION]`` (the ``Beta(2, 2)`` ceiling).

    Parameters
    ----------
    sense, total : np.ndarray
        Per-seed-node sense count ``K⁺`` and total gDNA-eligible unspliced count ``N``.
    gdna_weight : np.ndarray
        Per-seed-node count-clue gDNA fraction ``∈ [0, 1]`` (independent of the overdispersion).
    rna_sense_frac : float
        Library RNA sense fraction ``κ`` (the spliced-channel ``StrandModel`` mean).
    prior_overdispersion : float
        Prior overdispersion to shrink toward (the ``Beta(a, a)`` "floor"; see
        :func:`overdispersion_for_beta`). ``0`` with ``prior_weight = 0`` ⇒ pure MoM.
    prior_weight : float
        Prior strength in effective seed-node units. ``0`` ⇒ no shrinkage.

    Returns
    -------
    GdnaStrandModel
        ``gdna_strand_overdispersion`` in ``[0, _MAX_OVERDISPERSION]``; ``fallback_used`` when the
        seeds carry no gDNA strand signal (pooled denominator ≤ 0) ⇒ the prior is returned.
    """
    weight = np.clip(np.asarray(gdna_weight, dtype=np.float64), 0.0, 1.0)
    # gDNA component: the node's sense mean is the mixture ½·w + κ·(1−w); the component whose
    # shared rate inflates the variance is the gDNA fraction w, with its own mean ½ (⇒ scale ¼).
    node_mean = 0.5 * weight + float(rna_sense_frac) * (1.0 - weight)
    component_mean = np.full(node_mean.shape, 0.5, dtype=np.float64)
    od, n_seed_nodes, n_seed_fragments, fallback = _fit_overdispersion(
        sense,
        total,
        node_mean,
        weight,
        component_mean,
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )
    return GdnaStrandModel(
        gdna_strand_overdispersion=od,
        n_seed_nodes=n_seed_nodes,
        n_seed_fragments=n_seed_fragments,
        fallback_used=fallback,
    )


def _region_seeds(substrate, region_arrays, node_density):
    """``(sense, total, gdna_weight)`` from the count-observable contained regions.

    Intergenic (``TS_NONE``) and intron-only (``TS_POS``/``TS_NEG``) regions — i.e.
    ``node_density.region_count_observable`` — excluding ``TS_AMBIG`` (both strands, no defined
    sense). The weight is the count-clue gDNA fraction ``node_density.count_gdna_frac`` (=
    ``clip(density·eff_len / mass)``, density cleaned by the strand *mean* ½, not the dispersion).
    It reads the explicit count-prior MEAN (``count_gdna_frac``) directly — decoupled from the
    count-prior concentration, which carries the overdispersion-honest precision.
    """
    ts = np.asarray(region_arrays.strand_class)
    contained = substrate.contained
    pos = np.asarray(contained.n_unspliced_pos, dtype=np.float64)
    neg = np.asarray(contained.n_unspliced_neg, dtype=np.float64)
    total = pos + neg
    sense = np.where(ts == TS_NEG, neg, pos)  # orient to transcript sense (arbitrary for TS_NONE)
    weight = np.clip(np.asarray(node_density.count_gdna_frac, dtype=np.float64), 0.0, 1.0)
    seed = np.asarray(node_density.region_count_observable) & (ts != TS_AMBIG)
    return sense[seed], total[seed], weight[seed]


def fit_gdna_strand_from_substrate(
    substrate,
    region_arrays,
    node_density,
    boundary_side_eff_len,
    *,
    rna_sense_frac: float,
    prior_overdispersion: float = 0.0,
    prior_weight: float = 0.0,
) -> GdnaStrandModel:
    """Fit the global gDNA strand overdispersion from the calibration substrate.

    Pools two kinds of count-observable seed (the same seeds the density estimator trusts):

    * **contained regions** — intergenic + intron-only (:func:`_region_seeds`);
    * **boundary sides** — exon–intron / exon–intergenic seams
      (:func:`strand_deconv.boundary_side_seeds`), needed under hybrid capture, which depletes
      off-target intergenic/intronic gDNA.

    Both contribute ``(sense, total, gdna_weight)`` on the same footing — the weight being the
    overdispersion-free count-clue gDNA fraction — and the pooled estimator fits one global
    overdispersion, shrunk toward ``prior_overdispersion`` (strength ``prior_weight``).
    """
    r_sense, r_total, r_weight = _region_seeds(substrate, region_arrays, node_density)
    b_sense, b_total, b_weight = boundary_side_seeds(
        substrate, region_arrays, node_density, boundary_side_eff_len
    )
    sense = np.concatenate([r_sense, b_sense])
    total = np.concatenate([r_total, b_total])
    weight = np.concatenate([r_weight, b_weight])
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
    prior_overdispersion: float = 0.0,
    prior_weight: float = 0.0,
) -> RnaStrandModel:
    """Pooled method-of-moments fit of the global RNA strand overdispersion, with prior shrinkage.

    The twin of :func:`fit_gdna_strand_overdispersion`, on **pure-RNA spliced** nodes: spliced
    fragments are pure RNA, so each seed node is all-RNA — ``node_mean = rna_sense_frac`` (κ) and
    the component fraction is ``1`` (the whole node is the RNA component). The excess variance of
    the sense split beyond ``Binomial(N, κ)`` identifies the overdispersion. Shrinks toward
    ``prior_overdispersion`` by seed-node count and clamps to ``[0, _MAX_OVERDISPERSION]`` exactly
    like the gDNA fit, so under sparse data the two collapse to the same prior.

    Parameters
    ----------
    sense, total : np.ndarray
        Per-seed-node motif-relative spliced sense count and total spliced count ``N``.
    rna_sense_frac : float
        Library RNA sense fraction ``κ`` (the spliced-channel ``StrandModel`` mean) — the node mean.
    """
    total = np.asarray(total, dtype=np.float64)
    # Pure-RNA spliced node: mixture mean AND component mean are both κ; component fraction is 1.
    node_mean = np.full(total.shape, float(rna_sense_frac), dtype=np.float64)
    component_frac = np.ones(total.shape, dtype=np.float64)  # pure RNA
    od, n_seed_nodes, n_seed_fragments, fallback = _fit_overdispersion(
        sense,
        total,
        node_mean,
        component_frac,
        node_mean,  # component_mean = κ ⇒ scale κ(1−κ)
        prior_overdispersion=prior_overdispersion,
        prior_weight=prior_weight,
    )
    return RnaStrandModel(
        rna_strand_overdispersion=od,
        n_seed_nodes=n_seed_nodes,
        n_seed_fragments=n_seed_fragments,
        fallback_used=fallback,
    )


def fit_rna_strand_from_substrate(
    substrate,
    *,
    rna_sense_frac: float,
    prior_overdispersion: float = 0.0,
    prior_weight: float = 0.0,
) -> RnaStrandModel:
    """Fit the global RNA strand overdispersion from boundary-side spliced counts.

    Spliced fragments cross splice junctions, so they live on the **boundary sides** (the
    accumulator deposits a spliced jump as boundary flux), not in the region-contained view — and
    spliced ⇒ pure RNA. Each boundary side's motif-relative ``(n_spliced_sense, n_spliced_total)``
    is one seed node; the ``left`` and ``right`` sides are pooled. Sides with no spliced fragments
    drop out (the estimator filters ``total > 0``); orientation is motif-relative, so AMBIG regions
    are fine (no count-observable filter, unlike the gDNA fit).

    **Known approximation (v1):** a fragment spanning K junctions credits ~K boundary sides
    (``accumulator.cpp``), so multi-junction fragments are double-counted across sides. This
    inflates the effective count and mildly *under*-estimates the overdispersion (correlated
    repeats look more Binomial); accepted for v1, tracked as future work (count one side per
    boundary) — see ``docs/em_strand/05_rna_strand_bb_plan.md`` §9.
    """
    sense = np.concatenate(
        [
            np.asarray(substrate.left.n_spliced_sense, dtype=np.float64),
            np.asarray(substrate.right.n_spliced_sense, dtype=np.float64),
        ]
    )
    antisense = np.concatenate(
        [
            np.asarray(substrate.left.n_spliced_antisense, dtype=np.float64),
            np.asarray(substrate.right.n_spliced_antisense, dtype=np.float64),
        ]
    )
    total = sense + antisense
    return fit_rna_strand_overdispersion(
        sense,
        total,
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
    "fit_rna_strand_from_substrate",
    "overdispersion_for_beta",
]
