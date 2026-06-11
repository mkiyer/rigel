"""Per-node gDNA/RNA deconvolution by precision-weighted strand→count deference (decoupled).

Each **node** (every region and every boundary side, kept separate) is deconvolved into gDNA / RNA
mass by a **precision-weighted blend** of two independent estimators — never their product:

    g = w · g_strand  +  (1 − w) · g_count,     w = (2κ − 1)²   (κ = rna_sense_frac)
      g_strand : Beta-Binomial posterior median over the gDNA fraction (weak Beta(½,½) prior)
      g_count  : the count module's fraction (count_gdna_frac = density·eff_len / mass, raw)

The weight ``w = (2κ−1)²`` is the strand channel's **discriminability** — the per-fragment strand
Fisher information / the squared standardized separation between gDNA (sense ½) and RNA (sense κ). It
``→ 1`` at high strand specificity (trust the strand estimate) and ``→ 0`` at unstranded (defer fully
to count): a smooth deference with **no gate**. Because ``w`` is an *effect size* (a function of κ,
not of the spliced-read count), a near-unstranded library gets ``w ≈ 0`` regardless of depth — which
no significance threshold could deliver. A node with no defined sense (AMBIG / no shared strand) is
count-only (``w`` applies only to strand-observable nodes).

Why a blend, not their product: the strand estimator is **unbiased** (gDNA symmetric at ½, RNA at κ)
but noisy at low N/SS; the count estimator is **biased** under hybrid capture. We weight the unbiased
channel by its *signal strength* and lean on the biased channel only for what strand cannot resolve —
*not* naive inverse-variance weighting (which assumes both unbiased and penalises strand for
overdispersion variance that averages out at locus aggregation). See
``docs/calibration/count_channel_capture_design.md`` (the redesign plan, Phase 1) and
``decoupled_calibration_design.md``; the retired joint product is archived in
``archive/joint_deconvolution.md``.

The **strand module** is the Beta-Binomial posterior (the exact, clip-free robust strand
deconvolution): ``log_post(g) = log Beta(g; ½, ½) + strand_loglik(sense, antisense | g, κ, od)`` on a
grid, summarized by its median (point) and variance (width). The blended point estimate is then read
at the **FP-rate quantile** ``g(q) = clip(center + Φ⁻¹(q)·σ)`` (Phase 2; ``deconv_quantile`` default
0.5 ⇒ ``Φ⁻¹=0`` ⇒ bit-identical no-op). Amounts use the conserved **fractional mass** ``M``: gDNA =
``g·M``, RNA = ``(1−g)·M`` + the deterministic spliced mass. Regions and boundaries are combined into
loci only **after** calibration (``assemble_priors``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtri

from .run_fill import same_ref_left_right
from .signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from .strand_likelihood import strand_loglik

_GRID_EPS = 1.0e-4
_STRAND_PRIOR = 0.5  # Jeffreys Beta(½,½): the standard uninformative prior for the strand posterior


@dataclass(frozen=True, slots=True)
class NodeDeconv:
    """Per-node deconvolution result (regions or boundary sides; kept separate)."""

    gdna_mass: np.ndarray  # float64[K]
    rna_mass: np.ndarray  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac: np.ndarray  # float64[K] — reported gDNA fraction of the UNSPLICED mass
    gdna_frac_var: np.ndarray  # float64[K] — posterior variance (0 for count-routed nodes)


def _deconv_per_node(
    mass_unspl,
    mass_spliced,
    sense,
    antisense,
    strand_observable,
    count_gdna_frac,
    count_gdna_frac_var,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion,
    rna_strand_overdispersion,
    deconv_quantile,
    n_grid,
) -> NodeDeconv:
    """Blend the strand posterior and the count fraction per node, then read a quantile, then split mass.

    The **point estimate** (``center``) is the precision-weighted deference ``w·g_strand + (1−w)·g_count``
    with ``w = (2κ−1)²`` — the strand channel's discriminability (→1 at high SS, →0 at unstranded; a
    smooth, bias-robust deference with no gate). ``g_strand`` is the Beta-Binomial posterior median; a
    node with no defined sense (or κ=½) takes the count fraction.

    The **FP-rate quantile knob** (Phase 2) then reads ``g(q) = clip(center + Φ⁻¹(q)·σ)`` where ``σ`` is
    the combined per-node std: ``√(w²·σ²_strand + (1−w)²·σ²_count)`` (strand-observable) or ``σ_count``
    (count-routed). ``q=0.5 ⇒ Φ⁻¹=0 ⇒ g=center`` — a bit-identical no-op default. ``q>0.5`` is FP-averse
    (more gDNA); the shift is **uncertainty-aware** (wider posterior ⇒ larger shift). The count variance
    is used only to *widen* the quantile (safe), never to *sharpen* the blend (kept bias-robust at
    ``(2κ−1)²`` — the count σ is anti-calibrated under capture; see ``docs/calibration/phase2_design.md``).
    """
    k = mass_unspl.shape[0]
    grid = np.linspace(_GRID_EPS, 1.0 - _GRID_EPS, n_grid)
    log_grid = np.log(grid)
    log_1mgrid = np.log1p(-grid)
    log_prior = (_STRAND_PRIOR - 1.0) * log_grid + (_STRAND_PRIOR - 1.0) * log_1mgrid
    w_strand = (2.0 * float(rna_sense_frac) - 1.0) ** 2  # strand discriminability, in [0, 1]
    z = float(ndtri(min(max(float(deconv_quantile), 1e-6), 1.0 - 1e-6)))  # Φ⁻¹(q); 0 at q=0.5
    gdna = np.zeros(k)
    rna = np.zeros(k)
    gdna_frac = np.zeros(k)
    gdna_frac_var = np.zeros(k)
    for i in range(k):
        m = float(mass_unspl[i])
        if m <= 0.0:
            rna[i] = float(mass_spliced[i])  # only deterministic spliced RNA, if any
            continue
        g_count = min(max(float(count_gdna_frac[i]), 0.0), 1.0)
        var_count = max(float(count_gdna_frac_var[i]), 0.0)
        if w_strand > 0.0 and strand_observable[i] and (sense[i] + antisense[i]) > 0:
            # STRAND module: Beta-Binomial posterior (weak prior, no count term).
            log_post = log_prior + strand_loglik(
                grid,
                sense[i],
                antisense[i],
                rna_sense_frac,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
            )
            post = np.exp(log_post - log_post.max())
            post /= post.sum()
            mean = float(np.dot(grid, post))
            var_strand = float(np.dot((grid - mean) ** 2, post))
            g_strand = float(np.clip(np.interp(0.5, np.cumsum(post), grid), 0.0, 1.0))
            center = w_strand * g_strand + (1.0 - w_strand) * g_count
            var = w_strand**2 * var_strand + (1.0 - w_strand) ** 2 * var_count
        else:
            center = g_count  # no usable strand signal (unstranded / no sense) ⇒ count only
            var = var_count
        gdna_frac_var[i] = var  # combined per-node posterior variance (Phase 1/2)
        frac = center if z == 0.0 else min(max(center + z * var**0.5, 0.0), 1.0)
        gdna_frac[i] = frac
        gdna[i] = frac * m
        rna[i] = (1.0 - frac) * m + float(mass_spliced[i])
    return NodeDeconv(
        gdna_mass=gdna, rna_mass=rna, gdna_frac=gdna_frac, gdna_frac_var=gdna_frac_var
    )


def deconv_regions(
    substrate,
    region_arrays,
    node_density,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    deconv_quantile=0.5,
    n_grid=200,
) -> NodeDeconv:
    """Deconvolve each region's contained mass (a node) into gDNA / RNA by the strand/count blend.

    A region is strand-observable iff its transcript strand is defined (``TS_POS`` / ``TS_NEG``);
    ``TS_NONE`` (intergenic) and ``TS_AMBIG`` are count-only. The count fraction + its variance are the
    precomputed ``node_density.count_gdna_frac`` / ``count_gdna_frac_var``; the weight ``w=(2κ−1)²`` and
    the FP-rate quantile ``deconv_quantile`` are applied in ``_deconv_per_node``.
    """
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    sense = np.where(ts == TS_NEG, neg, pos)  # orient to transcript sense
    antisense = (pos + neg) - sense
    strand_observable = (ts == TS_POS) | (ts == TS_NEG)
    return _deconv_per_node(
        c.mass_unspliced,
        c.mass_spliced,
        sense,
        antisense,
        strand_observable,
        node_density.count_gdna_frac,
        node_density.count_gdna_frac_var,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        deconv_quantile=deconv_quantile,
        n_grid=n_grid,
    )


@dataclass(frozen=True, slots=True)
class _SideQuantities:
    """Per-region boundary-side quantities shared by the deconvolution and the strand-fit seeds."""

    sense: np.ndarray
    antisense: np.ndarray
    n_side: np.ndarray  # pos + neg (count evidence)
    mass: np.ndarray
    mass_spliced: np.ndarray
    count_gdna_frac: np.ndarray  # count module's gDNA fraction for this side (raw density ratio)
    count_gdna_frac_var: np.ndarray  # Phase-1 Poisson floor σ_g² = g²/n_side (cap g(1−g))
    strand_observable: np.ndarray
    count_observable: np.ndarray


def _left_right_neighbors(ts, ref_id, boundary_count_observable):
    """Per-region neighbour strand class + same-ref + boundary count-observability for both sides.

    Region ``r``'s LEFT side is the right side of boundary ``(r−1, r)`` (neighbour ``r−1``); its
    RIGHT side is the left side of boundary ``(r, r+1)`` (neighbour ``r+1``). ``boundary_count_observable[r]``
    describes boundary ``(r, r+1)``. Returns
    ``(left_same, ts_prev, left_observable, right_same, ts_next, right_observable)``.
    """
    r = ts.shape[0]
    left_same, right_same = same_ref_left_right(ref_id)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_observable = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        ts_prev[1:] = ts[:-1]
        left_observable[1:] = boundary_count_observable[:-1]
        ts_next[:-1] = ts[1:]
    return left_same, ts_prev, left_observable, right_same, ts_next, boundary_count_observable


def _compute_side(
    view, same, ts_self, ts_other, side_count_observable, eff, region_density
) -> _SideQuantities:
    """Per-side sense split, count-prior density, and strand-observability.

    A side is **strand-observable** iff its boundary's two regions define a single, consistent
    transcript sense (so 'sense' is defined). An **intergenic (TS_NONE) region is a strand wildcard**
    — it carries no transcript, so it is compatible with either strand: a gene-edge boundary
    (stranded exon ↔ intergenic) is therefore oriented by the gene side, ``{POS,NONE}→POS``,
    ``{NEG,NONE}→NEG`` (same as ``{POS,POS}`` / ``{NEG,NEG}``). Only a genuine conflict — opposite
    strands ``{POS,NEG}`` or a ``TS_AMBIG`` flank (overlapping +/− transcripts) — leaves the sense
    undefined (CALIBRATION_TODO Issue #3; treating gene-edge sides as strand-blind leaked ~half their
    pure-gDNA crossing mass). The **count** fraction is the raw crossing density ratio: a
    count-observable side (no shared exon) reads ``count_gdna_frac → 1`` from its own crossing mass;
    otherwise it borrows the swept region density. (No strand cleaning — the count module is raw.)

    The **count variance** is the Phase-1 Poisson floor (a boundary side is a *direct* measurement, not
    an imputed region — no variance~mean LOESS): the crossing count ``n_side`` is Poisson, so the
    fraction's relative variance is ``1/n_side``, giving ``σ_g² = g_count²·(1/n_side)`` capped at the
    Bernoulli maximum ``g_count(1−g_count)`` — the same form as an observable region's own-count floor
    (``density_model._count_fraction_variance``). It feeds the FP-rate quantile (Phase 2) only.
    """
    mass = view.mass_unspliced
    pos = view.n_unspliced_pos.astype(np.float64)
    neg = view.n_unspliced_neg.astype(np.float64)
    # Strand sense with TS_NONE (intergenic) as a wildcard: POS-sense if each flank is POS-or-NONE
    # and at least one is POS (mirror for NEG); a TS_AMBIG flank or an opposite-strand pairing leaves
    # it undefined. Reduces to {POS,POS}/{NEG,NEG} when neither flank is intergenic.
    either_ambig = (ts_self == TS_AMBIG) | (ts_other == TS_AMBIG)
    self_pos_or_none = (ts_self == TS_POS) | (ts_self == TS_NONE)
    other_pos_or_none = (ts_other == TS_POS) | (ts_other == TS_NONE)
    self_neg_or_none = (ts_self == TS_NEG) | (ts_self == TS_NONE)
    other_neg_or_none = (ts_other == TS_NEG) | (ts_other == TS_NONE)
    cons_pos = (
        same & ~either_ambig & self_pos_or_none & other_pos_or_none
        & ((ts_self == TS_POS) | (ts_other == TS_POS))
    )
    cons_neg = (
        same & ~either_ambig & self_neg_or_none & other_neg_or_none
        & ((ts_self == TS_NEG) | (ts_other == TS_NEG))
    )
    strand_observable = cons_pos | cons_neg
    sense = np.where(cons_neg, neg, pos)
    n_side = pos + neg
    antisense = n_side - sense
    with np.errstate(divide="ignore", invalid="ignore"):
        # Raw crossing density (no strand cleaning); count-observable sides use their own crossing
        # density, others borrow the swept region density.
        own = np.where((mass > 0.0) & (eff > 0.0), mass / np.maximum(eff, 1e-12), 0.0)
        density = np.where(side_count_observable, own, region_density)
        count_gdna_frac = np.clip(
            np.where(mass > 0.0, density * eff / np.maximum(mass, 1e-12), 0.0), 0.0, 1.0
        )
        # Phase-1 Poisson floor: crossing count n_side ~ Poisson ⇒ relative variance 1/n_side;
        # σ_g² = g²·(1/n_side), capped at the Bernoulli maximum g(1−g) (sides with no crossings → 0).
        v_rel = np.where(n_side > 0.0, 1.0 / np.maximum(n_side, 1e-12), 0.0)
        count_gdna_frac_var = np.minimum(
            count_gdna_frac**2 * v_rel, count_gdna_frac * (1.0 - count_gdna_frac)
        )
    return _SideQuantities(
        sense=sense,
        antisense=antisense,
        n_side=n_side,
        mass=mass,
        mass_spliced=view.mass_spliced,
        count_gdna_frac=count_gdna_frac,
        count_gdna_frac_var=count_gdna_frac_var,
        strand_observable=strand_observable,
        count_observable=np.asarray(side_count_observable, dtype=bool),
    )


def deconv_sides(
    substrate,
    region_arrays,
    node_density,
    boundary_side_eff_len,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    deconv_quantile=0.5,
    n_grid=200,
) -> tuple[NodeDeconv, NodeDeconv]:
    """Deconvolve each boundary **side** as an independent node by the strand/count blend.

    Region ``r`` owns the **right** side of its left boundary (``substrate.left[r]``) and the
    **left** side of its right boundary (``substrate.right[r]``); both lie inside region ``r`` and so
    use ``boundary_side_eff_len[r]``. Returns ``(left, right)`` per-region ``NodeDeconv`` results
    (zero where a side carries no mass — e.g. reference edges).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )

    def _deconv(view, same, ts_other, side_obs):
        sq = _compute_side(view, same, ts, ts_other, side_obs, eff, region_density)
        return _deconv_per_node(
            sq.mass,
            sq.mass_spliced,
            sq.sense,
            sq.antisense,
            sq.strand_observable,
            sq.count_gdna_frac,
            sq.count_gdna_frac_var,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            deconv_quantile=deconv_quantile,
            n_grid=n_grid,
        )

    left = _deconv(substrate.left, l_same, ts_prev, l_obs)
    right = _deconv(substrate.right, r_same, ts_next, r_obs)
    return left, right


def boundary_side_seeds(substrate, region_arrays, node_density, boundary_side_eff_len):
    """``(sense, total, gdna_weight)`` seed arrays from count- & strand-observable boundary sides.

    The exon–intron / exon–intergenic seam seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the contained-region seeds (needed under hybrid capture,
    which depletes off-target intergenic/intronic gDNA). The weight is the count-clue gDNA fraction
    ``count_gdna_frac`` (≈ 1 for a count-observable side, gDNA by signature) — the same quantity the
    contained-region seeds use (``node_density.count_gdna_frac``).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )
    sides = (
        _compute_side(substrate.left, l_same, ts, ts_prev, l_obs, eff, region_density),
        _compute_side(substrate.right, r_same, ts, ts_next, r_obs, eff, region_density),
    )
    senses, totals, weights = [], [], []
    for sq in sides:
        seed = sq.count_observable & sq.strand_observable & (sq.n_side > 0.0)
        senses.append(sq.sense[seed])
        totals.append(sq.n_side[seed])
        weights.append(sq.count_gdna_frac[seed])
    return np.concatenate(senses), np.concatenate(totals), np.concatenate(weights)


__all__ = ["NodeDeconv", "deconv_regions", "deconv_sides", "boundary_side_seeds"]
