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


@dataclass(frozen=True, slots=True)
class StrandSplit:
    """Standalone strand deconvolution of a node's UNSPLICED mass into gDNA vs RNA, + likelihood info.

    The strand module's per-node output (sequential redesign). All float64[R]:

    * ``gdna_frac`` — the gDNA fraction read at the FP-quantile where the node is strand-informative,
      else ``NaN`` (the strand makes no guess; the count module imputes those nodes from the field).
    * ``gdna_mass`` / ``rna_mass`` — ``mass_unspliced`` split by ``gdna_frac`` (``NaN`` where not split);
      conserved (``gdna_mass + rna_mass = mass_unspliced``) where split. Spliced mass is NOT included —
      it is deterministic RNA the count module folds in via the 3-term.
    * ``info`` — ``I = N·(2κ−1)²`` (N = unspliced count, κ = ``rna_sense_frac``): the **likelihood**
      precision the count module weights by; ``0`` where strand-uninformative (AMBIG / no sense) or
      unstranded (κ=½). NOT the posterior variance (which is the finite prior variance at κ=½).
    """

    gdna_frac: np.ndarray
    gdna_mass: np.ndarray
    rna_mass: np.ndarray
    info: np.ndarray


def _grid_posterior_quantile(post: np.ndarray, grid: np.ndarray, q: float = 0.5) -> np.ndarray:
    """Per-row posterior quantile — the batched ``np.interp(q, cumsum(post_row), grid)`` (median at q=½).

    Vectorises the quantile over the ``(n, n_grid)`` posterior: the CDF crossing of ``q`` on each row,
    linearly interpolated on ``grid`` with ``np.interp``'s exact arithmetic (``slope·(q−x_lo)+f_lo``)
    and its end-clamps (``q ≤ cdf[0] → grid[0]``; ``q ≥ cdf[-1] → grid[-1]``). ``q=0.5`` is the median.
    """
    cdf = np.cumsum(post, axis=1)
    ng = grid.shape[0]
    j = np.clip((cdf < q).sum(axis=1), 1, ng - 1)  # first index with cdf ≥ q (interval upper bound)
    rows = np.arange(post.shape[0])
    c_lo = cdf[rows, j - 1]
    c_hi = cdf[rows, j]
    g_lo = grid[j - 1]
    slope = (grid[j] - g_lo) / np.where(c_hi > c_lo, c_hi - c_lo, 1.0)
    out = slope * (q - c_lo) + g_lo
    out = np.where(cdf[:, 0] >= q, grid[0], out)  # q ≤ cdf[0] ⇒ grid[0]
    out = np.where(cdf[:, -1] <= q, grid[-1], out)  # q ≥ cdf[-1] ⇒ grid[-1]
    return out


def strand_posterior_gdna_frac(
    sense: np.ndarray,
    antisense: np.ndarray,
    rna_sense_frac: float,
    *,
    gdna_strand_overdispersion: float,
    rna_strand_overdispersion: float,
    n_grid: int,
    deconv_quantile: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Beta-Binomial **strand posterior** over the gDNA fraction, per node — ``(point, variance)``.

    The strand module's core: ``log_post(g) = log Beta(g; ½, ½) + strand_loglik(sense, antisense | g, κ,
    od)`` on the ``n_grid`` grid, read at the **FP-quantile** ``deconv_quantile`` (the point estimate;
    median at the ``0.5`` default), with the posterior **variance** alongside. ``sense``/``antisense`` are
    1-D arrays (length n, oriented to transcript sense); returns two length-n arrays. Shared by the
    per-region deconvolution (:func:`_deconv_per_node`, which uses the median + variance) and the
    standalone strand module (:func:`strand_deconvolve`, which reads its own quantile) so both share the
    *same* posterior machinery.
    """
    sense = np.asarray(sense, dtype=np.float64)
    antisense = np.asarray(antisense, dtype=np.float64)
    grid = np.linspace(_GRID_EPS, 1.0 - _GRID_EPS, n_grid)
    log_prior = (_STRAND_PRIOR - 1.0) * (np.log(grid) + np.log1p(-grid))
    log_post = log_prior[None, :] + strand_loglik(
        grid[None, :],
        sense[:, None],
        antisense[:, None],
        rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
    )  # (n, n_grid)
    post = np.exp(log_post - log_post.max(axis=1, keepdims=True))
    post /= post.sum(axis=1, keepdims=True)
    mean = post @ grid  # (n,)
    var = ((grid[None, :] - mean[:, None]) ** 2 * post).sum(axis=1)
    g_q = np.clip(_grid_posterior_quantile(post, grid, deconv_quantile), 0.0, 1.0)
    return g_q, var


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
    info_scale,
) -> NodeDeconv:
    """Blend the strand posterior and the count fraction per node, then read a quantile, then split mass.

    The **point estimate** (``center``) is the precision-weighted blend ``w·g_strand + (1−w)·g_count`` with
    a **per-node** weight ``w = I/(I+I₀)`` — the carried strand information ``I = N·(2κ−1)²`` (the
    per-fragment discriminability ``(2κ−1)²`` times the node's unspliced count ``N``), ``I₀ = info_scale``.
    This is the *gradient* of strand-trustworthiness (not a cliff): ``w→1`` at high information (confident
    strand → its own deconvolution governs, capture-invariant), ``w→0`` at ``κ=½`` or thin ``N`` (the count
    imputation governs). ``g_strand`` (the region's own strand) and ``g_count`` (the count module's spatial
    imputation, on raw contained + cleaned boundaries) are **independent** signals — strand *direction* vs
    count *magnitude* — so the blend does not double-count. A node with no defined sense (κ=½ / AMBIG)
    takes the count imputation (``w=0``).

    The **FP-rate quantile knob** then reads ``g(q) = clip(center + Φ⁻¹(q)·σ)`` with the combined per-node
    std ``√(w²·σ²_strand + (1−w)²·σ²_count)``. ``q=0.5 ⇒ Φ⁻¹=0 ⇒ g=center`` — the no-op default.
    """
    mass_unspl = np.asarray(mass_unspl, dtype=np.float64)
    mass_spliced = np.asarray(mass_spliced, dtype=np.float64)
    sense = np.asarray(sense, dtype=np.float64)
    antisense = np.asarray(antisense, dtype=np.float64)
    strand_observable = np.asarray(strand_observable, dtype=bool)
    discrim = (2.0 * float(rna_sense_frac) - 1.0) ** 2  # per-fragment strand discriminability, in [0, 1]
    n_node = sense + antisense  # unspliced count evidence per node
    info = n_node * discrim  # carried strand information I = N·(2κ−1)²
    w = info / (info + float(info_scale))  # per-node strand-trust gradient: 0 at I=0, → 1 at high I
    z = float(ndtri(min(max(float(deconv_quantile), 1e-6), 1.0 - 1e-6)))  # Φ⁻¹(q); 0 at q=0.5

    g_count = np.clip(np.asarray(count_gdna_frac, dtype=np.float64), 0.0, 1.0)
    var_count = np.maximum(np.asarray(count_gdna_frac_var, dtype=np.float64), 0.0)
    center = g_count.copy()  # count-routed default (no usable strand signal ⇒ count only)
    var = var_count.copy()

    # STRAND module (Beta-Binomial posterior, weak prior) for every strand-routed node at once: a node
    # uses strand iff the library is strand-identifiable (discrim>0), the node is strand-observable, it
    # carries unspliced mass, and it has a sense split. The grid posterior is one (n_use, n_grid) batch.
    active = mass_unspl > 0.0
    use_strand = active & strand_observable & (n_node > 0.0) & (discrim > 0.0)
    if use_strand.any():
        idx = np.flatnonzero(use_strand)
        g_strand, var_strand = strand_posterior_gdna_frac(
            sense[idx],
            antisense[idx],
            rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            n_grid=n_grid,
        )
        wi = w[idx]  # per-node strand-trust gradient
        center[idx] = wi * g_strand + (1.0 - wi) * g_count[idx]
        var[idx] = wi**2 * var_strand + (1.0 - wi) ** 2 * var_count[idx]

    # FP-rate quantile g(q)=clip(center+Φ⁻¹(q)·σ); q=½ ⇒ z=0 ⇒ frac=center (bit-identical no-op).
    frac = center if z == 0.0 else np.clip(center + z * np.sqrt(var), 0.0, 1.0)
    # Inactive nodes (no unspliced mass): no gDNA, only the deterministic spliced RNA; zero variance.
    frac = np.where(active, frac, 0.0)
    gdna_frac_var = np.where(active, var, 0.0)
    gdna = frac * mass_unspl
    rna = (1.0 - frac) * mass_unspl + mass_spliced
    return NodeDeconv(gdna_mass=gdna, rna_mass=rna, gdna_frac=frac, gdna_frac_var=gdna_frac_var)


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
    info_scale=1.0,
) -> NodeDeconv:
    """Deconvolve each region's contained mass (a node) into gDNA / RNA by the strand/count blend.

    A region is strand-observable iff its transcript strand is defined (``TS_POS`` / ``TS_NEG``);
    ``TS_NONE`` (intergenic) and ``TS_AMBIG`` are count-only. The count fraction + its variance are the
    precomputed ``node_density.count_gdna_frac`` / ``count_gdna_frac_var``; the per-node weight
    ``w=I/(I+info_scale)``, ``I=N·(2κ−1)²``, and the FP-rate quantile are applied in ``_deconv_per_node``.
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
        info_scale=info_scale,
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


def _side_strand_orientation(view, same, ts_self, ts_other):
    """Per-side sense/antisense split + strand-observability (the ``TS_NONE``-wildcard rule).

    A side is **strand-observable** iff its two flanks define a single consistent transcript sense:
    ``{POS,POS}``/``{NEG,NEG}``, or a gene-edge ``{POS,NONE}``/``{NEG,NONE}`` (intergenic ``TS_NONE`` is a
    wildcard, oriented by the gene side); an opposite-strand ``{POS,NEG}`` or a ``TS_AMBIG`` flank leaves
    the sense undefined. Returns ``(sense, antisense, n_side, strand_observable)`` (all length R). Shared
    by the count-side quantities (:func:`_compute_side`) and the standalone strand module
    (:func:`strand_deconvolve`).
    """
    pos = view.n_unspliced_pos.astype(np.float64)
    neg = view.n_unspliced_neg.astype(np.float64)
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
    return sense, n_side - sense, n_side, strand_observable


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
    sense, antisense, n_side, strand_observable = _side_strand_orientation(
        view, same, ts_self, ts_other
    )
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
    info_scale=1.0,
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
            info_scale=info_scale,
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


def strand_deconvolve(
    substrate,
    region_arrays,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    deconv_quantile=0.5,
    n_grid=200,
) -> tuple[StrandSplit, StrandSplit, StrandSplit]:
    """Standalone strand deconvolution of every node into gDNA / unspliced-RNA, with likelihood info.

    The strand module of the sequential redesign (``docs/calibration/sequential_calibration_redesign.md``):
    for each node — the region-contained view and each boundary side — the Beta-Binomial strand posterior
    (:func:`strand_posterior_gdna_frac`) is read at the FP-quantile ``deconv_quantile`` to split the
    UNSPLICED mass into gDNA vs RNA, and the **likelihood Fisher information** ``I = N·(2κ−1)²`` (N = the
    unspliced count, κ = ``rna_sense_frac``) is emitted as the precision the count module weights by.
    Strand-uninformative nodes (AMBIG / no defined sense, or κ=½) get ``I=0`` and ``gdna_frac=NaN`` — the
    strand makes no guess; the count module imputes them from the density field. The deterministic spliced
    mass is untouched (RNA the count module folds in via the 3-term). Returns ``(contained, left, right)``
    :class:`StrandSplit`. Has **no count input** — Jeffreys-regularised and independent (the redesign's
    "strand runs first, standalone"). Built alongside the live blend; not yet wired into ``calibrate``.
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    r = ts.shape[0]
    w_strand = (2.0 * float(rna_sense_frac) - 1.0) ** 2  # per-fragment strand info (discriminability)

    def _split(mass_unspliced, sense, antisense, n_side, observable) -> StrandSplit:
        mass_unspliced = np.asarray(mass_unspliced, dtype=np.float64)
        gdna_frac = np.full(r, np.nan)
        info = np.zeros(r)
        informative = np.asarray(observable, dtype=bool) & (n_side > 0.0) & (w_strand > 0.0)
        if informative.any():
            idx = np.flatnonzero(informative)
            g_q, _ = strand_posterior_gdna_frac(
                sense[idx],
                antisense[idx],
                rna_sense_frac,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
                n_grid=n_grid,
                deconv_quantile=deconv_quantile,
            )
            gdna_frac[idx] = g_q
            info[idx] = n_side[idx] * w_strand  # I = N·(2κ−1)²
        gdna_mass = gdna_frac * mass_unspliced
        return StrandSplit(
            gdna_frac=gdna_frac, gdna_mass=gdna_mass, rna_mass=mass_unspliced - gdna_mass, info=info
        )

    # Contained view: a region is strand-observable iff its own transcript strand is defined (POS/NEG).
    c = substrate.contained
    c_pos = c.n_unspliced_pos.astype(np.float64)
    c_neg = c.n_unspliced_neg.astype(np.float64)
    c_sense = np.where(ts == TS_NEG, c_neg, c_pos)
    c_n = c_pos + c_neg
    contained = _split(c.mass_unspliced, c_sense, c_n - c_sense, c_n, (ts == TS_POS) | (ts == TS_NEG))

    # Boundary sides: orient against the neighbour strand (the TS_NONE-wildcard rule, shared helper).
    left_same, right_same = same_ref_left_right(ref_id)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        ts_prev[1:] = ts[:-1]
        ts_next[:-1] = ts[1:]
    l_sense, l_anti, l_n, l_obs = _side_strand_orientation(substrate.left, left_same, ts, ts_prev)
    r_sense, r_anti, r_n, r_obs = _side_strand_orientation(substrate.right, right_same, ts, ts_next)
    left = _split(substrate.left.mass_unspliced, l_sense, l_anti, l_n, l_obs)
    right = _split(substrate.right.mass_unspliced, r_sense, r_anti, r_n, r_obs)
    return contained, left, right


def cleaned_gdna_count(split: StrandSplit, raw_count: np.ndarray, info_scale: float) -> np.ndarray:
    """Strand-cleaned gDNA count for the count module, degrading gracefully to the raw count.

    The count module consumes this in place of the raw unspliced count (``pos+neg``) when building a
    region's gDNA density. The cleaning fraction is

        ``w·g_strand + (1−w)·1``       with   ``w = info / (info + info_scale)``

    where ``info = N·(2κ−1)²`` is the strand deconvolution's information (Phase 1 ``StrandSplit.info``)
    and ``g_strand`` is its gDNA fraction (``StrandSplit.gdna_frac``). The fraction slides **continuously**
    from **1** (all-gDNA = no RNA removed = the raw count, a *no-op*) when ``info → 0`` (κ≈½ or thin) to
    ``g_strand`` (full strand clean) when ``info → ∞``.

    Robust by construction: ``w`` and the meaningfulness of ``g_strand`` share the same ``info``, so an
    uninformative split — bounded to ``[0,1]`` by the Jeffreys grid, at worst the prior median, *never*
    the wild unbounded MLE — is discarded (``w≈0``) rather than trusted. So the count module may trust
    the result at every strand specificity, κ=½ ± fit-noise included. See
    ``docs/calibration/redesign_phase2_plan.md`` and ``strand_first_plan.md``.
    """
    info = np.asarray(split.info, dtype=np.float64)
    raw = np.asarray(raw_count, dtype=np.float64)
    if info.shape != raw.shape:
        raise ValueError(f"info {info.shape} and raw_count {raw.shape} must match")
    if info_scale <= 0.0:
        raise ValueError(f"info_scale must be > 0, got {info_scale}")
    w = info / (info + float(info_scale))  # 0 at info=0, → 1 as info → ∞
    # gdna_frac is NaN where info==0 (no strand sense); pin it to 0 there — w==0 discards it anyway.
    g = np.where(info > 0.0, split.gdna_frac, 0.0)
    clean_frac = w * g + (1.0 - w)
    return clean_frac * raw


__all__ = [
    "NodeDeconv",
    "StrandSplit",
    "strand_posterior_gdna_frac",
    "strand_deconvolve",
    "cleaned_gdna_count",
    "deconv_regions",
    "deconv_sides",
    "boundary_side_seeds",
]
