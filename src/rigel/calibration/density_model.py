"""Phase 1 — the density model ("count clue"): per-region gDNA density from OBSERVED counts.

Acyclic by construction: the gDNA density is read **directly** from count-observable nodes
(where fragments are gDNA by construction) and **imputed locally** for the rest. It never consults
a global ``gdna_density_global * L`` product, so there is no density->deconv->density feedback loop.

Count-observability is a property of the region **signature** (4-bit exon/intron ± flags):

* **region** is observable ⇔ it has **no exon bit** (intergenic or intron-only). Its
  contained unspliced mass is gDNA (+ nascent RNA — an upper bias the strand clue removes);
  an exonic region's contained mass is contaminated by mature RNA.
* **boundary** is observable ⇔ **no exon bit is shared** across its two sides → no single
  exon-strand continues across it → no *unspliced mature RNA* crosses (a single-exon
  transcript spanning the seam would put unspliced mature RNA there). Its crossing
  unspliced mass is then gDNA(+nascent).

Raw counts (no strand cleaning): in the decoupled architecture the strand module owns the strand
channel; this count module is the fallback for strand-unobservable / unstranded nodes and works on
the raw unspliced count. The local imputation (no global sweep):

* **observable region** (intron/intergenic): its own contained ``count / region_eff_len``.
* **non-observable region** (exon/AMBIG): the gDNA density of each *observable boundary side*
  (``crossing count / fl_mean`` — the accumulator deposits the molecule's true span, so the
  one-side crossing flux is an unbiased density estimator), averaged over the available sides.
* **run interiors** (consecutive non-observable regions, no observable side): the nearest anchored
  neighbour, carried inward from the run's observable edges (forward + reverse, averaged).
* **no anchor in the whole reference**: the global gDNA density (the count-weighted mean of the
  count-observable regions' densities) — a rare fallback.

Counts → density via the gDNA FL effective lengths: contained ``count ÷ region_eff_len``, crossing
``count ÷ fl_mean``. For uniform gDNA at a given density both recover that density.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .run_fill import runfill_bidirectional, same_ref_left_right
from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class NodeDensity:
    """Per-region gDNA density (count clue) after local imputation."""

    density: np.ndarray  # float64[R] — local gDNA density (fragments per effective bp)
    count_gdna_frac: (
        np.ndarray
    )  # float64[R] — count module's gDNA fraction g_count = clip(density·region_eff_len /
    #   contained_mass): the gDNA fraction of the contained unspliced mass from the (raw) local gDNA
    #   density. Consumed by the count module (strand-unobservable / unstranded nodes) and the gDNA
    #   strand-fit seed weight.
    count_gdna_frac_var: np.ndarray  # float64[R] — Phase 1 count posterior variance σ_g² of
    #   count_gdna_frac: μ_g²·v_rel, v_rel = relative variance (CV²) from the NON-PARAMETRIC kNN
    #   variance~density curve (imputed), 1/N_own (observable), 1.0 (no anchor); floored by Poisson
    #   1/N_anchor; capped at μ_g(1−μ_g). Feeds the FP-rate quantile (Phase 2); NOT yet the blend.
    region_count_observable: np.ndarray  # bool[R] — count-observable region (non-exonic)
    boundary_count_observable: np.ndarray  # bool[R] — count-observable boundary right of region r
    n_region_count_observable: int
    n_boundary_count_observable: int


def count_observable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-observability for regions and (right-) boundaries.

    Returns ``(region_count_observable, boundary_count_observable)``, both ``bool[R]``. ``boundary_count_observable[r]``
    describes the internal boundary between region ``r`` and ``r+1`` (defined iff same ref).
    """
    sig = np.asarray(signature).astype(np.int64)
    ref = np.asarray(ref_id)
    r = sig.shape[0]
    region_count_observable = (sig & _EXON_BITS) == 0
    boundary_count_observable = np.zeros(r, dtype=bool)
    if r > 1:
        same = ref[:-1] == ref[1:]
        shared_exon = (sig[:-1] & sig[1:] & _EXON_BITS) != 0
        boundary_count_observable[:-1] = same & ~shared_exon
    return region_count_observable, boundary_count_observable


_LOESS_SPAN = 0.4  # LOESS bandwidth (fraction of points per local fit); CV-selectable later
_LOESS_MIN_FIT = 5  # need at least this many 2-anchor points to fit a local-linear curve
_BISQUARE_C = 6.0  # Cleveland's robustness constant (6·MAD); the standard LOWESS value, not tuned


def _loess(x: np.ndarray, y: np.ndarray, xq: np.ndarray, frac: float, robust_iters: int = 2):
    """Robust LOESS — local **linear** regression (tricube kernel) of y on x, evaluated at ``xq``.

    The standard Cleveland smoother: per query point, fit a weighted line over the ``frac·n`` nearest
    points (tricube distance weights, adaptive bandwidth = the k-th nearest distance), and read the line
    at the query. ``robust_iters`` passes reweight by the residual bisquare (``6·MAD``) to downweight
    outliers — important here because variance estimates are heavy-tailed. Local *linear* (not a local
    mean) corrects boundary bias; distance weighting (not k-by-rank) keeps a far cluster across a density
    gap from contaminating the fit. ~O(n_query·n) per pass — sub-ms at our region counts.
    """
    n = x.shape[0]
    k = int(np.clip(int(frac * n), 2, n))
    rob = np.ones(n)

    def fit_at(queries: np.ndarray) -> np.ndarray:
        out = np.empty(queries.shape[0])
        for i, xi in enumerate(queries):
            d = np.abs(x - xi)
            h = np.partition(d, k - 1)[k - 1]  # distance to the k-th nearest → local bandwidth
            w = np.clip(1.0 - (d / max(h, 1e-12)) ** 3, 0.0, 1.0) ** 3 * rob  # tricube × robustness
            sw = w.sum()
            if sw <= 0.0:
                out[i] = np.average(y, weights=rob) if rob.sum() > 0 else y.mean()
                continue
            mx = (w * x).sum() / sw
            my = (w * y).sum() / sw
            sxx = (w * (x - mx) ** 2).sum()
            b = (w * (x - mx) * (y - my)).sum() / sxx if sxx > 1e-12 else 0.0  # else local mean
            out[i] = my + b * (xi - mx)
        return out

    for _ in range(robust_iters):
        resid = y - fit_at(x)
        mad = np.median(np.abs(resid))
        if mad <= 0.0:
            break
        rob = np.clip(1.0 - (resid / (_BISQUARE_C * mad)) ** 2, 0.0, 1.0) ** 2
    return fit_at(xq)


def _count_fraction_variance(
    count_gdna_frac: np.ndarray,
    density: np.ndarray,
    *,
    own: np.ndarray,
    own_count: np.ndarray,
    d_left: np.ndarray,
    d_right: np.ndarray,
    n_anchor: np.ndarray,
    frac: float = _LOESS_SPAN,
) -> np.ndarray:
    """Phase-1 count posterior variance ``σ_g² = μ_g²·v_rel`` via a **non-parametric LOESS** variance~mean
    curve.

    Hybrid capture is bimodal (on/off-target chasm in both mean and variance), so no single parametric
    ``α·μ²`` law fits. Instead the imputation variance is a robust **LOESS curve in log-log space**:

    1. each 2-anchor node (both ``d_left``/``d_right`` finite, not observable) gives a raw density variance
       ``¼(d_L−d_R)²`` at mean density ``½(d_L+d_R)``;
    2. fit ``log σ²_density ~ log μ_ρ`` by robust local-linear LOESS (:func:`_loess`) — log-log because
       both span orders of magnitude (``var∝mean²`` is a line there) and it handles the chasm without a
       global law; local-linear corrects boundary bias and distance weighting avoids gap contamination;
    3. every region reads ``σ²_density`` off the curve at its own ``density`` → relative variance
       ``v_rel = σ²_density/density²`` (invariant under the linear density→fraction map).

    ``v_rel`` is then ``1/N_own`` for observable nodes (Poisson on own count), ``max(curve, 1/N_anchor)``
    for imputed nodes (LOESS curve floored by the anchors' Poisson count noise — Q1), and ``1.0`` for a
    node with no anchor anywhere (uninformative → defer). ``σ_g² = μ_g²·v_rel``, capped at ``μ_g(1−μ_g)``.

    ``frac`` is the LOESS span (the one smoothing knob; CV-selectable later — default ``_LOESS_SPAN``).
    """
    r = count_gdna_frac.shape[0]
    two_anchor = np.isfinite(d_left) & np.isfinite(d_right) & ~own
    mu_d = 0.5 * (d_left + d_right)
    raw_var = 0.25 * (d_left - d_right) ** 2

    fit_sel = two_anchor & (mu_d > _EPS) & (raw_var > _EPS)  # log-log fit drops exact-agreement zeros
    loess_v_rel = np.full(r, np.nan, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        if int(fit_sel.sum()) >= _LOESS_MIN_FIT:
            x_fit = np.log(mu_d[fit_sel])
            y_fit = np.log(raw_var[fit_sel])
            valid = density > _EPS
            sigma_d2 = np.exp(_loess(x_fit, y_fit, np.log(density[valid]), frac))
            loess_v_rel[valid] = sigma_d2 / density[valid] ** 2

        v_rel = np.full(r, 1.0, dtype=np.float64)  # default: no anchor → uninformative
        v_rel[own] = np.where(own_count[own] > 0.0, 1.0 / np.maximum(own_count[own], _EPS), 1.0)
        imputed = (~own) & (n_anchor > 0.0)
        floor = 1.0 / np.maximum(n_anchor, _EPS)  # Poisson count-noise floor (Q1)
        lv = np.where(np.isfinite(loess_v_rel), loess_v_rel, 0.0)
        v_rel[imputed] = np.maximum(lv[imputed], floor[imputed])
    return np.minimum(count_gdna_frac**2 * v_rel, count_gdna_frac * (1.0 - count_gdna_frac))


def node_gdna_density(
    substrate,
    region_arrays,
    region_eff_len: np.ndarray,
    fl_mean: float,
) -> NodeDensity:
    """Per-region gDNA density from the count clue via LOCAL boundary-anchored imputation.

    The count module estimates gDNA density from **raw** unspliced counts — no strand cleaning (in
    the decoupled architecture the strand module owns the strand channel; the count module is the
    fallback for strand-unobservable / unstranded nodes). The density is read directly from
    count-observable regions and imputed locally for the rest; a region with no local anchor anywhere
    takes the count-weighted-mean observable density as a global fallback. See the module docstring.
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    region_count_observable, boundary_count_observable = count_observable_masks(sig, ref_id)

    # Raw unspliced COUNT per node (no strand cleaning — the count module works on raw counts).
    def total_count(view) -> np.ndarray:
        pos = np.asarray(view.n_unspliced_pos, dtype=np.float64)
        neg = np.asarray(view.n_unspliced_neg, dtype=np.float64)
        return pos + neg

    contained_gdna = total_count(substrate.contained)
    left_gdna = total_count(substrate.left)  # right side of region r's LEFT boundary
    right_gdna = total_count(substrate.right)  # left side of region r's RIGHT boundary

    # Per-side boundary observability for region r: its LEFT side uses boundary (r−1, r); its RIGHT
    # side uses boundary (r, r+1). boundary_count_observable[k] describes boundary (k, k+1).
    left_same, right_same = same_ref_left_right(ref_id)
    left_anchor = np.zeros(r, dtype=bool)
    right_anchor = np.zeros(r, dtype=bool)
    if r > 1:
        left_anchor[1:] = boundary_count_observable[:-1] & left_same[1:]
        right_anchor[:-1] = boundary_count_observable[:-1] & right_same[:-1]

    density = np.full(r, np.nan, dtype=np.float64)
    # Observable region with a usable contained length → its own contained density.
    own = region_count_observable & (region_eff_len > _EPS)
    density[own] = contained_gdna[own] / region_eff_len[own]
    # Everything else: anchor from the available observable boundary sides (crossing count / fl_mean).
    inv_fl = 1.0 / fl_mean if fl_mean > 0.0 else 0.0
    # Per-node anchor bookkeeping for the Phase-1 count posterior variance (below): the two anchor
    # densities (d_L, d_R) and the total anchor crossing count behind the imputation.
    d_left = np.where(left_anchor, left_gdna * inv_fl, np.nan)
    d_right = np.where(right_anchor, right_gdna * inv_fl, np.nan)
    n_anchor = left_anchor * left_gdna + right_anchor * right_gdna  # crossing counts behind the impute
    for i in np.where(~own)[0]:
        est = []
        if left_anchor[i]:
            est.append(left_gdna[i] * inv_fl)
        if right_anchor[i]:
            est.append(right_gdna[i] * inv_fl)
        if est:
            density[i] = float(np.mean(est))

    # Run-fill: a region still unset (a run interior with no observable side) inherits the nearest
    # anchored neighbour, carried inward from both directions within its reference and averaged.
    density = runfill_bidirectional(density, ref_id)

    # A region still unset has NO local anchor anywhere in its reference. It takes the GLOBAL gDNA
    # density (the count-weighted mean of the count-observable regions' own densities — a sensible
    # baseline, not 0). A zero-gDNA library has ~no observable mass ⇒ baseline ≈ 0 ⇒ no manufactured
    # gDNA. (A rare fallback — local imputation handles the rest.)
    anchored = ~np.isnan(density)
    own_len = float(np.sum(region_eff_len[own]))
    global_density = float(np.sum(contained_gdna[own]) / own_len) if own_len > 0.0 else 0.0
    density = np.where(anchored, density, global_density)

    # Count-prior MEAN g_count = clip(density·eff_len / contained_mass): the gDNA fraction of the
    # contained unspliced mass implied by the local gDNA density. Consumed by the count module
    # (strand-unobservable / unstranded nodes) and by the gDNA strand-fit seed weight.
    contained_mass = np.asarray(substrate.contained.mass_unspliced, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(
            np.where(contained_mass > 0.0, density * region_eff_len / contained_mass, 0.0), 0.0, 1.0
        )

    # Phase 1: count posterior variance σ_g² of the gDNA fraction (see _count_fraction_variance).
    count_gdna_frac_var = _count_fraction_variance(
        count_gdna_frac, density, own=own, own_count=contained_gdna, d_left=d_left, d_right=d_right,
        n_anchor=n_anchor,
    )

    return NodeDensity(
        density=density,
        count_gdna_frac=count_gdna_frac,
        count_gdna_frac_var=count_gdna_frac_var,
        region_count_observable=region_count_observable,
        boundary_count_observable=boundary_count_observable,
        n_region_count_observable=int(region_count_observable.sum()),
        n_boundary_count_observable=int(boundary_count_observable.sum()),
    )


__all__ = [
    "NodeDensity",
    "count_observable_masks",
    "node_gdna_density",
]
