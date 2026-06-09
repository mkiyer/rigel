"""Phase 3 — the joint deconvolution: per-node gDNA/RNA deconvolution from count × strand.

Each **node** (every region and every boundary, kept separate) is deconvolved into gDNA / RNA
mass by combining the two **conditionally-independent** clues as a per-node posterior over the
gDNA fraction ``gdna_frac``::

    posterior(gdna_frac) ∝ Beta(gdna_frac ; gdna_frac^count, count_concentration)  ·  BB_strand(sense, antisense | gdna_frac)

- **count prior** ``Beta(a_c, b_c)``: mean ``count_gdna_frac = clip(density*eff_len / M, 0, 1)`` from
  the **strand-cleaned** global density (``density_model.gdna_frac``), concentration
  ``count_concentration`` = the expected gDNA count ``count_evidence = density*eff_len``.
  Jeffreys-floored, so a node expecting ~no gDNA gives Beta(1/2, 1/2) and defers to the strand clue.
- **strand likelihood** (Phase 2, Beta-Binomial): :func:`strand_likelihood.strand_loglik`. Absent
  (flat) for AMBIG / zero-flux nodes ⇒ count-only; degenerate for unstranded (rna_sense_frac = 0.5).

Because gdna_density_global is strand-cleaned upstream, the count clue is correct on every node — nascent introns
read as RNA (density ≫ clean gdna_density_global), silent gDNA genes as gDNA — and the two clues reinforce where
they agree and arbitrate where they disagree (the strand-specificity spectrum: unstranded ⇒
count-only, increasingly stranded ⇒ increasingly strand-informed). The clues are conditionally
independent given (gdna_density_global, rna_sense_frac); a node's ~1/N self-contribution to the global gdna_density_global is negligible.

The reported fraction is the posterior **median**, then shifted in LOG-ODDS by the FP-aversion
bias ``gdna_strand_llr_bias``: ``gdna_frac ← σ(λ + logit(median))`` — the calibration-stage twin of
the EM ``gdna_em_llr_bias`` and its exact analog: a flat per-node gDNA log-odds shift in nats, ``0`` =
neutral, ``>0`` ⇒ gDNA, ``<0`` ⇒ RNA, ``λ→+∞`` siphons all unspliced mass into gDNA. It is applied
to the point estimate (count-independent, so the nats are portable) rather than as a tilt on the
posterior (which a high-count node would drown out) or a quantile (which only moves a point on the
*unchanged* posterior). Amounts use the
conserved **fractional mass** ``M``; precisions use **discrete counts**. gDNA = ``gdna_frac·M``; RNA =
``(1−gdna_frac)·M`` + the deterministic spliced mass. Regions and boundaries are combined into loci only
**after** calibration (``assemble_priors``). Explicit grid posterior (numerically stable).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_model import strand_clean_gdna_frac
from .signature import TS_NEG, TS_POS
from .strand_likelihood import strand_loglik

_GRID_EPS = 1.0e-4
_JEFFREYS = 0.5  # principled prior floor (Beta(½,½) at zero count evidence), not a tunable


@dataclass(frozen=True, slots=True)
class JointDeconv:
    """Per-node deconvolution result (regions or boundaries; kept separate)."""

    gdna_mass: np.ndarray  # float64[K]
    rna_mass: np.ndarray  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac: np.ndarray  # float64[K] — reported gDNA fraction of the UNSPLICED mass
    gdna_frac_var: np.ndarray  # float64[K] — posterior variance


def _joint_per_node(
    mass_unspl,
    mass_spliced,
    sense,
    antisense,
    density,
    count_evidence,
    eff_len,
    strand_observable,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion,
    rna_strand_overdispersion,
    gdna_strand_llr_bias,
    n_grid,
) -> JointDeconv:
    k = mass_unspl.shape[0]
    grid = np.linspace(_GRID_EPS, 1.0 - _GRID_EPS, n_grid)
    log_grid = np.log(grid)
    log_1mgrid = np.log1p(-grid)
    gdna = np.zeros(k)
    rna = np.zeros(k)
    gdna_frac = np.zeros(k)
    gdna_frac_var = np.zeros(k)
    for i in range(k):
        m = float(mass_unspl[i])
        if m <= 0.0:
            rna[i] = float(mass_spliced[i])  # only deterministic spliced RNA, if any
            continue
        # Joint deconvolution: count prior Beta(a_c, b_c) (mean from the strand-cleaned global
        # density, precision = expected gDNA count count_evidence) x strand likelihood. The global
        # density is already strand-cleaned (density_model.gdna_frac), so the count clue is correct
        # on every node — nascent introns read as RNA (density >> clean global density), silent gDNA
        # genes as gDNA — and *reinforces* the
        # strand where they agree. Strand absent (flat) for AMBIG / zero-flux nodes ⇒ count-only;
        # unstranded (rna_sense_frac = 0.5) ⇒ count-only everywhere. count and strand are conditionally
        # independent given (gdna_density_global, rna_sense_frac); the node's ~1/N self-contribution to gdna_density_global is negligible.
        count_gdna_frac = min(max(density[i] * eff_len[i] / m, 0.0), 1.0)
        count_concentration = max(float(count_evidence[i]), 0.0)
        a_c = count_concentration * count_gdna_frac + _JEFFREYS
        b_c = count_concentration * (1.0 - count_gdna_frac) + _JEFFREYS
        log_post = (a_c - 1.0) * log_grid + (b_c - 1.0) * log_1mgrid
        if strand_observable[i] and (sense[i] + antisense[i]) > 0:
            log_post = log_post + strand_loglik(
                grid,
                sense[i],
                antisense[i],
                rna_sense_frac,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
            )
        w = np.exp(log_post - log_post.max())
        w /= w.sum()
        mean = float(np.dot(grid, w))
        var = float(np.dot((grid - mean) ** 2, w))
        # Unbiased posterior median (interp on the CDF — exact, no Gaussian approximation).
        frac = float(np.clip(np.interp(0.5, np.cumsum(w), grid), 0.0, 1.0))
        # gDNA FP-aversion LLR bias: shift the reported gDNA fraction's LOG-ODDS by λ nats —
        # frac ← σ(λ + logit(frac)) — the calibration-stage twin of the EM gdna_em_llr_bias, and its
        # exact analog: a flat per-node log-odds shift toward gDNA, independent of the node's
        # fragment count (a *tilt on the posterior* would instead be drowned out by high-count
        # nodes, so its nats would not be portable). λ=0 ⇒ unbiased median (neutral default);
        # λ→+∞ ⇒ frac→1 (all gDNA, even a confident-RNA node); λ→−∞ ⇒ frac→0 (all RNA), symmetric.
        if gdna_strand_llr_bias != 0.0 and 0.0 < frac < 1.0:
            logit_frac = np.log(frac) - np.log1p(-frac)
            frac = float(1.0 / (1.0 + np.exp(-(gdna_strand_llr_bias + logit_frac))))
        gdna_frac[i] = frac
        gdna_frac_var[i] = var
        gdna[i] = frac * m
        rna[i] = (1.0 - frac) * m + float(mass_spliced[i])
    return JointDeconv(
        gdna_mass=gdna, rna_mass=rna, gdna_frac=gdna_frac, gdna_frac_var=gdna_frac_var
    )


def deconv_regions(
    substrate,
    region_arrays,
    node_density,
    region_eff_len,
    *,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    gdna_strand_llr_bias=0.0,
    n_grid=200,
) -> JointDeconv:
    """Deconvolve each region's contained mass (a node) into gDNA / RNA."""
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    sense = np.where(ts == TS_NEG, neg, pos)  # orient to transcript sense
    antisense = (pos + neg) - sense
    strand_observable = (ts == TS_POS) | (ts == TS_NEG)
    return _joint_per_node(
        c.mass_unspliced,
        c.mass_spliced,
        sense,
        antisense,
        node_density.density,
        node_density.count_evidence,
        np.asarray(region_eff_len, dtype=np.float64),
        strand_observable,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        gdna_strand_llr_bias=gdna_strand_llr_bias,
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
    density: np.ndarray  # count-prior density: clean own crossing density, else swept region
    gdna_weight: np.ndarray  # strand-cleaned own gDNA fraction (mean-based; overdispersion-free)
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
    left_same = np.zeros(r, dtype=bool)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_observable = np.zeros(r, dtype=bool)
    right_same = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        left_same[1:] = ref_id[1:] == ref_id[:-1]
        ts_prev[1:] = ts[:-1]
        left_observable[1:] = boundary_count_observable[:-1]
        right_same[:-1] = ref_id[:-1] == ref_id[1:]
        ts_next[:-1] = ts[1:]
    return left_same, ts_prev, left_observable, right_same, ts_next, boundary_count_observable


def _compute_side(
    view, same, ts_self, ts_other, side_count_observable, eff, region_density, rna_sense_frac
) -> _SideQuantities:
    """Per-side sense split, strand-cleaned gDNA weight, and count-prior density.

    A side is **strand-observable** iff its boundary's two regions share a single transcript
    strand (so 'sense' is defined). Its own crossing mass is gDNA + nascent read-through; the
    closed-form strand-cleaned gDNA fraction ``(sense_frac − rna_sense_frac)/(½ − rna_sense_frac)``
    (mean-based, **overdispersion-free**) serves as both the count-prior density cleaner and the
    strand-fit seed weight. A side is **count-observable** iff its boundary has no shared exon
    (exon–intron / exon–intergenic seam) — those borrow their own clean crossing density; others
    borrow the swept region density.
    """
    mass = view.mass_unspliced
    pos = view.n_unspliced_pos.astype(np.float64)
    neg = view.n_unspliced_neg.astype(np.float64)
    both_pos = same & (ts_self == TS_POS) & (ts_other == TS_POS)
    both_neg = same & (ts_self == TS_NEG) & (ts_other == TS_NEG)
    strand_observable = both_pos | both_neg
    sense = np.where(both_neg, neg, pos)
    n_side = pos + neg
    antisense = n_side - sense
    with np.errstate(divide="ignore", invalid="ignore"):
        # Same closed-form strand clean as the count clue (one shared implementation); a side with
        # no defined sense (not strand-observable) keeps its full crossing mass (frac 1).
        clean_gdna_frac = strand_clean_gdna_frac(sense, n_side, rna_sense_frac)
        clean_gdna_frac = np.where(strand_observable, clean_gdna_frac, 1.0)
        own = np.where(
            (mass > 0.0) & (eff > 0.0), clean_gdna_frac * mass / np.maximum(eff, 1e-12), 0.0
        )
    density = np.where(side_count_observable, own, region_density)
    return _SideQuantities(
        sense=sense,
        antisense=antisense,
        n_side=n_side,
        mass=mass,
        mass_spliced=view.mass_spliced,
        density=density,
        gdna_weight=clean_gdna_frac,
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
    gdna_strand_llr_bias=0.0,
    n_grid=200,
) -> tuple[JointDeconv, JointDeconv]:
    """Deconvolve each boundary **side** as an independent node (R1/decision iii).

    The deconvolution unit is the boundary-side. Region ``r`` owns the **right** side of its
    left boundary (``substrate.left[r]``) and the **left** side of its right boundary
    (``substrate.right[r]``); both lie inside region ``r`` and so use
    ``boundary_side_eff_len[r] = E_FL[min(ℓ, L_r)]`` (R2/R3). A side is **count-observable** iff
    its boundary is (no shared exon) ⇒ ``count_gdna_frac → 1`` from its own crossing density; otherwise
    it borrows the swept region density. It is **strand-observable** iff the boundary's two
    regions share a single transcript strand. Returns ``(left, right)`` per-region JointDecodes
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
        sq = _compute_side(view, same, ts, ts_other, side_obs, eff, region_density, rna_sense_frac)
        return _joint_per_node(
            sq.mass,
            sq.mass_spliced,
            sq.sense,
            sq.antisense,
            sq.density,
            sq.n_side,
            eff,
            sq.strand_observable,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            gdna_strand_llr_bias=gdna_strand_llr_bias,
            n_grid=n_grid,
        )

    left = _deconv(substrate.left, l_same, ts_prev, l_obs)
    right = _deconv(substrate.right, r_same, ts_next, r_obs)
    return left, right


def boundary_side_seeds(
    substrate, region_arrays, node_density, boundary_side_eff_len, *, rna_sense_frac
):
    """``(sense, total, gdna_weight)`` seed arrays from count- & strand-observable boundary sides.

    The exon–intron / exon–intergenic seam seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the contained-region seeds (needed under hybrid capture,
    which depletes off-target intergenic/intronic gDNA). The weight is the strand-cleaned own
    gDNA fraction — the same overdispersion-free quantity the contained-region seeds use.
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )
    sides = (
        _compute_side(
            substrate.left, l_same, ts, ts_prev, l_obs, eff, region_density, rna_sense_frac
        ),
        _compute_side(
            substrate.right, r_same, ts, ts_next, r_obs, eff, region_density, rna_sense_frac
        ),
    )
    senses, totals, weights = [], [], []
    for sq in sides:
        seed = sq.count_observable & sq.strand_observable & (sq.n_side > 0.0)
        senses.append(sq.sense[seed])
        totals.append(sq.n_side[seed])
        weights.append(sq.gdna_weight[seed])
    return np.concatenate(senses), np.concatenate(totals), np.concatenate(weights)


__all__ = ["JointDeconv", "deconv_regions", "deconv_sides", "boundary_side_seeds"]
