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

The reported fraction is the posterior **median**, after an optional gDNA-favoring FP-aversion
prior ``Beta(1+λ, 1)`` (``log_post += λ·log gdna_frac``) whose strength ``λ`` comes from the
``gdna_strand_confidence_z`` z-score (0 = neutral/flat; λ→∞ ⇒ siphon all unspliced mass into gDNA,
overwhelming any finite strand evidence — unlike a quantile, which only reports a point on the
*unchanged* posterior). Amounts use the
conserved **fractional mass** ``M``; precisions use **discrete counts**. gDNA = ``gdna_frac·M``; RNA =
``(1−gdna_frac)·M`` + the deterministic spliced mass. Regions and boundaries are combined into loci only
**after** calibration (``assemble_priors``). Explicit grid posterior (numerically stable).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .signature import TS_NEG, TS_POS
from .strand_likelihood import strand_loglik

_GRID_EPS = 1.0e-4
_JEFFREYS = 0.5  # principled prior floor (Beta(½,½) at zero count evidence), not a tunable


def strand_confidence_z_to_lambda(z: float) -> float:
    """Map the strand-deconvolution gDNA confidence z-score to the prior pseudocount ``λ``.

    The FP-aversion prior on a node's gDNA fraction is ``Beta(1+λ, 1)`` with a-priori mean
    ``Φ(z)`` (the standard-normal CDF), so ``λ = (2·Φ(z) − 1)/(1 − Φ(z))``: ``z=0 ⇒ Φ=0.5 ⇒ λ=0``
    (flat, neutral); ``z>0 ⇒ λ>0`` (favor gDNA); ``z→∞ ⇒ λ→∞`` (all gDNA); ``z<0 ⇒ λ∈(−1,0)``
    (lean RNA). It is a strictly increasing function of ``z`` (more confidence ⇒ more gDNA).
    """
    if z == 0.0:
        return 0.0
    phi = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    phi = min(max(phi, 1e-12), 1.0 - 1e-12)  # keep λ finite at both extremes
    return (2.0 * phi - 1.0) / (1.0 - phi)


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
    strand_overdispersion,
    gdna_strand_lambda,
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
                strand_overdispersion=strand_overdispersion,
            )
        # gDNA FP-aversion prior Beta(1+λ, 1): log_post += λ·log(gdna_frac) pulls the fraction
        # toward gDNA (λ from gdna_strand_confidence_z; 0 = neutral/flat). It competes with the
        # strand likelihood on the same log scale, so λ→∞ siphons even a confident-RNA node into
        # gDNA, while λ=0 leaves the unbiased posterior.
        if gdna_strand_lambda != 0.0:
            log_post = log_post + gdna_strand_lambda * log_grid
        w = np.exp(log_post - log_post.max())
        w /= w.sum()
        mean = float(np.dot(grid, w))
        var = float(np.dot((grid - mean) ** 2, w))
        # Report the posterior MEDIAN of the (FP-aversion-shifted) posterior. interp on the
        # CDF — exact, no Gaussian approximation. λ=0 ⇒ the unbiased median (neutral default).
        frac = float(np.clip(np.interp(0.5, np.cumsum(w), grid), 0.0, 1.0))
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
    strand_overdispersion=0.0,
    gdna_strand_confidence_z=0.0,
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
        strand_overdispersion=strand_overdispersion,
        gdna_strand_lambda=strand_confidence_z_to_lambda(gdna_strand_confidence_z),
        n_grid=n_grid,
    )


def deconv_sides(
    substrate,
    region_arrays,
    node_density,
    boundary_side_eff_len,
    *,
    rna_sense_frac,
    strand_overdispersion=0.0,
    gdna_strand_confidence_z=0.0,
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
    r = ts.shape[0]
    gdna_strand_lambda = strand_confidence_z_to_lambda(gdna_strand_confidence_z)
    boundary_count_observable = (
        node_density.boundary_count_observable
    )  # boundary (b, b+1) observable, indexed by left region
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density

    def _side(view, same, ts_other, side_boundary_count_observable):
        mass = view.mass_unspliced
        pos = view.n_unspliced_pos.astype(np.float64)
        neg = view.n_unspliced_neg.astype(np.float64)
        both_pos = same & (ts == TS_POS) & (ts_other == TS_POS)
        both_neg = same & (ts == TS_NEG) & (ts_other == TS_NEG)
        strand_observable = both_pos | both_neg
        sense = np.where(both_neg, neg, pos)
        antisense = (pos + neg) - sense
        with np.errstate(divide="ignore", invalid="ignore"):
            # Strand-clean the observable side's own crossing density (the boundary analogue of
            # the global-density cleaning): an exon/intron seam is count-observable (no SHARED
            # exon), but its crossing mass is gDNA + nascent read-through. Multiply by the
            # closed-form gDNA fraction (sense_frac - rna_sense_frac)/(1/2 - rna_sense_frac) where
            # the side is strand-observable, so the count clue sees clean gDNA, not nascent.
            # Unstranded or non-oriented sides keep raw mass.
            n_side = pos + neg
            sense_frac = np.where(n_side > 0.0, sense / np.maximum(n_side, 1e-9), 0.5)
            if abs(0.5 - rna_sense_frac) > 1e-6:
                clean_gdna_frac = np.clip(
                    (sense_frac - rna_sense_frac) / (0.5 - rna_sense_frac), 0.0, 1.0
                )
            else:
                clean_gdna_frac = np.ones_like(mass)
            clean_gdna_frac = np.where(strand_observable, clean_gdna_frac, 1.0)
            own = np.where(
                (mass > 0.0) & (eff > 0.0), clean_gdna_frac * mass / np.maximum(eff, 1e-12), 0.0
            )
        density = np.where(
            side_boundary_count_observable, own, region_density
        )  # observable → clean own; else swept
        return _joint_per_node(
            mass,
            view.mass_spliced,
            sense,
            antisense,
            density,
            pos + neg,
            eff,
            strand_observable,
            rna_sense_frac=rna_sense_frac,
            strand_overdispersion=strand_overdispersion,
            gdna_strand_lambda=gdna_strand_lambda,
            n_grid=n_grid,
        )

    # LEFT view of r = right side of boundary (r-1, r); neighbour is r-1.
    left_same = np.zeros(r, dtype=bool)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_boundary_count_observable = np.zeros(r, dtype=bool)
    if r > 1:
        left_same[1:] = ref_id[1:] == ref_id[:-1]
        ts_prev[1:] = ts[:-1]
        left_boundary_count_observable[1:] = boundary_count_observable[:-1]
    left = _side(substrate.left, left_same, ts_prev, left_boundary_count_observable)

    # RIGHT view of r = left side of boundary (r, r+1); neighbour is r+1.
    right_same = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        right_same[:-1] = ref_id[:-1] == ref_id[1:]
        ts_next[:-1] = ts[1:]
    right = _side(substrate.right, right_same, ts_next, boundary_count_observable)
    return left, right


__all__ = ["JointDeconv", "deconv_regions", "deconv_sides"]
