"""Phase 3 — the joint decode: per-node gDNA/RNA deconvolution from count + strand.

Each **node** (every region and every boundary, kept separate) is deconvolved independently
into gDNA mass and RNA mass by combining the two **orthogonal** clues as a per-node posterior
over the gDNA fraction ``π_g``::

    posterior(π_g) ∝ Beta(π_g ; π_g^count, κ_c)  ·  BB_strand(sense, antisense | π_g)

- **count prior** ``Beta(a_c, b_c)``: mean ``π_g^count = clip(density·eff_len / M, 0, 1)``,
  concentration ``κ_c`` = the swept DISCRETE count evidence (Phase 1). Jeffreys-floored
  (``a_c, b_c ≥ ½``), so zero count evidence ⇒ Beta(½,½) and the strand decides.
- **strand likelihood** (Phase 2, Beta-Binomial): :func:`strand_decode.strand_loglik`.
  Absent (flat) for AMBIG / zero-flux nodes ⇒ count-only.

The reported fraction is the posterior quantile selected by ``confidence`` (z; 0 = mean).
Amounts use the conserved **fractional mass** ``M``; precisions use **discrete counts** —
the two are never conflated. gDNA = ``π_g·M``; RNA = ``(1−π_g)·M`` + the deterministic
spliced mass. Regions and boundaries are combined into loci only **after** calibration.

The implementation is the explicit grid posterior (numerically stable, §7 of the plan); the
O(R) Gaussian-product fast path is a Phase-6 optimization.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import TS_NEG, TS_POS
from .strand_decode import strand_loglik

_PI_EPS = 1.0e-4
_JEFFREYS = 0.5  # principled prior floor (Beta(½,½) at zero count evidence), not a tunable


@dataclass(frozen=True, slots=True)
class JointDecode:
    """Per-node deconvolution result (regions or boundaries; kept separate)."""

    gdna_mass: np.ndarray  # float64[K]
    rna_mass: np.ndarray  # float64[K]  (= (1−π_g)·M_unspliced + spliced mass)
    pi_g: np.ndarray  # float64[K] — reported gDNA fraction of the UNSPLICED mass
    pi_g_var: np.ndarray  # float64[K] — posterior variance


def _joint_per_node(
    mass_unspl, mass_spliced, sense, antisense, density, count_evidence, eff_len, strand_decodable,
    *, kappa_rna, strand_overdispersion, confidence, n_grid,
) -> JointDecode:
    k = mass_unspl.shape[0]
    grid = np.linspace(_PI_EPS, 1.0 - _PI_EPS, n_grid)
    log_grid = np.log(grid)
    log_1mgrid = np.log1p(-grid)
    gdna = np.zeros(k)
    rna = np.zeros(k)
    pi_g = np.zeros(k)
    pi_var = np.zeros(k)
    for i in range(k):
        m = float(mass_unspl[i])
        if m <= 0.0:
            rna[i] = float(mass_spliced[i])  # only deterministic spliced RNA, if any
            continue
        # count prior Beta(a_c, b_c): mean from fractional mass, concentration from discrete count
        pi_count = min(max(density[i] * eff_len[i] / m, 0.0), 1.0)
        kc = max(float(count_evidence[i]), 0.0)
        a_c = kc * pi_count + _JEFFREYS
        b_c = kc * (1.0 - pi_count) + _JEFFREYS
        log_post = (a_c - 1.0) * log_grid + (b_c - 1.0) * log_1mgrid
        if strand_decodable[i] and (sense[i] + antisense[i]) > 0:
            log_post = log_post + strand_loglik(
                grid, sense[i], antisense[i], kappa_rna, strand_overdispersion=strand_overdispersion
            )
        w = np.exp(log_post - log_post.max())
        w /= w.sum()
        mean = float(np.dot(grid, w))
        var = float(np.dot((grid - mean) ** 2, w))
        frac = min(max(mean + confidence * np.sqrt(var), 0.0), 1.0)  # z-quantile (Gaussian)
        pi_g[i] = frac
        pi_var[i] = var
        gdna[i] = frac * m
        rna[i] = (1.0 - frac) * m + float(mass_spliced[i])
    return JointDecode(gdna_mass=gdna, rna_mass=rna, pi_g=pi_g, pi_g_var=pi_var)


def decode_regions(
    substrate, region_arrays, node_density, region_eff_len, *,
    kappa_rna, strand_overdispersion=0.0, confidence=0.0, n_grid=200,
) -> JointDecode:
    """Deconvolve each region's contained mass (a node) into gDNA / RNA."""
    ts = np.asarray(region_arrays.ts_class)
    c = substrate.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    sense = np.where(ts == TS_NEG, neg, pos)  # orient to transcript sense
    antisense = (pos + neg) - sense
    strand_dec = (ts == TS_POS) | (ts == TS_NEG)
    return _joint_per_node(
        c.mass_unspliced, c.mass_spliced, sense, antisense, node_density.density,
        node_density.count_evidence, np.asarray(region_eff_len, dtype=np.float64), strand_dec,
        kappa_rna=kappa_rna, strand_overdispersion=strand_overdispersion,
        confidence=confidence, n_grid=n_grid,
    )


def decode_sides(
    substrate, region_arrays, node_density, boundary_side_eff_len, *,
    kappa_rna, strand_overdispersion=0.0, confidence=0.0, n_grid=200,
) -> tuple[JointDecode, JointDecode]:
    """Deconvolve each boundary **side** as an independent node (R1/decision iii).

    The deconvolution unit is the boundary-side. Region ``r`` owns the **right** side of its
    left boundary (``substrate.left[r]``) and the **left** side of its right boundary
    (``substrate.right[r]``); both lie inside region ``r`` and so use
    ``boundary_side_eff_len[r] = E_FL[min(ℓ, L_r)]`` (R2/R3). A side is **count-decodable** iff
    its boundary is (no shared exon) ⇒ ``π_count → 1`` from its own crossing density; otherwise
    it borrows the swept region density. It is **strand-decodable** iff the boundary's two
    regions share a single transcript strand. Returns ``(left, right)`` per-region JointDecodes
    (zero where a side carries no mass — e.g. reference edges).
    """
    ts = np.asarray(region_arrays.ts_class)
    ref_id = np.asarray(region_arrays.ref_id)
    r = ts.shape[0]
    bnd_dec = node_density.boundary_decodable  # boundary (b, b+1) decodable, indexed by left region
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density

    def _side(view, same, ts_other, side_bnd_dec):
        mass = view.mass_unspliced
        pos = view.n_unspliced_pos.astype(np.float64)
        neg = view.n_unspliced_neg.astype(np.float64)
        both_pos = same & (ts == TS_POS) & (ts_other == TS_POS)
        both_neg = same & (ts == TS_NEG) & (ts_other == TS_NEG)
        strand_dec = both_pos | both_neg
        sense = np.where(both_neg, neg, pos)
        antisense = (pos + neg) - sense
        with np.errstate(divide="ignore", invalid="ignore"):
            own = np.where((mass > 0.0) & (eff > 0.0), mass / np.maximum(eff, 1e-12), 0.0)
        density = np.where(side_bnd_dec, own, region_density)  # decodable → own; else swept
        return _joint_per_node(
            mass, view.mass_spliced, sense, antisense, density, pos + neg, eff, strand_dec,
            kappa_rna=kappa_rna, strand_overdispersion=strand_overdispersion,
            confidence=confidence, n_grid=n_grid,
        )

    # LEFT view of r = right side of boundary (r-1, r); neighbour is r-1.
    left_same = np.zeros(r, dtype=bool)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_bnd_dec = np.zeros(r, dtype=bool)
    if r > 1:
        left_same[1:] = ref_id[1:] == ref_id[:-1]
        ts_prev[1:] = ts[:-1]
        left_bnd_dec[1:] = bnd_dec[:-1]
    left = _side(substrate.left, left_same, ts_prev, left_bnd_dec)

    # RIGHT view of r = left side of boundary (r, r+1); neighbour is r+1.
    right_same = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        right_same[:-1] = ref_id[:-1] == ref_id[1:]
        ts_next[:-1] = ts[1:]
    right = _side(substrate.right, right_same, ts_next, bnd_dec)
    return left, right


__all__ = ["JointDecode", "decode_regions", "decode_sides"]
