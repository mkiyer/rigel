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


def decode_boundaries(
    substrate, region_arrays, node_density, mu_fl, *,
    kappa_rna, strand_overdispersion=0.0, confidence=0.0, n_grid=200,
) -> JointDecode:
    """Deconvolve each internal boundary's crossing mass (a node) into gDNA / RNA.

    A boundary node ``b`` (between regions ``b`` and ``b+1`` of the same reference) carries the
    crossing mass/flux from both sides. Strand-decodable only when both sides share a single
    transcript strand (POS-POS or NEG-NEG); otherwise count-only. Non-boundary rows (a
    reference's last region) are emitted with zero mass.
    """
    ts = np.asarray(region_arrays.ts_class)
    ref_id = np.asarray(region_arrays.ref_id)
    r = ts.shape[0]
    same = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    dens_next = np.zeros(r, dtype=np.float64)
    cev_next = np.zeros(r, dtype=np.float64)
    if r > 1:
        same[:-1] = ref_id[:-1] == ref_id[1:]
        ts_next[:-1] = ts[1:]
        dens_next[:-1] = node_density.density[1:]
        cev_next[:-1] = node_density.count_evidence[1:]

    lft, rgt = substrate.left, substrate.right
    mass = np.zeros(r)
    spliced = np.zeros(r)
    pos = np.zeros(r)
    neg = np.zeros(r)
    if r > 1:
        mass[:-1] = rgt.mass_unspliced[:-1] + lft.mass_unspliced[1:]
        spliced[:-1] = rgt.mass_spliced[:-1] + lft.mass_spliced[1:]
        pos[:-1] = rgt.n_unspliced_pos[:-1].astype(np.float64) + lft.n_unspliced_pos[1:].astype(np.float64)
        neg[:-1] = rgt.n_unspliced_neg[:-1].astype(np.float64) + lft.n_unspliced_neg[1:].astype(np.float64)
    mass = np.where(same, mass, 0.0)  # non-boundary rows carry no node

    both_pos = same & (ts == TS_POS) & (ts_next == TS_POS)
    both_neg = same & (ts == TS_NEG) & (ts_next == TS_NEG)
    strand_dec = both_pos | both_neg
    sense = np.where(both_neg, neg, pos)
    antisense = (pos + neg) - sense
    # A count-decodable boundary's crossings are gDNA-representative *by signature* (no exon
    # bit shared ⇒ no unspliced mature RNA crosses), so π_count → 1 there (= own crossing
    # density · μ_FL / M = 1); the strand then carves out any nascent RNA. Non-decodable
    # boundaries (an exon spans them) instead use the swept neighbour density.
    with np.errstate(divide="ignore", invalid="ignore"):
        own_density = np.where(mass > 0.0, mass / mu_fl, 0.0)
    swept_density = np.where(same, 0.5 * (node_density.density + dens_next), 0.0)
    density_b = np.where(node_density.boundary_decodable, own_density, swept_density)
    count_b = pos + neg  # the boundary's own discrete crossing count (statistical power)
    eff = np.full(r, mu_fl, dtype=np.float64)
    return _joint_per_node(
        mass, spliced, sense, antisense, density_b, count_b, eff, strand_dec,
        kappa_rna=kappa_rna, strand_overdispersion=strand_overdispersion,
        confidence=confidence, n_grid=n_grid,
    )


__all__ = ["JointDecode", "decode_regions", "decode_boundaries"]
