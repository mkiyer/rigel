"""The iterative propagation deconvolution — per-node {RNA+, RNA−, gDNA} by solving a node system.

Spec: ``docs/calibration/propagation_implementation_plan.md`` (validated by the prototypes
``scripts/debug/phaseC_propagation*.py``). Each node's contained/crossing unspliced mass is partitioned
into ``(f₊, f₋, f_g)``; the calibration reads ``f_g``. The partition is solved by **propagation from
seeds**:

* **seeds** (self-solvable from their own strand): intergenic (``f_g=1``); single-strand POS/NEG regions
  (the strand splits gDNA vs RNA_s). A seed emits its per-strand **nascent density** ``ρ_nasc_s`` (the RNA
  the strand sees, minus the mature the junctions see), tagged by **enrichment class** (exon vs intron —
  nascent is capture-enriched, so the two classes carry very different densities).
* **propagation**: ``ρ_nasc_s`` is propagated along each strand-``s`` gene body, **precision-weighted**
  (weight = the seed's strand information ``N·(2κ−1)²``) and **enrichment-matched** (exon density to exon
  regions, intron density to intron regions).
* **AMBIG solve** (the #1 resolution, validated): an AMBIG node has ≤2 strands and at least one is
  seedable. Subtract the propagated **seeded-strand** nascent + matures (κ-distributed) → the residual is
  ``gDNA + the unseeded strand's nascent`` → the node's **own strand** resolves it. Where **neither**
  strand is seedable, fall back to the count gDNA density.

**Regions** carry mature (within-exon); **boundary sides** are **mature-free** (a contiguous unspliced
crossing cannot span a junction — the user's physics: spliced crosses one-sided, unspliced both) so their
solve is the same with ``mature=0`` and the per-bp nascent density scaled by the *crossing* eff-length.

Not wired into ``calibrate`` by default (opt-in ``CalibrationConfig.use_propagation``); consolidates
``density_model`` + ``mature_density`` + the ``strand_deconv`` per-node blend + ``run_fill`` once it
replaces ``deconv_regions`` / ``deconv_sides``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mature_density import mature_density
from .signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS, TS_AMBIG, TS_NEG, TS_POS,
)
from .strand_deconv import NodeDeconv, strand_posterior_gdna_frac

__all__ = ["propagate_regions", "propagate"]

_STRAND_BITS = {TS_POS: (BIT_EXON_POS, BIT_INTRON_POS), TS_NEG: (BIT_EXON_NEG, BIT_INTRON_NEG)}
_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class _StrandField:
    """A per-strand propagated nascent-density field + which regions are seedable for that strand."""

    nasc_density: np.ndarray  # float64[R] — propagated per-bp nascent density on this strand
    seedable: np.ndarray  # bool[R] — a seed for this strand is reachable in the region's gene-body run


@dataclass(frozen=True, slots=True)
class _Fields:
    """The per-strand nascent fields + matures shared by the region and side solves."""

    nasc_density_pos: np.ndarray
    nasc_density_neg: np.ndarray
    seedable_pos: np.ndarray
    seedable_neg: np.ndarray
    mature_pos: np.ndarray
    mature_neg: np.ndarray


def _body_run_segments(in_body: np.ndarray, ref_id: np.ndarray) -> np.ndarray:
    """Segment id constant within each contiguous in-body run (per reference) — the propagation reach."""
    in_body = np.asarray(in_body, dtype=bool)
    ref = np.asarray(ref_id)
    r = in_body.shape[0]
    if r == 0:
        return np.zeros(0, dtype=np.int64)
    brk = np.zeros(r, dtype=bool)
    brk[1:] = (in_body[1:] != in_body[:-1]) | (ref[1:] != ref[:-1])
    return np.cumsum(brk)


def _propagate_strand(
    strand, *, sig, ref_id, strand_class, u_pos, u_neg, mature_s, e_rna, kappa, od_g, od_r, n_grid,
) -> _StrandField:
    """Propagate strand-``s`` nascent density from single-strand-``s`` seeds, precision-weighted + per class.

    A single-strand-``s`` region is a seed: its own strand posterior gives ``f_g``, so ``RNA_s=(1−f_g)·U``,
    ``nascent_s = RNA_s − mature_s``, density ``ρ = nascent_s / E_rna``, precision ``= N·(2κ−1)²``. Within
    each strand-``s`` gene-body run and **enrichment class** (exon / intron) the field is the
    precision-weighted mean of that class's seeds; a region with no same-class seed in its run is **not
    seedable** for ``s`` (its nascent is left to the AMBIG node-strand solve).
    """
    exon_bit, intron_bit = _STRAND_BITS[strand]
    r = sig.shape[0]
    sense = u_neg if strand == TS_NEG else u_pos
    antisense = u_pos if strand == TS_NEG else u_neg
    n = sense + antisense
    discrim = (2.0 * kappa - 1.0) ** 2

    is_seed = (strand_class == strand) & (n > 0.0)
    nasc_density = np.zeros(r, dtype=np.float64)
    precision = np.zeros(r, dtype=np.float64)
    idx = np.flatnonzero(is_seed)
    if idx.size:
        g_q, _ = strand_posterior_gdna_frac(
            sense[idx], antisense[idx], kappa,
            gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r, n_grid=n_grid,
        )
        rna_s = (1.0 - g_q) * n[idx]
        nasc = np.maximum(rna_s - mature_s[idx], 0.0)
        nasc_density[idx] = np.where(e_rna[idx] > _EPS, nasc / np.maximum(e_rna[idx], _EPS), 0.0)
        precision[idx] = n[idx] * discrim

    in_body = (sig & (exon_bit | intron_bit)) != 0
    run = _body_run_segments(in_body, ref_id)
    is_exon = (sig & exon_bit) != 0
    is_intron = (sig & intron_bit) != 0

    field = np.zeros(r, dtype=np.float64)
    seedable = np.zeros(r, dtype=bool)
    for cls_mask in (is_exon, is_intron):
        seed_here = is_seed & cls_mask & (precision > 0.0)
        for seg in np.unique(run[in_body & cls_mask]):
            in_seg = (run == seg) & cls_mask & in_body
            s_in = seed_here & in_seg
            wsum = float(precision[s_in].sum())
            if wsum > 0.0:
                wmean = float((nasc_density[s_in] * precision[s_in]).sum() / wsum)
                field[in_seg] = wmean
                seedable[in_seg] = True
    return _StrandField(nasc_density=field, seedable=seedable)


def _compute_fields(substrate, region_arrays, e_rna, rna_fl_mean, kappa, od_g, od_r, n_grid) -> _Fields:
    """Compute the per-strand nascent-density fields + matures (shared by the region and side solves)."""
    ts = np.asarray(region_arrays.strand_class)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    md = mature_density(substrate, region_arrays, e_rna, rna_fl_mean)
    common = dict(sig=sig, ref_id=ref_id, strand_class=ts, u_pos=u_pos, u_neg=u_neg, e_rna=e_rna,
                  kappa=float(kappa), od_g=od_g, od_r=od_r, n_grid=n_grid)
    fpos = _propagate_strand(TS_POS, mature_s=md.mature_pos, **common)
    fneg = _propagate_strand(TS_NEG, mature_s=md.mature_neg, **common)
    return _Fields(fpos.nasc_density, fneg.nasc_density, fpos.seedable, fneg.seedable,
                   md.mature_pos, md.mature_neg)


def _solve_regions(substrate, region_arrays, fields, e_rna, kappa, od_g, od_r, count_gdna_frac, n_grid):
    """Per-region f_g: intergenic = 1; single-strand = strand posterior; AMBIG = the propagation solve."""
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(c.mass_spliced, dtype=np.float64)
    nasc_pos = np.where(fields.seedable_pos, fields.nasc_density_pos * e_rna, 0.0)
    nasc_neg = np.where(fields.seedable_neg, fields.nasc_density_neg * e_rna, 0.0)

    f_g = np.full(U.shape, np.nan, dtype=np.float64)
    active = U > 0.0
    f_g[active & (ts != TS_POS) & (ts != TS_NEG) & (ts != TS_AMBIG)] = 1.0  # intergenic / no transcript
    for strand, sense_arr in ((TS_POS, u_pos), (TS_NEG, u_neg)):
        idx = np.flatnonzero(active & (ts == strand))
        if idx.size:
            g_q, _ = strand_posterior_gdna_frac(
                sense_arr[idx], (U - sense_arr)[idx], kappa,
                gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r, n_grid=n_grid,
            )
            f_g[idx] = g_q
    amb = active & (ts == TS_AMBIG)
    if amb.any():
        f_g[amb] = _solve_ambig(
            u_pos[amb], u_neg[amb], kappa, nasc_pos[amb], nasc_neg[amb],
            fields.seedable_pos[amb], fields.seedable_neg[amb],
            fields.mature_pos[amb], fields.mature_neg[amb], count_gdna_frac[amb],
        )
    f_g = np.where(active, np.clip(np.nan_to_num(f_g, nan=0.0), 0.0, 1.0), 0.0)
    return NodeDeconv(gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
                      gdna_frac=f_g, gdna_frac_var=np.zeros_like(f_g))


def _solve_sides(substrate, region_arrays, fields, rna_fl_mean, kappa, od_g, od_r, count_gdna_frac, n_grid):
    """Boundary sides (left, right): mature-free crossing → {gDNA, nascent±}.

    A side's unspliced crossing is gDNA + nascent (mature crosses spliced, one-sided). **Strand-observable**
    sides (a single consistent sense across the boundary) use their **own strand** directly — exactly like a
    single-strand region, and crucially *not* the propagated field (which carries the seed's small spurious
    nascent and would siphon gDNA→RNA at every observable side). Only **AMBIG** (strand-unobservable) sides
    use the propagation subtract-residual solve (mature=0, the per-bp nascent density × the **crossing**
    eff-length ``rna_fl_mean``). Region ``r``'s field applies to both of ``r``'s side views.
    """
    from .density_model import count_observable_masks
    from .strand_deconv import _left_right_neighbors, _side_strand_orientation

    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    _, bobs = count_observable_masks(sig, ref_id)
    l_same, ts_prev, _l_obs, r_same, ts_next, _r_obs = _left_right_neighbors(ts, ref_id, bobs)
    nasc_pos = np.where(fields.seedable_pos, fields.nasc_density_pos * rna_fl_mean, 0.0)
    nasc_neg = np.where(fields.seedable_neg, fields.nasc_density_neg * rna_fl_mean, 0.0)
    zeros = np.zeros_like(nasc_pos)

    def _solve(view, same, ts_other) -> NodeDeconv:
        u_pos = view.n_unspliced_pos.astype(np.float64)
        u_neg = view.n_unspliced_neg.astype(np.float64)
        U = u_pos + u_neg
        mass_unspl = np.asarray(view.mass_unspliced, dtype=np.float64)
        mass_spliced = np.asarray(view.mass_spliced, dtype=np.float64)
        sense, antisense, _n, observable = _side_strand_orientation(view, same, ts, ts_other)
        frac = np.zeros_like(U)
        # observable sides: their OWN strand resolves gDNA vs RNA (no propagated field — avoids the siphon).
        obs = np.flatnonzero(observable & (U > 0.0))
        if obs.size:
            g_q, _ = strand_posterior_gdna_frac(
                sense[obs], antisense[obs], kappa,
                gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r, n_grid=n_grid,
            )
            frac[obs] = g_q
        # AMBIG (unobservable) sides: the propagation subtract-residual solve (mature-free).
        amb = (~observable) & (U > 0.0)
        if amb.any():
            frac[amb] = _solve_ambig(
                u_pos[amb], u_neg[amb], kappa, nasc_pos[amb], nasc_neg[amb],
                fields.seedable_pos[amb], fields.seedable_neg[amb], zeros[amb], zeros[amb],
                count_gdna_frac[amb],
            )
        frac = np.where(U > 0.0, frac, 0.0)
        return NodeDeconv(gdna_mass=frac * mass_unspl, rna_mass=(1.0 - frac) * mass_unspl + mass_spliced,
                          gdna_frac=frac, gdna_frac_var=np.zeros_like(frac))

    return _solve(substrate.left, l_same, ts_prev), _solve(substrate.right, r_same, ts_next)


def propagate_regions(
    substrate, region_arrays, *, rna_region_eff_len, rna_fl_mean, rna_sense_frac,
    gdna_strand_overdispersion, rna_strand_overdispersion, count_gdna_frac, n_grid=200,
) -> NodeDeconv:
    """Per-region gDNA/RNA deconvolution by propagation (regions only; see :func:`propagate` for sides)."""
    e_rna = np.asarray(rna_region_eff_len, dtype=np.float64)
    fields = _compute_fields(substrate, region_arrays, e_rna, rna_fl_mean, rna_sense_frac,
                             gdna_strand_overdispersion, rna_strand_overdispersion, n_grid)
    return _solve_regions(substrate, region_arrays, fields, e_rna, rna_sense_frac,
                          gdna_strand_overdispersion, rna_strand_overdispersion, count_gdna_frac, n_grid)


def propagate(
    substrate, region_arrays, *, rna_region_eff_len, rna_fl_mean, rna_sense_frac,
    gdna_strand_overdispersion, rna_strand_overdispersion, count_gdna_frac, n_grid=200,
):
    """Propagation deconvolution for the region + both boundary sides → ``(regions, left, right)``.

    Computes the per-strand nascent fields once and applies the region solve (with mature) and the
    side solve (mature-free, crossing eff-length). Replaces ``deconv_regions`` + ``deconv_sides``.
    """
    e_rna = np.asarray(rna_region_eff_len, dtype=np.float64)
    fields = _compute_fields(substrate, region_arrays, e_rna, rna_fl_mean, rna_sense_frac,
                             gdna_strand_overdispersion, rna_strand_overdispersion, n_grid)
    regions = _solve_regions(substrate, region_arrays, fields, e_rna, rna_sense_frac,
                             gdna_strand_overdispersion, rna_strand_overdispersion, count_gdna_frac, n_grid)
    left, right = _solve_sides(substrate, region_arrays, fields, rna_fl_mean, rna_sense_frac,
                               gdna_strand_overdispersion, rna_strand_overdispersion, count_gdna_frac,
                               n_grid)
    return regions, left, right


def _solve_ambig(u_pos, u_neg, kappa, nasc_pos, nasc_neg, seedable_pos, seedable_neg,
                 m_pos, m_neg, count_gdna_frac):
    """f_g: subtract seeded nascent + matures, then the node strand resolves the unseeded residual.

    Seeded-strand RNA (nascent + mature) is distributed onto the observed counts by κ and subtracted; the
    residual = gDNA + the UNSEEDED strand's nascent. With one unseeded strand the residual tilt gives that
    nascent (gDNA = residual − |tilt|); with both seeded the residual is gDNA; with neither seedable we fall
    back to the count gDNA fraction.
    """
    U = u_pos + u_neg
    k = float(kappa)
    rna_pos_s = np.where(seedable_pos, nasc_pos, 0.0) + m_pos
    rna_neg_s = np.where(seedable_neg, nasc_neg, 0.0) + m_neg
    seeded_pos = k * rna_pos_s + (1.0 - k) * rna_neg_s
    seeded_neg = (1.0 - k) * rna_pos_s + k * rna_neg_s
    up2 = np.maximum(u_pos - seeded_pos, 0.0)
    un2 = np.maximum(u_neg - seeded_neg, 0.0)
    resid = up2 + un2

    both = seedable_pos & seedable_neg
    neither = ~seedable_pos & ~seedable_neg
    denom = 2.0 * k - 1.0
    nasc_resid = np.where(np.abs(denom) > _EPS, np.abs(up2 - un2) / np.abs(denom), 0.0)
    gdna = np.clip(resid - nasc_resid, 0.0, None)
    f_g = np.where(U > 0.0, gdna / np.maximum(U, _EPS), 0.0)
    f_g = np.where(both, np.clip(resid / np.maximum(U, _EPS), 0.0, 1.0), f_g)
    f_g = np.where(neither, np.clip(count_gdna_frac, 0.0, 1.0), f_g)
    # SAFETY FLOOR: the count gDNA density is a (capture-depleted) LOWER bound on the true gDNA, so the
    # propagated-nascent subtraction must not push f_g below it. This stops the catastrophic over-subtraction
    # when the seed nascent is spurious (mature_density misses the mature of single-exon/no-junction regions
    # → the mature RNA is mislabelled nascent → over-subtracted). Until the nascent source is made
    # mature-free (boundary-node nascent), this guarantees propagation never regresses below the count.
    return np.clip(np.maximum(f_g, np.clip(count_gdna_frac, 0.0, 1.0)), 0.0, 1.0)
