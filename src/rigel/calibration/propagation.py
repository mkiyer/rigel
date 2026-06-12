"""The iterative propagation deconvolution — per-region {RNA+, RNA−, gDNA} by solving a node system.

Spec: ``docs/calibration/propagation_implementation_plan.md`` (validated by the prototypes
``scripts/debug/phaseC_propagation*.py``). Each region's contained-unspliced mass is partitioned into
``(f₊, f₋, f_g)``; the calibration reads ``f_g``. The partition is solved by **propagation from seeds**:

* **seeds** (self-solvable from their own strand): intergenic (``f_g=1``); single-strand POS/NEG regions
  (the strand splits gDNA vs RNA_s). A seed emits its per-strand **nascent density** ``ρ_nasc_s`` (the RNA
  the strand sees, minus the mature the junctions see), tagged by **enrichment class** (exon vs intron —
  nascent is capture-enriched, so the two classes carry very different densities).
* **propagation**: ``ρ_nasc_s`` is propagated along each strand-``s`` gene body, **precision-weighted**
  (weight = the seed's strand information ``N·(2κ−1)²`` — a tiny/uninformative seed cannot poison the
  field) and **enrichment-matched** (exon density to exon regions, intron density to intron regions).
* **AMBIG solve** (the #1 resolution, validated): an AMBIG region has ≤2 strands and at least one is
  seedable. Subtract the propagated **seeded-strand** nascent + **both** matures (κ-distributed) → the
  residual is ``gDNA + the unseeded strand's nascent`` → the region's **own strand** resolves it (gDNA
  symmetric, the remaining nascent tilted). Propagation supplies one strand; the node strand finishes the
  other. Where **neither** strand is seedable (pathological), fall back to the count gDNA density.

This is **not** wired into ``calibrate`` yet — built + validated standalone first (Phase-A discipline). It
consolidates ``density_model`` + ``mature_density`` + the ``strand_deconv`` per-node blend + ``run_fill``;
those are retired once this replaces ``deconv_regions``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mature_density import mature_density
from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS, TS_NEG, TS_POS
from .strand_deconv import NodeDeconv, strand_posterior_gdna_frac

__all__ = ["propagate_regions"]

# Per-strand (exon bit, intron bit) — a region's enrichment class for strand s.
_STRAND_BITS = {TS_POS: (BIT_EXON_POS, BIT_INTRON_POS), TS_NEG: (BIT_EXON_NEG, BIT_INTRON_NEG)}
_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class _StrandField:
    """A per-strand propagated nascent-density field + which regions are seedable for that strand."""

    nasc_density: np.ndarray  # float64[R] — propagated per-bp nascent density on this strand
    seedable: np.ndarray  # bool[R] — a seed for this strand is reachable in the region's gene-body run


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
    strand: int,
    *,
    sig: np.ndarray,
    ref_id: np.ndarray,
    strand_class: np.ndarray,
    u_pos: np.ndarray,
    u_neg: np.ndarray,
    mature_s: np.ndarray,
    e_rna: np.ndarray,
    kappa: float,
    od_g: float,
    od_r: float,
    n_grid: int,
) -> _StrandField:
    """Propagate strand-``s`` nascent density from single-strand-``s`` seeds, precision-weighted + per class.

    A single-strand-``s`` region is a seed: its own strand posterior gives ``f_g``, so ``RNA_s=(1−f_g)·U``,
    ``nascent_s = RNA_s − mature_s``, density ``ρ = nascent_s / E_rna``, precision ``= N·(2κ−1)²``. Within
    each strand-``s`` gene-body run (contiguous ``E_s|I_s``) and **enrichment class** (exon / intron), the
    field is the precision-weighted mean of that class's seeds; a region with no same-class seed in its run
    is **not seedable** for ``s`` (its nascent is left to the AMBIG node-strand solve).
    """
    exon_bit, intron_bit = _STRAND_BITS[strand]
    r = sig.shape[0]
    sense = u_neg if strand == TS_NEG else u_pos
    antisense = u_pos if strand == TS_NEG else u_neg
    n = sense + antisense
    discrim = (2.0 * kappa - 1.0) ** 2

    # seed nascent density + precision on single-strand-s regions (the strand resolves them directly).
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
        precision[idx] = n[idx] * discrim  # strand information of the seed

    # propagate within (gene-body run × enrichment class): precision-weighted mean of same-class seeds.
    in_body = (sig & (exon_bit | intron_bit)) != 0
    run = _body_run_segments(in_body, ref_id)
    is_exon = (sig & exon_bit) != 0
    is_intron = (sig & intron_bit) != 0

    field = np.zeros(r, dtype=np.float64)
    seedable = np.zeros(r, dtype=bool)
    for cls_mask in (is_exon, is_intron):
        # group seeds of this class by run; assign each region of this class in the run the wtd mean.
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


def propagate_regions(
    substrate,
    region_arrays,
    *,
    rna_region_eff_len: np.ndarray,
    rna_fl_mean: float,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float,
    rna_strand_overdispersion: float,
    count_gdna_frac: np.ndarray,
    n_grid: int = 200,
) -> NodeDeconv:
    """Per-region gDNA/RNA deconvolution by propagation (see the module docstring).

    ``count_gdna_frac`` is the fallback gDNA fraction (the count density imputation) for regions where
    neither strand is seedable. Returns a :class:`NodeDeconv` (``gdna_mass``, ``rna_mass``, ``gdna_frac``,
    ``gdna_frac_var``) matching ``deconv_regions`` so it can replace it.
    """
    ts = np.asarray(region_arrays.strand_class)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    kappa = float(rna_sense_frac)
    e_rna = np.asarray(rna_region_eff_len, dtype=np.float64)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(c.mass_spliced, dtype=np.float64)

    md = mature_density(substrate, region_arrays, e_rna, rna_fl_mean)
    m_pos, m_neg = md.mature_pos, md.mature_neg

    # propagate each strand's nascent density (precision-weighted, enrichment-matched).
    common = dict(sig=sig, ref_id=ref_id, strand_class=ts, u_pos=u_pos, u_neg=u_neg, e_rna=e_rna,
                  kappa=kappa, od_g=gdna_strand_overdispersion, od_r=rna_strand_overdispersion,
                  n_grid=n_grid)
    fpos = _propagate_strand(TS_POS, mature_s=m_pos, **common)
    fneg = _propagate_strand(TS_NEG, mature_s=m_neg, **common)

    # nascent counts (density × eff-len) where seedable, else 0 (left to the node-strand residual solve).
    nasc_pos = np.where(fpos.seedable, fpos.nasc_density * e_rna, 0.0)
    nasc_neg = np.where(fneg.seedable, fneg.nasc_density * e_rna, 0.0)

    # ---- the per-region solve ----
    f_g = np.full(U.shape, np.nan, dtype=np.float64)
    active = U > 0.0

    # (a) intergenic / no transcript: all gDNA.
    none_mask = active & (ts != TS_POS) & (ts != TS_NEG) & ~_is_ambig(ts)
    f_g[none_mask] = 1.0

    # (b) single-strand seeds: the strand posterior gives f_g directly.
    for strand, sense_arr in ((TS_POS, u_pos), (TS_NEG, u_neg)):
        m = active & (ts == strand)
        idx = np.flatnonzero(m)
        if idx.size:
            g_q, _ = strand_posterior_gdna_frac(
                sense_arr[idx], (U - sense_arr)[idx], kappa,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion, n_grid=n_grid,
            )
            f_g[idx] = g_q

    # (c) AMBIG: subtract seeded-strand nascent + both matures -> residual -> node strand resolves.
    amb = active & _is_ambig(ts)
    if amb.any():
        f_g[amb] = _solve_ambig(
            u_pos[amb], u_neg[amb], kappa,
            nasc_pos[amb], nasc_neg[amb], fpos.seedable[amb], fneg.seedable[amb],
            m_pos[amb], m_neg[amb], count_gdna_frac[amb],
        )

    f_g = np.where(active, np.clip(np.nan_to_num(f_g, nan=0.0), 0.0, 1.0), 0.0)
    gdna = f_g * mass_unspl
    rna = (1.0 - f_g) * mass_unspl + mass_spliced
    return NodeDeconv(gdna_mass=gdna, rna_mass=rna, gdna_frac=f_g,
                      gdna_frac_var=np.zeros_like(f_g))


def _is_ambig(ts: np.ndarray) -> np.ndarray:
    from .signature import TS_AMBIG
    return np.asarray(ts) == TS_AMBIG


def _solve_ambig(u_pos, u_neg, kappa, nasc_pos, nasc_neg, seedable_pos, seedable_neg,
                 m_pos, m_neg, count_gdna_frac):
    """AMBIG f_g: subtract seeded nascent + matures, then the node strand resolves the unseeded residual.

    Seeded-strand RNA (nascent + mature) is distributed onto the observed counts by κ and subtracted; the
    residual = gDNA + the UNSEEDED strand's nascent. With one unseeded strand the residual tilt gives that
    nascent (gDNA = residual − |tilt|); with both seeded the residual is gDNA; with neither seedable we fall
    back to the count gDNA fraction.
    """
    U = u_pos + u_neg
    k = float(kappa)
    # seeded RNA per strand = nascent (if seedable) + mature; distribute by κ to observed pos/neg.
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
    # one unseeded strand: residual tilt is that strand's nascent; gDNA = residual − |nascent|.
    nasc_resid = np.where(np.abs(denom) > _EPS, np.abs(up2 - un2) / np.abs(denom), 0.0)
    gdna = np.clip(resid - nasc_resid, 0.0, None)
    f_g = np.where(U > 0.0, gdna / np.maximum(U, _EPS), 0.0)
    f_g = np.where(both, np.clip(resid / np.maximum(U, _EPS), 0.0, 1.0), f_g)  # residual is all gDNA
    f_g = np.where(neither, np.clip(count_gdna_frac, 0.0, 1.0), f_g)  # fallback: count density
    return np.clip(f_g, 0.0, 1.0)
