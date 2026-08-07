"""The aggregate DNA-background reference — a single genome-wide scalar ``(log ρ_bg, σ_bg)``.

The DNA
contamination level is NOT a per-region quantity: a region of effective length ``E`` resolves its DNA rate only
above ``~1/E`` (Fisher information ``ρ·E`` = the expected count), so a faint background — the very case that
matters under strong hybrid capture — is resolvable ONLY by pooling nodes into one aggregate support. This
module measures that pooled scalar; the sweep consumes it as a **one-sided log-floor** (never a scale or
denominator — see the derivation's strong-capture safeguard).

**Substrate (``include_introns``).** Two node sets, both purely intronic/intergenic (never exonic):

* ``include_introns=False`` (default here for *simulation* development) — **INTERGENIC only** (signature ``0``:
  no gene at all). Carries DNA only — no mature *and* no nascent RNA — so the pool is ``0`` for a DNA-free
  library even when nascent is abundant. The sim's ``nrna_present = 0.5×mature`` is wildly unrealistic, which
  would contaminate any intron included; intergenic-only is the safe set *there*.
* ``include_introns=True`` (the **real-data** path) — INTERGENIC + INTRON (non-exonic, ``sig & EXON == 0``).
  Real nascent RNA is *ridiculously sparse*, so a handful of transcribed introns are outliers that barely move
  a pooled mean — while the introns add an ENORMOUS aggregate span (``ΣE``), i.e. resolution, which it would be
  a shame to discard. The optional ``robust_trim_mad`` fence drops the rare nascent-contaminated intron (a high
  per-region density outlier) so it cannot move even the variance. This choice is to be explored **on real
  data**, not the sim.

Under capture the pool is the depleted off-target floor; without capture it is the uniform DNA level
(uniformity assumption; tumor CNV breaks it — see the derivation).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS

__all__ = ["BackgroundReference", "measure_background"]

_EPS = 1.0e-12
_MAD_TO_SIGMA = 1.4826  # MAD → σ for a Gaussian (a mathematical constant, not a tuned knob)
_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
_GENE_BITS = BIT_INTRON_POS | BIT_INTRON_NEG | _EXON_BITS


@dataclass(frozen=True, slots=True)
class BackgroundReference:
    """The pooled DNA-background scalar. ``log_rho_bg`` is ``log(Σg/ΣE)`` over the intergenic pool (natural
    log), ``−inf`` when the pool has zero counts (⇒ the consuming one-sided floor is dormant — "no lower
    bound", correct for a DNA-free or fully-depleted library). ``sigma_bg`` is the Poisson uncertainty on
    ``log ρ_bg`` (``1/√Σg``; ``+inf`` at zero counts). ``n_counts``/``eff_total`` are the aggregate diagnostics."""

    log_rho_bg: float  #: log(Σg/ΣE) over the intergenic pool; −inf when Σg=0 (a DIAGNOSTIC — the aggregate is
    #: still informative there: Σg=0 over a genome-scale ΣE is a *precise near-zero*, resolvable to ~1/ΣE)
    sigma_bg: float  #: √Var(log ρ_bg) ≈ 1/√Σg ; +inf at zero counts
    n_counts: (
        float  #: Σg — pooled DNA counts (the aggregate Poisson observation's count; works at Σg=0)
    )
    eff_total: float  #: ΣE — pooled gDNA effective length (the aggregate support; the source of its precision)
    n_regions: int = (
        0  #: number of pooled NODES — the aggregate cell's EM weight (the population it stands for)
    )
    # log ρ_floor — the DERIVED background-floor location. The pooled rate
    #: ``ln(Σg/ΣE)`` Fisher-blended with the per-region RESOLUTION WALL ``ρ_res = 1/harmmean(E of zero-count
    #: nodes)``: ``ln ρ_floor = (Σg·ln(Σg/ΣE) + n0·ln ρ_res)/(Σg + n0)``. EXACT limits: n0=0 ⇒ ln(Σg/ΣE) (the
    #: resolvable case, byte-identical to the old seed); Σg=0 ⇒ ln ρ_res. Replaces the ``1/ΣE`` seed (3 logs too
    #: low = the confident-FP seed) as the low-density mode the consuming NPMLE injects. ``−inf`` iff no support.
    log_rho_floor: float = -np.inf


def measure_background(
    substrate,
    region_arrays,
    node_eff_len,
    *,
    include_introns: bool = False,
    robust_trim_mad: float | None = None,
) -> BackgroundReference:
    """Pool the non-exonic NODES' contained unspliced counts over their gDNA effective length →
    ``(log ρ_bg, σ_bg)``. A pooled Poisson MLE — precise even when ρ_bg is faint (the aggregate ``Σg`` over a
    genome-scale ``ΣE`` resolves what no single node can), and honestly ``0`` (``σ_bg=∞``) when empty.

    ``node_eff_len`` is ``effective_length.contained_eff_length`` on the gDNA pmf — the count of start
    positions that place a whole molecule inside the node, which is the support the counts are a rate over.

    ``include_introns`` picks the node set (module docstring): ``False`` = intergenic-only (sim-safe);
    ``True`` = intergenic + intron (real-data resolution). ``robust_trim_mad`` (only with introns) is an upper
    MAD fence in log-density: pooled count-bearing nodes above ``median + robust_trim_mad·(MAD·1.4826)`` are
    dropped as nascent-contaminated outliers before the pool. ``None`` ⇒ no trim (plain pool)."""
    sig = np.asarray(region_arrays.signature)
    eff = np.asarray(node_eff_len, dtype=np.float64)
    exclude = _EXON_BITS if include_introns else _GENE_BITS
    pool = ((sig & exclude) == 0) & (eff > _EPS)
    # ⚠ The two columns are GENOME strand, and gDNA is strand-symmetric, so the pool is their sum.
    counts = np.asarray(substrate.node_contained.count, dtype=np.float64).sum(axis=1)

    if robust_trim_mad is not None:
        # Drop the rare high-density (nascent-contaminated) node from the pool. Count-bearing pooled regions
        # only — a count-0 region is never a HIGH outlier and stays (it is the faint-background bulk).
        nz = pool & (counts > _EPS)
        if int(nz.sum()) >= 8:
            ld = np.log(counts[nz]) - np.log(eff[nz])
            med = float(np.median(ld))
            mad = float(np.median(np.abs(ld - med))) * _MAD_TO_SIGMA
            if mad > _EPS:
                hi = med + float(robust_trim_mad) * mad
                dens = counts / np.maximum(eff, _EPS)
                pool = pool & ~((counts > _EPS) & (np.log(np.maximum(dens, _EPS)) > hi))

    sg = float(counts[pool].sum())
    se = float(eff[pool].sum())
    nr = int(pool.sum())
    if se <= _EPS:  # no pooled support at all — genuinely nothing to say
        return BackgroundReference(
            log_rho_bg=-np.inf,
            sigma_bg=np.inf,
            log_rho_floor=-np.inf,
            n_counts=sg,
            eff_total=se,
            n_regions=nr,
        )
    # THE DERIVED FLOOR. The old ``1/ΣE`` seed pools the background as ONE
    # genome-sized region ⇒ claims a resolution n× finer than any single node delivers (~3 logs too low = the
    # confident-FP seed). The honest floor is the density a TYPICAL background region still reads as ~zero — the
    # per-region Poisson resolution wall ``ρ_res = mean(1/E_i)`` over the ZERO-count regions (= 1/harmmean(E)),
    # Fisher-blended with the pooled rate ``Σg/ΣE`` (info Σg vs the n0 unresolved-region votes). EXACT limits:
    # n0=0 ⇒ ln(Σg/ΣE) (byte-identical to the old resolvable case); Σg=0 ⇒ ln ρ_res (the wall).
    zc = pool & (counts <= _EPS)  # the zero-count (unresolved) nodes
    n0 = int(zc.sum())
    s_recip = float((1.0 / eff[zc]).sum()) if n0 > 0 else 0.0  # Σ(1/E_i)
    rho_res = (s_recip / n0) if n0 > 0 else 0.0  # mean(1/E_i) = 1/harmmean(E of zero-count regions)
    if sg > _EPS and n0 > 0 and rho_res > _EPS:
        log_rho_floor = float((sg * (np.log(sg) - np.log(se)) + n0 * np.log(rho_res)) / (sg + n0))
    elif sg > _EPS:  # n0=0 fully resolved ⇒ the pooled rate EXACTLY (old resolvable behaviour)
        log_rho_floor = float(np.log(sg) - np.log(se))
    elif rho_res > _EPS:  # Σg=0 true zero ⇒ the resolution wall
        log_rho_floor = float(np.log(rho_res))
    else:
        log_rho_floor = -np.inf
    # Σg=0 is NOT dormant: it is a precise near-zero (log_rho_bg=−inf is only a diagnostic; the consuming fit
    # uses log_rho_floor as the low-density mode's location).
    log_rho_bg = float(np.log(sg) - np.log(se)) if sg > _EPS else -np.inf
    sigma_bg = float(1.0 / np.sqrt(sg)) if sg > _EPS else np.inf
    return BackgroundReference(
        log_rho_bg=log_rho_bg,
        sigma_bg=sigma_bg,
        log_rho_floor=log_rho_floor,
        n_counts=sg,
        eff_total=se,
        n_regions=nr,
    )
