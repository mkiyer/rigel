"""RNA density model — the per-strand RNA analog of :mod:`density_model` (the gDNA count clue).

The RNA prior of the per-node pie (PLAN v6 §8): each node carries same-strand RNA whose magnitude is
imputed from its chain neighbours' same-strand RNA — principally the **boundary** nodes, which own the
one-sided, motif-stranded **spliced** mass (definitive mature RNA). This module computes, per region and
per strand ``s ∈ {+, −}``, in **RNA density space** (FL-consistent: region ``E_rna[max(0,L−ℓ)]``, per-side
``E_rna[min(ℓ,L_side)]``):

* the region's CURRENT same-strand RNA density (``region_unspliced_s·(1−f_g)/E_rna_region``);
* each flanking boundary side's same-strand RNA density (``(side_unspliced_RNA_s + side_spliced_s)/E_rna_side``).

These feed two consumers off **one** assembly (DRY — :func:`rna_strand_densities`):

* :func:`fit_rna_imputation_varmean` — the RNA reliability ``var~mean`` (the generic
  :func:`variance_model.pair_imputation_points` over the single-strand-exon fit set, both strands pooled);
* :func:`node_rna_density` — the per-node imputed RNA-density **target** ``ρ̂±`` (the apply set: any node that
  can carry strand-``s`` RNA and has a same-strand RNA-observable flank), which :mod:`calibrate` converts to
  the pie fractions ``μ±`` and per-component precisions ``τ_rna±``.

RNA is **never chained globally** (it is spiky / isoform-driven — unlike genomically-smooth gDNA): a node
with no same-strand RNA-observable flank gets no RNA prior (``ρ̂_s = NaN`` → ``τ_rna_s = 0``, degrading to
the strand likelihood + the global foundation). The mature anchor is the flank's spliced (``spliced_s > 0``)
— the prior fires only where a genuine mature reference exists, so it cannot manufacture RNA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .density_model import count_observable_masks
from .run_fill import same_ref_left_right
from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_NEG, TS_POS
from .variance_model import MonotoneVarMean, pair_imputation_points

__all__ = ["StrandDensity", "rna_strand_densities", "fit_rna_imputation_varmean", "node_rna_density"]

_EPS = 1.0e-12


@dataclass(frozen=True, slots=True)
class StrandDensity:
    """Per-region RNA densities for one strand (all length R, RNA density space)."""

    strand: int  # TS_POS or TS_NEG
    region_density: np.ndarray  # the region's CURRENT same-strand RNA density
    left_density: np.ndarray  # the LEFT flanking boundary side's same-strand RNA density
    right_density: np.ndarray  # the RIGHT flanking boundary side's same-strand RNA density
    left_ok: np.ndarray  # bool: left flank is RNA-observable (count-observable side, same-strand spliced>0)
    right_ok: np.ndarray  # bool: right flank is RNA-observable
    fit_eligible: np.ndarray  # bool: single-strand exon dest on this strand (the var~mean LEARNING set)


def _side_rna_frac(side_gdna_frac, side_total, cleaned):
    """Per-side RNA fraction ``1 − side_gdna_frac``, with the strand-cleaned-count fallback where NaN."""
    out = 1.0 - np.where(np.isfinite(side_gdna_frac), side_gdna_frac, 0.0)
    if cleaned is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            fb = np.where(side_total > _EPS, (side_total - cleaned) / side_total, 0.0)
        out = np.where(np.isfinite(side_gdna_frac), out, np.clip(fb, 0.0, 1.0))
    return out


def rna_strand_densities(
    substrate,
    region_arrays,
    region_eff_len_rna,
    rna_boundary_side_eff_len,
    *,
    gdna_frac,
    left_gdna_frac,
    right_gdna_frac,
    cleaned_left=None,
    cleaned_right=None,
) -> dict[int, StrandDensity]:
    """Assemble the per-strand RNA densities (region + both flanking sides) — the shared source.

    ``gdna_frac`` (region, length R) is the region's CURRENT gDNA fraction; ``left/right_gdna_frac`` are the
    flanking boundary sides' gDNA fractions (the per-side ``StrandSplit.gdna_frac``, NaN where unsplit, with
    the ``cleaned_left/right`` strand-cleaned-count fallback). All densities are eff-len-consistent: region
    counts ÷ ``region_eff_len_rna`` (``E_rna[max(0,L−ℓ)]``), side mass ÷ ``rna_boundary_side_eff_len``
    (``E_rna[min(ℓ,L_side)]`` — the per-side density length, NOT ``E_rna[ℓ]``).
    """
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    inv_region = 1.0 / np.maximum(np.asarray(region_eff_len_rna, dtype=np.float64), _EPS)
    inv_side = 1.0 / np.maximum(np.asarray(rna_boundary_side_eff_len, dtype=np.float64), _EPS)

    rna_frac_region = 1.0 - np.clip(np.asarray(gdna_frac, dtype=np.float64), 0.0, 1.0)
    c, left, right = substrate.contained, substrate.left, substrate.right
    left_total = np.asarray(left.n_unspliced_pos + left.n_unspliced_neg, dtype=np.float64)
    right_total = np.asarray(right.n_unspliced_pos + right.n_unspliced_neg, dtype=np.float64)
    cl = None if cleaned_left is None else np.asarray(cleaned_left, dtype=np.float64)
    cr = None if cleaned_right is None else np.asarray(cleaned_right, dtype=np.float64)
    left_rna_frac = _side_rna_frac(np.asarray(left_gdna_frac, dtype=np.float64), left_total, cl)
    right_rna_frac = _side_rna_frac(np.asarray(right_gdna_frac, dtype=np.float64), right_total, cr)

    # Per-side count-observability of region r's flanks (boundary (r−1,r) / (r,r+1), same-ref).
    _rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(r, bool)
    rb = np.zeros(r, bool)
    if r > 1:
        la[1:] = bco[:-1] & ls[1:]
        rb[:-1] = bco[:-1] & rs[:-1]
    has_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0

    # A single-strand region's RNA is ALL on its one strand, regardless of genome orientation: the genome
    # +/− split of the UNSPLICED count is just the κ strand-noise the strand LIKELIHOOD handles, not two RNA
    # species. So the strand-s RNA is the TOTAL unspliced (pos+neg)·(1−f_g) + the motif-oriented spliced
    # (n_spliced_sense = the mature on the junction's own strand = the region's strand here). The two
    # StrandDensity entries therefore share the densities and differ only in the strand eligibility (which
    # selects the matching region set in the fit and the apply). Using a genome-strand-specific count would
    # read ≈0 RNA for a single-strand exon under a reverse-stranded library (κ≈0) — a catastrophic mis-call.
    c_total = np.asarray(c.n_unspliced_pos + c.n_unspliced_neg, dtype=np.float64)
    l_spl = np.asarray(left.n_spliced_sense, dtype=np.float64)
    r_spl = np.asarray(right.n_spliced_sense, dtype=np.float64)
    region_density = c_total * rna_frac_region * inv_region
    left_density = (left_total * left_rna_frac + l_spl) * inv_side
    right_density = (right_total * right_rna_frac + r_spl) * inv_side
    left_ok = la & (l_spl > 0.0)
    right_ok = rb & (r_spl > 0.0)
    return {
        strand: StrandDensity(
            strand=strand,
            region_density=region_density,
            left_density=left_density,
            right_density=right_density,
            left_ok=left_ok,
            right_ok=right_ok,
            fit_eligible=(ts == strand) & has_exon & ~_rco,
        )
        for strand in (TS_POS, TS_NEG)
    }


def fit_rna_imputation_varmean(densities: dict[int, StrandDensity], **kw) -> MonotoneVarMean:
    """The RNA imputation-reliability ``var~mean`` — both strands' node-pair points pooled into one fit.

    The fit set is the **single-strand-exon** dests with a same-strand RNA-observable flank
    (``fit_eligible & *_ok``). Per-strand points are concatenated (NOT the node arrays — so no spurious
    ``+``↔``−`` adjacency); the pooled points fit one monotone curve (Jensen dof=1).
    """
    means, raw = [], []
    for sd in densities.values():
        m, v = pair_imputation_points(
            sd.region_density,
            sd.left_density,
            sd.right_density,
            region_eligible=sd.fit_eligible,
            left_ok=sd.left_ok,
            right_ok=sd.right_ok,
        )
        means.append(m)
        raw.append(v)
    means = np.concatenate(means) if means else np.empty(0)
    raw = np.concatenate(raw) if raw else np.empty(0)
    return MonotoneVarMean.fit(means, raw, dof=np.ones(means.shape[0]), **kw)


def node_rna_density(
    densities: dict[int, StrandDensity], region_arrays
) -> tuple[np.ndarray, np.ndarray]:
    """Per-region imputed same-strand RNA-density **target** ``(ρ̂₊, ρ̂₋)`` — the RNA prior apply set.

    A node gets a strand-``s`` RNA prior where it is **single-strand** ``s`` (``ts == s``) AND has a
    same-strand RNA-observable flank; the target is the mean of the available flanking same-strand RNA
    densities (the boundary→region imputation, the mature flank as the anchor). ``NaN`` where no anchor
    exists (→ ``τ_rna_s = 0``, the prior is off). RNA is never chained globally; no run-fill fallback.

    **AMBIG nodes are excluded in Phase B (deferred to Phase C) — a scope choice, NOT a strand ambiguity.**
    A spliced fragment's strand is its **genomic (motif) strand** — the GT-AG motif is non-palindromic, so
    every junction has a definite strand, knowable independent of the read orientation. So an AMBIG node's two
    flanking junctions DO have definite (possibly opposite) strands, and its ``±`` mature could be assigned
    correctly. The reason to wait is **architectural**: the per-boundary junction strand is owned by the
    bipartite **boundary nodes** (Phase C), where each splice-junction boundary carries its one-strand spliced
    explicitly; the region-centric Phase B sweep does not wire per-flank junction strands. For a single-strand
    region the junction strand = the region's strand (trivially known), so Phase B covers that — the bulk of
    nodes — and an AMBIG node falls to the strand likelihood + global foundation until Phase C.
    """
    ts = np.asarray(region_arrays.strand_class)
    out = []
    for strand in (TS_POS, TS_NEG):
        sd = densities[strand]
        single_strand = ts == strand  # Phase B: single-strand only (AMBIG deferred to Phase C)
        lok = sd.left_ok & single_strand
        rok = sd.right_ok & single_strand
        k = lok.astype(np.float64) + rok.astype(np.float64)
        total = np.where(lok, sd.left_density, 0.0) + np.where(rok, sd.right_density, 0.0)
        out.append(np.where(k > 0.0, total / np.maximum(k, 1.0), np.nan))
    return out[0], out[1]
