"""Mature-RNA contained-unspliced density imputation — the RNA mirror of ``density_model``.

Per region and **per strand**, predict the contained-unspliced **mature** count from the spliced crossings
at the bounding eligible splice-junction boundaries. Structurally identical to ``density_model``'s gDNA
density imputation — anchor at observable boundaries, average the bounding anchors, run-fill interiors —
with three swaps:

* anchor boundary: an **eligible splice junction** (per strand) instead of a count-observable boundary;
* crossing count: the **spliced** flux (the junction-strand's mature) instead of the unspliced;
* effective lengths: the **RNA** FL — ``fl_mean_rna`` for the crossing, ``region_eff_length(L, rna_pmf)``
  (``E_mu``) for the contained mature.

A junction's spliced flux is an unbiased mature-rate estimator (``flux / fl_mean_rna = m̂``), validated in
``scripts/debug/phaseA_issueD_mature_density.py`` (recovers the contained-unspliced mature within ~5% across
exon sizes 400–1500 bp; predicts ~0 for sub-FL short exons where ``E_mu → 0``). The averaging of the two
bounding-junction densities reproduces the validated ``÷(n_junctions·fl_mean)`` divisor — no new geometry.

Run-fill stays **within a strand's exon**: an interior sub-region of a long exon (split by an overlapping
opposite-strand gene, so its own boundaries are exon↔exon, not exon↔intron, and carry no junction) inherits
the mature density carried across the *contiguous same-strand-exon run* — it never bridges across an intron
into a different exon. Regions with no junction anchor reachable in their run get ``0`` (no mature evidence
→ predict no mature; the safe direction — never invent mature to subtract). See
``docs/calibration/phaseA_mature_imputation_plan.md`` and ``deconvolution_roadmap.md`` (signal #2).

**Phase A scope:** this module *predicts* ``M_+`` / ``M_−``. The strand-aware subtraction (κ split) and the
residual deconvolution are Phases B/C.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .run_fill import runfill_bidirectional
from .signature import BIT_EXON_NEG, BIT_EXON_POS, TS_NEG, TS_POS
from .splice_junction import splice_junction_eligibility

__all__ = ["MatureDensity", "mature_density"]

_STRAND_EXON_BIT = {TS_POS: BIT_EXON_POS, TS_NEG: BIT_EXON_NEG}


@dataclass(frozen=True, slots=True)
class MatureDensity:
    """Per-region per-strand predicted contained-unspliced mature (the RNA mirror of ``NodeDensity``)."""

    m_pos: np.ndarray  # float64[R] — + strand mature density (rate, fragments / effective bp)
    m_neg: np.ndarray  # float64[R] — − strand mature density
    mature_pos: np.ndarray  # float64[R] — predicted contained-unspliced + mature count = m_pos · E_mu
    mature_neg: np.ndarray  # float64[R] — predicted contained-unspliced − mature count = m_neg · E_mu
    anchor_pos: np.ndarray  # bool[R] — region is the exon side of an eligible + splice junction
    anchor_neg: np.ndarray  # bool[R] — region is the exon side of an eligible − splice junction


def _exon_run_segments(is_exon: np.ndarray, ref_id: np.ndarray) -> np.ndarray:
    """Segment id constant within each contiguous same-strand-exon run (and per reference).

    Returns an int array usable as the ``ref_id`` argument to :func:`runfill_bidirectional`, so the carry
    fills only *within* a strand's exon (a split-up long exon) and never bridges an intron / a reference
    edge / a non-exon region. The boundary increments wherever ``is_exon`` flips or the reference changes;
    so each maximal exon run is one segment and every non-exon region is its own (carry-isolated) segment.
    """
    is_exon = np.asarray(is_exon, dtype=bool)
    ref = np.asarray(ref_id)
    r = is_exon.shape[0]
    if r == 0:
        return np.zeros(0, dtype=np.int64)
    brk = np.zeros(r, dtype=bool)
    brk[1:] = (is_exon[1:] != is_exon[:-1]) | (ref[1:] != ref[:-1])
    return np.cumsum(brk)


def mature_density(
    substrate,
    region_arrays,
    rna_region_eff_len: np.ndarray,
    fl_mean_rna: float,
) -> MatureDensity:
    """Per-region per-strand predicted contained-unspliced mature count (see module docstring).

    ``rna_region_eff_len`` = ``region_eff_length(region_size_bp, rna_pmf)`` (``E_mu``); ``fl_mean_rna`` =
    ``boundary_eff_length(rna_pmf)`` (the RNA crossing eff-length). Returns a :class:`MatureDensity`.
    """
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    e_mu = np.asarray(rna_region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    inv_fl = 1.0 / fl_mean_rna if fl_mean_rna > 0.0 else 0.0

    # Spliced flux into region r from each side (motif-oriented total = the junction-strand's mature, per
    # Issue A: one boundary aggregates ALL spliced crossings into the region). The junction's strand is
    # taken from the eligibility predicate, so the side total is attributed to that strand.
    left_spl = (substrate.left.n_spliced_sense + substrate.left.n_spliced_antisense).astype(np.float64)
    right_spl = (substrate.right.n_spliced_sense + substrate.right.n_spliced_antisense).astype(np.float64)

    def _field(strand: int) -> tuple[np.ndarray, np.ndarray]:
        """Mature density + anchor mask for one strand: anchor at eligible junctions, run-fill the run."""
        dens = np.full(r, np.nan, dtype=np.float64)
        anchor = np.zeros(r, dtype=bool)
        for i in range(r):
            ests: list[float] = []
            if i > 0 and ref_id[i] == ref_id[i - 1]:
                el = splice_junction_eligibility(int(sig[i - 1]), int(sig[i]))
                if el is not None and el.exon_side == "R" and strand in el.strands:
                    ests.append(left_spl[i] * inv_fl)
            if i < r - 1 and ref_id[i] == ref_id[i + 1]:
                el = splice_junction_eligibility(int(sig[i]), int(sig[i + 1]))
                if el is not None and el.exon_side == "L" and strand in el.strands:
                    ests.append(right_spl[i] * inv_fl)
            if ests:
                dens[i] = float(np.mean(ests))
                anchor[i] = True
        # Carry the anchored density across the contiguous same-strand-exon run (uniform mature over an
        # exon); never across an intron. Unreached regions (no junction in the run) → 0 (predict no mature).
        is_exon = (sig & _STRAND_EXON_BIT[strand]) != 0
        seg = _exon_run_segments(is_exon, ref_id)
        dens = runfill_bidirectional(dens, seg)
        dens = np.where(np.isnan(dens) | ~is_exon, 0.0, dens)
        return dens, anchor

    m_pos, a_pos = _field(TS_POS)
    m_neg, a_neg = _field(TS_NEG)
    return MatureDensity(
        m_pos=m_pos,
        m_neg=m_neg,
        mature_pos=m_pos * e_mu,
        mature_neg=m_neg * e_mu,
        anchor_pos=a_pos,
        anchor_neg=a_neg,
    )
