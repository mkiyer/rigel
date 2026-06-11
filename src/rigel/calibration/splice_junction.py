"""Splice-junction gDNA-fraction — eligibility predicate + boundary fraction.

Phase 4-mean (design: ``docs/calibration/count_mean_bias_design.md``). A non-observable exon region's
gDNA mass is biased ~2× low under hybrid capture when imputed from an *absolute* boundary density (the
boundary edge is depleted relative to the probe-enriched interior). The fix: at a **true splice-junction
boundary** the crossing fragments split cleanly *by construction* —

  * crossing **unspliced** = gDNA (+ nascent) — mature mRNA cannot contiguously span an exon→intron edge;
  * crossing **spliced**   = mature mRNA       — only mature carries the junction —

so the boundary yields a ground-truth-free gDNA *fraction*, which partitions the neighbouring region's
own (enriched) total. The hard part is **eligibility**: which adjacent-region 4-bit signature pairs are
genuine, unambiguous splice junctions. That predicate lives here; it is pure signature logic with **no
tunable parameters**. Realizability (which pairs a valid transcript layout can produce) is a *separate*,
empirically-tested property (``tests/calibration/test_splice_junction_realized.py``) — the predicate
classifies any pair it is handed; the region partition only ever feeds it realizable ones.

Where no boundary is eligible, calibration falls back to the absolute-density imputation
(``density_model``), which needs only an observable density and is always available.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_NEG,
    TS_POS,
    validate_signature,
)

__all__ = [
    "SpliceJunction",
    "splice_junction_eligibility",
    "boundary_gdna_fraction",
    "density_frac_to_count_frac",
    "region_splice_gdna_frac",
]

# (transcript-strand class, exon bit, intron bit) for each strand — the two annotation layers a
# signature carries. The predicate treats the strands independently and then requires them to agree
# (the matched-set rule below).
_STRAND_LAYERS: tuple[tuple[int, int, int], ...] = (
    (TS_POS, BIT_EXON_POS, BIT_INTRON_POS),
    (TS_NEG, BIT_EXON_NEG, BIT_INTRON_NEG),
)


@dataclass(frozen=True, slots=True)
class SpliceJunction:
    """An eligible splice-junction boundary.

    Attributes
    ----------
    exon_side:
        ``"L"`` or ``"R"`` — which adjacent region is the exon side (the region to impute; the other
        side is the intron reference).
    strands:
        The matched strand set, a sorted tuple of ``TS_NEG``/``TS_POS`` — the strands that both splice
        at this boundary *and* carry the exon side's mature body. Equals ``exon_strands(exon side)``.
    """

    exon_side: str
    strands: tuple[int, ...]


def _junction_side(sig_l: int, sig_r: int, exon_bit: int, intron_bit: int) -> str | None:
    """Return the exon side (``"L"``/``"R"``) of a clean exon↔intron transition on one strand, else None.

    "Clean" means a pure exon→intron (or intron→exon) flip with **no mixed state** on either side: the
    exon side has the exon bit and *not* the intron bit, and the intron side the reverse. A side that is
    both exon and intron on this strand (retained intron) is not a clean junction and yields None.
    """
    l_exon, l_intron = bool(sig_l & exon_bit), bool(sig_l & intron_bit)
    r_exon, r_intron = bool(sig_r & exon_bit), bool(sig_r & intron_bit)
    if l_exon and not l_intron and r_intron and not r_exon:
        return "L"  # exon → intron
    if l_intron and not l_exon and r_exon and not r_intron:
        return "R"  # intron → exon
    return None


def splice_junction_eligibility(sig_l: int, sig_r: int) -> SpliceJunction | None:
    """Classify the boundary between adjacent regions ``sig_l`` and ``sig_r``.

    Returns a :class:`SpliceJunction` if the boundary is an eligible splice-junction reference for its
    exon-side region, else ``None`` (→ caller falls back to absolute density).

    Eligibility (the *matched-set* rule, unstranded — the conservative rule that holds with no strand
    information): the set of strands that splice at the boundary must equal the set of strands whose
    exon body the imputed region carries, the splicing strands must agree on a single exon side, and the
    set is non-empty. Rationale: a crossing-spliced read's intron gap sits at *this* boundary, so the
    spliced reference covers exactly the strands that splice here; for it to reference the *whole* mature
    exon body in the region, those strands must be precisely the region's exon strands.
    """
    sig_l = validate_signature(sig_l)
    sig_r = validate_signature(sig_r)

    matched: list[int] = []
    exon_sides: set[str] = set()
    for strand, exon_bit, intron_bit in _STRAND_LAYERS:
        side = _junction_side(sig_l, sig_r, exon_bit, intron_bit)
        if side is not None:
            matched.append(strand)
            exon_sides.add(side)

    if not matched or len(exon_sides) != 1:
        return None  # no junction, or the two strands splice toward different regions

    side = exon_sides.pop()
    region_sig = sig_l if side == "L" else sig_r
    exon_strands = set()
    if region_sig & BIT_EXON_POS:
        exon_strands.add(TS_POS)
    if region_sig & BIT_EXON_NEG:
        exon_strands.add(TS_NEG)

    if exon_strands != set(matched):
        return None  # AMBIG-exon vs single-strand junction, or antisense contamination
    return SpliceJunction(exon_side=side, strands=tuple(sorted(matched)))


def boundary_gdna_fraction(
    *,
    unspliced_gdna: float,
    unspliced_rna: float,
    spliced: float,
    eff_gdna: float,
    eff_rna: float,
) -> float:
    """gDNA molecular fraction at a splice-junction boundary, from FL-normalised crossing densities.

    Each crossing count is divided by its origin's FL-marginal crossing effective length
    (``boundary_eff_length`` of the gDNA / RNA pmf) to recover a density, and the gDNA share is

        ``f = ρ_gDNA / (ρ_gDNA + ρ_RNA)``,
        ``ρ_gDNA = unspliced_gdna / eff_gdna``,  ``ρ_RNA = (unspliced_rna + spliced) / eff_rna``.

    Two forms share this one arithmetic (design §1, §5):

    * **2-term** (no prior deconvolution) — the crossing-unspliced is the gDNA+nascent lump; pass it all
      as ``unspliced_gdna`` with ``unspliced_rna = 0``. ``f`` is then the *not-mature* fraction.
    * **3-term** (a gDNA/RNA split is already available, from strand pre-cleaning or a sweep carry-over)
      — pass the deconvolved ``unspliced_gdna`` / ``unspliced_rna``; nascent moves to the RNA side and
      ``f`` is the *pure gDNA* fraction.

    Returns ``nan`` when there is no crossing evidence (zero total density) — the caller treats that as
    "no information from this boundary" and uses the other anchor or the fallback.
    """
    rho_gdna = unspliced_gdna / eff_gdna
    rho_rna = (unspliced_rna + spliced) / eff_rna
    total = rho_gdna + rho_rna
    if total <= 0.0:
        return float("nan")
    return rho_gdna / total


def density_frac_to_count_frac(
    density_frac: float, eff_gdna_region: float, eff_rna_region: float
) -> float:
    """Convert a molecular-**density** gDNA fraction to a region's contained-**count** fraction.

    ``boundary_gdna_fraction`` returns ``f_b = ρ_gDNA/(ρ_gDNA+ρ_RNA)`` — a *density* fraction (the
    crossing counts were divided by their **crossing** eff-lengths). Applying it to a region's contained
    *count* ``M_region`` requires re-weighting by the region's own gDNA/RNA effective lengths, because a
    contained count is ``ρ·E_region`` and the two species have different ``E_region`` when their FL
    distributions differ::

        g = (f_b·E^g_region) / ( f_b·E^g_region + (1−f_b)·E^r_region )

    This is the exact count fraction under uniform density (``= d_g·E^g/(d_g·E^g + d_r·E^r)``); it is
    the **identity** when ``E^g_region = E^r_region`` (matched FL), and corrects the otherwise-large
    short-exon bias when they differ (the absolute-density count fraction already carries this factor —
    ``density·E^g_region/M_region`` — so this aligns the splice fraction with it). See
    ``docs/calibration/fl_consistency_diagnostic.md``. A degenerate region with no contained mass
    (both eff-lengths ``0``) returns ``density_frac`` unchanged (the value is moot — ``M_region≈0``).
    """
    num = density_frac * eff_gdna_region
    den = num + (1.0 - density_frac) * eff_rna_region
    return float(num / den) if den > 0.0 else float(density_frac)


def region_splice_gdna_frac(
    substrate,
    region_arrays,
    fallback_frac: np.ndarray,
    *,
    eff_gdna: float,
    eff_rna: float,
    eff_gdna_region: np.ndarray,
    eff_rna_region: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Upgrade the count gDNA fraction for exon regions with an eligible splice-junction boundary.

    For each region, its left boundary ``(r−1, r)`` and right boundary ``(r, r+1)`` are checked with
    :func:`splice_junction_eligibility`; a boundary anchors region ``r`` when it is eligible *and* names
    ``r`` as the exon side (``"R"`` for the left boundary, ``"L"`` for the right). At each anchoring
    boundary the **2-term** :func:`boundary_gdna_fraction` is computed from the crossing counts on
    ``r``'s side (``substrate.left[r]`` / ``substrate.right[r]`` — the per-side flux), using only sides
    that carry a mature (spliced) reference. Eligible regions take the mean of their anchor *density*
    fractions (the same anchor combination the absolute-density imputation uses), then convert that to
    the region's contained-**count** fraction via :func:`density_frac_to_count_frac` using the region's
    own gDNA/RNA effective lengths (``eff_gdna_region`` / ``eff_rna_region`` = ``region_eff_length`` of
    each FL pmf). The boundary fraction is a *density* fraction; without this conversion it carries a
    large short-exon bias whenever the gDNA and RNA FL distributions differ (see
    ``docs/calibration/fl_consistency_diagnostic.md``). Regions with no usable splice anchor keep
    ``fallback_frac`` (the absolute-density imputation — itself already count-consistent — always
    available). Returns ``(count_gdna_frac, n_upgraded)``.

    First cut — the **2-term** form (crossing-unspliced lumped as gDNA+nascent), correct for libraries
    without nascent RNA. The 3-term / strand-resolved sweep (design §5, §5.1) is a later layer; this
    function is the per-region direct anchor only (no run-interior carry-over yet).
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    eff_gdna_region = np.asarray(eff_gdna_region, dtype=np.float64)
    eff_rna_region = np.asarray(eff_rna_region, dtype=np.float64)
    r = sig.shape[0]
    frac = np.array(fallback_frac, dtype=np.float64, copy=True)

    left, right = substrate.left, substrate.right
    left_unspl = (left.n_unspliced_pos + left.n_unspliced_neg).astype(np.float64)
    left_spl = (left.n_spliced_sense + left.n_spliced_antisense).astype(np.float64)
    right_unspl = (right.n_unspliced_pos + right.n_unspliced_neg).astype(np.float64)
    right_spl = (right.n_spliced_sense + right.n_spliced_antisense).astype(np.float64)

    def _anchor(unspl: float, spl: float) -> float | None:
        # A splice anchor needs a mature (spliced) reference; without one there is nothing to debias
        # against, so defer to the fallback.
        if spl <= 0.0:
            return None
        f = boundary_gdna_fraction(
            unspliced_gdna=unspl, unspliced_rna=0.0, spliced=spl, eff_gdna=eff_gdna, eff_rna=eff_rna
        )
        return None if math.isnan(f) else f

    n_upgraded = 0
    for i in range(r):
        anchors: list[float] = []
        if i > 0 and ref_id[i] == ref_id[i - 1]:
            elig = splice_junction_eligibility(int(sig[i - 1]), int(sig[i]))
            if elig is not None and elig.exon_side == "R":
                a = _anchor(left_unspl[i], left_spl[i])
                if a is not None:
                    anchors.append(a)
        if i < r - 1 and ref_id[i] == ref_id[i + 1]:
            elig = splice_junction_eligibility(int(sig[i]), int(sig[i + 1]))
            if elig is not None and elig.exon_side == "L":
                a = _anchor(right_unspl[i], right_spl[i])
                if a is not None:
                    anchors.append(a)
        if anchors:
            density_frac = float(np.mean(anchors))  # mean anchor DENSITY fraction for this region
            frac[i] = density_frac_to_count_frac(
                density_frac, eff_gdna_region[i], eff_rna_region[i]
            )
            n_upgraded += 1
    return frac, n_upgraded
