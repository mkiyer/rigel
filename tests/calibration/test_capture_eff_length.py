"""Regression tests for ``capture_eff_length`` — the exon→region incidence must cover every exon.

The bug (fixed 2026-06-14): ``_exon_region_incidence`` read the lower region index from region
**starts** (``searchsorted(starts, a, side='left')``), which skips the region that *contains* an
exon's left edge whenever that edge falls in a region's interior. Region boundaries do NOT always
coincide with exon edges — ``build_region_partition`` merges adjacent **same-signature** segments,
so a shorter alternative-transcript exon starts interior to a merged region. The off-by-one dropped
fully-contained exons/spans entirely (factor=1, no contraction) and produced ``len_t < exonic``,
which is geometrically impossible (a transcript's regions must contain its exons). The fix reads the
lower bound from region **ends** (containment semantics).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.capture_eff_length import (
    _exon_region_incidence,
    transcript_capture_eff_lengths,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from conftest import build_test_index


def _uniform_field_cal(region_arrays: RegionArrays, rho: float) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field over the region partition, built per the TRUE accumulator
    deposition law: every node (region contained + pooled seam) has density exactly ``rho``.
    ``contained[r] = rho·region_eff_len[r]``; the per-side density length ``boundary_len[r] =
    min(region_size_bp[r], 180)`` VARIES (short regions clip below E[ℓ]); each region's two sides carry
    ``rho·boundary_len[r]/2`` so the POOLED seam(r,r+1) = ``rho·(boundary_len[r]+boundary_len[r+1])/2`` =
    ``rho·S_s``. The factor-1-under-uniform invariant ⇒ every transcript's contraction factor is 1 ⇒
    eff_em == fl, EVEN for short exon flanks (the deposition-faithful averaged seam support)."""
    rel = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    bl = np.minimum(rel, 180.0)  # per-side density length E[min(ℓ,L)] — varies, clips on short regions
    ref_id = np.asarray(region_arrays.ref_id)
    n = rel.shape[0]
    contained = rho * rel
    side_right = np.zeros(n, dtype=np.float64)
    side_left = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # only genomically adjacent same-reference boundaries seam
        side_right[:-1] = np.where(same, rho * bl[:-1] / 2.0, 0.0)
        side_left[1:] = np.where(same, rho * bl[1:] / 2.0, 0.0)
    z = np.zeros(n, dtype=np.float64)
    return CalibrationResult(
        mass_gdna_contained=contained,
        mass_rna_contained=z.copy(),
        mass_gdna_left=side_left,
        mass_rna_left=z.copy(),
        mass_gdna_right=side_right,
        mass_rna_right=z.copy(),
        mass_rna_spliced=z.copy(),
        gdna_boundary_len=bl,
        gdna_region_eff_len=rel.copy(),
        gdna_density_global=rho,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )

# Two same-strand transcripts sharing gene g1. t1's first exon [150, 300) is a sub-interval of t0's
# [100, 300); both halves carry the identical EXON_POS signature, so the partition MERGES them into a
# single region [100, 300) whose interior holds t1's exon start (150) — the exact misalignment the bug
# mishandled (t1's first exon would otherwise be dropped, since 150 is interior to the region).
_MISALIGNED_GTF = """\
chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t151\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
"""


@pytest.fixture(scope="module")
def misaligned_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _MISALIGNED_GTF, genome_size=1000, name="misaligned")


def test_incidence_maps_every_transcript(misaligned_index):
    """No transcript (mature or nRNA span) is dropped by the exon→region incidence."""
    idx = misaligned_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    inc_t, _ = _exon_region_incidence(idx, ra)
    n_t = len(idx.t_df)
    mapped = set(int(t) for t in inc_t)
    dropped = set(range(n_t)) - mapped
    assert not dropped, f"transcripts dropped by incidence: {sorted(dropped)}"


def test_incidence_len_t_ge_exonic(misaligned_index):
    """Σ(region lengths over a transcript's incidence) ≥ its exonic/span length.

    Regions always *contain* the exons mapped to them, so ``len_t < exonic`` is impossible —
    the old searchsorted-on-starts off-by-one produced exactly that by skipping/dropping the
    region containing an interior exon edge.
    """
    idx = misaligned_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    inc_t, inc_r = _exon_region_incidence(idx, ra)
    n_t = len(idx.t_df)
    region_len = np.asarray(ra.region_size_bp, dtype=np.float64)
    len_t = np.zeros(n_t)
    np.add.at(len_t, inc_t, region_len[inc_r])
    exonic = idx.t_df["length"].to_numpy(dtype=np.float64)
    bad = np.flatnonzero(len_t < exonic - 1e-6)
    assert bad.size == 0, (
        f"len_t < exonic (geometrically impossible) for transcripts {bad.tolist()}: "
        f"len_t={len_t[bad].tolist()} exonic={exonic[bad].tolist()}"
    )


def test_merged_region_interior_exon_is_mapped(misaligned_index):
    """The transcript whose exon starts interior to a merged region maps to that region.

    Directly pins the fix: t1's first exon [150,300) lies inside the merged region [100,300); its
    incidence must include that region (the old code skipped to the next region, or dropped it)."""
    idx = misaligned_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    inc_t, inc_r = _exon_region_incidence(idx, ra)
    starts = np.asarray(ra.start)
    ends = np.asarray(ra.end)
    # the region covering position 150 (interior to a merged exon region)
    covering = np.flatnonzero((starts <= 150) & (ends > 150))
    assert covering.size == 1, f"expected one region covering pos 150, got {covering.tolist()}"
    region_at_150 = int(covering[0])
    # every transcript whose exonic footprint includes position 150 must map to that region
    t_df = idx.t_df
    for t in range(len(t_df)):
        a, b = int(t_df["start"].iloc[t]), int(t_df["end"].iloc[t])
        if a <= 150 < b:
            regions_for_t = set(int(inc_r[k]) for k in np.flatnonzero(inc_t == t))
            assert region_at_150 in regions_for_t, (
                f"transcript {t} spans pos 150 but its incidence {sorted(regions_for_t)} omits "
                f"the covering region {region_at_150}"
            )


# --- the bedrock invariant for the transcript path: factor = 1 under uniform gDNA -----------------


def test_transcript_factor_one_under_uniform_gdna(misaligned_index):
    """A uniform (unenriched) gDNA field contracts NO transcript's effective length: eff_em == fl.

    The density-correct effective-support divisor (gdna_region_eff_len for regions, the averaged
    per-side density length ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})]) for the pooled seams) makes every node's
    density ρ, so the Laplace-smoothed IPR over any transcript's region set returns its full effective
    support (factor 1). With the genomic region_size_bp divisor a short exon would fabricate a
    contraction even here — this pins that it does not."""
    idx = misaligned_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cal = _uniform_field_cal(ra, rho=0.02)
    n_t = len(idx.t_df)
    fl = np.linspace(800.0, 2000.0, n_t)  # arbitrary FL-marginal lengths; the factor must be exactly 1
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    np.testing.assert_allclose(eff, fl, rtol=1e-9)


def test_transcript_contracts_under_concentrated_gdna(misaligned_index):
    """Sanity: with gDNA concentrated on a SINGLE region (capture-like, non-uniform), the transcripts
    overlapping it contract below their FL-marginal length (factor < 1), while the contraction never
    exceeds fl. Confirms the method still does its job once the field is non-uniform."""
    idx = misaligned_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    n = ra.n_regions
    inc_t, inc_r = _exon_region_incidence(idx, ra)
    assert inc_t.size, "the index must map at least one transcript to a region"
    enriched_r = int(inc_r[0])  # a region genuinely inside some transcript's footprint
    contained = np.zeros(n, dtype=np.float64)
    contained[enriched_r] = 500.0  # all gDNA piled on one region — a sharply non-uniform (enriched) field
    z = np.zeros(n, dtype=np.float64)
    cal = CalibrationResult(
        mass_gdna_contained=contained,
        mass_rna_contained=z.copy(),
        mass_gdna_left=z.copy(),
        mass_rna_left=z.copy(),
        mass_gdna_right=z.copy(),
        mass_rna_right=z.copy(),
        mass_rna_spliced=z.copy(),
        gdna_boundary_len=np.full(n, 180.0, dtype=np.float64),
        gdna_region_eff_len=np.asarray(ra.region_size_bp, dtype=np.float64),
        gdna_density_global=0.01,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )
    n_t = len(idx.t_df)
    fl = np.full(n_t, 1500.0)
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    assert np.all(eff <= fl + 1e-9)  # contraction never expands
    overlaps = sorted({int(t) for t, r in zip(inc_t, inc_r, strict=True) if int(r) == enriched_r})
    assert overlaps, "expected at least one transcript overlapping the enriched region"
    assert np.any(eff[overlaps] < fl[overlaps] - 1e-6)  # at least one genuinely contracts
