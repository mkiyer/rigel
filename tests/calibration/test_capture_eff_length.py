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
    _transcript_node_incidence,
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
    inc_t, *_ = _transcript_node_incidence(idx, ra)
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
    inc_t, inc_r, *_ = _transcript_node_incidence(idx, ra)
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
    inc_t, inc_r, *_ = _transcript_node_incidence(idx, ra)
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
    inc_t, inc_r, *_ = _transcript_node_incidence(idx, ra)
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


# --- nascent<mature inversion guard (2026-07-09): splice-junction seams ---------------------------
# A multi-exon mRNA and a single-exon nascent parent covering the SAME genomic span. A nascent's genomic
# node set STRICTLY CONTAINS its mature child's, so its EM effective length can never be shorter. Before
# the junction-seam fix, a multi-exon mRNA's span_full (with splice junctions DROPPED) fell below its
# contiguous FL-marginal length, and the fl/span_full ratio (growing with exon count) inflated the mature's
# eff_em ABOVE its nascent parent's under capture — the physically impossible inversion. These pin the fix.

# six 500bp exons + a single-exon nascent covering the whole 1000..6500 span (genomic order per transcript
# so the incidence pairs adjacent exons into splice junctions).
_MULTIEXON_GTF = "".join(
    f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t+\t.\tgene_id "gm"; transcript_id "mrna";\n'
    for s in range(1000, 6500, 1000)
) + 'chr1\ttest\texon\t1001\t6500\t.\t+\t.\tgene_id "gm"; transcript_id "nasc";\n'


@pytest.fixture(scope="module")
def multiexon_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _MULTIEXON_GTF, genome_size=7000, name="multiexon")


def _tidx(idx, tid: str) -> int:
    tdf = idx.t_df
    return int(tdf.loc[tdf["t_id"] == tid, "t_index"].iloc[0])


def _field_cal(region_arrays: RegionArrays, density: np.ndarray, frag: float = 180.0) -> CalibrationResult:
    """Deposition-faithful CalibrationResult for an arbitrary per-region gDNA DENSITY field, with an
    FL-marginal region support (``region_eff_len = size − frag``) and crossing support (``boundary_len =
    min(size, frag)``). This makes a multi-exon mRNA's junction-dropped ``span_full`` fall BELOW its
    contiguous FL-marginal length — the exact gap the junction seams close. Uniform density ⇒ every node's
    m/S = density ⇒ factor 1 (the bedrock invariant), independent of the field values."""
    rel = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    region_eff = np.maximum(rel - frag, 1e-9)
    bl = np.minimum(rel, frag)
    ref_id = np.asarray(region_arrays.ref_id)
    n = rel.shape[0]
    d = np.asarray(density, dtype=np.float64)
    contained = d * region_eff
    side_left = np.zeros(n, dtype=np.float64)
    side_right = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]
        side_right[:-1] = np.where(same, d[:-1] * bl[:-1] / 2.0, 0.0)
        side_left[1:] = np.where(same, d[1:] * bl[1:] / 2.0, 0.0)
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
        gdna_region_eff_len=region_eff,
        gdna_density_global=float(d.mean()),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )


def test_junction_incidence_multiexon_only(multiexon_index):
    """A multi-exon mRNA yields one splice-junction seam per adjacent exon pair; a single-exon nRNA yields
    none, and each junction's flanking regions straddle the intron between the two exons."""
    idx = multiexon_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    _, _, _, _, jt, jl, jr = _transcript_node_incidence(idx, ra)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert (jt == mrna).sum() == 5, "a 6-exon mRNA must have 5 splice-junction seams"
    assert (jt == nasc).sum() == 0, "a single-exon nRNA must have NO splice junctions"
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    for k in np.flatnonzero(jt == mrna):
        assert ends[jl[k]] <= starts[jr[k]], "junction left flank must end at/before the right flank starts"


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """THE regression guard: under capture on a single exon, eff_em(nascent) >= eff_em(mature).

    Without the junction seams a 6-exon mRNA's fl/span_full ≈ 1.5 inflated its eff_em above its nascent
    parent's (an inversion, since the nascent's node set strictly contains the mature's). The imputed
    junction seams close the gap. Also asserts the mature genuinely CONTRACTS (the fix must not silently
    disable capture contraction)."""
    idx = multiexon_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    n = ra.n_regions
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    dens = np.full(n, 0.1)                       # depleted off-target
    dens[(ends > 1000) & (starts < 1500)] = 100.0  # capture the first exon [1000,1500)
    cal = _field_cal(ra, dens)
    frag = 180.0
    fl = np.maximum(idx.t_df["length"].to_numpy(dtype=np.float64) - frag, 1.0)  # contiguous FL-marginal
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert eff[nasc] >= eff[mrna] - 1e-6, (
        f"INVERSION: nascent eff_em {eff[nasc]:.2f} < mature eff_em {eff[mrna]:.2f}"
    )
    assert eff[mrna] < fl[mrna] - 1e-6, "the mature must genuinely contract under capture"


def test_spliced_factor_one_under_uniform(multiexon_index):
    """Capture-off bit-identity WITH junction seams: uniform gDNA ⇒ eff_em == fl for the multi-exon mRNA."""
    idx = multiexon_index
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cal = _field_cal(ra, np.full(ra.n_regions, 0.02))
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    np.testing.assert_allclose(eff, fl, rtol=1e-9)
