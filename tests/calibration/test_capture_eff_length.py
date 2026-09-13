"""``capture_eff_length``: the exon→region incidence must cover every region containing an exon.

``_transcript_region_incidence`` has to map an exon to every region that CONTAINS it, whatever
partition it is handed — reading the lower bound from region starts instead of region ends skips the
region holding an exon's left boundary whenever that boundary falls in a region's interior, which
drops fully contained exons and produces the geometrically impossible ``len_t < exonic``. The shipped
partition can no longer produce that geometry, because every exon boundary is a region interface in
the splice graph, so the coarse partition is built here by hand (:func:`_coarsened`) rather than taken
from an index: losing that fixture would silently retire the guard. The rest of the file holds what
the contraction itself must do — factor 1 under a uniform field, contraction but never expansion
under capture, no nascent/mature inversion once splice-junction boundaries are imputed, and a crossing
object reading its own density with no factor to get wrong.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.capture_eff_length import (
    _global_reference_density,
    _transcript_region_incidence,
    transcript_capture_eff_lengths,
)
from rigel.calibration.region_arrays import (
    RegionArrays,
    boundary_region_indices,
    region_right_boundary,
)
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from _index_builder import build_test_index

#: The gDNA crossing effective length every boundary carries. In production this is
#: ``crossing_eff_length(pmf, UNBOUNDED_REACH, UNBOUNDED_REACH) = mu_g − 1``, the SAME number at every
# boundary, because gDNA's template is the chromosome and takes no reach taper.
_CROSSING_EFF = 180.0


def _cal(region_arrays: RegionArrays, density, region_eff, boundary_eff) -> CalibrationResult:
    """THE fixture, and there is only one: a deposition-faithful result for an arbitrary per-region
    gDNA DENSITY field, where every object's mass is ``ρ × its own effective support`` on both axes.

    One builder rather than several is itself the guard. A contiguous boundary is a 0-bp boundary with
    ONE mass and ONE support, so there is no face, no half and nothing for two fixtures to disagree
    about — and a disagreement of exactly that shape, a per-side length halved in one place and not
    another, can cancel a spurious ½ and hide a factor of 2 from every assertion in this file.

    A boundary's density is its LEFT flank's. With a varying field the two flanks disagree, so the
    fixture must SAY which it means rather than average them into a number that is neither.
    """
    d = np.asarray(density, dtype=np.float64)
    region_eff = np.asarray(region_eff, dtype=np.float64)
    boundary_eff = np.asarray(boundary_eff, dtype=np.float64)
    lo, _hi = boundary_region_indices(np.asarray(region_arrays.ref_id))
    n = region_eff.shape[0]
    z = np.zeros(n, dtype=np.float64)
    ez = np.zeros(lo.shape[0], dtype=np.float64)
    return CalibrationResult(
        count_gdna_region=d * region_eff,
        count_rna_region=z.copy(),
        count_gdna_boundary=d[lo] * boundary_eff,
        count_rna_boundary=ez.copy(),
        count_rna_spliced_boundary=ez.copy(),
        # GEOMETRY, not a split: the mean conserved fragment-mass one crossing carries. 1.0 is the
        # identity — a boundary whose flanks both exceed every fragment length, where an incidence IS
        # a fragment — so a fixture that does not exercise K-inflation states it explicitly.
        boundary_mass_per_crossing=np.ones_like(ez),
        count_rna_sj=np.zeros(0, dtype=np.float64),
        boundary_spliced_mass_per_crossing=np.ones_like(ez),
        sj_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_region_eff_len=region_eff,
        gdna_boundary_eff_len=boundary_eff,
        rna_region_eff_len=region_eff,
        rna_boundary_eff_len=boundary_eff,
        gdna_frac_region=np.zeros(len(region_eff)),
        rna_pos_frac_region=np.zeros(len(region_eff)),
        rna_neg_frac_region=np.zeros(len(region_eff)),
        gdna_frac_boundary=np.zeros(len(boundary_eff)),
        rna_pos_frac_boundary=np.zeros(len(boundary_eff)),
        rna_neg_frac_boundary=np.zeros(len(boundary_eff)),
        gdna_density_global=float(d.mean()),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        n_boundaries=lo.shape[0],
        n_sj=0,
        config=CalibrationConfig(),
    )


def _uniform_field_cal(region_arrays: RegionArrays, rho: float) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field: every object's density is exactly ``rho``, with the region support
    the full genomic size. The factor-1-under-uniform invariant ⇒ every transcript's contraction factor
    is 1 ⇒ ``eff_em == fl``, even for exon flanks shorter than one fragment."""
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_boundaries = boundary_region_indices(np.asarray(region_arrays.ref_id))[0].shape[0]
    return _cal(
        region_arrays, np.full(size.shape[0], rho), size, np.full(n_boundaries, _CROSSING_EFF)
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


def _coarsened(idx) -> RegionArrays:
    """The index's regions with adjacent equal-signature neighbours MERGED — a deliberately coarse
    partition whose regions contain interior exon boundaries.

    This is the geometry the off-by-one mishandled. It is constructed here because no index emits it
    any more, and the function must still be correct on it: a caller may legitimately hand
    ``_transcript_region_incidence`` any partition that contains the exons.
    """
    n = idx.regions_df
    ref = n["ref_name"].astype(str).to_numpy()
    sig = n["signature"].to_numpy(np.uint8)
    start = n["start"].to_numpy(np.int64)
    end = n["end"].to_numpy(np.int64)
    keep = np.r_[True, (ref[1:] != ref[:-1]) | (sig[1:] != sig[:-1]) | (start[1:] != end[:-1])]
    i = np.flatnonzero(keep)
    last = np.r_[i[1:], ref.size] - 1
    merged = pd.DataFrame(
        {"ref_name": ref[i], "start": start[i], "end": end[last], "signature": sig[i]}
    )
    assert len(merged) < len(n), "the fixture must actually coarsen, or it tests nothing"
    return RegionArrays.from_frame(merged, idx.ref_name_to_id)


def test_incidence_maps_every_transcript(misaligned_index):
    """No transcript (mature or nRNA span) is dropped by the exon→region incidence."""
    idx = misaligned_index
    ra = _coarsened(idx)
    inc_t, *_ = _transcript_region_incidence(idx, ra)
    n_t = len(idx.t_df)
    mapped = set(int(t) for t in inc_t)
    dropped = set(range(n_t)) - mapped
    assert not dropped, f"transcripts dropped by incidence: {sorted(dropped)}"


def test_incidence_len_t_ge_exonic(misaligned_index):
    """Σ(region lengths over a transcript's incidence) ≥ its exonic/span length.

    Regions always *contain* the exons mapped to them, so ``len_t < exonic`` is impossible —
    the old searchsorted-on-starts off-by-one produced exactly that by skipping/dropping the
    region containing an interior exon boundary.
    """
    idx = misaligned_index
    ra = _coarsened(idx)
    inc_t, inc_r, *_ = _transcript_region_incidence(idx, ra)
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
    ra = _coarsened(idx)
    inc_t, inc_r, *_ = _transcript_region_incidence(idx, ra)
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


def test_incidence_is_correct_on_the_v8_partition_too(misaligned_index):
    """The LIVE path: production builds this geometry with ``from_index``, on a region partition where
    the alternative exon start at 150 is a region interface rather than a region interior.

    The three tests above pin the function against a partition that still has interior exon
    boundaries; this one pins it against the partition it actually runs on. Both properties must hold,
    and only one of them is reachable from an index.
    """
    idx = misaligned_index
    ra = RegionArrays.from_index(idx)
    assert ra.n_regions == len(idx.regions_df) > _coarsened(idx).n_regions  # it really does split
    assert 150 in set(ra.start.tolist())  # ...and it splits HERE

    inc_t, inc_r, *_ = _transcript_region_incidence(idx, ra)
    assert set(int(t) for t in inc_t) == set(range(len(idx.t_df))), "a transcript was dropped"

    len_t = np.zeros(len(idx.t_df))
    np.add.at(len_t, inc_t, np.asarray(ra.region_size_bp, dtype=np.float64)[inc_r])
    exonic = idx.t_df["length"].to_numpy(dtype=np.float64)
    bad = np.flatnonzero(len_t < exonic - 1e-6)
    assert bad.size == 0, f"len_t < exonic (geometrically impossible) for {bad.tolist()}"


# --- the bedrock invariant for the transcript path: factor = 1 under uniform gDNA -----------------


def test_transcript_factor_one_under_uniform_gdna(misaligned_index):
    """A uniform (unenriched) gDNA field contracts NO transcript's effective length: eff_em == fl.

    The density-correct effective-support divisor (gdna_region_eff_len for regions, the averaged
    per-side density length ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})]) for the pooled boundaries) makes every region's
    density ρ, so the Laplace-smoothed IPR over any transcript's region set returns its full effective
    support (factor 1). With the genomic ``region_size_bp`` divisor a short exon would fabricate a
    contraction even here — this pins that it does not."""
    idx = misaligned_index
    ra = RegionArrays.from_index(idx)
    cal = _uniform_field_cal(ra, rho=0.02)
    n_t = len(idx.t_df)
    fl = np.linspace(
        800.0, 2000.0, n_t
    )  # arbitrary FL-marginal lengths; the factor must be exactly 1
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    np.testing.assert_allclose(eff, fl, rtol=1e-9)


def test_transcript_contracts_under_concentrated_gdna(multiexon_index):
    """Under a realistic capture field (a subset of regions enriched, the rest depleted — a genuine bimodal
    region-density distribution the global detector can resolve), transcripts overlapping the DEPLETED regions
    contract below their FL-marginal length, and contraction never expands. (A single enriched region is NOT
    a detectable capture pattern under the global reference — that regime is correctly left uncontracted.)"""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    n = ra.n_regions
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    dens = np.full(n, 0.1)  # depleted (off-target) background
    dens[(ends > 1000) & (starts < 1500)] = 100.0  # capture the first exon (enriched on-target)
    cal = _field_cal(ra, dens)
    fl = np.maximum(idx.t_df["length"].to_numpy(dtype=np.float64) - 180.0, 1.0)
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    assert np.all(eff <= fl + 1e-9)  # contraction never expands
    assert np.any(eff < fl - 1e-6)  # at least one transcript genuinely contracts


# --- nascent<mature inversion guard: splice-junction boundaries ---------------------------------------
# A multi-exon mRNA and a single-exon nascent parent covering the SAME genomic span. A nascent's genomic
# region set STRICTLY CONTAINS its mature child's, so its EM effective length can never be shorter. With
# the splice junctions DROPPED a multi-exon mRNA's span_full falls below its contiguous FL-marginal
# length, and the fl/span_full ratio (growing with exon count) inflates the mature's eff_em ABOVE its
# nascent parent's under capture — a physically impossible inversion. Imputing the sj boundaries is what
# closes the gap, and these pin it.

# six 500bp exons + a single-exon nascent covering the whole 1000..6500 span (genomic order per transcript
# so the incidence pairs adjacent exons into splice junctions).
_MULTIEXON_GTF = (
    "".join(
        f'chr1\ttest\texon\t{s + 1}\t{s + 500}\t.\t+\t.\tgene_id "gm"; transcript_id "mrna";\n'
        for s in range(1000, 6500, 1000)
    )
    + 'chr1\ttest\texon\t1001\t6500\t.\t+\t.\tgene_id "gm"; transcript_id "nasc";\n'
)


@pytest.fixture(scope="module")
def multiexon_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _MULTIEXON_GTF, genome_size=7000, name="multiexon")


def _tidx(idx, tid: str) -> int:
    tdf = idx.t_df
    return int(tdf.loc[tdf["t_id"] == tid, "t_index"].iloc[0])


def _field_cal(
    region_arrays: RegionArrays, density: np.ndarray, frag: float = _CROSSING_EFF
) -> CalibrationResult:
    """An arbitrary per-region gDNA DENSITY field with an FL-MARGINAL region support
    (``region_eff = size − frag``). That makes a multi-exon mRNA's sj-dropped ``span_full`` fall
    BELOW its contiguous FL-marginal length — the exact gap the sj boundaries close. Uniform density ⇒
    every object's m/S = density ⇒ factor 1 (the bedrock invariant), independent of the field values."""
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_boundaries = boundary_region_indices(np.asarray(region_arrays.ref_id))[0].shape[0]
    return _cal(region_arrays, density, np.maximum(size - frag, 1e-9), np.full(n_boundaries, frag))


def test_sj_incidence_multiexon_only(multiexon_index):
    """A multi-exon mRNA yields one splice-junction boundary per adjacent exon pair; a single-exon nRNA yields
    none, and each sj's flanking regions straddle the intron between the two exons."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    _, _, _, _, jt, jl, jr = _transcript_region_incidence(idx, ra)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert (jt == mrna).sum() == 5, "a 6-exon mRNA must have 5 splice-junction boundaries"
    assert (jt == nasc).sum() == 0, "a single-exon nRNA must have NO splice junctions"
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    for k in np.flatnonzero(jt == mrna):
        assert ends[jl[k]] <= starts[jr[k]], (
            "sj left flank must end at/before the right flank starts"
        )


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """THE regression guard: under capture on a single exon, eff_em(nascent) >= eff_em(mature).

    Without the sj boundaries a 6-exon mRNA's fl/span_full ratio inflates its eff_em above its nascent
    parent's, which is an inversion because the nascent's region set strictly contains the mature's.
    The imputed sj boundaries close the gap. Also asserts the mature genuinely CONTRACTS, since a
    guard that silently disabled capture contraction would pass the inequality trivially."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    n = ra.n_regions
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    dens = np.full(n, 0.1)  # depleted off-target
    dens[(ends > 1000) & (starts < 1500)] = 100.0  # capture the first exon [1000,1500)
    cal = _field_cal(ra, dens)
    frag = 180.0
    fl = np.maximum(
        idx.t_df["length"].to_numpy(dtype=np.float64) - frag, 1.0
    )  # contiguous FL-marginal
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert eff[nasc] >= eff[mrna] - 1e-6, (
        f"INVERSION: nascent eff_em {eff[nasc]:.2f} < mature eff_em {eff[mrna]:.2f}"
    )
    assert eff[mrna] < fl[mrna] - 1e-6, "the mature must genuinely contract under capture"


def test_spliced_factor_one_under_uniform(multiexon_index):
    """Capture-off bit-identity WITH sj boundaries: uniform gDNA ⇒ eff_em == fl for the multi-exon mRNA."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    cal = _field_cal(ra, np.full(ra.n_regions, 0.02))
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    np.testing.assert_allclose(eff, fl, rtol=1e-9)


# --- the enriched-mode reference detector's core contract (locks the <5-region fallback + the enriched mode) ---


def test_global_reference_density_needs_five_gdna_regions():
    # Fewer than 5 gDNA-bearing regions ⇒ None (no reference ⇒ no contraction), even with a clean bimodal split.
    mass = np.array([100.0, 100.0, 1.0, 1.0])  # 4 gDNA-bearing regions
    support = np.full(4, 100.0)
    assert _global_reference_density(mass, support) is None


def test_global_reference_density_single_enriched_region_is_none():
    # A single enriched region among empties ⇒ only 1 gDNA-bearing region ⇒ None (not a detectable pattern).
    mass = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    support = np.full(6, 100.0)
    assert _global_reference_density(mass, support) is None


def test_global_reference_density_bimodal_returns_enriched_mode_snapped():
    # >=5 gDNA regions, bimodal (5 enriched at density 1.0 + 3 depleted at 0.01). The MASS-weighted KDE mode is
    # the enriched level, SNAPPED to an actual region density (exactly 1.0 here), not a grid value.
    mass = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 1.0, 1.0, 1.0])
    support = np.full(8, 100.0)
    rho = _global_reference_density(mass, support)
    assert rho == pytest.approx(1.0, rel=1e-9)


# ---------------------------------------------------------------------------
# The crossing-object density — the TRAPS: prefer-shares-to-differences factor-2 guard.
# ---------------------------------------------------------------------------
#
# The arithmetic that once made a factor of 2 available here is unrepresentable now: a contiguous
# boundary is a 0-bp boundary with one mass and one support, so there is no per-side length to halve,
# no pair of faces to sum, and no choice between a sum and an average. The PROPERTY that guarded is
# still real and is kept below — a crossing object under a uniform field must read ρ, and a boundary
# genuinely below the reference density must contract rather than clip.


def test_a_crossing_object_under_a_uniform_field_reads_RHO(multiexon_index):
    """One half of TRAPS: prefer-shares-to-differences. A boundary's mass over its own support is the
    true density — exactly, with no factor to get wrong — so a factor-of-2 regression in the crossing
    arithmetic reappears here first."""
    ra = RegionArrays.from_index(multiexon_index)
    rho = 0.037
    cal = _field_cal(ra, np.full(ra.n_regions, rho))
    # Read straight off the BOUNDARY axis: the incidence helper emits a boundary index, so a
    # boundary's density is read where it lives rather than through a region-shaped copy.
    boundary_mass = np.asarray(cal.count_gdna_boundary, dtype=np.float64)
    boundary_support = np.asarray(cal.gdna_boundary_eff_len, dtype=np.float64)
    live = boundary_support > 0.0
    assert live.any(), "the fixture produced no boundaries"
    np.testing.assert_allclose(boundary_mass[live] / boundary_support[live], rho, rtol=1e-12)


def test_a_boundary_below_the_reference_density_CONTRACTS_rather_than_clipping(multiexon_index):
    """TRAPS: prefer-shares-to-differences' other half, kept because the ``min(ρ/ρ_ref, 1)`` clip is
    still there and still hides anything that reads a density too HIGH.

    Put the true boundary density strictly inside ``(ρ_ref/2, ρ_ref)``: read correctly it is below the
    reference and must contract; read at any inflated multiple it lands above, clips to 1, and
    contributes no contraction at all — silent on exactly the boundaries the shrinkage exists to act on.

    The ρ_ref anchor sits on the LAST region, not the first. A boundary takes its LEFT flank's density
    (:func:`_cal`), so anchoring on region 0 would put boundary 0 itself at ρ_ref and the assertion would
    fail on the fixture rather than on the code.
    """
    ra = RegionArrays.from_index(multiexon_index)
    rho_ref, rho_boundary = 1.0, 0.7  # 0.7 ∈ (0.5, 1.0)
    dens = np.full(ra.n_regions, rho_boundary)
    dens[-1] = rho_ref  # the last region anchors ρ_ref and is no boundary's LEFT flank
    cal = _field_cal(ra, dens)

    boundary_mass = np.asarray(cal.count_gdna_boundary, dtype=np.float64)
    boundary_support = np.asarray(cal.gdna_boundary_eff_len, dtype=np.float64)
    band = (boundary_support > 0.0) & (boundary_mass > 0.0)
    assert band.any()
    ratio = (boundary_mass[band] / boundary_support[band]) / rho_ref
    assert np.all(ratio < 1.0 - 1e-9), (
        f"boundary density reads {ratio.max():.3f}×ρ_ref — at or above the reference it CLIPS and "
        "contributes no contraction, which is precisely how the factor-2 stayed invisible"
    )


# --- the boundary axis is a BOUNDARY index, and only a MULTI-reference index can prove it -------------

_TWO_REF_GTF = "".join(
    f'chrA\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "ga"; transcript_id "ta";\n'
    for s in (200, 800)
) + "".join(
    f'chrB\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "gb"; transcript_id "tb";\n'
    for s in (200, 800)
)


@pytest.fixture(scope="module")
def two_ref_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _TWO_REF_GTF, name="tworef", refs={"chrA": 2000, "chrB": 2000}
    )


def test_the_boundary_incidence_is_an_BOUNDARY_index_not_a_left_region_index(two_ref_index):
    """The gate a single-reference fixture cannot provide.

    ``region_right_boundary`` numbers boundaries over adjacent same-reference region pairs, so on ONE
    reference ``boundary(r) == r`` and a left-region index is indistinguishable from a boundary index.
    PERTURBATION: substituting one for the other changes nothing on any single-reference fixture here,
    which is how the conversion came to be untested (``TRAPS: perturb-every-gate``).

    On a second reference the two axes diverge by one per preceding reference boundary, and indexing a
    per-boundary array with a region index then reads THE WRONG BOUNDARY'S MASS — silently, since both
    are in range. This pins the axis: every emitted boundary index must be a valid boundary whose
    flanking regions are the ones the transcript actually crosses.
    """
    ra = RegionArrays.from_index(two_ref_index)
    _rt, _rr, bt, br, *_ = _transcript_region_incidence(two_ref_index, ra)
    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    assert br.size, "the fixture produced no interior boundaries"
    assert br.max() < lo.shape[0], "a boundary index outside the boundary axis"

    # the discriminating claim: each emitted boundary's flanks are same-reference neighbours, and the
    # transcript that emitted it overlaps BOTH of them.
    ref_id = np.asarray(ra.ref_id)
    np.testing.assert_array_equal(ref_id[lo[br]], ref_id[hi[br]])
    reg_t, reg_r = _rt, _rr
    for t, e in zip(bt.tolist(), br.tolist()):
        owned = set(reg_r[reg_t == t].tolist())
        assert lo[e] in owned and hi[e] in owned, (
            f"transcript {t} was given boundary {e} between regions {lo[e]},{hi[e]} — which it does not span"
        )

    # ...and the two axes genuinely differ here, so the assertions above are not vacuous
    right_boundary = region_right_boundary(ref_id)
    assert not np.array_equal(right_boundary[: lo.shape[0]], np.arange(lo.shape[0])), (
        "fixture is degenerate: the region and boundary axes coincide, so this test proves nothing"
    )
