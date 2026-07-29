"""Regression tests for ``capture_eff_length`` — the exon→region incidence must cover every exon.

The bug (fixed 2026-06-14): ``_exon_region_incidence`` read the lower region index from region
**starts** (``searchsorted(starts, a, side='left')``), which skips the region that *contains* an
exon's left edge whenever that edge falls in a region's interior. Region boundaries do NOT always
coincide with exon edges — ``build_region_partition`` merges adjacent **same-signature** segments,
so a shorter alternative-transcript exon starts interior to a merged region. The off-by-one dropped
fully-contained exons/spans entirely (factor=1, no contraction) and produced ``len_t < exonic``,
which is geometrically impossible (a transcript's regions must contain its exons). The fix reads the
lower bound from region **ends** (containment semantics).

⚠ **The shipped partition can no longer PRODUCE this geometry.** Every exon edge is a node
interface in the splice graph, so an index never yields a region with an interior exon edge. The
guarded property is nonetheless a property of the FUNCTION — ``_transcript_node_incidence`` must map
an exon to every region containing it, whatever the partition — and it is one that was found by
reasoning rather than by a test. So the coarse partition is built HERE, by hand
(:func:`_coarsened`), rather than taken from an index that would no longer supply one. Losing that
fixture would silently retire the regression.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.capture_eff_length import (
    _global_reference_density,
    _pooled_seam_arrays,
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
    # ⚠ D6 (2026-07-29): the STORED per-side density length is E[min(ℓ,L)]/2, and each face carries
    # ρ·gdna_boundary_len. This used to store the un-halved E[min(ℓ,L)] and deposit ρ·bl/2 — same face
    # mass, doubled length — which cancelled `_pooled_seam_arrays`' spurious ½ and hid the factor-2.
    bl = (
        np.minimum(rel, 180.0) / 2.0
    )  # per-side density length E[min(ℓ,L)]/2, clips on short regions
    ref_id = np.asarray(region_arrays.ref_id)
    n = rel.shape[0]
    contained = rho * rel
    side_right = np.zeros(n, dtype=np.float64)
    side_left = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]  # only genomically adjacent same-reference boundaries seam
        side_right[:-1] = np.where(same, rho * bl[:-1], 0.0)  # ρ·gdna_boundary_len per FACE
        side_left[1:] = np.where(same, rho * bl[1:], 0.0)
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


def _coarsened(idx) -> RegionArrays:
    """The index's nodes with adjacent equal-signature neighbours MERGED — a deliberately coarse
    partition whose regions contain interior exon edges.

    This is the geometry the off-by-one mishandled. It is constructed here because no index emits it
    any more, and the function must still be correct on it: a caller may legitimately hand
    ``_transcript_node_incidence`` any partition that contains the exons.
    """
    n = idx.nodes_df
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
    ra = _coarsened(idx)
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
    ra = _coarsened(idx)
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


def test_incidence_is_correct_on_the_v8_partition_too(misaligned_index):
    """The LIVE path (plan W1b): production builds this geometry with ``from_index``, on the v8 node
    partition, where the alternative exon start at 150 is a cut rather than a region interior.

    The three tests above pin the function against a partition that still has interior exon edges;
    this one pins it against the partition it actually runs on. Both properties must hold, and only
    one of them is now reachable from an index.
    """
    idx = misaligned_index
    ra = RegionArrays.from_index(idx)
    assert ra.n_regions == len(idx.nodes_df) > _coarsened(idx).n_regions  # it really does split
    assert 150 in set(ra.start.tolist())  # ...and it splits HERE

    inc_t, inc_r, *_ = _transcript_node_incidence(idx, ra)
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
    per-side density length ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})]) for the pooled seams) makes every node's
    density ρ, so the Laplace-smoothed IPR over any transcript's region set returns its full effective
    support (factor 1). With the genomic region_size_bp divisor a short exon would fabricate a
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
    node-density distribution the global detector can resolve), transcripts overlapping the DEPLETED regions
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


# --- nascent<mature inversion guard (2026-07-09): splice-junction seams ---------------------------
# A multi-exon mRNA and a single-exon nascent parent covering the SAME genomic span. A nascent's genomic
# node set STRICTLY CONTAINS its mature child's, so its EM effective length can never be shorter. Before
# the junction-seam fix, a multi-exon mRNA's span_full (with splice junctions DROPPED) fell below its
# contiguous FL-marginal length, and the fl/span_full ratio (growing with exon count) inflated the mature's
# eff_em ABOVE its nascent parent's under capture — the physically impossible inversion. These pin the fix.

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
    region_arrays: RegionArrays, density: np.ndarray, frag: float = 180.0
) -> CalibrationResult:
    """Deposition-faithful CalibrationResult for an arbitrary per-region gDNA DENSITY field, with an
    FL-marginal region support (``region_eff_len = size − frag``) and crossing support (``boundary_len =
    min(size, frag)``). This makes a multi-exon mRNA's junction-dropped ``span_full`` fall BELOW its
    contiguous FL-marginal length — the exact gap the junction seams close. Uniform density ⇒ every node's
    m/S = density ⇒ factor 1 (the bedrock invariant), independent of the field values."""
    rel = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    region_eff = np.maximum(rel - frag, 1e-9)
    bl = np.minimum(rel, frag) / 2.0  # D6: the STORED length is E[min(ℓ,L)]/2
    ref_id = np.asarray(region_arrays.ref_id)
    n = rel.shape[0]
    d = np.asarray(density, dtype=np.float64)
    contained = d * region_eff
    side_left = np.zeros(n, dtype=np.float64)
    side_right = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]
        side_right[:-1] = np.where(same, d[:-1] * bl[:-1], 0.0)  # ρ·gdna_boundary_len per FACE
        side_left[1:] = np.where(same, d[1:] * bl[1:], 0.0)
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
    ra = RegionArrays.from_index(idx)
    _, _, _, _, jt, jl, jr = _transcript_node_incidence(idx, ra)
    mrna, nasc = _tidx(idx, "mrna"), _tidx(idx, "nasc")
    assert (jt == mrna).sum() == 5, "a 6-exon mRNA must have 5 splice-junction seams"
    assert (jt == nasc).sum() == 0, "a single-exon nRNA must have NO splice junctions"
    starts, ends = np.asarray(ra.start), np.asarray(ra.end)
    for k in np.flatnonzero(jt == mrna):
        assert ends[jl[k]] <= starts[jr[k]], (
            "junction left flank must end at/before the right flank starts"
        )


def test_no_nascent_mature_inversion_under_capture(multiexon_index):
    """THE regression guard: under capture on a single exon, eff_em(nascent) >= eff_em(mature).

    Without the junction seams a 6-exon mRNA's fl/span_full ≈ 1.5 inflated its eff_em above its nascent
    parent's (an inversion, since the nascent's node set strictly contains the mature's). The imputed
    junction seams close the gap. Also asserts the mature genuinely CONTRACTS (the fix must not silently
    disable capture contraction)."""
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
    """Capture-off bit-identity WITH junction seams: uniform gDNA ⇒ eff_em == fl for the multi-exon mRNA."""
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    cal = _field_cal(ra, np.full(ra.n_regions, 0.02))
    fl = np.linspace(900.0, 2000.0, len(idx.t_df))
    eff = transcript_capture_eff_lengths(cal, ra, idx, fl)
    np.testing.assert_allclose(eff, fl, rtol=1e-9)


# --- the enriched-mode reference detector's core contract (locks the <5-node fallback + the enriched mode) ---


def test_global_reference_density_needs_five_gdna_nodes():
    # Fewer than 5 gDNA-bearing nodes ⇒ None (no reference ⇒ no contraction), even with a clean bimodal split.
    mass = np.array([100.0, 100.0, 1.0, 1.0])  # 4 gDNA-bearing nodes
    support = np.full(4, 100.0)
    assert _global_reference_density(mass, support) is None


def test_global_reference_density_single_enriched_region_is_none():
    # A single enriched region among empties ⇒ only 1 gDNA-bearing node ⇒ None (not a detectable pattern).
    mass = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    support = np.full(6, 100.0)
    assert _global_reference_density(mass, support) is None


def test_global_reference_density_bimodal_returns_enriched_mode_snapped():
    # >=5 gDNA nodes, bimodal (5 enriched at density 1.0 + 3 depleted at 0.01). The MASS-weighted KDE mode is
    # the enriched level, SNAPPED to an actual node density (exactly 1.0 here), not a grid value.
    mass = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 1.0, 1.0, 1.0])
    support = np.full(8, 100.0)
    rho = _global_reference_density(mass, support)
    assert rho == pytest.approx(1.0, rel=1e-9)


# ---------------------------------------------------------------------------
# D6 — the pooled-seam factor-2 (accumulator v5 §10.4; ledger W0.5).
# ---------------------------------------------------------------------------
#
# ⚠ These two tests are the FALSIFICATION tests for the fix: both MUST fail against the
# unpatched code, or the diagnosis is wrong. Verified failing at 3c293038 before the fix landed.
#
# The defect: ``gdna_boundary_len`` IS the halved per-side DENSITY length ``E[min(ℓ,L)]/2``
# (``effective_length.boundary_side_eff_length``, set at ``calibrate.py:226``), and the accumulator
# deposits ``ρ·gdna_boundary_len`` on EACH face — so a pooled seam holds ``ρ·(gbl_r + gbl_{r+1})``.
# Dividing that by the AVERAGE ``½·(gbl_r + gbl_{r+1})`` reads **2ρ**. The correct divisor is the SUM,
# which is what ``priors._gdna_region_node_arrays``'s own docstring says; the prose beside the code
# said "AVERAGE" and the code followed the prose.
#
# No test caught it because every fixture above builds ``gdna_boundary_len`` as the UN-halved
# ``E[min(ℓ,L)]`` and deposits ``ρ·bl/2`` per face — a doubled length that exactly cancels the
# spurious ½ — and because the ``min(ρ/ρ_ref, 1)`` clip rescues the uniform case regardless.


def _seam_faithful_cal(region_arrays: RegionArrays, density: np.ndarray, frag: float = 180.0):
    """A CalibrationResult built with PRODUCTION's ``gdna_boundary_len`` semantics.

    Unlike :func:`_field_cal` / :func:`_uniform_field_cal` above, ``gdna_boundary_len`` here is the
    genuinely halved per-side density length ``E[min(ℓ,L)]/2``, and each face carries
    ``ρ·gdna_boundary_len`` — exactly what ``boundary_side_eff_length``'s docstring specifies and what
    the accumulator deposits. A pooled seam therefore holds ``ρ·(gbl_r + gbl_{r+1})``.
    """
    rel = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    region_eff = np.maximum(rel - frag, 1e-9)
    gbl = np.minimum(rel, frag) / 2.0  # THE production semantic: E[min(ℓ,L)]/2
    ref_id = np.asarray(region_arrays.ref_id)
    n = rel.shape[0]
    d = np.asarray(density, dtype=np.float64)
    side_left = np.zeros(n, dtype=np.float64)
    side_right = np.zeros(n, dtype=np.float64)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]
        side_right[:-1] = np.where(same, d[:-1] * gbl[:-1], 0.0)  # ρ·gbl per FACE, not ρ·gbl/2
        side_left[1:] = np.where(same, d[1:] * gbl[1:], 0.0)
    z = np.zeros(n, dtype=np.float64)
    return CalibrationResult(
        mass_gdna_contained=d * region_eff,
        mass_rna_contained=z.copy(),
        mass_gdna_left=side_left,
        mass_rna_left=z.copy(),
        mass_gdna_right=side_right,
        mass_rna_right=z.copy(),
        mass_rna_spliced=z.copy(),
        gdna_boundary_len=gbl,
        gdna_region_eff_len=region_eff,
        gdna_density_global=float(np.mean(d)),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )


def test_pooled_seam_density_recovers_rho(multiexon_index):
    """⭐ D6: a pooled seam under a uniform field must read ρ, NOT 2ρ.

    Fails on the unpatched code with a ratio of exactly 2. Independently reproduced end-to-end by
    ``scratchpad/acc_seam_check.py``: 1.994 / 2.002 / 1.981 × truth at region lengths 2000 / 500 / 200.
    """
    ra = RegionArrays.from_index(multiexon_index)
    rho = 0.037
    cal = _seam_faithful_cal(ra, np.full(ra.n_regions, rho))
    seam_m, seam_S = _pooled_seam_arrays(cal, ra)
    live = seam_S > 0.0
    assert live.any(), "fixture produced no internal seams"
    np.testing.assert_allclose(seam_m[live] / seam_S[live], rho, rtol=1e-12)


def test_seam_inside_the_clip_band_still_contracts(multiexon_index):
    """⭐ D6, the reason no test caught it: the ``min(ρ/ρ_ref, 1)`` clip rescues the uniform case.

    Put the true seam density strictly inside ``(ρ_ref/2, ρ_ref)``. Read correctly it is BELOW the
    reference and must contract; read at 2ρ it lands ABOVE, clips to 1, and contributes no contraction
    at all — so the bug is silent on exactly the seams where the shrinkage was supposed to act.
    """
    idx = multiexon_index
    ra = RegionArrays.from_index(idx)
    n = ra.n_regions
    rho_ref, rho_seam = 1.0, 0.7  # 0.7 ∈ (0.5, 1.0): correct ⇒ contract; doubled ⇒ 1.4 ⇒ clipped
    dens = np.full(n, rho_ref)
    dens[1:] = rho_seam  # region 0 anchors ρ_ref; every seam among 1..n-1 sits in the band
    cal = _seam_faithful_cal(ra, dens)

    seam_m, seam_S = _pooled_seam_arrays(cal, ra)
    band = (seam_S > 0.0) & (seam_m > 0.0)
    assert band.any()
    ratio = (seam_m[band] / seam_S[band]) / rho_ref
    assert np.all(ratio < 1.0 - 1e-9), (
        f"seam density reads {ratio.max():.3f}×ρ_ref — at or above the reference it CLIPS and "
        "contributes no contraction, which is precisely how the factor-2 stayed invisible"
    )
