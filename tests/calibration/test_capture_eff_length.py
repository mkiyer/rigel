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
    _transcript_node_incidence,
    transcript_capture_eff_lengths,
)
from rigel.calibration.region_arrays import RegionArrays, edge_node_indices, node_right_edge
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from conftest import build_test_index

#: The gDNA crossing effective length every line carries. In production this is
#: ``crossing_eff_length(pmf, UNBOUNDED_REACH, UNBOUNDED_REACH) = mu_g − 1``, the SAME number at every
# line, because gDNA's template is the chromosome and takes no reach taper.
_CROSSING_EFF = 180.0


def _cal(region_arrays: RegionArrays, density, node_eff, edge_eff) -> CalibrationResult:
    """⭐ **THE fixture — there is only one now.** A deposition-faithful result for an arbitrary
     per-node gDNA DENSITY field: every object's mass is ``ρ × its own effective support``, on both axes.

     ⛔ **Three near-identical builders used to live here** — ``_uniform_field_cal``, ``_field_cal`` and
     ``_seam_faithful_cal`` — and they differed in exactly one thing: whether ``gdna_boundary_len`` was
     the halved per-side density length ``E[min(ℓ,L)]/2`` or the un-halved ``E[min(ℓ,L)]``. Two of them
     stored the un-halved length AND deposited half the mass, which cancelled a spurious ½ in
     ``_pooled_seam_arrays`` and hid an exact factor of 2 from every assertion in this file for months
    A contiguous edge is a 0-bp line with ONE mass and ONE support, so
     there is no face, no half, and nothing for three fixtures to disagree about.

     ⚠ A line's density is its **left flank's**. With a varying field the two flanks disagree, so the
     fixture must SAY which it means rather than average them into a number that is neither.
    """
    d = np.asarray(density, dtype=np.float64)
    node_eff = np.asarray(node_eff, dtype=np.float64)
    edge_eff = np.asarray(edge_eff, dtype=np.float64)
    lo, _hi = edge_node_indices(np.asarray(region_arrays.ref_id))
    n = node_eff.shape[0]
    z = np.zeros(n, dtype=np.float64)
    ez = np.zeros(lo.shape[0], dtype=np.float64)
    return CalibrationResult(
        mass_gdna_node=d * node_eff,
        mass_rna_node=z.copy(),
        mass_gdna_edge=d[lo] * edge_eff,
        mass_rna_edge=ez.copy(),
        mass_rna_spliced_edge=ez.copy(),
        # ⭐ GEOMETRY, not a split: the mean conserved fragment-mass one crossing carries. 1.0 is the
        # identity — a line whose flanks both exceed every fragment length, where an incidence IS
        # a fragment — so a fixture that does not exercise K-inflation states it explicitly.
        edge_mass_per_crossing=np.ones_like(ez),
        mass_rna_junction=np.zeros(0, dtype=np.float64),
        edge_spliced_mass_per_crossing=np.ones_like(ez),
        junction_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_node_eff_len=node_eff,
        gdna_edge_eff_len=edge_eff,
        rna_node_eff_len=node_eff,
        rna_edge_eff_len=edge_eff,
        gdna_density_global=float(d.mean()),
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_nodes=n,
        n_edges=lo.shape[0],
        n_junctions=0,
        config=CalibrationConfig(),
    )


def _uniform_field_cal(region_arrays: RegionArrays, rho: float) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field: every object's density is exactly ``rho``, with the node support
    the full genomic size. The factor-1-under-uniform invariant ⇒ every transcript's contraction factor
    is 1 ⇒ ``eff_em == fl``, even for exon flanks shorter than one fragment."""
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_edges = edge_node_indices(np.asarray(region_arrays.ref_id))[0].shape[0]
    return _cal(region_arrays, np.full(size.shape[0], rho), size, np.full(n_edges, _CROSSING_EFF))


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
    region_arrays: RegionArrays, density: np.ndarray, frag: float = _CROSSING_EFF
) -> CalibrationResult:
    """An arbitrary per-node gDNA DENSITY field with an **FL-marginal** node support
    (``node_eff = size − frag``). That makes a multi-exon mRNA's junction-dropped ``span_full`` fall
    BELOW its contiguous FL-marginal length — the exact gap the junction seams close. Uniform density ⇒
    every object's m/S = density ⇒ factor 1 (the bedrock invariant), independent of the field values."""
    size = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_edges = edge_node_indices(np.asarray(region_arrays.ref_id))[0].shape[0]
    return _cal(region_arrays, density, np.maximum(size - frag, 1e-9), np.full(n_edges, frag))


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
# The crossing-object density — what survives of the TRAPS: prefer-shares-to-differences factor-2 guard.
# ---------------------------------------------------------------------------
#
# ⛔ **TRAPS: prefer-shares-to-differences' two falsification tests are DELETED, and the reason is that the defect is now
# unrepresentable.** TRAPS: prefer-shares-to-differences was: ``gdna_boundary_len`` IS the halved per-side density length
# ``E[min(ℓ,L)]/2``, the accumulator deposits ``ρ·gdna_boundary_len`` on EACH face, so a pooled seam
# holds ``ρ·(gbl_r + gbl_{r+1})`` — and dividing that by the AVERAGE read **2ρ**. Every quantity in
# that sentence is gone with the faces. There is no per-side length to halve, no pair of faces to sum,
# and no choice between a sum and an average: a contiguous edge is a 0-bp line with one mass and one
# support ().
#
# ⚠ **The PROPERTY they protected is still real**, and is kept below: a crossing object under a uniform
# field must read ρ, and a line genuinely below the reference density must contract rather than clip.
# What is retired is the specific arithmetic that could break it.


def test_a_crossing_object_under_a_uniform_field_reads_RHO(multiexon_index):
    """⭐ The surviving half of TRAPS: prefer-shares-to-differences. A line's mass over its own support is the true density — exactly,
    with no factor to get wrong. Verified end-to-end at 1.994 / 2.002 / 1.981 × truth when the old
    average-vs-sum defect was live, so this is where that regression would reappear."""
    ra = RegionArrays.from_index(multiexon_index)
    rho = 0.037
    cal = _field_cal(ra, np.full(ra.n_regions, rho))
    # ⭐ Read straight off the EDGE axis. It used to go through `_left_keyed_edge_arrays`, a node-shaped
    # copy that existed only because the incidence helper emitted a left-node index; that helper now
    # emits an edge index and the copy is deleted, so a line's density is read where it lives.
    edge_mass = np.asarray(cal.mass_gdna_edge, dtype=np.float64)
    edge_support = np.asarray(cal.gdna_edge_eff_len, dtype=np.float64)
    live = edge_support > 0.0
    assert live.any(), "the fixture produced no lines"
    np.testing.assert_allclose(edge_mass[live] / edge_support[live], rho, rtol=1e-12)


def test_a_line_below_the_reference_density_CONTRACTS_rather_than_clipping(multiexon_index):
    """⭐ TRAPS: prefer-shares-to-differences' other half, kept because the ``min(ρ/ρ_ref, 1)`` clip is still there and still hides
    anything that reads a density too HIGH.

    Put the true line density strictly inside ``(ρ_ref/2, ρ_ref)``: read correctly it is below the
    reference and must contract; read at any inflated multiple it lands above, clips to 1, and
    contributes no contraction at all — silent on exactly the lines the shrinkage exists to act on.

    ⚠ The ρ_ref anchor sits on the LAST node, not the first. A line takes its LEFT flank's density
    (:func:`_cal`), so anchoring on node 0 would put line 0 itself at ρ_ref and the assertion would
    fail on the fixture rather than on the code.
    """
    ra = RegionArrays.from_index(multiexon_index)
    rho_ref, rho_line = 1.0, 0.7  # 0.7 ∈ (0.5, 1.0)
    dens = np.full(ra.n_regions, rho_line)
    dens[-1] = rho_ref  # the last node anchors ρ_ref and is no line's LEFT flank
    cal = _field_cal(ra, dens)

    edge_mass = np.asarray(cal.mass_gdna_edge, dtype=np.float64)
    edge_support = np.asarray(cal.gdna_edge_eff_len, dtype=np.float64)
    band = (edge_support > 0.0) & (edge_mass > 0.0)
    assert band.any()
    ratio = (edge_mass[band] / edge_support[band]) / rho_ref
    assert np.all(ratio < 1.0 - 1e-9), (
        f"line density reads {ratio.max():.3f}×ρ_ref — at or above the reference it CLIPS and "
        "contributes no contraction, which is precisely how the factor-2 stayed invisible"
    )


# --- the boundary axis is an EDGE index, and only a MULTI-reference index can prove it -------------

_TWO_REF_GTF = (
    "".join(
        f'chrA\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "ga"; transcript_id "ta";\n'
        for s in (200, 800)
    )
    + "".join(
        f'chrB\ttest\texon\t{s + 1}\t{s + 400}\t.\t+\t.\tgene_id "gb"; transcript_id "tb";\n'
        for s in (200, 800)
    )
)


@pytest.fixture(scope="module")
def two_ref_index(tmp_path_factory):
    return build_test_index(
        tmp_path_factory, _TWO_REF_GTF, name="tworef", refs={"chrA": 2000, "chrB": 2000}
    )


def test_the_boundary_incidence_is_an_EDGE_index_not_a_left_node_index(two_ref_index):
    """⭐⭐ **THE GATE A SINGLE-REFERENCE FIXTURE CANNOT PROVIDE, AND THE PERTURBATION THAT FOUND IT.**

    ``node_right_edge`` numbers edges over adjacent same-reference node pairs, so on ONE reference
    ``edge(r) == r`` and a left-node index is indistinguishable from an edge index. Every existing
    fixture here is single-reference, so substituting one for the other changed nothing and the
    conversion was effectively untested — found by injecting exactly that substitution
    (``TRAPS: perturb-every-gate``).

    ⛔ On a second reference the two axes diverge by one per preceding reference boundary, and indexing
    a per-edge array with a node index then reads **the wrong line's mass** — silently, since both are
    in range. This pins the axis: every emitted boundary index must be a valid edge whose flanking
    nodes are the ones the transcript actually crosses.
    """
    ra = RegionArrays.from_index(two_ref_index)
    _rt, _rr, bt, br, *_ = _transcript_node_incidence(two_ref_index, ra)
    lo, hi = edge_node_indices(np.asarray(ra.ref_id))
    assert br.size, "the fixture produced no interior boundaries"
    assert br.max() < lo.shape[0], "a boundary index outside the edge axis"

    # ⭐ the discriminating claim: each emitted edge's flanks are same-reference neighbours, and the
    # transcript that emitted it overlaps BOTH of them.
    ref_id = np.asarray(ra.ref_id)
    np.testing.assert_array_equal(ref_id[lo[br]], ref_id[hi[br]])
    reg_t, reg_r = _rt, _rr
    for t, e in zip(bt.tolist(), br.tolist()):
        owned = set(reg_r[reg_t == t].tolist())
        assert lo[e] in owned and hi[e] in owned, (
            f"transcript {t} was given line {e} between nodes {lo[e]},{hi[e]} — which it does not span"
        )

    # ...and the two axes genuinely differ here, so the assertions above are not vacuous
    right_edge = node_right_edge(ref_id)
    assert not np.array_equal(right_edge[: lo.shape[0]], np.arange(lo.shape[0])), (
        "fixture is degenerate: the node and edge axes coincide, so this test proves nothing"
    )
