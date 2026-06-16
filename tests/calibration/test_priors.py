"""assemble_priors — acyclic CalibrationResult → per-locus EM prior."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import _transport_boundary_flux, assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _calibration(
    *, mass_g, mass_d, gdna_geom_len, gdna_density_global=0.01, mass_spliced=None
) -> CalibrationResult:
    """Build a result with all mass in the contained node (left/right = 0).

    ``mass_g`` / ``mass_d`` are the per-region deconvolved gDNA / RNA mass;
    ``gdna_geom_len`` is the gDNA component's geometric length. ``gdna_density_global``
    is a required field but is **not** read by ``assemble_priors`` (it is already
    baked into the deconvolved mass). ``mass_spliced`` (default 0) is the spliced part of
    ``mass_d`` that ``assemble_priors`` withholds from ``rna_prior_count``.
    """
    mg = np.asarray(mass_g, dtype=np.float64)
    md = np.asarray(mass_d, dtype=np.float64)
    el = np.asarray(gdna_geom_len, dtype=np.float64)
    n = mg.shape[0]
    z = np.zeros(n, dtype=np.float64)
    ms = z.copy() if mass_spliced is None else np.asarray(mass_spliced, dtype=np.float64)
    return CalibrationResult(
        mass_gdna_contained=mg,
        mass_rna_contained=md,
        mass_gdna_left=z.copy(),
        mass_rna_left=z.copy(),
        mass_gdna_right=z.copy(),
        mass_rna_right=z.copy(),
        mass_rna_spliced=ms,
        gdna_geom_len=el,
        gdna_boundary_len=el,
        gdna_density_global=gdna_density_global,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )


def _regions(starts, ends) -> RegionArrays:
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    n = starts.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        signature=np.zeros(n, dtype=np.uint8),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=np.array([0, n], dtype=np.int32),
        order=np.arange(n, dtype=np.int64),
        n_refs=1,
    )


def _ml(locus_id, blocks) -> MultiLocus:
    loci = tuple(Locus(ref=str(rid), ref_id=rid, start=s, end=e) for rid, s, e in blocks)
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        gdna_span=sum(e - s for _, s, e in blocks),
        loci=loci,
    )


def test_single_locus_sums_region_nodes():
    # 3 regions fully inside one locus → φ = 1 each. gdna_prior_count = Σ G_r, rna_prior_count = Σ D_r,
    # gdna_eff_len = the Laplace-smoothed IPR (Σ G)²/Σ(G²/L) over the GENOMIC region lengths
    # region_size_bp (conservation-correct, PR-1). Region sizes [100, 200, 150] with gDNA mass ∝ size
    # ⇒ uniform per-bp density ⇒ the IPR equals the genomic span Σ region_size_bp = 450 (uniform mass →
    # full span; the smoothing is exact). gdna_geom_len is no longer read by assemble_priors.
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5],
        mass_d=[3.0, 4.0, 5.0],
        gdna_geom_len=[100.0, 200.0, 150.0],  # unused here now (capture_eff_length consumes it)
    )
    ra = _regions([0, 100, 300], [100, 300, 450])  # genomic sizes 100, 200, 150 (Σ = 450 bp)
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.rna_prior_count, [12.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [450.0])


def test_gdna_length_conservation_span_is_genomic():
    # CONSERVATION OF GENOMIC LENGTH (the PR-1 regression guard). The gDNA component's IPR span must be
    # the genomic span Σ region_size_bp, NOT the old FL-aware gdna_geom_len (which double-counts each
    # boundary seam + adds mean-FL per region ⇒ ~8% inflation, ≈2× on small exons). Here gdna_geom_len
    # is set 3× the genomic size to expose any leak of it into the length: with uniform per-bp gDNA the
    # eff-len must equal the GENOMIC span (300), never the inflated 900.
    cal = _calibration(
        mass_g=[1.0, 1.0, 1.0],
        mass_d=[0.0, 0.0, 0.0],
        gdna_geom_len=[300.0, 300.0, 300.0],  # 3× the genomic size — must NOT be used as the length
    )
    ra = _regions([0, 100, 200], [100, 200, 300])  # genomic sizes 100, 100, 100 (Σ = 300 bp)
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # uniform per-bp gDNA over genomic bp ⇒ eff_len = genomic span = 300 (the conserved length),
    # AND the cap is the genomic span: eff_len must never exceed Σ region_size_bp.
    np.testing.assert_allclose(priors.gdna_eff_len, [300.0])
    assert priors.gdna_eff_len[0] <= 300.0 + 1e-9


def test_gdna_mass_conservation_contained_plus_sides():
    # CONSERVATION OF MASS (pooled-seam node set, PR-A). The total gDNA = Σ contained + Σ pooled seam,
    # where each internal seam pools the boundary's two halves: seam(r,r+1) = mass_gdna_right[r] +
    # mass_gdna_left[r+1]. Every NON-terminal boundary side is counted exactly once via its seam; terminal
    # sides (region 0's left, region 2's right) carry zero mass on real data and are set zero here.
    #   contained = [2,3,1] (Σ=6); internal seams: (r0,r1) = right[0]+left[1] = 1+1 = 2;
    #   (r1,r2) = right[1]+left[2] = 2+1 = 3 ⇒ Σ seam = 5 ⇒ total gDNA = 11.
    cal = CalibrationResult(
        mass_gdna_contained=np.array([2.0, 3.0, 1.0]),
        mass_rna_contained=np.array([0.0, 0.0, 0.0]),
        mass_gdna_left=np.array([0.0, 1.0, 1.0]),  # left[0] = terminal (zero)
        mass_rna_left=np.array([0.0, 0.0, 0.0]),
        mass_gdna_right=np.array([1.0, 2.0, 0.0]),  # right[2] = terminal (zero)
        mass_rna_right=np.array([0.0, 0.0, 0.0]),
        mass_rna_spliced=np.array([0.0, 0.0, 0.0]),
        gdna_geom_len=np.array([100.0, 100.0, 100.0]),
        gdna_boundary_len=np.array([50.0, 50.0, 50.0]),
        gdna_density_global=0.01,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=3,
        config=CalibrationConfig(),
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # all 3 regions in one locus (φ=1) ⇒ gdna_prior_count = the total conserved gDNA mass.
    np.testing.assert_allclose(priors.gdna_prior_count, [11.0])


def _blind_seam_cal(contained_rna: float) -> CalibrationResult:
    """A locus whose gDNA contraction is driven by a FIXED pooled SEAM (no contained gDNA), with a tunable
    amount of CONTAINED RNA — to isolate the contained-evidence shrinkage (the seam fixes the IPR; the
    contained RNA drives the shrinkage weight w = C/(C+1))."""
    return CalibrationResult(
        mass_gdna_contained=np.array([0.0, 0.0, 0.0]),  # no contained gDNA — all signal is seam smear
        mass_rna_contained=np.array([contained_rna, 0.0, 0.0]),  # the only contained (unique-mapper) evidence
        mass_gdna_left=np.array([0.0, 1.0, 1.0]),  # seam: pooled = right[r]+left[r+1] = [2, 3, 0]
        mass_rna_left=np.array([0.0, 0.0, 0.0]),
        mass_gdna_right=np.array([1.0, 2.0, 0.0]),
        mass_rna_right=np.array([0.0, 0.0, 0.0]),
        mass_rna_spliced=np.array([0.0, 0.0, 0.0]),
        gdna_geom_len=np.array([100.0, 100.0, 100.0]),
        gdna_boundary_len=np.array([50.0, 50.0, 50.0]),  # FL-scale ⇒ the seam contracts the IPR to ≈121
        gdna_density_global=0.01,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=3,
        config=CalibrationConfig(),
    )


def test_contained_evidence_shrinkage_reverts_to_span_when_blind():
    # C = 0 (multimapper-blind locus: zero contained mass, only seam smear): w = 0 ⇒ eff_len → span (no
    # contraction), the smooth shrinkage's C→0 limit. The prior is honestly uninformative where calibration
    # cannot see. The seam mass is still counted in the prior pseudocount (only the eff-len is shrunk).
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(_blind_seam_cal(0.0), ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [300.0])  # span
    np.testing.assert_allclose(priors.gdna_prior_count, [5.0])  # pooled seams 2 + 3


def test_contained_evidence_shrinkage_is_smooth_not_a_cliff():
    # The shrinkage is SMOOTH in contained evidence (not a hard cliff at C=0). With the seam fixing the raw
    # IPR contraction (≈121 ≪ span=300), increasing CONTAINED RNA monotonically increases the trust w and
    # pulls eff_len from span (C=0) down toward the contraction (C≫1). Counts 1,2,3,... interpolate.
    ra = _regions([0, 100, 200], [100, 200, 300])
    span = 300.0

    def eff(c):
        return assemble_priors(_blind_seam_cal(c), ra, [_ml(0, [(0, 0, 300)])]).gdna_eff_len[0]

    e0, e1, e3, e_hi = eff(0.0), eff(1.0), eff(3.0), eff(1000.0)
    # smooth + monotone: span at C=0, strictly decreasing toward the contraction as C grows (no cliff).
    assert e0 == pytest.approx(span)
    assert e_hi < e3 < e1 < e0
    # C=1 lands strictly BETWEEN the contraction and span (a cliff would jump straight to the contraction).
    assert e_hi < e1 < span
    # high contained evidence ⇒ ≈ the earned IPR contraction (the leak win is preserved on real loci).
    # raw IPR = (G+1)²/(supp + (2G+1)/span), G=5 (pooled seams 2+3), supp = (2²+3²)/50 = 0.26.
    assert e_hi == pytest.approx(36.0 / (0.26 + 11.0 / 300.0), rel=0.02)
    assert e_hi < 0.5 * span


def test_spliced_mass_withheld_from_rna_prior():
    # A spliced fragment has no gDNA candidate in the EM (gDNA does not splice) → it is
    # guaranteed-RNA and the EM assigns it directly, so it must NOT load rna_prior_count. Only the
    # UNSPLICED RNA competes with gDNA. RNA mass [3,4,5] (Σ=12) of which spliced [1,1,2] (Σ=4) →
    # rna_prior = the unspliced remainder 8; the gDNA prior is unchanged (spliced never touches it).
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5],
        mass_d=[3.0, 4.0, 5.0],
        gdna_geom_len=[100.0, 200.0, 150.0],
        mass_spliced=[1.0, 1.0, 2.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.rna_prior_count, [8.0])  # 12 total RNA − 4 spliced
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])  # unchanged by the spliced withhold


def test_region_split_between_two_loci():
    # One region [0,100) straddling two adjacent loci ([0,50), [50,100)) →
    # overlap shares 0.5/0.5; each locus gets half of every projected quantity.
    cal = _calibration(mass_g=[5.0], mass_d=[10.0], gdna_geom_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(
        cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])]
    )
    np.testing.assert_allclose(priors.rna_prior_count, [5.0, 5.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [2.5, 2.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_evidence_free_region_gives_zero_gdna_prior():
    # Acyclic: no observed gDNA ⇒ zero gDNA pseudocount (no ρ_0·ω·L hallucination).
    # With no gDNA mass (G=0) the Laplace-smoothed IPR is (0+1)²/(1/span) = span exactly, so
    # the eff-len is the uniform geometric span (gdna_geom_len = 100) — never a tiny length.
    cal = _calibration(mass_g=[0.0], mass_d=[0.0], gdna_geom_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])  # geometric-span fallback


def test_eff_len_shrinks_toward_span_for_sparse_gdna():
    # gDNA concentrated on region 0 (region 1 empty): IPR = 100 (region-0 support),
    # geometric span = 200. The Laplace one-fragment uniform-per-base prior pulls the eff-len
    # from the concentrated IPR toward the uniform span in proportion to the gDNA count: a
    # sparse mass sits near the span, an abundant mass near the IPR. This is what stops the EM
    # from amplifying a tiny concentrated gDNA mass past the calibration's call.
    ra = _regions([0, 100], [100, 200])
    ml = [_ml(0, [(0, 0, 200)])]
    sparse = _calibration(mass_g=[2.0, 0.0], mass_d=[0.0, 0.0], gdna_geom_len=[100.0, 100.0])
    abundant = _calibration(mass_g=[100.0, 0.0], mass_d=[0.0, 0.0], gdna_geom_len=[100.0, 100.0])
    eff_sparse = assemble_priors(sparse, ra, ml).gdna_eff_len[0]
    eff_abundant = assemble_priors(abundant, ra, ml).gdna_eff_len[0]
    assert 100.0 < eff_abundant < eff_sparse < 200.0


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _calibration(
        mass_g=[1.0, 50.0],
        mass_d=[5.0, 99.0],
        gdna_geom_len=[100.0, 100.0],
    )
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.gdna_prior_count, [1.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])


def test_empty_multiloci_returns_empty():
    cal = _calibration(mass_g=[1.0], mass_d=[1.0], gdna_geom_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [])
    assert priors.rna_prior_count.shape == (0,)
    assert priors.gdna_prior_count.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _calibration(mass_g=[1.0, 2.0], mass_d=[1.0, 2.0], gdna_geom_len=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])


# --- boundary-flux transport (_transport_boundary_flux) -----------------------


def test_transport_smear_flows_to_exposed_exon():
    # A short captured exon (high density) beside a long intron with NO contained gDNA
    # (contained = 0) — its only mass is boundary smear. The transport must pull that
    # smear back onto the exon; the intron drains toward zero.
    contained = np.array([100.0, 0.0])  # exon has the real gDNA; intron has none
    left = np.array([0.0, 50.0])  # smear on the intron's exon-facing side
    right = np.array([50.0, 0.0])  # smear on the exon's intron-facing side
    length = np.array([1000.0, 14000.0])  # 1 kb exon, 14 kb intron
    e_cap = np.array([150.0, 150.0])
    ref_id = np.array([0, 0])
    g = _transport_boundary_flux(contained, left, right, length, e_cap, ref_id)
    assert g[0] > 199.0  # exon recaptures essentially all the smear
    assert g[1] < 1.0  # intron drained
    np.testing.assert_allclose(g.sum(), 200.0)  # mass conserved


def test_transport_equal_density_splits_evenly():
    # Two identical regions ⇒ the pooled boundary mass splits 50/50.
    contained = np.array([100.0, 100.0])
    left = np.array([0.0, 50.0])
    right = np.array([50.0, 0.0])
    length = np.array([1000.0, 1000.0])
    e_cap = np.array([150.0, 150.0])
    ref_id = np.array([0, 0])
    g = _transport_boundary_flux(contained, left, right, length, e_cap, ref_id)
    np.testing.assert_allclose(g, [150.0, 150.0])  # each keeps its 100 + half the 100 pool


def test_transport_does_not_cross_reference_boundary():
    # Regions on different references share no genomic boundary → no transport; each
    # region just sums its own three nodes.
    contained = np.array([100.0, 0.0])
    left = np.array([0.0, 50.0])
    right = np.array([50.0, 0.0])
    length = np.array([1000.0, 14000.0])
    e_cap = np.array([150.0, 150.0])
    ref_id = np.array([0, 1])  # cross-reference: boundary is not internal
    g = _transport_boundary_flux(contained, left, right, length, e_cap, ref_id)
    np.testing.assert_allclose(g, [150.0, 50.0])  # untouched node sums
