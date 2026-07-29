"""assemble_priors — acyclic CalibrationResult → per-locus EM prior.

The gDNA component eff-length is the density-correct, transport-free region + pooled-seam IPR
(``priors._gdna_region_node_arrays``): per-region CONTAINED node at effective support
``S_r = E[max(0,L−ℓ)] = gdna_region_eff_len`` + a POOLED SEAM node per internal boundary at the
averaged per-side density support ``S_s = ½·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])``. The bedrock invariant these tests pin: under a UNIFORM
gDNA field (every node mass = ρ·S_n) the Laplace-smoothed IPR returns ``span = ΣS_n`` exactly — an
unenriched library contracts NOTHING (factor 1). Using the genomic ``region_size_bp`` as the divisor
instead would understate short-region density and fabricate a contraction; the dedicated tests below
prove the method uses the EFFECTIVE support, not the genomic length.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import _gdna_region_node_arrays, assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _calibration(
    *, mass_g, mass_d, region_eff_len, fl_mean=150.0, gdna_density_global=0.01, mass_spliced=None
) -> CalibrationResult:
    """Build a CONTAINED-ONLY result (boundary sides = 0).

    ``mass_g`` / ``mass_d`` are the per-region deconvolved gDNA / RNA mass; ``region_eff_len`` is the
    per-region contained effective support ``E[max(0,L−ℓ)]`` (the density-correct IPR divisor);
    ``fl_mean`` seeds the per-side density length ``gdna_boundary_len`` (the seam support is the average
    of the two flanking values). ``gdna_density_global`` is required but NOT read by ``assemble_priors``.
    ``mass_spliced`` (default 0) is the spliced part of ``mass_d`` withheld from ``rna_prior_count``.
    """
    mg = np.asarray(mass_g, dtype=np.float64)
    md = np.asarray(mass_d, dtype=np.float64)
    rel = np.asarray(region_eff_len, dtype=np.float64)
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
        # ⚠ D6 (2026-07-29): the STORED per-side density length is E[min(ℓ,L)]/2, so for a long region
        # it is fl_mean/2 and a seam's support is the SUM of its two flanks = fl_mean — one whole
        # crossing opportunity per seam, which is exactly v5 §10.3's `ρ = E[flux]/fl_mean`.
        gdna_boundary_len=np.full(n, fl_mean / 2.0, dtype=np.float64),
        gdna_region_eff_len=rel.copy(),
        gdna_density_global=gdna_density_global,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        config=CalibrationConfig(),
    )


def _uniform_field_cal(region_eff_len, boundary_len, rho) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field over a single same-reference region chain, built per the TRUE
    accumulator deposition law: every node has density exactly ``rho``. ``contained[r] = rho·region_eff_len[r]``;
    each region's two boundary sides carry ``rho·boundary_len[r]/2`` (where ``boundary_len[r] = E[min(ℓ,L_r)]``
    is the per-side density length), so the POOLED seam(r,r+1) = ``rho·(boundary_len[r]+boundary_len[r+1])/2``
    = ``rho·S_s``. The factor-1-under-uniform invariant ⇒ ``gdna_eff_len = span = Σ region_eff_len +
    Σ_seam ½(boundary_len[r]+boundary_len[r+1])`` EXACTLY, and ``G/eff_len`` recovers the true ``rho``.

    Using ``boundary_len`` (NOT a constant ``E[ℓ]``) for the side masses is what makes this test a genuine
    guard on the seam-support choice: with VARIED / short ``boundary_len`` the field is uniform only if the
    helper divides the seam by the AVERAGED side density length — an ``E[ℓ]`` divisor would drift off 1."""
    rel = np.asarray(region_eff_len, dtype=np.float64)
    bl = np.asarray(boundary_len, dtype=np.float64) / 2.0  # E[min(ℓ,L)] → the STORED E[min(ℓ,L)]/2
    n = rel.shape[0]
    contained = rho * rel
    side_right = np.zeros(n, dtype=np.float64)
    side_left = np.zeros(n, dtype=np.float64)
    if n > 1:
        # ⚠ D6 (2026-07-29): ``gdna_boundary_len`` IS the halved per-side density length
        # ``E[min(ℓ,L)]/2``, and each face carries ``ρ·gdna_boundary_len``. This fixture used to store
        # the UN-halved ``E[min(ℓ,L)]`` and deposit ``ρ·bl/2``, which put the same mass on each face but
        # a DOUBLED length in ``gdna_boundary_len`` — exactly cancelling the spurious ½ that
        # ``_pooled_seam_arrays`` then applied, and hiding the factor-2 from every assertion here.
        side_right[:-1] = rho * bl[:-1]  # region r's right FACE ⇒ ρ·gdna_boundary_len[r]
        side_left[1:] = rho * bl[1:]  # region r+1's left FACE ⇒ ρ·gdna_boundary_len[r+1]
    return CalibrationResult(
        mass_gdna_contained=contained,
        mass_rna_contained=np.zeros(n, dtype=np.float64),
        mass_gdna_left=side_left,
        mass_rna_left=np.zeros(n, dtype=np.float64),
        mass_gdna_right=side_right,
        mass_rna_right=np.zeros(n, dtype=np.float64),
        mass_rna_spliced=np.zeros(n, dtype=np.float64),
        gdna_boundary_len=bl.copy(),
        gdna_region_eff_len=rel.copy(),
        gdna_density_global=rho,
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


# --- the bedrock invariant: factor = 1 under uniform gDNA -----------------------------------------


def test_factor_one_under_uniform_gdna():
    # THE correctness criterion. A uniform (unenriched) gDNA field — every node density = ρ — must
    # contract NOTHING: gdna_eff_len = span = Σ S_n EXACTLY, and the gDNA per-position rate G/eff_len
    # recovers the true ρ. region_eff_len=[120,200,80], per-side density length boundary_len=[180,60,180]
    # (region 1 is SHORT — its flanks deposit only ~60, not E[ℓ]≈180), ρ=0.02 over 3 same-ref regions.
    # Seam supports are the AVERAGES (180+60)/2=120 and (60+180)/2=120 ⇒ span = (120+200+80)+120+120 = 640.
    # The short middle region is the case an E[ℓ] seam divisor would get wrong (it would over-state the
    # support and drift the factor below 1) — so this pins the deposition-faithful averaged support.
    rel = [120.0, 200.0, 80.0]
    bl = [180.0, 60.0, 180.0]
    rho = 0.02
    span_eff = sum(rel) + 0.5 * (bl[0] + bl[1]) + 0.5 * (bl[1] + bl[2])  # 640
    cal = _uniform_field_cal(rel, bl, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 400)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span_eff], rtol=1e-9)
    # per-position gDNA rate = the true density ρ (the deconvolved gДНК mass spread over the eff-len)
    np.testing.assert_allclose(priors.gdna_prior_count[0] / priors.gdna_eff_len[0], rho, rtol=1e-9)


def test_factor_one_holds_for_any_density():
    # The factor-1 identity is exact for ANY ρ (the Laplace term cancels algebraically), so a 50000×
    # denser uniform library still contracts nothing. Guards against a ρ-dependent contraction.
    rel = [300.0, 150.0]
    bl = [200.0, 200.0]
    span_eff = sum(rel) + 0.5 * (bl[0] + bl[1])  # 650
    ra = _regions([0, 300], [300, 450])
    for rho in (1e-4, 0.01, 0.5, 5.0):
        priors = assemble_priors(_uniform_field_cal(rel, bl, rho), ra, [_ml(0, [(0, 0, 450)])])
        np.testing.assert_allclose(priors.gdna_eff_len, [span_eff], rtol=1e-9)


def test_eff_len_uses_effective_support_not_genomic_size():
    # PROOF the divisor is the EFFECTIVE support gdna_region_eff_len, NOT the genomic region_size_bp.
    # Genomic sizes are all 100 (Σ=300), but the contained effective support is [120,200,80] (Σ=400).
    # Build a uniform field against the EFFECTIVE support; the eff-len must equal the effective span,
    # never the genomic-based 300 + seams. If the method still used region_size_bp, the field would NOT
    # be uniform in its eyes and the factor would drift off 1.
    rel = [120.0, 200.0, 80.0]
    bl = [150.0, 150.0, 150.0]
    rho = 0.03
    span_eff = sum(rel) + 0.5 * (bl[0] + bl[1]) + 0.5 * (bl[1] + bl[2])  # 400 + 150 + 150 = 700
    cal = _uniform_field_cal(rel, bl, rho)
    ra = _regions([0, 100, 200], [100, 200, 300])  # genomic sizes 100,100,100 (Σ=300) ≠ rel
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span_eff], rtol=1e-9)
    assert not np.isclose(
        priors.gdna_eff_len[0], 300.0 + 0.5 * (bl[0] + bl[1]) + 0.5 * (bl[1] + bl[2])
    )


def test_node_arrays_uniform_density_is_constant():
    # The shared node-model helper itself: under a uniform field every node's density m_n/S_n = ρ
    # (regions AND pooled seams), the precondition for the per-node min() factor-1 identity. The varied/
    # short boundary_len exercises the averaged seam support (an E[ℓ] divisor would break the constancy).
    rel = np.array([120.0, 200.0, 80.0])
    bl = [180.0, 60.0, 180.0]
    rho = 0.02
    cal = _uniform_field_cal(rel, bl, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    gdna_region, support_len, pooled, seam_len = _gdna_region_node_arrays(cal, ra)
    # every node sits at density ρ — the folded region and the pooled-seam node alike:
    np.testing.assert_allclose(gdna_region / support_len, rho, rtol=1e-9)
    seam = seam_len > 0
    np.testing.assert_allclose(pooled[seam] / seam_len[seam], rho, rtol=1e-9)


# --- mass / projection (independent of the support choice) ----------------------------------------


def test_single_locus_sums_region_nodes():
    # gdna_prior_count = Σ contained gDNA, rna_prior_count = Σ RNA. These are MASS sums, unaffected by
    # the (effective-support) eff-len divisor. The eff-len magnitude is pinned by the factor-1 tests.
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5], mass_d=[3.0, 4.0, 5.0], region_eff_len=[100.0, 200.0, 150.0]
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.rna_prior_count, [12.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    # eff-len is bounded by the effective span and never expands beyond it:
    span_eff = (100.0 + 200.0 + 150.0) + 2 * 150.0  # Σ region_eff_len + 2 internal seams @ fl_mean
    assert 0.0 < priors.gdna_eff_len[0] <= span_eff + 1e-9


def test_gdna_mass_conservation_contained_plus_sides():
    # CONSERVATION OF MASS (pooled-seam node set). total gDNA = Σ contained + Σ pooled seam, where each
    # internal seam pools the boundary's two halves: seam(r,r+1) = mass_gdna_right[r] + mass_gdna_left[r+1].
    # Every NON-terminal boundary side is counted exactly once via its seam; terminal sides (region 0's
    # left, region 2's right) carry zero on real data and are set zero here. NO mass is moved.
    #   contained = [2,3,1] (Σ=6); seams: (r0,r1)=right[0]+left[1]=1+1=2; (r1,r2)=right[1]+left[2]=2+1=3
    #   ⇒ Σ seam = 5 ⇒ total gDNA = 11.
    cal = CalibrationResult(
        mass_gdna_contained=np.array([2.0, 3.0, 1.0]),
        mass_rna_contained=np.array([0.0, 0.0, 0.0]),
        mass_gdna_left=np.array([0.0, 1.0, 1.0]),  # left[0] = terminal (zero)
        mass_rna_left=np.array([0.0, 0.0, 0.0]),
        mass_gdna_right=np.array([1.0, 2.0, 0.0]),  # right[2] = terminal (zero)
        mass_rna_right=np.array([0.0, 0.0, 0.0]),
        mass_rna_spliced=np.array([0.0, 0.0, 0.0]),
        gdna_boundary_len=np.array([25.0, 25.0, 25.0]),  # D6: E[min]/2 ⇒ seam support 50
        gdna_region_eff_len=np.array([100.0, 100.0, 100.0]),
        gdna_density_global=0.01,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=3,
        config=CalibrationConfig(),
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # all 3 regions in one locus (φ=1) ⇒ gdna_prior_count = the total conserved gДНК mass.
    np.testing.assert_allclose(priors.gdna_prior_count, [11.0])
    # cross-check the helper's mass conservation directly: Σ gdna_region = Σ contained + Σ internal sides.
    gdna_region, _, _, _ = _gdna_region_node_arrays(cal, ra)
    internal_sides = (
        cal.mass_gdna_right[:-1].sum() + cal.mass_gdna_left[1:].sum()
    )  # every non-terminal side once
    np.testing.assert_allclose(gdna_region.sum(), cal.mass_gdna_contained.sum() + internal_sides)


def test_spliced_mass_withheld_from_rna_prior():
    # A spliced fragment has no gDNA candidate in the EM (gDNA does not splice) → it is guaranteed-RNA
    # and assigned directly, so it must NOT load rna_prior_count. RNA mass [3,4,5] (Σ=12) of which
    # spliced [1,1,2] (Σ=4) → rna_prior = the unspliced remainder 8; the gDNA prior is unchanged.
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5],
        mass_d=[3.0, 4.0, 5.0],
        region_eff_len=[100.0, 200.0, 150.0],
        mass_spliced=[1.0, 1.0, 2.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.rna_prior_count, [8.0])  # 12 total RNA − 4 spliced
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])  # unchanged by the spliced withhold


def test_region_split_between_two_loci():
    # One region [0,100) straddling two adjacent loci ([0,50), [50,100)) → overlap shares 0.5/0.5;
    # each locus gets half of every projected quantity. Single region ⇒ no seam ⇒ span = region_eff_len.
    cal = _calibration(mass_g=[5.0], mass_d=[10.0], region_eff_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0, 5.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [2.5, 2.5])
    # uniform single-region gDNA (5 over support 100) ⇒ each half-locus eff-len = 50 (factor 1).
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _calibration(mass_g=[1.0, 50.0], mass_d=[5.0, 99.0], region_eff_len=[100.0, 100.0])
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.gdna_prior_count, [1.0])
    assert priors.gdna_eff_len[0] > 0.0  # finite, positive (magnitude pinned by the factor-1 tests)


# --- Laplace shrinkage toward the (effective) span ------------------------------------------------


def test_evidence_free_region_gives_zero_gdna_prior():
    # Acyclic: no observed gDNA ⇒ zero gDNA pseudocount. With G=0 the Laplace-smoothed IPR is
    # (0+1)²/(1/span) = span exactly, so the eff-len is the uniform effective span (single region ⇒
    # span = region_eff_len = 100) — never a tiny length.
    cal = _calibration(mass_g=[0.0], mass_d=[0.0], region_eff_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])  # effective-span fallback


def _six_region_ra():
    return _regions([0, 100, 200, 300, 400, 500], [100, 200, 300, 400, 500, 600])


def _global_bimodal_cal(rna0: float, gdna0: float = 1.0) -> CalibrationResult:
    """6-region calibration for the CONTAINED-EVIDENCE SHRINKAGE under the GLOBAL reference. Five ENRICHED
    background regions (gDNA density 1.0) fix ρ_ref at the enriched mode; region 0 is DEPLETED (density
    gdna0/100 ≪ ρ_ref) so it contracts (elen ≪ span). ``rna0`` tunes the contained evidence C = gdna0 + rna0
    that drives the shrinkage weight w = C/(C+1). No crossing seam; region_eff_len = 100/region."""
    ng = 6
    mg = np.full(ng, 100.0)
    mg[0] = gdna0
    mr = np.zeros(ng)
    mr[0] = rna0
    z = np.zeros(ng)
    return CalibrationResult(
        mass_gdna_contained=mg,
        mass_rna_contained=mr,
        mass_gdna_left=z.copy(),
        mass_rna_left=z.copy(),
        mass_gdna_right=z.copy(),
        mass_rna_right=z.copy(),
        mass_rna_spliced=z.copy(),
        gdna_boundary_len=np.full(ng, 25.0),
        gdna_region_eff_len=np.full(ng, 100.0),
        gdna_density_global=0.5,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=ng,
        config=CalibrationConfig(),
    )


def test_eff_len_shrinks_toward_span_for_sparse_gdna():
    # Under the GLOBAL reference a DEPLETED locus (density ≪ ρ_ref, set by the enriched background) contracts
    # (elen ≪ span); the contained-evidence shrinkage pulls that contraction back toward span when the
    # unique-mapper evidence is sparse (abundant ⇒ the earned contraction, sparse ⇒ near span).
    ra = _six_region_ra()
    ml = [_ml(0, [(0, 0, 100)])]  # locus = region 0 (the depleted one)
    eff_sparse = assemble_priors(_global_bimodal_cal(rna0=0.3), ra, ml).gdna_eff_len[0]
    eff_abundant = assemble_priors(_global_bimodal_cal(rna0=100.0), ra, ml).gdna_eff_len[0]
    assert eff_abundant < eff_sparse  # more contained evidence ⇒ more of the earned contraction


def _blind_seam_cal(contained_rna: float) -> CalibrationResult:
    """A locus whose gDNA contraction is driven by a FIXED pooled SEAM (no contained gDNA), with a
    tunable amount of CONTAINED RNA — to isolate the contained-evidence shrinkage (the seam fixes the
    IPR; the contained RNA drives the shrinkage weight w = C/(C+1)). region_eff_len=100/region,
    boundary_len=50/region ⇒ seam support = ½(50+50) = 50 each ⇒ effective span = Σ region_eff_len +
    2 seams = 300 + 2·50 = 400."""
    return CalibrationResult(
        mass_gdna_contained=np.array(
            [0.0, 0.0, 0.0]
        ),  # no contained gDNA — all signal is seam smear
        mass_rna_contained=np.array(
            [contained_rna, 0.0, 0.0]
        ),  # the only contained (unique-mapper) evidence
        mass_gdna_left=np.array([0.0, 1.0, 1.0]),  # pooled seam = right[r]+left[r+1] = [2, 3, 0]
        mass_rna_left=np.array([0.0, 0.0, 0.0]),
        mass_gdna_right=np.array([1.0, 2.0, 0.0]),
        mass_rna_right=np.array([0.0, 0.0, 0.0]),
        mass_rna_spliced=np.array([0.0, 0.0, 0.0]),
        gdna_boundary_len=np.array([25.0, 25.0, 25.0]),  # D6: E[min]/2 ⇒ seam support 50
        gdna_region_eff_len=np.array([100.0, 100.0, 100.0]),
        gdna_density_global=0.01,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=3,
        config=CalibrationConfig(),
    )


def test_contained_evidence_shrinkage_reverts_to_span_when_blind():
    # C = 0 (multimapper-blind locus: zero contained mass, only seam smear): w = 0 ⇒ eff_len → span (no
    # contraction), the smooth shrinkage's C→0 limit. The prior is honestly uninformative where
    # calibration cannot see. The seam mass is still counted in the prior pseudocount.
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(_blind_seam_cal(0.0), ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [400.0])  # effective span = 300 + 2·50
    np.testing.assert_allclose(priors.gdna_prior_count, [5.0])  # pooled seams 2 + 3


def test_contained_evidence_shrinkage_is_smooth_not_a_cliff():
    # The shrinkage is SMOOTH in contained evidence (not a hard cliff). Under the global reference a DEPLETED
    # locus contracts (elen ≪ span); increasing the CONTAINED evidence C = gdna+rna monotonically raises the
    # trust w = C/(C+1) and pulls the eff-length from span (C→0) down toward the earned contraction (C≫1).
    # Counts 0,1,3,1000 interpolate strictly monotonically — no cliff.
    ra = _six_region_ra()
    ml = [
        _ml(0, [(0, 0, 100)])
    ]  # locus = region 0 (depleted vs the enriched background ⇒ contracts)

    def eff(c):
        return assemble_priors(_global_bimodal_cal(rna0=c), ra, ml).gdna_eff_len[0]

    e0, e1, e3, e_hi = eff(0.0), eff(1.0), eff(3.0), eff(1000.0)
    assert e_hi < e3 < e1 < e0  # smooth + monotone (no cliff)
    assert (
        e0 < 150.0
    )  # even C=0 carries the depleted contained gDNA as evidence ⇒ below the full span


def test_empty_multiloci_returns_empty():
    cal = _calibration(mass_g=[1.0], mass_d=[1.0], region_eff_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [])
    assert priors.rna_prior_count.shape == (0,)
    assert priors.gdna_prior_count.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _calibration(mass_g=[1.0, 2.0], mass_d=[1.0, 2.0], region_eff_len=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])


def test_gdna_eff_len_factor_one_under_uniform_gdna_with_kde_firing():
    # Priors-side factor-1: 6 regions at UNIFORM gDNA density (contained = ρ·region_eff, pooled seam =
    # ρ·seam_support) ⇒ the KDE detector FIRES (≥6 gDNA nodes) but is unimodal ⇒ ρ_ref = ρ ⇒ every node
    # min(m/ρ_ref, S) = S ⇒ gdna_eff_len == span EXACTLY (contraction factor 1). The existing gDNA-eff tests
    # run in the <5-region regime that returns None, so this locks the KDE code path on the priors side.
    from rigel.calibration.capture_eff_length import _global_reference_density

    ra = _six_region_ra()
    ml = [_ml(0, [(0, 0, 600)])]
    rho = 2.0
    cal = CalibrationResult(
        mass_gdna_contained=np.full(6, rho * 100.0),  # ρ·region_eff_len ⇒ contained density = ρ
        mass_rna_contained=np.zeros(6),
        mass_gdna_left=np.full(
            6, 50.0
        ),  # pooled seam m_s = right[r]+left[r+1] = 100 = ρ·seam_support(=50)
        mass_rna_left=np.zeros(6),
        mass_gdna_right=np.full(6, 50.0),
        mass_rna_right=np.zeros(6),
        mass_rna_spliced=np.zeros(6),
        gdna_boundary_len=np.full(6, 25.0),
        gdna_region_eff_len=np.full(6, 100.0),
        gdna_density_global=0.5,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=6,
        config=CalibrationConfig(),
    )
    # the KDE fires (unimodal ⇒ ρ_ref = ρ), NOT the <5-node None fallback
    assert _global_reference_density(
        cal.mass_gdna_contained, cal.gdna_region_eff_len
    ) == pytest.approx(rho)
    # span = Σ support_len = 5·(100 + ½·(50+50)) + 1·100 (terminal region carries no right seam)
    span = 5 * 150.0 + 100.0
    np.testing.assert_allclose(assemble_priors(cal, ra, ml).gdna_eff_len[0], span, rtol=1e-9)
