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
    # 3 regions fully inside one locus → φ = 1 each. alpha_gdna = Σ G_r,
    # alpha_rna = Σ D_r, gdna_eff_len = the IPR (Σ G)² / Σ(G²/L) over region lengths.
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5],
        mass_d=[3.0, 4.0, 5.0],
        gdna_geom_len=[100.0, 200.0, 150.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.rna_prior_count, [12.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    # gdna_eff_len is the Laplace-smoothed IPR over the FL-aware support gdna_geom_len =
    # [100, 200, 150]; here the gDNA density g/geom = 0.01 is uniform, so the IPR equals the
    # geometric span Σ gdna_geom_len = 450 (uniform mass → full span; the smoothing is exact).
    np.testing.assert_allclose(priors.gdna_eff_len, [450.0])


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
