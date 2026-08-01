"""assemble_priors — acyclic CalibrationResult → per-locus EM prior.

The gDNA component eff-length is the density-correct, transport-free node + line IPR
(``priors._component_node_arrays``): a per-node CONTAINED object at effective support
``S_r = E_f[(L_r − w + 1)+] = gdna_node_eff_len`` plus a per-line CROSSING object at
``S_e = E_f[w − 1] = gdna_edge_eff_len``. The bedrock invariant these tests pin: under a UNIFORM gDNA
field (every object's mass = ρ·S) the Laplace-smoothed IPR returns ``span = ΣS`` exactly — an
unenriched library contracts NOTHING (factor 1). Using the genomic ``region_size_bp`` as the divisor
instead would understate short-node density and fabricate a contraction; the dedicated tests below
prove the method uses the EFFECTIVE support, not the genomic length.

⭐ **THE UNIFORM-FIELD FIXTURE IS NOW THREE LINES, AND THAT IS THE EVIDENCE.** It used to carry a
ten-line note explaining that ``gdna_boundary_len`` was ALREADY the halved per-side density length
``E[min(ℓ,L)]/2``, that each face therefore deposited ``ρ·gdna_boundary_len``, and that an earlier
version of the fixture had stored the UN-halved length while depositing half the mass — cancelling
exactly, and hiding a factor of 2 from every assertion in this file for months
(`CARRY_FORWARD.md` §3 trap 2). A contiguous edge is a 0-bp line with one mass and one support, so a
uniform field is just ``mass = ρ·support`` on both axes and there is no ½ left to get wrong.

⚠ **Every span below is byte-identical to the pre-S5.f value** (640 / 650 / 700 / 400 / 850). The
schema changed; the geometry did not. A number that moved here would mean the port re-derived
something rather than re-keying it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import _component_node_arrays, assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _result(
    *,
    node_g,
    node_r,
    node_eff,
    edge_g=None,
    edge_r=None,
    edge_eff=None,
    edge_spliced=None,
    gdna_density_global=0.01,
    rna_node_eff=None,
    rna_edge_eff=None,
) -> CalibrationResult:
    """Build a result on the three axes. One reference with ``n`` nodes owns exactly ``n − 1`` lines.

    ``edge_*`` default to zeros, so a caller that cares only about contained mass writes only the node
    arrays — but the edge axis is still the RIGHT LENGTH, because an edge axis inconsistent with its
    own node axis is a mis-shaped fixture, not a "no edges" one.

    ⭐ **The RNA supports default to the gDNA ones**, so a test that is about projection, conservation
    or re-keying — not about the length tilt — keeps both components on one support and its ``g:r``
    ratio is unchanged from the pre-P1 fixture. The tilt itself is exercised where it belongs, in
    `test_prior_units.py`, by giving the two components genuinely different opportunities.
    """
    ng = np.asarray(node_g, dtype=np.float64)
    n = ng.shape[0]
    ne = max(n - 1, 0)
    ez = np.zeros(ne, dtype=np.float64)
    node_eff_arr = np.asarray(node_eff, dtype=np.float64)
    edge_eff_arr = ez.copy() if edge_eff is None else np.asarray(edge_eff, dtype=np.float64)
    return CalibrationResult(
        mass_gdna_node=ng,
        mass_rna_node=np.asarray(node_r, dtype=np.float64),
        mass_gdna_edge=ez.copy() if edge_g is None else np.asarray(edge_g, dtype=np.float64),
        mass_rna_edge=ez.copy() if edge_r is None else np.asarray(edge_r, dtype=np.float64),
        mass_rna_spliced_edge=(
            ez.copy() if edge_spliced is None else np.asarray(edge_spliced, dtype=np.float64)
        ),
        mass_rna_junction=np.zeros(0, dtype=np.float64),
        gdna_node_eff_len=node_eff_arr,
        gdna_edge_eff_len=edge_eff_arr,
        rna_node_eff_len=(
            node_eff_arr if rna_node_eff is None else np.asarray(rna_node_eff, dtype=np.float64)
        ),
        rna_edge_eff_len=(
            edge_eff_arr if rna_edge_eff is None else np.asarray(rna_edge_eff, dtype=np.float64)
        ),
        gdna_density_global=gdna_density_global,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_nodes=n,
        n_edges=ne,
        n_junctions=0,
        config=CalibrationConfig(),
    )


def _uniform_field(node_eff, edge_eff, rho) -> CalibrationResult:
    """A genuinely UNIFORM gDNA field: every object's mass is ``ρ × its own effective support``.

    That is the accumulator's deposition law stated directly — ``ρ·E_f[(L−w+1)+]`` contained,
    ``ρ·E_f[w−1]`` crossing — with no faces, no halving and nothing to cancel. The factor-1 invariant
    then says ``gdna_eff_len == span == Σ S`` exactly, and ``G/eff_len`` recovers the true ρ.
    """
    node_eff = np.asarray(node_eff, dtype=np.float64)
    edge_eff = np.asarray(edge_eff, dtype=np.float64)
    return _result(
        node_g=rho * node_eff,
        node_r=np.zeros_like(node_eff),
        node_eff=node_eff,
        edge_g=rho * edge_eff,
        edge_eff=edge_eff,
        gdna_density_global=rho,
    )


def _regions(starts, ends, signature=None) -> RegionArrays:
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    n = starts.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        signature=(
            np.zeros(n, dtype=np.uint8) if signature is None else np.asarray(signature, np.uint8)
        ),
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
    # THE correctness criterion. A uniform (unenriched) gDNA field — every object's density = ρ — must
    # contract NOTHING: gdna_eff_len = span = Σ S EXACTLY, and the gDNA per-position rate G/eff_len
    # recovers the true ρ. node_eff=[120,200,80] (node 1 is SHORT), edge_eff=[120,120], ρ=0.02 over
    # 3 same-ref nodes ⇒ span = 400 + 240 = 640.
    node_eff = [120.0, 200.0, 80.0]
    edge_eff = [120.0, 120.0]
    rho = 0.02
    span = sum(node_eff) + sum(edge_eff)  # 640
    cal = _uniform_field(node_eff, edge_eff, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 400)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    # ⭐ P1: the prior is a FRAGMENT COUNT — the recovered density times the GENOMIC span (400 bp, the
    # three nodes; the two lines are 0-bp and add none). Under a uniform field Σm/ΣS is ρ exactly, so
    # this asserts BOTH that the density is recovered and that the conversion uses the genomic span.
    # ⚠ It replaces `gdna_prior_count / gdna_eff_len == rho`, which held only while the prior was a sum
    # of per-object masses over the same ΣS the eff-len is built from — the units defect P1 removes.
    np.testing.assert_allclose(priors.gdna_prior_count, [rho * 400.0], rtol=1e-9)


def test_factor_one_holds_for_any_density():
    # The factor-1 identity is exact for ANY ρ (the Laplace term cancels algebraically), so a 50000×
    # denser uniform library still contracts nothing. Guards against a ρ-dependent contraction.
    node_eff = [300.0, 150.0]
    edge_eff = [200.0]
    span = 650.0
    ra = _regions([0, 300], [300, 450])
    for rho in (1e-4, 0.01, 0.5, 5.0):
        priors = assemble_priors(
            _uniform_field(node_eff, edge_eff, rho), ra, [_ml(0, [(0, 0, 450)])]
        )
        np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)


def test_eff_len_uses_effective_support_not_genomic_size():
    # PROOF the divisor is the EFFECTIVE support gdna_node_eff_len, NOT the genomic region_size_bp.
    # Genomic sizes are all 100 (Σ=300), but the contained effective support is [120,200,80] (Σ=400).
    # Build a uniform field against the EFFECTIVE support; the eff-len must equal the effective span,
    # never the genomic-based 300 + lines. If the method still used region_size_bp, the field would NOT
    # be uniform in its eyes and the factor would drift off 1.
    node_eff = [120.0, 200.0, 80.0]
    edge_eff = [150.0, 150.0]
    span = 700.0
    cal = _uniform_field(node_eff, edge_eff, 0.03)
    ra = _regions([0, 100, 200], [100, 200, 300])  # genomic sizes 100,100,100 (Σ=300) ≠ node_eff
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [span], rtol=1e-9)
    assert not np.isclose(priors.gdna_eff_len[0], 300.0 + sum(edge_eff))


def test_node_arrays_uniform_density_is_constant():
    # The shared object-model helper itself: under a uniform field every object's density m/S = ρ
    # (contained AND crossing), the precondition for the per-object min() factor-1 identity.
    node_eff = np.array([120.0, 200.0, 80.0])
    edge_eff = np.array([120.0, 120.0])
    rho = 0.02
    cal = _uniform_field(node_eff, edge_eff, rho)
    ra = _regions([0, 120, 320], [120, 320, 400])
    gdna_total, support_total, edge_mass, edge_support = _component_node_arrays(
        cal,
        ra,
        cal.mass_gdna_node,
        cal.mass_gdna_edge,
        cal.gdna_node_eff_len,
        cal.gdna_edge_eff_len,
    )
    np.testing.assert_allclose(gdna_total / support_total, rho, rtol=1e-9)
    live = edge_support > 0
    np.testing.assert_allclose(edge_mass[live] / edge_support[live], rho, rtol=1e-9)


# --- mass / projection (independent of the support choice) ----------------------------------------


def test_single_locus_projects_both_components():
    # ⭐ P1: the priors are DENSITIES INTEGRATED OVER THE SPAN, not mass sums.
    #   ΣS = Σ node_eff (450) + Σ edge_eff (300) = 750 ;  span_bp = 100+200+150 = 450
    #   gDNA: Σm 4.5  ⇒ ρ 0.006 ⇒ 0.006·450 = 2.7      (was 4.5, the raw incidence sum)
    #   RNA : Σm 12.0 ⇒ ρ 0.016 ⇒ 0.016·450 = 7.2      (was 12.0)
    # ⚠ The g:r RATIO is unchanged (0.375) because this fixture gives both components the same support;
    # the ratio moves only when their opportunities differ, which is `test_prior_units.py`'s job.
    cal = _result(
        node_g=[1.0, 2.0, 1.5],
        node_r=[3.0, 4.0, 5.0],
        node_eff=[100.0, 200.0, 150.0],
        edge_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.rna_prior_count, [7.2])
    np.testing.assert_allclose(priors.gdna_prior_count, [2.7])
    np.testing.assert_allclose(
        priors.gdna_prior_count[0] / priors.rna_prior_count[0], 4.5 / 12.0, rtol=1e-9
    )
    span = (100.0 + 200.0 + 150.0) + 2 * 150.0
    assert 0.0 < priors.gdna_eff_len[0] <= span + 1e-9


def test_gdna_mass_conservation_nodes_plus_lines():
    # CONSERVATION OF MASS. total gDNA = Σ node mass + Σ line mass, every object counted exactly once.
    #   nodes = [2,3,1] (Σ=6); lines = [2,3] (Σ=5)  ⇒  total gDNA = 11.
    # ⭐ The predecessor reached 11 by summing FOUR half-arrays — right[0]+left[1] and right[1]+left[2],
    # with the two terminal sides zeroed by hand. There are no halves and no terminals now.
    cal = _result(
        node_g=[2.0, 3.0, 1.0],
        node_r=[0.0, 0.0, 0.0],
        node_eff=[100.0, 100.0, 100.0],
        edge_g=[2.0, 3.0],
        edge_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    # ⭐ OBJECT conservation is what this test is named for and it is UNCHANGED: every object counted
    # exactly once, Σ = 11, asserted on `_component_node_arrays` below. ⚠ The PRIOR is no longer that
    # sum — a fragment deposits on max(K,1) objects, so 11 is an incidence count, not a fragment count:
    #   ΣS = 300 + 100 = 400 ;  span_bp = 300  ⇒  ρ = 11/400 = 0.0275  ⇒  0.0275·300 = 8.25
    np.testing.assert_allclose(priors.gdna_prior_count, [8.25])
    gdna_total, _, _, _ = _component_node_arrays(
        cal,
        ra,
        cal.mass_gdna_node,
        cal.mass_gdna_edge,
        cal.gdna_node_eff_len,
        cal.gdna_edge_eff_len,
    )
    np.testing.assert_allclose(
        gdna_total.sum(), cal.mass_gdna_node.sum() + cal.mass_gdna_edge.sum()
    )


def test_spliced_mass_withheld_from_rna_prior():
    # A spliced fragment has no gDNA candidate in the EM (gDNA does not splice) → it is guaranteed-RNA
    # and assigned directly, so it must NOT load rna_prior_count. Node RNA [3,4,5] (Σ=12) plus line RNA
    # [4,4] of which [1,3] is spliced ⇒ RNA mass = 12 + (4−1) + (4−3) = 16 (NOT 20).
    # ⭐ P1 converts that to a fragment count: ΣS = 450 + 300 = 750, span_bp = 300
    #   RNA : 16/750 · 300 = 6.4   ·   gDNA: 4.5/750 · 300 = 1.8
    # The WITHHOLDING is what this test pins, and it survives the units change intact: without it the
    # RNA mass would be 20 and the prior 8.0.
    cal = _result(
        node_g=[1.0, 2.0, 1.5],
        node_r=[3.0, 4.0, 5.0],
        node_eff=[100.0, 200.0, 150.0],
        edge_r=[4.0, 4.0],
        edge_spliced=[1.0, 3.0],
        edge_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.rna_prior_count, [6.4])
    np.testing.assert_allclose(priors.gdna_prior_count, [1.8])
    assert priors.rna_prior_count[0] < 20.0 / 750.0 * 300.0  # the spliced mass really is withheld


def test_the_junction_flux_does_NOT_enter_the_rna_prior():
    """⭐ The ruling, pinned. A junction fragment is certified RNA in exactly the sense a spliced
    crossing is withheld for — it has no gDNA candidate in the EM — so counting it would load the RNA
    side of a split that arbitrates only UNSPLICED fragments. A locus whose RNA is fully spliced SHOULD
    get a near-zero ``rna_prior_count``: its unspliced fragments really are gDNA or nascent.

    ⚠ The result carries the flux for QC (`test_calibrate`); ``assemble_priors`` must ignore it, and
    that is a deliberate asymmetry rather than an oversight.
    """
    base = _result(node_g=[1.0, 1.0], node_r=[2.0, 2.0], node_eff=[100.0, 100.0], edge_eff=[50.0])
    import dataclasses

    loud = dataclasses.replace(base, mass_rna_junction=np.array([10_000.0]), n_junctions=1)
    ra = _regions([0, 100], [100, 200])
    ml = [_ml(0, [(0, 0, 200)])]
    quiet_priors = assemble_priors(base, ra, ml)
    loud_priors = assemble_priors(loud, ra, ml)
    np.testing.assert_array_equal(quiet_priors.rna_prior_count, loud_priors.rna_prior_count)
    np.testing.assert_array_equal(quiet_priors.gdna_prior_count, loud_priors.gdna_prior_count)


def test_node_split_between_two_loci():
    # One node [0,100) straddling two adjacent loci ([0,50), [50,100)) → overlap shares 0.5/0.5;
    # each locus gets half of every projected quantity. Single node ⇒ no line ⇒ span = node_eff.
    cal = _result(node_g=[5.0], node_r=[10.0], node_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0, 5.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [2.5, 2.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_intergenic_node_dropped():
    # Node 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _result(node_g=[1.0, 50.0], node_r=[5.0, 99.0], node_eff=[100.0, 100.0])
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.gdna_prior_count, [1.0])
    assert priors.gdna_eff_len[0] > 0.0


def test_a_line_off_an_INTERGENIC_left_flank_is_rekeyed_to_its_right():
    """⛔ A locus's far-LEFT outer line has an intergenic left flank — a node the projection DROPS.

    Keying the line there loses its crossing gDNA silently: the locus's prior is under-counted and its
    eff-length inflated, with no shape error anywhere. The far-RIGHT line is already kept (its left
    flank is the locus's last node), so the re-key restores the symmetry.
    """
    from rigel.calibration.priors import edge_owner_nodes
    from rigel.calibration.signature import BIT_EXON_POS

    # node 0 intergenic, nodes 1-2 exonic ⇒ line 0 is the far-LEFT outer line, line 1 is interior.
    cal = _result(
        node_g=[0.0, 0.0, 0.0],
        node_r=[0.0, 0.0, 0.0],
        node_eff=[100.0, 100.0, 100.0],
        edge_g=[7.0, 3.0],
        edge_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300], signature=[0, BIT_EXON_POS, BIT_EXON_POS])
    np.testing.assert_array_equal(edge_owner_nodes(cal, ra), [1, 1])  # line 0 re-keyed off node 0
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 100, 300)])])  # the locus is nodes 1-2 only
    # ⭐ P1: mass 7+3 = 10 over ΣS = (100+50+50) + 100 = 300, times span_bp = 200 ⇒ 6.6667.
    # ⚠ The discrimination is UNCHANGED in strength: had line 0 stayed keyed to the dropped intergenic
    # node 0, the locus would see mass 3 over ΣS 250 ⇒ 2.4, not 6.667.
    np.testing.assert_allclose(priors.gdna_prior_count, [10.0 / 300.0 * 200.0])  # nothing lost


# --- Laplace shrinkage toward the (effective) span ------------------------------------------------


def test_evidence_free_node_gives_zero_gdna_prior():
    # Acyclic: no observed gDNA ⇒ zero gDNA pseudocount. With G=0 the Laplace-smoothed IPR is
    # (0+1)²/(1/span) = span exactly, so the eff-len is the uniform effective span (single node ⇒
    # span = node_eff = 100) — never a tiny length.
    cal = _result(node_g=[0.0], node_r=[0.0], node_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])
    np.testing.assert_allclose(priors.rna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])  # effective-span fallback


def _six_node_ra():
    return _regions([0, 100, 200, 300, 400, 500], [100, 200, 300, 400, 500, 600])


def _global_bimodal_cal(rna0: float, gdna0: float = 1.0) -> CalibrationResult:
    """6-node calibration for the CONTAINED-EVIDENCE SHRINKAGE under the GLOBAL reference. Five ENRICHED
    background nodes (gDNA density 1.0) fix ρ_ref at the enriched mode; node 0 is DEPLETED (density
    gdna0/100 ≪ ρ_ref) so it contracts (elen ≪ span). ``rna0`` tunes the contained evidence
    C = gdna0 + rna0 that drives the shrinkage weight w = C/(C+1). No crossing mass; node_eff = 100."""
    mg = np.full(6, 100.0)
    mg[0] = gdna0
    mr = np.zeros(6)
    mr[0] = rna0
    return _result(
        node_g=mg,
        node_r=mr,
        node_eff=np.full(6, 100.0),
        edge_eff=np.full(5, 50.0),
        gdna_density_global=0.5,
    )


def test_eff_len_shrinks_toward_span_for_sparse_gdna():
    # Under the GLOBAL reference a DEPLETED locus (density ≪ ρ_ref, set by the enriched background)
    # contracts (elen ≪ span); the contained-evidence shrinkage pulls that contraction back toward span
    # when the unique-mapper evidence is sparse.
    ra = _six_node_ra()
    ml = [_ml(0, [(0, 0, 100)])]  # locus = node 0 (the depleted one)
    eff_sparse = assemble_priors(_global_bimodal_cal(rna0=0.3), ra, ml).gdna_eff_len[0]
    eff_abundant = assemble_priors(_global_bimodal_cal(rna0=100.0), ra, ml).gdna_eff_len[0]
    assert eff_abundant < eff_sparse  # more contained evidence ⇒ more of the earned contraction


def _blind_line_cal(contained_rna: float) -> CalibrationResult:
    """A locus whose gDNA contraction is driven by FIXED crossing mass (no contained gDNA), with a
    tunable amount of CONTAINED RNA — to isolate the contained-evidence shrinkage. node_eff=100 each,
    edge_eff=50 each ⇒ effective span = 300 + 2·50 = 400."""
    return _result(
        node_g=[0.0, 0.0, 0.0],  # no contained gDNA — all signal is crossing
        node_r=[contained_rna, 0.0, 0.0],  # the only contained (unique-mapper) evidence
        node_eff=[100.0, 100.0, 100.0],
        edge_g=[2.0, 3.0],
        edge_eff=[50.0, 50.0],
    )


def test_contained_evidence_shrinkage_reverts_to_span_when_blind():
    # C = 0 (multimapper-blind locus: zero contained mass, only crossing): w = 0 ⇒ eff_len → span (no
    # contraction), the smooth shrinkage's C→0 limit. The prior is honestly uninformative where
    # calibration cannot see. The crossing mass is still counted in the prior pseudocount.
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(_blind_line_cal(0.0), ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_eff_len, [400.0])  # effective span = 300 + 2·50
    # ⭐ P1: lines 2+3 = 5 over ΣS = 300 + 100 = 400, times span_bp = 300 ⇒ 3.75. The point stands —
    # crossing-only mass IS still counted where calibration is contained-blind (it is > 0).
    np.testing.assert_allclose(priors.gdna_prior_count, [3.75])


def test_contained_evidence_shrinkage_is_smooth_not_a_cliff():
    # The shrinkage is SMOOTH in contained evidence (not a hard cliff). Counts 0,1,3,1000 interpolate
    # strictly monotonically.
    ra = _six_node_ra()
    ml = [_ml(0, [(0, 0, 100)])]

    def eff(c):
        return assemble_priors(_global_bimodal_cal(rna0=c), ra, ml).gdna_eff_len[0]

    e0, e1, e3, e_hi = eff(0.0), eff(1.0), eff(3.0), eff(1000.0)
    assert e_hi < e3 < e1 < e0  # smooth + monotone (no cliff)
    assert e0 < 150.0  # even C=0 carries the depleted contained gDNA as evidence


def test_empty_multiloci_returns_empty():
    cal = _result(node_g=[1.0], node_r=[1.0], node_eff=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [])
    assert priors.rna_prior_count.shape == (0,)
    assert priors.gdna_prior_count.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_node_count_mismatch_raises():
    cal = _result(node_g=[1.0, 2.0], node_r=[1.0, 2.0], node_eff=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 node vs calibration's 2
    with pytest.raises(ValueError, match="nodes"):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])])


def test_gdna_eff_len_factor_one_under_uniform_gdna_with_kde_firing():
    # Priors-side factor-1: 6 nodes at UNIFORM gDNA density ⇒ the KDE detector FIRES (≥6 gDNA objects)
    # but is unimodal ⇒ ρ_ref = ρ ⇒ every object min(m/ρ_ref, S) = S ⇒ gdna_eff_len == span EXACTLY.
    # The other gDNA-eff tests run in the <5-object regime that returns None, so this locks the KDE path.
    from rigel.calibration.capture_eff_length import _global_reference_density

    ra = _six_node_ra()
    rho = 2.0
    cal = _uniform_field(np.full(6, 100.0), np.full(5, 50.0), rho)
    assert _global_reference_density(cal.mass_gdna_node, cal.gdna_node_eff_len) == pytest.approx(
        rho
    )
    span = 6 * 100.0 + 5 * 50.0  # 850
    np.testing.assert_allclose(
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 600)])]).gdna_eff_len[0], span, rtol=1e-9
    )
