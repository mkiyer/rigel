"""assemble_priors — acyclic CalibrationResult → per-locus EM prior.

⭐⭐ **A NODE OWNS THE FRAGMENTS CONTAINED IN IT; AN EDGE OWNS THE FRAGMENTS THAT CROSS IT.** A locus
collects both — its nodes by genomic overlap, its edges by touching those nodes — and no line's mass is
ever folded into a node's total. The rule itself is gated in `test_edge_locus_projection.py`; this file
gates what ``assemble_priors`` builds on top of it.

The object set is a per-node CONTAINED object at effective support ``S_r = E_f[(L_r − w + 1)+] =
gdna_node_eff_len`` plus a per-line CROSSING object at ``S_e = E_f[w − 1] = gdna_edge_eff_len``. The
bedrock invariant these tests pin: under a UNIFORM gDNA field (every object's mass = ρ·S) every
object's ``min(m/ρ_ref, S)`` returns ``S``, so ``gdna_eff_len == span == ΣS`` exactly — an unenriched
library contracts NOTHING (factor 1). Using the genomic ``region_size_bp`` as the divisor instead would
understate short-node density and fabricate a contraction; the dedicated tests below prove the method
uses the EFFECTIVE support, not the genomic length.

⭐ **THE UNIFORM-FIELD FIXTURE IS NOW THREE LINES, AND THAT IS THE EVIDENCE.** It used to carry a
ten-line note explaining that ``gdna_boundary_len`` was ALREADY the halved per-side density length
``E[min(ℓ,L)]/2``, that each face therefore deposited ``ρ·gdna_boundary_len``, and that an earlier
version of the fixture had stored the UN-halved length while depositing half the mass — cancelling
exactly, and hiding a factor of 2 from every assertion in this file for months
A contiguous edge is a 0-bp line with one mass and one support, so a
uniform field is just ``mass = ρ·support`` on both axes and there is no ½ left to get wrong.

⚠ **Every span below is byte-identical to the pre-S5.f value** (640 / 650 / 700 / 400 / 850). The
schema changed; the geometry did not. A number that moved here would mean the port re-derived
something rather than re-keying it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import assemble_priors, contended_edges
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.calibration.signature import BIT_EXON_POS
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
    mass_per_crossing=None,
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
        # ⭐ GEOMETRY, not a split. 1.0 is the identity — a line whose flanks both exceed every
        # fragment length, where one crossing IS one fragment. A test exercising K-inflation overrides it.
        edge_mass_per_crossing=(
            np.ones_like(ez) if mass_per_crossing is None
            else np.asarray(mass_per_crossing, dtype=np.float64)
        ),
        mass_rna_junction=np.zeros(0, dtype=np.float64),
        edge_spliced_mass_per_crossing=np.ones_like(ez),
        junction_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_node_eff_len=node_eff_arr,
        gdna_edge_eff_len=edge_eff_arr,
        rna_node_eff_len=(
            node_eff_arr if rna_node_eff is None else np.asarray(rna_node_eff, dtype=np.float64)
        ),
        rna_edge_eff_len=(
            edge_eff_arr if rna_edge_eff is None else np.asarray(rna_edge_eff, dtype=np.float64)
        ),
        gdna_frac_node=np.zeros_like(ng),
        rna_pos_frac_node=np.zeros_like(ng),
        rna_neg_frac_node=np.zeros_like(ng),
        gdna_frac_edge=ez.copy(),
        rna_pos_frac_edge=ez.copy(),
        rna_neg_frac_edge=ez.copy(),
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
    # ⭐ The prior is a CONSERVED FRAGMENT COUNT read out of the bank: the contained mass ρ·Σnode_eff
    # (8.0) plus the crossing mass ρ·Σedge_eff (4.8) rescaled by q, which this fixture sets to the
    # identity 1.0 — flanks exceeding every fragment length, where one crossing IS one fragment.
    # ⚠ It replaces `ρ · span_bp` = 8.0, the retired density rule: that reached fragment units by
    # dividing the mass by its own opportunity and re-integrating, and dropped the 4.8 of crossing
    # fragments entirely, because a 0-bp line contributes no genomic span to integrate over.
    np.testing.assert_allclose(priors.gdna_prior_count, [rho * (sum(node_eff) + sum(edge_eff))], rtol=1e-9)
    assert not np.isclose(priors.gdna_prior_count[0], rho * 400.0)  # ⛔ not the retired ρ·span_bp


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


def test_every_OBJECT_has_the_same_density_under_a_uniform_field():
    """The precondition for the per-object ``min()`` factor-1 identity, asserted on the objects
    themselves rather than on a folded total.

    ⛔ It used to read through ``_component_node_arrays``, which summed each line's mass into a flank
    node before dividing — so it could only ever check the FOLD's density, never a line's own. Nodes and
    edges are peers now, so each axis is checked on its own axis.
    """
    node_eff = np.array([120.0, 200.0, 80.0])
    edge_eff = np.array([120.0, 120.0])
    rho = 0.02
    cal = _uniform_field(node_eff, edge_eff, rho)
    np.testing.assert_allclose(cal.mass_gdna_node / cal.gdna_node_eff_len, rho, rtol=1e-9)
    np.testing.assert_allclose(cal.mass_gdna_edge / cal.gdna_edge_eff_len, rho, rtol=1e-9)


# --- mass / projection (independent of the support choice) ----------------------------------------


def test_single_locus_projects_both_components():
    # ⭐ The priors are CONSERVED FRAGMENT COUNTS. No crossing mass here, so both are the contained
    # mass alone — one deposit per contained fragment, nothing to convert:
    #   gDNA: Σm = 4.5      (the retired ρ·span rule gave 4.5/750 · 450 = 2.7)
    #   RNA : Σm = 12.0     (the retired rule gave 7.2)
    # ⚠ The g:r RATIO is 0.375 either way, because a common divisor cancels from a ratio. That is
    # exactly why the ratio is NOT what discriminates the two rules — the totals are.
    cal = _result(
        node_g=[1.0, 2.0, 1.5],
        node_r=[3.0, 4.0, 5.0],
        node_eff=[100.0, 200.0, 150.0],
        edge_eff=[150.0, 150.0],
    )
    ra = _regions([0, 100, 300], [100, 300, 450])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 450)])])
    np.testing.assert_allclose(priors.rna_prior_count, [12.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
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
    # ⭐ OBJECT conservation is what this test is named for: the locus covers every node, so it collects
    # every node AND every line, and the prior is their total — nothing dropped, nothing double-counted.
    # ⛔⛔ AND HERE THE PRIOR EQUALS THE RAW SUM — 11.0 — WHICH IS NOT EVIDENCE THAT IT *IS* A RAW SUM.
    # This fixture's q is the identity 1.0, so incidence and fragment coincide by construction and this
    # test CANNOT tell the two rules apart. The discrimination lives in the q ≠ 1 test below and in
    # `test_prior_units.py`; asserting 11.0 here would otherwise read as a ruling that it is a raw sum.
    np.testing.assert_allclose(priors.gdna_prior_count, [11.0])
    np.testing.assert_allclose(cal.edge_mass_per_crossing, 1.0)  # ...the reason it coincides, pinned
    np.testing.assert_allclose(
        priors.gdna_prior_count.sum(), cal.mass_gdna_node.sum() + cal.mass_gdna_edge.sum()
    )
    assert contended_edges(ra, [_ml(0, [(0, 0, 300)])], 1).size == 0  # nothing double-claimed


def test_the_crossing_mass_is_rescaled_by_the_conserved_share():
    """⭐⭐ **THE ONE TEST IN THIS FILE THAT SEPARATES A CONSERVED COUNT FROM A RAW INCIDENCE SUM.**

    Every other fixture here leaves ``edge_mass_per_crossing`` at the identity 1.0, where one crossing
    IS one fragment and the two rules coincide — so they all pass under either. This one sets ``q`` to
    ``[0.5, 0.25]``: a fragment crossing line 0 deposited on 2 objects on average, line 1 on 4.

        gDNA = 3 (contained, one deposit each) + 4·0.5 + 8·0.25 = 7.0     raw sum would be 15.0
        RNA  = 6                               + 4·0.5 + 4·0.25 = 9.0     raw sum would be 14.0

    ⛔ The CONTAINED term must NOT be rescaled — a contained fragment touches exactly one node and is
    already a count. Rescaling it too would give 3·? and is the other wrong answer this pins out.
    """
    cal = _result(
        node_g=[1.0, 1.0, 1.0],
        node_r=[2.0, 2.0, 2.0],
        node_eff=[100.0, 100.0, 100.0],
        edge_g=[4.0, 8.0],
        edge_r=[4.0, 4.0],
        edge_eff=[50.0, 50.0],
        mass_per_crossing=[0.5, 0.25],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])])
    np.testing.assert_allclose(priors.gdna_prior_count, [7.0])
    np.testing.assert_allclose(priors.rna_prior_count, [9.0])
    assert not np.isclose(priors.gdna_prior_count[0], 15.0)  # ⛔ not the raw incidence sum
    assert not np.isclose(priors.rna_prior_count[0], 14.0)


def test_spliced_mass_withheld_from_rna_prior():
    # A spliced fragment has no gDNA candidate in the EM (gDNA does not splice) → it is guaranteed-RNA
    # and assigned directly, so it must NOT load rna_prior_count. Node RNA [3,4,5] (Σ=12) plus line RNA
    # [4,4] of which [1,3] is spliced ⇒ RNA mass = 12 + (4−1) + (4−3) = 16 (NOT 20).
    # ⭐ q is the identity here, so the conserved fragment count is that 16 unchanged; gDNA is its
    # contained 4.5 (no crossing mass). The retired ρ·span rule read 6.4 and 1.8.
    # The WITHHOLDING is what this test pins, and it survives the units change intact: without it the
    # RNA mass would be 20 and the prior 20.0.
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
    np.testing.assert_allclose(priors.rna_prior_count, [16.0])
    np.testing.assert_allclose(priors.gdna_prior_count, [4.5])
    assert priors.rna_prior_count[0] < 20.0  # the spliced mass really is withheld


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

    loud = dataclasses.replace(
        base,
        mass_rna_junction=np.array([10_000.0]),
        junction_mass_per_crossing=np.ones(1),
        n_junctions=1,
    )
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


def test_a_locus_keeps_the_outer_line_against_its_INTERGENIC_flank():
    """⭐ A locus's far-LEFT outer line has an intergenic left flank — a node the projection DROPS.

    A fragment crossing that line overlaps the locus, so it is one of its EM candidates and its mass
    must load the locus's prior. ⛔ The old assembler folded a line's mass into ONE flank node and so
    lost this line into the dropped intergenic flank — it needed an explicit intergenic RE-KEY to get it
    back. There is nothing to re-key now: the line touches node 1, so it is the locus's line.

    ⚠ The discrimination is UNCHANGED in strength: were the outer line dropped, the locus would see 3
    rather than 10 — a 70 % under-count with no shape error anywhere.
    """
    # node 0 intergenic, nodes 1-2 exonic ⇒ line 0 is the far-LEFT outer line, line 1 is interior.
    cal = _result(
        node_g=[0.0, 0.0, 0.0],
        node_r=[0.0, 0.0, 0.0],
        node_eff=[100.0, 100.0, 100.0],
        edge_g=[7.0, 3.0],
        edge_eff=[50.0, 50.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300], signature=[0, BIT_EXON_POS, BIT_EXON_POS])
    ml = [_ml(0, [(0, 100, 300)])]  # the locus is nodes 1-2 only
    priors = assemble_priors(cal, ra, ml)
    np.testing.assert_allclose(priors.gdna_prior_count, [10.0])  # 7 + 3, nothing lost
    assert contended_edges(ra, ml, 1).size == 0


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
    # ⭐ lines 2+3 = 5 crossing fragments (q is the identity). The point stands — crossing-only mass IS
    # still counted where calibration is contained-blind (it is > 0). The retired rule read 3.75.
    np.testing.assert_allclose(priors.gdna_prior_count, [5.0])


def _stray_on_a_dead_line_cal(stray: float) -> CalibrationResult:
    """7 nodes at a UNIFORM gDNA density of 1.0 (so ρ_ref = 1.0 and the KDE fires), 6 lines.

    Node 0 is intergenic, so ``edge_owner_nodes`` re-keys line 0 onto node 1 — which therefore owns
    BOTH line 0 (support **0**, carrying ``stray``) and line 1 (support 50, DEPLETED at mass 5). The
    depletion is what makes ``min(pooled/ρ_ref, seam_len)`` bind on the mass side, so stray mass there
    can actually move the answer.
    """
    return _result(
        node_g=[100.0] * 7,
        node_r=[0.0] * 7,
        node_eff=[100.0] * 7,
        edge_g=[stray, 5.0, 50.0, 50.0, 50.0, 50.0],
        edge_eff=[0.0, 50.0, 50.0, 50.0, 50.0, 50.0],
        gdna_density_global=1.0,
    )


def test_stray_mass_on_a_zero_opportunity_line_is_dropped_from_the_eff_len():
    """⭐ **THE P1e PERTURBATION, RELOCATED TO WHERE THE DIVISOR STILL LIVES.** Its other half is
    `test_prior_units.test_mass_on_a_zero_opportunity_object_STILL_COUNTS_because_a_count_has_no_divisor`,
    which asserts the opposite for the PRIOR — deliberately, because a count has nothing to divide by.
    The eff-length still divides by ρ_ref and still pools, so here the drop is load-bearing.

    ``mass > 0`` with ``support == 0`` is an ordinary configuration: ``contained_eff_length`` is exactly
    0 wherever an object is shorter than that component's shortest fragment (21.7 % of chr22 nodes for
    RNA, 18.7 % for gDNA), and the solver can still put mass there because ``f_g`` is an inference.

    ⛔ Two wrong answers this pins out, both measured by injecting them:

    * ``mass / max(support, 1e-9)`` — a density of ~1e9, which is how a "no data" default of 100 % gDNA
      once seeded false gDNA into neighbouring exons;
    * mass kept in the numerator with its support omitted from the denominator — ``ρ`` inflated with no
      exposure to pay for it. **Both sides of a pooled rate, or neither.**

    Either one moves the eff-length by **+19.97 bp** at ``stray = 20`` and **+44.93 bp** at
    ``stray = 5000``, where it pins against the 850 bp seam ceiling. With the drop it does not move at
    all, and this test sweeps 250× of stray mass to say so.
    """
    ra = _regions(list(range(0, 700, 100)), list(range(100, 800, 100)),
                  signature=[0] + [BIT_EXON_POS] * 6)
    ml = [_ml(0, [(0, 100, 700)])]  # the locus is nodes 1-6; node 0 is intergenic and dropped
    quiet = assemble_priors(_stray_on_a_dead_line_cal(0.0), ra, ml).gdna_eff_len[0]
    for stray in (20.0, 5000.0):
        loud = assemble_priors(_stray_on_a_dead_line_cal(stray), ra, ml)
        np.testing.assert_allclose(loud.gdna_eff_len, [quiet], rtol=1e-12)
        # ⛔ non-vacuity: the eff-len is genuinely contracted below the seam ceiling the undropped
        # mass would push it to, so "unchanged" is a real constraint and not both arms at the clamp.
        assert quiet < 850.0
        # ...and the stray mass is NOT silently discarded everywhere — the prior still counts it.
        assert loud.gdna_prior_count[0] == pytest.approx(
            assemble_priors(_stray_on_a_dead_line_cal(0.0), ra, ml).gdna_prior_count[0] + stray
        )


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


def test_the_contraction_is_applied_PER_OBJECT_not_over_a_folded_total():
    """⭐⭐ **THE GATE THE ``min()`` DOCSTRING CLAIMED AND NOTHING ENFORCED**, found by perturbation.

    ``elen`` contracts each object separately — ``Σ min(m_n/ρ_ref, S_n)`` — and the docstring has long
    said that folding them into one ``min()`` over the summed mass would UNDER-contract *"a captured
    exon whose line runs into a depleted intron"*. Nothing tested it.

    ⛔ **A single-density fixture cannot**: under a uniform field every object sits exactly at ``ρ_ref``,
    both forms return ``span``, and the locus-level clamp to ``span`` hides any excess anyway. The
    discriminating shape needs one object ABOVE ``ρ_ref`` and one BELOW, so the enriched object's excess
    would compensate the depleted one's deficit under a fold and cancel — while per object the excess is
    clipped and the deficit survives.

    6 nodes at density 1.0 fix ρ_ref; the locus is node 0 (density 0.1 — depleted) plus its two lines,
    one of which is ENRICHED at 5×. Per object: the enriched line clips to its support. Folded: the
    line's 5× mass pays for the node's shortfall and the locus reads as unenriched.
    """
    ra = _six_node_ra()
    mg = np.full(6, 100.0)
    mg[0] = 10.0  # the locus's node is DEPLETED (density 0.1 against ρ_ref = 1.0)
    eg = np.full(5, 50.0)
    eg[0] = 250.0  # ...and its right line is ENRICHED (density 5.0)
    cal = _result(
        node_g=mg,
        node_r=np.zeros(6),
        node_eff=np.full(6, 100.0),
        edge_g=eg,
        edge_eff=np.full(5, 50.0),
        gdna_density_global=1.0,
    )
    ml = [_ml(0, [(0, 0, 100)])]
    eff = assemble_priors(cal, ra, ml).gdna_eff_len[0]
    # span = node 100 + its ONE line (node 0 owns only its right line) = 150
    # per object: min(10/1, 100) + min(250/1, 50) = 10 + 50 = 60, shrunk toward span by w = C/(C+1)
    span, per_object = 150.0, 60.0
    contained_ev = 10.0
    w = contained_ev / (contained_ev + 1.0)
    np.testing.assert_allclose(eff, w * per_object + (1.0 - w) * span, rtol=1e-9)
    # ⛔ the folded form would read min((10+250)/1, 150) = 150 ⇒ eff == span ⇒ NO contraction at all
    assert eff < span - 1.0, "the enriched line paid for the depleted node — the fold is back"
