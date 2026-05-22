"""Bayesian-prior acceptance gate (v3 plan, Phase 2/3).

These tests pin the four invariants the Phase 2 redesign must satisfy:

1. ``gdna_prior_count == sum_loci eta_g(loc)`` from
   :func:`expected_gdna_count_global` — and *only* that.
2. ``enable_gdna`` is computed from per-unit ``is_spliced`` and
   ``gdna_log_liks``, **independent of prior strength** — in particular
    ``gdna_prior_count == 0`` does not force ``enable_gdna == 0``.
3. The global-only prior reads no target-locus payload evidence: the
   helper has no ``payload_arrays`` parameter, and perturbing the
    payload's per-region counts must not change ``gdna_prior_count``.

The high-nRNA toy-mini-genome sentinels now live as active regression tests in
``tests/scenarios/test_nrna_double_counting.py``; this file stays focused on the
global-prior assembly contract.
"""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._kappa import KappaEstimate
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.locus_prior import (
    ExpectedGdnaPriorParts,
    PriorTable,
    assemble_priors,
    enable_gdna_for_multilocus,
    expected_gdna_count_global,
)
from rigel.calibration.regions import RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModel
from rigel.locus import Locus, MultiLocus


# ---------------------------------------------------------------------------
# Fakes (mirroring tests/test_assemble_priors.py for self-contained gate)
# ---------------------------------------------------------------------------

def _kappa_zero() -> KappaEstimate:
    return KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason="")


def _delta_fl(length: int, *, max_size: int = 1024) -> FragmentLengthModel:
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _gdt(
    *,
    rho_ig: float = 0.0,
    rho_in: float = 0.0,
    rho_b: float = 0.0,
    fl_mean: int = 200,
) -> GlobalDensityTable:
    return GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC", rho=rho_ig, n_fragments=10,
            eff_length_bp=10_000.0, n_regions_used=1, kappa=_kappa_zero(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON", rho=rho_in, n_fragments=10,
            eff_length_bp=10_000.0, n_regions_used=1, kappa=_kappa_zero(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=rho_b, n_fragments=10,
            eff_length_bp=10_000.0, n_regions_used=1, kappa=_kappa_zero(),
        ),
        gdna_fl=_delta_fl(fl_mean),
    )


def _fake_index(
    region_rows: list[tuple[str, int, int, int, bool, bool]],
    transcripts: list[tuple[str, int, int]],
):
    refs = sorted({r[0] for r in region_rows} | {t[0] for t in transcripts})
    ref_name_to_id = {name: i for i, name in enumerate(refs)}
    n = len(region_rows)
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": [r[0] for r in region_rows],
            "start": np.array([r[1] for r in region_rows], dtype=np.int64),
            "end": np.array([r[2] for r in region_rows], dtype=np.int64),
            "type": np.array([r[3] for r in region_rows], dtype=np.uint8),
            "strand": np.zeros(n, dtype=np.uint8),
            "tx_pos_bp": np.zeros(n, dtype=np.int64),
            "tx_neg_bp": np.zeros(n, dtype=np.int64),
            "exon_pos_bp": np.zeros(n, dtype=np.int64),
            "exon_neg_bp": np.zeros(n, dtype=np.int64),
            "boundary_flux_left": np.array([r[4] for r in region_rows], dtype=bool),
            "boundary_flux_right": np.array([r[5] for r in region_rows], dtype=bool),
        }
    )
    t_df = pd.DataFrame(
        {
            "ref": pd.Categorical([t[0] for t in transcripts], categories=refs),
            "start": np.array([t[1] for t in transcripts], dtype=np.int64),
            "end": np.array([t[2] for t in transcripts], dtype=np.int64),
        }
    )
    return SimpleNamespace(
        region_df=region_df,
        t_df=t_df,
        ref_name_to_id=ref_name_to_id,
        ref_lengths={name: 1_000_000 for name in refs},
    )


def _make_payload(
    n_regions: int,
    counts_intergenic: list[int] | None = None,
    counts_intron: list[int] | None = None,
    u_left: list[int] | None = None,
    u_right: list[int] | None = None,
) -> CalibrationScanPayload:
    per_region = np.zeros((n_regions, MASK_N_STATES), dtype=np.int64)
    if counts_intergenic is not None:
        per_region[:, 0b100] = np.array(counts_intergenic, dtype=np.int64)
    if counts_intron is not None:
        per_region[:, 0b010] = np.array(counts_intron, dtype=np.int64)
    intron_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    intron_by_orient[:, 2] = per_region[:, 0b010]
    u_left_arr = np.array(u_left if u_left is not None else [0] * n_regions, dtype=np.int64)
    u_right_arr = np.array(u_right if u_right is not None else [0] * n_regions, dtype=np.int64)
    u_left_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    u_right_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    u_left_by_orient[:, 2] = u_left_arr
    u_right_by_orient[:, 2] = u_right_arr
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=per_region,
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=u_left_arr,
        u_right=u_right_arr,
        intron_counts_by_orient=intron_by_orient,
        exon_contained_counts_by_orient=np.zeros((n_regions, 3), dtype=np.int64),
        u_left_by_orient=u_left_by_orient,
        u_right_by_orient=u_right_by_orient,
        n_observed=0, n_excluded_multimap=0, n_excluded_chimera=0,
        n_excluded_artifact=0, n_unobserved=0, n_unannotated_ref=0,
    )


def _make_em(
    locus_t_indices: np.ndarray,
    *,
    n_units: int | None = None,
    is_spliced: np.ndarray | None = None,
    gdna_log_liks: np.ndarray | None = None,
) -> SimpleNamespace:
    if n_units is None:
        n_units = int(locus_t_indices.size + 32)
    return SimpleNamespace(
        locus_t_indices=locus_t_indices,
        is_spliced=(
            is_spliced if is_spliced is not None
            else np.zeros(n_units, dtype=np.uint8)
        ),
        gdna_log_liks=(
            gdna_log_liks if gdna_log_liks is not None
            else np.full(n_units, -1.0, dtype=np.float64)
        ),
    )


def _ml_single(
    multi_locus_id: int,
    transcript_indices: list[int],
    unit_indices: list[int],
    locus: Locus,
) -> MultiLocus:
    return MultiLocus(
        multi_locus_id=multi_locus_id,
        transcript_indices=np.array(transcript_indices, dtype=np.int32),
        unit_indices=np.array(unit_indices, dtype=np.int32),
        gdna_span=int(locus.span),
        loci=(locus,),
    )


# ---------------------------------------------------------------------------
# Acceptance gate
# ---------------------------------------------------------------------------

class TestPhase2GlobalOnlyContract:
    """Pin the four Phase 2 invariants from the v3 plan."""

    # --- (4) The global-only helper has no payload dependency. -----------

    def test_helper_signature_excludes_payload(self):
        """``expected_gdna_count_global`` must not accept payload args.

        This is a *type-level* guarantee that local counts cannot leak
        into the ordinary global-only prior (v3 plan §5 Phase 2).
        """
        sig = inspect.signature(expected_gdna_count_global)
        for forbidden in ("payload", "payload_arrays", "u_left", "u_right"):
            assert forbidden not in sig.parameters, (
                f"expected_gdna_count_global must not accept '{forbidden}'; "
                f"the global-only prior is payload-free by design."
            )

    def test_prior_invariant_to_payload_perturbation(self):
        """Perturbing per-region counts must not change ``gdna_prior_count``.

        The global-only prior is a function of locus geometry and the
        *global* densities only. Local payload counts feed only the
        diagnostic ``estimate_locus_gdna`` path.
        """
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
            ],
            transcripts=[("chr1", 100, 800)],
        )
        locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
        ml = _ml_single(0, [0], [0] * 10, locus)
        em = _make_em(np.zeros(10, dtype=np.int32))
        gdt = _gdt(rho_ig=1e-3, fl_mean=200)

        pt_a = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(n_regions=1, counts_intergenic=[5]),
            global_densities=gdt,
        )
        pt_b = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(
                n_regions=1, counts_intergenic=[500],
                u_left=[1234], u_right=[5678],
            ),
            global_densities=gdt,
        )
        assert pt_a.gdna_prior_count[0] == pytest.approx(pt_b.gdna_prior_count[0])
        # The global-only prior must be strictly positive here
        # (rho_ig > 0 and the locus has intergenic exposure).
        assert pt_a.gdna_prior_count[0] > 0.0

    # --- (1) gdna_prior_count == sum_loci eta_g(loc) ---------------------

    def test_gdna_prior_count_equals_sum_of_helper(self):
        """``gdna_prior_count`` must equal the helper totals summed."""
        # Two-locus multi-locus: each locus contributes its own eta_g.
        loc0 = Locus(ref="chr1", ref_id=0, start=0, end=1000)
        loc1 = Locus(ref="chr1", ref_id=0, start=2000, end=3500)
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
                ("chr1", 2000, 3500, int(RegionType.INTRON), True, True),
            ],
            transcripts=[("chr1", 100, 800), ("chr1", 2100, 3400)],
        )
        ml = MultiLocus(
            multi_locus_id=0,
            transcript_indices=np.array([0, 1], dtype=np.int32),
            unit_indices=np.array([0, 1, 2, 3, 4], dtype=np.int32),
            gdna_span=2500,
            loci=(loc0, loc1),
        )
        em = _make_em(np.zeros(5, dtype=np.int32))
        gdt = _gdt(rho_ig=1e-3, rho_in=2e-3, rho_b=5e-3, fl_mean=200)
        gdna_fl = _delta_fl(200)

        pt = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(n_regions=2),
            global_densities=gdt, gdna_fl=gdna_fl,
        )

        # Recompute the helper independently and compare.
        from rigel.calibration.locus_prior import (
            RegionArrays, RegionIndexPy,
        )
        region_arrays = RegionArrays.from_region_df(
            index.region_df, index.ref_name_to_id,
        )
        region_index = RegionIndexPy(arrays=region_arrays)
        eta_total = sum(
            expected_gdna_count_global(
                locus=loc, region_index=region_index,
                region_arrays=region_arrays, global_densities=gdt,
                gdna_fl=gdna_fl,
            ).total
            for loc in ml.loci
        )
        assert pt.gdna_prior_count[0] == pytest.approx(eta_total)
        assert pt.gdna_prior_count[0] > 0.0  # sanity: non-trivial setup

    # --- (2) enable_gdna independent of prior strength -------------------

    def test_enable_gdna_unspliced_unit_with_zero_prior(self):
        """``gdna_prior_count == 0`` does not force ``enable_gdna == 0``.

        The eligibility flag must be driven solely by per-unit
        ``is_spliced`` / ``gdna_log_liks``.
        """
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
            ],
            transcripts=[("chr1", 100, 800)],
        )
        locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
        ml = _ml_single(0, [0], [0, 1, 2, 3], locus)
        # Default _make_em: every unit unspliced + finite gDNA log-lik.
        em = _make_em(np.zeros(4, dtype=np.int32))
        # Zero global density => eta_g == 0.
        pt = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(n_regions=1),
            global_densities=_gdt(),
        )
        assert pt.gdna_prior_count[0] == pytest.approx(0.0)
        assert pt.enable_gdna[0] == 1, (
            "enable_gdna must be 1 when at least one unit is unspliced "
            "with a finite gDNA log-lik, even if gdna_prior_count == 0."
        )

    def test_enable_gdna_all_spliced_with_positive_prior(self):
        """All-spliced multi-locus is ineligible even with positive prior."""
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
            ],
            transcripts=[("chr1", 100, 800)],
        )
        locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
        n = 6
        ml = _ml_single(0, [0], list(range(n)), locus)
        em = _make_em(
            np.zeros(n, dtype=np.int32),
            n_units=n,
            is_spliced=np.ones(n, dtype=np.uint8),
            gdna_log_liks=np.full(n, -1.0, dtype=np.float64),
        )
        # Positive intergenic density => gdna_prior_count > 0.
        pt = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(n_regions=1),
            global_densities=_gdt(rho_ig=1e-3, fl_mean=200),
        )
        assert pt.gdna_prior_count[0] > 0.0
        assert pt.enable_gdna[0] == 0, (
            "enable_gdna must be 0 when no unit is unspliced with a "
            "finite gDNA log-lik, even if gdna_prior_count > 0."
        )

    def test_enable_gdna_helper_matches_unit_logic(self):
        """``enable_gdna_for_multilocus`` is the canonical decision."""
        n = 5
        # Unit 2 is unspliced with finite gdna log-lik; rest are spliced.
        is_spliced = np.array([1, 1, 0, 1, 1], dtype=np.uint8)
        gll = np.array([-np.inf, -np.inf, -2.0, -np.inf, -np.inf])
        em = _make_em(
            np.zeros(n, dtype=np.int32),
            n_units=n, is_spliced=is_spliced, gdna_log_liks=gll,
        )
        ml = _ml_single(
            0, [0],
            list(range(n)),
            Locus(ref="chr1", ref_id=0, start=0, end=1000),
        )
        assert enable_gdna_for_multilocus(ml, em) is True

        # Now hide unit 2 from the multi-locus' unit_indices.
        ml2 = MultiLocus(
            multi_locus_id=0,
            transcript_indices=np.array([0], dtype=np.int32),
            unit_indices=np.array([0, 1, 3, 4], dtype=np.int32),
            gdna_span=1000,
            loci=(Locus(ref="chr1", ref_id=0, start=0, end=1000),),
        )
        assert enable_gdna_for_multilocus(ml2, em) is False


# ---------------------------------------------------------------------------
# Output-schema gate (v3 plan §6: "output schema/goldens reflect the
# new prior semantics")
# ---------------------------------------------------------------------------

class TestPhase2PriorTableSchema:
    def test_prior_table_phase2_columns_present(self):
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
            ],
            transcripts=[("chr1", 100, 800)],
        )
        locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
        ml = _ml_single(0, [0], [0] * 3, locus)
        em = _make_em(np.zeros(3, dtype=np.int32))
        pt = assemble_priors(
            multi_loci=[ml], em_data=em, index=index,
            payload=_make_payload(n_regions=1),
            global_densities=_gdt(rho_ig=1e-3),
        )
        assert isinstance(pt, PriorTable)
        assert hasattr(pt, "gdna_prior_count")
        assert hasattr(pt, "enable_gdna")
        assert pt.enable_gdna.dtype == np.uint8
        assert pt.gdna_prior_count.dtype == np.float64
        assert not hasattr(pt, "rna_prior_count")
        assert not hasattr(pt, "alpha_rna")
        # And the legacy c_base scalar is gone.
        assert not hasattr(pt, "c_base_value")

    def test_expected_parts_decomposition_sums_to_total(self):
        """``ExpectedGdnaPriorParts.total`` must equal the sum of mass terms."""
        from rigel.calibration.locus_prior import (
            RegionArrays, RegionIndexPy,
        )
        index = _fake_index(
            region_rows=[
                ("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False),
                ("chr1", 1000, 2000, int(RegionType.INTRON), True, True),
                ("chr1", 1500, 1700, int(RegionType.EXON), True, False),
            ],
            transcripts=[("chr1", 100, 800)],
        )
        region_arrays = RegionArrays.from_region_df(
            index.region_df, index.ref_name_to_id,
        )
        region_index = RegionIndexPy(arrays=region_arrays)
        gdt = _gdt(rho_ig=1e-3, rho_in=2e-3, rho_b=5e-3, fl_mean=200)
        parts = expected_gdna_count_global(
            locus=Locus(ref="chr1", ref_id=0, start=0, end=2000),
            region_index=region_index,
            region_arrays=region_arrays,
            global_densities=gdt,
            gdna_fl=_delta_fl(200),
        )
        assert isinstance(parts, ExpectedGdnaPriorParts)
        recomputed = (
            parts.intergenic_contained
            + parts.intron_contained
            + parts.boundary_crossing_expected
            + parts.exon_contained_expected
        )
        assert parts.total == pytest.approx(recomputed)
