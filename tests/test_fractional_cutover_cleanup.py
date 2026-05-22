import numpy as np

from rigel.calibration._exposure import fractional_boundary_side_exposure
from rigel.calibration.fractional_evidence import FractionalEvidenceView
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import N_CHANNELS, N_FL_POOLS, N_SIGNATURES
from rigel.frag_length_model import FragmentLengthModel


def _scan_dict(*, n_observed: int, n_unannotated_ref: int = 0) -> dict[str, object]:
    region_mass = float(n_observed - n_unannotated_ref)
    region_counts = np.zeros((1, N_CHANNELS), dtype=np.float32)
    region_counts[0, 0] = region_mass

    channel_mass = np.zeros(N_CHANNELS, dtype=np.float64)
    channel_mass[0] = region_mass
    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    signature_mass[0] = region_mass

    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    fl_pool_mass[0, 10] = region_mass
    fl_pool_total = fl_pool_mass.sum(axis=1)

    return {
        "n_regions": 1,
        "region_counts": region_counts,
        "channel_mass": channel_mass,
        "signature_mass": signature_mass,
        "fl_pool_mass": fl_pool_mass,
        "fl_pool_total": fl_pool_total,
        "n_observed": n_observed,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_excluded_strand_ambig": 0,
        "n_unobserved": 0,
        "n_unannotated_ref": n_unannotated_ref,
        "n_fl_unavailable": 0,
        "resolver_splicing_anchor_tolerance": 3,
    }


def test_scan_payload_fl_bound_does_not_subtract_strand_ambig() -> None:
    scan_dict = _scan_dict(n_observed=10)
    scan_dict["n_excluded_strand_ambig"] = 7

    payload = CalibrationScanPayload.from_scan_dict(scan_dict, n_total=17)

    assert payload.fl_pool_total.sum() == 10.0


def test_scan_payload_region_mass_accounts_for_unannotated_ref_subset() -> None:
    payload = CalibrationScanPayload.from_scan_dict(
        _scan_dict(n_observed=10, n_unannotated_ref=2),
        n_total=10,
    )

    assert payload.region_counts.sum(dtype=np.float64) == 8.0


def test_fractional_evidence_view_masks_work_with_slots() -> None:
    payload = CalibrationScanPayload.from_scan_dict(_scan_dict(n_observed=1), n_total=1)

    view = FractionalEvidenceView.from_payload(payload, np.array([0], dtype=np.uint8))

    assert view.mask_intergenic.tolist() == [True]
    assert view.mask_intron_only.tolist() == [False]
    assert view.mask_exon_any.tolist() == [False]


def test_fractional_boundary_side_exposure_matches_exact_min_formula() -> None:
    counts = np.zeros(151, dtype=np.float64)
    counts[[50, 100, 150]] = [3.0, 5.0, 2.0]
    fl = FragmentLengthModel.from_counts(counts, max_size=150)
    lengths = np.array([0, 10, 80, 1000], dtype=np.int64)

    ell = np.arange(fl.pmf.size, dtype=np.float64)
    expected = np.array(
        [np.sum(fl.pmf[1:] * np.minimum((ell[1:] - 1.0) / 2.0, span / 2.0)) for span in lengths]
    )

    observed = fractional_boundary_side_exposure(lengths, fl)

    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)
