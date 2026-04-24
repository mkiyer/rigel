"""Tests for capture-class annotation and composite-Poisson calibration EM.

Phase 3a: per-bp ``(E_on, E_off)`` attribution and composite-Poisson
unmixing (``λ_G_on · E_on + λ_G_off · E_off``) via nested inner EM.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration import annotate_capture_class, calibrate_gdna, run_em
from rigel.calibration._stats import compute_region_stats


# ---------------------------------------------------------------------------
# Helpers (mirror tests/test_calibration.py shape)
# ---------------------------------------------------------------------------


def _make_region_df(n, *, ref="chr1", start=None, end=None, length=1000,
                    tx_pos=None, tx_neg=None,
                    exon_pos=None, exon_neg=None,
                    mappable_effective_length=None):
    if start is None:
        start = np.arange(n, dtype=np.int64) * length
    if end is None:
        end = start + length
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    lengths = end - start
    d = {
        "region_id": np.arange(n, dtype=np.int32),
        "ref": np.asarray([ref] * n) if isinstance(ref, str) else ref,
        "start": start,
        "end": end,
        "length": lengths.astype(np.int32),
        "tx_pos": tx_pos if tx_pos is not None else np.ones(n, dtype=bool),
        "tx_neg": tx_neg if tx_neg is not None else np.zeros(n, dtype=bool),
    }
    if exon_pos is not None:
        d["exon_pos"] = np.asarray(exon_pos, dtype=bool)
    else:
        d["exon_pos"] = np.ones(n, dtype=bool)
    if exon_neg is not None:
        d["exon_neg"] = np.asarray(exon_neg, dtype=bool)
    else:
        d["exon_neg"] = np.zeros(n, dtype=bool)
    if mappable_effective_length is None:
        mappable_effective_length = lengths.astype(np.float32)
    d["mappable_effective_length"] = np.asarray(
        mappable_effective_length, dtype=np.float32
    )
    return pd.DataFrame(d)


def _make_region_counts(n_unspliced_pos, n_unspliced_neg,
                        n_spliced_pos=None, n_spliced_neg=None):
    n = len(n_unspliced_pos)
    if n_spliced_pos is None:
        n_spliced_pos = np.zeros(n, dtype=np.float32)
    if n_spliced_neg is None:
        n_spliced_neg = np.zeros(n, dtype=np.float32)
    return pd.DataFrame({
        "region_id": np.arange(n, dtype=np.int32),
        "n_unspliced_pos": np.asarray(n_unspliced_pos, dtype=np.float32),
        "n_unspliced_neg": np.asarray(n_unspliced_neg, dtype=np.float32),
        "n_spliced_pos": np.asarray(n_spliced_pos, dtype=np.float32),
        "n_spliced_neg": np.asarray(n_spliced_neg, dtype=np.float32),
    })


def _write_bed(tmp_path, intervals, name="targets.bed"):
    path = tmp_path / name
    with path.open("w") as fh:
        for ref, start, end in intervals:
            fh.write(f"{ref}\t{start}\t{end}\n")
    return path


# ---------------------------------------------------------------------------
# annotate_capture_class — returns (E_on, E_off) float arrays
# ---------------------------------------------------------------------------


def test_annotate_returns_per_bp_attribution(tmp_path):
    """Full overlap → E_off=0; disjoint → E_on=0; sum == mappable_bp."""
    # Widely-spaced regions keep the BED coverage well under the 50% cap.
    n = 200
    region_df = _make_region_df(
        n, start=np.arange(n, dtype=np.int64) * 10_000,
        end=np.arange(n, dtype=np.int64) * 10_000 + 1000,
        mappable_effective_length=np.full(n, 1000, dtype=np.float32),
    )
    # 80 BED intervals that fully cover the first 80 regions.
    intervals = [(f"chr1", i * 10_000, i * 10_000 + 1000) for i in range(80)]
    # Add filler intervals on chr3 to pass the 100-interval floor.
    intervals += [(f"chr3", i * 10_000, i * 10_000 + 500) for i in range(100)]
    bed_path = _write_bed(tmp_path, intervals)

    E_on, E_off = annotate_capture_class(region_df, bed_path, pad=0, genome_size=3_000_000_000)
    assert E_on.dtype == np.float64
    assert E_off.dtype == np.float64
    assert E_on.shape == (n,)
    # Full probe coverage → E_on == mappable, E_off == 0 for covered regions.
    mappable = region_df["mappable_effective_length"].to_numpy(dtype=np.float64)
    np.testing.assert_allclose(E_on[:80] + E_off[:80], mappable[:80])
    np.testing.assert_allclose(E_on[:80], mappable[:80])
    np.testing.assert_allclose(E_off[:80], 0.0)
    # Far-away regions (chr1 indices >= 80) have no probe overlap → E_on = 0.
    assert (E_on[80:] == 0.0).all()
    np.testing.assert_allclose(E_off[80:], mappable[80:])


def test_annotate_partial_overlap_proportional(tmp_path):
    """Half-probed region should split mappable_bp 50/50."""
    n = 200
    region_df = _make_region_df(
        n, start=np.arange(n, dtype=np.int64) * 10_000,
        end=np.arange(n, dtype=np.int64) * 10_000 + 1000,
        mappable_effective_length=np.full(n, 1000, dtype=np.float32),
    )
    # First region: probe covers half of it [0, 500).  Rest: disjoint chr3.
    intervals = [("chr1", 0, 500)] + [
        ("chr3", i * 2000, i * 2000 + 200) for i in range(150)
    ]
    bed_path = _write_bed(tmp_path, intervals)
    E_on, E_off = annotate_capture_class(region_df, bed_path, pad=0, genome_size=3_000_000_000)
    mappable_0 = float(region_df["mappable_effective_length"].iloc[0])
    # Region 0 is 1000 bp long with 500 bp probe overlap → 50/50 split.
    assert E_on[0] == pytest.approx(0.5 * mappable_0, rel=1e-6)
    assert E_off[0] == pytest.approx(0.5 * mappable_0, rel=1e-6)


def test_annotate_sum_equals_mappable(tmp_path):
    """Invariant: E_on + E_off = mappable_bp for every region."""
    rng = np.random.default_rng(1)
    n = 300
    region_df = _make_region_df(
        n, start=np.arange(n, dtype=np.int64) * 10_000,
        end=np.arange(n, dtype=np.int64) * 10_000 + 1000,
        mappable_effective_length=np.full(n, 1000, dtype=np.float32),
    )
    # Random BED intervals on chr3 (disjoint from regions) → all E_on = 0
    # but the invariant must still hold.  Add a few probe-overlapping
    # intervals on chr1 for coverage.
    intervals = [
        ("chr1", int(s) * 10_000, int(s) * 10_000 + 200)
        for s in rng.integers(0, 50, size=40)
    ] + [
        ("chr3", int(s), int(s) + 200)
        for s in rng.integers(0, 300_000, size=200)
    ]
    bed_path = _write_bed(tmp_path, intervals)
    E_on, E_off = annotate_capture_class(region_df, bed_path, pad=0, genome_size=3_000_000_000)
    mappable = region_df["mappable_effective_length"].to_numpy(dtype=np.float64)
    np.testing.assert_allclose(E_on + E_off, mappable, atol=1e-9)
    assert (E_on >= 0).all()
    assert (E_off >= 0).all()


def test_annotate_padding_extends_on_target(tmp_path):
    """With pad=600 a probe 100 bp past the region reaches inside."""
    # Region at [1000, 2000); probe at [2100, 3100).  No pad: no overlap.
    # pad=600 extends probe to [1500, 3700) → overlaps 500 bp of region.
    n_far = 500
    starts = np.arange(n_far, dtype=np.int64) * 100_000
    ends = starts + 50_000
    filler = pd.DataFrame({
        "region_id": np.arange(n_far, dtype=np.int32),
        "ref": ["chr2"] * n_far,
        "start": starts,
        "end": ends,
        "length": (ends - starts).astype(np.int32),
        "tx_pos": np.ones(n_far, dtype=bool),
        "tx_neg": np.zeros(n_far, dtype=bool),
        "exon_pos": np.ones(n_far, dtype=bool),
        "exon_neg": np.zeros(n_far, dtype=bool),
        "mappable_effective_length": (ends - starts).astype(np.float32),
    })
    test_region = _make_region_df(1, ref="chr1", start=np.array([1000]),
                                  end=np.array([2000]))
    test_region["region_id"] = np.array([n_far], dtype=np.int32)
    region_df = pd.concat([filler, test_region], ignore_index=True)

    intervals = [("chr1", 2100, 3100)] + [
        ("chr3", i * 10_000, i * 10_000 + 500) for i in range(149)
    ]
    bed_path = _write_bed(tmp_path, intervals)

    E_on_nopad, _ = annotate_capture_class(region_df, bed_path, pad=0, genome_size=3_000_000_000)
    assert E_on_nopad[-1] == 0.0

    E_on_padded, _ = annotate_capture_class(region_df, bed_path, pad=600, genome_size=3_000_000_000)
    assert E_on_padded[-1] > 0.0


def test_annotate_rejects_trivial_bed(tmp_path):
    region_df = _make_region_df(100, length=1000)
    bed_path = _write_bed(tmp_path, [
        ("chr1", 0, 1000), ("chr1", 2000, 3000), ("chr1", 5000, 6000),
    ])
    with pytest.raises(ValueError, match="does not look like a capture panel"):
        annotate_capture_class(region_df, bed_path)


def test_annotate_warns_on_ref_mismatch(tmp_path, caplog):
    region_df = _make_region_df(500, ref="chr1", length=10_000)
    intervals = [("chrX", i * 1000, i * 1000 + 500) for i in range(150)]
    bed_path = _write_bed(tmp_path, intervals)

    import logging
    with caplog.at_level(logging.WARNING, logger="rigel.calibration._annotate"):
        E_on, E_off = annotate_capture_class(region_df, bed_path, genome_size=3_000_000_000)
    # No matching chromosomes → everything off-target (E_on ≡ 0).
    assert (E_on == 0.0).all()
    assert any("no matching intervals" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# End-to-end: composite-Poisson EM vs single-class
# ---------------------------------------------------------------------------


def test_em_no_bed_bit_identical_to_zero_e_on():
    """No-BED path (``capture_e=None``) must match ``E_on≡0``, ``E_off≡E``."""
    rng = np.random.default_rng(0)
    n = 400
    lengths = rng.integers(500, 3000, size=n)
    region_df = _make_region_df(
        n,
        start=np.cumsum(np.insert(lengths[:-1], 0, 0)).astype(np.int64),
        end=np.cumsum(lengths).astype(np.int64),
        mappable_effective_length=lengths.astype(np.float32),
    )
    n_u_pos = rng.poisson(lam=np.where(rng.random(n) < 0.8, 50, 5)).astype(np.float32)
    n_u_neg = rng.poisson(lam=2, size=n).astype(np.float32)
    n_s_pos = (rng.random(n) < 0.3).astype(np.float32) * 10
    counts = _make_region_counts(n_u_pos, n_u_neg, n_s_pos)
    stats = compute_region_stats(counts, region_df)

    fit_no = run_em(
        stats, fl_table=None, strand_specificity=0.9,
        mean_frag_len=250.0, max_iterations=20,
    )
    # Explicit composite with E_on≡0 should reproduce the single-rate path.
    E_on_zero = np.zeros(n, dtype=np.float64)
    E_off_full = stats["mappable_bp"].astype(np.float64)
    fit_comp = run_em(
        stats, fl_table=None, strand_specificity=0.9,
        mean_frag_len=250.0, max_iterations=20,
        capture_e=(E_on_zero, E_off_full),
    )
    assert fit_no.capture_class_mode is False
    assert fit_comp.capture_class_mode is True
    # λ_G (legacy effective rate) must match — σ_on·ΣE_on = 0 so the
    # coverage-weighted average degenerates to lam_G_off exactly.
    assert fit_comp.lam_G == pytest.approx(fit_no.lam_G, rel=1e-10)
    assert fit_comp.mu_R == pytest.approx(fit_no.mu_R, rel=1e-10)
    assert fit_comp.sigma_R == pytest.approx(fit_no.sigma_R, rel=1e-10)
    np.testing.assert_allclose(fit_comp.gamma, fit_no.gamma, atol=1e-12)
    # In the E_on≡0 branch lam_G_off equals the single-rate λ_G and
    # the lam_G_on rate is inert (its M-step sees zero weight).
    assert fit_comp.lam_G_off == pytest.approx(fit_no.lam_G, rel=1e-10)


def test_em_composite_recovers_distinct_lambdas():
    """Composite-Poisson unmixing recovers λ_on and λ_off from per-bp data."""
    rng = np.random.default_rng(42)
    n = 800
    lengths = np.full(n, 1000, dtype=np.int64)
    starts = np.cumsum(np.insert(lengths[:-1], 0, 0))
    ends = starts + lengths
    region_df = _make_region_df(
        n, start=starts, end=ends,
        mappable_effective_length=lengths.astype(np.float32),
    )
    # Per-region E_on fraction drawn uniformly in [0, 1]: mix of
    # fully-on, fully-off, and partial regions — the hardest case
    # for a binary-partition model.
    mappable = lengths.astype(np.float64)
    frac_on = rng.uniform(0.0, 1.0, size=n)
    E_on = mappable * frac_on
    E_off = mappable - E_on

    lam_on_true = 0.05
    lam_off_true = 0.001
    rate = lam_on_true * E_on + lam_off_true * E_off
    n_u = rng.poisson(lam=rate).astype(np.float32)
    n_u_pos = np.round(n_u * 0.5).astype(np.float32)
    n_u_neg = (n_u - n_u_pos).astype(np.float32)
    counts = _make_region_counts(n_u_pos, n_u_neg)
    stats = compute_region_stats(counts, region_df)

    fit = run_em(
        stats, fl_table=None, strand_specificity=0.5,
        mean_frag_len=250.0, max_iterations=50,
        capture_e=(E_on, E_off),
    )
    assert fit.capture_class_mode
    # Pure-gDNA regime: all γ → 1 so λ̂ → true rate per channel.  The
    # low-rate off-channel picks up some Poisson count variance and so
    # tolerates a wider relative band than the high-rate on-channel.
    assert fit.lam_G_on == pytest.approx(lam_on_true, rel=0.2)
    assert fit.lam_G_off == pytest.approx(lam_off_true, rel=0.6)
    assert fit.lam_G_on > 10.0 * fit.lam_G_off


def test_em_composite_zero_E_on_region_no_nan():
    """Region with E_on=0 must not produce NaN under composite M-step."""
    rng = np.random.default_rng(3)
    n = 200
    lengths = np.full(n, 1000, dtype=np.int64)
    starts = np.cumsum(np.insert(lengths[:-1], 0, 0))
    ends = starts + lengths
    region_df = _make_region_df(
        n, start=starts, end=ends,
        mappable_effective_length=lengths.astype(np.float32),
    )
    # Half the regions have zero E_on (fully off-target).
    mappable = lengths.astype(np.float64)
    E_on = mappable.copy()
    E_off = np.zeros_like(mappable)
    E_on[::2] = 0.0
    E_off[::2] = mappable[::2]
    # Random gDNA counts.
    n_u_pos = rng.poisson(lam=3, size=n).astype(np.float32)
    n_u_neg = rng.poisson(lam=3, size=n).astype(np.float32)
    counts = _make_region_counts(n_u_pos, n_u_neg)
    stats = compute_region_stats(counts, region_df)
    fit = run_em(
        stats, fl_table=None, strand_specificity=0.5,
        mean_frag_len=250.0, max_iterations=20,
        capture_e=(E_on, E_off),
    )
    assert np.isfinite(fit.lam_G_on)
    assert np.isfinite(fit.lam_G_off)
    assert np.all(np.isfinite(fit.gamma))


def test_calibrate_gdna_plumbs_capture_e():
    rng = np.random.default_rng(7)
    n = 300
    lengths = np.full(n, 1000, dtype=np.int64)
    starts = np.cumsum(np.insert(lengths[:-1], 0, 0))
    ends = starts + lengths
    region_df = _make_region_df(
        n, start=starts, end=ends,
        mappable_effective_length=lengths.astype(np.float32),
    )
    n_u_pos = rng.poisson(lam=5, size=n).astype(np.float32)
    n_u_neg = rng.poisson(lam=5, size=n).astype(np.float32)
    counts = _make_region_counts(n_u_pos, n_u_neg)
    fl = pd.DataFrame({"region_id": [], "frag_len": []}).astype(
        {"region_id": np.int64, "frag_len": np.int64}
    )
    mappable = lengths.astype(np.float64)
    E_on = np.where(np.arange(n) % 3 == 0, mappable, 0.0)
    E_off = mappable - E_on
    result = calibrate_gdna(
        counts, fl, region_df, strand_specificity=0.5,
        mean_frag_len=250.0, capture_e=(E_on, E_off),
    )
    assert result.capture_class_mode
    assert result.lam_G_on is not None
    assert result.lam_G_off is not None
    summary = result.to_summary_dict()
    assert summary["capture_class_mode"] is True
    assert "lambda_gdna_on" in summary
    assert "lambda_gdna_off" in summary
    assert "on_target_enrichment_ratio" in summary
