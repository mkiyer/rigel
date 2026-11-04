"""The oracle cache and the per-object scorer, shared by `calibration_vs_oracle.py` (the scorer) and
`calibration_oracle.py --build` (the cache). A helper module, so it has no row on the shelf.

T is the production accumulator run on the BAM split by origin, asserted to sum to the full payload in the drained
frame. `load_or_build_oracle` is the one place the per-origin caches are built.

Gate: ``tests/calibration/test_oracle_arms.py``.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))

from _oracle import (  # noqa: E402
    ORIGINS,
    RNA_STRAND_ORIGINS,
    OracleTruth,
    _split_bam,
    split_rna_by_transcript_strand,
)
from rigel.pipeline import scan_and_buffer  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402

#: The two axes ``CalibrationResult`` deconvolves; the sj axis is pure RNA and has nothing to score.
AXES = ("region", "boundary")


def object_fractions(gdna_mass, rna_mass) -> tuple[np.ndarray, np.ndarray]:
    """``(f_g, total)`` per object; ``f_g`` is NaN where the object carries no mass.

    NaN rather than 0 keeps "no data" out of the count of objects scored.
    """
    g = np.asarray(gdna_mass, np.float64)
    total = g + np.asarray(rna_mass, np.float64)
    frac = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(g, total, out=frac, where=total > 0.0)
    return frac, total


@dataclass(frozen=True, slots=True)
class AxisScore:
    """One arm's error over one axis. Every field except ``mwae`` is a mass, so it adds across a partition of
    the objects; the rate does not."""

    n_scored: int  #: objects with mass, never the object count of the axis
    mass: float  #: Σ total, the weight behind every number below
    net_err: float  #: Σ (gDNA_arm − gDNA_true). What the library-level figure sees
    abs_err: float  #: Σ |gDNA_arm − gDNA_true|. What the per-object answer is
    over_call: float  #: Σ (gDNA_arm − gDNA_true)+, gDNA claimed where there was RNA
    under_call: float  #: Σ (gDNA_true − gDNA_arm)+, gDNA missed
    mwae: float  #: mass-weighted mean |Δf_g| ≡ abs_err / mass. 0 = per-object perfect


def score_axis(arm_gdna, arm_rna, true_gdna, true_rna) -> AxisScore:
    """Score one arm against T on one axis.

    Refuses an arm and a truth on different bases: ``Σ w·|Δf_g| ≡ Σ|Δ gDNA mass|`` holds only when the two
    per-object totals agree.
    """
    arm_f, arm_total = object_fractions(arm_gdna, arm_rna)
    true_f, true_total = object_fractions(true_gdna, true_rna)
    if arm_total.shape != true_total.shape:
        raise ValueError(
            f"arm has {arm_total.shape[0]} objects, truth has {true_total.shape[0]} — these are "
            "different axes. A region mass scored against boundary truth is not a small error."
        )
    if not np.allclose(arm_total, true_total, rtol=1e-9, atol=1e-6):
        worst = int(np.argmax(np.abs(arm_total - true_total)))
        raise ValueError(
            f"arm and truth are on DIFFERENT BASES: per-object totals differ, worst at object "
            f"{worst} ({arm_total[worst]:.6g} vs {true_total[worst]:.6g}). Both must be "
            "contained count on the region axis and unspliced+spliced on the boundary axis."
        )

    live = np.isfinite(arm_f) & np.isfinite(true_f)
    d = (np.asarray(arm_gdna, np.float64) - np.asarray(true_gdna, np.float64))[live]
    mass = float(true_total[live].sum())
    abs_err = float(np.abs(d).sum())
    return AxisScore(
        n_scored=int(live.sum()),
        mass=mass,
        net_err=float(d.sum()),
        abs_err=abs_err,
        over_call=float(np.maximum(d, 0.0).sum()),
        under_call=float(np.maximum(-d, 0.0).sum()),
        mwae=abs_err / mass if mass > 0.0 else 0.0,
    )


def check_same_basis(name: str, arm, full_substrate) -> None:
    """Assert one ``CalibrationResult``-shaped arm's per-object totals are the payload's own totals.

    Checked per axis, so errors on the two cannot cancel. The region axis is ``region_contained``; the boundary
    axis is unspliced + spliced, matching the deconvolution's ``rna = (1−f_g)·unspliced + spliced``.
    """
    region_total = np.asarray(full_substrate.region_contained.count, np.float64).sum(axis=1)
    boundary_total = np.asarray(full_substrate.boundary_unspliced.count, np.float64).sum(
        axis=1
    ) + np.asarray(full_substrate.boundary_spliced.count, np.float64).sum(axis=1)
    for axis, expect in (("region", region_total), ("boundary", boundary_total)):
        got = np.asarray(getattr(arm, f"count_gdna_{axis}"), np.float64) + np.asarray(
            getattr(arm, f"count_rna_{axis}"), np.float64
        )
        if got.shape != expect.shape or not np.allclose(got, expect, rtol=1e-9, atol=1e-6):
            worst = int(np.argmax(np.abs(got - expect))) if got.shape == expect.shape else -1
            raise ValueError(
                f"{name}: the {axis} axis does not conserve the payload's own count"
                + (f" (worst at object {worst}: {got[worst]:.6g} vs {expect[worst]:.6g})"
                   if worst >= 0 else f" (shape {got.shape} vs {expect.shape})")
            )


def load_or_build_oracle(bam, index, pipeline_config, work_dir, tag, full_payload, cache_root,
                         lift=None, force=False):
    """T, from a per-origin cache when one is valid, otherwise split, scan, and populate it.

    The cache is keyed by the scan cache's own key, so ``read_scan_cache`` refuses one that does not describe this
    index and scan config; sum-to-full is re-run on every load. ``force`` rebuilds: the key hashes the deposit
    rule, not the scan's code, so a change to which fragments a scan offers needs it.
    """
    # ``full_payload`` is the drained whole and ``lift`` the drain's record; an empty ``lift`` drains nothing
    lift = lift or {}
    scan = pipeline_config.scan
    dirs = {k: Path(cache_root) / tag / k for k in ORIGINS}
    if not force:
        try:
            parts = {k: read_scan_cache(dirs[k], index, scan).payload for k in ORIGINS}
            truth = OracleTruth.from_cached_parts(full_payload, parts, lift)
        except (FileNotFoundError, KeyError, ScanCacheKeyError):
            pass  # no cache, or it does not describe this index/scan: rebuild it below
        else:
            # the three-way cache can be complete while the per-strand one is not
            ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root)
            return truth

    paths, read_counts = _split_bam(bam, Path(work_dir), tag)
    parts = {}
    for origin in ORIGINS:
        _stats, strand_model, _buf, payload = scan_and_buffer(paths[origin], index, scan)
        parts[origin] = payload
        # a part records the BAM it was split from, since the split file is scratch
        write_scan_cache(dirs[origin], payload=payload, strand_model=strand_model, index=index,
                         bam=bam, scan_config=scan)
    ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root, force=force)
    # the cache stays pass one (write_scan_cache refuses a drained payload); the truth is lifted as on a hit
    return OracleTruth.from_cached_parts(full_payload, parts, lift, read_counts)


def ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root, force=False) -> bool:
    """Cache the RNA reads split by transcript strand (``rna_pos`` / ``rna_neg``) beside the three ``ORIGINS``
    partitions, for the per-component truth. Skips the work when both load, unless ``force``; returns True if it
    built them.
    """
    dirs = {k: Path(cache_root) / tag / k for k in RNA_STRAND_ORIGINS}
    if not force:
        try:
            for k in RNA_STRAND_ORIGINS:
                read_scan_cache(dirs[k], index, scan)
            return False
        except (FileNotFoundError, KeyError, ScanCacheKeyError):
            pass
    paths, _counts = split_rna_by_transcript_strand(bam, Path(work_dir), tag)
    for key in RNA_STRAND_ORIGINS:
        _stats, strand_model, _buf, payload = scan_and_buffer(paths[key], index, scan)
        write_scan_cache(dirs[key], payload=payload, strand_model=strand_model, index=index,
                         bam=bam, scan_config=scan)  # the BAM it was split from
    return True


def library_f_gdna(result) -> float:
    """The library gDNA fraction, summed over both deconvolved axes.

    Applied to T as well as to the arm, so both count each fragment on every object it touched.
    """
    g = float(np.asarray(result.count_gdna_region).sum() + np.asarray(result.count_gdna_boundary).sum())
    r = float(np.asarray(result.count_rna_region).sum() + np.asarray(result.count_rna_boundary).sum())
    return g / (g + r) if (g + r) > 0 else 0.0
