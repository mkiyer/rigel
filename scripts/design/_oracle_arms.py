"""The oracle cache and the per-object scorer. A helper (no row on the shelf), shared by
`calibration_vs_oracle.py` (the scorer) and `calibration_oracle.py --build` (the cache), and gated by
`tests/calibration/test_oracle_arms.py`.

T is the production accumulator run on the BAM split by origin and asserted to sum to the full payload, in the
DRAINED frame (the partitions are lifted by replaying the whole's choices, `lift_drain_parts`).
`load_or_build_oracle` is the one place the per-origin caches are built — keyed by the shipped `read_scan_cache`,
so a stale cache is refused rather than reused, and the sum-to-full identity is re-run on every load.
"""

from __future__ import annotations

import dataclasses
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
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402

#: The two axes ``CalibrationResult`` deconvolves. The sj axis is pure RNA by construction, so
#: nothing is deconvolved there and there is nothing to score.
AXES = ("region", "boundary")


def object_fractions(gdna_mass, rna_mass) -> tuple[np.ndarray, np.ndarray]:
    """``(f_g, total)`` per object; ``f_g`` is NaN where the object carries no mass.

    NaN, never 0: "no data" must be inert, and a floored 0 reads as a confident "no gDNA here". The
    mass-weighted mean is blind to the difference (a zero-mass object carries zero weight), so the
    place it shows up is the count of objects scored, which is exactly where an inflated denominator is
    invisible.
    """
    g = np.asarray(gdna_mass, np.float64)
    total = g + np.asarray(rna_mass, np.float64)
    frac = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(g, total, out=frac, where=total > 0.0)
    return frac, total


@dataclass(frozen=True, slots=True)
class AxisScore:
    """One arm's error over one axis. Every field is a mass, not a rate, except
    ``mwae``, so the fields add across a partition of the objects and the rate does not."""

    n_scored: int  #: objects with mass, never the object count of the axis
    mass: float  #: Σ total, the weight behind every number below
    net_err: float  #: Σ (gDNA_arm − gDNA_true). What the library-level figure sees
    abs_err: float  #: Σ |gDNA_arm − gDNA_true|. What the per-object answer is
    over_call: float  #: Σ (gDNA_arm − gDNA_true)+, gDNA claimed where there was RNA
    under_call: float  #: Σ (gDNA_true − gDNA_arm)+, gDNA missed
    mwae: float  #: mass-weighted mean |Δf_g| ≡ abs_err / mass. 0 = per-object perfect


def score_axis(arm_gdna, arm_rna, true_gdna, true_rna) -> AxisScore:
    """Score one arm against T on one axis.

    Refuses arms and truths on different bases: ``Σ w·|Δf_g| ≡ Σ|Δ gDNA mass|`` holds only when the
    two per-object totals agree, and without that identity the mass-weighted mean is a weighted
    average of fractions over different denominators.
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

    Per axis, never pooled: ``n_regions`` and ``n_boundaries`` differ by only ``n_refs``, so an error
    on one axis cancelling an equal and opposite one on the other is a real class of mistake. The
    region axis holds no spliced molecule (``region_contained`` is credited only when the fragment
    used no sj); the boundary axis is unspliced + spliced, because the boundary deconvolution builds
    ``rna = (1−f_g)·unspliced + spliced`` and T must match that or the two are different quantities.
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
                         lift=None):
    """T, from a per-origin cache when one is valid, otherwise split, scan, and populate it.

    The oracle depends only on the accumulator and the index, so it is invariant across every solver
    change the debug loop makes: the cache is written once and hits for the rest of a campaign.
    Keyed by the scan cache's own key, never a new one: ``read_scan_cache`` refuses a payload whose
    ``graph_hash``, ``reach_digest``, ``payload_schema_digest`` or scan config does not describe the
    index it is loaded against (``reach`` in particular is covered by no other hash), so a stale
    oracle is refused loudly rather than silently feeding everything downstream. The sum-to-full
    identity is re-run over the loaded arrays regardless: the cache skips the scanning, never the
    validation.
    """
    # the drained frame: ``full_payload`` is the drained whole and ``lift`` the drain's box; the
    # partitions are drained by replaying the whole's choices. An empty ``lift`` (no held fragments)
    # keeps the pass-one identity unchanged.
    lift = lift or {}
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    dirs = {k: Path(cache_root) / tag / k for k in ORIGINS}
    try:
        parts = {k: read_scan_cache(dirs[k], index, scan).payload for k in ORIGINS}
        truth = OracleTruth.from_cached_parts(full_payload, parts, lift)
    except (FileNotFoundError, KeyError, ScanCacheKeyError):
        pass  # no cache, or it does not describe this index/scan: rebuild it below
    else:
        # the three-way cache can be complete while the per-strand one is not, so it is ensured on
        # the hit path too rather than only when the three are rebuilt
        ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root)
        return truth

    paths, read_counts = _split_bam(bam, Path(work_dir), tag)
    parts = {}
    for origin in ORIGINS:
        _stats, strand_model, _buf, payload = scan_and_buffer(paths[origin], index, scan)
        parts[origin] = payload
        write_scan_cache(dirs[origin], payload=payload, strand_model=strand_model, index=index,
                         bam=paths[origin], scan_config=scan)
    ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root)
    # the cache stays pass one (write_scan_cache refuses a drained payload); the returned truth is
    # lifted into the drained frame exactly as the cached path is.
    return OracleTruth.from_cached_parts(full_payload, parts, lift, read_counts)


def ensure_rna_strand_cache(bam, index, scan, work_dir, tag, cache_root) -> bool:
    """Cache the RNA reads split by transcript strand, ``rna_pos`` / ``rna_neg`` beside the three
    ``ORIGINS`` partitions, which is what a per-component truth needs (`_oracle.RNA_STRAND_ORIGINS`
    says why the payload's genome-strand columns cannot serve). Additive: the two partitions describe
    the same reads as ``mrna`` + ``nrna``, so ``calibration_oracle.py`` gates them against each other.
    Skips the work when both caches already load. Returns True if it built them.
    """
    dirs = {k: Path(cache_root) / tag / k for k in RNA_STRAND_ORIGINS}
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
                         bam=paths[key], scan_config=scan)
    return True


def library_f_gdna(result) -> float:
    """The library gDNA fraction, summed over both deconvolved axes, the reported deliverable.

    Printed beside the per-object answer because the gap between the two is the finding: the library
    figure is ``|Σ(g − t)|`` and the per-object answer is ``Σ|g − t|``. Applying the same functional
    to T rather than to the simulator's origin counts keeps the comparison on the accumulator's own
    basis (the origin-count fraction counts fragments once; this counts each fragment on every object
    it touched). Both axes, always: summing one axis reports a library's gDNA as a fraction of part
    of itself.
    """
    g = float(np.asarray(result.count_gdna_region).sum() + np.asarray(result.count_gdna_boundary).sum())
    r = float(np.asarray(result.count_rna_region).sum() + np.asarray(result.count_rna_boundary).sum())
    return g / (g + r) if (g + r) > 0 else 0.0
