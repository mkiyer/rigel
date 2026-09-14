#!/usr/bin/env python
"""How does pass-0 compare with the origin-split payload and two levered ceilings, per object and per class?

Pass-0 is calibration's prior-free first solve (``calib_refit_iters=0``). For each contaminated
condition this instrument scans the BAM once (or reads the cached scan), builds T, the truth, by
running the production accumulator on the BAM split by true origin (sum-to-full validated exactly,
in the drained frame: the whole drained at the production seed and every partition drained by
replaying its choices), runs P (pass-0) and the final solve, and scores each against T per object,
mass-weighted, on both deconvolved axes. The two ceilings are levers that already exist, never an
estimator: C_input runs the same solve with the simulator's own post-capture length pmfs, so
``C_input − P`` is how much of the error is wrong inputs rather than wrong solving; C_info is a
classification, never an estimate, of whether an object's own stored channels identify the
two-component split at all (`info_class_masks`), and it ignores neighbours, so ``C_info − P`` is not
a gap. Every object is also filed by the solver's own partition (`own_evidence` / `message_only` /
`struct_lock`), and the cross-tab "undetermined by C_info x answered by the messages" is the cell the
report exists to find. Scored per object and per class, never pooled; the directional split beside
the absolute; an object with no mass is absent, not ``f_g = 0``; zero-gDNA rows are held out because
a truth of exactly 0 is saturated. Its mass-weighted headline is the wrong yardstick for pass-0,
where honest ignorance reads as error; `solvability_audit.py` is that judge.

Also a library: `calibration_vs_oracle.py`, `solvability_audit.py` and `worst_objects.py`
import `measure_condition`, `score_axis`, `check_same_basis`, `object_fractions`,
`truth_length_pmf`, `truth_f_gdna`, `library_f_gdna`, the class tuples and the two defaults.

Usage::

    python scripts/design/pass0_vs_oracle.py --oracle-cache <dir>                    # the contaminated conditions
    python scripts/design/pass0_vs_oracle.py --conditions <name> --oracle-cache <dir>
    python scripts/design/pass0_vs_oracle.py --oracle-cache <dir> --jobs 4           # pre-warm the caches in parallel
    python scripts/design/pass0_vs_oracle.py --json out.json
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
import time
from dataclasses import dataclass, replace
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
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.effective_length import build_slot_moments  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.region_geometry import g1_locked  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag, scan_and_buffer  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
#: The current panel. Several instruments read this default from here; ``prior_vs_oracle.py``
#: carries its own copy, so a panel move must update both.
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: The solver's own "has own composition evidence" gate lives in production
#: (:func:`~rigel.calibration.region_init.has_own_composition_evidence`); this value exists only
#: because ``solver_slot_classes`` takes it as a movable argument so a gate can perturb the partition.
_EPS = 1.0e-9

#: The two axes ``CalibrationResult`` deconvolves. The sj axis is pure RNA by construction, so
#: nothing is deconvolved there and there is nothing to score.
AXES = ("region", "boundary")

#: Where did the answer come from? The solver's own three-way partition of a slot, reproducing
#: ``region_init``'s definitions. Mutually exclusive and exhaustive; a gate asserts the mass and the
#: error both decompose over them exactly.
SOLVER_CLASSES = ("own_evidence", "message_only", "struct_lock")

#: Is the answer there at all? The 2x2's identification status per object. ``absent`` is a class,
#: not a filter: an object with no mass has no answer to get right or wrong, and folding it into any
#: of the other three is how most of a genome reads as perfectly solved.
INFO_CLASSES = ("identified", "undet_no_separation", "undet_out_of_range", "absent")


# ── the per-object comparison ────────────────────────────────────────────────────────────────────


def object_fractions(gdna_mass, rna_mass) -> tuple[np.ndarray, np.ndarray]:
    """``(f_g, total)`` per object; ``f_g`` is NaN where the object carries no mass.

    NaN, never 0: "no data" must be inert, and a floored 0 reads as a confident "no gDNA here". The
    mass-weighted mean is blind to the difference (a zero-mass object carries zero weight), so the
    place it shows up is the count of objects scored and the class shares, which is exactly where an
    inflated denominator is invisible.
    """
    g = np.asarray(gdna_mass, np.float64)
    total = g + np.asarray(rna_mass, np.float64)
    frac = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(g, total, out=frac, where=total > 0.0)
    return frac, total


@dataclass(frozen=True, slots=True)
class AxisScore:
    """One arm's error over one selection of one axis. Every field is a mass, not a rate, except
    ``mwae``, so the fields add across a partition of the objects and the rate does not."""

    n_scored: int  #: objects with mass, never the object count of the axis
    mass: float  #: Σ total, the weight behind every number below
    net_err: float  #: Σ (gDNA_arm − gDNA_true). What the library-level figure sees
    abs_err: float  #: Σ |gDNA_arm − gDNA_true|. What the per-object answer is
    over_call: float  #: Σ (gDNA_arm − gDNA_true)+, gDNA claimed where there was RNA
    under_call: float  #: Σ (gDNA_true − gDNA_arm)+, gDNA missed
    mwae: float  #: mass-weighted mean |Δf_g| ≡ abs_err / mass. 0 = per-object perfect

    @property
    def cancellation(self) -> float:
        """``Σ|err| / |net|``: how much better the library-level number looks than the per-object
        answer, i.e. how much of an under-call is sitting next to an over-call."""
        return self.abs_err / abs(self.net_err) if self.net_err != 0.0 else float("inf")


def score_axis(arm_gdna, arm_rna, true_gdna, true_rna, select=None) -> AxisScore:
    """Score one arm against T on one axis, optionally restricted to ``select``.

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
    if select is not None:
        live &= np.asarray(select, bool)
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


# ── class 1: where did the answer come from? (the solver's own partition) ────────────────────────


def solver_slot_classes(capture, chain, eps: float = _EPS) -> dict[str, np.ndarray]:
    """Partition the chain's slots three ways, using ``region_init``'s own definitions.

    * ``struct_lock``: composition certain, on both axes
      (:func:`~rigel.calibration.region_geometry.g1_locked`). Neither RNA strand is admissible, so
      there is nothing to decide and ``f_g = 1`` is the pinned init. Locked is not the same as
      uninformed, and lumping the two reports a pure-gDNA intergenic region as a solver failure.
    * ``message_only``: no own composition evidence at all (``tau_lam`` at zero and not locked). Its
      gDNA/RNA split is decided entirely by neighbour messages and the population prior.
    * ``own_evidence``: everything else, where the strand Beta-Binomial or the intron factory's
      density deconvolution had something to say.

    ``eps`` is the solver's own gate (``own_composition_logvar`` tests ``tau > 1e-9``), so this
    partition answers "which mechanism did the solver use here", which is what the cross-tab and
    ``worst_objects.py`` need. It is deliberately not the question "should pass-0 be scored here": a
    fitted κ that misses ½ by a rounding step leaves a τ the solver treats as evidence that can
    resolve nothing, and that question is answered by ``solvability_audit``'s resolving-power curve
    over ``SD_LAMBDA_DECADES``. ``eps`` exists so a gate can move it and watch the partition move;
    production callers must not pass it.
    """
    tau = np.asarray(capture.tau_lam, np.float64)
    struct_lock = g1_locked(capture.free_pos, capture.free_neg)
    message_only = (tau <= eps) & (~struct_lock)
    return {
        "own_evidence": ~(struct_lock | message_only),
        "message_only": message_only,
        "struct_lock": struct_lock,
    }


def _project(slot_mask, chain, n_regions: int, n_boundaries: int) -> dict[str, np.ndarray]:
    """Scatter a per-slot boolean onto the region and boundary axes. The chain alternates
    REGION/BOUNDARY per reference, so every region and every contiguous boundary is exactly one slot
    and the map is a bijection; there is nothing to pool and nothing to drop."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    mask = np.asarray(slot_mask, bool)
    out = {"region": np.zeros(n_regions, bool), "boundary": np.zeros(n_boundaries, bool)}
    out["region"][obj[kind == REGION]] = mask[kind == REGION]
    out["boundary"][obj[kind == BOUNDARY]] = mask[kind == BOUNDARY]
    return out


def solver_class_masks(capture, chain, n_regions: int, n_boundaries: int) -> dict[str, dict]:
    """:func:`solver_slot_classes`, projected onto the two scored axes."""
    slots = solver_slot_classes(capture, chain)
    return {
        axis: {name: _project(m, chain, n_regions, n_boundaries)[axis] for name, m in slots.items()}
        for axis in AXES
    }


# ── class 2: is the answer THERE at all? (C_info) ────────────────────────────────────────────────


def info_class_masks(chain, region_arrays, substrate, gdna_pmf, rna_pmf) -> dict[str, dict]:
    """C_info: per object, is the two-component split recoverable from the stored channels?

    The 2×2 at one object is ``N = ρ_g·E_g + ρ_r·E_r`` and ``Σ1/L = ρ_g·D_g + ρ_r·D_r``, identified
    iff ``E_g/D_g ≠ E_r/D_r``, the two components' opportunity-weighted mean lengths. At a contiguous
    boundary the opportunity is ``(w−1)+`` and that reduces to the bare ``μ_g ≠ μ_r``; at a REGION it
    does not, because the opportunity is ``(ell − w + 1)+`` and ``1/w`` does not cancel it, so the
    moments are read in each slot's own frame from ``effective_length.build_slot_moments``, which
    already computes exactly this. Conditional on the object's count ``N``, the deposited weight has
    mean ``pi·m1_g + (1−pi)·m1_r`` where ``pi`` is the gDNA share of the landed fragments, so::

        pi_hat = (Σ(1/L)/N − m1_r) / (m1_g − m1_r)

    Four classes, and the middle two are the answer this function exists to give:

    * ``absent``: no count. There is no answer here to get right or wrong.
    * ``identified``: either only one component has any opportunity here (a region too short for any
      RNA fragment to fit inside can only contain gDNA), or both do, the moments separate, and
      ``pi_hat`` lands in ``[0, 1]``.
    * ``undet_no_separation``: ``m1_g == m1_r`` exactly. At equal opportunity-weighted means the
      channel carries zero information about composition, at any depth. Tested exactly, not to a
      tolerance: a difference computed from large floats comes out flat only to ~1e-11, and a 1e-11
      row would read as live and sell the grid's own width back as evidence. The degenerate "no
      opportunity for either component" case lands here without a special branch, because
      ``LandedMoments`` zeroes every moment at zero opportunity rather than flooring a division.
    * ``undet_out_of_range``: the moments separate but the solution is outside ``[0, 1]``, so the
      observation is not consistent with any mixture of these two components. This class also absorbs
      sampling noise (at one or two fragments a single draw of ``1/L`` easily falls outside the
      interval the two components span), which is why it is not pooled with the class above.

    This function ignores neighbours, which the sweep does not. It is a statement about one object's
    own channels and nothing else, and must not be read as a bound on the solver.
    """
    mg = build_slot_moments(chain, region_arrays, gdna_pmf)
    mr = build_slot_moments(chain, region_arrays, rna_pmf)

    kind = np.asarray(chain.kind)
    n_slots = kind.shape[0]
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    count = np.zeros(n_slots, np.float64)
    inv = np.zeros(n_slots, np.float64)
    for k, view in ((REGION, substrate.region_contained), (BOUNDARY, substrate.boundary_unspliced)):
        sel = kind == k
        count[sel] = np.asarray(view.count, np.float64).sum(axis=1)[obj[sel]]
        # two deposit rules, two names (TRAPS: two-masks-one-name): the REGION bank is the contained
        # rule's `inv_opportunity_sum`, the BOUNDARY bank the crossing rule's `inv_length_sum`.
        bank = view.inv_opportunity_sum if k == REGION else view.inv_length_sum
        inv[sel] = np.asarray(bank, np.float64)[obj[sel]]

    eff_g, eff_r = np.asarray(mg.eff, np.float64), np.asarray(mr.eff, np.float64)
    m1_g, m1_r = np.asarray(mg.m1, np.float64), np.asarray(mr.m1, np.float64)
    absent = count <= 0.0
    both = (eff_g > 0.0) & (eff_r > 0.0)
    one_only = (eff_g > 0.0) ^ (eff_r > 0.0)
    separates = both & (m1_g != m1_r)

    pi = np.full(n_slots, np.nan, np.float64)
    live = separates & ~absent
    np.divide(
        inv / np.where(absent, 1.0, count) - m1_r, m1_g - m1_r, out=pi, where=live
    )
    in_range = live & (pi >= 0.0) & (pi <= 1.0)
    slots = {
        "identified": (one_only & ~absent) | in_range,
        "undet_no_separation": ~absent & ~one_only & ~separates,
        "undet_out_of_range": live & ~in_range,
        "absent": absent,
    }
    n_regions, n_boundaries = int(substrate.n_regions), int(substrate.n_boundaries)
    return {
        axis: {name: _project(m, chain, n_regions, n_boundaries)[axis] for name, m in slots.items()}
        for axis in AXES
    }


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def calibrate_arm(payload, kwargs, config, *, gdna_pmf=None, rna_pmf=None, debug=None):
    """One ``calibrate`` run. ``gdna_pmf`` / ``rna_pmf`` override the fitted length models; that, and
    the solve depth in ``config``, are the only two things any arm varies."""
    call = dict(kwargs)
    if gdna_pmf is not None:
        call["gdna_fl_pmf"] = gdna_pmf
    if rna_pmf is not None:
        call["rna_fl_pmf"] = rna_pmf
    return calibrate(payload=payload, config=config, _debug=debug, **call)


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
    if cache_root is None:
        return OracleTruth.from_bam(
            bam, index, pipeline_config, Path(work_dir), tag, full_payload=full_payload,
            drain_with=((lift["undrained"], lift["choices"], lift["region_types"], lift["sj"])
                        if lift else None),
        )
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
    if cache_root is None:
        return False
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


@dataclass
class ConditionMeasurement:
    """Everything one condition produced, with the intermediates the gates interrogate."""

    condition: str
    payload: object
    oracle: OracleTruth
    truth: object  #: T, as a CalibrationResult
    arms: dict  #: name -> CalibrationResult
    calibrate_kwargs: dict
    debug_pass0: dict
    debug_final: dict
    scores: dict  #: arm -> axis -> class-or-"ALL" -> AxisScore
    info_scores: dict  #: arm -> axis -> info-class -> AxisScore
    info_shares: dict  #: axis -> info-class -> mass share
    #: The two classifications as boolean masks per axis, kept so a downstream instrument reads the
    #: same partition this one scored rather than recomputing its own. ``worst_objects.py`` consumes these.
    solver_masks: dict
    info_masks: dict
    cross: dict  #: axis -> (info class, solver class) -> AxisScore, for pass-0
    library_f_gdna: dict  #: "T" / arm name -> the library-level gDNA fraction (the thermometer)
    #: ``kind -> axis -> class -> objects in the class``, scored or not. Distinct from
    #: ``AxisScore.n_scored``, which counts only objects with mass; that distinction is the whole
    #: content of the ``absent`` class, whose scored count is 0 by definition.
    class_objects: dict
    #: Which length pmfs C_info was classified with, recorded rather than assumed: on a condition
    #: with no truth file it is the fitted ones.
    info_pmf_source: str
    seconds: float


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


def measure_condition(
    bam: str,
    index,
    pipeline_config,
    calibration_config,
    work_dir: Path,
    tag: str,
    *,
    truth_pmfs=None,
    oracle_cache=None,
) -> ConditionMeasurement:
    """Scan once (or read the cached scan), build T, run every arm, and score them per object and
    per class.

    ``truth_pmfs`` is a callable ``max_size -> (gdna_pmf, rna_pmf)``, or ``None`` to skip the C_input
    arms. A callable rather than two arrays because a pmf must be sized by the payload's own
    ``max_length``, which does not exist until the scan; handing ``calibrate`` two pmfs of different
    lengths is a silent frame mismatch rather than an error.
    """
    start = time.perf_counter()
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    # the main payload is cached too, by the same argument as the oracle cache: the scan depends only
    # on the BAM, the index and the scan config, never on calibration, so one cache serves every arm
    # of a campaign. Keyed by the shipped loader, never a home-made key: a refusal here is loud and
    # falls through to a rescan, where a home-made key would load a stale tally silently.
    _sc_dir = None if oracle_cache is None else Path(oracle_cache) / tag / "_main"
    payload = strand_model = None
    if _sc_dir is not None:
        try:
            _sc = read_scan_cache(_sc_dir, index, scan)
            payload, strand_model = _sc.payload, _sc.strand_model
        except (FileNotFoundError, KeyError, ScanCacheKeyError):
            payload = strand_model = None
    if payload is None:
        _stats, strand_model, _buffer, payload = scan_and_buffer(bam, index, scan)
        if _sc_dir is not None:
            write_scan_cache(_sc_dir, payload=payload, strand_model=strand_model, index=index,
                             bam=bam, scan_config=scan)

    # the drained frame: every arm below and T itself describe the tally production calibrates. The
    # drain replays at the production seed; the cache stays pass one.
    lift: dict = {}
    payload = _drain_side_buffer(
        payload, index, strand_model, seed=pipeline_config.second_pass_seed, _lift=lift
    )

    # T. Sum-to-full is validated on every bank exactly and raises if it does not hold, on the cached
    # path as well as the scanned one, so nothing below can run on an oracle that is not the
    # production payload split by origin. In the drained frame it is also the lift's identity gate.
    oracle = load_or_build_oracle(
        bam, index, pipeline_config, work_dir, tag, payload, oracle_cache, lift
    )

    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )

    from rigel.calibration.gdna_density import region_lengths_from_partition
    from rigel.calibration.splice_graph import build_region_partition_arrays

    max_size = int(payload.max_length)
    _fb, _fo, _frt = build_region_partition_arrays(index)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
        # production parity: the region args enable the two-pool contrast, exactly as pipeline.py
        region_lengths=region_lengths_from_partition(_fb, _fo, len(_frt)),
        region_types=_frt,
    )
    truth_gdna_pmf, truth_rna_pmf = truth_pmfs(max_size) if truth_pmfs is not None else (None, None)
    kwargs = dict(
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        sj=build_sj_geometry_arrays(index),
        boundary_flags=build_boundary_flags_array(index),
    )

    pass0_config = replace(calibration_config, calib_refit_iters=0)
    debug_pass0: dict = {}
    debug_final: dict = {}
    arms = {
        "pass0": calibrate_arm(payload, kwargs, pass0_config, debug=debug_pass0),
        "final": calibrate_arm(payload, kwargs, calibration_config, debug=debug_final),
    }
    if truth_gdna_pmf is not None and truth_rna_pmf is not None:
        for name, cfg in (("c_input_pass0", pass0_config), ("c_input_final", calibration_config)):
            arms[name] = calibrate_arm(
                payload, kwargs, cfg, gdna_pmf=truth_gdna_pmf, rna_pmf=truth_rna_pmf
            )

    truth = dataclasses.replace(arms["pass0"], **oracle.override_masses(ra))
    check_same_basis("T", truth, substrate)
    for name, arm in arms.items():
        check_same_basis(name, arm, substrate)

    # the classes come from pass-0's own run and are held fixed across every arm. ``tau_lam``
    # depends weakly on the incoming belief, so the final solve partitions slots slightly
    # differently, but "what evidence does this object have" is a property of the object, and a
    # class that moves between arms cannot be used to compare them.
    chain = debug_pass0["chain"]
    n_regions, n_boundaries = int(payload.n_regions), int(payload.n_boundaries)
    solver_masks = solver_class_masks(debug_pass0["capture"], chain, n_regions, n_boundaries)
    has_truth = truth_gdna_pmf is not None and truth_rna_pmf is not None
    info_pmf_source = "the simulator's own post-capture pmfs" if has_truth else "the FITTED pmfs"
    info_masks = info_class_masks(
        chain,
        ra,
        substrate,
        truth_gdna_pmf if has_truth else fl.gdna_pmf,
        truth_rna_pmf if has_truth else fl.rna_pmf,
    )

    def score_all(arm, masks, names):
        out = {}
        for axis in AXES:
            g = getattr(arm, f"count_gdna_{axis}")
            r = getattr(arm, f"count_rna_{axis}")
            tg = getattr(truth, f"count_gdna_{axis}")
            tr = getattr(truth, f"count_rna_{axis}")
            per = {"ALL": score_axis(g, r, tg, tr)}
            for name in names:
                per[name] = score_axis(g, r, tg, tr, select=masks[axis][name])
            out[axis] = per
        return out

    scores = {n: score_all(a, solver_masks, SOLVER_CLASSES) for n, a in arms.items()}
    info_scores = {n: score_all(a, info_masks, INFO_CLASSES) for n, a in arms.items()}
    info_shares = {
        axis: {
            name: (
                info_scores["pass0"][axis][name].mass / info_scores["pass0"][axis]["ALL"].mass
                if info_scores["pass0"][axis]["ALL"].mass > 0
                else 0.0
            )
            for name in INFO_CLASSES
        }
        for axis in AXES
    }

    # the cross-tab. Not a threshold on confidence: "undetermined by C_info × answered by the
    # messages" is a cell of a partition, and its error share is the statement.
    cross = {}
    p0 = arms["pass0"]
    for axis in AXES:
        g, r = getattr(p0, f"count_gdna_{axis}"), getattr(p0, f"count_rna_{axis}")
        tg, tr = getattr(truth, f"count_gdna_{axis}"), getattr(truth, f"count_rna_{axis}")
        cross[axis] = {
            (i, s): score_axis(g, r, tg, tr, select=info_masks[axis][i] & solver_masks[axis][s])
            for i in INFO_CLASSES
            for s in SOLVER_CLASSES
        }

    return ConditionMeasurement(
        condition=tag,
        payload=payload,
        oracle=oracle,
        truth=truth,
        arms=arms,
        calibrate_kwargs=kwargs,
        debug_pass0=debug_pass0,
        debug_final=debug_final,
        scores=scores,
        info_scores=info_scores,
        info_shares=info_shares,
        solver_masks=solver_masks,
        info_masks=info_masks,
        cross=cross,
        library_f_gdna={"T": library_f_gdna(truth)}
        | {name: library_f_gdna(arm) for name, arm in arms.items()},
        class_objects={
            kind: {
                axis: {name: int(np.asarray(m, bool).sum()) for name, m in masks[axis].items()}
                for axis in AXES
            }
            for kind, masks in (("solver", solver_masks), ("info", info_masks))
        },
        info_pmf_source=info_pmf_source,
        seconds=time.perf_counter() - start,
    )


# ── truth inputs ─────────────────────────────────────────────────────────────────────────────────


def truth_length_pmf(condition_dir: Path, kind: str, max_size: int):
    """The simulator's own post-capture length distribution for one origin class, as a pmf sized
    ``max_size + 1``; ``None`` when the condition has no truth file or no fragments of that kind.

    Post-capture empirical, not the configured ``frag_mean``: capture selects for length, so the
    configured parameters describe a library that was never sequenced.
    """
    path = condition_dir / "truth_fragment_lengths.tsv"
    if not path.is_file():
        return None
    pmf = np.zeros(max_size + 1, dtype=np.float64)
    with open(path) as handle:
        next(handle)
        for line in handle:
            row_kind, length_text, count_text, _fraction = line.rstrip("\n").split("\t")
            if row_kind != kind:
                continue
            length = int(length_text)
            if 0 <= length <= max_size:
                pmf[length] += float(count_text)
    total = pmf.sum()
    return pmf / total if total > 0 else None


def truth_f_gdna(condition_dir: Path):
    """The library's true gDNA fragment fraction, from the simulator's own origin counts; ``None``
    without a ``truth_summary.json``.

    From the file, never from the condition name: a name's gDNA token is a rate knob, not a fraction.
    """
    path = condition_dir / "truth_summary.json"
    if not path.is_file():
        return None
    counts = json.loads(path.read_text()).get("origin_counts", {})
    gdna = float(counts.get("gdna", 0.0))
    total = gdna + float(counts.get("mrna", 0.0)) + float(counts.get("nrna", 0.0))
    return gdna / total if total > 0 else None


# ── reporting ────────────────────────────────────────────────────────────────────────────────────

_ARM_LABEL = {
    "pass0": "P   pass-0 (prior-free)",
    "c_input_pass0": "C   pass-0, EXACT lengths",
    "final": "    final (3 refits)",
    "c_input_final": "C   final, EXACT lengths",
}
_ARM_ORDER = ("pass0", "c_input_pass0", "final", "c_input_final")


def _fmt(value: float) -> str:
    return f"{value:,.0f}"


def report(measurements: list[ConditionMeasurement]) -> None:
    print()
    print("═" * 112)
    print("⭐⭐ THE DELIVERABLE — pass-0, final, C_input and T, per object, on the CONTAMINATED rows")
    print("═" * 112)
    print(
        "   T is the production accumulator run on the BAM split by TRUE origin (sum-to-full "
        "validated, exactly).\n"
        "   C_input is the same solve handed the simulator's own post-capture length pmfs — a "
        "LENGTH-input ceiling.\n"
        "   ⚠ mass-weighted; Σ|err| is the per-object answer, net is what the library-level figure "
        "sees.\n"
        "   the DRAINED frame on every arm, T included: the partitions are lifted by replaying the whole's choices."
    )
    for m in measurements:
        true_contained = m.scores["pass0"]["region"]["ALL"].mass
        true_crossing = m.scores["pass0"]["boundary"]["ALL"].mass
        print()
        print(f"── {m.condition}    ({m.seconds:.0f} s)")
        print(
            f"   TRUE mass: region contained {_fmt(true_contained)}   "
            f"boundary crossing {_fmt(true_crossing)}   "
            f"true gDNA region {_fmt(float(np.asarray(m.truth.count_gdna_region).sum()))}"
        )
        for axis in AXES:
            print(f"   {axis.upper():<5} {'arm':<26} {'net err':>13} {'Σ|err|':>13} "
                  f"{'Σ|err|/|net|':>12} {'mwae':>8} {'under':>13} {'over':>13}")
            for name in _ARM_ORDER:
                if name not in m.scores:
                    continue
                s = m.scores[name][axis]["ALL"]
                print(
                    f"   {'':<5} {_ARM_LABEL[name]:<26} {s.net_err:>+13,.0f} {s.abs_err:>13,.0f} "
                    f"{s.cancellation:>12.1f} {s.mwae:>8.4f} {s.under_call:>13,.0f} "
                    f"{s.over_call:>13,.0f}"
                )
        t_lib = m.library_f_gdna["T"]
        print(f"   {'LIB':<5} {'arm':<26} {'f_gdna':>13} {'err vs T':>13}")
        print(f"   {'':<5} {'T   (the oracle)':<26} {t_lib:>13.4f} {'—':>13}")
        for name in _ARM_ORDER:
            if name not in m.library_f_gdna:
                continue
            value = m.library_f_gdna[name]
            print(f"   {'':<5} {_ARM_LABEL[name]:<26} {value:>13.4f} {value - t_lib:>+13.4f}")

    _report_classes(measurements)
    _report_info(measurements)
    _report_cross(measurements)


def _report_classes(measurements: list[ConditionMeasurement]) -> None:
    print()
    print("═" * 112)
    print("⭐ WHERE THE ERROR IS — pass-0, by the SOLVER's own three-way partition of each object")
    print("═" * 112)
    print("   own_evidence: strand / intron factory spoke here.  message_only: nothing did — "
          "the answer\n   came from neighbours and the population prior.  struct_lock: composition "
          "CERTAIN, nothing to decide.")
    for m in measurements:
        print(f"\n── {m.condition}")
        for axis in AXES:
            whole = m.scores["pass0"][axis]["ALL"]
            if whole.mass <= 0:
                continue
            print(f"   {axis.upper():<5} {'class':<16} {'in class':>10} {'w/ mass':>9} "
                  f"{'mass share':>11} {'Σ|err|':>13} {'err share':>10} {'mwae':>8} {'rate':>7}")
            for name in SOLVER_CLASSES:
                s = m.scores["pass0"][axis][name]
                mass_share = s.mass / whole.mass
                err_share = s.abs_err / whole.abs_err if whole.abs_err > 0 else 0.0
                # an empty class still prints its row: a report that only shows non-empty classes
                # can never show that a class it expected is missing.
                rate = f"{err_share / mass_share:.1f}x" if mass_share > 0 else "—"
                print(
                    f"   {'':<5} {name:<16} {m.class_objects['solver'][axis][name]:>10,} "
                    f"{s.n_scored:>9,} {mass_share:>10.1%} "
                    f"{s.abs_err:>13,.0f} {err_share:>9.1%} {s.mwae:>8.4f} {rate:>7}"
                )


def _report_info(measurements: list[ConditionMeasurement]) -> None:
    print()
    print("═" * 112)
    print("⭐ C_info — is the 2×2 IDENTIFIED at all, from this object's own stored channels?")
    print("═" * 112)
    print("   ⛔ NOT a bound on the solver: this ignores neighbour information, which the sweep "
          "uses. It says what\n   ONE OBJECT's own channels can and cannot resolve.")
    for m in measurements:
        print(f"\n── {m.condition}   (classified with {m.info_pmf_source})")
        for axis in AXES:
            whole = m.info_scores["pass0"][axis]["ALL"]
            if whole.mass <= 0:
                continue
            print(f"   {axis.upper():<5} {'class':<22} {'in class':>10} {'w/ mass':>9} "
                  f"{'mass share':>11} {'Σ|err| pass-0':>14} {'err share':>10} {'mwae':>8}")
            for name in INFO_CLASSES:
                s = m.info_scores["pass0"][axis][name]
                mass_share = s.mass / whole.mass
                err_share = s.abs_err / whole.abs_err if whole.abs_err > 0 else 0.0
                # ``absent`` scores 0 objects and 0 mass by definition; its content is the
                # "in class" column, which is why that column exists.
                print(
                    f"   {'':<5} {name:<22} {m.class_objects['info'][axis][name]:>10,} "
                    f"{s.n_scored:>9,} {mass_share:>10.1%} "
                    f"{s.abs_err:>14,.0f} {err_share:>9.1%} {s.mwae:>8.4f}"
                )


def _report_cross(measurements: list[ConditionMeasurement]) -> None:
    print()
    print("═" * 112)
    print("⭐⭐ THE CELL THIS SCRIPT EXISTS TO FIND — undetermined by C_info × answered by the MESSAGES")
    print("═" * 112)
    print("   An object with no own evidence AND no identifiable 2×2 has no answer of its own at "
          "all: whatever\n   pass-0 reports there came entirely from its neighbours and the "
          "population prior.")
    undetermined = ("undet_no_separation", "undet_out_of_range")
    for m in measurements:
        print(f"\n── {m.condition}")
        for axis in AXES:
            whole = m.scores["pass0"][axis]["ALL"]
            if whole.mass <= 0:
                continue
            print(f"   {axis.upper():<5} {'C_info × solver':<40} {'mass share':>11} "
                  f"{'Σ|err|':>13} {'err share':>10} {'mwae':>8}")
            for i in INFO_CLASSES:
                for s_name in SOLVER_CLASSES:
                    s = m.cross[axis][(i, s_name)]
                    if s.n_scored == 0:
                        continue
                    star = "⭐" if (i in undetermined and s_name == "message_only") else "  "
                    print(
                        f"   {'':<5} {star} {i + ' × ' + s_name:<37} "
                        f"{s.mass / whole.mass:>10.1%} {s.abs_err:>13,.0f} "
                        f"{(s.abs_err / whole.abs_err if whole.abs_err > 0 else 0):>9.1%} "
                        f"{s.mwae:>8.4f}"
                    )


def _to_json(measurements: list[ConditionMeasurement]) -> list[dict]:
    rows = []
    for m in measurements:
        rows.append(
            {
                "condition": m.condition,
                "seconds": m.seconds,
                "scores": {
                    arm: {
                        axis: {k: dataclasses.asdict(v) for k, v in per.items()}
                        for axis, per in by_axis.items()
                    }
                    for arm, by_axis in m.scores.items()
                },
                "info_scores": {
                    arm: {
                        axis: {k: dataclasses.asdict(v) for k, v in per.items()}
                        for axis, per in by_axis.items()
                    }
                    for arm, by_axis in m.info_scores.items()
                },
                "info_shares": m.info_shares,
                "library_f_gdna": m.library_f_gdna,
                "class_objects": m.class_objects,
                "info_pmf_source": m.info_pmf_source,
                "cross": {
                    axis: {f"{i}|{s}": dataclasses.asdict(v) for (i, s), v in cells.items()}
                    for axis, cells in m.cross.items()
                },
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")))
    ap.add_argument(
        "--oracle-cache",
        type=Path,
        default=None,
        help="persist the per-origin oracle payloads here and reuse them. ⭐ The oracle depends on "
        "the ACCUMULATOR and the index, never on calibration, so one cache serves an entire "
        "solver-debugging campaign. Keyed by the scan cache's own key, so a stale one is refused.",
    )
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument(
        "--jobs", type=int, default=1,
        help="PRE-WARM the per-origin oracle caches with this many worker processes. ⭐ It parallelises "
             "ONLY the cache build — the BAM split plus four scans per condition, which is ~95 %% of the "
             "wall clock — and leaves the measurement loop serial and in its original order, so the "
             "numbers are bit-identical to --jobs 1. ⚠ Core-bound, not memory-bound: one worker "
             "saturates ONE core at ~2 GB of real memory (a condition's ~45 GB RSS is reclaimable page "
             "cache, not a footprint).")
    ap.add_argument("--_prewarm", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if not args.suite.is_dir():
        print(f"no suite dir at {args.suite}", file=sys.stderr)
        return 2

    names = args.conditions or sorted(p.name for p in args.suite.iterdir() if p.is_dir())
    scored: list[str] = []
    zero_rows: list[str] = []
    for name in names:
        truth = truth_f_gdna(args.suite / name)
        if truth is None:
            continue
        # only an exactly-zero row is saturated: at truth = 0 any change that lowers the estimate
        # scores better. At truth = 0.01 it does not, and a 1 % row is a perfectly good scoring row;
        # a threshold that encodes one panel's value set would drop the next panel's low-gDNA end.
        if truth <= 0.0:
            zero_rows.append(name)
            continue
        scored.append(name)

    if zero_rows:
        print(f"  {len(zero_rows)} zero-gDNA row(s) held out as false-positive checks: "
              f"{', '.join(zero_rows)}")
    if not scored:
        print("no contaminated conditions found", file=sys.stderr)
        return 2

    index = TranscriptIndex.load(str(args.index))
    pipeline_config = PipelineConfig()
    calibration_config = CalibrationConfig()
    work_dir = args.work_dir / "rigel_pass0_oracle"

    # one condition's cache, then exit: the worker half of `--jobs`. Not part of the public CLI.
    if args._prewarm is not None:
        c = args._prewarm
        payload = read_scan_cache(Path(args.suite) / "scan_cache" / c, index).payload
        load_or_build_oracle(str(Path(args.suite) / c / "sim_oracle.bam"), index, pipeline_config,
                             work_dir / f"w_{c}", c, payload, args.oracle_cache)
        return 0

    # ── pre-warm in parallel, then measure serially off warm caches ──────────────────────────────
    # The split is what makes `--jobs` safe: `load_or_build_oracle` is a pure cache fill keyed by the
    # shipped loader, so filling it in another process cannot move a number (`measure_condition`
    # re-runs the sum-to-full gate over whatever it loads), and the scoring loop stays serial and in
    # order, so `--jobs N` is bit-identical to `--jobs 1`. A worker that fails is not fatal: the
    # serial loop rebuilds that condition itself, same code path.
    if args.jobs > 1 and args.oracle_cache is not None and len(scored) > 1:
        import concurrent.futures as _cf
        import subprocess as _sp
        todo = [c for c in scored
                if not all((Path(args.oracle_cache) / c / k / "payload.npz").is_file() for k in ORIGINS)]
        if todo:
            n = max(1, min(int(args.jobs), len(todo)))
            print(f"  pre-warming {len(todo)} oracle cache(s), {n} worker(s) …", flush=True)
            base = [sys.executable, str(Path(__file__).resolve()), "--suite", str(args.suite),
                    "--index", str(args.index), "--oracle-cache", str(args.oracle_cache),
                    "--work-dir", str(args.work_dir)]
            def _one(c):
                return c, _sp.run(base + ["--_prewarm", c], capture_output=True, text=True).returncode
            with _cf.ThreadPoolExecutor(max_workers=n) as ex:
                for c, rc in ex.map(_one, todo):
                    print(f"    {'✔' if rc == 0 else '⚠ serial loop will rebuild'} {c}", flush=True)

    measurements = []
    for name in scored:
        bam = str(args.suite / name / "sim_oracle.bam")
        print(f"  {name} …", flush=True)
        cond_dir = args.suite / name
        m = measure_condition(
            bam=bam,
            index=index,
            pipeline_config=pipeline_config,
            calibration_config=calibration_config,
            work_dir=work_dir,
            tag=name,
            truth_pmfs=lambda size, d=cond_dir: (
                truth_length_pmf(d, "gdna", size),
                truth_length_pmf(d, "rna", size),
            ),
            oracle_cache=args.oracle_cache,
        )
        measurements.append(m)
        print(f"  {name} done in {m.seconds:.0f} s", flush=True)

    report(measurements)
    if args.json:
        args.json.write_text(json.dumps(_to_json(measurements), indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
