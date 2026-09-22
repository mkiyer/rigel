"""The oracle arms — one condition scanned, split by true origin, calibrated at pass-0 and in full, and scored
per object against T. A helper (no row on the shelf), shared by `solvability_audit.py`, `calibration_vs_oracle.py`
and `calibration_oracle.py --build`, and gated by `tests/calibration/test_oracle_arms.py`. It was the instrument
`pass0_vs_oracle.py` until 2026-09-22; its own report — the C_input ceiling and the C_info classification, both
pricing the retired length channel — went with it.

T is the production accumulator run on the BAM split by origin and asserted to sum to the full payload, in the
DRAINED frame (the partitions are lifted by replaying the whole's choices, `lift_drain_parts`). The arms are
`pass0` (refits 0) and `final` (the shipped solve). `load_or_build_oracle` is the one place the per-origin caches
are built — keyed by the shipped `read_scan_cache`, so a stale cache is refused rather than reused, and the
sum-to-full identity is re-run on every load.
"""

from __future__ import annotations

import dataclasses
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
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.region_geometry import g1_locked  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag, scan_and_buffer  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402

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


def solver_slot_classes(capture, eps: float = _EPS) -> dict[str, np.ndarray]:
    """Partition the chain's slots three ways, using ``region_init``'s own definitions.

    * ``struct_lock``: composition certain, on both axes
      (:func:`~rigel.calibration.region_geometry.g1_locked`). Neither RNA strand is admissible, so
      there is nothing to decide and ``f_g = 1`` is the pinned init. Locked is not the same as
      uninformed, and lumping the two reports a pure-gDNA intergenic region as a solver failure.
    * ``message_only``: no own composition evidence at all (``tau_lam`` at zero and not locked). Its
      gDNA/RNA split is decided entirely by neighbour messages and the population prior.
    * ``own_evidence``: everything else, where the strand Beta-Binomial or the intron factory's
      density deconvolution had something to say.

    ``eps`` is the solver's own gate (``has_own_composition_evidence`` tests ``tau_lam > 1e-9``), so this
    partition answers "which mechanism did the solver use here", which is what the cross-tab and
    a dissection needs. It is deliberately not the question "should pass-0 be scored here": a
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
    slots = solver_slot_classes(capture)
    return {
        axis: {name: _project(m, chain, n_regions, n_boundaries)[axis] for name, m in slots.items()}
        for axis in AXES
    }


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
    arms: dict  #: "pass0" / "final" -> CalibrationResult
    calibrate_kwargs: dict
    debug_pass0: dict
    debug_final: dict
    scores: dict  #: arm -> axis -> solver class or "ALL" -> AxisScore
    #: The solver's classification as boolean masks per axis, kept so a downstream instrument reads
    #: the same partition this one scored rather than recomputing its own.
    solver_masks: dict
    seconds: float


def measure_condition(
    bam: str,
    index,
    pipeline_config,
    calibration_config,
    work_dir: Path,
    tag: str,
    *,
    oracle_cache=None,
) -> ConditionMeasurement:
    """Scan once (or read the cached scan), build T, run the two arms, and score them per object and
    per solver class."""
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
    from rigel.calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays
    from rigel.pipeline import library_fl_models

    fl = library_fl_models(payload, index)
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

    def score_all(arm):
        out = {}
        for axis in AXES:
            g = getattr(arm, f"count_gdna_{axis}")
            r = getattr(arm, f"count_rna_{axis}")
            tg = getattr(truth, f"count_gdna_{axis}")
            tr = getattr(truth, f"count_rna_{axis}")
            per = {"ALL": score_axis(g, r, tg, tr)}
            for name in SOLVER_CLASSES:
                per[name] = score_axis(g, r, tg, tr, select=solver_masks[axis][name])
            out[axis] = per
        return out

    return ConditionMeasurement(
        condition=tag,
        payload=payload,
        oracle=oracle,
        truth=truth,
        arms=arms,
        calibrate_kwargs=kwargs,
        debug_pass0=debug_pass0,
        debug_final=debug_final,
        scores={n: score_all(a) for n, a in arms.items()},
        solver_masks=solver_masks,
        seconds=time.perf_counter() - start,
    )


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
