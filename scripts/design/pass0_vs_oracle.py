#!/usr/bin/env python
"""⭐⭐ **STAGE B STEP 1** — pass-0 against the oracle, per object and per class.

**The question.** Pass-0 is calibration's PRIOR-FREE first solve, the first thing that happens after
initialisation and before any fitted prior exists. Does it find what the payload actually contains?
And where it does not, is that because the accumulator destroyed the information, or because the
solver missed it? A single library-level number cannot separate those two, and this project has
already spent sessions on the wrong one of them.

⭐ **THREE QUANTITIES PER OBJECT, AND THE DIFFERENCES BETWEEN THEM ARE THE WHOLE DIAGNOSTIC.**

===  =========================================================  ==========================================
 T   the truth: each object's real gDNA / RNA split             ``tests/calibration/_oracle.py``
 C   a ceiling: the best answer reachable under stated inputs   two LEVERS, below — never an estimator
 P   what pass-0 produces                                       ``calib_refit_iters=0``
===  =========================================================  ==========================================

``T − C`` is information the accumulator destroyed → Stage A work. ``C − P`` is solver gap → Stage B.

⛔ **THERE IS NO "BEST POSSIBLE PER-OBJECT ESTIMATOR" HERE, AND WRITING ONE WOULD BE THE MAGIC-NUMBER
FAILURE MODE WITH AN ESTIMATOR IN PLACE OF A CONSTANT.** So C is TWO things, each defined by a lever
that already exists, and they are reported separately:

* **C_input — the perfect-INPUT ceiling.** ``calibrate`` run with the simulator's own post-capture
  length distributions in place of the fitted ones (the same override ``calibration_truth_ab.py
  --ceiling`` uses), at both solve depths so each arm is compared against its own ceiling.
  ``C_input − P`` is "how much of the error is wrong INPUTS rather than wrong SOLVING?".
  ⚠ It is a *length*-input ceiling. The other library-level inputs (κ, the two strand
  overdispersions, the two Fisher sample sizes, the enrichment landscape, the two backgrounds) are
  injectable via ``InjectedCalibrationPriors`` — but the simulator writes no truth for them, so
  injecting anything there would be an A/B against a guess, not a ceiling. The one with a truth
  value is κ, and it is not a free lever either: protocol *fidelity* and a directional *sense
  fraction* are different quantities and a fitted κ of 0.0101 on a "0.99 stranded" library is
  correct. That is the next lever and it is a separate measurement.
* **C_info — the INFORMATION ceiling, reported as a CLASSIFICATION and never as an estimate.** Per
  object: is the two-component split recoverable from the stored channels at all? See
  :func:`info_class_masks`.

⚠ **C_info deliberately ignores neighbour information, which the sweep does use.** So C_info can be
"worse" than P on an object its neighbours rescued. ⛔ ``C_info − P`` is therefore NOT a gap and is
not reported as one. The useful statement is the reverse, and it is the point of the whole script:
**objects that C_info says are undetermined but that P answers anyway are objects whose answer came
entirely from the relay** — the class that carried 92 % of all error in the last full sweep.

SCORING RULES, NON-NEGOTIABLE
-----------------------------
* **Per object AND per class, never pooled.** One run once reported a −47.6 k error next to a +75.4 k
  error as a "nearly perfect" −1.9 k.
* **Mass-weighted**, because that is how the consumer weights it — a 1 bp node with 2 fragments must
  not count like an exon with 40,000.
* **The directional split alongside the absolute**, because here they nearly cancel.
* **The contaminated conditions only.** A zero-gDNA row is saturated at truth = 0 exactly, so
  anything that lowers the estimate "improves" it. Those rows are false-positive checks, nothing more,
  and this script refuses to average them in.
* **An object with no mass is ABSENT, not ``f_g = 0``.** Most of any real index carries no fragments.

⛔ **NO DRAIN, AND THAT IS FORCED, NOT AN OVERSIGHT.** T is the production accumulator run on the BAM
split by true origin, and the second pass's multinomial draw is scored against the payload's own
densities — so draining three partitions separately is not the same operation as draining the whole,
and the sum-to-full identity that makes the oracle trustworthy would not survive it. Both T and every
P arm here are therefore pass-one tallies. The drain is measured against per-fragment truth by
``second_pass_accuracy.py`` and against the deliverable by ``calibration_truth_ab.py``; it is a
different axis from this one.

⚠ **AND IT HAS ONE CONSEQUENCE WORTH STATING, BECAUSE IT BIASES EVERY P ARM HERE.** The fitted length
models are built from the UNDRAINED payload, and the drain is what fixes them: the RNA pool moves
**−4.0 % → −0.21 %** off capture and **−1.5 % → +1.11 %** under it once the held (systematically LONG)
fragments are returned. So the P arms are fed worse length models than production ships. ⭐ The direction
is conservative for the conclusion drawn from ``C_input − P``: a worse fitted input can only make
perfecting it look *more* valuable, and it still measures as worth ≈ 0. It would NOT be conservative for
any claim about how good pass-0 is in absolute terms.

Gates: ``tests/calibration/test_pass0_vs_oracle.py`` (9, each carrying its own perturbation).

Usage::

    python scripts/design/pass0_vs_oracle.py                      # the 4 contaminated conditions
    python scripts/design/pass0_vs_oracle.py --conditions NAME    # one
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

from _oracle import ORIGINS, OracleTruth, _split_bam  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.effective_length import build_slot_moments  # noqa: E402
from rigel.calibration.node_chain import EDGE, NODE  # noqa: E402
from rigel.calibration.node_geometry import g1_locked  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "pilot"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: ⭐ NOT DEFINED HERE ANY MORE. "Has own composition evidence" has ONE home in production —
#: :func:`~rigel.calibration.node_init.has_own_composition_evidence` — and every instrument imports
#: it. It used to be restated here and in ``composition_evidence_census.py``, each beside a comment
#: saying it must match the solver, which is precisely the arrangement TRAPS: a-test-that-redefines names: a change
#: to the solver would have moved neither. Kept as a name only because ``solver_slot_classes`` takes
#: it as a MOVEABLE argument so a gate can perturb the partition and watch it move.
_EPS = 1.0e-9

#: The two axes ``CalibrationResult`` deconvolves. The junction axis is pure RNA by construction —
#: nothing is deconvolved there, so there is nothing to score.
AXES = ("node", "edge")

#: ⭐ WHERE DID THE ANSWER COME FROM? The solver's own three-way partition of a slot, reproducing
#: ``node_init``'s definitions. Mutually exclusive and exhaustive — a gate asserts the mass and the
#: error both decompose over them exactly.
SOLVER_CLASSES = ("own_evidence", "relay_only", "struct_lock")

#: ⭐ IS THE ANSWER THERE AT ALL? The 2x2's identification status per object. ``absent`` is a class,
#: not a filter: an object with no mass has no answer to get right or wrong, and folding it into any
#: of the other three is how 80 % of a genome reads as perfectly solved.
INFO_CLASSES = ("identified", "undet_no_separation", "undet_out_of_range", "absent")


# ── the per-object comparison ────────────────────────────────────────────────────────────────────


def object_fractions(gdna_mass, rna_mass) -> tuple[np.ndarray, np.ndarray]:
    """``(f_g, total)`` per object — ``f_g`` is **NaN** where the object carries no mass.

    ⛔ NaN, never 0. "No data" must be inert; a floored 0 reads as a confident "no gDNA here" and its
    mirror reads as a confident "all gDNA here", and one of those was actively seeding false gDNA
    into neighbouring exons. The mass-weighted mean is blind to the difference (a zero-mass object
    carries zero weight), so the place it shows up is the COUNT of objects scored and the class
    shares — which is exactly where an inflated denominator is invisible.
    """
    g = np.asarray(gdna_mass, np.float64)
    total = g + np.asarray(rna_mass, np.float64)
    frac = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(g, total, out=frac, where=total > 0.0)
    return frac, total


@dataclass(frozen=True, slots=True)
class AxisScore:
    """One arm's error over one selection of one axis. Every field is a mass, not a rate, except
    ``mwae`` — so the fields ADD across a partition of the objects and the rate does not."""

    n_scored: int  #: objects with mass. ⚠ never the object count of the axis
    mass: float  #: Σ total, the weight behind every number below
    net_err: float  #: Σ (gDNA_arm − gDNA_true). What the library-level figure sees
    abs_err: float  #: Σ |gDNA_arm − gDNA_true|. What the per-object answer is
    over_call: float  #: Σ (gDNA_arm − gDNA_true)+ — gDNA claimed where there was RNA
    under_call: float  #: Σ (gDNA_true − gDNA_arm)+ — gDNA missed
    mwae: float  #: mass-weighted mean |Δf_g| ≡ abs_err / mass. 0 = per-object perfect

    @property
    def cancellation(self) -> float:
        """``Σ|err| / |net|`` — how much better the library-level number looks than the per-object
        answer. ⭐ This is the hook for the whole of Stage B: at 5.9x, the deliverable is reporting a
        near-miss over a large under-call sitting next to a large over-call."""
        return self.abs_err / abs(self.net_err) if self.net_err != 0.0 else float("inf")


def score_axis(arm_gdna, arm_rna, true_gdna, true_rna, select=None) -> AxisScore:
    """Score one arm against T on one axis, optionally restricted to ``select``.

    ⛔ **Refuses arms and truths on different bases.** ``Σ w·|Δf_g| ≡ Σ|Δ gDNA mass|`` holds *only*
    when the two per-object totals agree; without that identity the mass-weighted mean is a weighted
    average of fractions over different denominators, which is a number with no interpretation.
    """
    arm_f, arm_total = object_fractions(arm_gdna, arm_rna)
    true_f, true_total = object_fractions(true_gdna, true_rna)
    if arm_total.shape != true_total.shape:
        raise ValueError(
            f"arm has {arm_total.shape[0]} objects, truth has {true_total.shape[0]} — these are "
            "different axes. A node mass scored against edge truth is not a small error."
        )
    if not np.allclose(arm_total, true_total, rtol=1e-9, atol=1e-6):
        worst = int(np.argmax(np.abs(arm_total - true_total)))
        raise ValueError(
            f"arm and truth are on DIFFERENT BASES: per-object totals differ, worst at object "
            f"{worst} ({arm_total[worst]:.6g} vs {true_total[worst]:.6g}). Both must be "
            "contained count on the node axis and unspliced+spliced on the edge axis."
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

    ⚠ **Per axis, never pooled.** ``n_nodes`` and ``n_edges`` differ by only ``n_refs``, so an error
    on one axis cancelling an equal and opposite one on the other is not far-fetched — that is a real
    class of mistake a three-axis schema makes possible, and pooling the check invites it.

    The node axis holds no spliced molecule (``node_contained`` is credited only when the fragment
    used no junction); the edge axis is unspliced + spliced, because ``chain_edge_deconv`` builds
    ``rna = (1−f_g)·unspliced + spliced`` and T must match that or the two are different quantities.
    """
    node_total = np.asarray(full_substrate.node_contained.count, np.float64).sum(axis=1)
    edge_total = np.asarray(full_substrate.edge_unspliced.count, np.float64).sum(
        axis=1
    ) + np.asarray(full_substrate.edge_spliced.count, np.float64).sum(axis=1)
    for axis, expect in (("node", node_total), ("edge", edge_total)):
        got = np.asarray(getattr(arm, f"mass_gdna_{axis}"), np.float64) + np.asarray(
            getattr(arm, f"mass_rna_{axis}"), np.float64
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
    """Partition the chain's SLOTS three ways, using ``node_init``'s own definitions.

    * ``struct_lock`` — composition CERTAIN. A slot the solver could not move because there is
      nothing to decide: neither RNA strand admissible, so ``_type_belief`` pins ``{0,0,1}`` at
      ``Var(log f_g) = 0``. ⚠ Locked is not the same as uninformed, and lumping the two reports a
      pure-gDNA intergenic node as a solver failure.
      ⛔ **BOTH AXES** — :func:`~rigel.calibration.node_geometry.g1_locked`. This was
      ``(~solvable) & (kind == NODE)``, so every structurally-locked EDGE — an intergenic↔exon seam,
      where RNA cannot cross a gene boundary — was filed as ``relay_only``, i.e. as an object whose
      answer came from its neighbours, when nothing was ever asked of it and its ``f_g = 1`` is the
      pinned init.
      ⚠⚠ **AND IT IS NOT THE SAME MASK AS ``node_init.strand_evidence``'s ``struct_lock``**, which is
      node-only ON PURPOSE. That one answers "may this slot EMIT composition certainty into its
      messages?" and excludes G1 edges because a seam is structurally gDNA yet sits between
      RNA-carrying exons, so certainty there compounds into a phantom-gDNA emitter. Two questions, one
      word; ``g1_locked``'s docstring holds the distinction.
    * ``relay_only`` — no own composition evidence at all (``tau_lam`` at zero and not locked). Its
      gDNA/RNA split is decided entirely by neighbour messages and the population prior.
    * ``own_evidence`` — everything else: the strand Beta-Binomial, the intron factory's density
      deconvolution, or the length channel had something to say here.

    ⚠⚠ **``eps`` IS THE SOLVER'S OWN GATE AND MUST STAY THAT WAY.** ``own_composition_logvar`` tests
    ``tau > 1e-9``, so this partition answers "which mechanism did the solver USE here" — which is what
    the cross-tab and ``worst_objects.py`` need. It is deliberately **not** the question "should pass-0
    be SCORED here": a fitted κ that misses ½ by a rounding step leaves τ ~1e-7, which the solver does
    treat as evidence and which nonetheless has sd(λ) ~10³ nats and can resolve nothing. That second
    question needs a resolving-power floor and it lives in
    ``solvability_audit.own_evidence_tau_floor``. ⛔ Do not collapse the two — they gave different
    answers on the ladder's unstranded strata, and each is right for its own consumer.

    ``eps`` exists so a gate can move it and watch the partition move; production callers must not
    pass it. The three sources of ``tau_lam`` and their gating are ``node_init.build_node_init``.
    """
    tau = np.asarray(capture["_tau0_lam"], np.float64)
    struct_lock = g1_locked(capture["free_pos"], capture["free_neg"])
    relay_only = (tau <= eps) & (~struct_lock)
    return {
        "own_evidence": ~(struct_lock | relay_only),
        "relay_only": relay_only,
        "struct_lock": struct_lock,
    }


def _project(slot_mask, chain, n_nodes: int, n_edges: int) -> dict[str, np.ndarray]:
    """Scatter a per-slot boolean onto the node and edge axes. The chain is ``N E N E … N`` per
    reference, so every node and every contiguous edge is exactly one slot and the map is a
    bijection — there is nothing to pool and nothing to drop."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    mask = np.asarray(slot_mask, bool)
    out = {"node": np.zeros(n_nodes, bool), "edge": np.zeros(n_edges, bool)}
    out["node"][obj[kind == NODE]] = mask[kind == NODE]
    out["edge"][obj[kind == EDGE]] = mask[kind == EDGE]
    return out


def solver_class_masks(capture, chain, n_nodes: int, n_edges: int) -> dict[str, dict]:
    """:func:`solver_slot_classes`, projected onto the two scored axes."""
    slots = solver_slot_classes(capture, chain)
    return {
        axis: {name: _project(m, chain, n_nodes, n_edges)[axis] for name, m in slots.items()}
        for axis in AXES
    }


# ── class 2: is the answer THERE at all? (C_info) ────────────────────────────────────────────────


def info_class_masks(chain, region_arrays, substrate, gdna_pmf, rna_pmf) -> dict[str, dict]:
    """⭐ **C_info** — per object, is the two-component split recoverable from the stored channels?

    The 2×2 at one object is ``N = ρ_g·E_g + ρ_r·E_r`` and ``Σ1/L = ρ_g·D_g + ρ_r·D_r``, and it is
    identified iff ``E_g/D_g ≠ E_r/D_r`` — the two components' **opportunity-weighted** mean lengths.
    ⚠ At a contiguous edge the opportunity is ``(w−1)+`` and that reduces to the bare ``μ_g ≠ μ_r``;
    at a NODE it does **not**, because the opportunity is ``(ell − w + 1)+`` and ``1/w`` does not
    cancel it. Applying the edge form at a node would be scoring one frame's evidence against another
    frame's support, so this reads the moments in each slot's own frame.

    ⭐ **And it reads them from ``length_likelihood.build_slot_moments``, which already computes
    exactly this in exactly these frames.** Writing the algebra out again here would put two
    implementations of one quantity in the tree; the same argument is why ``LandedMoments`` carries
    ``eff`` at all. Conditional on the object's count ``N``, the deposited weight has mean
    ``pi·m1_g + (1−pi)·m1_r`` where ``pi`` is the gDNA share of the LANDED fragments — which is the
    quantity ``chain_node_deconv`` turns into mass — so::

        pi_hat = (Σ(1/L)/N − m1_r) / (m1_g − m1_r)

    Four classes, and the middle two are the answer this function exists to give:

    * ``absent`` — no count. There is no answer here to get right or wrong.
    * ``identified`` — either only ONE component has any opportunity here (a node too short for any
      RNA fragment to fit inside can only contain gDNA — determined, with no arithmetic at all), or
      both do, the moments separate, and ``pi_hat`` lands in ``[0, 1]``.
    * ``undet_no_separation`` — ``m1_g == m1_r`` **exactly**. At equal opportunity-weighted means the
      channel carries exactly zero information about composition, at any depth. ⚠ Tested exactly, not
      to a tolerance: a difference computed from large floats comes out flat only to ~1e-11, and a
      1e-11 row reads as live and then sells the grid's own width back as evidence. ⭐ The degenerate
      "no opportunity for EITHER component" case lands here without a special branch, because
      ``LandedMoments`` zeroes every moment at zero opportunity rather than flooring a division.
    * ``undet_out_of_range`` — the moments separate but the solution is outside ``[0, 1]``, so the
      observation is not consistent with ANY mixture of these two components. ⚠ This class also
      absorbs sampling noise: at one or two fragments a single draw of ``1/L`` easily falls outside
      the interval the two components span. That is why it is a separate class from the one above and
      not pooled into a single "undetermined" — the two mean different things.

    ⚠ This function ignores neighbours, which the sweep does not. It is a statement about one
    object's own channels and nothing else, and it must not be read as a bound on the solver.
    """
    mg = build_slot_moments(chain, region_arrays, gdna_pmf)
    mr = build_slot_moments(chain, region_arrays, rna_pmf)

    kind = np.asarray(chain.kind)
    n_slots = kind.shape[0]
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    count = np.zeros(n_slots, np.float64)
    inv = np.zeros(n_slots, np.float64)
    for k, view in ((NODE, substrate.node_contained), (EDGE, substrate.edge_unspliced)):
        sel = kind == k
        count[sel] = np.asarray(view.count, np.float64).sum(axis=1)[obj[sel]]
        inv[sel] = np.asarray(view.inv_length_sum, np.float64)[obj[sel]]

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
    n_nodes, n_edges = int(substrate.n_nodes), int(substrate.n_edges)
    return {
        axis: {name: _project(m, chain, n_nodes, n_edges)[axis] for name, m in slots.items()}
        for axis in AXES
    }


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def calibrate_arm(payload, kwargs, config, *, gdna_pmf=None, rna_pmf=None, debug=None):
    """One ``calibrate`` run. ``gdna_pmf`` / ``rna_pmf`` override the FITTED length models — that, and
    the solve depth in ``config``, are the only two things any arm below varies."""
    call = dict(kwargs)
    if gdna_pmf is not None:
        call["gdna_fl_pmf"] = gdna_pmf
    if rna_pmf is not None:
        call["rna_fl_pmf"] = rna_pmf
    return calibrate(payload=payload, config=config, _debug=debug, **call)


def load_or_build_oracle(bam, index, pipeline_config, work_dir, tag, full_payload, cache_root):
    """T, from a per-origin cache when one is valid — otherwise split, scan, and populate it.

    ⭐⭐ **WHY THIS IS WORTH CACHING AND THE FULL SCAN IS NOT.** The debug loop is "measure the whole
    panel → fix → measure again", and the fixes are in CALIBRATION. The oracle depends only on the
    accumulator and the index, so it is **invariant across every solver change the loop makes** — the
    cache is written once and hits for the rest of the campaign. Uncached, a 36-condition table costs
    hours of BAM splitting and re-scanning per iteration, which makes the loop impractical.

    ⛔ **KEYED BY THE SCAN CACHE'S OWN KEY, NOT A NEW ONE.** ``read_scan_cache`` refuses a payload
    whose ``graph_hash``, ``reach_digest``, ``payload_schema_digest`` or scan config does not describe
    the index it is loaded against — and ``reach`` in particular is covered by no other hash, so a
    rebuild that moved 38 % of contiguous reaches would verify clean against a home-made key. Reusing
    the shipped loader means a stale oracle is *refused*, loudly, rather than silently feeding
    everything downstream.

    ⚠ And the sum-to-full identity is re-run over the loaded arrays regardless
    (:meth:`OracleTruth.from_parts`): the cache can only skip the scanning, never the validation.
    """
    if cache_root is None:
        return OracleTruth.from_bam(bam, index, pipeline_config, Path(work_dir), tag,
                                    full_payload=full_payload)
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    dirs = {k: Path(cache_root) / tag / k for k in ORIGINS}
    try:
        parts = {k: read_scan_cache(dirs[k], index, scan).payload for k in ORIGINS}
        return OracleTruth.from_parts(full_payload, parts)
    except (FileNotFoundError, KeyError, ScanCacheKeyError):
        pass  # no cache, or it does not describe this index/scan — rebuild it below

    paths, read_counts = _split_bam(bam, Path(work_dir), tag)
    parts = {}
    for origin in ORIGINS:
        _stats, strand_model, _buf, payload = scan_and_buffer(paths[origin], index, scan)
        parts[origin] = payload
        write_scan_cache(dirs[origin], payload=payload, strand_model=strand_model, index=index,
                         bam=paths[origin], scan_config=scan)
    return OracleTruth.from_parts(full_payload, parts, read_counts)


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
    #: The two classifications as BOOLEAN MASKS per axis, kept so a downstream instrument reads the
    #: same partition this one scored rather than recomputing its own (two definitions of one class
    #: is how they drift). ``worst_objects.py`` consumes these.
    solver_masks: dict
    info_masks: dict
    cross: dict  #: axis -> (info class, solver class) -> AxisScore, for pass-0
    library_f_gdna: dict  #: "T" / arm name -> the LIBRARY-level gDNA fraction (the thermometer)
    #: ⭐ ``kind -> axis -> class -> objects IN the class``, scored or not. Distinct from
    #: ``AxisScore.n_scored``, which counts only objects with mass — and that distinction is the whole
    #: content of the ``absent`` class, whose scored count is 0 BY DEFINITION. Without this the class
    #: prints an empty row and the fact it exists to state ("most of the index carries no fragments")
    #: is invisible.
    class_objects: dict
    #: Which length pmfs C_info was classified with. ⚠ Recorded rather than assumed: the report says
    #: "with the simulator's own pmfs", and on a condition with no truth file it would be the FITTED
    #: ones — a caption describing something the code did not do.
    info_pmf_source: str
    seconds: float


def library_f_gdna(result) -> float:
    """The library gDNA fraction, summed over BOTH deconvolved axes — the reported deliverable.

    ⭐ Printed here beside the per-object answer for one reason: **the gap between the two is the
    finding.** The library figure is ``|Σ(g − t)|`` and the per-object answer is ``Σ|g − t|``, so a
    large under-call sitting next to a large over-call makes the first look an order of magnitude
    better than the second. Applying the same functional to T rather than to the simulator's origin
    counts keeps the comparison on the accumulator's own basis — the origin-count fraction is a
    different quantity (it counts fragments once; this counts each fragment on every object it
    touched), and subtracting one from the other would mix a basis change into the error.

    ⚠ Both axes, always: gDNA lives contained in a node AND crossing a line, and summing one axis
    reports a library's gDNA as a fraction of part of itself.
    """
    g = float(np.asarray(result.mass_gdna_node).sum() + np.asarray(result.mass_gdna_edge).sum())
    r = float(np.asarray(result.mass_rna_node).sum() + np.asarray(result.mass_rna_edge).sum())
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
    """Scan once, build T, run every arm, and score them per object and per class.

    ``truth_pmfs`` is a callable ``max_size -> (gdna_pmf, rna_pmf)``, or ``None`` to skip the C_input
    arms. ⚠ A callable rather than two arrays because a pmf must be sized by the PAYLOAD's own
    ``max_length`` — ``build_fl_models`` produces one of length ``max_length + 1``, and handing
    ``calibrate`` two pmfs of different lengths is a silent frame mismatch rather than an error. The
    payload does not exist until the scan below, and scanning a 1.6 GB BAM twice to learn one integer
    is four minutes per condition.

    ⚠ Four ``calibrate`` runs and four BAM scans. The oracle re-scans per partition and there is no
    cache for it; that is minutes, not seconds.
    """
    start = time.perf_counter()
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    # ⭐⭐ THE MAIN PAYLOAD IS CACHED TOO, and it is the same argument the oracle cache already makes:
    # the scan depends ONLY on the BAM, the index and the scan config — never on calibration — so one
    # cache serves every arm of a whole debugging campaign. Measured 8.3 s of a 24.5 s condition
    # (**34 %**), and it was being paid again for every arm of every A/B.
    # ⛔ Keyed by the SHIPPED loader, never a home-made key: ``read_scan_cache`` refuses a payload whose
    # ``graph_hash`` / ``reach_digest`` / ``payload_schema_digest`` / scan config does not describe this
    # index, and ``reach`` is covered by no other hash. A refusal here is loud and falls through to a
    # rescan; a home-made key would load a stale tally silently.
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

    # T. ⭐ Sum-to-full is validated on every bank exactly and RAISES if it does not hold — on the
    # cached path as well as the scanned one — so nothing below can run on an oracle that is not the
    # production payload split by origin.
    oracle = load_or_build_oracle(
        bam, index, pipeline_config, work_dir, tag, payload, oracle_cache
    )

    ra = RegionArrays.from_frame(index.nodes_df, index.ref_name_to_id)
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_edge_flags_array,
        build_junction_geometry_arrays,
    )

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        junction_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    truth_gdna_pmf, truth_rna_pmf = truth_pmfs(max_size) if truth_pmfs is not None else (None, None)
    kwargs = dict(
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        junctions=build_junction_geometry_arrays(index),
        edge_flags=build_edge_flags_array(index),
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

    # ⭐ THE CLASSES COME FROM PASS-0's OWN RUN, and are held fixed across every arm. ``tau_lam``
    # depends weakly on the incoming belief, so the final solve partitions slots slightly
    # differently — but "what evidence does this object have" is a property of the object, the
    # question Stage B asks is about the prior-free solve, and a class that moves between arms cannot
    # be used to compare them.
    chain = debug_pass0["chain"]
    n_nodes, n_edges = int(payload.n_nodes), int(payload.n_edges)
    solver_masks = solver_class_masks(debug_pass0["capture"], chain, n_nodes, n_edges)
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
            g = getattr(arm, f"mass_gdna_{axis}")
            r = getattr(arm, f"mass_rna_{axis}")
            tg = getattr(truth, f"mass_gdna_{axis}")
            tr = getattr(truth, f"mass_rna_{axis}")
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

    # ⭐ THE CROSS-TAB. Not a threshold on confidence — a threshold would be a magic number and this
    # needs none. "Undetermined by C_info × answered by the relay" is a cell of a partition, and its
    # error share is the statement.
    cross = {}
    p0 = arms["pass0"]
    for axis in AXES:
        g, r = getattr(p0, f"mass_gdna_{axis}"), getattr(p0, f"mass_rna_{axis}")
        tg, tr = getattr(truth, f"mass_gdna_{axis}"), getattr(truth, f"mass_rna_{axis}")
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
    """The simulator's own POST-CAPTURE length distribution for one origin class, as a pmf.

    ⚠ Post-capture empirical, not the configured ``frag_mean``: capture selects for length, so the
    configured parameters describe a library that was never sequenced. Same reader as
    ``calibration_truth_ab.py``'s ceiling arms, and it must stay the same reader — the whole value of
    a ceiling is that it is the consumer's own lever, not a second one that resembles it.
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
    """The library's TRUE gDNA fragment fraction, from the simulator's own origin counts.

    ⚠ From ``truth_summary.json``, never from the condition NAME: "gdna100" is a rate knob, not a
    fraction, and inferring one from the other is a rename away from silently wrong.
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
        "   ⛔ undrained on every arm, T included: three partitions cannot be drained independently."
    )
    for m in measurements:
        true_contained = m.scores["pass0"]["node"]["ALL"].mass
        true_crossing = m.scores["pass0"]["edge"]["ALL"].mass
        print()
        print(f"── {m.condition}    ({m.seconds:.0f} s)")
        print(
            f"   TRUE mass: node contained {_fmt(true_contained)}   "
            f"edge crossing {_fmt(true_crossing)}   "
            f"true gDNA node {_fmt(float(np.asarray(m.truth.mass_gdna_node).sum()))}"
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
    print("   own_evidence: strand / intron factory / length spoke here.  relay_only: nothing did — "
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
                # ⚠ An EMPTY class still prints its row. A report that only shows non-empty classes
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
                # ⚠ ``absent`` scores 0 objects and 0 mass BY DEFINITION — its content is the
                # "in class" column, which is why that column exists.
                print(
                    f"   {'':<5} {name:<22} {m.class_objects['info'][axis][name]:>10,} "
                    f"{s.n_scored:>9,} {mass_share:>10.1%} "
                    f"{s.abs_err:>14,.0f} {err_share:>9.1%} {s.mwae:>8.4f}"
                )


def _report_cross(measurements: list[ConditionMeasurement]) -> None:
    print()
    print("═" * 112)
    print("⭐⭐ THE CELL THIS SCRIPT EXISTS TO FIND — undetermined by C_info × answered by the RELAY")
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
                    star = "⭐" if (i in undetermined and s_name == "relay_only") else "  "
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
        # ⛔ **ONLY AN EXACTLY-ZERO ROW IS SATURATED**, and the distinction is not pedantry: at
        # truth = 0 any change that lowers the estimate scores better, which is the one-sidedness
        # that reversed a verdict here once. At truth = 0.01 it does not — a 1 % row is a perfectly
        # good scoring row and is where real libraries often sit.
        # ⚠ This test used to be ``truth <= 0.1``, which was correct for a panel whose only values
        # were 0 and 0.5 and became WRONG the moment the gDNA ladder added rungs at 1 %, 5 % and
        # 10 %: it would have silently dropped the entire low-gDNA end. A threshold that encodes one
        # panel's value set is a landmine for the next panel.
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
