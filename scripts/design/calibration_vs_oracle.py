"""IS THE SHIPPED CALIBRATION RIGHT, AND WHAT DOES ITS ERROR DO TO THE RULER THE EM DIVIDES BY?

⭐⭐⭐ **This is the 0.8.0 metric, scored one level ABOVE every other oracle instrument.** The
comparison is the calibration result against an ORACLE CALIBRATION — the same ``CalibrationResult``
object with only its deconvolved arrays replaced by the origin-split truth::

    P = calibrate(...)
    O = dataclasses.replace(P, **OracleTruth.from_parts(...).override_masses(region_arrays))

⛔ **AND IT IS THE ONLY INSTRUMENT THAT REACHES THE EFFECTIVE-LENGTH SHRINKAGE.** A calibration result
has two consumers and until this file existed every measurement arm in the tree patched the second one::

    calibrate(...)                          <- both consumers read THIS
      |-- transcript_capture_eff_lengths()  <- consumer A: `effective_lengths_em`, the EM's RULER
      |-- assemble_priors()                 <- consumer B: `LocusPriors`

``_setup_geometry_and_estimator`` builds consumer A **before** ``assemble_priors`` runs, so an arm that
wraps ``assemble_priors`` cannot see the ruler at all and every ceiling measured that way left the
shipped ruler installed. Substituting at the ``calibrate`` boundary reaches both.

⛔ **It deliberately does NOT re-score the prior.** ``prior_vs_oracle.py`` owns ``LocusPriors`` and a
second scorer is how a baseline and a ceiling drift apart. What is new here is the RULER.

⭐⭐ **THREE ARMS AND A NULL, and the null is what makes the ruler number attributable.**

=========  =====================================================================================
``P``      the shipped calibration
``O``      the oracle calibration — P with the six deconvolved arrays replaced by truth
``noop``   P with those six arrays replaced by THEMSELVES. ⛔ Must be BYTE-IDENTICAL to P, on the
           arrays *and* on the effective lengths derived from them; the gate runs before any table
``U``      ⭐ the NO-ENRICHMENT NULL: O's own total gDNA mass laid down at EXACTLY uniform density,
           ``m = rho_bar * support`` on both axes. Not a composition arm and never scored as one
=========  =====================================================================================

⛔⛔ **READ ``U`` IN TWO DIRECTIONS AND IT ANSWERS TWO DIFFERENT QUESTIONS.** Under uniform gDNA the
expected contained mass is exactly ``rho * gdna_region_eff_len``, so ``U`` is the same field as ``O``
with the SAMPLING NOISE removed and nothing else.

* At **capture-OFF** a uniform field is the physically correct configuration — there are no probes —
  so the shrinkage contract says the factor must be exactly **1.000**. Whatever ``U`` reads below 1 is
  contraction the estimator manufactured out of nothing, and it is reachable by no composition repair.
* At **capture-ON** ``U`` deliberately DESTROYS the real enrichment, so a factor near 1 there is the
  estimator working. It is the could-the-arm-have-fired check for the whole ruler column
  (TRAPS: could-the-arm-have-fired).

⭐ **Cheap enough to be the iteration loop: no solver, no EM, no BAM re-scan.** One ``calibrate`` off
the shipped scan cache per condition — measured 4-10 s, the whole 16-condition ladder in ~2 min. The
index and the region arrays load once, in 0.2 s.

⚠ **The four zero-gDNA conditions have no ``_main`` in the oracle cache and that is not a blocker.**
``_main`` is the UNDRAINED FULL PAYLOAD, which is the same scan as the plain scan cache, so it is read
from ``scan_cache/<condition>`` when absent — and which source was used is PRINTED per condition rather
than assumed, because a silent fallback is a truth source nobody checked.

⛔ **Score per stratum, never pooled**, and the 0.8.0 scope is stamped on every row: three strata are
the development target and unstranded x capture-ON is DEFERRED-but-REPORTED. A pooled total is a report
on the deferred stratum and nothing else (TRAPS: never-pool-the-strata).

Gates: ``tests/calibration/test_calibration_vs_oracle.py``. ``--self-test`` perturbs every comparator
with no I/O and is the same discipline in the other direction.

Usage::

    python scripts/design/calibration_vs_oracle.py                    # the whole ladder
    python scripts/design/calibration_vs_oracle.py --conditions NAME  # one condition
    python scripts/design/calibration_vs_oracle.py --jobs 4
    python scripts/design/calibration_vs_oracle.py --self-test        # no I/O
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


P0 = _sibling("pass0_vs_oracle.py")
PVO = _sibling("prior_vs_oracle.py")

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.capture_eff_length import (  # noqa: E402
    _global_reference_density,
    transcript_capture_eff_lengths,
)
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.frag_length_model import FragmentLengthModel  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402

#: ⭐ ONE HOME. The six deconvolved fields an oracle may substitute are ``prior_vs_oracle``'s list, not
#: a second copy — a set that drifted between two instruments would make their noop gates test
#: different things while both printed the word "identical".
OVERRIDE_FIELDS = PVO.OVERRIDE_FIELDS

DEFAULT_SUITE = PVO.DEFAULT_SUITE
DEFAULT_INDEX = PVO.DEFAULT_INDEX

#: The two axes ``CalibrationResult`` deconvolves; the sj axis is certified RNA and is never split.
AXES = P0.AXES


# ── the ruler: what the EM divides a transcript by ───────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class RulerScore:
    """One arm's transcript effective-length ruler. ⚠ ``factor`` is the OPPORTUNITY-WEIGHTED mean, not
    a mean of ratios: the per-transcript factors are weighted by the FL-marginal length they scale, so
    the aggregate is the ratio the EM's total opportunity actually moved by
    (TRAPS: a-mean-of-ratios-inherits-the-partition)."""

    rho_ref: float | None  #: the detected reference density; ``None`` = too little gDNA to detect one
    total_len: float  #: Σ eff_em over transcripts — the denominator the EM sums
    total_fl: float  #: Σ fl — the uncontracted FL-marginal length
    sum_ratio: float  #: Σ (eff_em / fl) — carried so the UNWEIGHTED mean is derivable
    n_transcripts: int

    @property
    def factor(self) -> float:
        """``Σ eff_em / Σ fl`` — 1.0 means no contraction. The aggregate the EM actually moves by."""
        return self.total_len / self.total_fl if self.total_fl > 0.0 else 1.0

    @property
    def factor_unweighted(self) -> float:
        """``mean(eff_em / fl)`` — every transcript counted once.

        ⚠ **Carried because the two disagree, and a reader comparing against a recorded number needs
        to know which one it was.** At the zero-gDNA control they read 0.095 and 0.345 on the same
        arrays: the contraction falls hardest on LONG transcripts, so weighting by the opportunity it
        scales makes it look far worse than counting transcripts equally does. ⛔ Neither is wrong —
        :attr:`factor` is the one that describes what the EM divides by, and it is the one to rank on.
        """
        return self.sum_ratio / self.n_transcripts if self.n_transcripts else 1.0


def ruler(calibration, region_arrays, index, fl_eff) -> tuple[RulerScore, np.ndarray]:
    """The shipped shrinkage, run on one arm. Returns its score and the per-transcript lengths.

    ⛔ **``transcript_capture_eff_lengths`` is called UNMODIFIED and that is the whole design.** The
    question this file asks is what a WRONG INPUT does to the shipped function, so re-deriving the
    contraction here would answer a different question and could be wrong in the same direction as the
    thing under test.
    """
    eff = transcript_capture_eff_lengths(calibration, region_arrays, index, fl_eff)
    return (
        RulerScore(
            rho_ref=_global_reference_density(
                np.asarray(calibration.mass_gdna_region, np.float64),
                np.maximum(np.asarray(calibration.gdna_region_eff_len, np.float64), 1e-9),
            ),
            total_len=float(np.asarray(eff, np.float64).sum()),
            total_fl=float(np.asarray(fl_eff, np.float64).sum()),
            sum_ratio=float(
                (np.asarray(eff, np.float64) / np.maximum(np.asarray(fl_eff, np.float64), 1e-9)).sum()
            ),
            n_transcripts=int(np.asarray(eff).size),
        ),
        np.asarray(eff, np.float64),
    )


#: ⭐ Signed per-object error buckets, in FRAGMENTS. Symmetric about an EXACT-ZERO bucket of its own,
#: because "this object is exactly right" and "this object is out by half a fragment" are different
#: statements and a histogram that merges them cannot show a solver that is right almost everywhere.
#: ⛔ Log-spaced, since per-object errors span five orders of magnitude on this panel; linear bins would
#: put 99 % of objects in the first bar and say nothing.
_HIST_EDGES = (1e-9, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


def signed_histogram(delta: np.ndarray) -> dict:
    """The full distribution of per-object signed error, zero at the centre.

    ``delta`` is ``gDNA_arm − gDNA_true`` per object, in fragments. Positive is an OVER-call — gDNA
    claimed where the truth was RNA — and negative is an under-call.

    ⭐ Returns counts AND the summed \\|error\\| carried by each bucket, because they answer different
    questions: a thousand objects out by one fragment and one object out by a thousand are the same
    ``Σ|err|`` and completely different defects.
    """
    d = np.asarray(delta, np.float64)
    edges = np.asarray(_HIST_EDGES, np.float64)
    out = {"n_objects": int(d.size), "buckets": []}
    labels: list[tuple[str, np.ndarray]] = []
    for i in range(edges.size - 1, 0, -1):
        labels.append((f"< -{edges[i - 1]:g}" if i == edges.size - 1 else
                       f"[-{edges[i]:g}, -{edges[i - 1]:g})",
                       (d <= -edges[i - 1]) & ((d > -edges[i]) if i < edges.size - 1 else True)))
    labels.append(("EXACT 0", np.abs(d) < edges[0]))
    for i in range(1, edges.size):
        labels.append((f"[{edges[i - 1]:g}, {edges[i]:g})" if i < edges.size - 1 else
                       f"> +{edges[i - 1]:g}",
                       (d >= edges[i - 1]) & ((d < edges[i]) if i < edges.size - 1 else True)))
    for name, sel in labels:
        out["buckets"].append(
            {"label": name, "n": int(sel.sum()), "abs_err": float(np.abs(d[sel]).sum())}
        )
    out["n_under"] = int((d < -edges[0]).sum())
    out["n_over"] = int((d > edges[0]).sum())
    out["n_exact"] = int((np.abs(d) <= edges[0]).sum())
    return out


#: ⭐ Buckets on the object's TRUE ``f_g``. The top one is the whole point: it is where the truth sits
#: at the simplex VERTEX, where ψ's Beta(½,½) reference is a fixed repulsion nothing swamps, and where
#: 57-86 % of this panel's calibration error lives. ⛔ Coarse away from the vertex on purpose — the
#: middle of the simplex is not where the defect is and a nine-column table hides the one column that
#: matters.
_VERTEX_EDGES = (0.0, 0.2, 0.8, 0.99, 0.999, 1.0000001)


def vertex_profile(p_arm, o_arm, axis: str) -> list[dict]:
    """Per-object error and composition CLOSURE, bucketed by the object's TRUE ``f_g``.

    ⭐⭐ **Two quantities per bucket, and they must be read together.**

    * ``mean_shortfall`` = ``mean(f_g_true − f_g_pred)``. Positive means the solver UNDER-calls gDNA.
      At the vertex bucket it is the repulsion, in units of composition.
    * ``mean_closure`` = ``mean(f_g + f_pos + f_neg)``. ⛔ **It is not a diagnostic of this arm — it is
      a CONFOUNDER for it**, and it is reported beside the shortfall so a reader cannot see one without
      the other.

    ⭐⭐⭐ **Why closure is a confounder rather than a curiosity.** ``f_g`` is published as a posterior
    MEDIAN and ``f_pos``/``f_neg`` as posterior MEANS, so for a single-strand object
    ``f_pos = 1 − E[f_g]`` and ``closure = 1 + (median − mean)`` exactly. **The composition closes iff
    the posterior is symmetric, and the deficit IS the skew.** Any change that moves posterior mass
    toward a vertex — which is what a composition prior does — moves the skew too, so a shortfall that
    improved and a closure that moved are not independent observations
    (TRAPS: one-thing-varied). Measured 2026-08-15: closure runs ~0.97-0.99 away from the vertex and
    **1.011-1.021 ABOVE 1 inside it**, i.e. it flips sign exactly where this arm operates.

    ⚠ Objects with no TRUE mass are excluded: they have no ``f_g`` to be right or wrong about.
    """
    gp = np.asarray(getattr(p_arm, f"mass_gdna_{axis}"), np.float64)
    rp = np.asarray(getattr(p_arm, f"mass_rna_{axis}"), np.float64)
    go = np.asarray(getattr(o_arm, f"mass_gdna_{axis}"), np.float64)
    ro = np.asarray(getattr(o_arm, f"mass_rna_{axis}"), np.float64)
    total = go + ro
    live = total > 0.0
    f_true = np.where(live, go / np.maximum(total, 1e-12), np.nan)
    f_pred = np.where(live, gp / np.maximum(gp + rp, 1e-12), np.nan)
    closure = (
        np.asarray(getattr(p_arm, f"gdna_frac_{axis}"), np.float64)
        + np.asarray(getattr(p_arm, f"rna_pos_frac_{axis}"), np.float64)
        + np.asarray(getattr(p_arm, f"rna_neg_frac_{axis}"), np.float64)
    )
    out = []
    for lo, hi in zip(_VERTEX_EDGES[:-1], _VERTEX_EDGES[1:], strict=True):
        sel = live & (f_true >= lo) & (f_true < hi)
        n = int(sel.sum())
        out.append(
            {
                "lo": float(lo),
                "hi": float(min(hi, 1.0)),
                "n": n,
                "mass": float(total[sel].sum()) if n else 0.0,
                "abs_err": float(np.abs(gp - go)[sel].sum()) if n else 0.0,
                "mean_shortfall": float(np.mean(f_true[sel] - f_pred[sel])) if n else 0.0,
                "mean_closure": float(np.mean(closure[sel])) if n else 0.0,
            }
        )
    return out


def pool_ledger(condition_dir: Path) -> dict:
    """The simulator's OWN starting fragment count per origin pool — the outer reference.

    ⛔ **Three pools on the truth side and TWO on the answer side, and that is structural rather than a
    limitation of this instrument.** ``calibrate`` deconvolves an object into ``(gDNA, RNA+, RNA−)`` and
    **cannot** split mature from nascent — that is the downstream EM's job, and
    ``quant_accuracy.py``'s pool table is where that split is scored. So ``nrna`` is reported here to
    keep the accounting honest and complete, and calibration's RNA answer is scored against
    ``mrna + nrna``.

    ⚠ From ``truth_summary.json``, never from the condition NAME.
    """
    summary = json.loads((Path(condition_dir) / "truth_summary.json").read_text())
    counts = summary["origin_counts"]
    return {k: float(counts.get(k, 0.0)) for k in ("gdna", "mrna", "nrna")}


def uniform_gdna_null(calibration):
    """⭐ ``U`` — the same gDNA TOTAL, laid down at exactly uniform density on both axes.

    Under uniform genomic gDNA at density ``rho`` the expected contained mass at an object is exactly
    ``rho * eff_len`` — the bedrock invariant the whole contraction rests on. So replacing each
    object's mass by ``rho_bar * eff_len``, with ``rho_bar = Σmass / Σeff_len``, removes the SAMPLING
    NOISE and changes nothing else: the library total is preserved exactly, per axis.

    ⛔ **Not a composition arm.** Its RNA arrays are untouched, so per-object conservation does not
    hold and it must never be handed to :func:`score_axis`. It exists only to answer "how much of the
    contraction survives when the field is noise-free", which is a question about the ESTIMATOR.
    """
    out = {}
    for axis in AXES:
        support = np.maximum(
            np.asarray(getattr(calibration, f"gdna_{axis}_eff_len"), np.float64), 1e-9
        )
        mass = np.asarray(getattr(calibration, f"mass_gdna_{axis}"), np.float64)
        total_support = float(support.sum())
        rho_bar = float(mass.sum()) / total_support if total_support > 0.0 else 0.0
        out[f"mass_gdna_{axis}"] = rho_bar * support
    return dataclasses.replace(calibration, **out)


# ── the gates, which run before any number is printed ────────────────────────────────────────────


def check_override_field_set(override: dict) -> None:
    """``override_masses`` must still write exactly :data:`OVERRIDE_FIELDS`.

    ⛔ Without this the ``noop`` arm would silently test a different set than the ``O`` arm — it would
    replace six fields with themselves while ``O`` replaced seven, and "byte-identical" would be a
    statement about the wrong six.
    """
    missing = set(OVERRIDE_FIELDS) - set(override)
    extra = set(override) - set(OVERRIDE_FIELDS)
    if missing or extra:
        raise SystemExit(
            f"⛔ override_masses no longer writes OVERRIDE_FIELDS: missing={sorted(missing)} "
            f"extra={sorted(extra)}. The noop gate would be testing a different set than the arm."
        )


def noop_differences(shipped, noop, eff_shipped, eff_noop) -> list[str]:
    """Every place the ``noop`` arm is not byte-identical to ``P``. Empty is the only pass.

    ⭐ **The effective lengths are compared too, and that is the point of doing it here.** Comparing
    only the six arrays would prove ``dataclasses.replace`` copies arrays, which nobody doubted. The
    claim under test is that the whole path from a substituted ``CalibrationResult`` down to the EM's
    ruler is inert when the substitution takes nothing — so the DERIVED quantity is what must match.

    ⚠ Byte-identity is reachable here and is NOT reachable for a ``quant_accuracy`` arm: there is no
    EM, no seed and no threaded scan in this path, so the two runs are the same arithmetic in the same
    order (TRAPS: the-deliverable-is-not-reproducible-by-default).
    """
    bad = [
        f for f in OVERRIDE_FIELDS
        if not np.array_equal(np.asarray(getattr(shipped, f)), np.asarray(getattr(noop, f)))
    ]
    if not np.array_equal(eff_shipped, eff_noop):
        worst = float(np.abs(eff_shipped - eff_noop).max())
        bad.append(f"effective_lengths_em (max|delta|={worst:.3e})")
    return bad


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def load_oracle(suite: Path, oracle_cache: Path, condition: str, index, drained_payload, lift):
    """The origin-split truth for one condition, **in the DRAINED frame** (the frame ruling of
    2026-08-31, `ISSUES: instruments-calibrate-undrained-cache`).

    ⭐ The FULL side is the very payload ``P`` calibrated — ``_main`` is byte-identical to the plain
    scan cache (both store pass one), so draining either yields the same frame and re-reading
    ``_main`` would be a second copy of the same quantity; the zero-gDNA rows, which have no
    ``_main``, therefore need no announced fallback any more. The cached PARTS are drained by
    replaying the whole's already-drawn choices (`from_cached_parts` → `lift_drain_parts`), and
    ``from_parts``' sum-to-full then validates the lift end to end on the drained frame.
    ⚠ The oracle may carry a nonzero ``gdna_spliced_leak`` — production's own drain behaviour
    (`ISSUES: drain-contaminates-certified-rna`) — and ``n_ambiguous`` bounds the lift's origin
    attribution; both are REPORTED on the row, never swallowed.
    """
    root = Path(oracle_cache) / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    return OracleTruth.from_cached_parts(drained_payload, parts, lift), "cached parts, drained frame"


def measure_condition(index, region_arrays, pipeline_config, suite: Path, oracle_cache: Path,
                      condition: str) -> dict:
    """calibrate once, build P / O / noop / U, gate them, and score. One JSON-able row."""
    start = time.perf_counter()
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    # ⭐ the DRAINED frame — everything below (P, the substrate, the oracle's full side) reads THIS
    # payload, never `cache.payload`, or the same-basis gates would be comparing two frames.
    payload = kw["payload"]
    p_arm = calibrate(config=pipeline_config.calibration, **kw)

    oracle, oracle_source = load_oracle(suite, oracle_cache, condition, index, payload, lift)
    override = oracle.override_masses(region_arrays)
    check_override_field_set(override)
    o_arm = dataclasses.replace(p_arm, **override)
    noop_arm = dataclasses.replace(p_arm, **{f: getattr(p_arm, f) for f in OVERRIDE_FIELDS})

    # ⛔ Both arms must be on the payload's own per-object totals, per axis. Without that identity a
    # mass-weighted mean of fractions is an average over different denominators.
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    P0.check_same_basis("P", p_arm, substrate)
    P0.check_same_basis("O", o_arm, substrate)

    # -- the ruler, on the FL-marginal lengths the pipeline builds it from --
    rna_fl = FragmentLengthModel.from_pmf(kw["rna_fl_pmf"], int(payload.max_length))
    fl_eff = rna_fl.compute_all_transcript_eff_lens(
        index.t_df["length"].values.astype(np.int64)
    )
    rulers, lengths = {}, {}
    for name, arm in (("P", p_arm), ("O", o_arm), ("noop", noop_arm),
                      ("U", uniform_gdna_null(o_arm))):
        rulers[name], lengths[name] = ruler(arm, region_arrays, index, fl_eff)

    bad = noop_differences(p_arm, noop_arm, lengths["P"], lengths["noop"])

    row = {
        "condition": condition,
        "stratum": list(PVO.stratum(condition)),
        "oracle_source": oracle_source,
        "gdna_spliced_leak": oracle.gdna_spliced_leak,
        "lift_n_ambiguous": oracle.n_ambiguous,
        "noop_differences": bad,
        "seconds": time.perf_counter() - start,
        "library_f_gdna_P": P0.library_f_gdna(p_arm),
        "library_f_gdna_O": P0.library_f_gdna(o_arm),
        "axes": {},
        "ruler": {k: dataclasses.asdict(v) for k, v in rulers.items()},
        # ⭐ Σ|Δ| over the ruler itself, in base pairs of opportunity — the quantity the EM divides by.
        "ruler_abs_err": float(np.abs(lengths["P"] - lengths["O"]).sum()),
        "ruler_n_moved": int(np.sum(lengths["P"] != lengths["O"])),
    }
    for axis in AXES:
        s = P0.score_axis(
            getattr(p_arm, f"mass_gdna_{axis}"), getattr(p_arm, f"mass_rna_{axis}"),
            getattr(o_arm, f"mass_gdna_{axis}"), getattr(o_arm, f"mass_rna_{axis}"),
        )
        row["axes"][axis] = dataclasses.asdict(s)
        # ⭐ THE FULL DISTRIBUTION, not just its sum. Objects with no mass at all are excluded: they
        # have no answer to get right, and folding them in puts 80 % of a genome in the EXACT-0 bar.
        gp = np.asarray(getattr(p_arm, f"mass_gdna_{axis}"), np.float64)
        go = np.asarray(getattr(o_arm, f"mass_gdna_{axis}"), np.float64)
        total = go + np.asarray(getattr(o_arm, f"mass_rna_{axis}"), np.float64)
        live = total > 0.0
        row["axes"][axis]["hist"] = signed_histogram(gp[live] - go[live])
        # ⛔ Δ_RNA ≡ −Δ_gDNA per object, because `check_same_basis` has just established that the two
        # arms carry the SAME per-object total. Recorded as a gate rather than left implicit: two
        # columns that are the same number would read as two independent measurements.
        rp = np.asarray(getattr(p_arm, f"mass_rna_{axis}"), np.float64)
        ro = np.asarray(getattr(o_arm, f"mass_rna_{axis}"), np.float64)
        row["axes"][axis]["rna_mirrors_gdna_max_dev"] = float(
            np.abs((rp - ro)[live] + (gp - go)[live]).max() if live.any() else 0.0
        )
        row["axes"][axis]["vertex"] = vertex_profile(p_arm, o_arm, axis)

    # ── the POOL LEDGER, in fragments ──
    truth = pool_ledger(Path(suite) / condition)
    row["pools"] = {
        "true_gdna": truth["gdna"],
        "true_mrna": truth["mrna"],
        "true_nrna": truth["nrna"],
        "true_rna": truth["mrna"] + truth["nrna"],
        # ⭐ the CONSERVED fragment counts — each axis converted by its own population's mass per
        # crossing, so these are fragments and not object incidences.
        "P_gdna": p_arm.library_gdna_fragments,
        "P_rna": p_arm.library_rna_fragments,
        "O_gdna": o_arm.library_gdna_fragments,
        "O_rna": o_arm.library_rna_fragments,
    }
    return row


# ── reporting ────────────────────────────────────────────────────────────────────────────────────

#: ⭐⭐ The 0.8.0 scope, stamped on the row rather than left to the reader. Neither
#: ``solvability_audit.py`` nor ``prior_vs_oracle.py`` marks it, so a reader ranking on a stratum table
#: had to remember which cell was the development target — and the deferred one carries most of the
#: error, so forgetting inverts the ranking. ⛔ DEFERRED IS REPORTED, never dropped.
_SCOPE = {
    ("stranded", "capture OFF"): "IN SCOPE",
    ("stranded", "capture ON"): "IN SCOPE",
    ("unstranded", "capture OFF"): "IN SCOPE",
    ("unstranded", "capture ON"): "DEFERRED",
}

#: Every selection every table prints, in order — one list, so a stratum cannot appear on some tables
#: and not others. ``None`` is a rule boundary.
_SELECTIONS = (
    *(
        (f"{' x '.join(st)}  [{_SCOPE[st]}]", (lambda c, st=st: PVO.stratum(c) == tuple(st)
                                               and not PVO.is_zero_gdna(c)))
        for st in _SCOPE
    ),
    (None, None),
    ("⛔ g00 ZERO-gDNA control (all strata)", PVO.is_zero_gdna),
)


def _agg_axis(scores: list[dict]) -> dict | None:
    """Sum a list of ``AxisScore`` dicts into one.

    ⛔ **The rate is RE-DERIVED from the summed totals, never averaged.** A mean of per-condition
    ``mwae`` over conditions of different depth is a number with no consumer
    (TRAPS: never-pool-the-strata). Every other field of ``AxisScore`` is a mass and adds.
    """
    scores = [s for s in scores if s is not None]
    if not scores:
        return None
    out = {
        k: sum(s[k] for s in scores)
        for k in ("n_scored", "mass", "net_err", "abs_err", "over_call", "under_call")
    }
    out["mwae"] = out["abs_err"] / out["mass"] if out["mass"] > 0.0 else 0.0
    return out


def _fmt_rho(x) -> str:
    return f"{'None':>10}" if x is None else f"{x:>10.4g}"


def report(rows: list[dict]) -> None:
    """The whole report from the per-condition JSON — **the only report path there is.**

    ⭐ It reads the serialised rows in the serial case too, so ``--jobs 1`` and ``--jobs 4`` print
    numbers produced by one code path rather than two that nothing compares.
    """
    # ── the gates first. A table read before its gate is a table nobody checked. ──
    print()
    print("=" * 118)
    print("  ⭐⭐⭐ CALIBRATION vs ORACLE CALIBRATION — the 0.8.0 metric, and the EM's RULER above it")
    print(f"  {len(rows)} conditions   no solver, no EM, no BAM re-scan")
    print("=" * 118)
    print()
    failed = [r for r in rows if r["noop_differences"]]
    if failed:
        for r in failed:
            print(f"  ⛔ {r['condition']}: noop is NOT byte-identical to P -> {r['noop_differences']}")
        raise SystemExit(2)
    print(f"  ✅ GATE  noop is byte-identical to P on all {len(OVERRIDE_FIELDS)} override fields AND on")
    print(f"           the effective lengths derived from them, on {len(rows)}/{len(rows)} conditions")
    print("  ✅ GATE  override_masses writes exactly the override field set")
    print("  ✅ GATE  P and O are on the payload's own per-object totals, both axes (check_same_basis)")
    fallbacks = sorted(r["condition"] for r in rows if r["oracle_source"] != "_main")
    if fallbacks:
        print(f"  ⚠  {len(fallbacks)} condition(s) read the full payload from the SCAN CACHE because the")
        print("           oracle cache holds no `_main` for a zero-gDNA row; sum-to-full passed on each:")
        for c in fallbacks:
            print(f"             {c}")

    def sel_rows(pred):
        return [r for r in rows if pred(r["condition"])]

    worst_mirror = max((r["axes"][a]["rna_mirrors_gdna_max_dev"] for r in rows for a in AXES),
                       default=0.0)
    print(f"  ✅ GATE  Δ_RNA ≡ −Δ_gDNA per object on both axes (max deviation {worst_mirror:.3e}), so")
    print("           ONE signed number per object describes the whole error — reported once, not twice")

    # ── ⓪ THE POOL LEDGER, in fragments ──
    print()
    print("  ⓪ ⭐⭐⭐ THE POOL LEDGER — FRAGMENT COUNTS, not ratios. What went in, and what calibration says.")
    print("     ⛔ THREE pools on the truth side, TWO on the answer side: `calibrate` splits gDNA from RNA")
    print("        and CANNOT split mature from nascent — that is the EM's job (`quant_accuracy.py`).")
    print(f"    {'condition / stratum':<44} {'true gDNA':>12} {'true mRNA':>12} {'true nRNA':>10} "
          f"{'P gDNA':>13} {'P RNA':>13} {'ΔgDNA vs O':>13} {'P/O gDNA':>9}")
    print("    " + "-" * 142)
    for title, pred in _SELECTIONS:
        if title is None:
            print("    " + "-" * 142)
            continue
        sub = sel_rows(pred)
        if not sub:
            continue
        t = {k: sum(r["pools"][k] for r in sub) for k in
             ("true_gdna", "true_mrna", "true_nrna", "P_gdna", "P_rna", "O_gdna")}
        d = t["P_gdna"] - t["O_gdna"]
        ratio = t["P_gdna"] / t["O_gdna"] if t["O_gdna"] > 0 else float("nan")
        print(f"    {title:<44} {t['true_gdna']:>12,.0f} {t['true_mrna']:>12,.0f} "
              f"{t['true_nrna']:>10,.0f} {t['P_gdna']:>13,.0f} {t['P_rna']:>13,.0f} "
              f"{d:>+13,.0f} {ratio:>9.4f}")
    if all(r["pools"]["true_nrna"] == 0.0 for r in rows):
        print("    ⛔ EVERY condition has true nRNA = 0 — this panel cannot resolve the nascent axis at")
        print("       all. A gDNA-vs-nascent-vs-annotated verdict needs nascent-bearing conditions.")
    print("    ⚠ `P/O gDNA` is the like-for-like comparison. The simulator's own count is the OUTER")
    print("      reference and differs from O by the accumulator's retention, which is a different axis.")

    # ── ① the composition, per axis ──
    for axis in AXES:
        print()
        print(f"  ① COMPOSITION on the {axis.upper()} axis — P against O, mass-weighted")
        print(f"    {'stratum':<38} {'objects':>9} {'mass':>14} {'Σ|Δ gDNA|':>14} "
              f"{'mwae':>8} {'net':>14} {'over':>13} {'under':>13}")
        print("    " + "-" * 130)
        for title, pred in _SELECTIONS:
            if title is None:
                print("    " + "-" * 130)
                continue
            a = _agg_axis([r["axes"][axis] for r in sel_rows(pred)])
            if a is None:
                continue
            print(f"    {title:<38} {a['n_scored']:>9,} {a['mass']:>14,.0f} {a['abs_err']:>14,.0f} "
                  f"{a['mwae']:>8.4f} {a['net_err']:>+14,.0f} {a['over_call']:>13,.0f} "
                  f"{a['under_call']:>13,.0f}")

    # ── ② the library figure ──
    print()
    print("  ② THE LIBRARY gDNA FRACTION — the one-number thermometer, per stratum")
    print(f"    {'stratum':<38} {'f_gdna P':>10} {'f_gdna O':>10} {'|Δ|':>10}")
    print("    " + "-" * 72)
    for title, pred in _SELECTIONS:
        if title is None:
            print("    " + "-" * 72)
            continue
        sub = sel_rows(pred)
        if not sub:
            continue
        # ⛔ A ratio of sums, not a mean of ratios — the conditions differ in depth.
        p = sum(r["library_f_gdna_P"] * r["axes"]["region"]["mass"] for r in sub)
        o = sum(r["library_f_gdna_O"] * r["axes"]["region"]["mass"] for r in sub)
        w = sum(r["axes"]["region"]["mass"] for r in sub)
        p, o = (p / w, o / w) if w > 0 else (float("nan"), float("nan"))
        print(f"    {title:<38} {p:>10.4f} {o:>10.4f} {abs(p - o):>10.4f}")

    # ── ③ THE RULER ──
    print()
    print("  ③ ⭐⭐⭐ THE RULER — `effective_lengths_em`, the transcript length the EM DIVIDES BY.")
    print("     No other instrument reaches this: it is built BEFORE `assemble_priors`, which is what")
    print("     every other arm patches. `factor` is Σ eff_em / Σ fl; 1.000 means no contraction.")
    print("     ⛔ At capture-OFF the contract says the factor is EXACTLY 1.000 — there are no probes.")
    print(f"    {'stratum':<38} {'factor P':>9} {'factor O':>9} {'factor U':>9} "
          f"{'P/O':>8} {'Σ|Δ len|':>15} {'moved':>9}")
    print("    " + "-" * 120)
    for title, pred in _SELECTIONS:
        if title is None:
            print("    " + "-" * 120)
            continue
        sub = sel_rows(pred)
        if not sub:
            continue
        fac = {}
        for arm in ("P", "O", "U"):
            tl = sum(r["ruler"][arm]["total_len"] for r in sub)
            tf = sum(r["ruler"][arm]["total_fl"] for r in sub)
            fac[arm] = tl / tf if tf > 0 else float("nan")
        print(f"    {title:<38} {fac['P']:>9.4f} {fac['O']:>9.4f} {fac['U']:>9.4f} "
              f"{fac['P'] / fac['O'] if fac['O'] else float('nan'):>8.3f} "
              f"{sum(r['ruler_abs_err'] for r in sub):>15,.0f} "
              f"{sum(r['ruler_n_moved'] for r in sub):>9,}")

    # ── ⑤ THE FULL SIGNED DISTRIBUTION ──
    for axis in AXES:
        print()
        print(f"  ⑤ ⭐⭐⭐ THE SIGNED ERROR DISTRIBUTION on the {axis.upper()} axis — zero IS the truth.")
        print("     Δ = gDNA_shipped − gDNA_true, in FRAGMENTS, per object. + is an OVER-call (gDNA")
        print("     claimed where the truth was RNA). Objects carrying no mass are excluded.")
        labels = [b["label"] for b in rows[0]["axes"][axis]["hist"]["buckets"]]
        print(f"    {'stratum':<38} " + " ".join(f"{lab:>13}" for lab in labels))
        print("    " + "-" * (38 + 14 * len(labels)))
        for title, pred in _SELECTIONS:
            if title is None:
                print("    " + "-" * (38 + 14 * len(labels)))
                continue
            sub = sel_rows(pred)
            if not sub:
                continue
            tot = [sum(r["axes"][axis]["hist"]["buckets"][i]["n"] for r in sub)
                   for i in range(len(labels))]
            print(f"    {title:<38} " + " ".join(f"{n:>13,}" for n in tot))
        print()
        print(f"    {'stratum':<38} {'under':>13} {'EXACT':>13} {'over':>13} {'objects':>13} "
              f"{'Σ|Δ| frags':>15}")
        print("    " + "-" * 110)
        for title, pred in _SELECTIONS:
            if title is None:
                print("    " + "-" * 110)
                continue
            sub = sel_rows(pred)
            if not sub:
                continue
            h = [r["axes"][axis]["hist"] for r in sub]
            print(f"    {title:<38} {sum(x['n_under'] for x in h):>13,} "
                  f"{sum(x['n_exact'] for x in h):>13,} {sum(x['n_over'] for x in h):>13,} "
                  f"{sum(x['n_objects'] for x in h):>13,} "
                  f"{sum(r['axes'][axis]['abs_err'] for r in sub):>15,.0f}")

    # ── ⑦ THE VERTEX PROFILE ──
    for axis in AXES:
        print()
        print(f"  ⑦ ⭐⭐⭐ WHERE THE ERROR SITS ON THE SIMPLEX — {axis.upper()} axis, bucketed by TRUE f_g.")
        print("     `shortfall` = mean(true − pred): + means the solver UNDER-calls gDNA.")
        print("     ⛔ `closure` = mean(f_g + f_pos + f_neg). It is a CONFOUNDER, not a diagnostic:")
        print("        f_g is a posterior MEDIAN and f_pos/f_neg are MEANS, so closure = 1 + (median −")
        print("        mean) exactly — the deficit IS the posterior skew, and a composition prior moves")
        print("        skew by construction. Read it beside the shortfall, never instead of it.")
        print(f"    {'stratum':<38} {'true f_g':<15} {'objects':>9} {'mass':>13} "
              f"{'shortfall':>10} {'Σ|Δ| frags':>12} {'of Σ|Δ|':>8} {'closure':>8}")
        print("    " + "-" * 128)
        for title, pred in _SELECTIONS:
            if title is None:
                print("    " + "-" * 128)
                continue
            sub = sel_rows(pred)
            if not sub:
                continue
            n_b = len(sub[0]["axes"][axis]["vertex"])
            grand = sum(r["axes"][axis]["abs_err"] for r in sub)
            for b in range(n_b):
                cells = [r["axes"][axis]["vertex"][b] for r in sub]
                n = sum(c["n"] for c in cells)
                if n == 0:
                    continue
                mass = sum(c["mass"] for c in cells)
                ae = sum(c["abs_err"] for c in cells)
                # ⛔ weighted by objects, never a mean of per-condition means
                sf = sum(c["mean_shortfall"] * c["n"] for c in cells) / n
                cl = sum(c["mean_closure"] * c["n"] for c in cells) / n
                lab = f"[{cells[0]['lo']:g}, {cells[0]['hi']:g}]"
                print(f"    {title if b == 0 else '':<38} {lab:<15} {n:>9,} {mass:>13,.0f} "
                      f"{sf:>+10.4f} {ae:>12,.0f} {ae / max(grand, 1e-9):>7.1%} {cl:>8.4f}")

    # ── ⑥ SCENARIOS RANKED WORST FIRST ──
    print()
    print("  ⑥ ⭐⭐⭐ SCENARIOS RANKED BY TOTAL Σ|Δ| IN FRAGMENTS (region + boundary) — worst first.")
    print("     ⛔ The scope tag is the ranking. The worst row overall is the DEFERRED stratum every")
    print("        time; take the worst IN-SCOPE row instead.")
    print(f"    {'#':>3} {'condition':<44} {'scope':<9} {'Σ|Δ| region':>13} {'Σ|Δ| bnd':>13} "
          f"{'Σ|Δ| TOTAL':>13} {'over':>12} {'under':>12}")
    print("    " + "-" * 126)
    ranked = sorted(
        rows, key=lambda r: -(r["axes"]["region"]["abs_err"] + r["axes"]["boundary"]["abs_err"])
    )
    for i, r in enumerate(ranked, 1):
        reg, bnd = r["axes"]["region"], r["axes"]["boundary"]
        scope = _SCOPE[tuple(r["stratum"])]
        print(f"    {i:>3} {r['condition']:<44} {scope:<9} {reg['abs_err']:>13,.0f} "
              f"{bnd['abs_err']:>13,.0f} {reg['abs_err'] + bnd['abs_err']:>13,.0f} "
              f"{reg['over_call'] + bnd['over_call']:>12,.0f} "
              f"{reg['under_call'] + bnd['under_call']:>12,.0f}")
    in_scope = [r for r in ranked if _SCOPE[tuple(r["stratum"])] == "IN SCOPE"]
    if in_scope:
        w = in_scope[0]
        print(f"    ⭐ WORST IN-SCOPE: {w['condition']}  "
              f"Σ|Δ| = {w['axes']['region']['abs_err'] + w['axes']['boundary']['abs_err']:,.0f} fragments")

    # ── ④ per condition ──
    print()
    print("  ④ PER CONDITION — rank within a stratum, never across one")
    print(f"    {'condition':<44} {'mwae reg':>9} {'mwae bnd':>9} {'lib P':>8} {'lib O':>8} "
          f"{'rho_ref P':>10} {'rho_ref O':>10} {'fac P':>7} {'fac O':>7} {'fac U':>7} {'s':>6}")
    print("    " + "-" * 146)
    for r in sorted(rows, key=lambda x: (tuple(x["stratum"]), x["condition"])):
        rl = r["ruler"]
        print(f"    {r['condition']:<44} {r['axes']['region']['mwae']:>9.4f} "
              f"{r['axes']['boundary']['mwae']:>9.4f} {r['library_f_gdna_P']:>8.4f} "
              f"{r['library_f_gdna_O']:>8.4f} {_fmt_rho(rl['P']['rho_ref'])} "
              f"{_fmt_rho(rl['O']['rho_ref'])} "
              f"{rl['P']['total_len'] / rl['P']['total_fl']:>7.4f} "
              f"{rl['O']['total_len'] / rl['O']['total_fl']:>7.4f} "
              f"{rl['U']['total_len'] / rl['U']['total_fl']:>7.4f} {r['seconds']:>6.1f}")
    print()
    print(f"  total wall clock {sum(r['seconds'] for r in rows):.1f} s over {len(rows)} conditions")


# ── the falsification: every comparator perturbed, with no I/O ───────────────────────────────────


def _toy_calibration(n_regions: int = 24, n_boundaries: int = 20, n_sj: int = 4, seed: int = 7):
    """A minimal well-formed ``CalibrationResult`` for the self-test. No index, no payload, no I/O."""
    from rigel.calibration.result import CalibrationResult
    from rigel.config import CalibrationConfig

    rng = np.random.default_rng(seed)
    ones_r, ones_b = np.ones(n_regions), np.ones(n_boundaries)
    gr = rng.uniform(1.0, 9.0, n_regions)
    gb = rng.uniform(1.0, 9.0, n_boundaries)
    return CalibrationResult(
        mass_gdna_region=gr,
        mass_rna_region=rng.uniform(1.0, 9.0, n_regions),
        mass_gdna_boundary=gb,
        mass_rna_boundary=rng.uniform(1.0, 9.0, n_boundaries),
        mass_rna_spliced_boundary=np.zeros(n_boundaries),
        boundary_mass_per_crossing=ones_b.copy(),
        count_rna_sj=np.zeros(n_sj),
        boundary_spliced_mass_per_crossing=ones_b.copy(),
        sj_mass_per_crossing=np.ones(n_sj),
        gdna_region_eff_len=rng.uniform(50.0, 500.0, n_regions),
        gdna_boundary_eff_len=np.full(n_boundaries, 120.0),
        rna_region_eff_len=ones_r.copy(),
        rna_boundary_eff_len=ones_b.copy(),
        gdna_frac_region=np.full(n_regions, 0.5),
        rna_pos_frac_region=np.full(n_regions, 0.25),
        rna_neg_frac_region=np.full(n_regions, 0.25),
        gdna_frac_boundary=np.full(n_boundaries, 0.5),
        rna_pos_frac_boundary=np.full(n_boundaries, 0.25),
        rna_neg_frac_boundary=np.full(n_boundaries, 0.25),
        gdna_density_global=0.1,
        rna_sense_frac=0.5,
        gdna_strand_overdispersion=0.0,
        rna_strand_overdispersion=0.0,
        n_regions=n_regions,
        n_boundaries=n_boundaries,
        n_sj=n_sj,
        config=CalibrationConfig(),
    )


def self_test() -> int:
    """⭐ Perturb every comparator and require it to FIRE. A green read proves nothing on its own —
    writing the gate is half the discipline and breaking the fixed code is the other half
    (TRAPS: perturb-every-gate)."""
    checks: list[tuple[str, bool]] = []

    def check(name: str, ok: bool):
        checks.append((name, bool(ok)))

    cal = _toy_calibration()

    # ① score_axis against ITSELF is exactly zero, and a ONE-ULP nudge makes it nonzero.
    s = P0.score_axis(cal.mass_gdna_region, cal.mass_rna_region,
                      cal.mass_gdna_region, cal.mass_rna_region)
    check("score_axis(P, P) is exactly 0", s.abs_err == 0.0 and s.mwae == 0.0)
    nudged = np.array(cal.mass_gdna_region, copy=True)
    nudged[3] = np.nextafter(nudged[3], np.inf)
    # ⚠ the RNA side moves the opposite way so the per-object TOTAL is preserved: score_axis refuses
    # two arms on different bases, and a basis refusal is not the perturbation under test.
    rna_n = np.array(cal.mass_rna_region, copy=True)
    rna_n[3] = np.nextafter(rna_n[3], -np.inf)
    s1 = P0.score_axis(nudged, rna_n, cal.mass_gdna_region, cal.mass_rna_region)
    check("score_axis resolves a ONE-ULP nudge", s1.abs_err > 0.0)

    # ② the noop comparator fires on a one-ULP nudge to an override array, and on the lengths alone.
    eff = np.linspace(100.0, 900.0, 16)
    check("noop comparator is clean on identical input",
          noop_differences(cal, cal, eff, eff) == [])
    bad_cal = dataclasses.replace(cal, mass_gdna_region=nudged)
    check("noop comparator resolves a ONE-ULP array nudge",
          noop_differences(cal, bad_cal, eff, eff) == ["mass_gdna_region"])
    eff2 = np.array(eff, copy=True)
    eff2[5] = np.nextafter(eff2[5], np.inf)
    check("noop comparator resolves a ONE-ULP LENGTH nudge, arrays identical",
          len(noop_differences(cal, cal, eff, eff2)) == 1)

    # ③ the override field-set gate refuses a dropped field and an extra one.
    full = {f: getattr(cal, f) for f in OVERRIDE_FIELDS}
    try:
        check_override_field_set(full)
        ok_full = True
    except SystemExit:
        ok_full = False
    check("field-set gate passes the real set", ok_full)
    for label, mutated in (("a DROPPED field", {k: v for k, v in list(full.items())[1:]}),
                           ("an EXTRA field", {**full, "gdna_density_global": 0.0})):
        try:
            check_override_field_set(mutated)
            fired = False
        except SystemExit:
            fired = True
        check(f"field-set gate refuses {label}", fired)

    # ④ ⭐ THE BEDROCK INVARIANT the U null rests on: an EXACTLY uniform field returns its own
    #    density, so `min(rho/rho_ref, 1)` is 1 everywhere and the contraction is exactly none.
    support = np.asarray(cal.gdna_region_eff_len, np.float64)
    rho_true = 0.037
    rr = _global_reference_density(rho_true * support, support)
    check("uniform field: rho_ref recovers its own density EXACTLY", rr == rho_true)
    # ⚠ Stated as "no contraction survives", not as bit-equality of the ratio: the realised density is
    # `(rho*S)/S`, and a multiply-then-divide does not round-trip, so demanding an exact 1.0 would test
    # float associativity rather than the invariant. The load-bearing claim is that `min(rho/rho_ref, 1)`
    # leaves nothing behind — measured below at 1e-16, against the 0.36 the real panel shows.
    weights = np.minimum((rho_true * support / support) / rr, 1.0)
    check("uniform field: the contraction it leaves is at float noise",
          float(1.0 - weights.min()) < 1e-12)
    # ⛔ and NOT vacuous: an object genuinely below the reference must still be contracted.
    depleted = np.array(rho_true * support, copy=True)
    depleted[0] *= 0.25
    check("a genuinely depleted object IS contracted",
          np.isclose(min((depleted[0] / support[0]) / rr, 1.0), 0.25, rtol=1e-9))
    # and it must NOT be blind: a bimodal field must return the ENRICHED mode, not the depleted one.
    bimodal = rho_true * support
    hot = np.arange(0, support.size, 3)
    bimodal[hot] *= 100.0
    rr_bi = _global_reference_density(bimodal, support)
    check("bimodal field: rho_ref finds the ENRICHED mode", rr_bi is not None and rr_bi > 10 * rho_true)

    # ⑤ the U null preserves the gDNA TOTAL per axis and flattens the field.
    u = uniform_gdna_null(cal)
    for axis in AXES:
        before = float(np.asarray(getattr(cal, f"mass_gdna_{axis}")).sum())
        after = float(np.asarray(getattr(u, f"mass_gdna_{axis}")).sum())
        check(f"U preserves the gDNA total on the {axis} axis", np.isclose(before, after, rtol=1e-12))
    rho_u = np.asarray(u.mass_gdna_region) / np.asarray(u.gdna_region_eff_len)
    check("U flattens the density field to a constant", np.allclose(rho_u, rho_u[0], rtol=1e-12))
    check("U leaves the RNA arrays untouched",
          np.array_equal(np.asarray(u.mass_rna_region), np.asarray(cal.mass_rna_region)))

    # ⑥ the aggregate is a RATIO OF SUMS, not a mean of ratios. Two scores of very different mass
    #    make the two answers differ, which is what makes this a test rather than a tautology.
    a = {"n_scored": 1, "mass": 1.0, "net_err": 0.0, "abs_err": 1.0,
         "over_call": 0.0, "under_call": 0.0, "mwae": 1.0}
    b = {"n_scored": 1, "mass": 999.0, "net_err": 0.0, "abs_err": 0.0,
         "over_call": 0.0, "under_call": 0.0, "mwae": 0.0}
    agg = _agg_axis([a, b])
    check("aggregate mwae is Σ|err| / Σmass", np.isclose(agg["mwae"], 1.0 / 1000.0))
    check("aggregate mwae is NOT the mean of the rates", not np.isclose(agg["mwae"], 0.5))
    check("aggregate returns None on an empty selection", _agg_axis([]) is None)

    # ⑦ ⭐ the SIGNED HISTOGRAM: it must place a known vector exactly, keep the sign, and never
    #    silently drop an object — the bucket counts must sum to the input size.
    known = np.array([0.0, 0.0, -5.0, -500.0, +2.0, +2000.0, -0.05, +0.05, +20000.0, -20000.0])
    h = signed_histogram(known)
    check("histogram conserves every object", sum(b["n"] for b in h["buckets"]) == known.size)
    check("histogram counts the exact zeros", h["n_exact"] == 2)
    check("histogram splits under/over correctly", (h["n_under"], h["n_over"]) == (4, 4))
    check("histogram Σ|err| equals Σ|input|",
          np.isclose(sum(b["abs_err"] for b in h["buckets"]), np.abs(known).sum()))
    # ⛔ and it must MOVE an object between buckets when the value crosses an edge — a histogram that
    # bins everything into one bar would pass every check above.
    lo = signed_histogram(np.array([9.99]))
    hi = signed_histogram(np.array([10.01]))
    check("histogram resolves a bucket EDGE crossing",
          [b["n"] for b in lo["buckets"]] != [b["n"] for b in hi["buckets"]])
    check("histogram is SIGN-sensitive: +x and -x land in mirrored buckets",
          [b["n"] for b in signed_histogram(np.array([+7.0]))["buckets"]]
          == [b["n"] for b in signed_histogram(np.array([-7.0]))["buckets"]][::-1])

    # ⑧ the pool ledger reads the FILE, and refuses to invent a pool it cannot find.
    import json as _json
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "truth_summary.json").write_text(
            _json.dumps({"origin_counts": {"gdna": 7.0, "mrna": 11.0, "nrna": 13.0}})
        )
        led = pool_ledger(Path(td))
        check("pool ledger reads all three origin pools",
              (led["gdna"], led["mrna"], led["nrna"]) == (7.0, 11.0, 13.0))
        (Path(td) / "truth_summary.json").write_text(_json.dumps({"origin_counts": {"gdna": 7.0}}))
        check("pool ledger reports a MISSING pool as 0, not as absent",
              pool_ledger(Path(td)) == {"gdna": 7.0, "mrna": 0.0, "nrna": 0.0})

    # ⑨ ⭐ THE VERTEX PROFILE: it must place objects in the right bucket by TRUE f_g, exclude the
    #    massless, conserve Σ|Δ| against the axis total, and RESOLVE a shortfall that is planted.
    # ⚠ its OWN toy, sized to the slices below: the shared one has 24 regions and a [20:30] slice would
    # silently take four. That is the defect this block exists to catch, met while writing it.
    cal_v = _toy_calibration(n_regions=40, n_boundaries=36, seed=11)
    n_obj = cal_v.n_regions
    g_true = np.zeros(n_obj)
    r_true = np.zeros(n_obj)
    g_true[:10], r_true[:10] = 0.0, 100.0  # truth f_g = 0.0  -> bottom bucket
    g_true[10:20], r_true[10:20] = 50.0, 50.0  # truth f_g = 0.5  -> middle bucket
    g_true[20:30], r_true[20:30] = 100.0, 0.0  # truth f_g = 1.0  -> ⭐ the VERTEX bucket
    # everything else stays massless and must be dropped entirely
    o_toy = dataclasses.replace(cal_v, mass_gdna_region=g_true, mass_rna_region=r_true)
    # a perfect arm: predictions equal to truth
    prof = vertex_profile(o_toy, o_toy, "region")
    check("vertex profile drops massless objects",
          sum(b["n"] for b in prof) == 30)
    check("vertex profile buckets by TRUE f_g",
          (prof[0]["n"], prof[1]["n"], prof[-1]["n"]) == (10, 10, 10))
    check("a perfect arm has zero shortfall in every bucket",
          all(abs(b["mean_shortfall"]) < 1e-12 for b in prof))
    # ⛔ plant a shortfall ONLY at the vertex and require it to appear THERE and nowhere else
    g_short = g_true.copy()
    r_short = r_true.copy()
    g_short[20:30], r_short[20:30] = 80.0, 20.0  # pred f_g = 0.8 where truth is 1.0
    p_toy = dataclasses.replace(cal_v, mass_gdna_region=g_short, mass_rna_region=r_short)
    prof2 = vertex_profile(p_toy, o_toy, "region")
    check("a planted vertex shortfall lands in the VERTEX bucket",
          abs(prof2[-1]["mean_shortfall"] - 0.2) < 1e-9)
    check("and leaves the other buckets untouched",
          all(abs(b["mean_shortfall"]) < 1e-12 for b in prof2[:-1]))
    check("vertex Σ|Δ| conserves the axis total",
          abs(sum(b["abs_err"] for b in prof2) - float(np.abs(g_short - g_true).sum())) < 1e-9)
    # ⛔ closure is read off the PREDICTING arm, and must move when its composition does
    open_arm = dataclasses.replace(p_toy, rna_pos_frac_region=np.zeros(n_obj))
    check("closure is read from the arm under test, not from truth",
          vertex_profile(open_arm, o_toy, "region")[-1]["mean_closure"]
          != prof2[-1]["mean_closure"])

    # ⑩ every stratum of the panel is covered exactly once by the scope table.
    check("the scope table covers all four strata", len(_SCOPE) == 4)
    check("three strata are IN SCOPE and one is DEFERRED",
          sorted(_SCOPE.values()) == ["DEFERRED", "IN SCOPE", "IN SCOPE", "IN SCOPE"])

    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"  {'✅' if ok else '⛔'} {name:<{width}}")
    n_ok = sum(1 for _n, ok in checks if ok)
    print(f"\n  {n_ok}/{len(checks)} self-test gates fired as specified")
    return 0 if n_ok == len(checks) else 1


# ── cli ──────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--work-dir", type=Path,
                    default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_cal_vs_oracle")
    ap.add_argument("--json", type=Path, default=None, help="write the per-condition rows and exit")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument(
        "--background-abundance",
        choices=("contained", "measured_total"),
        default=None,
        help="override CalibrationConfig.background_abundance for BOTH arms (P and O), so the run "
        "prices that estimator swap against oracle calibration. Default: the shipped config.",
    )
    ap.add_argument("--self-test", action="store_true", help="perturb every comparator; no I/O")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    if not names:
        raise SystemExit(f"⛔ no condition with a sim_oracle.bam under {args.suite}")
    cache = args.oracle_cache
    if cache is None and (args.suite / "oracle_cache").is_dir():
        cache = args.suite / "oracle_cache"
    if cache is None:
        raise SystemExit("⛔ --oracle-cache is required; this script refuses to invent a truth")

    if args.jobs > 1 and args.json is None:
        # ⭐ Shard by SUBPROCESS, exactly as prior_vs_oracle.py does — then merge and print through
        # the SAME report path, so a sharded run and a serial one cannot print different numbers.
        shards = [s for s in (names[i:: args.jobs] for i in range(args.jobs)) if s]
        tmp = args.work_dir / "_shards"
        tmp.mkdir(parents=True, exist_ok=True)
        procs = []
        for i, sh in enumerate(shards):
            out = tmp / f"shard{i}.json"
            cmd = [sys.executable, str(Path(__file__).resolve()),
                   "--suite", str(args.suite), "--index", str(args.index),
                   "--oracle-cache", str(cache), "--json", str(out), "--conditions", *sh]
            procs.append((subprocess.Popen(cmd), out))
        merged: list[dict] = []
        for proc, out in procs:
            if proc.wait() != 0:
                raise SystemExit("⛔ a shard failed; refusing to report a partial panel")
            merged += json.loads(out.read_text())
        report(merged)
        return 0

    t0 = time.perf_counter()
    index = TranscriptIndex.load(str(args.index))
    region_arrays = RegionArrays.from_index(index)
    pipeline_config = PipelineConfig()
    if args.background_abundance is not None:
        # ⛔ Applied to BOTH arms: P and O share one payload and differ only in the six deconvolved
        # arrays, so an override on one arm alone would compare two different tools rather than two
        # estimators.
        pipeline_config = dataclasses.replace(
            pipeline_config,
            calibration=dataclasses.replace(
                pipeline_config.calibration, background_abundance=args.background_abundance
            ),
        )
        print(f"⭐ background_abundance = {args.background_abundance!r} on BOTH arms")
    if args.json is None:
        print(f"index + region arrays loaded in {time.perf_counter() - t0:.2f} s", flush=True)

    rows = []
    for name in names:
        rows.append(
            measure_condition(index, region_arrays, pipeline_config, args.suite, cache, name)
        )
        if args.json is None:
            print(f"  scored {name}  ({rows[-1]['seconds']:.1f} s)", flush=True)

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows))
        return 0
    report(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
