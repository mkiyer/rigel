#!/usr/bin/env python
"""⭐⭐⭐ **THE TOOL'S ABSOLUTE ACCURACY, AND THE RE-SOLVE CEILING ABOVE IT** — one instrument, two jobs.

``--arm base`` alone is **tool-absolute-accuracy**: run the shipped pipeline end to end on a simulated
condition and score the transcript table against the simulator's own per-transcript truth.
``--arm oracle`` is **error-downstream-of-calibration**: run the identical pipeline with the ORACLE
prior injected in place of calibration's, and re-quantify. The difference between the two arms is the
whole of what a perfect calibration would buy, and it is the measurement that decides whether the work
belongs in calibration or in the EM.

⭐ **ONE SCRIPT BECAUSE IT IS ONE SCORER.** The two questions differ only in which ``LocusPriors`` the
EM is handed; scoring them with two implementations is how the ceiling and the baseline drift apart.

THE ARMS
--------
===============  ===============================================================================
``base``         the shipped pipeline, untouched
``noop``         the injection wrapper runs in full — the oracle is read, O is built and then
                 DISCARDED — and the shipped prior is returned. ⛔ MUST be byte-identical to
                 ``base``; this is the harness's own falsification (TRAPS: byte-identity-gate)
``oracle``       all three arrays replaced by O
``oracle_gdna``  only ``gdna_prior_count``      ⭐ the three single-array arms are what say
``oracle_rna``   only ``rna_prior_count``          WHICH of the prior's three numbers carries the
``oracle_efflen``only ``gdna_eff_len``             value — a joint ceiling cannot
===============  ===============================================================================

⛔ **WHAT "TRUTH" MEANS HERE, BECAUSE THE PANEL SHIPS TWO AND THEY ARE NOT INTERCHANGEABLE.** Each
condition carries its own ``truth_abundances.tsv`` whose ``mrna_abundance`` column is the **realised
observed fragment count** for that transcript in that condition (verified: it equals
``observed_mrna_fragments`` exactly and sums to the condition's RNA fragment total). The suite-level
``truth_abundances_nrna_none.tsv`` is the **pre-capture molar** abundance — a different quantity, and
under capture a very different one (log-space correlation 0.72 between them at
``g25 ss0.50 capture_on``). Scoring an observed-fragment estimate against a pre-capture molar truth
would charge the tool for hybrid capture, which it never claims to invert. ⭐ So the primary score is
**count against count**, with no length model anywhere in the comparison.

⚠ **The TPM rows are secondary and they measure ASSIGNMENT, not the length model.** Truth TPM is built
with the tool's OWN ``effective_length``, so the two sides share the length model and it cancels. That
is deliberate: the fragment-length models have their own instrument (``length_ceiling.py``) and mixing
the two questions is how a 14x split between two length pmfs stayed hidden
(TRAPS: price-the-halves-separately).

⛔ **THE ORACLE MASSES ARE UNDRAINED AND THE SHIPPED PIPELINE DRAINS, AND THAT BIAS IS MEASURED, NOT
ASSUMED.** ``prior_vs_oracle.py`` prices it on this panel at **0.153 %** of the gDNA prior and
**0.462 %** of the RNA prior. A drained oracle is inadmissible here (3.92 % of held fragments tie on
the deferred bank's key across origins, depositing spliced records into the gdna partition, which
``OracleTruth`` refuses). ⭐ The direction is conservative for either conclusion: a 0.2 % low oracle
cannot explain a large surviving error, and it can only understate a large removed one.

Gates: ``tests/calibration/test_quant_accuracy.py``. The byte-identity gate is
``scripts/design/arm_identity.py base noop`` — the rows are keyed ``(condition, axis)`` on purpose so
that already-falsified instrument can be pointed straight at them.

Usage::

    python scripts/design/quant_accuracy.py --arm base --out $RIGEL_ARMS/qa_base.jsonl --jobs 4
    python scripts/design/quant_accuracy.py --arm oracle --out $RIGEL_ARMS/qa_oracle.jsonl --jobs 4
    python scripts/design/arm_identity.py $RIGEL_ARMS/qa_base.jsonl $RIGEL_ARMS/qa_noop.jsonl
    python scripts/design/quant_accuracy.py --report $RIGEL_ARMS/qa_base.jsonl $RIGEL_ARMS/qa_oracle.jsonl
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "tests" / "calibration"))

from _oracle import ORIGINS, OracleTruth  # noqa: E402

import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, run_pipeline  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: ⭐⭐ ``warm_uniform`` seeds every component EQUALLY instead of by coverage-weighted share. It varies
#: ONE thing against ``base`` — where the EM is put down — so a difference is the seed's BASIN and
#: nothing else. ⛔ It exists because the true per-transcript split was measured to be a fixed point of
#: the prior-free objective (drifting 1.0 % of RNA mass over 200 iterations) while the shipped answer
#: sits 444,595 nats WORSE in data likelihood: the solver is not failing on ambiguity, it is landing
#: somewhere else, and this says whether the seed is what puts it there.
#: ⭐⭐⭐ ``oracle_alloc`` IS STAGE 5, AND IT IS A CAPABILITY PROOF RATHER THAN A CEILING. It hands the
#: EM the TRUE relative transcript abundances as the per-transcript allocation weights and zeroes the
#: coverage seed, so `theta` starts from the prior alone. The question is narrow and the answer is
#: binary: *given correct weights, does the machinery produce correct per-transcript abundances?* If it
#: does, a weighting function is worth designing; if it does not, no weighting function can help and the
#: work stops here.
#:
#: ⛔ **It is NOT a claim about reachable headroom**, and must not be quoted as one. Truth-derived
#: weights also hand over the true SUPPORT — a zero weight is exactly absorbing, so every silent
#: transcript is switched off for free, and silent transcripts carry 23.1 % of the transcript error at
#: g00. Pricing what a REAL weighting function could earn needs controls this arm deliberately omits.
ARMS = ("base", "base_reseed", "noop", "oracle", "oracle_gdna", "oracle_rna", "oracle_efflen",
        "warm_uniform", "oracle_alloc", "oracle_alloc_seed", "oracle_alloc_flip")

#: ⛔⛔ **THE SHIPPED PIPELINE IS NOT REPRODUCIBLE RUN TO RUN, AND EVERY ARM HERE PINS THE SEED THAT
#: MAKES IT SO.** ``EMConfig.seed`` defaults to ``None`` and ``assignment_mode`` to ``"sample"``, so
#: the EM's final hard assignment is an UNSEEDED categorical draw: measured on the gate toy, two
#: back-to-back runs of the identical pipeline differ on 4 transcripts by up to 43 fragments. A
#: byte-identical ``noop`` is therefore impossible on the default config — which is why this instrument
#: sets a seed rather than reporting a difference it cannot attribute. ⭐ ``base_reseed`` re-runs
#: ``base`` at ``seed + 1`` and is the NOISE FLOOR: any arm delta smaller than it is sampling.
DEFAULT_EM_SEED = 20260807

#: arm -> which ``LocusPriors`` fields come from O. ``noop`` takes NONE of them and still builds O,
#: which is what makes it a test of the wrapper rather than of an ``if``.
_ARM_FIELDS = {
    "base": (),
    "base_reseed": (),
    "noop": (),
    "oracle": ("gdna_prior_count", "rna_prior_count", "gdna_eff_len"),
    "oracle_gdna": ("gdna_prior_count",),
    "oracle_rna": ("rna_prior_count",),
    "oracle_efflen": ("gdna_eff_len",),
}


# ── the injection ────────────────────────────────────────────────────────────────────────────────


def load_oracle(bam: str, index, pipeline_config, cache_root: Path, tag: str) -> OracleTruth:
    """The origin-split truth for one condition, entirely from the shipped scan cache.

    ⛔ ``_main`` (the UNDRAINED full payload) is what sum-to-full is asserted against, and it must come
    from the same cache as the three partitions or the identity is checking two different scans
    against each other. ⚠ ``read_scan_cache`` refuses a payload whose ``graph_hash`` / ``reach_digest``
    / ``payload_schema_digest`` / scan config does not describe this index — ``reach`` is covered by no
    other hash — so a stale cache is a loud refusal rather than a silent wrong truth.
    """
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    root = Path(cache_root) / tag
    try:
        full = read_scan_cache(root / "_main", index, scan).payload
        parts = {k: read_scan_cache(root / k, index, scan).payload for k in ORIGINS}
    except (FileNotFoundError, KeyError, ScanCacheKeyError) as exc:
        raise SystemExit(
            f"⛔ {tag}: no valid oracle cache under {root} ({exc}). Build it first with "
            "pass0_vs_oracle.py --oracle-cache or prior_vs_oracle.py --oracle-cache; this script "
            "refuses to invent a truth."
        ) from exc
    return OracleTruth.from_parts(full, parts)


def install_arm(arm: str, oracle: OracleTruth | None):
    """Wrap :func:`~rigel.calibration.priors.assemble_priors` for one arm. Returns ``(restore, fired)``.

    ⭐ **The wrapper is IDENTICAL across every arm including ``noop``** — the oracle is read, O is
    built on the loci the run itself produced, and only then does the arm decide which fields to take.
    A ``noop`` that short-circuits before building O would prove the ``if`` works and nothing else
    (TRAPS: could-the-arm-have-fired); this one proves the whole path is inert when it takes no field.

    ⛔ **O is built on the run's OWN ``multi_loci``, never on a stored array.** ``build_multi_loci``
    unions transcripts linked by SCORED fragments, so the locus partition is a function of the run —
    an oracle prior keyed by ``multi_locus_id`` from another process is not portable and index-aligning
    it would compare locus 7 of one run with locus 7 of another.

    ⚠ ``fired`` is a mutable counter, not a bool, and the caller RAISES on zero: an override that never
    ran reads as "no effect", which is the most flattering possible failure (TRAPS: an-ablation-that-never-ran).
    """
    original = PRIORS.assemble_priors

    if arm in ("base", "base_reseed", "warm_uniform"):
        # ⚠ Counted as fired: ``base`` installs nothing by design, so the "did the override run?"
        # check must not fail on the one arm that has no override. The thing it guards — an
        # injection that silently did not happen — cannot occur here.
        return (lambda: None), {"n": 1}

    fields = _ARM_FIELDS[arm]

    fired = {"n": 0}

    def wrapper(calibration, region_arrays, multi_loci):
        shipped = original(calibration, region_arrays, multi_loci)
        o = original(
            dataclasses.replace(calibration, **oracle.override_masses(region_arrays)),
            region_arrays,
            multi_loci,
        )
        fired["n"] += 1
        if not fields:
            return shipped
        return dataclasses.replace(shipped, **{f: getattr(o, f) for f in fields})

    PRIORS.assemble_priors = wrapper

    def restore():
        PRIORS.assemble_priors = original

    return restore, fired


# ── the scoring ──────────────────────────────────────────────────────────────────────────────────


def _spearman_pearson(true_v: np.ndarray, est_v: np.ndarray) -> tuple[float, float]:
    from scipy.stats import pearsonr, spearmanr

    if true_v.size < 3:
        return float("nan"), float("nan")
    sp = float(spearmanr(true_v, est_v).statistic)
    pe = float(pearsonr(np.log2(true_v + 1.0), np.log2(est_v + 1.0)).statistic)
    return sp, pe


def score_transcripts(quant: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """Score one arm's transcript table against the condition's realised fragment truth.

    ⭐ **``count`` against ``mrna_abundance``, both fragment counts.** No effective length, no
    normalisation, no model between the two sides — so a difference is a fragment the tool put
    somewhere the simulator did not.

    ⛔ **The false-positive mass is reported separately and it is the number this tool exists for.**
    ``fp_mass`` is the estimate summed over transcripts the simulator gave ZERO fragments: every one of
    those is contamination the deconvolution failed to remove, with nothing to cancel it. Its mirror
    ``fn_mass`` is real RNA the tool called gDNA. A single ``Σ|Δ|`` cannot separate them and the two
    are different failures with different fixes.
    """
    m = truth.merge(
        quant[["transcript_id", "count", "count_em", "tpm", "effective_length"]],
        on="transcript_id",
        how="left",
    )
    for col in ("count", "count_em", "tpm", "effective_length"):
        m[col] = m[col].fillna(0.0)
    t = m["mrna_abundance"].to_numpy(np.float64)  # the realised OBSERVED fragment count
    e = m["count"].to_numpy(np.float64)
    d = e - t
    expressed = t > 0.0
    silent = ~expressed

    # TPM. ⚠ Truth TPM uses the TOOL's own effective length, so the length model is common to both
    # sides and cancels — this row is about ASSIGNMENT (see the module docstring).
    eff = np.maximum(m["effective_length"].to_numpy(np.float64), 1.0)
    rate = t / eff
    tpm_true = rate / rate.sum() * 1e6 if rate.sum() > 0 else np.zeros_like(rate)
    tpm_est = m["tpm"].to_numpy(np.float64)

    sp, pe = _spearman_pearson(t[expressed], e[expressed])
    # MARD over the expressed set — bounded in [0,1], finite at a zero estimate.
    denom = np.abs(e[expressed]) + np.abs(t[expressed])
    mard = float(np.mean(np.abs(d[expressed]) / np.where(denom > 0, denom, 1.0))) if expressed.any() else float("nan")
    rel = np.abs(d[expressed]) / t[expressed] if expressed.any() else np.array([])

    return {
        "n_tx": int(len(m)),
        "n_expressed": int(expressed.sum()),
        "n_detected": int((e > 0).sum()),
        "count_true": float(t.sum()),
        "count_est": float(e.sum()),
        "count_abs_err": float(np.abs(d).sum()),
        "count_net_err": float(d.sum()),
        "count_over": float(np.maximum(d, 0.0).sum()),
        "count_under": float(np.maximum(-d, 0.0).sum()),
        #: ⭐ contamination the deconvolution left on a silent transcript — a pure false positive
        "fp_mass": float(e[silent].sum()),
        "fp_n": int((e[silent] > 0).sum()),
        #: real RNA the tool did not assign to the transcript that produced it
        "fn_mass": float(t[expressed & (e == 0.0)].sum()),
        "fn_n": int((expressed & (e == 0.0)).sum()),
        "spearman": sp,
        "pearson_log2": pe,
        "mard": mard,
        "median_rel_err": float(np.median(rel)) if rel.size else float("nan"),
        "tpm_abs_err": float(np.abs(tpm_est - tpm_true).sum()),
        "tpm_fp": float(tpm_est[silent].sum()),
    }


def score_genes(quant: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """The SAME scorer, over genes — and the difference from the transcript row is the whole point.

    ⭐⭐ **THIS IS THE DISCRIMINATOR THAT SEPARATES "ISOFROM AMBIGUITY" FROM "SOMETHING ELSE", AND
    WITHOUT IT THE TRANSCRIPT NUMBER CANNOT BE ACTED ON.** Summing a gene's isoforms collapses exactly
    the error that comes from not knowing WHICH isoform a fragment came from — the difficulty every
    transcript quantifier has and that no prior can remove. What survives at gene level is error in
    deciding whether the fragment was RNA *from this gene* at all, which is the question Rigel is for.

    ⛔ Grouped on ``gene_id`` from the TRUTH table, not from the index, so the two sides are grouped by
    one definition. ⚠ Synthetic nRNA entities are absent from ``quant`` by construction
    (``get_counts_df`` drops them), so their mass is missing from the gene row too — it is reported on
    the ``library`` row as ``nrna_est`` and must be read there.
    """
    m = truth.merge(quant[["transcript_id", "count"]], on="transcript_id", how="left")
    m["count"] = m["count"].fillna(0.0)
    g = m.groupby("gene_id", observed=True)[["mrna_abundance", "count"]].sum()
    t = g["mrna_abundance"].to_numpy(np.float64)
    e = g["count"].to_numpy(np.float64)
    d = e - t
    expressed = t > 0.0
    silent = ~expressed
    sp, pe = _spearman_pearson(t[expressed], e[expressed])
    denom = np.abs(e[expressed]) + np.abs(t[expressed])
    mard = (
        float(np.mean(np.abs(d[expressed]) / np.where(denom > 0, denom, 1.0)))
        if expressed.any()
        else float("nan")
    )
    return {
        "n_tx": int(len(g)),  # genes, named to match the transcript row so one table renders both
        "n_expressed": int(expressed.sum()),
        "n_detected": int((e > 0).sum()),
        "count_true": float(t.sum()),
        "count_est": float(e.sum()),
        "count_abs_err": float(np.abs(d).sum()),
        "count_net_err": float(d.sum()),
        "count_over": float(np.maximum(d, 0.0).sum()),
        "count_under": float(np.maximum(-d, 0.0).sum()),
        "fp_mass": float(e[silent].sum()),
        "fp_n": int((e[silent] > 0).sum()),
        "fn_mass": float(t[expressed & (e == 0.0)].sum()),
        "fn_n": int((expressed & (e == 0.0)).sum()),
        "spearman": sp,
        "pearson_log2": pe,
        "mard": mard,
    }


def score_library(result, quant: pd.DataFrame, truth_summary: dict) -> dict:
    """The library-level split — the thermometer, not the target.

    ⛔⛔ **INTERGENIC FRAGMENTS ARE gDNA AND THE DENOMINATOR MUST SAY SO.** ``n_intergenic`` counts
    fragments that reached no locus at all; they never enter the EM, so ``estimator.gdna_em_count``
    excludes them. Off capture, gDNA is genome-uniform and **more than half of it is intergenic**
    (measured at `g50 ss0.50 capture_off`: 2,601,271 intergenic against 2,330,992 EM-assigned gDNA), so
    leaving them out of the numerator while RNA stays in the denominator understates the fraction
    enormously — it read **0.3151 against a truth of 0.5000**, and that fabricated a "systematic
    off-capture EM under-call" that does not exist. With them, the same run reads **0.4933**.
    ⭐ This is ``cli.py``'s ``gdna_fraction`` — ``(gdna_em + n_intergenic) / (rna + gdna_em +
    n_intergenic)`` — so the number here is the one a user reads. ⚠ ``gdna_em_frac_est`` is kept
    beside it because the two answer different questions and only one of them is the deliverable.

    The truth is ``origin_counts`` from the simulator, which counts each fragment exactly once.
    """
    est = result.estimator
    mrna = float(quant["count"].sum())
    nrna = float(est.nrna_em_count)
    gdna = float(est.gdna_em_count)
    intergenic = float(result.stats.n_intergenic)
    oc = truth_summary.get("origin_counts", {})
    t_g = float(oc.get("gdna", float("nan")))
    t_m = float(oc.get("mrna", 0.0))
    t_n = float(oc.get("nrna", 0.0))
    t_r = t_m + t_n
    total = mrna + nrna + gdna
    total_all = total + intergenic
    return {
        # ⭐ THE NASCENT ARM IS BROKEN OUT AND IT IS NOT COSMETIC. ``get_counts_df`` drops the
        # SYNTHETIC nRNA entities, so every fragment the EM parks on one is invisible in the
        # transcript table's ``fp_mass``. On an ``nrna_none`` panel the truth is exactly 0, which
        # makes ``nrna_est`` a THIRD false-positive channel with nothing to cancel it — and the
        # first measurement showed it larger than the annotated one.
        "mrna_est": mrna,
        "mrna_true": t_m,
        "nrna_est": nrna,
        "nrna_true": t_n,
        "rna_est": mrna + nrna,
        "rna_true": t_r,
        "gdna_est": gdna,
        "gdna_true": t_g,
        #: ⭐ THE DELIVERABLE — intergenic included, cli.py's denominator.
        "gdna_frac_est": (gdna + intergenic) / total_all if total_all > 0 else float("nan"),
        #: the EM's own view, intergenic in NEITHER side. A different question; not the deliverable.
        "gdna_em_frac_est": gdna / total if total > 0 else float("nan"),
        "gdna_frac_true": t_g / (t_g + t_r) if (t_g + t_r) > 0 else float("nan"),
        "n_intergenic": intergenic,
    }


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def seeded(pipeline_config, arm: str, em_seed: int):
    """The arm's pipeline config. ⭐ ``base_reseed`` differs from ``base`` in the SEED ALONE, so the
    gap between them is the sampling noise of the EM's own hard assignment and nothing else.

    """
    seed = em_seed + 1 if arm == "base_reseed" else em_seed
    warm = pipeline_config.em.warm_start
    if arm == "warm_uniform":
        warm = "uniform"
    elif arm == "oracle_alloc":
        # ⭐ the seed is zeroed so `theta` starts proportional to the prior ALONE — otherwise a
        # coverage-weighted seed is MULTIPLIED by an allocation from a different method and the result
        # is neither method's answer.
        warm = "prior"
    return dataclasses.replace(
        pipeline_config, em=dataclasses.replace(pipeline_config.em, seed=seed, warm_start=warm)
    )


def truth_weights(truth: pd.DataFrame, index) -> np.ndarray:
    """``float64[n_transcripts]`` — the TRUE realised fragment count per transcript, on the EM's axis.

    ⭐ This is the whole of stage 5's input: the relative abundances the simulator actually produced.
    Within a locus the EM only reads their RATIOS, so no normalisation is needed here.

    ⚠ Exact-duplicate transcripts are folded onto the twin the index kept — the truth table is keyed on
    the un-collapsed annotation, so without the fold their fragments would be dropped rather than
    attributed. ⛔ For such a pair the per-transcript truth is not merely awkward, it is UNDEFINED: the
    two are the same molecule and only the group total is a fact about the world.
    """
    col = next((c for c in ("observed_mrna_fragments", "mrna_abundance") if c in truth.columns), None)
    if col is None:
        raise SystemExit("⛔ truth table has neither observed_mrna_fragments nor mrna_abundance")
    t_index = dict(zip(index.t_df["t_id"].to_numpy(), index.t_df["t_index"].to_numpy(), strict=True))
    w = np.zeros(int(index.num_transcripts), dtype=np.float64)
    for tid, n in zip(truth["transcript_id"], truth[col], strict=True):
        i = t_index.get(str(tid))
        if i is not None:
            w[int(i)] += float(n)
    return w


def install_truth_weights(weights: np.ndarray):
    """Pass per-transcript weights into the EM, and COUNT THE ARM AT THE DEEPEST POINT IT CAN OBSERVE.

    ⭐ A wrapper rather than a config field: the weights are a per-run ARRAY, and an experiment that
    injects one should not add production surface.

    ⛔⛔ **THE COUNTER WATCHES THE ESTIMATOR, NOT THE INJECTION, AND THAT DISTINCTION COST A SESSION.**
    An earlier version counted nonzero entries in ``weights`` — a property of the array this function
    was handed, true whatever the pipeline then did with it. Meanwhile
    ``pipeline._run_locus_em_partitioned`` accepted ``rna_prior_weight`` and silently dropped it, so
    every allocation, however extreme, produced byte-identical output while the counter read healthy.
    ⭐ Counting where the SOLVER receives the array makes a dropped parameter impossible to miss, and
    ``--arm oracle_alloc_flip`` (a maximally WRONG allocation, which must move the answer) is the
    end-to-end half of the same check.
    """
    import rigel.pipeline as PL
    from rigel.estimator import AbundanceEstimator

    inner = PL._run_locus_em_partitioned
    inner_em = AbundanceEstimator.run_batch_locus_em_partitioned
    fired = {"n": 0}

    def wrapper(*args, **kw):
        kw["rna_prior_weight"] = weights
        return inner(*args, **kw)

    def em_wrapper(self, *args, **kw):
        w = kw.get("rna_prior_weight")
        if w is not None and int(np.asarray(w).size) and float(np.asarray(w).sum()) > 0.0:
            fired["n"] += 1
        return inner_em(self, *args, **kw)

    PL._run_locus_em_partitioned = wrapper
    AbundanceEstimator.run_batch_locus_em_partitioned = em_wrapper

    def restore():
        PL._run_locus_em_partitioned = inner
        AbundanceEstimator.run_batch_locus_em_partitioned = inner_em

    return restore, fired


def run_condition(arm: str, suite: Path, index, condition: str, pipeline_config,
                  oracle_cache: Path | None, em_seed: int = DEFAULT_EM_SEED) -> list[dict]:
    bam = str(suite / condition / "sim_oracle.bam")
    truth = pd.read_csv(suite / condition / "truth_abundances.tsv", sep="\t")
    summary = json.loads((suite / condition / "truth_summary.json").read_text())
    pipeline_config = seeded(pipeline_config, arm, em_seed)

    oracle = None
    if arm not in ("base", "base_reseed", "warm_uniform", "oracle_alloc", "oracle_alloc_seed", "oracle_alloc_flip"):
        if oracle_cache is None:
            raise SystemExit(f"⛔ arm {arm!r} needs --oracle-cache")
        oracle = load_oracle(bam, index, pipeline_config, oracle_cache, condition)

    if arm == "oracle_alloc_flip":
        w = truth_weights(truth, index)
        # ⛔ FALSIFICATION OF THE HARNESS, not a treatment: put the weight where the truth is NOT.
        # If a maximally wrong allocation moves nothing, the allocation never reached the solver.
        flip = np.zeros_like(w)
        nz = np.flatnonzero(w > 0)
        if nz.size:
            flip[nz] = w[nz][::-1]
        restore, fired = install_truth_weights(flip)
    elif arm.startswith("oracle_alloc"):
        # ⭐ `oracle_alloc_seed` keeps the SHIPPED coverage seed, so it varies ONE thing against `base`:
        # the allocation. `oracle_alloc` also zeroes the seed, which additionally removes RNA's evidence
        # advantage over calibration's gDNA prior — two changes, not one.
        restore, fired = install_truth_weights(truth_weights(truth, index))
    else:
        restore, fired = install_arm(arm, oracle)
    start = time.perf_counter()
    try:
        result = run_pipeline(bam, index, pipeline_config)
    finally:
        restore()
    if fired["n"] == 0:
        # ⛔ TRAPS: an-ablation-that-never-ran — an injection that never ran reads as "a perfect prior changes nothing".
        raise RuntimeError(
            f"{condition} [{arm}]: assemble_priors was never wrapped-and-called. This is not a "
            "measurement of zero effect."
        )
    seconds = time.perf_counter() - start

    quant = result.estimator.get_counts_df(index)
    common = {"arm": arm, "condition": condition, "seconds": seconds,
              "em_seed": int(pipeline_config.em.seed)}
    return [
        {**common, "axis": "transcript", **score_transcripts(quant, truth)},
        # ⭐ the SAME scorer over genes — isoform ambiguity summed away, see score_genes
        {**common, "axis": "gene", **score_genes(quant, truth)},
        {**common, "axis": "library", **score_library(result, quant, summary)},
    ]


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def stratum(cond: str) -> tuple[str, str]:
    return ("stranded" if "ss_0.99" in cond else "unstranded",
            "capture ON" if "capture_on" in cond else "capture OFF")


def is_zero_gdna(cond: str) -> bool:
    return "_g00_" in cond


_STRATA = (("stranded", "capture OFF"), ("stranded", "capture ON"),
           ("unstranded", "capture OFF"), ("unstranded", "capture ON"))


def _load(path: Path) -> dict:
    rows = [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]
    return {(r["condition"], r["axis"]): r for r in rows}


def report(paths: list[Path]) -> None:
    """One or more arms, per stratum. ⛔ Never pooled — the panel total hides a sign flip between
    strata, and on this panel one stratum carries almost all of the error."""
    arms = [(_load(p), Path(p).stem) for p in paths]
    keys = set(arms[0][0])
    for a, name in arms[1:]:
        if set(a) != keys:
            raise SystemExit(
                f"⛔ {name} has a different row set ({len(a)} vs {len(keys)}) — an arm that is "
                "missing conditions must not be aggregated against one that is not."
            )
    conds = sorted({c for c, _ax in keys})

    print()
    print("=" * 112)
    print("  ⭐⭐⭐ TRANSCRIPT-LEVEL ACCURACY, end to end, scored against the REALISED fragment truth")
    print(f"  {len(conds)} conditions   arms: {', '.join(n for _a, n in arms)}"
          "   messages OFF   length_likelihood OFF")
    print("=" * 112)

    def block(title, field, fmt="{:>14,.0f}", ratio=True, axis="transcript"):
        print()
        print(f"  {title}")
        head = f"    {'stratum':<26}"
        for _a, name in arms:
            head += f" {name[:16]:>16}"
        if ratio and len(arms) > 1:
            head += f" {'arm/base':>10}"
        print(head)
        print("    " + "-" * (26 + 17 * len(arms) + 11))

        def line(label, sel):
            vals = []
            for a, _n in arms:
                sub = [a[(c, axis)][field] for c in conds if sel(c) and (c, axis) in a]
                vals.append(sum(sub) if sub else float("nan"))
            row = f"    {label:<26}"
            for v in vals:
                row += " " + fmt.format(v).rjust(16)
            if ratio and len(arms) > 1:
                row += f" {vals[-1] / vals[0]:>10.3f}" if vals[0] else f" {'—':>10}"
            print(row)

        line("ALL (g00 excluded)", lambda c: not is_zero_gdna(c))
        for st in _STRATA:
            line(" x ".join(st), lambda c, st=st: stratum(c) == st and not is_zero_gdna(c))
        print("    " + "-" * (26 + 17 * len(arms) + 11))
        line("⛔ g00 ZERO-gDNA control", is_zero_gdna)

    block("① TOTAL MISASSIGNED FRAGMENTS  Σ|count_est − count_true|  ·  TRANSCRIPT level",
          "count_abs_err")
    if any((conds[0], "gene") in a for a, _n in arms):
        block("①g   … the SAME quantity at GENE level — isoform ambiguity summed away, so what is "
              "left is\n       error in deciding whether the fragment was RNA from this gene at all",
              "count_abs_err", axis="gene")
    block("② ⭐ FALSE-POSITIVE MASS — estimate on transcripts the simulator gave ZERO fragments",
          "fp_mass")
    block("③ FALSE-NEGATIVE MASS — real RNA the tool assigned nowhere", "fn_mass")
    block("④ TPM error, Σ|Δ| out of 1e6 (assignment only — see the docstring)", "tpm_abs_err")

    # rank statistics do not sum; report the mass-free mean and say so
    print()
    print("  ⑤ RANK / SHAPE (a mean over conditions — these do not add)")
    head = f"    {'stratum':<26}"
    for _a, name in arms:
        head += f" {name[:10] + ' sprmn':>17}{name[:10] + ' mard':>17}"
    print(head)
    print("    " + "-" * (26 + 34 * len(arms)))

    def rank_line(label, sel):
        row = f"    {label:<26}"
        for a, _n in arms:
            sub = [a[(c, "transcript")] for c in conds if sel(c)]
            sp = float(np.nanmean([x["spearman"] for x in sub])) if sub else float("nan")
            md = float(np.nanmean([x["mard"] for x in sub])) if sub else float("nan")
            row += f" {sp:>16.4f} {md:>16.4f}"
        print(row)

    rank_line("ALL (g00 excluded)", lambda c: not is_zero_gdna(c))
    for st in _STRATA:
        rank_line(" x ".join(st), lambda c, st=st: stratum(c) == st and not is_zero_gdna(c))
    print("    " + "-" * (26 + 34 * len(arms)))
    rank_line("⛔ g00 ZERO-gDNA control", is_zero_gdna)

    # ── ⑥ THE POOL LEVEL ─────────────────────────────────────────────────────────────────────────
    #
    # ⛔⛔ AXIOM 0 GUARD, AND IT IS THE WHOLE REASON THIS TABLE NEEDS A HEADER. The solver has THREE
    # populations — gDNA, RNA+, RNA− — and "nascent" is NOT one of them and never becomes one. This is
    # an ASSIGNMENT question, not a deconvolution one: `nrna_est` is mass the EM parked on nascent
    # entities, and reading it as a fourth component is the exact error `CLAUDE.md`'s AXIOM 0 exists to
    # prevent.
    # ⭐ On a `nrna_none` panel `nrna_true` is exactly 0, so every fragment in that column is a false
    # positive with NOTHING TO CANCEL IT — the same one-sided logic that makes the `g00` control
    # readable (`TRAPS: zero-target-guards-are-one-sided`), and it is why the column is worth printing
    # even on a panel that contains no nascent RNA at all.
    print()
    print("  ⑥ ⭐⭐⭐ POOL LEVEL — where the library's fragments went: gDNA vs NASCENT vs ANNOTATED")
    print("     ⛔ 'nascent' is a SIMULATOR input, never a population the solver has (AXIOM 0). This")
    print("        scores ASSIGNMENT. On a nrna_none panel nascent truth is 0, so that column is pure")
    print("        false positive and Δ% is undefined rather than large.")
    _POOL_FIELDS = ("gdna_est", "gdna_true", "nrna_est", "nrna_true", "mrna_est", "mrna_true",
                    "n_intergenic")
    for a, name in arms:
        print()
        print(f"    arm: {name}")
        # ⚠ An arm file written before this table existed has no pool fields. Say so and skip, rather
        # than dying on a KeyError halfway through a report whose other seven tables are fine.
        probe = next((a[(c, "library")] for c in conds if (c, "library") in a), None)
        if probe is None or any(f not in probe for f in _POOL_FIELDS):
            print("      ⚠ no pool-level fields in this arm — re-run it to get table ⑥")
            continue
        print(f"    {'stratum':<26} {'pool':<20} {'est':>15} {'true':>15} {'Δ':>15} {'Δ%':>8}")
        print("    " + "-" * 102)

        def pool_rows(label, sel, a=a):
            def tot(field):
                return sum(a[(c, "library")][field] for c in conds
                           if sel(c) and (c, "library") in a)

            # ⛔⛔ INTERGENIC FRAGMENTS ARE gDNA AND THE NUMERATOR MUST SAY SO. `gdna_est` is
            # `gdna_em_count`, which EXCLUDES fragments that reached no locus; `gdna_true` is the
            # simulator's origin count, which includes them. Scoring one against the other is the
            # documented mistake in `score_library`'s own docstring — it read 0.3151 against a truth of
            # 0.5000 and fabricated an off-capture EM under-call that does not exist. Measured here
            # before the fix: -50.7 % panel-wide, which is the intergenic pool and nothing else.
            # ⭐ The `of which` rows are informational and carry no truth of their own; only the TOTAL
            # is scored, and it is `cli.py`'s own `gdna_fraction` numerator.
            gdna_em, inter = tot("gdna_est"), tot("n_intergenic")
            rows = [
                ("gDNA (total)", gdna_em + inter, tot("gdna_true")),
                ("  of which EM", gdna_em, None),
                ("  of which intergenic", inter, None),
                ("nascent", tot("nrna_est"), tot("nrna_true")),
                ("annotated", tot("mrna_est"), tot("mrna_true")),
            ]
            first = True
            for pretty, est, tru in rows:
                if tru is None:
                    print(f"    {label if first else '':<26} {pretty:<20} {est:>15,.0f} "
                          f"{'—':>15} {'—':>15} {'—':>8}")
                else:
                    d = est - tru
                    # ⛔ a relative error against a TRUE ZERO is not a large number, it is undefined —
                    # printing `inf` or a huge % invites the misreading the header warns of.
                    pct = f"{100.0 * d / tru:>7.1f}%" if tru > 0 else "     n/a"
                    print(f"    {label if first else '':<26} {pretty:<20} {est:>15,.0f} "
                          f"{tru:>15,.0f} {d:>+15,.0f} {pct:>8}")
                first = False

        pool_rows("ALL (g00 excluded)", lambda c: not is_zero_gdna(c))
        for st in _STRATA:
            pool_rows(" x ".join(st),
                      lambda c, st=st: stratum(c) == st and not is_zero_gdna(c))
        print("    " + "-" * 102)
        pool_rows("⛔ g00 ZERO-gDNA control", is_zero_gdna)

    # ── ⑦ THE POOL LEVEL, PER CONDITION ──────────────────────────────────────────────────────────
    #
    # ⭐ The stratum roll-up above answers "which stratum is broken"; this answers "what happened in
    # THIS scenario", which is the row an owner reads when deciding whether a condition is usable. Same
    # three pools, same intergenic-inclusive gDNA, one line per condition.
    for a, name in arms:
        probe = next((a[(c, "library")] for c in conds if (c, "library") in a), None)
        if probe is None or any(f not in probe for f in _POOL_FIELDS):
            continue
        print()
        print(f"  ⑦ PER CONDITION — the three pools, est vs true (fragments)     arm: {name}")
        print(f"    {'condition':<44} {'gDNA est':>13} {'gDNA true':>13} {'nasc est':>11} "
              f"{'nasc true':>10} {'annot est':>13} {'annot true':>13}")
        print("    " + "-" * 122)
        for c in conds:
            r = a[(c, "library")]
            print(f"    {c:<44} {r['gdna_est'] + r['n_intergenic']:>13,.0f} {r['gdna_true']:>13,.0f} "
                  f"{r['nrna_est']:>11,.0f} {r['nrna_true']:>10,.0f} "
                  f"{r['mrna_est']:>13,.0f} {r['mrna_true']:>13,.0f}")
        tg = sum(a[(c, "library")]["gdna_est"] + a[(c, "library")]["n_intergenic"] for c in conds)
        print("    " + "-" * 122)
        print(f"    {'TOTAL (all ' + str(len(conds)) + ' conditions)':<44} {tg:>13,.0f} "
              f"{sum(a[(c, 'library')]['gdna_true'] for c in conds):>13,.0f} "
              f"{sum(a[(c, 'library')]['nrna_est'] for c in conds):>11,.0f} "
              f"{sum(a[(c, 'library')]['nrna_true'] for c in conds):>10,.0f} "
              f"{sum(a[(c, 'library')]['mrna_est'] for c in conds):>13,.0f} "
              f"{sum(a[(c, 'library')]['mrna_true'] for c in conds):>13,.0f}")

    # the library thermometer
    print()
    print("  ⑧ LIBRARY gDNA FRACTION — the thermometer")
    head = f"    {'condition':<44} {'truth':>8}"
    for _a, name in arms:
        head += f" {name[:14]:>14}"
    print(head)
    print("    " + "-" * (53 + 15 * len(arms)))
    for c in conds:
        row = f"    {c:<44} {arms[0][0][(c, 'library')]['gdna_frac_true']:>8.4f}"
        for a, _n in arms:
            row += f" {a[(c, 'library')]['gdna_frac_est']:>14.4f}"
        print(row)

    # per condition, the deliverable
    print()
    print("  ⑨ PER CONDITION — misassigned fragments and false-positive mass")
    head = f"    {'condition':<44}"
    for _a, name in arms:
        head += f" {name[:12] + ' Σ|Δ|':>18}{name[:12] + ' FP':>16}"
    print(head)
    print("    " + "-" * (44 + 34 * len(arms)))
    for c in conds:
        row = f"    {c:<44}"
        for a, _n in arms:
            r = a[(c, "transcript")]
            row += f" {r['count_abs_err']:>17,.0f} {r['fp_mass']:>15,.0f}"
        print(row)


# ── main ─────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", nargs="+", type=Path, default=None,
                    help="print the per-stratum tables from arm jsonl files and exit")
    ap.add_argument("--arm", choices=ARMS, default=None)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--em-seed", type=int, default=DEFAULT_EM_SEED,
                    help="⛔ pinned, because the shipped default is None and the EM's hard "
                         "assignment is an unseeded categorical draw — see DEFAULT_EM_SEED")
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()

    if args.report:
        report(args.report)
        return 0
    if args.arm is None or args.out is None:
        raise SystemExit("--arm and --out are required (or use --report)")

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    cache = args.oracle_cache
    if cache is None and (args.suite / "oracle_cache").is_dir():
        cache = args.suite / "oracle_cache"

    if args.jobs > 1 and len(names) > 1:
        # ⭐ Shards, not threads. Conditions share nothing but a read-only index and cache, so this
        # changes no number; and the measured path stays byte-for-byte the serial one.
        shards = [s for s in (names[i:: args.jobs] for i in range(args.jobs)) if s]
        tmp = args.out.parent / f".{args.out.stem}_shards"
        tmp.mkdir(parents=True, exist_ok=True)
        procs, outs = [], []
        for i, sh in enumerate(shards):
            o = tmp / f"{i}.jsonl"
            outs.append(o)
            cmd = [sys.executable, str(Path(__file__).resolve()), "--arm", args.arm,
                   "--suite", str(args.suite), "--index", str(args.index),
                   "--out", str(o), "--em-seed", str(args.em_seed),
                   "--conditions", *sh]
            if cache is not None:
                cmd += ["--oracle-cache", str(cache)]
            procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT, text=True))
        rc = 0
        for i, pr in enumerate(procs):
            out, _ = pr.communicate()
            if pr.returncode != 0:
                rc = pr.returncode
                print(f"  ⛔ shard {i} FAILED (rc={pr.returncode}):\n{out}", flush=True)
            else:
                print(f"  shard {i}: {len(shards[i])} conditions ok", flush=True)
        if rc:
            # ⛔ A short output file reads as a complete panel (TRAPS: an-ablation-that-never-ran's shape).
            raise SystemExit("a shard failed; refusing to concatenate a partial panel")
        with args.out.open("w") as fh:
            for o in outs:
                fh.write(o.read_text())
        print(f"  ⭐ {args.arm}: {sum(1 for _ in args.out.open())} rows -> {args.out}")
        return 0

    index = TranscriptIndex.load(str(args.index))
    pipeline_config = PipelineConfig()
    rows = []
    for name in names:
        print(f"  … {args.arm}  {name}", flush=True)
        rows += run_condition(args.arm, args.suite, index, name, pipeline_config, cache,
                              em_seed=args.em_seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"  ⭐ {args.arm}: {len(rows)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
