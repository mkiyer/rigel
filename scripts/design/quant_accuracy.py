#!/usr/bin/env python
"""How accurate is the tool end to end, and what is a perfect prior worth?

``--arm base`` runs the shipped pipeline on a simulated condition and scores its transcript table
against the simulator's own per-transcript truth. Every other arm runs the identical pipeline with
one thing substituted -- the oracle ``LocusPriors`` built from the origin-split truth (``oracle``, or
one of its two arrays alone: ``oracle_gdna``, ``oracle_efflen``), the ruler the EM
divides by (``oracle_ruler``: the SIMULATOR's capture-aware effective length in place of the shipped
ruler's, the only arm that reaches the lengths the EM divides by), the EM's seed (``warm_uniform``), or
the per-transcript
allocation weights (the ``oracle_alloc*`` arms, a capability proof and never a headroom claim) -- and
the difference from ``base`` is what that one thing is worth. One scorer serves every arm, so the
ceiling and the baseline cannot drift apart. The primary score is count against count: the truth is
each condition's realised observed fragment count, never the pre-capture molar abundance, and no
length model enters the comparison; the TPM rows share the tool's own effective length on both sides
and so measure assignment only. This is a thermometer above the 0.8.0 metric, never the target.
Reading rules: ``noop`` and ``oracle_ruler_noop`` must be byte-identical to ``base``
(gated in ``tests/calibration/test_quant_accuracy.py``); ``base_reseed`` is the sampling noise floor and any smaller delta
is noise; the oracle masses are undrained while the shipped pipeline drains, a small conservative
bias that cannot explain a large surviving error or hide a large removed one. The library-level
gDNA fraction counts intergenic fragments as gDNA, as ``cli.py`` does.

Gates: ``tests/calibration/test_quant_accuracy.py``.

Usage::

    python scripts/design/quant_accuracy.py --arm base --out $RIGEL_ARMS/qa_base.jsonl --jobs 4
    python scripts/design/quant_accuracy.py --arm oracle --oracle-cache DIR --out $RIGEL_ARMS/qa_oracle.jsonl
    python scripts/design/quant_accuracy.py --arm base --conditions COND --em-seed 1 --suite DIR --index INDEX
    python scripts/design/quant_accuracy.py --report $RIGEL_ARMS/qa_base.jsonl $RIGEL_ARMS/qa_oracle.jsonl
    python scripts/design/quant_accuracy.py --arm base --set em.assignment_mode=fractional --out F   # a config arm
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
sys.path.insert(0, str(_REPO / "scripts" / "design"))

from _oracle import ORIGINS, OracleTruth  # noqa: E402
from _shared import set_field  # noqa: E402


import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.config import EMConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, run_pipeline  # noqa: E402
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: ``warm_uniform`` seeds every component equally instead of by coverage-weighted share. It varies
#: one thing against ``base`` — where the EM is put down — so a difference is the seed's basin and
#: nothing else.
#:
#: ``oracle_alloc`` is a capability proof rather than a ceiling. It hands the EM the true relative
#: transcript abundances as the per-transcript allocation weights and zeroes the coverage seed, so
#: `theta` starts from the prior alone. The question is binary: given correct weights, does the
#: machinery produce correct per-transcript abundances? It is not a claim about reachable headroom:
#: truth-derived weights also hand over the true support (a zero weight is exactly absorbing, so every
#: silent transcript is switched off for free), and pricing what a real weighting function could earn
#: needs controls this arm deliberately omits.
#:
#: The ruler arms substitute the lengths the EM divides by. A transcript's capture-aware effective length
#: (``calibration.capture_eff_length.transcript_capture_eff_lengths``, the "ruler") is built by
#: ``pipeline._setup_geometry_and_estimator`` BEFORE ``assemble_priors`` runs, so no prior arm reaches it::
#:
#:     calibrate(...)
#:       |-- transcript_capture_eff_lengths()  <- the ruler arms substitute HERE: `effective_lengths_em`
#:       |-- assemble_priors()                 <- every prior arm wraps this: `LocusPriors`
#:
#: ``oracle_ruler`` hands the EM ``fl × factor``, the factor being the SIMULATOR's own capture factor per
#: transcript (``ruler_vs_truth.load_truth``: the sampler's partition, the truth the reads were drawn with),
#: anchored so the fully probed class reads 1 — the fully captured level the shipped efficiencies are read
#: against (``anchor_capture_factor``) — and never clipped: the anchor is a median, so the fully probed class
#: straddles 1. It varies one thing against ``base``: the priors, and the locus gDNA component's length with
#: them, stay as shipped, so the anchor also sets the transcripts' scale beside that component. Under capture a panel's probes can make capture ISOFORM-specific (a
#: probe across a junction captures only the isoforms holding it), and the split between isoforms is decided
#: by the ratio of their lengths, so this arm prices the ruler end to end; capture-OFF, the truth is the plain
#: length and the arm is the identity. It replaced an arm that swapped calibration's count arrays, which the
#: ruler stopped reading at ``c44fc306`` and so could not move it by one ulp.
#:
_RULER_ARMS = {"oracle_ruler": True, "oracle_ruler_noop": False}

ARMS = ("base", "base_reseed", "noop", "oracle", "oracle_gdna", "oracle_efflen",
        "warm_uniform", "oracle_alloc", "oracle_alloc_seed", "oracle_alloc_flip") + tuple(_RULER_ARMS)

#: The EM seed every arm pins: the shipped one, so ``base`` is the configuration that ships.
#: ``base_reseed`` re-runs ``base`` at ``seed + 1`` and is the noise floor. The seed reaches only the
#: ``sample`` assignment's draw, so under ``--set em.assignment_mode=fractional`` — how every arm is
#: benchmarked — the two differ by nothing but what varies from one run to the next.
DEFAULT_EM_SEED = EMConfig().seed

#: arm -> which ``LocusPriors`` fields come from O. ``noop`` takes NONE of them and still builds O,
#: which is what makes it a test of the wrapper rather than of an ``if``. There is no RNA arm:
#: calibration's RNA count does not reach the EM (``pipeline.em_pseudocounts`` reads the gDNA count
#: against the EM's own count of the locus), so an arm injecting it could not fire
#: (TRAPS: an-ablation-that-never-ran).
_ARM_FIELDS = {
    "base": (),
    "base_reseed": (),
    "noop": (),
    "oracle": ("gdna_count", "gdna_eff_len"),
    "oracle_gdna": ("gdna_count",),
    "oracle_efflen": ("gdna_eff_len",),
}


# ── the injection ────────────────────────────────────────────────────────────────────────────────


def load_oracle(bam: str, index, pipeline_config, cache_root: Path, tag: str) -> OracleTruth:
    """The origin-split truth for one condition, entirely from the shipped scan cache.

    ``_main`` (the UNDRAINED full payload) is what sum-to-full is asserted against, and it must come
    from the same cache as the three partitions or the identity is checking two different scans
    against each other. ``read_scan_cache`` refuses a payload whose ``graph_hash`` / ``reach_digest``
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
            "calibration_oracle.py --build (panel.py cache runs it); this script "
            "refuses to invent a truth."
        ) from exc
    return OracleTruth.from_parts(full, parts)


def install_arm(arm: str, oracle: OracleTruth | None):
    """Wrap :func:`~rigel.calibration.priors.assemble_priors` for one arm. Returns ``(restore, fired)``.

    The wrapper is IDENTICAL across every arm including ``noop`` — the oracle is read, O is
    built on the loci the run itself produced, and only then does the arm decide which fields to take.
    A ``noop`` that short-circuits before building O would prove the ``if`` works and nothing else
    (TRAPS: could-the-arm-have-fired); this one proves the whole path is inert when it takes no field.

    O is built on the run's OWN ``multi_loci``, never on a stored array. ``build_multi_loci``
    unions transcripts linked by SCORED fragments, so the locus partition is a function of the run —
    an oracle prior keyed by ``multi_locus_id`` from another process is not portable and index-aligning
    it would compare locus 7 of one run with locus 7 of another.

    ``fired`` is a mutable counter, not a bool, and the caller RAISES on zero: an override that never
    ran reads as "no effect", which is the most flattering possible failure (TRAPS: an-ablation-that-never-ran).
    """
    original = PRIORS.assemble_priors

    if arm in ("base", "base_reseed", "warm_uniform"):
        # Counted as fired: ``base`` installs nothing by design, so the "did the override run?"
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

    ``count`` against ``mrna_abundance``, both fragment counts. No effective length, no
    normalisation, no model between the two sides — so a difference is a fragment the tool put
    somewhere the simulator did not.

    The false-positive mass is reported separately and it is the number this tool exists for.
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

    # TPM. Truth TPM uses the TOOL's own effective length, so the length model is common to both
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
        #: contamination the deconvolution left on a silent transcript — a pure false positive
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

    This is the discriminator that separates isoform ambiguity from something else, and without it
    the transcript number cannot be acted on. Summing a gene's isoforms collapses exactly
    the error that comes from not knowing WHICH isoform a fragment came from — the difficulty every
    transcript quantifier has and that no prior can remove. What survives at gene level is error in
    deciding whether the fragment was RNA *from this gene* at all, which is the question Rigel is for.

    Grouped on ``gene_id`` from the TRUTH table, not from the index, so the two sides are grouped by
    one definition. Synthetic nRNA entities are absent from ``quant`` by construction
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

    Intergenic fragments are gDNA and the denominator must say so. ``n_intergenic`` counts
    fragments that reached no locus at all; they never enter the EM, so ``estimator.gdna_em_count``
    excludes them. Off capture, gDNA is genome-uniform and more than half of it is intergenic, so
    leaving them out of the numerator while RNA stays in the denominator understates the fraction
    badly. This is ``cli.py``'s ``gdna_fraction`` — ``(gdna_em + n_intergenic) / (rna + gdna_em +
    n_intergenic)`` — so the number here is the one a user reads. ``gdna_em_frac_est`` is kept
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
        # THE NASCENT ARM IS BROKEN OUT AND IT IS NOT COSMETIC. ``get_counts_df`` drops the
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
        #: THE DELIVERABLE — intergenic included, cli.py's denominator.
        "gdna_frac_est": (gdna + intergenic) / total_all if total_all > 0 else float("nan"),
        #: the EM's own view, intergenic in NEITHER side. A different question; not the deliverable.
        "gdna_em_frac_est": gdna / total if total > 0 else float("nan"),
        "gdna_frac_true": t_g / (t_g + t_r) if (t_g + t_r) > 0 else float("nan"),
        "n_intergenic": intergenic,
    }


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def seeded(pipeline_config, arm: str, em_seed: int):
    """The arm's pipeline config. ``base_reseed`` differs from ``base`` in the SEED ALONE, so the
    gap between them is the sampling noise of the EM's own hard assignment and nothing else.

    """
    seed = em_seed + 1 if arm == "base_reseed" else em_seed
    warm = pipeline_config.em.warm_start
    if arm == "warm_uniform":
        warm = "uniform"
    elif arm == "oracle_alloc":
        # the seed is zeroed so `theta` starts proportional to the prior alone — otherwise a
        # coverage-weighted seed is multiplied by an allocation from a different method and the result
        # is neither method's answer.
        warm = "prior"
    out = dataclasses.replace(
        pipeline_config, em=dataclasses.replace(pipeline_config.em, seed=seed, warm_start=warm)
    )
    return out


#: The truth table's RNA fragment columns, most realised first. Each entry is the MATURE and the
#: NASCENT column together, and that pairing is the point: a synthetic shadow entity emits no mature
#: RNA, so its whole realised mass sits in the nascent column of its OWN row, while a single-exon
#: annotated transcript carries both. A mature column read alone weights every shadow at zero.
_WEIGHT_COLUMNS = (
    ("observed_mrna_fragments", "observed_nrna_fragments"),
    ("mrna_abundance", "nrna_abundance"),
)


def truth_weights(truth: pd.DataFrame, index) -> np.ndarray:
    """``float64[n_transcripts]`` — the TRUE realised RNA fragment count per transcript, on the EM's axis.

    This is the whole of ``oracle_alloc``'s input: the relative abundances the simulator actually produced.
    Within a locus the EM only reads their RATIOS, so no normalisation is needed here.

    ⛔ THE WEIGHT IS MATURE **AND** NASCENT, because the RNA prior's recipients are every RNA component
    and a synthetic shadow entity is one of them (`EQUATIONS.md` §9b). A shadow's realised fragments are
    the NASCENT column of its own row and its mature column is identically zero, so reading the mature
    column alone hands every shadow a weight of exactly zero — which is not an allocation but the
    retired ``alpha = 0`` rule (`ISSUES: nascent-gets-no-rna-prior`, CLOSED) under an oracle's name, and
    it reads as a perfect allocation removing the nascent over-call. Gated by
    ``tests/calibration/test_quant_accuracy.py``. A table with no nascent column cannot express a
    shadow's allocation and is refused rather than silently weighted at zero.

    Exact-duplicate transcripts are folded onto the twin the index kept — the truth table is keyed on
    the un-collapsed annotation, so without the fold their fragments would be dropped rather than
    attributed. For such a pair the per-transcript truth is not merely awkward, it is UNDEFINED: the
    two are the same molecule and only the group total is a fact about the world.
    """
    pair = next((p for p in _WEIGHT_COLUMNS if all(c in truth.columns for c in p)), None)
    if pair is None:
        raise SystemExit(
            "⛔ truth table carries no mature+nascent RNA pair "
            f"({' or '.join('+'.join(p) for p in _WEIGHT_COLUMNS)}); a mature column alone cannot "
            "express a synthetic nascent entity's allocation"
        )
    total = truth[pair[0]].to_numpy(np.float64) + truth[pair[1]].to_numpy(np.float64)
    t_index = dict(zip(index.t_df["t_id"].to_numpy(), index.t_df["t_index"].to_numpy(), strict=True))
    w = np.zeros(int(index.num_transcripts), dtype=np.float64)
    for tid, n in zip(truth["transcript_id"], total, strict=True):
        i = t_index.get(str(tid))
        if i is not None:
            w[int(i)] += float(n)
    return w


def anchor_capture_factor(factor, probed_frac, plain_length) -> np.ndarray:
    """The simulator's capture factor on the reference's scale — a fully captured object reads 1: divided by
    the median over the fully probed class (probed fraction ≥ 0.9), every transcript when no such class exists (capture-OFF),
    and the identity where the factor is undefined (a transcript with no plain length). The anchor is the
    one ``ruler_vs_truth.score`` reads its errors against."""
    f = np.asarray(factor, dtype=np.float64)
    ok = np.isfinite(f) & (f > 0.0) & (np.asarray(plain_length, dtype=np.float64) > 0.0)
    logf = np.where(ok, np.log(np.where(ok, f, 1.0)), 0.0)
    probed = ok & (np.asarray(probed_frac, dtype=np.float64) >= 0.9)
    pool = probed if probed.any() else ok
    scale = float(np.median(logf[pool])) if pool.any() else 0.0
    return np.where(ok, np.exp(logf - scale), 1.0)


def capture_truth_factor(suite: Path, index, condition: str) -> np.ndarray:
    """``float64[n_transcripts]`` — the anchored capture factor on the index's axis, for the condition's
    capture label. It depends on the panel's capture config, the probe file, the simulated fragment-length
    law and the annotation, never on the condition, so it is computed once per label and cached at
    ``<suite>/oracle_cache/capture_truth_<label>.npz`` under a key over all four; a different key is
    recomputed, never read. Computing it integrates the sampler over every fragment length, minutes on the
    ladder."""
    import hashlib

    label = "on" if condition.endswith("_capture_on") else "off"
    manifest = json.loads((Path(suite) / "manifest.json").read_text())
    capture = next(c["config"] for c in manifest["capture_configs"] if c["label"] == label)
    probes = capture.get("probes")
    key_parts = {
        "capture": capture,
        "probes_sha": hashlib.sha256(Path(probes).read_bytes()).hexdigest() if probes else None,
        "fl": {k: manifest["simulation"][k] for k in ("frag_mean", "frag_std", "frag_min", "frag_max")},
        "t_id": [str(x) for x in index.t_df["t_id"]],
        "length": [int(x) for x in index.t_df["length"]],
    }
    key = hashlib.sha256(json.dumps(key_parts, sort_keys=True, default=str).encode()).hexdigest()[:16]
    path = Path(suite) / "oracle_cache" / f"capture_truth_{label}.npz"
    if path.exists():
        with np.load(path) as z:
            if str(z["key"]) == key:
                return anchor_capture_factor(z["factor"], z["probed_frac"], z["L_plain"])
    from _shared import sibling

    truth = sibling("ruler_vs_truth.py").load_truth(index, Path(index.index_dir), Path(suite), condition)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.stem}.{os.getpid()}.npz")
    np.savez(tmp, key=key, factor=truth.factor, probed_frac=truth.probed_frac, L_plain=truth.L_plain)
    tmp.replace(path)
    return anchor_capture_factor(truth.factor, truth.probed_frac, truth.L_plain)


def install_ruler_arm(arm: str, factor: np.ndarray):
    """Wrap the ruler so the EM divides by ``fl × factor``. Returns ``(restore, fired)``.

    ``transcript_capture_eff_lengths`` is patched as a MODULE ATTRIBUTE, which works because
    ``pipeline._setup_geometry_and_estimator`` imports it function-locally, at call time.

    The wrapper is identical for both arms: it runs the shipped ruler, builds the substitute, and records
    how far the substitute sits from the shipped lengths (``truth_delta``) and how far what it RETURNED
    sits from them (``max_abs_delta``). ``oracle_ruler_noop`` then hands back the shipped lengths' own
    object, so it is a test of the whole wrapper and not of an ``if`` (TRAPS: could-the-arm-have-fired).
    ``check_ruler_arm`` reads the counter: the arm must move the lengths wherever the truth differs from
    them, the noop must not move them at all.
    """
    from rigel.calibration import capture_eff_length as CEL

    substitute = _RULER_ARMS[arm]
    factor = np.asarray(factor, dtype=np.float64)
    orig_ruler = CEL.transcript_capture_eff_lengths
    fired = {"n": 0, "ruler": 0, "max_abs_delta": 0.0, "truth_delta": 0.0}

    def ruler_wrapper(calibration, region_arrays, index, fl_eff_lengths, rna_fl_pmf):
        shipped = orig_ruler(calibration, region_arrays, index, fl_eff_lengths, rna_fl_pmf)
        base = np.asarray(shipped, dtype=np.float64)
        truth = np.asarray(fl_eff_lengths, dtype=np.float64) * factor
        out = truth if substitute else shipped
        fired["n"] += 1
        fired["ruler"] += 1
        fired["truth_delta"] = max(fired["truth_delta"], float(np.abs(truth - base).max()))
        fired["max_abs_delta"] = max(
            fired["max_abs_delta"], float(np.abs(np.asarray(out, dtype=np.float64) - base).max())
        )
        return out

    CEL.transcript_capture_eff_lengths = ruler_wrapper

    def restore():
        CEL.transcript_capture_eff_lengths = orig_ruler

    return restore, fired


def check_ruler_arm(condition: str, arm: str, fired: dict) -> None:
    """Refuse a ruler arm that could not have measured what it names. The shrinkage must have run; the
    substituting arm must move the EM's lengths wherever the truth differs from the shipped ones (where it
    does not — capture-OFF — the arm is the identity and that is its answer); the noop must move nothing."""
    if fired["ruler"] == 0:
        raise RuntimeError(
            f"{condition} [{arm}]: the ruler never ran, so this arm substituted nothing the EM read."
        )
    if _RULER_ARMS[arm] and fired["truth_delta"] > 0.0 and fired["max_abs_delta"] == 0.0:
        raise RuntimeError(
            f"{condition} [{arm}]: the truth differs from the shipped ruler by {fired['truth_delta']:.3g} "
            "but the EM's lengths did not move — an arm that cannot move the quantity it names has not "
            "measured it."
        )
    if not _RULER_ARMS[arm] and fired["max_abs_delta"] != 0.0:
        raise RuntimeError(
            f"{condition} [{arm}]: the NOOP arm moved the EM's lengths by {fired['max_abs_delta']:.3e}."
        )


def install_truth_weights(weights: np.ndarray):
    """Pass per-transcript weights into the EM, and COUNT THE ARM AT THE DEEPEST POINT IT CAN OBSERVE.

    A wrapper rather than a config field: the weights are a per-run ARRAY, and an experiment that
    injects one should not add production surface.

    The counter watches the estimator, not the injection: counting nonzero entries in ``weights``
    would be true whatever the pipeline then did with the array, and a parameter dropped between the
    wrapper and the solver would read healthy. Counting where the solver receives the array makes a
    dropped parameter impossible to miss, and ``--arm oracle_alloc_flip`` (a maximally wrong
    allocation, which must move the answer) is the end-to-end half of the same check.
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
    if not (arm in ("base", "base_reseed", "warm_uniform") or arm.startswith("oracle_alloc")
            or arm in _RULER_ARMS):
        if oracle_cache is None:
            raise SystemExit(f"⛔ arm {arm!r} needs --oracle-cache")
        oracle = load_oracle(bam, index, pipeline_config, oracle_cache, condition)

    if arm == "oracle_alloc_flip":
        w = truth_weights(truth, index)
        # falsification of the harness, not a treatment: put the weight where the truth is NOT.
        # If a maximally wrong allocation moves nothing, the allocation never reached the solver.
        flip = np.zeros_like(w)
        nz = np.flatnonzero(w > 0)
        if nz.size:
            flip[nz] = w[nz][::-1]
        restore, fired = install_truth_weights(flip)
    elif arm.startswith("oracle_alloc"):
        # `oracle_alloc_seed` keeps the shipped coverage seed, so it varies one thing against `base`:
        # the allocation. `oracle_alloc` also zeroes the seed, which additionally removes RNA's evidence
        # advantage over calibration's gDNA prior — two changes, not one.
        restore, fired = install_truth_weights(truth_weights(truth, index))
    elif arm in _RULER_ARMS:
        # exact membership, and placed before the fall-through: `oracle_ruler` starts with "oracle"
        # and would otherwise land in `install_arm`, whose `_ARM_FIELDS[arm]` would raise — or worse,
        # would not have, had the name been one letter different.
        restore, fired = install_ruler_arm(arm, capture_truth_factor(suite, index, condition))
    else:
        restore, fired = install_arm(arm, oracle)
    start = time.perf_counter()
    try:
        result = run_pipeline(bam, index, pipeline_config)
    finally:
        restore()
    if fired["n"] == 0:
        # TRAPS: an-ablation-that-never-ran — an injection that never ran reads as "a perfect prior changes nothing".
        raise RuntimeError(
            f"{condition} [{arm}]: the override was never wrapped-and-called. This is not a "
            "measurement of zero effect."
        )
    if arm in _RULER_ARMS:
        # the ruler arms carry a second requirement: the EM's lengths are what they name (check_ruler_arm)
        check_ruler_arm(condition, arm, fired)
    seconds = time.perf_counter() - start

    quant = result.estimator.get_counts_df(index)
    common = {"arm": arm, "condition": condition, "seconds": seconds,
              "em_seed": int(pipeline_config.em.seed),
              # stamped on every row, so a fractional run cannot be reported as a sampled one
              "assignment_mode": pipeline_config.em.assignment_mode}
    if arm in _RULER_ARMS:
        # How far the arm moved the EM's lengths, and how far the truth sits from the shipped ruler, in
        # base pairs of opportunity on the worst transcript — recorded beside the score, so the arm's
        # REACH is readable off the output file.
        common["ruler_max_abs_delta"] = float(fired["max_abs_delta"])
        common["ruler_truth_delta"] = float(fired["truth_delta"])
    return [
        {**common, "axis": "transcript", **score_transcripts(quant, truth)},
        # the SAME scorer over genes — isoform ambiguity summed away, see score_genes
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


_GDNA_LEVEL = {"g00": "0 %", "g05": "5 %", "g50": "50 %", "g98": "98 %"}


def _level(cond: str) -> str:
    """The condition's designed gDNA level, from its name."""
    for k in _GDNA_LEVEL:
        if f"_{k}_" in cond:
            return k
    return "?"


def _pct(est: float, true: float) -> str:
    """Percent error against the realised truth. ⛔ UNDEFINED at a truth of zero, and printed as
    such rather than as a large number: at `g00` the true gDNA count IS zero, and a ratio there
    would invent a scale (`TRAPS: a-ratio-cannot-carry-zero`)."""
    if true <= 0.0:
        return "n/a"
    return f"{100.0 * (est - true) / true:+,.1f} %"


def _signed(x: float) -> str:
    return f"{x:+,.0f}"


def markdown_report(paths: list[Path], out: Path) -> None:
    """The full per-scenario accuracy report, as markdown — what a release is judged on.

    THREE POOLS AND THEY ARE NOT THE SAME QUESTION. Every fragment in the library is gDNA, SYNTHETIC
    nascent RNA, or ANNOTATED RNA, and the split is on ``is_synthetic`` — a span this index
    manufactured — never on ``is_nrna``. A single-exon ANNOTATED transcript carries ``is_nrna`` because
    it is at once the nascent and the mature form of a real gene, and it belongs to the ANNOTATED pool
    (`TRAPS: nrna-does-not-mean-synthetic`). The transcript and gene sections then score the annotated
    pool ALONE, because that is the table a user reads: synthetic entities are absent from it by
    construction and their mass is reported on the pool rows instead.

    ⛔ The gDNA estimate is ``gdna_em + n_intergenic``. Intergenic fragments reach no locus, so they
    never enter the EM, but they ARE gDNA and the truth counts them — comparing the EM's number alone
    against that truth would understate the estimate by more than half at capture-OFF.

    The first path is the arm reported; an arm whose stem contains ``reseed`` is used as the
    ATTRIBUTION FLOOR and printed beside every transcript row, because no delta below it is
    attributable (`TRAPS: the-deliverable-is-not-reproducible-by-default`).
    """
    arms = [(_load(p), Path(p).stem) for p in paths]
    primary, pname = arms[0]
    floor = next((a for a, n in arms if "reseed" in n), None)

    modes = {r.get("assignment_mode") for r in primary.values()}
    conds = sorted({c for c, _ax in primary})
    stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime(Path(paths[0]).stat().st_mtime))
    from rigel.config import CalibrationConfig as _CC

    L = []
    w = L.append
    w("# Rigel — end-to-end accuracy on the gDNA ladder")
    w("")
    w(f"**Arm** `{pname}` · **{len(conds)} conditions** · **assignment** "
      f"`{'/'.join(sorted(str(m) for m in modes))}` · **calibration** "
      f"`message_policy={_CC().message_policy!r}` · **scored** {stamp}")
    w("")
    w("Every number is scored against the simulator's own REALISED per-fragment truth — each fragment "
      "counted exactly once, never a pre-capture molar abundance — so a difference is a fragment the "
      "tool put somewhere the simulator did not.")
    w("")
    w("## How to read this")
    w("")
    w("- **Three pools, and they answer different questions.** `gDNA` is contamination the "
      "deconvolution must remove; `nascent RNA` is the SYNTHETIC entity this index manufactures for "
      "each multi-exon gene; `annotated RNA` is the transcript table a user reads. The split is on "
      "`is_synthetic`, never `is_nrna` — a single-exon annotated transcript is ANNOTATED RNA.")
    w("- **The gDNA estimate includes intergenic fragments**, which reach no locus and never enter "
      "the EM. Off capture they are more than half of all gDNA.")
    w("- **Percent error is undefined where the truth is zero** and is printed `n/a`. At `g00` the "
      "true gDNA count is exactly 0; the raw count beside it is the whole of the answer there.")
    w("- **Read per stratum, never pooled.** Three strata are in scope. `unstranded × capture ON` is "
      "DEFERRED — reported on every benchmark, never a development target — and it carries most of "
      "the error, so a pooled total would be its total.")
    w("- **Nothing below the attribution floor is attributable.** The floor is the same arm re-run "
      "under the next EM seed; it is printed beside every transcript row. Under the fractional "
      "assignment every arm runs, the seed reaches no number, so the floor is the run-to-run spread "
      "itself: the BAM scan's workers sum the tally's fractions in whatever batches each took, and the "
      "EM carries that last-bit difference into whole fragments — four runs of `g50 ss.99 OFF` spanned "
      "145 fragments (`TRAPS: the-deliverable-is-not-reproducible-by-default`). Treat every figure here "
      "as carrying that much noise.")
    w("- **`expressed` and `detected` are the scored sets, and the truth table is larger than "
      "either.** It carries one row per SYNTHETIC nascent entity as well, and those rows are zero on "
      "both sides — zero truth and zero estimate, since the transcript table drops them — so they "
      "enter no figure here. A raw row count would read as thousands of transcripts scored perfectly; "
      "it is not reported for that reason.")
    w("- **There is no false-negative column, because under fractional assignment it cannot fire.** "
      "A false negative needs an estimate of EXACTLY zero, and a fractional posterior essentially "
      "never is; it reads 0 on all 16 conditions and would look like a perfect score for something "
      "unmeasured. `under-assigned` is the live quantity and it is reported instead.")
    w("- **`MARD` is the SYMMETRIC mean absolute relative difference**, `|est−true| / (|est|+|true|)` "
      "averaged over the expressed set, so it is bounded in [0, 1] and finite at a zero estimate — "
      "0.5 is a 3× error, not a 50 % one. `median rel. err` is the ordinary `|est−true|/true`, over "
      "the same set. Both are unweighted by mass: a 1-fragment transcript counts as much as a "
      "25,000-fragment one.")
    w("- **One asymmetry in the pool rows.** The ESTIMATE splits on `is_synthetic` (the entity) while "
      "the TRUTH splits on the simulator's template kind, and a handful of annotated single-exon "
      "transcripts serve as their own nascent entity — their nascent fragments count in `nascent "
      "truth` while the tool's counts for them land in `annotated`. It is ~300 fragments of ~1.9 M "
      "and it cannot be removed by relabelling; it is recorded so the pool rows are not read as one "
      "partition measured twice.")
    w("- **Every stratum row includes its `g00` rung.** The terminal report excludes the zero-gDNA "
      "control from its roll-up; this one keeps it, so the two do not agree by construction.")
    w("")

    # ── pools ────────────────────────────────────────────────────────────────────────────────────
    w("## 1. Pool level — where the library's fragments went")
    w("")
    w("One row per scenario, grouped by stratum. `Δ` is `estimate − truth` in fragments.")
    w("")
    for pool, est_f, true_f, label in (
        ("gDNA", None, "gdna_true", "gDNA (EM + intergenic)"),
        ("nascent", "nrna_est", "nrna_true", "Nascent RNA — SYNTHETIC entities only"),
        ("annotated", "mrna_est", "mrna_true", "Annotated RNA — the transcript table"),
    ):
        w(f"### {label}")
        w("")
        w("| scenario | gDNA level | estimated | truth | Δ | % error |")
        w("|---|---|---:|---:|---:|---:|")
        for st in _STRATA:
            rows = [c for c in conds if stratum(c) == st]
            if not rows:
                continue
            tag = f"{st[0]} × {st[1]}" + ("  ⛔ DEFERRED" if st == ("unstranded", "capture ON") else "")
            w(f"| **{tag}** | | | | | |")
            for c in rows:
                r = primary[(c, "library")]
                est = (r["gdna_est"] + r["n_intergenic"]) if est_f is None else r[est_f]
                true = r[true_f]
                w(f"| `{c}` | {_GDNA_LEVEL[_level(c)]} | {est:,.0f} | {true:,.0f} | "
                  f"{_signed(est - true)} | {_pct(est, true)} |")
        w("")

    w("### Library gDNA fraction — the thermometer")
    w("")
    w("What `rigel quant` reports as the library's gDNA share. This is the number calibration exists "
      "to produce; the transcript table does not always keep it.")
    w("")
    w("| scenario | estimated | truth | Δ |")
    w("|---|---:|---:|---:|")
    for st in _STRATA:
        rows = [c for c in conds if stratum(c) == st]
        if not rows:
            continue
        w(f"| **{st[0]} × {st[1]}** | | | |")
        for c in rows:
            r = primary[(c, "library")]
            e, t = r["gdna_frac_est"], r["gdna_frac_true"]
            w(f"| `{c}` | {e:.4f} | {t:.4f} | {e - t:+.4f} |")
    w("")

    # ── transcript and gene ──────────────────────────────────────────────────────────────────────
    for axis, title, blurb in (
        ("transcript", "2. Transcript level — inside the annotated RNA pool",
         "Scored over the annotated transcripts alone: gDNA and the synthetic nascent entities are "
         "excluded, so what is left is how well the tool splits the RNA it kept. `Σ|Δ|` sums "
         "`|estimate − truth|` over every annotated transcript, SILENT ONES INCLUDED, and splits into "
         "`over-assigned` + `under-assigned`; `false-positive mass` is the part of the over-assignment "
         "that landed on a transcript the simulator gave ZERO fragments, so it is a subset of `Σ|Δ|` "
         "and the number this tool exists for."),
        ("gene", "3. Gene level — the same, with isoform ambiguity summed away",
         "The same scorer over genes, and the difference from the transcript row is the point: "
         "summing a gene's isoforms collapses exactly the error that comes from not knowing WHICH "
         "isoform a fragment came from. What survives is error in deciding whether the fragment was "
         "RNA from this gene at all, which is the question Rigel is for."),
    ):
        w(f"## {title}")
        w("")
        w(blurb)
        w("")
        unit = "expressed" if axis == "transcript" else "expressed genes"
        head = (f"| scenario | {unit} | detected | Σ\\|Δ\\| | Σ\\|Δ\\| as % of true | net Δ | "
                "over-assigned | under-assigned | false-positive mass | on n | MARD | Spearman |")
        if axis == "transcript":
            head = head.replace("| MARD |", "| median rel. err | MARD |")
        if floor is not None and axis == "transcript":
            head = head.replace("| net Δ |", "| net Δ | seed floor |")
        # ⛔ Count columns with the ESCAPED pipes removed. `Σ\|Δ\|` carries two literal `|`
        # characters that are cell CONTENT, not delimiters, and counting them put three phantom
        # columns in every separator row.
        ncols = head.replace("\\|", "").count("|") - 1
        w(head)
        w("|---" + "|---:" * (ncols - 1) + "|")
        for st in _STRATA:
            rows = [c for c in conds if stratum(c) == st]
            if not rows:
                continue
            tag = f"{st[0]} × {st[1]}" + ("  ⛔ DEFERRED" if st == ("unstranded", "capture ON") else "")
            w(f"| **{tag}** |" + " |" * (ncols - 1))
            for c in rows:
                r = primary[(c, axis)]
                share = 100.0 * r["count_abs_err"] / r["count_true"] if r["count_true"] > 0 else float("nan")
                cells = [f"`{c}`", f"{r['n_expressed']:,}", f"{r['n_detected']:,}",
                         f"{r['count_abs_err']:,.0f}", f"{share:.2f} %", _signed(r["count_net_err"])]
                if floor is not None and axis == "transcript":
                    f_r = floor.get((c, axis))
                    cells.append(f"±{abs(r['count_abs_err'] - f_r['count_abs_err']):,.0f}"
                                 if f_r else "—")
                cells += [f"{r['count_over']:,.0f}", f"{r['count_under']:,.0f}",
                          f"{r['fp_mass']:,.0f}", f"{r['fp_n']:,}"]
                if axis == "transcript":
                    cells.append(f"{r['median_rel_err']:.3f}")
                cells += [f"{r['mard']:.3f}", f"{r['spearman']:.4f}"]
                w("| " + " | ".join(cells) + " |")
        w("")

    # ── rollup ───────────────────────────────────────────────────────────────────────────────────
    w("## 4. Per stratum, summed over its four gDNA levels")
    w("")
    w("⛔ Summed WITHIN a stratum only. The four strata are never added together.")
    w("")
    w("| stratum | annotated truth | transcript Σ\\|Δ\\| | % | gene Σ\\|Δ\\| | % | nascent est | "
      "nascent truth | gDNA est | gDNA truth |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for st in _STRATA:
        rows = [c for c in conds if stratum(c) == st]
        if not rows:
            continue
        tx = sum(primary[(c, "transcript")]["count_abs_err"] for c in rows)
        gn = sum(primary[(c, "gene")]["count_abs_err"] for c in rows)
        tt = sum(primary[(c, "transcript")]["count_true"] for c in rows)
        ne = sum(primary[(c, "library")]["nrna_est"] for c in rows)
        nt = sum(primary[(c, "library")]["nrna_true"] for c in rows)
        ge = sum(primary[(c, "library")]["gdna_est"] + primary[(c, "library")]["n_intergenic"] for c in rows)
        gt = sum(primary[(c, "library")]["gdna_true"] for c in rows)
        tag = f"{st[0]} × {st[1]}" + (" ⛔ DEFERRED" if st == ("unstranded", "capture ON") else "")
        w(f"| {tag} | {tt:,.0f} | {tx:,.0f} | {100 * tx / tt:.2f} % | {gn:,.0f} | "
          f"{100 * gn / tt:.2f} % | {ne:,.0f} | {nt:,.0f} | {ge:,.0f} | {gt:,.0f} |")
    w("")
    w("---")
    w("")
    w(f"Generated by `scripts/design/quant_accuracy.py --markdown` from `{pname}.jsonl`"
      + (" with the seed floor from the `reseed` arm." if floor is not None else "."))
    w("")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L))
    print(f"  ⭐ markdown report -> {out}  ({len(L)} lines, {len(conds)} conditions)")


def report(paths: list[Path]) -> None:
    """One or more arms, per stratum. Never pooled — the panel total hides a sign flip between
    strata, and on this panel one stratum carries almost all of the error."""
    arms = [(_load(p), Path(p).stem) for p in paths]
    modes = {(name, r.get("assignment_mode")) for a, name in arms for r in a.values()}
    if len({m for _n, m in modes}) > 1:
        # `--set em.assignment_mode=fractional` moves every count, so a mix would report the mode's
        # effect as an arm's; a file written before the stamp existed is a different mode, not agreement.
        raise SystemExit(
            f"⛔ the arms were scored under different assignment modes {sorted(modes, key=str)} — "
            "compare arms run under one mode."
        )
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
    from rigel.config import CalibrationConfig as _CC

    _cc = _CC()
    print(f"  {len(conds)} conditions   arms: {', '.join(n for _a, n in arms)}"
          f"   calibration: the shipped defaults (message_policy={_cc.message_policy!r})")
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
    # Axiom 0 guard, and the reason this table needs a header. The solver has three
    # populations — gDNA, RNA+, RNA− — and "nascent" is not one of them and never becomes one. This is
    # an ASSIGNMENT question, not a deconvolution one: `nrna_est` is mass the EM parked on nascent
    # entities, and reading it as a fourth component is the exact error `CLAUDE.md`'s AXIOM 0 exists to
    # prevent.
    # On a `nrna_none` panel `nrna_true` is exactly 0, so every fragment in that column is a false
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
        print(f"    {'stratum':<26} {'pool':<20} {'est':>15} {'true':>15} {'Δ':>15} {'Δ%':>8}")
        print("    " + "-" * 102)

        def pool_rows(label, sel, a=a):
            def tot(field):
                return sum(a[(c, "library")][field] for c in conds
                           if sel(c) and (c, "library") in a)

            # Intergenic fragments are gDNA and the numerator must say so. `gdna_est` is
            # `gdna_em_count`, which EXCLUDES fragments that reached no locus; `gdna_true` is the
            # simulator's origin count, which includes them. Scoring one against the other is the
            # documented mistake in `score_library`'s own docstring — it read 0.3151 against a truth of
            # 0.5000 and fabricated an off-capture EM under-call that does not exist. Measured here
            # before the fix: -50.7 % panel-wide, which is the intergenic pool and nothing else.
            # The `of which` rows are informational and carry no truth of their own; only the TOTAL
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
                    # a relative error against a TRUE ZERO is not a large number, it is undefined —
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
    # The stratum roll-up above answers "which stratum is broken"; this answers "what happened in
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
    ap.add_argument("--markdown", type=Path, default=None,
                    help="with --report: also write the full per-scenario report as markdown. The "
                         "first --report file is the arm reported; an arm whose name contains "
                         "'reseed' becomes the attribution floor printed beside it")
    ap.add_argument("--arm", choices=ARMS, default=None)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--em-seed", type=int, default=DEFAULT_EM_SEED,
                    help="the shipped seed by default; it reaches only the sampled assignment's draw "
                         "— see DEFAULT_EM_SEED")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--set", dest="settings", action="append", default=[], metavar="SECTION.FIELD=VALUE",
                    help="a config value applied to every arm, repeatable (e.g. --set "
                         "em.assignment_mode=fractional); the same parser as every instrument's --set")
    args = ap.parse_args()

    if args.report:
        report(args.report)
        if args.markdown:
            markdown_report(args.report, args.markdown)
        return 0
    if args.markdown:
        raise SystemExit("--markdown needs --report FILES... (it renders arm jsonl, it runs nothing)")
    if args.arm is None or args.out is None:
        raise SystemExit("--arm and --out are required (or use --report)")

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    cache = args.oracle_cache
    if cache is None and (args.suite / "oracle_cache").is_dir():
        cache = args.suite / "oracle_cache"
    if args.arm in _RULER_ARMS:
        # the capture truth is minutes per label and shared by every condition: build it once, here,
        # before any shard could race to build it too
        truth_index = TranscriptIndex.load(str(args.index))
        for label in sorted({n.endswith("_capture_on") for n in names}):
            first = next(n for n in names if n.endswith("_capture_on") == label)
            capture_truth_factor(args.suite, truth_index, first)

    if args.jobs > 1 and len(names) > 1:
        # Shards, not threads. Conditions share nothing but a read-only index and cache, so this
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
            for spec in args.settings:
                cmd += ["--set", spec]
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
            # A short output file reads as a complete panel (TRAPS: an-ablation-that-never-ran's shape).
            raise SystemExit("a shard failed; refusing to concatenate a partial panel")
        with args.out.open("w") as fh:
            for o in outs:
                fh.write(o.read_text())
        print(f"  ⭐ {args.arm}: {sum(1 for _ in args.out.open())} rows -> {args.out}")
        return 0

    index = TranscriptIndex.load(str(args.index))
    pipeline_config = PipelineConfig()
    for spec in args.settings:
        pipeline_config = set_field(pipeline_config, spec)
        print(f"  ⭐ --set {spec}", flush=True)
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
