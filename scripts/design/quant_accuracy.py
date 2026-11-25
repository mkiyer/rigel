#!/usr/bin/env python
"""How accurate is the tool end to end, and what is a perfect prior worth?

``--arm base`` runs the shipped pipeline on a simulated condition and scores its transcript table against the
simulator's per-transcript truth. Every other arm runs the same pipeline with one thing substituted, and its
difference from ``base`` is what that one thing is worth:

- ``oracle``: the ``LocusPriors`` built from the origin-split truth (needs the oracle cache);
- ``oracle_ruler``: the simulator's capture-aware effective length in place of the shipped ruler's, the only arm
  that reaches the lengths the EM divides by, on every reference the condition drew gDNA from;
- ``oracle_alloc_seed``: the true per-transcript allocation weights, a capability proof, not a headroom claim;
- ``no_capture``: every capture efficiency read as 1, the line a capture reader must beat: under it a row is worth
  what reading no capture at all is worth. It substitutes no truth, so it scores an aligner's BAM too;
- ``noop`` and ``oracle_ruler_noop``: the full wrappers taking nothing; they must be byte-identical to ``base``.

One scorer serves every arm. The primary score is count against count: the truth is each condition's realised
fragment count, never the pre-capture molar abundance, and no length model enters it. The TPM rows use the tool's
own effective length on both sides, so they measure assignment only. The library gDNA fraction counts unspliced
intergenic fragments as gDNA, as ``cli.py`` does.

Reading rules: run an A/B pair with the scan pinned (``--set scan.total_threads=1``) so it is exactly reproducible,
judge an effect by its size, and read genes and pools beside the transcript table. The oracle masses are undrained
while the shipped pipeline drains, a small conservative bias.

Gate: ``tests/calibration/test_quant_accuracy.py``.

Usage::

    python scripts/design/quant_accuracy.py --arm base --out $RIGEL_ARMS/qa_base.jsonl --jobs 4
    python scripts/design/quant_accuracy.py --arm oracle --oracle-cache DIR --out $RIGEL_ARMS/qa_oracle.jsonl
    python scripts/design/quant_accuracy.py --arm base --conditions COND --suite DIR --index INDEX
    python scripts/design/quant_accuracy.py --report $RIGEL_ARMS/qa_base.jsonl $RIGEL_ARMS/qa_oracle.jsonl
    python scripts/design/quant_accuracy.py --arm base --set em.assignment_mode=fractional --out F   # a config arm
    python scripts/design/quant_accuracy.py --arm base --bam rustar_recommended.bam --out F   # an aligner's input

``--bam`` reads each condition from an aligner's BAM (``panel.py align``) instead of the simulator's alignments,
scored against the same truth so the two inputs' rows compare directly; the library row adds what the aligner left
out by origin (``unaligned_*``) and the gDNA fraction's truth over what reached Rigel (``gdna_frac_true_aligned``).
``--report`` prints the per-stratum tables of one or more arm files; ``--markdown`` also writes the per-scenario
report of the first.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _oracle import OracleTruth  # noqa: E402
from _shared import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_SUITE,
    DEFERRED,
    ORACLE_BAM,
    add_set_option,
    capture_config,
    capture_label,
    code_revision,
    config_overrides,
    configured,
    oracle_cache_dir,
    panel_conditions,
    run_shards,
    set_args,
    sibling,
    strata,
    stratum,
    stratum_selections,
)

import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import run_pipeline  # noqa: E402
from _scan_cache import ScanCacheKeyError, read_scan_cache  # noqa: E402

#: ``oracle_alloc_seed`` hands the EM the true relative abundances as allocation weights and keeps the shipped
#: coverage seed. It is no headroom claim: truth-derived weights also hand over the true support (a zero weight is
#: absorbing, so every silent transcript is switched off for free).
#:
#: The ruler arms substitute the lengths the EM divides by. The ruler (``transcript_capture_eff_lengths``) is built
#: by ``pipeline.transcript_geometry`` before ``assemble_priors`` runs, so no prior arm reaches it::
#:
#:     calibrate(...)
#:       |-- transcript_capture_eff_lengths()  <- the ruler arms substitute here: `effective_lengths_em`
#:       |-- assemble_priors()                 <- every prior arm wraps this: `LocusPriors`
#:
#: ``oracle_ruler`` hands the EM ``fl × factor``, the simulator's capture factor per transcript on the shipped scale,
#: where the uncaptured bulk reads 1 (``anchor_capture_factor``), never clipped; the priors and the locus gDNA length
#: stay as shipped. Probes can make capture isoform-specific and isoforms split by the ratio of their lengths, so this arm
#: prices the ruler end to end. Capture-OFF its factor is exactly 1, so it hands the EM RNA's own length, ``fl``; where
#: the shipped ruler returns another length with every efficiency at 1, the arm prices that difference too.
#:
#: The factor's scale is the kept gDNA length's, which reads capture off gDNA's own fragments. On a reference the
#: condition drew no gDNA from (an RNA-only spike-in, or every reference of a zero-gDNA condition) nothing measured
#: capture: the kept gDNA length there is the uncaptured one, and a probed transcript at ``fl × factor`` would stand its
#: whole capture factor above it, so the E-step would hand the gDNA component every fragment the two explain alike.
#: There both arms keep the shipped length (``gdna_free_transcripts``), and the arm prices nothing.
_RULER_ARMS = {"oracle_ruler": True, "oracle_ruler_noop": False}

ARMS = ("base", "noop", "oracle", "oracle_alloc_seed", "no_capture") + tuple(_RULER_ARMS)

#: The arms that substitute no truth, so they score any BAM of the condition's reads (``--bam``): every other arm
#: injects truth built from the simulator's own alignments.
_TRUTH_FREE_ARMS = ("base", "no_capture")

#: arm -> which ``LocusPriors`` fields come from the oracle O. ``noop`` takes none and still builds O, so it tests
#: the wrapper rather than an ``if``. There is no RNA arm: calibration's RNA masses reach the EM only inside
#: ``fragment_count``, each object's whole tally, which the oracle's split of it leaves where it was, so it could not
#: fire. ``oracle`` takes O whole, its fragment count with its gDNA count, so the share is one result's.
_ARM_FIELDS = {
    "base": (),
    "noop": (),
    "oracle": ("gdna_count", "fragment_count", "gdna_eff_len"),
}


# ── the injection ────────────────────────────────────────────────────────────────────────────────


def load_oracle(index, pipeline_config, suite: Path, cache_root: Path, tag: str) -> OracleTruth:
    """The origin-split truth for one condition, read from the scan caches.

    Sum-to-full is asserted against the condition's undrained scan cache, so a partition built from another scan
    is refused. ``read_scan_cache`` refuses a payload whose ``graph_hash``, ``reach_digest``,
    ``payload_schema_digest`` or scan config does not match this index, so a stale cache fails loudly.
    """
    scan = pipeline_config.scan
    try:
        full = read_scan_cache(Path(suite) / "scan_cache" / tag, index, scan).payload
        parts = sibling("_oracle_arms.py").oracle_parts(cache_root, tag, index, scan=scan)
    except (FileNotFoundError, KeyError, ScanCacheKeyError) as exc:
        raise SystemExit(
            f"⛔ {tag}: no valid scan or oracle cache ({exc}). Build both first with "
            "panel.py cache; this script refuses to invent a truth."
        ) from exc
    return OracleTruth.from_parts(full, parts)


def install_arm(arm: str, oracle: OracleTruth | None):
    """Wrap :func:`~rigel.calibration.priors.assemble_priors` for one arm. Returns ``(restore, fired)``.

    The wrapper is the same for every arm, ``noop`` included: it builds O on the loci the run produced, and only
    then does the arm pick which fields to take, so ``noop`` proves the whole path inert. O is built on the run's
    own ``multi_loci``, never a stored array: ``build_multi_loci`` unions transcripts linked by scored fragments,
    so the partition is a function of the run and a locus id from another process means nothing here.

    ``fired`` is a counter and the caller raises on zero, because an override that never ran reads as "no effect".
    """
    original = PRIORS.assemble_priors

    if arm == "base":
        # base installs nothing, so it counts as fired: the never-ran check must not fail on the one arm
        # with no override
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


def _count_scores(t: np.ndarray, e: np.ndarray) -> dict:
    """The count-against-count scores of a transcript or gene row: truth ``t`` and estimate ``e``, aligned."""
    d = e - t
    expressed = t > 0.0
    silent = ~expressed
    sp, pe = _spearman_pearson(t[expressed], e[expressed])
    # symmetric MARD over the expressed set: bounded in [0, 1], finite at a zero estimate
    denom = np.abs(e[expressed]) + np.abs(t[expressed])
    mard = float(np.mean(np.abs(d[expressed]) / np.where(denom > 0, denom, 1.0))) if expressed.any() else float("nan")
    return {
        "n_tx": int(t.size),
        "n_expressed": int(expressed.sum()),
        "n_detected": int((e > 0).sum()),
        "count_true": float(t.sum()),
        "count_est": float(e.sum()),
        "count_abs_err": float(np.abs(d).sum()),
        "count_net_err": float(d.sum()),
        "count_over": float(np.maximum(d, 0.0).sum()),
        "count_under": float(np.maximum(-d, 0.0).sum()),
        #: contamination the deconvolution left on a silent transcript
        "fp_mass": float(e[silent].sum()),
        "fp_n": int((e[silent] > 0).sum()),
        #: real RNA the tool did not assign to the transcript that produced it
        "fn_mass": float(t[expressed & (e == 0.0)].sum()),
        "fn_n": int((expressed & (e == 0.0)).sum()),
        "spearman": sp,
        "pearson_log2": pe,
        "mard": mard,
    }


def score_transcripts(quant: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """Score one arm's transcript table against the condition's realised fragment truth.

    ``count`` against ``observed_annotated_fragments``, both fragment counts with no length or normalisation between
    them, so a difference is a fragment the tool put somewhere the simulator did not. ``fp_mass`` is the estimate
    on transcripts the simulator gave zero fragments: contamination the deconvolution failed to remove, with
    nothing to cancel it. ``fn_mass`` is its mirror, real RNA the tool assigned nowhere. They are reported apart
    because one ``Σ|Δ|`` cannot separate two failures with different fixes. Every transcript of either table is
    scored (:func:`_joined`): one the truth does not hold is an estimate with truth 0, never dropped.
    """
    m = _joined(quant, truth, ["count", "count_em", "tpm", "effective_length"])
    t = m["observed_annotated_fragments"].to_numpy(np.float64)
    e = m["count"].to_numpy(np.float64)
    expressed = t > 0.0

    # truth TPM uses the tool's own effective length, so the length model cancels and the row measures assignment.
    # A transcript with no effective length has no rate: a shadow transcript the tool does not hold, or one it holds
    # at length 0. It carries no TPM on either side, so it is out of both vectors and both denominators; its
    # fragments are reported apart, never hidden, and the length-free count rows still score them.
    eff = m["effective_length"].to_numpy(np.float64)
    has_len = eff > 0.0
    rate = np.where(has_len, t / np.where(has_len, eff, 1.0), 0.0)
    tpm_true = _per_million(rate)
    tpm_est = _per_million(np.where(has_len, m["tpm"].to_numpy(np.float64), 0.0))
    rel = np.abs((e - t)[expressed]) / t[expressed] if expressed.any() else np.array([])
    return {
        **_count_scores(t, e),
        "median_rel_err": float(np.median(rel)) if rel.size else float("nan"),
        "tpm_abs_err": float(np.abs(tpm_est - tpm_true).sum()),
        "tpm_fp": float(tpm_est[~expressed].sum()),
        "n_no_length": int((~has_len).sum()),
        "no_length_fragments": float(t[~has_len].sum()),
    }


def _joined(quant: pd.DataFrame, truth: pd.DataFrame, estimate_cols: list[str]) -> pd.DataFrame:
    """Truth and estimate over the union of their transcripts, one row each: a transcript only one side holds reads 0
    on the other. Each side must name a transcript once (a repeated id would count its fragments twice), and every
    estimated fragment must survive the join: the scored total is reconciled against the table's."""
    m = truth.merge(
        quant[["transcript_id", *estimate_cols]], on="transcript_id", how="outer", validate="one_to_one"
    )
    m["observed_annotated_fragments"] = m["observed_annotated_fragments"].fillna(0.0)
    for col in estimate_cols:
        m[col] = m[col].fillna(0.0)
    scored, table = float(m["count"].sum()), float(quant["count"].sum())
    if not np.isclose(scored, table, rtol=1e-12, atol=0.0):
        raise ValueError(f"the scored estimate ({scored!r}) is not the table's ({table!r})")
    return m


def _per_million(x: np.ndarray) -> np.ndarray:
    """``x`` scaled to sum to one million, or zeros when it sums to nothing."""
    total = x.sum()
    return x / total * 1e6 if total > 0 else np.zeros_like(x)


def score_genes(quant: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """The same scorer over genes; its difference from the transcript row separates isoform ambiguity from the rest.

    Summing a gene's isoforms removes the error of not knowing which isoform a fragment came from. What remains is
    error in deciding whether a fragment was RNA from this gene at all, the question Rigel is for. Genes come from
    the truth table's ``gene_id``, so both sides share one grouping; a transcript the truth does not hold is grouped
    by the estimate's own ``gene_id``, which it must then carry. Synthetic nascent entities are absent from ``quant``
    (``get_counts_df`` drops them); their mass is on the library row as ``nascent_est``.
    """
    extra = ~quant["transcript_id"].isin(truth["transcript_id"])
    m = _joined(quant, truth, ["count"])
    if extra.any():
        if "gene_id" not in quant.columns:
            raise ValueError(f"{int(extra.sum())} estimated transcripts the truth does not hold carry no gene_id")
        estimate_gene = quant.loc[extra].set_index("transcript_id")["gene_id"].astype(str)
        gene = m["gene_id"].astype(object)
        gene = gene.where(gene.notna(), m["transcript_id"].map(estimate_gene))
        m["gene_id"] = gene
    g = m.groupby("gene_id", observed=True)[["observed_annotated_fragments", "count"]].sum()
    # `n_tx` counts genes here, named as on the transcript row so one table renders both
    return _count_scores(g["observed_annotated_fragments"].to_numpy(np.float64), g["count"].to_numpy(np.float64))


def score_library(result, quant: pd.DataFrame, truth_summary: dict) -> dict:
    """The library-level split by pool, and the gDNA fraction a user reads.

    ``n_intergenic`` counts fragments that reached no locus, all unspliced (the scan withholds a spliced hit no
    transcript holds); they never enter the EM, so ``estimator.gdna_em_count`` excludes them. Off capture, gDNA is
    genome-uniform and much of it is intergenic, so ``gdna_frac_est`` is ``cli.py``'s ``fraction_gdna``,
    ``(gdna_em + n_intergenic) / (rna + gdna_em + n_intergenic)``; ``gdna_em_frac_est`` is the EM's own share.
    The truth is the simulator's ``origin_counts``, which counts each fragment once.
    """
    est = result.estimator
    annotated = float(quant["count"].sum())
    nascent = float(est.synthetic_count)
    gdna = float(est.gdna_em_count)
    intergenic = float(result.stats.n_intergenic) + float(est.gdna_intergenic_em_count)
    oc = truth_summary.get("origin_counts", {})
    t_g = float(oc.get("gdna", float("nan")))
    t_m = float(oc.get("annotated", 0.0))
    t_n = float(oc.get("nascent", 0.0))
    t_r = t_m + t_n
    total = annotated + nascent + gdna
    total_all = total + intergenic
    return {
        # The nascent pool is broken out: ``get_counts_df`` drops the synthetic entities, so mass the EM
        # assigns them is invisible in the transcript table's ``fp_mass``. On an ``nascent_none`` panel its
        # truth is 0, which makes ``nascent_est`` a third false-positive channel.
        "annotated_est": annotated,
        "annotated_true": t_m,
        "nascent_est": nascent,
        "nascent_true": t_n,
        "rna_est": annotated + nascent,
        "rna_true": t_r,
        "gdna_est": gdna,
        "gdna_true": t_g,
        #: the reported fraction: unspliced intergenic fragments included, as cli.py computes it
        "gdna_frac_est": (gdna + intergenic) / total_all if total_all > 0 else float("nan"),
        #: the EM's own share, intergenic on neither side
        "gdna_em_frac_est": gdna / total if total > 0 else float("nan"),
        "gdna_frac_true": t_g / (t_g + t_r) if (t_g + t_r) > 0 else float("nan"),
        "n_intergenic_gdna": intergenic,
    }


def gdna_estimate(library: dict) -> float:
    """A library row's gDNA estimate, the one ``gdna_true`` is scored against: the EM's gDNA plus the
    unspliced intergenic fragments, which reach no locus and so never enter the EM (``score_library``)."""
    return library["gdna_est"] + library["n_intergenic_gdna"]


# ── an aligned input ─────────────────────────────────────────────────────────────────────────────


def check_input_bam(suite: Path, names: list[str], arm: str, bam_name: str) -> None:
    """Refuse an input no row could honestly come from. Only an arm that substitutes no truth (``base``,
    ``no_capture``) reads another BAM than the simulator's: every other arm substitutes truth built from the
    simulator's alignments. Every condition must hold the named BAM, or the panel would be scored short and read as
    complete."""
    if bam_name == ORACLE_BAM:
        return
    if arm not in _TRUTH_FREE_ARMS:
        raise SystemExit(
            f"⛔ --bam {bam_name} reads an aligner's output, which only {' or '.join(_TRUTH_FREE_ARMS)} can score: "
            f"arm {arm!r} substitutes truth built from the simulator's own alignments"
        )
    missing = [n for n in names if not (suite / n / bam_name).is_file()]
    if missing:
        raise SystemExit(f"⛔ no {bam_name} in {len(missing)} condition(s): {', '.join(missing)} (panel.py align)")


def aligner_loss(bam: Path, origin_counts: dict) -> dict:
    """What an aligner left out, by true origin: the templates of its BAM with no mapped record.

    The simulator's truth counts every fragment it drew, including gDNA drawn over N runs no aligner can place,
    so these counts say how much of a pool's gap is the aligner's. A template is lost only when none of its
    records is mapped; one placed mate reached Rigel. The BAM must hold each simulator template once, in one run
    of records (the read name gives its origin); anything else is refused."""
    import pysam

    from rigel.sim.read_name import parse_origin

    lost = {"gdna": 0, "annotated": 0, "nascent": 0}
    seen = dict(lost)
    name, mapped = None, False

    def close():
        kind = parse_origin(name).kind
        seen[kind] += 1
        lost[kind] += not mapped

    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.query_name != name:
                if name is not None:
                    close()
                name, mapped = r.query_name, False
            mapped = mapped or not r.is_unmapped
    if name is not None:
        close()
    want = {k: int(origin_counts.get(k, 0)) for k in seen}
    if seen != want:
        raise SystemExit(
            f"⛔ {bam}: its templates by origin {seen} are not the simulator's {want}: a template is missing, "
            "or its records are not in one run (the input is not collated)"
        )
    return {f"unaligned_{k}": v for k, v in lost.items()}


def aligned_library_truth(library: dict, bam: Path, summary: dict) -> dict:
    """The library row's additions for an aligner's input: what the aligner left out by true origin, and the gDNA
    fraction's truth over the templates that reached Rigel. The truth itself stays the simulator's."""
    loss = aligner_loss(bam, summary.get("origin_counts", {}))
    t_g = library["gdna_true"] - loss["unaligned_gdna"]
    t_r = library["rna_true"] - loss["unaligned_annotated"] - loss["unaligned_nascent"]
    return {**loss, "gdna_frac_true_aligned": t_g / (t_g + t_r) if (t_g + t_r) > 0 else float("nan")}


# ── one condition ────────────────────────────────────────────────────────────────────────────────


#: The truth table's realised RNA fragment columns, mature and nascent. A synthetic span emits no mature RNA, so
#: its whole realised mass sits in the nascent column of its own row; the mature column alone weights it at zero.
_WEIGHT_COLUMNS = ("observed_annotated_fragments", "observed_nascent_fragments")


def truth_weights(truth: pd.DataFrame, index) -> np.ndarray:
    """``float64[n_transcripts]``: the true realised RNA fragment count per transcript, on the EM's axis.

    This is ``oracle_alloc_seed``'s whole input. The EM reads only ratios within a locus, so no normalisation is
    needed. The weight is mature plus nascent, because every RNA component receives the RNA prior and a synthetic
    span is one of them; a mature-only weight would switch every synthetic span off. A table with no nascent
    column is refused rather than weighted at zero.

    Exact-duplicate transcripts fold onto the twin the index kept: the truth table is keyed on the un-collapsed
    annotation, and for such a pair only the group total is defined.
    """
    missing = [c for c in _WEIGHT_COLUMNS if c not in truth.columns]
    if missing:
        raise SystemExit(
            f"⛔ truth table lacks {'+'.join(missing)} of the realised mature+nascent RNA pair "
            f"({'+'.join(_WEIGHT_COLUMNS)}); a mature column alone cannot express a synthetic "
            "nascent entity's allocation"
        )
    mature, nascent = _WEIGHT_COLUMNS
    total = truth[mature].to_numpy(np.float64) + truth[nascent].to_numpy(np.float64)
    t_index = dict(zip(index.t_df["t_id"].to_numpy(), index.t_df["t_index"].to_numpy(), strict=True))
    w = np.zeros(int(index.num_transcripts), dtype=np.float64)
    for tid, n in zip(truth["transcript_id"], total, strict=True):
        i = t_index.get(str(tid))
        if i is not None:
            w[int(i)] += float(n)
    return w


def anchor_capture_factor(factor, probed_frac, plain_length) -> np.ndarray:
    """The simulator's capture factor on the shipped scale, where the uncaptured bulk reads 1.

    The shipped efficiencies are yields against the bulk, so the EM's lengths put an unprobed transcript at its plain
    length and a probed one above it, on the scale of the locus gDNA length the ruler arm leaves as shipped. The
    simulator's factor carries the off-target weight in both its counts, so it is already on that scale; it is divided
    by the median over the unprobed class (probed fraction 0), which reads 1 where the simulator's convention holds,
    and kept as it is when every transcript is probed. The identity where the factor is undefined (no plain length).
    ``ruler_vs_truth.score`` compares shapes, blind to a global scale, and anchors its errors on its own."""
    f = np.asarray(factor, dtype=np.float64)
    ok = np.isfinite(f) & (f > 0.0) & (np.asarray(plain_length, dtype=np.float64) > 0.0)
    logf = np.where(ok, np.log(np.where(ok, f, 1.0)), 0.0)
    bulk = ok & (np.asarray(probed_frac, dtype=np.float64) == 0.0)
    scale = float(np.median(logf[bulk])) if bulk.any() else 0.0
    return np.where(ok, np.exp(logf - scale), 1.0)


def capture_truth_factor(suite: Path, index, condition: str) -> np.ndarray:
    """``float64[n_transcripts]``: the anchored capture factor on the index's axis for the condition's capture label.

    It depends on the capture config, the probe file, the simulated fragment-length law and the annotation, never
    on the condition, so it is computed once per label and cached at
    ``<suite>/oracle_cache/capture_truth_<label>.npz`` under a key over all four; a mismatched key is recomputed.
    Computing it integrates the sampler over every fragment length and takes minutes."""
    import hashlib

    label = capture_label(condition)
    manifest = json.loads((Path(suite) / "manifest.json").read_text())
    capture = capture_config(manifest, condition)
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
    truth = sibling("ruler_vs_truth.py").load_truth(index, Path(index.index_dir), Path(suite), condition)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.stem}.{os.getpid()}.npz")
    np.savez(tmp, key=key, factor=truth.factor, probed_frac=truth.probed_frac, L_plain=truth.L_plain)
    tmp.replace(path)
    return anchor_capture_factor(truth.factor, truth.probed_frac, truth.L_plain)


def gdna_free_transcripts(manifest: dict, index, condition: str) -> np.ndarray:
    """``bool[n_transcripts]``: the transcripts on a reference *condition* drew no gDNA from, the RNA-only spike-ins and
    every reference of a zero-gDNA condition, read off the manifest (``calibration_oracle.genomic_ref_ids``), never off
    where gDNA landed. A condition the manifest does not list is refused."""
    refs = sibling("calibration_oracle.py").genomic_ref_ids(manifest, condition, index.ref_name_to_id)
    t_ref = index.t_df["ref"].map(index.ref_name_to_id).to_numpy(np.int64)
    return ~np.isin(t_ref, refs)


def install_ruler_arm(arm: str, factor: np.ndarray, held: np.ndarray):
    """Wrap the ruler so the EM divides by ``fl × factor``, and by the shipped length where ``held``. Returns
    ``(restore, fired)``.

    ``held`` marks the transcripts on a reference the condition drew no gDNA from (``gdna_free_transcripts``), where
    the factor has no gDNA length on its own scale to stand beside. ``transcript_capture_eff_lengths`` is patched as a
    module attribute, which works because ``pipeline.transcript_geometry`` imports it at call time. Both arms
    run the same wrapper: it runs the shipped ruler, builds the substitute, and records how far the substitute sits from
    the shipped lengths (``truth_delta``), how far the returned lengths do (``max_abs_delta``) and how many transcripts
    it held (``held``); ``oracle_ruler_noop`` returns the shipped object. ``check_ruler_arm`` reads the counter.
    """
    from rigel.calibration import capture_eff_length as CEL

    substitute = _RULER_ARMS[arm]
    factor = np.asarray(factor, dtype=np.float64)
    held = np.asarray(held, dtype=bool)
    orig_ruler = CEL.transcript_capture_eff_lengths
    fired = {"n": 0, "ruler": 0, "max_abs_delta": 0.0, "truth_delta": 0.0, "held": int(held.sum())}

    def ruler_wrapper(calibration, frame, fl_eff_lengths, junction_cap):
        shipped = orig_ruler(calibration, frame, fl_eff_lengths, junction_cap)
        base = np.asarray(shipped, dtype=np.float64)
        truth = np.where(held, base, np.asarray(fl_eff_lengths, dtype=np.float64) * factor)
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
    """Refuse a ruler arm that could not have measured what it names: the ruler must have run, the substituting
    arm must move the EM's lengths wherever the truth differs from the shipped ones (where it does not, the identity
    is the answer), and the noop must move nothing."""
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


def install_no_capture():
    """Read every capture efficiency as 1. Returns ``(restore, fired)``.

    ``calibrate._capture_fields`` is the one place the result takes its capture reading, so the wrapper replaces the
    reader's answer and nothing else: every region's and boundary's efficiency is 1 and the saturated level is 1, as a
    library with one gDNA level reads, at a capture weight of 0, and the bulk density is calibration's own. Everything downstream (the
    transcripts' lengths, the locus gDNA lengths, a multimapper's place weights, RNA's captured fraction) then reads no
    capture. The arm prices the capture reader as a whole: a row where ``base`` is worse than this arm is a row where
    reading capture cost more than it gave. ``fired`` counts the readings replaced; ``check_no_capture`` reads the
    result.
    """
    # the package re-exports the function `calibrate` under the module's name, so the module is imported by name
    CALIBRATE = importlib.import_module("rigel.calibration.calibrate")
    inner = CALIBRATE._capture_fields
    fired = {"n": 0}

    def wrapper(*args, **kwargs):
        fields = inner(*args, **kwargs)
        fired["n"] += 1
        return dict(
            fields,
            gdna_capture_weight=0.0,
            gdna_capture_saturation=1.0,
            gdna_capture_efficiency_region=np.ones_like(fields["gdna_capture_efficiency_region"]),
            gdna_capture_efficiency_boundary=np.ones_like(fields["gdna_capture_efficiency_boundary"]),
        )

    CALIBRATE._capture_fields = wrapper

    def restore():
        CALIBRATE._capture_fields = inner

    return restore, fired


def check_no_capture(condition: str, calibration) -> None:
    """The ``no_capture`` arm's claim is that the run read no capture; refuse a result that read any."""
    read = [
        np.asarray(calibration.gdna_capture_efficiency_region),
        np.asarray(calibration.gdna_capture_efficiency_boundary),
    ]
    if any(np.any(c != 1.0) for c in read) or calibration.gdna_capture_saturation != 1.0:
        raise RuntimeError(
            f"{condition} [no_capture]: the run read a capture efficiency other than 1; the arm never reached the "
            "result. This is not a measurement of reading no capture."
        )


def install_truth_weights(weights: np.ndarray):
    """Hand the EM per-transcript weights as the RNA prior's shares. Returns ``(restore, fired)``.

    A wrapper rather than a config field: the weights are a per-run array, and an experiment should not add
    production surface. The pipeline passes no ``rna_prior_weight``, so the wrapper sets it on the estimator's public
    batch call. The counter is not the injection: it counts a call only when the weights reach a transcript that call
    solves, so a run whose EM never ran, or a truth table that weighs none of the EM's transcripts, cannot read as
    healthy.
    """
    from rigel.estimator import AbundanceEstimator

    inner = AbundanceEstimator.run_batch_locus_em_partitioned
    fired = {"n": 0}

    def wrapper(self, partitions, locus_transcript_indices, *args, **kw):
        solved = np.concatenate([np.asarray(t, dtype=np.int64) for t in locus_transcript_indices])
        if float(weights[solved].sum()) > 0.0:
            fired["n"] += 1
        return inner(self, partitions, locus_transcript_indices, *args, rna_prior_weight=weights, **kw)

    AbundanceEstimator.run_batch_locus_em_partitioned = wrapper

    def restore():
        AbundanceEstimator.run_batch_locus_em_partitioned = inner

    return restore, fired


def install_condition_arm(arm: str, suite: Path, index, condition: str, pipeline_config, oracle_cache: Path | None,
                          truth: pd.DataFrame):
    """Install ``arm``'s one substitution for ``condition``; returns ``(restore, fired)``."""
    if arm == "oracle_alloc_seed":
        # the shipped coverage seed is kept, so the arm varies one thing against `base`: the allocation
        return install_truth_weights(truth_weights(truth, index))
    if arm in _RULER_ARMS:
        # checked before the fall-through: `install_arm` has no fields for a ruler arm
        manifest = json.loads((Path(suite) / "manifest.json").read_text())
        return install_ruler_arm(arm, capture_truth_factor(suite, index, condition),
                                 gdna_free_transcripts(manifest, index, condition))
    if arm == "no_capture":
        return install_no_capture()
    oracle = None
    if arm != "base":
        if oracle_cache is None:
            raise SystemExit(f"⛔ arm {arm!r} needs --oracle-cache")
        oracle = load_oracle(index, pipeline_config, suite, oracle_cache, condition)
    return install_arm(arm, oracle)


def row_stamp(arm: str, condition: str, pipeline_config, bam_name: str) -> dict:
    """What every row of one arm on one condition carries, so the config a number was scored under travels with it
    and the report reads it from the row, never from the shipped defaults at render time: the ``--set`` overrides
    (every field that differs from the shipped config), the message policy whatever set it, the assignment mode,
    the EM seed, the input BAM and the code revision."""
    return {
        "arm": arm,
        "condition": condition,
        "settings": config_overrides(pipeline_config),
        "message_policy": pipeline_config.calibration.message_policy,
        "em_seed": int(pipeline_config.em.seed),
        # so a fractional run cannot be reported as a sampled one
        "assignment_mode": pipeline_config.em.assignment_mode,
        # and the input it read, so an aligner's rows cannot be read as the simulator's
        "bam": bam_name,
        # provenance, never a refusal: an A/B of a `src/` change is two revisions by design
        "revision": code_revision(),
    }


def run_condition(arm: str, suite: Path, index, condition: str, pipeline_config,
                  oracle_cache: Path | None,
                  bam_name: str = ORACLE_BAM) -> list[dict]:
    """One arm on one condition: its transcript, gene and library rows."""
    check_input_bam(suite, [condition], arm, bam_name)
    bam = str(suite / condition / bam_name)
    truth = pd.read_csv(suite / condition / "truth_abundances.tsv", sep="\t")
    summary = json.loads((suite / condition / "truth_summary.json").read_text())

    restore, fired = install_condition_arm(arm, suite, index, condition, pipeline_config, oracle_cache, truth)
    start = time.perf_counter()
    try:
        result = run_pipeline(bam, index, pipeline_config)
    finally:
        restore()
    if fired["n"] == 0:
        # an injection that never ran would read as "a perfect prior changes nothing"
        raise RuntimeError(
            f"{condition} [{arm}]: the override was never wrapped-and-called. This is not a "
            "measurement of zero effect."
        )
    if arm in _RULER_ARMS:
        # the ruler arms carry a second requirement: the EM's lengths are what they name (check_ruler_arm)
        check_ruler_arm(condition, arm, fired)
    if arm == "no_capture":
        check_no_capture(condition, result.calibration)
    seconds = time.perf_counter() - start

    quant = result.estimator.get_counts_df(index)
    common = {**row_stamp(arm, condition, pipeline_config, bam_name), "seconds": seconds}
    if arm in _RULER_ARMS:
        # how far the arm moved the EM's lengths and how far the truth sits from the shipped ruler, in bp on the
        # worst transcript, and how many transcripts it kept at the shipped length (no gDNA on their reference), so
        # the arm's reach is readable from the output
        common["ruler_max_abs_delta"] = float(fired["max_abs_delta"])
        common["ruler_truth_delta"] = float(fired["truth_delta"])
        common["ruler_n_held"] = int(fired["held"])
    library = score_library(result, quant, summary)
    if bam_name != ORACLE_BAM:
        library.update(aligned_library_truth(library, Path(bam), summary))
    return [
        {**common, "axis": "transcript", **score_transcripts(quant, truth)},
        # the same scorer over genes, isoform ambiguity summed away
        {**common, "axis": "gene", **score_genes(quant, truth)},
        {**common, "axis": "library", **library},
    ]


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def _load(path: Path) -> dict:
    """An arm's rows by ``(condition, axis)``, each once: a repeated row (a shard appended twice) would otherwise win by
    its position in the file."""
    rows = [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]
    by_key: dict = {}
    for r in rows:
        key = (r["condition"], r["axis"])
        if key in by_key:
            raise SystemExit(f"⛔ {path} holds the row {key} twice — rebuild the arm rather than pick one.")
        by_key[key] = r
    return by_key


def _gdna_level(library_row: dict) -> str:
    """The condition's gDNA level, read off its own truth (``gdna_frac_true``, the simulator's origin
    counts) rather than mapped from the name's rung token, so no rung a panel carries can be missing."""
    f = float(library_row["gdna_frac_true"])
    return f"{100.0 * f:.3g} %" if np.isfinite(f) else "n/a"


def _pct(est: float, true: float) -> str:
    """Percent error against the realised truth; ``n/a`` at a truth of zero (the zero-gDNA control), where a
    ratio would invent a scale."""
    if true <= 0.0:
        return "n/a"
    return f"{100.0 * (est - true) / true:+,.1f} %"


def _signed(x: float) -> str:
    return f"{x:+,.0f}"


#: the "How to read this" section of the markdown report, one bullet a line
_MD_HOW_TO_READ = (
    "- **Three pools, and they answer different questions.** `gDNA` is contamination the "
    "deconvolution must remove; `nascent RNA` is the SYNTHETIC entity this index manufactures for "
    "each multi-exon gene; `annotated RNA` is the transcript table a user reads. The split is on "
    "`is_synthetic`, never `is_single_exon` — a single-exon annotated transcript is ANNOTATED RNA.",
    "- **The gDNA estimate includes intergenic fragments**, which reach no locus and never enter "
    "the EM. Off capture they are more than half of all gDNA.",
    "- **Percent error is undefined where the truth is zero** and is printed `n/a`. At `g00` the "
    "true gDNA count is exactly 0; the raw count beside it is the whole of the answer there.",
    "- **Read per stratum, never pooled.** Three strata are in scope. `unstranded × capture ON` is "
    "DEFERRED — reported on every benchmark, never a development target — and it carries most of "
    "the error, so a pooled total would be its total.",
    "- **An effect is judged by its size, with genes and pools read beside the transcript table.** "
    "An A/B pair runs with the scan pinned (`--set scan.total_threads=1`) and is exactly "
    "reproducible.",
    "- **`expressed` and `detected` are the scored sets, and the truth table is larger than "
    "either.** It carries one row per SYNTHETIC nascent entity as well, and those rows are zero on "
    "both sides — zero truth and zero estimate, since the transcript table drops them — so they "
    "enter no figure here. A raw row count would read as thousands of transcripts scored perfectly; "
    "it is not reported for that reason.",
    "- **`false-negative mass` is the truth on the expressed transcripts the estimate gives EXACTLY "
    "zero.** A component the EM left at its floor is assigned nothing, so the tool reports the "
    "transcript as absent, under fractional assignment too. It is the part of `under-assigned` "
    "the tool gave up whole, as `false-positive mass` is the part of `over-assigned` on a silent "
    "transcript; `detected` counts the estimates above zero.",
    "- **`MARD` is the SYMMETRIC mean absolute relative difference**, `|est−true| / (|est|+|true|)` "
    "averaged over the expressed set, so it is bounded in [0, 1] and finite at a zero estimate — "
    "0.5 is a 3× error, not a 50 % one. `median rel. err` is the ordinary `|est−true|/true`, over "
    "the same set. Both are unweighted by mass: a 1-fragment transcript counts as much as a "
    "25,000-fragment one.",
    "- **One asymmetry in the pool rows.** The ESTIMATE splits on `is_synthetic` (the entity) while "
    "the TRUTH splits on the simulator's template kind, and a handful of annotated single-exon "
    "transcripts serve as their own nascent entity — their nascent fragments count in `nascent "
    "truth` while the tool's counts for them land in `annotated`. It is a small share of the pool, "
    "and it cannot be removed by relabelling; it is recorded so the pool rows are not read as one "
    "partition measured twice.",
    "- **Every stratum row includes its `g00` rung.** The terminal report excludes the zero-gDNA "
    "control from its roll-up; this one keeps it, so the two do not agree by construction.",
)

#: the pool tables: (estimate field, or None for the intergenic-inclusive gDNA estimate; truth field; title)
_MD_POOLS = (
    (None, "gdna_true", "gDNA (EM + intergenic)"),
    ("nascent_est", "nascent_true", "Nascent RNA — SYNTHETIC entities only"),
    ("annotated_est", "annotated_true", "Annotated RNA — the transcript table"),
)

#: the transcript and gene sections: (axis, heading, blurb)
_MD_AXES = (
    ("transcript", "2. Transcript level — inside the annotated RNA pool",
     "Scored over the annotated transcripts alone: gDNA and the synthetic nascent entities are "
     "excluded, so what is left is how well the tool splits the RNA it kept. `Σ|Δ|` sums "
     "`|estimate − truth|` over every annotated transcript, SILENT ONES INCLUDED, and splits into "
     "`over-assigned` + `under-assigned`; `false-positive mass` is the part of the over-assignment "
     "that landed on a transcript the simulator gave ZERO fragments, so it is a subset of `Σ|Δ|` "
     "and the number this tool exists for. `false-negative mass` is its mirror, the truth on an "
     "expressed transcript the tool gave exactly zero, and a subset of `under-assigned`. Each is "
     "followed by the number of transcripts it sits on."),
    ("gene", "3. Gene level — the same, with isoform ambiguity summed away",
     "The same scorer over genes, and the difference from the transcript row is the point: "
     "summing a gene's isoforms collapses exactly the error that comes from not knowing WHICH "
     "isoform a fragment came from. What survives is error in deciding whether the fragment was "
     "RNA from this gene at all, which is the question Rigel is for."),
)


def _by_stratum(conds: list[str]):
    """``(stratum, its conditions)`` for every stratum that has any, in report order, the zero-gDNA control
    included."""
    for st in strata(conds):
        rows = [c for c in conds if stratum(c) == st]
        if rows:
            yield st, rows


def _stratum_name(st, gap: str = "  ") -> str:
    """``strand × capture``, marked when it is the deferred stratum."""
    return f"{st[0]} × {st[1]}" + (f"{gap}⛔ DEFERRED" if st == DEFERRED else "")


def _md_header(primary: dict, pname: str, conds: list[str], path: Path) -> list[str]:
    """The title, what was scored and when, and how to read the report."""
    modes = {r.get("assignment_mode") for r in primary.values()}
    policy, settings, revision = _stamps(primary.values())
    stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime(Path(path).stat().st_mtime))
    return [
        "# Rigel — end-to-end accuracy on the gDNA ladder",
        "",
        f"**Arm** `{pname}` · **{len(conds)} conditions** · **assignment** "
        f"`{'/'.join(sorted(str(m) for m in modes))}` · **calibration** `{policy}` · **--set** `{settings}` · "
        f"**revision** `{revision}` · **scored** {stamp}",
        "",
        "Every number is scored against the simulator's own REALISED per-fragment truth — each fragment "
        "counted exactly once, never a pre-capture molar abundance — so a difference is a fragment the "
        "tool put somewhere the simulator did not.",
        "",
        "## How to read this",
        "",
        *_MD_HOW_TO_READ,
        "",
    ]


def _md_pools(primary: dict, conds: list[str]) -> list[str]:
    """Section 1: the three pools per scenario, then the library gDNA fraction."""
    lines = ["## 1. Pool level — where the library's fragments went", "",
             "One row per scenario, grouped by stratum. `Δ` is `estimate − truth` in fragments.", ""]
    for est_f, true_f, label in _MD_POOLS:
        lines += [f"### {label}", "", "| scenario | gDNA level | estimated | truth | Δ | % error |",
                  "|---|---|---:|---:|---:|---:|"]
        for st, rows in _by_stratum(conds):
            lines.append(f"| **{_stratum_name(st)}** | | | | | |")
            for c in rows:
                r = primary[(c, "library")]
                est = gdna_estimate(r) if est_f is None else r[est_f]
                true = r[true_f]
                lines.append(f"| `{c}` | {_gdna_level(r)} | {est:,.0f} | {true:,.0f} | "
                             f"{_signed(est - true)} | {_pct(est, true)} |")
        lines.append("")
    lines += ["### Library gDNA fraction — the thermometer", "",
              "What `rigel quant` reports as the library's gDNA share. This is the number calibration exists "
              "to produce; the transcript table does not always keep it.", "",
              "| scenario | estimated | truth | Δ |", "|---|---:|---:|---:|"]
    for st, rows in _by_stratum(conds):
        lines.append(f"| **{st[0]} × {st[1]}** | | | |")
        for c in rows:
            r = primary[(c, "library")]
            e, t = r["gdna_frac_est"], r["gdna_frac_true"]
            lines.append(f"| `{c}` | {e:.4f} | {t:.4f} | {e - t:+.4f} |")
    return lines + [""]


def _md_axis(primary: dict, conds: list[str], axis: str, title: str, blurb: str) -> list[str]:
    """Section 2 or 3: the transcript or the gene table, one row per scenario."""
    unit = "expressed" if axis == "transcript" else "expressed genes"
    head = (f"| scenario | {unit} | detected | Σ\\|Δ\\| | Σ\\|Δ\\| as % of true | net Δ | "
            "over-assigned | under-assigned | false-positive mass | on n | false-negative mass | on n | "
            "MARD | Spearman |")
    if axis == "transcript":
        head = head.replace("| MARD |", "| median rel. err | MARD |")
    # count columns with the escaped pipes removed: the two `|` in `Σ\|Δ\|` are cell content, not delimiters
    ncols = head.replace("\\|", "").count("|") - 1
    lines = [f"## {title}", "", blurb, "", head, "|---" + "|---:" * (ncols - 1) + "|"]
    for st, rows in _by_stratum(conds):
        lines.append(f"| **{_stratum_name(st)}** |" + " |" * (ncols - 1))
        for c in rows:
            r = primary[(c, axis)]
            share = 100.0 * r["count_abs_err"] / r["count_true"] if r["count_true"] > 0 else float("nan")
            cells = [f"`{c}`", f"{r['n_expressed']:,}", f"{r['n_detected']:,}",
                     f"{r['count_abs_err']:,.0f}", f"{share:.2f} %", _signed(r["count_net_err"])]
            cells += [f"{r['count_over']:,.0f}", f"{r['count_under']:,.0f}",
                      f"{r['fp_mass']:,.0f}", f"{r['fp_n']:,}", f"{r['fn_mass']:,.0f}", f"{r['fn_n']:,}"]
            if axis == "transcript":
                cells.append(f"{r['median_rel_err']:.3f}")
            cells += [f"{r['mard']:.3f}", f"{r['spearman']:.4f}"]
            lines.append("| " + " | ".join(cells) + " |")
    return lines + [""]


def _md_rollup(primary: dict, conds: list[str]) -> list[str]:
    """Section 4: each stratum summed over its gDNA levels, never across strata."""
    lines = ["## 4. Per stratum, summed over its gDNA levels", "",
             "⛔ Summed WITHIN a stratum only. The strata are never added together.", "",
             "| stratum | annotated truth | transcript Σ\\|Δ\\| | % | gene Σ\\|Δ\\| | % | nascent est | "
             "nascent truth | gDNA est | gDNA truth |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for st, rows in _by_stratum(conds):
        tx = sum(primary[(c, "transcript")]["count_abs_err"] for c in rows)
        gn = sum(primary[(c, "gene")]["count_abs_err"] for c in rows)
        tt = sum(primary[(c, "transcript")]["count_true"] for c in rows)
        ne = sum(primary[(c, "library")]["nascent_est"] for c in rows)
        nt = sum(primary[(c, "library")]["nascent_true"] for c in rows)
        ge = sum(gdna_estimate(primary[(c, "library")]) for c in rows)
        gt = sum(primary[(c, "library")]["gdna_true"] for c in rows)
        lines.append(f"| {_stratum_name(st, ' ')} | {tt:,.0f} | {tx:,.0f} | {100 * tx / tt:.2f} % | {gn:,.0f} | "
                     f"{100 * gn / tt:.2f} % | {ne:,.0f} | {nt:,.0f} | {ge:,.0f} | {gt:,.0f} |")
    return lines + [""]


def markdown_report(path: Path, out: Path) -> None:
    """Write the per-scenario accuracy report as markdown, the report a release is judged on.

    Every fragment is gDNA, synthetic nascent RNA or annotated RNA, split on ``is_synthetic`` (a span this index
    made), never on ``is_single_exon``: a single-exon annotated transcript belongs to the annotated pool. The
    transcript and gene sections score the annotated pool alone, the table a user reads; synthetic mass is on the
    pool rows. The gDNA estimate is ``gdna_em + n_intergenic_gdna``, since intergenic fragments never enter the EM
    but the truth counts them.
    """
    primary, pname = _load(path), Path(path).stem
    conds = sorted({c for c, _ax in primary})
    lines = [
        *_md_header(primary, pname, conds, path),
        *_md_pools(primary, conds),
        *(line for axis, title, blurb in _MD_AXES for line in _md_axis(primary, conds, axis, title, blurb)),
        *_md_rollup(primary, conds),
        "---",
        "",
        f"Generated by `scripts/design/quant_accuracy.py --markdown` from `{pname}.jsonl`.",
        "",
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines))
    print(f"  ⭐ markdown report -> {out}  ({len(lines)} lines, {len(conds)} conditions)")


#: the library row's pool fields; an arm file without them has no per-condition pool table
_POOL_FIELDS = ("gdna_est", "gdna_true", "nascent_est", "nascent_true", "annotated_est", "annotated_true", "n_intergenic_gdna")


def _check_comparable(arms: list[tuple[dict, str]]) -> None:
    """Refuse arms scored under different assignment modes or config stamps, over different row sets, or on different
    inputs for one condition."""
    modes = {(name, r.get("assignment_mode")) for a, name in arms for r in a.values()}
    if len({m for _n, m in modes}) > 1:
        # `--set em.assignment_mode=fractional` moves every count, so a mix would report the mode's effect as
        # an arm's; a row with no stamp reads as its own mode, not as agreement
        raise SystemExit(
            f"⛔ the arms were scored under different assignment modes {sorted(modes, key=str)} — "
            "compare arms run under one mode."
        )
    stamps = {(name, json.dumps([r.get("settings"), r.get("message_policy")], sort_keys=True))
              for a, name in arms for r in a.values()}
    if len({st for _n, st in stamps}) > 1:
        # a pair is one difference, the arm: a second config under the same comparison would be reported as the
        # arm's effect. A row with no stamp reads as its own config, not as agreement.
        raise SystemExit(
            f"⛔ the arms carry different config stamps (--set overrides, message policy): "
            f"{sorted(stamps, key=str)} — compare arms scored under one config."
        )
    keys = set(arms[0][0])
    for a, name in arms[1:]:
        if set(a) != keys:
            raise SystemExit(
                f"⛔ {name} has a different row set ({len(a)} vs {len(keys)}) — an arm that is "
                "missing conditions must not be aggregated against one that is not."
            )
    # one input per condition: an aligner's rows beside the simulator's would report the alignment as the arm's
    # effect. A row with no stamp reads as its own input, not as agreement.
    inputs = {(r["condition"], name, r.get("bam")) for a, name in arms for r in a.values()}
    by_condition: dict = {}
    for condition, name, bam in inputs:
        by_condition.setdefault(condition, set()).add(bam)
    mixed = sorted(c for c, bams in by_condition.items() if len(bams) > 1)
    if mixed:
        raise SystemExit(
            f"⛔ the arms read different inputs for {len(mixed)} condition(s), e.g. {mixed[0]}: "
            f"{sorted(by_condition[mixed[0]], key=str)} — compare arms scored on one input."
        )


def _stamps(rows) -> tuple[str, str, str]:
    """``(message policy, --set overrides, revisions)`` as the rows stamp them, each value they carry once."""
    def one(field, render):
        vals = sorted({render(r.get(field)) for r in rows})
        return " / ".join(vals)

    return (
        one("message_policy", lambda v: f"message_policy={v!r}" if v is not None else "unstamped"),
        one("settings", lambda v: " ".join(f"{k}={x}" for k, x in sorted(v.items())) or "the shipped defaults"
            if v is not None else "unstamped"),
        one("revision", lambda v: str(v) if v is not None else "unstamped"),
    )


def _stamp_lines(arms) -> list[str]:
    """Per arm, the config its rows were scored under and the code they ran on, read from the rows."""
    out = []
    for a, name in arms:
        policy, settings, revision = _stamps(a.values())
        out.append(f"{name}: calibration {policy} · --set {settings} · revision {revision}")
    return out


def _print_sums(arms, conds, title: str, field: str, axis: str = "transcript") -> None:
    """One summed quantity per stratum and arm, with the last arm over the first beside it."""
    print()
    print(f"  {title}")
    head = f"    {'stratum':<26}"
    for _a, name in arms:
        head += f" {name[:16]:>16}"
    if len(arms) > 1:
        head += f" {'arm/base':>10}"
    print(head)
    print("    " + "-" * (26 + 17 * len(arms) + 11))
    for label, sel in stratum_selections(conds):
        vals = []
        for a, _n in arms:
            # a field a row never recorded reads nan, never 0
            sub = [a[(c, axis)].get(field, float("nan")) for c in conds if sel(c) and (c, axis) in a]
            vals.append(sum(sub) if sub else float("nan"))
        row = f"    {label:<26}"
        for v in vals:
            row += " " + f"{v:>14,.0f}".rjust(16)
        if len(arms) > 1:
            row += f" {vals[-1] / vals[0]:>10.3f}" if vals[0] else f" {'—':>10}"
        print(row)


def _print_rank(arms, conds) -> None:
    """⑤: the rank statistics, which do not sum, as a mean over conditions."""
    print()
    print("  ⑤ RANK / SHAPE (a mean over conditions — these do not add)")
    head = f"    {'stratum':<26}"
    for _a, name in arms:
        head += f" {name[:10] + ' sprmn':>17}{name[:10] + ' mard':>17}"
    print(head)
    print("    " + "-" * (26 + 34 * len(arms)))
    for label, sel in stratum_selections(conds):
        row = f"    {label:<26}"
        for a, _n in arms:
            sub = [a[(c, "transcript")] for c in conds if sel(c)]
            sp = float(np.nanmean([x["spearman"] for x in sub])) if sub else float("nan")
            md = float(np.nanmean([x["mard"] for x in sub])) if sub else float("nan")
            row += f" {sp:>16.4f} {md:>16.4f}"
        print(row)


def _print_pools(arms, conds) -> None:
    """⑥: where the library's fragments went, per stratum and arm.

    "Nascent" is a simulator input, not a solver population: `nascent_est` is mass the EM assigned to the synthetic
    entities, so this table scores assignment. On an `nascent_none` panel `nascent_true` is 0 and every fragment in that
    column is a false positive with nothing to cancel it, so the column is printed even there."""
    print()
    print("  ⑥ ⭐⭐⭐ POOL LEVEL — where the library's fragments went: gDNA vs NASCENT vs ANNOTATED")
    print("     ⛔ 'nascent' is a SIMULATOR input, never a population the solver has (AXIOM 0). This")
    print("        scores ASSIGNMENT. On a nascent_none panel nascent truth is 0, so that column is pure")
    print("        false positive and Δ% is undefined rather than large.")
    for a, name in arms:
        print()
        print(f"    arm: {name}")
        print(f"    {'stratum':<26} {'pool':<20} {'est':>15} {'true':>15} {'Δ':>15} {'Δ%':>8}")
        print("    " + "-" * 102)
        for label, sel in stratum_selections(conds):
            _print_pool_rows(a, conds, label, sel)


def _print_pool_rows(a: dict, conds, label: str, sel) -> None:
    def tot(field):
        return sum(a[(c, "library")][field] for c in conds if sel(c) and (c, "library") in a)

    # The gDNA total includes the intergenic fragments: `gdna_est` is `gdna_em_count`, which excludes fragments
    # that reached no locus, while `gdna_true` includes them. The `of which` rows carry no truth; only the total is
    # scored, and it is `cli.py`'s `count_gdna`, the `fraction_gdna` numerator.
    gdna_em, inter = tot("gdna_est"), tot("n_intergenic_gdna")
    rows = [
        ("gDNA (total)", gdna_em + inter, tot("gdna_true")),
        ("  of which EM", gdna_em, None),
        ("  of which intergenic", inter, None),
        ("nascent", tot("nascent_est"), tot("nascent_true")),
        ("annotated", tot("annotated_est"), tot("annotated_true")),
    ]
    for i, (pretty, est, tru) in enumerate(rows):
        first = label if i == 0 else ""
        if tru is None:
            print(f"    {first:<26} {pretty:<20} {est:>15,.0f} {'—':>15} {'—':>15} {'—':>8}")
            continue
        d = est - tru
        # a relative error against a true zero is undefined, not large
        pct = f"{100.0 * d / tru:>7.1f}%" if tru > 0 else "     n/a"
        print(f"    {first:<26} {pretty:<20} {est:>15,.0f} {tru:>15,.0f} {d:>+15,.0f} {pct:>8}")


def _print_pools_per_condition(arms, conds) -> None:
    """⑦: the same three pools, intergenic-inclusive gDNA, one line per condition: the stratum roll-up says which
    stratum is off, this says what happened in each scenario."""
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
            print(f"    {c:<44} {gdna_estimate(r):>13,.0f} {r['gdna_true']:>13,.0f} "
                  f"{r['nascent_est']:>11,.0f} {r['nascent_true']:>10,.0f} "
                  f"{r['annotated_est']:>13,.0f} {r['annotated_true']:>13,.0f}")


def _print_thermometer(arms, conds) -> None:
    """⑧: the library gDNA fraction per condition, truth beside every arm."""
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


def _print_errors_per_condition(arms, conds) -> None:
    """⑨: misassigned fragments and false-positive mass per condition, the deliverable."""
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


def report(paths: list[Path]) -> None:
    """Print one or more arms' tables per stratum, never pooled: a panel total hides a sign flip between strata."""
    arms = [(_load(p), Path(p).stem) for p in paths]
    _check_comparable(arms)
    conds = sorted({c for c, _ax in arms[0][0]})

    print()
    print("=" * 112)
    print("  ⭐⭐⭐ TRANSCRIPT-LEVEL ACCURACY, end to end, scored against the REALISED fragment truth")
    print(f"  {len(conds)} conditions   arms: {', '.join(n for _a, n in arms)}")
    for line in _stamp_lines(arms):
        print(f"  {line}")
    print("=" * 112)

    _print_sums(arms, conds, "① TOTAL MISASSIGNED FRAGMENTS  Σ|count_est − count_true|  ·  TRANSCRIPT level",
                "count_abs_err")
    if any((conds[0], "gene") in a for a, _n in arms):
        _print_sums(arms, conds, "①g   … the SAME quantity at GENE level — isoform ambiguity summed away, so what is "
                    "left is\n       error in deciding whether the fragment was RNA from this gene at all",
                    "count_abs_err", axis="gene")
    _print_sums(arms, conds, "② ⭐ FALSE-POSITIVE MASS — estimate on transcripts the simulator gave ZERO fragments",
                "fp_mass")
    _print_sums(arms, conds, "③ FALSE-NEGATIVE MASS — real RNA the tool assigned nowhere", "fn_mass")
    _print_sums(arms, conds, "④ TPM error, Σ|Δ| out of 1e6 (assignment only — see the docstring)", "tpm_abs_err")
    _print_sums(arms, conds, "④n   … the true fragments on transcripts with no effective length, which carry no TPM "
                "on either side (a shadow transcript the tool does not hold)", "no_length_fragments")
    _print_rank(arms, conds)
    _print_pools(arms, conds)
    _print_pools_per_condition(arms, conds)
    _print_thermometer(arms, conds)
    _print_errors_per_condition(arms, conds)


# ── main ─────────────────────────────────────────────────────────────────────────────────────────


def _parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", nargs="+", type=Path, default=None,
                    help="print the per-stratum tables from arm jsonl files and exit")
    ap.add_argument("--markdown", type=Path, default=None,
                    help="with --report: also write the full per-scenario report as markdown, of the "
                         "first --report file")
    ap.add_argument("--arm", choices=ARMS, default=None)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--bam", default=ORACLE_BAM,
                    help="the BAM each condition is read from, by name in its directory: the simulator's own "
                         "alignments by default, or an aligner's (panel.py align), which only base and no_capture "
                         "score")
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=1)
    add_set_option(ap, "a config value applied to every arm, repeatable (e.g. --set "
                       "em.assignment_mode=fractional); the same parser as every instrument's --set")
    return ap.parse_args()


def _build_capture_truth(suite: Path, index_dir: Path, names: list[str]) -> None:
    """The ruler arms' capture truth, once per capture label: it takes minutes per label and every condition
    shares it, so it is built here, before shards could race to build it."""
    truth_index = TranscriptIndex.load(str(index_dir))
    for label in sorted({capture_label(n) for n in names}):
        first = next(n for n in names if capture_label(n) == label)
        capture_truth_factor(suite, truth_index, first)


def _run_sharded(args, names: list[str], cache: Path | None) -> int:
    outs = run_shards(
        Path(__file__).resolve(), names, args.jobs, args.out.parent / f".{args.out.stem}_shards",
        lambda out, shard: ["--arm", args.arm, "--suite", str(args.suite), "--index", str(args.index),
                            "--out", str(out), "--bam", args.bam, "--conditions", *shard,
                            *([] if cache is None else ["--oracle-cache", str(cache)]),
                            *set_args(args.settings)],
    )
    with args.out.open("w") as fh:
        for o in outs:
            fh.write(o.read_text())
    print(f"  ⭐ {args.arm}: {sum(1 for _ in args.out.open())} rows -> {args.out}")
    return 0


def _run_serial(args, names: list[str], cache: Path | None) -> int:
    index = TranscriptIndex.load(str(args.index))
    pipeline_config = configured(args.settings)
    for spec in args.settings:
        print(f"  ⭐ --set {spec}", flush=True)
    rows = []
    for name in names:
        print(f"  … {args.arm}  {name}", flush=True)
        rows += run_condition(args.arm, args.suite, index, name, pipeline_config, cache,
                              bam_name=args.bam)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"  ⭐ {args.arm}: {len(rows)} rows -> {args.out}")
    return 0


def main() -> int:
    args = _parse_args()
    if args.report:
        report(args.report)
        if args.markdown:
            markdown_report(args.report[0], args.markdown)
        return 0
    if args.markdown:
        raise SystemExit("--markdown needs --report FILES... (it renders arm jsonl, it runs nothing)")
    if args.arm is None or args.out is None:
        raise SystemExit("--arm and --out are required (or use --report)")

    names = panel_conditions(args.suite, args.conditions)
    check_input_bam(args.suite, names, args.arm, args.bam)
    cache = oracle_cache_dir(args.suite, args.oracle_cache)
    if args.arm in _RULER_ARMS:
        _build_capture_truth(args.suite, args.index, names)
    if args.jobs > 1 and len(names) > 1:
        return _run_sharded(args, names, cache)
    return _run_serial(args, names, cache)


if __name__ == "__main__":
    raise SystemExit(main())
