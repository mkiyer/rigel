#!/usr/bin/env python
"""Is the ruler's formula right? The EM's effective length under capture, scored per transcript against the
simulator's own capture-aware effective length.

The truth is what generated the reads: a transcript's true captured length is ``Σ_w f_pre(w) · partition_t(w)``
(``CaptureSampler.partition_array``) and its plain length ``Σ_w f_pre(w) · off_target_weight · (L_t − w + 1)+``.
The truth factor is their ratio; an arm's factor is its effective length over the same plain lengths. The two are
compared in the log, anchored on the fully probed transcripts (a global scale is free), and summarised per probed
class and per kind (mRNA, annotated single-exon, synthetic span). A factor of exactly 0 is counted and dropped from
the quantiles.

Arms: ``shipped`` is the ruler in ``src/``; ``oracle_gdna`` feeds it the certified true gDNA counts per object
(``slot_truth.npz``) under a landscape fitted on the truth, at the capture weight the true counts give, so the gap
between the two is calibration's and what remains is the ruler's. ``--module FILE`` adds prototype rulers, ``ARMS = {name: ruler}`` with
``ruler(calibration, region_arrays, index, fl, **inputs)`` returning effective lengths, where ``inputs`` is
``calibrate_condition``'s debug bundle. ``--set SECTION.FIELD=VALUE`` prices a config value on every arm.

Read the classes apart: the probed class is exact when its witness is, the unprobed class is where the bulk bites,
and the partial classes carry the junction price, which reads a junction's capture from the gDNA objects beside it
(right on average, noisy per junction).

Usage::

    python scripts/design/ruler_vs_truth.py --panel test --condition gdna_g05_ss_0.99_nrna_file_capture_on
    python scripts/design/ruler_vs_truth.py --panel test                      # every capture-ON condition, one line each
    python scripts/design/ruler_vs_truth.py --panel test --panel-dir $RIGEL_RUNS/test_reference/scenarios_depth_d10
    python scripts/design/ruler_vs_truth.py --panel ladder --condition C --module proto.py --out table.tsv
    python scripts/design/ruler_vs_truth.py --panel test --condition C --set calibration.message_policy=silent
    python scripts/design/ruler_vs_truth.py --panel ladder --condition C --scale        # the one-scale read-out
    python scripts/design/ruler_vs_truth.py --self-test

``--panel-dir`` names a directory of conditions, never the index: pass ``--panel ladder`` with it for a
ladder-derived panel, or the test chromosome's index is loaded against the ladder's BAM and nothing deposits.

``--scale`` asks whether every hypothesis class's opportunities are on one scale: the E-step balances two hypotheses
only if their opportunities are the same multiple of the yield the simulator drew them with. It prints L / Y,
unanchored, for the locus gDNA component (``LocusPriors.gdna_eff_len`` against the gDNA-space yield over its
footprint), the synthetic spans and the annotated transcripts (the shipped lengths against ``partition_array``), by
probed class with the junction-probed transcripts apart, then the pairwise class-mean ratios, which must agree to
``ONE_SCALE_TOL``. The pipeline runs to the prior assembly and stops there: no EM.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pysam  # noqa: E402
from scipy.special import ndtr  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]

from _shared import (  # noqa: E402
    ORACLE_BAM,
    PANELS,
    add_set_option,
    capture_config,
    configured,
    refuse_unconsumed,
    refuses,
    sibling,
)
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.capture_eff_length import (  # noqa: E402
    fit_junction_cap,
    transcript_capture_eff_lengths,
)
from rigel.calibration.capture_efficiency import capture_efficiencies, visible_share, weighted  # noqa: E402
from rigel.calibration.diagnostics import CalibrationDiagnostics  # noqa: E402
from rigel.calibration.landscape import fit_landscape  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_geometry import region_gdna_geometry  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import length_tables, library_fl_models  # noqa: E402
from _scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from rigel.sim.capture import CaptureConfig, CaptureSampler  # noqa: E402
from rigel.sim.transcriptome import load_transcripts_from_index  # noqa: E402
from rigel.sim.wgs_config import GDNASimConfig  # noqa: E402


#: the transcript classes, by the probed fraction of the transcript's bases
CLASSES = (
    ("unprobed", 0.0, 1e-9),
    ("partial ≤ ½", 1e-9, 0.5),
    ("partial > ½", 0.5, 0.9),
    ("probed ≥ 0.9", 0.9, 1.01),
)
KINDS = ("mRNA", "single-exon", "entity")
#: the tolerance the per-class share is reported at: within ±0.1 nat of the truth, a 10 % effective length
WITHIN = 0.1
#: the keys a capture config carries from the panel manifest: every field of the simulator's config
_CAPTURE_KEYS = tuple(f.name for f in dataclasses.fields(CaptureConfig))


# ── the truth ────────────────────────────────────────────────────────────────────────────────────


def fl_pmf(sim) -> tuple[np.ndarray, np.ndarray]:
    """The exact fragment-length pmf the simulator draws from, for *sim*'s ``frag_mean``, ``frag_std``, ``frag_min``
    and ``frag_max``: `rigel.sim.sampling.truncated_normal_frag_lengths` rounds a Normal draw and keeps it in
    ``[frag_min, frag_max]``, so length ``k`` has the Normal's mass on ``[k − ½, k + ½)``, renormalised over the
    window. Each mass is taken on the side of the mean where it is a difference of two small tails, never of two
    numbers near 1."""
    widths = np.arange(int(sim.frag_min), int(sim.frag_max) + 1, dtype=np.int64)
    lo = (widths - 0.5 - float(sim.frag_mean)) / float(sim.frag_std)
    hi = lo + 1.0 / float(sim.frag_std)
    p = np.where(lo > 0.0, ndtr(-lo) - ndtr(-hi), ndtr(hi) - ndtr(lo))
    return widths, p / p.sum()


@dataclasses.dataclass(frozen=True, slots=True)
class Truth:
    """Per transcript: the plain (uncaptured) length, the true factor, the probed fraction of its bases,
    its kind, and the condition's observed RNA fragments and pre-capture abundance."""

    L_plain: np.ndarray
    factor: np.ndarray
    probed_frac: np.ndarray
    kind: np.ndarray
    frags: np.ndarray
    abundance: np.ndarray


def sampler_truth(sampler: CaptureSampler, lengths, widths, pmf, off_target: float) -> tuple:
    """``(L_plain, factor_true, probed_frac)`` from the sampler's own partition.

    ``L_plain`` counts the starts the simulator could draw without capture, so a transcript's unknown bases (masked,
    never drawn) leave its factor the capture's alone. The probed fraction is geometry, blind to how strongly each
    probe binds: the share of its ``n`` known bases that lie under a probe part (``CaptureSampler.probed_bases``),
    a split probe's pieces counted base by base."""
    L = np.asarray(lengths, dtype=np.int64)
    n = L.size
    rows = range(n)
    L_true = np.zeros(n)
    L_plain = np.zeros(n)
    for w, pw in zip(widths, pmf, strict=True):
        L_true += pw * sampler.partition_array("mrna", rows, L, int(w))
        L_plain += pw * off_target * sampler.admissible_counts("mrna", rows, L, int(w))
    ok = L_plain > 0
    factor = np.where(ok, L_true / np.maximum(L_plain, 1e-12), np.nan)
    binding = float(sampler.config.binding_per_base)
    if sampler.enabled and binding > 0.0:
        known = sampler.admissible_counts("mrna", rows, L, 1)
        probed = np.clip(sampler.probed_bases("mrna", rows) / np.maximum(known, 1), 0.0, 1.0)
    else:
        probed = np.zeros(n)
    return L_plain, factor, probed


@dataclasses.dataclass(frozen=True, slots=True)
class Sim:
    """The condition's simulator as the truth needs it: the sampler on the panel's probes, the RNA and gDNA generating
    fragment-length laws, the two weight constants, the reference lengths, the transcript rows, and the references
    gDNA is drawn from (the manifest's ``gdna.genomic_refs``; an RNA-only spike-in backbone is not among them)."""

    sampler: CaptureSampler
    widths: np.ndarray
    pmf: np.ndarray
    off: float
    bind: float
    ref_lengths: dict
    transcripts: list
    gdna_widths: np.ndarray
    gdna_pmf: np.ndarray
    gdna_refs: frozenset


def load_sim(index: TranscriptIndex, index_dir: Path, panel_dir: Path, condition: str) -> Sim:
    """The condition's simulator, read from the panel manifest: its capture config on the panel's probes, each probe
    weighed by its sequence in the panel's genome, masked where that genome is unknown, and its RNA and gDNA
    fragment-length laws. A captured condition whose recorded capture config lacks a field of the simulator's, or
    carries one it does not know, was simulated under another weight model, and is refused."""
    manifest = json.load(open(panel_dir / "manifest.json"))
    cap_cfg = capture_config(manifest, condition)
    missing = sorted(set(_CAPTURE_KEYS) - set(cap_cfg))
    unknown = sorted(set(cap_cfg) - set(_CAPTURE_KEYS))
    if cap_cfg.get("probes") and (missing or unknown):
        raise SystemExit(
            f"⛔ {condition}'s capture config records {missing or 'every field'} missing and {unknown or 'no'} unknown "
            "field(s): it was simulated under another capture weight model than this simulator's, so its truth "
            "cannot be rebuilt here; simulate the panel again"
        )
    sim = manifest["simulation"]
    transcripts = load_transcripts_from_index(index_dir)
    if len(transcripts) != len(index.t_df) or any(
        t.t_id != i for t, i in zip(transcripts, index.t_df.t_id, strict=True)
    ):
        raise SystemExit("⛔ the simulator's transcript rows do not match the index's")
    ref_lengths = {
        r["ref"]: int(r["length"])
        for r in pd.read_feather(index.index_dir + "/ref_lengths.feather").to_dict("records")
    }
    # the probes' strengths and the masks the reads were drawn under, from the panel's own genome
    if "genome" not in manifest:
        raise SystemExit("⛔ the manifest names no genome, so the truth cannot mask what the simulator never drew")
    with pysam.FastaFile(str(manifest["genome"])) as fasta:
        sampler = CaptureSampler.from_config(
            CaptureConfig(**{k: cap_cfg[k] for k in _CAPTURE_KEYS if k in cap_cfg}),
            transcripts,
            ref_lengths,
            fasta=fasta,
        )
        sampler.mask_unknown_sequence(fasta, manifest.get("gdna", {}).get("genomic_refs") or [])
    params = type("Sim", (), {k: sim[k] for k in ("frag_mean", "frag_std", "frag_min", "frag_max")})
    widths, pmf = fl_pmf(params)
    gdna_defaults = GDNASimConfig()
    gdna = manifest.get("gdna", {})
    gdna_params = type("GDNA", (), {
        k: gdna.get(k, getattr(gdna_defaults, k))
        for k in ("frag_mean", "frag_std", "frag_min", "frag_max")
    })
    gdna_widths, gdna_pmf = fl_pmf(gdna_params)
    return Sim(
        sampler,
        np.asarray(widths, dtype=np.int64),
        np.asarray(pmf, dtype=np.float64),
        float(cap_cfg["off_target_weight"]),
        float(cap_cfg.get("binding_per_base", 0.0)),
        ref_lengths,
        transcripts,
        np.asarray(gdna_widths, dtype=np.int64),
        np.asarray(gdna_pmf, dtype=np.float64),
        frozenset(manifest.get("gdna", {}).get("genomic_refs") or []),
    )


def transcript_kinds(t_df: pd.DataFrame) -> np.ndarray:
    """Per transcript, one of :data:`KINDS`: a synthetic span (``entity``), an annotated single-exon transcript, or
    an mRNA."""
    return np.where(
        t_df.is_synthetic.to_numpy(bool),
        "entity",
        np.where(t_df.is_single_exon.to_numpy(bool), "single-exon", "mRNA"),
    )


def load_truth(index: TranscriptIndex, index_dir: Path, panel_dir: Path, condition: str, sim: Sim | None = None) -> Truth:
    """The condition's per-transcript truth, from its simulator (``sim``, loaded here unless given) and its truth
    table."""
    if sim is None:
        sim = load_sim(index, index_dir, panel_dir, condition)
    L = np.array([int(t.length) for t in sim.transcripts], dtype=np.int64)
    L_plain, factor, probed = sampler_truth(sim.sampler, L, sim.widths, sim.pmf, sim.off)
    t = index.t_df
    truth = pd.read_csv(panel_dir / condition / "truth_abundances.tsv", sep="\t").set_index(
        "transcript_id"
    )
    frags = truth["observed_total_rna_fragments"].reindex(t.t_id).fillna(0).to_numpy(float)
    abund = truth["pre_capture_total_rna"].reindex(t.t_id).fillna(0).to_numpy(float)
    return Truth(L_plain, factor, probed, transcript_kinds(t), frags, abund)


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def calibrate_condition(index, region_arrays, panel_dir: Path, condition: str, config):
    """The shipped calibration on the condition's cached payload. Returns the result and the solve's debug bundle:
    the fields of ``CalibrationDiagnostics`` (the last refit's landscape under ``landscape``, ``None`` without a
    refit; the chain; the belief), the two fragment-length pmfs, the EM's length tables and the scan cache's
    provenance (``scan_provenance``), which keys the certified truth."""
    cache = read_scan_cache(panel_dir / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    diagnostics = CalibrationDiagnostics()
    cal = calibrate(config=config.calibration, diagnostics=diagnostics, **kw)
    debug = {f.name: getattr(diagnostics, f.name) for f in dataclasses.fields(diagnostics)}
    # the scan the solve read, which the certified truth must have been derived from
    debug["scan_provenance"] = cache.provenance
    debug["gdna_fl_pmf"] = kw["gdna_fl_pmf"]
    debug["rna_fl_pmf"] = kw["rna_fl_pmf"]
    # the EM's two length tables as production builds them on this payload: the shipped ruler reads them
    debug["length_tables"] = length_tables(library_fl_models(kw["payload"], index), cal, region_arrays)
    return cal, debug


def oracle_gdna_result(cal, region_arrays, panel_dir: Path, condition: str, inputs: dict):
    """The same result with the ruler's gDNA masses replaced by the certified truth per object, read by the
    shipped reader under a landscape fitted on the truth and raised to the weight the truth gives (the visible share
    of the true counts over the solve's chain, ``inputs`` being ``calibrate_condition``'s bundle): the ideal witness,
    with no estimator noise."""
    truth = sibling("_oracle_arms.py").certified_gdna(panel_dir / "oracle_cache" / condition, inputs["scan_provenance"])
    m_reg, m_bnd = truth["region"], truth["boundary"]
    chain, statics = inputs["chain"], inputs["statics"]
    at = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    booked = np.empty(at.shape[0])
    booked[at] = m_reg[idx[at]]
    booked[~at] = m_bnd[idx[~at]]
    both = np.asarray(statics.free_pos, dtype=bool) & np.asarray(statics.free_neg, dtype=bool)
    weight = visible_share(booked, region_gdna_geometry(inputs["geometry"])[0], ~both)
    S = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)
    rtype = coarse_type_array(np.asarray(region_arrays.signature))
    live = S > 1e-9
    anchor = live & (m_reg <= 1e-12) & (rtype != RegionType.EXON)
    ls = fit_landscape(
        m_reg[live],
        np.maximum(m_reg[live], 1.0),
        S[live],
        np.zeros(int(live.sum())),
        anchor=anchor[live],
    )
    # the reader takes the true supports, as calibrate does: a region with none reads 1
    reading = (
        None
        if ls is None
        else weighted(
            capture_efficiencies(
                ls,
                m_reg,
                np.asarray(cal.gdna_region_eff_len, float),
                m_bnd,
                np.asarray(cal.gdna_boundary_eff_len, float),
            ),
            weight,
        )
    )
    return dataclasses.replace(
        cal,
        count_gdna_region=m_reg,
        count_gdna_boundary=m_bnd,
        gdna_bulk_density=None if reading is None else reading.bulk_density,
        gdna_capture_weight=weight,
        gdna_capture_saturation=1.0 if reading is None else reading.saturation,
        gdna_capture_efficiency_region=(
            np.ones(m_reg.shape[0]) if reading is None else reading.efficiency_region
        ),
        gdna_capture_efficiency_boundary=(
            np.ones(m_bnd.shape[0]) if reading is None else reading.efficiency_boundary
        ),
    )


def load_arms(module_path):
    """``ARMS = {name: ruler}`` from a prototype file; each ruler takes
    ``(calibration, region_arrays, index, fl, **inputs)`` and returns the effective lengths."""
    if module_path is None:
        return {}
    spec = importlib.util.spec_from_file_location("ruler_vs_truth_arms", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    arms = getattr(module, "ARMS", None)
    if not isinstance(arms, dict) or not arms or not all(callable(v) for v in arms.values()):
        raise SystemExit(f"⛔ {module_path} must define ARMS = {{name: ruler}} with callable rulers")
    if "shipped" in arms or "oracle_gdna" in arms:
        raise SystemExit("⛔ a prototype arm may not be named `shipped` or `oracle_gdna`")
    return dict(arms)


# ── the score ────────────────────────────────────────────────────────────────────────────────────


def score(factor_arm, truth: Truth, min_frags: int) -> tuple[np.ndarray, np.ndarray, int]:
    """``(d, measured, n_zero)``: the anchored log error per transcript, the transcripts it is read on
    (a finite truth and error, at least ``min_frags`` observed RNA fragments), and how many measured
    transcripts read a factor of exactly 0. The anchor is the median error of the probed class (its own
    transcripts with at least ``min_frags``), the whole population when no probed transcript qualifies."""
    f = np.asarray(factor_arm, dtype=np.float64)
    ok = np.isfinite(truth.factor) & (truth.L_plain > 0) & (truth.frags >= min_frags)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.log(f) - np.log(truth.factor)
    n_zero = int(np.sum(ok & (f <= 0.0)))
    measured = ok & np.isfinite(d)
    anchor = measured & (truth.probed_frac >= 0.9)
    pool = anchor if anchor.any() else measured
    scale = float(np.median(d[pool])) if pool.any() else 0.0
    return d - scale, measured, n_zero


def class_rows(d, measured, truth: Truth) -> list[dict]:
    rows = []
    for cname, lo, hi in CLASSES:
        for k in KINDS:
            m = measured & (truth.probed_frac >= lo) & (truth.probed_frac < hi) & (truth.kind == k)
            if not m.any():
                continue
            q = np.quantile(d[m], [0.1, 0.5, 0.9])
            la = np.log(np.maximum(truth.abundance[m], 1e-12))
            slope = (
                float(np.polyfit(la, d[m], 1)[0]) if m.sum() >= 3 and la.std() > 0 else float("nan")
            )
            rows.append(
                dict(
                    cls=cname,
                    kind=k,
                    n=int(m.sum()),
                    med=float(q[1]),
                    q10=float(q[0]),
                    q90=float(q[2]),
                    within=float(np.mean(np.abs(d[m]) <= WITHIN)),
                    truth_factor=float(np.median(truth.factor[m])),
                    slope=slope,
                )
            )
    return rows


def headline(d, measured, truth: Truth) -> dict:
    """The one-line summary a condition sweep prints: the probed class's share within the tolerance
    and the unprobed and partial classes' median error, mRNA and single-exon transcripts pooled."""
    real = measured & (truth.kind != "entity")
    probed = real & (truth.probed_frac >= 0.9)
    unprobed = real & (truth.probed_frac < 1e-9)
    partial = real & (truth.probed_frac >= 1e-9) & (truth.probed_frac < 0.9)
    return dict(
        n_probed=int(probed.sum()),
        probed_within=float(np.mean(np.abs(d[probed]) <= WITHIN)) if probed.any() else float("nan"),
        n_unprobed=int(unprobed.sum()),
        unprobed_med=float(np.median(d[unprobed])) if unprobed.any() else float("nan"),
        n_partial=int(partial.sum()),
        partial_med=float(np.median(d[partial])) if partial.any() else float("nan"),
    )


def _fmt_capture(cal) -> str:
    if cal.gdna_bulk_density is None:
        return "None (factor 1 everywhere)"
    return f"bulk {cal.gdna_bulk_density:.3e}/bp, saturation {cal.gdna_capture_saturation:.4g}× the bulk"


def run_condition(index, region_arrays, index_dir, panel_dir, condition, config, arms, min_frags):
    """Every arm on one condition: ``(truth, {arm: (factor, calibration)})``."""
    truth = load_truth(index, index_dir, panel_dir, condition)
    cal, inputs = calibrate_condition(index, region_arrays, panel_dir, condition, config)
    out = {}
    tables = inputs["length_tables"]

    def shipped_ruler(result):
        # each result's lengths at its own junction cap, as the pipeline fits it
        cap = fit_junction_cap(result, region_arrays, index, tables)
        return transcript_capture_eff_lengths(
            result, region_arrays, index, truth.L_plain, tables, cap
        )

    out["shipped"] = (shipped_ruler(cal), cal)
    if (panel_dir / "oracle_cache" / condition / "slot_truth.npz").is_file():
        o = oracle_gdna_result(cal, region_arrays, panel_dir, condition, inputs)
        out["oracle_gdna"] = (shipped_ruler(o), o)
    for name, ruler in arms.items():
        out[name] = (
            np.asarray(ruler(cal, region_arrays, index, truth.L_plain, **inputs), dtype=np.float64),
            cal,
        )
    factors = {
        k: (np.where(truth.L_plain > 0, v[0] / np.maximum(truth.L_plain, 1e-12), np.nan), v[1])
        for k, v in out.items()
    }
    return truth, factors


def print_condition(condition, truth, factors, min_frags, elapsed):
    probed = truth.probed_frac >= 0.9
    unprobed = truth.probed_frac < 1e-9
    print(
        f"⭐ {condition}  ({elapsed:.0f} s) — {truth.L_plain.size:,} transcripts; truth factor: "
        f"unprobed {np.nanmedian(truth.factor[unprobed]) if unprobed.any() else float('nan'):.3g}, "
        f"probed {np.nanmedian(truth.factor[probed]) if probed.any() else float('nan'):.3g}"
    )
    for arm, (f, cal) in factors.items():
        d, measured, n_zero = score(f, truth, min_frags)
        anchor = int((measured & probed).sum())
        print(
            f"   arm {arm:<12} capture {_fmt_capture(cal)}; anchored on {anchor} probed transcripts "
            f"with ≥ {min_frags} fragments; factor exactly 0 (dropped): {n_zero}"
        )
        print(
            f"   {'class':<14}{'kind':<12}{'n':>6}   {'median log err':>14} {'10 %':>8} {'90 %':>8}"
            f"   {'within ±0.1':>11}  {'truth factor':>13}  {'slope vs log a':>14}"
        )
        for r in class_rows(d, measured, truth):
            print(
                f"   {r['cls']:<14}{r['kind']:<12}{r['n']:>6}   {r['med']:>+14.3f} {r['q10']:>+8.2f} "
                f"{r['q90']:>+8.2f}   {r['within']:>10.0%}  {r['truth_factor']:>13.3g}  {r['slope']:>+14.2f}"
            )


def _parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", choices=sorted(PANELS), default="test")
    ap.add_argument(
        "--panel-dir",
        type=Path,
        default=None,
        help="a scenario set on the panel's index (the depth ladder) instead of the panel's own",
    )
    ap.add_argument("--condition", default=None, help="one condition, the full per-class table")
    ap.add_argument("--module", type=Path, default=None, help="a .py defining ARMS = {name: ruler}")
    ap.add_argument("--min-frags", type=int, default=20, help="observed RNA fragments a scored transcript needs")
    add_set_option(ap, "override one PipelineConfig field on every arm; repeatable")
    ap.add_argument("--out", type=Path, default=None, help="the per-transcript table (one condition)")
    ap.add_argument("--scale", action="store_true", help="THE ONE-SCALE READ-OUT: L / Y per hypothesis class (one condition)")
    ap.add_argument("--width-step", type=int, default=1, help="--scale: thin the gDNA-space yield's width grid by this step")
    ap.add_argument("--self-test", action="store_true", help="perturb every comparator; no panel")
    return ap.parse_args()


def run_scale(args, index, index_dir, panel_dir, config) -> None:
    """``--scale``: the one-scale read-out for one condition, and its per-transcript and per-locus tables."""
    if not args.condition:
        raise SystemExit("⛔ --scale reads one condition: pass --condition")
    r = scale_readout(index, index_dir, panel_dir, args.condition, config, args.min_frags, args.width_step)
    print_scale(r)
    if args.out:
        pd.DataFrame(r["per_transcript"]).to_csv(args.out, sep="\t", index=False)
        pd.DataFrame(r["per_locus"]).to_csv(args.out.with_name(args.out.stem + "_loci" + args.out.suffix), sep="\t", index=False)
        print(f"   per-transcript and per-locus tables → {args.out} and its _loci twin")


def run_one(args, index, region_arrays, index_dir, panel_dir, config, arms) -> None:
    """``--condition``: every arm's per-class table on one condition, and its per-transcript table."""
    t0 = time.time()
    truth, factors = run_condition(
        index, region_arrays, index_dir, panel_dir, args.condition, config, arms, args.min_frags
    )
    print_condition(args.condition, truth, factors, args.min_frags, time.time() - t0)
    if args.out:
        table = pd.DataFrame(
            dict(
                t_id=index.t_df.t_id,
                kind=truth.kind,
                probed_frac=truth.probed_frac,
                frags=truth.frags,
                abundance=truth.abundance,
                L_plain=truth.L_plain,
                factor_true=truth.factor,
                **{f"factor_{k}": v[0] for k, v in factors.items()},
            )
        )
        table.to_csv(args.out, sep="\t", index=False)
        print(f"   per-transcript table → {args.out}")


def run_sweep(args, index, region_arrays, index_dir, panel_dir, config, arms) -> None:
    """Every cached capture-ON condition of the panel, one headline line per condition and arm."""
    conds = sorted(
        p.name
        for p in (panel_dir / "scan_cache").iterdir()
        if (p / "payload.npz").is_file() and p.name.endswith("_capture_on")
    )
    if not conds:
        raise SystemExit(f"⛔ no cached capture-ON condition under {panel_dir / 'scan_cache'}")
    names = None
    for c in conds:
        t0 = time.time()
        truth, factors = run_condition(
            index, region_arrays, index_dir, panel_dir, c, config, arms, args.min_frags
        )
        if names is None:
            names = list(factors)
            print(
                f"⭐ {panel_dir.name}: {len(conds)} capture-ON conditions — per arm, the probed class's share "
                f"within ±{WITHIN} nat and the unprobed / partial classes' median log error "
                f"(transcripts with ≥ {args.min_frags} fragments; entities excluded)"
            )
            print(f"{'condition':<44}{'bulk':>12} " + " ".join(f"{n:>30}" for n in names))
        cal = factors["shipped"][1]
        ref = (
            "None"
            if cal.gdna_bulk_density is None
            else f"10^{np.log10(cal.gdna_bulk_density):+.2f}"
        )
        cells = []
        for n in names:
            d, measured, _ = score(factors[n][0], truth, args.min_frags)
            h = headline(d, measured, truth)
            cells.append(
                f"{h['probed_within']:>5.0%}/{h['n_probed']:<4} {h['unprobed_med']:>+6.2f}/{h['n_unprobed']:<3} "
                f"{h['partial_med']:>+6.2f}/{h['n_partial']:<3}"
            )
        print(f"{c[-44:]:<44}{ref:>12} " + " ".join(f"{x:>30}" for x in cells) + f"   ({time.time() - t0:.0f} s)")
    print("   columns per arm: probed within ±0.1 / n · unprobed median / n · partial median / n")


def consumed_sections(scale: bool) -> frozenset[str]:
    """The `PipelineConfig` sections a run hands to the code it measures: the arms calibrate a cached scan, and
    ``--scale`` runs the pipeline from the BAM to the prior assembly, through the scan and the fragment scoring,
    and stops before the EM."""
    return frozenset({"scan", "scoring", "calibration"} if scale else {"calibration"})


def main() -> int:
    args = _parse_args()
    if args.self_test:
        return self_test()

    refuse_unconsumed(args.settings, consumed_sections(args.scale), "ruler_vs_truth.py")
    index_dir, panel_dir = PANELS[args.panel]
    if args.panel_dir is not None:
        panel_dir = args.panel_dir
    config = configured(args.settings)
    arms = load_arms(args.module)
    index = TranscriptIndex.load(str(index_dir))
    region_arrays = RegionArrays.from_index(index)
    if args.scale:
        run_scale(args, index, index_dir, panel_dir, config)
    elif args.condition:
        run_one(args, index, region_arrays, index_dir, panel_dir, config, arms)
    else:
        run_sweep(args, index, region_arrays, index_dir, panel_dir, config, arms)
    return 0


# ── the one-scale read-out ───────────────────────────────────────────────────────────────────────

#: the class-mean tolerance the read-out judges against: the synthetic pool is highly elastic to its own
#: opportunity, so even a 1 % disagreement between two classes' means moves that pool far more than 1 %
ONE_SCALE_TOL = 0.01


def em_loci(index: TranscriptIndex, panel_dir: Path, condition: str, config):
    """Run the pipeline to the prior assembly and stop there, with no EM. Returns the shipped calibration, the
    region arrays, the EM's multi-loci, the shipped ``LocusPriors``, and the ruler's inputs and output."""
    import rigel.calibration.capture_eff_length as CEL
    import rigel.calibration.fl as FL
    import rigel.calibration.priors as PRIORS
    from rigel.pipeline import run_pipeline

    cap: dict = {}
    orig_ruler = CEL.transcript_capture_eff_lengths
    orig_priors = PRIORS.assemble_priors
    orig_fl = FL.build_fl_models

    def fl_hook(*a, **k):
        models = orig_fl(*a, **k)
        cap["fl_models"] = models  # the last fit is calibration's
        return models

    def ruler_hook(calibration, region_arrays, index_, fl_eff_lengths, tables, junction_cap):
        out = orig_ruler(calibration, region_arrays, index_, fl_eff_lengths, tables, junction_cap)
        cap["fl"] = np.asarray(fl_eff_lengths, dtype=np.float64).copy()
        cap["rna_pmf"] = np.asarray(tables.rna, dtype=np.float64).copy()
        cap["length_tables"] = tables
        cap["ruler"] = np.asarray(out, dtype=np.float64).copy()
        return out

    def priors_hook(calibration, region_arrays, multi_loci):
        pri = orig_priors(calibration, region_arrays, multi_loci)
        cap.update(calibration=calibration, region_arrays=region_arrays, multi_loci=multi_loci, priors=pri)
        raise SystemExit(0)  # everything the read-out needs exists now; the EM is not run

    CEL.transcript_capture_eff_lengths = ruler_hook
    PRIORS.assemble_priors = priors_hook
    FL.build_fl_models = fl_hook
    try:
        run_pipeline(str(panel_dir / condition / ORACLE_BAM), index, config)
    except SystemExit:
        pass
    finally:
        CEL.transcript_capture_eff_lengths = orig_ruler
        PRIORS.assemble_priors = orig_priors
        FL.build_fl_models = orig_fl
    if "priors" not in cap or "ruler" not in cap:
        raise SystemExit("⛔ the pipeline never reached the prior assembly (TRAPS: an-ablation-that-never-ran)")
    return cap


def gdna_space_yield(sim: Sim, blocks: np.ndarray, n_loci: int, id2name: dict, width_step: int = 4) -> np.ndarray:
    """Per locus, ``Y_g = Σ_w f_g(w) Σ_{starts s whose fragment overlaps a block} (off + bind · best single-part
    overlap)``, read off the simulator's own gDNA-space landscape, the computation the reads were drawn with. A block
    on a reference the simulator never draws gDNA from (``sim.gdna_refs``) yields nothing, however it is probed.

    ``blocks`` is ``int64[n, 4]`` of ``(locus, ref_id, start, end)``. The width grid is thinned by ``width_step``,
    each kept width carrying its bin's pmf mass (the landscape is per width per reference)."""
    Y = np.zeros(n_loci, dtype=np.float64)
    if blocks.size == 0:
        return Y
    keep = np.arange(0, sim.gdna_widths.size, max(int(width_step), 1))
    wk = sim.gdna_widths[keep]
    pk = np.add.reduceat(sim.gdna_pmf, keep)
    live = sim.sampler.enabled and sim.bind > 0.0
    for rid in np.unique(blocks[:, 1]):
        name = id2name[int(rid)]
        if name not in sim.gdna_refs:
            continue
        RL = int(sim.ref_lengths[name])
        bb = blocks[blocks[:, 1] == rid]
        for w, pw in zip(wk, pk, strict=True):
            w = int(w)
            lo = np.maximum(bb[:, 2] - w + 1, 0)
            hi = np.minimum(bb[:, 3] - 1, RL - w)
            n_starts = np.maximum(hi - lo + 1, 0).astype(np.float64)
            # less the starts the simulator never drew, whose fragment covers an unknown base
            blocked = sim.sampler.blocked_starts("gdna", name, RL - w + 1, w)
            if blocked is not None:
                n_starts -= np.maximum(_blocked_below(blocked, hi + 1) - _blocked_below(blocked, lo), 0)
            sw = np.zeros(bb.shape[0])
            if live:
                pos, wt = sim.sampler.binding_landscape("gdna", name, RL, w)
                if pos.size:
                    cs = np.concatenate(([0.0], np.cumsum(wt)))
                    a = np.searchsorted(pos, lo, side="left")
                    b = np.searchsorted(pos, hi, side="right")
                    sw = cs[b] - cs[a]
            np.add.at(Y, bb[:, 0], pw * (sim.off * n_starts + sim.bind * sw))
    return Y


def _blocked_below(blocked: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Per ``x``, how many starts below it lie in the sorted, disjoint ``[start, end)`` runs ``blocked``: every run
    that starts below ``x`` lies wholly below it except the last, which may hold ``x``."""
    x = np.asarray(x, dtype=np.int64)
    length = blocked[:, 1] - blocked[:, 0]
    before = np.concatenate(([0], np.cumsum(length)))
    j = np.searchsorted(blocked[:, 0], x, side="left")
    last = np.maximum(j - 1, 0)
    partial = np.where(j > 0, np.minimum(x - blocked[last, 0], length[last]), 0)
    return (before[last] + partial).astype(np.float64)


def footprint_gdna_length(region_arrays, block: tuple, gdna_fl_pmf) -> float:
    """The shipped gDNA component's length at efficiency 1 over one footprint ``(ref_id, start, end)``:
    ``priors.assemble_priors`` itself, on a result carrying calibration's own geometry for this partition
    (the contained supports and ``effective_length.observable_cut_shares``) and no landscape."""
    from rigel.calibration.effective_length import contained_eff_length, observable_cut_shares
    from rigel.calibration.priors import assemble_priors
    from rigel.calibration.region_arrays import boundary_region_indices
    from rigel.calibration.result import CalibrationResult
    from rigel.locus import Locus, MultiLocus

    rid, s0, e0 = (int(v) for v in block)
    n = int(region_arrays.n_regions)
    ne = int(boundary_region_indices(np.asarray(region_arrays.ref_id))[0].shape[0])
    length = np.asarray(region_arrays.end, dtype=np.float64) - np.asarray(region_arrays.start, dtype=np.float64)
    z, ez = np.zeros(n), np.zeros(ne)
    cal = CalibrationResult(
        count_gdna_region=z, count_rna_region=z, count_gdna_boundary=ez, count_rna_boundary=ez,
        count_rna_spliced_boundary=ez, boundary_mass_per_crossing=np.ones(ne), count_rna_sj=np.zeros(0),
        boundary_spliced_mass_per_crossing=np.ones(ne), sj_mass_per_crossing=np.zeros(0),
        sj_src_region=np.zeros(0, dtype=np.int64), sj_dst_region=np.zeros(0, dtype=np.int64),
        gdna_region_eff_len=contained_eff_length(length, gdna_fl_pmf), gdna_boundary_eff_len=ez,
        gdna_boundary_conserved_len=observable_cut_shares(region_arrays, gdna_fl_pmf),
        gdna_frac_region=z, rna_pos_frac_region=z,
        rna_neg_frac_region=z, gdna_frac_boundary=ez, rna_pos_frac_boundary=ez, rna_neg_frac_boundary=ez,
        gdna_density_global=0.0, gdna_bulk_density=None, gdna_capture_weight=1.0,
        gdna_capture_saturation=1.0,
        gdna_capture_efficiency_region=np.ones(n), gdna_capture_efficiency_boundary=np.ones(ne),
        rna_sense_frac=0.5,
        n_regions=n, n_boundaries=ne, n_sj=0,
    )
    ml = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        span_bp=e0 - s0,
        loci=(Locus(ref=str(rid), ref_id=rid, start=s0, end=e0),),
    )
    return float(assemble_priors(cal, region_arrays, [ml]).gdna_eff_len[0])


def footprint_probed_fraction(sim: Sim, blocks: np.ndarray, n_loci: int, id2name: dict) -> np.ndarray:
    """Per locus: the fraction of the footprint's bases under a probe part, read off the sampler's gDNA-space probe
    map (``probe_intervals("gdna", ref)``, the union of its parts). It is the locus gDNA component's analogue of a
    transcript's probed fraction, so locus rows can be read by probed class beside transcript rows."""
    covered = np.zeros(n_loci, dtype=np.float64)
    length = np.zeros(n_loci, dtype=np.float64)
    if blocks.size == 0:
        return np.zeros(n_loci)
    for rid in np.unique(blocks[:, 1]):
        name = id2name[int(rid)]
        parts = sim.sampler.probe_intervals("gdna", name) if sim.sampler.enabled else []
        merged: list[list[int]] = []
        for iv in sorted(parts, key=lambda x: (int(x.start), int(x.end))):
            a, b = int(iv.start), int(iv.end)
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        ms = np.array([m[0] for m in merged], dtype=np.int64)
        me = np.array([m[1] for m in merged], dtype=np.int64)
        for key, _, s0, e0 in blocks[blocks[:, 1] == rid].tolist():
            length[key] += max(e0 - s0, 0)
            if ms.size == 0 or e0 <= s0:
                continue
            lo = int(np.searchsorted(me, s0, side="right"))
            hi = int(np.searchsorted(ms, e0, side="left"))
            if hi > lo:
                covered[key] += float(np.maximum(0, np.minimum(me[lo:hi], e0) - np.maximum(ms[lo:hi], s0)).sum())
    return np.where(length > 0, covered / np.maximum(length, 1.0), 0.0)


def junction_probed_flags(sim: Sim, index: TranscriptIndex) -> np.ndarray:
    """Per transcript: does a probe part in the transcript's own space strictly span one of its junctions?

    Read off the sampler's probe map (``probe_intervals("mrna", row)``) against the junction positions in transcript
    coordinates, 5′ to 3′. These are the junction-probed transcripts, whose junction placements bind a probe no
    gDNA witness can see."""
    from rigel.types import Strand

    t = index.t_df
    n = len(t)
    flags = np.zeros(n, dtype=bool)
    if not sim.sampler.enabled:
        return flags
    off, ex_s, ex_e, _ = index.build_exon_csr()
    off = np.asarray(off)
    ex_s = np.asarray(ex_s, dtype=np.int64)
    ex_e = np.asarray(ex_e, dtype=np.int64)
    strand = t["strand"].to_numpy().astype(np.int64)
    syn = t["is_synthetic"].to_numpy(dtype=bool)
    for i in range(n):
        if syn[i] or off[i + 1] - off[i] < 2:
            continue
        s, e = ex_s[off[i] : off[i + 1]], ex_e[off[i] : off[i + 1]]
        order = np.argsort(s)
        lengths = (e - s)[order]
        if int(strand[i]) != int(Strand.POS):
            lengths = lengths[::-1]
        junctions = np.cumsum(lengths)[:-1]
        intervals = sim.sampler.probe_intervals("mrna", i)
        if not intervals:
            continue
        for iv in intervals:
            if np.any((iv.start < junctions) & (junctions < iv.end)):
                flags[i] = True
                break
    return flags


def scale_stats(L: np.ndarray, Y: np.ndarray, weight: np.ndarray) -> dict:
    """``L / Y`` over the live members (both positive, weight positive): the weighted geometric mean, the weighted
    sd of the log, the median and the 10 % / 90 % quantiles, and the count."""
    L = np.asarray(L, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    w = np.asarray(weight, dtype=np.float64)
    live = np.isfinite(L) & np.isfinite(Y) & (L > 0) & (Y > 0) & (w > 0)
    if not live.any():
        return dict(n=0, geomean=float("nan"), sd_log=float("nan"), med=float("nan"), q10=float("nan"), q90=float("nan"))
    r = L[live] / Y[live]
    lx = np.log(r)
    ww = w[live]
    m = float(np.sum(ww * lx) / ww.sum())
    sd = float(np.sqrt(np.sum(ww * (lx - m) ** 2) / ww.sum()))
    q = np.quantile(r, [0.1, 0.5, 0.9])
    return dict(n=int(live.sum()), geomean=float(np.exp(m)), sd_log=sd, med=float(q[1]), q10=float(q[0]), q90=float(q[2]))


def scale_verdict(means: dict) -> tuple[float, str]:
    """The largest raw class-mean disagreement, as a fraction, and its pair.

    Reading it as a likelihood defect requires the component numerators to share a unit: rescaling both K and L
    for one component leaves K/L unchanged.
    """
    names = [k for k, v in means.items() if np.isfinite(v) and v > 0]
    worst, pair = 0.0, ""
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            d = abs(float(np.log(means[a] / means[b])))
            if d > worst:
                worst, pair = d, f"{a} / {b}"
    return float(np.expm1(worst)), pair


def locus_blocks(loci) -> tuple[int, np.ndarray]:
    """``(n_loci, blocks)``: the EM's loci as ``int64[n, 4]`` rows ``(locus, ref_id, start, end)``, one per
    non-empty footprint block."""
    n_loci = max((int(ml.multi_locus_id) for ml in loci), default=-1) + 1
    blocks = np.array(
        [(int(ml.multi_locus_id), int(b.ref_id), int(b.start), int(b.end)) for ml in loci for b in ml.loci if b.end > b.start],
        dtype=np.int64,
    ).reshape(-1, 4)
    return n_loci, blocks


def within_gene_sd(t_df: pd.DataFrame, live: np.ndarray, L: np.ndarray, Y: np.ndarray, weight: np.ndarray) -> float:
    """The weighted sd of ``log(L / Y)`` about each gene's own mean, over the ``live`` transcripts; NaN when none
    is."""
    if not live.any():
        return float("nan")
    g = t_df["g_index"].to_numpy()
    d = pd.DataFrame(dict(g=g[live], x=np.log(L[live] / Y[live]), w=weight[live]))
    d["x"] -= d.groupby("g")["x"].transform("mean")
    return float(np.sqrt(np.average(d["x"] ** 2, weights=d["w"])))


def scale_readout(index, index_dir, panel_dir, condition, config, min_frags: int, width_step: int):
    """The one-scale read-out for one condition: every hypothesis class's L / Y on the shipped lengths, the pairwise
    class-mean ratios and their verdict, the same by probed class, and the within-gene spread of the annotated
    transcripts."""
    t0 = time.time()
    sim = load_sim(index, index_dir, panel_dir, condition)
    truth = load_truth(index, index_dir, panel_dir, condition, sim)
    Y_t = np.where(np.isfinite(truth.factor), truth.L_plain * truth.factor, np.nan)
    cap = em_loci(index, panel_dir, condition, config)
    t = index.t_df
    syn = t["is_synthetic"].to_numpy(dtype=bool)
    single = (~syn) & t["is_single_exon"].to_numpy(dtype=bool)
    multi = (~syn) & ~single
    ruler, frags, pf = cap["ruler"], truth.frags, truth.probed_frac
    # the locus gDNA component against the gDNA-space yield, weighted by calibration's gDNA fragments
    id2name = {v: k for k, v in index.ref_name_to_id.items()}
    n_loci, blocks = locus_blocks(cap["multi_loci"])
    Y_g = gdna_space_yield(sim, blocks, n_loci, id2name, width_step)
    L_g = np.asarray(cap["priors"].gdna_eff_len, dtype=np.float64)
    gdna_count = np.asarray(cap["priors"].gdna_count, dtype=np.float64)
    jp = junction_probed_flags(sim, index)
    pf_loc = footprint_probed_fraction(sim, blocks, n_loci, id2name)

    def transcripts(m):
        return scale_stats(ruler[m], Y_t[m], frags[m])

    rows: list[tuple[str, dict]] = [
        ("locus gDNA component", scale_stats(L_g, Y_g, gdna_count)),
        ("synthetic span", transcripts(syn)),
        ("annotated single-exon", transcripts(single)),
        ("annotated multi-exon", transcripts(multi)),
    ]
    by_class: dict[str, dict] = {}
    for cname, lo, hi in CLASSES:
        m = multi & (pf >= lo) & (pf < hi)
        ml = (pf_loc >= lo) & (pf_loc < hi)
        spans = transcripts(syn & (pf >= lo) & (pf < hi))
        gdna = scale_stats(L_g[ml], Y_g[ml], gdna_count[ml])
        rows += [
            (f"  multi-exon {cname}, no junction probe", transcripts(m & ~jp)),
            (f"  multi-exon {cname}, JUNCTION-PROBED", transcripts(m & jp)),
            (f"  synthetic span {cname}", spans),
            (f"  locus gDNA footprint {cname}", gdna),
        ]
        by_class[cname] = {"isoforms": transcripts(m), "spans": spans, "gDNA": gdna}
    stats = dict(rows)
    means = {
        "gDNA": stats["locus gDNA component"]["geomean"],
        "spans": stats["synthetic span"]["geomean"],
        "isoforms": stats["annotated multi-exon"]["geomean"],
    }
    live = multi & np.isfinite(Y_t) & (Y_t > 0) & (ruler > 0) & (frags >= min_frags)
    return dict(condition=condition, stats=stats, means=means, by_class=by_class,
                within_gene_sd=within_gene_sd(t, live, ruler, Y_t, frags),
                n_junction_probed=int(jp[multi].sum()),
                n_multi=int(multi.sum()), n_loci=int(n_loci), elapsed=time.time() - t0,
                per_transcript=dict(t_id=t["t_id"].to_numpy(), kind=truth.kind, probed_frac=pf, junction_probed=jp, frags=frags,
                                    Y=Y_t, L=ruler),
                per_locus=dict(locus=np.arange(n_loci), Y=Y_g, L=L_g, gdna_count=gdna_count, probed_frac=pf_loc))


def _ratios(means: dict) -> str:
    return "  ".join(
        f"{a}/{b} {means[a] / means[b]:.3f}"
        for a, b in (("spans", "isoforms"), ("gDNA", "isoforms"), ("gDNA", "spans"))
        if np.isfinite(means[a]) and np.isfinite(means[b]) and means[b] > 0
    )


def print_scale(r: dict) -> None:
    print(
        f"⭐ {r['condition']}  ({r['elapsed']:.0f} s) — L / Y per hypothesis class, UNANCHORED (×10⁻³): weighted geometric "
        f"mean, sd of the log, median, 10 % / 90 %; {r['n_loci']:,} loci, {r['n_multi']:,} multi-exon transcripts of which "
        f"{r['n_junction_probed']:,} carry a probe part across one of their own junctions"
    )
    print(f"   {'class':<52}{'n':>7} {'geomean':>9} {'sd log':>7} {'median':>9} {'10 %':>9} {'90 %':>9}")
    for name, s in r["stats"].items():
        if s["n"] == 0:
            print(f"   {name:<52}{0:>7}")
            continue
        print(
            f"   {name:<52}{s['n']:>7} {s['geomean'] * 1e3:>9.4f} {s['sd_log']:>7.3f} {s['med'] * 1e3:>9.4f} "
            f"{s['q10'] * 1e3:>9.4f} {s['q90'] * 1e3:>9.4f}"
        )
    worst, pair = scale_verdict(r["means"])
    verdict = "ONE SCALE" if worst <= ONE_SCALE_TOL else f"NOT one scale ({worst:.1%} between {pair})"
    print(f"   class means: {_ratios(r['means'])}  →  {verdict}")
    print("   Cross-population scale requires matching likelihood-numerator units; compare K/L, not L alone.")
    print(f"   within-gene sd of log(L / Y), annotated multi-exon: {r['within_gene_sd']:.3f}   (tolerance {ONE_SCALE_TOL:.0%} on class means)")
    print("   WITHIN each probed class: class-mean L / Y ×10⁻³ and the pairwise ratios")
    print(f"   {'probed class':<16}{'n iso/span/loci':>17} {'isoforms':>9} {'spans':>9} {'gDNA':>9}  {'spans/iso':>9} {'gDNA/iso':>9} {'gDNA/spans':>10}")
    for cname, d in r["by_class"].items():
        iso, sp, g = (d[k] for k in ("isoforms", "spans", "gDNA"))
        m = [x["geomean"] if x["n"] else float("nan") for x in (iso, sp, g)]

        def rat(a, b):
            return f"{a / b:9.3f}" if np.isfinite(a) and np.isfinite(b) and b > 0 else f"{'…':>9}"

        print(f"   {cname:<16}{iso['n']:>5}/{sp['n']:>5}/{g['n']:>5} {m[0] * 1e3:>9.4f} {m[1] * 1e3:>9.4f} {m[2] * 1e3:>9.4f}  {rat(m[1], m[0])} {rat(m[2], m[0])} {rat(m[2], m[1]):>10}")


# ── the self-test ────────────────────────────────────────────────────────────────────────────────


@dataclasses.dataclass
class _Panel:
    """The self-test's synthetic panel: a tiny index with a fully probed (``p``), a half-probed (``h``) and an
    unprobed (``u``) gene, the probe file its sampler reads, and the fragment-length law."""

    tmp: Path
    index: TranscriptIndex
    transcripts: list
    sampler: CaptureSampler
    widths: np.ndarray
    pmf: np.ndarray
    off: float
    binding: float
    L: np.ndarray
    p: int
    h: int
    u: int

    def sim(self) -> Sim:
        """The panel as the yield reads it: RNA and gDNA under the same law."""
        widths = np.asarray(self.widths, dtype=np.int64)
        pmf = np.asarray(self.pmf, dtype=np.float64)
        return Sim(self.sampler, widths, pmf, self.off, self.binding, {"chr1": 16000}, self.transcripts, widths, pmf,
                   frozenset({"chr1"}))


#: the self-test genome's repeating unit: 20 bases holding ``20 * gc_optimum`` G or C, so every window a multiple of
#: 20 long, every probe of the self-test panel among them, holds exactly the optimum's GC and binds at strength 1
_UNIT = 20


def _self_test_genome(path: Path, length: int, at_strength: dict | None = None) -> Path:
    """A genome of ``length`` bases repeating the optimum's unit, with each ``(start, end)`` of ``at_strength``
    written as A and T instead (GC 0, a probe there binds weakly); indexed for ``pysam``."""
    n_gc = round(_UNIT * CaptureConfig.gc_optimum)
    unit = ("GC" * _UNIT)[:n_gc] + ("AT" * _UNIT)[: _UNIT - n_gc]
    bases = list(unit * (length // _UNIT))
    for lo, hi in at_strength or ():
        bases[lo:hi] = ("AT" * (hi - lo))[: hi - lo]
    path.write_text(">chr1\n" + "".join(bases) + "\n")
    pysam.faidx(str(path))
    return path


def _self_test_panel(tmp: Path) -> _Panel:
    # the self-test is test code, so it builds its index with the suite's one synthetic-index builder
    sys.path.insert(0, str(_REPO / "tests"))
    from _index_builder import build_index

    gtf = "".join(
        f'chr1\ttest\texon\t{s + 1}\t{s + 600}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n'
        for g, t, starts in (
            ("gP", "probed", (1000, 3000)),
            ("gH", "half", (6000, 8000)),
            ("gU", "unprobed", (11000, 13000)),
        )
        for s in starts
    )

    (tmp / "ruler_self_test").mkdir()
    index = build_index(tmp / "ruler_self_test", gtf, {"chr1": 16000})
    bed = tmp / "probes.bed"
    # the probed gene's one probe spans both exons as a two-block BED12 line, contiguous in cDNA, so every
    # start of every fragment lies under it at full scale; the half-probed gene's probe covers its first exon
    bed.write_text(
        "chr1\t1000\t3600\tp\t0\t.\t1000\t3600\t0\t2\t600,600\t0,2000\n"
        "chr1\t6000\t6600\th\t0\t.\t6000\t6600\t0\t1\t600\t0\n"
    )
    transcripts = load_transcripts_from_index(index.index_dir)
    off, binding = 1.0, 10.0
    # every probe binds at full strength on this genome, so the closed forms below read binding_per_base itself
    with pysam.FastaFile(str(_self_test_genome(tmp / "optimum.fa", 16000))) as fasta:
        sampler = CaptureSampler.from_config(
            CaptureConfig(
                probes=str(bed),
                probe_format="bed12",
                off_target_weight=off,
                binding_per_base=binding,
                min_overlap=1,
            ),
            transcripts,
            {"chr1": 16000},
            fasta=fasta,
        )
    widths, pmf = fl_pmf(type("Sim", (), dict(frag_mean=200, frag_std=50, frag_min=100, frag_max=300)))
    row = {tid: int(i) for i, tid in enumerate(index.t_df.t_id)}
    L = np.array([int(t.length) for t in transcripts], dtype=np.int64)
    return _Panel(tmp, index, transcripts, sampler, widths, pmf, off, binding, L,
                  row["probed"], row["half"], row["unprobed"])


def _check_manifest_loader(panel: _Panel, check) -> None:
    """The manifest loader keeps the RNA and gDNA laws apart and masks what the simulator never drew."""
    # its genome is known sequence but for two assembly gaps, one in the unprobed gene's intron and one in its
    # second exon
    bases = np.random.default_rng(0).choice(np.frombuffer(b"ACGT", dtype=np.uint8), 16000)
    bases[12200:12300] = ord("N")
    bases[13100:13120] = ord("N")
    genome = panel.tmp / "known.fa"
    genome.write_text(">chr1\n" + bases.tobytes().decode() + "\n")
    pysam.faidx(str(genome))
    rna_law = dict(frag_mean=200, frag_std=50, frag_min=100, frag_max=300)
    manifest = {
        "genome": str(genome),
        "simulation": rna_law,
        "gdna": dict(frag_mean=80, frag_std=10, frag_min=60, frag_max=100, genomic_refs=["chr1"]),
        "capture_configs": [{"label": "off", "config": {"off_target_weight": 1.0}}],
    }
    manifest_path = panel.tmp / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))

    def load():
        return load_sim(panel.index, Path(panel.index.index_dir), panel.tmp, "example_capture_off")

    loaded = load()
    gw, gp = fl_pmf(type("GDNA", (), {k: v for k, v in manifest["gdna"].items() if k != "genomic_refs"}))
    check("the manifest loader keeps distinct RNA and gDNA laws",
          np.array_equal(loaded.widths, panel.widths) and np.array_equal(loaded.pmf, panel.pmf)
          and np.array_equal(loaded.gdna_widths, gw) and np.array_equal(loaded.gdna_pmf, gp))
    # the loader masks the gDNA reference at its gaps, and the unprobed transcript at the exon gap in its own
    # coordinates (600 bases of first exon, then 100 into the second)
    u = panel.u
    check("the loader masks the gDNA reference at the genome's unknown runs",
          loaded.sampler.unknown("gdna", "chr1").tolist() == [[12200, 12300], [13100, 13120]])
    check("...and each transcript at its own", loaded.sampler.unknown("mrna", u).tolist() == [[700, 720]])
    plain, factor, _ = sampler_truth(loaded.sampler, panel.L, panel.widths, panel.pmf, 1.0)
    enumerated = sum(
        pw * sum(1 for st in range(int(panel.L[u]) - int(wv) + 1) if st + int(wv) <= 700 or st >= 720)
        for wv, pw in zip(panel.widths, panel.pmf, strict=True)
    )
    check("a masked transcript's plain length counts only the starts the simulator could draw",
          abs(plain[u] / enumerated - 1.0) < 1e-12)
    check("...so masking leaves an uncaptured transcript's factor exactly 1", factor[u] == 1.0)
    # the gDNA-space yield is the simulator's: on a reference the manifest draws gDNA from, a probed footprint yields
    # a captured gDNA mass; on one it never draws gDNA from (an RNA-only spike-in backbone) it yields none, probed or
    # not, so no locus there can be read against a gDNA yield
    capture_on = {**dataclasses.asdict(CaptureConfig()), "probes": str(panel.tmp / "probes.bed"),
                  "probe_format": "bed12", "off_target_weight": panel.off, "binding_per_base": panel.binding,
                  "min_overlap": 1}
    blocks = np.array([(0, 0, 11000, 13600), (1, 0, 1000, 3600)], dtype=np.int64)
    for refs in (["chr1"], []):
        manifest_path.write_text(json.dumps({
            **manifest, "gdna": {**manifest["gdna"], "genomic_refs": refs},
            "capture_configs": [{"label": "on", "config": capture_on}],
        }))
        captured = load_sim(panel.index, Path(panel.index.index_dir), panel.tmp, "example_capture_on")
        Y = gdna_space_yield(captured, blocks, 2, {0: "chr1"}, width_step=4)
        if refs:
            check("a probed footprint on a reference gDNA is drawn from yields a captured gDNA mass",
                  Y[0] > 0.0 and Y[1] > 2.0 * Y[0])
        else:
            check("...and on a reference gDNA is never drawn from it yields none, probed or not", not Y.any())
    # a captured record lacking a field of the simulator's capture config was simulated under another weight model
    manifest_path.write_text(json.dumps({
        **manifest, "capture_configs": [{"label": "on", "config": {
            k: v for k, v in capture_on.items() if k != "gc_half_distance"}}],
    }))
    check("a captured record lacking a capture-config field is refused",
          refuses(load_sim, panel.index, Path(panel.index.index_dir), panel.tmp, "example_capture_on"))
    manifest_path.write_text(json.dumps({
        **manifest, "capture_configs": [{"label": "on", "config": {**capture_on, "gdna_split_penalty": 0.2}}],
    }))
    check("...and so is one carrying a field this simulator does not know",
          refuses(load_sim, panel.index, Path(panel.index.index_dir), panel.tmp, "example_capture_on"))
    manifest_path.write_text(json.dumps({k: v for k, v in manifest.items() if k != "genome"}))
    fired = refuses(load)
    manifest_path.write_text(json.dumps(manifest))
    check("a manifest naming no genome is refused", fired)
    manifest["gdna"] = {"frag_mean": 100}
    manifest_path.write_text(json.dumps(manifest))
    loaded = load()
    gw, gp = fl_pmf(dataclasses.replace(GDNASimConfig(), frag_mean=100))
    check("missing gDNA fields use simulator defaults independently of RNA",
          np.array_equal(loaded.gdna_widths, gw) and np.array_equal(loaded.gdna_pmf, gp))


def _check_truth(panel: _Panel, check) -> Truth:
    """The truth factor and the probed fraction, read off the sampler; returns the panel's truth."""
    L_plain, factor, probed = sampler_truth(panel.sampler, panel.L, panel.widths, panel.pmf, panel.off)
    p, h, u = panel.p, panel.h, panel.u
    # an unprobed transcript's factor is exactly 1, a fully probed one's is 1 + binding·E[w]/off, a half-probed
    # one's lies strictly between
    check("an unprobed transcript's truth factor is exactly 1", factor[u] == 1.0)
    mean_w = float(np.dot(panel.widths, panel.pmf))
    check(
        "a fully probed transcript's truth factor is 1 + binding·E[w]/off",
        abs(factor[p] / (1.0 + panel.binding * mean_w / panel.off) - 1.0) < 0.02,
    )
    check("a half-probed transcript's factor lies between", 1.0 < factor[h] < factor[p])
    # the class is read off the sampler: probed 1, unprobed 0, half-probed ½
    check("probed fraction: fully probed reads 1", abs(probed[p] - 1.0) < 1e-9)
    check("probed fraction: unprobed reads 0", probed[u] == 0.0)
    check("probed fraction: half-probed reads ½", abs(probed[h] - 0.5) < 1e-9)
    n = panel.L.size
    return Truth(L_plain, factor, probed, transcript_kinds(panel.index.t_df), np.full(n, 100.0), np.full(n, 10.0))


def _check_strength(panel: _Panel, truth: Truth, check) -> None:
    """A probe binds at the strength its sequence gives: the same panel on a genome where the probed gene's probe is
    A and T alone captures that gene less, while its probed fraction, geometry, stays 1."""
    weak_genome = _self_test_genome(panel.tmp / "weak.fa", 16000, at_strength=[(1000, 1600), (3000, 3600)])
    with pysam.FastaFile(str(weak_genome)) as fasta:
        weak = CaptureSampler.from_config(panel.sampler.config, panel.transcripts, {"chr1": 16000}, fasta=fasta)
    _, factor, probed = sampler_truth(weak, panel.L, panel.widths, panel.pmf, panel.off)
    p, h = panel.p, panel.h
    check("an AT-only probe binds weaker: its gene's truth factor falls", 1.0 < factor[p] < truth.factor[p])
    check("...the other probe's gene is untouched", factor[h] == truth.factor[h])
    check("...and its probed fraction, geometry, stays 1", abs(probed[p] - 1.0) < 1e-9)


def _check_score(panel: _Panel, truth: Truth, check) -> None:
    """The anchored log error: zero on the truth, blind to a global scale, anchored on the probed class."""
    p, u = panel.p, panel.u
    factor = truth.factor
    d, measured, n_zero = score(factor, truth, min_frags=20)
    rows = class_rows(d, measured, truth)
    check("the truth scores exactly 0 on every class", all(abs(r["med"]) < 1e-12 for r in rows))
    check("every class is 100 % within the tolerance", all(r["within"] == 1.0 for r in rows))
    check("no zero factor on the truth", n_zero == 0)
    d3, _, _ = score(factor * 3.0, truth, min_frags=20)
    check("a global scale is anchored away", np.allclose(d3[measured], d[measured]))
    # the unprobed class moves by log 2 when its factor is doubled, the probed class not at all
    f2 = factor.copy()
    f2[u] *= 2.0
    d2, m2, _ = score(f2, truth, min_frags=20)
    check("a doubled unprobed factor reads +log 2 on the unprobed class", abs(d2[u] - np.log(2.0)) < 1e-9)
    check("...and leaves the probed class at 0", abs(d2[p]) < 1e-12)
    check("the headline carries the unprobed median", abs(headline(d2, m2, truth)["unprobed_med"] - np.log(2.0)) < 1e-9)
    # a factor of exactly 0 is counted and dropped, never a −inf in the quantiles
    f0 = factor.copy()
    f0[u] = 0.0
    d0, m0, z0 = score(f0, truth, min_frags=20)
    check("a zero factor is counted", z0 == 1)
    check("...and dropped from the measured set", not m0[u] and np.isfinite(d0[m0]).all())
    # the fragment floor: a transcript under it is not scored
    few = dataclasses.replace(truth, frags=np.where(np.arange(factor.size) == u, 5.0, 100.0))
    _, m5, _ = score(factor, few, min_frags=20)
    check("a transcript under the fragment floor is not measured", not m5[u] and m5[p])
    # the anchor is the probed class: with the probed class biased, everything shifts by its bias
    fb = factor.copy()
    fb[p] *= 5.0
    db, _, _ = score(fb, truth, min_frags=20)
    check("the anchor is the probed class's median", abs(db[p]) < 1e-12 and abs(db[u] + np.log(5.0)) < 1e-9)


def _check_module_loader(tmp: Path, check) -> None:
    """The prototype loader refuses a file without ``ARMS`` and a reserved arm name."""
    bad = tmp / "bad.py"
    bad.write_text("X = 1\n")
    check("a module without ARMS is refused", refuses(load_arms, bad))
    reserved = tmp / "reserved.py"
    reserved.write_text("ARMS = {'shipped': lambda *a: None}\n")
    check("a prototype arm named `shipped` is refused", refuses(load_arms, reserved))
    good = tmp / "good.py"
    good.write_text("def r(cal, ra, index, fl, **inputs):\n    return fl\nARMS = {'plain': r}\n")
    check("a well-formed prototype loads", list(load_arms(good)) == ["plain"])


def _check_one_scale(panel: _Panel, truth: Truth, check) -> None:
    """A class priced at its own yield reads 1; a class scaled by 1.185 is named as the pair, and a per-mille
    difference is inside the tolerance."""
    Y = truth.L_plain * truth.factor
    s_true = scale_stats(Y, Y, np.full(panel.L.size, 100.0))
    check("a class priced at its own yield reads L / Y = 1 with no spread", abs(s_true["geomean"] - 1.0) < 1e-12 and s_true["sd_log"] < 1e-12)
    worst, pair = scale_verdict({"gDNA": 1.0, "spans": 1.185, "isoforms": 1.0})
    check("a 1.185 scale on one class is named as its pair", abs(worst - 0.185) < 1e-9 and "spans" in pair)
    check("...and is outside the tolerance, while a per-mille difference is inside", worst > ONE_SCALE_TOL and scale_verdict({"a": 1.0, "b": 1.001})[0] < ONE_SCALE_TOL)


def _enumerated_yield(sampler: CaptureSampler, sim: Sim, blocks: np.ndarray) -> np.ndarray:
    """The gDNA-space yield by direct enumeration of every start whose fragment overlaps each block."""
    direct = np.zeros(2)
    for lid, _rid, lo, hi in blocks:
        for wv, pw in zip(sim.gdna_widths, sim.gdna_pmf, strict=True):
            direct[lid] += pw * sum(
                sampler.fragment_weight("gdna", "chr1", 16000, start, int(wv))
                for start in range(max(0, int(lo-wv+1)), min(int(hi), 16000-int(wv)+1))
            )
    return direct


def _check_gdna_space(panel: _Panel, check) -> None:
    """The junction-probed flag, the gDNA-space yield (under unequal laws and under a mask), and the footprint's
    length and probed fraction, each read off the sampler's own probe map."""
    sim, p, h, u, off = panel.sim(), panel.p, panel.h, panel.u, panel.off
    # the probed gene's two-block probe spans its junction, the half-probed gene's single-exon probe does not, the
    # unprobed gene has none
    jp = junction_probed_flags(sim, panel.index)
    check("a probe spanning a transcript's junction flags it, and only it", bool(jp[p]) and not jp[h] and not jp[u])
    # the yield counts every start whose fragment overlaps the footprint, so a block over the unprobed gene reads
    # exactly off · Σ_w f(w)(L + w − 1); a block over the probed gene reads more
    blocks = np.array([(0, 0, 11000, 13600), (1, 0, 1000, 3600)], dtype=np.int64)
    Yg = gdna_space_yield(sim, blocks, 2, {0: "chr1"}, width_step=1)
    plain = float(sum(pw * off * (2600 + int(wv) - 1) for wv, pw in zip(panel.widths, panel.pmf, strict=True)))
    check("an unprobed footprint's gDNA-space yield is its off-target plain length", abs(Yg[0] / plain - 1.0) < 1e-9)
    check("...and a probed footprint's yield exceeds it", Yg[1] > Yg[0] * 2.0)
    # unequal generating laws: the gDNA truth follows its own pmf, independent of RNA
    gap = dataclasses.replace(sim, gdna_widths=np.array([60, 80, 100]), gdna_pmf=np.array([0.25, 0.5, 0.25]))
    Ygap = gdna_space_yield(gap, blocks, 2, {0: "chr1"}, width_step=1)
    check("a gDNA length gap changes its footprint truth", not np.array_equal(Ygap, Yg))
    check("the unprobed gap footprint uses E_g[w]", abs(Ygap[0] - off * (2600 + 80 - 1)) < 1e-9)
    changed_rna = dataclasses.replace(gap, widths=np.array([400]), pmf=np.array([1.0]))
    check("changing only RNA widths leaves gDNA truth identical",
          np.array_equal(Ygap, gdna_space_yield(changed_rna, blocks, 2, {0: "chr1"}, width_step=1)))
    check("gap gDNA truth agrees with direct captured-start enumeration",
          np.allclose(Ygap, _enumerated_yield(panel.sampler, gap, blocks), rtol=1e-12, atol=1e-9))
    _check_footprints(panel, sim, plain, check)
    # under a mask: the yield drops every start whose fragment covers an unknown base, as the simulator's own
    # weights do, and a probed transcript's probed fraction is over its known bases
    sampler = panel.sampler
    sampler.mask_unknown("gdna", "chr1", [(12000, 12100), (2000, 2050)])
    Ym = gdna_space_yield(gap, blocks, 2, {0: "chr1"}, width_step=1)
    check("a masked gDNA footprint's yield agrees with direct enumeration of the drawable starts",
          np.allclose(Ym, _enumerated_yield(sampler, gap, blocks), rtol=1e-12, atol=1e-9) and (Ym < Ygap).all())
    sampler.mask_unknown("gdna", "chr1", [])
    sampler.mask_unknown("mrna", p, [(100, 140)])
    _, _, probed_m = sampler_truth(sampler, panel.L, panel.widths, panel.pmf, off)
    check("a fully probed transcript with unknown bases still reads probed fraction 1", abs(probed_m[p] - 1.0) < 1e-9)
    sampler.mask_unknown("mrna", p, [])


def _check_footprints(panel: _Panel, sim: Sim, plain: float, check) -> None:
    """The shipped gDNA component's length at efficiency 1 counts the yield's overlapping starts, and the
    footprint's probed fraction is read off the sampler's gDNA-space probe map."""
    # the footprint's outside neighbours hold every outside base, so each crossing fragment starts in a neighbour
    # and its whole unit lands on the footprint's objects: Σ_w f(w)(L + w − 1), as the yield counts
    ra_t = RegionArrays.from_index(panel.index)
    # the geometry indexes a pmf by width (calibration's convention); the simulator's pair is positional
    by_width = np.zeros(int(max(panel.widths)) + 1)
    by_width[np.asarray(panel.widths, dtype=np.int64)] = np.asarray(panel.pmf, dtype=np.float64)
    Lfp = footprint_gdna_length(ra_t, (0, 11000, 13600), by_width)
    check("the gDNA component's length at efficiency 1 is the yield's overlapping-start count", abs(Lfp / (plain / panel.off) - 1.0) < 1e-9)
    # the unprobed gene's block reads 0, a block exactly on one probe part reads 1, a block half on it reads ½
    parts = sim.sampler.probe_intervals("gdna", "chr1")
    a, b = int(parts[0].start), int(parts[0].end)
    pf_loc = footprint_probed_fraction(sim, np.array([(0, 0, 11000, 13600), (1, 0, a, b), (2, 0, a - (b - a), b)], dtype=np.int64), 3, {0: "chr1"})
    check("an unprobed footprint's probed fraction is 0", pf_loc[0] == 0.0)
    check("a footprint exactly on a probe part reads 1", abs(pf_loc[1] - 1.0) < 1e-12)
    check("...and one half on it reads ½", abs(pf_loc[2] - 0.5) < 1e-12)


def self_test() -> int:
    """Perturb every comparator on a synthetic panel with a known partition, and require each to fire.

    The panel is a tiny index with a fully probed, a half-probed and an unprobed gene and a probe file the sampler
    reads. The truth is the sampler's own; the arms are rulers written here."""
    import tempfile

    checks: list[tuple[str, bool]] = []

    def check(name: str, ok: bool):
        checks.append((name, bool(ok)))

    with tempfile.TemporaryDirectory() as tmp:
        panel = _self_test_panel(Path(tmp))
        _check_manifest_loader(panel, check)
        truth = _check_truth(panel, check)
        _check_strength(panel, truth, check)
        _check_score(panel, truth, check)
        _check_module_loader(panel.tmp, check)
        _check_one_scale(panel, truth, check)
        _check_gdna_space(panel, check)

    n = len(checks)
    failed = [name for name, ok in checks if not ok]
    for name, ok in checks:
        print(f"   {'✔' if ok else '⛔'} {name}")
    if failed:
        print(f"⛔ {len(failed)}/{n} self-test gates did NOT fire")
        return 1
    print(f"⭐ {n}/{n} self-test gates fired")
    return 0


if __name__ == "__main__":
    sys.exit(main())
