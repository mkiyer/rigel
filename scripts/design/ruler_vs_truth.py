#!/usr/bin/env python
"""IS THE RULER'S FORMULA RIGHT? — the EM's effective length under capture, scored per transcript against
the simulator's own capture-aware effective length.

The truth is what generated the reads: ``CaptureSampler.partition_array`` is the sampler's count of
admissible starts weighted by capture, so a transcript's true effective length under capture is
``Σ_w f_pre(w) · partition_t(w)`` over the pre-capture fragment-length pmf, and its plain length is
``Σ_w f_pre(w) · off_target_weight · (L_t − w + 1)+``. The truth FACTOR is their ratio; an arm's factor is
its effective length over the same plain lengths (``transcript_capture_eff_lengths(...) / fl``). The two
are compared in the log, ANCHORED on the fully probed transcripts (the reference makes them the unit on
both sides, so a global scale is free), and summarised per class — the probed fraction of the transcript's
bases, read off the sampler itself — and per kind (mRNA / annotated single-exon / synthetic nascent
entity). A factor of exactly 0 is an infinite abundance: counted, and dropped from the quantiles.

Arms: ``shipped`` is the ruler in ``src/`` on the shipped calibration; ``oracle_gdna`` is the same ruler
fed the CERTIFIED true gDNA counts per object (``slot_truth.npz``) with the reference read off a landscape
fitted on the truth — the ideal witness, so the gap between the two is the calibration's and what
remains under ``oracle_gdna`` is the ruler's (or the witness geometry's,
`ISSUES: ruler-witness-geometry-on-transcript-panels`); ``--module`` names a prototype file defining
``ARMS = {name: ruler}`` with ``ruler(calibration, region_arrays, index, fl, **inputs)`` returning the
effective lengths — ``inputs`` is ``calibrate``'s debug bundle (the last refit's landscape under
``gdna_hyperprior``, the chain, the belief) plus ``gdna_fl_pmf`` / ``rna_fl_pmf`` — run beside the
shipped one on the same calibration — the harness every ruler mechanism
is developed on before ``src/``; ``--set SECTION.FIELD=VALUE`` prices a config value on every arm.

⛔ Read the classes apart: the probed class is where the formula is exact when its witness is, the
unprobed class is where the reference bites, and the partial classes are where the junction price lives.
The junction price reads a junction's capture from the gDNA objects beside it, by conservation of bases,
right on average and noisy per junction (`ISSUES: the-junction-price-is-noisy-within-a-gene`); what no gDNA
object sees is declared, not repaired (`ISSUES: ruler-witness-geometry-on-transcript-panels`).

Usage::

    python scripts/design/ruler_vs_truth.py --panel test --condition gdna_g05_ss_0.99_nrna_file_capture_on
    python scripts/design/ruler_vs_truth.py --panel test                      # every capture-ON condition, one line each
    python scripts/design/ruler_vs_truth.py --panel test --panel-dir ~/Downloads/rigel_runs/test_reference/scenarios_depth_d10
    python scripts/design/ruler_vs_truth.py --panel ladder --condition C --module proto.py --out table.tsv
    python scripts/design/ruler_vs_truth.py --panel test --condition C --set calibration.message_policy=silent
    python scripts/design/ruler_vs_truth.py --panel ladder --condition C --scale        # THE ONE-SCALE READ-OUT
    python scripts/design/ruler_vs_truth.py --panel ladder --panel-dir ~/Downloads/rigel_runs/suite/ladder_nrna_lo --condition gdna_g50_ss_0.99_nrna_lo_capture_on --scale
    python scripts/design/ruler_vs_truth.py --self-test

⛔ ``--panel-dir`` names a directory of conditions, never the index: pass ``--panel ladder`` with it for a ladder-derived
panel (``ladder_nrna_lo``, the depth ladder), or the test chromosome's index is loaded against the ladder's BAM and no
fragment deposits.

``--scale`` — ARE EVERY HYPOTHESIS CLASS'S OPPORTUNITIES ON ONE SCALE? The E-step balances two hypotheses at a
fragment only if their opportunities are the same multiple of the yield the simulator actually drew them with, so
the quantity that decides the gDNA-versus-RNA split under capture is L / Y per hypothesis class, UNANCHORED: L the
length the EM divides by, Y the simulator's own expected captured yield per unit abundance — in mRNA space for a
transcript or a synthetic span (``partition_array``), in gDNA space for the locus gDNA component (the same weight
over the locus's footprint, ``_extra_landscape``). The read-out prints L / Y for the locus gDNA component
(``LocusPriors.gdna_eff_len``), the synthetic spans and the annotated transcripts on the shipped lengths, by probed
class with the JUNCTION-PROBED transcripts apart, then the pairwise class-mean ratios. ⛔ The tolerance is the
converged elasticity's: class means must agree to about 1 % (the synthetic pool moves ~20 % per 1 % of its own
opportunity at ``g50 ss.99 ON``). The pipeline is run to the prior assembly and stopped there — no EM.
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

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO / "scripts" / "design"))

from _shared import set_field, sibling  # noqa: E402
from rigel.calibration.abundance_landscape import located_enriched_mode  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.capture_eff_length import transcript_capture_eff_lengths  # noqa: E402
from rigel.calibration.capture_efficiency import capture_efficiencies  # noqa: E402
from rigel.calibration.landscape import fit_landscape  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from rigel.sim.capture import CaptureConfig, CaptureSampler  # noqa: E402
from rigel.sim.whole_genome import fl_pmf, load_transcripts_from_index  # noqa: E402

PB = sibling("policy_benchmark.py")
PANELS = PB.PANELS

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
#: the keys a capture config carries from the panel manifest
_CAPTURE_KEYS = (
    "probes",
    "probe_format",
    "off_target_weight",
    "binding_per_base",
    "min_overlap",
)


# ── the truth ────────────────────────────────────────────────────────────────────────────────────


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
    """``(L_plain, factor_true, probed_frac)`` from the sampler's own partition. The probed fraction is
    read off the partition at a one-base fragment: every start then carries exactly its own base, so
    ``(partition(1) − off·L) / binding`` is the count of probed bases (a split probe's pieces at their
    penalty), and no private probe map is consulted."""
    L = np.asarray(lengths, dtype=np.int64)
    n = L.size
    rows = range(n)
    L_true = np.zeros(n)
    L_plain = np.zeros(n)
    for w, pw in zip(widths, pmf, strict=True):
        L_true += pw * sampler.partition_array("mrna", rows, L, int(w))
        L_plain += pw * off_target * np.maximum(0.0, L - w + 1.0)
    ok = L_plain > 0
    factor = np.where(ok, L_true / np.maximum(L_plain, 1e-12), np.nan)
    binding = float(sampler.config.binding_per_base)
    if sampler.enabled and binding > 0.0:
        one = sampler.partition_array("mrna", rows, L, 1)
        probed = np.clip((one - off_target * L) / (binding * np.maximum(L, 1)), 0.0, 1.0)
    else:
        probed = np.zeros(n)
    return L_plain, factor, probed


@dataclasses.dataclass(frozen=True, slots=True)
class Sim:
    """The condition's simulator as the truth needs it: the sampler on the panel's probes, the pre-capture
    fragment-length law, the two weight constants, the reference lengths and the transcript rows."""

    sampler: CaptureSampler
    widths: np.ndarray
    pmf: np.ndarray
    off: float
    bind: float
    ref_lengths: dict
    transcripts: list


def load_sim(index: TranscriptIndex, index_dir: Path, panel_dir: Path, condition: str) -> Sim:
    manifest = json.load(open(panel_dir / "manifest.json"))
    label = "on" if condition.endswith("_capture_on") else "off"
    cap_cfg = next(c["config"] for c in manifest["capture_configs"] if c["label"] == label)
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
    sampler = CaptureSampler.from_config(
        CaptureConfig(**{k: cap_cfg[k] for k in _CAPTURE_KEYS if k in cap_cfg}),
        transcripts,
        ref_lengths,
    )
    params = type("Sim", (), {k: sim[k] for k in ("frag_mean", "frag_std", "frag_min", "frag_max")})
    widths, pmf = fl_pmf(params)
    return Sim(
        sampler,
        np.asarray(widths, dtype=np.int64),
        np.asarray(pmf, dtype=np.float64),
        float(cap_cfg["off_target_weight"]),
        float(cap_cfg.get("binding_per_base", 0.0)),
        ref_lengths,
        transcripts,
    )


def load_truth(index: TranscriptIndex, index_dir: Path, panel_dir: Path, condition: str) -> Truth:
    sim = load_sim(index, index_dir, panel_dir, condition)
    sampler, widths, pmf, transcripts = sim.sampler, sim.widths, sim.pmf, sim.transcripts
    L = np.array([int(t.length) for t in transcripts], dtype=np.int64)
    L_plain, factor, probed = sampler_truth(sampler, L, widths, pmf, sim.off)
    t = index.t_df
    kind = np.where(
        t.is_synthetic.to_numpy(bool),
        "entity",
        np.where(t.is_nrna.to_numpy(bool), "single-exon", "mRNA"),
    )
    truth = pd.read_csv(panel_dir / condition / "truth_abundances.tsv", sep="\t").set_index(
        "transcript_id"
    )
    frags = truth["observed_total_rna_fragments"].reindex(t.t_id).fillna(0).to_numpy(float)
    abund = truth["pre_capture_total_rna"].reindex(t.t_id).fillna(0).to_numpy(float)
    return Truth(L_plain, factor, probed, kind, frags, abund)


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def calibrate_condition(index, region_arrays, panel_dir: Path, condition: str, config):
    """The shipped calibration on the condition's cached payload; returns the result and the solve's
    debug bundle — the last refit's landscape under ``gdna_hyperprior`` (``None`` without a refit), the
    prior an expectation arm reads; the chain and the belief; and the two fragment-length pmfs."""
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    cache = read_scan_cache(panel_dir / "oracle_cache" / condition / "_main", index)
    kw = calibration_inputs(cache, index)
    debug: dict = {}
    cal = calibrate(
        payload=kw["payload"],
        config=config.calibration,
        region_arrays=region_arrays,
        strand_model=kw["strand_model"],
        gdna_fl_pmf=kw["gdna_fl_pmf"],
        rna_fl_pmf=kw["rna_fl_pmf"],
        sj=sj,
        boundary_flags=bflags,
        _debug=debug,
    )
    debug["gdna_fl_pmf"] = kw["gdna_fl_pmf"]
    debug["rna_fl_pmf"] = kw["rna_fl_pmf"]
    return cal, debug


def oracle_gdna_result(cal, region_arrays, panel_dir: Path, condition: str, gdna_fl_pmf):
    """The same result with the ruler's gDNA masses replaced by the CERTIFIED truth per object, the
    reference read off a landscape fitted on the truth and the efficiencies recomputed from both: the
    ideal witness, no estimator noise."""
    z = np.load(panel_dir / "oracle_cache" / condition / "slot_truth.npz", allow_pickle=True)
    kind, obj, n_g = z["kind"], z["obj"], z["n_gdna"]
    m_reg = np.zeros_like(np.asarray(cal.count_gdna_region, float))
    m_bnd = np.zeros_like(np.asarray(cal.count_gdna_boundary, float))
    r = kind == int(REGION)
    b = kind == int(BOUNDARY)
    m_reg[obj[r]] = n_g[r]
    m_bnd[obj[b]] = n_g[b]
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
    mode = located_enriched_mode(ls) if ls is not None else None
    ref = None if mode is None else float(np.exp(mode.mode.log_rho))
    if ref is None:
        efficiency = np.ones(m_reg.shape[0])
        efficiency_boundary = np.ones(m_bnd.shape[0])
    else:
        efficiency, efficiency_boundary = capture_efficiencies(
            ls,
            ref,
            m_reg,
            S,
            m_bnd,
            np.asarray(cal.gdna_boundary_eff_len, float),
        )
    return dataclasses.replace(
        cal,
        count_gdna_region=m_reg,
        count_gdna_boundary=m_bnd,
        gdna_reference_density=ref,
        gdna_reference_members=0 if mode is None else mode.n_members,
        gdna_capture_efficiency_region=efficiency,
        gdna_capture_efficiency_boundary=efficiency_boundary,
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


def _fmt_ref(cal) -> str:
    if cal.gdna_reference_density is None:
        return "None (factor 1 everywhere)"
    return f"{cal.gdna_reference_density:.3e}/bp from {cal.gdna_reference_members:,} located kernels"


def run_condition(index, region_arrays, index_dir, panel_dir, condition, config, arms, min_frags):
    """Every arm on one condition: ``(truth, {arm: (factor, calibration)})``."""
    truth = load_truth(index, index_dir, panel_dir, condition)
    cal, inputs = calibrate_condition(index, region_arrays, panel_dir, condition, config)
    out = {}
    rna_pmf = inputs["rna_fl_pmf"]
    out["shipped"] = (
        transcript_capture_eff_lengths(cal, region_arrays, index, truth.L_plain, rna_pmf),
        cal,
    )
    if (panel_dir / "oracle_cache" / condition / "slot_truth.npz").is_file():
        o = oracle_gdna_result(cal, region_arrays, panel_dir, condition, inputs["gdna_fl_pmf"])
        out["oracle_gdna"] = (
            transcript_capture_eff_lengths(o, region_arrays, index, truth.L_plain, rna_pmf),
            o,
        )
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
            f"   arm {arm:<12} reference {_fmt_ref(cal)}; anchored on {anchor} probed transcripts "
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


def main() -> int:
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
    ap.add_argument(
        "--set",
        dest="settings",
        action="append",
        default=[],
        metavar="SECTION.FIELD=VALUE",
        help="override one PipelineConfig field on every arm; repeatable",
    )
    ap.add_argument("--out", type=Path, default=None, help="the per-transcript table (one condition)")
    ap.add_argument("--scale", action="store_true", help="THE ONE-SCALE READ-OUT: L / Y per hypothesis class (one condition)")
    ap.add_argument("--width-step", type=int, default=1, help="--scale: thin the gDNA-space yield's width grid by this step")
    ap.add_argument("--self-test", action="store_true", help="perturb every comparator; no panel")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index_dir, panel_dir = PANELS[args.panel]
    if args.panel_dir is not None:
        panel_dir = args.panel_dir
    config = PipelineConfig()
    for spec in args.settings:
        config = set_field(config, spec)
    arms = load_arms(args.module)
    index = TranscriptIndex.load(str(index_dir))
    region_arrays = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)

    if args.scale:
        if not args.condition:
            raise SystemExit("⛔ --scale reads one condition: pass --condition")
        r = scale_readout(index, region_arrays, index_dir, panel_dir, args.condition, config, args.min_frags, args.width_step)
        print_scale(r)
        if args.out:
            pd.DataFrame(r["per_transcript"]).to_csv(args.out, sep="\t", index=False)
            pd.DataFrame(r["per_locus"]).to_csv(args.out.with_name(args.out.stem + "_loci" + args.out.suffix), sep="\t", index=False)
            print(f"   per-transcript and per-locus tables → {args.out} and its _loci twin")
        return 0
    if args.condition:
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
        return 0

    conds = sorted(
        p.name
        for p in (panel_dir / "oracle_cache").iterdir()
        if (p / "_main" / "payload.npz").is_file() and p.name.endswith("_capture_on")
    )
    if not conds:
        raise SystemExit(f"⛔ no cached capture-ON condition under {panel_dir / 'oracle_cache'}")
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
            print(f"{'condition':<44}{'reference':>12} " + " ".join(f"{n:>30}" for n in names))
        cal = factors["shipped"][1]
        ref = (
            "None"
            if cal.gdna_reference_density is None
            else f"10^{np.log10(cal.gdna_reference_density):+.2f}"
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
    return 0


# ── the one-scale read-out ───────────────────────────────────────────────────────────────────────

#: the class-mean tolerance the read-out judges against: the converged elasticity of the synthetic pool to its own
#: opportunity is about −20 at ``g50 ss.99 ON``, so a 1 % disagreement between two classes' means is a 20 % move
#: of that pool
ONE_SCALE_TOL = 0.01


def em_loci(index: TranscriptIndex, panel_dir: Path, condition: str, config):
    """The pipeline run to the prior assembly and stopped there: the shipped calibration, the region arrays, the
    EM's multi-loci, the shipped ``LocusPriors``, and the ruler's inputs and output — no EM."""
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
        cap["fl_models"] = models  # the last fit is calibration's: its gdna_pmf is the pairing arm's law
        return models

    def ruler_hook(calibration, region_arrays, index_, fl_eff_lengths, rna_fl_pmf):
        out = orig_ruler(calibration, region_arrays, index_, fl_eff_lengths, rna_fl_pmf)
        cap["fl"] = np.asarray(fl_eff_lengths, dtype=np.float64).copy()
        cap["rna_pmf"] = np.asarray(rna_fl_pmf, dtype=np.float64).copy()
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
        run_pipeline(str(panel_dir / condition / "sim_oracle.bam"), index, config)
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
    """Per locus, ``Y_g = Σ_w f(w) Σ_{starts s whose fragment OVERLAPS a block} (off + bind · best single-part
    overlap)``, read off the simulator's own gDNA-space landscape — the same computation the reads were drawn
    with. ``blocks`` is ``int64[n, 4]`` of ``(locus, ref_id, start, end)``. The width grid is thinned by
    ``width_step`` with each kept width carrying its bin's pmf mass (the landscape is per width per reference)."""
    Y = np.zeros(n_loci, dtype=np.float64)
    if blocks.size == 0:
        return Y
    keep = np.arange(0, sim.widths.size, max(int(width_step), 1))
    wk = sim.widths[keep]
    pk = np.add.reduceat(sim.pmf, keep)
    live = sim.sampler.enabled and sim.bind > 0.0
    for rid in np.unique(blocks[:, 1]):
        name = id2name[int(rid)]
        RL = int(sim.ref_lengths[name])
        bb = blocks[blocks[:, 1] == rid]
        for w, pw in zip(wk, pk, strict=True):
            w = int(w)
            lo = np.maximum(bb[:, 2] - w + 1, 0)
            hi = np.minimum(bb[:, 3] - 1, RL - w)
            n_starts = np.maximum(hi - lo + 1, 0).astype(np.float64)
            sw = np.zeros(bb.shape[0])
            if live:
                pos, wt = sim.sampler._extra_landscape("gdna", name, RL, w)
                if pos.size:
                    cs = np.concatenate(([0.0], np.cumsum(wt)))
                    a = np.searchsorted(pos, lo, side="left")
                    b = np.searchsorted(pos, hi, side="right")
                    sw = cs[b] - cs[a]
            np.add.at(Y, bb[:, 0], pw * (sim.off * n_starts + sim.bind * sw))
    return Y


def footprint_gdna_length(region_arrays, block: tuple, gdna_fl_pmf) -> float:
    """The shipped gDNA component's length at efficiency 1 over one footprint ``(ref_id, start, end)``:
    ``priors.assemble_priors`` itself, on a result carrying calibration's own geometry for this partition
    (the contained supports and ``calibrate._gdna_boundary_conserved_len``) and no reference."""
    from rigel.calibration.calibrate import _gdna_boundary_conserved_len
    from rigel.calibration.effective_length import contained_eff_length
    from rigel.calibration.priors import assemble_priors
    from rigel.calibration.region_arrays import boundary_region_indices
    from rigel.calibration.result import CalibrationResult
    from rigel.config import CalibrationConfig
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
        gdna_region_eff_len=contained_eff_length(length, gdna_fl_pmf), gdna_boundary_eff_len=ez,
        gdna_boundary_conserved_len=_gdna_boundary_conserved_len(region_arrays, gdna_fl_pmf),
        rna_region_eff_len=z, rna_boundary_eff_len=ez, gdna_frac_region=z, rna_pos_frac_region=z,
        rna_neg_frac_region=z, gdna_frac_boundary=ez, rna_pos_frac_boundary=ez, rna_neg_frac_boundary=ez,
        gdna_density_global=0.0, gdna_reference_density=None, gdna_reference_members=0,
        gdna_capture_efficiency_region=np.ones(n), gdna_capture_efficiency_boundary=np.ones(ne),
        rna_sense_frac=0.5, gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0,
        n_regions=n, n_boundaries=ne, n_sj=0, config=CalibrationConfig(),
    )
    ml = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        gdna_span=e0 - s0,
        loci=(Locus(ref=str(rid), ref_id=rid, start=s0, end=e0),),
    )
    return float(assemble_priors(cal, region_arrays, [ml]).gdna_eff_len[0])


def footprint_probed_fraction(sim: Sim, blocks: np.ndarray, n_loci: int, id2name: dict) -> np.ndarray:
    """Per locus: the fraction of the footprint's bases under a probe part, read off the sampler's own gDNA-space
    probe map (``_get_intervals("gdna", ref)``, the union of its parts) — the locus gDNA component's analogue of a
    transcript's probed fraction, so the locus rows can be read by probed class beside the transcript rows."""
    covered = np.zeros(n_loci, dtype=np.float64)
    length = np.zeros(n_loci, dtype=np.float64)
    if blocks.size == 0:
        return np.zeros(n_loci)
    for rid in np.unique(blocks[:, 1]):
        name = id2name[int(rid)]
        parts = sim.sampler._get_intervals("gdna", name) if sim.sampler.enabled else []
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
    """Per transcript: does a probe part in the transcript's OWN space strictly span one of its junctions? Read off
    the sampler's own probe map (``_get_intervals("mrna", row)``) against the junction positions in transcript
    coordinates, 5′ to 3′ — the panel's junction-probed transcripts, whose junction placements bind a probe no gDNA
    witness can see."""
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
        intervals = sim.sampler._get_intervals("mrna", i)
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
    """The largest pairwise disagreement among the class means, as a fraction, and its pair."""
    names = [k for k, v in means.items() if np.isfinite(v) and v > 0]
    worst, pair = 0.0, ""
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            d = abs(float(np.log(means[a] / means[b])))
            if d > worst:
                worst, pair = d, f"{a} / {b}"
    return float(np.expm1(worst)), pair


def scale_readout(index, region_arrays, index_dir, panel_dir, condition, config, min_frags: int, width_step: int):
    """The one-scale read-out for one condition: every hypothesis class's L / Y on the shipped lengths, the pairwise
    class-mean ratios and their verdict, the same by probed class, and the within-gene spread of the annotated
    transcripts."""
    t0 = time.time()
    truth = load_truth(index, index_dir, panel_dir, condition)
    sim = load_sim(index, index_dir, panel_dir, condition)
    Y_t = np.where(np.isfinite(truth.factor), truth.L_plain * truth.factor, np.nan)
    cap = em_loci(index, panel_dir, condition, config)
    loci, pri = cap["multi_loci"], cap["priors"]
    t = index.t_df
    syn = t["is_synthetic"].to_numpy(dtype=bool)
    single = (~syn) & t["is_nrna"].to_numpy(dtype=bool)
    multi = (~syn) & ~single
    frags = truth.frags
    # the locus gDNA component against the gDNA-space yield
    id2name = {v: k for k, v in index.ref_name_to_id.items()}
    n_loci = max((int(ml.multi_locus_id) for ml in loci), default=-1) + 1
    blocks = np.array(
        [(int(ml.multi_locus_id), int(b.ref_id), int(b.start), int(b.end)) for ml in loci for b in ml.loci if b.end > b.start],
        dtype=np.int64,
    ).reshape(-1, 4)
    Y_g = gdna_space_yield(sim, blocks, n_loci, id2name, width_step)
    L_g = np.asarray(pri.gdna_eff_len, dtype=np.float64)
    gdna_count = np.asarray(pri.gdna_count, dtype=np.float64)  # calibration's gDNA fragments: the component's weight
    ruler = cap["ruler"]
    jp = junction_probed_flags(sim, index)
    pf = truth.probed_frac
    pf_loc = footprint_probed_fraction(sim, blocks, n_loci, id2name)
    rows: list[tuple[str, dict]] = [
        ("locus gDNA component", scale_stats(L_g, Y_g, gdna_count)),
        ("synthetic span", scale_stats(ruler[syn], Y_t[syn], frags[syn])),
        ("annotated single-exon", scale_stats(ruler[single], Y_t[single], frags[single])),
        ("annotated multi-exon", scale_stats(ruler[multi], Y_t[multi], frags[multi])),
    ]
    by_class: dict[str, dict] = {}
    for cname, lo, hi in CLASSES:
        m = multi & (pf >= lo) & (pf < hi)
        ms = syn & (pf >= lo) & (pf < hi)
        ml = (pf_loc >= lo) & (pf_loc < hi)
        rows.append((f"  multi-exon {cname}, no junction probe", scale_stats(ruler[m & ~jp], Y_t[m & ~jp], frags[m & ~jp])))
        rows.append((f"  multi-exon {cname}, JUNCTION-PROBED", scale_stats(ruler[m & jp], Y_t[m & jp], frags[m & jp])))
        rows.append((f"  synthetic span {cname}", scale_stats(ruler[ms], Y_t[ms], frags[ms])))
        rows.append((f"  locus gDNA footprint {cname}", scale_stats(L_g[ml], Y_g[ml], gdna_count[ml])))
        by_class[cname] = {
            "isoforms": scale_stats(ruler[m], Y_t[m], frags[m]),
            "spans": rows[-2][1],
            "gDNA": rows[-1][1],
        }
    stats = dict(rows)
    means = {
        "gDNA": stats["locus gDNA component"]["geomean"],
        "spans": stats["synthetic span"]["geomean"],
        "isoforms": stats["annotated multi-exon"]["geomean"],
    }
    # within-gene spread of the annotated multi-exon transcripts: the sd of log(L / Y) about the gene mean
    live = multi & np.isfinite(Y_t) & (Y_t > 0) & (ruler > 0) & (frags >= min_frags)
    g = t["g_index"].to_numpy()
    within = float("nan")
    if live.any():
        d = pd.DataFrame(dict(g=g[live], x=np.log(ruler[live] / Y_t[live]), w=frags[live]))
        d["x"] -= d.groupby("g")["x"].transform("mean")
        within = float(np.sqrt(np.average(d["x"] ** 2, weights=d["w"])))
    elapsed = time.time() - t0
    return dict(condition=condition, stats=stats, means=means, by_class=by_class, within_gene_sd=within,
                n_junction_probed=int(jp[multi].sum()),
                n_multi=int(multi.sum()), n_loci=int(n_loci), elapsed=elapsed,
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


def self_test() -> int:
    """Perturb every comparator on a synthetic panel with a KNOWN partition — a tiny index with a fully
    probed, a half-probed and an unprobed gene, a probe file the sampler reads — and require each to
    fire (`TRAPS: perturb-every-gate`). The truth is the sampler's own; the arms are rulers written here."""
    import tempfile

    sys.path.insert(0, str(_REPO / "tests"))
    from _index_builder import build_test_index

    checks: list[tuple[str, bool]] = []

    def check(name: str, ok: bool):
        checks.append((name, bool(ok)))

    gtf = "".join(
        f'chr1\ttest\texon\t{s + 1}\t{s + 600}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n'
        for g, t, starts in (
            ("gP", "probed", (1000, 3000)),
            ("gH", "half", (6000, 8000)),
            ("gU", "unprobed", (11000, 13000)),
        )
        for s in starts
    )
    with tempfile.TemporaryDirectory() as tmp:

        class _Factory:
            def mktemp(self, name):
                p = Path(tmp) / name
                p.mkdir()
                return p

        index = build_test_index(_Factory(), gtf, genome_size=16000, name="ruler_self_test")
        bed = Path(tmp) / "probes.bed"
        # the probed gene's ONE probe spans both exons as a two-block BED12 line — contiguous in cDNA,
        # so every start of every fragment lies under it at full scale; the half-probed gene's probe
        # covers its first exon only
        bed.write_text(
            "chr1\t1000\t3600\tp\t0\t.\t1000\t3600\t0\t2\t600,600\t0,2000\n"
            "chr1\t6000\t6600\th\t0\t.\t6000\t6600\t0\t1\t600\t0\n"
        )
        transcripts = load_transcripts_from_index(index.index_dir)
        off, binding = 1.0, 10.0
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
        )
        params = type("Sim", (), dict(frag_mean=200, frag_std=50, frag_min=100, frag_max=300))
        widths, pmf = fl_pmf(params)
        L = np.array([int(t.length) for t in transcripts], dtype=np.int64)
        L_plain, factor, probed = sampler_truth(sampler, L, widths, pmf, off)
        t = index.t_df
        row = {tid: int(i) for i, tid in enumerate(t.t_id)}
        p, h, u = row["probed"], row["half"], row["unprobed"]

        # ① the truth: an unprobed transcript's factor is exactly 1, a fully probed one's is
        #    1 + binding·E[w]/off, a half-probed one's lies strictly between
        check("an unprobed transcript's truth factor is exactly 1", factor[u] == 1.0)
        mean_w = float(np.dot(widths, pmf))
        check(
            "a fully probed transcript's truth factor is 1 + binding·E[w]/off",
            abs(factor[p] / (1.0 + binding * mean_w / off) - 1.0) < 0.02,
        )
        check("a half-probed transcript's factor lies between", 1.0 < factor[h] < factor[p])
        # ② the class is read off the sampler: probed 1, unprobed 0, half-probed ½
        check("probed fraction: fully probed reads 1", abs(probed[p] - 1.0) < 1e-9)
        check("probed fraction: unprobed reads 0", probed[u] == 0.0)
        check("probed fraction: half-probed reads ½", abs(probed[h] - 0.5) < 1e-9)

        kind = np.where(
            t.is_synthetic.to_numpy(bool),
            "entity",
            np.where(t.is_nrna.to_numpy(bool), "single-exon", "mRNA"),
        )
        truth = Truth(L_plain, factor, probed, kind, np.full(L.size, 100.0), np.full(L.size, 10.0))

        # ③ a ruler that IS the truth scores 0 on every class, and a global scale on it is free
        d, measured, n_zero = score(factor, truth, min_frags=20)
        rows = class_rows(d, measured, truth)
        check("the truth scores exactly 0 on every class", all(abs(r["med"]) < 1e-12 for r in rows))
        check("every class is 100 % within the tolerance", all(r["within"] == 1.0 for r in rows))
        check("no zero factor on the truth", n_zero == 0)
        d3, _, _ = score(factor * 3.0, truth, min_frags=20)
        check("a global scale is anchored away", np.allclose(d3[measured], d[measured]))
        # ④ the unprobed class moves by log 2 when its factor is doubled, the probed class not at all
        f2 = factor.copy()
        f2[u] *= 2.0
        d2, m2, _ = score(f2, truth, min_frags=20)
        check("a doubled unprobed factor reads +log 2 on the unprobed class", abs(d2[u] - np.log(2.0)) < 1e-9)
        check("...and leaves the probed class at 0", abs(d2[p]) < 1e-12)
        h2 = headline(d2, m2, truth)
        check("the headline carries the unprobed median", abs(h2["unprobed_med"] - np.log(2.0)) < 1e-9)
        # ⑤ a factor of exactly 0 is counted and dropped, never a −inf in the quantiles
        f0 = factor.copy()
        f0[u] = 0.0
        d0, m0, z0 = score(f0, truth, min_frags=20)
        check("a zero factor is counted", z0 == 1)
        check("...and dropped from the measured set", not m0[u] and np.isfinite(d0[m0]).all())
        # ⑥ the fragment floor: a transcript under it is not scored
        few = dataclasses.replace(truth, frags=np.where(np.arange(L.size) == u, 5.0, 100.0))
        _, m5, _ = score(factor, few, min_frags=20)
        check("a transcript under the fragment floor is not measured", not m5[u] and m5[p])
        # ⑦ the anchor is the probed class: with the probed class biased, everything shifts by its bias
        fb = factor.copy()
        fb[p] *= 5.0
        db, _, _ = score(fb, truth, min_frags=20)
        check("the anchor is the probed class's median", abs(db[p]) < 1e-12 and abs(db[u] + np.log(5.0)) < 1e-9)
        # ⑧ the module loader refuses a file without ARMS, and a reserved arm name
        bad = Path(tmp) / "bad.py"
        bad.write_text("X = 1\n")
        try:
            load_arms(bad)
            fired = False
        except SystemExit:
            fired = True
        check("a module without ARMS is refused", fired)
        reserved = Path(tmp) / "reserved.py"
        reserved.write_text("ARMS = {'shipped': lambda *a: None}\n")
        try:
            load_arms(reserved)
            fired = False
        except SystemExit:
            fired = True
        check("a prototype arm named `shipped` is refused", fired)
        good = Path(tmp) / "good.py"
        good.write_text("def r(cal, ra, index, fl, **inputs):\n    return fl\nARMS = {'plain': r}\n")
        check("a well-formed prototype loads", list(load_arms(good)) == ["plain"])
        # ⑨ the one-scale read-out: a class priced at its own yield reads 1 and disagrees with nothing; a class
        #    scaled by 1.185 (the prototype rule's span:isoform ratio) is named as the pair, and a per-mille
        #    difference is inside the tolerance
        Y = L_plain * factor
        w = np.full(L.size, 100.0)
        s_true = scale_stats(Y, Y, w)
        check("a class priced at its own yield reads L / Y = 1 with no spread", abs(s_true["geomean"] - 1.0) < 1e-12 and s_true["sd_log"] < 1e-12)
        means = {"gDNA": 1.0, "spans": 1.185, "isoforms": 1.0}
        worst, pair = scale_verdict(means)
        check("a 1.185 scale on one class is named as its pair", abs(worst - 0.185) < 1e-9 and "spans" in pair)
        check("...and is outside the tolerance, while a per-mille difference is inside", worst > ONE_SCALE_TOL and scale_verdict({"a": 1.0, "b": 1.001})[0] < ONE_SCALE_TOL)
        # ⑩ the junction-probed flag is read off the sampler's own probe map: the probed gene's two-block probe
        #    spans its junction, the half-probed gene's single-exon probe does not, the unprobed gene has none
        sim = Sim(sampler, np.asarray(widths, dtype=np.int64), np.asarray(pmf, dtype=np.float64), off, binding, {"chr1": 16000}, transcripts)
        jp = junction_probed_flags(sim, index)
        check("a probe spanning a transcript's junction flags it, and only it", bool(jp[p]) and not jp[h] and not jp[u])
        # ⑪ the gDNA-space yield counts every start whose fragment OVERLAPS the footprint (the locus gDNA component
        #    takes its outer boundaries), so a block over the unprobed gene reads exactly off · Σ_w f(w)(L + w − 1):
        #    no probe, every start at the off-target weight; a block over the probed gene reads more
        blocks = np.array([(0, 0, 11000, 13600), (1, 0, 1000, 3600)], dtype=np.int64)
        Yg = gdna_space_yield(sim, blocks, 2, {0: "chr1"}, width_step=1)
        plain = float(sum(pw * off * (2600 + int(wv) - 1) for wv, pw in zip(widths, pmf, strict=True)))
        check("an unprobed footprint's gDNA-space yield is its off-target plain length", abs(Yg[0] / plain - 1.0) < 1e-9)
        check("...and a probed footprint's yield exceeds it", Yg[1] > Yg[0] * 2.0)
        # ⑫ the shipped gDNA component's length (`assemble_priors`) at efficiency 1 counts the SAME overlapping
        #    starts as the yield, Σ_w f(w)(L + w − 1), on a footprint whose outside neighbours hold every outside
        #    base: each crossing fragment then starts in a neighbour, so its whole unit lands on the footprint's
        #    objects
        from rigel.calibration.region_arrays import RegionArrays
        ra_t = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
        # the geometry indexes a pmf BY WIDTH (calibration's convention); the simulator's pair is positional
        by_width = np.zeros(int(max(widths)) + 1)
        by_width[np.asarray(widths, dtype=np.int64)] = np.asarray(pmf, dtype=np.float64)
        Lfp = footprint_gdna_length(ra_t, (0, 11000, 13600), by_width)
        check("the gDNA component's length at efficiency 1 is the yield's overlapping-start count", abs(Lfp / (plain / off) - 1.0) < 1e-9)
        # ⑬ the footprint's probed fraction is read off the sampler's own gDNA-space probe map: the unprobed gene's
        #    block reads 0, a block exactly on one probe part reads 1, a block half on it reads ½
        parts = sim.sampler._get_intervals("gdna", "chr1")
        a, b = int(parts[0].start), int(parts[0].end)
        pf_loc = footprint_probed_fraction(sim, np.array([(0, 0, 11000, 13600), (1, 0, a, b), (2, 0, a - (b - a), b)], dtype=np.int64), 3, {0: "chr1"})
        check("an unprobed footprint's probed fraction is 0", pf_loc[0] == 0.0)
        check("a footprint exactly on a probe part reads 1", abs(pf_loc[1] - 1.0) < 1e-12)
        check("...and one half on it reads ½", abs(pf_loc[2] - 0.5) < 1e-12)

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
