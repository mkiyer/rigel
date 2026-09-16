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
unprobed class is where the floor and the reference bite, and the partial classes are where the
junction rule lives. The gDNA witness under-reads a probe that spans a junction or a tiny exon by the
panel's own geometry, declared and not repaired.

Usage::

    python scripts/design/ruler_vs_truth.py --panel test --condition gdna_g05_ss_0.99_nrna_file_capture_on
    python scripts/design/ruler_vs_truth.py --panel test                      # every capture-ON condition, one line each
    python scripts/design/ruler_vs_truth.py --panel test --panel-dir ~/Downloads/rigel_runs/test_reference/scenarios_depth_d10
    python scripts/design/ruler_vs_truth.py --panel ladder --condition C --module proto.py --out table.tsv
    python scripts/design/ruler_vs_truth.py --panel test --condition C --set calibration.message_policy=silent
    python scripts/design/ruler_vs_truth.py --self-test
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
    "gdna_split_penalty",
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


def load_truth(index: TranscriptIndex, index_dir: Path, panel_dir: Path, condition: str) -> Truth:
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
    L = np.array([int(t.length) for t in transcripts], dtype=np.int64)
    L_plain, factor, probed = sampler_truth(sampler, L, widths, pmf, float(cap_cfg["off_target_weight"]))
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


def oracle_gdna_result(cal, region_arrays, panel_dir: Path, condition: str):
    """The same result with the ruler's gDNA masses replaced by the CERTIFIED truth per object and the
    reference read off a landscape fitted on the truth: the ideal witness, no estimator noise."""
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
    return dataclasses.replace(
        cal,
        count_gdna_region=m_reg,
        count_gdna_boundary=m_bnd,
        gdna_reference_density=ref,
        gdna_reference_members=0 if mode is None else mode.n_members,
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
    out["shipped"] = (
        transcript_capture_eff_lengths(cal, region_arrays, index, truth.L_plain),
        cal,
    )
    if (panel_dir / "oracle_cache" / condition / "slot_truth.npz").is_file():
        o = oracle_gdna_result(cal, region_arrays, panel_dir, condition)
        out["oracle_gdna"] = (transcript_capture_eff_lengths(o, region_arrays, index, truth.L_plain), o)
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
                gdna_split_penalty=0.2,
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
