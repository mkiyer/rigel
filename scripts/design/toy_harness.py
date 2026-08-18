#!/usr/bin/env python
"""⭐⭐ **THE TOY HARNESS — one transcript at a time, under a real library's global conditions.**

⛔ **THE PROBLEM THIS SOLVES.** Every calibration defect we have found so far was found by sifting a
36-condition, 10-million-fragment panel for objects that behave badly, and then trying to reason
backwards from a region id to a mechanism. That is slow, the examples are never quite the ones you
wanted, and a fix cannot be *demonstrated* — only measured in aggregate, where two errors can cancel.

⭐ **WHAT IT DOES.** A **mini chromosome** you define — a handful of genes with exactly the structure
you want to interrogate — simulated, scanned and calibrated in seconds, with the per-object answer
printed beside per-object truth. Small enough to read every row.

⭐⭐ **AND THE GLOBALS COME FROM A REAL CACHED CONDITION, NOT FROM DEFAULTS.** That is the whole
trick. A tiny toy cannot fit the library-level quantities calibration needs — a single transcript has
no population to estimate a strand balance, an enrichment landscape or an intergenic background from.
So one of the 36 ladder conditions acts as **donor**: it is calibrated once, and its fitted
library-level bundle is injected into the toy (:class:`~rigel.calibration.calibrate.InjectedCalibrationPriors`,
whose docstring specifies exactly this use). The toy then supplies only the controlled per-region
GEOMETRY, which is the thing under study.

The donor supplies, and the toy therefore does NOT have to invent:

===========================  ================================================================
gDNA + RNA length pmfs       the fitted models, passed as ``calibrate`` kwargs
strand balance kappa         plus both Beta-Binomial overdispersions and the Fisher noise-floor
                             sample sizes -- so the deadband behaves as it does on real data
enrichment landscape         the density NPMLE (an ABSOLUTE log-density model, hence the depth
                             matching below)
intergenic gDNA background   the intron-factory reference and the aggregate rho_bg
capture on/off + its knobs   reproduced in the toy's own SIMULATION, with toy probes
fragment-length simulation   frag mean/sd/min/max, read length, strand specificity
gDNA density per base        ⭐ see `_rate_from_capture` -- the toy is simulated to MATCH it
===========================  ================================================================

⭐⭐ **DEPTH MATCHING IS NOT OPTIONAL.** ``calibrate``'s own note on the injected enrichment prior says
it: the NPMLE is an *absolute* log-density landscape, so "the toy's densities must be generated at the
reference library's depth so its regions project onto the right cells". A toy at the wrong depth is not
a small version of the library, it is a different library. So the harness measures the donor's gDNA
density per base from the donor's OWN pure-gDNA population and simulates the toy to match it, rather
than taking a fragment count from the user.

⭐⭐ **TERMINOLOGY — one word per concept, owner-set 2026-08-04.**

=================  ==============================================================================
**counts**         discrete INTEGER fragment counts. What the accumulator stores and what the
                   solver's Poisson ``n`` is
**density**        counts per base. ⭐ **"abundance" is the SAME THING and the two are used
**abundance**      interchangeably** -- both mean counts/bp. ⛔ Not the simulator's molar
                   ``abundance=`` field, which is a per-transcript weight, not a density
=================  ==============================================================================

⚠ A region's stored *counts* are CONTAINED counts, so ``counts != density x bp`` -- it is
``density x effective_length``. The sweep therefore reports the density it ASKED the simulator for
and the counts each object actually RECEIVED, side by side, and never converts one into the other.

Usage::

    # list the built-in specs, then run one against a donor condition
    python scripts/design/toy_harness.py --list
    python scripts/design/toy_harness.py --spec two_exon --donor gdna_g50_ss_0.50_nrna_none_capture_off
    python scripts/design/toy_harness.py --spec all --donor <cond>     # the whole ladder of specs

Gates: ``tests/calibration/test_toy_harness.py``.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))

from _oracle import OracleTruth  # noqa: E402

from rigel.calibration.calibrate import InjectedCalibrationPriors, calibrate  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import (  # noqa: E402
    _drain_side_buffer,
    _native_detect_sj_tag,
    scan_and_buffer,
)
from rigel.scan_cache import index_derived_inputs  # noqa: E402
from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario  # noqa: E402

_EPS = 1.0e-12

DEFAULT_SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
DEFAULT_INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"

#: Coarse region type names, shared with every other instrument (`signature.coarse_type_array`).
TYPE_NAMES = {0: "intergenic", 1: "intron", 2: "exon"}


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# ⭐⭐⭐ THE RELAY'S DIAGNOSTIC BANKS — AND WHETHER THIS RUN HAS ANY
# ──────────────────────────────────────────────────────────────────────────────────────────────────
#
# ⛔⛔ **`_uni` HAS EXACTLY ONE WRITER AND IT IS BEHIND A CONFIG FLAG.** `messages/head.py` appends it,
# so it exists only under `HeadPolicy`; the shipped `CalibrationConfig.message_propagation` is `False`,
# which installs `SilentPolicy`. Five instruments read `capture["_uni"]` unconditionally and therefore
# died with `KeyError: '_uni'` the day that default flipped, while the SUITE stayed green because the
# test readers install the policy themselves. That is `TRAPS: a-green-suite-hid-five-dead-instruments`
# recurring, and this block is the shared repair so it cannot recur once per instrument.
#
# ⭐ The precedent is `ladder_arm_ab.py`: the message setting is PART OF THE ARM, it is STAMPED into the
# output, and an arm the policy cannot express is REFUSED up front rather than reported as a null result
# (`TRAPS: an-ablation-that-never-ran`). So every reader of these banks takes `--messages {off,on}`,
# prints `messages_stamp()`, and either degrades its relay-derived column with `relay_live()` or calls
# `require_relay()`.

#: ⭐ what the tool SHIPS, re-derived rather than written down — the honest default for an instrument
#: whose headline question does not involve the relay at all.
MESSAGES_SHIPPED: bool = bool(CalibrationConfig().message_propagation)

#: ⭐⭐ The `_uni_static` keys that are NOT the relay's to own: `head.py` copies them straight off its
#: `StepContext`, and `sweep.py` — the backbone, which runs under EVERY policy — already publishes the
#: same arrays at the top level of the capture under different names. So these four survive the mute,
#: and `relay_static()` serves them from the policy-independent source.
#: ⛔ The mapping is ASSERTED, never assumed: where both names are present (any `--messages on` run)
#: `relay_static()` requires byte-identity, so an alias that stops being an alias fails loudly instead
#: of quietly substituting a different quantity.
_STATIC_FROM_BACKBONE: dict[str, str] = {
    "M": "mass_global",  # head.py: M = ctx.mass  ·  sweep.py: mass=mass_global
    "E_g": "eff_global",  # head.py: E_g = ctx.eff_gdna_global  ·  sweep.py: eff_gdna_global=eff_global
    "E_r": "eff_rna",  # head.py: E_r = ctx.eff_rna  ·  sweep.py: eff_rna=ER
    "tau_own": "_tau0_lam",  # head.py: tau_own = ni.tau_lam  ·  sweep.py: _tau0_lam=own.tau_lam
}


def add_messages_flag(ap, *, default: bool) -> None:
    """Give an instrument the `--messages {off,on}` flag, with ITS OWN honest default.

    ⭐ `default=MESSAGES_SHIPPED` for an instrument whose measurement is policy-independent — it then
    runs the configuration the tool ships. `default=True` for one whose measurement IS the relay, which
    is a deliberate departure from the shipped config and is stamped as such on every run."""
    ap.add_argument(
        "--messages",
        choices=("off", "on"),
        default="on" if default else "off",
        help=f"belief propagation across objects. This instrument defaults to "
        f"{'on' if default else 'off'}; the SHIPPED CalibrationConfig is "
        f"{'on' if MESSAGES_SHIPPED else 'off'}. Stamped into the output either way.",
    )


def messages_on(args) -> bool:
    """The parsed flag as a bool."""
    return args.messages == "on"


def with_messages(config: CalibrationConfig, messages: bool) -> CalibrationConfig:
    """The same config with the relay switched. ⭐ No monkeypatching: `calibrate` selects
    `HeadPolicy` vs `SilentPolicy` off this one field."""
    return dataclasses.replace(config, message_propagation=bool(messages))


def messages_stamp(messages: bool) -> str:
    """The line every one of these instruments prints, so a reader can never mistake which
    configuration produced the numbers below it."""
    if messages:
        note = (
            "  ⛔ NOT the shipped config (shipped is OFF) — HeadPolicy installed"
            if not MESSAGES_SHIPPED
            else "  (the shipped config)"
        )
        return f"   messages = ON   ·   policy = HeadPolicy{note}"
    note = (
        "  (the shipped config)"
        if not MESSAGES_SHIPPED
        else "  ⛔ NOT the shipped config (shipped is ON)"
    )
    return (
        f"   messages = OFF  ·   policy = SilentPolicy{note}\n"
        f"   ⛔ ψ carries each slot's OWN evidence alone. Every relay-derived column below reads "
        f"'—' because the relay SENT NOTHING — not because it sent a zero."
    )


def relay_live(capture) -> bool:
    """Did a relay run in THIS capture? — read off the artifact, not off the config.

    ⭐ The artifact is the honest witness: `head.py` is the only writer of `_uni`, so its presence is
    the relay's own signature. An instrument that trusted its own `--messages` flag instead would lie
    the day a caller forgot to thread the config through (`TRAPS: an-ablation-that-never-ran`)."""
    return bool(capture.get("_uni"))


def relay_channels(capture) -> dict | None:
    """The LAST relay fuse's per-slot channels (`cg`, `cm_g`, `cm_p`, `cm_n`, `c_tau`, `mo_*`, …), or
    `None` when the relay was muted. ⛔ `None` means *no claim was made*; it does not mean zero."""
    return capture["_uni"][-1] if relay_live(capture) else None


def relay_static(capture) -> dict:
    """`_uni_static`, with the four backbone-owned entries backfilled so they survive the mute.

    ⛔ Every key `head.py` alone publishes (`og`, `pg_own`, `rho_lo`, `rho_hi`, `struct_lock`, …) is
    absent under `SilentPolicy` and stays absent — a caller that needs one must `require_relay()`."""
    st = dict(capture.get("_uni_static", {}))
    for head_key, backbone_key in _STATIC_FROM_BACKBONE.items():
        if backbone_key not in capture:
            continue
        alias = np.asarray(capture[backbone_key], np.float64)
        if head_key in st:
            # the alias, FALSIFIED on every run that has both names
            np.testing.assert_array_equal(
                np.asarray(st[head_key], np.float64),
                alias,
                err_msg=f"_uni_static[{head_key!r}] is no longer capture[{backbone_key!r}] — "
                f"`_STATIC_FROM_BACKBONE` in toy_harness.py is stale and must be re-derived",
            )
        else:
            st[head_key] = alias
    return st


def require_relay(capture, *, what: str) -> dict:
    """The measurement IS the relay, so refuse rather than report a null result.

    ⭐ Same contract as `ladder_arm_ab.py`'s up-front refusal of an arm the policy cannot express: a
    section that cannot run must say so in the reader's own vocabulary, not print zeros."""
    uni = relay_channels(capture)
    if uni is None:
        raise SystemExit(
            f"⛔ {what} is a MESSAGE-LAYER measurement and this run has message_propagation=False "
            f"(SilentPolicy), so the relay published no `_uni` banks at all.\n"
            f"   Re-run with `--messages on`. ⚠ That is NOT the shipped configuration "
            f"(CalibrationConfig.message_propagation = {MESSAGES_SHIPPED}) and the stamp will say so."
        )
    return uni


def relay_silent_note(what: str) -> str:
    """The one-line banner a DEGRADED section prints in place of its relay-derived numbers."""
    return (
        f"   ⛔ {what} NEEDS THE RELAY AND THE RELAY IS MUTED (SilentPolicy, the shipped config).\n"
        f"      Not measured here, and NOT a zero measurement. `--messages on` to measure it."
    )


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE DONOR — everything the toy cannot fit for itself
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class DonorGlobals:
    """The library-level conditions harvested from one real cached condition.

    ⚠ **Deliberately NOT cached to disk.** The bundle is a function of the calibration code that fit
    it, so a stored copy goes stale on exactly the changes this harness exists to test, and a stale
    global is invisible — it does not crash, it just quietly answers a different question. Harvesting
    costs one scan plus one calibrate (~30 s); harvest once per session and run many toys against it.
    """

    condition: str
    priors: InjectedCalibrationPriors
    gdna_fl_pmf: np.ndarray
    rna_fl_pmf: np.ndarray
    #: ⭐ gDNA molecules per base, measured on the donor's own structurally-pure-gDNA population.
    gdna_rate_per_base: float
    #: the donor's read-simulation settings, so the toy draws from the same library
    frag_mean: float
    frag_std: float
    frag_min: int
    frag_max: int
    read_length: int
    strand_specificity: float
    #: capture: whether the donor had it on, and the numeric knobs to reproduce it on toy probes
    capture_on: bool
    capture_knobs: dict = field(default_factory=dict)

    def describe(self) -> str:
        return (
            f"donor={self.condition}\n"
            f"  kappa={self.priors.rna_sense_frac:.6f}  n_rna_obs={self.priors.n_rna_obs:,.0f}  "
            f"n_gdna_obs={self.priors.n_gdna_obs:,.0f}\n"
            f"  od_rna={self.priors.rna_strand_overdispersion:.4g}  "
            f"od_gdna={self.priors.gdna_strand_overdispersion:.4g}\n"
            f"  gDNA rate = {self.gdna_rate_per_base:.6g} molecules/base   "
            f"FL {self.frag_mean:.0f}+-{self.frag_std:.0f} [{self.frag_min},{self.frag_max}] "
            f"read {self.read_length}\n"
            f"  strand_specificity={self.strand_specificity:.4g}  capture={'ON' if self.capture_on else 'off'}\n"
            f"  enrichment_prior={'yes' if self.priors.enrichment_prior is not None else 'NONE'}  "
            f"intron_background={'yes' if self.priors.intron_background is not None else 'NONE'}  "
            f"background={'yes' if self.priors.background is not None else 'NONE'}"
        )


def harvest(
    donor_dir: Path,
    index: TranscriptIndex,
    *,
    config: CalibrationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
    bam: str | None = None,
    name: str | None = None,
) -> DonorGlobals:
    """Calibrate one cached condition and keep everything the toy cannot fit for itself.

    ``bam`` / ``name`` override the defaults (``donor_dir/sim_oracle.bam`` and ``donor_dir.name``) so a
    gate can harvest from a scenario it built itself rather than from a 10 M-fragment panel condition
    that only exists on one machine.
    """
    config = config or CalibrationConfig()
    pipeline_config = pipeline_config or PipelineConfig()
    bam = bam or str(donor_dir / "sim_oracle.bam")
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, _buf, payload = scan_and_buffer(bam, index, scan)

    ra = RegionArrays.from_index(index)
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    debug: dict = {}
    calibrate(
        payload=payload,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=config,
        _debug=debug,
        **index_derived_inputs(index),
    )

    # the gDNA rate, from the pure population, in the SAME frame the toy will be measured in
    rate = _rate_from_capture(debug["capture"], debug["chain"], ra)

    name = name or donor_dir.name
    sim = _donor_sim_params(donor_dir, name)
    return DonorGlobals(
        condition=name,
        priors=debug["calibration_priors"],
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        gdna_rate_per_base=rate,
        capture_on="capture_on" in name,
        **sim,
    )


def _rate_from_capture(capture, chain, region_arrays) -> float:
    """``sum(count) / sum(eff_gdna)`` over the donor's INTERGENIC region slots.

    ⭐ Both arrays are the solver's own (``capture['count']`` and ``capture['eff_gdna']``), so the rate
    is in exactly the frame the toy's own regions will be measured in — no second implementation of an
    effective length, which is how two definitions of one quantity start to drift (TRAPS: two-docstrings-one-quantity).
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_region = kind == REGION
    pure = is_region.copy()
    pure[is_region] = rtype[obj[is_region]] == 0
    count = np.asarray(capture["count"], np.float64).sum(axis=1)
    eff = np.asarray(capture["eff_gdna"], np.float64)
    if not pure.any() or eff[pure].sum() <= 0.0:
        raise ValueError("no intergenic region slots with gDNA opportunity in the donor")
    return float(count[pure].sum() / eff[pure].sum())


def _donor_sim_params(donor_dir: Path, name: str) -> dict:
    """The donor's read-simulation settings.

    ⚠ Fragment lengths come from ``truth_summary.json``'s **post-capture** measurement where it
    exists, never from a configured ``frag_mean`` — capture selects for length, so the configured
    parameters describe a library that was never sequenced (TRAPS: capture-selects-for-length). Strand specificity is read
    off the condition name, which is where the panel encodes it.
    """
    ss = 0.5
    for part in name.split("_"):
        try:
            v = float(part)
        except ValueError:
            continue
        if 0.0 <= v <= 1.0:
            ss = v
            break
    summary = donor_dir / "truth_summary.json"
    frag_mean, frag_std, frag_min, frag_max = 206.0, 98.0, 50, 500
    if summary.is_file():
        d = json.loads(summary.read_text())
        allrow = d.get("fragment_lengths", {}).get("all", {})
        if allrow.get("mean") is not None:
            frag_mean = float(allrow["mean"])
            frag_std = float(allrow["std"])
            frag_min = int(allrow["min"])
            frag_max = int(allrow["max"])
    return {
        "frag_mean": frag_mean,
        "frag_std": frag_std,
        "frag_min": frag_min,
        "frag_max": frag_max,
        "read_length": 100,
        "strand_specificity": ss,
        "capture_knobs": {
            "off_target_weight": 1.0,
            "binding_per_base": 10.0,
            "gdna_split_penalty": 0.2,
            "min_overlap": 1,
            "probe_length": 120,
            "capture_fraction": 1.0,
        },
    }


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE TOY
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ToySpec:
    """A mini chromosome: what to put on it, and how much RNA to express.

    ⚠ **There is no gDNA knob.** The gDNA level is not a free parameter — it is set to match the
    donor's measured density per base, because the injected enrichment landscape is absolute. The RNA
    side is the experimental knob, and ``abundance`` on each transcript is what you vary.
    """

    name: str
    what_it_probes: str
    genome_length: int
    genes: list[dict]
    n_rna_fragments: int = 40_000
    nrna_abundance: float = 0.0
    seed: int = 7
    #: transcript ids to put a probe on when the donor is capture-ON. ``None`` = all of them.
    captured: tuple[str, ...] | None = None


@dataclass
class ToyResult:
    spec: ToySpec
    donor: DonorGlobals
    result: object  #: the CalibrationResult
    truth: OracleTruth
    payload: object
    region_arrays: RegionArrays
    chain: object
    capture: dict
    n_gdna_target: int
    seconds: float
    #: ⭐ the toy's OWN index, so a caller can rebuild anything index-derived — in particular
    #: `splice_graph.build_boundary_flags_array`, the TSS/TES/DONOR/ACCEPTOR bits per BOUNDARY.
    index: object = None


def _toy_probes(spec: "ToySpec", out: Path, knobs: dict) -> str:
    """Probes over the toy's transcripts, written from the SPEC's own exon coordinates.

    ⛔ **Deliberately not ``write_random_capture_probes``.** That draws a random subset of genes, which
    is right for a 36-condition panel and wrong for a toy: the whole point of a toy is that you decide
    what is captured, and a random draw would make the condition depend on a seed. Here every
    transcript named in ``spec.captured`` (default: all of them) gets ONE probe centred in its own
    transcript coordinates, which makes probe placement part of the controlled geometry.

    The file is the transcript-coordinate TSV the sampler already reads: ``transcript_id start end``,
    coordinates in TRANSCRIPT space (`sampler._load_transcript_probes`).
    """
    plen = int(knobs.get("probe_length", 120))
    lines = ["transcript_id\tstart\tend"]
    for gene in spec.genes:
        for t in gene["transcripts"]:
            if spec.captured is not None and t["t_id"] not in spec.captured:
                continue
            # ⭐⭐ TILED, not one centred probe. A single central probe leaves a 2 kb transcript's ENDS
            # uncovered, so its `intergenic|exon` boundaries stay at off-target density -- and a 0-bp line's
            # counts are `density x mean_FL` no matter how long the chromosome is, so those objects then
            # have ~0.26 counts and the chain has nothing to propagate. The donor panels are built with
            # `design_suite_probes.py` at probe_density 1.0, i.e. tiled, which is what enriches a
            # first/last exon's boundary in the first place. Match that.
            #
            # ⭐⭐⭐ AND TILED **PER EXON**, SO EVERY PROBE ABUTS THE intron|exon BOUNDARIES AND NONE STRADDLES
            # A SJ. Probes are written in TRANSCRIPT space, so a probe spanning an internal sj
            # offset has a genomic footprint in TWO blocks -- and `capture/sampler._split_scale` then
            # multiplies every gDNA fragment overlapping it by ``gdna_split_penalty``. Tiling across the
            # whole transcript put a split probe over each internal sj, which SUPPRESSED exactly the
            # population that spans an intron|exon BOUNDARY: measured on `spliced_exons` x
            # `g75 ss0.50 capture_on`, the two boundaries carried 2 and 5 gDNA counts while the exon interiors
            # carried 16 and 14. Per-exon tiling leaves every probe inside one exon, so it is unsplit, it
            # ends exactly ON the boundary, and an boundary-crossing fragment takes the full binding weight for
            # its exon-side overlap.
            # ⚠ It is also the honest geometry for the question: a probe boundary at an exon end is what a
            # real panel produces, and the split-probe case is a SEPARATE population worth its own rung
            # rather than an accident of how a tiling loop was written.
            off = 0
            for s, e in t["exons"]:
                elen = e - s
                for start in range(0, elen, plen):
                    lines.append(f"{t['t_id']}\t{off + start}\t{off + min(start + plen, elen)}")
                off += elen
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    return str(out)


def run_toy(
    spec: ToySpec,
    donor: DonorGlobals,
    work_dir: Path,
    *,
    config: CalibrationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
) -> ToyResult:
    """Simulate → scan → split by origin → calibrate with the donor's globals injected."""
    import time

    t0 = time.perf_counter()
    config = config or CalibrationConfig()
    pipeline_config = pipeline_config or PipelineConfig()
    wd = Path(work_dir) / spec.name
    wd.mkdir(parents=True, exist_ok=True)

    # ── the gDNA level is DERIVED, not chosen: match the donor's molecules per base ──────────────
    n_gdna = int(round(donor.gdna_rate_per_base * spec.genome_length))
    gdna_fraction = n_gdna / max(spec.n_rna_fragments, 1)

    sim_cfg = ReadSimConfig(
        frag_mean=int(round(donor.frag_mean)),
        frag_std=int(round(donor.frag_std)),
        frag_min=donor.frag_min,
        frag_max=donor.frag_max,
        read_length=donor.read_length,
        strand_specificity=donor.strand_specificity,
        seed=spec.seed,
    )
    gdna_cfg = GDNAConfig(
        abundance=0.0,
        frag_mean=int(round(donor.frag_mean)),
        frag_std=int(round(donor.frag_std)),
    )

    sc = Scenario(
        spec.name,
        genome_length=spec.genome_length,
        seed=spec.seed,
        work_dir=wd / "sim",
    )
    for gene in spec.genes:
        sc.add_gene(gene["gene_id"], gene["strand"], gene["transcripts"])

    capture_cfg = None
    if donor.capture_on:
        probes = _toy_probes(spec, wd / "probes.tsv", donor.capture_knobs)
        k = donor.capture_knobs
        capture_cfg = CaptureConfig(
            probes=probes,
            probe_format="transcript",
            off_target_weight=float(k["off_target_weight"]),
            binding_per_base=float(k["binding_per_base"]),
            gdna_split_penalty=float(k["gdna_split_penalty"]),
            min_overlap=int(k["min_overlap"]),
        )

    res = sc.build_oracle(
        n_rna_fragments=int(spec.n_rna_fragments),
        gdna_fraction=gdna_fraction,
        nrna_abundance=float(spec.nrna_abundance),
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        capture_config=capture_cfg,
    )

    bam = str(res.bam_path)
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, _buf, pass_one = scan_and_buffer(bam, res.index, scan)
    ra = RegionArrays.from_index(res.index)

    # ── ⭐⭐ THE SECOND PASS, on the WHOLE, exactly as production runs it ─────────────────────────
    # ⛔ Until this landed every number the toy reported was an UNDRAINED tally, and the population it
    # understates is the spliced one: a held fragment is held precisely because its unsequenced gap
    # admits more than one intron path. On `tes_readthrough` TA's and TB's sj share the donor at
    # 2,000 and differ only in acceptor, so 0.65 % of fragments defer — and they understate the
    # certified channel at @9,100 by a measured 13.5 %.
    # ⭐ `_lift` is what makes the ORACLE valid afterwards: score and draw ONCE on the whole, then replay
    # each fragment's chosen hypothesis inside whichever origin partition holds it (TRAPS: draining-breaks-the-oracle).
    lift: dict = {}
    payload = _drain_side_buffer(
        pass_one, res.index, strand_model, seed=pipeline_config.second_pass_seed, _lift=lift
    )
    truth = OracleTruth.from_bam(
        bam,
        res.index,
        pipeline_config,
        wd / "split",
        spec.name,
        # ⛔ DRAINED whole here (it is what calibration reads and what sum-to-full is asserted against);
        # UNDRAINED whole inside `drain_with` (the drained bank holds nothing, so it has no key pool).
        full_payload=payload,
        drain_with=(
            (lift["undrained"], lift["choices"], lift["region_types"], lift["sj"])
            if lift
            else None
        ),
    )

    debug: dict = {}
    out = calibrate(
        payload=payload,
        strand_model=strand_model,
        gdna_fl_pmf=donor.gdna_fl_pmf,
        rna_fl_pmf=donor.rna_fl_pmf,
        config=config,
        injected_priors=donor.priors,
        _debug=debug,
        **index_derived_inputs(res.index),
    )
    return ToyResult(
        spec=spec,
        donor=donor,
        result=out,
        truth=truth,
        payload=payload,
        region_arrays=ra,
        chain=debug["chain"],
        capture=debug["capture"],
        n_gdna_target=n_gdna,
        seconds=time.perf_counter() - t0,
        index=res.index,
    )


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE REPORT — every object, because there are few enough to read
# ──────────────────────────────────────────────────────────────────────────────────────────────────


def object_rows(r: ToyResult) -> list[dict]:
    """One row per region and per contiguous boundary, in genomic order along the chain."""
    chain, ra = r.chain, r.region_arrays
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    start = np.asarray(ra.start, np.int64)
    size = np.asarray(ra.region_size_bp, np.float64)

    cap = r.capture
    fg_loc = np.asarray(cap["fg_loc"], np.float64)
    fg = np.asarray(cap["f_g"], np.float64)
    var_g = np.asarray(cap["var_g"], np.float64)
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    count = np.asarray(cap["count"], np.float64).sum(axis=1)
    mature = np.asarray(cap["mature"], np.float64)
    spliced = np.asarray(cap["spliced"], np.float64)

    ov = r.truth.override_masses(ra)
    tg = {"region": np.asarray(ov["mass_gdna_region"], np.float64),
          "boundary": np.asarray(ov["mass_gdna_boundary"], np.float64)}
    tr = {"region": np.asarray(ov["mass_rna_region"], np.float64),
          "boundary": np.asarray(ov["mass_rna_boundary"], np.float64)}
    pg = {"region": np.asarray(r.result.mass_gdna_region, np.float64),
          "boundary": np.asarray(r.result.mass_gdna_boundary, np.float64)}
    pr = {"region": np.asarray(r.result.mass_rna_region, np.float64),
          "boundary": np.asarray(r.result.mass_rna_boundary, np.float64)}

    rows = []
    for s in range(int(chain.n_slots)):
        i = int(obj[s])
        axis = "region" if kind[s] == REGION else "boundary"
        if axis == "region":
            label = TYPE_NAMES[int(rtype[i])]
            where = f"{int(start[i]):,}–{int(start[i] + size[i]):,}"
            bp = int(size[i])
        else:
            lo, hi = s - 1, s + 1
            a = int(rtype[obj[lo]]) if lo >= 0 and kind[lo] == REGION else -1
            b = int(rtype[obj[hi]]) if hi < int(chain.n_slots) and kind[hi] == REGION else -1
            pair = "|".join(TYPE_NAMES.get(x, "?") for x in sorted((a, b)) if x >= 0)
            label = pair or "boundary"
            where = f"@{int(start[obj[hi]]):,}" if b >= 0 else f"#{i}"
            bp = 0
        t_tot = tg[axis][i] + tr[axis][i]
        p_tot = pg[axis][i] + pr[axis][i]
        rows.append({
            "slot": s,
            "axis": axis,
            "type": label,
            "where": where,
            "bp": bp,
            "n": float(count[s]),
            "spliced": float(spliced[s]),
            "sj": float(mature[s]),
            "true_fg": float(tg[axis][i] / t_tot) if t_tot > 0 else float("nan"),
            "fg_loc": float(fg_loc[s]),
            "pred_fg": float(fg[s]),
            "sd_fg": float(np.sqrt(max(var_g[s], 0.0))),
            "tau": float(tau[s]),
            "err": float(abs(pg[axis][i] - tg[axis][i])),
            "mass": float(t_tot if t_tot > 0 else p_tot),
        })
    return rows


def report(r: ToyResult) -> None:
    print()
    print("=" * 124)
    print(f"⭐⭐ TOY — {r.spec.name}:  {r.spec.what_it_probes}")
    print("=" * 124)
    print("   " + r.donor.describe().replace("\n", "\n   "))
    print(f"   toy: {r.spec.genome_length:,} bp · {len(r.spec.genes)} gene(s) · "
          f"{r.spec.n_rna_fragments:,} RNA fragments · {r.n_gdna_target:,} gDNA fragments "
          f"(DERIVED to match the donor's density) · nrna={r.spec.nrna_abundance:g}")
    print(f"   {r.seconds:.1f} s end to end")
    # ⭐ THE DRAIN, printed rather than assumed. ``held`` is what pass one could not resolve; ``ambig`` is
    # how many of those the origin lift could not attribute — it BOUNDS the truth error, so a run that
    # does not show it is a run whose oracle carries an unstated error bar (`second_pass.lift_choices`).
    d = getattr(r.payload, "drain", None)
    print(f"   second pass: {0 if d is None else d.offered:,} held · "
          f"{0 if d is None else d.deposited:,} deposited · "
          f"{0 if d is None else d.chose_spliced:,} chose a spliced path · "
          f"oracle lift ambiguous = {r.truth.n_ambiguous:,}")
    rows = object_rows(r)
    print()
    print(f"   {'slot':>4} {'kind':<6} {'type':<20} {'where':>17} {'bp':>7} {'counts':>8} "
          f"{'spl':>6} {'junc':>7} {'true':>6} {'loc':>6} {'pred':>6} {'sd':>7} {'Σ|err|':>9}")
    print("   " + "-" * 118)
    for row in rows:
        if row["mass"] <= 0.0:
            continue
        t = f"{row['true_fg']:.3f}" if np.isfinite(row["true_fg"]) else "  —  "
        print(f"   {row['slot']:>4} {row['axis']:<6} {row['type']:<20} {row['where']:>17} "
              f"{row['bp']:>7,} {row['n']:>8,.0f} {row['spliced']:>6,.0f} {row['sj']:>7,.0f} "
              f"{t:>6} {row['fg_loc']:>6.3f} {row['pred_fg']:>6.3f} {row['sd_fg']:>7.3f} "
              f"{row['err']:>9,.0f}")
    live = [x for x in rows if x["mass"] > 0 and np.isfinite(x["true_fg"])]
    if live:
        w = np.array([x["mass"] for x in live])
        e = np.array([x["err"] for x in live])
        d = np.array([abs(x["pred_fg"] - x["true_fg"]) for x in live])
        print("   " + "-" * 118)
        print(f"   {len(live)} objects with truth · Σ|err| = {e.sum():,.0f} fragments · "
              f"mass-weighted |Δf_g| = {float((d * w).sum() / w.sum()):.4f}")
    print("   ⭐ loc = the message-free local solve; pred = after the relay. If they differ, the")
    print("      messages moved it. sd is the declared sd of log f_g; tau is the own-evidence")
    print("      precision on the log-odds (0 = this object has no answer of its own).")


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE SPEC LADDER — start trivial, add one thing at a time
# ──────────────────────────────────────────────────────────────────────────────────────────────────


def _gene(gid, strand, exons, abundance, t_id=None):
    return {
        "gene_id": gid,
        "strand": strand,
        "transcripts": [{"t_id": t_id or f"{gid}_t1", "exons": exons, "abundance": abundance}],
    }


#: ⭐ Ordered simplest-first, and each one adds exactly ONE structure to the one before it. That is the
#: point: when a row goes wrong, the thing that changed is the thing to look at.
SPECS: dict[str, ToySpec] = {
    # ⚠ There is no gene-FREE rung: `TranscriptIndex` requires at least one transcript, so a
    # chromosome with no annotation cannot be indexed at all. A SILENT gene is the right first rung
    # anyway — it makes every object in the toy structurally pure gDNA, so any deviation from
    # f_g = 1 is a pure false positive with nothing to trade off against.
    "silent": ToySpec(
        name="silent",
        what_it_probes="⭐ ALL objects are pure gDNA (one silent gene) — every deviation from f_g = 1 "
        "is a false positive, with nothing to cancel against it",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 0.0)],
        n_rna_fragments=1,
    ),
    "TA_single_exon": ToySpec(
        name="TA_single_exon",
        what_it_probes="⭐⭐ OWNER'S SPEC. 5 kb chromosome, ONE single-exon transcript TA+ (1000,3000). "
        "NO introns and NO sj, so the exon can ONLY be solved through the two "
        "intergenic|exon BOUNDARIES: intergenic -> boundary -> exon -> boundary -> intergenic. It therefore "
        "tests exactly one thing — can an accurate intergenic gDNA level reach a single-stranded "
        "exon by message passing?",
        genome_length=5_000,
        genes=[_gene("TA", "+", [(1_000, 3_000)], 100.0, t_id="TA")],
        n_rna_fragments=1_000,
    ),
    "one_exon": ToySpec(
        name="one_exon",
        what_it_probes="a single-exon gene: one exon region, two intergenic|exon boundaries, no sj",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000)], 400.0)],
    ),
    "two_exon": ToySpec(
        name="two_exon",
        what_it_probes="ONE intron between two exons — the intron|exon boundaries and the sj flux",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0)],
    ),
    "spliced_exons": ToySpec(
        name="spliced_exons",
        what_it_probes="⭐⭐ OWNER'S SPEC. ONE two-exon transcript TA+ (1,000, 2,000) (9,000, 10,000) "
        "— so this is `nested_exons`'s TWIN at the same gene boundaries on the same 12 kb "
        "chromosome, with an INTRON and a SJ where the nesting was. FIVE REGIONS, FOUR "
        "contiguous BOUNDARIES and ⭐ the ladder's first SJ BOUNDARY:\n"
        "          REGION intergenic [0, 1000)        BOUNDARY @1,000   intergenic|exon, pure gDNA (TSS+)\n"
        "          REGION exon  [1000, 2000)   TA e1  BOUNDARY @2,000   intron|exon, the DONOR+ side\n"
        "          REGION intron [2000, 9000)  TA i1  BOUNDARY @9,000   intron|exon, the ACCEPTOR+ side\n"
        "          REGION exon  [9000, 10000)  TA e2  BOUNDARY @10,000  intergenic|exon, pure gDNA (TES+)\n"
        "          REGION intergenic [10000, 12000)\n"
        "          SJ BOUNDARY 2,000 → 9,000 (+), pure mature RNA, NOT a chain slot\n"
        "        ⭐⭐ What it adds over every rung before it, and why it is the hard one: the two "
        "exon↔intron BOUNDARIES. Mature RNA cannot cross an exon↔intron boundary contiguously, so their truth "
        "is pure gDNA — but the solver's own continuity gate says a strand IS admissible there "
        "(nascent RNA could cross), so they are NOT TRAPS: no-magic-numbers and the solver must *derive* what the structure "
        "already implies. ⛔ On an `nrna_none` donor that is the maximally-violated case of the "
        "intron↔exon imputation premise, so a nascent rung is the control this one needs.\n"
        "        ⛔⛔ CAPTURE-ON NEEDS `--genome-length 120000` ON THIS RUNG. At 12 kb the whole "
        "chromosome gets ~39 gDNA fragments and the two exon↔intron BOUNDARIES carry 2 and 5 — not a solve, an "
        "empty chromosome. At 120 kb they carry 20 and 36 at the EXON's own capture stratum (density "
        "0.079–0.142 against the exon interior's 0.158–0.162 and the intron interior's 0.00015), which is "
        "the regime the object actually matters in. ⭐ That is capture working as intended: the gDNA "
        "signal LEAVES the intergenic and intronic REGIONS and arrives at the BOUNDARIES abutting the exon.\n"
        "        ⭐ And unlike `nested_exons` there IS own evidence inside the gene: the 7,000 bp intron "
        "REGION is where the intron factory lives, so the gDNA level does not have to travel from the "
        "gene ends. The two exons each sit between a G1 gene-boundary BOUNDARY and an exon|intron BOUNDARY, and "
        "the sj's flux is the only measurement of their mature RNA.",
        genome_length=12_000,
        genes=[_gene("g1", "+", [(1_000, 2_000), (9_000, 10_000)], 300.0, t_id="TA")],
    ),
    "alt_splice": ToySpec(
        name="alt_splice",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — ALTERNATIVE SPLICING: several sj meeting "
        "at ONE boundary, which is the case a per-BOUNDARY sj total cannot represent.\n"
        "          TA+ (1,000, 2,000) (5,000, 6,000) (9,000, 10,000)   3 exons — the INCLUSION isoform\n"
        "          TB+ (1,000, 2,000) (9,000, 10,000)                  2 exons — the SKIPPING isoform\n"
        "        THREE sj over TWO shared sites:\n"
        "          j 2,000 -> 5,000   (TA's first intron)\n"
        "          j 6,000 -> 9,000   (TA's second intron)\n"
        "          j 2,000 -> 9,000   (TB's only intron — the exon-skipping jump)\n"
        "        ⭐⭐ SO THE SITES ARE SHARED, AND THAT IS THE POINT: the BOUNDARY @2,000 is the genomic-LOW "
        "end of TWO sj and the BOUNDARY @9,000 is the genomic-HIGH end of TWO. Both of @2,000's fluxes "
        "belong to its LOW flank and both of @9,000's to its HIGH flank, so each bank must POOL them as "
        "`Sum(count)/Sum(E)` — the ratio of sums, never the mean of ratios. ⛔ A single sj-inclusive "
        "total per BOUNDARY cannot express this at all, and neither can a per-sj rule that forgets the "
        "two share a line.\n"
        "        ⭐ It also adds a region that is exon AND intron on the SAME strand: [5,000, 6,000) is TA's "
        "middle exon and lies inside TB's intron. `splice_both_strands` had that contrast only ACROSS "
        "strands; here it is within one, so no strand bit can separate them and `coarse_type_array` calls "
        "it `exon`.\n"
        "        ⚠ What it does NOT cover: an BOUNDARY that is one sj's LOW end and another's HIGH end "
        "at once. That needs one transcript's intron to END where another's BEGINS, and it is gated in "
        "`tests/calibration/test_splice_flux_reframe.py` rather than simulated here.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "gA",
                "strand": "+",
                "transcripts": [
                    {
                        "t_id": "TA",
                        "exons": [(1_000, 2_000), (5_000, 6_000), (9_000, 10_000)],
                        "abundance": 300.0,
                    },
                    {"t_id": "TB", "exons": [(1_000, 2_000), (9_000, 10_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "tes_readthrough": ToySpec(
        name="tes_readthrough",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — the CERTIFIED-RNA CHANNEL AT A TERMINUS BOUNDARY, "
        "which is the case no other rung can produce at all.\n"
        "          TA+ (1,050, 2,000) (9,000,  9,100)\n"
        "          TB+ (1,000, 2,000) (9,050, 11,000)\n"
        "        Two sj from ONE shared donor: j 2,000 -> 9,000 (TA) and j 2,000 -> 9,050 (TB).\n"
        "        ⭐⭐⭐ **BOUNDARY @9,100 IS THE POINT.** It is TA's TES and NO sj touches it — yet "
        "transcription CONTINUES past it, because TB's second exon runs to 11,000. A TB fragment that "
        "USED TB's sj and reaches >50 bp past 9,050 crosses 9,100 **contiguously having spliced "
        "elsewhere**, so it lands in `boundary_spliced` at a line with no sj to price it against. That "
        "is exactly the population the new TSS/TES boundaries create, and if it is not binned as spliced it "
        "falls into the UNSPLICED pool and gets deconvolved — certified RNA fed to the gDNA solver.\n"
        "        ⛔ **No previous toy can make this fragment.** On every earlier rung the exons ARE the "
        "regions, so a spliced molecule never crosses an interior line contiguously and `boundary_spliced` is "
        "structurally zero everywhere (measured 0 on `alt_splice`, including at exons holding 68,000 RNA "
        "fragments).\n"
        "        ⭐ The other three structures, each a separate stress:\n"
        "          BOUNDARY @9,050 — TB's sj ACCEPTOR **and** a plain contiguity line for TA, whose "
        "exon 2 spans 9,000-9,100 unbroken. So one line carries sj flux for one transcript and an "
        "unspliced RNA crossing for another.\n"
        "          BOUNDARY @1,050 — TA's TSS, with TB already transcribing through it.\n"
        "          REGION [9,000, 9,050) — TA exon AND TB intron on the SAME strand, 50 bp wide, so it is "
        "also below one fragment length and has no resolvable density of its own (TRAPS: density-below-one-fragment-length).\n"
        "        ⚠ Both `abundance` values are meant to be SWEPT: the certified channel's strength at "
        "@9,100 is TB's alone, while the unspliced crossing there is gDNA + TB, so the TA/TB ratio moves "
        "the two independently. A single abundance pair tests one corner of that.",
        genome_length=13_000,
        genes=[
            {
                "gene_id": "gA",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_050, 2_000), (9_000, 9_100)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(1_000, 2_000), (9_050, 11_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "splice_both_strands": ToySpec(
        name="splice_both_strands",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — the rung the SPLICE-FLUX REFRAME logic must be "
        "derived against. FOUR transcripts, BOTH strands, overlapping exons AND overlapping introns, "
        "and TWO splice junctions pointing opposite ways:\n"
        "          TA+ (2,000, 3,000) (9,000, 10,000)      2 exons, + strand, intron 3,000-9,000\n"
        "          TB+ (2,000, 10,000)                     1 exon,  + strand, spans TA's intron\n"
        "          TC− (1,000, 11,000)                     1 exon,  − strand, spans everything\n"
        "          TD− (1,000, 2,500) (8,500, 11,000)      2 exons, − strand, intron 2,500-8,500\n"
        "        ⭐⭐ WHY THIS ONE. Every previous rung let an BOUNDARY answer 'is my neighbour an exon?' "
        "with a yes or a no. Here it cannot: THREE regions are simultaneously an INTRON on one strand and "
        "an EXON on the other — [2,500, 3,000), [3,000, 8,500) and [8,500, 9,000) — so 'exon' is not a "
        "property of a region at all, it is a property of (region, strand). And the two sj are on "
        "OPPOSITE strands, so an boundary can be the DONOR of one and sit beside the ACCEPTOR of the other.\n"
        "        ⛔ The question it exists to answer is per (BOUNDARY, side, strand, donor-or-acceptor, "
        "message direction): when this boundary reframes against that neighbour, does its splice flux belong "
        "in the total or not? The derivation is open.\n"
        "        ⚠ MY READING OF THE OWNER'S SPEC, flagged rather than assumed: the owner wrote the last "
        "two transcripts both as `TC-`. Two transcripts cannot share an id, so they are TC− and TD− "
        "here. If the intent was one transcript with two isoforms the ids change and nothing else does.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "gP",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(2_000, 3_000), (9_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 10_000)], "abundance": 300.0},
                ],
            },
            {
                "gene_id": "gM",
                "strand": "-",
                "transcripts": [
                    {"t_id": "TC", "exons": [(1_000, 11_000)], "abundance": 300.0},
                    {"t_id": "TD", "exons": [(1_000, 2_500), (8_500, 11_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "nested_exons": ToySpec(
        name="nested_exons",
        what_it_probes="⭐⭐ NESTED EXONS (owner's spec). THREE nested single-exon transcripts on one gene:\n"
        "                     TA+ (1,000, 10,000)   TB+ (2,000, 9,000)   TC+ (3,000, 8,000)\n"
        "        The partition is therefore SEVEN REGIONS and SIX BOUNDARIES, and there is NO INTRON "
        "anywhere:\n"
        "          REGION intergenic [0, 1000)        BOUNDARY @1,000   intergenic|exon, pure gDNA\n"
        "          REGION exon [1000, 2000)   TA      BOUNDARY @2,000   exon|exon\n"
        "          REGION exon [2000, 3000)   TA+TB   BOUNDARY @3,000   exon|exon\n"
        "          REGION exon [3000, 8000)   TA+B+C  BOUNDARY @8,000   exon|exon\n"
        "          REGION exon [8000, 9000)   TA+TB   BOUNDARY @9,000   exon|exon\n"
        "          REGION exon [9000, 10000)  TA      BOUNDARY @10,000  intergenic|exon, pure gDNA\n"
        "          REGION intergenic [10000, 12000)\n"
        "        ⭐ Every one of the five exon REGIONS carries RNA and the library is unstranded, so NONE "
        "of them has composition evidence of its own; and with no intron there is no object in the "
        "middle of the gene that can re-derive the gDNA level for itself. So the gDNA level has to "
        "travel from the two pure-gDNA BOUNDARIES at the gene ends through five evidence-free objects. "
        "⭐⭐ Because the transcripts are NESTED, the RNA density is a symmetric staircase (1, 2, 3, 2, 1 "
        "transcripts) while the gDNA density is UNIFORM — so the truth is known per object and any "
        "systematic drift in the delivered gDNA level is visible against it.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "g1",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 9_000)], "abundance": 300.0},
                    {"t_id": "TC", "exons": [(3_000, 8_000)], "abundance": 300.0},
                ],
            }
        ],
        n_rna_fragments=2_000,
    ),
    "nested_exons_neg": ToySpec(
        name="nested_exons_neg",
        what_it_probes="⭐ THE STRAND MIRROR of `nested_exons` — the same three nested transcripts on the "
        "MINUS strand. A − transcript runs right-to-left, so its TSS is at its HIGH coordinate and its "
        "TES at its low one: the FLAG_TSS_NEG / FLAG_TES_NEG bits must therefore implicate the OPPOSITE "
        "flank from their + counterparts. ⛔ That is a convention, not a derivation, and a convention "
        "that is assumed rather than pinned is how a sign gets flipped silently — so this rung exists to "
        "pin it against oracle truth.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "g1",
                "strand": "-",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 9_000)], "abundance": 300.0},
                    {"t_id": "TC", "exons": [(3_000, 8_000)], "abundance": 300.0},
                ],
            }
        ],
        n_rna_fragments=2_000,
    ),
    "deep_exon": ToySpec(
        name="deep_exon",
        what_it_probes="⭐ TRAPS: a-purity-filter-is-a-length-filter: a LARGE, RNA-rich exon beside an intron — the relay pins these to ~0.85",
        genome_length=80_000,
        genes=[_gene("g1", "+", [(20_000, 34_000), (40_000, 43_000)], 3000.0)],
        n_rna_fragments=120_000,
    ),
    "nascent": ToySpec(
        name="nascent",
        what_it_probes="two_exon plus NASCENT RNA — the only way an intron|exon boundary legitimately "
        "carries RNA, which the whole ladder panel cannot express",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0)],
        nrna_abundance=40.0,
    ),
    "two_genes": ToySpec(
        name="two_genes",
        what_it_probes="one expressed gene and one SILENT gene — the silent one's exons are pure gDNA",
        genome_length=90_000,
        genes=[
            _gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0),
            _gene("g2", "-", [(50_000, 53_000), (58_000, 61_000)], 0.0),
        ],
    ),
    "opposite_strands": ToySpec(
        name="opposite_strands",
        what_it_probes="two genes on OPPOSITE strands, well separated — the strand channel's control",
        genome_length=90_000,
        genes=[
            _gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0),
            _gene("g2", "-", [(50_000, 53_000), (58_000, 61_000)], 400.0),
        ],
    ),
}


def exon_bp(spec: ToySpec) -> int:
    """Total exonic bases across the spec's transcripts — the denominator for RNA density."""
    return int(sum(e - s for g in spec.genes for t in g["transcripts"] for s, e in t["exons"]))


def sweep_density(
    spec: ToySpec,
    donor: DonorGlobals,
    work_dir: Path,
    densities: list[float],
    *,
    config: CalibrationConfig | None = None,
) -> None:
    """⭐⭐ Vary the transcript's RNA **density** (counts/bp) and nothing else. One variable.

    The gDNA side is untouched by construction: :func:`run_toy` derives the gDNA count from the
    donor's density per base and the chromosome length, neither of which the sweep changes — so the
    gDNA background is **pinned** across every row while the RNA density moves. That is what makes the
    exon's true ``f_g`` sweep from mostly-gDNA to mostly-RNA with a single knob.

    ⚠ The densities are quoted as multiples of the donor's own gDNA density as well as absolutely,
    because what the solver has to resolve is the RATIO, not either level.
    """
    g_rate = donor.gdna_rate_per_base
    ebp = exon_bp(spec)
    print()
    print("=" * 126)
    print(f"⭐⭐ RNA DENSITY SWEEP — {spec.name}   ({spec.genome_length:,} bp, {ebp:,} exonic bp)")
    print("=" * 126)
    print("   " + donor.describe().replace("\n", "\n   "))
    print(f"   gDNA background is PINNED at the donor's {g_rate:.5g} counts/bp for every row;")
    print("   the ONLY thing varied is the transcript's RNA density.")
    print()
    print(f"   {'RNA density':>12} {'vs gDNA':>8} {'RNA':>8} {'gDNA':>6} | "
          f"{'intergenic':^14} | {'boundary@lo':^20} | {'EXON':^20} | {'boundary@hi':^20} | {'toy':>7}")
    print(f"   {'counts/bp':>12} {'x':>8} {'counts':>8} {'counts':>6} | "
          f"{'true':>6} {'pred':>6} | {'true':>6} {'loc':>6} {'pred':>6} | "
          f"{'true':>6} {'loc':>6} {'pred':>6} | {'true':>6} {'loc':>6} {'pred':>6} | {'mwae':>7}")
    print("   " + "-" * 121)
    starved = 0
    for dens in densities:
        n_rna = max(int(round(dens * ebp)), 1)
        one = dataclasses.replace(
            spec, name=f"{spec.name}_d{dens:g}".replace(".", "p"), n_rna_fragments=n_rna
        )
        r = run_toy(one, donor, work_dir, config=config)
        rows = object_rows(r)
        live = [x for x in rows if x["mass"] > 0 and np.isfinite(x["true_fg"])]
        w = np.array([x["mass"] for x in live]) if live else np.array([1.0])
        d = np.array([abs(x["pred_fg"] - x["true_fg"]) for x in live]) if live else np.array([0.0])
        ig = [x for x in rows if x["type"] == "intergenic" and x["mass"] > 0]
        ed = [x for x in rows if x["type"] == "intergenic|exon" and x["mass"] > 0]
        ex = [x for x in rows if x["type"] == "exon" and x["mass"] > 0]

        def cell(objs, keys=("true_fg", "fg_loc", "pred_fg"), pick=0):
            if not objs or pick >= len(objs):
                return " ".join(f"{'—':>6}" for _ in keys)
            o = objs[pick]
            out = []
            for k in keys:
                v = o[k]
                out.append(f"{v:>6.3f}" if np.isfinite(v) else f"{'—':>6}")
            return " ".join(out)

        if not ig or not ed or any(x["mass"] <= 0 for x in ig + ed):
            starved += 1
        print(f"   {dens:>12.5g} {dens / max(g_rate, _EPS):>8.2f} {n_rna:>8,} "
              f"{r.n_gdna_target:>6,} | {cell(ig, ('true_fg', 'pred_fg')):>13} | "
              f"{cell(ed, pick=0)} | {cell(ex)} | {cell(ed, pick=1)} | "
              f"{float((d * w).sum() / w.sum()):>7.4f}")
    print()
    if starved:
        print(f"   ⛔⛔ STARVED on {starved} of {len(densities)} rows: a pure-gDNA object had NO mass, so")
        print("      the chain has nothing to propagate and the exon's answer is meaningless there.")
        print(f"      The donor's gDNA density is {g_rate:.5g} counts/bp and the chromosome is only")
        print(f"      {spec.genome_length:,} bp. ⛔ The density is PINNED to the donor and must NOT be")
        print("      changed. Which lever helps depends on WHICH object is starved:")
        print("      • an intergenic REGION  — lengthen the chromosome (--genome-length): counts scale")
        print("        with bp at fixed density;")
        print("      • an BOUNDARY, capture OFF — ⛔ lengthening does NOTHING. A 0-bp line's counts are")
        print("        `density x mean_FL`, independent of the chromosome, so the only lever is depth.")
        print("      • ⭐⭐ an BOUNDARY, capture ON — LENGTHEN IT. This is the opposite of the capture-OFF")
        print("        case and the reason is the sampler's own mass split: the gDNA budget is")
        print("        `rate x genome_length` while the probe footprint is FIXED, and the on-probe share")
        print("        is `binding x overlap / (off_target x L + binding x overlap)`. So a longer")
        print("        chromosome hands capture a bigger budget to concentrate onto the same probes, and")
        print("        the boundary count grows with L until that ratio saturates. ⭐ Measured on")
        print("        `spliced_exons` x `g75 ss0.50 capture_on`, 12 kb -> 120 kb: the two intron|exon")
        print("        BOUNDARIES go 2 -> 20 and 5 -> 36 counts and the gene-boundary BOUNDARIES 1 -> 41 and")
        print("        2 -> 35, while the intron REGION stays at 1 and the intergenic REGION at ~0 density")
        print("        — i.e. the signal moves to the BOUNDARIES abutting the exon, which is what capture")
        print("        does. ⚠ Raising `binding_per_base` also works but un-matches the toy from the")
        print("        donor's chemistry; lengthening keeps every harvested global intact.")
        print()
    print("   ⭐ THE CHAIN under test:  intergenic → boundary → EXON → boundary → intergenic. With no intron")
    print("      and no sj, the exon's only route to an answer is the two boundaries, so a wrong")
    print("      exon here is a message-passing failure and nothing else.")
    print("   ⭐ loc = the message-free local solve; pred = after the relay. loc ≈ pred means the")
    print("      messages did nothing; loc far from truth with pred near it means they carried it.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="two_exon", help="a spec name, or 'all'")
    ap.add_argument(
        "--sweep-density",
        nargs="*",
        type=float,
        default=None,
        metavar="COUNTS_PER_BP",
        help="sweep the transcript's RNA density (counts/bp) instead of one fixed run; "
        "no values = a default decade ladder around the donor's own gDNA density",
    )
    # ⚠ `g25` until 2026-08-13, retired with the ladder rebuild. `g50` is the surviving mid rung and
    # is also `verify_toy_substrate.py`'s default, so the harness and its verifier agree by default.
    ap.add_argument("--donor", default="gdna_g50_ss_0.50_nrna_none_capture_off")
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_toy")
    ap.add_argument(
        "--genome-length",
        type=int,
        default=None,
        help="override the spec's chromosome length. ⭐ THE LEVER FOR CAPTURE: the gDNA DENSITY is "
        "pinned to the donor and must not change, so the only way to get a well-counted pure-gDNA "
        "background is more intergenic bp. A capture-ON donor's off-target density is ~24x lower "
        "than its capture-OFF twin, so a toy that works at 5 kb off capture needs ~120 kb on it",
    )
    ap.add_argument("--list", action="store_true", help="list the specs and exit")
    args = ap.parse_args()

    if args.list:
        print("specs (simplest first; each adds ONE structure to the one before):\n")
        for name, spec in SPECS.items():
            print(f"  {name:<18} {spec.genome_length:>7,} bp  {len(spec.genes)} gene(s)")
            print(f"  {'':<18} {spec.what_it_probes}")
        return 0

    names = list(SPECS) if args.spec == "all" else [args.spec]
    unknown = [n for n in names if n not in SPECS]
    if unknown:
        print(f"unknown spec(s): {unknown}. --list to see them.", file=sys.stderr)
        return 2

    donor_dir = args.suite / args.donor
    if not (donor_dir / "sim_oracle.bam").is_file():
        print(f"no sim_oracle.bam under {donor_dir}", file=sys.stderr)
        return 2

    index = TranscriptIndex.load(str(args.index))
    print(f"harvesting globals from {args.donor} …", flush=True)
    donor = harvest(donor_dir, index)
    print(donor.describe())

    if args.sweep_density is not None:
        if len(names) != 1:
            print("--sweep-density takes exactly one --spec", file=sys.stderr)
            return 2
        spec = SPECS[names[0]]
        if args.genome_length:
            spec = dataclasses.replace(spec, genome_length=int(args.genome_length))
        # ⭐ the default ladder is anchored on the DONOR's own gDNA density, so the rows are
        # 0.1x .. 1000x the background rather than absolute numbers that mean nothing on their own.
        ladder = args.sweep_density or [
            donor.gdna_rate_per_base * m for m in (0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1000.0)
        ]
        sweep_density(spec, donor, args.work_dir, list(ladder))
        return 0

    for name in names:
        one = SPECS[name]
        if args.genome_length:
            one = dataclasses.replace(one, genome_length=int(args.genome_length))
        report(run_toy(one, donor, args.work_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
