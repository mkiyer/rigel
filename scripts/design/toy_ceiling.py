#!/usr/bin/env python
"""⛔⛔ THE RE-SOLVE CEILING for the `intron|exon` BOUNDARY — NOT a substitution (TRAPS: substitution-understates-a-source).

A substitution replaces one object's answer with the truth and re-scores; that is honest for a SINK and
it UNDERSTATES a message SOURCE, whose whole value is what it carries to its neighbours. So every arm
here hands the BOUNDARY a different OWN BELIEF and then **re-solves the entire chain**, messages and all.

⭐ **The simulation is shared across arms.** Simulate + scan once per (condition, RNA rung), then
re-run only `calibrate`. The arms therefore differ by EXACTLY the solver input under test — same
fragments, same seed, same donor globals — and the extra arms are nearly free.

THE ARMS — a 2x2 on (what the BOUNDARY believes) x (may the BOUNDARY ORIGINATE a gDNA measurement?), plus the
handoff's literal level-transfer variant:

===============  ==========================  ==================  =====================================
arm              the BOUNDARY's own ``f_g``      ``struct_lock``     what the delta from `base` answers
===============  ==========================  ==================  =====================================
``base``         its own self-solve          off                 the baseline, re-recorded here
``intron_phi``   the flanking INTRON's       off                 ⭐ what face (I) as DERIVED can deliver
``intron_rho``   from the INTRON's DENSITY   off                 the LEVEL-transfer variant
``oracle_phi``   ORACLE TRUTH                off                 is the intron's VALUE good enough?
``lock_only``    its own self-solve          ⭐ ON               TRAPS: conservation-misses-mis-attribution: does the BOUNDARY need to ORIGINATE?
``oracle_lock``  ORACLE TRUTH                ⭐ ON               the absolute ceiling for this object
===============  ==========================  ==================  =====================================

⚠ ``struct_lock`` is what admits a slot to the MEASUREMENT stream (``messages.relay``'s ``mg_own``), and
`region_init.strand_evidence` scopes it to REGION slots — so an BOUNDARY can only ever RELAY a gDNA level, never
ORIGINATE one. That is `ISSUES: the-cancelling-pair` stated as code, and ``lock_only`` is its ceiling.

The override reuses `region_init`'s own `own_precision` / `own_composition_logvar`, so there is no second
implementation of the precision arithmetic to drift (TRAPS: two-docstrings-one-quantity).

⭐ **`--arms base noop` IS THE FALSIFICATION.** ``noop`` runs the whole wrapper with an empty target set
and must come back BYTE-IDENTICAL to ``base``; if it does not, every other arm is measuring the rebuild
of the densities rather than the override. Verified on `g50 ss0.50 capture_off`, 63 slots x 12 fields.
⚠ Two perturbations of that gate were run and only one fired: corrupting ``rho_g``'s divisor fires it
(138 values), while dropping the RNA liveness carry-forward is INERT on `spliced_exons` — recorded
rather than assumed, because a gate that cannot fail is not a gate.

⛔⛔ **WHAT IT MEASURED, 2026-08-04, on `spliced_exons` x 36 conditions x 7 rungs x 2 nascent arms.**
Face (I) IS a real mechanism and its ceiling is real but ONE-SIDED — and the ladder then refuted the
build anyway (`ladder_arm_ab.py`), which is why this file exists:

=====================  ==========  =============  ============  ============
stratum                base        intron_phi     oracle_phi    oracle_lock
=====================  ==========  =============  ============  ============
capture off  ss 0.50   0.0955      +0.003         **-0.000**    +0.001
capture off  ss 0.99   0.0803      -0.001         -0.002        +0.000
capture on   ss 0.50   0.2683      -0.027         -0.033        -0.047
capture on   ss 0.99   0.0959      -0.012         -0.012        -0.007
ALL                    0.1322      -0.009         -0.011        -0.012
=====================  ==========  =============  ============  ============

⛔ ``intron_rho`` (transfer the intron's DENSITY rather than its SHARE) is **+0.017 panel-wide and
+0.207 on capture-ON x unstranded** — capture depletes the intron ~1000x while enriching the BOUNDARY, so a
level transfer from it poisons the destination. ⛔ ``lock_only`` (let the BOUNDARY originate with its OWN
value) is **+0.010 / +0.19** — certainty without correctness is worse than no certainty.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str, path: Path):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, path / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


DESIGN = Path(__file__).resolve().parent
TH = _sibling("toy_harness.py", DESIGN)

from rigel.calibration import region_init as NI, sweep as SW  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_geometry import region_gdna_geometry  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer  # noqa: E402
from rigel.scan_cache import index_derived_inputs  # noqa: E402

sys.path.insert(0, str(DESIGN.parents[1] / "tests" / "calibration"))
from _oracle import OracleTruth  # noqa: E402

_EPS = 1.0e-9
#: ⭐ ``noop`` is the falsification arm (the wrapper with nothing to rewrite) and is not run by default.
ARMS = ("base", "noop", "intron_phi", "intron_rho", "oracle_phi", "lock_only", "oracle_lock")
DEFAULT_ARMS = tuple(a for a in ARMS if a != "noop")
LADDER = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
INTRON, EXON = 1, 2


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# simulate + scan ONCE, calibrate MANY times
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclasses.dataclass
class Substrate:
    """Everything downstream of the simulator that every arm shares."""

    spec: object
    donor: object
    payload: object
    strand_model: object
    index: object
    region_arrays: object
    truth: object
    n_gdna_target: int


def simulate(spec, donor, work_dir: Path, pipeline_config=None) -> Substrate:
    """`toy_harness.run_toy`'s first half — everything before `calibrate`."""
    pipeline_config = pipeline_config or PipelineConfig()
    wd = Path(work_dir) / spec.name
    wd.mkdir(parents=True, exist_ok=True)
    n_gdna = int(round(donor.gdna_rate_per_base * spec.genome_length))
    from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario

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
    sc = Scenario(spec.name, genome_length=spec.genome_length, seed=spec.seed, work_dir=wd / "sim")
    for gene in spec.genes:
        sc.add_gene(gene["gene_id"], gene["strand"], gene["transcripts"])
    capture_cfg = None
    if donor.capture_on:
        probes = TH._toy_probes(spec, wd / "probes.tsv", donor.capture_knobs)
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
        gdna_fraction=n_gdna / max(spec.n_rna_fragments, 1),
        nrna_abundance=float(spec.nrna_abundance),
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        capture_config=capture_cfg,
    )
    bam = str(res.bam_path)
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, _buf, payload = scan_and_buffer(bam, res.index, scan)
    truth = OracleTruth.from_bam(
        bam, res.index, pipeline_config, wd / "split", spec.name, full_payload=payload
    )
    return Substrate(
        spec=spec,
        donor=donor,
        payload=payload,
        strand_model=strand_model,
        index=res.index,
        region_arrays=RegionArrays.from_index(res.index),
        truth=truth,
        n_gdna_target=n_gdna,
    )


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# the override
# ──────────────────────────────────────────────────────────────────────────────────────────────────


def _slot_types(chain, region_arrays):
    """Per slot: its coarse type if it is a REGION, and for an BOUNDARY the pair of its two flanking REGIONS."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    n = int(chain.n_slots)
    is_region = kind == REGION
    region_type = np.where(is_region, rtype[np.clip(obj, 0, rtype.shape[0] - 1)], -1)
    left_t = np.full(n, -1, np.int64)
    right_t = np.full(n, -1, np.int64)
    left_t[1:] = np.where(is_region[:-1], region_type[:-1], -1)
    right_t[:-1] = np.where(is_region[1:], region_type[1:], -1)
    return is_region, region_type, left_t, right_t


def _target_slots(chain, region_arrays):
    """The `intron|exon` BOUNDARY slots, each paired with the slot of its flanking INTRON REGION."""
    is_region, region_type, left_t, right_t = _slot_types(chain, region_arrays)
    out = []
    for i in range(int(chain.n_slots)):
        if is_region[i]:
            continue
        pair = {int(left_t[i]), int(right_t[i])}
        if pair == {INTRON, EXON}:
            src = i - 1 if left_t[i] == INTRON else i + 1
            out.append((i, src))
    return out


def _true_fg_per_slot(chain, region_arrays, truth):
    """Oracle ``f_g`` at every slot, from the origin-split masses. NaN where the slot has no mass."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    ov = truth.override_masses(region_arrays)
    g = {"region": np.asarray(ov["mass_gdna_region"], np.float64),
         "boundary": np.asarray(ov["mass_gdna_boundary"], np.float64)}
    r = {"region": np.asarray(ov["mass_rna_region"], np.float64),
         "boundary": np.asarray(ov["mass_rna_boundary"], np.float64)}
    out = np.full(int(chain.n_slots), np.nan)
    for s in range(int(chain.n_slots)):
        ax = "region" if kind[s] == REGION else "boundary"
        i = int(obj[s])
        tot = g[ax][i] + r[ax][i]
        if tot > 0:
            out[s] = g[ax][i] / tot
    return out


def make_override(arm: str, region_arrays, truth):
    """A `build_region_init` wrapper that rewrites the two `intron|exon` BOUNDARIES' OWN belief, then lets the
    whole chain re-solve. Returns ``None`` for the untouched baseline.

    ⭐ ``noop`` runs the ENTIRE wrapper with an empty target set — the falsification arm. It must be
    byte-identical to ``base``; if it is not, the rebuild of the densities and precisions below is not
    faithful to `region_init` and every other arm is measuring the rebuild rather than the override."""
    if arm == "base":
        return None

    def wrapper(chain, statics, geometry, **kw):
        ni = NI.build_region_init(chain, statics, geometry, **kw)
        targets = [] if arm == "noop" else _target_slots(chain, region_arrays)
        f_g = np.array(ni.f_g, np.float64)
        f_pos = np.array(ni.f_pos, np.float64)
        f_neg = np.array(ni.f_neg, np.float64)
        tau = np.array(ni.tau_lam, np.float64)
        lock = np.array(ni.struct_lock, bool)
        M, E_g = region_gdna_geometry(geometry)
        M = np.asarray(M, np.float64)
        E_g = np.asarray(E_g, np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        n_region = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
        true_fg = _true_fg_per_slot(chain, region_arrays, truth)

        for boundary, src in targets:
            if arm in ("intron_phi", "intron_rho"):
                tau[boundary] = tau[src]
            if arm == "intron_phi":
                new_fg = float(f_g[src])
            elif arm == "intron_rho":
                rho_src = float(ni.rho_g[src])
                new_fg = (
                    float(np.clip(rho_src * E_g[boundary] / M[boundary], 0.0, 1.0))
                    if M[boundary] > _EPS
                    else float(f_g[boundary])
                )
            elif arm in ("oracle_phi", "oracle_lock"):
                new_fg = float(f_g[boundary]) if not np.isfinite(true_fg[boundary]) else float(true_fg[boundary])
            else:  # lock_only — the VALUE is untouched, only the certainty changes
                new_fg = float(f_g[boundary])
            if arm in ("lock_only", "oracle_lock"):
                lock[boundary] = True
            # the RNA side follows from the composition: keep the slot's own TILT, rescale to 1 − f_g.
            rna = max(0.0, 1.0 - new_fg)
            tot = f_pos[boundary] + f_neg[boundary]
            if tot > _EPS:
                f_pos[boundary] *= rna / tot
                f_neg[boundary] *= rna / tot
            else:
                fp_ok = bool(np.asarray(statics.free_pos, bool)[boundary])
                fn_ok = bool(np.asarray(statics.free_neg, bool)[boundary])
                k = float(fp_ok) + float(fn_ok)
                f_pos[boundary] = rna * (float(fp_ok) / k) if k else 0.0
                f_neg[boundary] = rna * (float(fn_ok) / k) if k else 0.0
            f_g[boundary] = new_fg

        # ── the densities + precisions, through `region_init`'s OWN arithmetic ──
        v_fg, v_fr = NI.own_composition_logvar(f_g, tau, lock)
        rho_g = np.where((M > _EPS) & (E_g > _EPS), f_g * M / np.maximum(E_g, _EPS), 0.0)
        rho_g = np.maximum(rho_g, 0.0)
        prec_g = NI.own_precision(n_region, v_fg, rho_g > _EPS)

        # ⚠ `region_init._rna` also gates on the local solve's posterior variance being finite, which is
        # not reconstructible from `RegionInit`. Its OUTCOME is: ``rho_old > 0`` exactly where that gate
        # passed — so carry the original liveness forward off the target slots and force it ON at the
        # slots being overridden, which is the whole point of the arm.
        touched = np.zeros(int(chain.n_slots), bool)
        for boundary, _ in targets:
            touched[boundary] = True

        def _rna(frac, free_s, rho_old):
            raw = np.where(
                (M > _EPS) & (E_r > _EPS) & np.asarray(free_s, bool),
                frac * M / np.maximum(E_r, _EPS),
                0.0,
            )
            live = (n_region > 0.0) & (raw > _EPS) & ((rho_old > 0.0) | touched)
            rho = np.where(live, raw, 0.0)
            return rho, NI.own_precision(n_region, v_fr, rho > _EPS)

        rho_pos, prec_pos = _rna(f_pos, statics.free_pos, np.asarray(ni.rho_pos, np.float64))
        rho_neg, prec_neg = _rna(f_neg, statics.free_neg, np.asarray(ni.rho_neg, np.float64))
        return NI.RegionInit(
            f_g=f_g,
            f_pos=f_pos,
            f_neg=f_neg,
            rho_g=rho_g,
            rho_pos=rho_pos,
            rho_neg=rho_neg,
            prec_g=prec_g,
            prec_pos=prec_pos,
            prec_neg=prec_neg,
            struct_lock=lock,
            tau_lam=tau,
        )

    return wrapper


def run_arm(sub: Substrate, arm: str, config) -> dict:
    """Calibrate the shared substrate under one arm and return the per-slot rows."""
    override = make_override(arm, sub.region_arrays, sub.truth)
    debug: dict = {}
    orig = SW.build_region_init
    if override is not None:
        SW.build_region_init = override
    try:
        out = calibrate(
            payload=sub.payload,
            strand_model=sub.strand_model,
            gdna_fl_pmf=sub.donor.gdna_fl_pmf,
            rna_fl_pmf=sub.donor.rna_fl_pmf,
            config=config,
            injected_priors=sub.donor.priors,
            _debug=debug,
            **index_derived_inputs(sub.index),
        )
    finally:
        SW.build_region_init = orig
    r = TH.ToyResult(
        spec=sub.spec,
        donor=sub.donor,
        result=out,
        truth=sub.truth,
        payload=sub.payload,
        region_arrays=sub.region_arrays,
        chain=debug["chain"],
        capture=debug["capture"],
        n_gdna_target=sub.n_gdna_target,
        seconds=0.0,
        index=sub.index,
    )
    rows = TH.object_rows(r)
    # ⭐ the channel state at every slot, so a null result can be ATTRIBUTED rather than guessed:
    # `cm_g` is the gDNA MEASUREMENT precision ψ receives (0 ⇒ the gDNA message is inert — TRAPS: conservation-misses-mis-attribution).
    # ⛔ These four are DECORATION on the ceiling, not the ceiling: `_uni` exists only under
    # `RelayPolicy`, and this file used to die here with `KeyError: '_uni'` under the shipped
    # `message_propagation = False`. Muted, they are NaN — the relay published no claim, and a 0 would
    # read as a measured inert channel (`TRAPS: an-ablation-that-never-ran`).
    uni = TH.relay_channels(debug["capture"])
    for s, row in enumerate(rows):
        for key in ("cm_g", "c_tau", "cg", "mo_g"):
            row[key] = float(uni[key][s]) if uni is not None else float("nan")
    return rows


# ──────────────────────────────────────────────────────────────────────────────────────────────────


def measure(spec_name, conditions, *, suite, index_path, work_dir, refit_iters, genome_length, nrna,
            arms, ladder=LADDER, abs_rna=0, messages=TH.MESSAGES_SHIPPED):
    index = TranscriptIndex.load(str(index_path))
    config = TH.with_messages(
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=int(refit_iters)), messages
    )
    print(TH.messages_stamp(messages), flush=True)
    base = TH.SPECS[spec_name]
    if genome_length:
        base = dataclasses.replace(base, genome_length=int(genome_length))
    if nrna:
        base = dataclasses.replace(base, nrna_abundance=float(nrna))
    ebp = TH.exon_bp(base)
    for cond in conditions:
        t0 = time.perf_counter()
        print(f"  harvesting {cond} …", flush=True)
        donor = TH.harvest(suite / cond, index, config=config)
        g_rate = float(donor.gdna_rate_per_base)
        for mult in ladder:
            # ⚠ ``abs_rna`` overrides the multiple with an ABSOLUTE count, and it is not a convenience:
            # on a ZERO-gDNA condition the donor's rate is 0, so ``m × rate × bp`` rounds to 1 fragment
            # and the row is degenerate rather than easy. A gDNA-free library is the sharpest
            # false-positive test there is and it needs a real RNA depth to be one.
            n_rna = int(abs_rna) if abs_rna else max(int(round(mult * g_rate * ebp)), 1)
            spec = dataclasses.replace(
                base,
                name=f"{base.name}_m{mult:g}".replace(".", "p"),
                n_rna_fragments=n_rna,
            )
            sub = simulate(spec, donor, work_dir)
            for arm in arms:
                for row in run_arm(sub, arm, config):
                    yield {
                        "condition": cond,
                        "arm": arm,
                        # ⭐ STAMPED INTO EVERY ROW, the `ladder_arm_ab.py` contract: `--out` rows outlive
                        # the run, and `cm_g`/`c_tau`/`cg`/`mo_g` are NaN under the muted relay — a row
                        # that does not carry its own message setting cannot be read later.
                        "messages": "on" if messages else "off",
                        "capture": "on" if donor.capture_on else "off",
                        "strand": f"{donor.strand_specificity:g}",
                        "gdna_rate": g_rate,
                        "mult": mult,
                        "n_rna": n_rna,
                        "nrna": float(base.nrna_abundance),
                        **{k: row[k] for k in (
                            "slot", "axis", "type", "where", "bp", "n", "spliced", "sj",
                            "true_fg", "fg_loc", "pred_fg", "sd_fg", "tau", "err", "mass",
                            "cm_g", "c_tau", "cg", "mo_g",
                        )},
                    }
        print(f"    {cond} done in {time.perf_counter() - t0:.0f} s", flush=True)


def _key(r):
    return f"{r['type']}@{r['where']}"


def _live(rows):
    return [r for r in rows if r["mass"] > 0 and np.isfinite(r["true_fg"])]


def _mwae(rows):
    live = _live(rows)
    if not live:
        return float("nan")
    w = np.array([r["mass"] for r in live])
    d = np.array([abs(r["pred_fg"] - r["true_fg"]) for r in live])
    return float((d * w).sum() / w.sum())


def report(rows, spec_name, refit_iters):
    arms = [a for a in ARMS if any(r["arm"] == a for r in rows)]
    conds = sorted({r["condition"] for r in rows})
    mults = sorted({r["mult"] for r in rows})
    strata = sorted({(r["capture"], r["strand"]) for r in rows})
    nrna = sorted({r["nrna"] for r in rows})
    print()
    print("=" * 136)
    print(f"⛔⛔ RE-SOLVE CEILING — {spec_name}, {len(conds)} conditions x {len(mults)} rungs x "
          f"{len(arms)} arms, refit={refit_iters}, nrna={nrna}")
    print("=" * 136)
    print("   Every arm RE-SOLVES the whole chain on the SAME simulated fragments; the arms differ only")
    print("   in what the two `intron|exon` BOUNDARIES believe about themselves. ⛔ Not a substitution (TRAPS: substitution-understates-a-source).")

    def _cells(fmt, cells, width=15):
        """arm 0 plain; every other arm with its DELTA from arm 0 beside it."""
        out = [f"{format(cells[0], fmt):>{width}}"]
        for c in cells[1:]:
            d = c - cells[0]
            out.append(f"{format(c, fmt) + f' ({d:+.3f})':>{width}}")
        return "".join(out)

    print()
    print("── 1. THE GENE'S mass-weighted |Δf_g|, per arm x stratum ──────────────────────────────────")
    print(f"\n   {'stratum':<22}" + "".join(f"{a:>15}" for a in arms))
    print("   " + "-" * (22 + 15 * len(arms)))

    def _gene_mwae(sel, a):
        per = [_mwae([r for r in sel if r["arm"] == a and r["condition"] == c and r["mult"] == m])
               for c in conds for m in mults]
        per = [v for v in per if np.isfinite(v)]
        return float(np.mean(per)) if per else float("nan")

    for cap, ss in strata:
        sel = [r for r in rows if r["capture"] == cap and r["strand"] == ss]
        print(f"   capture {cap:<3} ss {ss:<7}" + _cells(".4f", [_gene_mwae(sel, a) for a in arms]))
    print(f"   {'ALL':<22}" + _cells(".4f", [_gene_mwae(rows, a) for a in arms]))

    print()
    print("── 2. PER OBJECT — mean |Δf_g| per arm ─────────────────────────────────────────────────────")
    by = defaultdict(list)
    for r in rows:
        by[_key(r)].append(r)
    order = sorted(by, key=lambda k: -sum(r["err"] for r in by[k] if r["arm"] == "base"))
    print(f"\n   {'object':<26} {'n':>6}" + "".join(f"{a:>15}" for a in arms))
    print("   " + "-" * (33 + 15 * len(arms)))
    for k in order:
        g = _live(by[k])
        if not g:
            continue
        cells = []
        for a in arms:
            gg = [r for r in g if r["arm"] == a]
            cells.append(float(np.mean([abs(r["pred_fg"] - r["true_fg"]) for r in gg]))
                         if gg else float("nan"))
        print(f"   {k:<26} {np.mean([r['n'] for r in g]):>6.0f}" + _cells(".4f", cells))

    print()
    print("── 3. ⭐ THE ERROR MASS each arm removes (Σ|Δ gDNA mass|, the deliverable's own currency) ──")
    print(f"\n   {'object':<26}" + "".join(f"{a:>17}" for a in arms))
    print("   " + "-" * (26 + 17 * len(arms)))
    tot = {a: 0.0 for a in arms}
    for k in order:
        cells = [sum(r["err"] for r in by[k] if r["arm"] == a) for a in arms]
        for a, c in zip(arms, cells):
            tot[a] += c
        line = f"{cells[0]:>17,.0f}" + "".join(
            f"{f'{c:,.0f} ({(c - cells[0]) / max(cells[0], 1e-9):+.0%})':>17}" for c in cells[1:])
        print(f"   {k:<26}" + line)
    t0 = tot[arms[0]]
    print(f"   {'TOTAL':<26}" + f"{t0:>17,.0f}" + "".join(
        f"{f'{tot[a]:,.0f} ({(tot[a] - t0) / max(t0, 1e-9):+.0%})':>17}" for a in arms[1:]))

    print()
    print("── 4. ⭐⭐ THE CHANNEL: is the BOUNDARY→exon gDNA message ALIVE? (`cm_g` at the EXON slots) ────")
    print("   `cm_g` is the gDNA MEASUREMENT precision ψ receives. 0 ⇒ the level is carried only as a")
    print("   mode with no weight — `ISSUES: the-cancelling-pair`, a G1 BOUNDARY cannot ORIGINATE.")
    exon_keys = [k for k in by if k.startswith("exon@")]
    # ⛔ THE MUTED RELAY MUST NOT BE READ AS A DEAD CHANNEL. Rows carry their own `messages` stamp, and
    # a muted row has NaN in all four channel columns — so `% slots cm_g=0` would report 0 % (nothing is
    # `<= 0`) and read as "the channel is alive" on a run where no channel was ever asked.
    if any(r.get("messages", "on") == "off" for rs in by.values() for r in rs):
        print(TH.relay_silent_note("SECTION 4 (cm_g / c_tau at the exon slots)"))
    print(f"\n   {'object':<26} {'arm':<14} {'cm_g':>12} {'c_tau':>12} {'% slots cm_g=0':>16}")
    for k in exon_keys:
        for a in arms:
            g = [r for r in by[k] if r["arm"] == a]
            if not g:
                continue
            cm = np.array([r["cm_g"] for r in g])
            asked = np.isfinite(cm)
            zero = f"{(cm[asked] <= 0).mean():>15.0%}" if asked.any() else f"{'— not asked':>15}"
            print(f"   {k:<26} {a:<14} {cm.mean():>12.4g} "
                  f"{np.mean([r['c_tau'] for r in g]):>12.4g} {zero}")

    print()
    print("── 5. THE SWEEP, on the objects that carry the error ───────────────────────────────────────")
    for k in order[:4]:
        print(f"\n   {k}")
        print(f"      {'rung':>7} {'true':>7}" + "".join(f"{a:>13}" for a in arms))
        for m in mults:
            g = _live([r for r in by[k] if r["mult"] == m])
            if not g:
                continue
            cells = []
            for a in arms:
                gg = [r for r in g if r["arm"] == a]
                cells.append(np.mean([r["pred_fg"] for r in gg]) if gg else np.nan)
            print(f"      {'m=' + f'{m:g}':>7} {np.mean([r['true_fg'] for r in g]):>7.4f}"
                  + "".join(f"{c:>13.4f}" for c in cells))


def main() -> int:
    P0 = _sibling("pass0_vs_oracle.py", DESIGN)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="spliced_exons")
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE.parent / "ladder")
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--arms", nargs="*", default=list(DEFAULT_ARMS))
    ap.add_argument("--refit-iters", type=int, default=0)
    ap.add_argument("--genome-length", type=int, default=0)
    ap.add_argument("--nrna", type=float, default=0.0)
    ap.add_argument("--rungs", default=None,
                    help="comma-separated RNA multiples, overriding the default ladder. ⛔ `m` is a "
                         "multiple of the donor's OFF-TARGET gDNA density, and capture then "
                         "concentrates gDNA onto the exon ~65x — so the SAME m is a completely "
                         "different true composition on the two capture arms and they need "
                         "different rungs to reach the same operating point")
    ap.add_argument("--abs-rna", type=int, default=0,
                    help="an ABSOLUTE RNA fragment count, overriding the multiple. Needed on a "
                         "zero-gDNA condition, where the multiple degenerates to 1 fragment")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_toy_ceiling"))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--report", type=Path, default=None)
    # ⭐ the SHIPPED setting by default: the ceiling — hand one object class a different own belief and
    # re-solve — is policy-independent, and the four `_uni` columns are attribution decoration on it.
    TH.add_messages_flag(ap, default=TH.MESSAGES_SHIPPED)
    args = ap.parse_args()

    if args.report:
        rows = [json.loads(x) for x in args.report.read_text().splitlines() if x.strip()]
        report(rows, args.spec, args.refit_iters)
        return 0

    conds = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    rows = []
    fh = args.out.open("w") if args.out else None
    try:
        for row in measure(args.spec, conds, suite=args.suite, index_path=args.index,
                           work_dir=args.work_dir, refit_iters=args.refit_iters,
                           genome_length=args.genome_length, nrna=args.nrna, arms=args.arms,
                           ladder=(tuple(float(x) for x in args.rungs.split(","))
                                   if args.rungs else LADDER),
                           abs_rna=args.abs_rna, messages=TH.messages_on(args)):
            rows.append(row)
            if fh:
                fh.write(json.dumps(row) + "\n")
                fh.flush()
    finally:
        if fh:
            fh.close()
    report(rows, args.spec, args.refit_iters)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
