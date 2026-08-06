#!/usr/bin/env python
"""⭐⭐⭐ WHAT IS PERFECTING THE SIMPLEX VERTEX WORTH? — the re-solve ceiling, on the REAL 36-condition
ladder, plus the mechanism prototype in the same harness so the two are directly comparable.

⛔⛔ **THIS IS THE MEASUREMENT THAT DECIDES BUILD-VS-NOTE** (`TRAPS.md` B1, which has re-ranked this
project five times). A silent gene's objects are pure gDNA and the truth is ``f_g = 1.000`` exactly; a
zero-gDNA library's objects are pure RNA and the truth is ``0.000`` exactly. The solver lands 2–8 % inside
both. Before designing anything, hand those objects the exact answer, **re-solve the whole chain**, and
read what the headroom actually is.

⭐ **A RE-SOLVE, NOT A SUBSTITUTION** (`TRAPS.md` B17). A vertex object is mostly a message SOURCE — its
value is what it CARRIES — and substituting its answer after the fact does not propagate. So the pin goes
in at ``node_init.build_node_init``, the per-object message-free self-solve, and the relay then runs on
top of it. That is also the one place a mechanism can be expressed without touching either relay twin,
so a prototype cannot be gated in one twin and not the other (`TRAPS.md` B14).

**THE ARMS.**

| arm | what it does | what it is |
|---|---|---|
| ``base`` | nothing | the baseline, re-recorded from the current tree in the same run (`TRAPS.md` B8) |
| ``noop`` | pins the truth at ZERO objects | ⭐ the harness's own falsification — MUST be byte-identical to ``base`` |
| ``vertex_free`` | pins oracle truth at every vertex-truth object with **no own composition evidence** | ⭐⭐ **THE CEILING.** The population a vertex fix can reach |
| ``vertex_all`` | pins oracle truth at **every** vertex-truth object | a looser upper bound — includes objects that have their own evidence |
| ``ref_c=<C>`` | sets ψ's reference exponent to ``C`` instead of ½ | ⭐ the mechanism prototype (`TRAPS.md` B18 — the panel arm before ``src/``) |

⚠ ``vertex_free``'s "no own evidence" test is ``tau_lam <= 1e-4``, and that is a CLASSIFICATION FOR A
CEILING, never a production predicate — `TRAPS.md` B11 refused exactly this shape as a mechanism. Both
arms are reported side by side so the filter's effect is visible rather than assumed.

⛔ **A12 — the honesty columns are never quoted alone.** Every row carries ``mwae_all`` (denominator =
every object with mass) and ``abs_err_all`` (no denominator at all) beside the solvable-set numbers,
because an arm that changes what counts as solvable changes its own denominator.

⛔ **A10 — every arm counts its own firings and RAISES if it did not fire.** An override that never ran
reads as "no effect", which is publishable and wrong.

⚠ **A14 — the pin count per condition is printed.** An arm with zero opportunities is not a control.

---

⛔⛔⛔ **WHAT IT MEASURED, 2026-08-05 — and the verdict is NOT A BUILD. READ THIS FIRST.**

The number below is real and reproducible, and it is **the value of INFORMATION, not headroom for a
fix.** The objects it pins are HONEST: measured per-object, `|f_g - truth| / sd(f_g)` has median
**z = 0.5-0.6** on every arm and both simplex vertices, i.e. every wrong answer sits inside its own 1
sigma with a variance that is if anything conservative. And no prior-free solve can do better: every
PROPER prior on [0,1] has a median strictly inside (0,1), an object with zero composition information
has posterior = prior, so a vertex is unreachable there in ANY coordinate at ANY depth.
⛔ Pass-0 must stay prior-free — its purpose is to produce the substrate a prior is fitted ON — so
'fit a prior to fix this' is circular. ⭐ Use this number to SIZE the cost of missing information, and
look for the pass-0 defect in the **confidently-wrong** population instead, which is a different set of
objects (`SESSION_HANDOFF.md`).

``vertex_free``, against a ``base`` re-recorded in the same run, with ``noop`` byte-identical on all 36
rows of both axes (the harness's own falsification passing):

======================================  =========  ===========  ==================
the deliverable — library ``f_g``       base       vertex_free
======================================  =========  ===========  ==================
mean \\|error\\| at pass-0                0.1036     **0.0804**   −22.4 %
mean \\|error\\| on the SHIPPED solve      0.0538     **0.0407**   ⭐ **−24.4 %**
======================================  =========  ===========  ==================

⭐⭐ **For scale: perfecting BOTH fragment-length models is worth 2.6 % of the same deliverable.** The
vertex is ~9× the entire Stage-A length ceiling.

Per-object, pass-0 ``mwae`` over ALL objects (fixed denominator; ``solv%`` is byte-identical across every
arm, so none of this is a denominator move):

=========  ======  ============  ==============  ==============
axis       base    vertex_free   Σ\\|err\\| frags   better/worse
=========  ======  ============  ==============  ==============
node       0.1247  **0.0975**    −149,267        27 / 9
edge       0.1434  **0.1127**    −161,302        29 / 3
=========  ======  ============  ==============  ==============

⭐⭐⭐ **AND IT SPLITS ON EXACTLY ONE AXIS — STRAND — which is what the mechanism predicts.** Pass-0
``mwae`` delta, node axis: unstranded **−0.0188** (capture off, 9/0) and **−0.0963** (capture ON, 9/0);
stranded **−0.0003** and **+0.0064** (2/7). Every one of the 9 rows that got worse is ``ss_0.99``.
The strand channel's Fisher information is ``∝ (2κ−1)²`` and is EXACTLY zero at κ = ½, so on an
unstranded library ψ's reference is the only term left with a gradient at the vertex, while on a stranded
one the strand term supplies the λ information and the vertex is already reached. ⭐ The ceiling is
therefore entirely on unstranded data — which is also the panel's worst stratum
(``capture_ON × ss0.50``: base 0.3235 node / 0.2922 edge). Largest single row:
``gdna_g01_ss_0.50_capture_on``, **−0.2188**.

⛔⛔ **TWO WARNINGS THAT MUST TRAVEL WITH THE NUMBER.**

* **``vertex_all`` is WORSE than ``vertex_free``** — node pass-0 0.1076 vs 0.0975, and on the shipped
  solve it is worse than *base* (+0.0080). Pinning MORE truth hurts: the extra objects have their own
  evidence, and declaring them certain overrides it and propagates. That is `TRAPS.md` B20's shape
  reached with the TRUTH, so the harm is in the relay's dynamics and not in the values. ⭐ Quote
  ``vertex_free``, and note that a fix which hands out certainty broadly can lose even when it is right.
* **The honesty columns move the WRONG way** — confidently-wrong Σ\\|err\\| +9,175 (node) / +893 (edge),
  28 and 16 rows worse — while accuracy improves 22 %. `TRAPS.md` A12 exactly: certainty handed to an
  object moves it into the confident population. Read ``mwae_all`` and ``abs_err_all``, never these.

Usage::

    # one condition first, to check the levers are connected
    python scripts/design/vertex_ceiling.py --arm base       --conditions gdna_g50_ss_0.50_nrna_none_capture_off --out /tmp/v_base.jsonl
    python scripts/design/vertex_ceiling.py --arm vertex_free --conditions gdna_g50_ss_0.50_nrna_none_capture_off --out /tmp/v_free.jsonl
    # the whole ladder, with the oracle cache
    python scripts/design/vertex_ceiling.py --arm vertex_free \
        --oracle-cache ~/Downloads/rigel_runs/suite/ladder/oracle_cache --out /tmp/v_free.jsonl
    # and the comparison
    python scripts/design/vertex_ceiling.py --compare /tmp/v_base.jsonl /tmp/v_free.jsonl
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
sys.path.insert(0, str(DESIGN))


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, DESIGN / name)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[key] = mod
        spec.loader.exec_module(mod)
    return sys.modules[key]


SA = _sibling("solvability_audit.py")
P0 = _sibling("pass0_vs_oracle.py")

from rigel.calibration import bp_solver, node_init as NI  # noqa: E402
from rigel.calibration import simplex_logodds as SL  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.calibration.node_geometry import node_global_geometry  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

CAL = sys.modules["rigel.calibration.calibrate"]

_EPS = 1.0e-9
#: the ceiling's own classification of "this object has no composition evidence of its own". ⛔ NOT a
#: production predicate (`TRAPS.md` B11); the `vertex_all` arm exists so its effect is measured, not
#: assumed.
_TAU_FREE = 1.0e-4

#: filled by the wrappers, one call before `build_node_init` needs them.
_CTX: dict = {}
#: A10 — per-arm firing counters. A zero here RAISES.
_FIRED: dict = {"init": 0, "pinned": 0, "ref_g": 0, "ref_r": 0, "conditions": 0}


# ── the plumbing: get the oracle's per-object truth and the geometry to `build_node_init` ────────────


def _wrap_oracle():
    """Stash the origin-split oracle. It is built BEFORE any arm calibrates, so its per-object truth
    is available to the self-solve — which is what makes a re-solve ceiling possible at all."""
    orig = P0.load_or_build_oracle

    def wrapper(*a, **kw):
        oracle = orig(*a, **kw)
        _CTX["oracle"] = oracle
        return oracle

    P0.load_or_build_oracle = wrapper


def _wrap_node_sweep():
    """Stash `region_arrays` — `node_sweep` receives it and calls `build_node_init` after."""
    orig = CAL.node_sweep

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _CTX["region_arrays"] = region_arrays
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.node_sweep = wrapper


def _truth_fg_per_slot(chain):
    """The ORACLE's true ``f_g`` per SLOT, and the mass behind it.

    ⚠ ``NodeInit``'s arrays are indexed by SLOT (the chain's alternating NODE/EDGE sequence), while the
    oracle's masses are per OBJECT on two separate axes — so the mapping goes through
    ``chain.kind``/``chain.obj_idx`` rather than being assumed."""
    oracle = _CTX.get("oracle")
    ra = _CTX.get("region_arrays")
    if oracle is None or ra is None:
        return None, None
    ov = oracle.override_masses(ra)
    g = {
        "node": np.asarray(ov["mass_gdna_node"], np.float64),
        "edge": np.asarray(ov["mass_gdna_edge"], np.float64),
    }
    r = {
        "node": np.asarray(ov["mass_rna_node"], np.float64),
        "edge": np.asarray(ov["mass_rna_edge"], np.float64),
    }
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_node = kind == NODE
    n = int(chain.n_slots)
    tg = np.zeros(n)
    tr = np.zeros(n)
    for axis, msk in (("node", is_node), ("edge", ~is_node)):
        idx = np.flatnonzero(msk)
        if idx.size == 0:
            continue
        o = np.clip(obj[idx], 0, g[axis].shape[0] - 1)
        tg[idx] = g[axis][o]
        tr[idx] = r[axis][o]
    tot = tg + tr
    with np.errstate(invalid="ignore", divide="ignore"):
        fg = np.where(tot > 0.0, tg / np.maximum(tot, _EPS), np.nan)
    return fg, tot


def _install_vertex_pin(evidence_free_only: bool, force_empty: bool = False):
    """⭐⭐ THE CEILING ARM. Overwrite ``f_g`` with the ORACLE's exact answer at every object whose truth
    sits on a **vertex** of the composition simplex (``f_g`` exactly 0 or exactly 1), declare that belief
    certain, and let the relay re-solve on top of it.

    ⭐ Only the VERTEX population is pinned. An interior object keeps its own answer, so this prices the
    vertex and nothing else — which is the whole point of a channel ceiling.

    ⚠ ``evidence_free_only`` restricts the pin to objects with no own composition evidence
    (``tau_lam <= _TAU_FREE``). That is the population a vertex fix can actually reach; the unrestricted
    arm is the looser bound."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        _FIRED["init"] += 1
        true_fg, true_mass = _truth_fg_per_slot(chain)
        if true_fg is None:
            return ni
        f_g = np.array(ni.f_g, np.float64)
        f_pos = np.array(ni.f_pos, np.float64)
        f_neg = np.array(ni.f_neg, np.float64)
        tau = np.array(ni.tau_lam, np.float64)
        lock = np.array(ni.struct_lock, bool)

        at_vertex = np.isfinite(true_fg) & (true_mass > 0.0) & (
            (true_fg <= 0.0) | (true_fg >= 1.0)
        )
        if evidence_free_only:
            at_vertex &= tau <= _TAU_FREE
        if force_empty:
            # ⭐ the `noop` arm: the WHOLE wrapper runs — the oracle is read, the truth is mapped to
            #   slots, the classification is evaluated — and then nothing is pinned. Byte-identical to
            #   `base` is the assertion; anything else means the wrapper itself moves the answer (A5).
            at_vertex[:] = False
        tgt = np.flatnonzero(at_vertex)
        _FIRED["pinned"] += int(tgt.size)
        if tgt.size == 0:
            return ni

        # ── the pin: the exact composition, and the RNA half split across whichever strands are live.
        #    ⭐ A vertex is the one place this needs no share model: at f_g = 1 there is no RNA to
        #    split, and at f_g = 0 the split is the object's own strand freedom. ──
        new_fg = true_fg[tgt]
        rna = np.maximum(0.0, 1.0 - new_fg)
        fp_ok = np.asarray(statics.free_pos, bool)[tgt]
        fn_ok = np.asarray(statics.free_neg, bool)[tgt]
        tot = f_pos[tgt] + f_neg[tgt]
        k = fp_ok.astype(np.float64) + fn_ok.astype(np.float64)
        share_p = np.where(
            tot > _EPS, f_pos[tgt] / np.maximum(tot, _EPS), np.where(k > 0, fp_ok / np.maximum(k, 1.0), 0.0)
        )
        share_n = np.where(
            tot > _EPS, f_neg[tgt] / np.maximum(tot, _EPS), np.where(k > 0, fn_ok / np.maximum(k, 1.0), 0.0)
        )
        f_g[tgt] = new_fg
        f_pos[tgt] = rna * share_p
        f_neg[tgt] = rna * share_n
        # ⭐ CERTAIN, the same way a structurally pure-gDNA object is certain — via `struct_lock`, which
        #   `own_composition_logvar` already reads. A ceiling must hand over the answer AND the
        #   confidence, or the relay simply argues it back (`TRAPS.md` B17).
        lock[tgt] = True

        v_fg, v_fr = NI.own_composition_logvar(f_g, tau, lock)
        M, E_g = node_global_geometry(geometry)
        M = np.asarray(M, np.float64)
        E_g = np.asarray(E_g, np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        n_node = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
        rho_g = np.maximum(
            np.where((M > _EPS) & (E_g > _EPS), f_g * M / np.maximum(E_g, _EPS), 0.0), 0.0
        )
        prec_g = NI.own_precision(n_node, v_fg, rho_g > _EPS)
        touched = np.zeros(f_g.shape[0], bool)
        touched[tgt] = True

        def _rna(frac, free_s, rho_old):
            raw = np.where(
                (M > _EPS) & (E_r > _EPS) & np.asarray(free_s, bool),
                frac * M / np.maximum(E_r, _EPS),
                0.0,
            )
            live = (n_node > 0.0) & (raw > _EPS) & ((rho_old > 0.0) | touched)
            rho = np.where(live, raw, 0.0)
            return rho, NI.own_precision(n_node, v_fr, rho > _EPS)

        rho_pos, prec_pos = _rna(f_pos, statics.free_pos, np.asarray(ni.rho_pos, np.float64))
        rho_neg, prec_neg = _rna(f_neg, statics.free_neg, np.asarray(ni.rho_neg, np.float64))
        return NI.NodeInit(
            f_g=f_g, f_pos=f_pos, f_neg=f_neg,
            rho_g=rho_g, rho_pos=rho_pos, rho_neg=rho_neg,
            prec_g=prec_g, prec_pos=prec_pos, prec_neg=prec_neg,
            struct_lock=lock, tau_lam=tau,
        )

    bp_solver.build_node_init = wrapper


def _install_ref_exponent(c_value: float):
    """⭐ THE MECHANISM PROTOTYPE — ψ's reference exponent as a free number instead of ½.

    ψ writes ONE constant twice with opposite signs (``simplex_logodds._JEFFREYS_REF``):
    ``+C·log f_g`` bounds ``f_g → 0`` and ``+C·log(1−f_g)`` bounds ``f_g → 1``. Those two terms are the
    ONLY things holding the posterior median off the two vertices.

    ⛔ ``C = 0`` makes ψ improper on both sides (Beta(0,0) — `TRAPS.md` D5's Haldane, a vertex
    amplifier), so a small ``C`` is a PROTOTYPE that bounds what a derived fix could buy, never the fix
    itself. The point of putting it in this harness is that it is priced on the same panel and by the
    same instrument as the ceiling above."""
    def _gdna_arm(lam, global_logprior):
        _FIRED["ref_g"] += 1
        ref = float(c_value) * SL._log_fg(lam)[None, :]
        if global_logprior is None:
            return ref
        return ref + np.asarray(global_logprior, np.float64)

    def _rna_arm(lam):
        _FIRED["ref_r"] += 1
        return float(c_value) * SL._log1m_fg(lam)[None, :]

    SL._gdna_arm = _gdna_arm
    SL._rna_arm = _rna_arm


# ── the comparison ──────────────────────────────────────────────────────────────────────────────────


def _compare(paths: list[Path]) -> int:
    """Read two or more arm files and print the per-axis deltas, with the fixed-denominator columns
    first because those are the ones that cannot be gamed by knowing less (`TRAPS.md` A12)."""
    arms: dict[str, dict] = {}
    for p in paths:
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            arms.setdefault(row["arm"], {})[(row["condition"], row["axis"])] = row
    names = list(arms)
    if len(names) < 2:
        print(f"⛔ need >= 2 arms, got {names}")
        return 1
    base = names[0]
    # ⛔ A12 ORDER: the two FIXED-DENOMINATOR columns come FIRST, because they are the only two that
    #    cannot be gamed by the solver knowing less. `solvable_mwae` and the honesty columns follow.
    cols = [
        ("mwae_all", "mwae ALL", "lower"),
        ("abs_err_all", "Σ|err| ALL", "lower"),
        ("mwae_all_final", "mwae ALL (final)", "lower"),
        ("solvable_mwae", "mwae solvable", "lower"),
        ("solvable_mass_share", "solv% (mass)", "context"),
        ("conf_wrong_err", "confidently wrong", "lower"),
        ("conf_wrong_objects", "conf-wrong objects", "lower"),
        ("library_f_gdna_pass0", "library f_g pass0", "context"),
    ]
    for axis in ("node", "edge"):
        print(f"\n{'=' * 118}\n⭐ AXIS = {axis}\n{'=' * 118}")
        print(f"   {'metric':<22}{'arm':<16}{'mean':>12}{'vs base':>12}{'better':>9}"
              f"{'worse':>7}{'flat':>6}   rows")
        print("   " + "-" * 112)
        for key, label, _ in cols:
            bvals = {
                c: r.get(key) for (c, a), r in arms[base].items() if a == axis and r.get(key) is not None
            }
            if not bvals:
                continue
            print(f"   {label:<22}{base:<16}{np.mean(list(bvals.values())):>12.4f}"
                  f"{'—':>12}{'—':>9}{'—':>7}{'—':>6}   {len(bvals)}")
            for nm in names[1:]:
                vals = {
                    c: r.get(key)
                    for (c, a), r in arms[nm].items()
                    if a == axis and r.get(key) is not None
                }
                shared = sorted(set(vals) & set(bvals))
                if not shared:
                    continue
                b = np.array([bvals[c] for c in shared], float)
                v = np.array([vals[c] for c in shared], float)
                better = int(np.sum(v < b - 1e-12))
                worse = int(np.sum(v > b + 1e-12))
                flat = len(shared) - better - worse
                print(f"   {'':<22}{nm:<16}{v.mean():>12.4f}{v.mean() - b.mean():>+12.4f}"
                      f"{better:>9}{worse:>7}{flat:>6}   {len(shared)}")
        # ⛔ A byte-identical arm is NOT evidence of no change (`TRAPS.md` B4/A5) — EXCEPT for `noop`,
        #   where it is the assertion the arm exists to make. Say which of the two it is.
        for nm in names[1:]:
            shared = [c for (c, a) in arms[nm] if a == axis and (c, a) in arms[base]]
            same = sum(
                1
                for c in shared
                if abs(arms[nm][(c, axis)].get("mwae_all", 0.0)
                       - arms[base][(c, axis)].get("mwae_all", 0.0)) < 1e-12
            )
            if not shared or same != len(shared):
                continue
            if nm == "noop":
                print(f"   ✅ {nm} is byte-identical to {base} on all {len(shared)} rows — the harness's"
                      f" own falsification PASSES: the wrapper does not move the answer by itself.")
            else:
                print(f"   ⚠ {nm} is byte-identical to {base} on all {len(shared)} rows of this axis"
                      f" — if that arm was meant to CHANGE something, it did not fire (A10/A5).")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default=None,
                    help="base | noop | vertex_free | vertex_all | ref_c=<float>")
    ap.add_argument("--compare", nargs="*", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_vertex_ceiling"))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if args.compare:
        return _compare(args.compare)
    if not args.arm or not args.out:
        ap.error("--arm and --out are required unless --compare is given")

    _wrap_oracle()
    _wrap_node_sweep()
    arm = args.arm
    expect_fire: list[str] = []
    if arm == "vertex_free":
        _install_vertex_pin(True)
        expect_fire = ["pinned"]
    elif arm == "vertex_all":
        _install_vertex_pin(False)
        expect_fire = ["pinned"]
    elif arm == "noop":
        # ⭐ the harness's OWN falsification: the same wrapper, pinning nothing. Must be byte-identical
        #   to `base`, and if it is not, the wrapper itself is changing the answer (`TRAPS.md` A5).
        _install_vertex_pin(True, force_empty=True)
        expect_fire = ["init"]
    elif arm.startswith("ref_c="):
        _install_ref_exponent(float(arm.split("=", 1)[1]))
        expect_fire = ["ref_g", "ref_r"]
    elif arm != "base":
        ap.error(f"unknown arm {arm!r}")

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig()
    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    with args.out.open("w") as fh:
        for name in names:
            t0 = time.perf_counter()
            before = dict(_FIRED)
            cond = args.suite / name
            truth = P0.truth_f_gdna(cond) or 0.0
            m = P0.measure_condition(
                bam=str(cond / "sim_oracle.bam"), index=index, pipeline_config=PipelineConfig(),
                calibration_config=config, work_dir=args.work_dir / "rigel_pass0_oracle", tag=name,
                truth_pmfs=lambda size, d=cond: (
                    P0.truth_length_pmf(d, "gdna", size), P0.truth_length_pmf(d, "rna", size)
                ),
                oracle_cache=args.oracle_cache,
            )
            _FIRED["conditions"] += 1
            fired = {k: _FIRED[k] - before[k] for k in _FIRED}
            for k in expect_fire:
                if fired.get(k, 0) <= 0:
                    raise SystemExit(
                        f"⛔ A10: arm {arm!r} did not fire on {name} (counter {k} = 0). "
                        "An override that never ran reads as 'no effect'."
                    )
            for axis in ("node", "edge"):
                s = SA.summarise(SA.audit(m, axis=axis, config=config))
                sc = m.scores["pass0"][axis]["ALL"]
                s["mwae_all"] = float(sc.mwae)
                s["abs_err_all"] = float(sc.abs_err)
                s["mass_all"] = float(sc.mass)
                s["net_err_all"] = float(sc.net_err)
                fin = m.scores["final"][axis]["ALL"]
                s["mwae_all_final"] = float(fin.mwae)
                s["abs_err_all_final"] = float(fin.abs_err)
                s["library_f_gdna_pass0"] = float(m.library_f_gdna.get("pass0", float("nan")))
                s["library_f_gdna_final"] = float(m.library_f_gdna.get("final", float("nan")))
                s["library_f_gdna_truth"] = float(m.library_f_gdna.get("T", float("nan")))
                s["pinned"] = int(fired.get("pinned", 0))
                fh.write(json.dumps({"arm": arm, "condition": name, "axis": axis,
                                     "f_gdna": truth, **s}) + "\n")
                fh.flush()
            # ⚠ A14: the opportunity count, printed beside the result.
            print(f"  {name} {time.perf_counter() - t0:.0f}s   pinned={fired.get('pinned', 0)}"
                  f"  init_calls={fired.get('init', 0)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
