#!/usr/bin/env python
"""What is knowing the truth at the parameter-vertex objects worth, on the real ladder?

A silent gene's regions are pure gDNA (``f_g = 1`` exactly), the introns of a gene with no nascent
fragment are pure gDNA, and every counted object of a zero-gDNA library is pure RNA (``f_g = 0``
exactly). This harness hands those objects their exact answer from the certified slot truth and
RE-SOLVES the whole chain: under the two-phase transfer policy a node's neighbours receive its own
claim, so the pin is a delta at the truth installed as the node's claim (`messages.transfer._claims`)
and delivered as its ψ row (`PsiMessage.lam_rows`), never a substitution after the fact. The
population is the PARAMETER vertex (`_parameter_vertex`), never the realized one: an object whose few
fragments all happened to be gDNA is chance, and pinning chance as certainty prices nothing. Every arm
counts its own firings and raises if it did not fire; `noop` runs the whole wrapper and pins nothing,
so it must be byte-identical to `base`; `vertex_free` (no own composition evidence, the reachable
population) is the ceiling and `vertex_all` the looser bound. Two mechanism prototypes share the
harness so a ceiling and a mechanism are directly comparable: `ref_c=<a>[,<b>]` drives ψ's two Beta
reference exponents, `psi_mean` reports ``f_g`` as the posterior mean. The result prices missing
information, not headroom for a fix. Read ``mwae_all`` and ``abs_err_all`` (fixed denominators) and
never the honesty columns alone: an arm that changes what counts as solvable changes its own
denominator. Scoring is `solvability_audit.audit` over `pass0_vs_oracle.measure_condition`.

Usage::

    python scripts/design/vertex_ceiling.py --self-test
    python scripts/design/vertex_ceiling.py --arm base --conditions <cond> --out base.jsonl
    python scripts/design/vertex_ceiling.py --arm vertex_free --oracle-cache <suite>/oracle_cache --out free.jsonl
    python scripts/design/vertex_ceiling.py --arm ref_c=0.5,2.0 --out ref.jsonl
    python scripts/design/vertex_ceiling.py --compare base.jsonl free.jsonl
"""

from __future__ import annotations

import argparse
import contextlib
import inspect
import io
import json
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
sys.path.insert(0, str(DESIGN))


from _shared import sibling  # noqa: E402


SA = sibling("solvability_audit.py")
P0 = sibling("pass0_vs_oracle.py")

from rigel.calibration import region_init as NI, sweep as SW  # noqa: E402
from rigel.calibration import simplex_logodds as SL  # noqa: E402
from rigel.calibration.messages import PsiMessage  # noqa: E402
from rigel.calibration.messages import transfer as TR  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

CAL = sys.modules["rigel.calibration.calibrate"]

_EPS = 1.0e-9
#: the ceiling's own classification of "this object has no composition evidence of its own"; a
#: classification for a ceiling, not a production predicate (TRAPS: a-threshold-on-a-fitted-residue).
#: The `vertex_all` arm exists so the filter's effect is measured rather than assumed.
_TAU_FREE = 1.0e-4

#: filled by the wrappers, one call before `build_region_init` needs them.
_CTX: dict = {}
#: per-arm firing counters; an expected counter left at zero raises (TRAPS: an-ablation-that-never-ran).
_FIRED: dict = {
    "init": 0, "pinned": 0, "claimed": 0, "delivered": 0, "ref_g": 0, "ref_r": 0, "psi_mean": 0,
    "conditions": 0,
}
#: the pin's shape on the solve grid: everything but the cell nearest the truth is impossible
_PIN_WALL = -1.0e6


# ── the plumbing: get the oracle's per-object truth and the geometry to `build_region_init` ────────────


def _wrap_solve_chain():
    """Stash `region_arrays`: `solve_chain` receives it and calls `build_region_init` after.

    The wrapped name is a patch target the self-test checks for presence and identity, because a
    wrapped name that vanishes kills the instrument while the suite stays green."""
    orig = CAL.solve_chain

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _CTX["region_arrays"] = region_arrays
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.solve_chain = wrapper


def _install_psi_mean():
    """Report ``f_g`` as the posterior mean instead of the shipped ½-quantile, and price it.

    `_posterior_median_fg` is the single function both the single-strand and the AMBIG solve call, so
    one patch reaches both."""
    real = SL._posterior_median_fg

    def as_mean(post, lam, fg):
        _FIRED["psi_mean"] += 1
        return np.asarray(post, np.float64) @ np.asarray(fg, np.float64)

    SL._posterior_median_fg = as_mean
    del real


def _pin_row(lam, f_true: float) -> np.ndarray:
    """A delta at the true composition on the solve grid: zero at the cell nearest ``logit(f_true)``
    (a vertex lands on the grid's end cell) and a wall everywhere else, i.e. a certain claim in the
    currency every rule, lane and ψ read (a max-normalised log-profile over ``lam``)."""
    lam = np.asarray(lam, np.float64)
    f = float(np.clip(f_true, _EPS, 1.0 - _EPS))
    target = float(np.clip(np.log(f / (1.0 - f)), lam[0], lam[-1]))
    row = np.full(lam.shape[0], _PIN_WALL)
    row[int(np.argmin(np.abs(lam - target)))] = 0.0
    return row


class _PinnedPolicy(TR.TransferPolicy):
    """The shipped policy with the pinned nodes' ψ rows replaced by the delta at their truth: the
    node's own answer becomes the truth, whatever its evidence and its neighbours say. (What the node
    SENDS is pinned in `_claims`, see `_install_vertex_pin`.)"""

    name = "transfer"  # what the capture stamps: the instruments' "the arm ran" witness reads it

    def prepare(self, ctx, library):
        prepared = super().prepare(ctx, library)
        pins = _CTX.get("pins") or {}
        if not pins:
            return prepared
        lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), int(ctx.n_grid))
        n = int(ctx.n_slots)
        orig_solve = prepared.solve

        def solve(from_left, from_right):
            msg = orig_solve(from_left, from_right)
            rows = np.zeros((n, lam.shape[0])) if msg.lam_rows is None else np.array(msg.lam_rows, np.float64)
            for i, f_true in pins.items():
                rows[int(i)] = _pin_row(lam, f_true)
                _FIRED["delivered"] += 1
            return PsiMessage(lam_rows=rows, cube_rows=msg.cube_rows)

        prepared.solve = solve
        return prepared


def _parameter_vertex(chain, index, ra, slot_truth: dict, library_f_gdna: float):
    """The parameter-vertex population per slot, and its true composition: objects whose composition
    is a vertex by construction, never objects whose few fragments landed on one by chance (a boundary
    of a handful of crossings that all happened to be gDNA, pinned as certain, propagates chance).

    Three sources, each a parameter: (i) every region of a silent gene (no RNA fragment of any kind in
    any of its regions) that no expressed gene overlaps, ``f_g = 1``; (ii) every intron region of a
    gene with no nascent fragment in any of its introns that no nascent-bearing gene overlaps (mature
    RNA cannot be contained in an intron, so its only RNA is nascent), ``f_g = 1``; (iii) every live
    slot of a zero-gDNA library, ``f_g = 0``. A boundary joins (i)/(ii) when both its flanks do.
    Returns ``nan`` at every slot outside the population."""
    n = int(chain.n_slots)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region = kind == REGION
    n_rna = np.asarray(slot_truth["n_mrna"], np.float64) + np.asarray(slot_truth["n_nrna"], np.float64)
    n_nas = np.asarray(slot_truth["n_nrna"], np.float64)
    stratum = np.asarray(slot_truth["stratum"]).astype(str)
    count = np.asarray(slot_truth["count"], np.float64)
    f_true = np.full(n, np.nan)
    if library_f_gdna <= 0.0:
        f_true[count > 0.0] = 0.0
        return f_true
    starts = np.asarray(ra.start, np.int64)
    ends = np.asarray(ra.end, np.int64)
    ref_id = np.asarray(ra.ref_id)
    # per REGION slot: covered by an expressed gene? by a nascent-bearing gene? by any gene?
    region_slot = np.full(starts.shape[0], -1, np.int64)
    region_slot[obj[is_region]] = np.flatnonzero(is_region)
    covered = np.zeros(n, bool)
    expressed_cover = np.zeros(n, bool)
    nascent_cover = np.zeros(n, bool)
    genes = index.g_df
    if "is_synthetic" in genes.columns:
        genes = genes[~genes["is_synthetic"].astype(bool)]
    for g in genes.itertuples(index=False):
        rid = index.ref_name_to_id.get(g.ref)
        if rid is None:
            continue
        inside = np.flatnonzero((ref_id == rid) & (starts >= int(g.start)) & (ends <= int(g.end)))
        slots = region_slot[inside]
        slots = slots[slots >= 0]
        if slots.size == 0:
            continue
        covered[slots] = True
        if n_rna[slots].sum() > 0.0:
            expressed_cover[slots] = True
        if n_nas[slots].sum() > 0.0:
            nascent_cover[slots] = True
    silent = is_region & covered & ~expressed_cover
    dark_intron = is_region & (stratum == "R intron") & covered & ~nascent_cover
    at_one = (silent | dark_intron) & (count > 0.0)
    f_true[at_one] = 1.0
    left, right = np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)
    for b in np.flatnonzero(~is_region & (count > 0.0)):
        lo, hi = left[b], right[b]
        if lo >= 0 and hi >= 0 and at_one[lo] and at_one[hi]:
            f_true[b] = 1.0
    return f_true


def _install_vertex_pin(evidence_free_only: bool, force_empty: bool = False):
    """The ceiling arm. At every object whose truth sits on a vertex of the composition simplex
    (``f_g`` exactly 0 or exactly 1): make its own claim the oracle's exact answer (a delta on the
    solve grid, what its neighbours receive through every rule and lane) and deliver the same delta
    into ψ at the node (what it answers), then let the two passes and the solve run on top of it.
    An interior object keeps its own answer, so this prices the vertex and nothing else.

    ``evidence_free_only`` restricts the pin to objects with no own composition evidence
    (``tau_lam <= _TAU_FREE``), the population a vertex fix can actually reach; the unrestricted arm is
    the looser bound. ``force_empty`` runs the whole wrapper and pins nothing (the `noop` arm).

    Three patches, each a live target of the self-test: ``sweep.build_region_init`` (the
    classification, the truth mapped onto slots and the evidence filter, runs there once per sweep
    before the policy prepares); ``messages.transfer._claims`` (the pinned claims);
    ``calibrate.TransferPolicy`` (the pinned ψ rows)."""
    orig = NI.build_region_init
    orig_claims = TR._claims

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        _FIRED["init"] += 1
        _CTX["pins"] = {}
        slot_truth, index, ra = _CTX.get("slot_truth"), _CTX.get("index"), _CTX.get("region_arrays")
        if slot_truth is None or index is None or ra is None:
            return ni
        true_fg = _parameter_vertex(chain, index, ra, slot_truth, float(_CTX.get("library_f_gdna", 1.0)))
        tau = np.array(ni.tau_lam, np.float64)
        at_vertex = np.isfinite(true_fg)
        if evidence_free_only:
            at_vertex &= tau <= _TAU_FREE
        if force_empty:
            # the `noop` arm: the whole wrapper runs (the oracle is read, the truth is mapped to
            # slots, the classification is evaluated) and then nothing is pinned. Byte-identical to
            # `base` is the assertion; anything else means the wrapper itself moves the answer.
            at_vertex[:] = False
        tgt = np.flatnonzero(at_vertex)
        _FIRED["pinned"] += int(tgt.size)
        _CTX["pins"] = {int(i): float(true_fg[i]) for i in tgt}
        return ni

    def claims(c):
        own = orig_claims(c)
        for i, f_true in (_CTX.get("pins") or {}).items():
            own[int(i)] = _pin_row(c.lam, f_true)
            _FIRED["claimed"] += 1
        return own

    SW.build_region_init = wrapper
    TR._claims = claims
    CAL.TransferPolicy = _PinnedPolicy


def _install_ref_exponent(a_value: float, b_value: float | None = None):
    """ψ's two reference exponents as free numbers instead of the single shipped ½.

    They are pseudo-counts and the pair is a Beta: ``a·log f_g + b·log(1−f_g)`` on the λ grid is
    ``Beta(a, b)`` in ``f_g`` (the Jacobian ``|df/dλ| = f(1−f)`` turns ``f^{a−1}(1−f)^{b−1}`` into
    ``f^a (1−f)^b``), so the pair has a strength ``a+b`` and a mean ``a/(a+b)``, and the shipped
    ``a = b = ½`` fixes the mean at ½. ``a = b = 0`` makes ψ improper on both sides
    (TRAPS: no-prior-means-haldane), so small exponents bound what a derived rule could buy and are
    never the rule itself.

    ``b_value = None`` keeps the two equal. Both replacements take the shipped parameter names
    (``global_logprior``, ``rna_logprior``), which the self-test's arity block checks."""
    b_value = a_value if b_value is None else b_value

    def _gdna_arm(lam, global_logprior=None):
        _FIRED["ref_g"] += 1
        ref = float(a_value) * SL._log_fg(lam)[None, :]
        if global_logprior is None:
            return ref
        return ref + np.asarray(global_logprior, np.float64)

    def _rna_arm(lam, rna_logprior=None):
        _FIRED["ref_r"] += 1
        ref = float(b_value) * SL._log1m_fg(lam)[None, :]
        if rna_logprior is None:
            return ref
        return ref + np.asarray(rna_logprior, np.float64)

    SL._gdna_arm = _gdna_arm
    SL._rna_arm = _rna_arm


# ── the comparison ──────────────────────────────────────────────────────────────────────────────────


def _compare(paths: list[Path]) -> int:
    """Read two or more arm files and print the per-axis deltas, the fixed-denominator columns first
    because those cannot be gamed by knowing less (TRAPS: honesty-metrics-reward-ignorance)."""
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
    # the two fixed-denominator columns come first, because they are the only two that cannot be
    # gamed by the solver knowing less; `solvable_mwae` and the honesty columns follow.
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
    for axis in ("region", "boundary"):
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
        # a byte-identical arm is not evidence of no change (TRAPS: byte-identity-gate), except for
        # `noop`, where it is the assertion the arm exists to make; say which of the two it is.
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
                      f" — if that arm was meant to CHANGE something, it did not fire (TRAPS: an-ablation-that-never-ran/TRAPS: byte-identity-gate).")
    return 0


# ── --self-test: perturb every gate, with no I/O and no solver ──────────────────────────────────────
# Patch-target drift (a wrapped name that vanishes, a replacement whose arity no longer matches the
# shipped one) kills a monkey-patching harness while the suite stays green
# (TRAPS: a-green-suite-hid-five-dead-instruments), so those are the first two blocks and each is
# perturbed against the shape that would kill it.


#: every name this harness rebinds, as ``(module, attribute, definition_module_or_None)``. When the third
#: entry is given, the attribute must be the same object as its definition: patching a re-export only
#: reaches the solver while the re-export is live.
def _patch_targets():
    return (
        (CAL, "solve_chain", SW),
        (CAL, "TransferPolicy", TR),
        (TR, "_claims", None),
        (NI, "build_region_init", None),
        (SW, "build_region_init", NI),
        (SW, "CompositionPriors", SL),
        (SL, "CompositionPriors", None),
        (SL, "_posterior_median_fg", None),
        (SL, "_gdna_arm", None),
        (SL, "_rna_arm", None),
    )


def _target_live(mod, attr, definition) -> bool:
    """Is this patch target present, callable, and still the object it is a re-export OF?"""
    obj = getattr(mod, attr, None)
    if obj is None or not callable(obj):
        return False
    return definition is None or obj is getattr(definition, attr, None)


def _same_params(a, b) -> bool:
    """Same parameter names in the same order. Names rather than count, because the shipped caller
    passes `rna_logprior` by keyword."""
    return [p.name for p in inspect.signature(a).parameters.values()] == [
        p.name for p in inspect.signature(b).parameters.values()
    ]


def _try(fn):
    """Call ``fn`` and return its value, or ``None`` if it raised.

    Without this a broken replacement takes the whole self-test down with a traceback and the FAIL row
    for the check that already caught it never prints: the arity block diagnoses the failure, and the
    numeric block two lines later is what would crash."""
    try:
        return fn()
    except Exception:  # noqa: BLE001 — the self-test's job is to REPORT a broken arm, not to inherit it
        return None


def self_test() -> int:
    checks: list[tuple[str, bool]] = []
    saved = {(m.__name__, a): getattr(m, a, None) for m, a, _ in _patch_targets()}

    def restore():
        _CTX.clear()
        for m, a, _ in _patch_targets():
            v = saved[(m.__name__, a)]
            if v is not None:
                setattr(m, a, v)

    # ── ① PATCH TARGETS: present, callable, and the same object as their definition ──────────────────
    dead = [f"{m.__name__.rsplit('.', 1)[-1]}.{a}" for m, a, d in _patch_targets()
            if not _target_live(m, a, d)]
    checks.append((f"all {len(_patch_targets())} patch targets are live", not dead))
    if dead:
        print(f"  ⛔ DEAD PATCH TARGETS: {', '.join(dead)}", flush=True)

    # perturbation: the same predicate must refuse a name that does not exist (`CAL.region_sweep`).
    # A checker that cannot fail is not a check.
    checks.append(("the predicate REFUSES the dead `CAL.region_sweep` name",
                   not _target_live(CAL, "region_sweep", None)))

    # perturbation: break one re-export and the identity half must fire, not just the presence half.
    SW.build_region_init = lambda *a, **k: None
    checks.append(("a re-export rebound to a stranger => target reads DEAD",
                   not _target_live(SW, "build_region_init", NI)))
    restore()
    checks.append(("…and restoring it reads live again", _target_live(SW, "build_region_init", NI)))

    # ── ② ARITY: every replacement must take the shipped function's parameter names ──────────────────
    _install_ref_exponent(0.5)
    checks.append(("ref_c's two arms match the shipped signatures",
                   _same_params(SL._gdna_arm, saved[(SL.__name__, "_gdna_arm")])
                   and _same_params(SL._rna_arm, saved[(SL.__name__, "_rna_arm")])))
    # perturbation: a one-argument `_rna_arm` (one fewer than the shipped signature) must be rejected
    # by the same comparison.
    checks.append(("a one-argument `_rna_arm` is REJECTED by the same comparison",
                   not _same_params(lambda lam: None, saved[(SL.__name__, "_rna_arm")])))

    # ── ③ ref_c reproduces the shipped reference at ½, and moves off it elsewhere ────────────────────
    lam = np.linspace(-10.0, 10.0, 21)
    ship_g = saved[(SL.__name__, "_gdna_arm")](lam, None)
    ship_r = saved[(SL.__name__, "_rna_arm")](lam, None)
    before = dict(_FIRED)
    half_g = _try(lambda: SL._gdna_arm(lam, None))
    half_r = _try(lambda: SL._rna_arm(lam, None))
    checks.append(("ref_c=0.5 is BIT-IDENTICAL to the shipped ½ reference, both arms",
                   half_g is not None and half_r is not None
                   and np.array_equal(half_g, ship_g) and np.array_equal(half_r, ship_r)))
    checks.append(("…and both fire counters moved (TRAPS: an-ablation-that-never-ran)",
                   _FIRED["ref_g"] > before["ref_g"] and _FIRED["ref_r"] > before["ref_r"]))
    restore()
    # perturbation: a different exponent must not reproduce it, or the arm is inert.
    _install_ref_exponent(0.25)
    q_g = _try(lambda: SL._gdna_arm(lam, None))
    checks.append(("ref_c=0.25 DIFFERS from the shipped reference",
                   q_g is not None and not np.array_equal(q_g, ship_g)))
    restore()
    # perturbation: the pair is a Beta(a,b), so the two arms must be drivable independently.
    _install_ref_exponent(0.5, 2.0)
    ab_g = _try(lambda: SL._gdna_arm(lam, None))
    ab_r = _try(lambda: SL._rna_arm(lam, None))
    checks.append(("ref=A,B moves the RNA arm alone",
                   ab_g is not None and ab_r is not None
                   and np.array_equal(ab_g, ship_g) and not np.array_equal(ab_r, ship_r)))
    restore()

    # ── ④ psi_mean really is the MEAN and not the median ─────────────────────────────────────────────
    _install_psi_mean()
    post = np.array([[0.6, 0.4]])       # a deliberately skewed 2-point posterior …
    fg = np.array([0.0, 1.0])            # … whose mean (0.4) and median (0.0) differ
    got = SL._posterior_median_fg(post, lam, fg)
    checks.append(("psi_mean returns the posterior MEAN (0.4), not the median (0.0)",
                   float(np.ravel(got)[0]) == 0.4 and _FIRED["psi_mean"] > 0))
    restore()

    # ── ⑤ the vertex pin's `noop` shape: the wrapper runs and returns the init UNTOUCHED ─────────────
    sentinel = object()
    NI.build_region_init = lambda chain, statics, geometry, **kw: sentinel
    _install_vertex_pin(True, force_empty=True)
    before = dict(_FIRED)
    out = SW.build_region_init(None, None, None)
    checks.append(("noop: wrapper fires, pins nothing, returns the init object itself",
                   out is sentinel and _FIRED["init"] > before["init"]
                   and _FIRED["pinned"] == before["pinned"]))
    restore()

    # ── ⑥ the re-pointed pin: a delta on the grid, the claim pinned, the ψ row delivered ─────────────
    lam21 = np.linspace(-10.0, 10.0, 21)
    top, bottom, mid = _pin_row(lam21, 1.0), _pin_row(lam21, 0.0), _pin_row(lam21, 0.5)
    checks.append(("the pin row is a delta: at the top cell for f_g = 1, the bottom for 0, the centre for ½",
                   top[-1] == 0.0 and bottom[0] == 0.0 and mid[10] == 0.0
                   and np.sum(top == 0.0) == 1 and np.all(top[:-1] == _PIN_WALL)))
    # the claims wrapper pins exactly the slots the sweep wrapper classified, and nothing else
    from types import SimpleNamespace

    NI.build_region_init = lambda chain, statics, geometry, **kw: sentinel
    TR._claims = lambda c: [None] * c.n
    _install_vertex_pin(True)
    _CTX["pins"] = {2: 1.0}
    before = dict(_FIRED)
    own = TR._claims(SimpleNamespace(lam=lam21, n=4))
    checks.append(("the pinned slot's claim is the delta at its truth; every other claim is untouched",
                   own[2] is not None and np.array_equal(own[2], top)
                   and own[0] is None and own[1] is None and own[3] is None
                   and _FIRED["claimed"] == before["claimed"] + 1))
    # the pinned policy delivers the delta into ψ at the pinned slot, on top of a silent solve
    pol = CAL.TransferPolicy()
    prepared = pol.prepare(SimpleNamespace(n_slots=4, n_grid=21, logodds_window=10.0, factory_rows=None), None)
    msg = prepared.solve([None] * 4, [None] * 4)
    checks.append(("the pinned policy delivers the delta as the slot's ψ row (fires `delivered`)",
                   msg.lam_rows is not None and msg.lam_rows.shape == (4, 21)
                   and np.array_equal(msg.lam_rows[2], top) and not msg.lam_rows[0].any()
                   and _FIRED["delivered"] == before["delivered"] + 1))
    # perturbation: with no pins the policy is the shipped one, so silent stays silent
    _CTX["pins"] = {}
    checks.append(("with nothing pinned the policy's solve is untouched (silent stays silent)",
                   CAL.TransferPolicy()
                   .prepare(SimpleNamespace(n_slots=4, n_grid=21, logodds_window=10.0, factory_rows=None), None)
                   .solve([None] * 4, [None] * 4).is_silent))
    restore()
    checks.append(("…and every patch target is restored after the pin",
                   not [a for m, a, d in _patch_targets() if not _target_live(m, a, d)]))

    # ── ⑥b the population is the PARAMETER vertex ────────────────────────────────────────────────────
    import pandas as pd

    # regions 0..6 on one reference: A(silent gene: exon,intron,exon) | B(expressed, nascent on) | C(expressed, nascent OFF: exon,intron,exon)
    # chain: R0 B0 R1 B1 R2 B2 R3 B3 R4 B4 R5 B5 R6  (13 slots, regions at even indices)
    kind = np.array([REGION if i % 2 == 0 else 1 - REGION for i in range(13)])
    objx = np.array([i // 2 for i in range(13)])
    chain = SimpleNamespace(n_slots=13, kind=kind, obj_idx=objx,
                            left=np.array([-1] + list(range(12))), right=np.array(list(range(1, 13)) + [-1]))
    ra_ = SimpleNamespace(start=np.arange(7) * 100, end=np.arange(1, 8) * 100, ref_id=np.zeros(7, int))
    gdf = pd.DataFrame({"ref": ["c"] * 3, "start": [0, 300, 400], "end": [300, 400, 700],
                        "g_id": ["A", "B", "C"], "is_synthetic": [False] * 3})
    index_ = SimpleNamespace(g_df=gdf, ref_name_to_id={"c": 0})
    strat = np.array(["R exon", "B", "R intron", "B", "R exon", "B", "R exon", "B", "R exon", "B", "R intron", "B", "R exon"])
    count = np.ones(13) * 10.0
    n_mrna = np.zeros(13)
    n_nrna = np.zeros(13)
    n_mrna[6] = 5.0            # gene B's exon (slot 6) is expressed …
    n_nrna[6] = 1.0            # … with nascent
    n_mrna[8] = 4.0  # gene C's exons expressed …
    n_mrna[12] = 3.0  # … and its intron (slot 10) has no nascent
    st = dict(n_mrna=n_mrna, n_nrna=n_nrna, stratum=strat, count=count)
    f = _parameter_vertex(chain, index_, ra_, st, 0.5)
    checks.append(("silent gene A: every region and its interior boundaries are at f_g = 1",
                   all(f[i] == 1.0 for i in (0, 1, 2, 3, 4))))
    checks.append(("expressed gene B with nascent: nothing pinned",
                   np.isnan(f[6]) and np.isnan(f[5]) and np.isnan(f[7])))
    checks.append(("expressed gene C without nascent: its INTRON is at 1, its exons are not",
                   f[10] == 1.0 and np.isnan(f[8]) and np.isnan(f[12]) and np.isnan(f[9]) and np.isnan(f[11])))
    # perturbation: one nascent fragment in C's intron un-pins it
    st2 = dict(st, n_nrna=np.where(np.arange(13) == 10, 1.0, n_nrna))
    checks.append(("one nascent fragment in the intron un-pins it (chance is not the population)",
                   np.isnan(_parameter_vertex(chain, index_, ra_, st2, 0.5)[10])))
    # perturbation: an expressed gene overlapping the silent one un-pins the shared regions
    gdf2 = pd.concat([gdf, pd.DataFrame({"ref": ["c"], "start": [0], "end": [200], "g_id": ["D"], "is_synthetic": [False]})])
    st3 = dict(st, n_mrna=np.where(np.arange(13) == 0, 2.0, n_mrna))
    f3 = _parameter_vertex(chain, SimpleNamespace(g_df=gdf2, ref_name_to_id={"c": 0}), ra_, st3, 0.5)
    checks.append(("RNA inside a silent gene's span (an overlapping gene's) un-pins its EXONS — a slot's RNA "
                   "is not attributable, so the classification is conservative — while its nascent-free INTRON stays at 1",
                   np.isnan(f3[0]) and f3[2] == 1.0 and np.isnan(f3[4])))
    # the zero-gDNA library: every counted slot at 0
    f0 = _parameter_vertex(chain, index_, ra_, st, 0.0)
    checks.append(("a zero-gDNA library pins every counted slot at f_g = 0", bool(np.all(f0 == 0.0))))

    # ── ⑦ the comparator: byte-identical must be LABELLED, and differently for `noop` ────────────────
    def _rows(arm, bump=0.0):
        return "".join(
            json.dumps({"arm": arm, "condition": f"c{i}", "axis": ax,
                        "mwae_all": 0.10 + bump, "abs_err_all": 1000.0}) + "\n"
            for i in range(2) for ax in ("region", "boundary")
        )

    with tempfile.TemporaryDirectory() as td:
        base_p = Path(td) / "base.jsonl"
        base_p.write_text(_rows("base"))
        for arm, bump, want in (("noop", 0.0, "own falsification PASSES"),
                                ("vertex_free", 0.0, "did not fire"),
                                ("vertex_free", 0.05, None)):
            p = Path(td) / f"{arm}_{bump}.jsonl"
            p.write_text(_rows(arm, bump))
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                _compare([base_p, p])
            txt = buf.getvalue()
            if want is None:
                checks.append(("a MOVED arm is called neither identical nor unfired",
                               "own falsification PASSES" not in txt and "did not fire" not in txt))
            else:
                checks.append((f"identical `{arm}` is labelled {want!r}", want in txt))
        # perturbation: one arm alone is not a comparison and must be refused.
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = _compare([base_p])
        checks.append(("a single arm is REFUSED rather than compared with itself", rc == 1))

    width = max(len(name) for name, _ in checks)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}", flush=True)
    failed = [name for name, ok in checks if not ok]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} harness gates fire", flush=True)
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default=None,
                    help="base | noop | psi_mean | vertex_free | vertex_all | ref_c=<a>[,<b>]")
    ap.add_argument("--self-test", action="store_true",
                    help="perturb every harness gate; no I/O, no solver")
    ap.add_argument("--compare", nargs="*", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_vertex_ceiling"))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.compare:
        return _compare(args.compare)
    if not args.arm or not args.out:
        ap.error("--arm and --out are required unless --compare or --self-test is given")

    _wrap_solve_chain()
    arm = args.arm
    expect_fire: list[str] = []
    if arm == "psi_mean":
        # f_g as the posterior mean, which closes the simplex exactly.
        _install_psi_mean()
        expect_fire = ["psi_mean"]
    elif arm == "vertex_free":
        _install_vertex_pin(True)
        expect_fire = ["pinned", "claimed", "delivered"]
    elif arm == "vertex_all":
        _install_vertex_pin(False)
        expect_fire = ["pinned", "claimed", "delivered"]
    elif arm == "noop":
        # the harness's own falsification: the same wrapper, pinning nothing. Must be byte-identical
        # to `base`; if it is not, the wrapper itself is changing the answer.
        _install_vertex_pin(True, force_empty=True)
        expect_fire = ["init"]
    elif arm.startswith("ref_c="):
        # `ref_c=<a>` keeps both arms equal; `ref_c=<a>,<b>` drives them unequal (the Beta(a,b) design).
        _spec = arm.split("=", 1)[1]
        _install_ref_exponent(*(float(x) for x in _spec.split(",")))
        expect_fire = ["ref_g", "ref_r"]
    elif arm != "base":
        ap.error(f"unknown arm {arm!r}")

    index = TranscriptIndex.load(str(args.index))
    _CTX["index"] = index
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
            _CTX["library_f_gdna"] = float(truth)
            _CTX["slot_truth"] = None
            if arm in ("vertex_free", "vertex_all", "noop"):
                st = (args.oracle_cache or Path("/nonexistent")) / name / "slot_truth.npz"
                if not st.is_file():
                    raise SystemExit(
                        f"⛔ the vertex population is the PARAMETER vertex and needs the certified "
                        f"slot_truth: pass --oracle-cache (missing {st})"
                    )
                _CTX["slot_truth"] = dict(np.load(st, allow_pickle=True))
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
                        f"⛔ TRAPS: an-ablation-that-never-ran: arm {arm!r} did not fire on {name} (counter {k} = 0). "
                        "An override that never ran reads as 'no effect'."
                    )
            for axis in ("region", "boundary"):
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
            # the opportunity count, printed beside the result: an arm with zero opportunities is
            # not a control (TRAPS: could-the-arm-have-fired).
            print(f"  {name} {time.perf_counter() - t0:.0f}s   pinned={fired.get('pinned', 0)}"
                  f"  init_calls={fired.get('init', 0)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
