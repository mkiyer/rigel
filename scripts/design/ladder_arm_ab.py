#!/usr/bin/env python
"""⭐⭐ PROTOTYPE A SOLVER MECHANISM AS A PER-OBJECT OVERRIDE AND SCORE IT ON THE REAL LADDER — before
writing a line of it into ``src/``.

⛔ **THE WORKFLOW THIS EXISTS FOR.** A toy isolates a mechanism; it cannot rank one (`TESTING.md` §0b).
`toy_ceiling.py` re-solves the owner's two-exon toy under a set of arms and says what each is worth
THERE. This applies the SAME override to the 36-condition gDNA ladder and scores it with the project's
own acceptance instrument (`solvability_audit.py`), so the two arms differ by exactly the mechanism and
the answer is in the currency the project actually ships.

⭐ **The override goes in at `node_init.build_node_init`** — the per-object message-free self-solve —
which is the one place a mechanism can be expressed without touching either relay twin, so a prototype
cannot accidentally be gated in one twin and not the other (`TRAPS.md` B14). `region_arrays` reaches the
wrapper by way of a `node_sweep` wrapper, because `build_node_init` is not handed it.

⚠ **A prototype here is an UPPER bound on the built form, not a model of it**: it overwrites the
object's own belief unconditionally, with no licence and no inverse-variance fuse. If the upper bound
loses, the built form cannot win; if it wins, the built form still has to earn it.

⛔⛔ **THE ARM IT CARRIES, AND ITS VERDICT (2026-08-04).** ``intron_phi`` is `EQUATIONS.md` §3.6 face (I):
an `intron|exon` EDGE takes the flanking INTRON's COMPOSITION — never its level — and closes the level
with its own observed mass. On the toy it was worth −0.009 panel-wide and −0.027 on capture-ON ×
unstranded. On the ladder, 32 contaminated conditions, node axis:

======================  ==========  =============  ===================
metric                  base        intron_phi     better/worse/flat
======================  ==========  =============  ===================
solvable mwae           0.0413      **0.0426**     14 / 18 / 0
confidently-wrong       20,173      **22,336**      6 / 26 / 0
declared-precision      2.688       **2.821**       6 / 26 / 0
relay Δ                 +0.0021     **+0.0034**    14 / 18 / 0
======================  ==========  =============  ===================

⛔ **REFUTED.** And on the EDGE axis it is worse still (mwae 0.0216 → 0.1449) because it ADMITS those
EDGEs into the scored population (``solv%`` 43.1 → 48.2) carrying a wrong composition — `TRAPS.md` B13
inverted. ⭐ The one real signal in it: the split is by STRAND, not by capture. Every unstranded row is
worse and every stranded capture-OFF row is better, which is the DL mismatch test having nothing to
check an imported composition against when the destination has no self-solve.

⭐ C11 is the precedent that makes this run mandatory rather than optional: its toy delta was LARGER
(−0.035) and the ladder said it cost the panel nothing. `TRAPS.md` B1.

---

⭐⭐⭐ **THE ZERO-CLAIM DECOMPOSITION ARMS (`zc_*`).** The 39 % win of 2026-08-06 was three distinct
behavioural changes landed across two commits, and one coherent stratum — stranded × capture-ON at
``g10`` and above — got **20–34 % worse**. These arms revert **exactly one** of the three, holding the
other two at HEAD, so the regression can be attributed rather than guessed:

===================  ==============================================  ==========================
arm                  what it reverts                                 what it tests
===================  ==============================================  ==========================
``zc_own_count``     ``own_precision``'s count TERM back to the       CAN a zero-count object
                     ``n/(n·V+1)`` form, ``0`` at ``n = 0``           FORM a claim?
``zc_total_n``       the count ARGUMENT back to the slot TOTAL        does a component's
                     ``n_node``, from HEAD's COMPONENT count          precision come from its
                     ``f_c·M``                                       OWN count or the total?
``zc_live_count``    the ``live`` predicate back to COUNT-keyed —     CAN a zero-DENSITY
                     ``rho_g > 0`` on gDNA, ``n > 0 ∧ rho > 0``       component form one?
                     on each RNA strand                              (⭐ all three streams)
``zc_transfer``      ``composition_logvar``'s count term back to      CAN a zero-count object
                     ``1/n`` (``∞`` at 0) ⇒ ``σ²_transfer = ∞`` ⇒     TRANSMIT one?
                     every message it sends is annihilated
===================  ==============================================  ==========================

⚠ ``zc_total_n`` is the change that was NOT in the handoff's list and is the largest of the four by
construction: on a stranded, gDNA-rich exon the RNA arms' component count ``f_c·M`` is a few hundredths of
``n_node``, so their precision falls by orders of magnitude, while an unstranded object splits its RNA
evenly and barely moves.

⛔ Each arm asserts it FIRED (`TRAPS.md` A10 — an ablation that never ran reads as "no effect", and
``composition_logvar`` is bound as a module global in ``bp_solver`` as well as in its own leaf module, so
patching one name is exactly the miss A10 records). ``--arm base`` against itself is the byte-identity
falsification (`TRAPS.md` A5).
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import subprocess
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

import rigel.calibration.strand_balance  # noqa: E402,F401  (registers the module for patching)
from rigel.calibration import bp_solver, node_init as NI  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.calibration.node_geometry import g1_locked, node_global_geometry  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

CAL = sys.modules["rigel.calibration.calibrate"]

_EPS = 1.0e-9
INTRON, EXON = 1, 2
#: set by the `node_sweep` wrapper one call before `build_node_init` needs it.
_RA: dict = {}


def _install_kappa_half():
    """⭐⭐ THE ARM `ROADMAP.md` §1 step 3 ASKS FOR IN ITS OWN WORDS: "inject κ = 0.5 exactly and diff".

    On a genuinely unstranded library the strand channel carries **exactly** zero information about
    composition — the Fisher information is ``∝ (2κ−1)²`` and κ is ½. But κ is FITTED, so it lands on
    0.500689, and the tool reads that 0.000689 as evidence: ``disc = 4·max(0, (κ̂−½)² − σ²_d)`` survives
    the derived noise floor about one unstranded library in three, and a **boolean** composition licence
    is then flipped at full strength by it.

    Forcing κ to exactly ½ makes ``disc`` identically 0, which is what propagating ``Var(κ̂)`` would
    achieve by derivation — so on the UNSTRANDED conditions this is the ceiling on that fix.
    ⛔ On the STRANDED conditions it is not a fix at all, it is a DESTRUCTION TEST: κ̂ is ~97 standard
    errors from ½ there and the channel is carrying real signal, so those rows must get much worse. If
    they do not, the lever is not connected and no reading of the unstranded rows means anything.

    ⭐ Only κ moves. ``n_observations`` is carried through unchanged, so the noise floor ``σ²_d`` — which
    depends on the SAMPLE SIZE, not on κ — is identical in both arms; overriding it too would conflate
    "κ is exactly ½" with "there are no spliced reads"."""
    import dataclasses as _dc

    SB = sys.modules["rigel.calibration.strand_balance"]
    orig = CAL.fit_strand_balance

    def wrapper(strand_model):
        return _dc.replace(orig(strand_model), rna_sense_frac=0.5)

    CAL.fit_strand_balance = wrapper
    SB.fit_strand_balance = wrapper


#: ⛔ A10 — every ablation increments a counter and `main` RAISES if it did not fire. An arm that never
#: ran scores byte-identical to base, which reads as "this change is inert" and is publishable and wrong.
_FIRED: dict = {}


def _fire(name: str) -> None:
    _FIRED[name] = _FIRED.get(name, 0) + 1


def _install_zc_own_count():
    """⭐ REVERT 1 of 3 — ``own_precision``'s count term, so a ZERO COUNT cannot form a claim.

    HEAD: ``p = 1/(Var(log f) + trigamma(n + ½))`` — finite at every ``n``, including 0.
    Here: ``p = n/(n·Var(log f) + 1)`` gated on ``n > 0`` — algebraically ``1/(Var + 1/n)``, the
    large-count LIMIT of the same expression, which diverges at ``n = 0`` and so returns 0 there.

    ⚠ ``own_precision`` is called as a MODULE GLOBAL from inside ``build_node_init``, so rebinding
    ``NI.own_precision`` is sufficient and `bp_solver`'s re-export of ``build_node_init`` shares the same
    function object. The counter proves it rather than the reader assuming it."""
    orig = NI.own_precision

    def wrapper(n, v_log, live):
        _fire("zc_own_count")
        nn = np.asarray(n, np.float64)
        v = np.asarray(v_log, np.float64)
        ok = np.isfinite(v)
        v_fin = np.where(ok, v, 0.0)
        return np.where(
            (nn > 0.0) & ok & np.asarray(live, bool),
            nn / (nn * np.maximum(v_fin, 0.0) + 1.0),
            0.0,
        )

    NI.own_precision = wrapper
    return orig


def _rebuild_own(ni, statics, geometry, *, total_n: bool, count_live: bool):
    """Rebuild the three own-belief precisions from the arrays ``build_node_init`` itself was handed,
    varying ONE of its two count/gate decisions.

    ⚠ Done OUTSIDE ``build_node_init`` because both quantities are locals there. Every input comes from
    the geometry/statics the solver was given and the composition variance is HEAD's own
    ``own_composition_logvar``, so an arm differs from HEAD in its named decision and in nothing else.
    ⭐ The LOCATIONS (``rho_g``/``rho_pos``/``rho_neg``) are HEAD's untouched unless ``count_live``
    zeroes a component, which is exactly what the pre-fix predicate did."""
    M = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)  # ≡ node_global_geometry mass
    E_r = np.asarray(geometry.eff_rna, np.float64)
    f_g = np.asarray(ni.f_g, np.float64)
    v_fg, v_fr = NI.own_composition_logvar(f_g, ni.tau_lam, np.asarray(ni.struct_lock, bool))
    rho_g = np.asarray(ni.rho_g, np.float64)
    rho_p = np.asarray(ni.rho_pos, np.float64)
    rho_n = np.asarray(ni.rho_neg, np.float64)

    E_g = np.asarray(node_global_geometry(geometry)[1], np.float64)
    live_g = (rho_g > _EPS) if count_live else (E_g > 0.0)
    n_g = M if total_n else f_g * M
    prec_g = NI.own_precision(n_g, v_fg, live_g)

    def _rna(rho, free_s):
        admissible = np.asarray(free_s, bool)
        live = (
            admissible & (M > 0.0) & (rho > _EPS) if count_live else admissible & (E_r > 0.0)
        )
        out = np.where(live, rho, 0.0)
        n_c = M if total_n else rho * E_r  # rho·E_r ≡ f_c·M, the component's own count
        return out, NI.own_precision(n_c, v_fr, live)

    rho_pos, prec_pos = _rna(rho_p, statics.free_pos)
    rho_neg, prec_neg = _rna(rho_n, statics.free_neg)
    return NI.NodeInit(
        f_g=ni.f_g, f_pos=ni.f_pos, f_neg=ni.f_neg,
        rho_g=rho_g, rho_pos=rho_pos, rho_neg=rho_neg,
        prec_g=prec_g, prec_pos=prec_pos, prec_neg=prec_neg,
        struct_lock=ni.struct_lock, tau_lam=ni.tau_lam,
    )


def _install_rebuild(name: str, **flags):
    """Wrap ``build_node_init`` with one :func:`_rebuild_own` variant. ⚠ Both bindings are rebound —
    ``bp_solver`` re-exports the name (A10)."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        _fire(name)
        return _rebuild_own(ni, statics, geometry, **flags)

    NI.build_node_init = wrapper
    bp_solver.build_node_init = wrapper
    return orig


def _install_zc_transfer():
    """⭐ REVERT 3 of 3 — ``composition_logvar``'s count term back to ``1/n``, so a zero-count object
    cannot TRANSMIT a claim however precisely it holds one.

    ``Var(log ρ_tot) = ∞`` ⇒ ``σ²_transfer = logvar_tot[dst] + logvar_tot[src] = ∞`` ⇒
    ``1/(1/p + ∞) = 0`` on all three streams of every message that slot sends, the MEASUREMENT stream
    included. That is the annihilation `TRAPS.md` C0d records, and it is the half of the fix the 39 %
    was actually attributed to.

    ⛔⛔ **BOTH BINDINGS ARE PATCHED, and that is A10 verbatim.** ``bp_solver`` does
    ``from .enrichment_frame import composition_logvar``, which makes a SEPARATE module global; patching
    only the leaf module leaves the solver calling the original and the arm reads as inert.
    ⚠ ``transfer_logvar``'s ``~isfinite`` guard is NOT restored — it was deleted as treating the
    symptom, and restoring it here would vary two things."""
    EF = sys.modules["rigel.calibration.enrichment_frame"]
    orig = EF.composition_logvar

    def wrapper(f_g, E_g, E_r, var_fg, n):
        _fire("zc_transfer")
        nn = np.asarray(n, np.float64)
        with np.errstate(divide="ignore"):
            asym = np.where(nn > 0.0, 1.0 / np.maximum(nn, 1e-300), np.inf)
        return orig(f_g, E_g, E_r, var_fg, n) - EF.count_logvar(nn) + asym

    EF.composition_logvar = wrapper
    bp_solver.composition_logvar = wrapper
    return orig


def _zero_mass_locked(ni, geometry):
    """The slot set the two arms below act on: **zero observed mass AND structurally composition-certain.**

    ⭐ These are the ONLY slots whose own gDNA density is exactly 0 while their own precision is positive.
    ``rho_g = f_g·M/E_g`` is 0 iff ``M = 0`` (``f_g`` is a ψ grid point and never exactly 0), and a
    zero-mass slot has ``tau_lam = 0``, so ``own_composition_logvar`` returns ``∞`` and the precision is 0
    — **unless** ``struct_lock`` short-circuits the variance to 0. So "zero density with a live claim"
    is exactly ``struct_lock ∧ M = 0``: an empty intergenic NODE.
    ⚠ Printed as a count by both callers, because an arm acting on an empty set is not a control (A14)."""
    M = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
    return np.asarray(ni.struct_lock, bool) & (M <= 0.0), M


def _install_zc_anchor_mute():
    """⭐⭐⭐ THE SURGICAL ARM — mute the gDNA precision at EMPTY STRUCTURALLY-LOCKED anchors and change
    nothing else anywhere.

    ⛔ This is not a candidate fix; it is the **localisation**. ``zc_own_count`` and ``zc_live_count`` both
    remove the empty anchor's claim, but each also changes thousands of other slots, so neither can say the
    anchor is the site. This changes ``prec_g`` at ~122 slots at ``g75 capture_on`` and at 1,312 at ``g00``
    — so if the stranded × capture-ON regression goes and the ``g00`` win goes with it, the two are the SAME
    mechanism and a fix has to separate them by something other than muting."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        sel, _M = _zero_mass_locked(ni, geometry)
        _FIRED["zc_anchor_mute_slots"] = int(sel.sum())
        _fire("zc_anchor_mute")
        prec_g = np.array(ni.prec_g, np.float64)
        prec_g[sel] = 0.0
        return dataclasses.replace(ni, prec_g=prec_g)

    NI.build_node_init = wrapper
    bp_solver.build_node_init = wrapper
    return orig


def _install_zc_jeffreys_mean():
    """⭐⭐ THE CANDIDATE FIX — at a slot with **zero observed mass**, the gDNA density is the Jeffreys
    posterior MEAN ``½/E_g`` instead of the point value ``0``.

    ``a`` events over exposure ``E`` under the Jeffreys prior ψ is built on has posterior
    ``Gamma(a + ½, E)`` — the same posterior :func:`count_logvar` takes its variance from. Its mean is
    ``(a + ½)/E``, so at ``a = 0`` it is ``½/E``: **a long empty region says ~0 and a short empty region
    says "below about 1/E"**, which is `TRAPS.md` C0c's own sentence ("zero over 50.7 Mb and zero over
    200 bp are the same number to this code and opposite statements about the world") applied to the
    LOCATION rather than to the variance.

    ⭐ **The mass-pin objection does not reach ``M = 0``.** ``Σ_c ρ_c·E_c = M`` is what forbade the
    ``(a+½)/E`` location in general (three components each gaining ½ breaks it by 3/2 —
    `test_relay_mass_pin`'s ``R_own = 0.5 + 1/M``). At ``M = 0`` there is no mass to partition, the
    identity is vacuous, and the live components are independent Poisson RATES rather than shares of a
    total. So this is scoped to exactly the slots where the constraint that refused it does not apply.

    ⚠ It does NOT fix a level: an off-probe anchor's ``½/E_g`` is still far below an on-probe exon's true
    density. It fixes the CLAIM, and the arm exists to price that separately from the level."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        sel, _M = _zero_mass_locked(ni, geometry)
        _FIRED["zc_jeffreys_mean_slots"] = int(sel.sum())
        _fire("zc_jeffreys_mean")
        E_g = np.asarray(node_global_geometry(geometry)[1], np.float64)
        rho_g = np.array(ni.rho_g, np.float64)
        ok = sel & (E_g > 0.0)
        rho_g[ok] = 0.5 / E_g[ok]
        return dataclasses.replace(ni, rho_g=rho_g)

    NI.build_node_init = wrapper
    bp_solver.build_node_init = wrapper
    return orig


def _install_zc_reference_var():
    """⭐⭐⭐ THE CANDIDATE THE PERTURBATION PRODUCED — ``Var(f_g)`` at an EVIDENCE-FREE slot is the variance
    of ψ's OWN reference prior, ``Beta(½, ½)`` ⇒ exactly **1/8**, instead of ``f_g(1−f_g)`` evaluated at a
    default point estimate that sits on the corner.

    ⛔ **THE DEFECT.** ``node_sweep`` builds
    ``_var_fg = where(struct_lock, 0, where(τ > 0, min(fgfr²/τ, fgfr), fgfr))`` and hands it to
    `composition_logvar`, whose output IS ``σ²_transfer``. An evidence-free slot's default belief is
    ``f_g = 1`` **exactly** (measured on the ladder fixture: every zero-count slot), so ``fgfr = 0`` and the
    composition half of ``Var(log ρ_tot)`` VANISHES. ``logvar_tot`` is then ``count_logvar(n)`` alone — which
    is why replacing ``1/n`` (``∞`` at zero counts) left those hops with no damping whatsoever and turned a
    zero-mass slot from a relay BARRIER into a CONDUIT.

    ⭐ **NO TUNED CONSTANT.** ``f_g(1−f_g)`` is a BERNOULLI variance: the right scale away from the corner,
    and an assertion of exact knowledge at it. With ``τ_λ = 0`` there is no evidence at all, so the honest
    variance is the reference the solver is already built on — `simplex_logodds._JEFFREYS_REF` makes ψ's
    composition reference ``Beta(½, ½)``, whose variance is ``¼/(1·2) = 1/8`` and does not depend on where
    the point estimate sits.

    ⭐⭐ **AND IT IS SELF-TARGETING.** The coefficient it multiplies is ``[(1/E_g − 1/E_r)/B]²``, which
    ``composition_logvar``'s own docstring says "diverges as ``E_g`` collapses on a short region" — so the
    restored damping is LARGEST at the short, capture-depleted slots the regression lives on, and ~0.036 on a
    long contained region where the premise is fine.

    ⛔ ``struct_lock`` slots are untouched (they pass ``0.0`` exactly and keep it): a structurally pure-gDNA
    object IS composition-certain, so the zero-gDNA anchor's true claim is still undamped. That is what makes
    this different from every arm above, all of which cost the ``g00`` control.

    ⚠ Both bindings of ``composition_logvar`` are patched (`TRAPS.md` A10), and ``tau_lam``/``struct_lock``
    are captured from the ``build_node_init`` call that ``node_sweep`` makes immediately before."""
    orig_bni = NI.build_node_init
    orig_cl = sys.modules["rigel.calibration.enrichment_frame"].composition_logvar
    state: dict = {}

    def bni(chain, statics, geometry, **kw):
        ni = orig_bni(chain, statics, geometry, **kw)
        state["tau"] = np.asarray(ni.tau_lam, np.float64)
        state["lock"] = np.asarray(ni.struct_lock, bool)
        return ni

    def cl(f_g, E_g, E_r, var_fg, n):
        tau, lock = state.get("tau"), state.get("lock")
        v = np.asarray(var_fg, np.float64)
        if tau is not None and tau.shape == v.shape:
            free = (~lock) & (tau <= _EPS)  # no evidence, and not structurally certain
            _FIRED["zc_reference_var_slots"] = int(free.sum())
            _fire("zc_reference_var")
            v = np.where(free, 0.125, v)  # Var of Beta(½,½) — ψ's own reference
        return orig_cl(f_g, E_g, E_r, v, n)

    NI.build_node_init = bni
    bp_solver.build_node_init = bni
    sys.modules["rigel.calibration.enrichment_frame"].composition_logvar = cl
    bp_solver.composition_logvar = cl
    return orig_cl


def _discrepancy_logshift(phi_g, phi_r):
    """⭐⭐⭐ §5 OF `variance_model_notes.md` — THE DISCREPANCY INTEGRAL, as a pure function.

    ``phi_g`` is the share of the destination's OWN observed mass that the gDNA message accounts for, and
    ``phi_r`` the share the RNA arms account for. Both are observables: a measured message over a measured
    count. The rest of the mass is unexplained and could be EITHER capture enrichment of gDNA OR more RNA,
    and nothing prior-free distinguishes them.

    The enrichment ``γ`` is a SCALE, so its non-informative prior is ``1/γ`` — the same principle that gives
    ``ρ^{-1/2}`` for a rate. It is bounded above by the destination's own counts: gDNA cannot exceed the mass
    left after the RNA arms, so ``γ ≤ D`` with ``D = head/phi_g``, ``head = max(1 − phi_r, 0)``. Log-uniform
    on ``[1, D]`` means ``log γ`` is uniform on ``[0, log D]``, so::

        E[log f_g] = log(phi_g) + ½ log D          ⭐ the GEOMETRIC MEAN of the floor and the ceiling
        Var(log f_g) = (log D)² / 12               ⭐ the variance of a uniform of that width

    ⭐ **No tuned constant**: ``½`` and ``1/12`` are the mean and variance of a uniform distribution.
    ⭐ **Identically inert at ``D = 1``** — returns ``(0, 0)`` — which is the byte-identity falsification.
    ⛔ Returns a SHIFT and a VARIANCE, never a replacement, so the caller composes rather than overwrites.

    ⚠ ``head`` uses ``1 − phi_r``, not ``1``: the RNA arms may already account for part of the mass — notably
    a grafted spliced count, which is a MEASUREMENT in the destination's own frame and is not subject to
    ``γ``. Using ``1`` would hand gDNA room that is already spoken for."""
    pg = np.asarray(phi_g, np.float64)
    pr = np.asarray(phi_r, np.float64)
    head = np.maximum(1.0 - pr, 0.0)
    live = (pg > 0.0) & (head > pg)
    with np.errstate(divide="ignore", invalid="ignore"):
        logD = np.where(live, np.log(np.where(live, head / np.maximum(pg, 1e-300), 1.0)), 0.0)
    return 0.5 * logD, logD * logD / 12.0


def _install_zc_ref_prior():
    """⭐⭐⭐ **STEP 1 OF THE REBUILD** — an evidence-free slot's own belief is ψ's REFERENCE, not nothing.

    ⛔ **THE DEFECT.** ``own_composition_logvar`` returns ``∞`` whenever ``τ_λ = 0``, so ``own_precision``
    returns exactly 0, and ``_fuse(own, 0, msg, p) = msg`` for **any** ``p > 0``. The dissection found
    **100 % of the top error is `relay_only`** — own precision exactly zero — so on precisely the objects
    that carry the regression **no amount of damping can have any effect**. That is why `zc_disc_var`
    recovered 2 %: there was no own belief for a damped message to lose to.
    ⭐ And those objects are not opinionless — their ``fg_loc`` TRACKS TRUTH (measured: `worst_objects.py`
    on `g75 ss0.99 capture_on`, every one of the top 25). They hold a correct mode at a declared precision
    of zero, which makes them infinitely overridable.

    ⭐⭐⭐ **THE DERIVED VALUE, AND IT NEEDS NO NEW CONSTANT.** ``Var(log f_g) = f_r²·Var(λ)`` with
    ``Var(λ) = 1/τ``, so the only question is ``Var(λ)`` under ψ's own reference. For ``X ~ Beta(a,b)``,
    ``Var(logit X) = ψ₁(a) + ψ₁(b)``; ψ's composition reference is ``Beta(½,½)``
    (`simplex_logodds._JEFFREYS_REF`), so::

        Var(λ)|reference = 2·ψ₁(½) = π²  ≈ 9.8696        ⇒   τ_reference = 1/π² ≈ 0.1013

    ⭐ ``ψ₁(½) = π²/2`` is **already in the tree** — it is exactly what ``count_logvar(0)`` returns.

    ⭐⭐ **And the form is ADDITIVE, not a floor:** ``τ_eff = τ_λ + 1/π²``. A posterior precision is the prior's
    plus the evidence's, which is the same inverse-variance addition used everywhere else — so this is not
    `TRAPS.md` B11/D4f's refused threshold, it **deletes the ``τ = 0`` branch** instead of clipping it.
    Continuous, monotone, no discontinuity for a residue to flip.

    ⚠ ``struct_lock`` still short-circuits to 0: a structurally pure-gDNA object IS certain, and the
    reference must not weaken it — otherwise the zero-gDNA anchor loses its claim and the `g00` win with it.
    ⚠ ``has_own_composition_evidence`` is NOT touched, so the instruments' `relay_only`/`own_evidence`
    classification is unchanged and `solv%` cannot move (`TRAPS.md` B20 — admitting objects to the scored
    population is a cost that would confound the reading).
    ⚠ ``bp_solver``'s ``v_own_lam`` is a SEPARATE expression with its own ``∞`` branch inside a closure and
    is NOT reachable here, so the λ stream keeps ``∞`` while the per-component streams get the reference.
    Inconsistent, and deliberate for a prototype — it is one of the things the standalone sweep unifies."""
    ref_tau = 1.0 / (np.pi**2)  # = 1/(2·psi_1(1/2)); the Jeffreys reference's own precision on lambda
    orig = NI.own_composition_logvar

    def wrapper(f_g, tau_lam, struct_lock):
        _fire("zc_ref_prior")
        fg = np.clip(np.asarray(f_g, np.float64), 0.0, 1.0)
        fr = 1.0 - fg
        tau = np.asarray(tau_lam, np.float64) + ref_tau  # ⭐ additive: prior + evidence
        lock = np.asarray(struct_lock, bool)
        v_lam = 1.0 / tau
        return (
            np.where(lock, 0.0, fr * fr * v_lam),
            np.where(lock, 0.0, fg * fg * v_lam),
        )

    NI.own_composition_logvar = wrapper
    bp_solver.own_composition_logvar = wrapper
    return orig


def _install_zc_disc_var():
    """⭐⭐⭐ **S1a′ — THE VARIANCE ONLY, MODE UNTOUCHED.** The arm `zc_discrepancy`'s own refutation says the
    shift is the wrong half.

    ⛔ **WHY THE SHIFT IS WRONG, and it is a derivation error, not a tuning miss.** `zc_discrepancy` applies
    ``+½ log D`` and cost **+982 %** on the zero-gDNA control, with a monotone profile in gDNA amount that is
    the *mirror* of the lever it was meant to fix (+982 % at ``g00`` → −8.7 % at ``g98``). That is the
    signature of a **uniform bias toward gDNA**, not of a discrepancy-aware operator — because ``D`` is
    largest exactly where gDNA is ABSENT, so the shift is largest exactly where the message was already
    right. The geometric mean of a log-uniform on ``[1, D]`` is ``√D``, so the operator asserts an enrichment
    of ``√D`` **unconditionally**: at ``g00`` that is 1000× capture enrichment on a library with no gDNA in
    it at all.

    ⭐⭐⭐ **The error: ``γ = 1`` is not one point among many on ``[1, D]`` — it is the NULL, and the only value
    with physical backing.** Off capture gDNA is uniform and ``γ = 1`` exactly; capture makes ``γ > 1``
    *possible*, and nothing in pass-0 is evidence that it *happened*. A flat-in-log prior over ``[1, D]``
    silently converts "capture is possible" into "capture occurred, by ``√D``". `TRAPS.md` A12's shape: a
    channel whose correction tracks the answer rather than the evidence.

    ⭐ **So keep the mode at the null and charge only the width.** ``σ²_disc = (log D)²/12`` is still the
    honest statement — "I cannot rule out enrichment up to ``D``" — and it is now the *whole* statement:

        mode: unchanged (the message's own claim, i.e. γ = 1)
        variance: += (log D)²/12

    ⭐⭐ This is the owner's own proposal with my correction removed: *"the message doesn't need to change …
    the recipient converts the uncertainty into variance."* The dampening instinct was the right half.
    ⭐ And it should bite exactly where the regression is: on STRANDED data 82 % of the mass has its own
    evidence for the damped message to lose to, against 5 % unstranded (`SESSION_HANDOFF.md` §2)."""
    return _install_discrepancy(shift=False)


def _install_zc_discrepancy():
    """⭐⭐⭐ **S1a** — the geometric midpoint and its variance, applied to the delivered gDNA share.

    ⛔ **THE PATCH POINT IS CHOSEN SO THAT NOTHING IN `src/` MOVES** (`TRAPS.md` B18). The operator belongs
    at `bp_solver.py:1511`, inside `node_sweep`'s closure and therefore unreachable — but everything it needs
    is reconstructible one call later, at `simplex_logodds._solve_nodes_logodds_all`, which is module-level:

        gdna_imp_mode = log(cg·E_g/M) = log φ_g^msg      ⭐ ALREADY the delivered gDNA SHARE
        rna_imp_mode  = (log φ_p, log φ_n)

    so ``φ_g = exp(gdna_imp_mode)`` and ``φ_r = exp(mo_p) + exp(mo_n)`` — the two arguments §5 needs, with no
    access to any belief. **D4 holds**: the shift is built from the message and the destination's own count.

    ⚠ Patched on ``bp_solver``'s binding, which is a SEPARATE module global from the leaf module's
    (`TRAPS.md` A10). ``node_init`` and ``node_geometry`` also import the name, but their calls pass no
    imputation modes, so the operator cannot fire there — and the counter proves it rather than the reader
    assuming it.

    ⚠ Scoped to the **gDNA arm**. Every object in the regression has gDNA UNDER-called, and the dominant
    partial-message case is ``{g} → {g, R±}`` where the shared population set is exactly ``{g}``. The general
    multi-component shared case needs the population predicate and is S2."""
    return _install_discrepancy(shift=True)


def _install_discrepancy(*, shift: bool):
    """The shared installer. ``shift`` selects the mode treatment; the variance is always applied, so the
    two arms differ in exactly one term and nothing else."""
    SL = sys.modules["rigel.calibration.simplex_logodds"]
    orig = SL._solve_nodes_logodds_all
    tag = "zc_discrepancy" if shift else "zc_disc_var"

    def wrapper(*a, gdna_imp_mode=None, gdna_imp_prec=None, rna_imp_mode=None, **kw):
        if gdna_imp_mode is not None:
            mg = np.array(gdna_imp_mode, np.float64)
            phi_g = np.exp(mg)
            phi_r = np.zeros_like(phi_g)
            if rna_imp_mode is not None:
                for arm in rna_imp_mode:
                    phi_r = phi_r + np.exp(np.asarray(arm, np.float64))
            sh, var = _discrepancy_logshift(phi_g, phi_r)
            _FIRED[tag + "_slots"] = int((var > 0.0).sum())
            _fire(tag)
            gdna_imp_mode = mg + sh if shift else mg
            if gdna_imp_prec is not None:
                p = np.asarray(gdna_imp_prec, np.float64)
                gdna_imp_prec = np.where(p > 0.0, p / (1.0 + p * var), 0.0)
        return orig(
            *a,
            gdna_imp_mode=gdna_imp_mode,
            gdna_imp_prec=gdna_imp_prec,
            rna_imp_mode=rna_imp_mode,
            **kw,
        )

    SL._solve_nodes_logodds_all = wrapper
    bp_solver._solve_nodes_logodds_all = wrapper
    return orig


def _install_zc_struct_lock_g1():
    """⭐⭐⭐ ``struct_lock`` MEANS "STRUCTURALLY PURE gDNA" IN ITS OWN DOCSTRING AND MEANS "UNSOLVABLE"
    IN ITS CODE. This arm makes the mask match the prose: ``g1_locked(free_pos, free_neg) ∧ NODE``.

    ``strand_evidence`` is handed ``locked = ~solvable`` with
    ``solvable = (free_pos | free_neg) & (n_node > 0)``, so ``struct_lock = locked & is_region`` is true at
    **every NODE with zero counts** — an empty EXON and an empty INTRON included, not only at the "true
    intergenic NODE nodes" the docstring scopes it to. ``own_composition_logvar`` then returns
    ``Var(log f_g) = 0`` there: **composition CERTAIN**, on a slot whose ``f_g`` is a default belief and
    whose evidence is nothing at all.

    ⚠ It was INERT until 2026-08-06: ``own_precision``'s ``n > 0`` gate silenced every zero-count slot, so
    a certainty granted to 18,511 empty nodes could not leave them. Removing that gate un-masked it.
    ⭐ ``g1_locked`` is the ONE HOME for the predicate (`TRAPS.md` A11) and already exists; this arm simply
    calls it. At ``g00`` every intergenic anchor is G1, so the arm must KEEP the zero-gDNA win — that is
    what makes it different from muting the anchors."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        _fire("zc_struct_lock_g1")
        fp = np.asarray(statics.free_pos, bool)
        fn = np.asarray(statics.free_neg, bool)
        is_reg = np.asarray(chain.kind) == NODE
        want = g1_locked(fp, fn) & is_reg
        have = np.asarray(ni.struct_lock, bool)
        _FIRED["struct_lock_slots_HEAD"] = int(have.sum())
        _FIRED["struct_lock_slots_G1"] = int(want.sum())
        _FIRED["struct_lock_slots_REVOKED"] = int((have & ~want).sum())
        _FIRED["struct_lock_slots_ADDED"] = int((want & ~have).sum())
        M = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
        E_g = np.asarray(node_global_geometry(geometry)[1], np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        f_g = np.asarray(ni.f_g, np.float64)
        v_fg, v_fr = NI.own_composition_logvar(f_g, ni.tau_lam, want)
        prec_g = NI.own_precision(f_g * M, v_fg, E_g > 0.0)
        prec_pos = NI.own_precision(
            np.asarray(ni.rho_pos, np.float64) * E_r, v_fr, fp & (E_r > 0.0)
        )
        prec_neg = NI.own_precision(
            np.asarray(ni.rho_neg, np.float64) * E_r, v_fr, fn & (E_r > 0.0)
        )
        return dataclasses.replace(
            ni, struct_lock=want, prec_g=prec_g, prec_pos=prec_pos, prec_neg=prec_neg
        )

    NI.build_node_init = wrapper
    bp_solver.build_node_init = wrapper
    return orig


def _install_zc_logmean():
    """⭐⭐⭐ THE OTHER CONSISTENT LOCATION — ``exp(digamma(a+½))/E``, the **geometric** mean of the same
    ``Gamma(a + ½, E)`` posterior, at a zero-mass slot.

    ⛔ **THIS IS THE LOCATION THE EXISTING VARIANCE ACTUALLY BELONGS TO, and that is the argument for it
    over the arithmetic mean.** The relay carries a pair — a density and a precision on its LOG
    (``own_precision`` = ``1/(Var(log f) + Var(log ρ))``, ``count_logvar`` = ``trigamma(a+½)`` = the
    variance of ``log ρ`` under that posterior). A log-space variance describes a log-space location, and
    that location is ``E[log ρ] = digamma(a+½) − log E``. HEAD pairs ``Var(log ρ) = trigamma(½) = 4.93``
    with a location of **exactly 0**, whose log is ``−∞``: the two halves of the message describe
    different distributions, and a *relative* spread around zero is not a statement about anything.

    At ``a = 0`` this is ``exp(digamma(½))/E = 0.1404/E`` — five-eighths of a decade below the arithmetic
    mean ``0.5/E``, and both → ``a/E`` for ``a ≥ 10``, so like :func:`count_logvar` the whole change is at
    small counts. Scoped to ``M = 0`` for the same reason as ``zc_jeffreys_mean``: that is where
    ``Σ_c ρ_c·E_c = M`` places no constraint."""
    from scipy.special import digamma

    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        sel, _M = _zero_mass_locked(ni, geometry)
        _FIRED["zc_logmean_slots"] = int(sel.sum())
        _fire("zc_logmean")
        E_g = np.asarray(node_global_geometry(geometry)[1], np.float64)
        rho_g = np.array(ni.rho_g, np.float64)
        ok = sel & (E_g > 0.0)
        rho_g[ok] = float(np.exp(digamma(0.5))) / E_g[ok]
        return dataclasses.replace(ni, rho_g=rho_g)

    NI.build_node_init = wrapper
    bp_solver.build_node_init = wrapper
    return orig


def _targets(chain, region_arrays):
    """Vectorised: every `intron|exon` EDGE slot and the slot of its flanking INTRON NODE.

    ⭐ Genomic order IS slot order and the chain strictly alternates NODE/EDGE, so an EDGE's flanks are
    ``i-1`` and ``i+1`` and no traversal is needed. Breaks at a reference terminal are handled by
    `chain.left`/`chain.right` being −1 there."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_node = kind == NODE
    ntype = np.where(is_node, rtype[np.clip(obj, 0, rtype.shape[0] - 1)], -1)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    ok = (~is_node) & (left >= 0) & (right >= 0)
    lt = np.where(ok, ntype[np.clip(left, 0, ntype.size - 1)], -1)
    rt = np.where(ok, ntype[np.clip(right, 0, ntype.size - 1)], -1)
    hit = ok & (((lt == INTRON) & (rt == EXON)) | ((lt == EXON) & (rt == INTRON)))
    edges = np.flatnonzero(hit)
    srcs = np.where(lt[edges] == INTRON, left[edges], right[edges])
    return edges, srcs


def _wrap_node_sweep():
    """Stash `region_arrays` — `node_sweep` receives it and calls `build_node_init` after."""
    orig = CAL.node_sweep

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _RA["region_arrays"] = region_arrays
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.node_sweep = wrapper
    return orig


def _install_onesided_rna():
    """⭐⭐⭐ **THE CERTIFIED-RNA CLAIM IS A LOWER BOUND. DELIVER IT AS ONE.**

    ``bp_solver.py:439-447`` states the premise in the source's own words: *"what the graft actually knows
    is an INEQUALITY, ``rho_R(exon) >= rho_nu(B) + rho_mu(B)`` … and it uses it as an equality."* ψ then
    applies ``-1/2 p (log f_active - mo_p)^2`` in BOTH code paths — symmetric in the residual, no hinge
    anywhere in the file — so a destination holding MORE RNA than the bound, which the inequality
    explicitly permits, is penalised exactly as hard as one holding less.

    ⛔ **And an UNDER-claiming RNA message therefore drives ``f_g`` UP.** On an unstranded 1 %-gDNA library
    an exon's mass is essentially all RNA (``f_active`` near 1), so a message saying "the RNA share is 0.3"
    is read as "70 % of this is gDNA". Measured at ``g01 ss0.50 capture_on`` by ψ-boundary ablation with
    the A5 identity exact: HEAD's self-solve is 0.0086 against a truth of 0.0023, the message drives it to
    **0.3219** at precision 327, and muting the channel returns **0.0005**. Eight of HEAD's twelve worst
    slots behave that way.

    ⭐ **Why this is not the twelfth refused candidate.** Every one of the eleven in the graveyard was a
    rule for how to resolve DOUBT, and each was refused by the ``g00`` control because there the doubt must
    resolve to *no* gDNA. This adds doubt in **one direction only** — "at least this much RNA" — so at
    ``g00``, where the truth is ``f_active = 1``, every bound is satisfied and the channel goes **inert**
    rather than harmful. ⭐ And the message-precision sweep says the defect is a BIAS rather than a
    loudness (no plateau; the stranded optimum is exactly zero), which is `TRAPS.md` D1: a variance cannot
    fix a bias, and all three operators currently pricing this inequality are variances.

    ⚠ **Only the RNA channel.** The gDNA measurement is arguably a lower bound too — under capture
    ``gamma >= 1``, so the destination may hold MORE gDNA than a neighbour reports — but that is a second
    thing varied and a separate arm. ⚠ The bound also reaches ψ a SECOND time through ``tlam`` ->
    ``lam_imp`` (``bp_solver.py:1474-1479``), also two-sided; this arm does not touch that path, so a
    partial recovery is the expected shape rather than a full one.

    ⛔ The switch is `simplex_logodds.ONE_SIDED_RNA`, default False and byte-identical off — with the flag
    down `_rna_residual` returns its input difference unmodified, which is why ``--arm base`` must still
    reproduce the pre-refactor panel exactly (A5).
    """
    import rigel.calibration.simplex_logodds as SL

    SL.ONE_SIDED_RNA[0] = True
    _fire("onesided_rna")


def _install_msgscale(scale: float):
    """⭐⭐⭐ **DO MESSAGES ONLY NEED TO BE WEAK?** — the owner's hypothesis, as a one-parameter sweep.

    The owner's account, from having shipped the production tool: *"messages do not need to be confident.
    When messages become confident they overwrite the strand-specific data, and so they harm scenarios
    with strand specificity. When we have unstranded data none of the nodes has a solution — all have a
    precision of zero — so the weakest of weak messages is still going to work, because it is more
    precision than zero."*

    That is a quantitative claim and it has a shape: **multiply every message precision by one scalar and
    sweep it.** ``scale = 1`` is ``base``; ``scale = 0`` is ``msgfree_all``, which is already measured. If
    the account is right there is a PLATEAU at small scale where the stranded strata sit at their
    message-free optimum (because the messages have become negligible against real strand evidence) and
    unstranded x capture-ON still works (because any positive precision beats zero).

    ⛔⛔ **THIS IS A DIAGNOSTIC, NOT A PROPOSED FIX.** A global damping scalar is a tuned constant and
    `CLAUDE.md` G1 forbids one. What the sweep decides is which KIND of defect this is:

      * a PLATEAU exists  => the precision model is wrong by a roughly constant factor, the fix is a
                             derivation that produces that factor, and the backbone keeps the message
                             layer's structure;
      * NO plateau        => no single attenuation serves both regimes, the defect is structural rather
                             than a mis-calibration, and the backbone must change the messages themselves.

    ⚠ All four channels are scaled together — ``gdna_imp``, ``rna_imp``, ``lam_imp``, ``theta_imp`` — which
    is ONE thing varied and is exactly the claim as stated. Scaling them separately is a different, later
    experiment. The MODES are untouched, so nothing about what a message SAYS changes; only how loudly.
    ⛔ A10 — all three ``_solve_nodes_logodds_all`` bindings are patched and the arm raises if unfired.
    """
    import rigel.calibration.simplex_logodds as SL

    orig_solve = SL._solve_nodes_logodds_all
    k = float(scale)

    def _s(p):
        if p is None:
            return None
        if isinstance(p, tuple):
            return tuple(_s(q) for q in p)
        return np.asarray(p, np.float64) * k

    def solve(*a, **kw):
        _fire(f"msgscale_{scale:g}")
        for key in ("gdna_imp_prec", "rna_imp_prec", "lam_imp_prec", "theta_imp_prec"):
            if kw.get(key) is not None:
                kw[key] = _s(kw[key])
        return orig_solve(*a, **kw)

    for mod in (SL, NI, bp_solver):
        if hasattr(mod, "_solve_nodes_logodds_all"):
            mod._solve_nodes_logodds_all = solve


def _install_msgfree(where: str):
    """⭐⭐⭐ **HOW MUCH OF THE MESSAGE LAYER DOES THE SUBSTRATE ACTUALLY NEED?** — the consolidation arm.

    ⛔ **Pass-0's ONLY job is to be a training substrate for the gDNA hyperprior** (owner, 2026-08-07;
    ``variance_model_notes.md`` A2). It is not required to be accurate, and it is not the deliverable. So
    the question that sizes the whole backbone has never been asked: **if pass-0 sends no messages at all,
    is the fitted prior worse?**

    Two arms, because they answer two different questions and conflating them would be one experiment
    pretending to be two:

    ``msgfree_p0``   psi's four imputed channels are muted during the PRIOR-FREE sweep ONLY; every refit
                     sweep runs untouched. ⭐ This is the one that sizes the SUBSTRATE: it asks whether the
                     prior needs message-imputed values at the objects that have no evidence of their own.
    ``msgfree_all``  muted in every sweep. This asks whether the DELIVERABLE needs messages at all, and it
                     is the floor: whatever it scores is what the accumulator plus psi plus the fitted
                     prior are worth with no belief propagation whatsoever.

    ⚠ **What is muted is exactly the four ``*_imp`` arguments** — ``gdna_imp``, ``rna_imp``, ``lam_imp``,
    ``theta_imp``. Everything else psi receives is the slot's OWN evidence (its two strand counts, its
    spliced count, the reference, the fitted prior, the intron factory) and is untouched, so a difference
    is attributable to the message layer and to nothing else. ⭐ The whole relay still RUNS; only its
    delivery into psi is withheld. That keeps one thing varied and leaves the geometry identical.

    ⛔ A10 — ``_solve_nodes_logodds_all`` is bound as a module global in THREE places
    (``simplex_logodds``, ``node_init``, ``bp_solver``); all three are patched and the arm raises if it
    never fired. And the pass-0-vs-refit switch is read off ``gdna_prior is None``, which is exactly how
    ``calibrate`` itself distinguishes the two phases (``calibrate.py:528`` vs ``:572``).
    """
    import rigel.calibration.simplex_logodds as SL

    state = {"muted": where == "all"}
    orig_solve = SL._solve_nodes_logodds_all

    def solve(*a, **kw):
        if state["muted"]:
            _fire(f"msgfree_{where}")
            for k in ("gdna_imp_mode", "gdna_imp_prec", "rna_imp_mode", "rna_imp_prec",
                      "lam_imp_mode", "lam_imp_prec", "theta_imp_mode", "theta_imp_prec"):
                kw[k] = None
        return orig_solve(*a, **kw)

    for mod in (SL, NI, bp_solver):
        if hasattr(mod, "_solve_nodes_logodds_all"):
            mod._solve_nodes_logodds_all = solve

    if where == "p0":
        orig_sweep = bp_solver.node_sweep

        def sweep(chain, statics, geometry, belief, region_arrays, *a, **kw):
            # the prior-free sweep is the one `calibrate` calls with `gdna_prior=None`
            state["muted"] = kw.get("gdna_prior") is None
            try:
                return orig_sweep(chain, statics, geometry, belief, region_arrays, *a, **kw)
            finally:
                state["muted"] = False

        for mod in (CAL, bp_solver):
            if hasattr(mod, "node_sweep"):
                mod.node_sweep = sweep


def _install_eta(mute_level: bool = False):
    """⭐⭐⭐ THE ``eta`` REBUILD — the whole sweep replaced, not one term overridden.

    Every other arm in this file overrides a single expression inside the shipped sweep. This one swaps
    the sweep itself for `eta_node_sweep.eta_node_sweep`, which is the point: the derivation deletes the
    mass pin, the reframe ``r``, ``framed``, the flank pair and the graft/peel gating TOGETHER, and a
    sequence of partial edits would leave each new piece fighting the operators that remain.

    ⛔ A10 — ``node_sweep`` is bound as a module global in ``calibrate`` (``from .bp_solver import
    node_sweep``) as well as in ``bp_solver`` itself. BOTH are patched, and the arm asserts it fired: an
    arm that never ran scores byte-identical to base and reads as a clean "the rebuild is neutral", which
    is exactly the false negative this campaign must not publish.
    """
    ETA = _sibling("eta_node_sweep.py")
    ETA.MUTE_LEVEL[0] = bool(mute_level)

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _RA["region_arrays"] = region_arrays
        _fire("eta_nolevel" if mute_level else "eta")
        return ETA.eta_node_sweep(chain, statics, geometry, belief, region_arrays, *a, **kw)

    for mod in (CAL, bp_solver):
        if hasattr(mod, "node_sweep"):
            mod.node_sweep = wrapper


def _install_face_one():
    """⭐ FACE (I) AS DERIVED (`EQUATIONS.md` §3.6): the `intron|exon` EDGE takes the flanking INTRON's
    COMPOSITION — never its level — and closes the level with its OWN observed mass.

    ``T(INTRON) = {gDNA, nascent} = T(EDGE, unspliced only)``, so the two measure the same population and
    the share may cross. ⛔ The LEVEL may not: measured on the toy, transferring the intron's gDNA
    DENSITY instead is **+0.017 panel-wide and +0.207 on capture-ON × unstranded**, because capture
    depletes the intron ~1000× while enriching the EDGE."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        ra = _RA.get("region_arrays")
        if ra is None:
            return ni
        edges, srcs = _targets(chain, ra)
        if edges.size == 0:
            return ni
        f_g = np.array(ni.f_g, np.float64)
        f_pos = np.array(ni.f_pos, np.float64)
        f_neg = np.array(ni.f_neg, np.float64)
        tau = np.array(ni.tau_lam, np.float64)
        lock = np.asarray(ni.struct_lock, bool)
        M, E_g = node_global_geometry(geometry)
        M = np.asarray(M, np.float64)
        E_g = np.asarray(E_g, np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        n_node = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)

        new_fg = f_g[srcs]
        rna = np.maximum(0.0, 1.0 - new_fg)
        tot = f_pos[edges] + f_neg[edges]
        fp_ok = np.asarray(statics.free_pos, bool)[edges]
        fn_ok = np.asarray(statics.free_neg, bool)[edges]
        k = fp_ok.astype(np.float64) + fn_ok.astype(np.float64)
        share_p = np.where(tot > _EPS, f_pos[edges] / np.maximum(tot, _EPS),
                           np.where(k > 0, fp_ok / np.maximum(k, 1.0), 0.0))
        share_n = np.where(tot > _EPS, f_neg[edges] / np.maximum(tot, _EPS),
                           np.where(k > 0, fn_ok / np.maximum(k, 1.0), 0.0))
        f_pos[edges] = rna * share_p
        f_neg[edges] = rna * share_n
        f_g[edges] = new_fg
        tau[edges] = tau[srcs]

        v_fg, v_fr = NI.own_composition_logvar(f_g, tau, lock)
        rho_g = np.maximum(
            np.where((M > _EPS) & (E_g > _EPS), f_g * M / np.maximum(E_g, _EPS), 0.0), 0.0
        )
        prec_g = NI.own_precision(n_node, v_fg, rho_g > _EPS)
        touched = np.zeros(f_g.shape[0], bool)
        touched[edges] = True

        def _rna(frac, free_s, rho_old):
            raw = np.where(
                (M > _EPS) & (E_r > _EPS) & np.asarray(free_s, bool),
                frac * M / np.maximum(E_r, _EPS), 0.0,
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--arm",
        choices=(
            "base",
            "eta",
            "eta_nolevel",
            "msgfree_p0",
            "msgfree_all",
            "msgscale_0.001",
            "msgscale_0.01",
            "msgscale_0.1",
            "msgscale_0.5",
            "onesided_rna",
            "intron_phi",
            "kappa_half",
            "zc_noop",
            "zc_own_count",
            "zc_total_n",
            "zc_live_count",
            "zc_transfer",
            "zc_anchor_mute",
            "zc_jeffreys_mean",
            "zc_logmean",
            "zc_struct_lock_g1",
            "zc_reference_var",
            "zc_discrepancy",
            "zc_disc_var",
            "zc_ref_prior",
            "zc_ref_prior_damp",
        ),
        required=True,
    )
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_ladder_ceiling"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--jobs", type=int, default=1,
        help="run this many conditions CONCURRENTLY, by re-invoking this script on shards and "
             "concatenating. The conditions are independent, so this changes no number.",
    )
    args = ap.parse_args()

    # ── ⭐⭐ SHARDED PARALLELISM. Conditions are completely independent — separate BAMs, separate
    # calibrations, nothing shared but the read-only index and oracle cache — so this is a pure
    # wall-clock win with NO effect on any measurement. It re-invokes the SAME single-process path on
    # a subset rather than parallelising inside it, so the arm-installation logic every `zc_*` and the
    # `eta` arm depends on is byte-for-byte the code that was verified serially.
    # ⚠ ``OMP_NUM_THREADS=1`` is already forced at import, so the workers do not fight over threads.
    if args.jobs > 1:
        names = args.conditions or sorted(
            p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
        )
        shards = [names[i :: args.jobs] for i in range(args.jobs)]
        shards = [s for s in shards if s]
        tmp = args.out.parent / f".{args.out.stem}_shards"
        tmp.mkdir(parents=True, exist_ok=True)
        procs, outs = [], []
        for i, sh in enumerate(shards):
            o = tmp / f"{i}.jsonl"
            outs.append(o)
            cmd = [sys.executable, str(Path(__file__).resolve()), "--arm", args.arm,
                   "--suite", str(args.suite), "--index", str(args.index),
                   "--work-dir", str(args.work_dir / f"shard{i}"), "--out", str(o),
                   "--conditions", *sh]
            if args.oracle_cache is not None:
                cmd += ["--oracle-cache", str(args.oracle_cache)]
            procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                          text=True))
        rc = 0
        for i, pr in enumerate(procs):
            out, _ = pr.communicate()
            if pr.returncode != 0:
                rc = pr.returncode
                print(f"  ⛔ shard {i} FAILED (rc={pr.returncode}):\n{out}", flush=True)
            else:
                print(f"  shard {i}: {len(shards[i])} conditions ok", flush=True)
        if rc:
            # ⛔ A10's shape: a shard that died leaves a SHORT output file, and concatenating it
            # silently would publish a partial panel that reads like a complete one.
            raise RuntimeError(f"{args.arm}: a shard failed; refusing to concatenate a partial panel")
        with args.out.open("w") as fh:
            for o in outs:
                fh.write(o.read_text())
        n = sum(1 for _ in args.out.open())
        print(f"  ⭐ {args.arm}: {n} rows from {len(shards)} shards -> {args.out}", flush=True)
        return 0

    # ⚠ the eta arm REPLACES the sweep, so it installs its own region_arrays stash rather than wrapping
    #   the shipped one — wrapping first would leave the base sweep in the chain.
    if args.arm in ("eta", "eta_nolevel"):
        _install_eta(mute_level=args.arm == "eta_nolevel")
    elif args.arm in ("msgfree_p0", "msgfree_all"):
        _install_msgfree("p0" if args.arm == "msgfree_p0" else "all")
    elif args.arm.startswith("msgscale_"):
        _install_msgscale(float(args.arm.split("_", 1)[1]))
    elif args.arm == "onesided_rna":
        _install_onesided_rna()
    else:
        _wrap_node_sweep()
    if args.arm == "intron_phi":
        _install_face_one()
    elif args.arm == "kappa_half":
        _install_kappa_half()
    # ⭐ ``zc_noop`` re-derives HEAD's own two decisions through the SAME rebuild path as the three
    #   reverts, so it must come out byte-identical to ``base``. That is the falsification for the
    #   rebuild itself (`TRAPS.md` A5) — if noop ≠ base, no zc_* reading means anything.
    elif args.arm == "zc_noop":
        _install_rebuild("zc_noop", total_n=False, count_live=False)
    elif args.arm == "zc_own_count":
        _install_zc_own_count()
    elif args.arm == "zc_total_n":
        _install_rebuild("zc_total_n", total_n=True, count_live=False)
    elif args.arm == "zc_live_count":
        _install_rebuild("zc_live_count", total_n=False, count_live=True)
    elif args.arm == "zc_transfer":
        _install_zc_transfer()
    elif args.arm == "zc_anchor_mute":
        _install_zc_anchor_mute()
    elif args.arm == "zc_jeffreys_mean":
        _install_zc_jeffreys_mean()
    elif args.arm == "zc_logmean":
        _install_zc_logmean()
    elif args.arm == "zc_struct_lock_g1":
        _install_zc_struct_lock_g1()
    elif args.arm == "zc_reference_var":
        _install_zc_reference_var()
    elif args.arm == "zc_discrepancy":
        _install_zc_discrepancy()
    elif args.arm == "zc_disc_var":
        _install_zc_disc_var()
    elif args.arm == "zc_ref_prior":
        _install_zc_ref_prior()
    elif args.arm == "zc_ref_prior_damp":
        # ⭐ the PAIR: the own belief exists (step 1) AND the message is damped by the discrepancy
        #   (`zc_disc_var`). Neither can work alone — damping needs something to lose to, and an own
        #   belief without damping just competes with an undamped message. `TRAPS.md` D4j.
        _install_zc_ref_prior()
        _install_zc_disc_var()

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig()
    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    with args.out.open("w") as fh:
        for name in names:
            t0 = time.perf_counter()
            cond = args.suite / name
            truth = P0.truth_f_gdna(cond) or 0.0
            m = P0.measure_condition(
                bam=str(cond / "sim_oracle.bam"), index=index, pipeline_config=PipelineConfig(),
                calibration_config=config, work_dir=args.work_dir / "rigel_pass0_oracle", tag=name,
                # ⛔ NO truth pmfs. They exist only to build the two ``c_input_*`` arms, and this
                # harness reads ``pass0`` and ``final`` ONLY — so those two calibrate runs were pure
                # waste. Measured: 35.2 s -> 24.5 s per condition (**-30 %**), 10 node_sweep calls -> 5,
                # and all four scored fields BYTE-IDENTICAL on both axes. ⭐ Verified, not assumed
                # (`TRAPS.md` A5): dropping work that changes a number is a different change entirely.
                truth_pmfs=None,
                oracle_cache=args.oracle_cache,
            )
            for axis in ("node", "edge"):
                s = SA.summarise(SA.audit(m, axis=axis, config=config))
                # ⛔⛔ FIXED-DENOMINATOR COMPANIONS. `solvable_mwae`'s denominator is the SOLVABLE set,
                # and an arm that changes what counts as solvable changes its own denominator — so a
                # fall in it can be "solved better" or "stopped scoring the hard ones" and the number
                # cannot tell you which (`TRAPS.md` B13/B20). These three cannot move that way:
                #   * `mwae_all`  — every object with mass, so the denominator is the library.
                #   * `abs_err_all` — the raw error MASS, no denominator at all.
                #   * `library_f_gdna` — the deliverable itself, one number per condition.
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
                fh.write(json.dumps({"arm": args.arm, "condition": name, "axis": axis,
                                     "f_gdna": truth, **s}) + "\n")
                fh.flush()
            print(f"  {name} {time.perf_counter() - t0:.0f}s", flush=True)
    # ⛔⛔ A10 — an ablation that never ran scores byte-identical to base and reads as "inert".
    # ⚠ a COMPOSITE arm fires only its components' names, so require ANY zc_* firing rather than this
    #   arm's own name. The guard tripped on `zc_ref_prior_damp` AFTER a complete, valid run — it caught a
    #   bookkeeping gap, not a measurement one, and narrowing it here keeps A10's teeth for the real case
    #   (nothing fired at all).
    if (args.arm.startswith("eta") and not _FIRED.get(args.arm)) or (
        args.arm.startswith("zc_")
        and not any(k.startswith("zc_") and not k.endswith("_slots") for k in _FIRED)
    ):
        raise RuntimeError(
            f"arm {args.arm!r} NEVER FIRED — the patched name is not the one the solver calls. "
            f"fired: {_FIRED or '{}'} (TRAPS.md A10)"
        )
    if _FIRED:
        print(f"  ⭐ ablation fired {_FIRED} time(s) — {args.arm}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
