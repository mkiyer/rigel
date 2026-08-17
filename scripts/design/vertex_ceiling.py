#!/usr/bin/env python
"""⭐⭐⭐ WHAT IS PERFECTING THE SIMPLEX VERTEX WORTH? — the re-solve ceiling, on the REAL 36-condition
ladder, plus the mechanism prototype in the same harness so the two are directly comparable.

⛔⛔ **THIS IS THE MEASUREMENT THAT DECIDES BUILD-VS-NOTE** (TRAPS: measure-the-ceiling-first, which has re-ranked this
project five times). A silent gene's objects are pure gDNA and the truth is ``f_g = 1.000`` exactly; a
zero-gDNA library's objects are pure RNA and the truth is ``0.000`` exactly. The solver lands 2–8 % inside
both. Before designing anything, hand those objects the exact answer, **re-solve the whole chain**, and
read what the headroom actually is.

⭐ **A RE-SOLVE, NOT A SUBSTITUTION** (TRAPS: substitution-understates-a-source). A vertex object is mostly a message SOURCE — its
value is what it CARRIES — and substituting its answer after the fact does not propagate. So the pin goes
in at ``region_init.build_region_init``, the per-object message-free self-solve, and the relay then runs on
top of it. That is also the one place a mechanism can be expressed without touching either relay twin,
so a prototype cannot be gated in one twin and not the other (TRAPS: name-the-observable-per-site).

**THE ARMS.**

| arm | what it does | what it is |
|---|---|---|
| ``base`` | nothing | the baseline, re-recorded from the current tree in the same run (TRAPS: re-record-the-baseline) |
| ``noop`` | pins the truth at ZERO objects | ⭐ the harness's own falsification — MUST be byte-identical to ``base`` |
| ``vertex_free`` | pins oracle truth at every vertex-truth object with **no own composition evidence** | ⭐⭐ **THE CEILING.** The population a vertex fix can reach |
| ``vertex_all`` | pins oracle truth at **every** vertex-truth object | a looser upper bound — includes objects that have their own evidence |
| ``ref_c=<C>`` | sets ψ's reference exponent to ``C`` instead of ½ | ⭐ the mechanism prototype (TRAPS: panel-before-src — the panel arm before ``src/``) |
| ``ref_loc=noop`` | the full location path, taking NO location | ⭐ byte-identical to ``base`` on 16/16 conditions — the family's own falsification |
| ⭐ ``ref_loc=struct`` | ψ's reference MEAN pinned near 1 wherever mature RNA cannot be | **the candidate.** Final err-sum per stratum 0.381 / 0.659 / 0.363 / 0.800, control 1.000 |
| ``ref_loc=struct_grid`` | the same, capped at the highest grid point ``sigma(L)`` | ⭐⭐ the DERIVED strength — exactly ``L`` nats, no constant. Reproduces ``struct`` to 3 dp (0.384 / 0.660 / 0.366 / 0.800) |
| ``ref_loc=struct_soft`` | the same, floored at one pseudo-fragment AT THE OBJECT | ⛔ REFUSED — worse on every stratum (0.609 / 1.045 / 0.580 / 1.000) |
| ⛔ ``ref_loc={pooled,local}`` | the peel ``rho_g*E_g/M`` on the non-structural slots too | ⛔ REFUSED — **5.3x WORSE** on stranded × capture-ON, because the on-target gDNA density is under-read 2.6-3.6x pre-solve |
| ⭐⭐⭐ ``config_struct`` | ``CalibrationConfig(structural_reference=True)`` and **NO monkeypatch at all** | the same claim THROUGH ``calibrate``. ⛔ It is not a duplicate of ``ref_loc=struct_grid`` — that one injects at ψ's construction site and cannot reach `region_init`'s ``tau_lam``, which the shipped path now feeds too, so the two are DIFFERENT MEASUREMENTS and only this one is the deliverable |

⚠ ``vertex_free``'s "no own evidence" test is ``tau_lam <= 1e-4``, and that is a CLASSIFICATION FOR A
CEILING, never a production predicate — TRAPS: a-threshold-on-a-fitted-residue refused exactly this shape as a mechanism. Both
arms are reported side by side so the filter's effect is visible rather than assumed.

⛔ **TRAPS: honesty-metrics-reward-ignorance — the honesty columns are never quoted alone.** Every row carries ``mwae_all`` (denominator =
every object with mass) and ``abs_err_all`` (no denominator at all) beside the solvable-set numbers,
because an arm that changes what counts as solvable changes its own denominator.

⛔ **TRAPS: an-ablation-that-never-ran — every arm counts its own firings and RAISES if it did not fire.** An override that never ran
reads as "no effect", which is publishable and wrong.

⚠ **TRAPS: could-the-arm-have-fired — the pin count per condition is printed.** An arm with zero opportunities is not a control.

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
objects.

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
region       0.1247  **0.0975**    −149,267        27 / 9
boundary       0.1434  **0.1127**    −161,302        29 / 3
=========  ======  ============  ==============  ==============

⭐⭐⭐ **AND IT SPLITS ON EXACTLY ONE AXIS — STRAND — which is what the mechanism predicts.** Pass-0
``mwae`` delta, region axis: unstranded **−0.0188** (capture off, 9/0) and **−0.0963** (capture ON, 9/0);
stranded **−0.0003** and **+0.0064** (2/7). Every one of the 9 rows that got worse is ``ss_0.99``.
The strand channel's Fisher information is ``∝ (2κ−1)²`` and is EXACTLY zero at κ = ½, so on an
unstranded library ψ's reference is the only term left with a gradient at the vertex, while on a stranded
one the strand term supplies the λ information and the vertex is already reached. ⭐ The ceiling is
therefore entirely on unstranded data — which is also the panel's worst stratum
(``capture_ON × ss0.50``: base 0.3235 region / 0.2922 boundary). Largest single row:
``gdna_g01_ss_0.50_capture_on``, **−0.2188**.

⛔⛔ **TWO WARNINGS THAT MUST TRAVEL WITH THE NUMBER.**

* **``vertex_all`` is WORSE than ``vertex_free``** — region pass-0 0.1076 vs 0.0975, and on the shipped
  solve it is worse than *base* (+0.0080). Pinning MORE truth hurts: the extra objects have their own
  evidence, and declaring them certain overrides it and propagates. That is TRAPS: admitting-an-object-costs's shape
  reached with the TRUTH, so the harm is in the relay's dynamics and not in the values. ⭐ Quote
  ``vertex_free``, and note that a fix which hands out certainty broadly can lose even when it is right.
* **The honesty columns move the WRONG way** — confidently-wrong Σ\\|err\\| +9,175 (region) / +893 (boundary),
  28 and 16 rows worse — while accuracy improves 22 %. TRAPS: honesty-metrics-reward-ignorance exactly: certainty handed to an
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

from rigel.calibration import region_init as NI, sweep as SW  # noqa: E402
from rigel.calibration import simplex_logodds as SL  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_geometry import region_gdna_geometry  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

CAL = sys.modules["rigel.calibration.calibrate"]

_EPS = 1.0e-9
#: the ceiling's own classification of "this object has no composition evidence of its own". ⛔ NOT a
#: production predicate (TRAPS: a-threshold-on-a-fitted-residue); the `vertex_all` arm exists so its effect is measured, not
#: assumed.
_TAU_FREE = 1.0e-4

#: filled by the wrappers, one call before `build_region_init` needs them.
_CTX: dict = {}
#: TRAPS: an-ablation-that-never-ran — per-arm firing counters. A zero here RAISES.
_FIRED: dict = {"init": 0, "pinned": 0, "ref_g": 0, "ref_r": 0, "loc": 0, "conditions": 0}


# ── the plumbing: get the oracle's per-object truth and the geometry to `build_region_init` ────────────


def _wrap_oracle():
    """Stash the origin-split oracle. It is built BEFORE any arm calibrates, so its per-object truth
    is available to the self-solve — which is what makes a re-solve ceiling possible at all."""
    orig = P0.load_or_build_oracle

    def wrapper(*a, **kw):
        oracle = orig(*a, **kw)
        _CTX["oracle"] = oracle
        return oracle

    P0.load_or_build_oracle = wrapper


def _wrap_solve_chain():
    """Stash `region_arrays` — `solve_chain` receives it and calls `build_region_init` after.

    ⚠ It was ``CAL.region_sweep`` until the sweep was renamed, and this harness kept wrapping a name
    that no longer existed — an `AttributeError` on the first arm, so the instrument was DEAD while the
    suite stayed green. `TRAPS: a-green-suite-hid-five-dead-instruments`; `--self-test` is the gate that
    catches it now."""
    orig = CAL.solve_chain

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _CTX["region_arrays"] = region_arrays
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.solve_chain = wrapper


def _location_estimate(chain, statics, geometry, region_arrays, variant: str):
    """⭐⭐⭐ ψ's per-slot reference MEAN ``m_i``, PRIOR-FREE — the arm's whole payload.

    ⛔ **Nothing here reads a deconvolved array.** It reads the payload's own per-slot totals, the
    annotation, and the two opportunity models — which is what makes an estimate admissible at pass-0,
    where the gDNA landscape does not exist yet.

    ``m_i = rho_g,i * E_g,i / M_i``: the gDNA the object's own density predicts, as a share of what it
    actually holds. ⭐ **RNA is the RESIDUAL and is never predicted** — the asymmetry that makes this
    possible at all, since gDNA is near-uniform pre-solve while RNA spans six decades with no genomic
    autocorrelation (owner, 2026-08-16).

    ``variant`` selects where the density comes from, so the arms differ in ONE thing:

    ``struct``  only where mature RNA CANNOT be (``~mrna_active``) — the annotation-determined strata.
                psi's shipped 1/2 everywhere else, i.e. the landscape's population is left alone
    ``pooled``  the above, plus one pooled off-target density applied everywhere else
    ``local``   ⭐ the above, but an exon takes the density of its OWN flanking ``exon|intron`` splice
                boundaries when it has any. The boundary sits in the same probe footprint as its exon,
                so the capture enrichment CANCELS and never has to be estimated
    """
    import importlib.util as _il

    if "object_composition" not in sys.modules:
        _sp = _il.spec_from_file_location("object_composition", DESIGN / "object_composition.py")
        _m = _il.module_from_spec(_sp)
        sys.modules["object_composition"] = _m
        _sp.loader.exec_module(_m)
    OC = sys.modules["object_composition"]
    cls = OC.strata(chain, statics, geometry, region_arrays)
    label = cls["label"]
    # ⭐ ONE accessor for the pair, the same one the solver itself uses, so the estimate and the solve
    #   cannot end up on different bases: `region_contained` at a REGION, `boundary_unspliced` at a
    #   BOUNDARY, over the matching gDNA opportunity.
    M, eff_g = region_gdna_geometry(geometry)
    M = np.asarray(M, np.float64)
    eff_g = np.asarray(eff_g, np.float64)
    mature = np.asarray(statics.mrna_active_pos, bool) | np.asarray(statics.mrna_active_neg, bool)

    # the OFF-TARGET pooled density: intergenic + intron REGIONs, Sum mass / Sum E (never a mean of ratios)
    off = np.isin(label, list(OC.PURE_GDNA_STRATA))
    rho_off = OC.pooled_density(M, eff_g, off)
    rho = np.full(int(chain.n_slots), rho_off)

    if variant == "local":
        # ⭐ each slot's OWN flanking `exon|intron` splice boundaries, pooled over the one or two of them
        ei = (label == OC.ONTARGET_GDNA_STRATUM) & (
            np.asarray(geometry.eff_sj, np.float64).sum(1) > 0.0
        )
        num = np.zeros(int(chain.n_slots))
        den = np.zeros(int(chain.n_slots))
        for side in (np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)):
            ok = side >= 0
            q = np.zeros(int(chain.n_slots), bool)
            q[ok] = ei[side[ok]]
            b = side[q]
            num[q] += M[b]
            den[q] += eff_g[b]
        rho = np.where(den > 0.0, num / np.maximum(den, 1e-12), rho_off)

    # ⭐⭐⭐ THE STRUCTURAL CLAIM, WITH ITS OWN STRENGTH: where mature RNA CANNOT be, the expected gDNA
    #   count is ``rho_g * E_g,i`` and the RNA floor is exactly ONE PSEUDO-FRAGMENT **AT THIS OBJECT**::
    #
    #       m_i = E[g]_i / (E[g]_i + 1)
    #
    #   ⛔⛔ **The "one" must be one fragment HERE, not one over the genome's pooled opportunity.** A flat
    #   epsilon (``1/Sum eff_g`` = 1e-8) reads ``m = 1 - 1e-8``, at which **99.1 % of the reference's mass
    #   falls OUTSIDE the shipped L = 10 window** — worse than the (a,b) route this term exists to avoid,
    #   and the answer becomes a function of the grid. With the per-object floor the mass inside is
    #   **0.909-0.990** and ``m`` spans p5 0.788 / median 0.917 / p95 0.998 (measured, g50 ss0.99 off).
    #   ⭐ It is also the right MEANING: a big intergenic region may assert "pure gDNA" confidently, a
    #   tiny one may not, and the strength of the claim scales with the count the object expects.
    # ⛔⛔ **TWO FLOORS, AND THE PANEL PREFERS THE HARD ONE — kept as named variants because the choice
    #   is a MEASUREMENT, not a preference.**
    #   `struct`      pins `m -> 1 - eps` with a flat `eps = 1/Sum eff_g` (~1e-8).
    #   `struct_soft` uses the derived per-object floor `E[g]/(E[g]+1)` — ONE pseudo-fragment of RNA at
    #                 this object — which is the principled reading of the derivation.
    #   ⭐ Measured on the panel, final Sum|d| per stratum, ratio to base: hard **0.381 / 0.659 / 0.363 /
    #   0.800**, soft 0.609 / 1.045 / 0.580 / 1.000. The soft floor is WORSE on every stratum.
    #   ⭐⭐ The reason is the finding the whole investigation started from: on a structurally pure-gDNA
    #   object the truth IS `f_g = 1`, and objects at true `f_g ~ 1` carry **49-83 %** of in-scope error
    #   precisely because psi holds them off that vertex. A near-improper prior pushes them onto it; the
    #   soft floor pulls them back off (its median `m` is 0.917).
    #   ⚠ **The cost of the hard floor is real and this panel cannot see it**: at `m = 1 - 1e-8` only
    #   **0.94 %** of the reference's mass lies inside the shipped `L = 10` window, so the answer is a
    #   function of the grid (`TRAPS: a-clamp-at-the-closed-end-escapes-the-window`). ⛔ What is actually
    #   un-chosen here is the reference's STRENGTH on those slots — the second Beta degree of freedom,
    #   which every measurement in this project has held at `a + b = 1`. That is the open item, and
    #   picking either floor by its panel number alone would be tuning.
    expected_g = rho * eff_g
    if variant == "struct_soft":
        return np.where(mature, 0.5, expected_g / (expected_g + 1.0))
    if variant == "struct":
        e = 1.0 / max(float(np.sum(eff_g)), 1.0)
        return np.where(mature, 0.5, 1.0 - e)
    if variant == "struct_grid":
        # ⭐⭐⭐ **THE CAP IS THE GRID'S OWN RANGE, SO THERE IS NO CONSTANT TO CHOOSE** (owner,
        #   2026-08-16). A prior may not assert more than the lattice can represent, so the strongest
        #   admissible location is the highest grid point, ``m = sigma(L)``.
        #   ⭐ Its strength is then EXACTLY ``L`` nats — the term is worth ``log(1/eps)`` nats and
        #   ``eps = sigma(-L)`` — and the override budget is ``L / log(2*kappa)`` fragments: **14.6** at
        #   ``L = 10, kappa = 0.99``, against **27.0** for the flat clamp's 18.42 nats. ⛔ ``L`` is
        #   already ``CalibrationConfig.sweep_logodds_window``, gated and chosen; it is not a new number.
        from scipy.special import expit as _expit

        return np.where(mature, 0.5, float(_expit(_CTX["logodds_window"])))
    # the peel elsewhere: RNA is the RESIDUAL and is never predicted. ⚠ `pooled`/`local` carry the SOFT
    # floor on the structural slots, which is what their panel rows were measured with.
    m = np.clip(rho * eff_g / np.maximum(M, 1e-12), 0.0, 1.0)
    return np.where(mature, m, expected_g / (expected_g + 1.0))


def _count_structural_builder():
    """Count calls to the PRODUCTION location builder, without replacing it.

    ⛔ Every other arm here installs an override, so its fire counter is a property of the harness. This
    one measures the shipped path, so the counter must not become one: the real function is called and its
    result returned untouched, and only the tally is the harness's."""
    real = SL.structural_reference_location

    def counted(statics, logodds_window):
        _FIRED["loc"] += 1
        return real(statics, logodds_window)

    SL.structural_reference_location = counted
    SW.structural_reference_location = counted


def _install_reference_location(variant: str):
    """Inject ``m_i`` at ψ's ONE ``CompositionPriors`` construction site.

    ⭐ ``sweep`` holds the only construction; ``select``/``regrid`` rebuild through the real class inside
    ``simplex_logodds``, so patching the name in ``sweep`` intercepts exactly once and the slicing and
    regridding paths stay the shipped ones."""
    real = SL.CompositionPriors

    def factory(gdna=None, rna=None, location=None):
        _FIRED["loc"] += 1
        return real(gdna=gdna, rna=rna, location=_CTX.get("location"))

    SW.CompositionPriors = factory

    orig = CAL.solve_chain

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _CTX["region_arrays"] = region_arrays
        # ⛔ Read L from the CALL, never from a default: `struct_grid`'s whole claim is that its strength
        #   is the grid's own range, so silently assuming 10.0 while the solve ran on another window
        #   would make the arm's strength and its justification disagree with no tell.
        if "logodds_window" not in kw:
            raise SystemExit(
                "⛔ solve_chain was called without `logodds_window`; `ref_loc` cannot derive the grid cap."
            )
        _CTX["logodds_window"] = float(kw["logodds_window"])
        _CTX["location"] = (
            None
            if variant == "noop"
            else _location_estimate(chain, statics, geometry, region_arrays, variant)
        )
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.solve_chain = wrapper


def _truth_fg_per_slot(chain):
    """The ORACLE's true ``f_g`` per SLOT, and the mass behind it.

    ⚠ ``RegionInit``'s arrays are indexed by SLOT (the chain's alternating REGION/BOUNDARY sequence), while the
    oracle's masses are per OBJECT on two separate axes — so the mapping goes through
    ``chain.kind``/``chain.obj_idx`` rather than being assumed."""
    oracle = _CTX.get("oracle")
    ra = _CTX.get("region_arrays")
    if oracle is None or ra is None:
        return None, None
    ov = oracle.override_masses(ra)
    g = {
        "region": np.asarray(ov["mass_gdna_region"], np.float64),
        "boundary": np.asarray(ov["mass_gdna_boundary"], np.float64),
    }
    r = {
        "region": np.asarray(ov["mass_rna_region"], np.float64),
        "boundary": np.asarray(ov["mass_rna_boundary"], np.float64),
    }
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region = kind == REGION
    n = int(chain.n_slots)
    tg = np.zeros(n)
    tr = np.zeros(n)
    for axis, msk in (("region", is_region), ("boundary", ~is_region)):
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
    orig = NI.build_region_init

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
            #   `base` is the assertion; anything else means the wrapper itself moves the answer (TRAPS: byte-identity-gate).
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
        #   confidence, or the relay simply argues it back (TRAPS: substitution-understates-a-source).
        lock[tgt] = True

        v_fg, v_fr = NI.own_composition_logvar(f_g, tau, lock)
        M, E_g = region_gdna_geometry(geometry)
        M = np.asarray(M, np.float64)
        E_g = np.asarray(E_g, np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        n_region = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
        rho_g = np.maximum(
            np.where((M > _EPS) & (E_g > _EPS), f_g * M / np.maximum(E_g, _EPS), 0.0), 0.0
        )
        prec_g = NI.own_precision(n_region, v_fg, rho_g > _EPS)
        touched = np.zeros(f_g.shape[0], bool)
        touched[tgt] = True

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
            f_g=f_g, f_pos=f_pos, f_neg=f_neg,
            rho_g=rho_g, rho_pos=rho_pos, rho_neg=rho_neg,
            prec_g=prec_g, prec_pos=prec_pos, prec_neg=prec_neg,
            struct_lock=lock, tau_lam=tau,
        )

    SW.build_region_init = wrapper


def _install_ref_exponent(a_value: float, b_value: float | None = None):
    """⭐⭐ ψ's two reference exponents as free numbers instead of the single shipped ½.

    ⭐⭐⭐ **THEY ARE PSEUDO-COUNTS, AND THE PAIR IS A Beta.** ``a·log f_g + b·log(1−f_g)`` on the λ grid
    is exactly ``Beta(a, b)`` in ``f_g`` — the Jacobian ``|df/dλ| = f(1−f)`` turns ``f^{a−1}(1−f)^{b−1}``
    into ``f^a (1−f)^b``. So the pair has a STRENGTH ``a+b`` (one prior pseudo-fragment at the shipped
    value) and a MEAN ``a/(a+b)``, and the shipped ``a = b = ½`` fixes the mean at ½ — which asserts the
    library is half gDNA.

    ⛔ ``a = b = 0`` makes ψ improper on both sides (Beta(0,0) — TRAPS: no-prior-means-haldane's Haldane,
    a vertex amplifier), so small exponents are a PROTOTYPE that bounds what a derived rule could buy,
    never the rule itself.

    ⚠ ``b_value = None`` keeps the two equal, which is the historical one-knob behaviour; passing both
    is what an unequal-arm design needs. ⛔ The RNA arm takes ``rna_logprior`` exactly as its gDNA twin
    takes ``global_logprior``: this patch shadowed a one-argument ``_rna_arm`` after the arms were made
    symmetric, which raised on the first solve."""
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
    """Read two or more arm files and print the per-axis deltas, with the fixed-denominator columns
    first because those are the ones that cannot be gamed by knowing less (TRAPS: honesty-metrics-reward-ignorance)."""
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
    # ⛔ TRAPS: honesty-metrics-reward-ignorance ORDER: the two FIXED-DENOMINATOR columns come FIRST, because they are the only two that
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
        # ⛔ A byte-identical arm is NOT evidence of no change (TRAPS: hard-labels-miss-soft-change/TRAPS: byte-identity-gate) — EXCEPT for `noop`,
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
                      f" — if that arm was meant to CHANGE something, it did not fire (TRAPS: an-ablation-that-never-ran/TRAPS: byte-identity-gate).")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default=None,
                    help="base | noop | vertex_free | vertex_all | ref_c=<float> "
                         "| ref_loc={noop,struct,struct_grid,struct_soft,pooled,local}")
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
    _wrap_solve_chain()
    arm = args.arm
    expect_fire: list[str] = []
    structural_reference = False
    if arm == "config_struct":
        # ⭐⭐⭐ THE SHIPPED PATH, NOT AN OVERRIDE. Nothing is patched: the flag is set on the config and
        #   `calibrate` threads it to `solve_chain`'s ONE `CompositionPriors` site. ⚠ The counter wraps the
        #   PRODUCTION builder rather than replacing it, so TRAPS: an-ablation-that-never-ran is still
        #   satisfied without the measurement running through a wrapper (TRAPS: byte-identity-gate).
        structural_reference = True
        _count_structural_builder()
        expect_fire = ["loc"]
    elif arm == "vertex_free":
        _install_vertex_pin(True)
        expect_fire = ["pinned"]
    elif arm == "vertex_all":
        _install_vertex_pin(False)
        expect_fire = ["pinned"]
    elif arm == "noop":
        # ⭐ the harness's OWN falsification: the same wrapper, pinning nothing. Must be byte-identical
        #   to `base`, and if it is not, the wrapper itself is changing the answer (TRAPS: byte-identity-gate).
        _install_vertex_pin(True, force_empty=True)
        expect_fire = ["init"]
    elif arm.startswith("ref_loc="):
        # ⭐ the per-slot reference MEAN on the REAL panel, before any of it goes into `src/`
        #   (TRAPS: panel-before-src — four toy-positive changes have been panel-negative).
        _install_reference_location(arm.split("=", 1)[1])
        expect_fire = ["loc"]
    elif arm.startswith("ref_c="):
        # ⭐ `ref=C` keeps both arms equal; `ref=A,B` drives them unequal (the Beta(a,b) design).
        _spec = arm.split("=", 1)[1]
        _install_ref_exponent(*(float(x) for x in _spec.split(",")))
        expect_fire = ["ref_g", "ref_r"]
    elif arm != "base":
        ap.error(f"unknown arm {arm!r}")

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig(structural_reference=structural_reference)
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
            # ⚠ TRAPS: could-the-arm-have-fired: the opportunity count, printed beside the result.
            print(f"  {name} {time.perf_counter() - t0:.0f}s   pinned={fired.get('pinned', 0)}"
                  f"  init_calls={fired.get('init', 0)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
