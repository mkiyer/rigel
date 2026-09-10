#!/usr/bin/env python3
"""WHICH SLOTS TRAIN THE gDNA LANDSCAPE PRIOR, WITH WHAT EVIDENCE, AND HOW MUCH OF THAT TRAINING IS
FALSE? — the refit's training population censused per refit iteration, per node class and per
evidence class, against certified slot truth. No mechanism runs; nothing is proposed.

⭐⭐⭐ **THE QUESTION.** `calibrate._fit_gdna_hyperprior` fits `landscape.DensityLandscape` on the
previous sweep's solved belief (``f_g · mass`` per training region, weighted by
`landscape._reliability`) and the refit loop re-solves with it, ``calib_refit_iters`` times. The
prior therefore learns whatever the previous sweep believed, and a slot that had no evidence of
its own believes the reference and the prior — so a false positive at pass one can be taught back
to every blind slot at pass two (`ISSUES: gdna-landscape-trains-on-false-positives`). The owner's
ruling (2026-09-06) is that nodes whose only evidence is a BOUND do not train the prior; it LANDED
2026-09-10 as `RegionBelief.informed` (gate `tests/calibration/test_landscape_training_population.py`),
so the `delivered:bound` and `none` classes below are now OUTSIDE the training population and the
census shows them absent — run it at a commit before the landing to see what they trained at.
This instrument says WHO trains the prior, with WHAT evidence, at WHAT value, with WHAT weight — and
how much of that is false against the certified truth.

**What it records, per refit iteration.** It spies three things and changes nothing:

* every sweep's ``_capture`` (the belief that the NEXT fit trains on, ``tau_lam``, the intron
  factory's λ-factor) and the two held messages per slot at that sweep's solve
  (`_PreparedTransfer.solve`'s ``from_left`` / ``from_right``);
* every `fit_landscape` call's inputs — the exact ``(count, mass, eff, var, anchor)`` the estimator
  saw — and the landscape it returned.

The training selector is RE-DERIVED here from the same statics and GATED bit for bit against the
recorded inputs (a drifted re-derivation refuses to score), so the census is of the fit that ran.

**The evidence classes** — each slot's STRONGEST evidence at the sweep whose belief it trained
with, from the solver's own quantities and nothing invented:

    anchor                 the zero-count structural anchor (no unspliced mass; trains at 0, w = 1)
    locked                 structurally pure gDNA — neither RNA strand admissible
                           (`region_geometry.g1_locked`): trains at its whole mass, certain
    own:strand             the strand channel is live — `region_init.has_own_composition_evidence`
                           on ``tau_lam`` less the factory arm (`density_factor_precision`)
    own:factory            an ss intron's density-deconvolution factor is live (no strand)
    delivered:composition  no own channel; a COMPOSITION row arrived from a neighbour (a strand
                           profile through a face map, a splice-out row, a level-kept map)
    delivered:bound        no own channel, no composition; only a LEVEL arrived (the gDNA lane's
                           lower side, an RNA lane's ceiling, a cube row) — "a bound only"
    none                   no own channel, nothing delivered: the belief is the reference and the
                           prior alone — pure echo at every refit after the first

**How to read the tables.** Per refit iteration: ``Σw`` is what the estimator actually sums (each
region's Poisson kernel carries ``w``, not its mass), so ``Σw share`` is a class's say in the
landscape; ``gDNA trained`` is ``Σ f_g · M`` (the mass the class asserts is gDNA); ``FP w`` is the
share of the class's ``Σw`` on slots whose CERTIFIED gDNA is zero but which trained at ONE fragment or
more (below one the kernel centres at the resolution wall like the anchor's: the estimator's own
``max(count, 1)`` floor) — the false-positive training weight; ``centre`` is the median trained centre
``log10(max(count,1)/eff)`` and ``Δdec`` its median offset from the certified centre on slots with
true gDNA. ⛔ On a ``g00`` row EVERY slot is certified zero, so ``FP w`` there is simply the weight
share that trains at a fragment or more — the whole thing is false, and the zero controls are where
the entry's number comes from. ⚠ The training population is REGIONs only (boundaries are excluded by
the owner's ruling of 2026-07-27), so the node classes are the three region strata, exons split by
reach as `policy_benchmark.py --by-class` splits them.

**`--estimator` — THE METHODOLOGY AUDIT (owner, 2026-09-09).** For the last fit, the same
population re-fitted at the CERTIFIED values (``count = n_gdna``) with the weights as trained, and
again at ``var = 0``: the earth-mover distance in decades between the shipped landscape and each
tells apart the error the training VALUES put into the prior from what the estimator does with a
true population. It prices nothing downstream — `calibration_walk.py`'s C→E rung does that.

    python scripts/design/landscape_training_census.py --panel test --conditions gdna_g00_ss_0.50_nrna_file_capture_off
    python scripts/design/landscape_training_census.py --panel ladder --conditions gdna_g00_ss_0.50_nrna_mid_capture_off --estimator
    python scripts/design/landscape_training_census.py --panel ladder --g00
    python scripts/design/landscape_training_census.py --self-test
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


import importlib  # noqa: E402

# the package re-exports the FUNCTION `calibrate`, which shadows the submodule on attribute access
CAL = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration import landscape as LS  # noqa: E402
from rigel.calibration.density_deconv import density_factor_precision  # noqa: E402
from rigel.calibration.messages import Level, Message  # noqa: E402
from rigel.calibration.messages import transfer as TR  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_geometry import g1_locked  # noqa: E402
from rigel.calibration.region_init import has_own_composition_evidence  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

PB = _sibling("policy_benchmark.py")
PANELS = PB.PANELS
POLICIES = PB.POLICIES

EVIDENCE = (
    "anchor",
    "locked",
    "own:strand",
    "own:factory",
    "delivered:composition",
    "delivered:bound",
    "none",
)
_LN10 = np.log(10.0)


# ── the spy ──────────────────────────────────────────────────────────────────────────────────────────


class Spy:
    """Record every sweep's capture and held messages, and every landscape fit's inputs. Changes
    nothing: each wrapped function is called unchanged and its result returned as is."""

    def __init__(self):
        self.sweeps: list[dict] = []
        self.fits: list[dict] = []
        self._pending: dict = {}

    def __enter__(self):
        spy = self
        self._orig = (CAL.solve_chain, CAL.fit_landscape, TR._PreparedTransfer.solve)
        orig_solve_chain, orig_fit, orig_psolve = self._orig

        def solve_chain(*a, **k):
            out = orig_solve_chain(*a, **k)
            spy.sweeps.append(dict(capture=k.get("_capture"), belief=out, **spy._pending))
            spy._pending = {}
            return out

        def psolve(self_, from_left, from_right):
            msg = orig_psolve(self_, from_left, from_right)
            spy._pending = dict(from_left=list(from_left), from_right=list(from_right), msg=msg)
            return msg

        def fit(
            count, mass, eff, var, *, anchor, strength=1.0, knn_scale=LS._KNN_SCALE, domain=None, prev=None
        ):
            ls = orig_fit(
                count, mass, eff, var, anchor=anchor, strength=strength, knn_scale=knn_scale,
                domain=domain, prev=prev,
            )
            spy.fits.append(
                dict(
                    count=np.array(count, np.float64),
                    mass=np.array(mass, np.float64),
                    eff=np.array(eff, np.float64),
                    var=np.array(var, np.float64),
                    anchor=np.array(anchor, bool),
                    domain=domain,
                    prev=prev,
                    landscape=ls,
                )
            )
            return ls

        CAL.solve_chain, CAL.fit_landscape, TR._PreparedTransfer.solve = solve_chain, fit, psolve
        return self

    def __exit__(self, *exc):
        CAL.solve_chain, CAL.fit_landscape, TR._PreparedTransfer.solve = self._orig
        return False


# ── the pure pieces (self-tested) ───────────────────────────────────────────────────────────────────


def training_selector(
    kind, obj_idx, signature, free_pos, free_neg, mass_global, eff_global, informed=None
):
    """The refit's training population, re-derived from the statics exactly as
    `calibrate._fit_gdna_hyperprior` selects it: expressed REGIONs that are single-strand or
    structurally locked AND hold a composition (`RegionBelief.informed`, since 2026-09-10: a slot whose
    only evidence is a bound, or which has none, does not train), plus the zero-count anchor (an
    intergenic or intronic region with opportunity and no unspliced mass). Returns ``(sel, anchor)``;
    the caller GATES it against the recorded fit."""
    isr = np.asarray(kind) == REGION
    fp = np.asarray(free_pos, bool)
    fn = np.asarray(free_neg, bool)
    rtype = coarse_type_array(np.asarray(signature))
    ridx = np.clip(np.asarray(obj_idx, np.int64), 0, rtype.shape[0] - 1)
    mass = np.asarray(mass_global, np.float64)
    eff = np.asarray(eff_global, np.float64)
    expressed = isr & (eff > 1.0e-9) & (mass > 1.0e-12)
    anchor = isr & (eff > 1.0e-9) & (mass <= 1.0e-12) & (rtype[ridx] != RegionType.EXON)
    sel = expressed & ((fp ^ fn) | (~fp & ~fn))
    if informed is not None:
        sel &= np.asarray(informed, bool)
    sel |= anchor
    return sel, anchor


def gate_selector(sel, anchor, f_g, mass_global, eff_global, fit: dict) -> None:
    """Refuse unless the re-derived population reproduces the recorded fit inputs BIT FOR BIT."""
    count = np.asarray(f_g, np.float64)[sel] * np.asarray(mass_global, np.float64)[sel]
    checks = (
        ("count", count, fit["count"]),
        ("mass", np.asarray(mass_global, np.float64)[sel], fit["mass"]),
        ("eff", np.asarray(eff_global, np.float64)[sel], fit["eff"]),
        ("anchor", anchor[sel], fit["anchor"]),
    )
    for name, mine, theirs in checks:
        if mine.shape != theirs.shape or not np.array_equal(mine, theirs):
            raise AssertionError(
                f"the re-derived training selector does not reproduce the recorded fit's `{name}` "
                f"({mine.shape} vs {theirs.shape}) — the selector in `_fit_gdna_hyperprior` has "
                "drifted from this census's copy; fix the copy, do not score"
            )


def training_weights(fit: dict) -> np.ndarray:
    """Each recorded training slot's `_reliability` weight, exactly as the estimator applied it
    (its own ``live`` filter first; a slot it dropped weighs 0)."""
    count, eff, mass = fit["count"], fit["eff"], fit["mass"]
    live = np.isfinite(count) & np.isfinite(eff) & np.isfinite(mass) & (eff > LS._EPS)
    w = np.zeros(count.shape[0])
    w[live] = LS._reliability(np.maximum(count[live], 0.0), fit["var"][live], fit["anchor"][live])
    return w


def held_evidence(from_left, from_right, cube_slots, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Per slot: did a COMPOSITION row arrive on either side, and did a LEVEL (any lane) or a cube
    row — a bound — arrive on either side."""
    comp = np.zeros(n, bool)
    bound = np.zeros(n, bool)
    if from_left is None or from_right is None:
        return comp, bound
    for i in range(n):
        for m in (from_left[i], from_right[i]):
            if m is None:
                continue
            if m.composition is not None:
                comp[i] = True
            if any(getattr(m, lane) is not None for lane in Message.LANES[1:]):
                bound[i] = True
    for i in cube_slots:
        bound[int(i)] = True
    return comp, bound


def classify_evidence(tau_lam, tau_factory, comp, bound, anchor, locked=None) -> np.ndarray:
    """Each slot's STRONGEST evidence, by the order in :data:`EVIDENCE`. ``tau_lam`` is the solver's
    combined λ precision, ``tau_factory`` the factory arm read off the captured λ-factor; the strand
    arm is their difference; ``locked`` is `g1_locked` (structural certainty, above every channel)."""
    tau = np.asarray(tau_lam, np.float64)
    fac = np.zeros_like(tau) if tau_factory is None else np.asarray(tau_factory, np.float64)
    own = has_own_composition_evidence(tau)
    strand = has_own_composition_evidence(tau - fac)
    cls = np.full(tau.shape[0], "none", dtype=object)
    cls[np.asarray(bound, bool)] = "delivered:bound"
    cls[np.asarray(comp, bool)] = "delivered:composition"
    cls[own & ~strand] = "own:factory"
    cls[strand] = "own:strand"
    if locked is not None:
        cls[np.asarray(locked, bool)] = "locked"
    cls[np.asarray(anchor, bool)] = "anchor"
    return cls.astype(str)


def one_sided(row, eps: float = TR.EPS) -> bool:
    """Is a max-normalised log-profile a BOUND — its maximum reached at a grid end (a plateau that
    runs off the grid) — rather than a two-sided profile with an interior mode?"""
    r = np.asarray(row, np.float64)
    top = float(r.max())
    return bool(r[0] >= top - eps or r[-1] >= top - eps)


def census_rows(
    sel_idx, w, count, eff, true_gdna, node_class, evidence
) -> dict[tuple[str, str], dict]:
    """The census cells: ``(node class, evidence class) -> {n, sum_w, gdna, fp_w, centre, ddec}``.
    ``fp_w`` is the weight on certified-zero slots that trained at a positive count."""
    centre = np.log10(np.maximum(count, 1.0)) - np.log10(np.maximum(eff, LS._EPS))
    truth = np.asarray(true_gdna, np.float64)[sel_idx]
    t_centre = np.log10(np.maximum(truth, 1.0)) - np.log10(np.maximum(eff, LS._EPS))
    # a trained count below ONE fragment centres at the resolution wall exactly as the anchor does
    # (the estimator's own `max(count, 1)` floor), so it is not a false location; ≥ 1 is
    fp = (truth <= 0.0) & (count >= 1.0)
    out: dict[tuple[str, str], dict] = {}
    nc = np.asarray(node_class)[sel_idx]
    ev = np.asarray(evidence)[sel_idx]
    for key in sorted(set(zip(nc.tolist(), ev.tolist()))):
        m = (nc == key[0]) & (ev == key[1])
        pos = m & (truth > 0.0)
        out[key] = dict(
            n=int(m.sum()),
            sum_w=float(w[m].sum()),
            gdna=float(count[m].sum()),
            gdna_zero=float(count[m & (truth <= 0.0)].sum()),
            fp_w=float(w[m & fp].sum()),
            mean_w=float(w[m].mean()) if m.any() else float("nan"),
            centre=float(np.median(centre[m & (count > 0.0)])) if (m & (count > 0.0)).any() else float("nan"),
            ddec=float(np.median(centre[pos] - t_centre[pos])) if pos.any() else float("nan"),
        )
    return out


def emd_decades(a: LS.DensityLandscape, b: LS.DensityLandscape) -> float:
    """Earth-mover distance between two landscapes on ONE grid, in decades."""
    if a.log_rho.shape != b.log_rho.shape or not np.allclose(a.log_rho, b.log_rho):
        raise ValueError("the two landscapes are on different grids; EMD needs one axis")
    pa, pb = np.exp(a.logP), np.exp(b.logP)
    pa, pb = pa / pa.sum(), pb / pb.sum()
    step = float(a.log_rho[1] - a.log_rho[0]) / _LN10
    return float(np.abs(np.cumsum(pa) - np.cumsum(pb)).sum() * step)


# ── one condition ────────────────────────────────────────────────────────────────────────────────────


def run_condition(index, region_arrays, sj, boundary_flags, cache_dir: Path, policy: str) -> dict:
    """Calibrate one cached condition under ``policy`` with the spy on; return the per-fit census."""
    cache = read_scan_cache(cache_dir / "_main", index)
    kw = calibration_inputs(cache, index)
    payload = kw["payload"]
    slots = dict(np.load(cache_dir / "slot_truth.npz", allow_pickle=True))
    node_class = PB._slot_classes(slots, payload, boundary_flags)
    true_gdna = np.asarray(slots["n_gdna"], np.float64)
    debug: dict = {}
    with Spy() as spy:
        result = CAL.calibrate(
            payload=payload,
            config=dataclasses.replace(CalibrationConfig(), **POLICIES[policy]),
            region_arrays=region_arrays,
            strand_model=kw["strand_model"],
            gdna_fl_pmf=kw["gdna_fl_pmf"],
            rna_fl_pmf=kw["rna_fl_pmf"],
            sj=sj,
            boundary_flags=boundary_flags,
            _debug=debug,
        )
    chain = debug["chain"]
    n = int(chain.n_slots)
    if len(spy.sweeps) != len(spy.fits) + 1:
        raise AssertionError(
            f"{len(spy.sweeps)} sweeps for {len(spy.fits)} fits — expected one more sweep than fits "
            "(sweep, fit, re-sweep, …); the spy missed a call"
        )
    fits = []
    for k, fit in enumerate(spy.fits):
        sw = spy.sweeps[k]  # the sweep whose belief this fit trained on
        cap = sw["capture"]
        sel, anchor = training_selector(
            chain.kind,
            chain.obj_idx,
            region_arrays.signature,
            cap["free_pos"],
            cap["free_neg"],
            cap["mass_global"],
            cap["eff_global"],
            informed=sw["belief"].informed,
        )
        gate_selector(sel, anchor, cap["f_g"], cap["mass_global"], cap["eff_global"], fit)
        fg_grid = np.asarray(cap["solve_grid"], np.float64)  # the capture's grid is f_g = σ(λ)
        lam_grid = np.log(fg_grid) - np.log1p(-fg_grid)
        fac = density_factor_precision(cap.get("intron_prior"), lam_grid)
        msg = sw.get("msg")
        cube = () if msg is None or msg.cube_rows is None else tuple(msg.cube_rows)
        comp, bound = held_evidence(sw.get("from_left"), sw.get("from_right"), cube, n)
        locked = g1_locked(cap["free_pos"], cap["free_neg"])
        evidence = classify_evidence(cap["_tau0_lam"], fac, comp, bound, anchor, locked)
        rows_held = cap.get("lam_rows")
        onesided = None
        if rows_held is not None:
            onesided = np.array([one_sided(rows_held[i]) if np.ptp(rows_held[i]) > TR.EPS else False for i in range(n)])
        w = training_weights(fit)
        sel_idx = np.flatnonzero(sel)
        cells = census_rows(sel_idx, w, fit["count"], fit["eff"], true_gdna, node_class, evidence)
        fits.append(
            dict(
                iteration=k + 1,
                cells=cells,
                n_train=int(sel.sum()),
                sum_w=float(w.sum()),
                gdna=float(fit["count"].sum()),
                gdna_zero=float(fit["count"][true_gdna[sel_idx] <= 0.0].sum()),
                fp_w=float(w[(true_gdna[sel_idx] <= 0.0) & (fit["count"] >= 1.0)].sum()),
                onesided_share=_onesided_share(onesided, sel_idx, evidence),
                sel=sel,
                evidence=evidence,
                fit=fit,
                w=w,
            )
        )
    return dict(
        fits=fits,
        true_gdna=true_gdna,
        node_class=node_class,
        n_slots=n,
        policy=debug["capture"].get("policy_name"),
        landscape=debug.get("gdna_hyperprior"),
        controls=zero_controls(result, slots),
        sweeps=spy.sweeps,  # the raw per-sweep captures and held messages, for a dissection
    )


def zero_controls(result, slots: dict) -> dict:
    """THE TWO ZERO CONTROLS INSIDE ONE CONDITION, from the shipped final answer (the same estimate
    `policy_benchmark.py` scores, so ``total`` reproduces its row): invented gDNA on the slots whose
    certified gDNA is zero, and invented RNA on the slots that are certified PURE gDNA — the silent
    genes and the RNA-free introns and intergenic regions — each in fragments over its own
    population. The truth is a constant on both, so every fragment counted is a false positive."""
    kind = np.asarray(slots["kind"])
    obj = np.asarray(slots["obj"], np.int64)
    n_gdna = np.asarray(slots["n_gdna"], np.float64)
    n_rna = np.asarray(slots["n_nrna"], np.float64) + np.asarray(slots["n_mrna"], np.float64)
    mass = np.asarray(slots["count"], np.float64)
    is_r = kind == REGION
    est = np.zeros(kind.shape[0])
    est[is_r] = np.asarray(result.mass_gdna_region, np.float64)[obj[is_r]]
    est[~is_r] = np.asarray(result.mass_gdna_boundary, np.float64)[obj[~is_r]]
    err = np.abs(est - n_gdna)
    zg = (n_gdna <= 0.0) & (mass > 0.0)
    zr = (n_rna <= 0.0) & (n_gdna > 0.0)
    return dict(
        total=float(err.sum()),
        zero_gdna=(float(err[zg].sum()), int(zg.sum()), float(mass[zg].sum())),
        zero_rna=(float(err[zr].sum()), int(zr.sum()), float(mass[zr].sum())),
    )


def _onesided_share(onesided, sel_idx, evidence) -> float:
    """Among the training slots whose strongest evidence is a delivered composition row, the share
    whose held row is nonetheless one-sided (a plateau to a grid end) — `ISSUES: two-sided-exon-row`."""
    if onesided is None:
        return float("nan")
    m = np.asarray(evidence)[sel_idx] == "delivered:composition"
    return float(onesided[sel_idx][m].mean()) if m.any() else float("nan")


def estimator_audit(res: dict) -> dict:
    """The last fit re-run at the certified values on the same population: EMD in decades between
    the shipped landscape and (a) the truth values under the weights as trained, (b) the truth values
    at ``var = 0`` (every slot trusted), (c) the truth values with the slots whose only evidence is
    a bound or nothing removed — the owner's population — at the weights as trained."""
    last = res["fits"][-1]
    fit, sel = last["fit"], last["sel"]
    sel_idx = np.flatnonzero(sel)
    truth = res["true_gdna"][sel_idx]
    shipped = fit["landscape"]
    out = {}
    arms = {
        "truth values, weights as trained": (truth, fit["var"], np.ones(sel_idx.size, bool)),
        "truth values, var = 0": (truth, np.zeros_like(fit["var"]), np.ones(sel_idx.size, bool)),
    }
    keep = ~np.isin(last["evidence"][sel_idx], ("delivered:bound", "none"))
    arms["truth values, bound-only and evidence-free slots removed"] = (truth, fit["var"], keep)
    arms["SHIPPED values, bound-only and evidence-free slots removed"] = (fit["count"], fit["var"], keep)
    for name, (count, var, m) in arms.items():
        ls = LS.fit_landscape(
            count[m], fit["mass"][m], fit["eff"][m], var[m], anchor=fit["anchor"][m],
            strength=shipped.strength, domain=fit.get("domain"), prev=fit.get("prev"),
        )
        if ls is None:
            out[name] = (float("nan"), 0)
            continue
        # the grid is derived from (mass, eff): a removed slot can move it, so interpolate onto the
        # shipped grid before comparing
        if ls.log_rho.shape != shipped.log_rho.shape or not np.allclose(ls.log_rho, shipped.log_rho):
            p = np.exp(np.interp(shipped.log_rho, ls.log_rho, ls.logP, left=ls.logP[0], right=ls.logP[-1]))
            ls = LS.DensityLandscape(shipped.log_rho, np.log(p / p.sum()), ls.n_train, ls.strength)
        out[name] = (emd_decades(shipped, ls), int(m.sum()))
    return out


# ── printing ─────────────────────────────────────────────────────────────────────────────────────────


def print_condition(name: str, res: dict, *, by_node: bool) -> None:
    print(f"\n== {name}   policy {res['policy']}   {res['n_slots']:,} slots")
    c = res["controls"]
    print(
        f"   final answer: whole-library |err| {c['total']:,.0f} fragments; ZERO-gDNA control (certified-zero "
        f"slots) {c['zero_gdna'][0]:,.0f} invented gDNA on {c['zero_gdna'][1]:,} slots / {c['zero_gdna'][2]:,.0f} "
        f"fragments; ZERO-RNA control (certified pure-gDNA slots) {c['zero_rna'][0]:,.0f} invented RNA on "
        f"{c['zero_rna'][1]:,} slots / {c['zero_rna'][2]:,.0f} fragments"
    )
    for f in res["fits"]:
        sw = f["sum_w"]
        print(
            f"\n   refit {f['iteration']}: {f['n_train']:,} training regions, Σw {sw:,.1f}, "
            f"gDNA trained {f['gdna']:,.0f} fragments of which {f['gdna_zero']:,.0f} on certified-ZERO "
            f"slots ({(f['gdna_zero'] / f['gdna'] if f['gdna'] else 0):.1%}); false-positive weight "
            f"{f['fp_w'] / sw if sw else 0:.1%} of Σw; one-sided rows among delivered:composition "
            f"{f['onesided_share']:.0%}"
        )
        cells = f["cells"]
        # by evidence class, node classes summed (one condition, one fit — not a pooling across rows)
        print(f"   {'evidence':<24}{'n':>8}{'Σw':>10}{'w share':>9}{'FP w':>9}{'mean w':>8}{'gDNA':>12}{'on zero':>10}{'centre':>8}{'Δdec':>7}")
        for ev in EVIDENCE:
            ks = [k for k in cells if k[1] == ev]
            if not ks:
                continue
            n = sum(cells[k]["n"] for k in ks)
            s = sum(cells[k]["sum_w"] for k in ks)
            fp = sum(cells[k]["fp_w"] for k in ks)
            g = sum(cells[k]["gdna"] for k in ks)
            gz = sum(cells[k]["gdna_zero"] for k in ks)
            cs = [cells[k]["centre"] for k in ks if np.isfinite(cells[k]["centre"])]
            ds = [cells[k]["ddec"] for k in ks if np.isfinite(cells[k]["ddec"])]
            print(
                f"   {ev:<24}{n:>8,}{s:>10,.1f}{s / sw if sw else 0:>9.1%}{fp / s if s else 0:>9.1%}"
                f"{s / n if n else 0:>8.2f}{g:>12,.0f}{gz:>10,.0f}"
                f"{(np.median(cs) if cs else float('nan')):>8.2f}{(np.median(ds) if ds else float('nan')):>7.2f}"
            )
        if by_node:
            print(f"   {'node class × evidence':<62}{'n':>8}{'w share':>9}{'FP w':>9}{'gDNA':>12}{'on zero':>10}{'centre':>8}{'Δdec':>7}")
            for k in sorted(cells, key=lambda k: -cells[k]["sum_w"]):
                c = cells[k]
                print(
                    f"   {k[0] + '  ·  ' + k[1]:<62}{c['n']:>8,}{c['sum_w'] / sw if sw else 0:>9.1%}"
                    f"{c['fp_w'] / c['sum_w'] if c['sum_w'] else 0:>9.1%}{c['gdna']:>12,.0f}{c['gdna_zero']:>10,.0f}"
                    f"{c['centre']:>8.2f}{c['ddec']:>7.2f}"
                )


# ── self-test ────────────────────────────────────────────────────────────────────────────────────────


def _self_test() -> int:
    passed = 0

    def check(name, ok):
        nonlocal passed
        print(f"   {'✔' if ok else '✘'} {name}")
        passed += int(bool(ok))
        return ok

    K = 9
    lam = np.linspace(-4, 4, K)
    peaked = -0.5 * lam**2
    plateau = np.minimum(0.0, -0.5 * np.maximum(lam, 0.0) ** 2)  # flat below, falling above: a bound
    check("one_sided: an interior mode is two-sided", not one_sided(peaked))
    check("one_sided: a plateau to the grid end is a bound", one_sided(plateau))
    check("one_sided PERTURBED: the plateau tilted by one grid point's fall becomes two-sided",
          not one_sided(plateau - 0.01 * np.abs(lam)))

    # evidence classification on synthetic messages
    lv = Level(profile=plateau, n=3.0, a=100.0)
    fl = [None, Message(composition=peaked), Message(level_gdna=lv), None, Message(level_rna_pos=lv)]
    fr = [Message(), None, None, None, None]
    comp, bound = held_evidence(fl, fr, cube_slots=(3,), n=5)
    check("held_evidence: composition on slot 1 only", comp.tolist() == [False, True, False, False, False])
    check("held_evidence: bounds on slots 2 (gDNA level), 3 (cube), 4 (RNA level)",
          bound.tolist() == [False, False, True, True, True])
    tau = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    cls = classify_evidence(tau, None, comp, bound, np.zeros(5, bool)).tolist()
    check("classify: none / composition / bound / bound / bound",
          cls == ["none", "delivered:composition", "delivered:bound", "delivered:bound", "delivered:bound"])
    tau2 = np.array([1.0, 0.0, 5e-10, 0.0, 2.0])
    fac2 = np.array([0.0, 0.0, 0.0, 0.0, 2.0])
    cls2 = classify_evidence(tau2, fac2, comp, bound, np.array([False, False, False, True, False])).tolist()
    check("classify PERTURBED: strand beats none; τ under the guard is not a channel; a factory-only τ is "
          "own:factory over a bound; anchor beats everything",
          cls2 == ["own:strand", "delivered:composition", "delivered:bound", "anchor", "own:factory"])
    cls3 = classify_evidence(tau2, fac2, comp, bound, np.zeros(5, bool),
                             locked=np.array([True, False, False, False, False])).tolist()
    check("classify PERTURBED: a structurally locked slot outranks its own strand channel",
          cls3 == ["locked", "delivered:composition", "delivered:bound", "delivered:bound", "own:factory"])

    # the selector and its gate
    kind = np.array([REGION, 1, REGION, 1, REGION, 1, REGION])
    obj = np.array([0, 0, 1, 1, 2, 2, 3])
    sig = np.zeros(4, np.int64)  # every region's coarse type resolves through the same map
    fp = np.array([0, 0, 1, 1, 1, 1, 1], bool)
    fn = np.array([0, 0, 0, 0, 1, 1, 0], bool)
    mass = np.array([0.0, 5.0, 10.0, 5.0, 8.0, 5.0, 4.0])
    eff = np.full(7, 100.0)
    sel, anchor = training_selector(kind, obj, sig, fp, fn, mass, eff)
    rtype = coarse_type_array(sig)
    expect_anchor = rtype[0] != RegionType.EXON
    check("selector: the locked expressed region and the single-strand regions train; the AMBIG region "
          "does not; the empty region is the anchor iff it is not an exon",
          sel.tolist() == [bool(expect_anchor), False, True, False, False, False, True]
          and anchor.tolist() == [bool(expect_anchor), False, False, False, False, False, False])
    f_g = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    fit = dict(count=f_g[sel] * mass[sel], mass=mass[sel], eff=eff[sel], anchor=anchor[sel],
               var=np.zeros(int(sel.sum())))
    gate_selector(sel, anchor, f_g, mass, eff, fit)
    check("gate: the re-derived selector reproduces the recorded fit", True)
    bad = dict(fit, count=fit["count"] + np.eye(1, int(sel.sum()), 1)[0] * 1e-9)
    try:
        gate_selector(sel, anchor, f_g, mass, eff, bad)
        check("gate PERTURBED: a 1e-9 nudge in one recorded count is refused", False)
    except AssertionError:
        check("gate PERTURBED: a 1e-9 nudge in one recorded count is refused", True)

    # weights: the anchor weighs 1, a give-up slot collapses, a confident one keeps mass
    fit_w = dict(count=np.array([0.0, 50.0, 50.0]), mass=np.array([0.0, 100.0, 100.0]),
                 eff=np.array([100.0, 100.0, 100.0]), var=np.array([np.inf, 1e-4, 1e4]),
                 anchor=np.array([True, False, False]))
    w = training_weights(fit_w)
    check("weights: anchor 1.0, confident ≈ 1, give-up ≈ 0", w[0] == 1.0 and w[1] > 0.99 and w[2] < 0.01)

    # the census cells and the false-positive share
    sel_idx = np.array([0, 2, 4, 6])
    count = np.array([0.0, 3.0, 0.5, 4.0])  # 0.5 of a fragment centres at the wall: not a false location
    effc = np.full(4, 10.0)
    truth = np.array([0, 0, 0, 0, 0, 0, 4.0])
    node = np.array(["R intron"] * 7)
    ev = np.array(["anchor", "none", "none", "own:strand", "none", "none", "own:strand"])
    ww = np.array([1.0, 0.5, 0.5, 0.9])
    cells = census_rows(sel_idx, ww, count, effc, truth, node, ev)
    c_none = cells[("R intron", "none")]
    c_own = cells[("R intron", "own:strand")]
    check("census: the evidence-free cell's false-positive weight is the slot that trained at 3 on a zero",
          c_none["n"] == 2 and c_none["fp_w"] == 0.5 and c_none["gdna_zero"] == 3.5)
    check("census: the own cell is true (Δdec 0, no FP)", c_own["fp_w"] == 0.0 and c_own["ddec"] == 0.0)
    truth2 = truth.copy()
    truth2[2] = 3.0
    c2 = census_rows(sel_idx, ww, count, effc, truth2, node, ev)[("R intron", "none")]
    check("census PERTURBED: certifying that slot's 3 as true removes the FP weight", c2["fp_w"] == 0.0)

    # EMD on one grid
    g = np.linspace(-2, 2, 41) * _LN10
    p = np.exp(-0.5 * (g / _LN10 / 0.25) ** 2)  # narrow, so a rolled copy carries no wrapped tail
    a = LS.DensityLandscape(g, np.log(p / p.sum()), 1)
    q = np.roll(p, 5)
    b = LS.DensityLandscape(g, np.log(q / q.sum()), 1)
    check("emd: identical landscapes 0; a half-decade shift ≈ 0.5 decades",
          emd_decades(a, a) == 0.0 and abs(emd_decades(a, b) - 0.5) < 0.05)
    try:
        emd_decades(a, LS.DensityLandscape(g[:-1], a.logP[:-1], 1))
        check("emd PERTURBED: two grids refused", False)
    except ValueError:
        check("emd PERTURBED: two grids refused", True)

    print(f"\n{passed} passed")
    return 0 if passed == 17 else 1


# ── main ─────────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", choices=sorted(PANELS), default="test")
    ap.add_argument("--conditions", nargs="+", default=None, help="default: every certified condition")
    ap.add_argument("--g00", action="store_true", help="only the zero-gDNA rows (the zero controls)")
    ap.add_argument("--policy", choices=sorted(POLICIES), default="transfer")
    ap.add_argument("--by-node", action="store_true", help="also the node class × evidence table")
    ap.add_argument("--estimator", action="store_true", help="the methodology audit on the last fit")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return _self_test()

    index_dir, panel_dir = PANELS[args.panel]
    oracle = panel_dir / "oracle_cache"
    conditions = args.conditions or sorted(
        d.name for d in oracle.iterdir() if (d / "slot_truth.npz").exists()
    )
    if args.g00:
        conditions = [c for c in conditions if "_g00_" in c]
    if not conditions:
        raise SystemExit(f"no certified conditions under {oracle}")
    index = TranscriptIndex.load(str(index_dir))
    region_arrays = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    sj = build_sj_geometry_arrays(index)
    boundary_flags = build_boundary_flags_array(index)
    print(f"⭐ {args.panel} panel — the gDNA landscape prior's training population, per refit, against certified truth")
    print("   Σw is what the estimator sums; FP w is the share of it on certified-zero slots trained at ≥ 1 fragment.")
    for c in conditions:
        res = run_condition(index, region_arrays, sj, boundary_flags, oracle / c, args.policy)
        print_condition(c, res, by_node=args.by_node)
        if args.estimator and res["fits"]:
            print("\n   estimator audit — EMD in decades from the SHIPPED last fit (same population, same grid):")
            for name, (d, n) in estimator_audit(res).items():
                print(f"      {name:<62} {d:8.3f}   ({n:,} regions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
