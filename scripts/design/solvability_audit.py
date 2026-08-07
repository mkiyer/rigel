#!/usr/bin/env python
"""⭐⭐ **WHICH OBJECTS CAN BE SOLVED, AND OF THOSE, WHICH ARE SOLVED WRONG — AND CONFIDENTLY SO.**

⛔ **THE MISTAKE THIS REPLACES.** ``pass0_vs_oracle.py`` scores every object that carries mass and
mass-weights the result. In pass-0 that is the wrong question. Pass-0 is the **prior-free** solve, and
its job is not to be accurate everywhere — it is to produce a substrate the gDNA hyperprior can be
fitted against. A node with no own evidence that reports ``f_g ≈ ½`` at **zero precision** is not
making an error; it is correctly saying *I cannot be solved without a prior*, which is exactly true.
Counting that as error buries the thing that actually matters underneath it.

⭐⭐ **WHAT MATTERS IS THE OBJECT THAT IS SOLVABLE, IS SOLVED WRONG, AND IS SURE OF ITSELF.** That is
what corrupts the hyperprior fit, and through the prior it corrupts everything downstream. The whole
point of this instrument is to separate three populations that a mass-weighted error lumps together:

1. **UNDETERMINED** — no own evidence. ⛔ EXCLUDED from the error denominator. Its only failure mode
   is the opposite one: claiming a precision it has not earned.
2. **SOLVABLE and right.**
3. **SOLVABLE and wrong** — split by whether the solver was *confident*. Uncertainly wrong is a cost;
   **confidently wrong is a defect**, because a wrong value with a tight variance outvotes correct
   neighbours and anchors the prior.

⛔ **THE CONFIDENCE COMPARISON MUST BE IN LOG SPACE, AND THE FIELD NAME LIES ABOUT THIS.**
``var_gdna`` / ``gdna_frac_var`` is ``Var(log f_g)`` (`simplex_logodds`, required by TRAPS: two-gaussians-one-latent), while both
names read as the variance of the fraction. Writing ``|f_g − truth| / sqrt(var_gdna)`` is a linear
error over a log-space variance — the exact defect that, when corrected once before, moved a suite
total 0.046 → 1.007 and INVERTED the per-class ranking. So the standardised discrepancy here is::

    z = (log f_pred − log f_truth) / sqrt(var_gdna)

⚠ **and the truth is clipped to the solver's OWN λ-grid endpoints**, not to a chosen epsilon. Truth is
exactly 0 at a pure-RNA object and exactly 1 at a structurally-pure-gDNA one — both common — and
``log 0`` is not a number. The grid is the honest bound: the solver cannot express a fraction outside
it, so that is the best answer it could ever have given.

⭐ **NO TUNED THRESHOLD DECIDES ANYTHING HERE.** ``|z|`` bands are standard-deviation multiples, which
is a scale the quantity brings with it. And the headline is a **calibration curve** — realised RMS log
error against claimed sd, per precision decile — which needs no cut at all: a ratio of 1 means the
declared precision is earned, above 1 means overconfident. That curve is the answer to "is the solver
sure and wrong?" without anyone choosing a number.

⭐ **THE ABLATION LADDER IS FREE.** ``_debug["capture"]`` already carries the same solve at three
depths — ``fg_strand`` (strand likelihood alone), ``fg_loc`` (the message-free local self-solve: strand
+ intron factory + reference) and ``f_g`` (final, after the relay). Nothing is re-solved, so the rungs
are the solver's own arithmetic rather than a reimplementation of it. Comparing them says WHICH channel
moved an object off truth, and in particular recovers the class the old tooling called
``P1_OVERRULED``: strand had it right and confident, and the full solve overrode it.

⚠ On an unstranded library the strand rung is silent by construction (``I(f_g) ∝ (2κ−1)²`` is exactly
0 at κ = ½), so the FACTORY is the only own-evidence channel and is reported as its own rung. A ladder
without it would be blank on precisely the condition that matters most.

Usage::

    python scripts/design/solvability_audit.py --condition NAME [--suite DIR] [--top 20]
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


P0 = _sibling("pass0_vs_oracle.py")

from rigel.calibration.density_deconv import density_factor_precision  # noqa: E402
from rigel.calibration.node_chain import EDGE, NODE  # noqa: E402
from rigel.calibration.node_geometry import g1_locked  # noqa: E402
from rigel.calibration.simplex_logodds import _logodds_grid  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1.0e-9

#: ⚠ Standard-deviation multiples, not tuned constants: ``z`` carries its own scale, and these are
#: where a Gaussian's mass sits. Nothing branches on them — they are read points on a distribution.
Z_BANDS = ((0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, np.inf))

#: The own-evidence channels a slot can have in pass-0. ⚠ ``length`` is absent deliberately: the
#: per-node length likelihood is gated OFF, so it cannot contribute to what the SOLVER could get
#: right. C_info's identification test is a statement about the payload, not about the solver.
#:
#: ⛔ **THESE OVERLAP AND ARE NOT A PARTITION**, which a first draft of this instrument got wrong.
#: ``tau_lam = i_strand + tau_fac``, so a single-stranded INTRON node has BOTH — and that is the
#: best-evidenced kind of object there is, not a bookkeeping error. Only ``none`` is exclusive of the
#: rest. Mass summed over these rows therefore DOUBLE-COUNTS, and the report says so; the partition
#: that does add up is determined / undetermined.
CHANNELS = ("locked", "strand", "factory", "none")


#: ⭐⭐ **OWN-EVIDENCE STRENGTH IS A CONTINUUM AND MUST BE REPORTED AS ONE.** ``tau_lam`` is a Fisher
#: precision on ``λ = log(f_g/f_R)``, so a slot's own statement carries sd ``1/√τ``, in nats, against a
#: solver that can only represent ``λ ∈ [−L, +L]``. These are DECADE boundaries on that sd — read
#: points on a distribution, not a cut: nothing branches on them.
#:
#: ⛔ **A BINARY solvable/undetermined CUT AT ``τ > 1e-9`` WAS THE DEFECT, AND A BETTER THRESHOLD IS
#: NOT THE FIX.** The strand arm carries ``I(f_g) ∝ (2κ−1)²`` (`EQUATIONS.md` §5.2), *exactly* zero at
#: κ = ½; but κ is FITTED, so on a genuinely unstranded library it misses ½ by a few 1e-4 and τ lands
#: at ~1e-7 rather than at 0. ``τ > 1e-9`` then calls that a live channel and the object is scored as
#: SOLVABLE while its own statement has sd of 10³ nats. ⚠ A resolving-power floor at ``1/(2L)²`` was
#: derived, implemented and **REFUTED by its own insensitivity gate**: τ is CONTINUOUS across that
#: region on 4 of 5 ladder conditions (only the unstranded capture-OFF row is bimodal, and there the
#: two clusters are the silent strand arm and the live intron factory). With no empty interval, any
#: floor is a tuned constant. So the honest instrument reports the CURVE and lets the reader see how
#: much of the mass sits where nothing could be resolved.
SD_LAMBDA_DECADES = (1.0, 10.0, 100.0, 1000.0, np.inf)


def channel_masks(capture, chain, config) -> dict[str, np.ndarray]:
    """Which own-evidence channel(s) can speak at each slot — the solvability question, per channel.

    ⚠ **Overlapping capability flags, NOT a partition.** ``strand`` and ``factory`` are both live on a
    single-stranded intron node, because ``tau_lam`` is their SUM. Only ``none`` excludes the others.

    ⭐ Decomposed from ``tau_lam`` rather than re-derived: the factory arm is recoverable exactly by
    re-reading the captured ``intron_prior`` through the same ``density_factor_precision`` the solver
    used, and with the length channel off, whatever remains of ``tau_lam`` is the strand arm. So these
    are the solver's own numbers, split, not a second opinion about them.

    ⛔ **``locked`` IS THE TRAPS: no-magic-numbers CLASS ON BOTH AXES** (:func:`~rigel.calibration.node_geometry.g1_locked`).
    It used to be ``~solvable & (kind == NODE)``, which dropped every structurally-locked **edge** — an
    intergenic↔exon seam, where RNA cannot cross a gene boundary and ``_type_belief`` pins ``{0,0,1}``
    at ``Var(log f_g) = 0`` — into ``none``, i.e. excluded it from the scored population as honest
    ignorance. Those objects are the opposite of ignorant: they are structurally certain, and right.
    ⚠⚠ **This is NOT ``node_init.strand_evidence``'s ``struct_lock``, which is node-only ON PURPOSE** —
    that one governs whether a slot may EMIT certainty into its messages, and excludes G1 edges because
    a seam's crossing mass is RNA-contaminated. Two questions, one word; see ``g1_locked``.

    ⚠ **The τ tests here stay at the solver's own ``_EPS``**, so "has a channel" means what
    ``own_composition_logvar`` means by it. Strength is a separate question and it is reported as a
    curve over ``SD_LAMBDA_DECADES`` — see that constant for why a threshold cannot answer it.
    """
    tau = np.asarray(capture["_tau0_lam"], np.float64)
    # ⚠ ``chain`` is here to be CHECKED, not merely to be passed: a capture built against a different
    # partition would shift every mask by one slot, which is invisible in aggregate.
    if tau.shape != (int(chain.n_slots),):
        raise ValueError(
            f"capture['_tau0_lam'] has shape {tau.shape}; expected ({int(chain.n_slots)},), one per "
            f"chain slot. The capture and the chain describe different partitions."
        )
    # TRAPS: no-magic-numbers, from the ONE definition (`node_geometry.g1_locked`) — see there for why this is both axes
    # and why `node_init.struct_lock` is deliberately NOT.
    locked = g1_locked(capture["free_pos"], capture["free_neg"])
    lam_grid, _ = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
    fac = density_factor_precision(capture.get("intron_prior"), lam_grid)
    fac = np.zeros_like(tau) if fac is None else np.asarray(fac, np.float64)
    factory = (fac > _EPS) & ~locked
    strand = ((tau - fac) > _EPS) & ~locked
    return {
        "locked": locked,
        "strand": strand,
        "factory": factory,
        "none": ~(locked | strand | factory),
    }


def standardised_discrepancy(f_pred, f_true, var_log, fg_grid):
    """``z = (log f_pred − log f_truth) / sd``, both clipped to the solver's own grid support.

    ⛔ Log space, because ``var_gdna`` is ``Var(log f_g)`` and the names say otherwise. ⚠ ``sd == 0``
    is a slot the solver called CERTAIN; a wrong answer there is infinitely confident and is returned
    as ``inf`` rather than divided by zero — that is the honest encoding, and it puts structurally
    locked mistakes at the top of the ranking where they belong.
    """
    lo, hi = float(np.min(fg_grid)), float(np.max(fg_grid))
    lp = np.log(np.clip(np.asarray(f_pred, np.float64), lo, hi))
    lt = np.log(np.clip(np.asarray(f_true, np.float64), lo, hi))
    sd = np.sqrt(np.maximum(np.asarray(var_log, np.float64), 0.0))
    gap = lp - lt
    z = np.full(gap.shape, np.inf)
    np.divide(gap, sd, out=z, where=sd > _EPS)
    return np.where((sd <= _EPS) & (np.abs(gap) <= _EPS), 0.0, z), gap, sd


def audit(m, *, axis: str = "node", config=None) -> dict:
    """Partition one axis into undetermined / solvable-right / solvable-wrong and standardise."""
    config = config or CalibrationConfig()
    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    n_nodes, n_edges = int(m.payload.n_nodes), int(m.payload.n_edges)

    slots = channel_masks(cap, chain, config)
    per_axis = {
        name: P0._project(mask, chain, n_nodes, n_edges)[axis] for name, mask in slots.items()
    }

    g_p = np.asarray(getattr(m.arms["pass0"], f"mass_gdna_{axis}"), np.float64)
    r_p = np.asarray(getattr(m.arms["pass0"], f"mass_rna_{axis}"), np.float64)
    g_t = np.asarray(getattr(m.truth, f"mass_gdna_{axis}"), np.float64)
    r_t = np.asarray(getattr(m.truth, f"mass_rna_{axis}"), np.float64)
    total = g_t + r_t
    live = total > 0
    f_pred, f_true = P0.object_fractions(g_p, r_p)[0], P0.object_fractions(g_t, r_t)[0]

    # the solver's per-slot state, projected onto this axis
    def onto(values):
        out = np.zeros(n_nodes if axis == "node" else n_edges, np.float64)
        kind = np.asarray(chain.kind)
        obj = np.asarray(chain.obj_idx, np.int64)
        sel = (kind == NODE) if axis == "node" else (kind != NODE)
        out[obj[sel]] = np.asarray(values, np.float64)[sel]
        return out

    var_log = onto(cap["var_g"])
    _, fg_grid = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
    z, gap, sd = standardised_discrepancy(f_pred, f_true, var_log, fg_grid)

    determined = (per_axis["locked"] | per_axis["strand"] | per_axis["factory"]) & live
    undetermined = per_axis["none"] & live
    err = np.where(live, g_p - g_t, 0.0)
    # ⭐ the own-evidence STRENGTH each object earned, in the units the solver works in: sd(λ) = 1/√τ
    # nats. ``inf`` where there is no channel at all, 0 where the object is structurally certain.
    tau_axis = onto(cap["_tau0_lam"])
    with np.errstate(divide="ignore", invalid="ignore"):
        sd_lam = np.where(tau_axis > 0.0, 1.0 / np.sqrt(np.maximum(tau_axis, _EPS)), np.inf)
    sd_lam = np.where(per_axis["locked"], 0.0, sd_lam)
    return {
        "axis": axis,
        "live": live,
        "total": total,
        "err": err,
        "z": z,
        "gap": gap,
        "sd": sd,
        "sd_lam": sd_lam,
        "determined": determined,
        "undetermined": undetermined,
        "channels": per_axis,
        "ladder": {k: onto(cap[k]) for k in ("fg_strand", "fg_loc", "f_g")},
        "f_true": f_true,
        "f_pred": f_pred,
    }


#: ⭐ Read points on ``|f_pred − ½|``, which lives on [0, ½] by construction. Nothing branches on them.
OVERREACH_BANDS = (0.05, 0.15, 0.30, 0.50)


def undetermined_overreach_rows(a: dict) -> list[tuple]:
    """⛔⛔ **THE UNDETERMINED CLASS'S OWN GATE, AND IT WAS MISSING.**

    Excluding the undetermined population from the error denominator is right — an object with no own
    evidence reporting ``f_g ≈ ½`` at zero precision is stating a true fact about itself. But
    `SUCCESS.md` names the failure mode that exclusion leaves open, and nothing was checking it:

        "Its only failure mode is the opposite one: claiming a precision it has not earned."

    An undetermined object at ``f_g = 0.83`` is **not** honest ignorance. It is the relay and the
    population reference asserting an answer where the object had none, and because the class is
    excluded from every error total, the assertion is invisible. Measured cost of that blindness on the
    gDNA ladder: on ``gdna_g25_ss_0.50_nrna_none_capture_off``, **87 exon nodes are driven to
    ``f_g = 0.829`` against a truth of 0.009 — 395,251 fragments of error, 0.0 % of it scored**, and
    that is 3× the error on the row where the defect was first found.

    So: bucket the undetermined by how far from ½ they were moved, and report the error and the claimed
    precision in each bucket. ⭐ The correct answer for this class is ``½`` at ``sd = ∞``, so a row far
    from ½ with a FINITE sd is the defect, stated without any threshold deciding it.

    Returns ``(label, n, mass, Σ|err|, mean |f−½|, share claiming finite precision)``.
    """
    out = []
    und, fp, tot, err, sd = a["undetermined"], a["f_pred"], a["total"], np.abs(a["err"]), a["sd"]
    off = np.abs(np.asarray(fp, np.float64) - 0.5)
    lo = 0.0
    for hi in OVERREACH_BANDS:
        b = und & (off >= lo) & ((off < hi) if hi < 0.5 else (off <= hi))
        w = tot[b]
        out.append((
            f"{lo:.2f}–{hi:.2f}",
            int(b.sum()),
            float(w.sum()),
            float(err[b].sum()),
            float((off[b] * w).sum() / max(w.sum(), _EPS)),
            float(w[np.isfinite(sd[b]) & (sd[b] > _EPS)].sum() / max(w.sum(), _EPS)),
        ))
        lo = hi
    return out


def resolving_power_rows(a: dict, mask: np.ndarray) -> list[tuple]:
    """``(label, n, mass, Σ|err|, pred f_g, true f_g)`` per ``sd(λ)`` decade — the CURVE that replaces
    the solvable/undetermined cut for everything that is not structurally certain.

    ⭐ **This is the shape the measurement actually has.** ``sd(λ) ≫ L`` says the object's own evidence
    is flat across every λ the solver can represent, so whatever it reports came from its neighbours
    and the reference — the *substance* of "undetermined", stated as a magnitude instead of a class.
    ⛔ Nothing branches on the boundaries; they are decades, and a reader draws their own line.

    ⚠ **Pass a mask that EXCLUDES the structurally-locked slots.** Their ``sd(λ)`` is 0 because they
    are certain, not because their evidence is strong, and putting them in the first band would read
    as "0–1 nat: the best-evidenced population" when nothing was ever asked of them. They have their
    own row in the channel table.
    """
    out = []
    sd, err, total, fp, ft = a["sd_lam"], np.abs(a["err"]), a["total"], a["f_pred"], a["f_true"]
    lo = 0.0
    for hi in SD_LAMBDA_DECADES:
        # ⚠ the top band must be CLOSED at infinity: ``sd`` is exactly ``inf`` where no channel spoke,
        # and ``inf < inf`` is False, so a half-open top band silently drops that whole population —
        # the one the curve exists to make visible.
        b = mask & (sd >= lo) & ((sd < hi) if np.isfinite(hi) else np.ones_like(sd, bool))
        w = total[b]
        label = f"{lo:g}–{hi:g}" if np.isfinite(hi) else f">= {lo:g}"
        out.append((
            label,
            int(b.sum()),
            float(w.sum()),
            float(err[b].sum()),
            float((fp[b] * w).sum() / max(w.sum(), _EPS)),
            float((ft[b] * w).sum() / max(w.sum(), _EPS)),
        ))
        lo = hi
    return out


def _band_table(a: dict, mask: np.ndarray) -> list[tuple[str, int, float, float]]:
    """``(band, n, mass, Σ|err|)`` over ``|z|`` bands — the confidence profile of a population."""
    out = []
    z, err, total = np.abs(a["z"]), np.abs(a["err"]), a["total"]
    for lo, hi in Z_BANDS:
        b = mask & (z >= lo) & (z < hi)
        label = f"|z| {lo:g}–{hi:g}" if np.isfinite(hi) else f"|z| >= {lo:g}"
        out.append((label, int(b.sum()), float(total[b].sum()), float(err[b].sum())))
    return out


def report(m, a: dict, config=None) -> None:
    config = config or CalibrationConfig()
    _L = float(config.sweep_logodds_window)
    _kappa = float(m.arms["pass0"].rna_sense_frac)
    live, total, err = a["live"], a["total"], np.abs(a["err"])
    mass_all, err_all = total[live].sum(), err[live].sum()
    det, und = a["determined"], a["undetermined"]

    print()
    print("=" * 112)
    print(f"⭐⭐ SOLVABILITY AUDIT — {m.condition}   axis={a['axis']}   arm=pass-0 (prior-free)")
    print("=" * 112)
    print("   ⛔ Pass-0's job is to be a SUBSTRATE for the gDNA prior, not to be accurate everywhere.")
    print("      An object with no own evidence reporting f_g ~ 1/2 at zero precision is CORRECT.")
    print()
    print(f"   {'population':<34} {'objects':>9} {'mass':>14} {'share':>7} {'Σ|err|':>14} {'share':>7}")
    for label, mask in (("UNDETERMINED (excluded)", und), ("SOLVABLE (scored)", det)):
        print(f"   {label:<34} {int(mask.sum()):>9,} {total[mask].sum():>14,.0f} "
              f"{total[mask].sum() / max(mass_all, 1):>6.1%} {err[mask].sum():>14,.0f} "
              f"{err[mask].sum() / max(err_all, 1):>6.1%}")
    print()
    print("   by own-evidence CHANNEL (which can speak here at all). ⚠ THESE OVERLAP — a")
    print("   single-stranded intron node has both strand and factory — so mass DOUBLE-COUNTS here.")
    print("   The partition that adds up is the determined/undetermined split above.")
    for name in CHANNELS:
        c = a["channels"][name] & live
        print(f"     {name:<12} {int(c.sum()):>9,} objects  mass {total[c].sum():>13,.0f} "
              f"({total[c].sum() / max(mass_all, 1):>5.1%})   Σ|err| {err[c].sum():>13,.0f} "
              f"({err[c].sum() / max(err_all, 1):>5.1%})")

    print()
    print("   ⭐⭐ HOW STRONG IS THAT OWN EVIDENCE?  sd(λ) = 1/√τ nats, against a solver that can only")
    print(f"      represent λ ∈ [−{_L:g}, +{_L:g}].  ⛔ A row with sd(λ) far above {2 * _L:g} is scored")
    print("      as SOLVABLE and is not: its own evidence is flat over every λ the solver can express,")
    print("      so its answer came from neighbours and the reference. ⚠ NO threshold decides this —")
    print("      it is a curve, because τ is CONTINUOUS here and any cut would be a tuned constant.")
    print(f"   {'sd(λ) nats':<14} {'objects':>9} {'mass':>14} {'Σ|err|':>14} {'err share':>10} "
          f"{'pred f_g':>9} {'true f_g':>9}")
    _lock = a["channels"]["locked"] & live
    print(f"   {'CERTAIN (TRAPS: no-magic-numbers)':<14} {int(_lock.sum()):>9,} {total[_lock].sum():>14,.0f} "
          f"{err[_lock].sum():>14,.0f} {err[_lock].sum() / max(err[det].sum(), 1):>9.1%} "
          f"{'—':>9} {'—':>9}   structurally pure gDNA; nothing was asked of it")
    for label, n, mass, e, pred, true in resolving_power_rows(a, det & ~_lock):
        if n == 0:
            continue
        print(f"   {label:<14} {n:>9,} {mass:>14,.0f} {e:>14,.0f} "
              f"{e / max(err[det].sum(), 1):>9.1%} {pred:>9.4f} {true:>9.4f}")
    print(f"   ⚠ κ = {_kappa:.6f}, so the strand arm's information is scaled by (2κ−1)² = "
          f"{(2 * _kappa - 1) ** 2:.3e}.")
    print("      At κ = ½ it is EXACTLY zero (EQUATIONS §5.2); a fitted κ makes it merely tiny, and")
    print("      whether that lands above or below 1e-9 is what the old solvable/undetermined cut read.")

    print()
    print("   ⛔⛔ AND THE UNDETERMINED CLASS'S OWN FAILURE MODE — it is EXCLUDED from every error")
    print("      total above, so this is the only place it can be seen. Its correct answer is f_g = ½")
    print("      at sd = ∞; a row far from ½, and worse a row far from ½ CLAIMING precision, is the")
    print("      relay asserting an answer the object never had. ⚠ 0.0 % scored means 0.0 % reported.")
    print(f"   {'|f_pred − ½|':<14} {'objects':>9} {'mass':>14} {'Σ|err|':>14} {'mean |f−½|':>11} "
          f"{'claims sd':>10}")
    for label, n, mass, e, mean_off, prec_share in undetermined_overreach_rows(a):
        if n == 0:
            continue
        print(f"   {label:<14} {n:>9,} {mass:>14,.0f} {e:>14,.0f} {mean_off:>11.4f} "
              f"{prec_share:>9.1%}")

    print()
    print("   ⭐⭐ OF THE SOLVABLE OBJECTS — is the solver SURE when it is WRONG?")
    print(f"   {'band':<14} {'objects':>9} {'mass':>14} {'Σ|err|':>14} {'err share':>10}")
    det_err = err[det].sum()
    for label, n, mass, e in _band_table(a, det):
        print(f"   {label:<14} {n:>9,} {mass:>14,.0f} {e:>14,.0f} {e / max(det_err, 1):>9.1%}")
    print("   ⛔ The bottom rows are the defect: the solve is many sd from truth and SURE of it.")
    print("      A wrong value with a tight variance outvotes correct neighbours and anchors the prior.")

    print()
    print("   ⭐ IS THE DECLARED PRECISION EARNED?  realised RMS(log error) vs claimed sd, by decile")
    print("      of claimed sd. ⚠ Ratio ~1 = honest. >1 = OVERCONFIDENT. No threshold decides this.")
    sd, gap = a["sd"], a["gap"]
    ok = det & np.isfinite(gap) & (sd > _EPS)
    if int(ok.sum()) >= 10:
        q = np.quantile(sd[ok], np.linspace(0, 1, 11))
        print(f"     {'decile':<8} {'n':>8} {'claimed sd':>12} {'realised rms':>14} {'ratio':>8}")
        for i in range(10):
            b = ok & (sd >= q[i]) & (sd <= q[i + 1] if i == 9 else sd < q[i + 1])
            if not b.any():
                continue
            claimed = float(np.sqrt(np.mean(sd[b] ** 2)))
            realised = float(np.sqrt(np.mean(gap[b] ** 2)))
            print(f"     {i + 1:<8} {int(b.sum()):>8,} {claimed:>12.4f} {realised:>14.4f} "
                  f"{realised / max(claimed, _EPS):>8.2f}")
    else:
        print("     (too few solvable objects with a finite precision to form deciles)")

    print()
    print("   ⭐ THE ABLATION LADDER on the solvable set — which channel moved it off truth?")
    print("      strand-only -> local (strand+factory+reference) -> FINAL (after the relay)")
    lad, f_true = a["ladder"], a["f_true"]
    print(f"     {'rung':<28} {'mass-wtd |Δf_g| vs truth':>26}")
    for name, key in (("strand only", "fg_strand"), ("local (message-free)", "fg_loc"),
                      ("FINAL (with messages)", "f_g")):
        d = np.abs(lad[key] - f_true)
        w = total[det]
        print(f"     {name:<28} {float(np.sum(w * d[det]) / max(w.sum(), 1)):>26.4f}")
    print("   ⭐ A rung that is BETTER than the one below it means that channel HURT. Recovering the")
    print("      old `P1_OVERRULED` class: strand right and confident, then overruled by the full solve.")


#: ⭐ THE DEBUG CHAIN, in dependency order (owner, 2026-08-03). Each link is a precondition for the
#: next, so they are reported in this order and read top-down: the first one that is wrong explains
#: everything below it, and fixing anything lower before it is wasted work.
#:
#: 1. ``node/intergenic`` — structurally pure gDNA. If these are wrong, nothing else can be right.
#: 2. ``node/intron``     — the DENSITY PEEL: intron density against the intergenic background. On an
#:                          unstranded library this is the ONLY own-evidence channel there is.
#: 3. ``edge/intron|exon``      — must infer gDNA/RNA by propagation FROM the resolved intron node.
#: 4. ``edge/intergenic|exon``  — must impute from the resolved intergenic node.
CHAIN = (
    "node/intergenic",
    "node/intron",
    "node/exon",
    "edge/intergenic|intron",
    "edge/intron|exon",
    "edge/intergenic|exon",
    "edge/exon|exon",
    "edge/intron|intron",
    "edge/intergenic|intergenic",
)


def structural_classes(m, axis: str, config) -> dict[str, np.ndarray]:
    """Label each object by what it IS, structurally — a node's region type, or for a contiguous edge
    the PAIR of node types it separates.

    ⭐ The edge pair is the axis the debug chain turns on. An ``intron|exon`` line and an
    ``intergenic|exon`` line are the same kind of object to the solver and completely different
    problems: the first must inherit its answer from an intron node the density peel resolved, the
    second from a structurally-locked intergenic node. Lumping them as "edges" hides which
    propagation path is broken.

    ⚠ The pair is unordered (``intron|exon`` and ``exon|intron`` are one class) — the chain is
    genomic order, not a direction of inference.
    """
    from rigel.calibration.signature import coarse_type_array

    chain = m.debug_pass0["chain"]
    ra = m.debug_pass0["region_arrays"]
    names = {0: "intergenic", 1: "intron", 2: "exon"}
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n = int(m.payload.n_nodes) if axis == "node" else int(m.payload.n_edges)
    out = {c: np.zeros(n, bool) for c in CHAIN}

    if axis == "node":
        for t, nm in names.items():
            key = f"node/{nm}"
            if key in out:
                out[key][rtype == t] = True
        return out

    # a contiguous edge sits at chain slot s between node slots s-1 and s+1
    for s in np.flatnonzero(kind == EDGE):
        lo, hi = s - 1, s + 1
        if lo < 0 or hi >= kind.shape[0] or kind[lo] != NODE or kind[hi] != NODE:
            continue
        a, b = sorted((int(rtype[obj[lo]]), int(rtype[obj[hi]])))
        key = f"edge/{names[a]}|{names[b]}"
        if key in out:
            out[key][obj[s]] = True
    return out


def chain_report(m, a: dict, config) -> None:
    """The debug chain, in dependency order — the first broken link explains the ones below it."""
    struct = structural_classes(m, a["axis"], config)
    live, total, err = a["live"], a["total"], np.abs(a["err"])
    det, f_true, f_pred = a["determined"], a["f_true"], a["f_pred"]
    print()
    print("   ⭐⭐ THE DEBUG CHAIN — each link is a precondition for the next. Read TOP-DOWN: the first")
    print("      one that is wrong explains everything below it.")
    print(f"   {'structural class':<28} {'objects':>8} {'mass':>13} {'solv%':>6} "
          f"{'true f_g':>9} {'pred f_g':>9} {'mwae':>8} {'Σ|err|':>12}")
    for key in CHAIN:
        c = struct[key] & live
        if not c.any():
            continue
        w = total[c]
        tf = float(np.sum(w * f_true[c]) / w.sum())
        pf = float(np.sum(w * f_pred[c]) / w.sum())
        solv = float(total[c & det].sum() / w.sum())
        print(f"   {key:<28} {int(c.sum()):>8,} {w.sum():>13,.0f} {solv:>5.1%} "
              f"{tf:>9.4f} {pf:>9.4f} {err[c].sum() / w.sum():>8.4f} {err[c].sum():>12,.0f}")
    print("   ⚠ true/pred f_g are MASS-WEIGHTED means over the class, so they say whether the class is")
    print("     biased as a whole; mwae says whether its individual objects are right.")


def summarise(a: dict) -> dict:
    """The one-line-per-condition summary the debug loop's step 1-2 reads.

    ⭐ Every field is about the SOLVABLE population, because that is the only population pass-0 is
    accountable for. ``relay_delta`` is the one to watch across a panel: it is ``final − local`` on
    the solvable objects, so a POSITIVE value means the message layer moved objects that had their
    own answer AWAY from truth — the class the retired tooling called ``P1_OVERRULED``, and the one
    that matters out of proportion to its size because the hyperprior is fitted on exactly these
    objects.
    """
    det, total, err = a["determined"], a["total"], np.abs(a["err"])
    live = a["live"]
    mass_det = float(total[det].sum())
    conf = det & (np.abs(a["z"]) >= 2.0)
    sd, gap = a["sd"], a["gap"]
    ok = det & np.isfinite(gap) & (sd > _EPS)
    claimed = float(np.sqrt(np.mean(sd[ok] ** 2))) if ok.any() else float("nan")
    realised = float(np.sqrt(np.mean(gap[ok] ** 2))) if ok.any() else float("nan")
    w, f_true = total[det], a["f_true"]

    def rung(key):
        return float(np.sum(w * np.abs(a["ladder"][key][det] - f_true[det])) / max(w.sum(), 1))

    local, final = rung("fg_loc"), rung("f_g")
    # ⭐⭐⭐ THE TWO FIXED-DENOMINATOR FIELDS — TRAPS: honesty-metrics-reward-ignorance's own prescription, and the reason this
    # table can be ranked on at all. Every field above is defined over the DETERMINED set, whose size
    # the solver moves by declining to answer; these two are defined over the LIVE population, so
    # nothing the solver does to its own confidence can touch them.
    # ⛔ They exist because the boolean `determined` flips on FITTING NOISE at some conditions and not
    # others (TRAPS: deadband-from-the-wrong-sample): `g75 ss0.50 capture_off` reports solv% 96.2 %, weak% 65.4 %, conf-wrong 92,154
    # and calib 3.67 — and forcing κ = ½ takes solv% to 0.0 % while Σ|err| moves −0.6/−5.3/+4.6 %
    # across the three big classes. The row is a reporting artefact, and ranking the ladder on the
    # gameable columns picked it as the panel's worst. These two would not have.
    mass_live = float(total[live].sum())
    return_extra = {
        "all_mwae": float(err[live].sum()) / max(mass_live, 1.0),
        "abs_err": float(err[live].sum()),
    }
    # ⭐⭐ THE COMPANION COLUMN THAT MAKES ``solv%`` SAFE TO READ. ``solv%`` counts objects the SOLVER
    # treats as evidenced (``tau > 1e-9``), and that admits a strand arm whose own statement is 10³ nats
    # wide — statistically real, physically nil (TRAPS: a-threshold-on-a-fitted-residue). So report, beside it, the share of the
    # scored error that sits on objects whose own evidence cannot resolve one nat in ten. ⚠ 10 is a
    # DECADE off the curve, a read point, not a cut: nothing branches on it, and the curve above it in
    # the single-condition report is what a reader should actually consult.
    weak = det & (a["sd_lam"] >= 10.0)
    return {
        "solvable_mass_share": mass_det / max(float(total[live].sum()), 1.0),
        "weak_evidence_err_share": float(err[weak].sum()) / max(float(err[det].sum()), 1.0),
        "weak_evidence_mass_share": float(total[weak].sum()) / max(mass_det, 1.0),
        "solvable_err_share": float(err[det].sum()) / max(float(err[live].sum()), 1.0),
        "solvable_mwae": float(err[det].sum()) / max(mass_det, 1.0),
        "conf_wrong_objects": int(conf.sum()),
        "conf_wrong_err": float(err[conf].sum()),
        "conf_wrong_err_share": float(err[conf].sum()) / max(float(err[det].sum()), 1.0),
        "calibration_ratio": realised / claimed if claimed and claimed > 0 else float("nan"),
        "local_mwae": local,
        "final_mwae": final,
        "relay_delta": final - local,
        **return_extra,
    }


def panel_report(rows: list[tuple[str, float, dict]]) -> None:
    print()
    print("=" * 124)
    print("⭐⭐ THE PANEL, SCORED ON THE SOLVABLE POPULATION ONLY (pass-0)")
    print("=" * 124)
    print("   ⛔ Objects with no own evidence are EXCLUDED: in pass-0 they are correctly saying they")
    print("      cannot be solved without a prior. Scoring them buries everything that matters.")
    print("   ⭐ relay_delta = final − local on the solvable set. POSITIVE means the message layer")
    print("      moved objects that HAD their own answer away from truth.")
    print("   ⛔⛔ READ `weak%` BEFORE `mwae`. `solv%` counts what the SOLVER treats as evidenced, and")
    print("      that admits a strand arm whose own statement is 10³ nats wide against a ±10-nat grid")
    print("      (TRAPS a-threshold-on-a-fitted-residue). `weak%` is the share of the scored ERROR sitting on objects with")
    print("      sd(λ) ≥ 10 nats — i.e. on objects that had no answer of their own after all. A row")
    print("      with weak% near 100 is reporting the relay and the reference, not a solve.")
    print()
    print("   ⭐⭐ AND RANK ON THE LAST TWO, NOT ON `solv%`/`mwae`/`conf-wrong`/`calib`. Those four")
    print("      share a denominator the SOLVER moves — `determined` is a boolean on a continuous τ,")
    print("      and it flips on fitting noise (TRAPS deadband-from-the-wrong-sample). `mwae_all` and `Σ|err|` are over every")
    print("      LIVE object, so nothing the solver does to its own confidence can touch them.")
    print(f"   {'condition':<46} {'f_gdna':>7} {'solv%':>6} {'weak%':>6} {'mwae':>7} "
          f"{'conf-wrong':>11} {'calib':>6} {'local':>7} {'final':>7} {'relay Δ':>9} "
          f"{'mwae_all':>9} {'Σ|err|':>11}")
    print("   " + "-" * 125)
    for name, truth, s in rows:
        print(
            f"   {name:<46} {truth:>7.4f} {s['solvable_mass_share']:>5.1%} "
            f"{s['weak_evidence_err_share']:>5.1%} "
            f"{s['solvable_mwae']:>7.4f} {s['conf_wrong_err']:>11,.0f} "
            f"{s['calibration_ratio']:>6.2f} {s['local_mwae']:>7.4f} {s['final_mwae']:>7.4f} "
            f"{s['relay_delta']:>+9.4f} {s['all_mwae']:>9.4f} {s['abs_err']:>11,.0f}"
        )
    # ⛔ ZERO-gDNA ROWS ARE NEVER AVERAGED IN. Truth is 0 exactly there, so every log-space
    # discrepancy is measured against the grid floor and any change that lowers the estimate
    # "improves" the row — a one-sidedness that has already reversed a verdict in this project. They
    # are printed above as false-positive checks and excluded from every aggregate below.
    scored = [r for r in rows if r[1] > 0.0]
    zero = [r for r in rows if r[1] <= 0.0]
    if zero:
        print(f"\n   ⚠ {len(zero)} zero-gDNA row(s) shown above are FALSE-POSITIVE CHECKS and are")
        print("     excluded from the aggregates (truth = 0 exactly ⇒ the comparison is one-sided).")
    if not scored:
        return
    hurt = [r for r in scored if r[2]["relay_delta"] > 0]
    print(f"\n   ⭐ the relay HURTS the solvable set on {len(hurt)}/{len(scored)} CONTAMINATED "
          f"conditions   (mean Δ {np.mean([r[2]['relay_delta'] for r in scored]):+.4f})")
    over = [r for r in scored if r[2]["calibration_ratio"] > 1.0]
    print(f"   ⭐ the declared precision is NOT earned (ratio > 1) on {len(over)}/{len(scored)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--condition", default=None, help="one condition; omit for the whole panel")
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--axis", default="node", choices=("node", "edge", "both"))
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")))
    ap.add_argument("--oracle-cache", type=Path, default=None)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig()
    names = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file())
    )
    if not names:
        print(f"no conditions with a sim_oracle.bam under {args.suite}", file=sys.stderr)
        return 2

    panel: list[tuple[str, float, dict]] = []
    for name in names:
        cond = args.suite / name
        truth = P0.truth_f_gdna(cond) or 0.0
        print(f"  {name} …", flush=True)
        m = P0.measure_condition(
            bam=str(cond / "sim_oracle.bam"), index=index, pipeline_config=PipelineConfig(),
            calibration_config=config, work_dir=args.work_dir / "rigel_pass0_oracle", tag=name,
            truth_pmfs=lambda size, d=cond: (
                P0.truth_length_pmf(d, "gdna", size), P0.truth_length_pmf(d, "rna", size)
            ),
            oracle_cache=args.oracle_cache,
        )
        for axis in (("node", "edge") if args.axis == "both" else (args.axis,)):
            a = audit(m, axis=axis, config=config)
            if len(names) == 1:
                report(m, a, config)
                chain_report(m, a, config)
            if axis == ("node" if args.axis != "edge" else "edge"):
                panel.append((name, truth, summarise(a)))
    if len(panel) > 1:
        panel_report(panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
