#!/usr/bin/env python
"""Which objects are solvable, which are solved wrong, and which are confidently wrong?

Pass-0 is the prior-free solve, and its job is to produce the substrate the gDNA prior is fitted on,
not to be accurate everywhere: an object with no own evidence that reports ``f_g ≈ ½`` at zero
precision is correctly saying it cannot be solved without a prior, and counting that as error buries
what matters. So this audit takes `pass0_vs_oracle.measure_condition`'s truth and pass-0 arm and
separates three populations a mass-weighted error lumps together: undetermined (no own-evidence
channel; excluded from the error denominator, and checked only for the opposite failure, a value far
from ½ claiming a finite precision), solvable and right, and solvable and wrong, the last split by
confidence, because a wrong value with a tight variance outvotes correct neighbours and anchors the
prior. The standardised discrepancy is in log space, ``z = (log f_pred − log f_true) / sqrt(var)``,
because ``var_gdna`` is ``Var(log f_g)`` despite its name, with the truth clipped to the solver's own
λ-grid endpoints. No tuned threshold decides anything: ``|z|`` bands are sd multiples, own-evidence
strength is a curve over ``sd(λ) = 1/√τ`` decades, and the headline is a calibration curve of realised
RMS log error against claimed sd. The ablation ladder (``fg_strand`` / ``fg_loc`` / ``f_g``) is read
off the solver's own capture, so nothing is re-solved. Rank a panel on ``mwae_all`` and ``Σ|err|``
(fixed denominators), never on the columns whose denominator the solver moves; zero-gDNA rows are
printed as false-positive checks and never averaged in.

Also a library: `vertex_ceiling.py` and the calibration tests import `audit`, `summarise`,
`channel_masks`, `standardised_discrepancy`, `resolving_power_rows`, `undetermined_overreach_rows`
and the band constants.

Usage::

    python scripts/design/solvability_audit.py --condition <name> --oracle-cache <dir>   # one condition, in full
    python scripts/design/solvability_audit.py --oracle-cache <dir>                      # the whole panel, one row each
    python scripts/design/solvability_audit.py --condition <name> --axis both
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


from _shared import sibling  # noqa: E402


P0 = sibling("pass0_vs_oracle.py")

from rigel.calibration.calibrate import lattice_points  # noqa: E402
from rigel.calibration.density_deconv import density_factor_precision  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.region_geometry import g1_locked  # noqa: E402
from rigel.calibration.simplex_logodds import _logodds_grid  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1.0e-9

#: Standard-deviation multiples, not tuned constants: ``z`` carries its own scale, and these are
#: where a Gaussian's mass sits. Nothing branches on them; they are read points on a distribution.
Z_BANDS = ((0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, np.inf))

#: The own-evidence channels a slot can have in pass-0 (C_info's identification test is a statement
#: about the payload, not about the solver, and is not one of them). These overlap and are not a
#: partition: ``tau_lam = i_strand + tau_fac``, so a single-stranded intron region has both, and that
#: is the best-evidenced kind of object there is. Only ``none`` is exclusive of the rest. Mass summed
#: over these rows therefore double-counts, and the report says so; the partition that does add up
#: is determined / undetermined.
CHANNELS = ("locked", "strand", "factory", "none")


#: Own-evidence strength is a continuum and is reported as one. ``tau_lam`` is a Fisher precision on
#: ``λ = log(f_g/f_R)``, so a slot's own statement carries sd ``1/√τ``, in nats, against a solver that
#: can only represent ``λ ∈ [−L, +L]``. These are decade boundaries on that sd, read points and not
#: thresholds. A binary solvable/undetermined cut at ``τ > 1e-9`` is the solver's gate, not a
#: strength: the strand arm's information is ``∝ (2κ−1)²``, exactly zero at κ = ½, but κ is fitted,
#: so on an unstranded library τ lands near 1e-7 rather than 0 and the object is scored as solvable
#: while its own statement has an sd of thousands of nats. A resolving-power floor was refuted by its
#: own insensitivity gate (τ is continuous across that region, so any floor is a tuned constant); the
#: curve lets the reader see how much of the mass sits where nothing could be resolved.
SD_LAMBDA_DECADES = (1.0, 10.0, 100.0, 1000.0, np.inf)


def channel_masks(capture, chain, config) -> dict[str, np.ndarray]:
    """Which own-evidence channel(s) can speak at each slot: the solvability question, per channel.

    Overlapping capability flags, not a partition: ``strand`` and ``factory`` are both live on a
    single-stranded intron region, because ``tau_lam`` is their sum; only ``none`` excludes the others.
    Decomposed from ``tau_lam`` rather than re-derived: the factory arm is recovered exactly by
    re-reading the captured ``intron_prior`` through the same ``density_factor_precision`` the solver
    used, and whatever remains of ``tau_lam`` is the strand arm, so these are the solver's own numbers
    split, not a second opinion. ``locked`` is the G1 class on both axes
    (:func:`~rigel.calibration.region_geometry.g1_locked`): a structurally-locked boundary is certain
    and right, not ignorant, and must not fall into ``none``. The τ tests stay at the solver's own
    ``_EPS``, so "has a channel" means what ``region_init.has_own_composition_evidence`` means by it;
    strength is a separate question, reported as the curve over ``SD_LAMBDA_DECADES``.
    """
    tau = np.asarray(capture.tau_lam, np.float64)
    # ``chain`` is here to be checked, not merely passed: a capture built against a different
    # partition would shift every mask by one slot, which is invisible in aggregate.
    if tau.shape != (int(chain.n_slots),):
        raise ValueError(
            f"capture.tau_lam has shape {tau.shape}; expected ({int(chain.n_slots)},), one per "
            f"chain slot. The capture and the chain describe different partitions."
        )
    # G1, from the one definition (`region_geometry.g1_locked`), on both axes.
    locked = g1_locked(capture.free_pos, capture.free_neg)
    lam_grid, _ = _logodds_grid(
        lattice_points(config.sweep_logodds_window, config.sweep_logodds_step),
        float(config.sweep_logodds_window),
    )
    fac = density_factor_precision(capture.intron_prior, lam_grid)
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
    """``z = (log f_pred − log f_truth) / sd``, both clipped to the solver's own grid support; returns
    ``(z, gap, sd)``.

    Log space, because ``var_gdna`` is ``Var(log f_g)`` and the names say otherwise. ``sd == 0`` is a
    slot the solver called certain; a wrong answer there is infinitely confident and is returned as
    ``inf`` rather than divided by zero, which puts structurally locked mistakes at the top of the
    ranking where they belong.
    """
    lo, hi = float(np.min(fg_grid)), float(np.max(fg_grid))
    lp = np.log(np.clip(np.asarray(f_pred, np.float64), lo, hi))
    lt = np.log(np.clip(np.asarray(f_true, np.float64), lo, hi))
    sd = np.sqrt(np.maximum(np.asarray(var_log, np.float64), 0.0))
    gap = lp - lt
    z = np.full(gap.shape, np.inf)
    np.divide(gap, sd, out=z, where=sd > _EPS)
    return np.where((sd <= _EPS) & (np.abs(gap) <= _EPS), 0.0, z), gap, sd


def audit(m, *, axis: str = "region", config=None) -> dict:
    """Partition one axis into undetermined / solvable-right / solvable-wrong and standardise."""
    config = config or CalibrationConfig()
    cap, chain = m.debug_pass0["capture"], m.debug_pass0["chain"]
    n_regions, n_boundaries = int(m.payload.n_regions), int(m.payload.n_boundaries)

    slots = channel_masks(cap, chain, config)
    per_axis = {
        name: P0._project(mask, chain, n_regions, n_boundaries)[axis] for name, mask in slots.items()
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
        out = np.zeros(n_regions if axis == "region" else n_boundaries, np.float64)
        kind = np.asarray(chain.kind)
        obj = np.asarray(chain.obj_idx, np.int64)
        sel = (kind == REGION) if axis == "region" else (kind != REGION)
        out[obj[sel]] = np.asarray(values, np.float64)[sel]
        return out

    var_log = onto(cap.var_g)
    _, fg_grid = _logodds_grid(
        lattice_points(config.sweep_logodds_window, config.sweep_logodds_step),
        float(config.sweep_logodds_window),
    )
    z, gap, sd = standardised_discrepancy(f_pred, f_true, var_log, fg_grid)

    determined = (per_axis["locked"] | per_axis["strand"] | per_axis["factory"]) & live
    undetermined = per_axis["none"] & live
    err = np.where(live, g_p - g_t, 0.0)
    # the own-evidence strength each object earned, in the units the solver works in: sd(λ) = 1/√τ
    # nats. ``inf`` where there is no channel at all, 0 where the object is structurally certain.
    tau_axis = onto(cap.tau_lam)
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
        "ladder": {k: onto(getattr(cap, k)) for k in ("fg_strand", "fg_loc", "f_g")},
        "f_true": f_true,
        "f_pred": f_pred,
    }


#: Read points on ``|f_pred − ½|``, which lives on [0, ½] by construction. Nothing branches on them.
OVERREACH_BANDS = (0.05, 0.15, 0.30, 0.50)


def undetermined_overreach_rows(a: dict) -> list[tuple]:
    """The undetermined class's own gate.

    Excluding the undetermined population from the error denominator is right: an object with no own
    evidence reporting ``f_g ≈ ½`` at zero precision is stating a true fact about itself. Its one
    failure mode is the opposite one, claiming a precision it has not earned: an undetermined object
    far from ½ is the messages and the reference asserting an answer where the object had none, and
    because the class is excluded from every error total the assertion would otherwise be invisible.
    So the undetermined are bucketed by how far from ½ they were moved, with the error and the claimed
    precision in each bucket; the correct answer for this class is ``½`` at ``sd = ∞``, so a row far
    from ½ with a finite sd is the defect, stated without any threshold deciding it.

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
    """``(label, n, mass, Σ|err|, pred f_g, true f_g)`` per ``sd(λ)`` decade: the curve that replaces
    a solvable/undetermined cut for everything that is not structurally certain.

    ``sd(λ) ≫ L`` says the object's own evidence is flat across every λ the solver can represent, so
    whatever it reports came from its neighbours and the reference: the substance of "undetermined",
    stated as a magnitude instead of a class. Nothing branches on the boundaries. Pass a mask that
    excludes the structurally-locked slots: their ``sd(λ)`` is 0 because they are certain, not because
    their evidence is strong, and they have their own row in the channel table.
    """
    out = []
    sd, err, total, fp, ft = a["sd_lam"], np.abs(a["err"]), a["total"], a["f_pred"], a["f_true"]
    lo = 0.0
    for hi in SD_LAMBDA_DECADES:
        # the top band must be closed at infinity: ``sd`` is exactly ``inf`` where no channel spoke,
        # and ``inf < inf`` is False, so a half-open top band would silently drop that whole
        # population, the one the curve exists to make visible.
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
    """``(band, n, mass, Σ|err|)`` over ``|z|`` bands: the confidence profile of a population."""
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
    print("   single-stranded intron region has both strand and factory — so mass DOUBLE-COUNTS here.")
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
    print("      it is a curve, because τ is CONTINUOUS here and any region_bound would be a tuned constant.")
    print(f"   {'sd(λ) nats':<14} {'objects':>9} {'mass':>14} {'Σ|err|':>14} {'err share':>10} "
          f"{'pred f_g':>9} {'true f_g':>9}")
    _lock = a["channels"]["locked"] & live
    print(f"   {'CERTAIN (G1)':<14} {int(_lock.sum()):>9,} {total[_lock].sum():>14,.0f} "
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
    print("      whether that lands above or below 1e-9 is what the old solvable/undetermined region_bound read.")

    print()
    print("   ⛔⛔ AND THE UNDETERMINED CLASS'S OWN FAILURE MODE — it is EXCLUDED from every error")
    print("      total above, so this is the only place it can be seen. Its correct answer is f_g = ½")
    print("      at sd = ∞; a row far from ½, and worse a row far from ½ CLAIMING precision, is the")
    print("      messages asserting an answer the object never had. ⚠ 0.0 % scored means 0.0 % reported.")
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
    print("      strand-only -> local (strand+factory+reference) -> FINAL (after the messages)")
    lad, f_true = a["ladder"], a["f_true"]
    print(f"     {'rung':<28} {'mass-wtd |Δf_g| vs truth':>26}")
    for name, key in (("strand only", "fg_strand"), ("local (message-free)", "fg_loc"),
                      ("FINAL (with messages)", "f_g")):
        d = np.abs(lad[key] - f_true)
        w = total[det]
        print(f"     {name:<28} {float(np.sum(w * d[det]) / max(w.sum(), 1)):>26.4f}")
    print("   ⭐ A rung that is BETTER than the one below it means that channel HURT. Recovering the")
    print("      old `P1_OVERRULED` class: strand right and confident, then overruled by the full solve.")


#: The debug chain, in dependency order. Each link is a precondition for the next, so they are
#: reported in this order and read top-down: the first one that is wrong explains everything below
#: it, and fixing anything lower before it is wasted work.
#:
#: 1. ``region/intergenic``: structurally pure gDNA. If these are wrong, nothing else can be right.
#: 2. ``region/intron``: the density deconvolution, intron density against the intergenic
#:    background. On an unstranded library this is the only own-evidence channel there is.
#: 3. ``boundary/intron|exon``: must infer gDNA/RNA by propagation from the resolved intron region.
#: 4. ``boundary/intergenic|exon``: must impute from the resolved intergenic region.
CHAIN = (
    "region/intergenic",
    "region/intron",
    "region/exon",
    "boundary/intergenic|intron",
    "boundary/intron|exon",
    "boundary/intergenic|exon",
    "boundary/exon|exon",
    "boundary/intron|intron",
    "boundary/intergenic|intergenic",
)


def structural_classes(m, axis: str, config) -> dict[str, np.ndarray]:
    """Label each object by what it is structurally: a region's region type, or for a contiguous
    boundary the pair of region types it separates.

    The boundary pair is the axis the debug chain turns on: an ``intron|exon`` boundary and an
    ``intergenic|exon`` boundary are the same kind of object to the solver and completely different
    problems (the first inherits its answer from an intron region the density deconvolve resolved,
    the second from a structurally-locked intergenic region), and lumping them as "boundaries" hides
    which propagation path is broken. The pair is unordered (``intron|exon`` and ``exon|intron`` are
    one class); the chain is genomic order, not a direction of inference.
    """
    from rigel.calibration.signature import coarse_type_array

    chain = m.debug_pass0["chain"]
    ra = m.debug_pass0["region_arrays"]
    names = {0: "intergenic", 1: "intron", 2: "exon"}
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n = int(m.payload.n_regions) if axis == "region" else int(m.payload.n_boundaries)
    out = {c: np.zeros(n, bool) for c in CHAIN}

    if axis == "region":
        for t, nm in names.items():
            key = f"region/{nm}"
            if key in out:
                out[key][rtype == t] = True
        return out

    # a contiguous boundary sits at chain slot s between region slots s-1 and s+1
    for s in np.flatnonzero(kind == BOUNDARY):
        lo, hi = s - 1, s + 1
        if lo < 0 or hi >= kind.shape[0] or kind[lo] != REGION or kind[hi] != REGION:
            continue
        a, b = sorted((int(rtype[obj[lo]]), int(rtype[obj[hi]])))
        key = f"boundary/{names[a]}|{names[b]}"
        if key in out:
            out[key][obj[s]] = True
    return out


def chain_report(m, a: dict, config) -> None:
    """The debug chain, in dependency order: the first broken link explains the ones below it."""
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
    """The one-line-per-condition summary the panel table reads.

    Every field but the last two is about the solvable population, the only population pass-0 is
    accountable for. ``message_delta`` is ``final − local`` on the solvable objects, so a positive
    value means the message layer moved objects that had their own answer away from truth, which
    matters out of proportion to its size because the hyperprior is fitted on exactly these objects.
    ``all_mwae`` and ``abs_err`` are over the live population and are the fields to rank on.
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
    # the two fixed-denominator fields (TRAPS: honesty-metrics-reward-ignorance), the reason this
    # table can be ranked on at all. Every field above is defined over the determined set, whose size
    # the solver moves by declining to answer, and the boolean `determined` flips on fitting noise
    # (TRAPS: deadband-from-the-wrong-sample); these two are defined over the live population, so
    # nothing the solver does to its own confidence can touch them.
    mass_live = float(total[live].sum())
    return_extra = {
        "all_mwae": float(err[live].sum()) / max(mass_live, 1.0),
        "abs_err": float(err[live].sum()),
    }
    # the companion column that makes ``solv%`` safe to read: ``solv%`` counts objects the solver
    # treats as evidenced (``tau > 1e-9``), which admits a strand arm whose own statement is thousands
    # of nats wide (TRAPS: a-threshold-on-a-fitted-residue). So report, beside it, the share of the
    # scored error on objects whose own evidence cannot resolve one nat in ten. 10 is a decade off the
    # curve, a read point and not a threshold; the curve in the single-condition report is what a
    # reader should consult.
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
        "message_delta": final - local,
        **return_extra,
    }


def panel_report(rows: list[tuple[str, float, dict]]) -> None:
    print()
    print("=" * 124)
    print("⭐⭐ THE PANEL, SCORED ON THE SOLVABLE POPULATION ONLY (pass-0)")
    print("=" * 124)
    print("   ⛔ Objects with no own evidence are EXCLUDED: in pass-0 they are correctly saying they")
    print("      cannot be solved without a prior. Scoring them buries everything that matters.")
    print("   ⭐ message_delta = final − local on the solvable set. POSITIVE means the message layer")
    print("      moved objects that HAD their own answer away from truth.")
    print("   ⛔⛔ READ `weak%` BEFORE `mwae`. `solv%` counts what the SOLVER treats as evidenced, and")
    print("      that admits a strand arm whose own statement is 10³ nats wide against a ±10-nat grid")
    print("      (TRAPS a-threshold-on-a-fitted-residue). `weak%` is the share of the scored ERROR sitting on objects with")
    print("      sd(λ) ≥ 10 nats — i.e. on objects that had no answer of their own after all. A row")
    print("      with weak% near 100 is reporting the messages and the reference, not a solve.")
    print()
    print("   ⭐⭐ AND RANK ON THE LAST TWO, NOT ON `solv%`/`mwae`/`conf-wrong`/`calib`. Those four")
    print("      share a denominator the SOLVER moves — `determined` is a boolean on a continuous τ,")
    print("      and it flips on fitting noise (TRAPS deadband-from-the-wrong-sample). `mwae_all` and `Σ|err|` are over every")
    print("      LIVE object, so nothing the solver does to its own confidence can touch them.")
    print(f"   {'condition':<46} {'f_gdna':>7} {'solv%':>6} {'weak%':>6} {'mwae':>7} "
          f"{'conf-wrong':>11} {'calib':>6} {'local':>7} {'final':>7} {'msg Δ':>9} "
          f"{'mwae_all':>9} {'Σ|err|':>11}")
    print("   " + "-" * 125)
    for name, truth, s in rows:
        print(
            f"   {name:<46} {truth:>7.4f} {s['solvable_mass_share']:>5.1%} "
            f"{s['weak_evidence_err_share']:>5.1%} "
            f"{s['solvable_mwae']:>7.4f} {s['conf_wrong_err']:>11,.0f} "
            f"{s['calibration_ratio']:>6.2f} {s['local_mwae']:>7.4f} {s['final_mwae']:>7.4f} "
            f"{s['message_delta']:>+9.4f} {s['all_mwae']:>9.4f} {s['abs_err']:>11,.0f}"
        )
    # zero-gDNA rows are never averaged in. Truth is 0 exactly there, so every log-space discrepancy
    # is measured against the grid floor and any change that lowers the estimate "improves" the row.
    # They are printed above as false-positive checks and excluded from every aggregate below.
    scored = [r for r in rows if r[1] > 0.0]
    zero = [r for r in rows if r[1] <= 0.0]
    if zero:
        print(f"\n   ⚠ {len(zero)} zero-gDNA row(s) shown above are FALSE-POSITIVE CHECKS and are")
        print("     excluded from the aggregates (truth = 0 exactly ⇒ the comparison is one-sided).")
    if not scored:
        return
    hurt = [r for r in scored if r[2]["message_delta"] > 0]
    print(f"\n   ⭐ the messages HURT the solvable set on {len(hurt)}/{len(scored)} CONTAMINATED "
          f"conditions   (mean Δ {np.mean([r[2]['message_delta'] for r in scored]):+.4f})")
    over = [r for r in scored if r[2]["calibration_ratio"] > 1.0]
    print(f"   ⭐ the declared precision is NOT earned (ratio > 1) on {len(over)}/{len(scored)}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--condition", default=None, help="one condition; omit for the whole panel")
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--axis", default="region", choices=("region", "boundary", "both"))
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
        for axis in (("region", "boundary") if args.axis == "both" else (args.axis,)):
            a = audit(m, axis=axis, config=config)
            if len(names) == 1:
                report(m, a, config)
                chain_report(m, a, config)
            if axis == ("region" if args.axis != "boundary" else "boundary"):
                panel.append((name, truth, summarise(a)))
    if len(panel) > 1:
        panel_report(panel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
