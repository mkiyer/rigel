"""WHICH EXONS ARE SOLVED ACCURATELY — AND WHAT, WITHOUT TRUTH, TELLS YOU WHICH?

⛔⛔⛔ **UNDER REVIEW, AND ITS DOCSTRING VERDICT IS NOT YET REPRODUCIBLE FROM THIS FILE — 2026-08-20.**
An adversarial re-run of this instrument's own claims found, and reproduced, the following. Read them
BEFORE quoting any number from here:

* ⛔ **The headline pooled table was produced by a scratch script, not by this instrument.**
  :func:`report_summary` pools only ``flank structure`` and ``hops to a pure-gDNA slot``; every pooled
  strand / own-evidence / depth-quartile number came from an ungated helper. **Running this file as
  documented does not reproduce its own headline.**
* ⛔ **"depth is the only predictor monotone on all three in-scope strata" is FALSE** — the
  unstranded × capture-OFF quartiles are 6.9 / 18.6 / 18.2 / 37.1 and were printed with two of them
  transposed.
* ⛔ **"roughly independent of the strand bit" is FALSE and nothing measured it** — computed per
  quartile, the AMBIG penalty runs 2.5x → 10x → 15x → 1.4x. There is a strong interaction, and at the
  deepest quartile the strand bit is worth almost nothing, so "deep AND single-strand" discards the
  deepest AMBIG mass for no measured gain.
* ⛔ **The "phantom own-evidence mask, 7.6x worse" finding is a cross-condition artefact** — WITHIN
  `g98 ss0.50 capture_off` the evidence-bearing slots are BETTER (113.2 vs 168.3 per 1k). The
  comparator pooled two conditions whose mask is empty.
* ⛔ **Claims resting on the `g00` row pool ACROSS strata and reverse per condition** (both the
  ``sj only`` and the ``unreachable`` findings flip on at least one control condition).
* ⛔ **One self-test hole**: making ``quartile_bucket``'s edges RANK-based instead of MASS-based leaves
  ``--self-test`` green, so the equal-mass property the depth columns rest on is asserted but not gated.

⭐ **What survived verification**: the per-condition numeric columns reproduce exactly, and the
strand-bit finding itself — AMBIG exons carry ~18-20 % of exon mass and 36-50 % of exon misplacement —
holds per condition on the in-scope strata. ⚠ Everything else above must be re-derived before use.


⭐⭐⭐ **THE QUESTION IS THE TRAINING SUBSTRATE FOR THE gDNA LANDSCAPE PRIOR, AND IT IS CIRCULAR BY
CONSTRUCTION.** A per-object prior can only be trained on the exons pass-0 already solves, and pass-0
only solves an exon where something supplies its composition. So the deliverable is not "how big is the
exon error" — it is **which OBSERVABLE class of exon is safe to train on**, in FRAGMENTS, with every
predictor knowable WITHOUT truth, because on real data there is none.

⛔ **NO THRESHOLD ANYWHERE.** "Solved well" is the mass-weighted CURVE — the fraction of exon mass whose
``|f_g - true_f_g|`` is at most ``x``, over a grid of ``x`` — plus mass-weighted QUANTILES, which need no
grid at all, plus ``Sum|f_g - true_f_g|*M`` in FRAGMENTS. Nothing in this file branches on an error
value; the curve IS the answer (`solvability_audit.py`'s threshold-free calibration curve is the
precedent).

**THE FIVE TRUTH-FREE PREDICTORS, all read off the payload and the annotation:**

=============================  ==================================================================
**own composition evidence**   :func:`~rigel.calibration.region_init.has_own_composition_evidence`
                               on the slot's own ``_tau0_lam`` — the solver's OWN predicate,
                               imported rather than restated
**RNA strands admissible**     AMBIG (both admissible) vs single-strand — AXIOM 0's two bits
**own observed depth**         the slot's own count, bucketed by its mass-weighted QUARTILES.
                               Data-derived edges, so no constant enters — and the quartiles are
                               equal-MASS by construction, so their misplaced-fragment columns are
                               directly comparable
**hops to a pure-gDNA slot**   chain distance to the nearest slot where NEITHER RNA strand is
                               admissible, i.e. where gDNA is measured rather than imputed.
                               Reuses ``hop_currency.depth_to_measured``
**flank structure**            SPLICE SITE flank, TERMINUS flank, both, or a reference terminal
                               (``is_splice_site`` / ``is_terminus``, either strand)
=============================  ==================================================================

⭐⭐⭐ **THE VERDICT — 16/16 ladder conditions and 8/8 test-chromosome ones, 2026-08-20, shipped
``RelayPolicy`` and Stage-3 ``CurrencyPolicy``. Pooled per stratum over its 3 contaminated conditions
(the four ``g00`` rows are the control and are never folded in).** Rates are misplaced fragments per
1,000 exon fragments; every per-condition row is printed above the pooled one.

============================  ==============  ==========  =============================  ==========
stratum (relay)               exon fragments  misplaced   the SAFE class                 the UNSAFE
============================  ==============  ==========  =============================  ==========
stranded x capture OFF        5,708,389       80,722      single-strand 4,657,494 f       AMBIG
  [IN SCOPE]                  (14.1/1k)                   at **8.9/1k**                   1,050,895 f
                                                                                          at 37.5/1k
stranded x capture ON         14,215,613      363,200     single-strand 11,362,957 f      AMBIG
  [IN SCOPE]                  (25.5/1k)                   at **16.1/1k**                  2,852,656 f
                                                                                          at 63.2/1k
unstranded x capture OFF      5,702,758       115,051     single-strand 4,653,017 f       AMBIG
  [IN SCOPE]                  (20.2/1k)                   at **15.8/1k**                  1,049,741 f
                                                                                          at 39.4/1k
unstranded x capture ON       14,211,148      3,333,716   ⛔ NOTHING SEPARATES — every    —
  [DEFERRED]                  (234.6/1k)                  bucket 161-425/1k
============================  ==============  ==========  =============================  ==========

⭐⭐ **THE ANSWER IS THE STRAND BIT, AND IT IS THE SAME MASK AS "HAS OWN COMPOSITION EVIDENCE" ON EVERY
STRANDED CONTAMINATED CONDITION** — identical to the slot, which this file reports rather than assumes.
The mechanism is in ``region_init``: the strand lambda-term is credited to SINGLE-STRAND slots only (at
an AMBIG slot the Schur complement zeroes it, because the tilt is a free nuisance), and at an exon it is
the only live source. So **AMBIG exons are exactly the evidence-free exons**, they carry 18-20 % of exon
mass and **~50 % of all exon misplacement**, and they are the population a landscape prior must be
trained WITHOUT and then applied TO.

⭐⭐ **THE SECOND PREDICTOR IS THE EXON'S OWN DEPTH, AND IT IS THE ONLY ONE THAT RANKS MONOTONICALLY ON
ALL THREE IN-SCOPE STRATA** (relay, deepest to shallowest quartile, equal mass in each):
**4.6 / 7.7 / 14.5 / 29.7** per 1k (stranded x capture-OFF), **8.6 / 18.4 / 25.2 / 50.2** (stranded x
capture-ON), **6.9 / 18.2 / 18.6 / 37.1** (unstranded x capture-OFF). ⭐ It is also nearly independent of
the strand bit, so "deep AND single-strand" is the training set and it is knowable from the payload.

⛔⛔ **THE FLANK STRUCTURE IS NOT A SAFE CRITERION — IT REVERSES BETWEEN STRATA.** ``sj only`` is the
BEST class at the zero control (3.9/1k vs 7.0 for ``sj & term``) and the WORST on stranded x capture-ON
(49.6/1k vs 17.0 for ``term only``), on the same index. Its spread is 1.5-2.9x where the strand bit's is
3.9-4.2x, so it neither ranks nor separates.

⛔ **HOPS-TO-A-PURE-gDNA-SLOT RANKS THE RIGHT WAY BUT WEAKLY** (5.2 -> 17.3/1k on stranded x capture-OFF,
9.7 -> 31.1 on stranded x capture-ON) and it INVERTS on the deferred stratum, where depth 1 is the worst
bucket of all (424.6/1k). ⭐ Its one sharp finding is a class no other predictor names: the **47 ERCC
references that carry a live exon have NO pure-gDNA slot anywhere in them** — the exon's flanks are
reference terminals — so their 39,354 fragments are ``unreachable``/``reference end`` and cannot be
reached by any message. They are the worst class at the zero control (17.2/1k) and on unstranded x
capture-OFF (29.4/1k), and among the best under capture — a small, structurally distinct population that
a training rule should exclude by construction rather than by score.

⛔⛔ **THE CURRENCY POLICY IS A REGRESSION AT EXONS ON ALL THREE IN-SCOPE STRATA AND A WIN ONLY ON THE
DEFERRED ONE** (misplaced fragments, currency vs relay): stranded x capture-OFF **99,200 vs 80,722**
(1.23x), stranded x capture-ON **803,343 vs 363,200** (2.21x), unstranded x capture-OFF **124,196 vs
115,051** (1.08x), unstranded x capture-ON **2,667,047 vs 3,333,716** (0.80x), zero control **99,537 vs
84,747** (1.17x). ⭐ **What it does is MOVE accuracy from the evidence-bearing exons to the evidence-free
ones**: on stranded x capture-OFF the AMBIG class goes 39,401 -> 24,532 misplaced while single-strand
goes 41,321 -> 74,668. That is the recipient-decides architecture doing exactly what it says, priced.

⛔⛔ **AND THE TEST CHROMOSOME SAYS THE OPPOSITE, WHICH IS `TRAPS: panel-before-src` MET AGAIN.** On the
8-condition test chromosome at ``g50`` stranded x capture-OFF the currency arm reads **413** misplaced
against relay's **1,662** — a 4x WIN — where the ladder pools to a 1.23x LOSS. ⛔ Do not read an exon
verdict off the test chromosome; use it to iterate and the ladder to decide.

⚠ **THE ZERO-gDNA CONTROL IS ONE-SIDED AND IS STAMPED AS SUCH** — truth is exactly 0 at every exon, so
any RNA-favouring answer reads as correct and a class can look safe for a reason that does not survive
contamination.

⚠ **TWO SURPRISES WORTH A SECOND LOOK, REPORTED RATHER THAN EXPLAINED.** ① On an UNSTRANDED library the
own-evidence mask is empty at ``g05``/``g50`` — ``kappa = 1/2`` makes the strand discriminability
identically 0 — but at ``g98`` it opens on 305,606 exon fragments, and those slots are **7.6x WORSE**
than the evidence-free ones (113.2/1k vs 14.9/1k). An unstranded slot cannot have strand composition
evidence, so that is the deadband letting a phantom through where RNA is scarcest. ② At the zero control
the mask is empty on the STRANDED conditions too, which is the ``1/N_gdna`` term in the same deadband
doing its job (no gDNA ⇒ no discriminability) — correct, and worth knowing before reading a ``g00`` row.

⭐ Cost: ~8 s/condition on the test chromosome for three arms, ~2-4 min/condition on the ladder for two.

Usage::

    python scripts/design/exon_solvability.py --suite ... --index ...     # every cached condition
    python scripts/design/exon_solvability.py --condition NAME --arms relay currency
    python scripts/design/exon_solvability.py --refits 0                  # the PASS-0 reading
    python scripts/design/exon_solvability.py --self-test                 # no I/O, perturbed, 39/39
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


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


OC = _sibling("object_composition.py")
HC = _sibling("hop_currency.py")
PVO = OC.PVO

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import build_region_statics  # noqa: E402
from rigel.calibration.region_init import has_own_composition_evidence  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    FLAG_DONOR_POS,
    FLAG_TSS_POS,
    build_boundary_flags_array,
    is_splice_site,
    is_terminus,
)
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from rigel.types import Strand  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX

#: the EXON REGION population, named by `object_composition.strata` and carried in the certified
#: `slot_truth.npz`. ⛔ One name, one home — this file never re-derives the label.
EXON = "R exon"

#: ⭐ **THE CURVE'S X AXIS — A PRESENTATION GRID, NOT A THRESHOLD.** Nothing in this file branches on
#: it: every verdict is the whole curve, the mass-weighted quantiles (which need no grid) and
#: ``Σ err·M`` in fragments. It is here so two conditions print on the same axis and can be read side
#: by side. ``|f_g − true_f_g|`` is bounded by 1 because both are compositions, so the last column is
#: the bucket's whole mass BY CONSTRUCTION and a shortfall there is a bug, not a finding.
_CURVE_X = (0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 1.00)

#: mass-weighted quantile levels reported beside the curve. Also presentation: the quantile function
#: is threshold-free, and these are the points of it that get printed.
_QUANTILES = (0.50, 0.90)

#: ⭐ **THE ARMS.** ``key`` is the diagnostic key that PROVES the policy ran — the relay publishes
#: ``_uni`` and the currency policy ``_currency``; ``None`` means the arm must be message-free, which
#: is checked by ``f_g == fg_loc`` bit for bit (`TRAPS: an-ablation-that-never-ran`).
ARMS: dict[str, tuple[dict, str | None]] = {
    "relay": ({"message_propagation": True, "message_policy": "relay"}, "_uni"),
    "currency": ({"message_propagation": True, "message_policy": "currency"}, "_currency"),
    "silent": ({"message_propagation": False}, None),
}


# ── ① the threshold-free scorers ─────────────────────────────────────────────────────────────────


def mass_curve(err, mass, xs=_CURVE_X) -> np.ndarray:
    """FRAGMENTS whose ``|f_g − true_f_g|`` is at most each ``x``. Non-decreasing in ``x`` by
    construction; the last entry is the whole mass when the error axis is a composition."""
    err = np.asarray(err, np.float64)
    mass = np.asarray(mass, np.float64)
    return np.array([float(mass[err <= float(x)].sum()) for x in xs], np.float64)


def mass_quantiles(err, mass, qs=_QUANTILES) -> np.ndarray:
    """The mass-weighted quantiles of the error — the smallest error at which the cumulative FRAGMENT
    mass reaches ``q`` of the total. ⛔ No interpolation: an interpolated quantile invents an error
    value no slot has, and this axis is read as "what does a fragment at this rank actually see"."""
    err = np.asarray(err, np.float64)
    mass = np.asarray(mass, np.float64)
    total = float(mass.sum())
    if total <= 0.0:
        return np.full(len(qs), np.nan)
    order = np.argsort(err, kind="stable")
    cum = np.cumsum(mass[order])
    out = []
    for q in qs:
        j = int(np.searchsorted(cum, float(q) * total, side="left"))
        out.append(float(err[order][min(j, err.shape[0] - 1)]))
    return np.array(out, np.float64)


def misplaced_fragments(err, mass) -> float:
    """``Σ|f_g − true_f_g|·M`` — how many fragments the calibration puts on the wrong side of the
    gDNA/RNA split at these slots. The project's standard scalar, in FRAGMENTS, never a ratio."""
    return float((np.asarray(err, np.float64) * np.asarray(mass, np.float64)).sum())


# ── ② the truth-free predictors ──────────────────────────────────────────────────────────────────


def depth_to_pure_gdna(left, right, ref, free_pos, free_neg) -> np.ndarray:
    """Hops from each slot to the nearest STRUCTURALLY PURE-gDNA slot in the same reference — the
    slots where neither RNA strand is admissible, so gDNA is measured rather than imputed. ``0`` on
    such a slot, ``inf`` when the whole reference admits RNA everywhere.

    ⭐ The traversal is ``hop_currency.depth_to_measured`` — one home, so a chain-adjacency fix cannot
    reach one instrument and miss the other.
    """
    locked = ~np.asarray(free_pos, bool) & ~np.asarray(free_neg, bool)
    return HC.depth_to_measured(left, right, ref, ~locked)


def depth_bucket(depth) -> np.ndarray:
    """Label each slot by hops-to-pure-gDNA, on ``hop_currency``'s own bucket edges (exact to 4, then
    5-8, then 9+) — the same axis the hop census prints, so the two tables can be laid side by side."""
    d = np.asarray(depth, np.float64)
    out = np.full(d.shape[0], "unreachable", dtype=object)
    lo = 0
    for hi in (0, *HC._DEPTH_EDGES):
        sel = (d >= lo) & (d <= hi)
        out[sel] = f"{lo}" if lo == hi else f"{lo}-{hi}"
        lo = hi + 1
    out[(d >= lo) & np.isfinite(d)] = f">={lo}"
    return out.astype(str)


def flank_class(kind, left, right, flags) -> np.ndarray:
    """Per slot, what its two chain flanks ARE, from the splice graph's own flags on either strand:
    ``"sj only"``, ``"term only"``, ``"sj & term"``, ``"reference end"`` when a flank is a terminal.

    ⭐ **This is the split the seven object classes cannot see** — the same finding
    ``hop_currency.boundary_structure`` records at hop level, asked of a REGION about its
    neighbourhood: an exon reached across a SPLICE SITE receives a population that crossed into it,
    an exon reached across a TSS/TES does not.
    """
    kind = np.asarray(kind)
    flags = np.asarray(flags, np.uint16)
    term = is_terminus(flags, Strand.POS) | is_terminus(flags, Strand.NEG)
    spl = is_splice_site(flags, Strand.POS) | is_splice_site(flags, Strand.NEG)
    n = kind.shape[0]
    left = np.asarray(left, np.int64)
    right = np.asarray(right, np.int64)

    def gather(nbr, arr):
        ok = nbr >= 0
        return ok & arr[np.clip(nbr, 0, max(n - 1, 0))]

    has_sj = gather(left, spl) | gather(right, spl)
    has_term = gather(left, term) | gather(right, term)
    ends = (left < 0) | (right < 0)
    out = np.full(n, "reference end", dtype=object)
    out[has_sj & ~has_term] = "sj only"
    out[~has_sj & has_term] = "term only"
    out[has_sj & has_term] = "sj & term"
    out[ends & ~has_sj & ~has_term] = "reference end"
    return out.astype(str)


def strand_class(free_pos, free_neg) -> np.ndarray:
    """``"AMBIG"`` where both RNA strands are admissible, ``"single-strand"`` where one is,
    ``"no RNA"`` where neither (AXIOM 0's two bits, read as a label)."""
    fp = np.asarray(free_pos, bool)
    fn = np.asarray(free_neg, bool)
    out = np.full(fp.shape[0], "no RNA", dtype=object)
    out[fp ^ fn] = "single-strand"
    out[fp & fn] = "AMBIG"
    return out.astype(str)


def quartile_bucket(value, mass) -> tuple[np.ndarray, list[float]]:
    """Label each slot by which mass-weighted QUARTILE of ``value`` it falls in — edges derived from
    the data, so nothing is chosen. ⚠ Ties are kept together (the edges are VALUES, not ranks), so a
    bucket's mass is only approximately a quarter and the realised mass is printed rather than assumed.
    """
    v = np.asarray(value, np.float64)
    m = np.asarray(mass, np.float64)
    total = float(m.sum())
    if total <= 0.0 or v.size == 0:
        return np.full(v.shape[0], "quartile 1", dtype="<U12"), []
    order = np.argsort(v, kind="stable")
    cum = np.cumsum(m[order])
    edges = [float(v[order][min(int(np.searchsorted(cum, q * total, side="left")), v.size - 1)])
             for q in (0.25, 0.50, 0.75)]
    idx = np.searchsorted(np.asarray(edges, np.float64), v, side="right")
    labels = np.array([f"quartile {i + 1}" for i in range(4)])
    return labels[np.clip(idx, 0, 3)], edges


# ── ③ the tables ─────────────────────────────────────────────────────────────────────────────────


def bucket_rows(labels, err, mass, order=None) -> list[dict]:
    """One row per observable bucket: slots, FRAGMENTS, misplaced FRAGMENTS, the curve in FRAGMENTS
    and the mass-weighted quantiles. ⛔ The rows must PARTITION the population — a slot counted twice
    or not at all silently reweights the answer, so it is asserted rather than assumed."""
    labels = np.asarray(labels).astype(str)
    err = np.asarray(err, np.float64)
    mass = np.asarray(mass, np.float64)
    names = list(order) if order else sorted(set(labels.tolist()))
    rows = []
    seen = 0
    for nm in names:
        sel = labels == nm
        seen += int(sel.sum())
        if not sel.any():
            continue
        rows.append({
            "bucket": nm,
            "slots": int(sel.sum()),
            "mass": float(mass[sel].sum()),
            "misplaced": misplaced_fragments(err[sel], mass[sel]),
            "curve": mass_curve(err[sel], mass[sel]),
            "quantiles": mass_quantiles(err[sel], mass[sel]),
        })
    if seen != labels.shape[0]:
        missing = sorted(set(labels.tolist()) - set(names))
        raise AssertionError(
            f"the buckets do not partition the population: {seen:,} of {labels.shape[0]:,} slots "
            f"labelled, unhandled label(s) {missing}. A mass-weighted table over a non-partition is a "
            "mean over the wrong denominator, and an UNEXPECTED bucket is exactly the thing this "
            "instrument would otherwise report as a smaller population."
        )
    return rows


def run_arm(payload, kw, arm: str, refits: int | None = None,
            region_bank: str | None = None) -> dict:
    """One ``calibrate`` under one policy, with the ran/did-not-run assertion the repo requires: an
    arm that could not express its own policy scores byte-identically to another and looks like a
    finding (`TRAPS: an-ablation-that-never-ran`).

    ⭐ ``refits=0`` is the PASS-0 reading, and it is the one the training question actually wants: the
    landscape refits consume a prior fitted on the solved objects, so scoring the shipped belief and
    then training a prior on what it got right is circular. The default is the shipped setting, so both
    readings are available and neither is assumed.
    """
    overrides, key = ARMS[arm]
    cfg = dataclasses.replace(CalibrationConfig(), **overrides)
    if region_bank is not None:
        cfg = dataclasses.replace(cfg, region_abundance_bank=region_bank)
    if refits is not None:
        cfg = dataclasses.replace(cfg, calib_refit_iters=int(refits))
    debug: dict = {}
    calibrate(payload=payload, config=cfg, _debug=debug,
              **{k: v for k, v in kw.items() if k != "payload"})
    cap = debug["capture"]
    assert_arm_ran(cap, arm)
    return cap


def assert_arm_ran(cap: dict, arm: str) -> None:
    """The arm's own falsification, split out so the self-test can perturb it with no I/O."""
    key = ARMS[arm][1]
    if key is None:
        if not np.array_equal(np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])):
            raise AssertionError(f"arm {arm!r} is meant to be message-free, but its final belief "
                                 "differs from its own local solve — it is not muted")
        if any(k in cap for k in ("_uni", "_currency")):
            raise AssertionError(f"arm {arm!r} is meant to be message-free but a policy published state")
        return
    if key not in cap:
        raise AssertionError(
            f"arm {arm!r} did not run: the capture carries no {key!r}. An arm the config cannot "
            "express scores identically to another and reads as a finding.")


def load_condition(index, region_arrays, bflags, suite: Path, condition: str, arms,
                   refits: int | None = None, region_bank: str | None = None) -> dict:
    """Everything one condition needs. REFUSES a condition without a certified ``slot_truth.npz`` —
    this instrument derives no truth of its own, on purpose."""
    truth_path = Path(suite) / "oracle_cache" / condition / "slot_truth.npz"
    if not truth_path.is_file():
        raise FileNotFoundError(
            f"{truth_path} is missing — run `calibration_oracle.py --condition {condition}` first; "
            "this instrument scores against CERTIFIED truth only.")
    truth = np.load(truth_path, allow_pickle=True)
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    chain = build_region_chain(cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, bflags)

    kind = np.asarray(chain.kind)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    ref = np.asarray(truth["ref"])
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    stratum = np.asarray(truth["stratum"]).astype(str)
    mass = np.asarray(truth["count"], np.float64)
    exon = np.asarray(truth["live"], bool) & (stratum == EXON)
    if not np.all(kind[exon] == REGION):
        raise AssertionError(f"a {EXON!r} slot is not a REGION — the certified table and this chain "
                             "describe different scans")

    caps = {a: run_arm(cache.payload, kw, a, refits, region_bank=region_bank) for a in arms}
    any_cap = caps[arms[0]]
    obs = np.asarray(any_cap["count"], np.float64).sum(axis=1)
    n_bad = int(np.sum(obs != mass))
    gates = [{"gate": "gate truth-agrees", "ok": n_bad == 0, "n_bad": n_bad}]
    same = np.array_equal(np.asarray(any_cap["free_pos"], bool), fp) and np.array_equal(
        np.asarray(any_cap["free_neg"], bool), fn)
    gates.append({"gate": "gate statics-agree", "ok": same, "n_bad": 0 if same else 1})
    depth = depth_to_pure_gdna(left, right, ref, fp, fn)
    if exon.any() and float(np.nanmin(depth[exon])) < 1.0:
        raise AssertionError("an exon REGION reports depth 0 to a pure-gDNA slot — an exon admits at "
                             "least one RNA strand by construction, so it can never BE one")
    ev = {a: has_own_composition_evidence(np.asarray(caps[a]["_tau0_lam"], np.float64)) for a in arms}
    return {
        "condition": condition, "exon": exon, "mass": mass,
        "true_f_g": np.asarray(truth["true_f_g"], np.float64),
        "field_certified": bool(truth["field_certified"]),
        "f_g": {a: np.asarray(caps[a]["f_g"], np.float64) for a in arms},
        "fg_loc": {a: np.asarray(caps[a]["fg_loc"], np.float64) for a in arms},
        "own_evidence": ev, "depth": depth, "obs": obs,
        "flank": flank_class(kind, left, right, np.asarray(statics.boundary_flags)),
        "strand": strand_class(fp, fn),
        "gates": gates, "arms": list(arms),
        "refits": int(CalibrationConfig().calib_refit_iters) if refits is None else int(refits),
    }


def measure(d: dict) -> dict:
    """Every predictor's table, per arm. Pure — no I/O, so the self-test feeds it synthetic slots."""
    sel = np.asarray(d["exon"], bool)
    mass = np.asarray(d["mass"], np.float64)[sel]
    true_fg = np.asarray(d["true_f_g"], np.float64)[sel]
    dep = depth_bucket(np.asarray(d["depth"], np.float64))[sel]
    flank = np.asarray(d["flank"]).astype(str)[sel]
    strand = np.asarray(d["strand"]).astype(str)[sel]
    obs_q, obs_edges = quartile_bucket(np.asarray(d["obs"], np.float64)[sel], mass)
    out = {"mass": float(mass.sum()), "slots": int(sel.sum()), "obs_edges": obs_edges, "arms": {}}
    for a in d["arms"]:
        err = np.abs(np.asarray(d["f_g"][a], np.float64)[sel] - true_fg)
        ev = np.where(np.asarray(d["own_evidence"][a], bool)[sel], "own evidence", "NO own evidence")
        joint = np.char.add(np.char.add(flank.astype(str), "  x  "), strand.astype(str))
        # ⭐ **THE TWO PREDICTORS ARE OFTEN THE SAME MASK, AND SAYING SO IS THE FINDING.** The strand
        #   arm's lambda-term is credited to SINGLE-STRAND slots only (the Schur complement zeroes it at
        #   an AMBIG one), so at an exon "has own composition evidence" collapses onto "is single-strand"
        #   whenever the strand channel is the only live source. Reported, never asserted — an intron
        #   density factor can also feed tau, and then the two separate.
        collinear = bool(np.array_equal(ev == "own evidence", strand == "single-strand"))
        out["arms"][a] = {
            "own_evidence_is_single_strand": collinear,
            "all": bucket_rows(np.full(err.shape[0], "ALL exon REGIONs"), err, mass),
            "own evidence": bucket_rows(ev, err, mass, order=["own evidence", "NO own evidence"]),
            "hops to a pure-gDNA slot": bucket_rows(dep, err, mass),
            "flank structure": bucket_rows(
                flank, err, mass, order=["sj only", "sj & term", "term only", "reference end"]),
            "RNA strands admissible": bucket_rows(
                strand, err, mass, order=["single-strand", "AMBIG", "no RNA"]),
            "own observed depth": bucket_rows(
                obs_q, err, mass, order=[f"quartile {i + 1}" for i in range(4)]),
            "flank x strand": bucket_rows(joint, err, mass),
        }
    return out


# ── reporting ────────────────────────────────────────────────────────────────────────────────────

PREDICTORS = ("own evidence", "hops to a pure-gDNA slot", "flank structure",
              "RNA strands admissible", "own observed depth", "flank x strand")


def _hdr() -> str:
    cols = "".join(f"{'<=' + f'{x:g}':>9}" for x in _CURVE_X)
    qs = "".join(f"{'p' + f'{100 * q:.0f}':>7}" for q in _QUANTILES)
    return (f"   {'bucket':<26} {'slots':>7} {'fragments':>13} {'misplaced':>12}  |"
            f"{cols}  |{qs}")


def _row(r: dict) -> str:
    m = max(r["mass"], 1.0)
    cols = "".join(f"{100 * c / m:8.1f}%" for c in r["curve"])
    qs = "".join(f"{q:7.3f}" for q in r["quantiles"])
    return (f"   {r['bucket']:<26} {r['slots']:>7,} {r['mass']:>13,.0f} {r['misplaced']:>12,.0f}  |"
            f"{cols}  |{qs}")


def report_condition(d: dict, m: dict, verbose: bool) -> None:
    cond = d["condition"]
    scope = OC._scope(cond)
    ss, cap = PVO.stratum(cond)
    print(f"\n== {cond}   [{ss} x {cap} — {scope}]   "
          f"truth {'COMPOSITION+FIELD' if d['field_certified'] else 'COMPOSITION only'}   "
          f"calib_refit_iters={d['refits']}")
    for g in d["gates"]:
        print(f"   {'OK' if g['ok'] else 'FAIL'}  {g['gate']}" + ("" if g["ok"] else f"  {g}"))
    print(f"   exon REGIONs: {m['slots']:,} live slots, {m['mass']:,.0f} fragments"
          + (f"   own-depth quartile edges (fragments/slot): "
             f"{', '.join(f'{e:,.0f}' for e in m['obs_edges'])}" if m["obs_edges"] else ""))
    if PVO.is_zero_gdna(cond):
        print("   ⛔ ZERO-gDNA CONTROL: truth is exactly 0 everywhere, so any RNA-favouring answer "
              "reads as correct here. One-sided — never rank a class on this row.")
    for a in d["arms"]:
        print(f"\n   ARM {a}   — fraction of each bucket's exon MASS with |f_g − true_f_g| <= x; "
              "'misplaced' is Sum|f_g − true|*M in FRAGMENTS")
        if m["arms"][a]["own_evidence_is_single_strand"]:
            print("   note: 'own evidence' and 'single-strand' are the SAME mask here — the strand "
                  "lambda-term is the only live source, and it is credited to single-strand slots only")
        print(_hdr())
        print(_row(m["arms"][a]["all"][0]))
        for p in PREDICTORS:
            if p == "flank x strand" and not verbose:
                continue
            rows = m["arms"][a][p]
            if len(rows) < 2 and not verbose:
                print(f"   -- {p}: DEGENERATE, one bucket ({rows[0]['bucket']}) holds all the mass "
                      "— this predictor separates nothing here")
                continue
            print(f"   -- {p}")
            for r in sorted(rows, key=lambda r: -r["mass"]):
                print(_row(r))


def report_summary(results: list[tuple[dict, dict]], arms) -> None:
    """⛔ THE POOLED TABLE, AND IT COMES LAST BECAUSE IT IS A SUMMARY OF THE ROWS ABOVE IT (owner's
    rule). Pooled WITHIN a stratum only — the panel total hides a sign flip between strata."""
    print("\n\n" + "=" * 118)
    print("SUMMARY — pooled per STRATUM, and this is a summary of the per-condition rows above it, "
          "not a substitute for them.")
    print("⛔ The zero-gDNA control is its own row and is never folded into a stratum.")
    groups: dict[str, list[tuple[dict, dict]]] = {}
    for d, m in results:
        c = d["condition"]
        key = ("g00 ZERO-gDNA CONTROL" if PVO.is_zero_gdna(c)
               else f"{' x '.join(PVO.stratum(c))}  [{OC._scope(c)}]")
        groups.setdefault(key, []).append((d, m))
    for a in arms:
        print(f"\n   ARM {a}")
        print(f"   {'stratum':<38} {'conds':>5} {'predictor':<26} {'bucket':<16} "
              f"{'fragments':>13} {'misplaced':>12} {'misplaced/1k frag':>18}")
        for key, rows in sorted(groups.items()):
            for p in ("flank structure", "hops to a pure-gDNA slot"):
                acc: dict[str, list[float]] = {}
                for _d, m in rows:
                    for r in m["arms"][a][p]:
                        cur = acc.setdefault(r["bucket"], [0.0, 0.0])
                        cur[0] += r["mass"]
                        cur[1] += r["misplaced"]
                for b, (mass, mis) in sorted(acc.items(), key=lambda kv: -kv[1][0]):
                    rate = 1000.0 * mis / max(mass, 1.0)
                    print(f"   {key:<38} {len(rows):>5} {p:<26} {b:<16} {mass:>13,.0f} "
                          f"{mis:>12,.0f} {rate:>18,.1f}")
            print()


def write_tsv(path: Path, results: list[tuple[dict, dict]], arms) -> None:
    cols = ["condition", "strand", "capture", "scope", "arm", "predictor", "bucket", "slots",
            "fragments", "misplaced", *[f"frag_le_{x:g}" for x in _CURVE_X],
            *[f"p{100 * q:.0f}" for q in _QUANTILES]]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for d, m in results:
            c = d["condition"]
            ss, cap = PVO.stratum(c)
            for a in arms:
                for p in ("all", *PREDICTORS):
                    for r in m["arms"][a][p]:
                        fh.write("\t".join(str(v) for v in [
                            c, ss, cap, OC._scope(c), a, p, r["bucket"], r["slots"],
                            f"{r['mass']:.0f}", f"{r['misplaced']:.2f}",
                            *[f"{v:.0f}" for v in r["curve"]],
                            *[f"{v:.6f}" for v in r["quantiles"]]]) + "\n")
    print(f"\n   wrote {path}")


# ── the self-test: no I/O, every check falsified by perturbation ─────────────────────────────────


def self_test() -> int:
    passed = failed = 0

    def check(name, cond):
        nonlocal passed, failed
        ok = bool(cond)
        passed, failed = passed + ok, failed + (not ok)
        print(f"   {'PASS' if ok else 'FAIL'}  {name}")

    # ① the curve
    err = np.array([0.0, 0.02, 0.3])
    mass = np.array([10.0, 20.0, 30.0])
    c = mass_curve(err, mass, (0.001, 0.05, 1.0))
    check("curve counts FRAGMENTS at or below x", np.allclose(c, [10.0, 30.0, 60.0]))
    check("curve — perturbed mass moves it",
          not np.allclose(mass_curve(err, mass + np.array([1.0, 0, 0]), (0.001, 0.05, 1.0)), c))
    check("curve is non-decreasing in x", bool(np.all(np.diff(mass_curve(err, mass)) >= 0)))
    check("curve monotonicity — a reversed curve is caught",
          not bool(np.all(np.diff(mass_curve(err, mass)[::-1]) >= 0)))
    check("curve closes at x=1 (a composition error cannot exceed 1)",
          np.isclose(mass_curve(err, mass)[-1], mass.sum()))
    check("curve closure — an out-of-range error is caught",
          not np.isclose(mass_curve(np.array([0.0, 0.02, 1.5]), mass)[-1], mass.sum()))

    # ② the quantiles
    q = mass_quantiles(np.array([0.0, 1.0]), np.array([3.0, 1.0]), (0.5, 0.9))
    check("mass-weighted quantile follows the MASS, not the slot count", np.allclose(q, [0.0, 1.0]))
    check("quantile — perturbing the mass split moves it",
          not np.allclose(mass_quantiles(np.array([0.0, 1.0]), np.array([1.0, 3.0]), (0.5, 0.9)), q))
    check("quantile is non-decreasing in q",
          bool(np.all(np.diff(mass_quantiles(err, mass, (0.1, 0.5, 0.9, 1.0))) >= 0)))
    check("empty population gives nan, never 0", bool(np.all(np.isnan(mass_quantiles([], [])))))

    # ③ misplaced fragments
    check("misplaced is Sum err*M in fragments",
          np.isclose(misplaced_fragments([0.5, 0.25], [100.0, 200.0]), 100.0))
    check("misplaced — doubling the error doubles it",
          np.isclose(misplaced_fragments([1.0, 0.5], [100.0, 200.0]), 200.0))

    # ④ hops to a pure-gDNA slot, on a hand chain  R B R B R B R
    n = 7
    left = np.array([-1, 0, 1, 2, 3, 4, 5], np.int64)
    right = np.array([1, 2, 3, 4, 5, 6, -1], np.int64)
    ref = np.zeros(n, np.int32)
    fp = np.array([False, False, True, True, True, True, True])
    fn = np.zeros(n, bool)
    dep = depth_to_pure_gdna(left, right, ref, fp, fn)
    check("depth is 0 at a pure-gDNA slot and counts hops away from it",
          np.allclose(dep, [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]))
    check("depth — unlocking every slot makes it unreachable",
          bool(np.all(~np.isfinite(depth_to_pure_gdna(left, right, ref, np.ones(n, bool), fn)))))
    two = ref.copy()
    two[3:] = 1
    check("depth does not cross a reference boundary",
          bool(np.all(~np.isfinite(depth_to_pure_gdna(left, right, two, fp, fn)[3:]))))
    check("depth buckets label the hand chain",
          list(depth_bucket(np.array([0.0, 1.0, 4.0, 6.0, 12.0, np.inf]))) ==
          ["0", "1", "4", "5-8", ">=9", "unreachable"])
    check("depth buckets — a perturbed depth changes its bucket",
          depth_bucket(np.array([9.0]))[0] != depth_bucket(np.array([4.0]))[0])

    # ⑤ the flank class, on a hand chain  R B R B R  with a sj at slot 1 and a terminus at slot 3
    kind = np.array([REGION, BOUNDARY, REGION, BOUNDARY, REGION], np.int8)
    lf = np.array([-1, 0, 1, 2, 3], np.int64)
    rt = np.array([1, 2, 3, 4, -1], np.int64)
    sj_flag, tss_flag = np.uint16(int(FLAG_DONOR_POS)), np.uint16(int(FLAG_TSS_POS))
    flags = np.array([0, sj_flag, 0, tss_flag, 0], np.uint16)
    fc = flank_class(kind, lf, rt, flags)
    check("flank class reads the splice graph's own flags",
          list(fc[[0, 2, 4]]) == ["sj only", "sj & term", "term only"])
    check("flank class — clearing the sj flag reclassifies the exon between them",
          flank_class(kind, lf, rt, np.array([0, 0, 0, tss_flag, 0], np.uint16))[2] == "term only")
    check("flank class — a slot with no flagged flank is a reference end",
          flank_class(kind, lf, rt, np.zeros(5, np.uint16))[0] == "reference end")

    # ⑥ strand class
    sc = strand_class([True, True, False, False], [True, False, True, False])
    check("strand class names AXIOM 0's two bits",
          list(sc) == ["AMBIG", "single-strand", "single-strand", "no RNA"])
    check("strand class — flipping a bit moves the label",
          strand_class([True], [False])[0] != strand_class([True], [True])[0])

    # ⑦ the quartile bucket
    val = np.arange(8, dtype=np.float64)
    lab, edges = quartile_bucket(val, np.ones(8))
    check("quartile buckets are monotone in the value and cover all four",
          sorted(set(lab.tolist())) == [f"quartile {i + 1}" for i in range(4)]
          and list(lab) == sorted(lab.tolist()))
    check("quartile buckets — a mass-weighted edge follows the MASS",
          quartile_bucket(val, np.array([100.0, 1, 1, 1, 1, 1, 1, 1]))[1] != edges)

    # ⑧ the tables partition
    lab2 = np.array(["a", "b", "a", "b"])
    rows = bucket_rows(lab2, [0.1, 0.2, 0.3, 0.4], [1.0, 2.0, 3.0, 4.0])
    check("bucket rows partition the mass",
          np.isclose(sum(r["mass"] for r in rows), 10.0) and len(rows) == 2)
    try:
        bucket_rows(np.array(["a", "b"]), [0.1, 0.2], [1.0, 2.0], order=["a"])
        ok = False
    except AssertionError:
        ok = True
    check("bucket rows — a label the declared order does not cover is REFUSED", ok)
    check("bucket row totals equal the ungrouped total",
          np.isclose(sum(r["misplaced"] for r in rows),
                     misplaced_fragments([0.1, 0.2, 0.3, 0.4], [1.0, 2.0, 3.0, 4.0])))

    # ⑨ the arm assertions
    good = {"f_g": np.array([0.5]), "fg_loc": np.array([0.5])}
    try:
        assert_arm_ran(good, "silent")
        ok = True
    except AssertionError:
        ok = False
    check("a muted arm whose belief equals its local solve passes", ok)
    for nm, capd, arm in (("a muted arm that MOVED is refused",
                           {"f_g": np.array([0.4]), "fg_loc": np.array([0.5])}, "silent"),
                          ("a muted arm that published policy state is refused",
                           {"f_g": np.array([0.5]), "fg_loc": np.array([0.5]), "_uni": 1}, "silent"),
                          ("a relay arm that published nothing is refused", good, "relay"),
                          ("a currency arm that published nothing is refused", good, "currency")):
        try:
            assert_arm_ran(capd, arm)
            ok = False
        except AssertionError:
            ok = True
        check(nm, ok)
    check("the relay and currency arms are proven by DIFFERENT keys",
          ARMS["relay"][1] != ARMS["currency"][1])

    # ⑩ measure() end to end on synthetic slots — the assembled dict, no I/O
    d = {
        "condition": "gdna_g50_ss_0.99_nrna_mid_capture_off", "arms": ["relay"],
        "exon": np.array([True, True, True, False]),
        "mass": np.array([100.0, 200.0, 300.0, 400.0]),
        "true_f_g": np.array([0.5, 0.5, 0.5, 1.0]),
        "f_g": {"relay": np.array([0.5, 0.6, 0.9, 1.0])},
        "fg_loc": {"relay": np.array([0.5, 0.6, 0.9, 1.0])},
        "own_evidence": {"relay": np.array([True, True, False, True])},
        "depth": np.array([1.0, 2.0, 9.0, 0.0]),
        "obs": np.array([10.0, 20.0, 30.0, 40.0]),
        "flank": np.array(["sj only", "sj & term", "term only", "sj only"]),
        "strand": np.array(["AMBIG", "single-strand", "AMBIG", "no RNA"]),
        "field_certified": True, "gates": [],
    }
    mm = measure(d)
    check("measure() scores only the exon REGIONs",
          np.isclose(mm["mass"], 600.0) and mm["slots"] == 3)
    check("measure() reproduces the ungrouped misplacement",
          np.isclose(mm["arms"]["relay"]["all"][0]["misplaced"], 0.0 + 200 * 0.1 + 300 * 0.4))
    for p in PREDICTORS:
        tot = sum(r["mass"] for r in mm["arms"]["relay"][p])
        if not np.isclose(tot, 600.0):
            check(f"predictor {p!r} partitions the exon mass", False)
            break
    else:
        check("every predictor partitions the exon mass", True)
    check("the own-evidence / single-strand collinearity flag is FALSE when the masks differ",
          not mm["arms"]["relay"]["own_evidence_is_single_strand"])
    d_col = dict(d, own_evidence={"relay": np.array([False, True, False, True])})
    check("the collinearity flag — it is TRUE when they are the same mask",
          measure(d_col)["arms"]["relay"]["own_evidence_is_single_strand"])
    d2 = dict(d, exon=np.ones(4, bool))
    check("measure() — widening the population changes the mass",
          not np.isclose(measure(d2)["mass"], mm["mass"]))

    print(f"\n   {passed} passed, {failed} failed")
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None, help="one condition (default: every scan_cache entry)")
    ap.add_argument("--arms", nargs="+", default=["relay", "currency"], choices=sorted(ARMS),
                    help="which message policies to run (default: the shipped relay and the Stage-3 "
                         "currency policy)")
    ap.add_argument("--refits", type=int, default=None,
                    help="calib_refit_iters (default: the shipped value). 0 is the PASS-0 reading — "
                         "the one the training question wants, since the refits consume a fitted prior")
    ap.add_argument("--out", type=Path, default=None, help="write every row as TSV, in FRAGMENTS")
    ap.add_argument("--verbose", action="store_true",
                    help="print the joint flank x strand table and degenerate predictors in full")
    ap.add_argument("--region-bank", default=None, choices=("contained", "start"),
                    help="override CalibrationConfig.region_abundance_bank for every arm (None = "
                         "the shipped default). Only the currency arm can read it.")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    bflags = build_boundary_flags_array(index)
    conds = ([args.condition] if args.condition
             else sorted(p.name for p in (args.suite / "scan_cache").iterdir() if p.is_dir()))
    results, bad = [], 0
    for c in conds:
        d = load_condition(index, region_arrays, bflags, args.suite, c, args.arms, args.refits,
                           region_bank=args.region_bank)
        m = measure(d)
        report_condition(d, m, args.verbose)
        bad += sum(not g["ok"] for g in d["gates"])
        results.append((d, m))
    if len(results) > 1:
        report_summary(results, args.arms)
    if args.out:
        write_tsv(args.out, results, args.arms)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
