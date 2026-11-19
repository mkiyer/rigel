#!/usr/bin/env python
"""Does one arm forecast each fragment's source better than another, scored on identical read names?

The per-locus EM gives every fragment a posterior over its sources. ``--record`` runs one arm, the shipped pipeline
under its ``--set`` overrides on the build this process imports, on each condition, and keeps for every read pair of
the input BAM its forecast as four class masses (annotated RNA, nascent RNA, the locus's gDNA, intergenic gDNA, in
``native.EM_CLASSES`` order) beside the source its read name gives. ``--pair`` joins two arms' records by read name and
scores both on the fragments both forecast, so a difference cannot come from scoring different fragments. A read pair
one arm forecasts and the other does not is counted by source and by the other arm's status, never dropped. Two
``--set`` configurations of one build are two ``--record`` runs; a second build records under its own interpreter.
``--arm no_capture`` records under ``quant_accuracy.py``'s arm of that name, every capture efficiency read as 1, the
line every capture-ON row is read beside: through its own ``install_no_capture``, refused by its ``check_no_capture``.
It substitutes no truth, so it records from an aligner's BAM too. The record carries its arm as ``qa_arm``, and a pair
names it in what varies.

Where a forecast comes from:

- an EM unit: its posterior by class, the posterior the assignment keeps or draws from
  (``AbundanceEstimator.run_batch_locus_em_partitioned(emit_class_mass=True)``), so it is the same in either
  assignment mode;
- a fragment routing decides alone: a deterministic one is mass 1 on its transcript's class (annotated, or nascent
  for a synthetic span), a multimapper the gates left only intergenic places is mass 1 intergenic;
- a placed fragment that reaches no locus (the intergenic fast path): mass 1 intergenic. Only the annotated-BAM pass
  names these (``ZF`` intergenic on a fragment the annotation table does not hold), so each run writes that BAM to
  ``--work-dir``, reads it once for the read names, and deletes it.

Every other read pair has no forecast and is counted by status: ``unassigned`` (the EM or a gate holds it, with no
mass), ``withheld`` (it has a fragment id, the scan placed it nowhere) or ``excluded`` (no eligible record, so no
fragment id). The read names are the input BAM's; a template the simulator drew and the BAM lacks is counted per
source as ``not_in_bam``.

The scores are sums over fragments, in fragment units, split by the read name's source:

- the source forecast is ``p = P(gDNA) = (locus gDNA + intergenic gDNA) / (the row's sum)``. Its Brier score
  ``Σ (p − y)²`` is proper; the linear score ``Σ |p − y|``, the fractional source error, is not, because it rewards
  certainty: of 100 indistinguishable fragments, 40 of them gDNA, the calibrated forecast (0.4 each) costs linear 48
  and Brier 24, the certain one (RNA, all of them) linear 40 and Brier 40;
- the calibration check is CORP (Dimitriadis, Gneiting and Jordan 2021, PNAS 118:e2016191118), within one
  condition: ties pooled, the pool-adjacent-violators fit of the truth on the forecast, each block's value its exact
  count ratio. MCB, the Brier the monotone recalibration removes, is the miscalibration; DSC the discrimination; UNC
  the base rate's; Brier = MCB − DSC + UNC. The fit is monotone in the forecast alone, so biases of opposite sign at
  one forecast value in two places cancel inside it;
- OVER and UNDER split the forecasts' departure from the fit by direction (gDNA probability above, or below, the
  frequency at which fragments forecast there are gDNA); OVER − UNDER is the net gDNA pool error over the forecast
  fragments, so the pair shows what a net error cancels;
- the three-way Brier scores gDNA, annotated and nascent jointly. A single-exon annotated transcript that serves as
  its own nascent entity books its nascent fragments as annotated, so the three-way score carries that asymmetry as
  ``quant_accuracy.py``'s pool rows do;
- where the wrong-source mass lies: per source, the quartiles of the forecast weighted by the wrong-source mass,
  and that mass's share forecast at certainty (``p`` exactly 0 or 1).

A record is refused unless every check holds, each tolerance a rounding bound (Higham, Accuracy and Stability of
Numerical Algorithms, 2nd ed., §3.1 and §4.2; ``u = 2^-53``, ``γ(m) = m u / (1 − m u)``), never a chosen threshold:

- the read pairs, the excluded, the placed and the intergenic fragments equal the scan's own counts exactly, and the
  annotated BAM's ``ZF`` equals the annotation table's flag on every fragment the table holds: the join of fragment
  ids to read names;
- each EM unit's masses sum to 1 within ``γ(3n + 2)`` for its ``n`` candidates, or to exactly 0;
- the winner's class holds the winner's float32 posterior, so the masses are the annotation's rows;
- each locus's gDNA equals the solver's exactly (the same additions in the same order), its RNA and the intergenic
  term within two evaluations' bound;
- the closure: the four class sums equal the pools ``quant_accuracy.py`` reports for the same run (``score_library``:
  ``mrna_est``, ``nrna_est``, ``gdna_est``, ``n_intergenic_gdna``), each within two evaluations' bound at the
  check's own term count. This holds under fractional assignment.

The protocol is pinned and fractional (``scan.total_threads=1``, ``em.assignment_mode=fractional``), so a pair is
exactly reproducible and the closure applies; the instrument sets both and refuses a ``--set`` that undoes either.
Each record carries ``quant_accuracy.row_stamp``: the code revision, the settings and the input. A refused record
writes its JSON with ``accepted: false`` and no fragment file. ``--pair`` reads accepted records only, checks each
fragment file against the digest its JSON holds, and refuses two records with different inputs or conditions, a read
name held twice, different read names, or one read name with two sources. ``--report`` prints the tables per stratum
and rung (``g00``, contaminated), each summed over its conditions and never across strata.

A record holds every read pair of its condition (a 10 M-fragment condition writes about 140 MB): it is bulk, deleted
once paired; the pair's JSON lines are the receipt. Recording takes about two and a half times ``quant_accuracy.py``'s
run of the same condition (the annotated BAM, and one pass over it for the read names); pairing takes under a minute.

Gate: ``tests/test_paired_score.py``; ``--self-test`` plants each defect and requires its check to fire.

Usage::

    # 1. record each arm (one process per arm; a second build records under its own interpreter)
    python scripts/design/paired_score.py --record --conditions COND --out $RIGEL_ARMS/ps_base
    python scripts/design/paired_score.py --record --conditions COND --set calibration.message_policy=silent \\
        --out $RIGEL_ARMS/ps_silent
    python scripts/design/paired_score.py --record --conditions COND --arm no_capture --out $RIGEL_ARMS/ps_no_capture
    # 2. join them by read name, per condition, and print the tables
    python scripts/design/paired_score.py --pair $RIGEL_ARMS/ps_base $RIGEL_ARMS/ps_silent \\
        --out $RIGEL_ARMS/ps_base_silent.jsonl
    # 3. re-print a pair file, or render it as markdown
    python scripts/design/paired_score.py --report $RIGEL_ARMS/ps_base_silent.jsonl --markdown F.md
    python scripts/design/paired_score.py --self-test

``--bam`` records from an aligner's BAM of the same reads (``panel.py align``); both arms of a pair must read one
input.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import time
from array import array
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from _shared import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_SUITE,
    ORACLE_BAM,
    SCRATCH,
    add_set_option,
    configured,
    is_zero_gdna,
    panel_conditions,
    refuses,
    run_shards,
    scope,
    set_args,
    sibling,
    strata,
    stratum,
)

from rigel.native import EM_CLASSES, EM_WINNER_INTERGENIC  # noqa: E402

#: The class-mass columns, by the native names: annotated transcripts, synthetic spans, the locus's gDNA, intergenic.
ANN, NAS, GDNA, IG = (EM_CLASSES.index(c) for c in ("annotated", "nascent", "gdna", "intergenic"))

#: The read name's source, as the pools name it: gDNA, annotated RNA (a reference transcript's mature sequence),
#: nascent RNA. ``parse_origin`` spells them ``gdna``, ``mrna``, ``nrna``.
KINDS = ("gdna", "annotated", "nascent")
_KIND_OF_ORIGIN = {"gdna": 0, "mrna": 1, "nrna": 2}

#: A read pair's status in a record. Only ``forecast`` is scored.
STATUSES = ("forecast", "unassigned", "withheld", "excluded")
FORECAST, UNASSIGNED, WITHHELD, EXCLUDED = range(len(STATUSES))

#: The ``quant_accuracy.py`` arms a record can run under: the ones that substitute no truth, so they record from any
#: BAM of the condition's reads. ``base`` substitutes nothing.
ARMS = ("base", "no_capture")

#: The protocol every arm runs under: the scan pinned, so the pair is exactly reproducible, and the assignment
#: fractional, so the class sums close on the pools. Applied before the arm's own ``--set``.
PROTOCOL = ("scan.total_threads=1", "em.assignment_mode=fractional")

#: The SAM flag bits the scan's record filter reads (``bam_record.h``, ``is_filtered``): QC-failed, unmapped, and a
#: duplicate when duplicates are skipped.
_FLAG_UNMAPPED, _FLAG_QCFAIL, _FLAG_DUPLICATE = 0x4, 0x200, 0x400

#: The read pair's ZF when no eligible record carries one.
_NO_ZF = -1

#: Double and single precision's unit roundoff (round to nearest).
U64 = 2.0**-53
U32 = 2.0**-24

#: The quartiles reported for where the wrong-source mass lies.
QUARTILES = (0.25, 0.5, 0.75)
#: The wrong-source mass is also kept on a grid over the forecast, so a stratum's quartiles can be read from its
#: conditions pooled; the grid is a display resolution (the per-condition quartiles are exact).
GRID = 1000

#: The rows of every score dict that are sums over fragments, and so add across a stratum's conditions.
_SUMMED = ("brier", "linear", "brier3", "mcb", "dsc", "unc", "over", "under", "net_pool_error",
           *(f"{s}_{k}" for s in ("brier", "linear", "brier3") for k in KINDS))


# ── the rounding bounds ──────────────────────────────────────────────────────────────────────────


def gamma(m):
    """``γ(m) = m u / (1 − m u)``: a sum of nonnegative terms, evaluated by any order of additions in which no term
    passes through more than ``m`` of them, is within ``γ(m)`` of its exact value, relatively."""
    m = np.asarray(m, np.float64)
    return m * U64 / (1.0 - m * U64)


def two_order_bound(m, s1, s2):
    """The largest difference two evaluations of one nonnegative sum can show, each of height at most ``m``:
    ``2 γ(m) / (1 − γ(m)) · max(|s1|, |s2|)``."""
    g = gamma(m)
    return 2.0 * g / (1.0 - g) * np.maximum(np.abs(s1), np.abs(s2))


def conservation_bound(n_candidates):
    """How far one EM unit's four masses may sum from 1, for ``n`` candidates: the solver computes ``e_k``, their
    sequential sum ``S`` (height ``n − 1``), ``inv = fl(1 / S)`` and ``p_k = fl(e_k · inv)``, and adds each ``p_k``
    into its class from 0 (height at most ``n − 1``); the four classes add with height 3. The ``e_k`` cancel between
    numerator and denominator, so ``|row − 1| ≤ (1 + u)² (1 + γ(n + 2)) / (1 − γ(n − 1)) − 1 ≤ γ(3n + 2)``."""
    return gamma(3 * np.asarray(n_candidates, np.int64) + 2)


def winner_bound(held, n_candidates):
    """The most a winner's float32 posterior may exceed its class's double mass, which holds that same posterior:
    ``w32 ≤ (1 + u32) p_w ≤ (1 + u32) held / (1 − γ(n))``, and ``γ(3)`` more covers evaluating the bound."""
    return held * (1.0 + U32) / (1.0 - gamma(np.asarray(n_candidates, np.int64))) * (1.0 + gamma(3))


# ── the scores ───────────────────────────────────────────────────────────────────────────────────


def p_gdna(mass: np.ndarray) -> np.ndarray:
    """The source forecast of each row: its gDNA classes over its sum (NaN for a row with no mass)."""
    row = mass.sum(1)
    return np.divide(mass[:, GDNA] + mass[:, IG], row, out=np.full(row.size, np.nan), where=row > 0)


def scores(mass: np.ndarray, kind: np.ndarray):
    """``(scores, the CORP fit's blocks)`` of forecasts ``mass`` (``[n, 4]``, every row a forecast) against the read
    names' sources ``kind`` (``[n]``, 0 gDNA, 1 annotated, 2 nascent); ``(None, None)`` for no fragments."""
    from scipy.optimize import isotonic_regression

    kind = np.asarray(kind, np.int64)
    if kind.size == 0:
        return None, None
    row = mass.sum(1)
    if not np.all(row > 0.0):
        raise ValueError("a row with no forecast reached the scores")
    p = (mass[:, GDNA] + mass[:, IG]) / row
    q = np.stack([p, mass[:, ANN] / row, mass[:, NAS] / row], 1)  # gDNA, annotated, nascent
    y = (kind == 0).astype(np.float64)
    b = (p - y) ** 2
    onehot = np.zeros_like(q)
    onehot[np.arange(kind.size), kind] = 1.0
    b3 = ((q - onehot) ** 2).sum(1)
    out = {"n_forecast": int(kind.size)}
    for i, k in enumerate(KINDS):
        sel = kind == i
        out[f"n_{k}"] = int(sel.sum())
        out[f"linear_{k}"] = float(np.abs(p[sel] - y[sel]).sum())
        out[f"brier_{k}"] = float(b[sel].sum())
        out[f"brier3_{k}"] = float(b3[sel].sum())
    out["brier"] = float(b.sum())
    out["brier3"] = float(b3.sum())
    out["linear"] = float(np.abs(p - y).sum())

    # CORP within this condition: ties pooled, integer counts per forecast value, the fit's blocks from scipy, each
    # block's value its exact count ratio
    up, inv = np.unique(p, return_inverse=True)
    n = np.bincount(inv, minlength=up.size)
    s = np.bincount(inv[kind == 0], minlength=up.size)
    blk = isotonic_regression(s / n, weights=n.astype(np.float64), increasing=True).blocks
    n_b, s_b = np.add.reduceat(n, blk[:-1]), np.add.reduceat(s, blk[:-1])  # exact integers
    fit_b = s_b / n_b
    fit = np.repeat(fit_b, np.diff(blk))
    nf, sf = n.astype(np.float64), s.astype(np.float64)
    bs = float(np.sum(sf * (1.0 - up) ** 2 + (nf - sf) * up**2))
    bs_r = float(np.sum(s_b * (1.0 - fit_b) ** 2 + (n_b - s_b) * fit_b**2))
    n_all, y_all = float(n.sum()), float(s.sum())
    unc = n_all * (y_all / n_all) * (1.0 - y_all / n_all)
    out.update({
        "mcb": bs - bs_r,
        "dsc": unc - bs_r,
        "unc": unc,
        "over": float(np.sum(nf * np.clip(up - fit, 0.0, None))),
        "under": float(np.sum(nf * np.clip(fit - up, 0.0, None))),
        "net_pool_error": float(p.sum() - y_all),
        "n_unique_forecasts": int(up.size),
        "n_blocks": int(n_b.size),
    })
    # the decomposition's Brier, summed per forecast value, is the per-fragment Brier: a term takes at most four
    # roundings (a difference, a square, a product, the sum of the two outcomes) and each sum has height at most n
    out["check_brier"] = float(abs(bs - out["brier"]))
    out["check_brier_bound"] = float(two_order_bound(kind.size + 4, bs, out["brier"]))
    blocks = {"p_lo": up[blk[:-1]], "p_hi": up[blk[1:] - 1], "fit": fit_b, "n": n_b, "gdna": s_b}
    return out, blocks


def wquantiles(x: np.ndarray, w: np.ndarray, qs=QUARTILES) -> list:
    """The lower weighted quantiles of ``x``: the smallest ``x`` whose cumulative weight reaches ``q`` of the total."""
    if x.size == 0 or w.sum() <= 0:
        return [None for _ in qs]
    o = np.argsort(x, kind="stable")
    cw = np.cumsum(w[o])
    return [float(x[o][min(int(np.searchsorted(cw, q * cw[-1])), x.size - 1)]) for q in qs]


def where_they_lie(p: np.ndarray, kind: np.ndarray) -> dict:
    """Per source: the fragments forecast at certainty, and the quartiles of the forecast weighted by the
    wrong-source mass (gDNA called RNA, or RNA called gDNA). For a calibrated forecast ``E[(p − Y)² | p] = p(1 − p)``
    and ``E[|p − Y| | p] = 2p(1 − p)``, so a Brier of half the linear score holds at every ``p`` and locates nothing;
    this does."""
    out = {}
    for i, k in enumerate(KINDS):
        s = kind == i
        if not s.any():
            continue
        pk = p[s]
        wrong = (1.0 - pk) if i == 0 else pk
        out[k] = {
            "n": int(s.sum()),
            "at_0": int((pk == 0.0).sum()),
            "at_1": int((pk == 1.0).sum()),
            "wrong_mass": float(wrong.sum()),
            "wrong_mass_at_certainty": float(wrong[(pk == 0.0) | (pk == 1.0)].sum()),
            "wrong_mass_quartiles_of_p": wquantiles(pk, wrong),
            "wrong_mass_grid": np.bincount(np.minimum((pk * GRID).astype(np.int64), GRID - 1), weights=wrong,
                                           minlength=GRID).tolist(),
        }
    return out


# ── one arm's run: what the hooks capture ────────────────────────────────────────────────────────


@dataclasses.dataclass
class Run:
    """What one pipeline run hands its record: the EM's units, the annotation table's rows, the scan's counts and
    the pools ``quant_accuracy.score_library`` reports. Everything the record composes and checks, and nothing read
    from disk, so the self-test can plant a defect in any of it."""

    #: per EM unit, in the class masses' row order: fragment id, locus, the four masses, the native winner and its
    #: posterior, the candidate count
    em_fid: np.ndarray
    em_lid: np.ndarray
    em_mass: np.ndarray
    em_wtid: np.ndarray
    em_wpost: np.ndarray
    em_ncand: np.ndarray
    #: the annotation table's rows: fragment id, winner, ZF flag, posterior, class, candidate count
    t_fid: np.ndarray
    t_tid: np.ndarray
    t_flags: np.ndarray
    t_post: np.ndarray
    t_class: np.ndarray
    t_ncand: np.ndarray
    #: per transcript: whether it is a synthetic span
    is_synthetic: np.ndarray
    #: the solver's per-locus totals, ``{locus id: (gDNA, RNA)}``, and its intergenic gDNA
    locus_totals: dict
    intergenic_em: float
    #: the scan's own counts: ``n_read_pairs``, ``n_excluded``, ``n_fragments`` (placed), ``n_intergenic``
    ledger: dict
    #: ``quant_accuracy.score_library``'s pools: ``annotated``, ``nascent``, ``gdna``, ``intergenic``
    pools: dict
    #: the terms a pool's sum can pass through beyond the units' candidates and the routed fragments
    n_extra_terms: int


@dataclasses.dataclass
class ReadPairs:
    """The input BAM's read pairs in file order, as the annotated-BAM pass wrote them: each one's read-name digest, its
    source, whether it holds an eligible record (and so a fragment id), and the ZF its eligible records carry."""

    key: np.ndarray
    kind: np.ndarray
    eligible: np.ndarray
    zf: np.ndarray
    #: read pairs whose eligible records disagree on ZF, and filtered records carrying one (the writer passes a
    #: filtered record through untouched)
    zf_inconsistent: int = 0
    zf_on_filtered: int = 0


def name_key(name: str) -> int:
    """A read name's 64-bit BLAKE2b digest, the key two records are joined on (a digest repeated within one record is
    refused, so no collision passes unseen)."""
    return int.from_bytes(hashlib.blake2b(name.encode(), digest_size=8).digest(), "little")


def origin_code(name: str) -> int:
    """The read name's source as a ``KINDS`` index, by the simulator's own parser."""
    from rigel.sim.read_name import parse_origin

    return _KIND_OF_ORIGIN[parse_origin(name).kind]


def filtered(flag: int, skip_duplicates: bool) -> bool:
    """The scan's record filter (``is_filtered``): a filtered record gives its read pair no fragment id."""
    return bool(flag & (_FLAG_QCFAIL | _FLAG_UNMAPPED)) or bool((flag & _FLAG_DUPLICATE) and skip_duplicates)


def read_pairs(records, skip_duplicates: bool, kind_of=origin_code) -> ReadPairs:
    """The read pairs of ``records``, ``(name, flag, ZF or None)`` in file order, by the scan's rule: a read pair is a
    run of records sharing a name, and one with an eligible record gets the next fragment id."""
    keys, kinds, eligible, zfs = array("Q"), array("b"), array("b"), array("h")
    bad_zf = foreign_zf = 0
    name, elig, zf = None, False, _NO_ZF
    for qname, flag, tag in records:
        if qname != name:
            if name is not None:
                keys.append(name_key(name))
                kinds.append(kind_of(name))
                eligible.append(elig)
                zfs.append(zf)
            name, elig, zf = qname, False, _NO_ZF
        if filtered(flag, skip_duplicates):
            foreign_zf += tag is not None
            continue
        elig = True
        if tag is not None:
            bad_zf += zf not in (_NO_ZF, tag)
            zf = int(tag) if zf == _NO_ZF else zf
    if name is not None:
        keys.append(name_key(name))
        kinds.append(kind_of(name))
        eligible.append(elig)
        zfs.append(zf)
    return ReadPairs(
        key=np.frombuffer(keys, np.uint64).copy(),
        kind=np.frombuffer(kinds, np.int8).copy(),
        eligible=np.frombuffer(eligible, np.int8).astype(bool),
        zf=np.frombuffer(zfs, np.int16).astype(np.int32),
        zf_inconsistent=bad_zf,
        zf_on_filtered=foreign_zf,
    )


def _bam_records(path: Path):
    """``(name, flag, ZF or None)`` of every record of a BAM, in file order."""
    import pysam

    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            yield r.query_name, r.flag, (r.get_tag("ZF") if r.has_tag("ZF") else None)


# ── the record: compose every read pair's forecast, and check it ─────────────────────────────────


def compose(run: Run, pairs: ReadPairs):
    """``(status, mass, checks, closure, counts)`` per read pair of ``pairs``, in its order.

    A fragment id is a read pair's rank among the eligible ones. The EM's units take their class masses; the table's
    other rows are routing's (deterministic, gated, intergenic); a placed fragment the table does not hold and the
    writer stamped intergenic is the fast path. ``checks`` holds every acceptance check, ``closure`` the class sums
    against the pools."""
    from rigel.annotate import AF_ANNOTATED, AF_GDNA_INTERGENIC, AF_MULTIMAPPER_DROP, AF_NASCENT, AF_UNRESOLVED
    from rigel.buffer import FRAG_UNAMBIG

    checks: dict = {}
    pos = np.flatnonzero(pairs.eligible)  # the read pair of each fragment id
    n_frag = int(pos.size)
    fids = np.concatenate([run.em_fid, run.t_fid])
    checks["fragment_ids_out_of_range"] = int(np.sum((fids < 0) | (fids >= n_frag)))
    if checks["fragment_ids_out_of_range"]:
        return None, None, checks, {}, {}

    # -- the EM's units ----------------------------------------------------------------------------------------------
    in_em = np.zeros(n_frag, bool)
    in_em[run.em_fid] = True
    checks["em_units_repeated"] = int(run.em_fid.size - np.unique(run.em_fid).size)
    # -- routing's rows: the table's rows the EM does not hold -----------------------------------------------------------
    routed = ~in_em[run.t_fid]
    r_fid, r_tid, r_flags = run.t_fid[routed], run.t_tid[routed], run.t_flags[routed]
    r_post, r_class, r_ncand = run.t_post[routed], run.t_class[routed], run.t_ncand[routed]
    checks["table_rows_repeated"] = int(run.t_fid.size - np.unique(run.t_fid).size)
    checks["em_units_not_in_table"] = int(run.em_fid.size - int((~routed).sum()))
    det = (r_class == FRAG_UNAMBIG) & (r_tid >= 0) & (r_post == 1.0) & np.isin(r_flags, (AF_ANNOTATED, AF_NASCENT))
    gated = (r_flags == AF_UNRESOLVED) & (r_post == 0.0)
    inter = (r_flags == AF_GDNA_INTERGENIC) & (r_ncand == 0) & (r_post == 1.0)
    checks["routed_rows_of_no_kind"] = int((~(det | gated | inter)).sum())
    det_synth = run.is_synthetic[np.clip(r_tid[det], 0, None)]
    checks["deterministic_flag_not_its_class"] = int(np.sum((r_flags[det] == AF_NASCENT) != det_synth))

    mass = np.zeros((n_frag, len(EM_CLASSES)))
    has_row = np.zeros(n_frag, bool)
    mass[run.em_fid] = run.em_mass
    has_row[run.em_fid] = True
    mass[r_fid[det], np.where(det_synth, NAS, ANN)] = 1.0
    mass[r_fid[inter], IG] = 1.0
    has_row[r_fid] = True

    # -- the writer's ZF: the id join, and the fast path ------------------------------------------------------------------
    zf = pairs.zf[pos]
    checks["zf_not_the_tables_flag"] = int(np.sum(zf[run.t_fid] != run.t_flags.astype(np.int32)))
    checks["zf_inconsistent_within_a_read_pair"] = int(pairs.zf_inconsistent)
    checks["zf_on_a_filtered_record"] = int(pairs.zf_on_filtered)
    fast = ~has_row & (zf == AF_GDNA_INTERGENIC)
    mass[fast, IG] = 1.0
    placed = has_row | fast
    checks["withheld_with_another_zf"] = int(np.sum(~placed & ~np.isin(zf, (AF_UNRESOLVED, AF_MULTIMAPPER_DROP))))

    # -- the scan's own counts, exactly -----------------------------------------------------------------------------------
    n_pairs, n_excluded = int(pairs.key.size), int((~pairs.eligible).sum())
    ledger = {
        "read_pairs": (n_pairs, run.ledger["n_read_pairs"]),
        "excluded": (n_excluded, run.ledger["n_excluded"]),
        "placed": (int(placed.sum()), run.ledger["n_fragments"]),
        "intergenic": (int(fast.sum()) + int(inter.sum()), run.ledger["n_intergenic"]),
    }
    checks["ledger"] = {k: {"record": a, "scan": b} for k, (a, b) in ledger.items()}
    checks["ledger_mismatches"] = sum(a != b for a, b in ledger.values())

    # -- each EM unit conserves its mass ------------------------------------------------------------------------------------
    em_rows = run.em_mass.sum(1)
    zero = em_rows == 0.0
    resid = np.abs(em_rows[~zero] - 1.0)
    cb = conservation_bound(run.em_ncand[~zero])
    checks["em_units"] = int(em_rows.size)
    checks["em_units_with_no_mass"] = int(zero.sum())
    checks["em_units_with_no_mass_with_a_winner"] = int((zero & (run.em_wtid != -1)).sum())
    checks["conservation_over_bound"] = int((resid > cb).sum())
    checks["conservation_worst_ratio"] = float((resid / cb).max(initial=0.0))
    checks["negative_mass"] = int((run.em_mass < 0.0).sum())

    # -- the winner's class holds the winner's posterior: the masses are the annotation's rows ------------------------------
    won = run.em_wpost > 0
    wtid = run.em_wtid[won]
    wcls = np.where(wtid >= 0, np.where(run.is_synthetic[np.clip(wtid, 0, None)], NAS, ANN),
                    np.where(wtid == EM_WINNER_INTERGENIC, IG, GDNA))
    held = run.em_mass[np.flatnonzero(won), wcls]
    checks["winner_above_its_class"] = int(np.sum(run.em_wpost[won].astype(np.float64)
                                                  > winner_bound(held, run.em_ncand[won])))

    # -- per locus, the solver's totals ------------------------------------------------------------------------------------
    lids = np.array(sorted(run.locus_totals), np.int64)
    n_l = int(max(lids.max(initial=-1), run.em_lid.max(initial=-1))) + 1
    by_g = np.bincount(run.em_lid, weights=run.em_mass[:, GDNA], minlength=n_l)  # in unit order, as the solver adds
    by_r = np.bincount(run.em_lid, weights=run.em_mass[:, ANN] + run.em_mass[:, NAS], minlength=n_l)
    nat_g = np.array([run.locus_totals[i][0] for i in lids])
    nat_r = np.array([run.locus_totals[i][1] for i in lids])
    terms = (np.bincount(run.em_lid, weights=run.em_ncand, minlength=n_l)
             + np.bincount(run.em_lid, minlength=n_l) + 1)[lids]
    checks["locus_gdna_not_exact"] = int(np.sum(by_g[lids] != nat_g))
    checks["locus_rna_over_bound"] = int(np.sum(np.abs(by_r[lids] - nat_r) > two_order_bound(terms, by_r[lids], nat_r)))
    checks["units_outside_every_locus"] = int(np.sum(~np.isin(run.em_lid, lids)))
    ig_em = float(run.em_mass[:, IG].sum())
    ig_bound = float(two_order_bound(int(run.em_ncand.sum()) + lids.size + 1, ig_em, run.intergenic_em))
    checks["intergenic_em"] = {"mass": ig_em, "solver": float(run.intergenic_em),
                               "diff": abs(ig_em - float(run.intergenic_em)), "bound": ig_bound}

    # -- the closure: the class sums are the pools quant_accuracy.py reports -----------------------------------------------
    tot = mass.sum(0)
    t_pool = int(run.em_ncand.sum()) + int(routed.sum()) + int(fast.sum()) + int(run.n_extra_terms)
    closure = {}
    for c, name in enumerate(EM_CLASSES):
        m, p = float(tot[c]), float(run.pools[name])
        closure[name] = {"mass": m, "pool": p, "diff": abs(m - p), "bound": float(two_order_bound(t_pool, m, p))}
    checks["closure_over_bound"] = sum(v["diff"] > v["bound"] for v in closure.values())

    # -- every read pair's status and forecast ----------------------------------------------------------------------------
    status = np.full(n_pairs, EXCLUDED, np.int8)
    status[pos] = np.where(placed, np.where(mass.sum(1) > 0.0, FORECAST, UNASSIGNED), WITHHELD)
    full = np.zeros((n_pairs, len(EM_CLASSES)))
    full[pos] = mass
    counts = {
        "em_units": int(run.em_fid.size),
        "deterministic": int(det.sum()),
        "gated": int(gated.sum()),
        "routed_intergenic": int(inter.sum()),
        "fast_path_intergenic": int(fast.sum()),
        "withheld_by_zf": {str(int(z)): int(np.sum(~placed & (zf == z))) for z in np.unique(zf[~placed])},
    }
    return status, full, checks, closure, counts


def accepted(checks: dict) -> bool:
    """Every acceptance check of a record holds."""
    zero = ("fragment_ids_out_of_range", "em_units_repeated", "table_rows_repeated", "em_units_not_in_table",
            "routed_rows_of_no_kind", "deterministic_flag_not_its_class", "zf_not_the_tables_flag",
            "zf_inconsistent_within_a_read_pair", "zf_on_a_filtered_record", "withheld_with_another_zf",
            "ledger_mismatches", "conservation_over_bound", "negative_mass", "winner_above_its_class",
            "locus_gdna_not_exact", "locus_rna_over_bound", "units_outside_every_locus", "closure_over_bound",
            "read_names_repeated", "scores_check_over_bound")
    ig = checks.get("intergenic_em", {})
    return all(checks.get(k, 1) == 0 for k in zero) and ig.get("diff", 1.0) <= ig.get("bound", 0.0)


def counts_by(status: np.ndarray, kind: np.ndarray) -> dict:
    """Read pairs per status and source."""
    return {s: {k: int(np.sum((status == i) & (kind == j))) for j, k in enumerate(KINDS)}
            for i, s in enumerate(STATUSES)}


def build_record(run: Run, pairs: ReadPairs) -> tuple[dict, dict]:
    """``(arrays, summary)``: the record's per-read-pair arrays, and its checks, closure, counts and own scores."""
    status, mass, checks, closure, counts = compose(run, pairs)
    checks["read_names_repeated"] = int(pairs.key.size - np.unique(pairs.key).size)
    own = None
    if status is not None:
        f = status == FORECAST
        own, _ = scores(mass[f], pairs.kind[f])
        checks["scores_check_over_bound"] = int(own is not None and own["check_brier"] > own["check_brier_bound"])
    summary = {
        "accepted": accepted(checks),
        "checks": checks,
        "closure": closure,
        "counts": counts,
        "status": None if status is None else counts_by(status, pairs.kind),
        "in_bam": {k: int(np.sum(pairs.kind == j)) for j, k in enumerate(KINDS)},
        "scores": own,
    }
    arrays = {"key": pairs.key, "kind": pairs.kind, "status": status, "mass": mass}
    return arrays, summary


# ── one arm on one condition ─────────────────────────────────────────────────────────────────────


def capture_run(bam: Path, index, pipeline_config, annotated_bam: Path):
    """Run the pipeline once with the annotated BAM written to ``annotated_bam``; return ``(result, units, class
    masses, table rows)``: the EM's units (fragment id, locus, winner, its posterior, candidate count) in the order of
    their class masses (``emit_class_mass``), and the annotation table's rows. Each hook passes its call through and
    counts it; a run that fires the table hook other than once, or the EM's two hooks unequally, is refused."""
    import rigel.pipeline as PL
    from rigel.estimator import AbundanceEstimator

    got = {"table": [], "mass": [], "units": []}
    orig_table, orig_populate = PL._annotation_table, PL.populate_em_annotations
    orig_em = AbundanceEstimator.run_batch_locus_em_partitioned

    def table(buffer):
        t = orig_table(buffer)
        got["table"].append(t)
        return t

    def em(self, *args, **kw):
        if kw.get("emit_class_mass"):
            raise RuntimeError("the pipeline asks for the class masses itself; the hook would take them twice")
        out = orig_em(self, *args, **{**kw, "emit_class_mass": True})
        got["mass"].append(np.asarray(out[-1], np.float64))
        return out[:-1]

    def populate(batch_parts, out_winner_tid, out_winner_post, out_n_candidates, annotations, idx):
        got["units"].append({
            "fid": np.concatenate([p.frag_ids for p in batch_parts]).astype(np.int64),
            "lid": np.concatenate([np.full(p.frag_ids.shape[0], p.locus_id, np.int64) for p in batch_parts]),
            "wtid": np.asarray(out_winner_tid, np.int32).copy(),
            "wpost": np.asarray(out_winner_post, np.float32).copy(),
            "ncand": np.asarray(out_n_candidates, np.int64).copy(),
        })
        return orig_populate(batch_parts, out_winner_tid, out_winner_post, out_n_candidates, annotations, idx)

    PL._annotation_table, PL.populate_em_annotations = table, populate
    AbundanceEstimator.run_batch_locus_em_partitioned = em
    try:
        result = PL.run_pipeline(str(bam), index, dataclasses.replace(pipeline_config, annotated_bam_path=annotated_bam))
    finally:
        PL._annotation_table, PL.populate_em_annotations = orig_table, orig_populate
        AbundanceEstimator.run_batch_locus_em_partitioned = orig_em
    calls = {k: len(v) for k, v in got.items()}
    if calls["table"] != 1 or calls["mass"] != calls["units"] or calls["mass"] > 1:
        raise RuntimeError(f"the hooks fired {calls} times; a run fires the table once and the EM at most once, its "
                           "units with it. This is not a measurement.")
    units = got["units"][0] if got["units"] else None
    em_mass = got["mass"][0] if got["mass"] else np.zeros((0, len(EM_CLASSES)))
    if units is not None and em_mass.shape[0] != units["fid"].size:
        raise RuntimeError(f"{em_mass.shape[0]} class-mass rows for {units['fid'].size} EM units")
    t = got["table"][0]
    n = t.size
    table_rows = {k: np.asarray(getattr(t, k)[:n]).copy()
                  for k in ("frag_ids", "best_tid", "tx_flags", "posterior", "frag_class", "n_candidates")}
    return result, units, em_mass, table_rows


def record_condition(arm: str, suite: Path, index, condition: str, pipeline_config, bam_name: str,
                     work_dir: Path, qa_arm: str = "base") -> tuple[dict, dict]:
    """One arm on one condition of a panel: ``(arrays, row)``."""
    bam = suite / condition / bam_name
    if not bam.is_file():
        raise SystemExit(f"⛔ {condition}: no {bam_name}")
    summary = json.loads((suite / condition / "truth_summary.json").read_text())
    return record_bam(arm, bam, index, pipeline_config, work_dir, summary, condition, qa_arm)


def install_record_arm(qa, qa_arm: str):
    """Install the record's ``quant_accuracy.py`` arm (:data:`ARMS`) with that instrument's own code: nothing for
    ``base``, every capture efficiency read as 1 for ``no_capture``. Returns ``(restore, fired)``."""
    if qa_arm not in ARMS:
        raise SystemExit(f"⛔ arm {qa_arm!r} is not one a record runs under ({', '.join(ARMS)})")
    return qa.install_no_capture() if qa_arm == "no_capture" else qa.install_arm("base", None)


def record_bam(arm: str, bam: Path, index, pipeline_config, work_dir: Path, summary: dict,
               condition: str, qa_arm: str = "base") -> tuple[dict, dict]:
    """One arm on one BAM: ``(arrays, row)``, the row stamped as ``quant_accuracy.py`` stamps its rows, under the
    ``quant_accuracy.py`` arm ``qa_arm``. ``summary`` is the condition's ``truth_summary.json``, whose
    ``origin_counts`` are the simulator's templates by source. A ``no_capture`` run that read any capture is refused
    (``check_no_capture``), as is an arm whose wrapper never ran."""
    from rigel.splice import NUM_SPLICE_STRAND_COLS

    qa = sibling("quant_accuracy.py")
    bam_name = Path(bam).name
    work_dir.mkdir(parents=True, exist_ok=True)
    annotated = work_dir / f".{condition}.{os.getpid()}.annotated.bam"
    start = time.perf_counter()
    restore, fired = install_record_arm(qa, qa_arm)
    try:
        result, units, em_mass, t = capture_run(bam, index, pipeline_config, annotated)
        quant = result.estimator.get_counts_df(index)
        lib = qa.score_library(result, quant, summary)
        pairs = read_pairs(_bam_records(annotated), pipeline_config.scan.skip_duplicates)
    finally:
        restore()
        annotated.unlink(missing_ok=True)
    if fired["n"] == 0:
        raise RuntimeError(f"{condition} [{qa_arm}]: the arm's wrapper never ran, so the record is the shipped run's")
    if qa_arm == "no_capture":
        qa.check_no_capture(condition, result.calibration)
    est, stats = result.estimator, result.stats
    empty = np.zeros(0, np.int64)
    run = Run(
        em_fid=units["fid"] if units else empty,
        em_lid=units["lid"] if units else empty,
        em_mass=em_mass,
        em_wtid=units["wtid"] if units else empty.astype(np.int32),
        em_wpost=units["wpost"] if units else empty.astype(np.float32),
        em_ncand=units["ncand"] if units else empty,
        t_fid=t["frag_ids"].astype(np.int64),
        t_tid=t["best_tid"].astype(np.int64),
        t_flags=t["tx_flags"].astype(np.int64),
        t_post=t["posterior"].astype(np.float64),
        t_class=t["frag_class"],
        t_ncand=t["n_candidates"].astype(np.int64),
        is_synthetic=index.t_df["is_synthetic"].to_numpy().astype(bool),
        locus_totals={int(r["locus_id"]): (float(r["count_gdna"]), float(r["rna_total"]))
                      for r in est.locus_results},
        intergenic_em=float(est.gdna_intergenic_em_count),
        ledger={k: int(getattr(stats, k)) for k in ("n_read_pairs", "n_excluded", "n_fragments", "n_intergenic")},
        pools={"annotated": lib["mrna_est"], "nascent": lib["nrna_est"], "gdna": lib["gdna_est"],
               "intergenic": lib["n_intergenic_gdna"]},
        # every pool's other side also adds per transcript, per strand column, per locus, and once more
        n_extra_terms=int(index.num_transcripts) + len(est.locus_results) + NUM_SPLICE_STRAND_COLS + 1,
    )
    arrays, rec = build_record(run, pairs)
    truth = {k: int(summary.get("origin_counts", {}).get(o, 0)) for k, o in zip(KINDS, ("gdna", "mrna", "nrna"))}
    row = {
        **qa.row_stamp(arm, condition, pipeline_config, bam_name),
        "qa_arm": qa_arm,
        "seconds": time.perf_counter() - start,
        "stratum": list(stratum(condition)),
        "scope": scope(stratum(condition)),
        **rec,
        "truth": truth,
        "not_in_bam": {k: truth[k] - rec["in_bam"][k] for k in KINDS},
        "library": {k: lib[k] for k in ("mrna_est", "nrna_est", "gdna_est", "n_intergenic_gdna", "gdna_frac_est",
                                        "gdna_frac_true")},
    }
    return arrays, row


def protocol_config(settings):
    """The arm's ``PipelineConfig``: :data:`PROTOCOL`, then its own ``--set`` specs. A spec that undoes the protocol
    is refused: unpinned, the pair is not reproducible; sampled, the class sums cannot close on the pools."""
    cfg = configured([*PROTOCOL, *settings])
    if cfg.scan.total_threads != 1 or cfg.em.assignment_mode != "fractional":
        raise SystemExit(f"⛔ --set undoes the protocol ({', '.join(PROTOCOL)}): the arm would be unpinned or "
                         "sampled, and its pair not reproducible or its closure undefined")
    return cfg


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def write_record(out: Path, condition: str, arrays: dict, row: dict) -> None:
    """``<out>/<condition>.npz`` (an accepted record's arrays) and ``<out>/<condition>.json`` (its row, with the
    arrays' digest); a refused record writes only its row, and removes an earlier record's arrays."""
    out.mkdir(parents=True, exist_ok=True)
    npz, meta = out / f"{condition}.npz", out / f"{condition}.json"
    if row["accepted"]:
        tmp = out / f".{condition}.{os.getpid()}.npz"
        with tmp.open("wb") as fh:
            np.savez_compressed(fh, **arrays)
        tmp.replace(npz)
        row = {**row, "fragments_file": npz.name, "fragments_sha256": _sha256(npz)}
    else:
        npz.unlink(missing_ok=True)
    tmp = out / f".{condition}.{os.getpid()}.json"
    tmp.write_text(json.dumps(row, indent=1))
    tmp.replace(meta)


# ── the pair ─────────────────────────────────────────────────────────────────────────────────────


def load_record(rec_dir: Path, condition: str) -> dict:
    """A record's row and arrays; refused unless it was accepted and its arrays are the ones its row describes."""
    meta = json.loads((rec_dir / f"{condition}.json").read_text())
    if not meta.get("accepted"):
        raise SystemExit(f"⛔ {rec_dir / condition}: the record was refused ({meta.get('checks')})")
    npz = rec_dir / meta["fragments_file"]
    if _sha256(npz) != meta["fragments_sha256"]:
        raise SystemExit(f"⛔ {npz}: its digest is not the one its record holds; it was written by another run")
    with np.load(npz) as z:
        return {"meta": meta, **{k: z[k] for k in ("key", "kind", "status", "mass")}}


def _stamp(meta: dict) -> dict:
    return {k: meta.get(k) for k in ("arm", "qa_arm", "condition", "bam", "settings", "message_policy", "em_seed",
                                     "assignment_mode", "revision")}


def pair_records(a: dict, b: dict) -> dict:
    """Two records of one condition joined by read name: the forecast sets by source, both arms' scores on the
    fragments both forecast, each arm's on the fragments it alone forecasts, and where the wrong-source mass lies.
    Refused (``SystemExit``) unless both were accepted, they read one input, every read name is held once, the two
    hold the same read names, and each name has one source."""
    ma, mb = a["meta"], b["meta"]
    checks = {
        "a_accepted": bool(ma.get("accepted")),
        "b_accepted": bool(mb.get("accepted")),
        "same_condition": ma["condition"] == mb["condition"],
        "same_input": ma["bam"] == mb["bam"],
        "a_read_names_unique": bool(np.unique(a["key"]).size == a["key"].size),
        "b_read_names_unique": bool(np.unique(b["key"]).size == b["key"].size),
    }
    ka, kb = a["key"], b["key"]
    oa, ob = np.argsort(ka, kind="stable"), np.argsort(kb, kind="stable")
    checks["same_read_names"] = bool(ka.size == kb.size and np.array_equal(ka[oa], kb[ob]))
    ib = ob[np.searchsorted(kb[ob], ka)] if checks["same_read_names"] else None
    checks["same_source"] = bool(ib is not None and np.array_equal(a["kind"], b["kind"][ib]))
    failed = [k for k, v in checks.items() if not v]
    if failed:
        raise SystemExit(f"⛔ {ma['condition']}: the pair is refused: {', '.join(failed)}")

    kind = a["kind"].astype(np.int64)
    sa, sb = a["status"], b["status"][ib]
    xa, xb = a["mass"], b["mass"][ib]
    fa, fb = sa == FORECAST, sb == FORECAST

    def by_kind(sel):
        return {k: int(np.sum(sel & (kind == i))) for i, k in enumerate(KINDS)}

    def by_other_status(sel, other):
        return {s: by_kind(sel & (other == i)) for i, s in enumerate(STATUSES) if np.any(sel & (other == i))}

    both = fa & fb
    sets = {
        "both": by_kind(both),
        "a_only": by_kind(fa & ~fb),
        "b_only": by_kind(~fa & fb),
        "neither": by_kind(~fa & ~fb),
        "a_only_by_b_status": by_other_status(fa & ~fb, sb),
        "b_only_by_a_status": by_other_status(~fa & fb, sa),
        "neither_by_status": {f"{STATUSES[i]}__{STATUSES[j]}": by_kind(~fa & ~fb & (sa == i) & (sb == j))
                              for i in range(len(STATUSES)) for j in range(len(STATUSES))
                              if np.any(~fa & ~fb & (sa == i) & (sb == j))},
    }
    common = {side: scores(x[both], kind[both])[0] for side, x in (("a", xa), ("b", xb))}
    moved = both & (p_gdna(xa) != p_gdna(xb))
    row = {
        "condition": ma["condition"],
        "stratum": list(stratum(ma["condition"])),
        "scope": scope(stratum(ma["condition"])),
        "zero_gdna": is_zero_gdna(ma["condition"]),
        "arms": {"a": ma["arm"], "b": mb["arm"]},
        "stamps": {"a": _stamp(ma), "b": _stamp(mb)},
        "varies": _varies(ma, mb),
        "checks": checks,
        "sets": sets,
        "moved": by_kind(moved),
        "common": common,
        "delta": None if common["a"] is None else {k: common["b"][k] - common["a"][k] for k in _SUMMED},
        "own_only": {"a": scores(xa[fa & ~fb], kind[fa & ~fb])[0], "b": scores(xb[~fa & fb], kind[~fa & fb])[0]},
        "where_they_lie": {side: where_they_lie(p_gdna(x[both]), kind[both]) for side, x in (("a", xa), ("b", xb))},
        "closure": {"a": ma["closure"], "b": mb["closure"]},
        "not_in_bam": ma["not_in_bam"],
    }
    return row


def _varies(ma: dict, mb: dict) -> dict:
    """What differs between the two arms' stamps: the settings either sets (the other's value or the shipped one),
    the ``quant_accuracy.py`` arm, and the revision."""
    sa, sb = ma.get("settings") or {}, mb.get("settings") or {}
    out = {k: [sa.get(k), sb.get(k)] for k in sorted(set(sa) | set(sb)) if sa.get(k) != sb.get(k)}
    for k in ("qa_arm", "revision"):
        if ma.get(k) != mb.get(k):
            out[k] = [ma.get(k), mb.get(k)]
    return out


def recorded_conditions(rec_dir: Path) -> list[str]:
    return sorted(p.stem for p in rec_dir.glob("*.json") if not p.name.startswith("."))


def pair_dirs(dir_a: Path, dir_b: Path, conditions=None) -> list[dict]:
    """Every condition both directories recorded, or ``conditions``; refused when the two hold different conditions
    and none are named, since a pair missing a condition must not be read as a whole panel."""
    ca, cb = recorded_conditions(dir_a), recorded_conditions(dir_b)
    if conditions:
        missing = [c for c in conditions if c not in ca or c not in cb]
        if missing:
            raise SystemExit(f"⛔ not recorded in both: {missing}")
        names = list(conditions)
    else:
        if ca != cb:
            raise SystemExit(f"⛔ the two records hold different conditions (only A: {sorted(set(ca) - set(cb))}; "
                             f"only B: {sorted(set(cb) - set(ca))}); name them with --conditions")
        names = ca
    rows = []
    for c in names:
        row = pair_records(load_record(dir_a, c), load_record(dir_b, c))
        rows.append(row)
        delta = "—" if row["delta"] is None else f"{row['delta']['brier']:+,.1f}"
        print(f"  {c}: both {sum(row['sets']['both'].values()):,}; Brier {_fmt(row['common']['a'], 'brier')} → "
              f"{_fmt(row['common']['b'], 'brier')} (Δ {delta})", flush=True)
    return rows


# ── the report: per stratum and rung, never pooled ───────────────────────────────────────────────


def stratum_rung(condition: str) -> tuple[str, str, str]:
    """The group a condition is summed in: its stratum and its rung (``g00`` the zero-gDNA control, or contaminated)."""
    return (*stratum(condition), "g00" if is_zero_gdna(condition) else "contaminated")


def _add(dicts):
    present = [d for d in dicts if d]
    return {k: sum(d[k] for d in present) for k in _SUMMED} if present else None


def _sum_group(rows: list[dict]) -> dict:
    """The rows of one stratum and rung, summed. Refused when they span two: a total over strata mostly reports the
    deferred one (``never-pool-the-strata``)."""
    spans = {stratum_rung(r["condition"]) for r in rows}
    if len(spans) != 1:
        raise SystemExit(f"⛔ a total would pool {sorted(spans)}: scores are summed within one stratum and rung only")
    out = {"n_conditions": len(rows),
           "sets": {s: {k: sum(r["sets"][s][k] for r in rows) for k in KINDS}
                    for s in ("both", "a_only", "b_only", "neither")},
           "moved": {k: sum(r["moved"][k] for r in rows) for k in KINDS},
           "common": {side: _add(r["common"][side] for r in rows) for side in ("a", "b")},
           "lie": {}}
    for side in ("a", "b"):
        lie = {}
        for k in KINDS:
            w = [r["where_they_lie"][side][k] for r in rows if k in r["where_they_lie"][side]]
            if w:
                grid = np.sum([x["wrong_mass_grid"] for x in w], axis=0)
                lie[k] = {"wrong_mass": sum(x["wrong_mass"] for x in w),
                          "at_certainty": sum(x["wrong_mass_at_certainty"] for x in w), "grid": grid}
        out["lie"][side] = lie
    return out


def totals(rows: list[dict], key=stratum_rung) -> dict:
    """Each group's total, keyed by ``key(condition)``, every group held to one stratum and rung (:func:`_sum_group`)."""
    groups: dict = {}
    for r in rows:
        groups.setdefault(key(r["condition"]), []).append(r)
    return {k: _sum_group(g) for k, g in groups.items()}


def _grid_quartiles(grid) -> list:
    g = np.asarray(grid, np.float64)
    if g.sum() <= 0:
        return [None] * len(QUARTILES)
    c = np.cumsum(g)
    return [float(np.searchsorted(c, q * c[-1]) / GRID) for q in QUARTILES]


def _fmt(d, k) -> str:
    return "—" if d is None else f"{d[k]:,.0f}"


def _q3(q) -> str:
    return "—" if q[0] is None else " / ".join(f"{v:.3f}" for v in q)


def _group_label(g) -> str:
    deferred = " ⛔ DEFERRED" if scope(g[:2]) == "DEFERRED" else ""
    return f"{g[0]} × {g[1]}{deferred} · {g[2]}"


def _ordered_groups(rows) -> list:
    conds = [r["condition"] for r in rows]
    present = {stratum_rung(c) for c in conds}
    return [(*st, rung) for st in strata(conds) for rung in ("g00", "contaminated") if (*st, rung) in present]


def tables(rows: list[dict]) -> list[tuple[str, list[str], list[list[str]]]]:
    """The report's tables, ``(title, header, rows)``, from a pair file's rows: per stratum and rung, then per
    condition. One source for the terminal and the markdown."""
    tot = totals(rows)
    groups = _ordered_groups(rows)
    a, b = rows[0]["arms"]["a"], rows[0]["arms"]["b"]
    out = []
    sets_rows = []
    for g in groups:
        s = tot[g]["sets"]
        rna = {n: s[n]["annotated"] + s[n]["nascent"] for n in s}
        sets_rows.append([_group_label(g), *(f"{s[n]['gdna']:,}" for n in ("both", "a_only", "b_only", "neither")),
                          *(f"{rna[n]:,}" for n in ("both", "a_only", "b_only", "neither"))])
    out.append(("The forecast sets (fragments, by source)",
                ["stratum · rung", "gDNA both", f"gDNA {a} only", f"gDNA {b} only", "gDNA neither", "RNA both",
                 f"RNA {a} only", f"RNA {b} only", "RNA neither"], sets_rows))
    score_rows = []
    for g in groups:
        ca, cb = tot[g]["common"]["a"], tot[g]["common"]["b"]
        if ca is None:
            score_rows.append([_group_label(g)] + ["—"] * 9)
            continue
        ratio = f"{cb['brier'] / ca['brier']:.3f}" if ca["brier"] > 0 else "—"
        score_rows.append([_group_label(g), _fmt(ca, "brier"), _fmt(cb, "brier"), f"{cb['brier'] - ca['brier']:+,.0f}",
                           ratio, _fmt(ca, "mcb"), _fmt(cb, "mcb"), f"{_fmt(ca, 'over')} / {_fmt(ca, 'under')}",
                           f"{_fmt(cb, 'over')} / {_fmt(cb, 'under')}",
                           f"{sum(tot[g]['moved'].values()):,}"])
    out.append(("On the fragments both forecast: Brier, its miscalibration (MCB) and its direction",
                ["stratum · rung", f"{a} Brier", f"{b} Brier", "Δ", "ratio", f"{a} MCB", f"{b} MCB",
                 f"{a} OVER / UNDER", f"{b} OVER / UNDER", "moved"], score_rows))
    src_rows = []
    for g in groups:
        ca, cb = tot[g]["common"]["a"], tot[g]["common"]["b"]
        src_rows.append([_group_label(g), *(f"{_fmt(ca, f'brier_{k}')} → {_fmt(cb, f'brier_{k}')}" for k in KINDS),
                         f"{_fmt(ca, 'linear')} → {_fmt(cb, 'linear')}", f"{_fmt(ca, 'brier3')} → {_fmt(cb, 'brier3')}"])
    out.append(("By source: Brier, and the linear and three-way scores (A → B)",
                ["stratum · rung", "gDNA Brier", "annotated Brier", "nascent Brier", "linear", "three-way Brier"],
                src_rows))
    lie_rows = []
    for g in groups:
        for side, name in (("a", a), ("b", b)):
            lie = tot[g]["lie"][side]
            gd = lie.get("gdna")
            rna = [lie[k] for k in ("annotated", "nascent") if k in lie]
            r_mass = sum(x["wrong_mass"] for x in rna)
            r_cert = sum(x["at_certainty"] for x in rna)
            g_cells = ([f"{gd['wrong_mass']:,.0f}", _q3(_grid_quartiles(gd["grid"])),
                        f"{(gd['at_certainty'] / gd['wrong_mass'] if gd['wrong_mass'] else 0):.1%}"]
                       if gd else ["—", "—", "—"])
            r_cells = ([f"{r_mass:,.0f}", _q3(_grid_quartiles(np.sum([x["grid"] for x in rna], axis=0))),
                        f"{(r_cert / r_mass if r_mass else 0):.1%}"] if rna else ["—", "—", "—"])
            lie_rows.append([_group_label(g), name, *g_cells, *r_cells])
    out.append(("Where the wrong-source mass lies (on the fragments both forecast): quartiles of P(gDNA) weighted by "
                "it, and its share at certainty",
                ["stratum · rung", "arm", "gDNA→RNA", "its quartiles of P(gDNA)", "at certainty", "RNA→gDNA",
                 "its quartiles of P(gDNA)", "at certainty"], lie_rows))
    cond_rows = []
    for r in sorted(rows, key=lambda r: (_ordered_groups(rows).index(stratum_rung(r["condition"])), r["condition"])):
        ca, cb = r["common"]["a"], r["common"]["b"]
        worst = max((v["diff"] / v["bound"] if v["bound"] > 0 else (0.0 if v["diff"] == 0 else float("inf"))
                     for side in ("a", "b") for v in r["closure"][side].values()), default=float("nan"))
        cond_rows.append([r["condition"], f"{sum(r['sets']['both'].values()):,}",
                          f"{sum(r['sets']['a_only'].values()):,}", f"{sum(r['sets']['b_only'].values()):,}",
                          _fmt(ca, "brier"), _fmt(cb, "brier"),
                          "—" if r["delta"] is None else f"{r['delta']['brier']:+,.1f}",
                          _fmt(ca, "mcb"), _fmt(cb, "mcb"), f"{worst:.3f}"])
    out.append(("Per condition (closure: the worst |class sum − pool| / its rounding bound over both arms; ≤ 1 holds)",
                ["condition", "both", f"{a} only", f"{b} only", f"{a} Brier", f"{b} Brier", "Δ", f"{a} MCB",
                 f"{b} MCB", "closure"], cond_rows))
    return out


def _header_lines(rows) -> list[str]:
    st = rows[0]["stamps"]
    lines = [f"A = {rows[0]['arms']['a']}: --set {st['a']['settings']} · revision {st['a']['revision']}",
             f"B = {rows[0]['arms']['b']}: --set {st['b']['settings']} · revision {st['b']['revision']}",
             f"input {st['a']['bam']} · varies {rows[0]['varies'] or 'nothing'} · {len(rows)} conditions"]
    return lines


def check_pair_file(rows: list[dict]) -> None:
    """Refuse a pair file that mixes pairs: every row joins the same two arms on one input."""
    sig = {(r["arms"]["a"], r["arms"]["b"], r["stamps"]["a"]["bam"],
            json.dumps(r["varies"], sort_keys=True)) for r in rows}
    if len(sig) != 1:
        raise SystemExit(f"⛔ the rows join different arms or inputs: {sorted(sig)}")


def report(rows: list[dict], markdown: Path | None = None) -> None:
    """Print the tables; with ``markdown``, also write them there."""
    check_pair_file(rows)
    print()
    for line in _header_lines(rows):
        print(f"  {line}")
    for title, head, body in tables(rows):
        widths = [max(len(str(x)) for x in col) for col in zip(head, *body)]
        print(f"\n  {title}")
        print("    " + "  ".join(h.rjust(w) if i else h.ljust(w) for i, (h, w) in enumerate(zip(head, widths))))
        for row in body:
            print("    " + "  ".join(c.rjust(w) if i else c.ljust(w) for i, (c, w) in enumerate(zip(row, widths))))
    if markdown is None:
        return
    md = ["# Paired proper score", "", *(f"- {line}" for line in _header_lines(rows)), "",
          "Scored on the fragments both arms forecast, joined by read name; each stratum and rung summed over its "
          "conditions, never across strata. Lower Brier is better; MCB is the part a monotone recalibration "
          "removes.", ""]
    for title, head, body in tables(rows):
        md += [f"### {title}", "", "| " + " | ".join(head) + " |", "|---|" + "---:|" * (len(head) - 1),
               *("| " + " | ".join(row) + " |" for row in body), ""]
    markdown.parent.mkdir(parents=True, exist_ok=True)
    markdown.write_text("\n".join(md))
    print(f"\n  ⭐ markdown -> {markdown}")


# ── the self-test ────────────────────────────────────────────────────────────────────────────────


def _mass_of(p_gdna_values, annotated=None) -> np.ndarray:
    """Rows with ``P(gDNA) = p`` on the locus gDNA class, the rest annotated, or split as given."""
    p = np.asarray(p_gdna_values, np.float64)
    m = np.zeros((p.size, len(EM_CLASSES)))
    m[:, GDNA] = p
    m[:, ANN] = 1.0 - p if annotated is None else annotated
    if annotated is not None:
        m[:, NAS] = 1.0 - p - annotated
    return m


def _check_scores(check) -> None:
    """The scores against hand-worked values: the proper score prefers the calibrated forecast, the improper one the
    certain forecast; CORP's blocks, OVER, UNDER and the three-way score."""
    close = np.isclose
    kind = np.array([0] * 50 + [1] * 50)
    o, _ = scores(_mass_of(np.full(100, 0.5)), kind)
    check("50/50 at p = 0.5: linear 50, Brier 25, MCB 0, DSC 0, UNC 25",
          close(o["linear"], 50) and close(o["brier"], 25) and close(o["mcb"], 0) and close(o["dsc"], 0)
          and close(o["unc"], 25))
    o, _ = scores(_mass_of(np.zeros(100)), kind)
    check("50/50 called RNA: linear 50, Brier 50, MCB 25, UNDER 50, OVER 0",
          close(o["linear"], 50) and close(o["brier"], 50) and close(o["mcb"], 25) and close(o["under"], 50)
          and close(o["over"], 0))
    kind = np.array([0] * 40 + [1] * 60)
    cal, _ = scores(_mass_of(np.full(100, 0.4)), kind)
    cert, _ = scores(_mass_of(np.zeros(100)), kind)
    check("40/60: the calibrated forecast costs linear 48, Brier 24; the certain one linear 40, Brier 40",
          close(cal["linear"], 48) and close(cal["brier"], 24) and close(cert["linear"], 40)
          and close(cert["brier"], 40))
    check("the linear score prefers the certain forecast and the Brier the calibrated one",
          cert["linear"] < cal["linear"] and cal["brier"] < cert["brier"])
    o, _ = scores(_mass_of(np.r_[np.ones(30), np.zeros(70)]), np.array([0] * 30 + [1] * 70))
    check("a perfect forecast: Brier 0, MCB 0, DSC = UNC",
          o["brier"] == 0 and close(o["mcb"], 0) and close(o["dsc"], o["unc"]))
    p = np.r_[np.full(100, 0.9), np.full(100, 0.1)]
    kind = np.r_[np.zeros(60, int), np.ones(40, int), np.zeros(40, int), np.ones(60, int)]
    o, blocks = scores(_mass_of(p), kind)
    want_bs = 2 * (60 * 0.1**2 + 40 * 0.9**2)
    want_br = 2 * (60 * 0.4**2 + 40 * 0.6**2)
    check("overconfident both ways: OVER 30, UNDER 30, net 0, MCB = Brier − the refit's, two blocks",
          close(o["over"], 30) and close(o["under"], 30) and abs(o["net_pool_error"]) < 1e-9
          and close(o["brier"], want_bs) and close(o["mcb"], want_bs - want_br) and len(blocks["fit"]) == 2)
    p = np.r_[np.full(10, 0.2), np.full(10, 0.8)]
    kind = np.r_[np.zeros(8, int), np.ones(2, int), np.zeros(2, int), np.ones(8, int)]
    o, blocks = scores(_mass_of(p), kind)
    check("a violation the fit pools: one block at 0.5, DSC 0, OVER 3, UNDER 3",
          len(blocks["fit"]) == 1 and close(blocks["fit"][0], 0.5) and close(o["dsc"], 0) and close(o["over"], 3)
          and close(o["under"], 3))
    o, _ = scores(np.array([[0.3, 0.5, 0.2, 0.0]]), np.array([2]))
    check("three-way by hand: a nascent fragment at (gDNA 0.2, annotated 0.3, nascent 0.5)",
          close(o["brier3"], 0.2**2 + 0.3**2 + 0.5**2) and close(o["brier"], 0.2**2))
    o, _ = scores(np.array([[0.0, 0.0, 0.25, 0.75]]), np.array([0]))
    check("locus and intergenic gDNA are one source", o["brier"] == 0)
    check("a row with no forecast is refused", refuses(scores, np.zeros((1, 4)), np.array([0]), exc=ValueError))
    check("no fragments score nothing", scores(np.zeros((0, 4)), np.zeros(0, int)) == (None, None))


def _check_bounds(check) -> None:
    """The bounds' closed forms, and a solver-like normalisation never leaving its bound."""
    check("γ(1) = u / (1 − u)", float(gamma(1)) == U64 / (1 - U64))
    check("the conservation bound at 4 candidates is γ(14)",
          float(conservation_bound(4)) == 14 * U64 / (1 - 14 * U64))
    w = float(winner_bound(np.array([0.5]), np.array([2]))[0])
    check("the winner bound's closed form",
          w == 0.5 * (1 + U32) / (1 - float(gamma(2))) * (1 + float(gamma(3))))
    rng = np.random.default_rng(0)
    worst = 0.0
    for n in (1, 2, 3, 7, 50, 400):
        for _ in range(300):
            ll = rng.normal(0, 30, n)
            e = np.exp(ll - ll.max())
            s = 0.0
            for x in e:
                s += x
            post = e * (1.0 / s)
            row = np.zeros(4)
            for x, c in zip(post, rng.integers(0, 4, n)):
                row[c] += x
            worst = max(worst, abs(row.sum() - 1.0) / float(conservation_bound(n)))
    check("a solver-like normalisation stays within its conservation bound", worst <= 1.0)


def _toy_run():
    """A run and its read pairs that close: five read pairs, four of them fragments. Fragment 0 is an EM unit, 1 a
    deterministic annotated fragment, 2 the intergenic fast path, 3 withheld; the fifth read pair has no eligible
    record. Transcript 0 is annotated, 1 a synthetic span; the unit's locus is 7."""
    from rigel.annotate import AF_ANNOTATED, AF_GDNA_INTERGENIC, AF_UNRESOLVED
    from rigel.buffer import FRAG_UNAMBIG

    unit = np.zeros(len(EM_CLASSES))
    unit[[ANN, NAS, GDNA]] = (0.5, 0.125, 0.375)
    run = Run(
        em_fid=np.array([0]), em_lid=np.array([7]), em_mass=unit[None, :].copy(), em_wtid=np.array([0], np.int32),
        em_wpost=np.array([0.5], np.float32), em_ncand=np.array([3]),
        t_fid=np.array([1, 0]), t_tid=np.array([0, 0]), t_flags=np.array([AF_ANNOTATED, AF_ANNOTATED]),
        t_post=np.array([1.0, 0.5]), t_class=np.array([FRAG_UNAMBIG, 2], np.int8), t_ncand=np.array([1, 3]),
        is_synthetic=np.array([False, True]),
        locus_totals={7: (0.375, 0.625)}, intergenic_em=0.0,
        ledger={"n_read_pairs": 5, "n_excluded": 1, "n_fragments": 3, "n_intergenic": 1},
        pools={"annotated": 1.5, "nascent": 0.125, "gdna": 0.375, "intergenic": 1.0}, n_extra_terms=10,
    )
    names = ["t0:0-100:f:0", "t0:0-100:f:1", "gdna:chr1:0-100:f:2", "nrna_t0:0-100:f:3", "gdna:chr1:5-100:f:4"]
    # (name, flag, ZF): the unit's two mates, the deterministic one, the fast path, the withheld one, then an
    # unmapped pair
    records = [(names[0], 0x41, AF_ANNOTATED), (names[0], 0x81, AF_ANNOTATED), (names[1], 0x41, AF_ANNOTATED),
               (names[2], 0x41, AF_GDNA_INTERGENIC), (names[3], 0x41, AF_UNRESOLVED), (names[4], 0x4, None),
               (names[4], 0x4, None)]
    return run, records


def _check_record(check) -> None:
    """The record composes the toy, closes on its pools, and fires on each planted defect: two classes swapped, a
    fragment id joined to the wrong read name, a count the scan disagrees with, a unit that does not conserve."""
    from rigel.annotate import AF_ANNOTATED, AF_GDNA_EM

    run, records = _toy_run()
    pairs = read_pairs(records, skip_duplicates=True)
    check("the read-pair rule: five read pairs, the unmapped one with no fragment id",
          pairs.key.size == 5 and pairs.eligible.tolist() == [True] * 4 + [False]
          and pairs.kind.tolist() == [1, 1, 0, 2, 0])
    arrays, rec = build_record(run, pairs)
    check("the toy record is accepted and closes on every pool", rec["accepted"] and all(
        v["diff"] <= v["bound"] for v in rec["closure"].values()))
    check("each read pair gets its status: forecast ×3, withheld, excluded",
          arrays["status"].tolist() == [FORECAST, FORECAST, FORECAST, WITHHELD, EXCLUDED])
    check("the fast path is mass 1 intergenic, the deterministic fragment mass 1 annotated",
          arrays["mass"][2, IG] == 1.0 and arrays["mass"][1, ANN] == 1.0)

    swapped = dataclasses.replace(run, em_mass=run.em_mass[:, [GDNA, NAS, ANN, IG]])
    _, rec = build_record(swapped, pairs)
    check("PLANTED two classes swapped (annotated ↔ gDNA): the closure fires", rec["checks"]["closure_over_bound"] > 0
          and not rec["accepted"])
    check("PLANTED two classes swapped: the winner check fires", rec["checks"]["winner_above_its_class"] > 0)

    shifted = [(n, f | (0x4 if i == 2 else 0), z) for i, (n, f, z) in enumerate(records)]
    _, rec = build_record(run, read_pairs(shifted, skip_duplicates=True))
    check("PLANTED a fragment id joined to the wrong read name (an id rule off by one): the ZF join fires",
          rec["checks"].get("zf_not_the_tables_flag", 0) > 0 or rec["checks"].get("fragment_ids_out_of_range", 0) > 0)
    check("PLANTED the same: the record is refused", not rec["accepted"])

    uncollated = [(records[0][0] if i == 3 else n, f, z) for i, (n, f, z) in enumerate(records)]
    _, rec = build_record(run, read_pairs(uncollated, skip_duplicates=True))
    check("PLANTED one read name in two runs (an uncollated input): the repeated name fires",
          rec["checks"]["read_names_repeated"] == 1 and not rec["accepted"])

    fewer = dataclasses.replace(run, ledger={**run.ledger, "n_intergenic": 0})
    check("PLANTED a count the scan disagrees with: the ledger fires",
          build_record(fewer, pairs)[1]["checks"]["ledger_mismatches"] == 1)
    leaky = dataclasses.replace(run, em_mass=run.em_mass * (1 + 1e-9), pools={
        k: v * (1 + 1e-9) if k != "intergenic" else v for k, v in run.pools.items()})
    check("PLANTED a unit summing to 1 + 1e-9: conservation fires",
          build_record(leaky, pairs)[1]["checks"]["conservation_over_bound"] == 1)
    off = dataclasses.replace(run, pools={**run.pools, "annotated": 1.5 + 1e-9})
    check("PLANTED a pool one part in 1e9 away: the closure fires",
          build_record(off, pairs)[1]["checks"]["closure_over_bound"] == 1)
    gdna_won = dataclasses.replace(run, t_flags=np.array([AF_ANNOTATED, AF_GDNA_EM]))
    check("PLANTED a table flag the writer did not stamp: the ZF join fires",
          build_record(gdna_won, pairs)[1]["checks"]["zf_not_the_tables_flag"] == 1)


def _toy_records():
    """Two in-memory records of one condition over six read names (gDNA, annotated, nascent, twice)."""
    rng = np.random.default_rng(1)
    names = [f"gdna:chr1:{i}-{i + 100}:f:{i}" for i in range(2)] + [f"t{i}:0-100:f:{i}" for i in range(2)] + \
            [f"nrna_t{i}:0-100:f:{i}" for i in range(2)]
    key = np.array([name_key(n) for n in names], np.uint64)
    kind = np.array([0, 0, 1, 1, 2, 2], np.int8)

    def rec(arm):
        m = rng.dirichlet(np.ones(len(EM_CLASSES)), size=6)
        meta = {"accepted": True, "condition": "gdna_g50_ss_0.99_nrna_mid_capture_on", "bam": ORACLE_BAM,
                "arm": arm, "settings": {"em.assignment_mode": "fractional"}, "revision": "r", "closure": {},
                "not_in_bam": {k: 0 for k in KINDS}}
        return {"meta": meta, "key": key.copy(), "kind": kind.copy(), "status": np.zeros(6, np.int8), "mass": m}

    return rec("a"), rec("b")


def _check_pair(check) -> None:
    """The join: a mis-joined read name, a foreign one and a repeated one are refused; a fragment one arm forecasts
    and the other does not is counted, never scored; a shuffled record joins to the same scores."""
    a, b = _toy_records()
    row = pair_records(a, b)
    check("two records of one condition pair, every fragment in both", sum(row["sets"]["both"].values()) == 6)
    perm = np.array([1, 0, 3, 2, 5, 4])[::-1]
    shuffled = {**b, "key": b["key"][perm], "kind": b["kind"][perm], "status": b["status"][perm], "mass": b["mass"][perm]}
    check("a record in another order joins by name to the same scores",
          pair_records(a, shuffled)["common"]["b"]["brier"] == row["common"]["b"]["brier"])
    misjoined = {**b, "key": b["key"][[2, 1, 0, 3, 4, 5]]}
    check("PLANTED a mis-joined read name (two names exchanged across sources): refused",
          refuses(pair_records, a, misjoined))
    foreign = {**b, "key": np.r_[b["key"][:5], np.uint64(name_key("elsewhere"))]}
    check("PLANTED a foreign read name: refused", refuses(pair_records, a, foreign))
    repeated = {**b, "key": np.r_[b["key"][:5], b["key"][0]]}
    check("PLANTED a read name held twice: refused", refuses(pair_records, a, repeated))
    other = {**b, "meta": {**b["meta"], "bam": "rustar_recommended.bam"}}
    check("PLANTED two inputs: refused", refuses(pair_records, a, other))
    one_only = {**b, "status": np.r_[b["status"][:5], np.int8(UNASSIGNED)]}
    row1 = pair_records(a, one_only)
    check("a fragment only one arm forecasts is counted in its own set, by the other's status",
          row1["sets"]["a_only"]["nascent"] == 1 and sum(row1["sets"]["both"].values()) == 5
          and row1["sets"]["a_only_by_b_status"]["unassigned"]["nascent"] == 1)
    check("... and is not scored on the common set", row1["common"]["a"]["n_forecast"] == 5)
    same = pair_records(a, {**a, "meta": {**a["meta"], "arm": "a2"}})
    check("an arm against itself: Δ Brier exactly 0, nothing moved",
          same["delta"]["brier"] == 0.0 and sum(same["moved"].values()) == 0)
    arm = pair_records({**a, "meta": {**a["meta"], "qa_arm": "base"}}, {**b, "meta": {**b["meta"], "qa_arm": "no_capture"}})
    check("two records under one setting and two arms: the pair names the arm, and only it, as what varies",
          arm["varies"] == {"qa_arm": ["base", "no_capture"]})


def _check_strata(check) -> None:
    """The totals: one per stratum and rung, the deferred stratum apart, and a grouping that would pool strata is
    refused."""
    a, b = _toy_records()
    rows = []
    for cond in ("gdna_g50_ss_0.99_nrna_mid_capture_on", "gdna_g50_ss_0.50_nrna_mid_capture_on",
                 "gdna_g00_ss_0.99_nrna_mid_capture_on", "gdna_g25_ss_0.99_nrna_mid_capture_on"):
        ra = {**a, "meta": {**a["meta"], "condition": cond}}
        rb = {**b, "meta": {**b["meta"], "condition": cond}}
        rows.append(pair_records(ra, rb))
    tot = totals(rows)
    check("one total per stratum and rung: stranded ON contaminated (2), its g00 (1), unstranded ON (1)",
          sorted((k, v["n_conditions"]) for k, v in tot.items()) == sorted([
              (("stranded", "capture ON", "contaminated"), 2), (("stranded", "capture ON", "g00"), 1),
              (("unstranded", "capture ON", "contaminated"), 1)]))
    check("a stratum's total is the sum of its conditions",
          np.isclose(tot[("stranded", "capture ON", "contaminated")]["common"]["a"]["brier"],
                     rows[0]["common"]["a"]["brier"] + rows[3]["common"]["a"]["brier"]))
    check("PLANTED strata pooled (one key for every condition): refused", refuses(totals, rows, lambda c: "all"))
    check("PLANTED the zero control pooled with its stratum (rungs merged): refused",
          refuses(totals, rows, lambda c: stratum(c)))
    labels = [t[2] for t in tables(rows) if t[0].startswith("The forecast sets")][0]
    check("the deferred stratum is marked in every table it appears in",
          any("DEFERRED" in r[0] for r in labels) and sum("DEFERRED" in r[0] for r in labels) == 1)
    mixed = [rows[0], {**rows[1], "arms": {"a": "x", "b": "y"}}]
    check("PLANTED a pair file mixing two pairs: refused", refuses(check_pair_file, mixed))


def _check_protocol(check) -> None:
    """The protocol is applied, and a --set that undoes it is refused."""
    cfg = protocol_config(["calibration.message_policy=silent"])
    check("the protocol pins the scan and the assignment, and keeps the arm's own setting",
          cfg.scan.total_threads == 1 and cfg.em.assignment_mode == "fractional"
          and cfg.calibration.message_policy == "silent")
    check("PLANTED --set em.assignment_mode=sample: refused", refuses(protocol_config, ["em.assignment_mode=sample"]))
    check("PLANTED --set scan.total_threads=4: refused", refuses(protocol_config, ["scan.total_threads=4"]))


def self_test() -> int:
    """Every comparator and refusal, planted and required to fire; no I/O."""
    checks: list[tuple[str, bool]] = []

    def check(name: str, ok) -> None:
        checks.append((name, bool(ok)))

    for part in (_check_scores, _check_bounds, _check_record, _check_pair, _check_strata, _check_protocol):
        part(check)
    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"  {'✅' if ok else '⛔'} {name:<{width}}")
    n_ok = sum(ok for _, ok in checks)
    print(f"\n  {n_ok}/{len(checks)} self-test gates fired")
    return 0 if n_ok == len(checks) else 1


# ── cli ──────────────────────────────────────────────────────────────────────────────────────────


def _record(args) -> int:
    from rigel.index import TranscriptIndex

    cfg = protocol_config(args.settings)
    names = panel_conditions(args.suite, args.conditions)
    label = args.label or args.out.name
    qa_arm = args.arm or "base"
    if args.jobs > 1 and len(names) > 1:
        run_shards(Path(__file__).resolve(), names, args.jobs, args.work_dir / f".{label}_shards",
                   lambda _out, shard: ["--record", "--arm", qa_arm, "--suite", str(args.suite),
                                        "--index", str(args.index), "--bam", args.bam, "--out", str(args.out),
                                        "--label", label, "--work-dir", str(args.work_dir), "--conditions", *shard,
                                        *set_args(args.settings)])
        return 0
    index = TranscriptIndex.load(str(args.index))
    refused = []
    for name in names:
        print(f"  … {label}  {name}", flush=True)
        arrays, row = record_condition(label, args.suite, index, name, cfg, args.bam, args.work_dir, qa_arm)
        write_record(args.out, name, arrays, row)
        worst = max((v["diff"] / v["bound"] if v["bound"] > 0 else float(v["diff"] > 0))
                    for v in row["closure"].values()) if row["closure"] else float("nan")
        print(f"  {'⭐' if row['accepted'] else '⛔'} {name}: accepted {row['accepted']}; forecast "
              f"{sum(row['status']['forecast'].values()) if row['status'] else 0:,}; Brier "
              f"{_fmt(row['scores'], 'brier')}; closure worst diff/bound {worst:.3f}; "
              f"{row['seconds']:.0f} s", flush=True)
        if not row["accepted"]:
            refused.append(name)
    if refused:
        raise SystemExit(f"⛔ refused: {refused} (each record's checks are in its JSON)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--record", action="store_true", help="run one arm on each condition and write its records "
                                                            "to --out (a directory)")
    mode.add_argument("--pair", nargs=2, type=Path, metavar=("A", "B"),
                      help="join two record directories by read name and write the pair rows to --out (jsonl)")
    mode.add_argument("--report", type=Path, default=None, help="print the tables of a pair file")
    mode.add_argument("--self-test", action="store_true", help="plant every defect; no I/O")
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--bam", default=ORACLE_BAM,
                    help="the BAM each condition is read from, by name in its directory: the simulator's own "
                         "alignments by default, or an aligner's (panel.py align)")
    ap.add_argument("--out", type=Path, default=None,
                    help="--record: the arm's record directory; --pair: the pair file (jsonl)")
    ap.add_argument("--label", default=None, help="--record: the arm's name in the tables (default: --out's name)")
    ap.add_argument("--arm", choices=ARMS, default=None,
                    help="--record: the quant_accuracy.py arm to run under (default base; no_capture reads every "
                         "capture efficiency as 1)")
    ap.add_argument("--markdown", type=Path, default=None, help="--pair or --report: also write the tables here")
    ap.add_argument("--work-dir", type=Path, default=SCRATCH / "rigel_paired_score",
                    help="where each run's annotated BAM is written, read once, and deleted")
    ap.add_argument("--jobs", type=int, default=1, help="--record: conditions in parallel, one process each")
    add_set_option(ap, "--record: a config value for this arm, repeatable, after the protocol "
                       f"({', '.join(PROTOCOL)}); the same parser as every instrument's --set")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.report:
        report([json.loads(x) for x in args.report.read_text().splitlines() if x.strip()], args.markdown)
        return 0
    if args.out is None:
        raise SystemExit("⛔ --out is required")
    if args.record:
        return _record(args)
    if args.pair:
        if args.settings or args.arm:
            raise SystemExit("⛔ --set and --arm belong to --record: a pair reads the settings and the arm its records "
                             "were made under")
        rows = pair_dirs(*args.pair, args.conditions)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("".join(json.dumps(r) + "\n" for r in rows))
        print(f"  ⭐ {len(rows)} pairs -> {args.out}")
        report(rows, args.markdown)
        return 0
    raise SystemExit("⛔ one of --record, --pair, --report or --self-test")


if __name__ == "__main__":
    raise SystemExit(main())
