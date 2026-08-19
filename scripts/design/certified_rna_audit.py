#!/usr/bin/env python
"""⭐⭐⭐ IS THE CERTIFIED-RNA CHANNEL WIRED AT A TERMINUS BOUNDARY? — the audit, on a TA x TB abundance grid.

A **spliced** fragment cannot be gDNA, so ``boundary_spliced`` is the one observation in the tool that needs no
deconvolution: it is certified RNA. The solver's job is to deconvolve the **unspliced** population, and any
certified fragment that leaks into it is RNA handed to the gDNA solver.

⛔ **THE CASE THIS EXISTS FOR** (owner, 2026-08-05). The accumulator now creates an BOUNDARY at every
transcript start and end. **Those BOUNDARIES carry no splice junction** — nothing splices there — **yet spliced
fragments cross them**, because transcription of a longer overlapping transcript continues straight past.
Unless such a fragment is binned as ``boundary_spliced`` it lands in the unspliced pool and is deconvolved.

The rung is ``tes_readthrough``::

    TA+ (1,050, 2,000) (9,000,  9,100)      sj 2,000 -> 9,000
    TB+ (1,000, 2,000) (9,050, 11,000)      sj 2,000 -> 9,050

⭐ **BOUNDARY @9,100 is TA's TES and no sj touches it**, but TB's exon 2 runs to 11,000 — so a TB
fragment that used TB's sj and reaches past 9,100 crosses it CONTIGUOUSLY HAVING SPLICED ELSEWHERE.
⭐ **BOUNDARY @9,050 is TB's sj acceptor AND a plain contiguity boundary for TA**, whose exon 2 spans
9,000-9,100 unbroken. One boundary, sj flux for one transcript and an unspliced RNA crossing for another.

⛔⛔ **SCORE THE MASS, NOT ``f_g`` — AND THIS INSTRUMENT LEARNED THAT THE HARD WAY.** The solver's
``f_g`` is the gDNA fraction of the **unspliced** population; the oracle's per-object fraction is
spliced-INCLUSIVE. Comparing them directly reads a 1.5 %-of-mass defect as a half-unit error, which is
how the first run of this file mis-reported it. The columns that decide anything are ``Δgdna`` —
``|PRED gDNA − TRUE gDNA|`` in FRAGMENTS — and ``mass``, the object's own true total.

⚠ **Those two columns were MISSING from the printed table until 2026-08-17, so this paragraph's own rule
could not be followed off the output it describes**: every printed column was a fraction. The older
wording named ``TRUE gDNA`` / ``PRED gDNA`` columns, which this instrument does not have — `toy_harness`
publishes the pair as ``err`` and ``mass`` on each object row, and those are what is now printed.

**What is checked, per object, at every grid cell** — three things, because any one of them failing
silently produces the same symptom:

===  ====================================================================================
(a)  is the BANK populated?          ``spliced_count`` > 0 where the geometry says it must be
(b)  does it have a DIVISOR?         ``eff_rna`` > 0, since there is no ``eff_sj`` at a
                                     terminus BOUNDARY to price it against
(c)  is a PRECISION EMITTED?         the relay's own RNA measurement precision (``cm_p``/``cm_n``)
                                     — a bank with a divisor and no precision is inert
===  ====================================================================================

---

⛔⛔⛔ **THE VERDICT THAT USED TO BE HERE — "24 of 24 grid cells … (c) FAILS, ``cm_p = cm_n = 0``" — IS
WITHDRAWN, AND ITS NUMBERS ARE NOT REPLACED BY GUESSES (2026-08-17).** It cannot be reproduced by this
file and it is not honestly repairable into one, for three separate reasons:

* **The grid is not 24 cells.** The defaults are ``--ta 0 30 300 3000`` × ``--tb`` the same, with the
  silent ``(0, 0)`` corner skipped — **15** cells.
* **It named columns this instrument does not print.** ``TRUE gDNA`` / ``PRED gDNA`` per cell; see the
  MASS paragraph above.
* ⛔⛔ **It was measured BEFORE ``message_propagation`` was defaulted OFF, and check (c) is a RELAY
  quantity.** ``cm_p`` / ``cm_n`` are the relay's own RNA measurement precisions and exist only under
  ``RelayPolicy``. Reporting "(c) FAILS, precisions are 0" against a run in which the relay was never
  installed would be `TRAPS: an-ablation-that-never-ran`, so under the mute this file now reports (c)
  **UNANSWERABLE** and exits non-zero rather than passing OR failing.

⭐⭐ **WHAT A FRESH RUN SAYS, 2026-08-17, both policies, the full 15-cell default grid** — quoted from the
output, not remembered:

* ``--messages on`` (this file's default, and the only setting in which (c) has content): **rc 0**, and
  *"✅ every populated ``boundary_spliced`` bank has a divisor AND emits a precision."* A precision IS
  emitted — e.g. at TA = 3000 / TB = 3000, ``cm_p`` reads 0.1988 at @9,050 and 0.2026 at @9,100. ⚠
  ``cm_n`` is 0.0 at every row of that run, which is the NEGATIVE-strand channel on a two-plus-strand
  toy and is not evidence about the certified channel.
* ``--messages off`` (**the shipped config**): **rc 1**, and (c) is UNANSWERABLE on every populated bank
  in all 15 cells — ``cm_p`` / ``cm_n`` print as ``—`` because the relay published nothing.

⛔ **So the ONE claim the old verdict rested on — that the certified channel is inert — is now a claim
about the SHIPPED MUTE, not about the channel's wiring, and this file can no longer be used to make it.**
What survives unchanged is (a) and (b): ``boundary_spliced`` is populated at both terminus/contiguity
BOUNDARIES and carries ``eff_rna`` = 202.8 as its divisor on every populated row of both runs.

⭐ **@9,050 is still the diagnostic and @9,100 still its control, and the geometry is why.** Both carry
certified RNA; @9,100's flanks include the 1,900 bp exon [9,100, 11,000) which the relay can speak from.
@9,050's flanks are the two **50 bp** regions [9,000, 9,050) and [9,050, 9,100) — both below one mean
fragment length, so neither has a resolvable density (`TRAPS: density-below-one-fragment-length`) and a
muted relay has nothing to offer. ⚠ Under ``--messages off`` every one of those objects reads ``pred_fg``
0.4947–0.5012 — the uninformative reference — against ``true_fg`` 0.0006–0.0032; under ``--messages on``
the same objects move (0.0394 at TA = 3000 / TB = 30, 0.2603 at TA = 3000 / TB = 3000). ⛔ **The object's
own certified-RNA fragments are not what moves it in either case** — ψ has no spliced term at all — so
the answer is still set by whether a NEIGHBOUR happens to be informative.

⭐ **The grid is the experiment, not a sweep for its own sake.** At @9,100 the certified channel is TB's
alone while the unspliced crossing there is gDNA + TB, so the TA/TB ratio moves the two independently and a
single abundance pair tests one corner. TA also carries the *other* sj into @9,050, so the two
transcripts stress different boundaries.

⚠ **A CORRECT deferral to keep in view — and it is no longer a LOSS.** TA's and TB's sj share the donor
at 2,000 and differ only in acceptor (9,000 vs 9,050), so a fragment whose UNSEQUENCED mate gap could
hold either implies two different ``L`` values and is DEFERRED by pass one rather than deposited — the
owner's ruling, in `tests/native/_accumulator_reference.py`. ⛔⛔ **`run_toy` now runs the SECOND PASS on
the whole before it returns, so those fragments are PLACED and the old "13.5 % of the certified channel
at @9,100 is missing" figure describes an undrained tally this file no longer produces.** Measured
2026-08-17 at TA = 3000 / TB = 30: **36 held → 36 placed by the second pass, 0 still deferred**. ⚠ The
line printed per cell is the DRAIN's own ``offered``/``deposited`` pair for exactly that reason —
``ScanQC.deferred_undetermined_gap`` alone reads 0 on a drained payload and a reader takes that as "no
fragment was ever ambiguous".
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

DESIGN = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("toy_harness", DESIGN / "toy_harness.py")
TH = importlib.util.module_from_spec(_s)
sys.modules["toy_harness"] = TH
_s.loader.exec_module(TH)

from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"

#: ⛔⛔ **THE ORACLE THIS INSTRUMENT NEEDS CANNOT BE BUILT ON THE SHIPPED THREAD DEFAULT.**
#: `BamScanConfig.total_threads = 0` means all cores, and float addition is not associative across
#: worker threads — the same lesson `rescan_panels.py` records. `OracleTruth._validate`'s sum-to-full
#: budget is `max(|full|, 1) · n_cells · eps`, DERIVED and per cell, so a boundary whose own mass is
#: below 1.0 gets a floor budget of `n_cells · eps` (~1.6e-15 here) and cannot absorb thread dust.
#: ⭐ Measured on this rung, TA = 0 / TB = 30, four runs: at `total_threads=0` the residual is
#: NON-DETERMINISTIC — `boundary_spliced_mass` 2.7285e-12 then 8.1855e-12 on identical input — and at
#: `total_threads=1` it is deterministic AND both spliced banks go to EXACTLY ZERO
#: (`boundary_spliced_mass` and `sj_mass` drop out of the residual report entirely).
#: ⚠ So this is not a tolerance being widened: it is the instrument stopping the re-association that
#: manufactured the dust. `rename_identity.py` pins the same two knobs for the same reason.
DETERMINISM = dict(total_threads=1, bgzf_threads=1)


def _spec(base, ta: float, tb: float):
    """The rung with ONE thing changed per cell: the two abundances."""
    genes = [
        {
            "gene_id": g["gene_id"],
            "strand": g["strand"],
            "transcripts": [
                dict(t, abundance=(ta if t["t_id"] == "TA" else tb)) for t in g["transcripts"]
            ],
        }
        for g in base.genes
    ]
    return dataclasses.replace(base, genes=genes)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="tes_readthrough")
    ap.add_argument("--donor", default="gdna_g50_ss_0.50_nrna_none_capture_off")
    ap.add_argument("--ta", nargs="*", type=float, default=[0.0, 30.0, 300.0, 3000.0])
    ap.add_argument("--tb", nargs="*", type=float, default=[0.0, 30.0, 300.0, 3000.0])
    ap.add_argument("--n-rna", type=int, default=200_000)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_cert_audit"))
    # ⭐ DEFAULTS TO `on`, AND THAT IS THE HONEST DEFAULT HERE: check (c) — "is a precision emitted" —
    # reads `cm_p`/`cm_n`, which ARE the relay's own RNA measurement precisions. Under the shipped mute
    # they do not exist, so the instrument's headline question has no content and (c) can only be
    # reported UNANSWERED. The stamp says on every run that this is not the shipped configuration.
    TH.add_messages_flag(ap, default=True)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(INDEX))
    messages = TH.messages_on(args)
    cfg = TH.with_messages(
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=0), messages
    )
    # ⛔ ONE config drives the whole toy AND the oracle's four partition scans — `run_toy` hands this
    #   same object to `OracleTruth.from_bam`, so pinning here pins every scan sum-to-full compares.
    pipe = PipelineConfig()
    pipe = dataclasses.replace(pipe, scan=dataclasses.replace(pipe.scan, **DETERMINISM))
    donor = TH.harvest(SUITE / args.donor, index, config=cfg)
    base = dataclasses.replace(TH.SPECS[args.spec], n_rna_fragments=int(args.n_rna))

    print("=" * 148)
    print(f"⭐⭐⭐ CERTIFIED-RNA CHANNEL AUDIT — {args.spec} on a TA x TB abundance grid")
    print("=" * 148)
    print(f"   donor {args.donor}   gDNA {donor.gdna_rate_per_base:.6g}/bp   prior-free pass-0")
    print(f"   TA in {args.ta}   TB in {args.tb}   RNA budget {args.n_rna:,}")
    print(TH.messages_stamp(messages))
    print("   ⛔ (a) bank populated?  (b) divisor > 0?  (c) RNA precision emitted?"
          "  — any one failing gives the same symptom")

    fails: list[str] = []
    #: ⭐ banks where (a) and (b) pass and (c) was NEVER ASKED, because the relay is muted. Kept apart
    #: from `fails` so a config flag can never be reported as a defect in the certified-RNA channel.
    unasked: list[str] = []
    for ta in args.ta:
        for tb in args.tb:
            if ta == 0.0 and tb == 0.0:
                continue  # a silent gene has no certified channel to audit
            try:
                r = TH.run_toy(
                    _spec(base, ta, tb), donor, args.work_dir,
                    config=cfg, pipeline_config=pipe,
                )
            except AssertionError as exc:
                # ⛔ The oracle refused its own partition. With the scan pinned single-threaded that
                #   is NOT float dust any more, so it must not be swallowed — but a bare traceback
                #   out of `_oracle` reads as "this instrument is broken" when the diagnosis is one
                #   layer down. Name the layer and stop; do not continue with a truth-free grid.
                # ⚠ The verdict depends on whether the scan WAS pinned, so it is read off the config
                #   rather than asserted — a run with `DETERMINISM` edited must not claim the pin held.
                pinned = int(pipe.scan.total_threads) == 1
                print(f"\n⛔⛔ ORACLE REFUSED at TA = {ta:g} / TB = {tb:g}. Scan config: {DETERMINISM}")
                print(
                    "   ⭐ single-threaded, so this is NOT thread re-association."
                    if pinned
                    else "   ⚠ NOT pinned single-threaded — expect exactly this, non-deterministically."
                )
                print(f"   {exc}")
                print("   ⭐ The sum-to-full budget lives in `tests/calibration/_oracle.py::_validate`")
                print("     and is out of reach from `scripts/`. Re-derive it there, or report the")
                print("     grid cell above: the partition genuinely differs from the production scan.")
                return 1
            cap = r.capture
            st = TH.relay_static(cap)  # E_r survives the mute — the BACKBONE publishes it
            # ⛔ `_uni` exists only under `RelayPolicy`, and CHECK (c) IS a relay quantity: `cm_p`/`cm_n`
            # are the relay's own RNA measurement precisions. Muted, (c) is UNANSWERABLE, never FAILED —
            # see the `_UNASKED` note below. This file used to die here with `KeyError: '_uni'`.
            uni = TH.relay_channels(cap)
            spl = np.asarray(cap["spliced"], float)
            jun = np.asarray(cap["mature"], float)
            E_r = np.asarray(st["E_r"], float)
            nan_slots = np.full(int(r.chain.n_slots), float("nan"))
            cm_p = np.asarray(uni["cm_p"], float) if uni is not None else nan_slots
            cm_n = np.asarray(uni["cm_n"], float) if uni is not None else nan_slots
            fg = np.asarray(cap["f_g"], float)
            loc = np.asarray(cap["fg_loc"], float)
            cnt = np.asarray(cap["count"], float).sum(axis=1)
            rows = TH.object_rows(r)
            # ⚠ `ScanQC` is a dataclass, not a mapping — read the counter by attribute and
            #   fall back to 0 so a schema change degrades to "unknown" rather than crashing.
            # ⛔⛔ AND `deferred_undetermined_gap` ALONE IS THE WRONG COLUMN NOW. `run_toy` runs the
            #   SECOND PASS on the whole before it returns, so the payload it hands back is DRAINED and
            #   that counter reads 0 on every cell — which a reader takes as "no fragment was ever
            #   ambiguous" when what it means is "every ambiguous one has since been placed". The drain's
            #   own `offered`/`deposited` are what carry the population, so all three are printed.
            qc = getattr(r.payload, "qc", None)
            deferred = int(getattr(qc, "deferred_undetermined_gap", 0) or 0)
            drain = getattr(r.payload, "drain", None)
            held = int(getattr(drain, "offered", 0) or 0)
            placed = int(getattr(drain, "deposited", 0) or 0)

            print(f"\n{'─' * 148}")
            print(f"⭐ TA = {ta:<8g} TB = {tb:<8g}"
                  f"   ambiguous mate gap (CORRECT): {held:,} held → {placed:,} placed by the second"
                  f" pass, {deferred:,} still deferred")
            # ⭐⭐ `Δgdna` IS THE COLUMN THIS FILE'S OWN RULE DEMANDS AND THE TABLE HAD STOPPED PRINTING.
            #    "SCORE THE MASS, NOT `f_g`" was unfollowable off this table: every printed column was a
            #    fraction. `object_rows` publishes `err` = |PRED gDNA − TRUE gDNA| in FRAGMENTS and `mass`
            #    = the object's own true total, so both are printed and the rule is executable again.
            print(f"   {'object':<28}{'n':>8}{'spliced':>9}{'sj':>10}{'E_r':>8}"
                  f"{'cm_p':>10}{'cm_n':>8}{'true_fg':>9}{'fg_loc':>8}{'pred_fg':>9}{'Δ':>9}"
                  f"{'Δgdna':>10}{'mass':>10}   audit")
            print("   " + "-" * 146)
            for row in rows:
                s = row["slot"]
                if cnt[s] <= 0 and spl[s] <= 0:
                    continue
                S, J, er = float(spl[s]), float(jun[s]), float(E_r[s])
                pp, pn = float(cm_p[s]), float(cm_n[s])
                tf = row["true_fg"]
                d = fg[s] - tf if tf == tf else float("nan")
                # ⭐ the audit: only meaningful where the bank IS populated.
                if S <= 0:
                    note = ""
                elif er <= 0:
                    note = "⛔⛔ (b) certified RNA with NO DIVISOR — cannot become a density"
                    fails.append(f"TA={ta} TB={tb} {row['type']} {row['where']}: spliced>0, E_r=0")
                elif uni is None:
                    # ⛔⛔ NOT A FAILURE. `cm_p`/`cm_n` are the RELAY's precisions and the relay is muted,
                    # so no precision was ASKED FOR. Calling this "(c) FAILS" would manufacture a verdict
                    # against the certified-RNA channel out of a message-layer config flag — the exact
                    # false finding `TRAPS: an-ablation-that-never-ran` warns about.
                    # ⛔⛔ THE TEST IS `uni is None` — THE RELAY'S OWN ARTIFACT — AND NOT `cm_p != cm_p`.
                    #   Keying on NaN conflates two different verdicts: muted (`_uni` absent, nothing was
                    #   asked) and INSTALLED-BUT-PUBLISHED-A-NON-NUMBER, which is a (c) FAILURE. A NaN
                    #   under `RelayPolicy` would have printed "the relay is MUTED", which is false, and
                    #   filed a real wiring defect under `unasked` where nothing reads it. ⭐ Same
                    #   discipline as `toy_harness.relay_live`: trust the artifact, never the flag —
                    #   and here, never a sentinel value the artifact might legitimately carry.
                    note = "⚠ (c) UNANSWERABLE — the relay is MUTED, no precision was asked for"
                    unasked.append(f"TA={ta} TB={tb} {row['type']} {row['where']}")
                elif not (pp > 0.0 or pn > 0.0):
                    # ⚠ Written as `not (… > 0)` rather than `<= 0` so a NaN published by an INSTALLED
                    #   relay lands here — no usable precision is no usable precision.
                    note = "⛔⛔ (c) certified RNA, divisor OK, but NO PRECISION EMITTED — inert"
                    fails.append(f"TA={ta} TB={tb} {row['type']} {row['where']}: spliced>0, prec=0")
                else:
                    note = f"✅ certified RNA live (rho_R = {S / er:.4g})"
                f = lambda v, w=9, p=4: (f"{v:>{w}.{p}f}" if v == v else f"{'—':>{w}}")  # noqa: E731
                # ⛔ `cm_p`/`cm_n` MUST render as `—` and never as `nan` when the relay is muted. The
                #   run stamp promises "reads '—' because the relay SENT NOTHING — not because it sent
                #   a zero", and a printed `nan` says the relay produced one, which is a different and
                #   wrong story. This column went through a raw `.4g` until 2026-08-17.
                g = lambda v, w: (f"{v:>{w}.4g}" if v == v else f"{'—':>{w}}")  # noqa: E731
                print(f"   {row['type'] + ' ' + row['where']:<28}{cnt[s]:>8,.0f}{S:>9,.0f}"
                      f"{J:>10,.0f}{er:>8.1f}{g(pp, 10)}{g(pn, 8)}{f(tf)}{f(loc[s], 8)}"
                      f"{f(fg[s])}{f(d, 9, 4)}{row['err']:>10,.0f}{row['mass']:>10,.0f}   {note}")

    print(f"\n{'=' * 148}")
    print(TH.messages_stamp(messages))
    if fails:
        print(f"⛔⛔ {len(fails)} CERTIFIED-RNA WIRING FAILURE(S) — a bank populated but not usable:")
        for f_ in fails[:30]:
            print(f"     - {f_}")
    elif unasked:
        # ⛔ The exit code must NOT be 0 here: (a) and (b) passed but (c) — this instrument's whole
        # verdict — was never asked. A green exit would let a muted run stand in for an audit.
        print(f"⚠⚠ (a) and (b) PASS, and (c) IS UNANSWERED on {len(unasked)} populated bank(s): the")
        print("   relay is MUTED, so `cm_p`/`cm_n` do not exist. ⛔ This is NOT a pass and NOT a")
        print("   failure of the certified-RNA channel. Re-run with `--messages on` to audit (c).")
        for f_ in unasked[:30]:
            print(f"     - {f_}")
    else:
        print("✅ every populated `boundary_spliced` bank has a divisor AND emits a precision.")
        print("   ⚠ That is (a)+(b)+(c) only. It does NOT say the precision is the RIGHT SIZE, nor that")
        print("     ψ uses it for the object's OWN belief — ψ has no spliced term at all")
        print("     (`tests/calibration/test_vertex_reference.py`'s certified-RNA-blindness gate pins that).")
    print("=" * 148)
    return 1 if (fails or unasked) else 0


if __name__ == "__main__":
    raise SystemExit(main())
