#!/usr/bin/env python
"""⭐⭐⭐ IS THE CERTIFIED-RNA CHANNEL WIRED AT A TERMINUS EDGE? — the audit, on a TA x TB abundance grid.

A **spliced** fragment cannot be gDNA, so ``edge_spliced`` is the one observation in the tool that needs no
deconvolution: it is certified RNA. The solver's job is to deconvolve the **unspliced** population, and any
certified fragment that leaks into it is RNA handed to the gDNA solver.

⛔ **THE CASE THIS EXISTS FOR** (owner, 2026-08-05). The accumulator now creates an EDGE at every
transcript start and end. **Those EDGEs carry no splice junction** — nothing splices there — **yet spliced
fragments cross them**, because transcription of a longer overlapping transcript continues straight past.
Unless such a fragment is binned as ``edge_spliced`` it lands in the unspliced pool and is deconvolved.

The rung is ``tes_readthrough``::

    TA+ (1,050, 2,000) (9,000,  9,100)      junction 2,000 -> 9,000
    TB+ (1,000, 2,000) (9,050, 11,000)      junction 2,000 -> 9,050

⭐ **EDGE @9,100 is TA's TES and no junction touches it**, but TB's exon 2 runs to 11,000 — so a TB
fragment that used TB's junction and reaches past 9,100 crosses it CONTIGUOUSLY HAVING SPLICED ELSEWHERE.
⭐ **EDGE @9,050 is TB's junction acceptor AND a plain contiguity line for TA**, whose exon 2 spans
9,000-9,100 unbroken. One line, junction flux for one transcript and an unspliced RNA crossing for another.

⛔⛔ **SCORE THE MASS, NOT ``f_g`` — AND THIS INSTRUMENT LEARNED THAT THE HARD WAY.** The solver's
``f_g`` is the gDNA fraction of the **unspliced** population; the oracle's per-object fraction is
spliced-INCLUSIVE. Comparing them directly reads a 1.5 %-of-mass defect as a half-unit error, which is
how the first run of this file mis-reported it. The columns that decide anything are ``TRUE gDNA`` vs
``PRED gDNA`` in fragments.

**What is checked, per object, at every grid cell** — three things, because any one of them failing
silently produces the same symptom:

===  ====================================================================================
(a)  is the BANK populated?          ``spliced_count`` > 0 where the geometry says it must be
(b)  does it have a DIVISOR?         ``eff_rna`` > 0, since there is no ``eff_junction`` at a
                                     terminus EDGE to price it against
(c)  is a PRECISION EMITTED?         the relay's own RNA measurement precision (``cm_p``/``cm_n``)
                                     — a bank with a divisor and no precision is inert
===  ====================================================================================

---

⭐⭐⭐ **WHAT IT MEASURED, 2026-08-05 — 24 of 24 grid cells, and the verdict has two halves.**

✅ **(a) and (b) PASS, and the design intent is implemented.** ``edge_spliced`` is populated at both
terminus/contiguity EDGEs, it has ``eff_rna`` = 202.8 as its divisor, it is held OUT of the deconvolution,
and it is added back as RNA: ``PRED rna`` tracks ``TRUE rna`` to <2 % (11,205 vs 11,380 · 9,001 vs 9,018 ·
2,747 vs 2,784). Certified RNA is not being fed to the gDNA solver.

⛔ **(c) FAILS at every populated bank — ``cm_p = cm_n = 0`` in 24/24 cells.** The channel therefore
informs nothing about the UNSPLICED split at the very same line, and there the gDNA is over-called:

=================  =========  ==========  ==========  ==========  ===========
cell               object     unspl n     spliced     TRUE gDNA   PRED gDNA
=================  =========  ==========  ==========  ==========  ===========
TA 3000 / TB 30    @9,050     365         11,026      **11**      **186**  (17x)
TA 300  / TB 30    @9,050     288          8,739      9           26   (2.9x)
TA 30   / TB 30    @9,050     96           2,698      10          47   (4.7x)
TA 3000 / TB 30    @9,100     116            364      12          12   (exact)
=================  =========  ==========  ==========  ==========  ===========

⭐ **@9,050 is the diagnostic and @9,100 is its control.** Both carry certified RNA; @9,100's flanks
include the 1,900 bp exon [9,100, 11,000) which the relay can speak from, and it is exact. @9,050's flanks
are the two **50 bp** nodes [9,000, 9,050) and [9,050, 9,100) — both below one mean fragment length, so
neither has a resolvable density (TRAPS: density-below-one-fragment-length) and the relay has nothing to offer. ⛔ **So the answer at
@9,050 is set entirely by whether a neighbour happens to be informative, and never by the object's own
11,026 certified-RNA fragments.** ``f_g`` sits at 0.4902-0.5098 — the uninformative reference — across the
whole grid.

⚠ **It is abundance-driven in the direction that confirms the mechanism, not the object's own count.** At
TB = 0 there is no readthrough, @9,050 holds 11,649 certified fragments and reads 0.0725; at TB = 30 it
holds 11,026 and reads 0.5098. More of its own evidence, a worse answer — because what changed was the
neighbour, not the object.

⭐ **The grid is the experiment, not a sweep for its own sake.** At @9,100 the certified channel is TB's
alone while the unspliced crossing there is gDNA + TB, so the TA/TB ratio moves the two independently and a
single abundance pair tests one corner. TA also carries the *other* junction into @9,050, so the two
transcripts stress different lines.

⚠ **A known, CORRECT loss to keep in view.** TA's and TB's junctions share the donor at 2,000 and differ
only in acceptor (9,000 vs 9,050), so a fragment whose UNSEQUENCED mate gap could hold either implies two
different ``L`` values and is DEFERRED to the second pass rather than deposited — the owner's ruling, in
`tests/native/_accumulator_reference.py`. Measured on the pure-RNA arm: 1,298 of 200,000 fragments, which
is **13.5 % of the certified channel at @9,100** (9,627 truth vs 8,329 accumulated). ⛔ So the channel is
systematically UNDER-counted at a terminus EDGE wherever alternative acceptors exist, and that is a
property of the channel rather than a bug. It is printed per cell.
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

from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"


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
    args = ap.parse_args()

    index = TranscriptIndex.load(str(INDEX))
    cfg = dataclasses.replace(CalibrationConfig(), calib_refit_iters=0)
    donor = TH.harvest(SUITE / args.donor, index, config=cfg)
    base = dataclasses.replace(TH.SPECS[args.spec], n_rna_fragments=int(args.n_rna))

    print("=" * 132)
    print(f"⭐⭐⭐ CERTIFIED-RNA CHANNEL AUDIT — {args.spec} on a TA x TB abundance grid")
    print("=" * 132)
    print(f"   donor {args.donor}   gDNA {donor.gdna_rate_per_base:.6g}/bp   prior-free pass-0")
    print(f"   TA in {args.ta}   TB in {args.tb}   RNA budget {args.n_rna:,}")
    print("   ⛔ (a) bank populated?  (b) divisor > 0?  (c) RNA precision emitted?"
          "  — any one failing gives the same symptom")

    fails: list[str] = []
    for ta in args.ta:
        for tb in args.tb:
            if ta == 0.0 and tb == 0.0:
                continue  # a silent gene has no certified channel to audit
            r = TH.run_toy(_spec(base, ta, tb), donor, args.work_dir, config=cfg)
            cap = r.capture
            st = cap["_uni_static"]
            uni = cap["_uni"][-1]
            spl = np.asarray(cap["spliced"], float)
            jun = np.asarray(cap["mature"], float)
            E_r = np.asarray(st["E_r"], float)
            cm_p = np.asarray(uni["cm_p"], float)
            cm_n = np.asarray(uni["cm_n"], float)
            fg = np.asarray(cap["f_g"], float)
            loc = np.asarray(cap["fg_loc"], float)
            cnt = np.asarray(cap["count"], float).sum(axis=1)
            rows = TH.object_rows(r)
            # ⚠ `ScanQC` is a dataclass, not a mapping — read the counter by attribute and
            #   fall back to 0 so a schema change degrades to "unknown" rather than crashing.
            qc = getattr(r.payload, "qc", None)
            deferred = int(getattr(qc, "deferred_undetermined_gap", 0) or 0)

            print(f"\n{'─' * 132}")
            print(f"⭐ TA = {ta:<8g} TB = {tb:<8g}"
                  f"   deferred (ambiguous mate gap, CORRECT): {deferred:,}")
            print(f"   {'object':<28}{'n':>8}{'spliced':>9}{'junction':>10}{'E_r':>8}"
                  f"{'cm_p':>10}{'cm_n':>8}{'true_fg':>9}{'fg_loc':>8}{'pred_fg':>9}{'Δ':>9}   audit")
            print("   " + "-" * 126)
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
                elif pp <= 0 and pn <= 0:
                    note = "⛔⛔ (c) certified RNA, divisor OK, but NO PRECISION EMITTED — inert"
                    fails.append(f"TA={ta} TB={tb} {row['type']} {row['where']}: spliced>0, prec=0")
                else:
                    note = f"✅ certified RNA live (rho_R = {S / er:.4g})"
                f = lambda v, w=9, p=4: (f"{v:>{w}.{p}f}" if v == v else f"{'—':>{w}}")  # noqa: E731
                print(f"   {row['type'] + ' ' + row['where']:<28}{cnt[s]:>8,.0f}{S:>9,.0f}"
                      f"{J:>10,.0f}{er:>8.1f}{pp:>10.4g}{pn:>8.4g}{f(tf)}{f(loc[s], 8)}"
                      f"{f(fg[s])}{f(d, 9, 4)}   {note}")

    print(f"\n{'=' * 132}")
    if fails:
        print(f"⛔⛔ {len(fails)} CERTIFIED-RNA WIRING FAILURE(S) — a bank populated but not usable:")
        for f_ in fails[:30]:
            print(f"     - {f_}")
    else:
        print("✅ every populated `edge_spliced` bank has a divisor AND emits a precision.")
        print("   ⚠ That is (a)+(b)+(c) only. It does NOT say the precision is the RIGHT SIZE, nor that")
        print("     ψ uses it for the object's OWN belief — ψ has no spliced term at all")
        print("     (`tests/calibration/test_vertex_reference.py`'s certified-RNA-blindness gate pins that).")
    print("=" * 132)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
