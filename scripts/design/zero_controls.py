#!/usr/bin/env python
"""Does the tool hold at zero RNA and at zero gDNA? The two zero controls, which belong on every
experiment. Each arm runs a toy spec through the toy harness on a donor condition and prints every
object's answer beside a truth that is a constant, so every deviation is a false positive with
nothing to cancel against it: the zero-RNA arm silences every transcript (``abundance = 0``, a
silent gene) on a gDNA-rich donor and the truth is ``f_g = 1`` at every object; the zero-gDNA arm
runs the spec on the ``g00`` donor, whose measured gDNA rate is 0/bp, and the truth is ``f_g = 0``
at every object that carries RNA. The zero-RNA arm is the biologically dominant case, not a corner:
most annotated transcripts are off in any one sample, so their pure-gDNA objects are the modal
case. Per object it prints the three rungs — ``fg_strand`` (the strand likelihood alone),
``fg_loc`` (the message-free self-solve) and ``f_g`` (the final answer after the message passes) —
because a wrong ``fg_loc`` is an initialisation defect no message caused, while a right ``fg_loc``
and a wrong ``f_g`` is the messages; the localisation table reports the fraction of the gap the
messages closed beside whether the object had own evidence, since an evidence-free object is
supposed to return the reference. It runs the shipped message setting by default; muted, it prints
the measured rung-3 vs rung-2 separation rather than implying a message was sent. An object with
zero counts is flagged EMPTY and an arm whose every object is empty is reported as degenerate
(`TRAPS: could-the-arm-have-fired`): judge an arm by how many objects carried mass, never by the
total. Prior-free pass-0 (``calib_refit_iters = 0``), capture-OFF, unstranded throughout.

Usage::

    python scripts/design/zero_controls.py                                       # both arms, the default specs
    python scripts/design/zero_controls.py --specs silent spliced_exons TA_single_exon nested_exons
    python scripts/design/zero_controls.py --specs spliced_exons --arms rna       # zero-RNA only
    python scripts/design/zero_controls.py --arms gdna --n-rna 50000 --work-dir /path/to/work
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from _shared import sibling  # noqa: E402

TH = sibling("toy_harness.py")

from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"

#: capture-OFF x unstranded throughout: the simplest regime, no enrichment landscape and exactly zero
#: strand information, so nothing can mask a defect in the length or count channels. `g98` is the
#: rung with the most gDNA; the zero-RNA arm silences RNA by construction, so its thin RNA side costs
#: this control nothing.
DONOR_GDNA = "gdna_g98_ss_0.50_nrna_mid_capture_off"  # plenty of gDNA — the ZERO-RNA arm's substrate
DONOR_NONE = "gdna_g00_ss_0.50_nrna_mid_capture_off"  # zero gDNA — the ZERO-gDNA arm's substrate

FAIL: list[str] = []


def silence(spec):
    """Every transcript's abundance to 0 — a gene that is annotated and NOT expressed.

    ``n_rna_fragments`` is kept at 1, not 0: the simulator needs a nonzero RNA budget to run and a
    single fragment is the smallest thing that is still a library. `ToySpec.genes` is a list of dicts
    shared with the module-level `SPECS`, so it is deep-copied — mutating it in place would change
    every later run in the same process."""
    genes = copy.deepcopy(spec.genes)
    for g in genes:
        for t in g["transcripts"]:
            t["abundance"] = 0.0
            t["nrna_abundance"] = 0.0
    return dataclasses.replace(spec, genes=genes, n_rna_fragments=1, nrna_abundance=0.0)


def report(spec_name, arm, r, expect, messages):
    """Per object: the counts, the three rungs, and the deviation from a CONSTANT truth."""
    cap = r.capture
    fg = np.asarray(cap["f_g"], float)
    loc = np.asarray(cap["fg_loc"], float)
    strand = np.asarray(cap["fg_strand"], float)
    tau = np.asarray(cap["_tau0_lam"], float)
    cnt = np.asarray(cap["count"], float).sum(axis=1)
    rows = TH.object_rows(r)
    print(f"\n── {spec_name} · {arm} · truth f_g = {expect:.3f} at EVERY object ──────────────────")
    print(f"   {'slot':>4} {'type':<16} {'where':<16} {'n':>8} {'fg_strand':>10} {'fg_loc':>8} "
          f"{'f_g':>8} {'Δ':>8} {'tau':>9} {'err frags':>10}")
    print("   " + "-" * 112)
    live = tot_err = 0
    worst = (0.0, None, None)
    for row in rows:
        s = row["slot"]
        n = float(cnt[s])
        if n <= 0:
            print(f"   {s:>4} {row['type']:<16} {row['where']:<16} {0:>8} "
                  f"{'—':>10} {'—':>8} {'—':>8} {'—':>8} {'—':>9}  ⚠ EMPTY — not a control")
            continue
        live += 1
        d = float(fg[s]) - expect
        err = abs(d) * n
        tot_err += err
        if abs(d) > abs(worst[0]):
            worst = (d, row["type"], row["where"])
        flag = "" if abs(d) < 1e-6 else ("  ⛔" if abs(d) > 0.01 else "  ⚠")
        print(f"   {s:>4} {row['type']:<16} {row['where']:<16} {n:>8,.0f} {strand[s]:>10.4f} "
              f"{loc[s]:>8.4f} {fg[s]:>8.4f} {d:>+8.4f} {tau[s]:>9.3g} {err:>10.1f}{flag}")
    mass = float(cnt.sum())
    print(f"   {'':4} {'TOTAL':<16} {f'{live} live objects':<16} {mass:>8,.0f} {'':>10} {'':>8} "
          f"{'':>8} {'':>8} {'':>9} {tot_err:>10.1f}")
    if live == 0:
        print("   ⛔⛔ EVERY OBJECT IS EMPTY — this arm is DEGENERATE and tests nothing (TRAPS: could-the-arm-have-fired)")
        return
    print(f"   ⭐ worst object: {worst[1]} {worst[2]}  Δ = {worst[0]:+.4f}   "
          f"·   error share of mass = {tot_err / max(mass, 1):.4%}")
    # rung 3 vs rung 2, measured. Under `SilentPolicy` nothing is sent, so the third rung can only
    # repeat the second — but that is a claim about the code, and this instrument prints the NUMBER
    # instead: `max|f_g − fg_loc|` over the live objects, which is also the proof that the muted arm
    # is not silently doing something (`TRAPS: an-ablation-that-never-ran`).
    _live = np.asarray([float(cnt[row["slot"]]) > 0 for row in rows], bool)
    _sl = np.asarray([row["slot"] for row in rows], np.int64)[_live]
    _sep = float(np.max(np.abs(fg[_sl] - loc[_sl]))) if _sl.size else float("nan")
    print(f"   ⭐ RUNG 3 − RUNG 2, over the {live} live objects: max|f_g − fg_loc| = {_sep:.3g}   "
          f"— {'what the messages MOVED' if messages else 'muted: what no third rung MEASURES'}")
    # the localisation, which must NOT be a raw comparison of |Δ fg_loc| against |Δ f_g|: an object
    # with no own composition evidence is SUPPOSED to return psi's uninformative reference at zero
    # precision, so a rule that reads a large |Δ fg_loc| as "the self-solve is broken" mislabels every
    # evidence-free object. The honest quantity is what fraction of the gap the messages CLOSED,
    # reported per object beside whether it has own evidence, so the two kinds are read apart.
    print(f"\n   {'slot':>4} {'own evidence?':<14} {'|Δ| local':>10} {'|Δ| final':>10} "
          f"{'gap closed':>11}   reading")
    fin_bad = 0.0
    for row in rows:
        s = row["slot"]
        if cnt[s] <= 0:
            continue
        dl, df = abs(float(loc[s]) - expect), abs(float(fg[s]) - expect)
        fin_bad = max(fin_bad, df)
        has_own = float(tau[s]) > 1e-4
        closed = (1.0 - df / dl) if dl > 1e-12 else float("nan")
        if df < 1e-6:
            reading = "✅ exact"
        elif not has_own:
            # the wording is the measurement: muted, "the messages did not carry it" is trivially true
            # of a policy that was never asked and would read as a message-layer verdict, so say which.
            if not messages:
                reading = "⛔ NO own evidence, and the messages are MUTED — nothing COULD carry it"
            else:
                reading = ("⚠ NO own evidence — carried by messages, which stopped short"
                           if closed == closed and closed > 0.5
                           else "⛔ NO own evidence and the messages did not carry it")
        else:
            reading = "⛔ has own evidence and is still wrong — an INITIALISATION defect"
        cl = f"{closed:>10.1%}" if closed == closed else f"{'—':>10}"
        print(f"   {s:>4} {('yes' if has_own else 'no'):<14} {dl:>10.4f} {df:>10.4f} {cl:>11}"
              f"   {reading}")
    label = f"{spec_name} · {arm}"
    if fin_bad > 0.01:
        FAIL.append(f"{label}: worst |Δ| = {fin_bad:.4f} on a CONSTANT truth")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--specs", nargs="*", default=["silent", "TA_single_exon", "spliced_exons"])
    ap.add_argument("--arms", nargs="*", default=["rna", "gdna"], choices=["rna", "gdna"])
    ap.add_argument("--n-rna", type=int, default=200_000, help="the ZERO-gDNA arm's RNA depth")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_zero_controls"))
    # defaults to the SHIPPED setting: this is the admissibility control on every experiment, so it
    # must run the configuration the tool ships.
    TH.add_messages_flag(ap, default=TH.MESSAGES_SHIPPED)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(INDEX))
    messages = TH.messages_on(args)
    config = TH.with_messages(
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=0), messages
    )
    print("=" * 118)
    print("⭐⭐⭐ THE TWO ZERO CONTROLS — a CONSTANT truth, so every deviation is a false positive")
    print("=" * 118)
    print("   capture OFF · unstranded · PRIOR-FREE pass-0 (calib_refit_iters=0)")
    print(TH.messages_stamp(messages))

    if "rna" in args.arms:
        donor = TH.harvest(SUITE / DONOR_GDNA, index, config=config)
        print(f"\n{'=' * 118}\n⭐⭐ ARM 1 — ZERO RNA (silent genes).  donor {DONOR_GDNA}"
              f"\n   gDNA rate {donor.gdna_rate_per_base:.6g}/bp.  ⛔ Truth is f_g = 1.000 EVERYWHERE."
              f"\n   ⭐ This is the biologically dominant case: most annotated transcripts are OFF."
              f"\n{'=' * 118}")
        for name in args.specs:
            spec = silence(TH.SPECS[name])
            report(name, "ZERO RNA", TH.run_toy(spec, donor, args.work_dir / "rna", config=config), 1.0, messages)

    if "gdna" in args.arms:
        donor = TH.harvest(SUITE / DONOR_NONE, index, config=config)
        print(f"\n{'=' * 118}\n⭐⭐ ARM 2 — ZERO gDNA.  donor {DONOR_NONE}"
              f"\n   gDNA rate {donor.gdna_rate_per_base:.6g}/bp.  ⛔ Truth is f_g = 0.000 at every"
              f" object that carries RNA.\n{'=' * 118}")
        for name in args.specs:
            spec = dataclasses.replace(TH.SPECS[name], n_rna_fragments=int(args.n_rna))
            report(name, "ZERO gDNA", TH.run_toy(spec, donor, args.work_dir / "gdna", config=config), 0.0, messages)

    print("\n" + "=" * 118)
    if FAIL:
        print(f"⛔ {len(FAIL)} ZERO CONTROL(S) OFF BY MORE THAN 0.01 ON A CONSTANT TRUTH:")
        for f in FAIL:
            print(f"     - {f}")
    else:
        print("✅ every zero control is exact to 0.01")
    print("=" * 118)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
