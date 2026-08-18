#!/usr/bin/env python
"""⭐⭐⭐ THE TWO ZERO CONTROLS, ON EVERY RUNG — and they belong on every experiment.

⛔ **Owner's standing requirement, 2026-08-05:** *"We generally need zero gDNA and zero RNA controls
everywhere, generally, for almost any experiment and investigation that we run."*

Two arms, and the point of both is that **the truth is a CONSTANT**, so every deviation is a false
positive with nothing to cancel against it:

| arm | how | truth at every object | what a deviation means |
|---|---|---|---|
| ⭐⭐ **ZERO RNA** | every transcript's ``abundance = 0`` — a SILENT gene | ``f_g = 1.000`` exactly | invented RNA in a pure-gDNA library |
| ⭐⭐ **ZERO gDNA** | the ``g00`` donor, whose measured gDNA rate is 0/bp | ``f_g = 0.000`` exactly | invented gDNA — the phantom floor |

⭐⭐ **THE ZERO-RNA ARM IS THE BIOLOGICALLY DOMINANT CASE AND IT IS NOT A CORNER.** The annotation has
>50,000 genes and perhaps ~10,000 are expressed in any one sample, so **most annotated transcripts are
simply OFF**. Their objects are pure gDNA, they are the majority of the genome's objects, and they carry
real mass — so an error there is not an boundary case, it is the modal case. It should also be the easiest
thing the solver ever does: there is nothing to deconvolve.

**The three-rung ladder is printed per object, because it localises the defect with no guessing:**

    fg_strand   the strand likelihood ALONE  — on an unstranded library this MUST be psi's reference
    fg_loc      the message-free SELF-SOLVE  — strand + the intron factory + psi's own reference
    f_g         the FINAL answer, after the forward-backward relay

⭐ If ``fg_loc`` is already wrong the fault is in the per-object initialisation and no message caused it.
If ``fg_loc`` is right and ``f_g`` is not, it is the messages. ⛔ Reading only the final number cannot tell
those apart, and they have completely different fixes.

⛔⛔ **AND UNDER THE SHIPPED CONFIG THERE IS NO THIRD RUNG.** ``CalibrationConfig.message_propagation``
is ``False``, which installs ``SilentPolicy`` — no relay runs, so rung 3 can only repeat rung 2. This
instrument **prints the two rungs' measured maximum separation** rather than implying a relay ran, and
its relay-derived column reads ``—`` rather than a zero. ⭐ ``--messages on`` restores the third rung and
the stamp then says the configuration is not the shipped one. ⚠ It defaults to the SHIPPED setting on
purpose: this is the admissibility control the owner requires on **every experiment**, and an experiment
runs the configuration the tool ships.

⚠ **A zero arm can be DEGENERATE, and then it is not a control at all** (TRAPS: could-the-arm-have-fired). An object with
zero counts has no density, so the reframe is skipped and every message into it is inert — the arm then
reports "no error" while testing nothing. So this prints the COUNTS beside every answer and flags any
object that is empty. ⛔ Judge the arm by how many objects actually carried mass, not by the total.

Usage:

    python scripts/design/zero_controls.py --specs silent spliced_exons TA_single_exon nested_exons
    python scripts/design/zero_controls.py --specs spliced_exons --arms rna     # zero-RNA only
"""

from __future__ import annotations

import argparse
import copy
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

#: ⭐ capture-OFF × unstranded throughout: the simplest regime, no enrichment landscape and exactly zero
#: strand information, so nothing can mask a defect in the length or count channels.
# ⚠ `g75` until 2026-08-13, retired when the ladder was rebuilt to four rungs. `g98` is the surviving
# rung with the MOST gDNA (9.8 M fragments against g75's 7.5 M), which is what "plenty" asks for; the
# arm silences RNA by construction, so g98's thin RNA side costs this control nothing.
DONOR_GDNA = "gdna_g98_ss_0.50_nrna_none_capture_off"  # plenty of gDNA — the ZERO-RNA arm's substrate
DONOR_NONE = "gdna_g00_ss_0.50_nrna_none_capture_off"  # zero gDNA — the ZERO-gDNA arm's substrate

FAIL: list[str] = []


def silence(spec):
    """Every transcript's abundance to 0 — a gene that is annotated and NOT expressed.

    ⚠ ``n_rna_fragments`` is kept at 1, not 0: the simulator needs a nonzero RNA budget to run and a
    single fragment is the smallest thing that is still a library. ⛔ `ToySpec.genes` is a list of dicts
    shared with the module-level `SPECS`, so it is DEEP-copied — mutating it in place silently changed
    every later run in the same process, which is how this was first mis-measured."""
    genes = copy.deepcopy(spec.genes)
    for g in genes:
        for t in g["transcripts"]:
            t["abundance"] = 0.0
            t["nrna_abundance"] = 0.0
    return dataclasses.replace(spec, genes=genes, n_rna_fragments=1, nrna_abundance=0.0)


def report(spec_name, arm, r, expect):
    """Per object: the counts, the three rungs, and the deviation from a CONSTANT truth."""
    cap = r.capture
    # ⭐ read off the ARTIFACT, not off the flag: `head.py` is `_uni`'s only writer, so its presence is
    # the relay's own signature and a config that failed to thread through cannot lie about it.
    relay = TH.relay_channels(cap)
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
    # ⭐⭐ RUNG 3 vs RUNG 2, MEASURED. Under `SilentPolicy` the relay sends nothing, so the third rung can
    # only repeat the second — but that is a claim about the code, and this instrument's job is to print
    # the NUMBER instead. `max|f_g − fg_loc|` over the live objects is that number, and it is also the
    # honest proof that the muted arm is not silently doing something (`TRAPS: an-ablation-that-never-ran`).
    _live = np.asarray([float(cnt[row["slot"]]) > 0 for row in rows], bool)
    _sl = np.asarray([row["slot"] for row in rows], np.int64)[_live]
    _sep = float(np.max(np.abs(fg[_sl] - loc[_sl]))) if _sl.size else float("nan")
    if relay is None:
        print(f"   ⭐ RUNG 3 − RUNG 2, over the {live} live objects: max|f_g − fg_loc| = {_sep:.3g}   "
              f"— the relay is MUTED, so this is what 'no third rung' MEASURES")
    else:
        print(f"   ⭐ RUNG 3 − RUNG 2, over the {live} live objects: max|f_g − fg_loc| = {_sep:.3g}   "
              f"— what the relay MOVED")
    # ⛔⛔ THE LOCALISATION, and it must NOT be a raw comparison of |Δ fg_loc| against |Δ f_g|.
    # An object with no own composition evidence is SUPPOSED to return psi's uninformative reference
    # (~0.49) at zero precision — that is correct behaviour, not an error, so a rule that reads a large
    # |Δ fg_loc| as "the self-solve is broken" mislabels every evidence-free object. (It did, on the
    # first run of this file: it called an exon whose messages had moved it 0.49 → 0.93 a self-solve
    # defect.) ⭐ The honest quantity is what fraction of the gap the messages CLOSED, reported per
    # object beside its own precision, so an evidence-free object and an evidence-bearing one are read
    # differently rather than pooled.
    # ⭐⭐⭐ AND THE DECISIVE COLUMN: what the MESSAGES actually delivered, as an implied f_g. The relay
    # hands each slot a fused gDNA DENSITY ``cg``; the composition that density implies at this slot is
    # ``cg . E_g / M``. If that already implies f_g >= 1 and psi still returns less, the shortfall is
    # psi's SOLVE and not the level — a completely different fix from anything in the message layer, and
    # no amount of work on the reframe can reach it.
    # ⛔ `_uni` exists only under `HeadPolicy`; muted, there is NO delivered density to convert, and a 0
    # here would read as "the relay delivered nothing useful" rather than "the relay was not asked".
    cg = np.asarray(relay["cg"], float) if relay is not None else None
    st = TH.relay_static(r.capture)  # M / E_g survive the mute — the BACKBONE publishes them
    M = np.asarray(st["M"], float)
    E_g = np.asarray(st["E_g"], float)
    if relay is None:
        print(TH.relay_silent_note("THE 'msg implies f_g' COLUMN"))
    print(f"\n   {'slot':>4} {'own evidence?':<14} {'|Δ| local':>10} {'|Δ| final':>10} "
          f"{'gap closed':>11} {'msg implies f_g':>16}   reading")
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
            # ⛔ THE WORDING IS THE MEASUREMENT. Muted, "the messages did not carry it" is trivially true
            # of a relay that was never asked, and reads as a message-layer verdict. Say which it is.
            if relay is None:
                reading = "⛔ NO own evidence, and the relay is MUTED — nothing COULD carry it"
            else:
                reading = ("⚠ NO own evidence — carried by messages, which stopped short"
                           if closed == closed and closed > 0.5
                           else "⛔ NO own evidence and the messages did not carry it")
        else:
            reading = "⛔ has own evidence and is still wrong — an INITIALISATION defect"
        cl = f"{closed:>10.1%}" if closed == closed else f"{'—':>10}"
        imp = (cg[s] * E_g[s] / M[s] if M[s] > 0 else float("nan")) if relay is not None \
            else float("nan")
        ims = f"{imp:>16.4f}" if imp == imp else f"{'—':>16}"
        if expect == 1.0 and imp == imp and imp >= 0.999 and df > 1e-6:
            reading = "⛔⛔ the MESSAGE already implies f_g >= 1 — the shortfall is psi's SOLVE"
        print(f"   {s:>4} {('yes' if has_own else 'no'):<14} {dl:>10.4f} {df:>10.4f} {cl:>11} "
              f"{ims}   {reading}")
    label = f"{spec_name} · {arm}"
    if fin_bad > 0.01:
        FAIL.append(f"{label}: worst |Δ| = {fin_bad:.4f} on a CONSTANT truth")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--specs", nargs="*", default=["silent", "TA_single_exon", "spliced_exons"])
    ap.add_argument("--arms", nargs="*", default=["rna", "gdna"], choices=["rna", "gdna"])
    ap.add_argument("--n-rna", type=int, default=200_000, help="the ZERO-gDNA arm's RNA depth")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_zero_controls"))
    # ⭐ defaults to the SHIPPED setting: this is the admissibility control the owner requires on every
    # experiment, so it must run the configuration the tool ships.
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
            report(name, "ZERO RNA", TH.run_toy(spec, donor, args.work_dir / "rna", config=config), 1.0)

    if "gdna" in args.arms:
        donor = TH.harvest(SUITE / DONOR_NONE, index, config=config)
        print(f"\n{'=' * 118}\n⭐⭐ ARM 2 — ZERO gDNA.  donor {DONOR_NONE}"
              f"\n   gDNA rate {donor.gdna_rate_per_base:.6g}/bp.  ⛔ Truth is f_g = 0.000 at every"
              f" object that carries RNA.\n{'=' * 118}")
        for name in args.specs:
            spec = dataclasses.replace(TH.SPECS[name], n_rna_fragments=int(args.n_rna))
            report(name, "ZERO gDNA", TH.run_toy(spec, donor, args.work_dir / "gdna", config=config), 0.0)

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
