#!/usr/bin/env python
"""Score two ladder arms against each other, on the DELIVERABLE and on pass-0, per stratum.

⛔ `abs_err_all_final` is the deliverable (TRAPS: the-intermediate-is-not-the-deliverable — a −37.2 % pass-0 win was −3.9 % shipped);
`abs_err_all` is pass-0. BOTH are printed on every row, always.
⛔ `g00` is reported SEPARATELY, never inside the ALL row: its truth is exactly 0, so its relative
change is unbounded and would swamp every other condition in either direction.

⛔⛔ **THE PANEL SIZE IS DERIVED FROM THE ARM FILES, NEVER ASSERTED (repaired 2026-08-17).** The header
printed a hard-coded ``36 conditions x 2 axes`` and this docstring said ``the other 35`` — both were the
36-condition ladder, RETIRED on 2026-08-13 for the 16-condition `gdna_ladder.yaml`. ⚠ The lie was not
merely out of date: a shard of the panel scored *exactly* the same header, so the one line a reader uses
to check coverage confirmed a number the file had never counted. It is now ``len()`` of what was loaded.

⛔⛔ **AND `--messages` IS PART OF THE ARM, so two arms run under DIFFERENT policies are REFUSED BY
DEFAULT.** `ladder_arm_ab.py` stamps `messages` into every row for exactly this reason and
`arm_identity.py` already fails on it (`messages` is not in its `_NOT_A_MEASUREMENT` set); this
AGGREGATES, so pooling an `off` base with an `on` arm would have reported the config flip as the
mechanism's effect and nothing would have said so (TRAPS: an-ablation-that-never-ran).
⚠ An arm file written before 2026-08-11 carries no `messages` field at all; it is reported as
`unstamped` and refused against a stamped one rather than assumed to match.

⭐ **`--across-policies` is the deliberate exception, and it exists because the cross-policy comparison
is a REAL measurement** — `base --messages off` against `base --messages on` is how the mute was priced
at −58.3 % / −43.7 % / −32.1 % on the three in-scope strata and **+154.8 %** on the deferred one. It is
opt-in and it stamps `messages=off→on` across the header, following `ladder_arm_ab.py --allow-inert`:
refuse the accident, keep the experiment, and make the override visible in the output.

Usage::

    RIGEL_ARMS=<dir> python scripts/design/arm_score.py <base-arm> <arm> [--across-policies]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import os

#: where `ladder_arm_ab.py --out` wrote the arms. Override with $RIGEL_ARMS.
D = Path(os.environ.get("RIGEL_ARMS", Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_arms"))


def load(name):
    path = D / f"{name}.jsonl"
    if not path.is_file():
        raise SystemExit(
            f"no arm file at {path}\n"
            f"   `ladder_arm_ab.py --out {D}` writes them; $RIGEL_ARMS points at the directory."
        )
    rows = [json.loads(x) for x in path.read_text().splitlines() if x.strip()]
    if not rows:
        raise SystemExit(f"{path} is EMPTY — an arm with no rows must not score as agreement.")
    return {(r["condition"], r["axis"]): r for r in rows}


def policy_of(rows) -> str:
    """The `messages` setting every row of one arm file agrees on — or what disagrees.

    ⛔ Not defaulted. `ladder_arm_ab.py` stamps the setting on every row precisely so a consumer never
    has to guess it, and a file that predates 2026-08-11 carries no stamp at all, so it reads `unstamped`.

    ⚠ **`unstamped` compares equal to ANOTHER `unstamped`, and that is deliberate — it is a hole, and it
    is named rather than papered over** (measured 2026-08-17: two stripped copies of real arm files score
    `messages=unstamped`, rc=0, no refusal). Two pre-stamp files are the one pair this cannot separate;
    refusing them would break every comparison of two historical arms, which is the case the stamp was
    introduced too late to cover. ⛔ What it DOES catch is the mixed pair — `unstamped` against a stamped
    file is refused, because there the newer file's own stamp says a policy was chosen.
    """
    seen = {r.get("messages", "unstamped") for r in rows.values()}
    return seen.pop() if len(seen) == 1 else "MIXED" + repr(sorted(seen))


def stratum(cond):
    return ("stranded" if "ss_0.99" in cond else "unstranded",
            "capture ON" if "capture_on" in cond else "capture OFF")


#: ⛔ THE 0.8.0 SCOPE, stamped on the stratum rows. `unstranded x capture ON` is DEFERRED (owner,
#: 2026-08-14) and carries 64.5 % of transcript error, so it is the row a ranked read lands on by
#: accident every time. Stamped, not filtered — deferred is not dropped, and it must keep being reported.
DEFERRED_STRATUM = ("unstranded", "capture ON")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("base", help=f"the BASELINE arm's name, i.e. {D}/<base>.jsonl")
    ap.add_argument("arm", help="the arm under test, scored against the baseline")
    ap.add_argument(
        "--across-policies",
        action="store_true",
        help="compare two arms run under DIFFERENT --messages settings. ⛔ Only when the POLICY is the "
             "thing being measured; the difference is then the config flip, not the mechanism.",
    )
    args = ap.parse_args()

    a_name, b_name = args.base, args.arm
    A, B = load(a_name), load(b_name)
    keys = sorted(set(A) & set(B))
    assert len(keys) == len(A) == len(B), f"arms differ in shape: {len(A)} vs {len(B)} vs {len(keys)}"

    # ⛔ The policy is PART OF THE ARM. Refuse before a single number is printed.
    # ⚠ A MIXED file is refused whatever the flag says: `--across-policies` licenses ONE deliberate
    # step between two arms, never two policies smeared inside one of them.
    policy_a, policy_b = policy_of(A), policy_of(B)
    mixed = policy_a.startswith("MIXED") or policy_b.startswith("MIXED")
    if mixed or (policy_a != policy_b and not args.across_policies):
        raise SystemExit(
            f"⛔ REFUSED: these two arms were not run under the same message policy.\n"
            f"   {a_name}: messages={policy_a}\n   {b_name}: messages={policy_b}\n"
            f"   `--messages` is part of the arm (ladder_arm_ab.py stamps it on every row). Aggregating\n"
            f"   across it would report the config flip as the mechanism's effect.\n"
            + ("   ⛔ A file with MIXED stamps is refused outright — no flag licenses that.\n" if mixed
               else "   ⭐ If the POLICY is what you are measuring, pass --across-policies.\n")
        )
    policy_label = policy_a if policy_a == policy_b else f"⛔ {policy_a}→{policy_b} ACROSS POLICIES"

    def agg(field, pred):
        sa = sum(A[k][field] for k in keys if pred(k))
        sb = sum(B[k][field] for k in keys if pred(k))
        return sa, sb

    def line(label, pred, n):
        out = []
        for field in ("abs_err_all_final", "abs_err_all"):
            sa, sb = agg(field, pred)
            pct = (sb - sa) / sa * 100 if sa > 0 else float("nan")
            out.append((sa, sb, pct))
        better = sum(1 for k in keys if pred(k)
                     and B[k]["abs_err_all_final"] < A[k]["abs_err_all_final"])
        (fa, fb, fp), (pa, pb, pp) = out
        print(f"   {label:<26} {fa:>13,.0f} {fb:>13,.0f} {fp:>+8.1f}%   "
              f"{pa:>13,.0f} {pb:>13,.0f} {pp:>+8.1f}%   {better:>3}/{n}")

    def is_g00(k):
        return "_g00_" in k[0]

    # ⛔ COUNTED, never asserted — see the docstring. `len(keys)` is what the tables below actually sum.
    n_cond = len({k[0] for k in keys})
    axes = sorted({k[1] for k in keys})
    print()
    print(
        f"   {a_name}  vs  {b_name}      {n_cond} conditions x {len(axes)} axes "
        f"({', '.join(axes)}) = {len(keys)} rows      messages={policy_label}"
    )
    print(f"   {'':<26} {'':>13} {'DELIVERABLE':>13} {'':>9}   {'':>13} {'PASS-0':>13}")
    print(f"   {'stratum':<26} {'base':>13} {'arm':>13} {'':>9}   {'base':>13} {'arm':>13} "
          f"{'':>9}  better")
    print("   " + "-" * 112)
    n_ex = sum(1 for k in keys if not is_g00(k))
    line("ALL (g00 excluded)", lambda k: not is_g00(k), n_ex)
    def _in(k, s):
        return stratum(k[0]) == s and not is_g00(k)

    for s in (("stranded", "capture ON"), ("stranded", "capture OFF"),
              ("unstranded", "capture ON"), ("unstranded", "capture OFF")):
        def pred(k, s=s):
            return _in(k, s)
        # ⚠ ONE character, so the stamp cannot push the row out of the `{:<26}` column.
        tag = " ⛔" if s == DEFERRED_STRATUM else ""
        line(" x ".join(s) + tag, pred, sum(1 for k in keys if pred(k)))
    print("   ⛔ = DEFERRED for 0.8.0 (owner, 2026-08-14) — still measured and still reported;")
    print("      just not a development target. Rank on the three IN-SCOPE strata.")
    print("   " + "-" * 112)
    for ax in axes:
        def pred(k, ax=ax):
            return k[1] == ax and not is_g00(k)
        line(f"axis: {ax}", pred, sum(1 for k in keys if pred(k)))
    print("   " + "-" * 112)
    line("⛔ g00 ZERO CONTROL", is_g00, sum(1 for k in keys if is_g00(k)))

    # the library-level deliverable at the two controls, where truth is a constant
    print()
    print("   ⭐ library f_gdna at the CONTROLS (truth is a constant, so every deviation is signed)")
    print(f"      {'condition':<44} {'truth':>8} {'base':>8} {'arm':>8}")
    for cond in sorted({k[0] for k in keys}):
        if "_g00_" in cond or "_g98_" in cond:
            # ⚠ `library_f_gdna` is per CONDITION, so either axis carries it; take whichever is present
            # rather than hard-coding `region`, which an axis-sharded arm file need not have.
            k = (cond, next(ax for ax in axes if (cond, ax) in keys))
            print(f"      {cond:<44} {A[k]['library_f_gdna_truth']:>8.4f} "
                  f"{A[k]['library_f_gdna_final']:>8.4f} {B[k]['library_f_gdna_final']:>8.4f}")

    # the worst regressions and the best wins, on the deliverable
    delta = sorted(
        ((B[k]["abs_err_all_final"] - A[k]["abs_err_all_final"], k) for k in keys if not is_g00(k)),
    )
    # ⛔ HALF THE ROWS AT MOST, so the two lists cannot overlap. A fixed 6-and-6 is fine on the full
    # panel (24 non-g00 rows) and prints the SAME row as both a win and a regression on a shard.
    top = min(6, len(delta) // 2)
    print()
    print(f"   ⭐ the {top} biggest WINS and the {top} biggest REGRESSIONS on the deliverable (fragments)")
    for tag, rows in (("win", delta[:top]), ("regression", delta[len(delta) - top:][::-1])):
        for d, k in rows:
            base = A[k]["abs_err_all_final"]
            print(f"      {tag:<11} {k[0]:<44} {k[1]:<5} {base:>12,.0f} -> "
                  f"{B[k]['abs_err_all_final']:>12,.0f}  {d / base * 100 if base else 0:>+8.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
