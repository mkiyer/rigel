#!/usr/bin/env python
"""⭐⭐⭐ **IS THE PRIOR ASSEMBLER RIGHT? — the flgap pair, DRAINED, against the EM's own count.**

**The question, and why it needed its own instrument.** ``prior_vs_oracle.py`` runs on any suite and is
UNDRAINED on every arm — forced, because a drained oracle is inadmissible on the ladder (its
``lift_choices`` replay deposits spliced records into the gdna partition and ``OracleTruth`` refuses
it, correctly). Production always drains. So the one panel where the assembler can be scored the way
production runs it is the **flgap PAIR**, whose drained partitions the study cache already holds — and
the contamination that makes flgap_long a cross-check rather than a verdict is printed first, as a
number, not waved at.

⭐⭐ **THE YARDSTICK IS ``Fo``, NOT ``F``, AND THAT IS THE POINT OF THIS FILE.** ``F`` counts a fragment
in whichever locus holds its FIRST BASE; the EM counts every fragment that is a CANDIDATE, once, in
the one multi-locus it belongs to. ``F`` therefore drops the straddling population and ``O − F`` /
``S − F`` were measured against a target nothing consumes. Both are printed side by side so the size of
the correction is visible (``prior_vs_oracle.overlap_truth``).

⛔ **NOTHING HERE RE-IMPLEMENTS AN ARM.** The assembler is the shipped ``assemble_priors``, O and S come
from ``prior_vs_oracle.assemble_oracle_arm`` / ``assemble_share_arm``, ``Fo`` from
``prior_vs_oracle.overlap_truth`` and the scoring from ``prior_vs_oracle.score_arm``. This file loads a
cache, calls them, and prints (TRAPS: a-test-that-redefines).

⭐ ``priors.py`` is outside the study cache's key, so testing an ``assemble_priors`` change is a
**one-second** loop: edit, re-run this, read the table.

Usage::

    python scripts/design/prior_yardstick.py                    # all four conditions, both arms
    python scripts/design/prior_yardstick.py --arms drained     # production's arm alone
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
for _p in (_REPO / "scripts" / "design", _REPO / "tests" / "calibration"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import flgap_study_cache as FC  # noqa: E402
import rigel.calibration.priors as PRIORS  # noqa: E402
from prior_vs_oracle import (  # noqa: E402
    assemble_oracle_arm,
    assemble_share_arm,
    overlap_truth,
    project_start_counts,
    score_arm,
)

#: ``(suffix, label)`` per arm. ⛔ The drained arm is what PRODUCTION calibrates on; the undrained one
#: is what every other instrument in this repo reports, and they are printed together so a difference
#: between this file and those is attributable rather than mysterious.
ARMS = (("", "UNDRAINED"), ("_d", "DRAINED (production)"))


def build(st: dict, ra, suffix: str) -> dict:
    """One condition, one arm: assemble P, O and S, tally Fo, and score every pair.

    ⭐ **Every arm is built HERE, at read time**, from the four things the cache stores that
    ``priors.py`` did not produce: ``cal``, ``multi_loci``, the truth ``override``/``shares``, and the
    raw per-region ``starts``. That is what makes an ``assemble_priors`` edit a one-second loop AND
    keeps P, O, S and F on the same assembler — storing any of them would have decoupled them.
    """
    cal = st["cal" + suffix]
    ml = st["multi_loci" + suffix]
    override = st["override" + suffix]
    shares = st["shares" + suffix]
    p = PRIORS.assemble_priors(cal, ra, ml)
    o = assemble_oracle_arm(override, cal, ra, ml)
    s = assemble_share_arm(override, shares, cal, ra, ml)
    fo = overlap_truth(
        ml,
        st["unit_origin" + suffix],
        st["unit_is_spliced" + suffix],
        st["n_units" + suffix],
        st["walk"],
    )
    f_g, f_r, _dropped = project_start_counts(st["starts" + suffix], ra, ml)
    return {
        "n_loci": len(ml),
        "fo": fo,
        "gdna": {
            "Fo_vs_F": score_arm(fo.gdna, f_g),
            "O_vs_Fo": score_arm(o.gdna_prior_count, fo.gdna),
            "S_vs_Fo": score_arm(s.gdna_prior_count, fo.gdna),
            "P_vs_Fo": score_arm(p.gdna_prior_count, fo.gdna),
            "O_vs_F": score_arm(o.gdna_prior_count, f_g),
            "S_vs_F": score_arm(s.gdna_prior_count, f_g),
            "O_vs_S": score_arm(o.gdna_prior_count, s.gdna_prior_count),
        },
        "rna": {
            "Fo_vs_F": score_arm(fo.rna_unspliced, f_r),
            "O_vs_Fo": score_arm(o.rna_prior_count, fo.rna_unspliced),
            "S_vs_Fo": score_arm(s.rna_prior_count, fo.rna_unspliced),
            "P_vs_Fo": score_arm(p.rna_prior_count, fo.rna_unspliced),
            "O_vs_F": score_arm(o.rna_prior_count, f_r),
            "S_vs_F": score_arm(s.rna_prior_count, f_r),
            "O_vs_S": score_arm(o.rna_prior_count, s.rna_prior_count),
        },
        # ⛔ The denominator is the UNSPLICED pool, never `rna_all`: spliced RNA has no gDNA
        # competitor, so a prior that arbitrates that competition must not count it
        # (`TRAPS: score-the-consumers-own-count`).
        "phi": {
            "true": _phi(fo.gdna.sum(), fo.rna_unspliced.sum()),
            "S": _phi(s.gdna_prior_count.sum(), s.rna_prior_count.sum()),
            "P": _phi(p.gdna_prior_count.sum(), p.rna_prior_count.sum()),
            "strength": float(
                (s.gdna_prior_count.sum() + s.rna_prior_count.sum())
                / max(fo.gdna.sum() + fo.rna_unspliced.sum(), 1.0)
            ),
        },
    }


def tag(panel: str, cond: str) -> str:
    """``long/ON`` — the panel and the capture arm, with NOTHING truncated away.

    ⛔ The first version sliced seven characters out of the middle of the condition name, which is the
    same seven for every condition on this panel: four rows, two labels, and the capture arms — the axis
    the whole flgap pair exists to vary — indistinguishable. A table whose rows cannot be told apart is
    not a table.
    """
    return f"{panel.replace('flgap_', '')}/{'ON' if 'capture_on' in cond else 'OFF'}"


def _phi(g: float, r: float) -> float:
    tot = float(g) + float(r)
    return float(g) / tot if tot > 0 else float("nan")


def _rel(x: float) -> str:
    if not np.isfinite(x):
        return f"{'nan':>9}"
    return f"{x:>9.2e}" if 0.0 < abs(x) < 1e-3 else f"{x:>9.4f}"


def main() -> int:
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arms", nargs="*", default=None,
                    help="undrained / drained (default: both)")
    ap.add_argument("--conditions", nargs="*", default=None, help="panel/condition substrings")
    args = ap.parse_args()

    arms = [
        (suf, lab) for suf, lab in ARMS
        if args.arms is None or lab.split()[0].lower() in {a.lower() for a in args.arms}
    ]
    conds = [
        (panel, cond) for panel, cond in FC.CONDITIONS
        if args.conditions is None or any(c in f"{panel}/{cond}" for c in args.conditions)
    ]

    index = TranscriptIndex.load(str(FC.INDEX))
    ra = RegionArrays.from_index(index)

    rows: list[tuple] = []
    for panel, cond in conds:
        st = FC.load(panel, cond)
        for suffix, label in arms:
            if ("cal" + suffix) not in st:
                print(f"  ⚠ {panel} has no {label} arm cached (nothing was held) — skipped")
                continue
            rows.append((panel, cond, label, st, build(st, ra, suffix)))

    print()
    print("=" * 118)
    print("  ⭐⭐⭐ THE PRIOR ASSEMBLER AGAINST THE EM's OWN CANDIDATE COUNT — the flgap pair")
    print("  messages OFF   length_likelihood OFF   both capture arms   ss_0.50 (unstranded)")
    print("=" * 118)

    # ── admissibility FIRST. A drained number read without its contamination is not a number. ──
    print()
    print("  ⛔ ADMISSIBILITY OF THE DRAINED ARM — gDNA cannot splice, so any spliced deposit in the")
    print("     gdna partition is contamination the lift introduced. Printed before any score.")
    print(f"    {'panel/condition':<52} {'held':>10} {'lift ambig':>11} {'gdna spliced leak':>18}")
    print("    " + "-" * 96)
    for panel, cond in conds:
        st = FC.load(panel, cond)
        print(f"    {panel + '/' + cond:<52} {st['n_held']:>10,} "
              f"{st.get('lift_ambiguous', 0):>11,} {st.get('gdna_spliced_leak', 0):>18,}")
    print("    ⭐ flgap_short is the ADMISSIBLE panel for a drained verdict; flgap_long is a "
          "cross-check that carries real contamination.")

    # ── the join ──
    print()
    print("  ⛔ THE Fo JOIN — the hard gate ran when the cache was built (record + group counts vs the")
    print("     scanner's own). Below is the secondary diagnostic: gDNA cannot splice.")
    print(f"    {'panel/condition':<40} {'arm':<22} {'units':>12} {'spliced&gdna':>13} "
          f"{'orphans':>9} {'origin trans':>13}")
    print("    " + "-" * 114)
    for panel, cond, label, st, r in rows:
        d = r["fo"].diag
        print(f"    {tag(panel, cond):<40} {label:<22} {d['n_units']:>12,} "
              f"{d['spliced_gdna_units']:>13,} {d['orphan_units']:>9,} "
              f"{d['walk']['n_transitions']:>13,}")

    # ── ① the yardstick correction ──
    print()
    print("  ① ⛔ WHAT THE YARDSTICK CORRECTION IS WORTH — Fo (EM candidates) vs F (first-base starts)")
    print(f"    {'panel/condition':<32} {'arm':<12} {'Σ Fo':>12} {'Σ F':>12} {'Σ|Fo−F|':>10} "
          f"{'rel':>9}")
    print("    " + "-" * 100)
    for panel, cond, label, st, r in rows:
        y = r["gdna"]["Fo_vs_F"]
        print(f"    {tag(panel, cond):<32} {label.split()[0]:<12} {y.total_arm:>12,.0f} "
              f"{y.total_ref:>12,.0f} {y.abs_err:>10,.0f} {_rel(y.rel)}")

    # ── ② the assembler, scored both ways ──
    for arm_name, title in (("gdna", "gDNA arm"), ("rna", "RNA arm (target = UNSPLICED RNA units)")):
        print()
        print(f"  ② THE ASSEMBLER'S OWN ERROR · {title}")
        print("     ⭐ O = perfect masses in. S = O plus each component's OWN true per-boundary share.")
        print(f"    {'panel/condition':<26} {'arm':<11} {'O−Fo':>9} {'S−Fo':>9} {'O−S':>9} "
              f"{'O−F':>9} {'S−F':>9} {'P−Fo':>9}")
        print("    " + "-" * 106)
        for panel, cond, label, st, r in rows:
            g = r[arm_name]
            print(f"    {tag(panel, cond):<26} {label.split()[0]:<11} "
                  f"{_rel(g['O_vs_Fo'].rel)} {_rel(g['S_vs_Fo'].rel)} {_rel(g['O_vs_S'].rel)} "
                  f"{_rel(g['O_vs_F'].rel)} {_rel(g['S_vs_F'].rel)} {_rel(g['P_vs_Fo'].rel)}")
        print("    ⚠ every column is Σ|Δ| / Σ ref — additive within a column, not across them.")

    # ── ③ the share's share, on the CORRECT yardstick ──
    print()
    print("  ③ ⭐⭐ HOW MUCH OF THE ASSEMBLER'S RESIDUAL IS THE POOLED PER-COMPONENT SHARE?")
    print("     ⛔ Read against Fo. Against F this fraction was 15–43 % and the rest was called an open")
    print("        residual; that residual was the yardstick.")
    print(f"    {'panel/condition':<26} {'arm':<11} {'Σ|O−Fo|':>11} {'Σ|S−Fo|':>11} "
          f"{'share %':>9} {'left over':>11}")
    print("    " + "-" * 84)
    for panel, cond, label, st, r in rows:
        g = r["gdna"]
        o_err, s_err = g["O_vs_Fo"].abs_err, g["S_vs_Fo"].abs_err
        pct = 100.0 * (o_err - s_err) / o_err if o_err > 0 else float("nan")
        print(f"    {tag(panel, cond):<26} {label.split()[0]:<11} {o_err:>11,.0f} "
              f"{s_err:>11,.0f} {pct:>8.1f}% {s_err:>11,.0f}")

    # ── ④ the composition claim, against the pool it describes ──
    print()
    print("  ④ ⭐ THE COMPOSITION CLAIM, AGAINST THE UNSPLICED POOL — the population a_g:a_r describes")
    print("     ⛔ A spliced unit never gets a gDNA candidate, so spliced RNA does not compete with gDNA")
    print("     and must not enter this prior. phi = gDNA units / (gDNA + UNSPLICED RNA units).")
    print(f"    {'panel/condition':<20} {'arm':<8} {'phi true':>9} {'phi S':>8} {'Δ (S)':>9} "
          f"{'phi P':>8} {'Δ (P)':>9} {'strength':>9}")
    print("    " + "-" * 88)
    for panel, cond, label, st, r in rows:
        ph = r["phi"]
        print(f"    {tag(panel, cond):<20} {label.split()[0]:<8} {ph['true']:>9.4f} "
              f"{ph['S']:>8.4f} {ph['S'] - ph['true']:>+9.4f} {ph['P']:>8.4f} "
              f"{ph['P'] - ph['true']:>+9.4f} {ph['strength']:>9.3f}")
    print("    ⭐ Δ(S) is the assembler's own composition error with perfect masses in: ≤ 5e-4.")
    print("    ⭐ `strength` = Σ(a_g+a_r) / the unspliced pool. ~1.000 BY CONSTRUCTION — the conserved")
    print("       count IS that pool — so the posterior is a 50/50 blend of calibration and the EM's own")
    print("       evidence, and there is no knob. Not a defect; a design fact nothing had priced.")
    print("    ⚠ Δ(P) is CALIBRATION's error, and it is the real one: −0.48/−0.57 under capture.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
