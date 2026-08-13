#!/usr/bin/env python
"""⭐⭐ **ONE TOY SPEC × EVERY CACHED CONDITION × AN RNA-DENSITY LADDER, SCORED PER OBJECT.**

`toy_harness.py` answers "what does this structure do on ONE library at ONE RNA level". This answers
"where does it break across the whole space the cache can express", which is the question you ask once a
structure is the target rather than a probe. It is `toy_harness`'s own machinery — same spec, same
`run_toy`, same `object_rows`, same per-object truth — swept over two axes and aggregated.

⛔⛔ **IT MEASURES THE PRIOR-FREE PASS-0 SOLVE BY DEFAULT** (``--refit-iters 0``), because that is the
substrate the gDNA hyperprior is later fitted against: an object that is confidently wrong here anchors the
prior and corrupts everything downstream. ``--refit-iters 3`` gives the shipped solve for comparison, and
the two answer different questions — do not quote one for the other.

⭐ **The RNA density is quoted as a MULTIPLE of each donor's own gDNA density**, never absolutely. What the
solver has to resolve is the RATIO, and the left_boundaries' gDNA rates span ~100× across the ladder, so an absolute
density would put different left_boundaries at completely different true ``f_g`` on the "same" rung.

⚠ **What it cannot say.** Every cached condition is ``nrna_none``, so on any spec with an intron the
intron↔exon EDGEs have truth exactly 1.0 and the panel structurally cannot distinguish "no RNA crosses this
seam" from "no *mature* RNA crosses it". ⛔ Read `docs/TESTING.md` on the panels before concluding anything
about a seam.

⚠ **Harvesting is 30 s per donor and is deliberately not cached to disk** (`DonorGlobals`' own docstring:
a stored bundle goes stale on exactly the changes this harness exists to test). So a full 36-condition run
is ~20 min single-process. Shard it with ``--conditions``.

Usage::

    python scripts/design/toy_panel.py --spec spliced_exons
    python scripts/design/toy_panel.py --spec spliced_exons --conditions gdna_g25_ss_0.50_nrna_none_capture_off
    python scripts/design/toy_panel.py --report rows.jsonl        # re-aggregate, no re-measurement
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import os
import sys
from collections import defaultdict
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


TH = _sibling("toy_harness.py")

from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

#: ⭐ multiples of the DONOR's own gDNA density. At `m` the exon's true `f_g` is roughly `1/(1+m)`, so this
#: ladder spans ~0.91 → ~0.01 and brackets the crossover at 1x. ⚠ Read points, not tuned values.
LADDER = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)


def measure(spec_name: str, conditions: list[str], *, suite: Path, index_path: Path,
            work_dir: Path, refit_iters: int, ladder=LADDER, genome_length=0, nrna=0.0):
    """Yield one dict per (condition, density rung, chain slot)."""
    index = TranscriptIndex.load(str(index_path))
    config = dataclasses.replace(CalibrationConfig(), calib_refit_iters=int(refit_iters))
    spec = TH.SPECS[spec_name]
    if genome_length:
        spec = dataclasses.replace(spec, genome_length=int(genome_length))
    if nrna:
        spec = dataclasses.replace(spec, nrna_abundance=float(nrna))
    ebp = TH.exon_bp(spec)
    for cond in conditions:
        print(f"  harvesting {cond} …", flush=True)
        donor = TH.harvest(suite / cond, index, config=config)
        g_rate = float(donor.gdna_rate_per_base)
        for mult in ladder:
            n_rna = max(int(round(mult * g_rate * ebp)), 1)
            one = dataclasses.replace(
                spec,
                name=f"{spec.name}_m{mult:g}".replace(".", "p"),
                n_rna_fragments=n_rna,
            )
            r = TH.run_toy(one, donor, work_dir, config=config)
            for row in TH.object_rows(r):
                yield {
                    "condition": cond,
                    "capture": "on" if donor.capture_on else "off",
                    "strand": f"{donor.strand_specificity:g}",
                    "gdna_rate": g_rate,
                    "mult": mult,
                    "n_rna": n_rna,
                    "n_gdna": int(r.n_gdna_target),
                    **{k: row[k] for k in (
                        "slot", "axis", "type", "where", "bp", "n", "spliced", "junction",
                        "true_fg", "fg_loc", "pred_fg", "sd_fg", "tau", "err", "mass",
                    )},
                }


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE REPORT
# ──────────────────────────────────────────────────────────────────────────────────────────────────

#: ⭐ The chain slot's identity for aggregation: its object class plus, for a repeated class, WHICH one.
#: `intron|exon` appears twice on a two-exon gene and they are NOT interchangeable — one is the donor side
#: and one the acceptor side — so collapsing them would average away the asymmetry the graft creates.
def _key(row) -> str:
    return f"{row['type']}@{row['where']}"


def _live(rows):
    return [r for r in rows if r["mass"] > 0 and np.isfinite(r["true_fg"])]


def _mwae(rows) -> float:
    live = _live(rows)
    if not live:
        return float("nan")
    w = np.array([r["mass"] for r in live])
    d = np.array([abs(r["pred_fg"] - r["true_fg"]) for r in live])
    return float((d * w).sum() / w.sum())


def report(rows: list[dict], spec_name: str, refit_iters: int) -> None:
    conds = sorted({r["condition"] for r in rows})
    mults = sorted({r["mult"] for r in rows})
    arm = "PRIOR-FREE PASS-0" if refit_iters == 0 else f"{refit_iters}-refit (the shipped solve)"
    print()
    print("=" * 130)
    print(f"⭐⭐ {spec_name} — {len(conds)} cached conditions x {len(mults)} RNA rungs, {arm}")
    print("=" * 130)
    print("   ⛔ RNA density is a MULTIPLE of each donor's OWN gDNA density, so a rung means the same")
    print("      thing on every row: at m the exon's true f_g is roughly 1/(1+m).")

    # ── 1. the gene, per stratum x rung ───────────────────────────────────────────────────────────
    print()
    print("── 1. MASS-WEIGHTED |Δf_g| OVER THE WHOLE GENE, by stratum ─────────────────────────────")
    print(f"\n   {'stratum':<26}" + "".join(f"{'m=' + f'{m:g}':>9}" for m in mults) + f"{'mean':>9}")
    print("   " + "-" * (26 + 9 * (len(mults) + 1)))
    strata = sorted({(r["capture"], r["strand"]) for r in rows})
    for cap, ss in strata:
        cells = []
        for m in mults:
            per = [
                _mwae([r for r in rows if r["condition"] == c and r["mult"] == m])
                for c in conds
                if any(r["condition"] == c and r["capture"] == cap and r["strand"] == ss
                       for r in rows)
            ]
            per = [v for v in per if np.isfinite(v)]
            cells.append(float(np.mean(per)) if per else float("nan"))
        label = f"capture {cap:<3} ss {ss}"
        print(f"   {label:<26}" + "".join(f"{c:>9.4f}" for c in cells)
              + f"{float(np.nanmean(cells)):>9.4f}")

    # ── 2. per object, the error and who moved it ─────────────────────────────────────────────────
    print()
    print("── 2. PER OBJECT — mean |Δf_g|, its share of the gene's ERROR MASS, and loc vs pred ────")
    print("   ⭐ `loc` is the message-free local self-solve; `pred` is after the relay. `relay Δ` > 0")
    print("      means the messages moved this object AWAY from truth. ⛔ Rank by err MASS, not by mean.")
    by_obj = defaultdict(list)
    for r in rows:
        by_obj[_key(r)].append(r)
    tot_err = sum(r["err"] for r in rows)
    order = sorted(by_obj, key=lambda k: -sum(r["err"] for r in by_obj[k]))
    print(f"\n   {'object':<28} {'n':>7} {'true f_g':>9} {'pred':>7} {'|Δ|':>7} {'loc |Δ|':>8} "
          f"{'relay Δ':>8} {'err mass':>10} {'share':>7} {'sd':>6} {'tau':>9}")
    print("   " + "-" * 118)
    for k in order:
        g = by_obj[k]
        lv = _live(g)
        if not lv:
            continue
        d = np.array([abs(r["pred_fg"] - r["true_fg"]) for r in lv])
        dl = np.array([abs(r["fg_loc"] - r["true_fg"]) for r in lv])
        e = sum(r["err"] for r in g)
        print(f"   {k:<28} {np.mean([r['n'] for r in lv]):>7.0f} "
              f"{np.mean([r['true_fg'] for r in lv]):>9.4f} "
              f"{np.mean([r['pred_fg'] for r in lv]):>7.4f} {d.mean():>7.4f} {dl.mean():>8.4f} "
              f"{d.mean() - dl.mean():>+8.4f} {e:>10,.0f} {e / max(tot_err, 1):>6.1%} "
              f"{np.mean([r['sd_fg'] for r in lv]):>6.3f} "
              f"{np.mean([r['tau'] for r in lv]):>9.3g}")

    # ── 3. the sweep shape, for the objects that carry the error ──────────────────────────────────
    print()
    print("── 3. THE SWEEP SHAPE — does the answer TRACK the data? ─────────────────────────────────")
    print("   ⭐ A prediction that does not move when the truth moves is the tell (TRAPS: a-total-density-ratio).")
    for k in order[:4]:
        print(f"\n   {k}")
        print(f"      {'rung':>7} {'true':>7} {'loc':>7} {'pred':>7} {'|Δ|':>7} {'sd':>6} "
              f"{'n':>8} {'junction':>9}")
        for m in mults:
            g = _live([r for r in by_obj[k] if r["mult"] == m])
            if not g:
                continue
            print(f"      {'m=' + f'{m:g}':>7} {np.mean([r['true_fg'] for r in g]):>7.4f} "
                  f"{np.mean([r['fg_loc'] for r in g]):>7.4f} "
                  f"{np.mean([r['pred_fg'] for r in g]):>7.4f} "
                  f"{np.mean([abs(r['pred_fg'] - r['true_fg']) for r in g]):>7.4f} "
                  f"{np.mean([r['sd_fg'] for r in g]):>6.3f} "
                  f"{np.mean([r['n'] for r in g]):>8.0f} "
                  f"{np.mean([r['junction'] for r in g]):>9.0f}")

    # ── 4. confidently wrong: the population that corrupts a prior ────────────────────────────────
    print()
    print("── 4. ⛔ CONFIDENTLY WRONG — |Δf_g| against the DECLARED sd, per object ─────────────────")
    print("   ⭐ z = |pred − true| / sd(f_g). An object wrong by 3 sd anchors the hyperprior against")
    print("      correct neighbours. ⚠ sd here is sd(log f_g) as the solver stores it, so z is")
    print("      indicative rather than exact — the ranking is what matters.")
    print(f"\n   {'object':<28} {'|z|>1':>7} {'|z|>2':>7} {'|z|>3':>7} {'max |z|':>9} {'worst rung':>12}")
    print("   " + "-" * 78)
    for k in order:
        lv = [r for r in _live(by_obj[k]) if r["sd_fg"] > 0]
        if not lv:
            continue
        z = np.array([abs(r["pred_fg"] - r["true_fg"]) / r["sd_fg"] for r in lv])
        w = lv[int(np.argmax(z))]
        print(f"   {k:<28} {(z > 1).mean():>6.0%} {(z > 2).mean():>6.0%} {(z > 3).mean():>6.0%} "
              f"{z.max():>9.2f} {'m=' + f'{w['mult']:g}':>12}")


def main() -> int:
    P0 = _sibling("pass0_vs_oracle.py")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="spliced_exons")
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE.parent / "ladder")
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None, help="omit for every cached condition")
    ap.add_argument("--refit-iters", type=int, default=0, help="0 = the PRIOR-FREE pass-0 solve")
    ap.add_argument("--genome-length", type=int, default=0,
                    help="⛔ CAPTURE-ON NEEDS THIS. The gDNA budget is rate x genome_length while the "
                         "probe footprint is fixed, so a longer chromosome lets capture concentrate a "
                         "bigger budget onto the same probes and an EDGE's count grows with it. "
                         "`spliced_exons` needs ~120000 on a capture-ON donor (docs/TESTING.md 0b)")
    ap.add_argument("--nrna", type=float, default=0.0,
                    help="nascent abundance. ⭐ THE CONTROL every cached condition lacks: with "
                         "nrna_none an intron|exon EDGE's truth is exactly 1.000, so the intron-facing "
                         "identity is only testable non-trivially with this > 0")
    ap.add_argument("--work-dir", type=Path,
                    default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_toy_panel")
    ap.add_argument("--out", type=Path, default=None, help="write the per-object rows as JSONL")
    ap.add_argument("--report", type=Path, default=None, help="re-aggregate a JSONL, no measurement")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for name, spec in TH.SPECS.items():
            print(f"  {name:<18} {spec.genome_length:>7,} bp  {len(spec.genes)} gene(s)")
        return 0

    if args.report:
        rows = [json.loads(x) for x in args.report.read_text().splitlines() if x.strip()]
        report(rows, args.spec, args.refit_iters)
        return 0

    if args.spec not in TH.SPECS:
        print(f"unknown spec {args.spec!r}; --list to see them", file=sys.stderr)
        return 2
    conds = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    if not conds:
        print(f"no conditions with a sim_oracle.bam under {args.suite}", file=sys.stderr)
        return 2

    rows: list[dict] = []
    fh = args.out.open("w") if args.out else None
    try:
        for row in measure(args.spec, conds, suite=args.suite, index_path=args.index,
                           work_dir=args.work_dir, refit_iters=args.refit_iters,
                           genome_length=args.genome_length, nrna=args.nrna):
            rows.append(row)
            if fh:
                fh.write(json.dumps(row) + "\n")
                fh.flush()
    finally:
        if fh:
            fh.close()
    report(rows, args.spec, args.refit_iters)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
