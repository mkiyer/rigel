"""Exhaustive before/after net-flow analysis for the precision-state message fix (Phase 2).

Compares the post-EM 3-pool net fragment flow (`net_flow_3pool_per_condition.tsv`) AFTER the honest-send fix
vs the shipped BEFORE baseline (`.before_precision`). Reports: the full 16-condition table, regime rollups, and
an exhaustive residual-error budget (where the remaining error lives, by pool and direction).

Error axes (per condition):
  gdna_surplus   = gdna_assigned − gdna_true               (deconvolution over/under-call; for zero-gDNA = phantom)
  gdna↔rna leak  = net_gdna_to_nrna + net_gdna_to_mrna      (+ = gDNA→RNA leak, − = RNA→gDNA siphon)
  nrna↔mrna      = net_nrna_to_mrna                         (the EM maturity/isoform axis; calibration-adjacent)

    python scripts/debug/precision_benchmark_report.py [suite_dir]
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

SUITE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
BASE = sys.argv[2] if len(sys.argv) > 2 else "before_precision"  # baseline suffix (e.g. phase2 / before_precision)
NEW = SUITE / "net_flow_3pool_per_condition.tsv"
OLD = SUITE / f"net_flow_3pool_per_condition.tsv.{BASE}"


def load(p):
    df = pd.read_csv(p, sep="\t").set_index("condition")
    df["leak"] = df["net_gdna_to_nrna"] + df["net_gdna_to_mrna"]  # gDNA→RNA (+) / siphon (−)
    return df


def tag(r):
    cap = "CAP" if r["capture"] == "on" else "off"
    return f"{r['gdna_label']}/ss{r['ss']}/{r['nrna']}/{cap}"


def main():
    new, old = load(NEW), load(OLD)
    common = [c for c in new.index if c in old.index]
    rows = []
    for c in common:
        n, o = new.loc[c], old.loc[c]
        rows.append({
            "cond": tag(n), "gdna": n["gdna_label"], "cap": n["capture"], "ss": n["ss"], "nrna": n["nrna"],
            "gdna_true": n["gdna_true"],
            "surplus_old": o["gdna_net_surplus"], "surplus_new": n["gdna_net_surplus"],
            "leak_old": o["leak"], "leak_new": n["leak"],
            "nm_old": o["net_nrna_to_mrna"], "nm_new": n["net_nrna_to_mrna"],
        })
    df = pd.DataFrame(rows)

    print("=" * 110)
    print(f"PRECISION-STATE FIX — post-EM net fragment flow, BEFORE (.{BASE}) → AFTER (current)")
    print("=" * 110)
    print(f"\n{'condition':>26} {'gdna_true':>10} | {'gDNA surplus  old→new':>26} | "
          f"{'gDNA↔RNA leak  old→new':>28} | {'mark':>6}")
    for _, r in df.sort_values(["gdna", "cap", "ss", "nrna"]).iterrows():
        dl = r["leak_new"] - r["leak_old"]
        mark = "leak↓" if dl < -1000 else "leak↑" if dl > 1000 else "≈"
        print(f"  {r['cond']:>24} {r['gdna_true']:>10,.0f} | "
              f"{r['surplus_old']:>+12,.0f} → {r['surplus_new']:>+12,.0f} | "
              f"{r['leak_old']:>+13,.0f} → {r['leak_new']:>+13,.0f} | {mark:>6}")

    def block(title, mask):
        sub = df[mask]
        if sub.empty:
            return
        lo, ln = sub["leak_old"].abs().sum(), sub["leak_new"].abs().sum()
        so, sn = sub["surplus_old"].abs().sum(), sub["surplus_new"].abs().sum()
        d = 100 * (ln - lo) / max(lo, 1)
        print(f"  {title:<40} Σ|leak| {lo:>12,.0f} → {ln:>12,.0f} ({d:>+5.0f}%)   "
              f"Σ|surplus| {so:>11,.0f} → {sn:>11,.0f}")

    print("\n" + "-" * 110)
    print("REGIME ROLLUPS (Σ|gDNA↔RNA net flow| and Σ|gDNA surplus|; lower = better)")
    print("-" * 110)
    block("CAPTURE-ON, gDNA present (FLAGSHIP)", (df["cap"] == "on") & (df["gdna"] == "gdna300"))
    block("capture-off, gDNA present", (df["cap"] == "off") & (df["gdna"] == "gdna300"))
    block("zero-gDNA (pure false-positive test)", df["gdna"] == "none")
    block("  └ zero-gDNA + nascent (nrna=rnd)", (df["gdna"] == "none") & (df["nrna"] == "rnd"))
    block("  └ zero-gDNA, no nascent", (df["gdna"] == "none") & (df["nrna"] == "none"))
    block("stranded (ss0.99)", df["ss"] == 0.99)
    block("unstranded (ss0.50)", df["ss"] == 0.50)

    print("\n" + "-" * 110)
    print("RESIDUAL ERROR BUDGET (AFTER): the largest remaining |gDNA↔RNA net flow|, ranked")
    print("-" * 110)
    res = df.assign(absleak=df["leak_new"].abs()).sort_values("absleak", ascending=False)
    tot = df["leak_new"].abs().sum()
    cum = 0.0
    for _, r in res.head(8).iterrows():
        cum += abs(r["leak_new"])
        kind = "siphon RNA→gDNA" if r["leak_new"] < 0 else "leak gDNA→RNA"
        print(f"  {r['cond']:>24}  leak {r['leak_new']:>+13,.0f}  ({kind:>16})  "
              f"cum {100*cum/max(tot,1):>4.0f}% of Σ|leak|")
    print(f"\n  TOTAL Σ|gDNA↔RNA net flow|:  BEFORE {df['leak_old'].abs().sum():>13,.0f}  →  "
          f"AFTER {tot:>13,.0f}  ({100*(tot-df['leak_old'].abs().sum())/max(df['leak_old'].abs().sum(),1):+.0f}%)")
    print(f"  TOTAL Σ|gDNA surplus|:        BEFORE {df['surplus_old'].abs().sum():>13,.0f}  →  "
          f"AFTER {df['surplus_new'].abs().sum():>13,.0f}")
    # nascent↔mature axis (EM-side, calibration-adjacent) — report so we separate calibration from EM error
    print(f"  Σ|nascent↔mature net flow|:   BEFORE {df['nm_old'].abs().sum():>13,.0f}  →  "
          f"AFTER {df['nm_new'].abs().sum():>13,.0f}   (EM maturity axis)")


if __name__ == "__main__":
    main()
