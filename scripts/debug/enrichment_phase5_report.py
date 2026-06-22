"""Phase-5 before/after: post-EM net fragment-flow, enrichment-aware calibration vs the pre-enrichment
(.pre_enrichment) baseline, across the 16-condition quick_3to1_5mb suite.

Headline metric per condition: gdna_net_surplus (gDNA assigned − true; negative = gDNA→RNA leak) and the
gDNA→RNA net flow. Confirms (a) the capture conditions' leak drops, (b) the off-capture / no-gDNA controls are
unchanged.

    python scripts/debug/enrichment_phase5_report.py [suite_dir]
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

SUITE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
NEW = SUITE / "net_flow_3pool_per_condition.tsv"
OLD = SUITE / "net_flow_3pool_per_condition.tsv.pre_enrichment"


def main():
    new = pd.read_csv(NEW, sep="\t").set_index("condition")
    old = pd.read_csv(OLD, sep="\t").set_index("condition")
    common = [c for c in new.index if c in old.index]

    def short(c):
        r = new.loc[c]
        return f"{r['gdna_label']}/ss{r['ss']}/{r['nrna']}/{'CAP' if r['capture']=='on' else 'off'}"

    rows = []
    for c in common:
        n, o = new.loc[c], old.loc[c]
        leak_new = float(n["net_gdna_to_nrna"]) + float(n["net_gdna_to_mrna"])  # gDNA→RNA total
        leak_old = float(o["net_gdna_to_nrna"]) + float(o["net_gdna_to_mrna"])
        rows.append({
            "cond": short(c), "capture": n["capture"], "gdna": n["gdna_label"], "ss": n["ss"], "nrna": n["nrna"],
            "gdna_true": float(n["gdna_true"]),
            "surplus_old": float(o["gdna_net_surplus"]), "surplus_new": float(n["gdna_net_surplus"]),
            "leak_old": leak_old, "leak_new": leak_new,
        })
    df = pd.DataFrame(rows)

    def block(title, mask):
        sub = df[mask]
        if sub.empty:
            return
        print(f"\n{title}")
        print(f"  {'condition':>26} {'gdna_true':>10} | {'gDNA surplus old→new':>26} | {'gDNA→RNA leak old→new':>26}")
        for _, r in sub.iterrows():
            ds = r["surplus_new"] - r["surplus_old"]
            dl = r["leak_new"] - r["leak_old"]
            print(f"  {r['cond']:>26} {r['gdna_true']:>10,.0f} | "
                  f"{r['surplus_old']:>+11,.0f} → {r['surplus_new']:>+11,.0f} | "
                  f"{r['leak_old']:>+11,.0f} → {r['leak_new']:>+11,.0f}  "
                  f"({'leak↓' if dl < -500 else 'leak↑' if dl > 500 else '≈same'})")

    print("=== Phase-5 post-EM net fragment-flow: enrichment-aware vs pre-enrichment ===")
    block("CAPTURE-ON, gDNA present (the target — expect leak ↓):",
          (df["capture"] == "on") & (df["gdna"] == "gdna300"))
    block("CAPTURE-OFF, gDNA present (expect ≈ unchanged):",
          (df["capture"] == "off") & (df["gdna"] == "gdna300"))
    block("NO gDNA (pure false-positive test — expect ≈ unchanged, near 0):",
          (df["gdna"] == "none"))

    # aggregate
    cap = df[(df["capture"] == "on") & (df["gdna"] == "gdna300")]
    print("\n=== SUMMARY ===")
    print(f"  capture-on gDNA conditions: total gDNA→RNA leak  {cap['leak_old'].sum():+,.0f} → {cap['leak_new'].sum():+,.0f} "
          f"({100*(cap['leak_new'].sum()-cap['leak_old'].sum())/max(abs(cap['leak_old'].sum()),1):+.0f}%)")
    nog = df[df["gdna"] == "none"]
    print(f"  no-gDNA conditions (FP):    total gDNA→RNA leak  {nog['leak_old'].sum():+,.0f} → {nog['leak_new'].sum():+,.0f} "
          f"(must stay ≈ same)")
    off = df[(df["capture"] == "off") & (df["gdna"] == "gdna300")]
    print(f"  capture-off gDNA:           total surplus        {off['surplus_old'].sum():+,.0f} → {off['surplus_new'].sum():+,.0f} "
          f"(must stay ≈ same)")


if __name__ == "__main__":
    main()
