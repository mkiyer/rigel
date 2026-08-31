"""⭐⭐⭐ **ARE THE FOUR gDNA FRAGMENT-LENGTH POOLS ACTUALLY PURE gDNA, AND WHAT DOES THE SHIPPED LENGTH
MODEL SAY AGAINST TRUTH? — no solver, no model, straight off the origin-split oracle.**

`fl.py` asserts "Every pool is PURE BY CONSTRUCTION" and `scan_payload.py` labels two of them `pure gDNA`.
This measures that assertion. Each oracle condition holds a SEPARATE full scan per origin (`gdna`, `mrna`,
`nrna`), so every pool's composition and every component's length distribution are directly observable.

Two panels per row:

* **COMPOSITION** — per pool: the gDNA / nascent / mature fragment counts, each component's mean length,
  and the length bias the mixture imposes on that pool.
* **SHIPPED vs TRUTH** — `TRUE` (the gDNA partition's own lengths, what the model should estimate),
  `POOLED` (the four pools as they actually are) and `SHIPPED` (`build_fl_models(...).gdna_pmf`, i.e.
  opportunity-divided AND empirical-Bayes shrunk). ⭐ **Read `pool−true` against `ship−pool`**: the first is
  what CONTAMINATION costs, the second what the OPPORTUNITY DIVISOR and the SHRINKAGE cost, and confusing
  them has already produced one wrong attribution.

⛔⛔ **RUN IT ON A PANEL WHERE THE TWO COMPONENTS' FRAGMENT LENGTHS DIFFER, OR IT MEASURES NOTHING.** The
bias is `RNA_share × (len_RNA − len_gDNA)`, and the ladder and test chromosome give both components EQUAL
lengths deliberately (a forcing function for the EM). On those panels a 95 %-contaminated pool reads under
a bp of error; on the fl-gap arms the same pool reads over a hundred. That is why this defect shipped.

⚠ Needs the ORACLE cache with its origin partitions (`panel.py cache`), not just a scan cache.

    python scripts/design/fl_pool_purity.py --panel <scenarios dir> --index <index dir> [--conditions ...]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from rigel.calibration.fl import build_fl_models
from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
from rigel.calibration.sj_opportunity import crossing_probability_from_index
from rigel.index import TranscriptIndex
from rigel.scan_cache import read_scan_cache
from rigel.scan_payload import (
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)

#: The four pools `fl.build_fl_models` sums into the gDNA length model, in payload order.
GDNA_POOLS = (POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC, POOL_DNA_INTRON_EXON, POOL_DNA_INTERGENIC_EXON)
POOL_NAMES = {
    POOL_DNA_INTERGENIC: "0 intergenic contained",
    POOL_DNA_INTRONIC: "1 intronic contained",
    POOL_DNA_INTRON_EXON: "2 intron|exon crossing",
    POOL_DNA_INTERGENIC_EXON: "3 intergenic|exon crossing",
    POOL_RNA_SPLICED: "4 spliced (RNA)",
}


def mean_length(hist) -> float:
    """Mean of a length histogram, or ``nan`` when the component is absent (a real state, not a stub)."""
    v = np.asarray(hist, dtype=np.float64)
    n = v.sum()
    return float((v * np.arange(v.shape[0])).sum() / n) if n > 0 else float("nan")


def _pools(panel: Path, index, cond: str, part: str):
    d = panel / "oracle_cache" / cond / part
    if not d.exists():
        return None
    return np.asarray(read_scan_cache(d, index).payload.pool_lengths, dtype=np.float64)


def composition(panel: Path, index, cond: str) -> None:
    """Per pool: who is actually in it, and what the mixture does to that pool's mean length."""
    part = {p: _pools(panel, index, cond, p) for p in ("gdna", "mrna", "nrna")}
    if any(v is None for v in part.values()):
        print(f"══ {cond}: no origin partitions — run `panel.py cache`")
        return
    print(f"\n══ {cond}")
    print(
        f"   {'pool':<26}{'gDNA':>10}{'nascent':>10}{'mature':>9} | "
        f"{'RNA share':>10}{'len gDNA':>10}{'len RNA':>10}{'POOLED':>9}{'bias':>8}"
    )
    for p in GDNA_POOLS:
        g, n, m = part["gdna"][p], part["nrna"][p], part["mrna"][p]
        tot = g.sum() + n.sum() + m.sum()
        share = (n.sum() + m.sum()) / tot if tot > 0 else 0.0
        pooled = mean_length(g + n + m)
        print(
            f"   {POOL_NAMES[p]:<26}{g.sum():>10,.0f}{n.sum():>10,.0f}{m.sum():>9,.0f} | "
            f"{100 * share:>9.1f}%{mean_length(g):>10.1f}{mean_length(n + m):>10.1f}"
            f"{pooled:>9.1f}{pooled - mean_length(g):>+8.1f}"
        )
    g_all = sum(part["gdna"][p] for p in GDNA_POOLS)
    all_all = sum(part["gdna"][p] + part["nrna"][p] + part["mrna"][p] for p in GDNA_POOLS)
    rna_all = all_all - g_all
    print(
        f"   {'ALL FOUR gDNA POOLS':<26}{g_all.sum():>10,.0f}{rna_all.sum():>19,.0f} | "
        f"{100 * rna_all.sum() / max(all_all.sum(), 1):>9.1f}%{mean_length(g_all):>10.1f}"
        f"{mean_length(rna_all):>10.1f}{mean_length(all_all):>9.1f}"
        f"{mean_length(all_all) - mean_length(g_all):>+8.1f}"
    )
    sp = part["gdna"][POOL_RNA_SPLICED].sum()
    tot_sp = sp + part["nrna"][POOL_RNA_SPLICED].sum() + part["mrna"][POOL_RNA_SPLICED].sum()
    print(f"   {'4 spliced (the RNA pool)':<26}gDNA share {100 * sp / max(tot_sp, 1):.2f}%  ⭐ the one pool that IS pure")


def shipped_vs_truth(panel: Path, index, conds) -> None:
    """`TRUE` / `POOLED` / `SHIPPED`, so contamination and the divisor+shrinkage are attributed apart."""
    print(
        f"\n{'condition':<44}{'TRUE':>8}{'POOLED':>9}{'SHIPPED':>9}{'GLOBAL':>8} | "
        f"{'ship−true':>10}{'pool−true':>10}{'ship−pool':>10}{'n_gdna':>10}"
    )
    for cond in conds:
        gd = panel / "oracle_cache" / cond / "gdna"
        if not gd.exists():
            continue
        main = read_scan_cache(panel / "oracle_cache" / cond / "_main", index)
        true = np.asarray(read_scan_cache(gd, index).payload.deposited_lengths, dtype=np.float64)
        fl = build_fl_models(
            main.payload,
            sj_opportunity=crossing_probability_from_index(index, int(main.payload.max_length)),
            gdna_opportunity=gdna_opportunity_from_index(index, int(main.payload.max_length)),
        )
        pooled = np.asarray(main.payload.pool_lengths, dtype=np.float64)[list(GDNA_POOLS)].sum(axis=0)
        t, p, s, g = mean_length(true), mean_length(pooled), mean_length(fl.gdna_pmf), mean_length(fl.global_pmf)
        print(
            f"{cond:<44}{t:>8.1f}{p:>9.1f}{s:>9.1f}{g:>8.1f} | "
            f"{s - t:>+10.1f}{p - t:>+10.1f}{s - p:>+10.1f}{fl.n_gdna:>10,.0f}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, type=Path, help="a scenarios dir holding oracle_cache/")
    ap.add_argument("--index", required=True, help="the index the panel was scanned against")
    ap.add_argument("--conditions", nargs="*", default=None, help="default: every certified condition")
    ap.add_argument("--no-composition", action="store_true", help="the shipped-vs-truth table only")
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    conds = args.conditions or sorted(
        d.name for d in (args.panel / "oracle_cache").iterdir() if (d / "_main").exists()
    )
    if not args.no_composition:
        for cond in conds:
            composition(args.panel, index, cond)
    shipped_vs_truth(args.panel, index, conds)


if __name__ == "__main__":
    main()
