#!/usr/bin/env python
"""Are the four gDNA fragment-length pools actually pure gDNA, and what does the shipped length model
say against truth? No solver and no model: straight off the origin-split oracle.

`fl.py` treats every pool as pure by construction and `scan_payload.py` labels two of them pure gDNA;
this measures that claim. Each oracle condition holds a separate full scan per origin (`gdna`, `mrna`,
`nrna`), so every pool's composition and every component's length distribution are directly observable.
Two tables per condition. Composition: per pool, the gDNA / nascent / mature fragment counts, each
component's mean length, and the length bias the mixture imposes on that pool. Shipped against truth:
`TRUE` (the gDNA partition's own lengths, what the model should estimate), `POOLED` (the four pools as
they actually are), `SHIPPED` (`build_fl_models(...).gdna_pmf`, opportunity-divided and shrunk) and
`GLOBAL` (the unconditional anchor the shrinkage pulls toward, `deposited_lengths`), so that `pool-true`
(what contamination costs) and `ship-pool` (what the opportunity divisor and the shrinkage cost) are
attributed apart and never confused. Everything is read in the DRAINED frame the truth is certified in:
the whole drained as production drains it (`scan_cache.calibration_inputs`, which also builds `SHIPPED`
exactly as production builds it, with the two-pool contrast), and the three partitions lifted by
replaying the whole's choices so they sum to the drained whole (`OracleTruth.from_cached_parts`, whose
sum-to-full gate is the lift's own identity check). Until 2026-09-13 it priced pass one's model, which
production never builds. It measures nothing on an equal-length panel: the bias is the RNA share times
the length gap, and the ladder and the test chromosome give both components equal lengths deliberately,
so run it only where the two components' fragment lengths differ. It needs the oracle cache with its
origin partitions (`panel.py cache`), not just a scan cache.

Usage::

    python scripts/design/fl_pool_purity.py --panel <scenarios dir> --index <index dir>
    python scripts/design/fl_pool_purity.py --panel DIR --index DIR --conditions <cond> ...
    python scripts/design/fl_pool_purity.py --panel DIR --index DIR --no-composition   # the second table only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from rigel.index import TranscriptIndex
from rigel.scan_cache import calibration_inputs, read_scan_cache
from rigel.scan_payload import (
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402

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


def load(panel: Path, index, cond: str):
    """One condition in the DRAINED frame: production's calibration inputs for the whole (the payload and
    the two length models) and the three origin partitions lifted into the same frame; ``None`` when the
    condition has no origin partitions."""
    root = panel / "oracle_cache" / cond
    if not all((root / k).exists() for k in ORIGINS):
        return None
    lift: dict = {}
    kw = calibration_inputs(read_scan_cache(root / "_main", index), index, lift_out=lift)
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    return kw, OracleTruth.from_cached_parts(kw["payload"], parts, lift).parts


def composition(cond: str, parts: dict) -> None:
    """Per pool: who is actually in it, and what the mixture does to that pool's mean length."""
    part = {p: np.asarray(parts[p].pool_lengths, dtype=np.float64) for p in ORIGINS}
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


def shipped_row(cond: str, kw: dict, parts: dict) -> tuple:
    """`TRUE` / `POOLED` / `SHIPPED` / `GLOBAL` mean lengths and the pools' fragment count, one condition."""
    payload = kw["payload"]
    true = np.asarray(parts["gdna"].deposited_lengths, dtype=np.float64)
    pooled = np.asarray(payload.pool_lengths, dtype=np.float64)[list(GDNA_POOLS)].sum(axis=0)
    return (
        cond,
        mean_length(true),
        mean_length(pooled),
        mean_length(kw["gdna_fl_pmf"]),
        mean_length(payload.deposited_lengths),
        float(pooled.sum()),
    )


def shipped_vs_truth(rows: list[tuple]) -> None:
    """The second table, so contamination and the divisor+shrinkage are attributed apart."""
    print(
        f"\n{'condition':<44}{'TRUE':>8}{'POOLED':>9}{'SHIPPED':>9}{'GLOBAL':>8} | "
        f"{'ship−true':>10}{'pool−true':>10}{'ship−pool':>10}{'n_pooled':>10}"
    )
    for cond, t, p, s, g, n in rows:
        print(
            f"{cond:<44}{t:>8.1f}{p:>9.1f}{s:>9.1f}{g:>8.1f} | "
            f"{s - t:>+10.1f}{p - t:>+10.1f}{s - p:>+10.1f}{n:>10,.0f}"
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
    rows = []
    for cond in conds:
        loaded = load(args.panel, index, cond)
        if loaded is None:
            print(f"══ {cond}: no origin partitions — run `panel.py cache`")
            continue
        kw, parts = loaded
        if not args.no_composition:
            composition(cond, parts)
        rows.append(shipped_row(cond, kw, parts))
    shipped_vs_truth(rows)


if __name__ == "__main__":
    main()
