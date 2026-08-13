"""Verify a rebuilt index against the one it replaces: regions UNCHANGED, boundaries changed only in `reach`.

    TODO item 1   ·   Ledger: S1 (what changed and why)

⭐ WHY THIS CHECK IS THE WHOLE POINT. S1 changed `_contiguous_reaches` — the reach a contiguous boundary
reports — and **nothing else**. So a correct rebuild has an exactly predictable shape:

* `regions.feather` **byte-identical**. The partition did not move. ⚠ This is also the check that the rebuild
  used the RIGHT SOURCE: a different FASTA or GTF moves the cuts, and regions would differ immediately.
* `edges.feather` differing in the four `reach_*` columns of the **contiguous** rows and nowhere else.
  ⚠ Junction reach is deliberately unchanged — a junction boundary is only used by a molecule that spliced
  across it, so what remains either side is exonic, and `_junction_edges` stays on the exonic reach.

Anything else is a finding, not a rebuild. A "rebuild" that also moved the flags, the kinds, or the region
ids would be a different change wearing this one's clothes — and `partition_hash` would not notice, because
it covers `regions.feather` only.

    python scripts/design/verify_index_rebuild.py OLD_INDEX NEW_INDEX
"""

from __future__ import annotations

import argparse
import filecmp
import sys
from pathlib import Path

import pandas as pd

#: The columns S1 was allowed to move, and only on contiguous rows.
REACH_COLUMNS = ["reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg"]

EDGE_KIND_CONTIGUOUS = 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("old")
    ap.add_argument("new")
    args = ap.parse_args()
    old, new = Path(args.old), Path(args.new)

    failures: list[str] = []
    print(f"old  {old}\nnew  {new}\n")

    # ── regions: byte-identical, which also proves the source was right ─────────────────────────────────
    if filecmp.cmp(old / "nodes.feather", new / "nodes.feather", shallow=False):
        print("  OK    regions.feather byte-identical  (the partition did not move)")
    else:
        old_regions = pd.read_feather(old / "nodes.feather")
        new_regions = pd.read_feather(new / "nodes.feather")
        if old_regions.equals(new_regions):
            print("  OK    regions.feather differs on disk but is EQUAL as a table (compression only)")
        else:
            failures.append(
                f"regions.feather CHANGED: {len(old_regions):,} -> {len(new_regions):,} rows. The partition "
                f"moved, which S1 did not do — so this rebuild used a different FASTA or GTF."
            )
            print(f"  FAIL  regions.feather changed: {len(old_regions):,} -> {len(new_regions):,} rows")

    # ── boundaries: only the reach columns, only on contiguous rows ────────────────────────────────────────
    old_boundaries = pd.read_feather(old / "edges.feather")
    new_boundaries = pd.read_feather(new / "edges.feather")
    if len(old_boundaries) != len(new_boundaries):
        failures.append(f"edges.feather row count {len(old_boundaries):,} -> {len(new_boundaries):,}")
        print(f"  FAIL  edges.feather row count {len(old_boundaries):,} -> {len(new_boundaries):,}")
    elif list(old_boundaries.columns) != list(new_boundaries.columns):
        failures.append("edges.feather columns changed")
        print("  FAIL  edges.feather columns changed")
    else:
        contiguous = new_boundaries["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
        for column in old_boundaries.columns:
            same = old_boundaries[column].equals(new_boundaries[column])
            if column in REACH_COLUMNS:
                moved = int((old_boundaries[column].to_numpy() != new_boundaries[column].to_numpy()).sum())
                junction_moved = int(
                    (
                        old_boundaries.loc[~contiguous, column].to_numpy()
                        != new_boundaries.loc[~contiguous, column].to_numpy()
                    ).sum()
                )
                mark = "OK  " if junction_moved == 0 else "FAIL"
                print(
                    f"  {mark}  {column:<16} {moved:>10,} rows moved "
                    f"({100 * moved / len(old_boundaries):5.1f} %), of which junction rows: {junction_moved:,}"
                )
                if junction_moved:
                    failures.append(
                        f"{column}: {junction_moved:,} JUNCTION rows moved. Junction reach is exonic by "
                        f"design and S1 did not touch it."
                    )
            elif not same:
                changed = int((old_boundaries[column].to_numpy() != new_boundaries[column].to_numpy()).sum())
                failures.append(f"{column}: {changed:,} rows changed, but only reach may move")
                print(f"  FAIL  {column:<16} {changed:>10,} rows changed — NOT a reach column")

    print()
    if failures:
        print("⛔ NOT THE EXPECTED REBUILD")
        for line in failures:
            print(f"  {line}")
        sys.exit(1)
    print("✅ regions unchanged; boundaries changed only in contiguous reach — exactly what S1 did")


if __name__ == "__main__":
    main()
