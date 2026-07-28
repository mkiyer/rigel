"""Bit-identity check between two arms in /tmp/pass0_oracle_bench.tsv (exact string compare on every field)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

a, b = sys.argv[1], sys.argv[2]
rows: dict = {}
with Path("/tmp/pass0_oracle_bench.tsv").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.setdefault(r["arm"], {})[r["cond"]] = r
A, B = rows[a], rows[b]
conds = sorted(set(A) & set(B))
ident = 0
for c in conds:
    same = all(A[c][k] == B[c][k] for k in A[c] if k != "arm")
    ident += same
    if not same:
        d = [k for k in A[c] if k != "arm" and A[c][k] != B[c][k]]
        print(f"  DIFF {c}: {[(k, A[c][k], B[c][k]) for k in d]}")
print(f"{a} vs {b}: BIT-IDENTICAL on {ident}/{len(conds)} conditions")
