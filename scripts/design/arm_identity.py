#!/usr/bin/env python
"""Is this arm byte-identical to that one? The byte-identity gate for a restructure that must be a
numerical no-op. An aggregate answers "how much did this move?"; a no-op needs "did any field move
at all?", and a sum cannot answer that because two fields moving by +x and -x cancel. This compares
every scored field of every row of two arm files (the `.jsonl` rows `quant_accuracy.py` and
`prior_vs_oracle.py` write with `--out`), keyed by (condition, axis). Both failure modes it has
shown are gated rather than trusted (`TRAPS: byte-identity-gate`): the row-key sets must be equal,
not overlapping, and two empty arms fail rather than pass vacuously; and every differing field is
named with its own max |delta| and location, with both files' mtimes printed, because a baseline
recorded in an earlier session is the thing to re-record, never believe. Identity is bit-equality
or two `nan`s (a field not produced at a condition), never a tolerance. It measures nothing and
patches nothing in `src`; exit status is 0 only if every field of every row is identical, so it
composes into a shell gate.

Usage::

    python scripts/design/arm_identity.py qa_base qa_noop            # names under $RIGEL_ARMS
    python scripts/design/arm_identity.py a.jsonl b.jsonl            # or explicit paths
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

#: where an instrument's `--out` wrote the arms. Override with $RIGEL_ARMS.
D = Path(os.environ.get("RIGEL_ARMS", Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_arms"))

#: Not compared: fields that differ by construction and would report a difference on every row —
#: ``arm`` (the arm's own name) and ``seconds`` (wall clock, an observation of the machine). Every
#: other field is a measurement and is compared. A field belongs here only if it cannot be equal
#: between two identical arms: every name added is a way for a real difference to stop being reported.
_NOT_A_MEASUREMENT = frozenset({"arm", "seconds"})


def _resolve(name: str) -> Path:
    p = Path(name)
    return p if p.suffix == ".jsonl" else D / f"{name}.jsonl"


def load(name: str) -> tuple[dict, Path]:
    p = _resolve(name)
    if not p.is_file():
        raise SystemExit(f"⛔ no such arm file: {p}")
    rows = [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
    keyed = {(r["condition"], r["axis"]): r for r in rows}
    if len(keyed) != len(rows):
        raise SystemExit(f"⛔ {p}: {len(rows)} rows collapse to {len(keyed)} keys — duplicate rows")
    return keyed, p


def _same(a, b) -> bool:
    """Identity, not closeness. Two `nan`s are the same absence; anything else is compared as bits."""
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b):
            return True
    return a == b


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    a_name, b_name = sys.argv[1], sys.argv[2]
    A, pa = load(a_name)
    B, pb = load(b_name)

    print()
    print(f"   THE BYTE-IDENTITY GATE   {a_name}   vs   {b_name}")
    for nm, p in ((a_name, pa), (b_name, pb)):
        print(f"      {nm:<18} {p}   ({len(load(nm)[0])} rows, "
              f"written {time.strftime('%Y-%m-%d %H:%M', time.localtime(p.stat().st_mtime))})")

    fails: list[str] = []

    # ── the row sets must be EQUAL: a missing row is a failure, not a skipped row ──────────────────
    only_a, only_b = sorted(set(A) - set(B)), sorted(set(B) - set(A))
    if only_a or only_b:
        for k in only_a:
            fails.append(f"row present ONLY in {a_name}: {k[0]} {k[1]}")
        for k in only_b:
            fails.append(f"row present ONLY in {b_name}: {k[0]} {k[1]}")
    keys = sorted(set(A) | set(B))
    if not keys:
        fails.append("BOTH arms are empty — this gate would otherwise pass vacuously (TRAPS: byte-identity-gate)")

    # ── the field sets must match too, PER ROW KEY: a field added or dropped is a change to what was
    # measured. Rows of different axes may legitimately carry different schemas (`quant_accuracy.py`
    # writes a `transcript` row and a `library` row per condition), so a whole-file union of field
    # names would report spurious misses; per key, a field present in A's row and absent from B's
    # same row is still a failure and is attributed to that row.
    shared = sorted(set(A) & set(B))
    n_cmp = 0
    worst: dict[str, tuple[float, tuple[str, str]]] = {}
    n_diff: dict[str, int] = {}
    all_fields: set[str] = set()
    for k in shared:
        ra, rb = A[k], B[k]
        fa = set(ra) - _NOT_A_MEASUREMENT
        fb = set(rb) - _NOT_A_MEASUREMENT
        for f in sorted(fa - fb):
            fails.append(f"field {f!r} present ONLY in {a_name} at {k[0]} {k[1]}")
        for f in sorted(fb - fa):
            fails.append(f"field {f!r} present ONLY in {b_name} at {k[0]} {k[1]}")
        fields = sorted(fa & fb)
        all_fields |= set(fields)
        for f in fields:
            n_cmp += 1
            va, vb = ra[f], rb[f]
            if _same(va, vb):
                continue
            n_diff[f] = n_diff.get(f, 0) + 1
            try:
                d = abs(float(vb) - float(va))
            except (TypeError, ValueError):
                d = math.inf
            if f not in worst or d > worst[f][0]:
                worst[f] = (d, k)
    fields = sorted(all_fields)

    print()
    print(f"      compared {n_cmp:,} scored fields over {len(shared)} rows x {len(fields)} fields")

    if n_diff:
        print()
        print(f"      ⛔ {sum(n_diff.values()):,} FIELDS DIFFER")
        print(f"      {'field':<26} {'rows':>6}  {'max |delta|':>14}   where")
        for f in sorted(n_diff, key=lambda x: -n_diff[x]):
            d, k = worst[f]
            print(f"      {f:<26} {n_diff[f]:>6}  {d:>14.6g}   {k[0]} {k[1]}")

    if fails:
        print()
        for m in fails:
            print(f"      ⛔ {m}")

    print()
    if n_diff or fails:
        print(f"   ⛔ NOT IDENTICAL — {a_name} and {b_name} differ. TRAPS: byte-identity-gate: if one of these is the")
        print("      baseline and it was recorded in an EARLIER session, re-record it (TRAPS: re-record-the-baseline) before")
        print("      believing the delta.")
        return 1
    print(f"   ✅ BYTE-IDENTICAL on all {n_cmp:,} scored fields of all {len(shared)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
