#!/usr/bin/env python
"""⭐⭐⭐ **IS THIS ARM BYTE-IDENTICAL TO THAT ONE?** — the byte-identity gate, and it refuses to lie in either
direction.

`arm_score.py` AGGREGATES, so it answers "how much did this move?". A restructure gated on byte-identity
needs the other question — "did ANY field move at all?" — and a sum cannot answer it: two fields that move
by ±x cancel in a total and read as identical. This compares **every scored field of every row**.

⛔⛔ **TRAPS: byte-identity-gate SAYS THIS INSTRUMENT HAS LIED IN BOTH DIRECTIONS, so both failure modes are gated
here rather than trusted:**

* *"An arm with ZERO rows scored 32/32 IDENTICAL because the comparison looped over the new arm's rows."*
  ⇒ the row-key sets must be **EQUAL**, not overlapping, and the comparison runs over the UNION. A row
  present in one arm and absent from the other is a FAILURE, never a skip.
* *"A stored baseline went stale so unmodified HEAD no longer reproduced it."*
  ⇒ every field is named in the output with its own max |Δ|, so a mismatch says WHICH field and WHERE
  rather than only that the totals differ. And the two files' `mtime`s are printed, because a baseline
  recorded in an earlier session is exactly the thing TRAPS: byte-identity-gate warns about (`TRAPS: re-record-the-baseline` — re-record it).

⛔ **A `nan` equals a `nan` here, and that is deliberate.** `library_f_gdna_*` is `nan` at a condition
where the field was not produced, and `nan != nan` would report every such field as a difference on both
arms. Identity is `(both nan) or (bit-equal)` — compared as BITS (`==` on the float), never within a
tolerance, because a tolerance is what turns "byte-identical" into "close enough" and this gate exists to
refuse that.

Usage
-----
    python scripts/design/arm_identity.py base backbone_head       # names under $RIGEL_ARMS
    python scripts/design/arm_identity.py a.jsonl b.jsonl          # or explicit paths

Exit status is 0 only if every field of every row is identical, so it composes into a shell gate.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

#: where `ladder_arm_ab.py --out` wrote the arms. Override with $RIGEL_ARMS.
D = Path(os.environ.get("RIGEL_ARMS", Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_arms"))

#: ⛔ NOT compared: fields that differ BY CONSTRUCTION and would report a difference on every row.
#: Everything else in the row is a measurement and is compared.
#:
#: * ``arm`` — the arm's own name.
#: * ``seconds`` — wall clock. It is an observation of the MACHINE, not of the arm, and it never
#:   repeats. ⚠ Added when ``quant_accuracy.py`` became the first producer to record one; keeping the
#:   set this small is deliberate, because every name added here is a way for a real difference to
#:   stop being reported. A field belongs here only if it CANNOT be equal between two identical arms.
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

    # ── TRAPS: byte-identity-gate, failure mode 1: the row sets must be EQUAL. A missing row is not a skipped row. ────────────
    only_a, only_b = sorted(set(A) - set(B)), sorted(set(B) - set(A))
    if only_a or only_b:
        for k in only_a:
            fails.append(f"row present ONLY in {a_name}: {k[0]} {k[1]}")
        for k in only_b:
            fails.append(f"row present ONLY in {b_name}: {k[0]} {k[1]}")
    keys = sorted(set(A) | set(B))
    if not keys:
        fails.append("BOTH arms are empty — this gate would otherwise pass vacuously (TRAPS: byte-identity-gate)")

    # ── the field sets must match too: a field added or dropped is a change to what was measured ──────
    #
    # ⛔⛔ **PER ROW KEY, NEVER OVER THE WHOLE FILE.** This used to union the field names across every
    # row of an arm and then demand every row carry every name. That holds only while all rows share one
    # schema — true for `ladder_arm_ab`, whose two axes are `region` and `boundary` — and it reports a
    # spurious failure the moment an arm emits rows of DIFFERENT shapes on different axes:
    # `quant_accuracy.py` writes a `transcript` row and a `library` row per condition, and the global
    # union made each one look like it was missing the other's fields. **1,296 false "field missing"
    # failures on two arms that were in fact byte-identical.**
    # ⭐ The per-key comparison is also strictly STRONGER, not a relaxation: a field present in A's row
    # and absent from B's SAME row is still a failure, and it is now attributed to that row instead of
    # being masked by any other row that happened to carry the name.
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
