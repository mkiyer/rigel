# docs/dev — the sandbox

⭐ **Working notes live here and that is the point.** A feature write-up, a session handoff, a
half-finished derivation, a note to a collaborator. It is expected to be provisional, and nobody should
treat anything in this directory as settled.

## The two rules

1. ⛔ **Nothing outside `docs/dev/` may cite anything inside it.** Not the source, not a test, not one of
   the eight permanent docs. A citation is what turns a note into a dependency, and a dependency is what
   makes a temporary file permanent. `tests/test_docs_boundary.py` enforces this.

2. ⛔⛔ **When a finding settles, MOVE it out and delete it here, in the same edit.** Not copy — move.

   | the finding is… | its permanent home |
   |---|---|
   | a current measured number | `docs/ROADMAP.md` |
   | a mistake not to repeat | `docs/TRAPS.md` |
   | a decision that is settled | `docs/DESIGN.md` |
   | a derivation the code depends on | `docs/EQUATIONS.md` |
   | how performance is judged | `docs/SUCCESS.md` |
   | a panel, a harness, a gate | `docs/TESTING.md` |

## Why the rules exist

Two working docs once reached **1,181 lines between them — larger than DESIGN + ROADMAP + SUCCESS
combined** — and a new session was being told to read one of them as "⭐⭐⭐ THE STATE". They were deleted on
2026-08-07 and their content promoted. ⭐ **The failure was never that they existed.** It was that nothing
was ever moved out, so the provisional copy quietly became the authoritative one and the permanent docs
went stale beside it.

⚠ A stale dev doc is harmless. A dev doc that something else depends on is not.
