# docs/dev — the sandbox

Working notes live here: a session handoff, a half-finished derivation, a note to a collaborator.
Everything in this directory is provisional and nothing in it is settled.

## The two rules

1. **Nothing outside `docs/dev/` may cite anything inside it** — not the source, not a test, not a
   permanent doc. A citation turns a note into a dependency, and a dependency makes a temporary file
   permanent. `tests/test_docs_boundary.py` enforces this.

2. **When a finding settles, MOVE it out and delete it here, in the same edit.** Move, never copy: two
   homes diverge, and the provisional copy is the one that ends up being read.

   | the finding is… | its permanent home |
   |---|---|
   | an open problem, or a refusal with its measurement | `docs/ISSUES.md` (a named entry) |
   | a claim about the current state | `docs/ROADMAP.md`, naming the instrument that re-derives it |
   | a mistake not to repeat | `docs/TRAPS.md` (a named rule) |
   | a decision that is settled | `docs/DESIGN.md` |
   | a derivation the code depends on | `docs/EQUATIONS.md` |
   | how performance is judged | `docs/SUCCESS.md` |
   | a panel, a harness, a gate | `docs/TESTING.md` |

A stale dev doc is harmless. A dev doc that something else depends on is not.
