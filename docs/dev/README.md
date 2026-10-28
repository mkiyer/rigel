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

## What is here (2026-10-28)

Two session documents, each with its starting prompt: the splice artifacts' (`SPLICE_ARTIFACTS.md`, the next issue in
the tool's order: the problem, what changed since its design trail, the code today, what is measured, the open
questions) and the strand frame's (`STRAND_FRAME.md`, paused 2026-10-27).

Every other note was removed on 2026-10-28, its campaign landed or superseded, and its settled findings already in the
permanent docs: the owner's scope rulings of 2026-10-16/18 moved to `DESIGN.md` §1.0 that day. A tracked note is
recovered with `git show 84897128:docs/dev/<name>`: the splice-artifact trail (`GDNA_SPLICE_ARTIFACTS_*`,
`SPLICE_ARTIFACTS_*`), the paired-end roadmap and wave-1 plan, the observable genome's design, the region-length
handoff, and the owner's notes (`my_notes.md`, `TWO_PROBLEMS.md`). Copies of those, and the notes that were never
tracked, are in `$RIGEL_RUNS/archive/docs_dev_2026-10-28/`. Those notes are: the census length law's design, its
reviewer prompt, its five reviews and its handoff; the shared length law; the length laws' prior; the frame
correction; a separate agent's length-evidence plan; the strand frame's second-round reviewer prompt; the P2
multimapper review, with their evidence; and the multimapper intergenic alignments' document.

Earlier removals: the notes of every campaign settled by 2026-10-21 (`git show f763270c:docs/dev/<name>`; their
untracked review evidence in `$RIGEL_RUNS/archive/docs_dev_untracked_2026-10-21/`), and the strand frame's
rounds-one-to-three notebook (`$RIGEL_RUNS/archive/docs_dev_2026-10-27/`).
