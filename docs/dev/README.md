# docs/dev — the sandbox

Working notes live here while a piece of work is in progress. Nothing here is settled.

## The two rules

1. **Nothing outside `docs/dev/` may cite anything inside it**: not the source, a test or a permanent doc. A citation
   turns a note into a dependency. `tests/test_docs_boundary.py` enforces this.
2. **When a finding settles, MOVE it out and delete it here**, in the same edit. Move, never copy:

   | the finding is… | its home |
   |---|---|
   | an open problem | `docs/ISSUES.md` |
   | the current state and the next steps | `docs/ROADMAP.md` |
   | a mistake not to repeat | `docs/TRAPS.md` |
   | a settled decision | `docs/DESIGN.md` |
   | a derivation the code depends on | `docs/EQUATIONS.md` |
   | how performance is judged; panels and gates | `docs/SUCCESS.md`, `docs/TESTING.md` |

## What is here

`CAPTURE_IDEAS.md`: the working memory for the capture reader, kept at the 2026-11-03 reset. It holds the ideas from
the 2026-10-29 to 2026-11-03 campaigns judged most promising, the measured facts any new reader must respect, the owner's
standing rulings, and a short list of dead ends.

Everything else was retired at the reset. The design documents, reviews and session notes are in
`$RIGEL_RUNS/archive/cleanup_2026-11-03/docs/dev/`; earlier archives sit beside it.
