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

## What is here (2026-10-25)

The paired-end repair's roadmap and its wave-1 plan, the splice-artifact design trail (`GDNA_SPLICE_ARTIFACTS_*` and
`SPLICE_ARTIFACTS_*`), the region-length dissection and its design (`OBSERVABLE_OPPORTUNITY_DESIGN.md`, 2026-10-22), the
census length law's design (revision 4.2, the frozen candidate), its reviewer prompt and the five reviews
(`CENSUS_LAW_DESIGN.md`, `CENSUS_LAW_REVIEW_PROMPT.md`, `CENSUS_LAW_INDEPENDENT_REVIEW.md`,
`CENSUS_LAW_REVISION_REVIEW.md`, `CENSUS_LAW_APPROVAL_PATH.md`, `CENSUS_LAW_INTEGRATION_REVIEW.md`,
`CENSUS_LAW_WINDOW_REVIEW.md`, 2026-10-23), the census law's handoff with everything measured after revision 4.2 —
the STAR rung, real data, the strand-specificity spectrum, the refused arm CRr and the composition-signal feasibility
(`CENSUS_LAW_HANDOFF.md`, 2026-10-25), the diagnosis of why the census law did not converge and the shared-length-law
design that follows it, with its measured record (`SHARED_LENGTH_LAW.md`, 2026-10-25), a separate agent's
alternative — one reference law and a regularized log-contrast fitted to the raw strand counts, not yet fitted or run
(`LENGTH_EVIDENCE_PLAN.md`, 2026-10-25), the frame correction for gDNA's length table — derived, checked against the
simulator's truth, what the data can measure and the defects it exposed (`FRAME_CORRECTION.md`, 2026-10-26), the
strand frame's dedicated document — the conundrum, the measured tradeoffs, what was refused, the next estimator and a
starting prompt (`STRAND_FRAME.md`, rewritten 2026-10-27; its rounds-one-to-three notebook archived), its second-round reviewer prompt (`STRAND_FRAME_REVIEW_PROMPT.md`,
2026-10-26), the length laws' prior — each origin's law from its own counts instead of the all-fragment law, and the
strand frame's gDNA table as its other half (`LENGTH_LAW_PRIOR.md`, 2026-10-26), the
owner's notes (`my_notes.md`, `TWO_PROBLEMS.md`) and the session handoff (`SESSION_HANDOFF.md`). The notes of every settled
campaign were deleted on 2026-10-21, their findings already in the permanent docs: `git show f763270c:docs/dev/<name>`
recovers one. The review evidence that was never tracked is in `$RIGEL_RUNS/archive/docs_dev_untracked_2026-10-21/`.
The multimapper intergenic alignments' document left on 2026-10-28, when the campaign landed (`DESIGN.md` §3.1g): its
last revision is in `$RIGEL_RUNS/archive/docs_dev_2026-10-28/`.
