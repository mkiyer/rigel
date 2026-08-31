# NEXT SESSION — the state after the 2026-08-30 strand-overdispersion work

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md` §1. MOVE anything that settles into the permanent docs and DELETE this file.

## What the last session shipped

The gDNA strand overdispersion, refitted so it holds **regardless of the annotation**. The premise it
replaced — that some structural class of objects is pure gDNA — is false, and that is now a named trap
(`TRAPS: purity-is-a-property-of-the-annotation`). Three derived pieces, no new constant in any of them:
the **away-half** moment, **influence weighting**, and the two components **reconciling against each
other** instead of against a conjured `Beta(14,14)`. Derivations: `EQUATIONS.md` §6a/§6b/§6c. Rulings:
`DESIGN.md` §3.3a. State: `ROADMAP.md` §0.

Also landed: the clamp announces itself and the concentration is reported (`clamped_at_ceiling`,
`raw_overdispersion`, `effective_seeds`); the simulator gained a **shadow-transcript** mechanism and a
blank chromosome so the annotation premise stays testable (`TESTING.md` §0a); and the dead density
machinery went with the seed weight it existed to compute (`run_fill`, `RegionDensity`,
`region_gdna_density` deleted).

## The two live threads, in priority order

1. ⭐⭐⭐ **`ROADMAP.md` §1 rank 1b — the gDNA FRAGMENT-LENGTH model.** Diagnosed, not built. Same premise
   failure, in a second place: `fl.py` asserts "Every pool is PURE BY CONSTRUCTION" and the intronic pool
   measures **95 % nascent**. Read `docs/dev/DIAGNOSIS_gdna_fragment_length.md` then
   `DESIGN_gdna_fragment_length_fix.md` — the options are worked and **one is already refuted by
   measurement**, so start from what is left rather than from the obvious idea.
   ⛔ **The first move is a decision, not code**: the fix cannot be RANKED until the fl-gap arms are
   regenerated on the current sparse-nascent model, because every panel the tool is scored on gives the
   two components equal fragment lengths and therefore cannot see the defect.

2. ⭐⭐ **The message policy — the standing charter, and still at rung 0.** `MessagePolicy` remains
   byte-identical to `SilentPolicy` on all 30 test-chromosome conditions. The bar is unchanged: **win on
   unstranded, minimal harm on stranded, never pooled.** Before building any mechanism, answer rank 1a's
   question with no solver: **is the information a blind unstranded or AMBIG slot needs actually present
   in its neighbours?** The anchored twin block was designed for exactly that.
   ⚠ Re-baseline first — the policy numbers in `ROADMAP.md` predate this session's estimator change.

## Two cautions that will otherwise cost a day

* ⛔ **The benchmark resolves a LARGE effect, not a percent-level one.** A 1e−5 nudge to a nuisance
  parameter moves 1/30 rows by more than 0.5 % (worst 1.65 %, relay). Do not believe a single-row policy
  difference below ~2 %. `ROADMAP.md` §2 q8 and `TRAPS: a-constant-parked-a-value-off-a-knife-edge`.
* ⛔ **`docs/dev/` is a sandbox and nothing may cite into it.** Its own README carries the move-and-delete
  rule; a 727-line plan doc was deleted at the end of the last session precisely because its findings had
  been promoted.
