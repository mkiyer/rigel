# NEXT SESSION — the cleanup (owner, 2026-09-22)

The owner stopped the per-fragment-prior derivation line ("horribly overcomplicated; if we are thinking too hard we
are off track") and set the frame: `docs/dev/TWO_PROBLEMS.md` restated in two lines — a prior per transcript, and
the capture-contracted effective length per transcript — and BEFORE either is worked, a wholesale cleanup of the
code, prototypes, gates and infrastructure to a production-ready tree WITH THE SHIPPED INFRASTRUCTURE and no new
mechanism, as if 0.8.0 shipped what runs today.

The ledger is `docs/dev/CLEANUP.md`: what is done with its proof, the census, the owner's decisions each blocking
a deletion, and the stages ahead. Start there. The owner drives commits; the branch `calibrated-likelihood`
carries the previous campaign's inert refactors and findings plus today's retirements uncommitted — propose the
commit points, one kind of change each.

## What the next session does, in order

1. `python scripts/design/preflight.py --full`, then the suite; re-derive the count from `CLAUDE.md`'s table
   (today's retirements and the knob deletion move it; the ledger records the derivation).
2. Read the coverage table in `CLEANUP.md` and work stage 2: every function the suite never executes is read and
   either gated or deleted, one commit each, `rename_identity.py --check` after each.
3. Put the decisions in `CLEANUP.md` to the owner in plain words, one line each, and act on the answers.
4. Then stages 3–6 in order. The instrument retirements move `CLAUDE.md`'s table and the collected count by the
   table's rows.

## What not to do

No new mechanism, no new instrument, no derivation. Nothing that moves a number: the identity check is the proof
of every step. Do not touch the pseudocount, the ruler or the gDNA component's length. Keep `docs/dev/` to this
file, the ledger, `TWO_PROBLEMS.md` and the README.
