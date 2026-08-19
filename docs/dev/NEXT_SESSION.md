# NEXT SESSION — Stage 2, "THE MAP": measure the currency per hop type

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0/§1, the
    rulings `DESIGN.md`, the lessons `TRAPS.md`. ⛔ Delete this file when the plan it points at is executed.

    Written 2026-08-19, replacing `NEXT_SESSION_zero_control_pin.md` — which is SUPERSEDED, not lost: it
    queued a per-bug dissection of `RelayPolicy`, and the owner has since ruled that the relay is REBUILT
    rather than repaired. Its durable content was MOVED before deletion (see "what moved where", below).

## THE TREE

Clean, four commits ahead of where the day started. Suite **0 failed / 3,532 passed / 0 skipped /
9 xfailed, 3,541 collected**. `message_propagation = True`.

⭐ **The panel is fully rebuilt and certified**: 32/32 scan caches, 32/32 oracle caches, all 32
`slot_truth.npz` stamped **COMPOSITION + FIELD**, zero class-bias flags.
⚠ **Read that honestly**: the uniformity gate is VACUOUS on the 16 capture-ON conditions by design, so the
real verification is **16/16 capture-OFF, up from 8/16** — all eight failures were capture-OFF, so the win
is the same size, but "32/32" is a stamp count and not 32 verifications.
⚠ `panel.py status` still prints `oracle cache 11/32`; that is the DOCUMENTED quirk (it counts a fourth
part `_main`, written only for conditions `pass0_vs_oracle` SCORES, and redundant with `scan_cache`).

## WHAT THE NEXT SESSION DOES

⭐⭐⭐ **`PLAN_message_propagation_rebuild.md` §3a, in full. It is written to be executed.** In one line:
build `scripts/design/hop_currency.py` and run it on all 32 conditions, so the CURRENCY of each of the six
hop types is MEASURED before any policy code is written.

⛔ **What NOT to do, and this is the whole reason the previous handoff was replaced**: do not open a
per-bug dissection of `RelayPolicy`. Its zero-control gap is real, named and open — and its policy is
scheduled for deletion, so fixing it is work with no destination. The mechanism survives as a CONSTRAINT
the new policy must satisfy (`TRAPS: zero-the-precision-with-the-value`).

## WHY THIS IS THE NEXT STEP AND NOT SOMETHING ELSE

The owner's ruling (2026-08-18): *"the whole method from the top down is not simple, and it is not
elegant … it can be simple and elegant"*. `RelayPolicy` carries **19 switches**; the panel says the problem
has **six hop types**. `DESIGN.md` §0c.0 is the ruling that replaces the switch pile — a hop carries a
**LEVEL** or a **COMPOSITION**, and the two invariances are complementary, so the licence question IS the
currency question asked once per hop type.

⭐ Stage 2 is the step that makes that testable rather than argued. It costs no `src/` change and it
cannot break anything.

## THE ONE THING THAT WOULD OTHERWISE BE LOST

⛔⛔ **Every measurement in the plan's §2 is STALE and the plan now says so at the top of the section.**
It was taken against the PRE-FIX caches; the chimera repair moved boundary crossings by **2.4 %**, on
exactly the axis the currency oracle reads. Use §2 for the SHAPE and the METHOD, never as a result.
⚠ Two rows are additionally wrong on their own terms — the SPLICE IN row used the wrong source population,
and the `nrna_none` rows make a COMPOSITION "win" that is `nrna = 0` restated. §3a says how to fix both.

## WHAT MOVED WHERE (so nothing is hunted for)

| what | now lives in |
|---|---|
| the two-currency principle, the direction rule, the three settled points (forward-backward, three arms, no analogy) | `DESIGN.md` §0c.0 |
| compatibility-before-chimera | `DESIGN.md` §3.1b-ii |
| a refused claim must lose its precision with its value | `TRAPS: zero-the-precision-with-the-value` |
| a transcript predicate must not silently drop a molecule; make the fragment ledger CLOSE | `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule` |
| a class ratio needs a population that can supply its numerator | `TRAPS: a-ratio-needs-a-population-that-can-supply-its-numerator` |
| the state, the numbers, the ranked list | `ROADMAP.md` §0 and §1 ranks 1–2 |
| the two refused zero-control mechanisms (`struct_lock = g1_locked ∧ REGION`, the zero-mass rescale guard) | `ROADMAP.md` §4.1, the GRAVEYARD — with their named losses, so they are not rebuilt |

## THE QUEUE BEHIND STAGE 2

* Stage 3 — the third policy, rung by rung on toys, ending at the owner's TA–TF locus; Stage 4 — the panel
  in RAW COUNTS; Stage 5 — retire `RelayPolicy`.
* `arm` → `component`, the one token the vocabulary pass deliberately deferred: **THREE** senses, one of
  them `__ARM_NEON` in the C++ scanner, 1,358 sites. It needs its own `--sense` pass.
* The sj artifact BLACKLIST drops a whole fragment (`SPLICE_ARTIFACT` deposits nothing) — the chimera
  bug's twin, inert on the simulated panel (`sj_blacklist_size = 0`) and live only on real data. An owner
  decision, and the owner's chimera ruling arguably already answers it.
* `ROADMAP.md` §1 ranks 3+ are unchanged and independent of the rebuild.

## THE DISCIPLINE THAT PAID, AND THAT COST WHEN SKIPPED

* ⭐ **Run `rename_census.py --sense <token>` BEFORE renaming anything.** Skipping it corrupted the second
  sense of one token in ~50 places; the whole stage was reverted and redone scoped.
* ⭐ **A falsification test first, verified failing — then break the fix and watch each gate fire.** That
  second half found a hole: breaking the C++ to measure a GAP instead of a SPAN fired no native gate until
  one was written for it.
* ⭐ **Make the ledger close.** `n_fragments == deposited + Σ accounted drops` is what named the chimera bug
  after eleven other candidates had been eliminated.
* ⭐ **Score per gDNA-bearing reference, in FRAGMENTS.** Pooling references manufactured a 2 % deficit that
  was reported as a second bug before it was caught.
* ⭐ **Read `memory_pressure`, not RSS.** A condition's 45 GB RSS is reclaimable page cache; believing it
  would have forfeited a 6× rebuild speedup.
