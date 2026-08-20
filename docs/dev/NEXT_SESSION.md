# NEXT SESSION — the dissection loop, on the worst IN-SCOPE scenario

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0/§1, the
    rulings `DESIGN.md` §0c, the derivations `EQUATIONS.md` §3.5f–h, the lessons `TRAPS.md` §D, the
    substrate `TESTING.md` §0a/§0b. ⛔ Delete this file when the plan it points at is executed.

    Written 2026-08-20 at the end of the session that built the third policy.

## ⭐⭐⭐ FIRST COMMAND OF THE SESSION

    python scripts/design/preflight.py          # ~2 min, or --fast to skip the instrument sweep
    python -m pytest tests/ -q                  # the baseline is in CLAUDE.md; ANY failure is a regression

⛔ **Do not start work on a ✘.** A missing DERIVED artifact is a command you have not run, not damage.

## ⭐⭐⭐ THE GOAL, RESTATED SO IT IS NOT LOST IN THE MECHANISM

**Pass-0 must solve enough exons ACCURATELY that a gDNA LANDSCAPE can be trained on them. That trained
landscape is then the PRIOR for the second solve, and every region and boundary is refit and re-solved
with it.** Everything else — message propagation, the reference, the policies — is in service of that
one loop, and calibration is circular precisely because the prior cannot exist until some exons are
already solved (`DESIGN.md` §0c.0d). ⛔ So "which exons are solved, and what tells us so WITHOUT truth"
is the question the tool's trustworthiness on real data rests on.

## ⭐⭐⭐ THE METHOD — the owner's iterative protocol, and it is the whole plan

Owner, 2026-08-20. Follow it in order and do not skip the dissection:

1. **RUN THE FULL PANEL.** `relay_pool_ab.py --arms off on currency --out <tsv>` over all 16 ladder
   conditions, then `benchmark_report.py`. ⛔ Report every scenario; pool only at the end.
2. **PICK THE WORST SCENARIO WITHIN THE THREE IN-SCOPE STRATA.** ⛔ Never the deferred
   unstranded × capture-ON row, which is otherwise the worst every time.
3. **DISSECT IT.** Find the REGIONS and BOUNDARIES carrying the highest error MASS
   (`worst_objects.py`), then for a SAMPLE of them trace the error to its root cause: start at
   INITIALIZATION and follow every change to that object's belief — every message it received, every
   refit — comparing each step against the certified ground truth (`calibration_oracle.py`'s
   `slot_truth.npz`). ⭐ Do this for several objects; one object is an anecdote, a sample is a mechanism.
4. **DESIGN A FIX** for the root cause, gated first and falsified by perturbation.
5. **ADD THE OFFENDING TRANSCRIPTS TO THE TEST CHROMOSOME** (`TESTING.md` §0a). ⭐⭐ This is what makes
   the loop compound: the test chromosome becomes a repository of the difficult, complicated and
   error-prone transcript combinations the tool actually fails on, and eventually stands on its own as
   the benchmark for comparing calibration algorithms.
6. **RE-RUN THE PANEL** and see how the fix behaves across all 16 — not just where it was designed.
7. **REPEAT**, on the next highest source of error.

⛔ **AND THE TEST CHROMOSOME IS NOT A SUBSTITUTE FOR THE PANEL.** This session proved they can disagree
in RANK, not just magnitude (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`). Use the toy to isolate a
MECHANISM; use the panel to decide whether it is better.

## ⭐⭐ THREE POLICIES NOW EXIST, AND COMPARING THEM IS ITSELF A METHOD

`SilentPolicy` (no messages), `RelayPolicy` (the evolved one, 19 switches), `CurrencyPolicy` (this
session's, `message_policy = "currency"`). They fail differently, and the owner's instruction is to use
that: **contrast them to learn what a production policy must do.** What is measured so far — on the
16-condition ladder, whole chain, gDNA absolute error in FRAGMENTS, contaminated rows per stratum:

| stratum (contaminated rows summed) | Silent | Relay | currency |
|---|---|---|---|
| unstranded × capture-OFF ⭐ | **357,580** | 496,226 | 583,028 |
| stranded × capture-OFF ⭐ | **291,815** | 440,209 | 635,702 |
| stranded × capture-ON ⭐ | **564,678** | 961,174 | 848,458 |
| unstranded × capture-ON ⛔ deferred | 18,559,229 | **6,021,441** | 7,890,981 |
| **the four ZERO CONTROLS** | 1,381,959 / 416,652 / 68,431 / 31,692 | 54,852 / 122,976 / 15,527 / 42,378 | ⭐ **46,137 / 15,624 / 6,441 / 3,624** — best on all four |

⛔ **On this panel message propagation is net-harmful wherever the local solve has evidence, for BOTH
policies**, and its value is concentrated where the local solve is blind (the deferred stratum and the
capture-ON zero controls, where `currency` is the best arm). ⭐ These are the CURRENT numbers, with the weighted premise fit landed
(`<ladder>/benchmark/ladder_weighted_premise.tsv`).

## ⭐ THE OPEN LEAD IS CLOSED — measured, landed, and it did NOT close the gap

**The policy's premise variance was INERT on the panel.** The unweighted method-of-moments fit is
dominated by the sparsest hops: the median ladder slot holds **13** fragments, a quarter hold under five
and a fifth hold NONE, so `mean(v_r) = 2.07` swamped `Var(log r) = 0.76` and the premise floored at
**0.0** — against **0.992** on the test chromosome (median slot 354). The term that keeps an imputation
weak never fired on the panel and ran at full strength on the toy.

⭐ **A PRECISION-WEIGHTED fit returns 0.294 / 0.324 on the two substrates** — one number, as a
library-level property should be — and is landed and gated. The panel re-run moved **11 of 16 rows, up
to 2.5×**, and made `currency` the **best arm at all four zero controls** (7.9× and 11.7× better than
the relay on the two capture-ON ones).

⛔⛔ **AND IT DID NOT CHANGE THE RANK ON THE THREE IN-SCOPE CONTAMINATED STRATA — `SilentPolicy` still
wins all three.** So the gap is NOT one mechanism, and guessing at a second one is exactly what the
dissection protocol above exists to replace. ⭐ **Start there, on the worst IN-SCOPE scenario.**

## THE MEASUREMENT DISCIPLINE THAT FOUND EVERY DEFECT THIS SESSION

* **Split the error by whether the destination HAD ITS OWN composition evidence.** "The messages help"
  and "the messages trample a measurement" are different findings and a pooled number cannot tell them
  apart — it is how an 8.09× degradation of the measured half of a chain was found inside a modest
  total (`TRAPS: an-imputation-must-cost-something-every-hop`).
* **Verify a perturbation LANDED before believing a green result.** Four separate times this session an
  ablation silently did nothing: a `sed` that matched no line; a gate that read one slot past the
  stretch it meant to isolate; a patched module binding `region_init` does not use because it holds its
  own reference to the solver; and a TRAPS rule whose uppercase name escaped the file's own counting
  regex (`TRAPS: an-ablation-that-never-ran`).
* **Name the SUBSTRATE in the same sentence as the number.**

## STILL OPEN

* ⭐ **The capture-aware reference** (`ROADMAP.md` §1 rank 3 + §0's reference row): a per-object mean
  from the measured background takes both zero controls to EXACTLY 0 and improves every capture-OFF row
  on the ladder, and regresses every stranded × capture-ON row by the documented anchor under-read.
  Ranks 2 and 3 are ONE piece of work.
* ⚠ `scripts/design/exon_solvability.py` is UNTRACKED and stamped UNDER REVIEW: its per-condition
  columns reproduce, but its headline table is not produced by the instrument itself, one claim is
  transposed, its independence claim is refuted by direct measurement, and one self-test check survives
  a perturbation that breaks the property it guards. Repair or delete it — do not quote it.
* ⛔ Whether spliced fragments are enriched depends on PROBE PLACEMENT (owner) — a junction-spanning
  probe enriches them, a probe deep in an exon does not — so no fixed rule for the flux frame can be
  right and it must be learned. `docs/dev/currency_policy_design.md` holds the statement.
* `relay_pool_ab.py`'s docstring promises a `--table pipeline` that does not exist.
