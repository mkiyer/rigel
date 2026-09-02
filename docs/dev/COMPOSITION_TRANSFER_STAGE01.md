# THE COMPOSITION-TRANSFER THREAD — stages 0–3, 2026-09-01 (LADDER-CONFIRMED)

    ⚠ A DEV DOC — working results of the owner-approved plan in `message_notes.md`. Nothing here
    is settled; when a verdict lands it MOVES to its permanent home and this file shrinks.

**The mechanism under test (owner design, approved 2026-09-01):** one node type at a time,
starting with intron|exon BOUNDARIES; solve them by COMPOSITION TRANSFER from the adjacent
intron REGION (the shared-population pair: mature RNA cannot cross an intron|exon boundary, so
both objects' unspliced population set is {gDNA, unspliced RNA}); the exon-side message is OFF;
the open question is the transfer's HONEST PRECISION. Strictly one mechanism, zero `src/`
changes — the prototype runs through the shipped `message_policy = "message"` seam by patching
`rigel.calibration.calibrate.MessagePolicy` in-process, so every phase and refit runs it.

Substrate: the test chromosome, all 30 conditions, **100** intron|exon pairs. ⛔ Node-local
numbers below are NOT the 0.8.0 metric and NOT the panel — whole-library scoring and the ladder
confirmation are still owed before any shipping claim (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`).

## Stage 0 — the certified pair gap and the oracle-transfer ceiling (no solver)

* **The mature-cannot-cross check held exactly**: `n_mrna` at all pair boundaries, all 30
  conditions summed: **0.0**.
* **The composition gap is statistically ZERO off capture.** Noise-subtracted excess of the
  logit-composition gap (house trigamma law), per gene type: `nasc` OFF and `capnasc` OFF both
  0.000 (raw `E[d^2]` 0.36–0.45 fully explained by counting noise). **Under capture the only
  measurable excess is `capnasc` ON: sigma^2_pair = 0.138** (sigma = 0.372 nats, n = 102 pairs).
  The apparent large truth gaps seen in spot checks were counting noise — under capture an
  unprobed intron holds ~10–20 fragments while its exon-adjacent boundary holds 300–500.
* **The oracle-transfer ceiling** (hand every pair boundary its intron's realized true
  composition): tiny where the local solve is blind — `g50 ss0.50 ON` ceiling **783** against
  silent **19,832**; `g98 ss0.50 ON` ceiling 120 against silent 36,286. ⛔ And INVERTED at low
  gDNA stranded: `g05 ss0.99 ON` ceiling **1,039** against silent **148** — the neighbour is a
  misleading predictor there, matching `ISSUES: message-value-for-blind-slots`' first
  measurement. 12–21 pairs at `g05 ON` have an EMPTY intron (transfer must be silence).

## Stage 1 — the one-hop prototype, epsilon ladder + derived precision

Delivered claim, intron -> its boundary only: `lam_mode = clip(log(rho_g E_g) - log(rho_R E_r), ±L)`
read from the intron's OWN belief; `lam_prec` per arm. Arms: silent, relay, eps in
{0, 1e-10, 1e-6, 1e-3, 0.1, 1}, `drv` = 1/(1/tau_lam^intron + sigma2_pair), `flip` (mode
negated at drv precision — the perturbation falsification).

Node-local |gDNA err| at the 100 pair boundaries (fragments), the rows that decide:

| condition | silent | relay | eps0.1 | eps1 | drv | flip |
|---|---|---|---|---|---|---|
| g25 ss0.50 ON (deferred) | 10,036 | 1,785 | 6,681 | **1,706** | 10,033 | 9,828 |
| g50 ss0.50 ON (deferred) | 19,832 | 1,165 | 10,628 | **2,631** | 19,825 | 19,831 |
| g98 ss0.50 ON (deferred) | 36,286 | 78 | **1,761** | 2,519 | 4,590 | 37,750 |
| g05 ss0.99 ON (in scope) | 148 | 2,059 | 160 | 207 | **159** | 168 |
| g25 ss0.99 ON (in scope) | 323 | 1,694 | 369 | 569 | 423 | 508 |
| g50 ss0.99 ON (in scope) | 450 | 1,210 | 518 | 772 | 550 | 692 |
| g98 ss0.99 ON (in scope) | 753 | **99** | 953 | 1,379 | 1,078 | 1,777 |
| unstranded OFF g05..g98 (in scope) | 53–169 | 34–105 | 51–160 | 29–100 | 52–159 | 355–2,098 |
| g00 zero controls (six rows) | 0–100 | 5–289 | 1–76 | 10–40 | 0–89 | 17–2,083 |

What the ladder of arms says:

1. **eps <= 1e-6 moves NOTHING anywhere** — a blind boundary is not a vacuum (count likelihood +
   population prior compete), so the single tiny-epsilon probe would have measured nothing; the
   dose–response was the right instrument. Action begins at eps ~0.1–1.
2. **The mode is right and valuable where the solve is blind**: at unstranded × capture-ON the
   transfer at eps ~1 recovers most of the relay's win from ONE hop and ONE channel (vs the
   relay's whole operator stack), and `flip` destroys it — the falsification fires loudly.
3. **The harm profile is far better than the relay's on the in-scope stranded-ON rows** (relay
   2,059/1,694/1,210 vs drv 159/423/550 against silent 148/323/450) — the destination's own
   evidence survives because the declared precision is small and honest. ⚠ Exception: `g98 ss0.99
   ON`, where the relay's certified-flux anchor wins big (99) and every transfer arm is worse than
   silence — the certified-flux MEASUREMENT stream is orthogonal to this mechanism and will be
   wanted eventually (the `RowsSolve` seam exists).
4. ⛔ **The derived precision under-claims exactly where the value is.** `tau_lam^intron` at
   capture-ON introns is ~0.04 (density-vs-background reasoning is weak under capture, and the
   lambda-axis Jacobian crushes near-vertex information), so `drv` ~= silent on the blind
   capture-ON rows while eps = 1 wins. The honest tau is bracketed by measurement:
   **tau_lam (too weak at capture-ON) < honest < ~1 (too strong at stranded-ON)**. Stating the
   intron's composition posterior honestly in psi's lambda coordinate is the open derivation —
   note the family resemblance to `ISSUES: reference-prior-refuted-at-concept-level` (the
   information does not vanish in f-space; the lambda Jacobian is what vanishes). Do NOT reopen
   the refused solver-coordinate work; this is about the MESSAGE's declared precision only.
5. **Identity plumbing**: the stock rung-0 `message` policy is byte-identical to silent
   (re-confirmed in-process). Delivering a PRESENT lambda array with precision exactly 0 is NOT
   byte-inert: max |delta f_g| = 5.6e-16 (1 ULP) at exactly the delivered slots. A promoted
   implementation must deliver `None` (or claims only where the licence fires with prec > 0),
   never zero-filled channel arrays.

## Stage 2 — the honest precision DERIVED, two candidate laws priced (2026-09-01)

**The derivation.** The intron never observes its fragments' origins: it measures a TOTAL `n_I`
and knows its gDNA expectation `m_g = rho_bg * E_g` from the struct-locked intergenic anchor
(the background's reliability comes from the whole intergenic support — counts AND length). Its
RNA estimate is the SUBTRACTION `n_I − m_g`, so on the lambda axis

    Var(lambda_hat) = (n_I / n_r)^2 * ( v_bg + trigamma(n_I + 1/2) ),   n_r = (1 − f_hat) n_I

— the estimator variance is AMPLIFIED by how small the excess is relative to the total. Near the
pure-gDNA vertex `n_r -> 0`, the variance diverges, and that is not an artifact: the intron's
lambda-likelihood there is ONE-SIDED (flat top toward +lambda, a cliff below — "at least this
much gDNA"), and a symmetric `(mode, precision)` Gaussian cannot carry a cliff. `tau_lam` is the
factor's grid-wide second moment (`density_factor_precision`), so it reads the flat top's width
and collapses even when the cliff holds hundreds of nats. **The information is real; the Gaussian
is the wrong currency.** Both statements were then confirmed by measurement (below).

**The two arms.** `t_gauss` = the derived Gaussian (`tau = 1/(Var(lambda_hat) + sigma2_pair)`).
`t_rows` = the intron's actual lambda-likelihood delivered as a ROW FACTOR (`PsiMessage.lam_rows`
— psi's general evidence currency, final solve only, the certified-flux seam):
`mu(lambda) = m_g (1 + e^-lambda)`, `logL = n_I log mu − mu` (Poisson in the total the intron
saw), blurred along lambda by `N(0, v_bg + sigma2_pair)` — the transfer cost spent as a
coordinate BLUR, which weakens the cliff honestly instead of erasing it. An EMPTY intron with a
real background expectation legitimately claims "no RNA" (soft, `~m_g e^-lambda`), which is the
"zero counts over long support" intuition as a likelihood. `t_rows_flip` mirrors the factor —
the falsification.

**The verdict on the test chromosome (all 30 conditions, node-local at the 100 pair boundaries
+ whole-library, tables in the stage-2 harness output):**

* `t_gauss` behaves exactly as derived: honest, near-zero harm, and UNABLE to deliver the blind
  capture-ON value (13,253 at `g50 ss0.50 ON` node-local vs silent 19,832) — the cliff is the value.
* ⭐ **`t_rows` meets BOTH bars at this node type.** Node-local: beats silent on every unstranded
  row (e.g. `g05 ON` 556 vs 2,018; `g25 ON` 1,419 vs 10,036; `g98 OFF` 49 vs 160) and does **no
  stranded harm** — ties or WINS every `ss0.99` row (`g05 ON` 148 vs 148; `g98 ON` 714 vs 753;
  worst degradation anywhere +8 fragments at `g05 ss0.99 OFF`). The relay's stranded-ON harm at
  these nodes (2,059/1,694/1,210) simply does not occur. Whole-library: small consistent wins on
  11/12 in-scope contaminated conditions, one +8. `t_rows_flip` degrades everything — the gate fires.
* The remaining gap to the relay on the deferred rows is the EXON node type (the relay solves
  exons too; this mechanism touches only 100 boundaries): `g50 ss0.50 ON` whole-library rows
  68,532 vs relay 6,695 — that gap is the next rung's target, not this one's failure.
* ⚠ The one cost: small zero-control claims (`g00` node-local up to 147 vs silent's 100; up to
  +83 whole-library at `g00 ss0.70 ON`) — the intergenic background is not exactly empty at g00,
  so `m_g > 0` lets near-empty introns prefer gDNA. Real, small, to be priced by `zero_controls.py`
  standards before shipping.
* ⚠ Poisson (no od_g) in the factor and sigma2_pair carried from the test chromosome are the two
  prototype simplifications to revisit at promotion.

## Stage 3 — THE LADDER CONFIRMS (16 conditions, 19,610 pairs, silent/relay/t_rows/flip)

⭐ Credibility check first: the harness's whole-library scoring reproduces the recorded baseline
TO THE FRAGMENT (silent 298,597 / relay 456,838 at `g98 ss.99 ON`; 126,467 / 208,306 at OFF), so
every number below is on the official basis.

* **The WIN bar (unstranded), node-local**: t_rows beats silent on all 6 contaminated rows —
  `g05 ON` −53 % (107,541 → 50,053, ≈ relay), `g50 ON` −69 % (1,240,014 → 387,440), `g98 ON`
  −72 % (2,425,978 → 681,294), `g98 OFF` −54 % (45,910 → 21,322, relay 36,684). The relay stays
  ahead on the two deferred capture-ON rows (it also solves EXONS — the next rung's target).
* **The HARM bar (stranded): not merely minimal harm — t_rows WINS every `ss0.99` row**, node-local
  AND whole-library, and beats the relay on all of them: at the worst in-scope condition
  (`g98 ss.99 ON`) whole-library is **silent 298,597 / relay 456,838 / t_rows 281,835** — where
  the relay ADDS +158k (+53 %), the transfer REMOVES 16,762 (−5.6 %). `g98 ss.99 OFF`:
  109,963 vs silent 126,467 (relay 208,306). `g50 ss.99 ON`: 242,155 vs 260,629 (relay 409,168).
* **Whole-library, in scope: t_rows ≤ silent on ALL 12 contaminated conditions.**
* **Zero controls: t_rows ≡ silent EXACTLY on all four `g00` rows** (empty intergenic ⇒ no
  background expectation ⇒ every row inert) — the test chromosome's small g00 cost was a
  small-substrate artifact (its 20 kb spacers collect edge fragments; the ladder's intergenic is
  clean). ⚠ The relay's large `g00` whole-library wins (e.g. 1,254,145 → 58,840 at
  `g00 ss.50 OFF`) come from machinery this rung does not carry — that value is a LATER rung's
  target, and t_rows never goes below the silent floor there.
* `t_rows_flip` is catastrophic everywhere — the falsification fires on the ladder too.
* ⭐ The two substrates AGREE in sign on both halves (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`
  satisfied for this mechanism).

**Verdict: the first rung of the new policy — intron -> intron|exon boundary composition
transfer, delivered as the intron's blurred lambda-likelihood row factor — is measured, falsified,
and meets both bars on the shipping substrate.** Before src/: re-derive sigma2_pair on ladder
slot_truth (borrowed 0.138 here), decide the od_g/NB question in the factor, and take the owner's
promotion ruling (foundation-spec implementation + gates).

## Item ① — sigma2_pair on the LADDER's certified truth (stage0b, no solver)

Pooled per capture state (the deployable form): **capture-OFF 0.000** (noise explains everything,
n = 32,534 mixed pairs) and **capture-ON 0.200** (n = 3,155) — the borrowed 0.138 was close.
Two footnotes recorded so they are not rediscovered: `g98 OFF` alone shows ~0.40 excess (buried
in the pooled OFF value by the g05/g50 mass), and `g98 ON` has only ~22 mixed pairs (its 3.9–5.7
per-condition values are sample noise, not structure).

## Item ② — the od question, answered by construction + measurement (stage 4)

The hand-rolled Poisson-total factor duplicates a likelihood the tree already owns: the intron
factory's `density_lambda_factor` — NegBinom with `alpha_eff` fusing the per-region
overdispersion AND the background posterior's own width (`size >= 1/2` so an empty pool is wide,
never confident). So the production form is **the intron's own shipped factory row, blurred by
sigma2_pair, delivered one hop** (`t_fact`) — no second likelihood to keep in step, od handled
where it already lives, and the background width not double-counted (the blur is sigma2_pair
alone). One behavioural difference vs the Poisson-total form: a zero-count intron's factory row
is FLAT (no claim) where the Poisson-total form softly preferred gDNA — the factory is the more
conservative statement. Stage 4 races both on the ladder beside silent/relay + the flip
falsification.

## Stage 4 — the FINAL ladder A/B (ladder sigma2 + both factor laws; 16 conditions, 5 arms)

Both variants meet BOTH bars on the ladder, and they split cleanly:

* **`t_fact` (the shipped factory NB row, blurred, one hop) owns the ZERO CONTROLS**: node-local
  `g00` 4,155 / 4,559 / 685 / 1,635 against silent's 57,060 / 7,517 / 7,730 / 2,855 — it beats
  even the relay on 3 of 4 `g00` rows at these nodes, because a zero-background factory row
  truthfully claims "gDNA ≈ 0". Whole-library it improves every `g00` row. Its single in-scope
  degradation anywhere is **+169 fragments** (`g05 ss.99 ON`, +0.2 %); it beats silent on the
  other 11 contaminated in-scope rows.
* **`t_rows` (Poisson-total) owns the contaminated capture-ON rows**: node-local `g50 ss.50 ON`
  389,446 vs t_fact 880,317 (silent 1,240,014); `g98 ss.50 ON` 686,025 vs 1,452,996; in-scope
  stranded-ON ~20 % better node-local than t_fact. Off capture the two tie (t_rows slightly
  ahead everywhere). t_rows ≡ silent exactly at `g00` (its factor is inert with no background).
* `t_fact_flip` fires massively everywhere the claims are strong — the falsification holds.
* Whole-library at the worst in-scope condition (`g98 ss.99 ON`): silent 298,597 / relay 456,838 /
  t_rows 281,769 / t_fact 294,434.

**The attribution of the split**: one likelihood, two strengths. The factory's `alpha_eff` folds
the fitted od and the background posterior's width into the claim, widening it — honest on real
data, but the panel's gDNA is POISSON BY CONSTRUCTION (`synthetic_suite_is_poisson`), so the
sharper Poisson-total claim wins there partly by simulator fiat. The factory form is also the
structurally correct law at the zero corner (it scores `f*C` against the background, so `C > 0`
with an empty background forces `f -> 0`; the Poisson-total form degenerates there).

⭐ **RECOMMENDATION (awaiting the owner's ruling): promote the FACTORY form.** One shipped
likelihood, correct at every corner including `g00`, honest under real-data overdispersion, both
bars met with a worst in-scope cost of +169 fragments. The capture-ON upside of the sharper claim
(~2x at deferred rows, ~20 % node-local at in-scope stranded-ON) is RECORDED as the priced
subject of a follow-up derivation: which share of `alpha_eff` belongs in a one-hop COMPOSITION
claim (od is a rate-claim price; a composition is a ratio and part of the dispersion may cancel).
⛔ Not both behind a switch — one law.

## Stage 5 — the blur A/B, and THE CONSTANT IS DELETED

The three-arm race (no blur / uniform 0.200 / split by capture) on the full ladder: **the blur
does no measurable work** — every condition moves under half a percent between arms, the arms
tie exactly off capture, and the UNBLURRED row is slightly the best at the largest rows
(`g98 ss.50 ON` node-local 1,443,103 vs 1,452,996). `alpha_eff`'s own width dominates the
measured pair dispersion. So the hop cost is PRICED (0.000 off capture / 0.200 on, certified)
and paid at zero beyond what the row already carries — and **no dispersion constant ships**:
`transfer_pair_dispersion` was removed from the config, the blur from the policy, and its gate
from the suite. If the dispersion is ever re-priced as material, it re-enters as a
runtime-fitted widening (the `splice_in_premise_logvar` pattern), never as a constant.

## THE PROMOTION — landed in src/ (2026-09-01, uncommitted; the owner drives commits)

* `src/rigel/calibration/messages/transfer.py` — `TransferPolicy`: the intron's own factory
  rows (grid-keyed via calibrate's memoized `_intron_prior_at`, the flux-rows pattern),
  delivered VERBATIM as `PsiMessage.lam_rows` at the structurally derived intron|exon pair
  boundaries; `scan -> None` (one hop structural); silence — never zero-filled channels — when
  there is nothing to say. Declared in `_layers.py` (layer 6).
* `config.message_policy` gains `"transfer"`; the unknown-name refusal keeps its gate. No new
  config constant. `policy_benchmark.py` gains the `transfer` arm (named, not default).
* `tests/calibration/test_transfer_policy.py` — 4 gates, each watched FAILING first and each
  watched FIRING under a deliberate break (wrong source side, non-None scan, wrong policy
  installed; one gate was STRENGTHENED when an inverted blur kernel slipped a widening-only
  predicate — since deleted with the blur; the drop-exon-flank break is unfalsifiable on the
  gate toy and the docstring says so).
* Faithfulness: the installed policy reproduces the harness arms TO THE FRAGMENT through
  `policy_benchmark.py` (`g98 ss.99 ON`: silent 298,597 / transfer 294,676 = the no-blur arm).
* ⭐ **The suite baseline moves: 0 failed / 3,744 passed / 0 skipped / 8 xfail, 3,751
  collected** — accounted from 3,733/3,741: +3 the new calibration module, +2 the new
  tests/calibration file, +4 its own gates, +1 the docs/dev note. (CLAUDE.md's baseline line
  is updated by the commit that lands this, per its own rule.)

## Next, in order

1. The owner's commit of the landing (tree is clean-green, uncommitted).
2. **The exon-region rung** — the second node type, same discipline (derive -> prototype ->
   A/B -> ladder -> promote). Its targets, both measured here: the relay's remaining lead on
   the deferred unstranded-capture-ON rows (the relay solves exons; the transfer touches only
   boundaries) and the relay's large `g00` whole-library wins.
3. The prototype harnesses (`stage0`–`stage5`) live in the session scratchpad; the shipped
   gates and instruments carry everything load-bearing, so they die with the session by design.
