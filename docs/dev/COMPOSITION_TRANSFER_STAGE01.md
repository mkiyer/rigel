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

# RUNG 2 — EXON REGIONS: the derivation + design (2026-09-01, owner notes in message_notes.md)

## The route reformulation of the owner's licence, and the one sharp law

The owner's archetypes: equal populations (no TSS/TES, no strand change) -> composition
transfer; unequal -> not. The derivation sharpens this: at a licensed face the exon's population
is NOT equal to the boundary's — the exon additionally holds RNA that arrived by the SPLICED
route — but that difference is CERTIFIED-MEASURED (the face's route-summed flux). So the law is:

    ⭐ A population difference that is MEASURED at the face does not refuse the hop —
      it becomes part of the message. Only an UNMEASURED difference refuses
      (a terminus: molecules enter uncounted; a strand flip: membership changes).

The exon message is therefore composition transfer (the contiguous routes: gDNA + unspliced RNA,
continuous across a licensed face) PLUS the measured spliced route — never species, always
routes (AXIOM 0: the split is "how did RNA arrive at this face", not "what kind of RNA is it").

## ⛔ OWNER RULING (2026-09-01): the anchor's cross-locale assumption is REFUSED as rung 2's basis

The first draft of this section proposed migrating `rna_anchor`'s flux rows to exons. The owner
refused the load-bearing assumption — *"the RNA-frame boundary->exon ratio is capture-invariant"*
is a LEVEL claim across two locales, and it holds only under benign probe geometry. Two real
probe placements break it in OPPOSITE directions: a sparse mid-exon probe with short fragments
depletes the faces relative to the interior (flux under-reads the exon), and a junction-spanning
probe enriches spliced fragments over everything else (flux over-reads). An opposite-sign,
placement-dependent failure is a BIAS, not a variance — no global width can price it honestly.

⭐ **Why the measured 1.00–1.03 route identity was exact — and why that VINDICATES the concern**:
`test_probes.bed`'s own header says probes tile PER EXON, 8×125 bp SEAMLESS, abutting both exon
ends, **deliberately never spanning an sj**, and the sampler weights a fragment by its best
single-probe overlap — so on this substrate spliced, crossing and contained fragments at a face
share near-identical capture weight BY THE PANEL'S OWN DESIGN. The identity is a property of the
benign tiling, not of biology. (`TRAPS`-adjacent: real data is a test input, never a design
input — and so is a benign simulation.)

The owner also REFUSED two-flank disagreement as the width RULE: it is computable only where
both flanks are comparable, which is systematically the SIMPLE subpopulation (single isoform, no
terminus, one strand), so a width fitted there under-covers exactly the exons where the model is
most wrong. It survives only as a diagnostic on that subset. The replacement is a DERIVED width
(counting, opportunity/length conversions, the propagated split evidence) checked
population-unbiased by the calibration ratio (declared vs realised, `solvability_audit`), plus —
on stranded libraries — supervised validation against exons' own strand solves.

## Rung-2 stage 0 — the licence census and the two route identities (test chromosome, certified)

* **Licence census**: middle exons 25/25 BOTH faces licensed; outer exons 50/50 exactly ONE
  (the gene-end face carries TSS/TES). Archetype A covers every twin-block exon; outer exons
  are the one-message case.
* **Spliced route** (face route-rate vs exon true mature density): median ratio **1.00-1.03 at
  every gene type with flux, capture ON and OFF alike** (cap ON 1.010, capnasc ON 1.016) — the
  identity is essentially exact and capture-proof. `clean` ON reads 1.26 on 92 thin pairs
  (unprobed under capture — counting noise); `silent` is all-zero (the route cannot produce a
  false positive).
* **Contiguous route** (boundary vs exon true nascent density): **1.02 off capture**;
  **0.775 at `capnasc` ON** — the probe-shoulder gradient: the boundary sits at the enrichment
  shoulder while the exon interior is fully enriched. Real, directional, ~25 %; the honest cost
  of the contiguous transfer under capture (the anchor's recorded "+12 % mature offset stays in
  the width" is the same family).

## THE RUNG-2 DESIGN (owner direction): the FACE-COMPOSED COMPOSITION TRANSFER

**The whole face composition transfers** — spliced RNA + unspliced RNA + gDNA at the face, as
FRACTIONS, to the exon. Every ratio is formed WITHIN one locale, so the surviving assumption is
only that enrichment is COMPONENT-BLIND at each locale separately (each ratio cancels its own
locale's enrichment) — strictly weaker than the anchor's cross-locale equality. The one named
residual premise: spliced vs unspliced fragments at the SAME face must share capture affinity
(violated by deliberate junction probes); it is a library-measurable contrast (the sj-depletion
family) and can later feed a derived width — never a flank-pair fit.

**The construction** (all ingredients exist). At a licensed face `b` into exon `e`, with
densities g (gDNA), r (unspliced RNA), s (spliced route, route-summed):
`n_u = g·A_g + r·A_r` (the crossing count and its face opportunities), the crossing split
`lam_u = log(g·A_g) − log(r·A_r)` constrained by rung 1's intron row, and `s` measured with
counting width. The exon's composition claim is the MONOTONE map

    lam_e(lam_u) = log( n_u·sigma(lam_u)·E_g^e / A_g )
                 − log( ( n_u·(1−sigma(lam_u))/A_r + s ) · E_r^e )

and the exon row = the intron row PUSHED FORWARD through this map (change of variables on the
grid), composed with the flux count's marginal (quantile nodes, the honest-marginal pattern).
⭐ **The structural ceiling**: as `lam_u -> +inf` (crossing pure gDNA), `lam_e` SATURATES at
`log(n_u·E_g^e/A_g) − log(s·E_r^e)` — certified flux structurally CAPS the exon's claimable
`f_g`. No prior mechanism had this; it is the certified-RNA bound done in composition space.
gDNA needs no separate lane — it rides inside the transferred composition, and no level ever
crosses a locale.

**Derived width, per the owner's list** — count precision (`n_u`, `n_s` trigamma), length/
opportunity conversions (`A_g`, `A_r`, `A_sj`, `E^e` — exact geometry on the panel; the fl
model's priced uncertainty on real data), the split evidence's own width (the intron row,
propagated through the map's Jacobian — automatic in the pushforward), and the named premise
residual (recorded, initially unpriced, later a measured library contrast). No fitted pairs.

**Known-by-construction behaviour**: `g00` — intron row ≈ all-RNA, transfers, no false gDNA;
`silent` genes — no flux + pure-gDNA row -> pure-gDNA claim (true); empty face -> silence; two
licensed faces -> two rows summed by psi (left face uses the left intron, right the right —
independent sources; short-exon shared-fragment correlation across faces is a recorded residual);
outer exons -> one row (the terminus face refuses, archetype B deferred).

## THE ADVERSARIAL PROBE PANELS — built, certified, and the falsification FIRED (2026-09-01)

Two full 30-condition panels beside the benign one, drafted from the GTF (hand-edited-class
files pending owner review: `test_probes_sparse.bed` — one 125 bp probe centred per exon;
`test_probes_junction.bed` — one BED12 two-block probe per sj, 62+63 bp, contiguous only in
cDNA; configs `test_reference_probes_{sparse,junction}.yaml`, ONLY outdir+probes differ).
Simulated, cached, `slot_truth` certified 30/30 each (the recorded g00-prewarm + `_main`-copy +
`calibration_oracle.py` recipe). ⭐ The sampler expresses the real chemistry: block overlaps SUM
within a probe group, so a junction probe enriches spliced fragments over unspliced/gDNA at the
same face; `gdna_split_penalty` applies on top.

**The route identities on hostile probes (probed gene types, capture-ON; off-capture rows are
byte-identical across all three panels — the internal control):**

| panel | spliced-route ratio | contiguous ratio |
|---|---|---|
| benign (seamless tiling) | 1.00–1.03 | 1.02 (0.78 capnasc-ON) |
| SPARSE (mid-exon probe) | **0.003** (~300x under-read) | (faces too thin to read) |
| JUNCTION (sj-spanning) | **24–26x over-read** | **11x** at capnasc-ON |

⭐ The benign 1.00–1.03 was PROBE-PLACEMENT LUCK, as the owner suspected — the anchor's
cross-locale level assumption is off by 2.5 orders of magnitude in one direction on sparse
probes and 1.4 in the other on junction probes, on the same chromosome, same truth, same tool.
Unprobed gene types stay ~1.0 on both panels (the break is exactly where probes are).

## RUNG-2 v0 PROTOTYPE — the face-composed transfer, first sweep (benign panel, whole-library)

`t_face` = rung 1 + the exon rows (the intron row transported through the face map; soft flux
ceiling; silence on an empty face; likelihoods transport with NO Jacobian). Wins where it must:
`g25 ss.50 ON` 2,556 (silent 43,428, **relay 8,470**); `g50 ss.50 ON` 2,936 (silent 86,453,
relay 6,695) — face-local composition BEATS the anchor's cross-locale claim on the anchor's own
friendly ground. Relay keeps `g98`/`g00` ON (900 / 6,017 vs 8,333 / 26,435). Flip falsification
fires (10,072 -> 54,663 at `g00 ss.50 OFF`).
⚠ **v0's measured shadows, cause understood**: `g05 ss.50 ON` +8,909 vs silent (the recorded
misleading-neighbour regime), `g25 ss.50 OFF` +5,399 (mechanism NOT yet dissected), mild
stranded-ON harm (+150…+330). v0 transports the row at FULL strength through a map built from
noisy point ingredients (`n_u`, `s`) — **v1's owed derivation is the ingredient width**: the
counting uncertainty of `n_u` and `s` propagated through the map's Jacobian into a per-face
widening. Derived, per-face, no constants, no pairs.

## THE PRICING VERDICT ACROSS THREE PROBE WORLDS (whole-library, capture-ON rows)

**The anchor is dead, by measurement, on its own chromosome.** On SPARSE probes the relay
(carrying the anchor) is catastrophic: `g25 ss.99 ON` **140,424** vs silent 1,415 (**99x harm on
a stranded row**), `g05 ss.70 ON` 182,421 vs 1,233, `g50 ss.99 ON` 98,585 vs 6,914 — the 300x
flux under-read becomes a confident wrong level claim that destroys even rows the local solve
had right. On JUNCTION probes the relay is harmed on EVERY capture-ON row (21,663 vs silent
1,387 at the `g00 ss.50` zero control). It is good only on the benign panel it was measured on.

**`t_face` is best-or-near-best on nearly every capture-ON row of ALL THREE panels** — the
owner's robustness expectation confirmed, and stronger: even on the junction panel (the declared
unmodelable case) it WINS almost every capture-ON row (`g25 ss.50 ON` 2,200 vs silent 22,866;
`g50 ss.50 ON` 3,113 vs 45,011; `g98 ss.50 ON` 2,318 vs 88,385; every stranded-ON row a small
win). **Why the graceful degradation is structural**: a probe-affinity asymmetry enters a
COMPOSITION logarithmically (a 2x spliced-over-unspliced face asymmetry is ~0.7 nats of tilt in
one component of a ratio), while it enters a LEVEL claim in full (the 25x over-read is 3.2 nats
applied as a location); and interior depletion cancels inside the exon's own composition
entirely. Bounded damage vs unbounded damage — the owner's "compromise" is the right currency.

**v0's standing shadows** (each appears on every panel or is panel-shared): benign `g05 ss.50
ON` +8,909; `g25 ss.50 OFF` +5,399 (OFF rows are shared across panels — ONE defect, not three;
undissected); sparse `g25 ss.99 ON` +914; mild stranded-ON (+150…+330). The v1 ingredient-width
derivation targets all of these.

## THE `g25 OFF` DISSECTION — a mechanism worth keeping (2026-09-01)

Concentrated: FIVE exons in the two highest abundance blocks carried ~96 % of the +5.4k. The chain of
refuted hypotheses, each tested: refit feedback (REFUTED — final-sweep-only delivery is WORSE
everywhere, and the blind-row wins are legitimately bootstrap-mediated: true row information
trains the prior, which is where much of the whole-library win compounds); prior training on
message-shaped widths (REFUTED by the same arm). The REAL mechanism, named by a one-slot psi
factor replay at slot 258:

    landscape prior: peak at lam = −5.59 (truth), ~8 nats     row: cliff −18 nats AT that peak,
    flat above −4.6                                           SUM: flat to within 0.1 nat over
    [0.01, 1) — the factors ANNIHILATE, and psi's point estimate becomes the CENTROID of a
    vacuous plateau: f = 0.64 where truth is 0.003.

⭐ The lesson, stated generally: **two mildly-conflicting one-sided factors can cancel into a
flat posterior whose estimator wanders — the failure is not over-claiming but mutual
annihilation.** The row's cliff POSITION is a measurement (a 13-fragment face, 4-sigma lucky),
and its honesty is the ingredient width.

## v1 — THE INGREDIENT WIDTH (derived, per face, no constants)

The map's position variance by the delta method: `Var = trigamma(n_u+1/2) + trigamma(n_s+1/2)`
(the crossing count and the face's spliced count — the two measured map ingredients; the intron
row's own width is already inside the transported shape, so no double-count). Each face's
transported row is blurred by ITS OWN variance — a deep face blurs ~0, a 13-fragment face by
~0.28 nats, which softens the lucky cliff into a slope the prior's 8 nats overrule. Decisive-set
verdict: `g25 OFF` FIXED (12,969 → 8,070, slot 258 lands 0.005 vs truth 0.003); every win
retained (g50 ON 3,196; junction 2,351; sparse 2,714); the stranded-ON residue FLIPPED TO A WIN
(`g50 ss.99 ON` 1,544 vs silent 1,610). Standing shadow: `g05 ss.50 ON` (unchanged ~35.6k vs
silent 26.7k) — a DEFERRED-stratum row (unstranded × ON), the recorded misleading-neighbour
regime; not a 0.8.0 blocker, owed a dissection eventually.

## v1 CONFIRMED — three probe panels AND THE LADDER (2026-09-01)

Three-panel v1: every v0 win kept; `g25 OFF` residual +500; `g00 ss.50` zero controls now BEAT
silent on the benign panel; sparse panel worst degradation +582 (relay: 98k–188k harm on nine
rows); **junction panel: t_face_v1 is the best arm on every capture-ON row** — the bounded-damage
algebra held under the pathology built to break it.

**THE LADDER (silent / relay / rung1 / t_face_v1), the two halves apart:**

* HARM bar (ss.99): t_face_v1 ≤ silent on 7 of 8 rows (worst residue +406 = +0.5 % at
  `g05 ON`); at the WORST IN-SCOPE condition `g98 ss.99 ON`: silent 298,597 / relay 456,838 /
  rung1 294,676 / **t_face_v1 288,774**. `g98 ss.99 OFF`: 109,291 vs silent 126,467 (−14 %).
* WIN bar (ss.50): `g50 ON` **1,938,411** vs silent 6,141,095 (−68 %, within 8 % of the relay's
  entire machinery); `g98 ON` 3,077,782 vs 12,030,888 (−74 %; relay 2,433,908 keeps a 26 % lead);
  `g00 ss.50 OFF` (in-scope zero control) **303,826 vs silent 1,254,145** (−76 %; relay 58,840
  still ahead); in-scope `g50/g98 ss.50 OFF` −6 % / −15 %.
* Flip falsification catastrophic everywhere (3.9M at `g00 OFF`). Two-substrate sign agreement
  holds on both halves.

**Against the thread's original problem statement**: the shipped message layer ADDED +158,241
fragments (+53 %) at the worst in-scope condition; the two composition-transfer rungs REMOVE
9,823 (−3.3 %) there instead — while cutting the blind rows by 68–76 % and surviving hostile
probe designs the shipped relay cannot. The relay's remaining leads are the deferred-blind rows
and the benign-panel `g00`/`g98-ON` rows where its anchor still profits from probe-placement
luck (measured as such on the adversarial panels).

## RUNG 2 PROMOTED AND SHIPPED (2026-09-01, owner ruling "proceed, implement, ship")

Landed in `messages/transfer.py` as ONE mechanism beside rung 1: `face_is_licensed` (the
licence as a PURE predicate — extracted because the integration toy cannot falsify its terminus
branch, watched), `face_map_lambda` (monotone, flux-capped), `transport_row` (preimage read —
no Jacobian — then the delta-method ingredient width), and the exon-face block in
`TransferPolicy.prepare` (two licensed faces sum as independent witnesses; depleted/flat/
unlicensed faces deliver nothing). Rows ride EVERY sweep — the final-sweep-only citizenship was
built and REFUTED (true information training the prior bootstrap is where the blind-row value
compounds). Gates: 7, each written fail-first; watched perturbations: the ceiling dropped from
the map (fires), the width at half strength (fires), the terminus branch dropped from the
predicate (fires via the pure gate). Faithfulness: `policy_benchmark.py --policies silent
transfer` reproduces the harness TO THE FRAGMENT (`g98 ss.99 ON` 288,774; `g50 ss.50 ON`
1,938,411). Suite 3,747 passed / 8 xfail, 3,755 collected — accounted in CLAUDE.md (one case is
the owner's own `docs/dev/rename.md`).

# RUNG 3 — ARCHETYPE B case 1: the intergenic|exon EDGE (owner notes 2026-09-02)

Owner rulings (2026-09-02): the policy is INCOMPLETE and will not ship until terminus faces are
addressed; start with the intergenic|exon edge. Two design concerns RULED and absorbed: a biased
imputation mode cannot be rescued by a penalty (price bias vs benefit), and a single per-library
enrichment dispersion is the wrong model for a probed/unprobed/placement MIXTURE — the learned
enrichment structure is the LANDSCAPE's job, available only after a pass (the circularity the
refit bootstrap already resolves). Further rulings: `r ~ 1` rests on a DOCUMENTED TOOL-SCOPE
ASSUMPTION — probe panels do not target intergenic boundaries or trail off annotated transcript
ends — with discrete counting randomness around 1; and the substrate gains SINGLE-EXON
transcripts (owner-approved; Claude authored).

## Stage 0 (certified, no solver, all four panels)

* **(i) The structural premise holds EXACTLY**: certified RNA crossing intergenic|exon
  boundaries = 0.0 on every condition of every panel. (The recorded caution stands:
  `strand_evidence` refuses struct-lock CERTAINTY at the G1 boundary — ragged TSS/read-through
  are the real-data leak class — so the edge always speaks as a counted claim with width.)
* **(ii) The edge-level identity** (edge crossing density vs exon true gDNA density): **OFF
  capture 0.96–1.08 everywhere**; under capture a probe-design SPECTRUM — benign probed
  0.77–0.79, junction probed 0.16–0.19, sparse probed **0.004**, unprobed ~1.0–1.2, ladder ON
  pooled 0.43. The step is unknowable a priori; and the ratio never meaningfully exceeds 1 —
  **the bias has a KNOWN SIGN** (probes target exons, so an exon is at-least-as-enriched as its
  own edge crossing).
* **(iii) The census (ladder)**: 974 exons with edge+intron faces, **1,250 edge-only**
  (single-exon-like), 12,923 intron-only, 8,871 neither. 2,620 edges.

## THE REVISED DESIGN — the sign-certified INTERVAL claim, landscape-spanned

**Pass 1 (prior-free):** the edge delivers a soft ONE-SIDED row — the exon's observed gDNA
count is at least the edge's expectation (`sigma(lam)·n_e >= ~rho_edge·E_g^e`, cliff below,
flat above; softness = the two counting trigammas, which also carries the owner's r~1-by-
counting-randomness ruling). Sign-certified, so bias-immune under any probe design (an
under-read edge is a WEAKER true bound, never a false one); vacuous at `g00`; strong off
capture; harmlessly weak at probed exons under capture. The bias-cost half of the pricing
question drops to ~zero structurally.

**Refit sweeps (landscape fitted):** the claim becomes the MULTIPLICATIVE INTERVAL
`rho_e^obs ∈ [rho_edge·r_lo, rho_edge·r_hi]`, with `r_lo ~ 1` softened by counting and
`log r_hi` = the landscape's OWN learned span from the intergenic-anchored depleted level to
its support top (`log_rho[-1] − log rho_bg` — both derived per refit, no constants; the
landscape is a nonparametric `logP` over log-rate, so the span is read off the fitted object,
not off extracted "modes"). Off capture the fitted span is small -> tight intervals; under
capture it is decades -> weak upper sides at probed loci; a PASS-1-POISONED landscape (g00
unstranded false positives inflating the support) can only WIDEN `r_hi` — the safe direction.

⭐ **The `g00` kill returns through the MULTIPLICATIVE form, not through a learned narrow
span**: at `g00` the edge measures `rho_edge ~ 0` over real opportunity (structurally gDNA-only
crossing), and `[~0, r_hi·~0] = ~0` for ANY `r_hi` — the interval pins the exon at zero
regardless of enrichment, with softness = the edge count's own Gamma posterior tail. The same
immunity as before (`mu ~ 0` times any step is `~ 0`), now carried by both sides.

**Documented tool-scope assumption (to MOVE to `DESIGN.md` when this rung ships)**: probe
panels are assumed not to target intergenic sequence at gene boundaries nor to trail past
annotated transcript ends; `r ~ 1` (edge enrichment ≤ exon enrichment) rests on it, and a
panel violating it voids the edge claim's sign certificate.

## The substrate: THE MONO BLOCK (authored 2026-09-02, owner-approved)

20 single-exon transcripts appended to the twin block: one 1 kb exon, 20 kb intergenic, all
`+`, FOUR types × the five blocks — `mono` (unprobed, expressed), `capmono` (probed),
`monosilent` / `capmonosilent` (mrna = 0, the false-positive controls) — mrna by the standard
ladder, nascent 0 by ruling (nothing to model on a single exon). No sj anywhere, so the ONLY
message source is the gene edge. The generated chromosome was extended to **1.5 Mb** to hold
them (`GENOME_LENGTH`; ⚠ the blank contig also went to 1.5 Mb — an unintended-but-symmetric
side effect of the edit, kept; gDNA densities shift accordingly and the full rebuild
re-baselines). The adversarial probe generator gained the mono types (sparse 55 probes;
junction unchanged at 30 — a mono gene HAS no junction, so on the junction panel even "probed"
mono genes are uncapturable, itself a useful stress). ⛔ Every derived test-chromosome artifact
went stale with the GTF edit; the full rebuild (index + all seven panels + certification) ran
the same day.

## Next

1. Rung-3 prototype once panels certify: arms {off, pass-1 bound only, +landscape interval},
   flip; falsifiers: sparse probed exons (bound must stay harmless), `g00` kill by refit 1,
   mono-silent genes clean, stranded ladder rows unmoved.
2. Then the remaining terminus cases, opposite strands, overlapping isoforms — each a rung.
3. The deferred `g05 ss.50 ON` dissection stays parked.

## RUNG 3 CLOSED — the LOWER BOUND landed; the ceiling REFUSED (owner ruling, 2026-09-02)

The owner's ruling after the landscape-poison study: the enrichment-ceiling upper side is
over-engineering — "take the win, keep the policy clean and simple, accept the error." Landed in
`messages/transfer.py` as `edge_bound_row` + the edge block: the profile likelihood
`sup_{s>=1} Pois(n_b; c/s)` (0 wherever the exon's implied gDNA count covers the edge's; the
edge count's own one-sided Poisson tail below; IDENTICALLY VACUOUS at `n_b = 0` — the
near-zero-row artifact family is barred by law, not by threshold). 8 gates green, three
perturbations watched firing (two-sided bound, licence dropped, vacuity dropped). The A/B
(src-before vs src-after, the official instrument): `g50 ss.50 ON` 2,844 -> 2,734; stranded
1,772 -> 1,753; **`g00 ss.50 OFF` byte-unchanged at 14,602 — the accepted error stays visibly
accepted**. New-substrate sweep: unstranded 14/20 wins, stranded worst 1.02x. Suite 3,748 / 8
xfail (3,756 collected, accounted).

⚠ The harness's earlier `g00` "wins" (9,015; the flip-identical 33) are attributed to
harness-vs-src implementation divergence amplified by flat-posterior centroid sensitivity —
prototype artifacts, not mechanism; the src numbers above are the record.

## THE EXPOSED SYSTEMIC ISSUE — raised as `ISSUES: gdna-landscape-trains-on-false-positives`

The landscape-poison study's findings moved to the issue log (the MOVE rule): 100 % fiction
training at `g00 ss.50`; the annotation-vs-kappa gate mismatch; the reliability weight blind to
composition blindness; ⛔ the naive exclusion REFUTED by measurement (it starves the bootstrap
that generalizes message-delivered truth — `g50 ss.50 ON` 2,691 -> 56,422); and the separate
grid-vs-mass span-read defect. The open question and candidate mechanisms live in the entry.

# RUNG 4 — THE TERMINUS MECHANISM (holes ① and ②): derivation + stage 0 (2026-09-02)

## What the ladder census re-framed (chain of `g50 ss.50 OFF`, no solver)

Hole ① is two holes. Of the 12,811 exon|exon boundaries, **7,604 are internal TERMINI** (a TSS/TES of one
isoform inside another's exon; 950 k crossing mass) and 4,838 are alternative splice sites (licensed by
the rung-2 predicate; 458 k). 10,259 exons (1.05 M) are unreachable by ANY chain of licensed faces —
walled by terminus faces (11,630 of their faces). Hole ② (`exon|intron[term]`, 2,979 slots) is the same
mechanism at an intron flank. Strand-change faces are small (246 exon|exon + 888 exon|intron) and stay a
later rung, as the owner ruled (both-stranded loci last).

## Stage 0 — THE DIRECTED LICENCE, certified on ladder truth

A terminus boundary's terminating transcripts cover exactly one flank (the INSIDE); the direction is
read off the flag (TSS+ / TES− extend genomic-right, so outside = left flank; TES+ / TSS− extend left).
Noise-subtracted logit-composition gap (`gap² − Σ trigamma(n+½)`), mixed slots ≥ 20 fragments,
single-direction boundaries:

| condition | class | n | OUTSIDE excess | INSIDE excess |
|---|---|---|---|---|
| g50 ss.50 OFF | exon\|exon | 268 | −0.023 | 0.667 |
| g50 ss.50 OFF | exon\|intron | 176 | 0.003 | 2.716 |
| g50 ss.50 ON | exon\|exon | 312 | 0.021 | 1.802 |
| g98 ss.50 OFF / ON | exon\|exon | 109 / 187 | −0.016 / −0.021 | 0.354 / 0.786 |

**The OUTSIDE flank shares the crossing's composition to within counting noise on every condition; the
INSIDE flank does not.** The gDNA LEVEL across an `exon|intron[term]` face under capture is destroyed
(median log ratio 5.9, the cliff) while the composition survives.

## The derivation — three mechanisms, one currency, one source factor

Notation as rung 2: at boundary `b` the crossing count `n_b` with opportunities `A_g^b, A_r^b`; at a region
`R` the count `n_R` with `E_g^R, E_r^R`; λ = log-odds of the gDNA COUNT share at a slot; `s` a measured
spliced density.

**(i) The OUTSIDE transfer (composition, either direction).** `comp(b) = comp(O)` is the same
shared-population law as rung 1's intron→boundary pair, so the hop is the s = 0 face map (a pure
opportunity shift): `λ_b = λ_O − log(E_g^O/A_g^b) + log(E_r^O/A_r^b)` and its inverse. For hole ②
(`O` = intron) rung 1 ALREADY delivers the intron row to `b` — its pair detection never checked flags, so
`exon|intron[term]` boundaries are served. What is new is `O` = exon: `b` receives `O`'s composition —
and `O`'s composition is not `O`'s own belief (that would be a relay of a destination belief, the
foundation's reserved case) but **the source factor `O` was itself reached by, transported one map
further**. So rung 4 is built on

**(ii) THE COMPOSED TRANSPORT — one source factor through a chain of licensed maps.** The nearest
intron factory's row is carried region→boundary→region through licensed faces, each step a monotone
map, each measured ingredient adding its delta-method width (rung 2's `trigamma(n_u+½) + trigamma(n_s+½)`
per face; widths sum along the chain — a deep chain is honestly wide). Two step kinds:
* BOUNDARY → REGION: rung 2's face map (spliced route `s` joins on the entering side).
* REGION → BOUNDARY: the crossing carries the region's unspliced population MINUS the mature RNA that
  splices OUT at `b` in the direction of travel. Where nothing leaves — the sj at `b` is an ACCEPTOR
  w.r.t. travel (genomic-order flag `ACC_s` when travelling right, `DON_s` when travelling left; the
  boundary at the HIGH end of an intron is where spliced RNA JOINS rightward travel, at the LOW end it
  LEAVES, on either strand) — the step is the identity (s = 0 inverse map). Where mature RNA LEAVES, the
  step is a SUBTRACTION of the measured splice-out flux from the region's RNA density — the relay's
  classic trouble spot (`TRAPS: zero-the-precision-with-the-value`); rung 4 REFUSES leave steps
  (silence) and prices what that costs (below). ⭐ On `altss` every exon is reachable JOIN-ONLY from the
  appropriate side (the core at depth 2 from either intron), so the structure that indicts the licensed
  multi-hop needs no splice-out.

**(iii) THE INSIDE BOUND (level, one-sided).** `I`'s population = crossing + the terminating
transcripts' RNA (unmeasured), so composition cannot cross; gDNA CONTINUITY can: `ρ_g^I = s·ρ_g^b`,
`s` the unknown enrichment step. With `g_b(λ_b) = σ(λ_b)·n_b` the crossing's implied gDNA count and
`c(λ_I) = σ(λ_I)·n_I·A_g^b/E_g^I` the count `b` would show under continuity, the claim on `I` is the
PROFILE over the nuisance `s ≥ 1` AND over the crossing's own composition:

    L_I(λ_I) = max_{λ_b} [ L_b(λ_b) + sup_{s≥1} logPois( g_b(λ_b) ; c(λ_I)/s ) ]

— the source row's RUNNING MAXIMUM read through the level map, softened by the crossing count's
one-sided Poisson tail. It reduces EXACTLY to rung 3's `edge_bound_row` when `L_b` is the pure-gDNA
point (an intergenic crossing), so the edge rung is its special case; a peaked `L_b` becomes one-sided
(the upper side is dropped — the refused ceiling never enters); a cliff row keeps its cliff. ⚠ The sign
certificate `s ≥ 1` (destination at least as enriched as the crossing) holds for `exon|intron[term]` by the
edge's documented tool-scope assumption (probes target exons; the face is where they trail off) and is
NOT certified for `exon|exon[term]` under capture (both flanks are exons of the host; the ladder's median
log level ratio there is +0.14–0.17 ON, ≈ 0 OFF) — price it on the three probe panels before believing it.
Nothing goes from `I` to `O`: the gDNA level the other way is an UPPER bound on `O` (the refused
direction) and `O` is reached from its own side.

**What the rules reach on the ladder (chain census, `g50 ss.50 OFF`):** rungs 1–3 today leave 11,350
exons (1.11 M) unreached. Under (i)+(ii) join-only +(iii): **COMP 12,018 exons (701 k) · BOUND-only 5,392
(388 k) · UNREACHED 6,608 (841 k)**; non-edge boundaries COMP 21,524 · BOUND 4,071 · UNREACHED 6,826
(1.04 M crossing mass). Allowing LEAVE steps would move 1,741 exons (177 k) from unreached to
composition — the measured price of refusing splice-out, to weigh later. The remaining unreached are
exons between two termini facing away from each other (2,167; 190 k), leave-only sj chains (1,232;
165 k), and edge+terminus exons the edge rung already bounds (802; 63 k — not counted here).

**On the new test chromosome** the same census says the substrate indicts exactly these: every exon is
reached by composition except the mono block (edge-only, rung 3) and the 10 `nest` walls (BOUND-only,
3,840 mass); the `altss` cores sit at depth 2 (join-only); 30 exon|exon terminus boundaries and 10
`exon|intron[term]` faces exist, on both strands.

## The plan (DERIVE done → DESIGN → PLAN → PROTOTYPE → A/B → src)

1. Baseline on the rebuilt substrate FIRST (`policy_benchmark.py --panel test --policies silent
   transfer`, all seven panels certified), before any code.
2. Prototype through the `message` seam, ONE mechanism per arm, each with its flip falsification:
   arm A = (i)+(ii) join-only composed transport (reaches `altss` cores and the `altstart`/`nest`
   OUTSIDE pieces and their boundaries); arm B = A + (iii) the inside bound (reaches `instart`'s first
   exon, `altstart`'s inside piece, the `nest` walls). Score per type × strand × probe panel; the
   stranded rows must stay ≤ silent, the zero controls and `silent`/`monosilent` genes clean.
3. Confirm on the ladder, both halves apart, all three probe worlds; then the promotion ruling.

## RUNG 4 — the PROTOTYPE on the rebuilt benign panel (2026-09-02, same day)

**The new baseline first** (`policy_benchmark.py --panel test --policies silent transfer`, the
rebuilt 85-gene chromosome at 480 k fragments): unstranded `transfer` beats silence on **14/20**, worst
row **1.21×** (`g00 ss.50 ON`, the deferred zero control: 25,820 vs 21,422); stranded **9/10**, worst
**1.01×**. The in-scope unstranded OFF rows sit at 0.92–1.01× and are near the counting floor
(whole-library ~7–9 k fragments on 480 k) — the region axis is slightly WORSE than silence at
`g25/g50 ss.50 OFF` (7,808 vs 7,653; 6,805 vs 6,425) while the boundary axis is better. The blind
capture-ON rows keep the large wins (0.05–0.20×). Two rows (`g25 ss.70/.99 OFF`) certify at
COMPOSITION level only on every panel — the gdna-field-uniformity gate flags 1–2 of 981 z-scored slots
(a naive per-ref z reproduces no flag; the OFF simulations are shared across probe panels so it is one
event, not three). Recorded, not chased.

**The prototype** (`rung4_proto.py`, scratchpad; a `TransferPolicy` subclass installed by patching the
module attribute `rigel.calibration.calibrate.TransferPolicy` — ⚠ `from rigel.calibration import
calibrate` returns the FUNCTION, and patching it silently does nothing: every arm read identical to
`transfer` until the module was patched). Arms: `ident` (MULTI off, BOUND off — **byte-identical to
`transfer` on every condition and type**, the identity gate), `multi` (the composed join-only transport
+ the outside→terminus composition), `bound` (multi + the inside bound), `flip` (the inside bound with
the inequality reversed).

Whole-library |err| (fragments), `g50` rows:

| row | silent | transfer | multi | bound | flip |
|---|---|---|---|---|---|
| ss.50 OFF (in scope) | 6,983 | 7,064 | 7,161 | 7,269 | 7,133 |
| ss.50 ON (deferred) | 207,139 | 27,892 | **25,445** | 25,806 | 30,897 |
| ss.99 OFF | 6,209 | 6,050 | 6,068 | 6,078 | 6,064 |
| ss.99 ON | 4,155 | 4,119 | 4,115 | 4,112 | 4,359 |

Per type, the mechanism's own structures at `ss.50 ON`: `capaltss` 5,471 → 4,640 (multi), `capaltstart`
2,666 → 2,428, `capnest` 2,870 → 2,638, `capinstart` 8,979 → 8,732 (bound 8,795; **flip 14,901** — the
falsification fires where the structure is hole ②); stranded rows unmoved (±3). At `ss.50 OFF` multi
costs `altss` +44 and a few fragments elsewhere.

**Dissections (per slot, truth beside each arm):**
* `altss` OFF, block 1: the core (truth f_g 0.094, 363 fragments) reads 0.127 under `transfer` and 0.155
  under `multi` — the transported claim over-reads by ~20 fragments. The crossing's truth (0.108, 111
  fragments) and the T-only piece's (0.136) differ by counting noise; in per-opportunity terms the piece
  and core have the SAME RNA density and the composition gap is the gDNA count's Poisson noise (19 vs
  34). A 10-molar gene at 480 k is at the floor; the in-scope OFF cost is noise-level in absolute
  terms but it is a cost.
* `capaltstart` ON: the inside piece's `transfer` estimate already sits ABOVE truth (block 1: 0.708 vs
  0.635; block 2 [the − gene; the inside piece is the mirrored `[8000, 8500)`]: 0.769 vs 0.749) because
  rung 2 reaches it from its licensed intron face — so a LOWER bound is inert there and the flip's
  "win" is an accident of direction. ⭐ The directed licence reads correctly on the − gene: the inside
  piece carries more RNA, the crossing matches the outside flank (0.877 vs 0.873).
* `capnest` ON — the walled exon, the pure case: silence is blind (0.001) but **`transfer` already
  reads 0.666 vs truth 0.630 — the refit landscape prior, trained by the other rows, reaches the wall**.
  `bound` 0.665, `multi` 0.659, `flip` 0.460 (fires). The inside bound's marginal value is ~zero where
  the bootstrap prior already generalizes.

**Verdict so far**: the mechanisms behave exactly as derived (identity holds, both falsifications fire
where the structure is theirs, direction correct on both strands), the composed transport wins ~9 % on
the deferred blind row with ≤ 3 fragments of stranded movement, and the inside bound is INERT on this
substrate (every inside flank here has either a licensed witness or the prior). ⛔ Nothing here says
"ship": the in-scope unstranded OFF rows show a small cost, and the population that motivated the rung
(the ladder's 10,259 walled exons) is not on this chromosome in bulk. The LADDER run
(`rung4_proto.py ladder transfer,multi,bound,flip all`, ~1 h) is the decision.

**The adversarial probe panels (whole-library |err|, `g50` capture-ON rows + the `g00 ss.50 ON` zero control):**

| panel · row | silent | transfer | multi | bound | flip |
|---|---|---|---|---|---|
| SPARSE · ss.50 ON | 212,933 | 55,154 | **54,114** | 59,489 | 61,444 |
| SPARSE · ss.99 ON | 4,056 | 3,916 | 3,909 | 3,912 | 8,606 |
| SPARSE · g00 ss.50 ON | 55,416 | 56,101 | 56,101 | 56,101 | 56,101 |
| JUNCTION · ss.50 ON | 119,379 | 19,645 | **13,389** | 13,255 | 12,693 |
| JUNCTION · ss.99 ON | 4,510 | 3,987 | 3,993 | 4,001 | 4,252 |
| JUNCTION · g00 ss.50 ON | 103 | 102 | 102 | 102 | 102 |

⭐ `multi` is robust on both hostile designs (−2 % sparse, −32 % junction on the blind row; stranded
within 0.2 %; zero controls untouched). ⛔ **The inside bound HARMS on the sparse panel (+8 %)** — the
sign certificate `s ≥ 1` is not guaranteed at an exon|exon terminus (the derivation said so), and on
sparse probes the crossing can be MORE enriched than the inside piece, turning the "lower bound" into
an over-claim. So the inside bound is refuted as a shipping candidate by measurement: inert where a
witness or the prior exists, harmful where the enrichment sign flips. **The candidate for `src/` is the
composed join-only transport ALONE** — one mechanism — pending the ladder run.

## RUNG 4 — THE LADDER (16 conditions, silent / transfer / multi / bound / flip; 2026-09-02)

Ratios are against the SHIPPED `transfer` (rungs 1–3); silence for scale.

| unstranded (ss.50) | silent | transfer | multi | bound | flip |
|---|---|---|---|---|---|
| g00 OFF (in-scope zero control) | 1,254,145 | 303,826 | **285,785 (0.941×)** | 285,785 | 251,989 |
| g00 ON | 454,560 | 344,331 | 338,890 (0.984×) | 338,890 | 326,748 |
| g05 OFF | 50,435 | 51,074 | 51,149 (1.001×) | 51,313 (1.005×) | 51,145 |
| g05 ON (deferred) | 518,535 | 220,181 | 221,455 (1.006×) | 217,384 (0.987×) | 248,056 |
| g50 OFF | 155,660 | 147,736 | 148,661 (1.006×) | 149,613 (1.013×) | 148,901 |
| g50 ON (deferred) | 6,141,095 | 1,456,056 | 1,445,494 (0.993×) | **1,388,167 (0.953×)** | 1,811,445 |
| g98 OFF | 165,259 | 139,694 | 139,326 (0.997×) | 138,761 (0.993×) | 140,965 |
| g98 ON (deferred) | 12,030,888 | 2,522,236 | 2,479,621 (0.983×) | **2,385,564 (0.946×)** | 3,213,817 |

| stranded (ss.99) | silent | transfer | multi | bound | flip |
|---|---|---|---|---|---|
| g00 OFF | 30,606 | 18,960 | 18,305 (0.965×) | 18,305 | 17,775 |
| g00 ON | 20,787 | 18,197 | 18,117 (0.996×) | 18,117 | 17,974 |
| g05 OFF | 44,519 | 44,155 | 44,205 (1.001×) | 44,192 | 44,239 |
| g05 ON | 85,294 | 85,302 | 85,404 (1.001×) | 85,316 | 90,185 |
| g50 OFF | 125,634 | 117,679 | 118,614 (1.008×) | 118,551 (1.007×) | 118,864 |
| g50 ON | 260,629 | 244,163 | 244,035 (0.999×) | 243,677 | 333,588 |
| g98 OFF | 126,467 | 108,890 | 108,637 (0.998×) | 108,359 | 109,365 |
| g98 ON | 298,597 | 288,834 | 288,423 (0.999×) | 288,138 | 498,711 |

**Reading it, halves apart.** `multi` (the composed join-only transport) is NEUTRAL on the ladder: 5/8
rows better in each half, worst 1.006× unstranded and 1.008× stranded, one real move (−6 % at the
in-scope `g00 ss.50 OFF` zero control, which the shipped policy had left at 303,826). `bound` adds a
−5 % on the two large DEFERRED rows (`g50/g98 ss.50 ON`) and costs +1.3 % at the in-scope
`g50 ss.50 OFF`; on the sparse-probe panel it costs +8 %. `flip` is catastrophic where the bound acts
(deferred rows +25 %, `g98 ss.99 ON` +73 %) — the falsification fires on the ladder too. At `g00`
multi ≡ bound exactly (no background ⇒ flat intron rows ⇒ nothing to transport; the g00 gains come from
the outside→terminus and join hops feeding the prior bootstrap).

**Verdict.** The 1.05 M fragments of walled exon mass did not turn into whole-library error movement:
those exons are already served by the refit prior about as well as a transported row serves them. So
rung 4 CLOSES holes ① and ② structurally — every reachable terminus and exon|exon face now receives a
derived claim, join-only, no constants — at a measured cost of ≤ 0.8 % on any in-scope row and a −6 %
zero-control gain, but it does not buy a win. The inside bound is REFUTED for shipping (in-scope cost,
sparse-probe harm, its only gain on the deferred stratum). ⭐ The owner's call: ship `multi` as the
completion of the policy (structural coverage at neutral cost) or record it as measured-and-unneeded and
call the holes closed by the prior's reach.

# THE RESET (owner, 2026-09-02, late) — and ITEM 1: the exon → intron|exon boundary message (DERIVATION DRAFT)

The owner's rulings: rungs 1 and 2 are unfinished (the exon-side messages were nullified and are
owed), they come first, one message at a time; rung 4 grows one structure per step (`altstart` only
now); the tracker is `MESSAGE_RUNGS.md`. Every piece is derived AND taught.

## Item 1 — the setting, single isoform U on +, one intron|exon boundary b (an ACCEPTOR of U)

    intron I  |  b  |  exon E

Unspliced populations (spliced fragments are certified RNA and live in their own lane):

| slot | gDNA | U pre-mRNA (nascent) | U mature (fragments fully inside E) |
|---|---|---|---|
| I (contained) | yes | yes | no |
| b (crossing) | yes | yes | **no** — a mature molecule reaches E by the sj, never across b |
| E (contained) | yes | yes | **yes** |

Rung 1 (I → b) is exact because I and b share a population. Item 1 is the OTHER witness: E → b. E holds
one component more than b, so the message must REMOVE it. Two forms, in order of how much they assume.

**Form A — the inequality (assumes nothing beyond E's composition).** Mature RNA is non-negative, so the
crossing's RNA density is at most E's: `ρ_r^b ≤ ρ_r^E`, while gDNA is continuous, `ρ_g^b = ρ_g^E` (up to
the enrichment step between a face and an interior, which cancels in a ratio only if it is
component-blind at each locale — the same premise rung 2 already carries). Hence in log-odds

    λ_b ≥ λ_E + [ log(A_g^b/A_r^b) − log(E_g^E/E_r^E) ]

— a LOWER bound on the boundary's gDNA share: cliff below, flat above. It is the running maximum of E's
row read through the opportunity shift — the same one-sided-row currency rung 3 ships and the rung-4
inside bound used. It needs no spliced measurement at all.

**Form B — the subtraction (assumes the spliced route measures E's mature density).** Remove the mature
share: `ρ_r^b = ρ_r^E − s`, where `s` is the face's route-summed certified spliced density (rung 2's `s`).
Then `λ_b` is a point claim. ⚠ This is a SUBTRACTION of a quantity measured in a different locale (the
junction's) from one measured in E's, so a wrong `s` enters in FULL — the adversarial panels measured the
route identity at 0.003× (sparse probes) and 25× (junction probes) against 1.00–1.03 on the benign one.
Rung 2 survives that because it ADDS `s` (a wrong `s` only mis-caps a claim); form B would not. It is
also the relay's `TRAPS: zero-the-precision-with-the-value` family: near `ρ_r^E ≈ s` the estimate's
variance diverges, and the honest statement there is again one-sided.

**Where the message has value.** At a boundary whose intron is well measured, rung 1 already solves it
and E's witness adds little. The value is where the intron is thin — capture-ON unprobed introns hold
~10–20 fragments — and E is probed and, on stranded data, sharply solved by its own strand evidence. So
the first measurement is: per boundary, the intron row's width against E's, on capture-ON rows.

**The falsification** for either form: reverse the inequality (form A) or add instead of subtract (form
B); both must make the boundary's error jump on the rows where the message acts.

**Recommendation (for the owner's ruling):** derive and measure form A first — it is exact, needs no
`s`, and its residual (how far the crossing sits below E's share, i.e. the mature share of E) is a
measured number that then tells us what form B would be worth and where it is safe.

## THE BASELINE ON THE TRIMMED SUBSTRATE (55 genes: twin + mono + `altstart`; 480 k fragments; 2026-09-02)

`policy_benchmark.py --panel test --policies silent transfer`, the shipped rungs 1–3 against silence —
the number every item in `MESSAGE_RUNGS.md` is measured from. Unstranded: transfer beats silence on
**14/20**, worst **1.16×** (`g00 ss.50 ON`, the deferred zero control, 64,143 vs 55,290); stranded
**8/10**, worst **1.00×**. In-scope unstranded OFF rows: `g05` 1.00×, `g25` **1.06×** (12,738 vs
12,042 — the one in-scope row where the shipped policy costs), `g50` 0.99×, `g98` 0.92×. Blind
capture-ON rows: `g05` 0.08×, `g25` 0.76×, `g50` 0.05×, `g98` 0.07×. Full table in the session's
`baseline_trimmed_silent_transfer.txt`; every earlier test-chromosome number is on another substrate.

## Item 1, form B — THE OWNER'S ARITHMETIC, and what it reduces to (2026-09-02)

Owner design (`message_notes.md`): enrichment ratio `e = T_E / (U_b + S_b)` (exon total density over
the boundary's spliced + unspliced), rescale the boundary's spliced density into the exon's frame, subtract
it from the exon's RNA, rescale the leftover back. Done in symbols, with `f_E` the exon's gDNA share:

    mature at E's scale       m = e·S_b
    leftover RNA at E's scale (1 − f_E)·T_E − e·S_b
    boundary share            f_b = f_E·T_E / (T_E − e·S_b) = f_E / (1 − S_b/(U_b + S_b))

⭐ **`e` cancels.** `e·S_b/T_E = S_b/(U_b+S_b)` is the boundary's OWN spliced share, so

    f_b = f_E · (U_b + S_b) / U_b            (capped at 1)

— every ratio is formed within one locale (rung 2's law), no level ever crosses the face, and the
rescale/subtract/rescale is exactly the INVERSE of the shipped splice-in face map (`face_map_lambda`
solved for `f_b`). So at solve time the boundary evaluates the exon's row AT the splice-in map:
`L_b(λ_b) = L_E(M(λ_b))` — the same function, read the other way, no inversion code.

**Certified check (twin block, `g50 ss.99 OFF`, every splice-out face with ≥ 5 crossings):** the
prediction from the exon's TRUE share and the face's own counts tracks the boundary's true share across
the ladder — nascent genes 0.09/0.09, 0.20/0.21, 0.22/0.20, 0.62/0.61, 0.81/0.77; pure-gDNA faces
(`clean`/`cap`, true 1.000) read 0.65–1.00, the shortfall being the Poisson noise of `U_b` (7–18
fragments), the smallest count in the formula and its dominant variance term. Table in the session's
`splice_out_check.py`. So form B is UNBIASED in expectation and its honest width is the delta-method
variance of `S_b/(U_b+S_b)` through the map, which diverges as `U_b → 0` (the row saturates at "pure
gDNA" — one-sided, never a near-zero row). The premise it carries: spliced and unspliced fragments at
the SAME face share capture affinity (the junction-probe pathology inflates `S_b/U_b` and over-claims
gDNA at the boundary; to be priced on the adversarial panels, never assumed away).

On unstranded data the exon's own row is flat and the message is silence; its value is on stranded
rows where the exon is sharply solved and the intron is thin (capture-ON). Falsifiers: the nullified
form (`f_b = f_E`) and the reversed map must both lose where `S_b/U_b` is large.

## Item 1 — THE PROTOTYPE (form B through the `transfer` seam; `item1_proto.py`, 2026-09-02)

The exon publishes its OWN evidence: its strand log-likelihood over the grid (variance frozen at its
incoming belief, the solver's own count-zero-information freeze), gated by the solver's derived strand
DEADBAND (`own.tau_lam > 0`; inside the noise floor the exon says nothing — no constant). The boundary
reads that row at the splice-in face map, `L_b(λ_b) = L_E(M(λ_b))`, widened by the face's counting
variance `trigamma(U_b+½) + trigamma(S_b+½)`. Arms: `transfer` (src), `formb`, `null` (`f_b = f_E`,
splice-out ignored), `rev` (the exon claimed to hold LESS RNA than the crossing). ⚠ A first `rev` arm
accidentally re-implemented the correct map and tracked `formb` to the fragment — a falsifier that
cannot fail is not a falsifier; it was rewritten as the true reversal before anything was read.

**Benign panel, whole-library |err| (boundary axis in the last column):**

| row | silent | transfer | formb | null | rev | formb/transfer | boundary tr → fb |
|---|---|---|---|---|---|---|---|
| g05 ss.70 ON | 1,625 | 1,667 | **1,408** | 3,405 | 3,784 | 0.845 | 745 → 486 |
| g05 ss.99 ON | 1,433 | 1,367 | 1,357 | 3,488 | 3,818 | 0.993 | 237 → 227 |
| g25 ss.70 ON | 4,517 | 4,098 | **3,952** | 12,083 | 14,553 | 0.964 | 837 → 691 |
| g25 ss.99 ON | 2,403 | 2,216 | 2,247 | 10,122 | 12,553 | **1.014** | 522 → 553 |
| g50 ss.70 ON | 5,433 | 5,240 | **5,034** | 15,577 | 19,746 | 0.961 | 1,071 → 865 |
| g50 ss.99 ON | 3,356 | 3,288 | **3,189** | 14,160 | 18,909 | 0.970 | 613 → 514 |
| g98 ss.70 ON | 6,314 | 6,402 | **6,098** | 6,525 | 6,318 | 0.953 | 901 → 597 |
| g98 ss.99 ON | 3,763 | 3,719 | **3,086** | 4,227 | 5,002 | 0.830 | 1,082 → 448 |
| every capture-OFF row | — | — | ±4 | worse | worse | 1.000 | — |
| every ss.50 row (gated) | — | — | ≡ transfer | ≡ | ≡ | 1.000 | byte-identical |
| every g00 row | — | — | ≡ transfer | ≡ | ≡ | 1.000 | — |

The message acts exactly where derived — stranded and part-stranded capture-ON rows, where the probed
exon is sharp and the unprobed intron thin — and is silence off capture, at the zero controls and on
unstranded data. Both falsifiers fire wherever it acts (null 4–5×, rev 5–6× on the boundary axis).
Per boundary (`cap`, `g98 ss.99 ON`, ~1,000 crossings each, truth 1.000): residuals 12 → 2, 10 → 1.

**The adversarial panels (stranded capture-ON rows):** junction 0.926 / 0.883 / 0.961 (wins);
sparse **1.124** at `g50 ss.99 ON` (boundary 482 → 885), 0.993 at `g98`, 1.011 at `ss.70`.

**The one benign cost, dissected (`capnasc`, `g25 ss.99 ON`):** the nascent-bearing probed gene's
boundaries are OVER-claimed (block 3: 0.909 vs truth 0.814 where `transfer` read 0.816; 0.974 vs
0.851). Cause: at a probed face the unspliced crossing sits on the enrichment SHOULDER (rung-2 stage 0
measured the contiguous ratio at 0.775 for `capnasc` ON) while the spliced fragments, whose bodies lie
in two interiors, are fully enriched — so `S_b/U_b` over-reads the mature share and the subtraction
removes too much RNA. Nascent-free genes hide it (their boundaries are pure gDNA and the claim
saturates at the truth). The sparse panel is the same premise at full strength. ⭐ This is the
IMPUTATION COST the derivation named, now measured: the counting width alone under-states the
message, and the missing term is the premise "spliced and unspliced fragments at the same face share
capture affinity".

**Next (the honest width):** the premise must be FITTED at runtime, never a constant — the house
pattern is `messages/variance.premise_logvar` / `splice_in_premise_logvar` (method of moments: the
observed spread of a measured contrast minus the counting variance it was measured with, floored at
0). The natural contrast here is the TWO WITNESSES at a licensed intron|exon boundary: the intron's row
and the exon's form-B row are independent statements of the same composition, and the part of their
disagreement that counting cannot explain is the premise variance, added to the message's blur. Owed:
derive the estimator, measure it per condition (it should read ~0 off capture and grow with the
shoulder), then re-run the four arms.

## Item 1 — THE TRANSFER VARIANCE, derived fresh (owner: reuse nothing older; 2026-09-02)

**One measured ratio carries the message**, `ρ = S_b/U_b`; beyond the exon's own row, every error of the
message is an error in `log ρ`: (1) COUNTING, `Var(log ρ) = 1/S_b + 1/U_b`, shrinking with depth; (2) THE
PREMISE, `log ρ_obs = log ρ_true + log a` with `a` the spliced-vs-unspliced enrichment ratio AT THE FACE
(`a = 1` is the premise), which does not shrink with depth and is a property of probe geometry.

**The two-witness estimator of the premise** (fitted per library, no constant): at every licensed
intron|exon face where the exon AND the boundary both have a live strand row, the two strand solves
IMPLY a ratio, `ρ_imp = f_b(own)/f_E − 1`, and the face MEASURED one, `ρ_obs = S/U`; their difference is
`log a`. With each witness's counting variance carried by the delta method,
`v_prem = max(0, Var_w(d) − mean_w(v_count))`, `w = 1/v_count`, and the weighted mean of `d` is the
BIAS. ⚠ Two coordinates were tried and refused by their own numbers: logit-of-`f_b` (both witnesses'
variances explode at the vertex: counting terms of 40–65 logit²) and a first version against the INTRON
witness (one-sided cliff rows have no usable variance). `log ρ` is the coordinate where everything is
tame, and only faces with `ρ_imp > 0.05` and ≥ 5 fragments each side can speak about `a`.

**What it measured:**

| panel · row | faces | mean d = log a | Var_w(d) | mean_w(v_count) | v_prem |
|---|---|---|---|---|---|
| benign g50 ss.99 OFF | 95 | +0.12 | 0.48 | 1.13 | 0 |
| benign g98 ss.99 OFF | 60 | −0.00 | 0.10 | 0.29 | 0 |
| benign g25 ss.99 ON | 57 | **+0.27** | 0.016 | 0.010 | 0.005 |
| benign g50 ss.99 ON | 54 | **+0.29** | 0.020 | 0.006 | 0.014 |
| junction g50 ss.99 ON | 68 | **+0.78** | 0.031 | 0.020 | 0.011 |
| sparse g50 ss.99 ON | 8 | −0.15 | 0.008 | 0.020 | 0 |

⭐ Off capture the premise is clean (`a ≈ 1`, no spread). Under benign capture `a ≈ e^0.28 = 1.3` —
EXACTLY the shoulder rung 2's stage 0 measured from the other side (contiguous ratio 0.775 = 1/1.29):
spliced fragments, whose bodies lie in two fully-tiled interiors, are enriched ~1.3× over the crossers
on the shoulder. On junction probes `a ≈ 2.2`. And the SPREAD is tiny everywhere: on these panels the
premise error is a BIAS with a known sign, not a variance — so a fitted widening cannot fix it
(`formb_w` ≡ `formb` to a few fragments), and correcting the bias would be a level-like fudge: the
owner's ruling, not the derivation's. Recorded, not corrected.

**The width that was actually wrong: the coordinate of the blur.** The prototype applied the counting
width as a uniform blur in `λ_b`, but noise in `log ρ` moves `λ_b` by `σ_b/(1−f_b)` per unit — unbounded
near the pure-gDNA vertex — so thin faces near the vertex (the sparse panel's 1–5-crossing faces)
received a claim far more confident than their counts allow. The honest width is a MARGINAL over the
measured ratio: `log ρ ~ N(log ρ_obs, 1/S + 1/U)`, nine quantile nodes, the exon row read at each
node's map, likelihoods averaged (`formb_m`). Wide exactly where the map is sensitive, no wider
elsewhere, no constant (the node count is quadrature resolution like `n_grid`).

| row | transfer | formb (uniform blur) | **formb_m (marginal)** |
|---|---|---|---|
| benign g25 ss.99 ON | 2,216 | 2,247 | **2,171** (the one cost is gone) |
| benign g50 ss.99 ON | 3,288 | 3,189 | **3,168** |
| benign g98 ss.99 ON | 3,719 | 3,086 | **3,022** |
| benign g05 ss.70 ON | 1,667 | 1,408 | **1,407** |
| benign g50 ss.99 OFF | 8,495 | 8,498 | 8,496 |
| junction g50 ss.99 ON | 2,228 | 2,064 | 2,091 |
| sparse g50 ss.99 ON | 3,261 | 3,664 | 3,582 (boundary 482 → 885 → 802) |
| sparse g98 ss.99 ON | 3,182 | 3,161 | **3,147** |

**The sparse residual, dissected to one face, and what it taught (2026-09-02).** All of it sat in
`capaltstart`, at the inside piece's face toward the intron (offset 9000 on + genes, 8000 mirrored):
crossing 302 fragments, truth 1.000 (pure gDNA), `transfer` 0.972, form B 0.748. Not the premise: an
UNDER-claim. The face's own numbers say why: the exon-to-face opportunity ratio is 1.30 for gDNA but
0.69 for RNA there (`eff_rna` is capture-aware, `eff_gdna_global` is not; on the benign panel the two read
1.27 and 1.22 and nearly cancel). The owner's derivation needs both components to see ONE opportunity
ratio between the two locales — that is what lets the enrichment ratio cancel — and the shipped face
map's opportunities violate it under hostile probes. Two forms were then measured:
* the pure COUNT form (`f_b = f_E (U+S)/U`, every opportunity dropped) is CATASTROPHIC (benign
  `g50 ss.99 ON` 14,576 vs 3,168; junction 18,597): the placement geometry of a crossing versus a
  contained fragment (~1.27 between a face and a region) and of a spliced fragment versus a crosser is
  load-bearing and must stay;
* the GEOMETRIC form — both components on the capture-blind geometric opportunity `a_g` at both locales,
  the spliced density as `S / A_g^b` — keeps the geometry and drops only the capture asymmetry.

**The final table (item 1, `formb_g` = geometric form + the marginal width; `rev_g` its reversal):**

| panel · row | silent | transfer | formb_m | **formb_g** | rev_g |
|---|---|---|---|---|---|
| benign g25 ss.99 ON | 2,403 | 2,216 | 2,171 | **2,189** | 12,557 |
| benign g50 ss.99 ON | 3,356 | 3,288 | 3,168 | **3,210** | 18,722 |
| benign g98 ss.99 ON | 3,763 | 3,719 | 3,022 | **3,084** | 5,000 |
| benign g05 ss.70 ON | 1,625 | 1,667 | 1,407 | **1,408** | 3,784 |
| benign g50 ss.99 OFF | 8,572 | 8,495 | 8,496 | 8,496 | 9,418 |
| sparse g50 ss.99 ON | 3,300 | 3,261 | 3,582 | **3,269** | 4,430 |
| sparse g98 ss.99 ON | 3,431 | 3,182 | 3,147 | **3,137** | 3,287 |
| sparse g50 ss.70 ON | 10,651 | 11,798 | 11,891 | **11,744** | 12,566 |
| junction g50 ss.99 ON | 2,486 | 2,228 | 2,091 | **2,092** | 19,617 |
| junction g98 ss.99 ON | 2,778 | 2,432 | 2,213 | **2,200** | 2,390 |

`formb_g` wins every stranded / part-stranded capture-ON row on all three probe worlds, is neutral off
capture, exact silence on unstranded and zero-control rows (the deadband gate), and its reversal fires
everywhere it acts. It gives back a little of `formb_m`'s benign gain — the owner's compromise currency,
bounded damage over the last few fragments. ⭐ The lesson for every later message: a composition
transfer must convert counts to densities with ONE opportunity treatment for all components; a
capture-aware opportunity on one component alone re-introduces a level across locales.

**The full benign sweep, `formb_g` vs the shipped `transfer` (30 conditions):** IDENTICAL on 19 rows (every
unstranded row, every `g00` row, every capture-OFF row except two that differ by ONE fragment), BETTER on
the 9 stranded / part-stranded capture-ON rows — `g05 ss.70` 0.845, `g05 ss.99` 0.996, `g25 ss.70` 0.965,
`g25 ss.99` 0.988, `g50 ss.70` 0.966, `g50 ss.99` 0.976, `g98 ss.70` 0.956, `g98 ss.99` **0.829** — and
`rev_g` fires on every one of them (up to 5.7×). The ladder (16 conditions, run as 8 parallel
per-condition processes — the serial harness took ~4 min per ladder condition) is the shipping judgement.

## Item 1 — THE LADDER (16 conditions, 8 parallel per-condition processes; 2026-09-02)

| unstranded (ss.50) | transfer | formb_g | rev_g |
|---|---|---|---|
| all 8 rows | — | **byte-identical** (the deadband gate: exact silence) | identical |

| stranded (ss.99) | silent | transfer | formb_g | rev_g | boundary tr → g |
|---|---|---|---|---|---|
| g00 OFF / ON | 30,606 / 20,787 | 18,960 / 18,197 | identical | identical | — |
| g05 OFF | 44,519 | 44,155 | 44,163 (+8) | 45,824 | 19,769 → 19,777 |
| g05 ON | 85,294 | 85,302 | **84,886 (0.995×)** | 111,221 | 58,444 → 58,027 |
| g50 OFF | 125,634 | 117,679 | 117,712 (+33) | 134,036 | 48,514 → 48,547 |
| g50 ON | 260,629 | 244,163 | **242,803 (0.994×)** | 349,664 | 168,737 → 167,377 |
| g98 OFF | 126,467 | 108,890 | 108,917 (+27) | 109,995 | 44,953 → 44,979 |
| g98 ON | 298,597 | 288,834 | **285,139 (0.987×)** | 288,342 | 209,269 → 205,574 |

**Verdict.** The message is exact silence on every unstranded row and every zero control; on the
stranded half it wins every capture-ON row (−416 / −1,360 / −3,695 fragments; 0.5–1.3 %) and costs
8–33 fragments (< 0.03 %) on the capture-OFF rows, where the intron row is already sharp and the exon
adds nothing the boundary lacks. The reversal fires on every row where the message acts (up to 1.43×).
Smaller than on the test chromosome (3–17 %) because the ladder's whole-library error is dominated by
regions, not boundaries, and its faces are fewer per fragment of error; the direction and the sign
agreement between the two substrates hold (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank` satisfied).
⭐ Meets the bar for a stranded-data message: derived, checked on certified truth, falsified twice, wins
the rows it exists for on three probe worlds and the ladder, harmless elsewhere. READY FOR PROMOTION —
the owner's ruling — as ONE mechanism: the exon's own strand row, gated by the strand deadband, read at
the geometric-opportunity splice-in map, widened by the marginal over the measured ratio; delivered as a
`PsiMessage.lam_rows` row at the intron|exon boundary beside rung 1's intron row.

## Item 1 — SHIPPED (2026-09-02, owner: "land this")

Landed in `messages/transfer.py` as ONE mechanism beside rungs 1–3: `splice_out_row` (the geometric
splice-in map read backwards, marginalised over the face's measured ratio on equal-probability nodes
with the trigamma counting price rung 2 charges), the exon block on a shared `licensed_faces` helper
(rung 2 rewritten on the same helper — one licence, one flux-column rule for both directions of a
face), the strand triple passed by `calibrate`; `simplex_logodds.strand_row_logodds` is the public
row a region publishes about itself (the solver's own frozen-variance term). Ruling: `DESIGN.md`
§6b.4. Gates: 3 new, each watched failing first; perturbations watched firing: a capture-aware
opportunity in the map (2 gates), the deadband dropped (1), the subtraction reversed (2). Faithfulness:
`policy_benchmark.py --policies silent transfer` reproduces the harness's `formb_g` TO THE FRAGMENT
(`g98 ss.99 ON` 3,082, `g25 ss.99 ON` 2,206) — ⚠ only once the harness's parent message was switched
off: after a landing a subclass that calls the parent's `prepare` DELIVERS THE LANDED MESSAGE TOO and
double-counts it (the harness read 3,032 / 2,228 until then; three innocent candidates — node set,
width form, reference clip — were chased first). The recorded caution "compare src-vs-src across a
landing" has this second face: a harness built on the parent class is dead the moment the parent
gains the mechanism.
