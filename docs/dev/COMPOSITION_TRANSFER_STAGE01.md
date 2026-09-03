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

## RUNG 4 — prototyped, measured, SET ASIDE (2026-09-02)

The three mechanisms above were prototyped through the seam and measured on the benign panel, both
adversarial panels and the ladder. Verdict (`MESSAGE_RUNGS.md`, "recorded and set aside"): the composed
join-only transport is NEUTRAL on the ladder (5/8 better in each half, worst 1.006×/1.008×, −6 % at the
in-scope `g00 ss.50 OFF` zero control); the inside bound is REFUTED for shipping (in-scope +1.3 %,
sparse probes +8 % — the uncertified enrichment sign at an exon|exon terminus, as derived); the flip
falsifications fire. The 1.05 M fragments of walled-exon mass did not turn into whole-library error
movement: the refit prior already serves those exons (a `nest` wall reads 0.666 vs truth 0.630 under the
shipped policy where silence reads 0.001). The pieces re-enter one at a time as items 5–9 earn them.

# THE RESET (owner, 2026-09-02, late)

Rungs 1 and 2 were unfinished (the exon-side messages had been nullified) and come first, one message
at a time; rung 4 grows one structure per step (`altstart` only); `MESSAGE_RUNGS.md` is the tracker.

## THE BASELINE ON THE TRIMMED SUBSTRATE (55 genes: twin + mono + `altstart`; 480 k fragments; 2026-09-02)

`policy_benchmark.py --panel test --policies silent transfer`, the shipped rungs 1–3 against silence —
the number every item in `MESSAGE_RUNGS.md` is measured from. Unstranded: transfer beats silence on
**14/20**, worst **1.16×** (`g00 ss.50 ON`, the deferred zero control, 64,143 vs 55,290); stranded
**8/10**, worst **1.00×**. In-scope unstranded OFF rows: `g05` 1.00×, `g25` **1.06×** (12,738 vs
12,042 — the one in-scope row where the shipped policy costs), `g50` 0.99×, `g98` 0.92×. Blind
capture-ON rows: `g05` 0.08×, `g25` 0.76×, `g50` 0.05×, `g98` 0.07×. Full table in the session's
`baseline_trimmed_silent_transfer.txt`; every earlier test-chromosome number is on another substrate.

## Item 1 — SHIPPED (2026-09-02); the record moved

The exon → intron|exon boundary message: derivation, the certified-truth check, the prototype's
tables, the transfer variance (counting + premise, the two-witness estimator, the marginal width, the
geometric-opportunity form), the three-panel and ladder tables and the landing all MOVED to their
homes — the ruling and measurements to `DESIGN.md` §6b.4, the premise-bias decision to
`ISSUES: splice-out-premise-bias-uncorrected`, the harness double-count lesson to
`TRAPS: a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism`, the gates to
`tests/calibration/test_transfer_policy.py`. The session harness (`item1_proto.py`, scratchpad) was
promoted as `scripts/design/policy_prototype.py`.

## Item 2 — SHIPPED (2026-09-02); the record moved

The intron|exon boundary → intron message: the derivation (one shared population, the s = 0 map the
identity, the boundary's own strand row verbatim, the derived deadband at the boundary), stage 0 on
certified truth (zero excess variance over counting on all four substrates; the plug-in Poisson null
that exposed the digamma estimator's vertex artifact; the ladder's mate-gap mature crossings at
0.3–0.7 % of boundaries), the information budget, the node-local and whole-library tables on the three
test panels and the ladder, the reversed-row falsification and the prior-mediated `g05 ss.99 ON`
residue all MOVED to `DESIGN.md` §6b.5; the prior's second exposure to
`ISSUES: gdna-landscape-trains-on-false-positives`; the gates to `tests/calibration/test_transfer_policy.py`.
Two session instruments (the node-local scorer, the pair-gap census with its null) are described in
`NEXT_SESSION.md` for promotion. The prototype (`item2_proto.py`, scratchpad) subclassed the shipped
policy with `strand` passed THROUGH so item 1 stayed on in both arms — one thing varied — and is
retired with the landing.

# ITEM 5 — SHIPPED (2026-09-02); the record moved

The exon|exon terminus boundary from the outside flank: the orientation table, the certified
directed licence on both substrates, the census of the ladder's terminus boundaries by what lies
beyond the outside flank (half are chains of termini a dozen bases apart), the REFUTATION of the
verbatim own-row exchange and its diagnosis by component (the mature depletion of the unspliced
crossing), the corrected licence with the spliced crossing, the three messages, the node-local and
whole-library tables on all four substrates, the structural no-echo law and the two stage-0 lessons
all MOVED to `DESIGN.md` §6b.6; the gates to `tests/calibration/test_transfer_policy.py`; the chain-
of-termini derivation owed to the tracker's item 6b. The prototypes (`item5_proto.py` and its
refuted verbatim form `item5_proto_v0_verbatim.py`, scratchpad) are retired with the landing.

# ITEM 6 — SHIPPED (2026-09-02); the record moved

The owner's abundance-discrepancy rule into the inside flank of a terminus: the certified bracket
check, the three marginal forms (uniform, profile, fitted step) and their agreement, the falling-totals
sign error every prototype shared and its corrected support, the isolation of the ingredients (the
forwarded arrivals — held for the scan with a per-hop premise; the reverse direction — held), the
self-fitted hop premise, and the node-local, three-panel and ladder tables all MOVED to `DESIGN.md`
§6b.7; the gates to `tests/calibration/test_transfer_policy.py`. The prototype (`item6_proto.py`,
scratchpad, eleven arms) is retired with the landing.

## Item 7 — SHIPPED (2026-09-02); the record moved

The alternative splice site: the two licences (C through §6b.6's map with `S_b`, E through §6b.4's map
with `S_b + F`), the stage-0 certification on both substrates and by flank length, the component
dissection that named the capture taper, THE HOP PREMISE (the fitted step with its error, the excess,
the per-pair discrepancy), the refused flux-factor forms and every measurement now live in
`DESIGN.md` §6b.8 and `ISSUES.md` CLOSED/REFUSED (`the-flux-factor-hop-premise`); the substrate's
`altss` block and the REPLICATION RULE in `TESTING.md` §0a's YAML header. Instruments left in the
session scratchpad: `item7_stage0.py` / `item7_stage0_len.py` (the licences on certified truth, by type
and by flank length), `item7_components.py` (the by-component dissection), `item7_nodelocal.py` (the
node-local scorer at alt-ss objects), `item7_proto.py` (every arm tried: own rows, width-only fit,
step forms, flux forms, joint fit, per-pair discrepancy).
