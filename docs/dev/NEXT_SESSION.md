# NEXT SESSION — dissect the fan-out policy's residuals, stranded × capture-ON exons FIRST

    ⚠ **A DEV DOC, and it is a HANDOFF.** Read it, do the work, MOVE what settles into the permanent
    docs, and DELETE this file in the same session. (Precedent: five previous `NEXT_SESSION.md` files
    were each deleted per this instruction.)

## ① WHERE THE TREE IS

* Branch `fragment-length-gold-standard`, tree clean. **Suite: 0 failed / 3,719 passed / 9 xfail** —
  the baseline line in `CLAUDE.md` carries the full `+35` accounting for stages 3–4.
* **The first-pass redesign's stages 0–4 are ALL built, gated, priced and committed.** The working
  doc is `docs/dev/PLAN_first_pass_redesign.md`: §H the agreed design, §H.1 stage 3's verdicts,
  §H.2 stage 4's design (the completeness theorem, the three primitives, the owner's scope ruling),
  §H.3 the stage-4 verdicts, the residuals, and ⭐ **THE OWNER'S BINDING ANALYSIS FRAME** for this
  session's work — read that subsection before anything else.
* The shipped defaults are untouched: `message_policy` stays `"relay"`; the fan-out runs as
  `CalibrationConfig(message_policy="fanout")`. `measured_intron_reference` ships ON (stage 2).
* New since the last handoff: `rigel.calibration.structural_claims` (the stage-0 substrate, with the
  flank-completeness bits), `messages/fanout.FanOutPolicy` (the depth-2 pass-0 policy),
  `flank_to_exon_lambda` (the jointly-derived transfer), `PsiMessage.rna_one_sided` (per-message
  sidedness), and two instruments: `structural_claims_audit.py` (8/8) and `pass0_claimed_ab.py`
  (6/6) — the second is THIS session's main tool.

## ② THE WORK — the stranded × capture-ON claimed-exon regression, then the unstranded × OFF one

**The headline numbers (claimed slots only, `pass0_claimed_ab.py`, silent/relay/fanout):** the
fan-out wins everywhere it was designed to (deferred-stratum exons 169,659 vs relay 258,207 vs
silent 697,694; pooled ladder −37…−40 % vs the relay on all in-scope strata, 18/18 rows) — and
⛔ **stranded × capture-ON claimed exons regress 2.25× vs silence (40,528 vs 18,042 at
`g50 ss0.99 capture_on`)**, with a second, smaller residual at unstranded × capture-OFF exons
(+21 %, 7,755 vs 6,415). The owner: a resounding win overall, but the regressions get careful
attention, the stranded × capture-ON one first.

**THE OWNER'S FRAME (binding — the full text is `PLAN_first_pass_redesign.md` §H.3):**

1. **Split by gDNA RUNG** — `g05`/`g50`/`g98` × `ss_0.99` × capture-ON studied SEPARATELY, never
   pooled (`pass0_claimed_ab.py --condition …` per rung; add the DELIVER/REFUTE columns per rung).
2. **Survey the regressing REGIONS AND BOUNDARIES** — a slot-level dissection: which objects carry
   the regression's mass, and for each: the exon's OWN strand solve (mode + precision), the incoming
   message (mode + precision), the posterior, and certified truth. Bias and precision must be
   SEPARATELY visible per slot, because the owner's hypothesis is two-fold:
   * **① ACCURACY (bias)** — the imputed composition is good (it wins the unstranded scenarios) but
     less accurate than the strand model at stranded exons — anticipated. And the code is BRAND NEW:
     bugs may live in the implementation; the dissection doubles as debugging.
   * **② PRECISION (too high)** — the message's delivered precision evidently beats the strand
     model's. **The elegant fix is honest precision on both sides**: the strand model's exon
     precision honestly HIGH, the composed/transferred boundary precision honestly LOWER. The whole
     precision path needs evaluation theoretically AND at the regressed slots.
3. ⛔ **The owner holds a suspicion and will not act on it until the dissection data is on the
   table.** Do not land any mechanism first — the deliverable of step 1–2 is the DATA.

**The precision path to audit (the theory half of ②), file by file:**
* the exon's own side: `region_init.strand_evidence` → `tau_lam` (is the strand λ-term's precision
  at a deep stranded exon honestly high?) and `own_precision`'s count ceiling;
* the compose at the boundary: `messages/fanout._FanOutRelay.step` — intron `tau_lam` + the
  boundary's own `tau_lam`, precisions added (are those two honestly independent at a deep crossing
  count?);
* the transfer: `flank_to_exon_lambda`'s delta-method variance — Monte-Carlo-gated, but ONLY under
  the model's own assumptions: it prices sampling noise, NOT the premise ("the flank's composition
  is only locally-approximately the exon's"). ⚠ **"The drift" = that premise-variance term**, the
  measured excess when the TRUE value is transported (`hop_currency.py`: 5.1 % ladder / 10.2 %
  test-chr capture-ON, ~0 OFF) — one principled candidate inside branch ②, NOT presumed the fix;
* the one-sided branch's Jacobian `tau/f_g²` (`fanout.deliver`) — only where flanks are incomplete.

**Instruments for the slot survey:** `psi_channel_ablation.py` (every ψ-combine argument recorded,
then re-solved with one imputed channel nulled — exactly the message-vs-own question), `worst_objects.py`
(the concentration read: mechanism vs systematic), `toy_trace_error.py` on the test chromosome, and a
scratch per-slot table off `pass0_vs_oracle.measure_condition` (this session's `stage34_claimed`
scratch script pattern; consider promoting the per-slot version into `pass0_claimed_ab.py --dissect`
if it earns it). ⚠ Check each instrument's policy plumbing — some build their own
`CalibrationConfig`; the fan-out needs `message_policy="fanout"` reaching `calibrate`.

**Then the unstranded × capture-OFF exon cell** (+21 % on a 6.4k base): same discipline, own
dissection, no shared conclusions with the stranded one — different strata fail differently.

## ③ THE VALIDATION BAR for whatever fix eventually lands

One thing varied; per RUNG and per stratum, never pooled; claimed slots only
(`pass0_claimed_ab.py`), DELIVER/REFUTE split; both zero controls; test chromosome before the
ladder; falsification test first with the perturbation sweep (it has deleted dead clauses and found
fixture blind spots in every stage so far — trust it).

## ④ SESSION HANDOFF PROMPT — paste this to start the next session

```
Rigel — dissect the fan-out policy's residuals. Branch fragment-length-gold-standard, tree clean.

⛔ FIRST TWO COMMANDS:
     python scripts/design/preflight.py          (~2 s; --full is ~1 hour, don't)
     python -m pytest tests/ -q
   Baseline: 0 failed / 3,719 passed / 9 xfail. ANY failure is a regression.

READ, IN THIS ORDER:
  1. docs/dev/NEXT_SESSION.md   ← the handoff, written for this session. Delete it at the end per
     its own instruction, moving what settles into the permanent docs.
  2. docs/dev/PLAN_first_pass_redesign.md §H/§H.1/§H.2/§H.3 — the agreed design, the measured
     verdicts, and §H.3's OWNER'S ANALYSIS FRAME, which is binding for this session.
  3. CLAUDE.md's AXIOM 0 and the 0.8.0 SCOPE.

THE WORK — dissect and debug the STRANDED × CAPTURE-ON claimed-exon regression FIRST
(40,528 vs silent 18,042 at g50 ss0.99 capture_on, pass0_claimed_ab.py):
  - break the stratum by gDNA RUNG (g05/g50/g98 separately, never pooled), DELIVER/REFUTE split;
  - survey the regressing REGIONS AND BOUNDARIES slot by slot: own strand solve (mode, precision),
    incoming message (mode, precision), posterior, certified truth — bias and precision must be
    separately visible, because the hypothesis is two-fold (imputation bias — possibly a bug in the
    brand-new code — AND message precision too high against an honestly-precise strand model);
  - audit the message-precision path theoretically (the handoff lists it file by file);
  - ⛔ produce the DATA first: no mechanism lands before the dissection speaks — I have a suspicion
    and will rule on it when I see the survey.
Then the unstranded × capture-OFF exon cell (+21 %), its own dissection.
The owner drives commits.
```

## ⑤ OPEN ITEMS CARRIED, NOT THIS SESSION'S WORK

* **4g** — the landscape fit on SOLVED slots only (the estimate/bound completeness grades choose the
  training set) — after the residuals.
* The fl-gap side panels still carry the RETIRED uniform nascent model (re-simulation is an owner
  decision); `object_composition.PURE_GDNA_STRATA` still includes `R intron` (recorded, deliberate).
* The F-ledger's prospective deletions (`RelayPolicy`'s operators, the mass rescale, the global
  `ONE_SIDED_RNA`) stay gated on the full flow being priced — nothing is deleted on promise.
