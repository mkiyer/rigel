# The gDNA fragment-length model — diagnosis and root cause (2026-08-30)

Owner hypothesis: the gDNA fl estimate is contaminated by nascent RNA, the same structural error as the
strand overdispersion. ⭐ **CONFIRMED, quantified, and it is the dominant term off capture.** Two further
mechanisms dominate under capture. Instruments: `scratchpad/fl/pool_purity.py`, `shipped_vs_truth.py`
(origin-split oracle, no solver, no model).

## 1. The premise, and where it is written down

`fl.py`'s module docstring: **"Every pool is PURE BY CONSTRUCTION, and that purity is what removes the
circularity"**, and `scan_payload.py` labels pools 0 and 1 `# pure gDNA`. ⛔ The same file already
contradicts itself further down — `gdna_counts` is documented as "gDNA *plus* whatever nascent RNA sits in
an intron, **not a deconvolved pure-gDNA distribution**". The header asserts what the field docs deny.

## 2. MEASURED composition of the four "pure gDNA" pools (test chromosome, g05 ss0.99 capture-OFF, fl-gap arm)

| pool | gDNA | nascent | mature | RNA share | len gDNA | len RNA | pooled | bias |
|---|---|---|---|---|---|---|---|---|
| 0 intergenic contained | 7,804 | 0 | 8,881 | **53.2 %** | 259.4 | — | 181.2 | −78.2 |
| 1 intronic contained | 1,698 | 33,451 | 0 | **95.2 %** | 267.2 | 113.8 | 121.2 | **−146.0** |
| 2 intron\|exon crossing | 128 | 1,083 | 0 | **89.4 %** | 324.2 | 127.2 | 148.0 | −176.2 |
| 3 intergenic\|exon crossing | 65 | 0 | 0 | 0 % | 319.6 | — | 319.6 | +0.0 |
| **all four** | 9,695 | 43,415 | | **81.7 %** | **262.0** | | **140.9** | **−121.1** |

⛔ Pool 1 is **95 % nascent RNA**. ⛔ Pool 0 is **53 % MATURE** RNA — the shadow transcripts, unannotated,
so they read as intergenic: the contaminant is not only nascent, it is **any unannotated transcription**,
exactly as for the strand seeds. ⭐ The spliced pool is genuinely pure (gDNA share 0.00 %), so the asymmetry
is the same one the strand fit has: the RNA side is certified, the gDNA side is not.

## 3. ROOT CAUSE, and why no panel caught it

**The bias is the product of two factors: `bias ≈ RNA_share × (len_RNA − len_gDNA)`.** At g05 that is
0.817 × (113.3 − 262.0) = −121.5 against a measured −120.6. ⛔ **The ladder and the test chromosome give
gDNA and RNA EQUAL fragment lengths by design** (`DESIGN.md` §0b, the forcing function that stops the EM
splitting origins on length) — so the second factor is ZERO and a 95 %-contaminated pool reads a −0.8 bp
error. The defect is structurally invisible on every panel the tool is ranked on, and appears only on the
fl-gap SIDE arms, which are not a ladder rung.

Shipped model vs truth (mean L), capture-OFF, fl-gap arm — `ship−pool` is the shrinkage, `pool−true` the
contamination:

| condition | TRUE | POOLED | SHIPPED | ship−true | pool−true | ship−pool |
|---|---|---|---|---|---|---|
| g05 ss0.99 off | 261.5 | 140.9 | 140.3 | **−121.2** | −120.6 | −0.6 |
| g25 ss0.99 off | 261.1 | 200.1 | 199.0 | −62.0 | −60.9 | −1.1 |
| g50 ss0.99 off | 261.3 | 233.6 | 232.7 | −28.6 | −27.7 | −0.9 |

⭐ **The EB shrinkage is EXONERATED off capture (≤1.1 bp). The contamination is the whole error.**

## 4. Two further mechanisms, both dominant only UNDER capture

⚠ These are NOT the owner's hypothesis and must not be folded into it.

**(a) The opportunity divisor does not know about the probes.** On the EQUAL-LENGTH panel — where
contamination cannot move the mean — capture-ON still reads **+24.5 bp** (g05) and **+25.9 bp** (g25):
TRUE 233.6, POOLED 264.1, SHIPPED 258.1. Under capture the two CROSSING pools dominate, and a crossing
fragment must be long enough to span its boundary, so those pools are length-selected LONG. The divisor
(`gdna_opportunity_from_index`) is computed from the INDEX alone and is blind to the probe panel, so it
removes only ~6 bp of the ~30 bp selection.

**(b) The EB shrinkage anchor is the global MIXTURE, at a magic ESS.** `POOL_EB_PRIOR_ESS = 1000.0`
("Revisit the default on real-data pool sizes") shrinks the gDNA pmf toward `global_pmf` = every deposited
fragment — a mixture that is mostly RNA whenever gDNA is a minority. Inert while the pool is large; under
capture the gDNA pools collapse (n_gdna 53,110 → 5,106 at g05) and it bites: on the fl-gap arm at g05
capture-ON `ship−pool = −23.7` of a −31.7 total, i.e. **the shrinkage is then the DOMINANT term.** Same
shape as the deleted `Beta(14,14)`: a constant-weighted pull toward a contaminated anchor.

## 5. What a fix must satisfy (NOT a design — recorded so the next session starts here)

* ⭐ **The contaminant's length distribution is directly MEASURABLE**: nascent mean 113.8 against the pure
  spliced pool's 113.3 (and 214.7 vs 214.3 on the equal-length panel). So the pure RNA pool is a valid
  template for the contaminant and only the MIXING WEIGHT needs estimating — which the density
  deconvolution already estimates per object. An unmixing is therefore tractable; a purity assertion is not.
* ⛔ It may not assert any pool pure. Pool 3 (intergenic|exon) is the only one measured clean here, and
  that is a property of THIS annotation, not of the genome (`TRAPS: purity-is-a-property-of-the-annotation`).
* ⛔ It must be judged on the fl-gap arms, because the equal-length panels cannot see it — and those arms
  still carry the RETIRED uniform nascent model, so they need regenerating first (an owner decision).
* ⚠ The capture-ON mechanisms (4a, 4b) are separate defects with the opposite sign; fixing contamination
  alone will not close capture-ON, and pricing them together would confound them.

## 6. Three facts that raise the stakes (all verified at source, not inferred)

**(a) The contaminated pmf is not only an opportunity input — it is the EM's per-fragment length
likelihood.** `pipeline.py:817` turns `fl_models.gdna_pmf` into a `FragmentLengthModel` and hands it to the
scorer (`scoring.py` → the native LUT). ⚠ Not a contradiction of the retired length channel — `DESIGN.md`
retires it as a CALIBRATION COMPOSITION channel — but it means a −121 bp error in the gDNA length model is
applied per fragment in exactly the channel that separates origins.

**(b) A perfect gDNA length model is already PRICED at −5.90 %, and all of it is capture-ON.**
`length_ceiling.py`: the shipped solve's −6.31 % splits −5.90 % gDNA / −0.43 % RNA, "all of it capture-ON
(capture-OFF: −0.0000)". ⭐ That zero at capture-OFF is not evidence the model is right there — it is the
equal-length design making the contamination invisible, which is exactly §3. So the priced 5.90 % belongs to
mechanisms 4a/4b, and the contamination's own price has NEVER been measured on a panel that could see it.

**(c) There is already a recorded IMPOSSIBLE control failure that this diagnosis explains.**
`fl_anchor_gap.py`'s `G-gdna` is "⛔ THE CONTROL. `DNA_INTERGENIC` is pure gDNA and gDNA has NO INTRONS TO
MISS… If it MOVES, the fix reached fragments with no introns — which is impossible, and therefore a bug."
It holds at capture-OFF and reads **+6.0 % on all six capture-ON rows**, undiagnosed since 2026-08-17.
⭐ Mechanism 4a accounts for it: under capture the crossing pools dominate and their length selection is
corrected by a divisor computed from the INDEX, blind to the probe panel. ⚠ And the control's own premise —
"`DNA_INTERGENIC` is pure gDNA" — is the §2 assertion; measured here that pool is 53 % mature RNA once
unannotated transcription exists.
