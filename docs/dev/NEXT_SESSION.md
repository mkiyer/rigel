# NEXT SESSION — the state after the 2026-08-31 fragment-length deconvolution

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md` §1. MOVE anything that settles into the permanent docs and DELETE this file.

## What the last three sessions settled

1. **The gDNA strand overdispersion — SHIPPED.** Holds regardless of the annotation: the away-half
   moment, influence weighting, the two components reconciling. `EQUATIONS.md` §6a/§6b/§6c,
   `DESIGN.md` §3.3a. The premise it replaced is `TRAPS: purity-is-a-property-of-the-annotation`.
2. **The gDNA fragment-length model — PRICED**, and the pricing re-ranked it: not a calibration defect
   (a perfect pmf moves the 0.8.0 metric less than the benchmark resolves, in an inconsistent direction
   between the two fl-gap arms), but worth **8–32 % of the library gDNA bias** through the EM.
3. ⭐⭐⭐ **The gDNA fragment-length model — SHIPPED (2026-08-31).** `ROADMAP.md` §0 carries the state;
   `src/rigel/calibration/gdna_density.py` and `fl.py` carry the mechanism; the derivation and the
   refuted routes are `PLAN_gdna_fl_estimator.md` §3/§4/§7.1 and its §9.

## What the fl fix actually is, in three lines

The two CONTAINED pools are the same two shapes at different mixing weights, so a **contrast** cancels the
contaminant with **no template**: `g = [(1−a₁)f₀ − (1−a₀)f₁]/(a₀−a₁)`. The weights need one scalar, the
gDNA density, read off the **LOW side** of the per-object density where a contaminant that only ADDS
cannot reach — the root of a one-sided Poisson moment with a closed form (De Moivre). **No constant
anywhere**, and no accumulator or schema change.

## ✅ SHIPPED — the owner's call, 2026-08-31, with the equal-length cost accepted

A pre-registered criterion said "the estimator must be inert where nothing is wrong", and on the ladder's
equal-length capture-OFF rows it is not, quite: three degrade at **1.25-1.7x the reseed floor** (`g98`
inside it, and the `g00` zero control **exactly unchanged** because it declines). ⭐ That cost is the
expected variance and not a defect — with equal component lengths there is no length bias to remove, so
only the contrast's ~2x variance remains. **The owner accepted the trade**, on these grounds: the ladder
TOTAL favours it by ~9x (both capture-ON rows improve well outside the floor — 8.2 % of the standing bias
at `g50`, **51x the floor**), and wherever the two components' lengths actually differ, as on real cfRNA,
the gain is 20-250x larger.

## ⭐ Why it survives hybrid capture — SOLVED, and worth understanding before touching it

The shared-contaminant premise is FALSE under capture (TV between the two pools' contaminants is 0.95,
against 0.06-0.14 off capture), and the estimator is good there anyway. **Not luck, and no longer
unexplained**: under capture the intergenic pool is *depleted, not impure* (measured purity 0.48 -> 0.955
at `g05`), so `a_0 = rho*E_0/n_0` exceeds 1 and CLIPS. At `a_0 = 1` the contrast collapses algebraically
to `g = f_0` — the intronic pool's coefficient is exactly zero, so **the term the premise was needed for is
removed, not approximated**. Verified to 3e-17.

⛔ **So the safety rests on `a_0` clipping**, i.e. on intergenic really being near-pure under capture. A
probe panel that put RNA back into intergenic space would break it silently. That is the thing to watch,
and it is the reason the two refuted detectors are not missed: there is nothing to detect while the clip
holds.

⚠ Three smaller risks: the one-sided rate is biased **high** by 0.2-6.7 % and its scaling with the
contaminated fraction is unestablished; the decline rule does not fire at `g00` on the test chromosome
(harmless, measured — the derived repair is an SE carrying `rho`'s own relative error rather than dropping
it as common-mode); and **nothing has been run on real data** (`real-data-is-a-test-input`).

## ⭐⭐⭐ THE INTEGRATED fl ESTIMATOR — prototyped 2026-08-31, measured, NOT landed

The owner's design: capture is a SPECTRUM, so regions and boundaries must be ONE solve with no regime
switch. Prototyped in full (`integrated_fl_prototype.py.txt` beside this file is the exact code) and
measured on all 32 rows of four panels:

* **stage 1** — the shipped contained machinery, unchanged.
* **stage 2** — per-boundary composition, with the REGIONS calibrating the BOUNDARIES: probes bind
  indiscriminately, so at one boundary gDNA and nascent share the enrichment and it CANCELS from the
  composition: `R_b = (rho_r,adj/rho_off)·(mu_r−1)/(mu_g−1)`, `a_b = 1/(1+R_b)`, with `rho_r,adj` the
  adjacent region's own one-sided RNA excess. The crossing pair is then solved COMPLEMENT-FIRST
  (`r_hat` clipped, subtracted from the mixture) — algebraically the same 2×2, conditioned so the noise
  enters damped by `(1 − a_mix)`, small exactly when the separation is small.
* **stage 3** — the EXCESS-enrichment correction: the sampled union already carries every class's
  uniform part, so the unsampled on-target classes (contained-in-exon, spanning-exon) enter weighted by
  `(eps_e − 1)+` only — measured per exon from its own boundaries, signed means so noise cancels. At
  `eps = 1` the correction vanishes IDENTICALLY: no capture ⇒ byte-equal to the shipped estimator.

**Measured (pmf vs origin-split oracle): matches (±1–2 bp) or improves on 31/32 rows; every sign error
(+26 … +35 bp rows) eliminated; ladder capture-ON −24 → −3.1…−3.4 bp with TV 3× better; capture-OFF
inert everywhere (≤ +0.6 bp).** The one exception: gdna_long `g05` capture-ON −43.3 vs −38.4, TV better
by 29 % — a shared residual, not a regression.

⭐⭐⭐ **MEASURED END TO END, AND THE FIRST VERDICT WAS WRONG.** I scored it on `gdna_frac_est`, called it
worthless, and closed it. TRANSCRIPT accuracy — the fl model's actual consumer, the EM's per-fragment
assignment — says otherwise. ⛔ The mistake to learn from: `em_fl_ceiling.py` prints `gdna_frac_est` under
a "THE PRODUCT" banner, and that label is right for the CONTAINED contrast (a composition fix) and WRONG
for this one (a shape fix whose consumer is the EM).

**Capture-OFF — safe, and better.** Every test-chromosome row inside the reseed floor (`g00` EXACTLY
unchanged); the ladder's `g50` capture-OFF row **−15,544 fragments and −2,470 false-positive mass**, both
far outside it.

**Capture-ON — large wins, one large loss.** −7,669 (**−72 %**) on gdna_long `g50` unstranded with
false-positive mass 155 → 8; −6,970 on the ladder `g50`; −407, −172, −155 elsewhere. ⛔ **But the ladder's
`g05` capture-ON row is +150,809 WORSE, 86× the floor**, so the ladder net is +128k worse and **the
estimator does not meet the "matches or improves everywhere" bar.**

⭐⭐ **THAT ROW IS NOT THE ESTIMATOR'S FAULT, and the ceiling proves it**: fed the simulator's EXACT gDNA
pmf the same row costs **+186,333**. It REWARDS a wrong pmf — two errors partly cancelling, `ROADMAP.md`
§0's effective-length ruler finding in a second place. The integrated estimator lands at 81 % of the
perfect-pmf damage: it behaves like an accurate estimator on a row that punishes accuracy. ⚠ Counterweight,
stated honestly: on the row where perfect HELPS it captures only 13 % of the gain.

⭐⭐⭐ **THE ROBUSTNESS RESULT IS INDEPENDENT AND IS THE STRONGEST ARGUMENT.** Poisson-thinning the contained
fragments toward the infinite-capture limit — which is what stronger depletion IS — the SHIPPED estimator
**collapses at ~1,000 contained fragments to +31 bp and stays there SILENTLY** (its contrast declines and
it falls back to the raw four-pool sum), while the integrated one holds to **19 fragments** and then
declines honestly. **50× margin**, and a correctness property, not an accuracy trade. ⚠ Not immune at the
true limit: it reads `rho_off` and the adjacent RNA excess from the contained regions, needing a RATIO
where the shipped one needs a SHAPE — that buys the 50×, not immunity.

⭐⭐⭐ **THE +150k REGRESSION IS DIAGNOSED — a SIX-ARM channel split at ladder `g05` capture-ON, and the
answer is a FRAME error, not an accuracy one.** The gdna pmf feeds two consumers, and they want DIFFERENT
estimands:

| arm (true pmf injected into…) | transcript delta | fp_mass delta |
|---|---|---|
| the SCORER only (the EM's length term) | **−2,018** | −2,109 |
| GEOMETRY only (opportunity → prior) | **+188,208** | +55,141 |
| both | +186,355 | +54,733 |
| the RNA pmf, both routes | −18,228 | −8,123 |
| the fully-consistent TRUE PAIR, both routes | **+171,624** | +50,320 |

⛔ **Pair-consistency is REFUTED as the mechanism** — the true pair regresses too. ⭐ The mechanism: the
geometry/prior conversion assumes UNIFORM PLACEMENT (`EQUATIONS.md` §4.4 — the capture residual is a
placement model), so it wants the uniform-frame chemistry law (~216 on the ladder); the REALIZED pmf
(240, capture-tilted) answers the SCORER's question. Feeding the uniform-frame consumer the realized
quantity double-counts capture; the fitted pmf's bias had been accidentally keeping geometry in its own
frame. ⭐ Also real and separate: the true RNA pmf is worth **−18k** on this row — the sj-selection bias
of the RNA fit (−4.75 bp) has its own price, exactly the owner's original instinct. ⚠ The drain reads
`rna_pmf`/`global_pmf`, never `gdna_pmf`, so it is exonerated for g-arms but IS touched by r-arms.

⭐⭐⭐ **ROUTING IS VERIFIED — the regression vanishes.** Realized law to the SCORER, fitted (uniform-frame)
pair kept for GEOMETRY, measured on the ladder: `g05` capture-ON **+150,809 → +525** (inside the
attribution bar), `g50` capture-ON −6,004 kept, `g50` capture-OFF +723 (inside the bar). ⚠ Read every
delta against the ATTRIBUTION bar of a few thousand fragments (`ROADMAP.md` §0's attribution-floor row),
not the reseed floor alone.

⛔⛔ **AND THE HUNT'S DEEPEST FINDING IS GEOMETRY'S HYPERSENSITIVITY TO THE pmf TAIL, in both
directions**: +188k at capture-ON from the frame slip, and −15.5k at capture-OFF from a ~3 % tail change
(the full-integrated pmf in geometry) — the short-exon contained opportunity `E[(ell−L+1)+]` is on a
knife edge (−66 % at ell = 100), so the prior's unit conversion amplifies tail perturbations. That is the
effective-length-ruler territory §0 already carries, now with a mapped, reproducible 188k-fragment lever.
⚠ The −15.5k capture-OFF geometry win is REAL and attributable but knife-edge-shaped; do not chase it
without understanding the ruler first.

**THE LANDING PLAN (next session, gated as usual):** `FLModels` grows a NAMED second estimand — the
realized gDNA law (the integrated blend when its inputs are present, else the fitted pmf) — and
`pipeline.py` hands THAT to the scorer's `gdna_fl` while every geometry consumer keeps `gdna_pmf`
unchanged (byte-identical there). Also worth its own arm: the RNA fit's sj-selection bias is worth an
attributable **−18k** on the regression row (`r_both`), the owner's original instinct — a separate,
independently landable improvement.

⛔ **WHAT BLOCKS SHIPPING: the `g05` capture-ON compensating error**, which exists independently of this
work and is visible in the ceiling. Fix that and this lands; until then the ladder net is negative.

**Known residual, measured both ways**: the contained-in-exon class needs the within-exon hybridization
tilt `w(L)`; `g_B` as its proxy fixes big-exon substrates (test-chr ON −8 → +1) and breaks small-exon
ones (ladder ON +7 → +45), because small exons' crossers sample a different part of the binding curve.
`phi` is the conservative choice that meets the bar everywhere. The principled fix is measuring the
binding curve; a weighted-isotonic attempt is in the prototype's history and was refuted by tail noise.

**Also queued (owner-approved)**: a gene-edge-overlapping shadow transcript on `test_chr`, so pool 3's
contamination and the readthrough robustness are testable at all — requires re-simulation.

## The live thread

⭐⭐ **The message policy — the standing charter, still at rung 0.** `MessagePolicy` remains byte-identical
to `SilentPolicy` on all 30 test-chromosome conditions. The bar is unchanged: **win on unstranded, minimal
harm on stranded, never pooled.** Before building any mechanism, answer rank 1a's question with no solver:
**is the information a blind unstranded or AMBIG slot needs actually present in its neighbours?** The
anchored twin block was designed for exactly that.
⚠ Re-baseline first — the policy numbers in `ROADMAP.md` predate both the strand-estimator change and this
one.

## Cautions that will otherwise cost a day

* ⛔⛔ **COMPARE AGAINST WHAT SHIPS, NOT AGAINST YOUR PROTOTYPE'S BASELINE.** This session called capture-ON
  a blocker for two revisions of a plan on the strength of a prototype scored against the two CONTAINED
  pools; what ships is the FOUR-pool sum, and against that capture-ON is a large WIN. A wrong baseline
  produced a confident, wrong, blocking verdict that survived external review.
* ⛔ **The benchmark resolves a LARGE effect.** A 1e−5 nudge to a nuisance parameter moves 1/30 rows by
  more than 0.5 %. Do not believe a single-row policy difference below ~2 %.
* ⛔ **ALWAYS RUN BOTH fl-GAP SIGN ARMS AND THE EQUAL-LENGTH CONTROL.** The same change has been a win on
  `rna_long` and a loss on `gdna_long` at both the calibration and the transcript level.
* ⛔ **`np.diff(payload.region_bounds)` is WRONG** — bounds are concatenated per reference, so it
  manufactures a phantom region at every junction. Use `gdna_density.region_lengths_from_partition`. ⭐ No
  shipped code had this bug (audited); a prototype did, and it made an estimator look 6× worse than it is.
* ⛔ **`docs/dev/` is a sandbox and nothing may cite into it** — a permanent doc that does is a suite
  FAILURE (`tests/test_docs_boundary.py`), which caught exactly that in a `ROADMAP.md` edit.
