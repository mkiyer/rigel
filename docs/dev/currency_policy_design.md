# THE THIRD POLICY — design, as ruled by the owner 2026-08-19 (second pass)

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** Rulings will move to `DESIGN.md`
    when settled; derivations to `EQUATIONS.md`. ⛔ Delete when the policy ships.
    ⚠ The class is still spelled `CurrencyPolicy` — the owner ruled the name does not matter right now
    and the policy replaces `RelayPolicy` when finished.

## VOCABULARY (owner correction, this session)

⛔ **"level" is out.** The word is **ABUNDANCE** (counts/bp; *density* is the same quantity). The two
transfer strategies are the **ABUNDANCE strategy** and the **COMPOSITION strategy**. ⚠ `DESIGN.md`
§0c.0, `ROADMAP.md` §0 and `hop_currency.py` still say LEVEL — a scoped rename pass is owed once the
policy settles; new code in `messages/currency.py` says ABUNDANCE only.

## THE ARCHITECTURE (owner, verbatim rulings distilled)

**Sources.** Every slot may hold an INITIAL BELIEF at init: `g1_locked` slots (structurally pure gDNA,
abundance = its own measured rate); INTRONS via the density model against the intergenic background
(single-stranded introns deconvolve at init — rung-3 ruling); any single-stranded slot via the strand
model on stranded data. This is exactly `RegionInit` (`ctx.own`) — the policy invents no new source.

**The message** is the three abundances `{gDNA, RNA+, RNA−}` with per-component precisions.

**The sender sends blindly; the RECIPIENT decides** how to interpret (owner's semantics). Per hop the
destination asks ONE question of the static table: **is my transcript population the same as the
source's?**

**Population semantics** (from the owner's worked terminus example — TA+ (1000,10000), TB+
(5000,10000), TC− (500,20000), boundary B@5000 = TB's TSS):
* a REGION's population = the transcripts covering it; a BOUNDARY's population = the transcripts that
  CROSS it — a transcript STARTING at B is not in B's population;
* so `pop(B) = pop(left flank) ∩ pop(right flank)`, and a hop `R ↔ B` has equal populations **iff R's
  flank gains nothing at B** — one bit per (boundary, side), `terminus_flank_gain`'s own quantity,
  shared by both directions of the pair;
* in the example: R1(1000,5000) ↔ B: populations EQUAL ({TA+, TC−} both) ⇒ **COMPOSITION transfer is
  licensed even though B carries a TSS bit** — the raw terminus bit is NOT the licence, the SIDE is;
  B ↔ R2(5000,10000): R2 gains TB+ ⇒ ABUNDANCE strategy only.

**COMPOSITION strategy** (populations equal): rescale the message by the ratio of the two slots'
**measured total abundances — no belief enters the rescale** (owner ruling, answering the frame-debt
question). At an sj boundary the flux joins the flank-paired total and the SPLICE OUT / SPLICE IN
arithmetic applies (`EQUATIONS.md` §3.5e, worked numbers gated). ⭐ **the SAMPLING transfer variance — the
binomial up/down-sampling derivation** (owner: "the precision is still based on 2 fragments"):
rescaling transfers the source's SHARES, estimated from its raw counts; under a Dirichlet-Jeffreys
posterior over the K live components, `Var(log share_c) = trigamma(n_c + ½) − trigamma(n_tot + K/2)`
— every constant is Jeffreys, `n_c = 0` gives a finite wide claim, and multiplying by a 122× scale
factor adds NO precision.

**ABUNDANCE strategy** (populations differ): the abundances cross unscaled — they are claims about the
components the source could see; the destination's residual mass stays UNACCOUNTED (no rescale onto
gDNA — the mass pin's move is structurally absent). ⭐ **the DISAGREEMENT transfer variance**
(owner: "the greater the disagreement between total abundances, the less trustworthy the message"):
first derived form `σ² = (log(T_dst/T_src))²` — the squared mislift if the shared-rate premise is
wrong, the same form as `splice_in_frame_logvar`, no tuned constant. ⚠ The owner notes the honest
calibration eventually comes from a TOTAL-abundance landscape built BEFORE pass-0 (the Poisson-kernel
landscape is a candidate); that is a later rung, priced against this first form.

**Lifecycle** (owner's SPLICE-OUT walkthrough): a slot with NO belief RELAYS the (rescaled) running
message unmodified; a slot WITH a belief fuses precision-weighted (a strong local belief dominates —
"the message will likely dampen or die there"); at solve time each slot combines its two directional
messages precision-weighted, and a beliefless slot's combined message BECOMES its belief.

**Capture** (the 11.7× row's answer): the ABUNDANCE strategy's disagreement damping is what stops a
depleted rate from poisoning an enriched exon; the recovery of enriched exons is NOT the policy's job —
solve the solvable subset (single-stranded exons bordered by sj, via the boundary composition), fit the
BIMODAL gDNA landscape (Poisson kernel) on them, and the second pass pulls exons to the correct mode.
⭐ So pass-0 does not need to solve all exons — only enough to train the landscape.

## ⭐⭐⭐ THE KNOB — the two strategies are ONE mechanism (owner's intuition, 2026-08-20; built)

> *"I believe the absolute abundance transfer and the composition transfer are on two ends of a
> spectrum … there is actually a knob that connects them."*

There is, and it falls out of taking both claims seriously at once. ABUNDANCE says the enrichment is
1; COMPOSITION says it is exactly what the totals report. Both are point hypotheses about ONE unknown,
the log enrichment `eta`. Fusing them by inverse variance — the observed `log r` with its counting
variance `v`, against the no-enrichment premise whose error, if wrong, is `log r` itself — is a
SHRINKAGE estimator with no free constant:

    w = (log r)² / ((log r)² + v)      eta_hat = w · log r      Var(eta_hat) = w · v

so the delivered claim is `rho · r^w` — a continuous interpolation in log space, chosen PER HOP by the
data. `w → 0` at capture-OFF (the disagreement is counting noise) and `w → 1` under capture (it dwarfs
the noise). ⭐ And the residual `(1−w)·log r` is the ABUNDANCE strategy's transfer variance
`((1−w)·log r)²`, which vanishes at the composition end — ONE expression spans the continuum, so the
switch is gone rather than defaulted.

## ⭐⭐⭐ THE ENRICHMENT RATIO IS MODEL-FREE — the owner's fragment-length question, ANSWERED

> *"how are we trusting our enrichment ratios now? … It's a major problem when FL is highly
> different!"*

⛔ **A total abundance must never be `mass / effective_length`** — that divisor depends on the
composition being solved for (gDNA and RNA carry different length distributions), so the ratio is
circular and swings with the length gap. ⭐⭐ **The accumulator already deposits the exact quantity and
nothing was reading it**: the RECIPROCAL OPPORTUNITY per fragment — `1/(w−1)` crossing a BOUNDARY,
`1/(ell−w+1)` contained in a REGION — whose expectation is the density EXACTLY, for ANY length
distribution and ANY composition (`tests/native/_accumulator_reference.py`; `region_contained` and
`boundary_unspliced` carry it as `inv_length_sum`). It reached `CalibrationSubstrate` and stopped
there. Now `RegionGeometry.inv_abundance` → `StepContext.inv_abundance` → the knob's ratio.
⚠ So the equal-FL panel was masking a gap that is now closed on the enrichment path. What is NOT yet
model-free: the k-form's `M_dst` and the `E_c` divisors (composition-aware by construction, so exact
under the premise) — and the sj flux is not yet in the model-free total (the sj bank has its own
`inv_length_sum`; wiring it is the next small step).

## MEASURED (test chromosome, 8 scenarios, gDNA |err| in FRAGMENTS vs SilentPolicy)

| scenario | Silent | Relay | the knob |
|---|---|---|---|
| `g00 ss0.50 off` | 95,757 | 2 | **27** |
| `g00 ss0.99 off` | 24,836 | 2 | **41** |
| `g50 ss0.50 off` | 673 | 2,236 | **657** ✅ beats both |
| `g50 ss0.99 off` | 619 | 1,940 | **667** |
| ⛔ `g00 ss0.50 on` | 46,241 | 1,770 | 76,417 |
| ⛔ `g00 ss0.99 on` | 31,288 | 99 | 74,324 |
| `g50 ss0.50 on` | 65,707 | 36,158 | **45,408** |
| ⛔ `g50 ss0.99 on` | 5,632 | 19,086 | 25,942 |

⭐ **Capture-OFF is essentially solved** — the two zero controls fall from 26,663/26,760 (pure
abundance, no damping) to 27/41, a ~1000x reduction, and the two contaminated rows beat both shipped
policies. That is the disagreement variance doing exactly what the owner described.

⛔⛔ **THE OPEN DEFECT IS CAPTURE-ON, AND IT IS LOCALISED.** At `g00 capture_on` the worst slots read
`f_g = 1.000` where their own message-free solve reads 0.002–0.380 — the confident all-gDNA vertex.
Mechanism: under capture `k` is large, the rescale multiplies a SMALL gDNA claim by it, the delivered
share `rho_g·E_g/M` exceeds 1 and CLAMPS at the grid top, and ψ reads a confident "all gDNA".
⭐ The rescale is supposed to preserve the mass identity (shares summing to 1, no clamping possible),
so the identity is being broken between the rescale and the delivery — candidates, in order: the
two-direction combine averaging two claims each of which satisfies the identity alone; the
admissibility zeroing applied AFTER the rescale; and the sj flux missing from the model-free total.
⛔ Next session starts here, with `worst_objects.py` on `gdna_g00_ss_0.99_nrna_file_capture_on`.

## ⭐⭐⭐ WHERE THE ZERO-gDNA ERROR COMES FROM — MEASURED, and it is the REFERENCE'S MEAN

The owner's hypothesis, tested and confirmed on the local solve with messages OFF, scored against the
certified per-slot truth in FRAGMENTS. ⛔ The first run of this experiment read a NULL result because
`region_init` holds its OWN module reference to the solver, so patching the sweep's binding left the
local self-solve — the very thing under test — calling the original
(`TRAPS: an-ablation-that-never-ran`; the tell was pass-0 coming back byte-identical).

**① The Jeffreys ½ mean is the ENTIRE zero-gDNA local-solve error.** Setting the reference's mean per
object from the MEASURED background — `m_i = rho_bg·E_g,i / M_i`, with `rho_bg` the pooled gDNA density
of the structurally pure-gDNA slots, which are EMPTY at a zero-gDNA library — takes all four zero
controls to **exactly 0**:

| condition | shipped `m = ½` | anchor `m_i` |
|---|---|---|
| `g00 ss0.50 capture_off` | 95,757 | **0** |
| `g00 ss0.50 capture_on` | 46,241 | **0** |
| `g00 ss0.99 capture_off` | 24,836 | **0** |
| `g00 ss0.99 capture_on` | 31,288 | **0** |
| `g50 ss0.50 capture_off` | 673 | 711 |
| `g50 ss0.99 capture_off` | 619 | **597** |
| `g50 ss0.50 capture_on` | 65,707 | 65,917 |
| ⛔ `g50 ss0.99 capture_on` | 5,632 | 18,805 |

**⭐⭐⭐ AND THE 16-CONDITION LADDER AGREES, WITH A SHARPER SPLIT THAN THE TEST CHROMOSOME SHOWED**
(local solve, messages OFF, `nrna_mid` on every row so the nascent risk is genuinely exercised):

| stratum (0.8.0 scope) | rows | anchor `m_i` vs shipped ½ |
|---|---|---|
| both ZERO CONTROLS, capture-OFF **and** capture-ON | 4 | ⭐ **1,381,959 / 416,652 / 68,431 / 31,692 → EXACTLY 0** |
| unstranded × capture-OFF ⭐ IN SCOPE | 3 | **0.952 / 0.960 / 0.748** — every row better |
| stranded × capture-OFF ⭐ IN SCOPE | 3 | **0.956 / 0.987 / 0.756** — every row better |
| ⛔ stranded × capture-ON ⭐ IN SCOPE | 3 | **5.702 / 6.507 / 7.787** — every row worse |
| unstranded × capture-ON ⛔ DEFERRED | 3 | 1.000 / 0.999 / 1.441 — flat, then one loss at `g05` |

⭐ **Two of the three IN-SCOPE strata improve on EVERY row and both zero controls go to exactly zero;
the third regresses badly, and it is the one where the anchors are known to under-read.** ⛔ So this is
not a candidate to land as it stands — it is the measurement that says the reference's mean must be
per-object AND capture-aware, i.e. `ROADMAP.md` §1 rank 3, and that the two are one piece of work rather
than two.

⚠ **IT ALSO REFINES A RECORDED BLOCKER.** `ROADMAP.md` §1 rank 2 says *"`g05` REGRESSES 1.43x on both
strand settings"* — measured for a LIBRARY-WIDE mean. With a PER-OBJECT anchor mean, `g05` capture-OFF
**improves** (0.960 / 0.987) and only capture-ON regresses (1.441 / 6.507). The `g05` blocker is a
CAPTURE blocker, not a `g05` blocker.

⭐ This reproduces `object_composition.py`'s verdict (a per-object `m` is worth 25x, "exactly 0.000 at
both zero controls") **through the solver** rather than as a prior score, which is the step that was
missing. And it needs NO new constant: `_location_term`'s own derivation already says `location` IS the
object's prior expected composition and reduces to the shipped constant exactly at `m = ½`.

**② The one loss is the DOCUMENTED capture limitation, not a new defect.** Off-target anchors under-read
on-target gDNA by 2.6–3.6x under capture (`DESIGN.md` §0c's anchor ladder), so `m_i` is far too low at a
probed exon and the reference then under-calls gDNA where there IS gDNA. That is exactly the gap
`ROADMAP.md` §1 rank 3's spike-and-slab exists to close, reached from a new direction.

**③ THE DAMAGE IS AMPLIFIED BY THE REFIT LOOP, WHICH IS WHY THE STRENGTH KNOB LOOKS BETTER THAN IT IS.**
Weakening the reference's STRENGTH (`_JEFFREYS_REF` ½ → 0.05) makes pass-0 WORSE (99,115 → 107,991 at
`g00 ss0.50 capture_off`) and the FINAL answer perfect (94,366 → 1) — because a weak reference lets the
landscape learn the TRUE (zero) gDNA level, while the strong ½ teaches it "half gDNA" and the refits
compound it. ⛔ It is not a candidate: it costs `g50 ss0.99 capture_on` **3.6x**, where the library
genuinely IS ~half gDNA and the ½ mean is accidentally right. **The knob is the MEAN, not the strength.**

⚠ **What this does NOT resolve, and the owner named it first**: whether spliced fragments are enriched
depends on PROBE PLACEMENT — a probe spanning the junction enriches them, a probe deep in the exon
depletes them relative to the exon body — so no fixed rule for the flux frame can be right everywhere,
and it has to be learned. The policy's `(log r)²` frame damping is a "we do not know" treatment of
exactly that, not a solution to it.

## THE CONCEPT LADDER (one at a time, benchmark after each, both zero controls every time)

| concept | what lands | gate |
|---|---|---|
| **A** | the per-(boundary, side) population-equality table replaces the raw term-bit table | the owner's TA/TB/TC example as a unit test: R1→B licensed, B→R2 not |
| **B** | sources = `ctx.own`'s three abundances + precisions (g1 + intron factory + strand); running state relays unmodified through beliefless slots, fuses at believing ones; all three arms delivered to ψ | zero controls; the 8-scenario table |
| **C** | COMPOSITION transfer on equal hops: belief-free rescale (flank-paired totals incl. sj flux), §3.5e flux arithmetic wired, the Dirichlet-Jeffreys sampling precision | worked numbers end-to-end; benchmark |
| **D** | DISAGREEMENT transfer variance on unequal hops | the capture rows stop poisoning; benchmark |
| later | the landscape-calibrated disagreement variance; per-strand population equality; the λ/θ channels | — |

## Standing answers recorded (owner, this session)

* Rung 3: **yes** — introns deconvolve at INIT (density model vs intergenic background) and propagate
  that belief at pass-0; on stranded data every single-stranded node does.
* Naming: parked. The policy replaces `RelayPolicy` at the end of development.
* The test chromosome grows into a GOLD-STANDARD regression set — add a case whenever a new stressing
  situation is conceived. Queued candidates: a cassette (skipped) exon; a same-strand intron-retention
  isoform; divergent and convergent gene pairs; a probed exon shared by isoforms whose other exons are
  unprobed.
