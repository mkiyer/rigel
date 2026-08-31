# ROADMAP — where the tool is, and what to do next

Reading order for a new session: **`SUCCESS.md`** (how performance is measured) → **this file** (current
state and priority) → **`TRAPS.md`** (mistakes not to repeat) → `DESIGN.md` / `EQUATIONS.md` as needed.

⭐ **The map of this file:** the NUMBERS POLICY and the 0.8.0 scope block are the FRAME · **§0** where the
tool is · **§1** the ranked list and its per-rank detail, which is the point of the file · **§2** the open
questions · **§3** the rules · **§4** and **§4.1** what has already been priced and refused — read before
proposing a mechanism.

⛔⛔ **THIS FILE IS NOT A CHANGELOG — IT IS A CLEAR SET OF DIRECTIONS FOR THE NEXT STEPS** (owner ruling,
2026-08-14, and it REPLACES the rule this paragraph used to state, that a measurement is never deleted
here but only stamped). **A COMPLETED ITEM IS DELETED AND THE FILE IS COMPACTED.** No ✅ LANDED block, no
record of what was implemented, no report on finished work — the investigation stays in git.
⛔ **Before deleting one, check whether its durable content already has a permanent home; if it does not,
MOVE it there first** — a lesson to `TRAPS.md` as a NAMED rule, a ruling to `DESIGN.md`, a derivation to
`EQUATIONS.md`, and a number about where the tool IS to §0 below. ⛔ **MOVE, never COPY**: copying is what
creates two homes that then diverge.
⚠ **What is NOT a completed item and stays however old it is:** §0, including the ✅ rows whose job is to
say *do not work here*; §2's open questions; and §4 with its graveyard §4.1, which say *do not rebuild
these* and therefore point forward.

⭐⭐⭐ **THE FRAME IS THE 0.8.0 RELEASE SCOPE, IMMEDIATELY BELOW. READ IT BEFORE ANY RANKING IN THIS FILE**
— three strata are the target, one is deferred, the metric is calibration against oracle calibration, and
the length channel is retired until after 0.8.0. ⭐ **§1 is the ranked list and is the centre of this
file**; every other section exists to make it readable or to stop it being re-litigated.

## ⛔⛔⛔ THE NUMBERS POLICY — owner ruling, 2026-08-22, and it governs every section below

**THIS FILE NAMES INSTRUMENTS, NOT MEASUREMENTS.** The owner's instruction: *"avoid citing precise
numbers as it becomes quite cumbersome to keep these updated and maintained when we change things and
predisposed to stale docs and misinterpretation of the current tool behavior."* A figure written here
decays the moment the tree or the panel moves, and a decayed figure does not read as stale — it reads as
current behaviour, which is worse than no figure at all.

⭐ **So a claim in §0 or §1 states a DIRECTION and a MAGNITUDE CLASS, and names the instrument that
re-derives it.** "Message propagation costs the three in-scope contaminated strata and wins the deferred
one — re-derive with `policy_benchmark.py`" is durable

⭐ **A number may stay only if it is one of three kinds:** ① an EXACT structural identity (a ledger that
closes, a partition that sums to zero difference, a truth that is exactly zero, an inflation of exactly
1.0000×) — these are claims about correctness and do not decay; ② a DESIGNED PARAMETER of the panel or
the model (four gDNA rungs, two strand settings, the sparse model's own range and on-fraction); ③ a
figure a RANKING turns on, and then it must say what would move it.

⚠ **§4 and its graveyard are the deliberate exception.** Those rows are historical REFUSALS, each stamped
to the substrate it was measured on, and the evidence is what stops a refused mechanism being rebuilt — a
graveyard row without its number is an invitation. They are not descriptions of current behaviour and do
not decay in the way §0's figures do.

⛔ **THE SUBSTRATE MOVED TWICE AND BOTH MOVES INVALIDATE OLDER FIGURES.** On 2026-08-13 the ladder was
rebuilt 36 → 16 conditions (`pilot`, `flgap_short`, `flgap_long` deleted; the retired rungs `g01`/`g10`/
`g25`/`g75`/`g90` cannot be re-run without a further rebuild). On **2026-08-22** both panels were
regenerated with SPARSE nascent RNA (`DESIGN.md` §0b). **Assume any quantitative claim in this file
predates the sparse rebuild unless it says otherwise, and get a current number by running the named
instrument.** Re-derive rather than trust: a number that has moved is a result, not a documentation bug.

---

## ⭐⭐⭐ THE 0.8.0 RELEASE SCOPE — OWNER RULING 2026-08-14, AND IT IS THE FRAME FOR EVERY SECTION BELOW

Current version is **0.7.1** (`pyproject.toml`). The target is **0.8.0**.

⭐⭐ **WHAT 0.8.0 IS: A CALIBRATION RELEASE.** The focus of development is CALIBRATION, and **the metric
is the CALIBRATION RESULT scored against ORACLE CALIBRATION — not the end-to-end transcript number.** The
transcript number stays a thermometer (`SUCCESS.md`'s word), and it is measured and reported; it is not
what ranks the work.

**THE FOUR STRATA, AND THREE OF THEM ARE THE TARGET:**

| stratum | 0.8.0 |
|---|---|
| unstranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-ON | ⭐ **IN SCOPE** |
| ⛔ unstranded × capture-ON | ⛔⛔ **DEFERRED — not a development target for 0.8.0** |

⛔⛔ **THE DEFERRED STRATUM IS DEFERRED, NOT DROPPED.** unstranded × capture-ON **REMAINS in every
benchmark and every measurement and must keep being reported**. It is not a development target until the
other three are fully optimised. If it improves as a side effect of other work, that is a free win.

⛔⛔ **AND IT IS WHERE THE ERROR IS, WHICH IS EXACTLY WHY A PANEL TOTAL CANNOT RANK 0.8.0 WORK**
(measured 2026-08-13/14 on the rebuilt 16-condition ladder): unstranded × capture-ON carries **64.5 % of
transcript error and 90 % of gene-level error**. Rank on the three in-scope strata, per stratum, and read
the deferred one as a reported control — `TRAPS: never-pool-the-strata` now has a second, sharper reason.

⛔⛔⛔ **THE LENGTH CHANNEL IS RETIRED UNTIL AFTER 0.8.0.** The fragment-length / "length likelihood"
channel **as a CALIBRATION composition channel** is DEFERRED POST-0.8.0. **0.8.0 ships without it.**
⛔ Do not propose it, do not list it as a candidate, do not rank it.
⚠ It does **not** exist in `src/`: it was A/B'd on 2026-08-10 and never shipped, so this is a
documentation ruling, not a code removal. ⚠ `length_likelihood` in `src/rigel/second_pass.py` is a
**DIFFERENT thing** — per-fragment second-pass assignment — and is NOT affected by any of this.
⭐ **The record has permanent homes and is not repeated here:** the mechanism and its numbers are
`TRAPS: a-linear-likelihood-emits-a-sign` and `TRAPS: amplitude-fades-influence-does-not`; why a
single-level panel read it as a win is `TRAPS: a-single-level-panel-cannot-see-a-constant` (`0.0324 →
0.5222` against a truth of 0.507 on the flgap pair, and **54–57 % gDNA at the `g00` control, in a library
containing none**); the derivations are `EQUATIONS.md` §3d–§3e, which is explicit that deleting them again
is the wrong move. ⭐ **What the measurement CLOSED** is the search for a θ-independent channel to give an
AMBIG slot its own composition evidence: there is no such channel. ⛔ **What it does NOT close** is
unstranded × capture-ON, which stays blind and is the DEFERRED stratum.

⭐⭐ **WHY THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS, IN ONE SENTENCE: the EM ALREADY USES THE
FL DISTRIBUTION**, so a large gDNA-vs-RNA length difference lets it assign fragments on **LENGTH ALONE,
BYPASSING calibration entirely and MASKING bugs in it**. Equal lengths **FORCE the calibration phase to be
exercised**, which is the only way a calibration release can be measured at all — calibration itself uses
**strand and density**, plus belief propagation across objects, which is ON (`RelayPolicy`). ⭐ `TESTING.md` §0
and `SUCCESS.md` carry it in full; `TRAPS: a-length-gap-bypasses-calibration` is the name.

⭐ **Scenarios must be CACHED so calibration can be re-run extremely efficiently** (owner). `panel.py
cache` builds both caches, `build_scan_cache.py` is the scan half, and the pattern to copy is
the pattern to copy is a key made of the scan manifest plus a content hash of every producing source file,
with the module under test deliberately OUTSIDE the key so its edit loop is one second.

## §0 THE STATE

⭐⭐⭐ **THIS SECTION SAYS WHERE THE TOOL IS, QUALITATIVELY, AND NAMES THE INSTRUMENT THAT RE-DERIVES EACH
CLAIM.** Owner ruling, 2026-08-22: *"avoid citing precise numbers as it becomes quite cumbersome to keep
these updated and maintained when we change things and predisposed to stale docs and misinterpretation of
the current tool behavior."* ⛔ **So a figure that decays does not belong here.** A number survives in this
section only when it is (i) an exact identity or structural constant — a ledger closing, a partition
summing to zero difference, a truth that is exactly zero or exactly one, (ii) a **designed** panel
parameter, or (iii) load-bearing for a RANKING and unlikely to move — and where one survives, the row says
what would move it. Everything else is a DIRECTION plus a magnitude class plus the command.

⛔⛔⛔ **EVERY QUANTITATIVE CLAIM BELOW WAS MEASURED BEFORE THE SPARSE-NASCENT PANEL REBUILD OF 2026-08-22
UNLESS ITS ROW SAYS OTHERWISE, AND IS DUE A RE-MEASUREMENT.** ⭐ **The way to get a current number is to
run the instrument the row names.** That caveat is stated once, here, and is not repeated per row. ⛔ **Read
the stratum, never the total** — every row naming unstranded × capture-ON describes the **DEFERRED**
stratum, reported and not ranked (`TRAPS: never-pool-the-strata`).

| | |
|---|---|
| ⭐⭐ **THE PANEL — REBUILT ON THE SPARSE NASCENT MODEL, 2026-08-22** (`panel.py status`, certified by `calibration_oracle.py`) | The ladder is still **16 conditions** by design — four gDNA rungs (`g00/g05/g50/g98`) × two strand settings × capture off/on, one fixed fragment total, gDNA and RNA at deliberately EQUAL fragment lengths (`TRAPS: a-length-gap-bypasses-calibration`) — and nascent is now **SPARSE**: a per-gene-SPAN on/off draw at `on_fraction 0.50` with a log-uniform ABSOLUTE level independent of the mature one, so nascent above mature happens (the ruling is `DESIGN.md` §0b). ⚠ The nascent FRAGMENT share is EMERGENT rather than configured — price it in advance with `expected_rna_weights`, never read it off the config — and `0.50` is a **DEVELOPMENT STRESS** level, not real data (realistic is `0.10`), so a number measured here says what is ROBUST and never what is likely. ⭐⭐ **The one result that matters for reading the rest of this section**: at `g50` capture-OFF the mass-weighted intron true `f_g` is unchanged from the retired uniform panel, while a large share of live intron slots now read true `f_g` **exactly 1.0000** where the uniform model gave almost none. The aggregate held and the DISTRIBUTION moved, so a number that moved is attributable to the distribution — which is the contrast the rebuild was for. ⚠ Realised lengths stay near-equal (a few bp, widening slightly under capture; nascent runs a little long because a nascent span is not transcript-truncated), so the forcing function still holds. ⛔⛔ **THE fl-GAP SIDE PANELS WERE NOT REGENERATED and still carry the RETIRED uniform nascent model.** Their data on disk matches their own configs, so each is internally accurate, but the ladder and the side panel now use DIFFERENT nascent models and no claim may be carried across them. Re-simulating them is the owner's call |
| ⛔⛔ **THE POLICY BENCHMARK — THE STANDING MEASUREMENT** (`policy_benchmark.py`, gDNA absolute error in FRAGMENTS, per condition, never pooled) | ⚠ The three-way table that stood here was measured with `CurrencyPolicy`, DELETED 2026-08-27 — re-derive with the instrument rather than quoting it. ⭐⭐ **The verdict that survives and is what to carry**: propagation is net-harmful wherever the local solve HAS evidence and its value is concentrated where that solve is BLIND — which is why the bar is *win on unstranded, minimal harm on stranded*, never a pooled total |
| ⭐⭐⭐ **THE gDNA LENGTH MODEL HAS TWO ESTIMANDS, ROUTED AND COUPLED — 2026-08-31** (`FLModels.gdna_pmf` vs `gdna_realized_pmf`; gates in `tests/calibration/test_fl_realized.py`; priced by `em_fl_ceiling.py` and `quant_accuracy.py`) | "the gDNA fragment-length distribution" was ONE name over TWO quantities. **`gdna_pmf` is the UNIFORM-FRAME law** — what the chemistry makes, and what the opportunity/prior mathematics assumes, since those formulas count start positions under uniform placement and are applied where `rho` is fitted. **`gdna_realized_pmf` is the LIBRARY CENSUS** — capture's selection included, which is what the EM's per-fragment scorer conditions on. ⛔⛔ **The cost of confusing them is measured: the realized law fed to GEOMETRY is +188,208 misassigned transcripts on one ladder row**, while fed to the SCORER it helps; the fitted pmf's "bias" had been accidentally holding geometry in its own frame. ⭐⭐ **AND THE TWO ARE COUPLED, NOT SWITCHED BETWEEN.** Capture is a spectrum and at either end one stratum has no data — at zero capture the boundary pools are empty, at very strong capture the contained pools are — so the split between the estimands is carried explicitly and weighted by how well it is RESOLVED (`lam = S/(S+N)`, signal against its own sampling variance, no constant). ⭐ **`N = 1/m_C + 1/m_B` diverges when EITHER mass collapses**, so `lam -> 0` and the two laws become the SAME ARRAY from both ends; `lam = 1` is the uncoupled behaviour. ⛔ That replaced THREE binary thresholds, and removing them was measured **net −3,467 transcripts BETTER**, with every capture-OFF row byte-identical. ⚠ The on-target excess carries its OWN weight, not `lam` — `lam` asks whether the strata's LAWS differ, the excess asks whether an exon's ENRICHMENT differs from 1, and gating one on the other suppressed a resolved correction (share read 0.000 at 50× enrichment). Direction overall: the +150,809 regression → **+293**, net **−7,817** transcripts across eleven scenarios spanning both fl-gap sign arms and the whole capture spectrum, `g00` at +1. ⚠ **The RNA pmf does NOT mirror this** — its three consumers split scorer +156 / drain −822 / **geometry −18,208**, so RNA needs no second estimand, just a bias fix at `effective_lengths_em` (§2) |
| ⭐⭐⭐ **THE gDNA FRAGMENT-LENGTH MODEL IS DECONVOLVED — the purity assertion is gone, 2026-08-31** (`fl.build_fl_models` with `region_lengths`/`region_types`; the rate is `gdna_density.one_sided_rate`; priced end to end by `em_fl_ceiling.py`) | `fl.py` asserted "every pool is PURE BY CONSTRUCTION"; measured, the intronic pool is **95 % nascent** and the intergenic **53 % mature**. It now solves a **two-pool contrast** — the two CONTAINED pools are the same two shapes at different mixing weights, so the contaminant cancels with **no template** — and takes the weights from a gDNA density read off the **LOW side** of the per-object density, where a contaminant that only ADDS cannot reach (the away-half again; `TRAPS: purity-is-a-property-of-the-annotation`). ⭐ **Structural facts that will not decay**: no constant is introduced anywhere — the density is the root of a one-sided Poisson moment with a closed form (De Moivre), and the decline test compares the purity separation against **its own** standard error; it is the identity when the two pools agree; and it needs **no accumulator or schema change**. ⭐ Direction on the fl-gap arms, both sign directions, stranded and unstranded: the length pmf's error falls by **~90 %** and the deliverable `gdna_frac_est` captures **92–96 % of the perfect-pmf ceiling**, largest on an UNSTRANDED row. ⚠ Read the ceiling itself from `em_fl_ceiling.py`, never from here. ⭐⭐ **UNDER CAPTURE IT DEGENERATES TO THE INTERGENIC POOL ALONE, WHICH IS WHY THE FALSE PREMISE THERE IS HARMLESS.** The shared-contaminant premise fails badly under capture (TV 0.95 against 0.06-0.14 off it), but the intergenic pool is **depleted, not impure** there (purity 0.48 -> 0.955 at `g05`), so `a_0` exceeds 1 and clips; at `a_0 = 1` the contrast collapses ALGEBRAICALLY to `g = f_0` and the intronic pool's coefficient is exactly zero — verified to **3e-17**. The term the premise was needed for is not approximated, it is removed. ⛔ **The safety therefore rests on `a_0` clipping**: a probe panel that put RNA back into intergenic space would break it silently. ⭐⭐ **ON THE LADDER — the ranking panel — THE NET IS FAVOURABLE BY ROUGHLY 9×, and the sign differs by capture.** Ablated across seven conditions: the two capture-ON rows IMPROVE well outside the reseed floor (8.2 % of the standing bias at `g50`, 51× the floor; 30.5 % at `g05`, 6× the floor), while three capture-OFF rows degrade marginally at **1.25–1.7× the floor** and `g98` stays inside it. ⭐ **The `g00` zero control is EXACTLY unchanged** — the estimator declines there, as it must. ⛔ **Why capture-OFF costs anything at all**: the ladder gives both components EQUAL fragment lengths by design, so there is no length bias to remove and only the contrast's ~2× variance remains. ⭐ The same ablation on the equal-length TEST CHROMOSOME **improves**, and the difference is diagnostic rather than noise — that substrate carries shadow transcripts, so it has real contamination even at equal lengths, and the ladder does not. **The estimator helps wherever contamination exists and costs a little variance where it does not.** ⚠ **The marginal capture-OFF degradation still trips the pre-registered "must be inert where nothing is wrong" criterion, and that remains an OWNER CALL** — the criterion was written to stop exactly this being waved through, even though the panel total now favours the change |
| ⭐⭐⭐ **THE gDNA STRAND OVERDISPERSION — RE-DERIVED 2026-08-30, and the one calibration parameter whose fit is now robust to the ANNOTATION** (`gdna_strand.fit_gdna_strand_from_substrate`; QC line, or `GdnaStrandModel`) | The fit is the AWAY-HALF moment over every genic count- and strand-observable object, influence-weighted, and RECONCILED against the RNA component (`EQUATIONS.md` §6a/§6b/§6c) — no component shrinks toward a constant, and the CEILING is the only asserted constant left in the strand module. ⭐ **Structural facts that will not decay**: it is unbiased under ANY same-strand contamination of the seeds; at ρ = 0 it reduces EXACTLY to the previous pair-count estimator, which is why neither simulated panel moved; and its one blind spot is unannotated ANTISENSE. ⛔ **On real cfRNA a value AT the ceiling is a CLAMP, not a fit** — `clamped_at_ceiling` says so — and `effective_seeds` (a threshold-free participation ratio) says how many seeds the answer actually rests on. Direction on real data: the weighting moved libraries OFF the ceiling and raised effective seeds by two to three orders of magnitude; the gDNA CALL moved only ~1–2 %, so the value of this work is the integrity of the estimate, not its magnitude. ⚠ **No simulated panel can score it** — at a true ρ of 0 no seed can dominate, so the panels show no harm and never a benefit; the evidence is the blank chromosome (`test_shadow.gtf`) and the real-library census |
| ⚠ **the SUBSTRATE warning that goes with any policy claim** | Every claim about a policy must name its SUBSTRATE: a small substrate's objects mostly lack own evidence, so a message is nearly free there, and a three-way ranking has inverted between a toy and the panel (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`, `TRAPS: panel-before-src`) |
| ⭐⭐⭐ **AN ENRICHMENT RATIO IS MODEL-FREE** (`object_composition.py`) | The accumulator's reciprocal-opportunity banks reach the policy as `RegionGeometry.inv_abundance` plus the sj bank's own column, so a total never comes from `mass / effective_length` — a divisor that is a function of the composition being solved for. ⛔ **The exact-density claim holds at a BOUNDARY (`ρ·P(w≥2) = ρ`) and is FALSE at a REGION, where it is `ρ·P(w≤ℓ)`** — an order of magnitude at a short exon — so every REGION↔BOUNDARY hop carries a truncation factor (`TRAPS: a-cancellation-is-conditional-on-its-support`), and the per-component divisors `E_g`/`E_r` stay model-dependent. ⛔⛔ **AUDITED 2026-08-22: the capture-enrichment ANCHOR RATIO read on TOTALS — which is what the model-free detector actually sees — is roughly double the published off-capture figure at `g05` and over-reads the ON level.** It survives as a **BOOLEAN** (off and on are still separated by nearly two orders of magnitude) and is **wrong as a calibrated LEVEL** |
| ⛔⛔ **THE ZERO-gDNA ERROR IS ψ's REFERENCE MEAN, AND THE NAIVE PER-OBJECT FIX COSTS STRANDED × CAPTURE-ON** (`ladder_arm_ab.py` scored by `solvability_audit.py`; the reference harness is `vertex_ceiling.py --arm ref=A,B`) | Setting the reference mean per object from the MEASURED background (`m_i = rho_bg·E_g,i/M_i`) takes the local solve's error at **both zero controls to EXACTLY 0** on the test chromosome and on the ladder and improves **every** capture-OFF row — and regresses **every** stranded × capture-ON row by a large factor, which is the known in-gene anchor under-read. ⛔ So the reference must be per-object **and capture-aware**, which is why ranks 2 and 3 are ONE piece of work. ⭐ Under a corrected reference every policy including `SilentPolicy` gets all four zero controls right, so the relay's control win was largely compensating for the reference |
| ⭐⭐⭐ **THE EXON REFERENCE IS SOLVED AT CAPTURE-OFF BY ONE PER-OBJECT FORMULA — AND ITS ANCHOR IS NOW CONTAMINATED** (same harness; the anchor pool is `object_composition.py`) | `m_i = ρ_g,i·E_g,i / M_i` with `ρ_g` pooled prior-free over the anchors mature RNA cannot occupy, replacing BOTH shipped branches: **exactly zero phantom gDNA at both zero controls** — the test all eleven graveyard mechanisms failed — and an improvement on every capture-OFF rung, largest where gDNA is scarce. ⛔ **Not shippable**: stranded × capture-ON still regresses by a factor of about two after splitting the pooled anchor into the two whose ratio is the enrichment, and the residual is the in-gene anchor under-reading ON-TARGET gDNA — an OPPORTUNITY-MODEL repair, not an estimator one. ⛔ Applying the formula only at exonic slots is WORSE than base, so the UNIFICATION is the win. ⛔⛔ **RE-AUDITED 2026-08-22: `object_composition.PURE_GDNA_STRATA` includes `R intron`, which was exactly pure only while the panel held no nascent RNA and is now inflated — worst where gDNA is scarce.** Intergenic-only is inflated **exactly 1.0000×**, and the SHIPPED `fit_intron_background` pools intergenic only, so this is the INSTRUMENT's pool rather than shipped behaviour — but every figure in this row was produced through it and must be re-derived before it is quoted |
| ⭐⭐⭐ **ψ's λ BRACKET WAS TOO NARROW TO EXPRESS ITS OWN PRIOR — built, gated, priced, ships OFF** (`ladder_arm_ab.py`, swept by `arm_sweep.py`) | At the shipped bracket the reachable `f` floor sits orders of magnitude above the median density the prior points at, and `DensityLandscape.required_logodds_window` derives the correct bracket with **no chosen constant** (predicted matches measured on every stratum). Nearly every in-scope condition improves, one dense capture-ON rung regresses marginally, and `g05` improves on **both** strand settings — where every library-wide reference mean instead regressed it. ⚠ UNPRICED: memory at genome scale (a small multiple on `sweep_n_grid`) and the end-to-end thermometer |
| ✅ **ψ's COMPOSITION CLOSES — STRUCTURALLY** (gated by `tests/calibration/test_composition_closes.py`; accuracy from `solvability_audit.py`) | The composition is the IMAGE of ψ's two parameters, so **100.00 % of published objects close** on both axes and in every class, against roughly three quarters before. Accuracy per stratum sits inside the panel noise floor, so the win is CORRECTNESS rather than accuracy. ⛔ Taking MEANS everywhere also closes and is REFUSED — it is worse on every stratum, badly so under capture |
| ⛔⛔ **THE REFERENCE LOCATION IS DELETED** (owner refutation, 2026-08-24; `DESIGN.md` §6b.1) | `structural_reference`, `measured_intron_reference`, `_location_term` and both builders are GONE — the reference is the symmetric Jeffreys measure and asserts nothing. Measured price on this metric, same-session before/after: stranded × capture-ON +4–8 %, `g98` ss.99 OFF region +42 %, `g00` capture-OFF controls IMPROVE, `g00` capture-ON +30–45 % (the measured form's zero-control win — the density channel is the recorded path to winning it back). ⛔ Every §0 number above predates the deletion |
| ✅ **STAGE A — THE ACCUMULATOR IS DONE, AND THAT IS A MEASUREMENT** (`calibration_oracle.py` for the ledger, `length_ceiling.py` and `fl_anchor_gap.py` for the length models) | Perfecting BOTH fragment-length models is now worth a couple of percent of the deliverable, down from more than a fifth. The fragment LEDGER closes **exactly** (`n_chimeric` is **0** since the compatibility-before-chimera ruling; `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`) and every deposited fragment places **exactly one unit** across the objects it crosses, with **0 unaccounted, both origins**. ⚠ The one unfixed length defect is the gDNA pmf **UNDER CAPTURE**, which reads long and worsens with the mean; since it is exact OFF capture it is a PLACEMENT problem (`EQUATIONS.md` §4.4) and never a `gdna_opportunity` edit |
| ✅ **oracle FIELD certification** (`calibration_oracle.py`, stamped per condition) | **16/16 stamped** on the rebuilt ladder. ⚠ **That is a STAMP COUNT, not sixteen verifications**: the uniformity gate is VACUOUS BY DESIGN on the eight capture-ON rows and on the zero-gDNA rows, so it really verifies **6 of 16**. Read the stamp and the vacuity flag together, never the stamp alone |
| ✅ **THE SIMULATOR SIMULATES THE TOOL'S OWN TRANSCRIPTOME** (`simulator_gates.py`, `verify_capture.py`) | `rigel sim` builds a rigel index and draws annotated transcripts plus the synthetic nascent entities in ONE multinomial with `prob ∝ abundance × capture-aware effective length`, and a probe maps by GENOMIC OVERLAP to gDNA and to every transcript it touches on either strand. After that change gDNA and nascent are enriched under one probe at **the same rate**, where before nascent was enriched far less than gDNA |
| ✅ **THE HOP CURRENCY MAP, PER COMPONENT** (`hop_currency.py`, all 16 conditions, every gate green) | The hop key is `object class × {sj, term, sj+term}`, **not** the seven object classes — a `B exon\|intron` is a splice site or a terminus with OPPOSITE currencies (`TRAPS: an-object-class-does-not-see-a-terminus`) — and the components disagree, so the map is per component. COMPOSITION into an exon from a splice site on all three; at a terminus into an exon gDNA takes a LEVEL and the RNA components a COMPOSITION; an intron into its own boundary sits at floor. ⛔ **NEITHER currency serves `R exon <- B gene edge[term]` under capture**, and a substantial minority of imputed mass is nine or more hops from any measured gDNA |
| ⭐ **calibration, the THREE IN-SCOPE strata** (`prior_vs_oracle.py`, `solvability_audit.py`) | The library `f_gdna` error is small in the low single-digit thousandths to low hundredths, and the PRIOR the EM reads is within a few percent of a perfect one. ⚠ Measured on the 36-condition ladder retired 2026-08-13; **rank 9 is the re-derivation**, and until it runs these are the weakest numbers in this section |
| ⛔ **calibration, unstranded × capture-ON — the DEFERRED stratum** (`solvability_audit.py`, exon detail from `object_composition.py`) | **BLIND, not noisy.** It reports a near-zero gDNA fraction while truth spans the full range, hands the EM a gDNA prior short by the great majority of the mass, and at the exons answers near zero at every rung — so it looks acceptable at low gDNA **by coincidence**. ⭐ The mechanism is exact, not statistical: at `κ = ½` the strand λ-term is **identically 0** and no other channel reaches an AMBIG slot; the search for a θ-independent channel that would is CLOSED (the retired length channel) |
| ⛔⛔ **THE EFFECTIVE-LENGTH RULER — THREE FINDINGS THAT MUST BE READ TOGETHER** (`calibration_vs_oracle.py`; end to end with `quant_accuracy.py --arm oracle_ruler`) | ① At the zero-gDNA control the shipped shrinkage contracts the large majority of transcripts where the correct factor is **exactly 1.000**, capture-OFF included, with `rho_ref` fabricated from false-positive gDNA. ⛔ It is a SYMPTOM: feed the **SHIPPED** function correct composition arrays and it returns the correct factor at `g00` and the true factor on the deferred stratum, so repair the COMPOSITION it reads, not the function. ② ⛔⛔⛔ **A PERFECT RULER IS A *LOSS* ON THE TWO CAPTURE-OFF IN-SCOPE STRATA** — priced on all 16, more than twice base transcript error on each, with false-positive mass several times worse; the win is confined to the zero control and the DEFERRED stratum, and the single-condition reading invited exactly the wrong conclusion. ③ The reason is that **BOTH rulers are wrong**: the no-enrichment null `U` — the oracle's own gDNA total laid down at exactly uniform density, which is the physically correct capture-OFF field and must read 1.000 — reads essentially 1.000 with no fitting, while the shipped and oracle factors both sit well below it and the ORACLE one sits further away. ⭐ So `oracle_ruler` is an A/B between two wrong rulers and the real ceiling is the `U` ruler, which rank 6 has not yet run. ⚠ The two aggregates disagree several-fold because the contraction falls hardest on LONG transcripts — rank on `Σeff/Σfl` and read `ruler_n_moved`, never the aggregate |
| ⛔ **the PER-TRANSCRIPT prior lane is BUILT and NEVER PASSED, and it is the largest in-scope lever measured** (`quant_accuracy.py`'s oracle arms) | `rna_prior_weight` is complete end to end with its arms and falsifications, and the production call site in `pipeline.py` omits it, so the shipped EM carries **zero** per-transcript information. ⭐ A perfect PER-TRANSCRIPT prior takes the three in-scope strata to well under half their gene-level error and cuts false-positive mass by two orders of magnitude, while a perfect PER-LOCUS one does far less and neither rescues the deferred stratum — so the in-scope prize is the per-transcript one. ⚠ Passing it needs a weighting function; two are built and refused (§4), and what they settled is that the target is TOTAL abundance and that **the support problem is the whole problem** — no expressed transcript is wrongly zeroed, the error is that thousands of silent ones are not. The next candidate must be a SPARSITY mechanism, not a better mean |
| ✅ **the prior ASSEMBLER and its POPULATION** (`prior_vs_oracle.py`, `mass_prior_ab.py`) | With perfect masses in, the assembler's own relative error is a few parts in a thousand, and with perfect per-component shares it is a few parts in ten thousand; the composition claim `a_g:a_r` is **exact** against the unspliced pool. The conserved-count rewrite and then the right yardstick took it down by roughly two orders of magnitude (`TRAPS: score-the-consumers-own-count`) |
| **message propagation, and its priced side effect** (`relay_pool_ab.py`; the golden suite) | ON since 2026-08-18 (`message_propagation = True`, policy `relay`); its standing is the three-policy row above, and the old muted-versus-on percentage table is 2026-08-10 history predating the conserved-count rewrite. ⚠ The recorded price is real and in scope: two zero-gDNA golden AMBIG scenarios regress by three orders of magnitude |
| ✅ **end to end — the LIBRARY figure and the POOL split** (`relay_pool_ab.py`, rendered by `benchmark_report.py`) | Mean library `\|err\|` is in the low thousandths on the three in-scope strata and near half on the deferred one; at pool level the gDNA total is within a few percent in-scope and short by the great majority on the deferred stratum, where even a perfect prior leaves a few percent. ⛔ Never quote the pooled figure — it is almost entirely the deferred stratum. ⛔ Read `gDNA (total)`, never `gdna_est`, which excludes intergenic fragments |
| ⛔ **end to end — TRANSCRIPT assignment** (`quant_accuracy.py`; direction from `rigel.sim.net_flow`) | Roughly a third of RNA fragments are misassigned, and a perfect prior removes about a third of that. ⛔ **That ceiling is ONE STRATUM**: the deferred one improves by a large factor while the three in-scope strata are **neutral-to-worse**, and the gene-level split is sharper still. So calibration and assignment are two problems in two files, and for 0.8.0 this table is the thermometer |
| ⛔ **the ATTRIBUTION FLOOR — reproducibility has TWO sources** (`arm_identity.py`; the floor from `quant_accuracy.py --arm base_reseed`) | `EMConfig.seed` defaults to `None` **and** the multi-threaded scan re-associates the float banks per worker merge, so pinning the seed alone is insufficient (`TRAPS: the-deliverable-is-not-reproducible-by-default`). The signed-off spread is at the float64 noise level on `posterior_mean` and **exactly 0 on every integer column**. ⛔ Consequently `arm_identity.py` cannot gate a `quant_accuracy` arm — `--arm noop` differs on hundreds of fields against a reseed floor of over a thousand — so **no `quant_accuracy` delta below a few thousand fragments is attributable**. Re-derive the floor with `base_reseed` in the same session as the effect |
| ⛔ **the NASCENT channel** (`transcript_truth.py`) | Tens of millions of fragments are parked on entities whose transcript-level truth is **exactly 0**, invisible in every transcript table; the 2026-08-12 zero-prior rule for synthetic entities cut in-scope nascent misassignment by roughly an order of magnitude at a strand-gated gDNA price. ⛔ The "accurate per-entity nascent prior" successor is SUPERSEDED by the NASCENT SCOPE RULING (`DESIGN.md` §0b): `alpha = 0` — absent until proven — is the right default, and per-entity nascent accuracy is designing around nascent |
| ✅ **ONE NUMERIC CONVENTION, and the one place it is not honoured** (`rename_identity.py --check` for bit-identity) | A COUNT is an integer and a FRACTION is float64; nothing decodes a bank, and float64 is five to six orders of magnitude more accurate than the fixed point it replaced (`TRAPS: integer-channels-reproduce`). ⚠ The exception is the **f32 strand term at `κ = ½`**, which manufactures a spurious tilt growing linearly in depth where `w_pos` is **provably exactly ½** — negligible at panel scale, entirely pre-existing, and falsified decisively by evaluating only that term in float64, which returns ½ to nine decimals at every depth. Rank 10 owns the one-dtype repair |
| ✅ **the SUITE is green with no skips** (`python -m pytest tests/ -q`) | Zero failures, zero skips, and the xfails are **not one kind of thing** — some are the recorded price of a config default, the rest are executable records of proven defects whose fixes are panel-negative alone. ⛔ The collected count and its file-by-file accounting live in `CLAUDE.md` and are re-derived there, never carried here |

⭐⭐ **THE ONE SENTENCE, IN THREE PARTS.** ① On the LIBRARY gDNA fraction the tool is accurate on all three
IN-SCOPE strata and structurally blind on the DEFERRED one; that blindness is deferred, and the search for
a θ-independent channel that would fix it is closed. ② On TRANSCRIPT-LEVEL assignment the same sentence is
FALSE — a sixth to a fifth of fragments are misassigned even where calibration is perfect, and far more at
the zero-gDNA control — so these are two problems in two files, and for 0.8.0 the second is the
thermometer. ③ Two open defects reach the 0.8.0 metric: the effective-length shrinkage, which no ceiling
has priced because it is built BEFORE the function every measurement arm patches, and the per-transcript
prior lane production never passes.

⛔ **Read `mwae_all` / `Σ|err|`, never `solv%` / `mwae` / `conf-wrong` / `calib`** — those four share a
denominator the solver moves by declining to answer — and quote the SHIPPED column, never pass-0, where a
win once read an order of magnitude larger than it was on the deliverable
(`TRAPS: the-intermediate-is-not-the-deliverable`).

## §1 ⭐⭐⭐ THE 0.8.0 RANKED LIST — WHAT HAS TO HAPPEN, IN ORDER

⛔⛔ **RANKED BY WHAT MOVES THE 0.8.0 METRIC — the calibration result against oracle calibration, on the
THREE IN-SCOPE STRATA, per stratum.** ⛔ Never by the panel total: it is dominated by the DEFERRED stratum
(the scope block above carries that measurement), so ranking on the total ranks the thing 0.8.0 does not
ship.

⭐⭐ **THIS LIST CARRIES NO MEASUREMENT TABLE, AND THAT IS DELIBERATE.** A rank says WHAT to do, WHY it is
there, and THE ONE GATE OR CONSTRAINT that decides it. Where a figure would go, the row names the
**INSTRUMENT that re-derives it** — because a precise number written here goes stale at the next panel
rebuild and is then read as current behaviour, which is `TRAPS: re-record-the-baseline` committed in
documentation. ⭐ A number describing where the tool IS lives in §0, and even there it is re-derived rather
than trusted. ⛔ What stays exact in this section is an IDENTITY (a config flag, a predicate, a designed
parameter, a structural zero) and anything a RANKING turns on.

⛔ **NOT ON THIS LIST, DELIBERATELY:** the length channel (retired until after 0.8.0) and anything whose
only target is unstranded × capture-ON. §4 is everything already priced and refused — read it before
proposing a mechanism.

| rank | the work | why it is there, and the one gate that matters |
|---|---|---|
| **1** | ⭐⭐⭐ **THE DISSECTION LOOP, ON THE WORST IN-SCOPE SCENARIO — the owner's protocol (2026-08-20), and it is the METHOD rather than one task.** In order: ① run the full panel, all three arms, every scenario reported; ② pick the worst scenario **inside the three in-scope strata** (never the deferred one, which is otherwise worst every time); ③ rank its REGIONs and BOUNDARYs by error MASS and trace a SAMPLE of them from initialisation through every belief change, message and refit against the certified `slot_truth`; ④ design the fix, gated first; ⑤ **add the offending transcripts to the TEST CHROMOSOME**; ⑥ re-run the whole panel; ⑦ repeat on the next largest source | ⭐⭐ **Step ⑤ is what makes the loop compound** — the test chromosome becomes the repository of the combinations the tool actually fails on. ⭐ Step ⑤ concretely: edit the hand-authored sources under `scripts/sim/test_reference/`, including `test_abundances.tsv`, which is the ONE source of that reference's abundances and of its nascent pattern (its config is `abundance.mode: file`, so the `nrna:` block is dead by design); then re-simulate, re-certify with `calibration_oracle.py`, and confirm with `preflight.py`. ⭐ Step ③'s instruments are `worst_objects.py` then `calibration_walk.py`, and `exon_solvability.py` says which class is safe to train on. ⛔ **The last single-mechanism guess is closed as a lead**: the precision-weighted premise variance (the unweighted fit floored at zero on the panel's sparse hops) is landed and gated, moved most ladder rows by up to a small multiple, made `currency` best at all four zero controls — and left `SilentPolicy` winning all three in-scope contaminated strata. The gap is not one mechanism |
| **1a** | ⭐⭐ **THE MEASUREMENT THAT SURVIVED THE 2026-08-27 TEAR-DOWN** | Message propagation is net-harmful on every in-scope contaminated stratum, and its value concentrates where the local solve is BLIND (the deferred stratum, both capture-ON zero controls). ⭐ So the decisive question is whether a policy that speaks **only** into destinations with no own evidence beats one that speaks everywhere — a measurement, not a switch. ⛔ Score it split by whether the destination had its own composition evidence (`TRAPS: an-imputation-must-cost-something-every-hop`), on `policy_benchmark.py`, and name the substrate on every claim |
| **1c** | ⭐⭐ **EXPAND THE gDNA SPECTRUM — deliberately, not by multiplying the cross.** The panel holds `g00/g05/g50/g98`; the owner wants the spectrum filled (one, five, ten, twenty-five percent up past ninety) and has warned against dozens of benchmarks | ⛔ Levels are informative where behaviour CHANGES, not where it is smooth — interpolation is free, so each new level must be justified by a measured transition and should cross a REDUCED set of the other axes until an interaction is shown. ⚠ The cost is real (simulate, scan cache, oracle cache and certification per condition); pay it when the current scenarios are exhausted, which is what rank 1's characterisation establishes. ⛔⛔ **AND THE PANEL ESTATE NOW CARRIES TWO NASCENT MODELS.** The ladder and the test chromosome are SPARSE (`DESIGN.md` §0b); the two fl-gap side panels on disk were NOT re-simulated and still carry the RETIRED uniform fragment-share model. Each is internally consistent with its own config, so both remain valid on their own terms — but **a claim that spans the ladder and a side panel is not comparable**, and re-simulating them is the owner's call, not a prerequisite anyone should assert |
| **1d′** | ⭐⭐⭐ **THE REFERENCE PRIOR IS REFUTED AT THE CONCEPT LEVEL — do not repair it (owner + external review, 2026-08-24).** **The root cause is a COORDINATE mismatch, not a bad constant**: on the λ = logit(f_g) axis the data's information is `I ∝ N_eff·disc·[f_g(1−f_g)]²`, which vanishes at the vertices and is identically zero at κ = ½, while the reference's location tilt holds fixed nats there — so the prior-to-data ratio diverges exactly where truth lives. ⭐⭐⭐ **And the information does NOT vanish in `f_g`-space**: for `p = κ + f_g(½−κ)`, `I_{f_g} = N(½−κ)²/(p(1−p))` is bounded away from zero, so `[f_g(1−f_g)]²` is PURELY the Jacobian `df_g/dλ`. ⛔ Measured: the tilt is overturned at N = 3 fragments where the strand channel is alive and **NEVER at any depth at κ = ½** (0.7471 from N = 10 to N = 10⁶); and 82–95 % of the scored pass-0 error at unstranded conditions sits on slots the tilt decides, so the benchmark grades the tilt rather than the algorithm | ⛔ **What NOT to do, each measured:** a different CONSTANT (0.75 is optimal at none of 16; `σ(L)` is 3.3×/10.8× worse than nothing at the zero controls; the required location changes SIGN across the panel, certified 0.115 → 0.999); a data-DERIVED location (circular — it fits a prior on the data being deconvolved, and its exon form was built, priced and DELETED); or re-weighting the tilt (weakening helps only the zero controls and is up to 2.8× worse everywhere the library has gDNA). ⭐ **The candidates worth prototyping are a REPARAMETERISATION away from λ** (judged on `L`-invariance and the Berger–Bernardo cancellation before any panel number — the file already grids the nuisance tilt in an arcsine coordinate) and **extending the density CHANNEL** (`density_deconv.density_lambda_factor`, which already casts the background as a likelihood rather than a prior and ships at ss-intron REGIONs). ⛔ An information-weighted tilt `w = 1 − exp(−γI)` is DECLINED: it introduces `γ`, and at a slot whose truth genuinely IS a vertex it makes the tilt vanish and falls back to ½ |
| **1e** | ⭐ **THE UNSTRANDED × CAPTURE-OFF EXON CELL — a REFIT-versus-MESSAGE arbitration question, and much smaller than it was recorded as.** The "+21 %" is ONE RUNG of four: the other three are wins, one of them by an order of magnitude, and over the three contaminated rungs the stratum is a small net WIN. ⛔ **It is NOT a reference problem** — a `solvable_exon` never received a reference location even before the 2026-08-24 deletion (the location was neutral at every `mrna_active` slot), and the location concept is now deleted outright. ⭐ **What actually solves those exons is the REFITTED gDNA prior**: at κ = ½ the strand λ-term is identically zero and the intron factory does not reach exons, so a claimed exon's own λ evidence is exactly 0, and turning the refits off collapses the silent arm by more than an order of magnitude while the fan-out barely moves. The message is the ACCURATE voice there — it is simply out-performed by the refit and then displaces it | ⛔ **So the question to answer is a DESIGN one, not a bug**: a message and a refitted prior are two imputations arriving at the same slot with **nothing arbitrating them**. Neither knows the other exists. ⭐ Rank 1d's dampener does NOT touch this cell (byte-identical, by the `v_own = ∞` derivation), so this is genuinely separate work and should not be bundled with it. ⚠ Re-derive with `pass0_claimed_ab.py` (the deleted `--dissect` survey; its verdicts live in `DESIGN.md` §6b.2 and git), whose STAGE reading prints `calib_refit_iters` 0 beside the shipped count precisely so this comparison is visible; the arbitration belongs with rank 3/4's landscape prior, where the refit's own status is decided |
| **2** | ⭐⭐ **THE LAWS A MESSAGE POLICY MUST OBEY — they are RULES WITH NAMES, not a design queue** | The invariants earned by measurement live where a rule lives: a refused claim loses its precision WITH its value (`TRAPS: zero-the-precision-with-the-value`); an imputation must cost something on every hop (`TRAPS: an-imputation-must-cost-something-every-hop`); a claim's mode must be delivered in psi's own coordinate (`TRAPS: off-grid-message-mode`); and a single-source transform may only REDUCE a precision — only independent evidence may add — which the foundation spec now ENFORCES at runtime rather than asking anyone to remember. ⛔ The prescriptive list that stood here was a design queue for policy code DELETED on 2026-08-27; what was tried, measured and REFUTED that session is in git history and the session record, so check before re-proposing a mechanism |
| **3** | ⭐⭐⭐ **THE MEASURED PRIOR — the pre-solve TOTAL abundance, and the road back into the calibration solve.** ⭐⭐ **The goal in one line: ψ needs a proper prior, the correct one is FITTED from composition-free observables BEFORE any solve runs, and a TOTAL abundance is that observable** (`DESIGN.md` §3.1a-i–ii, `EQUATIONS.md` §2.3b). ✅ **Rung ① MEASURE** is done — `total_abundance`, START/END over `ℓ` under a derived wall rule, validated by the FIELD-FREE arm of `total_abundance_audit.py` on every condition, both zero controls included. ✅ **Rung ② INTEGRATE** is done — `fit_intron_background` behind `CalibrationConfig.background_abundance`, a tie off capture and a clear win under it. ✅ **Rung ③ FIT** is done — the `AbundanceLandscape`, censused by `abundance_landscape_census.py` and A/B'd by `landscape_head_to_head.py`, which is also what retired the NPMLE | ⭐ **Rung ④ consumes `rho_0`, the per-class enrichment responsibility `w`, and the ENRICHMENT DETECTOR — and NOT `span_R`.** `span_R` is grid-fragile: it moves by more than an order of magnitude with grid resolution because `split_basins` picks the enriched mode by basin MASS and over-resolution fragments the bulk into competitors, while `rho_0` and the anchor verdict are grid-robust (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`). Capture-OFF totals are bimodal from EXPRESSION, so the landscape alone cannot separate expression from enrichment — the detector is the other half. ⭐ It lands behind its own `CalibrationConfig` flag (`composition_reference`, default bit-identical), and the door is ψ's location term widening from one scalar to `(m_lo, m_hi, w)`. ⭐⭐⭐ **RUNG ④ NOW HAS A MEASURED REQUIREMENTS SPEC, and the derived REFERENCE is rung ④ rather than a separate track — but read rank 1d′ FIRST, which refutes the location-tilt concept outright.** Six requirements, each a measurement: ① the location's RANGE must span BOTH lattice ends (today it is `{½, 0.75}`, while the panel wants `σ(−L)` at `g00` — which beats asserting nothing by ~10⁴× — and `σ(+L)` at `g98`); ② `rho_bg` must be ON-RATE for the destination's enrichment class, because the intergenic pool carries EXACTLY ZERO probe bases and is therefore the UNPROBED rate — the single named cause of every capture-ON failure, and the reason behind the 2026-08-26 regions-only refusal; ③ exact at `g00` (the transported form structurally is; keep it, it is a free falsification); ④ strength stays one pseudo-fragment; ⑤ reference and message INTERACT and must be priced jointly; ⑥ both zero controls on every arm, the pure-gDNA/RNA-bearing split never pooled. ⛔⛔ **④a — widening the SHIPPED estimator's selector to `solvable_exon` — was BUILT, gated, priced on all 16 scenarios and then DELETED (2026-08-24), because rank 1d′ settles that a data-derived location is the circularity the whole concept is refuted for. Do not resurrect it as a flag.** ⭐ Its measurement stands as evidence about the PROBLEM: it wins all 8 capture-OFF rungs and takes all four `g00` rows to EXACTLY ZERO claimed-exon fragments — the `rho_bg = 0` branch, which does not care about capture — and it LOSES all 6 contaminated capture-ON rungs.** ⭐⭐ The loss is concentrated where the exon cannot argue back: catastrophic UNSTRANDED × capture-ON (`τ_λ ≡ 0` at κ = ½ leaves no own evidence to overturn a reference asserting near-zero gDNA), modest stranded. ⭐ **So a data-derived location wins exactly where the background rate is ON-RATE for the slot and fails catastrophically where it is not** — which is the cleanest statement of why the concept fails, and why rank 1d′ points at a reparameterisation rather than a better estimator. ⛔ **Rung ⑤ prices it per stratum, on both zero controls, against a SHUFFLE control (`TRAPS: attribution-must-survive-a-shuffle`), split by whether the destination had its own composition evidence.** ⛔⛔ **TWO CONSTRAINTS FROM THE 2026-08-22 AUDIT, BOTH ABOUT WHAT MAY BE BELIEVED HERE.** ⓐ **The intron-inclusive anchor pool is CONTAMINATED on this panel**: `object_composition.PURE_GDNA_STRATA` includes `R intron`, which was exactly pure only while the panel held no nascent, and nascent now inflates that pool UPWARD — several-fold where gDNA is scarce, negligibly where gDNA dominates — while the intergenic-only pool is inflated by nothing at all. ⭐ The SHIPPED `fit_intron_background` pools intergenic only, so shipped behaviour is clean and this is an INSTRUMENT defect; ⛔ but a measurement taken through that pool is not trustworthy, which retires the earlier reading that an intron-inclusive pool was licensed, and every §0 figure produced through it must be re-derived before it is quoted. ⓑ **The enrichment detector survives as a BOOLEAN and NOT as a calibrated LEVEL**: read on TOTALS, which is what the model-free detector actually sees, the off-capture ratio is well above the ~1 the older figure claimed and the ON level is over-read, because a ratio of totals is dominated by whichever component the claim is not about (`TRAPS: a-total-density-ratio`). The off-versus-on separation is still orders of magnitude, so it may decide capture-on-or-off and may never set a level; re-derive with `object_composition.py` before any level use. ⛔ Every cheap alternative is measured-refused: the library-wide truth mean, and BOTH single-mode per-object locations (pooled and flank-cancellation local), each of which damaged stranded × capture-ON and `B exon\|exon` identically |
| **4** | ⭐⭐ **CONFIRM THE LANDSCAPE THEN TRAINS ON A REAL SUBSTRATE** — nothing to build; this is rank 3's payoff. With exons BRACKETED rather than pinned, `_fit_gdna_hyperprior` trains on objects whose `f_g` reflects data | ⚠ The circularity to watch: it trains on `belief.f_g·mass` over expressed REGIONs INCLUDING exons, so too strong a reference means it learns the prior back. ⭐ Most of the training mass on the dense capture-ON rows already has its own composition evidence — measure that no-evidence share before and after with `composition_evidence_census.py`, and do not accept a change that moves only it |
| **5** | ⭐⭐⭐ **PERFORMANCE — MANDATORY BEFORE 0.8.0** (owner, 2026-08-17): the grid solve in memory-bounded parallel chunks, with an advanced CLI flag spanning one object at a time through many in parallel. Compute-versus-memory is the dial | ⛔⛔ **Optimise on HIGH-DEPTH REAL RNA-SEQ, not cfRNA** — this reverses the standing "profile on real cfRNA" instruction, which still appears in several places. ⚠ Calibration is depth-INDEPENDENT and scales with the INDEX while the EM scales with the DATA, which is why the panel profile and the genome-scale one are inverted and both correct (`TRAPS: toys-rank-hotspots-backwards`) |
| **6** | ⭐⭐⭐ **BUILD THE `U` RULER ARM** — the uniform-gDNA null ruler, the oracle's own gDNA total laid down at exactly uniform density, priced end to end | ⛔ **The oracle ruler is NOT the ceiling**: a perfect-composition ruler is a LOSS by a factor of a couple on the two capture-OFF in-scope strata, because there the correct factor is exactly 1.000 and the oracle ruler sits FURTHER from it than the shipped one — so `oracle_ruler` is an A/B between two wrong rulers. ⭐ The null reads essentially 1.000 with no fitting, and `calibration_vs_oracle.py` already carries the `U` column. ⚠ Admissible at capture-OFF only, and read `ruler_n_moved`, never the aggregate factor — the aggregate barely moves while most transcripts are redistributed |
| **7** | ⭐⭐⭐ **THE SHRINKAGE AT `g00`, REPAIRED UPSTREAM** — the shipped effective-length shrinkage contracts the great majority of transcripts where the correct factor is exactly 1.000, capture-OFF included, with its reference density fabricated from false-positive gDNA | ⛔ It is a SYMPTOM: feed the SHIPPED function correct composition arrays and it returns the correct factor at the zero control and on the deferred stratum alike. **Repair the composition, not the function** — and `priors.py` imports `_global_reference_density` from `capture_eff_length.py`, so one repair serves both consumers. ⭐ `calibration_vs_oracle.py` is the only instrument whose patch point is upstream of the ruler, which is why nothing else has ever priced this |
| **8** | ⭐⭐⭐ **PASS THE PER-TRANSCRIPT PRIOR LANE — and the weighting function IS the work.** `rna_prior_weight` is built end to end and the production call site in `pipeline.py` omits it, so the shipped EM carries zero per-transcript information | ⭐ The largest in-scope lever measured: a perfect per-transcript prior roughly halves in-scope gene-level error and cuts false-positive mass by orders of magnitude, and it does NOT rescue the deferred stratum — so the prize here is in scope. ⛔ Two candidate weighting functions are built and refused (`TRAPS: an-upper-bound-is-not-an-estimate`); what they settled is that the target is TOTAL abundance and that **the support problem is the whole problem** — no expressed transcript is wrongly zeroed, and the error is that thousands of silent ones are not. The next candidate must be a SPARSITY mechanism, not a better mean, and it is scored by `quant_accuracy.py`'s existing arms rather than by a second scorer |
| **9** | ⭐⭐ **RE-DERIVE THE IN-SCOPE NUMBERS ON THE CURRENT PANEL — and it now has TWO reasons, not one.** Most of §0 was measured on the 36-condition ladder retired 2026-08-13, and everything measured since predates the 2026-08-22 sparse-nascent rebuild | ⚠ Re-record the noise floor in the SAME session as the arm and quote it beside the effect: an arm inside the floor moved nothing, and an inherited floor is `TRAPS: re-record-the-baseline`. ⭐ This rank is also what makes the rest of the list quotable, so it is cheap insurance rather than bookkeeping |
| **10** | ⭐⭐ **THE f32 CUBE MANUFACTURES A STRAND TILT AT κ = ½ — one dtype, and the falsification is already written.** At κ = ½ the strand mean is ½ identically, so the tilt is provably unidentified and `w_pos` must be exactly ½; the AMBIG cube evaluates that sum in float32 and departs τ-dependently, giving a tilt that grows LINEARLY in depth | ⭐ Falsified decisively: the same cube with only the strand term in float64 returns `w_pos → 0.500000000` at every depth. ⚠ **Negligible at panel scale and entirely pre-existing**, so it would bite only a very deep library with a dominant locus. ⛔ Repair the strand term inside `_solve_ambig_logodds`, not the cube (its f32 storage is an authorised memory choice), and price it on the panel before landing (`TRAPS: panel-before-src`) |
| **11** | ⭐ **PRICE THE CANCELLING PAIR TOGETHER** — `struct_lock` rescoped to `g1_locked ∧ REGION` AND the `intergenic\|exon` boundary claiming its RNA-contaminated crossing mass as gDNA. Brief in §1.1 | ⛔ Neither half has an honest price alone: the `struct_lock` half by itself is catastrophic on the zero-gDNA control, by orders of magnitude (`TRAPS: a-cancelling-defect-pair`). Five xfails go green if and only if the pair lands. ⛔⛔ **It can only be priced with `--messages on`** — with them off it rewrites tens of thousands of slots and changes no scored field, which is `TRAPS: an-ablation-that-never-ran`. ⛔ **RE-PRICED 2026-08-26 WITH a local level channel in place (`ladder_arm_ab --arm stage1_pair{,_onesided}`, the measured intron reference as the replacement load) — STILL REFUSED**: marginal to the reference alone it worsens two of the three in-scope strata and the deferred one, and its wins are confined to `g00`, which rank 2 rules one-sided. ⭐ New fact for the brief: certified truth (`structural_claims_audit.py`) shows gene-edge BOUNDARIES carry EXACTLY ZERO RNA fragments on 32/32 conditions, so the mis-scoped mask's load is about zero-count slots emitting certainty, not contaminated counts |
| **12** | the index's duplicate map (an ALIAS MAP `dropped_t_id → kept_t_id`, not re-admission) · `mass_*_boundary` → `count_*_boundary` (they are crossing INCIDENCES) · restore the moment tests that went with the deleted length channel · the cheap ledger of dead surface and stale sibling references | correctness and hygiene, each its own commit (`TRAPS: one-thing-varied`), none of it moving the 0.8.0 metric. ⚠ The duplicate map needs an index rebuild but not a panel re-scan — verify with `rescan_panels.py` rather than assuming, since `reach` is covered by no other hash |
| **13** | ⏳ ⭐ **FINISH THE ORACLE-EFFECTIVE-LENGTH DIAGNOSTIC** — started, not finished; about half an hour, and it re-ranks the two ruler items above | ⛔ It needs a stashed pre-closure arm or it measures nothing — one arm is not an A/B. ⚠ Demoted because it ranks POST-calibration items while calibration's own path is ranks 1–4 |

⚠ **Tracked, not blocking:** two exact per-fragment mirrors of a PURE-RNA library deconvolve differently
in `mass_gdna_region`, by a few percent, and the effect is neither boundary-only nor monotone in
strandedness. That is a few percent of the zero-gDNA false-positive channel, and an R1-sense library is
now simulable, so it is measurable whenever someone wants it.

### §1.1 THE BRIEF FOR **the-cancelling-pair** — rank 11, kept whole because the analysis is complete

**THE TWO HALVES, AND WHY NEITHER HAS AN HONEST PRICE ALONE.** The certified-RNA channel is a LOWER BOUND
(`ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)` — the exon may hold molecules that never touch that boundary), and ψ
delivers it as a two-sided Gaussian, penalising a destination holding MORE RNA exactly as hard as one
holding less. Making it one-sided — `−½·p·max(0, mo − log f)²`, no new constant — is **the only mechanism
the zero-gDNA control has ever ENDORSED, on every one of its rows**, and yet the panel reads a small
REGRESSION. The reason is that the two-sided term was doing TWO jobs: the bound, and by its upward side a
de-facto gDNA LEVEL channel. Removing the accident exposes the real gap — worse on exactly the stratum
with no working level channel, better where the prior already suffices — which is why the two halves are
one experiment.

⛔⛔ **AND THE LEVEL CHANNEL IS STRUCTURALLY DISCONNECTED — a theorem, not a measurement.** The chain is
strictly REGION · BOUNDARY · REGION · …, so a one-slot-step channel is BIPARTITE; the only licensed
originators of a gDNA level are structurally pure-gDNA objects, which are REGIONs, so **no REGION can ever
receive a gDNA level.** Measured at `g00`, the receiver census is exactly that shape: thousands of
BOUNDARY receivers and **exactly zero** of every REGION class, with certified RNA the mirror image (zero
BOUNDARY, every live exon REGION). ⛔ Two repairs are already refused: letting a BOUNDARY originate a
level, which patches a symptom; and letting the level cross two steps, where the chain-fused level is
dominated by the intergenic anchors and hands every exon the off-probe floor against a true neighbouring
density orders of magnitude higher.

⚠ **What to run it against**: `--messages on`, at `g05 ss0.50 capture_on` (the nearest surviving rung to
the dropped `g01` one). There every single ψ-boundary channel ablation is small while the joint one is
large — `TRAPS: all-small-singly-large-jointly`, so the four share an upstream quantity — and at the worst
error slots the self-solve with the fitted prior is nearly correct while the messages destroy it, with the
certified-RNA channel alone recovering the truth at most of them.

## §2 THE OPEN QUESTIONS, ranked by how much they would change a decision

⚠ **Ranked as questions, not as work** — §1 is the ranked work list. ⛔ A question about the DEFERRED
stratum or a deleted panel cannot change a 0.8.0 decision until something else does. ⭐ Each entry names
what would ANSWER it and, where the sparse-nascent rebuild moved it, says so.

1. ⭐⭐⭐ **THE RNA LENGTH LAW IS BIASED −4 bp BY SJ OBSERVABILITY, AND THE FIX IS ALREADY BEING COMPUTED
   AND THROWN AWAY** (diagnosed 2026-08-31, unbuilt). ⭐ **What it is worth: an attributable −18,208
   transcripts** on the ladder's `g05` capture-ON row — measured by injecting the true RNA pmf, and the
   three consumers split **scorer +156 / drain −822 / GEOMETRY −18,208**, so this is a geometry problem
   and needs NO second estimand.
   ⭐⭐ **THE MECHANISM, decomposed by stage**: the raw spliced pool runs **+10.2 bp** long (longer
   fragments cross more junctions); the sj de-tilt removes **14.3**, overshooting by **4.1**; the EB
   shrinkage is inert (0.01). ⛔ The de-tilt divides by the probability a fragment CROSSES an sj, but the
   pool is selected on the splice being **OBSERVED** — and observability FALLS with length, since the
   junction must land in a sequenced read rather than the mate gap. `fl.py` already flags the exclusion of
   `sj_implicit`; this is its unpriced consequence. ⚠ A second, smaller term: `rna_pmf` estimates the
   MATURE law while geometry wants ALL RNA (mature + nascent), worth 0.3 bp at capture-ON and 1.7 bp off
   it — the same ESTIMAND mismatch as the gDNA frame error, a third instance.
   ⛔ **REFUTED, do not retry**: the sj bank's model-free identity (`count/inv_length_sum + 1`) is
   **equally** biased (−4.6, −5.3) because it samples the same observability-selected population. No
   estimator built on spliced fragments escapes this; the escape is a different population.
   ⭐⭐⭐ **AND THE TOOL ALREADY HAS ONE.** The two-pool contrast solves for BOTH components and discards
   the contaminant `r`. That discarded law is unspliced RNA from pools with **no observability selection**,
   and measured against truth it lands at **216.45 and 216.83 against 216.48 and 216.64 — inside 0.2 bp**.
   ⭐ So the fix is to keep it and combine it with the spliced-derived law by precision, reusing
   `_couple_estimands`; no read-length model and no new input. ⛔ Unbuilt and un-A/B'd: it must be judged
   on `quant_accuracy.py` transcript error, and the ceiling for the whole RNA side is the −18,208 above.

1. ⚠ **[SUPERSEDED by the row above, kept for its measurements] SHOULD THE gDNA LENGTH MODEL RUN A SECOND CONTRAST ON THE *CROSSING* POOLS? — the two pairs have
   OPPOSITE failure modes and the pool depths are anti-correlated, which is what makes this worth
   answering** (owner question, 2026-08-31; explored, not built). The crossing pair mirrors the contained
   pair exactly in what can contaminate each — mature RNA **cannot** cross an intron|exon boundary
   (it splices there), so pool 2's contaminant is nascent, the intronic region's contaminant; and pool 3's
   is unannotated transcription running past a gene edge, the intergenic region's contaminant.
   ⭐ **Measured against the origin-split oracle, with ORACLE weights, on BOTH fl-gap sign arms**: under
   capture the crossing route BEATS the shipped contained one and by a lot (`g50` capture-ON, TV vs truth
   **0.076 against 0.136** and **+2.8 bp against −9.8**; on the other arm **0.078 against 0.182** and
   **+1.6 against −23.9**), while OFF capture it is much worse (TV 0.31–0.41 against 0.016–0.098), because
   pool 3 is starved there — 29 to 630 fragments against 0.5–12 k under capture, exactly the mirror of the
   contained pools' collapse. ⛔ **Two things block it and both are named**: ① the weight estimator does
   NOT transfer — estimated `a_2`/`a_3` came back 0.002/0.002 against a truth of 0.980/1.000 under
   capture, because the crossing opportunity is computed from the INDEX and is blind to the probe panel
   (the same defect that makes `fl_anchor_gap.py`'s `G-gdna` control read +6 % on every capture-ON row);
   ② ⛔⛔ **the substrate cannot test the premise at all** — every shadow transcript lives on `test_blank`,
   which carries NO annotated genes and therefore NO intergenic|exon boundaries, so pool 3 measures
   **exactly 1.0000** pure on every condition and the contrast silently reduces to "pool 3 alone". On real
   data an unannotated 5'/3' extension crosses precisely that boundary
   (`TRAPS: purity-is-a-property-of-the-annotation`, a third time). ⭐⭐ **MEASURED 2026-08-31, and it reshapes the question**: each pool's
   enrichment above the off-target rate (`epsilon_p = true gDNA count / rho_off·E_p^index`) is **~1.0 for
   both CONTAINED pools even under capture** (the rate is fitted on the same off-target stratum, so the
   retention cancels) and **294–338× for both CROSSING pools** — the two-strata structure is the entire
   weight failure. ⭐ **And the two crossing pools' enrichments are EQUAL: `epsilon_2/epsilon_3` = 1.002
   and 1.034 under capture**, so the weight RATIO `a_2/a_3` survives capture computed from the index alone;
   only the common LEVEL (one scalar, the boundary-stratum rate) is missing. ⭐ **What would answer it**:
   ① a shadow transcript OVERLAPPING an annotated gene edge on `test_chr`, so pool 3's contamination is
   testable at all; ② the missing scalar, for which the candidates are a one-sided fit within the
   PROBED-boundary stratum (probe membership as a BOOLEAN, which the capture audit says is trustworthy
   where the calibrated LEVEL is not) or a §6c-style reconciliation that pins the level against the
   contained estimate and takes the shape from the crossing pair's ~30× larger counts. ⛔ Non-negativity
   identification is already refuted for the contained pair — do not re-try it here.

1. ⭐⭐⭐ **WHY IS PRIOR FIDELITY ANTI-CORRELATED WITH DELIVERABLE QUALITY? — the leading answer is "it is
   not the prior, it is the messages", and there is a second candidate that must be excluded first.** At
   the worst slots the self-solve WITH the fitted prior is nearly correct and the message layer then
   destroys it, and muting the certified-RNA channel alone recovers the truth at most of them. ⛔ That
   dissection ran at a rung dropped in 2026-08-13 (nearest surviving `g05`), so **confirm on a second
   stratum before closing**. ⛔ **The alternative to exclude first**: every prior-injection arm was
   measured with the effective-length shrinkage left unsubstituted, so part of "a better prior does not
   show up" may be the RULER rather than the messages — ranks 6 and 7 own that, and
   `calibration_vs_oracle.py` is the one instrument whose patch point reaches it.
   `TRAPS: the-intermediate-is-not-the-deliverable`.
2. ⭐⭐ **What is the joint price of the splice-flux reframe and the level defect?** The brief is §1.1 and
   the work is rank 11. Neither half has an honest price alone, and one of them is already in `src/`.
3. ⭐ **Do we solve ALTERNATIVE SPLICING correctly?** The `alt_splice` toy rung exists and is
   **unverified** (`toy_harness.py --list` carries it). Cheap, and it is the only structure where several
   splice junctions share a BOUNDARY — so it is in scope on all three shipping strata.
4. ⚠ **Does the message layer's transfer variance correctly price a ratio built on a handful of counts? —
   PARTLY ANSWERED, and what remains is narrower than the question was.** `transfer_variance_audit.py`
   named the mechanism: the shipped `transfer_logvar` is a counting term plus a composition term that
   vanishes at zero composition variance, so **every term shrinks as either slot is counted more deeply**
   and a deeply-counted transport arrives essentially undamped — a hole the retired between-mode term used
   to cover. ⛔ The same audit answered the obvious substitute with a NO in both directions at once: the
   landscape's per-slot posterior is TIGHTER than the counting variance (the wrong direction), and its
   population spread over-states a fitted premise by nearly an order of magnitude because it contains the
   real enrichment and expression structure the message is trying to transport. ⭐ What survives is the
   per-hop-type premise, and the honest seed for it is that under capture the posterior LOOSENS, which is
   mode-membership ambiguity reappearing on the right substrate. `EQUATIONS.md` §3.5d is the derivation.
5. ⚠ **PARKED, not open — the two capture-ON pilot rows that disagreed about the SIGN of every length
   correction.** Both panels were deleted in 2026-08-13, so answering it needs a rebuild; the correction
   it concerns is inside the retired length channel, so it cannot change a 0.8.0 decision either way. ⛔ If
   it is ever revisited, do not average the two rows — find which one is lying. ⚠ The fl-gap side panels
   are the nearest surviving substrate and they carry the RETIRED uniform nascent model, so they are not a
   drop-in replacement for either deleted panel.
6. ⭐⭐ **DOES ANY IN-SCOPE VERDICT DEPEND ON THE NASCENT STRESS LEVEL? — newly askable, and it is the one
   question the rebuild opened rather than closed.** The ladder runs at `on_fraction 0.50`, which is a
   DEVELOPMENT STRESS level and not real data; realistic is about a fifth of that (`DESIGN.md` §0b). A
   number measured on the ladder therefore says what is ROBUST, never what is likely. ⭐ This does **not**
   need another panel: re-simulate the single worst in-scope scenario at the realistic level and check
   whether any RANK moves. ⛔ If a verdict only holds at the stress level, it is a robustness finding and
   must be labelled as one.
7. ✅ **CLOSED as WON'T-FIX, 2026-08-22** — *should the simulator capture a pre-mRNA through every probe
   whose genomic blocks it spans?* It is nascent × capture fidelity, which the NASCENT SCOPE RULING
   (`DESIGN.md` §0b) puts out of scope: nascent is sparse in real data and dedicated nascent-capture is a
   different experiment. ⭐ The rebuild shrank the residual it explains as well, since capture depletes
   nascent by roughly an order of magnitude on the current ladder, and no verdict depends on it
   (`TRAPS: the-panel-enriches-nascent-by-its-own-probes`). Re-open only if a real library forces it.

8. ⭐⭐ **WHY IS THE RELAY DISCONTINUOUS IN `od_r` AT ~1e−5?** Measured 2026-08-30 on
   `g98 ss0.50 capture-OFF`: relay's error is 217,531 at `od_r` ∈ [0, 1e−7] and **212,581** at 1e−5, a 2.3 %
   jump across a 5e−6 parameter change, with everything else held fixed. That is a THRESHOLD in the
   relay/anchor path, not a response, and while it stands no relay comparison that crosses it is
   attributable — it is how a deleted constant looked like a regression (`TRAPS:
   a-constant-parked-a-value-off-a-knife-edge`). Same family as the κ̂ − ½ residue slope at unstranded g00.
   ⛔ It is in the anchor, not in the strand estimator. ⭐ **BOUNDED 2026-08-30** — a 1e−5 nudge to `od_r`
   with everything else held fixed, across all 30 test-chromosome conditions: **1/30 rows move more than
   0.5 % for each policy**, worst 0.59 % (silent) and 1.65 % (relay). So the benchmark RESOLVES a large
   effect and not a percent-level one: ⛔ **do not believe a single-row policy difference below ~2 %**, and
   a relay comparison that crosses the edge is not attributable. Driver:
   `scratchpad/awayhalf/stability.py`.

## §3 THE RULES THAT MADE THESE NUMBERS TRUSTWORTHY

⭐ These are `TRAPS.md`'s operational core, repeated here because they are what a new session must do on day
one, not read about later.

- **Measure the CEILING before building the correction** (TRAPS: measure-the-ceiling-first). It has re-ranked this project
  **five** times, and it is also how you learn a phase is finished — §0 is that story for Stage A.
  ⛔⛔ **AND CHECK WHAT THE CEILING'S ARM ACTUALLY SUBSTITUTES** — 2026-08-14: every prior-injection
  ceiling in this file patched `assemble_priors` and left `effective_lengths_em`, built earlier in the same
  pipeline, untouched. A ceiling measures the arm you wrote, not the quantity you named.
- ⭐⭐ **RANK PER STRATUM, AND AGAINST THE 0.8.0 SCOPE** — three strata are the target, unstranded ×
  capture-ON is deferred and reported, and it carries 64.5 % of transcript error, so a panel total ranks
  the deferred stratum (TRAPS: never-pool-the-strata).
- **Run the PANEL arm before writing a mechanism into `src/`** (TRAPS: panel-before-src). Three toy-positive changes have been
  panel-negative.
- **Zero-gDNA and zero-RNA controls on every experiment** (owner, 2026-08-05). The truth is a constant, so
  every deviation is a false positive with nothing to cancel it — and **re-anchor-on-the-deliverable** was found this way.
  ⚠ But check the arm could have fired at all (TRAPS: could-the-arm-have-fired).
- **Quote `mwae` over ALL objects and the raw Σ|err|**, never the honesty metrics alone (TRAPS: honesty-metrics-reward-ignorance).
- **One thing varied**, with the baseline re-recorded from the current tree in the same session (TRAPS: re-record-the-baseline).
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire**
  (TRAPS: perturb-every-gate). And name the observable for *each place* the change was made (TRAPS: name-the-observable-per-site).

---

## §4 ⛔ WHAT IS DELIBERATELY NOT NEXT — one boundary each, so it is not rebuilt

⛔⛔ **THIS SECTION AND ITS GRAVEYARD POINT FORWARD — *do not rebuild these* — so they are not completed
items and are never deleted.** Every row is a MEASUREMENT and stays exactly as recorded.
⚠ **Read the "its target" side of a row with the 0.8.0 scope in hand:** where a mechanism's only target
was **unstranded × capture-ON**, that target is now DEFERRED, so the row is moot as a 0.8.0 candidate on
top of already being refused. ⭐ **What is NOT moot is the `g00` zero-gDNA column** — the zero control is
owner-required on every experiment and applies to all four strata, so a row that failed there is refused
for 0.8.0 on grounds the scope does not touch.
⚠ **One row is about a different thing than its name suggests:** "the RNA fragment-length model" below is
the accumulator's FL *geometry*, which the EM consumes and which 0.8.0 ships. The length-channel
retirement is of a CALIBRATION COMPOSITION channel and does not reach it.
⚠ **PANEL STAMP:** every row whose `closed by` column says "all 36 conditions", or which quotes `g01` /
`g10` / `g25` / `g75` / `g90`, was measured on the ladder RETIRED 2026-08-13. The verdicts stand as
records; re-opening one means re-running it on the 16-condition panel, not re-reading the number.

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_rescale.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` BOUNDARY** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both BOUNDARIES the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | `length_ceiling.py`, one pmf at a time | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores sj *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. ⭐ Why an evidence-free object cannot reach the vertex at all — and why that is the value of missing information rather than headroom — is `EQUATIONS.md` §9a |
| **the Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`region_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to that fix, not to this arm. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| ⛔⛔ **A SPECIFIC per-transcript allocation RULE — soft-min over exclusive objects with a per-object Jeffreys half** | built end to end and A/B'd on all 36 conditions, seed pinned | ⛔ **REFUSED — worse on EVERY stratum and on the zero control**, transcript Σ\|err\| **57.5 M → 81.6 M (1.42×)**, and a length-proportional variant **2.10×**. ⭐ **The MECHANISM is not what failed** — the gDNA:RNA split moved **+0.2 %**, exactly as the conservation identity requires, so the A/B priced the ALLOCATION alone. Three defects, each measured: exclusivity hard-zeroed **38.7 %** of transcripts; estimating a density on a tiny exclusive region and extrapolating over the whole transcript amplified variance up to **6,534×** (44.6 % of weighted transcripts had their density from <200 bp); and a per-object `+½` revived the silent half of the annotation (`frac_expressed: 0.5`), taking false-positive mass **18.6 M → 41.6 M**. ⛔ The rule and its config flag were DELETED; the LANE it rode on was kept and is now proven (§1, rank 3). ⚠ The Jeffreys-mean half of it is §4.1 graveyard row one — see `TRAPS: a-trap-names-the-defect-not-the-repair` |
| ⛔⛔ **THE SOFT-MIN-ALONG-THE-PATH WEIGHTING FUNCTION — the owner's own theorem, built faithfully and REFUSED** | 12 arms (4 modes × 3 multipliers) on `g00 ss0.99 capture_off`, 3 of them re-run on the blind stratum `g50 ss0.50 capture_on`, base re-recorded in the same session | ⛔ **REFUSED. Worse than `base` at TRANSCRIPT level on every rung of every arm and on both strata** — 1.317–1.604× at `g00`, 1.262–1.331× on the blind stratum — with transcript false-positive mass 1.76–2.20× worse. ⭐ **The one encouraging number does not survive:** at `g00` the same weights took GENE error to **0.395–0.527×** (against `oracle_alloc`'s 0.128×), and on the blind stratum that collapses to **1.006–1.041×**, i.e. nothing. ⭐⭐ **The MECHANISM is `TRAPS: an-upper-bound-is-not-an-estimate`, and it is structural rather than a tuning failure:** the theorem bounds a transcript by the thinnest object on its path, but **3,644 of 4,839 silent transcripts (75.3 %) share an object with an expressed one** and inherit its bound. The zero-weight SET was byte-identical across all twelve arms, because a bound is zero only when every object is — a property of the data, not of the estimator. ⭐ Retreating to GENE granularity (split within the gene by effective length alone) did NOT rescue it — 1.340× against 1.317× — which proves the damage was never the within-gene split. ⚠ **Two things the experiment ESTABLISHED and which the next candidate should keep:** the dial is monotone in the theorem's favour on both axes and both strata (`min` < `harmonic` < `geometric` < `arithmetic`), so the pooled `Σmass/Σopportunity` control is the WORST rung and the soft min is doing real work; and **0.0 % of expressed transcripts were ever zeroed**, so the unrecoverable failure direction never fired. `scripts/design/transcript_weights.py`, `tests/calibration/test_transcript_weights.py` (31 gates) |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |

### §4.2 ⛔⛔ SEVEN MORE, FROM THE REFERENCE INVESTIGATION (2026-08-15/16). DO NOT REBUILD THESE EITHER.

⭐ Kept SEPARATE from §4.1 because that table is eleven rules for resolving DOUBT at an evidence-free slot
and is not to be edited. These are attempts to give ψ's reference a MEAN — plus two that came out of
SHIPPING it — and each was built, measured on the panel with both zero controls, and refused.
⛔ The form that survived is `DESIGN.md` §6b.1, and it now ships **ON** — see §0 for its measured numbers.

⛔⛔ **6 — GIVING `τ_λ` THE LOCATION TERM'S CURVATURE, "the asymmetry with the intron factory".** Built,
measured, refused on three counts, and the reasoning that motivated it was wrong at every step: the 3,227×
fall in `τ_λ` at a pinned slot is **~98 % the `[f(1−f)]²` Jacobian** (nothing was lost); the contribution
is a **BOOLEAN gate flip** that releases the full COUNT precision, not a ¾-unit increment (τ = 0.029 and
τ = 1e6 both return 850.44 of a 850.50 ceiling); and it carries no count, so it credits **data-free** slots
(`n = 0` ⇒ `prec_g` 0 → 0.2026) — the very population the structural reference's safety argument rests on
being empty. ⚠ Measured, it was **bit-identical on the deliverable on all 32 panel rows** and moved only
`has_own_composition_evidence`, which is 0.8.0's own denominator.
`TRAPS: a-priors-curvature-is-not-the-datas-information`.

⛔ **7 — SOFTENING THE PRIOR TO A PER-OBJECT ONE-PSEUDO-FRAGMENT FLOOR** (`m_i = E[g]_i/(E[g]_i+1)`), the
principled reading of the derivation. Worse on **every** stratum (0.609 / 1.045 / 0.580 / 1.000 against
0.381 / 0.659 / 0.363 / 0.800): on a structurally pure-gDNA object the truth IS `f_g = 1`, and a soft floor
pulls it back off the vertex it should sit on. ⭐ What replaced it introduces no constant at all — the
lattice's own top point `σ(L)` — see `EQUATIONS.md` §9c.1.

| # | mechanism | why it was refused |
|---|---|---|
| 1 | **a fitted RNA density `logP_r`**, the mirror of the gDNA landscape | ⛔ the only non-circular form (fit from the solver's own belief) reads **0.988 / 0.997 / 1.037** — nothing, then worse: feeding ψ a density fitted from ψ's own belief tells it what it already believes. ⭐ And the ORACLE version's gain did not survive a shuffle — a shape that is wrong on purpose BEAT the true one at `g98` (0.786 vs 0.854), so the attribution was never established (`TRAPS: attribution-must-survive-a-shuffle`) |
| 2 | **a library-wide Beta mean, `a = f_lib`** | ⛔ `g05` regresses **1.43×**; `f_lib` is calibration's own output so the loop has positive feedback with both vertices attracting; and moving `a`/`b` sets the TAILS as well as the location, so `b = 0.03` leaves **57 %** of the prior outside `L = 10` |
| 3 | **the OBJECT-weighted mean instead of `f_lib`** | ⛔ **the two split by STRAND and neither wins everywhere**: object-weighted 0.584 / 0.452 on the two stranded strata but **5.570×** on unstranded × capture-OFF. ⚠ The sweep that motivated it was `ss_0.99` on three of its four rows — `TRAPS: never-pool-the-strata`, met on a sweep rather than a panel |
| 4 | **a stratified ASSERTION** — pure-gDNA strata claim `f_g = 1`, reweighted by stratum size | ⛔ reads 1.000 at every condition **including `g00`**: an assertion cannot see a library with no gDNA in it. ⭐ What replaced it is a per-object DENSITY, which needs no reweighting at all — the strata select the training set, not the answer |
| 5 | **a pooled RNA density from sj flux** | ⛔ RNA spans six decades with no genomic autocorrelation, so a pooled flux is not a population parameter (owner). It scored well only by sitting on the mass-weighted centre — `TRAPS: a-mean-hits-the-mass-weighted-centre-by-luck`. ⭐ Replaced by RNA-as-residual, which predicts no RNA at all |

⚠ **One more that is a CAUTION rather than a refusal:** `f_g ≤ 1 − S/M` as an assumption-free bound from
certified RNA. `boundary_spliced` is a SEPARATE bank from `boundary_unspliced`, not a subset, so the bound
is simply false — the truth violated it by **302**. The correct statement is the identity
`ρ_r·E_r = unspliced_RNA + S`, i.e. S SUBTRACTS.

### §4.3 ⛔⛔ THE TRUNCATION-FREE REGION BANK INSIDE THE CURRENCY CHANNEL — BUILT, PRICED, REFUSED (2026-08-20). DO NOT REBUILD IT AS A DROP-IN.

The shipped REGION bank reads `rho·P(w<=ell)` (`TRAPS: a-cancellation-is-conditional-on-its-support`),
and the truncation-free repair — `region_start_count / ell`, the STARTS-IN relation, no fragment length
in the weight — was built behind `CalibrationConfig.region_abundance_bank`, gated by an absolute fill
gate fired three ways, proven inert outside the currency arm (16/16), and priced on all 16 conditions.
⛔ **REFUSED: it improves the currency arm's pass-0 exon solve (no-evidence ≤0.05 mass coverage
45.8 % → 66.1 % at `g98 ss0.99 ON`) and regresses its DELIVERABLE — the four zero controls 2.18×
(each 2.0–3.3×), the deferred stratum 1.84×, stranded × ON 1.20×; the only stratum win is 0.843× on
stranded × OFF.** ⭐ The named suspect (one leg measured): under the shipped bank **58.6 % of live exon
hops carry a REGION bank of exactly 0** and `enrichment_ratio` returns its 1.0 default — an ACCIDENTAL
MUTE the currency arm's standing partly rested on; the live bank un-mutes them and the policy's
machinery converts a better input into a worse answer
(`TRAPS: the-intermediate-is-not-the-deliverable`, measured inside one arm). ⭐ What survives: the
truncation algebra (`EQUATIONS.md` §2), the fill gate, the two-names rename, and the wall-exposure
numbers (FLUSH 3.7–9.1 % / BINDING 7.4–22.6 % of exonic starts by template population, spliced
coordinates, MAX collapse; the full table stays in the dev sandbox until a consumer lands).
⛔ The knob was DELETED after pricing (converge-and-delete); the full implementation is one commit
before the deletion (`a2b81b34`). Re-opening this requires a policy whose DELIVERABLE improves under a
better level channel — none exists today.

### §4.1 ⛔⛔⛔ THE GRAVEYARD — ELEVEN MECHANISMS PRICED, ELEVEN REFUSED. DO NOT REBUILD THESE.

⭐ Promoted from a working doc when it was deleted (2026-08-07). Every row is a
real build that was measured and refused. ⛔ **`g00` is the owner-required ZERO-gDNA control: its truth is
exactly 0, so every fragment there is a false positive with nothing to cancel it** — which is why a
mechanism can look good on its target and still be inadmissible.

⛔⛔ **THE TABLE IS UNCHANGED BY THE 0.8.0 SCOPE AND IS NOT TO BE EDITED — it is eleven measurements.**
⚠ How to read it now: the **`its target`** column is, for most rows, a win on **unstranded × capture-ON**,
which is the **DEFERRED** stratum — so those wins are moot as 0.8.0 arguments and were never admissible
anyway. ⭐ The **`g00`** column is the one that decided each row, it is in scope for all four strata, and
it is why the pattern paragraph below is the real result. ⚠ `zc_struct_lock_g1` is the one row still live
in §1 (rank 1's feeder ②, and rank 11) — as **half of a pair**, never alone, exactly as its row says. ⚠ Its `g00`
column was RE-PRICED 2026-08-18 on the licence + SPLICE IN-fix relay: **0.71×** on top of the SPLICE IN fix, with the four
`nrna_none` zero controls 1.8–15× worse — the sign on the total flipped, the verdict (half of a pair) did not; the table
row is the 2026-08-11 measurement and stays.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ REGION` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| the mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |
| `struct_lock = g1_locked ∧ REGION` **re-priced** | the standing strict xfail's own fix, on top of the SPLICE IN precision repair | ⭐ **0.71×** | in-scope **+2.5 / +2.1 / +0.5 %**, deferred +17 % | ⛔ **2026-08-18**, and the sign on `g00` FLIPPED from the 2026-08-11 row above once the relay was fixed underneath it — but the verdict did not: the four `nrna_none` zero controls go **1.8–15× WORSE**. The empty exons' "gDNA = 0 @ 0.2026" is LOAD-BEARING at AMBIG `exon\|exon` boundaries in an RNA-only library (an RNA+ level claim with no RNA− claim drifts ψ to 0.38 without it). ⭐ Still half of a pair |
| the mass rescale refuses a ZERO-MASS source (`pinM`) | `may_share_composition` additionally requires `M[src] > 0` | ⭐⭐ **0.31×** (134 k — better than the pre-licence relay's 154 k) | unstr × OFF 0.999, str × OFF 1.004, ⛔ **str × ON 1.032, six of six worse** | ⛔ **2026-08-18.** Under capture the empty slots between probe-covered stretches are the CONDUITS a relayed composition travels through (`TRAPS: the-divergence-was-a-barrier`), so the rescale at those hops is load-bearing. ⭐ The sharp predicate is NEITHER this nor the row above: refuse a source's OWN zero-count artefact, KEEP a relayed composition passing through an empty slot |

⭐⭐ **THE PATTERN, AND IT IS THE REAL RESULT: every one of the eleven was a rule for how to resolve DOUBT,
and at `g00` the doubt must resolve to NO gDNA.** A rule that lifts an evidence-free slot off zero is
inadmissible there however well it scores elsewhere. ⛔ The only candidate the control has ever ENDORSED is
the one-sided certified-RNA bound (−81.9 %, 8/8), and that one is panel-negative alone — it is half of
`TRAPS: a-cancelling-defect-pair` and is §1's rank 11, **the-cancelling-pair** (§1.1).

⚠ Four more `zc_*` arms exist as decomposition REVERTS used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.
