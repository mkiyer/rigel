# ROADMAP — where the tool is, and what to do next

Reading order for a new session: **`SUCCESS.md`** (how performance is measured) → **this file** (current
state and priority) → **`TRAPS.md`** (mistakes not to repeat) → `DESIGN.md` / `EQUATIONS.md` as needed.

⭐ **The map of this file:** the scope block is the FRAME · **§0** where the tool is · **§1** the ranked
list and its per-rank detail, which is the point of the file · **§2** the open questions · **§3** the
rules · **§4** and **§4.1** what has already been priced and refused — read before proposing a mechanism.

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

⛔ **Every number below was measured on the tree and panel current at its stamped date.** ⛔ **The panel changed on 2026-08-13**: `pilot`, `flgap_short` and `flgap_long` were DELETED and the ladder rebuilt 36 → 16 conditions (`g00/g05/g50/g98`), so **a number without a 2026-08-13-or-later stamp has not been re-derived on the panel now on disk** — and the retired rungs `g01`/`g10`/`g25`/`g75`/`g90` cannot be re-run at all without a further rebuild. Re-derive rather than
trust — a number that has moved is a result, not a documentation bug.

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
**strand and density**, plus belief propagation across objects, which is currently off. ⭐ `TESTING.md` §0
and `SUCCESS.md` carry it in full; `TRAPS: a-length-gap-bypasses-calibration` is the name.

⭐ **Scenarios must be CACHED so calibration can be re-run extremely efficiently** (owner). `panel.py
cache` builds both caches, `build_scan_cache.py` is the scan half, and the pattern to copy is
the pattern to copy is a key made of the scan manifest plus a content hash of every producing source file,
with the module under test deliberately OUTSIDE the key so its edit loop is one second.

## §0 THE STATE

⛔ **Re-derive rather than trust.** A number that has moved is a result, not a documentation bug.

⚠ **Twenty-four rows, re-counted 2026-08-20 rather than adjusted (`TRAPS: re-record-the-baseline`). If it
grows again, PRUNE rather than let it become a report** — this table says where the tool IS, so a ✅ row
earns its place only by saying *do not work here*, and a number nobody can attribute earns none.
⛔ **Read the stratum, never the total.** Every row that names unstranded × capture-ON is describing the
**DEFERRED** stratum, and it is reported rather than ranked.

| | | |
|---|---|---|
| ⛔⛔⛔ **THE THIRD POLICY IS BUILT — AND ON THE REAL LADDER IT IS THE WORST OF THE THREE ARMS ON ALL THREE IN-SCOPE STRATA. THE TEST-CHROMOSOME RESULT DID NOT TRANSFER.** | ⛔ **MEASURED 2026-08-20 on all 16 conditions** (`relay_pool_ab.py`, axis ALL, gDNA absolute error in FRAGMENTS, contaminated rows summed per stratum, artifact `<ladder>/benchmark/ladder_three_arms_2026-08-20.tsv`), Silent / Relay / **currency**: unstranded × capture-OFF **357,580** / 496,226 / 758,085; stranded × capture-OFF **291,815** / 440,209 / 774,909; stranded × capture-ON **564,678** / 961,174 / 2,158,432; DEFERRED unstranded × capture-ON 18,559,229 / **6,021,441** / 6,095,640. ⛔ **`SilentPolicy` wins all three IN-SCOPE strata and 8 of 16 rows** — on this panel message propagation is net-harmful wherever the local solve has evidence, for BOTH policies. ⚠ **On the test chromosome the same comparison made `currency` best on all three** (662 / 617 / 4,570 against Silent's 673 / 619 / 5,632): the toy and the panel DISAGREE IN RANK, which is `TRAPS: panel-before-src` for the fifth recorded time. ⛔ **Any claim about this policy must name its substrate.** | ⭐ **What currency IS best at, on the ladder**: both capture-ON zero controls (45,432 and 5,809 against Relay's 122,976 and 42,378) and the two dense unstranded capture-ON rows — i.e. exactly where the local solve is blind. ⭐ And an exon-level split says what it trades: it moves accuracy from evidence-BEARING exons to evidence-free ones, which is the wrong trade on a panel where evidence-bearing exons carry the mass. ⛔ It stays behind `message_policy` with `relay` the default; it is a DEVELOPMENT BASELINE with a named deficit, not a shippable improvement |
| ⭐⭐⭐ **AN ENRICHMENT RATIO IS NOW MODEL-FREE — and the quantity was already on disk, unread** | ⭐ **NEW 2026-08-20**, from the owner's question *"how are we trusting our enrichment ratios now?"* A total abundance must never be `mass / effective_length`: that divisor is a function of the composition being solved for, so 100 counts in 500 bp reads 0.25 as pure gDNA and 0.33 as pure RNA. The accumulator already deposits the RECIPROCAL OPPORTUNITY (`1/(w−1)` crossing, `1/(ℓ−w+1)` contained), whose expectation is the density EXACTLY for any length distribution and any composition — it reached `CalibrationSubstrate` and stopped. It now reaches the policy as `RegionGeometry.inv_abundance` (+ the sj bank's own column, per FACE and per TRANSCRIPT strand). `EQUATIONS.md` §3.5g | ⚠ **An equal-fragment-length panel is structurally blind to this**, which is why it had to be reasoned about rather than measured. ⛔ What is still model-dependent: the per-component divisors `E_g`/`E_r`, so a COMPONENT density still uses a length model — the TOTAL no longer does |
| ⛔⛔ **THE ZERO-gDNA ERROR IS THE REFERENCE'S MEAN, MEASURED THROUGH THE SOLVER — AND FIXING IT NAIVELY COSTS THE STRANDED × CAPTURE-ON STRATUM** | ⭐ **NEW 2026-08-20.** Setting ψ's reference mean per object from the MEASURED background (`m_i = rho_bg·E_g,i/M_i`, `rho_bg` from the structurally pure-gDNA slots, which are EMPTY at a zero-gDNA library) takes the LOCAL solve's error at **both zero controls to EXACTLY 0** — on the test chromosome (95,757 / 46,241 / 24,836 / 31,288 → 0) **and on the 16-condition ladder** (1,381,959 / 416,652 / 68,431 / 31,692 → 0). On the ladder it further improves **every** capture-OFF row (unstranded 0.952 / 0.960 / 0.748; stranded 0.956 / 0.987 / 0.756) and regresses **every** stranded × capture-ON row (5.702 / 6.507 / 7.787), which is the documented 2.6–3.6× anchor under-read. ⛔ So it is NOT a candidate as it stands: it says the reference's mean must be per-object **and capture-aware**, making rank 2 and rank 3 ONE piece of work | ⚠ **It refines a recorded blocker**: rank 2's *"`g05` regresses 1.43× on both strand settings"* was measured for a LIBRARY-WIDE mean; with a PER-OBJECT one `g05` capture-OFF **improves** and only capture-ON regresses — the `g05` blocker is a CAPTURE blocker. ⭐ And it re-prices the message layer: under a corrected reference **every** policy gets all four zero controls right, Silent included, so the relay's headline control win was substantially compensating for the reference |
| **Stage A — the accumulator** | ✅ **DONE**, and that is a measurement | perfecting BOTH fragment-length models is worth **2.6 %** of the deliverable, down from 22.2 % |
| ✅ **the fragment LEDGER CLOSES — no molecule is dropped without a counter** | ⭐ **NEW 2026-08-19.** `n_fragments == deposited + Σ accounted drops`, exactly: the panel read `10,000,000 − 9,988,687 = 11,313 = 7,226 deferred + 4,087 CHIMERIC`, and the chimeric ones were counted in NO qc field and were absent from the census. ⛔ Cause: `detect_chimera` joined a fragment's mates iff their TRANSCRIPT SETS intersect, so a contiguous genomic molecule spanning two unrelated transcripts was dropped — **4,087 gDNA fragments per condition, 0.04 % of fragments carrying 2.4 % of every boundary crossing**, because the predicate fires where transcripts are short and dense, which is where a fragment crosses many boundaries. ⭐ Repaired by the owner's COMPATIBILITY-BEFORE-CHIMERA ruling (`DESIGN.md` §3.1b-ii); `n_chimeric` **4,087 → 0** and the residual is now the deferred count exactly | ⭐⭐⭐ **FIELD CERTIFICATION IS **16/16 on the capture-OFF conditions, where the gate actually runs** (it was **8/16**; all eight failures were capture-OFF). ⚠ The other 16 are capture-ON, where `calibration_oracle.py` marks the uniformity gate **VACUOUS BY DESIGN** — capture makes the field deliberately non-uniform — so "32/32 stamped" is a stamp count and NOT 32 verifications** — re-derived 2026-08-19 on a fully rebuilt panel (32/32 scan caches, 32/32 oracle caches, all 32 `slot_truth.npz` stamped COMPOSITION + FIELD, **zero** class-bias flags anywhere). ⛔ **THE 2–3 % BOUNDARY-CROSSING DEFICIT WAS THIS, AND IT IS CLOSED.** FIELD certification failed on **8 of 32** conditions (the dense capture-OFF ones); on a fresh scan the gate now passes on every surviving capture-OFF partition, worst |z| per reference going **−11.93 → −2.78**, **−11.84 → −2.16**, **−8.23 → −2.11**, **−8.05 → 1.62**, **−6.59 → −1.07**. ⚠ Three of the eight have no surviving partition BAM and are verified only by their `nrna_mid` twins; ⚠ `g05_ss_0.50_nrna_mid` ref 0 moved the WRONG way (2.45 → 4.17, inside the band). ⛔ It hid because every dropped fragment was a CROSSER, so `region_contained` matched BAM truth EXACTLY on all 22,139 regions |
| ✅ **THE SIMULATOR SIMULATES THE TOOL'S OWN TRANSCRIPTOME — nascent RNA is the index's ENTITY, not a private pre-mRNA space** | ⭐ **NEW 2026-08-19, owner-ruled.** `rigel sim` now builds/loads a **rigel index** and simulates its transcript list — annotated transcripts PLUS `create_nrna_transcripts`'s synthetic single-exon nascent entities (TSS/TES clustered at 20 bp) — in ONE multinomial with `prob ∝ abundance × capture-aware effective length`, so the mature/nascent fragment split is realised rather than imposed. ⛔ **The defect it removes**: capture was applied inside each space separately, so a pre-mRNA saw only the probes of its OWN transcript — nascent spanning another transcript's probed exon read **0.16×** of "enriched alike" while the gDNA there was fully enriched, and unprobed sibling isoforms (360 of 4,154 in probed genes) were never enriched at all. A probe now maps by GENOMIC OVERLAP to gDNA and to every transcript whose exons it touches, **either strand** (ds-cDNA at capture). ⭐ Measured after: gDNA **1162×** and nascent **1003×** under one probe — the same rate (`TRAPS: the-panel-enriches-nascent-by-its-own-probes`; `TESTING.md`). ⛔⛔ **EVERY CACHED CONDITION IS NOW STALE** — mature capture-ON truth moves too, so all 32 need re-simulating before any capture-ON or nascent number is quoted | ⭐ Four perturbations fire: probe scoped to one transcript, nascent left on the contributor, the split penalty dropped, the GTF path restored (which now raises rather than silently simulating the old model). Suite **3,538 passed / 9 xfailed**; four golden scenarios regenerated with the diffs read (`effective_length` ≤ 0.7 %, `nrna_em_count` +6 to +20 %) |
| ⭐⭐⭐ **THE STANDING BENCHMARK — the scenario table in COUNTS, and the numbers the new policy must beat** | ⭐ **NEW 2026-08-19**, owner-specified: *"we need to track the full scenario table as a benchmark that we will gradually improve on … the table needs to report COUNTS (not ratios or fractions) … ground truth, estimated counts, and then different measures of error (net error, absolute error)."* The instrument is **`relay_pool_ab.py`** (the owner agreed it is the right one); the artifact is `<ladder>/benchmark/baseline_2026-08-19_shipped.log`. Σ **gDNA absolute error in FRAGMENTS** over the 16 rebuilt conditions, `SilentPolicy` → `RelayPolicy`: ⛔ **`g00` zero control 1,898,734 → 235,733 (0.124×)**; unstranded × capture-OFF 357,580 → 496,226 (**1.388×**); stranded × capture-OFF 291,816 → 440,209 (**1.509×**); stranded × capture-ON 564,678 → 961,174 (**1.702×**); ⛔ DEFERRED unstranded × capture-ON **18,559,229 → 6,021,441 (0.324×)**. ⭐⭐ **This RE-PRICES the relay's standing on the rebuilt panel** — the ~1.5–1.65× on the in-scope strata and the large control/deferred wins both survive, now measured with nascent RNA present in EVERY condition, so no row is `nrna = 0` restated. ⚠ `abs` is summed over OBJECTS (misplaced mass) and can exceed the 10 M library; `net` is the signed total and is printed beside it. | ⭐ Three truth pools are printed on every row (`gdna`/`nrna`/`mrna`), so a gDNA over-call names which RNA it ate. ⛔ `--donor` injects LIBRARY-level priors for a toy reference and REFUSES a donor whose strand/capture axes differ from the conditions — a mismatch is silent and once reported 82,581 false positives that a matching donor puts at 0 |
| ✅ **THE HOP CURRENCY MAP, PER ARM — Stage 2 of the message-propagation rebuild, MEASURED on all 16 (`hop_currency.py`)** | ⭐ **NEW 2026-08-19.** Every ordered adjacent pair scored with the source's TRUE value transported as a LEVEL and as a COMPOSITION, per component `{gDNA, RNA+, RNA−}`, in fragments, beside a Monte-Carlo noise floor; 64/64 gates green, 16 conditions in 23 s. ⛔⛔ **THE ARMS DISAGREE, which is why the owner ruled the map must be per arm**: at `R exon <- B exon\|exon[term]` off capture the gDNA arm wants a **LEVEL** (0.0 % excess over its floor) while both RNA arms want a **COMPOSITION** (5.2 % / 4.1 %) — the pooled first pass called that hop LEVEL outright and hid it. ⛔⛔ **AND THE HOP TYPES ARE NOT THE SEVEN OBJECT CLASSES** — a `B exon\|intron` is a SPLICE SITE (15,480) or a TERMINUS (3,867) and they have opposite currencies (`TRAPS: an-object-class-does-not-see-a-terminus`); the key is `object class × {sj, term, sj+term}`. **The map:** COMPOSITION into an exon from a SPLICE-SITE boundary, all three arms (0.4–1.1 % under capture, at floor off it); at a TERMINUS into an exon, gDNA takes a LEVEL and the RNA arms a COMPOSITION; out of an exon, gDNA a LEVEL and RNA a COMPOSITION under capture; an intron into its own boundary is at floor. ⛔ **NEITHER currency serves `R exon <- B gene edge[term]` UNDER CAPTURE — ~10 % excess on every arm**, the residual census and `ROADMAP.md` §1 rank 3's spike-and-slab, re-established rather than inherited. Depth: 23–29 % of imputed mass ≥ 9 hops from measured gDNA off capture. | ⭐ The per-arm truth is the oracle cache's two extra partitions `rna_pos`/`rna_neg` (the RNA reads by TRANSCRIPT strand), gated slot-by-slot against the mature/nascent split — **`n_rna_pos + n_rna_neg == n_mrna + n_nrna`, measured max gap 0** |
| ✅ **the tally CONSERVES a fragment count** | ⭐ every deposited fragment places **exactly one** unit across the objects it crosses — regions, contiguous boundaries and sj boundaries together. Measured on the origin-split oracle: **1.000× deposited, 0 unaccounted**, on BOTH origins | ⛔ RNA read **0.747×** until `sj_mass` landed (2026-08-11): a spliced fragment whose blocks cross no boundary deposited on no conserved bank at all — **1,222,375 of 4,830,713 (25.3 %)**. gDNA was always exact, because gDNA cannot splice |
| ✅ **ψ's COMPOSITION CLOSES — STRUCTURALLY** | ⭐ **NEW 2026-08-17.** `simplex_logodds._compose` builds the composition as the IMAGE of ψ's two parameters (`λ → f_g`, RNA total `:= 1 − f_g`, `θ →` the tilt share) instead of reading three coordinates independently, and `_posterior_median_fg` is a CONTINUOUS ½-quantile read on the λ lattice. Measured on real conditions: **100.00 % of published objects close**, both axes, every class, min/p5/p95 exactly 1.0000 — against **74.7 % / 77.2 %** before. Panel per stratum **0.999 / 0.995 / 0.996**, deferred 1.000, `g00` control **0.995**, pass-0 unstranded **0.993** | ⚠ Inside the 0.996–1.013 noise floor, so read it as neutral-to-marginally-better on accuracy; the WIN is correctness and simplicity (`DESIGN.md` §6c). ⛔ Taking MEANS everywhere also closes and is REFUSED at 1.352 / 1.573 / 3.756 and 1.801 on the control |
| ⭐⭐⭐ **ψ's λ BRACKET WAS TOO NARROW TO EXPRESS ITS OWN PRIOR — BUILT, GATED AND PRICED (`derive_logodds_window`, ships OFF)** | ⭐ **NEW 2026-08-17.** `landscape.logprior` is evaluated at `log ρ = log f + log M − log E` and ψ can only offer `f ∈ [σ(−L), σ(L)]`; at the shipped `L = 10` that floor is **4.540e-05**, **370×** above the median density the prior points at (`1.225e-07` over 16,885 slots carrying mass at `g00 ss_0.50 off`). ⭐⭐⭐ **THE BRACKET IS DERIVED AND NO CONSTANT IS CHOSEN — `DensityLandscape.required_logodds_window` is `max_i log((1−f_i)/f_i)` with `f_i = ρ_floor·E_i/M_i`, which collapses to the landscape's own log-dynamic range `log_rho[-1] − log_rho[0]`**: predicted **22.764 / 21.013 / 21.804 / 21.073** against a measured per-slot max of **22.76 / 21.01 / 21.80 / 21.07**. | ⭐⭐ **MEASURED THROUGH `calibrate` ON ALL 16 LADDER CONDITIONS, Σ\|Δ\| in object-incidence fragments, arm ÷ noop, PER STRATUM.** Stratum totals **0.4339** unstranded × OFF · **0.9075** stranded × OFF · **0.9674** stranded × ON · 0.9847 deferred — ⭐ **all four improve**. Per rung `g00 / g05 / g50 / g98`: unstranded × OFF **0.3449** / 0.9917 / 0.9840 / 0.9766; stranded × OFF **0.3242** / 0.9914 / 0.9820 / 0.9732; stranded × ON **0.3951** / 0.9972 / ⛔ **1.0077** / 0.9739; deferred 0.2796 / 0.9210 / 1.0004 / 1.0000. ⛔ **ELEVEN of the twelve IN-SCOPE conditions improve and ONE regresses** (`g50 ss_0.99 capture_on`). ⭐⭐ **`g05` is the tell: it regressed 1.43× under EVERY library-wide reference mean tried (§4.2 row 2) — the blocker that killed that family — and here it improves on both strand settings.** ⚠ On the DEFERRED stratum the arm is inert at `g50`/`g98` (1.0004 / 1.0000), which is that stratum's known blindness and not a property of the bracket. ⭐ Gates: the noop reproduces the shipped tree EXACTLY on all 8 conditions with a recorded reference; the bracket demonstrably FIRED on every condition (`TRAPS: an-ablation-that-never-ran`); and the resolution-only control moves the OTHER way (1.0147× / 1.0046× / 1.0036×), so the effect is the bracket and not the lattice. ⚠ **UNPRICED: memory at genome scale (~2.3× on `sweep_n_grid`, and profile on real cfRNA, never this panel) and the end-to-end thermometer.** ⛔ The intron factory's λ-factor is REBUILT at the widened bracket, never regridded — there is no map onto a domain it was never evaluated on, and it raised `IndexError` until that was seen |
| ⭐⭐⭐ **THE EXON REFERENCE IS SOLVED AT CAPTURE-OFF BY *ONE* PER-OBJECT FORMULA — AND THE ZERO CONTROL GOES TO EXACTLY 0** | ⭐ **NEW 2026-08-17, measured through `calibrate` on the widened bracket.** `DESIGN.md` §6b's ruling written literally — `m_i = ρ_g,i·E_g,i / M_i`, RNA as the residual, ρ_g pooled PRIOR-FREE over the anchors mature RNA cannot occupy — replacing BOTH shipped branches (the 0.75 four-class pin AND the neutral ½) with a single expression. Σ\|Δ\| arm ÷ base: `g00 ss_0.50 off` **654,786 → 0**, `g00 ss_0.99 off` **9,626 → 0** — ⭐⭐ **EXACTLY ZERO phantom gDNA at both zero controls, the test all ELEVEN §4.1 graveyard mechanisms failed** — then `g05` **0.870**, `g50` **0.745 / 0.731**, `g98` **0.628** on capture-OFF, both strand settings. | ⛔⛔ **IT IS NOT SHIPPABLE YET: stranded × capture-ON reads 2.33×, and the cause is NAMED, MEASURED AND NOT IN THE ESTIMATOR.** One pooled ρ_g reads **8.97×** there because probes deplete the off-target floor ~30×; splitting into the two anchors whose RATIO is the enrichment (in-gene `exon\|intron` for on-target, intergenic+intron REGIONs for off-target) takes it **8.97 → 2.33** and leaves capture-OFF untouched. ⭐ The residual matches `object_composition.py`'s recorded number: the in-gene anchor **under-reads on-target gDNA 2.6–3.6×** because it sits at the EDGE of the probe footprint while `eff_gdna` is built with unbounded reach. ⛔ **That is an OPPORTUNITY-MODEL repair (layer 2), not an estimator one** — fix the divisor, not the reference. ⭐⭐ **The win is the UNIFICATION, not a patch**: applying the formula ONLY at the exonic slots and keeping the pin elsewhere is WORSE than base (1.03–1.06× at `g05`/`g50`). ⚠ Priced on 7 conditions, ONE of them capture-ON — extend before landing. ⛔ It is only re-openable at all because the λ bracket now represents `m_i ≈ 1e-07`; the earlier refutation clamped it and measured the clamp |
| ✅ **the SUITE is green with no skips** | ⭐ **0 failed / 3,501 passed / 0 skipped / 10 xfail**, re-derived 2026-08-17 on a CLEAN TREE at `4cfd35d9`. ⚠ The `11 → 10` xfail is an xfail CLOSED BY REPAIRING THE THING: `test_scripts_index.py` xfails only while a script is BOTH listed in `BROKEN_ON_IMPORT` and genuinely broken, so repairing `prior_units_check.py` made the gate FAIL on its stale exemption — which is the behaviour it was written for. That dict is now EMPTY; it was **3,476** after `structural_reference` was defaulted ON, **3,459** before that, and **3,420** on 2026-08-15 after the object-composition work and **3,412** on 2026-08-14 after the ruler work and **3,404** earlier that day and **3,397** on 2026-08-13 and **3,235 / 0 / 7** on 2026-08-11, and every delta is ACCOUNTED file by file in `CLAUDE.md`, never adjusted. Two xfails and both skips were closed by BUILDING the missing thing, never by widening a bound | ⛔ The 7 interrogated on 2026-08-11 are not one kind of thing: 2 are the recorded price of `message_propagation = False`, 5 were WRITTEN as xfail as records of two proven defects whose fixes are panel-negative alone. `CLAUDE.md` carries that split; ⚠ the three added since have not been interrogated the same way |
| ✅ **ONE NUMERIC CONVENTION** | ⭐ a COUNT is an integer, a FRACTION is float64. No fixed point, no scale constant, nothing decodes a bank (owner, 2026-08-11) | ⭐ float64 is **1e5–7e5× MORE accurate** than the fixed point it replaced, measured against exact rational arithmetic. ⛔ The ~2.6 % once quoted against float is a **float32** number — `TRAPS: integer-channels-reproduce` |
| ⭐ **calibration, the THREE IN-SCOPE strata** | ✅ median library `f_gdna` error **0.005–0.012**, and ⭐ the PRIOR the EM reads is within **2.5–4.6 %** of a perfect one | stranded × capture-ON/OFF and unstranded × capture-OFF — **the 0.8.0 target**. ⚠ 36-condition ladder, RETIRED 2026-08-13 |
| ⛔ **calibration, unstranded × capture-ON — the DEFERRED stratum** | ⛔ **BLIND** — reports **0.033–0.058** while truth spans **0.00 → 0.98**, and hands the EM a gDNA prior **94.4 %** short. ⭐ Re-measured 2026-08-13/14 on the rebuilt ladder at the exon objects: `f_g` **0.040 / 0.0016 / 0.0021** at `g05 / g50 / g98` against truths **0.054 / 0.518 / 0.982** — a NEAR-ZERO answer regardless of truth, which **looks acceptable at low gDNA by coincidence** | not noisy; a flat boundary. ⛔⛔ **DEFERRED for 0.8.0** (owner, 2026-08-14) — reported on every panel, not a development target until the other three are optimised. The length channel was tried and refused, and is retired until after 0.8.0 (the scope block above) |
| ⛔⛔ **the EFFECTIVE-LENGTH SHRINKAGE at the ZERO-gDNA control** | ⛔ **NEW 2026-08-14.** The shipped shrinkage contracts every transcript by a mean factor **0.345** when the correct answer is exactly **1.000**, on **13,673 of 15,669** transcripts — **INCLUDING on capture-OFF, where the module contract says the factor must be 1**. `rho_ref` is fabricated entirely from false-positive gDNA | ⭐ **It is a SYMPTOM of the composition error, not an independent bug**: substituting ONLY the composition arrays and re-running the SHIPPED shrinkage function gives the correct factor (`g00` capture-off **0.345 → 1.000**; `g50` unstranded capture-on **0.834 → 0.401**, the truth). ⛔ So repair the composition it reads, not the function. One split, two consumers, one function — `priors.py` imports `_global_reference_density` **from** `capture_eff_length.py` |
| ⛔⛔⛔ **THE RULER IS NOW MEASURABLE, AND ON THE PANEL A PERFECT ONE IS A *LOSS* ON THE TWO CAPTURE-OFF IN-SCOPE STRATA** | ⭐ **PRICED ON ALL 16 CONDITIONS, 2026-08-14.** `oracle_ruler` substitutes at the `calibrate` boundary so it reaches BOTH consumers; `oracle` substitutes all three `LocusPriors` fields, so **`oracle_ruler` − `oracle` is the shrinkage and nothing else**. Transcript Σ\|err\| against `base`, per stratum — `oracle` / `oracle_ruler`: stranded × OFF **1.022 / 2.510**, stranded × ON **0.936 / 1.077**, unstranded × OFF **1.005 / 2.678**, ⛔ deferred **0.267 / 0.259**, `g00` control **0.928 / 0.207**. Gene level the same shape (1.697 / 0.816 / 2.047 / 0.104 / **0.060**); false-positive mass **2.955× and 3.340× WORSE** on the two capture-OFF strata and **0.074×** at the control | ⛔⛔ **THE ONE-CONDITION READING WAS WRONG AND THE PANEL IS WHY WE KNOW.** On `g00 ss_0.50 capture_off` alone a perfect ruler took transcript Σ\|err\| to **0.043×** where a perfect prior alone reached 0.802×, and that invited "every ceiling was measuring a fifth of the win". ⛔ It is not: the win is **confined to the zero control and the DEFERRED stratum**, and on the two in-scope capture-OFF strata a perfect ruler is **2.5–2.7× worse**. `TRAPS: panel-before-src` and `TRAPS: a-single-level-panel-cannot-see-a-constant`, met on a ceiling instead of a mechanism |
| ⭐⭐ **WHY A "PERFECT" RULER LOSES — the `U` null already said so and it was read too weakly** | ⭐ At capture-OFF there are no probes, so the physically correct contraction factor is **exactly 1.000**, and `calibration_vs_oracle.py`'s no-enrichment null confirms the estimator can reach it (**0.9963–0.9967**). But BOTH arms sit far below it: shipped **0.9318 / 0.9391** and *oracle* **0.9241 / 0.9193** — ⛔ **the shipped ruler is CLOSER to the truth than the oracle one**, so `oracle_ruler` substitutes one wrong ruler for a slightly worse one | ⛔ **So `oracle_ruler` is NOT a ceiling at capture-OFF; it is an A/B between two wrong rulers.** ⭐ The real ceiling arm is the **`U` ruler** (uniform gDNA field ⇒ factor ≈ 1.000), which is derivable with no fitting, and it has not been run. ⚠ Note also that a ~1–2 % aggregate gap costs 2.5× at transcript level, so `Σeff/Σfl` HIDES a large per-transcript redistribution — read `ruler_n_moved`, not the aggregate |
| ⭐ **the RULER's own error, per stratum** | ⭐ **NEW 2026-08-14**, `calibration_vs_oracle.py`, `factor = Σeff_em/Σfl` (1.000 = no contraction): **0.9318 / 0.0794 / 0.9391** on the three in-scope strata against O's 0.9241 / 0.0797 / 0.9193 — ⛔ and **0.0951 against a truth of exactly 1.000** at the `g00` control, 1.06 **billion** bp of opportunity on a library containing no gDNA | ⭐ **The no-enrichment null `U` settles that the shrinkage is a SYMPTOM, not a second defect**: O's own gDNA total laid down at exactly uniform density — the physically correct capture-OFF field, contract 1.000 — reads **0.9963–0.9967**, so perfect composition leaves only ~0.4 % of manufactured contraction. ⚠ **Two aggregates disagree by 3.6×**: the recorded **0.345** is `mean(eff/fl)`, this is `Σeff/Σfl`; the contraction falls hardest on LONG transcripts. Rank on `Σeff/Σfl` |
| ✅ **ψ's REFERENCE MEAN IS SET FROM THE ANNOTATION, AND IT IS ON** | ⭐ **NEW 2026-08-16, `CalibrationConfig.structural_reference = True`.** `m → (a+1)/(a+b+1) = 0.75` wherever `¬mrna_active`, neutral ½ elsewhere. Measured through `calibrate` — `vertex_ceiling.py --arm config_struct` vs a `base` re-recorded in the same session — final Σ\|Δ\| per stratum **0.930 / 0.908 / 0.925** on the three IN-SCOPE strata, 0.998 deferred, `g00` control **BIT-IDENTICAL on all 8 rows** at pass-0, final and `conf_wrong_err`. Better on all 12 contaminated conditions, worse on none; pass-0 **0.922 / 0.955 / 0.973**; confidently-wrong error **0.956 / 0.962 / 0.783 / 0.501**. Thermometer transcript **0.994 / 0.995 / 0.995**, deferred 1.002 | ⛔ **THE STRENGTH, NOT THE LOCATION, WAS THE WORK — and the ladder cannot rank a strength.** `strength = logit(m)`, so a location IS its strength in nats. The first form took it from the LATTICE (`m = σ(L)` ⇒ 9.31 nats ≈ 10,000:1) and scored **0.384 / 0.660 / 0.366** here — ⛔ but **worse than NO PRIOR at being refuted** (2.0247 vs 0.3946), because the ladder holds `nrna = 0` and therefore scores only the DELIVER obligation, where more nats is monotonically better. ⭐ One pseudo-observation on the reference's own exponents gives 0.75 exactly, a 3:1 claim refuted by ~1.5 stranded fragments and by `τ_fac = 161.4` unstranded. ⚠ Gene-level thermometer deltas (2.5–3.5 k) are BELOW the ~6,700-fragment noise floor and are not attributable; transcript-level (10–26 k) are |
| ⛔ **the PER-TRANSCRIPT prior lane is BUILT and NEVER PASSED** | ⛔ **NEW 2026-08-14.** `rna_prior_weight` is built end to end and the production call site in `pipeline.py` **omits it**, so the EM's default rule carries **zero per-transcript information** | ⚠ Passing it needs a weighting function, and the two built so far were REFUSED (§4). The lane, its arms and its falsifications are all in place |
| ⭐ **what a PERFECT prior is worth, by GRANULARITY** | ⭐ **NEW 2026-08-14.** Per-**LOCUS**: gene-level error on the deferred stratum **0.099**. Per-**TRANSCRIPT**: the three in-scope strata to **0.37–0.58** gene-level, with false-positive mass cut **100–200×** | ⛔ Per-transcript does **NOT** rescue unstranded × capture-ON (**0.641**) — which is the deferred stratum, so that is a reported fact and not a blocker. ⭐ The in-scope prize is the per-transcript one |
| ✅ **the prior ASSEMBLER and its POPULATION** | ✅ `rel` **0.0019–0.0027** with perfect masses in and **4.9e-4** with perfect per-component shares; the composition claim `a_g:a_r` is exact against the unspliced pool (`Δphi` **≤ 5e-4**) | ⭐ the assembler was **0.179**: the conserved-count rewrite took it to 0.0202 and the YARDSTICK took it to 0.0027 (`TRAPS: score-the-consumers-own-count`) |
| ✅ **the fragment-length models** | ✅ accurate — `pi(w)`'s de-tilt reads **211.77** against a true **212.20**, and the gDNA pmf is exact to **0.02 bp** off capture | ⛔ a claimed +10.7 % bias was a truth-parser bug (`TRAPS: a-truth-table-of-aggregates`) |
| ⚠ **the gDNA pmf under capture** | ⛔ the last unfixed length defect: **+13.6 bp** at a 330 bp gDNA mean and **+3.5 bp** at 120, drained per-bin `fit/true` **1.22 … 4.18** in the tail | ⭐ EXACT off capture and untouched by the drain, so it is a PLACEMENT problem — `EQUATIONS.md` §4.4 — and must not be attacked by editing `gdna_opportunity` |
| **message propagation** | ⛔ **OFF**, and it stays off until the tool is optimised end to end across all scenarios (owner, 2026-08-10) | net better on **all three IN-SCOPE strata** (−58 / −44 / −32 %, 16/16 on the two stranded ones); **+155 % on unstranded × capture-ON**, which is the **DEFERRED** stratum. ⭐ **The `/16` is SCORED ROWS, not conditions, and it does fit the 36-condition ladder** — settled 2026-08-17 by reading `arm_score.py`, whose rows are keyed `(condition, axis)` and whose every per-stratum line predicates on `not is_g00(k)`: 9 rungs − the control × 2 axes = exactly 16 rows per stratum. ⚠ On the 16-condition ladder the same arithmetic gives **6**, so a re-price reading `n/6` has not lost coverage. ⚠ **The 0.8.0 scope changes how that reads and does NOT re-open the ruling**: the whole recorded price lands on the stratum 0.8.0 does not target, so the honest question is now what it costs the ZERO controls. ⛔ RE-price on the rebuilt ladder, never inherit — every number in this row predates the conserved-count rewrite and the numeric convention |
| **the price of that, and it IS in scope** | ⚠ zero-gDNA golden scenarios go **0.029 → 89.93** and **0.005 → 9.58** | both AMBIG loci. ⛔ A zero control is owner-required on every experiment, so this price is in scope even though the stratum it sits on is not |
| ✅ **end to end — the LIBRARY figure** | ⭐ **RE-DERIVED 2026-08-11 on all 36 ladder conditions — ⚠ that ladder was RETIRED 2026-08-13 and this has not been re-derived on the 16-condition one**, messages OFF, `calibration_truth_ab.py --cache-subdir _main`. Mean `\|err\|` per stratum: **0.0025 / 0.0053 / 0.0056** (the three IN-SCOPE strata) and ⛔ **0.4040** (max **0.9151**) on unstranded × capture-ON, the DEFERRED one; `g00` control **0.0497** | ⛔ **NEVER quote the pooled 0.1043** — it is 97 % the deferred stratum (`TRAPS: never-pool-the-strata`). ⚠ The retired `0.1060` turns out to have been close, but it was unattributable and the stratified form is what ranks work. ⭐ The three in-scope strata are ~2× better than the "0.005–0.012" row above, measured the same day |
| ⛔ **end to end — TRANSCRIPT assignment** | ⛔ **31.3 %** of RNA fragments misassigned (**55,916,700** of 178,399,996, 32 contaminated conditions), and a perfect prior removes **32.2 %** of that, to 21.2 % | ⭐ **RE-DERIVED 2026-08-11 on the 36-condition ladder RETIRED 2026-08-13**, messages OFF, EM seed pinned — and it lands within 0.2 pp of the pre-rebuild 31.1 %, so the rebuild did not move it. ⛔ The ceiling is ONE STRATUM: unstranded × capture-ON goes **0.405**, the other three are **1.031 / 0.981 / 1.013** — a perfect prior is NEUTRAL-TO-WORSE on them. At GENE level the split is sharper still (**0.215** vs **1.704** on stranded × capture-OFF) |
| ⭐ **the pool level — gDNA vs NASCENT vs ANNOTATED** | ⭐ **2026-08-11**, `quant_accuracy.py --report` table ⑥ — ⚠ on the 36-condition ladder RETIRED 2026-08-13. gDNA total is **−0.6 % / −6.4 % / −0.5 %** on the three IN-SCOPE strata and ⛔ **−87.1 %** on the deferred unstranded × capture-ON; a perfect prior takes that last one to **−5.5 %** | ⛔ **Read `gDNA (total)`, never `gdna_est`** — the latter is `gdna_em_count` and EXCLUDES intergenic fragments, which are more than half of off-capture gDNA. Scoring it against the origin-count truth reads **−50.7 %** panel-wide that is entirely the missing pool; `score_library`'s docstring records the same mistake once fabricating an off-capture under-call (0.3151 vs a truth of 0.5000) |
| ⚠ **the `noop` BYTE-IDENTITY GATE CANNOT PASS — and the arm is SOUND anyway** | `--arm noop` differs from `base` on **376 of 2,016 fields**, max **399** on `count_abs_err`. ⭐ **The noise floor `--arm base_reseed` is 1,317 fields and 6,668** — so the arm's delta is **17× BELOW the floor** on the same field and touches 3.5× fewer fields. The injection wrapper is doing nothing, as designed | ⛔ **So the gate's SPECIFICATION is wrong, not the arm.** Byte-identity is unreachable while the scan is multi-threaded (`TRAPS: the-deliverable-is-not-reproducible-by-default`'s second source — the seed is pinned, the float banks still re-associate per run), so `arm_identity.py` is the wrong instrument for a `quant_accuracy` arm. ⭐ Either re-specify it as "below the `base_reseed` floor" or make the scan merge order-independent; ⚠ until then **no `quant_accuracy` arm delta below ~6,700 fragments is attributable** |
| ⛔ **the NASCENT channel** | ⛔ **20.2 M fragments** parked on entities whose truth is exactly **0** — 9.2 % of all true RNA, invisible in every transcript table | ⚠ measured 2026-08-11 on the ladder RETIRED 2026-08-13, and it PREDATES the 2026-08-12 rule that gives a **synthetic** nascent entity ZERO RNA prior (in-scope nascent **6,238,406 → 691,796**; the derivation, the `is_synthetic`-not-`is_nrna` distinction and the strand-gated +4.21 M gDNA price are `EQUATIONS.md` §9b). ⛔ Re-derive before quoting either number. ⭐ Next is an ACCURATE per-entity nascent prior rather than withholding from all (owner, 2026-08-12) — §1 rank 3 |
| ⛔ **reproducibility — TWO independent sources** | ⛔ `EMConfig.seed` defaults to `None`, **and** the scan is multi-threaded so the float banks are re-associated by the per-worker merge | `TRAPS: the-deliverable-is-not-reproducible-by-default`. Pinning the seed alone is no longer sufficient. ⭐ Owner signed this off (2026-08-11): the spread is **1.503e-15** on `posterior_mean` and **exactly 0** on every integer column, and tests validate the float banks within a DERIVED, bracketed tolerance |

⭐⭐ **THE ONE SENTENCE, AND UNDER THE 0.8.0 SCOPE IT HAS THREE PARTS.**

**① On the LIBRARY gDNA fraction, the tool is accurate on all three IN-SCOPE strata and blind on the
deferred one** — at unstranded × capture-ON `κ = ½` makes the strand λ-term identically 0 and no other
channel reaches an AMBIG slot. ⛔ **That blindness is now DEFERRED, so §1 is no longer "give that slot its
own evidence".** The search for a θ-independent channel is CLOSED (the scope block above), and the one
candidate that came out of it — the length channel — is RETIRED until after 0.8.0.

**② On TRANSCRIPT-LEVEL assignment that sentence is FALSE, and it was never tested until 2026-08-07.** The
tool misassigns 15.6–20.6 % of fragments even on the three strata calibration handles perfectly, and
42.9 % at `g00` where there is no gDNA at all to get wrong. A perfect prior removes ~32 % of the panel
total and makes two strata *worse*. **These are two different problems in two different files, and
conflating them is how "is calibration the bottleneck?" got asked for months without an answer.**
⚠ For 0.8.0 this number is the THERMOMETER, not the target — see the scope block above.

**③ The 0.8.0 metric is the CALIBRATION RESULT against ORACLE CALIBRATION on the three in-scope strata,
and two open defects reach it**: the effective-length shrinkage (which no ceiling has ever priced, because
it is built before the function every arm patches) and the per-transcript prior lane that production never
passes. Both are 2026-08-14 findings and both are in §1's top three.

⛔ **Read `mwae_all` / `Σ|err|`, never `solv%` / `mwae` / `conf-wrong` / `calib`** — those four share a
denominator the solver moves by declining to answer. And quote the SHIPPED column, not pass-0: a −37.2 %
pass-0 win was −3.9 % shipped (TRAPS: the-intermediate-is-not-the-deliverable).

## §1 ⭐⭐⭐ THE 0.8.0 RANKED LIST — WHAT HAS TO HAPPEN, IN ORDER

⛔⛔ **RANKED BY WHAT MOVES THE 0.8.0 METRIC — the calibration result against oracle calibration, on the
THREE IN-SCOPE STRATA, per stratum.** ⛔ Never by the panel total: the total is 64.5 % the DEFERRED
stratum, so ranking on it ranks the thing 0.8.0 does not ship.

| rank | the work | why it is there |
|---|---|---|
| **1** | ⭐⭐⭐ **THE THIRD POLICY IS BUILT AND IS THE BASELINE TO EVOLVE — what remains is CHARACTERISATION, not construction** | ⭐ **Stages 0–3 are DONE** (vocabulary, unblock, the map, the policy). `messages/currency.py` beats both shipped policies on all three IN-SCOPE contaminated strata (§0's row), on the same gated backbone, with `RelayPolicy` untouched and the A/B a config value. ⛔ **This is a baseline, not a finish** (owner, 2026-08-20: *"This is the beginning, not the end"*). What is next is to CHARACTERISE the residual: **which exons are solved accurately, by observable properties that exist on real data where there is no truth** — that set is the TRAINING SUBSTRATE the gDNA landscape needs, and the whole bootstrap (§0c.0d) turns on it. ⭐ Then: retire `RelayPolicy` once the panel says it is dead, fold the remaining operators down, and pursue a per-hop-type premise fit in place of the single library scalar | ⭐⭐ **The rulings are `DESIGN.md` §0c.0b–d and the derivations `EQUATIONS.md` §3.5f–h**; the lessons are `TRAPS.md`'s four new §D rules. ⚠ **The measurement discipline that found every defect here, and it is reusable**: split the error by whether the destination HAD ITS OWN composition evidence — "the messages help" and "the messages trample a measurement" are different findings, and a pooled number cannot tell them apart |
| **1b** | ⭐⭐ **EXPAND THE gDNA SPECTRUM — deliberately, and not by multiplying the cross** | The panel holds four gDNA levels (g00/g05/g50/g98); the owner has ruled the spectrum needs filling (*"one percent, five percent, ten percent, twenty five percent, all the way up to above ninety percent"*) and warned against blowing it up (*"we can blow this up and suddenly have dozens or hundreds of benchmarks"*). ⛔ Levels are informative where behaviour CHANGES, not where it is smooth — interpolation is free. So the choice must be justified by a measured transition, and the new levels should cross a REDUCED set of the other axes unless a measurement shows those axes interact with gDNA level | ⚠ The honest cost is real: simulate + scan-cache + oracle-cache + certify per condition. Pay it when the current scenarios are exhausted (owner), which the rank-1 characterisation is what establishes |
| **2** | ⭐⭐ **WHAT THE NEW POLICY MUST SATISFY — the per-bug loop's findings, kept as CONSTRAINTS rather than as a queue** | ⛔ **These are no longer tickets against `RelayPolicy`.** They were converted from mysteries into mechanisms by the dissection loop and they are what the rebuild must not re-acquire: ① **a refused claim must lose its precision with its value, and no later operator may hand it back** (`TRAPS: zero-the-precision-with-the-value`) — the SPLICE IN once added the sj count's precision onto an arm the licence had zeroed, and one hop later the mass rescale read that confident zero as "all your mass is gDNA"; ② the population of a message is **DIRECTION-DEPENDENT** — boundary→exon INCLUDES the spliced fragments that splice in, exon→boundary EXCLUDES them (`DESIGN.md` §0); ③ a boundary where the RNA populations differ may transfer a DENSITY, never a COMPOSITION; ④ both **ZERO CONTROLS** on every rung, because `g00` is ONE-SIDED — any RNA-favoring message reads as "right" there. | ⭐ **THE BASELINE THE NEW POLICY MUST BEAT is the relay as it stands**, measured on the 32-condition ladder against certified truth (`relay_pool_ab.py`): a large win on unstranded × capture-ON and at every zero control, a residual cost on the three in-scope contaminated strata. ⚠ **RE-PRICE IT — do not quote the old figures.** Every published relay number predates BOTH the SPLICE IN precision fix AND the chimera repair, and the panel has been rebuilt underneath them. ⛔ Still-open items belonging to the OLD policy — the `struct_lock` feeder, `transfer_var`, `strand_population` keep-or-cut — are deliberately NOT ranked: they are properties of a policy scheduled for deletion |
| **3** | ⭐⭐⭐ **THE EXON REFERENCE UNDER CAPTURE — SPIKE-AND-SLAB, the only gap left in the reference** | ⭐ capture-OFF is SOLVED (a per-object `m_i = rho_g,i*E_g,i/M_i`, locally imputed: **0.4037–0.4977** vs base on the contaminated capture-OFF rows, `g00` to 0) and off-target objects are SOLVED (`m_i` 0.9878–0.9916 against a truth of 1.0000). ⛔ **Exons under capture are the ONLY remaining gap**: probes enrich exons, so every gDNA anchor sits outside or at the edge of the footprint and reads **2.6–3.6x low**. Derivation: `EQUATIONS.md`'s capture-reference section. Ruling: `DESIGN.md`'s. | ⭐⭐⭐ **A Beta CANNOT put its median at `f_g = 0`; a SPIKE-AND-SLAB has an ATOM and CAN** — so this is the only proposed mechanism that reaches the vertex, and §9a's "value of missing information" is an artifact of the prior FAMILY, not of the data. ⛔ Gate: it must reduce EXACTLY to the shipped form at capture-OFF (anchor ratio 0.98 ⇒ the slab collapses onto the spike), plus both zero controls and a SHUFFLE control |
| **4** | ⭐⭐ **CONFIRM THE LANDSCAPE THEN TRAINS ON A REAL SUBSTRATE** | Nothing to build — the payoff, and why the exon matters at all. With exons BRACKETED rather than pinned, `_fit_gdna_hyperprior` trains on objects whose `f_g` reflects data and can fit the bimodal `P(log rho_g)` capture actually produces | ⚠ The circularity to watch: it trains on `belief.f_g*mass` over expressed REGIONs INCLUDING exons, so too strong a reference means it learns the prior back. ⭐ At `g50`/`g98` capture-ON **~80 % of the training mass has its own composition evidence** — real data to learn from, provided the reference does not overwrite it. Measure the no-evidence share before and after |
| **5** | ⭐⭐⭐ **PERFORMANCE — MANDATORY BEFORE 0.8.0** (owner, 2026-08-17) | The grid solve in memory-bounded parallel CHUNKS, with an advanced CLI flag spanning one-object-at-a-time through many-objects-in-parallel. Compute-vs-memory is the dial | ⛔⛔ **OPTIMISE ON HIGH-DEPTH REAL RNA-SEQ, NOT cfRNA** — cfRNA is too sparse and small to be the target. ⚠ This REVERSES the standing "profile on real cfRNA" instruction, which still appears in several places. ⚠ `scripts/profiling/` is the natural home, is covered by NO gate, and 2 of its 3 files were dead on 2026-08-17 |
| **6** | ⭐⭐⭐ **BUILD THE `U` RULER ARM — the oracle ruler is NOT the ceiling at capture-OFF** | the ruler arm is built and the panel is run (§0). Its verdict is that a *perfect-composition* ruler is **2.5–2.7× WORSE** on the two capture-OFF in-scope strata, because at capture-OFF the correct factor is **1.000** and the oracle ruler (0.919–0.924) is *further* from it than the shipped one (0.932–0.939). ⭐ The arm that would actually price the ruler is the **uniform-gDNA null**, measured at **0.9963–0.9967** and derivable with no fitting |
| **7** | ⭐⭐⭐ **THE SHRINKAGE AT `g00`, REPAIRED UPSTREAM** | 0.345 where the truth is exactly 1.000, on 13,673 of 15,669 transcripts, **including capture-OFF where the contract says 1** — and it is a SYMPTOM: feed the shipped function correct composition arrays and it returns the correct factor. ⛔ Repair the composition, not the function |
| **8** | ⭐⭐⭐ **PASS THE PER-TRANSCRIPT PRIOR LANE — and the weighting function is the work** | built end to end, omitted at the production call site, so the EM carries **zero** per-transcript information. The prize is the largest in-scope one measured: **0.37–0.58** gene-level on the three strata and false-positive mass **100–200×** down |
| **9** | ⭐⭐ **RE-DERIVE THE IN-SCOPE NUMBERS ON THE 16-CONDITION LADDER** | almost every number in §0 was measured on the ladder retired 2026-08-13. Noise floor **0.996–1.013** — an arm inside that band moved nothing |
| **10** | ⭐⭐ **THE f32 CUBE MANUFACTURES A STRAND TILT AT κ = ½ — one dtype, and the falsification is already written** | ⭐ **NEW 2026-08-17**, found while bug-hunting ψ's closure. At κ = ½ the strand mean is `½` identically, so the RNA tilt is PROVABLY UNIDENTIFIED and `w_pos` must be exactly ½; the AMBIG cube evaluates that sum in **float32**, where it departs from 1 by ~1e-7 τ-DEPENDENTLY, giving a spurious tilt that grows linearly in depth. Measured `\|w − ½\|`: 7e-6 at 1e4 fragments, 3.4e-4 at 1e6, **5.3e-3** at 1e7, **0.199** at 1e8, **0.425** at 1e9. ⭐ **Falsified decisively: recompute the same cube with the strand term in float64 — everything else unchanged, cube still stored f32 — and `w_pos → 0.500000000` at every depth.** ⛔ κ = 0.9/0.99 unaffected; specific to the degenerate-κ ridge, i.e. the `ss 0.50` half of the ladder | ⚠ **NEGLIGIBLE AT PANEL SCALE and entirely PRE-EXISTING**, which is why it is rank 6 and not higher: 10 M fragments over ~70 k slots is ~140 per slot (largest ~1e5–1e6), so ≤ 3.4e-4 there. It would bite a very deep library with a dominant locus. ⭐ The repair is one dtype on the strand term inside `_solve_ambig_logodds`, NOT the whole cube — the f32 cube is an authorised memory/perf choice and the reductions already accumulate in f64. ⛔ Price it on the panel before landing (`TRAPS: panel-before-src`) and check the AMBIG cube's memory does not move |
| **11** | ⭐ **PRICE THE CANCELLING PAIR TOGETHER** | five xfails go green iff it lands, and neither half has an honest price alone. ⛔ It can only be priced with `--messages on` |
| **12** | the index's duplicate map · `mass_*_boundary` → `count_*_boundary` · restore the moment tests · the cheap ledger | correctness and hygiene, each its own commit, `TRAPS: one-thing-varied` |
| **13** | ⏳ ⭐ **FINISH THE ORACLE-EFFECTIVE-LENGTH DIAGNOSTIC** | ⚠ **STARTED, NOT FINISHED (2026-08-17).** Half an hour; it re-ranks the two ruler items below it. Needs a `git stash` arm for the pre-closure tree or it measures nothing — one arm is not an A/B | ⚠ Demoted from rank 0 on 2026-08-17: it ranks POST-calibration items, and calibration's own path is now ranks 1–4 |

⛔ **NOT ON THIS LIST, DELIBERATELY:** the length channel (retired until after 0.8.0 — the scope block
above) and anything whose only target is unstranded × capture-ON (deferred). §4 is everything else
already priced and refused; read it before proposing a mechanism.

⛔ **The blocks below are the DETAIL, one per rank, in rank order.** A rank with no block is fully stated
by its table row.

⭐ **ZERO INSTRUMENTS ARE DARK.** The last two — `flgap_study_cache.py` and `prior_yardstick.py` — were
DELETED on 2026-08-17 along with the `flgap_short`/`flgap_long` configs, rather than left waiting for a
panel. ⛔ **Neither could be repointed at the ladder, and that is structural rather than incidental**: the
whole premise is a DRAINED oracle, which the ladder refuses. ⚠ `prior_yardstick.py` was the only place the
prior assembler was scored the way production runs it — its verdict is kept in `EQUATIONS.md` §3b and
`TRAPS: score-the-consumers-own-count`, and re-deriving it means **designing a drainable length-gap
panel**, not restoring a file. ⭐ `boundary_q_population.py` and `prior_vs_oracle.py` now run with no flgap
arm at all: the dead rows were removed and the surviving measurement reproduces its recorded numbers
exactly (`q_g` 0.6330 vs `q_r` 0.5233, all four gates green).

⭐⭐ **AND "EVERY OTHER INSTRUMENT UNDER `scripts/design/` RUNS" — WHICH THIS PARAGRAPH USED TO SAY — WAS
FALSE: TEN WERE DEAD. ALL TEN WERE REPAIRED ON 2026-08-17, each falsified by perturbation, and the suite
is 3,501 passed / 10 xfail / 0 failed.** The history below is kept because it is the evidence for the new
NAMED rule, not because anything is still broken. `certified_rna_audit.py`, `toy_dissect.py`, `toy_ceiling.py`,
`toy_trace_error.py`, `zero_controls.py` died on `KeyError: '_uni'` and `reframe_walk.py` on
`KeyError: 'rho_lo'`; the sweep then found four more (`backbone_parity.py`, `held_flux_census.py`,
`prior_units_check.py`, `boundary_q_population.py`) plus two dead `scripts/profiling/` files that NO gate
covers. `_uni` is written **only** at `messages/relay.py:1021`, i.e. under `RelayPolicy`, and
`CalibrationConfig.message_propagation` defaults **`False`** — none of them set the policy.
⛔⛔ **`zero_controls.py` is the control the owner requires on EVERY experiment**, so nothing in this
file could be honestly priced while it printed one table and then crashed.
⭐ **The suite stayed green throughout, and the mechanism is NEW**: the TEST readers install the policy
themselves (`tests/calibration/test_gdna_scale_rule.py:54`), so
`TRAPS: a-green-suite-hid-five-dead-instruments` — whose recorded trigger is a `src/` **deletion** — did
not cover this. The trigger here was a **CONFIG DEFAULT FLIP**.
⛔ `psi_channel_ablation.py` was the worst of them: all four arms are structurally inert under
`SilentPolicy` and it raised `IndexError` at `:132`. ⭐ Repairing it exposed TWO number-movers — its
scorer over-read the BOUNDARY axis **1.79×** (certified-spliced folded into an unspliced-only fraction:
2,182,619 spliced against 2,760,340 unspliced at `g05 ss0.50 capture_on`), and `combine[-1]` attributed
the FINAL arm's channels to `--arm pass0`, which disagree in SIGN (−9.5 % vs +102.5 % on `− gdna_imp`).
⚠ The `2.3×` first recorded here was superseded by that measurement.
⭐ **FIXED, and the cause was bigger than the symptom:** `module_census.py` printed
*"⛔ DEAD: nothing imports it, anywhere"* for the live `strand_summary` (imported at
`src/rigel/pipeline.py:164`) because its importer scan was TEXT, not AST — so every
`from .calibration.X import` site in `src/rigel/` was invisible and **nine** modules were understated,
one all the way to zero. It now walks the AST.

- ⏳ **OPEN · RANK 1 — THE ARM IS BUILT AND FALSIFIED; WHAT IS OPEN IS THE PANEL RUN AND THE RE-PRICING**
  (built 2026-08-14). `quant_accuracy.py --arm oracle_ruler` wraps `rigel.calibration.calibrate` as a
  module attribute — `run_pipeline` imports it function-locally, so the name resolves at call time — and
  the substituted `CalibrationResult` therefore reaches **both** consumers: `transcript_capture_eff_lengths`
  (the ruler) and `assemble_priors` (the prior).
  ⭐⭐ **`oracle_ruler` − `oracle` is the shrinkage and NOTHING else**, because `LocusPriors` has exactly
  three fields and `oracle` already takes all three from O. That is what makes it an attribution rather
  than one more composite.
  ⛔ **The counter watches the SHRINKAGE, not `calibrate`.** `calibrate` being called is necessary and not
  sufficient: hand `_setup_geometry_and_estimator` a `None` calibration and the arm silently degrades into
  `oracle` under another name. The arm additionally requires the substituting variant to have MOVED the
  ruler (measured `max\|Δ\| = 1,115,202` bp) and the noop not to have moved it at all (measured exactly
  `0.0`) — `TRAPS: an-ablation-that-never-ran`. ⚠ The noop reproduces `base` exactly on `count_abs_err` at
  both axes; `fp_mass` differs by **2 fragments**, ~4 orders below the `base_reseed` floor, which is the
  §0 row's point about byte-identity being unreachable for a `quant_accuracy` arm.
  ⭐ **A second thing had to be built for it to reach the control at all**: `quant_accuracy.load_oracle`
  now falls back to `scan_cache/<condition>` when the oracle cache holds no `_main`. `pass0_vs_oracle.py`
  holds every zero-gDNA condition out, so it never wrote one, and **no oracle arm could reach `g00` before
  this**. `_main` is the undrained full payload — the same quantity — and `from_parts` re-runs sum-to-full
  either way.
  ⛔⛔ **THE PANEL IS RUN, AND IT INVERTED THE ONE-CONDITION READING — §0 carries the table.** A perfect
  ruler is a large win at the `g00` control (**0.207×** transcript, **0.060×** gene) and on the DEFERRED
  stratum (0.259×), and a **2.5–2.7× LOSS** on the two capture-OFF in-scope strata. ⭐ The reason is in
  the instrument's own null column: at capture-OFF the correct factor is **1.000**, and the shipped ruler
  (0.932–0.939) is CLOSER to it than the oracle ruler (0.919–0.924) — so this arm is an A/B between two
  wrong rulers rather than a ceiling.
  ⭐⭐ **WHAT IS OPEN: the `U` ARM.** Substitute the uniform-gDNA null ruler — O's own gDNA total laid down
  at exactly uniform density, factor **0.9963–0.9967**, no fitting and no constant — and price THAT. It is
  a ceiling at capture-OFF in a way the oracle ruler is not. ⚠ It is admissible only at capture-OFF: with
  probes the field genuinely is not uniform.
  ⚠ **And read `ruler_n_moved`, never the aggregate factor.** A 1–2 % gap in `Σeff/Σfl` costs 2.5× at
  transcript level, so the aggregate is hiding a large per-transcript redistribution.
- ⏳ **OPEN · RANK 2 — ψ's COMPOSITION REFERENCE IS A Beta(a,b), AND ITS MEAN WAS NEVER CHOSEN**
  (2026-08-15). `a·log f_g + b·log(1−f_g)` on the λ grid **is exactly Beta(a,b) in `f_g`** — verified
  numerically to six decimals — so `a` and `b` are PSEUDO-COUNTS. A Beta has two degrees of freedom and
  `_JEFFREYS_REF = 0.5` fixes both: strength `a+b = 1` (one prior pseudo-fragment, correct) and **mean
  `a/(a+b) = ½`, which asserts the library is half gDNA**. That assertion is wrong on most libraries.
  ⭐ **THE ERROR IT CAUSES:** objects whose TRUE `f_g ≥ 0.999` carry **49–83 %** of all calibration error
  on the three in-scope strata, read **0.13–0.23 BELOW** the vertex; and the `[0, 0.2]` bucket is
  over-called on every stratum. One pull-toward-the-middle, both ends.
  ⭐⭐ **BUT ~¾ OF THE VERTEX SHORTFALL IS IRREDUCIBLE AND THAT MUST NOT BE FORGOTTEN.** Driving the solver
  on a synthetic pure-gDNA object, `1 − f_g ∝ n^(−1/2)` — ordinary statistical resolution, NOT a
  structural block, and the grid window is not binding (`L` 10 → 20 → 40 moves it ~2e-4).
  ⛔⛔ **BOTH OF THOSE CLAUSES ARE PRIOR-FREE AND STRANDED, AND NEITHER SURVIVES ON THE UNSTRANDED
  STRATUM WITH A FITTED LANDSCAPE LIVE (measured 2026-08-17).** The `n^(−1/2)` slope reproduces at
  **−0.5221** only in the `n ≤ 50` window it was measured on and at the STRANDED κ; at the real
  unstranded κ = 0.500369 it is **−0.0000 in every decade** over 5.4 decades (1.0 → 263,621
  fragments/slot), so the law makes no prediction there at all. And *"the grid window is not binding"*
  is false of the shipped pipeline: with the landscape live, `L` moves the `g00` answer by **3.2×**
  (1,898,257 → 614,587 at `dlam` fixed). ⭐ The bracket is rank 1's mechanism and §0 carries it — a
  session reading these two clauses would skip the one arm that works.
  The panel's
  vertex objects have **median 9–58 fragments**, and the `n^(−1/2)` law alone accounts for **72–83 %** of
  the observed shortfall. The reference is a MULTIPLIER on that law (13× at `n = 100`, 1.9× at `n = 1e6`),
  and that multiplier is the addressable quarter. `EQUATIONS.md` §9a's "value of missing information".
  ⭐ **THE CANDIDATE:** keep the strength at 1 and set the mean to the library's own composition,
  `a = f_lib`, `b = 1 − f_lib`. It reduces EXACTLY to the shipped constant at `f_lib = ½` — measured
  0.997 / 0.999 at `g50`, the null case returning the null answer. Panel **0.198×** over 8 capture-OFF
  conditions; `g00 ss0.50`, the worst IN-SCOPE condition at 1,906,953 fragments, falls to **0.002×**.
  L-robust at `L` = 10 / 20 / 40. ⚠ **That robustness claim is about THAT ARM and does not generalise
  to the shipped solve** — see the `L` correction above; on the deliverable, `L` is worth 3.2× at `g00`.
  ⛔⛔ **THREE THINGS BLOCK IT, ALL MEASURED:** `g05` REGRESSES **1.43×** on both strand settings;
  `f_lib` is calibration's own OUTPUT so the loop has POSITIVE feedback with both vertices attracting;
  and the optimum sits ABOVE the true `f_lib`, which says the target is the **OBJECT-weighted** mean
  composition rather than the fragment-weighted `f_lib`. ⭐ All three point at one diagnosis.
  ⛔⛔⛔ **THAT THIRD BLOCKER IS NOW MEASURED AND ITS DIAGNOSIS WAS WRONG — THE TWO CANDIDATE TARGETS
  SPLIT BY STRAND AND NEITHER WINS EVERYWHERE** (2026-08-15, `object_composition.py` +
  `vertex_ceiling.py --arm ref_c=a,b`, all 16 conditions, `base` re-recorded in session). The two
  quantities do differ by an order of magnitude — object-weighted **0.6497** against `f_lib` **0.0669**
  at `g05` capture-OFF, 0.7920 vs 0.5762 at `g50` — but fed to ψ from TRUTH on both sides, Σ\|Δ\| in
  FRAGMENTS per stratum reads **`f_lib` / object-weighted**: stranded × OFF 0.705 / **0.584**,
  stranded × ON 0.692 / **0.452**, ⛔ unstranded × OFF 0.907 / **5.570**, deferred 0.170 / 0.397.
  ⭐⭐ **The object-weighted mean wins on both STRANDED strata and loses 5.6× on unstranded ×
  capture-OFF.** ⛔ The sweep that produced "the optimum sits above `f_lib`" was `ss_0.99` on three of
  its four rows — `TRAPS: never-pool-the-strata`, reached on a SWEEP rather than a panel. ⭐ The
  mechanism: at `κ = ½` the strand λ-term is identically zero, so the reference's mean IS the answer
  rather than a nudge, and it is scored in fragments while it applies per object (`g05 ss0.50 off`
  library `f_g` final **0.1696** against a truth of 0.0479).
  ⭐⭐⭐ **WHAT SURVIVES IS THE ZERO CONTROL, AND IT IS THE LARGEST NUMBER IN THE STUDY: 0.003× —
  a 333-fold reduction — on ALL THREE arms.** The shipped mean of ½ is catastrophically wrong when the
  library contains no gDNA, and any estimate of the truth repairs it.
  ⭐⭐ **AND THE QUANTITY IS ESTIMABLE PRIOR-FREE, BY TWO ROUTES THAT FAIL IN OPPOSITE PLACES.** The
  pooled gDNA density of the intergenic + intronic REGIONs is **exact at `g00`** (the anchors are empty,
  so it reads 0) and within 0.06–0.12 at capture-OFF, and ⛔ collapses under capture, where probes
  deplete the off-target anchors ~30× (`rho` **0.000216** against a true **0.007303**) — which shows up
  in the arm as **1.731×** on stranded × capture-ON exactly as predicted. The certified sj flux is
  stable across capture and reads ⛔ **0.62–0.72 at a control containing no gDNA**.
  ⛔ **The owner's sj-boundary proposal is answered and the answer is ONE SLOT OVER**: `sj_count/eff_sj`
  overstates RNA at its own BOUNDARY by **7.1–10.0×** (the flux is a density on the SPLICED template;
  the boundary's RNA divisor is the unbounded-reach UNSPLICED-crossing one), and recovers the true RNA
  density at the ADJACENT EXON REGION to **1.020–1.281** off capture and **0.740–0.997** on it, on all
  16. ⚠ Intronic REGIONs are pure gDNA here only because `nrna = 0`, and they move the estimate by
  **< 0.002** — resolution, not accuracy.
  ⭐⭐⭐ **AND THE WHOLE SEARCH WAS IN THE WRONG SPACE — THE REFERENCE DOES NOT HAVE TO BE LIBRARY-WIDE**
  (owner, 2026-08-15). ψ solves one object at a time and the gDNA arm's fitted term is ALREADY per slot
  `(n_slots, K)`; the reference is the only scalar left in ψ. Derived from two Gamma rate priors, the
  induced prior on the λ grid has a THIRD term the code does not carry —
  `− (a+b)·log(f_g + r_i·(1−f_g))` — and **the shipped reference is exactly the degenerate case where
  the two components' scales match**. With `m_i = ρ_g,i·E_g,i / (ρ_g,i·E_g,i + ρ_r,i·E_r,i)` and Jeffreys'
  `a = b = ½` kept, it is `−log[(1−m_i)·f_g + m_i·(1−f_g)]`. ⭐ It reduces EXACTLY to the shipped constant
  at `m_i = ½`, keeps the `e^(−|λ|/2)` tails so `L`-invariance is untouched, is proper for every
  `m_i ∈ (0,1)`, introduces **no constant**, and makes `R_own = m_i` a closed form rather than the
  hard-coded ½.
  ⭐⭐ **MEASURED, prior-free, `Σ|m_i − f_g,i|·M_i` in FRAGMENTS against the shipped ½, per stratum:**
  **0.123 / 0.401 / 0.123 / 0.401** on the four strata and **0.096–0.103** at the `g00` control, against
  a truth-fed ceiling of 0.081 / 0.346 / 0.081 / 0.346 — ⭐ **8× at capture-OFF, 2.5× at capture-ON, 10×
  at the control, and within 1.5× of the ceiling everywhere.** `g05`, which regressed 1.43× under every
  library-wide mean, reads **0.191**.
  ⭐⭐⭐ **The per-object variation is carried by the OPPORTUNITY GEOMETRY, not by per-object densities**
  — swapping the pooled RNA density for the per-object sj flux is WORSE (0.515 vs 0.123). Two pooled
  scalars plus exact per-object geometry is the whole mechanism: no landscape, no substrate selection,
  no fit. ⭐ Four of six strata go to **0.000** (the ones the annotation determines); the entire residual
  is the two EXONIC strata.
  ⛔⛔ **THE BOUNDARY ANCHOR MUST BE CURATED TO `exon|intron`, NOT "has a sj"** (owner): mature RNA
  crosses an `exon|exon` boundary freely, which is why the uncurated pool read `f_g = 0.0000` over
  955,428 fragments at the zero control. Curated, it reads `ρ_g` **0.050249 against a true 0.050327**.
  ⭐ The same `mrna_active` predicate then gates the RNA density and takes `R intron` 0.498 → **0.000**
  and `B exon|intron` 0.485 → **0.000**.
  ⭐⭐ **Hybrid capture needs no detection step**: the ratio of the two anchors IS the enrichment —
  **0.98** without probes, **113–114** with, a 116× separation with no threshold and no flag. ⛔ The
  in-gene anchor is a detector and not yet a calibrated LEVEL: it under-reads on-target gDNA by
  **2.6–3.6×** under capture, because it sits at the EDGE of the probe footprint while `eff_gdna` is
  built with unbounded reach. That is an opportunity-model repair, not a new estimator.
  ⭐⭐⭐ **STAGE 0 IS RUN (2026-08-16) AND THE PRIOR-FREE `m_i` NOW READS 0.040 / 0.326 / 0.040 / 0.326
  AND EXACTLY 0.000 AT BOTH ZERO CONTROLS** — 25× better than the shipped ½ at capture-OFF and 3× at
  capture-ON, with `g05` at **0.009**. Three things moved it: dropping a spliced-count subtraction that
  over-corrected at `exon|exon`, shrinking the local sj flux toward the population rate by ONE
  pseudo-observation (raw local flux with a fallback is 0.121; shrunk is 0.040, and the naive
  `rho_r = 0`-where-absent form reads **0.088–0.130 at the zero control** where the shrunk one reads
  0.000), and letting the reference speak EVERYWHERE rather than only where the annotation determines the
  answer (restricting it is 10× worse — at pass-0 there is no landscape, so its exonic constant beats ½).
  ⛔⛔ **THE CLASS-POOLED TRUTH ARM IS NOT A CEILING**: the prior-free form beats it on every stratum,
  because that arm hands an on-target RNA density to objects mature RNA cannot occupy.
  ⛔ **AND THE WIN IS NOT PER-OBJECT RESOLUTION AT EXONS.** Within `R exon`, `m_i` has sd **0.0021**
  against a true `f_g` sd of **0.4441**, and within `B exon|exon` its sd is exactly 0 — it is a better
  CONSTANT there. The four ANNOTATION-DETERMINED strata go to exactly 0.000; exons, `exon|exon` and AMBIG
  are precisely the population the gDNA LANDSCAPE exists to serve, so the two terms partition the object
  universe rather than competing on it.
  ⭐⭐ **The post-solve on-target update DOES NOT RUN AWAY.** Iterated with each object's own likelihood
  removed — a strict upper bound on the feedback — four starts spanning three decades converge to the
  same fixed point to six decimals on all twelve contaminated conditions and to exactly 0.000 at the
  control. The per-object geometry damps what a library-wide scalar could not.
  ⭐⭐ **THE TERM IS IN `src/` AND IS INERT** (2026-08-16). `CompositionPriors.location` (a per-slot
  scalar) and `simplex_logodds._location_term` are shipped, added at BOTH ψ call sites; `None` ⇒ the term
  is not written ⇒ **bit-identical**, proven on real data (19/19 output arrays, with a perturbation that
  moves 10/19). 36 gates in `tests/calibration/test_reference_location.py`. ⛔ **Nothing sets it**, so
  behaviour is unchanged.
  ⭐⭐ **THE PANEL ARM IS RUN AND `struct` IS THE CANDIDATE** — the reference MEAN pinned near 1 wherever
  mature RNA cannot be, which is exactly the four annotation-determined classes (intergenic REGION,
  gene-edge BOUNDARY, intron REGION, exon\|intron BOUNDARY = 33,347 of 70,176 slots, true `f_g` **1.0000**,
  **zero** fragment cost, and EMPTY at `g00`). Final Σ\|Δ\| per stratum, ratio to base:
  **0.381 / 0.659 / 0.363 / 0.800**, control **1.000**, `noop` byte-identical 16/16.
  ⭐ Its strength needs no constant: the term is worth `log(1/eps)` nats, and capping at the highest grid
  point (`eps = sigma(-L)`) makes it exactly `L` — reproducing `struct` to three decimals (0.384 / 0.660 /
  0.366 / 0.800). ⛔ The deconvolve on non-structural slots is REFUSED: **5.3x worse** on stranded × capture-ON.
  ⛔⛔ **AND ONE BLOCKER STOPS THE DEFAULT BEING FLIPPED.** With the reference forced on, a hand-built
  pure-gDNA intron goes **0.9006 → 0.7661** — the wrong way. The local solve is RIGHT (`fg_loc`
  0.9858 → 0.9998); what inverts it is `tau_lam` collapsing to **exactly 0**, and `tau_lam > eps` is the
  predicate for *"this slot has no own composition evidence"*. ⭐ The strand term legitimately vanishes at
  the vertex (`c*a^2 ∝ f_g^2(1-f_g)^2`), and a FLOOR is already refused by
  `has_own_composition_evidence`'s own docstring. ⭐⭐⭐ **The repair is an asymmetry, not a rule**: the
  location term is a λ-FACTOR and must contribute its curvature via the shipped
  `density_factor_precision`, exactly as `intron_prior` already does. The working design lives in the dev
  sandbox.

- ⏳ **OPEN · RANK 3 — THE SHRINKAGE IS WRONG AT THE ZERO-gDNA CONTROL, AND IT IS A SYMPTOM** (measured
  2026-08-14). At `g00` the shipped shrinkage contracts **13,673 of 15,669** transcripts by a mean factor
  **0.345** where the correct factor is exactly **1.000** — ⛔ including on **capture-OFF**, where the
  module's own contract says a uniform enrichment must give exactly 1. `rho_ref` is fabricated **entirely**
  from false-positive gDNA.
  ⛔⛔ **AND THAT SAME FALSE-POSITIVE MASS IS ALREADY MEASURED AT THE PRIOR, ON AN IN-SCOPE STRATUM.** The
  shipped prior claims **2,067,637 gDNA fragments in libraries containing none**, and **1,707,321** of them
  are at unstranded × capture-**OFF** — a stratum that reads healthy (`rel` 0.046) on every contaminated
  row, and worse by **5.5×** than the blind one. ⭐ Only the `g00` zero control could have found it, and it
  is the same defect surfacing twice: once as a prior count, once as a wrong effective length.
  ⭐⭐ **The repair is upstream and that is measured, not assumed**: substitute ONLY the composition arrays
  and re-run the **SHIPPED** shrinkage function and the factor comes out right — `g00` capture-off
  **0.345 → 1.000**, `g50` unstranded capture-on **0.834 → 0.401** against a truth of 0.401. ⛔ So do not
  edit the shrinkage function; fix the composition it reads.
  ⭐ **One split, two consumers, one function**: `priors.py` imports `_global_reference_density` **from**
  `capture_eff_length.py`, so the transcript-level ruler and the per-locus gDNA effective length are the
  same derivation reached from two places. A repair lands once.
  ⭐⭐ **AND THE "SYMPTOM" CLAIM IS NOW MEASURED FROM THE OTHER SIDE TOO, NOT ONLY ARGUED.**
  `calibration_vs_oracle.py`'s `U` arm is O's own gDNA total laid down at EXACTLY uniform density — at
  capture-OFF the physically correct field, where the module contract demands a factor of exactly 1.000.
  It reads **0.9963–0.9967**, so under perfect composition the estimator manufactures only ~0.4 % of
  spurious contraction from sampling noise. ⛔ **A separate shrinkage repair therefore remains the wrong
  move**, and this is the number that says so rather than the reasoning.
  ⭐⭐⭐ **WHERE THE TARGET IS, AND WHETHER A CHANNEL EXISTS THERE — dissected 2026-08-14.** The largest
  IN-SCOPE composition error on the whole panel is **`g00 ss_0.50 capture_off`, region mwae 0.2558**: 6×
  the next in-scope number and **67× its own stranded twin** (0.0038) on identical geometry and identical
  (zero) gDNA. Split by population, truth 0 everywhere:

  | population | regions | FALSE gDNA at ss 0.50 | at ss 0.99 |
  |---|---|---|---|
  | `g1_locked` — RNA **inadmissible** | 1,312 | **0** (it holds no mass at all) | **0** |
  | TS_POS | 14,774 | 533,654 (`f_g` 0.2720) | 4,099 (0.0021) |
  | TS_NEG | 13,808 | 382,935 (0.2485) | 3,787 (0.0025) |
  | TS_AMBIG | 5,241 | 255,829 (0.2368) | 9,355 (0.0087) |

  ⛔ **NO channel reaches an RNA-admissible region there, and it is structural.** Strand is dead
  (`rna_sense_frac` fits 0.500369, so the λ-term is identically 0); the stranded twin gets 130× / 99× /
  27× lower error using precisely that channel, and its own worst residual is **TS_AMBIG** — the class
  strand cannot resolve even at ss 0.99, which is the mechanism confirming itself. The pure-gDNA anchors
  are EMPTY, and every rule for making an evidence-free slot speak is §4.1's graveyard.
  ⭐ **One thing found that is NOT one of the eleven:** the `g1_locked` zero is a Poisson observation over
  **50.8 Mbp** of support, so it bounds ρ from above with no fitting and no constant — 95 % UB **5.9e-08**
  against a fitted `gdna_density_global` of **1.9e-02**, a factor of **322,000**. ⭐ It falsifies cleanly:
  the same population predicts its own observed count to **5.6 %** at `g05` and **2.1 %** at `g50`, so it
  is informative rather than structurally empty. And it pulls the level DOWN, where all eleven refused
  mechanisms lifted an evidence-free slot UP.
  ⛔⛔ **BUT IT IS A DIAGNOSTIC AND NOT A LEVER, AND THAT MUST BE SAID BEFORE ANYONE BUILDS ON IT.**
  `derive.py` computes `gdna_density_global` **from** the per-object solve; its only consumers are
  `cli.py`, a log line and the QC report, and `density_model.py`'s docstring states there is no
  density→deconv→density feedback loop. Repairing the scalar would change a reported number and nothing
  else — a real lever has to sit inside the solve, which is exactly where the eleven died.
  ⚠ **And the bound is admissible only at capture-OFF**: `anchor_opportunity_census.py` measured the
  empty-anchor density claim false by **346×** under capture and true off it, because probes concentrate
  gDNA. That still covers two of the three in-scope strata.
- ⏳ **OPEN · RANK 4 — THE PER-TRANSCRIPT RNA PRIOR: THE FOUNDATION IS BUILT AND PROVEN; THE WEIGHTING
  FUNCTION IS THE WORK** (owner, 2026-08-12/13). The end goal is a prior that says WHICH transcript a
  locus's RNA pseudocounts belong to.
  ⛔⛔ **AND THE LANE IS NOT PASSED IN PRODUCTION** (found 2026-08-14): `rna_prior_weight` is built end to
  end and the production call site in `pipeline.py` **omits it**, so the shipped EM falls back to its
  default rule and carries **zero** per-transcript information. Passing it is not the work — having
  something admissible to put in it is.
  ⭐ **The prize, measured 2026-08-14 and stated per stratum:** a perfect PER-TRANSCRIPT prior takes the
  **three IN-SCOPE strata to 0.37–0.58** gene-level and cuts false-positive mass **100–200×**; a perfect
  PER-LOCUS prior takes the DEFERRED stratum to 0.099 and the per-transcript one does **not** rescue it
  (0.641). ⭐ So the per-transcript lane is the largest in-scope lever measured, and the locus-level one
  is a statement about the stratum 0.8.0 defers.
  ⭐⭐ **PROVEN 2026-08-13, `g00 ss0.99 capture_off`, true relative abundances as weights with the warm
  start zeroed:** transcript Σ|err| **4,275,988 → 696,900 (0.16×)**, gene **445,944 → 57,154 (0.13×)**,
  gene false-positive mass **23,529 → 0**. A deliberately FLIPPED allocation is **1.26×** worse. So the
  machinery converts good weights into good answers and a weighting function is worth designing.
  ⛔⛔ **THREE LIMITS TRAVEL WITH THAT NUMBER.** It is a CAPABILITY proof, not a ceiling — it does not
  price what a reachable weighting function could earn. True weights also hand over the true SUPPORT (a
  zero weight is EXACTLY absorbing, measured 0.0000 in all four mode × warm-start combinations), and
  4,579 of 8,750 annotated transcripts are silent at `g00`, so most of the fp_mass collapse is free
  rather than earned. And it is ONE condition on the stratum the tool already handles — the blind
  stratum is untested. ⚠ **Under the 0.8.0 scope the third limit is the mildest of the three**: the
  condition it was proven on is IN SCOPE and the untested stratum is the deferred one. The first two
  limits are untouched by any scope ruling.
  ⭐ **What is built:** `splice_graph.build_transcript_path` (the transcript → ordered
  REGION/BOUNDARY/SPLICE-JUNCTION walk, 0/15,669 disagreements against the verified incidence, 45,609
  splice steps, transcription order on 7,312/7,312 minus-strand); the per-transcript weight LANE into
  the EM; `EMConfig.warm_start` = `coverage`/`prior`/`uniform`; the three-way composition published per
  object; `scripts/design/transcript_truth.py` (true per-transcript counts split by splicedness, 41 s
  per 10 M-fragment condition, seven named gates all passing).
  ⛔⛔ **TWO CANDIDATE WEIGHTING FUNCTIONS HAVE BEEN BUILT AND REFUSED — §4 carries both rows and
  `TRAPS: an-upper-bound-is-not-an-estimate` the mechanism. Read them before designing a third.** ⭐ The
  owner's theorem is NOT what failed: the dial is monotone in its favour on every axis and both strata,
  so the pooled no-deconvolution control is the worst rung. What failed is that a bound shared with a
  loud neighbour is not evidence about the quiet one, and 75.3 % of silent transcripts have such a
  neighbour.
  ⭐⭐ **WHAT THOSE EXPERIMENTS SETTLED, AND IT NARROWS THE NEXT ATTEMPT SHARPLY:**
  * **The target is TOTAL abundance, not the unspliced count** — `oracle_alloc` (total truth weights)
     scores **0.163×** against `oracle_alloc_unspliced`'s **0.202×** at transcript level, even though the
     budget being split is literally unspliced pseudocounts. ⚠ At GENE level they invert (0.128 vs
     0.120), so this is a transcript-level ruling only.
  * **The density's units are settled and gated**: `density(r) = mass_rna_region[r]/ceff(r) = Σ_{t∋r} A_t/E_t`
     — the object's own opportunity CANCELS, which is what makes densities at objects of different sizes
     commensurate and what the `min` is a min OF. A transcript alone on its path recovers its own
     abundance exactly, on all four modes (`test_transcript_weights.py`).
  * **The support problem is the whole problem.** 0.0 % of expressed transcripts are ever zeroed; the
     error is entirely that 3,644 silent ones are not. Any next candidate is a SPARSITY mechanism, not a
     better mean.
  ⚠ The lane, the arms and the falsifications are all in place, so the next candidate is one function:
  `transcript_weights.build_weights` + one `alloc_*` arm.
  ⚠ **Identifiability, measured, because it was the objection:** 70.43 % of fragments are compatible
  with 2+ annotated transcripts, but only **0.29 %** of true mass lies in likelihood-flat directions
  (0.96 % by an independent null-space measure) and a census over 952 loci found **zero**
  exactly-exchangeable groups. An oracle per-transcript prior is not meaningfully fake. ⛔ The real
  difficulty is the middle tier: **23.99 %** of true mass sits on 4,277 transcripts with no fragment of
  their own.
  ⚠ **The shipped warm start is NOT the cause of the tool's error** — a UNIFORM start scores 0.976–1.024×
  against base on four conditions, at the noise floor.
- ⏳ **OPEN · RANK 6** — ⛔⛔⛔ **ψ's COMPOSITION DOES NOT CLOSE — a defect, and its own work path** (found 2026-08-12 while
  publishing the three-way composition; owner ruled it a separate branch, taken up after the
  per-transcript prior lands). `NodeDeconv` asserts *"posterior means; f_pos+f_neg+gdna_frac = 1"*.
  Measured on `g00 ss0.99 capture_off`, every object addressed by a chain slot:

  | axis | sums to 1 | sums into (0,1) | sums > 1 |
  |---|---|---|---|
  | REGION | 74.72 % | **25.25 %** — median 0.978, p5 **0.869** | 12 |
  | BOUNDARY | 77.24 % | **22.71 %** — median 0.979, p5 **0.850** | 16 |

  ⛔⛔ **THE MECHANISM RECORDED HERE WAS WRONG IN BOTH HALVES, AND IT IS NOW MEASURED (2026-08-17).**
  It read *"the three posterior means are `np.clip`-ed INDEPENDENTLY, and an unsolvable slot keeps an
  init instead"*, and reasoned *"by linearity of expectation three posterior means over ONE lattice
  should close"*. ⭐ **They are not three means.** `f_g` is the posterior **MEDIAN**
  (`_posterior_median_fg`) while `f_pos`/`f_neg` are posterior **MEANS** of `1 − f_g`, so

      SUM = 1 + median(f_g) − mean(f_g)     — the closure error IS the posterior's SKEW

  verified to **5.8e-15** on both ψ paths. ⚠ The clip never fires (the solver returns `[0, 0.987]` over
  48,000 slots × 3 components) and unsolvable slots close to **exactly 1.000000** — so neither recorded
  cause is real. ⚠ The asymmetry was never a decision: the median for `f_g` was argued for; the RNA
  fractions fall out as expectations of the grid quantity `1 − f_g`.

  ⭐ **TWO CONTRIBUTING TERMS, and the second only looks big on a toy.** (a) genuine posterior SKEW —
  plateaus at ~0.026 under grid refinement, shrinks with depth, **~95 % of it on real data**; (b) GRID
  QUANTISATION of the median — exactly the half-step (0.0423 at `n_grid` 60, 0.0098 at
  `n_grid_ss` 256), constant at every depth, 0.6–8 % on real data. ⛔ A float32 attribution for (b) was
  measured and **refused**: forcing the AMBIG cube to float64 leaves `|SUM−1|` at 0.042272 vs 0.042275.

  ⭐⭐ **SEVERITY: DORMANT, AND IT IS A LANDMINE IN FRONT OF RANK 1.** It reaches NEITHER the published
  masses (`f_g*count` / `(1-f_g)*count`, closed by construction), NOR the EM prior, NOR `derive` —
  `priors.py`/`derive.py`/`track.py` hold zero references to the strand split. `region_init`'s
  `rho_pos`/`rho_neg` carry it (`Σρ_c·E_c = M·SUM`, off by up to 9.5 % on a toy) and their only
  production consumer is `messages/relay.py`, with `message_propagation` **OFF**. ⛔ A claimed
  "+7.0 % library `f_gdna` inflation on unstranded data" was **REFUTED** — measured −0.0007…−0.0031,
  opposite sign, 11–52× smaller. ⚠ `sweep.py`'s own note names the trigger: *"a per-transcript prior
  reading composition per object"* — which is §1's top item.

  ⛔⛔ **AND THE OBVIOUS FIX IS MEASURED CATASTROPHIC.** `f_g` as the posterior MEAN closes the simplex
  exactly by linearity, and `vertex_ceiling.py --arm psi_mean` on all 16 conditions scores
  **1.352 / 1.573 / 3.756** on the three in-scope strata and **1.801** on the `g00` control. The median
  beats the mean at both vertices in **16/16** cases against a proper continuous quantile, so this is a
  real property and not the grid rounding (shipped-vs-proper bias −0.0017, i.e. noise).

  ⭐⭐⭐ **THE REPAIR IS TO ESTIMATE THE PARAMETERS AND MAP, WHICH MAKES CLOSURE STRUCTURAL.** The
  composition has TWO degrees of freedom, not three: `λ → f_g`, RNA total `:= 1 − f_g` (the
  parametrisation, not a second estimate), `θ →` the tilt SHARE splitting it (`simplex_logodds._compose`).
  Measured: closure exact to **2.22e-16 (1 ulp)** on both ψ paths at every κ, depth and strand split.

  ⛔⛔ **AND IT TOOK A SECOND CHANGE, WHICH THE FIRST ONE'S OWN FALLOUT FORCED.** Deriving the RNA total
  from `f_g` propagates whatever error `f_g` carries — and `_posterior_median_fg` SNAPPED to a lattice
  point. On an evidence-free object the posterior IS the reference, symmetric, mean exactly ½; the old
  read-out got `f_pos = 1 − mean = 0.5` exactly while the new one inherited the snap
  (`test_relay_mass_rescale`'s `R_own` read **0.51256**). So the quantile became CONTINUOUS — a histogram on
  the λ lattice, interpolated in the crossing bin.
  ⚠ **On λ, not on f_g**: σ's grid is highly non-uniform, and interpolating there biases a concentrated
  posterior toward ½ by **2.71e-03** at `n_grid` 60. On λ the same posterior returns its own grid point to
  **2.2e-16**, and σ being monotone makes the mapped quantile exact — median equivariance, which is the
  median's whole justification. Measured: `|median(1−f) − (1 − median(f))|` = **3.3e-15**, against
  **8.45e-02** for the snapped version.

  ⛔ **`f_g` IS THEREFORE *NOT* BIT-IDENTICAL** — an earlier draft of this row said it was, which was true
  of the composition change ALONE and false once the quantile moved: `f_g` shifts by up to the half-grid
  step it used to snap by (**9.8e-03** single-strand, **4.2e-02** AMBIG). ⚠ The panel must therefore be
  re-measured rather than argued from, and it was.
  ⛔ **This is NOT the renormalisation refused below**: nothing is rescaled to hide a residual — `f_g`
  and the RNA total are exact complements by construction, and a tilt is estimated as a share because a
  share is what it IS.

  ⛔⛔ **AND CLOSURE IS *NOT* THE VERTEX PROBLEM IN DISGUISE — HYPOTHESIS RAISED AND REFUTED
  (2026-08-17), WHICH SAVES AN EXPENSIVE EXPERIMENT.** The reading was that both are posterior
  ASYMMETRY, so fixing the vertex by structural certainty would collapse skew, closure and vertex bias
  together — and it would have been tested by rescoping `struct_lock`, i.e. by walking into
  `TRAPS: a-cancelling-defect-pair`. Measured against the origin-split oracle on three conditions, per
  band of the TRUTH's distance to a vertex, mass-weighted `|SUM − 1|`:

      truth band          g50 ss0.99   g50 ss0.50   g98 ss0.99 ON
      vertex   <0.02        0.0021       0.0074        0.0064
      centre   >0.30        0.0108       0.0282        0.0613      ⇐ 4-10x LARGER

  `corr(|closure|, closeness-to-vertex)` = **−0.069 / −0.149 / −0.080** — negative on all three. Closure
  fails **at the CENTRE of the simplex, not at the vertices**: a posterior pinned at a vertex is NARROW,
  so median ≈ mean; a low-information slot mid-composition is WIDE and skewed. ⭐ The driver is posterior
  WIDTH, not vertex repulsion. ⚠ The median's advantage over the mean is also not vertex-concentrated —
  `corr` with it is **−0.104 / −0.491 / +0.262**, inconsistent in sign, and by mass the median wins only
  46.4 % of the vertex band on `g50 ss0.99`. (That does not contradict the `psi_mean` panel arm: the arm
  changes the whole SOLVE — belief, refit, landscape — while this compares readouts on one solve.)
  ⭐⭐ **So closure needs its own fix and will not come free from the vertex work** — which is exactly why
  the parametrise-and-map repair above is worth having: it is structural, and independent of both.
  ⚠ The only population that closes exactly is the one that never reaches ψ: measured, **1,243 of 21,522
  slots close, and all 1,243 are intergenic** (0 introns, 0 exons) — the unsolvable slots that keep their
  simplex-valued init.

  ⚠⚠ **A SEPARATE, PRE-EXISTING DEFECT FOUND WHILE HUNTING BUGS IN THE ABOVE — RECORDED, NOT FIXED.**
  At `κ = ½` the strand mean is `p = ½·(f_g+f_pos+f_neg) = ½` identically, so the RNA tilt is PROVABLY
  UNIDENTIFIED and `w_pos` must be exactly ½. The AMBIG cube evaluates that sum in **float32**, where it
  departs from 1 by ~1e-7 in a τ-DEPENDENT way; times the slot count, that is a spurious log-likelihood
  tilt growing linearly in depth. Measured `|w − ½|`: **7e-6** at 1e4 fragments, 3.4e-4 at 1e6,
  **5.3e-3** at 1e7, **0.199** at 1e8, **0.425** at 1e9 (a 92.5/7.5 split where the truth is 50/50).
  Falsified decisively: recompute the same cube with the strand term in float64 — everything else
  unchanged, cube still stored f32 — and `w_pos → 0.500000000` at every depth. ⭐ κ = 0.9/0.99 are
  unaffected; it is specific to the degenerate-κ ridge, i.e. the `ss 0.50` half of the ladder.
  ⭐ **NEGLIGIBLE AT PANEL SCALE** — 10 M fragments over ~70 k slots is ~140 per slot, largest ~1e5–1e6,
  so ≤ 3.4e-4 — and **entirely pre-existing**: the retired read-out published `m_pos` directly, the
  parametrised one publishes `(1−f_g)·m_pos/(m_pos+m_neg)`, the SAME ratio with only the total rescaled.
  ⛔ It would matter on a very deep library with a dominant locus, and the repair is one dtype on the
  strand term, not a redesign.

  ⛔ **Do NOT repair it by renormalising all three at publication** — that makes a 15 %-short object
  indistinguishable from a solved one. `CalibrationResult` publishes the values as solved, bounds each
  to `[0,1]`, and asserts no closure; `test_a_composition_that_does_NOT_close_is_ACCEPTED_and_that_is_deliberate`
  pins that as a decision so the assertion is not added before ψ is fixed.
- ⏳ **OPEN · RANK 7 — Price the CANCELLING PAIR together** — `struct_lock` rescoped to
  `g1_locked & REGION` AND the
  `intergenic|exon` boundary claiming its RNA-contaminated crossing mass as gDNA. Five of the seven
  xfails interrogated on 2026-08-11 go green if and only if this lands. ⛔ Neither half may be priced
  alone: the `struct_lock` half is **+3,207 %** on the zero-gDNA control by itself.
  ⭐ **It is in scope**: its target is the false-positive gDNA channel, which the `g00` control measures
  on every stratum including the three 0.8.0 ships.
  ⛔⛔ **AND IT CAN ONLY BE PRICED WITH `--messages on`.** `zc_struct_lock_g1` is MEASURED inert with
  them off (2026-08-11): it rewrites the mask on 30,423 slots and changes no scored field, because
  `struct_lock` reaches the answer only through `own_composition_logvar` → a precision → a fusion that
  `SilentPolicy` never performs. Priced against a messages-off base it would read as "worth nothing",
  which is the exact reading `ladder_arm_ab.py`'s `--messages` stamp and up-front refusal exist to make
  impossible (`TRAPS: an-ablation-that-never-ran`).
- ⏳ **OPEN · RANK 8 — Re-price message propagation, and READ IT AGAINST THE 0.8.0 SCOPE.** Its recorded
  +155 % predates the conserved-count rewrite AND the numeric convention, and both of its recorded
  halves have moved scope: the −58 / −44 / −32 % is on the **three IN-SCOPE strata** and the +155 % is on
  the **DEFERRED** one. ⛔ So the question 0.8.0 asks of it is not "does it rescue the blind stratum" —
  that stratum is deferred — but **what it costs the ZERO controls**, where the recorded price is
  0.029 → 89.93 and 0.005 → 9.58. ⛔ RE-price, never inherit; the owner's 2026-08-10 ruling (off until
  the tool is optimised end to end) is not re-opened by a measurement.
  ⭐ It is one pair of commands — `ladder_arm_ab.py --arm base --messages off` against
  `--arm base --messages on`, `--jobs 8` — with the setting stamped into every row of both.
  ⚠ It also closes two of the xfails, which is a consequence and not a reason.
- ⏳ **OPEN · RANK 9 — CORRECTNESS AND HYGIENE.** ⛔ Each is its own commit (`TRAPS: one-thing-varied`),
  and none of them moves the 0.8.0 metric.
  * ⭐ **THE INDEX SHOULD EMIT ITS DUPLICATE MAP** (owner, 2026-08-13). The index drops EXACT duplicates
    — same ref, strand and sorted exon tuple — and is right to: they are *"mathematically unidentifiable
    from any read data"* (`index.py`'s own words). ⛔ But it records NOTHING about what it dropped, so a
    dropped id is unfindable through the index and every consumer re-derives the map by re-parsing the
    GTF — which defeats the index being freestanding. Measured on the suite reference: 8,755 GTF
    transcripts, 5 duplicate groups over 10 transcripts, 5 dropped → the 8,750 indexed; 2 of the 5 carry
    simulated fragments (329 at `g00`).
    ⭐ **The fix is an ALIAS MAP, not re-admission** — `dropped_t_id → kept_t_id` in the manifest, exposed
    as `TranscriptIndex.duplicate_map` defaulting to `{}` so an older index still loads. Re-admitting the
    transcripts would undo a collapse that is correct.
    ⭐⭐ **It needs an index rebuild but NOT a panel re-scan**, and that is worth stating because the two
    were conflated once: the scan cache is keyed on `graph_hash` / `reach_digest` /
    `payload_schema_digest` / `scan_config_digest`, and a metadata field on the same partition changes
    none of them. ⚠ A rebuild can still move `reach_digest` on its own (`reach` is covered by no other
    hash, and a rebuild once moved 38 % of human reaches), so verify with `rescan_panels.py` rather than
    assuming the caches survive. ⚠ Then DELETE `transcript_truth.py`'s local `duplicate_map` and read the
    index's — moving it, not copying it, so there is one home.
  * **`mass_*_boundary` → `count_*_boundary`.** `mass_gdna_boundary` and `mass_rna_boundary` on
    `CalibrationResult` are crossing INCIDENCES, not masses; the third one already landed as
    `count_rna_sj`.
  * **Restore the moment tests.** `contained_moments` / `crossing_moments` / `build_slot_moments` moved
    to `effective_length.py` without their tests, which went with the deleted length channel — untested
    geometry in a live layer-2 module. ⚠ This is TEST debt in a live geometry module and is **not** the
    length channel: retiring that channel until after 0.8.0 does not retire these.
  * **The cheap ledger:** dead substrate surface (`PopulationView.length_sum` / `.mean_length` /
    `.total_inv_length_sum` have no consumer in `src/`); `second_pass.py:443-445`'s comment, measured
    false in premise and harmless in effect; `module_census.py` re-run, since the purge left stale
    sibling references; and `UNDOCUMENTED_DEBT` — 8 scripts on disk and unlisted, in
    `tests/test_scripts_index.py`.
- ⏸ **OPEN · RANK 10 — PROFILE, AND THE PERFORMANCE SUBSTRATE IS A TRAP.** The one characterisation item
  never started; new code means slowdown and it has been a while. Measured 2026-08-07 on one 10 M-fragment
  ladder condition (a 35,135-REGION chr22 index): per-locus EM **15.9 s (47 %)**, native scan 6.5 s,
  calibration 6.5 s, second-pass drain 3.5 s, total 33.5 s. ⚠ Message propagation costs nothing measurable
  (33.7 s ON vs 33.5 s OFF) — "the relay is the expensive part" is REFUTED.
  ⛔⛔ **BUT THAT PROFILE IS UPSIDE DOWN FROM THE REAL ONE AND BOTH ARE CORRECT.** Calibration is
  depth-INDEPENDENT — every region in the index is solved regardless of read depth — so it scales with the
  INDEX while the EM scales with the DATA. On real cfRNA at genome scale (~1.5 M regions) calibration has
  measured ~66 s against the EM's ~24 s: the exact reverse. **Profile on real cfRNA**
  (`~/Downloads/rigel_runs/cfrna/`, genome index at `~/Downloads/rigel_runs/refs/rigel_index`), never on
  this panel — `TRAPS: toys-rank-hotspots-backwards`, which cost a whole analysis once.

⚠ **Tracked, not blocking:** two exact per-fragment mirrors of a PURE-RNA library give different
deconvolutions — `mass_gdna_region` +5.4 % at `ss=0.75`, +2.0 % at `ss=1.0`, 0 at `ss=0.90`. Not
boundary-only and not monotone. It is 5.4 % of the zero-gDNA false-positive channel, and an R1-sense
library is now simulable, so it is measurable.

### §1.1 THE BRIEF FOR **the-cancelling-pair** — everything that experiment needs, in one place

⭐ RANK 6's brief, and the only brief in this file. Kept whole because the analysis is complete and
re-deriving it would be waste.

**THE TWO HALVES, AND WHY NEITHER HAS AN HONEST PRICE ALONE.** The certified-RNA channel is a LOWER BOUND
(`ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)` — the exon may hold molecules that never touch that boundary) and ψ delivers it
as a two-sided Gaussian, penalising a destination holding MORE RNA exactly as hard as one holding less.
Making it one-sided — `−½·p·max(0, mo − log f)²`, no new constant — is **the only mechanism the zero-gDNA
control has ever ENDORSED: −81.9 %, 8/8.** The panel is nonetheless **+5.2 %**, and the reason is
`TRAPS: a-cancelling-defect-pair`: the two-sided term was doing TWO jobs, the bound *and* — by its upward
side — a de-facto gDNA LEVEL channel. Removing the accident exposes the real gap, +7.7 % on exactly the
stratum with no working level channel and −9.2 % where the prior already suffices.

⛔⛔ **AND THE LEVEL CHANNEL IS STRUCTURALLY DISCONNECTED, WHICH IS A THEOREM AND NOT A MEASUREMENT.**
Receivers of a one-slot-step channel at `g00`, by slot type:

| channel | BOUNDARY | intergenic REGION | intron REGION | exon REGION |
|---|---|---|---|---|
| gDNA level | 2,592 | 0 | 0 | **0** |
| certified RNA | **0** | 0 | 0 | 10,493 |

The chain is strictly REGION · BOUNDARY · REGION · BOUNDARY · …, so a one-slot-step channel is **BIPARTITE**: a REGION emitter can only
reach a BOUNDARY. The only licensed originators of a gDNA level are structurally pure-gDNA objects, which are
REGIONs — so **no REGION can ever receive a gDNA level.** It is not weak, it is disconnected from every object
carrying mass. ⛔ Two repairs are already refused: letting a BOUNDARY originate a level (a patch on a
symptom), and letting the level cross two steps (the chain-fused level is dominated by the 1,312 intergenic
anchors and hands every exon the off-probe floor, against a true neighbouring density **346×** higher).

⚠ **The measurement to run it against**: at `g01 ss0.50 capture_on` — ⛔ a rung the 2026-08-13 rebuild DROPPED; the nearest surviving condition is `g05 ss0.50 capture_on` — ψ-boundary ablation with the identity
exact, every single channel ablation is small (+2.3 / +6.8 / +4.9 / −0.1 %) and the joint one is +60.7 % —
`TRAPS: all-small-singly-large-jointly`, so the four share an upstream quantity. And at HEAD's top-12 error
slots the self-solve WITH the fitted prior is nearly correct while the messages destroy it; muting the
certified-RNA channel alone recovers the truth at 8 of the 12.

## §2 THE OPEN QUESTIONS, ranked by how much they would change a decision

⚠ **Ranked as questions, not as work** — §1 is the ranked work list. ⛔ A question about the DEFERRED
stratum or a deleted panel cannot change a 0.8.0 decision until something else does.

1. ⭐⭐⭐ **WHY IS PRIOR FIDELITY ANTI-CORRELATED WITH DELIVERABLE QUALITY? — LIKELY ANSWERED AS "IT IS NOT
   THE PRIOR, IT IS THE MESSAGES", AND THERE IS NOW A SECOND CANDIDATE ANSWER TO EXCLUDE.** At every one
   of HEAD's 12 worst slots the self-solve WITH the fitted prior is nearly correct (0.0043–0.0344 against
   truths of 0.0007–0.0136) and the message layer then destroys it (0.05–0.83); muting the certified-RNA
   channel alone recovers the truth at 8 of 12 (2026-08-06, ψ-boundary ablation at `g01 ss0.50 capture_on`
   — ⛔ a rung dropped 2026-08-13, nearest surviving `g05`). ⚠ **Confirm on a second stratum before
   closing** — and ⛔ **2026-08-14 adds the alternative that must be excluded first: every prior-injection
   arm was measured with the effective-length shrinkage left unsubstituted (§1, rank 1)**, so part of "a
   better prior does not show up" may be the RULER rather than the messages.
   `TRAPS: the-intermediate-is-not-the-deliverable`.
2. ⭐⭐ **What is the joint price of the splice-flux reframe and the level defect?** **the-cancelling-pair**. Neither has an
   honest price alone, and one of them is already in `src/`.
3. ⭐ **Do we solve ALTERNATIVE SPLICING correctly?** The `alt_splice` toy rung exists and is
   **unverified** (`toy_harness.py --list` carries it). Cheap, and it is the only structure where several
   sj share a BOUNDARY — so it is in scope on all three shipping strata.
4. ⚠ **Does the reframe's own `σ²_transfer` correctly price a ratio built on 19 counts?** `EQUATIONS.md`
   §3.5d says it is the right medicine for the noise and the wrong medicine for the share-weighting.
   Unmeasured.
5. ⚠ **The two capture-ON pilot rows disagreed about the SIGN of every length correction**, across two
   independently built panels. Unexplained, and ⛔ both panels were DELETED 2026-08-13, so answering it
   requires rebuilding them. Do not average them; find which one is lying.
6. ⚠ **Should the simulator capture a pre-mRNA through EVERY probe whose genomic blocks it spans, not
   only its own transcript's?** (2026-08-19, `TRAPS: the-panel-enriches-nascent-by-its-own-probes`.)
   Physically yes — every library molecule is ds-cDNA at capture — and the panel currently under-enriches
   nascent by up to 6× at shared or foreign probes, which puts a 5–37 % composition residual on the
   `nrna_mid × capture-ON` intron→boundary hops that is not a rule failure. It also touches the solver's
   intron deconvolution under capture with nascent. An OWNER call on the simulator; the map's verdicts
   do not depend on it.

---

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
