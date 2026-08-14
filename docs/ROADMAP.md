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
`flgap_study_cache.py`'s: keyed on the scan manifest plus a content hash of every producing source file,
with the module under test deliberately OUTSIDE the key so its edit loop is one second.

## §0 THE STATE

⛔ **Re-derive rather than trust.** A number that has moved is a result, not a documentation bug.

⚠ **Twenty-one rows, re-counted 2026-08-14 rather than adjusted (`TRAPS: re-record-the-baseline`). If it
grows again, PRUNE rather than let it become a report** — this table says where the tool IS, so a ✅ row
earns its place only by saying *do not work here*, and a number nobody can attribute earns none.
⛔ **Read the stratum, never the total.** Every row that names unstranded × capture-ON is describing the
**DEFERRED** stratum, and it is reported rather than ranked.

| | | |
|---|---|---|
| **Stage A — the accumulator** | ✅ **DONE**, and that is a measurement | perfecting BOTH fragment-length models is worth **2.6 %** of the deliverable, down from 22.2 % |
| ✅ **the tally CONSERVES a fragment count** | ⭐ every deposited fragment places **exactly one** unit across the objects it crosses — regions, contiguous boundaries and sj boundaries together. Measured on the origin-split oracle: **1.000× deposited, 0 unaccounted**, on BOTH origins | ⛔ RNA read **0.747×** until `sj_mass` landed (2026-08-11): a spliced fragment whose blocks cross no boundary deposited on no conserved bank at all — **1,222,375 of 4,830,713 (25.3 %)**. gDNA was always exact, because gDNA cannot splice |
| ✅ **the SUITE is green with no skips** | ⭐ **0 failed / 3,404 passed / 0 skipped / 10 xfail** (re-derived 2026-08-14); it was **3,397** on 2026-08-13 and **3,235 / 0 / 7** on 2026-08-11, and every delta is ACCOUNTED file by file in `CLAUDE.md`, never adjusted. Two xfails and both skips were closed by BUILDING the missing thing, never by widening a bound | ⛔ The 7 interrogated on 2026-08-11 are not one kind of thing: 2 are the recorded price of `message_propagation = False`, 5 were WRITTEN as xfail as records of two proven defects whose fixes are panel-negative alone. `CLAUDE.md` carries that split; ⚠ the three added since have not been interrogated the same way |
| ✅ **ONE NUMERIC CONVENTION** | ⭐ a COUNT is an integer, a FRACTION is float64. No fixed point, no scale constant, nothing decodes a bank (owner, 2026-08-11) | ⭐ float64 is **1e5–7e5× MORE accurate** than the fixed point it replaced, measured against exact rational arithmetic. ⛔ The ~2.6 % once quoted against float is a **float32** number — `TRAPS: integer-channels-reproduce` |
| ⭐ **calibration, the THREE IN-SCOPE strata** | ✅ median library `f_gdna` error **0.005–0.012**, and ⭐ the PRIOR the EM reads is within **2.5–4.6 %** of a perfect one | stranded × capture-ON/OFF and unstranded × capture-OFF — **the 0.8.0 target**. ⚠ 36-condition ladder, RETIRED 2026-08-13 |
| ⛔ **calibration, unstranded × capture-ON — the DEFERRED stratum** | ⛔ **BLIND** — reports **0.033–0.058** while truth spans **0.00 → 0.98**, and hands the EM a gDNA prior **94.4 %** short. ⭐ Re-measured 2026-08-13/14 on the rebuilt ladder at the exon objects: `f_g` **0.040 / 0.0016 / 0.0021** at `g05 / g50 / g98` against truths **0.054 / 0.518 / 0.982** — a NEAR-ZERO answer regardless of truth, which **looks acceptable at low gDNA by coincidence** | not noisy; a flat boundary. ⛔⛔ **DEFERRED for 0.8.0** (owner, 2026-08-14) — reported on every panel, not a development target until the other three are optimised. The length channel was tried and refused, and is retired until after 0.8.0 (the scope block above) |
| ⛔⛔ **the EFFECTIVE-LENGTH SHRINKAGE at the ZERO-gDNA control** | ⛔ **NEW 2026-08-14.** The shipped shrinkage contracts every transcript by a mean factor **0.345** when the correct answer is exactly **1.000**, on **13,673 of 15,669** transcripts — **INCLUDING on capture-OFF, where the module contract says the factor must be 1**. `rho_ref` is fabricated entirely from false-positive gDNA | ⭐ **It is a SYMPTOM of the composition error, not an independent bug**: substituting ONLY the composition arrays and re-running the SHIPPED shrinkage function gives the correct factor (`g00` capture-off **0.345 → 1.000**; `g50` unstranded capture-on **0.834 → 0.401**, the truth). ⛔ So repair the composition it reads, not the function. One split, two consumers, one function — `priors.py` imports `_global_reference_density` **from** `capture_eff_length.py` |
| ⛔⛔ **EVERY CEILING WAS MEASURED WITH A WRONG RULER INSTALLED** | ⛔ **NEW 2026-08-14.** `effective_lengths_em` is built by `_setup_geometry_and_estimator` **BEFORE** `assemble_priors` runs, and **every existing measurement arm patches `assemble_priors`** — so the shrinkage has **NEVER been priced by any ceiling** | ⛔ This invalidates the PRECISION of every prior-injection ceiling in this file, not their sign. It is why §1's rank 1 is an instrument and not a mechanism |
| ⛔ **the PER-TRANSCRIPT prior lane is BUILT and NEVER PASSED** | ⛔ **NEW 2026-08-14.** `rna_prior_weight` is built end to end and the production call site in `pipeline.py` **omits it**, so the EM's default rule carries **zero per-transcript information** | ⚠ Passing it needs a weighting function, and the two built so far were REFUSED (§4). The lane, its arms and its falsifications are all in place |
| ⭐ **what a PERFECT prior is worth, by GRANULARITY** | ⭐ **NEW 2026-08-14.** Per-**LOCUS**: gene-level error on the deferred stratum **0.099**. Per-**TRANSCRIPT**: the three in-scope strata to **0.37–0.58** gene-level, with false-positive mass cut **100–200×** | ⛔ Per-transcript does **NOT** rescue unstranded × capture-ON (**0.641**) — which is the deferred stratum, so that is a reported fact and not a blocker. ⭐ The in-scope prize is the per-transcript one |
| ✅ **the prior ASSEMBLER and its POPULATION** | ✅ `rel` **0.0019–0.0027** with perfect masses in and **4.9e-4** with perfect per-component shares; the composition claim `a_g:a_r` is exact against the unspliced pool (`Δphi` **≤ 5e-4**) | ⭐ the assembler was **0.179**: the conserved-count rewrite took it to 0.0202 and the YARDSTICK took it to 0.0027 (`TRAPS: score-the-consumers-own-count`) |
| ✅ **the fragment-length models** | ✅ accurate — `pi(w)`'s de-tilt reads **211.77** against a true **212.20**, and the gDNA pmf is exact to **0.02 bp** off capture | ⛔ a claimed +10.7 % bias was a truth-parser bug (`TRAPS: a-truth-table-of-aggregates`) |
| ⚠ **the gDNA pmf under capture** | ⛔ the last unfixed length defect: **+13.6 bp** at a 330 bp gDNA mean and **+3.5 bp** at 120, drained per-bin `fit/true` **1.22 … 4.18** in the tail | ⭐ EXACT off capture and untouched by the drain, so it is a PLACEMENT problem — `EQUATIONS.md` §4.4 — and must not be attacked by editing `gdna_opportunity` |
| **message propagation** | ⛔ **OFF**, and it stays off until the tool is optimised end to end across all scenarios (owner, 2026-08-10) | net better on **all three IN-SCOPE strata** (−58 / −44 / −32 %, 16/16 on the two stranded ones); **+155 % on unstranded × capture-ON**, which is the **DEFERRED** stratum. ⚠ **The 0.8.0 scope changes how that reads and does NOT re-open the ruling**: the whole recorded price lands on the stratum 0.8.0 does not target, so the honest question is now what it costs the ZERO controls. ⛔ RE-price on the rebuilt ladder, never inherit — every number in this row predates the conserved-count rewrite and the numeric convention |
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
| **1** | ⭐⭐⭐ **THE RULER FIRST — get the effective-length shrinkage inside a measurement arm** | `effective_lengths_em` is built BEFORE `assemble_priors`, and every arm patches `assemble_priors`, so **no ceiling has ever priced it**. Until an arm reaches it, every ranking below is provisional — including this one |
| **2** | ⭐⭐⭐ **THE SHRINKAGE AT `g00`, REPAIRED UPSTREAM** | 0.345 where the truth is exactly 1.000, on 13,673 of 15,669 transcripts, **including capture-OFF where the contract says 1** — and it is a SYMPTOM: feed the shipped function correct composition arrays and it returns the correct factor. ⛔ Repair the composition, not the function |
| **3** | ⭐⭐⭐ **PASS THE PER-TRANSCRIPT PRIOR LANE — and the weighting function is the work** | built end to end, omitted at the production call site, so the EM carries **zero** per-transcript information. The prize is the largest in-scope one measured: **0.37–0.58** gene-level on the three strata and false-positive mass **100–200×** down |
| **4** | ⭐⭐ **RE-DERIVE THE IN-SCOPE NUMBERS ON THE 16-CONDITION LADDER** | almost every number in §0 was measured on the ladder retired 2026-08-13. Noise floor **0.996–1.013** — an arm inside that band moved nothing |
| **5** | ⭐⭐ **ψ's COMPOSITION DOES NOT CLOSE** | a real defect on every stratum: 25.25 % of REGIONs and 22.71 % of BOUNDARYs sum into (0,1), p5 **0.869** / **0.850**. Owner ruled it its own branch, taken after the per-transcript prior |
| **6** | ⭐ **PRICE THE CANCELLING PAIR TOGETHER** | five xfails go green iff it lands, and neither half has an honest price alone. ⛔ It can only be priced with `--messages on` |
| **7** | ⚠ **RE-PRICE MESSAGE PROPAGATION** | ⛔ **its recorded price falls ENTIRELY on the deferred stratum**, so the in-scope question is what it costs the ZERO controls. The owner's 2026-08-10 ruling stands; this is a measurement, not a flip |
| **8** | the index's duplicate map · `mass_*_boundary` → `count_*_boundary` · restore the moment tests · the cheap ledger | correctness and hygiene, each its own commit, `TRAPS: one-thing-varied` |
| **9** | ⏸ **PROFILE — the one characterisation item never started** | new code means slowdown and it has been a while. ⛔ It does not move the 0.8.0 metric, and it has a trap of its own: profile on real cfRNA, never on this panel |

⛔ **NOT ON THIS LIST, DELIBERATELY:** the length channel (retired until after 0.8.0 — the scope block
above) and anything whose only target is unstranded × capture-ON (deferred). §4 is everything else
already priced and refused; read it before proposing a mechanism.

⛔ **The blocks below are the DETAIL, one per rank, in rank order.** A rank with no block is fully stated
by its table row.

⛔⛔ **FOUR INSTRUMENTS ARE DARK AND THE `CLAUDE.md` SCRIPTS TABLE STILL LISTS THEM AS RUNNABLE.** On
2026-08-13 `pilot`, `flgap_short` and `flgap_long` were deleted and the ladder rebuilt at 16 conditions, so
anything reading those panels no longer runs at all: **`flgap_study_cache.py`, `prior_yardstick.py`,
`boundary_q_population.py` and `prior_vs_oracle.py`'s flgap arms.** ⚠ `prior_yardstick.py` is the only
place the prior assembler is scored the way production runs it, and that — not the length channel — is
the reason rebuilding the flgap PAIR might be worth paying for (they only work as a pair). Every other
instrument under `scripts/design/` runs.

- ⏳ **OPEN · RANK 1 — THE RULER: THE EFFECTIVE-LENGTH SHRINKAGE HAS NEVER BEEN INSIDE A MEASUREMENT ARM**
  (found 2026-08-14). `effective_lengths_em` is built by `_setup_geometry_and_estimator`, which runs
  **before** `assemble_priors`; **every** existing measurement arm — every `oracle_*` arm of
  `quant_accuracy.py`, every ceiling — patches `assemble_priors`. So the ruler the EM measures transcripts
  with was never substituted, and **every ceiling number in this file was measured with a wrong ruler
  installed.**
  ⛔ **This ranks FIRST because it ranks everything else.** It does not change the SIGN of any recorded
  ceiling; it changes how much of the residual is attributable to the prior rather than to the length the
  prior is divided by. Until an arm can substitute both, "a perfect prior is worth X" is a statement about
  a composite.
  ⭐ **The build is an arm, not a mechanism**: an override that reaches `effective_lengths_em` at the same
  point the pipeline builds it, plus the `noop` gate the harness family already requires — and ⚠ the noop
  gate must be specified against the `base_reseed` noise floor, not byte-identity, for the reason §0's own
  row records.
- ⏳ **OPEN · RANK 2 — THE SHRINKAGE IS WRONG AT THE ZERO-gDNA CONTROL, AND IT IS A SYMPTOM** (measured
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
  ⚠ Ranks below the ruler only because without rank 1 its price cannot be read off any panel arm.
- ⏳ **OPEN · RANK 3 — THE PER-TRANSCRIPT RNA PRIOR: THE FOUNDATION IS BUILT AND PROVEN; THE WEIGHTING
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
- ⏳ **OPEN · RANK 5** — ⛔⛔⛔ **ψ's COMPOSITION DOES NOT CLOSE — a defect, and its own work path** (found 2026-08-12 while
  publishing the three-way composition; owner ruled it a separate branch, taken up after the
  per-transcript prior lands). `NodeDeconv` asserts *"posterior means; f_pos+f_neg+gdna_frac = 1"*.
  Measured on `g00 ss0.99 capture_off`, every object addressed by a chain slot:

  | axis | sums to 1 | sums into (0,1) | sums > 1 |
  |---|---|---|---|
  | REGION | 74.72 % | **25.25 %** — median 0.978, p5 **0.869** | 12 |
  | BOUNDARY | 77.24 % | **22.71 %** — median 0.979, p5 **0.850** | 16 |

  ⭐ **The mechanism is visible at `sweep.py`'s write-back**: the three posterior means are
  `np.clip(·, 0, 1)`-ed **INDEPENDENTLY**, and an unsolvable slot keeps an init instead — neither
  preserves a simplex. ⚠ By linearity of expectation three posterior means over ONE lattice should
  close, so the deficit is not inherent to taking means and the clip is a symptom rather than the whole
  cause. ⛔ It was invisible because nothing consumed the strand split: the crossing axis published
  `0` for both RNA strands, so its composition summed to `f_g` alone.
  ⛔ **Do NOT repair it by renormalising at publication** — that makes a 15 %-short object
  indistinguishable from a solved one. `CalibrationResult` publishes the values as solved, bounds each
  to `[0,1]`, and asserts no closure; `test_a_composition_that_does_NOT_close_is_ACCEPTED_and_that_is_deliberate`
  pins that as a decision so the assertion is not added before ψ is fixed.
- ⏳ **OPEN · RANK 6 — Price the CANCELLING PAIR together** — `struct_lock` rescoped to
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
- ⏳ **OPEN · RANK 7 — Re-price message propagation, and READ IT AGAINST THE 0.8.0 SCOPE.** Its recorded
  +155 % predates the conserved-count rewrite AND the numeric convention, and both of its recorded
  halves have moved scope: the −58 / −44 / −32 % is on the **three IN-SCOPE strata** and the +155 % is on
  the **DEFERRED** one. ⛔ So the question 0.8.0 asks of it is not "does it rescue the blind stratum" —
  that stratum is deferred — but **what it costs the ZERO controls**, where the recorded price is
  0.029 → 89.93 and 0.005 → 9.58. ⛔ RE-price, never inherit; the owner's 2026-08-10 ruling (off until
  the tool is optimised end to end) is not re-opened by a measurement.
  ⭐ It is one pair of commands — `ladder_arm_ab.py --arm base --messages off` against
  `--arm base --messages on`, `--jobs 8` — with the setting stamped into every row of both.
  ⚠ It also closes two of the xfails, which is a consequence and not a reason.
- ⏳ **OPEN · RANK 8 — CORRECTNESS AND HYGIENE.** ⛔ Each is its own commit (`TRAPS: one-thing-varied`),
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
- ⏸ **OPEN · RANK 9 — PROFILE, AND THE PERFORMANCE SUBSTRATE IS A TRAP.** The one characterisation item
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
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_pin.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
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
in §1 (rank 6) — as **half of a pair**, never alone, exactly as its row says.

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

⭐⭐ **THE PATTERN, AND IT IS THE REAL RESULT: every one of the eleven was a rule for how to resolve DOUBT,
and at `g00` the doubt must resolve to NO gDNA.** A rule that lifts an evidence-free slot off zero is
inadmissible there however well it scores elsewhere. ⛔ The only candidate the control has ever ENDORSED is
the one-sided certified-RNA bound (−81.9 %, 8/8), and that one is panel-negative alone — it is half of
`TRAPS: a-cancelling-defect-pair` and is §1's rank 6, **the-cancelling-pair** (§1.1).

⚠ Four more `zc_*` arms exist as decomposition REVERTS used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.
