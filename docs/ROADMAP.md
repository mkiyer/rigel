# ROADMAP — where the tool is, and what to do next

Reading order for a new session: **`SUCCESS.md`** (how performance is measured) → **this file** (current
state and priority) → **`TRAPS.md`** (mistakes not to repeat) → `DESIGN.md` / `EQUATIONS.md` as needed.

⛔⛔ **THIS FILE POINTS FORWARD. IT IS NOT A CHANGELOG.** A closed item earns **one line** in §3, and only
because that line is what stops it being rebuilt. The derivation goes to `EQUATIONS.md`, the lesson to
`TRAPS.md`, the ruling to `DESIGN.md`, and the investigation stays in git. ⚠ This file reached 1,030 lines
before its first prune (2026-08-05) by accumulating a section per campaign; if a section here reads like a
report on work already done, delete it.

⛔ **Every number below was measured on the current tree and the current panel.** Re-derive rather than
trust — a number that has moved is a result, not a documentation bug.

---

## §0 THE STATE

⛔ **Re-derive rather than trust.** A number that has moved is a result, not a documentation bug.

⚠ **Fourteen rows. If it grows again, PRUNE rather than let it become a report** — a closed item earns
one line, and a number nobody can attribute earns none.

| | | |
|---|---|---|
| **Stage A — the accumulator** | ✅ **DONE**, and that is a measurement | perfecting BOTH fragment-length models is worth **2.6 %** of the deliverable, down from 22.2 % |
| ✅ **the tally CONSERVES a fragment count** | ⭐ every deposited fragment places **exactly one** unit across the objects it crosses — nodes, contiguous edges and junction edges together. Measured on the origin-split oracle: **1.000× deposited, 0 unaccounted**, on BOTH origins | ⛔ RNA read **0.747×** until `sj_mass` landed (2026-08-11): a spliced fragment whose blocks cross no line deposited on no conserved bank at all — **1,222,375 of 4,830,713 (25.3 %)**. gDNA was always exact, because gDNA cannot splice |
| ✅ **the SUITE is green with no skips** | ⭐ **3,235 passed / 0 skipped / 7 xfail** (2026-08-11). Two xfails and both skips were closed by BUILDING the missing thing, never by widening a bound | ⛔ The 7 are not one kind of thing: 2 are the recorded price of `message_propagation = False`, 5 were WRITTEN as xfail as records of two proven defects whose fixes are panel-negative alone. `CLAUDE.md` carries the split |
| ✅ **ONE NUMERIC CONVENTION** | ⭐ a COUNT is an integer, a FRACTION is float64. No fixed point, no scale constant, nothing decodes a bank (owner, 2026-08-11) | ⭐ float64 is **1e5–7e5× MORE accurate** than the fixed point it replaced, measured against exact rational arithmetic. ⛔ The ~2.6 % once quoted against float is a **float32** number — `TRAPS: integer-channels-reproduce` |
| **calibration, 3 of 4 strata** | ✅ median library `f_gdna` error **0.005–0.012**, and ⭐ the PRIOR the EM reads is within **2.5–4.6 %** of a perfect one | stranded × on/off and unstranded × capture-OFF |
| ⛔ **calibration, unstranded × capture-ON** | ⛔ **BLIND** — reports **0.033–0.058** while truth spans **0.00 → 0.98**, and hands the EM a gDNA prior **94.4 %** short | not noisy; a flat line. **Still the whole open problem**, and nothing has fixed it — the length channel was tried and refused (§1.4) |
| ✅ **the prior ASSEMBLER and its POPULATION** | ✅ `rel` **0.0019–0.0027** with perfect masses in and **4.9e-4** with perfect per-component shares; the composition claim `a_g:a_r` is exact against the unspliced pool (`Δphi` **≤ 5e-4**) | ⭐ the assembler was **0.179**: the conserved-count rewrite took it to 0.0202 and the YARDSTICK took it to 0.0027 (`TRAPS: score-the-consumers-own-count`) |
| ✅ **the fragment-length models** | ✅ accurate — `pi(w)`'s de-tilt reads **211.77** against a true **212.20**, and the gDNA pmf is exact to **0.02 bp** off capture | ⛔ a claimed +10.7 % bias was a truth-parser bug (`TRAPS: a-truth-table-of-aggregates`) |
| ⚠ **the gDNA pmf under capture** | ⛔ the last unfixed length defect: **+13.6 bp** at a 330 bp gDNA mean and **+3.5 bp** at 120, drained per-bin `fit/true` **1.22 … 4.18** in the tail | ⭐ EXACT off capture and untouched by the drain, so it is a PLACEMENT problem — `EQUATIONS.md` §4.4 — and must not be attacked by editing `gdna_opportunity` |
| **message propagation** | ⛔ **OFF**, and it stays off until the tool is optimised end to end across all scenarios (owner, 2026-08-10) | net better on 3 of 4 strata (−58 / −44 / −32 %); **+155 % on unstranded × capture-ON**, which is the blind stratum. ⛔ It is the only remaining lever there, and its price must be RE-priced on the rebuilt ladder, never inherited |
| **the price of that** | ⚠ zero-gDNA golden scenarios go **0.029 → 89.93** and **0.005 → 9.58** | both AMBIG loci — the stratum above |
| ✅ **the LIBRARY figure is a FRAGMENT COUNT** | ⭐ it was an object-INCIDENCE sum in **all three** consumers that computed it, and the units did not match the truth it was scored against. Now `CalibrationResult.library_{gdna,rna}_fragments`, derived on read, each axis converted by its OWN `mass/count` | ⭐ On the pilot panel: better on **8 of 8** conditions, and the three contaminated ones the tool can see improve **73 / 81 / 98 %**. It RAISES the estimate on the contaminated arm and LOWERS it on the zero arm — both toward truth, which a merely-shifted number cannot do |
| ✅ **end to end — the LIBRARY figure** | ⭐ **RE-DERIVED 2026-08-11 on all 36 ladder conditions**, messages OFF, `calibration_truth_ab.py --cache-subdir _main`. Mean `\|err\|` per stratum: **0.0025 / 0.0053 / 0.0056** and ⛔ **0.4040** (max **0.9151**) on unstranded × capture-ON; `g00` control **0.0497** | ⛔ **NEVER quote the pooled 0.1043** — it is 97 % the one blind stratum (`TRAPS: never-pool-the-strata`). ⚠ The retired `0.1060` turns out to have been close, but it was unattributable and the stratified form is what ranks work. ⭐ The three good strata are ~2× better than the "0.005–0.012" row above, measured the same day |
| ⛔ **end to end — TRANSCRIPT assignment** | ⛔ **31.3 %** of RNA fragments misassigned (**55,916,700** of 178,399,996, 32 contaminated conditions), and a perfect prior removes **32.2 %** of that, to 21.2 % | ⭐ **RE-DERIVED 2026-08-11 on the rebuilt 36-condition ladder**, messages OFF, EM seed pinned — and it lands within 0.2 pp of the pre-rebuild 31.1 %, so the rebuild did not move it. ⛔ The ceiling is ONE STRATUM: unstranded × capture-ON goes **0.405**, the other three are **1.031 / 0.981 / 1.013** — a perfect prior is NEUTRAL-TO-WORSE on them. At GENE level the split is sharper still (**0.215** vs **1.704** on stranded × capture-OFF) |
| ⭐ **the pool level — gDNA vs NASCENT vs ANNOTATED** | ⭐ **NEW 2026-08-11**, `quant_accuracy.py --report` table ⑥. gDNA total is **−0.6 % / −6.4 % / −0.5 %** on three strata and ⛔ **−87.1 %** on unstranded × capture-ON; a perfect prior takes that last one to **−5.5 %** | ⛔ **Read `gDNA (total)`, never `gdna_est`** — the latter is `gdna_em_count` and EXCLUDES intergenic fragments, which are more than half of off-capture gDNA. Scoring it against the origin-count truth reads **−50.7 %** panel-wide that is entirely the missing pool; `score_library`'s docstring records the same mistake once fabricating an off-capture under-call (0.3151 vs a truth of 0.5000) |
| ✅ **the NASCENT leak — FIXED 2026-08-12** | ⭐ A synthetic nascent entity now receives **ZERO RNA prior**: it is one the index MANUFACTURED, nothing asserts it exists, so the null is that it is absent until the data proves otherwise. In scope: nascent **6,238,406 → 691,796** (0.111), gene Σ\|err\| **0.608**, `g00` gene **0.476** | ⛔ It was **NOT** the prior distributing uniformly — the shipped prior is a COMMON multiplicative factor and is provably neutral on the within-RNA split. The leak was the coverage-weighted warm start; `alpha = 0` is a new sparsity prior that counteracts it. ⚠ `is_synthetic`, NEVER `is_nrna` — a single-exon annotated transcript is flagged `is_nrna` and keeps its prior. ⚠ ~27 % of the freed mass goes to gDNA, not annotated |
| ⚠ **the price of that, and it is strand-gated** | ⚠ gDNA gains **+4.21 M** panel-wide (73 % of the freed mass returns to annotated). At the weakest in-scope strandedness the split is **86 % annotated / 14 % gDNA** | ⛔ On one nascent-ENRICHED toy (73 % nascent, a fully-nested antisense gene) it inverts to 83 % gDNA and costs 1721 → 1029 correct assignments at SS=0.65 — but **ZERO** at SS=1.0, because gDNA is strand-symmetric and RNA is not. Two strict xfails record it. ⭐ Next: an ACCURATE per-entity nascent prior instead of withholding from all (owner, 2026-08-12) |
| ⚠ **the `noop` BYTE-IDENTITY GATE CANNOT PASS — and the arm is SOUND anyway** | `--arm noop` differs from `base` on **376 of 2,016 fields**, max **399** on `count_abs_err`. ⭐ **The noise floor `--arm base_reseed` is 1,317 fields and 6,668** — so the arm's delta is **17× BELOW the floor** on the same field and touches 3.5× fewer fields. The injection wrapper is doing nothing, as designed | ⛔ **So the gate's SPECIFICATION is wrong, not the arm.** Byte-identity is unreachable while the scan is multi-threaded (`TRAPS: the-deliverable-is-not-reproducible-by-default`'s second source — the seed is pinned, the float banks still re-associate per run), so `arm_identity.py` is the wrong instrument for a `quant_accuracy` arm. ⭐ Either re-specify it as "below the `base_reseed` floor" or make the scan merge order-independent; ⚠ until then **no `quant_accuracy` arm delta below ~6,700 fragments is attributable** |
| ⛔ **the NASCENT channel** | ⛔ **20.2 M fragments** parked on entities whose truth is exactly **0** | 9.2 % of all true RNA, invisible in every transcript table |
| ⛔ **reproducibility — TWO independent sources** | ⛔ `EMConfig.seed` defaults to `None`, **and** the scan is multi-threaded so the float banks are re-associated by the per-worker merge | `TRAPS: the-deliverable-is-not-reproducible-by-default`. Pinning the seed alone is no longer sufficient. ⭐ Owner signed this off (2026-08-11): the spread is **1.503e-15** on `posterior_mean` and **exactly 0** on every integer column, and tests validate the float banks within a DERIVED, bracketed tolerance |

⭐⭐ **THE ONE SENTENCE, AND IT NOW HAS TWO HALVES BECAUSE THE TOOL HAS TWO DELIVERABLES.**

**On the LIBRARY gDNA fraction, the tool is accurate everywhere except unstranded × capture-ON, where it
cannot see gDNA at all** — `κ = ½` makes the strand λ-term identically 0 and no other channel reaches an
AMBIG slot. That is the whole of the calibration problem, and the oracle-prior injection proves it is the
*whole* problem for this deliverable: a perfect prior fixes it completely (`quant_accuracy.py --arm
oracle`). Everything in §2 is about giving that slot its own evidence.

⛔ **On TRANSCRIPT-LEVEL assignment, that sentence is FALSE and was never tested until 2026-08-07.** The
tool misassigns 15.6–20.6 % of fragments even on the three strata calibration handles perfectly, and
42.9 % at `g00` where there is no gDNA at all to get wrong. A perfect prior removes ~32 % of the panel
total and makes two strata *worse*. **These are two different problems in two different files, and
conflating them is how "is calibration the bottleneck?" got asked for months without an answer.**

⭐ **And the blindness now has a second, sharper statement, in the EM's own units.** The gDNA prior at
unstranded × capture-ON is not merely wrong, it is **PINNED**: `P/O = 0.040` at every one of
`g50 / g75 / g90 / g98`, with `P_g` sitting at 221 k → 428 k while the truth `O_g` runs 5.5 M → 10.8 M.
A flat line in the deliverable was already known; this is the same flat line one stage earlier.

⛔ **Read `mwae_all` / `Σ|err|`, never `solv%` / `mwae` / `conf-wrong` / `calib`** — those four share a
denominator the solver moves by declining to answer. And quote the SHIPPED column, not pass-0: a −37.2 %
pass-0 win was −3.9 % shipped (TRAPS: the-intermediate-is-not-the-deliverable).

## §1 ⭐⭐⭐ THE VALIDATION CAMPAIGN — items 1–3 DONE, item 4 open

⛔⛔ **OWNER RULING 2026-08-07, STILL IN FORCE: no new features until the tool is characterised end to
end**, with message propagation OFF and (until 2026-08-10) `length_likelihood` OFF, so every number is
attributable to the tool as it stands rather than to a feature landing mid-campaign.

| | the item | the question | state |
|---|---|---|---|
| **1** | calibration-prior-vs-oracle | calibration ships a PRIOR — how wrong is it? | ✅ `prior_vs_oracle.py` (+23 gates) |
| **2** | tool-absolute-accuracy | transcript-level accuracy, end to end? | ✅ `quant_accuracy.py --arm base` |
| **3** | error-downstream-of-calibration | inject the oracle prior — whose error is it? | ✅ `--arm oracle` |
| ⏸ **4** | performance | it has been a while, and new code means slowdown | **NOT STARTED** — §1.0 |

⭐⭐ **THE ANSWER TO 1–3, AND THE QUESTION HAD A HIDDEN ASSUMPTION.** "Calibration or the EM?" presumed
ONE deliverable. There are two and they answer differently:

* **the library gDNA/RNA separation** — calibration is the whole bottleneck. A perfect prior takes mean
  `|f_gdna − truth|` from **0.1060 to 0.0097**.
* **transcript assignment** — it is NOT calibration. **31.1 %** of fragments are misassigned and a
  perfect prior removes only **32 %** of that; on the two capture-OFF strata a perfect prior is **3–4 %
  WORSE**. 95.5 % of the stranded × capture-OFF error is ordinary isoform ambiguity (6.94 M → 311 k when
  isoforms are collapsed), which is not a Rigel-specific defect.

⛔ **Two things a perfect prior cannot touch, and they are §0's open rows**: the nascent channel (20.2 M
fragments on zero-truth entities, and at `g00` a perfect prior makes it WORSE, 4.30 M → 4.64 M), and
`g00`'s 7.05 M gene-level misassignment with no gDNA present at all.

⛔ **THE `g00` ZERO-gDNA CONTROL FAILS AT THE PRIOR, AND ITS WORST STRATUM IS NOT THE BLIND ONE.** The
shipped prior claims **2,067,637 gDNA fragments in libraries containing none**, and **1,707,321** of
them are at unstranded × capture-**OFF** — a stratum that reads healthy (`rel` 0.046) on every
contaminated row, and worse by **5.5×** than the blind one. Only a zero control could have found it.

⭐ **Where the live numbers are.** They are not repeated here: `prior_vs_oracle.py` regenerates the
`P−O` / `O−Fo` / `S−Fo` / `Fo−F` tables and the composition and strength diagnostics on demand
(~48 s/condition), and §0 carries the handful that rank the work. A closed item earns a line, not a
report.

### §1.4 ✅ `length_likelihood` — MEASURED, REFUSED AND DELETED 2026-08-10

Built, priced on the drained arm, and removed at `f470a570` along with every instrument that existed to
price it. The mechanism and the numbers are `TRAPS: a-linear-likelihood-emits-a-sign` and
`TRAPS: amplitude-fades-influence-does-not`; the panel lesson is
`TRAPS: a-single-level-panel-cannot-see-a-constant`. ⛔ **Do not rebuild it.** The opportunity-tilted
length moments it was built on survive in `effective_length.py`, where they were always a geometry rather
than a composition claim.

⚠ **What this closes and what it does not.** It closes the "give an AMBIG slot its own θ-independent
evidence" programme — there is no such channel and the search for one is over. It does NOT close
unstranded × capture-ON, which remains the one blind stratum; the surviving lever there is message
propagation, which is already measured to help exactly that stratum and is one config flag.

### §1.0 ⚠ ITEM 4 — AND THE PERFORMANCE SUBSTRATE IS A TRAP

Measured 2026-08-07 on one 10 M-fragment ladder condition (a 35,135-node chr22 index): per-locus EM
**15.9 s (47 %)**, native scan 6.5 s, calibration 6.5 s, second-pass drain 3.5 s, total 33.5 s.
⚠ Message propagation costs nothing measurable (33.7 s ON vs 33.5 s OFF) — the "the relay is the
expensive part" hypothesis is REFUTED.

⛔⛔ **BUT THAT PROFILE IS UPSIDE DOWN FROM THE REAL ONE AND BOTH ARE CORRECT.** Calibration is
depth-INDEPENDENT — every node in the index is solved regardless of read depth — so it scales with the
INDEX while the EM scales with the DATA. On real cfRNA at genome scale (~1.5 M nodes) calibration has
measured ~66 s against the EM's ~24 s: the exact reverse. **Profile on real cfRNA**
(`~/Downloads/rigel_runs/cfrna/`, genome index at `~/Downloads/rigel_runs/refs/rigel_index`), never on
this panel — `TRAPS: toys-rank-hotspots-backwards`, which cost a whole analysis once.

## §2 ⭐⭐ WHAT TO DO NEXT, IN ORDER

⭐ **The panels are rebuilt and gated** (2026-08-11: 184 payloads across ladder / flgap_short / flgap_long /
pilot, every non-target bank byte-identical). Every instrument under `scripts/design/` runs again — with
one exception, which is item 1.

1. ✅ **`ladder_arm_ab.py` — FIXED 2026-08-11, and the dead module path was the small half.**
   ⛔⛔ **THE FINDING: under the shipped `message_propagation = False`, 22 of the 26 arms could not move
   a number, and only ONE of them was supposed to be inert.** Six raised — four on
   `rigel.calibration.enrichment_frame`, deleted at `0d9d422b`. The other sixteen ran to completion,
   fired on tens of thousands of slots, and scored **byte-identical to `base`**: eight muted a relay that
   was already silent, and eight moved a PRECISION that nothing consumes, because `_fuse(own, p, msg)`
   returns `own` for every `p` when no message arrives. Only `intron_phi` and `kappa_half` were live.
   ⭐ **A fire counter cannot see this** — it answers "was my code called", and the question is "can this
   arm move a number". The repair is four things: `--messages {off,on}` stamped into **every row** so two
   arms run differently cannot be diffed silently; an up-front refusal for an arm the chosen policy
   cannot express; a fire assertion that names the arm (it covered `intron_phi` and `kappa_half` not at
   all, and the whole `zc_*` family only by "did ANY sibling fire"); and ⭐⭐ **`--self-test`, which runs
   every arm under BOTH policies and gates that each is INERT under one and MOVED under the other** —
   **25/25 in ~6 min**, falsified by three perturbations that each fired their own branch.
   ⚠ Two consequences: every arm file in `$RIGEL_ARMS` predates the `messages` field, so
   `arm_identity.py` now refuses to compare one against a fresh run — which is
   `TRAPS: re-record-the-baseline` doing its job, not a defect. And `anchor_opportunity_census.py` said
   in TWO places that its population was "exactly `strand_evidence`'s `struct_lock`"; measured, the
   solver's mask is **15–23× larger** (30,423 vs 1,312 at `g00`) because it also contains every
   zero-count NODE. The instrument's own mask is `g1_locked ∧ NODE` — the mask the standing xfail wants
   `struct_lock` rescoped TO — so its **346×** verdict stands but its SCOPE was overstated wherever it
   was quoted. Both sizes now print on every row.
2. ✅ **Re-derive the two end-to-end numbers on the 36-condition ladder — DONE 2026-08-11, BOTH halves.**
   `quant_accuracy.py --arm base / noop / oracle / base_reseed` on all 36 (7 min each at `--jobs 8`) and
   `calibration_truth_ab.py --pilot <ladder>/oracle_cache --cache-subdir _main`. The five §0 rows above
   carry the numbers; `$RIGEL_ARMS/qa_ladder_*.jsonl` carry the rows.
   ⭐⭐ **The verdict so far, and it is the one that ranks the next build:** a perfect prior is worth
   **32.2 %** of transcript misassignment and **every bit of it is unstranded × capture-ON**; on the
   other three strata it is neutral-to-worse, and at GENE level it makes stranded × capture-OFF **1.7×
   worse**. ⛔ Two things a perfect prior does NOT fix: the `g00` control (0.953) and half the nascent
   channel (15.86 M → 5.22 M, and at `g00` it goes the wrong way, 4.30 M → 4.64 M).
2b. ✅ **THE SIMULATION + BENCHMARKING WORKFLOW — DONE 2026-08-11** (owner ruling: build the workflow
   before the nascent panel). ⛔ **The diagnosis was THREE overlapping benchmarking stacks and no entry
   point**, with the working one undocumented: `scripts/design/*` (campaign-built, gated — the real one),
   `rigel.sim.analysis` (1,589 lines, a second scorer), `scripts/benchmarking/` (4,036 lines, external-tool
   comparison, 3–5 months stale, configs naming conditions that no longer exist).
   ⭐ **`scripts/sim/panel.py` is now the loop** — `status / build / simulate / cache / score / report`,
   every path derived from one panel YAML, adding NO measurement code. `TESTING.md` §2 is the full recipe;
   it previously stopped at "cache the scans" and **never mentioned running the tool or scoring it**, and
   the ORACLE cache had no step at all (it was a side effect of `pass0_vs_oracle.py`), so the documented
   recipe produced a panel every scorer refused.
   ⚠ **Deleted, and recoverable from `0d9d422b`:** `scripts/benchmarking/` entire; `rigel/sim/analysis.py`
   and `scripts/sim/evaluate_suite.py` (618 of 1,589 lines kept as `rigel/sim/net_flow.py` with its tests —
   the flow decomposition has no duplicate; the other 971 were a second scorer). `scripts/README.md` was
   rewritten: it documented three directories of which **two had not existed for months** and omitted
   `design/`, which is 56 files.
0. ⭐⭐⭐ **THE PER-TRANSCRIPT RNA PRIOR — THE FOUNDATION IS BUILT AND PROVEN; THE WEIGHTING FUNCTION IS
   THE WORK** (owner, 2026-08-12/13). The end goal is a prior that says WHICH transcript a locus's RNA
   pseudocounts belong to.
   ⭐⭐ **PROVEN 2026-08-13, `g00 ss0.99 capture_off`, true relative abundances as weights with the warm
   start zeroed:** transcript Σ|err| **4,275,988 → 696,900 (0.16×)**, gene **445,944 → 57,154 (0.13×)**,
   gene false-positive mass **23,529 → 0**. A deliberately FLIPPED allocation is **1.26×** worse. So the
   machinery converts good weights into good answers and a weighting function is worth designing.
   ⛔⛔ **THREE LIMITS TRAVEL WITH THAT NUMBER.** It is a CAPABILITY proof, not a ceiling — it does not
   price what a reachable weighting function could earn. True weights also hand over the true SUPPORT (a
   zero weight is EXACTLY absorbing, measured 0.0000 in all four mode × warm-start combinations), and
   4,579 of 8,750 annotated transcripts are silent at `g00`, so most of the fp_mass collapse is free
   rather than earned. And it is ONE condition on the stratum the tool already handles — the blind
   stratum is untested.
   ⭐ **What is built:** `splice_graph.build_transcript_path` (the transcript → ordered
   REGION/BOUNDARY/SPLICE-JUNCTION walk, 0/15,669 disagreements against the verified incidence, 45,609
   splice steps, transcription order on 7,312/7,312 minus-strand); the per-transcript weight LANE into
   the EM; `EMConfig.warm_start` = `coverage`/`prior`/`uniform`; the three-way composition published per
   object; `scripts/design/transcript_truth.py` (true per-transcript counts split by splicedness, 41 s
   per 10 M-fragment condition, seven named gates all passing).
   ⛔⛔ **STAGE 6 — THE WEIGHTING FUNCTION IS NOW BUILT, MEASURED AND REFUSED (2026-08-13), AND THE
   REFUSAL IS THE OWNER'S OWN THEOREM RATHER THAN A DEVIATION FROM IT.** The soft-min-along-the-path
   family — 12 arms, `min`/`harmonic`/`geometric`/`arithmetic` × three multipliers — is worse than
   `base` at transcript level on BOTH strata and its `g00` gene-level win (0.395–0.527×) collapses to
   1.006–1.041× on the blind stratum. §4 carries the full row; `TRAPS: an-upper-bound-is-not-an-estimate`
   carries the mechanism. ⭐ **The theorem is not what failed and must not be discarded**: the dial is
   monotone in its favour on every axis and both strata, so the pooled no-deconvolution control is the
   worst rung of the four. What failed is that a bound shared with a loud neighbour is not evidence about
   the quiet one, and 75.3 % of silent transcripts have such a neighbour.
   ⭐⭐ **WHAT THE EXPERIMENT SETTLED, AND IT NARROWS THE NEXT ATTEMPT SHARPLY:**
   * **The target is TOTAL abundance, not the unspliced count** — `oracle_alloc` (total truth weights)
     scores **0.163×** against `oracle_alloc_unspliced`'s **0.202×** at transcript level, even though the
     budget being split is literally unspliced pseudocounts. ⚠ At GENE level they invert (0.128 vs
     0.120), so this is a transcript-level ruling only.
   * **The density's units are settled and gated**: `density(r) = mass_rna_node[r]/ceff(r) = Σ_{t∋r} A_t/E_t`
     — the object's own opportunity CANCELS, which is what makes densities at objects of different sizes
     commensurate and what the `min` is a min OF. A transcript alone on its path recovers its own
     abundance exactly, on all four modes (`test_transcript_weights.py`).
   * **The support problem is the whole problem.** 0.0 % of expressed transcripts are ever zeroed; the
     error is entirely that 3,644 silent ones are not. Any next candidate is a SPARSITY mechanism, not a
     better mean.
   ⚠ The lane, the arms and the falsifications are all in place, so the next candidate is one function:
   `transcript_weights.build_weights` + one `alloc_*` arm.
   ⛔ An EARLIER first attempt was also built and REFUSED (§4) — soft-min over EXCLUSIVE objects with a
   per-object Jeffreys half. That one failed for three unrelated reasons and is a different row.
   ⭐ **Vocabulary ruling (owner):** REGION / BOUNDARY / SPLICE JUNCTION replace node / edge·line·cut /
   junction; `donor`/`acceptor` are banned as strand-dependent (the code is genomically ordered, so the
   name is backwards for 48.4 % of junctions). A CONVERGENCE, not an invention — all three axes already
   carry two live names each. `DESIGN.md` §0 is not updated until the bulk rename lands.
   ⚠ **Identifiability, measured, because it was the objection:** 70.43 % of fragments are compatible
   with 2+ annotated transcripts, but only **0.29 %** of true mass lies in likelihood-flat directions
   (0.96 % by an independent null-space measure) and a census over 952 loci found **zero**
   exactly-exchangeable groups. An oracle per-transcript prior is not meaningfully fake. ⛔ The real
   difficulty is the middle tier: **23.99 %** of true mass sits on 4,277 transcripts with no fragment of
   their own.
   ⚠ **The shipped warm start is NOT the cause of the tool's error** — a UNIFORM start scores 0.976–1.024×
   against base on four conditions, at the noise floor.
0b. ⛔⛔⛔ **ψ's COMPOSITION DOES NOT CLOSE — a defect, and its own work path** (found 2026-08-12 while
   publishing the three-way composition; owner ruled it a separate branch, taken up after the
   per-transcript prior lands). `NodeDeconv` asserts *"posterior means; f_pos+f_neg+gdna_frac = 1"*.
   Measured on `g00 ss0.99 capture_off`, every object addressed by a chain slot:

   | axis | sums to 1 | sums into (0,1) | sums > 1 |
   |---|---|---|---|
   | REGION (node) | 74.72 % | **25.25 %** — median 0.978, p5 **0.869** | 12 |
   | BOUNDARY (edge) | 77.24 % | **22.71 %** — median 0.979, p5 **0.850** | 16 |

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
0c0. ✅ **THE CONSERVED JUNCTION MASS IS PUBLISHED — DONE 2026-08-13, and it needed no re-scan.**
   `CalibrationResult.junction_conserved_mass` is `float64[n_junctions]`, the accumulator's `sj_mass`
   bank recovered exactly (12,758 of 13,482 elements bit-identical, worst 9.1e-13 on
   `g00 ss0.99 capture_off`). ⛔ It exists because the derivation was reachable only by knowing to
   perform it, while **`mass_rna_junction` sat next to it looking like the answer and reading 2.0719×
   high** — 5,668,526 incidences against 2,735,958.8 units of mass — so a weight that read it would be
   wrong in proportion to how spliced a transcript is. ⚠ **A PROPERTY, never a stored field, and that is
   forced rather than preferred:** `mass_rna_junction` is in `prior_vs_oracle.OVERRIDE_FIELDS`, so a
   cached array would survive an oracle arm's `dataclasses.replace` and describe the array it replaced
   (`TRAPS: a-hash-that-misses-its-artifact`, dataclass form — the same reason `library_rna_fragments`
   is derived). ⭐ `library_rna_fragments` now reads it, so the conversion has ONE home.
   ⭐ Six gates, four perturbations, and the sharpest is that a `where(count > 0, …)` fallback to the
   `1.0` identity would publish phantom mass on **4,636 of 13,482** zero-count junctions — a third of
   the axis, on the one axis that is certified RNA by construction. Exactly one gate catches that.
0c. ✅ **`sj_mass[2]` — THE PER-STRAND SPLICE-JUNCTION MASS, AND THE RULING IT REVERSES — LANDED
   2026-08-13** (owner, 2026-08-12), bundled with the two dead banks exactly as planned.
   ⭐ **The bundle is three parts pulling one way:** `JunctionEdge::mass` → `mass[2]` deposited at the
   same `col` the count uses (24 → 32 B); `Node::contained_length_sum` and
   `ContiguousEdge::unspliced_length_sum` deleted (24 → 16 B, 48 → 40 B). Net memory at genome scale is
   **down** ~19 MB. `sj_mass` moved `SINGLE_COLUMN_AXES` → `BANK_AXES`, and since the shape digest is
   DERIVED from those two tables it moved itself.
   ⭐⭐ **`substrate` folds the strand axis at the boundary**, so `PopulationView.mass` stays
   strand-agnostic and **no consumer below it changed at all**. The per-strand values are not re-exported
   there: their consumer is artifact filtering, which reads the payload.
   ⛔⛔ **THE TWO DELETED BANKS' STATED JUSTIFICATION WAS FALSE, and that is what made the deletion
   clean rather than lossy.** `scan_payload`'s docstring claimed `length_sum` removes the blind spot at
   `mu_g = mu_r`. It does not: there the third row is `(mu, mu)`, proportional to `(1, 1)` like the other
   two, so the 3x2 system is still rank one. It is an independent tilt only when the means already
   differ — precisely when the first two rows suffice. It survived because **nothing consumed it**, so no
   test could disagree. The retraction is written where the claim was, not merely deleted.
   ⚠ The lesson `mean_length` carried ("an object with no data must emit NOTHING, never a floored
   value") was MOVED onto `mass_per_crossing`'s 1.0-identity rather than deleted with its test.
0c1. ⛔⛔ **AND THE RE-SCAN GATE WAS BROKEN BEFORE ANY OF THIS — REPAIRED 2026-08-13.**
   `rescan_panels.py` compared every bank byte-exactly, on the stated grounds that *"these are integer
   tallies"*. The 2026-08-10 one-numeric-convention ruling made six of them float64, and float addition
   is not associative across worker threads: **two scans of the same BAM by the same binary differ on
   exactly those six banks**, by up to 3.5e-14 relative. The gate had therefore been UNSATISFIABLE since
   that ruling, and it failed this change on five banks the change never touched.
   ⭐ **The repair is a per-element budget `deposits · eps · value`**, with `deposits` read from the COUNT
   bank at the same object and the ×2 on the mass banks derived from the deposit rule — no tuned
   constant, and integers stay exact. `--self-test` 14/14, including that a ulp inside the budget is
   ACCEPTED, that a 1e-9 move is caught, and that a ZERO-count object gets a ZERO budget.
   ⛔ `TRAPS: a-stale-gate-accuses-the-newest-change`. What resolved it was a control on the INSTRUMENT
   — scan the same BAM twice — not on the change.
0c2. ⭐⭐⭐ **ORDERING RULING (2026-08-13): THE CACHES ARE REGENERATED BEFORE THE BULK RENAME, AND THE
   SECOND RE-SCAN THE RENAME COSTS IS ACCEPTED.** 13 of the payload's 62 schema names carry the old
   vocabulary (`node_*`, `edge_*`, `sj_*`, `cut_*`), and `payload_schema_digest` hashes NAMES — so the
   rename moves the digest again and refuses every cache a second time.
   ⛔ **Bundling them anyway would make the gate VACUOUS, and that is the whole argument.** Its power is
   *"31 shared banks, every one identical"*; rename the banks and there are almost no shared banks left,
   so the rebuild would be gated only by "the same symmetric difference appears everywhere" — which
   cannot see whether a renamed bank's VALUES survived. The gate has just demonstrated its worth twice in
   one session (it stopped an irreversible write, and it surfaced 0c1), and it only works with ONE change
   in flight (`one-thing-varied`).
   ⚠ A second consideration decided it as much as the first: **while every cache is refused, no
   instrument runs**, so a ~1,000-site mechanical rename would be made blind and then landed against a
   gate that cannot fail — the worst available combination.
   ⭐ **PREREQUISITE FOR THE RENAME, and it is easier to build while the caches are valid:** teach
   `rescan_panels.py` a `--renamed old=new` map, so a rename's rebuild asserts the new bank equals its
   predecessor within the same derived budget. Without it a pure rename cannot be gated at all. `accumulator.h` ruled that a mass bank is ONE value, not two, because *"nothing reads a
   mass per strand"*. ⭐ **That premise is now FALSE, which is what makes the reversal admissible:**
   artifactual splice junctions accumulate SYMMETRICALLY on both strands, like gDNA, so they are
   detectable by the strand model the tool already has — and a per-strand mass is the observable that
   model needs. ⚠ A second, structural reason: without it, artifact filtering needs TWO passes over the
   BAM (tally, filter, re-accumulate mass), which is the one thing the single-pass architecture exists
   to avoid.
   ⛔ **Record the reversal WITH the premise that changed**, or it gets re-litigated in both directions.
   ⛔ **Scope it deliberately.** The reversal is for the SJ axis, where the consumer is. Whether the
   BOUNDARY axis's `unspliced_mass`/`spliced_mass` also go per-strand is a SEPARATE question with no
   consumer named — `one-thing-varied` says do not bundle it.
   ⭐ **Price: the payload schema digest moves, so every scan cache is refused and rebuilds — 8.3 s per
   condition, ~6 min for the 36-condition ladder.** ⭐⭐ Bundle it with the two dead banks §2.7 already
   wants removed (`node_contained_length_sum`, `edge_unspliced_length_sum`), so the re-scan is paid once.
   ⛔ Drive it with `rescan_panels.py`, which gates every non-target bank on byte-identity and requires
   the SAME delta on every condition — that is only interpretable if ONE schema change lands at a time.
0d. ⭐ **THE INDEX SHOULD EMIT ITS DUPLICATE MAP** (owner, 2026-08-13). The index drops EXACT duplicates
   — same ref, strand and sorted exon tuple — and is right to: they are *"mathematically unidentifiable
   from any read data"* (`index.py`'s own words). ⛔ But it records NOTHING about what it dropped:
   `_collapse_duplicate_transcripts` logs a truncated summary and returns a filtered list. So a dropped
   id is unfindable through the index, and every consumer re-derives the map by re-parsing the GTF —
   which defeats the index being freestanding. Measured on the suite reference: 8,755 GTF transcripts, 5
   duplicate groups over 10 transcripts, 5 dropped → the 8,750 indexed; 2 of the 5 carry simulated
   fragments (329 at `g00`).
   ⭐ **The fix is an ALIAS MAP, not re-admission** — `dropped_t_id → kept_t_id` in the manifest, exposed
   as `TranscriptIndex.duplicate_map` defaulting to `{}` so an older index still loads. Re-admitting the
   transcripts would undo a collapse that is correct.
   ⭐⭐ **It needs an index rebuild but NOT a panel re-scan**, and that is worth stating because the two
   were conflated once: the scan cache is keyed on `graph_hash` / `reach_digest` /
   `payload_schema_digest` / `scan_config_digest`, and a metadata field on the same partition changes
   none of them. ⚠ A rebuild can still move `reach_digest` on its own (`reach` is covered by no other
   hash, and a rebuild once moved 38 % of human reaches), so verify with `rescan_panels.py` rather than
   assuming the caches survive.
   ⚠ Then DELETE `transcript_truth.py`'s local `duplicate_map` and read the index's — moving it, not
   copying it, so there is one home.
3. ⭐⭐ **Price the CANCELLING PAIR together** — `struct_lock` rescoped to `g1_locked & NODE` AND the
   `intergenic|exon` seam claiming its RNA-contaminated crossing mass as gDNA. Five of the seven
   remaining xfails go green if and only if this lands. ⛔ Neither half may be priced alone: the
   `struct_lock` half is **+3,207 %** on the zero-gDNA control by itself.
   ⛔⛔ **AND IT CAN ONLY BE PRICED WITH `--messages on`.** `zc_struct_lock_g1` is MEASURED inert with
   them off (2026-08-11): it rewrites the mask on 30,423 slots and changes no scored field, because
   `struct_lock` reaches the answer only through `own_composition_logvar` → a precision → a fusion that
   `SilentPolicy` never performs. Priced against a messages-off base it would read as "worth nothing",
   which is the exact reading item 1 was fixed to make impossible.
4. **Re-price message propagation on unstranded × capture-ON.** The only remaining lever on the one blind
   stratum, one config flag, and it closes two more xfails. Its recorded +155 % predates the
   conserved-count rewrite AND the numeric convention. ⛔ RE-price, never inherit.
   ⭐ It is now one pair of commands — `ladder_arm_ab.py --arm base --messages off` against
   `--arm base --messages on`, `--jobs 8` — with the setting stamped into every row of both.
5. **`mass_*_edge` → `count_*_edge`.** THREE of the five `mass_*` fields on `CalibrationResult` are
   incidences (`mass_gdna_edge`, `mass_rna_edge`, `mass_rna_junction`) and two are fragment counts. Its
   own commit, `one-thing-varied`.
6. **Restore the moment tests.** `contained_moments` / `crossing_moments` / `build_slot_moments` moved to
   `effective_length.py` without their tests, which went with the deleted length channel — untested
   geometry in a live layer-2 module.

7. **The cheap ledger** — moved here 2026-08-11 from a dev doc that was deleted, so it has one home:
   * dead substrate surface: `PopulationView.length_sum` / `.mean_length` / `.total_inv_length_sum` have
     no consumer in `src/`
   * the two dead banks `node_contained_length_sum` and `edge_unspliced_length_sum` ⛔ removing them
     invalidates every cache, so bundle with the next re-scan rather than paying it twice
   * `second_pass.py:443-445`'s comment — measured false in premise, harmless in effect
   * `module_census.py` re-run; the purge left stale sibling references
   * `UNDOCUMENTED_DEBT` — 8 scripts on disk and unlisted, in `tests/test_scripts_index.py`

⚠ **Tracked, not blocking:** two exact per-fragment mirrors of a PURE-RNA library give different
deconvolutions — `mass_gdna_node` +5.4 % at `ss=0.75`, +2.0 % at `ss=1.0`, 0 at `ss=0.90`. Not
boundary-only and not monotone. It is 5.4 % of the zero-gDNA false-positive channel, and an R1-sense
library is now simulable, so it is measurable.

⛔⛔ **The two briefs below are ON HOLD by owner ruling (2026-08-07)** and are kept only because the
analysis is complete and re-deriving it would be waste — not because they are next.

### §2.1 THE BRIEF FOR **the-cancelling-pair** — everything that experiment needs, in one place

⭐ Promoted from a working doc when it was deleted (2026-08-07). This is the next build.

**THE TWO HALVES, AND WHY NEITHER HAS AN HONEST PRICE ALONE.** The certified-RNA channel is a LOWER BOUND
(`ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)` — the exon may hold molecules that never touch that seam) and ψ delivers it
as a two-sided Gaussian, penalising a destination holding MORE RNA exactly as hard as one holding less.
Making it one-sided — `−½·p·max(0, mo − log f)²`, no new constant — is **the only mechanism the zero-gDNA
control has ever ENDORSED: −81.9 %, 8/8.** The panel is nonetheless **+5.2 %**, and the reason is
`TRAPS: a-cancelling-defect-pair`: the two-sided term was doing TWO jobs, the bound *and* — by its upward
side — a de-facto gDNA LEVEL channel. Removing the accident exposes the real gap, +7.7 % on exactly the
stratum with no working level channel and −9.2 % where the prior already suffices.

⛔⛔ **AND THE LEVEL CHANNEL IS STRUCTURALLY DISCONNECTED, WHICH IS A THEOREM AND NOT A MEASUREMENT.**
Receivers of a one-slot-step channel at `g00`, by slot type:

| channel | EDGE | intergenic NODE | intron NODE | exon NODE |
|---|---|---|---|---|
| gDNA level | 2,592 | 0 | 0 | **0** |
| certified RNA | **0** | 0 | 0 | 10,493 |

The chain is strictly `N E N E … N`, so a one-slot-step channel is **BIPARTITE**: a NODE emitter can only
reach an EDGE. The only licensed originators of a gDNA level are structurally pure-gDNA objects, which are
NODEs — so **no NODE can ever receive a gDNA level.** It is not weak, it is disconnected from every object
carrying mass. ⛔ Two repairs are already refused: letting an EDGE originate a level (a patch on a
symptom), and letting the level cross two steps (the chain-fused level is dominated by the 1,312 intergenic
anchors and hands every exon the off-probe floor, against a true neighbouring density **346×** higher).

⚠ **The measurement to run it against**: at `g01 ss0.50 capture_on`, ψ-boundary ablation with the identity
exact, every single channel ablation is small (+2.3 / +6.8 / +4.9 / −0.1 %) and the joint one is +60.7 % —
`TRAPS: all-small-singly-large-jointly`, so the four share an upstream quantity. And at HEAD's top-12 error
slots the self-solve WITH the fitted prior is nearly correct while the messages destroy it; muting the
certified-RNA channel alone recovers the truth at 8 of the 12.

### §2.2 ⛔⛔ THE LENGTH CHANNEL CANNOT BE PRICED ON THE LADDER — the substrate is BUILT AND CACHED

`length_likelihood` defaults **False** and is the **only** channel that can give an unstranded slot its OWN
composition evidence: it is θ-independent, so the Schur complement that zeroes the strand term at an AMBIG
node does not apply, and it enters `tau_lam` ungated. ⛔ But `TRAPS: equal-lengths-carry-no-composition` —
at equal component mean lengths the 2×2 is identified only through `μ_g − μ_r`, and the ladder was built
with equal configured lengths: the realised post-capture gap is **+1.5 %**. Enabling it there would measure
approximately nothing and that would be read as "the feature does not work".

⭐⭐ **THE PANELS ARE BUILT AND THE ORACLE CACHES EXIST (2026-08-07)**: `suite/flgap_short` (realised gap
**−41.0 %**) and `suite/flgap_long` (**+40.4 %**), 4 conditions each, 32 MB of cache each, verified
byte-identical on a cached re-run (2m39s cold → **21 s** warm). Nothing blocks the experiment but the
campaign in §1.

⭐ **And both panels reproduce the ladder's failure exactly, which is what makes them a clean A/B** — with
messages off and the length channel off, only unstranded × capture-ON is broken:

| | truth | reported | | | truth | reported |
|---|---|---|---|---|---|---|
| `flgap_short` ss0.50 OFF | 0.4854 | 0.4771 ✅ | | `flgap_long` ss0.50 OFF | 0.4933 | 0.4843 ✅ |
| ⛔ ss0.50 **ON** | 0.5288 | **0.0324** | | ⛔ ss0.50 **ON** | 0.5895 | **0.1048** |
| ss0.99 OFF | 0.4854 | 0.4785 ✅ | | ss0.99 OFF | 0.4932 | 0.4859 ✅ |
| ss0.99 ON | 0.5290 | 0.5147 ✅ | | ss0.99 ON | 0.5894 | 0.5764 ✅ |

A fix should move exactly those two rows and leave the other six alone.

## §3 ⛔ THE VERTEX PROBLEM — CLOSED, kept as one paragraph so it is not rebuilt

ψ lands 5–8 % short of `f_g = 1` on unexpressed genes, and that was ranked #1 for a campaign. It is
**not a build**: `vertex_ceiling.py` prices it and evidence-free objects sit at median `|Δ|/sd(f_g)` =
**z = 0.5–0.6**, inside their own 1σ — the 24.4 % it measures is *the value of missing information, not
headroom*. No prior-free solve can beat it: every proper prior on [0,1] has a median strictly inside
(0,1). ⭐ `test_certified_rna_licence.py`'s **zero-spliced-count** gate independently closed the one channel that might have
supplied vertex evidence — a zero certified count is consistent with `f_g = 1` too. ⚠ Its companion
hypothesis, the **phantom-gDNA floor**, turned out to be a DIFFERENT and much larger bug and is now
fixed: TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero, §0's 39 %.

## §4 ⛔ WHAT IS DELIBERATELY NOT NEXT — one line each, so it is not rebuilt

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_pin.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` EDGE** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both EDGEs the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | `length_ceiling.py`, one pmf at a time | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores junction *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. §3 below |
| **§2's Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`node_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to **§4**. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| ⛔⛔ **A SPECIFIC per-transcript allocation RULE — soft-min over exclusive objects with a per-object Jeffreys half** | built end to end and A/B'd on all 36 conditions, seed pinned | ⛔ **REFUSED — worse on EVERY stratum and on the zero control**, transcript Σ\|err\| **57.5 M → 81.6 M (1.42×)**, and a length-proportional variant **2.10×**. ⭐ **The MECHANISM is not what failed** — the gDNA:RNA split moved **+0.2 %**, exactly as the conservation identity requires, so the A/B priced the ALLOCATION alone. Three defects, each measured: exclusivity hard-zeroed **38.7 %** of transcripts; estimating a density on a tiny exclusive region and extrapolating over the whole transcript amplified variance up to **6,534×** (44.6 % of weighted transcripts had their density from <200 bp); and a per-object `+½` revived the silent half of the annotation (`frac_expressed: 0.5`), taking false-positive mass **18.6 M → 41.6 M**. ⛔ The rule and its config flag were DELETED; the LANE it rode on was kept and is now proven (item 0). ⚠ The Jeffreys-mean half of it is §4.1 graveyard row one — see `TRAPS: a-trap-names-the-defect-not-the-repair` |
| ⛔⛔ **THE SOFT-MIN-ALONG-THE-PATH WEIGHTING FUNCTION — the owner's own theorem, built faithfully and REFUSED** | 12 arms (4 modes × 3 multipliers) on `g00 ss0.99 capture_off`, 3 of them re-run on the blind stratum `g50 ss0.50 capture_on`, base re-recorded in the same session | ⛔ **REFUSED. Worse than `base` at TRANSCRIPT level on every rung of every arm and on both strata** — 1.317–1.604× at `g00`, 1.262–1.331× on the blind stratum — with transcript false-positive mass 1.76–2.20× worse. ⭐ **The one encouraging number does not survive:** at `g00` the same weights took GENE error to **0.395–0.527×** (against `oracle_alloc`'s 0.128×), and on the blind stratum that collapses to **1.006–1.041×**, i.e. nothing. ⭐⭐ **The MECHANISM is `TRAPS: an-upper-bound-is-not-an-estimate`, and it is structural rather than a tuning failure:** the theorem bounds a transcript by the thinnest object on its path, but **3,644 of 4,839 silent transcripts (75.3 %) share an object with an expressed one** and inherit its bound. The zero-weight SET was byte-identical across all twelve arms, because a bound is zero only when every object is — a property of the data, not of the estimator. ⭐ Retreating to GENE granularity (split within the gene by effective length alone) did NOT rescue it — 1.340× against 1.317× — which proves the damage was never the within-gene split. ⚠ **Two things the experiment ESTABLISHED and which the next candidate should keep:** the dial is monotone in the theorem's favour on both axes and both strata (`min` < `harmonic` < `geometric` < `arithmetic`), so the pooled `Σmass/Σopportunity` control is the WORST rung and the soft min is doing real work; and **0.0 % of expressed transcripts were ever zeroed**, so the unrecoverable failure direction never fired. `scripts/design/transcript_weights.py`, `tests/calibration/test_transcript_weights.py` (31 gates) |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |

### §4.1 ⛔⛔⛔ THE GRAVEYARD — ELEVEN MECHANISMS PRICED, ELEVEN REFUSED. DO NOT REBUILD THESE.

⭐ Promoted from a working doc when it was deleted (2026-08-07). Every row is a
real build that was measured and refused. ⛔ **`g00` is the owner-required ZERO-gDNA control: its truth is
exactly 0, so every fragment there is a false positive with nothing to cancel it** — which is why a
mechanism can look good on its target and still be inadmissible.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ NODE` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| §2's mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |

⭐⭐ **THE PATTERN, AND IT IS THE REAL RESULT: every one of the eleven was a rule for how to resolve DOUBT,
and at `g00` the doubt must resolve to NO gDNA.** A rule that lifts an evidence-free slot off zero is
inadmissible there however well it scores elsewhere. ⛔ The only candidate the control has ever ENDORSED is
the one-sided certified-RNA bound (−81.9 %, 8/8), and that one is panel-negative alone — it is half of
`TRAPS: a-cancelling-defect-pair` and is §1's **the-cancelling-pair**.

⚠ Four more `zc_*` arms exist as decomposition REVERTS used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.

---

---

## §5 THE OPEN QUESTIONS, ranked by how much they would change a decision

1. ⭐⭐⭐ ~~**WHY IS PRIOR FIDELITY ANTI-CORRELATED WITH DELIVERABLE QUALITY?**~~ ⭐⭐ **LIKELY ANSWERED,
   AND THE ANSWER IS THAT IT IS NOT THE PRIOR — IT IS THE MESSAGES** (2026-08-06, ψ-boundary ablation on
   HEAD at `g01 ss0.50 capture_on`, TRAPS: byte-identity-gate exact). At every one of HEAD's 12 worst slots the **self-solve WITH
   the fitted prior is nearly correct** (0.0043–0.0344 against truths of 0.0007–0.0136) and **the message
   layer then destroys it** (0.05–0.83). Muting the certified-RNA channel alone recovers the truth at 8 of
   the 12. So a better prior cannot show up in the deliverable: it is already right at the objects that
   carry the error, and something downstream overwrites it. ⛔ That also re-reads the earlier evidence
   rather than contradicting it — "the prior is degraded on unstranded where the solve improves" is what a
   prior looks like when it is not the binding constraint. ⚠ Confirm on a second stratum before closing.
   §1.1 above, TRAPS: the-intermediate-is-not-the-deliverable.
1b. ⭐⭐ ~~**Is the phantom-gDNA floor the same bug as the vertex problem?**~~ ⛔ **ANSWERED, and NO** — it was
   a different and much larger bug (§0's 39 %), and its residual is the CAPTURE LEVEL defect, not the vertex.
   TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero/TRAPS: the-divergence-was-a-barrier.
2. ⭐⭐ **What is the joint price of the splice-flux reframe and the level defect?** **the-cancelling-pair**. Neither has an
   honest price alone, and one of them is already in `src/`.
3. ⭐ **Do we solve ALTERNATIVE SPLICING correctly?** The `alt_splice` rung exists and is **unverified** — see
   the handoff. Cheap, and it is the only structure where several junctions share an EDGE.
4. ⚠ **The two capture-ON pilot rows disagree about the SIGN of every length correction**, across two
   independently built panels. Unexplained. ⛔ Do not average them; find which one is lying.
5. ⚠ **Does the reframe's own `σ²_transfer` correctly price a ratio built on 19 counts?** `EQUATIONS.md`
   §3.5d says it is the right medicine for the noise and the wrong medicine for the share-weighting.
   Unmeasured.

---

## §6 THE RULES THAT MADE THESE NUMBERS TRUSTWORTHY

⭐ These are `TRAPS.md`'s operational core, repeated here because they are what a new session must do on day
one, not read about later.

- **Measure the CEILING before building the correction** (TRAPS: measure-the-ceiling-first). It has re-ranked this project
  **five** times, and it is also how you learn a phase is finished — §0 is that story for Stage A.
- **Run the PANEL arm before writing a mechanism into `src/`** (TRAPS: panel-before-src). Three toy-positive changes have been
  panel-negative.
- **Zero-gDNA and zero-RNA controls on every experiment** (owner, 2026-08-05). The truth is a constant, so
  every deviation is a false positive with nothing to cancel it — and **re-anchor-on-the-deliverable** was found this way.
  ⚠ But check the arm could have fired at all (TRAPS: could-the-arm-have-fired).
- **Quote `mwae` over ALL objects and the raw Σ|err|**, never the honesty metrics alone (TRAPS: honesty-metrics-reward-ignorance).
- **One thing varied**, with the baseline re-recorded from the current tree in the same session (TRAPS: re-record-the-baseline).
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire**
  (TRAPS: perturb-every-gate). And name the observable for *each place* the change was made (TRAPS: name-the-observable-per-site).
