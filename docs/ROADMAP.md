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
| ⛔ **end to end — the LIBRARY figure** | ⛔ **NOT CURRENTLY KNOWN, and the old `0.1060 / 0.0097` pair must not be quoted.** It was unattributable: its `§1.3 ①` pointer named a section this file does not contain, and its two halves came from different instruments | ⚠ `calibration_truth_ab.py` reads mean `\|err\|` **0.1185** on the pilot panel — but that is 4 conditions on TWO gDNA levels, dominated by the blind stratum (`TRAPS: a-single-level-panel-cannot-see-a-constant`). ⭐ **The 36-condition ladder is rebuilt; re-derive there.** §2 item 1 |
| ⛔ **end to end — TRANSCRIPT assignment** | ⛔ **31.1 %** of fragments misassigned, and a perfect prior removes only **32 %** of that | 67.5 % survives, and on capture-OFF a perfect prior is 3–4 % WORSE. ⚠ Measured 2026-08-07 on the pre-rebuild panel; re-derive alongside the row above |
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

1. ⛔⛔ **`ladder_arm_ab.py` IS BROKEN and it blocks items 2–4.** Lines 264 / 385 / 406 do
   `sys.modules["rigel.calibration.enrichment_frame"]` on a module DELETED at `0d9d422b`, so the
   `zc_transfer` and `zc_reference_var` arms raise `KeyError` today. This is the panel harness — nothing
   below can be priced until it works. ⚠ Also check `anchor_opportunity_census.py:136`, whose comment
   claims its population is "exactly `strand_evidence`'s `struct_lock`" while it builds
   `g1_locked & INTERGENIC & NODE` — matching neither mask, in the instrument whose 346× verdict
   `CLAUDE.md` quotes.
2. ⭐⭐⭐ **Re-derive the two end-to-end numbers on the 36-condition ladder.** §0's LIBRARY row is
   unattributable and its TRANSCRIPT row predates the rebuild. `calibration_truth_ab.py` for the first,
   `quant_accuracy.py --arm base` / `--arm oracle` for the second. ⛔ The pilot panel has two gDNA levels
   and cannot carry either (`TRAPS: a-single-level-panel-cannot-see-a-constant`).
3. ⭐⭐ **Price the CANCELLING PAIR together** — `struct_lock` rescoped to `g1_locked & NODE` AND the
   `intergenic|exon` seam claiming its RNA-contaminated crossing mass as gDNA. Five of the seven
   remaining xfails go green if and only if this lands. ⛔ Neither half may be priced alone: the
   `struct_lock` half is **+3,207 %** on the zero-gDNA control by itself.
4. **Re-price message propagation on unstranded × capture-ON.** The only remaining lever on the one blind
   stratum, one config flag, and it closes two more xfails. Its recorded +155 % predates the
   conserved-count rewrite AND the numeric convention. ⛔ RE-price, never inherit.
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
