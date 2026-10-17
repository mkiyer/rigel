# Rigel 0.8.0 — independent verdict on the release reassessment (2026-10-12)

> **Superseded as a plan by `RNA_RELEASE_0_8_0_PLAN.md` (same day), after the owner's review of this
> verdict. Corrected here:**
> (1) Useful delivered information is not confined to pieces gDNA cannot see. At objects admitting both RNA
> strands, two strand counts cannot separate gDNA from balanced RNA at any strand specificity. Own evidence
> fails there: test g98 stranded ON goes 52.69 → 72.21 %, with three antisense genes overcalled by about
> 2,000 fragments each. §0 point 4, §6 and §9 understate this.
> (2) "Meets the tag at gene level on every in-scope stratum" holds for stratum aggregates only. Test g98
> stranded ON genes are 40.83 % against the tag's 26.07 %.
> (3) The per-length likelihood in §2 needs component-specific observed length laws, which the pipeline does
> not have. It is not interchangeable with the length-marginal form, even under flat capture.
> (4) The VCaP ±0.002 stop rule and the hybrid's opportunity cut are proposals, not gates or mechanisms.
> Offline merges of two arms' weights must first align their median units. The hybrid figures quoted below
> survive alignment; at ladder g05 stranded ON the units differ by 4.88 nat.
> The shared-field verdict and the release-bar table stand.

Sandbox, not authoritative. Answers `RNA_RELEASE_REASSESSMENT_PROMPT.md` against the packet
`RNA_RELEASE_REASSESSMENT_REVIEW.md`, the code at main `624dfade`, the tag `ff26d8e5` and the archived evidence.
No source change, no installation change, no commit, no pipeline run. The new measurements here are read-only
analyses of archived arrays (seconds each); scripts and outputs are in
`.cache/rigel_runs/2026-10-12_independent_review/`.

## 0. Verdict

1. **Abandon the shared capture field.** Its measured gain over shipped belongs to its new numerator/normalizer
   geometry, not to integration (g50: the fixed shipped field under the new geometry beats the integrated field on
   all four rows). At zero gDNA its model is unidentified along exact isoform/capture ridges, so integration turns
   the width of a broad prior into allocation error. Its field is also not the capture physics.
2. **The release bar is the tag, and main is one mechanism away from it.** On byte-identical panels, main beats
   tagged 0.7.1 on every in-scope stratum except two: the ladder's g00 stranded-OFF control (13.61 % vs 3.70 %)
   and the test chromosome's contaminated stranded-OFF rows (9.35 % vs 8.20 %). Both are capture-OFF rows where
   the reader invents capture.
3. **The single next mechanism is the own-evidence reader**: delete the delivered composition, DNA-level and
   RNA-level factors from the capture reader's likelihood. It is a deletion, and it is already measured on every
   panel. It clears both failures (2.79 %, 8.00 %). It meets the tag on every in-scope stratum at gene level. At
   transcript level it misses on two strata by 0.10 and 0.01 points, which I read as parity.
4. **Its cost is localised and known.** The delivered DNA level carried real information only on pieces that
   gDNA cannot see: exons shorter than gDNA fragments, whose bases are observed only by crossing fragments.
   Losing it costs RNA-short stranded ON 14.24 → 20.43 % (tag 20.33 %). Whether tag parity there is acceptable
   is the one decision the numbers leave to the owner.
5. **After that, in order and only if wanted:** admit neighbour gDNA levels where an object cannot see its own
   bases (derived, not cut); then the numerator/normalizer repair, isolated factor by factor under the shipped
   VB.

## 1. The release bar, re-measured from archives

Like-for-like checks: the five panel manifests hash identically to the October 5 tagged sweep's
(`fa4b7d31…`, `60b3c79f…`, `24fbaf69…`, `9ed0f915…`, `f3abaca1…`). Both scorers compute Σ|Δ| over Σ truth per
stratum, with g00 kept apart. Every run is fractional with seed 0. `src/` is identical between `0fdb9dce`, where
the 270-condition matrix was scored, and `624dfade`. The tag ran with CLI `--threads 1`; the matrix pinned only
the scan, which is the documented reproducibility condition.

Transcript % (gene %), ratio of sums. Tag = v0.7.1; Oct 5 = the detector-era tree; own = the own-evidence reader
(the matrix's `unary` arm).

| Panel · stratum | Tag | Oct 5 | Main | Own |
|---|---:|---:|---:|---:|
| Ladder str OFF, contaminated | 3.95 (0.44) | 2.02 (0.25) | 3.37 (0.33) | 2.47 (0.30) |
| Ladder str OFF, **g00** | 3.70 (0.47) | 1.93 (0.13) | **13.61 (0.83)** | 2.79 (0.16) |
| Ladder str ON, contaminated | 9.91 (2.52) | 6.14 (1.12) | 6.05 (1.37) | 6.37 (1.39) |
| Ladder str ON, g00 | 20.22 (3.47) | 6.47 (2.36) | 15.10 (2.61) | 6.61 (2.43) |
| Ladder unstr OFF, contaminated | 10.29 (1.49) | 2.29 (0.31) | 2.40 (0.32) | 2.34 (0.31) |
| Ladder unstr OFF, g00 | 20.09 (2.53) | 1.60 (0.14) | 1.59 (0.14) | 1.59 (0.14) |
| RNA-long str OFF / ON / unstr OFF | 5.23 / 7.38 / 5.29 | 2.26 / 6.17 / 2.13 | 3.52 / 6.09 / 2.55 | 3.47 / 6.54 / 2.21 |
| RNA-short str OFF / ON / unstr OFF | 6.31 / **20.33** / 6.00 | **42.30** / 19.82 / 3.78 | 4.99 / 14.24 / 3.95 | 4.56 / **20.43** / 3.94 |
| Test str OFF, contaminated | 8.20 (1.20) | 7.39 (1.08) | **9.35 (1.11)** | 8.00 (1.10) |
| Test str ON, contaminated | **8.58** (2.45) | 13.89 (4.01) | 8.17 (0.90) | **8.59** (1.28) |
| Test unstr OFF, contaminated | 9.25 (1.76) | 8.78 (1.33) | 9.18 (1.39) | 8.64 (1.33) |
| Test g00 str OFF / str ON / unstr OFF | 11.35 / 14.22 / 9.38 | 7.00 / 7.21 / 5.71 | 7.66 / 8.22 / 5.62 | 7.70 / 6.81 / 5.61 |

RNA-short stranded-ON genes are tag 3.52, main 1.50, own 1.74. Deferred (unstranded × ON, reported, out of
scope): ladder contaminated 33.21 / 12.23 / 62.96 / 52.76; RNA-long 28.99 / 8.24 / 34.11 / 14.75; RNA-short
29.76 / 18.84 / 68.06 / 63.71; test 18.94 / 18.15 / 13.66 / 21.33. Controls, main → own: single-probe stranded ON
11.46 → 8.67; junction-probed 28.21 → 27.96; test .70-specificity ON 7.88 → 10.34 (tag 10.28). VCaP gDNA fraction
is 0.2528 → 0.2532 against an approximate 0.2518.

Per rung, own trails the tag at test g05/g25 stranded OFF (7.20/7.73 vs 6.96/7.65) and test g05/g98 stranded ON
(7.24/72.21 vs 6.03/50.71). Main trails at six test rungs and at the ladder's g00 stranded OFF. The g98 rows
carry about 2 % RNA, and every arm's relative error there is large. All of this is in `bar_table.out`.

What this says: the October 5 tree shows what flat-at-g00 bought, and its 42.30 % shows what a library-level
reference costs. Main traded that catastrophe for false capture on flat stranded fields. No tree beats the tag
everywhere, but the own-evidence reader comes closest, and it is a deletion.

## 2. Estimand and likelihood (decision 1)

**Estimand.** The captured fragment count of each locus component (annotated transcripts, synthetic spans, the
gDNA component); TPM on plain lengths by ruling. Capture matters only through where each component's fragments
can land.

**Likelihood, conditional on mapping, annotation and observed length.** Component `k` with captured fraction
`a_k`, a fragment of length `w` on path `h`:

    p(f | k) = p_k^obs(w) · π_kh(f) / Z_k(w),        Z_k(w) = Σ_{s admissible for (k, w)} π_k(s, w)

Here `p_k^obs` is the component's observed length law and `π` the capture probability of the path's footprint. It
is exact: the observed law already holds capture's length selection, and the positional selection is normalised
per length, so no precursor law is needed. The production form uses one length-marginal normaliser
`L_k = Σ_w p(w) Z_k(w)` and drops `π_kh(f)` from the numerator. That equals the form above when
`π` does not depend on `w` and every candidate of `f` shares its footprint.

**What the conserved-share field is.** It is exact bookkeeping: each fragment's unit is partitioned over the
objects it deposits on, and lengths are conserved
([capture_eff_length.py:12–41](../../src/rigel/calibration/capture_eff_length.py)). As a selection model it is
first-order: `π(f) ≈ Σ_o d_o(f) c_o`, linear and independent of length. The simulator's selection is
`off_target + binding_per_base × best single-part overlap`
([sampler.py:12–19](../../src/rigel/sim/capture/sampler.py)). Three consequences, each a counterexample to
treating the field as physics:

- **Junctions are a max, not a sum.** The simulator's split-probe fixture has a true junction average of 44.56
  against an adjacent cap of 27.78 (follow-up review). The additive rule
  ([capture_eff_length.py:162–210](../../src/rigel/calibration/capture_eff_length.py)) can also return **zero**
  for a junction between captured exons when the intron pieces beside it carry imported levels. On RNA-short
  stranded ON, 2.6 % of junction cuts are exactly zero under the shipped reader (16 % of the junctions beside
  intron pieces with near-zero gDNA opportunity) and none under own evidence (`junction_dim.out`).
- **Capture depends on fragment length.** Inside a probed exon, binding is proportional to overlap, so a 250-bp
  gDNA fragment binds about 3.3× more than a 75-bp RNA fragment. The `c_o` read from gDNA is therefore not
  RNA's `c_o` under a length gap. The ruler shows a mirror image: isoform L/Y relative to gDNA is about 1.96 on
  RNA-short and about 0.45 on RNA-long, against about 1 on the equal-law ladder. This is a hypothesis, not proof
  of an allocation error: the unit audit shows raw L/Y can differ while K/L is right, so the test is K/L against
  the oracle's per-fragment probability.
- **Short pieces are invisible to gDNA.** A piece shorter than the gDNA fragments holds no contained gDNA; its
  bases are observed only by crossing fragments, which deposit on its boundaries. There are 5,829 such enriched
  regions on RNA-short stranded ON (§6).

**A frame inconsistency exists today.** The scorer reads the realized gDNA law
([pipeline.py:825](../../src/rigel/pipeline.py)). The gDNA normaliser is built from the uniform-frame law
([pipeline.py:944](../../src/rigel/pipeline.py) → `gdna_region_eff_len` →
[priors.py:334–341](../../src/rigel/calibration/priors.py)). Under capture, Σ K/L in a 100-bp window is 1.071
(ladder), 1.009 (RNA-long) and 1.024 (RNA-short) instead of 1. RNA uses one observed law in both places, which
is consistent.

**Verdict on decision 1.** The shared field integrates uncertainty over a misspecified linear field. That cannot
remove the field's bias, and the field's posterior is not identified where it matters (§3). Two parts of the
work are worth keeping, both deterministic: the frame consistency and candidate-path numerators where paths
differ.

## 3. Identification (decision 2)

**Exact ridge at zero gDNA** (`ridge_check.out`). Take T1 = {A} and T2 = {A, B}, with ℓ_A = 1000 and ℓ_B = 200,
flat truth and equal densities, so n_A = 2000 and n_B = 200. For every γ = c_B/c_A ≥ 0.5, the abundances
`n_T2 = n_B (1 + ℓ_A / (ℓ_B γ))` reproduce (n_A, n_B) exactly: log-likelihood −670.1994141683 at γ = 0.5, 1, 10,
10⁶. T2's count ranges over 200–2,200 (truth 1,200). Without gDNA, isoform expression and exon-specific capture
are not distinguishable, and the answer is whatever the γ prior says.

**What identifies what.**

| Claim | Supplied by |
|---|---|
| gDNA vs RNA at an object, stranded | likelihood (strand split) |
| gDNA vs RNA, unstranded | nothing local, except matched-footprint length laws when they differ (unbuilt in calibration) |
| Capture contrast where gDNA is present and separable | likelihood, via the physical assumption that gDNA is uniform before capture |
| Capture of a piece shorter than gDNA fragments | not its own data; its neighbours' crossings (a footprint assumption) |
| Capture contrast at zero gDNA, and the isoform split under it | prior only (the ridge) |
| DNA → RNA capture transfer | assumption; violated by overlap-length binding under gaps |
| Landscape, length laws, κ, grouped gDNA prior | the same library (empirical Bayes) |

The grouped prior's strength is one pseudo-fragment per gDNA-eligible unit, and DESIGN §3.1c says plainly that
it is not derived ([pipeline.py:625–656](../../src/rigel/pipeline.py)). It adds one pseudo-fragment of
calibration's share for every gDNA-eligible fragment the EM sees. That is the largest evidence reuse in the pipeline. It is also why in-scope gDNA pools are accurate,
so keep it for 0.8.0 and disclose it. The shared field's "one factor per fragment" removes none of these.

**What happens at zero gDNA under capture.** With the shipped geometry, nothing recovers the class scale. Ladder
g00 stranded ON span/isoform L/Y is 5.68 flat, 14.54 shipped and 17.84 own-evidence (`footprint_screen.out`).
End to end, own (6.61 %) ≈ the October 5 tree, which was flat there (6.47 %), ≪ shipped (15.10 %). This is an incompatibility, not a defect to engineer
away. "Capture-aware isoform quantification at zero gDNA, with no panel input" is not identified. The only
defensible default is the prior's centre: no differential capture. The alternative is a probe/target BED as an
input, which is a structural observation, not a detector. That is outside 0.8.0's interface.

**Conditioning (the prompt's question).** Yes, three times:

- **On a fragment's footprint.** `π` cancels among candidates that share it. That is why most production
  numerators may be dropped, and why a numerator is needed exactly where candidate paths differ.
- **On a locus's gDNA total.** Object gDNA counts become `Multinomial(G_L; E_o c_o / Σ_j E_j c_j)`. The gDNA
  level, and with it the library landscape and the median unit, drops out. The information about capture
  contrasts is the gDNA count itself, so G_L = 0 gives exactly the prior. This is the principled form of
  "weights need gDNA behind them" (§9).
- **On matched-footprint length composition.** This identifies gDNA:RNA where the laws differ. It is the
  deferred stratum's missing witness, and it is unbuilt.

## 4. Ranked findings (decision 3)

| # | Finding | Class | Evidence | Refuted by |
|---|---|---|---|---|
| 1 | The reader reuses count-solve **lower-bound** messages as factors for a **ratio** | model (estimand) | g00 str OFF: 696 objects with log w > 1 → 53 without the factors; 13.61 → 2.79 % | a factor-invariant reader that keeps the regression |
| 2 | Capture contrast is unidentified without gDNA | non-identifiability | the ridge; g00 ON class scale is wrong for every reader | an estimator recovering it without gDNA or a panel |
| 3 | gDNA kernel and normaliser use different length frames | code/model consistency | Σ K/L = 1.071/1.009/1.024; bundled in the g50 gain | the factorial shows no effect |
| 4 | No capture numerator where candidate paths differ | model approximation | multi-cut incidence 1.03–1.82; bundled in the g50 gain | the factorial shows no effect |
| 5 | DNA → RNA capture transfer is biased under length gaps | model limitation | ruler mirror image 1.96 vs 0.45 | per-fragment K/L matches the oracle |
| 6 | The additive junction rule can zero captured junctions | code fragility | 2.6 % zero cuts (shipped), 0 % (own) | zero cuts carry no spliced reads |
| 7 | The shared field regresses at g00 | model (ridge + prior) and inference (start dependence 1,965 fragments; bridge overlap ~1e-37) | test g00, all four rows | — |
| 8 | Free junctions add unidentified coordinates | model | 40.56 % str ON | — |
| 9 | Prototype gains crossed solvers | inference confound | shipped MAP 11.21 vs VB 9.24 (test g50 str OFF) | — |
| 10 | Residual own-strand false modes | model (binomial strand noise, od = 0 ruling) | 53 objects; 2.79 vs flat 1.93 % | a real gDNA-free stranded library shows none |
| 11 | Deferred stratum | non-identifiability + removed constraint + unbuilt channel | 62.96 % ladder | out of scope |

Finding 1 in one line: DESIGN.md's level-lane ruling (2026-09-05) says a level crossing a face means "at least
this much gDNA". The reader multiplies the intersection of the sides' levels into each object's posterior
([solve_kernel.cpp:668–680](../../src/rigel/native/solve_kernel.cpp),
[honest_reader.h:167–184](../../src/rigel/native/honest_reader.h)). When the object holds no gDNA, its own
likelihood is one-sided below an upper bound `U`. A false neighbour bound at `x_L < U` then sets the mode to at
least `x_L`, and the weight to at least `exp(x_L − median)`. That is unbounded as the library's gDNA → 0, even
though behind `x_L` stand a few strand-noise fragments. For counting the error is those fragments; as a capture
contrast it multiplies every RNA opportunity at the object.

## 5. The options (decision 4)

| Option | Promise | Discriminating experiment | Cost | Verdict |
|---|---|---|---|---|
| Own-evidence reader | removes the flat-field blocker | done: full matrix and VCaP; remains a production port that reproduces the arm's arrays exactly | port ≈ 1 day; matrix ≈ 45 min | **back, next** |
| Neighbour level admitted only where the object cannot see its own bases | recover RNA-short str ON without the flat-field cost | offline ruler on a derived visibility weight (seconds), then test chromosome (≈ 2 min), then matrix | ≈ 1 h | second, if the owner wants ≥ main on short captured exons |
| Numerator/normalizer repair | test g50 −2.4 to −3.3 points against shipped VB, genes +0.03 to +0.18 | 2³ factorial (path numerators × gDNA frame × membership), shipped VB, shipped field, test g00 and g50 | ≈ 15 min | third |
| Restore the Oct 5 reader / the factory | ladder 1.93 / 6.47 % | already measured: RNA-short str OFF 42.30 %, test str ON 13.89 %; detector ruling | — | refuse |
| Shared field | — | already measured | — | abandon |
| Locus-conditional contrast reader | flat at zero gDNA by construction; also suppresses finding 10 | enumerable fixture, then offline ruler, then matrix | days | document only (§9) |

## 6. The single next mechanism: the own-evidence reader

**Derivation.** Read each object's capture from its own observations once:
`x*_o = argmax_x [ℓ_o(x) + log P(x)]`, where `ℓ_o` is its strand-column likelihood with the RNA amount
integrated out. Publish `w_o = exp(x*_o − median x*)`. The delivered composition row and held levels stay in the
count solve, where they were designed and earn their keep. They leave the capture estimand, where finding 1 makes
their error unbounded. There is no new constant, prior or detector, and the landscape and median unit are as
shipped.

**What it buys and costs, measured.** The table in §1. Read-only localisation from the archived weights
(`false_modes.out`, `footprint_screen.out`):

- At ladder g00 stranded OFF, 543 of the shipped reader's 696 false weights sit on objects gDNA **can** see (eff_g ≥ 50;
  boundaries median depth 5), and own evidence removes 643 of the 696.
- At RNA-short g50 stranded ON, own evidence loses enrichment on **dim** regions (gDNA opportunity < 5 but
  positive, median 101 bp): 5,829 enriched regions become 956. Ruler within-gene log(L/Y) SD goes 0.095 → 0.137. On visible regions
  the two readers agree (2,008 vs 1,995).
- A diagnostic hybrid (shipped weights only on regions with gDNA opportunity < 5, own evidence elsewhere) recovers
  0.099 within-gene and span/isoform 1.16 on RNA-short (shipped 0.095 / 1.14). It also keeps ladder g00 stranded OFF at
  0.418, against 0.963 shipped. The useful neighbour information and the harmful one are separable by the object's own gDNA
  visibility. The cut is a localisation, not a proposal.
- Naive spatial pooling is refuted. Opportunity-weighted averaging of each region with its boundaries recovers
  within-gene ratios (0.137 → 0.105) but inflates span/isoform scale on every captured row (ladder g50 1.15 → 1.45,
  RNA-short 1.27 → 1.60). Enrichment contrasts of e⁵ let a 5 % share dominate a long intron.

**Falsifiers before `src/`.**

- **`port-identity`.** All six count arrays bit-identical to main (the reader touches no count), and weights
  identical to the archived `unary` arm wherever it saved them.
- **`delivered-row-invariance`.** Perturbing any delivered row changes no weight. The gate must fire on main.
- **`flat-field-fixture`.** On a zero-gDNA stranded chain with one false neighbour level, the weight is 1
  under own evidence. Re-admitting the level must fire the gate.
- **`release-table`.** Every in-scope stratum ≥ tag at gene level; transcripts ≥ tag except the two parity
  strata, reported as such; controls, deferred rows and VCaP beside them.
- **`full-suite`.** `preflight.py --full`, with the golden diff read and recorded before
  `--update-golden`.

**Why it earns its simplicity.** It deletes about 70 of `honest_reader.h`'s 575 lines and the reader block's
factor assembly ([solve_kernel.cpp:639–717](../../src/rigel/native/solve_kernel.cpp)). It also cuts the reader's
dependence on the message layer. That allows a later, numerically neutral refactor that moves the reader out of
the block solve into a per-object pass after the last refit.

## 7. Keep / delete / replace

| Mechanism | Decision | What it buys |
|---|---|---|
| Reader: own likelihood, landscape, mode, median unit | keep | the class scale under capture: test g50 str ON flat 32.11 % vs 6.58 % (fixed field); span/isoform 5.6 → 1.1–1.3 |
| Reader: delivered composition/DNA/RNA factors | **delete** | buy 6.2 points on RNA-short str ON (20.43 → 14.24 %); cost 10.8 on ladder g00 str OFF (2.79 → 13.61 %) and 1.35 on test str OFF contaminated (8.00 → 9.35 %) |
| Shared field, rate model, free junction, frozen M-step, transports, sampler | **delete from the release path** (archive) | nothing attributable to integration |
| Path export in the scorer (prototype patch) | hold | needed only if option 3 earns it |
| gDNA frame and path numerators | replace candidate (third) | unattributed g50 gain |
| Conserved-share length rule | keep | exact bookkeeping; first-order physics stated |
| Junction rule | keep as comparator, flagged (finding 6) | — |
| Count-solve level lane | keep | it serves counting; only its reuse is the defect |
| Grouped gDNA prior (strength = eligible units) | keep, disclose | in-scope pool accuracy |
| VB EM | keep | 0.4–2.1 points better than MAP at test g50 |
| Multimapper second pass | keep | VCaP 0.2401 → 0.2528 |
| Factory / deferred mechanisms | defer | out of scope |

## 8. Production boundary (decision 5)

**Capture reader.**

- **Inputs.** Per object: strand columns (u, v), gDNA opportunity, RNA opportunity, strand admissibility, κ; the
  last refit's landscape grid and log-density. Nothing from neighbours.
- **Output.** Per region and boundary: `w = exp(mode − median mode)`, and 1 where the reader did not read.
- **Assumptions, stated in the manual.** (i) gDNA is uniform before capture. (ii) DNA and RNA share capture per
  object; biased under fragment-length gaps. (iii) Pieces gDNA cannot see read the landscape's typical level;
  short probed exons are under-contracted. (iv) Without gDNA, within-locus capture contrast is not identified,
  and the weights stay near 1 except where own strand counts exceed the RNA leak.
- **Tests.** The existing curve-by-curve and mode-by-mode spec gates, trimmed to the own term, plus `delivered-row-invariance` and `flat-field-fixture`.
  Change the Python reference `tests/native/_honest_reader_reference.py` and the native code together.

## 9. Refused or deferred, and why

- **Shared field, any variant.** §0 and §3. Not salvageable by a better sampler: the target is not identified
  where it was meant to help, and its prior is not a capture prior.
- **Posterior averages (E[eˣ], E[x]), adjacency-coupled readers, opportunity-weighted footprints.** Each is
  measured and moves locations or spreads exon capture into introns.
- **The Oct 5 located-mode reader, or any detector or library reference.** Ruling, and RNA-short 42.30 %.
- **Flat weights in production.** Captured rows collapse (test g50 str ON 32.11 %, span/isoform 5.6).
- **The locus-conditional contrast reader** (§3, second bullet). It is the principled successor if finding 10
  matters on real data. It needs a new contrast prior and calibration's per-object gDNA posteriors as soft
  counts, a reuse to price. It is not for 0.8.0 without evidence from a real gDNA-free stranded library.

## 10. Plan, budgets, stop conditions

0. **Done, read-only, under 5 min.** The release-bar table and the localisation in §6.
1. **Port in a worktree after main is pushed.** Delete the factors from the native reader and its Python
   reference, then run `port-identity`, `delivered-row-invariance` and `flat-field-fixture`, and `pytest tests/native tests/calibration` (≈ 2 min). **Stop** if
   `port-identity` is not exact.
2. **The candidate arm, pinned and fractional, seed 0.** Ladder (≈ 11–15 min), both gap panels (≈ 6), test
   chromosome (≈ 1.5), single-probe and junction layouts (≈ 2), then VCaP alone (≈ 5–10). That is about 45 min.
   Main's arm is runtime-identical to the archived matrix, so it is reused. **Stop** if any in-scope stratum
   falls below the tag beyond the two known parity rows, or if VCaP's gDNA fraction moves more than 0.002.
3. **Suite and gates.** Full suite and `preflight.py --full` (≈ 10 min). Read the goldens' diff before updating.
   The owner commits.
4. **Owner decision.** Is tag parity on RNA-short stranded ON (20.43 vs 20.33 %) acceptable for 0.8.0? If not,
   go to step 5. If yes, release work moves to the manual's stated limitations.
5. **Visibility-weighted neighbour level.** Derive the visible fraction of a piece's gDNA observations from
   geometry: contained over overlapping placements under the gDNA law, with no constant. Then the offline ruler
   (seconds; reuse `footprint_screen.py`), then the test chromosome (≈ 2 min), then the matrix (≈ 45 min).
   **Stop** if span/isoform moves on any captured row or g00 stranded OFF rises above own-evidence's 2.79 %.
6. **The numerator/normalizer factorial.** 8 arms × 8 test conditions, shipped VB, shipped field (≈ 15 min).
   Promote only a factor that improves g00 and g50 rows without a gene regression, then the matrix.

## 11. Demonstrated, plausible, missing

- **Demonstrated.**
  - Main misses the tag on two in-scope strata, both from false capture on flat stranded fields.
  - The own-evidence reader clears both and is at parity or better elsewhere.
  - The shipped reader's flat-field false weights come mostly from delivered factors at gDNA-visible objects;
    their useful information sits at gDNA-blind objects.
  - The shared field loses to its own fixed-field control.
  - The zero-gDNA ridge is exact.
  - The gDNA kernel and normaliser disagree in frame.
  - Naive footprint pooling distorts class scale.
- **Plausible.**
  - Population-specific capture under length gaps (ruler mirror image; needs the K/L audit).
  - The g50 geometry gain comes from the gDNA frame and path numerators (needs the factorial).
  - Visibility-weighted neighbour levels recover RNA-short without the flat-field cost (ruler-level only).
- **Missing.**
  - A production port of the own-evidence reader and its receipt.
  - A real gDNA-free stranded library.
  - Per-fragment K/L against the oracle on the gap panels.
  - Plasma-library weight distributions under own evidence.
  - A transcript-truth real library.

## Appendix — this review's measurements

`.cache/rigel_runs/2026-10-12_independent_review/` holds the following. Each script is read-only over the
archived matrix (`2026-10-11_implementation_review`), the October 5 tagged sweep, the oracle caches and the
cached ruler truth. None reruns calibration or the EM.

| File | What it shows |
|---|---|
| `bar_table.py` / `.out` | per-rung and per-stratum tag, Oct 5, main and readers |
| `false_modes.py` / `.out` | where high weights sit, by kind, visibility and depth |
| `junction_dim.py` / `.out` | per-piece and junction differences between readers |
| `footprint_screen.py` / `.out` | the offline ruler: within-gene SD and span/isoform for shipped, own, flat, footprint and the diagnostic hybrid |
| `ridge_check.py` / `.out` | the exact ridge |

The screen approximates the RNA length law by the cached spliced pool. It reproduces the archived ruler (0.095 vs
0.096 shipped, 0.137 vs 0.141 own on RNA-short stranded ON) and reconstructs shipped lengths with median error
1e-3. It is a length-truth instrument; it does not stand in for the transcript table.
