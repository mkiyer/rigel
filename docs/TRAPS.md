# TRAPS — mistakes this project has already made

**Read this before designing anything.** Every entry is a real defect that shipped, or a real
measurement that was believed and was wrong. They are here because each one cost time, and several were
made twice.

⛔⛔⛔ **EVERY RULE HAS A NAME. NUMBERED LABELS ARE BANNED, AND A TEST ENFORCES IT** (owner, 2026-08-07;
gate: `tests/test_no_jargon_labels.py`). This file used to call its rules `A16`, `D4j`, `C0b`. Cite them as
**`TRAPS: off-grid-message-mode`**, **`TRAPS: a-cancelling-defect-pair`**,
**`TRAPS: frame-free-is-not-assumption-free`** — the name IS the identifier, so a citation says what it
means without a lookup, and it is still one greppable string with one home.

⚠ **The numbers were not merely opaque, they were AMBIGUOUS, and that is what forced this.** Measured
across the tree, `A1` meant a validation trap here *and* the FIDELITY criterion in `SUCCESS.md`; `C1` meant
a pool trap here *and* a moment variable in the length quadrature — then in `length_likelihood.py`, a module
PURGED at `b7ed7a0b`, the moments being mis-filed rather than dead: they live in `effective_length.py`,
which is where `tests/test_no_jargon_labels.py` scopes the `C1` exemption today; and **`G1` meant "no magic numbers"
here and "a structurally pure-gDNA object" in the solver — 201 occurrences, and a reader could not tell
which.** That is this file's own `TRAPS: two-masks-one-name` committed by the labelling scheme.

⭐ **The rename touched 980 citations across 114 files and 66 occurrences were correctly NOT touched**,
because they were the other meanings above. A name cannot collide that way: nothing else in the tree is
called `no-magic-numbers`.

⛔⛔ **THE SHAPE OF AN ENTRY, AND IT IS ENFORCED BY PRUNING.** A trap is **the mistake, the tell, and the
rule** — three or four lines. The investigation that produced it is in git; the number that describes
current state is in `ROADMAP.md`. ⚠ This file reached 741 lines and 91 entries before its first prune
(2026-08-05) because entries were being *appended to* as new instances arrived. ⭐ **Instances of ONE lesson
belong in a LIST inside that lesson, not as new entries** — see `TRAPS: a-message-from-the-destinations-belief`, which is nine costumes of one mistake.
⛔ A rule's NAME is cited from source docstrings and tests, so a name is never changed once written; a
merged rule survives as a row in its family's table, under its own name.

⛔ **This file holds LESSONS, not measurements.** A number appears only where it *is* the lesson (an
"exact factor of 2", a "12× over-represented"). Current measured state lives in `ROADMAP.md`; anything
here that reads like a status report is a bug in this file.

⚠ **AND A MEASUREMENT IS NEVER DELETED — ONLY STAMPED.** Where a number here was measured on the
**36-condition** gDNA ladder it predates that panel's rebuild to **16 conditions** (2026-08-13); it is kept
as a historical record rather than restated, because a lesson does not expire with its panel. Anything you
need as *current* comes from `ROADMAP.md`, never from a number in this file.
⛔ **The same stamp covers the `flgap_short` / `flgap_long` pair, which was DELETED in that same rebuild.**
Every flgap number below is a historical record and **cannot be re-run by anything on disk**: the panels
went on 2026-08-13, and their configs, their study cache and the two instruments that read them went on
2026-08-17. ⚠ Re-deriving one means DESIGNING a length-gap panel, not restoring a deleted file — and
`EQUATIONS.md` §3b is the argument any such design has to answer first.

---

## THE INDEX — every rule in one line, so you can scan instead of read

⭐ **144 rules, and every one has exactly one body — RE-DERIVED, never carried: the header said 105 while
the file held 115, the same way it once said 99 while the file held 101. Read the group that matches what
you are about to do, then open only those.**
⭐⭐ **RE-DERIVED 2026-08-19 AND THE INDEX IS COMPLETE: 144 entries, 144 bodies, and BOTH set differences
print empty.** (141/141 before Stage 2 of the message-propagation rebuild added three — one each to §B, §E, §F;
the header had read 138 while the file held 141, the drift this line exists to catch.) (137/137 before the relay-licence dissection added `a-scorer-scoped-to-the-mechanisms-targets`.) (124/124 before the composition-reference work added five to §D; 129/129 before ψ's
simplex-closure work added three more; 132/132 before the doc audit below added one; 133/133 before the
message-passing lock-down added four — two to §D and two to §A.)
⛔⛔ **AND THE +1 IS ITS OWN LESSON: `a-green-suite-hid-five-dead-instruments` HAD A NAME, FOUR CITATIONS
AND NO BODY.** `CLAUDE.md` and three instruments (`vertex_ceiling.py`, `toy_harness.py`,
`toy_trace_error.py`) cited it as `TRAPS: a-green-suite-hid-five-dead-instruments` while its text lived
only in `CLAUDE.md`'s prose — so the rule had two homes, one of which was not this file, and the citation
led nowhere. ⭐ **The set-difference derivation above cannot see that**, because it compares this file
against itself; a dangling citation is found only by grepping `TRAPS: ` across the tree and checking each
name resolves here.
⛔⛔ **The "120 index entries against 125 rule BODIES" this paragraph used to record was a COUNTING
ARTEFACT, not a gap — and that is the sharper lesson.** The body count came from matching `**<token>.` at
the start of a line, which also matches a bolded NUMBER where prose wraps: `**64.5 % of transcript
error…`, `**0.00 / 0.51 / 0.98**`, `**146.05**`, `**1.76–2.20× WORSE…`, `**3.4e9× finer**`. Exactly five,
and they were read as five missing rules. ⭐ **A rule name always begins with a LETTER**, so the
reproducible derivation is::

    index   ^- `([a-z0-9-]+)`          in the section above
    bodies  ^\*\*([a-z][a-z0-9-]*)\.\s in the sections below

⚠ **This is `TRAPS: re-record-the-baseline` committed one level up: the number was re-derived faithfully
and the METHOD was never validated, so a re-derivation reproduced the same wrong answer.** Re-derive with
the pattern above and check BOTH directions print empty; do not adjust the number. Cite a rule as
`TRAPS: <name>` — the name is the identifier and `tests/test_no_jargon_labels.py` keeps it that way.

**Validation and gates**

- `self-checking-validator` — A validator that calls the builder's own helper validates nothing.
- `perturb-every-gate` — A falsification test needs PERTURBATION, not just to be written first.
- `a-field-driven-gate-is-atomic` — A gate that enumerates the spec's fields fuses spec, implementation and schema into ONE commit.
- `waive-with-a-measurement` — AN ASSERTION THE SHIPPED CODE VIOLATES IS WAIVED WITH ITS MEASUREMENT, NEVER WIDENED — AND
- `a-docstring-that-misdescribes-the-graph` — A CLAIM ABOUT THE CODE, INSIDE THE CODE, THAT NOTHING
- `a-flat-pile-is-not-a-knot` — 35 MODULES WITH NO CYCLES AND ONE IMPORTER EACH IS NOT TANGLED — IT IS
- `the-rename-that-corrupted-a-diagram` — A MECHANICAL REWRITE OVER PROSE WILL HIT TOKENS THAT ARE NOT
- `a-gate-that-already-passed` — A gate that already passes is not a falsification.
- `right-conditional-wrong-marginal` — When a per-observation conditional is correct, no conditional comparison can detect a missing
- `byte-identity-gate` — A bit-identity gate has lied in both directions.
- `the-deliverable-is-not-reproducible-by-default` — THE SHIPPED PIPELINE DOES NOT REPRODUCE ITSELF, SO AN END-TO-END A/B ON
- `a-clip-hides-a-scale-error` — A min() clip hid an exact factor of 2 for months.
- `an-inverted-clip-is-a-constant` — a two-endpoint clip is silently a CONSTANT both when lo > hi and when
  lo == hi — and at a zero control that constant is the right answer, testing nothing.
- `an-ablation-that-never-ran` — AN ABLATION THAT NEVER RAN READS AS "NO EFFECT", AND TWO IMPORT HABITS MAKE THAT EASY.
- `an-inert-arm-reads-as-a-refutation` — ⭐⭐ AN INERT ABLATION READS "NO EFFECT"; AN INERT PROTOTYPE READS
  "YOUR IDEA DOES NOT WORK", AND THAT KILLS A DESIGN DIRECTION. Assert the patch target, count the firings,
  BEFORE any number is read.
- `a-green-suite-hid-five-dead-instruments` — AN INSTRUMENT DIES WHEN THE TREE MOVES UNDER IT, AND A GREEN
  SUITE SAYS NOTHING. Two mechanisms so far: a `src/` DELETION, and a CONFIG DEFAULT FLIP.
- `a-zero-count-is-a-measurement` — A ZERO COUNT IS A MEASUREMENT OF A DENSITY, NOT AN ABSENCE OF DATA — AND KEYING PRECISION
- `a-ratio-cannot-carry-zero` — THE RELAY'S CURRENCY IS A DENSITY RATIO, AND ZERO IS NOT IN THE MULTIPLICATIVE GROUP —
- `the-divergence-was-a-barrier` — THE DIVERGENCE WAS DOING A SECOND JOB, AND REMOVING IT TURNED A RELAY BARRIER INTO A
- `deadband-from-the-wrong-sample` — A NOISE DEADBAND WHOSE CUSHION IS SUPPLIED BY AN UNRELATED SAMPLE SIZE FAILS EXACTLY WHERE
- `honesty-metrics-reward-ignorance` — EVERY HONESTY METRIC IMPROVES AS THE SOLVER STOPS KNOWING ANYTHING — SO NONE IS READABLE
- `predicate-contradicts-its-docstring` — A PREDICATE CAN CONTRADICT ITS OWN DOCSTRING FOR MONTHS IF SOMETHING ELSE IS MASKING IT —
- `a-test-that-redefines` — A TEST THAT RE-DERIVES A DEFINITION CANNOT DETECT DRIFT IN IT.
- `a-gates-power-is-its-invariant-set` — A COMPARISON GATE'S POWER IS THE SIZE OF ITS INVARIANT SET; a change that
  shrinks that set disarms the gate, so it must not ride along with the change being gated.
- `a-gate-that-reconstructs` — A GATE THAT RECONSTRUCTS A VALUE IS VACUOUS WHEREVER THE CODE DECIDES NOT TO USE IT — GATE A
- `off-grid-message-mode` — A MESSAGE MODE OUTSIDE ITS GRID'S DOMAIN IS NOT A WEAK CLAIM, IT IS A PIN AT THE BOUNDARY —
- `a-comment-quoted-as-a-finding` — A WORD BORROWED FROM A COMMENT IS NOT A MEASUREMENT, AND ONE PROPAGATED INTO THE DOCS INVENTED
- `the-intermediate-is-not-the-deliverable` — AN INTERMEDIATE METRIC IS NOT THE DELIVERABLE UNTIL SOMEBODY MEASURES THE COUPLING — AND HERE
- `could-the-arm-have-fired` — "THE ARM CHANGED NOTHING" IS NOT A CONTROL UNTIL YOU CHECK THE ARM COULD HAVE CHANGED
- `prove-the-substrate` — Prove the SUBSTRATE before proving the code.
- `can-the-benchmark-resolve-it` — Before running a benchmark, prove it can resolve the axis you are changing.
- `toys-rank-hotspots-backwards` — A toy ranks performance hotspots BACKWARDS.

**Measurement and inference**

- `measure-the-ceiling-first` — MEASURE THE CEILING BEFORE BUILDING THE CORRECTION.
- `a-broad-population-carries-no-prior` — a population prior transfers information only if the pooled
  distribution is TIGHTER than the per-object uncertainty. RNA's is 3 decades wide: worth nothing.
- `attribution-must-survive-a-shuffle` — an oracle arm a SCRAMBLED oracle also wins has not priced
  the truth. Permute the truth you fed it and check the gain disappears.
- `score-against-truth` — Score against TRUTH, not against the previous run.
- `zero-target-guards-are-one-sided` — A zero-target guard is ONE-SIDED.
- `hard-labels-miss-soft-change` — Hard-label metrics are nearly blind to soft changes.
- `never-pool-the-strata` — THE DEFAULT INSTINCT IS A POOLED AVERAGE, AND IT IS WRONG HERE THREE WAYS.
- `a-scorer-scoped-to-the-mechanisms-targets` — A SCORER RESTRICTED TO THE POPULATION A MECHANISM TARGETS WILL REPORT A WIN
  THE WHOLE CHAIN REFUSES; the sign can flip. The scoped score is a diagnostic, never a verdict.
- `zero-the-precision-with-the-value` — A REFUSED CLAIM MUST LOSE ITS PRECISION IN THE SAME STATEMENT THAT ZEROES ITS
  VALUE, or a later operator hands the precision back and the "no claim" becomes a confident ZERO — which the
  mass pin then rescales into "all gDNA".
- `a-transcript-predicate-must-not-silently-drop-a-molecule` — A FRAGMENT REJECTED BY A TRANSCRIPT-LEVEL
  PREDICATE DEPOSITS NOTHING INTO CALIBRATION, AND THE REJECTED POPULATION IS NEVER RANDOM. Make the
  fragment ledger CLOSE.
- `a-ratio-needs-a-population-that-can-supply-its-numerator` — ADDING OPPORTUNITY THAT STRUCTURALLY CANNOT
  CARRY COUNTS MANUFACTURES A DEFICIT OF EXACTLY ITS OWN SHARE. Score over the population that can supply
  the numerator, which is why the field gate splits by reference.
- `a-threshold-on-a-fitted-residue` — A BINARY CUT ON A FITTED PARAMETER'S RESIDUE IS NOT A POPULATION TEST — AND A BETTER THRESHOLD I
- `excluding-a-population-hides-it` — EXCLUDING A POPULATION FROM THE DENOMINATOR WITHOUT GATING ITS OWN FAILURE MODE HIDES THE
- `name-the-observable-per-site` — IF EVERY GATE FOR A CHANGE READS ONE INSTRUMENT, THE CHANGE IS GATED IN ONE PLACE — AND
- `starved-is-not-depleted` — "STARVED" AND "DEPLETED" ARE NOT THE SAME DIAGNOSIS.
- `the-substrate-knob-fought-back` — A TILING LOOP CAN SUPPRESS THE VERY POPULATION UNDER STUDY.
- `key-on-a-realised-quantity` — AN RNA "LEVEL" THAT IS A MULTIPLE OF THE gDNA DENSITY IS NOT COMPARABLE ACROSS THE CAPTURE
- `price-the-halves-separately` — PRICE THE HALVES SEPARATELY — A DIAGNOSED DEFECT AND A VALUABLE ONE ARE NOT THE SAME DEFECT.
- `panel-before-src` — A TOY CEILING IS NOT A PANEL CEILING, AND THAT HAS NOW COST THREE BUILDS.
- `capture-inverts-the-counted-side` — "THE WELL-COUNTED SIDE" IS NOT A FIXED SIDE — CAPTURE INVERTS IT.
- `admitting-an-object-costs` — ADMITTING AN OBJECT TO THE SCORED POPULATION IS A COST, AND A MECHANISM CAN DO IT SILENTLY.
- `substitution-understates-a-source` — A CEILING BY SUBSTITUTION UNDERSTATES A MESSAGE SOURCE.
- `a-symptom-is-not-a-second-defect` — A WRONG NUMBER FED A WRONG INPUT IS A SYMPTOM. Substitute ONLY the
  upstream input, re-run the SHIPPED downstream function — and check the ceiling's injection point.
- `a-locked-object-is-not-a-control` — A STRUCTURALLY LOCKED OBJECT IS NOT A CONTROL.
- `draining-breaks-the-oracle` — A per-fragment-independent partition stops being one the moment a downstream step conditions on
- `an-equal-length-panel-defeats-the-lift` — THE REPAIR FOR THE RULE ABOVE HAS A SUBSTRATE IT CANNOT WORK ON, AND IT IS THE
- `a-length-gap-bypasses-calibration` — A LARGE gDNA-vs-RNA LENGTH GAP LETS THE EM ANSWER ON LENGTH ALONE
  AND MASKS CALIBRATION BUGS. The equal-length panel is what exposes them.
- `weight-it-like-the-consumer` — A bp-weighted mean and a fragment-weighted mean answer different questions.
- `a-support-ceiling-is-the-clamp` — A SUPPORT CEILING THAT MATCHES THE CLAMP IS NOT A MATCH — it is the clamp.
- `log-variance-is-not-linear` — Every "is the declared precision earned?" number written before 2026-07-28 compared a LOG-space
- `re-record-the-baseline` — A delta is only attributable if its baseline came from the same tree in the same session.
- `score-the-consumers-own-count` — TRUTH IS NOT A YARDSTICK UNTIL YOU NAME THE CONSUMER — a truth nothing reads measures your choice of target.
- `a-truth-table-of-aggregates` — A TRUTH TABLE'S LABEL COLUMN MAY HOLD NESTED AGGREGATES, AND SUMMING IT DOUBLE-COUNTS.
- `a-moment-match-is-not-sufficient` — MATCHING THE MODEL'S MOMENTS TO THE POPULATION'S CAN MAKE THE ANSWER WORSE.
- `a-single-level-panel-cannot-see-a-constant` — A PANEL THAT HOLDS THE ANSWER FIXED CANNOT TELL A GOOD ESTIMATOR FROM A CONSTANT.
- `the-floor-must-reproduce-the-selection` — A NOISE FLOOR DRAWN UNCONDITIONALLY UNDER A SCORER THAT CONDITIONS
  (on a non-empty destination) READS THE TRUNCATION AS A RULE ERROR — 10 % at ~1 expected count per object.

**Pools, selections and divisors**

- `two-divisors-opposite-sign` — TWO DIVISORS BUILT FROM ONE pmf CAN STILL DISAGREE — IF THEY RESPOND TO IT WITH OPPOSITE SIGN.
- `frame-free-is-not-assumption-free` — A TERM THAT IS FRAME-FREE IS NOT THEREFORE ASSUMPTION-FREE — LOOK AT THE NUISANCE YOU
- `a-purity-filter-is-a-length-filter` — A PURITY FILTER ON A LENGTH POOL IS A LENGTH FILTER.
- `pure-and-length-censored` — A POOL CAN BE COMPOSITION-PURE AND LENGTH-CENSORED AT THE SAME TIME, and the second is invisible
- `divide-by-a-probability` — A POOL DIVIDED BY ITS OPPORTUNITY MUST BE DIVIDED BY A PROBABILITY, NOT BY A COUNT.
- `opposite-tilts-must-not-pool` — POOLS WITH OPPOSITE TILTS MUST NOT BE POOLED RAW.
- `a-mean-of-ratios-inherits-the-partition` — AN ESTIMATOR THAT AVERAGES PER-OBJECT RATIOS IS A FUNCTION OF WHERE THE
  ANNOTATION DREW ITS BOUNDARIES.
- `a-trap-names-the-defect-not-the-repair` — CITING A TRAP DOES NOT LICENSE A MECHANISM. CHECK THE GRAVEYARD FOR THE
  MECHANISM ITSELF.
- `a-stale-gate-accuses-the-newest-change` — A GATE WHOSE PREMISE EXPIRED FAILS AND BLAMES WHATEVER IS IN
  FLIGHT. When a gate fails on something you did not touch, first ask whether it can still pass at all.
- `an-upper-bound-is-not-an-estimate` — A BOUND THAT IS SHARED IS NOT EVIDENCE ABOUT WHO SHARES IT. If the
  estimator's SUPPORT does not move as you vary it, you are varying the wrong thing.
- `a-gate-on-the-helper-is-not-a-gate-on-the-caller` — TESTING THE FUNCTION THAT COMPUTES A QUANTITY DOES NOT
  TEST THE CALL THAT ASKS FOR IT. Two holes, 29 green gates, both found only by perturbation.
- `fractional-mass-is-the-problem` — Fractional mass IS the partitioning problem.
- `conservation-misses-mis-attribution` — Mass conservation does not catch mis-attribution.
- `a-guard-outlives-its-divisor` — DELETE THE DIVISOR AND THE GUARD AGAINST IT GOES INERT — WHILE ITS TEST KEEPS PASSING.
- `a-fold-grows-a-heuristic` — A QUANTITY FOLDED ONTO ANOTHER AXIS TO FIT A CONSUMER'S INTERFACE WILL GROW A HEURISTIC TO REPAIR THE FOLD.

**Estimation and solver design**

- `we-keep-re-deriving-message-passing` — ⭐⭐⭐ **READ THIS ONE FIRST.** THREE OR FOUR sessions have
  re-derived message passing from scratch. *The tell:* you are reasoning about how an EXON gets its gDNA level and
  have not yet said the words. An exon's level cannot be measured — only imputed from neighbours — and
  imputation across the chain IS message passing. Carries the four STALE beliefs that keep sending
  sessions back to the start, and a pointer to the one home of the derivation.
- `one-hop-lifted-out-is-still-the-relay` — ⭐⭐ Lifting one hop out of `messages/` because it carries a
  MEASUREMENT rather than a BELIEF does not make it something else. A private copy of a framework you
  switched off is not a way to avoid the switch — fix the framework's known bug and re-price.
- `a-variance-cannot-fix-a-bias` — You cannot fix a biased mode with a variance.
- `two-gaussians-one-latent` — Never hand a solver two Gaussians built from one latent.
- `variance-fitted-on-the-belief` — Never fit a variance on the current, not-yet-solved belief.
- `a-message-from-the-destinations-belief` — A MESSAGE COMPUTED FROM THE DESTINATION'S OWN BELIEF CARRIES ZERO INFORMATION AND CONFIRMS THE
- `no-prior-means-haldane` — "No prior" does not exist on a grid — omitting a term lets the grid supply Haldane
- `prefer-shares-to-differences` — Sums are well conditioned; differences are not.
- `an-all-zero-factor-is-inert` — An all-zero factor is uninformative, not decisive.
- `density-below-one-fragment-length` — Density below one fragment length is not resolvable by ANY design.
- `identical-paralogs-are-bimodal` — An identical-paralog split is bimodal, and depth does not fix it.
- `compatibility-is-geometry-not-composition` — A PREDICATE BUILT FROM WHAT A FRAGMENT IS COMPATIBLE WITH
  CANNOT SEPARATE RNA FROM gDNA.
- `a-priors-curvature-is-not-the-datas-information` — A PRIOR MAY NOT CONTRIBUTE TO A FISHER PRECISION.
  Feeding one in is a gate flip that releases the full count precision, and it credits data-free slots.
- `deriving-one-coordinate-propagates-its-error` — closing a composition by deriving B from A makes B
  inherit A's defects, including ones B's own estimator did not have.
- `interpolate-on-the-axis-where-the-lattice-is-uniform` — a quantile interpolated on a non-uniform
  transform of the grid is biased toward the middle, by 2.7e-03 at n_grid 60.
- `read-the-whole-failure-list` — "21 golden failures" read off the last eleven lines; two were real.
  Derive the failure set with grep, never eyeball the tail.
- `a-four-decimal-print-is-not-a-zero` — 4.21e-05 printed as `0.0000` and a whole mechanism was built on
  it being zero. `repr` the number a diagnosis turns on.
- `a-refutability-test-needs-the-refuting-channel-in-the-fixture` — a prior measured with no evidence that
  could overturn it looked catastrophic; the channel was missing from the chain, not from the design.
- `a-strength-is-a-nat-a-prior-weight-is-a-count` — a prior's strength comes from what it CLAIMS, never
  from what the lattice can hold. The representation's limit is a cap, and a cap is not a choice.
- `a-constant-in-exact-arithmetic-is-not-constant-in-float64` — a term that "cancels" carried 2.2e-16 of
  asymmetry and tipped a median knife-edge by a full grid step. Return the vanishing value.

**Structure, indexes and plumbing**

- `one-reference-hides-refid-bugs` — Single-reference synthetic indexes hide reference-id-space mismatches.
- `annotated-is-not-genomic` — "HAS AN ANNOTATION" IS NOT "IS GENOMIC" — a classification must be an INPUT, not a proxy.
- `an-sj-is-not-a-gap` — A splice junction CANNOT be detected from a gap between deposited slices.
- `deposit-at-the-sj` — A splice deposit belongs at the SJ's coordinate, not the region's boundary.
- `splicing-makes-the-graph-cyclic` — Alternative splicing makes a region↔sj graph CYCLIC, and cycles are the common case — a
- `nrna-does-not-mean-synthetic` — ~is_synthetic & ~is_nrna as the "real transcript" filter deleted 52,104 real termini.
- `credit-exactly-one-sj` — A fragment crossing K sj must credit exactly ONE (the leftmost annotated).
- `strand-completes-the-sj-key` — (src, kind, dst) is not a total order for sj boundaries — two strand-coincident sj diffe
- `a-hash-that-misses-its-artifact` — Cache keys that do not cover the artifact they cache.
  ⚠ 2026-08-08: the payload digest hashed field NAMES, so collapsing a bank from `[n,2]` to `[n]` did not move it — and the cache-load path validates no shapes. It now hashes `name:columns`.
  ⚠ 2026-08-08, second form: an EXCLUSION from a key is a claim about what is stored, and the study cache stored two arrays the excluded file produced.
- `integer-channels-reproduce` — Integer channels are bit-identical across worker counts; float channels are not. ⚠ Its ~2.6 % price tag is a float32 number and does NOT carry to float64, which is 3.4e9x finer — and measured MORE accurate than the fixed point it replaced.
- `worktrees-run-the-wrong-code` — Worktrees silently run the wrong code — an editable install's meta-path finder beats
- `checkout-deletes-uncommitted-work` — git checkout -- <file> does not undo a perturbation when the work is uncommitted — it deletes th
- `two-masks-one-name` — TWO DIFFERENT MASKS SHARED THE WORD struct_lock, AND BOTH WERE RIGHT.
- `two-docstrings-one-quantity` — The prose next to the code said "the AVERAGE" and the code followed the prose, while a sibling
- `an-object-class-does-not-see-a-terminus` — THE OBJECT CLASSES CLASSIFY A BOUNDARY BY ITS FLANKS' EXON-NESS, AND A
  TSS/TES LYING INSIDE A GENE LOOKS EXACTLY LIKE A SPLICE SITE. A hop type is keyed on the boundary's structural
  flags, or the whole composition error into exons lands on one un-named sub-class.

**Domain facts that read like defects**

- `specificity-and-sense-are-complements` — Strand specificity is TWO different quantities and they are complements.
- `strand-measures-the-tilt` — Strand measures the TILT, not the gDNA fraction.
- `a-linear-likelihood-emits-a-sign` — A LIKELIHOOD THAT IS LINEAR IN THE PARAMETER HAS NO MODE, ONLY A DIRECTION — and on a bounded grid that is a saturated vote.
- `amplitude-fades-influence-does-not` — A TERM THAT NORMALISES AWAY CAN STILL DECIDE THE ANSWER, BECAUSE AN ARGMAX IS SCALE-FREE.
- `a-pooled-conversion-applied-per-component` — A RATIO MEASURED ON THE POOLED POPULATION AND APPLIED TO EACH COMPONENT IS POPULATION-BLIND.
- `an-identity-with-a-qualifier` — AN IDENTITY THAT HOLDS "OVER X" IS A MEASUREMENT WAITING TO BE MADE — PRICE THE COMPLEMENT.
- `equal-lengths-carry-no-composition` — At equal component mean lengths the length channel carries EXACTLY ZERO information about
- `capture-is-1000x-on-exons` — Hybrid capture is ~1000× on exons and only gDNA reads it cleanly.
- `capture-selects-for-length` — Capture SELECTS FOR LENGTH, and the post-capture distributions are the baseline.
- `on-target-by-start-is-geometry` — "On-target" defined by the START's territory is geometry, not capture efficiency.
- `real-data-is-a-test-input` — Real data is a TEST input, never a DESIGN input.
- `eff-lengths-do-not-cancel-at-an-end` — "Effective lengths cancel, so a region's marginal is just its length" is FALSE near any transcript
- `configured-lengths-are-not-realised` — EQUAL CONFIGURED FRAGMENT LENGTHS DO NOT GIVE EQUAL REALISED ONES.
- `mature-rna-never-crosses-a-boundary` — Mature RNA never crosses an exon↔intron boundary (0 of 1,146 boundaries, 7/7 conditions).
- `a-boundary-with-rna-is-not-an-sj` — But a boundary with RNA crossing it need not be a sj.
- `the-panel-enriches-nascent-by-its-own-probes` — A SIMULATOR THAT MODELS A POPULATION IN ITS OWN PRIVATE
  SPACE GIVES IT ITS OWN PRIVATE PHYSICS. Nascent RNA had a per-transcript pre-mRNA space, so a probe on
  another transcript's exon enriched the gDNA there and not the nascent. FOUND and REPAIRED 2026-08-19.

**Process**

- `no-magic-numbers` — No magic numbers.
- `one-thing-varied` — One thing varied per experiment, with the falsification test written first and verified failing,
- `converge-and-delete` — Converge and delete.
- `the-source-does-not-cite-docs` — The source does not reference the docs.
- `running-an-arm-is-a-fresh-process` — SIX OPERATIONAL TRAPS FROM RUNNING PANEL ARMS, each of which cost a
- `shard-an-arm-sweep-by-condition` — Shard a panel sweep by CONDITION, never by ARM over one condition: the instruments WRITE their caches and three writers on one directory truncated a payload.
- `a-mean-hits-the-mass-weighted-centre-by-luck` — A scalar that "wins" on a mass-weighted metric may only be sitting on the mass-weighted centre of a bimodal population; compare it against its own class mean.
- `a-clamp-at-the-closed-end-escapes-the-window` — A prior clamped at the closed end of its support can put almost all its mass outside the solve window; derive the floor at the OBJECT, not over the pool.
- `the-deconvolution-is-as-good-as-the-density-it-is-handed` — A residual estimator inherits the error of the density it subtracts; price that density against truth first.
- `no-enumeration-without-a-census` — Do not re-propose path or cell enumeration without a memory census.

---

## A. Validation and gates

**self-checking-validator. A validator that calls the builder's own helper validates nothing.** Deleting a negative-strand
coordinate swap left a 1,289-test suite green and the graph validator accepting; only re-deriving the
result by a *different algorithm* caught it. Worse than no check, because it reads as one.
*Corollaries:* a validator that only compares non-empty classes never sees a spurious flag — emit all
classes always; and prove every validator FIRES by deliberately corrupting its input.

**perturb-every-gate. A falsification test needs PERTURBATION, not just to be written first.** Writing the gate before
the fix is half the discipline; the other half is breaking the fixed code and watching each gate fire.
This has found holes in already-green gates repeatedly — once 7 of 9. Three recent examples, all blind
until perturbed: a capture A/B whose two arms did not share their random input (so it "passed" with the
feature removed); a de-tilt test whose pool had **one length bin**, which is invariant under any
divisor whatsoever; and a gate written in the same session as its own fix that did not fire on its own
named perturbation, because a redundant `var > 0` backstop was silently doing the guard's job — **two homes
for one predicate** (TRAPS: a-test-that-redefines), inside the falsification.

**a-field-driven-gate-is-atomic. A GATE THAT ENUMERATES THE SPECIFICATION'S FIELDS FUSES SPEC, IMPLEMENTATION AND SCHEMA INTO
ONE COMMIT — plan the change that way or the tree is red in between.** `test_accumulator_native_parity`
reads `dataclasses.fields(Tally)` and `test_accumulator_payload` does the same, deliberately: "add a
field to the specification and it joins this gate automatically". ⭐ That is the right design and it has
a consequence nobody had written down.

⚠ Measured 2026-08-08: adding **one** field to the specification's `Tally` — nothing else — turned **10
tests red** across the two files. The accumulator-prior plan had scheduled the work as four steps
(spec → C++ → digest → consumer) and warned that "steps 3 and 4 must land together"; the truth is that
**2, 3 and 4 are one commit**, and step 4 is not a step at all — `payload_schema_digest()` recurses over
the payload's field list, so the bump happens automatically, and a gate forbids keeping the field
accumulator-internal to dodge it.

⭐ **The cheap way to learn this is to try it**: add the field, run the suite, read the failure count,
revert. Ninety seconds, and it re-planned the commit before any C++ was written. ⛔ Do not discover it
half way through, with the extension rebuilt and the tree red.
⚠ The downstream cost is the real one: the digest bump invalidated **201 cached scans** across three
panels. Budget the rebuild as part of the commit, not as cleanup after it.

**waive-with-a-measurement. ⛔⛔ AN ASSERTION THE SHIPPED CODE VIOLATES IS WAIVED **WITH ITS MEASUREMENT**, NEVER WIDENED — AND
THE COUNT IS WORTH MORE THAN THE CRASH.** Five structural invariants were added to the solver's backbone,
each written from a defect that had already shipped. Two fired on HEAD immediately, and the instinct — loosen
the predicate until the tree is green — is exactly how TRAPS: a-gate-that-reconstructs' gate became vacuous. ⭐ **What was done instead:
each violated invariant is listed by name with the number beside it, a test asserts every waiver carries a
written reason, and every count is published with its ELIGIBLE set (TRAPS: could-the-arm-have-fired).** That turned five gates into six
new measurements on the FIRST run, the largest being that **61.1 % of live message packets claim the three
components together account for more fragments than the slot observed** — a fact consistent with a number
already in the source and which nothing could rank because nothing had made it checkable.
⚠ **The other half of the lesson is a tolerance, and it is derived rather than tuned.** A bare "the mode is
inside the grid" reported a *correct* answer as a defect at 63 slots, because the shipped tilt mode is a
convex mean of two messages that agree on `τ = ±1` and the division rounds **2 ULP** past `π/2`. Allowing one
**grid SPACING** — the coordinate's own resolution, no constant — separates that from TRAPS: off-grid-message-mode's real overshoot of
**57 spacings** by four orders of magnitude. ⛔ *The rule:* when an invariant fires on shipped code, **measure
the rate before touching the predicate**; then either the predicate was wrong about the coordinate's domain
(fix it, from the shipped builder) or the code is wrong (waive it, with the number).

**a-docstring-that-misdescribes-the-graph. ⛔⛔ A CLAIM ABOUT THE CODE, INSIDE THE CODE, THAT NOTHING
GATES — AND IT ROTS EXACTLY LIKE A STALE DOC CITATION, ONE LAYER DOWN.** `run_fill` said it was "shared by
`density_model`, `strand_deconv`, `priors`, and the `sweep` chain geometry" and had **one** importer;
`strand_likelihood` said "Used by the per-region strand module (`strand_deconv`)" and `strand_deconv` does not
import it at all. A developer who trusts either sentence goes looking for code that is not there. Measured
across the package: **14 module docstrings named a sibling with no import boundary in either direction, and 6
were genuinely stale.** ⭐ *The tell is free and nobody was reading it:* the import graph is in the AST, so
prose about it can be CHECKED — `scripts/design/module_census.py` does. ⚠ *And the instrument cannot finish
the job:* "fitted in X" is a DATA-FLOW claim, true of a value passed through a caller, and only a human can
tell it from an import claim. The count is a worklist, not a verdict.

**a-flat-pile-is-not-a-knot. ⭐⭐ 35 MODULES WITH NO CYCLES AND ONE IMPORTER EACH IS NOT TANGLED — IT IS
UNLAYERED, AND THE TWO NEED OPPOSITE FIXES.** The calibration package read as unmaintainable and the
instinct was to merge files. Measured from the AST it had **no import cycles at all** and **18 of 35 modules
had exactly one importer** — so nothing was entangled; there was simply no stated ORDER, and a flat pile of
peers is the one shape that cannot tell you where to add anything. ⭐ The layers were already in the boundaries:
naming them in `_layers.py` and gating them cost **zero behaviour change**, and the gate immediately found
two upward imports that were the same defect twice — `NodeDeconv`, the tool's central datum, defined in the
STRAND family and reached up for by three layers. ⛔ *The rule:* **before merging modules, ask whether the
problem is entanglement or missing order.** Merging a flat pile makes bigger files with the same problem.

**the-rename-that-corrupted-a-diagram. ⛔⛔ A MECHANICAL REWRITE OVER PROSE WILL HIT TOKENS THAT ARE NOT
PROSE, AND THE ONES IT HITS ARE THE ONES YOU DID NOT ENUMERATE.** Renaming 980 rule citations rewrote **20
slot ids inside ASCII chain diagrams** — `N0 E0 N1 E1 … E(k-2) N(k-1)` became a sentence — because `E1` is a
label in one namespace and a SLOT in another. The rename also stopped at `.py` and `.md`, so the C++ kept
its labels and a Python test asserting a deletion record ("DELETED by C2") failed against the header it
reads. ⭐ *The tell:* the gate written to forbid the labels is what found both, by flagging what remained.
⛔ *The rule:* **a mechanical rewrite needs the gate FIRST — then the residue it reports is the list of
things you did not think of**, and every one of them is a collision worth understanding rather than
exempting. Same family as `TRAPS: two-masks-one-name`, met while trying to fix it.

**a-gate-that-already-passed. A gate that already passes is not a falsification.** A simulator gate on "is on-target gDNA longer
than off-target?" passed *with the defect present*, because the engine's per-fragment **conditional** was
correct and only the **marginal** was being discarded. Before trusting a gate set, check which of them
would have failed before the fix. If none, the gate has not been written yet.

**right-conditional-wrong-marginal. When a per-observation conditional is correct, no conditional comparison can detect a missing
marginal.** The general form of TRAPS: a-gate-that-already-passed, and worth stating separately because it tells you *where* to look:
if `P(outcome | x)` is right for every `x` but the distribution of `x` is wrong, every conditional check
passes and only a check on the marginal fails.

**byte-identity-gate. A bit-identity gate has lied in both directions.** An arm with ZERO rows scored "32/32 IDENTICAL"
because the comparison looped over the new arm's rows; and a stored baseline went stale so unmodified
HEAD no longer reproduced it. **Re-record the baseline from the current tree in the same session**; if
HEAD-vs-baseline is not 100 %, the baseline is what is broken.

**the-deliverable-is-not-reproducible-by-default. ⛔⛔ THE SHIPPED PIPELINE DOES NOT REPRODUCE ITSELF, SO AN END-TO-END A/B ON
THE DEFAULT CONFIG MEASURES ITS EFFECT PLUS A SAMPLING DRAW.** `EMConfig.seed` defaults to **`None`** and
`assignment_mode` to **`"sample"`**, so the EM's final hard assignment is an *unseeded* categorical draw.
Two back-to-back runs of the identical pipeline on the identical BAM returned different transcript counts
— measured on the gate toy: **4 transcripts differing by up to 43 fragments**. Setting any fixed seed makes
it exact, and so does `assignment_mode="fractional"`.
⛔ **Found by the byte-identity gate above, on the arm that was supposed to be free.** The first
end-to-end `noop` arm — a wrapper that builds the oracle prior and then discards it — came back differing
on 80 % of transcripts, which reads as "the injection plumbing is broken" and is not. Every arm in
`quant_accuracy.py` now pins the seed, and `base_reseed` re-runs the baseline at `seed + 1` so the noise
floor is printed in the same units as the effect.
⚠ **This is a property of the tool, not of the harness, and it has not been ruled on.** A quantifier whose
output moves run to run is a defensible design (the draw is how a hard assignment is made from a posterior)
and it is also a surprise to a user diffing two runs. Whether the default should change is an owner call;
what is settled is that **no measurement of the deliverable is attributable without pinning it.**

**a-clip-hides-a-scale-error. A `min()` clip hid an exact factor of 2 for months.** Pooled boundary density read 1.994× truth
because the code SUMMED two faces' mass and divided by their AVERAGE. It survived 29 tests: all four
fixtures stored an un-halved length AND deposited half the mass (cancelling exactly), and `min(2S,S)=S`
under a uniform field, so the "contraction is 1 under a uniform field" bedrock invariant — written to
catch exactly this — passed with the bug present. *Lessons:* a clip can hide a scale error from a bedrock
test; repair fixtures, never relax assertions; **an exact algebraic 2 is never a modelling
approximation.**

**an-inverted-clip-is-a-constant. ⛔ A TWO-ENDPOINT CLIP HAS TWO SILENT DEGENERACIES AND BOTH TURN THE
ESTIMATOR INTO A CONSTANT — AND AT A ZERO CONTROL THAT CONSTANT IS THE RIGHT ANSWER.** **(a)**
`np.clip(m, lo, hi)` with `lo > hi` returns `hi`, with no warning and no error — verified directly.
**(b)** `lo == hi` returns that one value, which is not a defect in `np.clip` at all and is the harder of
the two to see. Measured 2026-08-17 while pricing an anchor-derived bound on an exon's `f_g`: at the `g00`
(zero-gDNA) condition **every anchor density is 0**, both endpoints collapse onto 0, and the arm asserts
`f_g = 0` unconditionally — the answer is right, and the zero control is testing nothing. ⭐ *The rule:*
**a zero control passed by a degenerate estimator may not be quoted as evidence that the estimator
DISCRIMINATES.** ⛔ *And the guard has to be the right one:* `assert lo <= hi` catches form (a) and
**PASSES form (b), which is the one that actually happened** — so the check that bites is to require the
SAME estimator to **MOVE at a contaminated condition** before either number is read, and to report the
interval's WIDTH at the control beside its verdict, because a zero-width interval is the tell. ⚠ Siblings,
each a different way the same thing hides: `TRAPS: a-clip-hides-a-scale-error` (a clip hiding a scale),
`TRAPS: an-all-zero-factor-is-inert`, `TRAPS: a-single-level-panel-cannot-see-a-constant` (a panel that
cannot tell a good estimator from a constant).

**an-ablation-that-never-ran. AN ABLATION THAT NEVER RAN READS AS "NO EFFECT", AND TWO IMPORT HABITS MAKE THAT EASY.** An arm
deleting a whole channel scored **byte-identical on every class of every condition** — a clean, publishable
"this channel is inert". It was inert only because the monkeypatch missed the module that binds the name
(`from .sweep import X` makes a separate module global), and the first repair missed it too because
`from pkg import mod` returns a re-exported *function* of that name which shadows the submodule.
⭐ *The tell:* the ablation was a strict SUPERSET of another and scored a smaller effect. ⭐ *The rule:*
every ablation increments a counter and the harness **raises** if it did not fire.

⚠ **AND IT FIRED TWICE MORE ON 2026-08-07, both times because a NAME MOVED rather than because an import
habit hid it.** (i) `import rigel.calibration.calibrate as CAL` binds the re-exported **function**, not the
module, so `CAL.solve_chain = spy` patched an attribute on a function object and the spy read as "never
called" — the fix is `sys.modules["rigel.calibration.calibrate"]`. (ii) When the sweep was renamed, three
instruments still patched `node_sweep`; patching a name the caller no longer calls raises **nothing**, sets
no `region_arrays`, and every override arm downstream would have silently overridden nothing. ⭐ *The rule
that generalises:* **a patch site must assert the name it is patching still EXISTS**, because a rename turns
a monkeypatch into a no-op with no error at all — and a no-op ablation is indistinguishable from an inert
mechanism.

⛔⛔ **AND ON 2026-08-11 THE COUNTER ITSELF WAS PROVED INSUFFICIENT — BY A CONFIG DEFAULT, NOT BY AN IMPORT.**
`ladder_arm_ab.py` carries 26 arms, every one of them counted and guarded by the rule above. Measured under
the shipped `message_propagation = False`: **22 could not move a number.** Six raised (four on a module
deleted at `0d9d422b`). The other sixteen **fired on tens of thousands of slots and scored byte-identical to
`base`** — eight muted a relay that was already silent, and eight moved a PRECISION nothing consumes, because
`_fuse(own, p, msg)` returns `own` for every `p` when no message arrives. The counter said *fired* for all
sixteen. ⭐ *The rule:* **a fire counter answers "was my code called"; the question is "can this arm move a
NUMBER", and only a diff against a base run under the same config answers it.** The harness now runs every arm
under both policies and requires INERT under one and MOVED under the other.
⚠ *And the deeper half:* **the config is part of the arm's definition.** `main` built a bare
`CalibrationConfig()`, so the day a shipped default flipped, twenty-two arms became no-ops with no error, no
warning and no diff — the files still parsed, still carried every row, and still scored. The setting is now
stamped into every row, which is `TRAPS: a-hash-that-misses-its-artifact` in its second form: **the artifact
must carry its own key.**

**an-inert-arm-reads-as-a-refutation. ⛔⛔ AN ABLATION THAT NEVER RAN READS AS "NO EFFECT"; A *PROTOTYPE*
ARM THAT NEVER RAN READS AS "YOUR IDEA DOES NOT WORK" — AND THAT ONE KILLS A DESIGN DIRECTION.** Same
mechanism as the rule above, strictly worse consequence, and it cost real measurements on 2026-08-17: three
arms testing a new exon reference were driven by `import rigel.calibration.calibrate as CAL` and
`CAL.solve_chain = ...`. `rigel/calibration/__init__.py` does `from .calibrate import calibrate`, so that
name binds the **FUNCTION**, not the module, and the patch set an attribute nothing reads. A **314-second**
run reported *all three arms byte-identical* — which reads as a clean, publishable NEGATIVE rather than as
an error, and it was believed.
⛔ **This is `TRAPS: an-ablation-that-never-ran` in a new costume, and its (i) records this EXACT import
from 2026-08-07** — so the mechanism and its fix (`sys.modules["rigel.calibration.calibrate"]`) live there
and are deliberately not repeated here. ⚠ Note what that means: the rule was already written, already
named, already cited, and the mistake was made anyway. A rule a reader has to remember is not a guard.
⭐ *What is new is the ORDER, and it is the rule:* **a monkeypatch arm ASSERTS ITS TARGET EXISTS and COUNTS
ITS OWN FIRINGS before any number is read** — not after the run, not when a result looks odd; a zero-firing
arm raises instead of scoring. ⚠ *Why the order matters more here than for an ablation:* an inert ABLATION
at least looks suspicious, because "no effect" invites a check. An inert PROTOTYPE looks like the negative
result you were braced for, and nobody audits a confirmation.

**a-green-suite-hid-five-dead-instruments. ⛔⛔ AN INSTRUMENT DIES WHEN THE TREE MOVES UNDER IT, AND A
GREEN SUITE SAYS NOTHING — BECAUSE NOTHING CHECKED THAT A SCRIPT STILL IMPORTS.** `tests/test_scripts_index.py`
checked only that a script was INDEXED and had a DOCSTRING, both true of a file that raises on line 1.
⭐ *The name records the FIRST mechanism, a `src/` DELETION:* `94d283c0` deleted the fixed-point layer
(`INV_LENGTH_SCALE`, `inv_length_quantum`) and killed three instruments, `0d9d422b` deleted
`enrichment_frame` and killed one, and one more died on `_component_node_arrays` — **five**, all found at
once. ⭐ *The rule as first written:* an import gate, and **run the suite after a `src/` deletion**.

⛔⛔ **AND IT RECURRED ON 2026-08-17 THROUGH A MECHANISM THAT RULE DOES NOT COVER: A CONFIG DEFAULT FLIP,
WITH NOTHING DELETED AND NOTHING RENAMED.** The diagnostic key `_uni` is written **only** by
`messages/relay.py` — i.e. only under `RelayPolicy` — and `CalibrationConfig.message_propagation` ships
`False`, which installs `SilentPolicy`. **Six more instruments** were dead: five raised
`KeyError: '_uni'` and `reframe_walk.py` raised `KeyError: 'rho_lo'`, another RelayPolicy-only capture.
⛔ **The import gate cannot see any of it** — every one of those files imports perfectly; they die on the
first read of a key the shipped configuration never writes. ⛔⛔ **And the suite stays green for a reason
worth stating on its own: the TEST readers install the policy THEMSELVES.**
`tests/calibration/test_gdna_scale_rule.py` binds `functools.partial(solve_chain, policy=RelayPolicy())`
and says why (its assertions would be vacuous otherwise) — so the tests exercise a configuration the
instruments do not, and no gate spans the gap. ⚠ A seventh, `psi_channel_ablation.py`, was found dead in
a way **no** config revives.

⭐ *The rule, extended:* **an instrument is alive only under the configuration it is RUN in, and a shipped
default is part of that configuration.** Run the instruments — not the suite — after a `src/` deletion,
after a rename, **and after flipping any shipped default**, because a default decides which diagnostic
keys exist at all. ⚠ Same family as `TRAPS: an-ablation-that-never-ran`'s third form (*"the config is part
of the arm's definition"*), one layer out: there a live arm could not move a number, here a whole
instrument cannot start.

**compatibility-is-geometry-not-composition. ⛔⛔ A PREDICATE BUILT FROM WHAT A FRAGMENT IS COMPATIBLE
WITH CANNOT SEPARATE RNA FROM gDNA — BECAUSE COMPATIBILITY IS GEOMETRY, AND gDNA IS COMPATIBLE WITH
EVERYTHING UNSPLICED.** Designed and refuted the same day, 2026-08-11/12, and the refutation is the
reusable part.

The proposal was to give a synthetic nascent entity warm-start mass only if it had **RNA-unambiguous
support** — at least one fragment whose candidate set, with the gDNA component removed, was exactly that
entity. The separating argument looked airtight: an intronic fragment is unspliced, no annotated
transcript admits it (the annotation splices that interval out), and only the nascent shadow span
contains it — so `C_R(u) = {nascent}` and the entity lives, while a zombie whose every fragment is also
admitted by a mature isoform has `|C_R| >= 2` everywhere and dies.

⭐ **Nothing in that argument requires the fragment to be RNA.** A genomic-DNA fragment landing in an
intron satisfies every clause. So the gate fires on intronic gDNA, and **the MORE gDNA a library
carries, the MORE nascent entities are revived** — the exact inverse of the intent. The predicate is
sound only at `f_gdna = 0`, which is the one library nobody runs.

⛔ *The rule:* **any RNA-vs-gDNA question resolved by SET MEMBERSHIP will be answered by geometry.** It
has to be resolved by LIKELIHOOD, where the two populations actually differ (strand, length, density,
splicing). ⭐ The replacement — give the entity no prior mass and let it earn its place — is immune to
the same counterexample precisely because it arbitrates the contested intronic fragment by likelihood
against a gDNA component that DOES hold prior mass, instead of by asking which sets contain it.

⚠ **And the measurement that motivated the refuted design was VACUOUS, which is its own warning.** It
reported synthetic nascent entities holding 1,037,811 fragments at **0.000 %** unambiguous support
against 15.551 % for annotated mRNA — devastating, and empty: `em_solver.cpp` appends a gDNA candidate
to *every* unspliced fragment, so an unspliced fragment can never be `unambig`, so a nascent entity's
`unambig_totals` is **structurally 0 on every input**. The statistic was never in a position to say
anything about the population it was being used to indict (`TRAPS: could-the-arm-have-fired`). It was
persuasive, and being persuasive is what made it dangerous.

**a-zero-count-is-a-measurement. ⛔⛔⛔ A ZERO COUNT IS A MEASUREMENT OF A DENSITY, NOT AN ABSENCE OF DATA — AND KEYING PRECISION
ON THE COUNT MAKES THE STRONGEST STATEMENT IN THE LIBRARY THE QUIETEST.** At `g00` (zero gDNA by
construction) pass-0 attributes **34–38 %** of unstranded mass to gDNA, Σ|err| **1.56 M / 1.95 M**
fragments, and **100 %** of it is `relay_only` with **0 %** `own_evidence` and **0 %** `struct_lock` on
both axes. The reason is a single census: all **1,298** intergenic regions hold **exactly zero** counts
over **50.7 Mb** of gDNA opportunity, and **all 1,298 emit nothing** — against 915–959 emitting at
`g25`/`g50`, so the anchor is healthy elsewhere and completely silent here. ⭐ Those regions are
structurally pure gDNA and `struct_lock`ed, i.e. composition-**CERTAIN** — and `own_precision`'s
`n > 0` still zeroes them, so they are *certain and silent at the same time*. ⛔ **The general defect is
that `p = n/(n·Var+1)` reads the COUNT as the evidence when the claim is a DENSITY: zero over 50.7 Mb
and zero over 200 bp are the same number to this code and opposite statements about the world.** ⭐ The
derived repair needs no new constant and reuses the reference ψ already carries: a Poisson rate with
`a` events over exposure `E` under the Jeffreys prior ψ is built on (`_JEFFREYS_REF`) has posterior
`Gamma(a+½, E)` — proper, finite precision at `a = 0`, mean `0.5/E`. At `ΣE = 50.7 Mb` that is
`9.9e-9 ± 1.4e-8`, which would take the worst object from `f_g = 0.510` to ~`7e-9` against a truth of
exactly 0. ⚠ The population prior papers over this downstream (that object's shipped answer is 0.054),
but pass-0 **is** the substrate the prior is fitted on, and a prior fitted on "everything is ½" cannot
be rescued per-object on a heterogeneous sample.

**a-ratio-cannot-carry-zero. ⛔⛔⛔ THE RELAY'S CURRENCY IS A DENSITY *RATIO*, AND ZERO IS NOT IN THE MULTIPLICATIVE GROUP —
SO "THERE IS NONE HERE" IS UNREPRESENTABLE BY CONSTRUCTION, NOT BY AN OVERSIGHT.** TRAPS: a-zero-count-is-a-measurement's repair was
built, gated (18 gates, 4 firing perturbations) and **measured at +0.0 % on the 36-condition panel, with
the four `g00` rows BYTE-IDENTICAL** — the conditions it was designed for did not move at all. The
two-step probe says exactly why: after the fix `NodeInit.prec_g > 0` at **1,298** intergenic slots
(value exactly `1/trigamma(½) = 0.2026`), and the relay publishes `prec_g > 0` at **0** of them. ⭐ The
reframe transports `r = ρ_tot(dst)/ρ_tot(src)` and every guard on that path is `rho > _EPS`
(the reframe's guards in `messages/relay.py`): a ratio with a zero denominator is undefined, so a source whose true
density is **0** cannot participate however precisely it knows it. ⛔⛔ **That is why the single most
informative statement available at a zero-gDNA library — 50.7 Mb of opportunity with no fragments — is
structurally unsayable.** ⭐ The lesson generalises past this solver: *before fixing a source that will
not speak, check whether the channel can carry the value it wants to send.* A multiplicative transport
cannot carry zero, and no amount of precision at the source changes that. ⚠ And the fix is not wasted —
it is a precondition: the source must be able to form the claim before the channel can be taught to
carry it. But it must be priced as what it measured, which is nothing.

**the-divergence-was-a-barrier. ⛔⛔⛔ THE DIVERGENCE WAS DOING A SECOND JOB, AND REMOVING IT TURNED A RELAY BARRIER INTO A
CONDUIT.** TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero replaced ``1/n`` with ``trigamma(n+½)`` in **two** places and won 39 % panel-wide — and
one coherent stratum, stranded × capture-ON at ``g10``+, got **20–34 % worse**. The whole regression is
owned by ONE of the two: ``composition_logvar``'s counting term, i.e. by ``σ²_transfer =
logvar_tot[dst] + logvar_tot[src]`` no longer being ``∞``. Reverting only that term reproduces the published
"before" column **to the fragment** (295,453 and 304,815), and a from-git pre-fix tree agrees to 0.17 %.
⭐⭐ **The mechanism is not the zero-count SOURCE, it is every hop that TOUCHES a zero-mass slot.**
``1/(1/p + ∞) = 0`` annihilated messages in both directions, so ``_fuse`` fell back to the destination's own
belief and **the chain was cut into segments at every empty slot**. Under capture the off-probe gaps between
the gaps between probe-covered stretches ARE the zero-mass slots, so a gDNA level now travels between capture strata that were
isolated — and it crosses UNSCALED, because ``framed`` needs ``ρ_tot > 0`` at both ends and forces ``r = 1``.
⛔ Muting every zero-mass emitter recovers **1.2–7.7 %**; restoring the barrier recovers **100 %**. So it is
PROPAGATION, not origination, and the previous session's hypothesis — an empty intergenic anchor asserting
``ρ_g ≈ 0`` — is refuted as the cause even though its premise is true and measured (the true gDNA density at
those anchors' neighbours is **346×** what a non-empty anchor reports). ⭐ **The lesson: before crediting a
divergence's removal, ask what the divergence was suppressing.** An ``∞`` in a damping term is a
STRUCTURAL gate wearing a variance's clothes, and it had been the only thing pricing a premise nothing else
prices — here "a gDNA level does not change across a BOUNDARY", which capture falsifies by ~350×.

**deadband-from-the-wrong-sample. ⛔⛔ A NOISE DEADBAND WHOSE CUSHION IS SUPPLIED BY AN UNRELATED SAMPLE SIZE FAILS EXACTLY WHERE
THAT SAMPLE GETS BIG — AND IT FAILS SILENTLY, INTO THE HONESTY COLUMNS.** `strand_evidence` gates the
strand channel on `disc = 4·max(0, (κ̂−½)² − σ²_d)` with `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`.
On unstranded data the TRUE κ is ½, so `(κ̂−½)²` is pure fitting noise — a χ²₁ with mean `Var(κ̂) ≈
¼/N_rna`, which a χ²₁ exceeds ~32 % of the time. ⭐ **The `¼/N_rna` half of the deadband is exactly the
right bias correction; the `¼/N_gdna` half is a different job** (gating a gDNA-free library) doing the
first job's arithmetic. So the cushion **shrinks as gDNA grows**, and the gate opens on **exactly 3 of 18**
unstranded ladder conditions — `g75`/`g90 ss0.50 capture_off` and `g98 ss0.50 capture_on` — measured
`(κ̂−½)²` = 4.74e-07 against `σ²_d` = 3.42e-07 at g75, a margin of 1.39×. ⛔ **And `I_strand = N_eff·disc·[…]`
then MULTIPLIES noise by the depth**, so the phantom information grows linearly in coverage. ⭐⭐ The damage
is not accuracy — forcing κ = ½ moves Σ|err| by −0.6 %/−5.3 %/**+4.6 %** across the three big classes, a
wash, reproducing the −0.2 % that `ROADMAP.md` §4's κ-residue row records. **The damage is that `solv%` goes 77–90 % → 0.0 %**: the panel
declared 92,154 objects confidently wrong and `calib` = 3.67 at a condition with *no solvable objects at
all*. ⛔⛔ **So the column used to pick the worst condition was inflated by the defect being looked for** —
TRAPS: honesty-metrics-reward-ignorance's lesson with the arrow reversed. ⭐ The repair is to propagate rather than gate:
`I = (½−κ̂)² / [p(1−p)/N_eff + (1−f_g)²·Var(κ̂)]`, which is smooth, needs no `max(0,·)`, and whose
`Var(κ̂)` denominator cancels the depth exactly. ⚠ Not built: `ROADMAP.md` §4's κ-residue row prices
the gate's ACCURACY cost at −0.2 %, so this repairs the honesty columns rather than the deliverable.

**honesty-metrics-reward-ignorance. ⛔⛔ EVERY HONESTY METRIC IMPROVES AS THE SOLVER STOPS KNOWING ANYTHING — SO NONE IS READABLE
WITHOUT A FIXED-DENOMINATOR ACCURACY NUMBER BESIDE IT.** A destruction control (force a 97-σ real strand
signal to ½) collapsed accuracy **0.0231 → 0.2064, +792 %, 16/16 worse** — while confidently-wrong Σ|err|
fell **87 %**, the declared-precision ratio went 4.53 → 0.67 and `weak%` went to 0, all 16/16 "better".
Three of those four are the Stage-B dashboard's headline columns. ⛔ **The rule: an A/B that moves `solv%`
has changed its own denominator, so quote `mwae` over ALL objects and the raw Σ|err| alongside — those two
cannot be gamed by knowing less.** ⚠ It also inverts single-condition readings: one row showed −64 % and
the same arm over 16 conditions is −0.2 %.

**predicate-contradicts-its-docstring. ⛔⛔ A PREDICATE CAN CONTRADICT ITS OWN DOCSTRING FOR MONTHS IF SOMETHING ELSE IS MASKING IT —
AND THE WRONG VERSION CAN BE LOAD-BEARING.** ``strand_evidence``'s ``struct_lock`` is documented as "scoped
to true intergenic REGIONs" — i.e. ``g1_locked ∧ REGION``, and `node_geometry.g1_locked` exists as the
designated ONE HOME for exactly that predicate (TRAPS: a-test-that-redefines). The code is handed ``locked = ~solvable`` with
``solvable = (free_pos|free_neg) & (n > 0)``, so it is true at **every zero-count REGION**: measured **19,709**
ladder slots against **1,312** that are actually structurally pure-gDNA objects, so **18,397** empty exons and introns declare their
composition CERTAIN. ⚠ It was INERT until 2026-08-06 because ``own_precision``'s ``n > 0`` gate silenced
every zero-count slot; removing that gate un-masked it. ⛔⛔ **And scoping it correctly is
PANEL-NEGATIVE**: the target row moved −1.2 %, ``g98`` went +0.4 % (worse), and the zero-gDNA control went
**+3,207 %** — the mis-scoped mask is what carries the zero-gDNA win. ⚠ **RE-PRICED 2026-08-18 on the
32-condition ladder with the per-strand licence + SPLICE IN fix live, and the sign on the control has
FLIPPED but the verdict has not**: on top of the SPLICE IN fix, `g1_locked ∧ REGION` takes the eight `g00` rows
**324 k → 232 k (0.71×)** — because the empty exons' own "RNA = 0 @ 0.2026" is one of the two feeders of the
pin's "all gDNA" (TRAPS: zero-the-precision-with-the-value) — while REGRESSING the four `nrna_none` zero
controls **1.8–15×** and the in-scope contaminated strata **+2.5 / +2.1 / +0.5 %** (deferred +17 %) — every
ratio here is against the SPLICE IN-fix tree. The
loss mechanism is named: at AMBIG `exon|exon` boundaries in an RNA-only library ψ was held at `f_g ≈ 0` by
the empty exons' "gDNA = 0 @ 0.2026" — a zero count over ~0.05 placements, right at `g00` for any reason —
and without it an RNA+ level claim with no RNA− claim drifts ψ to 0.38. Still half of a pair; the pair
is now sharper (see the dev-sandbox handoff). ⭐ *Lessons:* a gate that hands the
predicate in as an ARGUMENT cannot see the caller compute it wrong (`test_strand_evidence_struct_lock_regions_only`
did exactly that, and passed throughout) — gate it through the PRODUCTION path; and when a docstring and its
code disagree, find out which one the panel is relying on before "fixing" either.

**a-test-that-redefines. A TEST THAT RE-DERIVES A DEFINITION CANNOT DETECT DRIFT IN IT.** A gate existed precisely to keep
two instruments' shared class definition from diverging — and it recomputed that definition inline in
the test. Changing one instrument fired nothing. It is TRAPS: self-checking-validator's shape once more: the check and the thing
checked have to come from *different* places, and for a shared definition that means **one home**, with
every consumer importing it. Here the home is production code (`node_geometry.g1_locked`), because the
predicate is a production concept and `scripts/` is deliberately not importable.

**a-gates-power-is-its-invariant-set. ⭐⭐ A COMPARISON GATE'S POWER IS THE SIZE OF ITS INVARIANT SET, SO A
CHANGE THAT SHRINKS THAT SET MUST NOT RIDE ALONG WITH THE CHANGE IT IS MEANT TO GATE.** `rescan_panels.py`
gates a schema change by rebuilding every cache and demanding that **every bank the change did not touch
comes back identical** — its whole power is *"31 shared banks, every one identical"*. Bundling a ~1,000-site
mechanical rename into the same re-scan was proposed and refused (2026-08-13): renaming the banks leaves
almost no shared banks, so the rebuild would be gated only by "the same symmetric difference appears on
every condition" — which cannot see whether a renamed bank's **values** survived. The gate would still read
green, and would have proved nothing. ⭐ *The tell:* before bundling, ask **how many objects the gate will
still be comparing afterwards**; if the answer is "almost none", the gate has been disarmed rather than
satisfied. ⚠ **A second consideration decided it as much as the first: while every cache is refused, no
instrument runs** — so the rename would have been made BLIND and then landed against a gate that could not
fail, which is the worst available combination. Same family as `TRAPS: one-thing-varied`, one level up: not
"two changes confound each other's measurement" but "the second change deletes the first one's measuring
instrument".

**a-gate-that-reconstructs. ⭐⭐ A GATE THAT RECONSTRUCTS A VALUE IS VACUOUS WHEREVER THE CODE DECIDES NOT TO USE IT — GATE A
SYMMETRY THE CODE CANNOT FAKE INSTEAD.** A relay gate recomputed the delivered gDNA level and asked which
candidate the relay was closer to. It passed the honest run **and the deliberately mis-paired one**, because
on an unstranded fixture the composition licence is withheld, the level crosses unscaled, and the ratio
being reconstructed never enters the observable. Eleven green gates; the one perturbation they existed for
fired nothing. ⭐ The repair was to stop reconstructing: on a PALINDROMIC chain the relay's two passes must
be exact mirror images, which needs no assumption about which channels are live. ⛔ **When a gate needs a
value to flow through a conditional path, it is a gate on that conditional. Prefer an invariance.**

**off-grid-message-mode. ⛔⛔⛔ A MESSAGE MODE OUTSIDE ITS GRID'S DOMAIN IS NOT A WEAK CLAIM, IT IS A PIN AT THE BOUNDARY —
AND A CHANNEL ABSENT FROM THE CAPTURE IS A CHANNEL NO DISSECTION CAN RANK.** The `eta` prototype delivered
a raw log-odds `log(u₊/u₋)` into ψ's `theta_imp` slot. ψ's tilt grid is the **ANGLE** `arcsin(τ)` and
spans exactly `[−π/2, +π/2]`; the delivered modes were **±4.6**, 2.9× outside the whole domain. A Gaussian
`−½p(θ−m)²` with `m` off-grid is **monotone across every grid point**, so the tilt pinned at `τ = ±1`, the
AMBIG Schur protection that keeps the strand term out of `f_g` was destroyed, and the strand likelihood
explained the residue by calling the mass gDNA. **74 % of the `g00` error, one unit error.** ⭐⭐ *Three
lessons, and the second is the general one.* (i) **An out-of-range mode is the most confident statement a
channel can make, not the least** — the penalty has no interior minimum, so precision buys a corner rather
than a location. (ii) ⛔ **The tilt was the ONE message stream the prototype did not publish into
`_capture`, and five hypotheses died in the space it left.** The tell was free and unread the whole time:
the delivered `f_g` was the identical value `0.9898` at slot after slot — a **grid point**, i.e. a
posterior pinned into one cell. ⭐ *The rule:* **every channel a solver delivers must appear in its debug
capture, or no instrument can rank it and its absence reads as innocence.** TRAPS: could-the-arm-have-fired one level up — not "the
arm could not have fired" but "the arm could not be SEEN". (iii) ⚠ **Correcting the coordinate was NOT
the fix**: TRAPS: strand-measures-the-tilt gives `τ_obs = (2κ−1)(1−f_g)·τ`, so a raw COUNT tilt is the tilt of the *mass* and inverts
sign on an antisense protocol (`κ = 0.0101` ⇒ `2κ−1 = −0.98`). Build a message in the destination's own
coordinates — here the self-solve's RNA densities, which are ψ's `f_pos`/`f_neg` up to a factor that
cancels — and neither κ nor a belief is needed.

**a-comment-quoted-as-a-finding. ⛔⛔ A WORD BORROWED FROM A COMMENT IS NOT A MEASUREMENT, AND ONE PROPAGATED INTO THE DOCS INVENTED
A REGRESSION THAT NEVER HAPPENED.** the predecessor solver called its scan "a SEQUENTIAL Gauss-Seidel scan" — its
own shorthand for *in-place, therefore un-vectorisable*, written to justify the 15.7× do-not-merge twin. A
session read that phrase, wrote "HEAD is Gauss-Seidel, the rebuild is not" into `ROADMAP.md` and a
session handoff as if it were a structural finding, and a reviewer reading only the docs then
correctly derived a build plan from it — for a defect that does not exist. **Both sweeps run one forward
and one backward pass, both accumulate in place, and both combine at slot `i` from each neighbour's state
and never `i`'s own.** On a chain an in-place forward accumulation *is* the forward half of
forward-backward. ⭐ *The real difference the word was hiding is a good question:* HEAD relays **ten**
arrays (three component densities, three mode precisions, three MEASUREMENT precisions, `tau`); the
rebuild relays **four**. So the question is **which streams belong in the relayed state**, which is
answerable and already has one measurement against it (a chain-fused LEVEL is
dominated by the intergenic anchors and hands every exon the off-probe floor, 346×). ⛔ *The rule:*
**before repeating a source comment as a finding, check what the comment was written to justify** — and
never let a term of art cross from a code comment into a design doc without re-deriving it from the code.
⭐ **Settled by construction 2026-08-07:** the backbone's `_scan` is now literally `for i in seq: step(s, i)`
over one direction, called twice, so "one forward pass and one backward pass" is readable in six lines
instead of inferred from a 1,635-line function — and the word is gone from the tree.
Same family as TRAPS: two-docstrings-one-quantity, one layer up: there, two docstrings disagreed about one quantity; here, a docstring
and a document agreed on a word that meant different things in each.

**the-intermediate-is-not-the-deliverable. ⛔⛔⛔ AN INTERMEDIATE METRIC IS NOT THE DELIVERABLE UNTIL SOMEBODY MEASURES THE COUPLING — AND HERE
IT IS 10×.** Stage B was ranked, campaigned and celebrated on pass-0 `Σ|err|` for months. Measured
2026-08-06 on the same runs, in the same session, by reading a column that had been written all along
(`abs_err_all_final`): a **−37.2 %** pass-0 win is **−3.9 %** on the shipped solve, because the fitted prior
was already compensating for most of what pass-0 got wrong. ⛔⛔ **And the same change's REGRESSION grew
through the refit, +14.9 % → +30.0 %** — a pass-0 biased in the prior's own direction trains the prior to
repeat the bias, so the refit compounds instead of correcting. ⭐ *The rule:* when an intermediate stage
feeds a fitted stage, an improvement upstream can be absorbed and a bias upstream can be amplified —
**report both columns on every arm**, and never rank on the upstream one alone. ⚠ The tell was free the
whole time: the instrument already emitted the downstream number.

**could-the-arm-have-fired. ⛔⛔ "THE ARM CHANGED NOTHING" IS NOT A CONTROL UNTIL YOU CHECK THE ARM COULD HAVE CHANGED
SOMETHING.** A zero-gDNA condition was reported as a free falsification control because the new and old
reframe ratios agreed on every hop. They did, and the arm tested nothing: the composition licence was
withheld everywhere (no gDNA ⇒ no gDNA precision ⇒ no source can lend), and **6 of 8 hops per direction were
NO-FRAME** — the intron REGION had zero counts, so the reframe was skipped with `r` forced to 1 by the
divide-by-zero guard. The only framed hops were the ones the change does not touch.
⭐ **For every "no change" arm, count the opportunities the change had to fire and print that count beside
the result.** Same family as TRAPS: a-gate-that-reconstructs and TRAPS: can-the-benchmark-resolve-it: *the measurement was never in a position to say anything.*

⚠ **2026-08-10, and this form is a TEST FIXTURE rather than an arm.** The `sj_mass` deposit rule shares a
block's mass equally between the sj bounding it. Deleting that share left **18/18 gates passing** —
because the fixture carried ONE sj, so no block was ever bounded twice and `len(bounds)` was always
1. The rule was unfalsifiable on its own test set. A two-sj fixture plus a per-sj attribution
gate now catches it, and that gate is the *only* one that does: conservation cannot see the difference,
since giving a doubly-bounded block whole to each bound still sums to 1 when no such block exists.
⭐ **A fixture is an arm. Ask of it too: could this have failed?**

**prove-the-substrate. Prove the SUBSTRATE before proving the code.** For two milestones the simulated panel's
post-capture fragment-length distribution was byte-identical to its pre-capture one, and everything
measured against it inherited that. The tell was free the whole time: diff the two capture arms' truth
files. **When a simulated axis is the axis you are judging, gate the simulator on it.**

**can-the-benchmark-resolve-it. Before running a benchmark, prove it can resolve the axis you are changing.** A 32-condition suite
was used for months to judge a partition change it was structurally incapable of seeing: its fine region
set was row-for-row identical to its merged set. It also had `frag_std = 0`, was Poisson by construction,
and over-represented the terminus+splice-site boundary **12×**. `scripts/design/suite_resolves.py` is this
lesson made executable — every requirement scored against its *degenerate* value, no tuned thresholds.

**toys-rank-hotspots-backwards. A toy ranks performance hotspots BACKWARDS.** At 3.4 k vs 1.5 M regions: message relay 34 % → 81 %,
per-region solves 9 % → 29 %, the prior's EM **28 % → 0.7 %**. A whole analysis was spent on the toy's #1
hotspot. Profile on cached real data.

---

## B. Measurement and inference

**measure-the-ceiling-first. MEASURE THE CEILING BEFORE BUILDING THE CORRECTION.** Hand the consumer the *exact* answer for one
channel and see what perfecting it is worth. One channel was rank 0 for two sessions on the strength of a
coupling constant; the ceiling showed a **perfect** model of it bought ~1 %, while the channel nobody was
ranking was worth 21 %. An A/B tells you whether a change helped; a ceiling tells you whether the work is
worth starting — and it is available whenever the simulator writes truth.
*Instrument:* `scripts/design/calibration_truth_ab.py --ceiling`.

**a-scorer-scoped-to-the-mechanisms-targets. ⛔⛔ A SCORER RESTRICTED TO THE POPULATION A MECHANISM TARGETS WILL REPORT A WIN
THE WHOLE CHAIN REFUSES — and the sign can flip, not just the size.** Measured 2026-08-18: a relay arm
scored over the step-0 QUALIFYING exons read **0.916×** (a win) on unstranded × capture-OFF while the
same arm over EVERY live slot read **1.126×** (a loss) — the mechanism helped exactly the 9–18 % of exon
mass it was aimed at and hurt the rest. The scoped score is the right DIAGNOSTIC (it says whether the
mechanism does its own job) and the wrong VERDICT. This is `never-pool-the-strata`'s other axis: not
WHICH stratum, but WHICH SLOTS WITHIN one. Every verdict-level number is over all live objects; a scoped
number is labelled as a diagnostic or not reported.
*Sibling:* `never-pool-the-strata`. *Instrument:* `relay_pool_ab.py` scores all live slots for this reason.

**a-broad-population-carries-no-prior. ⛔⛔ A POPULATION PRIOR ONLY TRANSFERS INFORMATION IF THE POPULATION IS TIGHT.**
The gDNA landscape works because gDNA is near-uniform: pooling genuinely tells one object the library's level.
The same machinery fitted to RNA is worth **nothing** — measured 0.988 / 0.997 / 1.037 against base — because
the per-object RNA density spans ~3 decades, so the pooled distribution is wider than the per-object
uncertainty it was meant to reduce and the prior is flat where it matters. ⭐ *The rule:* before fitting a
population prior, measure the pooled distribution's WIDTH against the per-object uncertainty. If it is wider,
there is no prior to be had and the fit will read as neutral-to-worse. ⚠ The width was measured EARLY here
(2.4–3.0 decades) and read as "harmless"; it also means "uninformative", and that single reading would have
predicted the refusal before anything was built.

**attribution-must-survive-a-shuffle. ⛔⛔ AN ORACLE ARM THAT A SCRAMBLED ORACLE ALSO WINS HAS NOT PRICED THE TRUTH.**
An oracle RNA density fitted from the origin-split truth improved the panel 0.78–0.85×. Refitting it from the
SAME truth values SHUFFLED against their own opportunities — a shape wrong on purpose — **beat the true shape
at `g98` (0.786 vs 0.854)** and was neutral-to-worse elsewhere. So the gain was partly "any shape leaning
toward the vertex helps", not "the true density helps", and the ceiling did not mean what its name said.
⭐ *The rule:* before believing a ceiling built from truth, re-run it with that truth SCRAMBLED. If the gain
survives the scramble, the arm is measuring a property of the mechanism, not of the truth you fed it.
⚠ Cheap — one permutation — and it is the only check that separates the two.

**score-against-truth. Score against TRUTH, not against the previous run.** The simulator writes per-fragment ground truth
into the oracle BAM's read names. That is what found a score-annihilation bug, and what turned "the
deliverable improved" into "the deliverable got 23.9 % worse".

**zero-target-guards-are-one-sided. A zero-target guard is ONE-SIDED.** In a library with no gDNA, *any* change that lowers the
estimated gDNA fraction scores better. This reversed a published verdict once — a reported −13.1 % win
became +1.4 % worse over the full battery. Score on the contaminated conditions; use zero-target rows as
false-positive checks only.

**hard-labels-miss-soft-change. Hard-label metrics are nearly blind to soft changes.** A real change can move soft pool counts by
tens of thousands of fragments while a hard-label net is byte-identical. **Treat a byte-identical
hard-label result as no evidence, not as no change.** Related: never score pooled — one run had a −47.6 k
error and a +75.4 k error reading as a "nearly perfect" −1.9 k.
⭐ **And the cancellation is not bounded.** Scored per object against the oracle, `Σ|err| / |net|` reached
**274×** on one arm: the library-level figure read a +672 fragment error over a per-object Σ|err| of
184,271. **A net error is a lower bound on the real one, and it can be an arbitrarily weak one.** Report
`Σ|err|` beside every net, and the directional split beside both.

**never-pool-the-strata. THE DEFAULT INSTINCT IS A POOLED AVERAGE, AND IT IS WRONG HERE THREE WAYS.** (i) averaging error
over objects that cannot be solved buries the answer — but an exclusion is a promise to check that class
some other way (TRAPS: excluding-a-population-hides-it); (ii) a pooled mean hides a sign flip between strata; (iii) mass-weighting lets one
huge object set the number. ⭐ Report per stratum, with the denominator named.

⛔⛔ **AND SINCE 2026-08-14 THERE IS A FOURTH REASON, SHARPER THAN THE OTHER THREE: ONE STRATUM IS NOT A
DEVELOPMENT TARGET.** The 0.8.0 scope puts unstranded × capture-ON out of scope — DEFERRED, still
measured and still reported, but not what work is aimed at — while the three in-scope strata are
unstranded × capture-OFF, stranded × capture-OFF and stranded × capture-ON. The deferred cell carries
**64.5 % of transcript error and 90 % of gene-level error** on the rebuilt 16-condition ladder, so a
pooled total is mostly a report on the one stratum 0.8.0 does not ship, and an arm ranked on it is ranked
by the wrong number. ⚠ This is not the same objection as (ii): there the pooled mean *hid* a real signal,
here it faithfully reports a quantity nobody is optimising.

**a-threshold-on-a-fitted-residue. A BINARY CUT ON A FITTED PARAMETER'S RESIDUE IS NOT A POPULATION TEST — AND A BETTER THRESHOLD IS
NOT THE FIX.** `τ > 1e-9` promoted objects whose own statement was 1,377 nats wide into the scored
population, hiding a 1.06 M-fragment error. ⛔ A floor was implemented and refuted: τ is continuous across
the region, so any floor is a tuned constant. ⭐ The honest repair is to propagate the fitted parameter's
own variance, which drives the arm to ~0 by derivation (TRAPS: a-licence-with-no-floor).

**excluding-a-population-hides-it. ⛔⛔ EXCLUDING A POPULATION FROM THE DENOMINATOR WITHOUT GATING ITS OWN FAILURE MODE HIDES THE
LARGEST ERROR IN THE LIBRARY.** A condition reporting **mwae 0.0170** on a 39.5 %-solvable set had
**1,056,019 fragments** of error in the excluded class, with 6,425 objects moved a mean 0.414 from ½ and
100 % of them declaring a finite precision. The single worst mechanism was **0.0 % scored on every rung but
one**. ⭐ **The rule: every population you exclude needs its own gate, written at the same time as the
exclusion.** The tell is free — the *reason* to exclude a class is that its answer should be
uninformative, so measure how far from uninformative it actually is.

**name-the-observable-per-site. ⭐⭐ IF EVERY GATE FOR A CHANGE READS ONE INSTRUMENT, THE CHANGE IS GATED IN ONE PLACE — AND
`_relay`/`_transport` ARE TWO (now three, with `_flank_dom`).** Ten sound tests all read the combine's
publication; deleting the conjunct from `_relay` alone left **the entire calibration suite green**.
⭐ **For each place the change was made, name the observable that would move and confirm at least one gate
reads it.** A count of gates says nothing. ⚠ Re-opened by the very next change to the same pair, and again
by a third consumer nobody had listed (TRAPS: a-gate-that-reconstructs).

**starved-is-not-depleted. ⛔⛔ "STARVED" AND "DEPLETED" ARE NOT THE SAME DIAGNOSIS.** An object with almost no counts because
the experiment is too small is starved; one with almost no counts because the biology puts nothing there is
depleted. The first is fixed by more data, the second never is — and reading a depleted object as starved
threw away half a session chasing depth. ⭐ Decide which by asking whether the count grows with the lever
you have.

**the-substrate-knob-fought-back. A TILING LOOP CAN SUPPRESS THE VERY POPULATION UNDER STUDY.** Toy capture probes are written in
TRANSCRIPT space, so a probe spanning an internal sj offset has a genomic footprint in **two blocks**
— and `sim/capture/sampler._split_scale` then multiplies every gDNA fragment overlapping it by
`gdna_split_penalty`. Tiling across the whole transcript therefore put a split probe over every internal
sj and depressed exactly the fragments that span an `intron|exon` BOUNDARY. ⭐ The lesson is not "tile per
exon" (though that is the fix); it is that **a substrate knob can be adversarial to the object you are
measuring, and the way to find out is to read the sampler's weight for that object's own population** —
here, one `if len(blocks) > 1` in code nobody had reason to open.

**key-on-a-realised-quantity. ⛔⛔ AN RNA "LEVEL" THAT IS A MULTIPLE OF THE gDNA DENSITY IS NOT COMPARABLE ACROSS THE CAPTURE
AXIS.** At rung `m = 100` the exon's true `f_g` is **0.010 off capture and 0.31 on it**, because the
multiple is taken against the OFF-TARGET density while capture concentrates gDNA onto the exon ~65×. ⛔
**Key an operating point on a REALISED quantity you can see in the truth, never on a knob you set** — and
print the realised value so a row that cannot reach the target says so.

**price-the-halves-separately. ⛔⛔ PRICE THE HALVES SEPARATELY — A DIAGNOSED DEFECT AND A VALUABLE ONE ARE NOT THE SAME DEFECT.**
Pricing two length models together read −6.31 % and looked like one finding; split one pmf at a time it is
**−5.90 % gDNA / −0.43 % RNA**, and all of the gDNA half is capture-ON. The half that was fully diagnosed
was worth ~nothing and the half worth 14× more had not been looked at.

**panel-before-src. ⛔⛔ A TOY CEILING IS NOT A PANEL CEILING, AND THAT HAS NOW COST THREE BUILDS.** A re-solve ceiling
worth −0.009 on the toy carried to the 36-condition ladder as **+0.0013 (worse)**; a second change was
toy-positive and panel-negative the same way; a third improved its target rung 26 % and regressed the panel
axis that carries the mass. ⛔ **Run the panel arm before writing a mechanism into `src/`.** ⚠ And when the
panel disagrees, ask whether the defect is half of a cancelling pair before concluding (TRAPS: a-cancelling-defect-pair).

**capture-inverts-the-counted-side. ⛔ "THE WELL-COUNTED SIDE" IS NOT A FIXED SIDE — CAPTURE INVERTS IT.** Off capture an intron REGION
holds ~349 gDNA counts against a flanking BOUNDARY's 8–36; under capture the intron holds **1** and the BOUNDARY
holds 20–40, unchanged from 120 kb to 1.08 Mb. So a rule that transports from "the well-counted side"
silently reverses direction with the protocol, and transferring a DENSITY instead of a SHARE measured
**+0.207** on capture-ON × unstranded.

**admitting-an-object-costs. ⛔ ADMITTING AN OBJECT TO THE SCORED POPULATION IS A COST, AND A MECHANISM CAN DO IT SILENTLY.** A
prototype gave a class certainty it had not earned; `solv%` rose 43.1 → 48.2 and the edge-axis mwae went
0.0216 → **0.1449** — not because the answers got worse but because wrong answers started being counted.
⭐ Report `solv%` beside every accuracy number (TRAPS: excluding-a-population-hides-it inverted).

**substitution-understates-a-source. A CEILING BY SUBSTITUTION UNDERSTATES A MESSAGE SOURCE.** Replacing one object's answer with the
truth and re-scoring is the honest ceiling for a SINK — but a source's value is what it *carries*, and a
substitution does not propagate. Measured: substituting both `intron|exon` BOUNDARIES removed 9.1 % of the
gene's error, while the object they feed accounted for 82.7 %. ⭐ For a source the ceiling must be a real
arm — pin it and RE-SOLVE — and an instrument that offers substitution must say which of the two it is
doing (`toy_panel.py` now does).

**a-symptom-is-not-a-second-defect. ⛔⛔ A DOWNSTREAM NUMBER THAT IS WRONG BECAUSE ITS INPUT IS WRONG IS A
SYMPTOM, NOT A SECOND DEFECT — AND THE PROOF IS TO SUBSTITUTE ONLY THE UPSTREAM INPUT AND RE-RUN THE
*SHIPPED* DOWNSTREAM FUNCTION.** Measured 2026-08-13/14 on the effective-length shrinkage.

The shrinkage read as an independent bug, and a large one: at the ZERO-gDNA control the shipped factor
contracts **13,673 of 15,669** transcripts by a mean **0.345** where the correct answer is
exactly **1.000** — including on capture-OFF, where the module's own contract says the factor must be 1.
Every instinct says go and repair the shrinkage.

⭐ **It is not the shrinkage.** Substituting ONLY the composition arrays with truth and re-running the
SHIPPED function returns the right factor: `g00` capture-OFF **0.345 → 1.000**, and `g50` unstranded ×
capture-ON **0.834 → 0.401** against a truth of 0.401. `rho_ref` is fabricated entirely out of
false-positive gDNA. There is ONE split with two consumers and one function — `priors.py` imports
`_global_reference_density` from `capture_eff_length.py` — so repairing the composition repairs both, and
"repairing" the shrinkage would have installed a second error to cancel the first
(`TRAPS: a-cancelling-defect-pair`).

⛔ **The substitution is admissible here for a reason you must NAME before trusting the arm: the downstream
is a PURE FUNCTION of the substituted input, with no feedback.** Where the downstream feeds back — a message
source — the identical arm understates (`TRAPS: substitution-understates-a-source`) and the honest ceiling
is a pin-and-re-solve instead. Say which of the two you have before reading the number.

⛔⛔ **AND THE SECOND HALF IS WHY NOBODY HAD SEEN IT: A CEILING PRICES ONLY WHAT IS BUILT *AFTER* ITS
INJECTION POINT.** `effective_lengths_em` is built at `pipeline.py:816`, `assemble_priors` is called at
`pipeline.py:839`, and **every** measurement arm in the tree patches `assemble_priors`. The shrinkage was
therefore upstream of every ceiling ever run — identical in the baseline and in the ceiling, so it cancels
out of the difference and reads as worth nothing. **Every ceiling number in this project was measured with
a wrong ruler already installed.** ⭐ *The rule:* for every ceiling, write down WHERE the injection happens
and WHAT was already computed above it; anything above that point is not being priced, it is being assumed
correct. `TRAPS: measure-the-ceiling-first` says to price the ceiling — this says to check which stages the
ceiling can even see.

**a-locked-object-is-not-a-control. A STRUCTURALLY LOCKED OBJECT IS NOT A CONTROL.** An `intergenic|exon` boundary predicts
`f_g = 1.0000` exactly against a truth of 1.0, and that was read as the healthy twin of the broken
`intron|exon` boundary beside it — "the same object structurally, and the only difference is the sj
flux". It is not: it is **G1**, so it is never solved at all and keeps its pinned `{0,0,1}` init. It
cannot be wrong, so its correctness measures nothing. ⭐ *The control that works is the same class
split on the variable* — `intron|exon` boundaries **with** vs **without** sj flux — and it reversed
the conclusion: the class is 23 % low where **no** sj attaches, so the flux explains 26 % of its
error unstranded and **2.5–4 % on the strata where the class has own evidence at all.**

**draining-breaks-the-oracle. A per-fragment-independent partition stops being one the moment a downstream step conditions on
the whole tally.** The oracle is trustworthy because the accumulator deposits each fragment
independently, so splitting the BAM by origin and re-scanning must reconstruct the full payload exactly.
The second pass breaks that: its multinomial is scored against the payload's *own* densities, so three
partitions drained separately are not the whole drained. **The identity that makes a truth source valid
also fixes where in the pipeline it can be taken** — and the honest response is to measure the undrained
stage rather than to drain the parts and hope.

**an-equal-length-panel-defeats-the-lift. ⛔⛔ THE REPAIR FOR THE RULE ABOVE HAS A SUBSTRATE IT CANNOT WORK ON, AND IT IS THE
MAIN PANEL.** `lift_choices` restores the identity by replaying the whole library's already-drawn hypotheses
inside each origin partition, and its own docstring states the one ambiguity: records that tie on the
deferred bank's canonical key cannot be told apart, so the assignment is greedy and the count is returned
rather than swallowed. On **distinct-span** substrates that count is 0 and the lift is exact.
⛔ **The ladder is built with EQUAL configured fragment lengths, which makes span collisions the common case
by construction.** Measured over all 36 conditions: **265,781 of 6,774,503 held fragments (3.92 %)** are
ambiguous, and replaying deposits **22,518 spliced records inside the gdna partition** — which
`OracleTruth._validate` refuses, correctly, because gDNA does not splice.
⭐ **So the refusal is the finding, not an obstacle to route around.** The right response was to run both
sides UNDRAINED and then *price* the caveat instead of asserting it was small: re-deriving the shipped prior
on the drained payload moves it by **0.153 %** (gDNA arm) and **0.462 %** (RNA arm), against effects of
2.5–94 %. ⚠ A returned-and-reported ambiguity count is what made this a two-minute diagnosis instead of a
silently biased truth source — the same argument as `waive-with-a-measurement`, one layer down.
⚠ **Stamped 2026-08-14:** the counts above belong to the 36-condition ladder, since rebuilt to 16
conditions — a record, not a restatement. ⛔ **And equal lengths are NOT the panel's mistake**: read
`TRAPS: a-length-gap-bypasses-calibration` immediately below for what they buy and why the ladder keeps
them. The lift is what gets priced, not the substrate.

**a-length-gap-bypasses-calibration. ⛔⛔⛔ A LARGE gDNA-vs-RNA FRAGMENT-LENGTH GAP LETS THE EM ASSIGN
FRAGMENTS ON LENGTH ALONE, SO THE TOOL ANSWERS CORRECTLY WITH CALIBRATION BROKEN — AN EQUAL-LENGTH PANEL IS
WHAT EXPOSES THE BUG.** Owner, 2026-08-14. This is the gDNA ladder's real design reason, and the one that
had been written down is the weaker one.

**The EM ALREADY READS THE FRAGMENT-LENGTH DISTRIBUTION.** `second_pass.length_likelihood` scores every
hypothesis under the pmf appropriate to it, so a well-separated `μ_g` and `μ_r` separate the two populations
at the FRAGMENT level, before calibration's composition claim is consulted at all. The transcript numbers
then come out acceptable **however wrong that claim is**. Calibration is not being exercised on such a
panel; it is being carried.

⭐ **So the ladder's EQUAL fragment lengths are a FORCING FUNCTION, not a simplification.** With the means
equal the length channel carries exactly zero composition information
(`TRAPS: equal-lengths-carry-no-composition`), which leaves strand, density and — when it is switched on —
propagation across objects. Those are the code under test. ⛔ The weaker reason the docs used to give —
*"the length channel is neutralised, so the residual error is attributable to density and strand"* — is a
CONSEQUENCE, not the reason. Attribution is convenient; **not being lied to is the point.**

⛔ **The tell that a gap is carrying you, and it is free:** score the CALIBRATION result against ORACLE
CALIBRATION and print it beside the end-to-end number. A healthy end-to-end number sitting on top of an
unhealthy composition claim is the signature. Measured 2026-08-13/14: on unstranded × capture-ON the tool
emits a near-constant near-ZERO gDNA fraction regardless of truth — exon `f_g` **0.040 / 0.0016 / 0.0021**
against truths **0.054 / 0.518 / 0.982** — so at the lowest gDNA rung it looks correct by coincidence
(`TRAPS: a-single-level-panel-cannot-see-a-constant`, and it is why three gDNA levels is the minimum).

⭐ **The general form, and it is settled before the panel is built: for each axis a substrate varies, ask
which STAGE that axis lets a LATER stage answer WITHOUT.** A benchmark that hands the deliverable a
shortcut measures the shortcut. ⚠ It is the mirror of `TRAPS: an-equal-length-panel-defeats-the-lift`
directly above — equal lengths COST the drained oracle lift and BUY the only substrate on which calibration
can be caught being wrong. Both are true at once, and the panel keeps its equal lengths.

**weight-it-like-the-consumer. A bp-weighted mean and a fragment-weighted mean answer different questions.** An "11 % over-call"
was a bp-weighted geometric mean; the estimator is fragment-weighted, 89 % of the mass sat where the
effect is inert, and end-to-end it moved the answer by ≤ 0.0002. Both statements were true; only one
decided anything. **Weight the average the way the consumer weights it.**

**a-support-ceiling-is-the-clamp. A SUPPORT CEILING THAT MATCHES THE CLAMP IS NOT A MATCH — it is the clamp.** A distribution's
support "agreeing" with `max_frag_length` was recorded as a fix; the narrower estimate had been correct
and the "fix" was an uncut intron.

**log-variance-is-not-linear. Every "is the declared precision earned?" number written before 2026-07-28 compared a LOG-space
variance against a LINEAR-space error.** `Var(log f)` is not `Var(f)`; the delta method
(`Var(f) ≈ f²·Var(log f)`) converts, and at small `f` the two differ by orders of magnitude. Every
overconfidence figure from before that date is void, not merely imprecise.

**re-record-the-baseline. A delta is only attributable if its baseline came from the same tree in the same session.** Re-record
the before-picture, do not quote a stored one.

**a-truth-table-of-aggregates. ⛔⛔ A TRUTH TABLE'S LABEL COLUMN MAY HOLD NESTED AGGREGATES, AND
SUMMING IT DOUBLE-COUNTS — WITH THE CONTROL ARM UNAFFECTED, WHICH IS THE MOST PERSUASIVE SHAPE A WRONG
RESULT CAN TAKE.** Found 2026-08-09. `truth_fragment_lengths.tsv`'s `kind` column is
`{mrna, nrna, gdna, rna, all}`: the first three are POPULATIONS, `rna` re-counts (mrna+nrna) and `all`
re-counts the library. A parser written as `"gdna" if kind == "gdna" else "rna"` summed **mrna + rna +
all** into one bucket — mRNA twice plus the library once — and reported a "true RNA library" of 191.06 bp
where mRNA is **212.20**.

⛔ **The gDNA arm was untouched**, because only one kind matched it. So the control read exactly 1.000 and
the treatment looked broken, and the conclusion "the RNA pmf is +10.7 % long, `pi(w)` is 2× too flat"
survived TWO rounds of interpretation and two doc commits before a premise gate on an unrelated script
caught it. ⭐ Nothing about the analysis was sloppy; the input was.

⭐ **Two defences, and the second is the general one.** Enumerate the membership sets explicitly and RAISE
on an unrecognised label — guessing a bucket is the whole bug. And **put a premise gate on any truth
source**: the script that caught this compared read-name spans against the same table's mean and printed
both, which is why a 21 bp disagreement became visible instead of being absorbed.

**a-single-level-panel-cannot-see-a-constant. ⛔⛔ A PANEL THAT HOLDS THE QUANTITY OF INTEREST FIXED
CANNOT DISTINGUISH A GOOD ESTIMATOR FROM A CONSTANT THAT HAPPENS TO EQUAL IT.** Measured 2026-08-10, and
it nearly landed a feature.

`length_likelihood` was A/B'd on the flgap pair — 8 conditions varying fragment length (±40 %), capture
and strand specificity — and came back better on **7 of 8**, mean `|f_gdna − truth|` **0.133 → 0.0175
(−87 %)**, with the long-standing unstranded × capture-ON blindness apparently resolved (0.0324 → 0.5222
against a truth of 0.507). The one axis the flgap pair does NOT vary is the gDNA level: every condition is
`g50`, truth ≈ 0.5.

⛔ **On the `g00` zero-gDNA control the channel reports 54–57 % gDNA in a library containing none**, and at
`g98` (truth 0.980) it reports 0.287. Unstranded, its answer is **0.539 / 0.522 / 0.287** as truth goes
**0.00 / 0.51 / 0.98** — a near-constant near ½. The 87 % win was the panel agreeing with the constant.

⭐ **The tell is structural and available before the run: list the axes the panel varies, and check that
the DELIVERABLE is one of them.** A panel can be rich in every other dimension and still be blind, and
"7 of 8 conditions improved" reads as breadth when it is one point measured eight ways.
⚠ It is the mirror of `TRAPS: an-equal-length-panel-defeats-the-lift` — there the panel held the *input*
fixed, here it holds the *answer* fixed — and of the standing rule not to score on the zero-gDNA arm
ALONE. Both arms are required because each is one-sided on its own.

⛔⛔ **DEFERRED POST-0.8.0** (owner, 2026-08-14). The `length_likelihood` above is the fragment-length
COMPOSITION channel for CALIBRATION. It never shipped in `src/`, and 0.8.0 ships without it: **do not
propose it, do not rank it, do not list it as a candidate.** This entry is the record of why, not an
invitation to rebuild it. ⚠ It is a DIFFERENT thing from `second_pass.length_likelihood`, the per-fragment
assignment factor, which ships and is untouched by the ruling. ⭐ The PANEL lesson — a panel that holds the
deliverable fixed cannot tell an estimator from a constant — is general and is why this entry exists.

**a-moment-match-is-not-sufficient. ⛔ MATCHING A MODEL'S MOMENTS TO ITS POPULATION'S CAN MAKE THE ANSWER
WORSE, SO A MOMENT RATIO IS A NECESSARY AND NOT A SUFFICIENT GATE.** Measured 2026-08-10 on the length
channel. A derived opportunity took the RNA crossing moments from `realised/predicted` 0.925 to 0.971 off
capture and 0.948 to 0.983 under it — an improvement on BOTH arms. The channel's accuracy against truth
went 0.307 → 0.177 off capture (a large win) and **0.042 → 0.080 under it** (a regression). ⭐ The
likelihood reads a JOINT function of both components' moments and their covariance; improving one
component's marginal moment can move the pair apart. ⛔ So a table of moment ratios ranks candidate
repairs; it does not certify one. Score the deliverable.
⛔ **DEFERRED POST-0.8.0** (owner, 2026-08-14): the fragment-length COMPOSITION channel this was measured on
is not a 0.8.0 candidate — do not propose it, rank it or build it; it never shipped in `src/`. ⭐ The moment
lesson is general and stands without it.

**score-the-consumers-own-count. ⛔⛔ TRUTH IS NOT A YARDSTICK UNTIL YOU NAME THE CONSUMER — A TRUTH
NOTHING READS MEASURES YOUR CHOICE OF TARGET, AND IT WILL DO IT IN THE RIGHT UNITS.** Found and fixed
2026-08-08. `prior_vs_oracle.py`'s `F` arm was the per-locus count of gDNA fragments whose FIRST BASE
lands in the locus, taken from `node_start_count` — the accumulator's one exact invariant, one deposit
per accepted fragment, projected through the shipped projection. Impeccable provenance, in fragments,
and its docstring called it "EXACT on the gDNA arm ... with nothing to subtract".

⛔ **The EM counts something else.** `n_gdna` in `apply_grouped_prior_update` is the soft count over the
multi-locus's own EM UNITS, and a fragment becomes one by being a scored CANDIDATE — which a fragment
that starts in the intergenic flank and reaches into the locus is. `F` drops exactly that straddling
population: **0.35 %–3.35 %** of the total, per condition.

⭐⭐ **And the consequence was not a 3 % error in a number, it was a phantom research programme.** Scored
against `F`, the prior assembler's residual with perfect masses AND perfect per-component shares read
`rel` **0.0035–0.0302**, of which the shares explained only 15–43 % — so ~72 % was recorded as an open
residual with two candidate mechanisms and a place in the ranked plan. Re-scored against the EM's own
count it is **2.8e-5 … 2.0e-3**, the shares explain **82–99 %**, and *there is no residual to explain*.
**The 72 % was the yardstick.**

⚠ **AND IT RECURRED ON 2026-08-11, IN THE SAME FILE THAT ALREADY DOCUMENTED IT.** A new pool-level table
in `quant_accuracy.py` scored `gdna_est` against `gdna_true` and read **−50.7 %** panel-wide — a
spectacular, publishable gDNA under-call. `gdna_est` is `gdna_em_count`, which EXCLUDES fragments that
reached no locus; `gdna_true` is the simulator's origin count, which includes them. **Off capture more
than half of gDNA is intergenic**, so the "defect" was the missing pool and nothing else — with them the
same panel reads **−0.6 %**. ⭐ `score_library`'s own docstring, four lines above the field being
misused, records the identical mistake fabricating an off-capture EM under-call (0.3151 against a truth
of 0.5000). ⛔ *The lesson on top of the lesson:* **a documented trap does not protect the NEXT consumer
of the same field.** The columns not summing to the library was the tell, and it was visible on the
first render.

⛔ **The tell, and it is available before any measurement: write down which LINE of the consumer reads
the number.** If that line's population is a different set than the truth arm's population, the arm is
measuring the difference between the two sets and calling it error. Here the consumer was one C++
function and the answer took one grep. ⚠ Provenance is not the check — `F` had the best provenance
available and was still the wrong quantity; "exact" was a claim about the *bank*, not about the *target*.

⭐ **The repair generalises: take the target from the consumer's own bookkeeping.** `Fo` is
`MultiLocus.unit_indices` — the array `locus_partition` scatters by — labelled with each fragment's true
origin, so the locus assignment is the EM's own and nothing is re-derived. ⛔ Check it is not circular:
`assemble_priors` never reads a unit count, and that is gated behaviourally, not by grep.

⛔⛔ **AND THEN THE SAME COMMIT THAT WROTE THIS RULE BROKE IT AGAIN, ON THE OTHER ARM — WHICH IS THE
PART WORTH REMEMBERING.** The corrected instrument reported a "new defect": `a_r` allegedly tilted the
prior's composition claim by **+0.07…+0.10** because it withholds spliced RNA while the EM's `n_rna`
counts it. It had scored `a_g/(a_g+a_r)` against a truth counting **every** RNA unit. But a spliced unit
never receives a gDNA candidate (`em_solver.cpp`: `has_gdna = !is_spliced && isfinite(gdna_ll)`), so
spliced RNA does not compete with gDNA and a prior that arbitrates that competition must not count it —
doing so would penalise gDNA with fragments it could never have won. Against the population `a_g:a_r`
actually describes, the claim is exact to **≤ 5e-4**. ⚠ The right denominator was already a column in
the same table.

⭐ **Why it recurred, and the operational form of the rule.** Naming the consumer is not enough: the
consumer reads several arrays and each has its own population. `a_g` competes over the unspliced pool,
`a_r` over the unspliced pool, `n_rna` spans the locus — three populations in one function. ⛔ So the
check is **per array, not per instrument**: for each number the consumer reads, write down the SET of
fragments it is a count of, and score it against that set. An instrument that has got one arm right is
not thereby right about the next one.

---

**the-floor-must-reproduce-the-selection. ⛔⛔ A NOISE FLOOR DRAWN UNCONDITIONALLY UNDER A SCORER THAT
CONDITIONS READS THE CONDITIONING AS A RULE ERROR.** Measured 2026-08-19 in `hop_currency.py`. A hop is
scored only if its destination holds mass; at capture-OFF `g05` a boundary expects ~1.2 gDNA crossings,
so the scored destinations are a ZERO-TRUNCATED Poisson — realised mean 1.64, variance 0.71 — while a
perfect LEVEL rule claims the untruncated 1.17. The realised error was 972 fragments against a floor of
708 drawn with an unconditional redraw, and the 37 % gap printed as a 10-14 % "rule error" into every
small-boundary class at `g05`. The rule was exact; the null was wrong.
⭐ *The rule:* whatever selection the scorer applies to the data, the null must apply to its draws — here a
pure-gDNA destination is redrawn conditional on `n >= 1` (inverse CDF), and the same three classes read
"either" at their floors (972 vs 1,011). A floor that cannot be falsified by turning the conditioning off
is not a floor; `hop_currency.py --self-test` carries both arms.
*Sibling:* `a-ratio-needs-a-population-that-can-supply-its-numerator` (a selection on the denominator),
`excluding-a-population-hides-it`.


## C. Pools, selections and divisors

**two-divisors-opposite-sign. ⭐⭐ TWO DIVISORS BUILT FROM ONE pmf CAN STILL DISAGREE — IF THEY RESPOND TO IT WITH OPPOSITE SIGN.**
`E_J = E[w]−1` RISES with the mean fragment length while `E_r = e−E[w]+1` FALLS, so a length-model error
appears as a junction-vs-exon frame gap of `Δ·(1/E_J + 1/E_r)` = **0.62 %/bp**. The two are exactly
consistent (202.8 + 797.2 = 1000.0) and still disagree by 10 %. ⭐ **When two quantities built from one
model disagree, differentiate BOTH with respect to that model before looking for a bug in either.**
⛔⛔ And the second term added here was WRONG, which is the sharper lesson: a "finite-transcript placement"
factor of 1.024 was sound arithmetic about a generative model **the simulator does not use** — it reweights
the length marginal by the same opportunity, so the factor cancels exactly. ⭐ **Before pricing a selection
effect in a simulated substrate, read the simulator's own sampling code — not the docstring, the code.**

**frame-free-is-not-assumption-free. ⭐⭐⭐ A TERM THAT IS FRAME-FREE IS NOT THEREFORE ASSUMPTION-FREE — LOOK AT THE NUISANCE YOU
PROFILED OUT, NOT ONLY AT THE DIVISOR THAT CANCELLED.** `boundary_spliced` is certified RNA, so the obvious
move is a coefficient `S` on ψ's RNA arm: `E[S] = c·(1−f_g)·M`, and every opportunity ratio lives in `c`,
which multiplies the MEAN and is therefore an additive constant in log space. That reasoning is correct —
TRAPS: two-divisors-opposite-sign's opposite-sign trap is structurally absent, neither divisor reaches the retained term — and it is
**not enough**, because the *same* `c` contains the splice-visibility `q`, whose dropped term
`−(q/(1−q))(1−f_g)M` is the same size as the term kept. ⛔ With `q` free the profile likelihood in `f_g`
is **exactly flat on [0,1)**: the whole channel carries ONE BIT, "`f_g ≠ 1`". Measured on all 36 ladder
conditions against origin-split truth, the raw-count term is **worse than the uninformative reference on
12 of them**, worst **+0.4578** mwae. ⭐⭐ **The tell was free and general: its benefit tracked the ANSWER,
not the evidence** — −0.49 at `g00` where the truth is all-RNA, +0.45 at `g98` where it is all-gDNA. **A
channel whose sign follows the truth is a prior, not information** (TRAPS: honesty-metrics-reward-ignorance's shape, and the two zero controls
are what made it visible). ⭐ Ask of any "one-sided floor": *is the term I am dropping bounded, or is it
the dominant one?* `certified_q_census.py` · `test_certified_rna_licence.py`.

**a-purity-filter-is-a-length-filter. A PURITY FILTER ON A LENGTH POOL IS A LENGTH FILTER.** Barring fragments whose length was partly
*inferred* rather than sequenced selects exactly the ones whose mates sit far apart, so the pool became
**−9.6 % mean / −22.5 % sd** against truth where keeping them read +0.7 % / +2.4 %. **Before excluding a
population from a pool, ask what the exclusion criterion correlates with.** If it correlates with the
axis being measured, purity and accuracy point in opposite directions.
*The rule that replaced it:* a pool is keyed on **determinacy, not provenance** — a fragment enters when
exactly one hypothesis survived, however it got there.

**pure-and-length-censored. A POOL CAN BE COMPOSITION-PURE AND LENGTH-CENSORED AT THE SAME TIME, and the second is invisible to
every purity argument.** "gDNA contained in an intergenic or intronic region" is 100 % gDNA by construction
and, under hybrid capture, ~15 % short — because a long fragment beside a probe *reaches* the exon
boundary and stops being *contained*. It reads as contamination: the pools it spills into resemble the RNA
pool and were filed for two milestones as "not gDNA". They are gDNA. **Ask what a pool's selection rule
correlates with, not only what the selected fragments are.**

**divide-by-a-probability. A POOL DIVIDED BY ITS OPPORTUNITY MUST BE DIVIDED BY A PROBABILITY, NOT BY A COUNT.**
`count(w)/A(w)` recovers the distribution lengths were **drawn** from; every consumer needs the one the
library **realises**, which is the drawn one weighted by how many placements each length has. So the
divisor is `A(w)/T(w)`. The two forms differ in *shape*, and the ratio is also where an abundance
weighting cancels — swept to pathological regimes, the ratio form is never worse than not correcting and
the `A`-only form is.

**opposite-tilts-must-not-pool. POOLS WITH OPPOSITE TILTS MUST NOT BE POOLED RAW.** A *contained* pool's opportunity falls with
length; a *crossing* pool's rises. Summing the histograms and applying one divisor read a gDNA mean of
**146.05** where the contained pool alone said **88.0**. ⭐ Summing the counts **and** the matching
per-pool opportunities is a *different* operation and is correct — it is the opportunity-weighted average
of the per-pool estimates, and under Poisson counts those weights are exactly inverse-variance.

**a-mean-of-ratios-inherits-the-partition. ⭐⭐ AN ESTIMATOR THAT AVERAGES PER-OBJECT RATIOS IS A FUNCTION
OF WHERE THE ANNOTATION DREW ITS BOUNDARIES.** A region boundary appears wherever a signature changes — an
antisense feature overlapping on the other strand will split a uniform stretch in two — so how many
objects a transcript's path contains is an artefact, not biology. Averaging densities (harmonic,
geometric or arithmetic) weights each object by that artefact; `Σmass / Σopportunity` does not, because
splitting an object splits both sums together. Measured while designing the per-transcript RNA prior:
subdividing one 20 kb intron from **1 region to 4 to 10**, total mass and total opportunity held constant,
moved a shadow span's share of its locus's prior from **9.4 % → 18.6 % → 32.4 %**; the pooled form gives
the same answer at 1, 4 and 10. ⭐ `node_geometry` already states the rule for a different pool —
*"the ratio of sums, never the mean of ratios"* — and it generalises: **an estimator must be invariant to
a re-partition that leaves the data unchanged, and that is a property you can test directly.**

**a-trap-names-the-defect-not-the-repair. ⛔⛔ CITING A TRAP DOES NOT LICENSE A MECHANISM — THE TRAP NAMES
THE DEFECT, NOT WHICH QUANTITY TO REPAIR.** `TRAPS: a-zero-count-is-a-measurement` says a zero count over
a large exposure is a strong statement and must not be read as no data. From that, a per-transcript prior
was designed around the Jeffreys posterior MEAN as a density LOCATION, `(mass + ½)/opportunity`, and the
trap was cited as licensing it. ⭐ **It licenses no such thing: the same trap's shipped repair put the
half on the PRECISION (`trigamma(n+½)`), which measured +0.0 % with the zero-control rows byte-identical,
and the LOCATION form is in `ROADMAP.md` §4.1's graveyard at +7,269 % on `g00` (a second form at
+96,299 %).** `node_init.py` states the arithmetic in one line — three components each gaining `+½` breaks
`Σ_c ρ_c·E_c = M` by exactly 3/2 — and the sentence to remember is **"the half belongs to the rate's
VARIANCE, not to a share of a total."** ⛔ The rule: before adopting a mechanism because a trap seems to
motivate it, grep the GRAVEYARD for the mechanism. Eleven were priced and refused; a twelfth was
re-derived from first principles and was already row one.

**a-stale-gate-accuses-the-newest-change. ⛔⛔ A GATE WHOSE PREMISE EXPIRED DOES NOT GO QUIET — IT FAILS,
AND IT BLAMES WHATEVER IS IN FLIGHT.** `rescan_panels.py` gates an IRREVERSIBLE cache rebuild on
byte-identity, stating the reason in its own comparator: *"Byte-identity, never a tolerance: **these are
integer tallies** and the whole claim is that they did not move."* They were, under the fixed point. The
2026-08-10 owner ruling replaced it with one numeric convention — a COUNT is an integer, a FRACTION is
float64 — and float addition is not associative across worker threads, so six banks stopped being
bit-reproducible **that day**. `scan_payload`'s docstring recorded the consequence and even the remedy
("tests validate the float banks within a DERIVED tolerance"); this gate was never told.
⭐ **So on 2026-08-13 it failed a schema change on five banks that change had not touched** — 2,419 of
35,135 elements here, 3,930 of 13,482 there — and the natural reading, *"my change had a side effect"*,
was wrong. ⛔ The sibling `a-guard-outlives-its-divisor` is the INERT direction and is the gentler one: a
gate that stops testing merely lulls you, while one that cannot pass sends you hunting a defect that does
not exist, at the exact moment you most need to trust it.
⭐⭐ **WHAT RESOLVED IT WAS A CONTROL ON THE INSTRUMENT, NOT ON THE CHANGE**: scan the same BAM twice with
the same binary and compare. They differed on exactly those six banks, by at most 3.5e-14 relative — so
the gate was unsatisfiable for anyone, and had been for three days of commits. **When a gate fails on
something you did not touch, first ask whether it can still pass at all.**
⚠ And the exposure is structural: a gate that runs only at an irreversible moment is a gate nothing
exercises in between, so its rot is invisible until the one day it matters. `--self-test` existed and
passed 11/11 throughout — it perturbed the comparator against a fixture of its own, never against the
data the comparator now had to judge.

**an-upper-bound-is-not-an-estimate. ⛔⛔ A BOUND THAT IS SHARED IS NOT EVIDENCE ABOUT WHO SHARES IT.**
The per-transcript RNA prior's founding theorem is sound: mass at an object is shared by every transcript
covering it, so no transcript can be denser than the thinnest place on its own path, and its density is
bounded by the MINIMUM along that path. ⭐ The bound is true. It is also **useless wherever the minimum is
attained at an object the transcript does not own** — and measured on `g00 ss0.99 capture_off`, **3,644 of
4,839 silent transcripts (75.3 %) share at least one object with an expressed one**, so they inherit the
loud neighbour's bound and the prior asserts them into existence. Transcript false-positive mass went
**1.76–2.20× WORSE than base on all twelve rungs of the family**, while the *same weights* took gene-level
error to 0.40–0.53×: the bound is informative about the UNION and silent about the split.
⛔ **The tell is that the zero-weight SET was byte-identical across all twelve arms** — four soft-min modes
× three multipliers changed nothing about *which* transcripts were zeroed, because a bound is zero only
when every object is, which is a property of the data and not of the estimator. **When an estimator's
support does not move as you vary it, you are varying the wrong thing.** ⭐ Retreating to GENE granularity
did not rescue it (1.340× against 1.317×), which is the confirmation: the damage was never the within-gene
split. Priced and refused in `ROADMAP.md` §4.

**a-gate-on-the-helper-is-not-a-gate-on-the-caller. ⛔⛔ TESTING THE FUNCTION THAT COMPUTES A QUANTITY DOES
NOT TEST THE CALL THAT ASKS FOR IT.** Stage 6's weight builder shipped with 29 green gates and **two holes,
both found by perturbation and neither findable any other way**. (1) The re-partition gates called the
power-mean helper directly, so replacing the CALLER's opportunity weights with `ones` — the exact defect
`TRAPS: a-mean-of-ratios-inherits-the-partition` exists to prevent, and the whole reason the module
deviates from its spec — fired **nothing**. (2) The two opportunity arrays were compared in their own test,
so making `opportunity="total"` silently serve the UNSPLICED array fired **nothing** — because the
derivation gate used a SINGLE-exon transcript, where the two are equal by construction.
⭐ Both repairs are the same shape: a gate that goes through the PUBLIC entry point on a fixture where the
wrong answer is numerically different. ⚠ And (2) is `could-the-arm-have-fired` wearing a different hat —
the fixture was drawn from the majority population (single-exon transcripts are most of the annotation by
count) and the majority population is exactly where the two branches coincide.

**fractional-mass-is-the-problem. Fractional mass IS the partitioning problem.** A fragment spanning 4 regions writes six fractional
numbers whose values depend on region sizes, purely because a mass is conserved; the same fragment's
three crossing counts depend on nothing. *Corollary:* multimapper and ambiguous-path assignment must stay
INTEGRAL, or the non-integer observable returns and the count stops being a count.

**conservation-misses-mis-attribution. Mass conservation does not catch mis-attribution.** One fragment can credit the same boundary side
twice, and credit a boundary it never crossed, with total mass still exactly 1.0.

⭐⭐ **AND IT IS FAR WORSE THAN "a defect can slip through": AN ENTIRE ALTERNATIVE RULE CONSERVES.**
Measured 2026-08-08 while landing `boundary_unspliced_mass`. The shipped deposit shares a slice's
`slice_len / L` between the boundaries that bound it; the `1/K` rule an earlier design draft proposed gives
every crossed boundary an equal `1/K`. **Both sum to exactly one per fragment.** Injecting `1/K` into the
specification left *every* conservation gate green — both per-fragment laws, the exhaustiveness sum over
5,193 enumerated fragments, and even the closed form at a boundary whose flanks exceed every fragment, where
`K == 1` makes the two rules identical.

⭐ **What caught it was a re-derivation on a different axis**: attribute each fragment BASE to the boundaries
bounding its own region and divide by `L`. `1/K` is not expressible that way — a base carries no knowledge
of how many boundaries the whole fragment crossed — so the two rules separate immediately.
⛔ **The general form: if a property is invariant across the rules you are choosing between, a gate on
that property is not a gate on the choice.** Conservation was the obvious thing to test and it was the
one thing that could not decide anything. `tests/native/test_conserved_mass.py` states the four claims
separately for exactly this reason.

⚠ And a second, cheaper trap in the same file: the law is exact in RATIONALS, while the bank is a fixed
point. The first version of the conservation test asserted exact equality against a *quantised* array
and failed on its own correct code. Those are two claims — the rule conserves, and the representation
costs ≤ K quanta — and each needs its own assertion. Collapsing them either hides a real error inside a
rounding tolerance or reports rounding as a defect.

**a-guard-outlives-its-divisor. ⛔⛔ DELETE THE DIVISOR AND THE GUARD AGAINST IT GOES INERT — WHILE ITS TEST
KEEPS PASSING.** Measured 2026-08-08, re-basing `assemble_priors` onto the conserved fragment count.
`_mass_where_there_is_opportunity` exists because `rho = Σm/ΣS` is a rate: mass in the numerator with no
exposure in the denominator inflates it, and the floored variant `mass / max(S, 1e-9)` reaches ~1e9. The
rewrite removed the division — the prior now reads the count out of the bank — and **the guard's own
perturbation test went on passing, for a reason that had nothing to do with the guard**: its fixture
asserted `prior == 0`, and the mass really was zero on the object it was checking.

⭐ **Three separate things had quietly become true, and none was visible from the test result:**

1. The guard was **inert on the path the test exercised**. Injecting the floored divisor changed nothing
   there, so the test could not have failed.
2. The guard was **still load-bearing on a different path** — the eff-length, which still divides by
   `ρ_ref` — where nothing gated it at all. Injecting the floor there moved `gdna_eff_len` by **+19.97 bp**
   at 20 units of stray mass and **+44.93 bp** at 5,000, pinning against the boundary ceiling.
3. The correct behaviour on the first path had **reversed**: a count has nothing to divide by, so
   dropping zero-opportunity mass now silently loses fragments the accumulator really deposited. The
   test's assertion was not just unfalsifiable, it was backwards.

⛔ **The tell is structural, not statistical: a guard's justification names a specific operation** — here
a division — **so when that operation leaves the function, re-ask whether the guard still has a subject,
and where its subject went.** A green test is no evidence either way; only injecting the defect is.
⚠ The two halves are now a NAMED PAIR that cite each other
(`test_prior_units.test_mass_on_a_zero_opportunity_object_STILL_COUNTS_because_a_count_has_no_divisor`
and `test_priors.test_stray_mass_on_a_zero_opportunity_line_is_dropped_from_the_eff_len`), because either
one alone reads as a ruling about the whole file.

**a-fold-grows-a-heuristic. ⛔⛔ A QUANTITY FOLDED ONTO ANOTHER AXIS TO FIT A CONSUMER'S INTERFACE WILL GROW A
HEURISTIC TO REPAIR THE FOLD — AND THE HEURISTIC WILL OUTLIVE EVERYONE'S MEMORY OF WHY THE FOLD WAS THERE.**
Found 2026-08-08. A contiguous BOUNDARY is a 0-bp boundary; `_project_regions_to_loci` distributes by
`overlap / region_size_bp` and so cannot see an object with no extent. Rather than give the projection a
point path, each boundary's mass was **folded into one flank region** — and `_left_keyed_edge_arrays`' own
docstring says exactly that: *"all that remains is to hang it off a region so the genomic-overlap projection
can reach it."*

⭐ **Then the fold acquired a patch.** Keying left loses a locus's far-LEFT boundary into its intergenic flank,
which the projection drops — so an intergenic RE-KEY was added to send that one boundary right. The patch is
what made the attribution worst: it routes **100 %** of a boundary into the gene when a median
**8,066 bp** intergenic flank against a **211 bp** locus flank means ~61 % of that boundary's mass sits outside.

⛔ **Three tells that you are looking at one of these:**

1. A helper whose docstring justifies itself by a *consumer's* limitation rather than by the model.
2. A special case that exists to undo the general case (the re-key undoes left-keying).
3. Two callers of the same object for different purposes, where the fold is right for one and wrong for
   the other. Shipped v0.7.1 pooled a crossing's two sides because *"the two halves are one physical
   crossing event"* — true for the `min(m/ρ_ref, S)` CONTRACTION, false for the fragment COUNT, and one
   object served both.

⭐ **The repair is never a better heuristic; it is to give the consumer the axis it was missing.** Projecting
a boundary AS a boundary deleted `edge_owner_nodes`, the re-key, `_component_node_arrays`,
`_mass_where_there_is_opportunity` and `_left_keyed_edge_arrays` at once, needed no schema change, and was
numerically identical on real data — the fold had been obscure rather than wrong, which is precisely why it
survived. `DESIGN.md` §3.1b.

⚠ **And the fold hid a gate hole in its own axis**: `node_right_edge[r] == r` on a single reference, so every
fixture was blind to the difference between a boundary index and a left-node index. The replacement's index
conversion was *untested* until a perturbation said so (`TRAPS: perturb-every-gate`).

---

**a-ratio-needs-a-population-that-can-supply-its-numerator. ⛔⛔ ADDING OPPORTUNITY THAT STRUCTURALLY CANNOT
CARRY COUNTS MANUFACTURES A DEFICIT OF EXACTLY ITS OWN SHARE — and it reads as a defect in the estimator.**
Measured 2026-08-19, and it cost a false finding that was reported before it was caught. An ad-hoc A/B of
`calibration_oracle`'s field gate pooled every reference; the shipped gate calls `gate_uniformity` ONCE PER
gDNA-BEARING REFERENCE (`gdna_refs = {r : n_g > 0}`). Pooling folded the panel's **92 ERCC exon regions** —
real `eff_gdna`, structurally ZERO gDNA — into the `R exon` class: **2.03 %** of that class's opportunity
carrying **0** counts, which is precisely the "~2 % `R exon` deficit" that was written up as a second bug.
Scored per reference the class was never off (z **+0.63 / −0.14**).
⭐ *The rule:* a class ratio `n / (rho·E)` is meaningful only over objects that COULD supply the numerator.
Before reading one, ask what fraction of `E` sits on a population the numerator cannot reach — and if a
shipped gate stratifies, an A/B of that gate must stratify the same way or it is measuring something else.
⚠ This is `TRAPS: excluding-a-population-hides-it` run backwards: that one drops a population from the
denominator, this one adds one that cannot pay into the numerator.
*Sibling:* `never-pool-the-strata`, `a-scorer-scoped-to-the-mechanisms-targets`.

## D. Estimation and solver design

**we-keep-re-deriving-message-passing. ⭐⭐⭐ THREE OR FOUR SEPARATE SESSIONS HAVE INDEPENDENTLY RE-DERIVED
MESSAGE PASSING FROM SCRATCH, EACH TAKING MANY TURNS TO REACH THE SAME ENDPOINT, AND THE OWNER HAS HAD TO
EXPLAIN IT EVERY TIME** (the owner's own count, 2026-08-17; `DESIGN.md`'s message-passing section quotes
him). It is the most expensive recurring cost in this project, and it is not a
modelling error — nobody reasoned wrongly, everybody was slow. ⛔ **This entry exists so that arriving is
CHEAP; it deliberately does not contain the derivation.**

⭐⭐ *The tell, and a rule without a tell does not fire:* **you are reasoning about how an EXON gets its
gDNA level — from anchors, a ladder of rungs, a pooled reference, local imputation, run-fill, a throttle or
a bound — and you have not yet written the words "message passing".** Stop there. You are already inside
it, and the rest of this entry is what you were about to spend the session re-deriving.

⭐ *The endpoint, in two lines, so you can recognise having ARRIVED instead of re-deriving:* gDNA is
directly measurable ONLY where no mature transcript crosses — intergenic REGIONs, intron REGIONs, and the
BOUNDARIES against them. An EXON's unspliced mass is gDNA + RNA, which is the unknown itself, so an exon's
gDNA level **cannot be measured** and must be IMPUTED FROM ITS NEIGHBOURS along
`intron REGION ↔ intron|exon BOUNDARY ↔ exon REGION`. **Imputation across that chain IS message passing.**

⛔ *Where the settled derivation lives — ONE home each, and do not open a second:* the ruling is in
`DESIGN.md`, under the heading **"CALIBRATION IS MESSAGE PASSING, AND THE EXON IS WHY"**, and the maths in
`EQUATIONS.md` (both 2026-08-17). Find them with
**`grep -in 'message passing' docs/DESIGN.md docs/EQUATIONS.md`** rather than by section number: those two
files number by INSERTION (§3b, §6b, §9c are all inserts), so a section number quoted from here goes stale,
and this file has already been burned by exactly that.

⛔⛔ *And the reason a session re-derives rather than reads: FOUR stale beliefs, each of which reads as
"message passing was considered and rejected".* None is a reason to start over, and each is corrected in
one line:
- *"the relay was measured bad, so it is off"* — the **+154.8 %** price was measured with a NAMED,
  CONFIRMED bug live and on the **36-condition** ladder retired 2026-08-13. The bug: the composition
  licence knows about transcript TERMINI (`terminus_flank_gain`) and not about `mrna_active` flipping,
  which is the predicate that says the RNA population DIFFERS across an exon↔intron hop — so a correct
  pure-gDNA claim is relayed across a population change. It is a strict xfail in
  `tests/calibration/test_structural_reference.py`. ⛔ **Re-price, never inherit**
  (`TRAPS: re-record-the-baseline`).
- *"the mechanism does not exist yet"* — it is BUILT and switched off. `messages/relay.py`'s SPLICE IN
  (BOUNDARY → EXON) is exactly the sj-flux hop: only an EXON receives it, it carries a COUNT rather than
  an imputed belief, it has its own precision, and it is deliberately NOT tau-gated, so it survives
  unstranded data where the strand channel is dead. The switch is
  `CalibrationConfig.message_propagation = False`.
- *"a simplex vertex is unreachable, so most of the shortfall is irreducible"* — that theorem requires the
  prior to have a DENSITY, and it is correct and untouched FOR THAT FAMILY. It says nothing about a prior
  with an ATOM. `EQUATIONS.md`'s vertex section carries the scope; read it before quoting "irreducible".
- *"an exon can surely just borrow a density from somewhere"* — priced 2026-08-17 on stranded ×
  capture-ON, against `base`: a POOLED scalar reference **3.90× WORSE**, a NEAREST-RUNG one **1.27×**, a
  LOCALLY IMPUTED one **1.50×** — while that same local form is **0.4037–0.4977** (much BETTER than base)
  at capture-OFF. Hybrid-capture enrichment is PER EXON and arbitrary — it depends which probes the panel
  contains — and the anchors under-read a true exon by **2.6–3.6×**. ⛔ `capture_eff_length` cannot rescue
  it either: `transcript_capture_eff_lengths` takes a `CalibrationResult`, so it REQUIRES an
  already-solved system. ⚠ These five numbers are the RULING's, summarised here only so the shortcut reads
  as CLOSED; the anchor-ladder table beside them is their one home, and if the two ever disagree the ruling
  wins.

**one-hop-lifted-out-is-still-the-relay. ⛔⛔ LIFTING ONE HOP OUT OF `messages/` DOES NOT MAKE IT SOMETHING
OTHER THAN MESSAGE PASSING — PROPOSED THREE TIMES IN ONE SESSION** (2026-08-17). Each time the argument was
that the hop carries a MEASUREMENT — a count, an sj flux, a measured density — rather than a BELIEF, and so
"is not the relay". ⛔ *Owner's correction:* **it is message passing regardless of what it carries.** A
value computed at one object and consumed at another is a message; the payload's type is not the
definition, and "it is only one hop" is not either.

⭐ *The tell:* you are about to write a private copy of a mechanism that already exists behind a switch, and
the sentence justifying it begins *"this isn't really the relay, because…"*.

⭐ *The rule:* **rebuilding a private copy of a framework you have switched off is not a way of avoiding the
switch.** It is a second home for one mechanism, with none of the framework's licences, precisions, damping
or gates, and `TRAPS: converge-and-delete` forbids it. **The honest move is to fix the framework's KNOWN bug
and RE-PRICE it** — named in `TRAPS: we-keep-re-deriving-message-passing` above. ⚠ It is also the cheaper
move: the hop being reached for is usually already implemented (the SPLICE IN), so the private copy is a
re-implementation rather than a new capability.

**a-variance-cannot-fix-a-bias. You cannot fix a biased mode with a variance.** Established three times independently. Under
capture a counting estimate was systematically ~2× low but PRECISE — both flanking boundaries sat at the same
enriched boundary and agreed on the same biased-low density, so the bias was trusted. **A disagreement-based
variance model structurally cannot fix a bias.**

**two-gaussians-one-latent. Never hand a solver two Gaussians built from one latent.** A message on `log f` and one on
`log(1−f)` are rank-1 with correlation exactly −1, so adding their Fisher information is exactly **2×**
over-confident, rising to ~7× with deep spliced content.

**variance-fitted-on-the-belief. Never fit a variance on the current, not-yet-solved belief.** Adjacent WRONG regions agree, so the
variance collapses, the messages turn confident, and the error propagates. *Honesty measured against a
wrong belief is not honesty.* Same family: **any component trained on the solver's own output is
self-confirming** — refit iterations 1→5 went 0.0840 → 0.1056, monotonically worse.
**a-message-from-the-destinations-belief. ⭐⭐⭐ A MESSAGE COMPUTED FROM THE DESTINATION'S OWN BELIEF CARRIES ZERO INFORMATION AND CONFIRMS THE
DESTINATION.** **A message may use the destination's CONSTANTS — geometry, lengths — and its OBSERVATIONS;
never its BELIEFS.** Any "fix" that divides the destination's belief back out rebuilds the bug. ⭐ The check
is *"is the delivered value independent of the destination's own state?"*, and a prediction that does not
move when the data moves by four orders of magnitude is the tell.

⚠⚠ **THIS ONE LESSON HAS RECURRED NINE TIMES IN DIFFERENT COSTUMES.** Each sub-label below is cited from
source and tests, so the labels are kept; the investigations that produced them are in git, not here.

| | the costume it wore | the lesson that generalises |
|---|---|---|
| **TRAPS: a-message-from-the-destinations-belief** | delivered exactly `1/(1+f_own)`, reserving 33.6 % of the budget for imaginary gDNA so a zero-gDNA library read back 29.3 %; and at a gene end the rescale became `1/f(dst)`, 7× too big in the median, up to 190× | the rule above |
| **TRAPS: a-total-density-ratio** | ⛔ **a "TOTAL DENSITY" ratio is how TRAPS: a-message-from-the-destinations-belief comes back.** `r = ρ_tot(dst)/ρ_tot(src)` re-creates it whenever `ρ_tot(dst)` is dominated by a component the message is not about: a correct 0.0257 gDNA density delivered as **28.684**, and `f_g` then FLAT at 0.90 across a 10,000× RNA sweep | **a scale factor must be built from the component the claim is ABOUT, never from a total** |
| **TRAPS: substitute-the-definitions-first** | and the fix was not a better scale — substituting the definitions shows `ρ_c(src)·r ≡ φ_c(src)·ρ_tot(dst)`, a **pure composition imputation with no level in it**, so no corrective factor exists | ⭐⭐ **before correcting an operator, substitute its own definitions and read what it delivers.** Two sessions were spent hunting a better `r` for an expression that never carried a level |
| **TRAPS: the-pin-had-a-fixed-point** | the third TRAPS: a-message-from-the-destinations-belief was inside the operator built to defend against the second: the mass pin's `k = 1/(φ_msg + R_own)` has a fixed point at `R_own = ½`, so it drove the delivered gDNA FRACTION to ½ regardless of truth | ⚠ **and it hid from every aggregate** — the per-step rescales telescope back to 1 at a gene's far end, so only a per-object check away from a pure-gDNA object can see it |
| **TRAPS: no-belief-not-no-numbers** | gating the pin wherever it read the destination's own density looked clean and **broke the capture landscape** — the off-probe floor leaked into every exon | ⭐⭐ **state the licence as "no BELIEF may enter", not "the destination's numbers may not enter".** The first is TRAPS: a-message-from-the-destinations-belief and admits the structural case; the second is a superstition that costs a real mechanism |
| **TRAPS: a-licence-with-no-floor** | κ̂ = 0.500689 on a genuinely unstranded library leaves `I(f_g) ∝ (2κ−1)²` at ~5e-7, and a licence that tests a PRECISION with no floor is granted by it | a **structural** zero represented by a point estimate whose own variance covers it; the repair is to propagate `Var(κ̂)`, not to add a floor (TRAPS: a-threshold-on-a-fitted-residue refuted the floor) |
| **TRAPS: a-multiplication-gated-by-a-trace** | the same shape with a much larger trigger: an intron's **2.5 % phantom RNA** is the only nonzero RNA precision in a chain, and it alone unlocks a reframe that compounds a gDNA level to **2.16×** | ⭐ **a predicate that gates a MULTIPLICATION must be sized by how much density stands behind it, not by whether a precision is strictly positive** |
| **TRAPS: all-small-singly-large-jointly** | removing each ψ channel in turn moved the worst object by ≤0.016 of a 0.217 error while removing all of them moved it the whole way — because all three are built from the same relayed level | ⛔ **when single ablations are all small and the joint one is large, stop ablating consumers and go one stage upstream** |
| **TRAPS: recompute-from-the-oracle** | an "enrichment ratio" of 1.46 with capture OFF looked impossible and was correct: `ρ_tot` is a TOTAL and the destination held 55 RNA fragments where the source held 0. Compounded, 2.159 used vs 2.153 true | ⭐ **recompute the quantity from the ORACLE before assuming the formula is broken** — the alternative reading sends you to rewrite a correct function |
| **TRAPS: a-cancelling-defect-pair** | ⛔⛔ **fixing one of two errors that CANCEL is worse than fixing neither.** TRAPS: recompute-from-the-oracle's leak cancels across a two-hop pair; correcting one hop alone moved a toy's evidence-free exon 0.0107 → **0.0244** while the rung it was aimed at improved 26 % | ⭐ **when a fix is negative and its object sits on a multi-hop path, price it in the arm that also removes the other defect.** A pair of defects that cancel is one experiment, not two |
| **TRAPS: a-variance-cap-asserts-certainty** | ⛔⛔ **A VARIANCE CAP OF ``f(1−f)`` ASSERTS CERTAINTY AT THE CORNER, AND AN EVIDENCE-FREE SLOT'S DEFAULT BELIEF PARKS THERE.** ``node_sweep`` caps ``Var(f_g)`` at ``f_g(1−f_g)`` ("a fraction's max variance"), and the unsolved init is ``f_g = 1`` **exactly** — so the composition half of ``Var(log ρ_tot)`` is 0 at every slot with no evidence, and ``σ²_transfer`` there is the counting term alone. ⭐ Found by PERTURBING a gate for a different defect (TRAPS: perturb-every-gate): the gate still failed with ``struct_lock`` corrected, because the corner and not the mask is what zeroes it | ⭐ Replacing it with the reference prior's own variance — ``Beta(½,½)`` ⇒ ``⅛``, no tuned constant — is derivable, passes the zero control BYTE-IDENTICALLY, and is **inert (−0.3 %)**: at ``f_g = 1`` the coefficient ``[(1/E_g − 1/E_r)/B]²`` collapses to ``(1 − E_g/E_r)² ≲ 1`` too, so the term is bounded by ⅛ against ``trigamma(½) = 4.93``. ⛔ **The "coefficient diverges as ``E_g`` collapses" behaviour ``composition_logvar`` promises needs ``f_g < 1``** — any repair routed through that term is bounded by this |

**zero-the-precision-with-the-value. ⛔⛔ A REFUSED CLAIM MUST LOSE ITS PRECISION IN THE SAME STATEMENT THAT ZEROES ITS
VALUE — a value zeroed at one line and a precision added back at a later one is not "no claim", it is the
CONFIDENT ZERO the licence exists to forbid, and the relay's mass pin then converts it into its opposite.**
Found by the 2026-08-18 zero-control dissection (`g00 ss0.50 nrna_mid capture_on`, 62 k → 228 k under the
per-strand licence). The per-strand rule zeroed a refused RNA arm's value (with the sj SPLICE IN ``gp`` already
folded in) and its precisions, and the SPLICE IN block eleven lines later added the sj COUNT's precision back:
measured ``n = 0 @ pn = 280`` delivered into an empty exon. One hop on, the pin's licence read both arms as
SUPPLIED (``pg > 0``, ``pn > 0``), rescaled the budget onto the only non-zero component —
``k = M/(tg·E_g) = 467,000×`` — and delivered ``tg = M/E_g``, "all your mass is gDNA", at the gDNA arm's
unchanged precision, into a 26 k-fragment exon (true 0.000 → 0.9037). ⭐ *The rule:* an operator that
grants a precision (SPLICE IN, premise variance, conservation term) must run BEFORE the licence's zeroing, so
the zeroing is the last word — both twins, and the falsification fixture must carry the operator's input
(the first three-case gate had no sj, so its "precision zero" assertion passed on the SPLICE IN path it never
exercised). ⚠ *The same shape, from the other feeder:* an EMPTY exon whose ``struct_lock`` is the mis-scoped
``~solvable`` (TRAPS: predicate-contradicts-its-docstring) owns "RNA = 0 @ 1/trigamma(½)" on every arm, and once the
licence isolates that arm the same pin turns it into "all gDNA". ⛔ The scorer lesson beside it: a census
of the signature AT THE DELIVERED DESTINATION read 2.6 % of the error, because the damage lands one hop
past the confident zero — score the mechanism where it LANDS, not where it is written
(TRAPS: a-scorer-scoped-to-the-mechanisms-targets, second form).
*Repair landed:* the zeroing moved after the SPLICE IN in both twins (`test_gdna_scale_rule.py`, VECTORISED +
SCALAR gates, each fired by breaking its own twin). Panel: in-scope **0.931 / 0.999 / 1.000**, and
`g00` **0.753×** with all 8 zero-control rows improving, deferred 1.000 / 0.827.
*Sibling:* `a-licence-with-no-floor` (a boolean licence on a precision — here the pin's ``_may_share``).

**no-prior-means-haldane. "No prior" does not exist on a grid — omitting a term lets the grid supply Haldane**
(`p(x) ∝ 1/x`, improper, an amplifier toward the vertices). Posterior median spread over grid half-widths
4–20 is **0.045** at Jeffreys and **0.525** at Haldane.

**prefer-shares-to-differences. Sums are well conditioned; differences are not.** Subtracting across a sj gives
`Var(log ρ) = u²σ_T² + (u−1)²σ_μ²` with `u = 1/(continuing share)` — at the real median u = 2.3 already at
the boundary of validity, at p75 u = 5.3 hopeless at any depth. **Prefer shares.**

**an-all-zero-factor-is-inert. An all-zero factor is uninformative, not decisive.** In a multiplicative score, a factor that is
zero for every candidate used to annihilate the other factors and collapse the record to a coin toss.
Skip a flat-zero factor; do not multiply by it.

**density-below-one-fragment-length. Density below one fragment length is not resolvable by ANY design.** A 1 bp region has no
independently measurable density and never will (*composition* still does, since it depends on what the
fragments are, not where). *Corollaries:* an object with zero opportunity for a component must emit
nothing at zero precision, not a floored division; and "no data" must be inert, never "100 % gDNA" — that
default was actively seeding false gDNA into neighbouring exons.

**identical-paralogs-are-bimodal. An identical-paralog split is bimodal, and depth does not fix it.** Two sequence-identical
transcripts split either evenly or all-or-nothing, and which one is a coin flip on the fragment draw:
same code and seeds, `n_fragments` 3,000 → 171/0, 6,000 → 249/237, 12,000 → 679/0, 20,000 → 773/795.
⛔ **Do not "fix" this by moving a seed or a depth until it lands even** — that is tuning to green. It is
an unidentifiability, and the honest response is either to damp the degenerate direction or to score only
what is identifiable.

⚠ **THIS ENTRY SAID "the *total* is right in every case" AND THAT IS FALSE** (measured 2026-08-11). At
`gdna=100` the total is **171 against an expected ~146 — +17 %** — and the test file's own comment said so
in the same breath. Two defects live at this scenario, not one: the split is unidentifiable *and* the
total is over-called. Scoring "only what is identifiable" is right, but it must not be read as "the rest
is fine".

⭐ **The xfail here is now GREEN and the gate is stronger** (2026-08-11). `test_gdna_sweep[gdna_100]`
asserts the COLLAPSE structurally — `min == 0`, `max == total` — which still fires the day the tie is
broken, and additionally catches a PARTIAL break that a strict xfail could not distinguish from the
collapse. ⛔ The total is deliberately NOT pinned: it is a draw under an unpinned `EMConfig.seed` and moves
~14 counts across seeds, so pinning it would bake `TRAPS: the-deliverable-is-not-reproducible-by-default`
into a gate.

⭐ **And the mechanism is not what "coin flip" suggests.** With fractional assignment and only `iterations`
varied: 5 → 93.6/81.2, 12 → 147.1/27.7, 20+ → 174.5/0.0, stable to `convergence_delta=1e-9` and identical
under `mode=map`. The EM *converges* to the vertex on a real gradient. So there is a TILT entering from
outside the RNA likelihood terms, and it is a bug to localise rather than noise to tolerate.

---

**a-mean-hits-the-mass-weighted-centre-by-luck. A scalar that "wins" on a mass-weighted metric may only be
sitting on the mass-weighted centre.** Measured FOUR times in one investigation (2026-08-15/16), each time
looking like a result: the library-wide `f_lib` reference mean, the object-weighted mean, the pooled RNA
density, and the exonic class constant. ⭐ The population in each case was BIMODAL — objects at `f_g ≈ 0`
and at `f_g ≈ 1` with an empty valley — and a mean of a bimodal population is the LEAST likely value, so
the scalar carried no per-object information at all and scored well because `Σ|Δ|·M` weights by mass and
the mass sat at one mode. ⛔ **The tell: compare the candidate against ITS OWN CLASS MEAN.** If replacing
it by a constant costs nothing (measured 168,551 vs 164,074 fragments at `R exon`), it IS a constant, and
its sd against the truth's sd says so directly (0.0021 against 0.4441). ⭐ Report the sd and the
class-mean substitution beside any per-object claim.

**a-clamp-at-the-closed-end-escapes-the-window. A prior clamped at the closed end of its support can put
almost all its mass outside the solve window, and a sweep over the interior will never see it.** ψ's
reference location was clamped with a flat `eps = 1/Σeff_g = 1e-8`, i.e. `m = 1 − 1e-8`. Measured, **only
0.94 % of the reference's mass then lies inside the shipped `L = 10` window** — worse than the 57 % the
refused unequal-exponent design leaves outside, and the answer becomes a function of the grid. ⛔ The
`L`-invariance gate did not catch it because it sweeps `m ∈ [0.01, 0.99]` and the arm operated only at the
clamp. ⭐ **The repair is to derive the floor AT THE OBJECT, not over the pool**: one pseudo-fragment
*here* gives `m_i = E[g]_i/(E[g]_i + 1)`, which reads p5 0.788 / median 0.917 / p95 0.998 and keeps
0.909–0.990 of the mass inside the window. ⚠ A pooled floor (`1/Σeff_r`) does NOT fix it (0.0404 inside).
⭐ **Test a prior at the value the code will actually use, not only across its interior.**

**the-deconvolution-is-as-good-as-the-density-it-is-handed. A residual estimator inherits the error of the density
it subtracts, amplified.** `f_g = ρ_g·E_g/M` needs no RNA model — RNA is whatever the gDNA deconvolve leaves —
which is what makes it usable pre-solve. ⛔ But its accuracy is exactly the accuracy of `ρ_g`: measured
1.00× of truth at capture-OFF where the deconvolution scores **0.026**, and 0.28–0.38× under capture where the same
deconvolve scores **1.091**, i.e. worse than the uninformative reference. On the panel that became a **5.3×
LOSS** on stranded × capture-ON. ⭐ So price the DENSITY against truth before believing anything built on
it, and state the deconvolution's verdict per stratum of that density rather than pooled.

**deriving-one-coordinate-propagates-its-error. ⛔⛔⛔ CLOSING A COMPOSITION BY DERIVING ONE COORDINATE
FROM ANOTHER MAKES THE DERIVED ONE INHERIT EVERY DEFECT OF THE SOURCE — INCLUDING ONES THE OLD ESTIMATOR
DID NOT HAVE.** ψ's composition was closed by defining the RNA total as ``1 − f_g`` instead of reading it
independently. Correct, and it removed a real defect — but ``f_g`` was the GRID-SNAPPED median, and the
retired read-out had been the posterior MEAN, which on an evidence-free object is **exactly ½** because the
reference is symmetric. So the fix silently traded an exact number for a snapped one: `test_relay_mass_rescale`
read ``R_own = 0.51256`` where it had been ``0.5``. ⭐ *The rule:* when you make coordinate B a function of
coordinate A, enumerate what B used to be exact about and check A is exact about it too — the merge is only
sound if A is at least as good, and "at least as good on average" is not the test. ⚠ The repair was not to
undo the derivation but to fix the SOURCE (a continuous quantile), after which both were exact. Related:
`TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`.

**interpolate-on-the-axis-where-the-lattice-is-uniform. ⛔⛔ A QUANTILE INTERPOLATED ON A NON-UNIFORM
TRANSFORM OF THE GRID IS BIASED TOWARD THE MIDDLE, AND IT LOOKS RIGHT.** ψ's ½-quantile was interpolated in
``f_g``-space, where the σ lattice spacing runs ~1e-5 at the ends and ~0.085 in the middle. A bin's midpoint
in ``f`` is NOT the image of its midpoint in ``λ``, so a posterior concentrated on one grid point came back
**biased toward ½ by 2.71e-03** at ``n_grid`` 60 (1.48e-04 at 256) on every interior point. Interpolating on
``λ`` — where the lattice IS uniform — returns that posterior's own grid point to **2.2e-16**, and σ being
monotone makes the mapped quantile exact. ⭐ *The rule:* do quantile arithmetic on the axis the grid was
BUILT on, then map; never on the transformed axis. ⚠ It also restores the property the estimator exists for:
``|median(1−f) − (1 − median(f))|`` went **8.45e-02 → 3.3e-15**. ⛔ And the case is not synthetic — an
unsolved slot's fed-back belief produces a one-hot posterior, landing exactly on the bias.

**read-the-whole-failure-list. ⛔⛔ 21 FAILURES WERE REPORTED AS "21 GOLDEN FAILURES" FROM THE LAST ELEVEN
LINES OF THE OUTPUT. TWO WERE REAL.** pytest prints the tail; a long ``short test summary info`` scrolls the
rest away, and reading the visible names and generalising is how a genuine regression is filed as expected
churn. One was a broken caller (a changed return arity), the other was the composition defect above — the
single most informative failure of the session, nearly dismissed. ⭐ *The rule:* derive the failure set,
never eyeball it — ``pytest -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c`` prints one line per
FILE and fits on a screen at any failure count. ⚠ Same shape as `TRAPS: re-record-the-baseline` one level
up: the number was read faithfully and the METHOD of reading it was never checked.

**a-priors-curvature-is-not-the-datas-information. ⛔⛔⛔ A PRIOR MAY NOT CONTRIBUTE TO A FISHER PRECISION,
AND THREE SEPARATE MEASUREMENTS SAY SO.** `region_init`'s `tau_lam` is the DATA's information on the
composition axis. When ψ's reference gained an annotation-set MEAN, a slot it pinned sat at the λ grid's
edge where the strand term `c·a² ∝ f_g²(1−f_g)²` vanishes, and `tau_lam` fell **3,227×** (0.1358 →
4.21e-05) exactly as the slot's belief became most certain. The reading — *"the information moved into a
λ-factor whose curvature nothing reads; `density_factor_precision` already reads the intron factory's, so
let it read the location's too"* — was built and REFUSED:

* **the fall is the JACOBIAN, not a loss.** `[f(1−f)]²` between `f_g = 0.98576` and `0.99975` predicts
  **3,154×** on its own — ~98 % of it. The likelihood genuinely is that flat on λ at the point the prior
  chose. There was nothing to recover.
* **it is a BOOLEAN GATE FLIP, not a contribution.** At the vertex `Var(log f_g) = (1−f_g)²/τ ≈ 8e-08`, so
  `own_precision` saturates at the COUNT ceiling: τ = 0.029 and τ = 1e6 both return **850.44** against a
  ceiling of **850.50**. Only `τ > 0` does any work, and what it releases is 850 fragments of count
  precision — bought with a prior worth one.
* **it credits DATA-FREE slots.** The location carries no count, so it hands the same value to an empty
  slot and to a 10⁶-fragment one: `n = 0` goes `prec_g` 0 → 0.2026. `intron_prior` is NOT the precedent it
  looks like — its NegBinom curvature is count-derived and self-limits on thin data.

⭐ *The rule:* before adding any term to a precision, ask whether it scales with the DATA. If it does not,
it is a prior, and its place is in the posterior — never in the Fisher information. ⚠ Measured, the
contribution was **bit-identical on the deliverable across all 32 panel rows** and moved only
`solvable_mass_share` (13.16 → 11.91) and `weak_evidence_mass_share` (0.338 → 1.549) — i.e. it changed
0.8.0's own DENOMINATOR (`has_own_composition_evidence`) and nothing else. Related:
`TRAPS: a-four-decimal-print-is-not-a-zero`, which is how the wrong diagnosis was reached.

**a-refutability-test-needs-the-refuting-channel-in-the-fixture. ⛔⛔⛔ A PRIOR MEASURED WITHOUT THE
EVIDENCE THAT COULD OVERTURN IT IS NOT BEING TESTED — IT IS BEING ASKED TO ANSWER ALONE.** ψ's structural
prior (*presume gDNA where no annotated mature transcript crosses*) was measured against nascent RNA on a
toy chain and read **2.1–5.7× WORSE**, and the conclusion drawn was that the claim itself was unsafe. The
chain was `exon|intron|exon|intron|exon` with **no intergenic regions** — and `fit_intron_background` pools
INTERGENIC regions only, so it returned uninformative, `_build_intron_prior` returned `None`, and the
intron-vs-intergenic density mechanism **was not in the room**. At κ = ½ the strand mechanism is dead by
derivation too, so the fixture had NO refutation channel at all. Rebuilt with intergenic flanks — which
production always has — the same mechanism reads **τ_fac = 161.4** at every intron and the same prior
YIELDS to the same nascent, to within 0.02 of the no-prior answer.
⭐ *The rule:* before measuring whether a prior can be overturned, ASSERT that the overturning channel is
live in the fixture — `bg.informative`, a nonzero `τ`, a fired counter. ⚠ It is the same shape as
`TRAPS: an-ablation-that-never-ran` pointed the other way: there an override never ran, here the EVIDENCE
never ran, and both make the thing under test look like it has no competition.

**a-strength-is-a-nat-a-prior-weight-is-a-count. ⛔⛔ THE TWO ARE NOT THE SAME UNIT, AND EQUATING THEM IS AN
ANALOGY WEARING A DERIVATION'S CLOTHES.** ψ's location term has range `logit(m)` — a LOG-ODDS in nats — and
the reference it rides on is `Beta(a, b)` with a declared weight `a + b = 1` PSEUDO-COUNT. Setting
`m = σ(a+b)` reads as principled, gives 0.731, and is numerically fine; it is still dimensionally wrong.
⭐ The derivation that stays in one currency is *what mean would one pseudo-observation produce*:
`Beta(a,b) → Beta(a+1,b)`, mean `(a+1)/(a+b+1)` = **0.75** at `a = b = ½` — a composition mean from
pseudo-counts, tracking the exponents automatically.
⛔⛔ **And what it replaced is the sharper lesson: the strength had been taken from the LATTICE** (`m = σ(L)`
⇒ 9.31 nats ≈ 10,000:1), which is not a statement about the claim at all. Measured against being refuted it
was **worse than having no prior** (2.0247 vs 0.3946). ⭐ *The rule:* a prior's strength must come from what
it CLAIMS, never from what the representation can hold — the representation's limit is a CAP, and a cap is
not a choice. ⚠ Check the optimum is BROAD before trusting any derived value: here 0.69 → 1.50 nats lie
within 6 %, so the derivation is doing the work rather than a fit.

**a-four-decimal-print-is-not-a-zero. ⛔⛔ A DIAGNOSTIC PRINTED AT FOUR DECIMALS READ `0.0000`, AND A
SESSION'S WHOLE MECHANISM WAS BUILT ON IT BEING ZERO.** `tau_lam` under a sharp prior was recorded as
collapsing to **0.0000**, and the diagnosis followed: *"`tau_lam > eps` is `has_own_composition_evidence`,
so the strongest possible prior makes a slot read as evidence-FREE"*. The value was **4.21e-05**, the
predicate returns **True** on it, and the slot was never evidence-free. The real effect was a 3,227×
change of MAGNITUDE, which needs a different repair and — as it turned out — no repair at all.
⭐ *The rule:* a number that a diagnosis TURNS ON is `repr`'d, never formatted, and a claim of "exactly
zero" is asserted with `== 0.0` rather than read off a table. ⚠ The cost was a falsification gate written
against a non-defect: it demanded that a prior never LOWER the data's Fisher information at the posterior
mode, which is false for any nonlinear reparametrisation. It was verified failing, the code was changed to
make it pass, and both were wrong.

**a-constant-in-exact-arithmetic-is-not-constant-in-float64. ⛔⛔ "THIS TERM IS CONSTANT SO IT CANCELS" IS
A CLAIM ABOUT ℝ, AND THE GRID DOES NOT LIVE THERE.** ψ's location term at its NEUTRAL mean `m = ½` is the
constant `log 2` in exact arithmetic — its whole reduction property, and the reason setting a location on
some slots is supposed to leave every other slot alone. In float64 `logaddexp` leaves that row with
**`ptp = 2.22e-16` over three distinct values, at every grid size**. `_posterior_median_fg` reads the grid
point where the CDF first reaches ½, so wherever a posterior is exactly SYMMETRIC that rounding tips the
knife-edge: a balanced AMBIG slot at `κ = ½` moved `f_g` **0.5423 → 0.4577**, one full grid step, at a slot
the prior claimed to say nothing about — and there are ~37,000 of them, the exact population the gDNA
landscape is fitted to serve. ⭐ *The rule:* where a term is documented as vanishing, RETURN the vanishing
value rather than computing something that equals it; and gate the reduction on the ANSWER
(`solve(neutral) == solve(None)`, bit-for-bit), not on `allclose` of the term. ⚠ An `atol = 1e-15`
closeness gate on the term itself passed throughout.

## E. Structure, indexes and plumbing

**one-reference-hides-refid-bugs. Single-reference synthetic indexes hide reference-id-space mismatches.** A resolver assigned ref-ids
by first-seen interval order rather than the index's own order, so **476,719 of 476,732** real fragments
were silently dropped inside `deposit()` while every golden test passed. *And the fix is per-path:* extra
references make the space non-trivial for the paths those references exercise. RNA-only spike-ins do
nothing for the gDNA deposit path, which is a different branch — that needs **≥ 2 genomic references**.

**annotated-is-not-genomic. "HAS AN ANNOTATION" IS NOT "IS GENOMIC" — a classification must be an INPUT, not a proxy.** The
simulator chose its gDNA references with `{t.ref for t in transcripts}`. Every RNA-only spike-in carries
exactly one transcript, so every one qualified, and the panel filled with gDNA on synthetic RNA
references — on the very templates whose truth abundance of zero is the false-positive control. **And a
mis-stated classification must raise, not silently produce nothing.**

**an-sj-is-not-a-gap. A splice junction CANNOT be detected from a gap between deposited slices.** A contiguous spliced
read whose exon body straddles an internal region bound has no gap at all. The sj's identity is the cut-intron
coordinates, which the scanner already has — pass them through.

**deposit-at-the-sj. A splice deposit belongs at the SJ's coordinate, not the region's boundary.** Invisible for
annotated introns (their ends are region bounds); for unannotated sj the mass lands kilobases away.

**splicing-makes-the-graph-cyclic. Alternative splicing makes a region↔sj graph CYCLIC, and cycles are the common case** — a
cassette exon is a 4-cycle, and the human graph has ~404,000 independent cycles, one per sj.
Two-sweep forward-backward is exact only on a tree. **Never break a cycle by dropping a sj boundary** —
that re-isolates the exon the boundary exists for.

**nrna-does-not-mean-synthetic. `~is_synthetic & ~is_nrna` as the "real transcript" filter deleted 52,104 real termini.** On a
non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent", NOT "manufactured span". **ONE
filter: `~is_synthetic`.**

**credit-exactly-one-sj. A fragment crossing K sj must credit exactly ONE** (the leftmost annotated). Crediting all K
shifts the library sense fraction 21–34 % and creates between-side correlation that reads as
overdispersion (a zero-overdispersion simulator fit 0.092). A `1/K` split is provably biased by 4–12 σ.

**strand-completes-the-sj-key. `(src, kind, dst)` is not a total order for sj boundaries** — two strand-coincident sj differ
only in strand, so ordering becomes input-order-dependent and the duplicate check reads them as duplicates.
GENCODE has zero of them, so only a synthetic stress test finds it. Sort on `(src, kind, dst, strand)`.

**a-hash-that-misses-its-artifact. Cache keys that do not cover the artifact they cache.** A partition hash covered only the region file;
a flag fix rewrote every boundary file while leaving every region file byte-identical, so a stale cache would
verify CLEAN. **Never store a derived hash beside the data it describes; compute it on demand.**

⚠ **A SECOND FORM, and it is the one you build deliberately: AN EXCLUSION FROM A KEY IS A CLAIM ABOUT
WHAT IS STORED.** Measured 2026-08-08 on a study cache for the `flgap` pair, whose docstring said
`priors.py` was *deliberately* outside the key "so that editing the assembler does not invalidate a
5-minute scan" — a genuinely good decision that is sound only while nothing `priors.py` produces is in
the blob. Two things were: the assembled `p_arm` and the projected `f_gdna`. So an assembler edit served
a **fresh O beside a stale P and F**, and the comparison between them — the whole point of the cache —
was meaningless, with the key still reading `ok`. ⚠ That cache and its one consumer were deleted on
2026-08-17 with the panels they read, so this paragraph is the lesson's only home; the rule is about
CACHE KEYS, not about that file.

⭐ **The repair is not a better key, it is to store the INPUTS and derive at read time.** The blob was
changed to hold `cal`, `multi_loci`, the truth masses/shares and the raw per-region start counts, with
P, O, S and F all assembled on load, so they were always the same assembler. ⛔ And the two files that
*produced* the stored truth (`_oracle.py`, `prior_vs_oracle.py`) were missing from the key entirely, as
was the builder itself — a build that starts storing a new artifact must invalidate the blobs that lack
it.

⛔ **The tell: read the key's exclusion list as a sentence.** "This file is excluded *because* nothing it
makes is stored here" is checkable in one grep; "excluded because the loop is faster" is not a reason,
it is the cost of one. ⚠ And a stored artifact that is a PURE FUNCTION of stored inputs should carry a
build-time byte-identity gate proving the read-time derivation reproduces the production call, otherwise
dropping it from the blob is an unproven claim.

⚠ **A THIRD FORM, IN A DATACLASS, WITH NO CACHE ANYWHERE: A STORED FIELD SURVIVES
`dataclasses.replace`.** An oracle arm builds itself by swapping a `CalibrationResult`'s mass and count
arrays for truth with `dataclasses.replace` (`prior_vs_oracle.OVERRIDE_FIELDS`). Any field *computed at
construction* from those arrays comes through the swap intact and then silently describes **the arrays it
replaced** — the same staleness as a missed hash key, with no key to inspect and nothing to invalidate.
⭐ *The rule:* **if a value's inputs are in somebody's override list, that value may not be a stored
field.** `sj_conserved_mass` and `library_rna_fragments` are `@property` for exactly this reason, so the
oracle arm's number is the oracle's by construction — the same repair as the second form (store the
inputs, derive at read time) applied to an object rather than a blob.

**integer-channels-reproduce. Integer channels are bit-identical across worker counts; float channels are not.** Max relative
3.7e-7 per cell, which propagated to a ~2.6 % difference in the calibration output. Integer addition is
associative — that is the whole fix. *Corollary:* the one bank whose ORDER is observable (a list, not a
sum) must be sorted on its own content before it crosses the ABI.

⛔⛔ **AND THE PRICE TAG ABOVE IS A float32 NUMBER. IT WAS QUOTED AGAINST float64 FOR MONTHS AND IT DOES
NOT CARRY** (measured 2026-08-11). ``3.7e-7 relative per cell`` is ~3 float32 ulp; float64 is ~1e-16 —
**3.4e9× finer** — reaching the deliverable at ~1e-11, five orders below ``EMConfig.convergence_delta``.
The measured run-to-run spread of the shipped multi-threaded pipeline is **1.503e-15** on
``posterior_mean`` and **exactly 0** on every integer-derived column.

⭐⭐ **The bigger correction: the fixed point was LESS ACCURATE than the float it was defending against.**
Against exact rational arithmetic on the reciprocal-opportunity theorem — every admissible placement
deposits ``1/A``, so each length contributes exactly one density unit — the two representations miss the
integer answer by::

    node_len 151    fixed 7.0e-10    float64 5.8e-15      120,000x better
    node_len 400    fixed 1.7e-08    float64 1.0e-13      170,000x better
    node_len 1000   fixed 2.0e-07    float64 2.8e-13      714,000x better

⚠ And the "exactness" two gates asserted was a property of their FIXTURES: ``1/2 + 1/3 + 1/6`` lands back
on ``2^32`` because two rounding errors cancel, while ``1/3 + 1/3 + 1/3`` is one quantum short — and
float64 is exact on both. Both gates asserted the fixed point's OWN closed form, so they re-derived the
implementation and could not fail (`TRAPS: a-gate-that-reconstructs`).

⭐ **The rule's FIRST sentence is still true and is why the split survives**: a COUNT is an integer and
reproduces exactly; a FRACTION is float64 and does not. Owner ruling 2026-08-11: one numeric convention,
the tool is not bit-reproducible, and tests validate within a DERIVED tolerance — bracketed from both
sides, so a nudge just past it fails and one just inside it passes.

**worktrees-run-the-wrong-code. Worktrees silently run the wrong code** — an editable install's meta-path finder beats
`PYTHONPATH`, so an A/B inside a git worktree executes the main repo's source.

**checkout-deletes-uncommitted-work. `git checkout -- <file>` does not undo a perturbation when the work is uncommitted — it deletes the
work.** A perturbation harness must restore from a copy of the WORKING TREE. Cost one full
re-implementation.

**two-masks-one-name. TWO DIFFERENT MASKS SHARED THE WORD `struct_lock`, AND BOTH WERE RIGHT.** One answers "is this
belief pinned and certain?" (both axes); the other "may this slot EMIT composition certainty into its
messages?" (NODE-only on purpose, because a boundary is structurally gDNA but sits between RNA-carrying exons).
⭐ Two correct predicates under one name is worse than either being wrong: give each its own name and ONE
home, and let every consumer import it (TRAPS: a-test-that-redefines).

**two-docstrings-one-quantity. The prose next to the code said "the AVERAGE" and the code followed the prose,** while a sibling
module's docstring had the correct formula the whole time. Two docstrings disagreed about one quantity for
months and nobody diffed them.

---

**a-transcript-predicate-must-not-silently-drop-a-molecule. ⛔⛔⛔ A FRAGMENT REJECTED BY A TRANSCRIPT-LEVEL
PREDICATE DEPOSITS NOTHING INTO CALIBRATION — AND THE REJECTED POPULATION IS NEVER RANDOM, SO THE BIAS IS
CONCENTRATED WHERE THE PREDICATE FIRES.** Found 2026-08-19 as the cause of the 2–3 % boundary-crossing
deficit that failed FIELD certification on 8 of 32 conditions. `detect_chimera` joined a fragment's two
mates iff their TRANSCRIPT SETS intersect, so a contiguous genomic molecule spanning two transcripts that
share nothing was called `CHIMERA_CIS_STRAND_SAME` and dropped. Measured: **4,087 gDNA fragments per
condition — 0.04 % of fragments carrying 2.4 % of every boundary crossing**, because the predicate fires
exactly where transcripts are short and dense, which is exactly where a fragment crosses many boundaries
(3.66 crossings lost each against 1.65 for an average crosser).
⛔ **It hid because every dropped fragment was a CROSSER**, so `region_contained` matched the BAM's own truth
EXACTLY on all 22,139 regions and the conserved-mass gate closed. Only the crossing axis moved.
⭐ *The rule, owner-ruled:* **check genomic COMPATIBILITY before calling a fragment a chimera** — one
reference, mates facing inward, implied fragment length within `max_frag_length`. gDNA is contiguous and
routinely spans unrelated transcripts; that is what an annotation looks like, not a rearrangement, and such
fragments carry no junction to infer one from.
⭐⭐ *The general rule, and it is the transferable part:* **make the fragment ledger CLOSE.**
`n_fragments == deposited + Σ accounted drops`. Before the fix the panel read `10,000,000 − 9,988,687 =
11,313 = 7,226 deferred + 4,087 chimeric`, with the chimeric ones counted in NO qc field and absent from the
census; after it, `n_chimeric = 0` and the residual is the deferred count exactly. ⛔ A drop with no counter
is invisible to every downstream check.
⚠ **The same shape is still live elsewhere and is an owner decision**: a CIGAR sj rejected by the artifact
blacklist promotes the fragment to `SPLICE_ARTIFACT`, which also deposits nothing. It is inert on the
simulated panel (no blacklist artefact in the index, `sj_blacklist_size = 0`) and would fire only on real
data — which is how it has stayed unmeasured.

**an-object-class-does-not-see-a-terminus. ⛔⛔ THE SEVEN OBJECT STRATA CLASSIFY A BOUNDARY BY THE EXON-NESS OF
ITS TWO FLANKS, AND A TSS/TES LYING INSIDE ANOTHER TRANSCRIPT'S INTRON OR EXON IS INDISTINGUISHABLE FROM A
SPLICE SITE THERE.** Measured 2026-08-19 (`hop_currency.py`, Stage 2 of the message-propagation rebuild).
Classified by object class alone, `R exon <- B exon|intron` read COMPOSITION 221 k fragments against LEVEL 27 k
on `g50 ss0.50 nrna_mid capture_off` — the SPLICE-IN population fix had apparently not rescued the
composition into exons. Split by the splice graph's flags (`is_terminus` / `is_splice_site`, either
strand), **209 k of it sat on the 3,867 TERMINUS boundaries** (COMP 209 k vs LEVEL 12 k on 327 k) and the
15,480 SPLICE-SITE boundaries were near-exact both ways (11.7 k vs 15.0 k on 605 k). A `B exon|exon` splits
the same way: 4,838 `[sj]`, 7,850 `[term]`, 123 `[sj+term]`. The two sub-classes have OPPOSITE currencies —
at a terminus the RNA ORIGINATES and a composition cannot cross; at a splice site it ENTERS by the sj and
the composition does — so a table keyed on the class averaged a 64 % failure with a 2 % success and
named neither.
⭐ *The rule:* a hop type is `object class x {sj, term, sj+term}`, read off `RegionStatics.boundary_flags`; any
static hop table a policy builds in `prepare()` must key on it. ⚠ The classes themselves are right for what
they were built for (ψ's reference, `object_composition.py`) — a terminus boundary IS `exon|intron` by
flank; it is the HOP that needs the second bit.
*Sibling:* `two-masks-one-name`, `a-boundary-with-rna-is-not-an-sj`.


## F. Domain facts that read like defects

**specificity-and-sense-are-complements. Strand specificity is TWO different quantities and they are complements.** A simulator's
`strand_specificity` is protocol **fidelity** (direction-agnostic — an R1↔R2 swap with probability
`1 − ss`); a fitted `sense fraction` is **directional**. For an R1-antisense (dUTP) protocol — which real
cfRNA is — comparing them reads as a sign error. A fitted κ of 0.0101 on a "0.99 stranded" library is
**correct**, and forcing 0.99 measured 166× worse. The matching quantity is
`StrandModel.strand_specificity`, which recovers the knob directly.

⚠ **The recovered values are NOT restated here.** They lived in this paragraph *and* in
`test_strand_specificity_RECOVERS_the_simulated_parameter`, and both drifted: each said
`0.75 → 0.7701, 0.50 → 0.5020` while the tree measured **0.7494 / 0.5034**. A measured number belongs in
the test that measures it, where it cannot rot silently — one home, not two.

⭐ **And the protocol DIRECTION is now simulable in both orientations** (2026-08-11): `ReadSimConfig.r1_sense`
selects R1-antisense (TruSeq dUTP, the default) or R1-sense (KAPA), independently of the fidelity. That
closes the coverage gap this entry describes, and the per-fragment mirror between the two is gated.

**strand-measures-the-tilt. Strand measures the TILT, not the gDNA fraction.** With RNA tilt `d = f₊ − f₋`,
`p = ½ + (κ−½)·d` — the gDNA fraction **cancels identically**. Strand reaches gDNA only through the
triangle bound `f_g ≤ 1 − |d|`: tight on a single-strand region, slack on a both-strand region. And
`I(f_g) = 0` **exactly** at κ = ½, for any count and any overdispersion.

**a-linear-likelihood-emits-a-sign. ⛔⛔ A LIKELIHOOD THAT IS ASYMPTOTICALLY LINEAR IN ITS PARAMETER HAS
NO MODE, ONLY A DIRECTION — AND ON A BOUNDED GRID THAT DIRECTION SATURATES AT AN ENDPOINT.** Measured
2026-08-10 on the drained arm; it cost a feature, and the feature was purged at `b7ed7a0b` (`f470a570`
is the measurement checkpoint before it, which changed no `src/` and is where the instruments survive).

A fragment-length composition channel was built on a bivariate Gaussian in `(Σ1/w, Σw)` whose mean moves
with the composition `π` by `N·Δ`, `Δ = μ_g − μ_r`. The term quadratic in `π` is `O(Δ²)` while the linear
term is `O(Δ)`, so as the two fitted pmfs converge the row degenerates into a straight tilt whose maximum
can only be a grid endpoint. What it then emits is `sign(Δ' V⁻¹ (x − N·μ))` — the **direction of the
model's residual, with the magnitude divided out** — and when `Δ ≈ 0` that residual is mis-specification,
not composition.

⭐ **The measurement that proves the answer is not a function of `Δ` at all.** Interpolating `gdna_pmf`
toward `rna_pmf` across TWELVE orders of magnitude: at a gap of ~1e-9 bp the channel reported **0.72 /
0.59 / 0.72** on libraries whose truths were **0.00 / 0.00 / 0.57**, and growing the gap to its real value
made every one of them BETTER. ⛔ So closing the gap is strictly harmful and every "fix the pmf" repair is
refused. The fragment-length MODELS were exonerated by the same run.

⭐ **The tell, available before any measurement: expand the log-likelihood in the parameter and look at
the order of the leading term.** If it is linear, the argmax is a sign and the object is a vote, not an
estimate. The sufficient summary of such a channel is `(π̂, I)` — a location and an information — never a
row handed to a normaliser.
⛔ **DEFERRED POST-0.8.0** (owner, 2026-08-14): the channel was purged at `b7ed7a0b`, does not exist in
`src/`, and is not a 0.8.0 candidate — do not propose it, rank it or build it. ⭐ The lesson about a
likelihood that is linear in its parameter is general and applies to any channel.

**amplitude-fades-influence-does-not. ⛔ A TERM THAT NORMALISES AWAY CAN STILL DECIDE THE ANSWER, BECAUSE
AN ARGMAX IS SCALE-FREE.** Same campaign, and it is why fixing the precision was never going to be enough.

The channel's row amplitude faded correctly — linearly in `Δ`, `9.69e-12` over a `1e-12` gap ratio. Its
**participation** did not: the discrimination guard is an EXACT inequality, so `Δ = 0` gave 0 live slots
and `Δ = 5e-12 bp` gave **100.00 %** of the library mass. Its **declared precision** did not either:
`density_factor_precision` normalises the row and inverts the variance, so a near-flat row returns
`1/Var(uniform over the grid)` — the grid's own width sold as evidence, on 100 % of live slots.

⭐ **Three quantities, three different behaviours, and only the first one was ever checked.** When a term
enters a normalised sum, report its amplitude, its participation and its declared precision SEPARATELY;
they do not fade together, and it is the two that do not fade that decide the answer.
⛔ **DEFERRED POST-0.8.0** (owner, 2026-08-14): same channel, same ruling as
`TRAPS: a-linear-likelihood-emits-a-sign` — not a 0.8.0 candidate, so do not propose it. ⭐ The
three-quantities rule applies to ANY term entering a normalised sum and is why the entry is kept.

**a-pooled-conversion-applied-per-component. ⛔ A RATIO MEASURED ON THE POOLED POPULATION AND APPLIED TO
EACH COMPONENT SEPARATELY IS POPULATION-BLIND — AND THE BLINDNESS NEED NOT BE THE AXIS YOU EXPECT.**
Measured 2026-08-10, `boundary_q_population.py`.

`assemble_priors` converts a boundary's crossing INCIDENCE to a FRAGMENT count with `q = mass/count`, measured
per boundary on the whole population and applied to the gDNA and RNA parts separately. `q` is an explicit
function of fragment length, so the obvious worry was that a length gap breaks it.

⛔ **It does not. The equal-length null shows the LARGEST error.** At a 4.98 bp gap `q_g` = 0.6330 against
`q_r` = 0.5233 and the boundary term is over-counted by **+3.18 %** — five times the +0.63 % at a 17× larger
gap, and the sign does not track the gap either. ⭐ The driver is **PLACEMENT**: gDNA is genomically
uniform and crosses boundaries in long intergenic regions where `q → 1`; RNA is confined to transcripts, where
exon regions are short. The two populations occupy different parts of the partition, so they differ at
identical lengths.

⭐ **Bounded and left alone deliberately**: on the TOTAL prior the region term dilutes it to `Δphi`
**+0.00013 … +0.00596** — at most 0.6 pp, at or below calibration's own noise floor. ⚠ And it is not
repairable in production: the driver is *where* each population sits, so `q_c` cannot be modelled from the
fitted pmfs, and the pooled bank is the only per-line evidence there is. **Record the bound, build
nothing.**

**an-identity-with-a-qualifier. ⛔⛔ AN IDENTITY THAT HOLDS "OVER X" IS A MEASUREMENT WAITING TO BE MADE.
PRICE THE COMPLEMENT, OR THE QUALIFIER IS A HOLE NOBODY HAS SIZED.** Measured 2026-08-10; it was worth
25.3 % of the RNA library and had been documented, accurately, for months.

The conserved-mass rule's own specification said its identity was exact *"over deposited, **unspliced**,
annotated fragments"*. Every word was true and the scoping was deliberate. But nobody had ever asked how
big the complement was — and it was **1,222,375 of 4,830,713 RNA fragments** on one ordinary panel
condition. A spliced fragment whose every block lies inside one region crosses no boundary and is not
`contained` either, so it deposited on **no conserved bank at all**: it existed on the incidence axis and
on none of the conserved ones, which is why a library fragment count was not computable.

⭐ **The complement was cheap to price and the pricing is what made the fix obvious.** One census against
per-fragment truth: gDNA `1.000x` deposited, RNA `0.747x`. The asymmetry names its own cause — gDNA cannot
splice — and the repair followed in one step.

⛔ **The second half of the trap: a docstring had already promised the missing accounting existed.**
`boundary_spliced_mass` said a block with no interior boundary has *"their accounting on the sj axis"* —
and there was no `sj_mass`, so nothing kept the promise. **A cross-reference to a bank is not evidence the
bank exists.** Same family as TRAPS: two-docstrings-one-quantity: prose that describes an intended design reads
as prose that describes the code.

⭐ **The tell, and it is free:** grep your invariants for "over", "for", "assuming", "where". Each one
names a population; each population has a size; and until someone measures that size the invariant's
scope is an assumption wearing the clothes of a theorem.

**equal-lengths-carry-no-composition. At equal component mean lengths the length channel carries EXACTLY ZERO information about
composition, at any depth.** The 2×2 deconvolution is identified only through `μ_g − μ_r`. A claim that
one storage choice beats another was measured at a 4× mean separation and is **false** at every region
≥ 250 bp and reversed at equal means.
⭐⭐ **This fact is the gDNA ladder's FORCING FUNCTION, not merely a caveat about it** — zero length
information is exactly what stops the EM answering WITHOUT calibration
(`TRAPS: a-length-gap-bypasses-calibration`). ⛔ **DEFERRED POST-0.8.0** (owner, 2026-08-14): the
composition channel built on the non-equal case is not a 0.8.0 candidate — this entry states a domain fact,
never a candidate to build.

**capture-is-1000x-on-exons. Hybrid capture is ~1000× on exons and only gDNA reads it cleanly.** RNA's own 10⁴ expression range
hides the probe pattern; gDNA's uniform baseline does not. And **capture destroys the intron signal 75×**,
so nascent-vs-gDNA is fundamentally unidentifiable under capture.

**capture-selects-for-length. Capture SELECTS FOR LENGTH, and the post-capture distributions are the baseline.** Probes hybridise
to sequence, so a short fragment presents less sequence and is captured less efficiently. The pre-capture
parameters describe a library that was never sequenced — score against post-capture truth, for fragment
lengths exactly as for abundances. *And capture narrows the gDNA↔RNA length gap* whenever gDNA is the
shorter component, because the short tail it removes is disproportionately gDNA.

**on-target-by-start-is-geometry. "On-target" defined by the START's territory is geometry, not capture efficiency.** Conditioned on
being captured, a fragment whose start is in the intron is one that was **long enough to reach the
probe** — weight ~`w²/2` — while an exonic start carries ~`p²/2`, flat in `w`. So an intronic-start
population reads *longer* than an exonic-start one under any capture model. The population that
physically binds is the one that **overlaps a probe**.

**real-data-is-a-test-input. Real data is a TEST input, never a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
spectrum, not a sample of it. Sweep the plausible space, report the worst case, bring the domain call to
the owner. In particular: **never assume RNA fragments are longer than gDNA** — true for cfRNA, false
elsewhere.

**eff-lengths-do-not-cancel-at-an-end. "Effective lengths cancel, so a region's marginal is just its length" is FALSE near any transcript
end** — a mature fragment must fit in the remaining transcript; gDNA need not.

**configured-lengths-are-not-realised. EQUAL CONFIGURED FRAGMENT LENGTHS DO NOT GIVE EQUAL REALISED ONES.** Handing the simulator one
mean for gDNA and RNA still yields different realised means, because each pool's length marginal is
reweighted by its own template opportunity — a 2 kb transcript truncates the tail that a whole chromosome
does not. ⭐ So "equal lengths" is a claim about the CONFIG and never about the library; gate the axis on the
realised truth (TRAPS: prove-the-substrate).

**mature-rna-never-crosses-a-boundary. Mature RNA never crosses an exon↔intron boundary** (0 of 1,146 boundaries, 7/7 conditions). Exon↔exon boundaries
do. This is the hard empirical case that a contiguous boundary and a splice junction are physically different
objects — and it is what makes the two exon-crossing gDNA pools pure.

**a-boundary-with-rna-is-not-an-sj. But a boundary with RNA crossing it need not be a sj.** One position can be a splice donor for
transcript A and plain contiguous exon for transcript B; zero-gDNA libraries show boundaries with 44–55 k
unspliced fragments that are 100 % RNA.

**the-panel-enriches-nascent-by-its-own-probes. ⛔⛔ A SIMULATOR THAT MODELS A POPULATION IN ITS OWN
PRIVATE SPACE WILL GIVE IT ITS OWN PRIVATE PHYSICS — and the tool is then debugged against a library that
does not exist.** ⭐ **FOUND 2026-08-19 (`hop_currency.py`), REPAIRED THE SAME DAY.** Nascent RNA was a
parallel per-transcript "pre-mRNA space" beside the transcriptome, and hybrid capture was applied inside
each space separately: `sim/capture/sampler.py` gave a pre-mRNA only the probes declared on ITS OWN
transcript (`_add_transcript_probe` → `_add_nrna_blocks(t_idx, …)`; `_add_bed12_probe` skipped a
transcript the probe did not project onto, and filtered by strand). So nascent RNA spanning ANOTHER
transcript's probed exon — a nested or overlapping gene, an unprobed sibling isoform (360 of the 4,154
transcripts in the ladder's probed genes) — was NOT enriched while the gDNA at the same position WAS:
observed nascent at those boundaries **0.16x** of "enriched alike", and an 18 % composition residual on
the `nrna_mid x capture-ON` intron→boundary hops that is not the estimator's.
⭐ *The physics it violated:* at capture time every library molecule is ds-cDNA, so a probe binds anything
carrying its sequence — mature, nascent and gDNA alike, either strand.
⭐ *The repair (owner's ruling):* the simulator takes its transcriptome from a **rigel index**, so nascent
RNA is not a space at all — it is the index's own single-exon nascent ENTITY (`create_nrna_transcripts`,
TSS/TES clustered within `NRNA_MERGE_TOLERANCE`), an ordinary transcript in ONE multinomial with the
mature rows, and a probe maps by GENOMIC OVERLAP to gDNA and to every transcript whose exons it touches.
Measured after: gDNA and nascent under one probe are enriched **1162x and 1003x** — the same rate.
⛔ *The general rule:* a simulated population must be a first-class member of the same object set as the
others, weighted the same way and reached by the same mechanisms. A private space is where a private
physics hides, and the tool then inherits the error as if it were its own.
⚠ *And the reason it took a currency measurement to find it*: every gate the simulator had was WITHIN a
space (nascent totals, nascent lengths) and none compared ACROSS them at one genomic position.
*Sibling:* `real-data-is-a-test-input`, `a-single-level-panel-cannot-see-a-constant`, `prove-the-substrate`.


---

## G. Process

**no-magic-numbers. No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every
divisor must be derived from the deposit rule and unit-tested against brute-force enumeration.

**one-thing-varied. One thing varied per experiment**, with the falsification test written first and verified failing,
and a baseline re-recorded from the current tree in the same session.

**converge-and-delete. Converge and delete.** No legacy, no backwards compatibility, no speculative code. Code kept "for
comparison with the old version" is a defect. No version suffixes in file names.

**the-source-does-not-cite-docs. The source does not reference the docs.** Docs evolve and rot — 73 % of the citations that used to
be in the source pointed at documents that had already been deleted. A docstring may cite a **test** or
an executable specification, because code cannot rot silently; it may not cite a document.

**running-an-arm-is-a-fresh-process. ⚠ SIX OPERATIONAL TRAPS FROM RUNNING PANEL ARMS, each of which cost a
launch.** Promoted from a working doc when it was deleted (2026-08-07). They are not
about the model; they are about not wasting an hour.
(i) ⛔ **Never edit `src/` while an arm is running** — every arm is a fresh process that imports `src` at
start-up, so a mid-flight edit silently changes what half the shards measured. Adding a NEW arm to
`scripts/` is safe. (ii) ⛔ **zsh does not word-split an unquoted variable** — use an array and
`"${CONDS[@]}"`. (iii) ⛔ **A wait-loop whose `pgrep` pattern matches its own wrapper deadlocks** — wait on
a log marker instead. (iv) ✅ **`pass0_vs_oracle.DEFAULT_SUITE` pointed at the deleted PILOT panel until
2026-08-15** — three instruments read it from there and all three died on "Failed to open BAM" in their
no-argument form. Repaired to `ladder`; `prior_vs_oracle.py` carries a second copy that was already
correct, which is how one of two homes goes stale unnoticed. (v) ⚠ **Node-axis and region+boundary figures
differ by ~2×** — say which one, every time. (vi) ⚠ **A composite arm fires only its COMPONENTS' names**, so an
`TRAPS: an-ablation-that-never-ran` guard keyed on the arm's own name trips after a complete, valid run —
check WHY a guard fired before distrusting the data.

**shard-an-arm-sweep-by-condition. Shard a panel sweep by CONDITION, never by ARM over one condition —
the instruments WRITE their caches.** `pass0_vs_oracle.measure_condition` writes
`<oracle_cache>/<condition>/_main` whenever `--oracle-cache` is passed, and `load_or_build_oracle` writes
the three origin partitions when they are missing. Running three arms concurrently on the SAME condition
therefore races three writers on one directory: measured 2026-08-15, one `payload.npz` came out truncated
and every later read of it raised `BadZipFile`. ⭐ **The damage is loud rather than silent, and that is
worth knowing before panicking**: a torn zip raises, `read_scan_cache` re-verifies four digests, and
`OracleTruth.from_parts` re-runs sum-to-full on every bank — so a race cannot quietly feed a wrong truth,
it can only kill the run. The fifteen survivors were checked against `scan_cache/<condition>` and were
bit-identical on every INTEGER bank and **3.513e-14** apart on the six float ones, exactly the recorded
re-association floor. ⚠ Different arms on DIFFERENT conditions is safe and is the shard that parallelises;
so is one arm sequentially over the panel. ⛔ And an instrument that only READS should take its full
payload from `scan_cache/`, not from a directory another instrument writes.

**no-enumeration-without-a-census. Do not re-propose path or cell enumeration without a memory census.** Possible unspliced paths ≈
1 M regions × 3–6 reachable ends at ~100 B each = 0.3–0.6 GB, plus spliced paths. It was killed by memory,
and no consumer needs it.
