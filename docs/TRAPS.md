# TRAPS — mistakes this project has already made

⭐ **Read before designing anything.** Every rule below was paid for: something was built, measured,
and
came out wrong in a way that was not obvious in advance. They are LESSONS, not measurements — a
number
that says where the tool *is* lives in `ROADMAP.md` §0.

⛔ **CITE A RULE BY ITS NAME**, e.g. `TRAPS: panel-before-src`. The name is the identifier, so a
citation
says what it means without a lookup and is still one greppable string with one home.
`tests/test_no_jargon_labels.py` enforces it and refuses the numbered labels these rules used to
carry.

⚠ **A rule names a DEFECT, never a licence.** Citing one does not authorise a mechanism, and a rule is not
a veto. ⛔ Several were measured on a panel whose nascent RNA was a UNIFORM ratio on every multi-exon
transcript; since 2026-08-22 it is SPARSE (`DESIGN.md` §0b), and at a development STRESS level rather than
a realistic one — so such a rule says what is ROBUST, not what is likely. Where a rule's substrate matters
its entry says so; where it does not, treat the rule as general and the numbers inside it as history.

⛔ **The shape is load-bearing**: an entry starts at column 0 as `**name.` and nothing else may start
a
line with bold — a bold number there manufactures a phantom rule and has broken the uniqueness gate
three
times, most recently during this file's own compaction (2026-08-22).

## THE INDEX — every rule, by section

**A. Validation and gates** — `self-checking-validator` · `perturb-every-gate` · `a-field-driven-gate-is-atomic` · `waive-with-a-measurement` · `a-docstring-that-misdescribes-the-graph` · `a-flat-pile-is-not-a-knot` · `the-rename-that-corrupted-a-diagram` · `a-gate-that-already-passed` · `right-conditional-wrong-marginal` · `byte-identity-gate` · `the-deliverable-is-not-reproducible-by-default` · `a-clip-hides-a-scale-error` · `an-inverted-clip-is-a-constant` · `an-ablation-that-never-ran` · `a-green-suite-hid-five-dead-instruments` · `compatibility-is-geometry-not-composition` · `a-zero-count-is-a-measurement` · `a-ratio-cannot-carry-zero` · `the-divergence-was-a-barrier` · `deadband-from-the-wrong-sample` · `honesty-metrics-reward-ignorance` · `predicate-contradicts-its-docstring` · `a-test-that-redefines` · `a-gates-power-is-its-invariant-set` · `a-gate-that-reconstructs` · `off-grid-message-mode` · `a-comment-quoted-as-a-finding` · `the-intermediate-is-not-the-deliverable` · `could-the-arm-have-fired` · `prove-the-substrate` · `can-the-benchmark-resolve-it` · `toys-rank-hotspots-backwards`

**B. Measurement and inference** — `measure-the-ceiling-first` · `a-broad-population-carries-no-prior` · `attribution-must-survive-a-shuffle` · `score-against-truth` · `zero-target-guards-are-one-sided` · `hard-labels-miss-soft-change` · `never-pool-the-strata` · `a-threshold-on-a-fitted-residue` · `excluding-a-population-hides-it` · `name-the-observable-per-site` · `starved-is-not-depleted` · `the-substrate-knob-fought-back` · `key-on-a-realised-quantity` · `price-the-halves-separately` · `panel-before-src` · `capture-inverts-the-counted-side` · `admitting-an-object-costs` · `substitution-understates-a-source` · `a-symptom-is-not-a-second-defect` · `a-locked-object-is-not-a-control` · `draining-breaks-the-oracle` · `an-equal-length-panel-defeats-the-lift` · `a-length-gap-bypasses-calibration` · `weight-it-like-the-consumer` · `a-support-ceiling-is-the-clamp` · `log-variance-is-not-linear` · `re-record-the-baseline` · `a-truth-table-of-aggregates` · `a-single-level-panel-cannot-see-a-constant` · `score-the-consumers-own-count` · `the-floor-must-reproduce-the-selection`

**C. Pools, selections and divisors** — `a-cancellation-is-conditional-on-its-support` · `a-better-estimator-inside-a-weak-consumer-moves-nothing` · `a-pooled-rate-cannot-see-a-short-object-factor` · `two-estimators-of-one-rate-weight-the-field-differently` · `state-the-population-rule-do-not-inherit-it-from-a-table` · `two-divisors-opposite-sign` · `frame-free-is-not-assumption-free` · `a-purity-filter-is-a-length-filter` · `pure-and-length-censored` · `divide-by-a-probability` · `opposite-tilts-must-not-pool` · `a-mean-of-ratios-inherits-the-partition` · `a-trap-names-the-defect-not-the-repair` · `a-stale-gate-accuses-the-newest-change` · `an-upper-bound-is-not-an-estimate` · `a-gate-on-the-helper-is-not-a-gate-on-the-caller` · `fractional-mass-is-the-problem` · `conservation-misses-mis-attribution` · `a-guard-outlives-its-divisor` · `a-fold-grows-a-heuristic` · `a-ratio-needs-a-population-that-can-supply-its-numerator`

**D. Estimation and solver design** — `we-keep-re-deriving-message-passing` · `one-hop-lifted-out-is-still-the-relay` · `a-variance-cannot-fix-a-bias` · `two-gaussians-one-latent` · `variance-fitted-on-the-belief` · `a-message-from-the-destinations-belief` · `zero-the-precision-with-the-value` · `no-prior-means-haldane` · `prefer-shares-to-differences` · `an-all-zero-factor-is-inert` · `density-below-one-fragment-length` · `identical-paralogs-are-bimodal` · `a-mean-hits-the-mass-weighted-centre-by-luck` · `a-clamp-at-the-closed-end-escapes-the-window` · `the-deconvolution-is-as-good-as-the-density-it-is-handed` · `deriving-one-coordinate-propagates-its-error` · `interpolate-on-the-axis-where-the-lattice-is-uniform` · `read-the-whole-failure-list` · `a-priors-curvature-is-not-the-datas-information` · `a-refutability-test-needs-the-refuting-channel-in-the-fixture` · `a-strength-is-a-nat-a-prior-weight-is-a-count` · `a-four-decimal-print-is-not-a-zero` · `a-constant-in-exact-arithmetic-is-not-constant-in-float64` · `a-toy-and-a-panel-can-disagree-in-rank` · `a-rescale-that-reads-the-source-belief-is-unbounded` · `a-face-total-is-not-a-total-without-its-flux` · `an-imputation-must-cost-something-every-hop` · `a-floored-knob-is-not-the-bandwidth` · `a-mode-count-is-not-a-well-posed-quantity` · `measure-a-default-flip-before-you-write-it`

**E. Structure, indexes and plumbing** — `one-reference-hides-refid-bugs` · `annotated-is-not-genomic` · `an-sj-is-not-a-gap` · `deposit-at-the-sj` · `splicing-makes-the-graph-cyclic` · `nrna-does-not-mean-synthetic` · `credit-exactly-one-sj` · `strand-completes-the-sj-key` · `a-hash-that-misses-its-artifact` · `integer-channels-reproduce` · `worktrees-run-the-wrong-code` · `checkout-deletes-uncommitted-work` · `two-masks-one-name` · `two-docstrings-one-quantity` · `a-transcript-predicate-must-not-silently-drop-a-molecule` · `an-object-class-does-not-see-a-terminus`

**F. Domain facts that read like defects** — `specificity-and-sense-are-complements` · `strand-measures-the-tilt` · `a-linear-likelihood-emits-a-sign` · `amplitude-fades-influence-does-not` · `a-pooled-conversion-applied-per-component` · `an-identity-with-a-qualifier` · `equal-lengths-carry-no-composition` · `capture-is-1000x-on-exons` · `capture-selects-for-length` · `on-target-by-start-is-geometry` · `real-data-is-a-test-input` · `eff-lengths-do-not-cancel-at-an-end` · `configured-lengths-are-not-realised` · `mature-rna-never-crosses-a-boundary` · `a-boundary-with-rna-is-not-an-sj` · `the-panel-enriches-nascent-by-its-own-probes`

**G. Process** — `no-magic-numbers` · `one-thing-varied` · `converge-and-delete` · `the-source-does-not-cite-docs` · `running-an-arm-is-a-fresh-process` · `shard-an-arm-sweep-by-condition` · `no-enumeration-without-a-census`

---

## A. Validation and gates

**self-checking-validator. A validator that calls the builder's own helper validates nothing — worse
than no check, because it reads as one.** Only re-deriving the result by a *different algorithm*
caught a deleted coordinate swap that a 1,289-test suite waved through. Emit all classes always (a
validator comparing only non-empty classes never sees a spurious flag), and prove every validator
fires by corrupting its input.

**perturb-every-gate. Writing the falsification test first is half the discipline; the other half is
breaking the fixed code and watching each gate fire.** Perturbation has found holes in already-green
gates repeatedly — once 7 of 9, including a gate that did not fire on its own named perturbation
because a redundant backstop was silently doing the guard's job.

**a-field-driven-gate-is-atomic. A gate that enumerates the specification's fields fuses spec,
implementation and schema into one commit — plan the change that way or the tree is red in
between.** Adding one field to the specification turned 10 tests red across two files, and the
digest bump invalidated every cached scan — budget the rebuild as part of the commit. The cheap
probe: add the field, run the suite, read the failure count, revert.

**waive-with-a-measurement. An assertion the shipped code violates is waived WITH its measurement,
never widened.** Measure the rate before touching the predicate: either the predicate mis-stated the
coordinate's domain (fix it, from the shipped builder) or the code is wrong (waive it, with the
number and a written reason per waiver). Five such gates became six measurements on their first run,
the largest being that 61.1 % of live message packets over-claim the slot's observed fragments.
Derive any tolerance from the coordinate's own resolution — one grid spacing — never tune it.

**a-docstring-that-misdescribes-the-graph. A claim about the import graph, inside a docstring, that
nothing gates rots exactly like a stale doc citation, one layer down.** Measured: 14 module
docstrings named a sibling with no import boundary and 6 were genuinely stale. The graph is in the
AST, so the prose can be checked — `scripts/design/module_census.py` does; its count is a worklist,
not a verdict, because a data-flow claim only a human can judge.

**a-flat-pile-is-not-a-knot. Before merging modules, ask whether the problem is entanglement or
missing order — merging a flat pile makes bigger files with the same problem.** The calibration
package read as unmaintainable yet had zero import cycles and 18 of 35 modules with exactly one
importer; naming the layers already in the boundaries cost zero behaviour change and immediately
found two upward imports.

**the-rename-that-corrupted-a-diagram. A mechanical rewrite over prose hits tokens that are not
prose, and the ones it hits are the ones you did not enumerate — write the forbidding gate FIRST,
then treat its residue as the list of collisions worth understanding rather than exempting.** A
citation rename rewrote 20 slot ids inside ASCII chain diagrams because the same token was a label
in one namespace and a slot in another. *Sibling:* `TRAPS: two-masks-one-name`.

**a-gate-that-already-passed. A gate that already passes is not a falsification.** A simulator gate
passed with the defect present because the per-fragment conditional was correct and only the
marginal was being discarded. Before trusting a gate set, check which of them would have failed
before the fix; if none, the gate has not been written yet.

**right-conditional-wrong-marginal. When every per-observation conditional is correct, no
conditional comparison can detect a wrong marginal — only a check on the marginal fails.** The
general form of `TRAPS: a-gate-that-already-passed`; it tells you *where* to look.

**byte-identity-gate. A bit-identity gate has lied in both directions: an arm with ZERO rows scored
"identical" because the loop ran over the new arm's rows, and a stale stored baseline made
unmodified HEAD read as broken.** Require equal row-key sets, and re-record the baseline from the
current tree in the same session — if HEAD-vs-baseline is not 100 %, the baseline is what is broken.

**the-deliverable-is-not-reproducible-by-default. The shipped EM seed defaults to none with sampled
hard assignment, so an end-to-end A/B on the default config measures its effect plus a sampling
draw.** Two identical runs differed by up to 43 fragments on the gate toy. Pin the seed (or use
fractional assignment) on every measurement arm and print a reseeded noise floor beside the effect;
whether the default should change is an owner call.

**a-clip-hides-a-scale-error. A `min()` clip hid an exact factor of 2 for months** — the fixtures
cancelled it exactly, and under a uniform field the clip returned the unscaled value, so the bedrock
invariant written to catch exactly this passed with the bug present. Repair fixtures, never relax
assertions; an exact algebraic 2 is never a modelling approximation.

**an-inverted-clip-is-a-constant. A two-endpoint clip has two silent degeneracies — reversed
endpoints and equal endpoints — and both turn the estimator into a constant, which at a zero control
is the right answer.** A zero control passed by a degenerate estimator may not be quoted as evidence
the estimator DISCRIMINATES: require the same estimator to MOVE at a contaminated condition, and
report the interval's width at the control, because zero width is the tell. *Sibling:* `TRAPS:
a-clip-hides-a-scale-error`.

**an-ablation-that-never-ran. An ablation that never ran reads as "no effect" — and an inert
PROTOTYPE arm also reads as a refutation of the idea, a negative nobody audits because it looks like
the result you were braced for.** Three forms, each measured:
① a patch site must assert the name it patches still EXISTS — `from pkg import name` binds a
re-exported function, not the module (patch via `sys.modules[...]`), and a rename turns a
monkeypatch into a silent no-op;
② every ablation counts its own firings and the harness raises on zero, before any number is read;
③ the counter is not enough, because the config is part of the arm's definition: under one shipped
default, 22 of 26 arms fired on tens of thousands of slots and still scored byte-identical to base —
so require a diff against a base run under the SAME config, prove each arm INERT under one policy
and MOVED under the other, and stamp the setting into every row.

**a-green-suite-hid-five-dead-instruments. An instrument is alive only under the configuration it is
RUN in, and a green suite says nothing about it — the suite's own tests may install a policy the
instruments never see.** Five instruments died at once on a `src/` deletion; six more died on a
config default flip with nothing deleted or renamed, reading diagnostic keys the shipped default
never writes; a later flip broke a byte-identity arm because the flip changed which fields were
load-bearing. The trigger list: run the instruments — not just the suite — after a `src/` DELETION,
after a RENAME, and after flipping any SHIPPED DEFAULT; re-run every byte-identity gate after a
default flip; and re-derive a rebuilt quantity by the producer's own expression, never an algebraic
equivalent — one ULP is enough once something consumes the field.

**compatibility-is-geometry-not-composition. Any RNA-vs-gDNA question resolved by SET MEMBERSHIP
will be answered by geometry, because gDNA is compatible with everything unspliced — resolve it by
LIKELIHOOD, where the populations actually differ (strand, length, density, splicing).** A nascent
warm-start gate keyed on RNA-unambiguous support revived MORE nascent entities the MORE gDNA the
library carried — the exact inverse of the intent — and the statistic that motivated it was
structurally zero on every input. The refused gate aligns with the nascent scope ruling (`DESIGN.md`
§0b): nascent entities get no prior mass and earn their place by likelihood.

**a-zero-count-is-a-measurement. A zero count is a measurement of a density, not an absence of data
— keying precision on the count makes the strongest statement in the library the quietest.** Zero
fragments over 50.7 Mb and zero over 200 bp are the same number to `p = n/(n·Var+1)` and opposite
statements about the world. The derived repair needs no constant: a Poisson rate with `a` events
over exposure `E` under the Jeffreys prior has posterior `Gamma(a+½, E)` — proper, with finite
precision at `a = 0`.

**a-ratio-cannot-carry-zero. A multiplicative transport cannot carry zero — "there is none here" is
unrepresentable by a density RATIO, however precisely the source knows it.** The repair to the rule
above measured +0.0 % with the zero-gDNA rows byte-identical, because every guard on the relay path
requires a positive density, so a zero-density source cannot participate at all. Before fixing a
source that will not speak, check whether the channel can carry the value it wants to send; the
source fix is a precondition, priced as what it measured — nothing.

**the-divergence-was-a-barrier. Before crediting a divergence's removal, ask what the divergence was
suppressing — an infinity in a damping term is a structural gate wearing a variance's clothes.**
Making a transfer variance finite turned every zero-mass slot from a chain-cutting barrier into an
unscaled conduit, and one coherent stratum got 20–34 % worse: the infinity had been the only thing
pricing the premise "a gDNA level does not change across a boundary", which capture falsifies by
~350×.

**deadband-from-the-wrong-sample. A noise deadband whose cushion is supplied by an unrelated sample
size fails exactly where that sample gets big — and it fails silently, into the honesty columns.**
The strand gate's cushion shrank as gDNA grew, the phantom information then scaled linearly with
depth, and the damage was not accuracy (a wash) but the solvable fraction collapsing to zero —
inflating the very column used to pick the worst condition. Propagate a variance instead of gating
on it; the derived repair fixes the honesty columns, not the deliverable. *Sibling:* `TRAPS:
honesty-metrics-reward-ignorance`.

**honesty-metrics-reward-ignorance. Every honesty metric improves as the solver stops knowing
anything, so none is readable without a fixed-denominator accuracy number beside it.** A destruction
control made accuracy 792 % worse while three of the four headline honesty columns improved on every
condition. An A/B that moves the solvable fraction has changed its own denominator: quote the
mass-weighted error over ALL objects and the raw Σ|err| alongside — those two cannot be gamed by
knowing less.

**predicate-contradicts-its-docstring. A predicate can contradict its own docstring for months if
something else masks it, and the wrong version can be load-bearing — when docstring and code
disagree, find out which one the panel relies on before "fixing" either.** The mis-scoped lock
declared 18,397 empty exons and introns composition-certain, yet scoping it as documented was
panel-negative on every re-pricing. STAMP: those re-pricings were measured on rows carrying the
panel's 20 % stress nascent share — re-price before treating the verdict as binding at realistic
nascent levels. The universal half: a test that hands the predicate in as an ARGUMENT cannot see the
caller compute it wrong — gate it through the PRODUCTION path.

**a-test-that-redefines. A test that re-derives a definition cannot detect drift in it.** A gate
meant to keep two instruments' shared class definition aligned recomputed it inline, so changing one
instrument fired nothing. A shared definition needs ONE home — production code — with every consumer
importing it. *Sibling:* `TRAPS: self-checking-validator`.

**a-gates-power-is-its-invariant-set. A comparison gate's power is the size of its invariant set, so
a change that shrinks that set must not ride along with the change it is meant to gate.** Bundling a
~1,000-site rename into a schema re-scan would have left almost no shared banks for the identity
gate to compare — still green, proving nothing — while every refused cache also left the instruments
blind. The tell: before bundling, count how many objects the gate will still compare afterwards;
"almost none" means the gate has been disarmed, not satisfied.

**a-gate-that-reconstructs. A gate that reconstructs a value is vacuous wherever the code decides
not to use it — gate a symmetry the code cannot fake instead.** Eleven green relay gates passed the
deliberately mis-paired run, because on that fixture the reconstructed ratio never entered the
observable. When a gate needs a value to flow through a conditional path, it is a gate on that
conditional; prefer an invariance — here, mirror-image passes on a palindromic chain.

**off-grid-message-mode. A message mode outside its grid's domain is not a weak claim, it is a pin
at the boundary — the penalty has no interior minimum, so precision buys a corner, and an
out-of-range mode is the MOST confident statement a channel can make.** One such unit error caused
74 % of the zero-gDNA control's error. Two corollaries: every channel a solver delivers must appear
in its debug capture, or no instrument can rank it and its absence reads as innocence; and build a
message in the destination's OWN coordinates — correcting the coordinate alone was not the fix,
because a raw count tilt inverts sign on an antisense protocol.

**a-comment-quoted-as-a-finding. A word borrowed from a code comment is not a measurement — one
propagated into the docs invented a regression that never happened, and a build plan was correctly
derived from it.** Before repeating a source comment as a finding, check what the comment was
written to justify, and never let a term of art cross from a comment into a design doc without
re-deriving it from the code.

**the-intermediate-is-not-the-deliverable. An intermediate metric is not the deliverable until
somebody measures the coupling** — a −37.2 % pass-0 win was −3.9 % on the shipped solve, and the
same change's regression GREW through the refit, because a pass-0 biased in the prior's own
direction trains the prior to repeat the bias. When an intermediate stage feeds a fitted stage,
report both columns on every arm and never rank on the upstream one alone. The tell was free: the
instrument already emitted the downstream number.

**could-the-arm-have-fired. "The arm changed nothing" is not a control until you check the arm COULD
have changed something — count the opportunities the change had to fire and print that count beside
the result.** A celebrated free falsification tested nothing: 6 of 8 hops per direction were skipped
by a divide-by-zero guard, and the only live hops were ones the change did not touch. A fixture is
an arm too — a one-sj fixture left a deleted deposit rule fully green because no block was ever
bounded twice; ask of every fixture: could this have failed?

**prove-the-substrate. Prove the SUBSTRATE before proving the code — when a simulated axis is the
axis you are judging, gate the simulator on it.** For two milestones the panel's post-capture
fragment-length distribution was byte-identical to its pre-capture one, and everything measured
against it inherited that. The tell was free: diff the two capture arms' truth files.

**can-the-benchmark-resolve-it. Before running a benchmark, prove it can resolve the axis you are
changing.** A suite judged a partition change for months while its fine region set was row-for-row
identical to its merged set. `scripts/design/suite_resolves.py` is the lesson made executable —
every requirement scored against its degenerate value, no tuned thresholds.

**toys-rank-hotspots-backwards. A toy ranks performance hotspots backwards — profile on cached real
data.** At 3.4 k vs 1.5 M regions the prior's EM went from 28 % of runtime to 0.7 %, and a whole
analysis was spent on the toy's top hotspot.

---

---

## B. Measurement and inference

**measure-the-ceiling-first. Measure the ceiling before building the correction — hand the consumer
the exact answer for one channel and see what perfecting it is worth.** One channel was ranked first
for two sessions; its ceiling was worth ~1 % while an unranked channel was worth 21 %. An A/B says
whether a change helped; a ceiling says whether the work is worth starting.
*Instrument:* `scripts/design/calibration_truth_ab.py --ceiling`. *Sibling:*
`a-symptom-is-not-a-second-defect` (check what the ceiling can even see).

**a-broad-population-carries-no-prior. A population prior only transfers information if the
population is tight — before fitting one, measure the pooled distribution's width against the
per-object uncertainty, and if it is wider there is no prior to be had.** The gDNA landscape works
because gDNA is near-uniform; the same machinery fitted to RNA, whose per-object density spans ~3
decades, measured 0.988–1.037 against base — nothing.

**attribution-must-survive-a-shuffle. An oracle arm that a scrambled oracle also wins has not priced
the truth — re-run every truth-fed ceiling with the truth shuffled before believing it.** A fitted
oracle RNA density won 0.78–0.85×, but the same truth values shuffled against their opportunities
beat the true shape at g98 (0.786 vs 0.854): the gain was "any shape leaning toward the vertex", not
the truth. One permutation is the whole cost.

**score-against-truth. Score against TRUTH, not against the previous run.** The simulator writes
per-fragment ground truth into the oracle BAM's read names; scoring against it turned one "the
deliverable improved" into "the deliverable got 23.9 % worse".

**zero-target-guards-are-one-sided. A zero-target guard is ONE-SIDED: in a library with no gDNA, any
change that lowers the estimated gDNA fraction scores better.** Score on the contaminated conditions
and use zero-target rows as false-positive checks only. A reported −13.1 % win was +1.4 % worse over
the full battery.

**hard-labels-miss-soft-change. A byte-identical hard-label result is NO EVIDENCE, not no change —
and a net error is only a lower bound on the real one, arbitrarily weak.** `Σ|err| / |net|` reached
274× on one arm (a +672-fragment net over 184,271 misplaced). Report `Σ|err|` beside every net, and
the directional split beside both.

**never-pool-the-strata. Report per stratum with the denominator named — a pooled mean buries
unsolvable objects, hides a sign flip between strata, and lets one huge object set the number.**
Since the 0.8.0 scope there is a sharper reason: the deferred unstranded × capture-ON stratum
carries 64.5 % of transcript error and 90 % of gene-level error, so a pooled total mostly reports
the one stratum nobody is optimising. The same failure has a second axis within a stratum: a score
scoped to a mechanism's own targets is a diagnostic, only all-live-slots is a verdict — a scoped
relay score once read 0.916× (a win) while all live slots read 1.126× (a loss).
*Sibling:* `excluding-a-population-hides-it`.

**a-threshold-on-a-fitted-residue. A binary cut on a fitted parameter's residue is not a population
test, and a better threshold is not the fix — propagate the fitted parameter's own variance
instead.** A `tau > 1e-9` cut promoted objects whose own statement was 1,377 nats wide into the
scored population, hiding a 1.06 M-fragment error; any floor is a tuned constant because the residue
is continuous across the region.

**excluding-a-population-hides-it. Every population you exclude from the denominator needs its own
gate, written at the same time as the exclusion.** A condition reporting mwae 0.0170 on a 39.5
%-solvable set had 1,056,019 fragments of error in the excluded class. The tell is free: the reason
to exclude a class is that its answer should be uninformative, so measure how far from uninformative
it actually is.

**name-the-observable-per-site. For each place a change was made, name the observable that would
move and confirm at least one gate reads it — a count of gates says nothing when they all read one
downstream publication.** Ten sound tests all read the combine's publication, and deleting the
conjunct from one of the sites left the entire calibration suite green.

**starved-is-not-depleted. "Starved" and "depleted" are different diagnoses: starved (the experiment
is too small) is fixed by more data, depleted (the biology puts nothing there) never is.** Decide
which by asking whether the count grows with the lever you have; reading depleted as starved cost
half a session chasing depth.

**the-substrate-knob-fought-back. A substrate knob can be adversarial to the very population under
study — read the sampler's weight for that population's own objects before trusting the panel.** Toy
capture probes tiled in transcript space put a split probe over every internal sj, and the
split-probe gDNA penalty then suppressed exactly the boundary-spanning fragments being measured.

**key-on-a-realised-quantity. Key an operating point on a REALISED quantity visible in the truth,
never on a knob you set — and print the realised value so a row that cannot reach the target says
so.** An RNA "level" set as a multiple of the off-target gDNA density gave a true `f_g` of 0.010 off
capture and 0.31 under it at the same knob setting, because capture concentrates gDNA onto the exon
~65×.

**price-the-halves-separately. A diagnosed defect and a valuable one are not the same defect — price
the halves separately.** Two length models priced together read −6.31 % and looked like one finding;
split, they read −5.90 % gDNA / −0.43 % RNA, and the fully diagnosed half was the near-worthless one
(14× apart).

**panel-before-src. A toy ceiling is not a panel ceiling — run the panel arm before writing a
mechanism into `src/`, and make every claim about a mechanism name its substrate.** Recorded five
times now: toy-positive changes that were panel-negative, including a −0.009 toy win carrying to the
ladder as +0.0013 worse. When the panel disagrees, ask whether the defect is half of a cancelling
pair before concluding.
*Sibling:* `a-cancelling-defect-pair`.

**capture-inverts-the-counted-side. "The well-counted side" is not a fixed side — capture inverts
it, so a rule that transports from the well-counted side silently reverses direction with the
protocol.** Off capture an intron holds ~349 gDNA counts against a flanking boundary's 8–36; under
capture the intron holds 1 and the boundary 20–40.

**admitting-an-object-costs. Admitting an object to the scored population is a cost, and a mechanism
can do it silently — report the solvable fraction beside every accuracy number.** A prototype raised
solvability 43.1 → 48.2 % while the edge-axis error went 0.0216 → 0.1449 — not because answers got
worse but because wrong ones started being counted.
*Sibling:* `excluding-a-population-hides-it` (inverted).

**substitution-understates-a-source. A ceiling by substitution is honest for a SINK but understates
a message SOURCE, whose value is what it carries — for a source, pin the truth and RE-SOLVE.**
Substituting two boundaries removed 9.1 % of a gene's error while the object they feed accounted for
82.7 %. An instrument offering substitution must say which of the two it is doing.

**a-symptom-is-not-a-second-defect. A downstream number that is wrong because its input is wrong is
a symptom, not a second defect — the proof is to substitute only the upstream input and re-run the
SHIPPED downstream function.** The effective-length shrinkage read as an independent bug (factor
0.345 at the zero-gDNA control against a correct 1.000); fed truth composition, the shipped function
returned 1.000 exactly, so "repairing" it would have installed a second error to cancel the first.
The substitution is admissible only where the downstream is a pure function of the substituted input
with no feedback; otherwise pin-and-re-solve (`substitution-understates-a-source`). The second half
is why nobody had seen it: a ceiling prices only what is built AFTER its injection point — every
measurement arm patched `assemble_priors`, which runs after the shrinkage is built, so it cancelled
out of every ceiling ever run. For every ceiling, write down where the injection happens and what
was already computed above it.
*Sibling:* `a-cancelling-defect-pair`, `measure-the-ceiling-first`.

**a-locked-object-is-not-a-control. A structurally locked object keeps its pinned init and is never
solved, so it cannot be wrong and its correctness measures nothing.** The control that works is the
same class split on the variable — boundaries with vs without sj flux — and it reversed the
conclusion the "healthy twin" reading had supported.

**draining-breaks-the-oracle. A per-fragment-independent partition stops being one the moment a
downstream step conditions on the whole tally — the identity that makes a truth source valid also
fixes where in the pipeline it can be taken.** The second pass scores against the payload's own
densities, so three origin partitions drained separately are not the whole drained; measure the
undrained stage rather than draining the parts and hoping.

**an-equal-length-panel-defeats-the-lift. On an equal-length panel the drained-oracle lift's
span-tie ambiguity is the common case by construction, so run both sides UNDRAINED and price the
caveat instead of asserting it small.** Measured: 3.92 % of held fragments ambiguous, with spliced
records replayed into the gdna partition — a refusal the oracle is right to make. Equal lengths are
not the panel's mistake: see `a-length-gap-bypasses-calibration`, this rule's mirror.

**a-length-gap-bypasses-calibration. A large gDNA-vs-RNA fragment-length gap lets the EM assign
fragments on length alone, so the tool answers correctly with calibration broken — the EM already
reads the fragment-length distribution, and equal lengths FORCE calibration to be exercised (owner,
2026-08-14).** The tell is free: score calibration against oracle calibration beside the end-to-end
number; a healthy end-to-end number on top of an unhealthy composition claim is the signature. The
general form: for each axis a substrate varies, ask which stage that axis lets a LATER stage answer
without — a benchmark that hands the deliverable a shortcut measures the shortcut.
*Sibling:* `an-equal-length-panel-defeats-the-lift` (the cost side of the same design choice).

**weight-it-like-the-consumer. Weight the average the way the consumer weights it — a bp-weighted
and a fragment-weighted mean answer different questions.** An "11 % over-call" was bp-weighted; the
estimator is fragment-weighted, 89 % of the mass sat where the effect is inert, and end to end it
moved ≤ 0.0002.

**a-support-ceiling-is-the-clamp. A support ceiling that matches the clamp is not a match — it is
the clamp.** A distribution's support "agreeing" with `max_frag_length` was recorded as a fix; the
narrower estimate had been correct.

**log-variance-is-not-linear. `Var(log f)` is not `Var(f)` — convert with the delta method (`Var(f)
≈ f²·Var(log f)`); at small `f` the two differ by orders of magnitude.** Every overconfidence figure
computed without the conversion is void, not merely imprecise.

**re-record-the-baseline. A delta is only attributable if its baseline came from the same tree in
the same session — re-record the before-picture, never quote a stored one.**

**a-truth-table-of-aggregates. A truth table's label column may hold nested aggregates, and summing
it double-counts — with the control arm unaffected, the most persuasive shape a wrong result can
take.** Enumerate the membership sets explicitly, RAISE on an unrecognised label, and put a premise
gate on any truth source. A parser bucketing `{mrna, nrna, gdna, rna, all}` as "gdna, else rna"
summed mRNA twice plus the whole library once while the gDNA arm read exactly 1.000.

**a-single-level-panel-cannot-see-a-constant. A panel that holds the quantity of interest fixed
cannot distinguish a good estimator from a constant that happens to equal it — list the axes the
panel varies and check the DELIVERABLE is one of them.** A length channel won 7 of 8 conditions on
an all-g50 panel, then reported 54–57 % gDNA on the zero-gDNA control: the win was the panel
agreeing with a near-constant ½. Three gDNA levels is the minimum for this reason. (The channel it
was measured on is retired past 0.8.0; the panel lesson is general.)

**score-the-consumers-own-count. Truth is not a yardstick until you name the consumer — a truth
nothing reads measures your choice of target, in the right units — and the check is per array, not
per instrument: for each number the consumer reads, write down the SET of fragments it counts and
score against that set.** The first-base count `F` had impeccable provenance and was still not the
EM's target (the EM's soft count includes straddling candidates `F` drops, 0.35–3.35 % per
condition); scored against `F`, a 72 % "open residual" spawned a research programme, and against the
EM's own count it vanished. It recurred in the same file that documented it and again in the commit
that wrote this rule, on another array (spliced RNA never competes with gDNA, so a prior arbitrating
that competition must not be scored against a denominator counting it). The tell before any
measurement: write down which line of the consumer reads the number.

**the-floor-must-reproduce-the-selection. Whatever selection the scorer applies to the data, the
noise floor must apply to its draws — a floor drawn unconditionally under a conditioning scorer
reads the conditioning as a rule error.** Scored destinations were zero-truncated Poisson; an
unconditional floor printed a 10–14 % "rule error" on hops that were exact, and redrawing
conditional on `n >= 1` closed it (972 vs 1,011). A floor that cannot be falsified by turning the
conditioning off is not a floor.
*Sibling:* `a-ratio-needs-a-population-that-can-supply-its-numerator`,
`excluding-a-population-hides-it`.

---

## C. Pools, selections and divisors

**a-cancellation-is-conditional-on-its-support. ⭐⭐⭐ A RECIPROCAL-OPPORTUNITY DEPOSIT CANCELS ITS
OPPORTUNITY ONLY WHERE THAT OPPORTUNITY IS NON-ZERO.** Where `A(w) = 0` the fragment deposits
nothing, so
`E[Σ 1/A] = ρ · P(A > 0)` — and `P(A > 0)` is a functional of exactly the length distribution the
channel
claims independence from. ⛔ The tell is a comment that states the support as a REASSURANCE rather
than as a
factor. Measured on the shipped REGION bank `1/(ℓ−w+1)`: an **11.6× under-read at the panel's median
98 bp
exon**, and identically zero for the 7,123 of 24,018 exon REGIONs shorter than `frag_min`, at any
depth.
Write the expectation with its support factor everywhere and gate it with fragments placed outside
the
support; `EQUATIONS.md` §2 carries the corrected algebra.

**a-better-estimator-inside-a-weak-consumer-moves-nothing. ⭐⭐ A 1.8–4.3× REPAIR TO AN ESTIMATOR CAN
BE
WORTH ~1 % END TO END, AND THE CONSUMER'S STRENGTH IS WHY.** Swapping the pooled gDNA background
from
`count / E_contained` to the pmf-free START/END pair improved the estimator 1.8× on the intergenic
pool and
4.3× on the +introns pool under capture, while the deliverable against oracle calibration moved
0.9872 on
stranded × capture-ON and was a wash off capture. A background that enters as a weak one-sided floor
cannot
transmit a 2× better rate. ⭐ The rule is about SEQUENCING: price the consumer's sensitivity before
repairing its input, and spend the effort on whatever reads a level sharply.

**a-pooled-rate-cannot-see-a-short-object-factor. ⭐⭐⭐ A POOLED RATE IS A RATIO OF SUMS, SO IT IS
DOMINATED
BY THE LARGEST OBJECTS — AND A PER-OBJECT FACTOR THAT ONLY BITES AT SMALL ONES IS INVISIBLE IN IT,
AT ANY
SIZE.** The pooled gDNA background differs from the pmf-free form by `ρ_r·(E_r/E_g − 1)`, which
reaches
−99 % at `ℓ = 250` under a large fragment-length gap — yet the background pools did not move
(1.0007 / 1.0000 / 0.9994 across gaps of +169 / −4 / −171 bp, in both directions), because
`ΣE_r/ΣE_g` is
0.9966 over intergenic regions whose exposure is dominated by megabase-scale regions (ℓ-weighted
mean
**3 Mbp** against a median of 8 kb). Only exons are short enough to feel it (`ΣE_r/ΣE_g = 0.7294`,
median ℓ
143 bp) and the exon row moved exactly as the factor predicts. ⭐⭐ Before building a panel to exert a
per-object mechanism, compute the POOLED WEIGHT of the objects it acts on — and read the answer as
saying
which consumer the mechanism belongs to.

**two-estimators-of-one-rate-weight-the-field-differently. ⭐⭐⭐ TWO UNBIASED ESTIMATORS OF ONE SLOT'S
RATE
AGREE ONLY UNDER LOCAL UNIFORMITY; UNDER A PROBE-SHAPED FIELD THEY ARE DIFFERENT WEIGHTED AVERAGES
OF IT —
BOTH RIGHT, AND NOT INTERCHANGEABLE.** The START bank over `ℓ` and the CONTAINED reciprocal bank
reproduce
`P(w ≤ ℓ)` to 1.00025 / 1.00045 at capture-OFF over 14.4 M gDNA fragments; the same comparison under
capture reads **0.906 pooled and 0.034–0.077 on intron/intergenic regions**, because a START is
counted for
a fragment beginning anywhere in the region while a CONTAINED fragment must lie wholly inside it. ⛔
So
swapping one for the other under capture CHANGES THE QUANTITY rather than removing a bias, and may
never be
priced as a drop-in. ⭐ The escape is a pair that shares the weighting: START and END have the same
opportunity `ℓ`, so their agreement is field-free — both-exact regions read 1.0018–1.0057 while a
bound
side reads 0.712–0.879 or else 1.242–1.428, in the predicted direction (pre-rebuild figures; on the sparse-nascent panel 0.6974–0.9955 / 1.0216–1.4878, same direction). Build the field-free
comparison; do
not argue about the field-dependent one.

**state-the-population-rule-do-not-inherit-it-from-a-table. ⭐⭐ A POPULATION DEFINED BY WHICHEVER
ROWS A
FILE HAPPENS TO CONTAIN IS NOT DEFINED, AND TWO IMPLEMENTATIONS WILL REACH IT BY DIFFERENT DOORS.**
The
mature wall distances inherited "real transcripts only" from `intervals.feather` carrying no
synthetic
transcript; an independent enumeration through `index.get_exon_intervals`, which does return a
synthetic
span's own interval, disagreed on **57 distances**. ⛔ The tell is a filter you did not write: if you
cannot
point at the line that excludes a population, you are relying on a file's contents. Write the filter
and
the sentence saying why in BOTH implementations.

**two-divisors-opposite-sign. ⭐⭐ TWO DIVISORS BUILT FROM ONE pmf CAN STILL DISAGREE — IF THEY
RESPOND TO IT
WITH OPPOSITE SIGN.** `E_J = E[w]−1` rises with the mean fragment length while `E_r = e−E[w]+1`
falls, so a
length-model error appears as a junction-vs-exon frame gap of 0.62 %/bp: the two are exactly
consistent
(202.8 + 797.2 = 1000.0) and still 10 % apart. ⭐ When two quantities built from one model disagree,
differentiate both with respect to that model before looking for a bug in either. ⛔ And read the
simulator's own sampling code, not its docstring, before pricing a selection effect in a simulated
substrate — a 1.024 "finite-transcript placement" factor added here was sound arithmetic about a
generative
model the simulator does not use, and cancels exactly.

**frame-free-is-not-assumption-free. ⭐⭐⭐ A TERM THAT IS FRAME-FREE IS NOT THEREFORE ASSUMPTION-FREE
— LOOK
AT THE NUISANCE YOU PROFILED OUT, NOT ONLY AT THE DIVISOR THAT CANCELLED.** Putting certified-RNA
`boundary_spliced` on ψ's RNA arm as `E[S] = c·(1−f_g)·M` is correct about the frame and still not
enough:
the same `c` holds the splice-visibility `q`, whose dropped term is the same size as the term kept,
and
with `q` free the profile likelihood in `f_g` is exactly flat on [0,1) — the whole channel carries
one bit,
"`f_g ≠ 1`". The raw-count term scored worse than the uninformative reference on 12 of 36 ladder
conditions, worst +0.4578. ⭐⭐ The tell was free and general: its benefit tracked the ANSWER, not the
evidence (−0.49 where the truth is all-RNA, +0.45 where it is all-gDNA), and a channel whose sign
follows
the truth is a prior. Ask of any one-sided floor whether the term you drop is bounded or dominant.
*Sibling:* `two-divisors-opposite-sign`.

**a-purity-filter-is-a-length-filter. A PURITY FILTER ON A LENGTH POOL IS A LENGTH FILTER.** Barring
fragments whose length was partly *inferred* rather than sequenced selects exactly the ones whose
mates sit
far apart: the pool read **−9.6 % mean / −22.5 % sd** against truth where keeping them read +0.7 % /
+2.4 %.
Before excluding a population from a pool, ask what the exclusion criterion correlates with — if it
correlates with the axis being measured, purity and accuracy point in opposite directions. ⭐ The
rule that
replaced it: a pool is keyed on **determinacy, not provenance**.

**pure-and-length-censored. A POOL CAN BE COMPOSITION-PURE AND LENGTH-CENSORED AT THE SAME TIME, AND
THE
SECOND IS INVISIBLE TO EVERY PURITY ARGUMENT.** "gDNA contained in an intergenic or intronic region"
is
100 % gDNA by construction and, under hybrid capture, ~15 % short — a long fragment beside a probe
reaches
the exon boundary and stops being *contained*. The pools it spills into resemble the RNA pool, and
it was
filed for two milestones as "not gDNA". Ask what a pool's selection rule correlates with, not only
what the
selected fragments are.

**divide-by-a-probability. A POOL DIVIDED BY ITS OPPORTUNITY MUST BE DIVIDED BY A PROBABILITY, NOT
BY A
COUNT.** `count(w)/A(w)` recovers the distribution lengths were **drawn** from; every consumer needs
the one
the library **realises**, so the divisor is `A(w)/T(w)`. The two forms differ in shape, and the
ratio form
is also where an abundance weighting cancels — swept to pathological regimes it is never worse than
not
correcting, and the `A`-only form is.

**opposite-tilts-must-not-pool. POOLS WITH OPPOSITE TILTS MUST NOT BE POOLED RAW.** A *contained*
pool's
opportunity falls with length while a *crossing* pool's rises; summing the histograms and applying
one
divisor read a gDNA mean of **146.05** where the contained pool alone said **88.0**. ⭐ Summing the
counts
*and* the matching per-pool opportunities is a different operation and is correct — it is the
opportunity-weighted average of the per-pool estimates, and under Poisson counts those weights are
exactly
inverse-variance.

**a-mean-of-ratios-inherits-the-partition. ⭐⭐ AN ESTIMATOR THAT AVERAGES PER-OBJECT RATIOS IS A
FUNCTION OF
WHERE THE ANNOTATION DREW ITS BOUNDARIES.** A region boundary appears wherever a signature changes,
so how
many objects a transcript's path contains is an artefact rather than biology; `Σmass / Σopportunity`
is
invariant to it because splitting an object splits both sums together. Subdividing one 20 kb intron
from
1 region to 4 to 10 at constant total mass and opportunity moved a shadow span's share of its
locus's prior
9.4 % → 18.6 % → 32.4 %, while the pooled form gave the same answer at all three. ⭐ An estimator
must
be invariant to a re-partition that leaves the data unchanged, and that is a property you can test
directly.

**a-trap-names-the-defect-not-the-repair. ⛔⛔ CITING A TRAP DOES NOT LICENSE A MECHANISM — THE TRAP
NAMES
THE DEFECT, NOT WHICH QUANTITY TO REPAIR.** The rule that a zero count over a large exposure is a
measurement was cited to license a Jeffreys posterior MEAN as a density LOCATION, `(mass +
½)/opportunity`;
that same rule's shipped repair had put the half on the PRECISION (`trigamma(n+½)`, +0.0 % with the
zero-control rows byte-identical), and the location form sits in the graveyard at +7,269 % on `g00`.
Three
components each gaining `+½` breaks `Σ_c ρ_c·E_c = M` by exactly 3/2, and the sentence to remember
is that
the half belongs to the rate's VARIANCE, not to a share of a total. ⛔ Before adopting a mechanism
because a
trap seems to motivate it, grep the graveyard for the mechanism.

**a-stale-gate-accuses-the-newest-change. ⛔⛔ A GATE WHOSE PREMISE EXPIRED DOES NOT GO QUIET — IT
FAILS, AND
IT BLAMES WHATEVER IS IN FLIGHT.** `rescan_panels.py` gates an irreversible cache rebuild on
byte-identity
because "these are integer tallies"; the one-numeric-convention ruling made six banks float64, float
addition is not associative across worker threads, and the gate became unsatisfiable that day. It
then
failed a schema change on five banks that change had not touched, and the natural reading — "my
change had
a side effect" — was wrong. ⭐⭐ What resolved it was a control on the INSTRUMENT, not on the change:
scan
the same BAM twice with the same binary, and they differ on exactly those six banks by at most
3.5e-14
relative. **When a gate fails on something you did not touch, first ask whether it can still pass at
all.**
⚠ A gate that runs only at an irreversible moment is exercised by nothing in between, so its rot is
invisible until the one day it matters — its own `--self-test` passed throughout, against a fixture
rather
than against the data the comparator had to judge. *Sibling:* `a-guard-outlives-its-divisor`, the
inert and
gentler direction.

**an-upper-bound-is-not-an-estimate. ⛔⛔ A BOUND THAT IS SHARED IS NOT EVIDENCE ABOUT WHO SHARES
IT.** A
transcript's density is bounded by the minimum along its path — true, and useless wherever that
minimum is
attained at an object the transcript does not own: 3,644 of 4,839 silent transcripts (75.3 %) share
at
least one object with an expressed one, inherit the loud neighbour's bound, and are asserted into
existence.
Transcript false-positive mass came out 1.76–2.20× worse than base on all twelve rungs while the
same
weights took gene-level error to 0.40–0.53× — the bound is informative about the UNION and silent
about the
split. ⛔ The tell is that the zero-weight SET was byte-identical across all twelve arms: **when an
estimator's support does not move as you vary it, you are varying the wrong thing.**

**a-gate-on-the-helper-is-not-a-gate-on-the-caller. ⛔⛔ TESTING THE FUNCTION THAT COMPUTES A QUANTITY
DOES
NOT TEST THE CALL THAT ASKS FOR IT.** A weight builder shipped with 29 green gates and two holes,
both
found only by perturbation: the re-partition gates called the power-mean helper directly, so
replacing the
CALLER's opportunity weights with `ones` fired nothing; and the two opportunity arrays were compared
in
their own test on a single-exon transcript, where they are equal by construction, so silently
serving the
wrong one fired nothing. ⭐ Both repairs are one shape — a gate through the PUBLIC entry point, on a
fixture
drawn from a population where the two branches give numerically different answers, never from the
majority
population where they coincide. *Sibling:* `a-mean-of-ratios-inherits-the-partition`, the defect the
first
hole was meant to prevent.

**fractional-mass-is-the-problem. FRACTIONAL MASS IS THE PARTITIONING PROBLEM.** A fragment spanning
4 regions writes six fractional numbers whose values depend on region sizes, purely because a mass
is
conserved; the same fragment's three crossing counts depend on nothing. *Corollary:* multimapper and
ambiguous-path assignment must stay INTEGRAL, or the non-integer observable returns and the count
stops
being a count.

**conservation-misses-mis-attribution. MASS CONSERVATION DOES NOT CATCH MIS-ATTRIBUTION, AND IT IS
WORSE
THAN THAT: AN ENTIRE ALTERNATIVE DEPOSIT RULE CONSERVES.** One fragment can credit the same boundary
side
twice, and credit a boundary it never crossed, with total mass still exactly 1.0. Injecting the
rejected
`1/K` rule — an equal share to every crossed boundary — in place of the shipped `slice_len / L` left
every
conservation gate green: both per-fragment laws, the exhaustiveness sum over 5,193 enumerated
fragments,
and the closed form at a boundary where `K == 1` makes the two rules identical. ⭐ What separated
them was a
re-derivation on a different axis — attribute each fragment BASE to the boundaries bounding its own
region
and divide by `L`, which `1/K` cannot express. ⛔ **If a property is invariant across the rules you
are
choosing between, a gate on that property is not a gate on the choice.** ⚠ And the law is exact in
rationals while the bank is fixed point: "the rule conserves" and "the representation costs ≤ K
quanta" are
two claims, each needing its own assertion.

**a-guard-outlives-its-divisor. ⛔⛔ DELETE THE DIVISOR AND THE GUARD AGAINST IT GOES INERT — WHILE
ITS TEST
KEEPS PASSING.** `_mass_where_there_is_opportunity` existed because `rho = Σm/ΣS` is a rate;
re-basing
`assemble_priors` onto the conserved fragment count removed the division, and the guard's own
perturbation
test went on passing only because its fixture asserted `prior == 0` and the mass really was zero on
the
object it checked. Three things had quietly become true at once: the guard was inert on the path the
test
exercised; it was still load-bearing on the ungated eff-length path, where injecting the floored
divisor
moved `gdna_eff_len` by +19.97 bp at 20 units of stray mass and +44.93 bp at 5,000; and the correct
behaviour on the first path had REVERSED, because a count has nothing to divide by, so dropping
zero-opportunity mass now loses fragments the accumulator really deposited. ⛔ The tell is
structural: a
guard's justification names a specific operation, so when that operation leaves the function, re-ask
whether
the guard still has a subject and where its subject went. A green test is no evidence either way;
only
injecting the defect is.

**a-fold-grows-a-heuristic. ⛔⛔ A QUANTITY FOLDED ONTO ANOTHER AXIS TO FIT A CONSUMER'S INTERFACE
WILL GROW
A HEURISTIC TO REPAIR THE FOLD, AND THE HEURISTIC WILL OUTLIVE EVERY MEMORY OF WHY THE FOLD WAS
THERE.** A
contiguous BOUNDARY has no extent, so `_project_regions_to_loci`'s `overlap / region_size_bp` could
not see
it and each boundary's mass was folded into one flank region; keying left then lost a locus's
far-left
boundary, so an intergenic RE-KEY was patched in — and the patch is what made the attribution worst,
routing
100 % of that boundary into the gene where a median 8,066 bp intergenic flank against a 211 bp locus
flank
puts ~61 % of its mass outside. ⛔ Three tells: a helper whose docstring justifies itself by a
*consumer's*
limitation rather than by the model; a special case that exists to undo the general case; and two
callers of
one object where the fold is right for one and wrong for the other. ⭐ The repair is never a better
heuristic — it is to give the consumer the axis it was missing: projecting a boundary AS a boundary
deleted
five helpers at once, needed no schema change, and was numerically identical on real data, which is
why the
fold had survived. `DESIGN.md` §3.1b.

**a-ratio-needs-a-population-that-can-supply-its-numerator. ⛔⛔ ADDING OPPORTUNITY THAT STRUCTURALLY
CANNOT
CARRY COUNTS MANUFACTURES A DEFICIT OF EXACTLY ITS OWN SHARE, AND IT READS AS A DEFECT IN THE
ESTIMATOR.**
An ad-hoc A/B pooled every reference where the shipped field gate stratifies per gDNA-bearing
reference,
folding the panel's 92 ERCC exon regions — real `eff_gdna`, structurally zero gDNA — into the `R
exon`
class: 2.03 % of that class's opportunity carrying 0 counts, which is precisely the "~2 % `R exon`
deficit"
that was written up as a second bug and has been RETRACTED (scored per reference the class was never
off,
z +0.63 / −0.14). ⭐ A class ratio `n / (rho·E)` is meaningful only over objects that COULD supply
the
numerator; before reading one, ask what fraction of `E` sits on a population the numerator cannot
reach —
and if a shipped gate stratifies, an A/B of it must stratify the same way or it is measuring
something else.

---

## D. Estimation and solver design

**we-keep-re-deriving-message-passing. ⭐⭐⭐ FOUR SEPARATE SESSIONS HAVE INDEPENDENTLY RE-DERIVED
MESSAGE
PASSING FROM SCRATCH AND THE OWNER HAS HAD TO EXPLAIN IT EVERY TIME — the most expensive recurring
cost in
this project, and nobody reasoned wrongly, everybody was slow.** This entry exists so that arriving
is
cheap; it deliberately does not contain the derivation.

⭐⭐ *The tell:* you are reasoning about how an EXON gets its gDNA level — from anchors, a ladder of
rungs, a
pooled reference, local imputation, run-fill, a throttle or a bound — and you have not yet written
the words
"message passing". Stop there; you are already inside it.

⭐ *The endpoint, so you can recognise having ARRIVED:* gDNA is directly measurable only where no
mature
transcript crosses — intergenic REGIONs, intron REGIONs, and the BOUNDARIES against them. An exon's
unspliced mass is gDNA + RNA, which is the unknown itself, so an exon's gDNA level cannot be
measured and
must be IMPUTED from its neighbours along `intron REGION ↔ intron|exon BOUNDARY ↔ exon REGION`.
Imputation
across that chain IS message passing.

⚠ *And with strand-specific data most exons solve directly* (owner, 2026-08-22), so the imputation
problem
is weighted toward unstranded data and AMBIG slots rather than being the general case.

⛔⛔ *Four stale beliefs, each of which reads as "message passing was considered and rejected", each
false in
one line:* "the relay was measured bad, so it is off" — the +154.8 % price was measured with a
named,
confirmed licence bug live on a retired ladder, so re-price rather than inherit; "the mechanism does
not
exist yet" — it is built and switched off, `messages/relay.py`'s SPLICE IN (BOUNDARY → EXON) is
exactly that
hop, behind `CalibrationConfig.message_propagation`; "a simplex vertex is unreachable, so the
shortfall is
irreducible" — that theorem requires a prior with a DENSITY and says nothing about one with an ATOM;
"an
exon can just borrow a density from somewhere" — priced on stranded × capture-ON, a pooled scalar
reference
is 3.90× worse, a nearest-rung one 1.27× and a locally imputed one 1.50×, because capture enrichment
is per
exon and arbitrary.

⛔ Find the settled ruling and its maths with `grep -in 'message passing' docs/DESIGN.md
docs/EQUATIONS.md`
— both files number by insertion, so a section number quoted from here goes stale.

**one-hop-lifted-out-is-still-the-relay. ⛔⛔ LIFTING ONE HOP OUT OF `messages/` DOES NOT MAKE IT
SOMETHING
OTHER THAN MESSAGE PASSING.** Owner's correction: a value computed at one object and consumed at
another is
a message regardless of whether it carries a count, an sj flux or a belief, and "it is only one hop"
is not
a definition either.

⭐ *The tell:* the sentence justifying the new code begins *"this isn't really the relay, because…"*.
A
private copy of a framework you have switched off is not a way of avoiding the switch — it is a
second home
for one mechanism with none of the licences, precisions, damping or gates. Fix the framework's known
bug and
re-price it.

**a-variance-cannot-fix-a-bias. You cannot fix a biased mode with a variance.** Under capture a
counting
estimate was systematically ~2× low but PRECISE — both flanking boundaries sat at the same enriched
boundary
and agreed on the same biased-low density, so the bias was trusted. A disagreement-based variance
model
structurally cannot fix a bias.

**two-gaussians-one-latent. Never hand a solver two Gaussians built from one latent.** A message on
`log f`
and one on `log(1−f)` are rank-1 with correlation exactly −1, so adding their Fisher information is
exactly
2× over-confident, rising to ~7× with deep spliced content.

**variance-fitted-on-the-belief. Never fit a variance on the current, not-yet-solved belief.**
Adjacent
WRONG regions agree, so the variance collapses, the messages turn confident and the error
propagates.
Honesty measured against a wrong belief is not honesty. Same family: any component trained on the
solver's
own output is self-confirming — refit iterations 1→5 went 0.0840 → 0.1056, monotonically worse.

**a-message-from-the-destinations-belief. ⭐⭐⭐ A MESSAGE COMPUTED FROM THE DESTINATION'S OWN BELIEF
CARRIES
ZERO INFORMATION AND CONFIRMS THE DESTINATION.** A message may use the destination's CONSTANTS —
geometry,
lengths — and its OBSERVATIONS; never its BELIEFS, and any "fix" that divides the destination's
belief back
out rebuilds the bug. ⭐ The check is *"is the delivered value independent of the destination's own
state?"*, and a prediction that does not move when the data move by four orders of magnitude is the
tell.

⚠⚠ This one lesson has recurred nine times in different costumes. Each name below is cited from
source and
tests, so the names are kept; the investigations are in git.

| the costume | the lesson that generalises |
|---|---|
| **TRAPS: a-total-density-ratio** | `r = ρ_tot(dst)/ρ_tot(src)` re-creates the parent bug whenever the total is dominated by a component the message is not about — a correct 0.0257 gDNA density delivered as 28.684. A scale factor must be built from the component the claim is ABOUT, never from a total |
| **TRAPS: substitute-the-definitions-first** | substituting the definitions showed `ρ_c(src)·r ≡ φ_c(src)·ρ_tot(dst)` — a pure composition imputation with no level in it, so no corrective factor existed. Two sessions were spent hunting a better `r`. ⭐ Before correcting an operator, substitute its own definitions and read what it delivers |
| **TRAPS: the-pin-had-a-fixed-point** | the mass rescale's `k = 1/(φ_msg + R_own)` has a fixed point at `R_own = ½` and drove the delivered gDNA fraction to ½ regardless of truth. ⚠ It hid from every aggregate — the per-step rescales telescope back to 1 at a gene's far end, so only a per-object check away from a pure-gDNA object sees it |
| **TRAPS: no-belief-not-no-numbers** | gating the rescale wherever it read the destination's own density looked clean and broke the capture landscape — the off-probe floor leaked into every exon. ⭐ State the licence as "no BELIEF may enter", not "the destination's numbers may not enter" |
| **TRAPS: a-licence-with-no-floor** | κ̂ = 0.500689 on a genuinely unstranded library leaves `I(f_g) ∝ (2κ−1)²` at ~5e-7, and a licence that tests a PRECISION with no floor is granted by it. The repair is to propagate `Var(κ̂)`, not to add a floor |
| **TRAPS: a-multiplication-gated-by-a-trace** | an intron's 2.5 % phantom RNA is the only nonzero RNA precision in a chain and alone unlocks a reframe compounding a gDNA level to 2.16×. ⭐ A predicate gating a MULTIPLICATION must be sized by how much density stands behind it, not by whether a precision is strictly positive |
| **TRAPS: all-small-singly-large-jointly** | removing each ψ channel in turn moved the worst object by ≤0.016 of a 0.217 error while removing all of them moved it the whole way, because all three are built from the same relayed level. ⛔ All-small-singly plus large-jointly means stop ablating consumers and go one stage upstream |
| **TRAPS: recompute-from-the-oracle** | an impossible-looking enrichment ratio of 1.46 with capture OFF was correct: `ρ_tot` is a total and the destination held 55 RNA fragments where the source held 0. ⭐ Recompute the quantity from the ORACLE before assuming the formula is broken |
| **TRAPS: a-cancelling-defect-pair** | ⛔⛔ fixing one of two errors that CANCEL is worse than fixing neither — correcting one hop alone moved a toy's evidence-free exon 0.0107 → 0.0244 while the rung it targeted improved 26 %. Price such a fix in the arm that also removes the other defect; a cancelling pair is one experiment, not two |
| **TRAPS: a-variance-cap-asserts-certainty** | capping `Var(f_g)` at `f_g(1−f_g)` asserts certainty at the corner, and an evidence-free slot's `f_g = 1` init parks there, so the composition half of `Var(log ρ_tot)` is exactly 0 wherever there is no evidence. Replacing it with the reference prior's own `Beta(½,½) ⇒ ⅛` is derivable and inert (−0.3 %); any repair routed through `composition_logvar` is bounded by this |

**zero-the-precision-with-the-value. ⛔⛔ A REFUSED CLAIM MUST LOSE ITS PRECISION IN THE SAME
STATEMENT THAT
ZEROES ITS VALUE.** A value zeroed at one line and a precision handed back at a later one is not "no
claim",
it is the CONFIDENT ZERO the licence exists to forbid, and the relay's mass rescale converts it into
its
opposite: a delivered `n = 0 @ pn = 280` became `k = 467,000×`, "all your mass is gDNA", one hop on
in a
26 k-fragment exon whose truth was 0.000.

⭐ *The rule:* every operator that GRANTS a precision — SPLICE IN, premise variance, a conservation
term —
must run before the licence's zeroing, in both twins, so the zeroing is the last word; and the
falsification
fixture must carry that operator's input, or the assertion passes on a path it never exercised. ⛔
Score
such a mechanism where the damage LANDS, one hop past the confident zero, not where it is written —
a census
at the written site read 2.6 % of the error. *Sibling:* `a-licence-with-no-floor`.

**no-prior-means-haldane. "No prior" does not exist on a grid — omitting a term lets the grid supply
Haldane** (`p(x) ∝ 1/x`, improper, an amplifier toward the vertices). Posterior median spread over
grid
half-widths 4–20 is 0.045 at Jeffreys and 0.525 at Haldane.

**prefer-shares-to-differences. Sums are well conditioned; differences are not.** Subtracting across
a sj
gives `Var(log ρ) = u²σ_T² + (u−1)²σ_μ²` with `u = 1/(continuing share)` — at the real median u =
2.3,
already the boundary of validity, and at p75 u = 5.3, hopeless at any depth. Prefer shares.

**an-all-zero-factor-is-inert. An all-zero factor is uninformative, not decisive.** In a
multiplicative
score, a factor that is zero for every candidate annihilates the other factors and collapses the
record to a
coin toss. Skip a flat-zero factor; do not multiply by it.

**density-below-one-fragment-length. Density below one fragment length is not resolvable by ANY
design.** A
1 bp region has no independently measurable density and never will, though *composition* still does,
since
it depends on what the fragments are rather than where. *Corollaries:* an object with zero
opportunity for a
component must emit nothing at zero precision, not a floored division; and "no data" must be inert,
never
"100 % gDNA", which was actively seeding false gDNA into neighbouring exons.

**identical-paralogs-are-bimodal. An identical-paralog split is bimodal, and depth does not fix
it.** Same
code and seeds, `n_fragments` 3,000 → 171/0 but 6,000 → 249/237: which branch it lands on is a draw.
⛔ Do
not move a seed or a depth until it lands even — that is tuning to green. Assert the collapse
STRUCTURALLY
(`min == 0`, `max == total`), which still fires the day the tie breaks and additionally catches a
partial
break, and score only what is identifiable. ⚠ Two defects live here, not one: the total is also
over-called
by 17 %, so "score what is identifiable" must not be read as "the rest is fine". ⭐ And it is not
noise —
with fractional assignment and only `iterations` varied the EM converges to the vertex on a real
gradient
(5 → 93.6/81.2, 20+ → 174.5/0.0), so a tilt is entering from outside the RNA likelihood terms.

**a-mean-hits-the-mass-weighted-centre-by-luck. A scalar that "wins" on a mass-weighted metric may
only be
sitting on the mass-weighted centre.** Measured four times in one investigation — a library-wide
reference
mean, an object-weighted mean, a pooled RNA density, an exonic class constant — because each
population was
BIMODAL at `f_g ≈ 0` and `f_g ≈ 1` and `Σ|Δ|·M` weights by mass, which sat at one mode.

⛔ *The tell:* substitute the candidate by its own class mean. If that costs nothing (168,551 vs
164,074
fragments at `R exon`) it IS a constant, and its sd against the truth's says so directly — 0.0021
against
0.4441. Report both beside any per-object claim.

**a-clamp-at-the-closed-end-escapes-the-window. A prior clamped at the closed end of its support can
put
almost all its mass outside the solve window, and a sweep over the interior will never see it.** ψ's
reference was clamped at `m = 1 − 1e-8`, leaving only 0.94 % of its mass inside the shipped `L = 10`
window,
while the `L`-invariance gate swept `m ∈ [0.01, 0.99]` and the arm operated only at the clamp.

⭐ *The rule:* test a prior at the value the code will actually use, and derive its floor AT THE
OBJECT —
one pseudo-fragment here gives `m_i = E[g]_i/(E[g]_i + 1)`, which keeps 0.909–0.990 of the mass
inside the
window where a pooled floor keeps 0.04.

**the-deconvolution-is-as-good-as-the-density-it-is-handed. A residual estimator inherits the error
of the
density it subtracts, amplified.** `f_g = ρ_g·E_g/M` needs no RNA model — RNA is whatever the gDNA
deconvolve leaves — which is what makes it usable pre-solve, but its accuracy is exactly `ρ_g`'s: at
capture-OFF the density is 1.00× of truth and the deconvolution scores 0.026, under capture the
density is
0.28–0.38× and the same deconvolution scores 1.091, worse than the uninformative reference and a
5.3× panel
loss on stranded × capture-ON. ⭐ Price the DENSITY against truth before believing anything built on
it, and
state the verdict per stratum of that density rather than pooled.

**deriving-one-coordinate-propagates-its-error. ⛔⛔ CLOSING A COMPOSITION BY DERIVING ONE COORDINATE
FROM
ANOTHER MAKES THE DERIVED ONE INHERIT EVERY DEFECT OF THE SOURCE, INCLUDING ONES THE OLD ESTIMATOR
DID NOT
HAVE.** Defining the RNA total as `1 − f_g` removed a real defect, but `f_g` was the grid-SNAPPED
median
where the retired read-out had been a posterior mean that is exactly ½ on an evidence-free object,
so an
exact number was silently traded for a snapped one (0.5 → 0.51256).

⭐ *The rule:* when you make coordinate B a function of A, enumerate what B used to be exact about
and check
A is exact about it too; "at least as good on average" is not the test. ⚠ The repair was to fix the
SOURCE
— a continuous quantile — not to undo the derivation. *Sibling:*
`interpolate-on-the-axis-where-the-lattice-is-uniform`.

**interpolate-on-the-axis-where-the-lattice-is-uniform. ⛔⛔ A QUANTILE INTERPOLATED ON A NON-UNIFORM
TRANSFORM OF THE GRID IS BIASED TOWARD THE MIDDLE, AND IT LOOKS RIGHT.** ψ's ½-quantile was
interpolated in
`f_g`-space, where the σ lattice spacing runs ~1e-5 at the ends and ~0.085 in the middle, so a
posterior
concentrated on one grid point came back biased toward ½ by 2.71e-03; interpolating on `λ`, where
the
lattice IS uniform, returns that posterior's own grid point to 2.2e-16.

⭐ *The rule:* do quantile arithmetic on the axis the grid was BUILT on, then map. ⚠ The case is not
synthetic — an unsolved slot's fed-back belief is a one-hot posterior, landing exactly on the bias.

**read-the-whole-failure-list. ⛔⛔ 21 FAILURES WERE REPORTED AS "21 GOLDEN FAILURES" FROM THE LAST
ELEVEN
LINES OF THE OUTPUT, AND TWO WERE REAL** — one a broken caller, the other the most informative
defect of the
session, nearly filed as expected churn.

⭐ *The rule:* derive the failure set, never eyeball the tail —
`pytest -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c` prints one line per FILE and
fits on a
screen at any failure count.

**a-priors-curvature-is-not-the-datas-information. ⛔⛔⛔ A PRIOR MAY NOT CONTRIBUTE TO A FISHER
PRECISION.**
`region_init`'s `tau_lam` is the DATA's information on the composition axis; letting a prior's
location
curvature into it was built and refused on three independent counts — the observed 3,227× fall is
the
JACOBIAN (`[f(1−f)]²` alone predicts 3,154× of it, so the likelihood genuinely is that flat there),
the term
acts as a boolean gate flip rather than a contribution (τ = 0.029 and τ = 1e6 both return 850.44
against a
ceiling of 850.50), and it credits DATA-FREE slots (`n = 0` goes `prec_g` 0 → 0.2026).

⭐ *The rule:* before adding any term to a precision, ask whether it scales with the DATA. If it does
not,
it is a prior, and its place is in the posterior. *Sibling:* `a-four-decimal-print-is-not-a-zero`.

**a-refutability-test-needs-the-refuting-channel-in-the-fixture. ⛔⛔⛔ A PRIOR MEASURED WITHOUT THE
EVIDENCE
THAT COULD OVERTURN IT IS NOT BEING TESTED — IT IS BEING ASKED TO ANSWER ALONE.** ψ's structural
prior read
2.1–5.7× worse on a toy chain with no intergenic regions: `fit_intron_background` pools intergenic
only, so
it returned uninformative and the intron-vs-intergenic mechanism was never in the room, and at κ = ½
the
strand channel is dead by derivation too, so the fixture had no refutation channel at all. Rebuilt
with
intergenic flanks, which production always has, the same prior yields to the same evidence to within
0.02 of
the no-prior answer.

⭐ *The rule:* before measuring whether a prior can be overturned, ASSERT that the overturning
channel is
live in the fixture — `bg.informative`, a nonzero `τ`, a fired counter. ⚠ The measurement itself was
taken
on a nascent-heavy toy, so its magnitudes are a STRESS reading and not a design driver; the method
is
universal.

**a-strength-is-a-nat-a-prior-weight-is-a-count. ⛔⛔ A STRENGTH IN NATS AND A PRIOR WEIGHT IN
PSEUDO-COUNTS
ARE NOT THE SAME UNIT, AND EQUATING THEM IS AN ANALOGY WEARING A DERIVATION'S CLOTHES.** Setting
`m = σ(a+b)` reads as principled and is dimensionally wrong; the derivation that stays in one
currency asks
what mean one pseudo-observation would produce — `Beta(a,b) → Beta(a+1,b)`, mean 0.75 at `a = b =
½`.

⛔ The sharper half is what it replaced: taking the strength from the LATTICE (`m = σ(L)` ⇒ 9.31
nats) is
not a statement about the claim at all, and measured against being refuted it was worse than having
no prior
(2.0247 vs 0.3946). A prior's strength must come from what it CLAIMS, never from what the
representation can
hold — a cap is not a choice. ⚠ Check the optimum is BROAD before trusting any derived value.

**a-four-decimal-print-is-not-a-zero. ⛔⛔ A DIAGNOSTIC PRINTED AT FOUR DECIMALS READ `0.0000`, AND A
SESSION'S WHOLE MECHANISM WAS BUILT ON IT BEING ZERO.** The value was 4.21e-05, the predicate the
diagnosis
turned on returns True on it, and the real effect was a 3,227× change of MAGNITUDE that needed a
different
repair and, as it turned out, no repair at all.

⭐ *The rule:* a number a diagnosis TURNS ON is `repr`'d, never formatted, and "exactly zero" is
asserted
with `== 0.0` rather than read off a table. ⚠ The cost was a falsification gate written against a
non-defect: it was verified failing, and then correct code was changed to make it pass.

**a-constant-in-exact-arithmetic-is-not-constant-in-float64. ⛔⛔ "THIS TERM IS CONSTANT SO IT
CANCELS" IS A
CLAIM ABOUT ℝ, AND THE GRID DOES NOT LIVE THERE.** ψ's location term at its neutral `m = ½` is `log
2`
exactly — the whole reduction property — but `logaddexp` leaves that row with `ptp = 2.22e-16` at
every grid
size, and since the median reads the first grid point whose CDF reaches ½, a balanced AMBIG slot at
`κ = ½`
moved `f_g` 0.5423 → 0.4577, one full grid step, at ~37,000 slots the prior claimed to say nothing
about.

⭐ *The rule:* where a term is documented as vanishing, RETURN the vanishing value rather than
computing
something that equals it, and gate the reduction on the ANSWER (`solve(neutral) == solve(None)`, bit
for
bit), never on `allclose` of the term — an `atol = 1e-15` gate on the term passed throughout.

**a-toy-and-a-panel-can-disagree-in-rank. ⛔⛔⛔ A CHANGE MEASURED ON A TOY CAN NOT ONLY SHRINK ON THE
PANEL
BUT INVERT ITS ORDER** — the best arm of three on the test chromosome was the WORST of three on all
three
in-scope strata of the 16-condition ladder, same metric, same scorer, same axis.

⭐ *Why it inverts, which is the part worth carrying:* a toy's objects mostly have little or no own
evidence, so a message is nearly free information, while the panel's mass sits on objects that can
already
answer for themselves, where the same message competes with a measurement. Anything that trades
accuracy
from evidence-bearing to evidence-free objects therefore wins on the toy and loses on the panel by
construction.

⭐ *The rule:* a toy result is a MECHANISM check, never a RANKING — run the panel before any "this
beats
that", name the substrate in the same sentence as the number, and split the error by whether the
destination
had its own evidence to predict the inversion cheaply. *Sibling:* `panel-before-src`.

**a-rescale-that-reads-the-source-belief-is-unbounded. ⛔⛔⛔ A TRANSPORT WHOSE FACTOR IS DERIVED FROM
THE
SOURCE'S OWN CLAIM AMPLIFIES A WEAK CLAIM WITHOUT LIMIT, AND THE WEAKER THE SOURCE THE LARGER THE
FACTOR.**
The mass-identity rescale `k = M_dst / Σ_c ρ_c,src·E_c,dst` is exact under its premise, reproduces
the
worked numbers, and never touches the destination's belief — every property one would check — and it
still
handed a 23,889-fragment exon "all your mass is gDNA" at `k = 235,800` on a zero-gDNA library,
because the
source's claim is its DENOMINATOR.

⭐ *The rule:* transport by a quantity MEASURED at both ends — here the model-free abundance ratio —
never
by one the source's belief appears inside. The two agree exactly where the belief is right, which is
what
the worked numbers check, and diverge exactly where it is not. *Sibling:*
`zero-the-precision-with-the-value`.

**a-face-total-is-not-a-total-without-its-flux. ⛔⛔ COMPARING TWO OBJECTS' TOTALS WHEN ONE OF THEM
CANNOT
HOLD PART OF THE POPULATION REPORTS A DIFFERENCE THAT IS STRUCTURE, NOT ENRICHMENT.** Mature RNA
cannot
cross an `exon|intron` boundary contiguously, so it appears there as sj FLUX rather than as
crossings, and
an enrichment ratio built from the boundary's unspliced abundance against the exon's total read a
30.8×
"depletion" at a capture-OFF condition with no probes in it at all.

⭐ *The rule:* a face's total is its contained/crossing abundance PLUS the sj flux whose bodies lie
on that
side. ⚠ And that flux belongs to ONE face — pooling it over-counts the boundary against both flanks.

**an-imputation-must-cost-something-every-hop. ⛔⛔⛔ IF EVERY VARIANCE IN A MESSAGE LAYER SHRINKS WITH
COUNTS, THEN BETWEEN TWO DEEPLY-COUNTED SLOTS AN IMPUTATION CROSSES FOR FREE AND ARRIVES AT FULL
STRENGTH
BESIDE A REAL MEASUREMENT.** Owner: message propagation is an imputation and a weak predictor by
definition,
the strand model is a measurement, so propagation must be weak.

⭐ *The rule:* the PREMISE of a hop — "my neighbour's values apply here" — has a variance of its own
that
does NOT shrink with depth of counts, and it must be charged on every hop under every strategy. Fit
it by
method of moments so no constant enters: `premise = max(0, Var(log r) − mean(v_r))`, the observed
spread of
the hop ratios less the counting variance that spread already contains, which fits exactly 0 on a
substrate
whose ratios vary no more than Poisson predicts.

⭐ Measured by splitting each condition's error on whether the destination had its OWN composition
evidence, the policy was degrading the MEASURED half by 8.09× while leaving the imputed half neutral
at
0.98×; after the premise term, 0.81×. ⚠ The honest cost, and it is evidence the fix is real rather
than a
wash: unstranded × capture-ON — the one stratum with no strand channel — gets worse. ⭐ Split any
message
layer's error this way: "the messages help" and "the messages trample a measurement" are different
findings
and a pooled number cannot tell them apart.

**a-floored-knob-is-not-the-bandwidth. ⭐⭐⭐ A SMOOTHING CONSTANT THAT IS CLAMPED FOR 99 % OF THE DATA
IS NOT
THE BANDWIDTH, AND SWEEPING IT ANSWERS NOTHING.** `landscape.knn_widths` returns
`max(scale · d_k, grid_step)`, and on the ladder's ~30,700 regions 98.9 % of kernels sit at the
floor at the
shipped `_KNN_SCALE = 0.5` (92 % at scale 2.0), so the smoothing actually in force is the GRID STEP,
`span/(_N_GRID − 1)` ≈ 0.034 decades.

⭐ *The tell:* a sweep over a smoothing parameter returns the same score to three decimals at several
adjacent settings — that is a disconnected knob, not a flat optimum. Measure the share of kernels at
the
floor before reading any bandwidth sweep. ⚠ It is not a bug: the kernel asks for finer resolution
than the
render can represent and is correctly refused, which is why the same constant is live on a 70-region
toy and
dead on the panel.

⛔ The constant that IS in force was then selected honestly: over a 16× range the held-out predictive
likelihood's knee, the split-half mode reproducibility maximum (0.976, falling to 0.889) and the
anchor-gap
collapse (0.126 → 0.003 nats) all land on the shipped `_N_GRID = 260`.

**a-mode-count-is-not-a-well-posed-quantity. ⭐⭐⭐ HOW MANY MODES A FITTED DENSITY HAS IS A PROPERTY
OF THE
RENDER RESOLUTION, NOT OF THE FIELD — SO NO CONSUMER MAY DEPEND ON IT.** Across `_N_GRID` 65 → 1040
the
`AbundanceLandscape`'s mode count runs 3.7 → 17.7 and never converges, it tracks 1/step, it buys
nothing
(the held-out likelihood is flat over that range), and the modes reproduce worse across split halves
as the
count grows. So "the fit looks spiky" is answered with neither "signal" nor "wrong bandwidth": the
extra
maxima are the axis's own resolution, and a wiggle owning no basin mass is harmless.

⛔⛔ *And sweep the resolution before a consumer reads any SHAPE statistic off a fitted density.*
`rho_0`
and the anchor-consistency verdict survive, but `span_R` reads 58 → 95.6 → 1.9 on one row, because
`split_basins` picks the enriched mode by BASIN MASS and over-resolution fragments the bulk into
competitors. Keep only what does not move: a location and a containment verdict did, a ratio of two
selected
modes did not. *Sibling:* `a-floored-knob-is-not-the-bandwidth`.

**measure-a-default-flip-before-you-write-it. ⭐⭐⭐ A CONFIG DEFAULT AND A REFUSAL ARE ONE DESIGN, AND
FLIPPING THE FIRST CAN INVALIDATE THE SECOND.** `CalibrationConfig.abundance_landscape` was written
opt-in
and OFF and — correctly for an opt-in flag — refused to run without the wall arrays; flipping the
default
broke 65 callers (24 failed + 41 errors) on a single cause, that unit and toy fixtures have no wall
arrays
and never wanted the object at all.

⭐ *The tell:* you are about to change a default on a flag whose other branch raises. The refusal was
written for the population that opts IN; the default serves the population that never thought about
it, and
those are different populations with different correct answers.

⭐⭐ Re-read what the refusal was protecting against — it is usually narrower than the raise implies.
Here
the stated reason was "refusing rather than fitting on unmasked totals", and the alternative to
refusing was
never "fit unmasked", it is NOT TO FIT, so a missing input now skips LOUDLY and yields a `None` the
downstream reader already tolerates. ⛔ Keep the asymmetry explicit with a gate:
`background_abundance`
keeps its refusal because that pair feeds ψ, and the landscape's gate asserts that the first one
still
refuses.

---

## E. Structure, indexes and plumbing

**one-reference-hides-refid-bugs. Single-reference synthetic indexes hide reference-id-space
mismatches.** A resolver assigning
ref-ids by first-seen order silently dropped **476,719 of 476,732** real fragments while every
golden test passed.
The fix is per-path: RNA-only spike-ins do nothing for the gDNA deposit path — that needs at least
two genomic references.

**annotated-is-not-genomic. "Has an annotation" is not "is genomic" — a classification must be an
INPUT, not a proxy.** Choosing
gDNA references as `{t.ref for t in transcripts}` filled the panel with gDNA on the RNA-only
spike-ins whose zero truth
is the false-positive control. A mis-stated classification must raise, not silently produce nothing.

**an-sj-is-not-a-gap. A splice junction CANNOT be detected from a gap between deposited slices.** A
contiguous spliced
read whose exon body straddles an internal region bound has no gap at all. The sj's identity is the
cut-intron
coordinates, which the scanner already has — pass them through.

**deposit-at-the-sj. A splice deposit belongs at the SJ's coordinate, not the region's boundary.**
Invisible for
annotated introns (their ends are region bounds); for unannotated sj the mass lands kilobases away.

**splicing-makes-the-graph-cyclic. Alternative splicing makes a region↔sj graph CYCLIC, and cycles
are the common case** — the
human graph has ~404,000 independent cycles, one per sj, and two-sweep forward-backward is exact
only on a tree.
Never break a cycle by dropping a sj boundary — that re-isolates the exon the boundary exists for.

**nrna-does-not-mean-synthetic. `~is_synthetic & ~is_nrna` as the "real transcript" filter deleted
52,104 real termini** — on a
non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent", not "manufactured span". ONE
filter: `~is_synthetic`.

**credit-exactly-one-sj. A fragment crossing K sj must credit exactly ONE** (the leftmost
annotated). Crediting all K
shifts the library sense fraction 21–34 % and creates between-side correlation that reads as
overdispersion (a zero-overdispersion simulator fit 0.092). A `1/K` split is provably biased by 4–12
σ.

**strand-completes-the-sj-key. `(src, kind, dst)` is not a total order for sj boundaries** — two
strand-coincident sj differ
only in strand, so ordering becomes input-order-dependent and the duplicate check reads them as
duplicates.
GENCODE has zero of them, so only a synthetic stress test finds it. Sort on `(src, kind, dst,
strand)`.

**a-hash-that-misses-its-artifact. A cache key must cover the artifact it caches, in all three
forms.** ① A derived hash stored
beside its data verified a stale cache CLEAN when a fix rewrote files the hash did not cover —
compute hashes on demand, never store them.
② An exclusion from a key is a CLAIM that nothing the excluded file produces is stored in the blob;
when two of its products were,
the cache served one fresh arm beside stale ones with the key reading ok — read the exclusion list
as a sentence, and prefer storing
INPUTS and deriving at read time. ③ A dataclass field computed at construction survives
`dataclasses.replace` and then describes the
arrays it replaced: if a value's inputs are in somebody's override list, it must be a `@property`,
never a stored field.

**integer-channels-reproduce. A COUNT is an integer and reproduces bit-identically across worker
counts; a FRACTION is float64
and does not** — one numeric convention, with tests validating within a DERIVED tolerance bracketed
from both sides. The fixed-point
alternative was LESS accurate than the float it defended against: against exact rationals it missed
by up to 2.0e-7 where float64
missed by 2.8e-13 — **714,000× worse**. Corollary: the one bank whose ORDER is observable must be
sorted on its own content before crossing the ABI.

**worktrees-run-the-wrong-code. Worktrees silently run the wrong code** — an editable install's
meta-path finder beats
`PYTHONPATH`, so an A/B inside a git worktree executes the main repo's source.

**checkout-deletes-uncommitted-work. `git checkout -- <file>` does not undo a perturbation when the
work is uncommitted — it deletes the
work.** A perturbation harness must restore from a copy of the WORKING TREE. Cost one full
re-implementation.

**two-masks-one-name. Two different masks shared the word `struct_lock`, and BOTH were right** — one
meant "pinned and certain",
the other "may emit composition certainty". Two correct predicates under one name is worse than
either being wrong:
give each its own name and ONE home, and let every consumer import it.

**two-docstrings-one-quantity. The prose next to the code said "the AVERAGE" and the code followed
the prose,** while a sibling
module's docstring had the correct formula the whole time. Two docstrings disagreed about one
quantity for months and nobody diffed them.

**a-transcript-predicate-must-not-silently-drop-a-molecule. A fragment rejected by a
transcript-level predicate deposits NOTHING,
and the rejected population is never random — the bias concentrates exactly where the predicate
fires.** A chimera check that joined
mates only when their transcript sets intersect dropped **4,087 gDNA fragments per condition — 0.04
% of fragments carrying 2.4 % of
every boundary crossing** — and hid because every dropped fragment was a crosser, so the
contained-mass gate still closed.
Check genomic COMPATIBILITY (one reference, mates inward, plausible implied length) before calling a
fragment a chimera; gDNA
routinely spans unrelated transcripts. The transferable rule: make the fragment ledger CLOSE —
`n_fragments == deposited + Σ accounted drops`; a drop with no counter is invisible to every
downstream check.

**an-object-class-does-not-see-a-terminus. The seven object strata classify a boundary by the
exon-ness of its two flanks, so a
TSS/TES inside another transcript's intron is indistinguishable from a splice site there — and the
two have OPPOSITE currencies**:
at a terminus the RNA originates and a composition cannot cross; at a splice site it enters by the
sj and the composition does.
Measured: 209 k of a 221 k composition error into exons sat on 3,867 terminus boundaries while the
15,480 splice-site boundaries
were near-exact. A hop type is `object class x {sj, term, sj+term}`, read off
`RegionStatics.boundary_flags`; the classes themselves
remain right for ψ's reference — it is the HOP that needs the second bit.
*Sibling:* `two-masks-one-name`.

---

## F. Domain facts that read like defects

**specificity-and-sense-are-complements. Strand specificity is TWO different quantities and they are
complements** — a simulator's `strand_specificity` is direction-agnostic protocol fidelity; a fitted
sense fraction is directional, so on an R1-antisense (dUTP) protocol a fitted κ of 0.0101 on a "0.99
stranded" library is correct, and forcing 0.99 measured 166× worse. The matching quantity is
`StrandModel.strand_specificity`; the recovered values live in the test that measures them, not
here.

**strand-measures-the-tilt. Strand measures the TILT, not the gDNA fraction.** With RNA tilt `d = f₊
− f₋`, `p = ½ + (κ−½)·d` — the gDNA fraction cancels identically, and the information about it is
exactly zero at κ = ½ for any count. Strand reaches gDNA only through the triangle bound `f_g ≤ 1 −
|d|`.

**a-linear-likelihood-emits-a-sign. A likelihood that is asymptotically linear in its parameter has
no mode, only a direction — and on a bounded grid that direction saturates at an endpoint.** The
tell is free: expand the log-likelihood in the parameter; if the leading term is linear, the argmax
is a sign and the object is a vote, not an estimate — summarise it as a location plus an
information, never a row handed to a normaliser. The channel that taught this reported gDNA
fractions of 0.72 / 0.59 / 0.72 against truths 0.00 / 0.00 / 0.57 when its two pmfs were ~1e-9 bp
apart, and widening the gap improved every one — so "fix the pmf" repairs are refused. That channel
is purged from `src/` and is not a candidate; the lesson is general.

**amplitude-fades-influence-does-not. A term that normalises away can still decide the answer,
because an argmax is scale-free.** A term entering a normalised sum has three quantities —
amplitude, participation, declared precision — and they do not fade together: at a 5e-12 bp pmf gap
the amplitude was ~1e-11 while participation was 100 % of library mass and the declared precision
was the grid's own width sold as evidence. Report all three separately. *Sibling:*
`a-linear-likelihood-emits-a-sign` (same campaign).

**a-pooled-conversion-applied-per-component. A ratio measured on the pooled population and applied
to each component separately is population-blind — and the blindness need not be the axis you
expect.** The prior's crossing→fragment conversion `q = mass/count` errs worst at EQUAL lengths
(+3.18 % vs +0.63 % at a 17× larger gap): the driver is placement — gDNA crosses boundaries in long
intergenic regions where `q → 1`, RNA in short exons. Bounded at ≤ 0.6 pp of composition, at
calibration's noise floor: record the bound, build nothing.

**an-identity-with-a-qualifier. An identity that holds "over X" is a measurement waiting to be made
— price the complement, or the qualifier is a hole nobody has sized.** The conserved-mass identity's
accurate qualifier "over unspliced fragments" hid a complement of 1,222,375 of 4,830,713 RNA
fragments (25.3 % of the RNA library) depositing on no conserved bank at all. The tell: grep your
invariants for "over", "for", "assuming", "where" — each names a population with a size nobody has
measured.

**equal-lengths-carry-no-composition. At equal component mean lengths the length channel carries
exactly zero composition information, at any depth** — the deconvolution is identified only through
`μ_g − μ_r`. This fact is the gDNA ladder's forcing function, not a caveat about it: zero length
information is what stops the EM answering without exercising calibration.

**capture-is-1000x-on-exons. Hybrid capture is ~1000× on exons and only gDNA reads it cleanly** —
RNA's own 10⁴ expression range hides the probe pattern; gDNA's uniform baseline does not. Capture
also destroys the intron signal ~75×, so nascent-vs-gDNA is fundamentally unidentifiable under
capture — one more reason the nascent scope ruling puts dedicated nascent-capture experiments out of
scope.

**capture-selects-for-length. Capture selects for length, and the post-capture distributions are the
baseline.** A short fragment presents less sequence and is captured less efficiently, so the
pre-capture parameters describe a library that was never sequenced — score against post-capture
truth, for lengths exactly as for abundances. Capture narrows the gDNA↔RNA length gap whenever gDNA
is the shorter component.

**on-target-by-start-is-geometry. "On-target" defined by the START's territory is geometry, not
capture efficiency.** Conditioned on capture, an intronic-start fragment was long enough to reach
the probe (weight ~`w²/2`) while an exonic start is flat in `w`, so intronic-start populations read
longer under any capture model. The population that physically binds is the one that overlaps a
probe.

**real-data-is-a-test-input. Real data is a TEST input, never a DESIGN input.** The cfRNA on disk is
one far end of the RNA-seq spectrum, not a sample of it: sweep the plausible space, report the worst
case, bring the domain call to the owner. In particular, never assume RNA fragments are longer than
gDNA — true for cfRNA, false elsewhere.

**eff-lengths-do-not-cancel-at-an-end. "Effective lengths cancel, so a region's marginal is just its
length" is false near any transcript end** — a mature fragment must fit in the remaining transcript;
gDNA need not.

**configured-lengths-are-not-realised. Equal CONFIGURED fragment lengths do not give equal REALISED
ones** — each pool's length marginal is reweighted by its own template opportunity, and a 2 kb
transcript truncates the tail a whole chromosome does not. Gate any length axis on the realised
truth, never the config.

**mature-rna-never-crosses-a-boundary. Mature RNA never crosses an exon↔intron boundary** (0 of
1,146 boundaries, 7/7 conditions); exon↔exon boundaries it does cross. This is the hard empirical
case that a contiguous boundary and a splice junction are physically different objects, and it is
what makes the two exon-crossing gDNA pools pure.

**a-boundary-with-rna-is-not-an-sj. But a boundary with RNA crossing it need not be a splice
junction.** One position can be a splice donor for transcript A and plain contiguous exon for
transcript B; zero-gDNA libraries show boundaries carrying 44–55 k unspliced fragments that are 100
% RNA.

**the-panel-enriches-nascent-by-its-own-probes. A simulator that models a population in its own
private space gives it its own private physics** — the tool is then debugged against a library that
does not exist. The defect is FIXED (the index-driven simulator, 2026-08-19: nascent RNA is an
ordinary transcript in one multinomial, and gDNA and nascent under one probe now enrich at the same
rate, 1162× vs 1003×); the general rule survives: a simulated population must be a first-class
member of the same object set, weighted the same way and reached by the same mechanisms. Under the
nascent scope ruling, nascent-under-capture fidelity is doubly out of scope — dedicated
nascent-capture experiments are excluded, and any number at the panel's nascent share (20.2 % capture-OFF and 2.6 % capture-ON) is a
stress reading. *Sibling:* `real-data-is-a-test-input`.

---

---

## G. Process

**no-magic-numbers. No magic numbers.** Stop and discuss before adding any constant, heuristic or
tunable. Every
divisor must be derived from the deposit rule and unit-tested against brute-force enumeration.

**one-thing-varied. One thing varied per experiment**, with the falsification test written first and
verified failing,
and a baseline re-recorded from the current tree in the same session.

**converge-and-delete. Converge and delete.** No legacy, no backwards compatibility, no speculative
code. Code kept "for
comparison with the old version" is a defect. No version suffixes in file names.

**the-source-does-not-cite-docs. The source does not reference the docs.** Docs evolve and rot — 73
% of the citations that used to
be in the source pointed at documents that had already been deleted. A docstring may cite a **test**
or
an executable specification, because code cannot rot silently; it may not cite a document.

**running-an-arm-is-a-fresh-process. Six operational traps from running panel arms, each of which
cost a launch.**
(i) Never edit `src/` while an arm is running — every arm imports `src` at start-up, so a mid-flight
edit silently changes what half the shards measured; adding a NEW arm to `scripts/` is safe.
(ii) zsh does not word-split an unquoted variable — use an array and `"${CONDS[@]}"`.
(iii) A wait-loop whose `pgrep` pattern matches its own wrapper deadlocks — wait on a log marker
instead.
(iv) A default suite path stored in two homes goes stale in one of them unnoticed — three
instruments once died on "Failed to open BAM" because one copy still pointed at a deleted panel.
(v) Node-axis and region+boundary figures differ by ~2× — say which one, every time.
(vi) A composite arm fires only its COMPONENTS' names, so a guard keyed on the arm's own name trips
after a complete, valid run — check WHY a guard fired before distrusting the data. *Sibling:*
an-ablation-that-never-ran.

**shard-an-arm-sweep-by-condition. Shard a panel sweep by CONDITION, never by ARM over one condition
— the instruments WRITE their caches.** Concurrent arms on the same condition race writers on one
directory: one `payload.npz` came out truncated and every later read raised `BadZipFile`. The damage
is loud, not silent — digests and sum-to-full gates kill the run rather than feed a wrong truth; an
instrument that only READS should take its full payload from `scan_cache/`, not from a directory
another instrument writes.

**no-enumeration-without-a-census. Do not re-propose path or cell enumeration without a memory
census.** Possible unspliced paths ≈
1 M regions × 3–6 reachable ends at ~100 B each = 0.3–0.6 GB, plus spliced paths. It was killed by
memory,
and no consumer needs it.
