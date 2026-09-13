# TRAPS — mistakes this project has already made

Read before designing anything. Every rule below was paid for: something was built, measured, and came
out wrong in a way that was not obvious in advance. A rule is a lesson, stated as what looked right, what
was wrong with it, and what to do instead; it is not a measurement, and a number that says where the tool
*is* lives in `ROADMAP.md`, re-derived by the instrument it names. Open problems and refusals live in
`ISSUES.md`, rulings in `DESIGN.md`, derivations in `EQUATIONS.md`.

Cite a rule by its name, e.g. `TRAPS: panel-before-src`. The name is the identifier, so a citation says
what it means without a lookup and is one greppable string with one home. `tests/test_no_jargon_labels.py`
enforces it and refuses the numbered labels (`A16`, `D4j`, `C0b`) these rules used to carry.

The shape is load-bearing: a rule starts at column 0 as `**name.` followed by one paragraph, and nothing
else in this file may start a line with bold, because a bold token at column 0 manufactures a phantom rule.
A rule names a defect, never a licence: citing one authorises no mechanism and vetoes none. Rules measured
on the retired uniform-nascent panel, or at the current panel's stress-level nascent share, say what is
robust rather than what is likely.

## THE INDEX — every rule, by section

- **A. Validation and gates** — `self-checking-validator` · `perturb-every-gate` · `a-field-driven-gate-is-atomic` · `waive-with-a-measurement` · `a-docstring-that-misdescribes-the-graph` · `a-flat-pile-is-not-a-knot` · `the-rename-that-corrupted-a-diagram` · `a-constant-parked-a-value-off-a-knife-edge` · `a-gate-that-restates-the-implementation` · `a-gate-that-already-passed` · `right-conditional-wrong-marginal` · `byte-identity-gate` · `the-deliverable-is-not-reproducible-by-default` · `a-clip-hides-a-scale-error` · `an-ablation-that-never-ran` · `a-green-suite-hid-five-dead-instruments` · `compatibility-is-geometry-not-composition` · `a-zero-count-is-a-measurement` · `a-ratio-cannot-carry-zero` · `the-divergence-was-a-barrier` · `deadband-from-the-wrong-sample` · `honesty-metrics-reward-ignorance` · `predicate-contradicts-its-docstring` · `a-test-that-redefines` · `a-gates-power-is-its-invariant-set` · `a-gate-that-reconstructs` · `off-grid-message-mode` · `a-comment-quoted-as-a-finding` · `the-intermediate-is-not-the-deliverable` · `could-the-arm-have-fired` · `prove-the-substrate` · `can-the-benchmark-resolve-it` · `toys-rank-hotspots-backwards` · `an-identity-with-a-qualifier`

- **B. Measurement and inference** — `measure-the-ceiling-first` · `a-broad-population-carries-no-prior` · `attribution-must-survive-a-shuffle` · `score-against-truth` · `zero-target-guards-are-one-sided` · `hard-labels-miss-soft-change` · `never-pool-the-strata` · `a-threshold-on-a-fitted-residue` · `excluding-a-population-hides-it` · `name-the-observable-per-site` · `starved-is-not-depleted` · `the-substrate-knob-fought-back` · `key-on-a-realised-quantity` · `price-the-halves-separately` · `panel-before-src` · `admitting-an-object-costs` · `substitution-understates-a-source` · `a-symptom-is-not-a-second-defect` · `a-locked-object-is-not-a-control` · `draining-breaks-the-oracle` · `an-equal-length-panel-defeats-the-lift` · `a-length-gap-bypasses-calibration` · `weight-it-like-the-consumer` · `a-support-ceiling-is-the-clamp` · `log-variance-is-not-linear` · `re-record-the-baseline` · `a-truth-table-of-aggregates` · `a-single-level-panel-cannot-see-a-constant` · `score-the-consumers-own-count` · `the-floor-must-reproduce-the-selection`

- **C. Pools, selections and divisors** — `a-cancellation-is-conditional-on-its-support` · `a-better-estimator-inside-a-weak-consumer-moves-nothing` · `a-pooled-rate-cannot-see-a-short-object-factor` · `two-estimators-of-one-rate-weight-the-field-differently` · `state-the-population-rule-do-not-inherit-it-from-a-table` · `two-divisors-opposite-sign` · `frame-free-is-not-assumption-free` · `a-purity-filter-is-a-length-filter` · `pure-and-length-censored` · `divide-by-a-probability` · `opposite-tilts-must-not-pool` · `a-mean-of-ratios-inherits-the-partition` · `a-trap-names-the-defect-not-the-repair` · `a-stale-gate-accuses-the-newest-change` · `an-upper-bound-is-not-an-estimate` · `a-gate-on-the-helper-is-not-a-gate-on-the-caller` · `fractional-mass-is-the-problem` · `conservation-misses-mis-attribution` · `a-guard-outlives-its-divisor` · `a-fold-grows-a-heuristic` · `a-ratio-needs-a-population-that-can-supply-its-numerator`

- **D. Estimation and solver design** — `purity-is-a-property-of-the-annotation` · `pair-count-weighting-lets-one-seed-decide` · `we-keep-re-deriving-message-passing` · `one-hop-lifted-out-is-still-the-relay` · `a-variance-cannot-fix-a-bias` · `two-gaussians-one-latent` · `variance-fitted-on-the-belief` · `a-message-from-the-destinations-belief` · `a-total-density-ratio` · `substitute-the-definitions-first` · `the-pin-had-a-fixed-point` · `no-belief-not-no-numbers` · `a-licence-with-no-floor` · `a-multiplication-gated-by-a-trace` · `all-small-singly-large-jointly` · `recompute-from-the-oracle` · `a-cancelling-defect-pair` · `zero-the-precision-with-the-value` · `no-prior-means-haldane` · `prefer-shares-to-differences` · `an-all-zero-factor-is-inert` · `density-below-one-fragment-length` · `identical-paralogs-are-bimodal` · `a-mean-hits-the-mass-weighted-centre-by-luck` · `a-clamp-at-the-closed-end-escapes-the-window` · `the-deconvolution-is-as-good-as-the-density-it-is-handed` · `deriving-one-coordinate-propagates-its-error` · `interpolate-on-the-axis-where-the-lattice-is-uniform` · `read-the-whole-failure-list` · `a-priors-curvature-is-not-the-datas-information` · `a-refutability-test-needs-the-refuting-channel-in-the-fixture` · `a-strength-is-a-nat-a-prior-weight-is-a-count` · `a-four-decimal-print-is-not-a-zero` · `a-constant-in-exact-arithmetic-is-not-constant-in-float64` · `a-toy-and-a-panel-can-disagree-in-rank` · `a-rescale-that-reads-the-source-belief-is-unbounded` · `a-face-total-is-not-a-total-without-its-flux` · `an-imputation-must-cost-something-every-hop` · `a-floored-knob-is-not-the-bandwidth` · `a-mode-count-is-not-a-well-posed-quantity` · `measure-a-default-flip-before-you-write-it`

- **E. Structure, indexes and plumbing** — `one-reference-hides-refid-bugs` · `annotated-is-not-genomic` · `an-sj-is-not-a-gap` · `deposit-at-the-sj` · `splicing-makes-the-graph-cyclic` · `nrna-does-not-mean-synthetic` · `credit-exactly-one-sj` · `strand-completes-the-sj-key` · `a-hash-that-misses-its-artifact` · `integer-channels-reproduce` · `worktrees-run-the-wrong-code` · `checkout-deletes-uncommitted-work` · `two-masks-one-name` · `two-docstrings-one-quantity` · `a-transcript-predicate-must-not-silently-drop-a-molecule` · `an-object-class-does-not-see-a-terminus`

- **F. Domain facts that read like defects** — `specificity-and-sense-are-complements` · `strand-measures-the-tilt` · `a-linear-likelihood-emits-a-sign` · `a-pooled-conversion-applied-per-component` · `capture-inverts-the-counted-side` · `equal-lengths-carry-no-composition` · `capture-is-1000x-on-exons` · `capture-selects-for-length` · `on-target-by-start-is-geometry` · `eff-lengths-do-not-cancel-at-an-end` · `configured-lengths-are-not-realised` · `mature-rna-never-crosses-a-boundary` · `a-boundary-with-rna-is-not-an-sj`

- **G. Process** — `no-magic-numbers` · `one-thing-varied` · `a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism` · `converge-and-delete` · `the-source-does-not-cite-docs` · `real-data-is-a-test-input` · `running-an-arm-is-a-fresh-process` · `shard-an-arm-sweep-by-condition` · `no-enumeration-without-a-census`

---

## A. Validation and gates

**self-checking-validator. A validator that calls the builder's own helper validates nothing, and is worse
than no check because it reads as one.** Only re-deriving the result by a different algorithm caught a
deleted coordinate swap that a 1,289-test suite waved through. Emit all classes always (a validator
comparing only non-empty classes never sees a spurious flag), and prove every validator fires by
corrupting its input.

**perturb-every-gate. Writing the falsification test first is half the discipline; the other half is
breaking the fixed code and watching each gate fire.** Perturbation has found holes in already-green gates
repeatedly, once 7 of 9, including a gate that did not fire on its own named perturbation because a
redundant backstop was silently doing the guard's job.

**a-field-driven-gate-is-atomic. A gate that enumerates the specification's fields fuses spec,
implementation and schema into one commit, so plan the change that way or the tree is red in between.**
Adding one field to the specification turned 10 tests red across two files, and the digest bump invalidated
every cached scan; budget the rebuild as part of the commit. The cheap probe: add the field, run the suite,
read the failure count, revert.

**waive-with-a-measurement. An assertion the shipped code violates is waived with its measurement, never
widened.** Measure the rate before touching the predicate: either the predicate mis-stated the coordinate's
domain (fix it, from the shipped builder) or the code is wrong (waive it, with the number and a written
reason per waiver). Five such gates became six measurements on their first run. Derive any tolerance from
the coordinate's own resolution, one grid spacing, never tune it.

**a-docstring-that-misdescribes-the-graph. A claim about the import graph inside a docstring that nothing
gates rots exactly like a stale doc citation, one layer down.** Fourteen module docstrings named a sibling
with no import boundary and six were genuinely stale. The graph is in the AST, so the prose can be checked:
`scripts/design/module_census.py` does, and its count is a worklist rather than a verdict, because a
data-flow claim is one only a human can judge.

**a-flat-pile-is-not-a-knot. Before merging modules, ask whether the problem is entanglement or missing
order; merging a flat pile makes bigger files with the same problem.** The calibration package read as
unmaintainable yet had zero import cycles and most modules had exactly one importer. Naming the layers
already present in the boundaries cost zero behaviour change and immediately found two upward imports.

**the-rename-that-corrupted-a-diagram. A mechanical rewrite over prose hits tokens that are not prose, and
the ones it hits are the ones you did not enumerate.** Write the forbidding gate first, then treat its
residue as the list of collisions worth understanding rather than exempting. A citation rename rewrote 20
slot ids inside ASCII chain diagrams because the same token was a label in one namespace and a slot in
another. *Sibling:* `TRAPS: two-masks-one-name`.

**a-constant-parked-a-value-off-a-knife-edge. Removing a conjured constant can look like a regression when
the real defect is a discontinuity the constant happened to hide.** Deleting a shrinkage prior let one
row's overdispersion fall to its honest zero, and a downstream path with a threshold rather than a
response jumped. Bisect the parameter: if an infinitesimal move produces a finite jump, the defect is the
edge and no comparison across it is attributable. The tell is a change appearing in arms that inject the
value you altered and so cannot have been affected by it.

**a-gate-that-restates-the-implementation. A gate whose assertion is the code's own expression cannot fail
while the code is self-consistent, so it certifies the defect and fires on the repair.** An information
gate asserted the implementation's arithmetic while the docstring and the derivation said something else,
and the estimator's standard error was √2 too wide for as long as the gate stood. Assert the property (a
closed form, or a Monte-Carlo null the statistic must match), never the arithmetic the function performs,
and gate the wiring between helper and caller, where three more mutations passed unnoticed.

**a-gate-that-already-passed. A gate that already passes is not a falsification.** A simulator gate passed
with the defect present because the per-fragment conditional was correct and only the marginal was being
discarded. Before trusting a gate set, check which of them would have failed before the fix; if none, the
gate has not been written yet.

**right-conditional-wrong-marginal. When every per-observation conditional is correct, no conditional
comparison can detect a wrong marginal; only a check on the marginal fails.** The general form of
`TRAPS: a-gate-that-already-passed`, and it says where to look.

**byte-identity-gate. A bit-identity gate has lied in both directions: an arm with zero rows scored
"identical" because the loop ran over the new arm's rows, and a stale stored baseline made unmodified HEAD
read as broken.** Require equal row-key sets, and re-record the baseline from the current tree in the same
session; if HEAD-vs-baseline is not 100 %, the baseline is what is broken.

**the-deliverable-is-not-reproducible-by-default. The shipped EM seed defaults to none with sampled hard
assignment, so an end-to-end A/B on the default config measures its effect plus a sampling draw.** Two
identical runs differed by up to 43 fragments on the gate toy. Pin the seed (or use fractional assignment)
on every measurement arm and print a reseeded noise floor beside the effect; whether the default should
change is an owner call.

**a-clip-hides-a-scale-error. A clip hides errors on both sides of it.** A `min()` clip hid an exact
factor of 2 for months because the fixtures cancelled it, and a two-endpoint clip with reversed or equal
endpoints turns the estimator into a constant, which at a zero control is the right answer. Repair
fixtures rather than relaxing assertions, treat an exact algebraic 2 as a bug, and never quote a zero
control as evidence an estimator discriminates: require it to move at a contaminated condition.

**an-ablation-that-never-ran. An ablation that never ran reads as "no effect", and an inert prototype arm
reads as a refutation of the idea.** A patch site must assert the name it patches still exists (patch via
`sys.modules`; a rename makes a monkeypatch on an imported name a no-op); every arm counts its own firings
and raises on zero; and a counter is not enough, because an arm that fired can still be a hybrid or remove
more than it names. Diff against a base run under the same config, log what each arm applies rather than
what its flags say, and attribute source against source with a per-slot identity gate.

**a-green-suite-hid-five-dead-instruments. An instrument is alive only under the configuration it is run
in, and a green suite says nothing about it, because the suite's own tests may install a policy the
instruments never see.** Instruments have died on a `src/` deletion and on a default flip with nothing
deleted, reading keys the shipped default never writes. Run the instruments after a `src/` deletion, a
rename, or a shipped-default flip; re-run every byte-identity gate after a flip; and re-derive a rebuilt
quantity by the producer's own expression, never an algebraic equivalent, because one ULP is enough.

**compatibility-is-geometry-not-composition. Any RNA-vs-gDNA question resolved by set membership will be
answered by geometry, because gDNA is compatible with everything unspliced; resolve it by likelihood,
where the populations actually differ.** A warm-start gate keyed on RNA-unambiguous support revived more
nascent entities the more gDNA the library carried, the inverse of the intent. Nascent entities get no
prior mass and earn their place by likelihood (`DESIGN.md` §0b).

**a-zero-count-is-a-measurement. A zero count is a measurement of a density, not an absence of data;
keying precision on the count makes the strongest statement in the library the quietest.** Zero fragments
over 50 Mb and zero over 200 bp are the same number to `p = n/(n·Var+1)` and opposite statements about the
world. The repair needs no constant: a Poisson rate with `a` events over exposure `E` under the Jeffreys
prior has posterior `Gamma(a+½, E)`, proper, with finite precision at `a = 0`.

**a-ratio-cannot-carry-zero. A multiplicative transport cannot carry zero: "there is none here" is
unrepresentable by a density ratio, however precisely the source knows it.** A repair to the rule above
measured nothing, because every guard on the transport path required a positive density and a
zero-density source could not participate at all. Before fixing a source that will not speak, check
whether the channel can carry the value it wants to send.

**the-divergence-was-a-barrier. Before crediting a divergence's removal, ask what the divergence was
suppressing; an infinity in a damping term is a structural gate wearing a variance's clothes.** Making a
transfer variance finite turned every zero-mass slot from a chain-cutting barrier into an unscaled conduit
and one stratum got 20–34 % worse: the infinity had been the only thing pricing the premise that a gDNA
level does not change across a boundary, which capture falsifies by orders of magnitude.

**deadband-from-the-wrong-sample. A noise deadband whose cushion is supplied by an unrelated sample size
fails exactly where that sample gets big, silently, into the honesty columns.** The strand gate's cushion
shrank as gDNA grew, phantom information then scaled with depth, and the damage was not accuracy but the
solvable fraction collapsing to zero, inflating the column used to pick the worst condition. Propagate a
variance instead of gating on it. *Sibling:* `TRAPS: honesty-metrics-reward-ignorance`.

**honesty-metrics-reward-ignorance. Every honesty metric improves as the solver stops knowing anything, so
none is readable without a fixed-denominator accuracy number beside it.** A destruction control made
accuracy 792 % worse while three of four headline honesty columns improved on every condition. An A/B that
moves the solvable fraction has changed its own denominator: quote the mass-weighted error over all
objects and the raw Σ|err| alongside, which cannot be gamed by knowing less.

**predicate-contradicts-its-docstring. A predicate can contradict its own docstring for months if something
else masks it, and the wrong version can be load-bearing; when docstring and code disagree, find out which
one the panel relies on before "fixing" either.** A mis-scoped lock declared thousands of empty exons and
introns composition-certain, yet scoping it as documented was panel-negative on every re-pricing. The
universal half: a test that hands the predicate in as an argument cannot see the caller compute it wrong,
so gate it through the production path.

**a-test-that-redefines. A test that re-derives a definition cannot detect drift in it.** A gate meant to
keep two instruments' shared class definition aligned recomputed it inline, so changing one instrument
fired nothing. A shared definition needs one home, production code, with every consumer importing it.
*Sibling:* `TRAPS: self-checking-validator`.

**a-gates-power-is-its-invariant-set. A comparison gate's power is the size of its invariant set, so a
change that shrinks that set must not ride along with the change it is meant to gate.** Bundling a
thousand-site rename into a schema re-scan would have left almost no shared banks for the identity gate to
compare, still green and proving nothing. Before bundling, count how many objects the gate will still
compare afterwards; "almost none" means the gate has been disarmed, not satisfied.

**a-gate-that-reconstructs. A gate that reconstructs a value is vacuous wherever the code decides not to
use it; gate a symmetry the code cannot fake instead.** Eleven green gates passed a deliberately mis-paired
run because on that fixture the reconstructed ratio never entered the observable. When a gate needs a
value to flow through a conditional path, it is a gate on that conditional; prefer an invariance, such as
mirror-image passes on a palindromic chain.

**off-grid-message-mode. A message mode outside its grid's domain is not a weak claim but a pin at the
boundary: the penalty has no interior minimum, so precision buys a corner, and an out-of-range mode is the
most confident statement a channel can make.** One such unit error caused 74 % of the zero-gDNA control's
error. Two corollaries: every channel a solver delivers must appear in its debug capture, or its absence
reads as innocence; and build a message in the destination's own coordinates, because a raw count tilt
inverts sign on an antisense protocol.

**a-comment-quoted-as-a-finding. A word borrowed from a code comment is not a measurement; one propagated
into the docs invented a regression that never happened, and a build plan was correctly derived from it.**
Before repeating a source comment as a finding, check what the comment was written to justify, and never
let a term of art cross from a comment into a design doc without re-deriving it from the code.

**the-intermediate-is-not-the-deliverable. An intermediate metric is not the deliverable until somebody
measures the coupling.** A −37 % pass-0 win was −4 % on the shipped solve, and the same change's regression
grew through the refit, because a pass-0 biased in the prior's own direction trains the prior to repeat
the bias. When an intermediate stage feeds a fitted stage, report both columns on every arm and never rank
on the upstream one alone.

**could-the-arm-have-fired. "The arm changed nothing" is not a control until you check the arm could have
changed something: count the opportunities the change had to fire and print that count beside the
result.** A celebrated free falsification tested nothing, because most hops were skipped by a
divide-by-zero guard and the live ones were hops the change did not touch. A fixture is an arm too: a
one-sj fixture left a deleted deposit rule fully green because no block was ever bounded twice. Ask of
every fixture whether it could have failed.

**prove-the-substrate. Prove the substrate before proving the code: when a simulated axis is the axis you
are judging, gate the simulator on it.** For two milestones the panel's post-capture fragment-length
distribution was byte-identical to its pre-capture one, and everything measured against it inherited
that. The tell was free: diff the two capture arms' truth files.

**can-the-benchmark-resolve-it. Before running a benchmark, prove it can resolve the axis you are
changing.** A suite judged a partition change for months while its fine region set was row-for-row
identical to its merged set. `scripts/design/suite_resolves.py` is the lesson made executable: every
requirement scored against its degenerate value, no tuned thresholds.

**toys-rank-hotspots-backwards. A toy ranks performance hotspots backwards; profile on cached real data.**
Between a toy and the human index the prior's EM went from 28 % of runtime to under 1 %, and a whole
analysis was spent on the toy's top hotspot.

**an-identity-with-a-qualifier. An identity that holds "over X" is a measurement waiting to be made; price
the complement, or the qualifier is a hole nobody has sized.** The conserved-mass identity's accurate
qualifier "over unspliced fragments" hid a complement of a quarter of the RNA library depositing on no
conserved bank at all. Grep your invariants for "over", "for", "assuming" and "where": each names a
population with a size nobody has measured.

## B. Measurement and inference

**measure-the-ceiling-first. Measure the ceiling before building the correction: hand the consumer the
exact answer for one channel and see what perfecting it is worth.** One channel was ranked first for two
sessions; its ceiling was worth ~1 % while an unranked channel was worth 21 %. An A/B says whether a change
helped; a ceiling says whether the work is worth starting. *Instrument:*
`scripts/design/em_fl_ceiling.py`. *Sibling:* `TRAPS: a-symptom-is-not-a-second-defect`.

**a-broad-population-carries-no-prior. A population prior only transfers information if the population is
tight; before fitting one, measure the pooled distribution's width against the per-object uncertainty,
and if it is wider there is no prior to be had.** The gDNA landscape works because gDNA is near-uniform;
the same machinery fitted to RNA, whose per-object density spans ~3 decades, measured nothing.

**attribution-must-survive-a-shuffle. An oracle arm that a scrambled oracle also wins has not priced the
truth; re-run every truth-fed ceiling with the truth shuffled before believing it.** A fitted oracle RNA
density won, but the same truth values shuffled against their opportunities beat the true shape at high
gDNA: the gain was "any shape leaning toward the vertex", not the truth. One permutation is the whole cost.

**score-against-truth. Score against truth, not against the previous run.** The simulator writes
per-fragment ground truth into the oracle BAM's read names; scoring against it turned one "the deliverable
improved" into "the deliverable got 24 % worse".

**zero-target-guards-are-one-sided. A zero-target guard is one-sided: in a library with no gDNA, any change
that lowers the estimated gDNA fraction scores better.** Score on the contaminated conditions and use
zero-target rows as false-positive checks only. A reported −13 % win was slightly worse over the full
battery.

**hard-labels-miss-soft-change. A byte-identical hard-label result is no evidence, not no change, and a net
error is only a lower bound on the real one, arbitrarily weak.** `Σ|err| / |net|` reached 274× on one arm.
Report `Σ|err|` beside every net, and the directional split beside both.

**never-pool-the-strata. Report per stratum with the denominator named: a pooled mean buries unsolvable
objects, hides a sign flip between strata, and lets one huge object set the number.** The deferred
unstranded × capture-ON stratum carries most of the transcript error, so a pooled total mostly reports the
one stratum nobody is optimising. Within a stratum, a score scoped to a mechanism's own targets is a
diagnostic and only all-live-slots is a verdict. *Sibling:* `TRAPS: excluding-a-population-hides-it`.

**a-threshold-on-a-fitted-residue. A binary cut on a fitted parameter's residue is not a population test,
and a better threshold is not the fix; propagate the fitted parameter's own variance instead.** A tiny
`tau` cut promoted objects whose own statement was over a thousand nats wide into the scored population,
hiding a million-fragment error; any floor is a tuned constant because the residue is continuous.

**excluding-a-population-hides-it. Every population you exclude from the denominator needs its own gate,
written at the same time as the exclusion.** A condition reporting a fine error on a 40 %-solvable set had
a million fragments of error in the excluded class. The reason to exclude a class is that its answer
should be uninformative, so measure how far from uninformative it actually is.

**name-the-observable-per-site. For each place a change was made, name the observable that would move and
confirm at least one gate reads it; a count of gates says nothing when they all read one downstream
publication.** Ten sound tests all read the combine's publication, and deleting the conjunct from one of
the sites left the entire calibration suite green.

**starved-is-not-depleted. "Starved" and "depleted" are different diagnoses: starved (the experiment is too
small) is fixed by more data, depleted (the biology puts nothing there) never is.** Decide which by asking
whether the count grows with the lever you have; reading depleted as starved cost half a session chasing
depth.

**the-substrate-knob-fought-back. A substrate knob can be adversarial to the very population under study;
read the sampler's weight for that population's own objects before trusting the panel.** Toy capture
probes tiled in transcript space put a split probe over every internal sj, and the split-probe gDNA penalty
then suppressed exactly the boundary-spanning fragments being measured.

**key-on-a-realised-quantity. Key an operating point on a realised quantity visible in the truth, never on
a knob you set, and print the realised value so a row that cannot reach the target says so.** An RNA
"level" set as a multiple of the off-target gDNA density gave a true `f_g` of 0.01 off capture and 0.31
under it at the same knob setting, because capture concentrates gDNA onto the exon.

**price-the-halves-separately. A diagnosed defect and a valuable one are not the same defect; price the
halves separately.** Two length models priced together looked like one finding; split, the fully
diagnosed half was the near-worthless one, 14× smaller than the other.

**panel-before-src. A toy ceiling is not a panel ceiling: run the panel arm before writing a mechanism
into `src/`, and make every claim about a mechanism name its substrate.** Recorded five times: toy-positive
changes that were panel-negative. When the panel disagrees, ask whether the defect is half of a cancelling
pair before concluding. *Sibling:* `TRAPS: a-cancelling-defect-pair`.

**admitting-an-object-costs. Admitting an object to the scored population is a cost, and a mechanism can do
it silently; report the solvable fraction beside every accuracy number.** A prototype raised solvability
five points while the edge-axis error grew sevenfold, not because answers got worse but because wrong
ones started being counted. *Sibling:* `TRAPS: excluding-a-population-hides-it` (inverted).

**substitution-understates-a-source. A ceiling by substitution is honest for a sink but understates a
message source, whose value is what it carries; for a source, pin the truth and re-solve.** Substituting
two boundaries removed 9 % of a gene's error while the object they feed accounted for 83 %. An instrument
offering substitution must say which of the two it is doing.

**a-symptom-is-not-a-second-defect. A downstream number that is wrong because its input is wrong is a
symptom, not a second defect; substitute only the upstream input and re-run the shipped downstream
function.** The effective-length shrinkage read as an independent bug; fed truth composition it returned
exactly 1, so "repairing" it would have installed a second error to cancel the first. A ceiling prices only
what is built after its injection point, and every arm patched `assemble_priors`, which runs after the
shrinkage, so it cancelled out of every ceiling ever run. Write down where each ceiling injects.

**a-locked-object-is-not-a-control. A structurally locked object keeps its pinned init and is never solved,
so it cannot be wrong and its correctness measures nothing.** The control that works is the same class
split on the variable, boundaries with vs without sj flux, and it reversed the conclusion the "healthy
twin" reading had supported.

**draining-breaks-the-oracle. A per-fragment-independent partition stops being one the moment a downstream
step conditions on the whole tally; the identity that makes a truth source valid also fixes where in the
pipeline it can be taken.** The second pass scores against the payload's own densities, so three origin
partitions drained separately are not the whole drained; measure the undrained stage.

**an-equal-length-panel-defeats-the-lift. On an equal-length panel the drained-oracle lift's span-tie
ambiguity is the common case by construction, so run both sides undrained and price the caveat instead of
asserting it small.** About 4 % of held fragments were ambiguous, with spliced records replayed into the
gDNA partition, a refusal the oracle is right to make. Equal lengths are not the panel's mistake:
`TRAPS: a-length-gap-bypasses-calibration` is this rule's mirror.

**a-length-gap-bypasses-calibration. A large gDNA-vs-RNA fragment-length gap lets the EM assign fragments
on length alone, so the tool answers correctly with calibration broken; equal lengths force calibration
to be exercised.** The tell is free: score calibration against oracle calibration beside the end-to-end
number, and a healthy end-to-end number over an unhealthy composition claim is the signature. The general
form: for each axis a substrate varies, ask which stage that axis lets a later stage answer without.

**weight-it-like-the-consumer. Weight the average the way the consumer weights it; a bp-weighted and a
fragment-weighted mean answer different questions.** An "11 % over-call" was bp-weighted; the estimator is
fragment-weighted, most of the mass sat where the effect is inert, and end to end it barely moved.

**a-support-ceiling-is-the-clamp. A support ceiling that matches the clamp is not a match; it is the
clamp.** A distribution's support "agreeing" with `max_frag_length` was recorded as a fix; the narrower
estimate had been correct.

**log-variance-is-not-linear. `Var(log f)` is not `Var(f)`; convert with the delta method
(`Var(f) ≈ f²·Var(log f)`), because at small `f` the two differ by orders of magnitude.** Every
overconfidence figure computed without the conversion is void, not merely imprecise.

**re-record-the-baseline. A delta is only attributable if its baseline came from the same tree in the same
session; re-record the before-picture, never quote a stored one.**

**a-truth-table-of-aggregates. A truth table's label column may hold nested aggregates, and summing it
double-counts, with the control arm unaffected, the most persuasive shape a wrong result can take.**
Enumerate the membership sets explicitly, raise on an unrecognised label, and put a premise gate on any
truth source. A parser bucketing `{mrna, nrna, gdna, rna, all}` as "gdna, else rna" summed mRNA twice
plus the whole library once while the gDNA arm read exactly 1.

**a-single-level-panel-cannot-see-a-constant. A panel that holds the quantity of interest fixed cannot
distinguish a good estimator from a constant that happens to equal it; list the axes the panel varies and
check the deliverable is one of them.** A channel won 7 of 8 conditions on an all-g50 panel, then reported
over half gDNA on the zero-gDNA control: the win was the panel agreeing with a near-constant ½. Three gDNA
levels is the minimum for this reason.

**score-the-consumers-own-count. Truth is not a yardstick until you name the consumer: for each number the
consumer reads, write down the set of fragments it counts and score against that set.** The first-base
count had impeccable provenance and was still not the EM's target, whose soft count includes straddling
candidates it drops; scored against it, a 72 % "open residual" spawned a research programme, and against
the EM's own count it vanished. Before measuring, write down which line of the consumer reads the number.

**the-floor-must-reproduce-the-selection. Whatever selection the scorer applies to the data, the noise
floor must apply to its draws; a floor drawn unconditionally under a conditioning scorer reads the
conditioning as a rule error.** Scored destinations were zero-truncated Poisson; an unconditional floor
printed a 10–14 % "rule error" on hops that were exact, and redrawing conditional on `n >= 1` closed it. A
floor that cannot be falsified by turning the conditioning off is not a floor.

## C. Pools, selections and divisors

**a-cancellation-is-conditional-on-its-support. A reciprocal-opportunity deposit cancels its opportunity
only where that opportunity is non-zero.** Where `A(w) = 0` the fragment deposits nothing, so
`E[Σ 1/A] = ρ · P(A > 0)`, a functional of exactly the length distribution the channel claims independence
from; on the shipped REGION bank the omission was an 11.6× under-read at the median exon. Write the
expectation with its support factor and gate it with fragments placed outside the support; `EQUATIONS.md` §2.

**a-better-estimator-inside-a-weak-consumer-moves-nothing. A severalfold repair to an estimator can be
worth ~1 % end to end, and the consumer's strength is why.** Swapping the pooled gDNA background for a
pmf-free pair improved the estimator up to 4.3× under capture while the deliverable barely moved, because
a background that enters as a weak one-sided floor cannot transmit a sharper rate. Price the consumer's
sensitivity before repairing its input, and spend the effort on whatever reads a level sharply.

**a-pooled-rate-cannot-see-a-short-object-factor. A pooled rate is a ratio of sums, so it is dominated by
the largest objects, and a per-object factor that only bites at small ones is invisible in it at any
size.** A factor reaching −99 % on short objects did not move the background pools at all, because their
exposure is dominated by megabase regions; only the exon row moved, exactly as predicted. Before building
a panel to exert a per-object mechanism, compute the pooled weight of the objects it acts on.

**two-estimators-of-one-rate-weight-the-field-differently. Two unbiased estimators of one slot's rate
agree only under local uniformity; under a probe-shaped field they are different weighted averages of it,
both right and not interchangeable.** The START and CONTAINED banks agree off capture and disagree by up
to 30× on intron regions under it, so swapping one for the other changes the quantity rather than removing
a bias. START and END share one opportunity, so their agreement is field-free; build that comparison.

**state-the-population-rule-do-not-inherit-it-from-a-table. A population defined by whichever rows a file
happens to contain is not defined, and two implementations will reach it by different doors.** Wall
distances inherited "real transcripts only" from a table that carried no synthetic transcript; an
enumeration through the index, which does return the synthetic span, disagreed on 57 distances. If you
cannot point at the line that excludes a population, you are relying on a file's contents. Write the filter.

**two-divisors-opposite-sign. Two divisors built from one pmf can still disagree if they respond to it
with opposite sign.** `E_J = E[w]−1` rises with the mean fragment length while `E_r = e−E[w]+1` falls, so a
length-model error appears as a junction-vs-exon frame gap while the two remain exactly consistent.
Differentiate both with respect to the shared model before looking for a bug in either, and read the
simulator's own sampling code, not its docstring, before pricing a selection effect: a "finite-transcript
placement" factor added here was sound arithmetic about a generative model the simulator does not use.

**frame-free-is-not-assumption-free. A term that is frame-free is not therefore assumption-free; look at
the nuisance you profiled out, not only at the divisor that cancelled.** Putting certified RNA on ψ's RNA
arm as `E[S] = c·(1−f_g)·M` is correct about the frame and still not enough: the same `c` holds the
splice-visibility `q`, and with `q` free the profile likelihood in `f_g` is flat. The tell was general:
its benefit tracked the answer rather than the evidence, and a channel whose sign follows the truth is a
prior. Ask of any one-sided floor whether the term you drop is bounded or dominant.

**a-purity-filter-is-a-length-filter. A purity filter on a length pool is a length filter.** Barring
fragments whose length was partly inferred rather than sequenced selects exactly the ones whose mates sit
far apart: the pool read −9.6 % mean and −22.5 % sd against truth where keeping them was near-exact.
Before excluding a population from a pool, ask what the exclusion criterion correlates with; if it
correlates with the axis being measured, purity and accuracy point in opposite directions. A pool is
keyed on determinacy, not provenance.

**pure-and-length-censored. A pool can be composition-pure and length-censored at the same time, and the
second is invisible to every purity argument.** "gDNA contained in an intergenic or intronic region" is
gDNA by construction and, under hybrid capture, ~15 % short, because a long fragment beside a probe reaches
the exon boundary and stops being contained. The pools it spills into resemble the RNA pool, and it was
filed for two milestones as "not gDNA". Ask what a pool's selection rule correlates with, not only what
the selected fragments are.

**divide-by-a-probability. A pool divided by its opportunity must be divided by a probability, not by a
count.** `count(w)/A(w)` recovers the distribution lengths were drawn from; every consumer needs the one
the library realises, so the divisor is `A(w)/T(w)`. The two forms differ in shape, and the ratio form is
where an abundance weighting cancels: swept to pathological regimes it is never worse than not correcting,
and the `A`-only form is.

**opposite-tilts-must-not-pool. Pools with opposite tilts must not be pooled raw.** A contained pool's
opportunity falls with length while a crossing pool's rises; summing the histograms and applying one
divisor read a gDNA mean of 146 where the contained pool alone said 88. Summing the counts and the
matching per-pool opportunities is a different operation and is correct: it is the opportunity-weighted
average of the per-pool estimates, and under Poisson counts those weights are exactly inverse-variance.

**a-mean-of-ratios-inherits-the-partition. An estimator that averages per-object ratios is a function of
where the annotation drew its boundaries.** A region boundary appears wherever a signature changes, so
how many objects a path contains is an artefact; `Σmass / Σopportunity` is invariant to it because
splitting an object splits both sums together. Subdividing one intron from 1 region to 10 at constant mass
and opportunity tripled a shadow span's share of its locus's prior while the pooled form did not move. An
estimator must be invariant to a re-partition that leaves the data unchanged, and that is testable.

**a-trap-names-the-defect-not-the-repair. Citing a trap does not license a mechanism; the trap names the
defect, not which quantity to repair.** The rule that a zero count is a measurement was cited to license
a Jeffreys posterior mean as a density location, `(mass + ½)/opportunity`; the shipped repair had put the
half on the precision, and the location form sits in the graveyard at +7,269 % on the zero-gDNA control.
Before adopting a mechanism because a trap seems to motivate it, grep the graveyard for the mechanism.

**a-stale-gate-accuses-the-newest-change. A gate whose premise expired does not go quiet; it fails, and it
blames whatever is in flight.** `rescan_panels.py` gated an irreversible rebuild on byte-identity because
"these are integer tallies"; once six banks became float64 the gate was unsatisfiable, and it then failed a
schema change on banks that change had not touched. What resolved it was a control on the instrument: scan
the same BAM twice. When a gate fails on something you did not touch, first ask whether it can still pass.

**an-upper-bound-is-not-an-estimate. A bound that is shared is not evidence about who shares it.** A
transcript's density is bounded by the minimum along its path, which is useless wherever that minimum is
attained at an object the transcript does not own: three quarters of silent transcripts share an object
with an expressed one, inherit its bound, and are asserted into existence, so false-positive mass doubled
while gene-level error halved. The tell was a zero-weight set byte-identical across every arm: when an
estimator's support does not move as you vary it, you are varying the wrong thing.

**a-gate-on-the-helper-is-not-a-gate-on-the-caller. Testing the function that computes a quantity does not
test the call that asks for it.** A weight builder shipped with 29 green gates and two holes: the gates
called the helper directly, so replacing the caller's weights with ones fired nothing, and two arrays were
compared on a single-exon transcript, where they are equal by construction. Gate through the public entry
point, on a fixture from a population where the two branches give different answers.

**fractional-mass-is-the-problem. Fractional mass is the partitioning problem.** A fragment spanning four
regions writes six fractional numbers whose values depend on region sizes, purely because a mass is
conserved; the same fragment's three crossing counts depend on nothing. Multimapper and ambiguous-path
assignment must stay integral, or the non-integer observable returns and the count stops being a count.

**conservation-misses-mis-attribution. Mass conservation does not catch mis-attribution, and an entire
alternative deposit rule conserves.** Injecting the rejected `1/K` rule in place of the shipped
`slice_len / L` left every conservation gate green; what separated them was a re-derivation on a different
axis (each fragment base attributed to the boundaries bounding its own region), which `1/K` cannot express.
If a property is invariant across the rules you are choosing between, a gate on it does not gate the choice.

**a-guard-outlives-its-divisor. Delete the divisor and the guard against it goes inert while its test
keeps passing.** `_mass_where_there_is_opportunity` existed because `rho = Σm/ΣS` is a rate; re-basing the
prior onto a conserved count removed the division, and the guard's test kept passing because its fixture
asserted zero where the mass really was zero, while the guard's correct behaviour had reversed. A guard's
justification names an operation; when that operation leaves the function, re-ask whether the guard still
has a subject. Only injecting the defect tells you.

**a-fold-grows-a-heuristic. A quantity folded onto another axis to fit a consumer's interface will grow a
heuristic to repair the fold, and the heuristic will outlive every memory of why the fold was there.** A
contiguous BOUNDARY has no extent, so a bp-overlap projection folded each boundary's mass into one flank;
the patch that repaired the fold routed most of that mass into the wrong gene. Tells: a helper justified
by a consumer's limitation; a special case that undoes the general case; two callers where the fold is
right for one and wrong for the other. Give the consumer the missing axis instead. `DESIGN.md` §3.1b.

**a-ratio-needs-a-population-that-can-supply-its-numerator. Adding opportunity that structurally cannot
carry counts manufactures a deficit of exactly its own share, and it reads as a defect in the estimator.**
An A/B that pooled every reference where the shipped gate stratifies folded the panel's ERCC exons, real
`eff_gdna` with zero gDNA, into the `R exon` class, which was precisely the "2 % deficit" written up as a
second bug and retracted. A class ratio is meaningful only over objects that could supply the numerator,
and an A/B of a stratifying gate must stratify the same way.

## D. Estimation and solver design

**purity-is-a-property-of-the-annotation. No object class is pure gDNA.** "Intergenic" is whatever the
user's GTF leaves over, pervasive transcription is real, and most genes are off in any one sample but
nobody knows which, so an estimator that asserts a structural class clean is calibrated to the annotation,
not the genome; on a chromosome fed unannotated transcripts every purity-based fit moved and the away-half
moment did not. Where a contaminant can only push a statistic one way, use the half it cannot reach.

**pair-count-weighting-lets-one-seed-decide. Pooling a second moment weights each object by its pair
count, proportional to n², which is minimum-variance only when the parameter is zero.** Once it is not,
an object's information saturates with depth, so a handful of deep objects carry the estimate: on real
cfRNA one seed carried 78 % of a library's numerator. The fix is inverse-variance weighting
(`EQUATIONS.md` §6b), never trimming; and no simulated panel can expose this, because at a true zero no
seed can dominate.

**we-keep-re-deriving-message-passing. Several sessions have independently re-derived message passing
from scratch; nobody reasoned wrongly, everybody was slow.** The tell: you are reasoning about how an exon
gets its gDNA level (a pooled reference, local imputation, a throttle, a bound) and have not yet written
the words "message passing". Stop; you are inside it. gDNA is measurable only where no mature transcript
crosses, while an exon's unspliced mass is gDNA plus RNA, the unknown itself, so an exon's level must be
imputed along the chain, and that is message passing. It is built and ships as the `transfer` policy in
`calibration/messages/`; with stranded data most exons solve directly, so it matters most for unstranded
data and AMBIG slots. Find the ruling with `grep -in 'message passing' docs/DESIGN.md docs/EQUATIONS.md`.

**one-hop-lifted-out-is-still-the-relay. A hop implemented outside `calibration/messages/` is still message
passing, and belongs there.** A value computed at one object and consumed at another is a message whether
it carries a count, an sj flux or a belief, and "it is only one hop" is not a definition either. The tell
is a justification beginning "this isn't really a message, because": a private copy of the framework is
a second home for one mechanism with none of its pricing, licences or gates. Add the hop to the policy
and price it there.

**a-variance-cannot-fix-a-bias. You cannot fix a biased mode with a variance.** Under capture a counting
estimate was systematically ~2× low but precise: both flanking boundaries sat at the same enriched
boundary and agreed on the same biased-low density, so the bias was trusted. A disagreement-based
variance model structurally cannot fix a bias.

**two-gaussians-one-latent. Never hand a solver two Gaussians built from one latent.** A message on `log f`
and one on `log(1−f)` are rank-1 with correlation exactly −1, so adding their Fisher information is
exactly 2× over-confident, rising to ~7× with deep spliced content.

**variance-fitted-on-the-belief. Never fit a variance on the current, not-yet-solved belief.** Adjacent
wrong regions agree, so the variance collapses, the messages turn confident and the error propagates.
Any component trained on the solver's own output is self-confirming: refit iterations went monotonically
worse.

**a-message-from-the-destinations-belief. A message computed from the destination's own belief carries
zero information and confirms the destination.** A message may use the destination's constants and its
observations, never its beliefs, and any "fix" that divides the belief back out rebuilds the bug. The
check is whether the delivered value is independent of the destination's own state; a prediction that
does not move when the data move by four orders of magnitude is the tell. A corollary: capping `Var(f_g)`
at `f_g(1−f_g)` asserts certainty at the corner where an evidence-free init parks, so a composition
variance comes from the reference prior's own spread. The nine rules that follow are this lesson's other
costumes, kept because the source cites them.

**a-total-density-ratio. A scale factor must be built from the component the claim is about, never from a
total.** `r = ρ_tot(dst)/ρ_tot(src)` re-creates the parent bug whenever the total is dominated by a
component the message is not about: a correct gDNA density of 0.026 was delivered as 28.7.

**substitute-the-definitions-first. Before correcting an operator, substitute its own definitions and read
what it delivers.** Substituting showed `ρ_c(src)·r ≡ φ_c(src)·ρ_tot(dst)`, a pure composition imputation
with no level in it, so no corrective factor existed; two sessions had been spent hunting a better `r`.

**the-pin-had-a-fixed-point. A rescale of the form `k = 1/(φ_msg + R_own)` has a fixed point at
`R_own = ½` and drives the delivered fraction to ½ regardless of truth.** It hid from every aggregate
because per-step rescales telescope back to 1 at a gene's far end; only a per-object check away from a
pure-gDNA object sees it.

**no-belief-not-no-numbers. State the licence as "no belief may enter", not "the destination's numbers may
not enter".** Gating a rescale wherever it read the destination's own density looked clean and broke the
capture landscape, because the off-probe floor leaked into every exon.

**a-licence-with-no-floor. A licence that tests a precision with no floor is granted by it.** A fitted κ
within a millionth of ½ on a genuinely unstranded library leaves `I(f_g) ∝ (2κ−1)²` essentially zero yet
strictly positive. The repair is to propagate `Var(κ̂)`, not to add a floor.

**a-multiplication-gated-by-a-trace. A predicate gating a multiplication must be sized by how much density
stands behind it, not by whether a precision is strictly positive.** An intron's trace of phantom RNA was
the only nonzero RNA precision in a chain and alone unlocked a reframe compounding a gDNA level to 2.16×.

**all-small-singly-large-jointly. All-small-singly plus large-jointly means stop ablating consumers and go
one stage upstream.** Removing each ψ channel in turn moved the worst object by under a tenth of its
error while removing all of them moved it the whole way, because all three were built from one delivered
level.

**recompute-from-the-oracle. Recompute the quantity from the oracle before assuming the formula is
broken.** An impossible-looking enrichment ratio with capture off was correct: `ρ_tot` is a total and the
destination held RNA fragments where the source held none.

**a-cancelling-defect-pair. Fixing one of two errors that cancel is worse than fixing neither.** Correcting
one hop alone more than doubled a toy's evidence-free exon error while the rung it targeted improved.
Price such a fix in the arm that also removes the other defect; a cancelling pair is one experiment, not
two.

**zero-the-precision-with-the-value. A refused claim must lose its precision in the same statement that
zeroes its value.** A value zeroed at one line and a precision handed back at a later one is the confident
zero the licence forbids, and a downstream rescale converts it into its opposite: a delivered `n = 0` at
high precision became "all your mass is gDNA" one hop on. Every operator that grants a precision runs
before the licence's zeroing, the fixture carries that operator's input, and the mechanism is scored where
the damage lands, one hop past the confident zero. *Sibling:* `TRAPS: a-licence-with-no-floor`.

**no-prior-means-haldane. "No prior" does not exist on a grid; omitting a term lets the grid supply
Haldane** (`p(x) ∝ 1/x`, improper, an amplifier toward the vertices). Posterior median spread over grid
half-widths 4–20 is 0.045 at Jeffreys and 0.525 at Haldane.

**prefer-shares-to-differences. Sums are well conditioned; differences are not.** Subtracting across a sj
gives `Var(log ρ) = u²σ_T² + (u−1)²σ_μ²` with `u = 1/(continuing share)`; at the real median `u` is already
at the boundary of validity and at the upper quartile hopeless at any depth. Prefer shares.

**an-all-zero-factor-is-inert. An all-zero factor is uninformative, not decisive.** In a multiplicative
score, a factor that is zero for every candidate annihilates the other factors and collapses the record
to a coin toss. Skip a flat-zero factor; do not multiply by it.

**density-below-one-fragment-length. Density below one fragment length is not resolvable by any design.**
A 1 bp region has no independently measurable density and never will, though composition still does,
since it depends on what the fragments are rather than where. An object with zero opportunity for a
component must emit nothing at zero precision, not a floored division; and "no data" must be inert, never
"100 % gDNA", which was actively seeding false gDNA into neighbouring exons.

**identical-paralogs-are-bimodal. An identical-paralog split is bimodal, and depth does not fix it.** Which
branch it lands on is a draw, so do not move a seed or a depth until it lands even; that is tuning to
green. Assert the collapse structurally (`min == 0`, `max == total`), which still fires the day the tie
breaks, and score only what is identifiable. Two defects live here: the total was also over-called, and
the collapse follows a real gradient with iterations, so a tilt enters from outside the RNA likelihood.

**a-mean-hits-the-mass-weighted-centre-by-luck. A scalar that "wins" on a mass-weighted metric may only be
sitting on the mass-weighted centre.** Measured four times in one investigation, because each population
was bimodal at `f_g ≈ 0` and `f_g ≈ 1` and the metric weights by mass, which sat at one mode. The tell:
substitute the candidate by its own class mean; if that costs nothing it is a constant, and its sd
against the truth's says so directly. Report both beside any per-object claim.

**a-clamp-at-the-closed-end-escapes-the-window. A prior clamped at the closed end of its support can put
almost all its mass outside the solve window, and a sweep over the interior will never see it.** ψ's
reference was clamped at `m = 1 − 1e-8`, leaving under 1 % of its mass inside the shipped window, while the
invariance gate swept the interior and the arm operated only at the clamp. Test a prior at the value the
code will actually use, and derive its floor at the object: one pseudo-fragment gives
`m_i = E[g]_i/(E[g]_i + 1)`, which keeps most of the mass inside the window where a pooled floor keeps none.

**the-deconvolution-is-as-good-as-the-density-it-is-handed. A residual estimator inherits the error of the
density it subtracts, amplified.** `f_g = ρ_g·E_g/M` needs no RNA model, which makes it usable pre-solve,
but its accuracy is exactly `ρ_g`'s: off capture the density is right and the deconvolution scores well;
under capture the density is a third of truth and the same deconvolution is worse than the uninformative
reference. Price the density against truth before believing anything built on it, per stratum.

**deriving-one-coordinate-propagates-its-error. Closing a composition by deriving one coordinate from
another makes the derived one inherit every defect of the source, including ones the old estimator did
not have.** Defining the RNA total as `1 − f_g` removed a real defect, but `f_g` was the grid-snapped
median where the retired read-out had been exact at ½ on an evidence-free object. When you make coordinate
B a function of A, enumerate what B used to be exact about and check A is exact about it too; the repair
is to fix the source. *Sibling:* `TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`.

**interpolate-on-the-axis-where-the-lattice-is-uniform. A quantile interpolated on a non-uniform transform
of the grid is biased toward the middle, and it looks right.** ψ's ½-quantile was interpolated in
`f_g`-space, where the lattice spacing varies by four orders of magnitude, so a posterior concentrated on
one grid point came back biased toward ½; interpolating on `λ`, where the lattice is uniform, returns that
grid point exactly. Do quantile arithmetic on the axis the grid was built on, then map. The case is not
synthetic: an unsolved slot's fed-back belief is a one-hot posterior.

**read-the-whole-failure-list. Derive the failure set, never eyeball the tail.** A run reported as "21
golden failures" from the last screen of output held two real ones, one of them the most informative
defect of the session. `pytest -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c` prints one
line per file and fits on a screen at any failure count.

**a-priors-curvature-is-not-the-datas-information. A prior may not contribute to a Fisher precision.**
`region_init`'s `tau_lam` is the data's information on the composition axis; letting a prior's curvature
into it was refused on three counts: the observed fall is the Jacobian, the term acts as a boolean gate
flip rather than a contribution, and it credits data-free slots. Before adding any term to a precision,
ask whether it scales with the data; if not, it is a prior and its place is the posterior.

**a-refutability-test-needs-the-refuting-channel-in-the-fixture. A prior measured without the evidence
that could overturn it is not being tested; it is being asked to answer alone.** ψ's structural prior read
severalfold worse on a toy chain with no intergenic regions, where the background returned uninformative
and at κ = ½ the strand channel is dead too; rebuilt with intergenic flanks, which production always has,
it yields to the same evidence. Assert the overturning channel is live in the fixture before measuring.

**a-strength-is-a-nat-a-prior-weight-is-a-count. A strength in nats and a prior weight in pseudo-counts are
not the same unit, and equating them is an analogy wearing a derivation's clothes.** `m = σ(a+b)` reads as
principled and is dimensionally wrong; the derivation that stays in one currency asks what mean one
pseudo-observation would produce. Taking the strength from the lattice was worse than no prior: a prior's
strength comes from what it claims, never from what the representation can hold, and a cap is not a choice.

**a-four-decimal-print-is-not-a-zero. A diagnostic printed at four decimals read `0.0000`, and a session's
whole mechanism was built on it being zero.** The value was 4.2e-05, the predicate the diagnosis turned
on returns True on it, and the real effect was a change of magnitude needing a different repair and, as it
turned out, none. A number a diagnosis turns on is `repr`'d, never formatted, and "exactly zero" is
asserted with `== 0.0` rather than read off a table.

**a-constant-in-exact-arithmetic-is-not-constant-in-float64. "This term is constant so it cancels" is a
claim about ℝ, and the grid does not live there.** ψ's neutral location term is `log 2` exactly, but
`logaddexp` leaves that row one ULP from flat, and since the median reads the first grid point whose CDF
reaches ½, balanced AMBIG slots moved a full grid step. Where a term is documented as vanishing, return
the vanishing value, and gate the reduction on the answer (`solve(neutral) == solve(None)`, bit for bit).

**a-toy-and-a-panel-can-disagree-in-rank. A change measured on a toy can not only shrink on the panel but
invert its order.** The best arm of three on the test chromosome was the worst of three on every in-scope
stratum of the ladder. A toy's objects mostly have little own evidence, so a message is nearly free
information; the panel's mass sits on objects that can answer for themselves, where the same message
competes with a measurement. A toy result is a mechanism check, never a ranking: run the panel before any
"this beats that", and split the error by whether the destination had its own evidence.

**a-rescale-that-reads-the-source-belief-is-unbounded. A transport whose factor is derived from the
source's own claim amplifies a weak claim without limit, and the weaker the source the larger the
factor.** A mass-identity rescale `k = M_dst / Σ_c ρ_c,src·E_c,dst` is exact under its premise and never
touches the destination's belief, and still handed a large exon "all your mass is gDNA" on a zero-gDNA
library, because the source's claim is its denominator. Transport by a quantity measured at both ends,
never by one the source's belief appears inside. *Sibling:* `TRAPS: zero-the-precision-with-the-value`.

**a-face-total-is-not-a-total-without-its-flux. Comparing two objects' totals when one of them cannot
hold part of the population reports a difference that is structure, not enrichment.** Mature RNA cannot
cross an `exon|intron` boundary contiguously, so it appears there as sj flux rather than as crossings,
and an enrichment ratio built from the boundary's unspliced abundance read a 30× "depletion" with no
probes in the condition at all. A face's total is its contained/crossing abundance plus the sj flux whose
bodies lie on that side, and that flux belongs to one face; pooling it over-counts against both flanks.

**an-imputation-must-cost-something-every-hop. If every variance in a message layer shrinks with counts,
an imputation between two deeply-counted slots crosses for free and arrives at full strength beside a real
measurement.** The premise of a hop, "my neighbour's values apply here", has a variance that does not
shrink with depth and is charged on every hop, fitted by moments as `max(0, Var(log r) − mean(v_r))`.
Split a message layer's error by whether the destination had its own evidence: "the messages help" and
"the messages trample a measurement" are different findings.

**a-floored-knob-is-not-the-bandwidth. A smoothing constant that is clamped for almost all the data is not
the bandwidth, and sweeping it answers nothing.** `landscape.knn_widths` returns `max(scale · d_k,
grid_step)`, and on the ladder 99 % of kernels sit at the floor, so the smoothing in force is the grid
step. The tell: a sweep returns the same score at several adjacent settings, a disconnected knob rather
than a flat optimum. Measure the share of kernels at the floor before reading any bandwidth sweep.

**a-mode-count-is-not-a-well-posed-quantity. How many modes a fitted density has is a property of the
render resolution, not of the field, so no consumer may depend on it.** Across a 16× range of `_N_GRID`
the `AbundanceLandscape`'s mode count tracks 1/step and never converges, buys nothing in held-out
likelihood, and reproduces worse across split halves as it grows. Sweep the resolution before a consumer
reads any shape statistic off a fitted density: `rho_0` and the containment verdict survive, a ratio of
two selected modes does not. *Sibling:* `TRAPS: a-floored-knob-is-not-the-bandwidth`.

**measure-a-default-flip-before-you-write-it. A config default and a refusal are one design, and flipping
the first can invalidate the second.** `CalibrationConfig.abundance_landscape` was written opt-in and
refused to run without the wall arrays; flipping the default broke 65 callers on one cause, that unit and
toy fixtures have no wall arrays. The refusal serves the population that opts in; the default serves the
population that never thought about it. Re-read what the refusal was protecting against (here the
alternative to refusing was not to fit, so a missing input skips loudly) and keep the asymmetry gated.

## E. Structure, indexes and plumbing

**one-reference-hides-refid-bugs. Single-reference synthetic indexes hide reference-id-space mismatches.**
A resolver assigning ref-ids by first-seen order silently dropped nearly every real fragment while every
golden test passed. RNA-only spike-ins do nothing for the gDNA deposit path, which needs at least two
genomic references.

**annotated-is-not-genomic. "Has an annotation" is not "is genomic"; a classification must be an input,
not a proxy.** Choosing gDNA references as `{t.ref for t in transcripts}` filled the panel with gDNA on
the RNA-only spike-ins whose zero truth is the false-positive control. A mis-stated classification must
raise, not silently produce nothing.

**an-sj-is-not-a-gap. A splice junction cannot be detected from a gap between deposited slices.** A
contiguous spliced read whose exon body straddles an internal region bound has no gap at all. The sj's
identity is the cut-intron coordinates, which the scanner already has; pass them through.

**deposit-at-the-sj. A splice deposit belongs at the sj's coordinate, not the region's boundary.**
Invisible for annotated introns, whose ends are region bounds; for an unannotated sj the mass lands
kilobases away.

**splicing-makes-the-graph-cyclic. Alternative splicing makes a region↔sj graph cyclic, and cycles are
the common case.** The human graph has on the order of 400,000 independent cycles, one per sj, and
two-sweep forward-backward is exact only on a tree. Never break a cycle by dropping a sj boundary; that
re-isolates the exon the boundary exists for.

**nrna-does-not-mean-synthetic. `~is_synthetic & ~is_nrna` as the "real transcript" filter deleted tens
of thousands of real termini.** On a non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent",
not "manufactured span". One filter: `~is_synthetic`.

**credit-exactly-one-sj. A fragment crossing K sj must credit exactly one, the leftmost annotated.**
Crediting all K shifts the library sense fraction by tens of percent and creates between-side correlation
that reads as overdispersion on a zero-overdispersion simulator. A `1/K` split is provably biased.

**strand-completes-the-sj-key. `(src, kind, dst)` is not a total order for sj boundaries.** Two
strand-coincident sj differ only in strand, so ordering becomes input-order-dependent and the duplicate
check reads them as duplicates. GENCODE has none, so only a synthetic stress test finds it. Sort on
`(src, kind, dst, strand)`.

**a-hash-that-misses-its-artifact. A cache key must cover the artifact it caches, in all three forms.** A
hash stored beside its data verified a stale cache clean when a fix rewrote files it did not cover, so
compute hashes on demand. An exclusion from a key claims nothing the excluded file produces is in the
blob; when two of its products were, the cache served one fresh arm beside stale ones. A dataclass field
computed at construction survives `dataclasses.replace` and then describes the arrays it replaced: if a
value's inputs are in an override list, it must be a `@property`.

**integer-channels-reproduce. A count is an integer and reproduces bit-identically across worker counts;
a fraction is float64 and does not.** One numeric convention, with tests validating within a derived
tolerance bracketed from both sides. The fixed-point alternative was less accurate than the float it
defended against. The one bank whose order is observable must be sorted on its own content before
crossing the ABI.

**worktrees-run-the-wrong-code. Worktrees silently run the wrong code.** An editable install's meta-path
finder beats `PYTHONPATH`, so an A/B inside a git worktree executes the main repo's source.

**checkout-deletes-uncommitted-work. `git checkout -- <file>` does not undo a perturbation when the work
is uncommitted; it deletes the work.** A perturbation harness must restore from a copy of the working
tree. Cost one full re-implementation.

**two-masks-one-name. Two different masks shared the word `struct_lock`, and both were right.** One meant
"pinned and certain", the other "may emit composition certainty". Two correct predicates under one name
is worse than either being wrong: give each its own name and one home, and let every consumer import it.

**two-docstrings-one-quantity. The prose next to the code said "the average" and the code followed the
prose, while a sibling module's docstring had the correct formula the whole time.** Two docstrings
disagreed about one quantity for months and nobody diffed them.

**a-transcript-predicate-must-not-silently-drop-a-molecule. A fragment rejected by a transcript-level
predicate deposits nothing, and the rejected population is never random; the bias concentrates where the
predicate fires.** A chimera check that joined mates only when their transcript sets intersect dropped
0.04 % of fragments carrying 2.4 % of every boundary crossing, and hid because every dropped fragment was
a crosser. Check genomic compatibility (one reference, mates inward, plausible implied length) before
calling a fragment a chimera, and make the ledger close: `n_fragments == deposited + Σ accounted drops`.

**an-object-class-does-not-see-a-terminus. The object strata classify a boundary by the exon-ness of its
two flanks, so a TSS/TES inside another transcript's intron is indistinguishable from a splice site there,
and the two have opposite currencies.** At a terminus the RNA originates and a composition cannot cross;
at a splice site it enters by the sj and the composition does. A hop type is `object class × {sj, term,
sj+term}`, read off `RegionStatics.boundary_flags`; the classes stay right for ψ's reference, the hop
needs the second bit. *Sibling:* `TRAPS: two-masks-one-name`.

## F. Domain facts that read like defects

**specificity-and-sense-are-complements. Strand specificity is two different quantities and they are
complements.** A simulator's `strand_specificity` is direction-agnostic protocol fidelity; a fitted sense
fraction is directional, so on an R1-antisense (dUTP) protocol a fitted κ of 0.01 on a "0.99 stranded"
library is correct, and forcing 0.99 measured far worse. The matching quantity is
`StrandModel.strand_specificity`; the recovered values live in the test that measures them.

**strand-measures-the-tilt. Strand measures the tilt, not the gDNA fraction.** With RNA tilt `d = f₊ − f₋`,
`p = ½ + (κ−½)·d`: the gDNA fraction cancels identically, and the information about it is exactly zero at
κ = ½ for any count. Strand reaches gDNA only through the triangle bound `f_g ≤ 1 − |d|`.

**a-linear-likelihood-emits-a-sign. A likelihood that is asymptotically linear in its parameter has no
mode, only a direction, and on a bounded grid that direction saturates at an endpoint.** If the leading
term of the log-likelihood is linear, the object is a vote, not an estimate: summarise it as a location
plus an information, never a row handed to a normaliser. Its amplitude, participation and declared
precision do not fade together (a vanishing amplitude sat beside participation of the whole library and a
precision that was the grid's own width), so report all three. The channel that taught this is gone.

**a-pooled-conversion-applied-per-component. A ratio measured on the pooled population and applied to each
component separately is population-blind, and the blindness need not be the axis you expect.** The
prior's crossing→fragment conversion errs worst at equal lengths: the driver is placement, gDNA crossing
boundaries in long intergenic regions and RNA in short exons. Bounded at under a percentage point of
composition, at calibration's noise floor: record the bound, build nothing.

**capture-inverts-the-counted-side. "The well-counted side" is not a fixed side; capture inverts it, so a
rule that transports from the well-counted side silently reverses direction with the protocol.** Off
capture an intron holds hundreds of gDNA counts against a flanking boundary's dozens; under capture the
intron holds one and the boundary dozens.

**equal-lengths-carry-no-composition. At equal component mean lengths the length channel carries exactly
zero composition information, at any depth; the deconvolution is identified only through `μ_g − μ_r`.**
This is the gDNA ladder's forcing function, not a caveat about it: zero length information is what stops
the EM answering without exercising calibration.

**capture-is-1000x-on-exons. Hybrid capture is ~1000× on exons and only gDNA reads it cleanly.** RNA's own
expression range hides the probe pattern; gDNA's uniform baseline does not. Capture also destroys the
intron signal ~75×, so nascent-vs-gDNA is unidentifiable under capture, one more reason the nascent scope
ruling puts dedicated nascent-capture experiments out of scope.

**capture-selects-for-length. Capture selects for length, and the post-capture distributions are the
baseline.** A short fragment presents less sequence and is captured less efficiently, so the pre-capture
parameters describe a library that was never sequenced; score against post-capture truth, for lengths
exactly as for abundances. Capture narrows the gDNA↔RNA length gap whenever gDNA is the shorter component.

**on-target-by-start-is-geometry. "On-target" defined by the start's territory is geometry, not capture
efficiency.** Conditioned on capture, an intronic-start fragment was long enough to reach the probe
(weight ~`w²/2`) while an exonic start is flat in `w`, so intronic-start populations read longer under any
capture model. The population that physically binds is the one that overlaps a probe.

**eff-lengths-do-not-cancel-at-an-end. "Effective lengths cancel, so a region's marginal is just its
length" is false near any transcript end.** A mature fragment must fit in the remaining transcript; gDNA
need not.

**configured-lengths-are-not-realised. Equal configured fragment lengths do not give equal realised
ones.** Each pool's length marginal is reweighted by its own template opportunity, and a 2 kb transcript
truncates the tail a whole chromosome does not. Gate any length axis on the realised truth, never the
config.

**mature-rna-never-crosses-a-boundary. Mature RNA never crosses an exon↔intron boundary; exon↔exon
boundaries it does cross.** Measured at zero of every boundary on every condition, this is the hard
empirical case that a contiguous boundary and a splice junction are physically different objects, and it
is what makes the two exon-crossing gDNA pools pure.

**a-boundary-with-rna-is-not-an-sj. A boundary with RNA crossing it need not be a splice junction.** One
position can be a splice donor for transcript A and plain contiguous exon for transcript B; zero-gDNA
libraries show boundaries carrying tens of thousands of unspliced fragments that are entirely RNA.

## G. Process

**no-magic-numbers. No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable.
Every divisor must be derived from the deposit rule and unit-tested against brute-force enumeration.

**one-thing-varied. One thing varied per experiment**, with the falsification test written first and
verified failing, and a baseline re-recorded from the current tree in the same session.

**a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism. A prototype harness that
subclasses the shipped policy and calls the parent's `prepare` delivers the landed mechanism too, the
moment it lands, and then double-counts its own copy.** Three innocent candidates were chased before the
double delivery was found. A landing's faithfulness check compares `src` before against `src` after, and
a harness built on the parent class is retired the day the parent gains what it prototyped.
*Sibling:* `TRAPS: one-thing-varied`.

**converge-and-delete. Converge and delete.** No legacy, no backwards compatibility, no speculative code.
Code kept "for comparison with the old version" is a defect. No version suffixes in file names.

**the-source-does-not-cite-docs. The source does not reference the docs.** Docs evolve and rot: most of
the citations that used to be in the source pointed at documents that had already been deleted. A
docstring may cite a test or an executable specification, because code cannot rot silently; it may not
cite a document.

**real-data-is-a-test-input. Real data is a test input, never a design input.** The cfRNA on disk is one
far end of the RNA-seq spectrum, not a sample of it: sweep the plausible space, report the worst case,
bring the domain call to the owner. In particular, never assume RNA fragments are longer than gDNA; true
for cfRNA, false elsewhere.

**running-an-arm-is-a-fresh-process. Six operational traps from running panel arms, each of which cost a
launch.** (i) Never edit `src/` while an arm is running: every arm imports `src` at start-up, so a
mid-flight edit changes what half the shards measured. (ii) zsh does not word-split an unquoted variable;
use an array. (iii) A wait-loop whose `pgrep` pattern matches its own wrapper deadlocks; wait on a log
marker. (iv) A default path stored in two homes goes stale in one. (v) Node-axis and region+boundary
figures differ by ~2×; say which. (vi) A composite arm fires only its components' names, so a guard keyed
on the arm's own name trips after a valid run; check why a guard fired before distrusting the data.

**shard-an-arm-sweep-by-condition. Shard a panel sweep by condition, never by arm over one condition; the
instruments write their caches.** Concurrent arms on the same condition race writers on one directory:
one `payload.npz` came out truncated and every later read raised `BadZipFile`. The damage is loud, and an
instrument that only reads should take its full payload from `scan_cache/`, not from a directory another
instrument writes.

**no-enumeration-without-a-census. Do not re-propose path or cell enumeration without a memory census.**
Possible unspliced paths are about a million regions times a few reachable ends at ~100 B each, several
hundred MB, plus spliced paths. It was killed by memory, and no consumer needs it.
