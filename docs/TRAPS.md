# TRAPS — mistakes this project has already made

Read before designing anything. Each rule is a mistake already made here, stated as what to do and why. A rule is not
a measurement (where the tool stands is `ROADMAP.md`'s), and open problems live in GitHub issues, rulings and
derivations in `DESIGN.md` and `EQUATIONS.md`.

Cite a rule by its name, e.g. `TRAPS: panel-before-src`. A rule starts at column 0 as its bold name and a full stop,
and is one paragraph; no other line here starts with bold, because `tests/test_docs_boundary.py` reads every such line
as a rule. A rule names a defect, never a licence: citing one authorises no mechanism and vetoes none.

## Gates and validation

**self-checking-validator.** Validate a result by re-deriving it with a different algorithm, never through the
builder's helper or its arithmetic: such a gate cannot fail while the code is self-consistent, so it certifies the
defect and fires on the repair.

**perturb-every-gate.** After the falsification test passes, break the fixed code and watch each gate fire. A
perturbation that fires nothing names a missing gate, usually a fixture that cannot see the rule (one reference,
symmetric classes, pieces all longer than a fragment, one side of a clamp untouched): fix it and perturb again.

**could-the-arm-have-fired.** Count the opportunities a change had to fire and print the count beside its result.
Fixtures and controls are arms too: ask whether each could have failed (a structurally pinned object cannot) and holds
the channel that could overturn the claim; an oracle arm must reach every population its truth describes.

**an-ablation-that-never-ran.** Make every arm prove it ran: assert the patched name exists (patch through
`sys.modules`), count firings and raise on zero, and log what it applied rather than its flags. An inert arm reads as
"no effect", and an inert prototype as a refutation.

**a-gate-that-already-passed.** Check which gates would have failed before the fix; one that passes with the defect
present is no falsification.

**right-conditional-wrong-marginal.** When every per-observation conditional is right, only a check on the marginal
sees a wrong marginal.

**a-test-that-redefines.** Give a shared definition one home in production code, imported everywhere; a test that
re-derives it inline cannot see it drift.

**a-gate-that-reconstructs.** Gate a symmetry or invariance the code cannot fake, not a reconstructed value, which is
vacuous wherever the code decides not to use it.

**a-gate-on-the-helper-is-not-a-gate-on-the-caller.** Gate through the public entry point, on a fixture where the
branches differ, and name for every site a change touches the observable that would move and the gate that reads it.

**byte-identity-gate.** Before calling two arms identical, require equal row-key sets, count what the gate still
compares, and re-record the baseline from the current tree: an empty arm scores "identical", a stale baseline or an
expired premise blames whatever is in flight, and a gate comparing almost nothing is disarmed.

**a-clip-hides-a-scale-error.** Perturb both sides of every clip and treat an exact algebraic 2 as a bug: a clip hides
errors on both sides, and a range that "agrees" with a clamp is the clamp. A zero control proves no discrimination.

**a-constant-parked-a-value-off-a-knife-edge.** When removing a constant looks like a regression, bisect it: a finite
jump from an infinitesimal move is a downstream discontinuity the constant hid.

**the-divergence-was-a-barrier.** Before crediting the removal of an infinity, ask what it suppressed: an infinite
variance in a damping term is a structural gate, and a finite one can turn a barrier into a conduit.

**a-green-suite-hid-five-dead-instruments.** After a `src/` deletion, a rename or a shipped-default flip, run the
instruments, not only the suite (`preflight.py --full`): an instrument lives only under the configuration it runs in.
Re-derive a rebuilt quantity by the producer's own expression; one ULP breaks an identity.

**the-deliverable-is-not-reproducible-by-default.** Pin `scan.total_threads=1` wherever two runs must agree bit for
bit: the scan's workers sum float64 banks in finishing order, so above one scan thread the tally's last bits move and
reach the EM. The EM's seed is fixed (`EMConfig.seed = 0`); `em.n_threads=1` pins nothing.

**off-grid-message-mode.** Keep every message mode inside its grid's domain: an out-of-range mode pins the boundary, the
most confident statement a channel can make. Capture every delivered channel in the debug output, and build a message
in the destination's coordinates, since a raw count tilt inverts sign on an antisense protocol.

**honesty-metrics-reward-ignorance.** Quote a fixed-denominator accuracy (the mass-weighted error over all objects, and
Σ|err|) beside every honesty metric: each improves as the solver stops knowing anything.

**the-intermediate-is-not-the-deliverable.** Report the upstream and the shipped column on every arm and never rank on
the upstream one: a fitted stage downstream can learn the upstream bias and reverse the win.

**an-identity-with-a-qualifier.** Price the complement of every identity that holds "over", "for", "assuming" or
"where" something; the qualifier names a population nobody has sized.

**a-flat-pile-is-not-a-knot.** Before merging modules, ask whether the problem is entanglement or missing order:
naming the layers already present finds the upward imports, and merging a flat pile only makes bigger files.

## Measurement and inference

**score-against-truth.** Score against truth (the per-fragment truth in the oracle BAM's read names), never the
previous run, and raise on a truth label you do not recognise: nested aggregates double-count when summed.

**re-record-the-baseline.** Vary one thing per experiment and re-record its baseline from the current tree in the same
session; a delta against a stored baseline is not attributable.

**measure-the-ceiling-first.** Before building a correction, hand the consumer the exact answer for that channel and
see what perfecting it is worth (`quant_accuracy.py`'s injection arms); a far better input is worth nothing to a
consumer that reads it weakly. An A/B says whether a change helped; a ceiling says whether to start.

**a-symptom-is-not-a-second-defect.** Before calling a downstream number a second bug, substitute only its upstream
input and re-run the shipped downstream code; a ceiling prices only what is built after its injection point, so write
that point down. When each input alone moves an error little and all together move it fully, go one stage upstream.

**attribution-must-survive-a-shuffle.** Re-run every truth-fed arm with the truth shuffled; an oracle a scrambled
oracle also beats has not priced the truth.

**price-the-halves-separately.** Price each half of a diagnosis on its own; the fully diagnosed half need not be the
valuable one.

**never-pool-the-strata.** Report every score per stratum with its denominator named: a pooled mean buries unsolvable
objects, hides a sign flip between strata and lets the deferred stratum, which carries most of the error, set the
number. A score over a mechanism's own targets is a diagnostic; only all live slots is a verdict.

**excluding-a-population-hides-it.** Gate every population excluded from, or admitted to, a scored denominator, report
the solvable fraction beside every accuracy, and never cut on a fitted residue; a class ratio needs objects that can
supply its numerator.

**zero-target-guards-are-one-sided.** Use zero-gDNA rows as false-positive checks only: at truth zero, any change that
lowers the gDNA estimate scores better.

**hard-labels-miss-soft-change.** Report Σ|err| and the directional split beside every net error: an identical
hard-label result is no evidence, and a net only bounds the error from below.

**a-single-level-panel-cannot-see-a-constant.** Check the panel varies what the deliverable estimates, at three levels
at least, and substitute a candidate by its own class mean: if that costs nothing, the candidate is a constant.

**score-the-consumers-own-count.** Write down which line of the consumer reads a number and which fragments it counts,
score against that set, and apply the same selection to the noise floor's draws.

**weight-it-like-the-consumer.** Weight an average as the consumer weights it; a bp-weighted and a fragment-weighted
mean answer different questions.

**draining-breaks-the-oracle.** Take a per-fragment truth partition only before a step that conditions on the whole
tally: the second pass reads the payload's own densities, so parts drained separately are not the whole drained.

**an-equal-length-panel-defeats-the-lift.** On an equal-length panel run both sides undrained and price the drained
lift's span-tie ambiguity, which equal lengths make the common case.

**a-length-gap-bypasses-calibration.** On a length-gap panel score calibration against oracle calibration beside the
end-to-end number: the EM can assign by length alone, so a good end-to-end number can sit over broken calibration.

**judge-a-ruler-by-its-within-gene-spread.** Judge anything the EM divides by on its class scale and the within-gene
spread of its log error (`ruler_vs_truth.py --scale`), against the yield the reads were drawn with: the EM reads
ratios inside a locus, so a length closer to truth per transcript can still split isoforms worse.

**a-profile-share-is-a-ranking.** Use a profiler's share to choose what to open and an interleaved A/B pair to price
it: cProfile inflates many-call functions, and a toy ranks hotspots backwards, so profile on a deep real library.

**panel-before-src.** Run the panel arm before writing a mechanism into `src/`, and name the substrate of every claim:
a toy can invert a ranking, since its objects lack their own evidence and a message is nearly free there.

**prove-the-substrate.** When a simulated axis is the axis being judged, gate the simulator on it, and when an error
concentrates on one simulated mechanism, A/B the simulator's assumption before engineering against it. Read the
sampler's weight for the population under study, since a substrate knob can suppress it.

**can-the-benchmark-resolve-it.** Prove a benchmark resolves the axis being changed: score each requirement against
its degenerate value, which is what a structurally blind suite scores.

**a-broad-population-carries-no-prior.** Before fitting a population prior, measure the pooled width against the
per-object uncertainty; if it is wider, there is no prior to be had.

## Pools, selections and divisors

**a-cancellation-is-conditional-on-its-support.** Write a reciprocal-opportunity expectation with its support factor,
`E[Σ 1/A] = ρ · P(A > 0)`, and gate it with fragments placed outside the support, where nothing deposits.

**a-purity-filter-is-a-length-filter.** Before excluding fragments from a length pool, ask what the criterion
correlates with: barring the partly inferred lengths bars the long fragments. Key a pool on determinacy.

**pure-and-length-censored.** Ask what a pool's selection correlates with, not only what it selects: contained
intergenic or intronic gDNA is composition-pure and, beside a probe, length-censored.

**divide-by-a-probability.** Divide a pool by its opportunity as a probability, `A(w)/T(w)`, never a count:
`count(w)/A(w)` recovers the drawn distribution, and every consumer needs the realised one.

**opposite-tilts-must-not-pool.** Never apply one divisor to summed pools of opposite opportunity tilt; sum the counts
and the matching per-pool opportunities, which is the inverse-variance average of the per-pool estimates.

**two-divisors-opposite-sign.** Differentiate two divisors built from one pmf before hunting a bug in either; they can
respond with opposite sign. Read the simulator's sampling code, not its docstring, before pricing a selection.

**a-pooled-rate-cannot-see-a-short-object-factor.** Compute the pooled weight of the objects a factor acts on before
building a panel to exert it: a pooled rate is dominated by the largest objects, and blind when applied per component.

**a-mean-of-ratios-inherits-the-partition.** Use `Σmass / Σopportunity`, never a mean of per-object ratios, and test
invariance to a re-partition that leaves the data unchanged.

**a-zero-count-is-a-measurement.** A zero count measures a density with finite precision: under Jeffreys, `a` events
over exposure `E` give `Gamma(a+½, E)`, so zero over a megabase says far more than zero over 200 bp.

**a-ratio-cannot-carry-zero.** Before fixing a source that will not speak, check its channel can carry the value: a
multiplicative transport cannot say "there is none here".

**frame-free-is-not-assumption-free.** When a divisor cancels, examine the nuisance profiled out: a frame-free term can
still be flat in the parameter, and a channel whose benefit tracks the answer, not the evidence, is a prior.

**an-upper-bound-is-not-an-estimate.** Never read a shared bound as evidence about who shares it: a transcript inherits
its path's minimum density from objects it does not own.

**fractional-mass-is-the-problem.** Keep multimapper and ambiguous-path assignment integral: a conserved fractional
mass makes a fragment's deposits depend on region sizes, where its crossing counts depend on nothing.

**conservation-misses-mis-attribution.** Gate a deposit rule per object, on an axis the rejected rules cannot express;
mass conservation holds under every alternative.

**a-fold-grows-a-heuristic.** Give a consumer the axis it lacks rather than folding a quantity onto one it has; the fold
grows a repair heuristic that outlives its reason.

## Estimation and solver design

**purity-is-a-property-of-the-annotation.** Never assume a structural class is pure gDNA ("intergenic" is whatever the
annotation leaves over), and resolve RNA-vs-gDNA by likelihood, never membership: gDNA is compatible with everything
unspliced. Use the half of a statistic a contaminant cannot reach, and check nothing pushes the other way.

**one-hop-lifted-out-is-still-the-relay.** A value computed at one object and consumed at another is a message: put it
in the policy under `calibration/messages/` and price it there. Reasoning about how an exon gets its gDNA level (a
pooled reference, an imputation, a bound) is re-deriving the message passing that ships as `transfer`.

**a-message-from-the-destinations-belief.** Build a message from the destination's constants and observations, never
its belief, which it would only confirm; transport by the claim's own component measured at both ends, never by a
total or the source's belief.

**recompute-from-the-oracle.** Before correcting a formula, substitute its definitions and recompute it from the
oracle; an impossible-looking ratio may be right.

**variance-fitted-on-the-belief.** Never fit a variance, or any component, on the solver's own unsolved belief, and
never repair a bias with a variance: neighbours that share an error agree on it, and the confidence propagates it.

**two-gaussians-one-latent.** Never hand a solver two Gaussians built from one latent: `log f` and `log(1−f)` are
correlated exactly −1, so their summed Fisher information is over-confident.

**an-imputation-must-cost-something-every-hop.** Charge every hop a variance for its premise that does not shrink with
depth, and split a message layer's error by whether the destination had its own evidence.

**a-licence-with-no-floor.** Size a licence, or a predicate gating a multiplication, by the density behind it, never
by a precision being positive, and zero a refused claim's precision with its value: a confident zero inverts one
rescale on.

**a-cancelling-defect-pair.** Price a fix to one of two cancelling errors in the arm that removes the other too; one
alone is worse than neither.

**no-prior-means-haldane.** Write a prior for every grid term: omitting one lets the grid supply Haldane, `p(x) ∝ 1/x`,
an improper amplifier toward the vertices.

**a-priors-curvature-is-not-the-datas-information.** Add to a Fisher precision only terms that scale with the data; a
prior's curvature belongs in the posterior.

**an-all-zero-factor-is-inert.** Skip a factor that is zero for every candidate; multiplying by it makes the record a
coin toss.

**density-below-one-fragment-length.** An object shorter than a fragment has no measurable density: a component with
no opportunity emits nothing at zero precision, and "no data" is inert, never "100 % gDNA".

**identical-paralogs-are-bimodal.** Assert an identical-paralog split structurally (even, or `min == 0` and
`max == total`) and score only what is identifiable; never move a seed to make it land even.

**a-clamp-at-the-closed-end-escapes-the-window.** Test a prior at the value the code uses: clamped at the closed end of
its support, it can hold almost no mass inside the solve window while the interior looks fine.

**deriving-one-coordinate-propagates-its-error.** Before deriving coordinate B from A, list what B was exact about and
check A is too, and price A against truth per stratum; a derived quantity inherits its source's error, amplified.

**interpolate-on-the-axis-where-the-lattice-is-uniform.** Do quantile arithmetic on the grid's own axis (`λ`), then
map; interpolating on a non-uniform transform (`f_g`) biases toward the middle.

**a-constant-in-exact-arithmetic-is-not-constant-in-float64.** Where a term vanishes in ℝ, return the vanishing value
and gate it bit for bit (`solve(neutral) == solve(None)`); `repr` any number a diagnosis turns on, since a printed
`0.0000` need not be zero.

**a-face-total-is-not-a-total-without-its-flux.** Add a face's junction flux to its contained and crossing abundance
before comparing totals across a boundary; spliced RNA crosses as flux, so a raw total reads structure as enrichment.

**a-mode-count-is-not-a-well-posed-quantity.** Let no consumer read a shape statistic off a fitted density until a
resolution sweep shows it converges, and measure the share of kernels at the floor (`landscape.knn_widths`) before
reading a bandwidth sweep.

**a-linear-likelihood-emits-a-sign.** Summarise a likelihood asymptotically linear in its parameter as a location plus
an information, never a row handed to a normaliser: it has only a direction, which saturates at a grid endpoint.

## Structure, indexes and plumbing

**one-reference-hides-refid-bugs.** Build test indexes with at least two genomic references; one reference hides a
reference-id mismatch, and RNA-only spike-ins do nothing for the gDNA path.

**annotated-is-not-genomic.** Take a classification as an input and write the filter that defines a population,
raising when either is mis-stated: "has an annotation" is not "is genomic", and a file's rows are no definition.

**nrna-does-not-mean-synthetic.** Filter real transcripts with `~is_synthetic` alone: `is_single_exon` (once
`is_nrna`) means single-exon, not manufactured.

**an-sj-is-not-a-gap.** Pass a junction's cut-intron coordinates through and deposit at them; a junction is not a gap
between deposited slices, and a region's boundary can lie kilobases from an unannotated junction.

**splicing-makes-the-graph-cyclic.** Treat the region↔junction graph as cyclic and never break a cycle by dropping a
junction boundary, which re-isolates the exon it exists for.

**credit-exactly-one-sj.** Credit a fragment crossing K junctions to one, the leftmost annotated (all K shifts the
sense fraction and fakes overdispersion; `1/K` is biased), and key junctions on `(src, kind, dst, strand)`.

**a-transcript-predicate-must-not-silently-drop-a-molecule.** Check genomic compatibility before a transcript-level
predicate judges a fragment, and close the ledger (`n_fragments == deposited + Σ accounted drops`): a rejected
population is never random.

**an-object-class-does-not-see-a-terminus.** Read a hop type as `object class × {sj, term, sj+term}` from
`RegionStatics.boundary_flags`; a terminus and a splice site carry opposite currencies.

**two-masks-one-name.** Give each predicate its own name and each quantity one described home; two predicates under
one name, or two docstrings of one quantity, diverge silently.

**a-hash-that-misses-its-artifact.** Make a cache key cover what it caches: compute hashes on demand, check nothing
excluded from a key reaches the blob, and make an overridable input's derived value a `@property`, not a field.

**integer-channels-reproduce.** Keep counts integral, so they reproduce across worker counts, and hold float64
fractions to a derived tolerance; sort an order-observable bank on its own content before it crosses the ABI.

**a-kernel-edit-moves-ulps-under-fma-contraction.** Compare arms on one binary at one thread count against that
binary's shipped arm: `_em_impl` builds with `-ffp-contract=fast`, so an unrelated edit moves ulps the EM amplifies.
`--threads` pins every thread count; `--set scan.total_threads=1` pins only the scan.

**worktrees-run-the-wrong-code.** Import a prototype build first in every process, shards included (a
`sitecustomize.py` on `PYTHONPATH` reaches them): the editable install beats `PYTHONPATH`. Retire a harness built on
the shipped class the day that class gains what it prototyped, or it delivers the mechanism twice.

**a-rebuild-can-silently-no-op.** `touch` the sources before a rebuild and check the artifact, not the exit code
(`pip install -e` keys on mtime); restore a perturbation from a copy, since `git checkout` deletes uncommitted work.

## Domain facts that read like defects

**specificity-and-sense-are-complements.** Compare a simulator's `strand_specificity` with
`StrandModel.strand_specificity`, never a fitted sense fraction: on an R1-antisense protocol a κ near 0 is correct.

**strand-measures-the-tilt.** Strand measures the RNA tilt `d`, not the gDNA fraction: in `p = ½ + (κ−½)·d` it
cancels, and reaches strand only through `f_g ≤ 1 − |d|`.

**equal-lengths-carry-no-composition.** At equal component mean lengths the length channel carries no composition
information at any depth, which is why the ladder equalises them: the EM must exercise calibration.

**configured-lengths-are-not-realised.** Gate a length axis, and key an operating point, on the realised truth and
print it, never the config: each pool's lengths are reweighted by its template's opportunity, and capture moves gDNA.

**capture-selects-for-length.** Score against post-capture truth, for lengths as for abundances; the pre-capture
parameters describe a library that was never sequenced.

**capture-is-1000x-on-exons.** Read capture from gDNA, whose uniform baseline shows the probe pattern RNA's expression
hides; capture also erases the intron signal, so nascent RNA and gDNA blur under it.

**capture-inverts-the-counted-side.** Never fix which side of a transport is the well-counted one; capture inverts it.

**on-target-by-start-is-geometry.** Define on-target by overlap with a probe, not the start's territory: an
intronic-start fragment must be long enough to reach the probe, so it reads longer under any capture model.

**eff-lengths-do-not-cancel-at-an-end.** Near a transcript end an RNA fragment must fit in what remains and gDNA need
not, so effective lengths do not cancel there.

**mature-rna-never-crosses-a-boundary.** RNA spliced at an intron never crosses that intron's exon↔intron boundaries
contiguously, while it crosses exon↔exon ones, and one position can be a donor for one transcript and plain exon for
another.

## Process

**no-magic-numbers.** Stop and discuss before adding any constant, heuristic or tunable; derive every divisor from the
deposit rule and unit-test it against brute-force enumeration.

**real-data-is-a-test-input.** Use real data as a test input, never a design input: the cfRNA on disk is one far end
of the spectrum, so sweep the plausible space, report the worst case, and never assume RNA is longer than gDNA.

**read-the-whole-failure-list.** Derive the failure set, never eyeball the tail:
`pytest -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c`.

**running-an-arm-is-a-fresh-process.** Never edit `src/` while an arm runs; shard a sweep by condition, never by arm,
since the instruments write their caches; pass shell lists as zsh arrays; and wait on a log marker, not a `pgrep`.
