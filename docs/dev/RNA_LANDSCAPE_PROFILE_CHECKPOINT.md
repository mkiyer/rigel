# Landscape profiles: reference and next controlled experiment

*2026-10-07. Owner-authorized prototype. The numerical reference is outside the tree;
controlled calibration replays use frozen old curves. Production training is unchanged.*

## Why this is the next decision

The owner accepted the modest test-chromosome OFF tradeoff and requested larger validation.
The stranded results hold in both length-gap directions, but the factory-free candidate
regresses on RNA-long unstranded capture-ON. The owner explicitly keeps that stratum outside
the 0.8.0 release requirements while authorizing deeper robustness investigation. The accuracy tables
are in [RNA_FULL_PANEL_VALIDATION.md](RNA_FULL_PANEL_VALIDATION.md). The fixed-landscape
attribution and training census belong to
[ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice).

Most of that captured-row loss disappears when the old landscape is held fixed. Much of
the 98%-DNA OFF count loss remains. Population fitting is therefore a justified next
experiment, not a proven cure for both effects. Reusing the old factory-trained population
would reintroduce the removed assumption indirectly; it is an attribution control.

## Exactly what was implemented

`.cache/rigel_runs/2026-10-07_profile_fit/profile_fit.py` takes an object-by-density
log-likelihood array, object weights, a numerical objective tolerance and an iteration
budget. It returns fitted mass on the supplied density grid. It knows nothing about
introns, capture status, simulated truth or distant RNA abundance.

With weights normalized to sum to one, the fit is:

```text
F(p) = sum_i w_i log(sum_j p_j L_ij),   p_j >= 0, sum_j p_j = 1.

g_j = sum_i w_i L_ij / sum_k p_k L_ik
p_j <- p_j g_j.
```

The objective is concave. Since `sum_j p_j g_j = 1`, `max_j g_j - 1` bounds the remaining
objective improvement. The caller's tolerance controls numerical accuracy; the iteration
budget raises an error if convergence is not established. There is no pseudocount, density
threshold, fitted tuning constant or smoothing prior. Zero-weight rows are removed before
division, and per-row log scaling prevents overflow without changing the optimum.

This is not averaging normalized curves. A density-independent likelihood contributes a
constant to the objective, so an uninformative object cannot vote for a uniform population.
Whether the eventual data-derived profiles preserve an informative minority remains to be
measured; the formula alone is not that validation.

This is a dense small-grid reference, not a whole-genome implementation. It does not choose
grid support or regularization, calculate the input profiles, or certify local messages.
A density likelihood does not become a realized DNA count by multiplying its fitted rate
by opportunity.

## The distinction in plain language

The current admission rule asks whether the *answer after applying the prior* is precise
enough to teach that prior. An object just below the variance limit teaches it; one just
above the limit does not. The admitted object then contributes a curve built from its
estimated DNA count as though that count had been observed. A prior can consequently help
decide which observations get to correct it.

The proposed input instead asks how well each possible DNA density explains the observed
reads. Keep that whole answer: a precise observation rules out many densities; an uncertain
observation can remain compatible with several. A completely flat answer provides no
information about the population shape. Fit one population to those compatibility curves,
rather than average the curves as if each were a distribution of population members.

This removes a statistical decision rule, but does not remove all computation. The current
cutoff is shorter code; the proposed evidence contract is simpler reasoning. It needs
correct observation likelihoods, stable numerical integration and a verified population
optimizer. Existing posterior-derived weights remain fixed in the first A/B, so removing
admission alone is not a claim to have removed all prior feedback.

The population is still a shared **DNA-density** prior for counts. An object's evidence
uses its reads and justified local observations, not a distant gene's RNA expression.
Count estimates may still borrow information through that explicit DNA population, as
already authorized. The capture consumer and its locality contract remain separate.

## Tests and deliberate failures

The initial curve-average implementation fails three gates: an independent two-point
likelihood optimum, invariance to added flat rows, and the observed-DNA Poisson limit.
After replacement, a disjoint-support zero-weight example exposed division by zero; its
test failed before zero-weight rows were removed.

The initial six gates pass: the independent optimum; added uninformative objects; arbitrary
row scaling and common weight scaling; the observed-DNA Poisson-mixture limit; equivalence
between repetition and weighting; and a zero-weight object with disjoint support.

Four deliberate breakages are caught: averaging evidence, discarding weights, normalizing
columns instead of rows, and retaining zero-weight rows in the division. Every gate fails
under an appropriate mutation. Receipts are `before.log`, `zero_weight_before.log`,
`after.log`, `mutations.json` and `mutant_*.log` beside the reference.

### Extension after the admission approval

The reference now analytically cancels exactly constant likelihood rows before normalizing
weights for numerical convergence. Without this cancellation, an arbitrarily large number
of uninformative objects can dilute the stopping test even though they cannot change the
mathematical optimum. The new falsification fails before this change and passes afterward.
No approximate-flat threshold is used. With all rows constant, the reference returns a
uniform representative of an unidentified fit; this is not evidence for a uniform population.

Four added gates cover that cancellation, an exclusively supported rare population,
raw mixed-origin Poisson evidence against an independent two-point optimum, and the all-flat
limit. The likelihood reference now passes ten gates and six deliberate breakages.
The raw-observation case is a primitive test, not a proposed own-strand-only release reader.

## Keep the next contrasts separable

The previous plan said to change kernels first, but its proposed mixture fit also changes
the estimator from averaging kernels to fitting a population likelihood. Separate these:

1. Freeze the selected objects, weights, grid and rendered curves from one calibration.
   Verify that recombining those curves reproduces the old population. Compare its average
   with the likelihood fit on exactly those inputs. This prices the estimator change;
   it does not make the old curves honest evidence.
2. Replace those curves with certified observation likelihoods, holding selected objects,
   weights and grid fixed. Keep the factory absent. Do not relabel posterior-derived
   messages as likelihoods. Certify the relevant splice/level factor before using that
   channel on the unstranded examples.
3. **Owner authorized, 2026-10-07:** separately remove only the posterior-variance admission cutoff.
   Retain the composition-channel requirement, structural exclusions, zero-count anchors
   and current weights. Broad likelihoods then contribute uncertainty instead of being
   discarded because the count posterior is broad; flat curves provide no shape information.
   This is a separate experiment on DESIGN's location-admission rule, now explicitly approved.
4. If needed, price the posterior-derived reliability weights in another contrast. Do not
   combine unit weighting with admission. Admitting bound-only objects, boundaries or
   both-strand regions is a further question, not included in this request.
5. Only pair a successful partner with applying the DNA prior once. Preserve the RNA
   reference and count summary. The rejected prior-once-only arm is not a landing candidate.

First use cached inputs for the RNA-long unstranded ON failure, ladder g98 unstranded OFF,
a zero-DNA control and stranded controls. Require useful recovery without confidently
inventing DNA in controls. Check numerical support/refinement and informative-minority
retention before the full condition/class screens and per-stratum release A/B.

Do not add an unstranded factory switch, another class prior, capture detector, uncertainty
floor or panel-fitted multiplier. If the missing local information cannot be supplied under
the existing model, bring that limitation back to the owner. Population fitting must not
be credited with repairing the direct loss of local information without a measurement.

## Owner authorization

The owner approves replacing the **posterior-variance admission cutoff** with evidence-curve
contributions in the separate contrast above. Keep the composition requirement, structural
exclusions and weights fixed in the first comparison. This authorizes the prototype;
integration of the new capture reader remains a separate checkpoint.

## Implemented controls and what they establish

All new experiment files are in `.cache/rigel_runs/2026-10-07_profile_admission/`:

| File | Exact role |
|---|---|
| `admission.py` | Removes only the variance-selection line from the frozen Python fitter; an attribution control with old inputs |
| `screen.py` | Seven cached conditions, fresh current/candidate/control arms; freezes each refit's inputs and scores against slot truth by axis, class and admitted strand |
| `frozen_curves.py` | Reconstructs individual old count curves, including previous-prior weighting and population blur; retains reliability weights and existing pseudo-region |
| `fit_frozen.py` | Compares averaging with likelihood fitting of those same frozen curves |
| `optimize_zero.py` | Independent numerical check of the same zero-control objective, with a simplex-feasibility and global objective-gap certificate |
| `replay_frozen.py` | Replays the three frozen populations in calibration; mean controls reproduce original count arrays bit-for-bit |
| `test_admission.py`, `test_frozen_curves.py`, `test_evidence_contributions.py` | Admission invariants, frozen reconstruction, numerical neutrality and true observation-model limits |
| `mutate.py`, `mutate_curves.py`, `verification.json` | Deliberate failures and source/input identity checks; the likelihood mutations remain beside `profile_fit.py` |

Settled measurements, including the full seven-condition table and both estimator contrasts,
are in [ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice). Removing admission alone
barely affects the large failure. The estimator-only control repairs part of that failure,
but substantially worsens false DNA on the zero control and worsens intron counts on the
captured condition. Do not promote it or change a cutoff to conceal that tradeoff.

The combined reference/control checks comprise 22 passing gates and 16 caught deliberate
defects, with every gate exercised by a defect. These are scratch checks, not a substitute
for the standing production suite. No main source, existing tests, goldens or installed
native library changed; the owner's pre-existing count-repair diff is preserved.

## The next implementation checkpoint

The estimator comparison is now measured. Proceed to **observation inputs**, retaining the
same initial population, weights and domain for attribution. The exact own-observation
reference and one intron–boundary face are certified; splice/level factors and shared-data
accounting still need certification. Preserve uncertain RNA/DNA alternatives in those
curves rather than exporting a count median and a Gaussian width. Do not import the removed
factory or the target landscape through incoming messages.

Before a larger calibration candidate, require the raw-profile interface to pass the zero-rate,
zero-opportunity, small-count, strand-symmetry, prior-invariance, local-factor and repeated-
observation checks. Verify support extension and grid refinement independently of panel
error. The frozen old grid is an attribution device, not a justified density-rate support.
Neither more optimizer iterations nor narrow fitted population peaks repair incorrect inputs.

Then apply the authorized admission contrast with those verified curves. Price removal of
posterior reliability weights separately. Pair a successful population change with prior-once
only afterward. The unstranded ON regression remains a robustness investigation; this work
does not make completion of every evidence factor a new release requirement for the frozen
factory-free candidate. Its remaining in-scope end-to-end, real-library and landing checks
still determine production readiness.

The next input contrast and raw splice-face reference are now implemented in
[RNA_OBSERVATION_INPUT_CHECKPOINT.md](RNA_OBSERVATION_INPUT_CHECKPOINT.md). The partial
raw-input replacement restores the zero control but does not settle the mixed-object model.
Local route uncertainty is the next explicit assembly decision; the proposed profiling
experiment is specified there rather than silently adding another RNA reference measure.

**Current continuation, 2026-10-08:** the message-based density evaluator, separate-channel
correction and authorized evidence-curve admission contrast are now implemented outside
production. The admission result keeps the measured zero control stable but trades intron
improvement for worse exons and boundaries on RNA-long unstranded ON. It is not a landing
candidate. [RNA_EVIDENCE_ADMISSION_CHECKPOINT.md](RNA_EVIDENCE_ADMISSION_CHECKPOINT.md) records
the exact intervention, attribution controls and the subsequently approved and completed
equal-weight comparison. That comparison also trades class errors; no production weights
have been selected.
