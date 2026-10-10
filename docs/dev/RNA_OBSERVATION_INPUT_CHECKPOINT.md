# Observation inputs: implemented limits and the splice assembly checkpoint

*2026-10-07. Development checkpoint. Everything implemented here is outside production.
The owner authorized continuing the evidence-input work; no capture reader or prior-once
change has been integrated.*

## What was implemented

The first input contrast replaces only the curves for objects whose DNA count is structurally
known and objects with zero observations. Their own evidence has the exact Poisson form
`count log(rho Eg) - rho Eg`, with the zero-rate limit handled analytically. Integrating RNA
at an empty object changes only a density-independent constant. The replacement reads the
observed total, not a fitted DNA count, and takes no previous population or posterior variance.

This is deliberately partial. All other mixed-object curves, the selected population,
weights, density grid, fitted objective, existing pseudo-region and count solver are fixed.
It is not an own-strand-only capture reader. The two separate deletion controls retain either
the old population blur or the old previous-prior weighting, identifying their interaction.
Settled measurements have one home in
[ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice).

The same likelihood objective is solved with an independent constrained optimizer. A finite
objective, feasible population and global objective-gap bound certify its result; its generic
success flag does not. A deliberate zero-likelihood output exposed a missing finite-objective
guard, now fixed. Numerical tolerances are stated in each receipt and do not depend on panel
accuracy. The new reference is not the proposed whole-genome solver.

Files in `.cache/rigel_runs/2026-10-07_observation_inputs/`:

| File | Role |
|---|---|
| `limits.py` | Exact known-DNA and empty-object input limits |
| `fit_inputs.py` | Frozen population fits and separate blur/previous-prior contrasts |
| `replay.py` | Cached calibrations; old-curve controls reproduce their archived count arrays bit-for-bit |
| `splice_factor.py` | Raw conditional observations at one single-strand splice face |
| `test_limits.py`, `test_optimizer.py` | Independent Poisson/integration limits, untouched mixed rows and numerical certification |
| `test_splice_factor.py`, `test_splice_identification.py` | Poisson-conditioning identity, route opportunities, invariances and the unstranded identification family |
| `mutate.py` | Deliberate faults, including invented zero support, wrong component units and fixed unmeasured route shares |

The initial unchanged-input implementation fails seven of fourteen input gates. The provisional
equal-junction-opportunity approximation fails seven of sixteen splice gates. The corrected
components and optimizer checks total forty passing cases; thirteen deliberate faults exercise
every gate. These are scratch checks, not the production suite or release validation.

## The result and its limit

Raw observation limits remove the estimator-only arm's zero-control blow-up. The two partial
deletions do not; neither a smoothed curve nor one multiplied by the previous prior becomes
a likelihood by itself. The difficult captured condition remains close overall but loses
more intron accuracy, so the partial replacement is not a landing candidate.

This result supports the evidence contract, not a claim that the remaining population,
grid or mixed-object model is correct. Continue to read classes and strands apart. The
deferred unstranded ON comparison remains a robustness test, not a 0.8.0 requirement.

## What the splice factor adds

The source observations are the two unspliced columns and each junction's count, all at one
boundary. These are disjoint banks at that coordinate. Each junction retains its own
opportunity. Conditioning on their combined total cancels a common source abundance/capture
multiplier. The source factor remains inside the target RNA integral.

The exact derivation is in [EQUATIONS.md](../EQUATIONS.md), **Conditional observations at
one splice face**. It retains the existing same-capture premise at the face; it does not
solve the component-law frame mismatch or arbitrary route-specific capture.

The remaining unknown is local RNA routing: what share continued through this boundary,
and what share used each splice route. The current point-rate map hides this uncertainty
inside its measured rate and approximate width. A raw observation likelihood exposes it.
Simply fixing the shares at their count-based estimate would repeat the plug-in problem.

The unstranded reference makes the limit explicit: multiple DNA densities, including zero,
can explain the same source counts by changing the RNA route shares. A splice count certifies
RNA; it does not, by itself, certify positive DNA. Pure-DNA neighbors and the retained DNA
population prior remain possible information sources. No intron-background constraint is
silently restored.

## Authorized experiment: profile the local route uncertainty

The owner approved this experiment, emphasizing simplicity and robustness to changing
probe placement, fragment lengths and capture efficiency before panel scores. The scalar
one-route implementation and independent arithmetic gates are now complete. Its separate
generality condition fails: unequal capture averages across fragment banks can move the
source factor's preferred DNA density even on noiseless counts. The exact implementation
and proposed next step are in [RNA_ROUTE_PROFILE_CHECKPOINT.md](RNA_ROUTE_PROFILE_CHECKPOINT.md).
Settled derivations and findings have moved to EQUATIONS and ISSUES, respectively.

The original proposed assembly was:

For each candidate DNA density `rho` and target RNA amount `r`, maximize the conditional
source likelihood over the permitted route shares `z`, then perform the already specified
target RNA integral:

```text
H_profile(rho,r) = max_{z >= 0, sum z <= 1} H(rho,r,z)
L(rho) = integral own_observations(rho,r) H_profile(rho,r) r^(-1/2) dr.
```

This adds no route prior, distant RNA information, threshold or capture detector.
It keeps the target RNA reference explicit. It is an integrated profile likelihood, not
a claim to have marginalized a fully specified joint Bayesian model. It can be permissive
when route shares are weakly identified; this needs measurement rather than a correction
factor. No implementation of this profiling policy has been connected to calibration.

The local maximization is small and concave after a change of coordinates. At fixed `rho,r`,
each vertex assigns all local RNA to either the continuing route or one junction. Normalize
each vertex's source category means. Every allowed source probability vector lies in their
convex hull; fitting its mixture weights maximizes the multinomial likelihood. The existing
small-grid likelihood reference can supply an independent implementation, without a new
optimization framework or a population fit inside the capture reader.

The predeclared stop occurs before the RNA integral, multiple-route implementation,
population replay and panel scoring. Integrating the same incorrect opportunity assumption
instead of profiling it would not repair the observation model. Keep the tested arithmetic
as a diagnostic reference. Do not connect it to calibration pending the opportunity checkpoint.

## Production boundaries

Main source, existing tests, goldens and the installed native build are unchanged. Existing
uncommitted owner work is preserved. No commit or push was made. The full production suite,
preflight and release A/B have not been claimed for this partial prototype. Reader integration
and a successful prior-once partner still need their specified checkpoints.
