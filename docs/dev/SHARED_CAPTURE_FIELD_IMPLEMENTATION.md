# Shared capture field implementation campaign

Sandbox; 2026-10-11. The owner explicitly requested implementation after the follow-up review.
Baseline: `624dfade4afce37d795890ea8ab537086ea20a93`. No commits. Native worktree:
`<worktree>`; harness and evidence:
`<scratch>/rigel-shared-campaign-20261011`. Production remains the comparator.

## Current status — 2026-10-12

The isolated native prototype is implemented and truth-screened. It is **not shippable**: the current
samplers/optimizers have unresolved cross-start stability, and every tested g00 model regresses transcript
accuracy. Full shared-field ladder/gap/VCaP decisions remain outstanding. The table gives transcript error
percent, with gene error in parentheses, on the **test chromosome**, not the ladder.

| Test-chromosome g00 stratum | Shipped | Dependent field | Exposure rate | Free junction |
|---|---:|---:|---:|---:|
| Stranded OFF | 7.6610 (0.8350) | 10.9547 (0.8494) | 11.4295 (0.8475) | 20.8163 (7.9381) |
| Stranded ON | 8.2243 (0.6743) | 14.2810 (0.0498) | 15.4507 (0.0398) | 40.5560 (24.5732) |
| Unstranded OFF | 5.6250 (1.0137) | 14.9057 (4.9256) | 15.0008 (4.7157) | 25.5974 (11.8341) |
| Unstranded ON (deferred) | 9.3964 (2.4449) | 14.6496 (1.0378) | 15.5310 (1.0249) | 43.0029 (26.8839) |

Durable source, gates, receipts and reproduction notes:
`.cache/rigel_runs/2026-10-11_shared_field_implementation/README.md`.
The chronological derivations below deliberately retain refuted approaches; they are not current defaults.

## Derivation before implementation

For RNA, retain the scorer's chosen alignment and its conserved transcript-path deposit vector.
For gDNA, retain every eligible alignment with its existing likelihood-mixture weight. Let `w_Lo`
be the current region/boundary membership of object o in locus L. Define the gDNA coefficients as
`B_fgo = sum_h phi_fgh d_fgh,o w_Lo`, and its opportunity as
`A_go = w_Lo sum_(all genomic paths f,h) phi_fgh d_fgh,o`. RNA uses its transcript-path equivalent.
Then `sum_f B_fko = A_ko` and `K_fk(c) = sum_o B_fko c_o`, `L_k(c) = sum_o A_ko c_o` imply
`sum_f K_fk/L_k = 1`. All of a fragment's candidate rates enter one mixture likelihood.

This specifies the gDNA component by its deposit-object membership, exactly the domain of today's
length projection. A crossing fragment can have gDNA membership below one. Renormalizing its retained
shares would violate the identity. Consequently equality with the shipped flat-field scorer is required
for complete-member paths; incomplete-member paths are an explicit, measurable difference.

Use the scorer's **observed** gDNA length law for both B and A. This is a fixed empirical base measure,
not a claim to estimate an unobserved pre-capture law. RNA likewise keeps its observed input law.
This resolves algebraic normalization but does not resolve capture's selection of the empirical length
laws. That statistical limitation is a falsifier to examine on both length-gap directions. A separate
fixed-field arm will expose the measure change before attributing effects to field integration.
The finite path identity concerns the geometric fragment model; current alignment-error penalties,
likelihood pruning and compatibility approximations remain conditional scoring inputs, not a newly
derived full sequencing-error likelihood.

The free coordinates are the existing region and genomic boundary objects. Junction values initially
remain the current homogeneous clipped boundary/intron relation. Given assignments to component and
deposit object, including deterministic RNA paths, update each free coordinate from

    log p(c_o | z, c_-o) = log p_o(c_o)
        + sum_j m_j log c_j(c) - sum_k N_k log L_k(c) + constant.

Here j includes dependent junctions; `m_j` and `N_k` are INTEGER sampled sufficient statistics.
Terms with zero count contribute zero even at zero field. Each coordinate sees updated preceding
coordinates. No fragment is revisited inside this conditional. No reader posterior or delivered
message factor enters the field prior. The existing fitted landscape supplies grid probability masses;
a common log-unit translation changes both K and L and leaves those masses unchanged.

For abundance, maximize the same conditional-mixture likelihood with the retained aggregate penalty
`P_g log a_g + P_R log(sum_R a_k)`. At expected counts N its exact MAP M-step is
`a_g=(N_g+P_g)/T` and `a_k=(N_R+P_R) N_k/(T N_R)`, with `T=N+P_g+P_R`.
This is the shipped grouped MAP rule. The shipped VB update is a separate baseline, not relabeled as
exact inference for this objective. Use stochastic approximation EM with conditional responsibility
counts: `s_t=s_(t-1)+(N_RB(c_t,a_(t-1))-s_(t-1))/t`, followed by that M-step. The gain satisfies
`sum 1/t=infinity`, `sum 1/t^2<infinity`; field transitions target the current abundance's posterior.
This is a numerical route to a stationary point, not a guarantee of a global optimum or useful mixing.
Freeze abundance for final posterior sampling and average conditional responsibilities for output.
Different starts, doubled computational budgets, an independently sampled fixed-point residual and
exact finite-mixture EM are necessary diagnostics. No SQUAREM across a moving field objective.

## Falsifiers and budget ledger

Before a BAM arm: missing candidate paths, lost deterministic paths, duplicate multimapper observations,
normalized-away partial membership, inconsistent PMFs, stale coordinate lengths, substituted mean
assignments, unit-dependent priors, unsupported observations, and disagreement with exact finite inference.
Check the defects by mutation. Existing sampling counts are computational budgets, never fitted model
constants. Keep calibration unchanged and report its oracle scores beside allocation scores.

Initial enumerable numerical checks: two minutes, under 1 GB, one process. Native prototype build:
budget three minutes, under 4 GB. No timings claimed from compilation or concurrent accuracy work.
Subsequent run budgets and receipts are appended before advancing.

## First numerical correction, before the next prototype

Native collapsed posterior joint-state error 0.00274566 at 32,768 retained draws; stale-coordinate
and omitted-denominator mutations give 0.0753775 and 0.500024. The initial stochastic-approximation
schedule is too slow at a boundary optimum: fixed-field abundance error 0.00511191 after 32,768 steps.
On the integrated toy four starts miss exact abundance by 0.00181–0.00397. These are not convergence.

Replace the schedule by ordinary Monte Carlo EM: at outer iteration t, hold `a_t` fixed for a block
of collapsed Gibbs sweeps, average its conditional responsibility counts, then apply the grouped MAP
M-step to obtain `a_(t+1)`. Increasing the block length tends to exact integrated EM; finite blocks
do not promise monotone ascent. Different seeds and doubled block lengths measure the integration
error. At a point field the block is exactly one deterministic E-step, so the implementation must
reduce to conventional grouped MAP EM without the slow averaging bias. Keep the initial burn and
a final fixed-abundance burn/sample explicit in run receipts. Block counts are numerical budgets,
not additions to the statistical model. Repeat the finite audit before any BAM arm.

First native audit actual wall 1.34 s, excluding build; this is a budget receipt, not a speed comparison.
Rebuild plus expanded numerical gates budget three minutes, under 4 GB. Main inference untouched.

## Overhang-path support, derived before repair

The first BAM support gate fired on retained RNA overhang candidates: e.g. B2_chain_ab_U,
genomic start 2844797, transcript start 2845000, candidate width 209, coefficient 2.20e-194.
Clipping this to a six-base path invents an observation on an object with zero contained opportunity
under the 209-base kernel. Dropping small weights would hide the modeling defect behind a threshold.

The scorer already permits a misaligned/overhanging observation. When its proposed path lies outside
the template, retain its existing component likelihood but treat its position as unobserved. Marginalize
that missing path with the component's own baseline placement measure: `d_err,o=A_ko/L_k(1)`.
Then `K_err,k(c)=phi_err,k L_k(c)/L_k(1)` and `K_err,k(c)/L_k(c)=phi_err,k/L_k(1)` exactly.
It supplies component evidence but no spatial capture evidence. This uses neither a capture detector
nor a new prior/constant; it explicitly states the missing-location assumption for the existing error
channel. It is not a newly normalized full alignment-error model: the previously disclosed scorer
approximations remain. Fully contained candidate paths retain their actual shares. Test the cancellation
and its zero-opportunity support before reconnecting the arm. Record how many paths use it.

Pipeline hook integrity also needed a repair: Python silently prints a sitecustomize import error and
continues. The first attempted arm therefore ran the original allocator. It is archived as
`pilot_hook_failure`, never counted as an arm. The isolated loader now fails closed on any hook error.
Routing regression gates: 110 passed. Native path checks: 8,515 placements, max response error 3.55e-15;
packing verifies selected RNA, all gDNA hits, deterministic paths, and partial-member coefficients.

Next build, support gate and pilot: three minutes, under 4 GB per process, sequential.

## Continuous conditional reference, before implementation

The complete pipeline now runs. Fresh test g50/ss.99/OFF transcript errors: shipped 9.236288%,
flat-field 6.293893%, short uniform-start field 25.703244%. Increasing to 128 outer blocks of 32
sweeps gives 16.401991%; using the shipped coverage/grouped-prior initialization gives 14.704188%.
Neither is accepted as converged or an improvement. The wrapper initially used a uniform initialization
for every arm; it now honors coverage versus uniform and records it. The largest fixed-point residual
is still about 0.00180. Native path retention alone reproduces shipped output exactly on this condition.

A fixed 0.2-nat conditional grid is itself an approximation requiring a gate. For two objects of one
component with equal opportunity and counts N1,N2, the likelihood in `d=log(c1/c2)` is binomial logistic;
its local standard deviation is `sqrt(1/N1+1/N2)`. At 5,400/4,600 counts it is about 0.0201, far narrower
than that grid. More samples cannot recover an off-grid conditional. Test against numerical integration.

Use the existing landscape as a continuous piecewise-linear **log density in log field**, on its own
finite support. This avoids selecting a new model prior or an unobservable length law. Given sampled
integer assignments, the same collapsed coordinate log density is evaluated at any x=log(c_o):

    h_o(x)=log p_o(x)+sum_j m_j log c_j(c_-o,exp(x))-sum_k N_k log L_k(c_-o,exp(x)).

Draw a slice level `h_o(x_current)+log(U)`, start the interval at the prior's full bounded support,
and sample/shrink that interval about the current coordinate until accepted. There is no width or
step-size parameter. The full-range construction and reversible shrinkage are the bounded-support
case of Neal, *Slice sampling* (2003), sections 4.1–4.3,
https://glizen.com/radfordneal/slice.abstract.html .
Unlike an independence draw, this is a Markov transition preserving the exact conditional; mixing
still needs measurement. Sequential updates always use the latest field and lengths.

The common-scale orbit also has a continuous conditional. Given relative log fields x_o, allowable
shift s is the intersection of every object's `[lo-x_o,hi-x_o]`, with density proportional to
`exp(sum_o log p_o(x_o+s))`. The likelihood cancels by homogeneity, including the retained junction
relation. The same bounded slice transition samples this direction without a new tuning constant.

Retain the discrete kernel as an enumerable reference. Validate continuous responsibilities/field
ratios against quadrature, including the concentrated off-grid case, and preserve the earlier results.
Budget three minutes and under 4 GB for native rebuild and gates. No whole-genome run yet.

Continuous gates passed: dense two-object quadrature probability 0.539973956; discrete estimate
0.549833997 (error 0.00986004, misleading MCSE about 1e-17); continuous estimate 0.539955194
(error 0.00001876, batch MCSE 0.00003932). Broad two-component count error 0.000869; omitted
normalizer mutation 0.033207. Test g50/ss.99/OFF improves to 9.434294% at 128 outer blocks and
7.231464% at the existing 1,000-iteration limit, versus shipped 9.236288%. Gene error at 1,000
is 1.249322% versus 1.223986%. Largest abundance residual 0.000816; no convergence claim yet.
The 1,000-block run took 90.13 s and peaked at 896,416 KiB aggregate RSS; budget receipt only.

## Exact grouping and feasible dispersed starts, before implementation

Repeated observations with identical ordered component coefficients have identical responsibilities.
Replace their repeated categorical assignments by one multinomial draw with integer multiplicity m.
Its sufficient statistics have exactly the same distribution as m independent assignments; use
`m*r` for conditional expected component counts and expand the averaged probabilities back to each
original observation for delivery. Group only bit-identical coefficient vectors; no rounding, length
bins or approximation. This is the already-enumerated grouping identity, now applied in native code.
The sequence of random draws changes, so stochastic outputs need reference/convergence checks,
not a false promise of bit identity. Fixed-field probabilities must agree numerically.

A dispersed prior start can give a junction zero weight and an observed unit no support. Before
sampling, bisect its free log coordinates toward the common prior mode until all derived junctions
are positive. This only chooses a feasible starting point, preserves the target, and has no fitted
step-size or biological parameter. Keep the ordinary common-mode start as the comparator. Test
both starts and explicitly distinguish initialization failure from posterior mixing.

Next build/gates/repeated pilot: four minutes, under 4 GB per worker. No main-source promotion.

## Memory repair and analytical removal of missing positions

The four-stratum g50 test-chromosome screen improves transcript error (OFF/ON, stranded then
unstranded) from 9.236288/9.891189/10.547442/15.308427% to
7.261355/7.506631/9.347291/14.687100%. This is **not a pass**: unstranded OFF gene error rises
1.565813 -> 2.311910%, and several synthetic nascent allocations worsen. Numerical convergence
and attribution controls remain open. The run used 256 MCEM blocks of 32 sweeps, then 2,048 draws;
140.18 seconds, 4,704,256 KiB aggregate RSS with four workers (budget receipts, not speed claims).

The first RNA-short captured stranded pilot exceeded 180 seconds and reached 18,312,640 KiB RSS.
It has no accuracy result. Expand paths within each locus instead of keeping whole-library expanded
copies. Cache cumulative template coordinates once per template per deposit call; this reuses the
identical integer prefix sums. The 8,515-path identity remains within 3.55e-15. Preparation-only
gap diagnostic (zero fit/burn, one draw, **not an accuracy arm**) now finishes in 47.15 seconds,
8,786,224 KiB, including scoring 9,728,346 units and 103,093,845 candidates. It has 1,315,518
overhang candidates. Such candidates currently expand densely and sample positions that cancel.

For a missing-position candidate the earlier derivation gives a constant rate
`b_fk = phi_fk/L_k(1)`, independent of c. Integrate its latent object out analytically. Retain its
sampled component assignment in abundance counts, but not in the field's object counts or its
`-N_k log L_k` term. The latter now uses only assignments with observed spatial paths. The complete
field conditional is therefore `log prior + sum m_o log c_o - sum N_spatial,k log L_k`.
This is the marginal of the same augmented model, not a new error model or likelihood threshold.
For exact agreement at dependent-junction zeroes, a missing-position assignment still requires its
component length to be positive: preserve that support indicator using total sampled component
counts, even when its spatial count is zero. Keep the old dense expansion as a numerical comparator.
Tests must cover constant-rate responsibilities at several fields, mixed located/missing observations,
the integrated posterior against quadrature, and the zero-length support boundary before a BAM arm.

Budget the next native rebuild and numerical gates at three minutes / 4 GB; then a test-chromosome
identity/convergence pilot before another gap accuracy run. No model constants or capture detector.

Compact missing-position gates: fixed-field responsibilities match the dense representation exactly
at three fields, including the zero-junction face. Mixed spatial/missing-position quadrature gives
6.986212 RNA fragments; dense sampler 6.990196, compact 6.990457. Preserve zero numerical rates
as valid candidates, including empty vectors after floating-point underflow; never discard them with
an empirical threshold. An observation with zero total support still fails. The compact strOFF pilot
at 256 x 32 + 2,048 draws gives transcript 7.277881%, gene 1.254097%, 29.86 seconds / 740,736 KiB.

## First attribution controls (test chromosome, g50)

All pinned scan/fractional seed 0; percentage annotated-count absolute error. MAP-only retains the
shipped numerator and lengths. Fixed-reader uses the new candidate-specific numerator and coherent
length measure at the shipped reader's field. Flat uses that same new measure with c=1. Field32 is
the earlier 256 x 32 screen (before scalar marginalization); Field4 is the compact 256 x 4 screen
with 1,024 final draws. The latter is a **numerical screen, not a converged reference**.

| Stratum | Shipped VB | Shipped MAP | Flat MAP | Fixed-reader MAP | Field32 | Field4 |
|---|---:|---:|---:|---:|---:|---:|
| stranded OFF | 9.236288 | 11.205643 | 6.293275 | 6.866940 | 7.261355 | 7.419676 |
| stranded ON | 9.891189 | 10.305424 | 32.108723 | 6.578249 | 7.506631 | 7.631899 |
| unstranded OFF | 10.547442 | 11.793360 | 6.935230 | 7.386401 | 9.347291 | 9.634553 |
| deferred unstranded ON | 15.308427 | 17.421520 | 55.190869 | 11.984423 | 14.687100 | 14.951781 |

The fixed-reader control beats integrated field32 on transcript error on every row. This prevents
attributing the new model's gain over shipped VB to integration itself. Gene errors, same ordering:

| Stratum | Shipped VB | Shipped MAP | Flat MAP | Fixed-reader MAP | Field32 | Field4 |
|---|---:|---:|---:|---:|---:|---:|
| stranded OFF | 1.223986 | 1.248439 | 1.174021 | 1.248798 | 1.251237 | 1.252733 |
| stranded ON | 1.279939 | 1.335443 | 28.372790 | 1.311434 | 1.290556 | 1.292662 |
| unstranded OFF | 1.565813 | 1.604592 | 1.569625 | 1.606030 | 2.311910 | 2.353118 |
| deferred unstranded ON | 5.319915 | 5.463914 | 50.340285 | 5.499874 | 4.768604 | 4.782598 |

The short-block run finished in 35.88 seconds / 3,621,728 KiB for four workers; MAP/flat/fixed
controls took 9.61/7.53/7.50 seconds. These are budget receipts, not controlled speed comparisons.
Next RNA-short stranded ON pilot: 256 x 4, burn256, draws1,024, four-minute budget, expected12GB
with20GB hard worker guard. Hold all statistical assumptions fixed. This screen can expose large
defects but its 0.14–0.29 percentage-point block sensitivity forbids claiming fine convergence.

## Wide-support numerical falsifier, before repair

One missing-position observation has a constant likelihood, so with two independent log-uniform
fields on [-40,40], the posterior is exactly the prior. `P(|x1-x2|>40)=1/4` and
`P(x1<0,x2<0)=1/4`. The native incremental-length sampler at 65,536 draws gives
0.276748657 and 0.233261108. This falsifies the numerical kernel despite its moderate-range gates.
Subtracting `A*c_old` from a length dominated by that term can erase a small positive remainder;
the subsequent support check then rejects valid proposals. The common-scale orbit does not repair
that incorrect coordinate transition. Earlier BAM screens are retained but cannot establish full
sampler correctness. The active four-stratum RNA-short screen is stopped on this gate.

For each coordinate's affected component, compute once per coordinate the sum over **unchanged**
objects, using only nonnegative `A*c` terms. Evaluate each proposal as that base plus the nonnegative
sum over changed objects (including dependent junctions). Store the same canonical sum when accepting.
No subtraction, clipping floor, likelihood threshold, or model change is needed. The discrete reference
must likewise use canonical nonnegative sums, retaining the deliberate stale-length mutation only in
its named test branch. Recheck the exact prior-only distribution over a wide support and the prior
moderate-range, dependent-junction, missing-position and mutation gates before a new BAM arm.

## Exact factorization of equal candidate kernels, derived before prototype

The cached large locus has 189,126 observations, 8,298,415 candidates, 12,702,759 coefficient entries,
154 components and 537 objects. Its 66,363 observation groups spend 11.50 of 12.34 seconds in 192
E-steps, but only 0.294 seconds in field transitions. Reusing individual candidate numerators and
reciprocals increased this to 14.59 seconds and is reverted. The repeated matrix-vector work belongs
in an exact factorization, not another field approximation.

Within one observation f, partition candidates into classes Q whose **entire** coefficient vectors
`B_fk` are bit-identical, including object IDs and values. For a spatial class,

    v_k(c) = a_k/L_k(c)
    W_Q(c) = sum_(k in Q) v_k(c)
    H_fQ(c) = sum_o B_fQ,o c_o
    Z_f(c) = sum_Q W_Q(c) H_fQ(c)
    r_fQ = W_Q H_fQ/Z_f;  r_fk = r_fQ v_k/W_Q.

Membership sets Q recur across observations. Compute each W_Q once per sweep. Accumulate
`sum_f multiplicity_f r_fQ`, then distribute those expected counts over members with v_k/W_Q.
For a missing-position class with identical scalar b, use `H=b` and
`v_k=a_k 1{L_k>0}`; its assignments affect abundance and support, not spatial field counts.

Sample class/object counts for each observation group using its exact integer multinomial. Conditional
on its class, component identity is independent of the observation and its deposit object. Consequently
all integer assignments to the same membership class can be split by **one** multinomial with
probabilities v_k/W_Q. Add those integer component totals to the collapsed field statistics. This has
the same joint law for component/object sufficient statistics as the original per-observation draws;
no mean counts enter the field conditional. Keep spatial and missing-position classes distinct.

For fractional per-fragment delivery, still accumulate each candidate's `r_fQ v_k/W_Q` at each retained
draw. Factoring two separately averaged terms would be wrong. No guarantee that all kernels match is
assumed: unique kernels remain singleton classes. Equality is exact, without rounding, bins, tolerances,
depth cuts, or a changed likelihood. A shared-class audit must include unequal kernels at the same
objects, different objects with equal coefficients, a common class with unequal component lengths,
missing positions, and the dependent-junction support face. The finite posterior and fixed-field
reference must agree before any full-pipeline run. This is a representation of the same mechanism.

The factorization passes: enumerable joint-state maximum error 0.001997 (ordinary sampler 0.003751),
fixed-field fractional probabilities within 1.11e-16, and mixed missing/spatial counts within 0.0000663
fragments between samplers. The large cached locus profile is 4.34 seconds, with 2.39 in the E-step
and 1.08 in the now-stable field update. These cached profiles do not establish whole-pipeline timing.
Wide-support, finite-posterior/mutation, continuous-quadrature, MCEM, missing-position and dependent-
junction reference gates also pass with factorization enabled.

## Repaired first screen and an optimization hypothesis

The eight g00/g50 test-chromosome strata finished in 46.77 seconds / 5,368,816 KiB, four workers,
256 outer updates x4 sweeps, burn256, retained1,024. **Not a convergence or release pass.**
At g00, transcript errors (strOFF, strON, unstrOFF, deferredON) are
12.944025/17.364762/18.612468/18.179824%, versus shipped
7.661010/8.224251/5.624969/9.396448%. Gene errors are
0.855678/0.101941/5.614910/1.175434%, versus
0.834983/0.674347/1.013748/2.444940%. UnstrOFF synthetic RNA rises from43,327 to67,755 against43,034
truth; inferred EM gDNA remains only4.90 fragments. This is within-RNA confusion, not false gDNA.
At g50 the repaired field's transcript errors are7.285399/7.519984/9.393401/15.005010%, versus
9.236288/9.891189/10.547442/15.308427%. Gene regression on unstrOFF persists (1.565813->2.347278%).

Large g00 loci still have fixed-point count residuals of tens of fragments against approximately1–3
fragments of batch Monte Carlo error. A four-stratum1,000x16+2,048-draw convergence run is budgeted
at four minutes /6GB before attributing the regression to a converged model.

If ordinary MCEM remains the optimizer bottleneck, an exact alternative M-step can collapse the
component labels while retaining c as the EM latent variable. With
`q_t(c)=p(c|D,a_t)`, define

    Q_c(a|a_t) = E_qt [sum_f log(sum_k a_k K_fk(c)/L_k(c))] + log prior(a).

Then the change in the integrated log-posterior equals the change in Q_c plus
`KL(q_t(c)||p(c|D,a))`. Thus increasing Q_c increases the exact objective. For fixed posterior
field draws its Monte Carlo approximation is a concave mixture log-likelihood in simplex a, with the
same concave aggregate prior penalty. Its inner component-label EM is deterministic while those
draws are frozen. This permits ordinary SQUAREM **only inside that frozen-field M-step**, safeguarded
by the actual sampled log-sum objective; it does not permit extrapolating across moving field
posteriors. The sampled bound is not the exact bound, so seed/budget checks remain necessary.
This is derived but **not implemented or selected**. It must agree with finite integrated-EM optima
and reduce to ordinary fixed-field MAP before a pipeline comparison; do not silently replace the
objective with the shipped adaptive VB update, which is a different update rule.

The deeper g00 run finished in113.58 seconds /4,665,952 KiB. Transcript errors remain
10.954709/14.280981/14.905731/14.649606%; gene errors
0.849445/0.049779/4.925650/1.037791%. Largest fixed-point count residuals remain about6–11 fragments.
This improves on the short screen but neither converges convincingly nor repairs the zero-gDNA rows.

## Exposure-rate coupling: a separate model hypothesis, before prototype

The first arm holds **captured-fragment proportions** a fixed while integrating over c. A different
generative choice holds precursor positional rates lambda fixed. Then captured proportions are not
independent of c:

    W(lambda,c) = sum_k lambda_k L_k(c)
    a_k(lambda,c) = lambda_k L_k(c)/W(lambda,c)
    p(f|lambda,c) = sum_k lambda_k K_fk(c)/W(lambda,c).

The path-normalization identity gives sum_f p(f)=1 whenever W>0. A component with zero exposure has
zero captured probability; it does not need an undefined K_k/L_k. At a known field this is a bijective
reparameterization of the current mixture on positive exposures. Under an uncertain field the two
integrated models differ: this arm explicitly carries abundance/exposure covariance. This is a
**hypothesis**, not a declaration that it will resolve identification or the empirical length-law issue.

Retain the same empirical aggregate prior as a factor on captured origin probabilities,
`a_g(c)^P_g a_R(c)^P_R`, inside the field integral. Thus fixed-field optima must agree with the existing
grouped MAP control. The calibration quantities are unchanged and their evidence reuse remains disclosed.
There is no library reference or rate scale: multiplying every lambda by a common factor leaves all
probabilities unchanged. Use `sum lambda_k L_k(1)=1` only to choose numerical units.

For spatial observations, sample component/object assignments proportional to `lambda_k B_fko c_o`.
For a missing position, its already-derived kernel is `b_fk L_k(c)` with `b_fk=phi_fk/L_k(1)`;
sample its component proportional to `lambda_k b_fk L_k(c)` and integrate its object out. With
integer object counts m and missing-position component counts M, the collapsed field conditional is

    log prior(c) + sum_o m_o log c_o + sum_k M_k log L_k(c)
      + P_g log L_g(c) + P_R log W_R(lambda,c)
      - (N+P_g+P_R) log W(lambda,c) + constant.

Terms with zero coefficients are omitted, and all positive-log terms require positive arguments.
No per-fragment work enters that conditional. The background term is still spatially uniform before
capture under the common empirical length measure; no unobserved fragment-length law is invented.

An exact EM augmentation for the rate update uses **the analytic expectation**, not draws, of a
single auxiliary exposure `u|c,lambda ~ Gamma(T,W)`, `T=N+P_g+P_R`. This is not the earlier refused
per-component gamma sampler. Also analytically allocate the RNA-only prior term with
`h_k(c)=lambda_k L_k(c)/W_R(c)` for RNA components. The rate update is

    lambda'_k = [E N_k + P_g 1{k=g} + P_R E h_k(c)] / [T E(L_k(c)/W(c))].

All expectations use the same posterior at the current lambda, followed only by the unit normalization.
When no gDNA component exists, the RNA group prior is exactly1 and can be omitted algebraically.
The sufficient statistics still include only integer assignments of actual fragments. Fractional
empirical prior masses enter analytic factors, never fabricated integer observations.

Falsifiers before a BAM comparison: normalization, scale invariance, known-field grouped-MAP optima,
the analytic marginal gradient, monotone exact-enumeration EM, and posterior counts against an independent
optimizer. Then test the zero-gDNA chromosome first. Only if that survives should the separately named
rate arm receive the complete ladder/gap comparison; never mix it with the first arm's results.

The first native rate MCEM gate fails on the enumerable four-object fixture: objective -247.265184
versus -246.865089 from exact enumeration (gap0.400094). Fixed-rate fractional counts happen to be
field-invariant in this fixture, so their exact agreement hides a wrong exposure expectation. The
posterior puts probability0.6260 on field[3,1,10,1], but the coordinate chain stays near[1,1,3,1]; maximum
state-probability error is about0.64. This is correlated-field mixing, not permission to loosen the gate.
The first test harness also used bit equality after exp(log(grid)), corrected to numerical state
identification, and tried the log-uniform orbit on a nonuniform grid, correctly refused by the kernel.

An additional exact transition is derived before implementation: propose every free log-field coordinate
independently from its existing prior density, derive the dependent junctions, and accept using the
conditional **data-factor** likelihood ratio. The proposal prior cancels the target prior in Metropolis-
Hastings. Evaluate that data factor with integer m and missing counts, the existing origin prior factors,
and the same W denominator given above. For the continuous piecewise-linear log density, choose an
interval by its exact integrated exponential mass and invert that interval's exponential CDF; no proposal
width or empirical constant is introduced. Compose this full-support transition with the coordinate and
common-scale transitions once per sweep. High-dimensional acceptance may be negligible: report that
limitation, and do not promise it will solve general mixing. It specifically gives the small finite
instrument a route between modes which individual coordinates cannot traverse at useful cost.

With that additional transition, the finite posterior probability errors are0.000523 (no origin prior)
and0.001029 (fractional origin prior); deleting W gives0.312350. Four whole-model fits at128 outer blocks
x4096 sweeps have objective gaps0.0004713/0.0003563/0.0000219/0.0000054 from the exact reference.
The worst count difference is0.146 fragments out of300. Known-field grouped-MAP error2.75e-13;
missing-position/dependent-junction count6.924773 versus quadrature6.926928. Wide prior tail0.248779,
exact0.25, batchSE0.001687. Old captured-fraction reference/mutation gates still pass with the new arm off.
These are numerical gates, not evidence of a BAM accuracy improvement. A separately named four-stratum
g00 test run is budgeted240s/6GiB,1000x16,burn256,draw2048, four workers, pinned and fractional.

The exposure-rate g00 screen completes119.60s/4,640,304KiB. Transcript errors strOFF/strON/unstrOFF/
deferredON11.429524/15.450709/15.000810/15.530957%, versus shipped7.661010/8.224251/5.624969/9.396448%.
Gene errors0.847496/0.039809/4.715714/1.024916%. This arm does **not** repair the regression. Its largest
final log-rate score components are8.29/20.73/8.72/21.64, so it is not a fully converged model verdict.
A continuous prior-only analytic CDF check passes: exact0.6786302, samples0.6793060/0.6787415, SE0.0018242.

Return to the first captured-fraction model to isolate optimizer error before proposing another model.
The frozen-field M-step above is concave in a. Its deterministic EM map M has two iterates a1=M(a0),
a2=M(a1); let r=a1-a0 and v=a2-2a1+a0. A secant extrapolation is a(s)=a0+2s r+s^2 v. At s=1 it is
exactly a2. Start s=max(1,||r||/||v||); if the candidate is infeasible or its frozen-field log-sum
objective is below that of a2, halve its distance to1 until accepted (machine equality falls back to a2).
Optionally apply one ordinary M step and require objective nondecrease. This uses only the same frozen
posterior fields, so the acceptance objective is fixed throughout the inner solve; it does not extrapolate
across changing field posteriors. No clipping floor, artificial prior mass, or change to inference model.
The finite-reference and fixed-field gates must pass before any new BAM arm. First price on cached loci.

Frozen-field optimizer gates pass: fixed-field and finite posterior/mutation tests unchanged; four exact-
fixture abundance errors stay below0.001451. Cached strOFF g00 loci120/123/124 at64 outer blocks x32 retained
fields,128 equilibration sweeps per block,128 inner iterations, two seeds/dispersed starts, take110.46s
combined. Seed count-L1 differences122.39/1965.08/481.48 fragments. Final fixed-point count residuals
0.16/1.33,1.48/0.98,22.22/21.59 respectively. This does NOT pass a reliable-fitting gate. The old1000->10000
iteration check changed loci120/123/124 by1465.47/1142.20/2616.72 fragments; the manual probe was stopped
incomplete after six loci, so its declared120s was not a programmatic hard cap.

Before another BAM arm, distinguish local optimization from equivalent explanations. For two fitted
abundances a0,a1, draw fields separately under their posterior and evaluate delta(c)=logL(a1,c)-logL(a0,c).
The ratio of their integrated likelihoods is a bridge-sampling identity. With equal draw counts, d=logZ1/Z0
solves E_q0 logistic(delta-d)=E_q1 logistic(d-delta). Add the difference in grouped log-prior penalties to
d for the integrated MAP objective difference. Small overlap or disagreement between field-chain repeats
invalidates the estimate; do not select a winner from a tiny within-chain error. Verify the bridge against
finite enumeration before using the cached locus. This is an objective diagnostic, not a truth metric or
an added inference assumption.

Direct bridge sampling FAILS its overlap gate on locus123: overlap3.91e-37 and8.90e-38 in two repeats.
The apparent objective differences−18.80/−18.06 are therefore **not rankings**, despite within-chain
batch errors1.53/1.67. The finite bridge identity itself agrees with enumeration to5.55e-16.

A change of integration coordinates can improve this diagnostic without changing the model. In free
log-field coordinates use y=T(x)=x+d, with d fitted from a separate pilot half of each posterior chain.
The Jacobian is1. Compare unnormalised h0(x) to h1(T(x)); include BOTH field log priors and derive junctions
again from the transformed free coordinates. Outside the original bounded prior support, h=0 exactly.
For samples from q1, evaluate the inverse transform; samples mapping outside q0 have delta=+infinity,
not deletion or renormalisation. Samples from q0 mapping outside q1 have delta=−infinity. The same bridge
identity remains valid on the union of the supports; it still needs measured overlap and chain repeats.
No reference level, posterior clipping, prior refit, or likelihood approximation is introduced. This
transport diagnostic is the next isolated numerical measurement, not a new capture mechanism.

The additive transport remains inadequate: first overlap0.00452, and some diagnostic batches have no
support overlap at all (the original batch-SE code correctly cannot produce a finite estimate, but its
empty-array exception is a harness defect). Do not use its apparent +7.29 objective difference.

A support-preserving transport uses the existing continuous prior's CDF F. Set z=Phi^-1(F(x)); the prior
is standard normal in each z. Estimate means and Cholesky covariances from pilot samples only. The affine
map z1=mu1+L1 L0^-1(z0-mu0), then y=F^-1(Phi(z1)), is bijective on the SAME bounded field support.
Let H=L1 L0^-1. The log-Jacobian is sum log p(x)-log p(y) plus
log|det H|−(||z1||^2−||z0||^2)/2. Hence prior terms cancel analytically in the bridge log ratio; evaluate
the exact fragment log likelihood on each reconstructed field and add the normal-density correction and
log determinant. Invert the map for q1 samples. Refuse singular pilot covariance or nonfinite quantiles;
no jitter, density floor, clipped probability, new prior or detector. Compute log-CDF and log-survival
from exact piecewise-exponential interval integrals and invert by floating-point bisection, so small tails
need not be rounded to zero/one. Test inverse identity and Jacobian independently before comparison.

Prior-coordinate transport has inverse error2.84e-14 and Jacobian error4.69e-11. Its cached-locus bridge
still does not give a reliable ranking: overlaps0.0197/0.00195 and batch errors2.66/10.78 log units.
The apparent +7.15 difference is retained as a diagnostic, not selection evidence.

## Free positive junction object: separate geometry hypothesis, before prototype

The current comparator imposes c_J=clip(c_lo+c_hi−(c_il+c_ih)/2,0,max(...)). The earlier truth audit already
refuted this as a universal physical identity (truth44.5593 versus cap27.7797). It additionally permits
L_k=0 for a junction-only opportunity row while a captured mixture proportion remains positive: K/L is
then undefined. The first kernel used zero candidate rates on that face; this is not a general proof of
normalisation. Do not present the mixed missing-position face gate as such a proof. Audit actual rows
separately; no claim is made here that this face explains all BAM regressions.

Alternative: make each existing junction cut an ordinary positive shared field object, with the SAME
existing landscape prior density as the other objects. B and A are unchanged. Since A_k is nonnegative
and its row sum is positive, strict c_o>0 now implies L_k>0 for every component, with no floor. Every
fragment still enters once; all candidates sharing a cut use that same coordinate. Its integer field
sufficient statistic and collapsed conditional are the already-tested linear-object model. At a fixed
field this must be bit-identical to the dependent-junction arm supplied with the same field. No gDNA
opportunity is invented at the junction. Unobserved junctions carry negative evidence via the normalizer;
there is no claim of independent biological capture or unconditional abundance identification.

This changes one statistical relation, not an optimizer: use the original ordinary MCEM budgets on both
arms, not the unaccepted frozen-field optimizer or rate model. First truth screen g00 allfour strata,
1000x16,burn256,draw2048,240s/6GiB. Required falsifiers remain normalization, fixed-field identity,
seed/budget stability, and full ladder/both-gap truth A/B before any model recommendation.

The free-junction screen completes107.77s/4,790,832KiB. Transcript errors strOFF/strON/unstrOFF/
deferredON20.816344/40.555952/25.597426/43.002891%; gene errors7.938070/24.573163/11.834089/26.883946%.
Synthetic-RNA estimates122513/280656/142411/281075 versus truth43017/9476/43034/9293. Thus a completely
free independent junction field is not a usable repair at this budget. It increases within-RNA ambiguity
severely. Residuals9–15 fragments remain, so this is a failed implementation screen, not a converged proof
about every possible junction prior. Its known-field arm was exactly identical in probabilities, counts,
abundance, field and lengths. No new detector/reference/empirical constant was used.

Cached frozen-fit truth check (annotated transcripts only) further falsifies selection by small residual:
locus120 errors119.9998% for both seeds; locus123144.7078% versus89.2911%; locus12454.5502% versus56.3656%.
These are selected diagnostic loci, not genome-wide prevalence or aggregate accuracy. Calibration arrays
for rate and original dependent-junction g00 arms are identical, max delta0 across6arrays x4conditions.

## Current decision and unfinished release work

The native path-aware shared-field implementation exists outside main, with reference/mutation gates,
continuous posterior sampling, exact candidate-kernel factorization, and pinned fractional truth screens.
It is **not selected for production**. The ordinary captured-fraction model, exposure-rate alternative,
and independent free-junction variant all fail the g00 development screen. Numerical changes improved
correctness and cost but have not established a reliable optimizer or a release accuracy improvement.
No full ladder/both-gap field matrix or new field VCaP score has completed. Older calibration/readout
campaign A/Bs do not stand in for those missing arms. Existing primary calibration scores remain those
of the unchanged calibration stage. Do not mark 0.8.0 ready or the implementation objective complete.

The single next line of work is an **objective-checked joint-field fit**, retaining the path-specific
shared-field likelihood while resolving correlated posterior modes and abundance optimization before
another biological mechanism. The frozen-field M-step and posterior-coordinate transports above are
executable pieces, not a solution already proved. The required next falsifiers are: agreement with exact
integrated objectives; across-start posterior overlap and count stability on the bad cached loci; actual
normalization on every allowed junction-support state; test-chromosome truth improvement; then the whole
ladder, both gap directions, calibration companions, VCaP, and public output/annotation contracts. A raw
count-only Monte Carlo SE or small EM residual is explicitly insufficient. If adequate sampling leaves
truth worse, reject the model on that evidence; never pick an initial state by truth or tune a numerical
budget to improve a panel score.

Refused now: production/default promotion; claiming convergence from a fit cap; second-order broad-field
responsibilities; plug-in posterior averages already refuted in the prior campaign; single-object fragment
encoding; stale/subtractive length arithmetic; the hybrid variational objective; detector/reference/floor
repairs; independent free junctions as an already validated fix; reader deletion before a complete public
replacement; and any release/speed claim from these exploratory screens. These refusals concern present
implementation claims. They do not manufacture the missing full-panel scientific verdict.

The positive-free-opportunity support audit passes all four g00 test conditions in8.57s/3,921,568KiB,
with fitting disabled. This confirms that the newly identified junction-only zero-length face is not the
explanation for these g00 regressions. The adapter now refuses that unsupported geometry in future
dependent-field screens instead of silently treating an undefined component density as zero.
