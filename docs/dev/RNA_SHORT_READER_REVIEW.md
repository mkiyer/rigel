# Review of the count-inference capture reader

*2026-10-07. Sandbox review, not an approved model change or a release sign-off.
Reviews the count-inference reader recorded in §10 of plan v5, now preserved at
`.cache/rigel_runs/2026-10-06_plan_snapshots/RNA_SHORT_FIX_PLAN_v5_before_density_evidence.md`.
The current [RNA_SHORT_FIX_PLAN.md](RNA_SHORT_FIX_PLAN.md) incorporates these findings.*

**Recommendation: keep the reader outside production. Continue the counts track, but do not assume
that repairing point estimates alone will finish the reader.** The count defects reproduce. There
are also an overconfident evidence interface, a nonlocal prior ceiling, opportunity error, and
numerical edge failures. Each needs a separate intervention; none calls for tuning against a panel.

In plain language: calibration sometimes says “perhaps a few of these reads are DNA.” The reader
then treats those few reads as confirmed DNA observations. Against an almost empty DNA background,
that can become an enormous enrichment estimate and distort the RNA lengths. Better counts help,
but keeping the uncertainty attached to the evidence matters too.

## What this review actually ran

- Read the selected Step 3 implementation: `reader.py`, `capture_weight.py`, `hook.py`,
  `oracle_w.py`, the class analysis, arm runners, and the relevant calibration code and issue history.
- Reaggregated the archived 0.7.1, count-repair/current (“terminus”), and library-ceiling reader
  receipts. Condition/axis keys and truth totals must match before comparison. No pooled panel score.
- Ran three fresh cached test-chromosome calibrations: g00 ss0.99 ON, g25 ss0.70 ON, g50 ss0.99 ON.
  Each took 0.28–0.31 s excluding imports. Every object's observed total matched the oracle exactly.
- Ran small input perturbations, four falsification tests, and an independent numerical
  calculation on 127 frozen objects. The four falsification tests fail on the unmodified prototype.
- Replaced counts, then opportunities, on frozen inputs only. No scanning, EM, new synthetic
  libraries, human-genome debug dump, or full release benchmark.

Receipts and runnable audit code are in
`.cache/rigel_runs/2026-10-07_reader_review/`: `count_audit.py`, `counts.json`,
`counts_*.npz`, `receipts.py`, `strata.json`, `edge_checks.json`,
`test_reader_regressions.py`, `regressions.log`, `ci_integrals.py`, `integrals.json`,
and `counts_vs_opportunity.py/json`. The native module used by the fresh calibrations hashes to
`e7e0aba19c44947ac9977757c8c408ed262aadda9b8141b37a79f588e3b9a802`;
the Python prototype hashes are in `counts.json`.

No existing implementation, production source, production tests, or goldens were changed by this
review. The numerical alternative is an independent scratch reference, not a replacement reader.

## What holds, and what needs correcting

**The false capture reproduces exactly.** On test g00 ss0.99 ON, regions-then-boundaries object 2016:

| Quantity | Fresh value |
|---|---:|
| Observed unspliced reads | 4,930 |
| Inferred gDNA reads | 18.5657 |
| True gDNA reads | 0 |
| Median gDNA fraction | 0.00376586 |
| Variance of log gDNA fraction | 20.0737 |
| Reader capture weight | 23,831.03 |

The background is positive and uncertain (`mu = 3.65012e-6`, `alpha = 0.140174`), as intended.
This is not a division by an exact zero. Nor is a positive posterior median, by itself, proof that a
solver has broken: uncertainty near a nonnegative boundary can give a positive median. The harmful
step is presenting that median-derived estimate as independently observed DNA.

**The exon count bias is real, but not a universal 5% error.** The following is
`100 × (sum estimated / sum true − 1)`, with classes and admitted strands separate. It includes
every positive-opportunity object in the class, with no selection on the inferred count. Region
and boundary incidences are never added into a conserved-fragment pool.

| Object class | g25 ss0.70 ON | g50 ss0.99 ON |
|---|---:|---:|
| Positive-strand exons | +5.52% | +0.42% |
| Negative-strand exons | +3.60% | +0.13% |
| Both-strand exons | +10.74% | −10.89% |
| Positive-strand boundaries | −1.97% | −1.73% |
| Negative-strand boundaries | −4.85% | −1.22% |
| Both-strand boundaries | +25.52% | +4.04% |
| Boundaries admitting no RNA | 0.00% | 0.00% |

**Specific fix:** trace the largest signed errors separately through their own strand row,
received face rows, opportunity ratios, and density prior. Perturb one contributor with all other
inputs frozen. Do not introduce an exon multiplier, an ss0.70 branch, or a variance threshold.

**The remaining bias is not entirely count inference.** On g25 ss0.70 ON, 206 no-RNA boundaries
with expected capture factor above 20 have exactly correct gDNA counts. Their median capture-factor
error is nevertheless −0.05038 in log units, about 4.91% low. Substituting true counts leaves it
unchanged. Substituting the simulator's gDNA opportunities next moves it to −0.02212, about 2.19%
low. Their fitted opportunity is 1.02866 times the simulator opportunity. Background, prior,
ceiling and readout were held fixed. The remaining 2.19% has not been attributed.

**Specific fix:** retain a three-step frozen-input comparison: inferred counts → true counts →
true component opportunities. Score expected yields, not realized counts, for the resulting
weights. Use realized origin counts for the count census. Investigate the opportunity estimator
where the third step matters; do not replace a learned law with the simulator's law in production.

**Oracle final weights are not an isolation of count inference.** They bypass counts, their
Poisson interpretation, opportunity-to-density conversion, the background-relative readout, and
its prior simultaneously. They support the usefulness of the downstream length composition.
They do not establish that the reader is already as accurate as honest count inference permits.

| Test condition, capture ON | Current tree | CI reader | True final weights |
|---|---:|---:|---:|
| g05 ss0.99 | 5.23 / 0.35 | 7.82 / 0.50 | 5.36 / 0.31 |
| g25 ss0.70 | 6.49 / 1.34 | 8.93 / 2.37 | 7.00 / 1.02 |
| g25 ss0.99 | 7.38 / 0.80 | 8.06 / 1.15 | 6.54 / 0.86 |

These are transcripts / genes (%). “Reaches or beats the current tree” must be qualified by
metric and condition. The class-specific oracle runner also omits
`RIGEL_READER_CEILING=library`; unchanged classes therefore use the older own-content ceiling
unless the caller exported that variable. Its 9.01% baseline is not the latest 8.93% baseline.

**Specific fix:** put the ceiling explicitly in every arm's environment and receipt, regenerate
only the class contrasts needed for a decision, and compare to the same reader baseline.
The frozen comparisons above require no EM.

## The next statistical change should preserve evidence, not filter uncertainty

The selected path is `reader.object_weights` → `capture_weight(..., strands=NONE)`. Its likelihood is

```text
L(rho) ∝ (rho Eg)^k exp(−rho Eg),
k = calibration's median-derived gDNA count.
```

That is an appropriate likelihood for observed Poisson counts. It is not established for a
posterior point estimate produced using neighbours and the landscape. For the false boundary,
the implied Poisson log-rate variance is approximately `1/18.57 = 0.054`, versus the solver's
reported log-fraction variance of 20.07. These variances describe different distributions;
their mismatch exposes the discarded uncertainty, not a prescription to add them.

The two failed arms show that Gaussian blurring by the reported variance and discarding
“unlocated” objects are poor fixes. They do **not** show that propagating evidence is unnecessary.

**Specific fix to prototype, with an owner checkpoint before changing the model:**

1. Define one per-object density-evidence interface, including the neighbour information needed
   by unstranded and both-strand objects. It must represent a function of candidate density,
   including its low-density support and Poisson upper tail, rather than a centre and one variance.
2. Derive it from the observation model and the existing factor semantics. Test small examples
   against explicit latent-count enumeration. A posterior divided by the landscape is not enough:
   strand widths depend on the incoming belief, and the intron factory is itself a background prior.
3. With received rows and observations fixed, changing only the tested object's density prior must
   not change its evidence profile. Re-running neighbours under a different prior is a different
   intervention; do not claim global prior independence from this frozen-input test.
4. Use the same evidence interface to investigate the landscape's point-count kernels and the
   partner for applying the density prior once. Do not add a second population model to the reader.
   Keep the count and capture readouts' different prior assumptions explicit.
5. Test flat evidence, zero counts, vanishing opportunity, known pure DNA, weak real capture,
   both strands, and unstranded neighbour evidence. Do not accept a zero-control gain that comes
   from suppressing all weakly supported capture.

This is a bounded diagnostic project first. It does not authorize replacing the whole inference
engine, applying prior-once alone, reviving a class split, or returning to own-strand evidence as
the sole source. The existing prior-once experiment already demonstrates why a partner is needed.

## The library ceiling repairs truncation but introduces a new dependency

The implemented ceiling is the largest `(total reads + 0.5) / Eg` among objects with `Eg >= 1`.
It uses RNA reads as well as DNA reads.

A counterexample holds the local inferred count (10), opportunity (10), background
(`mu=0.1, alpha=20`), object count (`n=3500`), and all inferred gDNA counts fixed.
Increasing only another object's RNA-dominated total from 100 to 100,000 changes the local weight
from **8.15385 to 6.71449**. This is a prior-normalization effect, not a common unit conversion.
It violates the stated frozen-input locality property even after accounting for the known
`1/n` dependence.

The `Eg >= 1` rule also creates an opportunity cutoff and crashes when every positive opportunity
is fractional. This is in the ceiling calculation, not the previously refused landscape-training
guard, but it still lacks a model derivation.

**Specific fix:** separate numerical integration limits from the statistical prior. Remove the
scan of other objects' RNA densities and its whole-placement cutoff from the intended production
interface. A normalized, proper enrichment prior with unbounded support is one way to do this.
Simply increasing a normalized log-uniform ceiling until the answer “converges” is not a fix:
it changes the prior's odds.

A concrete candidate for owner discussion is a uniform prior on the background fraction
`b = 1/C`, giving slab density `p(C)=1/C²` for `C >= 1`. It has no upper cutoff, normalizes
exactly, and has `E[log C]=1`; the same counterarm construction remains defined. This is a **new
prior assumption**, not a consequence of the current evidence, and was **not implemented or
selected by this review**. It could penalize scarce true capture more than the current slab.
Any approved replacement must pass the single-probe/low-count controls as well as the locality
test. Do not pick a tail or its scale by optimizing the ladder.

The existing `p=1/n` choice is separately sensitive to annotation size. Its declared multiplicity
interpretation must be explicit; test scarce capture at genome-scale `n`, not only test-chromosome
`n`. Better counts do not settle either prior choice.

## A smaller and faster numerical implementation is available now

The current grid step depends on the ceiling-implied count, not the likelihood's observed count.
For the fresh g00 condition, the formula requests **357 million grid cells across 3,518 objects**,
including **3.87 million for one long region**. This was calculated from array sizes; the full
grid sweep was not run.

For the selected Poisson input, the entire background integration can be analytic. Let
`a0 = mu Eg`, `W = log(rho_max/mu)`, and drop the common count-factorial term. Define

```text
B = integral L(rho) background(rho) d rho
S = (1/W) integral from 0 to W of L(mu exp(t)) dt
T = (1/W) integral from 0 to W of t L(mu exp(t)) dt

log w = max(0, (T − B W/2) / ((n−1) B + S)).
```

This is the current counterarm formula after algebraic cancellation, in the continuous-integral
model. It avoids storing a posterior grid and avoids the `P(background)/(1-p)` division that
crashes at `n=1`. Its continuous `n=1` limit is defined; flat evidence gives exactly zero log
weight.

For a finite-shape Gamma background:

```text
log B =
    k log(a0/alpha)
  + log Gamma(alpha+k) − log Gamma(alpha)
  − (alpha+k) log(1+a0/alpha).
```

The independent scratch reference evaluates this in log space and integrates only the two
one-dimensional slab quantities. It preserves the current statistical assumptions and ceiling.
It also avoids taking a logarithm of a Gamma quantile that has underflowed to zero.

| Frozen condition | Objects | Original grid time | Integral reference time | Largest absolute log-weight difference |
|---|---:|---:|---:|---:|
| g00 ss0.99 ON | 42 | 0.435 s | 0.0160 s | 7.2e−15 |
| g25 ss0.70 ON | 43 | 0.0359 s | 0.0174 s | 5.8e−7 |
| g50 ss0.99 ON | 42 | 0.0317 s | 0.0169 s | 3.5e−13 |

These are back-to-back object evaluations on cached synthetic inputs, not a production or real-data
speed claim. The continuous integral is not bit-identical to the original grid approximation.
Before adoption, expand numerical equivalence tests across the supported input domain and record
back-to-back timing on a real library without whole-genome debug. A non-Poisson evidence profile
needs its own integration path; do not hard-code the point-count assumption into the new interface.

**Reproduced numerical failures and specific fixes:**

| Failure | Fix |
|---|---|
| `n=1` divides by zero, even with flat evidence | Use the cancelled component-integral formula and test its continuous limit. |
| All positive opportunities below 1 cause an empty maximum | Remove the ceiling's whole-placement selection with the approved prior redesign; do not add a cutoff fallback. |
| A valid highly dispersed fitted background gives `log(Gamma.ppf(...))=-inf` and crashes grid allocation | Use analytic background integration for the Poisson case; a general profile integrator must handle the zero-density tail without constructing an infinite log grid. |
| Remote RNA changes a fixed local weight | Remove the remote-density-dependent prior support; test with background, inferred evidence and `n` fixed. |

The dispersed-background failure is reachable from the existing fitter: 101 equal-opportunity
regions, 100 reads in one and zero in the rest. It is not an invalid manually invented Gamma shape.

**Code simplification at integration:** keep one selected evidence path, a background fit, a
readout, and the object-to-length composition. Archive the failed blur/located and own/CI-switch
arms outside production. Drop duplicate `library`/`noceiling` names. Remove strand and protocol
parameters from the pure CI numerical primitive. Pass explicit results instead of the hook's
global `_state`, fake reference/member fields, and `w/max(w)` adapter. Retain a common-rescale
test on the actual transcript and gDNA consumers. Update the docstrings that still describe the
own-read-only implementation.

## Release evidence, without pooling strata

“Current tree” below is the archived count-repair `terminus` arm used in §10, not the original
RNA-short-regressing tree. These are archived end-to-end results, reaggregated in this review.

Transcripts / genes (%); lower is better. Archived runs, reaggregated.

| Panel | Strand fraction | Capture | gDNA | 0.7.1 | Current tree | CI reader |
|---|---:|---|---|---:|---:|---:|
| test | 0.50 | OFF | nonzero gDNA | 9.25 / 1.76 | 8.66 / 1.33 | 9.26 / 1.39 |
| test | 0.50 | OFF | zero gDNA | 9.38 / 1.49 | 5.71 / 1.01 | 5.71 / 1.01 |
| test | 0.50 | ON | nonzero gDNA | 18.94 / 10.31 | 15.97 / 5.21 | 16.62 / 5.50 |
| test | 0.50 | ON | zero gDNA | 17.93 / 7.39 | 9.23 / 2.44 | 9.23 / 2.44 |
| test | 0.70 | OFF | nonzero gDNA | 9.44 / 1.54 | 8.62 / 1.25 | 8.40 / 1.25 |
| test | 0.70 | OFF | zero gDNA | 10.06 / 1.47 | 6.59 / 0.90 | 6.59 / 0.90 |
| test | 0.70 | ON | nonzero gDNA | 10.28 / 3.32 | 7.69 / 1.69 | 8.76 / 2.37 |
| test | 0.70 | ON | zero gDNA | 13.51 / 6.20 | 8.27 / 2.16 | 8.27 / 2.16 |
| test | 0.99 | OFF | nonzero gDNA | 8.20 / 1.20 | 7.67 / 1.08 | 7.61 / 1.08 |
| test | 0.99 | OFF | zero gDNA | 11.35 / 1.26 | 7.71 / 0.83 | 7.71 / 0.83 |
| test | 0.99 | ON | nonzero gDNA | 8.58 / 2.45 | 7.58 / 0.96 | 9.37 / 1.25 |
| test | 0.99 | ON | zero gDNA | 14.22 / 5.83 | 7.35 / 0.67 | 12.01 / 0.67 |
| ladder | 0.50 | OFF | nonzero gDNA | 10.29 / 1.49 | 2.29 / 0.31 | 2.32 / 0.31 |
| ladder | 0.50 | OFF | zero gDNA | 20.09 / 2.53 | 1.60 / 0.14 | 1.61 / 0.14 |
| ladder | 0.50 | ON | nonzero gDNA | 33.21 / 14.25 | 12.40 / 4.47 | 12.84 / 4.61 |
| ladder | 0.50 | ON | zero gDNA | 10.22 / 3.18 | 7.27 / 3.22 | 7.27 / 3.22 |
| ladder | 0.99 | OFF | nonzero gDNA | 3.95 / 0.44 | 2.02 / 0.25 | 2.04 / 0.25 |
| ladder | 0.99 | OFF | zero gDNA | 3.70 / 0.47 | 1.93 / 0.14 | 3.31 / 0.21 |
| ladder | 0.99 | ON | nonzero gDNA | 9.91 / 2.52 | 3.88 / 1.38 | 5.83 / 1.31 |
| ladder | 0.99 | ON | zero gDNA | 20.22 / 3.47 | 6.47 / 2.36 | 11.25 / 3.30 |
| flgap_rna_long | 0.50 | OFF | nonzero gDNA | 5.29 / 0.59 | 2.16 / 0.22 | 2.47 / 0.22 |
| flgap_rna_long | 0.50 | ON | nonzero gDNA | 28.99 / 5.84 | 6.55 / 0.94 | 7.07 / 0.96 |
| flgap_rna_long | 0.99 | OFF | nonzero gDNA | 5.23 / 0.46 | 2.27 / 0.19 | 2.27 / 0.19 |
| flgap_rna_long | 0.99 | ON | nonzero gDNA | 7.38 / 1.62 | 4.23 / 0.66 | 5.88 / 0.68 |
| flgap_rna_short | 0.50 | OFF | nonzero gDNA | 6.00 / 0.61 | 3.79 / 0.27 | 3.70 / 0.27 |
| flgap_rna_short | 0.50 | ON | nonzero gDNA | 29.76 / 6.19 | 17.57 / 2.23 | 21.90 / 2.57 |
| flgap_rna_short | 0.99 | OFF | nonzero gDNA | 6.31 / 0.44 | 3.53 / 0.23 | 3.55 / 0.23 |
| flgap_rna_short | 0.99 | ON | nonzero gDNA | 20.33 / 3.52 | 13.33 / 1.40 | 17.20 / 1.30 |

The major §10 numbers reproduce. Two qualifications matter for the strict release bar:

- Test unstranded OFF is 9.26 versus 9.25 in transcript error: approximately tied, but not a
  literal pass without an agreed repeatability tolerance.
- Ladder unstranded ON at zero gDNA has genes 3.22 versus 3.18. This is inherited from the current
  tree, but it still prevents saying the candidate passes 0.7.1 on every reported metric.
- Test ss0.99 ON remains a material transcript failure, 9.37 versus 8.58.

Do not change acceptance tolerances to make this candidate pass. Establish run-to-run variation
independently where a difference is tiny.

## Proposed next checkpoint

1. Keep the current count-repair tree and reader prototype frozen as separate baselines.
   Consolidate the plan's active specification: §§2–6 still describe own-read evidence and the
   superseded ceiling, while §10 runs a different interface. The standing suite count is 2,733.
2. Finish the ladder g00 attribution using one cached census. Rank false objects and trace their
   own evidence, received rows and length influence. Confirm the important object intervention
   on the test chromosome before spending another ladder EM run.
3. Prototype the density-evidence interface on frozen examples from the confirmed failure
   classes. Derive and test prior separation; separately price own-strand likelihood widths,
   received-face uncertainty, and the opportunity mismatch. No class- or condition-specific
   correction factor.
4. Resolve the capture prior's support/locality decision separately. Keep the analytic integration
   work as a numeric simplification with equivalence gates, not an accuracy arm.
5. Only then screen the test chromosome, both gap directions, all strata and zero controls.
   Use expected-yield checks across depth, continuous strand fidelity, scarce single probes,
   probe layouts, and regional DNA-rate changes. A local DNA-density increase can be copy number
   or mappability rather than capture; no estimator here has proved exact immunity.
   Fix choices before using these checks, and reserve independent noise realizations and real
   libraries for validation.
6. Integrate only after the owner checkpoint, with the move rule, falsification and perturbation
   gates, ruff, the full re-derived suite, and full preflight. Real libraries run one whole-genome
   process at a time. No new native implementation outside an isolated worktree.

**Decision to discuss:** make the next checkpoint honest density evidence and the demonstrated
opportunity/count separation, while keeping the new reader out of `src/`?
