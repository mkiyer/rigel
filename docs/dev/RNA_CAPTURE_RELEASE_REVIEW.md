# Rigel 0.8.0: capture reader and route to release

2026-10-09. External review brief. This file is a review snapshot, not a new
design ruling. Read [CLAUDE.md](../../CLAUDE.md) first. No source, production test,
golden, commit or push is part of this review checkpoint.

## The decision this review should help us make

The count foundation is integrated and validated. The detector-free capture
reader is not. Its mathematical reference is more defensible than treating an
inferred DNA count as an observed count, but its complete implementation is
expensive and has accumulated too much research machinery to copy into production.

Please identify the **smallest defensible completion**, including what to delete
or defer. Do not assume the proposed statistical model deserves to survive just
because we have written an exact evaluator. Conversely, do not replace uncertain
observations with confident estimates to make the evaluator cheap.

The owner wants robustness across fragment lengths, probe placements, capture
strengths, depths and real libraries, with roughly preserved accuracy. Improving
simulated accuracy is not the objective. Three release strata are required:
stranded capture-OFF, stranded capture-ON and unstranded capture-OFF. Unstranded
capture-ON is reported but deferred. Both-strand annotation objects also occur
in stranded libraries, so their numerical problem cannot simply be deferred.

## Plain-language description

Rigel answers two different questions about each region or boundary:

- How many fragments came from DNA versus RNA? The count solver combines the
  object's observations, local messages and a library-wide density prior.
- How much did capture enrich this object? That changes the opportunity used
  to quantify a transcript. Uncertainty in the first answer must not become
  confident evidence of enrichment in the second.

Different fragment lengths change the number of fragments that can fit inside
a region or cross its boundary. Counts therefore need to be divided by the
opportunity appropriate to their component. A count that is high merely because
more DNA fragments can fit must not be mistaken for capture. The original
RNA-short failure exposed this distinction.

The proposed reader asks which DNA densities are consistent with the actual two
strand counts and available local messages. It averages over possible RNA
amounts instead of subtracting a guessed RNA count. At overlapping-strand objects
it also averages over how RNA is divided between the strands. A broad range of
plausible DNA densities stays broad. This is what “honest evidence” means here;
it is not a claim that the inherited messages are exact independent observations.

The reader then compares that evidence with a background distribution and a
continuous range of enrichment. It does not decide whether the library is captured.
The background is shared; remote transcript expression and annotation size are
not reader inputs. Neighbour factors retain the existing propagation rules, so
conditional locality of this interface is weaker than complete pipeline locality.

## What is actually implemented

| Status | Implementation | Evidence / limitation |
|---|---|---|
| Integrated in the dirty main tree | Component-opportunity count repair; removal of the intron factory; observation-based strand messages; removal of unused plumbing. | Count-panel, real-library and landing receipts exist. Ambiguous weak-strand objects can over-call DNA; the owner accepts that conservative tradeoff. |
| Integrated separately | Local witness availability, positive numerical source coordinates, certified-source blur footprint and low-probability blur arithmetic repairs. | Isolated falsifications and deliberate native defects; the latest production suite has 2,814 passes, with lint and full preflight passing. |
| Still used in production | `capture_efficiency.py`, the landscape's enriched reference, clipping, and the library-wide `None` fallback. | This is the old reader. The detector-free release requirement is still unmet. |
| Prototype only | Density evidence with separate composition, DNA-level and RNA-level coordinates; independent profile support and local numerical units. | Frozen-input comparisons and independent references; source spacing remains an explicit failing gate. |
| Prototype only | Proper local capture prior and continuous correction score. | Mathematical limits, perturbations and small frozen-object screens. No full-panel or real-library result for this reader. |
| Not selected | Alternative landscape admission/weighting/population experiments. | They trade errors between classes. They are not dependencies of the proposed reader landing. |

Main remains at the owner's commit `a7b03103` with uncommitted work. Read the working
tree, not only HEAD. The count readiness audit describes its particular landing;
later repairs have their own receipts. Do not describe all archived “candidate”
results as results from today's complete tree.

The current source diff against that commit spans 15 files, with 306 added and
712 removed lines (net 406 fewer). That includes work before this checkpoint;
line count does not certify simplicity. The main complexity growth is in the
unselected research evaluators and harnesses. Their presence is not permission
to make them production dependencies.

The archived count-candidate report still contains individual high-DNA losses
against 0.7.1 and a severe deferred RNA-long regression. Its accepted foundation
status is not an all-condition pass or a detector-free release verdict.

## Exact candidate model and its assumptions

Let `rho` be DNA density, `Eg` its opportunity and `r` expected RNA count. With
one allowed RNA strand and protocol probability `q`, the two raw observations are

```text
u ~ Poisson(rho*Eg/2 + q*r)
v ~ Poisson(rho*Eg/2 + (1-q)*r).
```

The RNA nuisance reference is `r^(-1/2) dr`. The evidence integrates these
observations times the licensed composition factor at `log(rho*Eg/r)`, the DNA
factor at `rho`, and the RNA factor at `r/Er`. Their coordinates are distinct.
With both RNA strands, the existing reference gives equal mass to pure positive,
pure negative, and an arcsine continuum of strand shares; existing witness
exclusions remain. No fitted DNA count, DNA landscape prior or posterior variance
is inserted as a new observation. These modelling choices and the approximation
in the propagated factors deserve review independently of the quadrature.

For a normalized background `Q(rho)` with positive mean `mu`, the candidate has
equal reference odds for background and enrichment. Background has physical
capture label one. Enrichment has `rho=mu*C`, `C>=1`, and density `1/C²`.
This corresponds to a uniform retained background fraction `1/C`. It has no
upper capture ceiling. It scales the **mean** background; it does not multiply
an independently sampled background rate by capture.

```text
B = integral L(rho) Q(rho) d rho
S = integral from 1 to infinity L(mu*C) / C² dC
T = integral from 1 to infinity log(C) L(mu*C) / C² dC
log(weight) = max(0, (T-B)/(B+S)).
```

The last line is a **correction score**, not the posterior geometric mean of
physical capture. The background score is minus one and the enriched score is
`log C`, so flat evidence gives exactly weight one. The equal odds, proper tail
and loss/readout are choices, not conclusions forced by probability theory.
Ask whether this is a justified estimator of the quantity the length consumer
needs. The two controls using the old ceiling and old `1/n` odds are diagnostic
arms, not options proposed for production.

The reference currently takes a fixed Gamma background, including its point-mass
limit. It has not validated background estimation or propagated uncertainty in
the fitted parameters. Exact zero or unidentified background raises an error;
no positive floor or library capture switch is inserted. That explicit failure
is honest in a prototype and is **not a production resolution** for empty or
DNA-free inputs. A weight of one is the derived limit only under the stated
regularity conditions. See **Proper local capture reference and continuous
correction** in [EQUATIONS](../EQUATIONS.md).

## Numerical implementation under review

The scalar reference evaluates a complete evidence curve at each trial density.
At each density, a native integral averages over RNA amount; a Python outer
integral adds both RNA strands where needed. The capture calculation then
integrates that result against the background and enrichment priors. This
nested evaluation is the measured cost problem. The low-count difficult case
has only three reads: low depth does not make its message geometry cheap.

The authorized cheaper experiment changes only the outer evaluation. In
`x=log(rho/mu)`, the background, slab and log-moment integrands share the same
likelihood evaluation. The implementation accumulates the three components
together with `scipy.integrate.quad_vec` on each existing split interval. Input
factors, initial split points, inner tolerance, priors and readout stay fixed.
Its sole absolute tolerance is the smallest normal floating-point value, to
terminate fully underflowed intervals; it is not a biological rate floor.

The initial global-vector attempt missed about 16.6% of a narrow analytic
curve's slab mass while reporting convergence. The prewritten checks caught
this. Keeping the scalar reference's interval-local error checks resolves
that test. This is also why reported quadrature error is not proof of correctness.

**Outcome: reject this implementation for performance.** Four back-to-back
pairs agree within `7e-12` in log weight, but the shared calculation takes
2.4–2.6 times longer and increases expensive evidence calls. The old routine
already caches shared scalar nodes; the vector backend's subdivision work
outweighs any saving. No claim about all possible joint-integration methods
follows from this one negative result.

These are serial diagnostic Python/native-oracle timings, not whole-tool speed
claims. The two arms use identical input and native hashes, inner tolerance
`1e-8` and outer tolerance `1e-6`. The unstranded example uses a frozen stranded
background as a numerical control; it does not validate background fitting.
The detailed timing table has its permanent home under **SHARED CAPTURE
INTEGRATION** in [ISSUES](../ISSUES.md#the-gdna-prior-enters-psi-twice).
Both-strand timeout/cancellation receipts are reported in the artifact README.
Remaining jobs were stopped once the completed pairs rejected the implementation.
There is no new transcript/gene accuracy result and no production port.

This checkpoint passes 28 numerical tests, detects seven deliberately faulty
implementations and passes all 11 documentation-boundary tests. Hash checks
confirm 104 production source files and the installed native extension are
unchanged. Ruff passes on the prototype. The full 2,814-test production suite
was not rerun for this source-preserving experiment.

## Open issues that a fast evaluator would not resolve

| Concern | Required resolution or explicit release disposition |
|---|---|
| Unknown or tiny background | Derive a well-defined uncertainty-aware consumer for ordinary empty/zero-DNA inputs. Test the exact-zero limit; do not select a floor from simulations. |
| Background variation versus capture | Copy number, mappability and capture can produce similar local density changes. State the identification limit and test controlled regional shifts; do not promise immunity. |
| Sparse genuine capture and false correction | The current prior conservatively shrinks weak capture but permits small false weights. Evaluate expected sampling risk and downstream effect, not only selected successful objects. |
| Finite message support and narrow sources | Independent support/local units are prototypes. Three narrow-source interpolation checks still fail. Preserve the known source function through existing transport if necessary; do not globally refine all message tables without a measured need. |
| DNA-to-RNA transfer and length-law frame | DNA and RNA fragments overlap probes differently; the current RNA length law is capture-selected while DNA's law estimates the pre-capture frame. Component-specific divisors alone do not remove that mismatch or identify RNA capture opportunities. Price matched versus transferred expected yields across length gaps/layouts before claiming robustness. |
| Junctions and consumer scale | The old price is capped at one in a fully captured reference scale. A background-relative weight uses another scale. Derive the conversion and local junction rule; neither “max of neighbours” nor a new cap is automatically authorized. |
| Inherited local factors | They are approximations and can encode false confidence. Distinguish quadrature error from wrong evidence; do not add a second prior to hide the latter. |
| Real-data coverage | Existing real libraries are all stranded. This is a validation gap, not evidence that unstranded real data pass. |

The two captured weak exons in the proper-prior screen read about 444× and 478×
against expected capture 787×; the zero-DNA boundary reads 1.14× with messages.
Those are object diagnostics with frozen backgrounds, not release-accuracy
claims. The scalar screen's exact pure-DNA risk calculation also shows why one
DNA read cannot reliably recover a large capture factor.

## Complexity budget and proposed release path

Keep the accepted count solver and landscape fixed. Do not reopen population
training to chase the deferred stratum. Do not port the stacked prototype
integrators, import-hook harnesses, obsolete priors or diagnostic branch options.
Archive their receipts and use one ordinary worktree build for each proposed
production comparison.

The desired final shape is one typed local-evidence interface, one capture
consumer, and the existing length consumers. Remove the old reference reader,
clip and `None` machinery only after the replacement passes and the owner
approves integration. Do not keep both mechanisms behind a new feature flag.
Retain portable mathematical and mutation gates; archive research runners.

My recommendation is to keep the count foundation and the observation-based
evidence principle, reject this vector implementation, and obtain this review
before more numerical development. In particular, review whether the consumer
can compute its three required averages directly, without constructing an
entire marginalized density curve. Interchanging the nonnegative integrals is
mathematically valid; it does not establish a fast algorithm. Require a compact
implementation sketch and an independent small-case test before prototyping it.
If that requires a second general integration framework, reconsider the
consumer's approximation contract with the owner rather than silently accepting
that complexity. A faster implementation of the present formula would still
need the background, transfer and junction decisions listed above.

There are three finite delivery steps, in this order:

1. **Select a maintainable consumer.** Review the statistical target and the
   complete numerical contract, including both-strand evidence and zero
   background. Require bounded cost and equivalence to independent small-case
   references. If the cheaper experiment still fails, stop further quadrature
   work for review rather than commission an increasingly general integrator.
2. **Demonstrate robustness with the model frozen.** First expected-yield and
   observation-level checks varying length gaps, probe layout/count, capture
   strength, depth and background variation. Then test chromosome, both gap
   panels and ladder, every stratum separate and zero controls visible. Give
   every failing case a mechanism diagnosis; no retuning to individual panels.
3. **Land and consolidate after the owner checkpoint.** Portable failed-first
   tests and deliberate faults; reviewed goldens only; DESIGN/EQUATIONS/ISSUES
   updated together; full re-derived suite, Ruff and full preflight. Run real
   libraries serially last, assess runtime/memory, document residual limitations
   and check the publishing checklist. Commit/push/release stay owner actions.

These steps do not erase separate release-critical items already in
[ROADMAP](../ROADMAP.md). The reader is the present local bottleneck, not a claim
that every other release checklist item has closed.

## Minimal reading and reproduction

Read this brief, then the current sections of
[the active plan](RNA_SHORT_FIX_PLAN.md),
[the prior checkpoint](RNA_LOCAL_CAPTURE_PRIOR_CHECKPOINT.md), and
[the count landing audit](RNA_COUNT_READINESS_CHECKPOINT.md#final-engineering-landing).
Consult [all per-stratum count results](RNA_FULL_PANEL_VALIDATION.md) and
[real-library results](RNA_REAL_LIBRARY_VALIDATION.md) for accuracy context.
Their old-reader results cannot be advertised as new-reader performance.

Inspect these implementation paths:

- `src/rigel/calibration/capture_efficiency.py`, `calibrate.py`, `landscape.py`
  and `src/rigel/native/transfer_rows.h`: current reader, assembly, population
  and observation-based messages.
- `.cache/rigel_runs/2026-10-09_capture_prior/`: `local_prior.py`, `evidence.py`,
  `complete.py`, `test_local_prior.py`, `mutations.json`, `screen.json` and
  `verification.json`; six frozen typed input files.
- `.cache/rigel_runs/2026-10-09_shared_capture/`: bounded same-model experiment,
  its tests, paired receipts and reproduction instructions.
- `.cache/rigel_runs/2026-10-09_rna_spacing/`: narrow-source falsifications.
- `ISSUES: the-gdna-prior-enters-psi-twice`, especially **PROPER LOCAL CAPTURE
  PRIOR**, **NARROW-SOURCE INTERPOLATION** and the integration findings.

Prototype archives are local and ignored by git. A reviewer with only a clone
will need access to these artifacts; state any missing input rather than
assuming the brief proves the implementation. Do not require them to rebuild
whole libraries to give a design critique.

## Copyable prompt for the reviewer

```text
Review Rigel's path to a robust, maintainable 0.8.0 release. Start with CLAUDE.md
and docs/dev/RNA_CAPTURE_RELEASE_REVIEW.md, then inspect the named current
source and prototype files/receipts. This is a critique of the design and
implementation decisions, not a request to optimize synthetic accuracy.

Work read-only except for a new docs/dev/RNA_CAPTURE_EXTERNAL_REVIEW.md.
Do not change source/tests/goldens, commit, push or run full panels/real
libraries. Cheap independent calculations are welcome: first state their
expected time and check uptime. Identify unavailable artifacts explicitly.

Assess whether the proposed local-evidence model and correction score are
the simplest defensible way to satisfy the detector-free spectrum ruling.
Separate statistical assumptions, evidence/identification defects, numerical
error and avoidable engineering cost. Both-strand objects occur in stranded
data: preserving them matters even though unstranded capture-ON is deferred.
Review the new shared-integration result without assuming a native port will
solve a poor algorithm.

Answer these questions:
1. What should ship unchanged, be deleted, be deferred, or block release?
2. Is the quantity estimated by the proposed score the right input for the
   length consumer? If not, give the smallest replacement and its derivation.
3. Can the complete local evidence be consumed much more simply while
   retaining uncertainty, both RNA strands and local neighbour information?
4. How should background uncertainty/zero DNA, component opportunity
   transfer and junction scaling be resolved without a detector, hard cutoff,
   probe BED or simulation-specific constants?
5. Which remaining input/numerical defects materially block release, and
   what is the cheapest decisive test for each?

Every criticism must include a specific fix or explicit deferral, its
theoretical justification, a falsification test, and the cheapest experiment
that would reject it. Cite code and receipts, distinguish verified findings
from hypotheses, and do not treat prototype scores as production results.
End with one recommended implementation path, a short ordered release gate
list, stopping criteria, and only the owner decisions genuinely needed.
Favor deletion and a justified approximation over a second general solver.
Do not silently relax the no-detector, local continuous-weight or
fragment-length robustness requirements.
```
