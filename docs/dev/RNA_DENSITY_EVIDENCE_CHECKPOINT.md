# Density-evidence checkpoint and the next count-model decision

*2026-10-07. Development sandbox. Implements the approved diagnostic/reference checkpoint;
does not integrate a reader or change production count inference.*

## Completed work

The frozen count × opportunity comparison is complete on three test-chromosome capture-ON
conditions, RNA-short stranded OFF and RNA-long unstranded OFF. Every eligible object is
reported by class and admitted strand. Background, prior odds, ceiling and readout stay fixed
across substitutions. Fresh test arrays reproduce the earlier audit bit-for-bit.

The evidence primitive is implemented outside `src/`. It evaluates the two observed Poisson
strand columns while integrating unknown RNA with the existing total-RNA/tilt reference.
It covers no RNA, either RNA strand and both strands. Its public calculation takes no inferred
DNA count, posterior variance, landscape or distant RNA density. It is a small-case exact
reference, not the eventual whole-genome evaluation algorithm.

Validation:
- The old point-count reader fails the zero-DNA-support falsification before the replacement
  is written. The replacement passes it.
- All 113 observation-model tests pass, including independent direct integration for both-strand
  objects, zero observations/opportunity, rate tails, strand reversal and rate/opportunity units.
- Six deliberate mutations fail their targeted tests: the plug-in strand split, wrong RNA
  reference, omitted tilt atoms, flat zero-count evidence, a count-capped rate and wrong units.
- The oracle opportunity calculation passes an independent literal-placement test.
- Native local builders were replayed with only the incoming belief, factory contribution,
  coordinate origin or `split_live` bit changed. These are controlled diagnostics, not a new
  count-changing solver.

Settled measurements have moved to their issue homes under the move rule:
`ISSUES: calibration-detects-capture-on-a-capture-off-library` and
`ISSUES: the-gdna-prior-enters-psi-twice`.
The prior [reader review](RNA_SHORT_READER_REVIEW.md) retains the archived end-to-end tables.
No new transcript/gene accuracy result is claimed by these per-object diagnostics.

All code, arrays, hashes and logs are in
`.cache/rigel_runs/2026-10-07_density_checkpoint/`:

| Artifact | Purpose |
|---|---|
| `baseline.patch`, `baseline_head.txt`, `frozen_receipt.json` | Starting tree, source/native identity and cache provenance |
| `freeze_inputs.py`, five input NPZs | Fresh current calibration and oracle arrays |
| `factorial.py`, `factorial.json`, five factorial NPZs | Four-way count/opportunity substitution |
| `density_evidence.py`, `test_density_evidence.py` | Exact reference primitive and its independent tests |
| `legacy_falsification.log`, `mutations.json`, mutation logs | Failing old interface and deliberate-breakage evidence |
| `audit_messages.py`, `message_audit.json`, `message_local_inputs.npz` | Native local-message dependency audit |
| `false_boundary_evidence.json`, `test_oracle_geometry.py` | Actual ambiguous boundary and opportunity reference check |

## What the factory issue means

The factory's claim is effectively: “DNA here should look like the intergenic background.”
That is a prior. It can be useful, especially when the library has no strand information.
But transporting it to another object does not turn it into a new measurement.

If we put that claim inside the new evidence curve, train the landscape from the curve, and
then apply the landscape beside the factory, we have not achieved the intended separation.
If we simply delete the factory, we remove a constraint on which the existing unstranded
solver relies. The local audit makes both consequences concrete.

The own-observation primitive is therefore ready as a reference, but the entire neighbour
evidence interface is **not yet certified**. The population-partner experiment cannot honestly
start by treating today's delivered rows as likelihoods. Exporting a posterior and subtracting
one prior term is also insufficient because the strand widths already depend on that posterior.

## Authorized next prototype

The owner authorized removal of the intron constraint and continued investigation on
2026-10-07, after clarifying that the background-excess model is useful prior information,
not inherently invalid. This authorizes the prototype below; it does not establish that
removal improves accuracy or authorize reader integration. Start with a constraint-only
ablation to measure the information lost, then certify the local observation factors.

Test a **single landscape prior for count inference at every object**, with these exact terms:

1. Observation curves use the raw counts, declared RNA reference and independently certified
   local observation factors. The intron's own strand observations enter through the same
   model as other objects; its factory row is not an observation or a training kernel.
2. The count landscape is fitted from those curves. Its existing structural training exclusions
   remain initially; kernel replacement and admission/weighting changes remain separate contrasts.
3. Where a fitted landscape supplies the DNA prior, it is the only DNA prior at the target:
   no additional DNA Jeffreys tilt and no separate intron factory factor. The RNA reference stays.
4. The intergenic background remains measured and available to the eventual local capture
   consumer. It is not sent around the count graph disguised as extra observed DNA. The
   count landscape already has intergenic observations as training inputs.

This is a proposal to **test** the assembly, not permission to delete the production factory
or land prior-once alone. The no-fitted-landscape fallback is outside this first contrast and
must be specified separately before landing if it needs to change.

The smallest next work is the actual three-object face from the audit: derive its remaining
observational factor from counts/opportunities and compare to explicit latent enumeration.
Do not translate an old count-odds row into an intensity likelihood by renaming its axis.
Where the observations identify only a bound, retain a bound; do not invent a composition.
If the face requires a new spatial or capture-transfer assumption, stop and state it before
building the larger model.

After that reference passes, use an isolated native worktree and retain the existing count
path bit-identically in the diagnostic arm. Separate the observation-kernel, training-kernel
and prior-assembly contrasts. Only the successfully paired count candidate may progress
through the test chromosome and release panels. Never land the already-failed prior-once-only
arm as an intermediate release change.

The principal risk is explicit: ambiguous unstranded objects can lose composition information
that previously came from the factory. Neither the own-read reference nor a global capture
pattern is an acceptable substitute. The owner's subsequent 2026-10-07 clarification keeps
unstranded capture-ON deferred for 0.8.0 and authorizes investigating its regression for
robustness. Failure there is reported rather than treated as a release veto; it does not
authorize a new cutoff or another class prior.

## Boundaries that remain

- The capture prior and annotation-size dependence are still separate decisions.
- Honest counts do not fix component opportunities or DNA-to-RNA capture transfer.
- The sampled RNA-coordinate representations still need a proper numerical convergence gate;
  extreme tail log differences alone do not establish a practical error.
- `split_live` has a nonlocal input dependency, but this audit's triple did not show a delivered
  composition change. Do not report an unmeasured end-to-end effect.
- Exact own evidence is not guaranteed to identify capture when DNA is absent or strand
  information is missing. The implementation must express that uncertainty.

The factory decision is resolved for experimentation. Reader integration and any new local
biological assumption retain their owner checkpoints. Production source, existing prototypes
and goldens remain unchanged by this checkpoint.

The authorized continuation is implemented and reviewed in
[RNA_FACTORY_REMOVAL_CHECKPOINT.md](RNA_FACTORY_REMOVAL_CHECKPOINT.md). It includes the
isolated removal, native observational-message contrasts and one certified local face.
The remaining local interface and population partner still precede reader integration.
