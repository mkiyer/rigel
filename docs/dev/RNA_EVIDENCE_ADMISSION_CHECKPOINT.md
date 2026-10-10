# Evidence-curve admission and equal-weight comparisons

*2026-10-08. Prototype outside production. No commit or push.*

Removing the posterior-variance admission cutoff is safe on the measured zero-DNA
control, but it trades intron improvement for worse exon and boundary counts on the
deferred RNA-long captured condition. This population model is not ready to land.

The settled census, scores and smoothing attribution belong to
`ISSUES: the-gdna-prior-enters-psi-twice`, under **EVIDENCE-CURVE ADMISSION**. The preceding
interface implementation is [RNA_CHANNEL_EVIDENCE_CHECKPOINT.md](RNA_CHANNEL_EVIDENCE_CHECKPOINT.md).

## What changed

The only new statistical intervention admits expressed, single-strand regions with a
composition channel regardless of their posterior DNA variance. It retains zero-count
intron/intergenic anchors, the opportunity/mass predicates, boundary and both-strand
exclusions, and the existing reliability-weight formula. No variance floor, approximate-flat
threshold, capture detector, new population or fitted biological constant was added.

Each newly admitted object contributes the same density evidence already implemented:
its observed strand columns with the RNA amount integrated out, plus existing received
composition/DNA/RNA factors in their coordinates. The graph and message builders are unchanged.
Their known approximation limits remain; admission does not certify them.

This is a **frozen final-refit contrast**, not a self-consistent new calibration:

- The count-candidate observations, incoming beliefs, domain and received messages are frozen.
- Every previously admitted curve and weight is bit-identical to the control.
- Only newly admitted objects require new numerical evidence evaluations.
- The existing population objective, convergence certificate and smoothing formula remain.
- Earlier refits and capture-reference metadata are held fixed when counts are replayed.

The changed posterior is not fed back to create a different set of training inputs in this
experiment. A successful frozen result would need that follow-up before broader validation.

## Exact implementation and receipts

All new files are in `.cache/rigel_runs/2026-10-08_evidence_admission/`:

| File | Role |
|---|---|
| `population_inputs.py` | Removes the variance cut while retaining structural predicates and the existing weights. |
| `test_population_inputs.py`, `cutoff_control.py`, `mutate.py` | Falsification control, admission/weight/immutability gates and deliberate defects. |
| `freeze_inputs.py` | Rebuilds the unchanged count candidate; verifies old observations, factor tables, weights and final count arrays exactly. |
| `evaluate.py` | Evaluates added objects, reuses every old curve, and checks pilot results against the complete calculation. |
| `fit_population.py` | Reuses the earlier certified population objective and optimizer; no new estimator is introduced. |
| `replay.py` | Replays the final population in the frozen count solver and scores region/boundary and structural classes separately. |
| `floor_attribution.py` | Crosses old/new fitted population shape with old/new pseudo-region strength to separate their effects. |
| `influence.py` | Reports old/new objects' contributions to the change in the population objective. |

The admission control failed four of eight checks before implementation. All eight now
pass and all are exercised by six specific defects. The earlier ten population-reference
gates also pass. These scratch checks do not replace the production suite.

The numerical work completed within the planned screening budget. Both pilots reproduce
bit-for-bit in the complete calculations. Both population fits meet the unchanged objective-gap
criterion, and all matched count controls reproduce their earlier results exactly.

## What the experiment establishes

The measured zero control tolerates the newly admitted uncertain curves. Its apparent
small improvement mostly comes from diluting the existing uniform pseudo-region as total
training weight increases. Holding that strength fixed removes almost all the change.
This is not evidence that the new curves resolved its false-capture cases.

The RNA-long result has a real class tradeoff. Holding pseudo-region strength fixed leaves
that tradeoff almost intact. The added introns gain appreciable support under the changed
population; the many added exons gain very little under the retained weights. This does
not prove that the weights caused the failure. Weak or misspecified local factors, uncertain
population shape and the current count readout remain alternative causes.

No new end-to-end transcript/gene accuracy, capture-weight accuracy, release integration or
whole-genome throughput is claimed. The previously validated count candidate stays frozen.
The unstranded capture-ON stratum remains deferred for 0.8.0.

## Equal weights on the same evidence

**Approved and executed, 2026-10-08.** The owner separately authorized the equal-weight
prototype with curves, objects, smoothing strength and the count solver fixed. This
experiment does not select new production weights.

The statistical argument is straightforward. An ordinary likelihood gives each sampled
object one log-likelihood contribution. A broad evidence curve already conveys weak
information. Multiplying that curve's contribution by a weight calculated from the old
posterior can count uncertainty again and leak that old prior into the new population fit.
The current formula also weights introns and exons differently in aggregate. Equal weights
are therefore a theory-based control, not a class correction or a weight tuned to the panel.

The completed A/B:

1. Reused exactly the completed wide-admission objects, curves and density grid.
2. Compared the existing reliability weights with one unit per admitted object.
3. Held the pseudo-region mixing strength at the control value, separating relative
   evidence weighting with the already-measured smoothing effect.
4. Kept the objective, numerical tolerance, earlier refits, reference metadata and count
   readout fixed. No scan, new evidence calculation or source integration was needed.
5. Checked the weighted/unweighted optima against a closed-form two-density example,
   first on the old-weight control and then with three deliberate defects.
6. Certified every fit and replayed the zero control and RNA-long condition, reading introns,
   exons and boundaries separately. The runs stayed within the seconds-to-a-minute budget.

It improves introns while worsening exons and boundaries. The zero control is effectively
unchanged. The arm stops at this screen: no intermediate powers, class multipliers or
larger release runs. Settled numbers and receipts have one home under **EQUAL OBJECT
WEIGHTS** in the ISSUES entry above.

The comparison also found that two certified solutions of the same retained-weight
objective give different downstream counts. The current objective-gap check certifies
fit quality, not stability of the population or count readout. That sensitivity needs
separate numerical investigation. The completed follow-up largely attributes it to early
stopping: tighter fits from both starts give close class scores without a model change.
Exact measurements and the remaining per-object difference have one home under
**POPULATION FIT STOPPING** in ISSUES. The equal-weight tradeoff persists. The better-scoring
initializer is not a model choice, and no new production tolerance is selected here.

The implementation is `.cache/rigel_runs/2026-10-08_equal_weights/population.py`: a small
wrapper around the existing optimizer and smoother, with a copied unit-weight vector.
`test_population.py` and `mutate.py` exercise its four contracts. `run_fit.py` saves both
certified fits and hashes of their unchanged inputs; `replay.py` checks all incoming refit
arrays between arms. `optimizer_sensitivity.json` compares the fresh and archived controls.

Equal weights do not prove independence of overlapping neighbouring evidence and do not
repair incorrect opportunities. Those remain explicit limitations. Reader integration and
selection of capture's proper local prior remain separate owner checkpoints.

## Proposed input attribution: own observations in population training

**Held after owner feedback, 2026-10-08.** This remains a possible diagnostic ablation,
not the recommended next replacement or an authorized experiment. The owner points out
that unstranded captured exons need propagated information to teach the enriched population.
Spliced-fragment observations and structurally DNA-only boundaries are direct measurements
whose information reaches regions through the existing message system. Removing all
received factors would discard those measurements along with their transport assumptions.

The question is whether adding received neighbour factors to every training curve helps
the shared DNA population. Those factors may reintroduce another region's observations
at several training objects, and retain the known face approximations. The current
population objective treats each resulting curve as another contribution. Neither
no-echo propagation nor equal object weights establishes independent observations.

The alternative training input is the already-implemented own-observation integral:

```text
L_i(rho) = integral Pois(u_i; rho Eg_i/2 + q_i r)
                    Pois(v_i; rho Eg_i/2 + (1-q_i) r) r^(-1/2) dr.
```

Keep the declared RNA reference, component opportunities, admitted regions, reliability
weights, grid, smoothing strength and earlier refits fixed. Pure-DNA and zero-count
limits remain the existing ones. The scanner deposits each accepted contained path in at
most one region, so these own region columns do not copy that contained observation into
neighbouring training curves. This does not remove possible biological dependence between
regional rates or posterior dependence of the retained weights/admission.

Only the landscape-training evaluator changes: call the existing density evaluator with
the region's own observations and no received factors. Neighbour messages remain active
in the count solver and the proposed capture-evidence consumer. The global population
continues to inform both-strand count estimates. The experiment adds no capture detector,
class prior, capture probability, cutoff or new native kernel.

If this diagnostic is separately selected and approved:

1. Write gates for equality with the independent own-observation reference, insensitivity
   to changing received factors, sensitivity to changing own counts, and preservation of
   the unchanged training mask, weights, known-origin limits and smoothing. Verify the
   old input fails the isolation gate, then deliberately break each candidate contract.
2. Price a small pilot from the saved zero-control and RNA-long inputs. Reuse the existing
   native evaluator. Stop or reduce the measurement if complete evaluation would exceed
   the agreed screening budget.
3. Fit and replay matched own/received-input populations with the same numerical settings.
   Check readout refinement rather than relying only on an objective certificate.
4. Score each class separately. A gain would implicate the neighbour-inclusive training
   input, but would not by itself distinguish repeated observations from misspecified
   messages. A loss would price discarded information. Do not compensate with class
   weights or promote a tradeoff as a release improvement.
5. Advance a useful result to the in-scope stranded and OFF controls before considering
   self-consistent refits, the prior-once partner or broader release validation.

The prior-once and reader-integration checkpoints remain. No new production training
rule follows from this proposal.

## Discussion: a conservative first solve without the intron factory

The owner considered starting all intergenic and intronic reads as DNA on unstranded
data. After discussion, the owner agrees that using all intronic counts as known DNA
to train the first landscape would produce inaccurate calls, and authorizes continued
work with message-informed evidence. The ruling's home is **Unstranded bootstrap** in
[DESIGN.md](../DESIGN.md). No hard assignment or new initialization prior is selected.

The current count-only prototype first initializes beliefs, then performs a self-solve
and message propagation without a fitted landscape, then fits the first landscape.
Intergenic objects with no annotation-admitted RNA are already structurally DNA-only.
An intron that admits RNA is not structurally DNA-only. Its unstranded total measures the
sum of expected DNA and RNA counts; it cannot identify the split alone. At exactly
uninformative strand balance and with no other factors, the existing reference for one
admitted RNA strand has median one half. That is a reference choice, not a measured split.
Removing the intron factory removes its background-density assumption, not the RNA
population or the neighbouring observations.

Changing only an initial fraction does not implement the proposed training behavior:
the self-solve and propagation intervene before the landscape is fitted. Training on all
intronic counts as DNA would change the first training inputs or solve. It would also
allow real intronic RNA to raise the shared DNA distribution and influence other objects.
The weak-strand stress trace already demonstrates this direction of feedback; it is a
counterexample to guaranteed safety, not a reason to tune the model to that fixture.

Keep three concepts separate during this discussion: a conservative starting assignment,
the evidence that can overturn it, and the strength with which it teaches the population.
A prior assigning zero probability to any RNA cannot be overturned by finite likelihood
evidence in that update. A provisional DNA assignment with uncertainty is different, and
must not become a newly observed Poisson DNA count. No new prior odds, hard strand cutoff
or intron-specific confidence is selected.

The recommended evidence interface therefore continues to include local splice, boundary
and neighbour information. The demonstrated problem is loss of allocation uncertainty
when inferred counts become observations, with shared-observation dependencies in population
training a separate concern. It is not that propagated information is inherently unusable.

The subsequent saved-factor attribution is complete in the
[weak-strand checkpoint](RNA_WEAK_STRAND_EVIDENCE_CHECKPOINT.md). It traces the extreme
fixture's pressure to existing gene-edge messages, while the moderate fixture retains
pressure from other sources. This supports keeping the local evidence audit separate
from population weighting; it does not select a new admission or propagation rule.
