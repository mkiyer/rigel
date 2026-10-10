# RNA-short and capture: the 0.8.0 implementation plan

*Updated 2026-10-09 after the independent release-foundation review and owner acceptance.
Foundation integration and the reviewed golden updates are authorized. The owner drives
commits and release.*

## Objective and scope

Deliver a robust, maintainable release candidate at roughly the current accuracy,
including libraries whose RNA and DNA fragment-length distributions differ. Prefer a
small, general model over additional accuracy on the simulation panels. The owner's
rough accuracy guides are not fitted thresholds; their authoritative home is
**Release priority** in [DESIGN.md](../DESIGN.md).

First consolidate and verify the production foundation. Then finish these deliverables:

1. A compact, detector-free, continuous per-object capture reader that preserves
   uncertainty and uses the existing message graph.
2. Bounded validation across fragment lengths, probe layouts, capture strengths,
   zero-DNA controls, depth and real libraries.
3. One consolidated release candidate with complete landing checks.

Three strata remain required: stranded OFF, stranded ON and unstranded OFF. Report
unstranded ON separately; it remains deferred for 0.8.0. No pooling, capture detector,
probe BED, new fitted cutoff, unconditional landscape class split or restored intron
factory. Reader integration and a new statistical ruling still require owner discussion.
The owner drives commits and publication.

## What exists today

| Implementation | Location and status |
|---|---|
| Component-opportunity map repairs | Main working tree, uncommitted. Preserve the owner's existing source, tests and golden changes. |
| Factory-free count foundation | Integrated and verified in the main working tree on 2026-10-09. Removes the intron factory and supplies observed conditional-binomial strand messages, including eligible introns. |
| Existing capture reader | Still used in all full-panel and real-library count-candidate receipts. These results do not validate a new reader. |
| Detector-free point-count reader | Prototype only. Fails weak-evidence and zero-DNA controls by treating inferred counts as measured counts. |
| Density-evidence references | Separate Python/native prototypes. Own observations are certified; the complete neighbour factor and efficient both-strand implementation are unfinished. |
| Population, admission, equal-weight and edge-withdrawal alternatives | Measured research, not selected for production. None is part of the count candidate. |

The count worktree is `/private/tmp/rigel-count-clean-20261008`. Its current frozen package is
`.cache/rigel_runs/2026-10-09_foundation_landing/site`; the preceding cleanup control is
`.cache/rigel_runs/2026-10-08_release_foundation/site`. Main production source contains no
experimental density integrator or population alternatives to remove.

The useful count change is small: remove a duplicated intron assumption and express
own strand messages using observations rather than a posterior-frozen approximation.
It retains the existing graph, transfer maps, discrepancy rules, count solver and
landscape. It has not solved every uncertainty or capture-opportunity problem.

## Foundation integrated

**Pause new statistical experiments.** Preserve their code and findings, but do not
merge them into the selected candidate or optimize another panel average. The prior plan
is saved at `.cache/rigel_runs/2026-10-08_release_foundation/RNA_SHORT_FIX_PLAN.md.before`;
the checkpoint documents below preserve the research for review and resumption.

Only demonstrably unused code and interface baggage were removed. The completed cleanup removes
posterior beliefs and two dispersion arguments from message preparation, which no longer
reads them. The count solver still uses its own beliefs and dispersion parameters; those
stay. The Python policy now carries the single protocol parameter its messages use.
Obsolete fixture fields and comments were removed without compatibility paths. The final
landing also deletes the unused background fitter and Gaussian message helper, consolidates
duplicate row/protocol inputs, and removes fixed-zero assembly dispersion state. Low-level
dispersion parameters retain meaningful solver-test consumers.

The completed verification covered:

- The observation-only interface specification, verified failing on the old implementation.
- A native worktree build and a separately frozen, import-verified package.
- All calibration fields against fresh controls on the existing seven cached cases,
  and the existing end-to-end identity instrument. No output changed.
- The independent observation-law gates and actual compiled perturbations.
- The full integrated suite, lint and full preflight. The eight accepted golden cases were
  regenerated and checked against the archived candidate exactly; no tolerance was relaxed.
- The full-sweep belief-invariance gate and its compiled Gaussian restoration, a real-library
  identity across the entire cleanup chain, and ordinary-import end-to-end stratum checks.

The exact changes and receipts are in
[the readiness checkpoint](RNA_COUNT_READINESS_CHECKPOINT.md). Settled results belong
in `ISSUES: the-gdna-prior-enters-psi-twice`, under the count-cleanup entries.

## Review incorporated: integrate the foundation, then address the reader

The [independent review](RNA_RELEASE_FOUNDATION_REVIEW.md) supports integration of the
count foundation with the weak-strand tradeoff explicitly accepted. I agree. Do not
require a new solver or population model before integration. This is a recommendation
for a stable development baseline, not a claim that the release candidate is complete.

The owner has accepted the measured weak-strand limitation and authorized integration
plus the eight deliberate golden updates. The standing preference is
**Conservative allocation under uncertainty** in DESIGN. Goldens record expected outputs;
updating them does not certify accuracy. Preserve the truth-scored stress pair and the
reported deferred regression. Commit and push are separate owner instructions.

The review's frozen-row diagnosis is supported by a fresh array-only check. Gene-edge
information can say only that DNA exceeds a lower bound. With no opposing evidence,
the reference supplies the location and its posterior median trains the landscape.
This is a real modelling limitation under the retained rules, not broken map arithmetic.
The detailed measurement belongs in `ISSUES: the-gdna-prior-enters-psi-twice`, under
**WEAK-STRAND FLOOR DIAGNOSIS**. Its prevalence on unseen real libraries is not established.

### Completed foundation landing

1. Close the real-library identity gap against the validated count candidate, covering
   the full cleanup chain. Run one whole-genome process at a time, without debug capture.
   An earlier cleanup alone is not a sufficient control for the entire chain.
2. Add the full-sweep observation-only gate: vary valid incoming compositions, including
   endpoints, and compare delivered composition/level rows and presence masks with all
   observations and technical inputs fixed. The count posterior itself may change.
   Reintroduce the belief-dependent strand builder in a compiled copy to falsify the gate.
3. Audit the remaining cleanup suggestions and take only proved no-ops. Remove the
   redundant source-class predicate. Consolidate duplicate row inputs and protocol values
   where no live caller distinguishes them. Remove fixed-zero dispersion plumbing only
   after checking all solver, instrument and test callers. Do not add compatibility paths.
4. Apply the selected worktree source/test diff to main while preserving the owner's
   opportunity-map repairs and newer permanent docs. Rebuild; repeat identity checks on
   the actual integrated package. Do not copy stale worktree docs over main.
5. Regenerate the eight accepted golden outputs from that integrated package. Read the
   actual diff again, update DESIGN/EQUATIONS/ISSUES together under the move rule, and
   re-derive the suite count. Require no failures, lint and full preflight.
6. Reproduce a condition from each stratum, retaining deferred reporting. Present the
   concrete source, documentation, goldens and receipts for the owner's commit.

All six steps above are complete. Their receipts and the exact implementation summary are
in the readiness checkpoint. Commit and push remain owner actions.

All four available real libraries are stranded. The local cfRNA archive contains no genuine
unstranded real input. This coverage gap remains explicit; changing a protocol parameter or
scrambling a stranded library is not a substitute. Retain the synthetic unstranded checks.

### Qualifications to the review's proposed next changes

- **Reader scope:** accept a consumer-focused prototype with the count solver and landscape
  frozen. Retain licensed local neighbour evidence; the exact own-observation primitive
  alone repeats the refused own-strand-only approach at ambiguous objects. The primitive
  is a reference implementation, not an efficient complete production reader. Define its
  capture prior, normalization and no-evidence limit explicitly before implementation;
  removing the current reference and clip leaves those questions to answer.
- **Acceptance:** use roughly maintained in-scope accuracy, zero-control false enrichment,
  theoretical limits and input-variation robustness. A literal requirement for no change
  beyond A/B noise is inappropriate for a deterministic fractional comparison. A frozen
  count solver's unchanged false-DNA count is not a sufficient test of its new reader.
- **One-sided admission:** explicitly deferred by the owner until foundation integration.
  Keep the existing refusal in force. Any later experiment needs
  owner approval under `ISSUES: the-landscape-training-population-arms`. Do not assume
  that the frozen-landscape result bounds the outcome or that an endpoint maximum alone
  proves a mathematical bound; finite grid support can also produce an endpoint maximum.
- **Belief dependence:** the full-sweep gate is now present and passes at tested endpoint
  beliefs as well as interior values. The precision calculation clamps the DNA fraction
  before testing positivity; the proposed exact-vertex exception did not reproduce.
  Distinguish message availability
  from the separate posterior-based landscape admission/weighting, which this cleanup retains.
- **Background estimator:** `density_deconv` and its obsolete tests are removed because
  no selected consumer uses it. A future reader must justify its own background input;
  retaining an unused fitter did not settle that design.
- **Harnesses and documents:** after a reproducible baseline is established, stop using old
  constructor patches and import hooks for new A/Bs. Prefer ordinary isolated builds. Archive
  research with its dependencies, manifests and receipts before consolidating the active
  checkpoint documents into the plan and a ledger. Do not delete the only reproducible evidence.

The detector-free reader remains the main release requirement. Capture/length transfer is
an inherited limitation to measure with the existing cached expected-yield instrument, not
proof that the count candidate regressed. Report both absolute adequacy and changes versus
current; an equally biased baseline would not establish robustness. No new prior, numerical
range policy or reader integration is authorized by the review alone.

## Reader scope: one responsibility, complete local evidence

The count solver answers **how much DNA to allocate**. The reader answers **how capture
changes an object's contribution to effective length**. A conservative answer to the first
question is appropriate under the owner's ruling, but it is not strong evidence for the
second. If an object might contain a few DNA fragments, dividing that uncertain estimate
by a nearly zero background must not manufacture confident, enormous enrichment.

The recommended first prototype changes the reader while retaining the count model and
landscape training. This freezes the statistical model, not every helper that prepares its
inputs. It must use observations and licensed local neighbour evidence,
with their uncertainty, instead of treating a posterior point count as a fresh Poisson
observation. It reuses the existing region/boundary graph and splice information. It does
not assume equal expression at distant genes, introduce a second propagation system, or
select probes from a BED. Own-strand evidence alone remains insufficient.

| Approach | Advantage | Cost or limitation |
|---|---|---|
| Reader first, retain the count model | Small change with a clear cause; protects validated count behaviour; avoids reopening the population model to pursue panel accuracy. | Defective opportunity or message inputs still need isolated repairs; correct uncertainty cannot be recovered after the input has discarded it. |
| Change reader and count/population model together | Could repair an upstream limitation the reader cannot overcome. | More interacting assumptions, harder attribution, more validation and greater risk to already acceptable counts. |

Choose the first approach with broad robustness gates, not a restricted validation panel.
Start with analytic no-evidence, zero-opportunity and no-DNA limits, the archived false-capture
object, and a weak but real captured object. Then vary fragment laws in both directions,
probe placement/count, capture strength and depth, without retuning. Score capture against
expected yields, counts against origin truth, and transcripts/genes separately per stratum.
Report deferred unstranded ON and the lack of an unstranded real-library validation input.

If a failure persists because the local factor or opportunity is wrong, demonstrate that
cause with a frozen-input contrast and propose one upstream repair. Do not force a reader-only
solution through clipping, a detector, a threshold, or a new prior that conceals the failure.
The owner's robustness requirement decides whether to expand scope. Derive and review the
capture prior, normalization, numerical support and no-evidence readout before implementation;
the complete efficient density reader is still an open design, not a few lines ready to wire in.

### Scope made concrete after integration

A fresh ordinary-import counterexample confirms that fixing only the final consumer is too
restrictive. With all local observations, topology, opportunities and fitted protocol fixed,
an unrelated exon's expression can change an RNA message's interpretation. Separately, the
library coordinate can be zero or move the useful local profile outside its finite table.
Finer spacing does not restore a profile outside the table; a wider diagnostic range does.
These are evidence-preservation defects, not reasons to redesign the population prior.

The first correction is integrated independently: remove the counted-exon requirement on
the existing strand-protocol verdict. Source and native observation laws, coordinate origins
and hop formulas remain the same. The failed-first/mutation checks and seven calibration
identities are in `ISSUES: the-gdna-prior-enters-psi-twice`, **LOCAL-WITNESS INTEGRATION**.
The real-library pipeline is also identical. One tiny, truth-improving antisense golden move
is reviewed there; the final suite has 2,762 passes, with lint and full preflight green.

The zero-coordinate correction is integrated separately afterward: measured intron and
splice sources retain a positive numerical RNA origin when the exon-based reduction is
zero. Existing positive origins and native source rules stay unchanged. Its portable gates,
mutation coverage and seven complete calibration identities are recorded under
**ZERO-COORDINATE INTEGRATION** in the same issue. This does not implement the independent
density range or fix loss outside a finite table. The whole-library LBX0190 pipeline is also
bit-identical. After replacing the scheduler-dependent reorder gate with two deterministic
consumer checks, the actual main suite passes 2,788 tests; lint and full preflight pass.
No goldens changed in this repair.

The certified-source footprint correction follows separately. It evaluates the known
Poisson likelihood over the existing blur's footprint and retains the original output
grid. The fresh source and rounding falsifications, compiled defects, calibration census,
four test-chromosome strata and real-library comparison are recorded under
**SOURCE-FOOTPRINT INTEGRATION** in the same issue. It has a small downstream transcript
cost on the test chromosome; the moved stranded-ON assignments concentrate in ambiguous
isoform clusters and the contrast changes sign under a different existing EM start.
The shipped start stays. This is a numerical evidence repair, not a claimed accuracy gain.
The installed extension matches the tested candidate byte-for-byte; the final main suite
passes 2,796 tests, lint and full preflight pass, and the one negligible golden move is
reviewed before updating that fixture alone.

The low-probability arithmetic correction is then integrated independently, preserving
the same finite Gaussian convolution while removing its artificial likelihood floor.
Eighteen gates, five compiled defects, seven cached calibration comparisons, the changed
RNA-long condition end to end and the LBX0190 comparison are recorded under
**BLUR-TAIL INTEGRATION** in the same issue. The real-library output is bit-identical;
this correction does not resolve missing numerical support or select a capture prior.
The final integrated suite passes 2,814 tests, lint and full preflight pass, and no
goldens change. The installed native binary matches the tested candidate exactly.

The remaining range failure is freshly confirmed on that integrated build: the three
retained diagnostic specifications still fail, while wider controls reproduce the prior
convergence result. See **POST-INTEGRATION RANGE RECHECK** in the same issue and the
[representation proposal](RNA_MESSAGE_SUPPORT_DESIGN.md). This recheck changes no source
or production test. The owner subsequently authorized the isolated representation
prototype on 2026-10-09. Production integration remains a separate checkpoint.

That prototype is implemented in Python and as an isolated native numerical evaluator.
The owner also authorized local numerical units. The Python unit-selection prototype
removes the arbitrary library-coordinate dependence, with the count lattice fixed.
Convergence checks separately expose residual error from the existing density spacing.
The practical audit now exercises fixed, freshly fitted landscapes and the native count
readout: most checked counts change little, but two captured RNA-long exons have material
spacing sensitivity carried by their RNA messages. The landscape is not refitted and
the future density reader is not exercised. See the
[representation checkpoint](RNA_MESSAGE_SUPPORT_DESIGN.md) and **LOCAL-UNIT PROTOTYPE** /
**SPACING READOUT AUDIT** in the same issue. A native unit-selector port and production
consumer wiring remain undone; no numerical certification or accuracy gain is claimed.

Keep the remaining work in this order:

1. The bounded spacing audit and narrow-source diagnosis are complete. Keep the count
   foundation stable and stop broad refinement work. Retain the source-preserving
   correction as an input contract of the reader. Preserve local-unit covariance, independent
   support and the fixed count lattice. A small aggregate shift does not make those local
   failures disappear; a failed strict curve tolerance alone does not justify a general
   adaptive integrator. Do not substitute a larger fixed window, expression cap or
   fixture-selected spacing. Review any numerical contract change before a native port.
2. Review the smallest complete consumer before further numerical development. The bounded
   shared-quadrature calculation is complete and rejected for cost; its outcome is **SHARED
   CAPTURE INTEGRATION** in ISSUES. The owner's external-review packet is
   [RNA_CAPTURE_RELEASE_REVIEW.md](RNA_CAPTURE_RELEASE_REVIEW.md), including a copyable prompt.
   Keep both-strand and local neighbour evidence. Do not port the nested research oracle or
   start another integration framework because the first optimization failed.
3. Retain the owner-approved proper-prior prototype as the comparison reference, without
   reinstating the old library ceiling or annotation-dependent odds. Validate background
   estimation and uncertainty separately. Its mathematical screen does not select a
   production prior or certify robust performance across panels.
4. Run the bounded robustness screen, then the owner checkpoint before reader integration.

Any further upstream repair needs its own demonstrated defect and isolated comparison.
Do not broaden the count/population model to improve a simulated score.

This is the recommended scope for discussion: a reader with trustworthy inputs, while the
count/population model stays stable unless a separate causal diagnosis requires reopening it.
The independent-range and proper local capture-prior prototypes are now authorized;
production integration remains a separate decision. Bounds-training admission remains deferred.

### Local capture-prior prototype: mathematical screen passed, cost gate open

*Owner-authorized prototype completed, 2026-10-09. No production integration.*

The implemented reference consumes the object's density evidence and a fixed normalized
background distribution. The enrichment prior is uniform in retained background fraction,
with no upper ceiling; the continuous counterarm replaces neither the observation model nor
count inference. Its derivation, normalization, unit covariance and limits now live under
**Proper local capture reference and continuous correction** in EQUATIONS. The original
proposal's derivation has been moved there, rather than keeping a second mathematical home.
Equal odds remain a model assumption, not an estimate of probe prevalence.

The prior contrasts and complete-consumer measurements are recorded under **PROPER LOCAL
CAPTURE PRIOR** in `ISSUES: the-gdna-prior-enters-psi-twice`. The implementation and exact
reproduction commands are in [the local-prior checkpoint](RNA_LOCAL_CAPTURE_PRIOR_CHECKPOINT.md).
The proper tail and mixture odds were varied separately. Own observation curves and local
neighbour controls pass; the inherited message-support and narrow-source defects remain.
The old statistical dependencies on remote RNA and annotation count are removed from the
reader interface, conditional on fixed background and evidence. This does not certify the
locality of the complete production pipeline.

The complete nested evaluator fails the practical cost gate on both-strand objects. Stop
larger panel runs here. A correct scalar readout does not make an expensive evidence oracle
production-ready, and suppressing those objects would discard the evidence the owner wants
retained. Do not copy the nested research evaluator into production or weaken the prior to
make it faster.

The authorized cheaper calculation is now complete. Sharing the outer quadrature preserves
the selected answers but increases work; it does not repair the complete-consumer cost gate.
The numerical findings live under **SHARED CAPTURE INTEGRATION** in ISSUES. The owner has
requested an external review of the implementation decisions and the shortest release path;
the [review packet](RNA_CAPTURE_RELEASE_REVIEW.md) is ready for that purpose. Do not expand
this into another numerical engine while awaiting critique. A direct computation of the
three required averages is a question for that review, not an implemented or proven speedup.
Both RNA strands, witness exclusions, count uncertainty and distinct factor coordinates stay
in the contract. The hard source-precision cases remain mandatory input checks.

Background estimation and uncertainty remain unvalidated by this frozen-background contrast.
The exact-zero/unidentified background case is explicit in the reference; there is no rate
floor or capture-status fallback. Once the consumer passes numerical and cost checks, run
expected-yield screens across probe layouts and both fragment-length gaps, then per-stratum
transcript/gene A/Bs and serial real-library runs. Junction pricing and production integration
retain their separate owner checkpoints.

New production comparisons use main against an ordinary isolated
worktree build. Historical import hooks, packages and receipts remain archived only to reproduce
the earlier experiments. This plan is the active work list; the checkpoint documents below are
the research record, not parallel implementation plans.

## Research handoff

| Question | Preserved checkpoint and outcome |
|---|---|
| How much is count error versus opportunity/readout error? | [Reader review](RNA_SHORT_READER_REVIEW.md): archived contrasts and corrected attribution. True final weights bypass several mechanisms and do not prove count inference is the only limit. |
| What should honest density evidence mean? | [Density design](RNA_DENSITY_EVIDENCE_DESIGN.md) and [own-evidence checkpoint](RNA_DENSITY_EVIDENCE_CHECKPOINT.md): observations form the likelihood; population beliefs are not observations. |
| Can existing messages carry it? | [Message checkpoint](RNA_MESSAGE_EVIDENCE_CHECKPOINT.md) and [channel checkpoint](RNA_CHANNEL_EVIDENCE_CHECKPOINT.md): reuse the graph, keep composition and absolute-level coordinates distinct. |
| Does changing population inputs, admission or weights solve the problem? | [Admission checkpoint](RNA_EVIDENCE_ADMISSION_CHECKPOINT.md) and [weak-strand trace](RNA_WEAK_STRAND_EVIDENCE_CHECKPOINT.md): measured class tradeoffs; no selected replacement. |
| Is boundary enrichment always below the exon's average? | [Message-source checkpoint](RNA_MESSAGE_SOURCE_CHECKPOINT.md): the assumption can fail with probe placement. Blanket removal is not a robust substitute. |
| Is the full density reference ready to ship? | [Both-strand evidence](RNA_BOTH_STRAND_EVIDENCE_CHECKPOINT.md) and [support design](RNA_MESSAGE_SUPPORT_DESIGN.md): numerical support and cost remain unresolved; do not import the research implementation wholesale. |

The last frozen evidence attribution, in
`.cache/rigel_runs/2026-10-08_evidence_attribution/`, separates own observations, local
edge messages and population influence without fitting a new model. It is diagnostic
work, not an accuracy improvement or a selected production mechanism.

## Validation after the foundation is accepted

Start with frozen-input identities and analytically checkable limits, then the test
chromosome, both fragment-length gap panels and the ladder. Keep zero-DNA and weak-evidence
rows visible. Vary probe layout, fragment laws and capture strength without tuning constants
to each arm. Reuse the existing instruments; whole-genome libraries run last and serially.

Report transcript and gene errors per stratum against 0.7.1 and the current tree, with
calibration region and boundary errors beside them. The existing count-candidate coverage
is in [full-panel validation](RNA_FULL_PANEL_VALIDATION.md) and
[real-library validation](RNA_REAL_LIBRARY_VALIDATION.md). Those receipts are reusable only
for code proven numerically identical to the frozen candidate.

Production landing requires falsification and mutation coverage, resolved golden changes,
a fully passing re-derived suite, lint, full preflight and the move rule's coordinated
DESIGN/EQUATIONS/ISSUES updates. The reader's owner checkpoint precedes wiring consumers
and deleting the old reference/clip machinery. Commit, push and release remain owner actions.
