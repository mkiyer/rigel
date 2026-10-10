# Density-message support: proposed representation checkpoint

*Proposal opened 2026-10-08; rechecked against the integrated foundation and authorized
for isolated prototyping by the owner on 2026-10-09. Production integration remains a
separate checkpoint.*

I recommend prototyping a numerical range for density-message tables that is independent
of the composition/count range, at the same existing spacing. The single composition
lattice stays. The owner has authorized this representation experiment. It is not a new
capture or spatial-sharing model.

## Which table, which samples, and which algorithm

This is calibration's representation of neighbour evidence, before transcript assignment.
The scanner has already counted fragments in regions, boundaries and splice junctions.
Calibration builds evidence about the possible RNA/gDNA mixture at each object, sends
permitted messages between adjacent objects, and combines the delivered evidence with the
object's own observations and the landscape prior to infer its counts. The EM uses those
results later when assigning RNA to transcripts.

A density message asks: **how compatible is each possible density with the source's
observations?** Density means component count divided by that component's opportunity.
The opportunity accounts for the fragment lengths and geometry that allow a fragment to
be observed at that object. A message can convey a full curve or a one-sided constraint,
depending on the existing face rule.

The computer represents the curve as a row of log evidence scores at a finite set of
candidate densities. Those numerical evaluation points are the "samples". They are not
reads, patients, genomic positions or a new biological sampling procedure. Values between
them are obtained by interpolation. There are separate gDNA, RNA+ and RNA− density lanes.
This row is also distinct from the composition grid, whose candidates are RNA/gDNA
fractions, and from the shared landscape, which supplies a prior to the count solve.

The level coordinate is `log(density / numerical_unit)`. For an exact function the unit
only changes the axis labels and cancels when the recipient converts back to physical
density. The finite table does not automatically have that property: changing the unit
while keeping the numerical ticks fixed moves the physical densities that were evaluated.
Previously, observations elsewhere in the library could move that unit. That was a
numerical dependence, not an intended source of biological evidence for this object.

## The problem in plain language

A message records how plausible different RNA or DNA densities are. Its table axis is
possible density, not genomic position. At present that table is forced to use the same
finite numerical interval as the count-composition table. Changing the library's density
coordinate can move local evidence outside that interval. The table then repeats an edge
value, and the recipient can see a constant instead of the evidence.

Covering the original sources is necessary but insufficient. Each uncertainty blur needs
values beyond the requested output interval. Repeated propagation can exhaust the extra
range even when every source was represented correctly. This issue concerns the numerical
implementation of the existing propagation rules; it does not establish that those rules
are an exact biological model.

The settled diagnosis and measurements belong to **SOURCE-RANGE LIMIT** within
**MESSAGE-LOCALITY AUDIT** under
`ISSUES: the-gdna-prior-enters-psi-twice`. The operator identity belongs to
**Repeated blurs require the propagated footprint** in `../EQUATIONS.md`.

## What was tried and why it is not a release candidate

The scratch implementation is `.cache/rigel_runs/2026-10-08_shared_support/support.py`.
It derives an initial shared range from the existing source fraction limits, Poisson
likelihood tails, and certified flux's actual hop blur. It preserves the configured
spacing; no capture parameter or simulated truth enters the calculation.

| Check | Outcome |
|---|---|
| Original coordinate/unit/range specifications | Current range fails 13/13; source-derived range passes 13/13. |
| Certified sources at both junction ends and on both strands | Raw-source range fails 17/25, including an unnecessary no-source range change; adding the existing blur footprint and retaining inert ranges passes 25/25. |
| Seven frozen condition inputs | The shared range substantially increases every count table's allocation; no runtime benefit is established. |
| Repeated low-count propagation | Two new specifications fail the existing one-percent row target. Wider independent controls confirm the error. |

The last measurement uses an unstranded chain with a certified RNA source and balanced
low-count recipients. No capture mechanism or simulator yield is involved. Half-windows
27, 54, 108 and 216 are convergence controls, not proposed defaults. The fixed physical
query interval includes the recipient's useful RNA-density range.

The convergence table is recorded once in the ISSUES entry above. Its measurements are
differences between numerical representations of the same message, not transcript/gene
errors. No calibration or release benchmark was launched for the rejected range proposal.
A larger fixed window is not selected.

The post-integration recheck confirms that the separate source-footprint and blur-tail
repairs have not closed the range failures. All three saved specifications remain red
against the current native build; the convergence controls reproduce the earlier result.
The settled measurements are in **POST-INTEGRATION RANGE RECHECK** in the same ISSUES
entry. This check imports production code directly, without an experimental density
consumer or model hook. Receipts are in
`.cache/rigel_runs/2026-10-09_range_recheck/`.

## The proposed change

Keep the composition/count lattice, its reference measure and its configured spacing.
Give level-message rows enough density range for their sources, operations and consumers.
The first prototype should use intervals of the same uniform lattice, with explicit
origin and extent, and retain the existing graph, source rules, hop prices, intersections,
no-echo rules and propagation schedule.

For an output interval `[a,b]`, a finite blur with radius `h` cells requires input over
`[a-h*d,b+h*d]`, where `d` is the existing spacing. Required intervals can therefore be
propagated backward along the existing directional dependency chain before evaluating it
forward. Empty forwarding adds no blur. This numerical preparation creates no new
biological edges and no genomic radius. Its cost must be bounded to the relevant connected
neighborhood instead of widening every count solve in the library.

The reference implementation must additionally handle the lower-side prefix maximum,
row normalization, and the limits at zero density. A local crop cannot silently redefine
those operations. Their treatment must be derived and tested before a native port; the
simple finite-convolution interval identity alone is not a complete algorithm.

## Prototype outcome and the next decision

The authorized Python reference and isolated native numerical evaluator are implemented.
They repair the original missing-range examples with the count lattice fixed. The
normalization/prefix-maximum derivation is now **A finite core certifies level
normalization** in EQUATIONS. The 70 passing specifications, mutation checks and cost
measurements are recorded in **INDEPENDENT-RANGE PROTOTYPE** under
`ISSUES: the-gdna-prior-enters-psi-twice`. These are reference implementations; the generic
operation graph is not proposed as a second production message engine.

A broader falsification rejects **range alone** as the finished representation: moving
the numerical origin by a fraction of a cell changes physical sample positions. The new
three-test failure and the control that translates the knots are recorded in
**COORDINATE-PHASE LIMIT** in the same issue. The 70 green checks must be read beside those
three red checks. No calibration, transcript or gene improvement is established.

The owner authorized the following local-unit experiment on 2026-10-09. The Python
prototype is now implemented: one **local numerical unit per connected component of an
existing density lane**, shared by all its source and received rows:

- Use the lane's existing permitted faces to define connectivity, independently of block
  packing, message direction, query order and which consumers request output. No new
  neighbour relationship is introduced.
- Reduce the already admitted positive source count/exposure pairs once as `sum(n)/sum(E)`.
  This only supplies a positive numerical unit. It is not an RNA abundance estimate, an
  observation supplied to the count solve, a new prior, or an upper bound on density.
- Express table coordinates relative to that local unit. An unrelated library-origin
  change then only relabels the table; a change of exposure units relabels it covariantly.
  Keep the existing spacing, count lattice, biological sources, prices and propagation.
- With no admitted source, keep no row. Verify all-DNA/zero-RNA limits, empty forwarding,
  both strands, local-source conflicts, disconnected-component additions, block packing and
  opportunity-unit changes before any count A/B.

The advantage is removal of the arbitrary remote-origin dependence without extra
resolution levels or a capture parameter. The cost is a component label and a small local
reduction; changing local observations can still move its numerical sampling phase, so
this does not promise exact continuous-function evaluation. It needs convergence and
count checks, not a claim of improved simulated accuracy. A native production port should
use supported rows in the existing typed lanes, with one propagation implementation;
the Python operation graph remains the independent numerical reference.

**Outcome:** the arbitrary-library-coordinate tests now pass. The full local-unit
contracts and deliberate defects are recorded in **LOCAL-UNIT PROTOTYPE** under
`ISSUES: the-gdna-prior-enters-psi-twice`. The same checks pass through the existing
standalone C++ evaluator, while unit selection remains Python. No native unit-selector
port or production integration was done.

The convergence check exposes residual error from the existing density spacing; see
**SPACING SENSITIVITY** in that issue. The local unit removes the arbitrary remote-origin
dependence but does not make a finite table exact. A diagnostic through the current
one-object count solve finds a much smaller change in the inferred fraction. That check
has no fitted landscape and does not validate a future density/capture readout. Keep the
new convergence failures visible rather than calling the representation certified.

The owner-authorized bounded sensitivity audit is complete; its measurements and
limitations are recorded once under **SPACING READOUT AUDIT** in the same ISSUES entry.
It uses fresh fitted priors and preserves the native count calculation, verified by an
independent legacy-table replay. Most checked count readouts are insensitive. Two captured
RNA-long exons are exceptions: their coarse prototype allocations are materially wrong,
and refining one incoming RNA-level row accounts for the change. The finer reference is
checked again at those objects. This is a numerical prototype defect, not an observed
regression in the unchanged production build.

**Recommendation after measurement:** stop broad refinement experiments. Keep the count
foundation stable. The follow-up localizes the effect to interpolation of a narrow known
source between stored points; **NARROW-SOURCE INTERPOLATION** in ISSUES records the
fine-curve resampling contrast and independent depth-based falsifications. Preserve the
known source curve through its existing operations at the consumer as the next correction
to prototype; that correction is not yet implemented. Preserve the failures as consumer
checks. Do not turn the
diagnostic fine spacing into a new default or tune a tolerance to their truth. A small
aggregate count shift cannot certify the future capture readout; count-posterior density
moments also change, and the prior-free density reader was not exercised. The existing
strict convergence specifications remain visible rather than having their tolerances
relaxed to pass this sample. Production integration and capture-prior selection remain
separate decisions.

"Local" here means connected by the existing permitted lane, which can cover several
regions. The unit reduction does not pool those regions' biological abundances or require
equal densities. Connected observations can still move its sampling phase; disconnected
observations cannot, conditional on a fixed count grid, source admissions and fitted strand
protocol. The ordinary neighbour evidence and the separate shared landscape still have
their existing roles.

This changes the representation contract quoted below. It does not revive the refused
coarse/fine composition grids or interpolate between two composition readouts. There is
one composition lattice and no new density-resolution parameter.

## Authorized execution plan

1. Build a small Python reference for supported level rows and their existing operations.
   Make the source, prefix-maximum, normalization and tail contracts explicit. Use the
   existing spacing and numerical error budget; introduce no tunable biological limit.
2. Make the two propagation specifications pass beside the 38 source/coordinate checks.
   Extend the gates to both level types, empty relays, conflicting sources, zero-component
   limits and both opportunity-gap directions. Compare with independently widened
   evaluations, including checks that the reference itself has converged.
3. Deliberately mutate the representation and its wiring. Every new gate must reject an
   actual defect before native work begins.
4. Prototype the smallest native representation change in the existing isolated worktree.
   Keep the count solver's calculation unchanged when it receives identical factors.
   Require invariance to block packing and thread count.
5. Price working memory and time on cached test-chromosome inputs before larger panels.
   Then screen per-object counts and the weak-strand OFF stress case. Broaden to panels
   only after the numerical and count gates pass. Retain the separate reader-integration
   owner checkpoint and all production landing requirements.

## Open implementation concerns

- Each curve must carry its own coordinate metadata through the consumer. Today's
  both-strand cube shares one RNA origin; a future port must correctly convert both
  independently positioned RNA lanes before combining them. The prototype has not wired
  this interface into the count solver.
- Numerical unit covariance does not establish convergence at the existing spacing.
  The practical audit found a specific RNA-row failure alongside small aggregate count
  sensitivity. Coarse interpolation now accounts for that failure. Preserve the known
  source at the consumer and price the complete implementation before adding a
  refinement policy; the current audit does not refit the landscape or run EM.
- Variable row extent adds metadata and touches native interfaces. The prototype must
  establish that this is simpler and cheaper than enlarging all state tables.
- Prefix maxima and normalization need support outside a recipient's immediate queries.
  Their semantics must be preserved rather than approximated by a new endpoint rule.
- The exact density consumer integrates RNA amount over an unbounded domain. Its existing
  tail bound and the message's tail representation must agree, including at zero.
- Support is a numerical issue. Fixing it does not resolve the weak-strand OFF population
  feedback, the remaining likelihood approximations, capture-opportunity transfer, or the
  class tradeoff found in the separate landscape-weight experiment.

## Owner authorization and remaining checkpoint

`DESIGN.md` §6b.13 currently requires **“ONE representation everywhere — profiles on the
solve grid, no Gaussian summary anywhere in the transfer policy”**. The one-composition-
lattice ruling in §6b.15.10 remains binding. The owner authorized the isolated
representation prototype on 2026-10-09; see **Density-range prototype authorization** in
DESIGN. The same date's **Local density-unit prototype authorization** extends the
experiment to table placement. Keep the single composition/count lattice and existing spacing. A production
landing and the capture-prior decision still need their own checkpoints.

Receipts: `before.xml`, `source_only.xml`, `flux_before.xml`, `source_covered.xml`,
`range_census.json`, `propagation.json`, `propagation_curves.npz` and
`propagation_red.xml` in the scratch directory. The retained red specifications are
intentional; this proposal is not reported as a passing implementation.
