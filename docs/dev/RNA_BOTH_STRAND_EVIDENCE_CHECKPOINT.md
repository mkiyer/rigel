# Both-strand density evidence: numerical reference

*2026-10-08. Standalone prototype, with no production count call site.*

The density reader can now consider both RNA strands while preserving each strand's
delivered profile. For each possible split of RNA between strands, it combines the two
profiles in the RNA-amount coordinate and reuses the existing native integral over RNA
amount. It then integrates over the split and includes the two pure-strand hypotheses.
The total-RNA measure and hypothesis weights remain those already declared by the model.

This completes a numerical reference for the both-strand consumer. It does not complete
the evidence model, population partner, capture reader or release. The existing count
path does not call this code, and no calibration or transcript accuracy claim follows.

The derivation's home is **Both-strand density evidence retains the hypothesis support**
in `../EQUATIONS.md`. Settled verification findings are **BOTH-STRAND DENSITY REFERENCE**
in `ISSUES: the-gdna-prior-enters-psi-twice`.

## Exactly what was implemented

All files are under `.cache/rigel_runs/2026-10-08_both_evidence/`.

| File | Responsibility |
|---|---|
| `both_evidence.py` | Merge the two level tables at a fixed RNA split; reuse the isolated native one-amount integral; integrate the split and combine its three hypotheses. `delivered_evidence()` also passes the existing support bits. |
| `reference.py` | Independent SciPy integration in the opposite order: RNA split first, RNA total second. Uses the original two Poisson columns directly, without the native integral. |
| `test_both_evidence.py` | Observation-measure, coordinate, normalization, missing-factor, narrow-profile, zero-opportunity, support and input-preservation checks. |
| `mutate.py` | Alter the actual readout and verify targeted checks fail, then restore it. Includes separate challenges to hypothesis support and numerical integration. |
| `verify.py` | Validate receipts, executed-source snapshots and preserved main-tree/native identities. |

The initial reference used the frozen isolated native build in
`2026-10-08_level_support/site`. The subsequent exact-reuse optimization below changes
only the isolated inner implementation. Main source, production tests, goldens and the
installed library remain unchanged.

## A support rule the initial formula missed

The current count solver carries more than finite profile values. A delivered RNA
profile on one strand also rules out the hypothesis that all RNA is on the other
strand. That is the existing rule **The AMBIG tilt's hypothesis space is {pure +,
pure −, mixed} — the tilt atom**, not a newly introduced detector or cutoff.

The initial generic integral retained finite values at zero RNA and therefore did not
encode this exclusion. A separate falsification exposed the mismatch. The generic
numerical evaluator now accepts explicit witness bits; its delivery wrapper derives
them from the native profile-presence flags. The count solver itself is unchanged.
Generic finite factors can still be tested independently of a presence assertion.

Do not remove this inherited rule as an incidental numerical simplification. Its merits
and the status of uncertain own lower bounds are a separate model question. Preserving
the current rule is not proof that every delivered witness is an independent observation.

## Numerical checks and scope

The own-observation control fails the factor checks before implementation. The reference
then agrees with independent reverse-order integration and with the existing finite
latent-allocation calculation where profiles are neutral. The two RNA profiles retain
their separate coordinates; additive constants are preserved across all hypotheses.
Pure-strand, zero-opportunity, missing-profile and strand-mirror limits are explicit.

Narrow-profile tests vary their centers and widths independently of any simulated panel.
Integration uses the supplied profile knots, including where their RNA-amount ordering
changes with tilt. It does not select a fixed number of tilt nodes to fit panel accuracy.
Its tolerance allocates error to the existing inner quadrature, its tail and the outer
quadrature. Refinement tests check the requested density curve, including rates above
observed total divided by opportunity.

These are diagnostic numerical checks, not a global convergence theorem for every
possible table. The generic witness-free checks and the delivery-support checks are
distinct. Executed snapshots and the verification receipt preserve that distinction.

## Cost checkpoint and revised next step

The saved-table cost screen is complete. Its settled measurements live under
**BOTH-STRAND NUMERICAL COST** in the same ISSUES entry. The implementation and receipts
are in `.cache/rigel_runs/2026-10-08_both_cost/`:

- `census.py` reconstructs the existing native messages on frozen inputs and saves
  numerical factors, without running calibration or fitting a population.
- `time_one.py` applies a bounded diagnostic timer to one density at one object. It
  distinguishes native inner work from adapter work and records unfinished probes.
- `compact.py` and `run_compact.py` arm exact deletion of constant interior knots.
  The original reference is unchanged. The profile object and its witness status remain.
- `test_compact.py`, `mutate.py` and `audit_tables.py` check representation identity,
  deliberate breakages, actual saved tables and floating-point tilt intervals.

Compaction is a useful exact simplification, but it does not make the reference suitable
for production. Native inner evaluations dominate the remaining cost. The table-derived
tilt partition also contains intervals that cannot be subdivided in floating point, while
the outer reference asks for separate relative accuracy in every interval. A direct C++
port of that scheduling would preserve the wrong computational cost.

The shared-error prototype is now implemented and tested in
`.cache/rigel_runs/2026-10-08_tilt_budget/`. Its derivation is **Tilt variation can be
bounded without fitting RNA** in EQUATIONS. Its measurements, falsifications and rejected
single-call SciPy experiment are **BOTH-STRAND SHARED ERROR BUDGET** in the same ISSUES
entry. The retained code uses adaptive Simpson integration, including interval-edge
values, with one total error budget. It retains mass and a derived uncertainty for
intervals that floating-point arithmetic cannot subdivide. It uses no new capture
assumption, density cutoff, probe information or fitted RNA value.

| File | Responsibility |
|---|---|
| `bounds.py` | Bound tilt-dependent log-likelihood changes uniformly over unknown RNA amount, using observed columns and profile slopes/ranges. |
| `tilt_integral.py` | Adaptively refine the largest estimated error; account for indivisible intervals and re-sum before accepting the total. |
| `budget_evidence.py`, `run_budget.py` | Substitute only the outer numerical calculation into the frozen diagnostic reader. Retain the same compaction, conditional native integral, atoms and support bits. |
| `test_bounds.py`, `test_integral.py`, `mutate.py` | Explicit likelihood checks, narrow peaks, floating-point limits, normalization and actual implementation breakages. The frozen evidence checks run through the new adapter too. |
| `inspect_intervals.py`, `refine.py` | Inspect saved conditional values and tighten tolerance across representative density curves, without calibration or EM. |

The rejected global-SciPy implementation and its tests are retained as `.txt` snapshots.
Tests specific to its workspace API were replaced with numerical-behavior checks when
that implementation was rejected; the narrow-peak failure was retained. No failing
accuracy gate was weakened to accept the retained method.

An additional audit exposed cancellation of the outer Simpson error estimate. The
active diagnostic runner is now `2026-10-08_outer_refinement/run_guarded.py`; it retains
the inner integral's existing check of two successive refinements. `guarded_integral.py`
is the isolated outer replacement, `test_aliasing.py` supplies exact positive-polynomial
falsifications, and `screen.py` compares saved-object values and cost. The prior
`run_budget.py` remains a frozen comparison control. Settled results are **BOTH-STRAND
REFINEMENT CANCELLATION** in ISSUES, with the derivation **One embedded difference can
cancel** in EQUATIONS.

The mathematical envelope is distinct from ordinary quadrature's estimated error.
Independent integration and refinement provide additional checks; neither proves
convergence for every conceivable input profile. The improvement still leaves seconds
of work per density on the more expensive objects. Native inner evaluations dominate,
so a translation of this Python scheduler alone will not resolve production cost.
The read-only native call census in `2026-10-08_inner_cost/` identifies repeated inner
evaluations; **INNER-EVALUATION COST CENSUS** in ISSUES holds the measurements. The
follow-up is implemented in the isolated worktree: the recursive error check passes its
already evaluated quarter points into the children. It adds no cache or tunable and
changes no quadrature budget, tolerance or approximation. **INNER-SAMPLE REUSE** in
ISSUES holds the falsification, mutation, numeric-identity, cost and count-path results.

The reuse checkpoint uses `2026-10-08_outer_refinement/run_guarded.py`
with the native library in `2026-10-08_inner_reuse/site`, now a frozen comparison
control. `prepare.py` instruments the actual native header; `test_reuse.py`
compares values, unique sampled coordinates and complete refinement traces;
`mutate.py` recompiles deliberate breakages. `time_pairs.py` records serial scalar
comparisons, and `count_identity.py` checks the seven cached calibration conditions.
The single-strand and both-strand reference families run in separate Python processes
because their frozen oracles share the import name `reference`.

The reuse improves cost but does not make the consumer suitable for production. The
read-only interval census in `2026-10-08_inner_budget_cost/` preserves the existing
values and call counts while recording work, mass and estimated error separately.

The first shared-budget attempt allocated absolute error in proportion to interval
width. It passed the analytic cost gates but failed the independent narrow-profile
test: a nearly coincident interval demanded precision beyond floating-point subdivision.
That implementation and its failing receipt are retained as a rejected experiment.

The corrected total-error scheduler is implemented and tested in the native worktree.
Its derivation is **A total integral has a total numerical budget** in EQUATIONS;
**INNER TOTAL-ERROR BUDGET** in ISSUES owns the census, rejection and measured results.
The active combination is now the same guarded outer runner with the native library
in `2026-10-08_inner_budget/site`. The implementation replaces recursive local budgets
with a heap of active intervals, preserving the same integrand and error checks.
`prepare.py` builds analytic Poisson-integral and saved-factor probes from the actual
header. `test_budget.py` and `mutate.py` challenge cost, mass and the demonstrated narrow
case. Independent numerical families still run separately. `time_pairs.py`, `refine.py`
and `count_identity.py` keep runtime, numerical accuracy and count-path isolation distinct.

This comparison adds no tunable and does not assert a global quadrature convergence
theorem. Scalar timings alone still cannot establish whole-genome production cost.

The next contrast is now implemented in `2026-10-08_lazy_tilt/`: `lazy_integral.py`
reuses the uniform tilt envelopes before constructing each full Simpson stencil.
Its derivation is **Envelope-first quadrature retains mass** in EQUATIONS. The source
and receipts retain the initial stalled accounting variant; the correction sums live
intervals directly. `test_lazy.py` checks retained mass, narrow peaks, avoided samples
and cancellation, while `mutate.py` challenges the actual implementation. The independent
reference checks run through `run_lazy.py`.

The active diagnostic combination is now `2026-10-08_lazy_tilt/run_lazy.py` with the
unchanged native library in `2026-10-08_inner_budget/site`. `run_guarded.py` remains the
frozen comparison control. `2026-10-08_complete_curves/curve.py` evaluates an entire
existing landscape grid, records its actual integrator and compares with tighter
integration. **COMPLETE DENSITY-CURVE COST** in ISSUES owns the full-grid results,
scalar timing pairs and the rejected compression screen.

The complete curve is numerically consistent but remains too slow for a production
consumer. Avoid a large calibration run or a direct outer-loop port at this cost.
The next work returns to source-factor certification, including how source likelihood
tails survive row operations. Keep the pending numerical-support proposal separate;
this checkpoint does not authorize its representation change.
Any later reduction in tilt sampling must retain the narrow-peak challenge and be
controlled by the requested numerical error, with complete density-curve comparisons.

The first producer-side arithmetic check is now complete in
`2026-10-08_blur_tails/`. The native worktree retains deep likelihood tails in the
existing finite Gaussian convolution; **BLURRED-LIKELIHOOD TAILS** in ISSUES holds the
measured results and **A convolution preserves relative likelihood in its tails** in
EQUATIONS gives the identity. `test_blur.py` uses an independent discrete log-space
reference; `mutate.py` compiles actual altered headers, and `count_screen.py` records
fresh oracle-scored pairs. The prior native library remains frozen in `inner_budget/site`;
the corrected experimental library is in `blur_tails/site`. Existing consumer receipts
continue to refer to their original native library and original saved factors. This is
not a new full-curve timing or a claim that the producer model is certified.

The subsequent direct-edge audit is recorded in
[RNA_MESSAGE_SOURCE_CHECKPOINT.md](RNA_MESSAGE_SOURCE_CHECKPOINT.md). It isolates a
source already projected through the recipient's observed total, before final channel
separation. The small adapter corrects direct recipients only; native traces also find
forwarded copies. Do not substitute those partially corrected curves into a population
fit. Preserve source meaning before its count projection and certify the encountered
transfers independently of this numerical consumer.

## Remaining work before production

- Complete the numerical cost/error checkpoint above before selecting the native
  consumer. Certify that consumer against this independent reference, keeping the
  count path unchanged; this scalar Python reference is not a throughput claim.
- Complete finite-support repair after the pending representation decision. Integrating
  a truncated input accurately does not restore discarded evidence.
- Retain the separate audit of inherited face approximations and reused observations.
  The new integral preserves those factors; it does not make them a joint likelihood.
- Keep population admission, weights and prior-once separate. Both-strand objects are
  not silently added to population training by making their evidence computable.
- Integrate the detector-free capture reader only at its existing owner checkpoint,
  followed by the full per-stratum, zero-control, probe/length and real-library gates.

No new capture prior, probe rule, distant-expression input, density floor, fitted-DNA
Poisson observation or landscape-weight decision is introduced here.
