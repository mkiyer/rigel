# Existing messages, honest density readout

*2026-10-07. Development implementation and review checkpoint. Production source and the
frozen count candidate are unchanged. The capture reader is not integrated.*

## What was implemented

`message_evidence.py` is a Python reference that asks how well the observed strand columns
and the existing received composition row support each proposed DNA density. It sums over
the unknown RNA amount with the already specified RNA reference. The received row is read
inside that sum, at the proposed DNA/RNA odds.

This reuses Rigel's existing local transport, opportunity maps, discrepancy treatment and
directional passes. It adds no graph, spatial smoother, capture detector, neighbor radius,
new population or biological constant. The DNA landscape is not an argument. A positive
fitted DNA count is no longer substituted for an observed DNA count in this reference.

The current scope is one admitted RNA strand. Pure-DNA and absent-component limits are
explicit. The both-strand RNA cube was audited for prior dependence but is not yet integrated
by this evaluator. The existing landscape's structural exclusions make the single-strand
interface a useful first input contrast; this is not an own-strand-only release reader.

The statistical formula and numerical tail bound have moved to
[EQUATIONS.md](../EQUATIONS.md), **Density readout of an existing composition row**.
The numerical reference splits at message-table knots and bounds the discarded integration
tail by numerical tolerance. It does not cap a Poisson expectation at the observed total.
It is deliberately a certification reference, not the whole-genome evaluation strategy.

## What passed

The initial point-count control failed 13 of 21 checks. The completed reference passes
24 checks, with numerical warnings treated as errors. Eleven deliberate mutations/reversions
and the control exercise all 24 gates. These cover the independent latent-origin sum,
zero density and opportunity, the pure-DNA limit, high-rate tails, strand and unit changes,
row normalization, preserving input arrays, and rebuilding messages after belief changes.

An existing native intron-to-boundary message is passed directly to the evaluator and checked
against the independently derived two-object integral. The table is refined to isolate
its interpolation error. This certifies that face's adapter, not every splice or level rule.

Three fresh test-chromosome calibrations then rebuild the delivered messages with and without
the DNA prior on an identical grid. Composition and present cube rows are bit-identical even
though the count estimates change substantially. Replaying the original count solve with
diagnostics is also bit-identical. The measured findings and object-level trace have moved
to [ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice).

## What this clarifies

Keeping a curve makes a substantial difference to the strength of evidence. At the worst
zero-control boundary, treating the fitted count as observed DNA nearly rules out background;
the observation/message curve retains a much less decisive alternative. Neither calculation
by itself is a capture posterior.

The remaining pressure in that curve comes from observed neighboring strand imbalances.
The trace does not justify another message system or a rule that forces zero on simulated
RNA-only libraries. A realization can favor the wrong explanation. The intended repair
preserves uncertainty so an explicit prior and consumer can respond appropriately.

The profile is still a composite approximation wherever the existing messages are approximate
or reuse physical observations. The same source is not automatically independent because
it arrived along a different path. Prior independence is one requirement; it does not certify
capture-transfer assumptions, independence, finite-grid locality or calibration of uncertainty.

## Next implementation boundary

1. Finish the frozen count candidate's panel and serial real-data validation independently.
   Its reader remains the current one; the new evidence reference has changed no reported
   transcript result.
2. Implement bounded-memory evaluation of this same contract in the isolated native worktree,
   keeping count outputs identical while the new consumer is disabled. Check numerical
   convergence against the reference before a population replay. Do not replace the graph
   or add a new sharing-strength model to make this faster.
3. Replace the remaining mixed-object point-count inputs in the previously controlled
   landscape experiment. Keep training objects, weights, grid and fitting objective fixed
   for that contrast. Read zero controls, component length gaps and captured rows separately.
   A successful input/estimator partner precedes the separately reviewed prior-once change.
4. Retain the owner checkpoint for capture's proper prior and reader integration. Neither
   the annotation-dependent `1/n` prior nor a `1/C²` tail is selected by this work.

The component-opportunity frame, DNA-to-RNA capture transfer and shared-observation questions
remain open. Small panel improvements do not settle them. No new source landing, full-suite
pass, preflight pass or release readiness is claimed for this reference.

## Reproduction

Artifacts are in `.cache/rigel_runs/2026-10-07_message_evidence/`:

| Artifact | Contents |
|---|---|
| `message_evidence.py` | Single-strand reference and analytic component limits |
| `point_control.py`, `falsification.log` | Rejected point-count shortcut and its failing gates |
| `test_message_evidence.py`, `verification.log` | Independent checks, including existing native transport |
| `mutate.py`, `mutations.json`, `mutation_*.log` | Deliberate model, numerical-interface and native dependency reversions |
| `audit_delivery.py`, `delivery_audit.json`, `*_curves.npz` | Fixed-grid prior interventions, count identity and diagnostic density curves |
| `trace_zero.py`, `zero_trace.json` | Exact delivery reconstruction and channel attribution |

Use the frozen `site_exact` wheel through `run_native.py` from the factory-removal checkpoint.
The adapter's reference tests take about two seconds; the three cached delivery audits take
about six seconds on the measured machine. These are not whole-genome performance claims.
