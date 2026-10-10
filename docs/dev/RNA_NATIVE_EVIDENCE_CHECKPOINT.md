# Native density evidence: implementation checkpoint

*2026-10-07. Isolated native prototype. No main-tree source, tests, goldens or installed extension changed.*

The evidence formula now has a native evaluator. It answers how well each candidate DNA
density explains the observed strands and the existing received composition row, summing
over the unknown RNA amount. It accepts no fitted DNA count or DNA prior. This is the
same model as [the Python reference](RNA_MESSAGE_EVIDENCE_CHECKPOINT.md), implemented
for numerical testing and the next frozen population-input contrast.

## Exactly what is implemented

In `/private/tmp/rigel-density-20261007`:

- `src/rigel/native/density_evidence.h` implements the single-RNA-strand integral.
- `src/rigel/native/solve_kernel.cpp` adds an independent `density_evidence` binding.
  No existing count-solve call invokes it or has changed its arithmetic.

The Python adapter, independent tests, source snapshots, isolated package, build logs and
receipts are in `.cache/rigel_runs/2026-10-07_native_evidence/`.
The frozen count candidate remains separately installed in its original scratch package.

The evaluator takes observed integer strand incidences, component opportunities, strand
fidelity, a DNA-density vector, the existing log-odds message grid and a numerical tolerance.
It returns log evidence, retaining arbitrary density-independent message normalization.
RNA amount is integrated under the already-declared `r^(-1/2)` reference. RNA opportunity
enters the transported message maps; DNA opportunity converts density to expected counts.
A positive RNA opportunity is not an admission threshold.

With no RNA opportunity, the evaluator uses the analytic pure-DNA likelihood. At zero DNA
rate it uses the analytic RNA-only limit. It supports rates above the observed total, and
zero DNA opportunity gives a flat function of DNA density. Received rows are read-only.
The default scratch tolerance is a numerical request, not a model parameter.

For the mixed integral, the change of variable `r=t²` removes the reference singularity.
The integrator splits at the received row's interpolation knots and at interior likelihood
maxima. Adaptive Simpson integration checks two subdivisions before accepting an interval.
A derivative bound controls the remaining infinite tail; no biological upper density is
introduced. Log-rate ratios cancel the common large Poisson terms before subtraction.
The derivation is in `EQUATIONS.md`, **Density readout of an existing composition row**.

The native call releases the Python lock, so independent objects can be evaluated concurrently.
A two-worker comparison on twenty frozen objects is bit-identical to serial evaluation.
This is still a numerical prototype: no whole-genome evidence array or consumer is introduced.

## Falsification and verification

The first tests were run against the frozen module without the new interface: **38 failed,
1 passed**. The passing test was the already-existing rebuilt-message independence check.
The first implementation exposed two numerical defects before any consumer change:

- A single Simpson error difference could cancel. Two low-count, unstranded cases missed
  their declared numerical tolerance. Requiring the second subdivision fixes this.
- Direct subtraction of huge negative log probabilities exhausted subdivision at a DNA
  rate of `10^12`. Algebraically cancelling common terms fixes it without a rate cutoff.

The completed suite has **42 passing checks**. It compares against latent-origin sums,
independent adaptive integration and the separately derived two-object intron/boundary
integral. Cases include small and concentrated counts, both strand orientations, unstranded
observations, uneven or multimodal rows, normalization, absent opportunities and far tails.

**Twelve compiled C++ defects/reversions exercise every gate:** wrong density units;
forbidding RNA-only explanations; evaluating the message at a fixed RNA amount; losing
row normalization; integrating RNA when its opportunity is absent; wrong strand fidelity;
truncating DNA rates at observed counts; modifying the input row; reverting the numerical
error check; unstable Poisson subtraction; accepting invalid inputs; and restoring the
belief-frozen Gaussian strand claim. Every mutant is rejected; the source is restored.
`mutations.json` gives the failures and binary fingerprints.

Diagnostic count replays and frozen-training snapshots preserve the previous count-candidate
arrays bit-for-bit. The earlier fixed-grid prior-independence audit is also reproduced.
This checks the unchanged count path; it does not certify all message factors as independent
physical observations or remove their existing approximation assumptions.

## Consumer experiment and interpretation

The next input contrast is in `.cache/rigel_runs/2026-10-07_mixed_inputs/`:

1. Reproduce the frozen factory-free candidate's training selection, beliefs, grid and weights.
2. Keep the certified pure-DNA/empty-row input control unchanged.
3. Replace only the remaining nonempty, single-strand mixed rows with observation/message evidence.
4. Fit the same mixture objective, with the same weights and numerical optimum certificate.
5. Replay the final landscape and score region and boundary counts separately against truth.

The first decisive comparison uses only the final frozen refit: the count solve resets at
each refit, so earlier frozen landscapes cannot affect its final output. This saves redundant
integrations. It is an attribution experiment, **not a self-consistent new calibration**.
The zero-control final selection contains no mixed rows; its unchanged result is therefore
expected and is not proof that the eventual reader is safe on zero-DNA libraries.

The test-zero and RNA-long unstranded-ON snapshots preserve the original count outputs and
training arrays exactly. The larger final refit contains 1,159 mixed rows on 260 density
points. A small timing sample was used before starting that calculation. Production-scale
throughput remains to be established; scalar diagnostic timings are not a whole-tool claim.

## Frozen-input result and the completed coordinate follow-up

*Follow-up, 2026-10-08.* The frozen mixed-input comparison and its separate-channel
successor are complete. Neither qualifies the population model for production: the
aggregate count errors improve, but the intron regression remains. Their settled
measurements are in `ISSUES: the-gdna-prior-enters-psi-twice`, under **MIXED-INPUT AND
COORDINATE FOLLOW-UP** and **SEPARATE-CHANNEL IMPLEMENTATION**.

The held-density counterexample established that the fused composition row substitutes
an observed total into a density integral over expected RNA amount. Its derivation is
**Absolute level factors retain their coordinates** in `../EQUATIONS.md`.
[RNA_CHANNEL_EVIDENCE_CHECKPOINT.md](RNA_CHANNEL_EVIDENCE_CHECKPOINT.md) records the
implemented correction: keep existing composition, DNA-level and RNA-level factors in
their coordinates. Existing diagnostics already expose them. No new graph, propagation
pass, export subsystem or biological prior was needed; the count path stays unchanged.

## What this does not settle

The existing graph, message schedule, transport maps and discrepancy treatment are retained.
The evaluator does not make every inherited splice/level approximation exact. Shared-fragment
dependence and support/coordinate invariance remain audits, not assumed guarantees.
Both-RNA-strand cube integration is not implemented here; the frozen nonempty landscape
training selection already excludes those objects, but the eventual reader must address them.

The evidence-curve population partner is not yet accepted. Posterior-dependent selection and
weights are deliberately held fixed in this contrast, then tested separately. Applying the
DNA prior once remains a reviewed model change. Capture's proper local prior is still
unselected, and reader integration retains the owner's checkpoint.

The distinct count-candidate validation now includes [every panel stratum](RNA_FULL_PANEL_VALIDATION.md)
and [four real-library pairs](RNA_REAL_LIBRARY_VALIDATION.md). This numerical work claims no
new transcript accuracy. Production landing still needs the full suite, preflight, lint and
review of intended golden changes. Nothing was committed or pushed.
