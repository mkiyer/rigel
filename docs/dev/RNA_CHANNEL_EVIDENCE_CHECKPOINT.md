# Density evidence with the existing channels preserved

*2026-10-08. Implementation review, outside the production tree. No commit or push.*

The neighbour messages already contain composition evidence, DNA-density bounds and
RNA-density bounds. The former prototype first merged these through the object's observed
total, then integrated over its unknown expected total. Those two operations do not agree.
This change keeps each factor in its own coordinate throughout the density readout.

The derivation's home is **Absolute level factors retain their coordinates** in
`../EQUATIONS.md`. Settled measurements and remaining failures are in
`ISSUES: the-gdna-prior-enters-psi-twice`, under **SEPARATE-CHANNEL IMPLEMENTATION**.

## Exactly what was implemented

All scratch files and receipts are under
`.cache/rigel_runs/2026-10-08_channel_evidence/`. C++ work remains in the isolated worktree
`/private/tmp/rigel-density-20261007`; the main source, tests and goldens were preserved.

| File | Responsibility |
|---|---|
| `channel_evidence.py` | A small adapter: evaluate the DNA factor at candidate DNA density; pass RNA-level coordinates and values into the native integral. `Level` validates its coordinate origin and interpolation table. |
| `density_evidence.h` | Extend the existing integral with the RNA-level factor. Reuse adaptive quadrature, log-ratio arithmetic, peak partitioning and the mathematical tail bound. A single piecewise-linear factor helper serves both coordinates. |
| `solve_kernel.cpp` | Extend the standalone diagnostic binding with optional RNA factor arrays; absent arrays mean a neutral factor. Validate both tables. No existing count-solve call invokes it. |
| `delivered_channels.py` | Extract already-received factors. Add composition messages, intersect bounds, preserve the side exclusions and existing delivery admission. Reconstruct the old fused row as an identity control. |
| `reference.py`, `test_channels.py`, `test_delivery.py` | Independent SciPy integral and explicit coordinate/wiring gates, beside the previous evaluator's tests. |
| `mutate_channels.py`, `mutate_delivery.py` | Deliberately compile or apply defects, check that the gates fire, and restore the source. |
| `freeze_channels.py`, `evaluate_channels.py`, `fit_channels.py`, `replay.py` | Isolate the coordinate correction on frozen population inputs and replay its effect on counts. |
| `audit_channels.py` | Read the corrected evidence on previously diagnosed test objects; verify unchanged count outputs and fixed-grid prior independence. |

The C++ source snapshots and executable package are in this receipt directory. `executed/`
preserves scripts before formatting, matching their run-time hashes. `verification.json`
records final hashes and checks.

## How the readout works

For each possible DNA density:

1. Convert that density to expected DNA incidences using the object's DNA opportunity.
2. Consider possible RNA amounts, including the amount's uncertainty in the two observed
   strand columns. Evaluate composition evidence at each DNA/RNA ratio and RNA-level
   evidence at that possible RNA density.
3. Integrate those alternatives using the previously declared RNA reference measure.
4. Apply the already-delivered DNA-level factor at the candidate DNA density itself.

An RNA coordinate origin is converted back to absolute RNA amount using its opportunity.
It supplies no new expression prior. The method does not fix RNA amount to its estimated
count, cap DNA intensity at the observed total, or fit a capture factor in this checkpoint.

At zero DNA, a nonconstant RNA-level factor stays inside the integral. The former analytic
Gamma limit remains available for constant RNA factors. With no RNA opportunity, the
pure-DNA shape is retained. At zero DNA opportunity, own observations supply no DNA-density
information, but a separately delivered DNA-level factor is still meaningful.

## What the checks establish

The conversion control was run first and failed the coordinate tests. The empty extractor
was also run first and failed its wiring tests. The corrected evaluator agrees with an
independent integral; coordinate origins and component units cancel; normalization offsets
and input arrays are preserved. Every numerical and extraction gate was exercised by a
deliberate defect. Earlier numerical tests remain green.

The old assembly reconstructs bit-for-bit on captured tables. Frozen population inputs,
admission, weights and grids match the earlier candidate, as do its final counts. The
changed population contrast retains the same objective, numerical certificate and prior
application; only affected factor coordinates change. Its unchanged curves are reused
after exact identity checks. This avoids another unnecessary full calculation.

The difficult RNA-long captured condition improves in aggregate, but introns do not.
The zero-control final selection has no mixed objects, so its unchanged result cannot
validate the capture reader. Direct inspection of its problematic exon and boundary
also shows that this correction does not remove their DNA pressure. This is a completed
interface checkpoint, not an accepted population replacement or a capture-accuracy claim.

## Limits and the next execution boundary

- Earlier face builders can already encode absolute levels through observed totals.
  Preserving final delivery does not undo those approximations. Diagnose a particular
  factor before changing it; retain the measured count candidate as the control.
- Coordinate cancellation for a held table is established. The subsequent
  [locality audit](RNA_MESSAGE_LOCALITY_CHECKPOINT.md) demonstrates a remote-expression
  gate on `split_live` and tests its small correction. A separate zero-coordinate fallback
  preserves previously erased sources; finite physical support remains unresolved.
  Fixed-grid prior independence does not settle those questions.
- Both-RNA-strand integration is still absent. The current nonempty population selection
  excludes these objects, while the eventual capture reader must support them. The
  [both-strand delivery checkpoint](RNA_BOTH_STRAND_DELIVERY_CHECKPOINT.md) now certifies
  extraction of their existing native cube factors; the numerical consumer remains open.
- The inherited posterior-variance admission and reliability weights remain fixed here.
  The subsequent isolated admission experiment is complete and reports a class tradeoff;
  [RNA_EVIDENCE_ADMISSION_CHECKPOINT.md](RNA_EVIDENCE_ADMISSION_CHECKPOINT.md) records its
  implementation and the separately approved equal-weight comparison, also complete.
- Neighbour factors may share observations. Their product remains a composite model;
  this calculation does not assert statistical independence or exact shared capture.
- Whole-genome throughput is unmeasured. Scalar diagnostic timings are not a production
  speed claim. Optimize only demonstrated costs after numerical equivalence is proved.

The separate count candidate has completed the full-panel and real-library comparisons
linked from [RNA_SHORT_FIX_PLAN.md](RNA_SHORT_FIX_PLAN.md). Its
[readiness audit](RNA_COUNT_READINESS_CHECKPOINT.md) has now run the full suite, preflight,
lint and golden review; the candidate still has reference/golden failures and a documented
weak-strand OFF concern. It has not landed. The profile-based
population partner, applying the density prior once, capture's proper local prior and
reader integration retain their outstanding design and owner checkpoints.
