# Existing evidence at objects admitting both RNA strands

*2026-10-08. Prototype interface checkpoint. No production source, native worktree,
installed library, golden, commit or push changed in this checkpoint.*

The previous density adapter rejected an object admitting both RNA strands. Rigel already
delivers two RNA-density profiles at such objects. The new adapter reads those existing
profiles without flattening them into composition or discarding one strand. This removes
an interface gap; it does not yet implement the integral that consumes those profiles.

Settled findings and the condition-by-condition presence census live under **BOTH-STRAND
DELIVERY** in `ISSUES: the-gdna-prior-enters-psi-twice`. The wider implementation sequence
is [RNA_SHORT_FIX_PLAN.md](RNA_SHORT_FIX_PLAN.md).

## Exactly what changed

All new executable files are in `.cache/rigel_runs/2026-10-08_both_delivery/`.

| File | Responsibility |
|---|---|
| `both_delivery.py` | `BothDelivery` indexes native cube rows by slot once. Its `at()` returns composition, DNA, RNA+ and RNA− factors in their own coordinates. |
| `test_both_delivery.py` | Checks each channel, both own lower-side bounds, missing channels, row identity, density origins and input preservation. |
| `mutate.py` | Changes the actual adapter, exercises the checks and restores the original file. Retains the initial ineffective-fixture receipt. |
| `census.py` | Rebuilds the existing native messages on three saved test-chromosome contexts, reads every both-strand object and verifies delivery/input identity. It performs no calibration or EM fit. |
| `verify.py` | Checks falsification, candidate and mutation receipts, source snapshots, census inputs and preserved production/native identities. |

The adapter reuses `composition()` and `held_dna()` from the existing single-strand
diagnostic adapter. For RNA it reads `CubeRows`, already produced by native final message
assembly and exposed in sweep diagnostics. There is no second propagation implementation.
An eventual production caller can use that existing captured table; it need not rerun
message assembly as the small census instrument does.

The native cube already intersects both received RNA levels with the object's own
lower-side bound, independently for each RNA strand. The adapter preserves that result,
including the presence bits. It does not apply the single-strand rule that excludes
RNA from a side supplying composition. Missing factors are neutral; unused array storage
is never interpreted as an observed profile.

Each returned RNA profile retains the cube's axis and density origin. DNA retains its
own axis and origin. The adapter does not use the cube's observed total to turn an absolute
density profile back into a fraction profile. Invalid or duplicate cube slot identifiers
are rejected. This is an internal diagnostic interface, not a new user-facing mode.

## Validation and its limits

Falsification ran against a wrapper around the old single-strand extractor before the
new adapter was written. All tests were also exercised by actual implementation defects.
The first mutation attempt found that the test's neighbor bounds masked the own bounds;
the corrected fixture makes each contribute. Its initial receipt is retained rather than
represented as successful mutation coverage.

The saved-input audit compares only present native entries between repeated assemblies;
unmarked native buffers are uninitialized. It also hashes the same allocated input buffers
before and after extraction, proving the adapter did not mutate them. These checks certify
delivery identity, not block/thread invariance, a joint statistical likelihood, final
calibration identity under a new consumer, or transcript accuracy.

`*.py.txt` files preserve executed scripts before formatting. `verification.json` links
them to current source by syntax-tree identity and records the native and main-tree hashes.
The adapter uses the already-frozen isolated native build from the known-source blur
checkpoint. No new C++ work was needed. The documentation boundary gate and scratch-code
lint are the relevant final checks; the production suite was not rerun for a scratch
adapter that has no production call site.

## Next bounded step

Implement the both-strand density readout using the already-declared total-RNA measure
and tilt mixture, including its two pure-strand hypotheses. At each candidate DNA density
and possible RNA amount/split, evaluate the existing RNA+ and RNA− profiles at the
corresponding RNA densities. Composition remains a ratio factor and DNA remains an
absolute density factor.

The own-observation both-strand reference already exists; delivered RNA profiles make the
tilt integrand non-polynomial, so its former exact polynomial quadrature cannot simply be
reused with a fixed number of points. Start with independent integration and the neutral,
single-strand, mirror, zero and narrow-profile limits. Establish numerical convergence
before porting or benchmarking. Keep the pure-strand mixture weights and total-RNA measure
unchanged; two independent RNA priors would be a new model.

This work can use fixed, explicitly supplied factors without deciding the pending
message-extent or population-weight proposals. It will not, by itself, resolve those
issues, inherited factor approximations or shared-observation dependence. Do not advance
this interface alone to end-to-end benchmarks or integrate the capture reader.

**Subsequent checkpoint:** the [both-strand numerical reference](RNA_BOTH_STRAND_EVIDENCE_CHECKPOINT.md)
now implements that fixed-factor readout. Its audit also makes explicit the existing
delivered-witness support rule: each present RNA profile excludes the opposite pure-strand
atom in the production solver. Passing finite interpolated values alone would lose that
restriction. The new diagnostic consumer preserves it; the production count path remains
unchanged.
