# Existing-message locality: a small repair and explicit remaining failures

*2026-10-08 experiment; integration status updated 2026-10-09. No commit or push.*

An unrelated exon should not decide how a fixed object's neighbors are interpreted.
The measured strand protocol is shared technical information. RNA abundance elsewhere
should not become local biological evidence through an incidental implementation switch.

The audit found such a switch and separately found failures in the numerical RNA axis.
The switch and zero-coordinate source loss now have tested Python repairs. A separate native
repair corrects endpoint padding when blurring a known flux source. The broader finite-support
failure remains unresolved.
The settled findings live in **MESSAGE-LOCALITY AUDIT** under
`ISSUES: the-gdna-prior-enters-psi-twice`. The coordinate algebra lives in
**Changing a coordinate preserves physical support** in `../EQUATIONS.md`.

## Exactly what changed

The implementation is outside `src/`, in
`.cache/rigel_runs/2026-10-08_message_locality/local_witness.py`.
It subclasses the existing policy for isolated tests and provides the same method as a
process-local hook for calibration A/Bs. It replaces one Boolean reduction:

```python
# Existing availability rule
model_is_available and protocol_is_live and any_counted_single_strand_exon

# Prototype availability rule
model_is_available and protocol_is_live
```

The already-measured protocol decides whether a strand split can be interpreted.
The local counts still determine that split and its uncertainty. No new threshold,
statistical population, neighborhood, prior, message formula or native code is added.
The wrapper computes the original library facts and replaces this one field, keeping
both density-coordinate origins exactly as before.

This is a correction to the stated protocol rule. It does not certify the inherited
hop-width formula as an exact likelihood. In particular, with RNA on both strands,
the column difference reflects their contrast rather than either individual RNA amount.
That existing approximation is unchanged.

## Reproducible checks

| Artifact in the scratch directory | What it establishes |
|---|---|
| `audit.py`, `audit.json`, `curves.npz` | Rebuild native messages on fixed local observations; vary remote exon counts, coordinate support and grid spacing separately. |
| `test_witness.py`, `witness_before.xml`, `witness_candidate.xml` | Four of seven checks fail before the repair; all seven pass afterward. Both strand orientations, absent/dead protocol, unchanged origins and final local readouts are covered. |
| `mutate_witness.py`, `witness_mutations.json` | Reintroducing the remote gate, dropping either technical guard or changing the origin causes the relevant tests to fail. Every witness gate is exercised. |
| `actual_faces.py`, `actual_faces.json`, saved local face | Repeat the remote-exon intervention on annotation-derived neighborhoods and inspect the RNA bound delivered to a both-strand recipient. |
| `panel_census.py`, `panel_census.json`, `witness_panel.json` | Seven current-condition calibration A/Bs; every library input and final DNA-count array remains identical. No new transcript/gene error claim follows. |
| `test_coordinate_support.py`, `coordinate_unrepaired.xml` | The original three failing specifications for missing sources and lost grid support. The subsequent fallback below repairs two; one remains failing. |
| `verification.json` | Source/binary identities, complete test outcomes, mutation coverage and paired input/count checks. |

The native package is the already-tested separate-channel evaluator. No rebuild or
modification of the installed binary was required. The main count patch, test references
and goldens remain the owner's existing files. Earlier full-panel and real-library
measurements still describe the frozen count candidate; the new population and reader
are not production-ready.

## Separate zero-coordinate repair

This implementation is in
`.cache/rigel_runs/2026-10-08_rna_coordinate/positive_coordinate.py`. It calls the original
library reduction and returns its result unchanged whenever the RNA reference is positive.
Otherwise, it constructs a numerical scale from the available count/exposure pairs:

- Counted RNA-admitting objects contribute their total unspliced count and RNA opportunity.
- Certified route ends contribute their positive count and `count / route_rate` exposure.

The reference is the sum of those numerators divided by the sum of their exposures.
No pair is added when its count or opportunity/rate is absent. With no available pair,
the original result is retained. The coordinate's units transform with RNA opportunities;
DNA opportunities do not enter it. The original native source builders decide which
claims exist. This computes no RNA abundance estimate, prior or biological rate floor.

The fallback is tested alone, keeping the previous witness bit and every other library
field fixed. `test_coordinate.py` first fails 12/17; the repair passes 17/17. Six actual
defects exercise all gates. Both control and repaired calls pass the existing 17 native
RNA-lane tests. The seven archived panel contexts give identical library inputs and
unchanged arrays, so there is no new transcript/gene accuracy measurement.

The receipts are `before.xml`, `candidate.xml`, `mutations.json`, `existing_control.xml`,
`existing_candidate.xml` and `inputs.json`. `existing_locality.xml` re-runs the original
three locality specifications through the fallback: the two source-loss cases pass,
and the finite-support case still fails. The standalone native package is unchanged.

## Separate known-source blur repair

This implementation is in the isolated native worktree, with the exact header snapshot,
package and receipts in `.cache/rigel_runs/2026-10-08_level_support/`. Only
`transfer_rows.h` changes relative to the separate-channel prototype:

- Share the existing Gaussian kernel-radius calculation between the blur and its caller.
- Add reusable temporary axis/source/output buffers to the existing scratch object.
- In `flux_level`, evaluate the known Poisson source over the blur's footprint beyond
  the output table; blur there, retain the original cells and apply the same lower-side rule.
- Preserve the unblurred branch exactly and check buffer-size arithmetic before allocation.

The inference grids, source admission, hop widths, opportunity maps, propagation schedule,
message intersections and priors are unchanged. The extra source evaluations use the same
spacing; they are not a second inference lattice. Generic transported rows still use their
existing boundary treatment. The derivation is **A known source supplies the convolution
footprint** in `../EQUATIONS.md`; the settled measurement is in **MESSAGE-LOCALITY AUDIT**
under `ISSUES: the-gdna-prior-enters-psi-twice`.

| Receipt | What was checked |
|---|---|
| `audit.py`, `audit.json` | Independent log-sum-exp evaluation of the existing finite Gaussian convolution; retaining the source mode alone does not make endpoint padding valid. |
| `test_source_blur.py`, `before.xml`, `candidate_repaired.xml` | Control fails two of seven; repaired native passes all seven. Low/high rates, coordinate translation and exact zero-width behavior are checked. |
| `mutations.json` | Four compiled defects exercise all seven gates: restore padding, omit blur, discard the coordinate origin, and add width at zero. |
| `existing_repaired.xml` | All 17 existing RNA-lane checks pass. |
| `control_screen.json`, `candidate_screen.json` | Fresh seven-condition, calibration-only A/B; the control exactly reproduces frozen counts. Both length gaps and the zero control are included. |
| `control_stress_counts.json`, `candidate_stress_counts.json` | Both count-readiness stress cases, with factory kept/removed/frozen-landscape controls; all six paired region/boundary count arrays are bit-identical. |
| `existing_locality.xml` | The two repaired zero-coordinate specifications still pass. The remote-origin support-loss specification still fails. |
| `verification.json` | Exact source/binary identities, gate outcomes, perturbation coverage and archived-input comparisons. |

The fresh calibration screen reads **DNA-count error as percent of observed incidences**,
not transcript/gene error. Each cell below is frozen count candidate → source-blur repair.
No conditions are pooled.

| Condition | Regions (%) | Boundaries (%) |
|---|---:|---:|
| RNA-long, g50 ss0.50 ON, deferred | 41.675819 → 41.675809 | 32.913219 → 32.913214 |
| Ladder, g98 ss0.50 OFF | 2.464644 → 2.464645 | 9.123508 → 9.123508 |
| RNA-short, g50 ss0.99 OFF | 0.924558 → 0.924558 | 2.948237 → 2.948237 |
| RNA-long, g50 ss0.99 ON | 1.465722 → 1.465722 | 3.348012 → 3.348012 |
| Test, g00 ss0.99 ON | 0.003518 → 0.003518 | 0.023910 → 0.023910 |
| Test, g50 ss0.99 ON | 1.337062 → 1.337063 | 1.133639 → 1.133639 |
| Test, g50 ss0.50 ON, deferred | 3.688359 → 3.688369 | 2.542818 → 2.542821 |

Introns and exons are also scored separately in the receipts. The paired calibration
loops took about 20 seconds each; this is execution cost, not a measured speed-up.
The stress pair pins the scan and fractional EM, uses seed zero and default calibration/EM
thread counts. The known weak-strand OFF excess remains; this numerical repair does not
recover it. Existing full-panel transcript/gene tables still describe the frozen count
candidate. No new end-to-end or real-library claim is made.

Two implementation issues were resolved before these checks. The first signed-size
implementation crashed in the optimized native build and a standalone C++ reproduction;
the checked unsigned-size version passes normal release optimization and the address/
undefined-behavior sanitizer check. The small reproductions remain in `compiler_audit/`.
The first mutation campaign also caught a missed rebuild: the final two header edits
fell within Make's timestamp resolution. The corrected harness removes only its generated
translation-unit object before rebuilding and records every distinct native hash
(`TRAPS: a-rebuild-can-silently-no-op`). The successful receipts are from that corrected run.

## What remains and how to proceed

**A finite coordinate grid must cover the evidence it represents.** Refining points within
the wrong physical interval does not help. The widened-grid control diagnoses this failure;
it is not a recommendation to change the configured window to its diagnostic value.
Before changing the native range, write a numerical-support specification from the actual source
likelihoods and their limiting behavior. Test origin changes, large density ratios,
zero-component modes and both length-gap directions. Use convergence and the existing
numerical error budget; do not add a biological cutoff, an expression ceiling or a new
spatial-sharing rule. Keep this contrast separate from the fallback correction.

That source-support diagnostic has now been written and tested. It passes its 38 source
specifications but fails two subsequent repeated-propagation specifications, so it has
not advanced to calibration or release runs. The settled rejection and convergence
measurements are **SOURCE-RANGE LIMIT** under `ISSUES: the-gdna-prior-enters-psi-twice`.
[RNA_MESSAGE_SUPPORT_DESIGN.md](RNA_MESSAGE_SUPPORT_DESIGN.md) proposes the next numerical
representation checkpoint. It awaits the owner because independent density-message
extents would amend the current requirement that profiles use the solve grid. The
composition/count lattice and existing spacing would stay; implementation has not begun.

**Do not confuse robustness defects with the measured regression's cause.** The witness
bit already agrees with the protocol in all seven checked conditions. Fixing it cannot
recover their deferred RNA-long loss. The weak-strand OFF count stress case and the
landscape/zero-control limitations remain on the release path. The separate equal-weight
population experiment was subsequently approved and completed; its outcome is in
[the admission checkpoint](RNA_EVIDENCE_ADMISSION_CHECKPOINT.md). No population change
is part of this locality work.

The witness-availability repair was integrated separately on 2026-10-09, after the count
foundation. It has portable gates and fresh calibration comparisons; its settled receipts
are **LOCAL-WITNESS INTEGRATION** under `ISSUES: the-gdna-prior-enters-psi-twice`. The
zero-coordinate repair was then integrated with expanded source/unit gates; its settled
receipts are **ZERO-COORDINATE INTEGRATION** in that issue. The known-source blur repair
follows as **SOURCE-FOOTPRINT INTEGRATION**, with its additional spacing-rounding gate,
downstream transcript tradeoff and real-library comparison. Complete locality
certification, the profile population partner and capture-reader integration remain outstanding.
