# Weak-strand OFF: evidence-input checkpoint

*2026-10-08. Diagnostic extension of the frozen count candidate. No production change,
new model, commit or push.*

The weak-strand OFF regression needs more than removing the landscape's admission
cutoff: every eligible region already enters its training set. The count error is
present before the first landscape fit. Feeding those uncertain estimates back as
observed DNA also makes the population input too certain. The subsequent frozen-input
replay improves introns in both stress cases but trades errors between other classes.
It is not a sufficient population partner or a release candidate.

The settled measurements have one home: **WEAK-STRAND INPUT TRACE** under
`ISSUES: the-gdna-prior-enters-psi-twice`, followed by **FROZEN STRESS-INPUT COMPARISON**.
The earlier full count/control comparison is
**COUNT-CANDIDATE READINESS** in the same entry.

## Exactly what was implemented

The scripts are diagnostic wrappers in
`.cache/rigel_runs/2026-10-08_stress_evidence/`:

- `freeze.py` reuses the existing two stress fixtures and their factory-kept, removed
  and original-landscape controls. It captures each refit's objects, opportunities,
  beliefs, training selection, retained reliability weights and fitted landscape.
  For factory-removed runs it also extracts the existing composition, DNA-level and
  RNA-level factors through the already-tested adapter. The native assembly is
  reconstructed and checked for exact equality.
- `analyze.py` reads those saved inputs, joins their region/boundary counts to the
  certified read-name partition and compares three inputs at the same physical
  densities: raw inferred-count Poisson, own observed-strand evidence, and that
  evidence with existing delivered factors. It calls the previously tested native
  evaluator; it implements no new likelihood and fits no population.
- `fit_stress.py` reconstructs the old rendered curves and their mean, then reuses the
  certified population optimizer for three successive estimator/input contrasts.
  Pure-DNA/empty and mixed-object replacements are separate. It preserves the objects,
  weights, density grid, prior metadata and existing uniform pseudo-region strength.
- `replay_stress.py` regenerates each small fixture once and validates its read-name
  partition. A pipeline control and direct calibration control reproduce the earlier
  candidate exactly. The three contrasts reuse those calibration inputs and replace
  only the final landscape. Every refit's incoming beliefs and opportunities are
  checked against the frozen arrays; no candidate EM is run.

The density comparison is against the measured intergenic background, not a fitted
capture reference or exact per-object truth. The own-observation integral retains
the declared RNA reference and measured strand probability. The point-count comparison
precedes the old fitter's smoothing and previous-prior weighting. This experiment
does not certify the inherited message approximations or remove shared observations.

### Saved-factor attribution after the bootstrap discussion

The follow-up scripts are in `.cache/rigel_runs/2026-10-08_evidence_attribution/`:

- `attribute.py` evaluates the existing own, composition, DNA-level and RNA-level
  inputs at the same two saved densities. An independent concave optimization checks
  the older own-read profiling formula; this is a nuisance-treatment diagnostic, not
  a switch from the declared RNA integral to a profile.
- `trace.py` reconstructs the two small fixtures once each, verifies the archived
  controls, and stops before EM. Six saved contexts permit one EDGE producer at a time
  to be withdrawn with the existing face roles, parameters and all noncomposition
  factors held fixed. The model used by the candidate never changes.
- `verify.py` checks the receipts, input identities, unchanged production source/tests
  and installed native, the expected source reach and the full-control reconstruction.

Settled findings belong to **FACTOR ATTRIBUTION** in the same ISSUES entry. The extreme
case's pressure comes from existing local gene-edge messages. Their wholesale removal
has already been measured and is not a useful general replacement. The moderate case
does not share that complete explanation. This is progress in attribution, not a new
accuracy improvement or evidence that the population partner is ready.

## Verification and limits

All six final count controls exactly reproduce the preceding source-blur checkpoint.
The retained arrays also include the unchanged final populations and variances. Scan
threads are one, EM is fractional with seed zero, and calibration/EM thread counts
retain their defaults. The two 1,000-fragment fixtures use the existing simulator
inputs; no probe placement, length distribution or chemistry was optimized.

The freeze contains twelve refit snapshots. Its whole-array comparisons and typed
assembly checks pass. The evidence evaluation took seconds and did not run another
count solve. These checks validate a diagnostic instrument; there is no new inference
implementation to claim passing falsification or mutation gates for.

The final-refit replay now tests one consequence of profile training. It improves the
intron counts, while worsening extreme-case exons and moderate-case boundaries against
the frozen candidate. That class tradeoff prevents an aggregate improvement from being
treated as a release result. The six numerical fits pass their existing objective-gap
certificate. The complete table is in the ISSUES entry above.

This does not establish that local message bias has only one cause. The existing
delivered factor still favors excess DNA in the worst intron. The separate density-range
defect remains unresolved and its representation proposal awaits the owner.

## Disposition and remaining work

The planned frozen estimator/input sequence is complete. Do not rerun it as an untested
next step or advance this partial population model to large panels. The count candidate
remains frozen. The independent density-range proposal still awaits an owner decision.
The subsequently approved equal-weight comparison is complete and also trades class
errors; see [RNA_EVIDENCE_ADMISSION_CHECKPOINT.md](RNA_EVIDENCE_ADMISSION_CHECKPOINT.md).

Continue certification of the existing evidence factors and the prior/population
assembly. Both-strand evidence and the self-consistent population partner remain
unfinished. The current comparison does not license applying the DNA prior once without
its reviewed partner. Preserve both stress fixtures in future validation alongside the
zero control and the broader panels. Do not reinstate the intron constraint, tune a
strand threshold or introduce class weights to fit them. Their nascent share is a stress
reading, not the release design target.

Receipts: `freeze.json`, `diagnostic_stress_counts.json`, the per-refit/per-control
NPZ files, `analysis.json`, `*_fits.json`, `*_replay.json` and their logs. The native binary is the isolated
`2026-10-08_level_support/site` build; the installed library is unchanged.
