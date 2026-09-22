# NEXT SESSION — the capture-contracted effective length (owner, 2026-09-22)

The cleanup is done (`CLEANUP.md` records it). The next session's prompt:

> Read `CLAUDE.md`, then the sandbox's `TWO_PROBLEMS.md`. Problem 2: **correctly estimate the capture-contracted
> effective length for each transcript** from the calibration result. Use the six words of that note and no
> others; if the derivation needs more than a page, it is off track.
>
> Start from what ships and what is measured about it, all in `docs/ISSUES.md` and `docs/EQUATIONS.md` §11:
> transcripts and synthetic spans take the ruler (a sum over their own bases of capture efficiency times the
> fragment-end taper, `capture_eff_length.py`, `priors.assemble_priors`); the gDNA component takes a different
> rule (a sum over region and boundary objects, converted by the accumulator's `q`); against the simulator's own
> yield the two rules disagree by 14–16 % on capture, and this error cancels the pseudocount prior's bias, so
> neither can be repaired alone (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`,
> `ISSUES: the-pseudocount-prior-is-biased-toward-gdna`); within a probed class spliced and unspliced templates
> disagree by 3–5 % under any per-region efficiency, the junction problem
> (`ISSUES: ruler-witness-geometry-on-transcript-panels`). The instrument that judges it exists and runs with no
> EM: `ruler_vs_truth.py --scale` reads L / Y per hypothesis class against the simulator's yield, and the pass is
> that EVERY component class reads the same ratio to its truth within a stated tolerance — overall accuracy is not
> the metric; `ruler_vs_truth.py` (no `--scale`) reads the per-transcript length against the sampler's own
> capture-aware length, per probed class.
>
> Do, in order: (1) write the theory in one page — a base's capture efficiency is the probability a fragment
> covering it is captured, read from gDNA as its piece's density against the fully captured level; a
> component's capture-contracted length is the sum over its own bases of efficiency times taper; a junction is
> a coordinate, its fragments counted through the exon bases they cover — and state on that page which of the
> shipped rules (transcript, span, gDNA component) obeys it and which does not; (2) prototype ONE rule for every
> component outside `src/` (a `--module` ruler in `ruler_vs_truth.py`, or a worktree if the kernel is touched),
> A/B'd against the shipped rule on the same conditions with the identity references beside it; (3) read the
> per-class table on `g50 ss.99 ON`, `g05 ss.99 ON` and one capture-OFF row, the junction residual sized and
> reported with its number, not hidden; (4) only then the EM: `quant_accuracy.py` per stratum under fractional
> assignment, above the reseed floor, the pools per row before the transcript table, the zero rows beside — and
> remember that removing the length error alone makes the tool worse on capture because it was cancelling the
> pseudocount's bias, so the two are judged together, never one at a time. Stop and discuss before anything
> enters `src/`; the owner drives commits.
