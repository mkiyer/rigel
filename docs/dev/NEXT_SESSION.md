# NEXT SESSION — finish the cleanup (owner, 2026-09-22)

Two prompts. The first is the next session; the second is saved until the cleanup is done.

---

## Prompt 1 — the remaining cleanup

> Read `CLAUDE.md`, then `docs/dev/CLEANUP.md` — the ledger of the cleanup the owner ordered on 2026-09-22: a
> production-ready tree WITH THE SHIPPED INFRASTRUCTURE and no new mechanism, as if 0.8.0 shipped what runs today.
> The mechanical part is done and committed (the ledger's "Done" table): the knobs, the instrument shelf, the
> prototypes, the worktree, the branches, the docs that named them, the changelog, the manual, a wheel. What is
> left is THE REVIEW ITSELF, in this order, each step a deletion or a simplification proven a numeric no-op —
> `rename_identity.py --check` against `~/Downloads/rigel_runs/arms/review_identity_<condition>.json` for each of
> the two ladder references and the LBX0190 library after every source change (the default reference path does
> not exist; read the log, never the exit code) — and the suite's count re-derived from `CLAUDE.md`'s table at
> every commit, never adjusted.
>
> 1. **The top-level Python, read end to end for simplification**: `cli.py` (1,561 lines), `pipeline.py` (1,074),
>    `estimator.py` (978), `index.py` (1,789), `scan_payload.py`, `buffer.py`, `second_pass.py`, then the
>    calibration package's largest modules (`splice_graph.py` 1,597, `calibrate.py` 966, `fl.py` 881). For each:
>    what is dead (coverage found the whole functions; the branches inside live ones are not measured), what is
>    duplicated, what is a comment describing a tree that no longer exists, what parameter no caller varies. One
>    commit per file or per kind of change. The C++ (`em_solver.cpp` 2,982, `bam_scanner.cpp` 3,283,
>    `solve_kernel.cpp` 1,422) is coverage-blind: count its functions by name against their callers.
> 2. **The tests** (55,633 lines, 1.2× the source): the vacuous gate clauses the ledger names in stage 3 (read
>    each gate against its toy; a clause the toy cannot reach is deleted or given a substrate that reaches it);
>    tests of deleted code; duplicate gates of one invariant; test files whose docstring describes a mechanism
>    that is gone. Every deletion moves the collected count by an amount derived from the table.
> 3. **The permanent docs** (7,777 lines): `ISSUES.md` (93 entries — every open entry re-read: is it still open,
>    still ranked, still naming an instrument that exists?), `DESIGN.md` (2,015 lines — rulings, not a diary:
>    a section that narrates a campaign is condensed to the ruling and its measured price), `TRAPS.md`,
>    `TESTING.md`, `SUCCESS.md`, `CLAUDE.md` (a rule earns its place by changing what the next session does).
>    The move rule throughout: one home, delete where it was in the same edit; nothing permanent cites `docs/dev/`.
> 4. **Release readiness, the last mile**: `preflight.py --full`; the suite; the wheel; then the two commands of
>    `docs/PUBLISHING.md` walked once WITHOUT publishing, and the `Development Status` classifier put to the owner.
>
> Rules: no new mechanism, no new instrument, no derivation; do not touch the pseudocount, the ruler or the
> gDNA component's length; the owner drives commits — propose them, one kind of change each, and commit on the
> owner's go. Stop when the ledger's stages are closed and record the state in `CLEANUP.md`; then this file
> points at prompt 2.

---

## Prompt 2 — the capture-contracted effective length (saved until the cleanup is done)

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
