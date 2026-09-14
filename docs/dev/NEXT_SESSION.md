# NEXT SESSION — start here (2026-09-14, after L3, L5, the landscape's location floor and the port's prerequisites)

The whole picture and the reasoning are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` (THE LANES WORKLIST table,
every row DONE, sits before §C; §F is the port); the θ thread's derivations are `docs/dev/THETA_QUADRATURE.md`.
This file is only how to begin. Everything is committed: `81739da4` (L3), `41b55006` (L5, with the handoff),
`84923136` (the landscape's location floor); the port's prerequisites are fresh on that tree (below).

## The agreed order of the sessions (owner, 2026-09-14)

1. **THIS ONE — code review and cleanup.** Dead, stale, legacy, unused, rotten: source, instruments, tests,
   docs, artefacts on disk. Every step a numeric no-op, proven on the identity references and the replay; the
   suite's count re-derived; `preflight --full` first and last. The concrete list is below.
2. **The test chromosome's new structures**, then both panels remeasured: the long single-exon transcript
   encompassing a multi-exon transcript on the opposite strand; the single-exon antisense gene wholly inside
   a sense exon (the tilt atom's accepted limit, so the panel carries its number); the shared exon of two
   spliced genes on opposite strands (the deep stress); a single-exon gene inside the opposite strand's
   intron (the majority AMBIG class); head-to-head genes with overlapping UTRs. Edit `test_chr.yaml`, rebuild
   and re-certify by `docs/TESTING.md` §0a (`panel.py status` names each stage); then
   `calibration_vs_oracle.py` and `policy_benchmark.py --by-class` on both panels become the new standing
   numbers (`DESIGN.md` §7).
3. **The ruler at zero gDNA** (`ISSUES: g00-shrinkage-upstream-repair`, the largest in-scope number on the
   metric page, the modal real case) and **the flux price's witness in whole-strand units**
   (`ISSUES: flux-price-witness-units`, the one accuracy item inside the port's unit). One mechanism at a
   time, each its own commit.
4. **The performance re-baseline** (two back-to-back profiler pairs at 8 threads, fresh `port_identity_*`,
   `sweep_replay.py capture`, one hour) **and THE PORT** (plan §F).

Deferred by the owner: a real-data census of the atom's accepted limit (no representative real data yet).
Outside the port's unit and after it: the landscape estimator's remaining items (a)–(c) in
`ISSUES: gdna-landscape-trains-on-false-positives`, the prior assembler's entries (`eb-shrinkage-magic-ess`,
`capture-blind-gdna-divisor`, `per-transcript-prior-lane`, `u-ruler-arm`, `prior-fidelity-vs-deliverable`)
and the EM's assignment error (the thermometer).

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, at `84923136` or later). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md`. Run
> `scripts/design/preflight.py --full` and the suite before touching anything; `CLAUDE.md`'s baseline line is
> the count to reproduce (3,458 passed / 3 xfail / 3,461 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code; a
> restructure is proven bit-identical; no number moves in this session — every step is checked with
> `rename_identity.py --check` against the three `port_identity_*` references (two ladder conditions and
> `--bam` LBX0190) and `sweep_replay.py replay --dir ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step6
> --call 0..3`; the suite's count re-derived, never adjusted; each item its own commit, committed on my go;
> never `ruff format scripts/`; converge and delete — no legacy, no compatibility shims, no speculative code.
>
> THIS SESSION IS CODE REVIEW AND CLEANUP, nothing else lands. Census before you cut: a fresh coverage run
> of the suite and every `--self-test` (the W11 method, `ISSUES: hygiene-ledger`), then every never-executed
> statement ruled dead, rotten-but-live, or a gap. The list to start from is in `NEXT_SESSION.md` (source,
> instruments, tests, docs, disk). Run `preflight.py --full` first and last: a self-test broke unnoticed
> for a day this week. End with the handoff for the test-chromosome session.

## The cleanup list — what this week's sessions saw

Source (`src/rigel`, 35,244 lines; calibration 14,344):
* `scripts/profiling/sweep_replay.py`'s unpickling shim for captures taken before 2026-09-11 (`rows_at`);
  every capture before `sweeps_MO_3021_step6` is stale, so the shim is legacy.
* Docstrings that assert measurements: grep `src/` for "measured", "measurably", "refuted" and check each
  claim against the record — `landscape._reliability`'s claim contradicted the landed floor until 2026-09-14,
  and `region_init.has_own_composition_evidence`'s described a deadband that no longer exists.
* Dead parameter plumbing after L3 (`od_g` and `n_gdna_obs` left the strand channel; check every caller and
  every `_Sweep` / `_Strand` / `InjectedCalibrationPriors` field is read), after the θ landings (any
  `sweep_n_tilt` / tilt-lattice residue) and after the one-solver landing.
* The coverage census: W11 found 84 % of `src/rigel` executed and removed 2,191 statements' worth of dead
  code; L1–L6 and the θ thread changed the surface since. Re-run it; the kept GAPS are listed in
  `ISSUES: hygiene-ledger`.
* `landscape_training_census.py` prints `nan%` for an empty evidence class (cosmetic); its selector copy
  must keep reproducing `_fit_gdna_hyperprior` bit for bit (its gate refuses otherwise).

Docs:
* `DESIGN.md` §6b.15 is 264 lines and thirteen rulings under one heading, a running log of the pre-port
  work; split it into named sub-rulings (§6b.15.1 …) with no content change, so a ruling can be cited.
* `CLAUDE.md`'s baseline paragraph is 27 lines of history; the keep-it-lean ruling wants the count and the
  last change, one line (the history is git).
* `TRAPS.md`: 164 traps in 960 lines; a pass for duplicates and for traps whose substrate is gone (the
  length channel, the relay, the deadband), keeping every lesson that still changes what a session does.
* `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`: the DONE rows could be compressed to a ledger; `docs/dev/
  THETA_QUADRATURE.md` §§12–14 are landed and recorded in the permanent docs (a note at each says where).

Disk (not the repo; the owner's call): `~/Downloads/rigel_runs/arms/` holds 24 identity references of which
only `port_identity_*` describe the current tree (the `23a431a9` set is in a subfolder); `perf/` holds six
capture sets of which only `sweeps_MO_3021_step6` replays on the current tree.

Tests: 3 xfails, all executable records (`two-sided-exon-row`, `antisense-prior-assembly-casualty`, the
encompassing flank under `the-lower-bound-noise-ratchet`); 21 golden scenarios, the four with gDNA identical
through every landing this week.

## The two pricing questions inside the port's unit — what the owner decides

Both live in `messages/transfer_rows.py` and `messages/lanes.py`, the port's first step, so a fix after the
port is written twice; the decision is fix-before or accept-and-record.

* **`ISSUES: flux-price-witness-units`.** A junction's certified flux is a lower bound on the exon's RNA
  level, hence a ceiling on its gDNA. Its price — how much the bound is blurred — is `hop_price`, which
  compares the junction's route rate (the strand's RNA density in whole-strand units) with the exon's own
  count on the read column, and that column holds only `(1 − κ')` of the strand's RNA plus half the gDNA.
  So the two disagree by `log(1 − κ')` even when both are exactly right, and every flux ceiling pays
  `log(1 − κ')²` nats² of spurious disagreement: negligible at κ = 0.99 (1e−4), 0.13 at ss.65-style
  protocols, 0.48 on unstranded data — an in-scope stratum. The recorded symptom: the golden
  `strand_ss65_multi_iso`'s nested exon reads 0.152 gDNA from a ceiling of 0.38 where the flux says ≤ 0.
  The naive fix (the column split's asymmetry as the witness) was REFUSED — it reads zero at an
  equal-abundance overlap exon; the owed form is the strand's RNA count in whole-strand units at
  single-strand exons and a bounded one at both-stranded exons. Recommendation: derive and A/B it in
  session 3, before the port carves the price.
* **`ISSUES: the-lower-bound-noise-ratchet`** (the encompassing flank's xfail). A level is a lower bound
  delivered at the sender's SAMPLED density; the hop's price blurs its edge but cannot move it, so a
  neighbour whose sample ran 1.6σ high pins the flank above the truth (0.596 against 0.530 on 404
  fragments; on the ladder 0.8 % of one row). The ruling stands as "later, with the enrichment witness".
  Recommendation: accept and record; the port carves the level rule as it is.

## Where things stand

* **L3** — the strand channel's gate is a protocol decision (`DESIGN.md` §6b.15, `EQUATIONS.md` §5.2b).
* **L5** — the tilt atom (`DESIGN.md` §6b.15, `EQUATIONS.md` §9f); its cost where no witness can exist is an
  accepted limit (`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`).
* **The landscape's location floor** (`DESIGN.md` §7.1 rule 4): a slot trains only if it holds a composition
  AND its solve locates it, `Var(log f_g) ≤ 1 nat²`; the ladder's four zero controls 500 → 282, 211 → 194,
  550 → 265, 231 → 172, every stratum unchanged or better on both panels.
* **The encompassing gate** (`test_encompassing_locus.py`): the atom's population passes on both donors; the
  flank is the ratchet's xfail.
* **The port's prerequisites on `84923136`'s content**: `sweeps_MO_3021_step6` (four calls, each replays
  BIT-IDENTICAL: 30.0 / 12.9 / 13.0 / 13.1 s); `port_identity_{gdna_g05_ss_0.50_nrna_mid_capture_off,
  gdna_g05_ss_0.99_nrna_mid_capture_on, LBX0190}.json`, each `--check`ed BIT-IDENTICAL; the deep-library
  baseline `perf/baseline_2026-09-14/pair{1,2}_{a,b}.json` (18,568,456 fragments, 8 threads: wall 499.4 /
  498.8 and 502.0 / 498.4 s, peak 11.2 / 11.0 and 10.8 / 11.1 GB, `calibrate` 404–408 s of which the four
  sweeps 394 s, scan 32, second pass 22, quant 25; every stage 0.96–1.02 within a pair). The landscape floor
  landed after the captures were taken and moved numbers, so the replay's bit-identity is against the L5
  tree's ψ; re-capture in session 4.

## The session scratchpad (persists across sessions; nothing in the tree cites it)

`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4ec3e3a5-268d-4828-a6e7-91f9cd89647e/scratchpad/`:
`l3/` (the deadband arms and their oracle runs, the walk, the census, the identity references before L3),
`l5/` (the atom's prototype, the smoke test, the census, both stresses, the audit, the port's logs),
`lp/` (`dump_training.py` — the training rows per fit against truth with four candidate rules scored,
`arms.py`, the oracle runs of the four arms, `golden_l5/`, `commit_three.sh`). The previous session's
`…/25f7f3df-…/scratchpad/w12/` holds the θ thread's material.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no
  θ lattice anywhere (the tilt count is derived); the tilt's hypothesis space is {pure +, pure −, mixed}.
* The tilt measure stays the arcsine and the marginal exact (both flattenings REFUSED on the ladder).
* The tilt is not a message; presence per strand is what the atom adds and the lanes witness; no presence
  witness is built (the owner, 2026-09-14: a limit of the information, accepted).
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the two other readings of the floor are refused.
* `drain`, `row` and `face` stay; the arcsine coordinate stays refused; the vertex atom is parked.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
