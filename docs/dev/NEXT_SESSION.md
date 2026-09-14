# NEXT SESSION — start here (2026-09-14, after L3, L5, the landscape's location floor and the port's prerequisites)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
(THE LANES WORKLIST table — every row now DONE — sits before §C; §F is the port). The θ quadrature's
derivation, the tilt study and the encompassing-locus audit are `docs/dev/THETA_QUADRATURE.md` §§1–14. This
file is only how to begin.

## Three commits are waiting for the owner's go — in this order

The tree holds three landed items, kept apart: **L3 is STAGED** (the git index is exactly L3's change set, 136
files, most of them the 17 regenerated goldens); **L5** and **the landscape's location floor** are unstaged on
top, with this handoff. The scratchpad script makes the three commits in order — L3 from the index; L5 from
the working tree with its own golden set restored from `lp/golden_l5/` and the floor's four source files
held back; then the floor with its goldens — and the messages are beside it:

```bash
<scratchpad>/lp/commit_three.sh        # l3/COMMIT_MSG_L3.txt · l5/COMMIT_MSG_L5.txt · lp/COMMIT_MSG_LANDSCAPE.txt
```

The doc edits of the floor (DESIGN §7.1 rule 4, the two ISSUES entries, `CLAUDE.md`'s baseline) ride in the L5
commit, since the same files carry both; the source and the gate separate cleanly.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, with L3, L5 and the landscape's location floor committed). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md` and
> `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F. Run `scripts/design/preflight.py` and the suite before touching
> anything; `CLAUDE.md`'s baseline line is the count to reproduce (3,458 passed / 3 xfail / 3,461 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, BOTH zero controls), the panel, the suite, the
> encompassing-locus test and the shared-exon stress; a restructure is proven bit-identical; a change that
> moves numbers is judged with its magnitudes read first; one mechanism at a time; no magic numbers — every
> constant derived; a falsification test first, verified failing, then break the fix and watch each gate fire;
> each item its own commit, committed on my go; never `ruff format scripts/`; patch `calibrate` through
> `importlib.import_module`; any config value is a `--set` arm.
>
> THE PORT (plan §F): the unit is `sweep._solve_block`. Order: (i) the passes and `transfer_rows`, (ii)
> `prepare`'s builders, (iii) ψ (SIMD exp/log, the AMBIG cube — now `K × (K_t + 2)` with the tilt atom's two
> columns), (iv) threads over blocks. Every step behind the tolerance gate (`sweep_replay.py replay
> --tolerance` on `sweeps_MO_3021_step6`), the fresh `port_identity_*` references and the suite; the λ lattice
> stays a parameter (`sweep_logodds_step`). Read a timing only against the baseline pairs of 2026-09-14
> (`~/Downloads/rigel_runs/perf/baseline_2026-09-14/`, `profiler.py --compare`).

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## Where things stand

* **L3 — the strand channel's gate is a protocol decision** (`DESIGN.md` §6b.15, `EQUATIONS.md` §5.2b,
  `ISSUES: deadband-gates-a-gdna-free-library` CLOSED). `region_init.strand_discriminability(kappa, n_rna_obs)`
  is `4(κ̂−½)²` iff the spliced 2×2's Bayes factor favours a free κ over κ = ½ exactly; `n_gdna_obs` is gone
  everywhere. All four ladder `g00` rows had `N_gdna = 0`, so the stranded `g00` rows — and every `g00`-donor toy
  of the θ thread — had run with the channel dead. Ladder identical but unstranded OFF −0.22 %; the stranded
  zero controls 405 → 497 / 194 → 224 are the LANDSCAPE's vertex bias at the refit rung (filed as
  `gdna-landscape-trains-on-false-positives` (d)), not the gate's. The belief-read RNA level REFUSED as the relay.
* **L5 — the tilt atom** (`DESIGN.md` §6b.15, `EQUATIONS.md` §9f, `ISSUES: capture-on-strand-pure-ambig-undercall`
  CLOSED). Two atoms at `τ = ±1` in ψ's cube beside the continuum (`−log π` on its weights); a held level on a
  strand rules the other strand's atom out. Ladder stranded ON 470,862 → 427,046 (−9.3 %), every stratum
  better; stranded zero controls 497 → 550 / 224 → 231. ⚠ Its cost where a strand has NO possible witness is
  filed with its numbers as `ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` — the golden
  `antisense_contained` (a single-exon antisense gene inside a sense exon) reads its antisense transcript
  81 → 0 and 177.6 false gDNA of 1,000; the mono shared-exon stress doubles its false gDNA at 2–20 % minor.
  The owner's stance (2026-09-14): a limit of the information, accepted — no presence witness; the landscape
  prior decides on a real library, and its vertex bias (`gdna-landscape-trains-on-false-positives` (d)) is
  the lever. `ISSUES.md` was pruned the same day (done-records out, seven entries retired into the record).
* **The landscape's location floor** (`DESIGN.md` §7.1 rule 4; the owner's direction 2026-09-14: which nodes
  may teach the circular half). A slot trains only if it holds a composition AND its solve locates it,
  `Var(log f_g) ≤ 1 nat²` (`landscape._LOCATED_VAR`, the E-step's one-fragment floor through
  `Var(log c) = 1/c`); anchors regardless. The ladder's four zero controls 500 → 282 / 211 → 194 / 550 → 265
  / 231 → 172, every stratum unchanged or better on both panels (deferred −0.87 %); on `g00 ss.99 OFF` the
  training set from the second refit is the anchors alone. Two other readings of the floor refused with
  their numbers (the ruling). The tiny goldens move the other way (`antisense_contained` 177.6 → 200.8): a
  prior from four anchors pushes less.
* **The encompassing gate** (`test_encompassing_locus.py`): the exon∩exon slots and the region between TA+'s
  exons solve within 0.05 in every regime on both donors (un-xfailed); the one remaining miss is TB−'s
  shallow flank under the intergenic gDNA EDGE level, its own xfail under `the-lower-bound-noise-ratchet`.
* **The port's prerequisites on the landed tree** (this session, after L5): `sweeps_MO_3021_step6` captured
  (four calls) and replayed; `port_identity_{gdna_g05_ss_0.50_nrna_mid_capture_off,
  gdna_g05_ss_0.99_nrna_mid_capture_on, LBX0190}.json` re-frozen and `--check`ed (the `23a431a9` set moved
  to `arms/port_identity_23a431a9/`); the deep-library baseline (`mctp_vcap_rna20m_dna05m`, 8 threads, two
  back-to-back pairs) in `perf/baseline_2026-09-14/pair{1,2}_{a,b}.json`. The verdicts and timings are in the
  section below.
* `vertex_ceiling.py --self-test` had been broken by L1 (a bare stub context; `prepare` now reads the whole
  context) — fixed with a real empty block; `preflight --full` reads 10/10.

## The port's prerequisites — the record (2026-09-14, the L5 tree)

All taken on the L5 tree (`cd470d0e` + L3 + L5 uncommitted), the machine otherwise idle.

* **`sweep_replay.py capture`** → `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step6` (the MO_3021 library, 8
  threads): four calls; each replays **BIT-IDENTICAL** on this tree — 30.0 / 12.9 / 13.0 / 13.1 s (the step5
  capture on `23a431a9` read 30.0 / 12.7 / 12.7 / 12.7). `step5` and earlier describe trees before L2, L3 and L5.
* **`rename_identity.py --freeze`** → `~/Downloads/rigel_runs/arms/port_identity_{gdna_g05_ss_0.50_nrna_mid_capture_off,
  gdna_g05_ss_0.99_nrna_mid_capture_on, LBX0190}.json`, each `--check`ed **BIT-IDENTICAL** on this tree. The
  `23a431a9` set is in `arms/port_identity_23a431a9/` (stale: L2, L3 and L5 moved numbers).
* **The deep-library baseline** (`mctp_vcap_rna20m_dna05m`, 18,568,456 fragments, 8 threads, two back-to-back
  pairs, `profiler.py`) → `perf/baseline_2026-09-14/pair{1,2}_{a,b}.json`: wall 499.4 / 498.8 and 502.0 /
  498.4 s; peak RSS 11.2 / 11.0 and 10.8 / 11.1 GB; `calibrate` 404.3 / 405.1 and 407.6 / 404.8 s, the four
  sweeps 394 s of it (the ψ grid solves inside them ~102 s at 1,704 calls), the init ψ 4.5, the landscape fits
  4.3; scan 32, second pass 22, quant 25. Within each pair every stage reads 0.96–1.02 (`profiler.py --compare`),
  so the port's timings are read against these four reports, never against a number from another sitting.

## The session scratchpad (persists across sessions; nothing in the tree cites it)

`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4ec3e3a5-268d-4828-a6e7-91f9cd89647e/scratchpad/`:
`l3/` — `arms.py` (`run <arm> [instrument.py] -- args`: patches both bindings of `strand_discriminability`
and runs any instrument in-process), `run_oracle.sh` (one process per condition, merged), `strand_scalars.py`
(κ̂, N, the Bayes factor per condition), `deep_l3.py` / `audit_l3.py` (the θ thread's harnesses with L3 arms),
`perturb.sh`, the oracle jsons (`ladder_ref/drop/bf`, `test_ref/bf`), the walk, the census, the identity
references frozen before L3 (`identity/l3_pre_*`, both BIT-IDENTICAL after), the golden diff. `l5/` —
`tilt_atom.py` (the prototype), `smoke.py` (the strand-pure smoke test), `tilt_census.py`, `run_oracle.sh`,
the oracle jsons (`ladder_l3/atom_w/landed`, `test_*`), `census/` (the bands, L3 tree vs landed), `deep/`
(both stresses, before and after), `audit/`, `dissect_encompass.py` (the encompassing regimes slot by slot,
atoms on and off), `land_atom.py` / `adapt_gates.py` / `perturb.sh`, `port/` (the prerequisites' logs). `lp/` —
`dump_training.py` (the training rows per fit against truth, four candidate rules scored; `dump/*.json` on
seven conditions), `arms.py` (the admission arms through the block solve), the oracle jsons for the four
arms, `census/`, `walk/`, `golden_l5/` (the L5-state goldens the commit script restores), `commit_three.sh`. The
previous session's scratchpad (`…/25f7f3df-…/scratchpad/w12/`) still holds the θ thread's material.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no
  θ lattice anywhere (the tilt count is derived); the tilt's hypothesis space is {pure +, pure −, mixed}.
* The tilt measure stays the arcsine and the marginal exact (both flattenings REFUSED on the ladder).
* The tilt is not a message; presence per strand is what the atom adds and the lanes witness; the structural
  witness (the exon bits) adds nothing and is not written.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* `drain`, `row` and `face` stay; the arcsine coordinate stays refused; the vertex atom is parked.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
