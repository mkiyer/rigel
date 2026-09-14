# NEXT SESSION — start here (2026-09-13, after W12 and W13; the port is next)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C, and §8 is the design of the
remaining steps. The θ quadrature's derivation is `docs/dev/THETA_QUADRATURE.md`. This file is only how to
begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, the commit after `d4934c5c`). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md`,
> `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` and `docs/dev/THETA_QUADRATURE.md`. The PRE-PORT WORKLIST stands
> at: W1–W11 DONE and committed; W12 DERIVED, its first prototype REFUTED, and its DESIGN DECISION owed — the
> session's work; W13 (the re-capture of the deep library's sweeps and fresh identity references for the port
> thread) not started. Run `scripts/design/preflight.py` and the suite before touching anything; `CLAUDE.md`'s
> baseline line is the count to reproduce (3,419 passed / 2 xfail / 3,421 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel, the suite, and
> timing on back-to-back profiler pairs; a restructure is proven bit-identical against fresh references
> (`hygiene_identity_*` are current as of `d4934c5c`); a removal on the default rows of both substrates
> (`arms/2026-09-13_preport/oracle_{ladder,test}_default.json`, still current); one mechanism at a time; no
> magic numbers; each step's plan presented and gated; each item its own commit, committed on my go; never
> `ruff format scripts/`; patch `calibrate` through `importlib.import_module`; any config value is a `--set`
> arm. Before migrating or extending an instrument, census it and retire it if its question is closed.
>
> W12 FIRST: the design decision, then DERIVE its error bound, PROTOTYPE outside `src/`, A/B on both panels
> with both zero controls AND a deep stress (a toy with 500k-fragment pure-RNA AMBIG exons — the ladder is
> converged at `sweep_n_tilt` 60 and the defect is latent there), then `src/`, the ruling to `DESIGN.md`
> §6b.15, the derivation to `EQUATIONS.md`, the issue closed with its numbers. Then W13. Then the port.

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## W12 — where it stands

Step 1 LANDED 2026-09-13 (uncommitted, the owner's go given): the θ nodes follow the strand term's peak
(`simplex_logodds._tilt_window`, `_read_row_at`; `_T_NATS` and `_TILT_NODES` = 24 derived), the ruling in
`DESIGN.md` §6b.15, the derivation `EQUATIONS.md` §9e, `ISSUES: theta-quadrature-at-zero-gdna` CLOSED with its
numbers, `ISSUES: strand-marginal-volume-factor` OPENED (priority next, after W12 by the owner's word). Seven
gates in `test_vertex_reference.py` (four marginal-vs-adaptive-quadrature cases, the derived count converged, a
linear row read exactly, κ = ½ the whole domain); five perturbations fired; two antisense goldens regenerated
(transcript counts 5e−5 relative, a tiny toy's `em_effective_length` 2.4 % on an 11-bp entry); `preflight.py
--full` 10/10; the metric on both panels and the shared-exon stress reproduce the prototype. The instruments
(`tilt_census.py`, `deep_stress.py`, `quadrature_check.py`, `theta_window.py`, `oracle_summary.py`) are in the
session scratchpad `w12/` and cited from `docs/dev/THETA_QUADRATURE.md`; whether any becomes a
`scripts/design/` instrument (the census is the candidate: the owner's "measure where the tilt matters") is the
owner's call, after a census of what it would replace.

**Step 2 LANDED the same day (uncommitted):** the RNA level lanes deliver a row's ingredients
(`simplex_logodds.CubeRow` — the two held profiles, the slot's total and RNA opportunity, the lanes' reference
densities) and ψ evaluates them at its own nodes (`CubeRow.at`); `CalibrationConfig.sweep_n_tilt`,
`_tilt_grid`, `_read_row_at`, the sweep's `(K, K_t)` shape check and the cache's row arrays are gone. No tilt
count exists anywhere in the tool: the only θ count is the derived `_TILT_NODES` = 24.

W13 DONE the same day (`sweeps_MO_3021_step5`, four calls bit-identical; `port_identity_*` frozen and checked). The agreed order before the port (owner, 2026-09-13): ① the tilt-measure thread — DONE, REFUSED both forms (`ISSUES: strand-marginal-volume-factor`, CLOSED / REFUSED with the ladder's numbers); ② dissect `ISSUES: capture-on-strand-pure-ambig-undercall` (no code unless the dissection names it); ③ re-measure the deep library end to end as the port's baseline, then re-freeze the references if anything moved; ④ THE PORT (plan §F). `ROADMAP.md` rank 2 carries the accuracy items that remain (the owner asked to be taught it; the
teaching is in the session's closing message and `EQUATIONS.md` §9e's last paragraph).

## Where the worklist stands

W1–W10 DONE (W10: `0f3ca3e3`…`bdc1089c`). W11 DONE 2026-09-13 as thirteen commits (`97947ede`…`d4934c5c`):
the suite and all ten self-tests under coverage, 2,191 never-executed statements reviewed one by one; removed,
each its own commit, byte-identical on the default rows and the identity references — seven dead members,
the unreachable no-calibration path, the RNA reach taper's unfed switch, the RNA arm's unfed fitted-prior
socket (ψ takes one fitted arm, `gdna_logprior`; `CompositionPriors` is gone), and three dead simulator
features (`sim/locus_sweep`, `sim/net_flow`, the synthetic mini-genome suite path). Net −4,080 lines in 34
files. Kept as COVERAGE GAPS, not dead code (`ISSUES: hygiene-ledger` lists them): the five CLI command bodies,
the silent policy through `calibrate`, the simulator's sharded writers and whole-genome path (live in panel
builds, silent in the suite), the zarr splice blacklist. `preflight.py --full` 10/10 after the deletions.

Next: the port (plan §F); the identity references are `port_identity_*` and the captures `sweeps_MO_3021_step5`.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee.
* The two lattice knobs are the model of a kept knob: dimensionless, a guarantee, a recorded ladder, not on
  the CLI. The EM's `gdna_em_llr_bias` stays (the owner's).
* `drain`, `row` and `face` stay (2026-09-13). The arcsine coordinate stays refused; the vertex atom is parked.
* The deconvolved arrays are `count_<population>_<axis>` (2026-09-13): counts on their axis, never masses.
* ψ has ONE fitted composition arm, `gdna_logprior`; a fitted RNA arm, if ever built, lands as `_rna_arm`'s
  second argument (2026-09-13).
* Four instruments retired rather than migrated (2026-09-13, the owner's ruling after a census); the record and
  the reasons are in `ISSUES: hygiene-ledger`. Second-tier candidates the owner has not ruled on:
  `native_parity_on_real_data`, `prior_units_check`, `accumulator_cost`, `verify_index_rebuild`.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
