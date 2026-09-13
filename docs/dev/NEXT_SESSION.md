# NEXT SESSION — start here (2026-09-13, after W10 and W11; W12 is the next session's work, then W13 and the port)

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

## W12 — where it stands, and the decision owed

Derived (`docs/dev/THETA_QUADRATURE.md`): at fixed λ the frozen-variance strand term is an exact Gaussian in
the tilt τ, with centre `d/((1−f_g)(κ−½))` and width `σ_p/((1−f_g)|κ−½|)`; the θ-marginal is that Gaussian
against the arcsine weight; the uniform θ lattice is the trapezoid rule in `φ = θ + π/2` with both endpoint
nodes at weight 1 instead of ½; and beyond the strand-pure boundary the peak narrows as `√f_g`, so a fixed
lattice inflates the marginal more at higher gDNA — the gDNA bias is the λ-dependence of the resolution
error. Both recorded numbers are reproduced. Prototype A (the endpoint weights) is REFUTED: ≤ 0.7 fragments on
any `g00` row of either panel, the `K_t` 30 failure untouched (9,821 → 9,821 against 194 at 60).

The decision: **B**, Gauss–Hermite per λ centred on the Gaussian with the delivered cube rows and the arcsine
factor interpolated from the lattice (recommended: it keeps the delivered rows on their lattice and needs only
interpolation; the node count falls from 60 to ~8–16, the cube's cost lever), or **the analytic marginal** as
a precomputed table in `(τ̂, log σ_τ)` (removes the θ axis outright, but the delivered cube rows then need
the sifting approximation, exact only as `σ_τ → 0`). Either must be judged on a deep stress as well as the
panels, since the ladder is converged at 60.

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

Next: W12 (above), then W13, then the port (plan §F).

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
