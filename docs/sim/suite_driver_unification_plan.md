# Plan — unify the two parallel suite drivers (`suite.py` / `whole_genome.run_simulation`)

**Status:** plan for review (2026-06-07).

## The problem (why this exists)

Two modules each carry a full **condition-grid orchestration loop** — the cartesian product over
`nrna × gdna × strand-overdispersion × strand-specificity × capture`, the per-condition
simulate-and-write, and the manifest/truth emission:

- `whole_genome.run_simulation` — grid build `2055–2094`, loop `2096–2250`.
- `suite.main` — grid build `1155–1168`, loop `1179–1313`.

They are **near-duplicates**. Adding the gDNA strand-overdispersion axis required wiring the
*same* logic into both (suite `888–907, 1165, 1189–1191, 1227–1228, 1277, 1305`; whole_genome
`2065–2081, 2127–2129, 2172–2173, 2215, 2243`). Any future axis, manifest field, or seeding fix
must be made twice — the duplication this cleanup removes.

## What is NOT duplicated (stays put)

- **`suite.py` only:** synthetic-genome generation (`925–973`) and capture-**probe design**
  (`375–695` + `SuiteCaptureSpec`/`_suite_capture_specs` `230–372`). A suite-specific front-end.
- **`whole_genome.py` only:** the `WholeGenomeSimulator` engine (`889–1896`), the config
  dataclasses (`230–320`), abundance-file modes (`778–823`), and parallel sharding (`1902–1936`).

Key fact: **`suite.py` already imports `whole_genome`'s config dataclasses + simulator** (it builds
a `WholeGenomeSimConfig` at `1094`). The dependency is one-directional (suite → whole_genome), so
consolidation does not create a cycle.

## Target architecture

Make **`run_simulation` the single orchestrator.** Extract the grid+loop+manifest/truth body into
one reusable function in `whole_genome.py` (or a small `rigel/sim/orchestrator.py`), parameterized
by the two things that currently differ between the drivers:

1. **Abundance/fragment allocation** — whole_genome supports a `file` mode (`n_mrna = None`
   branch); suite uses a simple additive allocation (`n_mrna = n_rna`, `n_nrna = round(n_mrna·r)`,
   `n_gdna = round(rate·n_mrna)`). Reconcile by having the orchestrator take an already-built
   `WholeGenomeSimConfig` whose `AbundanceConfig`/`NRNAConfig` already encode both modes (suite
   just constructs the additive config it already builds today).
2. **Seeding** — suite uses `capture_paired_condition_seed` (capture variants of one base
   condition share a seed for paired comparison); confirm whole_genome's seeding and make the
   orchestrator use the **paired** scheme for both (it is a strict superset — non-capture suites
   are unaffected because there is one capture scenario per base condition).

`suite.main` then becomes a thin **front-end**: (1) generate genome, (2) design probes, (3) build
the `WholeGenomeSimConfig` (incl. the designed capture configs), (4) call `run_simulation`. The
condition loop + axis wiring + manifest/truth live **once**.

> **Open design question for review:** put the extracted orchestrator (a) as a top-level function
> in `whole_genome.py` that `run_simulation` and `suite.main` both call, or (b) in a new
> `rigel/sim/orchestrator.py`. (a) is less churn; (b) is cleaner separation if `whole_genome.py`
> (2343 lines) should shed orchestration. **Leaning (a)** — fewer moving parts, and whole_genome
> already owns `run_simulation`.

## Callers / blast radius (all must keep working)

- `scripts/sim/simulate_suite.py` → `suite.main` (unchanged signature).
- `scripts/sim/simulate_reads.py` → `whole_genome.main` (unchanged).
- Tests: `test_sim_capture.py` (imports `SuiteCaptureSpec`, `_suite_capture_specs`,
  `capture_paired_condition_seed`, `_load_suite_config` from suite — these **stay in suite**),
  `test_whole_genome_sim_config.py`, `test_sim_gdna_overdispersion.py`, `test_oracle_bam.py`
  (import whole_genome internals — unaffected). No skills / CLI / CI entry points invoke these.

## Regression-safety strategy (the gate)

This refactor **must be output-identical**. Before/after proof:

1. Run a small representative suite **on the current code** into a temp dir (a few conditions
   spanning gdna × ss × capture × overdispersion), capture `manifest.json` + each condition's
   `truth_*` + BAM record counts.
2. Refactor.
3. Re-run the same suite into a second temp dir; **diff** manifests + truth + counts. Require
   byte-identical manifest/truth (BAMs may differ only by allowed nondeterminism — pin seeds so
   they match exactly). Any diff is a regression to fix before merge.
4. Plus the existing `pytest tests/` (esp. `test_sim_capture.py`) green.

## Phased execution

- **P1** — extract the orchestrator function from `run_simulation` (pure move; `run_simulation`
  calls it). Prove `simulate_reads` output unchanged.
- **P2** — refactor `suite.main` to build the config + call the orchestrator instead of its own
  loop; delete suite's duplicate loop (`1179–1313`) and its duplicate axis parsing where it now
  comes from the shared config. Prove `simulate_suite` output unchanged (the gate above).
- **P3** — collapse the duplicated overdispersion-axis parse/validate/label into one helper used by
  both config paths; confirm the axis is wired in exactly one place.

## Effort / risk

**Medium-high effort, medium risk** — the logic is well-understood and one-directional, but suite
output identity (seeding, capture grouping, manifest schema) is the real risk; the before/after
diff gate is what makes it safe. ~300 LOC consolidated; no behavior change intended.
