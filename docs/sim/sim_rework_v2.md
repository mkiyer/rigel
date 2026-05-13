# Simulation Rework Plan

Status: implementation-ready plan, 2026-05-12.

This plan is ready to implement. Start with Phase 0 and do not begin the
additive nRNA behavior change until the behavior-preserving phases have pinned
the current read names, truth parsing, manifest fields, CIGAR projection, and
oracle BAM output.

## Scope

The rework covers only simulation, synthetic benchmark generation, manifests,
truth parsing, benchmark analysis, and the command-line entry points for those
tasks. It must not change Rigel quantification, calibration, EM behavior, output
schemas, or native C++ code.

## Goals

1. Move reusable simulator logic into `src/rigel/sim/`.
2. Keep `scripts/sim/` as thin command-line wrappers.
3. Preserve the existing no-nRNA synthetic suite as the regression baseline.
4. Add nRNA-positive synthetic scenarios where nascent RNA is additive relative
   to mature RNA.
5. Make truth, manifests, condition discovery, and analysis explicit enough that
   molecular truth and fragment truth cannot be confused.
6. Keep small in-memory test scenarios and whole-genome file-backed benchmarks
   as separate frontends that share the same core mechanics.
7. Make simulation deterministic under serial and parallel execution.
8. Replace obscure script names with entry points whose names state exactly what
   they generate, evaluate, or diagnose.

## Current Problems

### 1. Two frontends duplicate core mechanics

`ReadSimulator` / `OracleBamSimulator` in `src/rigel/sim/` and
`WholeGenomeSimulator` in `scripts/sim/sim.py` independently implement several
mechanics that must not diverge:

| Mechanic | Small scenario path | Whole-genome path |
|---|---|---|
| Reverse complement | `genome.py` | inline in `sim.py` |
| Transcript-space projection | `oracle_bam.py` | inline in `sim.py` |
| Pre-mRNA projection | `oracle_bam.py` | inline in `sim.py` |
| CIGAR construction | `oracle_bam.py` | inline in `sim.py` |
| BAM record construction | `oracle_bam.py` | inline in `sim.py` |
| SAM flag constants | `oracle_bam.py` | inline in `sim.py` |
| Truncated-normal fragment lengths | `reads.py` | inline in `sim.py` |
| Sequence extraction | `reads.py` | inline in `sim.py` |
| Read-name truth parsing | `scenario.py` | `run_rigel_analysis.py` |

This is the highest-risk architecture problem. A CIGAR, BAM flag, orientation,
or read-name fix in one path does not automatically reach the other.

### 2. `scripts/sim/sim.py` is a monolith

At roughly 1900 lines, `sim.py` mixes configuration dataclasses, YAML parsing,
abundance assignment, nRNA spike-in logic, genome caching, sequence extraction,
FASTQ writing, oracle BAM construction, fork-based workers, condition-grid
orchestration, manifest writing, and CLI parsing. It is difficult to test or
replace one behavior without moving unrelated pieces.

### 3. The synthetic runner uses an importlib workaround

`scripts/sim/run_synthetic_sim.py` loads `scripts/sim/sim.py` with
`importlib.machinery.SourceFileLoader` because the script name collides with the
`scripts/sim/` package-like directory. The package should own the implementation
so scripts can import normal `rigel.sim.*` modules.

### 4. Analysis is hard-coded to the old condition list

`scripts/sim/evaluate_suite.py` contains a fixed no-nRNA condition list and
assumes one truth-abundance file. New nRNA-positive conditions should require no
analysis-code edits. Analysis should discover conditions and truth files from
the manifest.

### 5. nRNA currently uses a closed-total fraction model

The whole-genome simulator's current nRNA spike-in partitions a fixed total RNA
pool:

```text
total_rna ~ base abundance
nrna_frac ~ configured fraction
mRNA = total_rna * (1 - nrna_frac)
nRNA = total_rna * nrna_frac
```

This is not the right control for the synthetic benchmark. Increasing nRNA also
reduces mature RNA truth, so mature RNA error and nascent contamination are
confounded.

### 6. Truth is too implicit

`truth_abundances.tsv` is molecular abundance truth, not fragment truth. True
fragment origin must come from oracle BAM / FASTQ read names, and generated pool
sizes must be stored in the manifest. The rework must make these distinctions
visible in the APIs and report output.

### 7. Runtime and parallelism are concerns, not the first correctness task

The existing simulator parallelizes within a condition. Running all expanded
conditions sequentially will be slower, but condition-level parallelism must be
added only after deterministic condition generation and manifest-driven analysis
are correct.

### 8. Script entry points are named by accident, not by purpose

The current `scripts/sim/` names do not communicate which layer each script owns:

| Current script | Current role | Problem |
|---|---|---|
| `sim.py` | Generic whole-genome read generation from config. | Too generic; sounds like the main script but is also a monolith. |
| `run_synthetic_sim.py` | Synthetic mini-genome benchmark-suite generation. | The `run_` prefix hides that this creates a whole suite and manifest. |
| `generate_synthetic_genome.py` | Standalone mini-reference construction. | Normal suite generation should own reference creation; this should not be a second public workflow. |
| `synthetic_sim_sweep.py` | Small locus-level diagnostic sweep using `Scenario`. | Name makes it sound like the synthetic mini-genome suite, but it is a separate diagnostic tool. |
| `run_rigel_analysis.py` | Runs Rigel on simulated conditions and compares outputs to truth. | It is not simulation generation; it is suite evaluation and should say so. |

The cleanup should not merely shrink these files while preserving the same
ambiguous names. The final CLI surface should distinguish read simulation,
suite generation, suite evaluation, and locus diagnostics.

## Design Decisions

### Shared core, separate frontends

Do not force one simulator class. Keep these two frontends:

- `ReadSimulator` / `OracleBamSimulator`: in-memory, small scenario, unit-test
  oriented, backed by `MutableGenome`.
- `WholeGenomeSimulator`: file-backed, high-throughput benchmark generator,
  backed by `pysam.FastaFile` and sequence caches.

They share projection, CIGAR, BAM, orientation, and truth helpers through
package modules. This gives one source of truth for mechanics without collapsing
two different use cases into an awkward inheritance hierarchy.

### Additive nRNA

Mature RNA abundance is fixed. Nascent RNA is added on top:

```text
mrna_abundance = base_mrna_abundance
nrna_abundance = base_mrna_abundance * nrna_ratio
```

Single-exon transcripts always have `nrna_abundance = 0.0`, because an unspliced
single-exon transcript is not identifiable as nascent RNA in this benchmark.

### Explicit truth layers

The simulator and analysis must keep three truth layers separate:

| Truth layer | Example location | Meaning |
|---|---|---|
| Molecular truth | `truth_abundances_nrna_<label>.tsv` | Transcript-level mRNA and nRNA abundance. |
| Generated pool truth | `manifest.json` condition fields | Requested and actual `n_mrna`, `n_nrna`, `n_gdna`. |
| Fragment origin truth | FASTQ / oracle BAM read names | Per-fragment mRNA, nRNA, or gDNA origin. |

Analysis should use molecular truth for abundance recovery and read-name truth
for leakage / assignment metrics.

### Determinism before speed

Every condition must be reproducible independent of process scheduling. The
suite config should carry one `base_seed`, and condition generation should derive
stable sub-seeds from string keys:

```text
abundance_seed = stable_seed(base_seed, "abundance")
condition_seed = stable_seed(base_seed, condition_name)
mrna_seed      = stable_seed(condition_seed, "mrna")
nrna_seed      = stable_seed(condition_seed, "nrna")
gdna_seed      = stable_seed(condition_seed, "gdna")
strand_seed    = stable_seed(condition_seed, "strand")
error_seed     = stable_seed(condition_seed, "errors")
```

Use a stable hash implementation, not Python's process-randomized `hash()`.
Parallel execution must consume the same seeds as serial execution.

### Package implementation, script wrappers

All reusable simulation, truth, manifest, suite, and analysis logic belongs in
`src/rigel/sim/`. Files in `scripts/sim/` should only parse arguments, call the
package API, and print high-level progress.

### Clear script names

Do not keep wrappers with the old names. The rework deletes the obsolete public
script names rather than preserving compatibility shims.

Script naming rules:

- No generic `sim.py`.
- No `run_` prefixes; they say nothing about the artifact being produced.
- Names should identify the product: reads, suite, evaluation, or locus sweep.
- Reference generation is part of suite construction unless a reference-only
   debugging option is explicitly requested.
- Locus-level diagnostics must be named separately from whole-suite generation.

## Target Package Layout

| Module | Responsibility |
|---|---|
| `rigel.sim.genome` | `MutableGenome`, random DNA helpers, FASTA I/O, `reverse_complement`. |
| `rigel.sim.annotation` | `GeneBuilder`, transcript dataclasses, splice-motif injection, GTF writing. |
| `rigel.sim.bam` | SAM flags, paired-end orientation helpers, BAM record creation, transcript/pre-mRNA projection, CIGAR creation. Export public helper names such as `make_bam_record`, `transcript_to_genomic_blocks`, and `blocks_to_cigar`. |
| `rigel.sim.reads` | Small FASTQ simulator for test scenarios. Uses `MutableGenome` and shared helpers. |
| `rigel.sim.oracle_bam` | Oracle BAM writer for test scenarios. Uses shared helpers. |
| `rigel.sim.scenario` | End-to-end small scenario orchestration. Uses shared truth parsing. |
| `rigel.sim.benchmark` | Small scenario accuracy summaries. Uses shared truth parsing where needed. |
| `rigel.sim.truth` | Read-name parsing, `Origin` dataclass, truth abundance TSV I/O, fragment-origin aggregation from FASTQ/BAM. |
| `rigel.sim.manifest` | Manifest dataclasses, condition records, condition naming, manifest read/write, condition discovery. |
| `rigel.sim.abundance` | Random abundance assignment, abundance-file loaders, truth writing, additive nRNA ratio assignment. |
| `rigel.sim.config` | Config dataclasses and YAML parsing for simulation, gDNA, nRNA, and suite profiles. |
| `rigel.sim.whole_genome` | `WholeGenomeSimulator`, genome cache, FASTQ buffers, shard helpers, high-throughput read generation. |
| `rigel.sim.synthetic_genome` | Synthetic mini-genome generator currently implemented in `generate_synthetic_genome.py`. |
| `rigel.sim.suite` | Condition grid orchestration, smoke/full profiles, deterministic seed derivation, manifest writing. |
| `rigel.sim.analysis` | Manifest-driven calibration, abundance, fragment-origin, and acceptance report generation. |

## Final Script Surface

The final public scripts in `scripts/sim/` should be small and intentionally
named:

| Final script | Purpose | Replaces |
|---|---|---|
| `simulate_reads.py` | Generate FASTQ and optional oracle BAM from an existing FASTA/GTF/config. This is the low-level whole-genome read simulator. | `sim.py` |
| `simulate_suite.py` | Generate a complete named simulation benchmark suite: synthetic reference as needed, condition grid, reads, oracle BAMs, truth files, and manifest. | `run_synthetic_sim.py`; normal use of `generate_synthetic_genome.py` |
| `evaluate_suite.py` | Run Rigel quantification on a simulated suite and produce manifest-driven accuracy, calibration, and fragment-origin reports. This is evaluation, not simulation generation. | `run_rigel_analysis.py` |
| `locus_sweep.py` | Run small locus-level combinatorial diagnostics using `Scenario`. This remains separate from whole-genome suite generation. | `locus_sweep.py` |

Do not keep a standalone public `generate_synthetic_genome.py` in the final
surface. The implementation should live in `rigel.sim.synthetic_genome`, and
`simulate_suite.py` should call it when the suite needs a synthetic reference.
If reference-only generation remains useful for debugging, expose it as
`simulate_suite.py --reference-only` or a subcommand, not as a second top-level
workflow.

Old names should not remain as public shims. Callers must move to the final
script names.

## nRNA Schema

### New schema

Checked-in benchmark configs should use additive ratios:

```yaml
nrna:
  mode: additive_ratio
  ratios: [0.0, 0.1, 0.5, 1.0]
  ratio_labels: [none, low, med, high]
```

`ratio_labels` must have the same length as `ratios`. Condition names use the
label, not the raw float.

Old `nrna.fracs` closed-total configs are removed. If `nrna.fracs` or
`nrna.frac_labels` appears in a config, parsing fails and the config must be
migrated to additive ratios.

## Additive Fragment Counts

For the synthetic mini-suite, fragment counts are explicit and additive:

```text
n_mrna = n_rna_fragments
n_nrna = round(n_mrna * nrna_ratio)
n_gdna = round(n_mrna * gdna_rate)
total  = n_mrna + n_nrna + n_gdna
```

The manifest records configured ratios and actual generated counts:

```json
{
  "name": "gdna_med_ss_0.99_nrna_low",
  "gdna_label": "med",
  "gdna_rate": 0.5,
  "strand_specificity": 0.99,
  "nrna_label": "low",
  "nrna_ratio": 0.1,
  "nrna_mode": "additive_ratio",
  "truth_abundances": "truth_abundances_nrna_low.tsv",
  "n_mrna": 1000000,
  "n_nrna": 100000,
  "n_gdna": 500000,
  "n_total": 1600000
}
```

During analysis, `n_mrna` is the mature RNA count field. `n_rna` is retained in
new manifests only as the explicit total RNA pool (`n_mrna + n_nrna`).

## Read-Name Truth Parser Requirements

`rigel.sim.truth.parse_origin(qname)` should return a structured `Origin` with
at least:

```python
Origin(
    kind="mrna" | "nrna" | "gdna",
    transcript_id=str | None,
    ref=str | None,
    start=int | None,
    end=int | None,
    strand="+" | "-" | "f" | "r" | None,
    index=int | None,
)
```

It must handle both FASTQ names with `/1` or `/2` suffixes and BAM query names
without suffixes. It must also handle both existing gDNA formats:

```text
<transcript_id>:<start>-<end>:<strand>:<index>
nrna_<transcript_id>:<start>-<end>:<strand>:<index>
gdna:<start>-<end>:<strand>:<index>
gdna:<ref>:<start>-<end>:<strand>:<index>
```

The parser should not infer biological truth from quantification output. It
should only decode simulator-origin names.

## Condition Grid

Preserve the current baseline family:

```text
gdna_{none,low,med,equal,high} x ss_{0.99,0.50} x nrna_none
```

Add nRNA-positive labels:

| Label | Ratio | Meaning |
|---|---:|---|
| `none` | 0.0 | Existing no-nRNA baseline. |
| `low` | 0.1 | 10% as many nRNA fragments as mature RNA fragments. |
| `med` | 0.5 | Moderate nascent contamination. |
| `high` | 1.0 | Equal nRNA and mature RNA fragment counts. Stress condition. |

The full grid is 4 nRNA labels x 5 gDNA rates x 2 strand-specificity settings =
40 conditions. Add profiles rather than changing semantics:

| Profile | Conditions |
|---|---|
| `smoke` | Existing 10 no-nRNA conditions. |
| `full` | All 40 conditions. |
| explicit condition list | User-selected subset from manifest/config names. |

## Parallelism Strategy

Parallelism is useful, but it is not part of the semantic refactor. Implement it
after additive nRNA and manifest-driven analysis are correct.

### Within-condition workers

Preserve the existing within-condition shard workers. This is the production
path for large simulations. The refactor should move it into
`rigel.sim.whole_genome` without changing default behavior.

### Condition-level workers

Add condition-level workers only in the final phase. Default behavior should be
conservative:

- If `n_condition_workers > 1`, default each condition to `n_workers = 1`.
- If `n_workers > 1`, default `n_condition_workers = 1`.
- Allow nested parallelism only with an explicit opt-in and a warning about CPU
  oversubscription.
- On macOS, use `spawn` for condition-level workers.
- On Linux, fork-based copy-on-write can remain available for the existing
  within-condition implementation.

### Shared memory

Do not make shared memory part of the implementation plan unless profiling shows
sequence extraction dominates runtime after the correctness refactor. Treat it
as a future optimization, not as a required design dependency.

## Implementation Phases

### Phase 0: Characterize current behavior

Add tests before moving code.

1. Add read-name parser tests for mRNA, nRNA, small-scenario gDNA, whole-genome
   gDNA, FASTQ `/1` and `/2` suffixes, and BAM names without suffixes.
2. Add tests for transcript-to-genomic projection, pre-mRNA projection, CIGAR
   creation, and BAM flag/orientation helpers using representative plus- and
   minus-strand transcripts.
3. Add a tiny manifest fixture with condition-specific truth files.
4. Add a small `WholeGenomeSimulator.simulate_and_write` test that uses explicit
   `n_mrna`, `n_nrna`, and `n_gdna` counts and verifies oracle read-name counts.
5. Save a no-nRNA smoke-suite baseline report for comparison during the rework.

Exit criteria:

- Current tests pass.
- New characterization tests pass against the current implementation.
- No production behavior has changed.

### Phase 1: Extract shared truth, manifest, and BAM helpers

Do this in small behavior-preserving moves.

1. Create `rigel.sim.truth` and `rigel.sim.manifest` with tests, but do not
   change callers yet.
2. Create `rigel.sim.bam` with projection, CIGAR, flags, orientation, and BAM
   record helpers. Prefer public helper names over copied underscored names.
3. Migrate `OracleBamSimulator` and `Scenario` to the shared modules.
4. Migrate `WholeGenomeSimulator` to the shared modules.
5. Remove duplicated helper code only after both callers use the package helper.

Exit criteria:

- Tests pass after each caller migration.
- Selected oracle BAM records have unchanged query names, flags, reference
  starts, CIGARs, mate starts, template lengths, and `NH` tags.
- Current no-nRNA synthetic suite analysis still passes.

### Phase 2: Move whole-genome implementation into the package

1. Create `rigel.sim.config` and move simulation, gDNA, nRNA, suite, and YAML
   parsing dataclasses into it.
2. Create `rigel.sim.abundance` and move random abundance assignment,
   abundance-file loaders, abundance format detection, and truth abundance I/O.
3. Create `rigel.sim.whole_genome` and move `WholeGenomeSimulator`, genome
   cache, FASTQ buffers, gzip writers, batch extraction helpers, and shard
   helpers into it.
4. Replace `scripts/sim/sim.py` with `scripts/sim/simulate_reads.py`, a thin
   CLI wrapper around `rigel.sim.whole_genome`.
5. Delete `scripts/sim/sim.py` after `simulate_reads.py` imports the package
   implementation directly.

Exit criteria:

- CLI wrappers still run checked-in configs.
- No-nRNA smoke output is unchanged except for intentional manifest metadata
  additions.
- The generic read-generation entry point is named `simulate_reads.py`; no
   implementation remains in `sim.py`.
- Ruff passes for `src/`, `tests/`, and `scripts/sim/`.

### Phase 3: Implement additive nRNA ratios

1. Add `apply_nrna_ratio(transcripts, ratio)` to `rigel.sim.abundance`.
2. Add `NRNAConfig.mode`, `ratios`, and `ratio_labels` handling to
   `rigel.sim.config`.
3. Reject `nrna.fracs` / `nrna.frac_labels` configs with a migration error.
4. Update suite orchestration to compute explicit `n_mrna`, `n_nrna`, and
   `n_gdna` counts.
5. Update truth writing so each nRNA label has its own molecular truth file.
6. Update manifests to record `nrna_mode`, `nrna_ratio`, `truth_abundances`,
   `n_mrna`, `n_nrna`, `n_gdna`, and `n_total`.
7. Convert checked-in synthetic benchmark configs to additive ratio schema.

Exit criteria:

- `nrna_none` has `n_nrna = 0` and preserves mature RNA truth.
- `nrna_low` has the same mature RNA truth as `nrna_none` plus additive nRNA.
- Single-exon transcripts have zero nRNA truth and zero nRNA fragments.
- Any `fracs` config fails; checked-in configs use additive ratios only.

### Phase 4: Make analysis manifest-driven

1. Create `rigel.sim.analysis` and move calibration, abundance, fragment-origin,
   and acceptance reporting out of `run_rigel_analysis.py`.
2. Replace hard-coded conditions with manifest condition discovery.
3. Load each condition's molecular truth file from the condition manifest entry.
4. Use `rigel.sim.truth.parse_origin` for fragment-origin classification.
5. Preserve current no-nRNA acceptance checks exactly.
6. Add nRNA-positive metrics as report-only sections.
7. Replace `run_rigel_analysis.py` with `scripts/sim/evaluate_suite.py`, a thin
   CLI wrapper around `rigel.sim.analysis`.
8. Delete `run_rigel_analysis.py` after `evaluate_suite.py` imports the package
   implementation directly.

Exit criteria:

- Existing no-nRNA acceptance checks pass from manifest discovery.
- Analysis works when conditions are added, removed, or selected by name.
- nRNA-positive sections report mature RNA recovery, nRNA recovery, and leakage
  metrics without strict thresholds.
- The evaluator name makes its role clear: it consumes a simulation suite,
   runs Rigel if needed, and compares outputs against simulator truth.

### Phase 5: Build the suite layer and synthetic genome package

1. Create `rigel.sim.suite` for condition naming, profile expansion,
   deterministic seed derivation, suite execution, and manifest writing.
2. Create `rigel.sim.synthetic_genome` and move the mini-genome generator out of
   `generate_synthetic_genome.py`. Reuse `rigel.sim.annotation` where it fits;
   keep specialized generator structures only where they express real generator
   behavior.
3. Create `scripts/sim/simulate_suite.py`, a thin wrapper around
   `rigel.sim.suite` and `rigel.sim.synthetic_genome`.
4. Fold normal synthetic-reference generation into `simulate_suite.py`; do not
   keep `generate_synthetic_genome.py` as a separate public workflow.
5. Add `--profile smoke`, `--profile full`, `--conditions`, and, if useful,
   `--reference-only` support.
6. Delete `run_synthetic_sim.py` and `generate_synthetic_genome.py` after
   wrappers and tests no longer import them.
7. Rename `locus_sweep.py` to `locus_sweep.py` so the diagnostic locus
   sweep is not confused with whole-genome suite generation.
8. Delete dead code only after wrappers and tests no longer import it.

Exit criteria:

- A clean smoke run produces the existing 10 no-nRNA conditions.
- A clean full run produces all 40 conditions.
- The manifest describes every condition and truth file.
- No importlib workaround remains.
- Final user-facing names are `simulate_reads.py`, `simulate_suite.py`,
  `evaluate_suite.py`, and `locus_sweep.py`.

### Phase 6: Add optional condition-level parallelism

1. Add `n_condition_workers` to suite config / CLI.
2. Enforce conservative defaults that avoid nested worker oversubscription.
3. Use deterministic condition and pool seeds so serial and parallel execution
   generate the same per-condition output.
4. Use `spawn` for condition-level workers on macOS.
5. Keep shared memory out of scope unless profiling justifies it.
6. Remove any remaining old-name files (`sim.py`, `run_synthetic_sim.py`,
   `generate_synthetic_genome.py`, `run_rigel_analysis.py`,
   `locus_sweep.py`).

Exit criteria:

- Serial and condition-parallel smoke runs produce equivalent manifests,
  condition metadata, and acceptance reports.
- Runtime improves on the synthetic mini-suite without destabilizing macOS
  development runs.
- `scripts/sim/` contains only clearly named entry points and config files.

## Acceptance Checks

Preserve existing checks and thresholds for no-nRNA / gDNA scenarios:

| Check | Scope | Threshold |
|---|---|---:|
| `rho_ex / rho_ig` | gDNA-positive conditions | `>= 0.95` |
| `nRNA in nrna_none` | nRNA-none conditions | `<= 10` |
| `gDNA -> RNA leak` | gDNA-positive conditions | `<= 0.015` |

Add nRNA-positive sections as report-only metrics until the first full additive
suite has been reviewed:

| Metric | Scope | Initial action |
|---|---|---|
| Mature RNA recovery under nRNA | nRNA-positive conditions | Report only. |
| nRNA recovery | nRNA-positive conditions | Report only. |
| nRNA -> mRNA leak | nRNA-positive conditions | Report only. |
| mRNA -> nRNA leak | nRNA-positive conditions | Report only. |
| nRNA -> gDNA leak | nRNA-positive conditions | Report only. |
| gDNA -> nRNA leak | gDNA-positive, nRNA-positive conditions | Report only. |

After the first full run, set strict guardrails only for metrics that are stable
and statistically interpretable.

## Validation Commands

Use the project environment for all commands:

```bash
conda activate rigel
```

After each Python-only phase:

```bash
pytest tests/test_sim.py tests/test_sim_analysis.py tests/test_oracle_bam.py -v
pytest tests/scenarios/test_nrna_double_counting.py -v
pytest tests/test_golden_output.py -v
ruff check src/ tests/ scripts/sim
```

Before declaring the rework complete:

```bash
pytest tests/ -v
ruff check src/ tests/ scripts/sim
```

End-to-end smoke validation:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /tmp/rigel_sim_smoke --profile smoke -j 4

python scripts/sim/evaluate_suite.py \
  --sim-base /tmp/rigel_sim_smoke
```

End-to-end full validation:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /tmp/rigel_sim_full --profile full -j 4

python scripts/sim/evaluate_suite.py \
  --sim-base /tmp/rigel_sim_full
```

Use fresh output directories for the first additive-ratio runs so old
`nrna_none` outputs are not mixed with new manifest semantics.

## Risks And Mitigations

| Risk | Mitigation |
|---|---|
| Refactoring changes BAM flags, CIGARs, or coordinates. | Phase 0 pins representative records; Phase 1 migrates one caller at a time and compares selected BAM fields. |
| nRNA ratio semantics invalidate old comparisons. | Preserve the no-nRNA smoke baseline; reject all `fracs` configs and require additive ratios. |
| Molecular truth is confused with fragment truth. | Separate APIs and manifest fields for molecular truth, generated pool counts, and read-name origin truth. |
| Seed drift hides regressions. | Derive stable per-condition and per-pool seeds from string keys; require serial and parallel equivalence for condition metadata and reports. |
| Parallelism destabilizes macOS or oversubscribes CPUs. | Defer condition-level workers to Phase 6; avoid nested workers by default; use `spawn` on macOS. |
| Runtime grows from 10 to 40 conditions. | Keep `smoke` profile for routine checks and `full` profile for complete benchmark runs. |
| nRNA-positive thresholds are arbitrary. | Report first, inspect the full additive run, then set guardrails based on observed stable behavior. |
| Circular imports appear during extraction. | Keep `truth`, `manifest`, and `bam` low-level; they should not import suite or analysis modules. |
| Old script names survive and keep the workflow murky. | Delete old public script names; keep only `simulate_reads.py`, `simulate_suite.py`, `evaluate_suite.py`, and `locus_sweep.py`. |

## Non-Goals

- No Rigel quantification, calibration, EM, or output changes.
- No new C++ code.
- No attempt to make the simulator biologically complete.
- No forced inheritance between `ReadSimulator` and `WholeGenomeSimulator`.
- No strict nRNA-positive acceptance thresholds before one manifest-driven full
  additive-ratio suite has been generated and reviewed.
- No shared-memory optimization unless profiling after the refactor justifies it.

## First Implementation Step

Begin with Phase 0. The first concrete patch should add tests around read-name
origin parsing and oracle BAM/projection behavior, then introduce
`rigel.sim.truth` and `rigel.sim.manifest` without changing callers. That gives
the later nRNA work a clean, tested surface and keeps the first review focused
on observable behavior rather than architecture churn.
