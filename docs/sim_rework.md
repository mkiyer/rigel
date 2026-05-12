# Simulation Rework Plan

Status: draft plan, 2026-05-12. This is the starting point for the simulation
cleanup before implementation. The intent is to make the simulator easier to
reason about before adding nRNA-positive synthetic benchmark scenarios.

## Goals

1. Move reusable simulator logic into `src/rigel/sim` and keep `scripts/sim` as
   thin command-line wrappers.
2. Add synthetic scenarios where nascent RNA is nonzero.
3. Treat nascent RNA as additive relative to mature RNA, not as a subtraction
   from a fixed total RNA pool.
4. Preserve the existing no-nRNA synthetic suite as a regression baseline.
5. Make manifests, truth files, condition naming, and analysis work for both
   nRNA-none and nRNA-positive conditions without hard-coded condition lists.

## Current Inventory

### `src/rigel/sim`

This package is the reusable test-scenario framework. It should remain the
home for simulation primitives.

| File | Current role | Plan |
|------|--------------|------|
| `genome.py` | `MutableGenome`, random DNA, FASTA writing, reverse complement. | Keep. Share sequence helpers with whole-genome code. |
| `annotation.py` | `GeneBuilder`, splice-motif injection, GTF writing. | Keep. Reuse style in synthetic reference generation. |
| `reads.py` | Small paired FASTQ simulator with mRNA, nRNA, and gDNA pools. | Keep for tests; align nRNA ratio semantics. |
| `oracle_bam.py` | Perfect-alignment BAM writer for small scenarios. | Keep; deduplicate projection and BAM helpers. |
| `scenario.py` | End-to-end small scenario orchestration. | Keep; add ratio helpers only if useful for tests. |
| `benchmark.py` | Pipeline-vs-truth comparison for `Scenario` outputs. | Keep; use shared truth parsing. |

### `scripts/sim`

This folder currently mixes CLIs with substantial implementation.

| File | Current role | Problem |
|------|--------------|---------|
| `generate_synthetic_genome.py` | Generates the 10 Mb synthetic mini-genome and GTF. | Reusable generator exists only as a script. |
| `sim.py` | Whole-genome simulator, config parser, abundance assignment, FASTQ/BAM writer, manifest writer, condition grid. | Too many responsibilities; duplicates logic from `src/rigel/sim`; current nRNA fraction reduces mature RNA. |
| `run_synthetic_sim.py` | Hard-coded 10-condition no-nRNA suite runner. | Uses an importlib workaround, bypasses generic config flow, hard-codes conditions. |
| `run_rigel_analysis.py` | Quantification and analysis report runner. | Hard-coded condition list and `truth_abundances.tsv`; only naturally supports `nrna_none`. |
| `synthetic_sim_sweep.py` | Small locus combinatorial sweep using `Scenario`. | Useful diagnostic tool, but separate semantics and output schema. |

## Main Problems

### nRNA is currently a closed total-RNA fraction

The whole-genome simulator currently models random nRNA as:

```text
total_rna ~ LogUniform(min, max)
nrna_frac ~ Uniform(min, max)
mRNA = total_rna * (1 - nrna_frac)
nRNA = total_rna * nrna_frac
```

For the synthetic benchmark, this is the wrong control. Increasing nRNA changes
the mature RNA truth. The new benchmark should keep mature RNA fixed and add
nascent RNA on top:

```text
mrna_abundance = base_mrna_abundance
nrna_abundance = base_mrna_abundance * nrna_ratio
```

For the mini-genome suite, the fragment-count control should also be additive:

```text
n_mrna_fragments = configured mature depth
n_nrna_fragments = round(n_mrna_fragments * nrna_ratio)
n_gdna_fragments = round(n_mrna_fragments * gdna_rate)
total_fragments = n_mrna_fragments + n_nrna_fragments + n_gdna_fragments
```

The existing `WholeGenomeSimulator.simulate_and_write(..., n_mrna=..., n_nrna=...)`
already has most of the low-level support for separate mRNA and nRNA pools. The
orchestration layer does not use it yet.

### Truth is not explicit enough

`truth_abundances.tsv` is molecular abundance truth, not fragment truth. Previous
analysis mistakes came from treating it as fragment truth. The rework should make
this distinction unavoidable:

- molecular truth: condition-specific `truth_abundances_nrna_<label>.tsv`
- generated pool truth: manifest fields such as `n_mrna`, `n_nrna`, `n_gdna`
- fragment assignment truth: oracle BAM read names

Analysis should load the condition-specific truth file named in `manifest.json`.

### Core mechanics are duplicated

These pieces appear in both the small simulator and the whole-genome simulator:

- reverse complement helpers
- transcript-space to genomic block projection
- pre-mRNA coordinate projection
- CIGAR creation
- paired-end orientation and BAM flags
- read-name truth conventions
- oracle BAM record construction

These should become shared package helpers so a CIGAR, strand, or read-name bug
cannot diverge between test scenarios and benchmark scenarios.

### Analysis is hard-coded to the old suite

`run_rigel_analysis.py` has a fixed `CONDITIONS` list and reads one truth file.
nRNA-positive scenarios need manifest-driven discovery and origin-aware analysis:

- `gdna:` read names are true gDNA.
- `nrna_<transcript_id>:` read names are true nascent RNA.
- `<transcript_id>:` read names are true mature RNA.

## Target Architecture

Use a layered design: package APIs do the work, scripts parse arguments.

### Proposed package modules

| Module | Responsibility |
|--------|----------------|
| `rigel.sim.config` | Dataclasses for simulation params, abundance config, gDNA config, nRNA ratio config, condition metadata, suite config. |
| `rigel.sim.abundance` | Random abundance assignment, abundance-file loaders, additive nRNA ratio assignment. |
| `rigel.sim.truth` | Read-name parsing, fragment truth counters, abundance truth writing/loading. |
| `rigel.sim.manifest` | Manifest schema, condition naming, manifest read/write, condition discovery. |
| `rigel.sim.projection` | Transcript/pre-mRNA coordinate projection and CIGAR helpers. |
| `rigel.sim.orientation` | Paired-end orientation and BAM flag helpers for mRNA, nRNA, and gDNA. |
| `rigel.sim.synthetic_genome` | Mini-genome gene/isoform generator currently in `generate_synthetic_genome.py`. |
| `rigel.sim.whole_genome` | Fast whole-genome simulator currently embedded in `scripts/sim/sim.py`. |
| `rigel.sim.suite` | Condition-grid orchestration for whole-genome simulation suites. |
| `rigel.sim.analysis` | Manifest-driven analysis helpers currently embedded in `run_rigel_analysis.py`. |

### Proposed script state

| Script | Target state |
|--------|--------------|
| `scripts/sim/generate_synthetic_genome.py` | Thin wrapper around `rigel.sim.synthetic_genome`. |
| `scripts/sim/sim.py` | Thin generic whole-genome simulation CLI. |
| `scripts/sim/run_synthetic_sim.py` | Thin mini-suite runner using `rigel.sim.suite`; no importlib workaround. |
| `scripts/sim/run_rigel_analysis.py` | Thin analysis CLI using `rigel.sim.analysis`; condition discovery from manifest by default. |
| `scripts/sim/synthetic_sim_sweep.py` | Keep as a diagnostic locus sweep initially; align terminology later. |

## nRNA Ratio Semantics

New configs should use `ratios`, not `fracs`:

```yaml
nrna:
  ratios: [0.0, 0.1, 0.5, 1.0]
  ratio_labels: [none, low, med, high]
  mode: fragment_ratio
```

For transcript abundance truth:

```text
mrna_abundance = base_mrna_abundance
nrna_abundance = base_mrna_abundance * nrna_ratio
```

For single-exon transcripts:

```text
nrna_abundance = 0
```

Single-exon nRNA is physically indistinguishable from mature RNA, so this should
remain zero.

For the synthetic mini-genome benchmark, use explicit fragment counts:

```text
n_mrna = n_rna_fragments
n_nrna = round(n_mrna * nrna_ratio)
n_gdna = round(n_mrna * gdna_rate)
```

The manifest should record both configured ratios and actual counts:

```json
{
  "nrna_label": "med",
  "nrna_ratio": 0.5,
  "n_mrna": 1000000,
  "n_nrna": 500000,
  "n_gdna": 2000000
}
```

Old `nrna.fracs` can be accepted temporarily with a warning, but checked-in
benchmark configs should move to the additive ratio schema.

## Synthetic Benchmark Augmentation

Preserve the current baseline condition family:

```text
gdna_{none,low,med,equal,high}_ss_{0.99,0.50}_nrna_none
```

Add nRNA-positive labels:

| Label | Ratio | Meaning |
|-------|-------|---------|
| `none` | `0.0` | Current baseline. |
| `low` | `0.1` | 10% as many nRNA fragments as mature RNA fragments. |
| `med` | `0.5` | Substantial nascent contamination. |
| `high` | `1.0` | Equal nRNA and mature RNA fragment counts. Stress condition. |

Crossing four nRNA labels with five gDNA rates and two strand-specificity values
gives 40 conditions. If runtime becomes awkward, add `--conditions` or a
`--profile smoke/full` option rather than changing the schema.

The analysis report should add nRNA-positive sections for:

- mature RNA recovery under nRNA contamination
- nRNA recovery from `nrna_quant` and locus-level `nrna`
- nRNA-to-mRNA, mRNA-to-nRNA, nRNA-to-gDNA, and gDNA-to-nRNA leakage from oracle
  read names
- `nrna_none` remains near zero
- calibration density behavior stratified by `nrna_label`, `gdna_label`, and
  strand specificity

## Migration Plan

### Phase 0: Characterize before moving code

Add focused tests before refactoring:

1. read-name parsing for mRNA, nRNA, and gDNA origins
2. `WholeGenomeSimulator.simulate_and_write` separate-pool mode with `n_mrna`,
   `n_nrna`, and `n_gdna`
3. a tiny manifest fixture with condition-specific truth files

Exit: current tests still pass, and the high-risk mechanics are pinned.

### Phase 1: Extract manifest and truth helpers

1. Create `rigel.sim.truth` and `rigel.sim.manifest`.
2. Move condition naming, manifest read/write, abundance truth writing, and
   read-name origin parsing into those modules.
3. Update `benchmark.py`, `scripts/sim/sim.py`, and `run_rigel_analysis.py` to
   use the shared helpers without changing behavior.

Exit: current no-nRNA suite analysis still works and tests pass.

### Phase 2: Implement additive nRNA ratios

1. Add `NRNAConfig.ratios`, `ratio_labels`, and mode handling.
2. Implement `apply_nrna_ratio(transcripts, ratio)` with mature abundance
   preserved.
3. Update suite orchestration to pass explicit `n_mrna`, `n_nrna`, and `n_gdna`
   to the whole-genome simulator.
4. Keep `fracs` as a temporary compatibility path with a warning.

Exit: a small generated condition has exact manifest counts for all pools;
single-exon transcripts still produce no nRNA reads; `nrna_none` is unchanged
apart from additive metadata.

### Phase 3: Make analysis manifest-driven

1. Replace hard-coded `CONDITIONS` with manifest condition discovery.
2. Load each condition's truth file from the manifest.
3. Extend fragment assignment analysis to classify true mRNA, true nRNA, and
   true gDNA.
4. Preserve current no-nRNA acceptance checks exactly.
5. Add nRNA-positive metrics, but report them before setting strict thresholds.

Initial acceptance policy:

| Check | Applies to | Initial threshold |
|-------|------------|-------------------|
| `rho_ex/rho_ig` | gDNA-positive conditions | `>= 0.95`, as today. |
| `nRNA in nrna_none` | nRNA-none conditions | `<= 10`, as today. |
| `gDNA->RNA leak` | gDNA-positive conditions | `<= 0.015`, as today. |
| `nRNA recovery` | nRNA-positive conditions | report first; gate after first full run. |
| `mature RNA recovery under nRNA` | nRNA-positive conditions | report first; gate after first full run. |
| `nRNA->gDNA leak` | nRNA-positive conditions | report first; gate after first full run. |

Exit: old acceptance checks remain visible and passing; reports group by nRNA,
gDNA, and strand condition.

### Phase 4: Extract the whole-genome simulator

1. Move `WholeGenomeSimulator` into `rigel.sim.whole_genome`.
2. Move projection and BAM helpers into shared package modules.
3. Move abundance loaders into `rigel.sim.abundance`.
4. Reduce `scripts/sim/sim.py` to a CLI wrapper.

Exit: checked-in configs still run; no generated no-nRNA output changes except
intentional metadata additions.

### Phase 5: Convert the synthetic mini-suite runner

1. Rewrite `run_synthetic_sim.py` as a wrapper around `rigel.sim.suite`.
2. Remove the importlib workaround.
3. Move the condition grid into YAML or a package dataclass.
4. Add `--conditions` and `--profile smoke/full` if useful.

Exit: a clean output directory produces the existing no-nRNA conditions plus
the expanded nRNA conditions, and the manifest describes every output.

### Phase 6: Clean stale paths and docs

1. Convert `scripts/sim/configs/synthetic_mini.yaml` to the additive ratio
   schema.
2. Update `sim_example.yaml` comments to remove closed-total nRNA language.
3. Decide which old `fracs` examples remain compatibility examples.
4. Keep `synthetic_sim_sweep.py` for locus diagnostics, but clarify abundance
   ratio versus fragment ratio terminology.
5. Delete dead code only after scripts and tests no longer import it.

## Tests To Add Or Update

- `test_nrna_ratio_preserves_mature_abundance`
- `test_nrna_ratio_zeroes_single_exon_transcripts`
- `test_condition_manifest_records_pool_counts`
- `test_truth_parser_classifies_mrna_nrna_gdna`
- `test_analysis_loads_condition_specific_truth`
- `test_run_rigel_analysis_discovers_manifest_conditions`
- small whole-genome run with explicit `n_mrna`, `n_nrna`, and `n_gdna`
- regression that no-nRNA condition names and acceptance behavior do not change

## Validation Commands

Use the project environment for all commands:

```bash
conda activate rigel
```

After Python-only simulator changes:

```bash
pytest tests/test_sim.py tests/test_sim_analysis.py -v
pytest tests/scenarios/test_nrna_double_counting.py -v
pytest tests/test_golden_output.py -v
ruff check src/ tests/ scripts/sim
```

Before declaring the rework complete:

```bash
pytest tests/ -v
ruff check src/ tests/ scripts/sim
```

Synthetic suite validation after implementation:

```bash
python scripts/sim/run_synthetic_sim.py \
  --outdir /Users/mkiyer/Downloads/rigel_runs/sim_synthetic

python scripts/sim/run_rigel_analysis.py \
  --sim-base /Users/mkiyer/Downloads/rigel_runs/sim_synthetic
```

For the first expanded 40-condition run, use a fresh output directory or a
clearly named sibling directory so old `nrna_none` outputs are not mixed with
new manifest semantics.

## Risks And Mitigations

| Risk | Mitigation |
|------|------------|
| nRNA semantics invalidate comparison to old runs. | Preserve no-nRNA conditions exactly; add explicit `nrna_ratio` and pool counts for new conditions. |
| Abundance ratio and fragment ratio are confused. | Benchmark suite uses explicit fragment-count control; general simulation can document abundance-ratio mode separately. |
| Refactoring changes BAM flags or CIGARs. | Add focused tests before moving code; compare oracle read-name counts and selected BAM tags. |
| nRNA-positive thresholds are arbitrary. | Report first full additive-ratio run, then set empirical guardrails. |
| Runtime grows from 10 to 40 conditions. | Add condition filters and smoke/full profiles. |
| Old configs break. | Temporarily accept `nrna.fracs` with a warning and convert checked-in configs deliberately. |

## Non-Goals

- Do not change Rigel quantification, calibration, EM, or output semantics as
  part of this simulator refactor.
- Do not introduce new C++ code.
- Do not make the simulator biologically complete. The target is controlled,
  interpretable synthetic benchmarking.
- Do not set strict nRNA-positive acceptance thresholds before one
  manifest-driven additive-ratio suite has been generated and analyzed.

## Recommended First Implementation Step

Start with a behavior-preserving extraction:

1. Add `rigel.sim.truth` and `rigel.sim.manifest`.
2. Move condition naming, manifest read/write, and read-name origin parsing into
   those modules.
3. Update `run_rigel_analysis.py` to use manifest condition discovery while
   still passing on the current no-nRNA suite.

That gives the nRNA work a clean surface: new conditions can be added through
simulation config and manifest fields, not by teaching every analysis function
another hard-coded condition pattern.