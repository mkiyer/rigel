# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

Rigel is a Bayesian RNA-seq transcript quantification tool that jointly models mRNA, nascent RNA (nRNA), and genomic DNA contamination (gDNA). It uses a single-pass C++ BAM scanner, a per-region **calibration** stage that deconvolves the library into gDNA vs RNA, and a locus-level EM solver. Python package name is `rigel-rnaseq` on PyPI; the import and CLI are `rigel`.

> **Calibration (WIP — NOT ready to ship).** ⚠ **Read
> `docs/calibration/ROADMAP.md` before doing any calibration work — it is the single entry point.** Most docs
> in `docs/calibration/` were archived to `archive/` on 2026-07-24; only the handful named in the ROADMAP are
> live. Do NOT reference anything in `archive/`.
>
> The calibration stage deconvolves each node's **unspliced** fragment mass into the 2-simplex **(f_rna₊,
> f_rna₋, f_g)** — sense-RNA / antisense-RNA / gDNA (calibration models **only RNA-vs-gDNA**; the per-locus EM
> separates nascent from mature downstream) — by a **belief-propagation sweep** over a **bipartite
> region↔boundary node chain** (`node_chain` builds the chain; `bp_solver` is the solver). This is the
> **prior-free first pass, "pass-0"** — an approximation whose result trains a **gDNA hyperprior** that is then
> required to **re-solve** (esp. AMBIG nodes). Per the **count-zero-information** principle
> (`docs/calibration/CALIBRATION_ARCHITECTURE.md`, authoritative): a fragment count carries no intrinsic
> gDNA/RNA information, so a node's composition is set by exactly three sources — the **strand likelihood**
> (Beta-Binomial tilt; constrains only the tilt, so AMBIG nodes get NO f_g from it), **cross-node imputation**
> (neighbour messages), and the **population gDNA prior** (the hyperprior, WIP).
>
> **⚠⚠ THE TWO STANDING DEBTS (P1e and `ω_graft`) — read `docs/calibration/P1D_P1E_DEBTS.md`.** Measured
> 2026-07-27 to be **REDUNDANT** (they price opposite directions of one failure: the graft's inequality
> used as an equality), and **both are RETAINED deliberately** until TSS/TES enter the region map (P1g),
> which is the structural fix both stand in for. ⚠ The 10 Mb synthetic **overstates P1e's damping by
> 26–300×** against real cfRNA (31.1 % of solver precision vs 0.10–1.18 %), and `ω_graft` fits **10×
> apart on two real samples** — so neither term's toy A/B is a production number.
>
> **⚠⚠ STANDING DEBT — P1e (`conservation surprise`) PRICES A BIAS AS A VARIANCE.** It damps a message by
> the unexplained part of `δ = log(M/S)`, but on a large share of its firing mass `δ` is systematic, not
> scatter (bias share 53–77 % on graft × one-component, 98.9–99.2 % at intergenic destinations), and
> **90–100 % of that damping lands on solvable nodes** — it is not inert. A variance cannot move a mode
> toward truth and it never shrinks. Landed anyway (the only change to improve accuracy AND honest precision
> together) and **scoped to the licensed direction `δ < 0`**; `RIGEL_P1E_OFF=1` ablates. The magnitude is
> also not what works — `δ` selects WHICH message to distrust and the calibration adds nothing.
> **When the bias strata are diagnosed, this term must SHRINK.** `variance_ledger.md` §6.
>
> **⚠⚠ STANDING DEBT — `ω_graft` (P1d) COMPENSATES FOR A STRUCTURAL FAILURE, AND MUST BE RE-DERIVED.**
> The region/boundary map has **no TSS/TES**, so the solver cannot tell a splice junction from a transcript
> terminus. `calibration/enrichment_frame.graft_premise_logvar` prices the graft's premise error with **one
> library-wide fitted scalar** standing in for a quantity measured to split **≥30×** on exactly that missing
> bit (`ω̂` 1.7–1.9 at termini vs 0.04–0.06 at junction-only seams). It works because over-charging a variance
> is cheap and under-charging is expensive — **not because it is right**, and it is **expected to be fragile
> on real data**. **When TSS/TES enter the region map (P1g), `ω` MUST be re-derived per structural class** —
> same equation, one scalar per class; the partial-pooling block in that function is the plug-in point. Do
> not treat the fitted value as a measured constant of the library.
>
> **PERFORMANCE (2026-07-27): calibration is ~4× faster at genome scale (266 s → 67 s), bit-identically.**
> `docs/calibration/PERFORMANCE.md` is authoritative; its §1 is the lesson — **the 10 Mb synthetic ranks
> the hotspots BACKWARDS**, so profile on the cached real cfRNA payloads, never the toy. The relay is
> scalar-native (`*_scalar` twins + Python-list operands) and both ψ solves are cache-blocked.
>
> **`bp_solver.py` holds ONE solve path, and as of 2026-07-26 it reads no environment variables** — every
> prototyping A/B gate has been removed and the production path has converged. The per-node
> **initialization self-solve** lives in `node_init.py` (`build_node_init` — the four sources:
> measured/intergenic, intron-factory, strand-deconv, unsolved-default). The **message-variance model is
> COMPLETE**, not a blocker: laws M1–M11 are derived, MC-validated in `scripts/debug/message_variance_mc.py`
> and A/B-won (`docs/calibration/message_variance_derivation.md`, audited in
> `docs/calibration/variance_ledger.md`). **The solver is now BP-legal** — all three third-node dependences
> (`_far`, the destination-own plug-in, and the two-ρ-iteration reframe loop) were removed on 2026-07-26. The
> **gDNA intron factory** is shipped (`intron_factory=True` default): it is the **intron special case** of the
> generic **density deconvolution** (`density_deconv.py`) — deconvolve counts into gDNA + RNA against a gDNA
> density prior (the intergenic node distribution, for introns), carrying its NB-derived precision
> (`density_deconv.density_factor_precision`).
> Fractional-accumulator spec in `docs/accumulator/00_design.md`. The full pipeline (scan → calibrate →
> **quant**) runs end-to-end: `quant_from_buffer` + `calibration.priors.assemble_priors` are wired in
> `run_pipeline`.

## Build & Development

> **ALL build, test, and lint commands MUST run inside the activated `rigel`
> conda environment.** It contains the full toolchain — htslib, C/C++
> compilers, `scikit-build-core`, `nanobind` — and the C++ build locates
> htslib via `$CONDA_PREFIX` (so the env must be active, not merely on `PATH`).
> In a non-interactive shell, prefix commands with:
> ```bash
> source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
> ```

```bash
# One-time environment setup (conda-forge + bioconda channels)
mamba env create -f mamba_env.yaml

conda activate rigel    # ← required before every build / test / lint

# Build and install (editable, C++ compiled via scikit-build-core + nanobind)
pip install --no-build-isolation -e .

# Install with dev extras (pytest)
pip install --no-build-isolation -ve ".[dev]"
```

After changing any C++ code in `src/rigel/native/`, you must re-run `pip install --no-build-isolation -e .` (in the activated `rigel` env) to recompile.

## Testing

```bash
pytest tests/ -v                          # all tests
pytest tests/calibration/ -v              # calibration unit tests
pytest tests/test_em_impl.py::test_name   # single test function
pytest tests/ --update-golden             # regenerate golden output files
pytest tests/ --cov=rigel                 # with coverage
```

Golden output regression tests live in `tests/golden/` (feather, TSV, JSON). Use `--update-golden` to regenerate after intentional output changes.

## Linting

Ruff is configured in `pyproject.toml`: Python 3.12 target, 100-char line length. `scripts/` is also
linted (with relaxed diagnostic-script idioms for `scripts/debug/**` — see `scripts/README.md`).

```bash
ruff check src/ tests/ scripts/
ruff format src/ tests/
```

## Architecture

### Pipeline (`pipeline.py`)

1. **BAM Scan** (`scan_and_buffer`): C++ htslib single-pass BAM reader. Resolves fragments against the reference index, trains strand/fragment-length models from unique mappers, buffers resolved fragments into a columnar `FragmentBuffer`, **and** deposits per-region/per-boundary fractional fragment mass into the C++ **accumulator** (4 channels: unspliced ±, spliced sense/antisense) → an `AccumulatorPayload`.

2. **Calibration** (`calibration.calibrate`): the **bipartite belief-propagation sweep** over the accumulator payload (see the Project Summary above). It deconvolves each node's unspliced mass into **gDNA vs RNA** and fits the library hyperparameters — `rna_sense_frac` (κ), the gDNA/RNA strand Beta-Binomial overdispersions, and `gdna_density_global` (a QC scalar) — plus the per-region and per-boundary-side deconvolved gDNA/RNA mass. Output: `CalibrationResult`.

3. **Quantification** (`quant_from_buffer`): bridges calibration → per-locus prior (`calibration.priors.assemble_priors`), scores fragments (`scan.FragmentRouter` → CSR `ScoredFragments`), builds loci via connected components (`locus.build_multi_loci`), partitions the CSR per locus (`locus_partition.partition_and_free`), and runs per-locus EM. Each locus is an independent subproblem with **`n_t + 1` components** — one per transcript row (annotated mRNA and nRNA spans alike) plus one gDNA component. The calibration prior enters as **two per-locus Dirichlet scalars** (`alpha_rna_add`, `alpha_gdna_add`) that set the gDNA-vs-RNA split; the EM distributes RNA mass among the compatible transcripts.

### Python Module Roles

- `cli.py` — CLI entry point with subcommands: `index`, `quant`, `sim`, `export`
- `pipeline.py` — orchestrator connecting scan → calibrate → quant
- `config.py` — frozen dataclasses (`EMConfig`, `BamScanConfig`, `PipelineConfig`, `CalibrationConfig`, `FragmentScoringConfig`)
- `index.py` — reference index build/load (feather files from GTF+FASTA); `region_df` defines the calibration region partition
- `scoring.py` — fragment likelihood scoring (`FragmentScorer`, RNA/gDNA FL models + strand/coverage/splice penalties)
- `scan.py` — `FragmentRouter` (scores the buffer → global CSR `ScoredFragments`)
- `scored_fragments.py` / `locus_partition.py` — CSR fragment container + per-locus partitioning
- `locus.py` — locus / `MultiLocus` construction (connected components), nRNA spans
- `estimator.py` — `AbundanceEstimator` (per-locus EM dispatch, output dataframes)
- `scan_payload.py` / `_accumulator.py` — `AccumulatorPayload` schema + Python `Accumulator` wrapper
- `buffer.py` — memory-efficient columnar fragment buffer
- `strand_model.py` / `frag_length_model.py` — model training from unique mappers
- `splice.py` / `splice_blacklist.py` — splice-type encoding + artifact blacklist
- `native.py` — public interface to the C++ extension modules

**Calibration package (`calibration/`):**

- `calibrate.py` — the calibrator orchestrator: RNA strand balance κ → gDNA/RNA strand overdispersions (seeded from `density_model.node_gdna_density` raw counts + `strand_deconv.boundary_side_seeds`) → build chain + geometry + statics → signature-binary `init_beliefs` → **`bp_solver.node_sweep`** (the single forward-backward pass; `σ²_imp` + every global-prior input fit once before the pass, no per-pass refit) → `chain_region_deconv` + `chain_boundary_side_deconv` → `gdna_density_global`. `sweep_n_grid`
- `node_chain.py` — builds the bipartite **region↔boundary** node chain (genomic B-R-B-…-B interleave + left/right adjacency) from the payload offsets; pure index arithmetic, no C++
- `bp_solver.py` — **the BP solver** (unified/composition, one path): the single forward-backward `node_sweep` (build the per-node self-solve → reframe/route/÷M relay + combine → ψ solve) + the region / boundary-side projections (`chain_region_deconv` / `chain_boundary_side_deconv`).
- `node_init.py` — **the pass-0 per-node INITIALIZATION**: `build_node_init` runs the message-free self-solve and assembles each node's own `(density, precision)` from the four sources of `docs/calibration/variance_model_concepts.md` — MEASURED (intergenic → Poisson precision), DENSITY DECONVOLUTION (`density_deconv`, `τ_λ` from NB curvature; intron factory = special case), STRAND DECONV (`strand_evidence`, the I_strand deadband), UNSOLVED default (100% gDNA, ZERO precision). Pure precision arithmetic (`own_composition_logvar`, `own_precision`); unit-tested per source in `tests/calibration/test_node_init.py`
- `density_deconv.py` — the generic **density deconvolution** primitive: fit a gDNA density prior NegBinom from a pure-gDNA pool (`fit_gdna_background`) + the per-node λ-factor `log NegBinom(f_g·C; ρ_bg·E_g, α_eff)` (`density_lambda_factor`) + its curvature precision (`density_factor_precision`). The **intron factory** is `fit_intron_background` (the gDNA prior = the intergenic node distribution)
- `node_geometry.py` — the per-node chain primitives beneath the solver: per-face geometry (`build_node_geometry`), the belief / density types, the static solver inputs (`build_node_statics`), the signature-binary `init_beliefs`, the pure geometry helpers `node_global_geometry` / `node_total_density`, and `_node_region_type`. Imports only lower layers (`node_chain`, `signature`, `effective_length`, `simplex_logodds`), so it sits below `bp_solver` / `node_init`
- `gdna_landscape.py` — **`GdnaLandscape` + `fit_gdna_landscape`: the Phase-2 gDNA-density hyperprior.** A weighted sum of zero-native per-node Poisson kernels rendered at a **population** resolution (nearest-neighbour spacing, `knn_widths`) — deterministic, no EM, no competition, so a capture-enriched minority cannot be competed away. Substrate selection lives in `calibrate._fit_gdna_hyperprior`; `strength` tempers the term (default 1). ⚠ Read `knn_widths` before touching the kernel: a *measurement* width (`1/√g`) collapses as ρ^(−1/2) and turns the enriched half of the landscape into a comb
- `npmle.py` — `DensityNPMLE`: **RETIRED from the hyperprior role** (superseded by `gdna_landscape`), kept for the Role-A enrichment fit + QC/diagnostics and toy injection. Do not delete before full production
- `substrate.py` — `CalibrationSubstrate` (payload → per-region contained / left / right views) + `BoundarySubstrate` (the per-boundary view: the bipartite graph's boundary nodes)
- `region_arrays.py` / `regions.py` — region geometry (`RegionArrays`, build partition from `index.region_df`)
- `signature.py` — region 4-bit strand/type signature + `strand_class` (POS/NEG/NONE/AMBIG)
- `density_model.py` — per-region gDNA **density** via local boundary-anchored imputation on the **raw** unspliced counts (`pos+neg`) + the signature-partition masks (`region_count_observable` / `boundary_count_observable`); returns `count_gdna_frac`, consumed only as the gDNA strand-overdispersion fit's seed selector + the partition masks the global prior uses
- `strand_balance.py` / `strand_summary.py` — **RNA** strand *mean*: `rna_sense_frac` (used by the decode). `StrandBalance.rna_strand_overdispersion` here is a QC-only thin-count power diagnostic (`1/(n_obs+3)`), distinct from the decode's RNA overdispersion (see `gdna_strand.py`)
- `gdna_strand.py` — **both** strand Beta-Binomial overdispersions (shared component-agnostic MoM core): `gdna_strand_overdispersion` (mean ½, fit from count-observable seed regions + boundary sides) and `rna_strand_overdispersion` (mean κ, fit from boundary-side spliced counts). Both applied symmetrically in `strand_likelihood` with the same default prior, so unstranded data is uninformative (see `docs/calibration/CALIBRATION_ARCHITECTURE.md`)
- `strand_likelihood.py` — `strand_loglik`: the two-component gDNA/RNA strand Beta-Binomial log-likelihood of one node over a `gdna_frac` grid (the only INTRINSIC gDNA/RNA signal; both overdispersions applied symmetrically)
- `strand_deconv.py` — the per-node `NodeDeconv` result type (the schema `bp_solver`'s projections + `simplex_logodds._solve_nodes_logodds_all` return) + `boundary_side_seeds` (the exon↔intron / exon↔intergenic boundary-side seeds for the gDNA strand-overdispersion fit, complementing the contained-region seeds under capture)
- `derive.py` — `gdna_density_global`: the library-average gDNA density (a QC scalar) from the deconvolved masses
- `effective_length.py` — the single owner of the FL-marginal effective lengths: `region_eff_length` (contained `E[max(0,L−ℓ)]`), `boundary_side_eff_length` (per-side crossing `E[min(ℓ,L)]`), the pooled-seam density divisor
- `capture_eff_length.py` — capture-aware **EM** effective lengths (`transcript_capture_eff_lengths`): contract each transcript's eff-len by the per-region gDNA-enrichment IPR over its region set, with an evidence-weighted shrinkage `w=G/(G+n_reg)` toward no-contraction on sparse gDNA
- `simplex.py` — per-node 2-simplex primitives (`_simplex_lattice`, `_mixture_strand_loglik`) shared by `simplex_logodds`
- `simplex_logodds.py` — the per-node ψ solve primitives (`_solve_nodes_logodds_all` + `_local_loglik_logodds` + the posterior-marginal helpers) over the 2-simplex `(f₊, f₋, f_g)` in **log-fraction** space: strand mixture + sided spliced floor + node-class prior (Jeffreys / global) + the gDNA & per-strand RNA imputation messages. `bp_solver` drives them over the chain
- `run_fill.py` — run/neighbour helpers (`same_ref_left_right`) for the boundary-side and chain geometry
- `fl.py` — gDNA / RNA / global fragment-length pmfs (empirical-Bayes smoothed)
- `result.py` — `CalibrationResult` schema + intrinsic invariants
- `priors.py` — `assemble_priors`: `CalibrationResult` → per-locus Dirichlet prior (+ the gDNA-component IPR effective length)
- `errors.py` — calibration exceptions

### C++ Extensions (`src/rigel/native/`)

Five nanobind modules compiled via CMakeLists.txt:

| Module | Source | Purpose |
|--------|--------|---------|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | Single-pass htslib BAM parser, fragment grouping, model training, **fractional accumulator** (per-region/boundary mass + FL pools) |
| `_em_impl` | `em_solver.cpp` | Per-locus EM (`n_t + 1` components), connected components, partition scatter helpers, OpenMP |
| `_scoring_impl` | `scoring.cpp` | Fragment likelihood scoring, bias correction (`-ffast-math`; no SIMD/`fast_exp`) |
| `_resolve_impl` | `resolve.cpp` | Fragment-to-transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges/cgranges_bind.cpp` | Vendored cgranges interval overlap library |

Key C++ details:
- C++17, compiled with `-O3`, LTO enabled
- SIMD: `fast_exp.h` provides AVX2/AVX-512/NEON code paths for exp() — used by `_em_impl` only (`scoring.cpp` has no SIMD)
- EM solver uses Kahan compensated summation and digamma for VBEM mode
- OpenMP for parallel BAM scanning and batch locus EM
- `_scoring_impl` uses `-ffast-math`; `_em_impl` uses `-ffp-contract=fast` (preserves Kahan)
- The accumulator must match the Python reference (`tests/native/_accumulator_reference.py`) byte-for-byte

### nRNA Architecture

nRNA components are not per-transcript shadows. Rigel builds a global table of unique nRNA spans keyed by `(ref, strand, start, end)` and shares each span across transcripts with matching genomic coordinates, reducing redundant states in isoform-rich loci. Each unique span is materialized as an ordinary **transcript row** in `index.t_df` (flagged `is_nrna` / `is_synthetic`), so the per-locus EM treats it like any other transcript component — hence `n_t + 1` components per locus (all transcript rows + one gDNA), not a separate mRNA/nRNA pair per transcript.

## Test Infrastructure

- `tests/` (top level) + `tests/calibration/` (calibration units) + `tests/native/` (accumulator spec / byte-for-byte) — ~60 test modules; fixtures in `conftest.py` (synthetic mini GTF/FASTA: 3 transcripts, 2 genes)
- Golden output tests (`test_golden_output.py`) compare against `tests/golden/` reference files
- Scenario BAM files in `tests/scenarios/` and `tests/scenarios_aligned/`
- CI runs on Ubuntu + macOS, Python 3.12 + 3.13

## CLI Subcommands

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/   # generate synthetic test scenarios
rigel export results/ -f tsv               # convert feather outputs to TSV/Parquet
```

Input BAM must be name-sorted with `NH` tag.
