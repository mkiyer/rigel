# Rigel architecture

The current code map and data flow. For the scientific method see `METHODS.md`; for calibration
theory see `calibration/CALIBRATION_ARCHITECTURE.md`; for usage see `MANUAL.md`.

Rigel is a Bayesian RNA-seq quantifier that jointly models **mRNA**, **nascent RNA (nRNA)**, and
**genomic DNA (gDNA) contamination**. A single-pass C++ BAM scanner feeds a per-region
**calibration** stage (deconvolve the library into gDNA vs RNA) and a locus-level **EM** solver.
PyPI package `rigel-rnaseq`; import and CLI are `rigel`.

## Pipeline (`pipeline.py`)

Three stages, one BAM pass:

```
BAM ──scan_and_buffer──▶ FragmentBuffer + AccumulatorPayload + trained models
                              │
                              ├── calibrate ──▶ CalibrationResult (per-region gDNA/RNA mass + scalars)
                              │
                              └── quant_from_buffer ──▶ AbundanceEstimator (per-locus EM)
```

1. **Scan** (`scan_and_buffer`): the C++ htslib reader resolves fragments against the reference
   index, trains the strand and fragment-length models from unique mappers, buffers resolved
   fragments into a columnar `FragmentBuffer`, **and** deposits per-region/per-boundary fractional
   fragment mass into the C++ **accumulator** (4 channels: unspliced ±, spliced sense/antisense)
   → an `AccumulatorPayload`.

2. **Calibrate** (`calibration.calibrate`): a **bipartite region↔boundary belief-propagation sweep**
   over the accumulator payload that deconvolves each node's unspliced mass into the simplex
   `(f_rna₊, f_rna₋, f_g)`. Per the **count-zero-information** principle, a fragment count carries no
   intrinsic gDNA/RNA information; a node's composition comes from exactly three sources — the strand
   likelihood (the Beta-Binomial tilt of the per-strand counts, the only intrinsic signal; the count
   enters only as overdispersed Fisher information), cross-node imputation (neighbour density messages
   at the belief-free Poisson disagreement variance σ²_imp, fit once), and the global gDNA prior
   (population baseline + a trained Phase-2 KDE). `node_chain` builds the chain; `bp_solver.node_sweep`
   runs a **single forward-backward pass** (exact on the chain, a forest of linear paths — no
   fixed-point loop). Fits the library hyperparameters (`gdna_density_global`, `rna_sense_frac`, the
   gDNA/RNA strand overdispersions) plus the per-region / per-boundary-side gDNA/RNA mass. Output:
   `CalibrationResult`. See `calibration/CALIBRATION_ARCHITECTURE.md`.

3. **Quantify** (`quant_from_buffer`): bridges calibration → per-locus prior
   (`calibration.priors.assemble_priors`), scores fragments (`scan.FragmentRouter` → CSR
   `ScoredFragments`), builds loci by connected components (`locus.build_multi_loci`), partitions
   the CSR per locus, and runs the per-locus EM. Each locus is an independent subproblem with
   **`n_t + 1` components** (one per transcript row + one gDNA). The calibration prior enters as two
   per-locus Dirichlet scalars (`gdna_prior_count`, `rna_prior_count`).

## Python modules

| module | role |
|---|---|
| `cli.py` | CLI: `index`, `quant`, `sim`, `export` |
| `pipeline.py` | orchestrator: scan → calibrate → quant |
| `config.py` | frozen config dataclasses (`EMConfig`, `BamScanConfig`, `PipelineConfig`, `CalibrationConfig`, `FragmentScoringConfig`) |
| `index.py` | reference index build/load (feather from GTF+FASTA); `region_df` = the calibration region partition |
| `scoring.py` | fragment likelihood scoring (RNA/gDNA FL models, strand/coverage/splice penalties) |
| `scan.py` | `FragmentRouter` — scores the buffer → global CSR `ScoredFragments` |
| `scored_fragments.py` / `locus_partition.py` | CSR fragment container + per-locus partitioning |
| `locus.py` | locus / `MultiLocus` construction (connected components), nRNA spans |
| `estimator.py` | `AbundanceEstimator` — per-locus EM dispatch, output dataframes |
| `scan_payload.py` / `_accumulator.py` | `AccumulatorPayload` schema + Python `Accumulator` wrapper |
| `buffer.py` | memory-efficient columnar fragment buffer |
| `strand_model.py` / `frag_length_model.py` | model training from unique mappers |
| `splice.py` / `splice_blacklist.py` | splice-type encoding + artifact blacklist |
| `annotate.py` | optional annotated-BAM writer (per-fragment tags ZT/ZG/ZF/ZC/ZL/ZW) |
| `native.py` | interface to the C++ extension modules |

### Calibration package (`calibration/`)

| module | role |
|---|---|
| `calibrate.py` | the calibrator orchestrator: κ → strand overdispersions → build chain + geometry + statics → signature-binary `init_beliefs` → `bp_solver.node_sweep` (the single forward-backward pass) → region / boundary-side projections → `gdna_density_global` |
| `node_chain.py` | builds the bipartite **region↔boundary** node chain (genomic interleave + left/right adjacency) from the payload offsets; pure index arithmetic |
| `bp_solver.py` | the BP solver: the anchored population gDNA prior + the belief-free message precision `adjacent_disagreement_variance` (fit once) + the single forward-backward `node_sweep` (gDNA + per-strand RNA messages) + the region / boundary-side projections (`chain_region_deconv` / `chain_boundary_side_deconv`) |
| `node_geometry.py` | per-node chain primitives beneath the solver: per-face geometry, the belief / density types + `node_densities` (the message currency ρ = f·M/E), the static solver inputs, the signature-binary `init_beliefs` |
| `simplex_logodds.py` / `simplex.py` | the per-node ψ solve over the 2-simplex `(f₊, f₋, f_g)` in log-fraction space (`_solve_nodes_logodds_all`): strand mixture + sided spliced floor + node-class prior + gDNA & per-strand RNA imputation messages, plus the lattice primitives |
| `strand_likelihood.py` | `strand_loglik`: the two-component gDNA/RNA strand Beta-Binomial log-likelihood of a node over a `gdna_frac` grid (the only intrinsic gDNA/RNA signal; both overdispersions applied symmetrically) |
| `gdna_density_prior.py` | the trained **Phase-2** gDNA-density KDE prior (`GdnaDensityPrior.fit` over gDNA-clean training nodes); consumed by `bp_solver._kde_logprior` |
| `strand_deconv.py` | the per-node `NodeDeconv` result type + `boundary_side_seeds` (the exon↔intron / exon↔intergenic boundary-side seeds for the gDNA strand-overdispersion fit) |
| `substrate.py` | `CalibrationSubstrate` (payload → per-region contained/left/right views) + `BoundarySubstrate` (the boundary-node view) |
| `region_arrays.py` / `regions.py` | region geometry (`RegionArrays`, partition from `index.region_df`) |
| `signature.py` | 4-bit exon/intron×strand signature + `strand_class` (POS/NEG/NONE/AMBIG) |
| `strand_balance.py` / `strand_summary.py` | RNA strand *mean* `rna_sense_frac` (κ); `rna_strand_overdispersion` here is a QC-only thin-count diagnostic |
| `density_model.py` | per-region gDNA density via local boundary-anchored imputation on the raw unspliced counts (`count_gdna_frac`): the gDNA strand-overdispersion fit's seed selector + the signature-partition masks the global prior uses |
| `gdna_strand.py` | both strand Beta-Binomial overdispersions (shared component-agnostic MoM core), applied symmetrically so unstranded data is uninformative |
| `effective_length.py` / `capture_eff_length.py` / `fl.py` | FL-marginal effective lengths; capture-aware EM effective lengths (gDNA-enrichment IPR contraction); gDNA/RNA/global FL pmfs |
| `run_fill.py` | run/neighbour helpers (`same_ref_left_right`) for the boundary-side and chain geometry |
| `derive.py` | `gdna_density_global` (a QC scalar) from the deconvolved masses |
| `priors.py` | `assemble_priors`: `CalibrationResult` → per-locus Dirichlet prior + gDNA eff-len |
| `result.py` / `errors.py` | `CalibrationResult` schema + invariants; exceptions |

## C++ extensions (`src/rigel/native/`)

Five nanobind modules (C++17, `-O3`, LTO; OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | single-pass htslib parser, fragment grouping, model training, **fractional accumulator** |
| `_em_impl` | `em_solver.cpp` | per-locus EM (`n_t+1` components), connected components, OpenMP, Kahan/digamma; uses the `fast_exp.h` SIMD exp (AVX2/AVX-512/NEON) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring, bias correction (`-ffast-math`; no SIMD/`fast_exp`) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges_bind.cpp` | vendored cgranges interval overlap |

The accumulator must match the Python reference (`tests/native/_accumulator_reference.py`)
byte-for-byte. After editing `src/rigel/native/`, re-run `pip install --no-build-isolation -e .`
in the activated `rigel` conda env.

## nRNA architecture

nRNA components are not per-transcript shadows. Rigel builds a global table of unique nRNA spans
keyed by `(ref, strand, start, end)` and shares each span across transcripts with matching
coordinates. Each unique span is materialized as an ordinary **transcript row** in `index.t_df`
(flagged `is_nrna`/`is_synthetic`), so the per-locus EM treats it like any transcript component —
hence `n_t + 1` components per locus, not a separate mRNA/nRNA pair per transcript.

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim [options]
rigel export results/ --format tsv    # feather → TSV / Parquet
```

Input BAM must be name-sorted with the `NH` tag.
