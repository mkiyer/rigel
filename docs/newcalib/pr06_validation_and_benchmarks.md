# PR 06 - Validation and Benchmarks

## Goal

Validate the rebuilt calibration system as a single production path for ordinary RNA-seq,
strand-specific RNA-seq, unstranded RNA-seq, and hybrid capture RNA-seq.

This PR removes temporary skips, deletes obsolete tests, updates goldens intentionally, and runs the
smallest benchmark set that can detect the known failure modes.

## Unit Test Matrix

### Latent Expression Model

- `p_states` has shape `(R, 2)` and rows sum to 1.
- Spliced evidence raises `p_expressed`.
- Strand-RNA lower-bound evidence raises `p_expressed`.
- Intergenic unspliced-only evidence remains mostly `p_unexpressed`.

### Native Support Payload

- One physical fragment spanning both region boundaries increments region support once.
- Two blocks from the same fragment hitting one region increment support once.
- Fractional mass identities still match `n_observed - n_unannotated_ref`.
- Support arrays sort with region arrays.

### Regional gDNA Mass

- `M_r + R_r = T_r` after float32 conversion.
- Strand deconvolution assigns antisense-compatible unspliced mass to gDNA when strand contrast is
  identifiable.
- Boundary projection imputes contained gDNA in unstranded cases.
- No-gDNA sentinel remains near zero.

### EB Exposure

- Uniform density learns low exposure variance and shrinks toward 1.0.
- Capture-like skew learns high exposure variance.
- Same raw ratio with low support shrinks more than high support.
- Extreme high support can produce `omega_r > 1000` without clipping to 1.
- Zero/near-zero mass produces finite exposure.

### Downstream EM

- `gdna_eff_len_em = gdna_eff_len_unweighted / omega_locus`.
- Native EM receives the adjusted denominator.
- Enriched gDNA loci become more competitive, not less.

## Integration Tests

Use focused synthetic cases before full goldens:

1. Non-capture, uniform gDNA, no RNA expression.
2. Non-capture, strand-specific RNA expression plus low gDNA.
3. Unstranded RNA expression plus boundary-crossing gDNA.
4. Hybrid capture locus: 100 kb locus with 2 kb captured subregion.
5. Tiny captured exon with many touching fragments but small fractional mass.
6. Tiny region with one or two touching fragments and a high raw ratio.
7. No-gDNA RNA-only sentinel.

The capture cases should assert relative behavior rather than exact counts:

- on-target `omega_r` above off-target `omega_r`,
- high-support on-target regions shrink less than low-support on-target regions,
- gDNA posterior increases after the denominator sign fix,
- mRNA/nRNA siphoning does not worsen in the no-gDNA sentinel.

## Golden Output Strategy

Golden files should be updated only after the unit and integration tests explain the expected
behavior shift. The output schema will change, so do not preserve old fields for compatibility.

Expected schema removals:

- capture state summaries,
- `gamma_r`,
- `p_captured`,
- `capture_enrichment_target`,
- exposure diagnostics whose sign implies effective-length multiplication.

Expected schema additions:

- regional support summaries,
- `omega` summaries,
- shrinkage summaries,
- `gdna_eff_len_unweighted`,
- `gdna_eff_len_em`,
- `gdna_exposure_factor`,
- `gdna_eff_len_adjustment_ratio`.

## Benchmark Checks

Use the existing benchmarking commands from project instructions when the implementation reaches an
end-to-end state:

```bash
conda activate rigel && python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml
conda activate rigel && python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml
conda activate rigel && python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml -o results/benchmark_report
```

For early iteration, run only targeted conditions first:

```bash
conda activate rigel && python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --conditions gdna_none_ss_1.00_nrna_none gdna_high_ss_0.90_nrna_none
```

## Acceptance Metrics

Track these in the report:

- Distribution of `omega_r` in non-capture conditions: p50 near 1.0 and narrow p95/p99.
- Distribution of `omega_r` in capture conditions: high p99/max without hard clipping.
- gDNA leakage into RNA in captured regions.
- RNA leakage into gDNA in highly expressed regions.
- nRNA siphon magnitude at complex loci.
- Low-support high-ratio regions: shrinkage weight visibly below high-support analogs.
- Runtime and memory change from native support arrays.

## Cleanup Checklist

- Remove obsolete capture-state tests instead of keeping them skipped.
- Remove temporary compatibility aliases for `A_r` if `omega` is the production name.
- Remove Python `BoundaryTable` production path after native boundary payload is consumed.
- Update docs that mention the four-state calibration model.
- Re-run native build after C++ changes.
- Run the smallest passing test set for each PR, then the full suite once the end-to-end path is
  restored.

## Done Means

- The redesigned calibration path is the only production path.
- Tests exercise ordinary, stranded, unstranded, and capture behavior.
- Goldens and summaries reflect the new model rather than compatibility with v6 capture states.
- Benchmark analysis shows whether remaining errors are biological identifiability issues or
  implementation defects.
