# Strand LLR: Binomial vs Beta-Binomial (mini-stress pilot)

**Date:** 2026-04-20
**Status:** Pilot complete — **Binomial retained as default**
**Toggle:** `CalibrationConfig.strand_llr_mode ∈ {"binomial" (default), "betabinom"}`

---

## Question

Does the v3-style shared-κ Beta-Binomial strand likelihood meaningfully
improve calibration accuracy vs the currently-shipping plain Binomial
when the strand channel is *actually active*?  (The original v3
implementation is preserved in `docs/calibration/calibration_v3_vs_v4_comparison.md`.)

## What was changed for this pilot

**Mini-genome redesigned** to let the aggregate-strand-z gate
actually open.  The previous 2 Mb / 150 bp-exon generator produced
near-zero *unspliced-annotated* RNA reads (fragments ≥ 200 bp always
spliced out), so the gate pooled gDNA-dominated counts and always
rejected the channel.

[scripts/benchmark/betabinom_pilot.py](../../scripts/benchmark/betabinom_pilot.py)
`build_realistic_mini_genome()`:

| Parameter         | Old    | New      |
|-------------------|--------|----------|
| genome size       | 2 Mb   | **5 Mb** |
| gene_pitch        | 80 kb  | **200 kb** |
| exon length       | 150 bp | **2 000 bp** |
| exons per tx      | 6      | 4        |
| exon footprint    | ~3.4%  | **~2.8%** (closer to real genomes ≈ 2 %) |

Exons 10× larger than `frag_mean=200` so unspliced RNA is plentiful;
intergenic is ≥ 95 % of the genome so gDNA does not swamp the
gene-strand bucket.  The strand gate now correctly fires at ss ≥ 0.9
and correctly stays silent at ss = 0.5.

## Grid

3 SS × 4 gDNA × 1 seed × 2 modes = **24 runs**,
50 000 RNA fragments per cell, gDNA fragments = `gdna_fraction` × RNA.

- `ss ∈ {0.5, 0.9, 1.0}`
- `gdna_fraction ∈ {0.0, 0.1, 0.5, 2.0}`
- `mode ∈ {binomial, betabinom}`

Output: `scripts/benchmark/results/betabinom_pilot_quick/results.tsv`.

## Results

### Strand gate behaviour (data-driven)

| ss  | gate fires? | z-score   | p̂_sense |
|-----|-------------|-----------|---------|
| 0.5 | **OFF** ✓   | 0.02–0.17 | 0.500   |
| 0.9 | **ON**      | 164–174   | ~0.90   |
| 1.0 | **ON**      | 205–216   | ~1.00   |

The z-test is doing exactly what it should — zero bias → off, strong
bias → on.

### λ̂ / λ_truth (both modes identical to ≥ 4 decimal places)

| ss  | gdna_frac | binomial | betabinom | κ (betabinom) |
|-----|-----------|----------|-----------|---------------|
| 0.5 | 0.0       | —        | —         | 0 (gate off)  |
| 0.5 | 0.1       | 0.9971   | 0.9971    | 0 (gate off)  |
| 0.5 | 0.5       | 1.0006   | 1.0006    | 0 (gate off)  |
| 0.5 | 2.0       | 0.9996   | 0.9996    | 0 (gate off)  |
| 0.9 | 0.0       | —        | —         | 500 (ceiling) |
| 0.9 | 0.1       | 0.9971   | 0.9971    | 500 (ceiling) |
| 0.9 | 0.5       | 1.0006   | 1.0006    | 500 (ceiling) |
| 0.9 | 2.0       | 0.9996   | 0.9996    | **116**       |
| 1.0 | 0.0       | —        | —         | 0.01 (floor)  |
| 1.0 | 0.1       | 1.0004   | 1.0004    | 500 (ceiling) |
| 1.0 | 0.5       | 1.0001   | 1.0001    | **115**       |
| 1.0 | 2.0       | 0.9996   | 0.9996    | **22**        |

Pairwise λ̂ agreement across all 24 cells is within 10⁻⁴ relative.

### κ behaviour

The MLE sweep returns a finite, sensible κ whenever there is
meaningful per-region overdispersion evidence (22, 115, 116 at
gdna_frac ≥ 0.5 and ss ≥ 0.9).  Two boundary failure modes:

- **Ceiling hits (κ → 500)** in 4/12 cells at low gDNA.  Here the
  binomial and beta-binomial are empirically indistinguishable, so the
  likelihood is monotone in κ and saturates against the bracket.
  Behaviour is statistically correct but suggests the upper bound
  should be widened or the optimiser should detect boundary-stuck
  fits.
- **Floor hit (κ ≈ 0.01)** in 1/12 cells at ss = 1.0 with *no* gDNA
  truth, where both components share the same data and κ is
  unidentifiable.  Degenerate.

## Interpretation

**On this benchmark the count channel dominates** — density alone
already recovers λ̂ to within 10⁻³, leaving no headroom for the strand
channel (or either of its two parametrisations) to differentiate
itself.  The synthetic grid is therefore not a discriminating test
between Binomial and Beta-Binomial.

The correct discriminator is real RNA-seq where:

- count-channel heterogeneity is large (biological expression
  variance, mappability artefacts, UTRs, highly-expressed loci),
- per-region strand observations are overdispersed (biological
  antisense, chimeric pairs, alignment quirks),
- κ → ∞ would over-weight per-region evidence in a way that the
  Binomial cannot temper.

That test is the Armis2 vcap pipeline, not this synthetic sweep.

## Decision

**Keep `binomial` as the shipping default.**  The Beta-Binomial
parametrisation introduces two extra tunables (κ-floor, κ-ceiling) and
a golden-section solver inside the EM loop for zero measurable gain on
the synthetic benchmark.  Simplicity wins in the absence of empirical
evidence for the additional complexity.

The toggle remains available as `strand_llr_mode="betabinom"` for
future head-to-head tests on real data.  Unit tests in
[tests/test_betabinom_strand.py](../../tests/test_betabinom_strand.py)
continue to guarantee:

- BetaBinom → Binomial as κ → ∞.
- LLR ≡ 0 at ss = 0.5.
- κ MLE recovers the true dispersion on overdispersed synthetic data
  to within a factor of ~4.
- `strand_llr_mode="binomial"` is bit-identical to pre-toggle
  behaviour (regression check).

## Related docs

- [Why the mini-genome needed redesigning](strand_llr_noise_floor_vcap_results.md)
- [v3 vs v4 structural comparison](calibration_v3_vs_v4_comparison.md)
- [Strand LLR noise-floor design (shipped layer (1))](strand_llr_noise_floor_design.md)
