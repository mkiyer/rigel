# Beta-Binomial strand-LLR pilot — mini-stress grid

*Date*: 2026-04-20

## TL;DR

The Beta-Binomial shared-κ strand LLR (v3 re-adoption, now a
`CalibrationConfig.strand_llr_mode = "betabinom"` toggle) is
implemented, unit-tested (9/9 pass), and causes **zero regression**
vs the default Binomial mode on the existing 105-test calibration
suite.

On the mini-stress grid, however, the mode switch has **no observable
effect**: all 24 quick-grid cells (ss ∈ {0.5, 0.9, 1.0} × gdna_frac ∈
{0.0, 0.1, 0.5, 2.0} × seed=1) hit the data-driven strand gate and
disable the strand pathway entirely. Binomial and BetaBinom therefore
produce bit-identical λ̂, π, and γ in every cell.

Empirical adoption decision has to come from a benchmark where the
strand channel **actually fires** (real stranded RNA-Seq, or a
synthetic scenario with strand-informative *unspliced* reads). The
mini-stress grid is currently not such a benchmark.

## Harness

- Script: [scripts/benchmark/betabinom_pilot.py](../../scripts/benchmark/betabinom_pilot.py)
- Output: `scripts/benchmark/results/betabinom_pilot_quick/results.tsv`
- Grid (quick): ss ∈ {0.5, 0.9, 1.0} × gdna_frac ∈ {0.0, 0.1, 0.5,
  2.0} × seed ∈ {1} × mode ∈ {binomial, betabinom} = **24 runs**.
- Uses the `build_mini_genome` generator from
  `calibration_sweep.py` (25 genes, 2 Mb, 6 short exons/gene, 50k
  mRNA fragments/cell).

## Results

### λ̂ recovery (λ̂ / λ_truth)

| ss   | gdna_frac | binomial ratio | betabinom ratio | Δ |
|-----:|----------:|---------------:|----------------:|---:|
| 0.50 | 0.10      | 0.9976         | 0.9976          | 0  |
| 0.50 | 0.50      | 0.9999         | 0.9999          | 0  |
| 0.50 | 2.00      | 0.9997         | 0.9997          | 0  |
| 0.90 | 0.10      | 0.9976         | 0.9976          | 0  |
| 0.90 | 0.50      | 0.9999         | 0.9999          | 0  |
| 0.90 | 2.00      | 0.9997         | 0.9997          | 0  |
| 1.00 | 0.10      | 0.9983         | 0.9983          | 0  |
| 1.00 | 0.50      | 0.9995         | 0.9995          | 0  |
| 1.00 | 2.00      | 0.9999         | 0.9999          | 0  |

(gdna_frac = 0.00 omitted; λ_truth = 0 by construction.)

### Diagnostics (identical across modes)

- `strand_used = False` in **all 24 cells**.
- `kappa = 0.0` in **all 24 cells** (never estimated — the BetaBinom
  branch never runs).
- `pi` ∈ [0.907, 1.000], consistent with the density-pathway-only
  identification of a sharp gDNA component.
- `em_n_iter` matches bit-for-bit between modes in all cells.

## Why the strand gate fires (root cause of the null)

The pre-EM aggregate z-test in
[_aggregate_strand_z](../../src/rigel/calibration/_em.py#L167) pools
sense / antisense counts over **unspliced** reads only
(`stats["n_unspliced"]`). In the mini-genome:

- Annotated-spliced mRNA fragments: **44,372** obs (these carry the
  strand signal but are *excluded* from this gate — they go to the
  separate spliced-strand trainer, not the EM strand channel).
- Unspliced fragments: ~6,870 (mostly exonic gDNA contamination plus
  a handful of short-fragment mRNA).
- Intergenic fragments: 94,072 (also excluded — no gene strand).

The unspliced pool is gDNA-dominated, so p̂ ≈ 0.497 (symmetric) and
z ≈ −0.5 for every cell — far below the `strand_z_threshold = 3.0`
cutoff. The gate correctly concludes "unspliced reads carry no
strand signal" and disables the strand channel.

**This is gate behaviour, not BetaBinom behaviour.** The BetaBinom
code never sees a fragment; its effect on calibration is null by
construction whenever the gate is shut.

## Unit-test evidence that the code works

All 9 unit tests in
[tests/test_betabinom_strand.py](../../tests/test_betabinom_strand.py)
pass. Key behavioural checks that rule out a silent no-op in the
implementation itself:

- `test_betabinom_large_kappa_approaches_binomial`: at κ = 1e6 the
  Beta-Binomial LLR matches Binomial LLR within 5%.
- `test_betabinom_smaller_kappa_tempers_llr`: at κ = 5 with
  overdispersed simulated draws (κ_true ≈ 10) the BetaBinom LLR is
  strictly smaller in magnitude than the Binomial LLR on the same
  counts.
- `test_kappa_mle_recovers_true_dispersion_roughly`: κ̂ ∈ (5, 80)
  when true dispersion is κ = 20.
- `test_run_em_betabinom_mode_smoke`: end-to-end EM with mode set to
  `"betabinom"` returns `fit.strand_llr_mode == "betabinom"` and
  `fit.kappa > 0` *when the strand gate opens*.

So: the code is live, correct, and off-by-default wiring is intact.
What the mini-stress grid shows is only that the gate upstream is
too conservative for this benchmark's fragment distribution, not
that BetaBinom is inactive.

## Recommendations

### 1. Keep BetaBinom as an opt-in mode (ship the toggle)

Zero regression risk confirmed. Users who know their library has
strong strand signal in unspliced reads (e.g. real stranded libs
with intron retention or single-exon transcripts) can flip
`strand_llr_mode: betabinom` and benefit from dispersion-aware
tempering. Default stays `binomial`.

### 2. Do not draw adoption conclusions from the mini-stress grid

The mini-stress grid cannot discriminate between the modes because
it doesn't deliver strand signal to the gate. Any future "is
BetaBinom better" decision needs a benchmark where `strand_used =
True` in a non-trivial fraction of runs.

Candidate benchmarks:
- **vcap** (`/scratch/…/rigel_benchmarks/vcap_sim/`): real stranded
  library, known to fire the strand channel.
- A synthetic scenario designed with intron retention or
  single-exon genes so the unspliced pool is strand-informative at
  ss ≥ 0.9.

### 3. Independent follow-up: audit the strand gate

Separate from the BetaBinom question, the observation that the
mini-stress grid **never** fires the strand channel — even at
ss = 1.0 with 44k spliced unique mappers — deserves scrutiny.
Two interpretations, not mutually exclusive:

- *Working as intended*: the EM-time strand channel was designed
  to extract information from unspliced reads specifically
  (because the spliced-strand pathway already trains a separate
  `ss_est`). In that framing, disabling it when unspliced reads
  are uninformative is correct.
- *Over-restrictive pool*: the pooled p̂ averages signal from
  oriented gene-proximal unspliced reads together with all the
  (symmetric) exonic-overlap gDNA unspliced reads, diluting the
  bias. A region-weighted stat might still detect signal where
  the pool does not.

If (2) is true in practice, BetaBinom may help more once the gate
lets signal through. This audit should be scoped as its own task
before re-running the mini-stress pilot.

## Files touched this pilot

- **New**: [scripts/benchmark/betabinom_pilot.py](../../scripts/benchmark/betabinom_pilot.py)
  (92 lines) — pilot harness reusing `build_mini_genome`, running
  each cell in both modes.
- **New**: [scripts/benchmark/results/betabinom_pilot_quick/results.tsv](../../scripts/benchmark/results/betabinom_pilot_quick/results.tsv)
  (24 rows, 39 cols).
- **Unchanged**: all production code paths — implementation
  itself was validated in the previous work unit (9/9 unit tests,
  105/105 calibration-suite pass).

## Status

- ✅ Implementation complete and covered by unit tests.
- ✅ Shipping toggle wired through `CalibrationConfig` →
  `calibrate_gdna` → `run_em` → `EMFit` → `CalibrationResult`.
- ⚠️ Mini-stress pilot: indistinguishable from baseline — strand
  gate closed in all cells.
- ⏳ Adoption decision: **deferred** pending a benchmark that
  actually exercises the strand channel (recommend vcap next).
