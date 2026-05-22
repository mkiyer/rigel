# RNA/gDNA FL Estimation Stress Test

Date: 2026-05-22
Script: `scripts/sim/fl_estimation_stress.py`
Primary output: `results/fl_estimation_stress_20260522/`

## Goal

Stress-test the RNA and gDNA fragment-length (FL) estimation path that is already
available before Phase 4 EM restoration:

1. Generate mini-genomes with arbitrary multi-exon transcripts.
2. Simulate oracle BAMs with independent RNA and gDNA FL distributions.
3. Run Rigel's native scan and calibration path.
4. Compare fitted RNA/gDNA FL models with true sampled fragment lengths encoded
   in oracle BAM read names.

The script intentionally stops after scan + calibration. `rigel quant` still
stops at `FractionalCutoverPending` until Phase 4 restores per-locus priors and
EM, but the FL path is fully exercised.

## Commands Run

Production-threshold five-case matrix:

```bash
conda activate rigel
python scripts/sim/fl_estimation_stress.py \
  --n-rna 6000 \
  --gdna-fraction 1.0 \
  --genome-length 80000 \
  --n-genes 12 \
  --threads 2 \
  --out-dir results/fl_estimation_stress_20260522
```

High-support single-case check:

```bash
python scripts/sim/fl_estimation_stress.py \
  --case rna150_gdna70 \
  --n-rna 60000 \
  --gdna-fraction 1.0 \
  --genome-length 80000 \
  --n-genes 12 \
  --threads 2 \
  --out-dir results/fl_estimation_stress_high_support_20260522
```

Low-support default-threshold check:

```bash
python scripts/sim/fl_estimation_stress.py \
  --case rna150_gdna70 \
  --n-rna 800 \
  --gdna-fraction 1.0 \
  --genome-length 30000 \
  --n-genes 6 \
  --threads 1 \
  --out-dir results/fl_estimation_stress_low_support_20260522
```

Low-support good-branch smoke check:

```bash
python scripts/sim/fl_estimation_stress.py \
  --case rna150_gdna70 \
  --n-rna 800 \
  --gdna-fraction 1.0 \
  --genome-length 30000 \
  --n-genes 6 \
  --threads 1 \
  --pool-quality-good 500 \
  --pool-quality-weak 50 \
  --out-dir results/fl_estimation_stress_smoke_20260522_good
```

## Main Matrix Results

All five production-threshold cases had `rna_quality=good` and
`gdna_quality=good`.

| Case | Truth RNA mean | Fitted RNA mean | RNA error | Truth gDNA mean | Fitted gDNA mean | gDNA error | RNA recovery | gDNA pool recovery |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `rna150_gdna70` | 149.45 | 149.45 | +0.00 | 70.86 | 70.80 | -0.05 | 1.000 | 0.937 |
| `rna70_gdna150` | 70.59 | 71.50 | +0.91 | 149.39 | 149.40 | +0.01 | 0.971 | 0.937 |
| `rna150_gdna150` | 149.89 | 149.89 | +0.00 | 149.66 | 149.66 | +0.00 | 1.000 | 0.938 |
| `rna250_gdna80` | 249.65 | 249.65 | +0.00 | 79.63 | 79.72 | +0.09 | 1.000 | 0.936 |
| `rna80_gdna250` | 80.08 | 80.53 | +0.45 | 250.48 | 249.88 | -0.60 | 0.989 | 0.923 |

Maximum absolute raw fitted mean error:

- RNA: 0.91 bp.
- gDNA: 0.60 bp.

Interpretation: the scanner/calibration FL histograms are accurate when the RNA
and gDNA channels have enough support to use the `good` branch.

## Original Caveat: Scoring PMF Add-One Smoothing

`FragmentLengthModel.mean` reports the raw-count mean in the `good` branch, but
`FragmentLengthModel.log_likelihood()` and `.pmf` use add-one smoothing across
all `max_frag_length + 1` bins. With `max_frag_length=600` and 6000
observations, this adds 601 uniform pseudocounts, which is not negligible for
short/narrow distributions.

Examples from the main matrix:

| Case | Raw fitted RNA mean | Scoring-PMF RNA mean | Raw fitted gDNA mean | Scoring-PMF gDNA mean |
| --- | ---: | ---: | ---: | ---: |
| `rna150_gdna70` | 149.45 | 163.15 | 70.80 | 92.93 |
| `rna70_gdna150` | 71.50 | 92.87 | 149.40 | 163.94 |
| `rna80_gdna250` | 80.53 | 100.72 | 249.88 | 254.78 |

At 60k observations in `rna150_gdna70`, the distortion is much smaller:

- RNA raw fitted mean: 149.47; scoring-PMF mean: 150.97.
- gDNA raw fitted mean: 70.64; scoring-PMF mean: 73.06.

This was meaningful for synthetic mini-runs and small calibration pools. It was
addressed after this stress test by replacing add-one-per-bin smoothing with one
total pseudo-observation spread across the full FL support. Finalized model
summary statistics now use the same predictive PMF used by scoring when a model
has observations or an EB prior.

## Original Low-Support Behavior

With the old production thresholds and only 800 RNA + 800 gDNA fragments, both
FL models were classified `weak` and shrank strongly toward the global mixed FL
histogram:

| Case | RNA quality | gDNA quality | Truth RNA mean | Fitted RNA mean | Truth gDNA mean | Fitted gDNA mean |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `rna150_gdna70` | weak | weak | 148.77 | 170.24 | 70.09 | 146.18 |

This was expected from the old EB policy, not a scanner bug. It showed that the
weak branch could not preserve sharply separated RNA/gDNA FL distributions at
very low support.

When the same low-support run is forced into the `good` branch by lowering
thresholds, the fitted means are accurate:

| Case | Truth RNA mean | Fitted RNA mean | Truth gDNA mean | Fitted gDNA mean |
| --- | ---: | ---: | ---: | ---: |
| `rna150_gdna70` | 148.77 | 148.77 | 70.09 | 70.11 |

## Policy Update Check

After the stress test, the weak-pool policy was changed so non-empty pools are
modeled by default and global shrinkage is capped logarithmically from the pool's
own evidence. A low-support rerun used the default adaptive weak branch:

```bash
python scripts/sim/fl_estimation_stress.py \
  --case rna150_gdna70 \
  --n-rna 800 \
  --gdna-fraction 1.0 \
  --genome-length 80000 \
  --n-genes 12 \
  --threads 2 \
  --out-dir results/fl_estimation_stress_after_fl_policy_20260522_low_support
```

| Case | RNA quality | gDNA quality | Truth RNA mean | Fitted RNA mean | Truth gDNA mean | Fitted gDNA mean | Scoring RNA mean | Scoring gDNA mean |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `rna150_gdna70` | weak | weak | 148.77 | 148.63 | 70.09 | 70.79 | 148.63 | 70.79 |

## Conclusions

1. The native scan + fractional gDNA FL pool accumulation + `build_fl_models`
   path is accurate under oracle conditions when channel-specific evidence is
   sufficient.
2. gDNA FL pool recovery is about 0.92-0.94 in these mini-genomes because exon
   and other non-gDNA-compatible gDNA fragments are intentionally excluded from
   the gDNA FL training pool. The retained pool is still representative enough
   to recover the mean accurately.
3. The original low-support EB shrinkage was too strong for divergent RNA/gDNA
  FLs. The updated adaptive policy preserves low-support channel-specific FL
  evidence while still borrowing a small amount of global shape.
4. The original add-one scoring PMF smoothing was too strong. The updated
  finalizer uses a one-pseudo-observation unseen-bin reserve, which keeps
  likelihoods finite without pulling short distributions toward the histogram
  midpoint.

## Suggested Follow-Ups

- Add a lightweight pytest or benchmark assertion around the good-branch stress
  matrix: max raw mean error < 1 bp under oracle conditions.
- Keep the new low-support weak-branch test in `tests/test_fl_models.py` as a
  guard against reintroducing global over-shrinkage.
- Re-run the full five-case matrix after Phase 4 EM restoration so FL scoring
  sensitivity is evaluated in the restored end-to-end inference path.
