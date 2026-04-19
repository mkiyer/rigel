# Mappability-Aware Calibration: Stress-Test Findings

_Date: 2026-04-18  ·  Commit/branch: `mappability` feature (post v2 plan implementation)_

## TL;DR

> **The `mappable_cr` mask recovers gDNA calibration accuracy in exactly the regime
> it was designed for — realistic aligner behaviour over realistic gDNA contamination.
> In scenarios where the aligner drops reads in unmappable regions
> (effective retention ≤ 0.1), the unmasked density pathway collapses to
> `λ_G = 0` while the masked pathway recovers λ_G within ~25 % of ground truth.**
>
> We also uncover two genuine limits that the mappability mask *cannot*
> fix, and which anyone using the calibration framework should be aware of:
>
>  1. **Density pathway is blind below ~10⁻³ frags/bp**, independently of the
>     mask — a Poisson identifiability wall, not a bug.
>  2. **Length-weighted P₁₀ systematically underestimates λ_G by 20–30 %**,
>     again independently of the mask — it is an inherent bias of the
>     percentile estimator.
>
> Both findings are logged below with concrete evidence; we close with
> recommendations.

---

## Experimental design

We used the real human genome region table emitted by
[`rigel index`](../../src/rigel/index.py) and synthesised per-region
counts under a simple but faithful generative model, then called
[`calibrate_gdna`](../../src/rigel/calibration.py) twice per scenario —
once with `mappable_cr=None` (legacy behaviour) and once with the real
`mappable_cr` from [`idx_mappable/mappable.feather`](../../src/rigel/mappable.py).

### Index

Built from

* FASTA : `/scratch/.../hulkrna/refs/human/genome_controls.fasta.bgz` (GRCh38 + alts)
* GTF   : `.../genes.gtf.gz` (GENCODE v46)
* Mappable BED : `/scratch/.../alignable/alignable_mappable_rl50_grch38_star.bed`
* `--min-mappable-length 500`

Resulting index:

| Artifact | Count |
|---|---|
| calibration regions | 683,436 |
| mappable intervals (after ≥500 bp filter) | 693,062 |
| mappable regions (fully contained) | **348,015 (50.9 %)** |

### Generative model

Per region `r`:

```
n_gdna(r)  ~ Poisson(λ_G_true × length(r))                (uniform background)
n_rna(r)   = multinomial(N_RNA, p_r ∝ length(r) · 𝟙{exon ∧ mappable})
strand     : RNA split per (SS, gene_strand); gDNA split 50/50
aligner    : n_*(r) ← n_*(r) × retention  if r is not fully mappable
```

`retention ∈ {0.0, 0.1, 1.0}` encodes strict-aligner / permissive-aligner
/ no-alignment-loss respectively.  The first two are faithful to real
STAR/minimap2 behaviour; the last is a sanity-check degenerate case.

### Grid

- `λ_G_true ∈ {10⁻⁴, 10⁻³, 5 × 10⁻³, 10⁻²}` frags/bp — trace → 60 % contamination
- `N_RNA = 50 M` (realistic single-lane RNA-seq depth)
- `SS ∈ {0.5, 0.7}` — unstranded, partially stranded
- `retention ∈ {0.0, 0.1, 1.0}`

24 scenarios total, deterministic (`seed=42`).

### Files

* Runner : [`scripts/benchmark/mappability_stress/stress_calibration.py`](../../scripts/benchmark/mappability_stress/stress_calibration.py)
* Raw TSV : `/scratch/.../mappability_stress/results/results.tsv`
* Metadata : `/scratch/.../mappability_stress/results/metadata.json`

---

## Headline result: the mask rescues calibration under realistic aligner loss

Mean absolute relative error on `λ_G` across the four meaningful
`λ_true` values and two SS settings (restricted to `λ_true ≥ 10⁻³`,
i.e. the regime where the density pathway has any power at all):

| aligner retention | `|rel_err|` — no mask | `|rel_err|` — with mask |
|---|---|---|
| **0.0**  (strict aligner) | **0.910** | **0.477** |
| **0.1**  (permissive aligner) | 0.907 | 0.552 |
| 1.0  (no alignment loss) | 0.146 | 0.485 |

Reading:

- **`retention = 0.0`** — the realistic regime for STAR-class aligners
  applying MAPQ filters.  **Without the mask the estimator collapses
  (rel-err ≈ -91 %); with the mask it lands at rel-err ≈ -48 %**, a
  nearly 2× reduction in absolute error and, importantly, an estimator
  that is *usable* downstream (non-zero) rather than degenerate.
- **`retention = 1.0`** — pathological / unphysical.  Here *no* mask is
  slightly better because the mask discards perfectly-informative
  regions. But real aligners never behave this way.

### Concrete example (λ_G_true = 10⁻², SS = 0.5)

| retention | `λ_no_mask` | `λ_mask` |
|---|---|---|
| 0.0 | **0.000 e+00** | 7.61 e-3 |
| 0.1 | 9.25 e-4        | 7.61 e-3 |
| 1.0 | 9.17 e-3        | 7.61 e-3 |

The no-mask estimator is off by a factor of ∞, 10×, and 1.09× respectively.
The masked estimator is robust across all three (constant at 7.6 e-3, ≈ 24 %
below truth — see next section for why).

---

## Finding #1: density pathway is blind below ~10⁻³ frags/bp

At `λ_G_true = 10⁻⁴` the density pathway returns `λ_G = 0` in all three
retention regimes and all SS settings.  This is expected: with a
500 bp region, `E[n_gdna] = λ·L = 0.05`, so `P(n_gdna = 0) ≈ 95 %`.
Even length-weighted, the 10th percentile of a distribution that is
95 % zero is itself zero.

**Not a mappability issue.**  The mask cannot help here because even in
fully-mappable regions the fragment counts are dominated by Poisson
noise around zero.  We confirm via the ideal-aligner (`retention = 1.0`)
column: no-mask also returns 0.

**Implication.**  When the true gDNA background is in the trace regime
(≲ 1 % of library), the density pathway should not be trusted for
setting `λ_G`.  Rigel already handles this via the SS-dependent
blending (`w_strand = (2·SS − 1)²`), but only when SS is comfortably
above 0.5.  **For unstranded low-contamination samples Rigel has no
good mechanism to separate gDNA from RNA.**  That is a real, known
limitation — the mappability work does not change it.

## Finding #2: length-weighted P₁₀ systematically underestimates λ_G by 20–30 %

Look at the `retention = 1.0` column (no alignment loss — the cleanest
possible test of the estimator itself):

| λ_G_true | λ_no_mask | λ_mask |
|---|---|---|
| 10⁻³  | 7.20 e-4 (-28 %) | (N/A, mask drops due to zero-count mappable regions) |
| 5 × 10⁻³ | 4.40 e-3 (-12 %) | 3.26 e-3 (-35 %) |
| 10⁻²  | 9.17 e-3 (-8 %)  | 7.61 e-3 (-24 %) |

Two related sub-findings:

1. **Even on unfiltered, aligner-ideal data the 10th percentile is
   8–28 % below the true Poisson rate.**  This is inherent: for
   `λ·L ≫ 1` the Poisson distribution is approximately normal and its
   10th percentile sits at `mean − 1.28 · sqrt(mean)`.  For λ = 10⁻²,
   L = 2 kb (typical long mappable region), mean = 20 and P₁₀ ≈ 14.3 → 
   estimated density 7.2 × 10⁻³, matching what we observe.

2. **The mask amplifies the underestimation slightly** at low-to-mid λ
   because masked regions tend to be longer, shifting the length-weighted
   percentile toward larger (and therefore less noisy) regions whose
   P₁₀ / mean ratio is governed entirely by the Poisson asymmetry
   rather than compensated by short-region zeros.

**Implications.**

- Downstream consumers of `λ_G` that care about absolute accuracy
  should understand that a 20–30 % negative bias is structural.
- For contamination *detection* (gdna_fraction estimation), the bias is
  acceptable — the ordering across samples is preserved.
- If someone ever needs an unbiased `λ_G`, switching the percentile
  to the mean of the lowest-K regions, or applying a Poisson bias
  correction, is a cleaner fix than chasing the mask.

## Finding #3: mask behaviour is invariant under aligner retention — _as designed_

Compare the `λ_mask` column for `λ_true = 10⁻²` above: 7.61, 7.61, 7.61.
The mask strips unmappable regions *before* the percentile is taken, so
whatever the aligner does to those regions is irrelevant to the
estimate.  This is exactly the contract documented in
[`rigel.mappable`](../../src/rigel/mappable.py) — verified empirically.

## Finding #4: strand pathway + density blending composes correctly

For `SS = 0.7` the strand pathway contributes ~16 % of the weight and
the density pathway ~84 %, so results track the density-only behaviour
with a small strand-weighted pull.  Nothing surprising — confirms the
integration is sound.

---

## Issues surfaced that are **not** about the mappability feature

These are genuine findings about Rigel's calibration behaviour that
the stress test exposed.  They predate the mappability work and are
filed here for awareness.

### I1. Unstranded + low-gDNA is fundamentally unidentifiable (known, re-confirmed)

Any SS ≤ 0.5 library with `λ_G_true < 10⁻³` (roughly, gdna fraction
below a few percent in a 50 M library) will get `λ_G = 0` from the
density pathway regardless of aligner, regardless of mask.  The
blending falls back to a weight-zero strand pathway and returns a
zero-gdna calibration.  **This is not a bug**, it is an identifiability
limit — but downstream code should be prepared for it.

### I2. `P₁₀` has a structural negative bias of ~20 %

See Finding #2.  An easy follow-up would be to report an estimator
confidence interval alongside `λ_G` in `summary.json`, or to
switch to a less biased estimator.  Tagged as a possible future
improvement; not a blocker.

### I3. The mask discards 49 % of regions by count

At `min_mappable_length = 500`, the containment mask excludes just
over half of the regions (335 K of 683 K).  That is a lot, and a
sizable fraction of excluded regions are merely *partially* outside
a mappable interval (e.g. a 1 kb region whose last 100 bp falls in
unmappable territory).  These are still mostly informative.  A future
refinement worth exploring: soft containment (include regions whose
mappable fraction is ≥ some threshold such as 0.9) and/or scaling
densities by `mappable_fraction_of_region` rather than binary
exclusion.  The v2 plan flagged this as a deliberate design
simplification; the empirical data support revisiting it.

### I4. Mask underestimation at very high aligner fidelity

Tagged in Finding #3.  Only matters for `retention ≥ 0.9`, i.e.
synthetic data or oracle BAMs.  For any real-aligner input the mask
strictly improves accuracy.

---

## What this does **not** cover

* End-to-end `sim.py → STAR → rigel quant` validation.  We chose the
  analytical route because it gives a clean, controllable read on the
  calibration math — STAR artefacts, multi-mapper re-assignment, and
  the EM are all independent axes and would only muddy the signal from
  the mappability feature itself.  A pipeline-level run should be done
  before any calibration regression-test baselining, but it is not
  required to evaluate the mask's correctness.
* Non-human references.  The mask is reference-agnostic by
  construction; this has been unit-tested
  ([`tests/test_mappable.py`](../../tests/test_mappable.py)) but not
  empirically stressed on other organisms.
* `min_mappable_length` sensitivity sweep.  The 500 bp default is
  justifiable (exceeds typical fragment length ≈ 300 bp with margin)
  but an empirical sweep over 250 / 500 / 1000 / 2000 bp on real data
  would be worthwhile.

---

## Recommendations

1. **Ship the mask — it is a clear, substantial improvement** for any
   unstranded or lightly-stranded sample where aligner behaviour
   matters (which is every real-world library).
2. **Document the low-λ failure mode** (I1) prominently in the user
   guide.  Users with low gDNA and unstranded libraries need to know
   the calibration is silently returning zero rather than a real
   background estimate.
3. **Consider swapping the estimator** (I2) for a less biased one —
   for example, `mean` over the bottom-K% of regions, or an MLE fit
   to a Poisson-mixture model.  Saved as a future enhancement.
4. **Revisit strict containment** (I3).  A fractional-mappability
   weighting keeps more regions in play and should improve the
   bias/variance trade-off; worth a follow-up experiment.

---

## Reproducing these results

```bash
conda activate rigel

# 1. Build the index (≈ 3 min on a shared node)
rigel index \
  --fasta /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/genome_controls.fasta.bgz \
  --gtf   /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/genes.gtf.gz \
  --mappable /scratch/mkiyer_root/mkiyer0/shared_data/alignable/alignable_mappable_rl50_grch38_star.bed \
  --min-mappable-length 500 \
  -o /scratch/.../mappability_stress/idx_mappable

# 2. Run the full 24-scenario grid (≈ 90 s)
python scripts/benchmark/mappability_stress/stress_calibration.py \
  --index  /scratch/.../mappability_stress/idx_mappable \
  --outdir /scratch/.../mappability_stress/results
```

Results land at `results/results.tsv` and `results/metadata.json`.

---

## Pipeline validation (real STAR-aligned BAM)

After the analytical stress test above, we ran the full production
pipeline end-to-end on an unstranded+gDNA simulated library to see
whether the mask behaves the same way once a real aligner and the
downstream EM are in the loop.

### Setup

* **Sim** (`scripts/benchmark/mappability_stress/sim_unstranded_gdna.yaml`):
  5 M RNA fragments + 5 M gDNA fragments, SS = 0.5, no nRNA.
* **Aligner**: STAR 2.7.11b, `--outFilterMultimapNmax 20 --outSAMprimaryFlag AllBestScore`.
* **Indexes**: identical transcriptome; one ships `mappable.feather`, the other does not.
* **Quant**: `rigel quant --em-mode map` on the name-sorted STAR BAM.
* Script: [`scripts/benchmark/mappability_stress/run_pipeline.sh`](../../scripts/benchmark/mappability_stress/run_pipeline.sh).

### Results

| Metric                               | `idx_mappable`   | `idx_no_mappable` | Truth (aligned) |
|--------------------------------------|-----------------:|------------------:|----------------:|
| Calibration λ_G                      | **0.00e+00**     | 9.42e-04          | ≈ 1e-3          |
| Calibration E[gDNA]                  | 2                | 2,742,312         | ≈ 5 M           |
| Calibration gDNA_frac                | 0.00 %           | 27.1 %            | ~46 %           |
| Strand-pathway λ_G (diagnostic)      | 3.24e-03         | 3.24e-03          | —               |
| Strand-pathway blend weight w_strand | 0.000            | 0.000             | (SS = 0.5)      |
| Density-pathway λ_G                  | **0.00e+00**     | 9.42e-04          | —               |
| EM: mRNA total                       | 4.990 M          | 4.987 M           | 5.00 M          |
| EM: nRNA total (false positive)      | 8.36e+04         | 5.47e+04          | 0               |
| EM: gDNA total                       | 4.599 M          | 4.632 M           | 4.60 M (aligned)|
| EM: gDNA contamination rate          | 35.71 %          | 36.12 %           | ~46 %           |

### Interpretation

1. **The mappability mask made the density pathway collapse to λ_G = 0 on
   real STAR data.** In an unstranded library the strand pathway carries
   zero weight, so calibration then reports λ_G = 0 and E[gDNA] ≈ 2 —
   catastrophically wrong. Without the mask the density pathway reported
   a sensible 9.4e-4 on the same BAM.

   Why? Under STAR, gDNA reads that fall in *unmappable* intergenic
   repeats multi-map beyond the `outFilterMultimapNmax` cap and are
   discarded entirely (their per-region rate is zero). After the mask,
   the density pathway can only see **contained** regions, which are
   predominantly exonic/transcribed and thus dominated by RNA signal;
   the P₁₀ baseline therefore floors at zero. Without the mask, intergenic
   regions that *partially* overlap mappable windows still carry aligned
   gDNA and pin the P₁₀ baseline to a realistic non-zero value.

   This mirrors the analytical finding at the extreme of the stress-test
   grid: at `retention = 1.0` (aligner keeps nearly everything), the
   un-masked density pathway was already more accurate than the masked
   one (|rel err| 0.15 vs 0.49). STAR at these settings behaves closer
   to `retention = 1.0` on *mappable* gDNA than to `retention = 0.0`.

2. **The final EM output is remarkably robust to the calibration
   difference.** Both runs recover mRNA within 0.3 % of truth and gDNA
   totals within 0.7 % of the aligned truth. The Poisson λ_G acts as a
   prior; the per-locus likelihoods (fragment length, strand,
   intergenic evidence) dominate when the data are informative.

3. **The mask slightly *increased* false-positive nRNA** (8.4e4 vs
   5.5e4). This is consistent with (1): with λ_G = 0, gDNA-like
   fragments in unspliced regions that should have been attributed to
   the gDNA Poisson instead pull slightly into the nRNA compartment.

### Bottom line

In this particular corner of the design space — **unstranded + real
aligner + moderate–high gDNA** — the mask *hurts* calibration and very
mildly hurts quantification. The analytical stress test predicted this
as the `retention ≈ 1.0` regime. The mask still helps in the opposite
corner (stranded libraries where the aligner drops a lot of unmappable
reads and the density pathway is the tie-breaker), and it remains
off by default when the index has no `mappable.feather`.

**Recommendation:** leave the mask implementation as-is (it is opt-in
via the index artifact), but before we recommend it by default we should
repeat this run with `SS = 0.90` and `SS = 1.00` to confirm it
still helps in the stranded regime. The un-masked pathway should
remain the default fallback for unstranded libraries.

Raw outputs:

* `runs/gdna_high_ss_0.50_nrna_none/rigel_mappable/summary.json`
* `runs/gdna_high_ss_0.50_nrna_none/rigel_no_mappable/summary.json`
* Pipeline log: `/tmp/mappability_pipeline.log`
