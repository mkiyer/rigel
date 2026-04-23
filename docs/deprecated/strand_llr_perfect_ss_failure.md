# Calibration v4 — Strand LLR Failure at Perfect Strand Specificity

**Date:** 2026-04-19
**Author:** analysis run from `scripts/calibration/stress_mini/`
**Status:** diagnostic report; no code changes applied.

---

## 1. One-paragraph summary

The v4 calibration mixture is **accurate to within ~9 %** across a 4×5
(ss × gdna_fraction) sweep on a controlled mini genome — **except** at
`strand_specificity = 1.0`, where `λ_G` is over-estimated by **1.45× – 2.68×**
depending on contamination level.  The failure mode is not in the mixture
structure but in one line of the strand log-likelihood-ratio: a symmetric clip
`ss_c = clip(ss, 1e-12, 1 − 1e-12)`.  At perfectly stranded input this maps
"one antisense read" to **+26.9 nats** of gDNA evidence per fragment —
~39× larger than the corresponding sense penalty of −0.69 nats.  A single
stray antisense observation in a low-coverage genic region is therefore
enough to flip γ from ≈0 to ≈1 and hijack the mixture.  The real vcap
simulations also all report `ss_est ≈ 1.0`, so the same pathology is
almost certainly driving the −11 % gDNA bias previously observed there.

---

## 2. The experiment

All artefacts under `/tmp/mini_stress/` (harness lives in
`scripts/calibration/stress_mini/`).

### Genome
A 1.27 Mb **mini genome** stitched from 500-bp blocks:

| kind  |  count | total bp | mappability |
|-------|-------:|---------:|:-----------:|
| U-tx  |     20 | 1.20 Mb  | 1.0 (one transcript per block) |
| U-int |      8 | 4 kb     | 1.0         |
| R2    |      4 | 2 kb     | 0.50        |
| R5    |     10 | 5 kb     | 0.20        |
| R10   |     20 | 10 kb    | 0.10        |
| R50   |    100 | 50 kb    | 0.02        |

Each U-tx block hosts **one** transcript with 1–10 exons (100–1000 bp each) and
500–10 000 bp introns.  Repeats are perfect 500-bp copies of the same motif so
reads multimap uniformly at N×.  STAR index + `alignable compute` (STAR,
`star_arriba_parameters.txt`) → rigel index with mappability.

### Sweep
For each `(ss, gdna_fraction)` in `{0.5, 0.75, 0.9, 1.0} × {0, 0.5, 1.0, 1.5, 2.0}`:

1. `OracleBamSimulator.write_bam(pool_split=(50 000, 0, int(50 000 · gdna_frac)))`
2. `scan_and_buffer` → trained strand/FL models
3. `calibrate_gdna` → `λ_est`, `π`, `π_soft`, γ, μ_R, σ_R

Truth: `λ_true = n_gdna / L_genome` (uniform scatter).

### Results

`λ_est / λ_true` grid:

| ss \ gdna_frac | 0.5  | 1.0  | 1.5  | 2.0  |
|---:|:---:|:---:|:---:|:---:|
| 0.50 | 1.08 | 1.09 | 1.08 | 1.08 |
| 0.75 | 1.08 | 1.09 | 1.08 | 1.08 |
| 0.90 | 1.08 | 1.09 | 1.08 | 1.08 |
| **1.00** | **2.68** | **1.86** | **1.59** | **1.45** |

All mixture parameters are healthy at ss ≤ 0.9 (π≈0.09, π_soft≈0.17, σ_R≈0.62).
At ss=1.0 they flip: **π_soft=0.99**, σ_R=0.79–0.81, γ distribution becomes
bimodal (99 regions γ<0.01 vs **140 regions γ≥0.99**).

---

## 3. Root cause — one line of `_strand_llr`

### What it should be doing

Under the mixture, each region has a log-odds that the reads there are gDNA.
One piece of evidence is the **sense/antisense imbalance**: under the RNA
model the fraction of sense reads is *p* = `strand_specificity`; under the
gDNA model it is *p* = 0.5.  The per-region strand LLR is

$$
\ell_i = \sum_{j \in \text{region } i} \log\frac{q_G(\mathrm{strand}_j)}{q_R(\mathrm{strand}_j)}
        = k_{\text{sense}} \log\!\frac{0.5}{ss} + k_{\text{anti}} \log\!\frac{0.5}{1-ss}
$$

This is **the right idea**: perfectly stranded libraries genuinely do give
the cleanest separation between RNA (p=ss≈1) and gDNA (p=0.5).  The issue
is calibration numerics, not the information structure.

### What it is actually doing

```python
# src/rigel/calibration/_em.py::_strand_llr  (line 219)
ss_c = float(np.clip(strand_specificity, _EPS, 1.0 - _EPS))     # _EPS = 1e-12
...
log_ratio_sense = math.log(0.5 / ss_c)
log_ratio_anti  = math.log(0.5 / (1.0 - ss_c))
```

The **symmetric 1 × 10⁻¹² floor** was intended as a numerical safety net.
At `ss = 1.0 − 10⁻¹²` it yields:

| term              | value     | interpretation                             |
|-------------------|-----------|--------------------------------------------|
| `log_ratio_sense` | −0.693    | "a sense read is 2× less likely under G"   |
| `log_ratio_anti`  | **+26.94**| "an antisense read is **5 × 10¹¹×** more likely under G" |
| asymmetry         | **39×**   | one antisense ≈ 39 sense reads             |

For comparison, at **`ss = 0.99`** (still essentially perfect):

| term              | value     |
|-------------------|-----------|
| `log_ratio_sense` | −0.683    |
| `log_ratio_anti`  | +3.912    |
| asymmetry         | **5.7×**  |

Concretely, a mostly-mRNA region with exactly **one** stray antisense read:

| n_sense | LLR @ ss=1−10⁻¹² | LLR @ ss=0.99 |
|--------:|----------------:|--------------:|
|       5 | **+23.5** (→G)  |  +0.50        |
|      20 | **+13.1** (→G)  |  −9.75        |
|     100 |  −42.4          | −64.4         |

So at ss=1.0 any RNA region with fewer than ~40 sense reads is flipped to
gDNA by a single antisense observation.  In the mini-genome sweep, most
genic regions are short exons with 10–30 reads each — **guaranteed to flip**
when the simulator injects any gDNA into the mix (gDNA is 50/50 strand, so
every small exon picks up a few antisense gDNA reads).

### Why the mixture then collapses

The M-step for π updates from the posterior γ:

$$
\pi_{\text{new}} = \frac{1}{N} \sum_i \gamma_i
$$

With 140 of 239 eligible regions pushed to γ≈1 by the runaway strand LLR,
π_soft shoots to 0.99.  The μ_R / σ_R anchor then refits against the
surviving 99 "clean" regions — which are large high-coverage exons where
the +26.9 nats per antisense read is finally outvoted.  This inflates
σ_R (more RNA variance "explained" by the few heavily-sampled survivors)
and biases μ_R downward, which in turn *further* widens the pool of
regions that look more gDNA-like than RNA-like.

Empirically this cascade is visible as a **σ_R increase from 0.62 → 0.79**
at ss=1.0 in the sweep, and σ_R=1.64–1.91 in the vcap oracle runs.

### Link to vcap

All three vcap conditions (pristine / gdna / gdna_nrna) report
`ss_est ∈ {1.0, 0.999998, 0.999995, 0.999991}` and π_soft ∈ [0.92, 0.99],
σ_R ∈ [1.64, 1.91] — the same pathological signature.  The −11.3 %
oracle gDNA bias we tracked earlier is therefore **not** a mystery of
deep human-genome biology; it is the same single-line numerical issue,
amplified by the more complex loci.

---

## 4. Why the symmetric clip was the wrong choice

There are **two different concerns** the `1-ε` clamp tries to handle with
one knob:

1. **Numerical**: keep `log(0)` out of the arithmetic.
2. **Modelling**: express disbelief that RNA can produce an antisense read.

The current code collapses these.  `_EPS = 1e-12` is appropriate for (1)
but absurd for (2): it asserts the probability of an antisense read under
the RNA model is 10⁻¹², i.e. once per 10¹² aligned reads.  Real libraries
(even PCR-free, top-scoring stranded kits) show **≥0.5 %** antisense from
readthrough, intron-containing transcripts, convergent transcription,
and simple alignment noise.  The clamp should reflect that floor.

The symmetric clamp is also wrong: the 0.5 floor is natural (we never
believe `ss < 0.5` for a "stranded" prep), but the **upper** bound should
encode the minimum antisense-read rate we are willing to ascribe to RNA.
That bound is biological, not numerical.

---

## 5. Proposed solutions

Ranked by how much I trust the reasoning, not by implementation cost.

### (a) Principled: Beta-Binomial strand likelihood — **recommended**

Replace the binomial with a Beta-Binomial that has a weakly-informative
Beta prior on the antisense rate under RNA:

$$
q_R(k_a \mid n) \;=\; \mathrm{BetaBinom}(k_a ;\; n,\, \alpha_R,\, \beta_R)
$$

with e.g. `α_R = 1 + n_total_train · ss_train` and `β_R = 1 + n_total_train · (1 − ss_train)` derived from the training corpus used to estimate `ss`.  Under gDNA, use a symmetric Beta(1, 1) ⇒ Binomial(n, 0.5).  **No clip needed**: the LLR is automatically bounded by the effective sample size of the prior.

- Pros: uses the actual measurement uncertainty on ss; no magic constants;
  smoothly interpolates between "weak evidence" and "hard evidence" regimes;
  tolerates a few stray antisense reads without flipping.
- Cons: slightly more compute per region (single Γ-function call pair);
  need a sensible rule for `α_R + β_R` (pseudo-count from train).

### (b) Pragmatic: biological lower bound on `1 − ss_c`

```python
STRAND_FLOOR = 0.005          # "never believe antisense rate < 0.5 % under RNA"
ss_c = min(strand_specificity, 1.0 - STRAND_FLOOR)
```

- Pros: one-line change, preserves all other logic, still lets us estimate
  `ss_est=0.99` from data; at ss=1.0 caps `log_ratio_anti` at
  `log(0.5/0.005) ≈ 4.6` nats — informative but non-catastrophic.
- Cons: a magic constant.  Picks one number for all library preps.

### (c) Belt-and-braces: per-fragment Huber-like clip on the LLR

```python
LLR_CAP = 5.0
log_ratio_sense = max(-LLR_CAP, min(LLR_CAP, math.log(0.5/ss_c)))
log_ratio_anti  = max(-LLR_CAP, min(LLR_CAP, math.log(0.5/(1-ss_c))))
```

- Pros: robust to any future clip mis-configuration; decouples the clip
  from the prior.
- Cons: distorts the likelihood for "legitimately informative" antisense
  observations; harder to reason about.

### (d) Meta-fix: always estimate `ss` with a sensible upper bound

The strand-model trainer currently reports `specificity` up to the precision
of the input (`0.9999…` ≈ `1 − 10⁻⁶`).  If the trainer were forced to pass
through a weakly-informative Beta-posterior mean with α,β ≥ 1, the reported
`ss_est` would never cross ~0.9995 for a 200 k-fragment library, which makes
the downstream LLR well-behaved even under the current clip.

- Pros: fixes `ss_est = 1.0` at its source (already fixes vcap oracle).
- Cons: doesn't help degenerate simulations or users who hand-set
  `strand_specificity = 1`.

### What I would do

Combine **(d) + (a)**: train `ss` with a Beta posterior and apply the
Beta-Binomial likelihood in `_strand_llr`.  As a cheap stop-gap until
then, land **(b)** with `STRAND_FLOOR = 0.005` — it is one line, it
immediately removes the failure from vcap, and it costs nothing in
information at real-world strand specificities.

---

## 6. What the fix should *not* break

- ss < 0.5 is still impossible.  Keep that clip on the low side.
- A library with a real 5 % antisense rate should still shift low-coverage
  genic regions toward gDNA when they look 50/50 — just not catastrophically
  so on a single stray read.
- The z-gate on strand evidence (`z ≥ 3`) stays: it's the right
  "activate strand LLR at all?" decision.

---

## 7. Validation plan after a fix

Re-run exactly the same 4×5 mini-stress sweep (`scripts/calibration/
stress_mini/run_sweep.py`) and require:

- `|λ_est/λ_true − 1| ≤ 0.15` across *every* ss row, including ss=1.0.
- `σ_R ≤ 0.8` across all cells (currently 0.62 at ss<1.0, 0.79–0.81 at ss=1.0).
- γ distribution stays unimodal in the "no-gDNA" cells
  (`gdna_fraction = 0`).

Then re-run vcap oracle: expect the −11 % bias to close (oracle σ_R should
drop from 1.85 → ≤ 1.0; λ should track the known truth within a few
percent).

---

## 8. Appendix — numeric scratch

```
ss         log_sense   log_anti    asymmetry
0.5        +0.000      +0.000      inf
0.75       -0.405      +0.693        1.7
0.9        -0.588      +1.609        2.7
0.95       -0.642      +2.303        3.6
0.99       -0.683      +3.912        5.7
0.999      -0.692      +6.215        9.0
1-1e-12    -0.693     +26.938       38.9   ← current behaviour
```

Vcap strand reports (from `summary.json` / `strand_model.strand_specificity`):

```
pristine    oracle_vbem  ss=1.0         π=0.84  π_soft=0.97  μ_R=-4.78  σ_R=1.85
pristine    star_vbem    ss=0.999998    π=0.77  π_soft=0.96  μ_R=-5.50  σ_R=1.82
gdna        oracle_vbem  ss=1.0         π=0.87  π_soft=0.99  μ_R=-4.72  σ_R=1.91
gdna        star_vbem    ss=0.999995    π=0.79  π_soft=0.99  μ_R=-5.54  σ_R=1.87
gdna_nrna   oracle_vbem  ss=1.0         π=0.80  π_soft=0.92  μ_R=-4.30  σ_R=1.79
gdna_nrna   star_vbem    ss=0.999991    π=0.72  π_soft=0.90  μ_R=-4.51  σ_R=1.66
```

All of these sit inside the pathological-ss regime described above.
