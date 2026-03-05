# Prior Architecture and Density-Weighting Analysis

## 1. Current Prior Architecture: Complete Trace

### Step 1: Base prior construction (`locus.py`, `build_locus_em_data`)

```python
prior = np.full(n_components, counter.em_prior, dtype=np.float64)
```

Every component (mRNA, nRNA, gDNA) starts with `prior[k] = em_pseudocount`.
Then certain components get zeroed:
- gDNA if no unspliced fragments at this locus
- nRNA if single-exon transcript
- nRNA if `nrna_init == 0` for that transcript

This is the **flat Dirichlet pseudocount** — call it α_base.

### Step 2: OVR is ADDED on top (`estimator.py`, `run_locus_em`)

```python
prior = prior + coverage_totals / n_ambiguous
```

This takes the coverage-weighted shares (how much "geometric attention"
each component received from ambiguous fragments), normalizes them to
sum to 1.0 across all components, and **adds** them to the existing
prior. The final prior per component is:

    α_k = α_base + coverage_totals_k / N_ambiguous

The OVR part sums to exactly 1.0 across all components. **It was never
a replacement — it was always additive.**

### Step 3: MAP-EM M-step (`estimator.py`, `_em_step`)

```python
theta_new = unambig_totals + em_totals + prior
theta_new /= theta_new.sum()
```

After each iteration, the unnormalized count for component k is:

    θ̃_k = unique_k + em_k + α_base + OVR_k

### Why the old default (α_base = 0.50) caused excessive false positives

Consider FGFR2 with 41 transcripts: that's 2 × 41 + 1 = 83 components.

With α_base = 0.50:
- **Total flat prior mass** = 83 × 0.50 = **41.5** virtual fragments
- **OVR mass** = 1.0 total (negligible by comparison)
- For a zero-truth transcript with no unambig counts, its prior alone is
  0.50 + tiny_OVR, which is enough to attract EM mass in the E-step
  and create a self-sustaining false positive

With α_base = 0.01:
- **Total flat prior mass** = 83 × 0.01 = **0.83** virtual fragments
- **OVR mass** = 1.0 (now *larger* than the flat prior!)
- The flat prior can no longer sustain false positives — only components
  that receive coverage-weighted OVR mass maintain meaningful prior
  support

**The OVR was always doing the right thing — distributing prior mass
according to geometric evidence. But at α = 0.50, the flat prior
completely dominated it (41.5 vs 1.0), making the OVR essentially
irrelevant.** At α = 0.01, the OVR becomes the primary source of prior
information, which is the design intent.

### Why VBEM still beats MAP-EM on FPs

VBEM uses ψ(α_k) − ψ(Σα) instead of log(θ_k) for E-step weights.
For small α_k, ψ(α_k) is extremely negative (e.g. ψ(0.01) ≈ −100.6),
which exponentially suppresses low-evidence components far more
aggressively than log(θ_k). This is an inherent mathematical property
of the digamma function — it acts as a built-in sparsity inducer.

### Phase C validation (10 regions × 4 conditions)

| Config    | MAE   | Spearman | FP  | Dropout |
|-----------|-------|----------|-----|---------|
| MAP 0.01  | 17.72 | 0.8077   | 110 | 178     |
| MAP 0.50  | 19.72 | 0.7760   | 462 | 49      |
| VBEM 0.25 | 18.26 | 0.8111   | 51  | 237     |

MAP 0.01 beats old default in 10/10 regions. VBEM still wins on FPs
(51 vs 110) but has more dropout (237 vs 178).

---

## 2. Density-Weighted EM: Analysis and Critique

### The idea

Currently, there is an asymmetry in how geometric information is used:

| Component       | Coverage-weighted? |
|-----------------|--------------------|
| unambig_totals   | No — raw fragment counts (1.0 per fragment) |
| em_totals       | No — posterior sums (each fragment contributes exactly 1.0 total) |
| OVR prior       | **Yes** — distributed by coverage weights |
| E-step routing  | **Partially** — log-likelihoods include per-fragment effective length |

The proposal: instead of counting fragments as integers, weight each
fragment by a density derived from its coverage probability under the
geometric model (e.g. the trapezoid). Convert raw counts N into
weighted density D, then solve EM on D.

### What this means concretely

Under a uniform-coverage model, the probability of observing a
fragment starting at position x on transcript t is:

    p(x | t) = trapezoid(x, L_t, frag_len) / eff_len_t

where the trapezoid accounts for the triangular ramp-up at the 5' end
and ramp-down at the 3' end. In the middle of the transcript,
p(x | t) = 1 / eff_len_t. At the edges, p(x | t) < 1 / eff_len_t.

The density weight for fragment f on transcript t would be:

    w(f, t) = 1 / [p(x_f | t) × eff_len_t]

For a fragment in the interior: w = 1.0 (rectangle region).
For a fragment near the edge: w > 1.0 (triangular ramp region).

A fragment observed at the tip of the 3' end is "more surprising" — it
landed in a region with lower coverage probability, so it provides
more information per observation. Weighting it up compensates for the
positional sampling bias.

### What changes in the EM

**Current M-step** (each fragment = 1.0):

    θ̃_k = Σ_i [1.0 × p(z_i = k | x_i, θ)] + unique_k + α_k

**Density-weighted M-step** (each fragment = w_ik):

    θ̃_k = Σ_i [w_ik × p(z_i = k | x_i, θ)] + Σ_{j ∈ unique_k} w_jk + α_k

The E-step posterior p(z_i = k) doesn't change — it still uses the
same log-likelihoods for routing. What changes is the *magnitude* of
each fragment's contribution after routing.

### Theoretical analysis

#### Relationship to effective length correction

The standard RNA-seq quantification model (salmon, kallisto) treats
effective length as a *summary statistic*:

    abundance_k ∝ count_k / eff_len_k

This is a single divisor applied after counting. It is the
maximum-likelihood estimator under a uniform-coverage Poisson model.

Per-fragment density weighting is a more granular version: instead of
one correction per transcript, it applies a correction per
fragment-transcript pair. Under the uniform model with infinite data,
these are **asymptotically equivalent** — the average of per-fragment
weights converges to 1/eff_len × eff_len = 1.0 for the correct
transcript.

The difference is in **finite samples**: per-fragment weighting has
higher variance but can reduce bias when the observed fragment
positions are not representative of the full coverage distribution.

#### When density weighting helps

1. **Short transcripts** — Most of the transcript length is in the
   ramp regions, so many fragments land in low-probability zones.
   Density weighting inflates these fragments, effectively saying
   "this short transcript is punching above its weight."

2. **Competing long vs short isoforms** — A short isoform's fragments
   are disproportionately in ramp zones. Without density weighting,
   they contribute 1.0 each, same as a long isoform's interior
   fragments. With weighting, the short isoform's ramp fragments
   contribute > 1.0, partially counteracting the length bias.

3. **Sparse loci** — With few fragments, the positional distribution
   may be skewed. Density weighting corrects for this.

#### When density weighting hurts

1. **Edge effects / high variance** — Fragments at the very tip of a
   transcript (within ~10 bp of the end) have very low coverage
   probability, so their density weight approaches infinity. This
   requires truncation/capping, introducing an arbitrary parameter.

2. **Model misspecification** — Real RNA-seq has 5'/3' bias, GC bias,
   and other non-uniformities. The trapezoid model assumes uniform
   coverage. If real coverage is non-uniform, density weighting based
   on the trapezoid may *increase* rather than decrease bias.

3. **Already handled by per-fragment effective length** — The current
   code already bakes per-fragment effective length corrections into
   each candidate's log-likelihood during the scan phase:

       log_lik[f, k] -= log(max(transcript_len - footprint + 1, 1))

   This means the E-step routing already "knows" about positional
   effects. The remaining gap is only in the M-step accumulation
   (fragment weight = 1.0 vs density-weighted).

4. **Interaction with OVR** — The OVR prior already uses coverage
   weights. If we also density-weight the data, we need to ensure
   consistency. The OVR's 1.0 total virtual read would need to be
   recalibrated against the total weighted density D rather than
   count N.

### Critical evaluation

The key question: **Does the current per-fragment effective length
correction in the log-likelihood achieve the same effect as density
weighting?**

**Not quite, but close.** The per-fragment correction adjusts the
*E-step routing* (which candidate gets the fragment), but not the
*M-step magnitude* (how much the fragment contributes). A fragment at
a transcript edge correctly gets routed away from long transcripts
(where its position is improbable) and toward short ones (where it's
more probable). But once routed, it still contributes 1.0.

Density weighting would additionally say: "this edge fragment, having
been correctly routed to the short transcript, should count as >1.0
because it was less likely to be observed there." This is an
inverse-probability-of-sampling correction.

**However**, in practice:
- The per-fragment likelihood correction already captures most of the
  positional signal (in the routing step)
- The remaining magnitude correction is a second-order effect
- The EM is already converging to good solutions (MAE ~17.7 with
  MAP 0.01)
- The dominant error source (110 FPs) comes from the prior structure,
  not from positional weighting of counts

### Recommendation: simpler path first

Before exploring density weighting (which adds complexity and new
parameters), there's a simpler intervention that addresses the same
root cause:

**Fold the flat prior into the OVR.**

Currently: α_k = α_base + OVR_k (two additive terms)

Proposed: α_k = ε + γ × OVR_k

where:
- ε ≈ 1e-10 is a numerical-stability floor (not a real prior)
- γ is a scale factor (default 1.0 = one virtual read total)

This eliminates the flat prior entirely. A component with zero
coverage weight gets α_k ≈ ε — effectively zero prior — and the EM
will drive it to zero. Only components with geometric evidence receive
meaningful prior support.

**Advantages over density weighting:**
1. Zero new parameters (γ = 1.0 is the natural default)
2. No truncation/capping issues
3. Directly addresses the FP mechanism (prior inflating dead
   components)
4. Preserves existing per-fragment effective length correction
5. Quick to implement and test

**Advantages of density weighting (for future exploration):**
1. More theoretically principled — consistent treatment of all
   fragments
2. Could reduce short-transcript bias
3. May improve accuracy in sparse loci
4. Not mutually exclusive with the OVR-only prior

### Suggested path forward

1. **Immediate (Phase D)**: Replace α_base with pure OVR + ε floor.
   Test on 10 regions. This directly tests whether eliminating the
   flat prior closes the FP gap with VBEM while avoiding VBEM's
   dropout penalty.

2. **If Phase D succeeds**: Ship MAP-EM with OVR-only prior as the
   default. VBEM remains available as an option.

3. **Future exploration**: Density weighting as a follow-on experiment
   for accuracy improvement on short transcripts and sparse loci.
   Requires careful design of the weight function (truncation,
   interaction with OVR, effect on convergence).
