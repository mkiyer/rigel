# Pool-Separated Likelihood-Based Pruning

**Date**: 2026-03-22
**Status**: Design (not yet implemented)
**Parent**: [implementation_plan.md](implementation_plan.md) — Problem 2
**Addresses**: Problem 2 (nRNA hard overhang gate bug) and partially Problem 1 (alignment-error rescue)

---

## 1. Motivation

### Why the current hard overhang gate is broken

The current hard overhang gate in `scoring.cpp` operates on a
single signal — **intronic overhang** — and applies a hard binary filter:

```
m_min_oh = min(oh for all candidates)
keep only candidates where oh == m_min_oh
```

This was adequate when all candidates were annotated mRNA transcripts
with similar exon structures. The introduction of synthetic nRNA
transcripts broke the gate because nRNA transcripts are single-exon
spans covering entire genomic loci (introns included), so they always
have `oh = 0`. The hard overhang gate therefore **always keeps nRNA and discards
any mRNA candidate with oh > 0** — even when the overhang is just 1 bp
from an alignment boundary error.

### Why the hard overhang gate was always suboptimal

Even without the nRNA bug, the overhang-only hard overhang gate is problematic.
The scoring engine computes a full log-likelihood from **four independent
signals**:

$$\ell = \ell_{\text{strand}} + \ell_{\text{FL}} + \text{oh} \cdot \lambda_{\text{oh}} + \text{nm} \cdot \lambda_{\text{nm}}$$

Gating on overhang alone discards the information from strand, fragment
length, and mismatch. A candidate with `oh = 1` but perfect strand
match and ideal fragment length is more likely correct than a candidate
with `oh = 0` but antisense strand and implausible fragment length. The
hard hard overhang gate cannot distinguish these cases.

### Design goals

1. **Pool separation**: nRNA and mRNA candidates must never compete in
   the same gate — they have fundamentally different expected scoring
   profiles (mRNA has exonic structure with potential overhang; nRNA is
   a single span with oh ≡ 0).

2. **Likelihood-based pruning**: Replace overhang-only gating with
   pruning based on the full log-likelihood $\ell$, which integrates
   all four scoring signals.

3. **No magic numbers**: Derive the pruning threshold from a principled
   probabilistic argument tied to the EM's posterior computation.

4. **Conservative default**: Err on the side of keeping more candidates.
   Excess candidates increase equivalence class sizes (computational
   cost) but do not degrade accuracy — the EM handles them correctly.

5. **Subsume Problem 1**: With pool-separated likelihood pruning, mRNA
   candidates with small overhang (1–5 bp alignment errors) are
   naturally retained alongside nRNA candidates, without any special
   rescue mechanism.

---

## 2. Design: Pool Separation

### Two pools

Every transcript in the index is classified into one of two pools:

| Pool | Composition | Characteristics |
|------|-------------|-----------------|
| **mRNA** | All annotated transcripts + annotated nascent-equivalents | Multi-exon or single-exon; has exonic structure; overhang is meaningful |
| **nRNA** | All synthetic nRNA transcripts (`is_synthetic_nrna = True`) | Single-exon spanning full locus; oh ≡ 0 always; fragment length = genomic span |

The pool assignment is determined at index build time and is immutable.

### Pool identity in the scorer

A new boolean array `t_is_nrna[n_transcripts]` is passed from Python to
the C++ `NativeFragmentScorer` constructor. This array is `true` for
synthetic nRNA transcripts and `false` for everything else.

The scorer uses this array to partition candidates during the pruning
step. Each pool is pruned independently using its own `pool_best_ll`.

### Why two pools, not three or more?

gDNA is already handled as a separate single candidate per unit
(computed at unit finalization, not part of the candidate CSR). There is
no need for a gDNA pool in the pruning gate.

Future extensions (e.g., separate pools for single-exon mRNA vs
multi-exon mRNA) are possible but not needed now. The mRNA vs nRNA
split addresses the critical structural asymmetry.

---

## 3. Design: Likelihood-Based Pruning

### Pruning criterion

Within each pool, after scoring all candidates:

1. Find `pool_best_ll` = max log-likelihood among candidates in this pool
2. **Keep** candidate $j$ if:
$$\text{pool\_best\_ll} - \ell_j \leq \Delta$$
3. **Prune** candidate $j$ otherwise

The union of survivors from both pools forms the final candidate set for
the EM unit.

### Why per-pool best (not global best)?

If we used a global `best_ll` across both pools, nRNA candidates (with
`oh = 0` → higher raw $\ell$) would set the bar, and mRNA candidates
with overhang > 0 would be pruned relative to nRNA's inflated
likelihood. This is exactly the current bug in a softer form.

Per-pool best ensures that mRNA candidates compete only with other mRNA
candidates, and nRNA candidates compete only with other nRNA candidates.
The EM then arbitrates between pools using its global mixing proportions.

---

## 4. Design: The Δ Threshold

### The principled derivation

The EM computes posterior probabilities for each candidate:

$$\gamma_j = \frac{\theta_j \cdot L_j}{\sum_k \theta_k \cdot L_k}$$

where $\theta_j$ is the mixing proportion (expression level) and
$L_j = \exp(\ell_j)$ is the likelihood.

**Question**: Can we safely prune candidate $j$ without affecting the
EM's convergence?

**Answer**: If $\gamma_j$ would be negligible regardless of the mixing
proportions $\theta$, then yes. The maximum possible $\gamma_j$ — under
the most favorable $\theta$ — is bounded by the likelihood ratio:

$$\gamma_j \leq \frac{L_j}{L_{\text{best}}} = \exp(\ell_j - \ell_{\text{best}}) = \exp(-\Delta_j)$$

where $\Delta_j = \ell_{\text{best}} - \ell_j \geq 0$.

This bound is achieved when all mixing weight is concentrated on
candidate $j$ (i.e., $\theta_j \to \infty$, all others zero). It is an
**upper bound** on the posterior, valid for any $\theta$ vector.

> **Interpretation**: If $\Delta_j > -\ln(\varepsilon)$, then the
> maximum possible posterior of candidate $j$ is less than
> $\varepsilon$, regardless of mixing proportions.

### Important caveat

The bound above assumes we are comparing one candidate to only the best
candidate. In reality, $\gamma_j$ depends on all candidates' mixing
proportions. The bound is tightest when we consider:

$$\gamma_j = \frac{\theta_j \cdot L_j}{\theta_j \cdot L_j + \theta_{\text{best}} \cdot L_{\text{best}} + \sum_{k \neq j, \text{best}} \theta_k \cdot L_k}$$

Ignoring the sum over other candidates (which only makes the bound
tighter):

$$\gamma_j \leq \frac{\theta_j \cdot L_j}{\theta_j \cdot L_j + \theta_{\text{best}} \cdot L_{\text{best}}}$$

This is maximized when $\theta_j / \theta_{\text{best}} \to \infty$,
giving $\gamma_j \to 1$. So in the absolute worst case, **no Δ
threshold can guarantee a negligible posterior**.

However, the EM constrains $\theta$ globally across all fragments.
If candidate $j$ has a massive likelihood deficit on most fragments, its
global $\theta_j$ will be tiny, making its posterior negligible
everywhere. The only scenario where $j$ could have large posterior
despite a large $\Delta_j$ is if $j$ is vastly more highly expressed
than all other candidates in the locus — and even then, $j$ would need
to be a strong competitor on *other* fragments to accumulate high
$\theta_j$.

### Practical implication

For the pruning threshold, we seek a $\Delta$ that is large enough that
**in practice**, pruned candidates never contribute meaningfully to the
EM. The posterior bound gives us a scale:

| $\varepsilon$ (max posterior) | $\Delta = -\ln(\varepsilon)$ | Interpretation |
|------|------|----------------|
| $10^{-2}$ | 4.6 | Aggressive — prunes 1-bp-overhang equivalents |
| $10^{-3}$ | 6.9 | Moderate |
| $10^{-4}$ | 9.2 | Conservative |
| $10^{-6}$ | 13.8 | Very conservative — retains candidates 6 bp overhang worse |
| $10^{-8}$ | 18.4 | Ultra-conservative — effectively no pruning |

### Recommended Δ

**Default**: $\varepsilon = 10^{-4}$, giving $\Delta \approx 9.21$.

Rationale:
- At the default overhang penalty (`log(0.01) ≈ -4.605` per bp), Δ = 9.21
  retains candidates with up to **2 bp overhang** beyond the best
  candidate's overhang, *even if all other signals are identical*. Since
  most alignment errors are 1–2 bp, this is sufficient.
- At $\varepsilon = 10^{-4}$, a pruned candidate's maximum posterior is
  0.01%. Even if its mixing proportion is 100× that of the best
  candidate, its posterior remains < 1%.
- The equivalence class expansion is bounded: candidates with
  $\Delta > 9.2$ are typically so far from the best that they represent
  entirely different transcript families (e.g., wrong strand, wrong
  locus). Pruning these has negligible computational cost.

### Parameter name and configuration

```python
@dataclass(frozen=True)
class FragmentScoringConfig:
    overhang_log_penalty: float = math.log(0.01)
    mismatch_log_penalty: float = math.log(0.1)
    gdna_splice_penalties: dict[int, float] | None = None
    # NEW:
    pruning_min_posterior: float = 1e-4  # ε
```

The user sets `pruning_min_posterior` (ε), and the code computes
$\Delta = -\ln(\varepsilon)$ internally. This is more interpretable
than exposing Δ directly — users think in terms of "what's the minimum
posterior I want to preserve?" rather than "what's the log-likelihood
ratio cutoff?"

**Invariants**:
- `pruning_min_posterior ∈ (0, 1)` — enforced at construction
- Setting `pruning_min_posterior = 0` → $\Delta = \infty$ → no pruning
  (keep everything)
- Setting `pruning_min_posterior = 1` → $\Delta = 0$ → keep only the
  best candidate per pool (maximally aggressive)

---

## 5. Design: Scoring Flow After Modification

### Current flow (unique mapper)

```
for each candidate k in [start, end):
    score → MrnaScored{t_idx, oh, log_lik, ...}
    track m_min_oh

for each scored:
    if oh == m_min_oh:
        emit candidate
```

### New flow (unique mapper)

```
for each candidate k in [start, end):
    score → MrnaScored{t_idx, oh, log_lik, ...}
    if t_is_nrna[t_idx]:
        add to nrna_pool, track nrna_best_ll
    else:
        add to mrna_pool, track mrna_best_ll

for each scored in mrna_pool:
    if mrna_best_ll - log_lik ≤ Δ:
        emit candidate

for each scored in nrna_pool:
    if nrna_best_ll - log_lik ≤ Δ:
        emit candidate
```

### Current flow (multimapper)

```
for each hit across all alignments:
    for each candidate:
        score → merge into merged_mrna[t_idx] (best per transcript)

mrna_min_oh = min(mc.overhang for mc in merged_mrna)
for each (t_idx, mc) in merged_mrna:
    if mc.overhang == mrna_min_oh:
        emit candidate
```

### New flow (multimapper)

```
for each hit across all alignments:
    for each candidate:
        score → merge into merged_mrna[t_idx] (best per transcript)

# Pool separation
mrna_best_ll = max(mc.log_lik for mc in merged_mrna where !t_is_nrna[t_idx])
nrna_best_ll = max(mc.log_lik for mc in merged_mrna where  t_is_nrna[t_idx])

for each (t_idx, mc) in merged_mrna:
    if t_is_nrna[t_idx]:
        if nrna_best_ll - mc.log_lik ≤ Δ:
            emit candidate
    else:
        if mrna_best_ll - mc.log_lik ≤ Δ:
            emit candidate
```

### Edge cases

| Condition | Behavior |
|-----------|----------|
| mRNA pool empty, nRNA pool non-empty | Emit nRNA survivors only |
| nRNA pool empty, mRNA pool non-empty | Emit mRNA survivors only |
| Both pools empty | Fragment gated out (stat_gated++) |
| All candidates in one pool pruned | Pool contributes zero candidates; other pool unaffected |
| Single candidate in pool | Always retained (Δ ≥ 0 trivially) |
| $\Delta = +\infty$ (ε → 0) | No pruning at all, every candidate retained |

---

## 6. Implementation: C++ Changes (scoring.cpp)

### 6a. New member: `t_is_nrna_`

```cpp
class NativeFragmentScorer {
    // ... existing members ...
    std::vector<uint8_t> t_is_nrna_;  // uint8[n_transcripts] — 1=nRNA, 0=mRNA
    double max_ll_delta_;             // Δ = -log(ε), pre-computed
```

Using `uint8_t` instead of `bool` because `std::vector<bool>` is a
space-optimized bitset with proxy references that complicates element
access in tight loops.

### 6b. Constructor changes

Add two new parameters to the `NativeFragmentScorer` constructor:

```cpp
NativeFragmentScorer(
    // ... existing parameters ...
    u8_1d  t_is_nrna_arr,           // NEW: uint8[n_transcripts]
    double pruning_max_ll_delta)     // NEW: Δ = -log(ε)
```

Copy `t_is_nrna_arr` into `t_is_nrna_` vector. Store
`pruning_max_ll_delta` as `max_ll_delta_`.

### 6c. Unique mapper scoring (score_chunk_impl)

Replace the `m_min_oh` hard overhang gate with pool-separated likelihood pruning.

**Before** (lines ~810-870):
```cpp
int32_t m_min_oh = 0x7FFFFFFF;
// ... score all candidates, track m_min_oh ...
for (int j = 0; j < m_n; ++j) {
    if (m_scored[j].oh == m_min_oh) {
        // emit
    }
}
```

**After**:
```cpp
// Score all candidates (unchanged scoring logic)
// ...

// Pool-separated likelihood pruning
double mrna_best = NEG_INF;
double nrna_best = NEG_INF;
for (int j = 0; j < m_n; ++j) {
    if (t_is_nrna_[m_scored[j].t_idx])
        nrna_best = std::max(nrna_best, m_scored[j].log_lik);
    else
        mrna_best = std::max(mrna_best, m_scored[j].log_lik);
}

for (int j = 0; j < m_n; ++j) {
    auto& s = m_scored[j];
    double pool_best = t_is_nrna_[s.t_idx] ? nrna_best : mrna_best;
    if (pool_best - s.log_lik <= max_ll_delta_) {
        // emit candidate
    }
}
```

The scoring logic for each candidate (strand, FL, overhang, NM) is
completely unchanged. Only the gate decision changes.

### 6d. Multimapper scoring (flush_mm_group)

Same transformation: replace `mrna_min_oh` gate with pool-separated
likelihood pruning.

**Before** (lines ~560-575):
```cpp
int32_t mrna_min_oh = 0x7FFFFFFF;
for (auto& [t_idx, mc] : merged_mrna)
    if (mc.overhang < mrna_min_oh)
        mrna_min_oh = mc.overhang;

for (auto& [t_idx, mc] : merged_mrna)
    if (mc.overhang == mrna_min_oh) { /* emit */ }
```

**After**:
```cpp
double mrna_best = NEG_INF;
double nrna_best = NEG_INF;
for (auto& [t_idx, mc] : merged_mrna) {
    if (t_is_nrna_[t_idx])
        nrna_best = std::max(nrna_best, mc.log_lik);
    else
        mrna_best = std::max(mrna_best, mc.log_lik);
}

for (auto& [t_idx, mc] : merged_mrna) {
    double pool_best = t_is_nrna_[t_idx] ? nrna_best : mrna_best;
    if (pool_best - mc.log_lik <= max_ll_delta_) {
        // emit candidate
    }
}
```

### 6e. Nanobind binding

Update the `NativeFragmentScorer` nanobind class definition to include
the two new constructor parameters. The Python side computes
`max_ll_delta = -math.log(pruning_min_posterior)` and passes it as a
`double`.

---

## 7. Implementation: Python Changes

### 7a. FragmentScoringConfig (config.py)

Add `pruning_min_posterior` parameter:

```python
@dataclass(frozen=True)
class FragmentScoringConfig:
    overhang_log_penalty: float = math.log(0.01)
    mismatch_log_penalty: float = math.log(0.1)
    gdna_splice_penalties: dict[int, float] | None = None
    pruning_min_posterior: float = 1e-4  # NEW
```

### 7b. FragmentScorer.from_models (scoring.py)

Pass `is_synthetic_nrna` array and computed Δ to the native scorer:

```python
t_is_nrna_arr = index.t_df["is_synthetic_nrna"].values.astype(np.uint8)
max_ll_delta = -math.log(max(scoring_config.pruning_min_posterior, 1e-300))

native_ctx = NativeFragmentScorer(
    # ... existing parameters ...
    t_is_nrna_arr=np.ascontiguousarray(t_is_nrna_arr),
    pruning_max_ll_delta=max_ll_delta,
)
```

### 7c. CLI (cli.py)

Expose `--pruning-min-posterior` as a command-line option for
`rigel quant`:

```
--pruning-min-posterior FLOAT   Minimum posterior threshold for candidate
                                pruning (default: 1e-4). Lower values are
                                more conservative (keep more candidates).
                                Set to 0 to disable pruning entirely.
```

---

## 8. Design Analysis: Why This Works

### The nRNA hard overhang gate bug is fixed

With pool separation, nRNA's `oh = 0` no longer establishes a baseline
that prunes mRNA candidates. Each pool competes internally:

- mRNA pool: candidates with `oh = 0, 1, 2, 3` compete among themselves.
  A candidate with `oh = 2` (penalty ≈ -9.2) is within
  $\Delta = 9.2$ of a candidate with `oh = 0`,
  assuming other signals are similar. **It survives**.
- nRNA pool: typically has one or a few candidates (per locus span).
  Almost always retained.

### Alignment-error rescue is built-in

Fragments with 1–3 bp overhang from minimap2 alignment errors:
- Under the old hard overhang gate: mRNA pruned (oh=1-3 vs nRNA oh=0). **BROKEN**.
- Under pool-separated pruning: mRNA with oh=1 has penalty ≈ -4.6.
  mRNA best candidate might have oh=0 with similar other signals.
  Δ = 4.6 < 9.2 → **mRNA candidate retained**. Problem 1 is solved
  without any special-case logic.

### Fragment length training is fixed

The FL model corruption cascade from the implementation plan:

```
Hard overhang gate prunes mRNA → nRNA wins → FL = genomic span → FL model corrupted
```

With pool separation, mRNA candidates survive hard overhang gate → fragment is
assigned to mRNA in the EM → FL computed against mRNA exonic structure
→ FL model trains on correct values.

### Equivalence class impact

The current hard overhang gate is maximally aggressive: it keeps only `oh == min_oh`
candidates. The new likelihood pruning is strictly more permissive.
Expected equivalence class growth:

- **Typical fragment** (oh=0 for all mRNA candidates, no nRNA in locus):
  Same as before — all candidates have similar $\ell$, all retained.
- **Fragment with overhang** (oh=1–3 on some mRNA, oh=0 on others):
  Now retains the oh>0 candidates that were previously pruned.
  Equivalence class grows by the number of overhang candidates.
  Bounded by the total number of transcripts overlapping the read.
- **Fragment in a locus with nRNA**: Now has nRNA candidates in addition
  to mRNA candidates. Equivalence class grows by 1 (or a few nRNA
  transcripts per locus).

For a typical locus with 5–20 transcripts, the equivalence class might
grow from ~5 to ~8 candidates. For mega-loci with thousands of
transcripts, the Δ threshold becomes important — candidates with very
poor FL or wrong strand (Δ > 9.2) are still pruned.

### Computational cost

The EM solver's time complexity per iteration scales as
$O(\sum_e n_e \cdot k_e)$ where $n_e$ = units in equivalence class $e$
and $k_e$ = components per class. Increasing $k_e$ by ~50% for
overhang-affected fragments adds ~50% to per-iteration cost for those
classes. Since overhang-affected fragments are ~7% of total, the
overall EM slowdown is approximately **3–5%**. This is negligible
compared to the accuracy improvement.

---

## 9. Invariants and Correctness Guarantees

### Monotonicity

As $\varepsilon \to 0$ (Δ → ∞), the pruning gate admits all candidates,
converging to the "no gating" baseline. Accuracy should be monotonically
non-decreasing as ε decreases (more candidates → EM has more information).

### Safety

The only way pruning can degrade accuracy is by removing the true source
transcript before the EM. With per-pool pruning at Δ = 9.2:
- The true mRNA source has $\ell$ at most ~4.6 below the best mRNA
  candidate (from 1 bp alignment error). Well within Δ.
- The true nRNA source (if the fragment is genuinely nascent RNA) has
  $\ell$ ≈ best nRNA candidate (since nRNA transcripts span the whole
  locus). Always retained.

### Determinism

Pool assignment is deterministic (from index). Likelihood computation is
deterministic. The pruning gate is a deterministic function of $\ell$
and $\Delta$. No new sources of non-determinism are introduced.

---

## 10. Testing Strategy

### Unit tests (regression gate)

All 998 existing tests must pass. Since the tests use the default
config, they will use the new pruning with $\varepsilon = 10^{-4}$.
Tests that depend on specific hard overhang gate behavior may need adjustment — but
only if they were testing the hard overhang gate directly, not just its downstream
effects.

### New unit tests

1. **Pool separation**: Create a fragment with mRNA candidates (oh=0,1)
   and nRNA candidate (oh=0). Verify that all three survive pruning.
   Under the old hard overhang gate, only oh=0 candidates would survive (nRNA + oh=0
   mRNA), pruning the oh=1 mRNA.

2. **Likelihood pruning boundary**: Create candidates with known
   log-likelihoods. Set ε to place the boundary exactly at one
   candidate. Verify it is pruned/retained correctly.

3. **Empty pool**: Fragment with no nRNA candidates → only mRNA pool
   active. Fragment with only nRNA candidates → only nRNA pool active.

4. **All pruned**: Both pools empty after pruning → fragment gated out.

5. **Conservative extreme**: Set ε = 0 (Δ = ∞) → all candidates
   retained. Verify equivalence class sizes match "no gating."

### Benchmark validation

Run the pristine benchmark (10M mRNA fragments, zero gDNA, zero nRNA):

| Metric | Target | Rationale |
|--------|--------|-----------|
| MAE | < 3.0 | Approaching oracle (1.85) and salmon (2.91) |
| nRNA siphon | < 10K | Oracle: 3,625 |
| gDNA false pos | < 5K | Oracle: 1,840 |
| FL mode (SPLICED_ANNOT) | 290–310 | Oracle: 299 |
| Spearman | > 0.85 | Oracle: 0.880, salmon: 0.894 |

### Sensitivity analysis (Δ sweep)

Run the pristine benchmark at multiple ε values to characterize the
accuracy–performance tradeoff:

| ε | Δ | Expected behavior |
|---|---|-------------------|
| 1.0 | 0.0 | Only best candidate per pool (most aggressive) |
| 0.01 | 4.6 | Aggressive — ~1 bp overhang tolerance |
| 1e-4 | 9.2 | Default — ~2 bp overhang tolerance |
| 1e-6 | 13.8 | Conservative — ~3 bp overhang tolerance |
| 1e-8 | 18.4 | Ultra-conservative — ~4 bp tolerance |
| 0 | ∞ | No pruning (all candidates retained) |

Report for each: MAE, equiv class size distribution, EM iterations,
wall-clock time.

---

## 11. Relationship to Other Problems

### Problem 1 (Alignment-Error Rescue)

**Subsumed**. With pool-separated likelihood pruning, mRNA candidates
with 1–3 bp overhang are retained alongside nRNA candidates. The EM
assigns fragments to the correct transcript based on global signal. No
additional "rescue" mechanism is needed.

If the benchmark shows residual MAE gap after this fix, the gap is
attributable to other sources (e.g., FL model quality, gap SJ correction)
rather than candidate pruning.

### Problem 3 (Robust FL Estimation)

**Partially addressed**. The FL model corruption was caused by nRNA
winning the hard overhang gate → FL computed against nRNA's genomic span → inflated FL
observations. With pool separation, mRNA candidates survive → FL
computed against mRNA exonic structure → correct observations.

However, the FL robustness enhancements (outlier rejection, cross-
validation) from the implementation plan should still be implemented as
defense-in-depth. They protect against other corruption vectors
(e.g., genuinely misaligned reads, structural variant artifacts).

### Gap SJ correction (overlap fix)

The overlap-based gap SJ correction currently in `resolve_context.h`
should be re-evaluated after this fix. Its benefit may be
complementary (fixing inter-mate gap FL computation) or redundant
(the FL benefit came from partially compensating for the hard overhang gate bug).
Re-run benchmark with and without the overlap fix after implementing
pool-separated pruning.

---

## 12. Implementation Checklist

### Phase 1: Core C++ changes

- [ ] Add `t_is_nrna_` (uint8 vector) and `max_ll_delta_` to
  `NativeFragmentScorer` class
- [ ] Extend constructor to accept new parameters, copy into member
  storage
- [ ] Modify unique mapper scoring path (`score_chunk_impl`): replace
  `m_min_oh` hard overhang gate with pool-separated likelihood pruning
- [ ] Modify multimapper scoring path (`flush_mm_group`): replace
  `mrna_min_oh` hard overhang gate with pool-separated likelihood pruning
- [ ] Update nanobind bindings to expose new constructor parameters

### Phase 2: Python interface changes

- [ ] Add `pruning_min_posterior` to `FragmentScoringConfig`
- [ ] Modify `FragmentScorer.from_models()` to extract `is_synthetic_nrna`
  from index and compute `max_ll_delta = -log(ε)`
- [ ] Pass new arrays/params to `NativeFragmentScorer` constructor
- [ ] Update CLI (`rigel quant`) to expose `--pruning-min-posterior`

### Phase 3: Testing

- [ ] Run existing 998 unit tests (must pass)
- [ ] Write new unit tests for pool separation and likelihood pruning
- [ ] Run pristine benchmark, compare against baseline and oracle
- [ ] Run Δ sensitivity sweep

### Phase 4: Evaluation

- [ ] Analyze benchmark results
- [ ] Decide on default ε value (may adjust from 1e-4)
- [ ] Re-evaluate gap SJ overlap fix
- [ ] Document results and update implementation plan status

---

## 13. Open Questions

1. **Det-unambig fast path**: Currently, FRAG_UNAMBIG + SPLICE_ANNOT
   fragments skip the scoring loop entirely and go straight to
   det-unambig tracking (no EM involvement). These fragments have a
   single candidate transcript and are deterministically assigned.
   **Question**: Do any of these fragments have nRNA as their single
   candidate? If so, they bypass the pruning gate entirely.
   **Resolution**: No — det-unambig requires SPLICE_ANNOT, meaning the
   fragment is spliced. Spliced fragments are incompatible with nRNA
   (single-exon synthetic transcripts cannot produce splice junctions).
   The det-unambig fast path is safe.

2. **Multimapper merge resolution**: In the multimapper path, candidates
   are merged across alignments with "lower overhang wins, ties broken
   by log_lik." After the pool-separated change, this merge still makes
   sense (we want the best representation per transcript regardless of
   pool). But should we verify that the merge doesn't create pathological
   interactions with pool pruning?
   **Status**: Needs further exploration after implementation. The merge
   logic is per-transcript (each transcript keeps its best alignment),
   which is independent of pool membership. Pool pruning operates on the
   merged results. No obvious pathology, but monitor empirically.

3. **Equivalence class explosion for mega-loci**: Mega-loci with
   thousands of transcripts might see significant equivalence class
   growth with permissive pruning. Should we monitor equiv class sizes
   in the benchmark and consider a hard cap (e.g., max 100 candidates
   per unit) as a safety valve?
   **Status**: Monitor equiv class size distribution in benchmarks. No
   hard cap for now — we want to observe the natural behavior first.

4. **nRNA pool pruning**: In many loci, there is only one nRNA
   transcript. If there are multiple (from overlapping loci with
   different strands), they might have very different strand likelihoods.
   **Resolution**: The per-pool Δ threshold handles this correctly.
   Opposite-strand nRNAs will have different strand likelihoods; the
   best nRNA sets the pool baseline and the other survives if within Δ.

5. **Interaction with FL training**: FL training occurs during the BAM
   scan (before scoring). The BAM scanner uses `get_unique_frag_length()`
   which requires all candidates to agree on FL.
   After this fix, the hard overhang gate in the *scoring* path changes, but FL
   training in the BAM scanner does NOT use the scorer — it uses the
   raw resolution from `resolve_context.h`. So this fix does NOT
   directly change FL training. However, indirectly, the hard overhang gate fix
   affects which fragments are classified as FRAG_UNAMBIG vs
   FRAG_AMBIG, which affects which fragments provide FL observations.
   **Status**: FL training needs separate focused work to ensure
   correctness in light of these changes. Planned as a follow-up task.
