# Positional Bias Model for Coverage-Weighted Quantification

## Motivation

The current effective-length correction assumes a **uniform** positional distribution
of fragment starts within each transcript:

```
P(f | t) = 1 / max(t_len − l_f + 1, 1)
```

Real RNA-seq exhibits systematic positional biases (5′/3′ enrichment, GC effects,
fragmentation non-uniformity). A positional bias model replaces the uniform
assumption with a per-base bias profile, improving quantification accuracy for
transcripts whose coverage deviates from uniform.

## Mathematical Framework

### Bias Profile

Each transcript `t` of length `L` has a base-level bias array:

```
b[0], b[1], ..., b[L−1]    (all > 0)
```

For the uniform (null) model, `b[i] = 1.0` for all `i`.

### Prefix Sum

Pre-compute the prefix sum for O(1) range queries:

```
B[0] = 0
B[i] = b[0] + b[1] + ... + b[i−1]
```

So `sum(b[s..e−1]) = B[e] − B[s]` in O(1).

### Fragment Weight

For a fragment mapped to transcript positions `[s, e)` (0-based, half-open):

```
W(f | t) = (B[e] − B[s]) / (e − s)
```

This is the **mean bias** over the fragment footprint. Under uniform bias
(`b[i] = 1.0`), `W = 1.0`.

### Bias-Aware Effective Length

For a given fragment length `l_f` and transcript length `L`:

```
                  L − l_f
L̃_t(l_f) = Σ     W(v, v + l_f | t)
               v=0
```

where `W(v, v + l_f | t) = (B[v + l_f] − B[v]) / l_f`.

This simplifies to:

```
L̃_t(l_f) = (1 / l_f) Σ_{v=0}^{L − l_f} (B[v + l_f] − B[v])
```

Under uniform bias: `L̃ = L − l_f + 1` (exactly the current behaviour).

### Fragment Likelihood

```
P(f | t) = W(f | t) / L̃_t(l_f)
```

In log space (what gets stored in `log_liks`):

```
log P(f | t) = log W(f | t) − log L̃_t(l_f)
```

Under uniform bias this collapses to:

```
log P(f | t) = log(1.0) − log(L − l_f + 1) = −log(max(t_len − l_f + 1, 1))
```

which is exactly the current effective-length correction baked into `log_liks`.

## Implementation Phases

### Phase 0 — CSR Extension + Uniform Bias Scaffold (this PR)

**Goal:** Extend the data pipeline to carry per-candidate transcript coordinates
and create the `bias.py` module with uniform profiles. Verify that the new code
path produces identical results to the current hard-coded effective-length
correction.

#### 0a. ScanData CSR Extension

- Add `tx_starts: list[int]` and `tx_ends: list[int]` to `ScanData` (estimator.py)
- These store the 0-based transcript-space `[start, end)` for each candidate
- Populated during scan in `_score_wta_mrna` and `_score_wta_nrna` (scan.py)
  - `tx_start` and `tx_end` are already computed for coverage weight — just store them
  - For nRNA: store `tx_start=0, tx_end=t_len` (full transcript)
- Propagate to `LocusEMInput` as `tx_starts` and `tx_ends` arrays

#### 0b. bias.py Module

```python
# src/hulkrna/bias.py

@dataclass
class BiasProfile:
    """Per-transcript positional bias profile with prefix-sum acceleration."""
    prefix_sum: np.ndarray   # float64, shape (L+1,) — B[0]=0, B[i]= Σ b[0..i-1]

    @staticmethod
    def uniform(length: int) -> "BiasProfile":
        """Uniform bias: b[i]=1.0 → prefix_sum = [0,1,2,...,L]."""
        return BiasProfile(prefix_sum=np.arange(length + 1, dtype=np.float64))

    def fragment_weight(self, tx_start: int, tx_end: int) -> float:
        """W(f|t) = (B[tx_end] - B[tx_start]) / (tx_end - tx_start)."""
        span = tx_end - tx_start
        if span <= 0:
            return 0.0
        return (self.prefix_sum[tx_end] - self.prefix_sum[tx_start]) / span

    def effective_length(self, frag_len: int) -> float:
        """L̃_t(l_f) = (1/l_f) Σ_{v=0}^{L-l_f} (B[v+l_f] - B[v])."""
        L = len(self.prefix_sum) - 1  # transcript length
        n_valid = L - frag_len + 1
        if n_valid <= 0:
            return 1.0  # degenerate; avoid log(0)
        B = self.prefix_sum
        eff = (B[frag_len:frag_len + n_valid] - B[:n_valid]).sum() / frag_len
        return max(eff, 1.0)
```

Helper:

```python
def build_uniform_profiles(t_lengths: np.ndarray) -> list[BiasProfile]:
    """Build uniform BiasProfile for each transcript."""
    return [BiasProfile.uniform(int(tl)) for tl in t_lengths]
```

#### 0c. Unit Tests — Uniform Equivalence

Verify for a range of `(t_len, frag_len)` combinations:

```python
profile = BiasProfile.uniform(t_len)
# Fragment weight is always 1.0
assert profile.fragment_weight(s, s + frag_len) == 1.0
# Effective length matches current formula
assert profile.effective_length(frag_len) == max(t_len - frag_len + 1, 1)
# Log-likelihood matches current baked-in correction
log_w = log(profile.fragment_weight(s, s + frag_len))
log_eff = log(profile.effective_length(frag_len))
assert log_w - log_eff ≈ -log(max(t_len - frag_len + 1, 1))
```

### Phase 1 — Wire Bias into EM (future PR)

- During scan, compute `log W(f|t)` and `log L̃_t(l_f)` via `BiasProfile`
- Store `log W − log L̃` instead of current `−log(eff_len)` in `log_liks`
- Verify all 605+ tests still pass (regression gate)
- Verify identical results with uniform profiles (no accuracy change)

### Phase 2 — Learn Bias Profiles (future PR)

- After EM convergence, reweight fragments by posterior assignment probabilities
- Estimate positional bias `b[i]` from the aggregate coverage across transcripts
- Iterate: re-run EM with updated profiles → re-estimate profiles → converge
- Add CLI flag `--bias-model` (none / uniform / learned)

### Phase 3 — Fragment-Length-Dependent Profiles (future PR)

- Extend `BiasProfile` to store per-fragment-length bias arrays
- Stratify the learning step by fragment length bins
- Evaluate on benchmarks for accuracy improvement

## Design Decisions

1. **Per-fragment-length** effective lengths (not mean fragment length) —
   matches current per-candidate correction granularity.
2. **Unified code path** — even the uniform model goes through `BiasProfile`,
   so there's a single code path to maintain.
3. **CSR storage of tx_starts/tx_ends** — these are already computed during scan
   for coverage weight; storing them avoids recomputation during EM/bias learning.
4. **Prefix-sum** for O(1) fragment weight queries — amortizes the cost of
   summing over bias arrays for each fragment.
5. **Phased rollout** — Phase 0 is a pure refactor (no accuracy change),
   Phase 1 wires it in, Phase 2+ learns profiles.
