# Nascent RNA Decoupling: Implementation Plan

## 1. Executive Summary

### Current State

Rigel models nascent RNA (nRNA) as a 1:1 shadow of every mature transcript. The EM component layout per locus is:

```
[0, n_t)        → mRNA (one per transcript)
[n_t, 2*n_t)    → nRNA (one per transcript, nascent shadow)
[2*n_t]         → gDNA (one per locus)
Total: 2*n_t + 1 components
```

This means a gene with 50 isoforms sharing the exact same genomic coordinates spawns 50 *identical* nascent RNA components. This bloats the Equivalence Class (EC) matrix, fragments intronic read evidence across redundant shadows, and relies on the EM prior to enforce symmetry.

### Future State

Nascent RNA is decoupled from individual transcripts. An nRNA is defined strictly by its genomic span: `(ref, strand, genome_start, genome_end)`. All transcripts sharing the same span map to a single unified nRNA component.

```
[0, T)          → mRNA (one per transcript, unchanged)
[T, T+N)        → nRNA (one per unique nascent RNA span)
T+N             → gDNA (one per locus)
Total: T + N + 1 components  (N ≪ T typically)
```

### Benefits

1. **Biological accuracy**: Pre-mRNA is defined by transcriptional start and termination, not splicing. Transcripts sharing the same genomic start/end originate from the same pre-mRNA molecule population.
2. **Computational speed**: Fewer EM components. Equivalence classes containing N identical nRNA shadows collapse to 1, shrinking EC width and accelerating the E-step kernel.
3. **Cleaner output**: A dedicated nRNA abundance file with physically meaningful entities. Downstream tools can map transcript ↔ nRNA ↔ gene without ambiguity.

### Example

```
Transcript T1 + strand: exons [(1000,2000), (5000,5500), (7000,7500), (9000,10000)]
Transcript T2 + strand: exons [(1000,2000), (9000,10000)]
Transcript T3 + strand: exons [(1000,2000), (5000,5500), (9000,10000)]
Transcript T4 + strand: exons [(4500,5500), (9000,10000)]

Current:  4 nRNA shadows (one per transcript), all scored against genomic span
Proposed: 2 unique nRNA spans:
  nRNA_0: (ref, +, 1000, 10000) ← T1, T2, T3
  nRNA_1: (ref, +, 4500, 10000) ← T4
```

---

## 2. Backward Compatibility

**There is NO backward compatibility requirement.** This is a clean break. Old methods, vestigial code paths, and the per-transcript nRNA shadow design will be deleted entirely. The implementation will produce production-quality, clean code.

---

## 3. Phase 1: Index Generation — Define nRNAs at Build Time

### 3.1 Goal

At `rigel index` time, compute the set of unique nascent RNA entities and persist them alongside the transcript index. This establishes the canonical nRNA definition for all downstream processing.

### 3.2 nRNA Extraction Algorithm

During `TranscriptIndex.build()`, after transcripts are parsed and sorted:

1. For each transcript, compute its genomic span: `(ref, strand, start, end)` where `start = min(exon_starts)` and `end = max(exon_ends)` (already the `Transcript.start` and `Transcript.end` properties).
2. Collect all unique `(ref, strand, start, end)` tuples across all transcripts.
3. Assign sequential `nrna_idx` values (0 to N-1) to these unique tuples.
4. For each transcript, record its corresponding `nrna_idx`.

**Uniqueness key**: The nRNA uniqueness key is strictly `(ref, strand, start, end)` — a purely genomic definition. Gene annotations are NOT part of the key because gene definitions are unreliable (overlapping annotations, readthrough transcripts, etc.). Nascent RNA is a physical molecule defined by transcription initiation and termination sites, independent of how genes are annotated. In rare cases, transcripts from *different* genes may share the same genomic span; these map to the same nRNA. The nRNA → gene relationship is therefore many-to-many in the general case (though 1:1 in the overwhelming majority of loci).

**Implementation location**: New function `compute_nrna_table()` in `index.py`, called from `TranscriptIndex.build()`.

```python
def compute_nrna_table(t_df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute unique nRNA entities and transcript-to-nRNA mapping.

    Parameters
    ----------
    t_df : pd.DataFrame
        Transcript table with columns: ref, start, end, strand, t_index.

    Returns
    -------
    nrna_df : pd.DataFrame
        Columns: nrna_idx, ref, strand, start, end, length.
    t_to_nrna : np.ndarray
        int32[n_transcripts] mapping t_index → nrna_idx.
    """
```

### 3.3 New Index File: `nrna.feather`

**Schema**:

| Column      | Type   | Description                                     |
|-------------|--------|-------------------------------------------------|
| `nrna_idx`  | int32  | Global unique identifier (0 to N-1)             |
| `ref`       | string | Chromosome / reference name                      |
| `strand`    | int8   | Strand (1=POS, 2=NEG)                            |
| `start`     | int32  | Genomic start of nRNA span                       |
| `end`       | int32  | Genomic end of nRNA span                         |
| `length`    | int32  | Genomic span length (`end - start`)              |

**Note**: `gene_idx` is intentionally absent from this table. The nRNA → gene relationship is derived from the transcript table (`transcripts.feather`) via `nrna_idx` ↔ `g_index`. Since the relationship is potentially many-to-many (rare but possible with overlapping gene annotations), it cannot be represented as a single column.

This file will be written as `nrna.feather` (with optional TSV mirror) alongside the existing index files.

### 3.4 Updated `transcripts.feather`

Add one column:

| Column      | Type   | Description                                     |
|-------------|--------|-------------------------------------------------|
| `nrna_idx`  | int32  | Foreign key to `nrna.feather`                    |

This requires updating:
- `Transcript.to_dict()` — add `nrna_idx` field
- `Transcript` dataclass — add `nrna_idx: int = -1` slot
- `transcripts_to_dataframe()` — will include `nrna_idx` automatically
- `TranscriptIndex.build()` — compute nRNA table, assign `nrna_idx` to transcripts, write `nrna.feather`

### 3.5 Updated `TranscriptIndex.load()`

Load `nrna.feather` into `self.nrna_df`. Derive:
- `self.t_to_nrna_arr`: int32[T] — fast transcript → nRNA lookup
- `self.num_nrna`: int — number of unique nRNA entities

**Not stored**: `nrna_to_gene` is NOT a stored array because the relationship is many-to-many. When gene information is needed for an nRNA (e.g., output files), it is derived at runtime from the transcript table by collecting all `g_index` values for transcripts sharing that `nrna_idx`.

The inverse mapping (nRNA → transcripts) is NOT stored on disk. It can be reconstructed trivially from `t_to_nrna_arr` at runtime when needed (CSR groupby).

### 3.6 Files Changed

| File | Changes |
|------|---------|
| `index.py` | Add `NRNA_FEATHER`, `NRNA_TSV` constants; add `compute_nrna_table()`; update `build()` to write nRNA table; update `load()` to read nRNA table and populate new arrays; update `_gen_cluster_unambig_intron_intervals()` to tag UNAMBIG_INTRON with `nrna_idx` and merge overlapping intervals per nRNA |
| `transcript.py` | Add `nrna_idx: int = -1` slot to `Transcript`; update `to_dict()` |

### 3.7 Validation

- Unit test: Given a GTF with transcripts sharing/differing in genomic span, verify correct nRNA deduplication.
- Unit test: Verify `t_to_nrna` mapping is correct.
- Index rebuild test: Build an index, load it, verify nRNA table integrity.
- Integration: Existing `test_index.py` tests must pass (with updated assertions for the new column).

---

## 4. Phase 2: Runtime Data Structures

### 4.1 Goal

Update `AbundanceEstimator`, `ScoredFragments`, and `LocusEMInput` to use the `T + N + 1` component layout instead of `2T + 1`.

### 4.2 AbundanceEstimator Changes

Currently the estimator is initialized with arrays of size `num_transcripts` for nRNA tracking:

```python
self.nrna_init = np.zeros(num_transcripts, ...)
self.nrna_em_counts = np.zeros(num_transcripts, ...)
self.nrna_frac_alpha = np.ones(num_transcripts, ...)
self.nrna_frac_beta = np.ones(num_transcripts, ...)
```

These must change:

- `nrna_init`: changes to **per-nRNA** (size N). Intronic evidence is accumulated per-nRNA during scanning (see Section 5 and 7.1). The unambiguous intron intervals are a global genome property — they don't belong to individual transcripts. All transcripts sharing the same nRNA span share the same intronic territory, so intronic evidence is naturally a per-nRNA quantity.
- `nrna_em_counts`: changes to **per-nRNA** (size N). EM assigns counts to the N unified nRNA components.
- `nrna_frac_alpha`, `nrna_frac_beta`: **stay per-transcript** (size T). The nascent fraction $\eta_t$ is a per-transcript property — it represents the fraction of transcript $t$'s total abundance that is in the nascent state.

New arrays needed:
- `self.t_to_nrna`: int32[T] — loaded from index
- `self.num_nrna`: int — number of unique nRNA entities
- `self.nrna_intronic_sense`: float64[N] — replaces `transcript_intronic_sense[T]`. Accumulated per-nRNA during scanning.
- `self.nrna_intronic_antisense`: float64[N] — replaces `transcript_intronic_antisense[T]`. Accumulated per-nRNA during scanning.
- `self.nrna_unspliced_sense`: float64[N] — replaces `transcript_unspliced_sense[T]`. Accumulated per-nRNA during scanning.
- `self.nrna_unspliced_antisense`: float64[N] — replaces `transcript_unspliced_antisense[T]`. Accumulated per-nRNA during scanning.
- `self.nrna_base_index`: changes from `num_transcripts` (currently) to `num_transcripts` (unchanged value, but semantics change — it's now the start of the N nRNA slots, not T nRNA slots)

The key semantic change: `nrna_base_index = T` and there are N nRNA components at positions `[T, T+N)` instead of T nRNA components at `[T, 2T)`.

### 4.3 ScoredFragments Changes

`ScoredFragments.nrna_base_index` stays at `num_transcripts`. But the `t_indices` array will now store:
- `[0, T)` for mRNA candidates
- `[T, T+n]` where `n = nrna_idx` for nRNA candidates (instead of `T + t_idx`)

This is the critical change in the scoring pass: when a fragment generates an nRNA candidate for transcript `t_idx`, the stored component index changes from `nrna_base + t_idx` to `nrna_base + t_to_nrna[t_idx]`.

### 4.4 LocusEMInput Changes

The component layout changes from `[0, n_t) + [n_t, 2*n_t) + [2*n_t]` to `[0, n_t) + [n_t, n_t + n_local_nrna) + [n_t + n_local_nrna]`.

Key changes in `build_locus_em_data()`:
- Compute `local_nrna_set`: the unique nRNA indices referenced by this locus's transcripts.
- `n_local_nrna = len(local_nrna_set)`.
- `n_components = n_t + n_local_nrna + 1`.
- `gdna_idx = n_t + n_local_nrna`.
- Build local mapping arrays for both transcripts (global → local) and nRNAs (global nrna_idx → local nrna position).
- Add `local_to_global_nrna`: int32[n_local_nrna] array mapping local nRNA index to global.
- `bias_profiles` (component lengths):
  - `[0, n_t)` → mRNA spliced exon lengths (unchanged)
  - `[n_t, n_t + n_local_nrna)` → nRNA genomic span length (from nrna table)
  - `[n_t + n_local_nrna]` → gDNA locus span (unchanged)
- `unambig_totals` layout changes to `[n_t + n_local_nrna + 1]`.
- Prior: nRNA prior gating changes — zero nRNA prior when `nrna_init[n] == 0` for that nRNA.
- `nrna_init` for each nRNA component: directly from the per-nRNA `nrna_init[n]` array (accumulated during scanning, see Section 7.1).
- `nrna_frac_alpha`, `nrna_frac_beta`: remain per-transcript (size n_t) and are passed through to the M-step.
- New: `local_t_to_local_nrna`: int32[n_t] mapping each local transcript index to its local nRNA index. Used by the C++ M-step.
- New: `nrna_to_t_offsets`, `nrna_to_t_indices`: CSR mapping from local nRNA → local transcript indices. Used by the C++ M-step for apportioning nRNA counts.

### 4.5 Files Changed

| File | Changes |
|------|---------|
| `estimator.py` | Update `AbundanceEstimator.__init__`, `ScoredFragments`, `LocusEMInput` dataclass |
| `locus.py` | Rewrite `build_locus_em_data()` for new component layout |

---

## 5. Phase 3: Scanner and Scoring

### 5.1 Goal

Update the C++ scoring pass to emit nRNA candidates indexed by nRNA entity rather than by transcript.

### 5.2 Scoring Changes

In `scoring.cpp`, the `FusedScoreBuffer` class currently stores nRNA candidates with index `nrna_base_ + t_idx`. This must change to `nrna_base_ + nrna_idx_for_transcript[t_idx]`.

The scorer needs access to the `t_to_nrna` mapping array. This must be passed into the C++ `FusedScoreBuffer` constructor.

**Key change in nRNA candidate emission** (both single-mapper and multi-mapper paths):

```cpp
// BEFORE:
st.ti[c] = nrna_base_ + t_idx;

// AFTER:
st.ti[c] = nrna_base_ + t_to_nrna_[t_idx];
```

**Equivalence class collapse**: When multiple transcripts share the same nRNA, their intronic fragments will now produce the *same* nRNA component index. The existing WTA (winner-takes-all) merge logic (`merged_nrna` map keyed by `t_idx`) must be rekeyed by `nrna_idx`:

```cpp
// BEFORE: key = t_idx (per transcript)
auto it = merged_nrna.find(t_idx);

// AFTER: key = nrna_idx (per unified nRNA)
int32_t nrna_idx = t_to_nrna_[t_idx];
auto it = merged_nrna.find(nrna_idx);
```

This naturally deduplicates nRNA candidates: if three transcripts all map to the same nRNA, only one candidate (the best-scoring one) is emitted per fragment. This is the source of the EC matrix shrinkage.

### 5.3 nRNA Candidate Scoring Details

After decoupling, nRNA scoring is almost entirely **transcript-independent**. An nRNA entity is a genomic interval `(ref, strand, start, end)` defined at index time — it does not depend on any transcript-level information.

**Coverage weight for nRNA**: Computed against the nRNA's own genomic span `(end - start)`, which is an invariant property of the nRNA entity. This is not derived from any transcript — it comes directly from the `nrna.feather` index table. The scorer should look up coverage weight via `nrna_span_[nrna_idx]` rather than `t_span_[t_idx]`.

**Per-fragment effective length**: Depends on the fragment's aligned extent relative to the nRNA span and the fragment length model, NOT on the transcript. Since the nRNA span is fixed at index time, this calculation is purely a function of the fragment and the nRNA geometry.

**tx_start / tx_end for nRNA candidates**: These represent the fragment's position relative to the nRNA's genomic start. Since the nRNA span is a fixed genomic interval, these are transcript-independent — just `frag_start - nrna_start` and `frag_end - nrna_start`.

**WTA merge**: Within the merge (now keyed by `nrna_idx`), we take the best-scoring candidate. Since all transcripts sharing the same nRNA produce the same nRNA geometry, the scores should be identical (barring minor differences from transcript-specific strand probability). The merge primarily serves to deduplicate: if three transcripts map to the same nRNA, only one nRNA candidate is emitted per fragment.

**Implication for scorer data structures**: The `FusedScoreBuffer` needs access to a per-nRNA span array (from the index) rather than deriving nRNA spans from per-transcript data. This is a new input to the scorer constructor.

### 5.4 Python Scoring Context

In `scan.py`, the `FragmentRouter` constructs the C++ `FusedScoreBuffer`. It must pass the `t_to_nrna` array from the index:

```python
# In FragmentRouter or equivalent setup:
scorer = FusedScoreBuffer(
    ...,
    t_to_nrna=index.t_to_nrna_arr,
    ...
)
```

### 5.5 Connected Components (Locus Builder)

In `build_loci()`, the C++ `connected_components()` function groups transcripts by shared fragments. The behavior is **mostly unchanged** by nRNA decoupling. An intronic fragment already matches all overlapping transcripts in the current code — each overlapping transcript generates an mRNA candidate (from exonic overlap) and/or an nRNA candidate. So transcripts sharing the same genomic region are already connected through their shared fragments.

After decoupling, when a fragment generates an nRNA candidate using `nrna_base + nrna_idx` (instead of `nrna_base + t_idx`), the union-find must map that single nRNA index back to ALL transcripts that share it. This is the nRNA → multi-transcript fan-out.

**Implementation**: Pass the nRNA → transcript CSR arrays (offsets + indices) to `connected_components()`. When an nRNA component `nrna_base + nrna_idx` is encountered, union all transcripts in `nrna_to_t[nrna_idx]`. This handles the edge case of deep-intronic fragments that have no mRNA candidates and only link transcripts through a shared nRNA.

### 5.6 Files Changed

| File | Changes |
|------|---------|
| `native/scoring.cpp` | Add `t_to_nrna_` member to `FusedScoreBuffer`; rekey nRNA WTA merge by `nrna_idx`; emit `nrna_base_ + nrna_idx`; change intronic + unspliced accumulation to per-nRNA; read `nrna_idx` from re-tagged UNAMBIG_INTRON intervals |
| `scan.py` | Pass `t_to_nrna` array to C++ scorer |
| `native/em_solver.cpp` | Update `connected_components()` to handle nRNA → multi-transcript mapping |
| `locus.py` | Update `build_loci()` call to pass nRNA→transcript mapping |

---

## 6. Phase 4: EM Solver — E-Step and M-Step

This is the mathematical core of the change.

### 6.1 Parameters

The optimization parameters are unchanged:

- $\theta_t$: Total abundance of transcript $t$ (mature + nascent combined)
- $\eta_t$: Nascent fraction of transcript $t$ (per-transcript, kept)

But the *component abundances* sent to the E-step change:

- $\alpha_t^{\text{mRNA}} = \theta_t \times (1 - \eta_t)$ — mature component for transcript $t$ (unchanged)
- $\alpha_n^{\text{nRNA}} = \sum_{t \in \mathcal{T}(n)} \theta_t \times \eta_t$ — nascent component for nRNA $n$ (NEW: aggregated across transcripts)

where $\mathcal{T}(n)$ is the set of transcripts mapping to nRNA $n$.

### 6.2 E-Step

The E-step operates on the `T + N + 1` component space (reduced from `2T + 1`). The math is standard — compute posteriors using log-weights:

$$\log w_t^{\text{mRNA}} = \log(\theta_t) + \log(1 - \eta_t) - \log(\text{eff\_len}_t^{\text{mRNA}})$$

$$\log w_n^{\text{nRNA}} = \log\left(\sum_{t \in \mathcal{T}(n)} \theta_t \cdot \eta_t\right) - \log(\text{eff\_len}_n^{\text{nRNA}})$$

$$\log w_{\text{gDNA}} = \log(\gamma) - \log(\text{eff\_len}_{\text{gDNA}})$$

The E-step yields expected counts:
- $E_{\text{mRNA}}[t]$ for each transcript $t$
- $E_{\text{nRNA}}[n]$ for each nRNA $n$
- $E_{\text{gDNA}}$

### 6.3 M-Step (The Key Change)

The M-step must apportion the nRNA expected counts $E_{\text{nRNA}}[n]$ back to individual transcripts proportional to their nascent mass share.

For transcript $t$ mapping to nRNA $n$:

1. **Compute transcript's share of nRNA mass:**

$$W_t = \frac{\theta_t^{(\text{old})} \cdot \eta_t^{(\text{old})}}{\alpha_n^{(\text{old})}}$$

where $\alpha_n^{(\text{old})} = \sum_{t' \in \mathcal{T}(n)} \theta_{t'}^{(\text{old})} \cdot \eta_{t'}^{(\text{old})}$.

If $\alpha_n^{(\text{old})} = 0$, distribute uniformly: $W_t = 1 / |\mathcal{T}(n)|$.

2. **Apportion nRNA counts:**

$$E_{\text{nRNA}}^{\text{apportioned}}[t] = E_{\text{nRNA}}[n] \times W_t$$

3. **Update transcript parameters:**

$$\theta_t^{(\text{new})} \propto E_{\text{mRNA}}[t] + E_{\text{nRNA}}^{\text{apportioned}}[t] + \text{unambig}[t] + \text{prior}[t]$$

4. **Update nascent fraction (MAP Beta, unchanged formula):**

$$\eta_t^{(\text{new})} = \frac{E_{\text{nRNA}}^{\text{apportioned}}[t] + \alpha_t - 1}{E_{\text{mRNA}}[t] + E_{\text{nRNA}}^{\text{apportioned}}[t] + \alpha_t + \beta_t - 2}$$

### 6.4 SQUAREM Acceleration

Currently, SQUAREM operates on `theta_t` (size `n_t + 1`: T transcript abundances + 1 gDNA). This remains the same. The nRNA component values are *derived* from theta_t and eta during each EM step — they are NOT extrapolated by SQUAREM (this is the existing design and it is correct).

### 6.5 C++ Implementation: `linked_map_em_step()`

The current function signature is:

```cpp
static void linked_map_em_step(
    const double* theta_t,       // [n_t + 1]: θ_1..θ_T, γ
    const double* nrna_frac,     // [n_t]
    const EmEquivClass& ec_data,
    const double* log_eff_len,   // [2*n_t + 1]
    ...
    int n_transcripts,
    int n_components)            // 2*n_t + 1
```

It changes to:

```cpp
static void linked_map_em_step(
    const double* theta_t,        // [n_t + 1]: θ_1..θ_T, γ
    const double* nrna_frac,      // [n_t]
    const EmEquivClass& ec_data,
    const double* log_eff_len,    // [n_t + n_nrna + 1]
    ...
    int n_transcripts,
    int n_nrna,                   // NEW: number of nRNA components
    int n_components,             // n_t + n_nrna + 1
    const int32_t* t_to_nrna,    // [n_t]: local transcript → local nRNA
    const int32_t* nrna_to_t_off, // [n_nrna + 1]: CSR offsets
    const int32_t* nrna_to_t_idx) // [sum]: CSR flat transcript indices
```

**Building log_weights:**

```cpp
// Compute nRNA component abundances by summing across transcripts
std::vector<double> nrna_abundance(n_nrna, 0.0);
for (int i = 0; i < n_t; ++i) {
    nrna_abundance[t_to_nrna[i]] += theta_t[i] * nrna_frac[i];
}

// mRNA components
for (int i = 0; i < n_t; ++i) {
    log_weights[i] = log(theta_t[i]) + log(1 - nrna_frac[i]) - log_eff_len[i];
}
// nRNA components
for (int n = 0; n < n_nrna; ++n) {
    log_weights[n_t + n] = log(nrna_abundance[n]) - log_eff_len[n_t + n];
}
// gDNA
log_weights[n_t + n_nrna] = log(theta_t[n_t]) - log_eff_len[n_t + n_nrna];
```

**M-step apportionment:**

```cpp
// After E-step yields em_totals[n_t + n_nrna + 1]

for (int n = 0; n < n_nrna; ++n) {
    double nrna_count = em_totals[n_t + n] + unambig_totals[n_t + n];
    double total_nrna_mass = nrna_abundance[n];  // from earlier

    // Apportion to transcripts
    int start = nrna_to_t_off[n];
    int end   = nrna_to_t_off[n + 1];
    for (int j = start; j < end; ++j) {
        int t = nrna_to_t_idx[j];
        double share = (total_nrna_mass > 0)
            ? (theta_t[t] * nrna_frac[t]) / total_nrna_mass
            : 1.0 / (end - start);
        double t_nrna_count = nrna_count * share;
        double t_mrna_count = unambig_totals[t] + em_totals[t];
        double total_t = t_mrna_count + t_nrna_count;

        theta_t_new[t] = total_t + prior[t] + prior[n_t + n] * share;
        // (normalization happens after all transcripts + gDNA)

        // Update nrna_frac via MAP Beta
        double nrna_frac_den = total_t + alpha[t] + beta[t] - 2.0;
        nrna_frac_new[t] = clamp((t_nrna_count + alpha[t] - 1.0) / nrna_frac_den);
    }
}
```

### 6.6 `linked_run_squarem()` Changes

- State size changes: still `n_t + 1` (theta_t + gamma). Unchanged.
- All inner calls to `linked_map_em_step()` pass the new nRNA mapping arrays.
- Decomposition of theta into full component vector: 
  - Before: `theta[n_t + i] = state0[i] * nrna_frac[i]`
  - After: `theta[n_t + n] = sum(state0[t] * nrna_frac[t] for t in nrna_to_t(n))`
- `alpha_out` size changes from `2*n_t + 1` to `n_t + n_nrna + 1`.
- Post-EM pruning operates at transcript level (unchanged). When a transcript is pruned, its nRNA contribution is zeroed. If ALL transcripts for an nRNA are pruned, the nRNA component naturally zeroes out.

### 6.7 `run_locus_em_native()` Changes

The nanobind entry point must accept the new mapping arrays:

```cpp
run_locus_em_native(
    ...,
    int n_transcripts,
    int n_nrna,                              // NEW
    i32_1d t_to_nrna,                        // NEW: [n_t]
    i32_1d nrna_to_t_offsets,                // NEW: [n_nrna + 1]
    i32_1d nrna_to_t_indices,                // NEW: [sum]
    f64_1d nrna_frac_alpha_arr,
    f64_1d nrna_frac_beta_arr)
```

### 6.8 `batch_locus_em()` Changes

The batched C++ entry point must propagate the nRNA mapping for each locus. Each locus has its own local nRNA subset, so the batch function must build per-locus local mappings from the global `t_to_nrna` array.

**Option A**: Build per-locus nRNA mappings in Python, pass flat CSR.
**Option B**: Build per-locus nRNA mappings inside C++.

**Recommendation**: Option B (C++ construction) for performance — the batch function already builds per-locus CSR data internally, and adding nRNA mapping construction is natural.

New parameters for `batch_locus_em`:
- `global_t_to_nrna`: int32[T] — global transcript → global nRNA mapping
- `num_global_nrna`: int — total nRNA count (for bounds checking)

Inside the per-locus loop, the C++ code will:
1. Extract the locus's transcript indices.
2. Compute the unique local nRNA set from `global_t_to_nrna[t_arr]`.
3. Build `local_t_to_local_nrna`, `nrna_to_t_offsets`, `nrna_to_t_indices`.
4. Build `n_components = n_t + n_local_nrna + 1`.

### 6.9 Files Changed

| File | Changes |
|------|---------|
| `native/em_solver.cpp` | Rewrite `linked_map_em_step()`, `linked_run_squarem()`, `run_locus_em_native()`, `batch_locus_em()`; update equivalence class builder for new component count; add nRNA aggregation/apportionment logic |
| `estimator.py` | Update `run_batch_locus_em()` to pass nRNA mapping arrays; update `LocusEMInput` fields |
| `locus.py` | Update `build_locus_em_data()` for new component layout |

---

## 7. Phase 5: nRNA Initialization and Prior Computation

### 7.1 nRNA Initialization

nRNA initialization is a **per-nRNA** quantity, not per-transcript. Unambiguous intron intervals (UNAMBIG_INTRON) are a global genome property: a genomic position is unambiguously intronic if it falls within an intron of at least one transcript and does NOT overlap any exon from any transcript in the cluster. These intervals don't conceptually "belong" to individual transcripts.

**Scanning changes**: During the scoring pass, intronic evidence (sense and antisense fragment counts) is accumulated **per-nRNA** rather than per-transcript. In `scoring.cpp`, the accumulation at lines ~878-879 changes from:

```cpp
// BEFORE (per-transcript):
if (has_ui) {
    if (is_anti) acc_ia[t_idx] += weight;
    else         acc_is[t_idx] += weight;
}

// AFTER (per-nRNA, nrna_idx read directly from interval):
if (has_ui) {
    int32_t nrna = nrna_ind[k];  // nrna_idx from re-tagged interval
    if (is_anti) acc_ia[nrna] += weight;
    else         acc_is[nrna] += weight;
}
```

The accumulator arrays `intronic_sense` and `intronic_antisense` change from size T to size N (num_nrna). The Python side (`scan.py`, `estimator.py`) passes `estimator.nrna_intronic_sense` and `estimator.nrna_intronic_antisense` (both size N) instead of the old per-transcript arrays.

**UNAMBIG_INTRON re-tagging**: The UNAMBIG_INTRON intervals in the index are re-tagged with `nrna_idx` (instead of `t_index`) at index build time (see Section 14.2). `_gen_cluster_unambig_intron_intervals()` groups intervals by `nrna_idx` and merges overlapping pieces. This eliminates redundant intervals for transcripts sharing the same nRNA and means the scanner reads `nrna_idx` directly from the interval — no `t_to_nrna` lookup needed during intronic accumulation.

**`compute_nrna_init()` changes**: This function now takes the per-nRNA intronic sense/antisense arrays (size N) and returns per-nRNA init values (size N):

$$\text{nrna\_init}[n] = f(\text{intronic\_sense}[n], \text{intronic\_antisense}[n])$$

where $f$ is the existing strand-excess formula.

### 7.2 nRNA Prior Gating

Currently, nRNA prior is zeroed when:
1. The transcript is single-exon (`transcript_span ≤ exonic_length`)
2. `nrna_init[t] == 0` for that transcript

After decoupling, the nRNA component prior should be zeroed when ALL transcripts mapping to it either:
1. Are single-exon (no introns), OR
2. The per-nRNA `nrna_init[n] == 0` (no intronic evidence for this nRNA)

Since nRNA init is now per-nRNA, this check is direct: if `nrna_init[n] == 0`, gate the prior for nRNA component `n`.

### 7.3 nrna_frac Prior Computation

The hierarchical nRNA fraction prior (`compute_hybrid_nrna_frac_priors()`) remains per-transcript. The hierarchy (global → locus-strand → TSS-group → transcript) is unchanged. Each transcript still has its own $\eta_t$ with Beta prior parameters $(\alpha_t, \beta_t)$.

**Interaction with per-nRNA intronic evidence**: The current prior computation uses per-transcript `transcript_intronic_sense/antisense` as part of the evidence feeding the hierarchy. With intronic evidence now accumulated per-nRNA, the prior function must read `nrna_intronic_sense[t_to_nrna[t]]` / `nrna_intronic_antisense[t_to_nrna[t]]` to get each transcript's intronic evidence (which is shared among all transcripts mapping to the same nRNA). The exonic evidence used in the prior (`transcript_exonic_sense/antisense`) remains per-transcript since exon structures differ across isoforms.

### 7.4 Files Changed

| File | Changes |
|------|---------|
| `locus.py` | Update `build_locus_em_data()` for nRNA prior gating using per-nRNA init |
| `estimator.py` | Replace `transcript_intronic_sense/antisense[T]` with `nrna_intronic_sense/antisense[N]`; update `compute_nrna_init()` for per-nRNA |
| `scan.py` | Pass per-nRNA intronic accumulator arrays to C++ scorer |
| `scoring.cpp` | Change intronic + unspliced accumulation to per-nRNA; read `nrna_idx` from re-tagged UNAMBIG_INTRON intervals |

---

## 8. Phase 6: Output Generation

### 8.1 New File: `nrna_quant.feather` / `nrna_quant.tsv`

A dedicated nascent RNA abundance output file.

**Schema**:

| Column            | Type    | Description                                              |
|-------------------|---------|----------------------------------------------------------|
| `nrna_idx`        | int32   | Unique nRNA index (matches index `nrna.feather`)         |
| `gene_id`         | string  | Gene identifier(s) — derived from transcript table; comma-delimited if multiple genes share this nRNA |
| `gene_name`       | string  | Gene name(s) — same derivation as `gene_id`              |
| `ref`             | string  | Chromosome                                               |
| `strand`          | string  | Strand (+/-)                                             |
| `start`           | int32   | Genomic start                                            |
| `end`             | int32   | Genomic end                                              |
| `effective_length` | float64 | Effective length (span - mean_frag + 1)                  |
| `est_counts`      | float64 | EM-assigned nascent RNA counts: $E_{\text{nRNA}}[n]$     |
| `tpm`             | float64 | Transcripts per million (nRNA-based, using nRNA eff len) |

**Note on `gene_id` / `gene_name`**: These are convenience columns derived at output time by collecting all distinct gene IDs among transcripts sharing this `nrna_idx`. In the vast majority of cases this is a single gene. For the rare case of overlapping gene annotations, values are comma-delimited.

**Implementation**: New method `AbundanceEstimator.get_nrna_counts_df(index)`.

### 8.2 Updated `quant.feather` (Transcript Abundance)

**Removed columns**: `nrna` (the per-transcript EM-assigned nascent count is no longer directly meaningful as a standalone column since nRNA is now a shared component)

**Added columns**:

| Column          | Type    | Description                                             |
|-----------------|---------|----------------------------------------------------------|
| `nrna_idx`      | int32   | Foreign key to `nrna_quant` — which nRNA produces this transcript |
| `nrna_fraction` | float64 | $\eta_t$ — nascent RNA fraction (0.0 to 1.0) from EM    |

**Retained** (with semantics clarified):
- `mrna`: mRNA counts = total assignment to the mature component
- `rna_total`: mrna + apportioned nascent counts for this transcript

### 8.3 Updated `gene_quant.feather` (Gene Abundance)

**Added columns**:

| Column             | Type    | Description                                          |
|--------------------|---------|------------------------------------------------------|
| `est_nascent_counts` | float64 | $\sum_{n \in G} E_{\text{nRNA}}[n]$ for all nRNAs in gene |
| `nascent_fraction`   | float64 | `est_nascent_counts / (mrna + est_nascent_counts)`   |

### 8.4 CLI Changes

Update `_write_quant_outputs()` in `cli.py` to write the new `nrna_quant.feather` (and `.tsv` mirror). Update `summary.json` to include nRNA-level statistics.

### 8.5 Files Changed

| File | Changes |
|------|---------|
| `estimator.py` | Add `get_nrna_counts_df()`; update `get_counts_df()` and `get_gene_counts_df()` |
| `cli.py` | Write `nrna_quant.feather`/`.tsv`; update summary.json |

---

## 9. Phase 7: Cleanup and Dead Code Removal

### 9.1 Remove Per-Transcript nRNA Shadow Vestiges

Delete all code that assumes `n_components = 2*n_t + 1`:
- Remove `nrna_base_index = num_transcripts` as a concept of "nRNA shadow per transcript"
- Clean up comments referencing the old layout
- Remove any `[n_t:2*n_t]` indexing patterns for nRNA

### 9.2 Update Tests

All tests referencing the `2*n_t + 1` layout must be updated. Key test files:
- `test_em_impl.py` — EM algorithm tests
- `test_estimator.py` — estimator output tests
- `test_pipeline_routing.py` — pipeline integration tests
- `test_buffer.py` — buffer/scoring tests
- `test_golden_output.py` — golden output comparison (will need new golden files)
- `test_bias.py` — bias correction tests

### 9.3 Documentation

Update:
- `docs/METHODS.md` — mathematical description of the EM
- `docs/parameters.md` — parameter descriptions
- `docs/MANUAL.md` — user-facing documentation
- CLI help strings

---

## 10. Implementation Sequence

### Stage 1: Index (Phase 1)

**Goal**: Define nRNAs at index time. Zero runtime impact.

1. Add `nrna_idx` slot to `Transcript` dataclass.
2. Implement `compute_nrna_table()` in `index.py`.
3. Update `TranscriptIndex.build()` to compute and write `nrna.feather`.
4. Update `_gen_cluster_unambig_intron_intervals()` to tag UNAMBIG_INTRON with `nrna_idx` and merge overlapping intervals per nRNA.
5. Update `TranscriptIndex.load()` to read `nrna.feather` and populate `t_to_nrna_arr`, `nrna_df`, `num_nrna`.
6. Write unit tests for nRNA extraction and index round-trip.

**Validation**: All existing tests pass. New tests verify nRNA deduplication.

### Stage 2: Scoring (Phase 3)

**Goal**: Emit nRNA candidates indexed by nRNA entity. EM still runs on old layout (temporary compatibility shim).

1. Pass `t_to_nrna` to C++ `FusedScoreBuffer`.
2. Change nRNA candidate emission to use `nrna_base + t_to_nrna[t_idx]`.
3. Rekey nRNA WTA merge by `nrna_idx`.
4. Change intronic + unspliced evidence accumulation to per-nRNA (read `nrna_idx` from re-tagged intervals).
5. Update `connected_components()` for nRNA → multi-transcript mapping.
6. Temporarily adapt `build_locus_em_data()` to handle the new nRNA indexing while still building the old `2*n_t + 1` layout (bridge code).

**Validation**: Scoring tests pass. Locus construction still works. EM results may differ slightly due to EC collapse but should be directionally correct.

### Stage 3: EM Solver (Phases 2, 4, 5)

**Goal**: Full `T + N + 1` component layout in the EM.

1. Update `LocusEMInput` for new component layout.
2. Rewrite `build_locus_em_data()` to build `n_t + n_local_nrna + 1` components with nRNA mapping CSR.
3. Rewrite C++ `linked_map_em_step()` with nRNA aggregation and M-step apportionment.
4. Rewrite C++ `linked_run_squarem()`.
5. Update `run_locus_em_native()` and `batch_locus_em()`.
6. Remove bridge code from Stage 2.

**Validation**: EM convergence verified. Total counts sum correctly: `sum(mRNA) + sum(nRNA) + gDNA == total_assigned`. EC matrix sizes verified to be smaller.

### Stage 4: Output (Phase 6)

**Goal**: New output schema.

1. Implement `get_nrna_counts_df()`.
2. Update `get_counts_df()` — remove `nrna` column, add `nrna_idx` and `nrna_fraction`.
3. Update `get_gene_counts_df()` — add `est_nascent_counts` and `nascent_fraction`.
4. Update CLI to write `nrna_quant.feather`/`.tsv`.
5. Update `summary.json`.

**Validation**: Output file schemas verified. TPM sums to 1M. Nascent counts are consistent between nRNA file and transcript file.

### Stage 5: Cleanup (Phase 7)

**Goal**: Remove all dead code and update tests/docs.

1. Delete old `2*n_t + 1` code paths.
2. Update all tests.
3. Regenerate golden output files.
4. Update documentation.

**Validation**: Full test suite green. No vestigial code.

---

## 11. Detailed File Change Matrix

| File | Stage | Nature of Change |
|------|-------|-----------------|
| `transcript.py` | 1 | Add `nrna_idx` slot |
| `index.py` | 1 | `compute_nrna_table()`, updated build/load, UNAMBIG_INTRON re-tagging |
| `native/scoring.cpp` | 2 | nRNA candidate indexing, WTA rekey, intronic accum per-nRNA |
| `scan.py` | 2 | Pass `t_to_nrna` to C++ |
| `native/em_solver.cpp` | 2,3 | `connected_components()`, then full EM rewrite |
| `locus.py` | 2,3 | `build_loci()` update, then `build_locus_em_data()` rewrite |
| `estimator.py` | 3,4 | `LocusEMInput`, `AbundanceEstimator`, output methods |
| `cli.py` | 4 | New nRNA output file |
| `config.py` | — | No changes expected |
| `pipeline.py` | 3 | Pass nRNA arrays through quant pipeline |
| `tests/*` | 5 | Update assertions, golden files |
| `docs/*` | 5 | Update methodology docs |

---

## 12. Risks and Mitigations

### Risk: Numerical stability of M-step apportionment

When $\alpha_n^{(\text{old})} \approx 0$, the share $W_t$ is ill-defined.

**Mitigation**: Fall back to uniform distribution ($W_t = 1/|\mathcal{T}(n)|$) when the denominator is below a small epsilon.

### Risk: SQUAREM divergence with new parameterization

The extrapolation operates on theta_t, which is unchanged. nrna_frac is updated via the M-step only, not extrapolated. This is the existing safe design.

**Mitigation**: None needed — the existing SQUAREM structure handles this correctly.

### Risk: Post-EM pruning logic

Pruning currently zeroes transcript pairs (mRNA + nRNA shadow). With decoupled nRNA, pruning a transcript should reduce its share of the shared nRNA but not zero the entire nRNA component.

**Mitigation**: Pruning operates on transcript theta_t. A pruned transcript's nrna_frac contribution naturally drops to zero. The nRNA component survives as long as any contributing transcript survives.

### Risk: Edge case — all transcripts of an nRNA are single-exon

If all transcripts mapping to an nRNA are single-exon, none has intronic evidence. The nRNA component should have zero prior and be effectively disabled.

**Mitigation**: Prior gating already handles this — nRNA prior is zeroed when all constituent transcripts have zero nrna_init.

### Risk: nRNA → gene many-to-many relationship

With `(ref, strand, start, end)` as the uniqueness key (no gene_idx), transcripts from different genes sharing the exact same genomic span map to the same nRNA. This creates a many-to-many nRNA → gene relationship.

**Mitigation**: This is rare in practice (requires different gene annotations on the same strand with identical start/end). Output files handle it via comma-delimited gene_id/gene_name. The EM is unaffected — it operates on transcript/nRNA components, not genes.

### Risk: UNAMBIG_INTRON interval deduplication (RESOLVED)

UNAMBIG_INTRON intervals are re-tagged with `nrna_idx` and merged/deduplicated at index build time (see Section 14.2). Double-counting is eliminated by construction.

### Risk: nRNA with a single transcript (degenerate case)

When an nRNA maps to exactly one transcript, the behavior should be mathematically identical to the old per-transcript shadow design (modulo float precision).

**Mitigation**: Verify with unit tests that single-transcript nRNAs produce identical results to the old code.

---

## 13. Verification Criteria

1. **Counts conservation**: `sum(mRNA[t]) + sum(nRNA[n]) + gDNA == total_ambiguous_fragments + total_unambig_fragments` per locus.
2. **EC collapse**: On a test locus with K transcripts sharing 1 genomic span, verify that ECs max out at width `K + 2` (K mRNA, 1 nRNA, 1 gDNA) instead of `2K + 1`.
3. **Degenerate equivalence**: For a single-transcript nRNA, EM results match the old per-transcript shadow within floating-point tolerance.
4. **Index round-trip**: Build index, load it, verify `nrna_df` and `t_to_nrna_arr` are consistent.
5. **TPM consistency**: nRNA TPM sums correctly. Transcript TPM sums correctly.
6. **Regression on real data**: Run on a standard test dataset and compare total mRNA/nRNA/gDNA fractions — should be similar (not identical due to EC collapse effects).

---

## 14. Resolved Design Decisions

The following decisions were resolved during plan review.

### 14.1 Unspliced accumulators → per-nRNA

**Decision**: Change `transcript_unspliced_sense/antisense[T]` to `nrna_unspliced_sense/antisense[N]` (per-nRNA). Unspliced reads within a genomic span are a per-nRNA property. This is conceptually simpler and consistent with the per-nRNA intronic accumulators. The gDNA EB prior aggregates to per-locus anyway, so the per-nRNA granularity is sufficient. In `scoring.cpp`, the unspliced accumulation changes from `acc_us[t_idx]` / `acc_ua[t_idx]` to `acc_us[t_to_nrna[t_idx]]` / `acc_ua[t_to_nrna[t_idx]]` (alongside the intronic accumulation change). The Python-side arrays in `estimator.py` change accordingly.

### 14.2 UNAMBIG_INTRON re-tagged with nrna_idx at index time

**Decision**: Re-tag UNAMBIG_INTRON intervals with `nrna_idx` (instead of `t_index`) at index build time. This:
- Eliminates redundant overlapping intervals for transcripts sharing the same nRNA
- Simplifies runtime accumulation — no `t_to_nrna` lookup needed during scanning; the scorer reads `nrna_idx` directly from the interval
- Requires updating `_gen_cluster_unambig_intron_intervals()` to group by `nrna_idx` and merge overlapping intervals within each nRNA

The `intervals.feather` schema changes: the `t_index` column on UNAMBIG_INTRON rows becomes `nrna_idx`. If the downstream code ever needs to map back to transcripts, it uses the separate `nrna → transcript` CSR mapping. This also resolves the UNAMBIG_INTRON deduplication risk (Section 12) — double-counting is eliminated by construction.

### 14.3 Both exonic and intronic unspliced fragments contribute to nRNA evidence

**Decision**: Confirmed. Both exonic unspliced and intronic unspliced fragments contribute evidence for nascent RNA. The nrna_frac prior computation (`_compute_hybrid_nrna_frac_vec()`) should use per-nRNA intronic evidence (`nrna_intronic_sense/antisense[n]`) alongside per-nRNA unspliced evidence. Transcripts sharing the same nRNA will have identical intronic and unspliced evidence. The exonic *spliced* evidence remains per-transcript (different exon structures produce different splicing patterns).

### 14.4 Exonic unspliced fragments are truly ambiguous (mRNA / nRNA / gDNA)

**Decision**: Unspliced reads in exonic regions (the `!has_ui` branch in `scoring.cpp`) are genuinely ambiguous — they could originate from mature mRNA, nascent RNA, or gDNA. They must NOT be treated as evidence for any single component. These fragments are correctly handled by the EM as ambiguous candidates with mRNA, nRNA, and gDNA components in their equivalence class. The exonic sense/antisense accumulators (`acc_es`, `acc_ea`) remain per-transcript since exon structures differ across isoforms, and they are used only for the nrna_frac prior hierarchy (not for nRNA initialization).
