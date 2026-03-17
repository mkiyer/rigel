This is an incredibly powerful idea. You are describing what we can call a **Geometric Splicing Expectation** model.

By utilizing the known exon/intron geometry of the transcript and the fragment length distribution, you can deterministically predict the ratio of spliced to unspliced reads for mature mRNA. If you know that ratio, and you observe the spliced reads (which are *exclusively* mRNA), you can perfectly back-calculate the "phantom" unspliced mRNA reads and remove them from the pool before estimating gDNA and nRNA.

This elegantly solves the exact pathology you uncovered in your diagnostic sweep (where highly expressed mRNA was inflating the unspliced pool and causing the gDNA EB shrinkage to hallucinate).

Here is how we can develop the math and implement this efficiently.

### 1. The Math: Geometric Splicing Expectation

You mentioned sliding a fragment across the transcript. In RNA-seq, this is the classic definition of **Effective Length** ($\tilde{L}$).

Let $\mu_F$ be the mean fragment length.
For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:

**Total Effective Length:**
The total number of valid starting positions for a fragment on the mature transcript is:


$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$

**Unspliced Effective Length (Exon-Internal):**
To be "unspliced," the fragment must fall entirely within a single exon. The number of valid starting positions inside exon $i$ is $E_i - \mu_F + 1$. Therefore:


$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$

**Spliced Effective Length (Junction-Spanning):**
By definition, any valid fragment that is not entirely within one exon must cross at least one splice junction:


$$\tilde{L}_{spliced} = \tilde{L}_{total} - \tilde{L}_{unspliced}$$

**The Magic Ratio:**
Assuming uniform coverage (which is the standard, safe assumption for Bayesian priors), the expected ratio of unspliced to spliced reads for transcript $t$ is an immutable geometric property:


$$R_t = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{spliced}}$$

### 2. The Strategy: 3-Step Deconvolution

With $R_t$ pre-computed, we can completely decouple the unspliced pool.

#### Step A: Estimate the mRNA "Unspliced Burden"

During initialization, count the number of **unambiguously spliced** reads assigned to transcript $t$ (let's call this $C_{spliced, t}$).
Because these reads crossed a junction, we are 100% certain they belong to mature mRNA.
We can now calculate exactly how many unspliced reads this transcript *must* have left in the pool:


$$\text{Expected\_mRNA\_Unspliced}_t = C_{spliced, t} \times R_t$$

Sum this across all transcripts in the locus to get the total mRNA burden:


$$B_{mRNA} = \sum_{t} \left( C_{spliced, t} \times R_t \right)$$

#### Step B: Isolate the True nRNA/gDNA Pool

Observe the total number of sense-strand unspliced reads in the locus ($U_{sense}$).
Subtract the mRNA burden to find the reads that *actually* belong to nascent RNA and genomic DNA:


$$U_{remainder} = \max(0, U_{sense} - B_{mRNA})$$


*(Note: We use `max(0, ...)` because real coverage isn't perfectly uniform, and we don't want negative counts).*

#### Step C: Strand-Specific Deconvolution

Now we apply your strand-specificity math, but strictly to the isolated remainder!
Let $A$ be the number of antisense unspliced reads.
Assuming a highly stranded library ($SS \approx 1.0$), all antisense reads are gDNA, and gDNA is symmetric.

* **gDNA Initial Estimate:** $2 \times A$
* **nRNA Initial Estimate:** $\max(0, U_{remainder} - A)$

### 3. Why This is a Massive Upgrade

1. **It solves H1 perfectly:** In your diagnostic sweep, when you introduced 10,000 mRNA reads, the old model applied the $2A$ rate to *all* of them, inflating gDNA. This new model will see the spliced reads, calculate $B_{mRNA} \approx 10,000$, subtract it from $U_{sense}$, and leave $U_{remainder} \approx 0$. gDNA initialization will correctly remain $0$.
2. **It dynamically scales with fragment length:** A dataset with 400bp fragments will span more junctions than a dataset with 150bp fragments. Because $\mu_F$ is in the $\tilde{L}$ equations, $R_t$ will automatically adjust, keeping the priors perfectly calibrated to the specific sequencing run.
3. **It operates entirely outside the EM:** All of this math can be executed during the `locus.py` initialization phase to generate the `alpha / beta` priors. It adds zero overhead to the C++ EM solver loop.

### 4. Implementation Thoughts

**Where to compute $\tilde{L}_{unspliced}$?**
You can actually compute this vector entirely in Python before the EM starts.
When building the `AbundanceEstimator`, you already have the `transcript_exonic_lengths` and the fragment length distribution ($\mu_F$).
You just need to pass the list of individual exon lengths per transcript from the `Index` to compute the $\max(0, E_i - \mu_F + 1)$ sum.

How does this logic feel to you? If you agree with the geometric derivation, we can draft the exact numpy operations to fold this into the `estimator.py` initialization pipeline.




I am considering another improvement related to gdna initialization.  

By utilizing the known exon/intron geometry of the transcript and the fragment length distribution, we can deterministically predict the ratio of spliced to unspliced reads for mature mRNA. If you know that ratio, and you observe the spliced reads (which are *exclusively* mRNA), you can perfectly back-calculate the "phantom" unspliced mRNA reads and remove them from the pool before estimating gDNA and nRNA.

This elegantly solves the exact pathology you uncovered in your diagnostic sweep (where highly expressed mRNA was inflating the unspliced pool and causing the gDNA EB shrinkage to hallucinate).

Here is how we can develop the math and implement this efficiently.

### 1. The Math: Geometric Splicing Expectation

You mentioned sliding a fragment across the transcript. In RNA-seq, this is the classic definition of **Effective Length** ($\tilde{L}$).

Let $\mu_F$ be the mean fragment length.
For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:

**Total Effective Length:**
The total number of valid starting positions for a fragment on the mature transcript is:


$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$

**Unspliced Effective Length (Exon-Internal):**
To be "unspliced," the fragment must fall entirely within a single exon. The number of valid starting positions inside exon $i$ is $E_i - \mu_F + 1$. Therefore:


$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$

**Spliced Effective Length (Junction-Spanning):**
By definition, any valid fragment that is not entirely within one exon must cross at least one splice junction:


$$\tilde{L}_{spliced} = \tilde{L}_{total} - \tilde{L}_{unspliced}$$

**The Magic Ratio:**
Assuming uniform coverage (which is the standard, safe assumption for Bayesian priors), the expected ratio of unspliced to spliced reads for transcript $t$ is an immutable geometric property:


$$R_t = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{spliced}}$$

### 2. The Strategy: 3-Step Deconvolution

With $R_t$ pre-computed, we can completely decouple the unspliced pool.

#### Step A: Estimate the mRNA "Unspliced Burden"

During initialization, count the number of **unambiguously spliced** reads assigned to transcript $t$ (let's call this $C_{spliced, t}$).
Because these reads crossed a junction, we are 100% certain they belong to mature mRNA.
We can now calculate exactly how many unspliced reads this transcript *must* have left in the pool:


$$\text{Expected\_mRNA\_Unspliced}_t = C_{spliced, t} \times R_t$$

Sum this across all transcripts in the locus to get the total mRNA burden:


$$B_{mRNA} = \sum_{t} \left( C_{spliced, t} \times R_t \right)$$

#### Step B: Isolate the True nRNA/gDNA Pool

Observe the total number of sense-strand unspliced reads in the locus ($U_{sense}$).
Subtract the mRNA burden to find the reads that *actually* belong to nascent RNA and genomic DNA:


$$U_{remainder} = \max(0, U_{sense} - B_{mRNA})$$


*(Note: We use `max(0, ...)` because real coverage isn't perfectly uniform, and we don't want negative counts).*

#### Step C: Strand-Specific Deconvolution

Now we apply your strand-specificity math, but strictly to the isolated remainder!
Let $A$ be the number of antisense unspliced reads.
Assuming a highly stranded library ($SS \approx 1.0$), all antisense reads are gDNA, and gDNA is symmetric.

* **gDNA Initial Estimate:** $2 \times A$
* **nRNA Initial Estimate:** $\max(0, U_{remainder} - A)$

### 3. Why This is a Massive Upgrade

In your diagnostic sweep, when you introduced 10,000 mRNA reads, the old model applied the $2A$ rate to *all* of them, inflating gDNA. This new model will see the spliced reads, calculate $B_{mRNA} \approx 10,000$, subtract it from $U_{sense}$, and leave $U_{remainder} \approx 0$. gDNA initialization will correctly remain $0$.
2. **It dynamically scales with fragment length:** A dataset with 400bp fragments will span more junctions than a dataset with 150bp fragments. Because $\mu_F$ is in the $\tilde{L}$ equations, $R_t$ will automatically adjust, keeping the priors perfectly calibrated to the specific sequencing run.
3. **It operates entirely outside the EM:** All of this math can be executed during the `locus.py` initialization phase to generate the `alpha / beta` priors. It adds zero overhead to the C++ EM solver loop.

### 4. Implementation Thoughts

**Where to compute $\tilde{L}_{unspliced}$?**
You can actually compute this vector entirely in Python before the EM starts.
When building the `AbundanceEstimator`, you already have the `transcript_exonic_lengths` and the fragment length distribution ($\mu_F$).
You just need to pass the list of individual exon lengths per transcript from the `Index` to compute the $\max(0, E_i - \mu_F + 1)$ sum.

How does this logic feel to you? If you agree with the geometric derivation, we can draft the exact numpy operations to fold this into the `estimator.py` initialization pipeline.. We should be able to estimate the number of estimate mature RNA fragments based on the fragment length distribution and the transcript structures. Geometric Splicing Prediction: The DerivationYour idea to predict $P(\text{unspliced} \mid \text{transcript}, \text{frag\_len})$ is completely viable. In RNA-seq, sliding a fragment across a transcript is formalized as computing its Effective Length ($\tilde{L}$).Let $\mu_F$ be the mean fragment length. For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:1. Total Valid Starting Positions (Total Effective Length):$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$2. Unspliced Starting Positions:For a read to be strictly "unspliced," it must fall entirely within the boundaries of a single exon.$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$3. The Splicing Ratio:Assuming uniform transcript coverage, the probability that a fragment drawn from mature mRNA $t$ appears unspliced is an immutable geometric property:$$P(\text{unspliced} \mid t) = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{total}}$$4. Executing the Deconvolution:During initialization, count the number of unambiguously spliced reads assigned to transcript $t$ ($C_{spliced}$). Because they cross a junction, we know with 100% certainty they are mature mRNA.We can calculate the "Phantom Unspliced Burden" that this mRNA must have left in the locus:$$Burden_{mRNA} = \sum_{t} \left( C_{spliced, t} \times \frac{P(\text{unspliced} \mid t)}{1 - P(\text{unspliced} \mid t)} \right)$$Subtract $Burden_{mRNA}$ from the total sense-strand unspliced pool before you ever calculate the gDNA or nRNA priors!  -- Can you do the following 1) Review and Explain the current gDNA initialization process. 2) Explain if and how gDNA can be overestimated when mRNA levels are very high. 3) Validate the accuracy of this theory, 4) After the BAM scanning pass, we could compute the P(unspliced | spliced, fragment length distribution) for each transcript. Ideally, we would initialize genomic DNA more conservatively based on 




# Geometric Splicing Expectation: Implementation Plan

## Problem

The gDNA initialization pipeline estimates per-locus gDNA counts from the
unspliced sense/antisense fragment pool using strand correction:

```
G = 2(A·SS - S·(1-SS)) / (2SS-1)
```

While this formula is **algebraically exact** (it perfectly cancels all RNA
contributions regardless of mRNA level), the unspliced pool is contaminated
by mature mRNA fragments that fall entirely within a single exon. These
"phantom unspliced" mRNA reads add **variance** to both S and A, which the
`1/(2SS-1)` amplification factor magnifies — especially at lower strand
specificity (SS=0.80 → amplification ≈ 2.8×).

This increased variance propagates through the 3-tier EB shrinkage. Noisy
per-locus estimates can miscalibrate the MoM κ estimation, and zero-gDNA
loci get pulled toward a global mean that itself has inflated variance.
The result: gDNA "hallucination" at high-expression loci.

## Solution: Geometric Splicing Expectation

For any multi-exon transcript, the ratio of unspliced to spliced reads is a
**deterministic geometric property** of the exon structure and fragment
length distribution. By observing spliced reads (which are 100% confirmed
mRNA — a spliced fragment cannot be gDNA), we can predict the expected
unspliced mRNA burden and subtract it from the pool before gDNA estimation.

### The Math (eCDF-Based)

Using a single mean fragment length μ_F is inaccurate for skewed fragment
length distributions. Instead, we integrate over the **empirical CDF** of
the observed fragment length distribution.

The key insight: for an exon of length E, the expected number of valid
fragment start positions that produce unspliced reads is:

```
L_unspliced_exon(E) = Σ_{f=1}^{E} P(F=f) × (E - f + 1)
                    = (E + 1) × (CDF(E) - P(0)) - CMOM(E)
```

where `CDF(k) = Σ_{l=0}^{k} P(l)` and `CMOM(k) = Σ_{l=0}^{k} l·P(l)`.

**This is exactly the formula used by `FragmentLengthModel.compute_all_transcript_eff_lens()`.**
Computing `eff_len(exon_length)` for each exon gives the per-exon unspliced
effective length. We reuse the existing infrastructure directly.

For transcript *t* with exon lengths E₁, E₂, …, Eₖ:

```
L_unspliced(t) = Σᵢ eff_len(Eᵢ)          (sum of per-exon effective lengths)
L_total(t)     = eff_len(Σᵢ Eᵢ)           (effective length of total exonic span)
L_spliced(t)   = L_total(t) - L_unspliced(t)
R_t            = L_unspliced(t) / L_spliced(t)
```

Under uniform coverage, if we observe C_spliced spliced reads from
transcript *t*, the expected unspliced mRNA burden is:

```
B_t = C_spliced(t) × R_t
```

### Why eCDF, Not Mean Fragment Length

The mean-based formula (`L_unspliced = Σ max(0, Eᵢ - μ_F + 1)`) is a
point estimate that ignores the shape of the distribution. The eCDF-based
formula integrates over the full fragment length distribution, properly
weighting each fragment length by its observed probability. This is
mathematically equivalent to computing per-fragment R_t values and taking
the probability-weighted mean — but computed once per exon length rather
than once per fragment.

For a 500 bp exon with a skewed distribution (mode=200, mean=250, tail
to 600), the eCDF naturally accounts for the fact that most fragments
(200 bp) fit easily (301 positions) while rarer long fragments (500 bp)
contribute only 1 position.

### Key Properties

- **R_t is an immutable geometric property** — it depends only on exon
  structure and the fragment length distribution, not on expression level.
- **Longer fragments → more junction-spanning → lower R_t.** A library
  with 400 bp fragments produces proportionally more spliced reads than
  one with 150 bp fragments. R_t automatically calibrates to the library.
- **Single-exon transcripts**: L_spliced = 0, so R_t = ∞. But since
  single-exon transcripts produce zero spliced reads (C_spliced = 0), the
  burden B_t = 0. No special-casing needed beyond guarding division by 0.
- **Very short exons** (all Eᵢ < min fragment length): L_unspliced ≈ 0,
  R_t ≈ 0. Nearly every fragment must span a junction. Correct.

## Counting Spliced Fragments

A spliced fragment is definitively mRNA regardless of transcript
ambiguity. A fragment that aligns to a unique genomic location but is
compatible with multiple overlapping transcripts is NOT a multimapper — it
is simply transcript-ambiguous. The splice junction proves RNA origin
regardless of which specific transcript it came from.

The rare exception is a fragment that multimaps (aligns to multiple genomic
locations) where some alignments are spliced and others are unspliced (e.g.,
processed pseudogenes). These are uncommon in practice and already handled
by the existing multimapper gating logic.

### EM Unit Semantics

**Each EM unit is one fragment** — not one fragment-transcript hit. The
`ScoredFragments` CSR has `is_spliced: bool[n_units]` with one entry per
fragment. Each unit then has a CSR row of candidate transcript indices in
`t_indices[offsets[u]:offsets[u+1]]`. So
`em_data.is_spliced[locus.unit_indices].sum()` counts spliced *fragments*,
not hits. No double counting.

### Available Spliced Counts

After the BAM scan pass, spliced fragment counts come from two sources:

1. **Deterministic unambiguous** (`FRAG_UNAMBIG + SPLICE_SPLICED_ANNOT`):
   bypass the EM entirely. Counted in `unambig_counts[t_idx, 4:6]`
   (SPLICED_ANNOT sense + antisense). These have a known transcript
   assignment.

2. **EM-entering spliced units**: all other spliced fragments (ambiguous,
   opposite-strand, multimapped, unannotated splice). Flagged by
   `em_data.is_spliced[unit_idx] == True`. These have multiple candidate
   transcripts but unknown assignment.

## Implementation

### Step 1: Compute R_t for Every Transcript (eCDF-Based)

**When**: After scan pass (fragment length model is trained).

**Where**: New function in `estimator.py`.

**Inputs**:
- `index._t_exon_intervals[t_idx]` → per-exon `[start, end)` intervals
- `frag_len_models.global_model` → `FragmentLengthModel` (provides the
  `compute_all_transcript_eff_lens()` method with eCDF infrastructure)

**Computation**:

```python
def compute_unspliced_to_spliced_ratios(
    index: TranscriptIndex,
    frag_len_model: FragmentLengthModel,
) -> np.ndarray:
    """Compute R_t = L_unspliced / L_spliced for every transcript.

    Uses the empirical fragment length distribution (eCDF) to compute
    effective lengths, properly weighting each fragment length by its
    observed probability.

    Returns inf for single-exon transcripts (L_spliced = 0).
    """
    n_t = index.num_transcripts

    # 1. Collect all per-exon lengths across all multi-exon transcripts
    #    into a flat array for vectorized eff_len computation.
    exon_lens_list = []
    exon_groups = []  # (t_idx, start_pos, n_exons) for reconstruction

    for t_idx in range(n_t):
        exon_ivs = index.get_exon_intervals(t_idx)
        if exon_ivs is None or len(exon_ivs) <= 1:
            continue
        lens = (exon_ivs[:, 1] - exon_ivs[:, 0]).astype(np.int64)
        exon_groups.append((t_idx, len(exon_lens_list), len(lens)))
        exon_lens_list.extend(lens)

    if not exon_lens_list:
        return np.full(n_t, np.inf, dtype=np.float64)

    all_exon_lens = np.array(exon_lens_list, dtype=np.int64)

    # 2. Vectorized eCDF-based effective length for all exons at once
    per_exon_eff = frag_len_model.compute_all_transcript_eff_lens(all_exon_lens)

    # 3. Also compute total effective lengths for each transcript
    total_exonic = index.t_df["length"].values.astype(np.int64)
    total_eff = frag_len_model.compute_all_transcript_eff_lens(total_exonic)

    # 4. Sum per-exon eff_lens by transcript → L_unspliced
    ratios = np.full(n_t, np.inf, dtype=np.float64)
    for t_idx, start, n_exons in exon_groups:
        L_unspliced = per_exon_eff[start : start + n_exons].sum()
        L_total = total_eff[t_idx]
        L_spliced = L_total - L_unspliced

        if L_spliced > 0:
            ratios[t_idx] = L_unspliced / L_spliced
        # else: keep inf (no spliced reads expected)

    return ratios
```

**Performance**: The Python loop iterates over multi-exon transcripts
(~100-200K), but all heavy computation (`compute_all_transcript_eff_lens`)
is vectorized numpy. The per-exon eff_len computation is a single
vectorized call over the flattened exon array. Expected runtime: <1s.

### Step 2: Compute Per-Locus mRNA Unspliced Burden

**When**: Inside `compute_eb_gdna_priors()`, per locus.

**Where**: `locus.py`, within the per-locus loop.

The burden has two components that must be handled differently because
they have different levels of transcript assignment certainty:

#### Component A: Unambiguous spliced reads (known transcript)

These have a deterministic transcript assignment. Use per-transcript R_t
directly:

```python
t_arr = locus.transcript_indices
burden_unambig = 0.0
for t_idx in t_arr:
    c_spliced = estimator.unambig_counts[t_idx, 4:6].sum()
    r_t = ratios[t_idx]
    if np.isfinite(r_t) and c_spliced > 0:
        burden_unambig += c_spliced * r_t
```

#### Component B: Ambiguous spliced EM units (unknown transcript)

Each spliced EM unit is one fragment with multiple candidate transcripts.
We don't know which transcript it came from (that's what the EM solves).
For each such unit, compute the expected burden as the mean R_t across
its candidate transcripts:

```python
n_transcripts = estimator.num_transcripts
burden_ambig = 0.0

spliced_mask = em_data.is_spliced[locus.unit_indices]
spliced_units = locus.unit_indices[spliced_mask]

for u in spliced_units:
    start = em_data.offsets[u]
    end = em_data.offsets[u + 1]
    candidates = em_data.t_indices[start:end]

    # Filter to mRNA transcripts (exclude nRNA shadows)
    tx_cands = candidates[candidates < n_transcripts]
    if len(tx_cands) == 0:
        continue

    cand_ratios = ratios[tx_cands]
    finite = np.isfinite(cand_ratios)
    if finite.any():
        burden_ambig += cand_ratios[finite].mean()
```

**Why mean R_t across candidates?** We lack transcript assignment
probabilities at this stage (pre-EM). Equal weighting is the
maximum-entropy prior. If one candidate has R_t=0.3 and another has
R_t=2.0, the mean (1.15) reflects our uncertainty. Once the EM runs,
the actual assignment will sort reads correctly — the burden subtraction
only needs to be a reasonable initialization prior.

**No double counting:** Each EM unit is one fragment. The inner loop
iterates over the unit's candidate transcripts to compute R_t, but adds
only ONE burden contribution per fragment (the mean R_t).

#### Total burden

```python
burden = burden_unambig + burden_ambig
```

### Step 3: Subtract Burden from Unspliced Pool

The mRNA unspliced burden is strand-biased (mostly sense). Distribute by
strand specificity:

```python
burden_sense = burden * strand_spec
burden_anti = burden * (1.0 - strand_spec)

adjusted_sense = max(0.0, locus_unspliced_sense - burden_sense)
adjusted_anti = max(0.0, locus_unspliced_anti - burden_anti)
```

Use `adjusted_sense` and `adjusted_anti` in place of raw counts when
computing `compute_gdna_density_hybrid()` for this locus.

### Step 4: Integrate into Pipeline

The insertion point is `pipeline.py` between the scan pass and
`compute_eb_gdna_priors()` (line ~542):

```python
# After scan, compute geometric splicing ratios (eCDF-based)
unspliced_to_spliced_ratios = compute_unspliced_to_spliced_ratios(
    index, frag_len_models.global_model,
)

# Pass to compute_eb_gdna_priors (new parameter)
gdna_inits = compute_eb_gdna_priors(
    loci, em_data, estimator, index, strand_models,
    unspliced_to_spliced_ratios=unspliced_to_spliced_ratios,
    ...
)
```

### Files Modified

| File | Change |
|------|--------|
| `src/rigel/estimator.py` | New `compute_unspliced_to_spliced_ratios()` using eCDF |
| `src/rigel/locus.py` | Modify `compute_eb_gdna_priors()`: accept ratios, compute per-locus burden (unambig + ambig), subtract from unspliced pool |
| `src/rigel/pipeline.py` | Compute ratios after scan, pass to EB prior computation |

No C++ changes. No new configuration parameters.

## Edge Cases

| Scenario | Handling |
|----------|----------|
| Single-exon transcript | R_t = ∞, C_spliced = 0 → burden = 0 |
| All exons < min frag len | L_unspliced ≈ 0, R_t ≈ 0 → no burden |
| No spliced reads in locus | burden = 0 (nothing to subtract) |
| Burden > unspliced pool | `max(0, ...)` clamp prevents negatives |
| Locus with only single-exon txs | All R_t = ∞, all finite checks fail → burden = 0 |
| Ambig unit with no mRNA candidates | `tx_cands` empty → skip (contributes 0) |
| Non-uniform coverage (5'/3' bias) | Burden is an expectation; `max(0, ...)` prevents over-subtraction. Coverage bias may cause actual unspliced mRNA to differ from predicted, but the clamp ensures safety. |

## Why This Works

1. **Eliminates the dominant noise source.** For a highly expressed gene
   with 10,000 reads and R_t = 0.6, the unspliced pool contains ~3,750
   mRNA phantom reads. Subtracting them reduces Poisson variance fed to
   the strand correction formula.

2. **Automatically scales with fragment length.** Longer fragments span
   more junctions (lower R_t), so the predicted burden adapts to the
   library. The eCDF integration handles skewed distributions correctly.

3. **Zero new parameters.** R_t is a deterministic geometric property
   computed from the empirical fragment length distribution and transcript
   annotation. No tuning, no κ, no pseudo-counts.

4. **Conservative by design.** Uses observed spliced counts (a lower bound
   on true mRNA) and clamps at zero. Cannot over-subtract.

5. **Pure Python, pre-EM.** Adds zero overhead to the C++ EM solver.
   Runs once during initialization.

6. **Reuses existing infrastructure.** The eCDF-based effective length
   computation (`FragmentLengthModel.compute_all_transcript_eff_lens`)
   already implements the exact CDF/CMOM table lookups needed. No new
   numerical machinery.

## Verification

1. **Unit test**: Construct synthetic transcripts with known exon
   structures. Verify R_t matches hand-computed values. Test with both
   uniform and skewed fragment length distributions to confirm the eCDF
   handles non-trivial cases.

2. **Golden test regression**: Run existing golden tests. The burden
   subtraction should produce slightly tighter gDNA estimates but should
   not break any test that passes today (gDNA estimates may decrease
   slightly at high-mRNA loci).

3. **Diagnostic sweep**: Run the benchmarking suite with the gDNA
   titration configs (`scripts/benchmarking/diag_gdna_titration.yaml`).
   Compare gDNA estimation accuracy before/after burden subtraction,
   especially at high mRNA expression levels and low SS.

4. **Targeted scenario**: Create a test locus with 10,000 mRNA reads,
   0 true gDNA, and verify that the burden subtraction drives gDNA init
   toward 0 (whereas without it, the EB shrinkage pulls it upward).
