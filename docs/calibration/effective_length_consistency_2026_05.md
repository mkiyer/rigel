# Effective length: audit + redesign

**Status:** problem statement, audit, proposed fix
**Date:** 2026-05-06
**Trigger:** Synthetic-sim post-R1 results show ~40 unexpressed transcripts get nonzero counts in *every* condition, including `gdna_none`. Top offenders are short transcripts (e.g. `GENE0013.1`, 1.7 kb spliced, true count 0, observed count 2,749, observed TPM 11,463). The leakage exists with **zero gDNA contamination**, so this is not a calibration bug — it is an **EM mass-allocation bug rooted in inconsistent length corrections**.

---

## 0. Notation key

| Symbol | Meaning |
|---|---|
| `L` | spliced length of a transcript (sum of exon spans, in bp) |
| `ℓ_f` | observed fragment length of a single fragment (insert size, bp) |
| `h(·)` | fragment-length probability mass function (PMF) — `h(ℓ)` = probability that a randomly drawn fragment has length `ℓ`. We have two such PMFs: `h_R` (RNA) and `h_G` (gDNA), trained from unique mappers. |
| **PMF** | "probability mass function" — the discrete probability distribution of fragment lengths |
| `(L − ℓ + 1)` | "raw" effective length: the number of distinct start positions a fragment of length `ℓ` can occupy inside a region of length `L`, assuming the fragment fits (`ℓ ≤ L`) |
| `L̃` | **FL-marginal effective length** — the average of `(L − ℓ + 1)` over the fragment-length distribution: `L̃ = Σ_ℓ h(ℓ) · max(L − ℓ + 1, 0)`. Used by salmon, kallisto, RSEM. |
| `θ_t` | EM "abundance" parameter for transcript `t` (a probability across all components) |
| `γ_{f, t}` | "responsibility" — posterior probability that fragment `f` came from transcript `t` |
| `α_t` | Dirichlet-prior pseudocount on component `t` |

Whenever the doc says "log h(ℓ_f)" it just means "the log of the trained fragment-length PMF evaluated at the observed fragment length of this fragment" — a single number per fragment, the *same* across all candidate transcripts for that fragment.

---

## 1. The four places we currently compute "effective length"

| # | File / function | Pseudocode | Used for |
|---|---|---|---|
| **A** | [src/rigel/native/em_solver.cpp::apply_bias_correction_uniform](src/rigel/native/em_solver.cpp) (lines 333–356) | `log_lik[i] -= log( max(L_t − ℓ_f + 1, 1) )` per candidate `i` | Per-fragment likelihood inside EM (mRNA) |
| **B** | [src/rigel/frag_length_model.py::compute_all_transcript_eff_lens](src/rigel/frag_length_model.py) (lines 365–425) | `L̃_t = Σ_ℓ h_R(ℓ) · max(L_t − ℓ + 1, 0)` | mRNA TPM denominator only |
| **C** | [src/rigel/calibration/density_global.py::l_eff_contained](src/rigel/calibration/density_global.py) (lines 41–63) | `Σ_ℓ h_G(ℓ) · max(L_region − ℓ + 1, 0)` | gDNA Poisson-rate denominator (calibration) |
| **D** | [src/rigel/native/scoring.cpp](src/rigel/native/scoring.cpp) (lines 530–540) | `log_lik -= log( max(t_span + flank − ℓ_f + 1, 1) )` per gDNA hit | Per-fragment likelihood inside EM (gDNA) |

Plus a fifth, [src/rigel/bias.py](src/rigel/bias.py) (`BiasProfile.effective_length`), which is **dead code** — defined but never imported by the live EM path.

A and D use the **observed** `ℓ_f` plugged into `(L − ℓ + 1)`, then floor at 1. B and C use the **FL-PMF integral** (no floor below 0, only floored at 1 in B for log-safety).

---

## 2. Audit findings (the live code paths)

### 2.1 The partitioned EM does not normalize by transcript length

In [`em_solver.cpp::run_em_with_partition`](src/rigel/native/em_solver.cpp) (around line 2195):

```cpp
// 4. Log effective lengths (all 1.0 → log = 0.0)
std::vector<double> log_eff_len(nc, 0.0);
```

**The vector that the EM uses for per-component length normalization is hard-coded to zero.** The only length-correction the EM sees is whatever was already baked into `log_lik[i]` for each candidate — and that is formula A applied **once**, before EM, to each candidate fragment. There is no per-tx normalizer in the M-step or warm-start.

This means rigel's EM responsibility ratio between two transcripts `t1, t2` for one fragment is:

```text
γ_{f, t1} / γ_{f, t2}
  = (θ_{t1} / θ_{t2}) · exp(log_lik_A(t1) − log_lik_A(t2))
  = (θ_{t1} / θ_{t2}) · ( max(L_{t2} − ℓ_f + 1, 1) / max(L_{t1} − ℓ_f + 1, 1) )
```

vs. what salmon does:

```text
γ_{f, t1} / γ_{f, t2}  =  (θ_{t1} / θ_{t2}) · ( L̃_{t2} / L̃_{t1} )
```

Same shape, different denominator. They agree when `L_t ≫ ℓ_f`. They diverge sharply when `L_t` is close to `ℓ_f`.

### 2.2 Floor at 1 is the smoking gun

When `ℓ_f > L_t` the formula-A correction collapses to `log(1) = 0` — the **maximum** per-fragment likelihood. Physically this candidate is impossible (the fragment is longer than the transcript), but rigel scores it as the best possible match. Short transcripts pay this penalty repeatedly during a normal RNA-seq run.

### 2.3 The mRNA scoring path already includes the FL log-prob

Confirmed in [`scoring.cpp`](src/rigel/native/scoring.cpp) line 560: every mRNA candidate gets `log_fl = frag_len_log_lik(flen)` added to `log_lik`. The gDNA branch (line 535) does the same with its own model. So the per-fragment FL probability `log h(ℓ_f)` is already plumbed; it just is not what controls per-tx allocation.

### 2.4 The Dirichlet prior inherits the bias

[`compute_ovr_prior_and_warm_start`](src/rigel/native/em_solver.cpp) (line 779) builds the per-component prior proportional to `coverage_totals[i]`, which itself comes from the soft-resolution of unambig + ambiguous candidates using the formula-A log-liks. **Any inflation in formula A flows directly into the prior**, then gets reused by the EM.

### 2.5 TPM uses a different denominator than the EM

[`estimator.py::get_transcript_counts_df`](src/rigel/estimator.py) (line 418) computes TPM from `count / L̃_t` with `L̃_t` from formula B. So a count that was inflated 3× by a formula-A bug is then divided by an `L̃` of 1.5 kb (instead of 5 kb), producing TPM inflation of 3× × 3.3× ≈ **10×**.

### 2.6 `bias.py` is dead code

`BiasProfile` exists, has the right docstring, has the right math, and has a uniform-bias fast path that *would* compute formula B if it were called. **Nothing imports it for EM purposes.** The C++ scorer/EM does not see it. It is intended as a future hook for biological bias profiles (GC, 3' degradation), but the structural piece — a single per-tx effective-length number — never made it across the language boundary.

### 2.7 Candidate validity

`f_len[k]` is the per-candidate spliced fragment footprint along the transcript. In principle a fragment that overhangs a transcript end could land in `_resolve_impl` with `flen > t_length`. There is **no check** that rejects such candidates today; they reach formula A and get the floor-at-1 max likelihood. (Worth a defensive filter regardless of the larger fix.)

---

## 3. What kallisto / salmon / RSEM do (in pseudocode)

All three converge on the same contract.

```text
# Per-fragment likelihood
log P(f | t) = log w(f, t)            # bias / position / GC term, often 1
             + log h(ℓ_f)             # FL probability, SAME across t
             − log L̃_t                # per-tx normalizer, FL-marginal

# EM update (M-step)
θ_{t} ∝ Σ_f γ_{f, t}                  # responsibilities computed from above

# TPM
TPM_t = (count_t / L̃_t) · 1e6 / Σ_t' (count_{t'} / L̃_{t'})
```

The key invariant:

> **The denominator inside the per-fragment likelihood (`L̃_t`) is the same number used by TPM.**

Rigel violates this. The EM uses raw `(L − ℓ_f + 1)` floored at 1; TPM uses `L̃`.

---

## 4. The big picture: one definition of "effective length"

Conceptually, an effective length is the **expected number of fragment-start positions where a fragment can land inside a region**, where the expectation is taken over the FL distribution. That's it. It's a property of `(L_region, FL-PMF)`. We currently have three implementations of this same concept:

* mRNA TPM uses formula B with `h_R` over transcript length `L_t`
* gDNA calibration uses formula C with `h_G` over region span
* mRNA EM uses *neither* — it uses formula A with the observed `ℓ_f`, no FL averaging
* gDNA EM uses *neither* — it uses formula D with the observed `ℓ_f`, no FL averaging

The cleanest reorganization: **a single `effective_length(L, fl_model)` function**, called everywhere.

```python
# proposed unified API (lives in rigel.frag_length_model)

def effective_length(
    region_lengths: np.ndarray,    # int64 array of region spans (transcript lengths,
                                   # locus spans, region spans — same function)
    fl_model: FragmentLengthModel, # the FL distribution to use (RNA or gDNA)
    *,
    floor: float = 1.0,            # 1.0 for log-safe scoring/TPM, 0.0 for Poisson denom
) -> np.ndarray:
    """Σ_ℓ fl_model.pmf(ℓ) · max(L − ℓ + 1, 0), floored at `floor`."""
```

This is what `compute_all_transcript_eff_lens` already does. Promote it from "transcript-only" to the canonical primitive; have `density_global.l_eff_contained` and any new gDNA-locus code call it. Delete the dead code in `bias.py` (or rebuild `BiasProfile` to wrap this primitive when we're ready for non-uniform biases).

---

## 5. Proposed definitive fix

Three coordinated changes, landed in this order. Each compiles and passes tests independently, but the real win comes when all three are in.

### Step 1 (defensive guard, 1-line change)

In [`apply_bias_correction_uniform`](src/rigel/native/em_solver.cpp), reject candidates that are physically impossible:

```cpp
if (frag_len > prof_len) {
    log_liks[i] = -INFINITY;   // or large negative; will be pruned
    continue;
}
```

Cost: zero. Benefit: removes the worst-case "floor-at-1 = max likelihood" pathology. Should be a no-op if upstream resolution is correct, but defends against any mis-resolved candidate.

### Step 2 (mRNA: salmon contract)

Replace the per-fragment formula A with a **per-transcript** correction using the FL-marginal `L̃_t`.

Pseudocode (current):

```text
# Pre-EM: for each candidate i
log_lik[i] += -log( max(L_{t(i)} − ℓ_{f(i)} + 1, 1) )
# EM step: log_eff_len[*] = 0   (no per-tx normalizer)
```

Pseudocode (proposed):

```text
# Pre-EM: for each candidate i — NO length correction here, just keep
# log_lik[i] = log_strand + log h_R(ℓ_f) + ohpen + log_nm

# At EM solver entry:
for each transcript t in locus:
    L̃_t = effective_length(L_t, rna_fl_model, floor=1.0)
    log_eff_len[t] = log(L̃_t)
log_eff_len[gdna_idx] = log(L̃_gdna)   # see Step 3

# EM step (already supports it):
log_weights[c] = log(θ[c] + ε) − log_eff_len[c]
```

This is **bit-for-bit the salmon contract**. The C++ EM already accepts `log_eff_len`; today it is hard-coded to zero. We just have to populate it.

Required code edits:

1. `src/rigel/frag_length_model.py` — already has `compute_all_transcript_eff_lens`; expose it (already exposed).
2. `src/rigel/pipeline.py` — already computes `effective_lengths` via formula B (line 386). Pass these into the EM as `log(L̃_t)`.
3. `src/rigel/native/em_solver.cpp::run_em_with_partition` — replace `std::vector<double> log_eff_len(nc, 0.0);` with the per-locus slice of the `log L̃_t` array provided by Python.
4. `src/rigel/native/em_solver.cpp::apply_bias_correction_uniform` — **remove** entirely (or keep as no-op for the legacy non-partitioned path). Per-fragment correction is gone.
5. `src/rigel/native/scoring.cpp` — keep `log_fl = frag_len_log_lik(flen)` as is (it's the `log h_R(ℓ_f)` term, unchanged).
6. Goldens regenerate.

### Step 3 (gDNA: symmetric treatment)

The gDNA component currently has per-hit length correction `−log(t_span + flank − ℓ_f + 1)` floored at 1. Same pathology as Step 2 in miniature: when `t_span` is small, the gDNA likelihood is artificially saturated.

Replace with a **per-locus** FL-marginal effective length for gDNA, computed from `h_G`:

```text
# Pre-EM (Python side):
for each locus L:
    L̃_gdna(L) = effective_length(locus_span(L) + flank, gdna_fl_model, floor=1.0)
log_eff_len[gdna_idx_in_locus] = log(L̃_gdna(L))

# Per-fragment gDNA score (scoring.cpp):
hit_log_lik = log h_G(ℓ_f) + log_strand + LOG_HALF + log_nm
              # ← no per-fragment length correction
              # ← log_eff_len[gdna_idx] applied inside EM as for mRNA
```

This unifies gDNA and mRNA: both use the same `effective_length(...)` primitive, both apply normalization once per component inside the EM, both retain `log h(ℓ_f)` per fragment.

Required edits:

1. Compute `L̃_gdna` in the locus prior path ([`locus_prior.py`](src/rigel/calibration/locus_prior.py)) and persist into the per-locus EM payload.
2. Strip the per-hit `− log(e_h)` from the gDNA branch in [`scoring.cpp`](src/rigel/native/scoring.cpp).
3. Populate `log_eff_len[gdna_idx]` accordingly inside `run_em_with_partition`.

For tiny loci where `L̃_gdna < 1`, the floor at 1 leaves `log_eff_len = 0`, which is the correct neutral behavior — the gDNA hypothesis simply gets no length boost over a single-position model. This is exactly the "robust to tiny regions" property we want.

### Step 4 (dead-code cleanup, optional)

Delete or rewire [`bias.py::BiasProfile`](src/rigel/bias.py). Either:

* Delete (it has no callers), or
* Promote it to wrap `effective_length(...)` and use it as the canonical container when we eventually add GC / position bias profiles.

---

## 6. Expected impact on the synthetic sim

Predictions (to be validated post-fix):

* False-positive transcript count: 39–48 → single digits (the floor-at-1 pathology is gone).
* TPM inflation on short txs (e.g. `GENE0013.1` 1.7 kb): drop from 11,463 to a small residual (purely from upstream multimapper resolution noise, not from length-correction asymmetry).
* Spearman / Pearson(log): unchanged or slightly improved.
* gDNA recovery numbers (already good post-R1): unchanged. Step 3 makes the per-locus gDNA score more honest at small loci but doesn't change global density estimates.
* MARD: should drop substantially because the low-expression bin (TPM 10–30, mean RE 180 %) is dominated by short-tx leakage.

---

## 7. Open questions before implementation

1. **Pre-EM normalization vs. inside-EM normalization.** Currently the EM does both: formula A is applied *once* in `apply_bias_correction_uniform` (modifying `log_lik` in place), and the EM uses `log_eff_len` (currently zero) every iteration. Salmon's `log L̃_t` is applied *every iteration* via the M-step. The difference matters numerically for SQUAREM stability — recommend matching salmon: keep `log_lik` clean (no length correction), apply `log L̃` inside the iteration. The C++ already supports this.

2. **Multi-locus components.** When a transcript appears in multiple loci (rare but possible via multimappers), should `L̃_t` be the same number everywhere? Yes — it is a property of the transcript, not the locus. Just look it up by global `t_idx`.

3. **`L̃` for very short transcripts.** A 50-bp single-exon "transcript" would get `L̃ ≈ small` even with the floor at 1. This is the correct geometric behavior (almost no fragments fit), but it produces gigantic TPMs from any sliver of count. Salmon's behavior: report TPM = 0 below a length threshold and emit a warning. We could adopt a similar floor (e.g., `L̃ ≥ mean(h_R) / 4`) at TPM time; not strictly part of this fix, but worth deciding.

4. **The `flank` term in gDNA.** Currently `L_h = t_span + gdna_flank` where `gdna_flank ≈ mean(h_G)`. Under the new contract, `L̃_gdna = Σ_ℓ h_G(ℓ) · max(t_span + flank − ℓ + 1, 0)`. The flank exists to allow fragments that straddle one boundary; under the FL-marginal formula the flank may need re-derivation. This is the deepest open question — propose deriving it from first principles in a follow-up note before locking in Step 3.
