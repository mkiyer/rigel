# gDNA Regional Exposure Plan v4.1

**Status**: design plan, not implemented.
**Date**: 2026-05-19
**Supersedes**: `gdna_exposure_plan_v4.md` (and all earlier plans).

v4.1 keeps v4's conceptual reset and corrects the following defects
identified in review:

- v4 conflated **overlap** geometry (gDNA locus) with **containment**
  geometry (transcript / nRNA span) under one weighted-length call.
  v4.1 separates validity from integration: one integration engine,
  two validity rules.
- v4 claimed the capture numerator `A(x_u)` always cancels in the
  E-step. It cancels only when competing candidates share the
  observed capture coordinate. v4.1 states the honest argument and
  the engineering reason we still drop the numerator.
- v4 demanded EM-output equality under a constant non-identity
  field `A = A_0`. That is only true with the prior held fixed. v4.1
  fixes the test.
- v4 had no performance plan for whole-transcriptome weighted mRNA
  lengths. v4.1 adds one.
- v4 left `AbundanceEstimator`'s EM-vs-output length split implicit.
  v4.1 makes it explicit and required.

Everything else from v4 — single global `ρ_ref`, denominator-only EM,
mandatory `η_g · L_g^A / L_g` prior rescaling, raw public TPM, no
native ABI change, no new CLI flags — stands unchanged.

---

## 0. Principle

One physical assay. One capture field. One integrator. Two validity
rules. Denominator-only in EM.

Hybrid capture is molecule-agnostic. A fragment from genomic position
`x` is retained with probability `A(x)`. The EM-visible consequence is
a single rule for every competing component:

> Replace each component's effective length `L_k` with its
> exposure-weighted form
> `L_k^A = Σ_ℓ pmf_k(ℓ) · ∫_{M_k(ℓ)} A(s) ds`,
> where `M_k(ℓ)` is the **genomic midpoint window set** for fragments
> of length `ℓ` that are valid for component `k`.

`M_k(ℓ)` is determined by `k`'s footprint and `k`'s validity rule:

- **Containment** (transcript-coordinate fit): used by mRNA
  (single- or multi-exon) and synthetic nRNA spans. A fragment is
  valid iff its `ℓ` bases lie entirely inside the footprint in
  transcript order. Multi-exon transcripts span genomic introns by
  construction; the fragment occupies `ℓ` transcript bases.
- **Overlap** (genomic-interval intersection): used by the gDNA
  multi-locus. A fragment is valid iff its genomic interval
  intersects the footprint's union. Fragments straddling locus
  boundaries are still gDNA evidence.

Both validity rules produce a set of genomic midpoint windows
`M_k(ℓ)`. A single integration engine integrates `A` over those
windows. The biology lives in the validity rule; the integrator is
shared.

---

## 1. Generative model and the role of the numerator

### 1.1 Capture-aware likelihood

For component `k`:

$$
f_k(x, \ell) \;=\; \frac{p_k(x, \ell)\,A(x)}{L_k^A},
\qquad
L_k^A \;=\; \sum_\ell h_k(\ell) \int_{M_k(\ell)} A(s)\,ds.
$$

### 1.2 When the numerator cancels, and when it doesn't

The E-step posterior at fragment `u` with candidate components `J(u)`:

$$
P(k\mid u) \;=\;
\frac{\theta_k\,p_k(u)\,A(x_u^{(k)})\,/\,L_k^A}
     {\sum_{j \in J(u)} \theta_j\,p_j(u)\,A(x_u^{(j)})\,/\,L_j^A}.
$$

The numerator cancels **exactly** only when every candidate `k ∈ J(u)`
assigns the same observed capture coordinate `x_u^{(k)} = x_u`. This
holds for the common case of multiple candidate molecule types
(mRNA / nRNA / gDNA) competing for one observed alignment, because
the observed fragment coordinates are fixed by the BAM record.

It does **not** hold across multimapper alternatives whose alignments
place the fragment at distinct genomic positions, and it does not hold
across competing spliced-vs-unspliced interpretations whose effective
midpoints differ.

Even where it does not cancel, v4.1 still drops the numerator from
production EM for two engineering reasons:

1. The regional exposure field is not a calibrated per-candidate
   capture-likelihood model. It is a smoothed sample-level field built
   from gDNA-typed conservative counts. Applying it as candidate-level
   evidence in EM under that smoothing was empirically harmful in
   v3-production (asymmetric leakage).
2. The denominator term `L_k^A` is identifiable and stable: it is a
   deterministic function of the component's footprint and the field,
   so it injects no per-fragment noise into EM.

The honest production rule is therefore: **use the denominator
correction; do not use the numerator.** Where the numerator would have
cancelled anyway it costs nothing; where it would not have cancelled,
we are conservatively declining to use a field we do not trust at the
per-candidate level.

### 1.3 What `θ_k` means

EM converges to `θ_k = (1/N) Σ_u P(k|u)`, the observed fragment
fraction routed to `k`. Up to the assumption that the numerator
cancels for the dominant ambiguity classes, `θ_k ∝ n_k · Ā_k` with
`Ā_k = L_k^A / L_k`, and `θ_k / L_k^A ∝ n_k / L_k`. The standard EM
rate-per-length output is capture-corrected.

---

## 2. The exposure field A(x)

Unchanged from v4. One global field, one global `ρ_ref`, no per-class
normalization.

### 2.1 Construction

1. **Per-region density.**
   `ρ̂_r = (Y_r + κ · ρ_global) / (E_r + κ)`,
   with `Y_r` the strand-corrected conservative gDNA count, `E_r` the
   class-appropriate FL-aware opportunity, `ρ_global = Σ Y_r / Σ E_r`,
   `κ` the existing EB shrinkage.
2. **Global reference.**
   `ρ_ref = weighted_quantile(ρ̂_r, weights = E_r, q = 0.95)` over all
   regions.
3. **Field.** `A_r = clip(ρ̂_r / ρ_ref, exp(LOG_A_FLOOR), 1.0)`.
4. **Auto-uniform.** If the global signal is below threshold,
   `A_r ≡ 1` and the pipeline takes the bit-exact uniform fast path.

### 2.2 Why one normalization

All three class evidence channels measure `ρ_g · A_local` at first
order; dividing by a single `ρ_ref ≈ ρ_g · A_typical` cancels `ρ_g`
and leaves the relative exposure on a common scale. The class
partition exists only to compute each region's opportunity `E_r`
correctly; it does not enter the normalization. See v4 §2.2 for the
table; unchanged here.

### 2.3 Sentinel test

A region whose `ρ̂` equals the global `ρ_ref` gets `A_r = 1`
regardless of class. A region whose `ρ̂` is 100× below `ρ_ref` gets
`A_r = 0.01` regardless of class.

---

## 3. The integration engine and the two validity rules

This is where v4.1 substantively differs from v4.

### 3.1 The integration engine

The shared primitive does one thing: given a set of per-reference
genomic midpoint windows for fragments of length `ℓ`, integrate `A`
over those windows.

```python
def integrate_exposure_over_midpoints(
    midpoint_windows: list[tuple[int, int, int]],  # (ref, g_lo, g_hi)
    exposure: ExposureField,
) -> float:
    """Sum of ∫A(s) ds over the given midpoint windows."""
    return sum(
        exposure.weighted_length_on_ref(ref, g_lo, g_hi)
        for (ref, g_lo, g_hi) in midpoint_windows
    )
```

`midpoint_windows` is required to be **disjoint** (the caller
de-duplicates / merges before passing in). The engine does no merging.

Under uniform `A ≡ 1`, this reduces to summing
`Σ (g_hi - g_lo)`, which is exactly the unweighted count of valid
midpoints (and therefore of valid starts).

### 3.2 The two validity rules

#### 3.2.1 Containment — transcript / nRNA span

Footprint: ordered exon blocks `[(ref, g_start_i, g_end_i), ...]` in
5' → 3' transcript order with strand `±1`. Single-exon and nRNA spans
are the N = 1 special case.

Let `L_t = Σ_i (g_end_i - g_start_i)`, `s_i = Σ_{j < i} block_length_j`.

For each `ℓ` with `pmf(ℓ) > 0`:
  `n_starts = L_t - ℓ + 1`; skip if `≤ 0`.
  `h = ℓ // 2`.
  Transcript-coordinate midpoint range: `[h, h + n_starts)`.
  For each block `i`:
    `t_lo = max(h, s_i)`
    `t_hi = min(h + n_starts, s_i + (g_end_i − g_start_i))`
    if `t_hi ≤ t_lo`: continue
    Project to genomic (strand-aware, see §3.4):
      strand +1: `g_lo = g_start_i + (t_lo − s_i)`,
                 `g_hi = g_start_i + (t_hi − s_i)`
      strand -1: `g_lo = g_end_i   − (t_hi − s_i)`,
                 `g_hi = g_end_i   − (t_lo − s_i)`
    emit `(ref_i, g_lo, g_hi)` for this `ℓ`'s midpoint window set.

Windows emitted within a single `ℓ` are disjoint by construction
(distinct blocks). Integrate, multiply by `pmf(ℓ)`, accumulate.

#### 3.2.2 Overlap — gDNA multi-locus

Footprint: a set of genomic intervals
`[(ref, g_start_i, g_end_i), ...]` (one or more per ref, possibly
overlapping or nearby). Validity: a fragment is gDNA-valid iff its
genomic `ℓ`-base interval intersects the footprint union.

For each `ℓ` with `pmf(ℓ) > 0`:
  For each interval, the valid genomic **start** range is
  `[g_start_i − ℓ + 1, g_end_i)` (overlap semantics).
  Convert to **midpoint** windows by shifting by `+h = +ℓ // 2`:
  `[g_start_i − ℓ + 1 + h, g_end_i + h)`
  (use the same midpoint convention as the existing
  `weighted_gdna_eff_len_for_loci`; v4.1 inherits its convention
  exactly to preserve bit-exactness in uniform mode).
  Group by `ref`, **sort and merge** overlapping/adjacent windows.
  Emit the merged disjoint window set; integrate; multiply by
  `pmf(ℓ)`; accumulate.

The merge step is the critical engineering fact: when a multi-locus
contains nearby blocks (common after locus expansion by FL), naive
per-block summation double-counts starts that fall in the overlap of
two expanded windows. The existing `weighted_gdna_eff_len_for_loci`
merges. v4.1 must preserve this. The contract is "uniform-mode
bit-exact against `gdna_eff_len_for_loci`."

Implementation note: v4.1 ports
`weighted_gdna_eff_len_for_loci` so that it constructs the merged
midpoint windows per `ℓ` and then calls
`integrate_exposure_over_midpoints`. The function's external behavior
does not change; only its internal structure is refactored to share
the engine with the containment path.

### 3.3 Public API

Two thin callers, one engine:

```python
def weighted_eff_len_contained(
    footprint: list[tuple[int, int, int]],   # ordered exon blocks
    fl: FragmentLengthModel,
    exposure: ExposureField,
    *,
    strand: int = 1,
    min_value: float = 1.0,
) -> float: ...

def weighted_eff_len_overlap(
    footprint: list[tuple[int, int, int]],   # genomic intervals, any order
    fl: FragmentLengthModel,
    exposure: ExposureField,
    *,
    min_value: float = 1.0,
) -> float: ...
```

Each builds the midpoint windows under its validity rule, ensures
disjointness (containment: by construction; overlap: by merge), calls
`integrate_exposure_over_midpoints`, weights by `pmf`, and clamps.

### 3.4 Strand

Containment must respect strand: the genomic projection mapping
differs for `+1` and `-1`. The exposure integral itself is symmetric
in `(g_lo, g_hi)`, so the choice only affects which genomic range a
given transcript-midpoint range maps to. A negative-strand transcript
with its exon blocks listed in transcript order (highest genomic
position first) must yield the same `L_t^A` whether expressed as
`strand=-1` with blocks in transcript order, or as `strand=+1` with
blocks in genomic order. This is testable (§9.1.6).

Overlap is strand-agnostic.

### 3.5 Properties

- **Uniform fast path bit-exact.** Under `A ≡ 1`, both rules return
  the existing unweighted effective length to floating-point
  precision. Containment matches `compute_all_transcript_eff_lens`;
  overlap matches `gdna_eff_len_for_loci`.
- **Constant non-identity field.** Under `A ≡ A_0`,
  `L_k^A = A_0 · L_k` for both rules within fp tolerance.
- **gDNA window merging is mandatory** in the overlap rule.

---

## 4. Applying the engine to components

### 4.1 gDNA, per multi-locus

```python
L_g_A = weighted_eff_len_overlap(
    footprint=[(loc.ref, loc.g_start, loc.g_end) for loc in M.loci],
    fl=gdna_fl,
    exposure=exposure,
)
```

Multi-ref multi-loci flow naturally because the merge step operates
per ref.

### 4.2 Synthetic nRNA, per unique span

```python
L_n_A[s] = weighted_eff_len_contained(
    footprint=[(s.ref, s.g_start, s.g_end)],
    fl=rna_fl,
    exposure=exposure,
    strand=s.strand,
)
```

### 4.3 Annotated mRNA, per transcript

```python
L_m_A[t] = weighted_eff_len_contained(
    footprint=exon_blocks_for(t),       # 5'→3' transcript order
    fl=rna_fl,
    exposure=exposure,
    strand=t.strand,
)
```

Single-exon transcripts are N = 1.

---

## 5. Performance plan for whole-transcriptome mRNA weighting

A naive Python loop over `~250 k` transcripts × FL pmf support
`~10²–10³` × exons-per-transcript will be too slow. v4.1 commits to
the following before §10 step 4 lands:

### 5.1 Per-reference exposure prefix sums

The exposure field is a step function on per-reference region
intervals. Build a per-reference cumulative integral so
`weighted_length_on_ref(ref, g_lo, g_hi)` is `O(log R_ref)` (binary
search) or `O(1)` (if `g_lo, g_hi` snap to region boundaries) per
call. The existing `RegionalGdnaExposure.weighted_length_on_ref`
already exists; v4.1 confirms its complexity and, if needed,
materializes per-ref cumulative arrays for `O(log R)` lookup.

### 5.2 Vectorize per transcript

For each transcript, vectorize the `ℓ` loop over the FL pmf support
rather than iterating in Python: build the per-ℓ midpoint ranges as
numpy arrays, call `weighted_length_on_ref` once per `(ℓ, block)`,
weight by `pmf` in one dot product.

### 5.3 FL truncation

Truncate the FL pmf to `pmf ≥ 1e-6` (existing convention in
`FragmentLengthModel`). Reduces inner loop from full support to the
effective support `~10²`.

### 5.4 CSR exon layout

Use the existing CSR exon table (`index._t_exon_intervals` or
equivalent) so per-transcript exon access is contiguous numpy slices,
not Python lists.

### 5.5 Optional native port

Containment weighted-length over the whole transcriptome is the most
expensive piece; gDNA overlap is bounded by locus count
(`~10⁴–10⁵`) and is not a concern. If §5.1–§5.4 do not hit the
benchmark gate, port `weighted_eff_len_contained` to
`_exposure_impl` C++ (no new ABI; called from Python).

### 5.6 Benchmark gate (mandatory)

On the synthetic 24-condition reference index
(`~250 k` transcripts), the regional-exposure-enabled pipeline must
add no more than **10%** to total wall-clock vs `--regional-exposure
off`. If the Python implementation alone clears this gate, the native
port is not required.

---

## 6. gDNA prior rescaling

Unchanged from v4. Mandatory:

$$
\eta_g^{(v4)}(M) \;=\; \eta_g(M)\,\cdot\,\frac{L_g^A(M)}{L_g(M)}.
$$

Diagnostic column `gdna_prior_count_unweighted` is emitted for one
release alongside the EM-consumed `gdna_prior_count`.

---

## 7. Output semantics and estimator plumbing

### 7.1 Public outputs stay raw

| File / column | Value |
|---|---|
| `quant.feather effective_length` | raw `L_t` |
| `quant.feather tpm`, `tpm_total_rna` | computed with raw `L_t` |
| `nrna_quant.feather effective_length`, `tpm` | raw `L_n` |
| `gene_quant.feather effective_length` | raw (current aggregation) |
| `loci.feather gdna_eff_len` | raw `L_g` |

### 7.2 Diagnostic columns

| Column | Where | Meaning |
|---|---|---|
| `em_exposure_weight` | `quant.feather` | `L_t^A / L_t` |
| `em_effective_length` | `quant.feather` | `L_t^A` |
| `em_exposure_weight` | `nrna_quant.feather` | `L_n^A / L_n` |
| `em_effective_length` | `nrna_quant.feather` | `L_n^A` |
| `gdna_em_exposure_weight` | `loci.feather` | `L_g^A / L_g` |
| `gdna_em_effective_length` | `loci.feather` | `L_g^A` |
| `gdna_prior_count_unweighted` | `loci.feather` | pre-rescale `η_g` (one release) |

### 7.3 Estimator plumbing (required)

`AbundanceEstimator` currently carries one transcript effective-length
array used for both EM and TPM output. v4.1 splits this explicitly:

```python
class AbundanceEstimator:
    effective_lengths_em: np.ndarray        # used by native EM
    effective_lengths_output: np.ndarray    # used by TPM and reports
```

Likewise per-locus:

```python
locus.gdna_eff_len_em: np.ndarray          # used by native EM
locus.gdna_eff_len_output: np.ndarray      # used by reports
```

In uniform / `off` mode, the two arrays are the same object (or
bit-exact copies). In regional mode they diverge. TPM and reports
must read `_output`; native EM must read `_em`. This is enforced by
type signature (two named fields) rather than convention.

---

## 8. Pipeline order

Same as v4 §7 with one rename: the unified call is now the appropriate
of `weighted_eff_len_overlap` / `weighted_eff_len_contained`.

1. Read `calibration.regional_exposure`.
2. Build raw `effective_lengths_output` and per-locus
   `gdna_eff_len_output` via existing paths.
3. If `exposure.mode == "uniform"`, set `_em` arrays to copies of
   `_output` arrays and take the bit-exact fast path.
4. Otherwise compute:
   - `L_g^A` per locus via `weighted_eff_len_overlap`.
   - `L_n^A` per nRNA span via `weighted_eff_len_contained` (N=1).
   - `L_m^A` per annotated transcript via `weighted_eff_len_contained`
     (N≥1, strand-aware).
   - Assemble `effective_lengths_em` and `gdna_eff_len_em`.
   - Rescale `η_g` per §6.
5. Construct `AbundanceEstimator` with both arrays distinguished.
6. Score fragments. No numerator weighting.
7. Build multi-loci, assemble priors, run EM with `_em` arrays.
8. Emit raw public abundance and §7.2 diagnostics.

---

## 9. Validation

### 9.1 Unit tests (CI)

1. **Uniform fast-path bit-exactness.** Whole-pipeline run with
   `exposure.mode == "uniform"` is bit-exact against
   `--regional-exposure off`.
2. **Overlap uniform parity.** `weighted_eff_len_overlap` under
   `A ≡ 1` equals `gdna_eff_len_for_loci` exactly on a battery of
   multi-locus shapes including overlapping and adjacent blocks
   (catches a missing merge).
3. **Containment uniform parity.** `weighted_eff_len_contained` under
   `A ≡ 1` equals `compute_all_transcript_eff_lens` for a fixture of
   single-exon and multi-exon transcripts.
4. **Constant `A_0` length law.** Under `A ≡ A_0` (with
   `A_0 ∈ {0.1, 0.5, 0.9}`):
   `L_k^A == A_0 · L_k` for each of gDNA, nRNA, mRNA within fp
   tolerance (relative error `< 1e-12`).
   *Do not* assert EM posterior equality; the prior rescaling
   changes the absolute pseudocount strength, which can shift
   prior-sensitive fixed points.
5. **Prior density preservation.**
   `η_g^{(v4)} / L_g^A == η_g / L_g` within fp tolerance for every
   multi-locus on a held-out fixture.
6. **Negative-strand exon-order test.** A transcript with
   `strand = -1` and exon blocks listed in transcript order yields
   the same `L_m^A` as the same transcript expressed with
   `strand = +1` and blocks reordered to genomic order, under a
   non-trivial field (relative error `< 1e-10`).
7. **Class-agnostic field test.** §2.3 sentinel.
8. **Window-merge stress.** Two overlapping locus blocks at
   distance `< FL_max` produce a `weighted_eff_len_overlap` value
   strictly less than the sum of the per-block values under any
   non-degenerate field (catches a regression that removes the
   merge).

### 9.2 Synthetic two-region locus

Deterministic locus, 5 kb high-exposure exon (`A_r = 1.0`), 95 kb
low-exposure intron (`A_r = 0.01`), 100 true gDNA + 100 true nRNA + 50
true mRNA fragments. Assertions:

- v4.1 recovers gDNA / nRNA / mRNA pool counts within ±10% of truth.
- No pool over-routed by more than 15%.
- Uniform-override control reproduces the pre-v4 baseline.
- Runtime < 2 s.

### 9.3 VCaP production gates

| Gate | Threshold |
|---|---|
| G1: gDNA → RNA leak | ≤ uniform baseline + 1 pp |
| G2: gDNA → nRNA leak | ≤ uniform baseline + 0.5 pp |
| G3: gDNA → mRNA leak | ≤ uniform baseline + 0.5 pp |
| G4: RNA → gDNA leak | within ±0.5 pp of uniform baseline |
| G5: mRNA → gDNA leak | ≤ uniform baseline + 0.5 pp |
| G6: synthetic 24-condition mRNA error | ≤ uniform baseline + 0.5% |
| G7: synthetic 24-condition nRNA error | ≤ uniform baseline + 0.5% |
| G8: `--regional-exposure off` regression suite | bit-exact |
| G9 (perf): wall-clock vs `off` | ≤ +10% on synthetic 24-cond. index |

Failure investigation order:

1. Field construction (§2): inspect `ρ_ref` and per-region `A_r`.
2. Validity rule (§3.2): which rule was applied; window-merge audit
   for gDNA loci.
3. Per-component weight-ratio diagnostics (§7.2).

Do not reach for a per-fragment numerator. If denominator-only does
not work, the field or the validity rule is wrong.

---

## 10. Implementation checklist

1. **Single-`ρ_ref` field.** Replace per-class normalization in
   `_regional_exposure.py` with §2.1 global construction. Add §9.1.7.
2. **Integration engine + overlap rule.** Implement
   `integrate_exposure_over_midpoints` and
   `weighted_eff_len_overlap` in `_exposure.py`. Port
   `weighted_gdna_eff_len_for_loci` to call them (preserving its
   external behaviour bit-exact in uniform mode). Add §9.1.2 and
   §9.1.8.
3. **Containment rule.** Implement `weighted_eff_len_contained` for
   N ≥ 1 blocks. Add §9.1.3 and §9.1.6.
4. **Performance** (§5). Add the §9.3 G9 benchmark in CI as a
   non-blocking metric initially; promote to blocking once stable.
5. **nRNA span bookkeeping.** Add `index.t_to_nrna_span_id` and
   `index.nrna_span_df` at `TranscriptIndex.load()`.
6. **`L_n^A` and `L_m^A`.** Wrappers around the containment rule.
   Uniform parity tests.
7. **Estimator plumbing** (§7.3). Split `effective_lengths_em` /
   `_output` and `gdna_eff_len_em` / `_output`. TPM/report sites must
   read `_output`; EM must read `_em`.
8. **`η_g` rescaling** (§6) and §9.1.5 invariant.
9. **Delete numerator path.** Remove `_apply_unit_gdna_weights` from
   `quant_from_buffer`. Remove `safe_unit`, cross-ref skip, and
   multimapper exposure scaffolding. `RegionalWeightApplicationStats`
   returns zeros for one release with a deprecation log line; then is
   removed.
10. **Diagnostic columns** (§7.2).
11. **§9.2 synthetic two-region locus** CI test.
12. **VCaP gates** (§9.3) and synthetic 24-condition regression.
13. **Ship** behind `--regional-exposure auto` on gate pass; mark all
    prior plans superseded.

No native EM ABI change. No new CLI flags. No internal modes. A
possible native port of `weighted_eff_len_contained` is a
performance-only addition; it does not appear in EM or scoring
interfaces.

---

## 11. What we delete

| Removed | Reason |
|---|---|
| `_apply_unit_gdna_weights` production call | denominator-only |
| `_apply_candidate_nrna_weights` (proposal) | same |
| Any mRNA per-candidate numerator (proposal) | same |
| `safe_unit` / cross-ref / multimapper exposure scaffolding | served the numerator |
| Per-class `ρ_ref` | §2 global field |
| `transcript_exposure_weights` / `Abar_t` average | §3.2.1 exact containment integral |
| Any "internal mode" / hidden CLI flags | one mode |

The native scoring `genomic_midpoint` field stays as a fragment-level
QC diagnostic and is not consumed by EM.

---

## 12. Open questions

1. **Field evidence sources.** v4.1 uses gDNA conservative count only.
   Multi-source field estimation is a follow-up.
2. **Midpoint vs full-footprint integration.** v4.1 integrates `A` at
   the genomic midpoint. For typical FL ≪ region size this is
   essentially exact. If diagnostics expose artefacts on long
   fragments crossing low-A regions, tighten to a per-base integral
   in a follow-up.
3. **Cross-sample comparability.** `L_k^A` is sample-internal; that
   is why §7.1 keeps public TPM raw.
4. **Native containment port.** Triggered only by the §5.6 benchmark
   gate.

---

## 13. Bottom line

One capture field. One integration engine. Two validity rules
(containment for transcripts and nRNA spans; overlap for gDNA loci).
Denominator-only in EM. Mandatory prior rescaling. Public outputs
stay raw. EM and output effective-length arrays are explicitly split
in `AbundanceEstimator`. No native ABI change.

The biology stays simple. The engineering admits that gDNA overlap
and transcript containment are different validity rules and shares
only what they actually share: the integrator over genomic midpoint
windows.

If the VCaP gates pass, ship. If they fail, investigate the field or
the validity rule, not the numerator.
