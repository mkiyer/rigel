# gDNA Regional Exposure Plan v4

**Status**: design plan, not implemented.
**Date**: 2026-05-19
**Supersedes**: `gdna_exposure_revised_plan_v3a.md`,
`gdna_exposure_revised_plan_v3.md`, `gdna_exposure_revised_plan_v2.md`,
`gdna_exposure_plan_v3.md`.

This document is written from scratch. It does not preserve sections,
modes, or branches from prior plans. Anything not stated here is
out of scope for v4.

---

## 0. Principle

One physical assay. One capture field. One correction. One algorithm.

Hybrid capture is molecule-agnostic. A fragment from genomic position
`x` is retained with probability `A(x)`, regardless of whether the
molecule that emitted it was chromosomal DNA, nascent RNA, or mature
mRNA. The EM-visible consequence is a single rule that applies to every
competing component in exactly the same way:

> Replace each component's effective length `L_k` with its
> exposure-weighted form
> `L_k^A = Σ_ℓ pmf_k(ℓ) · Σ_{valid start} A(midpoint_genomic)`.

That is the whole correction. No per-fragment numerator. No per-class
normalization of the field. No component-type-specific branch in the
weighted-length code. There are exactly two distinctions in this design:

1. **Numerator vs denominator.** Capture appears in both terms of the
   capture-aware likelihood, but the numerator factor `A(x_u)` is the
   same constant at fragment `u` for every component and cancels
   exactly in the EM posterior (§1.2). Production EM only computes
   the denominator term.
2. **Footprint geometry.** A component's footprint is an ordered list
   of genomic exon blocks. Contiguous-span components (gDNA per locus,
   nRNA span, single-exon mRNA) have N = 1 block. Multi-exon mRNA has
   N > 1. The unified weighted-length function handles both with one
   inner loop over the block list (§3). N = 1 is not a special case;
   it is the algorithm with N = 1.

Everything else in this document follows from these two distinctions
plus a faithful EM prior rescaling.

---

## 1. Generative model and the cancellation

### 1.1 Capture-aware likelihood

For component `k` with intrinsic per-position emission density
`p_k(x, ℓ)` (uniform over the component's footprint, weighted by its
FL pmf), the observed density of a fragment at `(x, ℓ)` after capture
is:

$$
f_k(x, \ell) \;=\; \frac{p_k(x, \ell)\,A(x)}{L_k^A},
\qquad
L_k^A \;=\; \sum_\ell h_k(\ell) \int_{\text{valid for }k} A(s)\,ds.
$$

`A: genome → (0, 1]` is the assay accessibility field. It is a single
function on the genome and applies to every molecule type.

### 1.2 Why the numerator cancels

The EM E-step at fragment `u`:

$$
P(k\mid u) \;=\; \frac{\theta_k\,p_k(u)\,A(x_u)\,/\,L_k^A}
                       {\sum_j \theta_j\,p_j(u)\,A(x_u)\,/\,L_j^A}
            \;=\; \frac{\theta_k\,p_k(u)\,/\,L_k^A}
                       {\sum_j \theta_j\,p_j(u)\,/\,L_j^A}.
$$

`A(x_u)` is identical in every term and cancels. The EM trajectory,
fixed point, and per-fragment posteriors depend only on `p_k(u)` and
`L_k^A`. Production does not compute the numerator.

### 1.3 What `θ_k` means after the change

EM converges to `θ_k = (1/N) Σ_u P(k|u)`, the observed fragment
fraction routed to `k`. The true molecule abundance `n_k` relates to
`θ_k` via `θ_k ∝ n_k · Ā_k` with `Ā_k = L_k^A / L_k`. So
`θ_k / L_k^A ∝ n_k / L_k`. EM produces the capture-corrected
abundance estimate for free.

---

## 2. The exposure field A(x)

### 2.1 One field, one normalization

There are no classes in the field's calibration. The region partition
(INTERGENIC, INTRON, EXON-INTRON) exists in calibration purely to
compute the right opportunity denominator `E_r` for each region's
evidence channel; after that, every region's `ρ̂_r` is on the same
scale (`ρ_g · Ā_local`) and feeds a single global normalization.

Construction:

1. **Per-region density.** For every region `r` (any class):
   `ρ̂_r = (Y_r + κ · ρ_global) / (E_r + κ)`
   where `Y_r` is the strand-corrected conservative gDNA count, `E_r`
   is the FL-aware opportunity (class-appropriate; this is where the
   class partition does its only useful work), `κ` is the existing EB
   shrinkage, and `ρ_global = (Σ_r Y_r) / (Σ_r E_r)` is the **single
   global** gDNA density.
2. **Global reference.**
   `ρ_ref = weighted_quantile(ρ̂_r, weights = E_r, q = 0.95)`
   taken over **all regions**.
3. **Field.**
   `A_r = clip(ρ̂_r / ρ_ref, exp(LOG_A_FLOOR), 1.0)`
4. **Auto-uniform.** If the global signal (observed log-spread vs the
   deterministic null spread, computed over all regions) is below
   threshold, the entire field collapses to `A_r = 1` and the
   pipeline takes the bit-exact uniform fast path.

### 2.2 Why no per-class normalization

Each region's `ρ̂_r` is an estimate of `ρ_g · A_local`. The constant
`ρ_g` is the same for all regions (it is the bulk gDNA molecule
density). Dividing by a single global `ρ_ref ≈ ρ_g · A_typical`
removes that constant and leaves `A_local / A_typical`. Class-specific
normalization in prior plans was a workaround for the asymmetric
numerator, not a real calibration requirement.

The first-order argument:

| Class | `ρ̂_r` expectation under uniform `ρ_g` |
|---|---|
| INTERGENIC (contained) | `ρ_g · Ā_intergenic_position` |
| INTRON (contained) | `ρ_g · Ā_intronic_position` |
| EXON-INTRON (boundary-crossing) | `ρ_g · Ā_boundary_position` |

Same units, same scale, same `ρ_g`. One normalization.

### 2.3 Sentinel test

Construct synthetic regions whose `ρ̂` spans three orders of magnitude
across classes. Assert:

- The global `ρ_ref` is consistent with the 95th percentile of the
  combined distribution weighted by `E_r`.
- A region whose `ρ̂` equals `ρ_ref` receives `A_r = 1.0` regardless of
  its class.
- A region whose `ρ̂` is 100× below `ρ_ref` receives `A_r = 0.01`
  regardless of its class.

There is no class-specific reference to test; there should be no
class-specific behaviour.

---

## 3. The unified weighted effective length

### 3.1 The function

```python
def weighted_effective_length(
    footprint: list[tuple[int, int, int]],
    fl: FragmentLengthModel,
    exposure: ExposureField,
    *,
    strand: int = 1,            # +1 or -1; only affects genomic projection
    min_value: float = 1.0,
) -> float:
    """Exposure-weighted FL-marginal effective length for a single
    component with the given genomic footprint.

    Parameters
    ----------
    footprint : list of (ref_id, g_start, g_end)
        Genomic exon blocks in transcript-coordinate order
        (5' → 3'). Single-block footprints describe contiguous-span
        components (gDNA, nRNA spans, single-exon transcripts).
    """
```

### 3.2 Algorithm

Let `L_t = Σ_i (g_end_i − g_start_i)` (transcript length).
Let `s_i = Σ_{j < i} (g_end_j − g_start_j)` (transcript start of block
`i`).

For each `ℓ` with `fl.pmf[ℓ] > 0`:
  Let `n_starts = L_t − ℓ + 1`. If `n_starts ≤ 0`, skip.
  Let `h = ℓ // 2`.
  Transcript-coordinate midpoint range: `[h, h + n_starts)`.

  For each exon block `i`:
    Transcript range of midpoints inside block `i`:
      `t_lo = max(h, s_i)`
      `t_hi = min(h + n_starts, s_i + (g_end_i − g_start_i))`
    If `t_hi ≤ t_lo`: continue.
    Project to genomic midpoint range:
      if `strand == +1`:
        `g_lo = g_start_i + (t_lo − s_i)`
        `g_hi = g_start_i + (t_hi − s_i)`
      else: # strand == -1
        `g_hi = g_end_i − (t_lo − s_i)`
        `g_lo = g_end_i − (t_hi − s_i)`
    Accumulate:
      `total += fl.pmf[ℓ] * exposure.weighted_length_on_ref(ref_id_i, g_lo, g_hi)`

Return `max(total, min_value)`.

### 3.3 Properties

- **N = 1 reduces to the existing gDNA / nRNA path** (one block, one
  midpoint range per `ℓ`).
- **Uniform-mode fast path is bit-exact** to the existing unweighted
  effective length when `A(x) ≡ 1`, because then
  `weighted_length_on_ref(ref, g_lo, g_hi) = g_hi − g_lo`.
- **Constant non-identity field** (`A ≡ A_0`) returns
  `A_0 · L_unweighted` to floating-point tolerance.
- **Strand only changes the projection mapping**, not the integration.
  This matters for transcripts; gDNA and nRNA spans are
  strand-symmetric.

### 3.4 Validity rule

A fragment is valid for component `k` iff it lies entirely within `k`'s
footprint *in transcript coordinates*. For multi-exon mRNA, this means
spliced fragments are valid (they cross genomic introns by construction
of the transcript), and we never project an intron into the genomic
midpoint sum. This is exactly the existing mRNA effective-length
convention; v4 inherits it without modification.

---

## 4. Applying the unified function to every component

There is no component-specific code path. Each component supplies a
footprint and an FL model, and calls the same function.

### 4.1 gDNA, per multi-locus

For each `MultiLocus` `M`, sum over its constituent `Locus` intervals.
Each `Locus` provides one block `(ref, g_start, g_end)`. Multi-ref
multi-loci pass blocks on different refs.

```python
L_g_A = sum(
    weighted_effective_length([block], gdna_fl, exposure, strand=+1)
    for block in multilocus_blocks(M)
)
```

The strand parameter is irrelevant here because gDNA blocks are
strand-symmetric; `+1` is conventional.

### 4.2 Synthetic nRNA, per unique span

For each unique nRNA span `s` with genomic interval
`(ref, g_start, g_end)` and strand `±1`:

```python
L_n_A[s] = weighted_effective_length(
    [(ref, g_start, g_end)], rna_fl, exposure, strand=s.strand,
)
```

Broadcast per-span values to all transcripts via
`index.t_to_nrna_span_id` (derived at index load from existing columns
— no index format bump).

### 4.3 Annotated mRNA, per transcript

For each annotated transcript `t` with exon blocks
`[(ref, e1_start, e1_end), ..., (ref, eN_start, eN_end)]` and strand:

```python
L_m_A[t] = weighted_effective_length(
    exon_blocks_for(t), rna_fl, exposure, strand=t.strand,
)
```

Single-exon annotated transcripts use the same call with N = 1.

### 4.4 Hand-off to EM

```python
t_eff_len_for_em = t_eff_len.copy()
t_eff_len_for_em[synthetic] = L_n_A[t_to_nrna_span_id[synthetic]]
t_eff_len_for_em[annotated] = L_m_A[annotated]
```

The per-locus gDNA effective length passed to EM is `L_g^A`. No native
EM ABI change.

---

## 5. gDNA prior rescaling

The v6 prior `η_g(M) = ρ_g · L_g(M)` was calibrated against an
unweighted opportunity. EM places prior mass proportional to
`η_g / L_g`. If `L_g` is replaced by `L_g^A = Ā_g L_g`, the un-rescaled
prior mass becomes `η_g / L_g^A = ρ_g / Ā_g`, which is too strong by
`1/Ā_g`. To preserve the calibrated rate per unit observed
opportunity, scale:

$$
\eta_g^{(v4)}(M) \;=\; \eta_g(M)\,\cdot\,\frac{L_g^A(M)}{L_g(M)}.
$$

This is mandatory in v4. The pre-rescaled `η_g` is emitted as
`gdna_prior_count_unweighted` for one release and then removed.

There is no analogous rescaling for RNA components because the existing
RNA pipeline does not place an explicit Dirichlet pseudocount per
transcript that was calibrated against unweighted `L_t`.

---

## 6. Output semantics

Public abundance outputs stay on **raw, unweighted** lengths so that
TPM and `effective_length` remain comparable across samples with
different capture designs. The weighted lengths live inside the EM
denominator only.

| File / column | Value |
|---|---|
| `quant.feather effective_length` | raw `L_t` (unweighted) |
| `quant.feather tpm`, `tpm_total_rna` | computed with raw `L_t` |
| `nrna_quant.feather effective_length`, `tpm` | raw `L_n` |
| `gene_quant.feather effective_length` | raw (current aggregation) |
| `loci.feather gdna_eff_len` | raw `L_g` |

New diagnostic columns (clearly named so users can't compare them
across samples by accident):

| Column | Where | Meaning |
|---|---|---|
| `em_exposure_weight` | `quant.feather` | `L_t^A / L_t` per transcript |
| `em_effective_length` | `quant.feather` | `L_t^A` per transcript |
| `em_exposure_weight` | `nrna_quant.feather` | `L_n^A / L_n` per span |
| `em_effective_length` | `nrna_quant.feather` | `L_n^A` |
| `gdna_em_exposure_weight` | `loci.feather` | `L_g^A / L_g` per locus |
| `gdna_em_effective_length` | `loci.feather` | `L_g^A` per locus |
| `gdna_prior_count_unweighted` | `loci.feather` | pre-rescale `η_g` (one release) |

If exposure-weighted TPM is ever useful for downstream analysis it can
be added later as `tpm_exposure_weighted`; do not overload existing
TPM semantics.

---

## 7. Pipeline order

`CalibrationResult.regional_exposure` is built by
`rigel.calibration._orchestrator.calibrate(...)` before
`quant_from_buffer(...)` begins.

Inside `quant_from_buffer`:

1. Read `calibration.regional_exposure` (constructed per §2).
2. Build raw `t_eff_len` (existing path) and `L_g` per locus
   (existing path).
3. If `exposure.mode == "uniform"`, set
   `t_eff_len_for_em = t_eff_len`, `L_g^A = L_g`, and take the
   bit-exact fast path.
4. Otherwise:
   - Compute `L_g^A` per locus using `weighted_effective_length`
     (§4.1).
   - Compute `L_n^A` per nRNA span (§4.2).
   - Compute `L_m^A` per annotated transcript (§4.3).
   - Assemble `t_eff_len_for_em` (§4.4).
   - Rescale `η_g` per §5.
5. Construct `AbundanceEstimator` with `t_eff_len_for_em` and
   `L_g^A`.
6. Score fragments. (No numerator weighting anywhere.)
7. Build multi-loci, assemble priors (with rescaled `η_g`),
   partition, run EM.
8. Emit raw public abundance and the §6 diagnostic columns.

No `safe_unit` mask, no cross-ref skip, no multimapper exposure
bookkeeping. None of those concepts exist in v4 because there is no
numerator to apply.

---

## 8. What we delete from prior plans

| Removed | Reason |
|---|---|
| `_apply_unit_gdna_weights` from production | denominator-only; numerator cancels |
| `_apply_candidate_nrna_weights` (was v2 proposal) | same |
| Any mRNA per-candidate numerator (v2 Mode 3b) | same |
| `safe_unit` mask, cross-ref skip, multimapper skip | only existed in service of the removed numerator |
| Per-class `ρ_ref` (`rho_ref_per_class`) | §2.2 argument: one global `ρ_ref` |
| `transcript_exposure_weights` / `Abar_t` approximation | §3 unified algorithm does the integral exactly with N exon blocks |
| `unspliced_symmetric`, `gdna_only_legacy`, etc. internal modes | only one mode exists |
| `--regional-mode-internal` hidden CLI | one mode; no internal flag |

The native scoring `genomic_midpoint` field stays in
`ScoredFragments` as a per-fragment QC diagnostic (e.g. low-A flag
on fragment-level reports) but is not consumed by EM.

---

## 9. Validation

### 9.1 Unit tests (CI)

1. **Uniform fast path bit-exactness.** Under `exposure.mode ==
   "uniform"`, the entire run is bit-exact to `--regional-exposure off`.
2. **Cancellation invariant.** Construct a hand-built constant
   non-identity field `A(x) = A_0` covering a small locus. Assert:
   - `L_k^A == A_0 · L_k` for each of gDNA, nRNA, mRNA within fp
     tolerance.
   - `η_g^{(v4)} == η_g · A_0`.
   - EM posteriors and assignment counts equal the unweighted-baseline
     run within fp tolerance.
3. **Class-agnostic field test.** §2.3 sentinel: a region at the
   global `ρ_ref` gets `A_r = 1` regardless of class; a region 100×
   below gets `A_r = 0.01` regardless of class.
4. **Unified algorithm parity.** `weighted_effective_length` with
   N = 1 block matches the existing single-locus weighted path.
   `weighted_effective_length` with N exon blocks for a transcript
   under `A ≡ 1` matches the existing unweighted mRNA effective
   length bit-exactly.
5. **Strand handling.** A transcript with `strand = -1` and a
   reversed list of exon blocks gives the same `L_m^A` as the
   `strand = +1` formulation (the integral is symmetric).

### 9.2 Synthetic two-region locus (CI)

`tests/test_regional_exposure_locus.py`: one deterministic locus,
5 kb high-exposure exon (`A_r = 1.0`), 95 kb low-exposure intron
(`A_r = 0.01`). 100 true gDNA, 100 true nRNA, 50 true mRNA fragments.

Assertions:

- v4 recovers gDNA / nRNA / mRNA pool counts within ±10% of truth.
- No pool over-routed by more than 15%.
- A control run with the field overridden to uniform reproduces the
  pre-v4 baseline.
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

If any gate fails, the investigation order is:

1. Field construction (§2): inspect `ρ_ref`, the per-region `A_r`
   distribution, and the auto-uniform signal.
2. Unified algorithm (§3): inspect `L_k^A / L_k` weight ratios per
   component across the failing loci.

Do not reach for a per-fragment numerator to recover the gates. If
denominator-only does not work, the field is wrong.

---

## 10. Implementation checklist

Ordered for minimal-blast-radius landings. Each step ends in a
green-test state.

1. **Single-`ρ_ref` field.** Replace per-class normalization in
   `_regional_exposure.py` with the §2.1 global construction. Add the
   §2.3 sentinel test.
2. **Unified `weighted_effective_length`.** Implement §3 in
   `src/rigel/calibration/_exposure.py` so it handles `N ≥ 1` blocks.
   Port the existing `weighted_gdna_eff_len_for_loci` to call it (N = 1
   per locus block, sum). Add §9.1 tests 4 and 5.
3. **nRNA span bookkeeping.** Add `index.t_to_nrna_span_id` and
   `index.nrna_span_df` at `TranscriptIndex.load()`.
4. **`L_n^A` per span and `L_m^A` per transcript.** Wrappers around
   `weighted_effective_length`. Bit-exact uniform parity tests.
5. **Pipeline hand-off.** Build `t_eff_len_for_em` and pass to
   `AbundanceEstimator`. Public outputs remain on raw lengths.
6. **`η_g` rescaling.** §5. Add prior-ratio invariant test.
7. **Delete numerator path.** Remove `_apply_unit_gdna_weights` from
   `quant_from_buffer`. Remove `safe_unit`-style scaffolding. Remove
   `RegionalWeightApplicationStats` plumbing (or keep it returning
   zeroes for one release with a deprecation log line).
8. **Diagnostic columns.** Add the §6 EM-weight columns.
9. **§9.2 synthetic two-region locus** CI test.
10. **VCaP gates** (§9.3) and synthetic 24-condition regression.
11. **Decide.** If G1–G8 pass: ship behind `--regional-exposure auto`.
    If any fail: §9.3 escalation path.
12. **Document.** Mark all prior plans superseded; link this document
    from the calibration README.

No native ABI change. No new CLI flags. No internal modes.

---

## 11. Open questions

1. **Field evidence sources.** v4 estimates `A(x)` from gDNA
   conservative count alone. RNA evidence could in principle also
   contribute (each RNA molecule was captured by the same probes),
   but the algebra requires factoring out RNA molecule abundance per
   region, which couples back to what EM is trying to estimate. v4
   stays gDNA-only; multi-source field estimation is a follow-up.
2. **Footprint approximation.** v4 looks up `A` at the genomic
   midpoint of the fragment. The exact form integrates `A` over the
   full fragment footprint. For typical FL (≲ few hundred bp) and
   region sizes (≳ kb), the midpoint approximation is essentially
   exact. If diagnostics show artefacts for long fragments crossing
   exon boundaries, this can be tightened in a follow-up.
3. **Cross-sample comparability.** Because `A` is sample-specific
   (capture protocol varies), `L_k^A` is not comparable across
   samples. This is exactly why §6 keeps public TPM and effective
   length on raw values; the EM-weight diagnostics are
   sample-internal.

---

## 12. Bottom line

There is one capture field, one normalization, one weighted-length
algorithm, and one EM correction. The native scoring path does not
change. The pipeline shrinks because the numerator scaffolding is
gone. The per-component code paths collapse into one function called
with different footprints.

If the VCaP gates pass, ship. If they fail, fix the field (§2) — the
EM machinery is doing one thing and doing it correctly by construction.
