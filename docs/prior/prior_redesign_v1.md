
# Prior Reformulation v1 — Group-Level gDNA/RNA Priors

## Problem Statement

The current per-component Dirichlet prior couples two distinct concerns:
1. **Aggregate calibration evidence** — calibration credibly estimates "how much
   gDNA mass exists in this locus" (numerator) and the corresponding RNA mass.
2. **Within-RNA distribution** — calibration says *nothing* about how RNA mass
   should be partitioned across competing isoforms.

Any prior scheme that places positive pseudocount α > 0 on individual transcripts
acts as a permanent gravitational pull toward non-zero expression for that isoform.
For loci containing many alternative isoforms (often dozens), this:

- **Boosts unexpressed transcripts** — α/Σα floor mass survives every EM iteration.
- **Replicates the EM's own job** — partitioning RNA mass across isoforms is
  exactly what likelihood-driven EM is designed to do, with no prior input needed.
- **Degrades sparsity** — VBEM's appeal is its ability to drive components to
  zero; a positive RNA prior breaks that.

Meanwhile, the gDNA component is *singular* (one per locus) — pseudocount mass
applied to it has no isoform-spreading pathology. The asymmetry is real and the
prior architecture should reflect it.

---

## Design Goals

1. **gDNA gets a calibrated prior** with smooth user-tunable strength.
2. **RNA aggregate (Σθ_RNA) gets a calibrated counter-prior** to balance the gDNA
   pull and prevent siphoning when calibration says "this locus is mostly RNA."
3. **Individual transcripts get zero pseudocount** — EM partitions RNA mass
   freely based on likelihood alone.
4. **VBEM Jeffreys asymmetry eliminated** — single code path for MAP and VBEM.
5. **Smooth, gate-free, single user knob.**

---

## Option Survey

### Option A — Per-component symmetric Dirichlet (rejected — current paired-prior proposal)

Place pseudocount on each RNA transcript: α_RNA,t = s·(1−w)/n_t.

**Verdict:** Rejected. Inflates unexpressed isoforms; user's primary concern.

### Option B — Hierarchical / two-level Dirichlet (recommended)

Decompose the simplex into an *outer split* and an *inner partition*:

- **Outer:**  (θ_gDNA, Σθ_RNA) ~ Beta(α_g, α_r) with α_g = s·w, α_r = s·(1−w)
- **Inner:**  (θ_t / Σθ_RNA)   ~ Dirichlet(0, 0, …, 0)  *(improper, no mass)*

This is *exactly* what we want: calibration controls the gDNA/RNA split; the
distribution within RNA is unrestrained.

**M-step closed form (MAP):**
```
n_g     = posterior expected gDNA mass (from E-step)
n_RNA,t = posterior expected mass for transcript t
S_RNA   = Σ_t n_RNA,t

# 1. Outer split — Beta MAP / VBEM mean (length-corrected analog)
p̂_g = (n_g + α_g) / (n_g + S_RNA + α_g + α_r)

# 2. Inner partition — pure likelihood, no prior mass
q̂_t = n_RNA,t / S_RNA           # if S_RNA > 0 else uniform fallback

# 3. Assemble
θ_g   = p̂_g
θ_t   = (1 − p̂_g) · q̂_t
```

**Properties:**
- **Sparsity preserved:** any transcript with n_RNA,t → 0 goes to θ_t = 0. No floor.
- **Aggregate balance enforced:** even if every isoform individually is unlikely
  given likelihood, the total RNA mass is held up by α_r against gDNA over-pull.
- **VBEM identical:** replace MAP point estimate with mean of Beta(n_g+α_g, S_RNA+α_r),
  same single code path. No Jeffreys baseline needed.
- **Limits are correct:**
  - α_g, α_r → 0  ⇒ pure likelihood EM (no prior).
  - α_r >> n_g  ⇒ θ_g shrinks toward zero (zero-gDNA libraries protected).
  - α_g >> S_RNA ⇒ gDNA captures locus (high-gDNA off-target loci).

### Option C — Transcript-group Dirichlet (extension of B)

Partition transcripts into K groups (by gene, by biotype, or by a "detection
class" learned during calibration). Outer prior becomes (K+1)-dim Dirichlet:

```
(θ_gDNA, Σθ_g1, Σθ_g2, …, Σθ_gK) ~ Dir(α_gDNA, α_g1, …, α_gK)
```

with improper Dirichlet within each group.

**When useful:**
- A gene has known evidence (e.g., spliced reads inside the locus) of being
  expressed, but a paralog gene in the same multi-locus does not. A group-level
  prior reflects that asymmetry without nominating individual isoforms.
- Calibration emits per-gene RNA evidence (which it does, via spliced fragment
  counts in `RegionCalibration`).

**Default group structure:** one group per gene. Single-gene loci collapse to
Option B exactly. Multi-gene loci get differentiated RNA-pool priors.

### Option D — Soft constraint via M-step penalty

Run unpriored EM, then in the M-step project θ onto the constraint

  KL( (θ_g, 1−θ_g) ‖ (w, 1−w) ) ≤ ε(s).

Equivalent to Option B at the optimum (Beta posterior mode) but requires a
projection solver. **Rejected** — Option B is the same math without iteration.

### Option E — Posterior-weighted RNA pseudocount

Distribute α_r across transcripts proportionally to *current* θ_RNA,t each
iteration: α_RNA,t = α_r · θ_RNA,t / Σθ_RNA.

**Verdict:** Functionally similar to Option B (when θ_t = 0 → no prior), but
adaptive priors complicate EM convergence proofs and SQUAREM acceleration.
**Rejected** in favor of B's closed form.

---

## Recommended Plan: Option B (Phase 1), Option C as Phase 2

### Phase 1 — Two-Level Prior (Locus-Aggregate RNA)

**Calibration side (`src/rigel/calibration/prior.py`):**
- Add `mu_rna: np.ndarray` to `RegionCalibration` (already implicit in 4-state E-step).
- Add to `PriorTable`:
  - `gdna_alpha: np.ndarray`  ≡  s · w
  - `rna_alpha:  np.ndarray`  ≡  s · (1 − w)
- Compute `w_l = μ_gDNA_l / (μ_gDNA_l + μ_RNA_l + ε)`.
- Single user knob `em_config.aggregate_prior_strength: float = 1.0` becomes `s`.
- Remove `enable_gdna` as policy gate; keep it only as the technical "is there an
  unspliced unit with finite gDNA log-lik?" flag.

**Native EM side (`src/rigel/native/em_solver.cpp`):**

Replace `compute_gdna_prior_and_warm_start` body (lines ~835–846) with:

```cpp
// Two-level prior: Beta(α_g, α_r) on the gDNA share; improper within RNA.
//
// We store the prior as additive pseudocounts compatible with the existing
// MAP/VBEM update kernels by rewriting only the gDNA component prior. RNA
// components stay at zero pseudocount (Jeffreys baseline removed).
//
// At every M-step, the aggregate Beta MAP is enforced by post-update rescaling:
//   θ_g_new = (n_g + α_g) / (n_g + S_RNA + α_g + α_r)
//   θ_RNA,t_new ∝ n_RNA,t * (1 - θ_g_new)
// This is mathematically equivalent to a per-component prior but avoids
// the per-transcript pseudocount inflation.
prior_out[gdna_idx] = locus_alpha_g;   // s · w
for (int i = 0; i < n_components; ++i) {
    if (i == gdna_idx) continue;
    prior_out[i] = (eligible[i] > 0.0) ? 0.0 : 0.0;   // RNA: zero pseudocount
}

// Warm start: gate ambiguous-row gDNA share by w to prevent zero-gDNA leak.
// In each ambiguous EC row, multiply gDNA's likelihood weight by w before
// the row-normalization. RNA components use weight (1−w)/n_rna_eligible
// as a *blending* factor (does NOT add pseudocount; only redirects ambiguous mass).
```

**The Beta MAP rescaling must be applied at the end of every M-step**, not just
once. Add a new helper invoked from `linked_map_em_step` (or its VBEM analog):

```cpp
static void apply_aggregate_beta_constraint(
    double* theta,
    int n_components,
    int gdna_idx,
    double alpha_g,
    double alpha_r)
{
    if (gdna_idx < 0) return;
    double n_g = theta[gdna_idx];
    double s_rna = 0.0;
    for (int i = 0; i < n_components; ++i)
        if (i != gdna_idx) s_rna += theta[i];
    if (n_g + s_rna <= 0.0) return;

    double p_g_new = (n_g + alpha_g) / (n_g + s_rna + alpha_g + alpha_r);
    double total   = n_g + s_rna;
    theta[gdna_idx] = p_g_new * total;
    double scale   = (s_rna > 0.0) ? ((1.0 - p_g_new) * total / s_rna) : 0.0;
    for (int i = 0; i < n_components; ++i)
        if (i != gdna_idx) theta[i] *= scale;
}
```

**Why this works:**
- Mass conservation: total Σθ is unchanged.
- Within-RNA ratios preserved: q̂_t = n_RNA,t / S_RNA untouched. Zero stays zero.
- gDNA share rescaled to Beta MAP. Smooth in (α_g, α_r).
- SQUAREM unaffected: rescale is a deterministic projection commuting with
  SQUAREM's affine combinations because it preserves the (gDNA, Σ_RNA) ratio
  as a function of one scalar.

**Tests:**
- `tests/test_em_impl.py`: assert θ_g / Σθ converges to Beta MAP at fixed α.
- test_batch_em_impl.py: zero-RNA-mass-in-locus → θ uniform with α_g
  saturating share.
- New `tests/test_aggregate_prior_constraint.py`: synthetic locus with one
  expressed isoform and 20 unexpressed isoforms — verify unexpressed θ_t = 0
  exactly (no floor), gDNA mass = Beta-MAP fraction.

### Phase 2 — Per-Gene Group Priors (Option C)

Extend `PriorTable`:
- `group_idx_of_transcript: np.ndarray[int]`  — group label per transcript.
- `group_alpha: np.ndarray[float]`            — aggregate α for each group.
- Pseudocount on gDNA = α_gDNA; on each group g = α_g; within group = 0.

M-step rescaling generalizes from Beta to Dirichlet over (K+1) aggregates:
identical math, K+1 cells instead of 2.

Default group assignment: one group per `gene_id` within the multi-locus.
Single-gene loci → identical to Phase 1. Multi-gene multi-loci get distinct
per-gene RNA pool priors driven by per-gene calibration evidence.

### User-Facing Knob

A single optional parameter in `EMConfig`:

```python
@dataclass(frozen=True)
class EMConfig:
    ...
    aggregate_prior_strength: float = 1.0  # s in α_g = s·w, α_r = s·(1-w)
```

Limits:
- `s = 0`  → pure likelihood (no prior at all; calibration ignored at EM step).
- `s = 1`  → "one calibrated observation per locus split between gDNA and RNA."
- `s >> 1` → "trust calibration; EM is constrained close to the Beta(α_g, α_r) mode."

Maps naturally to a UI "gDNA sensitivity" slider if desired:
`s = base * 10^(2·(0.5 − sensitivity))`.

---

## Migration Path

1. **PR 1** — Surface `μ_RNA` in `RegionCalibration` (already implicit in 4-state).
2. **PR 2** — Add `gdna_alpha`, `rna_alpha` to `PriorTable`; deprecate
   `gdna_prior_count_em`. Plumb through to native.
3. **PR 3** — Add `apply_aggregate_beta_constraint` to native EM; remove
   `use_vbem ? 0.5 : 0.0` baseline branch; gate warm-start ambiguous shares by w.
4. **PR 4** — Tests (zero-gDNA, zero-RNA, mixed; sparsity preservation in
   multi-isoform loci).
5. **PR 5** — Phase 2: per-gene Dirichlet aggregate priors.

---

## What This Achieves

| Concern | Resolution |
|---|---|
| Zero-gDNA libraries leak gDNA mass | α_r dominates when μ_gDNA → 0; Beta MAP shrinks θ_g toward 0 smoothly. |
| RNA prior inflates unexpressed isoforms | RNA prior is *aggregate* only; per-transcript pseudocount stays at 0. Sparsity preserved. |
| VBEM Jeffreys asymmetry | Baseline 0.5 removed entirely; MAP and VBEM differ only in posterior point estimate (mode vs mean). |
| Warm-start contamination | Ambiguous-row gDNA share gated by w (calibration confidence). |
| User control | One knob `aggregate_prior_strength`; smooth, monotone, principled. |
| Per-gene granularity (future) | Phase 2 Option C drops in cleanly: same M-step rescaling, K+1 aggregates. |

## Open Questions

1. Should `s` itself be **per-locus** (scaled by total locus exposure / unit count)?
   Otherwise small loci with few fragments are dominated by the prior. Probably yes —
   `s_locus = s_global · sqrt(N_observed_unspliced)` is a reasonable default but
   should be benchmarked.
2. Should we expose **per-gene priors** in the public API even in Phase 1, gated
   behind a flag, so calibration→EM contract is forward-compatible?
3. Interaction with SQUAREM acceleration step length: the aggregate constraint
   is affine in the (gDNA, Σ_RNA) marginal but nonlinear in the full simplex.
   Need empirical check that SQUAREM convergence rate is preserved.
