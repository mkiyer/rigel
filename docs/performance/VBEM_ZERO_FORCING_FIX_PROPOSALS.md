# VBEM Zero-Forcing Fix: Revised Proposals

**Date:** April 3, 2026  
**Status:** Planning — not yet implemented  
**Prerequisite:** Read [VCAP_VBEM_ANALYSIS.md](VCAP_VBEM_ANALYSIS.md) for benchmark findings  
**Related:** [vbem_convergence_plan.md](vbem_convergence_plan.md) (the original plan that introduced zero-forcing)

## Problem Recap

The VBEM convergence fix from `vbem_convergence_plan.md` introduced two changes:

1. **Change 1 (SQUAREM clamp):** Clamp extrapolated alpha to `prior[i]` instead of `EM_LOG_EPSILON`
2. **Change 2 (zero-forcing):** Collapse components to prior when `alpha_new ≤ prior * (1 + 1e-6)`

Together, these **successfully fix convergence** — VBEM no longer hits the iteration cap. But they **catastrophically degrade accuracy** on the VCaP benchmark (Pearson R drops from 0.986 to 0.815) by incorrectly zero-forcing highly expressed ribosomal protein genes in the mega-locus.

### Why the Original Plan's Assumptions Were Wrong

The plan stated: *"any component with genuine support has accumulated evidence far exceeding `prior × tol` within the first few iterations."*

This assumption fails because:

1. **The prior is vanishingly small.** With `C = 1.0` spread coverage-weighted over ~633K components, each RNA prior is ~$1.6 \times 10^{-6}$. The zero-forcing threshold `prior × (1 + 10^{-6})` is effectively `prior` itself — any component whose M-step output equals `prior + ε` passes, but any that rounds to `prior` is killed.

2. **SQUAREM triggers the absorbing barrier.** The real damage isn't from zero-forcing *per se* — it's from the combination of SQUAREM overshooting → clamping to prior → digamma creating an absorbing barrier:
   - SQUAREM extrapolation pushes a competing component's alpha below prior
   - Clamped to `prior ≈ 1.6e-6` (Change 1)
   - Stabilization E-step: `digamma(1.6e-6) ≈ -625,000` → effective weight `exp(-625,000) = 0`
   - M-step: `alpha_new = 0 + 0 + prior = prior`
   - Zero-forcing (Change 2): confirms `alpha ≤ prior * (1 + tol)` → locked at prior
   - **Even without Change 2**, the component is dead — `digamma(prior)` is so extreme that `em_totals = 0` self-reinforces

3. **The absorbing barrier is the real bug, not zero-forcing alone.** Once a component reaches `alpha ≈ prior ≈ 1.6e-6`, `digamma(1.6e-6) ≈ -625,000` makes recovery impossible regardless of whether we explicitly zero-force. Change 1 creates the absorbing barrier; Change 2 just makes it trigger faster.

### The Absorbing Barrier Mechanism

For a component at alpha value $\alpha$, the VBEM E-step weight is proportional to:

$$w_i \propto \exp\!\bigl(\psi(\alpha_i) - \psi(\Sigma) + s_i\bigr)$$

where $\psi$ is digamma, $\Sigma = \sum_j \alpha_j$, and $s_i$ is the score. For small $\alpha$:

$$\psi(\alpha) \approx -\frac{1}{\alpha} - \gamma_E \quad (\alpha \to 0^+)$$

| $\alpha$ | $\psi(\alpha)$ | $\exp(\psi(\alpha))$ | Recoverable? |
|-----------|----------------|-----------------------|--------------|
| 100 | 4.6 | 99.5 | ✓ Healthy |
| 1.0 | −0.58 | 0.56 | ✓ Marginally |
| 0.1 | −10.4 | $3 × 10^{-5}$ | ✓ Weak but nonzero |
| 0.01 | −100.6 | $3 × 10^{-44}$ | ✗ Effectively dead |
| $10^{-6}$ | $-10^{6}$ | 0 | ✗ Absorbing barrier |

**Recovery threshold:** A component needs $\alpha \gtrsim 0.1$ for `digamma` to yield a weight large enough to accumulate meaningful `em_totals`. Below ~0.01, the component is effectively dead in double precision.

**The current clamp floor (prior ≈ $10^{-6}$) is 4–5 orders of magnitude below the recovery threshold.** Any component pushed to this floor by SQUAREM is permanently eliminated, whether or not we explicitly zero-force it.

### Why MAP-EM Doesn't Have This Problem

MAP-EM uses $\log(\theta_i)$ instead of $\psi(\alpha_i)$. For two competing components:

- MAP: $\log(\theta)$ is smooth, monotone, and its gradient $1/\theta$ is gentle compared to digamma's $-1/\alpha^2$ near zero
- When $\theta$ gets small, $\log(\theta)$ decreases logarithmically — components retain nonzero weight and can recover
- MAP's normalization step ($\theta_i = \alpha_i / \Sigma$) ensures proportional sharing, not winner-take-all

For components with $\alpha \gg 1$, VBEM and MAP behave identically ($\psi(\alpha) \approx \log(\alpha)$ for large $\alpha$). The pathology is exclusively in the $\alpha \ll 1$ regime — exactly where SQUAREM clamping sends components.

---

## Proposal 1: Recoverable Clamping Floor (Recommended)

### Concept

Replace the prior-based clamp with a floor that keeps components in the *recoverable* regime of digamma. Instead of clamping SQUAREM extrapolation to `prior[i]` (which is in the absorbing regime), clamp to a value where `digamma(floor)` produces a small but nonzero weight.

### Changes

**Change 1 (revised):** Clamp SQUAREM extrapolation to `max(prior[i], recoverable_floor)`:

```cpp
static constexpr double VBEM_CLAMP_FLOOR = 0.1;

// In SQUAREM extrapolation loop:
double floor_i = std::max(prior[i], VBEM_CLAMP_FLOOR);
if (state_extrap[i] < floor_i)
    state_extrap[i] = floor_i;
```

**Change 2: Remove zero-forcing entirely.**

### Why VBEM_CLAMP_FLOOR = 0.1?

- $\psi(0.1) \approx -10.42$, so $\exp(\psi(0.1)) \approx 3 \times 10^{-5}$
- A component at the floor has weight $\sim 10^{-5}$ relative to a healthy component — small enough not to steal significant mass, but large enough to accumulate evidence and recover if it has genuine support
- For a ghost component with no EC entries: `em_totals = 0` always, so `alpha_new = 0 + 0 + prior = prior`. Next iteration, SQUAREM sees it stable → no extrapolation → clamped to floor → stable at floor. Ghost contributes exactly zero to L1 norm.
- For a "loser" component temporarily pushed down by SQUAREM: `em_totals > 0` (tiny but nonzero), so `alpha_new = tiny + prior`. Over subsequent iterations, if it has genuine read support, alpha grows back above floor → recovery.

### Steady-State Analysis

**Ghost component (no EC support):**
1. `alpha = floor = 0.1` → `digamma(0.1) = -10.42` → weight ≈ $3 × 10^{-5}$
2. But the component has no EC entries → `em_totals = 0`
3. `alpha_new = 0 + 0 + prior ≈ prior` → clamped back to floor
4. Stable at floor. L1 contribution = 0 per iteration.

**Ghost component that appears in ECs but has no real support:**
1. `alpha = floor = 0.1` → `digamma(0.1) = -10.42` → weight ≈ $3 × 10^{-5}$
2. The component appears in some ECs. Its tiny weight gives it a proportional share: `em_totals ≈ 3 × 10^{-5} × (fragment share)`
3. For a component competing against others with alpha >> 1: `em_totals` is negligible
4. `alpha_new ≈ em_totals + prior ≈ prior` → clamped to floor
5. Stable near floor. L1 contribution ≈ 0.

**Legitimate component pushed to floor by SQUAREM:**
1. `alpha = floor = 0.1` → weight ≈ $3 × 10^{-5}$
2. The component has genuine support in ECs
3. `em_totals > 0` (proportional to its weight × fragment share)
4. If it has, say, 1000 shared fragments: `em_totals ≈ 3 × 10^{-5} × 1000 / k_competitors ≈ 0.03`
5. `alpha_new = 0.03 + prior ≈ 0.03` → clamped to floor = 0.1
6. Next iteration: still at floor, gets another injection of evidence
7. Over time, evidence accumulates: `alpha = floor` each iteration, but the component stays in the game
8. If the competitors are also near the floor, proportional sharing occurs → correct partitioning

Wait — step 5 is key. If `alpha_new = 0.03 < floor`, it's clamped back to 0.1. But if `alpha_new = 0.5 > floor`, it escapes. Recovery speed depends on how much evidence the component accumulates per iteration from the floor.

### Alternative Floors

| Floor | $\psi(\text{floor})$ | $\exp(\psi)$ | Regime |
|-------|----------------------|--------------|--------|
| 0.5 | −1.96 | 0.14 | Strong recovery — component keeps ~14% relative weight |
| 0.1 | −10.42 | $3 × 10^{-5}$ | Weak recovery — needs many rounds of small evidence |
| 0.01 | −100.6 | $3 × 10^{-44}$ | Dead — no recovery in finite precision |

**Discussion:** A floor of 0.1 may be *too* conservative for fast recovery. Consider 0.5 or even 1.0. At floor = 1.0, `digamma(1.0) = -0.577`, `exp(-0.577) = 0.56` — the component retains ~56% relative weight and can compete fully. This is essentially saying "when SQUAREM overshoots, assume the component has ~1 virtual count."

A floor of 1.0 might over-protect ghost components (giving them too much weight), but since em_totals for ghosts is still ≈ 0, they won't steal significant mass. The main effect is on convergence speed — more components with nonzero weight means more "noise" in the E-step, potentially requiring more iterations.

### Recommendation

Start with `VBEM_CLAMP_FLOOR = 0.1` and test. If recovery is too slow, increase to 0.5 or 1.0. The key point: **any value ≥ 0.01 is infinitely better than the current floor of $10^{-6}$.**

### Tradeoffs

| Pro | Con |
|-----|-----|
| Simple — one constant change, one deletion | No Bayesian interpretation of the floor |
| Eliminates absorbing barrier | Ghost components with EC entries retain tiny but nonzero weight |
| Preserves SQUAREM acceleration | May need tuning (0.1 vs 0.5 vs 1.0) |
| No new data structures or passes | Convergence might be slightly slower than zero-forcing |

---

## Proposal 2: Ghost-Only Zero-Forcing

### Concept

The convergence problem comes from *ghost* components (those with no evidence), not from competing components. Pre-identify which components have structural support (appear in at least one EC) and only zero-force components with zero EC support.

### Changes

1. Before the SQUAREM loop, scan EC data to build a support mask:
```cpp
std::vector<bool> has_ec_support(nc, false);
for (const auto& ec : ec_data) {
    for (int j = 0; j < ec.k; ++j)
        has_ec_support[ec.comp_idx[j]] = true;
}
```

2. In the zero-forcing check:
```cpp
if (!has_ec_support[i]
    && state_new[i] <= prior[i] * (1.0 + VBEM_ZERO_FORCE_REL_TOL)) {
    state_new[i] = prior[i];
}
```

3. Keep Change 1 (clamp to prior).

### Why It Works

- Ghost components (no EC entries) → `em_totals = 0` always → oscillation source → zero-forcing eliminates the oscillation
- Real components (have EC entries) → never zero-forced → free to compete → correct mass partitioning
- The absorbing barrier from Change 1 still exists for components WITH EC support that SQUAREM pushes to prior, but these components have structural support and will tend to recover (they appear in ECs, so they get some em_totals every iteration)

### Problem: The Absorbing Barrier Still Exists

Even with ghost-only zero-forcing, Change 1 creates an absorbing barrier at `prior ≈ 1.6e-6` for ALL components. If SQUAREM pushes a real component (one with EC support) below prior, it gets clamped to prior, where `digamma(prior)` kills its weight. The fact that we don't *explicitly* zero-force it doesn't help — `digamma(1.6e-6) = -625,000` makes it dead anyway.

**This proposal must be combined with a raised clamp floor (Proposal 1) to be effective.** Ghost-only zero-forcing alone, with clamping to prior, has the same absorbing barrier problem.

### Combination: Proposals 1 + 2

The strongest variant:
- Clamp SQUAREM extrapolation to `max(prior[i], VBEM_CLAMP_FLOOR)` for all components
- Additionally zero-force components with no EC support (for faster convergence)
- Real components can recover; ghosts are cleanly eliminated

### Tradeoffs

| Pro | Con |
|-----|-----|
| Targeted — only zero-forces true ghosts | Requires O(EC) pre-scan (one-time, cheap) |
| Preserves competition dynamics for real components | Extra boolean array per locus |
| Fast convergence — ghosts eliminated immediately | Must be combined with Proposal 1 to avoid absorbing barrier |

---

## Proposal 3: Deferred Zero-Forcing (Burn-In)

### Concept

Don't zero-force for the first $B$ iterations. Let the posterior stabilize before making irreversible decisions. After burn-in, zero-force as currently implemented.

### Implementation

```cpp
// In SQUAREM loop:
bool allow_zero_force = (iter >= burn_in_iters);

if (allow_zero_force
    && state_new[i] <= prior[i] * (1.0 + VBEM_ZERO_FORCE_REL_TOL)) {
    state_new[i] = prior[i];
}
```

### Choosing B

| Strategy | Formula | Mega-locus value | Note |
|----------|---------|------------------|------|
| Fixed | $B = 50$ | 50 | Simple heuristic |
| Fraction of budget | $B = \text{max\_iters} / 4$ | 83 | Scales with iteration budget |
| Adaptive | Wait for $\delta < 10 \times \delta_{\text{conv}}$ | ? | Convergence-based |

### Analysis

After $B$ iterations, the posterior has had time to partition mass among competitors. Components with genuine support have accumulated significant evidence ($\alpha \gg \text{prior}$). Components with no support are at or near prior. Zero-forcing after burn-in should only kill true ghosts.

**Problem:** The absorbing barrier from Change 1 can still trap components *during* burn-in. If SQUAREM pushes RPL18A below prior in iteration 15 (within burn-in), it gets clamped to prior = $1.6 \times 10^{-6}$, `digamma` kills it, and it's dead by the time burn-in ends — even without zero-forcing. **This proposal must also be combined with a raised clamp floor.**

### Tradeoffs

| Pro | Con |
|-----|-----|
| Gives posterior time to stabilize | Still has absorbing barrier during burn-in |
| Simple to implement | B is a new hyperparameter |
| Reduces false zero-forcing | Doesn't address the root cause |

---

## Proposal 4: Adaptive Prior Budget

### Concept

Scale the total prior pseudocount $C$ with locus size so each component has a meaningful prior — one large enough that `digamma(prior)` doesn't create an absorbing barrier.

### Implementation

```cpp
double effective_C = std::max(
    total_pseudocount,
    n_rna_eligible * MIN_PER_COMPONENT_PRIOR
);
double rna_budget = (1.0 - locus_gamma) * effective_C;
```

### Choosing MIN_PER_COMPONENT_PRIOR

For no absorbing barrier, need `prior_i ≥ 0.01` per component (conservatively). Then:

| Locus size | $C$ needed | Effect |
|------------|-----------|--------|
| 100 components | max(1.0, 1.0) = 1.0 | No change |
| 10K components | max(1.0, 100) = 100 | 100× stronger prior |
| 633K components (mega-locus) | max(1.0, 6330) = 6,330 | 6330× stronger prior |

With `MIN_PER_COMPONENT_PRIOR = 1e-3`:

| Locus size | $C$ | Prior per component |
|------------|-----|-------------------|
| 100 | 1.0 | 0.005 |
| 10K | 10.0 | 0.0005 |
| 633K | 633 | $5 \times 10^{-4}$ |

Still below the recovery threshold for the mega-locus. Would need `MIN_PER_COMPONENT_PRIOR ≥ 0.01`.

### Analysis

**This changes the Bayesian model.** A prior of $C = 6330$ means we're asserting 6330 "virtual reads" of prior belief, vs the actual 61M real reads in the mega-locus. As a fraction, $6330/61M = 0.01\%$ — this is still a very weak prior and unlikely to bias results.

For small loci (say 50 transcripts, 5000 reads), keeping $C = 1.0$ preserves the current behavior.

**Key advantage:** This addresses the root cause directly. If each component has `prior ≈ 0.01`, then:
- `digamma(0.01) ≈ -100` — still strongly penalized but not absorbing
- Clamping to prior doesn't create a dead zone
- Zero-forcing at `prior × (1 + tol)` has a meaningful threshold

**Key disadvantage:** It reduces VBEM's sparsification in large loci. But the VCaP benchmark shows that VBEM's sparsification is *already harmful* in mega-loci — it zero-forces the *wrong* components. Reducing sparsification in mega-loci is a feature, not a bug.

### Tradeoffs

| Pro | Con |
|-----|-----|
| Addresses root cause (prior too small) | Changes the Bayesian model for large loci |
| No new mechanisms, just a scaled constant | Requires choosing MIN_PER_COMPONENT_PRIOR |
| Compatible with existing zero-forcing | May reduce desirable sparsification |
| Small loci unaffected | Prior is no longer "1 virtual read"semantically |

---

## Proposal 5: MAP → VBEM Two-Stage

### Concept

Use MAP-EM (which correctly handles the mega-locus) to identify the active component set. Then run VBEM only on the active set.

### Implementation

1. **Stage 1:** Run MAP-EM to convergence → obtain $\theta^{MAP}$
2. **Define active set:** $\mathcal{A} = \{i : \theta^{MAP}_i > \epsilon\}$ (e.g., $\epsilon = 10^{-8}$)
3. **Rebuild ECs:** Restrict equivalence classes to components in $\mathcal{A}$
4. **Stage 2:** Run VBEM on the reduced problem with the active set

### Analysis

**MAP identifies the truth well** (Pearson R = 0.986). VBEM on a reduced active set (~100K components instead of 633K) would have:
- Much larger per-component prior: $1.0 / 100K = 10^{-5}$ — still small but much better
- Fewer competitors for each read → faster convergence
- Reduced risk of SQUAREM overshooting into the absorbing regime

**VBEM adds value on the reduced set** through posterior uncertainty quantification and principled shrinkage, without the pathological sparsification that occurs in the full problem.

### Warm-Start Concern

The transition from MAP $\theta$ to VBEM $\alpha$ requires care:
- Initial $\alpha_i = \theta^{MAP}_i \times N_{frag}$ (scale by total fragment count)
- This gives a reasonable starting point for the Dirichlet posterior
- Zero components from MAP ($\theta = 0$) are excluded from the VBEM active set

### Tradeoffs

| Pro | Con |
|-----|-----|
| Most robust — MAP handles competition correctly | ~2× total EM runtime |
| VBEM on reduced set is well-conditioned | Added code complexity |
| Guaranteed accuracy ≥ MAP | Need to define active set threshold |
| Posterior uncertainty from VBEM | MAP→VBEM warm-start transition is non-trivial |

---

## Proposal 6: Per-Locus Mode Selection

### Concept

Use VBEM for small/medium loci where it provides useful sparsification and uncertainty. Use MAP-EM for mega-loci where VBEM's sparsification is pathological.

### Implementation

```cpp
bool use_vbem = (mode == "vbem") && (n_components < MEGA_LOCUS_THRESHOLD);
```

With `MEGA_LOCUS_THRESHOLD = 10000` or similar.

### Analysis

This is pragmatic but unprincipled. It acknowledges that VBEM and MAP have different failure modes:
- VBEM: excellent for small loci, pathological for mega-loci
- MAP: robust everywhere, loses VBEM's uncertainty benefits for small loci

The threshold would need to be determined empirically. The VCaP mega-locus has 424K components; the next-largest locus has 118. There's likely a natural gap.

### Tradeoffs

| Pro | Con |
|-----|-----|
| Simple, pragmatic | Unprincipled — no reason VBEM shouldn't work on large loci in principle |
| Preserves VBEM benefits for small loci | Threshold is a new parameter |
| MAP handles mega-locus correctly | Discontinuity at threshold |

---

## Comparison Matrix

| Proposal | Fixes convergence | Fixes accuracy | Addresses root cause | Complexity | New parameters |
|----------|:-:|:-:|:-:|:-:|:-:|
| 1. Recoverable floor | ✓ (likely) | ✓ | Partially | Low | VBEM_CLAMP_FLOOR |
| 2. Ghost-only ZF | ✓ | ✗ (alone) | No | Low | None |
| 3. Deferred ZF | ✓ (after burn-in) | ? | No | Low | burn_in_iters |
| 4. Adaptive prior | ✓ | ✓ | **Yes** | Low | MIN_PER_COMPONENT_PRIOR |
| 5. MAP→VBEM | ✓ | ✓ | Yes (by avoidance) | High | active_set_threshold |
| 6. Per-locus mode | ✓ (for MAP loci) | ✓ (for MAP loci) | No (avoidance) | Low | MEGA_LOCUS_THRESHOLD |

---

## Recommended Strategy

### Primary: Proposal 1 + 2 Combined (Recoverable Floor + Ghost-Only Zero-Forcing)

**Rationale:** This combination attacks the problem from both directions:

1. **Raised clamp floor** eliminates the absorbing barrier — SQUAREM can push components down, but they land at $\alpha = 0.1$ where `digamma(0.1) = -10.4` gives a small but recoverable weight. Legitimate components escape; undeserving ones are penalized but not killed.

2. **Ghost-only zero-forcing** provides fast convergence by immediately eliminating components with zero structural support (no EC entries). These are true ghosts that can never accumulate evidence.

Together:
- Ghosts: zero-forced immediately → L1 contribution = 0 → fast convergence ✓
- Losers WITH EC support: land at floor, retain tiny weight, can recover → no false zeros ✓
- Winners: unaffected (alpha >> floor) → correct quantification ✓

### Validation Plan

1. Implement Proposal 1 + 2
2. Recompile: `pip install --no-build-isolation -e .`
3. Run existing tests: `pytest tests/ -v`
4. Re-run VCaP benchmark with both `vbem` and `map` configs
5. **Success criteria:**
   - VBEM Pearson R ≥ 0.98 (vs current 0.815)
   - VBEM WARE ≤ 0.10 (vs current 0.252)
   - VBEM converges (doesn't hit iteration cap)
   - No ribosomal protein genes zero-forced
   - Runtime within 1.5× of MAP

### Fallback: Proposal 4 (Adaptive Prior)

If Proposal 1+2 doesn't fully resolve the accuracy issue (e.g., components at floor = 0.1 still accumulate enough mass to distort results), fall back to Proposal 4. Scaling the prior budget ensures every component has a meaningful prior, making the Dirichlet posterior well-conditioned even in the mega-locus.

The fallback could also be *combined* with Proposal 1+2: scale the prior AND raise the clamp floor AND ghost-only zero-force.

### Future: Proposal 5 (MAP→VBEM Two-Stage)

If VBEM continues to be problematic in mega-loci despite the above fixes, the MAP→VBEM two-stage approach is the nuclear option. It's the most robust but most expensive. Worth implementing if the single-pass VBEM fixes don't achieve MAP-level accuracy.

---

## Appendix: Detailed Digamma Analysis

### SQUAREM Overshoot Scenario

Consider RPL18A competing with pseudogene RPL18AP1. Both share 5000 reads.

**Iteration K:**
- RPL18A: $\alpha = 3000$, RPL18AP1: $\alpha = 2500$
- SQUAREM acceleration factor $step = 3.0$
- Extrapolation: $\alpha_{RPL18A}^{ext} = 3000 + 2(3)(-100) + 9(20) = 3000 - 600 + 180 = 2580$
- Both still healthy. No clamp triggered.

**Iteration K+5 (pseudogene winning):**
- RPL18A: $\alpha = 500$, RPL18AP1: $\alpha = 5000$
- r_vec for RPL18A: $-200$ (trending down strongly)
- v_vec: $+50$ (decelerating)
- $\alpha_{RPL18A}^{ext} = 500 + 2(3)(-200) + 9(50) = 500 - 1200 + 450 = -250$

**With current Change 1** (clamp to `prior ≈ 1.6e-6`):
- Clamped to $1.6 \times 10^{-6}$
- Stabilization: `digamma(1.6e-6) ≈ -625,000` → weight = 0 → dead
- Zero-forcing confirms: $\alpha_{new} = prior$ → permanent death
- **RPL18A incorrectly killed. Pseudogene inherits all 5000 reads.**

**With Proposal 1** (clamp to `floor = 0.1`):
- Clamped to 0.1
- Stabilization: `digamma(0.1) ≈ -10.4` → weight ≈ $3 \times 10^{-5}$
- em_totals: $3 \times 10^{-5} \times 5000 / 2 \approx 0.075$ (from shared reads)
- $\alpha_{new} = 0.075 + 1.6 \times 10^{-6} \approx 0.075$
- Clamped to floor = 0.1 for next SQUAREM iteration
- **RPL18A survives. Over subsequent iterations, if it has genuine support, it can compete back for shared reads. If the pseudogene truly has no independent support, RPL18A will eventually win.**

### Coverage-Weighted Prior Distribution

The prior isn't uniform — it's proportional to coverage (total EC weight). For the mega-locus:

```
rna_budget = (1 - 0.486) × 1.0 = 0.514
total_rna_coverage ≈ 61,000,000
prior_i = 0.514 × coverage_i / 61,000,000
```

A transcript with 1000 reads of coverage: `prior = 0.514 × 1000 / 61M ≈ 8.4 × 10^{-6}`
A transcript with 100,000 reads: `prior = 0.514 × 100,000 / 61M ≈ 8.4 × 10^{-4}`

Even the most highly covered transcript has `prior < 0.001` — well below the recovery threshold. **Coverage-weighting doesn't save us; the total budget $C = 1.0$ is simply too small for a 633K-component locus.**
