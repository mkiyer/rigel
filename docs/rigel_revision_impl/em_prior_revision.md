## 1. What's Happening Now (Current Architecture)

### The Prior is NOT What You Think

The Python prior array built in locus.py is **only a binary eligibility mask**. Here's the flow:

```
Python locus.py: prior = [1e-10, 1e-10, ..., 1e-10]  (all ε)
  ↓ zero prior[gdna] if gdna_init == 0
  ↓ zero prior[nrna] if all-single-exon
  ↓
C++ em_solver.cpp: eligible[i] = (prior[i] > 0) ? 1.0 : 0.0   ← BINARY MASK
  ↓
compute_ovr_prior_and_warm_start() OVERWRITES the prior:
  prior[i] = alpha_flat + gamma × coverage[i] / n_ambig   if eligible
  prior[i] = 0.0                                          if dead
```

The actual prior values entering the M-step are **always** the OVR-computed values — never the Python epsilon. The Python side only decides which components are alive or dead.

### The M-Step Formula

In MAP-EM (em_solver.cpp):

$$\theta_i^{\text{new}} = \frac{u_i + e_i + p_i}{\sum_j(u_j + e_j + p_j)}$$

where $u_i$ = unambiguous counts, $e_i$ = E-step assigned counts, $p_i$ = OVR prior.

This is the **Dirichlet posterior mean** (not the MAP/mode). It treats the prior as additive pseudocounts — every eligible component always gets some mass. **There is no mechanism to drive a component to zero once it's eligible.**

### Why nRNA FPs Are 83%

With default settings (`alpha_flat=0.01`, `gamma=1.0`):
- Every eligible nRNA component gets `prior[nrna] ≈ 0.01 + OVR_share`
- The OVR distributes ambiguous-fragment coverage proportionally — nRNA and gDNA components that can explain the data will attract some of these virtual reads
- In the M-step, even with zero unambiguous evidence, the OVR prior + leaked E-step mass produces non-zero theta
- The nRNA and gDNA components are **anti-correlated** (r = −0.997): fragments bounce between them. Whichever one the OVR seeds more heavily captures the leakage.

Bottom line: the current prior says "every eligible component should exist" — it can only *kill* components (prior=0), never *discourage* them. There's no graded skepticism.

---

## 2. The Dirichlet Sparsity Prior (α < 1)

### The Core Idea

A standard Dirichlet prior $\text{Dir}(\alpha_1, \dots, \alpha_K)$ has two regimes:

| $\alpha_i$ | Effect | Interpretation |
|---|---|---|
| $> 1$ | Concentrating — pushes mass toward the center | "All components should probably be present" |
| $= 1$ | Flat / uniform | "No opinion" |
| $< 1$ | **Sparsifying** — pushes mass toward corners | "Most components should be zero; only a few should be large" |

### Where It Fits: Replace the OVR Prior for nRNA

Currently the OVR prior treats all eligible components identically. The change:

| Component | Current Prior | New Prior |
|---|---|---|
| **mRNA** | $\alpha_{\text{flat}} + \text{OVR}_i$ | **Unchanged** — keep OVR |
| **nRNA** | $\alpha_{\text{flat}} + \text{OVR}_i$ | $\alpha_{\text{nrna}}$ (e.g. 0.9) — **replace** OVR |
| **gDNA** | $\alpha_{\text{flat}} + \text{OVR}_i$ | $\alpha_{\text{gdna}}$ from EB — **replace** OVR |

The prior is **not added** to the existing system — it **replaces** the OVR contribution for non-mRNA components. Think of it as: mRNA keeps its data-driven prior, while nRNA and gDNA get informative priors that encode our skepticism (nRNA) or regional estimate (gDNA).

### How It Crushes FPs: The Toll Mechanism

**MAP-EM formulation** — the true Dirichlet MAP (mode of the posterior) is:

$$\theta_i^{\text{MAP}} \propto \max\!\bigl(0,\; N_i + \alpha_i - 1\bigr)$$

For $\alpha_{\text{nrna}} = 0.9$:
- **Toll** = $1 - \alpha = 0.1$ fragments
- Any nRNA component where the total evidence (unambig + E-step) is less than 0.1 fragments gets **hard-zeroed**
- The M-step clamps $\theta_i = 0$, and the next E-step assigns it zero posterior weight
- Once zeroed, it stays zeroed — this is the sparsity "snap"

Concrete example with your sweep data:
- At SS=1.0 with NTA1=0 (no true nRNA): the spurious nRNA mean is 67.1 fragments across the whole locus, but this is spread across multiple nRNA components. If each component gets <0.1, they all snap to zero. Even if one gets 0.2 due to leakage, only that one survives — compared to all of them surviving today.

**VBEM formulation** (already implemented in your code):

$$\log w_i = \psi(\alpha_i) - \psi\!\left(\textstyle\sum_j \alpha_j\right) - \log L_i$$

For $\alpha = 0.9$: $\psi(0.9) \approx -1.07$ (negative!). This creates a **logarithmic penalty** in the E-step that a component must overcome with likelihood evidence. It's a softer version of the same sparsity — components don't hard-zero but become exponentially unlikely without strong data support.

### MAP-EM vs VBEM Trade-offs

| | MAP-EM with α−1 clamping | VBEM (already implemented) |
|---|---|---|
| Sparsity mechanism | Hard zero via max(0, ...) | Soft penalty via digamma |
| Code change | Need to modify M-step formula and add clamping | Just change prior values; solver works as-is |
| SQUAREM stability | Discontinuity at zero may cause extrapolation issues | Smooth everywhere |
| Existing support | **Not yet implemented** — current M-step uses posterior mean | **Already implemented** — just needs α < 1 |

**Recommendation**: Use VBEM. It already handles Dirichlet priors correctly, requires zero changes to em_solver.cpp, and is numerically smoother for SQUAREM. You'd only need to change the prior construction + set `mode="vbem"`.

---

## 3. The Spatial EB gDNA Prior

### What You Already Have

`compute_eb_gdna_priors()` in locus.py is a sophisticated three-level hierarchical estimator:

$$G_{\text{locus}} = w \cdot G_{\text{observed}} + (1-w) \cdot G_{\text{ref-level}}$$

where the ref-level estimate is itself shrunk toward a global estimate using MoM-estimated κ. The output `gdna_init` is in **fragment count** units (density × exonic_bp).

Currently this carefully computed value is used **only for binary gating**: if `gdna_init > 0`, the gDNA component is eligible; otherwise it's dead. The magnitude is thrown away.

### How to Use It as a Proper Prior

Instead of:
```python
prior[gdna_idx] = EM_PRIOR_EPSILON if gdna_init > 0 else 0.0
```

Use:
```python
prior[gdna_idx] = gdna_init  # fragment count as Dirichlet α
```

This tells the EM: "We expect approximately `gdna_init` fragments of gDNA at this locus, based on the regional background." The EB shrinkage ensures this is neither overconfident (raw locus estimate) nor uninformative (global average).

### Why This Fixes gDNA-nRNA Confusion

The current nRNA-gDNA correlation is −0.997 because both components compete for the same ambiguous fragments with **no anchoring** — either can absorb the leakage. With an informative gDNA prior:

1. The gDNA component is anchored to a regional expectation
2. It can still deviate if the likelihood strongly supports more/less gDNA
3. But it can't undergo the wild swings that produce the 190% gDNA overestimate at SS=0.50

Combined with α < 1 for nRNA, you get complementary forces:
- nRNA: sparsity prior pushes toward zero unless justified
- gDNA: informative prior anchors to regional background
- Result: the ambiguous fragment mass distributes according to the *informed* prior rather than bouncing between two equally-uninformed attractors

---

## 4. What to Keep, Change, and Remove

### Keep (no changes)

| What | Why |
|---|---|
| `log(0.5)` strand penalty for gDNA | Correct generative model; already in scoring.cpp |
| OVR warm-start for mRNA | Good multi-mapper resolution for expressed transcripts |
| Binary gating for single-exon nRNA | Sound: single-exon transcripts can't have intronic nRNA |
| EB gdna_init computation | Already good infrastructure; just need to use the value properly |
| nrna_init computation | The strand-corrected intronic evidence is useful signal |
| SQUAREM acceleration | Already working |
| The entire C++ EM solver | No M-step changes needed if using VBEM |

### Change

| What | From | To |
|---|---|---|
| nRNA prior | `alpha_flat + OVR_share` | `nrna_sparsity_alpha` (e.g. 0.9) — component-type-specific |
| gDNA prior | `alpha_flat + OVR_share` | `gdna_init` (EB fragment count) |
| EM mode | `"map"` (default) | `"vbem"` (default) — enables Dirichlet sparsity |
| `compute_ovr_prior_and_warm_start()` | Uniform treatment of all components | Apply OVR only to mRNA; nRNA/gDNA get their own α values |
| EMConfig | Single `prior_alpha` | Add `nrna_prior_alpha: float = 0.9` and `gdna_prior_mode: str = "eb"` |

### Remove (optional cleanup)

| What | Why |
|---|---|
| The OVR coverage contribution to nRNA/gDNA | Replaced by informative priors |
| `EM_PRIOR_EPSILON` as the Python-side prior | The Python prior array should carry the actual intended α values, not just a binary mask |

---

## 5. Concrete Next Steps

### Step 1: Make the Prior Component-Type-Aware (Python only)

Modify `compute_ovr_prior_and_warm_start()` in em_solver.cpp to accept per-component prior α values instead of a single `alpha_flat`. The Python side builds the full prior array:

```python
# In locus.py or estimator.py:
prior = np.full(n_components, alpha_flat, dtype=np.float64)  # mRNA: unchanged
for ln in range(n_nrna):
    if eligible[n_t + ln]:
        prior[n_t + ln] = nrna_sparsity_alpha   # e.g. 0.9
    else:
        prior[n_t + ln] = 0.0
prior[gdna_idx] = gdna_init if gdna_init > 0 else 0.0

# OVR only adds coverage to mRNA components
```

Then pass this pre-built prior into the solver, bypassing the OVR for non-mRNA.

### Step 2: Switch Default to VBEM

Change `EMConfig.mode` default from `"map"` to `"vbem"`. The existing VBEM implementation handles Dirichlet α < 1 correctly via digamma — no C++ changes needed.

### Step 3: Validate with the Sweep

Re-run the sweep with `nrna_sparsity_alpha` ∈ {0.5, 0.7, 0.9, 0.95} to find the optimal value. The key metrics to watch:
- nRNA FP rate (target: <10%, currently 83%)
- nRNA recall at SS=1.0 (shouldn't drop much)
- gDNA accuracy at SS=0.50 (should improve with anchoring)
- mRNA accuracy (should be unaffected)

### Step 4: Tune gDNA EB Prior Strength

The gdna_init value might need scaling. If the EB estimate says 10 fragments, using α_gdna=10 is a strong prior. You may want a scale factor: `α_gdna = scale × gdna_init` where `scale` is tuned via the sweep.

### Implementation Difficulty

- Step 1: Medium — modify the C++ OVR function to accept a pre-built prior array (or partition it by component type), plus Python-side prior construction changes
- Step 2: Trivial — one-line default change
- Step 3-4: Already have the infrastructure

Do you want me to start implementing Step 1 — making the prior component-type-aware? 

Completed: *Produce comprehensive explanation* (4/4)




Here is the concise, implementation-ready blueprint for the Tripartite Prior architecture. This completely decouples your prior policy (Python) from the EM engine (C++).

### Phase 1: Configuration Updates (Python)

Update your `EMConfig` class to support the new prior parameters and switch the default solver mode to ensure mathematical stability.

1. Add `nrna_sparsity_alpha: float = 0.9` (Must be $< 1.0$).
2. Add `gdna_prior_scale: float = 1.0` (To tune the strength of the EB prior).
3. Change default EM mode: `mode: str = "vbem"` (Crucial for SQUAREM stability when $\alpha < 1$).

### Phase 2: Prior Array Construction (Python)

In `locus.py` (or wherever you prepare the data for the C++ solver), build the complete `prior` array explicitly for all $T + N + 1$ components before passing it to C++.

1. **Initialize Array:** Create a float array of size `T + N + 1`.
2. **mRNA (Indices `0` to `T-1`):** * Set to `alpha_flat` (e.g., `0.01`). (The C++ engine will add the OVR to this).
3. **nRNA (Indices `T` to `T+N-1`):** * If eligible: Set to `nrna_sparsity_alpha` (e.g., `0.9`).
* If dead (e.g., single-exon): Set to `0.0`.


4. **gDNA (Index `T+N`):** * Set to `1.0 + (gdna_prior_scale * gdna_init)`.
* *(Adding 1.0 ensures it safely acts as a flat uniform prior even if `gdna_init` is exactly 0).*



### Phase 3: Engine Modification (C++)

Modify `compute_ovr_prior_and_warm_start()` in `em_solver.cpp` so it respects the pre-computed nRNA and gDNA priors and restricts OVR solely to mature transcripts.

1. **Restrict OVR:** Wrap the OVR coverage-sharing logic in a condition so it *only* adds mass to mRNA components (`if i < n_t`).
2. **Preserve Purity:** Ensure the $\alpha$ values passed from Python for nRNA and gDNA are completely untouched by the OVR calculations.
3. **Keep Warm Starts Unchanged:** You can continue using your intronic evidence (`nrna_init`) to seed the *initial* $\theta$ weights (the starting line) for the EM, but ensure they do not alter the `prior` array (the finish line).

### Phase 4: Execution & Cleanup (C++)

1. **Run VBEM:** Execute the EM using the `"vbem"` mode. The digamma function $\psi(\alpha)$ will automatically enforce the sparsity toll for nRNA and the continuous density anchor for gDNA.
2. **Post-EM Truncation:** Because VBEM pushes sparse components asymptotically close to zero (e.g., `1e-8`) rather than absolute zero, add a final cleanup loop after the EM converges:
* `if (theta[i] < 1e-5) theta[i] = 0.0;`
* Re-normalize the remaining $\theta$ array to sum to 1.0.



This plan requires zero changes to your actual EM math or Empirical Bayes logic. It strictly enforces the biological priors at the boundaries between Python and C++.