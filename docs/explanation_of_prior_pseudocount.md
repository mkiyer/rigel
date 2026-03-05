

 Let me trace the exact code path to give you a precise answer.  

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/locus.py#1-1), lines 1 to 50

Searched for text `def build_locus_em_data` (`**/src/hulkrna/locus.py`), 1 result

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/locus.py#147-147), lines 147 to 350

Now let me read the EM initialization in `run_locus_em` where the OVR is constructed: 

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/estimator.py#660-660), lines 660 to 780

Now let me also read the `_em_step` to show the full M-step formula:  

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/estimator.py#120-120), lines 120 to 155

Now I have the full picture. Here's the precise walkthrough:

---

## How the priors actually work — a complete trace

You're right that there's an accumulation problem. Let me walk through the exact data flow:

### Step 1: Base prior construction ([locus.py, line ~340](src/hulkrna/locus.py#L340))

```python
prior = np.full(n_components, counter.em_prior, dtype=np.float64)
```

Every component (mRNA, nRNA, gDNA) starts with `prior[k] = em_pseudocount` (was 0.50, now 0.01). Then certain components get zeroed: gDNA if no unspliced fragments, nRNA if single-exon or zero nrna_init. This is the **flat Dirichlet pseudocount** — call it $\alpha_{\text{base}}$.

### Step 2: OVR is ADDED on top ([estimator.py, line ~726](src/hulkrna/estimator.py#L726))

```python
prior = prior + coverage_totals / n_ambiguous
```

This takes the coverage-weighted shares (how much "geometric attention" each component got from ambiguous fragments), normalizes them to sum to 1.0, and **adds** them to the existing prior. So the final prior per component is:

$$\alpha_k = \alpha_{\text{base}} + \frac{\text{coverage\_totals}_k}{N_{\text{ambiguous}}}$$

The OVR part sums to exactly 1.0 across all components. **It was never a replacement — it was always additive.**

### Step 3: MAP-EM M-step ([estimator.py, line ~148](src/hulkrna/estimator.py#L148))

```python
theta_new = unambig_totals + em_totals + prior
theta_new /= theta_new.sum()
```

So after each iteration, the unnormalized count for component $k$ is:

$$\tilde{\theta}_k = \text{unique}_k + \text{em}_k + \underbrace{\alpha_{\text{base}} + \text{OVR}_k}_{\text{prior}_k}$$

### Why this caused excessive FPs

Consider a locus with 41 transcripts (like FGFR2): that's $2 \times 41 + 1 = 83$ components. With $\alpha_{\text{base}} = 0.50$:

- **Total flat prior mass** = $83 \times 0.50 = 41.5$ virtual fragments spread across all components
- **OVR mass** = 1.0 total (negligible by comparison)
- For a zero-truth transcript with no unambig counts, its prior alone is $0.50 + \text{tiny OVR}$, which is enough to attract EM counts in the E-step and create a self-sustaining false positive

With $\alpha_{\text{base}} = 0.01$:

- **Total flat prior mass** = $83 \times 0.01 = 0.83$ virtual fragments total
- **OVR mass** = 1.0 (now actually *larger* than the flat prior!)
- The flat prior can no longer sustain false positives — only components that receive coverage-weighted OVR mass maintain meaningful prior support

**The OVR was always doing the right thing — distributing prior mass according to geometric evidence. But at $\alpha = 0.50$, the flat prior completely dominated it (41.5 vs 1.0), making the OVR essentially irrelevant.** At $\alpha = 0.01$, the OVR becomes the *primary* source of prior information, which is what you originally intended.

### The conceptual fix

Your original intuition — "replace the flat pseudocount with a coverage-weighted virtual read" — is actually what happens *in practice* when $\alpha_{\text{base}}$ is small enough. At $\alpha = 0.01$, the prior is effectively:

$$\alpha_k \approx \text{OVR}_k$$

with just a tiny $0.01$ floor to prevent numerical degeneracy. The Phase C results confirm this works well: MAP 0.01 beats MAP 0.50 on every region, every metric except dropout.

### Why VBEM still beats MAP 0.01 on FPs (51 vs 110)

VBEM uses $\psi(\alpha_k) - \psi(\sum \alpha)$ instead of $\log(\theta_k)$ for E-step weights. For small $\alpha_k$, $\psi(\alpha_k)$ is extremely negative (e.g., $\psi(0.01) \approx -100.6$), which exponentially suppresses low-evidence components far more aggressively than $\log(\theta_k)$. This is an inherent mathematical property of the digamma function — it acts as a built-in sparsity inducer that MAP-EM's log can't match.

### What this means for the path forward

The simplest clean solution would be to **make the flat $\alpha_{\text{base}}$ even smaller** (or zero it entirely) and rely purely on OVR + a small numerical floor. This would make the prior *entirely* coverage-weighted with no flat component, which is exactly your original intent. Shall I test $\alpha = 0.001$ or even $\alpha = 0$ (with a tiny epsilon floor) to see if we can close the remaining FP gap with MAP-EM?