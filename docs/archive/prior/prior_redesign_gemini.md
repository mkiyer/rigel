Your critique identifies a fundamental hazard in mixture-model optimization for transcript quantification. If a locus contains $K$ candidate RNA transcripts (where $K$ can be large and expression is naturally sparse) and only one gDNA component, treating the paired prior problem by injecting a flat or uniform prior count across all individual RNA transcripts introduces a severe structural distortion:

1. **Sparsity Destruction:** Injecting even a small constant prior mass (e.g., $\alpha_k = 0.5$ or $\alpha_k = \alpha_{RNA}/K$) to every individual transcript exerts a persistent "gravitational pull" that prevents unexpressed isoforms from converging to zero.
2. **Scale Asymmetry:** If $\alpha_{RNA}$ is split uniformly among 40 transcripts, each receives a negligible count (e.g., $0.05$). Individually, these tiny counts fail to provide the necessary counter-weight against a concentrated gDNA prior count (e.g., $\alpha_{gDNA} = 2.0$) on a single component, allowing gDNA to inappropriately siphon reads from low-expression transcripts.

To resolve this tradeoff, the prior must operate at the **locus group level** for RNA, allowing the EM optimization to determine the internal distribution among transcripts without prior-induced distortion. Below is an enumeration of potential design strategies, followed by the ideal structural solution.

---

### Enumeration of Potential Design Plans

#### Plan A: Flat Uniform Allocation (The Naive Approach)

Divide the total calibrated RNA count $\alpha_{RNA}$ uniformly ($\alpha_k = \alpha_{RNA} / K$) or weight it by transcript effective length, adding it directly to each component's prior vector.

* **Pros:** Requires zero changes to the native C++ M-step solver; fits the existing `prior + posterior` array accumulation loop.
* **Cons:** Destroys expression sparsity. It introduces a baseline expression floor for all unexpressed isoforms, generating significant false-positive transcript calls and breaking down under SQUAREM acceleration.

#### Plan B: Evidence-Gated Prior Activation (The Heuristic Approach)

Pre-scan the locus fragments or run a single unregularized E-step iteration to identify transcripts with non-zero read support. Apply the prior count $\alpha_k$ only to active transcripts.

* **Pros:** Preserves absolute zero states for transcripts with no overlapping fragments.
* **Cons:** Introduces a discontinuous threshold that causes optimization instability. It fails to handle transcripts that lack unique mappers but possess valid multi-mapping evidence, replicating the responsibilities of the EM loop itself.

#### Plan C: Group Dirichlet / Hierarchical Prior (The Recommended Structural Solution)

Define a joint prior over the *aggregate* RNA expression mass versus the gDNA component. Instead of penalizing individual transcripts, the prior constrains the total locus allocation:


$$p(\theta) \propto \theta_{gDNA}^{\alpha_{gDNA} - 1} \left( \sum_{k \in RNA} \theta_k \right)^{\alpha_{RNA} - 1}$$

* **Pros:** Mathematically rigorous and perfectly preserves sparsity. If a transcript has zero evidence ($N_k = 0$), its prior allocation automatically drops to zero. The aggregate RNA mass exerts the exact required pull to balance gDNA, while the internal partitioning among isoforms remains entirely data-driven.
* **Cons:** Requires a minor modification to the M-step normalization logic within the C++ kernel.

---

### The Ideal Solution: Group Dirichlet Prior (Plan C)

By modeling the RNA transcripts as a cohesive group, the prior acts on the simplex boundary between total RNA and total gDNA.

#### Mathematical Formulation

Let $\theta_g$ represent the abundance of the single gDNA component, and $\theta_k$ represent the abundance of transcript $k$ (for $k = 1, \dots, K$). They satisfy the constraint $\theta_g + \sum_{k=1}^K \theta_k = 1$.

During the EM M-step, let $N_g$ be the accumulated posterior fragment count for gDNA, and $N_k$ be the accumulated posterior count for transcript $k$. The total accumulated RNA count is $N_{RNA} = \sum_{k=1}^K N_k$.

Incorporating the Group Dirichlet prior, the MAP objective function to maximize is:


$$\mathcal{L}(\theta) = N_g \log \theta_g + \sum_{k=1}^K N_k \log \theta_k + (\alpha_{gDNA} - 1) \log \theta_g + (\alpha_{RNA} - 1) \log \left(\sum_{k=1}^K \theta_k\right) + \lambda \left(1 - \theta_g - \sum_{k=1}^K \theta_k\right)$$

Taking the partial derivatives and solving the system of equations yields a clean closed-form update:

1. **Total Locus Mass Denominator ($M$):**

$$M = N_g + N_{RNA} + \alpha_{gDNA} + \alpha_{RNA} - 2$$


2. **gDNA Component Abundance ($\theta_g$):**

$$\theta_g = \frac{N_g + \alpha_{gDNA} - 1}{M}$$


3. **Individual Transcript Abundance ($\theta_k$):**

$$\theta_k = \frac{N_k}{M} \left( 1 + \frac{\alpha_{RNA} - 1}{N_{RNA}} \right), \quad \text{given } N_{RNA} > 0$$



#### The Sparsity-Preserving "Dynamic Pseudocount"

This formulation can be re-interpreted as applying a **dynamic, evidence-weighted pseudocount** ($\Delta_k$) to each transcript:


$$\theta_k = \frac{N_k + \Delta_k}{M} \quad \text{where} \quad \Delta_k = N_k \cdot \frac{\alpha_{RNA} - 1}{N_{RNA}}$$

* If an isoform has no read support ($N_k = 0$), then $\Delta_k = 0$. The transcript receives **zero prior boost**, allowing it to hit a true zero and maintain sparsity.
* If an isoform captures the majority of the expression evidence, it dynamically assumes the corresponding share of the protective RNA prior mass to defend itself against gDNA siphoning.
* The total prior mass added across all transcripts simplifies perfectly to: $\sum_{k=1}^K \Delta_k = \alpha_{RNA} - 1$, exerting the exact locus-level gravitational pull required.

---

### Implementation Blueprint

To wire this into the `rigel` architecture, changes are localized to the prior projection phase and the native M-step solver.

#### 1. Python Handoff Update (`src/rigel/calibration/prior.py`)

Modify `PriorTable` to carry paired locus-level prior scalars rather than a per-component dense vector:

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    gdna_prior_numerator: np.ndarray  # float64[L_loci] (corresponds to alpha_gdna)
    rna_prior_numerator: np.ndarray   # float64[L_loci] (corresponds to alpha_rna)
    gdna_eff_len: np.ndarray          # float64[L_loci]
    enable_gdna: np.ndarray           # uint8[L_loci]

```

These numerators are derived from the regional calibration arrays (`mu_gdna` and `rna_lower`) integrated across the locus blocks.

#### 2. Native C++ M-Step Update (`src/rigel/native/em_solver.cpp`)

In the M-step loop of `batch_locus_em_partitioned`, implement the Group Dirichlet scaling factor. This avoids altering the parallelized, performance-critical E-step kernel:

```cpp
// Within the serial/M-step portion of the locus solver loop:
double n_gdna = em_totals[n_transcripts]; // gDNA component is the last slot
double n_rna_total = 0.0;

for (int t = 0; t < n_transcripts; ++t) {
    n_rna_total += em_totals[t];
}

// Fetch locus-level paired prior values from the updated PriorTable
double alpha_gdna = locus_prior_gdna[locus_idx];
double alpha_rna  = locus_prior_rna[locus_idx];

// Target denominator
double M = n_gdna + n_rna_total + alpha_gdna + alpha_rna - 2.0;
if (M <= 0.0) { M = EM_LOG_EPSILON; }

// Update gDNA abundance
theta_new[n_transcripts] = std::max(n_gdna + alpha_gdna - 1.0, 0.0) / M;

// Update RNA transcript abundances using the dynamic scaling factor
if (n_rna_total > 0.0) {
    double rna_scale_factor = 1.0 + (alpha_rna - 1.0) / n_rna_total;
    for (int t = 0; t < n_transcripts; ++t) {
        theta_new[t] = std::max(em_totals[t] * rna_scale_factor, 0.0) / M;
    }
} else {
    // Boundary condition handling when no RNA fragments map to the locus
    handle_zero_rna_locus(theta_new, alpha_rna, n_transcripts, M);
}

```

#### 3. Handling Boundary Conditions ($N_{RNA} = 0$)

If a locus contains fragments mapping to gDNA but absolutely zero fragments mapping to any RNA transcript ($N_{RNA} = 0$), the scaling factor denominator encounters a division-by-zero. To remain numerically stable and well-behaved, `handle_zero_rna_locus` can distribute the total group prior mass ($\alpha_{RNA} - 1$) proportional to transcript effective lengths, or uniformly:


$$\theta_k = \frac{\alpha_{RNA} - 1}{K \cdot M}$$


Because there is no active read evidence to distort, distributing this minute baseline under a true-zero observation state will not trigger false-positive transcript calls; the entire group remains suppressed relative to the active gDNA component.