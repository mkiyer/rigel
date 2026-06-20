This is a phenomenal, clear-eyed diagnosis. You have perfectly isolated the core architectural vulnerability: the simplex solver operates as a precision-weighted democracy, and right now, certain terms are executing a bloodless coup by inflating their precision $\tau$ to astronomical levels.

Your diagnosis of the **AMBIG node vulnerability** (where the strand likelihood is flat along the $f_g$ axis) explains exactly why these nodes are the first to break when priors go rogue.

However, while your proposed fixes ($\rho_{src}$ querying and self-solvable global pooling) are excellent first-order band-aids that will patch the current leak, there is a **deeper structural flaw** hidden in the mathematics of your current formulation. If left unaddressed, this flaw will likely manifest as a new failure mode under different capture designs or depth regimes.

Here is a deep critique of the proposal, the hidden mathematical trap, and architectural alternatives for your working session.

---

## 1. The Hidden Structural Trap: The Quadratic vs. Linear Scaling Mismatch

In section 4.2, the document states: *"the $(M/E)^2$ jacobian (correct in itself) then magnifies the under-estimate..."* While the Jacobian is mathematically correct for a rigid change-of-variables, its interaction with the optimization objective creates a fundamental dimensional mismatch between your data terms and your propagation terms.

Consider how the terms scale relative to a massive increase in sequencing depth or capture enrichment (high mass $M$):

* **Strand Likelihood Information:** The Fisher Information $I(t)$ scales **linearly** with count $N$ (and therefore linearly with mass $M$).

$$I(t) \propto M$$


* **Imputation Message Precision:** The fraction precision $\tau_f$ is calculated via the Jacobian:

$$\tau_f = \tau_\rho \cdot \left(\frac{M}{E}\right)^2$$



If biological dispersion $\sigma^2_{bio}$ dominates $\tau_\rho$, then $\tau_\rho \approx 1/\sigma^2_{bio}$. This means $\tau_f$ scales **quadratically** with mass $M$.

$$\tau_f \propto M^2$$



### Why this is catastrophic

Biological spatial variance ($\sigma^2_{bio}$) represents how much the true template density fluctuates between adjacent genomic regions. **This fluctuation is an intrinsic property of biology, completely independent of how deep you sequence.** By multiplying $1/\sigma^2_{bio}$ by $(M/E)^2$, the model assumes that if an exon is highly enriched ($M$ is huge), the biological coupling between it and its neighbor magically becomes exponentially more rigid.

Because $M^2$ always outgrows $M$, **any highly enriched node will naturally allow its imputation messages to steamroll the direct strand evidence**, even if you fix the circularity bug. This is the root cause of the over-confidence.

---

## 2. Critique of Your Proposed Fixes

### 4.1 Global Prior Fix (Honest Mean & Density-Dependent Width)

* **The Good:** Switching to `self_solvable` regions correctly integrates the capture enrichment signal into the baseline population mean. Using the observed total density $\rho_{obs} = M/E$ completely breaks the circular feedback loop.
* **The Risk:** $\rho_{obs}$ can be highly noisy or skewed by local mapping artifacts (e.g., multi-mapping reads, ultra-short effective lengths). If $\rho_{obs}$ spikes artificially at a node, it will query a massive variance, causing the prior to instantly abandon the node. This is benign over-confidence turning into harmless under-confidence, but worth monitoring.

### 4.2 RNA Message Fix (Querying at $\rho_{src}$)

* **The Good:** Querying $\sigma^2_{bio}$ at $\rho_{src}$ instead of $\rho_{dst}$ stops the localized "sinkhole" effect where a slightly crushed density self-proves its own absolute certainty.
* **The Risk (Cascading Collapse):** While this fixes a single-node circularity loop, it does not fix a **directional chain collapse**. If an upstream `src` node happens to be crushed to zero density by a upstream neighbor, its $\rho_{src}$ will be $\sim 0$. When the downstream node queries $\sigma^2_{bio}(\rho_{src} \to 0)$, it will receive a tiny biological variance, generating an enormous $\tau$ that forces the downstream node to collapse to zero as well. The collapse will simply propagate down the graph like dominoes.

---

## 3. How to Radically Fix the Problem

To build a structurally bulletproof model, you need to prevent $\tau$ from escaping to infinity and ensure that technical/biological information scales on the same dimension.

### Alternative A: Linearize the Message Constraint (The Pseudo-Count Framework)

Instead of injecting neighbors' hints as rigid quadratic penalties on fractions ($-\frac{1}{2}\tau_f(f-\mu)^2$), translate the imputation message into **pseudo-counts** injected directly into the count-based likelihood.

If a neighbor says the imputed density should be $\mu_\rho$, translate this into an equivalent pseudo-count:


$$C_{pseudo} = \mu_\rho \cdot E_{dst}$$


Give this pseudo-count a weight/effective depth based entirely on the biological dispersion:


$$N_{pseudo} = \frac{1}{\text{dispersion}^2}$$


This forces the imputation message to compete on an even playing field with the strand reads. If a node has $16,000$ real physical reads, a weak biological hint worth $\sim 50$ pseudo-counts will naturally be outvoted. If a node has 2 reads, the message will gently guide it. Information will scale linearly ($M$ vs $M$) instead of quadratically ($M$ vs $M^2$).

### Alternative B: Move Optimization to Log-Ratio Space

The simplex boundary ($f_+ + f_- + f_g = 1$) is a dangerous place for Gaussian optimization because variances are absolute, but fractions are relative.

If you map your fractions to **Additive Log-Ratios (ALR)** using gDNA as the baseline:


$$y_+ = \log\left(\frac{f_+}{f_g}\right), \quad y_- = \log\left(\frac{f_-}{f_g}\right)$$


The boundaries $[0, 1]$ open up to $[-\infty, +\infty]$. In log-ratio space, biological variance $\sigma^2_{bio}$ acts as a constant dispersion factor across all expression levels (handling the $10^4$-fold dynamic range of RNA natively). A crushed RNA density simply means $y \to -\infty$ with a modest, stable variance, eliminating the hyper-inflated precision spikes.

### Alternative C: The Defensive Floor (The Immediate Practical Implementation)

If you must keep the current solver architecture for the upcoming milestone, you must introduce a non-zero **variance floor** to structurally guarantee that no message can ever achieve an infinite veto power.

1. Calculate the global dynamic range of your RNA library (e.g., the 90th percentile of variance among your self-solvable nodes).
2. Enforce a hard ceiling on precision:

$$\tau_{\rho} = \frac{1}{\max\left(\sigma^2_{bio}(\rho_{src}), \, \sigma^2_{floor}\right) + V_p}$$


3. This guarantees that $\tau_f$ can never climb to $1e11$, preserving the local strand likelihood's ability to fight back and win when it has thousands of reads.

---

## 4. Recommendations for Your Working Session

1. **Implement the floor immediately:** Add a defensive maximum precision cap to `honest_message_tau` to break the failure mode while you experiment.
2. **Test `ρ_src` vs `ρ_obs` for RNA:** Use your `dissect_regions.py` harness to test if querying at $\rho_{src}$ triggers the cascading chain collapse described in Section 2. If it does, using $\rho_{obs}$ (total observed mass density) is much safer because it anchors the variance query to raw data rather than moving graph states.
3. **Unify the logic:** Do not stage this purely for RNA. Build a unified `compute_honest_precision()` helper for the Global Prior, RNA, and gDNA messages. If a structural design rule is true for RNA, it is true for gDNA—the only reason gDNA isn't breaking right now is because capture enrichment targets exons, masking the problem.

Given that switching to pseudo-counts or log-ratio space would require reworking the solver's underlying grid/lattice structure, would you like to sketch out the exact mathematical formulation for a safe, non-circular biological variance floor ($\sigma^2_{floor}$) that can be dropped directly into your existing P-spline model?



=================

Switching from continuous Gaussian penalties in fraction space to **pseudo-counts in log-fraction space** (effectively treating messages as a local Dirichlet/Multinomial update) is an exceptionally elegant architectural pivot.

By communicating in count space, you treat information from adjacent nodes exactly like real physical sequencing fragments. This unifies your objective function into a single coherent dimension: **log-likelihood of discrete counts**.

Crucially, it resolves the dangerous quadratic scaling explosion ($M^2$) introduced by the fraction Jacobian. If the biological variance scales with the square of the density (a constant coefficient of variation $CV_{bio}$), the total pseudo-counts a message can deliver hits a hard structural ceiling:


$$N_{pseudo} = \frac{\rho^2}{\sigma^2_{bio}(\rho)} = \frac{1}{CV_{bio}^2}$$


Even if an exon has millions of reads, its biological neighbor can never input more than $N_{pseudo}$ (e.g., $\sim 100$ counts) of structural confidence into the simplex local grid.

Below is the blueprint for refactoring your codebase to implement this count-space message framework.

---

### 1. The Core Transformation: Modifying `simplex_sweep.py` & `simplex.py`

Currently, `_local_loglik` evaluates your messages using a Gaussian penalty: `−½ · τ · (f − μ)²`. You need to replace this with a Multinomial log-likelihood term: $\sum \alpha_c \log(f_c)$.

Since your simplex lattice points (`f_pos_g`, `f_neg_g`, `f_g_g`) are fixed grids generated by `_simplex_lattice`, you can pre-compute their logarithms to maintain high performance and prevent $\log(0)$ boundary exceptions.

#### Changes in `simplex_sweep.py`:

Modify the `_local_loglik` signature and implementation to swap out the fraction/precision parameters for pseudo-count ($\alpha$) parameters.

```python
# In rigel/calibration/simplex_sweep.py

def _local_loglik(
    u_pos, u_neg, spliced_pos, spliced_neg, allow_pos, allow_neg, kappa,
    od_g, od_r, lattice,
    strand_obs=None, global_mu=0.0, global_tau=0.0,
    # REPLACE THE OLD IMPUTATION ARGUMENTS WITH DIRECT PSEUDO-COUNTS (ALPHAS)
    alpha_pos=None, alpha_neg=None, alpha_g=None
):
    """
    Computes the local log-likelihood over the 2-simplex lattice.
    alpha_pos, alpha_neg, alpha_g should have shape (nodes, 1) to broadcast 
    against the lattice of shape (1, P).
    """
    f_pos_g, f_neg_g, f_g_g = lattice
    
    # 1. Compute baseline intrinsic strand likelihood
    psi = _mixture_strand_loglik(
        u_pos, u_pos + u_neg, f_g_g, f_pos_g, f_neg_g, kappa, od_g, od_r
    )
    
    # Pre-calculate protected log-lattice positions to prevent log(0)
    # The grid points are static, so this is safe and incredibly fast
    log_f_pos = np.log(np.maximum(f_pos_g, 1e-12))
    log_f_neg = np.log(np.maximum(f_neg_g, 1e-12))
    log_f_g   = np.log(np.maximum(f_g_g, 1e-12))

    # 2. Inject Imputation Messages as Multinomial Pseudo-Counts
    # This replaces the old continuous quadratic penalties entirely!
    if alpha_pos is not None:
        psi += alpha_pos * log_f_pos
    if alpha_neg is not None:
        psi += alpha_neg * log_f_neg
    if alpha_g is not None:
        psi += alpha_g * log_f_g

    # 3. Keep your global prior / boundary constraints intact
    if global_tau > 0.0:
        psi -= 0.5 * global_tau * (f_g_g - global_mu) ** 2
        
    # [Rest of your existing spliced lower-bound / masking logic remains unchanged]
    ...
    return psi

```

---

### 2. Updating Message Generation in `bp_solver.py`

This is where the massive simplification happens. In your current pipeline, you take the neighbor's imputed density, compute the target fraction, query the variance model, and then perform a dangerous Jacobian correction $(M/E)^2$ to get a fraction-space precision.

**You can throw away the Jacobian correction completely.** Instead, calculate how many discrete pseudo-counts the source node's density scale is authorized to speak for.

#### Changes in `bp_solver.py`:

Locate where messages are prepped for the next node's simplex solve inside your directional sweep loops:

```python
# In rigel/calibration/bp_solver.py (Inside your pass loops)

# 1. Grab the source node's converged density components from the current snapshot
rho_src_pos = density_pn[src_nodes] * f_pos[src_nodes]
rho_src_neg = density_pn[src_nodes] * f_neg[src_nodes]
rho_src_g   = density_pn[src_nodes] * f_g[src_nodes]
rho_src_total = density_pn[src_nodes]

# 2. Query the biological variance at the source density scale
# This isolates the intrinsic biological predictability floor
var_bio = var_mean_model.predict(rho_src_total)
var_bio = np.maximum(var_bio, 1e-8) # Defensive floor to prevent division by zero

# 3. Compute the Total Pseudo-Count Capacity (N_pseudo) allowed by biology
# If var_bio track overdispersion (phi * rho^2), then n_pseudo simplifies cleanly to 1/phi
n_pseudo = (rho_src_total ** 2) / var_bio

# 4. Partition the total pseudo-counts into component alpha channels
# Reshape to (nodes, 1) for proper downstream simplex lattice broadcasting
alpha_pos = (n_pseudo * f_pos[src_nodes]).reshape(-1, 1)
alpha_neg = (n_pseudo * f_neg[src_nodes]).reshape(-1, 1)
alpha_g   = (n_pseudo * f_g[src_nodes]).reshape(-1, 1)

# 5. Pass these alphas straight into your updated _solve_nodes routine
# No fractions, no matrix inversion, no Jacobians.

```

---

### 3. Adapting the Variance Model Fitting (`variance_model.py`)

Your current `MonotoneVarMean.fit_offset` method is already exceptionally well-designed for this change. It already decomposes the total observed cross-node density variance into a biological component and a Poisson sampling floor:


$$\mathbb{E}[R^2] = \sigma^2_{bio}(\mu) + V_p$$

Because it intentionally subtracts out $V_p$ during optimization, `predict(mean)` cleanly isolates the **pure biological excess uncertainty** ($\sigma^2_{bio}$). This is exactly the value needed to establish your pseudo-count capacity in step 2 above.

No changes are strictly required for the internal P-spline math of `variance_model.py`. However, you should update the documentation and the docstrings to reflect that it now exports biological variance to scale pseudo-counts rather than calculating fraction precisions:

```python
# In rigel/calibration/variance_model.py

# Update the main description to reflect the new paradigm:
"""
...
bp_solver fits this each sweep pass on cross-node density residuals.
The model isolates biological excess variance σ²_bio. 
Instead of calculating fraction-space precisions, bp_solver uses this 
to scale the effective sample size of messages: N_pseudo = mean² / σ²_bio.
"""

```

---

### Summary of Benefits to Your Working Session

1. **Information Symmetry:** Real reads entering the solver provide information scaling linearly with read count. Imputation messages entering via the Dirichlet terms now *also* scale linearly with count, up to the saturation point $N_{pseudo}$. The two terms can now compete fairly.
2. **Elimination of Boundary Spikes:** In your old approach, if a fraction approached 0, the quadratic parameterization could force the precision to spike, artificially pinning values against the simplex walls. In log-fraction space, a zero-count message ($\alpha_c = 0$) simply means $0 \cdot \log(f_c) = 0$, dropping out of the equation natively without exerting any toxic pulling force.
3. **Simpler Code:** You completely sidestep the maintenance and debugging of derivative transformations and coordinate transitions. A message is now just an array of pseudo-counts.


=======================

This hits on the exact mathematical beauty of moving to count space.

The short answer is: **The recipient node handles the effective length mismatch automatically because its own local evidence (the strand likelihood) scales naturally with its own effective length.** By *not* normalizing the incoming pseudo-counts by the destination's effective length, you are allowing the destination node to perform a perfect, text-book Bayesian update where the **sample size** of local data naturally determines how much weight to give the neighbor's hint.

Here is exactly how the destination node integrates these un-normalized counts, why the units line up perfectly, and why this natively solves your scaling crisis.

---

### 1. The Unified Objective Function

When the destination node runs its local simplex solve, it maximizes its log-posterior lattice $\log \psi(f)$. Let's look at the mathematical units of the terms inside that objective function under the count-space model:

$$\log \psi(f) = \underbrace{L_{\text{strand}}(f, U_{\text{dst}})}_{\text{Scaled by Real Counts } N_{\text{dst}}} + \underbrace{\sum_{c} \alpha_{c} \log(f_c)}_{\text{Scaled by Virtual Counts } n_{\text{pseudo}}}$$

Both terms are now expressed in the exact same currency: **log-likelihood of discrete fragment counts**.

* **The Local Strand Term** ($L_{\text{strand}}$) is derived from the real, physical fragments captured at the destination node ($U_{\text{pos}}, U_{\text{neg}}$). Its absolute magnitude and curvature scale **linearly** with the destination node's total physical count $N_{\text{dst}}$.
* **The Message Term** ($\sum \alpha_c \log f_c$) is driven by the virtual pseudo-counts ($\alpha_c = n_{\text{pseudo}} \cdot f_{c,\text{src}}$) sent by the neighbor. Its absolute magnitude scales **linearly** with $n_{\text{pseudo}}$.

---

### 2. How Effective Length Automatically Balances the Vote

Let's trace your exact example. Suppose the source node tells the destination node that the biological ratio should be $f_{\text{src}} = [0.4, 0.4, 0.2]$ (Sense, Antisense, gDNA), and based on biological variance, it is authorized to send a total of $n_{\text{pseudo}} = 100$ virtual counts.

This means the destination node receives a message array: $\alpha = [40, 40, 20]$.

Now let's look at how two different destination nodes (one huge, one tiny) with the exact same true underlying density ($\rho = 0.1 \text{ fragments/bp}$) integrate this identical message:

#### Regime A: The Destination Node is HUGE ($E_{\text{dst}} = 10\text{kb}$)

* **Real Counts:** Because it is long, it catches a ton of fragments: $N_{\text{dst}} = 0.1 \times 10,000 = 1,000$ real physical counts.
* **The Integration:** Inside the solver, the local strand term has a commanding weight of $\sim 1,000$ counts, while the neighbor's incoming message only provides $100$ virtual counts.
* **The Outcome:** The local data holds a **91% majority vote** ($1000 / 1100$). The local data easily outvotes the neighbor's message. This is statistically correct: a node with 1,000 real counts has massive local precision ($\text{CV} \approx 3\%$) and does not need to be bossed around by a neighbor.

#### Regime B: The Destination Node is TINY ($E_{\text{dst}} = 1\text{kb}$)

* **Real Counts:** Because it is short, it catches very few fragments: $N_{\text{dst}} = 0.1 \times 1,000 = 100$ real physical counts.
* **The Integration:** Inside the solver, the local strand term only has a weight of $\sim 100$ counts. The neighbor's incoming message *also* provides $100$ virtual counts.
* **The Outcome:** The local data and the neighbor's message have an **equal 50/50 vote** ($100 / 200$). The neighbor's message acts as a powerful guide, smoothing out the local assignment. This is *also* statistically correct: a node with only 100 real counts is highly vulnerable to Poisson sampling noise ($\text{CV} \approx 10\%$). It *should* heavily rely on its spatial neighbors to stabilize its partition.

---

### 3. Why This Extinguishes the "Over-Confidence Coup"

In your original fraction-space design, the message precision $\tau_f$ was multiplied by the Jacobian $(M/E)^2$. If a destination node had a massive mass $M$, the penalty skyrocketed to $1e11$, meaning the neighbor message forced a **100% tyrannical veto** over local data, regardless of how much local evidence existed.

By communicating a fixed pool of pseudo-counts ($n_{\text{pseudo}}$), you place a hard, structural ceiling on how much influence a neighbor can ever wield.

If your variance model learns that the biological coefficient of variation between adjacent nodes plateaus at $CV_{\text{bio}} = 0.10$ (10% biological fluctuation), then:


$$n_{\text{pseudo}} = \frac{1}{CV_{\text{bio}}^2} = 100 \text{ counts}$$

No matter how deeply you sequence the library, and no matter how massive the effective length of the source or destination nodes becomes, **that neighbor can never cast more than 100 virtual votes in the destination node's local election.** If the destination node has 5,000 real reads, it will effortlessly override a wrong neighbor. If the destination node has 2 reads, the neighbor will gracefully protect it from dropping to 0.

### Summary

The recipient node doesn't need to normalize the pseudo-counts by effective length because **counts are the natural currency of evidence magnitude.** A larger effective length yields more real physical counts, which automatically dilutes the relative power of the incoming virtual pseudo-counts. The graph balances itself natively.


================

