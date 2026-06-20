# TODO list for entire package

## Tiny nodes (count or size)

Nodes will zero (or low) counts kill propagation

Nodes (regions) with tiny length can also kill propagation. Very small regions with length much less than fragment length are not large enough to hold contained fragments. The regions would therefore have zero counts (correctly).

The problem is that these nodes will kill any message propagation.

A caveat with message propagation is the node size. Boundaries have a fixed size and depend upon fragment length. Short fragments will be more likely contained and less likely boundary crossing. Long fragments are more likely to cross boundaries and less likely to be contained.

Our method is currently not robust to a variety of fragment lengths. It is especially vulnerable to long fragment lengths. If a region is much smaller than the fragment length, we end up with zero (or very low) counts. The region kills message propagation.


We may need to improve the solution for message propagation. Instead of killing the message, nodes with 0 counts or tiny size should act as a "relay" or "router". They have no information to contribute to support or refute the propagating message. They shouldn't kill the message and they shouldn't change the message.

The "belief propagation" idea is critically important. As a belief (message) propagates, each node has the opportunity to dampen the message, amplify the message, or relay the message. Nodes with zero information should act as relays.



### Answer from Gemini 6/18/2026

This is a fascinating architecture you're building! Combining a bipartite DAG with a bidirectional message-passing sweep to partition RNA and gDNA is an elegant way to handle the ambiguities of transcriptomic alignment.

You hit the nail on the head regarding the problem. This is a classic **structural zero** caused by a mismatch between your feature size (node length) and your assay's physical constraints (fragment length). If a node is 40bp and your fragments are 200bp, it is physically impossible to have a fragment *contained* within that node.

If your message-passing math multiplies incoming signals by the node's internal contained count (or a probability directly tied to it), a zero there acts as a total circuit breaker.

Here are a few ways to fix this, ranging from structural graph modifications to algorithmic adjustments.

---

## 1. The Algorithmic Fix: Decouple "Contained Mass" from "Node Abundance"

The cleanest fix doesn't change your graph at all; it changes how the node calculates its transmitting power. Currently, your node updates are likely treating "0 contained fragments" as "0 expression." Instead, the node's true hidden variable—let's call it Abundance ($\theta_i$)—should be driven by both its internal counts *and* the flux of the edges passing through it.

* **Effective Length Conditioning:** For any node $i$, the effective length for a *contained* fragment is $\tilde{L}_i = \max(0, L_i - \mu + 1)$, where $\mu$ is the fragment length.
* **The Adjustment:** If $\tilde{L}_i = 0$, you must omit or bypass the node's internal containment likelihood term during the sweep. Instead, the node acts strictly as a **pure router (or transit hub)**. Its abundance is estimated entirely by conservation of flow from the incoming edges:

$$\theta_i = \sum_{j \in \text{In}(i)} \text{Message}_{j \to i}$$

During the sweep, instead of multiplying the incoming message by the internal count, the node normalizes and distributes the incoming edge mass to the outgoing edges based on transition probabilities, completely ignoring the fact that its internal count is zero.

---

## 2. The Structural Fix: Topological Contraction (With Conditions)

You mentioned merging nodes with surrounding nodes but worried about creating too many paths. Your intuition is correct: blind merging in a complex DAG can cause a combinatorial explosion of alternative paths, wrecking your graph's resolution.

However, you *can* safely collapse nodes if you apply **strict topological rules**:

* **The "Linear Bridge" Exception:** If a tiny node has an $\text{in-degree} = 1$ and an $\text{out-degree} = 1$, it is a simple linear bridge. You can safely collapse this node into its upstream and downstream neighbors (or absorb it entirely into a single continuous edge) without changing the global path topology.
* **The Junction Rule:** If the tiny node has an $\text{in-degree} > 1$ or an $\text{out-degree} > 1$ (e.g., it's an alternative splice site choice or a transcription start site), **do not merge it**. Merging a junction node alters the logic of the DAG and creates phantom paths that don't exist biologically. For these nodes, you *must* use the algorithmic fix in Option 1.

---

## 3. The Pseudo-Mass / Prior Distribution Approach

If modifying the core sweep logic feels too invasive, you can use a regularized Bayesian prior to prevent the math from crashing into a hard zero.

* **Fractional Pseudo-counts:** Instead of initializing a node with 0 mass, assign it a tiny fractional pseudo-count based on its capacity.
* **Edge-to-Node Imputation:** Before you start the main bidirectional sweep, run a quick local preprocessing step. Look at the fragments sitting on the *edges* connected to the tiny node. If Edge A $\to$ Tiny Node has 50 fragments, and Tiny Node $\to$ Edge B has 50 fragments, you know structurally that roughly 50 fragments traversed this node. You can artificially seed the node with an "implied transit mass" scaled by its length relative to the fragment size, preventing the zero-multiplication problem entirely.

---

## Summary Recommendation

If your pipeline can handle it, **Option 1 (Flux-based routing for $\tilde{L}=0$)** is the most mathematically rigorous because it honors the true physical constraints of RNA-seq without altering your graph topology. If you want a quick structural fix, implement **Option 2** but restricted *only* to nodes where $\text{in-degree} = 1$ and $\text{out-degree} = 1$.

How are your transition probabilities on the edges currently calculated—are they derived entirely from the edge-spanning fragment counts, or do they rely on the relative abundances of the connected nodes?


## Simplex solver normalization

## Calibrating lattice grid dimension

don't allow solver to have ridiculously high precision
- recommended grid_size * grid_size as cap, justification?
- recommended capping at max strand precision


## Linear graph -> DAG

docs/accumulator/dag.md

## FL distribution in calibration

may help discriminate gdna from rna during calibration
requires per-node FL distributions
deferred, potentially low ROI