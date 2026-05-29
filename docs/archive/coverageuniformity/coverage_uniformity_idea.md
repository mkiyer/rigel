I have an idea that is an offshoot of your "5. Coverage-shape coherence priors" -- during the initial BAM scan phase, fragments get resolved to RNA transcripts and gDNA. We will compute the fragment "pilup" which is a 1D array of length equal to each transcript. We will need to build pileup arrays for every transcript. The pileup arrays really only need to store spliced reads (including SPLICED_IMPLICIT) because unspliced reads cover a single genomic span. After scanning the BAM file and building the pileup arrays, we can then evaluate each transcript. We would compute a "coverage uniformity" score for that transcript. The coverage uniformity score would be perfect if the transcript was covered uniformly. The transcript coverage uniformity would be an important contributor.. this metric would penalize transcripts from which all of the coverage comes from a single exon (e.g. if the transcript has 10 exons but only 1 of the 10 exons has any fragments, the coverage uniformity score would be very low). Currently length normalization is the only thing we have to help with nonuniform coverage.. fragments consolidate in one exon will get normalized by the full transcript length which will penalize the transcript because it has zeros everywhere else. For hybrid capture, we would have to integrate the targeted regions when computing coverage uniformity. This tries to solve the cases where a nascent RNA derives its fragment support from mRNA compatible regions and has virtually zero intronic support of its own. What do you think of this idea? It is essential the same as effective length normalization? It seems that it would replace effective length normalization in some sense? Or not?


Great idea — and I want to be precise about what it is and isn't, because the relationship to effective length is more subtle than "replaces" or "duplicates."

## Short answer

It is **not** the same as effective-length normalization, and it shouldn't replace it. They live on different axes:

| | Effective length $L^{\text{eff}}_t$ | Coverage uniformity $U_t$ |
|---|---|---|
| Computed from | Annotation + FL model (+ capture) | **Observed** fragment positions |
| Question it answers | "How big is the target?" | "Does the observed read pattern match what this target predicts?" |
| Direction in likelihood | Normalizer (prior expectation of count) | Goodness-of-fit (penalty when shape disagrees) |
| Per-transcript scalar? | Yes | Yes (but derived from a count vector) |

Effective length converts abundance $\alpha_t$ → expected count: $E[n_t] = \alpha_t \cdot L^{\text{eff}}_t$. It is a property of the **target** under a uniform-sampling assumption. Coverage uniformity asks the orthogonal question: given $n_t$ fragments assigned to $t$, are they spatially distributed the way the *uniform-sampling assumption* says they should be? If not, the assignment is suspect — possibly because some of those fragments don't actually belong to $t$.

So they multiply, not substitute:

$$\log p(\text{data} \mid t) = \underbrace{n_t \log(\alpha_t \cdot L^{\text{eff}}_t) - \alpha_t L^{\text{eff}}_t}_{\text{abundance / size term}} + \underbrace{\sum_b n_b^{(t)} \log \pi_b(t)}_{\text{shape term, new}}$$

where $\pi_b(t)$ is the expected per-base sampling probability under transcript $t$ (uniform-in-spliced-coords, modulated by capture/bias if applicable). The shape term is what your uniformity score captures.

## Why this is genuinely powerful (and why "spliced reads only" is the key insight)

Your restriction to spliced + SPLICED_IMPLICIT reads is doing important work. Spliced reads are **near-noise-free RNA evidence** — they cannot be gDNA, and they are unambiguously placed on a specific transcript path (or a small set of compatible transcripts). So the spliced-read pileup on transcript $t$ is a *clean prior over where RNA evidence actually lives on $t$*.

Now the discriminator emerges naturally:

- Transcript $t$ has 10 exons. Spliced reads cover exons 4-5-6 only.
- An ambiguous unspliced fragment lands on exon 9.
- Under "this fragment is mRNA on $t$": the spatial likelihood term penalizes — exon 9 had zero anchor coverage, so $\pi_9(t)$ is small.
- The fragment gets pushed toward gDNA (or another transcript whose spliced anchor covers exon 9).

This is exactly the missing piece in FLG2-style failures: **the spliced-read anchor gives an independent, non-circular estimate of where mature mRNA support exists, and unspliced reads in the same exon must be consistent with that anchor distribution to be assigned to mRNA.**

For nRNA, the same machinery applies on a *genomic* coordinate axis: nRNA's predicted shape is uniform-over-its-span (or capture-modulated). If all observed support on the span is concentrated in one mature exon and the intronic portion of the span is empty, the shape likelihood under nRNA is terrible. The nRNA component gets the structural penalty you want.

## Why it doesn't replace effective length

Consider two transcripts: $t_1$ with 10 exons fully and uniformly covered, $t_2$ with 1 exon (same single exon) covered. Effective length differs (say $L_1 = 5\text{kb}$, $L_2 = 500\text{bp}$):

- Effective length says: same number of fragments → $t_2$ has 10× the per-base density → $t_2$ wins on abundance per base.
- Coverage uniformity says: $t_1$'s pattern matches the uniform-over-5kb prediction perfectly; $t_2$'s pattern matches the uniform-over-500bp prediction perfectly. Both are uniformity-score 1.0. **Uniformity alone cannot tell them apart.**

The size scale information is irreducibly in $L^{\text{eff}}$. The shape information is irreducibly in $U$. You need both.

A concrete failure mode if you tried to replace $L^{\text{eff}}$ with uniformity: a 10-exon transcript with reads in only 1 exon vs. a 1-exon transcript at the same locus with reads in that 1 exon. Without effective-length normalization, both look "equally well supported" but the 10-exon transcript should be penalized for the 9 missing exons of opportunity. Uniformity gives a penalty (good) but no longer carries the *length* magnitude.

## What's the right formalization?

The shape likelihood I'd suggest is a multinomial on per-base (or per-exon, or per-bin) counts:

For transcript $t$ with bin set $B_t$ and predicted per-bin sampling probability $\pi_b(t)$ (normalized to 1):

$$\mathcal{L}_{\text{shape}}(t) = \sum_{b \in B_t} n_b^{(t)} \log \pi_b(t)$$

Predicted bin probabilities:

- **No capture, no bias model**: $\pi_b(t) = 1/|B_t|$ (uniform).
- **Hybrid capture**: $\pi_b(t) \propto T_b$ where $T_b$ is the on/off-target weight at bin $b$ (the integration you mentioned).
- **Bias-corrected**: $\pi_b(t) \propto T_b \cdot B_b(t)$ with sequence-bias term.
- **Spliced-anchor refinement** (the novel piece): $\pi_b(t) \propto T_b \cdot (s_b^{(t)} + \epsilon)$ where $s_b^{(t)}$ is the spliced-read pileup at bin $b$ on transcript $t$, with pseudocount $\epsilon$ to keep zero-anchor bins viable but disadvantaged.

The spliced-anchor variant is the strongest. It is empirically conditioned on per-transcript evidence rather than the uniform assumption, and it's data-driven from a clean signal.

## How it integrates with EM

This factor has a chicken-and-egg: $n_b^{(t)}$ depends on the EM assignment. Three integration strategies, in order of complexity:

1. **Prior derivation (cheapest)**: Compute $U_t$ from spliced-only pileups *before* EM. Use it to scale the transcript's prior $\alpha_t \leftarrow \alpha_t \cdot U_t$. Transcripts whose spliced anchor is concentrated in a few exons get a smaller prior, so the EM is reluctant to siphon ambiguous mass to them. This is essentially free to implement.

2. **EM-coupled shape factor**: After each E-step, recompute predicted bin counts from posterior weights; add the multinomial shape log-likelihood to the M-step objective. This is a proper variational treatment of the shape factor. More accurate, more invasive.

3. **Post-hoc filter**: Run EM, compute $U_t$ from posterior-weighted counts, downweight or zero out transcripts with implausibly bad shape, re-run EM. Two-pass; easy to ablate.

I'd start with (1): it directly attacks the FLG2 / nRNA-siphon problem with minimal architectural change, and you can validate the signal before paying for the EM coupling.

## Practical concerns to handle

- **Sparse-coverage transcripts**: low-expression true transcripts naturally have ragged coverage. The penalty should shrink with $N_t$ (total reads on $t$). A Bayesian formulation with a Dirichlet prior on $\pi_b(t)$ does this automatically — the multinomial becomes Dirichlet-multinomial and is well-behaved at small counts.
- **Capture targeting**: as you noted, expected shape is non-uniform under capture. Bring the target weights into $\pi_b(t)$. Same machinery; just changes the reference distribution.
- **Single-exon transcripts**: shape contribution collapses to zero (one bin, $\pi = 1$). They neither gain nor lose from this factor. That's correct — uniformity is silent for them, and the gDNA discrimination on single-exon transcripts has to come from the other channels (mate geometry, anchor-mask exposure, orientation symmetry).
- **Memory for pileup arrays**: spliced-only pileups are very sparse. Store as `(transcript_id, bin, count)` records during scan, materialize to dense arrays only for transcripts that have spliced support. For a typical RNA-seq sample most transcripts will have few or zero spliced anchors, so memory should be modest.
- **Bin resolution**: per-exon is probably enough — finer resolution adds variance without much signal, and exon-level is exactly the granularity at which the failure presents itself.

## Bottom line

It's a real, complementary signal — specifically the *shape* likelihood that the current model is missing. The "spliced reads only" anchor idea is the part that makes it identifying rather than circular: it grounds the expected coverage in evidence that cannot itself be gDNA. Don't think of it as replacing effective length; think of it as adding a second axis (shape) on top of the existing scale axis (size). This is the most attractive of the seven ideas I proposed earlier because (a) you've spotted a clean non-circular anchor, (b) it generalizes naturally to nRNA spans and to hybrid capture, and (c) the cheapest version (prior scaling) requires no EM surgery.