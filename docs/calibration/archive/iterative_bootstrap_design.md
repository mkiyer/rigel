# Iterative calibration bootstrap — resolving the variance-model circularity under capture

**Status:** design (2026-06-15). Promotes the Stage-3 EB loop from "optional enhancement" to **required
architecture**: under hybrid capture the variance models *cannot* be trusted on a single pass, so
calibration must bootstrap and iterate. Builds on `per_node_deconv_hierarchy_design.md` §16-18 (the 3-tier
Bayesian per-node solve + the SCAM `var~mean` fitter in `variance_model.py`).

## 1. The circularity (the problem)

The IMPUTATION variance model must answer "how well do a region's boundaries predict the region?" — which
needs **triplets where the region itself is observable** (the region is the known answer). But under
hybrid capture:

- the **enriched** regions are **exons** (the high-mean regime), which are **not** count-observable (no
  contained gDNA signal to anchor) — so they yield **no** trainable triplets;
- the trainable triplets come almost entirely from **depleted** intergenic/intron regions (the low-mean
  regime).

So a variance model fit at initialization has seen only the **depleted** mean and variance range. Asking
it to predict at **enriched** means is pure extrapolation — and a confident extrapolation there is
**dangerous**: it would manufacture false certainty about the count at exactly the regions where the count
is least trustworthy. The circularity:

> trustworthy imputation model → needs enriched-region training examples → needs the enriched regions
> solved → needs a trustworthy imputation model.

## 2. The resolution — commit to iteration, with a *self-aware* variance model

Three mechanisms, together, break the circularity gracefully:

### (a) The variance model carries its own confidence (the safety lynchpin)

A `var~mean` prediction must come with a **confidence** = the local density of *training* data near the
query mean. Where the model **extrapolates** (a mean regime it has not been trained on), confidence → 0 and
the **count is down-weighted to the global-prior level** — so an untrained model can never give false
certainty. Concretely the per-node count precision becomes

```
τ_count(μ)  =  confidence(μ) · 1/var(μ)          (confidence ∈ [0,1]; →0 where unseen ⇒ count yields)
```

This is what makes a deliberately-uncertain first pass *automatic*: in Pass 0 the model is trained only on
depleted seeds, so `confidence(enriched μ) ≈ 0` ⇒ enriched nodes ignore the count and defer to the strand
(if stranded) and the global prior. No hand-set "make it uncertain" — the data density *is* the confidence.

### (b) Strand breaks the circularity (Pass 0 = the strand-only solve)

Strand-resolvable nodes need **no** count/imputation model: single-strand (TS_POS/NEG) nodes and
intergenic seeds are solved by the strand likelihood + Jeffreys alone. Crucially, in **stranded** data this
includes single-strand **exons** — so Pass 0 produces reliable region estimates **across the enriched
range**, giving the IMPUTATION model its first enriched-region training triplets (the strand-solved region
is the "answer"). This is the bootstrap that the count-only world cannot provide.

### (c) Enriched-mean *variance* is observable even pre-solve (your boundary-pair point)

An intron→exon→intron pattern gives a region anchored by **two observable boundaries**. The boundary-pair
disagreement is a **variance observation at an enriched mean** — without the region truth. So the DIRECT /
own-variance `var~mean` can be trained across the enriched mean range from Pass 0, even though the
imputation-*correctness* (boundary→region accuracy) still needs the region truth from (b)/iteration. This
widens the initial model's support so (a)'s confidence is non-zero across more of the range.

## 3. The bootstrap + iteration loop

```
PASS 0  (bootstrap — no trustworthy count model yet):
  • solve STRAND-RESOLVABLE nodes (single-strand + intergenic) by strand alone → reliable region
    estimates spanning depleted AND (in stranded data) enriched exons.
  • fit the DIRECT var~mean from observable own-variance + intron-exon-intron boundary-pair disagreement
    (covers the enriched mean range for VARIANCE).
  • fit the IMPUTATION var~mean from the now-solved triplets (strand-solved region = the answer).
  • everywhere the var~mean is extrapolating (low confidence), the count yields to strand + global.

PASS k≥1  (refine):
  • re-fit DIRECT + IMPUTATION var~mean from ALL solved nodes (far more data, now spanning the full mean
    range; the prior pass's estimates are the noisy "answers").
  • re-solve every node with the refined var~mean (confidence now high across the trained range),
    the propagation, and the global prior.
  • repeat until the var~mean curves and the deconvolved masses stabilize (a few passes; seed- and
    strand-anchored, so bounded — not open EM).
```

## 4. Why this resolves the circularity (and degrades gracefully)

- **Stranded data (the common case):** strand (b) solves the enriched exons in Pass 0 → the imputation
  model gets enriched training data immediately → one or two refine passes converge. The confidence (a)
  guards Pass 0 against the count's depleted-only training.
- **Unstranded data (ss≈½, the hard case):** strand can't break it; Pass 0 has no enriched region truth.
  Then (a) keeps the count uncertain at enriched regions (confidence→0), the global prior governs, and the
  boundary-pair variance (c) + iteration refine what they can. This is the honest limit — without strand
  *or* a trustworthy count, an enriched AMBIG node is genuinely under-determined, and deferring to the
  global gDNA baseline is the correct behavior, not a confident guess.
- **No false certainty, ever:** the confidence mechanism (a) is the invariant — the count is trusted only
  where the variance model has actually seen data at that mean. This is the danger you named, designed out.

## 5. Building blocks — status

Ready: the SCAM `var~mean` fitter (`variance_model.MonotoneVarMean`), the simplex/propagation solver, the
region↔boundary↔boundary partition + accumulator, the strand likelihood + Jeffreys (Pass-0 solver),
`count_observable_masks` (which nodes are strand/count observable). **To build:**

1. **Confidence in `MonotoneVarMean`** (this turn): `predict` returns `(var, confidence)` where confidence
   = the local density of training means at the query (∈[0,1], →0 on extrapolation). This is (a).
2. **The DIRECT/IMPUTATION builders by the right cases** (`variance_model`): DIRECT from observable
   own-variance + intron-exon-intron boundary pairs; IMPUTATION from solved triplets (region = answer).
3. **The bootstrap + iteration control flow** in `calibrate`: Pass-0 strand-only solve → fit var~mean →
   Pass-k refine → converge. Hooks: the strand-only solver path already exists (the node-class prior at
   `info`-high nodes); the re-fit reuses the builders.
4. **Convergence + stopping** (a few passes; track the var~mean / mass deltas).

## 6. Open questions

- **Stopping rule / max passes** — likely 2-3; measure the mass delta per pass.
- **Pass-0 strand-resolvable predicate** — exactly which nodes (info `I=N(2κ−1)²` above a threshold? or all
  single-strand?); the cleaner choice is "solve where the strand posterior is confident," reusing `w`.
- **Confidence bandwidth** — the kernel width for the local-data-density confidence (tie it to the SCAM
  knot spacing so it's consistent with the fit's own resolution).
- **Unstranded bootstrap** — is the boundary-pair variance (c) enough to seed anything useful, or do we
  just accept the global-prior fallback there? (Probably the latter.)
