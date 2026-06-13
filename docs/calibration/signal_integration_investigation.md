# Signal integration — methodological investigation (strand × count trust)

**Status:** investigation, 2026-06-13. Written after a gutting mistake (replacing `solve_node` with a flat
density average — reverted) and a design objection: the integration of the strand and count signals is the
crux, it must be done *inside* the elegant solver, and the lever is **calibrated uncertainty / trust**,
not a new algorithm. This doc establishes the framework, grounds it in a measurement, and proposes a
principled, calibratable trust scheme that generalizes the old `w = I/(I+I₀)`.

## 0. The grounding measurement (flagship, ss 0.99, capture-on)

Per-node estimate error vs oracle gDNA fraction, by node class:

| node class | n | oracle | strand bias / MAE | count bias / MAE |
|---|---|---|---|---|
| single-strand (both avail) | 940 | 0.925 | −0.123 / 0.134 | −0.258 / 0.263 |
| ↳ EXON (on-target) | 524 | 0.865 | −0.150 / 0.172 | **−0.459 / 0.468** |
| ↳ INTRON (off-target) | 416 | 1.000 | −0.097 / 0.097 | **−0.005 / 0.005** |
| AMBIG (count only) | 141 | 0.967 | — | −0.383 / 0.383 |

**Findings.** (a) Strand is ~unbiased and uniformly decent (MAE ≈ 0.13). (b) The count's reliability is
**node-dependent**: excellent where count-observable (introns/intergenic, MAE 0.005), catastrophic where
*imputed across the capture enrichment discontinuity* (exons, MAE 0.47, bias −0.46). (c) The count error is
a **systematic bias**, not sampling noise. (d) AMBIG nodes have no strand and a biased count (−0.38) — the
leak source.

## 1. Are we operating with log-likelihoods? Yes.

`simplex.solve_node` maximises a log-posterior over the pie `(f₊, f₋, f_g)` that is a **sum of
log-likelihood terms** (strand mixture BB + a Gaussian count term `−½·c·(f_g − g_count)²` + the
spliced lower bound + the gDNA prior). Integrating evidence = **adding log-likelihoods** = multiplying
likelihoods — the standard Bayesian fusion. The **trust** of each signal is the **curvature (precision)**
of its log-likelihood at the optimum. So "trust strand over count" is, mathematically, "give the strand
term more curvature than the count term." Nothing about the algorithm needs to change — only the
precisions.

## 2. The old `w = I/(I+I₀)` IS Gaussian log-likelihood fusion

`g = w·g_strand + (1−w)·g_count`, `w = I/(I+I₀)`, `I = N·(2κ−1)²` is **exactly** the posterior mean of
fusing two Gaussian signals:

```
strand:  mean g_strand,  precision ∝ I  = N·(2κ−1)²     (grows with depth × discriminability)
count:   mean g_count,   precision ∝ I₀ = const (10)     (FIXED — independent of the count's own data)
fused mean = (I·g_strand + I₀·g_count)/(I+I₀) = w·g_strand + (1−w)·g_count
```

The decisive design choice hidden in the old method: **the count's precision is a fixed cap `I₀`, never its
own data variance.** That is *why* it was robust — it never let a confident-but-wrong count dominate. It's
a crude, one-number acknowledgement that the count is misspecified.

## 3. Why naive inverse-variance regressed (the lesson)

Setting `count precision = 1/var_count` (the `var~mean` sampling variance) lets a count with small sampling
variance dominate — even where it is badly **biased** (exons, −0.46). Variance ≠ bias. A misspecified
estimator can be *precise and wrong*. The flat-average gut built on this and regressed (pool gDNA 0.686 →
0.628 on the flagship; further from the true 0.75). **The count's trust must reflect its bias /
misspecification, not just its noise.**

## 4. The principled framing: the count likelihood is MISSPECIFIED

The count is a biased estimator under capture (the gDNA-density imputation crosses the enrichment
discontinuity). The literature tool for fusing a misspecified likelihood is the **tempered / power
likelihood** (generalized Bayes): raise it to a power `β ∈ [0,1]`,

```
log L_total = log L_strand  +  β · log L_count  +  log prior
```

`β = 1` → full trust; `β → 0` → ignore the count; `β` between → "strongly penalize the count so the strand
is much more likely," exactly the requested lever. Tempering multiplies the count term's curvature by `β`,
so the **effective count precision = β · (1/var_count)**. The old `I₀` cap is the special case
`β = I₀·var_count` (i.e. choose `β` so the count contributes exactly `I₀` worth of precision). Tempering is
the *general* knob; the cap is one calibration of it. (Equivalent dual view: **variance inflation**
`var_count → var_count/β`, adding the bias to the uncertainty.)

## 5. The key generalization — `β` (count trust) must be PER-NODE

§0 shows the count's reliability is **not uniform**: trust it on introns (MAE 0.005), distrust it on
imputed exons (MAE 0.47). A single global `β` (or `I₀`) is a compromise that under-trusts good introns and
over-trusts bad exons. The generalization that beats the old method:

```
β_node  high  where the count is a DIRECT, reliable measurement   (count-observable: introns, intergenic)
β_node  low   where the count is IMPUTED across the capture jump   (exons, AMBIG exons)
```

Signals available to set `β_node` *without oracle* (all already computed):
- **count-observability** (`region_count_observable` / `boundary_count_observable`) — a direct vs imputed flag.
- the **`var~mean` curve** — larger imputation variance ⇒ lower `β` (captures *some* of the unreliability).
- **capture-class crossing** — whether the node's density was imputed from a different enrichment class
  (the exon-imputed-from-intron case, the −0.46 bias) ⇒ low `β`.

This is the same shape as the old `w` (a strand-trust gradient) but now the **count side carries its own
data-driven reliability**, so the fusion is honest both ways.

## 6. Where propagation earns its keep

The integration above is per node. **Propagation** adds the missing ingredient for AMBIG exons: they have
*no strand* and a *biased count* (−0.38). But their single-strand EXON neighbours have a **strand-derived**
gDNA density that is far better (bias −0.15 vs the count's −0.46). The RTS propagates that strand-derived
density (a signal carrying a value **and** a precision) along the exon-class chain into the AMBIG node,
where — *because `β_node` for the AMBIG count is low* — it out-trusts the local biased count. So:

- propagation **carries strand-quality information to nodes that have no strand of their own**;
- the **per-node count trust `β`** is what lets that propagated signal win over the local biased count.

Both are needed: propagation without the count-trust calibration lets the biased count resist (the AMBIG
under-call we saw); the calibration without propagation leaves AMBIG with nothing better than its count.

## 7. Recommended scheme (to review, then implement inside `solve_node`)

1. **Keep `solve_node`** as the integrator (log-likelihood sum on the pie). No new solver.
2. **Strand term**: precision `I = N·(2κ−1)²`, unbiased — unchanged.
3. **Count term**: tempered precision `β_node · (1/var_count)`, with `var_count` from the `var~mean` curve
   and `β_node` a **data-driven reliability** (count-observable → high; imputed-across-capture → low).
   Calibrate the `β_node` mapping on the benchmark using §0-style strand-vs-count error by node class
   (we *measured* the reliabilities — they become the calibration targets, not magic numbers).
4. **Propagation**: the RTS carries the gDNA-density signal (value + precision) along enrichment-class
   chains; the strand-derived seed densities reach AMBIG nodes and, with low `β` on the AMBIG count, govern.
5. The pie keeps the no-over-subtraction safety; the gDNA prior keeps the no-evidence default at gDNA.

## 8. Open questions (research, before implementing)

1. **Calibrating `β_node`.** Map (observability, `var~mean`, capture-class-crossing) → `β`. Is a 2-level
   `β` (direct vs imputed-across-capture) enough, or a continuous function? Fit on the benchmark; validate
   that the fused estimate's MAE beats both signals alone.
2. **Tempering vs variance-inflation vs cap** — three equivalent encodings; pick the one that composes
   cleanly with the propagation precision (likely tempering, so the count term's curvature is just scaled).
3. **Does the propagated strand-density actually reach the AMBIG exons?** §0 says single-strand exon
   neighbours carry −0.15-quality density; confirm the exon-class chain connects them to the AMBIG exons
   (vs being seed-starved) — measure the propagated density's error at AMBIG nodes directly.
4. **Avoid double-counting a node's own strand** when it is both a seed (its strand-derived density enters
   the field) and re-uses the field as count evidence in its own `solve_node`. The exclude-self rule (BP
   message excludes the target's own contribution) is the clean fix.
5. **`I₀` vs `β`** — is the old global `I₀=10` already near-optimal for single-strand nodes, with the
   per-node `β` only needed to *raise* trust on count-observable introns and *lower* it on imputed exons?
   Quantify the headroom.

## 9. What this preserves

The simplex (the pie, no over-subtraction) and the propagation (BP on a chain, order-independent) are
kept — they generalize the old per-node combine to a spatial one. The only thing being *designed* is the
**count's trust**, which the old method hard-capped at `I₀` and which §0 shows we can calibrate per node
to beat. No solver is replaced.
