# Effective-length shared-reference fix (FIX2) — roadmap & proof

**File:** `src/rigel/calibration/capture_eff_length.py`
**Prerequisite:** the junction-seam fix (landed 2026-07-09) — closes the `fl/span_full` inflation and adds
the splice-junction seams FIX2 reuses.
**Status:** designed + proven; implement in a fresh session, gated by one A/B.
**Published:** https://claude.ai/code/artifact/abb1c6da-388f-407c-a121-5b35cbcc1a7f
**Proof harness:** `scratchpad/efflen_fix_proof.py` (reproduces §4 end to end).

## At a glance

- **Problem:** each transcript is contracted against its **own** reference density `ρ* = G_c/E_c` (read on
  its own nodes). A mature on depleted exons reads a low `ρ*` (no contraction, `eff = fl`); its nascent
  parent, spanning an enriched region, reads a high `ρ*` (contracts) → mature ends up **longer** than
  nascent — the physically impossible inversion, now driven by the reference mismatch (not the junction
  geometry the first fix repaired, and not the calibration mass-attribution collapse).
- **Evidence (flagship, after junction fix):** ~2.5% of pairs still invert (synthetic geometry-only field),
  with `ρ*_mature ≈ 0.1` vs `ρ*_nascent ≈ 93.6`; `fl/span_full ≈ 0.93–0.99` confirms the junction fix
  works and the residual is *reference*, not geometry.
- **Fix (FIX2):** contract every transcript in a locus against a **single shared per-locus reference**
  `ρ_ref` — the same `G_c/E_c` the gDNA component already uses in `assemble_priors`.
- **Bonus:** puts transcripts on the same density scale as the gDNA component (removes a latent
  gDNA-vs-RNA arbitration inconsistency).

## 1. The problem — a per-transcript reference

The contraction writes each transcript's EM effective length as its enrichment-weighted footprint: node
supports scaled by enrichment relative to a reference `ρ_ref`. Today `ρ_ref = ρ* = G_c/E_c` is computed
**per transcript** on its own contained (exon) nodes, so a mature and its nascent parent are measured
against different rulers. Under **targeted** capture (only some exons probed):

- a **mature** with all exons on *depleted* regions reads a low, uniform `ρ*` ⇒ no contraction ⇒ `eff = fl`;
- its **nascent parent** spans the same locus but crosses an *enriched* region ⇒ high `ρ*` ⇒ contracts.

⇒ `eff_mature (uncontracted) > eff_nascent (contracted)`. Survives perfect (oracle) calibration — it is
purely the per-transcript reference choice.

## 2. Ramifications

- **Residual inversions** (~2.5%), concentrated on low-count depleted matures — small EM impact, but the
  invariant is violated and the guard test cannot be made strict.
- **Scale inconsistency with the gDNA component:** `assemble_priors` gives gDNA a *per-locus* reference;
  transcripts use a *per-transcript* one, so they compete on subtly different density scales.

## 3. The fix — one shared per-locus reference

```
Λ_t = Σ_{n∈t}  S_n · min(ρ_n / ρ_ref, 1)      ρ_n = m_n/S_n,  ρ_ref shared across the locus
   mature nodes  = exon regions + interior + splice-junction seams (imputed, from the landed fix)
   nascent nodes = all regions + genomic seams
   ρ_ref[locus]  = G_c/E_c over the locus's contained nodes  (== the gDNA component's reference)
```

The junction seams already make `Σ S_n = fl` for a spliced mRNA; the only change is that `ρ_ref` becomes
per-locus instead of per-transcript. **Capture-off stays bit-identical:** uniform field ⇒ `ρ_n = ρ_ref` ⇒
`w_n = 1` ⇒ `Λ_t = Σ S_n = fl`.

## 4. The proof

**4.1 Capture-off bit-identity.** Uniform ⇒ `w_n = 1` ⇒ `Λ_t = Σ S_n = fl_t` (junction seams restore
contiguity). Bedrock invariant preserved.

**4.2 Monotonicity theorem.** For a mature `m` and its nascent parent `n` sharing exons, shared `ρ_ref`:

```
Λ_nascent − Λ_mature = Σ_{introns I} (ρ_I / ρ_ref) · |I|  ≥ 0
```

The gap is exactly the enrichment-weighted intron content the nascent has and the mature lacks; → 0 under
strong capture (depleted introns), matching the biology (equal effective length under infinite capture on
a shared exon).

*Derivation, per intron between exons e, e+1:*
```
nascent: seam(e→I) + intron I + seam(I→e+1);   mature: seam(e→e+1)
gap = S_j·w(½(ρ_e+ρ_I)) + (|I|−L)·w(ρ_I) + S_j·w(½(ρ_I+ρ_{e+1})) − S_j·w(½(ρ_e+ρ_{e+1}))
```
Shared exon/junction terms cancel (identical `ρ_ref`); half-sum weights are in `w`'s linear regime so
`w(½(ρ_e+ρ_I)) + w(½(ρ_I+ρ_{e+1})) − w(½(ρ_e+ρ_{e+1})) = ρ_I/ρ_ref`, giving
`gap = (ρ_I/ρ_ref)·L + (ρ_I/ρ_ref)·(|I|−L) = (ρ_I/ρ_ref)·|I| ≥ 0`. ∎

**Robust to the `ρ_ref` choice:** if a captured exon saturates the cap (`ρ_e > ρ_ref`), the shared
exon/junction terms still cancel pairwise and the remaining nascent terms are ≥ 0 — the gap stays
non-negative for *any* shared reference. Sharing the reference is what matters, not its exact definition.

**4.3 Numeric verification** (`scratchpad/efflen_fix_proof.py`):

| scenario | Λ mature | Λ nascent | gap | Σ(ρ_I/ρ_ref)·\|I\| | exact? |
|---|--:|--:|--:|--:|:--:|
| 6-exon, capture exon1 | 502.9 | 512.9 | 10.00 | 10.00 | ✓ |
| 6-exon, capture 1,3,5 | 1701.7 | 1711.7 | 10.00 | 10.00 | ✓ |
| 6-exon, capture all | 3400.0 | 3410.0 | 10.00 | 10.00 | ✓ |
| micro-exons, capture exon1 | 701.9 | 707.9 | 6.00 | 6.00 | ✓ |

Every case: gap = enrichment-weighted intron length to machine precision, no inversions. The per-transcript
reference inverts on the same inputs.

## 5. Implementation plan

1. **Thread the per-locus reference in.** `transcript_capture_eff_lengths` runs per-transcript over the
   region arrays with no locus structure. Compute `ρ_ref[locus] = G_c/E_c` once (already computed for the
   gDNA component in `assemble_priors`) and pass a per-region→per-locus reference map into the contraction;
   unify so gDNA and transcripts share one reference (don't compute twice).
2. **Rewrite the contraction** as `eff_em = Σ S_n · min(ρ_n/ρ_ref, 1)` over the node set (regions +
   interior + junction seams, all already assembled), replacing `fl · (θ/ρ*)/span_full`.
3. **Preserve the multimapper shrinkage** `w = C/(C+1)` (toward no contraction on sparse contained gDNA)
   through the rewrite.
4. **A/B before shipping** — the per-transcript "contained ρ*" won a prior binding-sweep A/B (memory
   `efflen_node_mapping_contained_rhostar`). A/B the shared reference on the soft 3-pool surplus across the
   flagship set; confirm no abundance regression on captured matures. One change at a time.
5. **Make the invariant test strict** — assert `eff_em(nascent) ≥ eff_em(mature)` for *all* pairs plus the
   exact-gap identity on a synthetic field; regenerate goldens.
