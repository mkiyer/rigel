# Capture-aware effective lengths for all EM components (design)

**Status:** design for review. 2026-06-10 (rev 2). Fixes a confirmed **major bug**: the hybrid-capture
IPR effective-length contraction is applied to the **gDNA component only**; mRNA and nRNA components use
the capture-blind full FL-marginal length, so under capture the gDNA component is artificially
concentrated and out-competes nRNA (and mature) for capture-enriched reads.

> **Rev 2 correction.** An earlier draft proposed deriving a new per-region enrichment `e_r` from total
> coverage density. That is wrong: transcript expression spans ~10⁴-fold dynamic range, which would swamp
> the probe signal. **gDNA is already the correct per-region enrichment** — it is source-uniform within a
> locus, so its deconvolved per-region mass reflects *only* the probe pattern, carrying no expression
> structure. The IPR machinery that already contracts the gDNA component is the whole solution; we simply
> apply it to each transcript's region set. No new readout, no new constants.

## 1. The bug (confirmed, code + data)

- **gDNA component** (`priors.assemble_priors`): `gdna_eff_len` = the Laplace-smoothed IPR of the
  per-region deconvolved gDNA mass over the locus's regions — capture-aware, contracts under capture.
- **mRNA / nRNA components** (`pipeline._setup_geometry`): `effective_lengths` =
  `rna_fl.compute_all_transcript_eff_lens(exonic_length)` — capture-**blind** full length;
  `effective_lengths_em = None`, so the EM uses these for every transcript row.

In the EM a component's per-position rate ∝ `abundance / eff_len`. Under capture, reads concentrate on
probed regions; the gDNA component (short, contracted) wins them, the nRNA component (full genomic span,
uncontracted) is diluted and loses. Empirically, on the siphon loci the nRNA effective length is **5–12×
the gDNA effective length** (locus 18: 38,448 vs 3,894), and of 190,588 nascent reads only 4,620 (2.4%)
reach nRNA — 96,520 siphon to gDNA, ~89,000 to mRNA.

## 2. The principle — extend the existing gDNA IPR to every component

The per-region deconvolved gDNA mass `m_r` (`gdna_region`, transport-corrected) and its FL-aware support
`L_r` (`gdna_geom_len`) are **already computed** and are the clean enrichment readout (gDNA uniform ⇒
expression-immune). The gDNA component's effective length is already the Laplace-smoothed IPR of `(m_r,
L_r)` over the **locus's** region set. The fix: compute the same IPR over **each transcript's** region
set:

| component | region set `S_t` |
|---|---|
| gDNA | all regions in the locus (unchanged) |
| nRNA | all regions its (single, full-span) exon overlaps — exons **and** introns |
| mRNA | the regions its exons overlap — exon regions only |

Both transcript cases fall out of one routine — *the regions a transcript's exon-blocks overlap* —
because nRNA rows are single-exon spans (covering introns too) while mRNA exon-blocks skip introns.

## 3. The formula — identical to the gDNA IPR, restricted to `S_t`

```
              ( G_t + 1 )²
 eff_len_t = ───────────────── ,  capped at  span_t
             supp_t + (2·G_t + 1)/span_t

 G_t   = Σ_{r∈S_t}  m_r            (gDNA mass over the transcript's regions)
 supp_t= Σ_{r∈S_t}  m_r² / L_r
 span_t= Σ_{r∈S_t}  L_r            (the FL-aware genomic/exonic span)
```

This is `priors.assemble_priors`'s existing gDNA formula with the sum restricted from the locus to `S_t`.
Behaviour (all inherited, no new logic):
- **No gDNA over the transcript** (`G_t = 0`): `eff_len_t = span_t` — *no contraction*, the genomic/exonic
  span. (Exactly what you described: zero coverage ⇒ span; the Laplace `+1` gives this with no special
  case.)
- **gDNA present & concentrated** (capture): contracts toward the probed regions, ∝ the gDNA evidence.
- An nRNA whose span ≈ the gene gets `S_t` ≈ the locus regions ⇒ `eff_len_nRNA ≈ eff_len_gDNA` → the two
  compete on equal footing (the bug's asymmetry is gone). mRNA contracts over its exon regions only.

Sparse/zero, normalization, and the coupling of coverage→contraction are **already handled by the IPR +
Laplace smoothing** — there is nothing new to design there (this was my confusion in rev 1).

## 4. The O(incidence) algorithm — one region pass per locus, all transcripts at once

`assemble_priors` already iterates regions × overlapping loci to project `m_r`/`m_r²/L_r`/`L_r` to loci.
Extend that single pass to also accumulate per **transcript**: for each region `r` and each transcript `t`
whose `S_t` includes `r`, add `m_r → G_t`, `m_r²/L_r → supp_t`, `L_r → span_t`. After the pass,
`eff_len_t = (G_t+1)² / (supp_t + (2G_t+1)/span_t)`, capped at `span_t`. Cost = `O(Σ incidences)` =
`O(R · transcripts-per-region)` per locus — one pass, every transcript simultaneously (your design). The
incidence list (`region → transcripts whose S_t includes it`) is **annotation-only** (transcript exons vs
the region partition), precomputed once at index/geometry build; only `m_r` is per-sample.

## 5. Implementation plan

1. **Transcript→region incidence + `L_r`** — precompute at geometry build from `index.t_df` exons vs the
   region partition (one overlap pass; reuses `priors._project_regions_to_loci`'s math). Stored on the
   geometry/index.
2. **Per-transcript contracted eff-len** — a routine taking the calibration per-region gDNA mass
   (`gdna_region` transported, `gdna_geom_len`) + the incidence, returning `eff_len_t` for every
   transcript via §3/§4. Naturally lives beside `assemble_priors` (it already holds `gdna_region`).
3. **`pipeline`** — set the geometry's `effective_lengths_em` to these (currently `None`); the gDNA
   component keeps its `gdna_eff_len` (the same formula, locus region set).
4. **EM** — already consumes `effective_lengths_em`; **no C++ change**.

## 6. Remaining open details (minor; not blockers)

- **6.1 mRNA exon support.** `L_r = gdna_geom_len` is the *genomic* FL-aware support; a mature fragment is
  spliced, so its true per-exon-region support uses the spliced neighbourhood, not the genomic one. Using
  the genomic `L_r` over exon regions is a v1 approximation (good when region≈exon); a spliced-aware `L_r`
  for mRNA is a later refinement. (Optimize later, per your note.)
- **6.2 Partial overlap.** Region boundaries align with exon/intron edges, so a transcript's exons
  usually *fully* cover the regions they touch ⇒ summing `m_r` over `S_t` is exact. Where a region is only
  partially covered (a boundary set by another transcript), the exact form weights `m_r` by the overlap
  fraction; v1 uses whole-region `m_r` (negligible error). Note for later.
- **6.3 Output/TPM length.** The EM uses the contracted `effective_lengths_em`; the reported abundance/TPM
  can keep the FL-marginal `effective_lengths` (convention) — separable, decide at wiring.
- **6.4 Acyclicity.** `m_r` is the calibration's per-region gDNA mass (pre-EM), so the contraction is a
  fixed EM input — no feedback loop. Confirm in the pass order.

## 7. Validation

- **Primary:** the nascent siphon collapses — `gdna_none ss0.99 capture_on`: nRNA absorbs its nascent
  (gdna_observed 96,520 → ≈0; nRNA 4,620 → ≈190,588).
- Capture-off ≈ inert (~1e-7 relative shift — the residual real gDNA-density variation; the Laplace `+1`
  keeps the factor ≈1, so it is not exactly bit-identical but practically so).
- **Validated (quick suite):** the nascent siphon collapses (nRNA 2.4%→93% of nascent), but with fair
  eff-lengths nRNA and gDNA — both unspliced-genomic — become **confounded**, so real gDNA also leaks into
  the resurrected nRNA component (the gDNA→RNA net leak rises under capture). Mature RNA stays controlled.
  This is a **real identifiability challenge uncovered** (nRNA↔gDNA), not a regression: the eff-length
  asymmetry had been masking it by handing gDNA every unspliced read. The separate discriminator
  (FL / strand) between nRNA and gDNA is the follow-up; report leak at the **3-pool** (gDNA/nRNA/mRNA)
  level, watching the mature-RNA false-positive rate.
- Complex QUICK suite: gDNA→RNA leak and nascent siphon both improve; pure-RNA FP stays ≈0.
- Golden changes for captured/nascent scenarios (intended); regenerate + spot-check direction.

## 8. Scope & risk

Cross-cutting (precompute incidence → per-transcript IPR from `gdna_region` → `effective_lengths_em` → EM)
but **no C++ change**, **no new constants**, **no new readout** — it reuses the existing gDNA IPR + the
existing overlap projection. It changes mRNA *and* nRNA EM eff-lengths under capture (the intended fix),
should improve the nascent siphon *and* mature quantification under capture, and capture-off is unchanged.

**Ready to implement:** the `e_r` question that gated rev 1 is resolved (it is the existing gDNA mass).
The only open items (§6) are minor v1 approximations, not blockers.
