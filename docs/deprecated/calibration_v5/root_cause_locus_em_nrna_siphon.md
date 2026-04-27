# Root-Cause Analysis: "Calibration 10× failure" is actually nRNA→gDNA siphon in the locus EM

**Status:** Diagnosed; ready to discuss fix.
**Date:** 2026-04-22
**Input data:** `/scratch/.../runs/human/mctp_vcap_rna20m_dna10m/rigel/annotated.bam` quantified by current (Option-B) build; outputs in `/scratch/.../runs/human_optionb/mctp_vcap_rna20m_dna10m/`.
**Per-region counts dump:** `/scratch/.../runs/human_optionb_diag/dna10m/region_counts.feather`
**Diagnostic scripts:** [scripts/debug/dump_region_counts.py](scripts/debug/dump_region_counts.py)

---

## Executive summary

The "10× calibration gap" in `vcap_mixture_p1_corrected.md` is **not** a calibration failure.  The calibration prior is approximately right.  The gap opens **inside the locus EM**, which inflates `gdna_prior` into `gdna_rate` by factors of 10×–330× per locus.  The mechanism is an **nRNA siphon**: the per-transcript nRNA component captures ≤3% of fragments everywhere, so unspliced reads from expressed genes (pre-mRNA, retained introns, intronic STAR alignments) have nowhere to go in the competing-component EM — the gDNA component absorbs them.

The three hypotheses from `vcap_mixture_p1_corrected.md` are ranked as follows in light of the data:
- **H1 (clip removal):** ruled out as primary — the clip removes ≈60% of the non-clipped expectation but that only moves calibration ≈2×, not 12×.
- **H2 (class-R absorption in calibration EM):** *partially true* (calibration λ_G is ≈3.6× higher than the clean intergenic estimate), but this shift is in the wrong direction for the reported failure and is not what matters downstream.
- **H3 (locus-EM over-allocation):** **confirmed primary.**  This was listed as a secondary concern.  It should be the primary target.

---

## Evidence

### 1. Intergenic gives a clean ground-truth λ_G reference

Intergenic regions (no known transcript overlap) are by construction pure gDNA.  Pooling per-region fragment counts for dna10m:

| region class | n_reg | Σ mappable (Gb) | Σ n_u | naive λ_G |
|---|---:|---:|---:|---:|
| intergenic | 30,566 | 1.102 | 246,910 | **2.24e-4 / bp** |
| intronic | 247,898 | 1.597 | 3,113,281 | 1.95e-3 / bp |
| exonic | 117,333 | 0.125 | 8,083,018 | 6.46e-2 / bp |

Intergenic λ_G is homogeneous across chromosomes (per-chrom median 2.0e-4, range 1–3e-4 for the 20 largest chroms — 40% of total mappable bp).  Treat **λ_G_true ≈ 2.2e-4/bp** as the reference.

### 2. Calibration λ_G is 3.6× high (intronic contamination pulls it up)

```
Calibration-fit λ_G              = 8.02e-4 / bp   (3.6× the intergenic reference)
Σ (λ_G_cal · mbp) across genome  = 2.29 M fragments
Σ min(λ_G_cal · mbp, n_total)    = 0.886 M fragments   ← summary.total_expected_gdna
```
The calibration EM's "expressed-class" LogNormal absorbs some low-rate intronic regions (which contain genuine intronic RNA plus gDNA) as class R, so λ_G converges to a value between the intergenic truth and the intronic rate.  This is what `p1_corrected.md` called out, and it is a real ≈3× bias.

### 3. The locus EM turns 0.886M into 8.66M (13×)

```
calibration total_expected_gdna  =   0.886 M
Σ (gdna_prior × n_locus)         =   0.699 M   ← per-locus prior expectation
Σ (posterior gdna)               =   8.663 M   ← quant.gdna_total
Σ posterior / Σ prior_E_gdna     ≈   12.4×
```

By locus-size bucket (dna10m):

| size | n_loci | mean_prior | mean_rate(post) | mRNA% | nRNA% | gDNA% |
|---|---:|---:|---:|---:|---:|---:|
| <1 kb | 3,712 | 0.217 | 0.894 | 43% | 0.7% | 56% |
| 1–10 kb | 6,332 | 0.311 | 0.745 | 60% | 2.3% | 37% |
| 10–100 kb | 11,527 | 0.269 | 0.642 | 62% | 2.1% | 36% |
| 0.1–1 Mb | 4,412 | 0.236 | 0.644 | 53% | 2.8% | 44% |
| 1–10 Mb | 106 | 0.251 | 0.732 | 35% | 2.7% | 62% |
| >10 Mb | 2 | 0.081 | 0.633 | 47% | 0.08% | 53% |

The blow-up is universal across locus sizes — no size escapes it.  Extreme per-locus examples:

| locus | span | n_tx | total | mrna | **nrna** | gdna | prior | rate | blow-up |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 (mega) | 53.4 Mb | 6,883 | 296,718 | 146,190 | **225** | 150,303 | 0.053 | 0.507 | 9.6× |
| 18122 | 468 kb | 242 | 50,236 | 15,786 | **19** | 34,431 | 0.005 | 0.685 | **145×** |
| 12781 | 310 kb | 178 | 26,207 | 263 | **57** | 25,887 | 0.009 | 0.988 | 110× |
| 3985 | 28 kb | 160 | 12,186 | 5,126 | **73** | 6,987 | 0.002 | 0.573 | **331×** |

**The signature is identical everywhere:** priors in the 0.01–0.1 range are driven to posteriors of 0.5–0.99, and the nRNA component captures <3% of fragments.

### 4. Global quant vs ground-truth gDNA across the titration

| lib | cal λ_G | π_soft | quant gDNA | quant gDNA_frac | implied true gDNA (λ_G_intergenic × mbp_all) | quant / true |
|---|---:|---:|---:|---:|---:|---:|
| dna00m | 1.05e-5 | 0.920 | 0.19 M | 1.4% | ≈0.03 M | ≈6× |
| dna10m | 8.02e-4 | 0.863 | 8.66 M | 37.5% | ≈0.64 M | ≈13× |
| dna80m | 5.45e-3 | 0.865 | 59.54 M | 76.8% | ≈(scaled) 5 M | ≈12× |

(Intergenic λ_G scales approximately linearly with DNA input in a spike-in titration, so the "implied true gDNA" column extrapolates the dna10m intergenic estimate.)

The **quant gDNA is 10–13× higher than even a generous upper-bound estimate of true gDNA.**

---

## Mechanism: why the locus EM inflates gDNA

Per-locus, rigel's EM has three competing fragment classes for each read:
1. **mRNA**: transcript isoform × coverage model × spliced-RNA FL model.
2. **nRNA**: one "nRNA entity" per unique intronic/gene-body span × genomic coverage × RNA FL model.
3. **gDNA**: flat rate λ_G × locus span × gDNA FL model.

The observed symptoms:
- nRNA globally captures **2.4%** of fragments (0.51 M of 21.9 M) — implausibly low for a prostate cancer cell line where nascent/retained-intron signal is typically 10–20% of polyA-selected reads and far more for total RNA.
- mean `gdna_rate` posterior per locus is 0.28–0.89 across size buckets, driven by the locus EM adding mass that was not in the prior.
- The allocation gradient is `(gDNA ↑, nRNA ≈ 0)` — consistent with **unspliced reads being routed to gDNA because nRNA has too low a data likelihood** (or too few entities, or too tight a coverage profile, or too low a prior weight).

This matches the mega-locus observation in `p1_corrected.md`: 53 Mb, 6,883 tx, 225 nRNA out of 297 K reads is absurd — a locus this size at a VCaP gene-rich region should have ≥30 K nRNA/pre-mRNA fragments.

### Plausible proximate causes (ranked)

1. **nRNA FL model is being chosen equal to the RNA FL model rather than unspliced-specific** → unspliced reads score poorly against mRNA (no junction match) and equally against nRNA and gDNA, so the flat gDNA rate dominates.  Worth checking whether Option B's scorer corrections introduced an FL-model mismatch for nRNA.
2. **nRNA entity coverage model is over-specified** (requires the read to be inside a specific intron span) → a read that spans an intron boundary or falls in a gene-body region not covered by any enumerated nRNA entity has `P(read | nRNA_i) ≈ 0`, leaving only gDNA.
3. **gDNA likelihood is too generous** — `λ_G × locus_span_bp × P(FL | gdna_fl_model)` can exceed the nRNA likelihood even inside a transcribed gene, because it's a flat rate over the whole locus span with no penalty for landing inside an annotated gene body.
4. **Calibration's per-locus prior γ_ℓ is inflated via clipping** in regions where `λ_G × mbp > n_total` (lots of intronic regions) — but this makes γ_ℓ *smaller*, which would pull the posterior toward nRNA, not gDNA.  So this probably isn't the driver.

---

## Implications for the user's stated priorities

> "Calibration appears to be dramatically failing, and now that we have corrected the gdna EM issues, we need to resurrect calibration and make it work properly."

The diagnostic is that **calibration is not the primary failure mode right now.**  Its λ_G bias is ≈3.6×; the locus-EM over-allocation is ≈12×.  Fixing calibration alone will move the quant by at most ≈3×, not 12×.

Two orthogonal fixes are needed:
- **P0 (locus EM / nRNA):** Diagnose why the nRNA component is collecting <3% of fragments globally, despite representing gene-body / intronic real-RNA signal.  Most likely a scoring / entity-enumeration issue that slipped through the Option B changes.
- **P1 (calibration):** Tighten the "class R absorbs intergenic-like regions" bleed.  The cleanest proposal is to **seed λ_G from intergenic-only, high-mappable regions** and use the EM only to refine, not re-estimate, it.  The intergenic pool is genuinely pure gDNA and the per-chrom estimate is stable — so an intergenic-anchored prior for λ_G would remove the 3.6× bias.

If we only have bandwidth for one, P0 has ≥4× the leverage of P1.

---

## Proposed next step (for user to direct)

**Option A — diagnose the nRNA siphon (P0) first.**  Write a locus-level diagnostic that, on a representative locus (e.g., locus 18122: 242 tx, 50 K fragments, 145× blow-up), dumps:
- per-fragment `(mRNA_i, nRNA_j, gdna)` likelihoods and final posterior,
- per-nRNA-entity count and coverage statistics,
- the RNA vs gDNA FL model comparison at that locus.
This should reveal in minutes whether (a) nRNA entities are missing, (b) nRNA likelihood is zero-floored, or (c) gDNA dominates by flat-rate volume.

**Option B — stabilize calibration (P1) via intergenic anchor.**  Replace `_robust_initial_lam_G` with a pooled intergenic-region MLE, and add an intergenic-anchored prior term to `_m_step` to keep λ_G near the pure-gDNA rate.  Quick to implement; effect is bounded to ≈3×.

**Option C — both in parallel.**  P0 is a diagnostic pass (no code changes) so it doesn't block P1.

Recommendation: **Option C** — start P0 diagnostic now (no risk to code), queue P1 fix pending P0 findings.
