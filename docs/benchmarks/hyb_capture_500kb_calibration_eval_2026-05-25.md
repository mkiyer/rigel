# Hybrid-Capture × Strand-Specificity Synthetic Benchmark — Calibration Diagnosis

**Date:** 2026-05-25
**Suite:** `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb`
**Rigel version:** 0.4.1 (post calibration overhaul)
**Scope:** 8 conditions = 2 capture (off/on) × 2 strand specificity (ss=0.50, 0.99) × 2 gDNA (none, high=1.0)
**Library composition:** 100k mRNA + 100k gDNA fragments (high) or 100k mRNA only (none). nRNA = 0 everywhere.

This report focuses on **calibration health** (gDNA prior, FL models, density anchors, locus-level priors) and the new **hybrid-capture exposure problem** surfaced by these scenarios. Transcript-level abundance is summarised but is *downstream* of the calibration findings and is not the primary subject.

Pipeline used: full rigel analysis driver (`python -m rigel.sim.analysis --sim-base …`) → quant on `sim_oracle.bam` with `--annotated-bam` per condition. Drill-down script: [scripts/debug/capture_eval_inspect.py](../../scripts/debug/capture_eval_inspect.py). Raw report: `analysis_report.txt` in the suite directory.

---

## TL;DR

| Capture | gDNA | ss | Spearman | gDNA→RNA leak | gDNA frac est (truth=0.5) | gDNA FL mean (truth=150) |
|---|---|---|---|---|---|---|
| off | none | 0.99 | **0.943** | — | 0.000 ✓ | fallback (249) |
| off | none | 0.50 | **0.938** | — | 0.000 ✓ | fallback (249) |
| off | high | 0.99 | **0.956** | 1.0% ✓ | **0.510 ✓** | **152.1 ✓** (+1.4%) |
| off | high | 0.50 | **0.939** | 1.5% ✓ | **0.516 ✓** | **151.9 ✓** (+1.3%) |
| on  | none | 0.99 | 0.397 ✗ | — | 0.000 ✓ | fallback |
| on  | none | 0.50 | 0.306 ✗ | — | 0.000 ✓ | fallback |
| on  | high | 0.99 | 0.287 ✗ | **13.1% ✗** | 0.449 (-10%) | 168.3 (+12%) |
| on  | high | 0.50 | 0.189 ✗ | **22.0% ✗** | 0.407 (-19%) | 168.6 (+12%) |

**Bottom line:** Without capture the new calibration is healthy across all dimensions tested here. With capture turned on, **every** calibration surface degrades in the same correlated way, and transcript quant collapses regardless of contamination or strand. The failure pattern points to a single missing model component: **per-transcript exposure / capture probability**.

---

## 1. What is working

### 1.1 gDNA fragment-length recovery (no capture)
Capture-off `gdna_high` runs recover the simulated gDNA FL mean (150 bp) within +1.4%: `gdna_fl_mean = 152.1` (ss=0.99) and `151.9` (ss=0.50). Quality gate reports `good`. RNA FL mean is +2.7% (256 vs 250), also in tolerance.

### 1.2 Global gDNA fraction (no capture)
The 1-D FL mixture (`gdna_fraction` in `summary.json/quantification`) returns 0.510 / 0.516 vs truth 0.500 — within 3% on 200k-fragment libraries.

### 1.3 Locus-level gDNA prior (no capture)
For `gdna_high_ss_0.99_capture_off` the per-locus `gdna_prior` tracks truth gDNA rate per locus tightly (locus 0 prior=0.93 vs estimate 0.98; locus 2 prior=0.96, est 1.00; locus 4 prior=0.78 est 0.82). Posterior `gdna_rate` recovers truth to within a few percent on all 10 loci.

### 1.4 RNA quantification (no capture)
Spearman ≈ 0.94 across all four capture-off conditions, independent of strand or gDNA contamination. Captured transcripts are not even relevant here, so quality is determined by the EM alone — which behaves.

### 1.5 nRNA stays at zero
Loci report `nrna = 0` everywhere; no false nRNA inflation in any nrna_none condition (post-fix acceptance check passes 8/8).

### 1.6 Pure gDNA classification (no capture)
Annotated-BAM fragment audit shows **98.5–99.0%** of gDNA fragments correctly identified as gDNA. The boundary-density / intron-anchor pathway is healthy when probes are absent.

---

## 2. What is failing — and why

### 2.1 Capture collapses the intron-based gDNA density anchor (root cause #1)

The synthetic suite has **no intergenic anchor** (the 500 kb chromosome is fully populated by gene loci), so calibration falls back to **INTRON** density. With capture off this works (mean_density ≈ 0.20, n_frags ≈ 27 k). With capture on the picture changes drastically:

| Condition | INTRON `mean_density` | INTRON n_frags | ratio to capture_off |
|---|---|---|---|
| gdna_high_ss_0.99_capture_off | 0.2016 | 27 588 | 1.00 |
| gdna_high_ss_0.99_capture_on  | 0.0029 | **394** | **0.014** |
| gdna_high_ss_0.50_capture_off | 0.1996 | 27 309 | 1.00 |
| gdna_high_ss_0.50_capture_on  | 0.0029 | **395** | **0.014** |

Capture removes ≈ 99% of intronic gDNA fragments. The estimated *intronic* gDNA density is then ~70× lower than capture-off, which is *factually correct* but is then applied as the per-locus gDNA prior on exonic regions where gDNA actually piles up at probe sites. Result: per-locus `gdna_prior` shrinks from 0.92–0.96 (capture_off) to 0.03–0.30 (capture_on) for the same loci. EM still rescues some of it through likelihood, but ~13–22% of true gDNA leaks into RNA.

**Diagnosis:** the intron anchor under-predicts exonic gDNA density when probes pull gDNA fragments onto captured regions. The calibration has no signal that distinguishes "globally low gDNA contamination" from "high gDNA contamination redistributed onto captured exons."

### 2.2 Capture biases the gDNA fragment-length mean upward (+12%)

Capture-on libraries report `gdna_fl_mean ≈ 168 bp` vs truth 150 bp. Only fragments that overlap a probe survive enrichment, and longer gDNA fragments are more likely to overlap a probe → length-biased survival. The FL quality gate still says `good`, so this bias is silent.

Effect on downstream: the gDNA likelihood factor (FL density) over-weights longer fragments, mildly favouring RNA assignment for shorter exonic-gDNA fragments. Secondary to issue 2.1 but reinforces the leak.

### 2.3 Per-transcript capture exposure is not modelled — quant collapses (root cause #2)

This is the dominant transcript-level failure mode. Aggregating expressed transcripts by whether they have any probe coverage:

| Condition | captured truth TPM | captured rigel TPM | uncaptured truth TPM | uncaptured rigel TPM |
|---|---|---|---|---|
| gdna_none_ss_0.99_capture_off | 646 607 | 648 051 | 353 393 | 351 834 |
| **gdna_none_ss_0.99_capture_on**  | 646 607 | **999 196** | 353 393 | **670** |
| **gdna_none_ss_0.50_capture_on**  | 646 607 | 998 402 | 353 393 | 860 |
| gdna_high_ss_0.99_capture_on  | 646 607 | 961 948 | 353 393 | 18 022 (mostly gDNA leak) |
| gdna_high_ss_0.50_capture_on  | 646 607 | 930 973 | 353 393 | 11 371 |

Under capture, **all real signal from uncaptured transcripts is gone from the BAM** — the simulator depleted them by ~500× via the `binding_per_base=10`/`off_target_weight=1.0` enrichment. Rigel has no exposure model, so it converts the surviving captured fragments straight into TPM. Captured TPMs inflate ~+55% to absorb the missing mass; uncaptured TPMs go to ≈0. This degrades Spearman from 0.94 → 0.19–0.40 even in the *no-gDNA* capture conditions, proving the failure is capture-driven, not contamination-driven.

Examples of within-locus capture distortion (`gdna_none_ss_0.99_capture_on`):
- `GENE0002.4` (cov_frac=1.0, truth TPM = 291 k) → rigel TPM = 454 k (+56%).
- `GENE0006.3` (cov_frac=0.99, truth = 269 k) → rigel = 412 k (+53%).
- `GENE0008.4` (cov_frac=0.0, truth = **219 k**) → rigel = 79 (-99.96%, signal annihilated).
- `GENE0010.4` (cov_frac=0.0, truth = 37 k) → rigel = 31.

### 2.4 Capture creates new false positives at captured non-expressed transcripts

Top FPs in `gdna_high_ss_0.99_capture_on`:
- `GENE0003.3` (captured, truth = 0) → 783 counts, TPM 7 818.
- `GENE0008.1` (captured, truth = 0) → 110 counts, TPM 2 667.
- `GENE0006.2` (uncaptured, truth = 0) → 526 counts (gDNA leak).
- `GENE0001.2` (uncaptured, truth = 0) → 122 counts (gDNA leak).

The captured non-expressed transcripts pick up gDNA reads at their exonic probe sites (the intronic-anchored prior is too low; EM cannot rule them out). The uncaptured non-expressed transcripts get gDNA reads that landed in their loci. Both classes of FP follow directly from §2.1 and §2.3.

### 2.5 No INTERGENIC anchor on the synthetic genome (evaluation gap)

In all 8 conditions `density_evidence.priors.INTERGENIC = None` — the synthetic 500 kb chromosome contains no genuine intergenic region. The intron anchor is the only signal driving the gDNA prior. This is a property of the *simulator* (the synthetic genome is gene-packed), not a rigel bug, but it amplifies §2.1. **For the production benchmark suite we need conditions with realistic intergenic padding** so we can characterise the intergenic-vs-intron path under capture.

---

## 3. Calibration-component scorecard

| Component | capture off | capture on | Comment |
|---|---|---|---|
| RNA FL mixture | ✓ | ✓ | +2.7–3.1% bias in both modes, independent of capture |
| gDNA FL mixture | ✓ (152) | ✗ (168, +12%) | capture-induced length-biased survival |
| Pool fraction π_pool | ✓ (0.51) | partial (0.41–0.45) | global FL mixture under-counts capture-relocated gDNA |
| INTRON density anchor | ✓ (0.20) | "correct but misleading" (0.003) | depleted of gDNA under capture, but the density family then mis-priors exonic regions |
| INTERGENIC anchor | n/a (genome) | n/a (genome) | evaluation gap |
| Per-locus gDNA prior | ✓ (0.13–0.96) | ✗ (0.02–0.65) | downstream of intron-anchor collapse |
| nRNA states stay off | ✓ | ✓ | acceptance passes 8/8 |
| Per-locus posterior gdna_rate | ✓ | partial — EM rescues some loci to 0.78–0.98 | likelihood overrules the bad prior in high-gDNA loci, but the leak still costs ~13–22% |
| Capture exposure model | n/a | **absent** | no per-transcript or per-locus exposure factor |

---

## 4. Where to focus effort

Prioritised by impact on the four-library matrix (capture {on,off} × ss {SS, US}):

### Priority 1 — Per-transcript / per-region exposure model for capture

Capture-on quant is broken even without gDNA. Without an exposure model TPMs are *capture coverage × truth*, not truth. Options worth scoping:

1. **Data-driven exposure detection.** Capture is detectable from the data: enormous (≥100×) ratio between per-base coverage of high-coverage exons vs introns/uncaptured regions; bimodal per-transcript coverage distribution. We can flag capture mode automatically.
2. **Per-transcript exposure factor inferred from data.** Empirical Bayes prior over an exposure scalar `s_t` per transcript (or per region), initialised by per-base coverage relative to a robust baseline. For uncaptured transcripts `s_t → 0`, the transcript is excluded from quant with an explicit flag rather than being assigned spurious near-zero TPM.
3. **Capture-region exposure family.** Mirror the existing INTRON/INTERGENIC density evidence with a `CAPTURED_EXON` family, fit from the bimodal coverage signature. Use it both for the gDNA prior (§2.1) *and* the RNA exposure.
4. **Optional user-supplied probe BED.** Cheap and unambiguous when available. Should not be the only path; pure data-driven detection is the strategic target ("rigel learns capture from data").

### Priority 2 — Capture-aware gDNA prior

Even if we fix exposure, the intron-anchor collapse (§2.1) still mis-priors gDNA at probe sites. Once we have a `CAPTURED_EXON` density family (Priority 1.3) we can:

- Fit gDNA density separately on intronic vs captured-exon regions.
- Use the **global FL π_pool** (which is robust under capture; off by only ~10%) as the per-library anchor when intron/intergenic anchors disagree by >10×.
- Add a *capture-detected* fallback that boosts the per-locus `gdna_prior` from the captured-exon family rather than the intron family.

### Priority 3 — Capture-induced gDNA FL bias

Once Priority 1 is in place, fit a probe-aware FL mixture: among capture-on libraries, the gDNA component is length-survival-biased (truncated at the short end). A simple capture-mode flag plus a re-weighted FL mixture should pull `gdna_fl_mean` from 168 back toward 150.

### Priority 4 — Simulator / evaluation gaps

- Add an INTERGENIC region to the synthetic 500 kb genome so we can validate the intergenic anchor under capture.
- Add nRNA × capture conditions (currently nRNA=none everywhere) once the above land.
- Extend the eval script to print `cov_frac` × TPM-error and captured/uncaptured aggregates automatically (currently in [scripts/debug/capture_eval_inspect.py](../../scripts/debug/capture_eval_inspect.py)).
- Acceptance thresholds currently treat 1.5% gDNA leak as the limit; for capture-on we need a separate threshold pegged to whatever the post-Priority-1/2 system can reach.

---

## 5. Quantitative appendix

### 5.1 Locus-level prior vs estimate, `gdna_high_ss_0.99_capture_on`
Captured & expressed loci 1, 5 — EM converges to gdna_rate=0.17, 0.20 with priors 0.03, 0.02 (recovered).
Uncaptured / gDNA-dominated loci 0, 2, 3, 4, 8 — truth gdna_rate ≈ 1.0; rigel reaches 0.93–0.98 with priors only 0.11–0.30 (under-prior, but EM rescues).
Locus 9 — truth gdna_rate ≈ 0.8; rigel 0.78 with prior 0.03 (under-prior).

### 5.2 Calibration anchors (capture on vs off)
| Anchor | metric | cap_off / cap_on | drop |
|---|---|---|---|
| INTRON `mean_density` | 0.2016 / 0.0029 | 70× |
| INTRON n_frags | 27 588 / 394 | 70× |
| ALL `mean_density` | 0.2002 / 0.0029 | 69× |
| intergenic_fraction (quant) | 0.315 / 0.005 | 70× |
| gdna_em_fraction (quant) | 0.195 / 0.444 | inverted — EM picks up the slack |
| total gdna_fraction | 0.510 / 0.449 | -12% |
| gdna_fl_mean | 152.1 / 168.3 | +12% |

### 5.3 Capture/uncaptured aggregate TPM
See §2.3 table.

### 5.4 Quant accuracy (Spearman / Pearson / MARD)
| Condition | Spearman | Pearson | MARD |
|---|---|---|---|
| gdna_none_ss_0.99_capture_off | 0.943 | 0.867 | 62 |
| gdna_none_ss_0.99_capture_on  | 0.397 | 0.532 | 49 |
| gdna_none_ss_0.50_capture_off | 0.938 | 0.877 | 53 |
| gdna_none_ss_0.50_capture_on  | 0.306 | 0.392 | 46 |
| gdna_high_ss_0.99_capture_off | 0.956 | 0.920 | 64 |
| gdna_high_ss_0.99_capture_on  | 0.287 | 0.384 | 293 |
| gdna_high_ss_0.50_capture_off | 0.939 | 0.818 | 72 |
| gdna_high_ss_0.50_capture_on  | 0.189 | 0.213 | 486 |

### 5.5 Fragment-level allocation accuracy (annotated-BAM audit, mRNA exact-transcript precision)
| Condition | mRNA exact | mRNA→gDNA | gDNA correct | gDNA→RNA |
|---|---|---|---|---|
| gdna_none_ss_0.99_capture_off | 0.813 | 0.000 | — | — |
| gdna_none_ss_0.99_capture_on  | 0.835 | 0.000 | — | — |
| gdna_high_ss_0.99_capture_off | 0.798 | 0.030 | **0.990** | 0.010 |
| gdna_high_ss_0.99_capture_on  | 0.797 | 0.028 | 0.869 | **0.131** |
| gdna_high_ss_0.50_capture_off | 0.782 | 0.046 | **0.985** | 0.014 |
| gdna_high_ss_0.50_capture_on  | 0.784 | 0.035 | 0.780 | **0.220** |

mRNA exact-transcript precision is **insensitive to capture** (~0.79–0.84) — i.e. when rigel decides a fragment is mRNA it picks the right transcript well regardless of capture. The failure is entirely in (a) routing exonic gDNA fragments to RNA, and (b) producing meaningful TPM in a library with capture-induced exposure heterogeneity.

---

## 6. Recommended next step

Implement a data-driven **capture-mode detector + per-region exposure family** (Priorities 1.1 + 1.3). It is the single change that unblocks both calibration (§2.1) and quant (§2.3). The synthetic suite already encodes the test signature — once exposure is modelled we expect:

- INTRON gDNA prior to remain low *and* a new CAPTURED_EXON gDNA anchor to come online with `mean_density` comparable to the capture_off intron anchor (≈ 0.2).
- Per-locus gdna_prior at captured loci to climb from 0.03 → ~0.5+.
- gDNA→RNA leak to drop from 13–22% back toward the 1–1.5% capture_off baseline.
- Uncaptured-transcript TPM to be explicitly flagged as low-exposure rather than silently produced as near-zero noise.
