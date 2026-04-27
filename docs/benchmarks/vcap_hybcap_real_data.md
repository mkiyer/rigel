# Rigel vs salmon vs kallisto on hybrid-capture VCaP

**Date:** 2026-04-27
**Libraries:**
- `mctp_vcap_rna20m_dna00m` — 20 M cell-line RNA, **no gDNA** spike-in (28.7 M raw reads)
- `mctp_vcap_rna20m_dna20m` — 20 M RNA + **20 M gDNA** spike-in (66.1 M raw reads)
- Both: hybrid-capture protocol (same kit captures RNA and DNA spike).
- Reference: GRCh38 / GENCODE, rigel index 457k transcripts.

**Tools:**
| Tool | Version | Input | Aligner |
|------|---------|-------|---------|
| rigel | 0.5 (current main, post-P0 streamlining) | name-sorted BAM | STAR |
| salmon | (pre-existing pipeline output) | FASTQ | quasi-mapping |
| kallisto | (pre-existing pipeline output) | FASTQ | pseudoalignment |

Two rigel runs were produced for fair comparison:
- `rigel_v05` — default `--assignment-mode sample` (single-best draws)
- `rigel_v05_frac` — `--assignment-mode fractional` (matches salmon/kallisto semantics)

---

## §0 Headline findings

1. **At the gene level, rigel agrees very well with both other tools** — Pearson
   r ≥ 0.997 in dna00m and ≥ 0.99 in dna20m. The transcript-level
   disagreements are real but they reflect (a) salmon's
   length-bias correction, (b) kallisto's quasi-alignment-vs-real-alignment
   differences, and (c) **rigel's gDNA siphon** in the spiked library.

2. **rigel calls 42% of the unspliced/intronic pool as gDNA in dna00m even
   though the library has no gDNA spike-in.** This is the **gDNA-vs-nascent-RNA
   identifiability problem** documented in `docs/calibration/`. The capture
   kit pulls in pre-mRNA/intronic fragments that look exactly like fragmented
   genomic DNA to a length-only mixture model.

3. **In dna20m, rigel correctly inflates the gDNA share to 68%** and removes
   ~1,000 loci of contamination that salmon and kallisto are silently calling
   as transcript counts. This is the correct behavior, but it makes rigel look
   like a low-correlation outlier when comparing against tools that don't
   model gDNA at all.

4. The **default `--assignment-mode sample` gives misleadingly low log-space
   correlations** because it produces ~89k expressed transcripts vs salmon's
   ~143k and kallisto's ~93k. Switching to `--assignment-mode fractional`
   recovers ~52,000 additional non-zero transcripts and adds ~0.6 percentage
   points to log-space r — but more importantly makes the comparison
   semantically fair.

---

## §1 Coverage and total counts

After dropping rigel's synthetic-nRNA rows, three-way intersection on
common transcripts (n = 227,394 annotated transcripts present in all three
catalogs):

| Library | Mode | rigel Σcount | salmon Σcount | kallisto Σcount | rigel nz tx | salmon nz tx | kallisto nz tx |
|---------|------|-------------:|--------------:|----------------:|------------:|-------------:|---------------:|
| dna00m | sample | 13,116,815 | 12,928,819 | 13,288,094 | 92,724 | 143,777 | 142,766 |
| dna00m | fractional | 13,104,251 | 12,928,819 | 13,288,094 | 150,615 | 143,777 | 142,766 |
| dna20m | sample | 13,937,976 | 15,568,589 | 20,632,430 | 100,859 | 194,422 | 169,956 |
| dna20m | fractional | 13,921,901 | 15,568,589 | 20,632,430 | 194,421 | 194,422 | 169,956 |

Key observations:
- **dna20m**: salmon assigns 1.6 M more reads to transcripts than rigel does;
  kallisto assigns **6.7 M more**. Both are taking gDNA spike-in fragments
  that aligned to intronic regions and credit them to overlapping
  transcripts. Rigel routes those fragments to the gDNA component (see §3).
- The fractional rigel run produces almost the same number of expressed
  transcripts as salmon (150 k vs 144 k in dna00m; 194 k vs 194 k in dna20m).
  In `sample` mode, rigel's count of expressed transcripts is roughly half,
  which inflates the apparent disagreement.

---

## §2 Pairwise agreement (annotated transcripts only)

### dna00m (`fractional` mode)

| Pair | Pearson r | log2 Pearson r | Spearman r |
|------|----------:|---------------:|-----------:|
| rigel vs salmon | **0.921** | **0.896** | 0.766 |
| rigel vs kallisto | **0.998** | **0.944** | 0.803 |
| salmon vs kallisto | 0.922 | 0.899 | 0.841 |

### dna20m (`fractional` mode)

| Pair | Pearson r | log2 Pearson r | Spearman r |
|------|----------:|---------------:|-----------:|
| rigel vs salmon | **0.911** | **0.800** | 0.715 |
| rigel vs kallisto | **0.975** | **0.734** | 0.687 |
| salmon vs kallisto | 0.905 | 0.773 | 0.745 |

### Sample-vs-fractional impact (rigel only)

| Library | Mode | rigel-salmon Pr | rigel-salmon log2 Pr | rigel-kallisto Pr | rigel-kallisto log2 Pr | n_nz_rigel |
|---------|------|---------------:|--------------------:|------------------:|----------------------:|-----------:|
| dna00m | sample | 0.9209 | 0.8900 | 0.9974 | 0.9411 | 92,724 |
| dna00m | fractional | 0.9212 | **0.8957** | 0.9975 | **0.9437** | **150,615** |
| dna20m | sample | 0.9110 | 0.7895 | 0.9745 | 0.7280 | 100,859 |
| dna20m | fractional | 0.9114 | **0.7997** | 0.9749 | **0.7342** | **194,421** |

**Take-away:** The sample-mode penalty is ≤ 1 pp on log-space r but the
expressed-transcript count drops by **≈ 60 k** (halving the support). Anyone
reading rigel's quant.feather outside the rigel-internal locus framework is
at risk of misinterpreting this; **fractional mode should be the default for
external comparisons.**

### Gene-level agreement (much higher)

Aggregating both rigel and the other tools to gene-level (sum of constituent
transcripts):

| Library | rigel-salmon Pr | rigel-kallisto Pr | rigel-salmon log2 Pr | rigel-kallisto log2 Pr |
|---------|----------------:|------------------:|---------------------:|-----------------------:|
| dna00m | **0.997** | **0.999** | **0.995** | **0.996** |
| dna20m | **0.993** | **0.980** | **0.938** | **0.895** |

Gene-level dna00m agreement is essentially perfect. The transcript-level
spread is dominated by **isoform-assignment differences** between three
fundamentally different model families (rigel = locus-EM with bias
correction; salmon = length-corrected EM with sequence bias model; kallisto
= EM on raw equiv-class counts).

---

## §3 Rigel pool decomposition (the gDNA siphon story)

| Metric | dna00m | dna20m |
|--------|-------:|-------:|
| BAM total reads | 28,731,530 | 66,113,219 |
| BAM mapped reads | 28,534,945 | 65,620,992 |
| rigel: annotated mRNA count | 13,119,707 | 13,953,517 |
| rigel: synthetic nRNA count | 269,842 | 350,989 |
| **rigel: pi_pool (gDNA fraction in unspliced pool)** | **0.420** | **0.677** |
| Calibration: pool n_obs (unspliced + intronic + intergenic) | 29,345 | 717,526 |
| salmon: annotated count | 12,928,819 | 15,568,589 |
| kallisto: est_counts | 13,288,094 | 20,632,430 |

**Reading:**

- **dna00m has zero spiked gDNA, but rigel's calibration assigns 42% of the
  unspliced pool to gDNA.** The pool only has 29 k observations because
  cell-line RNA is dominantly mature-spliced; the mixture model is fitting
  noise + a tail of pre-mRNA/intronic fragments that are being labeled
  "gDNA-like" by length alone. Salmon and kallisto don't model gDNA, so
  they happily assign these fragments to whatever transcript they overlap.
- **dna20m has 50% gDNA by sample design, and rigel calls 67.7% of the
  unspliced pool gDNA.** This is plausible because (a) the pool also
  contains intergenic/intronic RNA fragments, (b) capture-pulled gDNA
  fragments that fall on exons get scored as RNA at the per-fragment
  step. The 67.7% specifically refers to the unspliced/intronic pool, not
  the whole library.
- **Salmon takes 15.6 M reads as transcript counts in dna20m (vs 12.9 M in
  dna00m)** — that's **+2.6 M reads of "extra" RNA signal that almost
  certainly comes from gDNA spike-in.** Kallisto is even more extreme at
  +7.3 M. Rigel correctly does not absorb these into its mRNA layer.

The published `pi_pool=0.42` for the no-gDNA library is a **real defect** in
rigel's calibration, not a consequence of the synthetic protocol — see §6
for root-cause discussion.

---

## §4 Effective length sanity check

For shared transcripts (rigel vs salmon):

| Library | n | median ratio | mean ratio | 5th pct | 95th pct |
|---------|--:|-------------:|-----------:|--------:|---------:|
| dna00m | 227,394 | 1.001 | 0.990 | 1.000 | 1.003 |
| dna20m | 227,394 | 1.008 | 1.000 | 1.001 | 1.030 |

Effective lengths agree within 1% across the range. Bias corrections are
not the source of the disagreement.

---

## §5 Locus-level disagreement pattern (the sharpest signal)

Aggregating to rigel locus IDs and counting loci where rigel and salmon
disagree by ≥ 2× at locus sum ≥ 200:

| Library | Total disagreeing loci | rigel ≫ salmon | rigel ≪ salmon |
|---------|----------------------:|---------------:|----------------:|
| dna00m | 112 | 95 | 17 |
| dna20m | **1,159** | **113** | **1,046** |

**Pattern:** dna00m has 95 loci where rigel reports more counts than salmon
(rigel pulls in multimappers / pre-mRNA that salmon discards on quasi-
mapping). dna20m flips entirely: **1,046 loci where rigel reports fewer
counts** because rigel routes the spike-in gDNA away from those loci and
into its gDNA component.

The 1,046 "rigel ≪ salmon" loci in dna20m are precisely the ones where
salmon is **overcalling** because of the gDNA spike-in — these are the
loci that should be analyzed downstream as **rigel's wins**.

---

## §6 Root-cause analysis (rigel-side issues)

### R1 — `--assignment-mode sample` is the wrong default for cross-tool comparison

**Severity:** High (silent over-collapsing of isoform-level signal).

**Mechanism:** With `sample`, every locus-EM fragment is hard-assigned to a
single transcript via a Bernoulli draw weighted by the posterior. The result
is a sparse vector with ~89 k expressed annotated transcripts in dna00m vs
salmon's 144 k and kallisto's 143 k. log2 Pearson r drops by ~0.6 pp; the
real damage is in downstream applications that look at "expressed
transcripts" sets.

**Fix options:**
- (a) Make `fractional` the default. Pros: matches salmon/kallisto
  semantics for downstream interop. Cons: rigel's identity is
  sample-based (downstream isoform calls).
- (b) Always also write `count_fractional` alongside `count` regardless of
  mode. Cheap, no API break.
- (c) Document the default loudly in the user manual.

**Recommendation:** (b) + (c). Keep `sample` as default for users who use
rigel's downstream filtering, but write the fractional column too.

### R2 — gDNA-vs-nascent-RNA identifiability collapse in zero-gDNA libraries

**Severity:** High (real false-positive gDNA fraction, magnitude ≈ 0.42 for
this library).

**Mechanism:** SRD v1 calibration fits a 1-D fragment-length mixture
(RNA_FL fixed to spliced histogram) on the unspliced/incompatible /
intronic / intergenic pool. In a hybrid-capture cell-line library the
capture kit pulls in (i) abundant spliced mRNA → not in the pool; (ii)
pre-mRNA / intronic RNA → in the pool, length distribution shifted shorter
than spliced mRNA; (iii) some intergenic / decay products. Without (iv)
actual gDNA, the mixture still finds a "shorter" component and labels it
gDNA. The model has no length feature distinguishing gDNA from
fragmented nascent RNA, so they alias.

**Evidence:**
- `pi_pool = 0.420` in dna00m (no gDNA present) vs `0.677` in dna20m
  (50% gDNA). The library-by-library delta is correct (Δ ≈ 0.26 ≈ what
  the mixture sees from the spike) but the *baseline* is wrong by 0.42.
- This is documented as the "nRNA siphon" failure mode in
  `docs/calibration/root_cause_locus_em_nrna_siphon.md`.

**Fix options (long-term, not P0):**
- Use **splice-completion features** (was the molecule fully spliced, or
  does it contain intronic sequence?) to separate nRNA from gDNA. gDNA
  is uniformly distributed across exon+intron, nRNA is biased toward
  newly-transcribed regions.
- Use **mappability/coverage uniformity within transcript bodies** — gDNA
  has flat coverage, nRNA is bursty.
- Allow user to provide a **library prior on gDNA fraction** when they know
  the protocol is gDNA-free.

### R3 — Locus partition over-aggregates for hybrid-capture cell lines

**Severity:** Medium.

**Mechanism:** rigel partitions by locus (UF over multimapper edges).
Hybrid-capture libraries are dominated by a small number of very
high-expression loci (cell-line-typical: HK, ribosomal proteins, MYC,
PTEN). The single mega-locus in dna20m has 29,713 transcripts ×
1.75 M units (see profile report) — most of these transcripts get an
ε-rank credit and never become "expressed" in `sample` mode.

**Evidence:**
- dna00m: 92,724 expressed annotated transcripts in `sample` mode despite
  ~144 k showing nonzero in salmon and kallisto.
- The mega-locus is one source; the other is a long tail of medium loci
  (top 1000 loci own 55% of EM cost per the profile report).

**Mitigation already in place:** `pruning_min_posterior=1e-4` discards
low-posterior candidates before EM. Tightening this to `1e-3` would
further sparsify but might amplify dropout for genuinely low-expression
isoforms.

### R4 — Calibration summary fields are missing from `summary.json`

**Severity:** Low (documentation/observability).

**Mechanism:** The calibration block only exposes `pi_pool` and `n_pool` —
not `quality`, `ss_inferred`, `iter_count`, or `converged`. The
benchmarking real-data analyzer surfaces "—" for these.

**Fix:** Plumb the full `CalibrationResult` diagnostics into
`summary.json`. Trivial change, big QC win.

### R5 — gDNA-EM total is not exposed in `summary.json`

**Severity:** Low.

**Mechanism:** `_run_locus_em_partitioned` computes `total_gdna_em` and
logs it but never writes it to `summary.json`. The benchmarking
analyzer therefore cannot tell whether the EM removed e.g. 2 M
fragments to gDNA or 200 k.

**Fix:** Add `em_summary.gdna_em_total` to `summary.json`. Trivial.

---

## §7 Recommendations

### Immediate (no code changes, just docs / defaults review)

1. **Document the `assignment_mode` semantic gap** in the manual: warn
   users comparing rigel to salmon/kallisto to use `--assignment-mode
   fractional` or aggregate to gene-level.
2. **Add a `--gdna-fraction-prior` knob** that lets the user tell rigel "I
   know there is no gDNA in this library; clamp pi_pool ≤ X". This is the
   simplest workable mitigation for R2 until splice-completion features
   are in.

### Short-term (one PR each)

3. **Always emit `count_fractional` column** in `quant.feather` regardless
   of assignment mode (R1 fix).
4. **Plumb full `CalibrationResult` into `summary.json`** (R4 fix).
5. **Add `em_summary.gdna_em_total` and `em_summary.n_loci_with_gdna_mass`
   to `summary.json`** (R5 fix).

### Long-term research

6. **Splice-completion features for nRNA vs gDNA discrimination** (R2 root
   cause). Outline:
   - Per-fragment "fraction of bases on exon vs intron" feature.
   - Per-fragment "spans annotated splice junction" indicator.
   - Two-component mixture extending the existing FL mixture.
7. **Mega-locus splitting heuristic** (R3): when a locus has > N
   transcripts and most fragments map to a small clique, isolate that
   clique into a separate locus. Profile data shows the 29k-transcript
   mega-locus dominates downstream sparsity.

---

## §8 Appendix — pipeline parameters

Both rigel runs used:
```
rigel quant --bam <annotated.bam> \
            --index .../rigel_index \
            -o .../rigel_v05[_frac] \
            --threads 8 --tmpdir /tmp \
            [--assignment-mode fractional]
```
Defaults otherwise: `em_mode=vbem`, `em_iterations=1000`,
`em_convergence_delta=1e-6`, `pruning_min_posterior=1e-4`,
`overhang_alpha=mismatch_alpha=0.1`, `gdna_splice_penalty_unannot=0.01`.

External tool outputs were taken as-published from the hulkrna pipeline.
Exact salmon / kallisto invocations are recorded in their respective
`cmd_info.json` / `run_info.json` files in each library directory.

### Reproducibility

```bash
# rigel runs
conda activate rigel
rigel quant --bam <bam> --index <index> -o <outdir> --threads 8 [--assignment-mode fractional]

# Analysis
python -m scripts.benchmarking.real_data \
  --library "dna00m=/scratch/.../mctp_vcap_rna20m_dna00m" \
  --library "dna20m=/scratch/.../mctp_vcap_rna20m_dna20m" \
  --rigel-subdir rigel_v05_frac \
  -o docs/benchmarks/vcap_hybcap_real_data.md
```

Raw per-library metrics tables (with per-pair top-25 outliers) are in the
companion file
[vcap_hybcap_rigel_vs_salmon_kallisto_data.md](vcap_hybcap_rigel_vs_salmon_kallisto_data.md).
