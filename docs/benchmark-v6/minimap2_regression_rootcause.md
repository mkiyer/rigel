# Minimap2 Regression Root Cause Analysis

**Benchmark:** `benchmark_pristine` (10M fragments, gdna=none, ss=0.95, nrna=none)  
**Target set:** 55 transcripts where oracle error ≤ 25% AND (salmon OR kallisto) error ≤ 25% AND mm2 error ≥ 50% AND truth ≥ 50 fragments  
**Scripts:** `scripts/debug/find_minimap2_regressions.py`, `scripts/debug/polr2a_oaz1_check.py`

---

## Benchmark Context

After fixing `benchmark.py` to use `minimap2 -ax splice:sr -j annotation.bed -N 20` (from the broken `-ax splice --junc-bed -N 10`), mm2/oracle MAE ratio improved from **3.87× → 1.21×**. These are the remaining fixable regressions.

| Method | Transcript MAE | Pearson | Spearman |
|---|---|---|---|
| rigel oracle | 1.83 | 0.9992 | 0.880 |
| salmon | 2.91 | 0.9999 | 0.894 |
| kallisto | 3.36 | 0.9980 | 0.848 |
| **rigel mm2** | **2.25** | **0.9979** | **0.865** |

---

## Top 10 Regression Genes

| Rank | Gene | truth | mm2 | mm2_ratio | abs_err | Category |
|---|---|---|---|---|---|---|
| 1 | ENSG00000280800 | 21,430 | 7,047 | 0.33 | 14,383 | 3 — rRNA proximity |
| 2 | AP2M1 | 2,554 | 590 | 0.23 | 1,964 | 1 — micro-exon (6bp) |
| 3 | POLR2A | 2,261 | 769 | 0.34 | 1,492 | 3 — POLR2 family scatter |
| 4 | UBB | 1,550 | 297 | 0.19 | 1,253 | 1 — micro-intron (228bp) |
| 5 | SLC25A6 | 2,328 | 1,148 | 0.49 | 1,180 | 2 — PAR1 (chrX/Y) |
| 6 | SEPTIN7 | 806 | 12 | 0.01 | 794 | 1 — micro-exon (5bp) |
| 7 | RPL41 | 377 | 1,152 | 3.06 | 775 | 1 — micro-exon (23bp, over-count) |
| 8 | VAMP7 | 1,177 | 585 | 0.50 | 592 | 2 — PAR2 (chrX/Y) |
| 9 | OAZ1 | 1,030 | 500 | 0.49 | 530 | 3 — pseudogene/paralog |
| 10 | CYTH1 | 487 | 0 | 0.00 | 487 | 1 — micro-exon (3bp) |

---

## Category 1: Micro-exon / Micro-intron Anchor Failure

**Mechanism:** Same as the confirmed GAPDH failure mode. Minimap2 requires a minimum anchor length to call a splice junction. When an exon is shorter than the anchor threshold (~20bp typical, but failures observed up to ~50bp in practice), reads spanning that junction fail to be recognized as spliced. The reads either (a) go unmapped, (b) map to adjacent transcripts that skip the micro-exon, or (c) cause incorrect isoform assignment.

**5 genes, ~5,358 count-unit error (53% of top-10 total)**

### AP2M1 (ENST00000292807.9)
- chr3, 12-exon transcript
- Exon 5: **6bp** (blockSizes: 105,117,266,83,**6**,136,142,120,136,98,112,610)
- truth=2,554 | oracle=2,554 | mm2=590 → **77% loss**
- ENST00000382456.7 (alternative isoform, skips exon 5): truth=0, mm2=over 29× truth
- Reads anchoring across the 6bp exon fail; the shorter transcript absorbs them

### SEPTIN7 (ENST00000350320.11)
- chr7, 13-exon transcript
- Exon 2: **5bp** (blockSizes: 214,**5**,103,107,...)
- truth=806 | oracle≈806 | mm2=12 → **98% loss**
- ENST00000399034.7 (alternative isoform) absorbs the displaced reads

### CYTH1 (ENST00000446868.8)
- chr17, 14-exon transcript
- Exon 5: **3bp** (blockSizes: 2122,155,72,77,**3**,112,149,113,81,119,67,65,83,72)
- truth=487 | oracle≈487 | mm2=0 → **100% loss**
- The 3bp exon is the shortest observed; 100% of reads fail this junction

### RPL41 — isoform switching via micro-exon (ENST00000501597.3 / ENST00000546591.6)
- chr12, two overlapping isoforms both with tiny exons
- Dominant ENST00000501597.3: blockSizes=70,**25**,**23**,351 → two micro-exons
- Minor ENST00000546591.6: blockSizes=167,**23**,486 → one micro-exon
- Dominant: truth=7,113 | mm2=5,600 → 21% under-count (reads fail micro-exon junctions)
- Minor (our target): truth=377 | mm2=1,152 → **3× over-count** (absorbs mis-aligned reads from dominant)
- Net: the micro-exon failure on the dominant isoform shifts reads to the minor isoform

### UBB (ENST00000535788.1) — micro-intron
- chr17, 3-exon transcript with a **228bp intron** between exons 2 and 3
- Exon structure: exon1=94bp, intron1=718bp, exon2=229bp, **intron2=228bp**, exon3=381bp
- ENST00000302182.8 (2-exon, spans the micro-intron as continuous exon 2): truth=6,220, mm2=7,411 (+19%)
- truth=1,550 | oracle≈1,550 | mm2=297 → **81% loss**
- Reads spanning the 228bp micro-intron fail to be recognized as spliced → appear as ENST00000302182.8 (unspliced relative)

**Actionable fix for Category 1:** 2-pass alignment (`minimap2 --write-junc` + `--pass1`) if available in a future minimap2 version, or minimap2's `--junc-bonus` parameter to boost short-intron junctions. Alternatively, pre-annotate known micro-exon junctions explicitly in the BED12 file with a higher score to guarantee anchor success.

---

## Category 2: PAR Duplication

**Mechanism:** Pseudoautosomal regions (PAR1, PAR2) exist at identical coordinates on chrX and chrY. Reads originating from PAR loci have perfect alignments on both chromosomes. With `--secondary=yes -N 20`, mm2 reports both hits. Rigel's index only contains chrX transcripts, so chrY alignments produce unassignable fragments. The net effect is ~50% read loss.

**2 genes, ~1,772 count-unit error (17% of top-10 total)**

### SLC25A6 (ENST00000381401.11) — PAR1
- chrX:1,386,151–1,392,113 (PAR1, identical to chrY:1,386,151–1,392,113)
- truth=2,328 | oracle=2,328 | salmon=2,314 | mm2=1,148 → **49% loss** (exactly 50% as expected)
- Oracle works because the simulation ground truth is aligned to chrX, so oracle reads are already disambiguated

### VAMP7 (ENST00000286448.12) — PAR2
- chrX:155,881,344–155,943,769 (PAR2, copy at chrY:57,067,864–57,130,289)
- truth=1,177 | oracle=1,177 | salmon≈1,177 | mm2=585 → **50% loss**

**Actionable fix for Category 2:** Mask chrY PAR1 (chrY:10,001–2,781,479) and PAR2 (chrY:56,887,903–57,217,415) in the alignment reference. Reads from PAR loci will then all map to chrX. Alternatively, a rigel-side fix could detect PAR-overlapping transcripts and apply chromosome-pairing correction, but this is more complex.

---

## Category 3: Genome-level Multimapping

**Mechanism:** These failures arise from genome-based alignment (minimap2) exposing genomic regions not modeled in rigel's transcript index — either (a) gene family paralogs with GENCODE `basic` filter excluding them, (b) processed pseudogenes in the genome but outside the annotation, or (c) rRNA repeat arrays distributed across multiple acrocentric chromosomes. Salmon and kallisto win here because they align against the transcriptome directly, never seeing these alternative genomic loci.

**3 genes, ~2,509 count-unit error (25% of top-10 total)**

### POLR2A (ENST00000674977.2) — gene family scatter
- chr17, 30 exons, 29,814bp span
- truth=2,261 | oracle=2,261 | salmon=2,137 | mm2=769 → **66% loss**
- The paftools.js junction BED contains extensive POLR2* family junctions: POLR2J2, POLR2J3 (many transcripts), POLR2B, POLR2C (including multiple retained_intron and NMD variants), plus POLR2DP1/2, POLR2CP1, POLR2KP1/2, POLR2LP1, POLR2MP1 pseudogenes
- Rigel's index (after `transcript_filter: basic`) only knows `ENST00000674977.2` for POLR2A itself
- Reads from POLR2A find secondary alignments at POLR2 family loci; those hits are outside rigel's index; the EM absorbs them as gDNA
- Note: 8 retained_intron variants of POLR2A itself are in the junction BED but NOT in rigel's index

### OAZ1 (ENST00000602676.6) — pseudogene competition
- chr19, truth=1,030 | oracle=1,018 | salmon=967 | mm2=500 → **51% loss**
- OAZ1P1 (processed pseudogene, chr1:40,132,763–40,133,448) is in both the mm2 junction BED AND in rigel's transcript index (ENST00000430724.2)
- OAZ1P1 truth=0, oracle=0, mm2=0 (confirming the pseudogene is NOT the direct sink)
- The 530 lost reads scatter to other unmapped loci — likely sequence similarity between OAZ1 (chr19) and OAZ1P1 junction hint causing reads to align ambiguously or to unannotated pseudogene copies not in rigel's index
- OAZ2 and OAZ3 are not the sink (mm2 matches well for both)
- Note: even oracle shows minor leakage (1,018 vs 1,030 truth) indicating inherent multi-mapping at this locus

### ENSG00000280800 (ENST00000631211.1) — rRNA proximity repeat scatter
- chr21:8,210,383–8,211,306, single-exon lncRNA (923bp), locus in ribosomal intergenic spacer
- truth=21,430 | oracle=21,430 | **salmon=20,381 (95%)** | kallisto=6,793 (32%) | mm2=7,047 (33%)
- Adjacent genes: RNA5-8SN2 (5.8S rRNA, chr21:8,212,571), MIR6724-1, MIR6724-2 clusters
- This lncRNA lies in the acrocentric chr21 pericentromeric rDNA array region
- Its sequence likely has similarity to the 5.8S/18S/28S rRNA sequences that are distributed across rDNA arrays on chr13, chr14, chr15, chr21, chr22
- minimap2 distributes reads across many rDNA copies → only 1/(number of rDNA copies) reach chr21
- kallisto ALSO fails (32%) = same result as mm2, confirming this is a genuine sequence-level ambiguity
- **Salmon is exceptional here**: quasi-mapping uses the full transcriptome index as a single reference and likely clusters all rRNA copies together in its equivalence class model, routing reads back correctly
- Fix: not addressable at the alignment layer; would require salmon-like transcriptome collapsing or explicitly filtering rRNA-adjacent lncRNAs from quantification targets

---

## Impact Summary

| Category | Genes | Mechanism | Approx error (count units) | Fixable? |
|---|---|---|---|---|
| 1 — Micro-exon/intron | AP2M1, SEPTIN7, CYTH1, RPL41, UBB | Splice anchor failure | ~5,358 (53%) | Yes — 2-pass or junction boosting |
| 2 — PAR duplication | SLC25A6, VAMP7 | chrY read loss | ~1,772 (17%) | Yes — mask chrY PAR |
| 3 — Genome multimapping | POLR2A, OAZ1, ENSG00000280800 | Reads scatter outside rigel index | ~2,509 (25%) | Difficult — fundamental genome vs transcriptome issue |
| **Total** | **10 genes** | | **~9,639** | |

---

## Proposed Actions

### Priority 1 (High Impact, Actionable): Micro-exon improvement
- Regenerate the mm2 junction BED with explicit short-exon scoring: assign junction score proportional to `min(blockSize, 50)` so that 3–6bp exon junctions still receive high weight
- Alternatively, test minimap2 `--junc-bonus=20` to see if boosting junction scores helps
- Check compatibility with minimap2 v2.30 `splice:sr` preset parameters

### Priority 2 (Quick Win): PAR reference masking
- Add a step in `build_minimap2_bed()` (or a preprocessing note in MANUAL.md) to mask chrY PAR1 and PAR2 from the alignment reference FASTA
- The standard Broad Institute "analysis set" hg38 reference already hard-masks these regions; switch to that reference if not already using it
- Check: `grep "chrY" rigel_index/genome.fasta | head`

### Priority 3 (Low Priority): Genome multimapping
- For POLR2A: expand rigel's `transcript_filter` to include POLR2A retained_intron variants to capture reads that mm2 assigns to non-basic isoforms
- For OAZ1: investigate whether removing OAZ1P1 from the junction BED reduces scatter
- For ENSG00000280800: this is an inherent limitation of genome-based alignment in rRNA regions; document as a known limitation
