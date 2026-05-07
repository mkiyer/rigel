# Synthetic Sim Suite — Re-run & Deep RCA

## Headline finding (corrected)

**Per-gene RNA quantification is essentially perfect.** Median per-gene relative error ≈ 0.0005; total absolute gene error ≈ 4–5K out of 1M fragments (<0.5%) across all conditions, including high-gDNA + low strand-specificity.

Two methodological gotchas in earlier reports:

1. The `truth_abundances.tsv` `mrna_abundance` column is in **simulated molecules** (sum ≈ 112K), not fragments (sum = 1M). All the "MARD ~39 / gdna_leak ~22K" headline numbers in the existing benchmark report were computed against an incorrectly-scaled truth. Using fragment-level oracle truth from the BAM read-name (`{tid}:{frag_id}`), accuracy is dramatically better than reported.
2. Paired-end BAMs have 2 records per fragment — must filter to read1 when counting.

Once corrected, the picture is clean enough to enumerate a small number of distinct issues in priority order.

---

## P0 — Exact-duplicate isoforms in the GTF cause unresolvable EM splits

### Evidence
Stratified per-transcript accuracy (gdna_none_ss_0.99):

| Stratum | n | Σtruth | Σpred | abs_err | RE |
|---|---:|---:|---:|---:|---:|
| Structurally unique exon set | 95 | 486 925 | 485 758 | 9 731 | **2.0 %** |
| Near-dup (same intron chain, UTR-only diff) | 151 | 513 074 | 514 241 | 81 193 | **15.8 %** |
| Exact-duplicate exon set | 28 (in 14 groups) | 43 856 | 44 091 | 41 817 | **95.4 %** |

### Root cause
The synthetic GTF contains 14 transcript groups (28 transcripts) whose exon coordinates are bit-identical:

- `GENE0013.1 ≡ GENE0013.4` (all 11 exons identical)
- `GENE0048.1 ≡ GENE0048.2` (all 12 exons identical)
- `GENE0001.1 ≡ GENE0001.9`, `GENE0005.1 ≡ GENE0005.2`, `GENE0008.1 ≡ GENE0008.5 ≡ GENE0008.7`, `GENE0010.1 ≡ GENE0010.2`, etc.

For these pairs there is **no fragment that can carry information distinguishing the isoforms** — they are mathematically unidentifiable. Annotated-BAM tracing confirms `count_unambig = 0` for every FP transcript and EM correctly converges to a 50/50 split (`ZW ≈ 0.49`). The salmon-style L̃<sub>t</sub> fix from this session does not (and cannot) help because L̃ values are equal for identical structures.

A second, quieter cohort (151 transcripts) shares an intron chain but differs in UTR start/end by tens of bases. These also degrade (15.8% RE) because only fragments overlapping the differential UTR window can break the tie.

### Research areas (priority-ordered)
1. **Index-time isoform collapse / equivalence-class output** — detect identical exon signatures during `rigel index`; either (a) collapse duplicates to a single canonical representative and emit a `transcript_id_alias` column in `quant.feather`, or (b) keep both and report posterior at the equivalence-class level with a per-class uncertainty annotation. This is the single most impactful change.
2. **Sparsity/concentration prior in EM** — a Dirichlet `α < 1` (or an L1 penalty) for transcripts whose distinguishing-fragment count is below threshold pushes mass onto the dominant isoform. Need to verify it does not break the `Issue 5` calibration.
3. **UTR-tolerance equivalence classes** — collapse near-duplicates whose intron chain matches and whose UTR boundary differs by < N bases (configurable, e.g. 50 bp), unless distinguishing fragments are observed.
4. **Per-class identifiability flag in summary.json** — surface the unresolvable-pairs explicitly so users do not mis-attribute zero-vs-low expression calls.

---

## P1 — gDNA → RNA leakage at high contamination

### Evidence
At 67% gDNA fraction (gdna_high), of the 2M oracle gDNA fragments that reach the locus pipeline, 23 K (1.1%) get assigned a gene tag. Concentrated in the few longest, highest-expressed loci (`GENE0019`, `GENE0001`, `GENE0040`, `GENE0052`). Persists at both ss=0.99 and ss=0.50.

### Root cause (hypothesis, partially confirmed)
gDNA fragments falling on long single-exon stretches (or large terminal exons) of the brightest genes look indistinguishable from unspliced exonic RNA. The fragment-length signal (gDNA mean ≈ 350 vs RNA mean ≈ 250) is informative on average but per-fragment overlap of the two FL distributions allows ~1% leak. The locus-level prior assigns most exonic mass to RNA when RNA dominates, which is the right Bayesian behavior but the wrong frequentist behavior.

### Research areas
1. **Length-aware gDNA prior per gene** — currently the locus-level π_gdna applies uniformly across all candidates in the locus. A per-position gDNA density (using the genome-wide intergenic ρ_ig estimate) gives a more informative prior on long single-exon stretches.
2. **Stricter splice-junction-bearing weighting** — fragments that have a junction can be near-deterministically classified RNA; only unspliced fragments should ever participate in the gDNA-vs-RNA call. Verify this is already enforced.
3. **FL likelihood ratio cap** — diagnostic to confirm gDNA-typed fragments at high abundance are dominated by unspliced reads with FL ≈ 350.

---

## P2 — nRNA component absorbs gDNA at high contamination

### Evidence
`nrna_none` truth (zero nRNA) → 11–15 K fragments assigned to nRNA in `gdna_high_*` conditions. Effect scales with gDNA fraction, not strand specificity.

### Root cause
nRNA components share intronic span with gDNA. With heavy gDNA contamination, intronic gDNA fragments fit the nRNA likelihood (long unspliced intronic reads) almost as well as the gDNA likelihood. EM splits the ambiguous mass.

### Research areas
1. **Sparser nRNA prior when no spliced-RNA evidence** — gate nRNA component activation on observation of any junction-bearing fragment that overlaps the nRNA span.
2. **gDNA-vs-nRNA discriminator** using FL: nRNA tends to be processed-length distributed while genomic gDNA is shotgun-fragmented. The `_fl_mixture.py` calibration should already separate these — verify it actually does for the high-gDNA condition.

---

## P3 — Exonic-vs-intergenic gDNA density slightly under-allocates exons

### Evidence
Calibration `analyze_deep` reports `ρ_ex / ρ_ig ≈ 0.86–0.91` across conditions; theoretical expectation ≈ 1.0. Under-allocation of gDNA mass to exonic positions inflates the residual RNA estimate slightly.

### Root cause (hypothesis)
Either (a) the exonic effective length under-counts mappable positions vs. intergenic — possibly because the FL-PMF containment correction subtracts edge effects more aggressively in short exons, or (b) a mappability differential between exonic and intergenic regions in the synthetic genome.

### Research areas
- Add a unit test that the SRD calibration recovers ρ_ex/ρ_ig ≈ 1.0 on a uniform-coverage gDNA-only synthetic.
- Cross-check exonic effective-length computation against a brute-force Monte-Carlo on a few exemplar exons.

---

## P4 — Wrong-gene assignment at SS=0.50

### Evidence
At ss=0.50 (unstranded), antisense-overlapping or convergent loci show ~0.06–0.10% wrong-gene assignments. Magnitude small (~600–1000 reads at 1M scale).

### Root cause
With unstranded reads the strand-LLR component contributes zero discrimination for antisense-overlapping pairs (e.g. `Gene48` and `Gene49as` share locus 34). The fragment-length prior plus exon-coverage break the tie correctly most of the time, but ~0.1% slips.

### Status
Acceptable; not worth optimizing in isolation. Likely will ride along with any P0/P1 work.

---

## Summary of recommendations

| Priority | Issue | Suggested next step |
|---|---|---|
| **P0** | Exact-duplicate isoforms unidentifiable | Implement equivalence-class collapse at `rigel index` time + reporting at class level. Largest win by far. |
| P1 | gDNA→RNA leakage at high contamination | Per-position gDNA prior on long exons; verify junction-only RNA-vs-gDNA gating. |
| P2 | gDNA→nRNA siphon | Gate nRNA activation on observed splice evidence. |
| P3 | ρ_ex/ρ_ig ≈ 0.88 (slight) | Audit exonic effective-length math; add unit test. |
| P4 | Wrong-gene at ss=0.50 | Defer; minor and likely improves with P1/P2. |

The earlier "nRNA siphon" / "GENE0013.1 phantom" stories are now mostly explained by P0. Once duplicate-isoform handling is in place, the residual error budget will be small enough that P1–P3 become individually visible and easier to attack.

Diagnostic scripts kept at `/tmp/{inspect_fp_origin,diff_isoforms,oracle_truth,oracle_truth_v2,stratified_accuracy}.py`. Findings saved to repo memory `synthetic-sim-rca-2026-05-07.md`.