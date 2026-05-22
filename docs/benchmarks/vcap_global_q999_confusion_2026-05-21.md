# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-21

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/global_q999_v1/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 48.10% gDNA and 50.60% RNA across all counted fragments, with 1.31% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 48.74% gDNA and 51.26% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 95.83% and gDNA at 3.28%. gDNA-source fragments are called gDNA at 82.49% and RNA at 15.88%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 95.83% | 3.28% | 0.89% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 15.88% | 82.49% | 1.63% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,406,465 (95.83%) | 459,112 (3.28%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 2,894,671 (15.88%) | 15,037,544 (82.49%) | 296,462 (1.63%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,860,850 (91.93%) | 545,615 (3.90%) | 459,112 (3.28%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 1,085,039 (5.95%) | 1,809,632 (9.93%) | 15,037,544 (82.49%) | 296,462 (1.63%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 1,224,147 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 848,745 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 462,095 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 373,766 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 201,730 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 145,027 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 63,104 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_implicit | 42,163 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 34,686 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 19,914 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 11,966 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 11,703 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_annot | 11,084 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 9,574 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,491 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_annot | 9,159 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | spliced_implicit | 6,867 |
| H7MFFDSXY | gdna | rna | mrna | multimapper | unspliced | 6,387 |

## Artifacts

- Flowcell summary: `results/vcap_global_q999_confusion_2026-05-21/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_global_q999_confusion_2026-05-21/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_global_q999_confusion_2026-05-21/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_global_q999_confusion_2026-05-21/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_global_q999_confusion_2026-05-21/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_global_q999_confusion_2026-05-21/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_global_q999_confusion_2026-05-21/summary.json`
