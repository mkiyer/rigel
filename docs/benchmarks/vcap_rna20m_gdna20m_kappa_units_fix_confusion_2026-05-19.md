# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-19

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 49.81% gDNA and 48.89% RNA across all counted fragments, with 1.31% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 50.47% gDNA and 49.53% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 95.48% and gDNA at 3.64%. gDNA-source fragments are called gDNA at 85.24% and RNA at 13.13%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 95.48% | 3.64% | 0.89% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 13.13% | 85.24% | 1.63% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,356,978 (95.48%) | 508,599 (3.64%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 2,393,558 (13.13%) | 15,538,657 (85.24%) | 296,462 (1.63%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,860,866 (91.93%) | 496,112 (3.55%) | 508,599 (3.64%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 1,062,195 (5.83%) | 1,331,363 (7.30%) | 15,538,657 (85.24%) | 296,462 (1.63%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 877,209 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 830,287 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 412,914 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 342,908 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 197,129 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 145,027 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 71,368 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_implicit | 42,010 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 28,111 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 15,016 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 11,966 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 11,703 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 11,117 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_annot | 10,799 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,491 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_annot | 9,159 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | spliced_implicit | 6,821 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_implicit | 6,509 |

## Artifacts

- Flowcell summary: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_rna20m_gdna20m_kappa_units_fix_confusion_2026-05-19/summary.json`
