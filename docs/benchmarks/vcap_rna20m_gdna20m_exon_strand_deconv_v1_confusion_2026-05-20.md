# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-20

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 49.35% gDNA and 49.35% RNA across all counted fragments, with 1.31% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 50.00% gDNA and 50.00% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 95.77% and gDNA at 3.34%. gDNA-source fragments are called gDNA at 84.66% and RNA at 13.72%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 95.77% | 3.34% | 0.89% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 13.72% | 84.66% | 1.63% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,398,759 (95.77%) | 466,818 (3.34%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 2,500,290 (13.72%) | 15,431,925 (84.66%) | 296,462 (1.63%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,910,790 (92.29%) | 487,969 (3.49%) | 466,818 (3.34%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 1,171,473 (6.43%) | 1,328,817 (7.29%) | 15,431,925 (84.66%) | 296,462 (1.63%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 918,169 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 854,615 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 378,530 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 361,080 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 217,716 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 145,027 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 64,449 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_implicit | 41,919 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 30,380 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 14,926 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 11,966 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 11,703 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_annot | 10,662 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 10,644 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,491 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_annot | 9,159 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | spliced_implicit | 6,821 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_implicit | 6,600 |

## Artifacts

- Flowcell summary: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion/summary.json`
