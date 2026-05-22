# FLG2 Hotspot EM Failure Diagnosis - 2026-05-21

Window: `chr1:152,350,001-152,360,000`

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1/annotated.bam`

The hotspot is solved inside MultiLocus/locus `3`, a chromosome-scale component with 42,639 transcript components, 22,701 nRNA entities, 2,823,632 EM fragments, and a gDNA effective length of 92,407,095 bp. The same component contains 1,446,008 assigned gDNA fragments, but its gDNA mass is spread over a very large effective opportunity.

## Local Confusion

Fragments overlapping the 10 kb window: 4,699. True RNA-source fragments in this window are only 23, while true gDNA-source fragments are 4,676. Rigel assigns 3,357 true-gDNA fragments to RNA (71.79%).

| Truth | Pred pool | Detail | Count |
| --- | --- | --- | --- |
| gdna | rna | nrna | 1737 |
| gdna | rna | mrna | 1620 |
| gdna | gdna | gdna | 1307 |
| gdna | unresolved | unresolved | 12 |
| rna | rna | nrna | 8 |
| rna | rna | mrna | 7 |
| rna | unresolved | unresolved | 7 |
| rna | gdna | gdna | 1 |

## Local Strand Split

`read1_strand` is the genomic strand of R1. Since this run inferred an R1-antisense protocol, `rna_strand_if_rna` is the transcript strand implied if the fragment were RNA. In the true-gDNA population, the two orientations are both present, which lets the EM explain the window as two RNA components: negative-strand FLG2 mRNA plus positive-strand antisense/nRNA.

| Truth | Pred | R1 | RNA strand if RNA | Count | Fraction |
| --- | --- | --- | --- | --- | --- |
| gdna | mrna | + | - | 1,600 | 98.77% |
| gdna | nrna | - | + | 1,517 | 87.33% |
| gdna | gdna | - | + | 749 | 57.31% |
| gdna | gdna | + | - | 558 | 42.69% |
| gdna | nrna | + | - | 220 | 12.67% |
| gdna | mrna | - | + | 20 | 1.23% |
| rna | nrna | - | + | 8 | 100.00% |
| gdna | unresolved | - | + | 7 | 58.33% |
| gdna | unresolved | + | - | 5 | 41.67% |
| rna | mrna | + | - | 5 | 71.43% |
| rna | unresolved | - | + | 5 | 71.43% |
| rna | unresolved | + | - | 2 | 28.57% |
| rna | mrna | - | + | 2 | 28.57% |
| rna | gdna | + | - | 1 | 100.00% |

## False-RNA Targets

| Pred | Target | Count | False-RNA frac |
| --- | --- | --- | --- |
| mrna | mRNA FLG2 ENST00000388718.5 | 1,581 | 47.10% |
| nrna | nRNA chr1:152,335,378-152,364,080(+) | 1,504 | 44.80% |
| nrna | nRNA chr1:152,348,734-152,360,006(-) | 212 | 6.32% |
| mrna | mRNA SMCP ENST00000368765.4 | 1 | 0.03% |
| mrna | mRNA CHRNB2 ENST00000637900.1 | 1 | 0.03% |
| nrna | nRNA chr1:1,702,880-1,721,655(-) | 1 | 0.03% |
| nrna | nRNA chr1:154,479,480-154,502,412(-) | 1 | 0.03% |
| mrna | mRNA F5 ENST00000367796.3 | 1 | 0.03% |
| mrna | mRNA AHDC1 ENST00000642416.1 | 1 | 0.03% |
| mrna | mRNA MROH9 ENST00000367758.7 | 1 | 0.03% |
| mrna | mRNA IPO9 ENST00000361565.9 | 1 | 0.03% |
| mrna | mRNA HMCN1 ENST00000271588.9 | 1 | 0.03% |

## False-RNA Categories

| Pred | ZC | ZS | R1 | RNA strand if RNA | Count | False-RNA frac |
| --- | --- | --- | --- | --- | --- | --- |
| mrna | ambig_opp_strand | unspliced | + | - | 1,552 | 46.23% |
| nrna | ambig_opp_strand | unspliced | - | + | 1,470 | 43.79% |
| nrna | ambig_opp_strand | unspliced | + | - | 209 | 6.23% |
| mrna | multimapper | unspliced | + | - | 18 | 0.54% |
| nrna | multimapper | unspliced | - | + | 15 | 0.45% |
| nrna | multimapper | spliced_unannot | - | + | 11 | 0.33% |
| mrna | ambig_same_strand | spliced_implicit | + | - | 11 | 0.33% |
| nrna | ambig_opp_strand | spliced_unannot | - | + | 10 | 0.30% |
| mrna | ambig_same_strand | spliced_implicit | - | + | 10 | 0.30% |
| mrna | ambig_opp_strand | spliced_unannot | + | - | 9 | 0.27% |
| mrna | ambig_opp_strand | spliced_implicit | + | - | 8 | 0.24% |
| mrna | ambig_opp_strand | spliced_implicit | - | + | 8 | 0.24% |
| nrna | ambig_same_strand | spliced_implicit | - | + | 6 | 0.18% |
| nrna | ambig_opp_strand | spliced_implicit | - | + | 5 | 0.15% |
| nrna | ambig_opp_strand | spliced_implicit | + | - | 4 | 0.12% |
| nrna | ambig_same_strand | spliced_implicit | + | - | 4 | 0.12% |

## Overlapping Transcript Models

The reference index presents FLG2, positive-strand annotated antisense transcripts, and synthetic positive/negative nRNA spans as simultaneous candidates in this window.

| Transcript | Gene | Strand | Kind | Exon bp in window | Span bp in window | Run count | TPM |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ENST00000388718.5 | FLG2 | - | annotated | 7,852 | 10,000 | 1,583 | 17.501 |
| ENST00000655109.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000630125.3 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000669062.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000664213.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000669830.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000653548.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000659844.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000666686.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000392688.7 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000658099.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000629331.1 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| ENST00000445097.2 | CCDST | + | annotated | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152189296_152365652 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152189296_152404507 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152189332_152365557 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152189332_152365652 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152189332_152405050 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152313458_152365652 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152331475_152365652 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_1_152335378_152364080 |  | + | synthetic nRNA | 0 | 10,000 | 1,513 | 5.014 |
| RIGEL_NRNA_chr1_1_152337067_152366692 |  | + | synthetic nRNA | 0 | 10,000 | 0 | 0.000 |
| RIGEL_NRNA_chr1_2_152348734_152360006 |  | - | synthetic nRNA | 0 | 10,000 | 212 | 1.813 |

## Calibration Regions In The Window

| Region | Type mask | Strand | tx+ | tx- | exon+ | exon- |
| --- | --- | --- | --- | --- | --- | --- |
| 152,348,735-152,357,647 | 2 | - | 8,913 | 8,913 | 0 | 8,913 |
| 152,357,648-152,358,746 | 1 | ? | 1,099 | 1,099 | 0 | 0 |
| 152,358,747-152,358,906 | 2 | - | 160 | 160 | 0 | 160 |
| 152,358,907-152,359,955 | 1 | ? | 1,049 | 1,049 | 0 | 0 |
| 152,359,956-152,360,006 | 2 | - | 51 | 51 | 0 | 51 |

## Nearby Antisense Evidence Check

This table compares the hotspot to the positive-strand antisense transcript span outside the hotspot. It is a read-level way to test the browser observation that the antisense transcript model has little support outside the FLG2-overlapping exon/intron region.

| Feature | Truth | Pred | Count | Within-truth frac |
| --- | --- | --- | --- | --- |
| antisense_exons_outside_window | gdna | nrna | 173 | 36.50% |
| antisense_exons_outside_window | gdna | mrna | 155 | 32.70% |
| antisense_exons_outside_window | gdna | gdna | 142 | 29.96% |
| antisense_exons_outside_window | gdna | unresolved | 4 | 0.84% |
| antisense_exons_outside_window | rna | nrna | 21 | 55.26% |
| antisense_exons_outside_window | rna | unresolved | 8 | 21.05% |
| antisense_exons_outside_window | rna | mrna | 7 | 18.42% |
| antisense_exons_outside_window | rna | gdna | 2 | 5.26% |
| antisense_introns_outside_window | gdna | nrna | 1,986 | 35.10% |
| antisense_introns_outside_window | gdna | gdna | 1,847 | 32.64% |
| antisense_introns_outside_window | gdna | mrna | 1,810 | 31.99% |
| antisense_introns_outside_window | gdna | unresolved | 15 | 0.27% |
| antisense_introns_outside_window | rna | nrna | 548 | 85.89% |
| antisense_introns_outside_window | rna | gdna | 68 | 10.66% |
| antisense_introns_outside_window | rna | unresolved | 14 | 2.19% |
| antisense_introns_outside_window | rna | mrna | 8 | 1.25% |
| flg2_exons | gdna | nrna | 1,736 | 37.13% |
| flg2_exons | gdna | mrna | 1,620 | 34.65% |
| flg2_exons | gdna | gdna | 1,307 | 27.96% |
| flg2_exons | gdna | unresolved | 12 | 0.26% |
| flg2_exons | rna | nrna | 8 | 34.78% |
| flg2_exons | rna | mrna | 7 | 30.43% |
| flg2_exons | rna | unresolved | 7 | 30.43% |
| flg2_exons | rna | gdna | 1 | 4.35% |
| hotspot_window | gdna | nrna | 1,737 | 37.15% |
| hotspot_window | gdna | mrna | 1,620 | 34.64% |
| hotspot_window | gdna | gdna | 1,307 | 27.95% |
| hotspot_window | gdna | unresolved | 12 | 0.26% |

## Approximate EM Log-Odds

For one unspliced fragment, the native VBEM E-step compares components roughly as:

```text
log posterior score_k = log_lik_k(fragment) + digamma(alpha_k) - log(effective_length_k) + constant
```

So for gDNA versus a transcript component:

```text
log odds(gDNA / RNA_k)
  ~= [log_lik_gDNA - log_lik_RNA_k]
   + [log(alpha_gDNA) - log(L_gDNA)]
   - [log(alpha_k) - log(L_k)]
```

Using the converged locus outputs as an approximation, the median fragment length in the window is 316 bp. The saved summary for this run does not include the RNA/gDNA FL histograms, so the table below excludes the FL likelihood ratio and shows the density plus strand terms. Negative total log-odds means the RNA component wins for a strand-compatible fragment before any additional FL-model effect.

| Component | Count | L_eff | Density log-odds | FL log-odds | Strand log-odds | Total log-odds |
| --- | --- | --- | --- | --- | --- | --- |
| mRNA FLG2 ENST00000388718.5 | 1,583 | 8,879 | -2.268 | n/a | -0.692 | -2.961 |
| nRNA RIGEL_NRNA_chr1_1_152335378_152364080 | 1,513 | 20,745 | -1.374 | n/a | -0.692 | -2.067 |
| nRNA RIGEL_NRNA_chr1_1_152168124_152332686 | 2,328 | 54,784 | -0.834 | n/a | -0.692 | -1.527 |

## Interpretation

The browser intuition is right: the local data look like gDNA, but the model can explain the same 50/50 strand mixture as the sum of two strand-specific RNA components. The key mathematical issue is that Rigel's gDNA competitor is one component for the entire MultiLocus, with `L_gDNA=92,407,095` bp in this run. FLG2 has an EM effective length on the order of kilobases, and the synthetic antisense nRNA spans are still far shorter than the MultiLocus gDNA opportunity.

That creates a strong per-base density penalty against gDNA. Even though the locus has many gDNA-assigned fragments overall, the gDNA abundance must be divided by the huge MultiLocus effective length before it competes with FLG2 or the antisense nRNA spans. For a strand-compatible FLG2 read, RNA also gets a strand-likelihood advantage over gDNA: the RNA component gets probability near 1 on its matching orientation, while gDNA pays a 0.5 strand-symmetry factor. The opposite read orientation is not evidence against RNA globally because a positive-strand nRNA/antisense component can absorb it.

This means the failure is both a likelihood-model problem and a prior/geometry problem:

1. Likelihood problem: there is no coupled local test that says `positive and negative unspliced reads at the same exonic interval should be explained by one symmetric gDNA source unless the positive-strand transcript has independent support elsewhere`.
2. Prior/effective-length problem: the gDNA component is normalized over the whole chromosome-scale MultiLocus, not over the local captured opportunity that generated these reads.
3. Model-structure problem: nRNA components are allowed to claim intronic overlap locally without a penalty for absent support in their other exons/introns.

## What We Need

The direct fix is not just another scalar kappa tweak. We need local evidence coupling in the model. The most actionable design is a local gDNA-versus-ambiguous-unspliced guard:

1. For an exon-contained unspliced pileup, aggregate both orientations over a small genomic segment before hard assignment or before EM prior construction.
2. Require independent RNA support for the antisense/nRNA explanation: splice junctions, transcript-unique exons, or reads in other introns/exons of the same nascent span.
3. Give gDNA a local captured-opportunity effective length for that segment, not the entire chromosome-scale MultiLocus denominator.
4. Penalize nRNA components whose only support is a local overlap with a high-depth opposite-strand mature exon and whose remaining span has near-zero support.

Artifacts are in `results/flg2_hotspot_diagnostics_2026-05-21`.
