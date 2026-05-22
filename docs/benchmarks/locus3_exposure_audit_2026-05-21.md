# Locus 3 Exposure-Weighted Opportunity Audit - 2026-05-21

Input run: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`

## Saved Denominator

Exposure weighting is active, but it is weak for this MultiLocus.

| Metric | Value |
| --- | --- |
| locus_id | 3 |
| locus_span_bp | 353,317,609 |
| gdna_eff_len_unweighted | 353,762,179 |
| gdna_em_exposure_weight | 0.261212 |
| gdna_eff_len passed to EM | 92,407,095 |
| reduction factor | 3.83x |
| weight needed for 10 kb opportunity | 2.827e-05 |
| current / 10 kb opportunity | 9240.7x |

Important correction: `92.4M` is not the raw locus size. It is already the exposure-weighted denominator. The unweighted FL-marginal gDNA denominator is `353.8M`, close to the saved `locus_span_bp` plus the FL expansion term.

## EM-Relevant Thresholds

These thresholds ask: how small would `L_gDNA` need to be for a strand-compatible gDNA fragment to break even against this RNA component, assuming the current alpha values and the gDNA `0.5` strand factor?

| Component | Current L_gDNA | Break-even L_gDNA | Current / break-even | Break-even weight |
| --- | --- | --- | --- | --- |
| mRNA FLG2 ENST00000388718.5 | 92,407,095 | 4,781,453 | 19.33 | 1.352e-02 |
| nRNA chr1_1_152335378_152364080 | 92,407,095 | 11,688,385 | 7.91 | 3.304e-02 |
| nRNA chr1_1_152168124_152332686 | 92,407,095 | 20,063,223 | 4.61 | 5.671e-02 |

The current denominator is about 19.3x too large for gDNA to beat FLG2 mRNA on the strand-compatible half of the pileup.

## Regional Exposure Summary

The exposure model itself is in regional mode. However, the reference density is the global 95th percentile (`rho_ref=3.609e-04`), and weights are capped at 1. Regions above that reference all saturate. In this run, ordinary exon composite density is already above the reference median, so many exon/capture-like regions are not distinguished from the FLG2 hotspot.

| Class | rho q50 | rho q95 | Opportunity | Evidence |
| --- | --- | --- | --- | --- |
| INTERGENIC | 7.571e-06 | 8.506e-05 | 1,308,005,751 | 10,370 |
| INTRON | 5.238e-05 | 4.863e-04 | 1,541,887,073 | 301,316 |
| EXON-COMPOSITE | 6.098e-04 | 6.281e-03 | 14,175,229 | 1,469,756 |
| EXON-INTRON | 0.000e+00 | 1.931e-01 | 180,896,006 | 5,407,683 |
| EXON-CONTAINED | 0.000e+00 | 1.857e-03 | 96,253,082 | 1,578,425 |

## RNA Component Exposure Weights

| Set | N | Min | P05 | P50 | P95 | Max | Count-weighted mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| annotated RNA components in locus 3 | 21,233 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| nRNA entities in locus 3 | 21,440 | 0.0188 | 0.1300 | 0.5841 | 1.0000 | 1.0000 | 0.7040 |

Top locus-3 RNA/nRNA components:

| Component | Kind | Count | EM effective length | EM weight |
| --- | --- | --- | --- | --- |
| ENST00000618132.2 | annotated | 25,056 | 4,738 | 1.0000 |
| ENST00000512805.6 | annotated | 19,963 | 895 | 1.0000 |
| ENST00000710862.1 | annotated | 16,710 | 6,450 | 1.0000 |
| ENST00000710935.1 | annotated | 11,160 | 1,609 | 1.0000 |
| ENST00000426088.5 | annotated | 11,134 | 7,189 | 1.0000 |
| ENST00000367545.8 | annotated | 11,070 | 12,673 | 1.0000 |
| chr2_2_178525988_178804642 | nRNA | 15,151 | 272,987 | 0.9805 |
| chr2_1_178523417_178776534 | nRNA | 7,746 | 248,693 | 0.9835 |
| chr1_1_203795653_203802173 | nRNA | 2,774 | 6,275 | 1.0000 |
| chr1_1_152168124_152332686 | nRNA | 2,328 | 54,784 | 0.3334 |
| ENST00000534336.3 | nRNA | 1,931 | 10,189 | 1.0000 |
| chr5_2_141136682_141241744 | nRNA | 1,869 | 75,103 | 0.7165 |

## Footprint Reconstruction Check

Using saved annotated transcript rows plus parseable nRNA span ids, the reconstructed merged footprint is 353,318,962 bp across 1,497 merged intervals. This is close enough to the saved 353,317,609 bp to confirm that locus 3 is a large connected component, not a 10 kb local problem.

## Code-Path Finding

The production prior path computes:

```text
L_gDNA_EM = gdna_eff_len_for_loci(ml.loci) * footprint_exposure_weight(ml.loci)
```

That is an unweighted FL-expanded denominator multiplied by one scalar bp-weighted average over the unexpanded MultiLocus footprint. The repository also contains `weighted_gdna_eff_len_for_loci(...)`, which directly integrates the exposure field over the FL-expanded midpoint windows, but `assemble_priors(...)` does not use it.

This explains the saved number mechanically: exposure weighting is enabled, but for locus 3 it collapses to a broad scalar `0.261`, not a local/capture opportunity near kb scale.
