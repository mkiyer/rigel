# VCaP Problem Regions For Hybrid-Capture Simulation - 2026-05-20

Inputs:

- Annotated BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix/annotated.bam`
- Locus output: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix/loci.feather`
- Region analysis: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/`
- Hotspot analysis: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/`

Coordinates below are 1-based inclusive for readability. The matching 0-based
half-open interval is `[start_1based - 1, end)`.

## Selected Targets

| Target | Region | 0-based half-open | Locus / MultiLocus | Error signal | Why this is useful |
| --- | --- | --- | --- | --- | --- |
| T1 | chr2:178,559,311-178,576,831 | chr2:178559310-178576831 | MultiLocus 3 | 4,607/6,838 true-gDNA fragments called RNA (67.37%); 4,266 nRNA, 341 mRNA | Strong mostly-gDNA-to-nRNA error in an ambiguous-strand EXON region with very little true RNA support. Good for isolating exon-union/nRNA siphon behavior. |
| T2 | chr11:65,497,640-65,508,073 | chr11:65497639-65508073 | MultiLocus 3 | 4,396/5,443 true-gDNA fragments called RNA (80.76%); 2,239 nRNA, 2,157 mRNA | Highest-rate large hotspot. Ambiguous-strand EXON, 70 exon-overlapping transcripts, 5 exon genes, and abundant true RNA nearby. Good stress test for multi-gene ambiguous exons. |
| T3 | chrX:67,544,021-67,546,762 | chrX:67544020-67546762 | Locus 21523 | 2,793/6,173 true-gDNA fragments called RNA (45.25%); 2,711 mRNA, 82 nRNA | Compact AR locus with strong gDNA-to-mRNA error. Useful non-mega comparison where the locus is small enough to simulate and inspect carefully. |

## Supporting Details

### T1: chr2 ambiguous EXON, MultiLocus 3

- Calibration-region row: `region_id=99068`, `region_type=EXON`, `region_strand=AMBIG`, length `17,521 bp`.
- Boundary flags: left and right are both eligible.
- Complexity: `16` exon-overlapping transcripts, `2` exon genes, `40` transcript spans, `2` span genes.
- Dominant-region evidence: `6,838` true-gDNA source fragments dominated by this region; `4,607` called RNA.
- Error split: `4,266` nRNA, `341` mRNA.
- Hotspot target: `nRNA chr2:178,525,988-178,804,642(-)`.
- Nearby 10 kb windows show the same behavior, e.g. `chr2:178,560,001-178,570,000` has `2,665/3,949` true-gDNA fragments called RNA.

### T2: chr11 ambiguous EXON, MultiLocus 3

- Calibration-region row: `region_id=374859`, `region_type=EXON`, `region_strand=AMBIG`, length `10,434 bp`.
- Boundary flags: neither side is boundary-flux eligible, which makes this a good test of how exon-contained expected gDNA is handled without local boundary evidence.
- Complexity: `70` exon-overlapping transcripts, `5` exon genes, `70` transcript spans, `5` span genes.
- Dominant-region evidence: `5,443` true-gDNA source fragments dominated by this region; `4,396` called RNA.
- Error split: `2,239` nRNA, `2,157` mRNA.
- Hotspot target: `nRNA ENST00000698129.1` plus mRNA calls in the same window.
- Nearby 10 kb window `chr11:65,500,001-65,510,000` has `3,698/4,206` true-gDNA fragments called RNA.

### T3: chrX AR-region compact locus, Locus 21523

- Calibration-region row: `region_id=640049`, `region_type=EXON`, `region_strand=POS`, length `2,742 bp`.
- Boundary flags: right side is boundary-flux eligible.
- Complexity: `8` exon-overlapping transcripts, `1` exon gene.
- Dominant-region evidence: `6,173` true-gDNA source fragments dominated by this region; `2,793` called RNA.
- Error split: `2,711` mRNA, `82` nRNA.
- Locus 21523 summary: `13` transcripts, `1` gene, `59,500` EM fragments, `gdna_rate=0.230607`, `gdna_prior_count_em=123.71`, `gdna_em_exposure_weight=0.230402`.
- Transcript-coordinate span for locus 21523 on chrX: `67,544,021-67,730,619`.
- Nearby 10 kb window `chrX:67,540,001-67,550,000` has `2,875/6,585` true-gDNA fragments called RNA.

## Notes For Simulation Design

- T1 and T2 both belong to MultiLocus 3, a very large cross-reference mega-component (`42,639` transcripts, `4,005` genes, `2,823,632` EM fragments). For simulation, use local windows around the listed coordinates rather than attempting to reproduce the full MultiLocus at first.
- T1 is the cleanest nRNA-siphon target: high gDNA-to-nRNA, low true RNA in the exact dominant region.
- T2 is the strongest ambiguous-exon stress target: very high false-RNA rate, multi-gene exon complexity, and substantial nearby true RNA.
- T3 is the best compact control target: strong error, a biologically recognizable locus, and a much smaller Locus/MultiLocus footprint.