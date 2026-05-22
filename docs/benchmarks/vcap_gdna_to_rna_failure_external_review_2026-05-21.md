# VCaP gDNA-to-RNA Leakage Failure Brief For External Review

Date: 2026-05-21

Primary benchmark: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m`

Purpose: this document is intended as a handoff to an outside reviewer. It describes
the current Rigel failure mode, the fixes and diagnostics we have already tried,
why the results are disappointing despite reasonable model ideas, and where we
think external advice would be most valuable.

## Executive Summary

Rigel is currently stuck on a recurring failure mode: true genomic DNA fragments
are being assigned to RNA, either as mature mRNA or as nascent RNA (nRNA). In the
VCaP real RNA/gDNA mixture, flowcell labels let us treat RNA and gDNA source as
known truth. The current best-performing recent run still calls roughly 13-16%
of true-gDNA fragments as RNA, depending on the calibration variant.

The central symptom is not a small numerical tuning problem. It is a model
identifiability problem exposed by high-depth, hybrid-capture-like, exon-rich
regions. Many true gDNA fragments are unspliced, genic, and compatible with both
gDNA and one or more RNA states. The EM can explain the two genomic-strand
orientations of a true gDNA pileup as a combination of strand-specific mRNA and
opposite-strand nRNA components. Once this happens, local evidence that looks
obviously like gDNA in a genome browser can be mathematically competitive with,
or even dominated by, RNA states.

We have tried several serious fixes:

1. **Regional/capture exposure weighting.** This helped after a kappa-units bug
   was fixed, but it did not solve the problem. It reduced one bad VCaP run from
   20.44% gDNA->RNA to 13.13%, but still left millions of false RNA calls.
2. **EXON-contained strand deconvolution and EXON-COMPOSITE density.** The
   theory was plausible: use stranded exon-contained counts to estimate exonic
   gDNA density instead of relying only on boundary flux. In practice, it made
   true-gDNA->RNA worse versus the immediate kappa-fix baseline: 2,393,558 to
   2,500,290 false RNA calls.
3. **Exact exposure-weighted gDNA denominator.** This was implemented, measured,
   and rolled back. It barely changed the FLG2 mega-locus denominator
   (`92,407,095` -> `92,459,136`, only 1.00056x), so it did not address the
   root failure.
4. **Changing the exposure density reference quantile from Q95 to Q99.9/Q99.95.**
   Offline diagnostics suggested this might shrink the gDNA denominator toward a
   useful local scale. Full VCaP quantification did the opposite of what we
   wanted: Q99.9 increased gDNA->RNA by 394,381 fragments, and Q99.95 increased
   it by 428,138 fragments versus Q95.

The ideas are not obviously foolish. They are attempts to put the same physical
quantity, local gDNA opportunity, into both the numerator and denominator of the
gDNA likelihood. The bad results suggest that the missing ingredient is not just
another scalar exposure normalizer. We likely need a different model structure:
local evidence coupling, stronger requirements for nRNA support, a real
hybrid-capture exposure model, or a local gDNA state/guard that prevents
unspliced exon-compatible gDNA pileups from being split into mRNA+nRNA by EM.

## What Rigel Is Trying To Do

Rigel jointly models three sources of RNA-seq-like fragments:

- mature RNA (`mRNA`) assigned to annotated transcript states;
- nascent RNA (`nRNA`) assigned to transcript-derived genomic/nascent spans;
- genomic DNA contamination (`gDNA`) assigned per locus.

The pipeline is roughly:

1. A native C++ scanner reads a name-sorted BAM, resolves fragments against a
   transcript/reference index, trains strand and fragment-length models, and
   buffers fragment-level candidate information.
2. Python calibration estimates global and regional gDNA density, fragment-length
   distributions, and locus priors.
3. A locus-level EM solver assigns ambiguous fragments across mRNA, nRNA, and a
   per-locus gDNA component.
4. Output BAM annotation writes tags that can be compared to external truth. For
   the VCaP benchmark, query-name flowcell IDs provide RNA/gDNA truth.

An important design principle in this repository is that genes are output
summaries only. Internal inference uses transcripts and transcript-derived
states. That matters here because the failure is not gene-level summarization;
it happens at the fragment-to-state level inside large transcript/nRNA/gDNA
loci.

## Benchmark And Truth Definition

The main real-mixture benchmark is:

```text
/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam
```

Truth source is derived from flowcell IDs in the query name:

| Flowcell | Truth source |
| --- | --- |
| `C6EL5ANXX` | RNA |
| `H7MFFDSXY` | gDNA |

The confusion analyses count one primary read1 record per fragment. Secondary,
supplementary, and read2 records are skipped. All recent confusion scans counted
`32,218,601` fragments:

| Truth source | Fragments |
| --- | ---: |
| RNA | 13,989,924 |
| gDNA | 18,228,677 |

The central metric is true-gDNA fragments predicted as RNA:

```text
gDNA -> RNA = gDNA -> mRNA + gDNA -> nRNA
```

This metric is the failure mode because it creates false transcript/nascent RNA
signal from contaminating DNA.

## Current Failure Mode

The current failure is broad, but it has a very recognizable anatomy.

### Symptom 1: Millions Of True gDNA Fragments Are Called RNA

The current EXON-COMPOSITE run (`exon_strand_deconv_v1`) has:

| Metric | Value |
| --- | ---: |
| true gDNA fragments | 18,228,677 |
| gDNA -> RNA | 2,500,290 |
| gDNA -> mRNA | 1,171,473 |
| gDNA -> nRNA | 1,328,817 |
| gDNA recall | 84.66% |
| RNA -> gDNA | 466,818 |

This is already better than the worst v4.3 regional-exposure run, but it is not
good enough for a decontamination/quantification tool. A false-RNA rate of
13.72% among true-gDNA fragments means that a large amount of DNA contamination
is being presented as biological RNA.

### Symptom 2: The Errors Are Mostly Unspliced, Genic, Ambiguous Fragments

In the current EXON-COMPOSITE run, 96.12% of gDNA->RNA false positives are
`ZS=unspliced`, and 96.97% are either `ambig_same_strand` or
`ambig_opp_strand`. The dominant categories are:

| Predicted pool | Category | Count | Fraction of false RNA |
| --- | --- | ---: | ---: |
| mRNA | `ambig_same_strand / unspliced` | 918,169 | 36.72% |
| nRNA | `ambig_same_strand / unspliced` | 854,615 | 34.18% |
| nRNA | `ambig_opp_strand / unspliced` | 361,080 | 14.44% |
| mRNA | `ambig_opp_strand / unspliced` | 217,716 | 8.71% |

This is not mostly a spliced-read artifact. It is not mostly an alignment
secondary/supplementary artifact. It is the core ambiguity of unspliced genic
sequence in an RNA/gDNA mixture.

### Symptom 3: The Problem Is Diffuse, With Strong Hotspots

The top 10 EM loci account for only 13.20% of all gDNA->RNA false positives, and
the top 10 genomic 10 kb windows account for only 1.13%. This means the problem
is recurrent across many genic ambiguous regions rather than one bad locus.

There are still striking hotspots. Examples from the current run:

| Window | gDNA false RNA | False rate among gDNA-source fragments | Dominant target |
| --- | ---: | ---: | --- |
| `chr11:65,500,001-65,510,000` | 3,809 | 90.56% | nRNA / mRNA ambiguity |
| `chr1:152,350,001-152,360,000` | 3,297 | 71.52% | FLG2 mRNA and antisense/nRNA |
| `chrX:67,540,001-67,550,000` | 2,890 | 43.89% | AR mRNA |
| `chr2:178,560,001-178,570,000` | 2,833 | 71.74% | nRNA span |

### Symptom 4: FLG2 Shows Why Browser Intuition And EM Diverge

The FLG2 hotspot at `chr1:152,350,001-152,360,000` is the clearest local example.
Within the 10 kb window:

| Local quantity | Value |
| --- | ---: |
| fragments overlapping window | 4,699 |
| true RNA-source fragments | 23 |
| true gDNA-source fragments | 4,676 |
| true gDNA assigned RNA | 3,357 |
| local true-gDNA false-RNA rate | 71.79% |

The local data are overwhelmingly gDNA by flowcell truth. However, Rigel can
explain the two genomic orientations as two RNA stories:

- negative-strand FLG2 mature mRNA explains one orientation;
- positive-strand antisense/nRNA explains the other orientation.

The local false-RNA targets are concentrated:

| Predicted target | Count | False-RNA fraction |
| --- | ---: | ---: |
| mRNA FLG2 `ENST00000388718.5` | 1,581 | 47.10% |
| nRNA `chr1:152,335,378-152,364,080(+)` | 1,504 | 44.80% |
| nRNA `chr1:152,348,734-152,360,006(-)` | 212 | 6.32% |

This is the model failure in miniature: a true symmetric gDNA pileup can be
split into strand-specific RNA states if both strand-compatible RNA candidates
exist locally.

## Why The EM Chooses RNA In These Cases

For an unspliced fragment, the native VBEM E-step approximately compares:

```text
log posterior score_k
  = log_lik_k(fragment)
  + digamma(alpha_k)
  - log(effective_length_k)
  + constant
```

For gDNA versus an RNA component:

```text
log odds(gDNA / RNA_k)
  ~= [log_lik_gDNA - log_lik_RNA_k]
   + [log(alpha_gDNA) - log(L_gDNA)]
   - [log(alpha_RNA_k) - log(L_RNA_k)]
```

In FLG2/locus 3, the gDNA competitor is one component for a huge MultiLocus.
The run has:

| Locus 3 quantity | Value in Q95 / EXON-COMPOSITE run |
| --- | ---: |
| merged footprint | 353,317,609 bp |
| unweighted FL-marginal gDNA denominator | 353,762,179 bp |
| exposure-weighted gDNA denominator passed to EM | 92,407,095 bp |
| exposure weight | 0.261212 |
| EM fragments in locus | 2,823,632 |
| transcript components | 42,639 |
| nRNA entities | 22,701 |

The denominator is already exposure-weighted, but it is still tens of millions
of bases. Local RNA components have effective lengths on the order of kilobases
to tens of kilobases. Even with many gDNA-assigned fragments in the locus, the
gDNA abundance is spread over a huge opportunity. For strand-compatible RNA, the
RNA component also gets a strand-likelihood advantage while gDNA pays a
strand-symmetry factor.

The result is that local evidence which biologically looks like DNA can lose to
RNA in the model because the model lacks a coupled local test that says:

```text
opposite orientations at the same exon-contained interval are more parsimonious
as one symmetric gDNA source unless the opposite-strand RNA/nRNA explanation has
independent support elsewhere.
```

## Chronology Of What We Tried

### 1. gDNA Effective-Length Contract

We changed the gDNA scoring/EM contract so that gDNA uses a component-level
effective length analogous to RNA. The native scorer emits gDNA fragment-length
terms and non-length terms, while the EM receives a per-MultiLocus gDNA
effective length. This was necessary to make gDNA/RNA comparisons dimensionally
consistent.

This was a foundational correction, not by itself a solution to gDNA->RNA
leakage. It made later regional exposure work possible and auditable.

### 2. Regional Exposure Weighting

The regional exposure model introduces a spatial opportunity field:

```text
A_r in (0, 1]
```

The intended physical model is:

```text
local gDNA density ~= global gDNA density * local opportunity
```

The same exposure field is supposed to enter:

1. the gDNA numerator, by adding `log A(midpoint)` to a gDNA candidate fragment;
2. the gDNA denominator, by weighting the effective length/opportunity.

This is conceptually sound. It is exactly the kind of correction needed for
hybrid-capture data, where the on-target gDNA opportunity can be much smaller
than a whole gene or mega-locus footprint.

However, the first v4.3-style regional exposure run regressed badly:

| Run | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA |
| --- | ---: | ---: | ---: | ---: | ---: |
| v4.3 with multimappers | 77.94% | 3,725,095 (20.44%) | 1,346,710 | 2,378,385 | 383,809 |

The root cause was a kappa-units error. The empirical-Bayes kappa returned by
`estimate_kappa` is a negative-binomial/Gamma-Poisson shape parameter, not a bp
opportunity length. Using it directly made empty or low-evidence regions shrink
toward nearly zero exposure instead of toward the global mean. The correction is
to convert alpha to opportunity units:

```text
beta = alpha / rho_global
```

After fixing this kappa-units issue, the VCaP confusion improved substantially:

| Run | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA |
| --- | ---: | ---: | ---: | ---: | ---: |
| v4.3 with multimappers | 77.94% | 3,725,095 (20.44%) | 1,346,710 | 2,378,385 | 383,809 |
| kappa-units fix | 85.24% | 2,393,558 (13.13%) | 1,062,195 | 1,331,363 | 508,599 |

This is the one effort that clearly helped gDNA->RNA leakage. But it introduced
or exposed a tradeoff: RNA->gDNA worsened from 383,809 to 508,599. More
importantly, it still left 2.39 million true-gDNA fragments called RNA.

### 3. EXON-Contained Strand Deconvolution And EXON-COMPOSITE Density

Boundary flux alone is a weak proxy for exonic gDNA density. The EXON-contained
strand-deconvolution feature tried to estimate gDNA density inside mature exon
regions using strand information. In a stranded library, mature RNA is strongly
oriented, while gDNA should be symmetric. The idea was to subtract the
strand-asymmetric RNA component and obtain a cleaner exonic gDNA density.

The implementation added:

- native payload counts for exon-contained fragments by orientation;
- global `EXON-CONTAINED` density;
- precision-weighted fusion of boundary and contained evidence into an
  `EXON-COMPOSITE` density;
- regional exposure rows using separate boundary and exon-contained evidence;
- diagnostics for boundary-only, contained-only, and composite channels.

The feature did what it was mathematically designed to do:

| Density channel | Value |
| --- | ---: |
| EXON-INTRON boundary rho | 0.02990 |
| EXON-CONTAINED rho | 0.01977 |
| EXON-COMPOSITE rho | 0.02381 |

The contained channel estimated lower exonic gDNA density than boundary flux,
and the composite moved downward. But the benchmark result got worse relative to
the kappa-units baseline:

| Run | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA |
| --- | ---: | ---: | ---: | ---: | ---: |
| kappa-units fix | 85.24% | 2,393,558 (13.13%) | 1,062,195 | 1,331,363 | 508,599 |
| EXON-COMPOSITE v1 | 84.66% | 2,500,290 (13.72%) | 1,171,473 | 1,328,817 | 466,818 |

Interpretation: EXON-contained deconvolution helped one side of the tradeoff
(RNA->gDNA improved by 41,781), but it made gDNA->RNA worse by 106,732. The
regression was almost entirely gDNA->mRNA. The likely mechanism is that the
contained exon channel correctly sees strong strand-asymmetric RNA signal and
therefore lowers exonic gDNA density, but that also weakens gDNA protection for
true captured gDNA fragments lying inside mature exons.

This is a key lesson: a correction that is locally reasonable for removing RNA
from the gDNA density estimate can still be harmful if ambiguous gDNA fragments
need that exonic prior to avoid being siphoned into RNA states.

### 4. Exact Exposure-Weighted gDNA Effective-Length Trial

Production was using a scalar approximation:

```text
L_gDNA_EM = gdna_eff_len_for_loci(locus) * footprint_exposure_weight(locus)
```

The repository had previously contained an exact FL-aware weighted denominator
helper that integrates the exposure field over fragment-length-expanded midpoint
windows. We suspected that the scalar approximation might be why the FLG2
mega-locus denominator remained enormous.

The exact denominator trial was implemented, measured on VCaP, then removed
after it did not materially change the problematic denominator:

| Locus 3 metric | Scalar run | Exact trial | Interpretation |
| --- | ---: | ---: | --- |
| `gdna_eff_len_unweighted` | 353,762,179 | 353,762,179 | raw opportunity unchanged |
| `gdna_eff_len` | 92,407,095 | 92,459,136 | exact was only 1.000563x scalar |
| `gdna_em_exposure_weight` | 0.261212 | 0.261360 | essentially unchanged |
| `gdna_prior_count_em` | 259,499 | 259,645 | prior scaled consistently |

This ruled out exact projection geometry as the main cause of the 92M FLG2
denominator. The problem is the exposure field and model structure, not the
small fragment-length shoulder difference between scalar and exact integration.

### 5. Cap-At-1 And Density Reference Quantile Diagnostics

The exposure field maps regional density estimates to exposure weights roughly as:

```text
A_r = clip(rho_hat_r / rho_ref, floor=1e-4, upper=1)
```

The default reference was global Q95. In FLG2/locus 3, this produced a mean
exposure of `0.2612`, so the gDNA denominator remained `92.4M`. Diagnostics
showed that raising the reference scale could shrink the denominator:

| Normalizer | Mean exposure | L_gDNA |
| --- | ---: | ---: |
| current global Q95 | 0.261213 | 92,407,283 |
| global Q99 | 0.101606 | 35,944,243 |
| global Q99.9 | 0.026853 | 9,499,710 |
| class Q99.9 | 0.015386 | 5,442,962 |
| locus Q99.9 | 0.014838 | 5,249,150 |
| global max | 0.000118 | 41,737 |

This looked promising because hotspot math suggested that FLG2 gDNA needed a
denominator closer to a few million or lower to become competitive with local
RNA components. We therefore made the regional exposure reference quantile
configurable and ran full VCaP quantification at Q99.9 and Q99.95.

Full validation was bad:

| Run | rho_ref | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Q95 baseline | 0.0003609 | 84.66% | 2,500,290 (13.72%) | 1,171,473 | 1,328,817 | 466,818 |
| Q99.9 | 0.006698 | 82.49% | 2,894,671 (15.88%) | 1,085,039 | 1,809,632 | 459,112 |
| Q99.95 | 0.01010 | 82.31% | 2,928,428 (16.06%) | 1,087,668 | 1,840,760 | 455,500 |

The higher reference quantiles reduced gDNA->mRNA slightly, but they greatly
increased gDNA->nRNA:

| Comparison | Delta gDNA -> RNA | Delta gDNA -> mRNA | Delta gDNA -> nRNA |
| --- | ---: | ---: | ---: |
| Q99.9 - Q95 | +394,381 | -86,434 | +480,815 |
| Q99.95 - Q95 | +428,138 | -83,805 | +511,943 |

At FLG2/locus 3, the denominator and prior collapsed:

| Run | Exposure weight | L_gDNA | Prior EM | EM gDNA | EM mRNA | EM nRNA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Q95 | 0.261212 | 92,407,095 | 259,499 | 1,446,008 | 1,162,014 | 288,677 |
| Q99.9 | 0.026882 | 9,510,006 | 26,706 | 1,378,282 | 1,140,289 | 378,128 |
| Q99.95 | 0.018469 | 6,533,723 | 18,348 | 1,372,478 | 1,139,974 | 384,247 |

The high-tail quantiles moved fragments from gDNA into nRNA rather than fixing
false-RNA assignment. This is a particularly important negative result because
it shows that simply making the gDNA denominator smaller is not sufficient under
the current prior/EM coupling and state structure.

### 6. Other Diagnostic Efforts

We also performed several targeted diagnostic passes:

- **Hotspot analysis.** Identified that false RNA calls are mostly unspliced and
  genic, with FLG2 and other hotspots showing local gDNA truth but RNA EM wins.
- **Local RNA evidence stratification.** False RNA calls occur both in windows
  with high RNA-source coverage and in windows with little or no RNA-source
  evidence. This separates true expression pull from self-seeding ambiguous RNA
  components.
- **Cap-at-1 decomposition.** Showed that Q95 saturation is not the only issue:
  saturated regions contribute heavily, but broad moderately high exposure in
  the `[0.1,1)` bin also keeps the denominator large. A soft cap at Q95 still
  gives a 65.1M FLG2 denominator.
- **Class-aware normalization theory.** Considered per-class exposure references
  because region classes have different physical density scales. We became
  concerned that class-aware normalization has a weak identifying basis and may
  hide model errors rather than solve them. It remains a diagnostic idea, not a
  validated production fix.

## What We Think Is Fundamentally Wrong

We do not think this is now primarily a coding bug. The code paths are tested,
the changes are doing what they were mathematically designed to do, and the
benchmarks are consistently telling us that the current model is missing
structure.

### 1. Local gDNA Is Not A Local State

The gDNA competitor is one component per MultiLocus. In mega-loci, especially
those created by multimapper connectivity, this means gDNA is normalized over a
huge footprint. Local RNA states can be much shorter and therefore much more
competitive per effective base. We tried exposure weighting to shrink the gDNA
opportunity, but the state is still a broad locus-level object rather than a
local source.

### 2. nRNA States Are Too Easy To Invoke Locally

nRNA components can absorb unspliced signal in a local overlap even if the rest
of the nRNA span has little independent support. In FLG2, a positive-strand nRNA
span helps explain the opposite orientation of a true gDNA pileup. The model has
no strong penalty for this local-only explanation.

### 3. Opposite-Strand Evidence Is Not Coupled

gDNA is symmetric; mature/nascent RNA is strand-specific. But the model compares
individual fragments to individual states. It does not aggregate a local pileup
and ask whether the two orientations are jointly better explained by one
symmetric gDNA source. As a result, a 50/50 gDNA pileup can be split into two RNA
components.

### 4. Exposure Is Not Identified From Ambiguous Exonic Counts

Hybrid capture makes local gDNA opportunity highly nonuniform. But exonic counts
are confounded by mature RNA expression and nRNA. If ambiguous exon-contained
counts influence the exposure field too freely, we can either under-protect true
gDNA in exons or create circular self-rescue. We need a stronger identifying
assumption or external information, such as a capture target BED or control data.

### 5. Scalar Normalizers Change Tradeoffs, Not The Root Cause

Kappa, Q95/Q99.9, class-relative references, and caps all move mass between
gDNA, mRNA, and nRNA. They have not removed the false-RNA pressure. Recent
changes often improve one side of the confusion matrix while worsening another:

- kappa-units fix improved gDNA recall but worsened RNA->gDNA;
- EXON-COMPOSITE improved RNA->gDNA but worsened gDNA->RNA;
- high-tail quantiles reduced gDNA->mRNA but greatly worsened gDNA->nRNA.

This pattern is why we feel stuck. The model has enough flexibility to move the
error around, but not enough structure to identify the true source reliably.

## Current State Of The Tool

Rigel is in a technically strong but scientifically unresolved state.

What is working:

- The scanner, calibration payloads, EM, and output annotation run at full VCaP
  scale.
- The codebase has focused tests for the implemented calibration and exposure
  components.
- The output diagnostics are now rich enough to localize failures by truth
  source, predicted pool, locus, genomic window, fragment category, and target.
- Several previous bugs or weak assumptions have been found and corrected,
  especially the gDNA effective-length contract and kappa-units issue.

What is not working:

- The tool is not yet reliable for separating gDNA from RNA in high-depth,
  exon-rich, hybrid-capture-like mixtures.
- The remaining gDNA->RNA leakage is too large to treat as a tolerable edge
  case.
- The current model can confidently assign true gDNA fragments to RNA with high
  winner weights in biologically implausible situations.
- The fixes we have tried are increasingly moving the error between mRNA and
  nRNA rather than reducing total false RNA.

For users, the practical warning is: Rigel's RNA abundance estimates should be
treated cautiously in unspliced ambiguous exonic regions, especially in large
mega-loci, high-expression/capture hotspots, or loci with antisense/nRNA
candidate states. The tool can still be useful, but the current output should
not be considered a solved decontamination result for VCaP-like data.

## Recommendations For What To Do Next

### Recommendation 1: Add A Local gDNA-Versus-RNA Guard

The most direct missing structure is local coupling. Before allowing a set of
unspliced exon-contained fragments to be assigned to mRNA/nRNA, aggregate a small
genomic segment and evaluate whether the local orientation balance and support
pattern are more consistent with gDNA.

A candidate rule/model:

1. Segment ambiguous exon-contained unspliced evidence into local windows or
   region-partition intervals.
2. Count both orientations together.
3. If both orientations are present and there is weak independent RNA evidence,
   boost or floor the local gDNA posterior/prior.
4. Require independent RNA support before the opposite-strand/nRNA explanation
   can win: splice junctions, transcript-unique exons, coverage outside the
   local overlap, or coherent nascent-span support.

This could be implemented first as a conservative guard or diagnostic, then
folded into the generative model later. It directly targets the FLG2 anatomy.

### Recommendation 2: Make nRNA Support Nonlocal

nRNA should not be able to win solely because it overlaps a local exonic pileup.
For each nRNA candidate, introduce a support score that looks outside the local
mature-exon overlap:

- intronic or pre-mRNA body support away from mature exon-only sequence;
- continuity across the nascent span;
- strand-consistent evidence in other parts of the span;
- absence of support in unrelated span regions as a penalty.

This is not just a heuristic. It encodes a biological distinction: nascent RNA is
a molecule with genomic extent, not an arbitrary local absorber for whichever
orientation gDNA needs.

### Recommendation 3: Build A Real Hybrid-Capture Exposure Model

The exposure model should explicitly represent assay target/off-target states.
If a capture BED or bait design is available, use it. If not, infer a low-
complexity exposure field from conservative gDNA-proxy evidence.

Important constraints:

- Do not learn the target map from the same ambiguous exon-contained fragments
  that need to be classified.
- Treat mature RNA coverage as a confounder, not target evidence.
- Use intergenic, intronic, boundary-crossing, and high-confidence gDNA posterior
  evidence as anchors.
- Keep the exposure field low-complexity: a small number of shared exposure
  states, spatial smoothing, and shrinkage to uniform when evidence is weak.
- Use the same `A(x)` in the numerator and denominator.

The current Q95/Q99.9 scalar reference experiments are not a substitute for this
model. They are only global normalization conventions.

### Recommendation 4: Consider Local gDNA States Or Locus Splitting

One gDNA component per huge MultiLocus may be too coarse. Options:

- split mega-loci into local gDNA subcomponents even when RNA transcript states
  remain connected by multimapping;
- introduce tiled/local gDNA states with hierarchical shrinkage to a locus-level
  gDNA abundance;
- perform EM in a way that lets gDNA compete locally while keeping global
  contamination calibrated.

This is a bigger architectural change, but it may be necessary. The FLG2 failure
is partly caused by forcing a local captured gDNA pileup to compete through a
chromosome-scale gDNA denominator.

### Recommendation 5: Add An Ambiguous/Abstain Output State

When gDNA and RNA are not identifiable, forcing a hard mRNA/nRNA/gDNA assignment
can create false biological signal. A practical improvement may be to emit an
`ambiguous_gdna_rna` or low-confidence category for fragments/loci where:

- the winning RNA state is unspliced and exon-contained;
- local opposite-strand balance is compatible with gDNA;
- independent RNA evidence is weak;
- posterior is sensitive to exposure normalizer or prior settings.

This would not solve quantification alone, but it would prevent overconfident
false RNA output and provide a safer downstream interpretation.

### Recommendation 6: Build A Focused Synthetic Benchmark For This Failure

The VCaP benchmark is valuable but complex. We need a synthetic or semi-synthetic
benchmark that isolates the current failure:

- one mature exon with true RNA on one strand;
- true gDNA sampled symmetrically across the exon;
- an antisense or nRNA candidate overlapping the same interval;
- configurable capture on/off exposure;
- known ground truth for every fragment.

This benchmark should be part of the acceptance suite for any future fix. The
acceptance metric should be total gDNA->RNA leakage, split by mRNA and nRNA, not
only aggregate EM likelihood or total gDNA mass.

## Directions We Should Avoid As Primary Fixes

The recent evidence argues against spending more time on these as standalone
solutions:

1. **More scalar quantile tuning.** Q99.9 and Q99.95 made leakage worse.
2. **Unconstrained class-aware normalization.** It may produce attractive local
   denominators but lacks a clean identifying assumption unless tied to a real
   assay exposure model.
3. **More kappa-only tuning.** The kappa-units bug was real and fixing it helped,
   but the remaining leakage is not just shrinkage strength.
4. **Exact denominator geometry alone.** We measured it; it does not move FLG2.
5. **EXON-contained deconvolution as an unbounded density-lowering correction.**
   It can protect RNA from gDNA overcalling but under-protects true exonic gDNA.

## Specific Questions For An External Reviewer

We would benefit most from advice on these points:

1. Is the local coupling diagnosis correct? In particular, should opposite-strand
   unspliced reads in the same exon-contained segment be modeled jointly before
   assigning to separate RNA states?
2. What is the most principled way to let gDNA compete locally without exploding
   the EM state space in mega-loci?
3. How should nRNA candidate states be penalized or validated when their support
   is local and overlaps a highly expressed mature exon?
4. Can a no-BED hybrid-capture exposure model be identified from conservative
   gDNA-proxy evidence, or do we need target BED/control data for this problem?
5. Would a hierarchical model with local gDNA tiles and locus-level/global
   shrinkage be preferable to the current one-gDNA-state-per-MultiLocus design?
6. Should ambiguous gDNA/RNA fragments be allowed to abstain from transcript
   quantification, and if so, how should that uncertainty propagate to output?
7. Are there known RNA-seq/gDNA contamination methods that handle this exact
   ambiguity better, especially in capture or high-gDNA mixtures?

## Key Artifacts For Review

Primary validation reports:

- `docs/benchmarks/vcap_exon_strand_deconv_validation_2026-05-20.md`
- `docs/benchmarks/vcap_global_high_tail_exposure_validation_2026-05-21.md`
- `docs/benchmarks/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20.md`
- `docs/benchmarks/flg2_hotspot_diagnostics_2026-05-21.md`
- `docs/benchmarks/locus3_exposure_audit_2026-05-21.md`
- `docs/benchmarks/vcap_locus3_cap_at_1_diagnostics_2026-05-21.md`

Supporting design/implementation notes:

- `docs/calibration/hybrid_capture_exposure_model_2026-05-18.md`
- `docs/calibration/exact_gdna_eff_len_exposure_plan_2026-05-21.md`
- `docs/calibration/exon_strand_deconv_plan_v1.md`
- `docs/calibration/exon_strand_deconv_implementation_log_v1.html`

Key result tables:

- `results/vcap_exon_strand_deconv_validation_2026-05-20/run_metrics.tsv`
- `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion_detail_compare.tsv`
- `results/vcap_global_high_tail_exposure_2026-05-21/run_metrics.tsv`
- `results/vcap_global_high_tail_exposure_2026-05-21/confusion_detail_compare.tsv`
- `results/flg2_hotspot_diagnostics_2026-05-21/`
- `results/vcap_cap_at_1_diagnostics_2026-05-21/`

Recent run directories:

- Q95 / EXON-COMPOSITE: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`
- global Q99.9: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/global_q999_v1`
- global Q99.95: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/global_q9995_v1`

## Bottom Line

We are no longer short on plausible scalar corrections. We have tried several,
and the best recent ideas either helped only modestly or made the primary
gDNA->RNA leakage worse. The current evidence points to a deeper identifiability
and model-structure problem: local, symmetric gDNA evidence is being compared
against flexible strand-specific RNA/nRNA states without enough local coupling or
independent support requirements.

The next useful step is not another Q value or kappa tweak. It is external review
of the generative model and a redesign around local gDNA competition,
nonlocal nRNA support, and capture-aware exposure identification.