# gDNA Regional Exposure Revised Plan v1

**Status**: design plan, not implemented  
**Date**: 2026-05-18  
**Supersedes directionally**: `gdna_exposure_plan_v3.md` production strategy  
**Related diagnosis**: `gdna_exposure_class_normalization_fix_plan_2026-05-18.md`

## 1. Executive Summary

The implemented v3 regional exposure model has the right instinct but the wrong
competitive structure. Hybrid capture changes the probability that a fragment
from a genomic location is observed. That is an assay observation-process term,
not a gDNA-specific molecule-identity term.

The v3 implementation applies regional exposure only to gDNA:

```text
gDNA numerator     += log A(x_u)
gDNA denominator    = weighted L_gDNA
nRNA numerator      = unchanged
nRNA denominator    = unchanged
```

This makes regional exposure a direct likelihood-ratio penalty against gDNA and
in favor of nRNA for low-exposure intronic fragments. That is unfair in a
hybrid-capture RNA-seq assay because nascent RNA molecules from the same low
exposure intron are also hard to observe.

The revised principle is:

```text
Regional exposure must be applied to every candidate component whose molecule
would occupy the exposed genomic location.
```

For the immediate gDNA-to-nRNA regression, that means symmetric regional
exposure for unspliced genomic hypotheses:

```text
gDNA: weighted numerator and weighted effective length
nRNA: weighted numerator and weighted effective length
mRNA: unchanged in v1 unless a later validation shows it is required
```

This plan starts from a clean slate. The class-normalization bug remains real,
but it is no longer the conceptual center of the fix. The center is fairness:
the capture exposure field should affect gDNA and nRNA symmetrically whenever
they explain the same unspliced genomic fragment.

## 2. What Went Wrong In v3

### 2.1 Cross-Class Normalization Bug

v3 computes per-class exposure references and then collapses them into one
global maximum:

```python
rho_ref = max(rho_ref_per_class.values())
```

In the VCaP hybrid-capture mix, the EXON-INTRON boundary channel dominates that
global reference:

```text
EXON-INTRON rho_ref_class  = 1.928126e-01
INTRON rho_ref_class       = 3.792302e-04
INTERGENIC rho_ref_class   = 1.065998e-09
```

This incorrectly makes ordinary intronic regions look extremely depleted. A
class-specific reference is necessary.

### 2.2 Incompatible Evidence Units

Even class-specific normalization does not fully solve the problem because the
three current evidence channels are not physically exchangeable.

| Class | Exposure denominator | Conservative count | Physical meaning |
|---|---:|---:|---|
| `INTERGENIC` | contained-fragment opportunity | contained intergenic fragments | broad off-target genomic sampling |
| `INTRON` | contained-fragment opportunity | contained intronic fragments | unspliced intron-contained sampling |
| `EXON-INTRON` | boundary-crossing side opportunity | exon-intron boundary fragments | capture shoulder / exon-edge flux |

Boundary-crossing density near an exon is not the same quantity as contained
intronic density. It can inform local capture exposure, but it should not be
treated as a directly comparable per-base density without a model that maps it
onto a common assay-accessibility field.

### 2.3 Structural Asymmetry Against gDNA

The deeper failure is the asymmetric use of `A(x)`. The EM score is:

```text
score_k(u) = log_lik_k(u) + log theta_k - log L_k
```

v3 changes only the gDNA side:

```text
score_gDNA(u) = base_gDNA(u) + log A(x_u) + log theta_g - log L_g^A
score_nRNA(u) = base_nRNA(u)              + log theta_n - log L_n
```

Relative to nRNA, the gDNA exposure contribution is:

```text
Delta_u = log A(x_u) - log(L_g^A / L_g)
        = log A(x_u) - log Abar_locus
        = log(A(x_u) / Abar_locus)
```

Where `Abar_locus` is the FL-weighted average gDNA exposure over the locus.

In a hybrid-capture locus with high-exposure exons and low-exposure introns,
`Abar_locus` tracks the exon-enriched opportunity. A true gDNA fragment in a
low-exposure intron then receives a large negative shift relative to nRNA. The
observed regression is exactly this failure: true unspliced gDNA fragments move
to nRNA.

This failure persists even if class-specific normalization is fixed. The
penalty may shrink, but any high-variance within-locus exposure field can still
penalize gDNA against unweighted nRNA.

## 3. Revised Design Principle

Regional exposure is an assay accessibility term. It should answer:

```text
Given that a molecule generated sequence from genomic position x, how likely is
this assay to recover that molecule?
```

It should not answer:

```text
Is this molecule gDNA rather than nRNA?
```

Therefore, when two hypotheses explain the same observed unspliced genomic
fragment at the same location, they should receive the same local exposure
factor.

The corrected component score is:

```text
score_k(u) = base_log_lik_k(u)
           + log A_k(u)
           + log theta_k
           - log L_k^A
```

Where `A_k(u)` is defined by the genomic footprint of the molecule under
component `k`, and `L_k^A` is the corresponding exposed effective length.

## 4. Component-Specific Exposure Rules

### 4.1 gDNA

gDNA is an unspliced genomic molecule. Its exposure is defined on the genomic
footprint of the fragment.

Numerator:

```text
base_gDNA_log_lik[u] += log A(midpoint_or_footprint_u)
```

Denominator:

```text
L_gDNA^A(M) = sum_ell h_G(ell) * integral over valid gDNA midpoint windows A(x) dx
```

This is the same denominator family implemented by
`weighted_gdna_eff_len_for_loci()`.

### 4.2 nRNA

nRNA is also an unspliced genomic molecule. It should receive the same assay
exposure treatment as gDNA for the same genomic location.

Numerator:

```text
base_nRNA_log_lik[u, t] += log A(midpoint_or_footprint_u)
```

for synthetic nRNA candidate transcript `t`, provided the candidate explains
the fragment as unspliced genomic sequence on the same reference.

Denominator:

```text
L_nRNA^A(t) = sum_ell h_RNA(ell) * integral over valid nRNA midpoint windows A(x) dx
```

For a synthetic nRNA span `[start_t, end_t)`, the valid windows are the same
fragment-length-valid windows used for an unspliced genomic transcript span.
This can be computed from `index.t_df.ref`, `start`, `end`, `length`, and the
synthetic/nRNA flags.

Important: nRNA is unspliced by definition. Spliced alignments should not be
routed to nRNA and should not use nRNA exposure logic.

### 4.3 mRNA

mRNA is a mature spliced molecule. It does not physically contain introns. In
v1, leave mRNA on the existing transcript effective-length model.

Rationale:

1. The catastrophic regression is gDNA to nRNA on unspliced fragments.
2. Mature mRNA exposure requires projection through exon blocks and splice
   structure, which is a larger modeling problem.
3. Applying symmetric exposure to gDNA and nRNA should remove the unfair
   low-intron penalty without destabilizing mature transcript quantification.

Future work may add mRNA exposure over exonic transcript footprints if VCaP or
other hybrid-capture truth data show a remaining mature-RNA bias.

## 5. Revised Exposure Field `A(x)`

The exposure field should be a single assay-accessibility field, not three
unrelated class-specific density scales stitched together.

### 5.1 Short-Term Field For v1

For a first implementation, keep the existing `RegionalGdnaExposure` table shape
but repair the most dangerous construction choices:

1. Use class-specific references when converting `rho_hat` to weights.
2. Treat the output as an approximate assay field, not as gDNA-only evidence.
3. Apply the field symmetrically to gDNA and nRNA.

This is an incremental bridge, not the final exposure learner.

### 5.2 Longer-Term Field

Replace the class-specific density stitching with a smoothed genomic exposure
estimator. The estimator should transform all conservative evidence sources
onto one target:

```text
A(x) = relative assay accessibility at genomic position x
```

Possible evidence sources:

- contained intergenic fragments;
- contained intronic fragments;
- exon-intron boundary flux, projected onto local exon shoulders and nearby
  intronic windows;
- high-confidence gDNA posterior fragments from a preliminary uniform pass;
- optional target BED intervals in a later supervised mode.

The emitted field should be normalized robustly, for example with a high
exposure quantile or central on-target state set to one. Per-class diagnostics
should remain available, but downstream EM should consume only the single
`A(x)` field.

## 6. Candidate Modes To Test

Before choosing the final production default, run a controlled sweep on the same
VCaP input. These modes are diagnostic, not all intended for long-term CLI
surface.

### Mode 0: Uniform Baseline

Current `--regional-exposure off` behavior.

Purpose: same-version control.

### Mode 1: Current gDNA-Only Regional

Current v3 behavior.

Purpose: failed comparator; should reproduce the known gDNA-to-nRNA regression.

### Mode 2: Denominator-Only gDNA Regional

```text
gDNA numerator: unchanged
gDNA denominator: weighted L_gDNA^A
nRNA numerator: unchanged
nRNA denominator: unchanged
```

Purpose: isolate the harm caused by local numerator penalties. This mode should
not be considered fully generative, but it is conservative because it cannot
create local low-A penalties against gDNA relative to nRNA.

Expected result: gDNA-to-nRNA leak should be much lower than Mode 1. If this
happens, it confirms that the numerator asymmetry is load-bearing.

### Mode 3: Symmetric Unspliced Regional

```text
gDNA numerator: weighted by A
gDNA denominator: weighted L_gDNA^A
nRNA numerator: weighted by A
nRNA denominator: weighted L_nRNA^A
mRNA: unchanged
```

Purpose: production candidate for revised v1.

Expected result: removes the unfair gDNA-versus-nRNA local exposure penalty
while retaining exposure-aware opportunity normalization.

### Mode 4: Mean-Centered Numerator Diagnostic

```text
numerator factor = A(x) / Abar_component_or_locus
denominator returns to unweighted or cancels analytically
```

Purpose: diagnostic only. It tests whether within-locus exposure variation is
useful after removing average shifts. This is not the preferred production path
because it discards across-locus exposure information.

## 7. Implementation Sketch

### 7.1 Data Model

Keep `RegionalGdnaExposure` as the lookup table initially, but consider renaming
it later to `RegionalExposure` once the API is no longer gDNA-specific.

Required additions:

1. A function to compute weighted effective lengths for synthetic nRNA
   transcripts:

   ```python
   weighted_nrna_eff_lens(
       index: TranscriptIndex,
       rna_fl: FragmentLengthModel,
       exposure: RegionalGdnaExposure,
   ) -> np.ndarray
   ```

2. A function to apply per-unit exposure to nRNA candidates in the global CSR
   before `partition_and_free()`:

   ```python
   _apply_unit_unspliced_exposure_weights(
       em_data: ScoredFragments,
       exposure: RegionalGdnaExposure,
       index: TranscriptIndex,
       *,
       apply_gdna: bool,
       apply_nrna: bool,
   ) -> RegionalWeightApplicationStats
   ```

3. A way to pass exposed nRNA effective lengths into EM. The native EM already
   consumes one `t_eff_lens` array for all transcript-like components, so the
   Python `TranscriptGeometry.effective_lengths` can be adjusted before
   constructing `AbundanceEstimator` or before EM dispatch.

### 7.2 Effective Length Handling

Current EM uses:

```text
log_weights[c] = log(theta[c]) - log_eff_len[c]
```

for every transcript-like component and for gDNA. Therefore, symmetric nRNA
weighting requires the nRNA transcript effective lengths in `t_eff_lens` to be
the exposed lengths.

Proposed rule:

- annotated mature transcripts: keep existing RNA effective lengths;
- synthetic nRNA transcripts: replace with `L_nRNA^A(t)` when regional exposure
  is enabled;
- gDNA: keep per-locus `L_gDNA^A(M)`.

Uniform mode must remain bit-identical to the current path.

### 7.3 Numerator Application

The current scorer emits:

- per-candidate transcript log likelihoods in CSR arrays;
- one per-unit `gdna_log_liks` array;
- per-unit `genomic_midpoint`.

For gDNA, the current application path mutates `gdna_log_liks[u]`.

For nRNA, mutate only candidate entries whose transcript is synthetic nRNA and
whose unit is unspliced:

```text
em_data.log_liks[j] += log A(genomic_midpoint[u])
```

where candidate `j` belongs to unit `u` and `t_indices[j]` is synthetic nRNA.

Do not apply this to spliced units. Do not apply this to annotated mature mRNA
in v1.

Cross-reference behavior should match the current conservative gDNA behavior:

- if candidate transcripts in a unit span multiple references, skip the
  per-unit exposure term for that unit;
- count the skip in summary diagnostics;
- leave denominator computation unchanged because denominators are component
  specific and reference-aware.

### 7.4 Config Surface

Do not expose all experimental modes as permanent user-facing flags yet. Use
internal/debug modes while validating.

Production-facing target:

```text
--regional-exposure off
--regional-exposure auto
```

Where `auto` eventually means symmetric unspliced regional exposure if the
validation gate passes.

For development, add an internal enum or debug-only config value, for example:

```text
off
gdna_only
gdna_denominator_only
unspliced_symmetric
```

The CLI can keep exposing only `off` and `auto` until the model is validated.

## 8. Validation Plan

### 8.1 Oracle/Truth Confirmation

For the VCaP flowcell-truth run, continue using:

```text
C6EL5ANXX = RNA
H7MFFDSXY = gDNA
```

If an oracle BAM or true-origin tag is available for a simulated or controlled
subset, confirm that the implicated gDNA-to-nRNA transitions are truly gDNA
origin. The flowcell truth is already strong for this real mix, but an oracle
confirmation is useful for synthetic reproduction.

### 8.2 VCaP Mode Sweep

Run the same input BAM and seed for:

1. uniform/off;
2. current gDNA-only regional;
3. denominator-only gDNA regional;
4. symmetric unspliced regional.

Metrics:

- full truth-by-predicted pool confusion matrix;
- gDNA to mRNA;
- gDNA to nRNA;
- RNA to gDNA;
- transition matrix from uniform baseline to each regional mode;
- per-locus `gdna_eff_len_weight_ratio`;
- for false gDNA-to-nRNA fragments, `log A(x)`, `log Abar_locus`, and
  `log(A(x) / Abar_locus)`.

Acceptance criteria for the production candidate:

1. gDNA-to-nRNA must not exceed the same-version uniform baseline by a material
   amount.
2. RNA-to-gDNA should not regress materially relative to uniform.
3. Any residual gDNA-to-RNA movement should be concentrated in biologically
   plausible regions, not globally across low-A introns.
4. Synthetic uniform scenarios must remain near bit-identical or within known
   stochastic assignment variation.

### 8.3 Synthetic 24-Condition Validation

Rerun the existing synthetic 24-condition framework. Synthetic uniform data is
not sufficient to validate hybrid-capture behavior, but it is necessary to prove
the revised model does not regress ordinary use cases.

Acceptance criteria:

- all 24 conditions complete;
- strand-specific scenarios retain prior performance;
- no new nRNA siphon or gDNA leakage beyond known baseline behavior;
- `--regional-exposure off` remains exactly the same as current uniform.

### 8.4 Small Locus Simulation

Create a targeted synthetic locus with:

- one high-exposure exon;
- one low-exposure intron;
- true gDNA fragments in both regions;
- true nRNA fragments in both regions;
- equal or controlled underlying molecule abundance.

The expected behavior:

- gDNA-only numerator weighting misroutes low-A gDNA to nRNA;
- denominator-only reduces that local penalty;
- symmetric unspliced weighting preserves gDNA-vs-nRNA ratios across high-A and
  low-A positions, modulo component effective-length differences.

This test should be small enough to run in CI once implemented.

## 9. Risks And Open Questions

### 9.1 Does nRNA Effective Length Need RNA FL Or gDNA FL?

nRNA is RNA, but unspliced nRNA spans genomic sequence. Its fragment-length
distribution should use the RNA fragment-length model, not the gDNA model,
unless existing nRNA handling already intentionally uses another model. The
implementation must audit current nRNA effective-length construction before
patching.

### 9.2 What Is The Right Footprint For `A(u)`?

The current implementation uses a midpoint lookup. That is efficient and likely
adequate as a first approximation. A more exact model would integrate `A(x)`
over the fragment footprint.

Revised v1 should keep midpoint weighting for parity with the current gDNA path
unless validation shows boundary artifacts.

### 9.3 Should Mature mRNA Eventually Be Weighted?

Probably yes for a fully generative hybrid-capture model. But mRNA weighting is
not required to fix the unfair gDNA-to-nRNA competition and should be deferred
until the unspliced-symmetric model is validated.

### 9.4 Can Existing Class-Derived `A(x)` Be Trusted?

Only as a bridge. The long-term model should learn one physical exposure field
from smoothed assay evidence. The current class-derived field can still answer
the immediate question: does symmetric application eliminate the catastrophic
regression?

## 10. Recommended Path Forward

1. Freeze the current gDNA-only regional mode as experimental and unsafe for
   production.
2. Implement the class-reference fix only as hygiene, not as the main solution.
3. Add development modes for denominator-only and symmetric unspliced exposure.
4. Implement weighted effective lengths for synthetic nRNA components.
5. Apply midpoint `log A` to both gDNA and synthetic nRNA candidates for
   unspliced units.
6. Run the VCaP four-mode sweep and synthetic 24-condition validation.
7. If symmetric unspliced exposure passes, make it the behavior behind
   `--regional-exposure auto`.
8. Start a second design phase for a unified smoothed assay exposure field.

## 11. Bottom Line

The fair model is not "penalize intronic gDNA in hybrid capture." The fair model
is:

```text
Penalize every hypothesis whose molecule would be hard to observe from that
same intronic location.
```

For the current regression, the relevant competition is gDNA versus nRNA. Both
are unspliced genomic hypotheses. Therefore, regional exposure must be symmetric
across them. That is the conceptual overhaul this plan adopts.
