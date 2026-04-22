# gDNA Length Model — Problem Statement and Solution Options

**Date:** 2026-04-22
**Status:** Design document. **Option B selected for implementation.**
**Related:** [`gdna_locus_prior_deep_diagnostic.md`](gdna_locus_prior_deep_diagnostic.md)

---

## 0. Decision summary

After a second review (including a critique from Gemini), we pivot to
**Option B (per-fragment scalar gDNA effective length via harmonic
mean of per-Locus lengths)** as the chosen implementation. Option A
(per-hit gDNA candidates) was initially recommended but is rejected
for two substantive reasons:

1. **VBEM sparsification hazard.** Option A introduces K independent
   gDNA parameters per MultiLocus (one per contiguous region).
   Rigel's VBEM digamma penalty will prune low-count candidates to
   zero. gDNA, as a single biophysical contamination process, must
   not be fragmented across disjoint parameters — any region that
   gets pruned can no longer accept gDNA posterior, which breaks the
   gDNA↔RNA competition in precisely the regions that need it most.
2. **Option B is the closed-form marginal, not an approximation.**
   Under a uniform-gDNA-over-compatible-positions prior, the exact
   per-fragment likelihood is $p(F\mid \text{gDNA}) = \sum_h (1/\text{NH})\cdot (1/e_h)$.
   Baking this into a single per-fragment scalar via the harmonic
   mean of per-hit $e_h$ values preserves unified $\theta_{\text{gDNA}}$
   statistical strength and is structurally identical in likelihood
   terms to Option A, without Option A's parameter fragmentation.

Option A's information about per-region gDNA density variation (which
is real, due to CNV) is better recovered under **Option C
(hierarchical per-region rates)** as a future workstream, not under
a regional EM that overfits noise.

**Implementation plan:** see Section 6.

---

---

## 1. The problem

### 1.1 What Rigel does today

Rigel partitions the quantification problem into **loci**: connected
components of transcripts linked by **shared fragments**. A transcript
$t_i$ and a transcript $t_j$ land in the same locus if any fragment has
compatible alignment hits on both of them (directly, or transitively
through other transcripts/fragments). This is a pure graph operation;
union-find over `(transcript, transcript)` edges induced by the
multi-alignment structure of the BAM.

Inside each locus, the EM solves for posteriors across a set of
candidates. Those candidates come in three flavors:

- **mRNA transcript candidates** — one per annotated transcript in
  the locus. Candidate $t$ has `eff_len = spliced_length(t)`.
- **nRNA transcript candidates** — one per synthetic nRNA entity
  in the locus. Candidate $n$ has `eff_len = genomic_span(n)`
  (single-exon synthetic transcript).
- **a single gDNA candidate per locus** — `eff_len = gdna_span + gdna_flank`,
  where `gdna_span` is the *merged genomic footprint* of every
  transcript in the locus (intron-inclusive).

The per-candidate likelihood for a fragment $F$ is then adjusted in-place
by `log_liks[i] -= log(eff_len - frag_len + 1)` before EM iteration
(see [`em_solver.cpp::apply_bias_correction_uniform`](../../src/rigel/native/em_solver.cpp#L330)).
This is the uniform-position prior familiar from Salmon-style isoform
quantification.

### 1.2 The asymmetry

| Aspect | mRNA | nRNA | gDNA (current) |
|---|---|---|---|
| # candidates per fragment in locus | one per compatible transcript hit | one per compatible nRNA hit | **one per locus** |
| Each candidate's `eff_len` | spliced length of *that* transcript | genomic span of *that* nRNA | **summed genomic span across ALL regions in the component** |
| Locally defined? | Yes | Yes | **No — globally inflated by union-find** |

RNA candidates are **local to the hypothesis**: "if this fragment came
from transcript $T$, its support set (positions it could have originated
from) is the transcript body itself." The eff_len is a property of the
hypothesis alone.

gDNA is different. One gDNA "hypothesis" per locus currently asserts
"if this fragment came from gDNA, its support set is the entire
union of every genomic region any multimapper touches in this
component." **This is true only if a gDNA fragment could have originated
anywhere in the locus with equal probability — which is false whenever
multimappers link physically disjoint regions.**

### 1.3 The multimapper explosion

Consider a single fragment with NH = 10 alignments landing on 10
different genomic regions. That fragment alone yokes 10 regions into
a single connected component (because its presence forces all 10
transcripts to be jointly solved for). After union-find:

- `locus.gdna_span` now includes **all 10 regions** summed.
- **Every other fragment** in the component — including fragments
  that have only a single, unique alignment in one of those 10
  regions — now has its gDNA likelihood divided by the 10-region
  total.

Under the current model, a single multimapper with NH=10 effectively
shrinks the gDNA per-position probability by ~10× **for every other
fragment in the component**, including the unique mappers. One
hypothesis's support-set width is being imposed on a completely
different hypothesis's fragment.

On the human genome with permissive aligners (minimap2 `-N 20`),
this cascades. In the dna80m benchmark we observe a **569 Mb**
connected component spanning large chunks of multiple chromosomes,
containing 10.5 M EM fragments (41% of the EM pool). Per-fragment
gDNA log-likelihood deficit vs typical RNA is ~12 nats. The gDNA
posterior inside that mega-locus collapses to 0.50 even though the
true contamination fraction is 0.80 and the locus-local prior
$\gamma_\ell \approx 0.28$ warm-starts it reasonably.

### 1.4 What is the *right* answer?

The Bayesian question is: *"if a fragment $F$ came from gDNA, what is
the probability of observing its alignments?"* Under a uniform-over-
positions model of gDNA origin within the genomic support set:

$$
p(F \mid \text{gDNA}) \;=\; \sum_{h=1}^{\text{NH}} \frac{1}{\text{NH}} \cdot \frac{1}{L_h - \ell_F + 1}
$$

where:

- $L_h$ is the length of the **local** genomic support region
  containing hit $h$ (a calibration region, a contiguous mappable
  interval, whichever convention we pick),
- $\ell_F$ is the fragment's genomic footprint,
- $1/\text{NH}$ distributes the gDNA hypothesis uniformly across
  the fragment's hit positions.

The essential property: **$L_h$ is a property of hit $h$ itself, not
of the locus $h$ happens to end up in.** Two fragments that each have
unique hits in the same 50 kb region should see the same
$L = 50\text{ kb}$ in their gDNA likelihood, regardless of whether
they ended up in a 50 kb locus, a 500 kb locus, or a 569 Mb
mega-locus.

The locus-level eff_len = sum-over-regions model **conflates**:

- the fragment's actual alignment footprint (what we can observe),
- the connectivity of the transcript graph (purely an artifact of
  how other fragments happened to align).

These are independent things. Only the first should enter the
per-fragment likelihood.

---

## 2. Solution options

### 2.1 Option A — Per-hit gDNA candidates (recommended)

#### 2.1.1 The core idea

**Treat each alignment hit of a fragment as potentially originating
from gDNA at that specific location**, rather than treating "gDNA"
as a single locus-wide abstract state.

Today, when a fragment has NH=10 hits, we emit:
- Up to 10 **transcript candidates** (one per hit that resolves to
  an annotated transcript; RNA is already per-hit),
- **1 gDNA candidate** — a single state shared by all hits.

Option A emits:
- Up to 10 transcript candidates (unchanged),
- **Up to 10 gDNA candidates**, one per alignment hit, each with
  its own local eff_len $L_h$ (e.g., the calibration region
  containing that hit).

#### 2.1.2 Are we "breaking apart" multimapping fragments?

**No.** A multimapping fragment remains a single fragment with a
single total posterior mass of 1.0. We are not assigning the fragment
to multiple places and summing duplicates.

What we are doing is **enumerating all alternative origin hypotheses
and letting the EM weight them**. Today a fragment with NH=10 has
10 RNA hypotheses and 1 collapsed gDNA hypothesis. Option A gives it
10 RNA hypotheses and 10 gDNA hypotheses (one per hit). The
fragment's posterior still sums to 1. The EM now has 20 candidates
competing for that 1.0 instead of 11, and each gDNA candidate is
spatially specific.

Concretely, a fragment with NH=3 hits in regions $R_1, R_2, R_3$
would contribute candidates (with $\ell$ = fragment length):

| Candidate | component | eff_len | log_lik scale |
|---|---|---|---|
| $t_A$ at hit 1 | mRNA | $L_{t_A}$ | $\log \tfrac{1}{L_{t_A}-\ell+1}$ |
| $t_B$ at hit 2 | mRNA | $L_{t_B}$ | $\log \tfrac{1}{L_{t_B}-\ell+1}$ |
| $t_C$ at hit 3 | mRNA | $L_{t_C}$ | $\log \tfrac{1}{L_{t_C}-\ell+1}$ |
| gDNA at hit 1 | gDNA | $L_1$ | $\log \tfrac{1}{L_1-\ell+1}$ |
| gDNA at hit 2 | gDNA | $L_2$ | $\log \tfrac{1}{L_2-\ell+1}$ |
| gDNA at hit 3 | gDNA | $L_3$ | $\log \tfrac{1}{L_3-\ell+1}$ |

All six candidates carry the fragment's 1/NH = 1/3 alignment-uniform
prior already embedded in the scoring stage, plus fragment-length and
splice-compatibility terms identical to current practice.

#### 2.1.3 Why this is the same pattern as RNA

For RNA, a fragment with multiple hits gets one RNA candidate **per
hit** because each hit is a physically distinct alignment position on
a specific transcript. The EM uses the transcript's length as that
candidate's support set because the transcript is the biological unit
that the fragment supposedly originated from.

The gDNA analogue is: each hit is a physically distinct alignment
position on the **genome**; the local genomic region containing that
hit is the biological unit the fragment supposedly originated from,
so the local region's length is that candidate's support set.

Rigel's current per-locus collapse breaks this symmetry by erasing
the spatial specificity of each hit. Option A restores it.

#### 2.1.4 What "local region" means

Candidate definitions of $L_h$ for hit $h$:

1. **The calibration region containing $h$** — we already index
   `regions.feather` with `region_cr` (cgranges). A single overlap
   query per hit gives us the region, and we have both `length` and
   `mappable_effective_length` available. Using
   `mappable_effective_length` simultaneously resolves the P2
   mappable-bp/genomic-bp asymmetry noted in the diagnostic doc.
2. **A local contiguous mappable interval** — query alignable for a
   window around the hit and compute mappable mass locally. More
   physically precise but more expensive; probably unnecessary given
   that calibration regions are already tuned to be locally uniform.

Recommendation: **start with (1), use `mappable_effective_length`**.
If a hit happens to span a calibration region boundary, use the
region with greater overlap or sum the contribution of both regions'
lengths (weighted by overlap).

#### 2.1.5 Interaction with mRNA and nRNA

**mRNA is unchanged.** Each mRNA candidate already carries a
transcript-local eff_len (spliced length). The fact that the locus
now contains more gDNA candidates does not change any RNA candidate's
likelihood.

**nRNA is unchanged at the candidate level** but its competitive
position changes. Today nRNA is a frequent siphon for misassigned
gDNA fragments because:
- nRNA's eff_len ≈ pre-mRNA genomic span (10–100 kb typically),
- gDNA's eff_len ≈ locus span (up to 569 Mb),
- per-position probability of nRNA » per-position probability of gDNA,
- fragments that should go to gDNA end up at nRNA.

Under Option A, gDNA candidates now have eff_len $L_h$ comparable to
nRNA (both are local genomic regions, same order of magnitude). The
competition becomes **physically fair** — both hypotheses are scored
against roughly-equal-sized genomic supports, and the prior
($\gamma_\ell$ for gDNA, separate for nRNA) drives the final
allocation.

This is the strongest argument for A over alternatives: the mega-locus
pathology and the nRNA siphon are **the same underlying bug**, and
A fixes both simultaneously.

#### 2.1.6 Locus-level rollups

`get_loci_df` today reports `gdna` as the EM-summed posterior to the
single gDNA state. Under Option A:

- Per-locus `gdna` = sum over all hits whose containing region is
  inside the locus of the per-hit gDNA posteriors.
- Natural definition of "the locus's gDNA abundance" because each
  hit has a specific genomic location, and the locus is still a
  well-defined set of regions.
- Cross-chromosome mega-locus: `gdna` splits naturally across its
  sub-regions without any special handling.

#### 2.1.7 Cost

- **Candidate count:** per fragment, roughly doubles (NH gDNA
  candidates instead of 1). In the mega-locus (10.5 M fragments) this
  adds ~10.5 M × mean_NH extra rows. Mean NH in the EM pool is
  typically 2–3, so this is ~20–30 M rows — tractable.
- **Memory:** each extra candidate is a few floats + int. Peak EM
  memory grows proportional to mean NH on gDNA-heavy data.
- **CPU:** EM per-iteration cost scales linearly with candidate
  count. Expect a ~20–30% increase. Minor.

#### 2.1.8 Implementation sketch

1. Scoring stage (`scoring.py` / `scoring.cpp`): for each fragment
   hit, emit a gDNA candidate with `t_index = GDNA_SENTINEL` and
   `bias_profile_hint = L_h`. Per-hit splice penalty and fragment-
   length terms identical to current gDNA scoring.
2. Partitioning: gDNA candidates follow the same locus assignment as
   their hit's transcript, so they land in the same locus as the
   RNA candidates they're competing with.
3. Bias correction: remove the special-case `bias_profiles[gdna_idx]`
   path; gDNA candidates now have their own `bias_profiles[i]`
   entries populated by the scorer.
4. EM: no change to the solver itself beyond recognizing that each
   fragment can now have multiple gDNA rows. The E-step already
   handles arbitrary candidate layouts.
5. Rollups: sum gDNA posteriors per locus as today, but with the
   understanding that locus is defined by spatial regions.

---

### 2.2 Option B — Per-fragment scalar gDNA eff_len

Keep one gDNA row per fragment per locus (one candidate), but make its
eff_len a **per-fragment scalar** computed at scoring time:

$$
\text{eff\_len}_F \;=\; \frac{\text{NH}}{\sum_h 1/L_h}
$$

This is the harmonic-mean collapse of Option A. Each fragment's gDNA
likelihood gets a scalar baked into `gdna_log_liks[F]` before EM,
and the per-locus `bias_profiles[gdna_idx]` mechanism is deprecated
or ignored for gDNA.

**Pros:** smallest code footprint; no extra candidates; directly
addresses the locus-size blowup.

**Cons:** loses the per-hit spatial specificity. When two gDNA hits
in different sub-regions of a mega-locus have very different
densities (e.g., one in a repeat region, one in a clean region),
the collapsed scalar averages them. Option A handles that
heterogeneity naturally.

Useful as a **minimal first patch** to validate the direction before
committing to Option A's candidate structure overhaul.

---

### 2.3 Option C — Genome-wide joint EM with a single gDNA component

#### 2.3.1 The idea

Abandon locus partitioning for the gDNA component only. Solve one
giant EM that includes:

- All mRNA candidates genome-wide (transcript isoform EM),
- All nRNA candidates genome-wide,
- **One gDNA component** that accepts fragments from anywhere on the
  genome, including intergenic fragments that currently bypass EM
  entirely,
- **One intergenic rate parameter** (or equivalently $\lambda_G$
  estimated jointly with everything else).

This extends the analogy you drew: *"the way we fold spliced reads
into RNA components, we fold intergenic reads into the gDNA
component."* The EM would concurrently estimate isoform abundances
and a global gDNA density parameter.

As you noted, in the current data the mega-locus is *already*
approximately a genome-wide EM — it has 41% of the EM fragments,
70k transcripts, and spans 569 Mb. Option C formalizes this.

#### 2.3.2 Why it's attractive

- **Clean semantics.** gDNA is a global property; $\lambda_G$ has
  direct physical meaning; no per-locus γ aggregation needed.
- **No locus-boundary artifacts.** The P1 pathology disappears because
  there is no "locus denominator" for gDNA.
- **Intergenic fragments contribute.** Currently they're accumulated
  into per-region counts for calibration but never pass through EM.
  Giving them a gDNA candidate lets them influence $\lambda_G$
  consistently with everything else.

#### 2.3.3 Why the naïve version is wrong

A single global gDNA rate $\lambda_G$ assumes gDNA contamination is
**spatially uniform** across the genome. This is biologically false.
Sources of true per-locus gDNA variation include:

- **Copy-number amplifications** (CNAs): a cancer cell line with a
  10× amplification of chr8 has 10× more genomic DNA per cell in
  that region. Mapped to the reference, this projects as 10× gDNA
  coverage in that locus — not a factor-of-1 global rate.
- **Chromosome-level amplifications and polyploidy:** common in
  tumor samples.
- **Deletions and losses of heterozygosity:** the opposite — some
  regions have zero or reduced gDNA.
- **Technical biases:** GC bias, fragmentation bias, capture bias in
  targeted libraries all create spatially non-uniform gDNA
  depositions.

**This is the same identifiability problem that intratumor
heterogeneity poses to CNV callers.** The EM cannot distinguish
"spatially non-uniform gDNA" from "spatially non-uniform RNA" without
some structural constraint.

#### 2.3.4 Implication

The model needs to **allow per-locus (or per-region) gDNA rate
variation** while still respecting the global rate structure. Naïve
Option C fights the data. A more sophisticated version would:

- Use $\lambda_G$ as a **prior mean** on per-region gDNA density, not
  a hard constraint.
- Allow data-driven per-region deviation from the prior (hierarchical
  Bayesian gDNA density) — exactly what Rigel's per-locus $\gamma_\ell$
  already does for the prior. The difference is that Option C would
  have the EM itself update per-region gDNA density, rather than
  fixing it from calibration.

#### 2.3.5 Bottom line on C

A correctly-formulated global-gDNA EM is **not incompatible with
Option A**. In fact, Option A already approximates this: each
per-hit gDNA candidate has $L_h = $ local region length, which is
the per-region support. What Option A doesn't include is
per-region rate variation beyond the locus-level prior. That's a
future extension.

---

## 3. Recommended path (revised)

**Option B is the implementation target.**

1. **Implement Option B directly** (no separate prototype step).
   Per-fragment scalar gDNA eff_len $L_{\text{eff},F} = \text{NH}/\sum_h (1/e_h)$
   where $e_h = \max(L_h + \text{gdna\_flank} - \ell_F + 1, 1)$ and
   $L_h$ is the length of the **contiguous merged Locus interval**
   that hit $h$ lands in. Bake into `gdna_log_liks[F]` at scoring /
   partitioning time. Disable the current per-locus gDNA bias
   correction inside the C++ EM.
2. **Option A is rejected.** (See Section 0 decision summary.)
3. **Option C (hierarchical per-region gDNA rates) is a future
   workstream** to address CNV-driven spatial variation properly.
   Not compatible with Option B's "one gDNA candidate per MultiLocus"
   structure but can extend it: Option C would generalize the
   locus-level $\gamma_\ell$ prior to a per-Locus rate parameter
   that the EM updates jointly.

---

## 4. Summary table

| Aspect | A (per-hit candidates) | **B (per-frag scalar)** ✅ | C (global EM) |
|---|---|---|---|
| Fixes mega-locus pathology | ✅ | **✅** | ✅ |
| Fixes nRNA siphon | ✅ | **✅** | ✅ |
| Safe under VBEM pruning | ❌ (fragments $\theta$) | **✅** (unified $\theta$) | ✅ |
| Exact marginal likelihood | ✅ | **✅** (closed-form collapse) | depends |
| Handles CNV variation | nominal (overfits noise) | via prior only | via hierarchical extension |
| Symmetric with RNA model | natively | **✅** (same bias-correction pathway) | rearchitected |
| Code footprint | medium | **small** | large |
| Candidate count growth | ~2× | **none** | none |
| CPU cost increase | ~20–30% | **negligible** | architecture-dependent |

---

## 5. Open questions (resolved for Option B)

- **Boundary cases for $L_h$:** $L_h$ = length of the contiguous
  merged Locus interval (entry in `locus.merged_intervals`) that
  contains hit $h$. Not calibration regions. This is the natural
  "region of the genome where this fragment can compete" definition.
- **Intergenic hits:** remain excluded from the EM. Folding them in
  requires Option C's architectural changes. Option B solves the
  mega-locus problem without touching the routing logic.
- **gDNA prior inside the locus:** unchanged — exactly one gDNA
  candidate per MultiLocus, one $\gamma_\ell$ per MultiLocus, Jeffreys
  correction untouched.
- **Flank treatment:** applied per Locus interval, not per MultiLocus.
  Each merged interval has two gDNA-accessible boundaries so the
  flank is added to each interval's length individually.

---

## 6. Implementation plan for Option B — **FINAL (in-scorer harmonic mean)**

### 6.0 Guiding principle

Compute the exact Option B harmonic-mean length correction **inside
the scorer**, at the same point where every other per-fragment
physical quantity (fragment-length density, strand LLR, nm penalty)
is computed. One pass, no back-references to locus topology, no
per-hit state propagated to Python.

The scorer already iterates every multimapper alignment inside
`score_mm_alignment` and has access to `t_span_[t_idx]` — the
transcript's genomic span (including introns) — for free. Adding the
harmonic-mean accumulator is a per-hit ALU op. The resulting
`gdna_log_liks[u]` array is emitted **already carrying the full
Option B length correction**; everything downstream (buffer,
partition, EM) is ignorant of gDNA length physics entirely.

This is a net-deletion refactor:

| Deleted | Location |
|---|---|
| `gdna_flank` threading through pipeline | `pipeline.py` (3 sites) |
| `gdna_spans` parameter in batch EM wrapper | `pipeline.py`, `estimator.py` |
| `gdna_spans` C++ kwarg and plumbing | `em_solver.cpp` (5 sites) |
| `gdna_span` parameter in `extract_locus_sub_problem_from_partition` | `em_solver.cpp` |
| `sub.bias_profiles[sub.gdna_idx] = gdna_span` | `em_solver.cpp` |
| `genomic_footprints` per-unit array (entire pipeline) | `scored_fragments.py`, `scan.py`, `partition.py`, `pipeline.py`, `em_solver.cpp`, `scoring.cpp` |
| `mm_best_gdna_fp`, `mm_first_gfp` accumulators | `scoring.cpp` |

The only new data on the wire is `gdna_flank` passed once to the
scorer constructor alongside `gdna_fl_log_prob_obj` it is already
paired with.

### 6.1 Mathematical specification

For each EM unit $u$ built from a fragment with multimapper hits
$\mathcal{H}(u)$ (all unspliced hits; spliced hits exclude the
fragment from the gDNA hypothesis as today):

- Per hit $h \in \mathcal{H}(u)$ with transcript $t_h$ and genomic
  footprint $\ell_u$ (footprint is a fragment property, constant
  across hits):
  $$L_h = \texttt{t\_span\_}[t_h] + \texttt{gdna\_flank}$$
  $$e_h = \max(L_h - \ell_u + 1,\; 1)$$
  $$\log p_h^{\text{gDNA}} = \texttt{gdna\_fl}(\ell_u) + \texttt{gdna\_log\_sp} + \log\tfrac{1}{2} + \texttt{nm}_h \cdot \texttt{mm\_log\_pen} - \log e_h$$
- Fragment-level gDNA log-likelihood (marginalizing uniformly over
  the NH equally-likely alignment choices):
  $$\texttt{gdna\_log\_liks}[u] = -\log|\mathcal{H}(u)| + \operatorname{logsumexp}_{h \in \mathcal{H}(u)} \log p_h^{\text{gDNA}}$$
- If any hit of the fragment is spliced: $\texttt{gdna\_log\_liks}[u] = -\infty$
  (unchanged behavior).
- Unique-mapper fallthrough: $|\mathcal{H}(u)| = 1$, reduces to
  $\log p_1^{\text{gDNA}}$.

Interpretation: `t_span_[t_h]` is the local genomic window
consistent with an alignment at hit $h$. The harmonic mean is the
exact closed form of $\log p(\text{fragment} \mid \text{gDNA})$
under a uniform-position prior restricted to the union of the
fragment's NH plausible origins, marginalized over equally likely
hit selection.

**Why `t_span_` and not merged-interval length:** `t_span_` answers
the right question per-hit (*"how big is the gDNA window that
contains this specific alignment?"*). Merged-interval length couples
transcripts whose only connection is a third-party multimapper
elsewhere, re-introducing the very coupling we are trying to break.
For overlapping transcript clusters, `t_span_` values differ
negligibly; for paralog/pseudogene multimappers, `t_span_`
per-hit is strictly more informative than one lumped interval.

### 6.2 C++ scorer changes — the core of the fix

#### 6.2.1 `FillState` accumulator

Replace the scalar `mm_best_gdna_ll` / `mm_best_gdna_fp` /
`mm_first_gfp` triplet with a numerically stable online
logsumexp accumulator:

```cpp
// In struct FillState:
// -----  REMOVED  -----
// double  mm_best_gdna_ll;
// int32_t mm_best_gdna_fp;
// int32_t mm_first_gfp;

// -----  ADDED  -----
double  mm_gdna_lse_max;   // running max of per-hit gDNA log-liks
double  mm_gdna_lse_sum;   // Σ exp(x_h − mm_gdna_lse_max)
int32_t mm_nh_gdna;        // # of unspliced hits contributing
```

`reset_mm_group()`:
```cpp
mm_gdna_lse_max = -std::numeric_limits<double>::infinity();
mm_gdna_lse_sum = 0.0;
mm_nh_gdna = 0;
```

Helper (static, file-local):
```cpp
static inline void lse_update(double& max_v, double& sum_v, double x) {
    if (x == -std::numeric_limits<double>::infinity()) return;
    if (x > max_v) {
        if (max_v > -std::numeric_limits<double>::infinity())
            sum_v *= std::exp(max_v - x);
        max_v = x;
        sum_v += 1.0;
    } else {
        sum_v += std::exp(x - max_v);
    }
}
```

#### 6.2.2 `score_mm_alignment` — per-hit contribution

Replace the existing gDNA block (lines ~478–495) with:

```cpp
// gDNA: accumulate per-hit log-lik with Option B length correction.
// Only unspliced hits enter the gDNA hypothesis.
if (stype != SPLICE_SPLICED_ANNOT &&
    stype != SPLICE_SPLICED_UNANNOT)
{
    int32_t gfp_val = cp.g_fp[row];
    double gdna_fl = gdna_frag_len_log_lik(gfp_val);

    // Harmonic-mean length term: pick the "best" RNA candidate's
    // transcript as the spatial anchor for this hit (all candidates
    // at a given hit overlap the same genomic region, so their
    // t_span values are near-identical; we use the first one for
    // speed — no sensitivity to this choice).
    int32_t anchor_t = (n_cand > 0) ? cp.t_ind[start] : -1;
    int32_t t_span   = (anchor_t >= 0) ? t_span_[anchor_t] : 0;
    int64_t L_h      = static_cast<int64_t>(t_span) + gdna_flank_;
    int64_t e_h      = L_h - static_cast<int64_t>(gfp_val) + 1;
    if (e_h < 1) e_h = 1;

    double hit_log_lik = gdna_fl + gdna_log_sp + LOG_HALF + log_nm
                       - std::log(static_cast<double>(e_h));

    lse_update(st.mm_gdna_lse_max, st.mm_gdna_lse_sum, hit_log_lik);
    ++st.mm_nh_gdna;
}
```

Note: `anchor_t` is the first RNA candidate of the current hit
(from `cp.t_ind[cp.t_off[row]]`). For MM hits with zero candidates
(unresolved), `anchor_t = -1` and the hit is skipped — consistent
with its inability to contribute to gDNA sampling.

#### 6.2.3 `flush_mm_group` — emit one scalar per unit

Replace the gdna_ll emission block (lines ~646–652) with:

```cpp
double final_gdna_ll;
if (st.mm_is_any_spliced || st.mm_nh_gdna == 0) {
    final_gdna_ll = NEG_INF;
} else {
    final_gdna_ll = st.mm_gdna_lse_max
                  + std::log(st.mm_gdna_lse_sum)
                  - std::log(static_cast<double>(st.mm_nh_gdna));
}
st.v_gdna_ll->push_back(final_gdna_ll);
```

The `v_gfp->push_back(...)` lines are **deleted** (array no longer
exists; see 6.4).

#### 6.2.4 Non-MM path — NH=1 special case

Replace the non-MM gdna_ll block (lines ~902–915) with:

```cpp
bool is_spl = (stype == SPLICE_SPLICED_ANNOT
            || stype == SPLICE_SPLICED_UNANNOT);
st.v_is_spliced->push_back(is_spl ? 1 : 0);

if (!is_spl && best_t >= 0) {
    double gdna_fl = gdna_frag_len_log_lik(genomic_footprint);
    int64_t L   = static_cast<int64_t>(t_span_[best_t]) + gdna_flank_;
    int64_t e_h = L - static_cast<int64_t>(genomic_footprint) + 1;
    if (e_h < 1) e_h = 1;
    st.v_gdna_ll->push_back(
        gdna_fl + gdna_log_sp + LOG_HALF + log_nm
      - std::log(static_cast<double>(e_h)));
} else {
    st.v_gdna_ll->push_back(NEG_INF);
}
```

Here `best_t` is the already-computed mRNA-winner for this
fragment; for NH=1 that's the correct (and only) anchor.

#### 6.2.5 Scorer constructor — accept `gdna_flank`

Add one scalar to `NativeFragmentScorer`:

```cpp
// In member section:
int32_t gdna_flank_;     // NEW

// In constructor signature — add after gdna_fl_tail_base:
int32_t gdna_flank,

// In initializer list:
gdna_flank_(gdna_flank),
```

Nanobind binding: add one int kwarg to the constructor. All other
scorer state is unchanged.

### 6.3 `apply_bias_correction_uniform` guard

```cpp
// In em_solver.cpp
static void apply_bias_correction_uniform(
    double*        log_liks,
    const int32_t* t_indices,
    const int32_t* tx_starts,
    const int32_t* tx_ends,
    const int64_t* profile_lengths,   // indexed [0..n_t)
    int32_t        n_t,                // NEW
    size_t         n_candidates)
{
    for (size_t i = 0; i < n_candidates; ++i) {
        int32_t local = t_indices[i];
        if (local >= n_t) continue;          // gDNA row — pre-corrected
        // ... existing body unchanged ...
    }
}
```

One new call-site kwarg (`n_t`) added at each of the two
invocations (lines ~1142 and ~2105). The `bias_profiles` vector is
now sized `n_t` instead of `nc`; `sub.bias_profiles.resize(n_t)`.

### 6.4 `genomic_footprints` / `v_gfp` deletion

With no downstream consumer of per-unit genomic footprint (it was
used solely to size the gDNA row's bias correction), the entire
plumbing is removed:

**Files and exact sites:**
- `src/rigel/scored_fragments.py`: remove `genomic_footprints`
  field from `ScoredFragments`, `ScoredFragmentsView`, and
  `LocusPartition` dataclasses (3 sites).
- `src/rigel/scan.py` line 157, 226: drop from `StreamingScorer`
  output tuple and `ScoredFragments` constructor call.
- `src/rigel/partition.py` line 75, 110: drop from scatter map
  and `LocusPartition` constructor call.
- `src/rigel/pipeline.py` line 558: drop from `_call_batch_em`
  tuple; shift subsequent tuple indices (see 6.5.1).
- `src/rigel/native/em_solver.cpp`:
  - `PartitionView::genomic_footprints` field removed (line 1720).
  - `pv.genomic_footprints` binding read removed (line 1971);
    subsequent `tup[N]` indices shifted.
  - `extract_locus_sub_problem_from_partition`: the
    `sort_buf[k++] = {sub.gdna_idx, gdna_ll, 1.0, 0, footprint, 0}`
    line (1832) becomes
    `sort_buf[k++] = {sub.gdna_idx, gdna_ll, 1.0, 0, 0, 0}` —
    `tx_start` and `tx_end` are ignored by the gDNA-guarded bias
    correction (6.3).
- `src/rigel/native/scoring.cpp`: remove `v_gfp` from `FillState`
  (field + constructor + push_back sites in both MM and non-MM
  flush paths). Remove `mm_best_gdna_fp` and `mm_first_gfp`.

This is the "ripple" that scales the net deletion from "several
lines" to "a clean, visible improvement" in diff volume.

### 6.5 Pipeline + estimator cleanup

#### 6.5.1 `src/rigel/pipeline.py`

- Line ~815: `gdna_flank = int(calibration.gdna_fl_model.mean)`
  moves to a site where the `NativeFragmentScorer` is
  constructed. Grep for that construction (likely in
  `scoring.py::FragmentScorer.__init__`); pass `gdna_flank` as
  a new constructor kwarg.
- Lines 604, 639: delete `locus.gdna_span + gdna_flank` span
  construction.
- Line 824: delete `gdna_flank=gdna_flank` kwarg from
  `_run_locus_em_partitioned` call.
- `_call_batch_em`: drop the `batch_spans` parameter; drop the
  `p.genomic_footprints` slot from `partition_tuples`. Re-count
  tuple slots consumed on the C++ side.
- `_run_locus_em_partitioned` signature: drop `gdna_flank`
  parameter.

#### 6.5.2 `src/rigel/estimator.py`

- `run_batch_locus_em_partitioned` wrapper: drop `gdna_spans`
  parameter; drop from the `_em_impl.batch_locus_em_partitioned`
  call signature.

#### 6.5.3 `src/rigel/scoring.py`

- `FragmentScorer.__init__`: accept `gdna_flank: int`; pass to
  `_scoring_impl.NativeFragmentScorer(...)` binding.
- Callers of `FragmentScorer` — grep `FragmentScorer(` — get
  `gdna_flank` plumbed from the same calibration source.

### 6.6 `src/rigel/native/em_solver.cpp` signature changes

```cpp
// extract_locus_sub_problem_from_partition — drop gdna_span
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,
    double alpha_gdna,
    // int64_t gdna_span,           <-- REMOVE
    const double* unambig_row_sums,
    const int32_t* all_t_lengths,
    const int32_t* local_map,
    int32_t local_map_size);

// batch_locus_em_partitioned — drop gdna_spans kwarg
// (the `i64_1d gdna_spans` parameter at line ~1930 is removed)
// Callers: estimator.py binding; pipeline.py tuple construction.
```

### 6.7 Touch-point summary

| File | Net delta | Reason |
|---|---:|---|
| `src/rigel/native/scoring.cpp` | +~25, −~15 | logsumexp accumulator; per-hit length correction; kill `v_gfp`/mm_best_fp |
| `src/rigel/native/em_solver.cpp` | +~3, −~25 | `n_t` guard; drop `gdna_span` / `gdna_spans` / `genomic_footprints` / `bias_profiles[gdna_idx]` |
| `src/rigel/scoring.py` | +~3 | accept + forward `gdna_flank` |
| `src/rigel/scan.py` | −~3 | drop `genomic_footprints` from output |
| `src/rigel/partition.py` | −~3 | drop `genomic_footprints` scatter |
| `src/rigel/scored_fragments.py` | −~6 | drop field from 3 dataclasses |
| `src/rigel/pipeline.py` | +~2, −~12 | wire `gdna_flank` into scorer; drop batch `gdna_spans` / per-call span construction / tuple slot |
| `src/rigel/estimator.py` | −~4 | drop `gdna_spans` wrapper kwarg |
| **Net** | **≈ −35 LOC** | |

Plus: ~60 LOC of new unit tests and regression assertions.

### 6.8 Rollout sequence (implementation order)

**Phase 1 — Scorer-side physics (self-contained, compiles alone):**
1. `scoring.cpp`: add `gdna_flank_` member, constructor param, nanobind kwarg.
2. `scoring.cpp`: introduce `lse_update` helper; add logsumexp accumulator state; delete `mm_best_gdna_ll`/`mm_best_gdna_fp`/`mm_first_gfp`.
3. `scoring.cpp`: rewrite per-hit gDNA block in `score_mm_alignment`.
4. `scoring.cpp`: rewrite flush emission in `flush_mm_group`.
5. `scoring.cpp`: rewrite non-MM gDNA emission.
6. `scoring.py`: forward `gdna_flank` to the scorer constructor.
7. `pipeline.py`: pass `int(calibration.gdna_fl_model.mean)` into scorer construction.
8. **Compile + run full test suite.** Expect many golden failures; do NOT regenerate yet.

**Phase 2 — EM solver simplification:**
9. `em_solver.cpp`: drop `gdna_span` param from extraction; resize `bias_profiles` to `n_t`; delete `bias_profiles[gdna_idx] = gdna_span`.
10. `em_solver.cpp`: add `n_t` param to `apply_bias_correction_uniform`; update both call sites.
11. `em_solver.cpp`: drop `gdna_spans` kwarg from `batch_locus_em_partitioned` binding; drop `gs_ptr` use.
12. `estimator.py`: drop `gdna_spans` wrapper kwarg.
13. `pipeline.py`: drop `batch_spans` from `_call_batch_em`; drop `np.array([locus.gdna_span + gdna_flank], ...)` at both call sites; drop `gdna_flank` param from `_run_locus_em_partitioned`.
14. **Compile + test.** Golden failures persist (output is changing); still do not regenerate.

**Phase 3 — `genomic_footprints` rip-out:**
15. `scoring.cpp`: delete `v_gfp` field + all push sites.
16. `scan.py`, `scored_fragments.py`, `partition.py`, `pipeline.py`: drop `genomic_footprints` from output tuple / dataclass / scatter / C++ input tuple.
17. `em_solver.cpp`: drop `genomic_footprints` from `PartitionView` + binding; adjust `sort_buf` gDNA entry to use zero footprint.
18. **Compile + test.**

**Phase 4 — Validation and golden regeneration:**
19. Add unit test `test_gdna_harmonic_length.py` (6.9).
20. Run `pytest tests/` — all tests should pass except goldens.
21. Human-review diffs of representative goldens
    (`loci_df.tsv`, `quant.tsv`) on 2–3 fixtures.
22. `pytest tests/ --update-golden` once diffs look sane.
23. Run `mctp_vcap_rna20m_dna80m` benchmark. Confirm mega-locus
    gdna_rate ≥ 0.75 and size-binned flatness.
24. Run all 13 synthetic benchmark conditions via
    `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml`.
25. Commit in phases 1–4 as separate commits for bisectability;
    open single PR.

### 6.9 Validation

**Unit test plan (new file `tests/test_gdna_harmonic_length.py`):**

```python
def test_nh1_unique_mapper():
    # Fragment with NH=1, single transcript of t_span=50kb, gfp=300
    # Expected: gdna_ll = gdna_fl + gdna_log_sp + LOG_HALF + 0
    #                   - log(50000 + flank - 300 + 1)
    ...

def test_nh2_similar_intervals():
    # Fragment with NH=2, both hits on t_span=50kb transcripts
    # Expected: harmonic mean = 50kb; gdna_ll = ... - log(50k+flank-gfp+1)
    # (logsumexp of two equal terms minus log(2))
    ...

def test_nh2_dissimilar_intervals():
    # Fragment with NH=2, t_span values 5kb and 500kb
    # Expected: harmonic mean of e_h values ≈ 2e(5kb)/(e(5kb)+e(500kb))
    # Confirm value matches hand-computed logsumexp exactly.
    ...

def test_spliced_kills_gdna():
    # Fragment with NH=2, one hit spliced, one unspliced
    # Expected: gdna_log_liks[u] = -inf
    ...

def test_all_spliced_nh_gdna_zero():
    # Fragment with NH=2, both spliced
    # Expected: gdna_log_liks[u] = -inf, mm_nh_gdna = 0
    ...
```

**Mega-locus regression:** extend `tests/test_golden_output.py`
(or add `tests/test_mega_locus_regression.py`) to load the
`dna80m`-like synthetic fixture and assert
`loci_df.gdna_rate >= 0.75` for the largest connected component.

**Benchmark targets** on `mctp_vcap_rna20m_dna80m`:

1. Mega-locus `gdna_rate`: 0.500 → ≥ 0.75.
2. Size-binned `gdna_rate` variance (max − min over size bins):
   today ≈ 0.25; target ≤ 0.05.
3. Overall `mrna` count vs truth: no regression (±1 %).
4. Overall `nrna` count: expected reduction in siphon.

### 6.10 Risk assessment

| Risk | Likelihood | Mitigation |
|---|---|---|
| `t_span_` is not already loaded in scorer | none — verified at `scoring.cpp:151` | — |
| `gdna_flank` is unavailable at scorer-construction time | low — `calibration.gdna_fl_model.mean` is computed before scorer build; one grep confirms | Verify in step 6 of rollout; fall back to plumbing if needed |
| NH=1 non-MM path uses `best_t` of fragment which may be −1 | real — fragments with zero candidates skipped by `new_cands > 0` guard; gdna emission is inside that guard | Covered by existing guard |
| Logsumexp numerical overflow | low — inputs are bounded log-probs with typical range ±30 | `lse_update` uses stable online algorithm |
| Spliced+unspliced mixed MM group with current "any-spliced → -inf" rule is too harsh | possible but unchanged vs status quo | Out of scope for this PR; file issue if needed |
| Downstream consumer of `genomic_footprints` was overlooked | low — grep came back clean outside the plumbing path | Phase 3 gated on Phase 2 passing; easy to revert |
| VBEM convergence changes due to gDNA log-lik dynamic range | low — shrinks effective length ~4 orders of mag in mega-locus (log gain ~+9 nats) which is within the solver's RNA candidate range | Monitor `em_n_iter` in benchmark summary.json; flag if >2× baseline |

### 6.11 Follow-up: locoregional $\lambda_G(x)$ (the "Simulated Option C")

Unchanged from previous draft. With Option B owning the
**likelihood** side of gDNA physics (local denominator), a future
PR can refine the **prior** side by replacing the global
$\lambda_G$ in `locus.py::compute_locus_priors` with a windowed
estimator over intergenic reads. The two multiply cleanly. No
code changes in scoring, EM, partition, or buffer would be
required for that follow-up.

---

## 6. Implementation plan for Option B (ARCHIVED — representative-hit draft)

The representative-hit approximation below was superseded by the
in-scorer harmonic mean above after review identified that tying
the gDNA effective length to the RNA-winning transcript
reintroduces the RNA↔gDNA coupling Option B is designed to
break. Kept for design-trail traceability.

### 6.0 Guiding principle: simplify, don't bolt-on

This change is a **net code deletion**. The current locus-level
gDNA bias correction is a special case wedged into RNA's per-candidate
bias-correction pathway, complete with a sentinel `gdna_idx` row in
`bias_profiles`, a per-locus `gdna_span` parameter threaded through
pipeline.py into the C++ EM, and a fragile `locus.gdna_span +
gdna_flank` computation at two call sites. All of that is conceptually
*one thing* — "the denominator of the uniform-position prior for
gDNA" — but it is currently spread across four files and two
languages.

Option B's natural home for this logic is **a single Python function
in `locus.py`** that bakes the gDNA length correction into
`ScoredFragments.gdna_log_liks` once, immediately after `build_loci()`
returns. Every line of the existing gDNA bias-correction plumbing
becomes redundant and is deleted:

| Deleted | Location |
|---|---|
| `gdna_flank = int(calibration.gdna_fl_model.mean)` threading | `pipeline.py` (2 sites) |
| `gdna_spans` parameter in batch EM signature | `pipeline.py`, `em_solver.cpp` |
| `sub.bias_profiles[sub.gdna_idx] = gdna_span` | `em_solver.cpp::extract_locus_sub_problem_from_partition` |
| Implicit gDNA row inside `apply_bias_correction_uniform` | `em_solver.cpp` (guard added) |

What remains is a **clean, symmetric model**: per-RNA-candidate bias
correction inside the EM (unchanged); gDNA log-liks arrive
pre-corrected from Python (new, single source of truth). The
`Locus.gdna_span` field stays but becomes purely informational —
used only by `get_loci_df` for `loci.feather`, never by the
likelihood.

### 6.1 Mathematical specification

For each EM unit $u$ with best transcript $t^*(u)$ and genomic
footprint $\ell_u$:

- $I(u)$ := the merged interval in unit $u$'s locus that contains
  $t^*(u)$. (By construction of `build_loci`, every transcript is
  contained in exactly one merged interval.)
- $L(u) := \text{len}(I(u)) + \text{gdna\_flank}$ where
  $\text{gdna\_flank} = \lfloor \mathbb{E}[\ell_{\text{gDNA}}] \rfloor$.
- Per-unit gDNA length correction:
  $$\text{e}(u) \;=\; \max\!\big(L(u) - \ell_u + 1,\; 1\big)$$
  $$\texttt{gdna\_log\_liks}[u] \mathrel{-}= \log \text{e}(u)$$

For unique-mapper units this is **exact**: NH=1, and $t^*$ is the
sole alignment. For multimapper units, this uses $t^*(u)$ as a
representative hit from the fragment's NH alignments. Since all NH
hits land in the *same MultiLocus* by union-find construction, and
since gDNA log-lik contribution is dominated by $\ell_u$ (identical
across hits — it's the fragment's own length) and nm (near-identical
across hits for reasonable aligners), the representative-hit
approximation differs from the full harmonic-mean formula only when
the NH hits land in merged intervals of dramatically different
sizes. This can be upgraded to full per-hit later if benchmarks
justify it (see Section 6.7); the API is designed so this upgrade
is a one-function change.

### 6.2 The new `locus.py` function

```python
def apply_gdna_length_correction(
    sf: ScoredFragments,
    loci: list[Locus],
    gdna_flank: int,
    num_transcripts: int,
) -> None:
    """Bake per-unit gDNA length correction into sf.gdna_log_liks.

    For each unit u, subtract log(max(L(u) - gfp(u) + 1, 1)) from
    gdna_log_liks[u], where L(u) = len(merged_interval(t*(u))) +
    gdna_flank.

    This replaces the former per-locus ``bias_profiles[gdna_idx] =
    gdna_span`` mechanism.  It is the *only* place in the codebase
    where the gDNA effective-length model is defined.

    Called once, in pipeline.quant_from_buffer, after build_loci and
    before partition_and_free.  Operates in-place on sf.
    """
    t_to_local_len = np.zeros(num_transcripts, dtype=np.int64)
    for locus in loci:
        # Each merged interval's effective length
        # (cached once per locus during build_loci; see 6.3)
        for iv_idx, (ref, start, end) in enumerate(locus.merged_intervals):
            L = (end - start) + gdna_flank
            t_to_local_len[locus.interval_t_indices[iv_idx]] = L

    best_t = sf.locus_t_indices                    # int32[n_units]
    gfp    = sf.genomic_footprints                 # int32[n_units]
    L_loc  = t_to_local_len[best_t]                # int64[n_units]
    e      = np.maximum(L_loc - gfp.astype(np.int64) + 1, 1)
    finite = np.isfinite(sf.gdna_log_liks)
    sf.gdna_log_liks[finite] -= np.log(e[finite])
```

That is the full algorithm. ~25 lines including docstring.

### 6.3 Augmenting `Locus` with `interval_t_indices`

We need a mapping from each merged interval to its member transcripts,
built once during `build_loci`. Extend the `Locus` dataclass:

```python
@dataclass(slots=True)
class Locus:
    locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray
    gdna_span: int
    merged_intervals: list  # list[tuple[ref, start, end]]
    interval_t_indices: list  # NEW: list[np.ndarray], one per merged_interval
```

Populated inside `build_loci`'s existing merge loop by sorting
member transcripts by `(ref_code, start)` and partitioning by
interval membership. Incremental O(n_transcripts_in_locus) cost
per locus, no new passes.

### 6.4 Pipeline wiring

One new call in `pipeline.quant_from_buffer`, one deletion of the
flank-threading in `_run_locus_em_partitioned`:

```python
# pipeline.py — after build_loci
loci = build_loci(em_data, index)
apply_gdna_length_correction(
    em_data, loci,
    gdna_flank=int(calibration.gdna_fl_model.mean),
    num_transcripts=index.num_transcripts,
)
alpha_gdna, alpha_rna = _compute_priors(...)
partitions = partition_and_free(em_data, loci)
# NOTE: `gdna_flank` no longer threaded into _run_locus_em_partitioned
_run_locus_em_partitioned(
    estimator, partitions, loci, index,
    alpha_gdna, alpha_rna, em_config,
    # removed: gdna_flank=gdna_flank,
    emit_locus_stats=emit_locus_stats,
    annotations=annotations,
)
```

Two delete-ops:
- `pipeline.py` line ~815: `gdna_flank = int(...)` moves into
  `apply_gdna_length_correction` call; original binding and
  downstream threading go away.
- `pipeline.py` line ~545: `_build_locus_meta` no longer has
  access to (nor need for) `gdna_flank`. `locus.gdna_span` still
  used for `locus_span_bp` reporting.

### 6.5 C++ EM solver simplification

The solver currently stores a per-component `bias_profiles[nc]`
array where `bias_profiles[gdna_idx] = gdna_span`. Under Option B
there is no per-locus gDNA profile. Two micro-changes to
`em_solver.cpp`:

**(a)** `extract_locus_sub_problem_from_partition` signature drops
its `gdna_span` parameter:

```cpp
// BEFORE
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,
    double alpha_gdna,
    int64_t gdna_span,              // REMOVE
    ...);

// AFTER
static void extract_locus_sub_problem_from_partition(
    LocusSubProblem& sub,
    const PartitionView& pv,
    double alpha_gdna,
    ...);
```

Body: delete `sub.bias_profiles[sub.gdna_idx] = gdna_span;`.
`sub.bias_profiles` is still sized `nc` (RNA + gDNA) to keep the
indexing convention, but the gDNA slot is unused. Alternatively,
size it to `n_t` — mechanically a little more work but
philosophically cleaner. Prefer sizing to `n_t` and guarding the
bias-correction loop with `if (local_comp < n_t)`.

**(b)** `apply_bias_correction_uniform` gains an `n_t` parameter
and guards the gDNA row:

```cpp
static void apply_bias_correction_uniform(
    double*        log_liks,
    const int32_t* t_indices,
    const int32_t* tx_starts,
    const int32_t* tx_ends,
    const int64_t* profile_lengths,   // now indexed [0..n_t)
    int32_t        n_t,                // NEW
    size_t         n_candidates)
{
    for (size_t i = 0; i < n_candidates; ++i) {
        int32_t local_comp = t_indices[i];
        if (local_comp >= n_t) continue;   // gDNA row, pre-corrected
        ...
        log_liks[i] -= std::log((double)eff_len);
    }
}
```

**(c)** `_batch_locus_em_partitioned` Python-side signature drops
the `gdna_spans` argument; C++-side binding drops the parameter;
all call sites updated.

### 6.6 Touch-point summary

| File | Net change |
|---|---|
| `src/rigel/scored_fragments.py` | +1 field on `Locus` (`interval_t_indices`) |
| `src/rigel/locus.py` | +1 function `apply_gdna_length_correction`. Minor addition in `build_loci` to populate `interval_t_indices`. |
| `src/rigel/pipeline.py` | +1 call site for `apply_gdna_length_correction`. −4 lines threading `gdna_flank` / `gdna_spans` into the batch EM. |
| `src/rigel/estimator.py` | −1 arg (`gdna_spans`) from `_batch_locus_em_partitioned` wrapper. |
| `src/rigel/native/em_solver.cpp` | −`gdna_span` parameter from extraction. −1 line setting `bias_profiles[gdna_idx]`. +1 guard in `apply_bias_correction_uniform`. −1 kwarg in nanobind binding. |
| `tests/` | +1 test for `apply_gdna_length_correction`. +1 mega-locus regression test. Regenerate goldens once validated. |

**Rough LOC budget:** +80 lines of new code (function + test + dataclass field). −50 lines of deleted plumbing. **Net: +30 lines**, with the conceptual surface reduced (one Python function owns all gDNA length physics).

### 6.7 Validation and upgrade path

**Validation targets** on `mctp_vcap_rna20m_dna80m`:

1. Mega-locus `gdna_rate`: 0.500 → ~0.78.
2. Size-binned `gdna_rate` (weighted mean): flatten from
   `0.75 / 0.66 / 0.63 / 0.53 / 0.50` to ~0.78 across all bins.
3. Overall `nrna` count: should drop meaningfully as gDNA absorbs
   its fair share.
4. Overall `mrna` count vs truth: should not regress.

**Unit tests** for `apply_gdna_length_correction`:

- NH=1 unit in a 50 kb interval → correction = −log(50000 + flank
  − gfp + 1). Closed form.
- Two-interval locus with a unit whose best transcript sits in
  the smaller interval → uses the small-interval length.
- Spliced unit (gdna_log_liks = −inf) → unchanged.
- Unique mapper vs multimapper: same code path, identical output
  for matched transcripts.

**Upgrade path to full per-hit harmonic mean (Section 1.4's exact
formula):**

If benchmarks show residual error on mega-loci where multimapper
hits span merged intervals of dramatically different sizes,
replace Section 6.2's representative-hit lookup with a per-hit
harmonic mean. The upgrade is **localized to
`apply_gdna_length_correction`** — no scoring changes, no C++
changes, no partition changes — and requires threading per-hit
genomic positions through the scorer. That in turn requires a new
CSR-ish pair of arrays (`gdna_hit_offsets`, `gdna_hit_positions`)
emitted alongside `gdna_log_liks`. Left for a follow-up PR only
if needed.

### 6.8 Future: locoregional $\lambda_G(x)$ (the "Simulated C")

With Option B in place, the per-fragment gDNA likelihood uses the
local merged-interval length — a spatial denominator. The per-locus
$\gamma_\ell$ prior is the spatial numerator. Together they define
a local Poisson density expectation for gDNA.

**Current**: $\gamma_\ell$ is computed from a global $\lambda_G$
applied uniformly across mappable bp, weighted by observed region
totals. Fine for uniform contamination; blind to CNVs.

**Upgrade**: replace the global $\lambda_G$ with a windowed
$\lambda_G(x)$ estimated from intergenic reads in sliding genomic
windows (2–10 Mb). Loci falling inside amplified regions get a
higher $\gamma_\ell$ prior; loci in deleted regions get a lower
one. Calibration-only change — no EM, scoring, or locus plumbing
affected.

**Important:** this is an enhancement orthogonal to Option B.
Option B is a correctness fix for the mathematical formula of gDNA
likelihood; locoregional $\lambda_G(x)$ is a refinement of the
*prior* feeding into that formula. They multiply cleanly.

### 6.9 Rollout sequence

1. Add `Locus.interval_t_indices` field and populate in `build_loci`.
   Existing consumers are unaffected (additive field). Tests still
   pass.
2. Add `apply_gdna_length_correction` function with unit tests.
3. Remove `gdna_flank` / `gdna_spans` threading from pipeline.py
   and estimator.py. Remove C++ parameter.
4. Adjust `apply_bias_correction_uniform` with `n_t` guard.
5. Recompile (`pip install --no-build-isolation -e .`) and run
   full test suite. Regenerate golden outputs (`pytest --update-golden`)
   after human review of before/after `loci_df.tsv` deltas.
6. Run `mctp_vcap_rna20m_dna80m` benchmark. Compare mega-locus
   `gdna_rate` to expected 0.78–0.80.
7. Run all 13 synthetic benchmark conditions; analyze.
8. Single PR, no feature flag. The old behavior is mathematically
   incorrect for multi-interval loci; Option B is strictly better.

---
