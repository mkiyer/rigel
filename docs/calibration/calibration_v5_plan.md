Calibration v5 — Consolidated Plan & Implementation
Status: AUTHORITATIVE. Single source of truth for the v5 calibration system. Consolidates and supersedes:

calibration_v5_plan_v4.md (settled design + Phase 0 contract)
calibration_v5_plan_v5.md (Phase 1 layout-sweep refactor)
calibration_v5_plan_v6.md (gDNA-only boundary flux narrowing)
v5_phase4_handoff.md (Phase 4 native implementation handoff)
Date: 2026-05-04 Predecessors archived to docs/archive/calibration/.

This document is the ONLY authoritative reference. Phases 0–5 are implemented; Phases 6–11 remain. Re-opening any locked design decision requires a written amendment appended to this file.

1. Verified codebase facts
All facts in this section have been verified against the live tree and inform every phase.

em_solver.cpp (Phase 7a target) fans alpha_rna over all eligible non-gDNA components (mRNA + nRNA together). A native fix is required.
bam_scanner.cpp::AssembledFragment.genomic_footprint() is the bbox of aligned exons only. The gap between paired mates is structurally absent. ⇒ Fragment masks 110/101/111 are observable; they are NOT errors.
pipeline.py calls calibration before quant_from_buffer; result is plumbed through to quant.
quant_from_buffer raises if calibration.gdna_fl_model is None. ⇒ Phase 0 supplies a stub with a populated gdna_fl_model.
cli.py _PARAM_SPECS registry: new config knobs go through it.
locus.py merged_intervals is a union of disconnected spans.
pyproject.toml — pandas + pyarrow only. No Polars.
index.py::load_reference_lengths keys by string ref_name. ref_lengths.feather is the canonical ref table.
No estimate_kappa exists yet; Phase 6 implements it.
cli.py summary.json has a single "calibration" key. Reuse with version: "v5" discriminator.
resolve_context.h ResolvedFragment carries splice_type, chimera_type, chimera_gap, num_hits, nm, genomic_footprint, genomic_start, plus per-T arrays. This is the right struct for the native accumulator hook.
read_length in resolve_context.h is the sum of aligned exon block spans (M/D/=/X), NOT the full query sequence length. Soft-clipped bases are excluded.
2. Settled design (LOCKED)
These decisions are closed. Re-opening any requires a written amendment.

2.1 Region partition
Cluster all transcripts on either strand into genic spans (maximal runs of overlapping transcripts; strand-agnostic).
Complement of genic spans on each reference = intergenic spans.
Within each genic span, merge all-strand exons → EXON regions; the interstices become INTRON regions.
EXON wins at every position with any-strand exon coverage.
Synthetic nRNAs are excluded from geometry.
Tiny regions are retained; filtered downstream when too small.
The partition tiles the entire genome with no gaps and no overlaps.
Naming note: the codebase uses genic / intergenic, not island / gap.

2.2 Region attributes (on-disk schema)
regions.feather schema:

Column	Dtype	Notes
region_id	int64	Monotonic across genome; row index ≡ region_id
ref_name	string	Matches ref_lengths.feather; canonical on disk
start	int64	0-based inclusive
end	int64	0-based exclusive
type	uint8	INTERGENIC=0, INTRON=1, EXON=2
strand	uint8	NONE=0, POS=1, NEG=2, AMBIG=3 (strict bitwise OR)
tx_pos_bp	int64	bp of (+)-strand tx overlap
tx_neg_bp	int64	bp of (−)-strand tx overlap
exon_pos_bp	int64	bp of (+)-strand exon overlap
exon_neg_bp	int64	bp of (−)-strand exon overlap
boundary_flux_left	bool	left boundary is internal (touches an INTRON of the same genic span)
boundary_flux_right	bool	right boundary is internal (touches an INTRON of the same genic span)
Invariants validated on load:

region_df.index.equals(region_df["region_id"]).
(region_df["end"] - region_df["start"] > 0).all().
For each ref: regions are non-overlapping, sorted, and cover [0, ref_lengths[ref]) with no gaps.
Every ref_name exists in ref_lengths.feather.
region_cgr (cgranges) is built at load time, keyed by integer ref_id. Native code receives ref_id (translated from string at the single load point — see §2.10).

2.3 Fragment categorization (8-state mask)
3-bit mask {has_intergenic, has_intron, has_exon}:

Mask	Name	gDNA FL pool?	Density numerator?
100	INTERGENIC_ONLY	yes	intergenic
010	INTRON_ONLY	yes	intron
001	EXON_ONLY	no	—
011	EXON_INTRON	yes	exon-intron
110	INTRON_INTERGENIC	no (annotation-gap)	annotation-gap (QC)
101	EXON_INTERGENIC	no (annotation-gap)	annotation-gap (QC)
111	EXON_INTRON_INTERGENIC	no (annotation-gap)	annotation-gap (QC)
000	NO_OVERLAP	no	impossible if partition complete; warning counter
Bit order in C++ / numpy arrays (canonical): bit 0 = EXON, bit 1 = INTRON, bit 2 = INTERGENIC. The Python type encoding {INTERGENIC=0, INTRON=1, EXON=2} is translated to the bitmask via type_mask = 1 << (2 - type) at the C++/Python boundary (pipeline.py ~line 254).

Annotation-gap exclusion (locked): include_annotation_gap_in_fl_pool defaults to False. Annotation-gap fragments are often unannotated splicing or retained introns (i.e., RNA, not gDNA); including them biases the gDNA FL distribution toward RNA shape. Per-mask QC counts are surfaced in the result; if n_pool_annotation_gap >> 0, that is a data-quality signal worth investigating.

2.4 Counting unit (locked)
A fragment contributes at most +1 to a given region for category counts, region counts, and boundary-flux counters. Implementation uses SmallRegionSet (§4.2) for amortized O(1) dedup with no per-fragment heap traffic.

2.5 Boundary-flux counters — gDNA-only
Per EXON region 
R
R we accumulate two int64_t counters during the BAM scan, gated on splice_type == SPLICE_UNSPLICED:

u_L(R) — count of unspliced fragments straddling the left boundary.
u_R(R) — count of unspliced fragments straddling the right boundary.
Tests (with boundary_tol, default 0):

u_L : frag_start  < region_start - tol   AND  frag_end > region_start
u_R : frag_start  < region_end           AND  frag_end > region_end + tol
frag_start = min(exons[*].start), frag_end = max(exons[*].end) across the fragment's aligned blocks. A fragment may contribute to both u_L and u_R of the same region only when it spans the entire region.

These counters are inputs, not estimators. They feed the global EXON-INTRON gDNA density (§2.6), which is then propagated to per-Locus EXON gDNA mass via §2.7. There is no per-EXON-region gDNA density: per-region EXON unspliced counts are structurally ambiguous (mRNA, nRNA, and gDNA all contribute), so we never compute or publish such a quantity.

At scan time, counters are recorded for every EXON boundary regardless of boundary_flux_left/right flag — eligibility filtering happens at aggregation time (§2.6.2). This Phase-4 contract lets us swap eligibility policies without re-scanning the BAM.

Eligibility rules (locked)
Internal (non-terminal) EXON boundary that abuts an INTRON of the same genic span → flag = True.
Terminal EXON boundaries (TSS, TES) → flag = False. Rationale: hybrid-capture probes typically do not tile across TSS/TES; unspliced flux at terminal boundaries is dominated by capture-edge artefacts rather than genuine gDNA. Excluding terminal boundaries makes the estimator capture-aware by construction.
Single-exon transcripts → both flags = False.
The Phase-1 builder (src/rigel/calibration/regions.py) implements these rules; no schema change is required.

Why no spliced (RNA) boundary flux
Under the locked v5.v1 policy prior_weight_rna[nrna] = 0 (§2.9), every non-mRNA fragment in an exon is attributed to gDNA in the prior. The unspliced boundary flux is the gDNA estimator. A spliced/unspliced ratio would estimate the non-mRNA fraction — collapsing to the same answer the §2.7 density-propagation route already gives. Two counters, one estimator path, one diagnostic surface.

2.6 Three global gDNA densities (locked)
The library is summarized by exactly three scalar gDNA densities in units of fragments / bp:

Density	Numerator (Σ over all regions)	Denominator (Σ over all regions)
ρ
^
INTERGENIC
gDNA
ρ
^
​
  
INTERGENIC
gDNA
​
 	per_region_counts[R, 0b100] over INTERGENIC R	
L
eff
(
R
)
L 
eff
​
 (R) over INTERGENIC R
ρ
^
INTRON
gDNA
ρ
^
​
  
INTRON
gDNA
​
 	per_region_counts[R, 0b010] over INTRON R	
L
eff
(
R
)
L 
eff
​
 (R) over INTRON R
ρ
^
EXON-INTRON
gDNA
ρ
^
​
  
EXON-INTRON
gDNA
​
 	
u
L
(
R
)
 
1
L
(
R
)
+
u
R
(
R
)
 
1
R
(
R
)
u 
L
​
 (R)1 
L
​
 (R)+u 
R
​
 (R)1 
R
​
 (R) over EXON R	
(
1
L
(
R
)
+
1
R
(
R
)
)
⋅
L
ˉ
gDNA
(1 
L
​
 (R)+1 
R
​
 (R))⋅ 
L
ˉ
  
gDNA
​
  over EXON R
where the gDNA effective length of a region of length 
L
(
R
)
L(R) is

L
eff
(
R
)
  
=
  
max
⁡
 ⁣
(
0
,
  
L
(
R
)
−
L
ˉ
gDNA
+
1
)
L 
eff
​
 (R)=max(0,L(R)− 
L
ˉ
  
gDNA
​
 +1)

i.e., the number of distinct fragment start positions of mean gDNA fragment length that fit inside 
R
R. 
L
ˉ
gDNA
L
ˉ
  
gDNA
​
  is the mean of the gDNA FL model from Phase 8 (during Phase 6 we use the bootstrap-stub mean; see §5.1.2).

There is no fourth (EXON) density. Per-EXON-region unspliced fragment counts (mask 0b001) cannot be attributed to gDNA without ambiguity vs nRNA and pre-mRNA isoforms, so we never compute or publish a per-EXON gDNA density. EXON gDNA mass is derived from 
ρ
^
EXON-INTRON
gDNA
ρ
^
​
  
EXON-INTRON
gDNA
​
  via §2.7 — boundary flux measures gDNA density on exonic territory unambiguously, and that is the only EXON-side density we need.

The legacy pi_pool scalar is a fragment-count-weighted summary (summary.json only); it does not feed prior assembly.

2.6.1 Two-level EB shrinkage
A two-level Negative-Binomial empirical-Bayes hierarchy:

global density (per type, computed above)
        ↓ shrink
locoregional density (per Locus interval, per type)
There is no per-reference (per-chromosome) level. The locoregional level is computed on demand at locus-prior assembly time (§2.7); it is not a pre-computed table.

Standard NB shrinkage:

ρ
^
loco
(
L
,
t
)
  
=
  
N
(
L
,
t
)
+
κ
t
⋅
ρ
^
global
(
t
)
L
eff
(
L
,
t
)
+
κ
t
ρ
^
​
  
loco
​
 (L,t)= 
L 
eff
​
 (L,t)+κ 
t
​
 
N(L,t)+κ 
t
​
 ⋅ 
ρ
^
​
  
global
​
 (t)
​
 

where 
N
(
L
,
t
)
N(L,t) and 
L
eff
(
L
,
t
)
L 
eff
​
 (L,t) are the type-
t
t fragment count and effective-length sum aggregated over all type-
t
t regions inside Locus 
L
L (clipped to 
L
L's contiguous interval). 
κ
t
κ 
t
​
  is the global NB overdispersion estimate per type (Phase 6.1).

2.6.2 Counters → density conversion
EXON-INTRON numerator/denominator are computed after filtering EXON regions by boundary_flux_left/right flags (terminal exclusion, §2.5):

for R in exon_regions:
    n_L = u_L[R] if R.boundary_flux_left  else 0
    n_R = u_R[R] if R.boundary_flux_right else 0
    sides = int(R.boundary_flux_left) + int(R.boundary_flux_right)
    if sides == 0:
        continue
    numerator   += n_L + n_R
    denominator += sides * L_gdna_mean
The sides * L̄_gDNA denominator follows the boundary-flux geometry: each eligible boundary contributes one mean-fragment-length window in which an unspliced gDNA fragment can be observed straddling it.

2.7 Per-Locus gDNA mass (locked)
Recall (§2.16): a MultiLocus is a connected component of transcripts; a Locus is one of its contiguous genomic intervals. A MultiLocus contains one or more Locus intervals. gDNA mass is estimated per Locus, independently, then handed to the EM either as a single combined MultiLocus prior (v5.v1) or as per-Locus priors (v5.v2; see §2.7.3).

2.7.1 Algorithm (per Locus)
For each Locus interval 
L
=
(
ref
,
start
,
end
)
L=(ref,start,end):

Collect overlapping calibration regions via region_cgr:

region_ids = region_cgr.overlap(ref_id, L.start, L.end)
rows = region_df.loc[region_ids]   # invariant: index == region_id
Per type 
t
∈
{
INTERGENIC
,
INTRON
,
EXON-INTRON
}
t∈{INTERGENIC,INTRON,EXON-INTRON}, compute the locoregional EB-shrunk density. Numerator and denominator aggregate inside 
L
L:

INTERGENIC / INTRON: as in the global formula (§2.6) but restricted to regions of that type that overlap 
L
L (clip to 
L
L for the effective-length sum).
EXON-INTRON: aggregate 
u
L
,
u
R
u 
L
​
 ,u 
R
​
  over EXON regions in 
L
L that are flag-eligible (§2.6.2).
Then EB-shrink each to the global density of that type (§2.6.1).

Convert each locoregional density to a predicted gDNA fragment count by multiplying by the corresponding effective length sum inside 
L
L:

N
^
t
gDNA
(
L
)
  
=
  
ρ
^
loco
(
L
,
t
)
  
⋅
  
L
eff
(
L
,
t
)
N
  
t
gDNA
​
 (L)= 
ρ
^
​
  
loco
​
 (L,t)⋅L 
eff
​
 (L,t)

For 
t
=
EXON-INTRON
t=EXON-INTRON the effective length is summed over all EXON regions in 
L
L — not just flag-eligible ones. The density already encodes capture-aware bias by construction; the predicted count is the gDNA mass we expect on all exonic territory inside 
L
L:

L
eff
EXON-INTRON
(
L
)
  
=
  
∑
R
∈
EXON
(
L
)
max
⁡
 ⁣
(
0
,
  
L
(
R
)
−
L
ˉ
gDNA
+
1
)
L 
eff
EXON-INTRON
​
 (L)= 
R∈EXON(L)
∑
​
 max(0,L(R)− 
L
ˉ
  
gDNA
​
 +1)

Sum the three predicted counts to get the Locus gDNA mass:

N
^
gDNA
(
L
)
  
=
  
N
^
INTERGENIC
gDNA
(
L
)
+
N
^
INTRON
gDNA
(
L
)
+
N
^
EXON-INTRON
gDNA
(
L
)
N
  
gDNA
 (L)= 
N
  
INTERGENIC
gDNA
​
 (L)+ 
N
  
INTRON
gDNA
​
 (L)+ 
N
  
EXON-INTRON
gDNA
​
 (L)

Compute the Locus gDNA fraction:

π
^
gDNA
(
L
)
  
=
  
N
^
gDNA
(
L
)
 
/
 
N
obs
(
L
)
π
^
  
gDNA
 (L)= 
N
  
gDNA
 (L)/N 
obs
​
 (L)

where 
N
obs
(
L
)
N 
obs
​
 (L) is the exact fragment count assigned to 
L
L from the buffer (an integer). 
π
^
π
^
  is clipped to 
[
0
,
1
]
[0,1].

This procedure is per Locus; each contiguous interval has its own EB-shrunk density triple, its own gDNA mass estimate, and its own gDNA fraction. It is never run once over the union of a MultiLocus's intervals — that would mix calibration regions from disjoint genomic neighborhoods and inflate the apparent locoregional sample size.

2.7.2 MultiLocus → EM prior assembly (v5.v1)
For Phase 7b (v5.v1), the EM still receives a single gDNA prior mass per MultiLocus. We aggregate per-Locus estimates by simple sum:

def assemble_multilocus_prior(multi_locus, per_locus_estimates):
    n_gdna = sum(e.n_gdna for e in per_locus_estimates)
    n_obs  = sum(e.n_obs  for e in per_locus_estimates)
    n_rna  = max(0.0, n_obs - n_gdna)
    return MultiLocusPrior(
        n_gdna=n_gdna, n_rna=n_rna,
        per_locus=per_locus_estimates,   # retained for diagnostics + v5.v2
    )
The MultiLocus EM then splits n_rna across mRNA and nRNA components via prior_weight_rna (§2.9).

2.7.3 Per-Locus gDNA components in EM (v5.v2 spike, deferred)
The per-Locus estimates are retained on MultiLocusPrior.per_locus. A v5.v2 EM extension can emit one gDNA component per Locus interval rather than a single shared gDNA component, assigning fragment 
f
f to the gDNA component whose Locus interval contains 
f
f's genomic footprint. Design spike scheduled in Phase 7d.

2.7.4 Sanity checks (raised at assembly time)
A Locus with no overlapping calibration regions → raise RuntimeError (BAM reference mismatch with index).
π
^
gDNA
(
L
)
>
1
π
^
  
gDNA
 (L)>1 before clipping → log a warning with 
L
L's coordinates, the three component counts, and 
N
obs
(
L
)
N 
obs
​
 (L).
A MultiLocus with 
π
^
<
0
π
^
 <0 or NaN → raise (defensive guard against EB underflow).
2.8 c_base(ℓ) formula
OPEN — Phase 7c design spike. Candidates: fixed / coverage-scaled / EB-shrinkage-derived / hybrid. Decision recorded in this document after Phase 7c.

2.9 nRNA suppression — continuous prior weight in C++
Native EM solver receives per-locus prior_weight_rna: std::vector<float>, one entry per component, in [0, 1]. The fan-out becomes:

double total_weighted_rna_coverage = 0.0;
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] > 0.0 && i != gdna_idx) {
        total_weighted_rna_coverage += prior_weight_rna[i] * coverage_totals[i];
    }
}
// later, in the prior-assignment loop:
prior_out[i] = baseline + std::max(
    alpha_rna * prior_weight_rna[i] * coverage_totals[i]
        / total_weighted_rna_coverage,
    EM_LOG_EPSILON);
v5.v1: prior_weight_rna[mrna] = 1.0, prior_weight_rna[nrna] = 0.0. v5.v2 may substitute any continuous weight with no ABI change.

2.10 Reference identity contract
On disk: all DataFrames key references by ref_name (string).
In native arrays / cgranges payloads: integer ref_id derived from the canonical ref_lengths.feather order.
Translation point: exactly one place — TranscriptIndex.load() builds ref_name_to_id: dict[str, int] and ref_names: list[str] from ref_lengths.feather. All callers consume one or the other; no ad hoc translation elsewhere.
Calibration overlap engine: the native RegionIndex (§2.15) is the single source for region-overlap queries on the calibration path. Its set_regions rejects unknown ref_id rather than auto-creating refs, which would silently shadow a ref-table mismatch.
2.11 Library / style decisions
pandas + pyarrow only. No Polars.
C++17, -O3, LTO. Hot-path code uses inline storage; no per-fragment heap allocation in CalibrationAccumulator::observe.
All C++ public types live under namespace rigel::calibration.
2.12 Diagnostic surfaces explicitly removed
The following per-region diagnostics are not part of the v5 result schema:

s_L, s_R per-EXON spliced-block counters.
ρ
^
RNA
(
R
)
ρ
^
​
  
RNA
​
 (R) per-region RNA boundary density.
π
^
gdna
(
R
)
π
^
  
gdna
​
 (R) per-region boundary-flux ratio.
Per-region (per-EXON) gDNA density of any kind. EXON gDNA mass is a derived quantity from the global EXON-INTRON density propagated through the per-Locus EB step (§2.7); it never lives at region granularity.
Per-reference (per-chromosome) EB level. The hierarchy is global → locoregional only.
density_by_type carries exactly three entries: INTERGENIC, INTRON, EXON-INTRON — each a single global density.

2.13 Calibration observation policy (locked)
Exactly one policy decides which fragments enter the calibration accumulator. The policy is enforced at the scanner observation site (not downstream); the payload is therefore "usable-only plus exclusion counters" — there is no implicit mixture and no post-hoc filtering.

Fragment class	v5.v1 policy	Counter
Unique mapper, non-chimeric, non-artifact, resolved	observed	n_calibration_observed
Unique mapper, intergenic-only (resolved-as-intergenic)	observed	(same)
Multi-mapper (num_hits > 1), any class	excluded	n_calibration_excluded_multimap
Chimeric (any chimera_type != CHIMERA_NONE)	excluded	n_calibration_excluded_chimera
Artifact (splice_type == SPLICE_ARTIFACT)	excluded	n_calibration_excluded_artifact
Out-of-region / no-overlap (mask == 0b000)	observed (counted in global_counts[0]); flagged in QC	n_calibration_oor
Unresolved (no transcript hits and no intergenic resolution)	not observed (never reaches the call site)	—
Counters are integers; the payload remains integer-typed.

2.13.1 Multimapper rationale
Multimapper inclusion would require either (a) per-hit +1 increments (over-counting by NH, biasing intergenic and repeat-region densities upward) or (b) fractional 1/NH increments (requires float counters, mass-conservation tests, and complicates κ-estimation since the underlying observations are no longer Poisson/NB). v5.v1 chooses exclusion.

This causes systematic under-estimation of gDNA density on multi-mappable regions (notably intergenic repeat tracts and pseudogenes). The compensation is deferred to v5.v2 via mappability-adjusted effective lengths: 
L
eff
map
(
R
)
=
∫
R
m
(
x
)
 
d
x
L 
eff
map
​
 (R)=∫ 
R
​
 m(x)dx where 
m
(
x
)
∈
[
0
,
1
]
m(x)∈[0,1] is per-bp mappability under the read-length configuration. Then density estimates rebalance correctly without re-introducing multimapper observations into the numerator. Tracked in §9 open questions.

2.13.2 Chimera & artifact rationale
Chimeric fragments have ambiguous genomic footprints (paired mates on different references / strands / large gaps). Letting them inflate per_region_counts would corrupt every downstream density. Artifact splice junctions are pre-classified as gDNA-derived junk by the splice blacklist; including them in the calibration sample double-counts that signal in the boundary-flux numerator while simultaneously biasing the gDNA FL pool toward the artifact length distribution.

2.13.3 Test contract
tests/test_v5_observation_policy.py (Phase 5.5, ≥ 6 cases) — pins each row of the table above using a synthetic mini-BAM with one fragment per class. Failure of this test is the canonical indicator that the observation contract has drifted.

2.14 Per-Locus n_obs derivation (locked)
Per-Locus observed fragment counts (N_obs(L) in §2.7) are derived at calibration time from FragmentBuffer geometry:

def per_locus_n_obs(buffer, multi_locus, ref_name_to_id):
    """For each Locus interval in `multi_locus.loci`, count fragments
    whose (ref_id, genomic_start, genomic_start + genomic_footprint)
    falls inside the interval."""
    counts: dict[Locus, int] = {l: 0 for l in multi_locus.loci}
    for f_idx in buffer.fragments_for_multi_locus(multi_locus.id):
        ref_id = buffer.ref_id[f_idx]
        gs     = buffer.genomic_start[f_idx]
        ge     = gs + buffer.genomic_footprint[f_idx]
        for l in multi_locus.loci:
            if l.ref_id == ref_id and gs >= l.start and ge <= l.end:
                counts[l] += 1
                break    # disjoint Locus intervals → unique assignment
    return counts
buffer.fragments_for_multi_locus already exists implicitly via the CSR fragment_router mapping in scan.py; Phase 7b adds a public view (MultiLocusFragmentIndex) so calibration can iterate without duplicating router internals.

The denominator 
N
obs
(
L
)
N 
obs
​
 (L) in 
π
^
gDNA
(
L
)
π
^
  
gDNA
 (L) is the post-policy observed count (§2.13 policy). The same exclusion policy applies to both numerator (predicted gDNA from regions) and denominator (Locus-assigned fragments) — i.e., the buffer view used for per_locus_n_obs is filtered through the same multimap / chimera / artifact mask. Otherwise 
π
^
π
^
  would be biased by having different fragment populations in the numerator and denominator.

2.15 Region-lookup contract (locked)
Calibration does region overlap lookups in two contexts:

Context	Lookup engine	Owner
Hot path (BAM scan)	rigel::calibration::RegionIndex (binary-search, native)	bam_scanner.cpp
Phase 7b per-Locus prior	same RegionIndex, via a Python nanobind wrapper	_calibration_impl.RegionIndex
The Python region_cgr member of TranscriptIndex is legacy for calibration use as of Phase 6.5. It remains available for resolver debugging; calibration code paths must not depend on it. The native RegionIndex is exposed to Python via rigel.native.RegionIndex (Phase 6.5 deliverable) with a single overlap method:

class RegionIndex:
    def overlap(self, ref_id: int, start: int, end: int) -> np.ndarray:
        """Return int32 array of region_ids overlapping [start, end) on ref_id."""
set_regions raises on unknown ref_id (no auto-creation; the contract from §2.10 is enforced at the seam).

2.16 Nomenclature — Locus vs MultiLocus (locked)
The current Locus class in src/rigel/locus.py represents a connected component of transcripts that share at least one multi-mapping fragment, possibly spanning multiple disconnected genomic intervals (chimeric multimappers, tandem-duplicate paralogs, etc.). This makes "Locus" a misleading name — it can refer to a set of disjoint genomic regions.

v5 codifies the distinction:

Term	Definition
Locus	A single contiguous genomic interval (ref, start, end). The atomic unit of gDNA-density estimation (§2.7).
MultiLocus	A connected component of transcripts (the current Locus); contains one or more Locus intervals. The atomic unit of EM execution.
MultiLocus.loci	The ordered list of contiguous Locus intervals composing the MultiLocus (≥ 1). For most genes this list has length 1.
A Phase 6.5 refactor renames the existing Locus class to MultiLocus and introduces a new lightweight Locus dataclass for contiguous intervals (see §5.1.5). All call sites that consume Locus.merged_intervals are migrated to iterate MultiLocus.loci explicitly. This refactor is mechanical (no semantic change to EM behaviour) and lands as its own phase to keep the diff reviewable.

Why the rename matters for calibration
§2.7 estimates per-Locus gDNA mass independently. If we ran the recipe over MultiLocus.merged_intervals directly, we would mix calibration regions from disjoint genomic neighborhoods (e.g., a multi-mapping pseudogene cluster on chr3 and chrX), inflating the locoregional sample size with statistically independent signal and biasing the EB shrinkage. Iterating MultiLocus.loci and estimating per Locus is the correct unit.

3. Phase status tracker
Phase	Title	Status
0	Cleanup + transitional contract	completed (2026-04-28)
1	Region partition (layout sweep + two emitters)	completed (2026-05-04)
2	Region index persistence + load validation	completed (2026-05-04)
3	C++ scanner refactor (bit-identical sink seam)	completed (2026-05-04)
4	Online fragment categorization (8-state mask + per-region counts + per-mask FL hist)	completed (2026-05-04)
4.5	RNA boundary-flux scrub (v6 lock-in)	completed (2026-05-04)
5	Unspliced boundary-flux counters (gDNA-only)	completed (2026-05-04)
5.5	Observation policy + payload promotion	next — not-started
6	Global gDNA densities + estimate_kappa (3 types)	not-started
6.5	Native RegionIndex Python wrapper + Locus→MultiLocus rename	not-started
7	Per-Locus gDNA mass + nRNA suppression + c_base spike	not-started
7d	Per-Locus gDNA components in EM (v5.v2 spike)	not-started (deferred)
8	FL pool refresh + CalibrationResultV5 schema	not-started
9	Pipeline integration	not-started
10	Validation harness + benchmarks	not-started
11	Documentation + release	not-started
4. Implemented phases (0–5)
This section is the operational reference for code that already exists. Phases 6–11 are spec-only (§5).

4.1 Phase 0 — Cleanup + transitional contract
Deleted v1 calibration code (_simple.py, _categorize.py, _result.py, tests/test_calibration_simple.py, tests/test_categorize.py, tests/test_gdna.py, tests/test_gdna_harmonic_length.py).

Moved general FL utilities to src/rigel/frag_length_mixture.py and src/rigel/frag_length_eb.py.

Public surface of rigel.calibration package:

from rigel.calibration import (
    CalibrationStub,                # transitional dataclass
    bootstrap_fl_calibration,       # builds CalibrationStub from buffer
    calibrate_v5,                   # raises NotImplementedError until Phase 9
)
pipeline.run_pipeline calls bootstrap_fl_calibration instead of the legacy calibrate_gdna. The stub's gdna_fl_model is the global FL shrunk toward the RNA FL with a heavy pseudocount, logged as WARNING: using bootstrap FL fallback (calibration v5 not yet active). summary.json reports calibration.active == False.

compute_locus_priors_from_partitions was removed; the temporary replacement compute_default_locus_priors returns uniform Dirichlet pseudocounts and is the only locus-prior path between Phase 0 and Phase 7.

4.2 Phase 1 — Region partition (layout sweep + two emitters)
src/rigel/index.py was refactored so a single per-reference sweep over sorted transcripts feeds both intervals.feather (resolver cgranges payload) and regions.feather (calibration partition):

@dataclass(frozen=True)
class _IntergenicSpan:
    start: int
    end: int

@dataclass(frozen=True)
class _GenicSpan:
    start: int
    end: int
    transcripts: list[Transcript]

def _iter_reference_layout(
    ref_length: int,
    transcripts: list[Transcript],
) -> Iterator[_IntergenicSpan | _GenicSpan]:
    """Walk one reference in sorted-transcript order, yielding the
    interleaved sequence of intergenic and genic spans that tiles
    [0, ref_length) exactly."""
Two consumers (_emit_genomic_intervals, _emit_regions) operate on the same stream. build_index_artifacts(transcripts, ref_lengths) -> (intervals_df, regions_df) is the single driver.

intervals.feather is byte-identical to pre-refactor; regions.feather matches the §2.2 schema. region_id is assigned globally in genomic order so regions_df.index == regions_df["region_id"].

4.3 Phase 2 — Region index persistence + load validation
regions.feather is written next to transcripts.feather. TranscriptIndex.load() reads it, validates the §2.2 invariants, builds region_cgr keyed by integer ref_id (from ref_name_to_id), and exposes regions (DataFrame) and region_cgr properties. Stale indexes fail-fast with explicit "rebuild" guidance. INDEX_FORMAT_VERSION was bumped.

4.4 Phase 3 — C++ scanner refactor (bit-identical sink seam)
A CalibrationSink abstract header was introduced in src/rigel/native/ with a NoOpSink implementation. The Phase 3 seam is now retired in favor of the Phase 4 in-place accumulator (see §4.5); the header remains as inert scaffolding.

4.5 Phase 4 — Online fragment categorization
4.5.1 Native data structures
src/rigel/native/calibration/region_index.h — rigel::calibration::RegionIndex:

class RegionIndex {
public:
    RegionIndex(const int32_t* ref_ids,
                const int64_t* starts,
                const int64_t* ends,
                const uint8_t* type_masks,
                int32_t n_regions,
                int32_t n_refs);

    void overlap_into(int32_t ref_id, int64_t start, int64_t end,
                      SmallRegionSet& out) const noexcept;

    uint8_t type_mask(int32_t region_id) const noexcept;
    int64_t start    (int32_t region_id) const noexcept;  // Phase 5
    int64_t end      (int32_t region_id) const noexcept;  // Phase 5
    int32_t n_regions() const noexcept;
    int32_t n_refs()    const noexcept;
};
Sorted-array storage; overlap_into uses binary search + forward-walk on the per-ref slice. No interval tree, no cgranges in the hot path — the partition tiles each reference, so binary search suffices.

src/rigel/native/calibration/small_region_set.h — rigel::calibration::SmallRegionSet:

class SmallRegionSet {
public:
    static constexpr int kInline = 16;
    void clear() noexcept;
    bool insert(int32_t v) noexcept;     // linear-scan dedup
    int32_t size() const noexcept;
    template <class F> void for_each(F&& f) const noexcept;
private:
    int32_t inline_[kInline];
    int32_t size_ = 0;
    std::vector<int32_t> spilled_;       // empty unless size_ > kInline
};
Owned per-worker; reused across all fragments. Spilled vector is reused, never freed mid-scan.

src/rigel/native/calibration/accumulator.h — rigel::calibration::CalibrationAccumulator:

inline constexpr int kMaskStates = 8;
inline constexpr int kFLBins     = 1024;
inline constexpr uint8_t kExonBit = 0x01;

class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int32_t n_regions);

    void observe(int32_t splice_type,
                 int32_t genomic_footprint,
                 const ExonBlock* exons,
                 size_t n_exons,
                 const RegionIndex& regions,
                 int32_t boundary_tol = 0) noexcept;

    void merge_from(const CalibrationAccumulator& other) noexcept;

    int32_t n_regions() const noexcept;
    const std::array<int64_t, kMaskStates>& global_counts() const noexcept;
    const std::vector<int64_t>& per_region_counts() const noexcept;
    const std::array<std::array<int64_t, kFLBins>, kMaskStates>&
        fl_hist() const noexcept;
    const std::vector<int64_t>& u_L() const noexcept;   // Phase 5
    const std::vector<int64_t>& u_R() const noexcept;   // Phase 5

private:
    int32_t n_regions_;
    std::array<int64_t, kMaskStates> global_counts_{};
    std::vector<int64_t> per_region_counts_;            // [rid*8 + mask]
    std::array<std::array<int64_t, kFLBins>, kMaskStates> fl_hist_{};
    std::vector<int64_t> u_L_;                          // Phase 5, n_regions
    std::vector<int64_t> u_R_;                          // Phase 5, n_regions
    SmallRegionSet scratch_;
};
4.5.2 observe() algorithm — four passes
void CalibrationAccumulator::observe(...) noexcept {
    if (n_exons == 0) return;
    scratch_.clear();

    // Pre-compute fragment footprint (used by passes A and D).
    int32_t frag_ref_id = exons[0].ref_id;
    int64_t frag_start  = exons[0].start;
    int64_t frag_end    = exons[0].end;
    bool    same_ref    = true;
    for (size_t i = 1; i < n_exons; ++i) {
        if (exons[i].ref_id != frag_ref_id) same_ref = false;
        frag_start = std::min<int64_t>(frag_start, exons[i].start);
        frag_end   = std::max<int64_t>(frag_end,   exons[i].end);
    }

    // (A) Splice-class geometry dispatch.
    if (splice_type == SPLICE_UNSPLICED && same_ref) {
        regions.overlap_into(frag_ref_id, frag_start, frag_end, scratch_);
    } else {
        // SPLICE_SPLICED_*, SPLICE_IMPLICIT, or mixed-ref unspliced.
        for (size_t i = 0; i < n_exons; ++i)
            regions.overlap_into(exons[i].ref_id, exons[i].start,
                                 exons[i].end, scratch_);
    }

    // (B) mask = OR of region.type_mask over touched regions.
    uint8_t mask = 0;
    scratch_.for_each([&](int32_t rid) noexcept {
        mask |= regions.type_mask(rid);
    });

    // (C) global, per-region, FL counters with FINAL mask.
    global_counts_[mask] += 1;
    scratch_.for_each([&, mask](int32_t rid) noexcept {
        per_region_counts_[size_t(rid) * kMaskStates + mask] += 1;
    });
    if (genomic_footprint > 0 && genomic_footprint < kFLBins)
        fl_hist_[mask][genomic_footprint] += 1;

    // (D) Phase 5: unspliced boundary-flux counters.
    if (splice_type == SPLICE_UNSPLICED && same_ref) {
        scratch_.for_each([&](int32_t rid) noexcept {
            if ((regions.type_mask(rid) & kExonBit) == 0) return;
            const int64_t rs = regions.start(rid);
            const int64_t re = regions.end(rid);
            if (frag_start < rs - boundary_tol && frag_end > rs)
                u_L_[size_t(rid)] += 1;
            if (frag_start < re && frag_end > re + boundary_tol)
                u_R_[size_t(rid)] += 1;
        });
    }
}
merge_from() sums global_counts_, fl_hist_, per_region_counts_, u_L_, u_R_ element-wise. Worker merge is associative and unit-tested.

4.5.3 Scanner integration
WorkerState (in bam_scanner.cpp) carries rigel::calibration::CalibrationAccumulator cal_acc, sized to n_regions. BamScanner::scan():

Builds RegionIndex once at the top from region_df columns (passed through nanobind) plus the resolver's ref_name → ref_id table.
Constructs each WorkerState with n_regions.
In process_qname_group_threaded, calls ws.cal_acc.observe(...) at the three append sites (intergenic-only, chimeric, resolved non-chimeric).
After workers join, folds each cal_acc into cal_acc_merged_.
Exports the payload (§4.5.4).
4.5.4 Returned scan payload
result["calibration"] (nanobind dict, returned by BamScanner::scan()):

{
    "global_counts":      np.ndarray, shape (8,),                int64
    "per_region_counts":  np.ndarray, shape (n_regions, 8),      int64
    "fl_hist":            np.ndarray, shape (8, kFLBins),        int64
    "u_L":                np.ndarray, shape (n_regions,),        int64    # Phase 5
    "u_R":                np.ndarray, shape (n_regions,),        int64    # Phase 5
    "n_regions":          int,
    "fl_bins":            int,                                    # == kFLBins
}
Arrays are zero-copy numpy views over the merged accumulator's storage; nanobind owns lifetime via nb::capsule. The 3-bit mask order on the trailing axis: bit 0 = EXON, bit 1 = INTRON, bit 2 = INTERGENIC.

4.5.5 Python wrapper
src/rigel/calibration/scan_payload.py:

@dataclass(frozen=True)
class CalibrationScanPayload:
    global_counts: np.ndarray       # (8,)            int64
    per_region_counts: np.ndarray   # (n_regions, 8)  int64
    fl_hist: np.ndarray             # (8, fl_bins)    int64
    u_L: np.ndarray                 # (n_regions,)    int64
    u_R: np.ndarray                 # (n_regions,)    int64

    @classmethod
    def from_scan_dict(cls, d: dict) -> "CalibrationScanPayload": ...
from_scan_dict validates shapes against the declared n_regions and fl_bins and takes defensive copies (the native module's capsule lifetime is tied to the scan call).

4.6 Phase 4.5 — RNA boundary-flux scrub (v6 lock-in)
Audit confirmed that no spliced-flux scaffolding (s_L, s_R, spliced_left/right, boundary_flux_rna) leaked into the codebase. tests/test_v6_boundary_flux_eligibility.py pins the eligibility contract:

Two-exon (+) → exon[0] (False, True); exon[1] (True, False).
Three-exon (+) → exon[0] (False, True); exon[1] (True, True); exon[2] (True, False).
Single-exon either strand → (False, False).
Opposite-strand transcripts sharing TSS/TES → terminal flags False.
Regression sweep over a 4-transcript GTF.
Docstring of RegionRecord.boundary_flux_left/right updated with the v6 (gDNA-only, hybrid-capture-aware) interpretation.

4.7 Phase 5 — Unspliced boundary-flux counters
Implemented as Pass D inside CalibrationAccumulator::observe() (§4.5.2). Two int64_t vectors of length n_regions — u_L_ and u_R_ — gated on splice_type == SPLICE_UNSPLICED AND single-ref. Counters recorded for all EXON regions including terminals; eligibility filtering via boundary_flux_left/right is a Phase 6 concern.

tests/test_v6_phase5_boundary_flux.py (6 tests):

Payload shape/dtype ((n_regions,) int64, non-negative).
mRNA-only sim → u_L.sum() == u_R.sum() == 0 (spliced or within-exon fragments cannot be SPLICE_UNSPLICED and straddle an EXON boundary).
gDNA contamination → u_L.sum() > 0, u_R.sum() > 0, and non-zero counts only on EXON regions.
Per-fragment dedup bound: max(u_L) ≤ total_fragments.
1-thread vs 4-thread merge equality.
Terminal exons not silently filtered (Phase 5 records counts regardless of eligibility flags).
5. Remaining phases (5.5–11)
5.0 Phase 5.5 — Observation policy + payload promotion
Scope. Lock the §2.13 observation policy at the scanner; widen the scan payload with exclusion counters; promote CalibrationScanPayload onto PipelineResult so Phase 6+ are pure transforms. Required prerequisite for Phase 6 — Phase 6 densities are statistically meaningless on the current payload because chimeras and multimappers inflate the numerator.

5.0.1 Native scanner patch
src/rigel/native/bam_scanner.cpp — three observation sites currently violate §2.13:

Site	Line region	Current behavior	v5.v1 fix
Intergenic resolved	~1401	observed (already gated is_unique_mapper)	unchanged; verify via test
Chimeric resolved	~1424	observed unconditionally	gate on !result.get_is_chimeric() (skip + bump n_excluded_chimera)
Genic non-chimeric resolved	~1503	observed once per hit (multimapper inflation)	gate on is_unique_mapper (skip + bump n_excluded_multimap once per fragment); also gate on splice_type != SPLICE_ARTIFACT (bump n_excluded_artifact)
Add to WorkerState:

struct CalibrationExclusionCounters {
    int64_t multimap = 0;
    int64_t chimera  = 0;
    int64_t artifact = 0;
    int64_t oor      = 0;   // mask == 0 after observe()
    int64_t observed = 0;   // sanity: == sum(global_counts) excluding mask 0
};
Counters merge element-wise across workers (merge_from). Exported in result["calibration"] alongside the existing arrays.

5.0.2 Scanner ABI extension
result["calibration"] gains:

"n_excluded_multimap": int64
"n_excluded_chimera":  int64
"n_excluded_artifact": int64
"n_oor":               int64
"n_observed":          int64
CalibrationScanPayload (scan_payload.py) gains the same five fields. from_scan_dict validates n_observed == global_counts.sum() − global_counts[0] (the OR-zero NO_OVERLAP row).

5.0.3 Payload promotion onto PipelineResult
src/rigel/pipeline.py:

@dataclass(frozen=True)
class PipelineResult:
    ...
    calibration_payload: CalibrationScanPayload    # NEW
run_pipeline stops discarding the payload returned by scan_and_buffer; calibration consumers in Phase 6+ accept the payload as a pure input.

5.0.4 Tests
tests/test_v5_observation_policy.py (≥ 6 cases — pin the §2.13.3 contract):

Unique mapper, non-chimeric → observed; counters all zero.
Multimapper (NH=2) → excluded; n_excluded_multimap == 1, global_counts.sum() == 0.
Chimeric (cis-strand-diff) → excluded; n_excluded_chimera == 1.
Artifact (SPLICE_ARTIFACT) → excluded; n_excluded_artifact == 1.
Implicit-spliced unique mapper → observed; mask = 0b001.
Unknown-ref fragment → never reaches the call site (resolver drops it); counters zero. Sanity assertion only.
Plus a balance test:

n_observed + n_excluded_multimap + n_excluded_chimera \
           + n_excluded_artifact + n_unresolved == n_total_fragments
5.0.5 Exit gate
tests/test_v5_observation_policy.py green.
All existing tests still green (the policy fix changes calibration counts but no other code path).
PipelineResult.calibration_payload populated end-to-end on tests/scenarios_aligned/.
§2.13 balance equation holds on a 1k-fragment sample.
5.1 Phase 6 — Global gDNA densities + estimate_kappa
Scope. Compute the three global gDNA densities (§2.6) and one NB-overdispersion estimate 
κ
t
κ 
t
​
  per type. No per-Locus work in this phase — that is Phase 7. No per-reference EB level.

Artifacts:

src/rigel/calibration/density_global.py   # NEW: compute_global_densities()
src/rigel/calibration/_kappa.py           # NEW: estimate_kappa()
tests/test_v5_density_global.py           # NEW
tests/test_v5_estimate_kappa.py           # NEW
5.1.1 compute_global_densities — public surface
@dataclass(frozen=True)
class GlobalGdnaDensity:
    type: Literal["INTERGENIC", "INTRON", "EXON-INTRON"]
    rho: float                  # fragments / bp
    n_fragments: int64
    eff_length_bp: int64
    n_regions_used: int64       # for EXON-INTRON: # eligible boundaries
    kappa: KappaEstimate        # per-type NB overdispersion (Phase 6.2)


@dataclass(frozen=True)
class GlobalDensityTable:
    intergenic:    GlobalGdnaDensity
    intron:        GlobalGdnaDensity
    exon_intron:   GlobalGdnaDensity
    gdna_fl_mean:  float        # source: Phase 8 result (or stub during Phase 6)

    def by_type(self) -> dict[str, GlobalGdnaDensity]: ...


def compute_global_densities(
    region_df: pd.DataFrame,
    payload: CalibrationScanPayload,
    gdna_fl_mean: float,
) -> GlobalDensityTable: ...
Implementation is a pure pass over region_df + payload:

INTERGENIC: filter type == 0; numerator Σ per_region_counts[R, 0b100]; denominator Σ max(0, L(R) − L̄_gDNA + 1).
INTRON: filter type == 1; numerator Σ per_region_counts[R, 0b010]; denominator as above.
EXON-INTRON: filter type == 2; for each EXON region apply the flag-eligibility filter (§2.6.2): n_L = u_L if flag_left else 0, similarly n_R; numerator Σ (n_L + n_R); denominator Σ (sides * L̄_gDNA) where sides ∈ {0, 1, 2}.
If a denominator is zero (no eligible regions / boundaries), 
ρ
^
ρ
^
​
  is set to 0.0 and KappaEstimate.fallback_used is forced True.

5.1.2 gdna_fl_mean source
During Phase 6 (before Phase 8 lands the real gDNA FL model), use bootstrap_calibration.gdna_fl_model.mean() from Phase 0. Plumbing this through is a one-line change. Phase 9 swaps it for the calibration-derived mean.

5.1.3 estimate_kappa — per-region NB MoM with heterogeneous exposures
Naïve MoM that collapses all regions into one mean exposure 
L
eff
‾
L 
eff
​
 
​
  conflates true overdispersion with exposure heterogeneity, biasing 
κ
^
κ
^
  low whenever region lengths span more than ~½ of an order of magnitude. v5.v1 uses a per-region Method-of-Moments NB estimator that respects exposures.

Model: 
N
R
∣
ρ
R
∼
Poisson
(
ρ
R
⋅
L
eff
(
R
)
)
N 
R
​
 ∣ρ 
R
​
 ∼Poisson(ρ 
R
​
 ⋅L 
eff
​
 (R)), 
ρ
R
∼
Γ
(
κ
,
κ
/
ρ
^
)
ρ 
R
​
 ∼Γ(κ,κ/ 
ρ
^
​
 ) (mean 
ρ
^
ρ
^
​
 , shape 
κ
κ). Marginally 
N
R
∼
NegBin
(
μ
R
=
ρ
^
⋅
L
eff
(
R
)
,
shape
=
κ
)
N 
R
​
 ∼NegBin(μ 
R
​
 = 
ρ
^
​
 ⋅L 
eff
​
 (R),shape=κ) with 
V
a
r
(
N
R
)
=
μ
R
+
μ
R
2
/
κ
Var(N 
R
​
 )=μ 
R
​
 +μ 
R
2
​
 /κ.

Sum-of-squares moment match (over 
n
n regions):

∑
R
(
N
R
−
μ
R
)
2
  
=MoM
  
∑
R
μ
R
  
+
  
1
κ
∑
R
μ
R
2
R
∑
​
 (N 
R
​
 −μ 
R
​
 ) 
2
  
=
MoM
  
R
∑
​
 μ 
R
​
 + 
κ
1
​
  
R
∑
​
 μ 
R
2
​
 

Solving for 
κ
κ with 
μ
R
=
ρ
^
⋅
L
eff
(
R
)
μ 
R
​
 = 
ρ
^
​
 ⋅L 
eff
​
 (R):

κ
^
  
=
  
∑
R
μ
R
2
 
∑
R
(
N
R
−
μ
R
)
2
  
−
  
∑
R
μ
R
 
κ
^
 = 
∑ 
R
​
 (N 
R
​
 −μ 
R
​
 ) 
2
 −∑ 
R
​
 μ 
R
​
 
∑ 
R
​
 μ 
R
2
​
 
​
 

KAPPA_DEFAULT       = 100.0   # neutral fallback
KAPPA_MIN           = 1.0
KAPPA_MAX           = 1.0e6
MIN_REGIONS_FOR_MOM = 5

@dataclass(frozen=True)
class KappaEstimate:
    value: float
    fallback_used: bool
    n_regions: int
    excess_variance: float    # Σ(N − μ)² − Σ μ; negative ⇒ underdispersion

def estimate_kappa(
    counts: np.ndarray,        # (n_regions,) int64
    eff_lengths: np.ndarray,   # (n_regions,) float64
    rho_hat: float,            # global density of this type
) -> KappaEstimate:
    """Per-region NB Method-of-Moments. Falls back to KAPPA_DEFAULT
    when:
      - n_regions < MIN_REGIONS_FOR_MOM
      - rho_hat == 0 or all eff_lengths == 0 (μ identically zero)
      - excess_variance ≤ 0 (Poisson-or-tighter than Poisson)
      - resulting κ outside [KAPPA_MIN, KAPPA_MAX] → clipped, not flagged
    """
    if counts.size < MIN_REGIONS_FOR_MOM or rho_hat == 0.0:
        return KappaEstimate(KAPPA_DEFAULT, True, int(counts.size), 0.0)
    mu = rho_hat * eff_lengths
    if mu.sum() == 0.0:
        return KappaEstimate(KAPPA_DEFAULT, True, int(counts.size), 0.0)
    excess = float(((counts - mu) ** 2).sum() - mu.sum())
    if excess <= 0.0:
        return KappaEstimate(KAPPA_DEFAULT, True, int(counts.size), excess)
    kappa = float((mu * mu).sum() / excess)
    kappa = min(max(kappa, KAPPA_MIN), KAPPA_MAX)
    return KappaEstimate(kappa, False, int(counts.size), excess)
This uses rho_hat as the global mean (not the per-region MoM mean), matching the §2.6.1 EB shrinkage target. The estimator is consistent when 
L
eff
(
R
)
L 
eff
​
 (R) varies arbitrarily — no implicit "single mean exposure" assumption.

5.1.4 Tests
tests/test_v5_estimate_kappa.py (degenerate cases):

n_regions = 0 → fallback.
n_regions = 3 < threshold → fallback.
All counts equal, all eff_lengths equal → underdispersion fallback.
Underdispersion (Poisson-like) → fallback.
rho_hat == 0 → fallback.
Healthy NegBin generative test (uniform exposures) → recovers κ within 25%.
Heterogeneous-exposure consistency. Generate n_regions = 200 regions with L_eff log-uniform over 
[
100
,
100
 
000
]
[100,100000] (3 orders of magnitude). Sample 
ρ
R
∼
Γ
(
κ
∗
,
κ
∗
/
ρ
)
ρ 
R
​
 ∼Γ(κ 
∗
 ,κ 
∗
 /ρ) with 
κ
∗
∈
{
2
,
20
,
200
}
κ 
∗
 ∈{2,20,200}, then 
N
R
∼
Poisson
(
ρ
R
⋅
L
eff
,
R
)
N 
R
​
 ∼Poisson(ρ 
R
​
 ⋅L 
eff,R
​
 ). Recover 
κ
^
κ
^
  within 
±
30
%
±30% for all three target values — this test would fail under the old single-mean-exposure formula for 
κ
∗
=
2
κ 
∗
 =2.
tests/test_v5_density_global.py:

Hand-counted 3-region scenario: pin the three densities to their exact values.
EXON-INTRON with all flags False → density = 0, fallback κ.
INTRON-only sample (no INTERGENIC fragments) → INTERGENIC density = 0, others non-zero.
Pure-mRNA sample (no gDNA) → all three densities < 1e-3.
50/50 dna20m sample → INTERGENIC ≈ INTRON ≈ EXON-INTRON within 30%.
5.1.5 Exit gate
11+ tests green.
Protected suite green.
Integration smoke: compute_global_densities runs end-to-end on tests/scenarios_aligned/ and produces a GlobalDensityTable with three finite, non-negative densities.
5.1.5 Phase 6.5 — Native RegionIndex Python wrapper + Locus → MultiLocus rename
Scope. Two coupled refactors that together deliver the Phase 7b prerequisites: (a) a Python-callable wrapper around the existing native RegionIndex so calibration overlap queries go through one engine, and (b) the mechanical Locus → MultiLocus rename. No semantic change to EM behavior.

5.1.5.1 Native RegionIndex Python wrapper
src/rigel/native/_calibration_impl.cpp (existing module): expose rigel::calibration::RegionIndex as a nanobind class.

# rigel.native (re-export)
class RegionIndex:
    def __init__(
        self,
        ref_ids:    np.ndarray,   # int32, (n_regions,)
        starts:     np.ndarray,   # int64, (n_regions,)
        ends:       np.ndarray,   # int64, (n_regions,)
        type_masks: np.ndarray,   # uint8, (n_regions,)
        n_refs:     int,
    ) -> None: ...

    def overlap(self, ref_id: int, start: int, end: int) -> np.ndarray: ...
    def n_regions(self) -> int: ...
TranscriptIndex.load() builds and caches a RegionIndex instance alongside the existing region_cgr. A docstring on region_cgr marks it legacy for calibration use (the resolver continues to use it until a separate cleanup).

5.1.5.2 Locus → MultiLocus rename
Scope. Mechanical rename + minimal API extraction. No semantic change to EM behavior. Required before Phase 7 because §2.7 operates per Locus, not per MultiLocus.

Artifacts:

src/rigel/locus.py
    - rename class Locus → class MultiLocus
    - add @dataclass(frozen=True) class Locus:   # NEW
          ref: str; ref_id: int; start: int; end: int
    - add MultiLocus.loci: tuple[Locus, ...]
      (built from existing merged_intervals at construction)
    - keep MultiLocus.merged_intervals as a property → Locus list
      adapter for one release; deprecate in Phase 11.

src/rigel/estimator.py, src/rigel/pipeline.py, src/rigel/scan.py,
tests/test_locus_partition.py, tests/test_pipeline_routing.py, ...
    - update type hints `Locus` → `MultiLocus`
    - call sites that iterate merged_intervals migrate to .loci

tests/test_v5_locus_rename.py            # NEW
Implementation steps:

vscode_renameSymbol Locus → MultiLocus across the workspace.
Reintroduce Locus as a small frozen dataclass (4 fields).
Build MultiLocus.loci once at construction by walking the existing merged_intervals (already a list of contiguous intervals).
Update tests.
Confirm EM golden outputs unchanged.
Tests (tests/test_v5_locus_rename.py):

MultiLocus.loci is non-empty for every connected component in a 1-ref synthetic GTF.
For a multi-ref scenario (pseudogene paralogs), MultiLocus.loci has length > 1 and len(set(l.ref for l in m.loci)) > 1.
Σ (l.end - l.start) == sum of len(merged_intervals) parity check.
Exit gate: EM golden outputs bit-identical; protected suite green.

5.2 Phase 7 — Per-Locus gDNA mass + nRNA suppression + prior strength
Phase 7 is the consumer phase: it turns the global gDNA densities from Phase 6, the calibration regions from Phase 1–2, and the canonical Locus geometry from Phase 6.5 into a per-MultiLocus Dirichlet prior that the EM consumes. It also lands the long-planned continuous nRNA-suppression knob in the native EM.

Subphases (3 active + 1 deferred):

7a — Native EM prior_weight_rna (continuous nRNA suppression ABI). Self-contained C++/Python change; unblocks the v5.v1 nRNA policy (weight=0 for nRNA components).
7b — Per-Locus gDNA mass estimator (§2.7) and MultiLocus prior aggregation (§2.7.2). The bulk of the new Python surface.
7c — c_base(ℓ) design spike + decision. Determines how Phase-7b intensities (π̂(L)) become EM Dirichlet pseudocounts. Locks §2.8.
7d — (deferred to v5.v2) Per-Locus gDNA components in EM (§2.7.3).
Sequencing & dependencies. 7a is independent of 7b/c and may land first (or in parallel), since its only Python contract is "EM accepts an extra prior_weight_rna array; pass None ≡ all-ones for backwards compatibility." 7b consumes Phase 6 (GlobalDensityTable), Phase 6.5 (Locus, MultiLocus.loci, native RegionIndex), and the existing FragmentBuffer. 7c consumes 7b output and is purely a Python calibration / scoring exercise — no native changes. 7c writes the chosen c_base formula back into §2.8 of this document.

End-to-end pipeline integration is deferred to Phase 9. Phase 7 ships the artifacts (estimator, schema, c_base formula); Phase 9 swaps compute_default_locus_priors for the new MultiLocusPrior → (α_gdna, α_rna) adapter.

Artifact inventory:

src/rigel/native/em_solver.cpp                     # MODIFIED: prior_weight_rna fan-out  (7a)
src/rigel/native/_em_impl.cpp                      # MODIFIED: nanobind signature        (7a)
src/rigel/locus.py                                 # MODIFIED: prior_weight_rna helper   (7a)
src/rigel/calibration/density_loco.py              # NEW: shrink_to_loco() (~30 LOC)     (7b)
src/rigel/calibration/_arrays.py                   # NEW: pre-flattened region+payload   (7b)
src/rigel/calibration/locus_prior.py               # NEW: estimator + aggregation        (7b)
src/rigel/calibration/_locus_n_obs.py              # NEW: §2.14 n_obs derivation         (7b)
src/rigel/calibration/c_base.py                    # NEW: c_base candidates + dispatch   (7c)
src/rigel/calibration/__init__.py                  # MODIFIED: export 7b/7c surface
scripts/calibration/c_base_spike.py                # NEW: 5-scenario spike harness       (7c)
tests/test_v5_em_prior_weight.py                   # NEW                                  (7a)
tests/test_v5_density_loco.py                      # NEW                                  (7b)
tests/test_v5_per_locus_gdna_mass.py               # NEW                                  (7b)
tests/test_v5_multilocus_aggregation.py            # NEW                                  (7b)
tests/test_v5_locus_n_obs.py                       # NEW                                  (7b)
tests/test_v5_c_base.py                            # NEW                                  (7c)
docs/calibration/calibration_v5_plan.md            # UPDATED: §2.8 amended                (7c)
5.2.0 Architecture & data flow
The Phase 7 entry point is a single pure function:

def assemble_priors_v5(
    multi_loci:        list[MultiLocus],
    em_data:           ScoredFragments,           # for unit→fragment
    buffer:            FragmentBuffer,            # for genomic_start / footprint
    index:             TranscriptIndex,           # region_df, region_index, ref_lengths
    payload:           CalibrationScanPayload,    # Phase 5/5.5 output
    global_densities:  GlobalDensityTable,        # Phase 6 output
    gdna_fl_mean:      float,                     # Phase 8 output (Phase 0 stub during 7)
    c_base:            CBaseFormula,              # Phase 7c choice
) -> PriorTable:
    ...
PriorTable carries everything Phase 9 needs:

@dataclass(frozen=True)
class PriorTable:
    multi_locus_priors: tuple[MultiLocusPrior, ...]   # one per MultiLocus
    alpha_gdna:          np.ndarray  # (n_multi_loci,) float64 — EM input
    alpha_rna:           np.ndarray  # (n_multi_loci,) float64 — EM input
    prior_weight_rna:    list[np.ndarray]   # one per MultiLocus, float32
Three cleanly separable layers (top → bottom):

┌────────────────────────────────────────────────────────────────┐
│  LAYER A — orchestration                                        │
│    locus_prior.assemble_priors_v5() : list[MultiLocus] → PriorTable│
└──────────────────────────────┬─────────────────────────────────┘
                               │
┌──────────────────────────────▼─────────────────────────────────┐
│  LAYER B — per-Locus estimator (pure functions)                 │
│    _arrays.RegionArrays      : pre-flatten region_df once        │
│    _arrays.PayloadArrays     : pre-flatten payload columns once  │
│    _locus_n_obs.partition()  : MultiLocus.unit_indices → per-Locus│
│    locus_prior.estimate_locus_gdna()  : Locus → LocusGdnaEstimate│
│    density_loco.shrink_to_loco()      : EB shrinkage scalar      │
└──────────────────────────────┬─────────────────────────────────┘
                               │
┌──────────────────────────────▼─────────────────────────────────┐
│  LAYER C — c_base dispatch (pure)                                │
│    c_base.evaluate(formula, ctx)  : LocusGdnaEstimate × ctx → α │
└────────────────────────────────────────────────────────────────┘
Layer B is the hot path (called once per Locus interval, ~50 K times per pipeline run). It must be allocation-free in the inner loop — hence the _arrays.py factor-out (one pandas → numpy hop at the top of assemble_priors_v5, then pure numpy fancy-indexing per Locus).

Key design decisions affecting all three subphases:

pi_gdna(L) is intensive; mass scaling is c_base's job. §2.7.1 defines π̂(L) = N̂_gDNA(L) / N_obs(L) as a fraction. Phase 7b's estimator outputs (π̂, n_gdna_predicted, n_obs) per Locus. Phase 7c's c_base(L) then converts the intensity into a Dirichlet pseudocount strength — α_gdna = c_base(L) · π̂(L) etc. This decouples the estimator from the prior-strength policy and makes the c_base spike a pure post-hoc computation that does not require re-running the estimator.

N_obs(L) is derived from EM units, not the raw BAM. §2.14 originally specified iterating the buffer with the calibration exclusion policy. We reinterpret this as: N_obs(L) is the count of EM units assigned to Locus L. This is what the Dirichlet prior will actually be applied against, so it is the natural denominator. (If we kept §2.14's "post-policy buffer count," we would have to rescale α_gdna to match the EM input population — an extra source of error.) Amendment to §2.14 recorded inline (see §5.2.2.5).

One pandas → numpy hop, then no pandas in the inner loop. A pandas .loc[region_ids] per Locus is ~50× slower than a numpy fancy-index. The _arrays.py module materializes the four region columns we touch (start, end, type, boundary_flux_*) and the three payload reductions (per_region[:, 0b100], per_region[:, 0b010], (u_L, u_R)) once at entry.

Single global gDNA component per MultiLocus (v5.v1). A MultiLocus with len(.loci) > 1 (paralog clusters across refs) still has one gDNA component in the EM. Per-Locus diagnostics are retained on MultiLocusPrior.per_locus to enable the v5.v2 per-Locus-component spike (§5.2.4) without re-estimating.

5.2.1 Subphase 7a — native EM continuous prior_weight_rna
Goal. Replace the binary "include / exclude" nRNA suppression hack (zero-coverage shadowing) with a per-component float weight in [0, 1] honored by the native fan-out. v5.v1 uses {1.0 for mRNA, 0.0 for nRNA}; v5.v2 may set arbitrary continuous weights with no ABI change.

5.2.1.1 Native change (em_solver.cpp::compute_ovr_prior_and_warm_start)
The existing fan-out (em_solver.cpp lines 821–849) already computes coverage_totals[i] per component during the warm-start E-step. Patch is minimal and confined to two short blocks:

//  --- Existing signature gains one parameter ---
static void compute_ovr_prior_and_warm_start(
    ...,
    const float*  prior_weight_rna,   // [n_components], NEW; nullptr ≡ all-ones
    ...);

// --- Replace the total_rna_coverage accumulator (current line ~824) ---
double total_weighted_rna_coverage = 0.0;
int    n_rna_eligible              = 0;
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] > 0.0 && i != gdna_idx) {
        const double w = prior_weight_rna ? double(prior_weight_rna[i]) : 1.0;
        total_weighted_rna_coverage += w * coverage_totals[i];
        if (w > 0.0) ++n_rna_eligible;     // count ELIGIBLE-AND-WEIGHTED
    }
}

// --- Replace the per-component fan-out (current line ~838) ---
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] <= 0.0) {
        prior_out[i] = 0.0;
    } else if (i == gdna_idx) {
        prior_out[i] = baseline + std::max(alpha_gdna, EM_LOG_EPSILON);
    } else {
        const double w = prior_weight_rna ? double(prior_weight_rna[i]) : 1.0;
        if (total_weighted_rna_coverage > 0.0) {
            prior_out[i] = baseline + std::max(
                alpha_rna * w * coverage_totals[i] / total_weighted_rna_coverage,
                EM_LOG_EPSILON);
        } else if (n_rna_eligible > 0) {
            prior_out[i] = baseline + std::max(
                alpha_rna * w / n_rna_eligible, EM_LOG_EPSILON);
        } else {
            prior_out[i] = baseline;   // no weighted RNA mass → fall through
        }
    }
}
Invariants the patch preserves:

Passing nullptr (or a numpy None from Python) gives bit-identical output to the current implementation.
prior_weight_rna[gdna_idx] is ignored by construction (the branch is gated on i != gdna_idx). Callers may set it to any value; we recommend 0.0 for self-documentation.
A weight of exactly 0.0 collapses that component's prior to the baseline (0.0 for MAP, 0.5 for VBEM Jeffreys correction). It does not zero out theta_init — the warm-start is an empirical measurement and stays untouched.
Performance. The hot loop already iterates n_components. The patch adds a single prior_weight_rna ? float(...) : 1.0 ternary per iteration; on a typical locus (n_components ≤ 200) this is < 1 µs per warm-start. The branch predicts perfectly because the pointer is constant for the call.

5.2.1.2 Nanobind signature change
run_em_locus and run_em_batch (the Python entry points in em_solver.cpp lines 1761 / 1946) accept a new optional argument:

# Single-locus
run_em_locus(
    ...,
    alpha_gdna: float,
    alpha_rna:  float,
    prior_weight_rna: np.ndarray | None = None,  # float32, (n_components,)
    ...
)

# Batch
run_em_batch(
    ...,
    locus_alpha_gdna: np.ndarray,                # (n_loci,) float64
    locus_alpha_rna:  np.ndarray,                # (n_loci,) float64
    locus_prior_weight_rna: list[np.ndarray] | None = None,  # NEW
    ...
)
Both default to None for trivial backwards compatibility. The batch form takes a list (not a 2-D array) because per-locus n_components varies.

5.2.1.3 Python helper (locus.py)
def build_prior_weight_rna(
    multi_locus: MultiLocus,
    em_data:     ScoredFragments,
    nrna_weight: float = 0.0,        # v5.v1 default
) -> np.ndarray:
    """Return float32[n_components] in component-index order.

    Components per MultiLocus EM are laid out as
    ``[mRNA_t for t in transcripts] + [nRNA_t for t in transcripts] + [gDNA]``.
    This helper depends only on the layout, not on EM internals.
    """
5.2.1.4 Tests (tests/test_v5_em_prior_weight.py, ≥ 5 cases)
Bit-identical regression. Call run_em_locus with prior_weight_rna=None and again with prior_weight_rna=ones(n) on a 3-component scenario; outputs must be identical to machine precision.
Total suppression. weights = [1, 1, 0] on a [mRNA, nRNA, gDNA] locus → the nRNA component receives only baseline (0 in MAP, 0.5 in VBEM); mRNA / gDNA priors unchanged from the weight=[1,1,1] case after re-normalization.
Half-share. weights = [1, 1, 0.5] → nRNA prior is exactly half the contribution it would have under uniform weighting (after re-normalization).
Degenerate input. All-zero RNA weights with non-zero alpha_rna → fan-out skips RNA components entirely; output is well-formed (no NaN, no division by zero).
VBEM baseline preserved. Same as case 2 but with mode="vbem"; suppressed component receives exactly 0.5, not 0.
5.2.1.5 Exit gate (7a)
5 tests green; protected suite green; EM golden outputs bit-identical when caller passes None.
5.2.2 Subphase 7b — Per-Locus gDNA mass
The numerical core of Phase 7. Implements §2.7 verbatim, with performance and numerical-safety hardening described below.

5.2.2.1 Schemas
# density_loco.py — pure scalar EB shrinkage
def shrink_to_loco(
    n_loco:     int,        # observed count inside Locus
    leff_loco:  float,      # effective-length sum inside Locus
    rho_global: float,      # global density of this type
    kappa:      float,      # per-type NB overdispersion (Phase 6)
) -> float:
    """Section 2.6.1 closed-form NB EB shrinkage.

    Returns ``(n + κρ) / (L + κ)`` with all-zero / NaN guards.
    Edge cases:
        * leff_loco = 0 ⇒ returns rho_global (prior-only)
        * κ = ∞         ⇒ returns rho_global (full shrinkage)
        * κ = 0         ⇒ returns n / L      (no shrinkage)
    """


# locus_prior.py — value types
@dataclass(frozen=True, slots=True)
class LocusGdnaEstimate:
    locus:                Locus
    n_obs:                int          # EM-unit count inside Locus
    n_gdna_intergenic:    float
    n_gdna_intron:        float
    n_gdna_exon_intron:   float
    n_gdna:               float        # = sum of the three above
    pi_gdna:              float        # = n_gdna / max(n_obs, 1), clipped [0, 1]
    rho_loco:             tuple[float, float, float]   # (intergenic, intron, exon-intron) — diagnostic
    leff_loco:            tuple[float, float, float]   # (intergenic, intron, exon-intron eff-lengths)
    n_eligible_boundaries: int          # EXON-INTRON sides count
    fallback_flags:       int           # bitmask: see §5.2.2.7

@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    multi_locus_id:  int
    n_obs:           int                # = Σ per_locus n_obs
    n_gdna:          float
    n_rna:           float              # = max(0.0, n_obs − n_gdna)
    pi_gdna:         float              # = n_gdna / max(n_obs, 1)
    per_locus:       tuple[LocusGdnaEstimate, ...]
LocusGdnaEstimate lives on MultiLocusPrior.per_locus and survives into the Phase 8 per_locus_gdna_df. No information is discarded.

5.2.2.2 Pre-flattened arrays (_arrays.py)
A single pandas → numpy hop, owned per assemble_priors_v5 call:

@dataclass(frozen=True, slots=True)
class RegionArrays:
    starts:      np.ndarray   # int64 [n_regions]
    ends:        np.ndarray   # int64 [n_regions]
    types:       np.ndarray   # uint8 [n_regions]
    flag_left:   np.ndarray   # bool  [n_regions]
    flag_right:  np.ndarray   # bool  [n_regions]

    @classmethod
    def from_region_df(cls, region_df: pd.DataFrame) -> "RegionArrays": ...

@dataclass(frozen=True, slots=True)
class PayloadArrays:
    n_intergenic_per_region:  np.ndarray  # int64 [n_regions] — payload[:, 0b100]
    n_intron_per_region:      np.ndarray  # int64 [n_regions] — payload[:, 0b010]
    u_L:                      np.ndarray  # int64 [n_regions]
    u_R:                      np.ndarray  # int64 [n_regions]

    @classmethod
    def from_payload(cls, payload: CalibrationScanPayload) -> "PayloadArrays": ...
These are constructed once and threaded through every per-Locus call. They reference (do not copy) the underlying numpy arrays where possible; the from_payload slice path is payload.per_region_counts[:, 4] which is a view, not a copy, in C-order.

5.2.2.3 Per-Locus algorithm (estimate_locus_gdna)
def estimate_locus_gdna(
    locus:            Locus,
    region_index:     RegionIndex,
    region_arrays:    RegionArrays,
    payload_arrays:   PayloadArrays,
    global_densities: GlobalDensityTable,
    gdna_fl_mean:     float,
    n_obs_locus:      int,
) -> LocusGdnaEstimate:
    # 1. Native overlap (binary search; allocation-free in C++).
    rids = region_index.overlap(locus.ref_id, locus.start, locus.end)
    if rids.size == 0:
        raise RuntimeError(
            f"Locus {locus} has no overlapping calibration regions — "
            "BAM reference does not match index. Rebuild index against "
            "the BAM's reference."
        )

    # 2. Pull region geometry (numpy fancy-index — no pandas).
    rs = region_arrays.starts[rids]
    re = region_arrays.ends[rids]
    rt = region_arrays.types[rids]
    cl = np.minimum(re, locus.end) - np.maximum(rs, locus.start)        # clipped length
    leff = np.maximum(0, cl - gdna_fl_mean + 1).astype(np.float64)      # gDNA-effective

    fallback = 0      # bitmask collected as we go (see §5.2.2.7)

    # 3a. INTERGENIC contribution.
    is_inter = rt == 0
    n_inter      = int(payload_arrays.n_intergenic_per_region[rids[is_inter]].sum())
    leff_inter   = float(leff[is_inter].sum())
    rho_inter    = shrink_to_loco(
        n_inter, leff_inter,
        global_densities.intergenic.rho,
        global_densities.intergenic.kappa.value,
    )
    n_gdna_inter = rho_inter * leff_inter
    if leff_inter == 0.0 and is_inter.any():
        fallback |= FLAG_INTERGENIC_ZERO_LEFF

    # 3b. INTRON contribution.
    is_intron = rt == 1
    n_intron      = int(payload_arrays.n_intron_per_region[rids[is_intron]].sum())
    leff_intron   = float(leff[is_intron].sum())
    rho_intron    = shrink_to_loco(
        n_intron, leff_intron,
        global_densities.intron.rho,
        global_densities.intron.kappa.value,
    )
    n_gdna_intron = rho_intron * leff_intron
    if leff_intron == 0.0 and is_intron.any():
        fallback |= FLAG_INTRON_ZERO_LEFF

    # 3c. EXON-INTRON contribution (boundary flux).
    is_exon = rt == 2
    leff_exon = float(leff[is_exon].sum())     # multiplier for predicted count (§2.7.1 step 3)
    if is_exon.any():
        exon_rids = rids[is_exon]
        fl = region_arrays.flag_left [exon_rids]
        fr = region_arrays.flag_right[exon_rids]
        n_L = np.where(fl, payload_arrays.u_L[exon_rids], 0).sum()
        n_R = np.where(fr, payload_arrays.u_R[exon_rids], 0).sum()
        sides = int(fl.sum() + fr.sum())
        n_ei         = int(n_L + n_R)
        leff_ei_eb   = sides * gdna_fl_mean         # EB-shrinkage denom (§2.6 EXON-INTRON)
        if sides == 0:
            fallback |= FLAG_EXON_INTRON_NO_ELIGIBLE
        rho_ei = shrink_to_loco(
            n_ei, leff_ei_eb,
            global_densities.exon_intron.rho,
            global_densities.exon_intron.kappa.value,
        )
        n_gdna_ei = rho_ei * leff_exon              # multiply by FULL exonic L_eff (§2.7.1 step 3)
    else:
        rho_ei, n_gdna_ei, sides = 0.0, 0.0, 0

    # 4. Sum & finalize.
    n_gdna = n_gdna_inter + n_gdna_intron + n_gdna_ei
    if not math.isfinite(n_gdna):
        raise RuntimeError(
            f"Locus {locus}: non-finite predicted gDNA count "
            f"({n_gdna_inter}, {n_gdna_intron}, {n_gdna_ei}). "
            "EB κ underflow or shape mismatch."
        )

    pi_clip_warned = n_gdna > n_obs_locus
    pi_gdna = min(max(n_gdna / max(n_obs_locus, 1), 0.0), 1.0) if n_obs_locus > 0 else 0.0
    if pi_clip_warned:
        fallback |= FLAG_PI_CLIPPED
        logger.warning(
            "Locus %s: π̂_gdna > 1 before clipping (predicted=%.1f, n_obs=%d). "
            "Components: intergenic=%.1f, intron=%.1f, exon-intron=%.1f.",
            locus, n_gdna, n_obs_locus,
            n_gdna_inter, n_gdna_intron, n_gdna_ei,
        )

    return LocusGdnaEstimate(
        locus=locus, n_obs=n_obs_locus,
        n_gdna_intergenic=n_gdna_inter, n_gdna_intron=n_gdna_intron,
        n_gdna_exon_intron=n_gdna_ei, n_gdna=n_gdna, pi_gdna=pi_gdna,
        rho_loco=(rho_inter, rho_intron, rho_ei),
        leff_loco=(leff_inter, leff_intron, leff_exon),
        n_eligible_boundaries=sides, fallback_flags=fallback,
    )
Numerical hardening notes.

The EB shrinkage denominator leff + κ is always > 0 because κ ≥ KAPPA_MIN = 1.0 (Phase 6 contract). No divide-by-zero possible.
n_obs_locus = 0 is silent: π̂ = 0 (the prior collapses to all-RNA). This is the correct behavior for a Locus with no observed fragments.
A π̂ > 1 before clipping is logged at WARN, not raised: it is diagnostic of either (a) a Locus with strong coverage from unannotated genomic territory leaking past the boundary or (b) an EB-shrunk density that doesn't match the locus, which is exactly the failure mode the per-Locus diagnostics are designed to expose. Raising would prevent quant from completing on real-world edge cases.
Non-finite output is a hard failure (defensive guard against upstream EB pathologies).
5.2.2.4 MultiLocus aggregation (assemble_multilocus_prior)
def assemble_multilocus_prior(
    multi_locus:      MultiLocus,
    per_locus:        tuple[LocusGdnaEstimate, ...],
) -> MultiLocusPrior:
    n_obs  = sum(e.n_obs  for e in per_locus)
    n_gdna = sum(e.n_gdna for e in per_locus)
    return MultiLocusPrior(
        multi_locus_id=multi_locus.locus_id,
        n_obs=n_obs,
        n_gdna=n_gdna,
        n_rna=max(0.0, n_obs - n_gdna),
        pi_gdna=(n_gdna / max(n_obs, 1)) if n_obs > 0 else 0.0,
        per_locus=per_locus,
    )
For a typical 1-Locus MultiLocus, this is a four-line reduction. For multi-ref MultiLoci, n_rna = max(0, ...) is a defensive guard against the (rare) case where one Locus's predicted gDNA exceeds its observed fragments while another's is near zero — the sum can in principle exceed total obs without the per-Locus clip catching it.

5.2.2.5 Per-Locus n_obs derivation (amends §2.14)
§2.14 originally specified iterating the FragmentBuffer with the calibration exclusion policy. Phase 7b reinterprets this: n_obs(L) is the count of EM units assigned to Locus L. This is the denominator that the EM will use the prior against, so it is the correct denominator for π̂. (Using the calibration-policy buffer count would require rescaling α_gdna by n_em_units / n_calibration_obs, which introduces a free error term whenever the policy fractions differ between the global library and the local Locus.)

The implementation lives in _locus_n_obs.py:

def partition_units_to_loci(
    multi_locus: MultiLocus,
    em_data:     ScoredFragments,
    buffer:      FragmentBuffer,
    t_to_ref_id: np.ndarray,           # int32 [n_transcripts] — index.t_df["ref"].cat.codes
                                       # remapped via canonical ref_lengths order
) -> tuple[np.ndarray, ...]:
    """Return per-Locus arrays of EM unit indices.

    Fast path for ``len(multi_locus.loci) == 1`` (≈99 % of MultiLoci):
    returns ``(multi_locus.unit_indices,)`` directly.
    """
    if len(multi_locus.loci) == 1:
        return (multi_locus.unit_indices,)

    units = multi_locus.unit_indices                                  # int32
    fids  = em_data.frag_ids[units]                                   # int64
    t0    = em_data.locus_t_indices[units]                            # int32 — best transcript per unit
    rid   = t_to_ref_id[t0]
    gs    = buffer.genomic_start[fids].astype(np.int64)
    ge    = gs + buffer.genomic_footprint[fids].astype(np.int64)

    # For each unit, find the unique Locus interval whose
    # (ref_id, start, end) contains it. Loci are disjoint by
    # construction, so the assignment is unique (or the unit is
    # an orphan — see fallback below).
    out: list[list[int]] = [[] for _ in multi_locus.loci]
    orphans:   list[int] = []
    for k in range(units.size):
        for li, l in enumerate(multi_locus.loci):
            if l.ref_id == rid[k] and gs[k] >= l.start and ge[k] <= l.end:
                out[li].append(int(units[k]))
                break
        else:
            orphans.append(int(units[k]))

    if orphans:
        # Route orphans to the largest Locus (defensive; logged).
        # Realistic cause: chimeric multimapper whose best transcript's
        # ref doesn't match its genomic_start (synthesis artifact).
        biggest = max(range(len(multi_locus.loci)),
                      key=lambda i: multi_locus.loci[i].end - multi_locus.loci[i].start)
        out[biggest].extend(orphans)
        logger.debug("MultiLocus %d: %d orphan units routed to Locus %d",
                     multi_locus.locus_id, len(orphans), biggest)

    return tuple(np.array(x, dtype=np.int32) for x in out)
The fast path covers the vast majority of MultiLoci (the slow path runs only for paralog-cluster MultiLoci with > 1 Locus interval, which we estimate at < 1 % of MultiLoci on the human transcriptome).

The t_to_ref_id array is built once, outside the per-MultiLocus loop, in assemble_priors_v5:

canonical = {n: i for i, n in enumerate(index.ref_lengths.keys())}
t_to_ref_id = np.fromiter(
    (canonical[r] for r in index.t_df["ref"].cat.categories[
        index.t_df["ref"].cat.codes]),
    dtype=np.int32, count=index.num_transcripts,
)
5.2.2.6 Orchestration (assemble_priors_v5)
def assemble_priors_v5(
    multi_loci, em_data, buffer, index,
    payload, global_densities, gdna_fl_mean, c_base,
) -> PriorTable:
    region_arrays  = RegionArrays.from_region_df(index.regions)
    payload_arrays = PayloadArrays.from_payload(payload)
    t_to_ref_id    = _build_t_to_ref_id(index)

    multi_priors:    list[MultiLocusPrior] = []
    alpha_gdna_list: list[float]           = []
    alpha_rna_list:  list[float]           = []
    weight_list:     list[np.ndarray]      = []

    for ml in multi_loci:
        per_locus_units = partition_units_to_loci(ml, em_data, buffer, t_to_ref_id)
        per_locus_est = tuple(
            estimate_locus_gdna(
                locus=l, region_index=index.region_index,
                region_arrays=region_arrays, payload_arrays=payload_arrays,
                global_densities=global_densities, gdna_fl_mean=gdna_fl_mean,
                n_obs_locus=int(units.size),
            )
            for l, units in zip(ml.loci, per_locus_units, strict=True)
        )
        prior = assemble_multilocus_prior(ml, per_locus_est)
        multi_priors.append(prior)

        # Phase 7c: convert pi_gdna → Dirichlet pseudocounts.
        a_g, a_r = c_base.evaluate(prior, ml)
        alpha_gdna_list.append(a_g)
        alpha_rna_list.append(a_r)
        weight_list.append(build_prior_weight_rna(ml, em_data))

    return PriorTable(
        multi_locus_priors=tuple(multi_priors),
        alpha_gdna=np.array(alpha_gdna_list, dtype=np.float64),
        alpha_rna=np.array(alpha_rna_list,   dtype=np.float64),
        prior_weight_rna=weight_list,
    )
5.2.2.7 Performance contract & fallback bitmask
Per-Locus call budget (target on the human transcriptome, ~50 K MultiLoci × ~1.05 Loci/MultiLocus ≈ 53 K Locus calls):

Step	Budget per Locus	Justification
RegionIndex.overlap	< 5 µs	binary search + small forward walk in C++
numpy fancy-index	< 20 µs	4 arrays × O(k) where k ≈ 5 typical
EB shrinkage (3 scalar)	< 1 µs	pure scalar ops
Total	< 50 µs	per-Locus
End-to-end target: ≤ 3 s for the full human run. Above this we add the C++ port (deferred to v5.v3 unless profiling demands it).

fallback_flags bitmask (defined once in locus_prior.py):

FLAG_INTERGENIC_ZERO_LEFF   = 1 << 0   # had INTERGENIC region but L_eff = 0 (region too small)
FLAG_INTRON_ZERO_LEFF       = 1 << 1   # ditto INTRON
FLAG_EXON_INTRON_NO_ELIGIBLE = 1 << 2  # had EXON regions but no eligible boundaries (all terminal)
FLAG_PI_CLIPPED             = 1 << 3   # n_gdna > n_obs before clipping
These are diagnostic, not blocking. Phase 8 surfaces them in per_locus_gdna_df; Phase 10 validates aggregate counts on each benchmark scenario.

5.2.2.8 Tests
tests/test_v5_density_loco.py (≥ 6 cases, all pure-scalar, no fixtures):

n=0, L=0 → rho_global (prior-only edge).
n=100, L=10000, ρ=0.001, κ=10 → exact closed-form check.
Strong signal n=10000, L=10000 with small κ → ≈ n / L.
Weak signal n=2, L=10 with large κ → ≈ ρ_global.
Monotone in n for fixed L, ρ_global, κ.
Float64 stability: n=10⁹, L=10⁹ does not overflow.
tests/test_v5_locus_n_obs.py (≥ 5 cases):

Single-Locus MultiLocus → fast path returns (multi_locus.unit_indices,) (identity).
Two-Locus MultiLocus, units cleanly split between the two refs → each Locus gets exactly its share.
Orphan unit (synthetic: t0's ref doesn't match buffer genomic_start) → routed to the larger Locus, debug log emitted.
Empty MultiLocus → returns one empty array per Locus.
Σ per-Locus n_obs == len(multi_locus.unit_indices).
tests/test_v5_per_locus_gdna_mass.py (≥ 8 cases — pin §2.7.1 end-to-end with hand-built fakes for RegionIndex, RegionArrays, PayloadArrays, GlobalDensityTable):

Pure-mRNA Locus (zero per-region intergenic/intron counts, zero u_L/u_R) → n_gdna ≈ 0 (only the EB prior contributes; quantify bound).
Pure-gDNA Locus (counts ≈ density × L_eff) → pi_gdna ≈ 1 and each per-type predicted count recovers ground-truth within 20 %.
INTRON-only signal (no INTERGENIC, no boundary flux) → INTRON contribution dominates; other two are EB-prior only.
EXON-INTRON-only signal → EXON-INTRON contribution dominates; n_gdna_exon_intron = ρ̂_loco_ei × Σ L_eff(EXON in L).
No eligible EXON-INTRON boundaries (all terminal) → FLAG_EXON_INTRON_NO_ELIGIBLE set; rho falls back to global.
Locus partially overlapping a region → clipped L_eff is used (assert exact value).
n_obs = 0 → pi_gdna = 0, no warning.
n_gdna > n_obs → pi_gdna = 1.0 after clip, FLAG_PI_CLIPPED set, WARN logged.
No-overlap Locus raises with the canonical "rebuild index" message.
tests/test_v5_multilocus_aggregation.py (≥ 5 cases):

Single-Locus aggregation → totals == per-Locus values.
Two-Locus aggregation → totals == sum.
n_rna clamping when one Locus's predicted gDNA > its n_obs but the MultiLocus sum is still positive.
Per-Locus diagnostics retained on per_locus.
NaN per-Locus (synthetic) → assemble_multilocus_prior raises.
Performance test (tests/test_v5_per_locus_gdna_mass.py::test_perf):

Build a synthetic 5 000-MultiLocus, 200 000-region scenario; call assemble_priors_v5 once; assert wall-time < 1.0 s on CI. Marked @pytest.mark.benchmark; runs only with -m benchmark.
5.2.2.9 Exit gate (7b)
All 7b tests green; protected suite green.
Performance test passes (< 1 s for 5 K MultiLoci).
assemble_priors_v5 runs end-to-end on tests/scenarios_aligned/ with no RuntimeError.
Aggregate fallback-flag counts logged at INFO and recorded in the Phase 10 validation matrix.
5.2.3 Subphase 7c — c_base(ℓ) design spike + decision
Question. Phase 7b produces an intensity π̂(L) ∈ [0, 1]. The EM consumes Dirichlet pseudocounts (α_gdna, α_rna) with absolute mass interpretation. The conversion is:

α
gdna
(
L
)
=
c
base
(
L
)
⋅
π
^
(
L
)
,
α
rna
(
L
)
=
c
base
(
L
)
⋅
(
1
−
π
^
(
L
)
)
α 
gdna
​
 (L)=c 
base
​
 (L)⋅ 
π
^
 (L),α 
rna
​
 (L)=c 
base
​
 (L)⋅(1− 
π
^
 (L))

c_base(L) is the prior strength in fragment-equivalents. The Phase 0 transitional default was c_base = 10.0 (uniform 50/50 split). Phase 7c picks the right functional form.

5.2.3.1 Candidate formulas
All four candidates share the signature c_base(prior, multi_locus, em_data) -> float:

ID	Formula	Intuition
C1 — constant	c_base = C₀ (default 10.0)	Phase 0 baseline; uniform across loci
C2 — sqrt_n	c_base = C₀ + α · sqrt(n_obs)	Prior strength grows as √N → negligible vs data at high coverage; Bayesian reference-prior intuition
C3 — eb_inverse_var	c_base = κ̄ · n_obs / (n_obs + κ̄) where κ̄ = mean(κ_INTERGENIC, κ_INTRON, κ_EXON-INTRON)	Pseudocount = effective sample size of the EB prior at this Locus; falls naturally out of the NB hierarchy
C4 — coverage_weighted	c_base = clip(α · L_eff_locus, C_min, C_max)	Strength scales with effective coverage capacity; saturates to keep megaloci stable
C1 is the trivial baseline. C2/C4 are heuristic but parameter-free once α is fixed. C3 is the principled answer if you trust the Phase 6 κ estimates — and is automatically self-calibrating.

5.2.3.2 Spike harness (scripts/calibration/c_base_spike.py)
5 synthetic scenarios, all built from existing tests/scenarios/ infrastructure:

Scenario	What it stresses
pure_mrna	π̂ = 0 ⇒ α_gdna ≈ 0 must NOT siphon mRNA into nRNA
pure_gdna	π̂ = 1 ⇒ α_gdna must dominate; verify EM converges to all-gDNA
ta1_siphon	The canonical TA1 nRNA ratio failure (NTA1 ≫ TA1)
low_coverage	n_obs ≈ 5; the prior should dominate, not the data
mixed	Heterogeneous π̂ across loci; verify no one formula degrades
For each formula × scenario, run the full quant pipeline and record:

Mean / max relative error vs ground truth on (mRNA, nRNA, gDNA) per-pool.
Per-locus α_gdna / α_rna distribution (quartiles).
nRNA siphon magnitude (the TA1 metric).
Output: a TSV scoreboard + a markdown summary auto-appended to this plan as §2.8 once the decision is made.

5.2.3.3 Decision criteria
A candidate wins if it satisfies all of:

TA1 siphon < 5 % (the historical regression we are fixing).
Pure-mRNA error < 2 % (no spurious gDNA mass).
Pure-gDNA error < 2 % (prior doesn't under-claim either).
Low-coverage behaves: mean α_gdna among loci with n_obs < 10 sits between C1's value (10) and 0.5×n_obs (no silent prior blow-up).
Mixed scenario does not regress vs C1 on any pool by > 1 %.
If two candidates tie, pick the one with fewer tunable parameters (C3 > C2 > C4 > C1).

5.2.3.4 §2.8 amendment template
After the spike, replace §2.8 with:

### 2.8 `c_base(ℓ)` formula (RESOLVED — Phase 7c)

**Chosen:** <C-id>

**Formula:** <python-1-liner>

**Rationale:** <2–3 paragraphs from spike report>

**Diagnostics:** spike output committed to
`docs/calibration/c_base_spike_results.md`.
5.2.3.5 Tests
tests/test_v5_c_base.py (≥ 6 cases):

Each candidate is callable, returns a finite positive float on a trivial input.
C2: monotone non-decreasing in n_obs.
C3: returns 0 when n_obs = 0; returns ≈ κ̄ when n_obs ≫ κ̄.
C4: respects C_min, C_max clipping.
Dispatch (CBaseFormula.from_name) maps {"constant", "sqrt_n", "eb_inverse_var", "coverage_weighted"} to the right callable.
Unknown name raises ValueError with the list of valid names.
5.2.3.6 Exit gate (7c)
6 unit tests green.
Spike harness runs all 5 scenarios × 4 formulas successfully.
Decision recorded; §2.8 amended.
Selected formula passes all 5 decision criteria on at least 4 of 5 scenarios; documented carve-outs for any scenario that fails one criterion.
5.2.4 Subphase 7d — Per-Locus gDNA EM components (deferred to v5.v2)
A v5.v2 EM extension may emit one gDNA component per Locus interval inside a MultiLocus (rather than one shared component). Each gDNA component receives only fragments whose genomic footprint falls inside its Locus.

Why deferred: v5.v1 with a single shared gDNA component is correct (and competitive) for the > 99 % of MultiLoci with one Locus. The marginal benefit is concentrated in paralog clusters across refs — exactly the regime where Phase 5.5's multimapper exclusion already biases us. v5.v2 should pair this with the mappability-adjusted effective lengths (§9 open question) for a coherent "repeat-aware" release.

Design questions that the v5.v2 spike must resolve:

How are ambiguous fragments (no unique Locus assignment, e.g., chimera-spanning) routed? Plan: route to the gDNA component of the Locus containing the first mate; record a chimera_routed counter.
Do per-Locus gDNA components couple via a shared Dirichlet prior (one α_gdna split by Locus L_eff) or are they fully independent (each Locus gets its own α)? Independence is simpler but loses the global prior strength; coupling preserves it but requires the EM to do the split internally.
What does identifiability look like when one Locus has zero fragments? An independent-prior gDNA component with α_gdna > 0 and n_obs = 0 is identifiable (collapses to its prior mean) but uninformative; a coupled prior would let the other Loci absorb it.
The Phase 7b MultiLocusPrior.per_locus field carries enough information to support either design without re-estimating gDNA mass.

5.2.5 Phase 7 exit gate
Criterion	Subphase	Verification
Continuous prior_weight_rna ABI lands; bit-identical with None	7a	tests/test_v5_em_prior_weight.py + golden EM regression
Per-Locus estimator produces correct EB-shrunk densities	7b	tests/test_v5_per_locus_gdna_mass.py
MultiLocusPrior aggregates correctly	7b	tests/test_v5_multilocus_aggregation.py
No-overlap Locus raises with actionable message	7b	dedicated test case
Performance: 5 K MultiLoci in < 1 s	7b	@pytest.mark.benchmark test
Per-Locus diagnostics surfaced on MultiLocusPrior	7b	schema test
c_base formula chosen and §2.8 amended	7c	spike report; doc diff
Selected c_base meets 4/5 decision criteria	7c	spike scoreboard
TA1 nRNA siphon < 5 %	7b + 7c	spike ta1_siphon scenario
Pure-mRNA / pure-gDNA error < 2 %	7b + 7c	spike scenarios
Protected test set green	all	CI
Phase 7 does not require pipeline integration; that lands in Phase 9. The spike harness in 7c invokes the new estimator directly (behind a --use-v5-priors script flag) so Phase 9 can swap compute_default_locus_priors → assemble_priors_v5 mechanically.

5.2.6 Open questions & known limitations
Region double-counting between Loci within one genic span. A single calibration region (e.g., a long INTRON shared by a readthrough-transcript pair) can overlap two adjacent MultiLoci if no fragment links them. Current behavior: both Loci attribute the region's full per-region count to themselves. Realistic incidence: < 1 % of regions on the human transcriptome. v5.v2 fix: build a region → Locus assignment table at calibration time and split fractionally by clipped L_eff.
Single global gDNA component for multi-ref MultiLoci (§5.2.4 — Phase 7d).
c_base interaction with VBEM Jeffreys baseline. VBEM adds +0.5 to every component; the c_base spike is run in MAP mode. Re-run on VBEM scenarios to verify the chosen formula doesn't pathologically interact with the +0.5 baseline. Tracked in Phase 10.
Orphan units in partition_units_to_loci. A unit whose best transcript's ref disagrees with the buffer's genomic_start is routed to the largest Locus and logged. Realistic incidence: chimeric multimapper artifacts; expected to be < 0.01 % of units. If this exceeds 0.1 % on real data, escalate to an investigation.
gdna_fl_mean source during Phase 7 development. The bootstrap stub mean from Phase 0 is RNA-biased. Phase 8 swaps in the calibration-derived mean. Phase 7c spike must be re-run after Phase 8 lands to confirm the c_base decision survives the FL refresh; deviations > 5 % on any scenario reopen the formula choice.
5.3 Phase 8 — gDNA FL pool refresh + CalibrationResultV5 schema
Phase 8 is the consolidation phase. It does not change any algorithm that ships in v5.v1; it produces (a) a real gDNA fragment-length model derived from the calibration pool and (b) the canonical CalibrationResultV5 schema that Phase 9 will plumb through the pipeline. The schema is the single, frozen carrier for everything any downstream consumer (EM, summary JSON, diagnostic notebooks, future post-hoc audits) needs to know about a calibration run.

Phase 8 is split into three subphases:

8a — compute_pool_fl_models(): build the three FL models (gDNA, RNA, global) from CalibrationScanPayload.fl_hist plus a conservative quality classifier.
8b — CalibrationResultV5 dataclass + builder: assemble the schema from Phase 6 / Phase 7b / Phase 8a outputs, including serialization to summary.json plus optional feather emission of the per-MultiLocus / per-Locus diagnostic dataframes.
8c — Tests + skip-policy bookkeeping. No pipeline wiring yet (Phase 9 owns the wiring).
5.3.0 Architecture & data flow
CalibrationScanPayload      Phase 6 result          Phase 7b result
    ↓ (fl_hist)                 ↓                        ↓
compute_pool_fl_models   GlobalDensityTable      PriorTable + per-Locus
    ↓                         ↓                        ↓
gdna_fl_model            global_gdna_densities    multi_locus_prior_df
rna_fl_model                                       per_locus_gdna_df
global_fl_model
gdna_fl_quality
                           ↓ all three ↓
                  build_calibration_result_v5(...)
                              ↓
                     CalibrationResultV5  ← Phase 9 plumbs into the pipeline
compute_pool_fl_models and build_calibration_result_v5 are pure functions over their inputs — no pipeline knowledge, no I/O. They live in src/rigel/calibration/_fl_pool.py and src/rigel/calibration/_result.py respectively.

5.3.1 Subphase 8a — Pool FL models (_fl_pool.py)
Pool definition (canonical, locked in §2.3):

Pool	Source mask(s) on fl_hist	Semantics
gDNA pool	mask 2 (INTRON_ONLY) ∪ mask 3 (EXON|INTRON) ∪ mask 4 (INTERGENIC_ONLY)	The boundary-flux + intron-only + intergenic-only set; all three are gDNA-dominated by design
RNA pool	mask 1 (EXON_ONLY)	Strictly exonic; the only mask that is RNA-dominated
Global pool	sum over all 8 masks	Already produced by the BAM scan as frag_length_models.global_model (re-derived here for self-containment)
Annotation-gap masks (mask 5/6/7, all of which include INTERGENIC) are excluded from the gDNA pool by default (§2.3) because INTERGENIC overlaps under the current annotation pass tend to be either chimeric or artifactual; their counts go into n_pool_annotation_gap for diagnostics.

Schema:

@dataclass(frozen=True, slots=True)
class PoolFLModels:
    """Three FL models built from the calibration pool."""
    gdna_fl_model:   FragmentLengthModel
    rna_fl_model:    FragmentLengthModel
    global_fl_model: FragmentLengthModel
    gdna_fl_quality: Literal["good", "weak", "fallback"]
    n_pool:          int                       # total fragments in gDNA pool
    n_rna:           int                       # total fragments in RNA pool
    n_global:        int                       # total fragments globally
    n_pool_annotation_gap: dict[str, int]      # per-mask leakage counts (5/6/7)
Function signature:

def compute_pool_fl_models(
    fl_hist:        np.ndarray,                # shape (8, fl_bins), int64
    *,
    max_size:       int,                       # carries through to FragmentLengthModel
    prior_ess:      float = 1000.0,            # for low-count EB shrinkage
    quality_threshold_good:     int = 5_000,   # n_pool ≥ this → "good"
    quality_threshold_weak:     int = 200,     # n_pool ≥ this → "weak"; else "fallback"
) -> PoolFLModels:
    """Build (gdna, rna, global) FragmentLengthModel triple from the
    Phase 4 calibration ``fl_hist`` matrix.

    The three pools are built by simple mask-row summation; no shrinkage
    is performed when the pool count is healthy. When ``n_pool`` falls
    below ``quality_threshold_good`` the gDNA FL is built by EB-shrinking
    the pool histogram toward the global histogram with ``prior_ess``;
    when it falls below ``quality_threshold_weak`` the gDNA FL falls back
    to the global histogram entirely (and a warning is logged).
    """
Quality-classifier behavior (locked):

n_pool	gdna_fl_quality	gDNA FL build path
≥ 5000	"good"	Direct FragmentLengthModel.from_counts(pool_counts)
200–4999	"weak"	build_gdna_fl(pool_counts, global_counts, prior_ess=1000)
< 200	"fallback"	Direct FragmentLengthModel.from_counts(global_counts) (i.e., gDNA ≡ global)
The thresholds are conservative defaults; both are tunable via PipelineConfig once Phase 9 wires the orchestrator. The boundaries are picked so that the existing bootstrap behavior on the protected test scenarios maps to "weak" (the bootstrap currently uses prior_ess=1000), keeping the Phase 8 → Phase 9 transition smooth.

Tests (tests/test_v5_pool_fl.py, ≥ 8 cases):

Pool sums equal expected mask-row sums (round-trip on a hand-built fl_hist).
Annotation-gap masks (5/6/7) appear in n_pool_annotation_gap, not in n_pool.
n_pool ≥ quality_threshold_good ⇒ gdna_fl_quality == "good", no shrinkage applied (FL counts are exactly the pool counts).
n_pool ∈ [200, 5000) ⇒ gdna_fl_quality == "weak", EB shrinkage applied (FL mean is between pool mean and global mean).
n_pool < 200 ⇒ gdna_fl_quality == "fallback", gDNA FL counts exactly equal global counts.
Empty fl_hist (all zeros) ⇒ "fallback" + uniform FL.
max_size is propagated correctly to all three returned models.
Custom thresholds via kwargs override defaults.
5.3.2 Subphase 8b — CalibrationResultV5 schema + builder
Schema (frozen, slots, full-documented):

@dataclass(frozen=True, slots=True)
class CalibrationResultV5:
    """Calibration v5 frozen result. The single, canonical carrier of
    every artifact a downstream consumer may need."""

    # Discriminator
    version: str = "v5"
    active:  bool = True

    # FL models (consumed by FragmentScorer in quant_from_buffer)
    gdna_fl_model:   FragmentLengthModel
    rna_fl_model:    FragmentLengthModel
    global_fl_model: FragmentLengthModel
    gdna_fl_quality: str            # "good" | "weak" | "fallback"

    # Library-level signal (read by downstream EM as a diagnostic)
    strand_specificity: float

    # Region partition + global gDNA densities (Phase 6)
    region_df:             pd.DataFrame                   # passthrough from index.region_df
    global_gdna_densities: GlobalDensityTable             # 3 entries

    # Phase 7b prior tables — populated lazily after build_loci
    multi_locus_prior_df: pd.DataFrame | None = None      # one row per MultiLocus
    per_locus_gdna_df:    pd.DataFrame | None = None      # one row per Locus
    c_base_formula:       str = "constant"                # name from c_base.VALID_NAMES

    # Pool diagnostics
    pi_pool:                float = 0.0                   # n_pool / (n_pool + n_rna)
    n_pool:                 int = 0
    n_rna:                  int = 0
    n_global:               int = 0
    n_pool_annotation_gap:  dict[str, int] = field(default_factory=dict)

    # Scan-time exclusion counters (from CalibrationScanPayload)
    n_multimap_excluded:    int = 0
    n_chimera_excluded:     int = 0
    n_artifact_excluded:    int = 0
    n_oor_excluded:         int = 0

    # Category counts from the scan (8-element mask vector)
    category_counts:        np.ndarray = field(
        default_factory=lambda: np.zeros(8, dtype=np.int64))

    # Per-region boundary-flux summary (Phase 5 u_L/u_R aggregates)
    boundary_flux_gdna_summary: dict = field(default_factory=dict)

    # Per-type κ diagnostics (Phase 6)
    kappa_diagnostics: dict[str, KappaEstimate] = field(default_factory=dict)

    # Reserved for future use (e.g., post-hoc audit tags)
    extra: dict = field(default_factory=dict)

    def to_summary_dict(self) -> dict[str, Any]:
        """JSON-serializable summary for ``summary.json``. Skips the
        bulky pandas dataframes; emits only scalars + the two FL-mean
        and quality scalars."""
        ...

    def with_priors(
        self,
        multi_locus_prior_df: pd.DataFrame,
        per_locus_gdna_df:    pd.DataFrame,
        c_base_formula:       str,
    ) -> "CalibrationResultV5":
        """Return a new CalibrationResultV5 with the prior tables
        backfilled. Used by Phase 9 after build_loci + assemble_priors_v5."""
        return dataclasses.replace(
            self,
            multi_locus_prior_df=multi_locus_prior_df,
            per_locus_gdna_df=per_locus_gdna_df,
            c_base_formula=c_base_formula,
        )
Schema design notes:

Lazy prior tables. multi_locus_prior_df and per_locus_gdna_df default to None. They are filled in by with_priors() after build_loci runs in quant_from_buffer. This keeps the global Phase 8 calibration step decoupled from EM-time work — calibrate_v5 can return a usable result before any locus has been constructed.
region_df is a passthrough, not a copy. The dataframe is already in index.region_df; storing the same reference here lets downstream consumers ask the calibration result for the schema-locked view without an extra index argument.
No per-region density columns (§2.12 — explicitly removed).
c_base_formula is the name from c_base.VALID_NAMES, not the callable; the callable is reconstructable via from_name.
with_priors uses dataclasses.replace, not mutation: the schema stays frozen.
The schema deliberately does NOT include:
Any reference to CalibrationStub / pi_pool == 0.5 — this is real calibration, not a placeholder.
mixture_converged / mixture_iterations — these belong to a different mixture-EM design that is not part of v5.v1. The plan's earlier mention of these fields is removed in this rewrite.
Builder signature:

def build_calibration_result_v5(
    pool_models:           PoolFLModels,
    global_densities:      GlobalDensityTable,
    payload:               CalibrationScanPayload,
    region_df:             pd.DataFrame,
    strand_specificity:    float,
) -> CalibrationResultV5:
    """Compose Phase 6 + Phase 8a outputs into a CalibrationResultV5.

    Per-MultiLocus prior tables are not populated here — they are
    backfilled by Phase 9's quant_from_buffer via ``with_priors``.
    """
The builder is a thin assembler: it sums payload.global_counts into category_counts, computes pi_pool = n_pool / (n_pool + n_rna), populates the κ diagnostics dict from global_densities.by_type(), aggregates payload.u_L / payload.u_R per-mask into boundary_flux_gdna_summary, and constructs the dataclass.

Dataframe schemas (locked):

multi_locus_prior_df — one row per MultiLocus, in EM iteration order:

Column	Dtype	Source
multi_locus_id	int64	MultiLocus.locus_id
n_loci	int32	len(MultiLocus.loci)
n_obs	int64	MultiLocusPrior.n_obs
n_gdna	float64	MultiLocusPrior.n_gdna
n_rna	float64	MultiLocusPrior.n_rna
pi_gdna	float64	MultiLocusPrior.pi_gdna
alpha_gdna	float64	PriorTable.alpha_gdna[i]
alpha_rna	float64	PriorTable.alpha_rna[i]
c_base	float64	alpha_gdna + alpha_rna (sanity / introspection)
per_locus_gdna_df — one row per Locus (i.e., Σ len(ml.loci) rows, multiple per MultiLocus when a MultiLocus has > 1 Locus):

Column	Dtype	Source
multi_locus_id	int64	parent MultiLocus
locus_idx	int32	index of this Locus inside multi_locus.loci
ref	str	Locus.ref
start	int64	Locus.start
end	int64	Locus.end
n_obs	int64	LocusGdnaEstimate.n_obs
n_gdna	float64	LocusGdnaEstimate.n_gdna
n_gdna_intergenic	float64	per-type predicted count
n_gdna_intron	float64	per-type predicted count
n_gdna_exon_intron	float64	per-type predicted count
pi_gdna	float64	per-Locus clipped ratio
rho_loco_intergenic	float64	EB-shrunk locoregional density
rho_loco_intron	float64	EB-shrunk locoregional density
rho_loco_exon_intron	float64	EB-shrunk locoregional density
leff_intergenic	float64	per-type effective length sum
leff_intron	float64	per-type effective length sum
leff_exon	float64	per-type effective length sum
n_eligible_boundaries	int32	EXON-INTRON eligible-side count
fallback_flags	int32	bitmask (see §5.2.2.7)
Both dataframes are constructed by helper functions (_build_multi_locus_prior_df, _build_per_locus_gdna_df) in src/rigel/calibration/_result.py. They live next to the dataclass because the dataclass owns the schema contract.

5.3.3 Subphase 8c — Tests
tests/test_v5_calibration_result.py (≥ 10 cases):

Smoke: build_calibration_result_v5 on hand-built fakes returns a CalibrationResultV5 with version == "v5", active is True, prior dataframes None.
to_summary_dict returns JSON-serializable dict (round-trip via json.dumps); contains version, gdna_fl_quality, pi_pool, n_pool, n_rna, gdna_fl_mean, the four exclusion counters, and the κ diagnostics.
to_summary_dict does NOT serialize the dataframes (size guard).
with_priors(...) returns a new instance with both dataframes populated; original instance is unchanged (frozen contract).
pi_pool = n_pool / (n_pool + n_rna); zero-input case (n_pool == n_rna == 0) yields pi_pool == 0.0, no division by zero.
category_counts is the 8-element sum from payload.global_counts.
Schema-locked DataFrame columns: every column listed in §5.3.2 appears with its declared dtype after _build_multi_locus_prior_df on a hand-built PriorTable.
Same for _build_per_locus_gdna_df.
Dataframe row counts: len(multi_locus_prior_df) == n_multi_loci and len(per_locus_gdna_df) == sum(len(ml.loci) for ml in ...).
kappa_diagnostics is keyed by "INTERGENIC", "INTRON", "EXON-INTRON" and each value is a KappaEstimate.
Combined with test_v5_pool_fl.py (≥ 8 from 5.3.1), Phase 8 ships ≥ 18 new unit tests.

5.3.4 Exit gate (Phase 8)
#	Criterion	Verification
1	compute_pool_fl_models callable; ≥ 8 unit tests green	pytest
2	CalibrationResultV5 callable; ≥ 10 unit tests green	pytest
3	to_summary_dict round-trips through json.dumps	unit test
4	with_priors returns a new instance; original is frozen	unit test
5	Both dataframe schemas exactly match §5.3.2 column tables	unit test
6	Public exports added to src/rigel/calibration/__init__.py; test_v5_phase0_cleanup.py::test_package_imports_cleanly updated	regression
7	Protected suite green (no new regressions)	full pytest
8	No pipeline wiring changes (Phase 9 owns wiring)	code review
Phase 8 must not modify pipeline.py, quant_from_buffer, or bootstrap_fl_calibration. The new schema is built but not yet consumed.

5.4 Phase 9 — Pipeline integration
Phase 9 is the wiring step. It replaces the Phase 0 bootstrap path inside quant_from_buffer with the real calibrate_v5() orchestrator, threads the per-MultiLocus prior table through to the EM (closing the loop on Phase 7), and re-enables every test currently skipped behind until_phase="9".

Phase 9 is split into four subphases:

9a — calibrate_v5() orchestrator (composing 8a + Phase 6 + Phase 8b builder).
9b — quant_from_buffer rewiring: bootstrap → real calibration, flat priors → assemble_priors_v5, with the prior tables backfilled via with_priors.
9c — Backwards-compat policy: when index.region_df is None (older indexes built before the v5 region table), fall back to the bootstrap stub with a clear log line.
9d — Test re-enablement, golden regeneration, and tests/skip_during_v5_dev.txt cleanup.
5.4.0 Architecture & data flow
            scan_and_buffer()  ──┐
                                  ├──► (buffer, payload, fl_models, strand_models)
                index             │
                  │               │
                  ▼               ▼
            calibrate_v5(buffer, index, payload, fl_models, strand_models)
                              │
                              ├── compute_pool_fl_models(payload.fl_hist) ─┐
                              ├── compute_global_densities(...)            │
                              └── build_calibration_result_v5(...) ────────┘
                                       │
                                       ▼
                              CalibrationResultV5  (priors=None)
                                       │
                                       ▼
              quant_from_buffer(buffer, index, calibration=cal, ...)
                  ↓
              _score_fragments → em_data
                  ↓
              build_loci → multi_loci
                  ↓
              assemble_priors_v5(multi_loci, em_data, buffer, index,
                                 payload, global_densities,
                                 c_base=cal.c_base_formula or default)
                  ↓
              prior_table : PriorTable
                  ↓
              cal2 = cal.with_priors(_build_multi_locus_prior_df(...),
                                     _build_per_locus_gdna_df(...),
                                     c_base_formula=name)
                  ↓
              run_batch_locus_em_partitioned(
                  ..., prior_table.alpha_gdna, prior_table.alpha_rna,
                  prior_weight_rna_per_locus=prior_table.prior_weight_rna)
Notable structural decisions:

calibrate_v5 is global-only. It returns a result whose prior tables are None. Per-locus work happens after build_loci inside quant_from_buffer; the result is then immutably backfilled via with_priors. This keeps calibrate_v5 testable in isolation without needing locus or EM-data dependencies.
quant_from_buffer becomes cal-aware. It accepts a CalibrationResultV5 or a CalibrationStub (legacy). When it sees the v5 result it routes through the new prior path; when it sees the stub it falls back to compute_default_locus_priors. This lets us switch the default in run_pipeline without breaking any external script that constructs CalibrationStub directly.
compute_default_locus_priors is preserved, not deleted, for the backwards-compat path. It is annotated as deprecated and a future Phase 11 can remove it once external callers migrate.
5.4.1 Subphase 9a — calibrate_v5() orchestrator
Replace the NotImplementedError stub in src/rigel/calibration/__init__.py with a real implementation. Function lives in a new file src/rigel/calibration/_orchestrator.py to keep __init__.py thin.

def calibrate_v5(
    buffer:             FragmentBuffer,         # not used in v5.v1; reserved for v5.v2
    index:              TranscriptIndex,
    payload:            CalibrationScanPayload,
    frag_length_models: FragmentLengthModels,
    strand_models:      StrandModels,
    *,
    pool_quality_thresholds: tuple[int, int] = (5_000, 200),
) -> CalibrationResultV5:
    """Phase 9 calibration orchestrator (Calibration v5).

    Composes:
      1. compute_pool_fl_models(payload.fl_hist)            (Phase 8a)
      2. compute_global_densities(index.region_df, payload, gdna_fl_mean)  (Phase 6)
      3. build_calibration_result_v5(...)                   (Phase 8b)

    Per-MultiLocus prior tables are NOT populated here; they are
    backfilled by quant_from_buffer after build_loci.

    Raises ValueError if index.region_df is None (no calibration
    artifacts in this index — caller should fall back to bootstrap).
    """
    if index.region_df is None:
        raise ValueError(
            "calibrate_v5 requires an index built with calibration "
            "artifacts (index.region_df is None). Either rebuild the "
            "index, or pass a CalibrationStub via "
            "bootstrap_fl_calibration."
        )

    pool = compute_pool_fl_models(
        payload.fl_hist,
        max_size=frag_length_models.max_size,
        quality_threshold_good=pool_quality_thresholds[0],
        quality_threshold_weak=pool_quality_thresholds[1],
    )
    global_dens = compute_global_densities(
        index.region_df, payload,
        gdna_fl_mean=pool.gdna_fl_model.mean,
    )
    return build_calibration_result_v5(
        pool_models=pool,
        global_densities=global_dens,
        payload=payload,
        region_df=index.region_df,
        strand_specificity=strand_models.strand_specificity,
    )
Buffer is reserved for v5.v2 (mappability-adjusted L_eff) but plumbed through now to lock the signature.

Tests (tests/test_v5_orchestrator.py, ≥ 6 cases):

Smoke: callable on a synthetic minimal payload returns a CalibrationResultV5 with prior_tables == None.
index.region_df is None ⇒ ValueError with actionable message.
The returned result's gdna_fl_model.mean matches pool_models.gdna_fl_model.mean (no aliasing surprise).
The returned result's global_gdna_densities.gdna_fl_mean is pool_models.gdna_fl_model.mean (Phase 6 / Phase 8a coupling).
pool_quality_thresholds kwarg propagates to compute_pool_fl_models.
strand_specificity is taken from strand_models.strand_specificity.
5.4.2 Subphase 9b — quant_from_buffer rewiring
Three precise edits to src/rigel/pipeline.py:

Edit 1 — type-relax the calibration parameter:

# In quant_from_buffer, the calibration parameter is widened.
# Both CalibrationStub and CalibrationResultV5 expose
# `gdna_fl_model` and a `version` attribute; the dispatch is by
# isinstance check at the prior-construction site below.
calibration: CalibrationStub | CalibrationResultV5
Edit 2 — replace the prior-construction block:

# Old (Phases 0–8):
alpha_gdna, alpha_rna = compute_default_locus_priors(loci)

# New (Phase 9):
prior_weight_rna_per_locus: list | None = None
if isinstance(calibration, CalibrationResultV5):
    # Real calibration path — build per-MultiLocus priors.
    prior_table = assemble_priors_v5(
        multi_loci=loci,
        em_data=em_data,           # Note: still alive at this point
        buffer=buffer,
        index=index,
        payload=calibration_payload,   # NEW kwarg; see Edit 3
        global_densities=calibration.global_gdna_densities,
        # Phase 7c: c_base formula. Default to "constant" until 7c
        # spike picks a winner; the choice is recorded on the result.
        c_base=C_BASE_DEFAULT,
        nrna_weight=1.0,           # Phase 7a knob; default = legacy
    )
    alpha_gdna = prior_table.alpha_gdna
    alpha_rna  = prior_table.alpha_rna
    prior_weight_rna_per_locus = prior_table.prior_weight_rna

    # Backfill the prior dataframes onto the calibration result.
    calibration = calibration.with_priors(
        multi_locus_prior_df=_build_multi_locus_prior_df(
            loci, prior_table),
        per_locus_gdna_df=_build_per_locus_gdna_df(
            loci, prior_table),
        c_base_formula="constant",
    )
else:
    # Legacy stub path — symmetric flat priors.
    alpha_gdna, alpha_rna = compute_default_locus_priors(loci)
Edit 3 — thread payload and the per-locus weight vector through the EM call:

_run_locus_em_partitioned(
    estimator,
    partitions,
    loci,
    index,
    alpha_gdna, alpha_rna,
    em_config,
    emit_locus_stats=emit_locus_stats,
    annotations=annotations,
    prior_weight_rna_per_locus=prior_weight_rna_per_locus,  # NEW
)
_run_locus_em_partitioned already accepts the kwarg via Phase 7a; this is a pure passthrough.

run_pipeline switch — the default calibration construction moves from bootstrap_fl_calibration to calibrate_v5:

# Old (line ~884):
calibration = bootstrap_fl_calibration(frag_length_models)

# New:
if cal_payload is not None and index.region_df is not None:
    calibration = calibrate_v5(
        buffer=buffer, index=index, payload=cal_payload,
        frag_length_models=frag_length_models,
        strand_models=strand_models,
    )
else:
    # Older indexes without region_df → bootstrap fallback.
    logger.warning(
        "[CAL] No calibration payload or region_df available "
        "(index built before v5); falling back to bootstrap stub."
    )
    calibration = bootstrap_fl_calibration(frag_length_models)
PipelineResult.calibration type widens to CalibrationStub | CalibrationResultV5.

Tests (tests/test_v5_pipeline_integration.py, ≥ 8 cases):

End-to-end on a synthetic scenario: PipelineResult.calibration is CalibrationResultV5, version == "v5", active is True.
multi_locus_prior_df is non-None and has one row per MultiLocus.
per_locus_gdna_df is non-None and has the expected row count.
The EM ran with the new (alpha_gdna, alpha_rna) (i.e., not the flat 5/5 split).
Backwards-compat: pass a CalibrationStub directly to quant_from_buffer ⇒ legacy path runs (alpha_gdna/alpha_rna come from compute_default_locus_priors).
Older index (no region_df) ⇒ run_pipeline warns and falls back to bootstrap; result has version == "stub".
The summary JSON for a v5 run contains "version": "v5" and a non-zero pi_pool on a contaminated scenario.
n_pool in the result equals the sum of the gDNA-pool masks on the calibration payload.
5.4.3 Subphase 9c — Backwards-compat policy
Lock down explicit acceptance criteria for the legacy paths:

Any caller (script, notebook, external pipeline) that constructs CalibrationStub directly continues to work.
Any caller that calls compute_default_locus_priors directly continues to work (the function is undeleted).
Older indexes without region_df (built by rigel index versions before the Phase 1 region partition landed) keep working via the bootstrap fallback. A clear [CAL] No calibration payload warning is logged exactly once per run_pipeline invocation.
bootstrap_fl_calibration retains its docstring, but a DeprecationWarning is added: "Direct use of bootstrap_fl_calibration is supported only for legacy indexes. Prefer rebuilding the index with v5 calibration artifacts and using calibrate_v5."
5.4.4 Subphase 9d — Test re-enablement & golden regeneration
tests/skip_during_v5_dev.txt lists four test files currently blocked by until_phase="9":

Test file	Blocking reason	Phase 9 action
tests/test_golden_output.py	Goldens were captured against deleted SRD v1 outputs	Run the suite with pytest --update-golden; manually inspect the diff; commit new goldens with a git add -p review pass
tests/scenarios/test_diagnostics.py	Unspliced / low-strand assertions need real calibration	Re-enable; verify all assertions pass
tests/scenarios/test_non_overlapping.py	Negative-control gDNA accuracy assertions need real calibration	Re-enable; verify gDNA pi_pool lands in the expected range
tests/scenarios/test_nrna_double_counting.py	nRNA recovery requires real calibration	Re-enable; verify TA1 siphon < 5 % per the criterion in §5.2.3.3
tests/skip_during_v5_dev.txt is emptied of the v5 entries (or deleted entirely if no other entries remain).

The golden-regeneration step is delicate. The procedure:

Run the full protected suite once to confirm pre-Phase-9 state.
Apply Phase 9 wiring.
Run pytest tests/test_golden_output.py --update-golden.
git diff --stat tests/golden/ to review the magnitude of the diffs.
Read the diffs file-by-file. Every numeric change must trace back to a calibration design decision, not a bug. Document the net direction (e.g., "mRNA estimates ↑ ~3 % across all scenarios because real calibration suppresses gDNA siphon") in the commit message.
Commit goldens.
5.4.5 Exit gate (Phase 9)
#	Criterion	Verification
1	calibrate_v5() returns a valid CalibrationResultV5; ≥ 6 unit tests	pytest
2	quant_from_buffer end-to-end on a v5 index produces CalibrationResultV5 with both prior tables populated	pytest
3	Backwards compatibility: CalibrationStub path still works	pytest
4	Older-index fallback: index.region_df is None ⇒ bootstrap path runs with a single warning	pytest
5	tests/skip_during_v5_dev.txt v5 entries removed	grep
6	All four previously-skipped test files green	pytest
7	Golden outputs regenerated with documented diffs	git review
8	Full suite green (no skips except those documented for non-v5 reasons)	pytest
9	Smoke benchmark runs on tests/scenarios/ without regression vs the bootstrap baseline on RNA-only scenarios; demonstrates ≥ 80 % reduction in TA1 siphon on the contaminated scenario	benchmark
10	summary.json contains the v5 schema fields (version, pi_pool, gdna_fl_quality, etc.)	unit test
Phase 9 is the final wiring phase. Phase 10 is validation (real data, hybrid capture, full benchmarks); Phase 11 is documentation and release. After Phase 9 lands, the calibrator's behavior on the protected suite is locked.

5.4.6 Known follow-ups punted to later phases
c_base_formula selection. Phase 9 ships with c_base="constant" (the Phase 0 baseline). Promoting the Phase 7c spike's chosen formula into the orchestrator is a one-line change in quant_from_buffer and should be done immediately after the spike decision is recorded in §2.8.
nrna_weight selection. Phase 9 ships with nrna_weight=1.0 (legacy all-ones weights). The TA1 siphon mitigation that justified the Phase 7a knob is a Phase 10 follow-up: validate against the TA1 scenario in tests/scenarios/test_nrna_double_counting.py and pick the value that drives siphon < 5 %.
PipelineConfig exposure of pool_quality_thresholds and nrna_weight. These are constructor kwargs in Phase 9; Phase 11 promotes them to fields on PipelineConfig for end-user control.
Removal of compute_default_locus_priors and CalibrationStub. Targeted for Phase 11 once all external callers migrate.
5.5 Phase 10 — Validation harness + benchmarks
Validation matrix:

Pristine RNA-seq: pi_pool ≈ 0.
50/50 dna20m: pi_pool ∈ [0.45, 0.55].
TA1 siphon: residual < 5%.
Hybrid-capture sample: gDNA flux density on internal exons matches synthetic ground truth within ±15% on a low-contamination scenario.
Full benchmark vs salmon / kallisto on the cluster matrix.
5.6 Phase 11 — Documentation + release
Update docs/MANUAL.md, docs/METHODS.md, docs/parameters.md.
Bump version to 0.6.0.
CHANGELOG entry summarizing v5 design changes since v4.
6. Code-quality contracts (all phases)
6.1 C++ style
C++17, namespaces under rigel:: and rigel::calibration::.
No raw new/delete; smart pointers for ownership; references for non-owning views.
No heap allocation in CalibrationAccumulator::observe() after warm-up.
Public C++ types documented with doxygen-style comments.
All hot-path counters are int64_t.
6.2 SmallRegionSet contract
Inline storage kInline = 16.
Linear-scan dedup (N small).
Spilled vector reused across fragments.
Unit-tested for inline, overflow, dedup, clear.
6.3 Python style
Frozen dataclasses for all result objects.
Type hints on all public functions.
pandas DataFrames; no Polars.
No assert in production paths; raise explicit exceptions with actionable messages.
pathlib.Path for all file paths.
6.4 Diagnostics
Every phase emits structured INFO logs summarizing key counts.
All counts surface in CalibrationResultV5 and through summary.json.
Fallback paths are explicit and recorded (e.g., KappaEstimate.fallback_used).
7. Protected test set
The following must pass after every phase boundary; they cover end-to-end pipeline behavior independent of calibration intelligence:

tests/test_index.py
tests/test_index_integrity.py
tests/test_buffer.py
tests/test_resolution.py
tests/test_em_impl.py
tests/test_estimator.py
tests/test_pipeline_smoke.py
tests/test_oracle_bam.py
tests/test_golden_output.py
tests/test_cli.py
Tests legitimately blocked by missing v5 calibration get pytest.mark.skip(reason="...", until_phase="9"). CI fails on any skip lacking an until_phase= annotation or pointing at a phase ≤ the current phase.

8. Risk register
#	Risk	Phase	Likelihood	Mitigation
1	Phase 0 bootstrap FL too inaccurate; quant outputs unusable during phases 0–8	0	Medium	Bootstrap uses global FL shrunk toward RNA — biased toward calling fragments mRNA. Acceptable for protected smoke tests; flagged loudly in logs
2	Native nRNA suppression introduces bit-level diffs in non-nRNA loci	7a	Medium	Bit-identical regression on TA1-free scenario before merge
3	c_base decision regresses high-coverage loci	7c	Medium	5-locus design spike + decision recorded in §2.8
4	Annotation-gap masks more frequent than expected on real data	4 (done)	Medium	Per-mask QC counters; revisit if > 1% of pool
5	estimate_kappa MoM hits degenerate input on real data	6	Medium	Explicit fallback + diagnostics
6	SmallRegionSet overflows hot path on dense annotations	4 (done)	Low	Spilled vector reused; unit-tested
7	Stale-index migration confuses users	2 (done), 11	Low	Clear error + CHANGELOG + MANUAL update
8	prior_weight_rna ABI change for v5.v2	future	Low	Already float; no future ABI change required
9	Dropping spliced flux underestimates gDNA on capture-edge / FL-miscalibrated loci	5 (done), 10	Low	Boundary-flux extrapolation is FL-divisor-based; if a regression appears, re-introduce 
s
L
s 
L
​
  / 
s
R
s 
R
​
  behind a config flag — additive change
10	Terminal-exon exclusion costs effective boundary count on whole-genome libraries	5 (done)	Low	Loss is at most 2 boundaries per transcript; quantify on the validation matrix in Phase 10
11	Memory pressure from per-worker per_region_counts_ (~256 MB × n_workers worst case)	4 (done)	Low	Acceptable at planned n_workers ≤ 8. If real-world traces show pressure, switch to a streaming-merge (workers emit chunks; single shared buffer accumulates under a mutex). Not blocking
12	EXON-INTRON global density biases toward genes with many internal exons (more eligible boundaries)	6	Medium	Each eligible boundary contributes equally to numerator AND denominator (sides * L̄_gDNA); per-gene weighting cancels in the ratio. Quantify on Phase 10 hybrid-capture validation row
13	Locus rename in 6.5 breaks downstream callers (scripts, notebooks, tests)	6.5	Medium	LSP-driven vscode_renameSymbol; keep MultiLocus.merged_intervals as a deprecated property for one release; CI grep guard against bare Locus( constructions outside locus.py
14	Locoregional EB inflates κ-shrinkage when a Locus has zero overlapping regions of a given type	7b	Low	The numerator is 0 + κ·ρ_global and denominator 0 + κ → density collapses to ρ_global; degenerate but correct (uses the prior). Validated in test_v5_density_loco.py
15	Single MultiLocus-level gDNA component (v5.v1) misallocates gDNA across disjoint Locus intervals (e.g., pseudogene paralogs on different chromosomes)	7b	Medium	Per-Locus diagnostics surfaced; Phase 7d v5.v2 spike addresses the underlying identifiability question
16	Excluding multimappers under-estimates gDNA density on multi-mappable regions (intergenic repeats, pseudogenes)	5.5 (locked)	High	Documented + accepted for v5.v1; v5.v2 mappability-adjusted effective lengths 
L
eff
map
(
R
)
=
∫
R
m
(
x
)
d
x
L 
eff
map
​
 (R)=∫ 
R
​
 m(x)dx recover unbiased density without re-introducing multimappers (§2.13.1)
17	Per-region kappa MoM instability when one region has runaway count (Var dominated by one term)	6	Medium	Excess-variance is non-robust by construction; if instability observed on real data, swap to per-region NB MLE (single-parameter Newton on 
∑
log
⁡
L
∑logL). Tracked as a Phase 10 follow-up
18	region_cgr legacy property silently used by new calibration code	6.5+	Medium	Calibration code paths import rigel.native.RegionIndex only; CI grep guard against region_cgr references in src/rigel/calibration/**
9. Open questions reserved for live iteration
Phase 7c — c_base(ℓ) formula. Empirical decision after spike; amend §2.8.
Phase 7d — per-Locus EM gDNA components. Design spike for v5.v2; default off; resolution gates on the three identifiability questions in §5.2.4.
v5.v2 — mappability-adjusted effective lengths. Compensate for multimapper exclusion (§2.13.1) by replacing 
L
(
R
)
L(R) with 
∫
R
m
(
x
)
d
x
∫ 
R
​
 m(x)dx everywhere 
L
eff
L 
eff
​
  is computed. Requires a per-bp mappability track keyed to the read-length configuration; precompute at index time.
Annotation-gap inclusion in FL pool. Default False; flip to True only if QC counts justify and oracle validation supports it.
Mean RNA FL source for boundary-flux normalization. Default gdna_fl_model.mean() (Phase 6 input); fall back to global_fl_model.mean() if gdna_fl_quality == "fallback".
10. References
Earlier design notes (informative, not authoritative): /memories/repo/rigel-calibration-design.md.
User TODO original sketch: docs/TODO.md §"Calibration".
Archived predecessors:
docs/archive/calibration/calibration_v5_plan_v{1,2,3,4,5,6}.md
docs/archive/calibration/v5_phase4_handoff.md
docs/archive/calibration/srd_v2_phase6_handoff.md
docs/archive/calibration/srd_v3_phase2_oracle_findings.md
External reviews folded in: GPT 5.4 second-pass contracts review; Gemini critique (docs/calibration/gemini_critique.md).