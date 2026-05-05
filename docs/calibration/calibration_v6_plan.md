# Rigel Calibration — Plan

**Status:** AUTHORITATIVE. This is the single source of truth for the
calibration subsystem. There is no "v5", no "v6", no "phase X" terminology
inside the implementation: this document is named `calibration_v6_plan.md`
purely as filesystem provenance. Symbols, modules, tests, and comments
inside the codebase use the unmarked names (`calibrate`, `CalibrationResult`,
`assemble_priors`, etc.). If a future redesign supersedes this plan, that
redesign will retire the entire surface, not annotate it with a version
suffix.

**Date:** 2026-05-05
**Repo state at plan inception:** HEAD `a0dac2d` ("working on new calibration
system"). The `src/rigel/calibration/` package contains the legacy v1 SRD
code (`_simple.py`, `_categorize.py`, `_result.py`, `_fl_mixture.py`,
`_fl_empirical_bayes.py`). All Phase 0–5 work from the previous attempt
was lost in the 2026-05-01 data-loss incident
(`/memories/repo/data-loss-incident-2026-05-01.md`).

**Predecessor (informative, not authoritative):**
`docs/calibration/calibration_v5_plan.md`. Most design decisions in this
plan are inherited verbatim from that document. The differences are noted
in §0.

---

## 0. Differences from the predecessor plan

| Topic | Predecessor (v5 plan) | This plan |
|---|---|---|
| Versioning in code | `calibrate_v5`, `CalibrationResultV5`, `version="v5"` field, `active=True` field | All unmarked: `calibrate`, `CalibrationResult`. No `version`/`active` fields. |
| Phase numbering | Phases 0–11 with subphases | Single linear ordering of milestones (M1–M9). Each milestone is one PR-sized commit. |
| `L_eff` geometry | `max(0, L − FL + 1)` (containment) | `L + FL − 1` (overlap). Locked from day 1. |
| `partition_units_to_loci` | Buffer-materialization fast path | Anchor-transcript membership (`em_data.locus_t_indices[unit]` → containing `Locus`). No buffer geometry pass. |
| Bootstrap stub (`CalibrationStub`, `bootstrap_fl_calibration`, `compute_default_locus_priors`) | Transitional scaffolding kept for legacy indexes | NOT introduced. We accept that intermediate commits between M1 and M6 will not produce statistically meaningful priors; the protected suite uses synthetic mini-fixtures that don't depend on calibration intelligence. |
| `tests/skip_during_v5_dev.txt` | Tracked file with phase gate enforcement | NOT introduced. Tests are either green at every milestone or marked `pytest.mark.skip(reason=...)` with a TODO referencing the milestone. No registry. |
| Backwards-compat with old indexes | Bootstrap-fallback path in `quant_from_buffer` | Hard error: "rebuild your index". The index format version is bumped once, at M2. |
| `c_base(ℓ)` formula | Open question with 4-candidate spike (Phase 7c) | Constant `c_base = 10.0` from M6 forward. Tunable empirical formula deferred to a future plan once we have real benchmark evidence. |
| `prior_weight_rna` ABI knob (continuous nRNA suppression) | Phase 7a | Kept; shipped at M5. v5.v1 default `nrna_weight = 0.0`. |

---

## 1. Goals

The calibration subsystem has three jobs. In order of priority:

1. **Estimate the gDNA fragment-length (FL) distribution** from a high-purity
   gDNA sub-pool of the input library, so the EM's per-fragment likelihood
   under the gDNA component matches the actual gDNA shape (not the bootstrap
   global FL biased toward RNA).
2. **Estimate the global gDNA contamination level** (one density per
   region type), so per-locus priors are anchored to a library-wide baseline.
3. **Provide per-`MultiLocus` Dirichlet priors** that split the locus's
   prior mass between gDNA and total RNA (= mRNA + nRNA), so the EM has
   informative starting points and does not siphon mRNA into nRNA on
   short or coverage-poor loci.

A non-goal of this plan: separating mRNA from nRNA in calibration. That
split is left entirely to the EM's per-fragment likelihoods. Calibration
delivers a `(gDNA, total RNA)` split per `MultiLocus`; the EM further
distributes total RNA across mRNA and nRNA components.

---

## 2. Locked design

### 2.1 Region partition

Per-reference, walk transcripts in sorted-start order and emit a tiling
of `[0, ref_length)`:

1. **Genic spans** = maximal runs of transcripts that overlap on either
   strand (strand-agnostic). The complement on each reference =
   **intergenic spans**.
2. Within each genic span, merge all-strand exons into `EXON` regions.
   Interstices = `INTRON` regions.
3. `EXON` wins at every position with any-strand exon coverage.
4. Synthetic nRNAs are excluded from this sweep so they cannot
   coalesce real genic spans.
5. Tiny regions are retained at this stage; coverage filtering is a
   downstream concern.

The partition tiles every reference exactly: no gaps, no overlaps.
Naming: **genic / intergenic** (not "island / gap").

### 2.2 Region attributes (on-disk schema)

`regions.feather` schema:

| Column | Dtype | Notes |
|---|---|---|
| `region_id` | int64 | Monotonic across genome; row index ≡ `region_id` |
| `ref_name` | string | Matches `ref_lengths.feather` |
| `start` | int64 | 0-based inclusive |
| `end` | int64 | 0-based exclusive |
| `type` | uint8 | `INTERGENIC=0`, `INTRON=1`, `EXON=2` |
| `strand` | uint8 | `NONE=0`, `POS=1`, `NEG=2`, `AMBIG=3` (strict bitwise OR) |
| `tx_pos_bp` | int64 | bp of (+)-strand transcript overlap |
| `tx_neg_bp` | int64 | bp of (−)-strand transcript overlap |
| `exon_pos_bp` | int64 | bp of (+)-strand exon overlap |
| `exon_neg_bp` | int64 | bp of (−)-strand exon overlap |
| `boundary_flux_left` | bool | left boundary is internal (touches an `INTRON` of the same genic span) |
| `boundary_flux_right` | bool | right boundary is internal |

**Eligibility rules for `boundary_flux_*` (locked):**
- Internal `EXON` boundary that abuts an `INTRON` of the same genic span → `True`.
- Terminal `EXON` boundaries (TSS, TES) → `False`. Rationale: hybrid-capture
  probes typically don't tile across TSS/TES, so unspliced flux at terminal
  boundaries is dominated by capture-edge artefacts. Excluding terminal
  boundaries makes the estimator capture-aware by construction.
- Single-exon transcripts → both flags `False`.
- `INTRON` and `INTERGENIC` regions: both flags always `False`.

**Invariants validated on load:**
- `region_df.index.equals(region_df["region_id"])`.
- `(region_df["end"] - region_df["start"] > 0).all()`.
- For each ref: regions are non-overlapping, sorted, and cover
  `[0, ref_lengths[ref])` with no gaps.
- Every `ref_name` exists in `ref_lengths.feather`.

A native `RegionIndex` (sorted-array, binary-search keyed by integer `ref_id`)
is built at load time. Native code receives `ref_id`; the
`ref_name → ref_id` translation happens once, at `TranscriptIndex.load()`.

### 2.3 Fragment categorization (8-state mask)

3-bit mask `{has_intergenic, has_intron, has_exon}`. Bit ordering in C++ /
numpy is **bit 0 = EXON**, bit 1 = INTRON, bit 2 = INTERGENIC. The Python
`type` column encoding `{INTERGENIC=0, INTRON=1, EXON=2}` is converted
via `type_mask = 1 << (2 - type)`.

| Mask | Name | gDNA FL pool? | Density numerator? |
|---|---|---|---|
| `100` | `INTERGENIC_ONLY` | yes | intergenic |
| `010` | `INTRON_ONLY` | yes | intron |
| `001` | `EXON_ONLY` | no | — |
| `011` | `EXON_INTRON` | yes | exon-intron (boundary flux) |
| `110` | `INTRON_INTERGENIC` | no (annotation gap) | annotation-gap (QC) |
| `101` | `EXON_INTERGENIC` | no (annotation gap) | annotation-gap (QC) |
| `111` | `EXON_INTRON_INTERGENIC` | no (annotation gap) | annotation-gap (QC) |
| `000` | `NO_OVERLAP` | no | impossible if partition complete |

Annotation-gap fragments (masks 5/6/7) are tracked in QC counters but
excluded from the gDNA FL pool and density numerators. They are typically
unannotated splicing or retained introns (i.e., RNA), so including them
biases the gDNA FL toward RNA shape.

### 2.4 Counting unit

A fragment contributes **at most +1** to a given region for category counts,
per-region counts, and boundary-flux counters. Implementation: per-fragment
`SmallRegionSet` (inline-16 + spill vector reused across fragments) for
amortized O(1) dedup with no per-fragment heap traffic.

### 2.5 Boundary-flux counters — gDNA-only

For each `EXON` region we accumulate two `int64_t` counters during the BAM
scan, gated on `splice_type == SPLICE_UNSPLICED AND single-ref fragment`:

- `u_L(R)`: count of unspliced fragments straddling the left boundary.
- `u_R(R)`: count straddling the right boundary.

Tests (with `boundary_tol`, default 0):
```
u_L : frag_start  < region_start − tol   AND  frag_end > region_start
u_R : frag_start  < region_end           AND  frag_end > region_end + tol
```

A fragment may contribute to both `u_L` and `u_R` of the same region only
when it spans the entire region.

Counters are recorded for **every** `EXON` boundary regardless of the
`boundary_flux_*` flag. Eligibility filtering happens at aggregation
time (§2.6.2) so we can swap eligibility policies without rescanning.

There are no spliced (RNA) boundary-flux counters. Under the locked policy
`prior_weight_rna[nrna] = 0`, every non-mRNA fragment in an exon is
attributed to gDNA in the prior; `u_L`/`u_R` is the gDNA estimator.

### 2.6 Three global gDNA densities

Library is summarized by exactly three scalar gDNA densities (fragments / bp):

| Density | Numerator (Σ over regions) | Denominator (Σ over regions) |
|---|---|---|
| ρ̂_INTERGENIC | `per_region_counts[R, 0b100]` over INTERGENIC R | `L_eff(R)` over INTERGENIC R |
| ρ̂_INTRON | `per_region_counts[R, 0b010]` over INTRON R | `L_eff(R)` over INTRON R |
| ρ̂_EXON-INTRON | `u_L(R)·1_L(R) + u_R(R)·1_R(R)` over eligible EXON R | `(1_L(R) + 1_R(R)) · L̄_gDNA` over eligible EXON R |

**Effective length (LOCKED):**
$$L_{\mathrm{eff}}(R) \;=\; L(R) + \bar{L}_{\mathrm{gDNA}} - 1$$

This is the **overlap** geometry — the number of fragment **start positions**
of mean gDNA length that **overlap** region `R` at all (anywhere from a
single bp at the right boundary to a full containment to a single bp at the
left). It must match the resolver's inclusion criterion (`any_overlap`) so
that calibration counts and EM exposure stay consistent. The previous
`max(0, L − L̄ + 1)` (containment) silently zeroed any region shorter than
`L̄_gDNA` (~350 bp on real libraries), causing per-locus gDNA priors to
collapse to zero on short loci. See
`/memories/repo/calibration-overlap-leff-2026-04.md`.

**Why no fourth (EXON) density:** per-`EXON`-region unspliced counts cannot
be attributed to gDNA without ambiguity vs nRNA and pre-mRNA isoforms. EXON
gDNA mass is derived from ρ̂_EXON-INTRON via §2.7.

#### 2.6.1 Two-level EB shrinkage

```
global density (per type)
        ↓ NB-shrink with κ_t
locoregional density (per Locus, per type)
```

Closed-form NB shrinkage:
$$\hat\rho_{\mathrm{loco}}(L,t) \;=\; \frac{N(L,t) + \kappa_t \cdot \hat\rho_{\mathrm{global}}(t)}{L_{\mathrm{eff}}(L,t) + \kappa_t}$$

There is no per-reference (per-chromosome) level. The locoregional level is
computed on demand at locus-prior assembly time (§2.7); never pre-tabulated.

#### 2.6.2 Counters → density conversion (EXON-INTRON eligibility)

```python
for R in exon_regions:
    n_L = u_L[R] if R.boundary_flux_left  else 0
    n_R = u_R[R] if R.boundary_flux_right else 0
    sides = int(R.boundary_flux_left) + int(R.boundary_flux_right)
    if sides == 0:
        continue
    numerator   += n_L + n_R
    denominator += sides * L_gdna_mean
```

#### 2.6.3 NB overdispersion κ — per-region MoM

Naïve MoM that collapses all regions into one mean exposure conflates
true overdispersion with exposure heterogeneity. We use a per-region NB
Method-of-Moments estimator:

$$\hat\kappa = \frac{\sum_R \mu_R^2}{\sum_R (N_R - \mu_R)^2 - \sum_R \mu_R}$$

where $\mu_R = \hat\rho_{\mathrm{global}} \cdot L_{\mathrm{eff}}(R)$.

```python
KAPPA_DEFAULT       = 100.0
KAPPA_MIN           = 1.0
KAPPA_MAX           = 1.0e6
MIN_REGIONS_FOR_MOM = 5
```

Fallback to `KAPPA_DEFAULT` when:
- `n_regions < MIN_REGIONS_FOR_MOM`, or
- `rho_hat == 0`, or all `eff_lengths == 0`, or
- `excess_variance ≤ 0` (Poisson or tighter — underdispersion).

Result is clipped to `[KAPPA_MIN, KAPPA_MAX]`.

### 2.7 Per-`Locus` gDNA mass

A `MultiLocus` is a connected component of transcripts (the EM-execution
unit). It contains one or more `Locus` intervals (contiguous genomic
spans; the calibration-estimation unit). Most `MultiLocus`es have exactly
one `Locus`; paralog clusters across refs have several.

#### 2.7.1 Algorithm (per `Locus`)

For each `Locus = (ref, start, end)`:

1. **Overlap query:** `region_ids = region_index.overlap(ref_id, start, end)`.
   If empty → raise `RuntimeError` (BAM reference doesn't match index).

2. **Per-type contributions** (clipping each region to the locus):
   - `INTERGENIC` and `INTRON` follow §2.6 numerator/denominator restricted
     to overlapping regions of that type.
   - `EXON-INTRON` aggregates `u_L`/`u_R` over flag-eligible `EXON` regions
     in the locus.
   - EB-shrink each to the global density of that type via §2.6.1.

3. **Per-type predicted gDNA count:**
   $$\hat{N}^{\mathrm{gDNA}}_t(L) = \hat\rho_{\mathrm{loco}}(L,t) \cdot L_{\mathrm{eff}}(L,t)$$

   For `t = EXON-INTRON`, the multiplier is the **full** exonic
   $L_{\mathrm{eff}}$ inside `L` (not just flag-eligible boundaries):
   $$L_{\mathrm{eff}}^{\mathrm{EXON-INTRON}}(L) = \sum_{R \in \mathrm{EXON}(L)} \left(L(R) + \bar L_{\mathrm{gDNA}} - 1\right)$$

   The density already encodes capture-aware bias (eligibility filtering at
   the density level); the predicted count is the gDNA mass we expect on
   all exonic territory inside `L`.

4. **Total locus gDNA mass and fraction:**
   $$\hat N^{\mathrm{gDNA}}(L) = \sum_t \hat N^{\mathrm{gDNA}}_t(L)$$
   $$\hat\pi^{\mathrm{gDNA}}(L) = \mathrm{clip}\big(\hat N^{\mathrm{gDNA}}(L) / N_{\mathrm{obs}}(L),\ 0,\ 1\big)$$

   `N_obs(L)` is the count of EM units assigned to `L` (see §2.8).

5. **Sanity:** π̂ > 1 before clipping → `WARN` (diagnostic of capture leakage
   or EB pathology, but not blocking). Non-finite output → hard error.

#### 2.7.2 `MultiLocus` aggregation

For the EM, each `MultiLocus` receives a single combined gDNA prior:

```python
def assemble_multilocus_prior(ml, per_locus_estimates):
    n_obs  = sum(e.n_obs  for e in per_locus_estimates)
    n_gdna = sum(e.n_gdna for e in per_locus_estimates)
    n_rna  = max(0.0, n_obs - n_gdna)
    return MultiLocusPrior(
        n_obs=n_obs, n_gdna=n_gdna, n_rna=n_rna,
        per_locus=per_locus_estimates,   # retained for diagnostics
    )
```

Per-`Locus` estimates are retained on `MultiLocusPrior.per_locus` so a
future per-`Locus`-component EM extension can use them without
re-estimating.

### 2.8 `N_obs(L)` — anchor-transcript routing

`N_obs(L)` is the count of EM units routed to `Locus L`. Routing uses
**anchor-transcript membership**, not per-fragment genomic coordinates:

```python
def partition_units_to_loci(ml, em_data, t_to_locus_idx):
    """Route EM units in MultiLocus ml to its Locus intervals.

    Each unit carries an anchor transcript em_data.locus_t_indices[unit].
    Each transcript belongs to exactly one Locus in ml.loci (build_loci
    invariant). The mapping unit → anchor → Locus is therefore a pure
    transcript-table lookup with no buffer access.
    """
```

Fast path (`len(ml.loci) == 1`, ≈99% of `MultiLoci`): return
`(ml.unit_indices,)` directly.

Slow path: pre-compute `t_idx → local_locus_idx` for all transcripts in
`ml.transcript_indices` by mask-and-assign over `ml.loci` (typically ≤ 3
loci); then `np.searchsorted` against
`em_data.locus_t_indices[ml.unit_indices]` to bin units.

Membership is total by construction. A transcript that does not lie in any
`Locus` is a `build_loci` invariant violation → `RuntimeError`.

See `/memories/repo/multilocus-partition-design-2026-04.md` for the
anti-pattern (buffer-materialization) we are explicitly **not** doing.

### 2.9 nRNA suppression — continuous prior weight in C++

Native EM solver receives a per-component float weight:
`prior_weight_rna: std::vector<float>` (one entry per EM component, in
`[0, 1]`).

Fan-out:
```cpp
double total_weighted_rna_coverage = 0.0;
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] > 0.0 && i != gdna_idx) {
        total_weighted_rna_coverage += prior_weight_rna[i] * coverage_totals[i];
    }
}
// later, per-component prior:
prior_out[i] = baseline + std::max(
    alpha_rna * prior_weight_rna[i] * coverage_totals[i]
        / total_weighted_rna_coverage,
    EM_LOG_EPSILON);
```

**Default policy:** `prior_weight_rna[mrna] = 1.0`, `prior_weight_rna[nrna] = 0.0`.

`nullptr` (Python `None`) is a backwards-compatible all-ones default.

### 2.10 Calibration observation policy

Exactly one policy decides which fragments enter the calibration accumulator,
**at the scanner observation site** (not downstream). The payload is
"usable-only plus exclusion counters":

| Fragment class | Policy | Counter |
|---|---|---|
| Unique mapper, non-chimeric, non-artifact, resolved | observed | `n_observed` |
| Unique mapper, intergenic-only | observed | (same) |
| Multi-mapper (`num_hits > 1`), any class | excluded | `n_excluded_multimap` |
| Chimeric (`chimera_type != CHIMERA_NONE`) | excluded | `n_excluded_chimera` |
| Artifact (`splice_type == SPLICE_ARTIFACT`) | excluded | `n_excluded_artifact` |
| Out-of-region (mask `0b000` after observe) | observed (in `global_counts[0]`); QC flagged | `n_oor` |

**Multimapper rationale:** including them would require either per-hit `+1`
(over-counts by `NH`) or `1/NH` fractional increments (requires float
counters and complicates κ estimation). Excluding them systematically
under-estimates gDNA density on multi-mappable regions; this is documented
and accepted. A future mappability-adjusted $L_{\mathrm{eff}}^{\mathrm{map}}(R) = \int_R m(x) dx$
can recover unbiased density without re-introducing multimappers.

**Chimera/artifact rationale:** chimeric fragments have ambiguous
footprints; artifact junctions are pre-classified as gDNA-derived junk by
the splice blacklist. Including either corrupts every downstream density.

### 2.11 Reference identity contract

- On disk: all DataFrames key references by `ref_name` (string).
- In native arrays / cgranges payloads: integer `ref_id` derived from the
  canonical `ref_lengths.feather` order.
- Translation point: exactly **one** place — `TranscriptIndex.load()` builds
  `ref_name_to_id: dict[str, int]` and `ref_names: list[str]` from
  `ref_lengths.feather`. All callers consume one or the other.
- Calibration overlap engine: the native `RegionIndex` is the **single
  source** for region-overlap queries on the calibration path. Its
  `set_regions` rejects unknown `ref_id` rather than auto-creating refs.

### 2.12 Diagnostic surfaces explicitly removed

The following per-region quantities are **not** part of the calibration
result schema:
- `s_L`, `s_R` per-`EXON` spliced-block counters.
- ρ̂_RNA(R) per-region RNA boundary density.
- π̂_gdna(R) per-region boundary-flux ratio.
- Per-region (per-`EXON`) gDNA density of any kind. EXON gDNA mass is
  derived (§2.7), never published at region granularity.
- Per-reference (per-chromosome) EB level.

`density_by_type` carries exactly three entries: `INTERGENIC`, `INTRON`,
`EXON-INTRON`.

### 2.13 Library / style decisions

- pandas + pyarrow only. **No Polars.**
- C++17, `-O3`, LTO. Hot path is allocation-free.
- All C++ public types live under `namespace rigel::calibration`.

---

## 3. Implementation discipline (LESSONS FROM DATA LOSS)

The 2026-05-01 incident destroyed ~14 days of uncommitted Phase 6–9 work.
Root causes that this plan must prevent:

1. **A bulk regex `re.sub(r'  +(?=\\S)', ' ', source)` collapsed Python
   indentation across 192 files in seconds.** Never apply whitespace tidies
   to source code without operating line-by-line with leading-whitespace
   preservation. Comment-only edits should be filtered to comment lines.
2. **VS Code Local History only snapshots files at editor save events;
   in-conversation agent edits to brand-new files were not captured.** This
   means untracked files have **no recoverable backup** if an automated
   edit corrupts them.
3. **Multi-day untracked work piled up because there was no commit
   discipline tied to milestone boundaries.**

**Mandatory rules for this rebuild:**

- **Commit at every milestone boundary**, even if subsequent work would
  squash. Each milestone in §4 below is a single PR-sized commit.
- **Commit all new files immediately** after they pass their own unit tests
  (do not let untracked `.py`/`.cpp` files accumulate).
- **Never run a bulk text transformation on `src/` or `tests/` without
  first committing.** If a transformation is wrong, the recovery path is
  `git checkout -- .` — which only works if the prior state is committed.
- **Whitespace and `Phase X`-style comment scrubs are not part of any
  milestone.** Leave new code clean from the start; if cleanup is needed,
  it is a dedicated commit done via deterministic `replace_string_in_file`
  edits, not regex passes.

---

## 4. Milestones

Each milestone produces one commit with a focused diff: ~one new module,
one new test file, and any wiring changes required to keep the protected
suite green. Total: 9 milestones.

### M1 — Region partition (Python)

**Scope:** add the region table to the index build; do not yet load it
back, do not yet wire it into the scanner.

**Code:**
- `src/rigel/calibration/regions.py` (new):
  - `RegionType(IntEnum)`: `INTERGENIC=0`, `INTRON=1`, `EXON=2`.
  - `RegionStrand(IntFlag)`: `NONE=0`, `POS=1`, `NEG=2`, `AMBIG=3`.
  - `emit_regions(ref, layout) -> Iterator[RegionRecord]` consumer of the
    layout iterator.  Per-genic-span partition uses an **endpoint event
    sweep** (O(E log E) time, O(E) memory in the number of intervals;
    NEVER allocates per-base arrays) so that mega-loci such as HLA — where
    a single genic span can exceed a few Mb — never trigger a transient
    multi-hundred-MB allocation.
- `src/rigel/index.py`:
  - `_GenicSpan`, `_IntergenicSpan` frozen dataclasses.
  - `_iter_reference_layout(ref_length, transcripts)` — single sweep
    yielding interleaved spans.
  - Refactor existing `_gen_genomic_intervals` to consume the same iterator.
  - `build_index_artifacts(transcripts, ref_lengths) -> (intervals_df, regions_df)`.
  - `TranscriptIndex.build()` writes both `intervals.feather` and
    `regions.feather`.

**Tests:** `tests/test_layout_iter.py`, `tests/test_region_partition.py`.

**Exit gate:** `intervals.feather` byte-identical to pre-refactor baseline;
new tests green; protected suite green.

### M2 — Region index persistence + load validation

**Scope:** load and validate `regions.feather` on `TranscriptIndex.load()`;
attach `region_df`, `ref_lengths`, and `ref_name_to_id` to the index;
bump `INDEX_FORMAT_VERSION`.

The load-time invariant — per-reference, regions are sorted, non-
overlapping, and tile `[0, ref_lengths[ref])` exactly — is what lets M3
use a **native per-reference CSR** (`region_index.h`, `upper_bound` plus
linear scan) for overlap queries.  The calibration path **does not**
introduce a second cgranges instance; reusing the existing resolver-side
cgranges would cost a per-fragment dynamic-allocation hit that the
per-ref CSR avoids.

**Code:**
- `src/rigel/calibration/regions.py`:
  - `load_regions(path) -> pd.DataFrame`.
  - `validate_against_ref_lengths(region_df, ref_lengths)`.
- `src/rigel/index.py`:
  - Bump `INDEX_FORMAT_VERSION` (current = 2 → 3).
  - Load `regions.feather`, validate invariants, attach to `TranscriptIndex`.
  - `TranscriptIndex.region_df`, `TranscriptIndex.ref_lengths` slots.
  - **Hard error** on stale index — no fallback path.

**Tests:** `tests/test_region_persist.py` (load happy-path + 3 invariant
violations).

**Exit gate:** new tests green; protected suite green; `tests/test_index.py`
fixture extended with `ref_lengths.feather` + minimal `regions.feather`.

### M3 — C++ scanner: in-place calibration accumulator

**Scope:** scanner builds a per-reference CSR `RegionIndex` from
`region_df`; per-worker `CalibrationAccumulator` tracks the 8-state mask,
per-region counts, per-mask FL histogram, and `u_left`/`u_right`
boundary-flux counters.  Worker results merged after `workers.join()`.
Payload exported via nanobind.

**Code:**
- `src/rigel/native/calibration/region_index.h` — per-reference CSR over
  sorted region intervals, `upper_bound` + linear scan overlap.  Public
  accessors: `set(ref_ids, starts, ends, type_masks, n_regions, n_refs)`,
  `overlap_into<N>(ref_id, qstart, qend, out_inline, out_spill,
  already_inline)`, `type_mask`, `start`, `end`, `n_regions`, `n_refs`.
- `src/rigel/native/calibration/accumulator.h` and `.cpp` —
  `CalibrationAccumulator(int32_t n_regions)` with `observe(splice_type,
  ref_id, frag_start, frag_end, exons, n_exons, region_index)` and
  `merge_from(other)`.
- `src/rigel/native/bam_scanner.cpp`:
  - `BamScanner::set_regions(ref_ids, starts, ends, type_masks, n_refs)`
    — Python-callable, takes `int32_t` ref-ids (NOT strings; the
    Python-side `_wire_calibration_regions` translates ref-name → ref-id
    via `index.resolver.get_ref_to_id()` before calling).  Must be called
    before `scan()`; raises on length mismatch or double-set.
  - `WorkerState::cal_acc(n_regions)` per worker.
  - One observation site at the end of `process_qname_group_threaded`,
    plus an early dispatch for `is_multimap` that calls
    `cal_acc.note_multimap()` and returns.  Chimeric and artifact
    excluded classes use `cal_acc.note_chimera()` /
    `cal_acc.note_artifact()`.  Resolved fragments with no observable
    hit call `cal_acc.note_unobserved()` so the balance assertion stays
    exact.
  - After `workers.join()`, fold worker accumulators into the merged
    payload.
  - `build_result()` exports `result["calibration"]` dict (see below).
- `src/rigel/calibration/scan_payload.py`:
  - `@dataclass(frozen=True) class CalibrationScanPayload` with fields:
    `global_counts (8,)`, `per_region_counts (n_regions, 8)`,
    `fl_hist (8, 1024)`, `u_left (n_regions,)`, `u_right (n_regions,)`,
    `n_observed`, `n_excluded_multimap`, `n_excluded_chimera`,
    `n_excluded_artifact`, `n_unobserved`, `n_oor`.
  - `from_scan_dict(d, *, n_total=None)` validator + balance assertion
    `n_observed + n_excluded_multimap + n_excluded_chimera +
    n_excluded_artifact + n_unobserved == n_total` where the balance
    basis `n_total` is `stats.n_read_names` (not `stats.n_fragments`).

**Pass D (boundary-flux) inside `observe()`** (§2.5) — fires for
`splice_type == SPLICE_UNSPLICED` and only against exon-typed regions
that the fragment overlaps; **no `boundary_tol`** parameter (the closed-
box `frag_start < rs && frag_end > rs` test is exact for half-open
regions and the integer fragment endpoints htslib delivers):

```cpp
if (splice_type == SPLICE_UNSPLICED) {
    auto check_one = [&](int32_t rid) {
        if ((regions.type_mask(rid) & kExonBit) == 0) return;
        const int64_t rs = regions.start(rid), re = regions.end(rid);
        if (frag_start < rs && frag_end > rs) u_left[rid]++;
        if (frag_start < re && frag_end > re) u_right[rid]++;
    };
    /* iterate the inline+spill region set */
}
```

**Tests (`tests/test_calibration_accumulator.py`):**
- `TestSetRegions` (3): binding contract, length-mismatch rejection,
  double-set rejection.
- `TestPipelinePayload` (2): scan returns payload, `run_pipeline`
  attaches `PipelineResult.calibration_payload`.
- `TestPayloadValidation` (7): shape, dtype, balance, `n_oor`
  accounting, `None`-dict rejection.
- `TestWorkerMergeEquality` (1): 1-vs-4-worker payload byte-equal.
- `TestMaskCorrectness` (5): hand-built BAM, one unspliced fragment per
  intended mask (`EXON_ONLY`, `INTRON_ONLY`, `INTERGENIC_ONLY`,
  `EXON+INTRON`, `INTERGENIC+EXON`); assert exactly one bin populated.
- `TestBoundaryFlux` (4): single-fragment left-straddle, right-straddle,
  full-span (both flags), interior (neither flag) against a known exon
  region in `mini_index`.
- `TestObservationPolicy` (2 + 1 deferred): NH=2 → `n_excluded_multimap`,
  trans-chromosomal pair → `n_excluded_chimera`.  `SPLICE_ARTIFACT` is
  deferred to M9 once `splice_blacklist.feather` is wired into the test
  fixtures (skip with explicit reason).

**Exit gate:** new tests green; protected suite green;
`PipelineResult.calibration_payload: CalibrationScanPayload | None`
populated end-to-end.

### M4 — Global gDNA densities + κ

**Scope:** compute the three densities from §2.6 and one κ estimate per type.

**Code:**
- `src/rigel/calibration/density_global.py`:
  - `@dataclass(frozen=True) class GlobalGdnaDensity` with
    `type, rho, n_fragments, eff_length_bp, n_regions_used, kappa`.
  - `@dataclass(frozen=True) class GlobalDensityTable` with
    `intergenic, intron, exon_intron, gdna_fl_mean`.
  - `compute_global_densities(region_df, payload, gdna_fl_mean) -> GlobalDensityTable`
    — pure pass over the region table + payload.
- `src/rigel/calibration/_kappa.py`:
  - `KappaEstimate` dataclass + `estimate_kappa(counts, eff_lengths, rho_hat)`
    per §2.6.3.

**Tests:** `tests/test_density_global.py` (hand-counted 3-region scenario;
EXON-INTRON with all flags False → density 0; INTRON-only sample;
pure-mRNA → all densities < 1e-3; 50/50 sample → three densities within 30%).
`tests/test_estimate_kappa.py` (degenerate cases + heterogeneous-exposure
consistency at $\kappa^* \in \{2, 20, 200\}$ recovering κ within ±30%).

**Exit gate:** ≥ 11 new tests green; protected suite green.

### M5 — EM `prior_weight_rna` ABI + `Locus`/`MultiLocus` rename

**Scope:** native EM accepts continuous nRNA suppression; locus
nomenclature codified.

**Code (a) — EM ABI:**
- `src/rigel/native/em_solver.cpp`:
  - `compute_ovr_prior_and_warm_start` gains `const float* prior_weight_rna`
    (nullptr ≡ all-ones for backwards-compat).
  - Fan-out per §2.9.
- `src/rigel/native/_em_impl.cpp`:
  - `run_em_locus(..., prior_weight_rna: np.ndarray | None = None)`.
  - `run_em_batch(..., locus_prior_weight_rna: list[np.ndarray] | None = None)`.

**Code (b) — Rename:**
- `src/rigel/locus.py`:
  - Rename existing `class Locus` → `class MultiLocus`.
  - New `@dataclass(frozen=True) class Locus` with fields `(ref, ref_id, start, end)`.
  - `MultiLocus.loci: tuple[Locus, ...]` built from existing
    `merged_intervals` at construction.
  - Drop `merged_intervals` entirely (no deprecation alias — strict cut).
  - `build_prior_weight_rna(multi_locus, em_data, nrna_weight=0.0) -> np.ndarray`
    helper.
- All call sites that consumed `Locus.merged_intervals` migrate to
  `multi_locus.loci`.

**Tests:**
- `tests/test_em_prior_weight.py` (5 cases: bit-identical with `None` /
  full suppression / half-share / VBEM baseline preserved / degenerate
  all-zero RNA weights).
- `tests/test_locus_rename.py` (single-Locus MultiLocus, multi-ref
  MultiLocus, `Σ (l.end - l.start)` parity).

**Exit gate:** ≥ 8 new tests green; EM golden outputs bit-identical when
caller passes `None`; protected suite green.

### M6 — Per-`Locus` gDNA mass + `MultiLocus` priors

**Scope:** the numerical core of calibration. §2.7 verbatim.

**Code:**
- `src/rigel/calibration/density_loco.py`:
  - `shrink_to_loco(n_loco, leff_loco, rho_global, kappa) -> float` per §2.6.1.
- `src/rigel/calibration/_arrays.py`:
  - `RegionArrays.from_region_df(region_df)` — pre-flatten the 5 columns
    we touch (starts, ends, types, flag_left, flag_right) into numpy arrays.
  - `PayloadArrays.from_payload(payload)` — pre-extract the 4 reductions
    (intergenic per region, intron per region, `u_L`, `u_R`).
- `src/rigel/calibration/_locus_n_obs.py`:
  - `partition_units_to_loci(ml, em_data, t_to_locus_idx)` per §2.8.
- `src/rigel/calibration/locus_prior.py`:
  - `LocusGdnaEstimate` (frozen, slots) — see schema below.
  - `MultiLocusPrior` (frozen, slots) — see schema below.
  - `PriorTable` (frozen) carrying `multi_locus_priors`, `alpha_gdna`,
    `alpha_rna`, `prior_weight_rna`.
  - `estimate_locus_gdna(...)` per §2.7.1.
  - `assemble_multilocus_prior(...)` per §2.7.2.
  - `assemble_priors(multi_loci, em_data, index, payload, global_densities,
    gdna_fl_mean, c_base) -> PriorTable` orchestrator.
  - `C_BASE_DEFAULT = 10.0`.
  - Fallback bitmask constants:
    `FLAG_INTERGENIC_ZERO_LEFF`, `FLAG_INTRON_ZERO_LEFF`,
    `FLAG_EXON_INTRON_NO_ELIGIBLE`, `FLAG_PI_CLIPPED`.

**Schemas:**
```python
@dataclass(frozen=True, slots=True)
class LocusGdnaEstimate:
    locus:                Locus
    n_obs:                int
    n_gdna_intergenic:    float
    n_gdna_intron:        float
    n_gdna_exon_intron:   float
    n_gdna:               float
    pi_gdna:              float        # clipped [0, 1]
    rho_loco:             tuple[float, float, float]
    leff_loco:            tuple[float, float, float]
    n_eligible_boundaries: int
    fallback_flags:       int

@dataclass(frozen=True, slots=True)
class MultiLocusPrior:
    multi_locus_id:  int
    n_obs:           int
    n_gdna:          float
    n_rna:           float
    pi_gdna:         float
    per_locus:       tuple[LocusGdnaEstimate, ...]
```

**Critical implementation details:**
- `L_eff(R)` uses **overlap geometry**: `cl + L̄_gDNA - 1` where `cl` is the
  region length clipped to the locus interval. **Never** use
  `max(0, cl - L̄_gDNA + 1)`.
- `partition_units_to_loci` slow path uses **anchor-transcript membership**
  via `em_data.locus_t_indices[unit]`. **Never** re-introduce a
  `buffer.materialize_frag_id_indexed_geometry()` shortcut.

**Tests:**
- `tests/test_density_loco.py` (≥ 6 scalar cases).
- `tests/test_partition_units.py` (≥ 5 cases including transcript-anchor
  routing on a 2-`Locus` paralog `MultiLocus`).
- `tests/test_per_locus_gdna_mass.py` (≥ 8 cases pinning §2.7.1 end-to-end
  with hand-built fakes — including the **short-locus overlap-L_eff** case
  that the previous containment formula silently zeroed).
- `tests/test_assemble_priors.py` (≥ 5 cases for `MultiLocus` aggregation
  + orchestration).

**Exit gate:** ≥ 24 new tests green; protected suite green; performance
budget < 1s for 5K MultiLoci on a synthetic stress test.

### M7 — Pool FL models + `CalibrationResult` schema

**Scope:** build the gDNA pool FL model from `payload.fl_hist`; carry the
scan-trained RNA/global FL models through unchanged; assemble the
canonical `CalibrationResult` carrier.

**Critical FL semantics — read this before touching M7.**  The M3
accumulator records `fl_raw = frag_end - frag_start` (genomic span) into
`fl_hist[mask, :]` for every observed fragment.  For unspliced fragments
this equals the true fragment length; for **spliced** fragments
(`SPLICE_SPLICED`, which dominates `mask = EXON_ONLY = 0b001`) this is
the genomic span, NOT the fragment length, and is on average several
kilobases too large.  Therefore:

- **`payload.fl_hist[gdna_pool_masks, :]` IS a valid fragment-length
  histogram** — gDNA fragments are never spliced, so their genomic span
  equals their fragment length by construction.
- **`payload.fl_hist[1, :]` (`EXON_ONLY`) is NOT a fragment-length
  histogram** — it is dominated by spliced fragments whose genomic span
  bears no relationship to their FL.
- M7 therefore uses the payload **only** for the gDNA pool FL.  The
  RNA-pool and global FL models continue to come from the scan-trained
  `frag_length_models` carried out of `scan_and_buffer` (which work in
  transcript space via `_resolve_core` and so are correct for spliced
  fragments by construction).  This is what `pipeline.py` already does;
  M7 codifies it.

  *Future option (out of scope here):* widen `fl_hist` to a `(8, 2,
  kFlBins)` tensor that splits each mask into spliced vs unspliced strata
  using the resolved `splice_type`.  That would let an `EXON_ONLY +
  unspliced` slice supply a payload-derived RNA FL.  Tracked as a
  follow-up; not required for the calibration to be correct.

**Code:**
- `src/rigel/calibration/_fl_pool.py`:
  - `PoolFLModels` (frozen, slots): `gdna_fl_model, rna_fl_model,
    global_fl_model, gdna_fl_quality, n_pool, n_rna, n_global,
    n_pool_annotation_gap`.
  - `compute_pool_fl_models(
      fl_hist, *, scan_rna_fl_model, scan_global_fl_model,
      max_size, prior_ess=1000.0,
      quality_threshold_good=5_000, quality_threshold_weak=200,
    ) -> PoolFLModels`.
  - Pool definition:
    - gDNA pool = mask 2 ∪ mask 3 ∪ mask 4 (`INTRON_ONLY`, `EXON_INTRON`,
      `INTERGENIC_ONLY`) — all unspliced by construction, so genomic
      span equals FL.
    - RNA pool model = `scan_rna_fl_model` (passed through; no payload
      use).
    - Global model = `scan_global_fl_model` (passed through; no payload
      use).
  - Quality classifier (gates only the gDNA model):
    - `n_pool ≥ 5000` → `"good"`, no shrinkage.
    - `200 ≤ n_pool < 5000` → `"weak"`, EB-shrink toward global with
      `prior_ess=1000`.
    - `n_pool < 200` → `"fallback"`, gDNA FL ≡ global.
- `src/rigel/calibration/_result.py`:
  - `CalibrationResult` (frozen, slots) — schema below.
  - `build_calibration_result(pool_models, global_densities, payload,
    region_df, strand_specificity) -> CalibrationResult`.
  - `_build_multi_locus_prior_df`, `_build_per_locus_gdna_df` helpers.

**`CalibrationResult` schema:**
```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    # FL models (consumed by FragmentScorer)
    gdna_fl_model:   FragmentLengthModel
    rna_fl_model:    FragmentLengthModel
    global_fl_model: FragmentLengthModel
    gdna_fl_quality: str            # "good" | "weak" | "fallback"

    # Library-level diagnostic
    strand_specificity: float

    # Region partition + global densities
    region_df:             pd.DataFrame
    global_gdna_densities: GlobalDensityTable

    # Per-MultiLocus priors (lazily backfilled by `with_priors`)
    multi_locus_prior_df: pd.DataFrame | None = None
    per_locus_gdna_df:    pd.DataFrame | None = None
    c_base_value:         float = C_BASE_DEFAULT

    # Pool diagnostics
    pi_pool:                float = 0.0
    n_pool:                 int = 0
    n_rna:                  int = 0
    n_global:               int = 0
    n_pool_annotation_gap:  dict[str, int] = field(default_factory=dict)

    # Scan-time exclusion counters
    n_multimap_excluded:    int = 0
    n_chimera_excluded:     int = 0
    n_artifact_excluded:    int = 0
    n_oor_excluded:         int = 0

    # 8-element mask vector + boundary-flux summary + κ diagnostics
    category_counts:        np.ndarray
    boundary_flux_gdna_summary: dict
    kappa_diagnostics:      dict[str, KappaEstimate]

    def to_summary_dict(self) -> dict[str, Any]: ...
    def with_priors(self, multi_locus_prior_df, per_locus_gdna_df) -> "CalibrationResult":
        return dataclasses.replace(
            self,
            multi_locus_prior_df=multi_locus_prior_df,
            per_locus_gdna_df=per_locus_gdna_df,
        )
```

**No `version` field. No `active` field.** This is THE calibration result.

**Tests:** `tests/test_pool_fl.py` (≥ 8), `tests/test_calibration_result.py`
(≥ 10 including frozen contract via `with_priors` and dataframe schema
locks).

**Exit gate:** ≥ 18 new tests green; protected suite green.

### M8 — `calibrate()` orchestrator + pipeline integration

**Scope:** compose M4 + M7 into a top-level `calibrate()` and migrate
the pipeline off the live SRD-v1 surface.  This is **NOT** a hard cut
from a stub; the migration target is the in-tree calibration code
(`_simple.py`, `_categorize.py`, `_fl_mixture.py`,
`_fl_empirical_bayes.py`, `_result.py`, the SRD-v1 `CalibrationConfig`
fields, the SRD-v1 CLI flags, and the SRD-v1 `summary.json` keys).
Splitting into three focused commits keeps each diff reviewable and lets
us bisect any benchmark regression to a single switchover.

**Migration table (SRD-v1 → v6):**

| Surface | SRD-v1 (current) | v6 (target) |
|---|---|---|
| Module entrypoint | `from .calibration import calibrate_gdna` | `from .calibration import calibrate` |
| Result type | `CalibrationResult` (in `_result.py`, v1 schema) | `CalibrationResult` (in `_result.py`, v6 schema below) |
| Config dataclass | `CalibrationConfig.exon_fit_tolerance_bp`, `fl_prior_ess`, `max_iter`, `tol` | `CalibrationConfig.pool_quality_thresholds=(5000,200)`, `prior_ess=1000.0`, `nrna_weight`, `c_base` |
| CLI flags | `--cal-exon-fit-tolerance-bp`, `--cal-fl-prior-ess`, `--cal-max-iter`, `--cal-tol` | `--cal-quality-good`, `--cal-quality-weak`, `--cal-prior-ess`, `--cal-nrna-weight`, `--cal-c-base` |
| `summary.json` keys | `calibration.pi_pool`, `gdna_fl.{mu,sigma,quality}`, `srd.*` | `calibration.pi_pool`, `calibration.gdna_fl.{mu,sigma,quality}`, `calibration.kappa_diagnostics`, `calibration.boundary_flux_gdna_summary`, `calibration.global_densities` |
| `quant_from_buffer` arg | `calibration: CalibrationResult` (v1) | `calibration: CalibrationResult` (v6 — same name, new schema) |
| Per-locus prior | `compute_locus_priors_from_partitions(...)` | `assemble_priors(...)` (M6) → `calibration.with_priors(...)` |

**Three-commit split:**

#### M8a — Introduce v6 surface alongside SRD-v1 (additive)

- New `src/rigel/calibration/_orchestrator.py` exporting `calibrate(...)`.
- New `src/rigel/calibration/_fl_pool.py` (per M7).
- New `src/rigel/calibration/_result_v6.py` exporting
  `CalibrationResultV6` (renamed back to `CalibrationResult` in M8c when
  the v1 class is deleted; the temporary suffix avoids collision).
- New `CalibrationConfig.v6_*` fields added alongside the v1 fields; v1
  fields untouched.
- New CLI flags added (`--cal-quality-good`, etc.); v1 flags untouched.
- No change to `pipeline.run_pipeline` or `quant_from_buffer` yet.
- Tests: `tests/test_calibrate_orchestrator.py` (6 cases) verifying the
  new surface produces a sensible `CalibrationResultV6` from a synthetic
  payload + index.
- Exit gate: new tests green; SRD-v1 path still wired and still passes
  the protected suite untouched.

#### M8b — Pipeline switchover (replace v1 wiring)

- `pipeline.run_pipeline` calls `calibrate(...)` instead of
  `calibrate_gdna(...)`.
- `pipeline.quant_from_buffer` consumes the v6 result shape;
  `assemble_priors(...)` (M6) backfills via `with_priors(...)`.
- `_run_locus_em_partitioned` accepts `prior_weight_rna_per_locus`.
- CLI: v1 flags become deprecation-warning shims that map onto the
  closest v6 flag (or are ignored with a one-line warning if no
  equivalent exists).
- `summary.json` writer emits the v6 keys; v1 keys are temporarily
  emitted in parallel to avoid breaking external dashboards.
- Golden outputs regenerated; diffs documented in the commit message.
- Tests: `tests/test_pipeline_integration.py` (8 cases) — end-to-end on
  synthetic scenarios + `π̂_pool > 0` on a contaminated scenario.
- Exit gate: full suite green; benchmark sweep (Armis2 13-condition
  matrix, see `.github/copilot-instructions.md`) shows no regression vs
  the M8a baseline.

#### M8c — Legacy deletion (subtractive)

- Delete `_simple.py`, `_categorize.py`, `_fl_mixture.py`,
  `_fl_empirical_bayes.py`, the v1 `_result.py`.  Move the two general
  FL utilities to `src/rigel/frag_length_mixture.py` and
  `src/rigel/frag_length_eb.py` (they are not calibration-specific).
- Rename `CalibrationResultV6` → `CalibrationResult`; rename
  `_result_v6.py` → `_result.py`.
- Drop the v1 `CalibrationConfig` fields and v1 CLI flags (after one
  release cycle of the M8b deprecation warnings).
- Drop the v1 keys from `summary.json`.
- Delete the v1 tests: `tests/test_calibration_simple.py`,
  `tests/test_categorize.py`, `tests/test_gdna.py`,
  `tests/test_gdna_harmonic_length.py`.
- Exit gate: full suite green; legacy modules absent (verified via
  `git grep`); CHANGELOG entry written.

**Hard cut** — no bootstrap fallback, no `CalibrationStub`.

**Code:**
- `src/rigel/calibration/_orchestrator.py`:
  ```python
  def calibrate(
      buffer:             FragmentBuffer,         # reserved
      index:              TranscriptIndex,
      payload:            CalibrationScanPayload,
      frag_length_models: FragmentLengthModels,
      strand_models:      StrandModels,
      *,
      pool_quality_thresholds: tuple[int, int] = (5_000, 200),
  ) -> CalibrationResult:
      if index.region_df is None:
          raise RuntimeError(
              "Index has no region table. Rebuild the index "
              "(rigel index --fasta ... --gtf ...). Older indexes "
              "are not supported."
          )
      pool = compute_pool_fl_models(
          payload.fl_hist,
          max_size=frag_length_models.max_size,
          quality_threshold_good=pool_quality_thresholds[0],
          quality_threshold_weak=pool_quality_thresholds[1],
      )
      global_dens = compute_global_densities(
          index.region_df, payload, gdna_fl_mean=pool.gdna_fl_model.mean,
      )
      return build_calibration_result(
          pool_models=pool,
          global_densities=global_dens,
          payload=payload,
          region_df=index.region_df,
          strand_specificity=strand_models.strand_specificity,
      )
  ```
- `src/rigel/calibration/__init__.py` exports the public surface listed below.
- `src/rigel/pipeline.py`:
  - `scan_and_buffer` returns `(buffer, frag_length_models, strand_models,
    n_total, calibration_payload)`.
  - `run_pipeline` calls `calibrate(...)` once after the scan.
  - `quant_from_buffer` accepts `calibration: CalibrationResult` (no
    `Stub` union). After `build_loci`, calls `assemble_priors(...)` and
    backfills the priors via `calibration.with_priors(...)`.
  - `_run_locus_em_partitioned` accepts `prior_weight_rna_per_locus`.

**Public package surface (`src/rigel/calibration/__init__.py`):**
```python
from ._orchestrator import calibrate
from ._fl_pool import (
    PoolFLModels, compute_pool_fl_models,
    DEFAULT_PRIOR_ESS, DEFAULT_QUALITY_THRESHOLD_GOOD,
    DEFAULT_QUALITY_THRESHOLD_WEAK,
)
from ._kappa import KappaEstimate, estimate_kappa, KAPPA_DEFAULT, KAPPA_MIN, KAPPA_MAX
from ._result import CalibrationResult, build_calibration_result
from .density_global import GlobalGdnaDensity, GlobalDensityTable, compute_global_densities
from .density_loco import shrink_to_loco
from .locus_prior import (
    C_BASE_DEFAULT,
    LocusGdnaEstimate, MultiLocusPrior, PriorTable,
    estimate_locus_gdna, assemble_multilocus_prior, assemble_priors,
    build_prior_weight_rna,
    FLAG_INTERGENIC_ZERO_LEFF, FLAG_INTRON_ZERO_LEFF,
    FLAG_EXON_INTRON_NO_ELIGIBLE, FLAG_PI_CLIPPED,
)
from .scan_payload import CalibrationScanPayload
```

**Legacy deletion** is performed in M8c (see above), not in M8a or M8b.
The v1 modules (`_simple.py`, `_categorize.py`, `_fl_mixture.py`,
`_fl_empirical_bayes.py`, the v1 `_result.py`) and their tests
(`tests/test_calibration_simple.py`, `tests/test_categorize.py`,
`tests/test_gdna.py`, `tests/test_gdna_harmonic_length.py`) are deleted
in M8c.

**Combined exit gate (M8a + M8b + M8c):** ≥ 14 new tests green; legacy
tests deleted (not skipped); golden outputs regenerated with documented
diffs; full suite green; benchmark matrix shows no regression.

### M9 — Validation + documentation

**Scope:** confirm calibration meets the design goals on representative
benchmarks; update user-facing documentation.

**Validation matrix:**
- Pristine RNA-seq (`gdna_none`): π̂_pool ≈ 0 (< 0.05).
- 50/50 dna20m: π̂_pool ∈ [0.45, 0.55].
- TA1 nRNA-siphon scenario: residual mRNA error < 5%.
- Hybrid-capture sample: gDNA flux density on internal exons matches
  synthetic ground truth within ±15% on a low-contamination scenario.
- Full Armis2 cluster benchmark vs salmon / kallisto on the 13-condition
  matrix (see `.github/copilot-instructions.md` §"Armis2 SLURM Cluster").

**Documentation:**
- `docs/MANUAL.md` — calibration section (1 page).
- `docs/METHODS.md` — full mathematical exposition (this plan condensed).
- `docs/parameters.md` — `pool_quality_thresholds`, `nrna_weight`,
  `c_base`.
- `CHANGELOG.md` — major release entry.

**Exit gate:** validation matrix passes; docs published; CHANGELOG entry
written; tag a release.

---

## 5. Code-quality contracts

### 5.1 C++

- C++17, `-O3`, LTO; namespaces under `rigel::` and `rigel::calibration::`.
- No raw `new`/`delete`; smart pointers for ownership.
- No heap allocation in `CalibrationAccumulator::observe()` after warm-up.
- All hot-path counters are `int64_t`.
- `SmallRegionSet`: inline storage `kInline = 16`, linear-scan dedup,
  spilled vector reused across fragments.

### 5.2 Python

- Frozen dataclasses with `slots=True` for all result objects.
- Type hints on all public functions.
- pandas + pyarrow only.
- No `assert` in production paths — raise explicit exceptions with
  actionable messages (e.g., the canonical "Rebuild the index" message).
- `pathlib.Path` for all file paths.
- `dataclasses.replace` for "modify" operations on frozen results.

### 5.3 Diagnostics

- Every milestone emits structured `INFO` logs summarizing key counts.
- All counts surface in `CalibrationResult` and through `summary.json`.
- Fallback paths are explicit and recorded
  (e.g., `KappaEstimate.fallback_used`, `LocusGdnaEstimate.fallback_flags`).

---

## 6. Protected test set

The following must pass after every milestone boundary; they cover end-to-end
pipeline behaviour independent of calibration intelligence:

- `tests/test_index.py`
- `tests/test_index_integrity.py`
- `tests/test_buffer.py`
- `tests/test_resolution.py`
- `tests/test_em_impl.py`
- `tests/test_estimator.py`
- `tests/test_pipeline_smoke.py`
- `tests/test_oracle_bam.py`
- `tests/test_cli.py`

Tests that legitimately cannot pass at a given milestone (e.g.,
`tests/test_golden_output.py` between M5 and M8 because outputs change
shape) get `pytest.mark.skip(reason="re-enabled at M<n>: ...")` with a
clear TODO. **No registry file**; the skip reason is the only enforcement.

---

## 7. Risk register

| # | Risk | Milestone | Mitigation |
|---|---|---|---|
| 1 | Bulk text transformation corrupts source files | all | §3 commit discipline; never apply whitespace regex to source; stage-then-commit before any global edit |
| 2 | Untracked work pile-up enables next data-loss | all | Each milestone is a single commit; no inter-milestone uncommitted scaffolding |
| 3 | `c_base = 10.0` constant is wrong for high/low-coverage regimes | M6+ | Single config knob; M9 benchmarks measure sensitivity; future plan can change formula |
| 4 | Hard-error on stale indexes alienates users with old indexes | M2, M8 | One-line rebuild instruction in error message; M9 docs flag the migration in CHANGELOG |
| 5 | Excluding multimappers under-estimates gDNA on multi-mappable territory | M3 | Documented; future mappability-adjusted `L_eff` recovers without re-introducing multimappers |
| 6 | Per-region NB MoM hits degenerate input on real data | M4 | Explicit fallback to `KAPPA_DEFAULT`, recorded on `KappaEstimate.fallback_used`; M9 monitors fallback rate |
| 7 | Single MultiLocus-level gDNA component misallocates across disjoint Locus intervals | M6 | Per-Locus diagnostics surfaced on `MultiLocusPrior.per_locus`; future per-Locus EM extension can use them |
| 8 | EM `prior_weight_rna` ABI change requires native re-spec | M5 | Bit-identical regression test with `None` argument; default is backwards-compat |
| 9 | nRNA siphon (TA1 scenario) not fixed by `nrna_weight=0` alone | M9 validation | If M9 benchmark fails the < 5% siphon criterion, escalate to a follow-up plan; do not silently widen tolerances |
| 10 | π̂(L) > 1 (capture leakage) | M6 | `WARN` + `FLAG_PI_CLIPPED`; clipped to 1.0; not blocking; M9 quantifies frequency |

---

## 8. Open questions deferred to future plans

- **`c_base(ℓ)` formula.** Constant `10.0` for now. A future plan can run
  the 4-candidate spike (`constant`, `sqrt_n`, `eb_inverse_var`,
  `coverage_weighted`) once we have real-data evidence on which scenario
  it matters for.
- **Per-`Locus` EM gDNA components.** Currently one shared component per
  `MultiLocus`. A future plan can emit one component per `Locus` using
  the diagnostics already stored on `MultiLocusPrior.per_locus`.
- **Mappability-adjusted $L_{\mathrm{eff}}$.** Multimapper exclusion biases
  density estimates downward on repeat tracts. A future plan can replace
  `L(R)` with $\int_R m(x) dx$ per the read-length configuration.
- **Annotation-gap inclusion in FL pool.** Default `False`. Future plans
  can flip if QC counts justify and oracle validation supports.
- **Per-region NB MLE.** Replace MoM with a single-parameter Newton on
  $\sum \log L$ if instability is observed on real data.

---

## 9. References

- Predecessor: `docs/calibration/calibration_v5_plan.md`.
- User TODO original sketch: `docs/TODO.md` §"Calibration".
- Memory notes informing this plan:
  - `/memories/repo/data-loss-incident-2026-05-01.md` — what we are
    rebuilding and why.
  - `/memories/repo/calibration-overlap-leff-2026-04.md` — `L_eff` overlap
    geometry.
  - `/memories/repo/multilocus-partition-design-2026-04.md` —
    anchor-transcript routing for `partition_units_to_loci`.
  - `/memories/repo/rigel-calibration-design.md` — original design
    rationale (intergenic-only is not enough; expression-based purity
    scoring; etc.).
  - `/memories/repo/gdna-siphon-root-cause-2026-03.md` — the TA1 nRNA
    siphon failure mode this plan is designed to fix.
  - `/memories/repo/calibration-fallback-fix-2026-04.md` — why
    `prior_weight_rna[nrna] = 0` matters.
