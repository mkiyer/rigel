# SPLICED_IMPLICIT Misclassification - Bug Report and Fix Plan

**Date:** 2026-05-09  
**Author:** Copilot, with diagnosis driven by user investigation  
**Status:** Revised plan (rev 4) - reviewed against current code; ready for implementation

---

## Revision history

- **rev 1 (2026-05-09):** initial plan; aggregate `intron_in_gaps` discriminant; new `implicit_splice_tolerance` parameter.
- **rev 2 (2026-05-09):** per-intron whole-containment discriminant; reused `BamScanConfig.boundary_tolerance`; iterated introns directly from the resolver exon CSR.
- **rev 3 (2026-05-09):** fixes several implementation errors in rev 2 and standardizes the single user-facing knob as **splicing anchor tolerance**:
  - canonical CLI flag becomes `--splicing-anchor-tolerance`;
  - canonical config field becomes `BamScanConfig.splicing_anchor_tolerance`;
  - legacy `--boundary-tolerance` / YAML `boundary_tolerance` are compatibility aliases only;
  - resolver setter becomes `FragmentResolver::set_splicing_anchor_tolerance(int K)`, not a separate implicit-splice tolerance;
  - the same resolved K is passed to both implicit-splice classification and calibration boundary-flux clearance;
  - the implicit-splice reference C++ no longer uses an invalid `return;` inside `_resolve_core`, preserves the existing no-chimera behavior, sorts/merges blocks defensively, and requires positive intron-gap overlap before applying K slack;
  - native ZS label for `SPLICE_ARTIFACT` is aligned with Python (`"splice_artifact"`, not `"spliced_artifact"`);
  - summary/config persistence is corrected to the actual `summary.json` structure.
- **rev 5 (2026-05-09):** drops back-compat for `boundary_tolerance`. The
  CLI/YAML/config/payload/summary keys are renamed unconditionally to
  `splicing_anchor_tolerance`. There is no `--boundary-tolerance` alias,
  no YAML alias, and no read-side fallback for legacy summaries. Tests
  that referenced the old name are updated.

- **rev 4 (2026-05-09):** tightens the implementation for performance and maintainability:
  - implicit-splice detection uses reusable `ResolverScratch` buffers and small helper functions instead of a large inline `_resolve_core` block;
  - transcript introns are queried by binary search over the existing per-transcript exon CSR;
  - scanner-sorted exon blocks take a no-copy fast path, with defensive copy/sort only for unsorted Python-facing calls;
  - no new cgranges intron index is planned unless profiling shows the CSR path is a bottleneck.

---

## 1. Summary

The C++ `_resolve_core` function classifies fragments as `SPLICE_IMPLICIT`
(splice junction inferred from a transcript-space projection shorter than the
genomic span) using a discriminant that **does not distinguish intron content
inside aligned blocks from intron content inside the unsequenced paired-end
gap**. The intended semantic is:

> An annotated intron of at least one candidate transcript is contained in the
> paired-end gap, allowing a small splice-anchor tolerance K.

The implemented test is instead:

> The total intron content of at least one candidate transcript inside the
> fragment outer span `[gstart, gend]` is positive.

For a true gDNA fragment whose contiguously aligned reads overlap a region that
some annotated transcript calls intronic, the implemented test trips. The
fragment is then:

- hard-gated out of the gDNA likelihood in
  [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp)
  because the gDNA per-hit branch requires `splice_type == SPLICE_UNSPLICED`;
- excluded from boundary-flux fan-out in
  [src/rigel/native/calibration/accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp)
  because `flux_eligible = (splice_type == SPLICE_UNSPLICED)`;
- forced to compete only as RNA, almost always landing in nRNA because true
  gDNA fragments lack mRNA support.

Empirically, the misallocation is small as a fraction of the gDNA pool
(~0.15% of total true-gDNA fragments per condition in
`scripts/debug/implicit_splice_antisense_probe.py`), but it is fully one-sided:
the misclassified true-gDNA fragments are dumped into RNA, mostly nRNA, and
they are exactly the unspliced boundary-crossing fragments whose absence
depresses `rho_ex / rho_ig`. This makes the bug a plausible contributor to the
residual deficit observed in the boundary-tolerance investigation.

The annotated-BAM `ZS` tag also lacks native labels for `SPLICE_IMPLICIT` and
`SPLICE_ARTIFACT`, which both currently render as `"unknown"` from the C++ BAM
writer. Python's `_splice_type_label` already handles these through
`SpliceType(code).name.lower()`.

---

## 2. Root cause

### 2.1 The current discriminant

In [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h),
the current implicit-splice block is:

```cpp
// --- SRD v2: SPLICED_IMPLICIT detection ---
if (cr.splice_type == SPLICE_UNSPLICED && exons.size() >= 2) {
    int32_t genomic_span = cr.genomic_footprint;
    for (const auto& kv : cr.frag_length_map) {
        int32_t fl = kv.second;
        if (fl > 0 && (genomic_span - fl) > 0) {
            cr.splice_type = SPLICE_IMPLICIT;
            break;
        }
    }
}
```

`cr.frag_length_map[t]` is filled by `compute_frag_lengths`, which projects only
the outermost endpoints `(gstart, gend)` into transcript space:

```cpp
for (int32_t t : t_inds) {
    int32_t tx_s = genomic_to_tx_pos(gstart, t);
    int32_t tx_e = genomic_to_tx_pos(gend,   t);
    int32_t fl   = std::abs(tx_e - tx_s);
    if (fl > 0) result[t] = fl;
}
```

`genomic_to_tx_pos` collapses every transcript intron contained in the outer
span, regardless of whether the intron lies in an aligned block or in the
paired-end gap. Therefore:

```text
genomic_span - fl_t = total transcript-t intron content in [gstart, gend]
```

That is not the same quantity as intron content inside the PE gap. It includes
introns that the CIGAR explicitly contradicts by aligning through them.

### 2.2 Empirical fingerprint

`scripts/debug/implicit_splice_antisense_probe.py` was run across synthetic
conditions with K=3:

| condition | true gDNA | inferred-implicit true-gDNA | assigned to RNA | mRNA | nRNA |
|---|---:|---:|---:|---:|---:|
| `gdna_low_ss_0.99_nrna_none` | 100,000 | 155 | 155 | 1 | 154 |
| `gdna_med_ss_0.99_nrna_none` | 500,000 | 739 | 739 | 2 | 737 |
| `gdna_equal_ss_0.99_nrna_none` | 1,000,000 | 1,593 | 1,593 | 4 | 1,589 |
| `gdna_high_ss_0.99_nrna_none` | 2,000,000 | 3,008 | 3,008 | 5 | 3,003 |

Geometry of the misclassified true-gDNA fragments:

- strict full-intron-contained-in-PE-gap among true-gDNA inferred implicit: **0**;
- any intron overlap across the full outer span among inferred implicit: **100%**;
- true-RNA inferred-implicit fragments do meet the strict intron-in-gap definition.

The false positives are ordinary gDNA fragments by length (`gdna_high` truth FL
median 422 bp, q90 537 bp, q99 629 bp, max 740 bp), not tail artifacts.

---

## 3. Fix design

### 3.1 Canonical parameter: splicing anchor tolerance

Use one conceptual and user-facing parameter everywhere:

```text
splicing_anchor_tolerance = K bp
```

It governs two related aligner-anchor tolerances:

1. implicit-splice classification: how far an annotated intron boundary may
   protrude outside a paired-end gap while still being considered contained;
2. calibration boundary-flux clearance: how much aligned sequence must clear an
   exon-intron boundary before a fragment contributes to boundary-flux counts and
   EXON|INTRON region-mask observations.

Canonical public names:

- CLI: `--splicing-anchor-tolerance K`
- YAML/config key: `splicing_anchor_tolerance`
- Python dataclass field: `BamScanConfig.splicing_anchor_tolerance`
- Native resolver setter: `set_splicing_anchor_tolerance(K)`
- Summary/config output key: `splicing_anchor_tolerance`

This is a hard rename. There is no `--boundary-tolerance` alias, no YAML
alias, and no read-side fallback for old serialized payloads or
`summary.json` files. All occurrences of `boundary_tolerance` in the
code base, tests, debug scripts, and docs are updated in one mechanical
pass.

Important semantic detail: the **implicit-splice** predicate uses raw K, so
`K = 0` means strict containment. The **calibration boundary-flux** code keeps
its existing integer-coordinate compatibility rule `q(K) = max(K, 1)`, so
`K = 0` still reproduces pre-2026.05 strict-crossing behavior bit-for-bit.
Those are two consumers of the same parameter, not two parameters.

### 3.2 Per-intron whole-containment discriminant

Replace the current fragment-length-shortening discriminant with an explicit
per-candidate-intron, per-PE-gap containment test.

For each fragment with sorted aligned blocks `{[b_s^k, b_e^k)}`:

1. Merge overlapping/adjacent aligned blocks into non-overlapping blocks
   `{B^j}`. Sort defensively by `(ref_id, start, end)` before merging; the BAM
   scanner already sorts `AssembledFragment.exons`, but the Python-facing
   `resolve()` / `resolve_fragment()` path should not rely on caller order.
2. Derive PE gaps between consecutive merged blocks on the same reference:
   `G^j = [B_e^j, B_s^{j+1})` when `B_s^{j+1} > B_e^j`.
3. Query transcript introns directly from the resolver exon CSR (`I_t^i = [exon_ends_[i], exon_starts_[i + 1])`).
4. Require positive half-open overlap between intron and PE gap before applying tolerance; this prevents a large K from classifying a nearby-but-disjoint intron as implicit.
5. Apply one-sided boundary slack:

   ```text
   overlaps(I, G) && (G.start - I.start) <= K && (I.end - G.end) <= K
   ```

6. Classify as `SPLICE_IMPLICIT` as soon as any candidate transcript has any
   intron-gap pair satisfying the predicate.

### 3.3 Hot-path implementation choice

Use the existing per-transcript exon CSR and binary search. Do **not** add a
new cgranges intron index for the first implementation.

Rationale:

- `_resolve_core` already knows the candidate transcripts in `cr.t_inds`; a
  global intron cgranges query would have to rediscover interval hits and then
  intersect transcript sets back against `cr.t_inds`.
- The exon CSR is already built at index-load time for every transcript, with
  exons in genomic order. Consecutive exon ends/starts are monotone intron
  starts/ends, so `std::lower_bound` / `std::upper_bound` can jump directly to
  introns near the PE gap.
- The expected workload is small and local: most paired-end fragments have one
  unsequenced gap, and each candidate transcript contributes only the introns
  whose starts/ends fall near that gap.
- A new cgranges index would add memory, build-time, binding, and scratch-buffer
  complexity to a path that can be implemented as cache-friendly pointer math
  over arrays already resident in `FragmentResolver`.

Complexity target for the implicit-splice check:

```text
O(n_gaps * n_candidate_transcripts * (log n_introns_t + n_local_introns))
```

not `O(n_gaps * n_candidate_transcripts * n_introns_t)` and not a global
interval-index query per gap. If profiling later shows this helper is material
in mega-loci, the fallback design is an optional intron cgranges index keyed by
intron interval with transcript-set labels; that should be a measured second
step, not the default design.

Implementation style:

- Add tiny helpers in [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h), rather than growing `_resolve_core` with another monolithic block.
- Add reusable scratch vectors to `ResolverScratch`, for example
  `implicit_blocks` and `implicit_gaps`, so the defensive sort fallback and gap
  collection do not allocate on every fragment after warm-up.
- Fast-path the common BAM-scanner input, which is already sorted/merged: only
  copy and sort into scratch when `std::is_sorted` says the caller supplied
  unsorted blocks.
- Keep the helper read-only with respect to candidate sets and exon CSR. It
  should only return `true` / `false`; the caller alone mutates `cr.splice_type`.

Reference C++ shape for [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h):

```cpp
struct GapBlock {
  int32_t ref_id;
  int32_t start;
  int32_t end;
};

inline bool transcript_has_implicit_intron_in_gap(
  int32_t t, int32_t gap_start, int32_t gap_end, int32_t K) const
{
  if (t < 0 || t + 1 >= static_cast<int32_t>(exon_offsets_.size())) return false;

  const int32_t begin = exon_offsets_[t];
  const int32_t end = exon_offsets_[t + 1];
  const int32_t n_introns = end - begin - 1;
  if (n_introns <= 0 || gap_end <= gap_start) return false;

  const int32_t* intron_starts = exon_ends_.data() + begin;
  const int32_t* intron_ends = exon_starts_.data() + begin + 1;

  const int64_t min_start = static_cast<int64_t>(gap_start) - K;
  const int64_t max_end = static_cast<int64_t>(gap_end) + K;

  const int32_t by_start = static_cast<int32_t>(
    std::lower_bound(intron_starts, intron_starts + n_introns, min_start)
    - intron_starts);
  const int32_t by_end = static_cast<int32_t>(
    std::upper_bound(intron_ends, intron_ends + n_introns, gap_start)
    - intron_ends);

  for (int32_t i = std::max(by_start, by_end); i < n_introns; ++i) {
    const int32_t is = intron_starts[i];
    if (is >= gap_end) break;        // no positive overlap with later introns
    const int32_t ie = intron_ends[i];
    if (static_cast<int64_t>(ie) > max_end) break;
    if (ie > is) return true;        // all containment/overlap tests hold
  }
  return false;
}

inline bool has_implicit_splice_gap(
  const std::vector<ExonBlock>& exons,
  const std::vector<int32_t>& candidate_t,
  ResolverScratch& scratch) const
{
  if (exons.size() < 2 || candidate_t.empty() || !has_exon_index()) return false;

  auto& gaps = scratch.implicit_gaps;
  gaps.clear();

  auto exon_block_less = [](const ExonBlock& a, const ExonBlock& b) {
    if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
    if (a.start != b.start) return a.start < b.start;
    if (a.end != b.end) return a.end < b.end;
    return a.strand < b.strand;
  };

  const std::vector<ExonBlock>* ordered = &exons;
  if (!std::is_sorted(exons.begin(), exons.end(), exon_block_less)) {
    auto& blocks = scratch.implicit_blocks;
    blocks.assign(exons.begin(), exons.end());
    std::sort(blocks.begin(), blocks.end(), exon_block_less);
    ordered = &blocks;
  }

  int32_t cur_ref = (*ordered)[0].ref_id;
  int32_t cur_end = (*ordered)[0].end;
  for (size_t k = 1; k < ordered->size(); ++k) {
    const auto& block = (*ordered)[k];
    if (block.ref_id != cur_ref) {
      cur_ref = block.ref_id;
      cur_end = block.end;
    } else if (block.start > cur_end) {
      gaps.push_back({cur_ref, cur_end, block.start});
      cur_end = block.end;
    } else {
      cur_end = std::max(cur_end, block.end);
    }
  }

  const int32_t K = splicing_anchor_tolerance_;
  for (int32_t t : candidate_t) {
    for (const GapBlock& gap : gaps) {
      if (transcript_has_implicit_intron_in_gap(t, gap.start, gap.end, K)) {
        return true;
      }
    }
  }
  return false;
}

// In _resolve_core, after cr.t_inds and chimera_type are known:
if (cr.splice_type == SPLICE_UNSPLICED &&
  cr.chimera_type == CHIMERA_NONE &&
  has_implicit_splice_gap(exons, cr.t_inds, scratch)) {
  cr.splice_type = SPLICE_IMPLICIT;
}
```

Key corrections relative to rev 2:

- Do **not** `return;` when `gaps.empty()`. `_resolve_core` is a `bool`
  function and still must set `genomic_start`, clean scratch, and return `true`.
- Preserve current no-chimera behavior with `cr.chimera_type == CHIMERA_NONE`.
- Require positive intron-gap overlap before applying K slack.
- Use the resolver's single `splicing_anchor_tolerance_` member, not a separate
  `implicit_splice_tolerance_`.
- Avoid per-fragment heap churn by putting temporary blocks/gaps in
  `ResolverScratch`.
- Avoid scanning all introns in every candidate transcript; use binary search to
  inspect only introns that can overlap the PE gap and satisfy K slack.

### 3.4 Why per-intron, not aggregate

A true gDNA fragment whose 200 bp PE gap overlaps a small slice of a 50 kb
annotated intron would satisfy an aggregate `intron_in_gap > K` test even though
no specific intron is contained in the gap. The per-intron rule rejects that
class because at least one intron boundary is far outside the PE gap. This is
the empirical false-positive class from the diagnostic.

### 3.5 Why use the exon CSR instead of cgranges

`genomic_to_tx_pos` is correct for fragment-length projection because endpoint
bases in introns must count toward physical fragment length. It is the wrong
primitive for deciding which intron intervals lie in which gaps. The exon CSR
already provides authoritative genomic intron coordinates via consecutive exon
ends/starts.

The CSR path also keeps the implementation maintainable: no new index file, no
new native builder, no transcript-set payload for intron intervals, and no new
thread-safety rules beyond the existing `ResolverScratch` convention. It should
be faster for the common case because it restricts work to already-known
candidate transcripts and performs local monotone-array searches.

cgranges remains a reasonable contingency only if profiling shows pathological
mega-loci spend meaningful time in `has_implicit_splice_gap`. In that case, add
a measured phase with a global intron cgranges index labelled by transcript-set
offsets and compare it against the CSR binary-search helper before committing
the added complexity.

`cr.frag_length_map` remains unchanged. It still represents the projected
transcript-space fragment length conditional on a spliced interpretation. Only
the `SPLICE_IMPLICIT` classification changes.

### 3.6 ZS tag completeness

Update native `splice_type_label` in
[src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) to
match Python label semantics:

```cpp
case SPLICE_UNSPLICED:       return "unspliced";
case SPLICE_SPLICED_UNANNOT: return "spliced_unannot";
case SPLICE_SPLICED_ANNOT:   return "spliced_annot";
case SPLICE_IMPLICIT:        return "spliced_implicit";
case SPLICE_ARTIFACT:        return "splice_artifact";
default:                     return "unknown";
```

Note the artifact label is `"splice_artifact"` because Python currently emits
`SpliceType.SPLICE_ARTIFACT.name.lower()`. Do not introduce the inconsistent
label `"spliced_artifact"` unless the Python enum is renamed in the same change.

Also expose `SPLICE_IMPLICIT` and `SPLICE_ARTIFACT` from `_resolve_impl` in
[src/rigel/native/resolve.cpp](../../src/rigel/native/resolve.cpp); the native
constants already exist in [src/rigel/native/constants.h](../../src/rigel/native/constants.h).

---

## 4. Implementation plan

### Commit 1 - Native splice labels and constants

**Files:**

- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp)
- [src/rigel/native/resolve.cpp](../../src/rigel/native/resolve.cpp)
- [tests/test_annotate.py](../../tests/test_annotate.py) or a native annotation test if one exists

Changes:

- Add `SPLICE_IMPLICIT` and `SPLICE_ARTIFACT` cases to native `splice_type_label`.
- Use `"splice_artifact"`, matching Python.
- Expose `SPLICE_IMPLICIT` and `SPLICE_ARTIFACT` as module-level `_resolve_impl`
  attrs.
- Add/extend a test that native annotated BAM labels no longer render these
  splice classes as `"unknown"`.

Acceptance:

- Recompile with `conda activate rigel && pip install --no-build-isolation -e .`.
- Run focused annotation/splice tests, then the full suite if time permits.
- No quantification output changes are expected from this metadata-only commit.

### Commit 2 - Rename/converge the tolerance parameter

**Files:**

- [src/rigel/config.py](../../src/rigel/config.py)
- [src/rigel/cli.py](../../src/rigel/cli.py)
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py)
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp)
- [src/rigel/native/calibration/accumulator.h](../../src/rigel/native/calibration/accumulator.h)
- [src/rigel/native/calibration/accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp)
- [src/rigel/calibration/_exposure.py](../../src/rigel/calibration/_exposure.py)
- [src/rigel/calibration/density_global.py](../../src/rigel/calibration/density_global.py)
- [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py)
- [src/rigel/calibration/scan_payload.py](../../src/rigel/calibration/scan_payload.py)
- [src/rigel/calibration/_orchestrator.py](../../src/rigel/calibration/_orchestrator.py)
- [src/rigel/calibration/_result.py](../../src/rigel/calibration/_result.py)
- debug/benchmark scripts that construct `BamScanConfig(boundary_tolerance=...)`
- [docs/parameters.md](../parameters.md)
- [CHANGELOG.md](../../CHANGELOG.md)

Changes:

- Rename `BamScanConfig.boundary_tolerance` to
  `BamScanConfig.splicing_anchor_tolerance` (validation, docstrings).
- Rename CLI flag to `--splicing-anchor-tolerance`. **Remove** the old
  `--boundary-tolerance` flag entirely.
- Update YAML config loaders to accept only `splicing_anchor_tolerance`.
  **Remove** any `boundary_tolerance` key handling.
- Update `_PARAM_SPECS` to map `splicing_anchor_tolerance` to
  `scan.splicing_anchor_tolerance`.
- Update `summary.json` keys to `splicing_anchor_tolerance` in
  `configuration.scan`, `command.arguments`, and `calibration`.
- Rename calibration payload/result fields to `splicing_anchor_tolerance`
  with no fallback.
- Rename Python kwargs throughout calibration math (`density_global`,
  `_exposure`, `locus_prior`, `_orchestrator`) to
  `splicing_anchor_tolerance`.
- Rename native accumulator/scanner variable names and binding docs to
  `splicing_anchor_tolerance`. Preserve the existing `q(K)=max(K,1)`
  calibration-side semantic.
- Update tests, debug scripts, and benchmark configs to use the new name.

Acceptance:

- CLI tests cover the canonical flag and the YAML canonical key.
- A grep for `boundary_tolerance` (any case) in the code base returns no
  hits.
- Summary/config tests assert the canonical key is written.
- Calibration tests still prove `K = 0` reproduces strict-crossing
  behavior.

### Commit 3 - Per-intron implicit-splice discriminant

**Files:**

- [src/rigel/native/resolve_context.h](../../src/rigel/native/resolve_context.h)
- [src/rigel/native/resolve.cpp](../../src/rigel/native/resolve.cpp)
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py)
- [tests/test_resolution.py](../../tests/test_resolution.py)

Changes:

- Add `int32_t splicing_anchor_tolerance_ = 0;` to `FragmentResolver`.
  The class currently keeps resolver state public, so do not describe this as a
  private member unless the class is deliberately refactored.
- Add `set_splicing_anchor_tolerance(int K)` with `K >= 0` validation and bind it
  in `_resolve_impl`.
- Add reusable `ResolverScratch` storage for implicit-splice gap detection, such
  as `implicit_blocks` and `implicit_gaps`, and reserve a small initial capacity
  in the constructor. These vectors must be cleared and reused per fragment.
- Implement the implicit-splice geometry as small helpers:
  `collect_pe_gaps`, `transcript_has_implicit_intron_in_gap`, and
  `has_implicit_splice_gap`. Keep `_resolve_core` to a short gate plus one
  setter of `cr.splice_type`.
- Use binary search over per-transcript intron starts/ends derived from the exon
  CSR. Do not add a cgranges intron index in this commit.
- In `scan_and_buffer`, call
  `resolve_ctx.set_splicing_anchor_tolerance(int(scan.splicing_anchor_tolerance))`
  before `_NativeBamScanner` construction, alongside the existing resolver setup.
- In `_resolve_core`, replace the current fragment-length-shortening
  implicit-splice block with the per-intron/per-gap predicate in Section 3.2.
- Do not change `compute_frag_lengths`.

New unit tests in [tests/test_resolution.py](../../tests/test_resolution.py):

Use function-scoped indexes or reset the resolver tolerance after each test.
Avoid mutating a session-scoped `TranscriptIndex` fixture without cleanup,
because `FragmentResolver` state persists across tests.

| case | blocks | candidate intron geometry vs PE gap | K | expected splice type |
|---|---|---|---:|---|
| (a) intron strictly inside PE gap | 2 | `G_s <= I_s`, `I_e <= G_e` | 0 | `SPLICED_IMPLICIT` |
| (b) intron strictly inside aligned block | 2 | no positive intron-gap overlap | 0 | `UNSPLICED` |
| (c) intron protrudes 2 bp left of gap | 2 | positive overlap, left protrusion 2 | 3 | `SPLICED_IMPLICIT` |
| (d) intron protrudes 5 bp left of gap | 2 | positive overlap, left protrusion 5 | 3 | `UNSPLICED` |
| (e) 200 bp slice of a 50 kb intron lies in gap | 2 | both intron boundaries far outside gap | 3 | `UNSPLICED` |
| (f) 4 bp microintron entirely in gap | 2 | strict containment | 3 | `SPLICED_IMPLICIT` |
| (g) 4 bp microintron entirely in block | 2 | no positive intron-gap overlap | 3 | `UNSPLICED` |
| (h) two-intron transcript, one intron contained | 2 | any one intron-gap pair satisfies predicate | 3 | `SPLICED_IMPLICIT` |
| (i) single-block fragment | 1 | no PE gap | 3 | `UNSPLICED` |
| (j) CIGAR-spliced annotated fragment | spliced | implicit gate not entered | 3 | `SPLICED_ANNOT` |
| (k) CIGAR-spliced unannotated fragment | spliced | implicit gate not entered | 3 | `SPLICED_UNANNOT` |
| (l) two candidates, only second has contained intron | 2 | any-candidate semantics | 3 | `SPLICED_IMPLICIT` |
| (m) nearby but disjoint intron within K bp of gap | 2 | no positive overlap | 3 | `UNSPLICED` |
| (n) chimeric multi-block fragment | 2+ | would otherwise satisfy geometry | 3 | not promoted to `SPLICED_IMPLICIT` |

Acceptance:

- Recompile after C++ changes.
- Focused resolver tests pass.
- Add a lightweight performance assertion or benchmark note for a synthetic
  high-candidate locus: the helper should scale with local introns near the gap,
  not total introns in the locus.
- Full suite passes after any intentional golden updates.
- If goldens shift, inspect before regenerating. Expected shifts are fewer
  fragments routed to nRNA in paired-end gDNA contamination scenarios; scenarios
  without paired-end gDNA contamination should be unchanged or explainable.

### Commit 4 - Synthetic validation and documentation

Re-run [scripts/debug/implicit_splice_antisense_probe.py](../../scripts/debug/implicit_splice_antisense_probe.py)
with default `--splicing-anchor-tolerance 3` against:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic
```

Hard correctness expectations:

1. Inferred-implicit true-gDNA counts drop from 155 / 739 / 1,593 / 3,008 to a
   small residual explained by genuine intron-in-gap geometry. Any non-zero
   residual should be auditable per fragment.
2. Total nRNA assignments in `gdna_high_ss_0.99_nrna_none` drop by roughly the
   corrected count, allowing that some corrected fragments may compete and land
   elsewhere.
3. CIGAR-spliced fragment classification is unchanged.

Hypothesis-framed targets:

4. Re-run [scripts/debug/boundary_side_geometry_probe.py](../../scripts/debug/boundary_side_geometry_probe.py).
   The previous headline `rho_ex / rho_ig` was 0.924. A monotone improvement is
   positive evidence; `>= 0.95` is a useful target but not a correctness gate.
5. Run `python -m scripts.benchmarking analyze` and confirm transcript-level
   relative error is no worse, ideally modestly improved in high-gDNA conditions.
6. Document headline results in [CHANGELOG.md](../../CHANGELOG.md) under
   `[Unreleased] - Fixed` and update [docs/parameters.md](../parameters.md) to
   describe `--splicing-anchor-tolerance`.

---

## 5. Why this fix is correct and bounded

1. **Defers to the alignment.** An intron inside an aligned block is
   contradicted by the CIGAR and should not by itself promote a fragment to
   `SPLICE_IMPLICIT`.
2. **Preserves the strict K=0 semantic.** At `K=0`, implicit means at least one
   candidate transcript has a positive intron interval wholly contained in a PE
   gap.
3. **Uses one tolerance knob.** The same K handles small splice-anchor placement
   errors for both implicit-splice detection and calibration boundary flux.
4. **Does not change scoring gates.** gDNA scoring and calibration already key
   off `splice_type == SPLICE_UNSPLICED`; fixing classification lets downstream
   behavior self-correct.
5. **Keeps fragment-length semantics stable.** `frag_length_map` remains the
   outer-endpoint transcript-space projection conditional on a spliced
   interpretation.
6. **Has bounded cost.** The new work is roughly `n_gaps * n_candidates * (log n_introns_t + n_local_introns)` for unspliced, non-chimeric multi-block fragments only. Typical paired-end fragments have one gap, and the binary-search window avoids scanning all introns in large transcripts.

---

## 6. Risk register

1. **Parameter rename blast radius.** `boundary_tolerance` appears in config,
   CLI, pipeline, calibration payloads, density denominators, native scanner
   bindings, docs, debug scripts, and changelog text. The rename is hard:
   one mechanical pass plus a final `grep -ri boundary_tolerance src tests
   scripts docs` that must return zero hits.
2. **Legacy configs.** Pre-1.0 project; no compatibility window. Users who
   regenerate runs with old YAMLs containing `boundary_tolerance` will get
   a clear unknown-key error from the YAML loader. CHANGELOG documents the
   rename.
3. **Chimeras.** Current implicit detection effectively skips chimeras because
   `frag_length_map` is only filled for `CHIMERA_NONE`. The new CSR-based code
   must keep that gate explicitly.
4. **Near-but-disjoint introns.** K slack must not classify introns that do not
   overlap the PE gap at all. The explicit positive-overlap check is required.
5. **Microintrons.** No `min_implicit_intron_bp` is introduced here. A
   microintron wholly inside the PE gap remains a valid implicit splice.
6. **Native/Python label consistency.** The native C++ label for artifact must
   be `"splice_artifact"`, matching Python's enum-derived label.
7. **Session-scoped resolver state in tests.** The new resolver tolerance setter
   mutates resolver state. Tests must isolate or reset it.
8. **Hot-path regression.** This check runs inside `_resolve_core`, so avoid allocation and global interval queries. Temporary vectors belong in `ResolverScratch`, and the intron lookup should use binary search over the existing exon CSR. Add cgranges only after profiling shows a real regression.

---

## 7. Acceptance summary

A successful implementation will:

- add native ZS labels for `spliced_implicit` and `splice_artifact`;
- expose native `SPLICE_IMPLICIT` and `SPLICE_ARTIFACT` constants from
  `_resolve_impl`;
- standardize the CLI/config parameter as `splicing_anchor_tolerance`
  (hard rename; no `boundary_tolerance` aliases anywhere);
- pass exactly one resolved K into both the resolver and calibration scanner;
- replace the implicit-splice discriminant with the per-intron/per-PE-gap
  containment test using raw K slack and positive intron-gap overlap;
- implement the containment test as reusable scratch-backed helpers using binary
  search over the existing exon CSR, not a new cgranges intron index;
- preserve chimeric and CIGAR-spliced classification behavior;
- add the geometry tests in Section 4 Commit 3;
- regenerate goldens only after manual inspection of intentional output shifts;
- validate on the synthetic sweep that true-gDNA to nRNA misallocation drops to
  near zero and boundary-density ratios improve or remain explainable.
