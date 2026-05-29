# Accumulator v2 — Phase 1 Audit Memo

Date: 2026-05-29
Scope: gate artifact for Phase 1 of
[`01_implementation.md`](01_implementation.md). Documents the **exact
current behavior** of the native calibration accumulator and the
Python-layer boundary derivation, locks the v2 spec, and identifies
any legacy bugs that must be fixed before the Phase 5 parity gate.

**Bottom line.** Legacy code is correct. **No bugs identified, no
legacy fixes required.** Phase 1 proceeds with Python reference +
xfail spec tests only. The v2 work is a *substrate extension* (adds
per-boundary flux + per-boundary FL histograms) and a *storage
reshape* (region-keyed boundary slots → boundary-keyed) on top of
already-correct mass accounting.

---

## 1. Files audited

| File | Role |
|---|---|
| [`src/rigel/native/calibration/accumulator.h`](../../src/rigel/native/calibration/accumulator.h) | Native accumulator struct + API |
| [`src/rigel/native/calibration/accumulator.cpp`](../../src/rigel/native/calibration/accumulator.cpp) | `CalibrationAccumulator::observe` implementation |
| [`src/rigel/native/calibration/region_index.h`](../../src/rigel/native/calibration/region_index.h) | Region overlap-query interface |
| [`src/rigel/native/calibration/region_signature.h`](../../src/rigel/native/calibration/region_signature.h) | Channel/signature index constants |
| [`src/rigel/native/bam_scanner.cpp`](../../src/rigel/native/bam_scanner.cpp) | Per-fragment finalize block → `cal_acc.observe(...)` |
| [`src/rigel/calibration/region_count_ledger.py`](../../src/rigel/calibration/region_count_ledger.py) | Python ledger view |
| [`src/rigel/calibration/boundaries.py`](../../src/rigel/calibration/boundaries.py) | Python-side boundary-keyed derivation |
| [`src/rigel/calibration/scan_payload.py`](../../src/rigel/calibration/scan_payload.py) | Payload + payload-level invariants |

Tests inventoried: `tests/test_boundaries.py` (will be `tests/test_boundary_model.py`/`tests/test_boundary_sweep.py` per design doc's older naming — confirm in Phase 6), `tests/test_region_count_ledger.py`, `tests/test_region_partition.py`, `tests/test_region_persist.py`, plus calibration integration tests.

---

## 2. Per-fragment vs per-block normalization

**Verdict: correct, per-fragment.** Confirmed at three sites in
`accumulator.cpp`:

- Lines 48–53: `total_aligned_bp = Σ block_len_k` computed once
  outside the per-block loop.
- Line 76: `const double inv_total = 1.0 / static_cast<double>(total_aligned_bp);`
- Line 133 (per-block deposit): `w = overlap_bp * inv_total` — overlap of *this block* with *this region*, divided by per-fragment $L$.

Sum-of-deposits per fragment = $\sum_\text{blocks} \sum_\text{regions overlapping block} (\text{overlap\_bp} / L) = L / L = 1.0$.

**No bug.** The v2 reference implementation must reproduce this
exactly.

---

## 3. Boundary-mass-mirroring inventory

Today's `RegionCountLedger` stores boundary mass **region-keyed**:
`region_counts[r, channel]` where the 12 channels encode
(compartment × splice × strand), compartment ∈ {CONTAINED,
BOUNDARY_LEFT, BOUNDARY_RIGHT}.

For a fragment crossing the boundary at position $p$ that separates
regions $R_a$ (left) and $R_b$ (right):

- The block in $R_a$ deposits into `region_counts[R_a, BOUNDARY_RIGHT, ...]`
  with mass $\ell_a / L$ — because the fragment span exceeds $R_a$'s right edge.
- The block in $R_b$ deposits into `region_counts[R_b, BOUNDARY_LEFT, ...]`
  with mass $\ell_b / L$ — because the fragment span begins before $R_b$'s left edge.

The two slots are **partitioned, not duplicated** — they receive
*different* mass amounts (the per-block overlap shares). Boundary
$p$'s aggregate "evidence" is reconstructed by
[`boundaries.py`](../../src/rigel/calibration/boundaries.py) joining
`region_counts[R_a, BOUNDARY_RIGHT, ...]` with
`region_counts[R_b, BOUNDARY_LEFT, ...]` into one record.

**v2 mapping.** The legacy `boundary_right_*` of $R_a$ = v2
`boundary[p].mass_left_*` (block was on the *left* of the boundary).
The legacy `boundary_left_*` of $R_b$ = v2 `boundary[p].mass_right_*`
(block was on the *right* of the boundary). This sign-flip is the
single most important mapping; the parity test in Phase 4 must apply
it exactly.

---

## 4. Fully-spans-region case

**Verdict: implemented as documented in CLAUDE.md, w/2 split.**

`accumulator.cpp` lines 159–178 explicitly handle the case where a
single block's overlap with region $R$ has both `cross_left` and
`cross_right` true (i.e. the fragment span exceeds *both* of $R$'s
edges). Behavior:

- `region_counts[R, BOUNDARY_LEFT, ...]` += `0.5 * w`
- `region_counts[R, BOUNDARY_RIGHT, ...]` += `0.5 * w`
- FL pool += `w` once (via the LEFT compartment slot only, to avoid double-counting).

Where `w = overlap_bp_of_block_in_R / total_aligned_bp_of_fragment`.

**v2 reproduction.** In v2 §4.5.4, this case is handled by
"decompose the straddling block into per-region slices and apply
the §4.3 rule as if they were adjacent block-fragments with no
inter-block gap." For a single block straddling exactly one whole
region $R$:

- Left slice (in $R_{prev}$) contributes $\ell_{prev}/L$ to $R_{prev}$.right_boundary.mass_left.
- Middle slice (in $R$) contributes $\ell_R/L$ to both $R$.left_boundary.mass_right AND $R$.right_boundary.mass_left.
- Right slice (in $R_{next}$) contributes $\ell_{next}/L$ to $R_{next}$.left_boundary.mass_right.

**Important compatibility note.** This v2 attribution is *different*
from the legacy "w/2 of $R$'s block-overlap-mass goes to each of
$R$'s left/right boundary slots". Under v2, the middle slice of
length $\ell_R$ deposits *all* $\ell_R/L$ on each side (not $\ell_R/(2L)$),
because both sides have observed evidence. The legacy w/2 was a
heuristic for the "I don't know which side to credit" problem; v2
resolves the ambiguity by crediting both sides fully whenever the
block actually straddles a boundary.

**Implication for Phase 5 parity gate.** The legacy and v2 outputs
will *not* match bit-for-bit in any locus that contains a
fully-spans-region fragment. The Phase 5 parity test must (a)
identify all loci where fully-spanning fragments occur, (b) exclude
those boundary slots from the equality assertion or (c) apply a
custom comparator that reproduces the legacy w/2 from v2's full-credit
deposits. Decision deferred to Phase 4.

Likelihood of impact: **low** — fully-spans-region requires
$\text{region length} < \text{fragment insert length}$, which for typical
RNA-seq fragments (~200–400 bp) only happens in very short exonic
regions (< 200 bp). Will measure prevalence in the Phase 1 reference
test corpus.

---

## 5. Splice-flag attribution

**Verdict: per-fragment.** `accumulator.cpp` line 70 reads:

```cpp
const int splice_idx = (splice_type == SPLICE_UNSPLICED)
                          ? region_signature::kSpliceUnspliced
                          : region_signature::kSpliceSpliced;
```

The `splice_type` parameter is passed once per fragment from
`bam_scanner.cpp` line 1648, set upstream by the resolver based on
whether any junction in the fragment is recognized as spliced.

**v2 reproduction.** Design doc §4.3 sketched per-junction splice
flags. **Phase 1 decision: inherit the legacy per-fragment flag.**
Per-junction would require the resolver to expose per-junction
classifications, which is a separate refactor outside this scope.
The reference implementation and tests use a single per-fragment
spliced flag matching legacy behavior. Update
[`00_design.md`](00_design.md) §4.3 to reflect this.

---

## 6. Strand attribution

**Verdict: per-fragment, strand-ambiguous rejected.**
`accumulator.cpp` lines 45–56:

```cpp
if (fragment_strand == STRAND_POS)      strand_idx = kChannelStrandPos;
else if (fragment_strand == STRAND_NEG) strand_idx = kChannelStrandNeg;
else { ++payload_.n_excluded_strand_ambig; return; }
```

The strand_idx is set once per fragment and used for all per-block
deposits. Strand-ambiguous fragments are dropped (not routed to
either +/- channel).

**v2 reproduction.** Match exactly: per-fragment strand sign,
ambiguous → drop.

---

## 7. Dedup of per-region per-compartment support

`accumulator.cpp` lines 182–219 sort+unique the `touched_contained_`,
`touched_left_`, `touched_right_` vectors and increment integer
support counters per (region, compartment). This guarantees a
region is counted **at most once per compartment per fragment**,
even if multiple blocks of the fragment touch the same region.

**v2 reproduction.** The v2 `Region.contained_*_*` counter is
incremented exactly once per fragment classified as "contained"
(by `+= 1`). The v2 `Boundary.flux_*_*` counter is incremented at
most once per boundary per fragment (§4.4 flux dedup). The v2 fl
histograms are incremented at most once per region/boundary per
fragment (§5).

---

## 8. Strand and splice channels — index parity

Legacy 12 channels: 3 compartments × 2 splice × 2 strands.
- channel index = `compartment * 4 + splice * 2 + strand` (verify in `region_signature.h`).

v2 4 region channels: spl×strand (compartment fixed to CONTAINED).
v2 16 boundary mass channels: side × spl × strand × 2 (left/right side).
v2 4 boundary flux channels: spl×strand.

Mapping table (legacy → v2):

| Legacy slot | v2 slot |
|---|---|
| `region_counts[r, CONTAINED, spl, strand]` | `regions[r].contained_{spl}_{strand} += w` (but v2 stores integer count, not float mass; equivalence holds because for contained fragments $\sum w = 1$) |
| `region_counts[r, BOUNDARY_RIGHT, spl, strand]` | `boundaries[right_of_r].mass_left_{spl}_{strand}` |
| `region_counts[r, BOUNDARY_LEFT, spl, strand]` | `boundaries[left_of_r].mass_right_{spl}_{strand}` |

The integer-vs-float region storage is the **second** known v2/legacy
divergence (the first being the fully-spans-region behavior, §4
above). For a fully-contained fragment, $\sum_\text{regions overlapped by any block}
w = 1$ — but the regions overlapped by a contained fragment are *all
the same single region*, so the legacy float adds $w_1 + w_2 + \ldots = 1$
to one slot. v2's `+= 1` is therefore exactly equivalent. **No
divergence in practice** for contained fragments.

---

## 9. Found bugs

None. The legacy implementation is well-designed and matches its
documented invariants. The accumulator/payload sanity checks in
`scan_payload.py` lines 195–350 confirm mass conservation, channel
totals, signature totals, and FL pool consistency on every payload
produced.

---

## 10. Decisions locked for Phase 2 reference implementation

1. **Per-fragment normalization.** $L = \sum_k \ell_k$ computed once
   per fragment; each block-region overlap deposits
   $\text{overlap\_bp} / L$ to the relevant slot.
2. **Per-fragment splice flag.** Inherit the upstream resolver's
   single per-fragment classification. Update
   [`00_design.md`](00_design.md) §4.3 to drop the
   "per_junction_spliced[K-1]" framing; replace with a single bool.
3. **Per-fragment strand.** Inherit; ambiguous fragments dropped.
4. **Region storage.** `uint32×4` per region (contained channels
   only). Contained fragment increments by exactly +1 regardless of
   block count.
5. **Boundary storage.** `float32×8` (mass_left × 4 + mass_right × 4)
   + `uint32×4` flux per boundary. Per-fragment-junction-event
   deposits per §4.3 of [`00_design.md`](00_design.md).
6. **FL bins deferred indefinitely.** The calibration substrate
   does **not** accumulate per-region or per-boundary FL histograms.
   Spliced fragments are FL-ambiguous (compatible with many isoform
   lengths), and the calibrator does not need FL as a likelihood
   signal — its discrimination comes from count (NB), strand
   (Beta-Binomial), and spliced (Poisson) channels. FL continues to
   be consumed downstream by the EM scorer
   ([`scoring.py`](../../src/rigel/scoring.py),
   [`frag_length_model.py`](../../src/rigel/frag_length_model.py))
   per-fragment at full resolution; that path is unaffected.
7. **Fully-spans-region.** Decompose block into per-region slices;
   apply §4.3 rule to consecutive slices. Accept that v2 and legacy
   diverge in this rare case; Phase 5 parity test will exclude or
   custom-compare these boundaries.
8. **Strand-ambiguous, chimeric, multimap fragments.** All dropped
   from calibration substrate per the existing
   `bam_scanner.cpp:1640` gate.

---

## 11. Action items for Phase 1 completion

- [x] Audit memo (this document).
- [ ] Python reference implementation
      (`tests/native/_accumulator_reference.py`).
- [ ] 11 xfail spec tests
      (`tests/native/test_accumulator_spec.py`).
- [ ] Update [`00_design.md`](00_design.md) §4.3 to drop the per-junction splice framing.
- [ ] `pytest tests/ -v` confirms legacy tests still green and v2 tests xfail.

No legacy code changes required.
