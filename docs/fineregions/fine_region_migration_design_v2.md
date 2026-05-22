# Fine-Grained Region Migration Design — v2 (Implementation-Ready)

Date: 2026-05-22
Supersedes: [fine_region_migration_design.md](fine_region_migration_design.md)

## 0. Scope and Outcome

This document is the implementation-ready version of the fine-grained region
migration. A developer should be able to start Phase 0 from this doc alone.

**Outcome.** Replace the persisted 3-type / 4-strand coarse region partition
with a non-overlapping per-reference partition keyed by a 4-bit fine signature
`{intron_pos, intron_neg, exon_pos, exon_neg}`. All current calibration math
continues to function bit-for-bit through a deterministic fine→coarse adapter
("parity bridge"), and downstream phases progressively replace coarse reads
with fine-aware reads.

**Non-goals (this migration).**

- No new calibration math is enabled in the first PR. Parity first.
- No reintroduction of multimapper calibration evidence.
- No change to EM/locus/scoring code paths other than what calibration outputs
  feed them.
- No backward compatibility for old indexes or scan payloads.

---

## 1. Critique of v1

v1 ([fine_region_migration_design.md](fine_region_migration_design.md))
captured the right *direction* but is not implementation-ready. The
substantive issues:

1. **Underestimated current sophistication.** The existing
   `CalibrationAccumulator` already tracks:
   - per-region per-mask integer counts (`per_region_counts[R, 8]`),
   - per-region orientation channels for INTRON contained and EXON contained
     (`{intron,exon_contained}_counts_by_orient[R, 3]`),
   - per-region left/right boundary flux with orientation
     (`u_{left,right}[_by_orient]`),
   - fragment-length histogram stratified by 8-state mask (`fl_hist[8, 1024]`),
   - splicing-anchor-tolerance `K` qualification with boundary predicate
     `frag_start + q <= b && frag_end >= b + q`.

   v1's "replace integer mask counts with a fractional 12-channel array"
   would *delete* the orientation, contained-vs-flux split, and the
   endpoint-based boundary semantics that current densities depend on,
   and would require simultaneously rewriting
   [density_global.py](../../src/rigel/calibration/density_global.py),
   [_regional_exposure.py](../../src/rigel/calibration/_regional_exposure.py),
   and
   [locus_prior.py](../../src/rigel/calibration/locus_prior.py).
   That coupling is the opposite of "seamless".

2. **Conflated two orthogonal changes.** Fine *partition* is a strict
   improvement to the M2 invariant table. Fractional 12-channel *accumulation*
   is a different decision about how to bin per-fragment evidence. v2 keeps
   them separate, and only the partition change is mandatory.

3. **No native ABI spec.** v1 said "pass signatures to native code" without
   specifying the new `RegionIndex::set(...)` signature, the new payload
   dict keys, or the C++ namespace constants.

4. **No interaction with `splicing_anchor_tolerance`.** The current builder
   sets `boundary_flux_left/right` to drive the boundary-flux density's
   eligibility filter, and the accumulator uses `q = max(K, 1)` for both
   per-block overlap qualification and boundary-flux endpoint predicates.
   For a *fine* partition, a region's left/right boundary semantics
   depend on the neighboring fine signature (e.g. intron→exon vs exon→exon
   on opposite strand vs exon→intergenic). v1 did not address this.

5. **No test plan.** "Add tests" is not enough. v2 enumerates the test matrix,
   the golden update protocol, and the parity gates that must pass before
   each phase is mergeable.

6. **No performance / memory budget.** A fine partition can grow the region
   count by ~1.5–2x on GENCODE-class annotations. Per-worker `R*12*8` bytes
   for a float64[R,12] accumulator becomes nontrivial at high thread counts.
   v2 sets a budget and documents the worst case.

7. **No rollout / risk register.** Since we are explicitly *not* preserving
   compatibility, the cutover protocol must be explicit (index rebuild,
   golden refresh, downstream tool re-runs).

8. **Bit-layout typo unresolved.** v1 noted the TODO mislabel of `b0100`
   but treated it as prose; v2 ships a `signature.py` module whose
   constants and `pack_signature(...)` make the layout unambiguous, plus
   a test that exercises all 16 states.

9. **Region ref-id mapping inconsistency.** v1 noted but did not commit to
   the fix. `_wire_calibration_regions(...)` currently drops regions whose
   `ref_name` is not in the resolver's `ref_name_to_id`. With a
   reference-complete fine partition this is a quiet correctness bug
   because intergenic-only contigs will silently disappear. v2 fixes it.

---

## 2. Design Invariants and Definitions

### 2.1 Canonical signature

The persisted fine-region signature is a `uint8` with the canonical layout:

```text
bit 3 (0x8): intron_pos
bit 2 (0x4): intron_neg
bit 1 (0x2): exon_pos
bit 0 (0x1): exon_neg

signature = (intron_pos << 3) | (intron_neg << 2)
          | (exon_pos   << 1) | (exon_neg   << 0)
```

The TODO's prose example `b0100 = intron_pos` is a swapped-bit typo. The
coarse-state truth table in the TODO is consistent with the canonical
layout above and is the authoritative reference (reproduced in §2.3).

A `signature.py` module is the single source of truth for the layout in
Python; in C++ the same constants live in
`src/rigel/native/calibration/region_signature.h`.

### 2.2 Derived coarse fields

Coarse fields are **derived, not persisted as primary semantics**. They
exist in the on-disk schema solely so existing readers continue to work
during the parity bridge, and so `pandas` queries are ergonomic.

```python
has_pos     = intron_pos | exon_pos
has_neg     = intron_neg | exon_neg
coarse_type = EXON       if (exon_pos | exon_neg)   else (
              INTRON     if (intron_pos | intron_neg) else
              INTERGENIC)
coarse_strand = NONE | POS if has_pos | NEG if has_neg   # POS|NEG = AMBIG=3
coarse_ambig  = has_pos & has_neg
```

`coarse_type`'s exon-wins-over-intron rule preserves the current
[regions.py](../../src/rigel/calibration/regions.py) semantics
(emit_regions exon override). This is required for parity.

### 2.3 Coarse derivation table

| signature | bits                          | coarse_type | coarse_strand | ambig |
| --------: | ----------------------------- | ----------- | ------------- | :---: |
| `0x0`     | intergenic                    | INTERGENIC  | NONE          |       |
| `0x1`     | exon_neg                      | EXON        | NEG           |       |
| `0x2`     | exon_pos                      | EXON        | POS           |       |
| `0x3`     | exon_pos+exon_neg             | EXON        | AMBIG         |   y   |
| `0x4`     | intron_neg                    | INTRON      | NEG           |       |
| `0x5`     | intron_neg+exon_neg           | EXON        | NEG           |       |
| `0x6`     | intron_neg+exon_pos           | EXON        | AMBIG         |   y   |
| `0x7`     | intron_neg+exon_pos+exon_neg  | EXON        | AMBIG         |   y   |
| `0x8`     | intron_pos                    | INTRON      | POS           |       |
| `0x9`     | intron_pos+exon_neg           | EXON        | AMBIG         |   y   |
| `0xA`     | intron_pos+exon_pos           | EXON        | POS           |       |
| `0xB`     | intron_pos+exon_pos+exon_neg  | EXON        | AMBIG         |   y   |
| `0xC`     | intron_pos+intron_neg         | INTRON      | AMBIG         |   y   |
| `0xD`     | intron_pos+intron_neg+exon_neg| EXON        | AMBIG         |   y   |
| `0xE`     | intron_pos+intron_neg+exon_pos| EXON        | AMBIG         |   y   |
| `0xF`     | all four flags                | EXON        | AMBIG         |   y   |

### 2.4 `boundary_flux_left/right` policy

For coarse parity, `boundary_flux_{left,right}` is recomputed as a *derived
property of the fine partition*:

- `boundary_flux_left[i]  = (coarse_type[i] == EXON) AND
                             (i > 0)                AND
                             (ref_name[i-1] == ref_name[i])`
- `boundary_flux_right[i] = (coarse_type[i] == EXON) AND
                             (i+1 < R)              AND
                             (ref_name[i+1] == ref_name[i])`

This is the strict fine analog of the current builder: an EXON region's
edge is flux-eligible iff there is a same-reference neighbor on that side.
Note that on a fine partition this *does* include adjacent EXON–EXON edges
(e.g. `b0010 → b0011` opposite-strand entry); these were collapsed into a
single EXON run under the coarse builder. We accept this delta and address
it in Phase 4 by deriving boundary eligibility from adjacent fine
signatures rather than this coarse-parity rule.

### 2.5 Splicing-anchor-tolerance interaction

`splicing_anchor_tolerance K` is unchanged in semantics. The fine
partition does not alter:

- the per-block overlap qualification `overlap_bp >= q = max(K, 1)`,
- the boundary-flux endpoint predicate
  `frag_start + q <= b && frag_end >= b + q`,
- the `n_below_tolerance` accountability.

These remain in `CalibrationAccumulator::observe`. The fine signature only
changes (a) which regions exist on disk, and (b) what neighbor relationship
is encoded by `boundary_flux_{left,right}`.

### 2.6 No legacy load path

`INDEX_FORMAT_VERSION` is bumped from 3 to **4**. On encountering a
version < 4 (or a `regions.feather` lacking the `signature` column), the
loader raises `ValueError` with `Rebuild the index (rigel index ...)`.
There is no migration shim.

---

## 3. Phase 0 — Cleanup, Signature Helpers, Ref-id Fix

**Goal.** Land low-risk infrastructure that the rest of the migration
builds on, without changing any persisted formats or downstream math.

### 3.1 Deliverables

- `src/rigel/calibration/signature.py` (new) — single source of truth for
  the 4-bit layout and derivation helpers.
- `src/rigel/native/calibration/region_signature.h` (new) — mirror C++
  constants. Build-fail if Python/C++ disagree (compile-time
  `static_assert` via a small Python-generated header check, or a runtime
  test in `tests/test_signature_layout_native.py`).
- `tests/test_signature.py` (new) — exhaustive 16-state derivation test.
- `_wire_calibration_regions(...)` in
  [pipeline.py](../../src/rigel/pipeline.py): unify region ref-id mapping
  to `index.ref_name_to_id` (BAM tid space), with explicit handling of
  references that are present in the index but absent from the BAM header.
- Remove dead code:
  - `set_regions_legacy(...)` if present in
    `src/rigel/native/bam_scanner.cpp` and its nb binding,
  - any optional/missing-region fallback in `scan_and_buffer(...)`.

### 3.2 `signature.py` skeleton

```python
# src/rigel/calibration/signature.py
"""Canonical fine-region 4-bit signature layout and derivations."""
from __future__ import annotations

import numpy as np

# Bit positions (canonical layout — DO NOT REORDER).
BIT_INTRON_POS = 3
BIT_INTRON_NEG = 2
BIT_EXON_POS   = 1
BIT_EXON_NEG   = 0

MASK_INTRON_POS = 1 << BIT_INTRON_POS
MASK_INTRON_NEG = 1 << BIT_INTRON_NEG
MASK_EXON_POS   = 1 << BIT_EXON_POS
MASK_EXON_NEG   = 1 << BIT_EXON_NEG

N_SIGNATURES = 16

# Coarse derivation enums (must match the values used by
# rigel.calibration.regions.RegionType / RegionStrand).
COARSE_INTERGENIC = 0
COARSE_INTRON     = 1
COARSE_EXON       = 2

COARSE_STRAND_NONE  = 0
COARSE_STRAND_POS   = 1
COARSE_STRAND_NEG   = 2
COARSE_STRAND_AMBIG = 3   # POS | NEG


def pack_signature(
    intron_pos: bool | np.ndarray,
    intron_neg: bool | np.ndarray,
    exon_pos:   bool | np.ndarray,
    exon_neg:   bool | np.ndarray,
) -> np.ndarray:
    """Vectorized packing into the canonical uint8 layout."""
    ip = np.asarray(intron_pos, dtype=np.uint8)
    in_ = np.asarray(intron_neg, dtype=np.uint8)
    ep = np.asarray(exon_pos,   dtype=np.uint8)
    en = np.asarray(exon_neg,   dtype=np.uint8)
    return ((ip << BIT_INTRON_POS) | (in_ << BIT_INTRON_NEG)
          | (ep << BIT_EXON_POS)   | (en  << BIT_EXON_NEG)).astype(np.uint8)


def unpack_signature(sig: np.ndarray) -> dict[str, np.ndarray]:
    s = np.asarray(sig, dtype=np.uint8)
    return {
        "intron_pos": ((s >> BIT_INTRON_POS) & 1).astype(bool),
        "intron_neg": ((s >> BIT_INTRON_NEG) & 1).astype(bool),
        "exon_pos":   ((s >> BIT_EXON_POS  ) & 1).astype(bool),
        "exon_neg":   ((s >> BIT_EXON_NEG  ) & 1).astype(bool),
    }


def coarse_type_from_signature(sig: np.ndarray) -> np.ndarray:
    s = np.asarray(sig, dtype=np.uint8)
    has_exon   = ((s & (MASK_EXON_POS   | MASK_EXON_NEG))   != 0)
    has_intron = ((s & (MASK_INTRON_POS | MASK_INTRON_NEG)) != 0)
    out = np.full(s.shape, COARSE_INTERGENIC, dtype=np.uint8)
    out[has_intron] = COARSE_INTRON
    out[has_exon]   = COARSE_EXON   # exon wins
    return out


def coarse_strand_from_signature(sig: np.ndarray) -> np.ndarray:
    s = np.asarray(sig, dtype=np.uint8)
    has_pos = ((s & (MASK_INTRON_POS | MASK_EXON_POS)) != 0)
    has_neg = ((s & (MASK_INTRON_NEG | MASK_EXON_NEG)) != 0)
    out = np.full(s.shape, COARSE_STRAND_NONE, dtype=np.uint8)
    out[has_pos] = COARSE_STRAND_POS
    out[has_neg & ~has_pos] = COARSE_STRAND_NEG
    out[has_pos & has_neg]  = COARSE_STRAND_AMBIG
    return out


def coarse_ambig_from_signature(sig: np.ndarray) -> np.ndarray:
    s = np.asarray(sig, dtype=np.uint8)
    has_pos = ((s & (MASK_INTRON_POS | MASK_EXON_POS)) != 0)
    has_neg = ((s & (MASK_INTRON_NEG | MASK_EXON_NEG)) != 0)
    return has_pos & has_neg
```

### 3.3 C++ mirror

```cpp
// src/rigel/native/calibration/region_signature.h
#pragma once
#include <cstdint>

namespace rigel::calibration::sig {
constexpr uint8_t BIT_INTRON_POS = 3;
constexpr uint8_t BIT_INTRON_NEG = 2;
constexpr uint8_t BIT_EXON_POS   = 1;
constexpr uint8_t BIT_EXON_NEG   = 0;

constexpr uint8_t MASK_INTRON_POS = 1u << BIT_INTRON_POS;
constexpr uint8_t MASK_INTRON_NEG = 1u << BIT_INTRON_NEG;
constexpr uint8_t MASK_EXON_POS   = 1u << BIT_EXON_POS;
constexpr uint8_t MASK_EXON_NEG   = 1u << BIT_EXON_NEG;

constexpr int     N_SIGNATURES    = 16;

inline constexpr uint8_t coarse_type_from(uint8_t s) {
    const bool exon   = (s & (MASK_EXON_POS   | MASK_EXON_NEG))   != 0;
    const bool intron = (s & (MASK_INTRON_POS | MASK_INTRON_NEG)) != 0;
    // 0 = INTERGENIC, 1 = INTRON, 2 = EXON  (matches Python COARSE_*)
    return exon ? 2u : (intron ? 1u : 0u);
}
inline constexpr uint8_t coarse_strand_from(uint8_t s) {
    const bool pos = (s & (MASK_INTRON_POS | MASK_EXON_POS)) != 0;
    const bool neg = (s & (MASK_INTRON_NEG | MASK_EXON_NEG)) != 0;
    return static_cast<uint8_t>((pos ? 1u : 0u) | (neg ? 2u : 0u));
}
}  // namespace rigel::calibration::sig
```

### 3.4 Ref-id mapping fix

In `pipeline._wire_calibration_regions(...)`, replace any region ref-name
→ id lookup that uses the resolver's transcript-ref space with
`index.ref_name_to_id`. Add a debug-log line:

```python
n_dropped = (region_df["ref_name"].map(index.ref_name_to_id).isna()).sum()
if n_dropped:
    log.warning("calibration: dropping %d region rows (%s) — refs not in BAM header",
                n_dropped, sorted_unique_dropped_refs)
```

References in the index but missing from the BAM header are dropped with a
warning. References in the BAM header but missing from the index never had
regions to begin with.

### 3.5 Tests

- `tests/test_signature.py`:
  - exhaustively check all 16 states against the §2.3 table for `pack`,
    `unpack`, `coarse_type_from_signature`, `coarse_strand_from_signature`,
    `coarse_ambig_from_signature`;
  - randomized round-trip `pack(unpack(s)) == s` for `s in 0..15`.
- `tests/test_signature_layout_native.py`:
  - import the native module and compare the C++ `sig::*` constants and
    `coarse_type_from`/`coarse_strand_from` outputs against the Python
    table over all 16 states.
- `tests/test_pipeline_wiring.py` (extend):
  - assert that when the index has an intergenic-only reference, that
    reference's regions still reach the native scanner.

### 3.6 Acceptance gates

- All existing tests pass unchanged.
- New tests pass.
- No persisted-format change. No `INDEX_FORMAT_VERSION` bump in Phase 0.

---

## 4. Phase 1 — Fine Index Builder + Parity Bridge

**Goal.** Build the fine region table and persist it. Continue emitting
the coarse columns the rest of the code currently consumes, all
*derived from `signature`*. Native code and scan payload schema are
**unchanged in this phase**.

### 4.1 Deliverables

- `build_fine_region_table(transcripts, ref_lengths) -> pd.DataFrame` in
  [regions.py](../../src/rigel/calibration/regions.py).
- New on-disk `regions.feather` schema (see §4.3).
- `INDEX_FORMAT_VERSION = 4`. Hard-fail loaders on version < 4.
- `RegionRecord` retired as a tuple-of-arrays emitter; replaced by a
  vectorized builder. Old `emit_regions(...)` and per-base / event-sweep
  EXON/INTRON code removed.
- `validate_against_ref_lengths(...)` extended for new fields.
- `RegionArrays.from_region_df(...)` and `PayloadArrays.from_payload(...)`
  carry the new fields, but downstream code (density_global,
  _regional_exposure, locus_prior) reads only the *derived coarse*
  columns. Bit-for-bit parity required.

### 4.2 Builder algorithm

Per reference:

1. Collect events from every non-synthetic transcript:
   - exon: `(start, exon_<strand>, +1)`, `(end, exon_<strand>, -1)`.
   - intron (between consecutive exons of the same transcript):
     `(start, intron_<strand>, +1)`, `(end, intron_<strand>, -1)`.
2. Sort events by position. Sweep left-to-right with four running
   counters `(ip, in_, ep, en)`.
3. Between consecutive distinct event positions `[p_k, p_{k+1})`, the
   active 4-tuple is constant. Pack it as a uint8 signature.
4. Emit a region only when the signature changes from the previous
   emitted segment (adjacent-merge semantics, identical to the old
   builder's behavior).
5. Reference boundaries `0` and `ref_length` are implicit start/end
   events.
6. Empty references (no events) emit exactly one `signature == 0`
   region covering `[0, ref_length)`.

Per-reference cost: `O(E log E)` where `E` is the number of transcript
exon+intron interval edges on that reference. No `O(L)` allocation.

### 4.3 New `regions.feather` schema

Columns (locked dtypes):

| Column                | dtype     | Notes                                            |
| --------------------- | --------- | ------------------------------------------------ |
| `region_id`           | int64     | Globally sequential = row index.                 |
| `ref_name`            | string    | FASTA reference name.                            |
| `start`               | int64     | 0-based inclusive.                               |
| `end`                 | int64     | 0-based exclusive.                               |
| `length`              | int64     | `end - start`. Persisted for diagnostics.        |
| `signature`           | uint8     | Canonical 4-bit fine state.                      |
| `intron_pos`          | bool      | Explicit flag (== `signature & 0x8 != 0`).       |
| `intron_neg`          | bool      | Explicit flag (== `signature & 0x4 != 0`).       |
| `exon_pos`            | bool      | Explicit flag (== `signature & 0x2 != 0`).       |
| `exon_neg`            | bool      | Explicit flag (== `signature & 0x1 != 0`).       |
| `type`                | uint8     | *Derived* coarse type; `RegionType` values.      |
| `strand`              | uint8     | *Derived* coarse strand; `RegionStrand` values.  |
| `boundary_flux_left`  | bool      | §2.4 fine-derived rule.                          |
| `boundary_flux_right` | bool      | §2.4 fine-derived rule.                          |

Removed from the schema:
`tx_pos_bp`, `tx_neg_bp`, `exon_pos_bp`, `exon_neg_bp`. These were never
consumed past the v6 plan; no current downstream code reads them
(verified via grep). They are not replaced.

### 4.4 Validation

Extend `validate_against_ref_lengths(...)`:

- `signature` values are in `[0, 15]`.
- `intron_pos/neg`, `exon_pos/neg` boolean flags are *exactly*
  `unpack_signature(signature)`.
- `type` == `coarse_type_from_signature(signature)`.
- `strand` == `coarse_strand_from_signature(signature)`.
- `boundary_flux_left[i]`, `boundary_flux_right[i]` match §2.4.
- All existing tiling checks (coverage, sorted, contiguous, region_id =
  row index).

### 4.5 Loader cutover

```python
def load_regions(path):
    df = pd.read_feather(path)
    if "signature" not in df.columns:
        raise ValueError(
            f"regions.feather at {path} is missing 'signature' (pre-v4 "
            f"index). Rebuild the index (rigel index --fasta ... --gtf ...)."
        )
    ...
```

`index.py`'s metadata loader checks `INDEX_FORMAT_VERSION >= 4` and
raises the same actionable error otherwise.

### 4.6 `RegionArrays` / `PayloadArrays`

[_arrays.py](../../src/rigel/calibration/_arrays.py):

- `RegionArrays` gains a `signature: np.ndarray  # uint8, (R,)` field.
  All other fields keep current dtypes. `type` and `strand` continue to
  be carried (derived).
- `PayloadArrays` unchanged in Phase 1.

### 4.7 Tests

`tests/test_regions.py` — replace existing tests with the fine-builder
matrix. Each is a synthetic transcript layout assertion on the emitted
DataFrame:

| Case                                                              | Signatures expected (in order)                |
| ----------------------------------------------------------------- | --------------------------------------------- |
| Empty reference                                                   | `[0x0]`                                       |
| Single (+) transcript with 2 exons + 1 intron, intergenic flanks  | `0x0, 0x2, 0x8, 0x2, 0x0`                     |
| Two (+) transcripts with overlapping exons                        | merges identical adjacent signatures          |
| (+) and (-) transcripts overlapping on exons                      | includes `0x3` (`exon_pos+exon_neg`)          |
| (+) exon overlapping (-) intron                                   | includes `0x6` (`intron_neg+exon_pos`)        |
| (+) exon overlapping (+) intron of another transcript             | includes `0xA` (`intron_pos+exon_pos`)        |
| Single-exon transcript                                            | `0x0, 0x{2 or 1}, 0x0`                        |
| Read-through transcript spanning two genes                        | no signature collapse where bits differ       |
| Multi-reference                                                   | `region_id` is globally sequential            |
| Intergenic-only reference                                         | single `0x0` row, propagated through pipeline |

`tests/test_region_parity.py` (new) — for each scenario above, assert:

- `signature == pack_signature(intron_pos, intron_neg, exon_pos, exon_neg)`,
- `type == coarse_type_from_signature(signature)`,
- `strand == coarse_strand_from_signature(signature)`,
- `boundary_flux_*` matches §2.4.

`tests/test_index_format_version.py`: building an index writes
`INDEX_FORMAT_VERSION == 4`; loading a hand-written feather with
version 3 raises `ValueError`.

### 4.8 Golden refresh

All golden outputs in `tests/golden/` are recomputed once in this PR
with `pytest tests/ --update-golden`. We expect deltas only where:

- region counts change (more rows on indexes that have overlapping
  transcripts on opposite strands, or read-through structures), which
  is observable in `loci_df.feather` only through region-driven
  exposure terms;
- nothing else, because the calibration math reads derived coarse
  columns identically.

Treat any non-trivial change to `transcript_df`, `gene_df`, or
`scalars.json` outside known explanations as a regression and root-cause
before merging.

### 4.9 Acceptance gates

- New `regions.feather` schema validates on every test fixture.
- Coarse-parity assertion: for every fixture used in calibration tests
  (`tests/test_density_global.py`, `tests/test_regional_exposure.py`,
  `tests/test_locus_prior.py`), the *derived* coarse columns
  (`type`, `strand`, `boundary_flux_*`) produce bit-identical
  `CalibrationScanPayload`-derived densities to the v3-builder
  baseline.
  - Mechanism: a temporary in-tree `tests/_v3_builder_snapshot/`
    captures the v3 builder output on every fixture before PR1 lands;
    the parity test reloads those snapshots and compares row-aligned
    derived columns.
- All existing tests pass with golden refresh.
- `INDEX_FORMAT_VERSION == 4` enforced.

---

## 5. Phase 2 — Native Region Index Carries Signature

**Goal.** The native `RegionIndex` learns about `signature` *additively*.
The accumulator does not change behavior, but it can now query
`regions.signature(rid)` if it wants to. We keep the integer
`per_region_counts[R, 8]`, orientation channels, and boundary flux
exactly as today.

### 5.1 Deliverables

- `RegionIndex::set(...)` adds a `signatures` parameter and exposes
  `signature(rid)`. `type_masks` and `strands` continue to be passed
  in (derived in Python from `signature`) so the rest of the C++ code
  is byte-identical.
- nanobind binding updated; Python `set_regions(...)` API gains a
  `signatures` argument.
- `pipeline._wire_calibration_regions(...)` derives and passes
  `signatures` from `region_df["signature"]` (no Python-side derivation
  of `type_masks`/`strands` — read them from the new derived columns
  to avoid drift).

### 5.2 New native `set` signature

```cpp
void set(const int32_t* ref_ids,
         const int64_t* starts,
         const int64_t* ends,
         const uint8_t* signatures,   // NEW — required
         const uint8_t* type_masks,   // kept for now; trivially derivable
         const uint8_t* strands,      // kept for now; trivially derivable
         int64_t n_regions,
         int32_t n_refs);

inline uint8_t signature(int32_t rid) const { return signatures_[rid]; }
```

A future PR may remove `type_masks` and `strands` parameters and derive
them inside `set(...)` via `sig::coarse_type_from`/`coarse_strand_from`,
once Phase 4 has eliminated direct coarse-mask reads on the C++ side.
We do not do this now to minimize Phase 2's diff surface.

### 5.3 Wiring (Python side)

In `pipeline._wire_calibration_regions(...)`:

```python
sig_arr = region_df["signature"].to_numpy(dtype=np.uint8, copy=False)
type_arr = region_df["type"].to_numpy(dtype=np.uint8, copy=False)
strand_arr = region_df["strand"].to_numpy(dtype=np.uint8, copy=False)
# Sanity: derived columns must match.
assert np.array_equal(type_arr, coarse_type_from_signature(sig_arr))
assert np.array_equal(strand_arr, coarse_strand_from_signature(sig_arr))

# 3-bit type mask the accumulator currently consumes.
# COARSE_*  -> 3-bit obs_mask bit
#  INTERGENIC -> 0b100
#  INTRON     -> 0b010
#  EXON       -> 0b001
type_mask_arr = np.choose(type_arr, [0b100, 0b010, 0b001]).astype(np.uint8)

native.set_regions(ref_ids, starts, ends, sig_arr, type_mask_arr, strand_arr, ...)
```

### 5.4 Tests

- `tests/test_region_index_native.py` (extend): construct a fixture
  with all 16 signatures and assert `RegionIndex.signature(rid)` returns
  the expected value.
- Verify that the Python `assert` on derived columns matches before the
  call into native (i.e. parity check is a runtime invariant, not just a
  Python-side test).

### 5.5 Acceptance gates

- All existing scanner tests pass unchanged.
- `pytest tests/` is green; goldens unchanged (Phase 1 already refreshed
  them).

---

## 6. Phase 3 — Fine Payload Channels (Additive)

**Goal.** Add a new per-region count block keyed by fine signature
without removing any existing channel. This is the first calibration
math change: downstream code may *prefer* fine reads but the parity
bridge from Phase 1 stays in place.

### 6.1 Deliverables

- `CalibrationPayload` gains:
  - `per_region_counts_by_sig: int64[R, 16]` — per-region per-signature
    counts of qualified hits (the same `qualified_hits_` loop in the
    accumulator). Indexing: `rid * 16 + signature(rid)`. This is
    sparse-by-construction (only `signature(rid)` is incremented for
    each `rid`), so the only useful information it adds over the
    existing 3-bit mask is the *cross-region distribution* of qualified
    fragments — see §6.3 for why we keep it anyway.
  - `per_sig_global_counts: int64[16]` — global per-fragment fine
    `obs_sig`, where `obs_sig = OR over qualified hits of signature(rid)`.
    Replaces `global_counts[8]` as the primary fragment-level
    classification; the old 8-state counter is kept and derived from
    `per_sig_global_counts` via `coarse_type_from_signature` (sum over
    fine states that map to each coarse type) so legacy readers continue.
  - `fl_hist_by_sig: int64[16, 1024]` — fragment-length histogram
    stratified by `obs_sig`. The old `fl_hist[8, 1024]` is kept and
    derived analogously.
- Update Python `CalibrationScanPayload` to validate and expose the new
  fields; new fields are *opt-in* readers (`obs_sig_*`), old fields
  remain canonical.
- No change to orientation channels or boundary flux in this phase.

### 6.2 Accumulator changes

`observe(...)` keeps the existing logic. After computing `qualified_hits_`,
also compute:

```cpp
uint8_t obs_sig = 0;  // OR of fine signatures over qualified hits
for (int32_t rid : qualified_hits_) {
    obs_sig |= regions.signature(rid);
}
// Existing path: obs_mask, global_counts, fl_hist, per_region_counts...
// New: per-fine-state aggregates.
payload_.per_sig_global_counts[obs_sig]++;
payload_.fl_hist_by_sig[obs_sig * 1024 + fl_idx]++;
for (int32_t rid : qualified_hits_) {
    payload_.per_region_counts_by_sig[rid * 16 + regions.signature(rid)]++;
}
```

Cost: one extra `OR` per qualified hit and one extra 16x1024 histogram
increment per observed fragment. The 16x1024 histogram is 128 KiB per
worker — fits comfortably in L2.

### 6.3 Why keep `per_region_counts_by_sig` if it is sparse-by-rid?

A multi-block fragment whose blocks fall in regions of *different*
signatures will increment multiple slots in `per_region_counts_by_sig`
(one per (rid, signature(rid)) pair). The *row-summed* per-fine-state
distribution then tells calibration "how many qualified hits across
the genome had fine state s", which is not derivable from
`per_sig_global_counts[s]` alone (the latter is per-fragment OR'd).
This is the key new statistic that Phase 4's per-signature density
estimator will consume.

### 6.4 Payload schema validation (Python)

`CalibrationScanPayload.from_scan_dict`:

- Adds shape checks for the new arrays.
- Adds derived-equivalence assertion at construction time:
  - `global_counts == coarse_collapse(per_sig_global_counts)`
  - `fl_hist == coarse_collapse(fl_hist_by_sig)`
  where `coarse_collapse(x)` sums fine-state slots whose
  `coarse_type_from_signature` matches each of the 3 coarse bits and
  packs into the 8-state `obs_mask` (treating `obs_mask` as the OR of
  bits, which means a fine state with both `exon` and `intron` flags
  contributes to *both* coarse bits — match the current accumulator
  behavior; see §6.5).

### 6.5 Coarse-fine equivalence

Subtlety: the current 8-state `obs_mask` is the OR of per-region
3-bit `type_mask`s. A region with `coarse_type == EXON` contributes
`0b001`. A region with `coarse_type == INTRON` contributes `0b010`.
A fragment that touches both an EXON region (signature `0x5`,
intron_neg+exon_neg) and an INTRON region (signature `0x4`,
intron_neg) has `obs_mask = 0b011`. The same fragment under fine
accounting has `obs_sig = 0x5 | 0x4 = 0x5`. The coarse collapse of
`obs_sig = 0x5` is `coarse_type(0x5) = EXON` → `0b001`. These
*disagree*. We resolve this by defining `coarse_collapse(obs_sig)`
not as `coarse_type_from_signature(obs_sig)` but as the bitwise OR of
per-fine-state coarse bits computed pointwise:

```python
# Build a length-16 lookup table.
COARSE_OBS_MASK_FROM_SIG = np.zeros(16, dtype=np.uint8)
for s in range(16):
    bits = 0
    if s & (MASK_EXON_POS   | MASK_EXON_NEG):   bits |= 0b001
    if s & (MASK_INTRON_POS | MASK_INTRON_NEG): bits |= 0b010
    if s == 0:                                   bits |= 0b100
    COARSE_OBS_MASK_FROM_SIG[s] = bits
```

A fragment's coarse `obs_mask` is then the OR of
`COARSE_OBS_MASK_FROM_SIG[regions.signature(rid)]` across qualified
hits, which the accumulator already computes via the existing
`type_mask(rid)` path (since `type_mask` was derived from
`signature` upstream, agreement is automatic). The derived-equivalence
assertion in §6.4 therefore reads:

```python
assert global_counts == fold_obs_mask(per_sig_global_counts)
```

where `fold_obs_mask` accumulates each fragment's per-fine state into
a coarse 8-state bin using `COARSE_OBS_MASK_FROM_SIG` *as if* the
fragment's qualified hits were drawn from a single fine state — which
is only an approximation. **In practice we do not assert pointwise
equality**; we assert that the *totals* match (`per_sig_global_counts.sum()
== global_counts.sum() == n_observed`) and document the deliberate
non-equivalence on multi-region fragments. The fine block is *new
evidence*, not a redundant view of the coarse block.

### 6.6 Tests

- `tests/test_accumulator_fine.py` (new):
  - Single fragment in pure-intergenic region: `per_sig_global_counts[0] = 1`.
  - Single fragment in exon_pos region (`0x2`):
    `per_sig_global_counts[2] = 1`, `per_region_counts_by_sig[rid, 2] = 1`.
  - Two-block spliced fragment crossing `0x2 → 0x8 → 0x2`: assert
    `obs_sig == 0xA`, `per_region_counts_by_sig` has 3 incremented
    slots, FL histogram updates `[0xA, fl_idx]`.
- `tests/test_scan_payload_fine.py`: shape/dtype checks; total-count
  identity assertions.
- Re-run all existing calibration tests; they must pass unchanged
  because they read coarse channels only.

### 6.7 Acceptance gates

- Existing tests pass unchanged.
- Total-count identities hold on every synthetic fixture.
- No golden updates required (no consumer of the new arrays yet).

---

## 7. Phase 4 — Calibration Consumers Move to Fine Reads

**Goal.** Replace direct reads of the coarse 8-state `per_region_counts`
and 4-state `RegionStrand` with reads of `per_region_counts_by_sig` and
`signature`. The parity bridge from Phase 1 is *removed* at the end of
this phase: `type`, `strand`, `boundary_flux_*` become diagnostic-only
columns or are dropped from `RegionArrays`.

### 7.1 Deliverables

- `density_global.compute_global_densities(...)` rewritten to compute
  numerators from `per_region_counts_by_sig` rather than
  `per_region_counts[:, MASK_*]`. Boundary-flux density continues to
  use `u_{left,right}` (the integer channels), but the eligibility
  filter `boundary_flux_{left,right}` is replaced by a per-side
  *adjacent-signature* lookup table:

  ```python
  # left side is eligible iff the same-ref neighbor on the left has
  # an intron flag matching at least one of this region's exon flags
  # (i.e. it is an actual exon-intron transition on a shared strand).
  ```

- `_regional_exposure.RegionalGdnaExposure.build(...)` similarly
  switches `type_arr` reads to `signature` queries, and uses the
  fine per-signature density when present.
- `locus_prior.assemble_priors(...)` and friends switch reads.
- `RegionArrays` loses its `type`, `strand`, `bf_left`, `bf_right`
  fields once nothing reads them (verified by removing the fields and
  running `pytest`).

### 7.2 Fine-aware density classes

Define 5+ fine-aware density classes:

| Class                    | Eligibility                                                      |
| ------------------------ | ---------------------------------------------------------------- |
| `INTERGENIC`             | `signature == 0x0`                                               |
| `INTRON-PURE`            | `signature in {0x4, 0x8, 0xC}` (intron only, any strand combo)   |
| `EXON-CONTAINED-SOLO`    | `signature in {0x1, 0x2}` (single-strand exon, no overlap)       |
| `EXON-CONTAINED-AMBIG`   | `signature in {0x3, 0x6, 0x9, 0xA}` (mixed exon/intron overlap)  |
| `EXON-INTRON-BOUNDARY`   | qualified left/right boundaries where adjacent signature differs |

Each class derives its own `rho`, `kappa`, and `eff_length_bp`
exactly as today, but from the fine numerator/denominator.

### 7.3 Coarse-parity adapter retired

Once §7.1 lands, `type`, `strand`, `bf_left`, `bf_right`:

- removed from `RegionArrays`,
- moved to "diagnostic only" columns in `regions.feather` (still
  emitted for `pandas` ergonomics) OR dropped entirely — TBD by PR
  reviewers.

### 7.4 Tests

- All `tests/test_density_global.py` cases rewritten to use signature
  fixtures rather than `RegionType` / `RegionStrand` literals.
- `tests/test_regional_exposure.py` rewritten similarly.
- `tests/test_locus_prior.py` extended for read-through and
  opposite-strand-overlap fixtures.
- Golden refresh expected and reviewed: this is the first phase that
  *should* shift downstream numbers materially. Each material delta
  must be explained in the PR description.

### 7.5 Acceptance gates

- Synthetic locus sweep (`scripts/sim/locus_sweep.py`) against
  `scripts/benchmark/configs/locus_simple_baseline.yaml` shows no
  regression in mRNA/nRNA/gDNA relative error.
- Full-scale benchmarks on Armis2
  (`scripts/benchmarking/configs/default.yaml`) over the 13 available
  conditions: no condition regresses mRNA error by more than 5 %
  relative; gDNA error improves on at least the
  `gdna_high_ss_*_nrna_none` conditions (where the fine partition is
  expected to help).
- Runtime: BAM scan wall-clock within 10 % of Phase 3 on the largest
  benchmark BAM.

---

## 8. Phase 5 — Validation, Benchmarks, and Documentation

### 8.1 Test matrix (all phases, run on every PR)

- Unit:
  - `test_signature.py`, `test_regions.py`, `test_region_parity.py`,
    `test_signature_layout_native.py`,
    `test_region_index_native.py`, `test_accumulator_fine.py`,
    `test_scan_payload_fine.py`,
    `test_density_global.py`, `test_regional_exposure.py`,
    `test_locus_prior.py`.
- Integration:
  - `test_cli.py` (`rigel index`, `rigel quant` smoke).
  - `test_golden_output.py` against refreshed goldens.
- Synthetic benchmark:
  - `scripts/sim/locus_sweep.py` on `locus_simple_baseline.yaml`
    diffed against
    `scripts/benchmark/golden/locus_simple_baseline/`.
- Full-scale (run manually before merging Phase 4):
  - `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml`
  - `python -m scripts.benchmarking analyze -c .../default.yaml -o results/fine_regions_phase4/`.

### 8.2 Golden update protocol

- Phase 0: no goldens touched.
- Phase 1: refresh all goldens in a *single commit*, with the commit
  message listing every fixture whose region count changed.
- Phases 2, 3: no goldens touched (regression evidence).
- Phase 4: refresh goldens after acceptance gates pass; document each
  material delta in the PR.

Always use `pytest tests/ --update-golden`; never hand-edit golden
files.

### 8.3 Cross-cutting concerns

**Performance budget.** Region count growth on GENCODE-class
annotations is expected to be ~1.3–1.8x (estimated from overlap rates
in `tests/test_regions.py` synthetic fixtures and confirmed against
the current `regions.feather` of a benchmark run). The native
accumulator adds:
- one OR per qualified hit (negligible),
- one 16×1024 FL histogram (128 KiB per worker),
- one R×16 fine count vector (128 bytes × R per worker; e.g. 128 MiB
  per worker for 1 M regions — emit a warning when total exceeds
  2 GiB across workers).

**Concurrency.** Per-worker accumulators with deterministic merge
order (already implemented). Float channels are not introduced in
this migration, so Kahan summation is not required.

**Observability.** Emit a `regions_summary.json` next to
`regions.feather` at index-build time:

```json
{
  "index_format_version": 4,
  "n_regions": 4731283,
  "per_signature_count": {"0x0": 1234567, "0x1": 89012, ...},
  "per_coarse_type": {"INTERGENIC": ..., "INTRON": ..., "EXON": ...},
  "per_reference_n_regions": {"chr1": ..., ...}
}
```

Add a corresponding `regions_observed_summary.json` in scan output
showing per-signature qualified-hit counts.

**Decoy / contig handling.** References not in the index (and present
in the BAM header) are not touched. References in the index but not in
the BAM header are dropped with a warning. Both behaviors are tested
in `test_pipeline_wiring.py`.

**nRNA span table interaction.** Unaffected. nRNA spans are derived
from per-transcript exon coordinates, not from the calibration region
partition.

### 8.4 Risk register

| Risk                                                              | Likelihood | Severity | Mitigation                                                                                       |
| ----------------------------------------------------------------- | ---------- | -------- | ------------------------------------------------------------------------------------------------ |
| Region count growth breaks memory budget on mega-loci             | medium     | high     | Emit explicit warning at index build when R > 5e6 per ref; add CI fixture with overlapping txs.  |
| Bit-layout typo reintroduced by hand-rolled signature packing     | low        | high     | `signature.py` is the only writer; Python/C++ layout test in CI.                                 |
| `boundary_flux_*` semantic drift in Phase 4 changes EXON-INTRON ρ | medium     | medium   | Compare ρ pre/post on every golden fixture; require PR explanation for material deltas.          |
| Goldens accepted without inspection                               | medium     | high     | Phase 1 PR explicitly lists every changed golden; reviewers must check off.                      |
| `_wire_calibration_regions` ref-id fix regresses BAMs missing refs| low        | medium   | Add `test_pipeline_wiring.py` covering "BAM header subset of index".                             |
| Native ABI drift between Python `signature` constants and C++     | low        | high     | `tests/test_signature_layout_native.py` runs in CI.                                              |

### 8.5 Rollout

Single repository — no feature flag. The migration cuts over via a
sequence of PRs (Phase 0..4). Each PR is independently mergeable and
leaves the codebase green; consumers must rebuild their indexes after
the Phase 1 PR merges.

CHANGELOG entry (added in Phase 1 PR):

```
## Unreleased

### Breaking
- regions.feather format bumped to INDEX_FORMAT_VERSION = 4. Old
  indexes will fail to load with an explicit "Rebuild the index"
  message. The new schema persists a 4-bit fine signature
  (intron_pos, intron_neg, exon_pos, exon_neg) instead of the
  previous 3-type / 4-strand encoding. Rebuild with
  `rigel index --fasta ... --gtf ...`.
```

---

## 9. Open Questions Resolved from v1

1. **Fractional 12-channel accumulator?** Deferred indefinitely. The
   current integer accumulator already provides contained vs flux
   separation via endpoint predicates, and `splicing_anchor_tolerance`
   guarantees per-hit qualification. A fractional rewrite is not
   required to support fine partition.
2. **Persist booleans and signature both?** Yes — both, for
   auditability and so loaders can validate one against the other.
3. **Multimapper calibration evidence?** Not in this migration.
4. **Spliced boundary flux?** Existing rules preserved (unspliced only,
   `flux_eligible = (splice_type == SPLICE_UNSPLICED)`).
5. **Uncertainty into EM prior?** Out of scope here; Phase 4 surfaces
   per-region precision in summaries only.

---

## 10. Implementation Checklist (Phase 0 → 4)

### Phase 0
- [ ] Create `src/rigel/calibration/signature.py` per §3.2.
- [ ] Create `src/rigel/native/calibration/region_signature.h` per §3.3.
- [ ] Add `tests/test_signature.py` (16-state coverage).
- [ ] Add `tests/test_signature_layout_native.py`.
- [ ] Fix `_wire_calibration_regions(...)` ref-id mapping per §3.4.
- [ ] Remove `set_regions_legacy(...)` and missing-region fallback if
      present.
- [ ] CI green; no golden changes.

### Phase 1
- [ ] Implement `build_fine_region_table(...)` per §4.2.
- [ ] Update `regions.feather` schema per §4.3; drop `*_bp` columns.
- [ ] Bump `INDEX_FORMAT_VERSION` to 4 in `index.py`.
- [ ] Update `validate_against_ref_lengths(...)` per §4.4.
- [ ] Update `load_regions(...)` to require `signature` column.
- [ ] Update `RegionArrays.from_region_df(...)` to carry `signature`.
- [ ] Snapshot pre-PR coarse-derived columns into
      `tests/_v3_builder_snapshot/`; add parity test.
- [ ] Refresh all goldens in a single commit; PR description lists
      changed golden files and counts.
- [ ] Add CHANGELOG entry per §8.5.

### Phase 2
- [ ] Add `signatures` parameter to native `RegionIndex::set(...)`.
- [ ] Expose `signature(rid)` on `RegionIndex`.
- [ ] Update nb binding and Python `set_regions(...)` call site.
- [ ] Add `tests/test_region_index_native.py` for `signature(rid)`.

### Phase 3
- [ ] Add `per_region_counts_by_sig`, `per_sig_global_counts`,
      `fl_hist_by_sig` to `CalibrationPayload`.
- [ ] Update `CalibrationAccumulator::observe(...)` and `merge_from(...)`.
- [ ] Update `CalibrationScanPayload.from_scan_dict(...)`.
- [ ] Add `tests/test_accumulator_fine.py` and
      `tests/test_scan_payload_fine.py`.
- [ ] Verify all coarse-consuming tests still pass.

### Phase 4
- [ ] Rewrite `compute_global_densities(...)` to read fine channels.
- [ ] Rewrite `RegionalGdnaExposure.build(...)` to read fine channels.
- [ ] Rewrite `assemble_priors(...)` reads.
- [ ] Drop `type`, `strand`, `bf_left`, `bf_right` from `RegionArrays`.
- [ ] Run synthetic locus sweep and Armis2 benchmarks; document deltas.
- [ ] Refresh goldens and document material changes per fixture.

---

## 11. Files Touched (Reference)

Python:
- [src/rigel/calibration/regions.py](../../src/rigel/calibration/regions.py)
- [src/rigel/calibration/signature.py](../../src/rigel/calibration/signature.py) (new)
- [src/rigel/calibration/scan_payload.py](../../src/rigel/calibration/scan_payload.py)
- [src/rigel/calibration/_arrays.py](../../src/rigel/calibration/_arrays.py)
- [src/rigel/calibration/density_global.py](../../src/rigel/calibration/density_global.py)
- [src/rigel/calibration/_regional_exposure.py](../../src/rigel/calibration/_regional_exposure.py)
- [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py)
- [src/rigel/index.py](../../src/rigel/index.py)
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py)
- [src/rigel/native.py](../../src/rigel/native.py)

Native:
- [src/rigel/native/calibration/region_signature.h](../../src/rigel/native/calibration/region_signature.h) (new)
- [src/rigel/native/calibration/region_index.h](../../src/rigel/native/calibration/region_index.h)
- [src/rigel/native/calibration/accumulator.h](../../src/rigel/native/calibration/accumulator.h)
- [src/rigel/native/calibration/accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp)
- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp)

Tests (new or substantially rewritten):
- [tests/test_signature.py](../../tests/test_signature.py) (new)
- [tests/test_signature_layout_native.py](../../tests/test_signature_layout_native.py) (new)
- [tests/test_regions.py](../../tests/test_regions.py) (rewrite)
- [tests/test_region_parity.py](../../tests/test_region_parity.py) (new)
- [tests/test_region_index_native.py](../../tests/test_region_index_native.py) (extend)
- [tests/test_accumulator_fine.py](../../tests/test_accumulator_fine.py) (new)
- [tests/test_scan_payload_fine.py](../../tests/test_scan_payload_fine.py) (new)
- [tests/test_density_global.py](../../tests/test_density_global.py) (rewrite)
- [tests/test_regional_exposure.py](../../tests/test_regional_exposure.py) (rewrite)
- [tests/test_locus_prior.py](../../tests/test_locus_prior.py) (extend)
- [tests/test_pipeline_wiring.py](../../tests/test_pipeline_wiring.py) (extend)
- [tests/test_index_format_version.py](../../tests/test_index_format_version.py) (new)
