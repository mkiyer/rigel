# Accumulator junction-strand carry + strand-convention unification — design

**Status:** ✅ IMPLEMENTED 2026-06-26. The annotation-derived `boundary_junction_strands` hack is
removed. The splice junction's genomic strand is carried on the boundary node, recorded "for free"
from the motif strand the scanner already computes. Validated: native byte-for-byte spec (incl. new
junction-strand tests) + 232 native/calibration tests green; full suite 1046 non-golden green;
re-scan parity = the hack (over-call +25,624 vs hack +25,968 — the motif catches 2 junctions the
annotation lookup missed). Golden quant outputs shift (the spliced fix) → regenerate with
`--update-golden` once the gDNA-present over-reach (message precision) is settled.

## 0. Why

A boundary node's spliced mass is the **mature-RNA floor** that anchors the adjacent exon's RNA
fraction. By the biology it is **single-strand**: a splice junction's strand is fixed by its GT/AG
motif (axiom 1), and at most one motif-stranded junction occupies a genomic position (axiom 2 —
verified: **0 +/− conflicts** across all junction positions in the benchmark). But the accumulator
stores spliced mass only as motif-**relative** SENSE/ANTISENSE (ch2/ch3), summed into
`mass_spliced`; the junction's **genomic** strand label is discarded. The calibrator then cannot
route the one-sided spliced mass to the correct strand's exon at AMBIG / exon↔exon junctions
(region signatures recover the strand for only 19% of AMBIG-exon spliced sides).

The SENSE/ANTISENSE split is **kept** — it measures library strand-specificity (sense-vs-antisense
RNA fraction; spliced reads are known RNA, so their motif-relative split is the cleanest strandedness
signal). We **add** the junction's genomic strand as a separate per-boundary label.

The scanner already computes the motif strand at deposit (`cr.sj_strand` / the implicit introns'
transcript strand) to orient SENSE/ANTISENSE, so recording it is free — no new computation, just a
field that rides the payload like every other boundary quantity. This deletes the
annotation-lookup plumbing (re-reading `sj.feather`, threading through `calibrate`/`pipeline`) and
handles **novel** junctions (absent from the annotation) for real data.

## 1. Strand-convention unification (do FIRST, independently)

Three names, two value conventions, one of them inconsistent:

| name (file) | NONE | POS | NEG | AMBIG |
|---|---|---|---|---|
| `Strand` (`types.py`, IntEnum, bitwise) | 0 | 1 | 2 | 3 |
| `RegionStrand` (`signature.py`, IntFlag) | 0 | 1 | 2 | 3 |
| `TS_*` (`signature.py`) | 0 | 1 | **−1** | **2** |

`TS_NEG = −1` was justified by "select the sense channel by the sign" — **but no code does sign
arithmetic**; every usage is an equality test against the named constant (`ts == TS_NEG`,
`np.where(ts == TS_NEG, …)`). `strand_class` is never persisted (recomputed from `signature` via
`transcript_strand_class` on load), and there are no literal `== −1` / `== 2`-means-AMBIG
comparisons. So the values are free to change.

**Fix:** redefine `TS_NEG = 2`, `TS_AMBIG = 3` so `TS_* ≡ RegionStrand ≡ Strand` (one value
convention everywhere). Make the source single: `TS_NONE/POS/NEG/AMBIG = int(RegionStrand.NONE/…)`.
`transcript_strand_class` then emits the same values as `coarse_strand_from_signature` (they become
the same map — `transcript_strand_class` stays as the vectorized form `RegionArrays` uses). Delete
the "by the sign" comment. (`Strand` vs `RegionStrand` remain duplicate enums with identical values
— harmless, same convention; consolidating them is an optional later nicety, out of scope here.)

**Blast radius (all equality, all safe):** `bp_solver` (`free_pos/neg`, `spliced_pos/neg`,
`_spliced_faces`), `gdna_strand`, `strand_deconv`, `region_arrays`, `substrate`. Plus any test that
hardcodes a TS value — audit `tests/`. Run the full calibration suite; values changing 1:1 under the
named constants ⇒ green.

## 2. The junction-strand carry (accumulator)

**Storage — one `int8` per boundary** (axiom 2 ⇒ a boundary hosts ≤1 junction, one strand). Values
in the **`Strand`** convention: `NONE=0`, `POS=1`, `NEG=2` (never AMBIG — a junction has a definite
strand). A separate parallel array keeps the 64-byte `Boundary` struct (mass/flux) **byte-for-byte
untouched**:

```cpp
// accumulator.h — Accumulator private state
std::vector<std::int8_t> boundary_junction_strand_;   // size = n_boundaries, 0/1/2 (Strand)
```
with an accessor `const std::int8_t* boundary_junction_strand_data() const`.

**Deposit — record the motif strand on each junction boundary.** Add one parameter:

```cpp
void deposit(const std::int64_t* starts, const std::int64_t* ends, std::size_t n,
             bool spliced, bool primary, std::int32_t strand);   // strand = genomic motif strand
```
In the **crossing** path, whenever a **spliced** fragment deposits its exon-side slice to a boundary
side (channel ≥ 2), set `boundary_junction_strand_[b] = strand`. (Unspliced ⇒ no junction ⇒ field
left 0; contained spliced ⇒ region channel, no boundary.) A multi-junction fragment carries one
transcript strand, applied to all its junction boundaries — consistent. By axiom 2 repeated/merged
sets agree; defensively, set-if-unset and keep the existing value (a later mismatch is a partition
artifact, not biology — optionally count it as a stat).

**merge_from:** `if dst[b]==0: dst[b]=src[b]` per boundary (NONE is the identity; nonzero values
agree by axiom 2).

**Scanner (`deposit_to_accumulator`)** passes the strand it already has:
```cpp
const std::int32_t strand = spliced ? motif_strand : cr.align_strand;  // genomic strand
…deposit(starts, ends, n, spliced, primary, strand);
```
`motif_strand` is `cr.sj_strand` (explicit) or the OR of `cr.implicit_introns` strands (implicit) —
both already computed for the `primary`/orientation logic above the deposit call.

**Python reference (`tests/native/_accumulator_reference.py`)** mirrors all of the above
byte-for-byte: `deposit(..., strand)` sets `boundaries[b].junction_strand` on spliced crossings; a
new `Boundary.junction_strand: int` field; `merge` identical.

## 3. Payload → substrate → consumer (plumbing the field, not re-deriving it)

- **`bam_scanner.cpp build_result`**: emit `cal["boundary_junction_strand"]` — an `int8[B_obj_total]`
  packed in the same per-ref boundary loop that fills `boundary_mass_left` (one extra line reading
  `a.boundary_junction_strand_data()[i]`). Also add a `boundaries_junction_strand` nanobind property
  on `Accumulator` (parallel to `boundaries_mass_left`) so the byte-for-byte test can read it.
- **`scan_payload.py`**: add `boundary_junction_strand: np.ndarray` (int8[B], default zeros for
  backward-compat if a stale payload lacks the key — but a recompiled scanner always emits it).
- **`substrate.py BoundarySubstrate`**: add `junction_strand: np.ndarray` (int8[B]) from the payload,
  aligned with `left_region`/`right_region`.
- **`bp_solver.build_node_geometry`**: read `boundary_substrate.junction_strand` (no new function
  parameter — it rides the substrate). `_spliced_faces` routes a side's `mass_spliced` to its exon
  flank iff `junction_strand == Strand.POS` (resp. `NEG`) — the unified value now equals `TS_POS` /
  `TS_NEG`. Keep the "side is an exon" gate. This replaces the signature-heuristic AND the
  `~other_exon` AMBIG-exclusion (so AMBIG exons get their correct-strand anchor).

## 4. Remove the annotation hack

Delete: `bp_solver.boundary_junction_strands`; the `boundary_junction_strand` parameter on
`calibrate` and `build_node_geometry`'s *call from calibrate*; the `sj.feather` block in
`pipeline.py`; the `from ..types import Strand` added for the hack stays only if still used by the
unified routing (it is — for `Strand.POS/NEG`). The junction strand now flows in the payload, so
calibrate needs nothing extra.

## 5. Byte-for-byte + tests

1. **Accumulator reference parity** (`tests/native/`): extend the table-driven deposit tests with
   the `strand` arg; assert `boundaries_junction_strand` (C++) == reference on: a + spliced crossing
   (donor + acceptor boundaries get `POS`), a − spliced crossing (`NEG`), an unspliced crossing (0),
   a contained spliced fragment (no boundary set), a multi-junction fragment (all its junctions get
   the one strand), and `merge_from` consistency.
2. **Convention**: a unit test pinning `TS_NEG==2`, `TS_AMBIG==3`, `TS_* == int(RegionStrand.*)`,
   and `transcript_strand_class == coarse_strand_from_signature` elementwise.
3. **Calibration suite** green after the value change.
4. **End-to-end**: re-scan one zero-gDNA condition; assert the routed spliced mass and AMBIG-exon
   residual match the (validated) annotation-hack numbers — the accumulator strand == the annotation
   strand for annotated junctions, so results are identical.

## 6. Rebuild / cache note

C++ changed ⇒ `pip install --no-build-isolation -e .` (in the `rigel` env). The cached dissect
payloads (`/tmp/dissect_cache/*.pkl`) were scanned by the old accumulator and **lack the field** —
they must be re-scanned to test the production path (delete the cache; the small 5 Mb zero-gDNA
conditions are fast, gdna300 slower). `scan_payload` defaulting the field to zeros lets old caches
load without crashing (they then behave as "no junction strand" = the signature fallback), but the
fix is only exercised after a re-scan.

## 7. Decisions / open points

- **Per-boundary single label** (not per-side): faithful to axiom 2, minimal, elegant. Revisit only
  if a real +/− co-occurrence at one position ever appears (none observed).
- **Value convention = `Strand`** (POS=1/NEG=2), the scanner/C++ canonical; `TS_*` unified to it.
- The junction strand is **observed from the motif**, not the annotation — so it covers novel
  junctions and needs no index tables in calibration (the user's elegance requirement).
- Empirically on the synthetic suite this is **neutral** (the AMBIG residual was already shrunk by
  the relaxed-face fix; the dominant residual is message *precision* on anchored exons). It is
  adopted for **correctness/robustness** (real data, novel junctions, AMBIG) and to **remove the
  AMBIG-exclusion latent bug** + the annotation hack, not for a benchmark delta.
