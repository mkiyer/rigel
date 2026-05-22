# Fine-Grained Region Migration Design — v3 (Implementation-Ready)

Date: 2026-05-22
Supersedes: [fine_region_migration_design_v2.md](fine_region_migration_design_v2.md)
            (which deferred the fractional accumulator — this version
            commits to it as the primary design)

---

## 0. Scope and Outcome

Replace the persisted coarse 3-type / 4-strand region partition with a
reference-exhaustive per-reference partition keyed by a 4-bit fine
signature `{intron_pos, intron_neg, exon_pos, exon_neg}`, and replace
the integer mask accumulator with a **fractional float32 accumulator**
that emits, per region, twelve channels:

```text
(contained, boundary_left, boundary_right) × (unspliced, spliced) × (pos, neg)
```

The fractional accumulator preserves a **one-fragment-one-unit
invariant** that the current integer accumulator cannot maintain under
a fine partition (a long fragment touching many small regions
double-counts in the integer formulation). The fractional design
unblocks the per-signature gDNA density estimator and makes
boundary-flux evidence available everywhere, not just at coarse EXON
edges.

**Non-goals.**

- No new gDNA density estimator algorithm in this migration; Phase 4
  *defines the consumer API*, not the new statistical model.
- No reintroduction of multimapper calibration evidence.
- No backward compatibility for pre-v4 indexes or pre-v4 scan payloads.
- No change to EM/locus/scoring code paths beyond what calibration
  produces.

---

## 1. Why v3 differs from v2

v2 deferred the fractional accumulator on three grounds. All three were
wrong on closer inspection:

1. **"Fractional accumulation deletes orientation and contained-vs-flux
   channels."** Incorrect. The 12-channel layout is
   `(3 compartments × 2 splice × 2 strand) = 12`, which is strictly
   richer than the current integer channels: it adds a splice axis,
   preserves strand orientation, and preserves contained-vs-flux.

2. **"The fine partition can be supported by an additive Phase 3 that
   keeps the integer accumulator."** This is true mechanically, but
   the integer accumulator suffers a *fragment-length bias* on a fine
   partition. Under coarse regions, a 200 bp fragment crosses ≤2
   regions; on fine partition the same fragment may touch 5–10. The
   integer `per_region_counts` is then implicitly length-weighted:
   long fragments contribute more "evidence" to density estimators
   than short fragments. This bias is dataset-dependent and grows
   with the fineness of the partition. The fractional formulation
   eliminates it by enforcing `Σ_regions Σ_channels contribution = 1`
   per fragment.

3. **"Float64 is needed for precision."** Reconsidered. Typical RNA-seq
   libraries are 20M–200M fragments, and per-region maxes are far
   below the float32 exact-integer ceiling (`2^24 ≈ 16M`). float32
   halves the dense-merge memory and is precise enough for the
   regime we operate in. We add a debug-build invariant
   (`max_per_region_observed < 2^23`) to catch outliers in unusual
   libraries.

Additionally, v3 **drops the accumulator-side boundary-flux use of
`splicing_anchor_tolerance K`**. The public scan configuration stays:
the resolver still uses `K` for implicit-splice gap classification.
What goes away is the calibration accumulator's hard per-hit
`overlap_bp < K` threshold and the matching boundary-flux denominator
clearance. Under fractional accumulation, micro-overhang aligner errors
are naturally down-weighted by `overlap_bp / total_aligned_bp`, provided
the density consumer uses the small-overlap side of an exon-intron
boundary (§3.6, §9.2). This removal simplifies the hot path and
tightens the one-fragment-one-unit invariant from `approximately 1.0`
to **exactly 1.0**.

The remaining v2 content (signature layout, coarse derivation table,
ref-id mapping fix, INDEX_FORMAT_VERSION = 4, no legacy load path) is
carried forward verbatim.

---

## 2. Design Invariants and Definitions

### 2.1 Canonical signature

Bit layout (uint8, single source of truth in
[signature.py](../../src/rigel/calibration/signature.py) and
[region_signature.h](../../src/rigel/native/calibration/region_signature.h)):

```text
bit 3 (0x8): intron_pos
bit 2 (0x4): intron_neg
bit 1 (0x2): exon_pos
bit 0 (0x1): exon_neg
```

### 2.2 Coarse derivation table

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

Coarse fields are derived, not persisted as primary semantics. The
parity bridge in Phase 1 emits them as columns for downstream code that
has not yet been rewritten.

### 2.3 No legacy load path

`INDEX_FORMAT_VERSION` is bumped from 3 to **4**. Pre-v4 indexes or
`regions.feather` lacking the `signature` column raise `ValueError`
with `Rebuild the index (rigel index ...)`. No migration shim.

---

## 3. Fractional Accumulation — The Algorithm

This is the core new content of v3.

### 3.1 Channel layout

Each region accumulates a **12-vector of float32**:

```text
channel index = compartment * 4 + splice_idx * 2 + strand_idx

compartment in {0=contained, 1=boundary_left, 2=boundary_right}
splice_idx  in {0=unspliced,  1=spliced}
strand_idx  in {0=pos,        1=neg}
```

This packing places strand as the innermost dimension so consumers that
want `pos + neg` summed (strand-agnostic queries) can SIMD-add adjacent
float32 pairs.

### 3.2 Per-fragment preconditions

A fragment is **eligible for calibration accumulation** iff (matches
today's gate, plus one new exclusion):

- not chimeric;
- not a splice artifact;
- single reference across all aligned blocks;
- resolved single strand — **strand-ambiguous fragments are excluded
  with a new `n_excluded_strand_ambig` counter** (today's chimera gate
  already catches most of these; the explicit counter makes the policy
  visible);
- not used as multimapper calibration evidence (same policy as today).

Per-fragment metadata at the hot path entry:
`aligned_blocks: [(start, end)]`, `splice_idx`, `strand_idx`,
`fl_idx` (clipped fragment-length bin).

`splicing_anchor_tolerance` remains an upstream resolver setting for
implicit-splice classification. It is **not** passed into the
fractional accumulator and does not affect mass distribution.

### 3.3 Mass distribution rule

Let `total_aligned_bp = Σ (end - start)` over the fragment's aligned
blocks; `inv_total = 1.0f / total_aligned_bp`.

For each aligned block `(bs, be)`, query `RegionIndex.overlap(ref, bs, be)`
to get overlapping rids. For each rid with region `[rs, re)`
(`RegionIndex.overlap` guarantees `overlap_bp ≥ 1`):

```text
overlap_bp = min(be, re) - max(bs, rs)
w = overlap_bp * inv_total               # float32

cross_left  = (bs < rs)   # block extends past region's left edge
cross_right = (be > re)   # block extends past region's right edge

if !cross_left && !cross_right:
    accum[rid][CONTAINED]      += w     # block fully inside region
elif cross_left && cross_right:
    accum[rid][BOUNDARY_LEFT]  += 0.5f * w   # block fully spans region
    accum[rid][BOUNDARY_RIGHT] += 0.5f * w   # → flux on BOTH sides
elif cross_left:
    accum[rid][BOUNDARY_LEFT]  += w     # block extends past LEFT edge
else: # cross_right
    accum[rid][BOUNDARY_RIGHT] += w     # block extends past RIGHT edge
```

Channel slots use the `(splice_idx, strand_idx)` quadrant computed once
per fragment.

No tolerance gate is applied at the per-hit level. The natural
`overlap_bp / total_aligned_bp` weighting handles micro-overhangs
(§3.6).

### 3.4 Semantic notes

- **Fully spanned = flux on both sides.** A block that fully spans a
  small region crosses both edges; the fragment continues past both.
  This is flux, not containment. (Mass is split 50/50 between the two
  edges; the per-region total is `w`, preserving the one-unit
  invariant.) This is intentional and matches the user's gDNA density
  estimator design: a fully-spanning fragment is the strongest possible
  evidence of contiguous coverage across that region.

- **Splice axis is per-fragment, not per-block.** A fragment is either
  spliced (has at least one N CIGAR op) or not; both blocks of a
  two-block spliced fragment write into `splice_idx = 1` channels.

- **One-fragment-one-unit invariant (exact).**
  `Σ channels Σ rids contribution = 1.0` per fragment, exactly
  (up to float32 rounding). The fine partition is reference-exhaustive
  and `RegionIndex.overlap` returns every overlapping region, so every
  aligned base contributes to exactly one (rid, channel) pair. Debug
  builds assert `|Σ - 1.0| < 1e-5` per fragment.

### 3.5 Worked examples

Block `[80, 140)` over regions `[0,100)`, `[100,150)`,
`total_aligned_bp = 60`:

- R1: overlap_bp = 20, cross_left=false, cross_right=true →
  `boundary_right += 20/60`.
- R2: overlap_bp = 40, cross_left=true, cross_right=false →
  `boundary_left  += 40/60`.

Total mass = 60/60 = 1.0. ✓

Block `[50, 250)` over regions `[0,100)`, `[100,150)`, `[150,300)`,
`total_aligned_bp = 200`:

- R1: overlap_bp = 50, bs(50) > rs(0) and be(250) > re(100) →
  cross_left=false, cross_right=true → `boundary_right += 50/200`.
- R2: overlap_bp = 50, bs(50) < rs(100), be(250) > re(150) →
  cross_left=true, cross_right=true → `boundary_left += 25/200`,
  `boundary_right += 25/200`.
- R3: overlap_bp = 100, bs(50) < rs(150), be(250) < re(300) →
  cross_left=true, cross_right=false → `boundary_left += 100/200`.

Total mass = (50 + 50 + 100) / 200 = 1.0. ✓

Two-block spliced fragment with blocks `[0, 100)` and `[200, 300)`,
both fully contained in distinct exon regions R1 and R3,
`total_aligned_bp = 200`:

- R1 (block 1): overlap_bp = 100, contained →
  `contained += 100/200 = 0.5`.
- R3 (block 2): overlap_bp = 100, contained →
  `contained += 100/200 = 0.5`.

Total mass = 1.0, in `splice_idx=1` channels of both regions. ✓

### 3.6 Why accumulator-side `splicing_anchor_tolerance K` is removed

The historical `K` (default 3 bp) has two current roles:

1. **Implicit-splice resolver slack.** `resolve_context.h` uses `K` to
   decide whether a paired-end gap wholly contains an annotated intron
   with one-sided slack. This remains unchanged.
2. **Calibration boundary-flux slack.** The integer accumulator uses
   `K` to drop region hits whose overlap is below `q=max(K,1)` and to
   require `K` bp clearance on both sides of a boundary before
   incrementing `u_left/u_right`. This role is removed.

The accumulator-side `K` existed because integer boundary flux was a
binary event. A 2 bp intronic overhang at an exon→intron boundary and a
50 bp genuine gDNA crossing both contributed `+1` to the old boundary
counter unless `K` filtered the former out.

Fractional accumulation changes the signal, but only if consumers use
the correct side of the boundary. Let exon `E=[0,100)`, intron
`I=[100,...)`, and an unspliced block of length `L=100` extend from
the exon into the intron by `o=2` bp:

```text
E.boundary_right += (L - o) / L = 0.98
I.boundary_left  += o / L       = 0.02
```

For a genuine gDNA fragment with 50 bp on each side:

```text
E.boundary_right += 0.50
I.boundary_left  += 0.50
```

Therefore:

| Consumer of EXON_INTRON flux      | RNA 98/2 overhang | gDNA 50/50 crossing | Useful? |
| --------------------------------- | ----------------- | ------------------- | :-----: |
| Exon-side flux only               | 0.98              | 0.50                | no      |
| Intron-side flux only             | 0.02              | 0.50                | yes     |
| Sum of exon side + intron side    | 1.00              | 1.00                | no      |

This is the key implementation constraint: **do not estimate the
EXON_INTRON boundary density by summing both sides, and do not use the
exon side for the RNA-suppression argument.** The Phase 3 compatibility
consumer uses the intronic side of EXON_INTRON boundaries (or,
equivalently, the non-exonic side for future EXON_INTERGENIC channels),
with a denominator that mirrors the same fractional side mass (§9.2).

Fractional accumulation does not make RNA micro-overhang noise vanish.
If RNA coverage is enormous, aggregate `N_R * o/L` noise can still be
nontrivial. It does, however, turn an unbounded integer false-positive
count into a small proportional mass that the downstream estimator can
model using splice-axis, strand-symmetry, and neighbor-signature
features. A hard accumulator threshold would destroy that information.

**Removed from the implementation:**

- Per-hit `overlap_bp < q → skip` gate in
  `CalibrationAccumulator::observe(...)`.
- Boundary-endpoint predicate `frag_start + q <= b && frag_end >= b + q`
  (superseded by `cross_left`/`cross_right`).
- `n_below_tolerance` payload counter.
- `CalibrationAccumulator` constructor/member storage for `K`.
- `BamScanner.set_regions(..., splicing_anchor_tolerance)` plumbing
  into the accumulator.
- Use of `boundary_crossing_exposure(..., splicing_anchor_tolerance=K)`
  as the denominator for fractional boundary-flux channels.

**Kept from the implementation:**

- `BamScanConfig.splicing_anchor_tolerance` and the
  `--splicing-anchor-tolerance` CLI flag, now documented as
  implicit-splice resolver slack only.
- `resolve_ctx.set_splicing_anchor_tolerance(...)` in the scan setup.
- Optional scan-payload/result provenance recording the configured `K`,
  explicitly marked as resolver provenance, not accumulator behavior.

### 3.7 Strand-ambig and other exclusions

Each excluded category gets its own counter, all of which sum into the
existing `n_read_names` accountability identity:

| Counter                          | Meaning                                              |
| -------------------------------- | ---------------------------------------------------- |
| `n_excluded_chimera`             | mixed-ref or chimeric fragments                      |
| `n_excluded_splice_artifact`     | failed splice sanity                                 |
| `n_excluded_strand_ambig`        | resolved strand is AMBIG (new explicit counter)      |
| `n_excluded_multimapper`         | multimapper not used for calibration                 |
| `n_observed_for_calibration`     | the float32 quadrant accumulators count these        |

Accountability identity: every input fragment lands in exactly one of
the above counters. `n_below_tolerance` is gone (§3.6).

---

## 4. Payload Schema

### 4.1 New `CalibrationPayload` (native)

```cpp
struct CalibrationPayload {
    // Fractional per-region channels — the primary new evidence block.
    // SoA on the Python boundary; layout choice below.
    float* region_counts;        // dense float32[R, 12], row-major
    int64_t n_regions;

    // Fragment-length mass table stratified by 12-channel quadrant
    // (compartment × splice × strand). Entries are fractional mass,
    // not integer fragment counts.
    double fl_mass_by_channel[12][1024];

    // Per-signature global mass and FL mass table — used by the new
    // gDNA density estimator in Phase 4.
    double per_sig_global_mass[16];      // double for safe global sum
    double fl_mass_by_sig[16][1024];     // fractional mass by sig/bin

    // Accountability counters (int64).
    int64_t n_observed_for_calibration;
    int64_t n_excluded_chimera;
    int64_t n_excluded_splice_artifact;
    int64_t n_excluded_strand_ambig;     // NEW
    int64_t n_excluded_multimapper;

    // Debug-only invariant trackers (compiled out in release).
    #ifdef RIGEL_DEBUG_ACCUM
    double total_mass_emitted;           // expected ≈ n_observed_for_calibration
    int64_t n_unit_invariant_violations;
    float   max_per_region_observed;     // float32 saturation watchdog
    #endif
};
```

`region_counts` is **float32, dense `[R, 12]`** at the Python boundary.
Per-worker accumulation uses a sparse representation (see §5).

`per_sig_global_mass` uses `double` because it accumulates across the
whole genome and is a useful summary that benefits from extra
precision; it is O(16) so cost is trivial.

`splicing_anchor_tolerance` is not a `CalibrationAccumulator` field in
Phase 3. The native scanner may still include it in the Python scan
dict as run provenance because the resolver used it upstream for
implicit-splice classification.

### 4.2 Removed from the payload

The following are removed from the native `CalibrationPayload`. They
have no role under fractional accumulation:

- `per_region_counts[R, 8]` (replaced by `region_counts[R, 12]`)
- `intron_contained_counts_by_orient[R, 3]` (orientation now lives in
  the `strand` axis of `region_counts`)
- `exon_contained_counts_by_orient[R, 3]` (ditto)
- `u_left[R]`, `u_right[R]` (replaced by
  `region_counts[r, BOUNDARY_LEFT, *, *]` summed)
- `u_left_by_orient[R, 3]`, `u_right_by_orient[R, 3]` (ditto)
- `global_counts[8]` (replaced by `per_sig_global_mass[16]`)
- `fl_hist[8, 1024]` (replaced by `fl_mass_by_channel[12, 1024]` and
  `fl_mass_by_sig[16, 1024]`)

This is a hard cut. The integer accumulator infrastructure goes away
in Phase 3.

### 4.3 Python schema (`CalibrationScanPayload`)

```python
@dataclass(frozen=True)
class CalibrationScanPayload:
    # Fractional region evidence (SoA on Python side).
    contained_unspliced_pos:    np.ndarray   # float32[R]
    contained_unspliced_neg:    np.ndarray
    contained_spliced_pos:      np.ndarray
    contained_spliced_neg:      np.ndarray
    boundary_left_unspliced_pos:  np.ndarray
    boundary_left_unspliced_neg:  np.ndarray
    boundary_left_spliced_pos:    np.ndarray
    boundary_left_spliced_neg:    np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray
    boundary_right_spliced_pos:   np.ndarray
    boundary_right_spliced_neg:   np.ndarray
    # — or equivalently, a single float32[R, 12] view backed by the
    # same buffer, exposed via named-property accessors.

    fl_mass_by_channel: np.ndarray             # float64[12, 1024]
    per_sig_global_mass: np.ndarray            # float64[16]
    fl_mass_by_sig:     np.ndarray             # float64[16, 1024]

    # Accountability
    n_observed_for_calibration: int
    n_excluded_chimera: int
    n_excluded_splice_artifact: int
    n_excluded_strand_ambig: int
    n_excluded_multimapper: int

    # Scan-config provenance only. The resolver used this value for
    # implicit-splice classification; the fractional accumulator did
    # not use it for mass distribution.
    splicing_anchor_tolerance: int = 0

    # Convenience derived (lazy properties, not stored)
    @property
    def contained(self) -> np.ndarray:        # float32[R]; sum over splice×strand
        ...
    @property
    def boundary_left(self) -> np.ndarray:    # float32[R]
        ...
    @property
    def boundary_right(self) -> np.ndarray:   # float32[R]
        ...
    @property
    def region_total_mass(self) -> np.ndarray:  # float32[R]; contained + bL + bR
        ...
```

Schema validation at `from_scan_dict` time:

- All 12 SoA arrays are float32, shape `(R,)`, finite, non-negative.
- `fl_mass_by_channel.shape == (12, 1024)`, dtype float64,
  non-negative.
- `per_sig_global_mass.shape == (16,)`, dtype float64, non-negative,
  sum approximately equals `n_observed_for_calibration`.
- `fl_mass_by_sig.shape == (16, 1024)`, dtype float64, non-negative.
- `splicing_anchor_tolerance >= 0`; this is provenance for the
  resolver, not an accumulator tolerance.
- Total-mass identity (exact): `region_counts.sum() ≈
  n_observed_for_calibration` within float32 tolerance
  (`|Σ − N_obs| / N_obs < 1e-5`). There is no dropped-mass term —
  with no tolerance gate, every observed fragment contributes exactly
  1.0 to the total.

---

## 5. Performance Architecture

A dense per-worker `float32[R, 12]` is too large at scale. For
`R ≈ 10M` (post-fine-partition GENCODE), dense is 480 MB per worker;
at 16 workers, 7.7 GB.

Real per-worker hit sets are tiny. A name-sorted BAM chunk processed
by a single worker typically touches `10^4 – 10^5` distinct regions
(the contiguous mapped portion of the genome chunk).

### 5.1 Per-worker accumulator: sparse hash, AoS

```cpp
struct WorkerAccum {
    // Open-addressing flat hash map: rid -> slot index in `slots`.
    // robin_hood::unordered_flat_map<int32_t, uint32_t> idx;
    // OR: vendored simple linear-probing table; we already vendor cgranges.
    HashMap<int32_t, uint32_t> idx;

    // Dense-on-touch slot vector. Each slot is a 12-float row (AoS).
    // Pre-reserve to 64 K slots; double on demand.
    std::vector<std::array<float, 12>> slots;
};
```

Per-worker hot-path write:

```cpp
inline void worker_add(WorkerAccum& w, int32_t rid, int channel, float mass) {
    auto [it, inserted] = w.idx.try_emplace(rid, w.slots.size());
    if (inserted) w.slots.emplace_back();   // zero-initialized
    w.slots[it->second][channel] += mass;
}
```

Per-worker memory: `O(distinct_regions_touched × 12 × 4 B + hashmap overhead)`,
typically a few MB. Independent of `R`.

### 5.2 Global merge: dense SoA

After all workers finish, allocate the dense `float32[R, 12]` once and
sum workers into it in deterministic worker order (matches current
behavior). Then write to the Python payload as **12 separate
`float32[R]` arrays** (SoA), backed by either:

- 12 independent allocations (simplest, no aliasing issues), or
- a single `float32[R, 12]` buffer with strided views per channel
  (less alloc churn).

SoA on the Python boundary is the right cache layout for consumers
that read one channel at a time across many rids (which is what
density estimators do).

### 5.3 Hot-loop micro-optimizations

These apply to the per-fragment, per-block, per-rid inner loop in
`CalibrationAccumulator::observe`:

1. **Precompute `inv_total` once per fragment** (no division per
   rid).
2. **Branchless compartment selection.** Encode
   `(cross_left, cross_right)` as a 2-bit index into a 4-entry
   small lookup:
   ```cpp
   //  cl cr | action
   //   0  0 | contained  += w           (weight_left = 0,    weight_right = 0,   weight_contained = w)
   //   0  1 | bR         += w           (0,   w,    0)
   //   1  0 | bL         += w           (w,   0,    0)
   //   1  1 | bL+=w/2, bR+=w/2          (w/2, w/2,  0)
   ```
   Implement as three multipliers indexed by the 2-bit code:
   ```cpp
   static constexpr float kMultBL[4] = {0.0f, 0.0f, 1.0f, 0.5f};
   static constexpr float kMultBR[4] = {0.0f, 1.0f, 0.0f, 0.5f};
   static constexpr float kMultC [4] = {1.0f, 0.0f, 0.0f, 0.0f};
   const uint32_t code = (cross_left << 1) | cross_right;
  const int q = splice_idx * 2 + strand_idx;
  slot[q + 0] += w * kMultC [code];  // contained
  slot[q + 4] += w * kMultBL[code];  // boundary_left
  slot[q + 8] += w * kMultBR[code];  // boundary_right
   ```
   Two of the three adds are zero each iteration; modern OOO cores
   pipeline these effectively and the codegen avoids branch
   mispredictions.

3. **Precompute `region_starts[R]`, `region_ends[R]` as flat int64
   arrays** at `RegionIndex::set(...)` time. Inner loop reads
   these by `rid` — fully predictable, vectorizable for the
   overlap_bp computation across contiguous hit ranges.

4. **Pre-compute per-fragment `quadrant = splice_idx * 2 +
  strand_idx`**; channels are `quadrant + {0, 4, 8}` for
  contained, boundary_left, and boundary_right respectively. This
  matches §3.1's `compartment * 4 + splice_idx * 2 + strand_idx`
  layout.

5. **Per-fragment debug invariant check** is `#ifdef RIGEL_DEBUG_ACCUM`
   and compiled out in release. The cost in debug is one float sum
   accumulated across the fragment's inner loop.

6. **No atomics on the hot path.** Per-worker accumulators are private.

7. **Hash-map sizing.** Reserve `WorkerAccum::idx` to a power-of-two
   initial bucket count (e.g. 1 << 16 = 65 536) and grow geometrically.
   For typical chunk sizes this avoids any reallocation.

### 5.4 Estimated cost

The hot path adds per overlap hit (no tolerance branch):

- 2 int64 loads (region_starts, region_ends),
- 2 min/max + 1 sub for overlap_bp,
- 1 mul for `w = overlap_bp * inv_total` (float32),
- 1 hashmap probe (amortized O(1)),
- 3 float adds to a 48 B AoS slot (all in one cache line).

Compared to the current integer accumulator's per-hit cost (sorted-merge
deduplication, then int64 ++ into 8-channel array + orientation channel
+ FL histogram update + tolerance gate + boundary endpoint predicate),
this is fewer instructions and one fewer mispredictable branch. The
AoS slot fits in one cache line.

The 12×1024 channel FL mass table (96 KiB) and 16×1024 by-signature FL
mass table (128 KiB) are both per-worker private and both fit in L2.

---

## 6. Phase 0 — Cleanup, Signature Helpers, Ref-id Fix

Identical to v2 §3 except cross-references are updated. Verbatim
deliverables:

- `src/rigel/calibration/signature.py` (new) — single source of truth
  for the 4-bit layout and derivation helpers; also defines the
  channel layout constants (`CHAN_CONTAINED`, `CHAN_BOUNDARY_LEFT`,
  `CHAN_BOUNDARY_RIGHT`, `N_CHANNELS = 12`,
  `channel_index(compartment, splice_idx, strand_idx)`).
- `src/rigel/native/calibration/region_signature.h` (new) — C++ mirror,
  including `channel_index(...)` constexpr.
- `tests/test_signature.py` — exhaustive 16-state derivation test;
  channel-index round-trip test.
- `tests/test_signature_layout_native.py` — Python↔C++ constants
  agree.
- `_wire_calibration_regions(...)` ref-id fix in
  [pipeline.py](../../src/rigel/pipeline.py): use
  `index.ref_name_to_id`; warn on dropped refs.
- Remove `set_regions_legacy(...)` and any missing-region fallback in
  `scan_and_buffer(...)`.

**Acceptance gate.** All existing tests pass unchanged. No persisted
format change. No `INDEX_FORMAT_VERSION` bump.

---

## 7. Phase 1 — Fine Region Index Builder + Parity Bridge

### 7.1 Deliverables

- `build_fine_region_table(transcripts, ref_lengths) -> pd.DataFrame`
  in [regions.py](../../src/rigel/calibration/regions.py).
- New on-disk `regions.feather` schema (see §7.3).
- `INDEX_FORMAT_VERSION = 4`. Hard-fail loaders on `version < 4`.
- Retire `emit_regions(...)` / `RegionRecord`; replace with the
  vectorized builder.
- `validate_against_ref_lengths(...)` extended for new fields.
- `RegionArrays` gains a `signature: uint8[R]` field. Coarse fields
  (`type`, `strand`, `boundary_flux_*`) remain *as derived columns*
  through Phase 3.

### 7.2 Builder algorithm

Per reference:

1. Collect events from every non-synthetic transcript:
   - exon: `(start, exon_<strand>, +1)`, `(end, exon_<strand>, −1)`.
   - intron (between consecutive exons of the same transcript):
     `(start, intron_<strand>, +1)`, `(end, intron_<strand>, −1)`.
2. Sort events by position. Sweep left to right with four running
   counters `(ip, in_, ep, en)`.
3. Between consecutive distinct event positions `[p_k, p_{k+1})`, the
   active 4-tuple is constant. Pack as a uint8 signature.
4. Emit a region only when the signature changes from the previously
   emitted segment (adjacent-merge semantics).
5. Reference boundaries `0` and `ref_length` are implicit start/end.
6. Empty references emit a single `signature == 0` row covering
   `[0, ref_length)`.

Per-reference cost: `O(E log E)`, `E` = number of transcript event
edges on that reference.

### 7.3 `regions.feather` schema

Columns (locked dtypes):

| Column                | dtype     | Notes                                       |
| --------------------- | --------- | ------------------------------------------- |
| `region_id`           | int64     | Globally sequential = row index             |
| `ref_name`            | string    | FASTA reference name                        |
| `start`               | int64     | 0-based inclusive                           |
| `end`                 | int64     | 0-based exclusive                           |
| `length`              | int64     | `end - start`                               |
| `signature`           | uint8     | Canonical 4-bit fine state                  |
| `intron_pos`          | bool      | == `signature & 0x8 != 0`                   |
| `intron_neg`          | bool      | == `signature & 0x4 != 0`                   |
| `exon_pos`            | bool      | == `signature & 0x2 != 0`                   |
| `exon_neg`            | bool      | == `signature & 0x1 != 0`                   |
| `type`                | uint8     | derived coarse (Phase 1 bridge)             |
| `strand`              | uint8     | derived coarse (Phase 1 bridge)             |
| `boundary_flux_left`  | bool      | coarse-parity rule, Phase 1 bridge          |
| `boundary_flux_right` | bool      | coarse-parity rule, Phase 1 bridge          |
| `left_signature`      | uint8     | adjacent-region signature on the left       |
| `right_signature`     | uint8     | adjacent-region signature on the right      |
| `boundary_kind_left`  | uint8     | precomputed left boundary classification    |
| `boundary_kind_right` | uint8     | precomputed right boundary classification   |

The `left_signature`/`right_signature` columns store sentinel `0xFF`
for reference-end edges (no neighbor on that side). They are computed
in the builder via a single linear pass; cost is negligible relative
to the build itself.

`boundary_kind_{left,right}` is a uint8 enum with values:

```text
0 = NONE                    # no neighbor (reference end)
1 = EXON_INTRON             # exon bit on one side, intron bit on other, same strand
2 = EXON_INTERGENIC         # exon side ↔ intergenic
3 = INTRON_INTERGENIC       # intron side ↔ intergenic
4 = EXON_EXON_STRAND_SWITCH # exon flag set on both, but only opposite-strand
5 = OTHER                   # any other transition
```

`boundary_kind` is **additive**; nothing in Phase 1–3 reads it, but
it's cheap to compute now and saves a future migration. Phase 4's
density estimator may consume it.

### 7.4 Boundary-flux parity rule (Phase 1 only)

`boundary_flux_{left,right}` uses the same coarse-parity rule as v2 §2.4:

- `boundary_flux_left[i] = (coarse_type[i] == EXON) AND
                            (i > 0) AND (ref_name[i-1] == ref_name[i])`
- `boundary_flux_right[i] = (coarse_type[i] == EXON) AND
                             (i+1 < R) AND (ref_name[i+1] == ref_name[i])`

These columns exist only to keep the integer accumulator and current
density consumers bit-identical through Phase 2. They are deleted in
Phase 3 along with the integer accumulator.

### 7.5 Removed schema fields

`tx_pos_bp`, `tx_neg_bp`, `exon_pos_bp`, `exon_neg_bp` are not
persisted. No current downstream code reads them (verified by grep).

### 7.6 Loader cutover

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

`index.py` metadata check raises the same actionable error on
`INDEX_FORMAT_VERSION < 4`.

### 7.7 Tests

`tests/test_regions.py` — fine-builder matrix (same as v2 §4.7):

| Case                                                             | Signatures expected (in order)              |
| ---------------------------------------------------------------- | ------------------------------------------- |
| Empty reference                                                  | `[0x0]`                                     |
| Single (+) tx, 2 exons + 1 intron, intergenic flanks             | `0x0, 0x2, 0x8, 0x2, 0x0`                   |
| Two (+) tx with overlapping exons                                | merges identical adjacent signatures        |
| (+) and (−) tx overlapping on exons                              | includes `0x3` (`exon_pos+exon_neg`)        |
| (+) exon overlapping (−) intron                                  | includes `0x6` (`intron_neg+exon_pos`)      |
| (+) exon overlapping (+) intron (other tx)                       | includes `0xA` (`intron_pos+exon_pos`)      |
| Single-exon tx                                                   | `0x0, 0x{1 or 2}, 0x0`                      |
| Read-through tx spanning two genes                               | no signature collapse where bits differ     |
| Multi-reference                                                  | `region_id` globally sequential             |
| Intergenic-only reference                                        | single `0x0` row                            |

`tests/test_region_parity.py` (new) — for each scenario assert:

- `signature == pack_signature(intron_pos, intron_neg, exon_pos, exon_neg)`,
- `type == coarse_type_from_signature(signature)`,
- `strand == coarse_strand_from_signature(signature)`,
- `boundary_flux_*` matches §7.4,
- `left_signature`/`right_signature` are 0xFF at reference edges and
  match neighboring `signature` otherwise.

`tests/test_index_format_version.py` — building an index writes
`INDEX_FORMAT_VERSION == 4`; pre-v4 feather raises `ValueError`.

### 7.8 Golden refresh (single commit)

Refresh all `tests/golden/` with `pytest tests/ --update-golden`.
Expected deltas: only `regions.feather` content and any downstream
fields that scale with region count (e.g. region-driven exposure
sums). Treat unexplained changes to `transcript_df`, `gene_df`,
`scalars.json` as regressions.

### 7.9 Acceptance gates

- New schema validates on all fixtures.
- Snapshot-based parity: capture v3-builder outputs of the *derived*
  coarse columns from a temporary `tests/_v3_builder_snapshot/`
  taken before this PR lands; new builder produces matching
  derived columns row-aligned.
- All tests pass with refreshed goldens.
- `INDEX_FORMAT_VERSION == 4` enforced.

---

## 8. Phase 2 — Native Region Index Carries Signature

Identical in spirit to v2 §5. Small, low-risk phase. Native
accumulator behavior is **unchanged**.

### 8.1 Deliverables

- `RegionIndex::set(...)` gains a `signatures` parameter and exposes
  `signature(rid)`, `left_signature(rid)`, `right_signature(rid)`,
  `boundary_kind_left(rid)`, `boundary_kind_right(rid)` accessors.
- nanobind binding updated; Python `set_regions(...)` API gains the
  new arrays.
- `pipeline._wire_calibration_regions(...)` derives and passes
  signatures from `region_df["signature"]`; passes neighbor
  signatures and boundary_kinds from the new columns.

### 8.2 New native `set` signature

```cpp
void set(const int32_t* ref_ids,
         const int64_t* starts,
         const int64_t* ends,
         const uint8_t* signatures,            // NEW required
         const uint8_t* left_signatures,       // NEW required (0xFF sentinel)
         const uint8_t* right_signatures,      // NEW required
         const uint8_t* boundary_kind_left,    // NEW required
         const uint8_t* boundary_kind_right,   // NEW required
         const uint8_t* type_masks,            // kept; will be removed in Phase 3
         const uint8_t* strands,               // kept; removed in Phase 3
         int64_t n_regions,
         int32_t n_refs);
```

### 8.3 Tests

- `tests/test_region_index_native.py` (extend): construct fixture with
  all 16 signatures + each `boundary_kind_*` value; assert accessors.

### 8.4 Acceptance gates

- Existing scanner tests pass unchanged.
- Goldens untouched (refreshed in Phase 1).

---

## 9. Phase 3 — Fractional Accumulator Cutover

This is the largest single phase. It replaces the integer
accumulator and rewrites every direct consumer of integer per-region
counts in one PR. Phases 0–2 deliberately set this up so it is a
mechanical translation, not a research effort.

### 9.1 Deliverables

- **Native:**
  - `CalibrationPayload` replaces integer arrays with the §4.1
    float32 layout.
  - `CalibrationAccumulator::observe(...)` rewritten per §3.
  - `CalibrationAccumulator` no longer accepts or stores
    `splicing_anchor_tolerance`; `BamScanner.set_regions(...)` no
    longer forwards it to the accumulator.
  - Per-worker sparse hash map per §5.1; deterministic global merge
    per §5.2.
  - Channel index helpers from `region_signature.h`.
  - Remove `intron_contained_counts_by_orient`,
    `exon_contained_counts_by_orient`, `u_left*`, `u_right*`,
    `per_region_counts[R, 8]`, `global_counts[8]`, `fl_hist[8, 1024]`.
- **Python:**
  - `CalibrationScanPayload` rewritten per §4.3.
  - `density_global.compute_global_densities(...)` rewritten to read
    `region_counts[R, 12]` SoA fields directly.
  - Boundary denominators rewritten from `boundary_crossing_exposure(K)`
    to the exact fractional side-exposure helper in §9.2.
  - `_regional_exposure.RegionalGdnaExposure.build(...)` rewritten.
  - `locus_prior.assemble_priors(...)` rewritten.
  - `BamScanConfig.splicing_anchor_tolerance` and
    `--splicing-anchor-tolerance` retained, but documentation/help text
    narrowed to implicit-splice resolver slack.
  - `RegionArrays` loses `type`, `strand`, `bf_left`, `bf_right`
    (verify by removal + `pytest`).
  - `regions.feather`: `type`, `strand`, `boundary_flux_*` remain in
    the schema for `pandas` ergonomics only. Hot code reads
    `signature` directly.

### 9.2 New density consumers — minimal API

Phase 3 keeps the density estimator's *shape* close to today so this
phase is a translation, not a redesign. The new gDNA density
estimator algorithm lives in Phase 4.

For each existing coarse class, the Phase 3 implementation maps to a
fractional read:

| Coarse class today        | Phase 3 read                                                            |
| ------------------------- | ----------------------------------------------------------------------- |
| INTERGENIC numerator      | `contained.sum()` over rids with `signature == 0x0`                     |
| INTRON-PURE numerator     | `contained.sum()` over rids with intron-only signatures `{0x4,0x8,0xC}` |
| EXON-CONTAINED numerator  | `contained.sum()` over rids with any exon flag                          |
| EXON-INTRON BOUNDARY num. | **Intron-side flux only**: `boundary_left` or `boundary_right` on the region side whose `boundary_kind_* == EXON_INTRON` and whose signature has intron flag(s) but no exon flag. Do not sum both sides of the same boundary. |

Each `.sum()` is a strand-marginalized sum (`pos + neg`) and a
splice-marginalized sum (`unspliced + spliced`) unless the consumer
specifically wants the stratification.

Boundary denominators must mirror the fractional numerator, not the old
integer `boundary_crossing_exposure(K)` event count. Add a helper such
as `fractional_boundary_side_exposure(length, gdna_fl)` that computes
the expected mass emitted to **one side** of a region under the exact
§3.3 routing rule. For a region `[0,S)` and a left-boundary crossing
fragment `[a,a+ℓ)` with `a < 0 < a+ℓ`, the per-start contribution is:

```text
overlap = min(a + ℓ, S) - 0
w = overlap / ℓ
side_mass = 0.5 * w   if a + ℓ > S   # fragment fully spans the region
          = w         otherwise      # crosses only the left boundary
```

Sum `side_mass` over all integer starts `a ∈ [1-ℓ, -1]`, then average
over the gDNA FL PMF. The right-boundary formula is symmetric. For
large regions (`S ≥ ℓ` for all relevant ℓ), this reduces to
`Σℓ h(ℓ) * (ℓ - 1) / 2`, i.e. half the old strict-crossing exposure at
`K=0`. For small fine regions, the exact helper is required because
fully-spanned regions split mass 50/50 across both sides.

Containment denominators continue to be computed from region `length`
and the same eligibility filters.

### 9.3 Tests (Phase 3)

`tests/test_accumulator_fine.py` (new):

- **Unit invariant**: single fragment whose blocks total `L` bp;
  assert `sum(all channels, all rids) ≈ 1.0` within float32
  tolerance.
- **Worked examples from §3.5** as explicit assertions.
- **Splice axis routing**: spliced and unspliced fragments otherwise
  identical land in different quadrants.
- **Strand axis routing**: pos and neg fragments otherwise identical
  land in different quadrants.
- **Strand-ambig exclusion**: a fragment whose strand resolves AMBIG
  increments `n_excluded_strand_ambig` and does *not* contribute to
  `region_counts`.
- **Micro-overhang routing (no tolerance gate)**: a fragment with a
  small intronic overhang `o` at an exon→intron boundary contributes
  exactly `(L − o) / L` to exon `boundary_right` and `o / L` to intron
  `boundary_left`; total per-fragment mass is exactly 1.0. Test for
  `o = 1, 2, 3` and `L = 100` against expected float32 values.
- **Boundary consumer guardrail**: the EXON_INTRON density numerator
  uses intron-side flux only. A 98/2 RNA overhang must contribute
  `0.02`, not `0.98` and not `1.0`, to that numerator.
- **Fractional side exposure denominator**: for a large region and a
  degenerate `L=100` gDNA FL model,
  `fractional_boundary_side_exposure == 49.5`; for a short region,
  assert the exact enumerated value including 50/50 fully-spanned
  routing.
- **Multi-block spliced fragment**: two blocks, each fully contained
  in distinct regions, sum to 1.0 across both regions' `contained`
  channel.
- **Resolver K still wired**: setting `BamScanConfig.splicing_anchor_tolerance`
  still changes implicit-splice classification in the existing
  resolver test matrix; it does not change accumulator mass routing.

`tests/test_scan_payload_fractional.py` (new):

- Schema, dtype, shape, non-negativity, finiteness.
- Global identity (exact): `Σ region_counts ≈ n_observed_for_calibration`
  within float32 tolerance.
- `per_sig_global_mass.sum() ≈ n_observed_for_calibration`.

`tests/test_density_global.py` (rewrite):

- All cases reworked against fractional reads.
- Add a parity test only where coarse and fractional results are
  *expected* to match in degenerate cases (e.g. all fragments fully
  contained in single coarse-EXON regions on a non-overlapping
  partition); use this as a sanity floor.
- Add a non-parity test for exon-intron micro-overhangs: old integer
  `u_*` would count 1, fractional intron-side flux counts `o/L`.

`tests/test_regional_exposure.py`, `tests/test_locus_prior.py`:
rewritten similarly.

### 9.4 Goldens

Goldens **will move materially** in this phase. This is expected
and is the point of the migration. Update protocol:

- Run `pytest tests/ --update-golden`.
- In the PR description, table every changed golden file with a
  one-line explanation grounded in the algorithm (length-bias
  removal, splice axis stratification, fractional flux on
  fully-spanned regions, etc.).
- Reviewers spot-check at least 3 fixtures by hand against the
  algorithm in §3.

### 9.5 Performance acceptance

On the largest available benchmark BAM:

- BAM scan wall-clock within **±15 %** of Phase 2 baseline.
- Peak RSS across all workers within **+25 %** of Phase 2 baseline.
  (Fractional accumulator's sparse maps should keep total memory
  *lower* than the dense integer per-region arrays the integer
  accumulator allocated for orientation channels; the budget allows
  for measurement noise.)

If wall-clock regresses beyond 15 %, treat as a Phase 3 blocker and
profile before merging. The likely culprit if it happens is hash
map contention or grow-realloc; tune `WorkerAccum::idx` capacity.

### 9.6 Correctness acceptance

- Synthetic locus sweep (`scripts/sim/locus_sweep.py`,
  `locus_simple_baseline.yaml`) shows no regression in mRNA/nRNA/gDNA
  relative error.
- Full-scale Armis2 benchmarks
  (`scripts/benchmarking/configs/default.yaml`) over the 13
  conditions: no condition regresses mRNA error by more than **5 %
  relative**. We expect (and should observe) gDNA improvement on
  `gdna_high_ss_*_nrna_none` conditions.

---

## 10. Phase 4 — New gDNA Density Estimator (Forward)

Phase 4 is **out of scope for this migration's implementation work**;
it is the research direction unlocked by Phases 0–3. v3 specifies the
*interface* the new estimator will consume so Phase 3's deliverables
satisfy its needs.

### 10.1 Interface contract

A future per-signature density estimator consumes:

- `signature[R]`, `length[R]`, `left_signature[R]`, `right_signature[R]`,
  `boundary_kind_left[R]`, `boundary_kind_right[R]`.
- `region_counts[R, 12]` SoA (the 12 channels).
- `per_sig_global_mass[16]`.
- `fl_mass_by_channel[12, 1024]` and `fl_mass_by_sig[16, 1024]`.

The estimator's output schema (`gdna_density_by_signature`,
`gdna_density_by_boundary_kind`, per-region shrinkage weights) is
TBD and not specified here.

### 10.2 Properties the estimator can exploit

- **Universal boundary flux.** Every fine region has a left and a
  right boundary kind; flux evidence is available everywhere, not
  only at coarse EXON edges.
- **Splice/unspliced separability.** gDNA cannot be spliced; the
  spliced quadrants are a near-pure RNA channel and can be used to
  bound the RNA contribution to ambiguous regions.
- **Strand symmetry of gDNA.** For genuinely gDNA-dominated regions,
  pos and neg channels should be roughly equal; strand imbalance is
  evidence against the gDNA hypothesis at exonic signatures.
- **Per-region total mass `M_r = contained + boundary_left + boundary_right`**
  is a well-defined "fragments per region" exposure that scales
  correctly under fine partitioning (the integer accumulator could
  not provide this).

### 10.3 Phase 4 deliverables (placeholder)

To be specified in a follow-up design after Phase 3 lands and we have
real fractional payloads to inspect on benchmark BAMs.

---

## 11. Test Matrix Summary (all phases)

Unit:
- `test_signature.py`, `test_signature_layout_native.py`,
- `test_regions.py`, `test_region_parity.py`,
- `test_region_index_native.py`,
- `test_accumulator_fine.py`, `test_scan_payload_fractional.py`,
- `test_density_global.py`, `test_regional_exposure.py`,
- `test_locus_prior.py`, `test_index_format_version.py`,
- `test_pipeline_wiring.py`.

Integration:
- `test_cli.py` (`rigel index`, `rigel quant` smoke).
- `test_golden_output.py` against refreshed goldens.

Synthetic benchmark:
- `scripts/sim/locus_sweep.py` on `locus_simple_baseline.yaml`
  diffed against `scripts/benchmark/golden/locus_simple_baseline/`.

Full-scale (Phase 3 gate):
- `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml`
- `python -m scripts.benchmarking analyze -c .../default.yaml -o results/fractional_phase3/`

---

## 12. Golden Update Protocol

| Phase | Goldens touched? | Notes                                                         |
| ----- | ---------------- | ------------------------------------------------------------- |
| 0     | No               |                                                               |
| 1     | Yes, single PR   | Region-count-driven changes only; document per fixture.       |
| 2     | No               |                                                               |
| 3     | Yes, single PR   | Material; per-fixture explanation required.                   |
| 4     | TBD              | Future PR.                                                    |

Always `pytest tests/ --update-golden`; never hand-edit golden files.

---

## 13. Risk Register

| Risk                                                              | Likelihood | Severity | Mitigation                                                                                  |
| ----------------------------------------------------------------- | ---------- | -------- | ------------------------------------------------------------------------------------------- |
| float32 precision loss at extreme depth                           | low        | medium   | Debug invariant `max_per_region_observed < 2^23`; documented escape hatch to float64.       |
| Region count growth blows memory budget                           | medium     | high     | Sparse per-worker hash; dense merge once; alert at index build when R > 5e6 / ref.          |
| Goldens accepted without inspection                               | medium     | high     | Phase 3 PR lists every changed golden; reviewers spot-check 3+ fixtures.                    |
| Hash-map contention or rehashing dominates hot path               | medium     | medium   | Pre-reserve buckets; benchmark before merging Phase 3.                                      |
| Boundary semantics drift surprises downstream gDNA model          | medium     | high     | Phase 4 is forward; Phase 3 keeps API compatible with current consumers.                    |
| Bit-layout typo reintroduced by hand-rolled signature packing     | low        | high     | `signature.py` is the only writer; Python↔C++ layout test in CI.                            |
| Ref-id fix regresses BAMs missing refs                            | low        | medium   | Test "BAM header subset of index" in `test_pipeline_wiring.py`.                             |
| Removing accumulator-side `K` increases boundary-flux contamination on real data | medium     | medium   | Use intron-side EXON_INTRON flux only; strand/splice diagnostics in Phase 3 benchmarks; if still needed, add an explicit small-overhang side channel in a later estimator PR rather than silently dropping mass in the accumulator. |
| Fragment-length bias proves smaller than expected, Phase 3 looks  | low        | low      | Synthetic locus sweep + Armis2 benchmarks gate the PR; if no improvement, document and     |
| like a no-op on real data                                         |            |          | accept (the structural cleanup still has value).                                            |

---

## 14. Rollout

Single repository — no feature flag. Phases land as four PRs (0, 1,
2, 3). Each PR is mergeable and leaves the codebase green.
Consumers must rebuild indexes after the Phase 1 PR merges.

CHANGELOG entry (added in Phase 1 PR; expanded in Phase 3 PR):

```text
## Unreleased

### Breaking
- regions.feather format bumped to INDEX_FORMAT_VERSION = 4. Old
  indexes will fail to load with an explicit "Rebuild the index"
  message. The new schema persists a 4-bit fine signature
  (intron_pos, intron_neg, exon_pos, exon_neg) per region. Rebuild
  with `rigel index --fasta ... --gtf ...`.
- Calibration scan payload schema is replaced: integer 8-state
  per-region counts and orientation channels are removed; a
  fractional float32 12-channel per-region array
  (contained/boundary_left/boundary_right × spliced/unspliced ×
  pos/neg) becomes the primary evidence block. Downstream gDNA
  density estimators consume this directly.
```

---

## 15. Implementation Checklist

### Phase 0
- [ ] Create `src/rigel/calibration/signature.py` (signature + channel
      layout helpers).
- [ ] Create `src/rigel/native/calibration/region_signature.h` (C++
      mirror).
- [ ] Add `tests/test_signature.py` (16-state coverage; channel index
      round-trip).
- [ ] Add `tests/test_signature_layout_native.py`.
- [ ] Fix `_wire_calibration_regions(...)` ref-id mapping per §6.
- [ ] Remove `set_regions_legacy(...)` and missing-region fallback.

### Phase 1
- [ ] Implement `build_fine_region_table(...)` per §7.2.
- [ ] Update `regions.feather` schema per §7.3 (incl. `left_signature`,
      `right_signature`, `boundary_kind_{left,right}`).
- [ ] Bump `INDEX_FORMAT_VERSION` to 4.
- [ ] Update `validate_against_ref_lengths(...)` per §7.7.
- [ ] Update `load_regions(...)` to require `signature`.
- [ ] Update `RegionArrays.from_region_df(...)` to carry `signature`.
- [ ] Snapshot pre-PR coarse-derived columns into
      `tests/_v3_builder_snapshot/`; add parity test.
- [ ] Refresh all goldens (single commit); list changed files in PR.
- [ ] Add CHANGELOG entry per §14.

### Phase 2
- [ ] Add `signatures`, neighbor signatures, `boundary_kind_*` to
      native `RegionIndex::set(...)`.
- [ ] Expose accessors on `RegionIndex`.
- [ ] Update nanobind binding and Python wiring.
- [ ] Extend `tests/test_region_index_native.py`.

### Phase 3
- [ ] Rewrite `CalibrationPayload` (native) per §4.1; remove integer
      arrays.
- [ ] Rewrite `CalibrationAccumulator::observe(...)` per §3 (no
      tolerance gate; `cross_left`/`cross_right` routing).
- [ ] Remove `splicing_anchor_tolerance` from `CalibrationAccumulator`
      constructor/member state and from `BamScanner.set_regions(...)`;
      keep `resolve_ctx.set_splicing_anchor_tolerance(...)` in scan
      setup.
- [ ] Keep `BamScanConfig.splicing_anchor_tolerance` and
      `--splicing-anchor-tolerance`, but update docs/help text to state
      that it controls implicit-splice resolver slack only.
- [ ] Implement per-worker sparse hash + global dense merge per §5.
- [ ] Add debug-only per-fragment unit invariant check.
- [ ] Rewrite `CalibrationScanPayload` per §4.3 (SoA float32 channels).
- [ ] Implement `fractional_boundary_side_exposure(...)` and remove
  `boundary_crossing_exposure(K)` from fractional boundary density
  consumers.
- [ ] Rewrite `compute_global_densities(...)`.
- [ ] Rewrite `RegionalGdnaExposure.build(...)`.
- [ ] Rewrite `assemble_priors(...)`.
- [ ] Drop `type`, `strand`, `bf_left`, `bf_right` from
      `RegionArrays`.
- [ ] Add `tests/test_accumulator_fine.py`,
      `tests/test_scan_payload_fractional.py`.
- [ ] Rewrite `tests/test_density_global.py`,
      `tests/test_regional_exposure.py`, `tests/test_locus_prior.py`.
- [ ] Run synthetic locus sweep + Armis2 benchmarks; document deltas.
- [ ] Refresh goldens (single commit); per-fixture explanation in PR.
- [ ] Expand CHANGELOG entry.

---

## 16. Files Touched (Reference)

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
- [tests/test_scan_payload_fractional.py](../../tests/test_scan_payload_fractional.py) (new)
- [tests/test_density_global.py](../../tests/test_density_global.py) (rewrite)
- [tests/test_regional_exposure.py](../../tests/test_regional_exposure.py) (rewrite)
- [tests/test_locus_prior.py](../../tests/test_locus_prior.py) (extend)
- [tests/test_pipeline_wiring.py](../../tests/test_pipeline_wiring.py) (extend)
- [tests/test_index_format_version.py](../../tests/test_index_format_version.py) (new)
