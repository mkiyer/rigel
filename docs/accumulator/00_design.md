# Fractional Accumulator — Region/Boundary Split Design

Status: design proposal. Drives a native rewrite of the fractional
accumulator. Calibration v6 ([`../caljointmodel/`](../caljointmodel/))
inherits the new substrate.

## 1. Goals

1. **Conceptual clarity.** A reference is partitioned into $N$
   non-overlapping `Region`s and $N+1$ `Boundary`s. Region $i$ borders
   boundary $i$ (left) and boundary $i+1$ (right). The first and last
   boundaries on each reference are terminal (no neighbour); all
   interior boundaries are shared by exactly two regions.
2. **Single-object semantics.** Each boundary is one object owning
   left and right fractional masses + a single integer flux count.
   No more "is `regions[r].boundary_right` redundant with
   `regions[r+1].boundary_left`?" — they were never duplicates, but
   the question keeps coming up because the storage is region-keyed.
3. **Per-fragment accounting.** A fragment contributes exactly 1.0
   total mass to the substrate, distributed across at most one region
   (contained) or a set of boundaries (crossing). Multi-block
   fragments normalize by the **sum of aligned block lengths**, not
   per-block.
4. **Build the calibration substrate in the scan.** The python
   `BoundaryTable` derivation in
   [`../../src/rigel/calibration/boundaries.py`](../../src/rigel/calibration/boundaries.py)
   becomes unnecessary; the native accumulator emits the new layout
   directly.

## 2. Current behavior audit

### 2.1 Storage today

`bam_scanner.cpp` + `accumulator.cpp` write into `RegionCountLedger`
which is region-keyed:

```
per region r:
  contained_{un,sp}_{pos,neg}      // 4 channels
  boundary_left_{un,sp}_{pos,neg}  // 4 channels (mass from fragments that crossed
                                   //               the LEFT edge of region r)
  boundary_right_{un,sp}_{pos,neg} // 4 channels (mass from fragments that crossed
                                   //               the RIGHT edge of region r)
```

Total: 12 float counters per region. A boundary-crossing fragment
deposits mass into **two** region slots (`regions[r].boundary_right`
and `regions[r+1].boundary_left`) — partitioned by the bp-split, not
duplicated. The Python helper
[`boundaries.py`](../../src/rigel/calibration/boundaries.py) rebuilds
a boundary-keyed view from these region slots; it is correct but
costs an extra pass and an extra copy.

### 2.2 Per-block vs per-fragment

**Open finding from this audit.** Need to confirm in code review
that the current accumulator deposits mass per-fragment (with
block-length normalization summing to 1.0) and not per-block.
Existing C++ comments suggest per-fragment but several call sites in
`accumulator.cpp` iterate aligned blocks and call deposit helpers
once per block — risk of double-counting when a block-iteration loop
is not gated by a per-fragment "already-touched" set. This audit is
**part of the implementation work**, not a precondition; the new accumulator makes
per-fragment semantics explicit at the API level so the bug surface
shrinks.

### 2.3 The "fully spans" quirk

When a fragment is longer than its containing region, it **encompasses**
that region — overlapping both its left and right boundaries. The region's
slice mass splits `w/2` into each side (the RIGHT side of its left boundary
and the LEFT side of its right boundary), conserving fragment mass. The new
accumulator handles this with the general per-slice attribution rule (§4.3):
an interior slice crosses two boundaries, so each side gets half — the
`n_cross = 2` branch, no bespoke case.

### 2.4 Why ad-hoc Python boundary derivation must end

[`boundaries.py`](../../src/rigel/calibration/boundaries.py) builds
~250 LoC of derivation logic with `_validate_inputs`,
`build_boundary_table`, and `validate_boundary_table`. The
validation function is a self-test that the derivation is correct —
which it is, but the existence of a "validate this derivation"
function is a smell. Move the substrate into the scan and the
derivation + validator + tests disappear.

## 3. Data structures

```cpp
// src/rigel/native/accumulator.h

struct Region {
  uint32_t contained_unspl_pos;        // contained, unspliced, sense
  uint32_t contained_unspl_neg;        // contained, unspliced, antisense
  uint32_t contained_spl_pos;          // contained, spliced, sense
  uint32_t contained_spl_neg;          // contained, spliced, antisense
};
// 16 bytes per region.

struct Boundary {
  // Fractional mass: left side / right side, broken down by channel.
  float    mass_left_unspl_pos,  mass_left_unspl_neg;
  float    mass_left_spl_pos,    mass_left_spl_neg;
  float    mass_right_unspl_pos, mass_right_unspl_neg;
  float    mass_right_spl_pos,   mass_right_spl_neg;
  // Integer flux: count of fragment-events at this boundary,
  // by channel only (NOT split left/right).
  uint32_t flux_unspl_pos, flux_unspl_neg;
  uint32_t flux_spl_pos,   flux_spl_neg;
};
// 8*4 + 4*4 = 48 bytes per boundary.

struct Accumulator {
  // Per reference: N_r regions, N_r + 1 boundaries.
  std::vector<Region>   regions;       // size = N
  std::vector<Boundary> boundaries;    // size = N + n_refs
};
```

> **Note.** Per-region and per-boundary FL histograms are
> **deliberately omitted** — see audit_phase1.md decision #6. The
> calibrator does not consume FL as a likelihood signal; the
> downstream EM scorer
> ([`../../src/rigel/scoring.py`](../../src/rigel/scoring.py),
> [`../../src/rigel/frag_length_model.py`](../../src/rigel/frag_length_model.py))
> continues to use full-resolution FL per-fragment as it does today.

### 3.1 Region/boundary indexing

For reference $f$ with regions $[r_{f,0}, r_{f,1})$:
- Boundary indices on reference $f$ are $[b_{f,0}, b_{f,1}) = [r_{f,0} + f, r_{f,1} + f + 1)$ (offset by 1 per prior ref).
- For local region $i \in [0, N_f)$: left boundary = $b_{f,0} + i$, right boundary = $b_{f,0} + i + 1$.
- $b_{f,0}$ and $b_{f,1} - 1$ are terminal (no left/right neighbour respectively).

A small `ref_region_offsets[F+1]` and `ref_boundary_offsets[F+1]` map
references to their region/boundary index ranges, parallel to the
existing `RegionArrays.ref_offsets`.

### 3.2 Memory

At $N = 200{,}000$ regions across $\sim 25$ references:
- Regions: 16 × $N$ = **3.2 MB**
- Boundaries: 48 × ($N$ + 25) ≈ **9.6 MB**
- **Total: ~13 MB.** An order of magnitude smaller than the
  legacy `RegionCountLedger` + FL histograms (which were ~75 MB
  before the FL deferral). FL is now consumed exclusively by the
  downstream EM scorer per-fragment at the existing point of use.

## 4. Accumulation algorithm

A fragment is observed as an aligned-block list $B_1, \ldots, B_K$
(sense strand, splice flag). Define $\ell_k$ = aligned bp of block
$k$ (CIGAR M/D/=/X), $L = \sum_k \ell_k$ (total observed aligned bp).
The fragment carries a single strand $s \in \{+, -\}$ and a single
spliced flag (set by the upstream resolver based on whether any
junction in the fragment is recognized as spliced). The strand and
spliced flag are per-fragment; they apply to every deposit the
fragment makes. Per-junction splice classification is **out of
scope** for fractional accumulator — see Phase 1 audit memo decision #2.

### 4.1 Decide: contained vs crossing

Per-fragment classification:

```
region_of(B_k) = unique region containing all of block B_k's bp
                 (or -1 if B_k straddles a boundary — rare; see §4.4)
fragment_regions = { region_of(B_k) for k }
```

- If `len(fragment_regions) == 1` and no block straddles a boundary →
  **contained**. The fragment increments exactly one region's count.
- Otherwise → **crossing**. The fragment deposits mass into one or
  more boundaries; the regions are not touched.

### 4.2 Contained: integer increment

A contained fragment with strand $s$ and any splice flag adds
**exactly +1** to the corresponding contained channel of its region:

```cpp
regions[r].contained_{spl|unspl}_{pos|neg} += 1;
```

Block count and block lengths are irrelevant — the fragment is one
biological observation. The integer-only design is sufficient because
$\sum_\text{mass contributions for one contained fragment} = 1$ exactly.

A spliced fragment is one whose CIGAR contains an N (or whose
block-pair gap crosses an annotated junction) **and** stays within
one region (which can happen for short introns or wide regions).

### 4.3 Crossing: per-slice boundary attribution

A fragment not contained in one region is decomposed into per-region
**slices** $s_1, \ldots, s_M$ (one per region it touches, in genomic order;
§4.5.4). Each slice deposits its mass onto the boundaries it crosses,
**conserving the fragment's mass**.

**Boundary side convention (the foundation).** A boundary $B$ separates
two regions. `B.mass_left` is the mass contributed by the slice in the
region **to the left of** $B$; `B.mass_right` by the slice **to the right
of** $B$. `mass_left` is *not* "the left half of $B$"; it is *the side
from which the supporting slice came*.

**Attribution rule.** Slice $s_i$ (in region $R_i$, length $\ell_i$):

1. It crosses its **left** boundary iff $i > 1$ (a slice precedes it) and
   its **right** boundary iff $i < M$ (a slice follows it). Let
   $n_i = [\,i>1\,] + [\,i<M\,] \in \{1, 2\}$ be its crossing count; its
   per-crossing mass share is $\ell_i / (L\,n_i)$.
2. If it crosses right, deposit the share into
   $(\text{right boundary of } R_i).\text{mass\_left}$ (the slice is to the
   *left* of that boundary) and increment that boundary's `flux_left` by 1.
3. If it crosses left, deposit the share into
   $(\text{left boundary of } R_i).\text{mass\_right}$ (the slice is to the
   *right* of that boundary) and increment that boundary's `flux_right` by 1.

**Encompassed regions split 50/50.** An *interior* slice ($1 < i < M$) — a
region the fragment fully traverses, overlapping **both** its boundaries —
has $n_i = 2$, so each of its two boundary sides receives **half** its mass.
End slices ($i = 1$ or $i = M$) have $n_i = 1$ and deposit full mass on their
one crossed side. This is what makes per-fragment mass conserve (§6): the sum
over all slice-sides is $\sum_i \ell_i / L = 1$.

**Adjacent vs. spliced — no special case.** If $R_i, R_{i+1}$ are adjacent,
$s_i$'s right boundary *is* $s_{i+1}$'s left boundary — the one shared
boundary receives `mass_left` from $s_i$ and `mass_right` from $s_{i+1}$ (a
contiguous crossing, both sides credited). If non-adjacent (a spliced
intron-skip), the interior boundaries between them get **nothing**: a skipped
region produces no slice, so it is never an $s_i$ — the **only observed
events** principle. Flux is the integer per-side count (§4.4), never split.

Mnemonic: slice-on-left → `mass_left`; slice-on-right → `mass_right`.
Never the other way around.

### 4.4 Flux increment

Each junction event also increments the **integer flux** counter of
each touched boundary, by channel:

```cpp
boundaries[right_boundary_of(R_k)].flux_{spl|unspl}_{pos|neg} += 1;
boundaries[left_boundary_of(R_{k+1})].flux_{spl|unspl}_{pos|neg} += 1;
```

If the two boundaries are identical (adjacent regions), the single
boundary gets +1 (not +2): a junction event is **one** observation
per boundary. Implementation: deduplicate by boundary index before
incrementing flux.

Per fragment-with-multiple-junctions: each junction is an independent
event, so the middle boundary in a 3-block fragment may get two
flux increments (one for each junction adjacent to it). Per the
worked example in §5.2, this is the intended semantics.

### 4.5 Worked examples

#### 4.5.1 Two-block, non-adjacent gap (user-verified walkthrough)

Exons (1000,2000), (5000,6000). Region partition includes exon-1
$R_1 = (1000,2000)$, intronic gap $R_2 = (2000,5000)$, exon-2
$R_3 = (5000,6000)$. Boundaries at positions 1000, 2000, 5000, 6000
(call them $B_{1000}, B_{2000}, B_{5000}, B_{6000}$).

Fragment with blocks $B_1 = [1800, 1950)$ ($\ell_1 = 150$) in $R_1$,
$B_2 = [5050, 5950)$ ($\ell_2 = 900$) in $R_3$. $L = 1050$. Strand
$+$, implicit splice between block-1 and block-2.

Applying the §4.3 rule:

- Block-1 lies in $R_1$ and implicitly exits via $R_1$'s right
  boundary = $B_{2000}$. Block-1 is to the **left** of $B_{2000}$ →
  $B_{2000}.\text{mass\_left}.\text{spl\_pos} \mathrel{+}= 150/1050$;
  $B_{2000}.\text{flux.spl\_pos} \mathrel{+}= 1$.
- Block-2 lies in $R_3$ and implicitly enters via $R_3$'s left
  boundary = $B_{5000}$. Block-2 is to the **right** of $B_{5000}$ →
  $B_{5000}.\text{mass\_right}.\text{spl\_pos} \mathrel{+}= 900/1050$;
  $B_{5000}.\text{flux.spl\_pos} \mathrel{+}= 1$.
- Intronic gap region $R_2 = (2000, 5000)$ has no block support —
  skipped. No contained increment.
- $B_{2000}.\text{mass\_right}$ and $B_{5000}.\text{mass\_left}$ are
  **not** touched: there is no block in $R_2$ to support them.
- Terminal-side boundaries $B_{1000}$ and $B_{6000}$: unchanged.

Mass deposited: $150/1050 + 900/1050 = 1.0$. ✓

#### 4.5.2 Three-block, all-adjacent-exon (middle slice splits 50/50)

Three exonic regions $R_1, R_2, R_3$ with shared interior boundaries
$B_a$ (between $R_1$ and $R_2$) and $B_b$ (between $R_2$ and $R_3$).
Fragment slices $s_1$ in $R_1$ ($\ell_1 = 100$), $s_2$ in $R_2$
($\ell_2 = 80$), $s_3$ in $R_3$ ($\ell_3 = 120$). $L = 300$. Strand
$+$, both inter-block junctions are spliced.

- $s_1$ (end slice, $n_1 = 1$): crosses right only → $B_a.\text{mass\_left}.\text{spl\_pos} \mathrel{+}= 100/300$; $B_a.\text{flux\_left.spl\_pos} \mathrel{+}= 1$.
- $s_2$ ($R_2$ is **encompassed**, $n_2 = 2$): each side gets half → $B_a.\text{mass\_right}.\text{spl\_pos} \mathrel{+}= 40/300$ ($B_a.\text{flux\_right}\mathrel{+}=1$) and $B_b.\text{mass\_left}.\text{spl\_pos} \mathrel{+}= 40/300$ ($B_b.\text{flux\_left}\mathrel{+}=1$).
- $s_3$ (end slice, $n_3 = 1$): crosses left only → $B_b.\text{mass\_right}.\text{spl\_pos} \mathrel{+}= 120/300$; $B_b.\text{flux\_right.spl\_pos} \mathrel{+}= 1$.

Total mass: $(100 + 40 + 40 + 120)/300 = 300/300 = 1.0$. ✓ The middle
slice's mass ($80/300$) splits 50/50 between its two boundary sides, so the
fragment conserves mass (§6). Flux is unaffected — each boundary side still
records one integer crossing event.

#### 4.5.3 Single-region, multi-block (contained spliced fragment)

Long exonic region containing both blocks of a spliced fragment.
$\ell_1 = 60$, $\ell_2 = 40$, $L = 100$. Strand $+$, spliced.

`regions[r].contained_spl_pos += 1`. No boundary touched. No
per-block accounting. The fragment is one biological observation.

#### 4.5.4 Fully-spans-region (the encompassing case)

Fragment with one block straddling a small region $R_k$ entirely
(implies fragment longer than the region). The block decomposes into
slices $\ell^{(k-1)}, \ell^{(k)}, \ell^{(k+1)}$ across $R_{k-1}, R_k,
R_{k+1}$; $L = $ their sum. $R_k$ is **encompassed** — the fragment
overlaps both its boundaries.

Apply the §4.3 per-slice rule: the two end slices ($\ell^{(k-1)},
\ell^{(k+1)}$) deposit full mass on their single crossed side; the
interior slice $\ell^{(k)}$ has $n = 2$ and splits 50/50 —
$\ell^{(k)}/(2L)$ to the RIGHT side of $R_k$'s left boundary and
$\ell^{(k)}/(2L)$ to the LEFT side of its right boundary. Total mass
conserves to $1.0$.

This is the only case where a single CIGAR-block straddles region
boundaries; it falls out of the per-slice rule with no special case.

## 5. Reserved

*Reserved.* Earlier drafts placed FL histograms here; see
`audit_phase1.md` decision #6 — the calibration substrate does not
accumulate FL. Section number preserved so cross-references remain
stable.

## 6. Mass conservation invariants

> **Per-fragment mass conservation.** Every fragment deposits total mass
> exactly $1.0$ — across `contained` (single-region fragment) or across the
> boundary `mass_left`/`mass_right` slots (crossing fragment). A region the
> fragment **encompasses** (an interior slice, §4.3) splits its slice mass
> 50/50 between its two boundary sides, so no mass is double-counted. This
> restores the legacy `w/2` behavior (§2.3) and the
> $\int \text{overlap}/\ell\,dx = L$ identity the calibration exposure relies on.

Aggregate invariants:

1. **Per-fragment mass = 1.0.** For a fragment with slices $s_1, \ldots, s_M$
   (after the §4.5.4 block decomposition), each slice $s_i$ deposits
   $\ell_{s_i}/L$ total, split across the $n_i \in \{1,2\}$ boundary sides it
   crosses. The sum is
   $$\frac{1}{L} \sum_{i=1}^{M} \ell_{s_i} \;=\; 1.0$$
   for every $M$ — contained ($M=1$), simple crossing ($M=2$), and spanning
   ($M \geq 3$) alike.
2. **Contained counts.** $\sum_r \mathrm{contained}_r = N_\text{contained}$
   (each contained fragment increments exactly one region by +1).
3. **Flux totals.** $\sum_b (\mathrm{flux\_left} + \mathrm{flux\_right})_{b,c}$
   counts the (slice, crossed-boundary-side) events in channel $c$. Flux is
   the **integer** per-side crossing count and is *not* split by the 50/50
   rule — only mass is. (Power from flux, density from mass.)
4. **No terminal-boundary deposits in typical references.** Terminal
   boundaries receive mass only when fragments extend beyond a region
   partition's edge — zero for well-formed reference partitions.

**Implication for calibration v6.** Per-region mass now conserves, so the
substrate's $M^{(\text{cont})} + M^{(L)} + M^{(R)}$ is a faithful per-region
fragment-mass total (D1). Statistical *power* still comes from the integer
flux; fractional *mass* from these conserved deposits.

## 7. Public API (native + python)

### 7.1 C++ surface

```cpp
// accumulator.h
struct Accumulator {
  std::vector<Region>   regions;
  std::vector<Boundary> boundaries;
  std::vector<int64_t>  ref_region_offsets; // [F+1]
  std::vector<int64_t>  ref_boundary_offsets;

  void deposit(const FragmentEvent& fe);
};
```

`bam_scanner.cpp` is unchanged in its outer loop; the per-fragment
deposit path swaps `RegionCountLedger.deposit_*` for
`Accumulator.deposit`.

### 7.2 nanobind binding

Exposes the same fields as numpy arrays (zero-copy or with explicit
copy depending on lifetime). Replaces today's
`RegionCountLedger`/`BoundaryTable` pair.

### 7.3 Python adapter

```python
@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    # Region-keyed
    region_contained_unspl_pos: np.ndarray  # uint32[N]
    region_contained_unspl_neg: np.ndarray
    region_contained_spl_pos:   np.ndarray
    region_contained_spl_neg:   np.ndarray
    region_fl_hist:             np.ndarray  # uint32[N, FL_BINS]

    # Boundary-keyed
    boundary_mass_left_unspl_pos:  np.ndarray  # float32[B]
    boundary_mass_left_unspl_neg:  np.ndarray
    boundary_mass_left_spl_pos:    np.ndarray
    boundary_mass_left_spl_neg:    np.ndarray
    boundary_mass_right_unspl_pos: np.ndarray
    boundary_mass_right_unspl_neg: np.ndarray
    boundary_mass_right_spl_pos:   np.ndarray
    boundary_mass_right_spl_neg:   np.ndarray
    boundary_flux_unspl_pos: np.ndarray  # uint32[B]
    boundary_flux_unspl_neg: np.ndarray
    boundary_flux_spl_pos:   np.ndarray
    boundary_flux_spl_neg:   np.ndarray
    boundary_fl_hist:        np.ndarray  # uint32[B, FL_BINS]

    # Topology
    ref_region_offsets:   np.ndarray  # int64[F+1]
    ref_boundary_offsets: np.ndarray  # int64[F+1]
```

The `compute_locus_priors_from_partitions` consumer (and the
calibration v6 calibrator) read directly from `AccumulatorPayload`.
No intermediate `BoundaryTable` derivation.

## 8. Implementation phases

1. **Audit & spec lock.** Read `accumulator.cpp` line-by-line and
   confirm the per-fragment normalization semantics; fix any
   per-block double-counts found. Write the per-fragment deposit
   reference test (synthetic fragments → expected `Accumulator`
   contents).
2. **Native rewrite.** Implement `accumulator.h` + `.cpp`.
   `bam_scanner.cpp` switches its inner deposit calls. Keep the old
   `RegionCountLedger` compiled-but-unused for one commit so tests
   can compare.
3. **nanobind + python adapter.** Expose `AccumulatorPayload`.
   Update `scan_payload.py`.
4. **Consumer migration.** Rewrite `boundaries.py` consumers to read
   from `AccumulatorPayload` directly. Delete `boundaries.py` and
   its tests.
5. **Calibration substrate.** The calibrator's substrate becomes a
   thin view over `AccumulatorPayload` (native, not a Python extension).
6. **Cleanup.** Delete `RegionCountLedger`-only fields. Update
   `CLAUDE.md` and `.github/copilot-instructions.md`.

## 9. Test plan

Per-fragment deposit tests in `tests/native/test_accumulator.py`:

| Test | Fragment | Expected deposit |
|---|---|---|
| `test_contained_single_block` | block in $R_5$ | `regions[5].contained_unspl_pos += 1`; nothing else |
| `test_contained_multi_block_spliced` | 2 blocks both in $R_5$ | `regions[5].contained_spl_pos += 1`; nothing else |
| `test_two_block_adjacent_regions` | blocks in $R_5$, $R_6$ (boundaries 5, 6, 7) | `boundaries[6].mass_left += $\ell_1/L$; mass_right += $\ell_2/L$; flux += 1` |
| `test_two_block_non_adjacent_regions` | blocks in $R_5$, $R_8$ | `boundaries[6].mass_left += $\ell_1/L$; boundaries[8].mass_right += $\ell_2/L$; intermediate boundaries unchanged` |
| `test_three_block_all_adjacent` | spec example §4.5.2 | per spec |
| `test_fully_spans_region` | one block straddling $R_4, R_5, R_6$ | per spec §4.5.4 |
| `test_mass_conservation` | random N=10000 fragments | total mass == N exactly |
| `test_flux_dedup` | adjacent-region 2-block fragment | one flux increment per touched boundary, not two |

Plus golden tests preserved against today's `RegionCountLedger`
behavior on existing scenario BAMs — `test_accumulator_matches_legacy.py`
runs both accumulators on the same fragment stream and asserts that
the boundary view, when projected back to region-keyed slots,
matches the legacy `boundary_left`/`boundary_right` to bit-exact.

## 10. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Hidden per-block double-counts in legacy `accumulator.cpp` discovered during §8 audit | Medium | Medium | Fix during the rewrite; document each bug fixed in CHANGELOG |
| Boundary index math off-by-one across reference seams | Low | High | Property tests over multi-reference synthetic fixtures |
| Memory pressure on very large reference sets (high $N$) | Low | Low | Memory is smaller than today's layout; bounded |
| Downstream `boundaries.py` consumers have hidden dependencies on the validator | Low | Low | Audit `validate_boundary_table` callsites before deletion |
| Implicit-splice classification (which inter-block gaps are "spliced") inconsistent with legacy | Medium | Medium | Reuse the existing classifier code path; do not reinvent here |

## 11. Open questions deferred to implementation

1. **Block-straddles-boundary case detail.** §4.5.4 sketches the
   rule but the C++ implementation will need to enumerate "virtual
   block slices" per region intersection. Decide whether to do this
   in the deposit function or in a pre-processing step.
2. **Implicit-splice detection.** Currently inferred from CIGAR N or
   from inter-block-gap-crosses-annotated-junction. The accumulator
   This design inherits this — no new rule.
3. **Strand assignment for ambiguous fragments.** Inherits today's
   logic from `bam_scanner.cpp`.
4. **Soft-clip handling.** Soft-clipped bases excluded from $\ell_k$
   (consistent with current `read_length` definition per
   [`../../CLAUDE.md`](../../CLAUDE.md)).
