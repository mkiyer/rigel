# Accumulator — implementation plan

This document says *what code to write*, in what order, with real names and line numbers. Line numbers
are from `4bb4d191`.

**Read in this order:** this state block · `ACCUMULATOR_DESIGN.md` (the design) · `LEDGER.md` (what has
landed, and why each thing is the way it is) · `CARRY_FORWARD.md` (measured facts and traps from the
project's deleted documentation).

---

## 0. STATE — start here

| | |
|---|---|
| HEAD | `4a7a78d7`. S1, S2 and the doc corrections are **committed**; the tree is clean |
| done | **S1** (index: reach + junction CSR) · **S2** (reference accumulator, spec matrix, scan profiler, real-data shim, adversarial review, two owner rulings) |
| next | ⭐ **S3 — the C++** (§3.1–§3.4). The doc corrections are done |
| suite | **1303 pass**, 53 in the accumulator spec matrix; ruff + `ruff format` clean |
| bench | r0 `0.079005` / r3 `0.046675`, 32/32 flat — unchanged since `3c293038` |
| deposit budget | **~357 ns/fragment** end-to-end (`scripts/design/scan_profile.py`), plus ~0.108 s fixed for the 1.04 M-node partition |

⚠ Line numbers in §1–§3 are from `4bb4d191` and are still valid: nothing committed since has touched
`src/rigel/native/`. S3 is the first step that does.

### ⭐ NAMING — one word per concept, and they are the codebase's own words

Fixed 2026-07-29 after the owner rejected three inventions. Recorded because the C++ will follow these.

| ⛔ invented | ✅ use | why |
|---|---|---|
| `genome_channel`, `CHANNEL_PLUS/MINUS`, `kNChannels` | `STRAND_COLUMNS`, `Strand::POS`/`NEG` | "channel" was leftover vocabulary from the 4-way `spliced × primary` axis being deleted. The axis IS the strand; say strand |
| `N_STRANDS = 2` | `N_STRAND_COLUMNS = 2` | there are **four** strands (`Strand` has OR semantics: `POS\|NEG == AMBIGUOUS`, plus `NONE`). Only two name a column |
| `motif_strand`, then `sj_motif_strand` | **`sj_strand`** | `sj` is the codebase's abbreviation already — `sj.feather`, `sj_map`, `SJKey`, `sj_lookup_into`, `cr.sj_strand`. "motif" does not say *which* motif, and a third name for one quantity is worse than an ambiguous one |
| `lo` / `hi` for a fragment | **`start` / `end`** | the codebase uses `start`/`end` everywhere for coordinates; `lo`/`hi` appears only in `reach_lo_pos`/`reach_hi_pos`, which is a different quantity |

⚠ **The shipped C++ is where `motif_strand` came from**: `bam_scanner.cpp:1468` declares
`int32_t motif_strand = cr.sj_strand;` — a local alias whose name is worse than the field it copies. S3
rewrites that adapter, so **delete the alias and use `cr.sj_strand` directly.**

⚠ A junction edge's **annotated** strand and a fragment's **observed** strand are the same quantity from two
sources, so both are `sj_strand`, disambiguated by scope. The comparison then reads as what it is:
`if p.sj_strand[k] != sj_strand` — annotated versus observed.

### ✅ THE DOC CORRECTIONS — LANDED (2026-07-29)

The plan and the design are what the C++ follows, so they were corrected before any C++ was written. Four
were known and listed; **three more of the same class were found while making them.** All seven are
documentation only — no code changed, so the suite (1293 at the time) and the bench were untouched by
construction. S2.2 then added 5 tests, taking the suite to 1298.

| | correction | where |
|---|---|---|
| 1 | **§5.1 rewritten for the ONE strand convention.** Was `spliced × {sense, antisense}`; now `{genome +, genome −}` in every bank without exception, with sense/antisense stated as **derived, never stored**. Records what the channel physically is: the fragment's **read-1 genome orientation** (`bam_scanner.cpp:661-686`), which makes declared strandedness a header field rather than a channel convention | design §5.1, §5.2 |
| 2 | **§3.2's pseudocode replaced by a transcription of the reference**, with the four bugs it used to encode listed and each attributed: `crossed == 0` for CONTAINED · `L` from the **unclipped** span by a second formula · `lo`/`hi−1` instead of the path's first and last **covered** base · unguarded `quantum(L−1)`. Plus the ordering constraint the old block violated: the junction lookup must run **before** the crossing loop, because `spliced` picks the edge bank | plan §3.2 |
| 3 | **§3.3 — `junction_edge_id_` deleted, not documented.** The id **is** the CSR slot. ⭐ Measured on both indexes: the slot equals the dense rank within `flatnonzero(kind == JUNCTION)` for all 404,168 (and 537), so `edges_df.edge_row` is purely a payload→index join key and **never crosses the ABI** | plan §3.3 |
| 4 | **The `L == 1` guard stated, with its residue**: an `L = 1` path using an annotated junction books a junction count against density 0 — the schema's only count/density co-support violation, and correct | plan §3.2 bug 4 |
| **5** | ⭐ **found while making the others — §8 said "two structurally pure pools", §8.2 then said "a third", and the settled ruling is five.** §8 now carries the five-pool table with its measured means (gDNA 88.0/88.5, splash 138.8/211.6, RNA 220.7) and the note that this moves the shipped gDNA mean length by ~40 % | design §8, §8.2 |
| **6** | ⭐ **found while making the others — the plan's S3 gate said "no regression against the 124 ns baseline", but the ledger's S3 budget is ~357 ns/fragment end-to-end.** Two different measurements of two different things, one of them named as the gate in two places with different values. §6 now distinguishes them explicitly and records that S2 built the harness (§6 still claimed none existed) | plan §5 S3 row, §6 |
| **7** | ⛔ **found while making the others, and the worst of the three — `CARRY_FORWARD.md` §4 is a PRE-SETTLEMENT idea list, and six of its bullets contradict the settled deposit rule.** It told a C++ author to weight an edge `1/L`, to treat "crossed nothing" as contained, that `L` is "exonic bases only", to deposit unannotated junctions on the **spliced** channel, to implement three conservation identities the design rejected as tautologies, and to leave a 1 bp spanned node untouched. §4 now opens with a banner and each item is marked `⛔ SUPERSEDED` in place with what replaced it | `CARRY_FORWARD.md` §4 |

Also landed in the design: **§12.A**, which called the minimum intron length the one item blocking
implementation, now records both owner rulings — **every CIGAR `N` is an intron, no floor** (a floor would
be a magic number, and no floor is what ships) and **a zero-length exon is impossible, so ABUTTING introns
are malformed and merge** (reversing the strict-overlap decision taken during the S2 review). The artifact
catalogue (§3.3) gained rows for overlapping and abutting introns; §7.2's dangling `§12.A` pointer now
points at the calibration list where the capture-model question actually lives. **Nothing blocks
implementation.**

### TODO — deferred work, none of it blocking S3

Each item is deferred deliberately, with the reason. This is the project's TODO list; add to it rather than
starting a new file.

* ⛔ **`detect_chimera` is blind to two real populations, because of one gate.** It considers only blocks
  with a **non-empty transcript set** (`constants.h:339-343`) and needs ≥ 2 mutually disjoint transcript
  components, so it returns `CHIMERA_NONE` for both of these:
  1. **Same-reference-orientation mates.** Paired-end sequencing reads the two ends of one double-stranded
     molecule, so R1 and R2 must map to **opposite** reference strands; `build_fragment` encodes that by
     flipping R2, so a normal FR pair collapses to one `(ref, strand)` group. Two groups means the mates
     agree in reference orientation, which is aberrant — consistent with an inversion placing both ends on
     the same reference strand. It is evidence of a rearrangement, not a molecule. Today it becomes
     `align_strand = STRAND_AMBIGUOUS` (`POS|NEG`) and is **silently dropped** at `bam_scanner.cpp:1474-1480`.
  2. **The multi-megabase spans.** ⭐ 95 % carry a supplementary record, and intergenic blocks have empty
     transcript sets, so the gate never sees them. Same root cause, second symptom.

  The vocabulary already exists: `CHIMERA_CIS_STRAND_SAME` / `CHIMERA_CIS_STRAND_DIFF`.
  **Why deferred:** widening the gate reclassifies currently-accepted fragments as chimeras, which moves the
  bench. That is a change to *what counts as a fragment*, not to how one is tallied — so it is its own arm
  with its own before/after measurement, after S3 is byte-identical. Owner-agreed 2026-07-29.
* **The 8 indexes on disk carry stale `reach`.** S1 changed the builder, not the artifacts. Nothing reads
  those columns so nothing is wrong today, but every index must be rebuilt before S5/S6 makes reach
  load-bearing. `partition_hash` will not notice — it covers `nodes.feather` only.
* **The shipped C++ accumulator has no reference test** until S3 replaces it (the old spec matrix tested
  the fractional rule, which no longer exists). The oracle bench is the guard in that window.
* **`sj.feather` duplicates junction coordinates.** It is one row per *(junction, transcript)* and feeds
  `sj_map` at load (`index.py:1351`), so it cannot merge into `edges.feather` (one row per distinct
  junction). But its `(ref, start, end, strand)` columns duplicate what `(src, dst, strand)` already
  determines, and nothing enforces they agree. **Available simplification: re-key it to
  `(junction_edge_id, t_index)`.** Owner's idea; log it for after the accumulator lands.
* **Goldens moved at the index change** (`gdna_em_count` fell 16–52 %) and that movement has never been
  validated against truth. The new benchmark suite should adjudicate it, not this sequence.

---

## 1. The verdict: refactor in place, do not replace

`accumulator.cpp` is 356 lines in eight functions. Only one of them is genuinely rewritten.

| function | line | fate |
|---|---|---|
| `Accumulator::Accumulator` | 30 | **keep** — already takes a sorted cut array; cuts are cuts |
| `Accumulator::region_of_pos` | 74 | **keep**, rename `node_of_pos`. `upper_bound − 1`, byte-compatible with `np.searchsorted(side='right')−1` |
| **`Accumulator::deposit`** | **86–232** | **rewrite** — the only substantial one |
| `Accumulator::merge_from` | 234 | adapt — same element-wise shape, new fields |
| `AccumulatorSet::AccumulatorSet` | 275 | **keep** — per-reference slicing of the flat partition |
| `AccumulatorSet::at` × 2 | 328, 337 | **keep** |
| `AccumulatorSet::merge_from` | 346 | **keep** |

The 60 % that survives is the part that is easy to get subtly wrong — offset arithmetic and the
per-worker merge, which is exactly where a ref-id mismatch once silently dropped 476,719 of 476,732
fragments while every golden test passed. **No file is created and none is deleted wholesale.** Same
files, same class names, same binding shape, no version suffixes anywhere.

## 2. What the tree already gives us

Four findings that shrink the work substantially. Each is the reason a step below is small.

**2.1 The scanner already computes the design's path.**
`fragment_genomic_spans` (`bam_scanner.cpp:750-794`) returns, per reference, `[min block start, max
block end]` **with the mate gap filled and the cut introns excised**. That is §3.1's path segments
exactly, and `Σ segment lengths` is already `span − Σ intron` = **`L`**. No new span code is needed.

**2.2 Calibration's read surface is four call sites.** Everything calibration knows about the payload
enters through:

| | |
|---|---|
| `CalibrationSubstrate.from_payload` | `substrate.py:109, 117-120` |
| `BoundarySubstrate.from_payload` | `substrate.py:195-200` |
| `build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)` | `calibrate.py:305` |
| `gdna_fl_mass(payload)` | `fl.py:130` |

Nothing else in `src/` touches a payload field, and no consumer names one as a string. So the rewiring
is a boundary-layer rewrite plus a cascade through **one** derived type, not a diffuse edit.

**2.3 The reach columns have no consumer in `src/`.** Grep-verified: `reach_lo_pos` and friends appear
only in `splice_graph.py` and `tests/calibration/test_splice_graph.py`. Changing their meaning (§P1)
breaks five tests and nothing else. **Cheap now, expensive later.**

**2.4 The contiguous rows of `edges.feather` are redundant.** `kind == 0` always has `dst == src + 1`
and `strand == 0`; there are exactly `n_nodes − n_refs` of them. The accumulator derives contiguous
edges from the cut array it already receives. Only the **junction** rows carry information the
accumulator needs.

## 3. The target code

### 3.1 C++ types — `accumulator.h`

`Region` and `Boundary` are replaced. `channel_idx`, `kNChannels` and the FL-pool block go.

```cpp
namespace rigel::accumulator {

//: The two array columns ARE the two genome strands, Strand::POS and Strand::NEG. `Strand` has FOUR
//: values (OR semantics make POS|NEG == AMBIGUOUS, and NONE means no strand), so only these two name a
//: column and a fragment carrying neither is REJECTED, never filed under one of them.
inline constexpr std::size_t kNStrandColumns = 2;

//: Densities accumulate as round(kDensityScale / placements). Integer addition is associative, so the
//: merge is bit-identical at any worker count -- which the float32 mass channels are not.
inline constexpr std::uint64_t kDensityScale = 1ull << 32;

inline std::uint64_t density_quantum(std::int64_t placements) noexcept;   // half away from zero

struct Node {                                   // 48 B exactly, no padding
    std::uint32_t contained_count[kNStrandColumns];
    std::uint32_t spanning_count[kNStrandColumns];
    std::uint64_t contained_density[kNStrandColumns];
    std::uint64_t spanning_density[kNStrandColumns];
};

struct ContiguousEdge {                         // 48 B
    std::uint32_t unspliced_count[kNStrandColumns];
    std::uint32_t spliced_count[kNStrandColumns];
    std::uint64_t unspliced_density[kNStrandColumns];
    std::uint64_t spliced_density[kNStrandColumns];
};

struct JunctionEdge {                           // 24 B
    std::uint32_t count[kNStrandColumns];
    std::uint64_t density[kNStrandColumns];
};
```

`static_assert` each size, as the current header does for `Region`/`Boundary`.

**`start_count` lives in its own flat array**, `std::vector<std::uint32_t> node_start_count_`, not in
`Node`. Three reasons: `Node` then measures 48 B *exactly* with no padding member; design §10.2's invariant
array becomes a single `memcpy` at build-result time instead of a strided gather; and the hot struct
carries only what the deposit's inner loops touch.

**Per-worker memory at human scale:** `1,043,881 × 48 + 1,043,595 × 48 + 404,168 × 24 + 1,043,881 × 4`
= **110 MB**, against **~92 MB** today. Both are per worker and both are dwarfed by the 8.6 GB peak RSS.
Not a blocker; recorded so nobody is surprised.

### 3.2 The deposit — allocation-free

The current body builds a `std::vector<Slice>` per fragment (`accumulator.cpp:105-106`). ⭐ Measured
this session: that allocation is **22.8 ns, 18 % of the 124 ns per-fragment deposit**, and it is
invisible to any profiler that samples by function because it is attributed to `malloc`. The new body
allocates nothing.

```cpp
struct FragmentPath {
    std::int32_t        ref_id;
    const std::int64_t* block_starts;   // merged aligned blocks, sorted, one reference
    const std::int64_t* block_ends;
    std::size_t         n_blocks;
    const IntronBlock*  introns;        // deduplicated on (ref, start, end), strands OR-ed
    std::size_t         n_introns;
    std::int32_t        align_strand;   // cr.align_strand -- the READ-1 genome orientation. THE channel
    std::int32_t        sj_strand;   // cr.sj_strand. ⚠ NEVER selects a channel; §3.3 disambiguation only
    bool                introns_inferred;   // SPLICE_IMPLICIT -- bars it from the pure-RNA pool
};

DepositOutcome deposit(const FragmentPath& path);   // outcome: DEPOSITED / TOO_LONG / EMPTY
```

The algorithm, per fragment. ⚠ **This is a transcription of the reference's `deposit`, which is the
specification. If the two ever disagree, the reference wins and this block is the thing that is wrong.**

```
// ── clip to the reference, then derive L exactly ONCE ─────────────────────────────────────
cuts = cut_positions_ + ref_cut_offsets_[ref];  n = cuts on this reference
if n < 2:                      return EMPTY                  // reference contributes no nodes
lo = max(path.lo, cuts[0]);  hi = min(path.hi, cuts[n-1])    // ⚠ CLIP FIRST, then measure
if hi <= lo:                   return EMPTY

introns  = normalise(path.introns, lo, hi)   // sorted, DISJOINT, clipped. Overlapping AND ABUTTING
                                             //   introns merge; each merge bumps introns_absorbed
segments = [lo, hi) with those introns cut out
L        = Σ (b - a) over segments           // ⭐ L is DEFINED here. There is no second formula.
if L <= 0:                     return EMPTY
if L > max_fragment_length_:    return TOO_LONG              // a QC counter, not a silent drop

// ── which annotated junctions does the path use? this also selects the edge bank ──────────
for each intron:  jid = sj_edge_id(ref, intron, path.sj_strand)   // §3.3; -1 if unannotated
spliced = (at least one jid >= 0)
unannotated_introns_ += n_introns - n_annotated               // deposits nothing across the gap (design §4.1)

channel = (path.align_strand == STRAND_POS) ? kChannelPlus : kChannelMinus   // ⭐ genome strand, always

// ⚠ the path's own first and last COVERED base -- NOT lo and hi-1
first_base = segments.front().a;   last_base = segments.back().b - 1

node_start_count_[node_base + node_of(cuts, first_base)] += 1     // always, exactly once

// ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────
q_edge = (L >= 2) ? density_quantum(L - 1) : 0     // ⚠ the L == 1 guard; see bug 4 below
q_node = density_quantum(L)
edge   = spliced ? spliced_bank : unspliced_bank
for each segment [a, b):
    j0 = upper_bound(cuts, a) - cuts       // first line strictly inside the segment
    j1 = lower_bound(cuts, b) - cuts       // first line at or past its end
    for j in [j0, j1):      edge[edge_base + j - 1]      += (1, q_edge)
    for j in [j0, j1 - 1):  node_spanning[node_base + j] += (1, q_node)   // interior of THIS range
    n_crossed += j1 - j0

for each jid >= 0:  junction[jid] += (1, q_edge)

// ── contained: the WHOLE path lies inside ONE node ────────────────────────────────────────
if no annotated junction and node_of(first_base) == node_of(last_base):
    node_contained[node_base + node_of(first_base)] += (1, q_node)

pool_lengths_[pool(spliced, inferred, contained_node, sole_line)][L] += 1     // §3.4
```

Two binary searches per segment, then straight loops. No allocation, no set, no hash.

⚠ **The junction lookup runs BEFORE the crossing loop**, because `spliced` chooses which of the two edge
banks the crossings land in. This block used to resolve junctions afterwards.

**⛔ The four things this block got wrong before S2, each a real bug the reference now pins.** They are
listed because a C++ author following the old text would have reproduced all four — and byte-identity, the
only gate S3 has, would then have been byte-identity to the wrong answer.

1. **CONTAINED was `crossed == 0`.** It must be *the whole path lies inside one node*. An unannotated
   intron can swallow every line between two blocks, leaving a fragment that crosses **nothing** yet
   straddles two nodes; the old rule credited its whole length to the first of them. Such a fragment
   deposits on **no object at all** but still increments `node_start_count`, so the loss is visible
   rather than silent.
2. **`L` came from the UNCLIPPED span**, by a second formula (`(hi−lo) − Σ intron`) contradicting the
   segment total the code builds anyway. ⭐ Measured: the two disagree by up to 1.5× as soon as two
   introns overlap, and a wide overlap drives `L` **negative** (−290 against a true 20), filing a good
   fragment as `dropped_empty`. Live on real data at 1 in 875,670 — and blocking far beyond that rate,
   because a C++ author builds the segments first and will naturally sum them.
3. **The contained test and `start_count` used `lo` / `hi−1`**, the fragment's extent. With a leading
   intron the molecule does not begin at `lo`: `introns=[(150,480)]` over `[150,500)` has its whole path
   in `[480,500)`, a different node — so the invariant was credited to a node the fragment never touched.
4. **`quantum(L−1)` was called unguarded** — integer division by zero, reachable from
   `deposit(ref, 500, 501)`: `L = 1`, accepted, contained. A length-1 molecule cannot cross a 0-bp line,
   so `q_edge` is **0** there. ⚠ **Its residue is the schema's only count/density co-support violation**:
   an `L = 1` path using an annotated junction books a junction *count* against density 0. That is
   correct — one observation, zero measured density — but a consumer dividing count by density must not
   assume the two share a support.

⚠ **"Spanned by the *same* segment" is a correction to the design's wording**, and it is load-bearing.
Under "both lines crossed" a spliced fragment would be counted as spanning a node it *jumps over*. It
is also the cheaper rule — the spanned nodes are simply the interior of each crossed range. It stays
self-consistent because a node lying wholly inside a segment has equal genomic and transcript length,
so the `(L − ℓ − 1)` placement count is right in both spaces.

⚠ ⭐ **A bucketed / windowed search is wrong if written naively.** A first attempt this session reported
0.650 crossings per fragment against a correct 0.064, because it breaks when one node spans many
buckets — and human intergenic nodes reach megabases (chr1 mean node 2,552 bp against a median of 160).
A plain forward scan measured the same speed and is correct. **Whatever ships must be equality-tested
against a full binary search over the real chr1 cut array, not a toy.**

### 3.3 The junction lookup — a CSR keyed by cut index

⭐ Every annotated intron has **both** endpoints as partition cuts (measured: 100.00 % of 404,168). So
the "is this intron annotated?" test is the binary search the deposit is already doing:

```cpp
std::vector<std::int32_t> sj_offsets_;   // size n_cuts + 1, CSR over donor cut index
std::vector<std::int32_t> sj_acceptor_cut_;  // acceptor CUT INDEX, not a coordinate
std::vector<std::int8_t>  sj_strand_;    // disambiguates a strand-coincident pair -- nothing else
```

`sj_edge_id(donor_cut, acceptor_cut, sj_strand)` scans `sj_offsets_[d] .. [d+1]` — one to three
entries at human scale (⭐ 70.4 % of cuts are not a donor at all; over those that are, mean fan-out
**1.31**, max 25). If the intron's start is not a cut it is unannotated and the search never happens. No
hash map in the hot loop; the arrays are built once per reference from `edges.feather`'s `kind == 1` rows.

### ⛔ The junction-edge id IS the CSR slot. There is no third array.

The scan returns `k`, its position in `sj_acceptor_cut_`, and `sj_count[k]` is indexed by it
directly. **The draft's `junction_edge_id_` does not exist** — it was an indirection to
`edges_df`'s row, and taking that row as the bank index is a memory-safety bug: ⭐ the highest junction row
is **1,447,755** in a **404,168**-row bank, so `junction[edge_row]` runs off the end, and even in bounds it
permutes every surviving row.

⭐ **Measured, on both the human and the `ambig_dense_10mb` index:** the CSR slot equals the dense rank
within `flatnonzero(kind == JUNCTION)` for **all** junctions (404,168/404,168 and 537/537), because
`edges_df`'s junction rows are already in `(ref, src, dst, strand)` order and the CSR's
`lexsort((strand, acceptor_cut, donor_cut))` is the same order under the per-reference shift. So
`JunctionEdgeArrays.edge_row` is a strictly ascending map **CSR slot → `edges_df` row**, needed *only* to
join a payload row back to the index. **It does not cross the ABI into the accumulator.**

⚠ Two things follow, and both are load-bearing:

* **The builder's sort order is part of the contract.** `build_junction_edge_arrays` and the reference's
  `Partition.from_cuts` must sort identically, because the id *is* the rank. They disagreed once during S2
  (`(acceptor, donor)` against `(strand, acceptor, donor)`), which permutes every row and breaks
  byte-identity — reachable only via a strand-coincident junction pair, which is exactly what one spec
  test constructs.
* ⚠ `JunctionEdgeArrays.edge_row`'s docstring in `splice_graph.py` still calls itself *"the junction-edge
  id"*. It is not; it is the join key. **Fix the docstring in S3** — it is the one place a C++ author would
  look for the id's definition.

### 3.4 The fragment-length pools — five, each structurally pure

The current axis (3 region types × 2 compartments = 6 pools, `accumulator.h:57-63`) is replaced, not
renamed. The new pools are chosen so each is **pure by construction**, which is what removes the
circularity: nothing is estimated from the fragments it will later explain.

```cpp
enum class FragmentPool : std::uint8_t {
    kDnaIntergenic       = 0,   // contained in an intergenic node
    kDnaIntronic         = 1,   // contained in an intronic node
    kDnaIntronExon       = 2,   // crossing a contiguous edge with {intron, exon} flanks
    kDnaIntergenicExon   = 3,   // crossing a contiguous edge with {intergenic, exon} flanks
    kRnaSpliced          = 4,   // on an ANNOTATED junction edge
};
inline constexpr std::size_t kNFragmentPools = 5;
```

The two crossing pools are the "splash" reads — fragments that straddle a capture-probe edge — and they
are the only gDNA population that is *on target*, so a length model fitted from the intergenic pool
alone is mis-centred for exactly the fragments that leak. Having them as named pools makes that
comparison a QC output rather than a guess.

**No pool for an exonic contained fragment or an exon↔exon crossing**: those are gDNA/RNA mixtures, and
an impure pool is worse than a missing one.

⭐ **Bin at `L`, not at the genomic span.** The current code deliberately bins at span
(`accumulator.cpp:144-150`) because its `L` was the *covered* length, which drops the mate gap and
collapses the gDNA distribution to a spike at twice the read length. Our `L` already includes the mate
gap and excludes introns, so it *is* the molecule length for both components — one rule, no comment
needed to defend it.

The flank types of a contiguous edge come from `region_types` (the `coarse_type_array` the scanner
already receives), so the pool index is `f(type[j-1], type[j])` — no new input crosses the ABI.

⚠ A fragment with an **inferred** intron (`SPLICE_IMPLICIT`) must not enter `kRnaSpliced`: its splice is
a model inference, not an observation. That is what `FragmentPath::introns_inferred` is for.

### 3.5 The payload — `scan_payload.py`

Same file, same class name, no version suffix. `@dataclass(frozen=True, slots=True)`, fields typed bare
`np.ndarray` with an inline `# uint32[n, 2] — meaning` comment, all coercion and validation inside
`from_scan_result`, derived quantities as `@property`.

```
node_contained_count      uint32[n_nodes, 2]
node_contained_density    uint64[n_nodes, 2]
node_spanning_count       uint32[n_nodes, 2]
node_spanning_density     uint64[n_nodes, 2]
node_start_count          uint32[n_nodes]
edge_unspliced_count      uint32[n_edges, 2]
edge_unspliced_density    uint64[n_edges, 2]
edge_spliced_count        uint32[n_edges, 2]
edge_spliced_density      uint64[n_edges, 2]
sj_count            uint32[n_junctions, 2]
sj_density          uint64[n_junctions, 2]
ref_node_offsets          int64[n_refs + 1]
ref_edge_offsets          int64[n_refs + 1]
cut_positions             int64[n_cuts]
ref_cut_offsets           int64[n_refs + 1]
fragment_length_pools     int64[5, max_fragment_length + 1]        # §3.4 -- the pool axis changes
qc                        design §10.3's denominators, as fields
```

⚠ **Document the ownership**, because it is a live footgun: `np.ascontiguousarray(x, dtype=D)` is a
**no-op** when C++ already emitted `D`, so the payload holds **views over capsule-owned C++ heap** and
the payload object is the keep-alive. Someone "adding a cast for safety" would silently double peak
memory. Say so in the docstring.

⚠ **Delete, do not port:** `scan_payload.py:143-149`'s *"zeros for a stale payload that predates the
field"* — a backwards-compat path the rules forbid. And `src/rigel/_accumulator.py` **whole** (192
lines): a row-view façade that exists only so the old spec tests can write
`acc.regions[i].contained[ch]`.

⭐ If the new storage is **SoA flat vectors** rather than AoS structs, the per-ref → flat copy at
`bam_scanner.cpp:1947-1975` collapses from an element-wise per-channel loop into one `memcpy` per array
per reference, and the nanobind row views stop needing strides. Worth doing; it is the same edit.

### 3.6 The Python reference — `tests/native/_accumulator_reference.py`

Rewritten in place; still the byte-for-byte specification the C++ is gated against. ~300 lines. Model
it on `calibration/node_chain.py`, which this codebase treats as its best-written module: a module
docstring that draws the structure and states the addressing scheme before any code, one frozen+slots
dataclass with an inline `# int64[n] — meaning` per field, one loop over references with everything
inside it vectorised, and error messages that state the invariant, derive the expected value, and
pre-empt the wrong fix.

## 4. What dies in calibration

**`SubstrateView` survives** — it is a fine per-object statistics record.

**`CalibrationSubstrate` and `BoundarySubstrate` do not.** Their entire content beyond `.contained` is
`.left`/`.right` (`substrate.py:111-120`) and the same numbers re-keyed by boundary
(`substrate.py:170-173`). Two classes holding one set of numbers in two keyings exist *solely* because
a boundary has two sides. Target: **one substrate type with two arrays of objects (nodes, edges)** — no
left/right, no re-keying identity, no `_make_view`. **Keep the alignment check** (`substrate.py:130-148`)
— that one is load-bearing.

**Ranked by risk, highest first:**

| | what | why |
|---|---|---|
| **R1** | `build_node_geometry` + `NodeGeometry` (`node_geometry.py:61-315`) | 213 lines, **18 per-face arrays**, plus the junction-strand routing (`_spliced_faces`, :180-187), the exon-bit gating (:165-178) and the `_continues`/`_eff_spl_face` far-boundary machinery. Every "face" concept dissolves |
| R2 | `bp_solver`'s per-face consumers | the `(left, right)` tuples threaded through six functions, and the deliberate scalar/vector twin pairs that **must not be merged** (a measured 15.7×/op) |
| R3 | `fl.py`'s pool axis | the pools are re-keyed to the five of §3.4 |
| R4 | `effective_length.py` | the three old divisors go; the one `placements` formula replaces them |

## 5. The sequence — a clean cut

**Owner ruling: clean cut.** No additive-then-subtractive migration, no temporary duplicate fields. The
old objects are replaced in one pass and the pipeline is red until S5 closes.

That decision buys concision and costs attribution, so the plan must replace what attribution would have
given us. It does, in one move: **the Python reference (S2) is a working, correct implementation of the
whole deposit rule before any C++ or any consumer changes.** During the red period it is the oracle for
both halves — "is the C++ right?" is answered by byte-identity against it, and "is calibration rewired
right?" is answered by running calibration on a reference-produced payload. Neither question waits on
the other, and neither is answered by staring at a bench delta.

Which makes **S2 the most important step in this plan, not a formality.** It is written first, its tests
are written before it and verified failing, and nothing downstream starts until it is green.

| step | what | gate |
|---|---|---|
| **S1** | **Index.** `splice_graph._contiguous_reaches` → genomic span reach (§P1). Add `build_junction_edge_arrays(index)` returning the §3.3 CSR | bench **bit-identical 32/32** — nothing reads reach; 5 named tests in `test_splice_graph.py` updated; I1–I13 pass at human scale. **Independent of everything else; land it alone.** |
| **S2** | **The Python reference**, rewritten in place, plus design §10.4's test matrix written first and verified failing. Plus the scan-profiling harness (§6) | every test green; the reference runs end to end on a **real cfRNA payload** through a shim, and its QC denominators are sane |
| **S3** | **C++**: types, `deposit`, `merge_from`, the five pools, `set_junctions(...)` as a **second** method (the one-shot guard at `bam_scanner.cpp:1118-1122` throws today); dedupe introns on `(ref,start,end)` with strands OR-ed (§8 trap 1) | **byte-identical to S2's reference** on real data; **bit-identical at 1/2/4/8 workers** — newly achievable now every channel is integer; **no regression against ~357 ns/fragment end-to-end** (§6 — *not* the 124 ns microbenchmark) |
| **S4** | **Payload** rewritten; `_accumulator.py` deleted; SoA + `memcpy` at build-result | schema tests; `Σ node_start_count == n_accepted` on real data |
| **S5** | **Consumers, in one pass**: substrate collapsed to one type (`CalibrationSubstrate`/`BoundarySubstrate` deleted), `build_node_geometry` rewritten (R1), `effective_length` reduced to the one `placements` formula, `fl.py` re-keyed to the five pools | calibration runs end to end; the delta is **measured and recorded**, not required to be zero — there is no comparable baseline across this change |
| **S6** | **Delete.** `ruff check src/ tests/ scripts/`; undefined-name failures are the authoritative list | suite green; goldens regenerated **once** |

⚠ Goldens already moved at the index change (`gdna_em_count` fell 16–52 %) and that movement has never
been validated against truth. The new benchmark suite — not this sequence — is what should adjudicate it.

⚠ **The one thing a clean cut forbids:** starting S5 before S3 is byte-identical to S2. If both halves
are in flight at once the oracle is gone and the red period becomes unbounded.

### P1 — the one index change, in detail

`splice_graph._contiguous_reaches` (`splice_graph.py:501-534`) computes **exonic** reach over **mature**
transcripts: `lo = ma.before[...] + (p − ma.start[...])`, maximised per strand. It is nonzero only where
an exon strictly contains the position, so it is **zero throughout every intron**.

The design (§7.1) needs the opposite. Nascent RNA is a single-exon transcript spanning the gene, so RNA
reach is the **genomic distance to the ends of the widest transcript covering the position**, per
strand — and it is *not* zero in an intron, because that is where nascent RNA lives. A mature-only
reach declares zero RNA opportunity across every intron in the genome.

Replace the cumulative-exonic-offset arithmetic with a per-strand max over transcript **spans**. Same
four columns, same dtypes, same file. Add an invariant: reach is zero outside every transcript span and
monotone within one.

## 6. Performance

⭐ **Measured this session** against the shipped `accumulator.cpp`, driven by the **real chr1 v8 cut
array** (97,537 nodes, median 160 bp, mean 2,552 bp):

| | |
|---|---|
| current deposit, C++ microbenchmark | **124 ns / fragment** |
| of which the per-fragment `std::vector<Slice>` | **22.8 ns (18 %)** |

⭐ **S2 built the harness** — `scripts/design/scan_profile.py`. There had been none in the repository;
`scratchpad/perf_*.py` profile *calibration*, and every recorded scan number came from somewhere that no
longer exists. It isolates the accumulator by scanning each BAM twice, with and without `set_regions`, so
everything else cancels.

⚠ **The S3 gate is `~357 ns/fragment` end-to-end, NOT the 124 ns above.** Wiring the accumulator adds two
costs that scale differently and must not be averaged: `O(fragments)` is the deposit, `O(partition)` is
allocating and zeroing a per-worker accumulator, merging workers, and copying into the flat payload.
Against the 1.04 M-node human partition the second dominates a shallow BAM, so a single-BAM
`delta / n_fragments` is not a per-fragment cost. Regressing `delta(n) = fixed + c·n` over three real
cfRNA libraries gives **c ≈ 357 ns/fragment** (treat as ~350–400; three points, visible scatter) and
**fixed ≈ 0.108 s**. The end-to-end figure is ~2.9× the microbenchmark because the A/B also contains
`fragment_genomic_spans` and the scan-side gate in `deposit_to_accumulator` — and it is the one a user
pays.

⭐ Projected single-thread deposit cost, partition excluded: **3.6 s at 10⁷ fragments, 35.7 s at 10⁸.**
That is the regime this tool is for, and it is why the 22.8 ns per-fragment allocation is worth removing.

**Non-negotiables for the new deposit:** no allocation in the per-fragment path (`WorkerState` already
keeps reusable span scratch vectors at `bam_scanner.cpp:391-397` — follow that pattern); the crossed-line
enumeration stays a contiguous index range, never a container; junction lookup stays a short CSR scan.

## 7. How the code must read

Observed conventions, not preferences — this is what the tree does everywhere.

* Module docstring is `"""rigel.<pkg>.<mod> — one-line purpose."""` followed by a body that **is** the
  schema or the derivation (`scan_payload.py:1-40`, `splice_graph.py:1-66`).
* `from __future__ import annotations`; stdlib / third-party / relative import groups; `__all__` either
  straight after the imports or as the last statement — one or the other per module.
* `@dataclass(frozen=True, slots=True)`; arrays are plural snake_case nouns carrying their axis
  (`ref_cut_offsets`); offset arrays are `*_offsets` of length `n_refs + 1`, CSR-sliced.
* Every field gets an inline `# uint32[n, 2] — meaning` comment. Derived values are `@property`, never
  stored.
* Errors are `ValueError` naming **both** the observed and the expected value.
* No Greek in identifiers. No magic numbers — every divisor derived from the deposit rule and unit-tested
  against brute-force enumeration.
* C++: `static_assert` every struct size; `noexcept` on the hot path; `nb::class_` with `def_prop_ro`
  zero-copy views.

## 8. Traps found while reading — each needs a test

1. ⛔ **Introns can appear twice with different strands.** `parse_bam_record` reads `XS`/`ts` once per
   *record* (`bam_scanner.cpp:618`) and stamps it on all that record's junctions; `build_fragment`
   inserts into a `set<tuple<ref,start,end,strand>>`. A pair where R1 has `XS` and R2 does not yields
   the same intron twice. `fragment_genomic_spans` survives this by accident (`pos = max(pos, c1)`); a
   junction deposit iterating `f.introns` would **double-credit**. **De-duplicate on `(ref, start, end)`
   with the strands OR-ed before crediting.**
2. ⛔ **`build_fragment` keys blocks by `(ref, strand)`** (`:653-656`). Mates whose orientations disagree
   produce two groups and `align_strand = AMBIGUOUS`, which the gate at `:1474-1480` silently drops.
   Count it in QC rather than losing it.
3. **Comment says zero-copy; the code copies** (`bam_scanner.cpp:1902-1904` vs `:1929-1975`). The copy is
   correct — delete the comment.
4. **Stale stride comment**: `bam_scanner.cpp:2666` says a `Boundary` is 48 B / 12 floats; it is 64 B /
   16. The code is right because it computes `sizeof(Boundary)/sizeof(float)`; only the comment lies.
   It repeats at :2682, :2697, :2712.
5. **FL pooling disables itself silently** when `region_types` is empty or `max_fl <= 0`
   (`accumulator.cpp:30-72`). Make it explicit or make it throw.
6. **`test_float_mass_channels_are_NOT_bit_identical_across_shardings`** (`test_accumulator_spec.py:652`)
   asserts the *current* nondeterminism. Once every channel is integer this becomes false by
   construction — **delete it, do not port it**, and replace it with the positive assertion.

## 9. Settled since draft

| | ruling |
|---|---|
| **naming** | **`spanning`**, everywhere, including the design doc |
| **sequencing** | **clean cut** (§5). S2 becomes the oracle that replaces per-step attribution |
| **fragment-length pools** | **five, structurally pure** (§3.4): DNA intergenic / intronic / intron-exon crossing / intergenic-exon crossing, and RNA annotated-spliced. Binned at `L` |
| **`start_count`** | its **own flat array**, not a `Node` field. `Node` is 48 B with no padding |
| **strand** | ⭐ **ONE convention: genome strand everywhere** (design §5.1). Sense/antisense derived, never stored |
| **minimum intron length** | **none — every CIGAR `N` is an intron** (design §12.A) |
| **malformed introns** | overlapping **and abutting** introns merge, and the merge is counted (design §12.A) |
| **junction-edge id** | the **CSR slot**. `edges_df.edge_row` is a join key and does not cross the ABI (§3.3) |

Nothing is open. S1 and S2 are done; the next step is **S3**.
