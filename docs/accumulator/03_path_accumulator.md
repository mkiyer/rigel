# The path accumulator — design and implementation plan

**v1, 2026-07-28.** Consumes the splice graph (`../index/00_splice_graph_design.md`), produces the path
counts that `../calibration/PATH_MARGINALIZATION.md` turns into node/edge evidence. Derivation of *why*:
`02_redesign_derivation.md`.

---

## 1. CONTRACT

> Convert each deposited fragment into a **path** through the splice graph and count paths. One fragment
> increments exactly one path. Nothing is split, nothing is fractional, nothing is double-counted.

Inputs: the graph's per-reference cut arrays + junction edge table; the fragment's contiguous genomic
spans (already produced by `fragment_genomic_spans`, shipped); its channel.
Output: a path table with per-channel counters (§3).

---

## 2. PATH REPRESENTATION

A path is the ordered node sequence a fragment occupies. Between junctions the walk is contiguous, so the
node sequence is fully determined by the endpoints and the junctions used. The stored key is therefore:

```
    key = ( first_node : int32,  last_node : int32,  junction_edge_ids : int32[n_j] )
```

* `n_j == 0` → **unspliced**: the path is `first .. last` contiguous. Key packs into a single `uint64`
  (`first << 32 | last`) — the fast path, no allocation, no hashing beyond an integer mix.
* `n_j > 0` → **spliced**: hash the key bytes; compare exactly on collision.

`spliced` is implied by `n_j > 0` and is not stored separately. Strand is carried by the **channel**, not
by the key — unchanged encoding, `ch = (spliced ? 2 : 0) + (primary ? 0 : 1)`, so a path row holds four
channels exactly as `region_contained` does today.

⚠ **Do not key paths on "skipped node ids"** as a way to infer splicing. It is tempting (a junction skips
ids) but wrong at a **bookended** junction — two exons abutting with a zero-length intron do not skip an
id — and wrong for a fragment that legitimately crosses many contiguous nodes. The junction id list is the
only correct discriminator, and it is needed anyway for the effective length.

---

## 3. WHAT IS ACCUMULATED — three counters per path per channel

```
    N_p   uint32    fragment COUNT on this path                     (the integer atom)
    D_p   float64   Σ_f 1 / L_f      over fragments on this path     (the density sum)
    F_p   float64   Σ_f o_first,f / L_f                             (the first-node overlap sum)
```

`L_f` is the fragment's **path length** — the number of bases it occupies along the path (exonic only;
introns are not part of `L`). `o_first,f` is its overlap with the path's first node.

Why exactly these three, and why `F_p` is not optional: `PATH_MARGINALIZATION.md` §3 derives it, and the
short version is that a path's **interior** nodes are always *fully* covered (the graph cuts at every exon
endpoint, so a junction donor is always a node boundary), so an interior node's share is `|n| · D_p` —
`D_p` alone. Only the two **end** nodes have variable overlap, and they are constrained by
`o_first + o_last = L_f − Σ_interior |n|`, so one extra accumulator determines both:

```
    interior node n :   |n| · D_p
    first node      :   F_p
    last node       :   N_p − (Σ_interior |n|) · D_p − F_p
```

**This is exact per fragment, not an approximation.** Three numbers replace an unbounded per-fragment
record.

⚠ **M-acc2 must measure whether `F_p` earns its 8 bytes.** The alternative is to substitute the
*expected* first-node overlap given the path geometry and the FL pmf. That is an approximation; it may be
a good one. Measure before dropping it.

⚠ **Variance.** The node marginal is a weighted sum, so its variance is `Σ_f (o/L)²`, not its mean. Inside
a path the weights vary only over `L_f ∈ [Σ_interior + 2, Σ_interior + |first| + |last|]`, which bounds the
dispersion, but the honest quantity is a Kish `n_eff`. Adding `Q_p = Σ 1/L_f²` makes it exact for interior
nodes at 8 more bytes. **Deferred pending M-acc3** — do not add it on aesthetics.

---

## 4. DEPOSIT ALGORITHM

```
deposit(spans, channel):
    # 1. walk spans -> node path
    first = node_of(spans[0].start)
    last  = node_of(spans[-1].end - 1)
    juncs = []
    for gap between span[i].end and span[i+1].start:
        e = junction_edge(donor = span[i].end, acceptor = span[i+1].start, strand)
        if e is None:  -> §7 unannotated
        juncs.append(e)
    # 2. per-fragment quantities
    L       = Σ span lengths
    o_first = min(end(first), spans[0].end) - spans[0].start
    # 3. one increment
    pid = intern(first, last, juncs)
    row = table[pid][channel]
    row.N += 1;  row.D += 1.0 / L;  row.F += o_first / L
```

`node_of` is a `searchsorted` into the reference's cut array — the same primitive
`Accumulator::region_of_pos` already uses. `junction_edge` is a lookup into a
`(donor, acceptor, strand) → edge_id` hash built once at load.

**Per-fragment cost**: `O(n_spans · log n_nodes)` for the lookups + one hash insert. The current deposit
is `O(nodes touched)` with two array increments per touched node; the new one is `O(1)` increments and one
hash probe. **It should be faster on multi-node fragments and marginally slower on contained ones.**
M-acc1 measures it; do not assume.

---

## 5. NATIVE DESIGN

* **Per-worker open-addressing hash map** `key → path_id`, plus a growable counter array. Mirrors today's
  per-worker `Accumulator` + `merge_from`, so the threading model is unchanged.
* **Merge**: for each worker map, probe the master and add the four channels' `(N, D, F)`.
* ⭐ **Accumulate `D` and `F` in FIXED POINT, not floating point.** A float accumulation makes the result
  merge-order dependent — the known `calibrate_cross_process_nondeterminism` (~2.6 %, from C++ parallel FP
  reduction), which today's float32 mass already suffers. Integer accumulation is **associative**, so the
  result is **bit-exact at any thread count**:

  ```
      D_fixed += round(2^40 / L)          uint64
      F_fixed += round(2^40 * o_first / L)
  ```

  Range check: `L ∈ [~50, ~2000]` ⇒ each term ≤ `2^40/50 ≈ 2.2e10`; at `10^8` fragments the sum is
  `≈ 2.2e18` against a `uint64` ceiling of `1.8e19`. Fits with an order of magnitude to spare, and the
  relative precision of a single term is `2^-40 · L ≈ 2e-9`. **This removes a long-standing nondeterminism
  source rather than extending it** — a strict improvement over the status quo, not merely a wash.
* **Fast path**: `n_j == 0` keys are a `uint64`; keep them in a separate map from spliced keys so the
  common case never touches variable-length comparison.
* **No `float32`.** `D` accumulates ~10⁷ terms of order `1/200`; `float32` has 24 bits of mantissa and
  would lose ~4 decimal digits. `float64`, and Kahan if M-acc4 shows drift.

**Memory** (human, gated by M-acc1): key 8–24 B + 4 channels × 20 B = **~100 B per distinct path**.

⚠ **Do not read that as "small".** The *possible* unspliced path count is roughly
`n_nodes × (reachable last nodes)` ≈ 992 K × 3–6 ≈ **3–6 M**, and a deeply sequenced library saturates a
large share of it. At ~100 B that is **0.3–0.6 GB, plus spliced paths.** This is a first-order design
pressure, not a footnote — M-acc1 is the go/no-go and it should run before Y2, not after.

⭐ **And the schema size depends on a decision that is currently deferred.** `D_p` and `F_p` exist only for
the **coverage** marginal; the **anchored-cell** marginal needs `N_p` alone (16 B + key, a **4× reduction**).
`PATH_MARGINALIZATION.md` §4 currently says "ship Option E, measure Option A later" — but the memory
consequence means **the attribution decision should be taken before the accumulator schema is frozen.**

Fallbacks if M-acc1 fails, in order of preference:
1. store only the channels that occur (most paths are single-channel) — ~4× on the counters;
2. drop `F_p` (M-acc2) — 32 B;
3. cap the contiguous span `last − first` and fold the tail into a coarser key, **logging what was
   folded** (`no silent caps`).

---

## 6. AMBIGUOUS FRAGMENTS — the small side buffer

Owner decisions D4/D5. `FragmentBuffer` is **not** touched.

* A fragment whose path is **unique** deposits directly during the scan. This is the overwhelming
  majority.
* A fragment with several compatible paths — an **implicit splice** whose mate gap matches more than one
  annotated intron, or a **multimapper** — is appended to a small side buffer holding only
  `(candidate path keys, channel, L, o_first)`. Nothing else.
* **Pass 2** runs over that buffer alone, after the RNA FL model is trained, and assigns each fragment to
  its **maximum-likelihood** path under the fitted FL pmf and the paths' own counts.

⚠ **Assignment must stay integral.** A fractional split would reintroduce exactly the non-integer
observable this redesign removes, and `N_p` would stop being a count. An EM with fractional weights is a
later, separately measured arm; if it is ever adopted the variance model must change with it.

**M-acc5** measures the ambiguous fraction on real cfRNA before the buffer is sized. Multimapper recovery
(D5) is out of scope for this phase and is enabled by, not part of, the buffer.

---

## 7. ⚠ UNANNOTATED JUNCTIONS — the open hole

A fragment may cross a junction that is not in the annotation. Today it deposits spliced mass on whatever
boundary it touches. **On a graph it has no edge to walk**, and this design must say what happens.

| option | effect |
|---|---|
| (a) hold out, like `SPLICE_ARTIFACT` | honest, loses real RNA evidence, biases against novel splicing |
| (b) add the edge at runtime | the graph stops being an index artifact; every consumer must handle a growing edge set |
| (c) ⭐ deposit on a path keyed `(first, last, UNANNOTATED)` | keeps the RNA evidence at both ends, needs no graph mutation, and the effective length falls back to the two end nodes' geometry |

**Recommend (c), gated on M-acc6**: measure the unannotated-junction fragment fraction on real cfRNA
first. If it is negligible, (a) is simpler and more honest; if it is not, (c) is required and its
effective-length treatment must be specified before it ships. **Do not decide this from the synthetic
suite** — it has no novel junctions by construction.

---

## 8. TESTS

Reference implementation `tests/native/_path_accumulator_reference.py`, byte-for-byte against the C++, in
the same shape as today's `_accumulator_reference.py`.

| # | test |
|---|---|
| A1 | contained fragment → path `(n, n)`, `N=1`, `D=1/L`, `F=1` |
| A2 | two-node contiguous crossing → `N=1`, `D=1/L`, `F=o_first/L` |
| A3 | fragment spanning 4 contiguous nodes → **one** path, interior nodes fully covered |
| A4 | simple spliced fragment → path with 1 junction id |
| A5 | multi-junction fragment (the `(B→D→F→I)` tiny-exon case) → distinct path from `(B→I)` |
| A6 | bookended junction (zero-length intron) → junction id present although no node is skipped |
| A7 | **count conservation**: `Σ_p N_p == n_deposited`, exactly, integers |
| A8 | ⭐ **density conservation**: for every fragment, the marginalised node shares sum to `1` (`Σ_n o_n/L = 1`) — the §3 identity, checked per fragment |
| A9 | merge determinism: N workers in any order → identical `N`; `D`/`F` identical at `total_threads=1` |
| A10 | opposite-strand junctions at the same coordinates → distinct paths |
| A11 | fragment running off the reference end → clipped as today, still one path |
| A12 | unannotated junction → whatever §7 decides, asserted explicitly |
| A13 | ambiguous fragment → side buffer, **not** the main table; pass 2 assigns exactly one path |

---

## 9. MEASUREMENTS THAT GATE THIS DESIGN

```
M-acc1  distinct-path count + resident memory, cfrna:LBX0190          -> §5, the go/no-go number
M-acc2  is F_p needed, or does the geometric expectation suffice?     -> §3
M-acc3  weight dispersion within a path: does the variance need Q_p?  -> §3
M-acc4  float64 drift in D over ~10^8 fragments (Kahan or not)        -> §5
M-acc5  ambiguous-fragment fraction on real cfRNA                     -> §6
M-acc6  unannotated-junction fragment fraction on real cfRNA          -> §7
M-acc7  deposit throughput vs today, at genome scale                  -> §4
```

⚠ Every one of these on **real cfRNA**. The toy has 1,698 regions, no novel junctions, and is Poisson by
construction.

---

## 10. PHASES

```
Y0  M-acc1..7, then spec lock
Y1  reference implementation (Python) + tests A1-A13 written to target, failing
Y2  C++ path accumulator + nanobind; per-worker map + deterministic merge
    Gate: byte-for-byte vs the reference; A1-A13 green; M-acc7 budget met
Y3  dual-write today's AccumulatorPayload from the same spans
    Gate: legacy payload BIT-IDENTICAL 32/32 (it adds a table nobody reads yet)
Y4  side buffer + pass 2 (§6)
Y5  hand off to PATH_MARGINALIZATION.md
```
