> ## ⛔ AMENDED 2026-07-29 — the graph SHIPPED, and three gates below did not survive contact
>
> * **P2 is unsatisfiable** (the merge deletes cuts) and **P3 is retired with its comparator.** The v7
>   partition was deleted at accumulator plan W1b-clean, so "merging v8 reproduces `regions.feather`"
>   has nothing to compare against. Its content — an independent check on the signature — is now
>   **I3b**: `validate_graph` recomputes each node's signature from its midpoint by direct interval
>   containment and asserts equality. No v7 needed, runs at build and load, proven to fire.
> * **§9's option (b)** ("hold v8 behind the loader until the accumulator lands") was **not** taken.
>   Option (a) shipped: the graph IS the partition as of W1b, measured. ⚠ Its warning that "shipping
>   v8 alone changes calibration's numbers" is right in general and **false on the 32-condition
>   bench**, whose annotation has zero mergeable adjacencies. See ledger W1b.
> * **§8's scan-deposit budget is MEASURED**, on real cfRNA at the +38.7 % human partition:
>   **+9.7 %** (LBX0190) / **+6.9 %** (MO_3021). Measurable, small, and now on record for D5/W6.
>
> **Still current (2026-07-28).** The graph model, builder, invariants and test matrix below are otherwise
> unchanged by the accumulator's v5 revision. The consumer is
> [`../accumulator/05_accumulator_v5.md`](../accumulator/05_accumulator_v5.md).
> ⭐ **One addition v5 requires:** junction edges also store `reach_donor` / `reach_acceptor` (`int32`
> each) — the transcript's cumulative exonic length either side of the junction, over the annotation. They
> feed the junction divisor `E_J` (v5 §4.2b) and are pure index-time arithmetic — the exonic length
> remaining either side of the junction, maximal over the annotation's isoforms (v5 §4.4).

# The splice graph index — design and implementation plan

**v1, 2026-07-28.** The foundation of the accumulator rework. Scoped as its own effort: the index is
built once and everything downstream inherits its bugs, so this document specifies the model, the on-disk
schema, the build algorithm, the invariants, and a test matrix that must be green before anything consumes
it.

Companion: `../accumulator/05_accumulator_v5.md` — what runs on this graph, and why.
(`02_redesign_derivation.md`, `03_path_accumulator.md`, `04_accumulator_v3.md` and
`../calibration/PATH_MARGINALIZATION.md` are superseded drafts; do not cite them.)

---

## 1. SCOPE

Replace the calibration region/boundary partition with a **splice graph** persisted in the index.

| | today (v7) | v8 |
|---|---|---|
| partition rule | maximal intervals of **constant signature** (adjacent equal-signature segments **merged**) | every distinct annotation event is a cut; **no merging** |
| regions (human) | 752,654 | **992,068** (+32 %) |
| median region | 305 bp | **145 bp** |
| transcript termini visible | 40.5 % | **100 %** |
| boundaries | `k+1` positional rows per ref, `[boundary_id, ref_name, position]` | **edges**: `k−1` contiguous + 404,168 junction (human) |
| splice junctions in the map | absent — inferred at runtime from the observed motif | **first-class directed edges** |
| TSS/TES | absent | **typed, per strand, per side** |

**Why no merging is the whole point.** A terminus is always an exon endpoint, so it is always in the event
set — but if the signature does not change there (an alternative TSS inside another transcript's exon), the
merge deletes the cut. Measured on the real human index: **59.5 % of transcript termini fall strictly inside
a merged region.** Dropping the merge preserves all of them with no extra machinery. The cut set does not
change; only the merge step is removed.

**Non-goals.** No change to `transcripts.feather`, `intervals.feather`, `sj.feather`, the resolver, the
scorer, the EM, or the locus builder. This document covers `nodes.feather` + `edges.feather` and their
builder, loader and validators only.

---

## 2. THE GRAPH MODEL

### 2.1 Nodes

A **node** is a half-open genomic interval `[start, end)` on one reference. Nodes tile each reference
exactly and are numbered globally in genomic order (reference order = `ref_lengths` iteration order, then
`start`), so **node ids are contiguous and increasing within a reference** — the property every downstream
CSR and the topological order depend on.

```
node_id  int64    global, genomic order
ref_id   int32    index into ref_names
start    int64    half-open
end      int64
length   int64    end - start   (derived, stored for convenience/validation)
signature uint8   4-bit {intron_neg, intron_pos, exon_neg, exon_pos} — unchanged encoding
```

⚠ **The `neighbour-differs` invariant is DELETED.** Adjacent nodes may share a signature; that is the
point. It is replaced by a stronger, more meaningful invariant (§6, I4): *every interior interface is an
annotation event.*

### 2.2 Edges

An **edge** is a directed transition between two nodes, always `src < dst` in genomic order — so **genomic
order is a topological order of the whole graph**, and no BFS or explicit toposort is needed anywhere.

```
edge_id     int64   global, sorted by (src, kind, dst)
src         int64   node id
dst         int64   node id, dst > src, same ref
kind        uint8   0 = CONTIGUOUS, 1 = JUNCTION
strand      int8    JUNCTION only: 1 = POS, 2 = NEG  (0 on CONTIGUOUS)
flags       uint16  CONTIGUOUS only: the structural bits at this position (§2.3)
```

* **CONTIGUOUS**: `dst == src + 1`, same reference. It *is* the genomic position `p = end(src) = start(dst)`.
  There are `k − 1` of them for a `k`-node reference. Undirected in physics, stored `src → dst` for
  ordering.
* **JUNCTION**: one per distinct `(ref, intron_start, intron_end, strand)`. `src` = the node **ending** at
  `intron_start` (the donor), `dst` = the node **starting** at `intron_end` (the acceptor). Directed by the
  splice motif's strand. Human: **404,168** distinct junctions, and **zero** of them are strand-coincident,
  so `strand` is a clean 2-value field (the old `genomic_sj_strand == 3` "both" case does not occur in
  GENCODE and must still be *validated*, not assumed).

⚠ The v7 terminal boundaries (`B0`, `Bk` at each reference edge) **disappear**. They are not transitions.
They never carried deposits (`00_design.md` §6 invariant 4) and the node chain can be rebuilt without them.

### 2.3 Contiguous-edge flags — the structural bits

8 meaningful bits, one byte, per strand `s ∈ {pos, neg}`:

| bit | meaning at position `p` |
|---|---|
| `TSS_s` | some strand-`s` transcript **starts** at `p` (its 5′ terminus) |
| `TES_s` | some strand-`s` transcript **ends** at `p` (its 3′ terminus) |
| `DONOR_s` | `p` is the `intron_start` of some strand-`s` junction |
| `ACCEPTOR_s` | `p` is the `intron_end` of some strand-`s` junction |

These are **not mutually exclusive** — that is exactly the case v1 of the calibration work showed the
signature is blind to (73 toy positions that are *both* terminus and junction; at human scale the majority).

⭐ **Set from `~is_synthetic` — ONE filter, the same as the event set** (ledger W1c). A brief attempt
at `~is_synthetic & ~is_nrna` deleted the TSS/TES of 26,475 real single-exon transcripts, because on
a non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent", not "manufactured span".
Every manufactured span is `is_synthetic`. Invariant **I13** now pins it.
Two derived predicates matter downstream and are computed, not stored:

```
    rna_crosses_contiguously(p, s)  =  exon_s(src) AND exon_s(dst)  AND NOT (TSS_s or TES_s at p)
    is_terminus(p, s)               =  TSS_s(p) OR TES_s(p)
```

`rna_crosses_contiguously` replaces the observational `accept_l/accept_r = (SP+SN) > 0` in `bp_solver`
(`P1G_SCOPE` C3) with structure.

### 2.4 What the graph is NOT

It is not a per-transcript path store. A transcript's path through the graph is recoverable
(`transcripts.feather` + `intervals.feather` are unchanged) but is not materialised — nothing downstream
needs it, and materialising 228 K paths of ~10 nodes would add 2.3 M rows for no consumer.

---

## 3. THE PARTITION RULE

### 3.1 The event set

For each reference:

```
    cuts = sorted(unique( {exon.start, exon.end : all exons of all NON-SYNTHETIC transcripts}
                          ∪ {0, ref_length} ))
    nodes = [ (cuts[i], cuts[i+1]) for i in 0 .. len(cuts)-2 ]
```

Every intron endpoint is already an exon endpoint of the same transcript, so introns contribute no new
cuts; every TSS/TES is a first/last exon endpoint, so termini contribute no new cuts either. **The cut set
is exactly today's cut set. Only the merge is removed.** That is what makes this a low-risk change to a
high-blast-radius artifact.

⚠ Synthetic / nRNA transcript rows are excluded, exactly as today (`build_region_partition` skips
`tx.is_synthetic`). Their spans are derived from real transcripts, so they add no events; including them
would change the partition and is a behaviour change nobody asked for.

### 3.2 Signature

Unchanged semantics: a node's 4 bits record whether any real transcript places an exon / an intron on the
`+` / `−` strand over that node. Because nodes are now sub-intervals of constant coverage by construction,
the signature is constant over a node — the same guarantee as today, obtained the same way.

### 3.3 Degenerate cases, and what each must do

| case | rule |
|---|---|
| two transcripts share an exon endpoint | one cut (the set dedups) |
| a 1 bp node (two events 1 bp apart) | **kept**. Human has **15,315** of them. They are legal nodes with `length = 1`; nothing downstream may assume `length > 1`. |
| an exon of length 0 in the GTF | rejected at parse time (already) — must raise, not silently drop |
| reference with no transcripts | one node spanning `[0, ref_length)`, zero edges |
| reference of length 0 | zero nodes, zero edges (already handled) |
| an exon endpoint at 0 or `ref_length` | no duplicate cut (the set dedups); no zero-length node |
| a junction whose donor/acceptor is not a node interface | **impossible by construction — assert it** (I5) |
| coincident opposite-strand junctions at the same `(start,end)` | two JUNCTION edges, one per strand. Zero occurrences in GENCODE; must still work. |
| a transcript on strand `.` / ambiguous | skipped for signature and junctions, exactly as today |

### 3.4 Scale (measured, human GENCODE, `refs/rigel_index_v7` inputs)

| | |
|---|---|
| real transcripts | 227,844 |
| real exon rows | 1,640,987 |
| **nodes** | **992,068** (median 145 bp; p25 59, p75 730) |
| nodes < 10 bp | 8.5 % (0.01 % of bp) |
| nodes < 200 bp | 58.1 % (1.45 % of bp) |
| CONTIGUOUS edges | ≈ 992,068 − n_refs |
| **JUNCTION edges** | **404,168** |
| total graph | ≈ 992 K nodes, ≈ 1.40 M edges |

---

## 4. ON-DISK SCHEMA (INDEX FORMAT v8)

```
nodes.feather    node_id:int64  ref_name:string  start:int64  end:int64  length:int64  signature:uint8
edges.feather    edge_id:int64  src:int64  dst:int64  kind:uint8  strand:int8  flags:uint16
```

* `regions.feather` / `boundaries.feather` are **removed**. `nodes.feather` is the strict superset of
  `regions.feather`; `edges.feather` replaces `boundaries.feather`.
* Row order is the id order and is part of the contract: `nodes` by `(ref order, start)`; `edges` by
  `(src, kind, dst)` — so the out-edges of a node are contiguous and a CSR is a single `searchsorted`.
* `manifest.json` gains `n_nodes`, `n_edges`, `n_junctions` for cheap validation and for the loader to
  pre-allocate.
* `INDEX_FORMAT_VERSION` 7 → **8**, with the history note appended in `index.py`. Loaders refuse v7
  (they already refuse mismatched versions), so **every index must be rebuilt** — see §9.

**Derived at load, not stored** (cheap, and storing them invites them going stale):

```
    ref_node_offsets   int64[n_refs+1]     CSR: nodes of reference f
    out_edge_offsets   int64[n_nodes+1]    CSR: edges leaving node i
    in_edge_offsets / in_edge_ids          the reverse CSR (one counting sort)
```

Load cost at human scale: two feather reads (992 K + 1.4 M rows) plus three `O(n)` passes — budget
**< 1 s**, gate in §8.

**Size on disk**: nodes ≈ 992 K × 33 B ≈ 33 MB raw, ~12 MB compressed; edges ≈ 1.4 M × 28 B ≈ 39 MB raw,
~10 MB compressed. Comparable to today's `regions` + `boundaries`.

---

## 5. BUILD ALGORITHM

Today's builder is a Python event sweep appending `dict` rows (`_build_reference_rows`). At 992 K nodes
that is the wrong shape. **The v8 builder is fully vectorised**; there is no per-node Python object.

Per reference `f`:

```python
# 1. cut set                                                  O(E log E)
cuts = np.unique(np.concatenate([exon_starts, exon_ends, [0, ref_len]]))
starts, ends = cuts[:-1], cuts[1:]                          # the nodes

# 2. signature: 4 difference arrays over the node index      O(E + n)
sig = np.zeros(n_nodes, np.uint8)
for bit, (iv_start, iv_end) in the four (exon±, intron±) interval sets:
    d = np.zeros(n_nodes + 1, np.int32)
    np.add.at(d, np.searchsorted(cuts, iv_start), 1)         # cuts are exact hits
    np.add.at(d, np.searchsorted(cuts, iv_end),  -1)
    sig |= np.where(np.cumsum(d)[:-1] > 0, bit, 0)

# 3. contiguous edges                                        O(n)
src = np.arange(n_nodes - 1); dst = src + 1

# 4. junction edges                                          O(J log n)
j_src = np.searchsorted(cuts, intron_start) - 1              # node ENDING at the donor
j_dst = np.searchsorted(cuts, intron_end)                    # node STARTING at the acceptor
assert (cuts[j_src + 1] == intron_start).all()               # I5
assert (cuts[j_dst]     == intron_end).all()

# 5. flags on contiguous edges                               O(T + J)
np.bitwise_or.at(flags, np.searchsorted(cuts, tss_pos) - 1, TSS_POS)   # etc.
```

**Complexity** `O(E log E)` per reference, dominated by one sort of ~3.3 M int64 events genome-wide.
**Budget: < 10 s and < 1 GB peak for human GENCODE** (gate in §8). Today's `dict`-per-row builder is the
thing being replaced, so this must be *measured*, not assumed — it is a required acceptance number, not a
nice-to-have.

**Determinism.** The output is a pure function of `(transcripts, ref_lengths)`. `np.unique` sorts, ids are
assigned by position, edges by `(src, kind, dst)` — no dict iteration order, no hashing, no parallel
reduction. **Byte-identical across runs and platforms** is an invariant, tested (§7, T-D1).

**Parallelism.** Per-reference independent; genome-wide the build is seconds, so ship it serial and
revisit only if §8's budget fails. Adding threads to a 10 s one-time build is not worth the determinism
risk.

---

## 6. VALIDATION INVARIANTS

Every one of these is a `validate_graph()` assertion run at build time and again at load time (the load
check is cheap and has already caught real breakage — the current loader's region/boundary cross-check is
what surfaces stale indexes).

| # | invariant |
|---|---|
| **I1** | nodes tile each reference exactly: `start[0] == 0`, `end[i] == start[i+1]`, `end[-1] == ref_length`, all `length > 0` |
| **I2** | `node_id` is `0..n-1` in genomic order; `ref_node_offsets` is non-decreasing and ends at `n_nodes` |
| **I3** | `signature` ∈ `[0, 16)`; a node's signature equals the recomputed coverage of its midpoint |
| **I4** | ⭐ **every interior interface is an annotation event** — `{end[i] : i interior}` ⊆ `{exon endpoints}` and ⊇ it. *(This replaces `neighbour-differs`.)* |
| **I5** | every JUNCTION edge lands on interfaces: `end(src) == intron_start`, `start(dst) == intron_end`, `src < dst`, same ref |
| **I6** | every JUNCTION `(ref, donor, acceptor, strand)` appears exactly once; the set equals the distinct real-transcript intron set |
| **I7** | CONTIGUOUS edges are exactly `{(i, i+1)}` within each reference — count `= n_nodes − n_refs_with_nodes` |
| **I8** | `src < dst` on **every** edge ⇒ genomic order is a topological order |
| **I9** | a JUNCTION's `src` carries the `exon_s` bit and its `dst` carries the `exon_s` bit for that junction's strand `s` |
| **I10** | `TSS_s`/`TES_s` flag counts equal the distinct terminus-position counts derived independently from `transcripts.feather` |
| **I11** | every transcript's exon path is walkable: consecutive exons are joined by a JUNCTION edge on that strand, and each exon's interior interfaces are CONTIGUOUS edges ⭐ *the end-to-end structural check* |
| **I12** | edge rows are sorted by `(src, kind, dst, strand)`; `edge_id` is the row index. ⚠ **strand is part of the key** — without it two strand-coincident junctions collide and read as a duplicate |
| **I3b** | ⭐ **the signature, recomputed independently**: each node's 4 bits, re-derived from its MIDPOINT by direct interval containment, equal the stored value. A different algorithm from the builder's cumulative-difference sweep over the same interval sets, so the two can only agree by both being right. *(Replaces the retired P3, which pinned the signature against the deleted v7 partition.)* |
| **I13** | ⭐ **the flags ARE the events**, both directions: each structural bit is set at exactly the positions that generate it — no event without its flag, no flag without its event. Subsumes I10. ⚠ **One exemption:** an interior interface carries NO flag when adjacent exons of one transcript are **bookended** (a zero-length intron) — an exon endpoint that is neither a terminus nor a splice site. Zero occurrences in GENCODE; constructed by G14. This is why I13 is *"the flags are the events"* and not the tempting *"every interface carries a flag"* |

**I11 is the one that catches real bugs.** It walks every real transcript through the built graph and
fails if any transcript is not representable. It is `O(total exons)` and it subsumes most of I5/I6/I9.

---

## 7. TEST PLAN

Existing hooks to build on: `tests/conftest.py::build_test_index(tmp_path_factory, gtf_text, genome_size)`
already builds a full index from a GTF string. The toy matrix below is a table-driven parametrisation of
that helper. New file `tests/calibration/test_splice_graph.py`; `tests/calibration/test_regions.py` and
`test_boundary_partition.py` are rewritten against the new artifact.

### 7.1 The toy GTF matrix — one graph assertion per case

| # | case | what it must prove |
|---|---|---|
| **G1** | single-exon transcript | 3 nodes (upstream / exon / downstream), 2 contiguous edges, 0 junctions, TSS+TES flags on the two edges |
| **G2** | 2-exon transcript | 1 junction edge donor→acceptor; the intron node exists and is not on the transcript's path |
| **G3** | 3-exon transcript | 2 junction edges; I11 walk succeeds |
| **G4** | ⭐ **alternative TSS inside another transcript's exon** | the interior cut **survives** (this is the case v7's merge deletes); both nodes carry the same signature; `TSS` flag on the interior edge |
| **G5** | ⭐ **alternative TES inside another transcript's exon** | mirror of G4 |
| **G6** | ⭐ **a position that is BOTH a terminus and a junction** (one transcript ends where another splices) | `TES_s` **and** `DONOR_s` both set on the same edge — the case the signature is structurally blind to |
| **G7** | **exon skipping** (A-B-C and A-C) | two junction edges out of A's last node; ⚠ the undirected graph now has a **cycle** — see `PATH_MARGINALIZATION.md` |
| **G8** | **mutually exclusive exons** | disjoint junction pairs; no spurious edge between the two alternatives |
| **G9** | **retained intron** (one isoform splices, one does not) | the intron node carries `exon_s` **and** `intron_s`; the junction edge spans it |
| **G10** | ⭐ **overlapping transcripts on OPPOSITE strands** | AMBIG signature nodes; junction edges keep distinct `strand`; no strand leakage in flags |
| **G11** | **nested / contained transcript** (one wholly inside another's intron) | correct signature; both walkable |
| **G12** | **shared exon endpoint across transcripts** | exactly one cut, one edge |
| **G13** | **1 bp node** (two events 1 bp apart) | node with `length == 1` is emitted and walkable |
| **G14** | **adjacent exons with no intron** (bookended) | contiguous edge, no junction edge |
| **G15** | **transcript at reference position 0** and **at `ref_length`** | no zero-length node, no duplicate cut |
| **G16** | **reference with no transcripts** | 1 node, 0 edges |
| **G17** | **two references** | node/edge ids contiguous per ref; no cross-reference edge |
| **G18** | **coincident opposite-strand junctions** (same donor/acceptor, both strands) | two junction edges, distinct `strand`; the v6 `sj_strand == 3` case never needed |

### 7.2 Property tests (`hypothesis`-style, or seeded random GTFs)

* **P1** — random transcript sets: I1–I12 all hold.
* **P2** — the cut set of the v8 builder equals the cut set of the v7 builder (`np.union1d` of v7 region
  starts/ends). ⭐ **This is the migration gate: v8 must differ from v7 only by the removed merge.**
* **P3** — merging v8's adjacent equal-signature nodes reproduces `v7 regions.feather` **exactly**. This is
  the strongest possible statement that nothing else changed, and it is cheap.
* **P4** — determinism: build twice, assert byte-identical feathers (T-D1).
* **P5** — every real transcript walks (I11) on a full human build, not just toys.

### 7.3 Regression / integration

* **R1** — `build_region_partition_arrays`-equivalent output feeds the C++ scanner unchanged (§9).
* **R2** — full pipeline smoke on `tests/scenarios_aligned/` with a v8 index.
* **R3** — the human build: run I1–I12 + P5 on `refs/rigel_index_v7`'s GTF and record the §3.4 census as an
  asserted regression (node count, junction count, 1 bp node count).
* **R4** — goldens: **LAST**, and expected to move (the partition changes).

---

## 8. PERFORMANCE BUDGET (acceptance numbers, measured on human GENCODE)

| stage | budget | how |
|---|---|---|
| build nodes + edges | **< 10 s**, **< 1 GB** peak | vectorised §5; no `dict` rows |
| `validate_graph` (I1–I12 incl. the I11 walk) | **< 5 s** | vectorised; I11 is `O(total exons)` |
| write feathers | < 3 s | LZ4 |
| load + CSR build | **< 1 s** | two reads + three `O(n)` passes |
| ⚠ **BAM scan deposit** | **no measurable regression** | see below |

⚠ **The scan is a budget item, not just the build.** The partition is 32 % finer, so `region_of_pos`
binary-searches a larger cut array and every crossing fragment touches more nodes. On a 992 K-node
reference the cut array is ~8 MB — well past L2 — so the per-fragment walk's cache behaviour changes.
This is a **required acceptance measurement on real cfRNA before X5 lands**, not an afterthought.

⚠ Profile on the **real** GTF. `PERFORMANCE.md` §1's standing lesson — the 10 Mb synthetic ranks costs
backwards — applies to the builder too: the toy has 1,698 regions and will not show a `dict`-per-row
regression at all.

---

## 9. WHAT DOWNSTREAM MUST CHANGE, AND WHAT MUST NOT

**Must not change in this phase.** `transcripts.feather`, `intervals.feather`, `sj.feather`, the resolver,
`scoring.py`, the EM, `locus.py`. This phase ships a new artifact and an adapter; it does not touch the
calibration model.

**The complete consumer set** (grepped, `src/**`): `calibration/region_arrays.py`, `index.py`,
`pipeline.py`, `calibration/{substrate,priors,capture_eff_length,track,calibrate,node_chain,signature}.py`.
That is the entire blast radius — ten files, all inside calibration + the pipeline wiring, none in the
resolver / scorer / EM. `track.py` uses coordinates only; `locus.py`'s hit is a comment.
⚠ `scripts/debug/**` has ~10 more consumers; they are diagnostics and may lag, but
`pass0_oracle_bench.py` and `calib_cache.py` are on the critical path for every A/B and must be migrated
with the rest.

**Adapter (the whole point of shipping this phase alone).** `build_region_partition_arrays(index)` returns
`(boundary_positions, ref_pos_offsets, region_types)` to `BamScanner.set_regions`. From v8:
`boundary_positions` = the node cut array per reference, `region_types` = `coarse_type_array(signature)`.
**The C++ ABI is unchanged.** So the scanner keeps working on the v8 graph with no C++ change — it simply
sees a finer partition.

⚠ **Consequence, and it must be stated plainly: shipping v8 alone changes calibration's numbers**, because
the partition is finer (median region 305 → 145 bp). It is *not* a bit-identical arm. Two options, and the
choice belongs to the owner:

* **(a) Ship v8 + the adapter first**, accept a measured calibration delta, then land the accumulator.
  Pro: small, verifiable steps. Con: one arm with a known-degraded intermediate (finer regions on the
  *old* deposit rule move evidence from `EV_OWN` into `EV_IMPUTED`).
* **(b) ⭐ Recommended: hold v8 behind the loader until the accumulator lands**, and gate the pair
  together. The index work is still developed, tested and reviewed on its own (§7 is entirely
  self-contained and needs no accumulator), it just does not become the *default* alone.

**Rebuild cost.** Every index: `~/Downloads/rigel_runs/*/rigel_index`, `refs/rigel_index_v7`, the four cfRNA
payload caches, every scenario suite, and the test fixtures. ⚠ `_selfsolve_cache` payloads are keyed to an
index — **namespace `/tmp/rigel_selfsolve` by index hash as part of this phase**; it is one line and it has
already cost two full rebuilds.

---

## 10. IMPLEMENTATION PHASES

```
X0  SPEC LOCK — this document reviewed; §3.3 degenerate-case table agreed
X1  TESTS FIRST — tests/calibration/test_splice_graph.py with G1–G18 + P1–P5,
    written to the target behaviour so they FAIL against today's code.
    Gate: the whole matrix runs and fails for the right reasons.
X2  BUILDER — calibration/splice_graph.py: build_splice_graph(transcripts, ref_lengths)
    -> (nodes_df, edges_df); vectorised per §5.
    Gate: G1–G18 green; P2/P3 green (v8 differs from v7 ONLY by the merge)
X3  VALIDATORS — validate_graph() with I1–I12; wired at build AND load
    Gate: P1, P4, P5 green; human build R3 green
X4  SCHEMA + LOADER — nodes/edges feather, INDEX_FORMAT_VERSION 8, CSR at load,
    manifest counts; regions/boundaries removed
    Gate: §8 budgets met on the human GTF; R1 green
X5  ADAPTER — build_region_partition_arrays from the graph; C++ ABI untouched
    Gate: R2 green; scanner runs end to end on a v8 index
X6  REBUILD every index + namespace /tmp/rigel_selfsolve by index hash
X7  hand off to the accumulator (05_accumulator_v5.md)
```

**Standing gates on every phase:** ruff; full test suite; no C++ change; determinism (P4) re-checked after
any builder edit.

---

## 11. OPEN QUESTIONS

1. **Do nRNA/synthetic transcript rows stay excluded from the event set?** They are today. Excluding them
   keeps the cut set identical to v7's (P2), which is what makes the migration checkable. But `P1G_SCOPE`
   §5 measured that 211 of 211 signature false-positives were nRNA-span edges — worth re-checking once the
   graph exists. **Recommend: keep excluded, revisit with a measurement.**
2. **Should JUNCTION edges carry the observed motif strand as well as the annotated strand?** The
   accumulator observes a motif at deposit today (`boundary_junction_strand`). With annotated junctions the
   two should agree; a disagreement is a QC signal worth surfacing rather than silently discarding.
   **Recommend: annotated strand in the index, observed strand compared in the accumulator, mismatch
   counted and reported.**
3. **Unannotated junctions.** A fragment may cross a junction not in the annotation. Today it deposits as
   spliced mass at whatever boundary it touches. On a graph it has **no edge to walk**. This is a real hole
   and it belongs to the accumulator design — `03_path_accumulator.md` §7 must answer it, and the answer
   affects whether `edges.feather` needs to grow at runtime.
4. **1 bp and sub-FL nodes** — 8.5 % of nodes are < 10 bp. Nothing in the graph objects, but every
   downstream divisor must tolerate `length = 1`. Tracked as an invariant, not a special case.
5. ⚠ **NASCENT RNA HAS NO WELL-DEFINED 3′ END, so `reach_V` is not the TES.** `02_redesign_derivation.md`
   §2.2 defines the nascent support as "the pre-mRNA span" and its reach as the distance to that span's
   end. But a nascent molecule is a *partially transcribed* pre-mRNA — its 3′ end is wherever polymerase
   had reached, which is a distribution, not a point. Three candidates: (i) `reach_V = ∞` inside the gene
   body (no taper — simplest, and defensible because nascent fragments are not size-selected against a
   transcript end); (ii) taper at the TES like mature; (iii) model a 5′→3′ nascent gradient. **This is
   undecided and it changes `E_V` near every TES — exactly the region P1g is about.** Recommend (i) with
   an explicit note, and a measurement against the oracle on `nrna_present` scenarios.

---

## 12. ⭐ A SEQUENCING REFINEMENT — take the P1g win off the index alone

The three-document plan as first written lands the index, then the accumulator, then the calibration
migration, and **nothing measurable improves until the last of those.** That is a long unrewarded stretch on
a branch that already has an uncommitted change set.

But the P1g prize does **not** need the accumulator change at all. `gdna_reframe_terminus.md`'s 190× mode error and
`ω_graft`'s ≥30× split are both fixed by the **structural flags alone** (§2.3), which this document delivers
at X4. So:

```
X0-X6   the graph index, as specified                       (foundation, no behaviour change yet)
X6b ⭐   C1/C2/C3 from P1G_SCOPE, on the new graph:
          C1  the reframe at terminus seams  (prize: capture-OFF -0.0178 / -0.0124)
          C2  omega_graft re-fitted per structural class
          C3  accept_l/accept_r -> rna_crosses_contiguously (structural, coverage-independent)
        gate: the pre-registered A/B of P1G_SCOPE §8, incl. its falsification test
        (junction-only edges must be UNMOVED — they are already exact at 1.0x)
X7      hand off to the accumulator (05_accumulator_v5.md)
```

⚠ The partition also changes at X4 (median region 305 → 145 bp), so X6b's A/B measures **both** effects at
once unless the two are separated. Separate them: record the baseline on the v8 partition *before* wiring
C1, so the arm varies one thing.

**This is the recommended order.** It gets the highest-value measured win out of the foundation work
immediately, and it de-risks the accumulator by proving the graph is correct under a real solver first.
