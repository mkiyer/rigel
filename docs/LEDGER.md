# Ledger

**One row per step, appended as it lands, never retroactively.** Plan: `IMPLEMENTATION_PLAN.md`.
Design: `ACCUMULATOR_DESIGN.md`.

The point of this file is attribution. The accumulator rework moves several entangled quantities, and a
delta is only attributable if it is recorded against a baseline taken from the same tree in the same
session.

**Standing baseline** (32-condition `ambig_dense_10mb`, mass-weighted mean absolute error vs the oracle,
`OMP_NUM_THREADS=1`): **r0 `0.079005` / r3 `0.046675`**, recorded from the working tree at `3c293038`
and reproduced 32/32 exactly at every step below.

---

## S1 — the index: reach semantics + the junction CSR (2026-07-29)

**One thing varied:** what a **contiguous** edge's four `reach_*` columns mean, plus a new read-only
accessor. Nothing in `src/` reads either, so the gate is bit-identity — and it held.

| gate | result |
|---|---|
| bench, both refits | ⭐ **32/32 flat.** r0 `0.079005`, r3 `0.046675` vs arm `w2a_pre` |
| suite | **1287 pass** (from 1280: one reach test became two, six new junction-CSR tests) |
| ruff + `ruff format` | clean |
| human-scale graph gate | build **3.1 s** / peak RSS 1.63 GB; `validate_graph` I1–I13 **2.8 s** (§8 budget < 5 s) |
| census | unchanged — 1,043,881 nodes / 1,043,595 contiguous / 404,168 junction edges |

### What changed, and why

**1. `_contiguous_reaches` now measures the GENOMIC distance to a covering transcript's span ends**,
per strand, maximised independently per side (D2). It measured the *exonic* distance over *mature*
transcripts, which is nonzero only where an exon strictly contains the position — i.e. **zero throughout
every intron.**

Since nascent RNA is an ordinary transcript that happens to be single-exon and to span its gene, and a
genomic distance is never shorter than the exonic distance inside the same span, the widest RNA molecule
covering a position is the nascent one. The old columns therefore declared **zero RNA opportunity across
every intron in the genome** — backwards, because the intron is exactly where nascent RNA lives.

Verified directly on a fixture with exons `[500,700)` and `[1200,1500)`:

| position | 50 | 100 | 200 | 300 | 450 | 500 | 650 | 700 | 800 |
|---|---|---|---|---|---|---|---|---|---|
| region | outside | exon | **intron** | **intron** | **intron** | exon | exon | outside | outside |
| `reach_lo` | 0 | 0 | 100 | 200 | 350 | 400 | 550 | 600 | 0 |
| `reach_hi` | 0 | 600 | 500 | 400 | 250 | 200 | 50 | 0 | 0 |

**2. ⚠ A JUNCTION edge's reach is deliberately NOT changed.** A junction edge is used only by a molecule
that spliced across it, so what remains either side is exonic. `_junction_edges` stays on the exonic
reach. The three junction-reach tests passing unmodified is the check on that split.

**3. New: `build_junction_edge_arrays(index) -> JunctionEdgeArrays`** — the lookup the deposit rule needs
to answer *"is this observed intron an annotated junction, and which edge?"* in the hot loop.

A CSR keyed by the **donor cut index**, which is the index the deposit already computes while locating
the lines its path crosses. It is cheap because of a measured property: every annotated intron has both
endpoints as partition cuts, so "is it annotated?" reduces to "are both endpoints cuts, and is the pair
registered?". If the start is not a cut the table is never consulted.

⭐ **Verified against `sj.feather`** — a derivation that never passes through `edges_df`'s node ids — and
the two agree exactly on both indexes: **537/537** (`ambig_dense_10mb`) and **404,168/404,168** (human).

⭐ **Hot-loop cost confirmed:** 70.4 % of cuts have no outgoing junction at all (identical on both
indexes); over the cuts that do, mean fan-out **1.31**, max 25. A one-to-two iteration scan.

### ⚠ Outstanding — the indexes on disk carry STALE reach

The change is in the *builder*. All 8 existing indexes still hold the old exonic reach in
`edges.feather`. Nothing reads those columns, so nothing is wrong today — but **every index must be
rebuilt before S5/S6 makes reach load-bearing**, and a rebuilt index will not compare byte-for-byte with
an old one on `edges.feather`. `partition_hash` will not notice: it covers `nodes.feather` only, and
nodes are unchanged.

### Files touched (S1)

* `src/rigel/calibration/splice_graph.py` — `_contiguous_reaches` rewritten; `JunctionEdgeArrays` +
  `build_junction_edge_arrays` added; `dataclass` import; `__all__`; the module docstring's "Reaches"
  section, which described the old semantics and would otherwise have become the "two docstrings
  disagreeing about one quantity" trap.
* `tests/calibration/test_splice_graph.py` — `test_reach_is_zero_where_no_mature_crosses` replaced by
  `test_contiguous_reach_is_NONZERO_INSIDE_AN_INTRON` and
  `test_reach_is_zero_outside_a_span_and_on_a_strand_with_no_transcript`.
* `tests/calibration/test_junction_edge_arrays.py` — new, 6 tests.

---

## S2 — the reference accumulator, the spec matrix, and the scan profiler (2026-07-29)

**Written tests first and verified failing** (`ImportError` — the API did not exist), then implemented
until green. The reference is the executable specification; the native accumulator will be gated on
byte-identity to it.

| gate | result |
|---|---|
| spec matrix | **33 pass** (`tests/native/test_accumulator_spec.py`) |
| suite | **1280 pass** |
| ruff + `ruff format` | clean |
| ⭐ `E[Σ 1/(L−1)] = ρ` **with no length model supplied** | within [0.98, 1.02] at **50 / 200 / 1000 bp** node spacing |
| ⭐ `E[count] = ρ·(E[L] − 1)` | within [0.98, 1.02] |
| ⭐ order-independence | the tally is byte-identical under a random permutation of arrival order |

### ⛔ A real bug the matrix caught, in the first draft of the reference

*Contained* was implemented as **"crossed no line"**. That is wrong: an unannotated intron can swallow
every line between two blocks, leaving a fragment that crosses nothing yet **straddles two nodes** — and
the code then credited its whole length to the first of them.

**Contained now means the whole path lies inside ONE node.** Such a straddling fragment deposits on **no
object** but still increments `node_start_count`, so the loss is visible rather than silent. Pinned by
`test_a_fragment_straddling_two_nodes_without_crossing_a_line_is_NOT_contained`.

Two further errors, both mine and both in the *tests*, caught the same way: a claim that block
`[120,200)` crosses the line at 100 (it is to the left of the block), and a fixture where
`deposit(0, 150, 300)` crosses **two** lines, not one.

### ⭐ THE SCAN PROFILER — and why the obvious number is wrong

`scripts/design/scan_profile.py`. There was **no scan-profiling harness in the repository**; the only
recorded figure was *"the scan is ~2 % of runtime"*, measured on a 39 MB BAM. The deposit is isolated by
scanning each BAM twice — with and without `set_regions` — so everything else cancels.

⚠ **Wiring the accumulator adds two costs with different scaling, and they must not be averaged.**
`O(fragments)` is the deposit; `O(partition)` is allocating and zeroing a per-worker accumulator, merging
the workers, and copying per-reference structs into the flat payload. Against the 1.04 M-node human
partition the second dominates on a shallow BAM, so a single-BAM `delta / n_fragments` is **not** a
per-fragment cost. Measured, same index, three real cfRNA libraries:

| library | fragments | whole scan | accumulator | naive delta/fragment |
|---|---|---|---|---|
| LBX0190 | 155,352 | 1.036 s | 0.143 s | **917.8 ns** |
| MO_3021 | 875,670 | 4.485 s | 0.476 s | **543.1 ns** |
| LBX0588 | 1,299,234 | 2.896 s | 0.536 s | **412.9 ns** |

The naive figure falls as depth rises — it is the fixed partition cost being amortised. Regressing
`delta(n) = fixed + c·n`:

* ⭐ **c ≈ 357 ns/fragment** — the deposit, the number that generalises
* ⭐ **fixed ≈ 0.108 s** — the price of the 1.04 M-node partition itself

⚠ Three points with visible scatter, so treat `c` as **~350–400 ns/fragment**, not 356.9. And note it is
~2.9× the 124 ns that a pure C++ microbenchmark of `Accumulator::deposit` reports: the A/B additionally
contains `fragment_genomic_spans` and the scan-side gate in `deposit_to_accumulator`. **The S3 budget is
the ~357 ns end-to-end figure**, because that is the one a user pays.

Projected single-thread deposit cost, partition excluded: **3.6 s at 10⁷ fragments, 35.7 s at 10⁸.** That
is the regime this tool is for, and it is why the per-fragment allocation in the current deposit (22.8 ns,
18 %) is worth removing.

### ⚠ Outstanding in S2

The reference has **not yet been run on a real cfRNA payload through a shim** — the other stated S2 gate.
It is currently exercised only on hand-built fixtures of 2–6 nodes, so a scaling defect would not have
been caught. This should close before the C++ is written against it.

### Files touched (S2)

* `tests/native/_accumulator_reference.py` — rewritten in place (`Partition`, `Tally`, `Accumulator`,
  `density_quantum`, `FragmentPool`, `DepositOutcome`). No version suffix; the old fractional reference
  is gone.
* `tests/native/test_accumulator_spec.py` — rewritten in place, 33 tests.
* `scripts/design/scan_profile.py` — new.

⚠ The shipped C++ accumulator now has **no reference test** until S3 replaces it (the old spec matrix
tested the fractional rule, which no longer exists). The oracle bench is the only guard in that window.
It is one step, on code scheduled for deletion, and it is stated here rather than discovered later.

### ⭐ S2 CLOSED — the reference validated on real data at full human scale

`scripts/design/reference_on_real_data.py`. The spec matrix exercises fixtures of 2–6 nodes; this drives
the reference from a real cfRNA BAM against the **1,043,881-node / 404,168-junction** human partition,
and re-derives every reported total **by a different method** — the reference locates crossed lines with
two `searchsorted` calls and takes the index range, this walks each segment's cuts with `bisect` one at a
time.

60,000 read groups of LBX0190 → 57,393 accepted:

| cross-check | reference | independent re-derivation |
|---|---|---|
| `Σ start_count == accepted` | 57,393 | 57,393 |
| edge crossings | 54,051 | 54,051 |
| **edge density (fixed point)** | **1,165,766,162,679** | **1,165,766,162,679** |
| junction crossings | 41,133 | 41,133 |

The density total matching *exactly as an integer* is the strongest of the four — it depends on the
crossing set, the quantum, and the rounding mode all being right simultaneously.

Sanity, all consistent with figures measured independently earlier: `dropped_too_long` **4.3 %** (against
5.45 % from `bam_spans.py` on the full library); unannotated introns **2.5 %** of intron observations
(against 2.33 % measured from CIGAR ops); **36.6 %** of fragments contained, as expected with median
`L` = 157 against a median node of 151.

### ⭐⭐ AND THE POOLS SEPARATE BY LENGTH, MEASURED WITH NO MODEL

| pool | fragments | mean `L` |
|---|---|---|
| `DNA_INTERGENIC` | 2,260 | **88.0** |
| `DNA_INTRONIC` | 4,432 | **88.5** |
| `DNA_INTRON_EXON` (splash) | 884 | 138.8 |
| `DNA_INTERGENIC_EXON` (splash) | 65 | 211.6 |
| `RNA_SPLICED` | 27,328 | **220.7** |

⭐ The two pure gDNA pools agree with each other at ~88 bp and the pure RNA pool sits at ~221 bp — a
**2.5× separation, and `μ_g − μ_r` is the determinant of the 2×2** that separates gDNA from RNA by length
alone. This is the design's central claim, measured directly from structurally pure populations with no
fitted model anywhere.

⚠ **And it disagrees with what ships today.** The shipped `gdna_fl_pmf` has mean **146.05** on this
library, because it sums four differently-tilted pools (intergenic and intronic × contained and
boundary). The pure intergenic-contained pool says **88.0**. An independent measurement during the code
recon put the corrected pure estimate at **86.49** — so the new pool axis is not a re-labelling, it moves
the gDNA mean length by ~40 %, and the old value was biased long by pooling the boundary compartments.

The two splash pools behaving *between* the two extremes (139, 212) is the expected signature of
on-target gDNA plus leaked RNA, and is why they are named pools rather than folded in.

⚠ The reference costs **503 µs/fragment** (pure Python) — ~1,400× the native path. Expected for a
specification; it is why the shim takes a `--limit`.

### ⛔ S2 REVIEW — 5 adversarial lenses, 4 blocking findings, all reproduced and fixed

The reference is the only gate the C++ will have, so it was reviewed before a line of C++ was typed.
**Verdict: BUGGY.** Every finding below was reproduced by running code, then fixed, then pinned.

**B1 — two contradictory formulas for `L`, 43 lines apart.** `length = (hi−lo) − Σ(e−s)` sums intron
lengths; `_segments` takes their **union**. They disagree the moment introns overlap:

| introns on `[150,500)` | naive `L` | true (segment union) |
|---|---|---|
| `[(210,260),(240,300)]` | 240 | **260** |
| `[(200,400),(250,300)]` | 100 | **150** |
| `[(150,480),(160,470)]` | **−290** | **20** |

A negative `L` filed a good fragment as `dropped_empty` — invisible to the start-count invariant, since a
rejected fragment never reaches it. **Live on real data**: MO_3021 `chr8:138206290-138206943`, two mates
disagreeing about an acceptor. 1 in 875,670.

⭐ Blocking far beyond that rate: **a C++ author builds the segments first and will naturally sum them**,
so byte-identity — the only gate — would have reduced to which of two contradictory rules the file
happened to carry. **Fix:** `_normalise_introns` returns a sorted disjoint set, `L` is *defined* as
`Σ segment lengths`, and the absorbed count is a QC field. One definition, no second formula possible.

⚠ **Strict overlap only.** Two introns that merely *abut* are two splices with a zero-length exon
between them, and a fragment can legitimately use both — merging them destroyed a junction silently
(caught by `test_a_fragment_using_TWO_junctions_credits_BOTH` going red).

**B2 — one array holding two strand conventions.** `node_spanning_*` was indexed by the same channel as
the edge banks, i.e. **sense-relative-to-motif for spliced fragments and genome strand for unspliced
ones**. The design's justification for having no spliced node channel ("an annotated-spliced fragment can
never be *contained*") is true but covers only containment; a spliced fragment's blocks routinely **span**
a node. ⭐ Measured: **65–69 %** of all `node_spanning` deposits on real cfRNA come from spliced
fragments, **40–44 %** of them landing in the opposite column from their unspliced neighbours.
**Fix:** node banks are **always** genome strand. Edge and junction banks keep the spliced/sense
convention. Pinned by two tests, including the discriminating case (a spliced fragment *sense* to a
`−` junction must still book its node as genome-minus).

**B3 — the 33-test matrix could not tell the reference from a materially different one.** The reviewer
copied it, changed three behaviours, and the unmodified matrix passed **33/33**. Now 45 tests, and the
three behaviours are each pinned.

**B4 — a third bug, found while fixing the others.** The contained test and `start_count` used the
fragment's *extent* (`lo`, `hi−1`) rather than the path's first and last **covered** base. With a leading
intron the molecule does not begin at `lo`: `introns=[(150,480)]` over `[150,500)` has its path entirely
in a different node, and the invariant was being credited to a node the fragment never touches.

**Also fixed:** `build_junction_edge_arrays` sorted `(acceptor, donor)` while the reference's
`Partition.from_cuts` sorted `(strand, acceptor, donor)`. The junction id **is** the rank in that order,
so the two disagreeing would permute every row and break byte-identity — reachable only via a
strand-coincident junction pair, which is exactly what one test constructs.

| gate after the fixes | result |
|---|---|
| spec matrix | **45 pass** (from 33) |
| suite | **1292 pass** |
| ruff + format | clean |
| real-data shim, LBX0190 + MO_3021 | all four cross-checks exact on both; `dropped_empty` now 0, `introns_absorbed` 2 |

⭐ MO_3021's pools, second library, confirming the separation is real and sample-specific:
gDNA **118.4 / 124.1**, RNA **203.3** — a 1.7× separation where LBX0190 gave 2.5×.

### ⚠ Still open from the review, for the plan rather than the code

* `IMPLEMENTATION_PLAN.md` §3.2's pseudocode still encodes CONTAINED as `crossed == 0` and computes `L`
  from the **unclipped** span — both now wrong. Correct the plan before the C++ follows it.
* §3.3's `junction_edge_id_` must be documented as the **dense CSR rank**, not `edges_df`'s row (which
  reaches 1,447,755 for a 404,168-row bank).
* The `L == 1` guard on `quantum(L−1)` is load-bearing and the plan calls it unguarded — integer division
  by zero. `deposit(0, 500, 501)` reaches it.
* `ACCUMULATOR_DESIGN.md` §12.A still lists the minimum intron length as blocking; the owner's ruling
  (every CIGAR `N` is an intron, no floor) needs recording there.

### S2 — two owner rulings landed after the review (2026-07-29)

**1. Zero-length exons are impossible, so ABUTTING introns are malformed.** This REVERSES the strict-overlap
decision taken during the review fix. Two introns sharing an endpoint imply a zero-length exon between
them — a transcript with one is molecularly identical to a transcript without it, so a single molecule can
never legitimately use both, and the pair is an alignment artifact. `_normalise_introns` merges them and
counts the merge.

⭐ The index cannot produce the pair either, so the rule holds on both sides: a zero-length exon is dropped
when `_Exons` is built, which fuses its two flanking introns into one. Pinned by
`test_ABUTTING_introns_are_MALFORMED_and_merge`; the test that legitimately uses two junctions now has a
real exon between them and its own partition.

**2. ⭐ ONE STRAND CONVENTION, and it deleted a concept rather than fixing one.** Owner ruling: store
everything relative to the **genome** strand; *sense* / *antisense* name the transcript-relative notion and
are never stored.

```python
CHANNEL_PLUS, CHANNEL_MINUS = 0, 1        # every count and density array, without exception
def genome_channel(strand): ...
```

Sense is **derived**: a junction edge carries its own genomic strand, so a consumer computes
`sense = (fragment strand == junction strand)`. Nothing is lost — including the wrong-stranded-junction
aligner diagnostic — and the `spliced_primary` branch is gone entirely. `motif_strand` is now needed *only*
to disambiguate strand-coincident junctions in the lookup, never to select a channel.

This is what made B2 possible: a schema with two conventions in it. There is now one.

| gate | result |
|---|---|
| spec matrix | **46 pass** |
| suite | **1293 pass** |
| ruff + `ruff format` | clean |
| real-data shim, LBX0190 + MO_3021 | all four cross-checks exact on both |

### ⛔ CONSEQUENCE — `ACCUMULATOR_DESIGN.md` §5 is now STALE

Its channel table still says `spliced × {sense, antisense}`. That is the schema the C++ implements, so it
must be rewritten before S3. Listed as the first task in `IMPLEMENTATION_PLAN.md`'s state block.

---

## S2.1 — the doc corrections (2026-07-29)

**Documentation only; no code touched, so the suite (1293) and the bench are untouched by construction.**
The plan and the design are what the C++ follows, and statements in them contradicted the reference — which
is the C++'s only gate. Four were known and listed in the plan's state block; **three more of the same class
surfaced while making them.**

| | correction | was | now |
|---|---|---|---|
| 1 | design §5.1 channel axis | `spliced × {sense, antisense}` | `{genome +, genome −}`, every bank, no exception |
| 2 | plan §3.2 deposit pseudocode | four bugs (below) | a transcription of the reference |
| 3 | plan §3.3 junction id | `junction_edge_id_` → `edges_df` row | the **CSR slot**; the third array deleted |
| 4 | plan §3.2 `quantum(L−1)` | unguarded | guarded, with its residue stated |
| **5** | design §8 pool count | "two", then "a third" in §8.2 | **five**, with the measured means |
| **6** | plan S3 performance gate | "no regression against 124 ns" | **~357 ns/fragment end-to-end**, and the two are distinguished |
| **7** | `CARRY_FORWARD.md` §4 | six bullets contradicting the settled deposit rule | banner + `⛔ SUPERSEDED` in place, each with its replacement |

**On (1):** the rewrite records what the channel physically *is*, which no document had stated —
`align_strand` is the OR of the exon blocks' strands with read 1 keeping its reference orientation and read
2 flipped to agree (`bam_scanner.cpp:661-686`), i.e. the fragment's **read-1 genome orientation**. That is
protocol-independent, which is exactly why declared strandedness is a *header field* and not a channel
convention: whether read-1 orientation equals the RNA molecule's strand is a property of the library.
⭐ The new schema **cannot express** the shipped `primary = (align_strand == motif_strand)` bug, because no
channel is labelled *sense*, so no label can be wrong.

**On (2):** the four bugs are each attributed in place — `crossed == 0` for CONTAINED · `L` from the
**unclipped** span *by a second formula* · `lo`/`hi−1` instead of the path's first and last **covered**
base · unguarded `quantum(L−1)`. Plus an ordering constraint the old block violated: **the junction lookup
must run before the crossing loop**, because `spliced` selects which edge bank the crossings land in.

**On (3) — measured, and it resolved more cleanly than the state block anticipated.** The block asked for
`junction_edge_id_` to be *documented* as the dense CSR rank. It should be **deleted**: the reference
indexes `junction_count[k]` by the CSR slot directly, so the id *is* the slot and no indirection exists.

⭐ Verified on both indexes — human and `ambig_dense_10mb`:

| | human | `ambig_dense_10mb` |
|---|---|---|
| junctions | 404,168 | 537 |
| `edges_df` rows / highest junction row | 1,447,763 / **1,447,755** | 2,234 / 2,230 |
| CSR slot == dense rank within `flatnonzero(kind == JUNCTION)` | **404,168 / 404,168** | **537 / 537** |
| `edge_row` strictly ascending | yes | yes |

So `JunctionEdgeArrays.edge_row` is a monotone map **CSR slot → `edges_df` row**, needed only to join a
payload row back to the index, and it **does not cross the ABI**. Writing `junction[edge_row]` would run
1.04 M rows off the end of a 404,168-row bank. ⚠ `edge_row`'s docstring in `splice_graph.py` still calls
itself *"the junction-edge id"* — the one place a C++ author would look for the definition. **Fix in S3.**

**On (5):** ⭐ the design asserted two pure pools, then §8.2 asked for "a third", while the settled ruling
is five. §8 now carries the five-pool table with the S2 measurements and the consequence that ships today:
`gdna_fl_pmf`'s mean is **146.05** on LBX0190 against **88.0** for the pure intergenic-contained pool, so
the new axis moves the gDNA mean length by **~40 %** rather than re-labelling it.

**On (6):** the plan named the S3 performance gate in two places with two different numbers — 124 ns (a
pure C++ microbenchmark of `deposit`) and ~357 ns (end-to-end, from the S2 regression). Both are correct
measurements of different things; only one is a gate. §6 now separates them and states that the gate is the
end-to-end figure, *because that is the one a user pays*. §6 also still claimed **"there is no scan
profiling harness in this repo"** — S2 built it.

**On (7) — the worst of the three, and it was hiding in plain sight.** `CARRY_FORWARD.md` is one of the four
documents in the reading order, and its §4 is a **pre-settlement** proposal list distilled from the 278
deleted docs. Six of its bullets tell a C++ author to do something the design subsequently rejected:

| §4 said | settled rule |
|---|---|
| `+1/L` on each crossed edge | **`1/(L−1)`** — a 0-bp line has `L−1` placements; `1/L` reads 0.54 % low |
| "crossing none, it is contained" | **the whole path lies inside one node** — the exact bug the S2 matrix caught |
| "`L` is the fragment's PATH length — **exonic bases only**" | `L` = span − introns, so the **mate gap counts**; exonic-only collapses the distribution to a spike at 2× read length |
| unannotated junctions deposit "in the **spliced** channel" | **unspliced** (owner ruling) — they are disproportionately artifactual, so they must compete with gDNA |
| "three exact integer conservation identities" | **tautologies, rejected** (design §10.2); replaced by `Σ start_count == n_accepted` |
| 1 bp node spanned → "node untouched" | **its `spanning` counter is incremented**; only `contained` is untouched |

⚠ Four of the six are in the *same* bullet list a C++ author would read as the deposit rule, and the first
three are individually sufficient to produce a wrong tally that byte-identity could not catch — because
byte-identity is against the reference, and someone who implemented §4 would not have consulted it.

§4 now opens with a banner stating that it is history and that the design plus the reference are the
authority, and every superseded item is marked in place with what replaced it. Nothing was deleted: the
rejected ideas and the reasons they were rejected are worth more together than either alone.

Also landed: design **§12.A** turned from "the one item blocking implementation" into the two owner rulings
(no minimum intron length; abutting introns merge), the artifact catalogue gained rows for overlapping and
abutting introns, §7.2's dangling `§12.A` pointer was repaired, and the plan's duplicate `### 3.4` heading
was renumbered (pools §3.4, payload §3.5, reference §3.6).

**Nothing now blocks implementation.** Next: S3.

---

## S2.2 — the `edge_row` docstring, and the ordering contract it turned out nothing tested (2026-07-29)

**One thing varied:** the junction-id semantics, in the one place a C++ author would look them up. What was
meant to be a one-line docstring fix uncovered a missing gate.

| gate | result |
|---|---|
| suite | **1298 pass** (from 1293; +5 tests) |
| ruff + `ruff format` | clean |
| bench | not re-run — no `src/` behaviour changed, only a docstring |

### The docstring

`JunctionEdgeArrays.edge_row` described itself as *"row in index.edges_df, **i.e. the junction-edge id**"*.
It is not the id — the id is the **CSR slot** — and using it to index a junction bank writes 1.04 M rows
past the end of a 404,168-entry array. It is now labelled a **join key**, with the ordering contract stated
beside it. `_lookup` in the test module, which advertises itself as the deposit's worked example, returned
`edge_row[k]`; it now returns `k`.

### ⛔ And the contract had never been tested

The ledger's S2 entry records that `build_junction_edge_arrays` and the reference's `Partition.from_cuts`
disagreed once on sort key, and that agreement is what the S3 byte-identity gate rests on. **Nothing tested
it.** The spec matrix exercises only `from_cuts`; `reference_on_real_data.py` builds its `Partition`
*straight from* `build_junction_edge_arrays` (`scripts/design/reference_on_real_data.py:59-68`), so the two
sorts were never compared to each other. Now they are, by a route that shares only the definition of a
junction: the builder walks node ids and applies a per-reference `cut_base − node_base` shift, the check
names each junction genomically and lets `from_cuts` resolve both endpoints with its own `searchsorted`.

### ⭐ The instructive part: the first two versions of that test had NO TEETH

Written, green, and worthless — the B3 failure again, caught this time by perturbing the code rather than
trusting the green.

| perturbation of the builder's key | v1 fixtures | + nested fixture |
|---|---|---|
| `strand` promoted to primary | **passes** ⛔ | **fails** ✅ |
| donor/acceptor priority swapped | **passes** ⛔ | **fails** ✅ |
| `strand` dropped from the key | **passes** ⛔ | **passes** — see below |

⭐ **Why:** on every fixture in the module *and on both real indexes*, donor order, acceptor order and
strand order happen to agree, so every permutation of the sort key yields the identical answer. The fix is
one fixture — a **nested** intron pair with the **outer** intron on the minus strand, `[400,1400)` NEG
enclosing `[600,1000)` POS. Nesting breaks donor-vs-acceptor (`400 < 600` but `1400 > 1000`); the strand
assignment breaks donor-vs-strand (the smaller donor is the NEG one). Two permutations become observable
immediately.

⚠ **The third stays unobservable, and that is a finding, not a gap left open.** ⭐ Measured: `edges_df`
emits a strand-coincident pair **POS-before-NEG whichever order the GTF lists the genes in** — the graph
builder sorts — and `np.lexsort` is stable, so both routes start from an already-correct tie order. So the
**builder's** `strand` key is *defensive*: keep it, because it makes the contract explicit rather than
resting on `edges_df`'s internal sort, but no test can defend it. **`from_cuts`'s `strand` key is a
different matter and IS load-bearing**, because a caller may list junctions in any order — now pinned by
`test_a_junction_id_is_a_function_of_the_PARTITION_not_of_argument_order`, verified to fail when the key is
dropped, and only on the reversed-argument arm.

⚠ The pre-existing `test_opposite_strand_junctions_at_the_same_coordinates_are_DISTINCT_edges` does **not**
pin the strand key either: it passes the pair POS-first, so a stable sort keeps NEG at id 1 regardless. It
is a good test of a different property, and it was reasonable to think it covered this. It does not.

### Files touched (S2.2)

* `src/rigel/calibration/splice_graph.py` — `JunctionEdgeArrays` docstring: the id is the slot, `edge_row`
  is a join key, the ordering contract, and the memory-safety consequence of confusing them.
* `tests/calibration/test_junction_edge_arrays.py` — `_lookup` returns the slot; nested chr3 and a
  strand-coincident fixture added; `test_the_csr_slot_order_matches_the_reference_accumulator`
  (×2 fixtures) and `test_a_strand_coincident_pair_is_two_distinct_slots` new.
* `tests/native/test_accumulator_spec.py` — `test_a_junction_id_is_a_function_of_the_PARTITION_not_of_argument_order`
  (×2 arms) new.

---

## S2.3 — naming: one word per concept, and they are the codebase's own words (2026-07-29)

**Owner rejected four inventions in the spec. All four were mine, and one of them had CAUSED a bug** — I
coined a word, then reasoned about the coinage instead of the concept it stood for.

| ⛔ invented | ✅ now | why it was wrong |
|---|---|---|
| `genome_channel`, `CHANNEL_PLUS/MINUS` | `STRAND_COLUMNS`, `Strand.POS`/`NEG` | "channel" is leftover vocabulary from the 4-way `spliced × primary` axis this rework **deletes**. The axis is the strand |
| `N_STRANDS = 2` | `N_STRAND_COLUMNS = 2` | there are **four** strands (`POS\|NEG == AMBIGUOUS`, plus `NONE`); only two name a column |
| `motif_strand` → `sj_motif_strand` | **`sj_strand`** | already the codebase's name (`sj.feather`, `sj_map`, `SJKey`, `cr.sj_strand`). My "fix" added a *third* name for one quantity |
| `align_strand` | **`strand`** | it is the universal case; qualify the exception, not the norm |
| `lo` / `hi` for a fragment | `start` / `end` | the codebase's coordinate names; `lo`/`hi` belongs to `reach_lo_pos`, a different quantity |

### ⭐ The owner's ruling that resolved it: TWO strands, and they are independent

| | what it is | who has one | what it does |
|---|---|---|---|
| `strand` | the genomic strand the read **aligned** to | every read | selects the column, nothing else |
| `sj_strand` | a splice junction's strand, from its genomic **motif** (`GT..AG` → `+`, `CT..AC` → `−`); `XS` for STAR, `ts` for minimap2, auto-detected | spliced reads only | resolves an intron against the annotation, nothing else |

Neither constrains the other. Comparing them yields *sense* vs *antisense*, which is **derived** and never
stored. ⭐ The shipped code collapses the comparison into a third concept — a bool named `primary`
(`bam_scanner.cpp:1493`) — and that is how a dUTP library put 0.6 % of its spliced fragments in the column
labelled *sense*. `primary` is deleted in S3, not renamed.

### ⛔ AND THE RENAME EXPOSED A REAL BUG — which is the whole argument for doing it

S2.2 had made an `AMBIGUOUS` `sj_strand` match a junction *"on coordinates alone"*. **Wrong.** `sj_strand`
is the OR of a per-RECORD tag, so AMBIGUOUS means **the two mates disagreed about one molecule** —
contradictory evidence, the same family as mates agreeing in reference orientation. It must credit no
junction. The rule is three-way, and both my earlier versions had only two:

| `sj_strand` | means | rule |
|---|---|---|
| `NONE` | no strand tag in the BAM at all | **no information** → match on coordinates alone |
| `POS` / `NEG` | one definite observed strand | must agree with the junction edge's own `sj_strand` |
| `AMBIGUOUS` | the two mates' tags **disagree** | **contradictory** → trust no splice; own QC counter |

⚠ `NONE` staying permissive is load-bearing, and S2.2 had nearly broken it. Aligners differ — STAR writes
`XS`, minimap2 writes `ts`, some write neither — so on an untagged BAM **every** spliced fragment arrives
with `NONE`. Demanding a strand there deletes 100 % of that aligner's annotated junctions, and it would
present as a stale annotation rather than a convention bug. Pinned by
`test_a_MISSING_sj_strand_MATCHES_on_coordinates_alone`, which did not exist before.

⚠ AMBIGUOUS increments **`contradictory_sj_strand`**, not `unannotated_introns`. That counter measures
annotation coverage; feeding it alignment disagreements makes it report a stale annotation whenever the
aligner is inconsistent.

| gate | result |
|---|---|
| suite | **1303 pass**, 53 in the spec matrix |
| ruff + `ruff format` | clean |
| harness self-check | 6/6 |
| real-data cross-checks | all four exact |

### ⚠ Two process mistakes of mine, recorded so they are not repeated

* I ran `ruff format` over `scripts/`, reformatting ~130 unrelated files. **The project's format scope is
  `src/ tests/`**; `scripts/` is linted, not formatted. Reverted.
* A `sed`-style bulk rename left `_normalise_introns` using the *fragment's* bounds where it meant the
  *intron's* — the two are both in scope there. Hand-written instead. **Bulk-rename a file whose locals
  share a name with its parameters and you get a silent semantic change**, which is exactly the failure the
  spec matrix exists to catch and did.

### Deferred to its own arm

`align_strand` → `strand` **tree-wide**: 101 occurrences remain in `src/` and `tests/`, and **seven are
string keys** (`buffer.py:193,342,381` — a parquet column in the spill path — plus `resolve.cpp:38`,
`resolve_context.h:307,352`, `tests/test_buffer.py:45`). A string key survives compilation and fails at
runtime, on the buffer→EM path, which the accumulator rework does not touch. Bundling it into S3 would put
an unrelated runtime-failure mode inside the one step whose gate is byte-identity.

---

## S3 (WIP) — `accumulator.{h,cpp}` rewritten; on branch `s3-accumulator` (2026-07-29)

**`main` stays GREEN at `c9fc7c1c`. The rewrite is on branch `s3-accumulator` at `e1114d48`** because
`bam_scanner.cpp` still calls the old API, so the extension does not build. The two files compile
standalone (`g++ -std=c++17 -fsyntax-only`), and the deposit is a direct transcription of the spec.

### What landed

`Region` / `Boundary` / the 4-channel `spliced × primary` axis are gone. In their place:

| | size | holds |
|---|---|---|
| `Node` | 48 B | `contained_count` / `spanning_count` `uint32[2]` + both densities `uint64[2]` |
| `ContiguousEdge` | 48 B | `unspliced` / `spliced` count + density |
| `JunctionEdge` | 24 B | count + density |

`node_start_count_` is its **own flat `uint32` array**, which is what keeps `Node` at 48 B with no padding
member. `static_assert` on all three sizes.

* `kNStrandColumns = 2`, and **`strand_column()` returns −1** for a strand that names no column. ⚠ The −1 is
  the whole point: `align_strand == STRAND_POS ? 0 : 1` compiles, is shorter, and silently books an
  undefined strand into the minus column — the bug this rework exists to delete.
* `density_quantum()` — `round(2^32 / placements)`, halves away from zero, contract-pinned.
* `FragmentPath` — one fragment as the scanner has it. ⚠ Its introns **need not** be sorted, disjoint or
  de-duplicated; `deposit` normalises them, which is what lets `L` be *defined* as the segment total
  instead of computed by a second formula.
* `DepositScratch` — reusable buffers, so the deposit allocates nothing. ⭐ The shipped one spent **22.8 ns,
  18 % of the deposit**, in a per-fragment `std::vector`, invisible to a sampling profiler because the time
  is attributed to `malloc`.
* `DepositCounters` — the QC denominators, mergeable.
* `set_junctions()` as a **second** method, because `set_regions` throws if called twice.

The deposit follows the spec's order, and every ordering decision carries a comment saying why it is that
order: strand column **first** → clip to the reference → normalise introns (overlapping **and** abutting
merge) → `L` = **sum of the path's segments**, one formula → the length gate → resolve junctions **before**
the crossing loop, because `spliced` picks the bank → per-segment crossings via two binary searches over a
contiguous index range → spanning as the **interior** of that range → junction deposits → contained iff the
whole path lies in one node → the pool.

⭐ **`max_length < 1` now throws**, with a message that says why, instead of silently disabling itself as
the shipped code does (`accumulator.cpp:30-72`). It gates `L` *and* sizes the pool histograms, so at 0
every real fragment is dropped as too long and the entire tally is silently empty. Owner-confirmed.

### What remains

`IMPLEMENTATION_PLAN.md` §0 carries the ordered checklist, the four verified landmines, and the gates. In
brief: the scanner's `WorkerState`, `set_junctions` plumbing, the deposit adapter (delete `primary`, delete
the `motif_strand` alias), `build_result`'s keys, the nanobind bindings, the byte-identity gate driven
**through the bindings** rather than `AccumulatorPayload`, and the deletion of `fragment_genomic_spans`
plus its parallel spec.
