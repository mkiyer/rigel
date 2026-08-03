# Ledger archive — the accumulator rework, S1 to the deletion (2026-07-29 / 2026-07-30)

**Moved verbatim from `docs/WIP.md`, never rewritten.** The ledger is append-only because attribution is
its whole point; splitting it preserves that, editing it would not. `docs/WIP.md` keeps the index of every
entry (including these) and the full text of everything from the index rebuild onward.

⛔ Every benchmark number below refers to the deleted `ambig_dense_10mb` suite and is VOID — see the
deletion entry at the foot of this file. The accumulator results (byte-identity, worker determinism) are
statements about two implementations agreeing and remain reproducible against any index.

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

> ✅ **CLOSED 2026-07-30** — see the *Index rebuild* entry at the end of this file. ⭐ It turned out to be
> **7** distinct indexes, not 8, and **~38 % of human contiguous edges** had the wrong reach.

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

* `docs/accumulator/IMPLEMENTATION_PLAN.md` §3.2's pseudocode still encodes CONTAINED as `crossed == 0` and computes `L`
  from the **unclipped** span — both now wrong. Correct the plan before the C++ follows it.
* §3.3's `junction_edge_id_` must be documented as the **dense CSR rank**, not `edges_df`'s row (which
  reaches 1,447,755 for a 404,168-row bank).
* The `L == 1` guard on `quantum(L−1)` is load-bearing and the plan calls it unguarded — integer division
  by zero. `deposit(0, 500, 501)` reaches it.
* `docs/accumulator/DESIGN.md` §12.A still lists the minimum intron length as blocking; the owner's ruling
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

### ⛔ CONSEQUENCE — `docs/accumulator/DESIGN.md` §5 is now STALE

Its channel table still says `spliced × {sense, antisense}`. That is the schema the C++ implements, so it
must be rewritten before S3. Listed as the first task in `docs/accumulator/IMPLEMENTATION_PLAN.md`'s state block.

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
| **7** | `docs/SESSION_HANDOFF.md` §4 | six bullets contradicting the settled deposit rule | banner + `⛔ SUPERSEDED` in place, each with its replacement |

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

**On (7) — the worst of the three, and it was hiding in plain sight.** `docs/SESSION_HANDOFF.md` is one of the four
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

`docs/accumulator/IMPLEMENTATION_PLAN.md` §0 carries the ordered checklist, the four verified landmines, and the gates. In
brief: the scanner's `WorkerState`, `set_junctions` plumbing, the deposit adapter (delete `primary`, delete
the `motif_strand` alias), `build_result`'s keys, the nanobind bindings, the byte-identity gate driven
**through the bindings** rather than `AccumulatorPayload`, and the deletion of `fragment_genomic_spans`
plus its parallel spec.

---

## S3 — the C++ wired, and byte-identical to the specification (2026-07-30)

**S3 is closed.** `bam_scanner.cpp` now speaks the new API, the extension builds, and the native
accumulator reproduces `tests/native/_accumulator_reference.py` **byte for byte on real data at full human
scale**. The pipeline is red from here until S4 rewrites the payload — by design (§5's clean cut), and every
one of those failures is attributed below rather than assumed.

| gate | result |
|---|---|
| parity, fixtures | **6 pass** (`test_accumulator_native_parity.py`): 45 named cases + 10,000 random fragments + merge bit-identity at 2/4/8 shards |
| ⭐ parity, **real data** | **BYTE-IDENTICAL** on LBX0190 (60,000 fragments, 44 refs) and on the **whole of MO_3021** (875,670 fragments, 104 refs) — every array, every dtype, every QC counter |
| ⭐ **bit-identical at 1/2/4/8 workers** | pass, and every bank asserted non-empty first |
| `Σ node_start_count == deposited` | holds on a real scan |
| spec matrix | **56 pass** (from 53) |
| suite | **1048 pass**; 240 fail + 17 error, **all 257 from one line** — see below |
| ruff + `ruff format` | clean |
| bench | not run: the pipeline cannot complete a scan until S4, so there is no comparable number |

### ⭐ The gate, and the two holes it had

The parity harness was written **before** the bindings and verified red (5/5, against the old constructor —
it does not `importorskip`, so a stale extension is a hard error). It reads its field list off
`dataclasses.fields(Tally)`, so a new specification field joins the gate automatically.

⛔ **Then 23 perturbations of the C++ were run to see whether it actually bites. Two came back GREEN, and
both were real gaps.**

| green perturbation | why the gate could not see it | fixed by |
|---|---|---|
| **only the LEFTMOST annotated junction is credited** | the battery had **one** annotated junction, so no fragment could use two — and this is a rule the design deliberately REVERSED (`docs/SESSION_HANDOFF.md` §3 trap 21 still recommends the old one) | three junctions in the partition, and 8 cases using two at once |
| `node_of_pos` upper clamp removed | unreachable from `deposit` — the path is clipped to the reference first — but the method is **bound**, so a caller can reach it | `test_node_of_pos_agrees_everywhere_including_outside_the_reference` |

Both were re-run after the fix and now fail. The other 21 are caught, including: edge weight `1/L` instead
of `1/(L−1)` · CONTAINED as "crossed no line" · `L` from the raw introns (B1) · `first_base` from the
fragment extent · spanning as "both lines crossed" · the spliced bank never used · an undefined strand
booked into the minus column · truncating instead of rounding half-away-from-zero · the `sj_strand` filter
applied when NONE · abutting introns not merging · the node quantum using `L−1` · outcome precedence ·
introns filtered instead of clipped · `start_count` charged to a rejection · a multi-line crossing entering
a splash pool · the raw introns unsorted · `exact_cut` by `upper_bound` · the junction acceptor ignored ·
`merge_from` dropping `node_start_count`.

⚠ **Two perturbations are green for reasons that are NOT gate holes, and both are worth knowing.**

1. **`L = (end − start) − Σ intron` over the NORMALISED introns is an algebraic identity**, not a bug —
   once the introns are sorted, disjoint and clipped, that sum *is* the complement of the segment total. B1
   is only reachable from the **raw** list, which is what the corrected perturbation used, and that one is
   caught. ⭐ So the reference's fix ("`L` is *defined* as the segment total") is robust rather than merely
   correct: normalising first makes the wrong formula unreachable.
2. **Removing the `L == 1` edge-quantum guard is undetectable on this machine.** ARM64 `udiv` by zero
   yields 0, which is exactly what the guard produces, so the observable answer is identical. On x86 the
   same code raises SIGFPE. **The gate cannot defend that guard here**; it is defended by the comment and by
   review, and this is stated rather than left as an apparent pass.

### What the wiring changed, beyond transcription

**`set_junctions` is a second method, and `scan` refuses to run without it.** A missing junction table is
invisible — every observed intron simply reads as unannotated, so all 404,168 junction edges and both
spliced banks come back empty from a scan that looks well-formed. Pass empty arrays to declare an
annotation with none. `AccumulatorSet::set_junctions` does the per-reference slicing of the flat CSR
(`offsets[c0..c1] − offsets[c0]`, `acceptor_cut − c0`), which is also why a reference's junction ids are
`slot − j0` and the payload's junction axis is the flat slot order concatenated in reference order.
⭐ `native_parity_on_real_data.py` re-derives that slicing **in numpy** and compares, so the two routes to
the same arithmetic are checked against each other rather than one calling the other.

**⛔ A live bug the new API cannot express, found while writing the adapter.** For a fragment with blocks on
more than one reference the shipped code computed a span *per reference* and deposited **all of them** onto
`exons.front().ref_id` — chr7 coordinates landing on chr1's cut axis. `ws.span_ref` recorded which
reference each span belonged to and **nothing ever read it**, which is how it survived. Such fragments now
deposit nothing. ⚠ Deliberately narrow: the test is multi-reference, **not** `cr.chimera_type`, which is
also set for single-reference *cis* chimeras that the intergenic path deposits today — stopping those is a
change to *what counts as a fragment* and belongs in its own arm (the plan's TODO). ⭐ Measured: the two
gates admit **identical** fragments on LBX0190 and MO_3021, so the narrow choice costs nothing today.

**`primary` is deleted, not renamed**, and nothing replaces it. The `align_ok`/`motif_ok` gate is gone too:
the deposit rejects an undefined `align_strand` itself and *counts* it (158 / 615 / 246 on the three
libraries — a population that used to vanish).

**The `memcpy` §3.5 promised is not available**, and S3 took the other side of the trade. The accumulator
stores AoS (a 48 B `Node` interleaving two `uint32` pairs and two `uint64` pairs) and the payload is SoA, so
the copy is a strided transpose and every nanobind view needs explicit strides
(`sizeof(struct)/sizeof(element)`). It is O(partition) once per scan against a deposit that is
O(fragments). ⚠ The binding this replaced passed **no strides at all** and was correct only because the old
`Region` happened to be four contiguous `uint32`.

**The payload keys are the specification's `Tally` field names, character for character.** §3.5 had called
one of them `fragment_length_pools` where the reference says `pool_lengths` — one quantity, two names, in
the two files that must agree, which is the `motif_strand` mistake again. `ref_sj_offsets` was missing
outright and is required: without it a consumer cannot locate a reference's junction rows.

**Three silent-skip branches removed, two of them in code S3 itself wrote.** `node_types` is now required
whenever a reference owns a node; a pool-histogram width mismatch throws in `merge_from` and in
`build_result` instead of quietly not summing. All three would have dropped a whole output with nothing to
notice it by. ⚠ Found in the closing review, not while writing — the habit reproduces itself, so it is worth
grepping for `if (a.size() == b.size())` deliberately.

**The silent pool-disable is gone.** `node_types` is now required whenever a reference owns a node, and
throws otherwise. The shipped code disabled the length pools on an empty array — a whole output going
missing with nothing to notice it by — and the specification has no such state, so a C++ mode the spec
cannot express could only ever be a way to disagree with it.

### ⭐ One word for the implicit splice: `sj_implicit` (owner ruling)

`introns_inferred` / `inferred_intron_fragments` were a **second** word for a concept that already had one
(`SPLICE_IMPLICIT`), which is the mistake `sj_strand` was renamed to fix. `sj_` is already the codebase's
prefix for a splice-junction attribute, so a `FragmentPath` now carries `sj_strand` and `sj_implicit` — two
attributes of one splice under one prefix. *Inferred* is gone from the vocabulary, in the reference, the
spec matrix, both C++ files and the two design documents.

### ⭐ `path_ambiguous` — a new rejection outcome, and the measurement that corrected it twice

Owner ruling (design §9.1): an implicit splice **deposits**, with `sj_strand` inferred from the implying
transcript — but only if the candidates all imply the **same intron set**. Otherwise `L`, both quanta, the
pool bin and the set of crossed lines are all undetermined, and it cannot fall back to an unspliced deposit
either, because `L` involves an intron. So it deposits on nothing, leaves `start_count` untouched, and is
counted on its own denominator until the second-pass side buffer drains it.

Three tests were written first and verified failing, then the reference was perturbed three ways to prove
they bite (check deleted · precedence swapped · the two concepts conflated). **Precedence is part of the
contract:** a fragment that is both strand-undefined and path-ambiguous is filed under the **strand**,
because `dropped_ambiguous_path` sizes the population the second pass can **recover** and a fragment with
no genome strand is not recoverable.

⛔ **And the first measurement of it looked like a bug, because it was one.** The rule deferred **100 % of
implicit fragments on all three real libraries and deposited none.** Cause: every gene carries a synthetic
nascent shadow transcript spanning its locus as one exon, so that row implies no intron in any gap, ever —
counting it as a dissenting candidate makes the test unsatisfiable. It also conflates *which intron did this
molecule splice* with *was this molecule nascent*, which is a **component** the accumulator does not decide.
The candidates are the real transcripts, `~is_synthetic` alone (never `is_nrna`; §3 trap 3).

⚠ **A single-isoform falsification test would have caught it immediately and was not written first.** It is
now (`test_implicit_splice_deposit.py`), with a vacuity guard — the first version of that test reported "no
implicit fragment at all", which is how the bug surfaced.

⭐ **Corrected, and the number the side buffer needs:**

| library | implicit splices | deposited | deferred | deferred as a share of everything accepted |
|---|---|---|---|---|
| LBX0190 | 1,684 | 403 (23.9 %) | 1,281 (76.1 %) | **1.03 %** |
| MO_3021 | 4,728 | 1,427 (30.2 %) | 3,301 (69.8 %) | **0.48 %** |
| LBX0588 | 632 | 159 (25.2 %) | 473 (74.8 %) | **0.04 %** |

Consistent across three libraries: about a quarter to a third of implicit splices have a determined path.
⚠ The deferred share is small as a fraction of everything, so the side buffer is not urgent — but it is
**70–76 % of the implicit population**, so that population is mostly unavailable until it lands.

### ⚠ The suite: 240 failures + 17 errors, ALL from one line

`scan_payload.py:102` reads `cal["n_channels"]`, which S3 deleted. Every affected module runs a scan through
`scan_and_buffer`. ⭐ **Verified rather than assumed:** with payload construction temporarily skipped, every
remaining failure is a *consumer* reading a payload attribute (`fl_pool_mass`, `region_contained`,
`r_total`, `ref_boundary_offsets`, `boundaries`). **No failure comes from the deposit rule.** The first
group is S4, the second S5.

### ⚠ The performance gate is deferred to after S4, by owner ruling

`~357 ns/fragment` is measured through `scripts/design/scan_profile.py` → `scan_and_buffer` → the payload,
which cannot be constructed until S4, so the gate is unmeasurable now. Owner: *"i'm not worried about the
performance gate. we will make it fast eventually."* Recorded here rather than papered over with a shim in
the profiler. ⚠ The deposit is allocation-free per fragment (`DepositScratch` + the adapter's intron
buffer), which was the one performance requirement with a measured motive — 22.8 ns, 18 % of the shipped
deposit, in a per-fragment `std::vector`.

### Files touched (S3)

* `src/rigel/native/calibration/accumulator.{h,cpp}` — `sj_implicit`; `path_ambiguous` +
  `DepositOutcome::kAmbiguousPath` + `outcome_key()`; `AccumulatorSet::set_junctions`; `node_types` required.
* `src/rigel/native/bam_scanner.cpp` — `WorkerState` scratch; `set_regions` renamed parameters and
  `max_length >= 1`; new `set_junctions` + the missing-table guard in `scan`; the deposit adapter; the
  payload; the bindings. **Deleted:** `fragment_genomic_spans`, its binding, `ACCUMULATOR_N_CHANNELS`,
  `accumulator_channel_idx`.
* `src/rigel/native/resolve_context.h` — `collect_implicit_splice_introns` reports `out_ambiguous` over the
  real candidates; the emission and `SPLICE_IMPLICIT` classification are unchanged by construction.
* `src/rigel/native/constants.h` — `RawResolveResult::implicit_ambiguous`.
* `src/rigel/pipeline.py` — `_wire_calibration_regions` installs the junction CSR.
* `src/rigel/native.py` — `Accumulator` comes straight from `_bam_impl`; `fragment_genomic_spans` gone.
* `tests/native/` — `test_accumulator_native_parity.py`, `test_accumulator_worker_determinism.py`,
  `test_implicit_splice_deposit.py` new; 3 tests added to the spec matrix.
* `scripts/design/` — `native_parity_on_real_data.py`, `implicit_splice_census.py` new.
* **Deleted:** `src/rigel/_accumulator.py`, `tests/native/_fragment_spans_reference.py`,
  `tests/native/test_fragment_spans_spec.py`.

### What remains — S4

Rewrite `AccumulatorPayload` against the new keys (they are the `Tally` field names) and delete
`scan_payload.py:143-149`'s backwards-compat path. Then delete the two test-local monkeypatches that stand
in for it (`test_accumulator_worker_determinism.py`, `implicit_splice_census.py`) and re-record the
performance number.

---

## S4 — the payload, typed against the specification's own vocabulary (2026-07-30)

**S4 is closed.** `AccumulatorPayload` is rewritten against the new keys, the scan completes end to end
again, and the two test-local monkeypatches that stood in for the payload are deleted. The red has moved
from *"the payload cannot be built"* to *"consumers read fields that no longer exist"* — which is S5, and
is what the S3 probe predicted it would be.

| gate | result |
|---|---|
| schema | **15 pass** (`test_accumulator_payload.py`), written first and verified red |
| ⭐ perturbations | 6 run, **2 holes found and closed** — see below |
| parity + workers + implicit | **80 pass**, now driven through the real payload rather than a monkeypatch |
| `Σ node_start_count == qc.deposited` | holds on a real scan, **through the payload** |
| suite | **1031 pass**; 266 fail + 15 error, all of them consumers reading old attributes |
| ruff + `ruff format` | clean |
| ⚠ deposit cost | **410 ns/fragment**, fixed **0.348 s** — a regression; see below |

### The schema is checked against `Tally`, not against a list

The payload's field names, the C++ payload keys and `_accumulator_reference.Tally`'s field names are now
the same strings. So `test_accumulator_payload.py` reads its field list off `dataclasses.fields(Tally)`
rather than writing one out: a hand-written list would be exactly the mapping table this naming convention
exists to abolish, and it would be free to drift.

⛔ **Two of six perturbations came back green, and both were real holes in the tests.**

| green perturbation | why it passed | closed by |
|---|---|---|
| **dtype coerced instead of asserted** | the test checked the payload's OUTPUT dtype, and `ascontiguousarray(x, dtype=uint32)` happily narrows an int64 — so a coercing payload still reports the right dtype and passes a check on itself. What that hides is a C++ side that has stopped agreeing with the schema | `test_a_WRONG_dtype_is_REJECTED_rather_than_coerced` — the wrong dtype now has to arrive and be refused |
| **a missing QC denominator** | nothing was feeding an incomplete `qc` block | `test_a_MISSING_qc_denominator_is_REJECTED` |

The other four are caught: offsets trusted rather than re-derived · `qc` left as a plain dict · arrays
copied instead of viewed · the row-divisibility check removed.

### What the payload does beyond reshaping

**`qc` is a typed `ScanQC`, not a dict.** Design §10.3 requires these denominators to be emitted so that
every conservation statement can name what it excluded. A dict answers a typo with a `KeyError` at the
call site, possibly far away and possibly never; a frozen dataclass answers at the boundary. Field names
are the specification's own `Tally.qc` keys, checked against them.

**The per-reference offsets are RE-DERIVED from the cut axis, not trusted.** An offset array of the right
*length* can still be inconsistent, and every consumer slices by these — a per-reference offset that drifts
is the defect class that once dropped 476,719 of 476,732 fragments while every golden test passed.

⭐ **Provenance now covers nodes AND edges: `TranscriptIndex.graph_hash`.** `partition_hash` covers
`nodes.feather` only, and its own docstring says so — it keys a cached *scan*, which sees the cut array.
This payload is **edge-keyed by construction**: its junction axis is meaningless against a different
junction CSR. ⚠ The two genuinely differ — the 2026-07-29 flag fix rewrote every `edges.feather` while
leaving every `nodes.feather` byte-identical, so a nodes-only key would have verified **clean** against a
stale payload. `graph_hash` is `partition_hash` plus the junction CSR (offsets, acceptor cuts, strands),
computed on demand and never stored (§3 trap 25). ⚠ `edge_row` is deliberately excluded: it is a join key
that never crosses the ABI, and hashing it would invalidate a good payload whenever an unrelated edge row
moved.

**`strandedness`, design §5.2's other header field, is deliberately NOT here.** It is derivable from
`strand_observations` (the agreement of read-1 orientation with the splice-motif strand over spliced
fragments) and belongs with the strand model rather than half-built into the payload. Logged in `TODO.md`.

### ⚠ `fl.py` owns its own stale pool count now, and that is not compatibility-keeping

`calibration/fl.py` imported `N_FL_POOLS` from `scan_payload`, which cascaded into **four collection
errors** — worse than test failures, because a collection error stops the whole suite and destroys failure
attribution. The constant describes fl.py's OLD six-pool grid (3 region types × 2 compartments), which the
accumulator no longer emits; re-keying that module onto the five structurally-pure pools is S5's R3. So the
constant moved to the module that still indexes by it, marked stale, rather than being kept alive in a
module that has deleted the concept.

### ⚠ The suite: 266 failures + 15 errors, and every one is a consumer

| root cause | count | step |
|---|---|---|
| `'AccumulatorPayload' object has no attribute 'fl_pool_mass'` | 243 | S5 R3 |
| `AccumulatorPayload.__init__() got an unexpected keyword 'boundaries'` | 28 | S5 — test fixtures building the old payload |
| `region_contained` / `r_total` / `ref_boundary_offsets` / `boundaries` | 10 | S5 R1, R2 |

⭐ Nothing fails at payload *construction* any more, which was the whole of S4. The count rose from 240 to
266 because four modules that previously died at collection now collect and fail individually — more
failures, more information.

### ⚠ PERFORMANCE REGRESSED, and both numbers moved

Now measurable again (it runs through the payload), so re-recorded on the same three real cfRNA libraries:

| | before (S2, shipped accumulator) | now | |
|---|---|---|---|
| per-fragment deposit `c` | ~357 ns (treat as 350–400; 3 points, visible scatter) | **410.0 ns** | ⚠ at or just above the top of the band |
| fixed partition cost | 0.108 s | **0.348 s** | ⛔ **3.2×** |

The per-fragment figure is close enough to the old band's edge that scatter could account for some of it;
**the fixed cost is not** — it is `O(partition)` work that tripled, which is consistent with the payload
now materialising eleven arrays over 1.04 M / 1.04 M / 404 k rows against six before, plus a larger
per-worker accumulator. It is paid once per scan, so it amortises on a deep BAM and dominates a shallow
one. ⚠ Recorded as a regression rather than reported as "roughly unchanged": owner ruling is that speed
comes later, but the number that comes later has to be measured against an honest one now. Logged in
`TODO.md`.

### Files touched (S4)

* `src/rigel/scan_payload.py` — rewritten. `ScanQC`; `N_STRAND_COLUMNS` / `N_FRAGMENT_POOLS` replace
  `N_CHANNELS` / `N_FL_POOLS`; offsets re-derived; dtypes asserted rather than coerced; `graph_hash`.
* `src/rigel/index.py` — `graph_hash`.
* `src/rigel/pipeline.py` — passes `index.graph_hash` into the payload.
* `src/rigel/calibration/fl.py` — owns its stale `N_FL_POOLS`.
* `tests/test_accumulator_payload.py` — rewritten, 15 tests.
* `tests/native/test_accumulator_worker_determinism.py`, `tests/native/test_implicit_splice_deposit.py`,
  `scripts/design/implicit_splice_census.py` — monkeypatches deleted; all three read the payload.

---

## Index rebuild — the artifacts now carry the S1 reach (2026-07-30)

S1 changed `_contiguous_reaches`, the **builder**. The indexes on disk still held the old exonic reach, and
S5 is what makes reach load-bearing — so they were rebuilt before it, not after. ⚠ The hazard was never
that something read a wrong number today; it is that `partition_hash` covers `nodes.feather` only, so a
stale `edges.feather` verifies **clean** (`docs/SESSION_HANDOFF.md` §3 trap 25).

⭐ **There were 7 distinct indexes, not 8.** `rigel_runs/rigel_index` and `refs/rigel_index_v7` are
byte-identical — the same index in two places, despite the name. Both were rebuilt from the same output so
they stay identical. (Two further `rigel_index_v5_old` directories are format 5; nothing loads a v5 index.)

### The verification is what makes a rebuild safe, and it is a re-derivation

`scripts/design/verify_index_rebuild.py`. A correct rebuild has an exactly predictable shape, because S1
changed one function:

* **`nodes.feather` byte-identical** — the partition did not move. ⚠ This doubles as the check that the
  rebuild used the **right source**: a different FASTA or GTF moves the cuts, and nodes differ immediately.
* **`edges.feather` differing in the four `reach_*` columns of CONTIGUOUS rows and nowhere else.** Junction
  reach is deliberately unchanged — a junction edge is only used by a molecule that spliced across it, so
  what remains either side is exonic, and `_junction_edges` stays on the exonic reach.

All 7 rebuilt clean:

| index | contiguous rows whose reach moved | junction rows moved |
|---|---|---|
| human (×2 paths) | **37.5 – 39.2 %** of 1,447,763 | **0** |
| `ambig_dense_10mb` | 26.6 – 41.3 % | 0 |
| the four 5 Mb scenarios | 25.5 – 34.7 % | 0 |
| `quick_3to1_5mb` | 31.1 – 34.7 % | 0 |

⭐ **And ~38 % of human contiguous edges had the wrong reach**, which is the size of what S5 would have been
measured against had this been left.

### ⭐ Both hashes are unchanged, and that is the right answer

`partition_hash` **514f422a85c17163** and `graph_hash` **7856cc59463c903d**, identical before and after.
Not a null result — it is the provenance design being correct: neither key covers `reach`, because neither
the scan nor the accumulator reads it. **A reach change genuinely does not invalidate a cached payload.**
⚠ The corollary is a gap, and it is logged: reach IS consumed by calibration, so anything caching a
*calibration* output needs a key that includes it. Nothing does today.

Byte-identity re-confirmed against the rebuilt index: LBX0190 (60,000 fragments, 44 refs) and the whole of
MO_3021 (875,670, 104 refs). Suite and ruff unchanged from S4.

### ⛔ An index cannot be rebuilt from its own manifest

`manifest.json` records `format_version`, `rigel_version` and `mappability` — **not the FASTA, the GTF, or
the flags**. The human source had to be inferred (by matching node and transcript counts), and
`--collapse-duplicate-transcripts` was discovered from a build failure: the GENCODE GTF has 141
duplicate-transcript groups totalling 283 transcripts, which are mathematically unidentifiable and must be
collapsed. ⚠ The rebuild was only safe because `nodes.feather` byte-identity would have caught a wrong
source. Logged in `TODO.md`: record the source paths, their hashes and the build flags in the manifest.

### Files touched

* `scripts/design/verify_index_rebuild.py` — new.
* Seven index directories on disk (outside the repo): `edges.feather` replaced; everything else, including
  `splice_blacklist.feather`, preserved.

---

## ⛔ Benchmarks and indexes DELETED — the standing baseline is void (2026-07-30)

**Owner decision: every simulated benchmark suite and every built index was deleted**, and both are being
rebuilt from scratch along with a new benchmark suite. This entry exists so that the numbers above are read
correctly from here on, rather than silently reused.

⚠ **This file is append-only, so nothing above has been rewritten.** Every measurement in the earlier
entries stood when it was recorded. What changes is which of them can still be *reproduced*.

### What is VOID — do not quote, compare against, or attempt to reproduce

| | |
|---|---|
| **The standing baseline** `r0 0.079005` / `r3 0.046675` | It was the 32-condition `ambig_dense_10mb` suite, which no longer exists. Every "32/32 flat" gate above refers to it |
| **The goldens** | Same reason. ⚠ Their unexplained 16–52 % `gdna_em_count` move at the index change was **never validated against truth**; the new suite is what should adjudicate it, and that question is still open |
| **Cached payloads** | Built against a deleted index *and* the pre-S3 payload schema |
| **The `calibration-benchmark` skill** | Deleted — it pointed exclusively at deleted paths, and a skill that runs a nonexistent command is worse than an absent one. ⭐ Its methodology was preserved into `docs/testing/BENCHMARKING.md`: net fragment flow, why hard per-fragment labels are the wrong target, and the three caveats |

### What SURVIVES, and why

* **The specification and its gates.** `tests/native/_accumulator_reference.py` is code, and the
  byte-identity result is a statement about **two implementations agreeing**, not about any dataset. It is
  reproducible against any index.
* **Every real-cfRNA measurement** — the fragment-length pools (88.0 / 88.5 / 138.8 / 211.6 / 220.7), the
  implicit-splice census, the multimapper shares, the crossing statistics. These are observations about
  real libraries, and ⚠ **the cfRNA BAMs are not part of this deletion**.
* **The whole test suite.** ⭐ Grep-verified: **no test depended on the deleted data.** The four modules
  that name a suite do so only in comments recording where a number was measured.
* **The index builder and its invariants** (I1–I13), and `scripts/design/verify_index_rebuild.py`.

### ⚠ What becomes unmeasurable until the rebuild

The **S4 performance numbers** (410 ns/fragment, 0.348 s fixed) were measured against the human index and
the real cfRNA BAMs. The BAMs survive; the index does not, so those figures cannot be re-derived until an
index is rebuilt. They are recorded as what they were, not as a live gate.

⚠ **And the census numbers that describe the human graph** — 1,043,881 nodes, 1,043,595 contiguous edges,
404,168 junction edges, median node 151 bp, 15,687 nodes of length 1 — are properties of **the annotation
they were built from**, not constants of the tool. A rebuild from a different GTF moves all of them. They
should be re-derived and re-recorded rather than assumed, and `nodes.feather` byte-identity is the check
that says whether the source is the same one.

### The immediate consequence for sequencing

The open question in `TODO.md` — *build the benchmark suite before S5, or run S5 blind* — is now settled by
circumstance: **there is no suite to run S5 against.** The suite comes first.

---

---

# Moved from `WIP.md` on 2026-08-03 — the index, the suite, S5 and S6

⭐ **Settled work whose conclusions are no longer live inputs.** The cut is at S5.g-2: everything above
it predates the fragment-length critical path. ⛔ Moved verbatim, never rewritten — attribution is the
point of a ledger, and a rewritten entry is a claim about the past rather than a record of it.
## I1 — an index can now be rebuilt from its own manifest (2026-07-30)

`TODO.md` item 1's sub-item, done **before** the rebuild because that is the cheapest moment: the manifest
is written by the build, so recording provenance after the fact would mean recording what I *believed* the
source was — which is exactly the failure being fixed. The previous manifest held `format_version` and
`rigel_version` and nothing else.

⚠ **The manifest was thinner than `TODO.md` said.** That entry lists `mappability` as a recorded key. There
is no such key and never was: per-region mappability was removed in v0.5.0 (`index.py:955`) and the
alignable store now contributes only `splice_blacklist.feather`. Two keys, not three.

### What it records

```json
"sources": { "fasta": {path, bytes, sha256}, "gtf": {...}, "alignable_zarr": null },
"build_flags": { collapse_duplicate_transcripts, feather_compression, gtf_parse_mode,
                 nrna_tolerance, splice_blacklist_min_count, write_tsv }
```

Defaults are written **explicitly**: a rebuilder must not have to know this rigel version's defaults to
reproduce the artifact. `alignable_zarr` is `null` when absent rather than omitted — `null` says "no store
was supplied", a missing key would say "this manifest predates the field".

⚠ **This is not trap 25, and the docstring says so** so that nobody "fixes" it later. Trap 25 forbids
storing a hash of an artifact *beside* that artifact, because the two drift and the stale hash then
verifies clean. These are digests of **inputs the index cannot recompute from itself** — the genome and the
annotation are not in the index. Provenance, not a cache key. `partition_hash` and `graph_hash` are
untouched and still computed on demand.

⭐ `sha256` rather than the `blake2b` the two partition keys use, deliberately: those are internal keys
nobody types, this one is meant to be checked against `shasum -a 256` by whoever is hunting for the source.

### ⭐ The gate had TWO holes, and only perturbation found them

13 tests, written first and verified failing (7 failed / 1 passed — the 1 was the regression guard on the
two pre-existing keys, which is the correct answer). Then the implementation was deliberately broken ten
ways. **Four perturbations were caught by the tests as first written; two were not**, and both would have
been live on the real sources:

| | perturbation | caught by the first 8? |
|---|---|---|
| P1 | digest the **path** instead of the content | ✅ 2 tests |
| P2 | drop one build flag (`nrna_tolerance`) | ✅ |
| P3 | omit `alignable_zarr` when absent instead of writing `null` | ✅ |
| P4 | record the flag **defaults** rather than the values used | ✅ |
| **P6** | **read only the FIRST chunk of the file** | ⛔ **NO — all 8 passed** |
| **P8** | **record the raw input string, not the resolved path** | ⛔ **NO — all 8 passed** |
| P9 | directory digest ignores file **names** | (branch had no test at all) |
| P10 | directory digest walks only the top level | (branch had no test at all) |

**Why they hid:** every fixture is ~2 KB against a 1 MB read chunk, so the chunk loop never iterates; and
every fixture path is already absolute, so `resolve()` is a no-op. On the real inputs P6 would digest one
megabyte of a **1.6 GB** GTF — two annotations differing after byte 1,048,576 would record the *same*
digest, which is precisely the failure this exists to prevent. P8 would record `../refs/genes.gtf`.

Closed by `TestTheDigestItself`: two files with a byte-identical first chunk must digest differently; the
chunk size is monkeypatched to 7 and must not move the answer (**if it does, it is not an I/O buffer**);
a relative path must be recorded absolute; and the directory branch — reachable because an alignable store
may be an unpackaged `.zarr` directory — gets its first tests. All ten perturbations now fail the gate.

⭐ **The flag list is not written out in the test.** It is read off `inspect.signature(TranscriptIndex.build)`
minus the source parameters, so a **new** build knob that never reaches the manifest fails this file rather
than silently escaping provenance. Same device as S4's schema tests reading their field list off `Tally`:
a hand-written list is how the manifest got thin in the first place.

### Gates

* 13 new tests green; **the full suite is 1044 passed / 266 failed / 15 errors** against a baseline of
  **1031 / 266 / 15** re-recorded from this tree in this session. ⭐ The failure and error counts are
  **unchanged** — every one is still an S5 consumer, and I1 added exactly its 13 passes.
* `ruff check src/ tests/ scripts/` clean; `ruff format src/ tests/` clean (⚠ never over `scripts/`).

### Files touched

* `src/rigel/index.py` — `source_record`, `_sha256_of_file`, `_sha256_of_directory`,
  `_DIGEST_CHUNK_BYTES`; the manifest block in `build`.
* `tests/test_index_provenance.py` — new, 13 tests.

---

## I2 — the human index rebuilt from scratch, and it is the SAME index (2026-07-30)

`TODO.md` item 1. Every index on disk was deleted on 2026-07-30; this is the replacement, built from the
sources that survived the deletion:

| | |
|---|---|
| FASTA | `rigel_runs/refs/genome_controls.fasta.bgz` — 286 refs (GRCh38 primary + 169 scaffolds + 92 ERCC) |
| GTF | `rigel_runs/refs/genes_controls.sorted.gtf` — GENCODE v46 / Ensembl 112, 1.6 GB |
| alignable | `~/Downloads/alignable/alignable_grch38_star.zarr.zip`, 26.6 GB |
| flags | `--collapse-duplicate-transcripts --no-tsv`, and they are now **in the manifest** (I1) |
| out | `rigel_runs/refs/rigel_index` |

⭐ **69 seconds, 5.5 GB peak RSS.** The rebuild was assumed to be the expensive part of this session and it
is not — which is why the determinism gate below was affordable.

### ⛔ The named gate could not run, and what replaced it

`scripts/design/verify_index_rebuild.py` diffs a rebuild against **its predecessor**. There is no
predecessor: the deletion removed all seven. The substitute is the one that *is* available from scratch and
is strictly stronger for this purpose — **build twice and diff the two builds**:

* `verify_index_rebuild.py` on the pair: nodes byte-identical, **0 reach rows moved** in all four columns,
  0 junction rows. The script runs unmodified; 0/0 is the from-scratch shape.
* and then every artifact compared byte-for-byte, which the script does not do (it *permits* reach to
  move): `nodes` `edges` `intervals` `sj` `transcripts` `ref_lengths` `splice_blacklist` — **all seven
  byte-identical**, and the two manifests identical. The build is deterministic end to end.

The second build was then deleted. It is a derived artifact and the rules do not keep those.

### ⭐ THE CENSUS RE-DERIVED — and it reproduces the deleted index EXACTLY

`scripts/design/index_census.py`, new. The deletion entry warned that these numbers describe an annotation
rather than the tool and must be re-derived. They were, and **every one matches**:

| quantity | re-derived | previously recorded |
|---|---|---|
| nodes | **1,043,881** | 1,043,881 ✅ |
| contiguous edges | **1,043,595** | 1,043,595 ✅ |
| junction edges | **404,168** | 404,168 ✅ |
| median node length | **151 bp** | 151 bp ✅ |
| nodes of length 1 | **15,687** | 15,687 ✅ |
| seam census (terminus / splice / BOTH / neither) | **40.70 / 58.31 / 0.99 / 0.00 %** | identical ✅ |
| BOTH, absolute | **10,337** | 10,337 ✅ |
| merged regions (v7 partition) | **752,654** | 752,654 ✅ |
| termini the merge hid | **232,451 = 53.42 %** | 232,451 / 53.4 % ✅ |
| transcripts parsed | **254,319** | — |
| single-exon transcripts | **26,475** | 26,475 (trap 3) ✅ |
| mean junction fan-out / max | **1.31 / 25** | 1.31 / 25 ✅ |

⭐ **That agreement is the provenance result, not a formality.** `nodes.feather` byte-identity was the only
thing that would have caught a wrong source last time, and there is no old `nodes.feather` to compare to —
so the census *is* the source check, and it says the surviving GTF is the one the deleted index was built
from. The manifest now records its sha256 so nobody has to infer it again.

⚠ **Two numbers that were quoted with the wrong scope, corrected here.** *Mean node length 2,552 bp* is
**chr1**; genome-wide it is **2,970**. And `docs/SESSION_HANDOFF.md` §1 fact 2's denominator, *"232,451 of
435,291"*, re-derives as 232,451 of **435,107** terminus-flagged contiguous seams — the numerator and the
53.4 % are exact, the denominator is 184 out and depends on how a terminus is counted.

Also new, from the census and load-bearing for the deposit: **591,442 nodes (56.7 %) are shorter than
200 bp**, i.e. shorter than one RNA fragment — `docs/SESSION_HANDOFF.md` §0 C5's unpriced correlation risk applies
to the majority of the partition, not a corner of it.

### Gates, all on the rebuilt index

* **I1–I13 at human scale**: run *inside* the build (`index.py:819` passes `transcripts`, which is what
  gates I3b/I4/I11/I13); a violation raises, and the build exited 0.
* **The manifest's digests verified against the system tool**: `shasum -a 256` on the 798 MB FASTA returns
  `aa1f257f…64a`, character-for-character what the manifest recorded. The interoperability claim in I1's
  docstring is a measurement, not an intention.
* **The S3 gate re-run on real cfRNA**: ✅ byte-identical on LBX0190 (60,000 fragments / 44 refs) and on
  **the whole of MO_3021 — 875,670 fragments across 104 references**, every `Tally` field and every QC
  denominator. The specification and the C++ still agree, now against an index built from scratch.
* **Suite and lint** unchanged from I1: 1044 / 266 / 15, ruff clean.

### ⚠ Performance re-recorded — the S4 figures were unmeasurable and now are not

`scripts/design/scan_profile.py`, three real cfRNA libraries, `OMP_NUM_THREADS=1`, 3 repeats, minimum used:

| | pre-rework | S4 (deleted index) | **now** |
|---|---|---|---|
| per-fragment deposit | ~357 ns | 410.0 ns | **367.8 ns** |
| fixed partition cost | 0.108 s | 0.348 s | **0.315 s** |

⚠ **Read this as a re-measurement, not as an improvement.** It is a 3-point regression on a different
machine state against a rebuilt index, and 367.8 sits inside the 350–400 band the S2 baseline was already
quoted with. What does *not* move is the finding S4 recorded: **the fixed `O(partition)` cost is ~3× the
pre-rework 0.108 s** and that is still unexplained by scatter. `TODO.md` keeps it, with the same reason.

### Files touched

* `scripts/design/index_census.py` — new.
* `~/Downloads/rigel_runs/refs/rigel_index/` — the index itself (outside the repo).

---

## I3 + B0 — the suite backbone, and the gate that says whether it can judge anything (2026-07-30)

`TODO.md` item 2. Two things landed together because the second is what licenses the first.

### I3 — a real human backbone: chr22 + 92 ERCC controls

`scripts/sim/build_suite_reference.py`, new. Carves references out of the real genome+annotation,
preserving the source `.fai` order so the suite's reference-id space is a **subsequence** of the
production one rather than a reshuffle. 93 references, 50,901,224 bp, 71,566 GTF lines. The index builds
in **1.4 s** — which is the whole point of `testing_plan.md`: calibration cost is depth-independent
(`docs/SESSION_HANDOFF.md` §1 fact 22), so a genome-scale index taxes every future iteration at ~66 s a run.

⚠ **The ERCC controls are not filler.** `docs/SESSION_HANDOFF.md` §3 trap 20: a single-reference synthetic index
hid a reference-id-space mismatch that silently dropped **476,719 of 476,732** real fragments inside
`deposit()` while every golden test passed. 92 extra references cost ~83 kb and make that space non-trivial.

| | human | **chr22 + ERCC** | ⛔ deleted `ambig_dense_10mb` |
|---|---|---|---|
| nodes | 1,043,881 | **22,231** | 1,698 |
| median node | 151 bp | **147 bp** | — |
| **nodes / merged regions** | 1.387 x | **1.386 x** | ⛔ **1.000 x** |
| **termini the merge hid** | 53.42 % | **52.72 %** | ⛔ ~0 |
| seam census term/splice/BOTH | 40.70/58.31/0.99 % | **40.82/57.96/1.22 %** | ⛔ 11.8 % BOTH, 12x over |
| 1 bp nodes | 15,687 | **332** | — |

⭐ **The backbone is structurally human at 2 % of the size.** That is the argument for a real chromosome
over a generated mini-genome, and it is now a measurement rather than an assertion.

### B0 — `suite_resolves.py`: prove it can resolve the axis BEFORE quoting a number from it

`docs/SESSION_HANDOFF.md` §3 trap 15 made executable, and **written before the suite exists** so it cannot be
retrofitted to whatever the suite turned out to be.

⭐ **THERE ARE NO TUNED THRESHOLDS IN IT.** A pass mark would be a magic number. Every requirement is
scored against its **degenerate value** — the number a structurally blind suite scores — and those are
boundaries, not choices: a partition that cannot be resolved has `nodes / merged` **exactly 1.000**; a
suite with no length variation has variance **exactly 0**; a Poisson simulator has `omega` **exactly 0**.
A requirement passes iff it is strictly on the non-degenerate side. The human index prints beside it for
scale, never as a threshold.

**Structural half (needs no reads): (d) (e) (g) — all green on the backbone.** The empirical half —
(a) density step, (b) fragment-length variance, (c) non-Poisson counts, (f) low-gDNA x strong-capture —
says `SKIP` until a suite exists. That is the honest state, not a stub.

#### ⭐ Its teeth, proven on TWO degenerate references, each failing for its own reason

| | deleted-suite shape | starved toy |
|---|---|---|
| construction | one single-exon transcript per gene, one strand | + and − transcripts at IDENTICAL coordinates, nothing else |
| (g) nodes/merged | ⛔ **121 nodes vs 121 merged = 1.000 x** | ⛔ 1.000 x |
| (g) hidden termini | ⛔ 0 | ⛔ 0 |
| (d) termini inside an exon | ⛔ 0 | ⛔ 0 |
| (e) single-stranded nodes | OK 60 | ⛔ **0** |
| (e) both-stranded nodes | ⛔ **0** | OK 50 |
| exit code | **1** | **1** |

⭐ The first reproduces the deleted suite's exact defect — **121 == 121**, the same shape as its
`1,698 == 1,698`. The two are complementary: (e) decomposes into "has single-stranded nodes" and "has
both-stranded nodes for them to support", and each degenerate reference fails a different half. A single
pooled (e) check would have passed on both.

### ⛔ AND THE PILOT TURNED UP A STRUCTURAL COUPLING NOBODY HAD RECORDED

Isolated by bisection, same reference, same 20,000 fragments, one condition:

| | capture OFF | capture ON |
|---|---|---|
| wall clock | **1.3 s** | **92.1 s** — ~70 x |

`CaptureSampler._extra_landscape` runs several numpy sorts and `reduceat`s **per (transcript, fragment
length)** pair, and `_mass_cache` is keyed the same way. Distinct fragment lengths drawn, measured:

| | distinct lengths |
|---|---|
| ⛔ deleted suite, `frag_std: 0` | **1** |
| pilot RNA, 206 +/- 98 measured from MO_3021 | **538** |
| pilot gDNA, 157 +/- 125 measured from MO_3021 | **596** |

So the capture path is **O(on-panel transcripts x distinct fragment lengths)** in time *and* in cache
size — 5,326 x 538 = 2.9 M entries. ⭐ **The zero fragment-length variance that made the deleted suite
unable to test anything length-dependent is also what made it fast.** Requirement (b) and the capture
axis are coupled through the simulator, and that coupling is a plausible reason `frag_std: 0` was ever
chosen. Recorded here; the fix (the per-start weight profile is a sliding window over probe coverage and
is derivable for all lengths at once) is its own arm with its own gate.

### Files touched

* `scripts/sim/build_suite_reference.py`, `scripts/sim/design_suite_probes.py`,
  `scripts/design/suite_resolves.py`, `scripts/sim/configs/chr22_pilot.yaml` — all new.
* `~/Downloads/rigel_runs/suite/{reference,rigel_index}` — outside the repo.

---

## B1 — the capture simulator made fast, exactly (2026-07-30)

Owner: *"We need to simulate a full fragment length distribution accurately. We need to optimize the
capture simulation performance. It is abysmal."* Overdispersion deferred by the same instruction.

### The measurement that framed it

Bisected on identical inputs, one condition, 20,000 fragments:

| | capture OFF | capture ON |
|---|---|---|
| wall clock | 1.3 s | **92.1 s** |

`cProfile` on the capture-on run, and it is not ambiguous:

| | |
|---|---|
| `partition_array` | **123.6 s of 127.9 s = 96.6 %** |
| scalar `partition` calls | **1,052,172** |
| `_local_overlap_weights` calls | **4,210,311** |
| peak RSS at 200 k fragments | **18 GB**, still climbing |

The arithmetic was never the cost. `_accumulate_pool` loops over every **distinct fragment length** and
the old body looped over every transcript inside that, doing ~10 numpy operations on ~300-element arrays
each time. Distinct lengths drawn, measured:

| | distinct lengths |
|---|---|
| ⛔ deleted suite, `frag_std: 0` | **1** |
| pilot RNA, 206 +/- 98 (measured from MO_3021) | **538** |
| pilot gDNA, 157 +/- 125 (measured from MO_3021) | **596** |

⭐ **The zero fragment-length variance that made the deleted suite unable to test anything
length-dependent is also what made it fast.** Requirement (b) and the capture axis are coupled through
the simulator, and that is a plausible reason `frag_std: 0` was ever chosen.

### The two changes, both exact

1. **`partition_array` computes every key at once.** Probes are flattened once per space (the layout does
   not depend on fragment length) and ranked within their key, so probes sharing a rank belong to
   different keys and therefore write to disjoint slices of one scratch buffer — which is what lets the
   max be a plain read-max-write instead of `np.maximum.at`. ⚠ It falls back to the old per-key loop when
   any probe **group** holds more than one interval: a probe split across exons must have its pieces
   *summed* before the max across probes, and reordering that into a max would silently change the model.
2. **`_accumulate_pool` computes effective lengths only for nonzero-abundance transcripts.** With
   `frac_expressed: 0.5` that is half the annotation. ⚠ `eff` is still assembled at full length with
   zeros in the dead rows, so `weights`, `probs` and therefore the `rng.choice` draw are unchanged — a
   speed change, not a behaviour change.

⚠ **Nothing is cached in the partition path any more.** `_mass_cache` was 100 % miss-then-store there:
`partition_array` visits each `(key, fragment length)` pair exactly once, so every entry it wrote was
never read. That was the 18 GB. `TestNoUnboundedCache` pins it, and failed with **12,000 stray entries**
before the fix.

### Result

| | before | after |
|---|---|---|
| 20 k fragments, capture on, 1 worker | 127.9 s | **15.4 s** — 8.3 x |
| peak RSS at that size | 18 GB (200 k) | **1.26 GB** |
| **1 M fragments, capture on, 8 workers** | (would not complete) | **24.9 s, 797 MB** |

⭐ **The output is content-identical to the pre-optimisation simulator** — decompressed FASTQ compared
byte for byte across all three implementations. ⚠ The compressed `.fq.gz` files differ, but only in the
gzip mtime header; comparing them directly is the wrong test and briefly read as a regression.

### ⛔ THE GATE HAD NO TEETH AT FIRST, AND ONLY PERTURBATION SHOWED IT

131 tests written first, gated on **brute-force enumeration over every start position** — an explicit
Python loop sharing no helper with the implementation (`docs/SESSION_HANDOFF.md` §3 trap 1). They passed against
the old code, establishing it as the oracle, and the cache test failed, proving the memory claim.

Then the new code was broken nine ways. **Seven of nine perturbations PASSED.** The cause: one layout in
the shared pool was a probe split across two pieces, which makes `_flat_probes` return `None` for the
**entire space** — so every case silently ran the fallback and the batched code was never executed.

Fixed by running the matrix through `partition_array` rather than the scalar `partition`, and by
`assert_batched_path_is_live`, which asserts the path under test is the path being taken. After that:

| | perturbation | caught |
|---|---|---|
| P1 | last rank wins instead of max | 25 failures |
| P2 | sum across probes instead of max | 43 |
| P3 | drop the `min_overlap` gate | 9 |
| P4 | drop the per-probe scale | 10 |
| P5 | support starts at the probe, not one window before it | 76 |
| P6 | no right-hand template clipping | 46 |
| P7 | drop the zero-width-segment mask (the `reduceat` quirk) | 2 |
| P8b | two probes of one key share a rank -> duplicate indices, silent loss | 25 |
| P9 | ignore the key offset, all keys share slice 0 | 6 |

⚠ Two attempted perturbations (rank by probe end; globally-assigned ranks) are **not** defects — rank
assignment only has to be a bijection per key, and max is order-independent. They are recorded so nobody
reads their green as a hole.

### Gates

* 250 tests in `tests/test_sim_capture_partition.py`, all green; nine real perturbations caught.
* Full suite **1294 passed / 266 failed / 15 errors** against 1031/266/15 at session start — failures and
  errors unchanged, every one still an S5 consumer.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.

### Files touched

* `src/rigel/sim/capture/sampler.py` — `_flat_probes`, batched `partition_array`, `_flat_probe_cache`.
* `src/rigel/sim/wgs_engine.py` — `_accumulate_pool` restricted to live abundances.
* `tests/test_sim_capture_partition.py` — new, 250 tests.

---

## B1a — the gDNA arm: root cause found, and it is NOT where three hypotheses said (2026-07-30)

Appended rather than folded into B1, because B1's numbers were measured on the RNA arm and stand. The
gDNA arm did not move with them, and this records why — including the wrong guesses, since each cost a
measurement.

### The RNA arm, re-verified after every change

| | |
|---|---|
| 1 M fragments, capture on, 8 workers | **29.0 s, 1.13 GB** |
| content vs the pre-optimisation simulator | **IDENTICAL** |

And the full 8-condition pilot at 5 M fragments each completed in **1101.7 s (18.4 min)**, 8.3 GB on disk.

### ⛔ Three hypotheses, three measurements, three misses

The gDNA arm at 200 k fragments / 1 worker started at **267.6 s / 37.2 GB** against the RNA arm's
**31.2 s / 2.9 GB** at the same depth.

| | hypothesis | change made | result |
|---|---|---|---|
| 1 | the scratch buffer spans whole chromosomes | buffer compacted to merged probe runs | 268.8 s / 37.2 GB — **no move** |
| 2 | the scalar `partition` bypasses the batched path | `partition` routed through `partition_array` | 268.7 s / 37.2 GB — **no move** |
| 3 | `_flat_probes` rejects the gDNA space (split probe groups) | slots re-scoped from KEY to merged RUN, so split groups batch | 354.4 s / 38.5 GB — **worse** |
| 4 | `_accumulate_gdna_counts` calls `partition` per chromosome | batched to one call per fragment length | 319.7 s / 38.6 GB — **no move** |

⚠ **Every one of those changes is kept**: each is correct, each is content-identical, and 2–4 are what
make the gDNA space eligible for the batched path at all. None of them was the bottleneck.

### ⭐ The actual cause, from the profile

`sample_starts` — **22,838 calls, 231.1 s cumulative**, driving `_extra_landscape` (216.1 s) and
**6,409,160** `_local_overlap_weights` calls, plus 21.0 s of `argsort`.

It rebuilds the capture landscape of a **whole chromosome** for each `(reference, fragment length)` pair,
and `_mass_cache` retains those arrays — that is the 38 GB. The decisive number is the ratio: 200,000
gDNA fragments spread over 93 references x ~596 distinct lengths means each call materialises a
~2.9 M-element landscape **to draw a mean of about 4 fragments.** Enormous setup per tiny draw.

⚠ **This is why `partition` was the top of the FIRST profile and is not the top now.** Both were always
present; batching `partition` simply uncovered the one underneath. A single profile taken before any
change would have mis-ranked the work — `docs/SESSION_HANDOFF.md` §3 trap 16, in a new costume.

### What the fix requires, and why it is not taken here

The landscape exists only to sample a start position from `off_target_weight + binding_per_base * w(s)`.
That draw can be done in two stages — choose on-probe versus off-target using the per-key totals
`partition_array` **already computes**, then place within the chosen region — with no per-start array at
all. Distributionally identical.

⛔ **But it consumes RNG draws in a different order, so the output stops being bit-identical.** Every
other change in B1 and B1a was gated on content identity against the original simulator, and this one
cannot be. That makes it its own arm with its own gate, and an owner decision rather than mine.

### Gates

* `tests/test_sim_capture_partition.py` **262 tests** green (was 250; the split-group case moved from
  asserting a fallback to asserting the batched path is live).
* Full suite **1306 passed / 266 failed / 15 errors** — failures and errors unchanged from 1031/266/15
  at session start, all S5 consumers.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.
* RNA arm content-identical; gDNA arm content-identical through all four changes.

### Files touched

* `src/rigel/sim/capture/sampler.py` — `_flat_probes` now carries `probe_group`; `partition_array`
  slots on `(group, piece)` within a merged run; `partition` delegates to it.
* `src/rigel/sim/wgs_engine.py` — `_accumulate_gdna_counts` batched over references.

---

## B1b — the gDNA arm fixed, and it needed no RNG change after all (2026-07-30)

Owner authorised a sampling rewrite that would reorder the RNG draws. It was not spent: the cost was
structural, not algorithmic, and removing it left the output **content-identical**.

### What B1a's profile actually pointed at

`sample_starts` -> `_extra_landscape`: **22,838 calls, 231 s cumulative, 6,409,160
`_local_overlap_weights` calls, 21 s of `argsort` alone.** It looped in Python over every probe on the
key and lexsorted the concatenation — and in the gDNA space a key is a whole chromosome carrying every
probe, so one call built a 76 k-element landscape **to place a mean of about four fragments**.

### The change

The run/slot computation that `partition_array` already used was extracted into `_run_landscape`, and
`_extra_landscape` now reads out that same buffer instead of having its own implementation. ⭐ Two
consequences beyond speed: the Python probe loop and the lexsort are gone, and **the two consumers can no
longer disagree about what the capture landscape is** — they were separate implementations of one
quantity, which is the shape of `docs/SESSION_HANDOFF.md` §3 trap 27.

⚠ **`_mass_cache` is deleted outright.** It never hit in either path: `partition_array` visits each
`(key, fragment length)` pair once, and `sample_starts` is driven by a counts dict keyed by exactly that
pair. Every entry it ever stored was written and never read. That was the 18 GB (RNA) and 38 GB (gDNA).
Two tests now assert that **no** per-call dict grows across 300 fragment lengths — pinning the property
rather than the one attribute that happened to hold it.

### Result — the gDNA arm, 200 k fragments, 1 worker

| | |
|---|---|
| start of session | 267.6 s / **37.2 GB** |
| after B1a's four changes | 319.7 s / 38.6 GB (none was the bottleneck) |
| **after this** | **161.9 s / 2.72 GB** |

**13.7x less memory**, and R1 and R2 content-identical to the original simulator.

### ⚠ What is still there, measured not guessed

`_run_landscape` is now 47.8 s tottime over 23,788 calls, and `scatter` **511,711 calls / 44.8 s** — the
`(group slot x piece slot)` loop, averaging 21.5 iterations per call because at large fragment lengths
probe neighbourhoods merge into big runs. The remaining lever is the rejection sampler the owner
authorised: propose from the separable per-probe trapezoid (closed form, no per-start array), accept with
`max_g / sum_g` evaluated at the one drawn position. That one **would** reorder the RNG draws.

### Gates

* `tests/test_sim_capture_partition.py` **263 tests** green.
* Full suite **1307 passed / 266 failed / 15 errors** — failures and errors unchanged from 1031/266/15 at
  session start, all S5 consumers.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.
* Content identity verified on both arms against the pre-optimisation simulator.

### Files touched

* `src/rigel/sim/capture/sampler.py` — `_run_landscape` extracted and shared; `_extra_landscape` rebuilt
  on it; `_mass_cache` removed.
* `tests/test_sim_capture_partition.py` — cache tests generalised to "no per-call state accumulates".

---

## B2 — the pilot suite exists, and B0 says what it can and cannot judge (2026-07-30)

### The pilot

The owner's 8-condition minimum: `gDNA {none, 100 %} x strand {0.50, 0.99} x capture {off, on}`, a clean
2^3, nascent held at `none`. 5,000,000 RNA fragments per condition on the chr22 + 92 ERCC backbone,
oracle BAM per condition. `scripts/sim/configs/chr22_pilot.yaml`.

⚠ **Every fragment-length parameter is measured, not chosen** — the count-weighted mean and sd of
MO_3021's own pools (`rigel/fragment_lengths.feather`, the library that is also the S3 byte-identity
gate): RNA 206.1 +/- 98.3, gDNA 156.5 +/- 124.6.

| | before the capture work | **after** |
|---|---|---|
| 8 conditions x 5 M, 8 workers | 1101.7 s / **43.1 GB** | **741.4 s / 4.09 GB** |

⭐ **10.5x less memory on the actual deliverable**, and 8.3 GB on disk.

### B0's verdict on it — 6 structural requirements resolve, 2 empirical ones do not

| | requirement | measured | blind | human |
|---|---|---|---|---|
| (g) | nodes / merged regions | **1.386 x** | 1.000 | 1.387 |
| (g) | termini the merge hid | **4,906** (52.72 %) | 0 | 232,451 (53.42 %) |
| (d) | termini strictly inside an exon | **4,958** (53.28 %) | 0 | 234,400 |
| (e) | single-stranded nodes | **18,127** (81.5 %) | 0 | 839,200 |
| (e) | both-stranded nodes | **3,313** | 0 | 171,600 |
| (e) | single-stranded per both-stranded | **5.47** | 0 | 4.89 |
| (a) | capture density step, p90/p10 | **2,123 x** | 1.0 | — |
| (b) | fragment-length sd, worst pool | **87.4 bp** | 0 | — |
| (c) | replicate pairs for omega | ⛔ **0** | 0 | — |
| (f) | conditions at 0 < gDNA <= 10 % with capture | ⛔ **0** | 0 | — |

⭐ **The suite is structurally human at 2 % of the size** — every structural ratio tracks the full index.
⛔ **And it is honest about the two it cannot judge.** (c) is the deferred overdispersion item: with one
draw per parameter set there is nothing to estimate dispersion against, so B0 reports *estimability*
rather than a misleading number. (f) needs gDNA in the 1-10 % band, which the owner's 8-condition
pilot deliberately does not have. **Neither is a defect in the suite; both are named work.**

### ⭐ The empirical half's teeth, on the deleted suite's own shape

Truth files were hand-written to the `ambig_dense_10mb` configuration — every fragment exactly 200 bp,
capture off, one condition, no replicates — and fed to B0. **All four empirical requirements fail**, each
for its own reason: (a) no capture-on condition, (b) sd exactly 0.0, (c) 0 replicate pairs, (f) 0
low-gDNA conditions. Exit code 1. Together with the two degenerate references from I3/B0, every check in
the file has now been shown to fire.

### ⛔ A simulator defect found while building the falsification case

`sampling.truncated_normal_frag_lengths` (`sampling.py:30`) rejection-samples in a `while` loop **with no
termination guard**. An impossible truncation — `Normal(206, 0)` into `[200, 200]`, which is exactly what
"reproduce the deleted suite's frag_std = 0" produces — never draws a valid value and the simulator hangs
forever instead of raising. Logged in `TODO.md`; it cost a 10-minute timeout here.

### Files touched

* `scripts/design/suite_resolves.py` — empirical half added: (a) (b) (c) (f), all read from the suite's
  own truth files rather than a BAM scan, which matters because calibration is red until S5.
* `scripts/sim/configs/chr22_pilot.yaml` — the pilot.
* `~/Downloads/rigel_runs/suite/pilot/` — 8 conditions, 8.3 GB, outside the repo.

---

## B3 — the scan cache: scan once, calibrate many times (2026-07-30)

`TODO.md` item 2's cached substrate, and `docs/testing/testing_plan.md`'s objective 1. Owner's four steps:
determine what calibration needs · cache it · load and resume · seed a toy from it. **Steps 1–3 landed;
step 4 is stubbed with a strict-xfail that names its blocker.**

### Step 1 — what calibration needs, read off the call site

`calibrate()` takes eight arguments from three different places, and only one is expensive:

| input | origin | cached? |
|---|---|---|
| `payload` (21 fields) | the scan | ⭐ **yes** |
| `strand_model` | the scan | ⭐ **yes**, including the per-junction table |
| FL histograms | the scan | ⭐ **yes, RAW** |
| `region_arrays` | `RegionArrays.from_index` — **0.11 s** | no |
| `boundary_flags` | `build_boundary_flags_array` — **0.04 s** | no |
| `gdna_fl_pmf` / `rna_fl_pmf` | `build_fl_models(...)` | no — derived |
| `config` | the thing being varied | no |
| `injected_priors` | fitted BY `calibrate` | step 4 |

⚠ **Anything derivable from the index is rebuilt, never stored** — 0.15 s together against an 8.45 s index
load that happens anyway, and a stored copy is how a cache goes stale against what it describes (trap 25).
⭐ **The FL histograms are cached RAW, not as pmfs**, so `build_fl_models` stays the single source of
truth; freezing its output would mean a change to the FL model silently misses every cached scan.

### ⭐ Step 2 — the key needed THREE parts, and the middle one was a logged gap coming due

* `graph_hash` — the payload already carries it, so the tally self-keys.
* ⭐ **a REACH digest.** `TODO.md` logged that `reach` is consumed by calibration and covered by
  **neither** `partition_hash` **nor** `graph_hash`, and that the gap "becomes live the moment something
  caches a calibration output". **A cache loaded against an index is that moment.** The 2026-07-30 rebuild
  moved ~38 % of human contiguous reaches with *both* existing hashes byte-identical — a reach-blind key
  would have verified clean against a moved index.
* **the scan config**, plus a manifest self-consistency check (the recorded config must re-digest to the
  recorded digest, so a tampered manifest fails loudly).

⚠ **Not pickle.** The old `scripts/debug/calib_cache.py` pickled `calibrate()`'s kwargs; a pickle of
numpy-holding dataclasses is fragile exactly across schema changes — and the schema changes at S5, so a
pickle written today would not load after it. Arrays to `.npz`, scalars and provenance to JSON.

### Result, on the real pilot

`scripts/design/build_scan_cache.py`, all 8 conditions against the chr22 index:

| | |
|---|---|
| scan all 8 once | **45.7 s** |
| reload all 8 | **0.16 s** — ⭐ **291x** |
| per condition on disk | **0.3 – 0.9 MB** |

The writer reads each cache straight back before reporting success: a cache that cannot be loaded against
the index it was just written from should fail then, not three days later.

### ⛔ TWO S5 BLOCKERS FOUND, and neither was in the plan

1. **`gdna_fl_mass` still reads `payload.fl_pool_mass`** — the field S4 replaced with the five pure pools
   (`pool_lengths`). So `calibration_inputs` cannot build the FL pmfs today. Split out
   `index_derived_inputs`, which works now and is tested; `calibration_inputs` is strict-xfail. This is
   plan §4's **R3**, and it is now pinned by a test rather than a list entry.
2. **Step 4's prior extraction needs `calibrate` to run.** `InjectedCalibrationPriors` already exists and
   already carries exactly the population quantities a toy cannot fit — `rna_sense_frac`, `n_rna_obs`,
   `n_gdna_obs`, both strand overdispersions, the enrichment NPMLE, the intron background, the aggregate
   background — and `calibrate` already stashes the fitted bundle in `_debug["calibration_priors"]`.
   ⭐ **The seed mechanism is not new work; `scripts/debug/toy_inject.py` is prior art** — it just points
   at the deleted `ambig_dense_10mb` and its `_selfsolve_cache`, so it is dead code against dead paths.

### ⛔ THE GATE HAD TWO HOLES, AND ONLY PERTURBATION FOUND THEM

14 tests written first and verified failing. Then eight deliberate breakages; **two passed**:

| | perturbation | first version | now |
|---|---|---|---|
| P1 | reach digest covers only the POS columns | ⛔ **passed** | 2 failures |
| P2 | drop the per-junction strand table | caught | caught |
| P3 | cache the derived pmf instead of the raw histogram | caught | caught |
| P4 | also store `region_arrays` in the cache dir | caught | caught |
| P5 | skip the reach check on load | caught | caught |
| P6 | skip the **index** graph_hash check | ⛔ **passed** | caught |
| P7 | drop one payload array on write | caught (7) | caught |
| P8 | lose the dtypes | caught | caught |
| P9 | skip the manifest self-consistency check | — | caught |

**Why they hid.** P1: the digest test perturbed only `reach_lo_pos`, so half the reach was unguarded —
now parametrised over all four columns. P6: the graph_hash test tampered with the *manifest*, which trips
the self-consistency check **before** the index comparison is ever reached — now there is a second test
that moves the manifest and the payload together, producing a cache that is internally consistent but
describes a different index. ⚠ One further perturbation (dropping the column name from the digest) is
**not** a defect — the columns are a fixed tuple of equal length, so the separator cannot change a
collision — and is recorded so its green is not read as a hole.

### Gates

* `tests/test_scan_cache.py`: **14 passed, 2 xfailed** (both strict, both naming an S5 dependency).
* Full suite **1321 passed / 266 failed / 2 xfailed / 15 errors** — failures and errors unchanged from
  1031/266/15 at session start.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.

### Files touched

* `src/rigel/scan_cache.py` — new.
* `scripts/design/build_scan_cache.py` — new.
* `tests/test_scan_cache.py` — new, 16 tests.

---

## Cleanup — docs and scripts, and one coupling the sweep nearly broke (2026-07-30)

Owner request at session end. Recorded because a deletion of 160 files deserves its reasoning written
down, and because it turned up a dependency that a grep had missed.

### The rule applied

**A script stays if it can run today and answers a question we still have.** Everything else goes to git
history, which is the archive — the standing rule is *converge and delete*, not *keep for comparison*.

| | before | after | why |
|---|---|---|---|
| `scripts/debug/` | **113 files, 18,192 lines** | **gone** | 76 referenced deleted suites/indexes or deleted payload fields and could not run at all; the rest were calibration-internals sandboxes that S5 invalidates |
| `scratchpad/` | 31 tracked scripts, 8.9 MB | gone | one-off prototypes for a design that is now implemented and tested. Their *results* are `docs/SESSION_HANDOFF.md` §1; the citations now point at git history |
| `scratch/` | 1 file | gone | a `PROVENANCE.txt` under an empty tree |
| `archive/calibration_legacy_2026_05/` | 3 files | gone | ⛔ legacy code kept for comparison, which the standing rules forbid outright |
| `scripts/sim/` | 11 scripts, 12 configs | 8 scripts, 3 configs | dropped the generator for the deleted suite, a deprecated shim, and 9 configs pointing at deleted outdirs |
| `scripts/calibration/analyze_srd_v2.py` | — | gone | "SRD v2 Phase 4" is a retired calibration generation |

### ⛔ AND THE SWEEP NEARLY DELETED A LIVE DEPENDENCY

`scripts/debug/oracle.py` is **the calibration oracle** — "the oracle IS the production accumulator,
partitioned by true fragment origin", the thing that proves a partition sums to the full payload. It was
in `debug/` and `tests/calibration/test_oracle.py` imports it via a **`sys.path` hack**, which my
coupling check (`grep 'from scripts|import scripts'`) did not match.

Caught by the suite refusing to collect. Restored and **moved to `tests/calibration/_oracle.py`**, which
is where it belonged: it is test infrastructure, and the repo already has that exact pattern in
`tests/native/_accumulator_reference.py`.

⚠ **The first fix for the import was also wrong and the suite said so.** Adding
`tests/calibration/__init__.py` to allow `from ._oracle import …` turned the directory into a package and
broke **five** sibling modules that rely on pytest's implicit `sys.path` insertion (`import _synthetic`).
Reverted; the bare `from _oracle import …` matches the directory's own convention.

⭐ **Lesson: a grep for import statements is not a coupling check.** `sys.path` manipulation, plugin
registries and string-keyed lookups all evade it. Running the suite is the coupling check.

### Documentation

* **`docs/WIP.md` split at 1,715 lines.** Entries up to the deletion moved **verbatim** to
  `docs/WIP_ARCHIVE.md` — relocated, never rewritten, so the append-only property that makes attribution
  meaningful survives. `docs/WIP.md` keeps an index of *every* entry plus the full text from `I1` on.
* **`docs/SESSION_HANDOFF.md` §4 deleted** (441 → 310 lines). It was a pre-settlement proposal list whose six
  worst bullets already carried `⛔ SUPERSEDED` markers; a proposal list that contradicts the shipped
  deposit rule is a hazard, not history. ⭐ **The rest of the file is NOT dead** — §3's traps caught real
  defects repeatedly this session and §2 is the derivation reference; the header now says so and gives a
  read order.
* **`docs/testing/BENCHMARK_SUITE.md` is new** — what the suite is, how to build it, and what it can and cannot judge.
* **`TODO.md` rewritten** against the actual state; **`docs/accumulator/IMPLEMENTATION_PLAN.md` §0** replaced with the
  session handoff; **`CLAUDE.md`** doc list, test counts and tooling table corrected.

### Gates

* Full suite **1321 passed / 266 failed / 2 xfailed / 15 errors** — **identical to before the cleanup**.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.
* ⚠ Nothing committed. `git checkout` recovers any file the owner disagrees with.


---

## S5.0 — S5 stopped, and the observable question answered (2026-07-30)

**Owner ruling.** S5 was written as "rewire the four consumers". Asked which observables the rewiring
should consume, the design had no answer — so wiring it would have frozen a schema nobody had checked.
S5 became a derivation first. `docs/calibration/S5_DESIGN_LOG.md` is the running log;
`docs/accumulator/NODE_DENSITY_DERIVATION.md` is the write-up.

**Two tools built, both re-runnable.**

| | |
|---|---|
| `scripts/design/observable_efficiency.py` | the full length histogram is the sufficient statistic and its Fisher information is closed form, so `efficiency = Var(full) / Var(stored set)` is exactly "what fraction of the information does this storage choice keep?". Swept over 756 gDNA×RNA fragment-length pairs (means 50–300 bp × cv 0.15–1.00, **both directions**) × 3 shape families × 4 compositions |
| `scripts/design/node_density_derivation.py` | T0–T6, each perturbed. Opportunity formulas verified by **exact enumeration** of integer starts |

**What was found.**

* ⭐ **The shipped `(count, Σ1/L)` pair has an exact blind spot.** At an edge the count row is
  `(μ_g−1, μ_r−1)` and the reciprocal row is `(1,1)`, so the determinant is `μ_g − μ_r`: when the two
  components share a mean fragment length the pair carries **literally zero** information about the
  split, at any depth, in every shape family, at every composition. Efficiency minimum **0.000**.
* ⭐ **The deposit weight is the reciprocal OPPORTUNITY, not the reciprocal length** — `h = 1/A(w)`, and
  it is unique (0 violations against 266 for `1/L`). The edge rule is its `A = w−1` case; an **edge is
  the `ℓ→0` limit of a node**, verified identically. `contained + spanning` at reciprocal opportunity
  gives `ρ·(1 − f(ℓ+1))` to 4e-16.
* ⛔ **A model-free channel provably carries no composition information** — model-free means both
  components share a coefficient, so its row is `(K,K)`. The level and the split are two jobs.
  Measured: `Σ1/A` alone scores 0.113 at a 151 bp node, 0.000 at 1000 bp.
* **`Σ L` is the tilt worth adding**: minimum over the grid `0.000 → 0.078` at an edge, `0.078 → 0.188`
  at a 151 bp node; median at a 3000 bp node `0.182 → 0.692`.

⚠ **Two results were wrong before they were right, and both were caught by discipline, not luck.** An
assumed gamma fragment-length shape reversed the ranking (a measured histogram with the same moments
disagreed) — hence three shape families. And the first Monte Carlo came out **2× better than the theory
it existed to test**, because it clipped the estimate and used the wrong solver. The harness now asserts
**monotonicity** (a superset of channels can never carry less information), which is what caught a
covariance-conditioning failure reporting `0.000` where a subset scored `0.832`.

⚠ **`docs/SESSION_HANDOFF.md` §1 fact 5 does not reproduce** and needs correcting or retiring. It was measured
at a 4× mean separation; over the full grid it is false at every node ≥ 250 bp and at the edge.

**Rulings taken** (owner, all recorded in `docs/calibration/S5_DESIGN_LOG.md` §4): every population stores
`(count, inv_length_sum, length_sum)`, nodes **and** edges; the node deposit weight is unchanged; `Σ1/A`
is not stored; `CalibrationResult`'s left/right split becomes per-edge.

## S5.a — `length_sum` on every population, and `density` renamed (2026-07-30)

Two sub-steps, deliberately separated so each has its own gate.

**S5.a-1, the rename — behaviour-neutral, and proven so.** `density` → `inv_length_sum`,
`DENSITY_SCALE` → `INV_LENGTH_SCALE`, `density_quantum` → `inv_length_quantum`, across 9 files.
⭐ `Σ1/L` **is** an exact density at an edge and is **not** one at a node — the design's own §6 says so —
so the old name put one word on two concepts (`docs/SESSION_HANDOFF.md` §3 trap 27). Gate: the full suite
returned **1321 passed / 266 failed / 2 xfailed / 15 errors**, byte-identical to the pre-rename baseline.

**S5.a-2, the new channel.** Five arrays — `node_contained_length_sum`, `node_spanning_length_sum`,
`edge_unspliced_length_sum`, `edge_spliced_length_sum`, `sj_length_sum` — through the reference, the C++
structs, the deposit, the per-worker merge, `build_result`, the zero-copy bindings and the payload.
The deposit change is **purely additive**: one `+= L` per object touched, no quantum changed.

| gate | result |
|---|---|
| 6 spec tests, written first | ✅ **verified failing first**; one failure was a bug in the TEST (wrong edge id), fixed there, not in the code |
| 5 reference perturbations | ✅ all caught — `L−1`, skipped spanning, genomic span instead of `L`, a `L/n` share, skipped contained |
| 2 C++ perturbations | ✅ caught by the right gate — a wrong deposit fails byte-identity, a dropped merge term fails only the shard test |
| real cfRNA, human scale | ✅ **byte-identical on 60,000 fragments over 44 references**, all 5 new fields, all 9 QC denominators |
| ⭐ the gate refused to skip | the real-data script **failed loudly** on an unregistered field rather than silently omitting it from the comparison |
| full suite | **1327 passed / 266 failed** — failures unchanged at exactly the baseline 266 |
| deposit cost | **339.7 ns/fragment** (S4: 367.8), fixed **0.356 s** (S4: 0.315). Regressed over 3 real cfRNA libraries |

⚠ **Per-worker memory grew as designed**: `Node` 48 → 80 B, `ContiguousEdge` 48 → 80 B, `JunctionEdge`
24 → 40 B ⇒ **187 MB/worker** at human scale, against 114 MB before and a measured 8.6 GB peak RSS. The
fixed `O(partition)` cost rose with it (0.315 → 0.356 s), which is the zeroing and the AoS→SoA transpose;
the per-fragment cost did not.

**S5.a-3, the scan-cache schema key.** `payload_schema_digest()` — the accumulator's own field list.
⛔ None of the three existing keys moves when the accumulator changes, so a cache written before this step
would have been **accepted** and then died inside `_payload_from_parts` with a bare `KeyError`. 3 tests,
3 perturbations, all caught — including "a missing key must not read as agreement", which is
`docs/SESSION_HANDOFF.md` §3 trap 25 exactly.

⚠ **Every existing scan cache under `~/Downloads/rigel_runs/suite/pilot/scan_cache` is now correctly
refused and must be re-scanned** (45.7 s for all 8 conditions).


## S5.b — `fl.py` re-keyed to the five pure pools (2026-07-30)

`gdna_fl_mass` read `payload.fl_pool_mass`, a field S4 deleted. Re-keyed onto `payload.pool_lengths`,
and the pool axis moved to `scan_payload.py` beside `N_FRAGMENT_POOLS` — it is the accumulator's own
enum in a third language, and a consumer holding a private copy is how three files come to disagree
about which row is which.

**Three accessors, and the partition between them is the point.**

| | pools | fitted? |
|---|---|---|
| `gdna_fl_mass` | `DNA_INTERGENIC + DNA_INTRONIC` — **contained only** | yes |
| `rna_fl_mass` | `RNA_SPLICED` — annotated junction, splice OBSERVED | yes |
| `splash_fl_mass` | the two `*_EXON` crossing pools | **no — QC only** |

⭐ **Both component models now come from the PAYLOAD, not the scanner's category models.** The scanner's
spliced histogram is transcript-space and requires a UNIQUE transcript; the accumulator's pool is a
structural rule over a larger population, is binned at the same `L` as everything else, and already
excludes `sj_implicit` fragments — whose splice was never sequenced, so certifying them as RNA would make
the pool depend on the model it is used to fit. One quantity, one source. Only the unconditional
`global` histogram still comes from the scan, because no pool is unconditional.

| gate | result |
|---|---|
| 6 tests, written first | ✅ verified failing (ImportError, then per-assertion) |
| 4 perturbations | ✅ all caught — splash folded into gDNA, RNA leaking into gDNA, a swapped index in the SCHEMA, a pool reachable twice |
| the three-way pool-axis contract | ✅ pinned against the executable specification's own `FragmentPool` enum, not a written-out list |
| ⭐ the strict xfail | ✅ **XPASS(strict)** — `calibration_inputs` works; marker removed |
| full suite | **1335 passed / 266 failed / 1 xfailed / 15 errors.** Failures unchanged at the baseline 266 |

**⭐ Real cfRNA, LBX0190 — the first model-free fragment-length measurement this project has, and it
reproduces `docs/accumulator/DESIGN.md` §8's table.** (A *test*, not a design input.)

| pool | n | mean L | design §8 |
|---|---|---|---|
| `DNA_INTERGENIC` | 1,594 | **79.9** | 88.0 |
| `DNA_INTRONIC` | 2,873 | **79.1** | 88.5 |
| `DNA_INTRON_EXON` | 2,092 | 144.1 | 138.8 |
| `DNA_INTERGENIC_EXON` | 123 | 231.1 | 211.6 |
| `RNA_SPLICED` | 63,566 | **227.3** | 220.7 |

The two pure gDNA pools agree with each other (79.9 / 79.1) with no fitted model anywhere, the splash
pools land *between* the two extremes as §8 predicts of on-target gDNA plus leaked RNA, and the
separation is **2.86×** — the determinant of §6's 2×2. ⭐ Folding the splash pools in would move the gDNA
mean **79.4 → 102.4**, confirming on real data the bias §8 attributes to the shipped four-pool axis.

⚠ **NEW FINDING, logged to `TODO.md`: the EB shrink is now doing much more work than it was tuned for.**
`POOL_EB_PRIOR_ESS = 1000` against a pure gDNA pool of only **4,467** fragments pulls the fitted mean
**79.4 → 100.3**, a **+26 %** shift toward the global FL. The old four-pool aggregate was several times
larger, so the same constant shrank far less. The constant has not been changed — it is a magic number
and its value is the owner's call.


## S5.c — `effective_length.py` reduced to the one placements formula (2026-07-30)

Two frames, one family, and nothing else::

    contained   E_f[ (node_len − w + 1)+ ]
    crossing    E_f[ max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)) ]

⭐ **The crossing formula covers both edge kinds and both components.** Mean fragment length is its
large-reach limit, not a separate case: gDNA's template is the chromosome, so `UNBOUNDED_REACH` on both
sides collapses it to `mu − 1` exactly. Each of the four `min` terms is one binding constraint and none
is droppable.

| gate | result |
|---|---|
| **89 tests**, written first | ✅ verified failing (ImportError, then per-assertion) |
| every formula **enumerated** | ✅ 49 contained × 30 crossing cells against integer start positions — *not* a restatement of the closed form. The previous file computed its "brute force" from the same algebra as the implementation, so the pair agreed while both were off by one |
| ⭐ independent cross-check | ✅ reproduces `docs/SESSION_HANDOFF.md` §2's published taper table: **160.095 against 160.1** at R=200, and 88.2/19.8/198.5 against 87.8/19.6/199.0 |
| ⭐ a documentation defect found | that table **mixes two conventions**: its first four entries are symmetric (both reaches = R) and the last, captioned "at a first exon", is one-sided. Recomputing it is what revealed this; the test now pins both readings |

**Deleted, not ported:** `boundary_side_eff_length` (`E[min(l,R)]/2`), `spliced_side_eff_length`
(`E[min²/2l]`), `boundary_side_crossing_count_eff_length`, `boundary_eff_length`, `region_eff_length`.
The first three divided a per-FACE mass and a contiguous edge has no faces; the last two have direct
successors (`contained_eff_length`, `fl_mean`).

⚠ **THE CLEAN CUT COST 23 TESTS, AND THIS IS THE ONE NUMBER TO LOOK AT.**
Suite: **1379 passed / 289 failed / 1 xfailed / 15 errors**, against 266 failures at session start.

| | |
|---|---|
| +23 failures | `test_bp_solver` (16), `test_node_init` (6), `test_node_chain` (1) |
| what they were | tests of `build_node_geometry`'s **per-face** model, driven by hand-built substrates — so they passed while that model existed. It is deleted; they are S5.e's to rewrite |
| ⭐ how they fail | **20 of the 23 raise a NotImplementedError naming S5.e or S5.f explicitly.** A failing test that says which step owns it is worth more than one that passes against a deleted model |

⛔ **An import-time break was rejected in favour of a runtime one, deliberately.** Removing the deleted
symbols outright made `rigel.calibration` unimportable and turned 266 attributable runtime failures into
**36 collection errors — the whole suite stopped running**, destroying the signal every other step is
gated on. `build_node_geometry` and `calibrate()` therefore raise with a message naming their step, and
`test_bp_solver` binds a shim for the two divisors with no successor (never reached: the tests hit
`build_node_geometry` first). All three are deleted by S5.e/S5.f.

⛔ **`tests/calibration/test_message_frames.py` DELETED** (5 tests). Every one asserted a *relationship
between* the per-face divisors — the half-triangle against the two-sided face. That model does not
exist, so the tests could not be ported, only rewritten from scratch against a different geometry. Git
retains them.

⚠ **A7 is still open and S5.c is written so it stays open:** reach constrains an *unspliced* edge
**per component** — the RNA part is bounded by its transcript, the gDNA part is not — so an edge has no
single reach. Production has ignored reach entirely; "keep ignoring it" is expressible as
`UNBOUNDED_REACH` rather than as a second code path, so the decision is not pre-empted.


## S5.d — one substrate, and a chain with no terminals (2026-07-30)

**`substrate.py` — two classes became one.** `CalibrationSubstrate` held three per-region views
(contained / left / right); `BoundarySubstrate` held the same numbers re-keyed by boundary. Two classes,
one set of numbers, two keyings — and they existed **solely because a boundary had two sides**. An edge
is a 0-bp line with one set of numbers, so the second class, the left/right axis, the re-keying identity
and `_make_view` all dissolved together.

What replaces them: one `PopulationView` (`count`, `inv_length_sum`, `length_sum`, each `[n, 2]`) held
five times — node contained / node spanning / edge unspliced / edge spliced / junction — on three axes.

⭐ **The fixed point is decoded HERE and only here.** `inv_length_sum` arrives as
`round(2^32 / placements)` and leaves as a real number, so no downstream module needs to know the scale.
`length_sum` is a plain integer sum and is deliberately **not** scaled.
⚠ `mean_length` is **NaN at zero count, not 0** — zero would read as "the fragments here are infinitely
short" (`docs/SESSION_HANDOFF.md` §3 trap 23).

**`node_chain.py` — `B R B R … R B` became `N E N E … E N`.** A reference with `k` nodes owns `k − 1`
interior lines, not `k + 1` boundary slots, so the chain **starts and ends with a NODE and there are no
terminal slots**. `2k − 1` slots, not `2k + 1`. An **edge therefore always has a node on both sides** —
an invariant the old shape could not state, because its terminal boundaries had one flank.

| gate | result |
|---|---|
| **24 tests**, written first | ✅ verified failing |
| **7 perturbations** | ✅ all caught — undecoded fixed point, `length_sum` scaled too, `mean_length` floored to 0, the per-reference offset check dropped, two populations swapped, the **old `k+1` terminal shape** (9 tests), references bleeding into each other |
| the alignment check | ✅ kept, and it checks per-reference offsets as well as the total — matching totals are not sufficient evidence, that is the defect that dropped 476,719 of 476,732 fragments with every golden green |
| ⚠ two fixture bugs, caught by the tests | two banks shared a total (so "no population aliases another" was untestable), and the fixture floored where the spec rounds half-away-from-zero |

Suite: **1384 passed / 291 failed / 1 xfailed / 15 errors** (S5.c left it at 1379 / 289).

**Renamed, not aliased: `REGION` → `NODE`, `BOUNDARY` → `EDGE`, `ref_idx` → `obj_idx`.** These are the
same concepts under the new vocabulary — a chain slot is an interval or the interface between two — so
the S5.e modules were renamed mechanically rather than given a compatibility alias. An alias would have
been exactly the two-names-one-concept defect the naming rules exist to prevent.

⛔ **Two test modules DELETED, both testing classes that no longer exist.** Git retains them.

* `test_boundary_substrate.py` (11 tests) — every one asserts the left/right re-keying identity of
  `BoundarySubstrate`.
* (S5.c already deleted `test_message_frames.py`, 5 tests, for the same reason.)

`test_density_model.py` binds a `SubstrateView` shim so it still imports; `node_gdna_density` is an S5.e
consumer and its tests now fail naming that step. Deleted by S5.e.


## S5.e — the faces dissolve (2026-07-30)

**A7 was ruled first, before any code**, so it could not be retrofitted to whatever the implementation
turned out to be (`docs/calibration/S5_DESIGN_LOG.md` §1 A7, §4). ⭐ The decision reduces to **one call site**, and that
is what made it rulable: gDNA at any edge is `taper_g = 1` by physics and the two contained frames take
no reach argument, so the only question was the RNA component at a *contiguous* edge.

| edge kind | component | reach | why |
|---|---|---|---|
| contiguous | gDNA | `UNBOUNDED_REACH` | its template is the chromosome — physics, not a ruling |
| **contiguous** | **RNA** | ⭐ **`UNBOUNDED_REACH`** | keeps S5.e varying ONE thing; A7 proper is **S5.g**, the only ordering where it gets an end-to-end A/B against S5.f's first baseline |
| junction | RNA | ⭐ **the real exonic per-strand reach** | a BRAND-NEW population — the predecessor had no junction divisor at all — so wiring it regresses nothing, and it keeps the one code path exercised rather than a dead branch |

⚠ **The first baseline will therefore carry a known bias and S5.f's entry must say so**:
`docs/SESSION_HANDOFF.md` §1 fact 6 — ignoring the taper over-calls gDNA by **11.0 %** genome-wide and by
**+0.36** in the last node before a polyA site.

### 18 per-face arrays became 5

`NodeGeometry` held a `_left`/`_right` pair for the unspliced mass, its integer flux, three effective
lengths and four spliced channels. The pairs existed because a boundary's two sides lay in
differently-sized flanks and so had **different divisors**. A contiguous edge is a 0-bp line. What went
with them:

* **`mass` versus `n`.** The old accumulator split one fragment's mass across objects, so the density
  numerator was fractional and a separate integer flux was carried for the Poisson variance. The new one
  deposits `+1` on every object touched — one `count`, and `Var(log ρ) = 1/n` is honest against it.
  ⭐ That alone dissolves 12 of the 18.
* **The junction-strand routing and the exon-bit flank gating**, plus the `_continues`/`_eff_spl_face`
  far-boundary machinery. All of it existed to *guess* which flank a spliced deposit belonged to, because
  the old accumulator attributed a splice to the node's edge rather than the junction's own coordinate
  (trap 6). The v8 index states `(src, dst, strand)`, so the guess is replaced by the chain's own
  adjacency: a junction's donor is `chain.right[src]`, its acceptor `chain.left[dst]`.
* **`BoundarySubstrate`'s `left_region`/`right_region` with their `-1` holes**, and the `max(left, right)`
  that de-duplicated the straddle. An edge always has a node on both sides (S5.d), so the branch has no
  cases and the `max` has nothing to de-duplicate.
* **`node_global_geometry`'s two-face sum**, and with it the long note about a `½` here silently
  cancelling a `½` missing from the per-face length — which is why that frame read the correct ρ while
  every per-face MESSAGE read ρ/2.
* In `bp_solver`: the `df`/`sf` face parameters through `_relay`, `_transport`, `_peel_share`,
  `_peel_share_scalar`, `_seam_pair` and `_flank_dom`; the `accept_l`/`accept_r` acceptor test; and
  `_rho_faces`' (node, left face, right face) triple, now one ρ_tot per slot.
  ⚠ **The scalar/vector twins are NOT merged** — `bp_solver` documents a measured **15.7×/op** reason.

### ⭐ The naming, after the owner's review

The first cut called them `count` / `mature_count`, and the owner asked whether `unspliced_count` /
`spliced_count` would be clearer. It would — and the question exposed a collision: `NodeStatics` already
held a `spliced_count` (`edge_spliced`), so renaming the junction flux to `spliced_count` would have put
**two fields of that name on two objects consumed side by side**, holding quantities that differ by two
orders of magnitude at the same line. That is trap 27 exactly.

Resolved by taking the accumulator's own three bank names, so one vocabulary runs from the executable
specification to the solver. At one line, of the molecules that touched it:

| field | population | strand axis | what it is |
|---|---|---|---|
| `unspliced_count` | `edge_unspliced` / `node_contained` | GENOME | crossed CONTIGUOUSLY, spliced nowhere — **the mixture being deconvolved** |
| `spliced_count` | `edge_spliced` | GENOME | crossed CONTIGUOUSLY, spliced elsewhere — certified RNA, the **floor** on this line's own population |
| `junction_count` | `sj_count`, gathered onto donor + acceptor | **TRANSCRIPT** | never crossed it, it **jumped** — the flux the graft hands an exon and the peel measures against |

⛔ **"mature" names nothing now.** It fits both spliced populations and so distinguishes neither.
Worked example, at a donor seam with 40 nascent read-throughs, 2 spliced-elsewhere crossings and 251
junction fragments: the peel needs 251 (using 2 puts the continuing share at 0.95 against a true 0.14),
while the strand floor needs 2 (using 251 asserts 251 certified-RNA fragments inside a population of 42).
Pinned by `test_spliced_count_and_junction_count_are_DIFFERENT_POPULATIONS` and
`test_the_word_MATURE_names_no_field`.

⭐ **`NodeStatics` is now structure only** — all three populations live on `NodeGeometry`, where their
difference is visible, and `build_node_statics` no longer takes a substrate at all. The two strand axes
stay named apart: `unspliced_count`/`spliced_count` are GENOME strand (where the read aligned, the
accumulator's one convention); `junction_count`/`eff_junction` are TRANSCRIPT strand, joined from each
junction's annotated strand — which *is* "sense derived, never stored".

⚠ **The rename is behaviour-neutral and was proven so**: the suite returned byte-identical counts
(266 / 1428 / 1 / 15) across it, and all 15 perturbations were re-run and still caught.

⭐ **A real defect fixed, not merely ported: the divisor is no longer floored.** The predecessor wrapped
every effective length in `max(·, _EPS)`. `effective_length`'s own contract and trap 23 say an object
with no opportunity must emit *nothing*; the floor is what produced densities of ~1e9 on the **12.4 %**
of fine-partition nodes where the contained effective length collapses to exactly 0.

### ⛔ THE GATE HAD A HOLE, AND ONLY PERTURBATION FOUND IT

23 tests in `tests/calibration/test_node_geometry.py`, written first and **verified failing**, gated on
**brute-force enumeration of integer start positions** in an explicit Python loop sharing no helper with
the implementation. Then the implementation was broken fifteen ways. **Fourteen were caught; P15 was
not:**

| | perturbation | caught |
|---|---|---|
| P1–P2 | node/edge populations swapped · genome-strand columns transposed | ✅ |
| P3 | **A7 accidentally ON** — the contiguous RNA divisor tapers | ✅ (the ruling is PINNED, so turning it on in S5.g must break this test) |
| P4 | the junction divisor ignores its reach | ✅ |
| P5 | mature flux keyed by the ALIGN column, not the junction's own strand | ✅ |
| P6 · P11 | deposits only on the donor line · acceptor taken from `src` | ✅ |
| P7 | the pooled junction divisor is a MEAN of ratios, not a ratio of sums | ✅ |
| P8 · P12 · P13 | the divisor floored to `_EPS` · zero divisor yields `inf` · `eff_gdna` floored | ✅ |
| P9 | the contained divisor collapses the pmf to its mean | ✅ |
| P10 · P14 | mature written onto NODE slots · junction flux takes one align column | ✅ |
| **P15** | **`slot_of_node = 2 * arange` — the LAYOUT ASSUMPTION instead of reading the chain** | ⛔ **NO — all 22 passed** |

**Why it hid.** Slot ids run `N E N E N` per reference, so *within the first reference* node `i` really
does sit at slot `2i`. Every fixture had one reference, and even the two-reference fixture put its
junction on chr1. chr2's node 3 sits at slot **5**, not 6 — so the assumption misplaces every junction
on every reference after the first, which at human scale is 285 of 286. Closed by
`test_a_junction_on_a_LATER_REFERENCE_lands_on_that_references_own_line`; all fifteen now fail the gate.

### Gates

* **`tests/calibration/test_node_geometry.py` — 23 new tests**, 15 perturbations caught.
* `tests/calibration/test_edge_flags.py` — 5 tests, replacing `test_boundary_flags.py`. The `k+1`
  boundary axis became a plain per-contiguous-edge array (`docs/calibration/S5_DESIGN_LOG.md` §4's ruling), so the
  "two spaces off by one per reference in opposite directions" problem is gone rather than handled.
* The whole solver **runs end to end** on the collapsed geometry — `init_beliefs` → `node_sweep` →
  `chain_node_deconv`/`chain_edge_deconv`, including every `_capture` diagnostic branch.
* ⭐ **The factor-1 bedrock invariant holds on the new geometry**: lay down a uniform ρ = 0.5 field and
  the strand-locked intergenic anchors read it back to 1e-9, and the AMBIG-recovery guard (the live
  9.0 %-low bound) still passes. That is the strongest evidence the port is faithful rather than merely
  green.
* Full suite **1440 passed / 256 failed / 1 xfailed / 15 errors** against a baseline of
  **1384 / 291 / 1 / 15** re-recorded from this tree in this session: **−35 failures, +56 passes**.
  ⭐ **Every remaining calibration failure is now S5.f's** — `test_calibrate` (9),
  `test_gdna_strand_integration` (4), `test_spliced_boundary_onesidedness` (3),
  `test_region_index_alignment` (3), `test_accumulator_span_unbiased` (3),
  `test_substrate_conservation` (2), `test_oracle` (2), `test_ambig_scenario` (2). Nothing in the tree
  still fails on the per-face model.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.

### ✅ `test_bp_solver.py` fully ported — 19 failures → 0, and the port found two hidden confounds

Three tests were of the per-face GEOMETRY itself and were **deleted**, not ported: their model does not
exist, so they could only be rewritten from scratch, and `test_node_geometry.py` is that rewrite (the
call S5.c made on `test_message_frames.py`). The other 16 are sweep-behaviour tests — the mature
graft/peel, the τ-gag fix, the intron relay, sweep determinism — and all 16 now run on
`_synthetic.make_chain_parts`. The `_pending_s5e` shim is deleted.

⭐ **Porting the fixtures exposed two artefacts the old shape had been hiding, and both are now gone
rather than commented around:**

* **`_mature_exon_chain` had to grow from 3 nodes to 5.** It placed mature flux on two intron↔exon
  *boundaries* by hand, which only worked because the old accumulator attributed a splice to the node's
  edge. A junction now states its own `(src, dst)`, so it must HAVE endpoints — the chain became
  `exon | intron | EXON | intron | exon` and the flux is *derived* from the graph. ⭐ The fixture is
  also more honest for it: `edge_spliced` is 0 on every seam, which is a measured fact rather than a
  convenience (mature RNA never crosses an exon↔intron seam — 0 of 1,146, `docs/SESSION_HANDOFF.md` §1
  fact 13).
* **`test_gdna_sweep_zero_gdna_pin_and_monotone` carried a ten-line comment explaining that its introns
  read ~0.44 instead of ~0.22 because the chain's two TERMINAL boundary slots were G1-locked and emitted
  structural all-gDNA into their neighbours** — "an ARTEFACT OF THIS ARTIFICIAL CHAIN". Those slots do
  not exist any more: `k` nodes own `k−1` lines. The invariant the test protects is unchanged and the
  artefact is gone with the shape that caused it.

### Files touched

* `src/rigel/calibration/node_geometry.py` — `NodeGeometry` (18 arrays → 5), `build_node_geometry`,
  `node_global_geometry`, `node_total_density`, `build_node_statics` (the region/boundary twin pair
  collapsed into one slot-keyed pass), `init_beliefs`, `_rate`, `_check_edge_flags`.
* `src/rigel/calibration/bp_solver.py` — the six per-face consumers; `chain_region_deconv` /
  `chain_boundary_side_deconv` → `chain_node_deconv` / `chain_edge_deconv` (per-EDGE, per the §4 ruling
  that killed the split-then-re-pool trap-2 pattern); `chain.order` → `arange`.
* `src/rigel/calibration/splice_graph.py` — `build_boundary_flags_array` → `build_edge_flags_array`;
  `JunctionGeometry` + `build_junction_geometry_arrays`, joining through `build_junction_edge_arrays`'s
  slot order rather than recomputing it (one implementation of a byte-identity contract).
* `src/rigel/calibration/node_init.py` · `density_model.py` (two arguments dissolved: `region_eff_len`
  and `fl_mean` are both the geometry's own `eff_gdna` now) · `calibrate.py` (names only; S5.f owns it) ·
  `src/rigel/pipeline.py` · `src/rigel/scan_cache.py`.
* `tests/calibration/_synthetic.py` — `make_chain_parts`, the shared node/edge/junction fixture.
* `tests/calibration/test_node_geometry.py`, `test_edge_flags.py` — new.
  `test_boundary_flags.py` deleted.

---

## S5.f — `calibrate()` runs, and there is a number (2026-07-30)

⭐ **THE PIVOT.** `calibrate()` had not run since the accumulator rewrite began; nothing downstream was
measurable, so nothing downstream could be judged. It runs now, end to end, on all eight chr22 pilot
conditions — and **its numbers are THE FIRST BASELINE**, recorded below, against which every deferred
derivation (S5.g, S5.a2, A6/A3) is to be A/B'd.

**Baseline re-recorded from this tree in this session** (trap 17, never a stored one):
**1440 passed / 256 failed / 1 xfailed / 15 errors**, matching the S5.e handoff exactly.
**After: 1729 passed / 22 failed / 1 xfailed / 0 errors** — **−234 failures, +289 passes, −15 errors.**
All 22 remaining are the golden files (21) plus one EM-nondeterminism finding, both recorded below.

### ⛔ The owner's ruling, taken FIRST because it decides the schema

> **Does `CalibrationResult` gain a third (junction) axis?** — **Yes: exported, QC only.**

The accumulator has three populations at a line: `edge_unspliced` (the mixture being deconvolved),
`edge_spliced` (crossed contiguously, spliced elsewhere) and `sj_count` (jumped). The junction flux was
already *consumed* — `build_node_geometry` gathers it onto donor/acceptor lines for the graft/peel — but
never *exported*, so the calibration's own output was silent about the population that at a donor seam is
251 of 253 fragments. It is now `mass_rna_junction`, carried verbatim (a junction edge is pure mature RNA
by construction, so there is nothing to deconvolve).

⚠ **`assemble_priors` does NOT read it, and that asymmetry is deliberate.** The prior arbitrates only the
*unspliced competing* mass: a spliced fragment has no gDNA candidate in the EM, which is exactly why
`mass_rna_spliced_edge` is withheld. A junction fragment is certified RNA in the same sense, so feeding it
in would load the RNA side of a split it does not belong to — and a locus whose RNA is fully spliced
*should* get a near-zero `rna_prior_count`, because its unspliced fragments really are gDNA or nascent.
Pinned by `test_the_junction_flux_does_NOT_enter_the_rna_prior`.

### ⭐ THE FIRST BASELINE — chr22 pilot, 8 conditions

`scripts/design/build_scan_cache.py --force` (66.6 s for all eight), then `calibrate()` per condition.
Axes: **N = 22,231 nodes · E = 22,138 edges · J = 8,471 junctions**. ⚠ **Bit-identical on re-run** —
every digit below reproduces, so this is a baseline rather than a sample.

| condition | s | `rho_g` | κ | `od_g` | `od_r` | gDNA | RNA | junction | **f_gdna** |
|---|---|---|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 6.4 | 0.102275 | 0.5000 | 0.0000 | 0.0000 | 5,370,056 | 5,373,548 | 2,117,198 | **0.4998** |
| gdna100 ss0.50 capture_on | 5.2 | 0.084227 | 0.4990 | 0.0000 | 0.0000 | 4,384,623 | 7,295,911 | 1,315,320 | **0.3754** |
| gdna100 ss0.99 capture_off | 6.4 | 0.095879 | 0.0101 | 0.0000 | 0.0000 | 5,034,150 | 5,704,854 | 2,114,778 | **0.4688** |
| gdna100 ss0.99 capture_on | 5.2 | 0.109805 | 0.0101 | 0.0000 | 0.0000 | 5,714,602 | 5,955,867 | 1,316,066 | **0.4897** |
| **none** ss0.50 capture_off | 3.7 | 0.008417 | 0.5002 | 0.0345 | 0.0000 | 443,277 | 5,145,770 | 2,117,047 | **0.0793** |
| **none** ss0.50 capture_on | 2.8 | 0.001066 | 0.4998 | 0.0345 | 0.0000 | 56,141 | 5,502,246 | 1,314,146 | **0.0101** |
| **none** ss0.99 capture_off | 3.6 | 0.000318 | 0.0100 | 0.0345 | 0.0000 | 16,745 | 5,568,579 | 2,118,177 | **0.0030** |
| **none** ss0.99 capture_on | 2.8 | 0.000171 | 0.0099 | 0.0345 | 0.0000 | 9,021 | 5,554,579 | 1,315,272 | **0.0016** |

Truth, from the suite manifest: `gdna_rate` 1.0 ("gdna100", ⇒ **f_gdna ≈ 0.5**) or 0.0 ("none", ⇒
**f_gdna = 0 exactly**); `gdna_strand_overdispersion = 0.0`; `strand_specificity` 0.5 or 0.99.

**What it says, read against that truth:**

* **`od_g` recovers 0.0000 exactly** on every gdna100 condition — the simulator generated gDNA with zero
  strand overdispersion and the fit returns zero, not a floor. On the zero-gDNA conditions it falls back
  to `_PRIOR_OVERDISPERSION` = 0.0345 on all four, which is correct: there is no gDNA to fit from.
* **`od_r` = 0.0000 everywhere**, consistent with the simulator being Poisson by construction
  (`docs/SESSION_HANDOFF.md` §1 fact 18) — so this baseline can say nothing about dispersion, by design.
* **f_gdna lands within 6 % of truth on 3 of 4** gdna100 conditions.
* ⚠ **The zero-gDNA false positive is 7.9 % in the worst corner** (unstranded, capture off) and falls to
  0.16 % when stranded and captured. That ordering is what §2 predicts — at κ = ½ the strand channel
  carries *exactly zero* information about composition, so the solver has only count and density left.
  ⛔ Quote it **only alongside the gdna100 column**: a zero-gDNA guard is ONE-SIDED (trap 19) and any
  change that merely lowers gDNA scores better on it.

### ⚠ TWO ANOMALIES THIS BASELINE SURFACES — neither is S5.f's, both are now measurable

1. ⛔ **The fitted κ is `1 − truth`.** A library simulated at `strand_specificity = 0.99` calibrates to
   **κ = 0.0101**; at 0.50 it calibrates to 0.4990–0.5002 (where the flip is invisible). This is
   `docs/SESSION_HANDOFF.md` §0 **C4** — "sense fraction is 0.002–0.012 on all four real cfRNA libraries
   (nearly fully antisense) — possibly a read-orientation convention bug" — now answered against ground
   truth for the first time: it is a **convention flip, not biology**. S5.f does not touch
   `fit_strand_balance`, so this is pre-existing and merely became visible. **It needs its own step.**
2. ⚠ **`gdna100 ss0.50 capture_on` reads f_gdna 0.3754, 25 % low**, and its RNA total (7.30 M) is 21–36 %
   above every other condition's. The other three gdna100 conditions are within 6 %. Unexplained; it is
   the both-worst-case corner (no strand information *and* capture), so it is the natural first target.

### ⚠ THE BASELINE CARRIES A KNOWN BIAS, AND IT IS NAMED

**A7 is OFF at contiguous edges by owner ruling** (`docs/calibration/S5_DESIGN_LOG.md` §1 A7): the RNA half of an unspliced
crossing takes `UNBOUNDED_REACH` rather than its transcript's real remaining length. Measured cost,
already on record (`docs/SESSION_HANDOFF.md` §1 fact 6): an **11.0 % genome-wide gDNA over-call**, contiguous
seams worse than junctions (0.750 vs 0.886), and the gDNA fraction off by **+0.36** in the last node
before a polyA site. Junction edges *do* take their real exonic reach — a brand-new population with no
predecessor divisor, so wiring it regresses nothing. **S5.g turns the contiguous taper on and A/Bs it
against the table above**; that ordering is the entire reason A7 was deferred.

### ⭐ WHAT DISSOLVED, AND THE TRAP IT TAKES WITH IT

`CalibrationResult` carried **six** per-region mass arrays and a per-side length; it now carries three
axes and two supports.

| gone | successor |
|---|---|
| `mass_{gdna,rna}_{left,right}` | ⭐ **`mass_{gdna,rna}_edge`** — one number per 0-bp line |
| `gdna_boundary_len` (`E[min(ℓ,L)]/2`) | ⛔ **none.** Its ½ existed only because a face's mass was half a crossing. The per-edge divisor is `crossing_eff_length`, and there is no ½ in it |
| `priors._gdna_region_node_arrays` → `_pooled_seam_arrays` | ⭐ the pooling **disappears**, it is not re-derived |
| `region_arrays.{region_boundary,boundary_region}_indices` | `edge_node_indices` / `node_right_edge` — the last `k+1` axis in the tree |
| `strand_deconv.boundary_side_seeds` (+ `_SideQuantities`, `_compute_side`, `_left_right_neighbors`) | `edge_seeds` — **ONE seed per line, not two per boundary** |

⛔ **The split-then-re-pool was a no-op with a history.** `assemble_priors` summed
`mass_gdna_right[r] + mass_gdna_left[r+1]` straight back together, and `_pooled_seam_arrays` did the
identical thing — putting back a split the calibrator had just made. That exact sum-then-halve pattern
hid an **exact factor of 2 for months** (`docs/SESSION_HANDOFF.md` §3 trap 2). It cannot recur: there is no face
to halve and no pair to sum.

⭐ **The evidence that it is really gone is the FIXTURES.** `test_capture_eff_length.py` carried **three
near-identical result builders** that differed in exactly one thing — whether `gdna_boundary_len` was
halved — two of which stored an un-halved length *and* deposited half the mass, cancelling exactly and
hiding the factor-2 from every assertion in the file. They are now **one** builder. `test_priors.py`'s
uniform-field fixture carried a ten-line note explaining the same convention; it is now three lines
(`mass = ρ · support`, both axes). And **every span in `test_priors.py` is byte-identical to its
pre-S5.f value** — 640 / 650 / 700 / 400 / 850 — so the schema moved and the geometry did not.

⭐ **Two further dissolutions, each with its derivation:**

* **The edge seed weight is exactly 1.** It was `clip(density · eff / mass)` where a count-observable
  side read its own crossing density `mass / eff` — algebraically 1 on every seed the mask admits, with
  the borrowed-density branch unreachable. With `count` and `mass` the same number that is exact, so the
  ratio, its effective length, and the `boundary_side_eff_len` argument that carried it all go.
* **`derive.gdna_density_global`'s terminal masking.** It masked each boundary side with
  `same_ref_left_right` because a reference terminal has nothing on its far side. `E = N − n_refs` and
  every entry is a real object, so there is nothing to exclude.

### ⛔ THE GATE HAD A HOLE, AND ONLY PERTURBATION FOUND IT

Nine perturbations across the two new pieces; **eight caught, one not**:

| | perturbation | caught |
|---|---|---|
| P1 | edge id = LEFT NODE id — correct on reference 0 only (**the P15 shape**) | ✅ |
| P2 · P4 | the grouped-`ref_id` check removed · every adjacent pair an edge (lines straddle references) | ✅ |
| P3 | `hi_node = lo_node` — an edge collapses onto one node | ✅ |
| P5 | `sense` always takes the POS column (the NEG orientation ignored) | ✅ |
| P7 | count-observability ignored — a mid-exon line seeds the gDNA fit | ✅ |
| P8 · P9 | the shared-exon test becomes an OR · edge observability read at the NODE index | ✅ |
| **P6** | **an AMBIG flank no longer blocks the seed** | ⛔ **NO — all 29 passed** |

**Why it hid, and what it exposed.** `TS_AMBIG` is a fourth distinct value (3), so
`(ts == TS_POS) | (ts == TS_NONE)` is already `False` on it — the `~either_ambig` guard the predecessor
carried on both clauses was **dead code reading as the rule**. Deleting it changed nothing, which is
worse than absent: it invites the belief that the rule lives there. The guard is removed and the rule is
now pinned by `test_an_AMBIG_flank_cannot_seed`, whose fixture uses INTRON bits on both strands so the
flank stays *count*-observable and the test can only fail on strand. Re-armed
(`TS_AMBIG` admitted as POS-compatible), the new test fires.

### ⭐ CALIBRATION IS DETERMINISTIC; SOMETHING DOWNSTREAM IS NOT

`tests/scenarios/test_nrna_double_counting.py::test_full_sweep[g20_n0_s100]` fails its negative-control
limit (25 counts) — and, run in isolation five times, returns **29 / 29 / 32 / 29 / 32**. Fixing
`PYTHONHASHSEED` does not stabilise it. Localised:

| stage | two runs, same input | verdict |
|---|---|---|
| **scan** — all five banks + `node_start_count` | `np.array_equal` on every one | ✅ bit-identical |
| **calibrate** — same payload, two calls | `max|diff| = 0.000e+00` on every array; `gdna_density_global` identical to the last digit | ✅ bit-identical |
| **calibrate** — two INDEPENDENT scans | `max|diff| = 0.000e+00` | ✅ bit-identical |

⭐ So **the first baseline is reproducible** and the nondeterminism is downstream of calibration, in the
per-locus EM / quantification. ⚠ **It bounds every future end-to-end A/B**: a transcript-count comparison
is only meaningful above a noise floor of ~3 counts on a ~30-count control, and that floor has not been
characterised. Calibration-level A/Bs (S5.g's) are unaffected. **Not S5.f's, not chased here** — it is a
pre-existing property that only became observable once `calibrate()` ran.

### Gates

* **`calibrate()` runs end to end on all 8 chr22 pilot conditions**, 2.6–6.4 s each, and the numbers
  above are bit-identical on re-run. ⭐ Every number in `test_calibrate.py` is *arithmetic on the
  fixture's own banks* — node totals 12/21/15, edge totals 5/11 (unspliced + spliced), junction 13,
  divisors 100−50+1 = 51 contained and 50−1 = 49 crossing — not a recorded output.
* **`tests/calibration/` is fully green: 543 passed, 0 failed** (was 25 failed at the start of the port).
* `test_region_arrays.py` — the node↔edge mapping re-derived by a **different algorithm**
  (`build_node_chain` walking the payload's CSR offsets), on a deliberately **multi-reference** topology
  with a 3-node, a 2-node and a 1-node reference. Trap 1: a validator calling the builder's own helper
  validates nothing.
* `test_result_schema.py` — **31 tests**, every array pinned to its own axis against all three lengths,
  because `E = N − n_refs` makes a mis-keyed array a *plausible* shape that nothing else would fault on.
* **The overdispersion end-to-end recovery still holds**: a planted BetaBinom(½, od) over 400 intergenic
  nodes is recovered at od ∈ {0.05, 0.10, 0.20} through the rewired path, and the Binomial limit still
  floors below 0.02.
* **The crossing-vs-contained unbiasedness identity still holds** on a real simulated scan
  (`docs/SESSION_HANDOFF.md` §1 fact 4) — now as a ratio of SUMS on both axes rather than a mean of per-boundary
  ratios, and with a line's ½-averaged two faces replaced by its one count.
* **The oracle validates every bank EXACTLY, on all three axes.** Sum-to-full used to be byte-exact on
  two arrays and *approximate* on the two float32 mass arrays; every bank is an integer count now, so
  `boundary_mass_tol` is **deleted rather than carried at zero** — a tolerance that must be zero is a
  claim that some comparison is approximate, and none is. "gDNA is never spliced" is now checked on the
  junction axis too, which the old 4-channel layout had no room for.
* `ruff check src/ tests/ scripts/` and `ruff format` clean.

### Deleted rather than ported

* **`tests/calibration/test_spliced_boundary_onesidedness.py`** — its entire subject is that spliced mass
  sits on the exon FACE of a boundary and zero on the intron-facing one. There are no faces, and a splice
  now deposits on its **junction edge only, never on the lines it splices over** — a stronger statement,
  already gated at the specification level (`test_a_spliced_jump_deposits_NOTHING_on_the_lines_it_splices_over`,
  verified passing before the deletion). The one part that survives — the two edge banks are disjoint
  populations — is re-pinned in `test_substrate_conservation.py`. Same call S5.c made on
  `test_message_frames.py` and S5.e made on three `test_bp_solver.py` tests.
* **`test_capture_eff_length.py`'s two D6 tests** — falsification tests for the pooled-seam factor-2, in
  which *every quantity named* is gone. The property they protected is kept as
  `test_a_crossing_object_under_a_uniform_field_reads_RHO` and
  `test_a_line_below_the_reference_density_CONTRACTS_rather_than_clipping`; what is retired is the
  specific arithmetic that could break it.
* **`test_substrate_conservation.py` was rewritten, not ported** — it asserted the region→boundary map
  attributes "the non-terminal boundary mass, no double count, no loss", a property of the `k+1` axis
  whose two data-free terminals were the only thing that could lose mass. What is kept is what made it
  worth having: it runs on a **real scan**, and its invariants are re-derived from the FRAGMENT BUFFER
  and the INDEX — sources independent of the accumulator. ⚠ Its scenario gained a gDNA arm, because
  without one every fragment is mature RNA that either sits inside an exon node or jumps, so
  `edge_unspliced` — the one population calibration deconvolves — was identically zero and the
  "conservation" test passed over an empty bank.

### Files touched

* `src/rigel/calibration/result.py` — the three-axis schema, per-axis validation.
* `calibrate.py` — the body; the chain + geometry now built FIRST, so the geometry owns every divisor and
  the four length models the predecessor passed around (`region_eff_length`, `boundary_side_eff_length`,
  `boundary_eff_length`, a scalar `fl_mean`) collapse to one projection off `eff_gdna`.
* `priors.py` (`_gdna_node_arrays`, `edge_owner_nodes`) · `capture_eff_length.py`
  (`_left_keyed_edge_arrays`) · `derive.py` · `strand_deconv.py` (`edge_seeds`,
  `edge_strand_orientation`) · `gdna_strand.py` (`_node_seeds`) · `density_model.py`
  (`count_observable_masks` → the edge axis) · `background_reference.py` · `density_deconv.py` ·
  `region_arrays.py` · `track.py` · `src/rigel/pipeline.py` · `src/rigel/cli.py`.
* Tests: `test_result_schema` `test_calibrate` `test_priors` `test_capture_eff_length` `_oracle`
  `test_oracle` `test_substrate_conservation` `test_gdna_strand` `test_gdna_strand_integration`
  `test_region_arrays` `test_region_index_alignment` `test_background_reference`
  `test_accumulator_span_unbiased` `test_ambig_scenario` `_synthetic`
  `tests/test_scanner_accumulator_integration.py` `tests/test_summary_report.py`.

### ⚠ Left for S6, and one tool defect found on the way

* **21 golden files** are the only remaining failures besides the EM finding. They move NUMERICALLY, not
  merely start running, and are regenerated **once**, at S6 — `pytest tests/ --update-golden`.
* ⛔ **`scripts/design/build_scan_cache.py`'s skip logic does not cover the payload schema digest.** It
  printed `cached / skip` for all eight pilot conditions that `read_scan_cache` then **refused**
  (`payload_schema_digest None != '66b41ea0b645209d'`). `--force` is the workaround; the check belongs in
  the skip test. Same family as `docs/SESSION_HANDOFF.md` §3 trap 25 — a cache key that does not cover the
  artifact it caches.

---

## S5.f-addendum — the κ mirror is CONSISTENT, so the inference is correct (2026-07-30)

**Not a step; a measurement made against S5.f's baseline, recorded here rather than by editing S5.f's
entry.** The S5.f entry reported the fitted κ as `1 − truth` and said the flip "needs its own step"
without establishing how severe it is. It is now established, and it is **much less severe than it
looks**.

**The question.** κ enters the only intrinsic gDNA/RNA signal — `p = ½·f_g + κ·(1−f_g)`
(`docs/SESSION_HANDOFF.md` §2). If κ were mirrored while the per-node sense counts were not, the strand
likelihood would score observations against a mirrored `p` and **the one channel that separates gDNA
from RNA would be read backwards**. If both are mirrored, the mirror cancels and only the label is wrong.

**The experiment.** Inject κ via `InjectedCalibrationPriors.rna_sense_frac` on the two **zero-gDNA**
stranded pilot conditions, where truth is `f_gdna = 0` **exactly**:

| κ used | `gdna_none_ss_0.99_capture_off` | `gdna_none_ss_0.99_capture_on` |
|---|---|---|
| **0.0100 — the FITTED value** | **0.0030** | **0.0016** |
| 0.99 — the simulated "truth" | ⛔ **0.4992** | ⛔ **0.4822** |
| 0.50 — uninformative (the control) | 0.0792 | 0.0100 |

⭐ **Forcing the nominal truth is 166× worse than the fitted value, and 6× worse than supplying no strand
information at all.** So κ and the per-node sense columns are mirrored the same way, the mirror cancels
inside the solve, and **the deconvolution is correct**. The defect is confined to the exported scalar.

⚠ **The control is what makes this conclusive.** Without the κ = 0.50 row, "0.99 is worse" could have
been a mis-specified-prior effect of any kind; landing *between* the two — 0.0030 correct, 0.0792
uninformative, 0.4992 anti-informative — is the signature of a sign error rather than a noisy fit.

**Consequences recorded:** `TODO.md` §6 (ranked 6, with the gate any fix must pass — κ must move to 0.99
*while* `f_gdna` stays at 0.0030); `docs/SESSION_HANDOFF.md` §0 C4 (answered) and §1 fact 17 (its κ column reads
backwards — the four real cfRNA libraries are ordinary **highly stranded** libraries, not near-purely
antisense ones, which is a materially different claim about the data).

⚠ **What this does NOT settle:** where the mirror lives — the scanner's genome-strand convention or the
simulator's read orientation. That is `TODO.md` §6's first step.

---

## S6 — the goldens regenerated, and two defects the regeneration exposed (2026-07-30)

**Gate: suite green; goldens regenerated ONCE.** Suite **1729 passed / 22 failed → 1752 passed /
1 failed**. The one remaining failure is **not** S6's and **not** flakiness — see below.

### The delete half was already done, and that is the point

S6 is specified as "delete; `ruff check` undefined-name failures are the authoritative list"
(`docs/accumulator/IMPLEMENTATION_PLAN.md` §469). That list is **empty**: `ruff check --select F821` passes, there is not
one `NotImplementedError` shim left in `src/` or `tests/`, and no `_pending_*` markers. S5.a–S5.f each
deleted as they went, which is what "converge and delete" is supposed to produce — a step whose own
deletion phase finds nothing.

⭐ **So the dead-code hunt went past what ruff can see**: an AST sweep for module-level symbols defined
in `src/rigel/` and never referenced anywhere in `src/`, `tests/` or `scripts/`. **One hit**, and it was
not dead code — it was an **uncalled safety check** (below).

### ⛔ TWO DEFECTS, BOTH INTRODUCED BY S5.f, BOTH INVISIBLE TO A GREEN SUITE

**1. `scan_cache.calibration_inputs` was broken by S5.f.** Its docstring says "Exactly the keyword
arguments `calibrate` needs"; `calibrate` gained `junctions` at S5.f and this helper was never updated.
Any caller would trip S5.f's own alignment guard. ⚠ **The existing test could not catch it** — it checked
a *hand-kept* `required` set, and a hand-kept list cannot notice an argument nobody added to it. Replaced
with the set of every parameter of `calibrate` that has no usable default, read off `inspect.signature`,
so the next added argument fails here instead of at a call site. A second test now **actually calls**
`calibrate(**calibration_inputs(...))`, because a signature check sees a missing NAME and never a
mis-sized ARRAY — and a junction axis of the wrong length places every splice on the wrong line.

**2. `check_scan_config` had no caller, while `read_scan_cache`'s docstring promised what it does.**
`read_scan_cache` said it would "REFUSE it unless it describes this index **and this scan
configuration**". It checked the manifest's `scan_config_digest` against **the config the manifest
itself records** — self-consistency — and never against the config the caller intends. So a cache
scanned under `sj_strand_tag="XS"` loaded silently for a caller who wanted `"auto"`. Two statements
about one contract, disagreeing (`docs/SESSION_HANDOFF.md` §3 trap 27). ⭐ Resolved by **wiring it, not deleting
it**: `read_scan_cache` takes an optional `scan_config` and calls the check when given one. Deleting a
safety check because nobody calls it is the wrong way to resolve a promise the prose already made.

⭐ **A strict xfail was xfailing for the WRONG REASON.**
`test_population_priors_can_be_extracted_from_a_cached_scan` was `xfail(strict=True)` naming S5 as the
blocker. `calibrate()` has run since S5.f — so it should have started passing, and it did not, because
of defect 1. **Strict is what made this visible at all** (a non-strict xfail would have stayed quietly
red forever), but strictness proves only that *something* fails, never that the recorded cause is still
the cause. The xfail is removed and the test is live. ⭐ This closes `TODO.md` §4 (the scan cache's
population-prior seed) as a side effect.

### ⭐ THE TWICE-AND-DIFF EARNED ITS PLACE, AND THEN DISPROVED ITS OWN PREMISE

The goldens were regenerated **twice** and the two outputs diffed, on the reasoning that the EM's default
`assignment_mode="sample"` draws from the posterior and a flaky expectation baked into a golden is
permanent. **2 of 21 scenarios differed between regenerations** — `wide_intron` and
`higher_frag_count_3k`, exactly the two that had been failing.

⭐ **But the difference is 1–2 ulp of float64, not a resampled assignment.** The COUNTS are identical
(`331.0`, `0.0`); what moves is `687.9988063212445` vs `687.998806321244` and
`623458.9808348428` vs `623458.9808348427` — ~**1e-15** relative, float accumulation order, the same
family as `docs/SESSION_HANDOFF.md` §1 fact 11. Against `RTOL = 1e-6` that is **nine orders of magnitude of
headroom**, and the regenerated goldens pass three consecutive runs. ⚠ Worth recording *because* it came
out negative: the check cost two minutes and converted "the goldens are probably fine" into a measured
statement, on precisely the two scenarios where it mattered.

### ⛔ THE ONE REMAINING FAILURE IS A MODELLING QUESTION, NOT A TEST DEFECT

`tests/scenarios/test_nrna_double_counting.py::test_full_sweep[g20_n0_s100]` — a silent transcript
(truth 0) against a limit of 25. Measured across all three assignment modes:

| `assignment_mode` | `t_ctrl` counts | `tests/scenarios/` (120 tests) |
|---|---|---|
| `sample` (the default) | 29 / 29 / 32 / 29 / 33 | 1 failed |
| `fractional` | 30 | 1 failed |
| **`map`** | ≤ 25 | ⭐ **120 passed** |

⭐ **It does not fail because the answer wobbles — every draw is over the limit.** The variation is real
but incidental; the leak is ~30 counts under two of three modes.

⛔ **`map` was NOT adopted, and adopting it would be a mistake.** MAP assigns each fragment to its single
highest-posterior component, so it is the mode that most suppresses low-posterior assignment — and a
negative control is a **one-sided** metric (`docs/SESSION_HANDOFF.md` §3 trap 19: in a library with none of X,
any change that lowers X scores better). Going green under MAP would be fitting the test to the mode.
Whether ~30 counts of leak onto a silent transcript is acceptable is the **owner's modelling call**,
recorded in `TODO.md` §7.

### ⭐ And the practical consequence for every A/B after this

Owner, 2026-07-30: the EM samples from the posterior **by design**. Measured over four runs of one
scenario at a fixed seed, varying only the mode: `sample` spreads **~0.5 %** (5440.7 / 5458.8 / 5465.3 /
5455.8), `map` **~1e-10** (6002.246879554673 / …736 / …669 / …744), `fractional` **~1e-11**. ⭐ **So an
end-to-end A/B has a switch** — run both arms under `map` or `fractional` and the noise drops eight
orders of magnitude. ⚠ The three modes give materially different totals (5441 / 6002 / 6277), so they are
different estimators: hold the mode fixed across both arms, and never compare a number from one mode to
one from another.

⚠ **A seeded `sample` run still does not reproduce**, and that is a genuine defect — the seed is plumbed
(`estimator.py:187` → `:369` → `em_solver.cpp:2165`'s `SplitMix64(rng_seed ^ li·φ)`) and
`tests/scenarios/conftest.py:92` already sets one. `TODO.md` §7 ranks it and lists the candidates.

### Files touched

* `tests/golden/` — **147 files regenerated** (21 scenarios × transcript/gene/loci/scalars, feather+tsv).
* `src/rigel/scan_cache.py` — `index_derived_inputs` gains `junctions`; `read_scan_cache` gains an
  optional `scan_config` and calls `check_scan_config`; `check_scan_config` exported; the module
  docstring's step-4 block corrected (it still said the seed path was blocked on S5).
* `tests/test_scan_cache.py` — the hand-kept `required` set replaced by `inspect.signature`; a new test
  that actually drives `calibrate`; the obsolete strict xfail removed.
* `docs/TODO.md` (§6, §7 + the ranked critical path) · `docs/SESSION_HANDOFF.md` (§0 C4, §1 fact 17).

---

## S5.g-1 — the per-contiguous-edge RNA reach, built and gated (2026-07-30)

**The input A7 needs, landed on its own.** S5.g turns on the contiguous-edge RNA taper; this is the
array that makes it possible, split out because it has an unambiguous correct answer (it must match
`edges_df` on the payload's edge axis) while the divisor change does **not** — see the open question
below. Suite **1752 → 1758 passed / 1 failed** (+6 tests, the failure unchanged and unrelated).

`splice_graph.build_contiguous_edge_reach_arrays(index) → (reach_lo, reach_hi)`, each `float64[E, 2]`,
column 0 the POS-strand transcript's reach and column 1 the NEG's. Keyed by ``src`` and laid out per
reference exactly as `build_edge_flags_array` is, so the two are **the same axis element for element**.

**Three properties that are not arbitrary, each with its citation:**

* **PER STRAND and per SIDE** — `docs/SESSION_HANDOFF.md` §2: reach is "maximised over transcripts
  independently per side AND per strand". This is what forces `eff_rna` to gain a strand axis; a single
  averaged reach describes neither transcript.
* **GENOMIC, not exonic** — unlike a junction's. A junction is used only by a spliced molecule so what
  remains either side is exonic; a contiguous line is also crossed by *nascent* RNA, which is genomic.
  Taking the exonic reach here would declare an intronic nascent fragment impossible.
* ⚠ **A reach of 0 is the ANSWER, not a missing value.** No template of that strand at that line ⇒ zero
  opportunity ⇒ divisor 0 ⇒ the consumer emits nothing (`docs/SESSION_HANDOFF.md` §3 trap 23). Measured on the
  chr22 pilot index: **40.6 %** of contiguous edges have no POS template, **42.9 %** no NEG.

### ⭐ It independently reproduces fact 6's magnitude

On the chr22 pilot index with the suite's own fragment-length distribution (mean 206, sd 98), over the
lines where the strand has any template at all:

| | unbounded divisor | tapered mean | ratio |
|---|---|---|---|
| POS | 216.6 | 182.5 | **0.8425** |
| NEG | 216.6 | 186.1 | **0.8592** |

`docs/SESSION_HANDOFF.md` §1 fact 6 records a bp-weighted **0.8904** over 3,000 genes, and notes contiguous
seams are the worse half (0.750 against junctions' 0.886). Different index, different weighting,
different FL — and the same magnitude. ⭐ That is an independent reproduction of the number A7 exists to
remove, from an array built for a different purpose.

### Gates

* `tests/calibration/test_edge_reach.py` — **6 tests**, written first and verified failing, with the
  edge coordinate derived from the PARTITION rather than from the reach builder (trap 1).
* **Six perturbations, all caught**: strand columns transposed · keyed by `dst` not `src` · `lo`/`hi`
  swapped · the `k−1` slice becomes `k` (a terminal slot returns) · zero reach floored to
  `UNBOUNDED_REACH` · reach clipped to a one-fragment minimum.
  ⚠ The flooring perturbation was **placed wrongly the first time** — after the slices were already
  taken — and passed. That was a no-op, not a gate hole; re-armed where it bites, it fails 4 of 6. A
  perturbation that changes nothing proves nothing, and checking *that* is part of the discipline.
* The fixture is deliberately **multi-reference** with a POS and a NEG transcript that are **disjoint**,
  so a transposed strand axis moves a non-zero to the wrong line rather than permuting a pair.

### ⛔ OPEN — the one question that gates the rest of S5.g, and it is the owner's

Turning the taper on means `NodeGeometry.eff_rna` becomes `[n_slots, 2]`. Two of its three consumers
take that without a decision:

* `node_init._rna` is **already called once per strand** — it just needs the matching column.
* `node_geometry.node_total_density` sums the two strands, exactly as it already does for junctions.

⛔ **`bp_solver.py:334` does not.** There `E_r` feeds the transfer variance
`Var(log ρ_tot) = 1/n + [(1/E_g − 1/E_r)/B]²·Var(f_g)` — a **strand-agnostic scalar**. Collapsing two
strand divisors into one for that term requires a reduction rule, and the honest candidates disagree:

| candidate | argument for | argument against |
|---|---|---|
| the belief-weighted divisor `M·(f_pos+f_neg) / ρ_r_total` | the only one that reproduces the actual total RNA density | belief-dependent, and this term is evaluated at the INPUT belief — a self-reference of the trap-11 family |
| max over strands (the loosest reach) | describes the strand with template, and at the 40 % of lines with only one live strand it IS that strand's reach | ⛔ **no monotone safety argument** — see the correction below |
| leave it at `UNBOUNDED_REACH` | minimal change; it is a damping heuristic, not a density | ⛔ puts two meanings on `eff_rna` — trap 27, the exact defect S5.e removed |

⛔ **CORRECTION, made before this entry was committed.** An earlier draft justified "max over strands"
as *"conservative: least damping, never over-confident"*. **That is false**, and the arithmetic says so:
the term is `(1/E_g − 1/E_r)²`, which is **not monotonic in `E_r`** — it is ZERO at `E_r = E_g` and grows
on both sides. On the suite's own pools (`E_g = mu_g − 1 = 156`, `E_r` unbounded `= 205`):

| `E_r` | 205 (unbounded) | 190 | 170 | **156** | 140 | 100 | 50 |
|---|---|---|---|---|---|---|---|
| `(1/E_g − 1/E_r)²` | 2.35e-6 | 1.32e-6 | 2.79e-7 | **0** | 5.37e-7 | 1.29e-5 | **1.85e-4** |

The zero sits **inside the plausible tapered range**, so tapering first *reduces* the damping to nothing
and then increases it sharply — at a deep terminus (`E_r = 50`) it is **79×** the unbounded value. No
"larger is safer" argument survives that.

⭐ **And the zero is correct physics, not a pathology**: when the two components convert mass to density
identically, composition uncertainty cannot move the total density, so its variance genuinely vanishes.

**This is a new heuristic, and the standing rule is to stop and discuss before adding one.** It is
recorded rather than guessed. `TODO.md` §1.

---

## S5.g-2 — the A7 taper A/B: MEASURED, and it does essentially nothing (2026-07-30)

⛔ **The owner challenged the concept of reach, and the measurement backs the challenge.** A7 was
deferred past S5.f specifically so it could be A/B'd against a baseline. It has been, and **turning the
contiguous-edge RNA taper on moves the library gDNA fraction by ≤ 0.0002.**

### The gate first: arm A reproduces the committed baseline EXACTLY

`edge_rna_reach` was added to `build_node_geometry`/`calibrate` as ONE optional argument — `None` keeps
`UNBOUNDED_REACH`, a `(reach_lo, reach_hi)` pair turns the taper on — so the two arms share every line of
code and differ in one value. **All 8 conditions of arm A are bit-identical to S5.f's recorded baseline**
(trap 17: if HEAD-vs-baseline is not 100 %, the baseline is what is broken).

Arm B uses the **max over strands** reach. That is the physically right single-number approximation — an
unspliced crossing could belong to either strand's transcript, so its opportunity is the better of the
two — and it is also the *least* aggressive taper, hence a **lower bound** on the per-strand version's
effect. If the generous version moves nothing, the refined one cannot justify its machinery.

### The result

| condition | truth | ARM A | ARM B | Δ |
|---|---|---|---|---|
| gdna100 ss0.50 capture_off | ~0.50 | 0.4998 | 0.4999 | +0.0001 |
| gdna100 ss0.50 capture_on | ~0.50 | 0.3754 | 0.3756 | +0.0002 |
| gdna100 ss0.99 capture_off/on | ~0.50 | 0.4688 / 0.4897 | 0.4688 / 0.4897 | +0.0000 |
| **none** ss0.50 capture_off | **0.0000** | 0.0793 | 0.0793 | −0.0000 |
| **none** ss0.50 capture_on | **0.0000** | 0.0101 | 0.0101 | −0.0000 |
| **none** ss0.99 capture_off/on | **0.0000** | 0.0030 / 0.0016 | 0.0030 / 0.0016 | −0.0000 |

Direction is *toward* truth on all four zero-gDNA conditions — by less than 1e-4.

### ⭐ THE NULL IS REAL — the taper did reach the solve

A null result is worthless if the change never arrived. It did:

* `eff_rna` changed on **5,030 of 22,138 edges (22.7 %)**; mean edge divisor **233.7 → 204.2**; max
  single change 233.7 (some edges go to 0).
* Per-object, **1,297 nodes** and **2,053 edges** moved by more than 1e-9.

### ⛔ AND THE REASON IS THE WEIGHTING — which is what fact 6 got wrong

| taper ratio | value |
|---|---|
| **position**-weighted (fact 6's convention) | **0.8738** — reproduces the recorded 0.8904 |
| **mass**-weighted (what the estimator sees) | ⭐ **0.9596** |

Fragments concentrate mid-transcript, which is exactly where the taper is inert:

| lines | slots | edge mass |
|---|---|---|
| taper inert (~1.00×) | 17,559 | ⭐ **89.17 %** |
| 0.50–0.99× | 1,837 | 7.21 % |
| 0.25–0.50× | 392 | 0.97 % |
| **< 0.25× (severe)** | 2,350 | **2.66 %** |

**A bp-weighted mean gives every terminal position equal say with every mid-transcript position. The
calibration does not.** That single substitution turns an "11.0 % over-call" into a 4.0 % geometric
effect — and then the solver turns 4.0 % into ≤ 0.0002, because composition is anchored by far more than
one divisor.

### ⭐ "+0.36 in the last node" IS TRUE, AND IRRELEVANT — both at once

Max per-node `|Δf_gdna|` = **0.3961**, so the original claim is right about the worst node. But:

| `|Δf_gdna| >` | nodes | fragments they carry | share of node mass |
|---|---|---|---|
| 0.30 | **6** | **17** | **0.0002 %** |
| 0.10 | 57 | 454 | 0.0064 % |
| 0.01 | 662 | 20,398 | 0.29 % |

Total gDNA mass moved: **0.021 %** on nodes, 0.52 % on edges. ⚠ This is the shape trap 19 warns about in
reverse — the per-object number is dramatic and the pooled number is nothing, and here the pooled number
is the one that decides. Both were measured rather than one being assumed from the other.

### ⚠ THE CONDITION UNDER WHICH THIS NULL WOULD FAIL, AND THE ONE-NUMBER TEST

The null holds **because mass is mid-transcript**. It would not hold for a 3′-biased or heavily degraded
library — **cfRNA is degraded** — or where transcripts are short relative to the fragment length. ⛔ So
this is not "reach never matters"; it is "reach does not matter when mass is mid-transcript".

⭐ **The screening test is the mass-weighted taper ratio** (0.9596 here). It is one number, needs only the
index and a payload, and costs nothing next to wiring the taper. Compute it before assuming.

### What this decides

`docs/SESSION_HANDOFF.md` §1 fact 6 corrected — the geometry stands, the "11.0 % gDNA over-call" label does not.
`CLAUDE.md` and `docs/calibration/S5_DESIGN_LOG.md` §0 updated: A7 was the headline expected gain of S5.g and it is not
there. ⚠ **The `edge_rna_reach` switch and `build_contiguous_edge_reach_arrays` are now machinery that
buys ≤ 0.0002 and costs a strand axis on `eff_rna` plus a new reduction heuristic at `bp_solver:334`.**
Recommendation on the table: delete both under "converge and delete", keep the finding and the screening
test. Owner's call.

---

