# TODO — the project's deferred work, ranked

**This is the one list.** Add here rather than starting another file, and delete an item when it lands
rather than marking it done — `LEDGER.md` is where finished work is recorded, with its gates and its
reasoning. Every item states **why it is deferred**, because "we'll get to it" is how a list stops being
read.

Sequence: `IMPLEMENTATION_PLAN.md` §0 is the live handoff and says what is happening *now*. This file says
what is queued behind it.

---

## ⛔ 2026-07-30 — every benchmark and every index was DELETED

Owner decision. Both are being rebuilt from scratch. `LEDGER.md`'s deletion entry says exactly what that
voided; the short version is that **the standing baseline and the goldens refer to nothing** and must not be
quoted, and that **no test depended on the deleted data** (grep-verified), so the suite still runs.

⭐ **It also settles the sequencing question that was open here:** build the benchmark suite before S5, or
run S5 blind? There is now no suite to run S5 against, so the suite comes first.

## The critical path

| | item | why now / why not yet |
|---|---|---|
| **1** | ⭐ **Index building** | Rebuilt from scratch. ⚠ Fix the manifest-provenance gap below while doing it — it is the cheapest moment |
| **2** | ⭐ **The new benchmark suite** | Now blocking S5. Requirements below, and they are a requirements document, not a wish list |
| **3** | **S5 — the consumers** | The substrate collapse, `build_node_geometry`, `effective_length`, `fl.py`'s pool axis. `IMPLEMENTATION_PLAN.md` §4 ranks the four by risk |
| **4** | **S6 — delete** | `ruff check` undefined-name failures are the authoritative list; goldens regenerated once |
| **5** | **A new benchmark skill** | Only once the suite exists. The old one was deleted; writing a replacement against a suite that does not exist yet would be speculative |

---

## ⭐ 2. The new benchmark suite — the requirements

⛔ **The 32-condition `ambig_dense_10mb` suite cannot judge this design.** Measured: fragment-length
variance is **zero** (every fragment exactly 200 bp), it is Poisson by construction, and its fine node set
is **row-for-row identical** to its merged region set, so it cannot see a partition change at all. It also
ranks performance hotspots backwards. Standing rule: *before running a benchmark, prove it can resolve the
axis you are changing.*

**So S5's stated gate — "the delta is measured and recorded, not required to be zero" — is unfalsifiable
without one.**

The design (§11.3) and `docs/testing/testing_plan.md` agree on the shape: **real human genes as the
backbone**, so calibration sees real fragment-length distributions and a real strand model; **cache the
whole BAM-scan intermediate** (a 24-condition run went ~13 min → ~9 s, bit-identically); then piggyback a
synthetic stress chromosome onto the cached objects. It **must** contain, or the design's own failure modes
are invisible:

* a **density step**, not just a uniform background — over a run of flat nodes a relayed message decays
  geometrically per hop, so a uniform scenario cannot distinguish "the relay works" from "the global prior
  reached it";
* **fragment-length variance** — the deleted suite had **exactly zero** (every fragment 200 bp), so nothing
  length-dependent was exercised at all, and length is now the core of the design;
* **non-Poisson counts** — ⚠ the simulator draws multinomial at fixed abundance (`sim/wgs_engine.py:473`),
  so overdispersion is **0 by construction**. This has to be built IN or the new suite inherits the same
  blind spot;
* **alternative TSS/TES inside exons** — 53.4 % of real transcript termini fall strictly inside an old
  merged region, which is the entire reason for the v8 partition;
* **ample single-stranded nodes** in any both-strand stress test — the population prior trains on them, and
  an isolated AMBIG node with no single-stranded neighbours is a starved toy, not a hard case;
* a **low-gDNA × strong-capture corner** (1–10 % gDNA) — real libraries live there, and 0/100/300 %
  conditions cannot see it;
* and it must be able to **resolve a partition change at all**: the deleted suite's fine node set was
  row-for-row identical to its merged one.

⭐ **How to evaluate it** is already written down — `docs/BENCHMARKING.md` now carries the net-fragment-flow
methodology preserved from the deleted skill, including why hard per-fragment label recovery is the wrong
target and why a byte-identical hard-label result is *no evidence* rather than *no change*.

⚠ Also outstanding and waiting on this: **the goldens moved at the index change** (`gdna_em_count` fell
16–52 %) and that movement has **never been validated against truth**. The new suite should adjudicate it;
the accumulator sequence should not.

---

## ⚠ Performance regressed at S4, and it is measured

Re-recorded on the same three real cfRNA libraries the S2 baseline used:

| | before (shipped accumulator) | after S4 |
|---|---|---|
| per-fragment deposit | ~357 ns (350–400; 3 points, visible scatter) | **410.0 ns** |
| fixed partition cost | 0.108 s | **0.348 s** — ⛔ 3.2× |

The per-fragment number sits at the edge of the old band and scatter could account for part of it; **the
fixed cost cannot be explained that way.** It is `O(partition)` work, consistent with the payload
materialising eleven arrays over 1.04 M / 1.04 M / 404 k rows against six before, plus a larger per-worker
accumulator. Paid once per scan: amortised on a deep BAM, dominant on a shallow one.

**Why deferred:** owner ruling — *"i'm not worried about the performance gate. we will make it fast
eventually."* ⚠ The non-negotiable part already holds: the deposit allocates nothing per fragment. The
obvious first look is `build_result`'s AoS→SoA transpose and the per-worker zeroing, both O(partition).

## `strandedness` is not on the payload

Design §5.2 lists it as a header field — *measure it, declare it, assert it in QC*. It is derivable from
`strand_observations` (the agreement of read-1 orientation with the splice-motif strand over spliced
fragments) and belongs with the strand model, so S4 left it out rather than half-building it. ⚠ The new
schema **cannot express** the bug it was meant to guard against — no channel is labelled *sense*, so no
label can be wrong — which is why this is a diagnostic to add rather than a hole to plug.

## Nothing keys a cached CALIBRATION result against `reach`

⭐ Surfaced by the index rebuild: `partition_hash` covers `nodes.feather`, `graph_hash` covers that plus the
junction CSR — and **neither covers `reach`**, correctly, because neither the scan nor the accumulator reads
it. Both hashes were byte-identical across a rebuild that moved ~38 % of human contiguous reaches. But reach
IS consumed by calibration, so any cache of a *calibration* output needs a third key that includes it.
Nothing caches one today; this becomes live the moment something does.

## ⛔ An index cannot be rebuilt from its own manifest

`manifest.json` records `format_version`, `rigel_version` and `mappability` — **not the FASTA, the GTF, or
the flags**. Rebuilding the human index meant guessing the source and then discovering
`--collapse-duplicate-transcripts` from a build failure (the GENCODE GTF has 141 duplicate-transcript
groups, 283 transcripts, which are mathematically unidentifiable and must be collapsed). ⚠ An artifact that
cannot be reproduced from its own provenance is one nobody can safely re-derive later — and the check that
saved this rebuild was **`nodes.feather` byte-identity**, which is the only reason a wrong source would have
been caught rather than silently swapped in.

**Fix:** record the source paths, their content hashes, and the build flags in `manifest.json`.

## Deferred deliberately, with the reason

### ⛔ `detect_chimera` is blind to two real populations, because of one gate

It considers only blocks with a **non-empty transcript set** (`constants.h:339-343`) and needs ≥ 2 mutually
disjoint transcript components, so it returns `CHIMERA_NONE` for both of:

1. **Same-reference-orientation mates.** Paired-end sequencing reads the two ends of one double-stranded
   molecule, so R1 and R2 must map to **opposite** reference strands; `build_fragment` encodes that by
   flipping R2, so a normal FR pair collapses to one `(ref, strand)` group. Two groups means the mates
   agree in reference orientation — evidence of a rearrangement, not a molecule. It becomes
   `align_strand = STRAND_AMBIGUOUS` and, since S3, is **counted** as `dropped_strand_undefined` (158 / 615
   / 246 on the three cfRNA libraries) rather than vanishing.
2. **The multi-megabase spans.** ⭐ 95 % carry a supplementary record, and intergenic blocks have empty
   transcript sets, so the gate never sees them. Same root cause, second symptom.

The vocabulary already exists: `CHIMERA_CIS_STRAND_SAME` / `CHIMERA_CIS_STRAND_DIFF`.
**Why deferred:** widening the gate reclassifies currently-accepted fragments as chimeras, which moves the
bench. That is a change to *what counts as a fragment*, not to how one is tallied — so it is its own arm
with its own before/after measurement. Owner-agreed 2026-07-29.

### Single-reference *cis* chimeras still deposit on the intergenic path

Related to the above but distinct, and raised by S3. `cr.chimera_type` is set for single-reference cis
chimeras, and the resolved path honours it (`get_is_chimeric()` → `continue`) while the **intergenic path
does not** — it has no chimera check at all. S3's adapter deliberately gates only on **multi-reference**,
which is what a `FragmentPath` cannot express, and left this alone because stopping it is the same
"what counts as a fragment" change as the item above.
⭐ Measured: the narrow and broad gates admit **identical** fragments on LBX0190 and MO_3021, so nothing is
lost by having deferred it.

### The side buffer: multimappers and `path_ambiguous`

Design §9 has the recipe — group by `(candidate set, column)`, apportion each group's integer count across
its candidates in proportion to their densities, rounding by **largest remainder**. Deterministic, integral,
single-pass, no RNG. ⚠ The assignment must be **integral**: a fractional `1/NH` split re-creates the
non-integer observable this redesign exists to delete.

⭐ **Now quantified rather than speculative**, which is what makes it worth scheduling:

| population | share |
|---|---|
| implicit splices with an undetermined path (`path_ambiguous`) | **70–76 %** of all implicit splices; 0.04–1.03 % of everything accepted |
| multimappers, as a share of **intergenic** fragments | **59.3 %** (LBX0190) / 39.5 % (MO_3021) / 20.2 % (vcap) / 16.7 % (LBX0588) |
| multimappers, intronic | ~53.6 % |

The multimapper gate thins the **one pure gDNA pool** non-randomly, toward repeats — so this is not only a
recovered-fragment count, it is a bias in the length model.

### `contradictory_sj_strand` has never fired on real data

⭐ Measured 0 on all three cfRNA libraries. Expected for STAR, which writes a consistent `XS` per record —
but it means one branch of the three-way `sj_strand` rule is untested outside the spec matrix. Worth one
look before it is trusted; a deliberately mixed-tag BAM would settle it.

### `sj.feather` duplicates junction coordinates

One row per *(junction, transcript)*, feeding `sj_map` at load (`index.py:1351`), so it cannot merge into
`edges.feather` (one row per distinct junction). But its `(ref, start, end, strand)` columns duplicate what
`(src, dst, strand)` already determines, and **nothing enforces they agree**.
**Available simplification: re-key it to `(junction_edge_id, t_index)`.** Owner's idea.

### `align_strand` → `strand` tree-wide

101 occurrences in `src/` and `tests/`, and ⚠ **seven are string keys** (`buffer.py:193,342,381` — a parquet
column in the spill path — plus `resolve.cpp:38`, `resolve_context.h:307,352`, `tests/test_buffer.py:45`).
A string key survives compilation and fails at **runtime**, on the buffer→EM path. **Why deferred:** that
path is untouched by the accumulator rework, so bundling it in would put an unrelated runtime-failure mode
inside a step whose gate is byte-identity.

### `CLAUDE.md`'s document list is wrong

It says four documents are the whole story and that "a doc path not listed above does not exist". There are
in fact **eight**: the four, plus `BENCHMARKING.md`, `MANUAL.md`, `PUBLISHING.md` and
`docs/testing/testing_plan.md` — the last of which is the owner's plan for the cached-substrate harness and
is load-bearing for item 3 above. Fix the list rather than the files.
