# Ledger

**One row per step, appended as it lands, never retroactively.** Plan: `IMPLEMENTATION_PLAN.md` §0.
Design: `ACCUMULATOR_DESIGN.md`. Deferred work: `TODO.md`.

The point of this file is **attribution**: a delta is only attributable if it is recorded against a
baseline taken from the same tree in the same session.

⛔ **The old benchmark baseline is VOID.** `r0 0.079005 / r3 0.046675` was the 32-condition
`ambig_dense_10mb` suite, deleted 2026-07-30. Do not quote it, compare against it, or try to reproduce it.

⭐ **THE FIRST BASELINE EXISTS as of S5.f (2026-07-30)** — eight chr22 pilot conditions, in the S5.f entry
below, bit-identical on re-run. ⚠ It carries A7's known **11.0 %** genome-wide gDNA over-call by
deliberate ruling, and it says nothing about dispersion (the simulator is Poisson by construction).

## The index — every entry, oldest first

⚠ **Split 2026-07-30 because this file passed 1,700 lines.** Entries up to and including the deletion
moved **verbatim** to `LEDGER_ARCHIVE.md` — relocated, never rewritten, so the append-only property that
makes attribution meaningful survives. Everything from `I1` on is below in full.

| entry | what | where |
|---|---|---|
| S1 | index: reach semantics + the junction CSR | `LEDGER_ARCHIVE.md` |
| S2 · S2.1 · S2.2 · S2.3 | the Python reference, the spec matrix, the doc corrections, the naming | `LEDGER_ARCHIVE.md` |
| S3 | the C++ wired, byte-identical to the specification | `LEDGER_ARCHIVE.md` |
| S4 | the payload, typed against the specification's vocabulary | `LEDGER_ARCHIVE.md` |
| Index rebuild | the artifacts carry the S1 reach | `LEDGER_ARCHIVE.md` |
| ⛔ Deletion | every benchmark and every index deleted; the baseline voided | `LEDGER_ARCHIVE.md` |
| **I1** | an index can be rebuilt from its own manifest | below |
| **I2** | the human index rebuilt from scratch — and it is the SAME index | below |
| **I3 + B0** | the suite backbone, and the gate that says what it can judge | below |
| **B1 · B1a · B1b** | the capture simulator made fast, exactly | below |
| **S5.0** | ⭐ S5 STOPPED and turned into a derivation; the observable question answered | below |
| **S5.a** | ⭐ `length_sum` added to every population; `density` renamed `inv_length_sum` | below |
| **S5.b** | `fl.py` re-keyed to the five pure pools; the scan cache unblocked | below |
| **S5.c** | `effective_length.py` reduced to the ONE placements formula | below |
| **S5.d** | substrate collapsed to ONE type; the chain re-keyed to nodes/edges | below |
| **S5.e** | ⭐ A7 ruled; the 18 per-face arrays became 5; the faces dissolve | below |
| **S5.f** | ⭐⭐ **`calibrate()` RUNS — the FIRST BASELINE**; the third axis exported; the pooling dissolves | below |
| **S5.f-addendum** | the κ mirror is CONSISTENT — the inference is correct, only the label is wrong | below |
| **S6** | ⭐ goldens regenerated; two S5.f defects exposed; the EM's sampling mode measured | below |
| **S5.g-1** | the per-contiguous-edge RNA reach built + gated; ⛔ one open question gates the divisor | below |
| **S5.g-2** | ⛔ **the A7 taper A/B: it moves the gDNA fraction by ≤0.0002.** Fact 6's 11.0 % was bp-weighted geometry | below |
| **B2** | the pilot suite, and B0's verdict on it | below |
| **B3** | the scan cache — scan once, calibrate many times | below |
| **C2.6** | ⭐ **gap introns cut on EVERY fragment; the anchor's impossible tail is 85 % gone and the residual is LOCALISED to D3.** ⚠ D1 measured to cost more than it buys | below |

---


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
**chr1**; genome-wide it is **2,970**. And `CARRY_FORWARD.md` §1 fact 2's denominator, *"232,451 of
435,291"*, re-derives as 232,451 of **435,107** terminus-flagged contiguous seams — the numerator and the
53.4 % are exact, the denominator is 184 out and depends on how a terminus is counted.

Also new, from the census and load-bearing for the deposit: **591,442 nodes (56.7 %) are shorter than
200 bp**, i.e. shorter than one RNA fragment — `CARRY_FORWARD.md` §0 C5's unpriced correlation risk applies
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
(`CARRY_FORWARD.md` §1 fact 22), so a genome-scale index taxes every future iteration at ~66 s a run.

⚠ **The ERCC controls are not filler.** `CARRY_FORWARD.md` §3 trap 20: a single-reference synthetic index
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

`CARRY_FORWARD.md` §3 trap 15 made executable, and **written before the suite exists** so it cannot be
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
Python loop sharing no helper with the implementation (`CARRY_FORWARD.md` §3 trap 1). They passed against
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
change would have mis-ranked the work — `CARRY_FORWARD.md` §3 trap 16, in a new costume.

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
quantity, which is the shape of `CARRY_FORWARD.md` §3 trap 27.

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
| `scratchpad/` | 31 tracked scripts, 8.9 MB | gone | one-off prototypes for a design that is now implemented and tested. Their *results* are `CARRY_FORWARD.md` §1; the citations now point at git history |
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

* **`LEDGER.md` split at 1,715 lines.** Entries up to the deletion moved **verbatim** to
  `LEDGER_ARCHIVE.md` — relocated, never rewritten, so the append-only property that makes attribution
  meaningful survives. `LEDGER.md` keeps an index of *every* entry plus the full text from `I1` on.
* **`CARRY_FORWARD.md` §4 deleted** (441 → 310 lines). It was a pre-settlement proposal list whose six
  worst bullets already carried `⛔ SUPERSEDED` markers; a proposal list that contradicts the shipped
  deposit rule is a hazard, not history. ⭐ **The rest of the file is NOT dead** — §3's traps caught real
  defects repeatedly this session and §2 is the derivation reference; the header now says so and gives a
  read order.
* **`BENCHMARK_SUITE.md` is new** — what the suite is, how to build it, and what it can and cannot judge.
* **`TODO.md` rewritten** against the actual state; **`IMPLEMENTATION_PLAN.md` §0** replaced with the
  session handoff; **`CLAUDE.md`** doc list, test counts and tooling table corrected.

### Gates

* Full suite **1321 passed / 266 failed / 2 xfailed / 15 errors** — **identical to before the cleanup**.
* `ruff check src/ tests/ scripts/` and `ruff format --check src/ tests/` clean.
* ⚠ Nothing committed. `git checkout` recovers any file the owner disagrees with.


---

## S5.0 — S5 stopped, and the observable question answered (2026-07-30)

**Owner ruling.** S5 was written as "rewire the four consumers". Asked which observables the rewiring
should consume, the design had no answer — so wiring it would have frozen a schema nobody had checked.
S5 became a derivation first. `docs/S5_DESIGN_LOG.md` is the running log;
`docs/NODE_DENSITY_DERIVATION.md` is the write-up.

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

⚠ **`CARRY_FORWARD.md` §1 fact 5 does not reproduce** and needs correcting or retiring. It was measured
at a 4× mean separation; over the full grid it is false at every node ≥ 250 bp and at the edge.

**Rulings taken** (owner, all recorded in `S5_DESIGN_LOG.md` §4): every population stores
`(count, inv_length_sum, length_sum)`, nodes **and** edges; the node deposit weight is unchanged; `Σ1/A`
is not stored; `CalibrationResult`'s left/right split becomes per-edge.

## S5.a — `length_sum` on every population, and `density` renamed (2026-07-30)

Two sub-steps, deliberately separated so each has its own gate.

**S5.a-1, the rename — behaviour-neutral, and proven so.** `density` → `inv_length_sum`,
`DENSITY_SCALE` → `INV_LENGTH_SCALE`, `density_quantum` → `inv_length_quantum`, across 9 files.
⭐ `Σ1/L` **is** an exact density at an edge and is **not** one at a node — the design's own §6 says so —
so the old name put one word on two concepts (`CARRY_FORWARD.md` §3 trap 27). Gate: the full suite
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
`CARRY_FORWARD.md` §3 trap 25 exactly.

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
reproduces `ACCUMULATOR_DESIGN.md` §8's table.** (A *test*, not a design input.)

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
| ⭐ independent cross-check | ✅ reproduces `CARRY_FORWARD.md` §2's published taper table: **160.095 against 160.1** at R=200, and 88.2/19.8/198.5 against 87.8/19.6/199.0 |
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
short" (`CARRY_FORWARD.md` §3 trap 23).

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
turned out to be (`S5_DESIGN_LOG.md` §1 A7, §4). ⭐ The decision reduces to **one call site**, and that
is what made it rulable: gDNA at any edge is `taper_g = 1` by physics and the two contained frames take
no reach argument, so the only question was the RNA component at a *contiguous* edge.

| edge kind | component | reach | why |
|---|---|---|---|
| contiguous | gDNA | `UNBOUNDED_REACH` | its template is the chromosome — physics, not a ruling |
| **contiguous** | **RNA** | ⭐ **`UNBOUNDED_REACH`** | keeps S5.e varying ONE thing; A7 proper is **S5.g**, the only ordering where it gets an end-to-end A/B against S5.f's first baseline |
| junction | RNA | ⭐ **the real exonic per-strand reach** | a BRAND-NEW population — the predecessor had no junction divisor at all — so wiring it regresses nothing, and it keeps the one code path exercised rather than a dead branch |

⚠ **The first baseline will therefore carry a known bias and S5.f's entry must say so**:
`CARRY_FORWARD.md` §1 fact 6 — ignoring the taper over-calls gDNA by **11.0 %** genome-wide and by
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
  boundary axis became a plain per-contiguous-edge array (`S5_DESIGN_LOG.md` §4's ruling), so the
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
  convenience (mature RNA never crosses an exon↔intron seam — 0 of 1,146, `CARRY_FORWARD.md` §1
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
  (`CARRY_FORWARD.md` §1 fact 18) — so this baseline can say nothing about dispersion, by design.
* **f_gdna lands within 6 % of truth on 3 of 4** gdna100 conditions.
* ⚠ **The zero-gDNA false positive is 7.9 % in the worst corner** (unstranded, capture off) and falls to
  0.16 % when stranded and captured. That ordering is what §2 predicts — at κ = ½ the strand channel
  carries *exactly zero* information about composition, so the solver has only count and density left.
  ⛔ Quote it **only alongside the gdna100 column**: a zero-gDNA guard is ONE-SIDED (trap 19) and any
  change that merely lowers gDNA scores better on it.

### ⚠ TWO ANOMALIES THIS BASELINE SURFACES — neither is S5.f's, both are now measurable

1. ⛔ **The fitted κ is `1 − truth`.** A library simulated at `strand_specificity = 0.99` calibrates to
   **κ = 0.0101**; at 0.50 it calibrates to 0.4990–0.5002 (where the flip is invisible). This is
   `CARRY_FORWARD.md` §0 **C4** — "sense fraction is 0.002–0.012 on all four real cfRNA libraries
   (nearly fully antisense) — possibly a read-orientation convention bug" — now answered against ground
   truth for the first time: it is a **convention flip, not biology**. S5.f does not touch
   `fit_strand_balance`, so this is pre-existing and merely became visible. **It needs its own step.**
2. ⚠ **`gdna100 ss0.50 capture_on` reads f_gdna 0.3754, 25 % low**, and its RNA total (7.30 M) is 21–36 %
   above every other condition's. The other three gdna100 conditions are within 6 %. Unexplained; it is
   the both-worst-case corner (no strand information *and* capture), so it is the natural first target.

### ⚠ THE BASELINE CARRIES A KNOWN BIAS, AND IT IS NAMED

**A7 is OFF at contiguous edges by owner ruling** (`S5_DESIGN_LOG.md` §1 A7): the RNA half of an unspliced
crossing takes `UNBOUNDED_REACH` rather than its transcript's real remaining length. Measured cost,
already on record (`CARRY_FORWARD.md` §1 fact 6): an **11.0 % genome-wide gDNA over-call**, contiguous
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
hid an **exact factor of 2 for months** (`CARRY_FORWARD.md` §3 trap 2). It cannot recur: there is no face
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
  (`CARRY_FORWARD.md` §1 fact 4) — now as a ratio of SUMS on both axes rather than a mean of per-boundary
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
  the skip test. Same family as `CARRY_FORWARD.md` §3 trap 25 — a cache key that does not cover the
  artifact it caches.

---

## S5.f-addendum — the κ mirror is CONSISTENT, so the inference is correct (2026-07-30)

**Not a step; a measurement made against S5.f's baseline, recorded here rather than by editing S5.f's
entry.** The S5.f entry reported the fitted κ as `1 − truth` and said the flip "needs its own step"
without establishing how severe it is. It is now established, and it is **much less severe than it
looks**.

**The question.** κ enters the only intrinsic gDNA/RNA signal — `p = ½·f_g + κ·(1−f_g)`
(`CARRY_FORWARD.md` §2). If κ were mirrored while the per-node sense counts were not, the strand
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
*while* `f_gdna` stays at 0.0030); `CARRY_FORWARD.md` §0 C4 (answered) and §1 fact 17 (its κ column reads
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
(`IMPLEMENTATION_PLAN.md` §469). That list is **empty**: `ruff check --select F821` passes, there is not
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
about one contract, disagreeing (`CARRY_FORWARD.md` §3 trap 27). ⭐ Resolved by **wiring it, not deleting
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
family as `CARRY_FORWARD.md` §1 fact 11. Against `RTOL = 1e-6` that is **nine orders of magnitude of
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
negative control is a **one-sided** metric (`CARRY_FORWARD.md` §3 trap 19: in a library with none of X,
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
* `docs/TODO.md` (§6, §7 + the ranked critical path) · `docs/CARRY_FORWARD.md` (§0 C4, §1 fact 17).

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

* **PER STRAND and per SIDE** — `CARRY_FORWARD.md` §2: reach is "maximised over transcripts
  independently per side AND per strand". This is what forces `eff_rna` to gain a strand axis; a single
  averaged reach describes neither transcript.
* **GENOMIC, not exonic** — unlike a junction's. A junction is used only by a spliced molecule so what
  remains either side is exonic; a contiguous line is also crossed by *nascent* RNA, which is genomic.
  Taking the exonic reach here would declare an intronic nascent fragment impossible.
* ⚠ **A reach of 0 is the ANSWER, not a missing value.** No template of that strand at that line ⇒ zero
  opportunity ⇒ divisor 0 ⇒ the consumer emits nothing (`CARRY_FORWARD.md` §3 trap 23). Measured on the
  chr22 pilot index: **40.6 %** of contiguous edges have no POS template, **42.9 %** no NEG.

### ⭐ It independently reproduces fact 6's magnitude

On the chr22 pilot index with the suite's own fragment-length distribution (mean 206, sd 98), over the
lines where the strand has any template at all:

| | unbounded divisor | tapered mean | ratio |
|---|---|---|---|
| POS | 216.6 | 182.5 | **0.8425** |
| NEG | 216.6 | 186.1 | **0.8592** |

`CARRY_FORWARD.md` §1 fact 6 records a bp-weighted **0.8904** over 3,000 genes, and notes contiguous
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

`CARRY_FORWARD.md` §1 fact 6 corrected — the geometry stands, the "11.0 % gDNA over-call" label does not.
`CLAUDE.md` and `S5_DESIGN_LOG.md` §0 updated: A7 was the headline expected gain of S5.g and it is not
there. ⚠ **The `edge_rna_reach` switch and `build_contiguous_edge_reach_arrays` are now machinery that
buys ≤ 0.0002 and costs a strand axis on `eff_rna` plus a new reduction heuristic at `bp_solver:334`.**
Recommendation on the table: delete both under "converge and delete", keep the finding and the screening
test. Owner's call.

---

## P0 — the composition-evidence census: HOW MUCH MASS REACHES THE SOLVER BLIND (2026-07-31)

    Plan: `docs/SOLVER_OBSERVABLES_PLAN.md` §4   ·   Tool: `scripts/design/composition_evidence_census.py`
    Prediction recorded BEFORE the run: the ss0.50 conditions carry materially more no-evidence mass
    than the ss0.99 ones, concentrated on AMBIG slots.

⭐ **CONFIRMED, and the finding is larger than the prediction.** P0 exists to falsify P2's premise before
P2 is built. It did not falsify it.

**What is measured.** The mass-weighted share of unspliced fragment mass on chain slots with
`tau_lam <= 1e-9` and **no structural lock** — i.e. slots whose own composition belief carries
`own_precision = 0`, so their gDNA/RNA split is decided entirely by neighbour messages and the
population prior. ⚠ Structurally locked (intergenic pure-gDNA) nodes are EXCLUDED: they are composition
*certain*, not uninformed, and lumping the two would report a correct answer as a failure.

⚠ **`tau_lam = 0` is not "the strand likelihood was ignored".** The local ψ solve still sets the
*location* `f_g` from strand; what is zero is the *precision*, so the node emits nothing and its own
answer carries no weight in the fusion. The census measures blindness of the message system, not of the
per-node read-out.

### The result — 8 chr22 pilot conditions, no production code changed

| condition | κ | **no-ev** | node | edge | AMBIG | 1-str | lock |
|---|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 0.5000 | **49.4 %** | 28.5 | 20.9 | 12.4 | 36.7 | 28.7 % |
| gdna100 ss0.50 capture_on | 0.4990 | **98.2 %** | 58.7 | 39.5 | 29.0 | 68.1 | 1.0 % |
| gdna100 ss0.99 capture_off | 0.0101 | **12.7 %** | 6.2 | 6.5 | 12.4 | 0.0 | 28.7 % |
| gdna100 ss0.99 capture_on | 0.0101 | **30.0 %** | 16.3 | 13.7 | 29.0 | 0.0 | 1.0 % |
| **none** × 4 (both κ, both capture) | — | **100.0 %** | 61.6–67.2 | 32.8–38.4 | 26.6–40.1 | 59.9–73.4 | 0.0 % |

Unstranded 86.9 % against stranded 60.7 % overall — but the pooled ratio understates it, because the
`none` arm is saturated at 100 % on both. **On the gdna100 arm, where the comparison is clean, it is
49.4 % vs 12.7 % and 98.2 % vs 30.0 % — 3.3× and 3.9×.**

### ⭐ Four findings, three of them not predicted

1. **The strand channel's on/off is exact, and visible in one column.** `no-evidence | single-strand` is
   **0.0 %** on both κ = 0.0101 conditions and **63.6 % / 98.9 %** on the κ = ½ ones. Every single-strand
   slot in a stranded library has own evidence; most have none in an unstranded one. This is
   `CARRY_FORWARD.md` §2's `I(f_g) ∝ (2κ−1)²` read off production output for the first time.
2. ⭐ **AMBIG mass is blind in EVERY condition, stranded or not**: `no-evidence | AMBIG` is **93.3 % –
   100.0 %** on all eight. That is the Schur gate (`node_init`, approach E) doing exactly what it was
   designed to do — a both-strand node's tilt is a free nuisance, so strand cancels out of `f_g` — and it
   means **13.3 – 40.1 % of library mass has never had own composition evidence in any library Rigel has
   ever run.** ⚠ This is not an unstranded-data problem. It is every both-strand locus, always.
3. ⭐ **A zero-gDNA library is 100 % blind at any κ.** `strand_evidence`'s derived noise floor is
   `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`, so `N_gdna = 0 ⇒ σ²_d → ∞ ⇒ disc = 0`. Documented and
   deliberate, but the consequence had not been quantified: on the `none` arm **nothing** gives any node
   own evidence, and `none ss0.50 capture_off` nevertheless reports `f_gdna = 0.0793` against a truth of
   exactly 0 — 443,277 phantom fragments, produced entirely by messages and the prior.
4. ⭐ **Capture roughly doubles the blindness**, by removing the anchor rather than by adding noise:
   structurally-locked mass collapses **28.7 % → 1.0 %** because intergenic nodes are depleted, and
   no-evidence goes 49.4 % → 98.2 % (ss0.50) and 12.7 % → 30.0 % (ss0.99). The pure-gDNA anchor the
   prior-free pass leans on is a *capture-off* asset.

### The falsification — the census fires, and lands on the natural experiment

⛔ A census that does not move when the thing it claims to measure is switched off measures nothing
(`falsification_needs_perturbation`). Re-run `gdna100 ss0.99 capture_off` with **κ = 0.5 injected and
nothing else changed** (`--inject-kappa 0.5`, via `dataclasses.replace` on the condition's own fitted
`InjectedCalibrationPriors`):

| | no-evidence | no-ev \| single-strand |
|---|---|---|
| as fitted, κ = 0.0101 | 12.7 % | 0.0 % |
| **κ → 0.5 injected** | **49.4 %** | **63.6 %** |
| `gdna100 ss0.50 capture_off`, the natural experiment | **49.4 %** | **63.6 %** |

⭐ **The injected arm reproduces the natural experiment to the last digit** — same payload, one variable,
landing on the independently-simulated condition's value. The census tracks the strand channel and
nothing else.

### What this decides

* **P2 is justified**: the length likelihood is the only proposed source that is independent of κ, of
  `N_gdna`, and of the single-strand/AMBIG gate. Finding 2 is the strongest case for it — the Schur
  argument that silences strand on AMBIG nodes does not apply to a channel that has no dependence on the
  tilt `θ` at all.
* ⚠ **P2's gate must not be scored on the `none` arm alone.** All four zero-gDNA conditions are saturated
  at 100 % blind, so any change scores "better" there (trap 19, one-sidedness). The clean comparison is
  the gdna100 arm, and `gdna100 ss0.50 capture_on` (98.2 % blind, `f_gdna` 0.3754 against ~0.50 — the
  row S5.f recorded as "unexplained") is the sharpest single target in the suite.
* **`SOLVER_OBSERVABLES_PLAN.md` §4 updated** with the result. P1 (the units fix) is unaffected by this
  measurement and remains next.

---

## P1 — the EM prior is a FRAGMENT COUNT, not an object-incidence sum (2026-07-31)

    Plan: `docs/SOLVER_OBSERVABLES_PLAN.md` §5   ·   Unit gates: `tests/calibration/test_prior_units.py`
    End-to-end gate: `scripts/design/prior_units_check.py`

⭐ **A UNITS BUG, WITH A PROOF — not a modelling change.** ``assemble_priors`` handed the EM two additive
pseudocounts that are added straight to its own fragment counts (``G = n_gdna + a_g``,
``em_solver.cpp:apply_grouped_prior_update``), but built them by **summing per-object masses**. A fragment
deposits on ``max(K, 1)`` objects, so for a partition of spacing ``s``::

    incidences(w) = max( 1 , (w-1)/s )

Counts are conserved exactly where every node is longer than every fragment and become a
**length-weighted** count where they are not — and the weight is the fragment's own length, so it does
not cancel between two components with different mean lengths.

**The fix, in one line of model:** density is intensive, so pool it as a ratio of sums and integrate it
over the span. ``rho_c = Σ share·mass_c / Σ share·A_c``, then ``prior_c = rho_c · span_bp`` — the **same
genomic span for both components**, so the ratio carries no length tilt. ``A_r`` did not exist:
``CalibrationResult`` gained ``rna_node_eff_len`` / ``rna_edge_eff_len``, **projected** off
``NodeGeometry.eff_rna`` (never recomputed — trap 27) exactly as their gDNA twins are.

### The gates, written first and verified failing

⭐ The failure pattern was **self-confirming**, which is why it was worth writing them deterministically
rather than end to end:

| gate | before | after |
|---|---|---|
| partition invariance, 1200 bp re-tiled | ✅ passes at 1200 and 400 bp nodes, ⛔ **fails at 100 bp** | passes at all four tilings |
| …the RNA prior across that re-tiling | **50.05 → 109.45, a 2.19× swing** with the library unchanged | invariant |
| prior == true fragment count | 34.53 against 36.0 | exact |
| g:r flat in ``mu_g/mu_r`` (swept 0.25×–2×, both directions) | ⛔ fails everywhere **except exactly ``mu_g = mu_r``** | flat |

Invariance holding at 1200/400 bp and breaking at 100 bp is precisely the ``max(1, (w-1)/s)`` crossover;
the ratio test passing at its own null (equal means ⇒ no tilt) identifies the mechanism as the length
tilt and nothing else. gDNA is *partition-invariant in that fixture* because its 50 bp fragments fit
inside every node — the swing is entirely on the component whose fragments do not.

### ⛔ Perturbation P1e found a real defect, and the test was the weak one

Five perturbations. Four were caught (RNA÷gDNA opportunity, dropped ``span``, reverted raw sums, mean of
ratios). **P1e — flooring the divisor to ``1e-9`` instead of testing ``support > 0`` — left all 15 tests
green**, because the only zero-support fixture also had zero mass.

⚠ **``mass > 0`` with ``support == 0`` is an ordinary configuration.** ``contained_eff_length`` is exactly
0 wherever a node is shorter than that component's shortest fragment; measured on the chr22 pilot index
against its own measured pure pools, that is **21.7 % of nodes for RNA and 18.7 % for gDNA**
(`CARRY_FORWARD.md` §1 fact 9 records 12.4 % genome-wide). The solver can still put mass there — ``f_g``
is an inference, not a fact. The floor would turn one stray fragment on a 40 bp node into a density of
~1e9 (trap 23, the mechanism that once seeded false gDNA into neighbouring exons).

⭐ **And behind the weak test was a defect in the fix itself**: leaving that mass in the pooled numerator
while it contributes nothing to the denominator inflates ``rho`` with no exposure to pay for it. **Both
sides of a pooled rate, or neither** — ``_mass_where_there_is_opportunity`` now drops it from both, the
same contract ``node_geometry._rate`` already keeps on the solver side. New test
``test_mass_on_a_zero_opportunity_object_is_dropped_from_BOTH_sides``; P1e now fails it.

### ⛔ T3 caught a hard violation: the old prior EXCEEDED the whole library

`prior_units_check.py` computes both arms from ONE ``CalibrationResult`` — no re-scan, no second code
path, no simulation noise:

| condition | fragments | A: raw sum | **B: P1** | A/frag | **B/frag** | A f_gdna | **B f_gdna** |
|---|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 9,634,502 | 6,386,844 | 5,660,367 | 0.663 | 0.588 | 0.4390 | 0.4616 |
| gdna100 ss0.50 capture_on | 9,775,761 | 10,466,755 | 9,190,222 | ⛔ **1.071** | **0.940** | 0.4085 | 0.4285 |
| gdna100 ss0.99 capture_off | 9,633,129 | 6,382,577 | 5,650,154 | 0.663 | 0.587 | 0.3863 | 0.4080 |
| gdna100 ss0.99 capture_on | 9,774,382 | 10,457,211 | 9,268,515 | ⛔ **1.070** | **0.948** | 0.5360 | 0.5556 |
| none ss0.50 capture_off | 4,636,413 | 3,799,206 | 3,247,412 | 0.819 | 0.700 | 0.1167 | 0.1196 |
| none ss0.50 capture_on | 4,791,669 | 4,459,761 | 3,784,949 | 0.931 | 0.790 | 0.0126 | 0.0130 |
| none ss0.99 capture_off | 4,636,640 | 3,797,474 | 3,239,915 | 0.819 | 0.699 | 0.0044 | 0.0045 |
| none ss0.99 capture_on | 4,792,251 | 4,460,570 | 3,784,966 | 0.931 | 0.790 | 0.0020 | 0.0021 |

⭐ **On both `capture_on` gdna100 conditions the pre-P1 prior totalled MORE pseudocounts than the library
has accepted fragments** (1.071× and 1.070×) — for a quantity that arbitrates only the *unspliced*
subset. That is structurally impossible and nothing was checking it. Every arm-B row is now below 1.

⚠ **The f_gdna column is the LOCUS-PROJECTED prior ratio, not `LEDGER`'s library-wide scalar** (intergenic
nodes are dropped by the projection, and they are gDNA-rich — hence 0.4390 here against the S5.f table's
0.4998 on the same condition). The two are not comparable and must not be quoted against each other.

### ⚠ The honest scoring, including where it gets worse

The correction moves gDNA **up** on every condition, as §5.4 predicted before the run. On the zero-gDNA
arm — truth exactly 0, and the one arm whose truth survives the projection — that is **slightly away from
truth**: 0.1167 → 0.1196 on `none ss0.50 capture_off`, and +0.0001 to +0.0004 on the other three.

⛔ **This is trap 19 in reverse and must not be read as a regression.** On a library with no gDNA, any
bias that under-calls gDNA scores better; the raw sum's under-call was accidentally flattering exactly
that arm. P1 is a units error with a proof, not a knob tuned to a scenario, and the residual 0.1196 on a
zero-gDNA library is **P2's target** — it is the unstranded condition P0 measured at **100 % blind**.

### Suite

**1757 passed / 22 failed**, against HEAD's **1758 / 1** re-measured in the same session by `git stash`
(trap 17). Accounting is exact: `1758 − 21 golden + 20 new = 1757`. The 21 are `test_golden_output.py`
moving numerically, which is expected and deliberate — ⚠ **the goldens regenerate ONCE, after P2.** The
1 pre-existing failure (`test_nrna_double_counting`, `TODO.md` §7) is byte-identical.

`tests/calibration/` **569 passed** (549 baseline + 16 new + 4 new schema parametrizations, the latter
because `rna_*_eff_len` were added to the shape/sign/dtype guard lists rather than left ungated).

### ⚠ Unrelated finding: commit 49e9b456 is incomplete

`S6 + S5.g-1` committed `splice_graph.build_contiguous_edge_reach_arrays` and `test_edge_reach.py` but
**left `node_geometry.build_node_geometry`'s `edge_rna_reach` parameter uncommitted**. At HEAD the
builder therefore has no consumer — the A7 switch exists only in the working tree. Not caused by this
work and no data was lost; flagged because a builder whose only caller is uncommitted reads as dead code
to the next reader.

---

## P2 — the fragment-length likelihood: BUILT, MEASURED, and ⛔ **GATED OFF — it is blocked upstream** (2026-07-31)

    Plan: `docs/SOLVER_OBSERVABLES_PLAN.md` §6   ·   Module: `rigel.calibration.length_likelihood`
    Gates: `tests/calibration/test_length_likelihood.py` (43)   ·   A/B: `scripts/design/length_likelihood_ab.py`
    Switch: `CalibrationConfig.length_likelihood`, **default False** — the False arm is byte-identical to P1

⭐ **THE MECHANISM WORKS, EXACTLY AS DESIGNED.** P0's headline was that 13.3–40.1 % of library mass
reaches the solver with no own composition evidence in every condition, and **100 %** on an unstranded or
zero-gDNA one. Turning this channel on:

| no-evidence mass share | OFF | **ON** |
|---|---|---|
| gdna100 ss0.50 capture_off | 49.4 % | **0.0 %** |
| gdna100 ss0.50 capture_on | 98.2 % | **3.9 %** |
| gdna100 ss0.99 capture_on | 30.0 % | **1.1 %** |
| **none** × 4 | 100.0 % | **0.0 %** |

The blind mass is gone. The Schur argument that silences strand on AMBIG nodes does not apply to a term
with no ``θ`` dependence, and the measurement confirms it.

⛔ **AND THE ANSWER IT GIVES IS WRONG — 15× worse on the only arm with an unambiguous truth.**

| condition | truth | OFF | **ON** |
|---|---|---|---|
| none ss0.99 capture_off | **0.0000** | 0.0030 | **0.1269** |
| none ss0.99 capture_on | **0.0000** | 0.0016 | **0.2664** |
| none ss0.50 capture_off | **0.0000** | 0.0793 | **0.4689** |
| none ss0.50 capture_on | **0.0000** | 0.0101 | **0.6156** |
| gdna100 ss0.50 capture_off | ~0.52 | 0.4998 | 0.6064 (overshoots) |
| gdna100 ss0.99 capture_on | ~0.52 | 0.4897 | 0.4264 (worse) |

Zero-gDNA arm mean ``|f_gdna|``: **0.0235 → 0.3695**.

### ⭐ THE CAUSE IS UPSTREAM, AND IT IS NOT IN THE LIKELIHOOD

⚠ **CORRECTION, 2026-07-31, same session — the first diagnosis below was half wrong and the owner caught
it.** "The pmf is FABRICATED" is not what happens. ``build_fl_models`` EB-shrinks each pool toward the
global histogram, and at ``n_gdna = 0`` it shrinks **all the way**: verified, ``gdna_pmf == global_pmf``
byte-identically on both zero-gDNA conditions. That is exactly the intended behaviour, and the owner's
follow-on is also right — **identical pmfs for the two components make this channel exactly inert**, which
is proven byte-identically end to end (`test_identical_components_are_EXACTLY_inert_not_nearly_inert`).

⭐ **So the real defect is narrower and sharper: on a zero-gDNA library ``global`` and ``rna_fl_pmf``
describe the SAME population — every fragment is RNA — and they disagree.** Measured on
`none ss0.99 capture_off`: global **210.1 ± 85.5** against the RNA pool **234.5 ± 146.2**, i.e. **+11.6 %
in the mean and +71.1 % in the sd**, against `S5_DESIGN_LOG.md` §3.6's independently-predicted +14 % / +50 %
from junction opportunity. The length likelihood then reads that disagreement as composition, correctly —
there is no defect in the likelihood at all. ⛔ **The channel is a comparator, and it is being handed two
rulers of different length.** The original text follows.

**The gDNA fragment-length pmf reads a fallback on a zero-gDNA library.** Measured:

| condition | gDNA pool n | gdna_fl_pmf mean ± sd | rna_fl_pmf mean ± sd |
|---|---|---|---|
| gdna100 ss0.50 capture_off | 4,526,536 | 194.4 ± 97.4 | 234.9 ± 146.8 |
| **none** ss0.50 capture_off | **0** | **210.1 ± 85.5** | 234.7 ± 146.7 |
| **none** ss0.50 capture_on | **0** | **211.0 ± 85.9** | 242.3 ± 148.0 |

⛔ **Zero observations, and `build_fl_models` still returns a confident pmf** — a fallback that happens to
differ from the RNA pmf, so the likelihood dutifully discriminates against a component that does not
exist. **Every other consumer uses `gdna_fl_pmf` as a DIVISOR, where a wrong pmf is a scale error. This is
the first consumer that uses it as a DISCRIMINANT, where a fabricated pmf manufactures composition out of
nothing.**

⭐ **`strand_evidence` already guards exactly this case and the length channel does not.** Its derived
noise floor is `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`, so `N_gdna = 0 ⇒ disc = 0` — the strand
channel refuses to speak when there is no gDNA to calibrate against. The length channel needs the same
principle, and it is a **derivation**, not a flag: the discriminant is the separation between two fitted
pmfs, and it must be gated by the sampling uncertainty of that separation.

⚠ **And a second upstream defect, which the gdna100 overshoot points at.** `ACCUMULATOR_DESIGN.md` §8.1(b)
records that **every pool histogram is length-biased by its own opportunity** and must be divided by it
before use; `S5_DESIGN_LOG.md` §3.6 measured the effect on this very suite (the RNA pool reads 234.9
against a configured 206.1 — **+14 % in the mean, +50 % in the sd** — and the gDNA intergenic pool reads
195.4 against 156.5, *unexplained*). That correction is listed under §4 **"Not yet decided"**. So the two
pmfs are already tilted, by different and partly unknown opportunities, and `length_likelihood` tilts them
a second time by ``A(w)``. ⛔ **§8.1(b) is no longer a tidy-up. It BLOCKS P2.**

### ⭐ Three defects found by the falsification discipline, two of them in my own work

| | found by | what |
|---|---|---|
| **P2d** | perturbation | deleting the ``−½ log det`` term left all 41 tests green. Its grid-variation is ``O(1)`` while the quadratic's is ``O(N)``, so it is negligible at depth and **decisive at the median node**: the peak moves **0.32 at N=1**, 0.05 at N=5, 0.004 at N=50. ⚠ It also documents the model's own limit — the Gaussian is asymptotic in ``N`` and at ``N=1`` is not a trustworthy *location*, because the true single-draw mixture likelihood is bimodal |
| **the equal-pmf null** | A/B | forcing ``rna_fl_pmf = gdna_fl_pmf`` must be byte-identical and was not (0.5008 vs 0.4924) |
| **⭐ the grid-width leak** | chasing that null | rows that are flat *to 1e-11* rather than to 0 pass `density_factor_precision`'s ``ptp > 1e-12`` gate; the near-uniform posterior then returns ``tau = 1/Var(uniform over λ) = 0.029016`` — **the grid's own width sold as composition evidence**, which is exactly what that function's docstring says must never happen. Measured: **689 slots, max tau 0.02902**. Fixed by gating on the MOMENTS structurally, so the null is now *exactly* zero, and pinned by an **exact-equality** assertion (``assert_allclose`` passes at 1e-11 and would miss it) |

⚠ **The grid-width leak is latent in `density_factor_precision` itself, not only in this caller.** Any
factor that is nearly-but-not-exactly flat gets 0.029 of free precision from the grid. The intron factory
is safe today only because its inactive rows are written as exact zeros. A threshold-free fix exists —
report ``max(0, 1/Var_post − 1/Var_uniform)``, the precision the factor ADDS over the grid's own width —
but it changes a shipped, default-on feature and needs its own A/B. **Recorded, not smuggled in.**

### What landed, and what did not

**Landed** (behind `length_likelihood=False`, byte-identical when off): `calibration/length_likelihood.py`
(the tilted moments in closed form — ``O(n_nodes)`` cumulative sums, no ``n_nodes × max_len`` array, which
would be 8 GB at human scale); the two length channels gathered onto `NodeGeometry`; a `length_loglik`
argument threaded through `simplex_logodds`'s 1-D and 2-D solves, `node_init` (as **ungated** source 4) and
`bp_solver`'s local *and* final solve; `_build_length_loglik` in `calibrate`.

**Verified**: U1 against **exact enumeration over integer start positions** (a literal loop, not the
module's own closed form — trap 1) over 4 pmfs × 6 node lengths, both frames; ``q12 ≡ 1`` at a node as a
structural identity of the ``1/L`` deposit weight; ``moments.eff`` byte-identical to
`contained_eff_length`/`crossing_eff_length` (trap 27); 5 perturbations, all now caught.

**Suite**: **1800 passed / 22 failed** — the same 21 golden moves from P1 plus the 1 pre-existing
`TODO.md` §7 failure; nothing new. `tests/calibration/` 654.

### ⛔ THE VERDICT, AND IT IS NOT MINE TO MAKE

P2 is **implemented and inert**. It must not be switched on until the fragment-length pools are trustworthy
as a *discriminant*:

1. **`build_fl_models` must not fabricate a pmf from an empty pool.** A zero-observation pool has no
   distribution, and the honest output is "no gDNA length model", which must make the length channel inert
   — the same statement `strand_evidence`'s `1/N_gdna` already makes for strand.
2. **`ACCUMULATOR_DESIGN.md` §8.1(b)'s per-pool opportunity correction must land** — divide each pooled
   histogram by its own opportunity before normalising. Listed as undecided since S5.b; it is now the
   blocker.
3. Only then re-run this A/B. The before-picture (P0's blindness census) and the after-picture (this
   entry's no-evidence collapse to 0.0 %) are both recorded, so the re-run is a two-command comparison.

⚠ **Do not "fix" this by damping the channel.** The channel is doing what it was built to do; it is being
fed a length model that does not describe the library. Damping would hide the upstream defect behind a
tuned constant, which is the failure mode `CARRY_FORWARD.md` §3 trap 12 records three times over.

---

## C0 — ⭐ the accumulator's `L` is PROVEN, and it is fit to be the gold standard (2026-07-31)

    Audit: `docs/FRAGMENT_LENGTH_AUDIT.md`   ·   Gate: `tests/native/test_fragment_length_proof.py`
    Owner precondition, 2026-07-31: *"the new fragment length computation is the newest implementation
    and I'm not sure how rigorously it has been tested in all cases; we need to prove that very
    carefully if it's going to become the gold standard."*

⛔ **The precondition was right, and the existing coverage did not meet it.** `test_accumulator_spec`
pinned **six hand-picked malformed intron lists**, chosen by the same author as the code they check —
the coverage pattern that finds what the author thought of and nothing else.

### The oracle is a different ALGORITHM, not a different spelling

`CARRY_FORWARD.md` §3 trap 1. The oracle is **integer set arithmetic**::

    covered = set(range(start, end))  −  ∪ set(range(a, b)) for each intron

No sorting, no merging, no cursor walk, no ``searchsorted``. Every malformed case the production path
handles by explicit logic — overlapping, nested, abutting, duplicated, zero-length and out-of-range
introns — the oracle handles **by construction**, because set subtraction is idempotent and order-free.

⭐ **And all four deposit populations fall out of that one set**, which is what makes this a proof of the
**geometry** and not only of ``L``:

| | oracle |
|---|---|
| ``L`` | ``len(covered)`` |
| line ``p`` crossed | ``p−1 ∈ covered and p ∈ covered`` — `ACCUMULATOR_DESIGN.md` §2's definition, verbatim |
| node spanned | ``cuts[i]−1 … cuts[i+1]`` all covered |
| node contained | no junction used, and ``min``/``max`` of ``covered`` fall in one node |

⚠ §3.1 requires that *"whatever counts toward ``L`` must also count as coverage for crossing"*. **Nothing
tested the two against each other before.** They now come from one set, so they cannot disagree.

### Coverage

| | |
|---|---|
| **exhaustive**, 0 introns | 78 configurations |
| **exhaustive**, 1 intron | 7,098 |
| **exhaustive**, 2 introns | **326,508** |
| randomised, ≤4 introns, realistic coordinates (9 kb ref, 3 kb spans, 900 bp introns), fixed seed | 4,000 |

**333,684 exhaustive configurations** — the complete space over a 12 bp reference whose 4 nodes include a
**1 bp** node — plus a randomised sweep at a scale exhaustive enumeration cannot reach (a 300 bp molecule
spanning a long intron, the case §3.2's length limit turns on). **All pass.**

### ⭐ The proof was proven to FIRE — 7 perturbations

| | perturbation | |
|---|---|---|
| **L1'** | ``L`` from ``span − Σ(RAW intron lengths)`` — the formula §3.3 says goes **negative** on a wide overlap | ✅ caught |
| **L2'** | introns not clipped to the fragment | ✅ caught |
| **L4** | fragment not clipped to the reference | ✅ caught |
| **L5** | crossing boundary ``searchsorted`` side flipped right→left | ✅ caught |
| **L6** | spanning loop credits one node too many | ✅ caught |
| **L7** | containment keyed on the fragment EXTENT, not its first/last COVERED base | ✅ caught |
| **L3** | ABUTTING introns no longer merge | ⚠ **not caught** |

⚠ **L3 is a correct non-failure.** Merging abutting introns changes only the ``introns_absorbed`` QC
counter — set subtraction removes ``[10,20) ∪ [20,30)`` and ``[10,30)`` identically, so ``L`` and all four
populations are untouched. The counter is pinned separately by ``test_accumulator_spec``.

⭐ **And a first attempt at L1 was ALSO a no-op**, which is worth recording because it is the design's own
claim demonstrated: replacing ``Σ segments`` with ``span − Σ introns`` **after** ``_normalise_introns``
changes nothing, because normalisation is precisely what makes the two formulas agree. The bug §3.3 warns
about only exists on the RAW list — and that version is caught.

⚠ **An error the oracle caught in ME first**: my first oracle reported crossed lines by *cut* index where
the deposit writes ``edge_base + line − 1``. Writing the oracle independently forced the offset to be
re-derived rather than copied, which is the entire point of trap 1.

### ⛔ Scope, stated rather than implied

The fixture carries **no annotated junctions**, so every fragment is unspliced and the **spliced routing**
(``edge_spliced`` vs ``edge_unspliced``, junction credit, the containment block) is not exercised here;
``test_accumulator_spec`` covers it. That is the right scope for C1: the unconditional histogram bins by
``L``, and ``L`` does not depend on which population the fragment is routed to.

⚠ Proving the Python reference proves the C++ too — `test_accumulator_native_parity` gates the C++ on
byte-identity to it, so a defect here would have been reproduced faithfully in both.

### Verdict

✅ **`L` is fit to be the tool's one definition of fragment length. C1 may proceed.** Suite **1805 / 22**
(the 21 golden moves + the pre-existing `TODO.md` §7 failure; +5 new proof tests).

---

## C1 — ⭐ the accumulator gets an unconditional length histogram, and the frame mismatch is GONE (2026-07-31)

    Audit: `docs/FRAGMENT_LENGTH_AUDIT.md` §4 C1   ·   Gates: `tests/native/test_fragment_length_proof.py`,
    `tests/test_accumulator_payload.py`   ·   Precondition: C0 (`L` proven) — done first, by owner ruling

**What.** One `uint32[max_length + 1]` row, `deposited_lengths`, incremented on the same line as
`node_start_count` and the `DEPOSITED` counter. The five pure pools are untouched — `_pool()` still
returns `None` for a mixture, because an impure pool is worse than a missing one (`ACCUMULATOR_DESIGN.md`
§8). ⭐ **This is a separate tally, not a sixth pool**, and that distinction is the whole point.

**Why.** The accumulator binned only *conditioned* pools, so it had no unconditional histogram; the
empirical-Bayes shrinkage needs one; so the anchor was taken from the **scanner**, which measures length
by two other rules over another population. C1 removes the reason that borrowing existed.

⭐ **It is "unconditional GIVEN DEPOSIT", and the name says so.** It excludes what the accumulator
rejects — too long, ambiguous path, strand-undefined, empty — each of which is counted in `qc`. That is
**exactly** the population the pools are drawn from, which is what makes it the right anchor rather than
merely a convenient one; an anchor over a *wider* population would re-create the frame mismatch in a new
place.

### ⭐ THE MEASUREMENT: C1 separated the two mechanisms the audit could only assert were separable

On the zero-gDNA conditions every fragment is RNA, so the anchor and the RNA pool describe **one
population** and any gap between them is bias:

| | mean | sd |
|---|---|---|
| **the old scanner anchor** | 210.1 | 85.5 |
| ⭐ **C1, the accumulator's own** | **218.0** | **111.2** |
| the RNA_SPLICED pool | 234.7 | 146.7 |

| gap to the RNA pool | before C1 (scanner anchor) | **after C1** | closed |
|---|---|---|---|
| mean | **+11.6 %** | **+7.7 %** | **34 %** |
| sd | **+71.1 %** | **+32.0 %** | **55 %** |
| support ceiling | 713 vs 1000 | **1000 vs 1000** | ⭐ **entirely** |

⭐ **The frame component is gone and the opportunity component is now isolated and measurable.** §2 of the
audit named two mechanisms and could only argue they were separable; C1 separated them. The residual
**+7.7 % / +32 %** is the junction-opportunity tilt, and it is C3's (§8.1(b)'s) target — now against a
same-frame reference instead of a confounded one.

⚠ The old anchor's sd was too *small* (85.5) precisely because definition **B** is transcript-space and
unanimity-gated and its support stopped at 713 bp. An anchor narrower than the thing it anchors is the
worst direction for an EB prior: it shrinks the pools toward a distribution that cannot produce them.

> ⛔ **CORRECTION, added 2026-08-01 by C2.6 — THIS ENTRY READ THE SUPPORT CEILING BACKWARDS.**
> The "713 vs 1000 → 1000 vs 1000, ⭐ entirely" row above is **wrong, and in the reassuring direction**.
> **713 bp is the library's TRUE maximum**, read from `truth_fragment_lengths.tsv`; **1000 is
> `max_frag_length`, the clamp**. Definition **B** took its length from the *transcript*, so it could not
> produce an uncut intron — on the ceiling it was right and `L` was wrong, by an intron that was never
> sequenced. ⚠ The paragraph immediately above ("the old anchor's sd was too small... its support
> stopped at 713 bp") is the same error stated a second way. ⭐ The residual "+7.7 % / +32 %" was
> therefore **not one mechanism**: the mean is the junction tilt (C3's) and the sd was the uncut intron,
> now fixed — **the anchor's sd against truth is +1.98 %, from +26.97 %.** See the C2.6 entry.
> ⚠ None of this undoes C1 or C2, whose other reasons all still hold.

### Gates

| | |
|---|---|
| **G1 byte-identity** | the parity gate reads its field list off `dataclasses.fields(Tally)`, so `deposited_lengths` joined it **automatically**, dtype and shape included. ⭐ Proven to fire: three reference perturbations (stop binning · bin one short · **bin at the SPAN instead of `L`**) all caught |
| **G2 the invariant** | `Σ deposited_lengths == Σ node_start_count == qc.deposited`, incremented on one line so the three cannot drift by construction. Same externally-checkable form as §10.2's start count and a **different statement**: that one says every fragment was located in space, this one that every fragment was binned by length |
| **at the payload door** | `from_scan_result` refuses a histogram that does not sum to `qc.deposited`, or is the wrong length. ⚠ It must live at the boundary and not only in the accumulator's tests, because the payload is what a **cached** scan is rebuilt from |
| **superset** | every pooled fragment is also binned here at the same `L` — the property that makes it a usable anchor. Strict superset in general: an exonic contained fragment and a multi-line crossing enter no pool at all |
| **rejected ≠ binned** | too-long, ambiguous-path and empty fragments bin nothing |

⭐ **6 perturbations, and the 6th found a hole in the tests rather than the code**: deleting the payload's
door-check left everything green, because nothing asserted the check *existed*.
`test_a_deposited_lengths_HISTOGRAM_THAT_DOES_NOT_BIN_EVERY_FRAGMENT_IS_REJECTED` closes it and now fails
without the check.

### Blast radius, all as designed

* `payload_schema_digest` **66b41ea0b645209d → b7d29676c58b2c65**, so every existing scan cache is refused
  at load. ⭐ That is exactly what that key exists for (`scan_cache` docstring: S5.a's `length_sum` is the
  precedent). The 8 pilot caches were rebuilt — **56.0 s**, and all 8 satisfy G2.
* Cache serialisation is driven by `dataclasses.fields(AccumulatorPayload)`, so the new row persists with
  **no change to `scan_cache`**.
* ⚠ **A pre-existing defect fixed in passing**: `Tally` declared `edge_spliced_length_sum` and
  `sj_length_sum` **twice each**. Harmless (the dataclass de-duplicates to 18 fields) but it is a wart in
  the file that IS the executable specification, and it was committed.

**Suite 1809 / 22** — the 21 golden moves + the pre-existing `TODO.md` §7 failure; +4 C1 tests. Native
96/96 including parity and worker determinism.

⛔ **`build_fl_models` still reads the scanner anchor.** C1 makes the correct anchor *exist*; **C2** is
what switches the consumer over and deletes the scanner histogram. Nothing downstream has moved yet, by
design — this step is additive and every number in the tool is unchanged.

---

## C2 — ⭐ ONE definition of fragment length, and the scanner's histogram is DELETED (2026-08-01)

    Audit: `docs/FRAGMENT_LENGTH_AUDIT.md` §4 (C2.0–C2.5)
    Gates: `tests/test_one_fragment_length_definition.py` (the standing grep), `tests/test_d7_transcript_eff_lengths.py`,
           `TestSpliceCensus` + `TestFragmentLengthAnchor` in `tests/test_scanner_accumulator_integration.py`
    Measurement: `scripts/design/fl_anchor_gap.py` (new)   ·   Report A/B: this tree vs a worktree at `d045d820`

⭐ **The tool now has ONE definition of fragment length — the accumulator's `L` — measured in ONE place,
over ONE stated population.** D1, D2, D3 and D5 are closed; D6 is closed by deletion (the population that
was silently dropped no longer exists as a concept); D7 is **verified rather than assumed**, which the
audit demanded in those words.

### ⛔ THE OWNER'S RULING REJECTED THE PREMISE OF BOTH OPTIONS, AND IT WAS RIGHT

C2 was blocked on one decision: `rigel report`'s five per-fragment splice-type counts had no accumulator
equivalent. Option **(a)** grew five counters *in the accumulator*; option **(b)** reshaped the report and
lost two of five categories. The recommendation was (a). The owner ruled **neither**:

> *"There are some QC counts that the scanner must be responsible for. They truly live in the scanner …
> If the QC counts are generated in one place and have no algorithmic use, there's no need to pass them.
> … I don't see a reason to propagate artifacts into the accumulator, it's the scanner's job to identify
> and filter these out."*

⭐ **The principle the audit was missing.** §4's "ONE definition, ONE place, ONE population" governs
**model inputs**. The five splice counts are not a model input — they are a **census of the scanner's own
classification decisions**, with no consumer but the report. They were entangled with the histogram only
because they were *derived from it* (`category_models[stype].n_observations`), an implementation accident.
So C2.0 is not "move them to the accumulator"; it is **sever them from the histogram and leave them where
they are generated**.

⛔ **And option (a) would have been wrong on its own terms.** `SPLICE_ARTIFACT` fragments never reach the
accumulator — the deposit adapter returns early, because a blacklisted junction's span derives from an
alignment the scanner has already refused to believe — so (a)'s stated gate, *"the five counts sum to
`qc.deposited`"*, **cannot hold**. That was found by pricing the option, not by building it.

⭐ **What the ruling bought, measured**: `payload_schema_digest` stays **`b7d29676c58b2c65`**. No accumulator
change, no reopened S3 byte-identity gate, and **all 8 pilot scan caches remained valid** — against option
(a), which would have moved the digest and rebuilt every one.

### C2.0 — the splice census, and an identity that closes the books across two subsystems

One `std::array<int64_t, NUM_SPLICE_TYPES>` in `BamScanStats`, incremented at **one** site (the top of the
deposit adapter, before the artifact hold-out), plus **one** counter `n_deposit_not_offered`.

⚠ **An array, not five named fields**, so `census[st]++` is one statement, the merge one loop, the export
one loop over the *existing* `splice_type_label` table — and a category added to `SpliceType` cannot slip
past any of the three. ⭐ **No name table anywhere**: `splice_type_label`'s strings are exactly the
`SpliceType` member names lower-cased, so `rigel.splice.census_field` derives the same key on the Python
side. The correspondence is **tested**, not asserted — which matters because `_apply_scan_stats` copies with
`dict.get(key, 0)`, so drift reads as *zero*, not as an error.

⭐ **The gate as written in §4 was unsatisfiable; the replacement is stronger.** In C1's G2 form:

    Σ census − census[SPLICE_ARTIFACT] == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered

⭐ **And the artifact count is derived a SECOND way, by a subsystem that has never heard of an artifact.**
Scan one BAM against the same index with and without a splice blacklist: blacklisting relabels fragments
`SPLICED_ANNOT → SPLICE_ARTIFACT`, so `qc.deposited` — which knows nothing of splice types or blacklists —
must fall by **exactly** the artifact census. Measured 200 → 40 against a census of 160. `CARRY_FORWARD.md`
§3 trap 1: a check that re-derives the number by the same route checks nothing.

**6 perturbations, all caught — but TWO needed fixtures built for them, and that is the finding:**

| | perturbation | |
|---|---|---|
| **PC1** | one category never censused | ✅ caught by 3 tests |
| **PC2** | census moved AFTER the artifact hold-out | ✅ caught — ⚠ **only by the blacklist fixture** |
| **PC3** | artifacts counted as `not_offered` as well | ✅ caught |
| **PC4** | one C++ export key renamed (the `.get(key, 0)` silent-zero) | ✅ caught by 3 tests |
| **PC5** | the multi-reference hold-out made silent again | ⛔ **NOT caught at first** |
| **PC6** | census narrowed to only what deposits | ✅ caught by 3 tests |

⛔ **PC5 is the C1-perturbation-6 pattern again**: `n_deposit_not_offered` was zero on both sides of the
identity in every fixture, so deleting it left the suite green. Unlike C0's L3 that is **not** a proven
no-op — it is a hole. Closed by `multi_reference_bam`, a hand-written two-contig BAM whose mates land on
different references and which reaches the adapter *by the intergenic path* — precisely the route the
shipped defect took. PC5 now fires.

### C2.1 — the anchor moved, and the fix is STRUCTURAL

`build_fl_models(payload)` — **the payload is the only argument**. All three histograms are read off one
object in one frame, so the mixed-frame call is **unrepresentable**, not merely discouraged. The EB kernel
over three free histograms survives as `_fl_models_from_histograms`, which production never calls.

⚠ **A value gate alone would have been vacuous.** On the plain oracle fixture the scanner's histogram and
`deposited_lengths` are **byte-identical** — a perfect BAM with no ambiguity makes definitions A/B and C
agree — so the test would have passed whichever anchor was wired in. A byte-identical result is no
evidence (`benchmark_ab_methodology_cautions`). The blacklist fixture separates them (200 vs 40).

⭐ **§2's table re-measured on all 8 pilot conditions — the gate's stated numbers, to the decimal:**

| `gdna_none ss0.50 capture_off` | anchor | RNA pool | gap |
|---|---|---|---|
| mean | **218.0** | 234.7 | **+7.7 %** |
| sd | **111.2** | 146.7 | **+32.0 %** |
| support ceiling | 1000 | 1000 | ⭐ matched |

C1's numbers exactly, against the pre-C1 **+11.6 % / +71.1 %**. The frame component is gone from the
*shipped* code, not merely available to it.

⚠ **NEW, and C3 should expect it: the tilt is capture-dependent.** `capture_on` reads **+11.0 % / +43.4 %**
against `capture_off`'s +7.7 % / +32.0 %. §8.1(b) is aimed at the junction-opportunity tilt; there is a
second effect riding on it.

⚠ **Three test harnesses had already drifted from production** — `_oracle.py`, `test_ambig_scenario`,
`test_accumulator_span_unbiased` still built FL models from the scanner's `category_models`, which
production stopped using at S5.d. They were calibrating against a model the tool does not ship. Converged.

### C2.2 / C2.3 — the deletion, and the report

Deleted: sites **A** and **B**, `FragLenObservations`, `frag_length_observations`,
`_replay_fraglen_observations`, `FragmentLengthModels` (**plural**), `n_frag_length_{un,}ambiguous`, and
`get_unique_frag_length_mrna` — **definition B itself**. `scan_and_buffer` returns a 4-tuple; 14 files
followed.

⚠ **`FragmentLengthModel` (SINGULAR) stays** and has its own survivor test. The two names differ by one
character and the singular is the scorer; a search-and-delete that caught it would have removed the
fragment-length term from scoring and reported it as numerical drift.

⭐ **The gate is a standing search over `src/`, not a one-off command** — `tests/test_one_fragment_length_definition.py`.
A partial delete that still compiles is the failure mode. ⚠ It found **three genuine stale references**
that a build could not: a nanobind docstring still advertising `'frag_length_observations'` as a returned
key, a comment in `resolve_context.h`, and a `config.py` docstring. It exempts an explicit
`DELETED by C2` tombstone and nothing else — a comment that merely *mentions* a deleted symbol is how the
next reader concludes it still exists.

**The report A/B against a worktree at `d045d820`, same scenario, same seed, `assignment_mode="map"`:**

| splice key | `d045d820` | **C2** |
|---|---|---|
| `unspliced` | 1482 | **1547** |
| `spliced_annotated` | 412 | 412 |
| `spliced_unannotated` | 0 | 0 |
| `spliced_implicit` | **8** | **41** |
| `splice_artifact` | 0 | 0 |
| **total** | **1902** | ⭐ **2000** |

⭐ **The key set is IDENTICAL — nothing was lost — and the old numbers were missing 4.9 % of the library.**
2000 is exactly the simulated fragment count. That reproduces G6's 4.6 % scanner-drop rate on a completely
different fixture, and it is the first time the report's splice section has summed to the library.
⚠ **`spliced_implicit` was under-counted 5×**: an implicit splice is exactly the case whose transcript-space
length most often fails the unanimity gate, so it was the worst-hit class.

`fragment_lengths.feather`'s categories change by design: the per-`SpliceType` histograms give way to
`global` + `gdna`/`rna` + the **five pure pools**. ⭐ That surfaces `splash_fl_mass` — the two crossing
pools, the only **on-target** gDNA population — for the first time; its own docstring asked to be reported
("makes that comparison an output instead of an assumption") and nothing was doing it.

### C2.4 / C2.5 — the dead field, and D7 CHECKED

`ScanCache.fl_global_counts`, `fl_rna_counts` (**D5** — written, read back, consumed by nothing) and
`fl_max_size` (a duplicate of `payload.max_length`) are gone, with `fl.npz` itself. ⚠ **Caches written
before C2 still load** — an extra file on disk is not a key — and all 8 pilot caches were verified loading
after the change.

⭐ **D7 is verified, not assumed** — the audit's own instruction. The shipped per-transcript effective
lengths are asserted **exactly equal** to those rebuilt from the payload alone, end to end through
`run_pipeline`. Two supporting gates: the quantity is *sensitive* to the pmf (shifting it 50 bp shortens
every effective length — an equality against a constant array would prove nothing), and perturbation
**PD1** (feeding `global_pmf` instead of `rna_pmf`) is caught.

⚠ **A finding, pinned rather than fixed**: `pipeline`'s guard `if rna_fl.n_observations > 0` **cannot
fire**. `from_pmf` sets `_total_weight = 1.0`, so `n_observations` is exactly 1 whatever went in. The
reachable behaviour is correct — an empty RNA pool EB-shrinks to the unconditional anchor, which beats a
hard-coded 200 bp mean — but the guard *looks* like a data-presence fallback, and that appearance is what
would let a future empty-pool bug hide behind it.

### ⚠ Unrelated: `scripts/sim/fl_estimation_stress.py` was ALREADY DEAD and is deleted

It called `calibrate(index=…, scan_trained=…, fl_prior_ess=…, pool_quality_good=…)` — **not one of those
is a parameter of the current `calibrate`**. Dead since well before this branch; the 2026-07-30 sweep
missed it. Recoverable: `git checkout d045d820 -- scripts/sim/fl_estimation_stress.py`.

### Suite

**1824 passed / 21 failed**, from a re-measured baseline of **1809 / 22**. Accounting is exact:
`1809 + 7 (C2.0) + 2 (C2.1 net) + 10 (C2.2 gate) + 2 (C2.3) + 3 (C2.5) − 9 (deleted container tests) = 1824`.

⭐ **The failure count went DOWN by one: `TODO.md` §7's `test_nrna_double_counting[g20_n0_s100]` now
PASSES.** The silent negative control reads **0 counts against a limit of 25**, where it leaked ~30 before
— not a marginal pass, and stable across three re-runs. ⚠ **Do not close §7 on this.** A negative control
is one-sided (trap 19) and this was not the change's target; the honest statement is that correcting the
fragment-length models removed the leak on this condition, and §7 should be re-read against the other
modes before it is retired.

⚠ **The 21 remaining are `test_golden_output` and they have now moved TWICE** — P1's units fix and C2's FL
models. ⛔ Still **do not regenerate**: they move again at C3. Regenerate **once**, after C3, **twice, and
diff**.

### ⚠ The standing check's INVOCATION matters, and it cost a false alarm

⛔ Use **`python -m pytest`**, not bare `pytest`. Bare `pytest` does not put the repo root on `sys.path`,
so `tests/calibration/test_fl.py`'s `import tests.native._accumulator_reference` raises and the suite reads
**1808 / 23**. That 23rd failure is an artefact of the invocation, not a regression — and under the rule
"a 23rd failure is a regression" it reads as one. `CLAUDE.md`'s command should say `python -m pytest`.

---

## C2.6 — an intron in the mate gap is cut on EVERY fragment, not only unspliced ones (2026-08-01)

    Spec: `docs/SPEC_GAP_INTRONS.md`  ·  Cause and evidence: `docs/JUNCTION_OPPORTUNITY.md` §4
    Owner ruling, 2026-08-01: *"we should be searching for gap introns within every fragment."*
    Baseline: e2f41cd0, re-measured in this session — **1824 passed / 21 failed**, all 21 goldens.

C0 proved the accumulator's `L` correct **given its inputs**. This entry is about an input that was
incomplete. `resolve_context.h` ran implicit-splice detection only on fragments the resolver had already
classified `SPLICE_UNSPLICED`, so a fragment carrying an observed CIGAR-N splice never had its
**unsequenced mate gap** examined, and an annotated intron sitting in that gap stayed inside `L`.

⚠ `UNSPLICED` never meant "one aligned block" — an unspliced paired-end fragment already has two blocks
and a mate gap, and that case always worked. The missed population is **spliced fragments that also have
a gap intron**, necessarily long because they span two or more introns. That is exactly the tail.

### The change, in three places

| | | |
|---|---|---|
| **§2.1** | `collect_implicit_splice_introns` takes the fragment's **observed** introns and drops every gap that matches one by **exact `(start, end)` equality** | ⛔ the trap: the gap finder walks consecutive aligned blocks, so a CIGAR-N intron is a "hole" too |
| **§2.2** | detection is **unconditional**; the `SPLICE_IMPLICIT` **promotion** stays unspliced-only | `splice_type` feeds scoring, the buffer, strand training and the report's census. This work is about `L`, not classification |
| **§2.3** | the deposit adapter **unions** the two intron lists (sorted, then de-duplicated), `sj_implicit = !implicit_introns.empty()`, `path_ambiguous = implicit_ambiguous` | the two lists stopped being mutually exclusive |

⭐ **The exact-equality rule is the whole of §1.** `transcript_has_implicit_intron_in_gap` accepts any
annotated intron inside the gap **± K** (K = 3 by default), so an overlap-based filter — or none — lets a
**different** nearby intron answer for one the CIGAR already stated. The two then normalise into one
wider interval and `L` comes out too **SHORT**. Measured under X2 on the hand-written fixture: the
near-match fragment reads **500 bp where the molecule is 502**.

### ⭐ THE HEADLINE — `gdna_none ss0.99 capture_off`, scored on `truth_fragment_lengths.tsv`

Truth: mean **217.13**, sd **87.41**, ceiling **713 bp**. Every target is read from that file or is a
control; **nothing here is tuned**.

| the anchor (`deposited_lengths`) | before | after | |
|---|---|---|---|
| **sd vs truth** | **+26.97 %** | ⭐ **+1.98 %** | **G-sd — 92.7 % of the excess removed** |
| mean vs truth | +0.38 % | −1.88 % | ⚠ see the selection effect below |
| mass **> 713 bp** (truth: **0**) | 0.00909 | ⭐ **0.00137** | **G-tail — 85 % removed, ⛔ not 0** |
| mass ≥ 700 bp | 0.00965 | 0.00143 | |
| `qc.dropped_too_long` | **280,558** (5.71 %) | ⭐ **38,309** (0.80 %) | 86 % collapse |
| `qc.deposited` | 4,636,640 | 4,770,233 | |
| `qc.dropped_ambiguous_path` | 82,802 | 191,458 | +108,656 — the new deferrals |
| `qc.sj_implicit_fragments` | 46,700 | 323,654 | 6.8× — the mixed population, seen for the first time |
| `qc.introns_absorbed` | 0 | **0** | the two lists really are disjoint |

### ⛔ G-gdna — THE CONTROL, AND IT DID NOT MOVE ONE DIGIT

`DNA_INTERGENIC` is pure gDNA and gDNA **has no introns to miss**, so a change there would mean the fix
reached fragments with no introns — impossible, therefore a bug. On `gdna100 ss0.99 capture_off`:

| | before | after |
|---|---|---|
| mean vs truth gdna | −0.5976 % | **−0.5976 %** |
| sd vs truth gdna | −0.3424 % | **−0.3424 %** |
| ≥600 / ≥700 bp | 0.00023 / 0.00001 | **0.00023 / 0.00001** |
| ceiling · n | 785 · 4,527,474 | **785 · 4,527,474** |

⭐ Bit-identical on every statistic **including the fragment count**. The control is clean.

### ⛔ G-tail DID NOT REACH 0 — and the residual is MEASURED, not excused

The gate says the mass above the true ceiling must be 0. It is **0.00137**. Per the spec's §6 the
remaining mass was measured rather than closed with a constant, by an experiment (**M3**) that emits
**every** annotated intron inside a gap instead of the first:

| | shipped | M3 (all introns in the gap) |
|---|---|---|
| mass > 713 bp | 0.00137 | ⭐ **0.00002** |
| `dropped_too_long` | 38,309 | ⭐ **389** |
| anchor sd vs truth | +1.98 % | −1.28 % |

⭐ **D3 is 98.5 % of the residual, and it is now a measurement rather than a hypothesis.**
`transcript_has_implicit_intron_in_gap` returns the **first** matching intron and stops, so a mate gap
spanning two annotated introns keeps only one cut. ⚠ M3 is a MEASUREMENT ONLY — it emits multiple
introns without extending the per-gap unanimity test to compare intron *sets*, which is the real work.
`TODO.md`-grade, and it is now the only known mechanism left in the tail.

### ⛔ D1 WAS MEASURED AND IT COSTS MORE THAN IT BUYS — the owner's call to revisit

D1 (recommended in the spec, implemented as ruled) removes a mixed fragment from `RNA_SPLICED`, because a
length partly inferred from the annotation is a product of the model that pool is used to fit. The spec
required the two effects to be reported **separately**; they were, by an experiment (**M1**) that cuts the
intron but leaves the fragment in the pool.

| `RNA_SPLICED` vs truth | before | ⭐ M1: intron cut, fragment KEPT | shipped: intron cut, fragment REMOVED |
|---|---|---|---|
| **mean** | +8.00 % | ⭐ **+0.67 %** | ⛔ **−9.58 %** |
| **sd** | +67.35 % | ⭐ **+2.40 %** | ⛔ **−22.46 %** |
| n | 1,475,626 | 1,609,219 | 1,332,265 |
| mass > 713 bp | 0.02812 | 0.00364 | 0.00000 |

⭐ **Cutting the intron very nearly fixes the pool (+8.00 → +0.67 mean, +67.35 → +2.40 sd). Removing the
mixed fragments then breaks it again in the opposite direction**, because the fragments removed are
exactly the ones whose mates sit far apart — a **length-selection bias**, and the pool is what the
fragment-length model is fitted from. ⚠ The purity argument for D1 is unchanged and still real; what is
new is its price, and it is 10 % of the pool's mean. Reversing it is one line
(`path.sj_implicit = implicit_only`), which is literally M1. **Owner decision.**

### ⚠ The anchor's mean went +0.38 % → −1.88 %, and about half of that is SELECTION

`dropped_ambiguous_path` rose by 108,656 (2.3 % of deposits) because a gap whose candidates disagree now
defers the fragment even when its *other* splice was observed. Those fragments are long by construction,
so the surviving anchor is biased short. Experiment **M2** lets them deposit again:

| anchor | shipped | M2 (ambiguous still deposits) |
|---|---|---|
| mean vs truth | −1.88 % | **−0.83 %** |
| sd vs truth | +1.98 % | +4.44 % |
| mass > 713 bp | 0.00137 | 0.00183 |

⭐ So ~1.05 pp of the −1.88 % is the deferral, not a measurement error — and deferring is the right call
(`L` genuinely is undetermined), so this is a cost recorded, not a defect.

### Gates

| | | |
|---|---|---|
| **U1** | a spliced fragment's mate-gap intron is found, and `L` excludes **both** cuts | ✅ resolver-level and end-to-end (`deposited_lengths` reads exactly `{300, 500, 502}`) |
| **U2** | the observed CIGAR-N intron is **not** re-derived | ✅ |
| **U3** | a near match within ±K is **not** substituted, and `L` does not shrink | ✅ |
| **U4** | `splice_type` does not move; the census is unchanged | ✅ 3 `SPLICED_ANNOT`, 1 `SPLICED_UNANNOT`, **0** implicit |
| **U5** | a MIXED fragment whose candidates disagree is rejected **and counted** | ✅ |
| **D1** | only the fully-observed fragment stays in `RNA_SPLICED` | ✅ two-sided |
| **S3** | byte-identity to `_accumulator_reference.py`, bit-identical at 1/2/4/8 workers | ✅ automatic, unchanged |

⚠ **`_accumulator_reference.py` did NOT move, and that is correct.** The specification describes the
**accumulator**, whose contract is `FragmentPath` in and a tally out. Which introns a fragment *has* is
decided upstream, in the resolver and the deposit adapter. `sj_implicit`'s meaning — "this `L` depends on
an intron that was never sequenced" — is unchanged; what changed is which fragments satisfy it. S3 never
reopened, and `payload_schema_digest` never moved.

### The perturbations

| | perturbation | result |
|---|---|---|
| **X1** | restore the `splice_type == SPLICE_UNSPLICED` gate | ⭐ **6 gates fail** |
| **X2** | drop the §2.1 observed-gap filter | ⭐ **7 fail**, including U2 and U3. `near` reads **500 bp against 502** — the "L too SHORT" defect, exactly as §1 predicts. ⚠ Also caught by the worker-determinism module's own non-vacuity guard: every spliced fragment re-derives its own observed intron, so **`pool_lengths` sums to zero** |
| **X3** | `sj_implicit` back on the splice class | ⭐ **D1 alone fails** — precisely the predicted gate, nothing else |
| **X4** | `path_ambiguous` gated on implicit-only | ⭐ **U5 fails** (+2 collateral) |
| **X5** | leave the union **unsorted** | ⛔ **NOTHING FAILS — the spec's prediction was wrong** |

### ⛔ X5: the spec predicted S3 byte-identity would break. It does not, and here is why

`normalise_introns` **sorts the intron list itself**, in the C++ accumulator and in
`_accumulator_reference.py` alike. So an unsorted union cannot change `L`, the segments, the junction ids
or any tally — the adapter's sort protects exactly one quantity, `qc.introns_absorbed`, by keeping the
adjacent-pair de-duplication a de-duplication.

⚠ And the case where that matters is close to unreachable: an implied intron can only duplicate an
observed one under K-slack, and it must then be **non-adjacent** in the concatenation, which needs a
fragment with ≥2 observed introns *and* a coincident K-slack match. The sort is kept because it costs
nothing on lists of ≤4 and restores the invariant the de-duplication was written against — but it is
**defensive, not load-bearing**, and this entry records that rather than claiming a gate it does not have.
`qc.introns_absorbed == 0` is now asserted on the fixture, which pins §2.1's disjointness claim instead.

### ⚠ Instrumentation added, and why it is not test-only scaffolding

`ResolvedFragment` now exposes `implicit_introns` and `implicit_ambiguous`. Without it **U2 is
unobservable**: the adapter de-duplicates and the accumulator normalises, so an intron re-derived from an
observed splice is merged away and "never emitted" is indistinguishable from "emitted and absorbed". §1
is precisely about the case where they differ, so the emission itself has to be visible. ⛔ It is
**copied, not moved**, out of `RawResolveResult` — `bam_scanner.cpp` builds the `ResolvedFragment` before
calling `deposit_to_accumulator(frag, cr)`, so a move would hand the adapter an empty list and silently
stop cutting every gap intron.

### ⚠ Cost: +8.9 % scan time

`gdna100 ss0.99 capture_off`, 10,000,000 fragments, `OMP_NUM_THREADS=1`, best of three:
**845.4 ns/fragment against 775.6** with detection gated back to unspliced-only. Detection is now a
binary search per candidate per gap on every fragment rather than on the unspliced ones.

### Blast radius, as predicted

`payload_schema_digest` **did not move** (no new field), so the 8 pilot caches stayed loadable — ⛔ but
their contents were wrong and all 8 were rebuilt (45.6 s). The goldens move a **third** time. ⛔ Still
regenerate **once**, after C3, **twice, and diff**.

### Suite

**1835 passed / 21 failed**, from a re-measured baseline of **1824 / 21** in the same session.
Accounting is exact: `1824 + 4 (U1–U4, resolver level) + 7 (end-to-end module) = 1835`. All 21 failures
are `test_golden_output`, unchanged in identity from the C2 entry.

### Files touched

`src/rigel/native/resolve_context.h` · `src/rigel/native/bam_scanner.cpp` · `src/rigel/native/resolve.cpp`
· `scripts/design/fl_anchor_gap.py` (the truth panel: G-tail, G-sd, G-gdna, D1, D3) ·
`tests/test_implicit_splice.py` · `tests/native/test_gap_introns_are_cut.py` (new)

### ⚠ A correction to the record, carried from `JUNCTION_OPPORTUNITY.md` §4.4

**C1 read the support ceiling backwards.** It recorded the anchor's ceiling moving from the scanner's
**713** to the accumulator's **1000** as "the support ceiling mismatch closed entirely". ⛔ **713 was the
library's true maximum**, to the base pair; 1000 is `max_frag_length`, the clamp. The deleted definition
**B** took its length from the *transcript*, so it could not produce an uncut intron — on this one
question it was right and `L` was wrong. ⚠ That does not undo C2, whose reasons all still hold.

---

## S1 + S2 — the accumulator ARBITRATES, and the side buffer exists (2026-08-01)

    Plan: `docs/PLAN_TWO_PASS.md` §4 (S1, S2)   ·   Spec: `docs/SPEC_GAP_PATHS.md` §0–§4
    Baseline: 750cc8ee — **44 failed / 1824 passed**, of which 21 are the stale goldens and
    **23 were a half-finished interface change**. After: **21 failed / 1863 passed**, all 21 goldens.
    ⛔ NOT judged by a calibration A/B, by design — see "what this deliberately does not measure".

`750cc8ee` moved the SPECIFICATION (`tests/native/_accumulator_reference.py`) to the hypothesis
interface and left the C++ behind it. This entry finishes that migration and lands the side buffer.

> A fragment arrives at the accumulator with its **hypothesis set**, not with one path and a flag. Exactly
> one survivor **deposits**; two or more and `L` is not one number, so the fragment is **held WHOLE** for
> the second pass. The empty hypothesis is the genomic one and needs no flag.

### What landed

| | | |
|---|---|---|
| **the C++ caught up** | `DeferredFragments` gains `ref`, is int64 throughout, and gains a canonical `sort` | the bank is the ONE tally whose order is observable; everything else is a sum of integers |
| **`Accumulator` knows which reference it IS** | `ref_id` is a required ctor argument, stamped into every held record | the drain replays a record through `deposit` onto **that** reference's cut axis |
| **the binding exposes both new banks** | `Accumulator.deferred` and `.gap_resolution`, through the same two exporters `build_result` uses | one set of key strings on this side of the ABI, so the parity surface and the payload cannot disagree |
| **the payload carries them** | `AccumulatorPayload.deferred` (a validated CSR) and `.gap_resolution` (a typed census) | |
| **S2: the cache round-trips them** | nested banks split by sub-field type: arrays to the `.npz`, counters to the manifest | ⛔ `dataclasses.asdict` on the bank yields ndarrays, and the manifest is written with `json.dumps(..., default=str)` — which would stringify each array to a **truncated repr**, silently |

### ⭐ THE OWNER'S KNOWN DEFECT IS FIXED, AND IT IS NOW GATED

`payload_schema_digest()` hashed `AccumulatorPayload`'s top-level field names **alone**, so a change
inside a nested bank was invisible to it — a renamed `ScanQC` field let a stale cache be accepted by the
key and then fail deep in the loader with a bare `TypeError`, exactly the failure the digest exists to
prevent. It recurses one level now. ⛔ S1 raised the stakes rather than lowering them: `DeferredFragments`
puts **thirteen** array names inside one payload field and every one is an `.npz` key.

⚠ **The fix had no gate until a perturbation said so** — see X8 below.

### ⭐ FINDING 1 — `GapResolution.RESOLVED_UNSPLICED` COULD NOT BE ENTERED, and it is deleted

The census had a class documented as *"the genomic hypothesis survived alone because every spliced path
was longer than `max_fragment_length`"*. **That is impossible.** A spliced hypothesis CUTS bases the
genomic one keeps, so `L_spliced <= L_genomic` always; the one arbitration filter is
`L <= max_fragment_length`; therefore if the genomic path survives the filter **every** spliced path does
too, and the survivor set can never be `{genomic}` while a spliced path was offered — which is the
condition for being in the census at all.

⭐ **Found by writing the non-vacuity gate**, not by review: the parity battery asserts every census class
is reached, and this one could not be. Checked before deleting — 200,000 random hypothesis sets, the
genomic path was the longest in **every one**.

⭐ **The ordering is now pinned as the REASON rather than the class as a consequence**
(`test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST`, 20,000 randomised sets). A future filter that broke
the ordering fails there instead of quietly needing a deleted class back.

### ⛔ FINDING 2 — THE REFERENCE STAMP WAS UNTESTABLE, AND TWO FIXTURES HAD TO MOVE

Two perturbations **passed the entire suite of 1860 tests**:

| | perturbation | why nothing caught it |
|---|---|---|
| **X2** | `deferred_.append(fragment, 0, …)` — a constant instead of the accumulator's own id | every accumulator-level fixture was reference 0 |
| **X3** | `AccumulatorSet` gives every reference the id `0` | every scan fixture was single-contig, or deferred only on reference 0 |

⛔ A constant is indistinguishable from a correct value until two references both hold something, and the
consequence is not cosmetic: the second pass replays each held record onto `ref`'s cut axis, so a constant
stamp drains one chromosome's coordinates onto another's partition — **the defect the predecessor deposit
adapter actually had**. Two fixtures changed rather than one assertion added:

* `test_accumulator_native_parity.py` now sits on **reference 3**, with three leading references that
  contribute no cuts (legal, and it exercises the per-reference offset arithmetic too);
* `test_gap_introns_are_cut.py` is now **two contigs**, with the disagreeing-isoform gene on chr2.

Both perturbations fire after the change. ⚠ The extent check is the half that matters most: `[66100,67000)`
is a plausible interval on chr1 as well, so "stamped 1" and "inside reference 1's cut range" are different
statements and only the second fails against the real partition.

### Behaviour that MOVED, deliberately, and where it is recorded

| | change | why |
|---|---|---|
| **D1 is deleted** | the RNA pool is keyed on **determinacy**, not provenance | a fragment enters a pool only when exactly ONE hypothesis survived, so its `L` is not in doubt however it was arrived at. C2.6 measured the alternative: **+0.67 % / +2.40 %** against truth on determinacy, **−9.58 % / −22.46 %** on provenance. ⭐ A purity filter on a length pool is a length filter |
| **D3 / C2.7 is solved** | every annotated intron in a gap is cut, and candidates group by their **whole-fragment path** | two transcripts differing only in their SECOND intron are now two hypotheses. ⛔ Its measurement cannot be re-read until the drain — see below |
| ⚠ **a SPLICED_UNANNOT fragment with a gap intron now DEFERS** | `SPEC_GAP_PATHS.md` §2: ∅ is available whenever no **annotated** junction was sequenced | the `near` fixture's junction is 2 bp off the annotation, so the molecule is not certified RNA and the genomic explanation survives. Its `L` gate moved from "deposits at 502" to "the enumeration produced the intron set that yields 502 and never the near-miss that yields 500" — asserted on **coordinates**, which is stronger |
| ⚠ **a transcript the read CONTRADICTS is no longer a candidate** | `SPEC_GAP_PATHS.md` §3, concern C1 | `test_h_any_intron_satisfies`'s old geometry ran a block contiguously across the transcript's FIRST intron, so the read has bases that transcript splices out. Cutting its *other* intron on that transcript's authority is a path no molecule has. The case keeps its meaning on a geometry that does not contradict, and the contradiction is now its own gate (`test_h2_…`) |
| ⚠ **the `max_fragment_length` prefilter changes CLASSIFICATION** | `SPEC_GAP_PATHS.md` §8, concern C3 | the same annotation defers or resolves depending only on how far apart the exons sit — a 300 bp intron holds the fragment, an 1800 bp one rules the genomic path out by length and deposits it. ⭐ Now gated **two-sided on a real scan** rather than assumed |

### Gates

| | | |
|---|---|---|
| **byte-identity** | C++ ↔ the specification, over every `Tally` field including the deferred CSR and the census, read off `dataclasses.fields(Tally)` | ✅ 54 named cases + 10,000 random fragments with random hypothesis sets |
| **non-vacuity** | the battery reaches **every** census class, and defers > 100 of the 10,000 | ✅ — this is what found FINDING 1 |
| **conservation** | `deposited + deferred + dropped_* == offered`, and the bank holds the fragments the counter claims | ✅ at the reference level, at the payload door, and on a real scan against the scanner's **independent** splice census |
| **determinism** | bit-identical at 1/2/4/8 workers, deferred bank included | ✅ at the accumulator (2/4/8 shards) and end-to-end through the scan |
| **S2 round trip** | the side buffer survives the cache byte-identically and comes back **typed** | ✅, and the fixture is asserted non-empty first |

### The perturbations

| | perturbation | result |
|---|---|---|
| **X1** | the canonical sort does nothing | ⭐ **5 fail** — parity and worker determinism |
| **X2** | the `ref` stamp is a constant | ⛔ **passed 1860 tests**; after the fixture moved, **7 fail** |
| **X3** | `AccumulatorSet` stamps every reference 0 | ⛔ **passed 1860 tests**; after the fixture moved, **1 fails** — the one that names it |
| **X4** | the census swaps `rna_or_gdna` and `which_introns` | ⭐ **6 fail** |
| **X5** | the fragment is COUNTED as deferred but not HELD | ⭐ **12 fail** — including the payload's door check |
| **X6** | the length filter keeps only the first hypothesis when the filter empties the set | ⭐ **1 fails** — precisely the case written for it |
| **X7** | `lex_compare` treats a PREFIX as equal | ⭐ **3 fail** |
| **X8** | `payload_schema_digest` stops recursing | ⛔ **not caught** — the fix had no gate; one added, then **1 fails** |
| **X9** | the cache writes a nested bank's arrays into the manifest JSON | ⭐ **2 fail** |

### ⛔ WHAT THIS DELIBERATELY DOES NOT MEASURE

**No calibration A/B, and that is the plan's instruction** (`PLAN_TWO_PASS.md` §4). Between S1 and S3 the
tally is deliberately **thinner**: the ambiguous mass is held, not yet deposited. And the held population
is the **long** one — a longer gap admits more hypotheses — so the surviving anchor is biased short.

⛔ **Every junction-opportunity number in the docs predates the drain and must be re-measured after S3,
not before.** So must D3's residual, which is a statistic of the deposited set.

### Files touched

`src/rigel/native/calibration/accumulator.{h,cpp}` · `src/rigel/native/bam_scanner.cpp` ·
`src/rigel/scan_payload.py` · `src/rigel/scan_cache.py` · `tests/native/_accumulator_reference.py` ·
`tests/native/test_accumulator_native_parity.py` · `tests/native/test_gap_hypothesis_arbitration.py` ·
`tests/native/test_gap_introns_are_cut.py` · `tests/native/test_fragment_length_proof.py` ·
`tests/native/test_implicit_splice_deposit.py` · `tests/native/test_accumulator_worker_determinism.py` ·
`tests/test_accumulator_payload.py` · `tests/test_scan_cache.py` · `tests/test_implicit_splice.py` ·
`tests/calibration/_synthetic.py` · `tests/calibration/test_gdna_strand_integration.py` ·
`scripts/design/fl_anchor_gap.py` · `scripts/design/implicit_splice_census.py`

### Blast radius

⛔ **`payload_schema_digest` MOVED** — it recurses now, and the payload gained two fields. Every scan cache
is invalidated **by design**; the 8 pilot caches rebuild (~46 s). The goldens move a **fourth** time.
⛔ Still regenerate **once**, at the very end, **twice, and diff**.

---

## S2.1 — the pipeline is REPRODUCIBLE, and it was never the seed (2026-08-01)

    Gate: `tests/test_scan_order_independence.py`   ·   Retires: `TODO.md` rank 3
    Owner question: *"the RNG seed issue seems minor — is there an easy solution?"*
    Baseline before: 21 failed / 1863 passed. After: **21 failed / 1868 passed**, all 21 goldens.

`TODO.md` rank 3 recorded that `EMConfig.seed` *"does not make `assignment_mode="sample"` reproducible"*
and filed it under the seed. ⛔ **The seed was never the problem**, and the recorded diagnosis sent the
investigation to the wrong subsystem.

### What was measured, before anything was changed

| | |
|---|---|
| a seeded `sample` run, twice at the default thread count | ⛔ counts differ by up to **22** |
| the `rng_seed` handed to the C++ EM on those two runs | ⭐ **bit-identical** — `9008465316566228445` |
| the EM's own thread count (`EMConfig.n_threads` 1 / 2 / 4 / auto) | ⚠ **irrelevant** — fails at all of them |
| scan `total_threads` = 1, 2, 4 | ✅ reproducible |
| scan `total_threads` = 8, 16 (**the default** = `os.cpu_count()`) | ⛔ **not** reproducible |
| the fragment buffer, two 16-thread scans | ⛔ **different row order** (1 chunk → identical; 6 chunks → not) |
| the `frag_id → fragment` mapping, same two scans | ⭐ **byte-identical**, ids contiguous `0..n-1` in BAM order |

⭐ **So: the scanner streams finalized chunks in worker-COMPLETION order, the fragment buffer inherits
that, and the EM's sampled assignment consumes its per-locus RNG in buffer-row order.** Permuting the
units permutes which fragment gets which draw, and units in different equivalence classes have different
posteriors, so the assignment genuinely moves. ⚠ The 1- and 4-thread passes were an artifact of scale —
one chunk, nothing to reorder. On a real library, which emits many chunks at any thread count above one,
the shipped default was **not reproducible at all**.

### The fix — an ORDERING, not a seed, in one place

```python
unit_indices=_in_bam_order(comp_u_indices[lo:hi], em_data.frag_ids)
```

`build_multi_loci` orders each locus's units by `frag_id`. Every per-locus array is scattered by that
index list (`locus_partition.partition_and_free`), so **all of them inherit the canonical order and no
consumer has to know**. ⚠ `unit_indices` is used nowhere else except for its length, so nothing else moves.

### ⛔ WHY NOT SEED EACH DRAW FROM THE FRAGMENT — the plausible wrong fix, and it is caught

Hashing the fragment's own content would also be order-independent, and would be **wrong**: identical
fragments hash alike and so draw alike, collapsing a 60/40 posterior to **100/0** for every group of
duplicates. That is the owner's own D-D ruling (`PLAN_TWO_PASS.md` §5.3) arrived at from the other side.

⭐ **And it fails even on its own terms.** A content key **ties exactly on the duplicates**, so the tie
falls back to buffer order — measured (perturbation Y2): ordering by `gdna_log_liks` instead of `frag_id`
leaves **both** `sample` and `fractional` thread-count dependent. `frag_id` is an identity, so it never
ties; ordering by identity and keeping one stream per locus is what preserves the multinomial spread.

### ⭐ THE CODE GOT SMALLER — 46 lines deleted

`em_solver.cpp`'s `build_equiv_classes` sorted the ROWS within each equivalence class on a
log-likelihood fingerprint, ~45 lines, explicitly *"because multi-threaded BAM scanning produces
fragments in non-deterministic order"*. Fixed at the source, that reason is gone: units now arrive
canonical, so `unit_list` is already canonical. **Deleted, and proven a no-op** — the golden
`effective_length` reads `501.82302350785244` before and after, to the last digit.

⚠ The **class** sort is KEPT. It guards a different mechanism — `unordered_map` iteration order, which is
a property of the hash table, not of the scan — and ⛔ it has **no gate**: perturbation Y3 removes it and
nothing fails, because a single build's hash order is stable. Pre-existing; recorded, not fixed.

### ⭐ A SECOND DEFECT FELL OUT, and it was never filed

`fractional` mode was thread-count dependent too — it scatters float posteriors into shared accumulators,
so a permutation reorders the summation and the answer drifts by ULPs. Nobody had noticed because nobody
had compared *across* thread counts. The same ordering fixes it, and the gate covers all three modes.

### Cost

One `argsort` per locus over an int64 key. Measured: **0.15 s at 1 M fragments, 0.91 s at 5 M** — against
a ~66 s calibration. ⭐ An earlier draft densified the ids into a global rank first; it is exactly
equivalent (rank is monotone in `frag_id`) and cost **5×** more — 731 ms of that 906 ms — so it was
deleted too.

### Gates

| | | |
|---|---|---|
| **the contract** | one BAM, one seed, one answer — byte-identical across scan thread counts {1, 4, 16} × 2 repeats, in **all three** assignment modes | ✅ |
| **non-vacuity** | the buffer's row order really does differ between those thread counts | ✅ retried 3×, because it is a race |
| **the assumption** | no two EM units share a `frag_id` — otherwise the sort ties and falls back to buffer order | ✅ (2,118 unique ids from 3,000 fragments; a subset, so gapped but unique) |

⚠ **Two repeats per thread count, not one.** The permutation is random: a single pair can agree by luck,
and did — two 16-thread runs disagreed while 1 vs 4 vs 16 matched on the same tree.

### The perturbations

| | perturbation | result |
|---|---|---|
| **Y1** | revert to buffer-row order | ⭐ **fires** |
| **Y2** | order by CONTENT (`gdna_log_liks`) instead of identity | ⭐ **fires on both `sample` and `fractional`** — the wrong fix is caught |
| **Y3** | delete the equiv-CLASS sort as well | ⛔ **not caught** — pre-existing gap, recorded above |

⚠ Y1 fired via `fractional` rather than `sample` on that run. That is the honest shape of this gate:
`sample` is the mode that matters and the flaky detector; `fractional` is the reliable one, because ULP
drift is deterministic given a permutation. Together they hold.

### What this changes for the owner

⭐ **`assignment_mode="sample"` is now a legitimate A/B arm.** `CARRY_FORWARD.md` and `TODO.md` §7 both say
to A/B under `map`/`fractional` because the default wobbles by ~0.5 %. With a fixed seed all three modes
now reproduce exactly. ⚠ The three remain **different estimators** — 5441 / 6002 / 6277 on one scenario —
so an A/B must still hold the mode fixed across both arms.

⭐ **And it answers D-D ahead of S3.** The second pass needs no new seed machinery: S1 already gives the
deferred queue a canonical order, so `(global_seed, index in that queue)` is well defined for exactly the
reason established here. The rule is *order by identity, then one stream* — not *hash the content*.

---

## P0 — there was never a strand sign bug (2026-08-02)

    Spec: `docs/SPEC_SECOND_PASS.md` §3.3, P0   ·   Gate: `tests/calibration/test_strand_sense_convention.py`
    Retires: `TODO.md` rank 6   ·   Corrects: S5.f, S5.f-addendum, `CARRY_FORWARD.md` §1 fact 17 / §0 C4
    Suite: **21 failed / 1873 passed**, all 21 goldens.

`TODO.md` rank 6 recorded that the fitted κ is `1 − truth` and that *"only the exported scalar is
mis-labelled"*. ⛔ **Both statements are wrong, and the item had been filed twice.** The owner approved
fixing the sign; the audit found there is no sign to fix.

### The collision: two quantities, both called "strand specificity"

| | |
|---|---|
| `ReadSimConfig.strand_specificity` | *"probability an RNA fragment preserves correct read orientation … an R1↔R2 swap with probability 1 − ss"* — protocol **fidelity**, direction-agnostic |
| `StrandModel.p_r1_sense` | `P(align_strand == the junction's strand)` — **directional**. Its own docstring already said it: *"High (≈0.95) for R1-sense libraries (KAPA). Low (≈0.05) for R1-antisense (Illumina TruSeq dUTP)."* |

⭐ **For an R1-antisense protocol these are complements**, so comparing one against the other reads as a
sign error and is not one. The simulator emits R1-antisense — explicit in the gDNA path,
`r2_seqs, r1_seqs = _batch_extract_reads(...)`, the first extracted read becoming **R2** — which is the
most common real protocol. A dUTP library at 99 % fidelity genuinely has a sense fraction near 0.01.

### ⭐ THE MEASUREMENT — the tool already exposes the matching quantity

`StrandModel.strand_specificity = max(p_r1_sense, p_r1_antisense)` is direction-agnostic, like the
simulator's knob. Two genes, one per strand, 4,000 zero-gDNA fragments:

| simulated `ss` | `p_r1_sense` | `strand_specificity` | `rna_sense_frac` |
|---|---|---|---|
| 1.00 | 0.0000 | ⭐ **1.0000** | 0.0008 |
| 0.99 | 0.0156 | **0.9844** | 0.0164 |
| 0.75 | 0.2299 | **0.7701** | 0.2303 |
| 0.50 | 0.4980 | **0.5020** | 0.4980 |

⭐ It recovers the simulated parameter directly. The comparison the record should always have made.

⚠ **Opposite-strand genes are load-bearing in that fixture.** On a single-strand locus a convention error
that swapped the comparison's operands would move both together and be invisible.

### ⭐ And the 166× measurement now has an explanation rather than a rationalisation

S5.f-addendum forced κ to the nominal 0.99 and a zero-gDNA library read `f_gdna = 0.4992` against the
fitted value's 0.0030, and concluded *"the mirror cancels"*. ⛔ It does not cancel — **there is no mirror**.
`0.0101` is the right answer for an R1-antisense library, and 0.99 is a different quantity substituted for
it. The numbers stand; the explanation changes, and it is the simpler one.

### ⭐ What this unblocks

`SPEC_SECOND_PASS.md` §3.3's strand term, immediately and with no fix. `rna_sense_frac` is exactly the
`P(align_strand agrees | RNA)` the second pass needs to score an unspliced fragment's competing strand
hypotheses.

⚠ **But mind the direction.** On an R1-antisense library the value is ≈ 0.01, so the RNA hypothesis whose
implied strand *disagrees* with `align_strand` is the LIKELY one. ⛔ A scorer written as "agreement ⇒
multiply by `rna_sense_frac`" is exactly backwards on every real cfRNA library. Recorded in the spec.

### ⚠ A real gap found while closing a phantom one

**The simulator can only emit R1-ANTISENSE libraries.** `strand_specificity` is a swap probability about a
*fixed* orientation, never a choice of orientation, so **no simulated condition exercises the R1-sense
branch** — and real R1-sense libraries (KAPA-style) exist. Filed as `TODO.md` rank 6b and marked in place
by a **strict xfail** that deletes itself when the simulator gains the switch. ⚠ Low urgency: the branch is
a `max()` and a comparison, and real cfRNA is dUTP.

### Gates

| | | |
|---|---|---|
| **the recovery** | `StrandModel.strand_specificity` recovers the simulated knob at 1.00 / 0.75 / 0.50 | ✅ within 0.03 — ~4σ of the binomial standard error at 4,000 fragments, loose enough never to flake and far tighter than the 0.5-scale error a real flip would give |
| **the direction** | a perfectly stranded simulated library has `p_r1_sense ≈ 0`, and `read1_sense` reports the protocol as antisense | ✅ |
| **one convention** | `rna_sense_frac` is the Beta posterior mean of `p_r1_sense` and does not diverge from it | ✅ — a divergence would mean a second strand convention had appeared between the model and the balance fit |

⛔ **No code changed.** The deliverable is the gate and the corrected record: the absence of that gate is
what let one non-defect be filed twice, and it is the only thing that stops a third time.

---

## P1 — strand in the gap enumeration: two behaviours gated, one defect fixed (2026-08-02)

    Spec: `docs/SPEC_SECOND_PASS.md` §3.3 and D-5   ·   Gate: `tests/native/test_gap_hypothesis_strand.py`
    Suite: **21 failed / 1877 passed / 1 xfailed**, all 21 goldens.

The owner's 2026-08-02 ruling: *"any fragment with a splice junction, immediately, we can constrain the
gap [hypotheses] to that strand … fragments without a splice junction are unspliced and could be either
strand."* The audit found the first half **already implemented and tested nowhere**; this entry gates both
halves and fixes the one case that was wrong.

### ⭐ Gated: an observed junction pins the strand

`tP` (+) and `tM` (−) cover the same sequenced blocks and share the observed junction's coordinates, so
neither annotation nor overlap can separate them — only the sequenced **motif** can. Measured through a
real scan:

| observed XS motif | hypotheses | outcome |
|---|---|---|
| `+` | tP's 400 bp gap intron alone | deposits at **L = 450** |
| `−` | tM's 200 bp gap intron alone | deposits at **L = 650** |

⭐ The deposited `L` says which transcript was believed, with no need to read the hypothesis set. ⛔ And it
fires on the count too: without the pin both would survive and the fragment would be **deferred instead of
deposited**.

### ✅ Gated: an unspliced fragment offers both strands

Same locus, no motif: three hypotheses — `(1600,2000)` on `+`, `(1700,1900)` on `−`, and ∅ — with the two
spliced ones carrying **opposite** implied strands. ⚠ This is the case that makes the second pass's strand
term necessary rather than decorative: drop it and these two are separated by length alone.

### ⛔ Fixed: D-5 — one path claimed by both strands

`add_hypothesis` grouped paths by intron **coordinates only**, so two opposite-strand transcripts implying
the same intron merged into one hypothesis carrying the **first-seen** strand — an answer that flips when
two GTF lines are swapped. Grouping is right (it *is* one path); the strand is now set to **AMBIGUOUS**
when supporters disagree, which is what that value already means everywhere else and which `deposit`
already refuses to credit a junction on. Idempotent, three lines.

⚠ Unreachable on human data — **0 of 404,168** junction coordinates are annotated on both strands, and the
index warns that it is biologically impossible. Fixed because the alternative is order-dependent.

### ⚠ THE FIXTURE TRAP, recorded because it cost an hour

Driving `FragmentResolver` directly leaves `t_strand_arr_` empty, so **every hypothesis's implied strand
silently reads NONE** and a strand gate written that way passes against nothing. ⛔ All of these gates go
through the scan.

### The perturbations

| | perturbation | result |
|---|---|---|
| **Z1** | the merged hypothesis keeps the first supporter's strand | ⭐ **fires** — D-5 |
| **Z2** | the implied strand is dropped; every hypothesis reports NONE | ⭐ **fires on two gates** |
| **Z3** | `!certified_rna` removed from the ∅ condition | ⛔ **not caught** — and the reason is worth having |

### ⛔ Z3 DOES NOT FIRE, AND THAT IS NOT A WEAK FIXTURE

Measured: the `open` fragment's candidates are `t_inds = [0, 1, 4, 5]`, where **4 and 5 are the nRNA
shadows** (`RIGEL_NRNA_chr1_{1,2}_1000_2200`, `is_synthetic`). A shadow is single-exon, so it implies
nothing in the gap, so `any_candidate_implies_nothing` is already true and ∅ is emitted without the clause.

⭐ **That is `SPEC_GAP_PATHS.md` §2's own ruling arriving by a different route** — *"the nascent shadow IS
the ∅ hypothesis"* — and the two mechanisms agree rather than conflict:

* **spliced** fragment: the observed CIGAR-N intron falls inside the shadow's single exon, so the shadow
  cannot explain the read and drops out of `t_inds`. ∅ is correctly absent — the certified-RNA `mixed`
  fragment in `test_gap_introns_are_cut.py` has exactly one hypothesis.
* **unspliced** fragment: the shadow survives and supplies ∅, which is what §2's table requires.

⚠ So `!certified_rna` is **redundant, not wrong**, and it is KEPT: it states the rule directly rather than
depending on the shadow mechanism continuing to exist, and a single-exon gene has no separate shadow row
at all (`is_nrna` on a non-synthetic row means mature ≡ nascent). ⛔ Do not delete it on the strength of
this perturbation.

### Files touched

`src/rigel/native/resolve_context.h` (the D-5 fix) · `tests/native/test_gap_hypothesis_strand.py` (new)

---

## P2 — the scorer is BUILT and UNGATED (2026-08-02) ⛔ PARTIAL, READ THIS BEFORE USING IT

    Spec: `docs/SPEC_SECOND_PASS.md` §3 and P2   ·   Suite: **21 failed / 1877 passed / 1 xfailed**
    ⛔ **`src/rigel/second_pass.py` HAS NO TESTS.** It is inert — nothing imports it — but it is
    production code in `src/` that no gate covers, which is exactly the state where the next reader
    assumes it works. It does not yet have the right to be believed.

### What landed

| | |
|---|---|
| `Accumulator.length_under` | ⭐ `L` under ONE hypothesis without depositing. **Exposed, not reimplemented**: a Python scorer computing its own `L` would be a second definition of the quantity C0/C2 unified, and the one the drain would then disagree with. Verified: spliced 101, genomic 800, clipping correct |
| `MarshalledFragment` | the binding's hypothesis marshalling factored out, shared by `deposit` and `length_under`, so the two cannot drift about how a hypothesis set crosses the ABI |
| `src/rigel/second_pass.py` | the scorer: all three terms of §3, each kept separately on `HypothesisTerms` so a regression is attributable to one rather than to "the score moved" |

⚠ **D-1 is a parameter (`min` vs `geometric`) BECAUSE it must be measured.** It is not a tunable to be
left in: the docstring says it is expected to be deleted once the measurement is in. ⛔ If it is still
there in a month, that is the failure mode `CARRY_FORWARD.md` warns about.

### ⭐ THE SMOKE RUN ALREADY SHOWS D-3, AND IT IS NOT SUBTLE

Four-fragment fixture, both held records scored **uniform**:

```
record 0: [62100,63000)   L=502  rho=2.00e-03  f=0.000   -> undecided
record 1: [66100,67000)   L=500  rho=0.000     f=0.500   -> undecided
```

On record 1 **both junctions have zero flux** — the only fragment that would have used them is the held
one. That is exactly the systematic hole the spec predicted: *the held fragments are excluded from the
tally they are scored against*. On a toy it is 100 %. ⛔ **The real number is unknown and it is the fact
that decides D-3.**

### ⛔ WHAT IS DELIBERATELY NOT DONE, AND WHY

**The P2 gate is not written.** It says *"on a hand-built locus where the truth is known, the correct
hypothesis takes the larger share"* — and the available fixture has no depth, so every score is uniform
and the gate would pass or fail for reasons that have nothing to do with the scorer. ⚠ Writing it against
this fixture would have produced a green gate over an undecided scorer, which is worse than no gate.

Three things remain, and the order matters:

1. ⭐ **Measure D-3 on the pilot first** — pure measurement, data already on disk
   (`~/Downloads/rigel_runs/suite/pilot`, 8 conditions with `truth_fragment_lengths.tsv`). It needs no new
   fixture, and its answer may change what the gate should even assert. ⛔ If the zero-flux share is large
   the score needs a fallback, and **that is an owner decision — no pseudocount will be invented**.
2. **A discriminating fixture, then the gate.** Deep enough that rho and f actually separate.
3. **Measure D-1**, then delete the parameter.
