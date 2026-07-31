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
| **B2** | the pilot suite, and B0's verdict on it | below |
| **B3** | the scan cache — scan once, calibrate many times | below |

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
