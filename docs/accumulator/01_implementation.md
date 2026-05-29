# Accumulator v2 — Phased Implementation Plan

Companion to [`00_design.md`](00_design.md). This document is the
step-by-step playbook. Each phase has explicit file lists,
acceptance gates, and rollback plans. **No phase begins until the
previous phase's acceptance gate is green and is reviewed.**

---

## Cross-cutting conventions

- **Build between every phase.** `pip install --no-build-isolation -e .`
  must succeed before any test command runs.
- **Tests scoped per phase.** Each phase adds its own test file
  under `tests/native/` or `tests/` and is mergeable as a stand-alone
  commit. Never delete a legacy test until the parity gate (Phase 5)
  signs off.
- **No premature deletions.** Legacy code (`RegionCountLedger`
  boundary fields, [`boundaries.py`](../../src/rigel/calibration/boundaries.py))
  remains compiled-and-callable through Phase 6. Deletion happens
  only in Phase 7.
- **FL is not accumulated.** Per audit_phase1.md decision #6, the
  v2 substrate does not collect per-region or per-boundary FL
  histograms. FL stays in the downstream EM scorer
  ([`scoring.py`](../../src/rigel/scoring.py),
  [`frag_length_model.py`](../../src/rigel/frag_length_model.py))
  where it can be applied per-fragment without binning.
- **Naming.** Old type: `RegionCountLedger`. New type:
  `Accumulator`. Old Python adapter: `RegionCountsPayload`
  (whatever it's called today). New: `AccumulatorPayload`. The two
  coexist until Phase 7.

---

## Phase 1 — Audit & spec lock

**Goal.** Document the *exact* current per-fragment accounting
semantics in `src/rigel/native/accumulator.cpp` and
`bam_scanner.cpp`. Fix any per-block double-counting bugs found
**inside the legacy implementation** so the parity gate in Phase 5
compares two correct accumulators. Lock the v2 spec by writing the
reference test corpus.

### 1.1 Code review tasks (read-only)

Files to read top-to-bottom and annotate with a written memo
(`docs/accumulator/audit_phase1.md`):

1. `src/rigel/native/accumulator.h`
2. `src/rigel/native/accumulator.cpp`
3. `src/rigel/native/bam_scanner.cpp` — focus on the per-fragment
   loop that calls accumulator deposit helpers.
4. `src/rigel/native/region_count_ledger.{h,cpp}` (or equivalent)
5. `src/rigel/calibration/boundaries.py`
6. Tests that touch the accumulator:
   `tests/test_region_count_ledger.py`,
   `tests/test_region_partition.py`,
   `tests/test_region_persist.py`,
   `tests/test_boundary_model.py`,
   `tests/test_boundary_sweep.py`.

### 1.2 Memo contents (`audit_phase1.md`)

Required sections:

1. **Per-fragment vs per-block.** Quote each line in
   `accumulator.cpp` where a deposit occurs. State whether it is
   inside a per-block loop or a per-fragment loop, and whether the
   normalization divides by per-block $\ell_k$ or per-fragment $L$.
2. **Boundary-mass-mirroring inventory.** Today's
   `RegionCountLedger` stores `boundary_left_*` on region $r$ and
   `boundary_right_*` on region $r-1$ for the same boundary. Confirm
   that the two slots are never both incremented independently.
3. **Fully-spans-region path.** Document the current logic for
   fragments where one block straddles a region. Today this triggers
   the `w/2` split rule (see [`00_design.md`](00_design.md) §2.3).
4. **Spliced-flag attribution.** State exactly when a fragment is
   marked spliced vs unspliced in the boundary deposit. Compare to
   the v2 spec (per-junction-event spliced flag).
5. **Strand attribution.** State how strand $s \in \{+, -\}$ is set
   per fragment.
6. **Found bugs.** Enumerate any discrepancies between the
   documented intent and the code. Each bug gets a ticket-style
   entry: file:line, observed behavior, expected behavior, fix
   proposal.

### 1.3 Fix any bugs found in legacy code

If §1.2.6 surfaces bugs, fix them **inside the legacy
`accumulator.cpp`** with their own commits and tests in this phase
(not v2). The point is to make the parity gate in Phase 5 meaningful.

### 1.4 Lock the v2 spec via reference test corpus

Create `tests/native/test_accumulator_spec.py` with **expected
outputs only** (the v2 implementation does not exist yet; tests are
xfail-marked until Phase 2 lands). Each test is a tiny synthetic
fragment list paired with the expected `AccumulatorPayload`
contents, computed by hand or by a Python reference implementation
in the test file itself.

Tests required (from [`00_design.md`](00_design.md) §9):

| ID | Name |
|---|---|
| T1 | `test_contained_single_block` |
| T2 | `test_contained_multi_block_spliced` |
| T3 | `test_two_block_adjacent_regions` |
| T4 | `test_two_block_non_adjacent_regions` (user-verified §4.5.1 walkthrough) |
| T5 | `test_three_block_all_adjacent` (user §4.5.2 walkthrough) |
| T6 | `test_fully_spans_region` |
| T7 | `test_mass_conservation_random` |
| T8 | `test_flux_dedup_adjacent_regions` |
| T10 | `test_strand_attribution_negative_strand` |
| T11 | `test_spliced_flag_per_junction` |

Plus a pure-Python reference implementation in
`tests/native/_accumulator_reference.py` (~150 LoC) that any
human can read and confirm matches the spec. Phase 2's native
implementation will be tested against this Python reference.

### 1.5 Acceptance gate

- [ ] `docs/accumulator/audit_phase1.md` exists with all 6 required sections.
- [ ] All bugs found in §1.2.6 either fixed (with tests) or filed as deferred (with justification).
- [ ] `tests/native/test_accumulator_spec.py` exists with all 11 tests, currently `pytest.mark.xfail(reason="accumulator not yet implemented")`.
- [ ] `tests/native/_accumulator_reference.py` is a working Python implementation.
- [ ] `pytest tests/ -v` passes (legacy tests still green; v2 tests xfail).

### 1.6 Rollback

This phase is read-only / docs / fix-legacy-bugs. Rollback = revert
the (small) legacy-bug-fix commits if they regress something.

---

## Phase 2 — Native rewrite (`accumulator.{h,cpp}`)

**Goal.** Land the C++ implementation. Compile alongside legacy
accumulator (both are live in the binary). No call sites switched
yet.

### 2.1 Files added

1. `src/rigel/native/accumulator.h` — full struct + class API per
   [`00_design.md`](00_design.md) §3, §7.1.
2. `src/rigel/native/accumulator.cpp` — implementation.

### 2.2 Files modified

- `src/rigel/native/CMakeLists.txt` — add the new `.cpp` to
  `_bam_impl` target sources. Bind module unchanged in this phase.

### 2.3 Implementation outline (`accumulator.cpp`)

```cpp
void Accumulator::deposit(const FragmentEvent& fe) {
  // fe carries: ref_id, blocks[K] (start, end), strand_sign,
  //             per_block_region_id[K], per_junction_spliced[K-1].

  // 1. Compute L = sum(block lengths).
  // 2. Classify contained vs crossing:
  //      crossing = any (per_block_region_id[k] != per_block_region_id[0])
  //                 OR any block straddles a boundary (rare; §4.5.4).

  if (!crossing) {
    int r = per_block_region_id[0];
    bool spliced_any = std::any_of(...);
    // Increment contained channel.
    if (spliced_any) {
      if (strand_pos) ++regions_[r].contained_spl_pos;
      else            ++regions_[r].contained_spl_neg;
    } else {
      if (strand_pos) ++regions_[r].contained_unspl_pos;
      else            ++regions_[r].contained_unspl_neg;
    }
    return;
  }

  // Crossing path.
  // For each junction (k, k+1):
  std::array<int, MAX_JUNCTIONS> touched_boundaries{};
  int n_touched = 0;
  for (int k = 0; k + 1 < K; ++k) {
    int Rk   = per_block_region_id[k];
    int Rkp1 = per_block_region_id[k+1];
    int b_out = right_boundary_of(Rk);        // = ref_boundary_offset + local_region_index_k + 1
    int b_in  = left_boundary_of(Rkp1);       // = ref_boundary_offset + local_region_index_kp1
    bool spliced = per_junction_spliced[k];
    float mass_k   = float(block_len[k])   / float(L);
    float mass_kp1 = float(block_len[k+1]) / float(L);

    // block-on-left (B_k is left of b_out) -> mass_left
    deposit_left (b_out, strand_pos, spliced, mass_k);
    // block-on-right (B_{k+1} is right of b_in) -> mass_right
    deposit_right(b_in,  strand_pos, spliced, mass_kp1);

    // Flux dedup: each touched boundary gets one flux per fragment-junction.
    flux_increment_dedup(touched_boundaries, n_touched, b_out, strand_pos, spliced);
    flux_increment_dedup(touched_boundaries, n_touched, b_in,  strand_pos, spliced);
  }
}
```

`flux_increment_dedup` checks whether the boundary index is already
in `touched_boundaries[0..n_touched]`; if not, appends and
increments. The "channel" granularity (spl_pos/spl_neg/etc.) is
finer than the dedup key, but per [`00_design.md`](00_design.md)
§4.4 the dedup is at boundary level — a single junction event over
adjacent regions is one observation, one flux increment.

### 2.4 The fully-spans-region case (§4.5.4)

Block $B$ straddles boundaries $B_a, B_b, \ldots$. Decompose $B$
into per-region intersections $\ell^{(R)}$ for each region $R$ the
block overlaps. Treat the resulting sequence of $(R, \ell^{(R)})$
slices **as if they were adjacent blocks with no inter-block gap**
and apply the §4.3 rule. Mark each implicit boundary-crossing as
unspliced (it's a contiguous CIGAR M run, not an N).

Implementation: a helper `expand_block_to_region_slices(block,
out)` runs before the main per-block loop. After expansion, the
deposit loop is identical to the multi-block case.

### 2.5 Tests added

- `tests/native/test_accumulator_spec.py` — remove `xfail` markers.
- `tests/native/test_accumulator_against_python_reference.py` —
  property tests: 1000 random synthetic fragment streams, assert
  the native `Accumulator` output matches the Python reference
  byte-for-byte (after the float32→bit equality dance).

### 2.6 Acceptance gate

- [ ] `pip install --no-build-isolation -e .` succeeds.
- [ ] All 11 spec tests pass (no `xfail`).
- [ ] Python-reference property tests pass on 1000 random fragments.
- [ ] Legacy tests still pass (we haven't switched callers yet).

### 2.7 Rollback

Drop the new source files from CMakeLists.txt. Nothing else is
touched.

---

## Phase 3 — nanobind binding + Python adapter

**Goal.** Expose `Accumulator` to Python via the existing
`_bam_impl` nanobind module. Add `AccumulatorPayload` dataclass.
Still no call sites switched.

### 3.1 Files added

1. `src/rigel/scan_payload.py` — augmented (or new module if missing)
   with the `AccumulatorPayload` dataclass from
   [`00_design.md`](00_design.md) §7.3.

### 3.2 Files modified

- `src/rigel/native/bam_scanner.cpp` — nanobind binding block at the
  bottom: expose `Accumulator::regions`, `boundaries`, and offset
  arrays as numpy arrays (zero-copy with `nb::ndarray`).
- `src/rigel/native.py` — re-export the new bindings.

### 3.3 Numpy view contract

Each of the 16 boundary mass/flux fields and 4 region contained
fields is exposed as a separate 1-D numpy array. Lifetime is tied
to the `Accumulator` C++ object via nanobind's owner-reference
machinery.

### 3.4 Tests added

`tests/test_accumulator_payload.py`:

- Construct a tiny synthetic accumulator from Python (using a thin
  Python factory exposed by the binding).
- Confirm field dtypes (`uint32`, `float32`), shapes, and contiguity.
- Confirm zero-copy: writes to the numpy view mutate the C++ state
  (or, if we choose to expose read-only views, writes raise).
- Confirm garbage collection: the numpy arrays keep the underlying
  C++ object alive until the last reference is dropped.

### 3.5 Acceptance gate

- [ ] Build succeeds.
- [ ] `pytest tests/test_accumulator_payload.py -v` passes.
- [ ] No regression in existing tests.

### 3.6 Rollback

Revert the nanobind binding block in `bam_scanner.cpp` and remove
`AccumulatorPayload` from `scan_payload.py`.

---

## Phase 4 — Wire v2 into `bam_scanner.cpp` (dual-write)

**Goal.** While scanning a BAM, deposit each fragment into **both**
the legacy `RegionCountLedger` and the new `Accumulator`. Emit
both payloads. No downstream consumer changes yet.

### 4.1 Files modified

- `src/rigel/native/bam_scanner.cpp` — in the per-fragment finalize
  block, after calling the legacy deposit, call
  `accumulator_.deposit(fe)`. Wrap the v2 call in a runtime
  on/off flag (env var `RIGEL_ACCUMULATOR=1` default-on for
  developer builds; flip the default after Phase 5 parity passes).
- `src/rigel/scan.py` — surface the `AccumulatorPayload` alongside
  whatever payload `scan_and_buffer` already returns.

### 4.2 Tests added

`tests/test_dual_write_parity.py`:

- Run the scanner on one scenario BAM with both accumulators enabled.
- Project the v2 boundary view back to region-keyed `boundary_left` /
  `boundary_right` channels using the indexing in
  [`00_design.md`](00_design.md) §3.1.
- Assert the projection equals the legacy `RegionCountLedger`'s
  boundary slots within float32 tolerance (`abs_tol = 1e-5`).
- Assert contained channels match exactly (integer compare).

### 4.3 Acceptance gate

- [ ] Build succeeds.
- [ ] `pytest tests/test_dual_write_parity.py -v` passes on at least
      3 scenario BAMs (`tests/scenarios_aligned/`).
- [ ] No regression in any existing test.
- [ ] Scan-time overhead measured: report `time pytest tests/test_pipeline_smoke.py` before and after.

### 4.4 Rollback

Flip `RIGEL_ACCUMULATOR=0`. The legacy path is fully intact.

---

## Phase 5 — Parity gate against real data

**Goal.** Run the full benchmark suite with dual-write enabled and
prove the v2 substrate matches legacy on real-scale BAMs. This is
the single critical-path gate of the project.

### 5.1 Tasks

1. Run `python -m scripts.benchmarking run -c
   scripts/benchmarking/configs/default.yaml --dry-run` first to
   confirm conditions are discoverable.
2. Re-run rigel on 3 representative conditions:
   - `gdna_none_ss_1.00_nrna_none`
   - `gdna_high_ss_0.50_nrna_rand`
   - `gdna_low_ss_0.90_nrna_none`
   Use the new dual-write build. Output goes to
   `runs/<condition>/rigel/accumulator_phase5/`.
3. Compare to the existing golden `runs/<condition>/rigel/default/`
   transcript-level quant.feather. Allowed: floating-point noise at
   ~1e-6 relative tolerance.

### 5.2 Diagnostics (`scripts/debug/accumulator_parity.py`)

A new script that:
- Loads both `AccumulatorPayload` and the legacy `RegionCountLedger`
  payload from a single dual-write scan.
- Computes per-boundary, per-region, per-channel diffs.
- Reports the top-N divergent boundaries with synthetic-context
  metadata (which transcripts overlap, what the locus shape is).

### 5.3 Acceptance gate

- [ ] All 3 conditions match golden quant.feather within 1e-6 relative tolerance for `mrna_abundance`.
- [ ] Per-boundary parity report shows zero diffs above 1e-5 absolute.
- [ ] If any diff exceeds threshold: do NOT proceed. Root-cause first.
       Most likely culprits: implicit-splice classification mismatch
       (§4.3 spliced flag), or fully-spans-region attribution
       (§4.5.4). File a bug, fix in Phase 2, return here.

### 5.4 Rollback

If the gate fails and the cause is structural (i.e. the v2 design
itself is wrong), pause and re-open the design phase. Otherwise:
fix in the appropriate phase and re-run.

---

## Phase 6 — Consumer migration

**Goal.** Switch every downstream consumer from `RegionCountLedger`
+ `BoundaryTable` to `AccumulatorPayload`. Legacy code still
compiled but no longer read.

### 6.1 Consumers to migrate (audit before starting)

Run `grep_search` for these symbols and migrate each callsite:

1. `BoundaryTable` — primary consumer in
   `src/rigel/calibration/`.
2. `build_boundary_table`, `validate_boundary_table` — from
   `src/rigel/calibration/boundaries.py`.
3. `RegionCountLedger` boundary-slot reads (anything reading
   `boundary_left_*` or `boundary_right_*` off the ledger).
4. `compute_locus_priors_from_partitions` — confirm whether it reads
   from boundary table or directly from ledger.
5. Any test fixtures that construct synthetic `BoundaryTable`s
   (they get replaced with synthetic `AccumulatorPayload`s; reuse
   the Python reference impl from Phase 1).

### 6.2 Files modified (expected list — confirm via grep before edits)

- `src/rigel/calibration/_simple.py`
- `src/rigel/calibration/_categorize.py`
- `src/rigel/locus.py` (if it reads boundary fields)
- `src/rigel/priors.py` (if it reads boundary fields)
- Test fixtures: `tests/test_boundary_model.py`, `tests/test_boundary_sweep.py`

### 6.3 Strategy

One commit per consumer migration. Each commit:
1. Switches one consumer from `BoundaryTable` to
   `AccumulatorPayload`.
2. Updates that consumer's tests.
3. Runs full `pytest tests/ -v` to confirm no other consumer broke.

### 6.4 Acceptance gate

- [ ] No remaining import of `boundaries.py` outside `boundaries.py` itself.
- [ ] No remaining read of `RegionCountLedger.boundary_left_*` or `.boundary_right_*` fields.
- [ ] Full `pytest tests/ -v` passes.
- [ ] Re-run Phase 5 parity gate: still passes (consumer migration must not regress numerics).

### 6.5 Rollback

Revert the offending consumer-migration commit. The parallel dual-write
infrastructure means the legacy path is still alive and callable.

---

## Phase 7 — Legacy deletion

**Goal.** Remove the legacy boundary path. The legacy contained
counts may stay if `RegionCountLedger` still serves other purposes;
otherwise, delete the whole ledger.

### 7.1 Files deleted

- `src/rigel/calibration/boundaries.py`
- `tests/test_boundary_model.py`, `tests/test_boundary_sweep.py`
  (or rewrite as tests of `AccumulatorPayload` if their coverage is
  still valuable).
- `RegionCountLedger::boundary_left_*` / `::boundary_right_*` fields
  in `src/rigel/native/region_count_ledger.h`.
- Dual-write call in `bam_scanner.cpp` (delete the legacy deposit
  call; keep only the v2 deposit).
- `RIGEL_ACCUMULATOR` env-var flag.

### 7.2 Files modified

- `src/rigel/native/region_count_ledger.cpp` — strip boundary
  accounting helpers.
- `tests/test_region_count_ledger.py` — strip boundary tests.

### 7.3 Documentation updates

- `CLAUDE.md` — update the "Python Module Roles" table: remove
  `boundaries.py` row, update `scan_payload.py` description.
- `.github/copilot-instructions.md` — same.
- `docs/caljointmodel/04_interface_contract.md` §2 — remove the
  "Status (2026-05-29)" pending banner; mark substrate as live.
- `CHANGELOG.md` — entry for the rewrite.

### 7.4 Acceptance gate

- [ ] Full `pytest tests/ -v` passes.
- [ ] `ruff check src/ tests/` clean.
- [ ] Benchmark re-run on 3 conditions still matches golden.
- [ ] No grep hits for `BoundaryTable`, `boundaries.py`,
      `boundary_left_`, `boundary_right_`.

### 7.5 Rollback

Revert the deletion commit if any test surface lights up. The
dual-write deletion is the most likely source of regression.

---

## Phase 8 — Calibration v6 substrate hookup

**Goal.** Re-open [`../caljointmodel/`](../caljointmodel/) and wire
the calibrator to consume `AccumulatorPayload` directly per docs
01–03. This phase is *out of scope* for the accumulator v2 plan
itself but is the immediate downstream consumer; flag it here so
the v2 work doesn't drift away from its purpose.

### 8.1 Trigger

Phase 7 acceptance gate green.

### 8.2 First task

Read [`../caljointmodel/05_implementation_plan.md`](../caljointmodel/05_implementation_plan.md)
Phase 3 ("substrate adapter") and re-scope: the substrate is now
the live `AccumulatorPayload`, not an extension. Update that doc to
mark it as completed-by-accumulator-v2.

---

## Phase 9 — Cleanup & perf check

**Goal.** Final polish. Optional perf measurements to confirm v2 is
no slower than legacy (it should be slightly faster: no
mass-mirroring across two regions).

### 9.1 Tasks

- Profile a full scan with `scripts/profiler.py` on the largest
  available scenario BAM.
- Report wall-clock and memory delta vs the Phase 5 measurement.

### 9.2 Acceptance gate

- [ ] Scan-time within ±5% of legacy baseline.
- [ ] Peak RSS within ±10% of legacy baseline.

---

## Risk register summary

| Phase | Top risk | Mitigation |
|---|---|---|
| 1 | Audit memo skipped | Required artifact at gate |
| 2 | Fully-spans-region case is subtle | Pre-decomposition helper + dedicated test T6 |
| 3 | Nanobind lifetime bugs | Explicit owner-reference + GC test |
| 4 | Dual-write doubles scan time | Profile in 4.3 gate; if >2× slowdown, redesign |
| 5 | Numerical divergence on real data | Diagnostics script + per-boundary report |
| 6 | Hidden consumer found mid-migration | Grep-driven audit before any edits |
| 7 | Deletion regresses unrelated tests | Revert commit; tests are rebuildable |
| 8 | Out of scope | Documented as future work |
| 9 | Perf regression | Profile then optimize |

---

## File-touch matrix

| File | P1 | P2 | P3 | P4 | P5 | P6 | P7 |
|---|---|---|---|---|---|---|---|
| `accumulator.{h,cpp}` (legacy) | read+fix | — | — | — | — | — | delete boundary fields |
| `accumulator.{h,cpp}` | — | **add** | — | — | — | — | — |
| `bam_scanner.cpp` | read | edit (CMake) | edit (binding) | edit (dual-write) | — | — | edit (drop legacy call) |
| `scan.py` | — | — | — | edit | — | — | — |
| `scan_payload.py` | — | — | **add** | edit | — | — | — |
| `native.py` | — | — | edit | — | — | — | — |
| `calibration/boundaries.py` | read | — | — | — | — | — | **delete** |
| `calibration/_simple.py` | read | — | — | — | — | edit | — |
| `calibration/_categorize.py` | read | — | — | — | — | edit | — |
| `locus.py` | read | — | — | — | — | edit if needed | — |
| `priors.py` | read | — | — | — | — | edit if needed | — |
| `tests/native/test_accumulator_spec.py` | **add (xfail)** | unxfail | — | — | — | — | — |
| `tests/native/_accumulator_reference.py` | **add** | — | — | — | — | — | — |
| `tests/test_accumulator_payload.py` | — | — | **add** | — | — | — | — |
| `tests/test_dual_write_parity.py` | — | — | — | **add** | — | — | delete |
| `tests/test_boundary_model.py` | read | — | — | — | — | rewrite or — | delete |
| `tests/test_boundary_sweep.py` | read | — | — | — | — | rewrite or — | delete |
| `scripts/debug/accumulator_parity.py` | — | — | — | — | **add** | — | — |
| `docs/accumulator/audit_phase1.md` | **add** | — | — | — | — | — | — |
| `docs/caljointmodel/04_interface_contract.md` | — | — | — | — | — | — | edit (drop banner) |
| `CLAUDE.md` | — | — | — | — | — | — | edit |
| `.github/copilot-instructions.md` | — | — | — | — | — | — | edit |
| `CHANGELOG.md` | — | — | — | — | — | — | edit |

---

## Sign-off

This plan covers Phase 1 through Phase 9. Approval to begin Phase
1 = approval of the plan. Each phase gate is its own checkpoint.
