# hulkrna Test Suite Audit

**Date:** 2025-01-XX  
**Scope:** All 21 files in `tests/` (including `conftest.py`)  
**Source modules:** 17 in `src/hulkrna/` + `sim/` subpackage (4 modules)

---

## Executive Summary

| Metric | Value |
|--------|-------|
| Test files | 21 (incl. conftest) |
| Approx. test cases (after parametrization) | ~600 |
| Lines of test code | ~8,800 |
| Stale / broken imports | 0 |
| Duplicate helper functions | 1 pair (`_make_locus_em_data` in 2 files) |
| Source modules without dedicated test file | 4 (`locus.py`, `scan.py`, `scoring.py`, `transcript.py`) |
| Primary trim candidate | `test_scenarios.py` (~219 parametrized end-to-end tests) |

### Top Recommendations

1. **Trim `test_scenarios.py`** – reduce ~219 parametrized tests to ~80 by cutting redundant sweep levels and consolidating similar scenarios.
2. **Merge `test_diagnostic.py` into `test_scenarios.py`** – they test the same thing (simulation → pipeline) with identical imports and patterns.
3. **Extract shared helper `_make_locus_em_data`** into `conftest.py` – duplicated verbatim in `test_gdna.py` and `test_estimator.py`.
4. **Redirect private-API tests in `test_gdna.py`** – imports `_score_gdna_candidate` / `_GDNA_SPLICE_PENALTIES` (aliases that may break on refactor).
5. **Expand `conftest.py` fixtures** – current fixtures (`mini_gtf_file`, `mini_fasta_file`) are used by only 2 files; shared mocks could serve 4+ files.

---

## Per-File Analysis

### 1. `conftest.py` (75 lines)

| Attribute | Detail |
|-----------|--------|
| Purpose | Shared fixtures: `MINI_GTF` string, `mini_gtf_file`, `mini_fasta_file` |
| Test count | 0 (fixture-only) |
| Consumers | `test_index.py`, `test_gtf.py` |
| Staleness | None |
| Quality | Clean but under-utilized — many tests build their own identical fixtures |

**Recommendation: Expand.** Add shared helpers (`_make_locus_em_data`, mock `HulkIndex`, `FragmentLengthModels.uniform()`) here to DRY up `test_gdna.py`, `test_estimator.py`, `test_pipeline_routing.py`, and `test_buffer.py`.

---

### 2. `test_cli.py` (26 lines, 2 tests)

| Attribute | Detail |
|-----------|--------|
| Tests | `test_index_requires_args`, `test_count_requires_args` |
| Imports | `hulkrna.cli.build_parser` |
| Staleness | None |
| Quality | Minimal — only checks argparse raises on missing args |

**Recommendation: Keep.** Consider adding a test for `get_version()` and a smoke test for `main()` with `--help`.

---

### 3. `test_gtf.py` (50 lines, 4 tests)

| Attribute | Detail |
|-----------|--------|
| Tests | `test_strict_valid`, `test_strict_invalid_raises`, `test_warn_skip_valid`, `test_warn_skip_drops_bad` |
| Imports | `hulkrna.gtf.GTF`, `hulkrna.transcript.Transcript` |
| Staleness | None |
| Quality | Good — focused, uses `conftest` fixtures |

**Recommendation: Keep as-is.**

---

### 4. `test_stats.py` (84 lines, ~11 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestPipelineStats`, `TestBamStatsProxy` |
| Imports | `hulkrna.stats.PipelineStats`, `_BamStatsProxy` |
| Staleness | None |
| Quality | Good — covers initialization, accumulation, repr, `_BamStatsProxy` summation |

**Recommendation: Keep as-is.**

---

### 5. `test_categories.py` (90 lines, ~9 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestSpliceType`, `TestSpliceStrandCol`, `TestColumnSubsets` |
| Imports | `hulkrna.categories.SpliceType`, `SpliceStrandCol`, `NUM_SPLICE_STRAND_COLS` |
| Staleness | None |
| Quality | Good — validates enum values, column indexing, `NUM_SPLICE_STRAND_COLS` invariant |

**Recommendation: Keep as-is.**

---

### 6. `test_index.py` (98 lines, 2 tests)

| Attribute | Detail |
|-----------|--------|
| Tests | `test_load_roundtrip`, `test_load_validates_version` |
| Imports | `hulkrna.index.HulkIndex`, `load_reference_lengths`, `read_transcripts`, etc. |
| Staleness | None |
| Quality | Good — tests serialization round-trip and version validation |

**Recommendation: Keep.** Could add tests for `build_splice_junctions`, `build_genomic_intervals`, `gtf_to_bed12`.

---

### 7. `test_types.py` (244 lines, ~30 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestStrand` (6), `TestInterval` (5), `TestGenomicInterval` (3), `TestIntervalType`, `TestRefInterval` (3), `TestMergeCriteria` (4), `TestMergeResult` (4), `TestEmptyMerge` (2) |
| Imports | All from `hulkrna.types` |
| Staleness | None |
| Quality | Excellent — thorough enum, dataclass, and edge-case coverage |

**Recommendation: Keep as-is.** Model test file for the project.

---

### 8. `test_strand_model.py` (222 lines, ~20 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestObserve` (4), `TestPosterior` (5), `TestLikelihood` (3), `TestSerialization` (3), `TestProperties` (2), `TestFallback` (3) |
| Imports | `hulkrna.strand_model.StrandModel/StrandModels`, `hulkrna.types.Strand` |
| Staleness | None |
| Quality | Excellent — Bayesian correctness, serialization, edge cases |

**Recommendation: Keep as-is.**

---

### 9. `test_annotate.py` (249 lines, ~8 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestAnnotationTable` (unit, ~4 tests), `TestAnnotatedBamIntegration` (integration, ~4 tests) |
| Imports | `hulkrna.annotate.AnnotationTable`, `write_annotated_bam` |
| External deps | minimap2, samtools (integration tests) |
| Staleness | None |
| Quality | Good — unit and integration nicely separated |

**Recommendation: Keep as-is.** Mark integration tests with `@pytest.mark.integration` if not already.

---

### 10. `test_frag_length_model.py` (375 lines, ~30 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestBasic`, `TestStatistics`, `TestLikelihood`, `TestSerialization`, `TestModelsContainer`, `TestTranscriptEffectiveLength`, `TestEffectiveLengthCaching` |
| Imports | `hulkrna.frag_length_model.FragmentLengthModel/FragmentLengthModels`, `hulkrna.categories.SpliceType` |
| Staleness | None |
| Quality | Excellent — statistical validation, caching, effective-length calculations |

**Recommendation: Keep as-is.**

---

### 11. `test_pipeline_routing.py` (393 lines, 6 tests)

| Attribute | Detail |
|-----------|--------|
| Tests | 6 functions testing `_scan_and_build_em_data` |
| Imports | `hulkrna.pipeline._scan_and_build_em_data` (private), `hulkrna.buffer.*`, `hulkrna.estimator.*`, etc. |
| Staleness | None (verified `_scan_and_build_em_data` exists at pipeline.py:589) |
| Quality | Good — uses mocks for `HulkIndex`/`FragmentBuffer`/`BufferedFragment` to test routing logic in isolation |
| Issue | Builds its own mock classes that partially overlap with `test_buffer.py` mocks |

**Recommendation: Keep.** Consider moving shared mock classes (`_MockIndex`, `_MockBuffer`) into `conftest.py`.

---

### 12. `test_fragment.py` (400 lines, ~25 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestFromReads` (5), `TestBasics` (5), `TestMultiTag` (5), `TestSupplementaryMerging` (5), `TestNMExtraction` (5) |
| Imports | `hulkrna.fragment.Fragment`, `hulkrna.types.GenomicInterval/Strand` |
| Staleness | None |
| Quality | Good — uses mock `pysam.AlignedSegment` objects extensively |

**Recommendation: Keep as-is.**

---

### 13. `test_sim.py` (538 lines, ~30 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestMutableGenome`, `TestGeneBuilder`, `TestReadSimulator`, `TestGeneBuilderHulkIndex` |
| External deps | minimap2, samtools, pysam |
| Staleness | None |
| Quality | Good — validates the simulation infrastructure that `test_scenarios.py` and `test_diagnostic.py` depend on |

**Recommendation: Keep as-is.** This validates the test infrastructure itself.

---

### 14. `test_diagnostic.py` (551 lines, ~9 tests) ⚠️

| Attribute | Detail |
|-----------|--------|
| Classes | `TestIsoformCollapse` (2), `TestUnsplicedAnnihilation` (2), `TestGDNAOverAbsorption` (2), `TestAntisenseDrift` (3) |
| External deps | minimap2, samtools |
| Imports | `hulkrna.pipeline.run_pipeline`, `hulkrna.sim.Scenario/SimConfig/GDNAConfig/run_benchmark` |
| Staleness | None |
| Redundancy | **High** — uses identical simulation→pipeline pattern as `test_scenarios.py` |
| Quality | Well-motivated (each class targets a specific known bug), but the scenarios overlap with `test_scenarios.py` classes |

**Recommendation: Merge into `test_scenarios.py` as additional scenario classes.** The "diagnostic" framing is valuable but doesn't justify a separate file — these are just named regression scenarios. After merging, mark them with `@pytest.mark.regression`.

---

### 15. `test_oracle_bam.py` (599 lines, ~25 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestCoordinateProjection` (4), `TestCigar` (3), `TestTakeFromLeft` (3), `TestTakeFromRight` (3), `TestOracleBamBasic` (3), `TestOracleBamGDNA` (2), `TestOracleBamNRNA`, `TestOracleBamStrand`, `TestOracleBamCoordSorted`, `TestOracleBamTruthCounts` |
| External deps | minimap2, samtools, pysam |
| Staleness | None |
| Quality | Good — unit tests for coordinate math + integration tests for BAM simulation |

**Recommendation: Keep as-is.** Tests sim infrastructure.

---

### 16. `test_bam.py` (769 lines, ~35 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestParseRead` (8), `TestPairReads` (5), `TestGroupRecordsByHit` (6), `TestParseBamFile` (6), `TestComplexReadGroups` (3) |
| External deps | pysam |
| Staleness | None |
| Quality | Very thorough — covers SAM/BAM parsing edge cases, strand tags, multi-mapping, chimeras |

**Recommendation: Keep as-is.**

---

### 17. `test_buffer.py` (775 lines, ~40 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestBufferedFragment`, `TestFragId`, `TestFragmentBufferBasic`, `TestFragmentClasses`, `TestAmbiguousMask`, `TestDiskSpill`, `TestBufferSummary`, `TestIterChunks`, `TestChimeraHandling` |
| Imports | `hulkrna.buffer.*`, `hulkrna.resolution.ResolvedFragment`, `hulkrna.types.*`, `hulkrna.categories.SpliceType` |
| Staleness | None |
| Quality | Excellent — covers disk spill, chunk iteration, ambiguity masking, chimera handling |

**Recommendation: Keep as-is.**

---

### 18. `test_gdna.py` (877 lines, ~35 tests) ⚠️

| Attribute | Detail |
|-----------|--------|
| Classes | `TestScoreGDNACandidate` (6), `TestLocusGDNATheta` (3), `TestGDNAAssignment`, `TestGDNAAttribution`, `TestGDNAProperties`, `TestGDNAOutput`, `TestGDNAStats`, `TestGDNARate`, `TestNRNAInit`, `TestLocusBehavior` |
| Imports | `hulkrna.pipeline._score_gdna_candidate`, `_GDNA_SPLICE_PENALTIES` (private aliases) |
| Duplicate | `_make_locus_em_data` helper duplicated from `test_estimator.py` |
| Staleness | None — `_score_gdna_candidate` exists as alias at pipeline.py:80 |
| Quality | Good tests, but fragile coupling to private API |

**Recommendation: Trim.**
1. Extract `_make_locus_em_data` to `conftest.py`.
2. Import `score_gdna_standalone` from `hulkrna.scoring` instead of the private alias `_score_gdna_candidate` from `hulkrna.pipeline`.
3. Similarly, import `GDNA_SPLICE_PENALTIES` from `hulkrna.scoring` instead of `_GDNA_SPLICE_PENALTIES`.

---

### 19. `test_estimator.py` (1,006 lines, ~50 tests) ⚠️

| Attribute | Detail |
|-----------|--------|
| Classes | `TestIsAntisense`, `TestAssignUnique`, `TestScanData`, `TestLocusEM`, `TestLocusAssignment`, `TestPosteriorMean`, `TestSimultaneousResolution`, `TestMultimapperEM`, `TestCountsOutput`, `TestDetailOutput`, `TestGDNASummaryOutput`, `TestGDNAInLocusEM` |
| Duplicate | `_make_locus_em_data` — same helper as in `test_gdna.py` |
| Staleness | None |
| Quality | Thorough — convergence, determinism, posterior correctness, output formats |
| Overlap | `TestGDNAInLocusEM` class overlaps with `test_gdna.py:TestLocusBehavior`; both create loci with gDNA and verify EM behavior |

**Recommendation: Keep with cleanup.**
1. Extract `_make_locus_em_data` to `conftest.py`.
2. Consider merging `TestGDNAInLocusEM` into `test_gdna.py` or vice versa to avoid duplicate gDNA-in-EM coverage.

---

### 20. `test_resolution.py` (966 lines, ~45 tests)

| Attribute | Detail |
|-----------|--------|
| Classes | `TestMergeSetsWithCriteria` (6), `TestFragmentLength` (5), `TestResolvedFragment` (5), `TestDetectIntrachromosomalChimera` (5), `TestResolveFragment` (6), `TestComputeOverlapProfile` (7), `TestFilterByOverlap` (6), `TestOverlapFiltering` (3), `TestFragLengthDiscrimination` (5) |
| Staleness | None |
| Quality | Excellent — covers set merging, chimera detection, overlap profiles, fragment-length discrimination |

**Recommendation: Keep as-is.**

---

### 21. `test_scenarios.py` (1,361 lines, ~219 parametrized tests) ⚠️⚠️

| Attribute | Detail |
|-----------|--------|
| Classes | 10 scenario classes |
| External deps | minimap2, samtools (every test) |
| Cost per test | Full simulation → alignment → pipeline (several seconds each) |
| Total runtime | Estimated **15–40 minutes** for all 219 tests |

**Scenario inventory:**

| Class | Geometry | Parametrized cases |
|-------|----------|--------------------|
| `TestSingleExon` | 1 single-exon gene + ctrl | baseline + 4 sweeps = **21** |
| `TestSplicedGene` | 1 spliced gene + ctrl | baseline + 4 sweeps = **21** |
| `TestNonOverlappingGenes` | 2 non-overlapping genes + ctrl | baseline + 3 sweeps = **15** |
| `TestTwoIsoforms` | 2 isoforms of 1 gene + ctrl | 5 sweeps = **25** |
| `TestOverlappingAntisense` | 2 overlapping antisense genes + ctrl | 5 sweeps = **25** |
| `TestContainedAntisense` | Antisense gene in another's intron + ctrl | 5 sweeps = **25** |
| `TestParalogMultimapping` | Identical paralogs + helper + ctrl | 2 fixed + 4 sweeps = **19** |
| `TestDistinguishableParalogs` | Shared exon paralogs + ctrl | 4 sweeps + stress = **22** |
| `TestTwoExonWithControl` | 2-exon gene + ctrl | 8 method groups = **33** |
| `TestWideIntronInsertPenalty` | Wide-intron gene + ctrl | 4 methods = **13** |

**Issues:**
1. **Excessive parameter grid** — Each scenario sweeps gDNA ∈ {0,5,20,50,100}, nRNA ∈ {0,10,30,50,70}, strand ∈ {0.65,0.8,0.9,0.95,1.0}. Most sweeps test the same assertion pattern (alignment OK, accountability OK, negative ctrl OK). The middle grid points add runtime but rarely catch bugs that edge points miss.
2. **Heavy duplication of assertion patterns** — Nearly every test body follows the exact same pattern: `_build_and_run` → `_assert_alignment` → `_assert_accountability` → `_assert_negative_control` → optional extra check. The sweeps differ only in parameter values.
3. **Scenario class duplication** — `TestSingleExon` and `TestSplicedGene` are structurally identical (same sweeps, same assertions), differing only in gene geometry. `TestOverlappingAntisense` and `TestContainedAntisense` are extremely similar.
4. **`TestTwoExonWithControl` is a superset** — Tests all 3 axes plus all 2-way and 3-way interactions (8 methods, 33 cases). This alone provides strong parameter-space coverage, making the simpler scenarios' sweeps redundant for regression detection.

**Recommendation: Aggressive trim.**
1. **Reduce sweep grid sizes** — Use {0, 50, 100} for gDNA, {0, 30, 70} for nRNA, {0.65, 0.9, 1.0} for strand. This cuts each sweep from 5 to 3 levels, reducing total cases by ~40%.
2. **Drop stress tests from simpler scenarios** — `TestSingleExon`, `TestSplicedGene`, `TestNonOverlappingGenes` don't need 3-parameter stress combos; `TestTwoExonWithControl` already covers interactions.
3. **Merge `TestOverlappingAntisense` + `TestContainedAntisense`** into one class parametrized on geometry.
4. **Merge `TestSingleExon` + `TestSplicedGene`** into one class parametrized on `spliced=True/False`.
5. **Target: ~80 parametrized tests** (from 219), with the same scenario diversity but fewer redundant parameter combos.

---

## Source Modules Without Dedicated Tests

| Module | Tested Indirectly By |
|--------|---------------------|
| `locus.py` | `test_gdna.py`, `test_estimator.py` (via `LocusEMInput`/`Locus` creation) |
| `scan.py` (`EmDataBuilder`) | `test_pipeline_routing.py` (via `_scan_and_build_em_data`) |
| `scoring.py` | `test_gdna.py` (via `_score_gdna_candidate` alias) |
| `transcript.py` | `test_gtf.py`, `test_sim.py`, `test_oracle_bam.py` (as data carrier) |

These modules have reasonable indirect coverage. Adding dedicated test files is lower priority than trimming the scenario suite.

---

## Redundancy Map

```
test_gdna.py  ←─ _make_locus_em_data ─→  test_estimator.py   (DUPLICATE)
test_gdna.py:TestLocusBehavior   ↔   test_estimator.py:TestGDNAInLocusEM   (OVERLAP)
test_diagnostic.py   ↔   test_scenarios.py   (same pattern: Scenario → sim → pipeline → assert)
test_pipeline_routing.py mocks   ↔   test_buffer.py mocks   (similar mock classes)
```

---

## Staleness Check

**No stale imports found.** Every `from hulkrna.X import Y` was verified against the source module. Two imports reference private aliases (`_score_gdna_candidate`, `_GDNA_SPLICE_PENALTIES` in `test_gdna.py`) that exist but are fragile — they should be redirected to the public `hulkrna.scoring` API.

---

## Action Items (Priority Order)

| # | Action | Files | Impact |
|---|--------|-------|--------|
| 1 | Trim `test_scenarios.py` parameter grids | `test_scenarios.py` | −139 tests, −10–25 min runtime |
| 2 | Merge `test_diagnostic.py` → `test_scenarios.py` | `test_diagnostic.py`, `test_scenarios.py` | −1 file, cleaner organization |
| 3 | Extract `_make_locus_em_data` to `conftest.py` | `conftest.py`, `test_gdna.py`, `test_estimator.py` | DRY, single source of truth |
| 4 | Fix private imports in `test_gdna.py` | `test_gdna.py` | Future-proofing |
| 5 | Deduplicate gDNA-in-EM tests | `test_gdna.py`, `test_estimator.py` | −1 overlapping class |
| 6 | Share mock classes via `conftest.py` | `conftest.py`, `test_pipeline_routing.py` | Minor DRY improvement |
| 7 | Add `@pytest.mark.integration` to external-tool tests | Multiple files | Enables `pytest -m "not integration"` for fast dev loop |
