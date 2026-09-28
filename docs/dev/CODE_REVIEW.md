# Code review — the whole tree, 2026-09-27

A read-only review of `src/rigel/` (Python and C++), the build files, `src/rigel/sim/` and `scripts/`, done by nine
parallel reviewers (one per area) and consolidated here. The owner's brief: find bugs, dead code, stale and legacy
code, rotten comments and logical issues, and make the code simpler, clearer and more maintainable. Nothing was
changed. File:line pointers are as of `1ccdaa4c` and drift.

**Verification.** Every finding was checked by its reviewer against the code and a grep over `src/`, `tests/` and
`scripts/`. Findings marked ✔ were re-checked by the consolidating session. None was benchmarked: where a fix
changes behaviour it is marked **moves numbers**, and per `CLAUDE.md` it needs its own A/B before `src/`.

**Tooling.** `ruff check` is clean. `vulture` (run from a scratch directory, not installed) flagged 43 candidates;
most were false (fields the native side fills, or the report reads by name) and the true ones are below.

---

## 1. Bugs

### 1a. Results are wrong today (moves numbers; each is its own A/B)

1. ✔ **The realized gDNA length law is fed RNA COUNTS where it expects a pmf** (`calibration/fl.py:537`, called at
   `:781`). `g_C` is normalised at `:520`, `rna` is not: `rna_counts` comes from `detilt_pool`, which preserves the
   spliced-fragment total, and `contained_opportunity` is linear in its input. So the RNA rate `rho_r` is divided by
   N_spliced, `R_b → 0`, and every exon-flanking boundary reads as pure gDNA. The reviewer measured the intron|exon
   share at 0.026 / 0.963 / 1.000 for the same fixture at scale 1 / ×1e3 / ×1e6. It reaches both the realized law and,
   through `_couple_estimands`, the uniform gDNA law. The tests pass a normalised pmf, so they cannot see it.
   **Fix:** normalise at the top of `_realized_gdna_counts`. Ladder A/B required.
2. ✔ **The influence-weighted strand overdispersion is computed and then discarded** (`calibration/calibrate.py:417`
   then `:426-433`). Production overwrites `gdna_od` with `reconcile_overdispersions(raw_overdispersion, …)`, the
   ρ = 0 pair-count moment, paired with an `information` evaluated at the fitted ρ. The RNA side passes null
   information, which `gdna_strand.py:211-213` says reconcile must not be fed. **Owner question:** was the raw moment
   intended? If not, pass the fitted value. Stale docstrings beside it: `gdna_strand.py:46-48, :111, :148, :352, :524`.
3. ✔ **Nascent-span contributor count is overwritten** (`index.py:493`): `s_tx.nrna_n_contributors = …` inside a loop
   over merged spans, so a single-exon transcript covering two spans keeps the last count. It becomes
   `n_contributing_transcripts` in the nascent table. **Fix:** `+=`. Output changes only in that case.
4. **`lift_choices` keys records too coarsely** (`second_pass.py:535-569`). It keys on (ref, start, end,
   align_strand, sj_strand); the C++ canonical key adds the observed introns and the hypotheses
   (`native/calibration/accumulator.cpp:121-127`). Records tying on the short key but differing in introns are
   pooled, so a hypothesis index can land on the wrong record (a raise in `drain`, or a wrong replay in an oracle
   partition). **Fix:** key on the full record. Frequency on the panels unmeasured.
5. **Two policy benchmarks pool ss 0.70 into "unstranded"** (`scripts/design/policy_benchmark.py:176`,
   `_shared.py:50`, `quant_accuracy.py:666`): anything not `ss_0.99` counts as unstranded, so the test panel's 0.70
   rows land in the "must WIN" half against the never-pool rule. Ladder output unchanged.

### 1b. Broken features and crashes

6. ✔ **The HTML report's Vega-Lite marks were renamed** (`report/specs.py:62, :220, :270`). Rename commit
   `ffa05dd2` (line → boundary) rewrote Vega-Lite's own `"line"` keyword: `"type": "boundary"` fails to compile, so
   the gDNA-density genome chart never renders; `"boundary": {"strokeWidth": 2}` on two area marks is ignored, so
   their outlines are gone. No test compiles the specs. **Fix:** restore `"line"`, and add a spec-compile gate.
7. ✔ **A C++ exception in a scan worker aborts the Python process** (`native/bam_scanner.cpp:1336-1364`): the worker
   lambda has no try/catch (the reader thread has one), so e.g. `checked_u16_buffer_value`'s throw or `bad_alloc`
   calls `std::terminate`. **Fix:** catch, store an `exception_ptr`, abort both queues, rethrow after the join.
   Related: `thread_pool.h` has the same no-throw hazard (a throw on thread 0 returns before the workers finish).
8. ✔ **`--tsv` and `--emit-locus-stats` in a config YAML are silently ignored** (`cli.py:1146-1152, 1339-1348,
   943-960`): both are `store_true, default=False`, and the resolver treats any non-None CLI value as "set on the
   command line". `config.yaml`'s "rerun with" therefore drops them. **Fix:** `default=None` (e.g.
   `BooleanOptionalAction`).
9. ✔ **`FragmentBuffer`'s safety-net finalizer keeps the buffer alive** (`buffer.py:491`):
   `weakref.finalize(self, …, self)` passes the object to its own finalizer, which the Python docs warn prevents
   collection. Any buffer dropped without `cleanup()` holds its chunks (up to 2 GiB), spill files and writer thread
   until exit — in loops in `build_scan_cache.py`, `_oracle_arms.py`, `prior_vs_oracle.py`, `profiler.py` and tests.
   **Fix:** finalize on a small state object.
10. **A latent use-after-free in the nanobind casts** (`native/solve_kernel.cpp:752-756, 829-831`;
    `native/em_solver.cpp:2303-2320, 182-191, 1880, 1943-1944, 2000`). `nb::cast<A>(h).data()` converts into a
    temporary when dtype or order mismatch, and the pointer dangles (with the GIL released); output arguments without
    `.noconvert()` write into a discarded copy. Safe today only because every caller passes exact dtypes.
    `StreamingScorer` keeps borrowed arrays without `nb::keep_alive`. **Fix:** `nb::cast<…>(obj, false)`,
    `.noconvert()` on outputs, `keep_alive`.
11. **`--overhang-alpha 0` (documented as a hard gate) yields NaN** (`native/scoring.cpp:580-581, 870-871`):
    `oh * log(0)` is `0 × −inf` for zero-overhang candidates, in a `-ffast-math` unit. Mismatch is guarded, overhang
    is not. **Fix:** require alpha in (0, 1] or carry the gate as a bool. Default path unchanged.
12. **The scan cache can accept a stale or torn cache** (`scan_cache.py`):
    - `deposit_digest` (`:247-274`) hashes only ndarray outputs of a fixture whose fragments offer no hypotheses, so
      a change to gap arbitration, deferral or the qc classification leaves the key unchanged;
    - `_payload_from_parts` (`:490-516`) bypasses `from_dict` and both conservation checks that
      `scan_payload.py:773-777` justifies by "a cache can be truncated";
    - `write_scan_cache` is not atomic, so an interrupted `--force` rebuild leaves new files beside the old manifest;
    - `_WORK_ONLY_FIELDS` (`:290`) omits `spill_dir`, `buffer_size_bytes` and three other non-tally fields, so a
      cache is refused for a different spill directory.
13. **GTF parsing edge cases** (`gtf.py:42`, `transcript.py:154, 160`): the attribute regex needs a trailing `;`, so
    a last attribute without one is dropped; a missing `transcript_id` raises a bare `KeyError`, which `warn-skip`
    does not catch; an exon on another chromosome or strand than its transcript is merged in silently.
14. **Simulator: a capture entry enabled by default comes back disabled** (`sim/whole_genome.py:131`): the probe
    check fires only for an explicit `enabled: true`, so `{label: on}` with no probes yields a `_capture_on`
    condition with no capture.
15. **`panel.py --index` is not forwarded to the simulator** (`scripts/sim/panel.py:203-208`), so the panel is
    simulated on the config's index but cached and scored on the other.

### 1c. Test fixtures that do not test what they say

16. ✔ **`build_oracle` ignores `gdna_fraction` unless `n_rna_fragments` is given** (`sim/scenario.py:512-540`): the
    fixed-total branch takes gDNA from `gdna_config.abundance` (0 without a config). Six fixtures pass
    `gdna_fraction` alone and simulate NO gDNA: `test_second_pass_pipeline.py:75`, `test_d7_transcript_eff_lengths.py:50`,
    `test_summary_report.py:164`, `test_scan_order_independence.py:68, 168, 379` (whose docstring at `:344` says it
    puts genomic fragments in introns). **Fix:** raise on the combination, then repair the fixtures — their
    expectations will move.

### 1d. Latent (not reachable on today's inputs)

17. `calibration/total_abundance.py:133-134`: `IndexError` when the LAST reference owns no regions (reproduced).
18. `calibration/splice_graph.py:219, 430-443`: `_ref_slices` assumes exon rows grouped by reference; unsorted input
    silently drops exons (reproduced). Production input is sorted by contract only.
19. `calibration/messages/transfer.py:183`: the RNA-coordinate fallback fires only at `den <= 0`; with single-strand
    exons all at zero count, `rho_rna = 0` and no RNA lane is built — the failure the gDNA fallback prevents.
20. `calibration/calibrate.py:634-645`: a failed later refit overwrites the last landscape with `None` before the
    `break`, so the reference density, every efficiency and the debug record are lost.
21. `calibration/calibrate.py:376-378`: injecting `rna_sense_frac` without `n_rna_obs` kills the strand channel.
22. `native/fast_exp.h:77-99`: an all-`−inf` row gives `x = NaN`, and `static_cast<int64_t>(NaN)` is undefined.
23. `native/em_solver.cpp:2806-2816`: an early return yields a 2-tuple where callers unpack 5.
24. `native/resolve_context.h:300-324`: `finalize()` leaves `t_offsets_` without its leading 0.
25. `native/calibration/accumulator.cpp:699-702` (and the spec at `:964-965`): a first intron starting exactly at the
    fragment start shifts every junction end by one. Fix the spec first.
26. `sim/sampling.py:33`: `.astype(int)` truncates, so the simulated fragment-length mean is `frag_mean − 0.5` while
    `fl_pmf` (the ruler's truth) evaluates the density at integers. An empty `[frag_min, frag_max]` window loops
    forever (`sim/wgs_engine.py:314-315`).
27. `native/bam_scanner.cpp:87-94`: an unknown `--sj-strand-tag` spec silently becomes `XS_TS`.
28. `pipeline.py`: the auto-detected strand tag is never logged or recorded in `summary.json` / `config.yaml`
    (`:918-923, 1067-1071`); `second_pass_seed=None` draws OS entropy (`config.py:392`).
29. `cli.py:1197`: `--em-iterations 0` is documented as "unambiguous-only", but still runs one SQUAREM iteration and
    assigns every unit.
30. `truth.py:142`, `whole_genome.py:952`: `genomic_span` is 0 for a transcript starting at position 0.

---

## 2. Dead code (behaviour-preserving to delete unless noted)

The large ones, each an end-to-end chain:

1. ✔ **`gdna_locus_counts`** (`em_solver.cpp:1836-1850`; `estimator.py:221, 383`): filled, never read. It is the
   only consumer of the per-unit `locus_t_indices` / `locus_count_cols` chain (scorer → scan → `ScoredFragments` →
   scatters → partition tuple slots 7-8 → `LocusSubProblem`). Delete the chain; the tuple shrinks to 7.
2. **`genomic_midpoint`** (`scoring.cpp:416, 450-452, 538-546, 686-690, 939-952, 1164`; `scan.py:167, 223`;
   `scored_fragments.py:66-70`): 8 bytes per EM unit kept alive through the EM; its named consumer does not exist.
3. **`gdna_splice_penalties`** (`config.py:122` → `scoring.py` → `pipeline.py:553` → `scan.py:115-119` → the C++
   argument): the only value is `{UNSPLICED: 1.0}`, so `log(1) = 0` is added to every gDNA likelihood.
4. **`exon_bp_pos/neg`, `tx_bp_pos/neg`** (`constants.h:202-205`, `resolve_context.h:1117-1153`,
   `resolve.cpp:45-48`): a per-transcript strand lookup on every interval hit in the hot path; nothing reads them
   and tests assert they are absent from the buffer.
5. **`chimera_gap`** (`constants.h:446-463`, `resolve_context.h`): an O(components² × blocks²) loop computing a
   value that is never exported.
6. **The EM's float64 payload path** (`em_solver.cpp:149-197, 2306-2312`, `locus_partition.py:29-48`): the scorer
   emits float32 only; the float64 dispatch and `scatter_*_f64` exist for two test harnesses.
7. **Calibration geometry fields** (`calibration/region_geometry.py:146-163, 259-328`; `substrate.py:73-78,
   186-189`): `inv_sj_lo/hi` unread; `inv_abundance`, `eff_sj*` test-only; the substrate's `inv_length_sum` /
   `inv_opportunity_sum` channels feed only these. About 60 lines of comments go with them.
8. **`splice_graph.build_transcript_path`** and its types (`calibration/splice_graph.py:1401-1597`, ~200 lines):
   used only by `tests/calibration/test_transcript_path.py`.
9. **`AbundanceLandscape.w_slot`** (`calibration/abundance_landscape.py:296-306`): a dense (n_regions × 260) float64
   matrix built on every capture-ON run and read only by tests — exactly what `_render` tiles to avoid. Also
   `rho_0` (unread), `span_R`, `anchor_*` (test-only).
10. **`CalibrationResult.rna_region_eff_len` / `rna_boundary_eff_len`** (`result.py:161-178`, `calibrate.py:849-851`)
    and the second `_project_eff` call.
11. **`_IntronFactory`'s per-grid cache** (`calibrate.py:172-222, 438-481`): `FactoryRows.fg` / `.shape` unread, and
    `kernel()` returns only grid-independent inputs, so the cache rebuilds identical objects. Bit-identical.
12. **`annotate.frag_id_to_row`** (`annotate.py:232, 273, 316, 349-365`): a Python dict with one entry per annotated
    fragment (many GB on a deep library with `--annotated-bam`); only a test reads it.
13. **The Python interval index on every load** (`index.py:1328-1348, 1577-1614`): `cr` (cgranges) and `_iv_type`
    are built and never freed; the only reader, `query()`, is test-only and raises in production. The
    `_cgranges_impl` extension exists only for it.
14. **`buffer.py`'s test-only surface**: `BufferedFragment`, `append`/`finalize`, `__iter__`, `summary` and its
    counters; the `sj_strand` and `merge_criteria` columns are stored and spilled for every fragment but never read.
15. **The simulator's retired paths**: the `fragment_share` nascent mode (`wgs_config.py:91-113`,
    `whole_genome.py:239-245, 291-295, 739-782, 1009-1015`; live only in the stale fl-gap configs; the 11
    `test_reference*.yaml` carry inert blocks), `transcript_filter` (only "all" is legal), the orchestrator's unused
    parameters and always-None manifest fields, the `whole_genome` re-export shim, `pgzip`.
16. **`calibration/strand_summary.py`**: its identifiability test (a 99 % z-test and a 1e-3 floor) duplicates the
    Bayes factor calibration actually uses, so `pipeline.py:163-192` can warn differently from what calibration
    does, and names an estimator that no longer exists. No tests.

Smaller ones:

- Native: `rows.lgamma`, `rows.pairwise_sum` bindings; `solve_block`'s `written`; `Arena::fp_loc/fn_loc`;
  `ChainArrays::n`; `Params::od_g`; the gDNA `SolveLane`'s `own_rows/own_mask/flux_*`; `ResolverScratch::tmp_b/tmp_out`;
  `ParsedAlignment::mate_ref_*`; `parse_cigar`'s `sj_strand` parameter; `sj_base`; six `AF_*` constants duplicated
  from `annotate.py`; `n_chimeric` (exported, never incremented); `squarem_grouped_fallback_used` (always false),
  `bias_us` (times nothing); `_resolve_core`'s bool return; unused includes; `using FRAG_AMBIG_OPP_STRAND`.
- Three of the five sweep assertions cannot fire (`solve_kernel.cpp:680-685, 712-714`).
- Python: `index.t_to_ref_arr`; `transcripts.feather`'s `abundance`, `nrna_abundance`, `n_exons` and `sj.feather`'s
  `interval_type` (an index-format change); `GdnaDensityFit.total_counts`; `RegionWallMask.w_max`; `CubeRows.select/
  shifted/nbytes`; two `_EPS`; `substrate.strand_class`, `region_span_count`; `blocks.SweepCapture`'s six unread
  fields; `FLModels.gdna_contrast/gdna_realized`; `ChainView.n_grid/logodds_window`; `density_deconv.include_introns`;
  `result.config`; `scoring.LOG_HALF`; `scan.py`'s `n_chim` and `_native_ctx is None`; `HeldScores.n_hypotheses`;
  `native.Accumulator`; `frag_length_model.pmf` and the raw-count smoothing branch; eight report-model fields
  report.js never reads; the `pileup` extra in `pyproject.toml`; `track.capture_summary` (a third capture census, with
  unexplained constants and a swallowed exception).
- Used only by tests or scripts (move to tests, or keep deliberately): `InjectedCalibrationPriors`,
  `_project_regions_to_loci`, `contended_boundaries`, `psi_cube`, `gdna_arm`, `posterior_median_fg`, `compose`,
  `fl_mean`, `split_basins`, `lattice_points`, `load_manifest`, eight `types.py` / `splice.py` members.

---

## 3. Performance (behaviour-preserving)

1. `estimator.get_loci_df` (`estimator.py:835-887`) is O(loci × transcripts): measured 0.48 ms per locus at 450k
   transcripts, about 19 s at 40k loci, in the CLI writer. **Fix:** `np.bincount`.
2. `scan.py:177-208`: with `--annotated-bam`, a Python loop calls `annotations.add` once per deterministic fragment;
   the batch API already exists.
3. ✔ `em_solver.cpp:2094-2106`: every locus clears `local_map` up to its largest global transcript index, although
   the trailing loop restores every entry it set.
4. `resolve_context.h:176-212`: the scan path copies `gap_introns` / offsets per resolved hit (a malloc each) into a
   test-only binding type.
5. `sim/capture/sampler.py:269-273`: the partition memo holds one population, so the mRNA and gDNA spaces evict each
   other; the promised reuse across conditions almost never happens. `_transcripts_overlapping` is O(T) per probe.
6. `pipeline.py:354, 406, 462`: `build_region_partition_arrays` runs four times per run, `build_sj_arrays` twice.
7. `solve_kernel.cpp:527-540`: slots that are in-slots but not solvable build the whole cube and discard it; the
   thread pool is re-created on every call with fewer tasks than cores (`:103-117`).
8. `em_solver.cpp:405, 574-588`: a heap vector per equivalence class per E-step; the task list rebuilt every E-step.

---

## 4. Duplication and simplification

1. **The gates test a copy of the production wiring**: `solve_one_block` and `transfer_prepare`
   (`solve_kernel.cpp:500-608` vs `924-1000`) wire the chain and build the lanes twice, and the gates exercise the
   second. One shared builder.
2. **`scoring.cpp` scores twice**: candidate scoring, pruning and unit finalisation are written once for
   multimappers and once for single fragments (`:555-697` vs `:845-964`); the deterministic and EM antisense tests
   disagree at p_sense = 0.5 exactly (`:788-800` vs `:856-864`).
3. **SQUAREM's VBEM and MAP branches** (`em_solver.cpp:1161-1373`) differ only in the weight, the floor and the
   normalisation: one template, about 150 lines fewer.
4. **`BamAnnotationWriter::write` re-implements `parse_bam_record`** (`bam_scanner.cpp:2649-2687` vs `599-653`), and
   the multimapper detection twice.
5. **`TranscriptIndex.load`** (~250 lines) converts arrays to frozensets and dicts and back to CSR for the resolver,
   with the same lexsort-and-diff grouping twice. Build the CSR directly; serve exons from the cached exon CSR (the
   simulator must pass `retain_test_structures=True` today).
6. **One annotation walk**: `intervals.feather` exons are re-read and re-expanded in `splice_graph.py:1296-1314`,
   `capture_eff_length.py:95-131` and `_Exons`.
7. **One ramp primitive**: `gdna_density.contained_opportunity(pmf, L)` equals `effective_length.contained_eff_length`,
   and `sj_opportunity._ramp_sum` equals `gdna_opportunity.contained_opportunity(L, W)` — two different functions
   share the name `contained_opportunity` in one layer.
8. **`FragmentScorer`** (`scoring.py:78-227`): 13 frozen fields copied into the native scorer, `_native_ctx`
   attached with `object.__setattr__`. Return the native scorer directly; store alphas, not log penalties, and the
   `cli.py:869-887` special case goes.
9. **One validating constructor for `AccumulatorPayload`**, shared by the scan, the cache and the drain.
10. **`pipeline._apply_scan_stats`** hand-lists 33 keys with `.get(key, 0)`, so a renamed C++ key reads 0 silently.
    Iterate the dataclass fields and index strictly.
11. **The simulator's truth**: ten `ground_truth_*` methods and four file passes in `run_benchmark`; one pass would do.
    `tests/calibration/_oracle.py` is imported by five instruments through `sys.path` hacks; move it to
    `scripts/design/`. The runs directory is spelled five times despite `_shared.py`.
12. **The numpy bit-matching machinery** (`pairwise_sum`'s copy of numpy's order, the `lgamma` binding, named
    temporaries) no longer has a Python counterpart. Retiring it moves results by ulps: an owner decision.

---

## 5. Rotten comments and rule violations

- **Rename damage** from the automated renames: `line → boundary` in Vega-Lite (a bug, above), "cut → REGION_BOUND"
  (`config.py:111-115`, `scan_payload.py:453`), "junction → SpliceJunction" (`scan_payload.py:21`,
  `strand_model.py:145`, `index.py:536`), `num_boundaries` for exon features (`transcript.py:150`), "BED12 boundary",
  "splice sj" (three files). The renames were proven numerically identical, which cannot see prose or JSON keywords.
- **Source cites the docs** in about ten places (`calibrate.py:875`, `abundance_landscape.py:210, 223`,
  `priors.py:328`, `simplex_logodds.py:401`, `capture_efficiency.py:15`, `psi_kernel.h:80, 231, 258`,
  `transfer_rows.h:214`). `test_docs_boundary.py` polices only citations into `docs/dev`, so the rule has no gate.
- **References to deleted code**: ports of `bam.py`, `fragment.py`, `resolution.py`; `transfer_rows.EPS`,
  `faces.py`, `strand_row_logodds`, `_mixture_strand_loglik`, `density_deconv._log_negbinom`,
  `region_init.strand_evidence`, `sweep._gdna_arm`, `logprior`, `parse_bam_file`, `mappable_effective_length`,
  "count-clue density", "invariant I2", `ec.scratch`, `rigel.categories`, `merge_sets_with_criteria`,
  `synthetic_genome.inject_splice_sites`, `suite.main`, `--no-fastq`, `gdna_total`, `shadow_init`.
- **Message-cache leftovers** (deleted in `25511abe`): `sweep.py:11-15`, `solve_kernel.cpp:400-401`,
  `transfer.py:153-155`.
- **Wrong claims**: "integer addition is associative, so identical at any worker count" above float64 merges
  (`accumulator.cpp:891`, `accumulator.h:323-325`); "nothing in src/ reads them" (`native.py:60-61`,
  `solve_kernel.cpp:20, 1164`); the assignment weights "mirror the E-step" under VBEM (`em_solver.cpp:1500-1504`);
  "three assignment modes" (two since `map` was deleted); `fast_exp.h:54-55` on underflow; `HeldScores.n_undecided`
  "all scored zero"; the ZS tag listing four labels (`annotate.py:73-76`); `strand_model.py:184-186` on caches that
  do not exist; `strand_model.py:328-358`'s promised closed form; schema "v2" in `report/substrate.py`.
- **History narrated in source** (belongs in git): accumulator headers, "Phase 0/2/3", "was named … until",
  "replaces Python …", stale measured numbers (`calibrate.py:253-254, 624-625`, `transfer_rows.h:214`).
- **Unexplained constants**: queue sizes and reserves in the scanner; `1e-9`/`1e-12` in `transfer_kernel.h:287`,
  `transfer_rows.h:255, 263`; `track.capture_summary`'s five; the hand-typed `z = 1.959964`; `_DEFAULT_L = 10.0`
  duplicating `CalibrationConfig.sweep_logodds_window`; hand-typed `RegionType` codes in four places.
- **Build**: the SIMD block in `CMakeLists.txt` sets `RIGEL_ARCH`, which nothing reads, and claims runtime dispatch
  that does not exist; `-ffast-math` implies `-ffinite-math-only` while `scoring.cpp` uses `−inf` sentinels.

---

## 6. Decisions for the owner

1. The strand-overdispersion estimator (bug 2): fitted or raw moment into `reconcile_overdispersions`. **Decided 2026-09-28: one shared value** (`ISSUES: strand-overdispersion-one-shared-value`).
2. `cli.py:522-534`: `gdna_fraction` in `summary.json` counts SPLICED intergenic fragments as gDNA, against Axiom 0.
3. The per-transcript `rna_prior_weight` lane and `warm_start="prior"`: production plumbing whose only producer is
   `quant_accuracy.py`'s oracle-allocation arm. Keep as an instrument hook, or delete.
4. `region_span_count` (`accumulator.cpp:661-676`): deposited, exported, read by nothing in `src/`; it came in with
   the START/END/SPAN taxonomy.
5. The `fragment_share` nascent mode: retiring it needs the fl-gap panels re-simulated or retired.
6. Dropping dead index columns bumps `INDEX_FORMAT_VERSION`.
7. Retiring the numpy bit-matching machinery (ulp-level moves; a golden and replay re-record).

## 7. A proposed fix order

One case at a time, each verified before the next (`CLAUDE.md`'s working rules):

1. **Behaviour-preserving, no number can move**: the Vega-Lite marks, worker exception handling, the finalizer, the
   YAML booleans, the dead chains in §2, the comments in §5. Gate: the suite, plus `rename_identity.py --check` on
   a frozen reference to prove each is a numeric no-op. **DONE 2026-09-27** in thirteen packages, each proven
   bit-identical on three frozen references (`g05 ss.50 OFF`, `g05 ss.99 ON`, the real LBX0190); what they left is §8.
2. **Performance**, each proven a numeric no-op the same way.
3. **Bugs that move numbers**, each with a falsification test first and its own ladder A/B: the realized law (bug 1),
   the contributor count (3), `lift_choices` (4), then the owner's call on the overdispersion (2).
4. **The fixtures** (16), whose expectations will move.
5. **The larger simplifications** in §4, one per change.

## 8. Step 1 follow-ups — DONE 2026-09-28 (step 1b)

All but item 17 landed in eleven packages, each proven bit-identical on the three frozen references, the
transcript table unchanged. The only arrays gone are the dead ones deleted:
- step 1's two `CalibrationResult` arrays;
- step 1b's `region_contained_inv_opportunity_sum` (a payload-schema change, so every scan and oracle cache was
  rebuilt).

Outcomes that are not self-evident from git:
- **Item 5 (test-only code).** `split_basins`, `lattice_points` and `_project_regions_to_loci` turned out to be
  live, so they stay. `InjectedCalibrationPriors` stays by the owner's ruling: it is the seam the toy harness
  injects priors through.
- **Item 13.** The scanner's per-reference intron filter is not redundant (a CIGAR with only an `N` reaches the
  deposit); the comment now says so.
- **Items 20 and 21.** `-fno-finite-math-only` and `_solve_impl`'s optimisation flags both proved bit-identical and
  stayed.
- **Item 25 (owner, 2026-09-28).** The `t_index` presence check is gone. The row-alignment check stays and is now
  tested, because a reordered transcript table would silently mis-map every per-transcript array. `rigel index`
  refuses a GTF with no transcripts. A build deletes the old manifest before its first write, so an interrupted
  rebuild is refused.
- **Item 26.** Confirmed by the owner.
- **Item 17** waits for the strand-overdispersion work.

## 9. Step 1c — what step 1b left, with the owner's rulings (2026-09-28) — 9a–9d DONE 2026-09-28

Each is its own change, verified like steps 1 and 1b unless it is marked as moving numbers.

### 9a. Instruments (first: the strand-overdispersion A/Bs need them) — DONE

Outcome: all five fixed, each with a gate verified failing and then fired by a deliberate break.
- `solvability_audit.py` reads the truth from the simulator's ledger (`_shared.pool_ledger`) and raises when it is
  missing. `--oracle-cache` defaults to `<suite>/oracle_cache`.
- `quant_accuracy.py` labels a rung from its own `gdna_frac_true`, not from a map of known rungs.
- `_shared.strandedness` is the one rule, and it raises on a name with no `_ss_` token. Any strand specificity other
  than 0.50 or 0.99 is reported APART, never pooled.
- `ruler_vs_truth.py` and `policy_benchmark.py` put `src` on `sys.path` only when no rigel is importable.
- Left: `scripts/sim/build_test_reference.py` has the same unconditional insertion. `solvability_audit.py` has no
  `--self-test`; its gate is a suite test.


1. `solvability_audit.py` is broken on the current tree: it calls the deleted `_oracle_arms.truth_f_gdna`
   (recorded in ISSUES since 2026-09-24). Repair it.
2. `quant_accuracy.py --markdown` crashes on panels with a g25 rung (`_GDNA_LEVEL` has no g25 entry).
3. `calibration_vs_oracle.py` ignores `--jobs` when given `--json`.
4. `ruler_vs_truth.py` puts `<repo>/src` first on `sys.path` unconditionally, so it breaks in a worktree
   (`policy_benchmark.py` guards the same insertion).
5. `policy_benchmark.py`, `_shared.py` and `quant_accuracy.py` fold ss 0.70 into "unstranded" (§1a bug 5).

### 9b. The report — DONE

Outcome:
- `summary.json`'s `calibration.capture` block, its census and the density chart are deleted. The report's capture
  tile reads `gdna_reference_density` and `gdna_reference_members`.
- MANUAL lists exactly the nine keys written, and the CHANGELOG records the removal.
- Open (owner): a fold (on-target over off-target gDNA density) would need `split_basins`' depleted mode as a new
  `CalibrationResult` field.


6. Delete `track.capture_summary`'s separate KDE census (owner). `summary.json`'s capture block and the report read
   calibration's own answer (`located_enriched_mode`, `split_basins`). That removes its five unexplained constants.
   The report's capture numbers change; quantification does not.
7. `report.js` reads `c.n_nodes`, which the rename made `n_regions`, so every report with a capture block shows
   "NaN nodes". MANUAL's `summary.json` capture section is stale: it documents `n_nodes`, `separation_nats` and
   `enrichment_factor`, and omits the keys actually written. Fixed together with item 6.

### 9c. The index — DONE

Outcome:
- `build()` removes an earlier blacklist.
- `load()` applies the blacklist only when the manifest records `sources.alignable_zarr`. A recorded blacklist
  whose file is missing loads with detection off (`sj_blacklist_loaded: false`), with no refusal.
- Three tests that placed a blacklist file by hand now rebuild through `tests/_index_builder.py`'s
  `rebuild_with_splice_blacklist`.
- The cluster's production index must record its store before its next quant (GDNA_SPLICE_ARTIFACTS_PLAN.md).


8. A stale `splice_blacklist.feather` survives a completed rebuild without the alignable store, and `load()`
   applies it whatever the manifest records (owner: an index rebuild must address the blacklist). The blacklist
   should be used only when the manifest records its source. Test: build with a blacklist, rebuild without one,
   load, and expect no blacklist.
9. Hardening from the index-integrity package's review:
   - a test that an empty-GTF rebuild into an existing index leaves that index loadable;
   - `_mini_sources` reused by the two inlined fixtures in `test_index_integrity.py`.

### 9d. The EM and the length laws — item 10 DONE (already at 3bb35bdd); item 11 analysed

Outcome:
- Item 10 had already landed at 3bb35bdd. That commit removed three `locus_stats` columns, not the two its message
  names.
- Item 11 is derived in `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/06_fl.md`. The guards that never bind
  go first, bit-identical. The realized-law fix, the empty pools and a single refresh follow, each after its own A/B
  and the owner's rulings.


10. Delete `squarem_extrapolation_clamp_count` (owner). Since the SQUAREM backtracking fix it counts only
    components the EM's own step had already floored.
11. A focused analysis of `fl.py`'s constants (owner), taken with the realized-law bug (§1a bug 1), which lives in
    the same function (`_realized_gdna_counts`):
    - four zero-guards never bind on a reachable input (`max(·, 1e-30)` three times, `max(μ − 1, 1e-9)` twice). They
      can go, or become explicit branches, with identical output.
    - the 0.25 bp refresh test in the boundary-stratum loop is a genuine tunable with no derivation, and it moves
      numbers.
    - `transfer_rows.h`'s `1e-12` never binds on a reachable input.

### 9e. Watch

12. `test_scan_order_independence::test_THE_FIXTURE_REALLY_DOES_REORDER_THE_BUFFER` failed once on thread timing, a
    dependence it documents, and passed on three reruns.

### 9f. Real-data findings from the strand-overdispersion measurement, the same under every arm (for ISSUES)

Outcome (investigated 2026-09-28; pages in `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/pages/`):
- 13 is `ISSUES: a-pure-gdna-library-reads-as-nascent-rna`.
- 14 is the gDNA landscape collapsing at low depth: the refits own 69–86 % of the drift.
- 15 is not a defect. It is the `test_blank` control's unannotated transcription, locked to gDNA by structure.
  Beside it, the intergenic background's dispersion is fitted on a pool that transcription contaminates.
- Whether 14 and that dispersion fit become ISSUES entries is the owner's call.


13. The pure-DNA VCaP exome half is read as about 93 % RNA (the deferred unstranded × capture-ON stratum).
14. LBX0588's gDNA share per deposited fragment moves with depth: 0.11 → 0.45 → 0.83 at 1 % / 10 % / full.
15. `full_lowg` stranded × capture-OFF over-calls gDNA 2.3×.
