## Scope read

- **Read in full:** every file in `src/rigel/calibration/`: `_layers`, `__init__`, `errors`, `blocks`, `sweep`, `calibrate`, `landscape`, `simplex_logodds`, `region_geometry`, `region_init`, `region_chain`, `capture_efficiency`, `capture_eff_length`, `effective_length`, `gdna_density`, `gdna_opportunity`, `sj_opportunity`, `fl` (1–846), `priors`, `result`, `region_arrays`, `substrate`, `signature`, `strand_balance`, `derive`, `density_model`, `track`, and `messages/{__init__,silent,transfer}`.
- **`splice_graph.py`:** read 1–65, 115–416, 845–880 and 886–1163; the helpers and validator in between were skimmed.
- **Native, read in full:** `solve_kernel.cpp` (1–1224), `transfer_kernel.h` (1–814), `transfer_rows.h` (1–404), `psi_kernel.h` (1–348).
- **Cross-checks outside the slice:** `test_layering.py`, `test_sweep_backbone.py` 455–495, `test_region_init.py` 285–310, `_transfer_harness.py` 595–610, `test_capture_efficiency.py` 295–310, `ruler_vs_truth.py` 240–265, `pipeline.py` 379–406, `scan_cache.py` 99–118, `config.py` 785–797. Every test name the slice cites exists.
- **Clean on four rules:**
  - Layering: no upward import (every sibling import is relative, so the gate sees it). `capture_eff_length`'s layer-7 import is only inside `TYPE_CHECKING`.
  - The slice cites no doc (no `docs/` and no `§` anywhere in it).
  - No Greek identifiers.
  - Nothing in the slice repairs alignments.

## Findings

**1. Two "own composition evidence" predicates under one name**
- **Where:** `native/solve_kernel.cpp:420-421, 607`
  ```
  as_bool(A.has_own)[i] = tau > 0.0;
  as_bool(A.own_ev)[i] = tau > OWN_EPS;
  out.has_comp[start + i] = as_bool(A.own_ev)[i] || (!fp[i] && !fn[i]) || as_bool(A.held)[i];
  ```
  and `region_init.py:38,58` ("THE ONE DEFINITION…", `return ... > _EPS`).
- **Category:** duplicate / naming. **Severity:** medium.
- **What is wrong:** τ > 0 gates the transfer policy's own claims, the level rule's strand term and the splice-site width (`transfer_kernel.h:223,295,330`). τ > 1e-9 gates `has_composition`, which decides landscape training. The Python "one definition" is the second form, but the gate harness builds its `has_own_composition` as τ > 0 (`_transfer_harness.py:609`). A slot with 0 < τ ≤ 1e-9 sends a claim but is not counted as having a composition. The "divide-by-zero guard" rationale is vestigial: nothing divides by τ.
- **Evidence:**
  - `rg -n "OWN_EPS|has_own|own_ev" src/rigel/native` shows `OWN_EPS` used only for the clamps and this threshold.
  - At κ = 0.95, τ ≈ 4.3·n·f², so τ < 1e-9 iff f < 1.5e-5/√n. That is below the L = 10 lattice floor (4.5e-5), so the window opens only inside a refit's widened bracket or on a weakly stranded live library.
  - `test_the_evidence_predicate_agrees_with_the_kernel_liveness_test` uses τ ∈ {0, 1e-4, 1, 850}, never the window.
- **Proposed change:** one predicate (τ > 0), with `OWN_EPS` kept and renamed as the clamp only; extend the gate into (0, 1e-9].
- **Moves numbers:** yes, only for slots in the window (count unmeasured).
- **Confidence:** code confirmed; how often the window is reached is suspected.

**2. Contiguous-boundary reaches are computed, stored, gated and hashed, but read by nothing**
- **Where:** `splice_graph.py:317` `c_reach = _contiguous_reaches(ex, rows, pos)`, against `region_geometry.py:152-155`:
  ```
  BOUNDARY, RNA     ``UNBOUNDED_REACH`` by ruling ⇒ ``mu_r − 1``. ...
                    tenth. A per-boundary taper would be an A/B, never run
  ```
- **Category:** dead (speculative). **Severity:** medium.
- **What is wrong:** the only reach consumer, `build_sj_geometry_arrays` (`splice_graph.py:1100`), reads sj rows only. `scan_cache.reach_digest` ("the part of the index calibration reads") also hashes the contiguous rows, so changing values nobody reads invalidates every scan cache.
- **Evidence:** `rg "reach_(lo|hi)_(pos|neg)" src scripts` hits only `splice_graph.py` and `scan_cache.py`. Nothing in `src/rigel/native` reads them. `test_splice_graph.py:338-518` gates their values.
- **Proposed change:** either give the taper A/B an ISSUES entry that owns this code, or zero the contiguous reaches, digest sj rows only, and trim the module docstring.
- **Moves numbers:** no (cache keys change). **Confidence:** confirmed.

**3. `SweepCapture` has five fields nothing reads; one is named `mature`, one costs a ψ solve per slot**
- **Where:** `blocks.py:57-59` ("Every field has a reader … a field nothing reads is dead surface"); `sweep.py:357-365`:
  ```
  var_g=out["var_gdna"].copy(),
  spliced=view.spliced_slot,
  mature=np.asarray(view.sj_count, np.float64).sum(axis=1),
  eff_gdna=np.asarray(view.eff_gdna, np.float64),
  ```
  and `solve_kernel.cpp:597-603`, the strand-only solve behind `fg_strand`.
- **Category:** dead + naming. **Severity:** medium.
- **What is wrong:**
  - `fg_strand`, `var_g`, `spliced`, `mature` and `eff_gdna` have no reader in `tests/` or `scripts/`.
  - `fg_strand` adds a ψ solve per solved slot on every captured sweep. That includes every `ruler_vs_truth.py` calibration, whose arms read none of the capture.
  - `mature` names the sj-flux sum, which `test_the_word_MATURE_names_no_field` forbids on `RegionGeometry` for exactly this reason.
- **Evidence:** a field-by-field `rg "\.f\b|'f'"` sweep finds 0 readers. The getattr gate at `test_sweep_backbone.py:477` reads only `f_g`, `fg_loc`, `tau_lam`, `solvable`, `mass_global` and `count`.
- **Proposed change:** delete the five fields, `Outputs.fg_strand` and the strand-only solve; fix both docstrings.
- **Moves numbers:** no. **Confidence:** confirmed.

**4. `build_fl_models`' optional divisors are production-dead; scoring tests use gDNA laws production never builds**
- **Where:** `fl.py:753-758`:
  ```
  if sj_opportunity is not None:
      rna_counts = detilt_pool(rna_counts, sj_opportunity)
  realized_counts = None
  if gdna_opportunity is None:
      gdna_counts = gdna_contained_fl_mass(payload)
  ```
- **Category:** dead + test-gap. **Severity:** medium.
- **What is wrong:** the only production caller, `pipeline.library_fl_models` (`pipeline.py:400-406`), always passes all four inputs, and its own docstring warns that a different length model "would be measuring a tool that does not ship". The None branches (no de-tilt, the contained-pair fallback, no contrast, no realized law) are reached only by tests.
- **Evidence:** `build_fl_models(payload)` is called at `test_second_pass_scoring.py:234`, `test_second_pass_placements.py:104`, `test_second_pass_strand.py:135`, `test_scanner_accumulator_integration.py:313`, `test_scan_cache.py:205` and `test_d7_transcript_eff_lengths.py:91,112`.
- **Proposed change:** make the four inputs required, route those tests through `library_fl_models`, and delete `gdna_contained_fl_mass`.
- **Moves numbers:** production no; test expectations may move. **Confidence:** confirmed.

**5. ψ's signal mask solves spliced-only slots and throws the answer away**
- **Where:** `simplex_logodds.py:342-345`
  ```
  signal = (
      u_pos + u_neg + np.asarray(mass_unspl, np.float64) + np.asarray(mass_spliced, np.float64)
  ) > 0.0
  ```
  against `psi_kernel.h:316-317` (`active = (u_pos + u_neg) > 0`, otherwise it writes zeros).
- **Category:** dead. **Severity:** low.
- **What is wrong:** a slot with spliced crossings and no unspliced count gets a full (λ, θ) cube and then reports the zeros it would have had unsolved. The sweep's `in_slots` (`solve_kernel.cpp:400-403`) admits it too, and then `solvable` (n_u > 0) discards the result: two wasted solves per such slot per sweep, three under a capture. `mass_unspl` duplicates `u_pos + u_neg` at the one production call (`region_geometry.py:523`). The docstring's "`mass_spliced` enters only the signal mask" describes an input with no effect.
- **Proposed change:** mask on `u_pos + u_neg > 0`, drop both mass arguments, and prove it a no-op with `rename_identity.py`.
- **Moves numbers:** no (the wasted time is unmeasured). **Confidence:** confirmed.

**6. Production-dead surface not on the ledger**
- **Where:** `gdna_density.py:124` `def pooled_log_rate(...)`; `priors.py:36-38` `_RNA_SIGNATURE_BITS`; `splice_graph.py:1124,1130` `is_terminus` / `is_splice_site`; `capture_eff_length.py:121`.
- **Category:** dead. **Severity:** low.
- **What is wrong:**
  - `pooled_log_rate` is called only by `test_gdna_density.py`. `_layers.py:57-59` says both rate estimators are "shared by their consumers", and `one_sided_rate` spells the pooled rate inline anyway (`gdna_density.py:173`).
  - `_RNA_SIGNATURE_BITS` is unused in `src`; its own comment names a test as its one consumer.
  - `is_terminus` / `is_splice_site` have no production caller (the C++ carries its own bits, `transfer_kernel.h:38-48`), yet `region_geometry.py:399-401` tells consumers to compose with them.
  - `capture_eff_length.py:121`'s `tdf is not None and "is_synthetic" in tdf.columns` is a format guard contradicted by its own line 151 (`tdf.shape[0]`) and by `sj_opportunity.py:157`.
- **Proposed change:** move these into tests (or ledger item (c)) and correct both docstrings.
- **Moves numbers:** no. **Confidence:** confirmed.

**7. Rules written twice (or three times)**
- **Where:** `priors.py:95-96` against `region_arrays.py:205-206`:
  ```
  first = lo + np.searchsorted(ends[lo:hi], b_start, side="right")
  last = lo + np.searchsorted(starts[lo:hi], b_end, side="left")
  ```
- **Category:** duplicate. **Severity:** low.
- **What is wrong:**
  - `_region_locus_shares` inlines `overlapping_region_runs`.
  - `gdna_opportunity_from_index` (217-223) re-implements `region_lengths_from_partition`.
  - "Intergenic region" is spelled three times: `transfer_kernel.h:148-150`, `transfer.py:149`, and `sweep.py:319` (`terminal`).
- **Proposed change:** call the shared helpers; one predicate.
- **Moves numbers:** no. **Confidence:** confirmed.

**8. The landscape's tiling reaches its bits; bare literals beside it**
- **Where:** `landscape.py:291-297`:
  ```
  rows = _block_rows(grid.size, 8)
  ...
  sums[b] += (weights[sl][m][:, None] * kernels[m]).sum(0)
  ```
- **Category:** bare-constant + contract. **Severity:** low.
- **What is wrong:**
  - Tile partial sums are added across tiles, so the float association depends on the tile size. `config.py:794-796` says the block size "cannot reach the arithmetic", and the landscape has no tiling gate (capture efficiency has one).
  - `_block_rows`' second argument is a bare 8 here and 16 at `capture_efficiency.py:139`. It is called `itemsize` in the source and `n_arrays` in `test_capture_efficiency.py:306`, and neither value counts the 3–5 temporaries per cell.
  - Also bare: `landscape.py:133` `lo - 0.5, lo + 0.5` and `:189` `max(..., 2)`.
- **Proposed change:** make the accumulation tile-invariant (or add a gate and correct the docstring); derive the per-cell bytes; document or derive the two literals.
- **Moves numbers:** no at the current value; ULPs if the constant changes. **Confidence:** confirmed by reading, not run.

**9. Comments describing old structure**
- **Where:** `region_geometry.py:23-25` and `264-265`; `region_chain.py:74-77` and `99`; `simplex_logodds.py:9-11` and `89-90`; `result.py:202-203`.
- **Category:** stale-comment. **Severity:** low.
- **What is wrong:**
  - "densities … computed inline in `sweep.solve_chain`": the sweep is one native call.
  - `RegionStatics` is said to hold "per-strand counts"; it holds none.
  - `region_gdna_geometry` claims one shared definition, but the kernel rebuilds `n_u` itself (`solve_kernel.cpp:392`).
  - "the strand family (layer 4)" no longer exists.
  - The masses are described as what "calibrate/derive read ONLY", but `_result` also publishes the fractions.
  - `simplex_logodds` mentions a fitted "RNA arm" / `logP_r` that no code has.
  - The od fields are called "fitted", but they are always 0.0 and only tests read them.
- **Proposed change:** rewrite these lines; relabel or drop the od fields.
- **Moves numbers:** no. **Confidence:** confirmed.

**10. "Mature" and "nascent" used as populations in comments**
- **Where:** `simplex_logodds.py:48` ("at a sj mature RNA *splices*"); `splice_graph.py:36` ("crossed by gDNA, by nascent RNA and by mature RNA alike"), `:49`, `:56`, `:849/862/876`; `fl.py:22`, `:481`; `sj_opportunity.py:153`; `region_geometry.py:78`. The field `blocks.SweepCapture.mature` is covered in finding 3.
- **Category:** naming. **Severity:** low.
- **Proposed change:** reword as "RNA that splices / has not spliced there" and "annotated introns".
- **Moves numbers:** no. **Confidence:** confirmed.

## Already recorded

- **`hygiene-ledger`:**
  - (c) `has_own_composition_evidence`, `_project_regions_to_loci`, `GdnaContrast` / `GdnaRealized`.
  - (d) the vertex claim, `sweep`'s UNIDENTIFIED deferral, "no per-region spliced floor", and `fl`'s "gDNA by structure" (`fl.py:16-17,101`).
  - (e) `density_model.count_observable_masks` (its mature/nascent wording goes with it) and the log-odds window gate.
  - (f) three unread scalars in `solve_kernel.cpp`, the `-g` build, the lgamma machinery.
  - (h) `solve_one_block` / `transfer_prepare` wiring, the two `contained_opportunity` functions and the ramp written twice, `intervals.feather` read three ways.
  - (i) native constants and inline C++ literals, e.g. `blur_radius`'s 4σ. **Correction to (i):** `transfer_rows.h:16-17` `EPS` / `TINY` carry no comment, and `EPS` serves as a variance floor, a flatness threshold in nats, and a fraction clamp.
- **`latent-defects`:** `level_bound_row`'s 1e-12, nanobind `opt()`, `run_parallel`, `_ref_slices`, `fast_exp`.
- **`the-transfer-price-arithmetic`:** empty-piece witness units, `f_pred`, `hop_price` frames.
- **`the-gdna-prior-enters-psi-twice`:** includes `simplex_logodds`' "replaces" vs "added".
- **`strand-overdispersion-one-shared-value`:** `strand_variance`'s (n·f)².
- **`psi-reads-kappa-where-the-strand-channel-is-dead`.**
- **`debug-capture-memory-is-unbounded`:** every sweep is captured, and `ruler_vs_truth` always passes `_debug`. Finding 3 is the unread-field half of this.
- **`an-anchor-only-landscape-reads-as-capture`:** `_render` / `census`.
- **TRAPS:** `mature-rna-never-crosses-a-boundary` (the domain wording behind `fl`'s comments).
- **Also mention files in this slice, as design entries rather than hygiene:** `capture-blind-gdna-divisor`, `the-junction-price-is-noisy-within-a-gene`, `the-crossing-apportionment`, `the-efficiency-posterior-floor-on-empty-pieces`, `the-reader-reads-false-gdna-as-capture`, `the-fl-boundary-inversion-*`, `the-realized-gdna-length-law-reads-rna-counts`, `the-gdna-length-law-falls-back-at-identical-purities`, `theta-quadrature-at-zero-gdna`, `f32-strand-tilt-at-half`, `transfer-variance-premise`, `flux-source-skipped-at-an-empty-exon-piece`, `the-pooled-q-in-the-gdna-count`.
