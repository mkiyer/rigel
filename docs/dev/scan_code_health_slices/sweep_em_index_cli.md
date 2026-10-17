## Scope read

I read every file in the slice in full: `native/em_solver.cpp` 1–2950, `estimator.py` 1–888, `index.py` 1–1780, `cli.py` 1–1575, `config.py` 1–878, `frag_length_model.py`, `strand_model.py`, `gtf.py`, `transcript.py`, `types.py`, `junction_fit.py`, `native.py`, `report/{__init__,build,substrate,model,specs,html}.py`, and the headers `simd_detect.h`, `thread_pool.h`, `thread_queue.h`, `ndarray_util.h`. I read `fast_exp.h` 1–80 closely and skimmed the rest. In `report/assets/report.js` I read 262–282 and 385–410.

To check claims I also read: `pipeline.py` 540–575, 715–935, 1000–1016; `scan_cache.py` 296–310; `native/resolve_context.h` 670–860; `calibration/fl.py` 78–214; `locus.py` 150–230; `CMakeLists.txt` (the flags); and the test fixtures I cite. From the docs I read ISSUES `hygiene-ledger`, `latent-defects`, `per-transcript-prior-lane`, `overlapping-synthetic-shadows`, `sj-strand-tag-chosen-from-the-first-reads`, `the-em-stops-before-it-converges` and `the-format-changes-to-batch-before-release`, plus DESIGN §5 and TRAPS `nrna-does-not-mean-synthetic`.

I found no doc citations in the slice and no Greek letters in any identifier. Every config field and every `CONSTANTS` field has at least one `src/` reader, except `unseen_smoothing_ess` (finding 11).

## Findings

**1. The outputs use "nascent" for two different sets of fragments, and the report mislabels its own row.** `src/rigel/report/model.py:298-304`
```
        ["Mature mRNA", mrna, q.get("mrna_fraction", 0.0), "spliced + EM-assigned transcript mass"],
        [
            "Nascent RNA",
            nrna,
            q.get("nrna_fraction", 0.0),
            "synthetic + annotated single-exon spans",
```
- **Category:** contract · **Severity:** medium · **Confidence:** confirmed · **Moves numbers:** no (labels and column names; the goldens' headers change if renamed).
- **What:** `nrna_total` is `estimator.nrna_em_count`, which counts synthetic spans only (`estimator.py:229-234`). Annotated single-exon transcripts are counted in `mrna_total = quant_df["count"].sum()` (`cli.py:433-434`), so the row's description is false. `gene_quant` does the reverse: `nascent_count` holds exactly the annotated single-exon (`is_nrna`) transcripts and leaves synthetic spans out (`estimator.py:564-575`). So one fragment on a single-exon gene is "Mature mRNA" in `summary.json`, `loci.feather` and the report, and "nascent" in `gene_quant`. Both splits also treat mature and nascent as populations, against Axiom 0: `estimator.py:460` ("they are both mature and nascent"), `cli.py:342`, `html.py:146`, `report.js:277-278`.
- **Evidence:** read the two masks. No test asserts the row text: `rg "synthetic + annotated"` finds only `model.py`.
- **Change:** fix the row description now ("synthetic spans"). The names (`mature_count`/`nascent_count`, `mrna_*`/`nrna_*`, the "Mature mRNA · nascent RNA" header) are an output contract recorded in the MANUAL and the goldens, so renaming them needs an owner ruling.

**2. The `--em-mode map` help describes an update the code does not run.** `src/rigel/cli.py:1319-1321`
```
        help=f"EM algorithm variant (default: {qd['em_mode']}). "
        "'vbem' uses Variational Bayes EM with digamma-based soft "
        "updates. 'map' uses MAP-EM with hard max(0, n+a-1) updates.",
```
- **Category:** stale-comment · **Severity:** medium · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** `map_em_step` (`em_solver.cpp:774-791`) adds the pseudocounts to `unambig + em_totals` through `apply_grouped_prior_update`, then normalises. The MAP path has no `−1` and no `max(0,·)`, so it never does the sparse Dirichlet-MAP update the help promises. `docs/MANUAL.md:193` repeats the formula.
- **Evidence:** `rg "n\+a-1|max\(0" src/rigel/native/em_solver.cpp` finds nothing.
- **Change:** describe MAP as EM with the calibration pseudocounts added to the counts, with no −1.

**3. Most index load failures crash `rigel quant` with a traceback instead of the rebuild message.** `src/rigel/index.py:1323-1326`
```
        for name, path in ((REGIONS_FEATHER, regions_path), (BOUNDARIES_FEATHER, boundaries_path)):
            if not os.path.exists(path):
                raise RuntimeError(
                    f"Index at {index_dir} is missing {name} (the splice graph). "
```
- **Category:** contract · **Severity:** low–medium · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** The `IndexLoadError` docstring (169-171) promises a rebuild message for "a missing file", and `quant_command` catches only `IndexLoadError` (`cli.py:132-135`). But:
  - a missing `regions.feather` or `edges.feather` raises a plain `RuntimeError`;
  - a missing `transcripts`, `intervals` or `sj` feather raises pyarrow's `FileNotFoundError` (unguarded reads at 1276, 1335, 1473);
  - a `t_index` mismatch raises `ValueError` (1277-1281).
- **Evidence:** delete `regions.feather` from any index and run `rigel quant`.
- **Change:** check every required file up front and raise `IndexLoadError` for each failure.

**4. The docs for `splicing_anchor_tolerance` say calibration is unaffected; it is not.** `src/rigel/config.py:219-221`, `cli.py:1452-1455`
```
    The calibration accumulator does not read this value: compartment,
    splice and strand are recorded directly into the per-region
    channels. It tunes the resolver alone.
```
- **Category:** stale-comment · **Severity:** low–medium · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** K decides which transcripts count as contradicted and which introns a mate gap may hold (`resolve_context.h:797-818`). Those gap hypotheses are what the accumulator arbitrates and deposits (720-721). So K changes calibration's tally, which is why `scan_cache._scan_config_digest` (`scan_cache.py:301-309`) includes it in the cache key as a field that "can change a tally". The CLI help also says "the fractional calibration accumulator", a name that appears nowhere else.
- **Change:** say that K shapes the hypotheses the accumulator arbitrates, so it changes calibration and the scan-cache key.

**5. Possibly out of scope: the anchor slack exists to absorb aligner overrun.** `src/rigel/native/resolve_context.h:694-696`, reached from `--splicing-anchor-tolerance` and `BamScanConfig.splicing_anchor_tolerance`
```
    /// SPLICING ANCHOR TOLERANCE K before it counts -- the same slack the gap match uses, because it is
    /// the same physical thing: an aligner routinely runs a read a base or two past a donor. Without
    /// this the predicate rejects exactly the fragments K exists to accept.
```
- **Category:** out-of-scope · **Severity:** low · **Confidence:** suspected · **Moves numbers:** yes, if it is ever removed (the default 3 would become 0).
- **What:** By its own comment, K's purpose is to accept up to K aligned bases that lie inside an annotated intron, as if they had been spliced out. That overrides what the aligner wrote. Whether that counts as alignment repair under the 2026-10-16 scope ruling is the owner's call. I do not propose deleting it.

**6. The strand trainer's log line prints two different estimates of the same rate side by side.** `src/rigel/strand_model.py:363-366`
```
        alpha = self.n_minor + 1.0
        beta = (n - self.n_minor) + 1.0
        # UCL on minor-orientation rate r = 1 − ss at the given confidence.
        r_ucl = float(betaincinv(alpha, beta, confidence))
```
- **Category:** duplicate · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no (log only).
- **What:** `n_minor` comes from the pooled 2×2 table (252-255), while `strand_specificity` comes from the genuine-junction κ (266-295). The docstring calls the result an upper limit on `1 − strand_specificity`, but it is the limit on the pooled rate, splice artifacts included. `pipeline.py:1011-1016` logs `ss_est` and `ε_CI` together. On cfRNA-B (pooled 0.0635, genuine 0.0030, `ISSUES: strand-overdispersion-one-shared-value`) that prints ss ≈ 0.997 next to an "upper limit" of about 0.065.
- **Change:** compute the Beta posterior from `genuine_n_same`, or label the value as the pooled rate.

**7. The EM's mode codes are written out on both sides of the binding.** `estimator.py:351,372`; `em_solver.cpp:55-57,880`
```
        _ASSIGNMENT_MODE_MAP = {"fractional": 0, "sample": 1}
            {"coverage": 0, "prior": 1, "uniform": 2}[self.em_config.warm_start],
// Assignment mode constants (must match Python _ASSIGNMENT_MODE_MAP)
static constexpr int ASSIGN_FRACTIONAL = 0;
```
- **Category:** duplicate · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** `config.py:13-14` says constants shared with Python are defined once in C++ and exported, and hygiene-ledger (i) says such values are never restated. These two enums are restated in Python. In C++ the warm-start codes are bare `1`/`2` (901, 906). The binding checks neither: an assignment code other than 0 or 1 indexes the empty `choice` vector at `em_solver.cpp:1618`.
- **Change:** export the codes from `_em_impl` as `EM_LOG_EPSILON` is, read them in Python, and refuse unknown codes in the binding.

**8. The `summary.json` schema version is defined twice.** `cli.py:208` (`SUMMARY_SCHEMA_VERSION = 3`) and `report/substrate.py:25` (`EXPECTED_SCHEMA_VERSION = 3`).
- **Category:** duplicate · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** Bumping the writer's value alone puts the drift banner on every new report. The only gate (`tests/test_summary_report.py:192`) checks the writer's value.
- **Change:** have the report read the writer's constant.

**9. `IndexConfig` is never constructed, so its validation never runs.** `config.py:439-451`, `cli.py:1164-1173`, `index.py:981`
- **Category:** dead · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** `rg "IndexConfig\("` over `src/`, `tests/` and `scripts/` finds nothing; only class-attribute defaults are read. So `rigel index --nrna-tolerance -5` builds without complaint (it clusters exactly as tolerance 0 would) and records −5 in the manifest. The same setting also has two names: `nrna_merge_tolerance` in the config and `nrna_tolerance` in the CLI, the build call and the manifest.
- **Change:** build an `IndexConfig` from the CLI arguments and pass it to `build` under one name.

**10. Members that no production code reaches.**
- `TranscriptIndex.g_to_strand_arr` (`index.py:854,1310`): read only by `tests/test_index_integrity.py:105-108` and a test mock.
- `AbundanceEstimator(is_nrna=…)`, stored as `_nrna_mask` (`estimator.py:163,171-175`): read only by `tests/test_opp_strand_order.py:152`.
- `IntervalType.SJ` (`types.py:111`): nothing reads it. `sj.feather` rows carry no interval type, and the native side defines only EXON and TRANSCRIPT. Its docstrings (`types.py:122-125`, `index.py:600`) describe rows that are never built.
- `index.write_bed12` (`index.py:610-655`): called only by `sim/scenario.py:566-582` and `tests/test_sim.py`.

**Category:** dead · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no. **Change:** delete the first three; move `write_bed12` into `rigel.sim`.

**11. `FragmentLengthModel` has an unreachable smoothing path and a redundant flag.** `frag_length_model.py:251-257`
```
        total = self.total_weight
        if total == 0:
            n = self.max_size + 1
            return np.full(n, 1.0 / n, dtype=np.float64)
        n = self.max_size + 1
        ess = CONSTANTS.fragment_length.unseen_smoothing_ess
        return (self.counts + ess / n) / (total + ess)
```
- **Category:** dead / duplicate · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.
- **What:**
  - Every caller of `pmf` and `compute_all_transcript_eff_lens` (`pipeline.py:556,561` and every test) passes a model built by `from_pmf`. The raw view that `fl._empirical` builds is only ever serialised with `to_dict`. So this branch and `CONSTANTS.fragment_length.unseen_smoothing_ess` are reached by nothing.
  - `from_pmf` sets `counts[:] = prob` and `_total_weight = 1.0` (119-125), so `_stats_use_prob`, `_stats_probs` and `_stats_weights` choose between two equal vectors.
  - Two docstrings describe the dead path: 144-145 and 233-234.
- **Change:** make the raw view statistics-only; delete the branch, the constant and the flag.

**12. The estimator's output tables recompute values they already hold.** `estimator.py:861-863`
```
            rna_total = float(r["rna_total"])  # EM: annotated mRNA + synthetic nRNA
            mrna_em = max(rna_total - nrna, 0.0)  # EM-only annotated mRNA
            mrna = mrna_em + unambig_mrna  # total mRNA = EM + unambig
```
- **Category:** duplicate · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** yes, by ulps, in `loci.feather` only (item a).
- **What:**
  - (a) The loci table rebuilds `mrna` by subtraction, behind a clamp that hides any drift, although `t_total[annot_in_locus].sum()` is the same quantity.
  - (b) `get_nrna_counts_df` recodes the TPM denominator inline (743-745); `_tpm(t_total, eff)[idx]` gives the same, bit-identical result.
  - (c) The `is_synthetic` fallback at 608-612 exists only for the test `MockIndex` (`tests/test_estimator.py:83-92`). The comment at 604-607 says `annotated_mask` zeroes synthetic rows in "all aggregates above", which is false for `count_unambig`, `count_spliced` and the length sums.
- **Change:** read the values already held; add `is_synthetic` to the mock and index the column strictly.

**13. QC and display constants live outside `CONSTANTS.qc`.**
- `report/model.py:58,66-68,91,106`: the thresholds that colour the headline QC tiles (mapping rate, strandedness, gDNA fraction, genic fraction); also `max_rows=20000` at 335.
- `strand_model.py:160-161` and `488-489`: the junction-depth cuts 100 and 1000, written twice.
- `index.py:152,322-323`: error messages list 3 examples, while `CONSTANTS.qc.report_examples = 5` serves the same purpose at 252.
- `cli.py:353-355,583`: `round(…, 2)` and `round(…, 8)` beside `SUMMARY_DECIMALS`.
- `report/model.py:7-9` says the module "emits raw numbers", yet it rounds at 275, 374-376 and 453-455.

**Category:** bare-constant · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.

**14. Simulator inputs sit on the production `Transcript`.** `transcript.py:86-89`
```
    nrna_abundance: float = 0.0
    """Nascent RNA (pre-mRNA) molecular abundance.

    Independent of ``abundance`` (mature mRNA molecular abundance).
```
- **Category:** naming · **Severity:** low · **Confidence:** confirmed · **Moves numbers:** no.
- **What:** `abundance` (set from the GTF score column at 109) and `nrna_abundance` are read only by `sim/scenario.py:313-316`. Outside `rigel.sim`, their docstrings use mature and nascent as populations.
- **Change:** move them to the simulator's own record.

**15. Stale comments** (all stale-comment · low · confirmed · moves numbers: no)
- `types.py:81` and `index.py:1666` refer to a `Fragment` class that no longer exists anywhere in the repository.
- The file list in `index.py`'s module docstring omits `transcript_aliases.feather`, and `index.py:840-842` says only transcripts, intervals and junctions are stored on disk.
- `em_solver.cpp:1890-1892,2047-2049,2081` call the EM's lengths "FL-marginal"; they are the capture-contracted lengths.
- `em_solver.cpp:141-147`: digamma's `x <= 0` branch is unreachable, because both callers clamp at `EM_LOG_EPSILON`; its two returns are identical and its comment is wrong.
- Smaller residue: an orphaned comment at `cli.py:878-881`; "25+ individual keyword arguments" at `config.py:463-464`; `config.py:471` opens a doc comment with `#` instead of `#:`.

## Already recorded

- `hygiene-ledger`:
  - (c) `_cgranges_impl` is kept for the test-only `TranscriptIndex.query()`.
  - (e) The report test's v3 summary omits the blacklist keys.
  - (f) `--no-mappability` naming.
  - (g) `quant_command`'s docstring omits the annotated-BAM pass.
  - (h) The SQUAREM VBEM/MAP duplication, `load`'s frozenset round-trip, and `intervals.feather` read three ways.
  - (i) `SQUAREM_BUDGET_DIVISOR` and `em_solver.cpp`'s unaudited inline literals (`MAX_K_STACK`, the restated `1e-300`/`-1e300`); `rigel sim` restating `Scenario`'s defaults; `FragmentLengthModel.to_dict`'s inline rounding.
- `latent-defects`:
  - The ~1e-300 posterior a dead component keeps in `assign_posteriors`.
  - `summary.json`'s `n_transcripts`, `read_groups` and `mapped_reads`, and the report's read-fate bar.
  - Sample mode's `posterior_mean` is Σp²/Σp.
  - A `gene_id` with transcripts on two references (`_build_gene_table`).
  - `run_parallel` and `ThreadPool` throws.
  - `nb::cast<A>(h).data()` in `em_solver.cpp`.
  - `fast_exp` on an all −inf row.
- `per-transcript-prior-lane`: `rna_prior_weight` and `warm_start='prior'`, ruled kept.
- `sj-strand-tag-chosen-from-the-first-reads`: the `auto` detection described in the `--sj-strand-tag` help.
- `overlapping-synthetic-shadows`: per-span `nrna_quant` rows and `nrna_id`.
- `the-format-changes-to-batch-before-release`: the report's capture note for older (schema-2) summaries, the one backward-compatibility path in `report/`.
- `the-em-stops-before-it-converges`: the EM iteration budget.
- `TRAPS: nrna-does-not-mean-synthetic`: the meaning of `is_nrna` (background to finding 1).
- Not in ISSUES or TRAPS but pinned by a test: the `n_observations <= 0` guard can never fire on a `from_pmf` model (`tests/test_d7_transcript_eff_lengths.py::test_the_n_observations_GUARD_ON_THAT_PATH_CANNOT_FIRE`).
