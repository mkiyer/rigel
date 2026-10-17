## Scope read
Read in full: `src/rigel/native/scoring.cpp` (1–1158), `src/rigel/scoring.py` (1–161), `src/rigel/scored_fragments.py` (1–120), `src/rigel/locus.py` (1–233), `src/rigel/locus_partition.py` (1–109).

Read for context:
- `scan.py` (all); `pipeline.py` 560–960
- `em_solver.cpp`: 170–295 (equivalence classes), 870–935 (warm start), 1880–2010 (sub-problem), 2112–2135, 2560–2700 (connected components)
- `constants.h` 1–120 and 300–400; `buffer.py` 120–190; `frag_length_model.py` 1–140; `config.py` 111–157 and 700–737; `splice.py` 1–100; `strand_model.py` 284–300; `annotate.py` 60–110 and 425–460; `index.py` 1715–1790
- `resolve_context.h`: 180–195, 920–985, 1048–1100, 1280–1365
- `bam_scanner.cpp`: 1355–1400, 1840–2030, 2408–2430, 2600–2700, 2780–2860
- the length filter in `calibration/accumulator.{h,cpp}`
- tests: `tests/test_pipeline_routing.py` (all), `tests/test_locus.py` 1–110
- docs: ISSUES entries hygiene-ledger, latent-defects, scoring-penalties-are-underived-constants, multimapper-intergenic-alignments, binary-cuts-on-continuous-quantities, the-minus-strand-coverage-weight-sat-one-fragment-length-off; DESIGN's nascent ruling and §3.1d; the TRAPS rule nrna-does-not-mean-synthetic.

Checks were pure-Python transcriptions, with no rigel import, in the scratchpad (`<scratch>/`): `lse_check.py`, `cov_check.py`, `cov_weight.py`. I edited nothing in the repo.

## Findings

### 1. Under `--mismatch-alpha 0`, a multimapper can lose its gDNA term depending on the order of its hits
`scoring.cpp:57-73, 519-527, 646-653`
```
    } else {
        sum_v += std::exp(x - max_v);
...
                st.mm_gdna_lse_sum > 0.0) {
```
- **Category / severity:** defect · medium
- **Moves numbers:** yes, only at the gate
- **Confidence:** confirmed
- **What is wrong:** At the gate, `log_nm = nm * -inf` gives a mismatched hit a gDNA term of −inf. If the first two gDNA-explainable hits are both −inf, the update computes exp(−inf − −inf) = NaN. NaN then persists, `NaN > 0` fails, and the group is certified RNA, even though a later NM = 0 hit is a valid gDNA reading. The comment at 169–170 ("nm * -inf is -inf … the EM path needs no flag") is false here.
- **Evidence:** In `lse_check.py`, the order (clean, gated, gated) gives −7.10; the order (gated, gated, clean) drops gDNA. No test applies the mismatch gate to a multimapper.
- **Proposed change:** Skip a gated hit in both the RNA merge and the gDNA term, as the overhang gate already does. Make `lse_update` ignore −inf. Add a test over the three hit orders.

### 2. The multimapper merge keeps each transcript's least-overhang hit, not its most likely one
`scoring.cpp:569-582`
```
                if (oh < prev.overhang ||
                    (oh == prev.overhang &&
                     log_lik > prev.log_lik))
```
- **Category / severity:** defect / duplicate · medium
- **Moves numbers:** yes
- **Confidence:** confirmed
- **What is wrong:** Overhang is already inside `log_lik` (`oh * oh_log_pen_`), and it is used again as the primary key. At the default alphas, a hit with overhang 0 and NM 3 beats one with overhang 1 and NM 0, so the transcript is scored 4.6 nats below its best hit. Under the mismatch gate, the transcript is lost outright: the −inf entry is never replaced by the clean, more-overhanging hit, and is then pruned. `count_col` and the coverage weight follow the chosen hit. No doc rules on which hit a transcript takes.
- **Evidence:** Every multimapper test in `test_pipeline_routing.py` puts its hits on distinct transcripts, so the replace branch (573–581) is never exercised.
- **Proposed change:** Choose the hit by `log_lik`, or state and derive another rule. Add a test with two hits on one transcript, at the default alpha and at the gate.

### 3. Pruning uses three different references
`scoring.cpp:843-868` (single fragments), `604-629` (multimappers), `883-888` (gDNA)
```
                double pool_best = t_is_nrna_[s.t_idx]
                                 ? nrna_best : mrna_best;
...
                    has_gdna = gdna_competes(gdna_ll, best_ll);
```
- **Category / severity:** defect / duplicate · medium
- **Moves numbers:** yes; size not measured
- **Confidence:** confirmed by reading
- **What is wrong:** An annotated candidate is pruned against the best annotated candidate, a synthetic span against the best span, and gDNA against the best of both. The span and gDNA read the same footprint; with equal length laws they differ only by log(2·p_sense). Yet the span survives any distance below the annotated reading, while gDNA is cut at 9.21 nats. So an IMPLICIT unit with a long gap can go to the annotated transcript or the span, but never to gDNA. DESIGN §3.1d's "competes exactly as an RNA candidate does" holds for neither pool, and `FragmentScoringConfig`'s "its fragment's best" is inexact.
- **Evidence:** Take one unit with the annotated transcript at A, the synthetic span at A−12 and gDNA at A−12.7. The span is emitted, gDNA is dropped, and `is_spliced` is 1. The routing tests set `is_synthetic` to False everywhere, so pool separation is untested.
- **Proposed change:** Use one reference, the unit's best candidate of any kind, for every candidate. A/B it per stratum.

### 4. The coverage weight's plateau is wrong for transcripts shorter than two fragment lengths, and its gate restates the formula
`scoring.cpp:93-96, 111`; `tests/test_pipeline_routing.py:1004`
```
    double w = std::min(static_cast<double>(f_len), L / 2.0);
```
- **Category / severity:** defect + test-gap · medium
- **Moves numbers:** yes, through the EM warm start, so only where the EM depends on its start
- **Confidence:** confirmed
- **What is wrong:** Enumerating uniform placements gives coverage min(x, f, L−f, L−x), so the plateau is min(f, L−f). This equals min(f, L/2) only when L ≥ 2f. For L = 300 and f = 200 at offset 50, the code gives 1.50 against 1.14 enumerated; for L = 260 and f = 250 it gives 1.93 against 1.02. `_trapezoid_weight` uses the same `min(f, length / 2)`, so the `_SHORT` cases (L = 300, f = 200) pass on the wrong model. Separately, w is taken from the cut interval rather than the fragment's length.
- **Evidence:** `cov_check.py` and `cov_weight.py`.
- **Proposed change:** Derive w from the placement count, with weight 1 when f ≥ L. Gate against enumeration rather than against the formula.

### 5. `is_spliced` really means "no gDNA candidate", duplicates `isfinite(gdna_log_liks)`, and its docstring states a replaced rule
`scored_fragments.py:48-51`
```
        bool[n_units] — True if gDNA cannot explain this unit: certified RNA (the scorer's
        ``gdna_can_explain``). An artifact or an implicit splice within the maximum fragment length is
        False; a sequenced junction, or an implicit splice beyond that length, is True.
```
- **Category / severity:** naming + stale-comment + duplicate · medium
- **Moves numbers:** no
- **Confidence:** confirmed
- **What is wrong:** The field is set to `!has_gdna`, so it is also True for an UNSPLICED or ARTIFACT unit whose gDNA term was pruned. `test_an_ARTIFACT_whose_footprint_the_law_cannot_produce_is_pruned` asserts exactly that, contradicting the docstring. The "maximum fragment length" clause describes the rule DESIGN §3.1d says pruning replaced. The field is the same bit as `gdna_log_liks == -inf`, which is why both readers check both: `em_solver.cpp:1959` and `pipeline.py:770`. Separately, the `log_liks` docstring ("log(P_strand × P_insert)") omits the overhang and NM terms.
- **Proposed change:** Derive the field from `isfinite(gdna_log_liks)` or rename it, and fix both docstrings.

### 6. A multimapper's reported splice type drops IMPLICIT and ARTIFACT
`scoring.cpp:505-511`
```
        if (stype == SPLICE_SPLICED_ANNOT) {
            st.mm_best_stype = SPLICE_SPLICED_ANNOT;
        } else if (stype == SPLICE_SPLICED_UNANNOT &&
```
- **Category / severity:** defect (reporting) · low
- **Moves numbers:** no; only the ZS tag
- **Confidence:** confirmed
- **What is wrong:** A group made only of IMPLICIT or ARTIFACT hits reports UNSPLICED. The writer stamps that type as ZS on every hit (`bam_scanner.cpp:2652`). `n_calibrated` excludes multimappers by class, so no count moves.
- **Evidence:** `test_a_multimapper_of_IMPLICIT_and_ARTIFACT_alignments_is_a_gdna_candidate` builds exactly this group and asserts nothing about `splice_type`.
- **Proposed change:** Rank all five splice types and assert the result in that test.

### 7. `pruning_min_posterior` is never validated
`scoring.py:129-130`
```
        _eps = max(pruning_min_posterior, CONSTANTS.scoring.pruning_floor)
        max_ll_delta = -math.log(_eps) if _eps < 1.0 else 0.0
```
- **Category / severity:** contract · low
- **Moves numbers:** no, for valid input
- **Confidence:** confirmed
- **What is wrong:** Both alphas are range-checked; this CLI flag is not. A value ≥ 1, or NaN (`max(nan, x)` returns nan), silently sets the margin to 0, so each unit keeps only its tied-best candidate. A negative value disables pruning. The name is also misleading: it bounds a likelihood ratio, not a posterior.
- **Proposed change:** Check 0 ≤ value < 1 in `__post_init__`.

### 8. Mature/nascent vocabulary, and an array named `t_is_nrna` that holds `is_synthetic`
`scoring.py:132-136`; `scoring.cpp:181, 398, 533, 605-611, 793-800`; `locus.py:109`
```
        t_is_nrna_arr = index.t_df["is_synthetic"].values.astype(np.uint8)
```
- **Category / severity:** naming · low
- **Moves numbers:** no
- **Confidence:** confirmed
- **What is wrong:** Identifiers and comments treat mRNA and nRNA as populations (`MergedMrna`, `MrnaScored`, `mrna_best`, `nrna_best`, "1=nRNA, 0=mRNA", "mRNA + nRNA candidates"). The array named for nRNA actually holds the synthetic flag, which is the confusion `TRAPS: nrna-does-not-mean-synthetic` records.
- **Proposed change:** Rename to `t_is_synthetic`, `annotated_best` and `synthetic_best`, then run `rename_identity.py --check`.

### 9. `MultiLocus.gdna_span` is just a reported span, and its floor does nothing
`locus.py:228` (`gdna_span=max(span, 1)`); `tests/test_locus.py:4, 102`
- **Category / severity:** naming + stale-comment · low
- **Moves numbers:** no
- **Confidence:** confirmed
- **What is wrong:** The only reader is `pipeline.py:830` (`locus_span_bp`), and nothing divides by it, so "floored at one so nothing downstream divides by zero" is stale. Any transcript with end > start gives a span ≥ 1. The name also suggests a gDNA length, which the EM takes from `priors.gdna_eff_len` instead.
- **Proposed change:** Rename to `span_bp` and drop the floor and its assertion.

### 10. Residue
All low, no numbers move, confirmed.
- **Dead:** `stat_chim_` is returned by `finish()` and discarded as `_` (`scan.py:129`).
- **Bare constants:** `exon_str == 1 || exon_str == 2` at `scoring.cpp:501, 758, 785`, although `STRAND_POS` and `STRAND_NEG` exist.
- **Stale comments:**
  - "Track member count for NH bookkeeping" (`:530`): the count is only a non-empty flag.
  - "called from the Python scoring loop" (`:4-5`).
  - "g_offsets must stay alive" (`locus_partition.py:52`): no such name exists.
  - "replaces a 457K-iteration Python loop" (`scoring.py:102-103`).
- **Vacuous test:** `test_multimapper_spliced_annot_skips_shadows` only asserts `t_indices < num_transcripts`, which the scorer cannot violate.

## Contracts
- **`NativeFragmentScorer` ctor:** takes `n_transcripts` from `t_strand_arr`. It assumes the length, synthetic-flag and exon-offset arrays, and a length table of `fl_max_size + 1` entries, without checking any of them; a short array is an out-of-bounds read.
- **`StreamingScorer` ctor:**
  - `unambig_counts` must have at least 6 columns, or deterministic counts are dropped silently (`col < ua_ncols`).
  - The transcript strand arrives twice: once in the scorer for EM units, once as this pointer for deterministic counts.
- **`score_chunk`:**
  - N comes from the splice-type array. Nothing checks that `t_off` has N+1 int32 entries or that the candidate arrays have `t_off[N]`.
  - A multimapper's rows must be contiguous under one frag_id. This holds today because the scanner cuts chunks between read-name groups.
- **Chimeric branch:**
  - A chimeric hit of a multimapper is annotated as chimeric under the group's frag_id while its other hits form an EM unit. Only the writer's last-wins `frag_to_row` hides the duplicate.
  - If that group is then gated, the fragment stays marked chimeric, against the scanner's rule that a fragment is chimeric only if all its hits are.
- **Length terms:**
  - A length ≤ 0, or a missing length table, scores 0 (log 1), which is the best possible score.
  - Production never reaches this, but the routing tests' unfinalized models rely on it, as they do on `getattr(..., "_tail_base", 0.0)`.
- **`coverage_weight`:**
  - Assumes each transcript's exons are sorted and disjoint and that `t_length` equals the sum of its exon lengths; neither is checked.
  - It restates the resolver's `tx_frag_length` start rule.
- **`stat_gated`:** also counts units whose every candidate has `exon_bp == 0` at any alpha, which the MANUAL's definition of `n_gated_out` does not cover.
- **`build_multi_loci`:** an out-of-range candidate index silently leaves the unit out of every locus. The BAM-order sort assumes frag_ids are unique.
- **`partition_and_free`:**
  - Keys partitions and `locus_id` by list position, while `pipeline` looks them up by `multi_locus_id` and stamps ZL from `locus_id`. They agree only because `build_multi_loci` numbers loci by position.
  - `em_data` is unusable afterwards.
- **Hash-map order:** the order of each multimapper unit's candidates comes from an `unordered_map` and so depends on history. The EM extractor re-sorts it, so only another reader of the global CSR would see it.

## Already recorded
- `ISSUES: hygiene-ledger`:
  - (h) scoring is written twice and disagrees on antisense at p = ½; `FragmentScorer` copies 13 fields into native.
  - (i) translation-unit constants such as `SCORED_STACK_CAPACITY`.
  - (e) the ZS gate checks set membership only.
- `ISSUES: latent-defects`:
  - the p = ½ antisense disagreement;
  - the nanobind `cast(...).data()` dangle and `StreamingScorer`'s missing keep_alive. The entry names `em_solver.cpp` and `solve_kernel.cpp`; `scoring.cpp`'s twelve `score_chunk` casts and its stored constructor pointers are the same pattern.
- `ISSUES: multimapper-intergenic-alignments`: the `n_cand > 0` gDNA condition and the gDNA term being a mean over hits.
- `ISSUES: scoring-penalties-are-underived-constants`: the 0.1 alphas. Its line "a zero alpha gives NaN today" points at a latent-defects item that is no longer listed; finding 1 is the live remainder.
- `ISSUES: binary-cuts-on-continuous-quantities`: `TAIL_DECAY`, `max_frag_length = 1000`, and pruning on a k·ln10 lattice.
- `ISSUES: the-minus-strand-coverage-weight-sat-one-fragment-length-off` (closed): covers the weight's strand handling, not finding 4.
- `ISSUES: the-em-prior-share-excludes-drawn-multimappers`.
- `TRAPS: nrna-does-not-mean-synthetic`.
