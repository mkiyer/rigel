<!-- title: Spliced-fragment-length deposit — fix the mature-projection bug (handoff spec) -->
# Deposit spliced fragments by fragment length — handoff spec

**Status:** audited + scoped, NOT started. This is the precise spec to resume in a new session. The fix lives in
the **accumulator** (C++ `src/rigel/native/calibration/accumulator.cpp` + its byte-for-byte Python reference
`tests/native/_accumulator_reference.py`) and the **scanner** (`src/rigel/native/bam_scanner.cpp`), so it is best
done as its own piece of work — separate from the calibration-side precision fix already shipped (main@f95b953d).

## 0. TL;DR — what's wrong and what to do

The spliced-derived mature estimate over-attributes **~63% of mature RNA to gDNA** because spliced mass is
*deposited* at **covered/read length** scale but *normalized* by a **fragment-length** eff-len. The two scales must
agree. **The fix: deposit spliced fragments using their transcript FRAGMENT length** (genomic span minus the
spliced-out introns), so the spliced mass lands on the same fragment-length scale as the contained/unspliced mass
and the existing `boundary_side_eff_length` recovers the density correctly. The hard part: a spliced alignment's
transcript fragment length is **often ambiguous** (the inner mate gap can cross junctions compatible with multiple
isoforms → multiple compatible lengths).

## 1. The bug

The unstranded gDNA estimate (`bp_solver.fit_enrichment_transfer`, the spliced-derived `ρ_g` blend, weight
`w_spl = 1−(2κ−1)²`, so it dominates as κ→½) subtracts mature RNA from the unspliced mass:

```
ρ_mature   = M_spliced(B) / E_rna_crossing(B)            # boundary spliced density
M_mature(R)= ρ_mature · E_rna_contained(R)               # projected contained mature
ρ_g(R)     = clip(M_unspliced(R) − M_mature(R), 0) / E_gdna(R)
```

`E_rna_crossing(B)` is `effective_length.boundary_side_eff_length` = `E[min(ℓ, L)]` over the RNA **fragment-length**
pmf (≈ `fl_mean` ≈ 257 in the toy). But `M_spliced(B)` is deposited at **covered length** scale (see §3), whose
per-side support is ≈ 110 (≈ read_length). Dividing a covered-scale mass by a fragment-scale eff-len under-states
`ρ_mature` by **~2.3×**, so only ~37% of the mature is subtracted; the other ~63% reads as gDNA.

## 2. Evidence (reproducible)

Diagnostic: `scripts/debug/mature_projection_check.py` (uses `dissect_regions.build_or_load_cache`; nrna_none ⇒
`payload_rna` = mature only, so "true contained mature" is exact). On `gdna_none_ss_0.50_nrna_none_capture_on`:

```
Σ TRUE contained mature (oracle)        = 644,358
Σ PROJECTED mature  ρ_mature·E_rna_cont  = 238,290   (projected/true = 0.370)
Σ residual (→ phantom gDNA)              = 409,148
E_rna_crossing USED  ≈ 257   |   EMPIRICAL spliced support S_emp = M_spliced/ρ_true ≈ 110   (ratio 0.43)
```

Stable across capture/no-capture and stranded/unstranded (0.35–0.37) ⇒ a pure normalization-scale bug, **not**
capture- or strand-specific. **Acceptance target for the fix: `projected/true → 1.0`.**

Post-EM impact (`scripts/debug/precision_benchmark_report.py`, suite `~/Downloads/rigel_runs/quick_3to1_5mb`):
this is the **non-nascent half of the unstranded residual** — the `none/ss0.5/none/CAP` +152K (prior-level)
phantom, and a contributor to the gDNA-present unstranded leak. **Masked in stranded data** (`w_spl≈0`), so it only
moves the unstranded conditions. The nascent half of the unstranded residual is a *separate*, FL-limited problem
(gDNA-vs-nascent) — out of scope here.

## 3. Root cause — the deposit geometry

`_accumulator_reference.py::deposit()` (the spec the C++ mirrors byte-for-byte):

- A fragment is a list of aligned `blocks` (genomic). A **spliced** fragment has ≥2 blocks with the intron gap(s)
  between them (the gap is *not* in `blocks`).
- Blocks are expanded into per-region `slices`; `L = Σ(slice_end − slice_start)` = the **covered (aligned)
  length** — for a spliced fragment this is ~the read coverage (≈ 2×read_length for paired), **excluding the
  intron**.
- Crossing mass is distributed `share = (slice_len)·(1/L)/n_cross` per boundary side.

So for a spliced fragment the per-side mass ≈ `exon_block_len / covered_L` — a **covered/read-scale** quantity. An
*unspliced* crossing fragment has `covered_L = genomic_span = fragment_length ℓ`, so its deposit *is* fragment-scale
and `boundary_side_eff_length` (a fragment-FL eff-len) recovers it correctly. The spliced fragment breaks this:
`covered_L ≠ ℓ` (the intron is missing), so its support is read-scale while the eff-len is fragment-scale.

`effective_length.py` has only `region_eff_length` (`E[max(0,L−ℓ)]`, contained) and `boundary_side_eff_length`
(`E[min(ℓ,L)]`, crossing) — both over the **fragment-length** pmf. There is **no spliced/covered-length eff-len**,
and the spliced channel is (correctly, for FL pooling) excluded from the gDNA FL pool (`fl_on = not spliced`).

## 4. The fix — deposit spliced by transcript fragment length

Put the spliced mass on the **fragment-length scale** so the existing fragment-FL eff-len matches. Conceptually: a
spliced fragment is *contiguous in transcript space* with length `ℓ_tx`; it should deposit across the splice
junction exactly like an unspliced crossing fragment of length `ℓ_tx` does across a region boundary — per-side
support `E[min(ℓ_tx, exon_side)]`. Then `boundary_side_eff_length` (unchanged) recovers `ρ_mature` correctly and
`projected/true → 1`.

`ℓ_tx` (transcript fragment length) = `genomic_span − Σ(spliced-out intron lengths between the mates)`.

## 5. The challenge — ambiguous transcript fragment length

- **Observed introns** (a read's CIGAR has `N`): the intron length is known; subtract it. No ambiguity for that
  read's own splice.
- **Inner-mate-gap introns** (paired-end; read1 ends in exon A, read2 starts in exon B, the junction(s) fall in the
  *unsequenced* inner gap): which intron(s) the gap spans is **inferred from the annotation**, and if several
  isoforms place different introns in the gap, `ℓ_tx` is a **set** of compatible lengths, not a point. This is the
  "fragment length of a spliced alignment is not known" case.

**Handling options (decide in the new session):**
1. **Compatible-set + FL model.** Enumerate the introns compatible with the mate placement (from the index's
   junction set), form the candidate `ℓ_tx` set, and weight by the RNA FL pmf (the maximum-likelihood / posterior
   `ℓ_tx`). Most principled; needs the scanner to reach the annotation junctions at deposit time.
2. **Observed-only.** Use only CIGAR-`N` introns; for inner-gap junctions, fall back (e.g. assume the shortest
   compatible intron, or the annotated primary). Simpler, biased for gap-spanning fragments.
3. **Per-fragment expected `ℓ_tx`.** Deposit a fractional mass spread over the compatible `ℓ_tx` set (mass-conserving),
   so the eff-len integrates the ambiguity. Cleanest statistically, most work.

The accumulator must stay **byte-for-byte** between C++ and the Python reference, and **mass-conserving** (each
fragment deposits total mass 1) — see `accumulator_mass_conservation` / `docs/accumulator/00_design.md`.

## 6. Implementation plan (file by file)

1. **Scanner — `src/rigel/native/bam_scanner.cpp`** (and any fragment-resolution helper): for a spliced fragment,
   compute `ℓ_tx` (genomic span − CIGAR `N` introns − inner-gap introns via the index junctions, per the §5 choice).
   Surface `ℓ_tx` (point, or the compatible set/weights) to the deposit call. This is where the index's junction
   annotation must be reachable.
2. **Accumulator — `accumulator.cpp` + `tests/native/_accumulator_reference.py::deposit()`:** for spliced
   fragments, deposit on the **transcript** scale: treat as a contiguous `ℓ_tx` fragment crossing the junction, so
   the per-side mass integrates to `ρ_mature·E[min(ℓ_tx, exon_side)]`. Keep the C++ and reference identical; keep
   mass conservation. (Unspliced path unchanged.)
3. **`effective_length.py`:** likely **no change** — the goal is for the *deposit* to be fragment-scale so the
   existing `boundary_side_eff_length` is correct. (If a separate spliced eff-len is cleaner, add it here and own
   it next to the others.)
4. **`bp_solver.fit_enrichment_transfer`:** **no change** if the deposit fix lands — `ρ_mature = M_spliced/E_rna_crossing`
   becomes correct automatically. (Double-check `build_node_geometry`'s `spliced_*`/`eff_rna_*` wiring still lines up.)
5. **Tests:** update `tests/native/` accumulator byte-for-byte tests for the new spliced deposit; add a
   uniform-mature **factor-1 invariant** test (analogous to the gDNA factor-1-under-uniform) asserting projected
   contained mature == true under uniform density.

## 7. Validation

- **Unit:** the new accumulator factor-1-under-uniform-mature invariant.
- **Diagnostic:** `OMP_NUM_THREADS=1 python scripts/debug/mature_projection_check.py
  gdna_none_ss_0.50_nrna_none_capture_on …` → `projected/true → ~1.0` (from 0.37). (Re-build the dissect cache
  after the accumulator change: it's keyed in `/tmp/dissect_cache/`; delete to force a re-scan.)
- **Benchmark:** re-quant the suite and compare with `precision_benchmark_report.py <suite> before_precision`
  (and `… phase2`). Expect the unstranded zero-gDNA-no-nascent phantom (`none/ss0.5/none/*`) and the unstranded
  gDNA leak to drop; stranded should be ≈ unchanged (the spliced term is `w_spl≈0` there). Regenerate goldens.

## 8. Code map / references

- **Bug consumer:** `src/rigel/calibration/bp_solver.py::fit_enrichment_transfer` (the `rho_spliced` block); the
  geometry it reads: `build_node_geometry` (`spliced_pos/neg_left/right`, `eff_rna_left/right`).
- **Eff-len:** `src/rigel/calibration/effective_length.py` (`boundary_side_eff_length`, `region_eff_length`).
- **Deposit:** `tests/native/_accumulator_reference.py::deposit` (lines ~140–260; the spliced/crossing path) +
  `src/rigel/native/calibration/accumulator.cpp` (mirror).
- **Scanner / FL:** `src/rigel/native/bam_scanner.cpp`; `src/rigel/frag_length_model.py` (RNA FL = `SPLICED_ANNOT`,
  bins at fragment length — confirm what it should bin for the spliced covered vs transcript length).
- **Diagnostics:** `scripts/debug/mature_projection_check.py` (the projected/true validator),
  `scripts/debug/precision_benchmark_report.py` (before/after net-flow; arg2 = baseline suffix).
- **Related:** `docs/calibration/precision_state_design.md` (the shipped precision fix that exposed this),
  `docs/accumulator/00_design.md` (deposit/mass-conservation spec), memory `precision_state_message_fix` /
  `accumulator_mass_conservation`.
