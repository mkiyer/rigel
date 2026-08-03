# SPEC — search every fragment's unsequenced gaps for introns, not just unspliced ones

    Status: ✅ **IMPLEMENTED 2026-08-01.** What LANDED, what it MEASURED and what it did NOT close are
            in `docs/WIP.md`'s **C2.6** entry, which supersedes this file wherever they differ.
            This document is kept as the record of the DESIGN, not of the outcome.
    Cause and evidence: `docs/accumulator/JUNCTION_OPPORTUNITY.md` §4.
    Owner ruling, 2026-08-01: *"we should be searching for gap introns within every fragment."*
    Blocked: C3 (the RNA opportunity work) — ⭐ now UNBLOCKED.

    ⭐ THE THREE THINGS THE MEASUREMENT CHANGED, so this file is not read as the final word:
      1. G-tail did NOT reach 0. It went 0.00909 -> 0.00137 above the true 713 bp ceiling and
         `dropped_too_long` 280,558 -> 38,309. ⭐ The residual is **D3**, and that is measured, not
         assumed: emitting EVERY intron in a gap takes it to 0.00002 / 389. §3's D3 was right to be
         deferred and is now the ONLY known mechanism left.
      2. ⛔ **D1 costs more than it buys.** Cutting the intron takes `RNA_SPLICED` from +8.00 % to
         +0.67 % against truth; then REMOVING the mixed fragments takes it to −9.58 %, because they are
         the long ones. §3's D1 asked for the two effects to be reported separately — they were, and the
         answer inverts the recommendation. **Owner decision, one line to reverse.**
      3. ⛔ **X5's prediction was wrong.** Leaving the union unsorted fails NOTHING, because
         `normalise_introns` sorts internally in both the C++ and the specification. The sort is
         defensive, not load-bearing. §4's perturbation table is corrected in the LEDGER.

---

## §0 THE BUG IN ONE PARAGRAPH

`resolve_context.h:1442` runs implicit-splice detection **only when the resolver has already classified
the fragment `SPLICE_UNSPLICED`**. A paired-end fragment whose reads carry an observed CIGAR-N splice
therefore never has its **unsequenced mate gap** examined. If an annotated intron lies inside that gap it
is never cut, so the accumulator's `L` — the molecule length, and since C2 the tool's *one* definition of
fragment length — over-counts by the whole intron.

⚠ **`UNSPLICED` is not "single aligned block".** An unspliced PE fragment already has two blocks and a
mate gap; that is the case the detector was built for. The gate is on the **splice type**, not on block
count. The population being missed is *spliced fragments that also have a gap intron* — necessarily
long, because they span two or more introns.

### The evidence, restated (`docs/accumulator/JUNCTION_OPPORTUNITY.md` §4 has the derivation)

The pilot's simulated mRNA support is **50–713 bp**. Measured on the shipped tally:

| series | ≥600 bp | ≥700 | ≥800 |
|---|---|---|---|
| `DNA_INTERGENIC` pool — gDNA, **no introns to miss** | 0.00023 | 0.00001 | 0.00000 |
| ⭐ **truth, gdna** | **0.00024** | **0.00001** | **0.00000** |
| `RNA_SPLICED` pool — RNA, **always has an observed splice** | **0.04182** | **0.03002** | **0.01853** |
| the anchor `deposited_lengths` — mostly unspliced, so mostly *does* get detection | 0.00663 | 0.00467 | 0.00288 |

⭐ **`L` is exact to five decimals when the fragment has no introns and badly wrong when it has them**,
and the pool where the gate *guarantees* detection is skipped is **6× worse than the anchor**. A further
**6.1 % of deposits** (`qc.dropped_too_long = 280,558`) are discarded as longer than 1000 bp on a library
whose longest molecule is 713 bp. Truncating the anchor at 600 bp removes **95.5 %** of its sd excess.

⚠ **C0's proof is untouched.** `L` is correct *given its inputs*; the **intron list is incomplete**.

---

## §1 ⛔ THE TRAP: A CIGAR-N INTRON IS ALSO A "GAP"

`collect_implicit_splice_introns` builds its gap list by walking **consecutive aligned blocks** and
emitting every hole between them (`resolve_context.h:858-871`). A read `50M200N50M` produces **two
blocks**, so the 200 bp `N` is a hole exactly like a mate gap.

⛔ **Deleting the `splice_type == SPLICE_UNSPLICED` gate and nothing else is wrong.** The detector would
then re-derive introns the CIGAR already stated, and:

* an **exact** re-derivation is harmless — `_normalise_introns` de-duplicates;
* a **near** match is not. `transcript_has_implicit_intron_in_gap` accepts any annotated intron inside
  the gap **± K** (`splicing_anchor_tolerance_`), so a *different* nearby intron can match. Two
  overlapping introns then normalise into one wider interval and `L` comes out **too short** — the same
  class of defect in the other direction.

⭐ **The rule: only search gaps that no observed intron already explains.** A hole between blocks is
either a CIGAR-N intron — in which case it equals an entry in `f.introns` exactly, both being derived
from the same CIGAR — or it is unsequenced. Only the latter is a candidate.

---

## §2 THE CHANGE

### 2.1 `resolve_context.h` — `collect_implicit_splice_introns`

Add an input: the fragment's **observed** introns. Skip any gap that exactly matches one.

```
bool collect_implicit_splice_introns(
        const std::vector<ExonBlock>&   exons,
        const std::vector<IntronBlock>& observed_introns,   // NEW
        const std::vector<int32_t>&     candidate_t,
        ResolverScratch&                scratch,
        std::vector<IntronBlock>&       out_introns,
        bool&                           out_ambiguous) const
```

* after the gap list is built, drop every `gap` for which some `observed_introns[i]` has
  `start == gap.start && end == gap.end` on the same `ref_id`;
* ⚠ **exact equality, not overlap.** Overlap would also drop a genuine gap intron that happens to abut
  an observed one, and that is the population this work exists to find.
* everything downstream of the gap list — the unanimity scan, the synthetic-row filter, `out_ambiguous`
  — is **unchanged**.

### 2.2 `resolve_context.h:1442` — the gate

```
-       if (cr.splice_type == SPLICE_UNSPLICED &&
-           cr.chimera_type == CHIMERA_NONE &&
-           collect_implicit_splice_introns(exons, cr.t_inds, scratch,
-                                           cr.implicit_introns, cr.implicit_ambiguous)) {
-           cr.splice_type = SPLICE_IMPLICIT;
-       }
+       if (cr.chimera_type == CHIMERA_NONE &&
+           collect_implicit_splice_introns(exons, introns, cr.t_inds, scratch,
+                                           cr.implicit_introns, cr.implicit_ambiguous) &&
+           cr.splice_type == SPLICE_UNSPLICED) {
+           cr.splice_type = SPLICE_IMPLICIT;      // ⚠ the PROMOTION is still unspliced-only
+       }
```

⭐ **Detection is unconditional; the `SPLICE_IMPLICIT` *promotion* stays unspliced-only.** `splice_type`
is consumed by `scoring.cpp`, the buffer, the strand training and — since C2 — `rigel report`'s census.
Re-labelling an observed `SPLICED_ANNOT` fragment as `SPLICED_IMPLICIT` would silently move mass between
reported categories and change what the strand model trains on. **This work is about `L`, not about
classification.** Keep them separate.

⚠ `chimera_type == CHIMERA_NONE` is **kept**: a chimeric fragment's blocks are not one molecule and its
gaps are not introns.

### 2.3 `bam_scanner.cpp` — the deposit adapter, ~line 1506

The two intron lists stop being mutually exclusive, so **union** them:

```
-       const bool implicit = (st == SPLICE_IMPLICIT);
-       // Mutually exclusive — implicit is only detected when no CIGAR-N exists.
-       const std::vector<IntronBlock>& cut_introns = implicit ? cr.implicit_introns : f.introns;
+       // ⭐ A fragment may have BOTH: an observed CIGAR-N splice AND an annotated intron inside its
+       // unsequenced mate gap. Both are cut from L. They are disjoint by construction — §1.
+       (concatenate f.introns and cr.implicit_introns into ws.deposit_introns)
```

⚠ The existing per-`ref_id` filter and the adjacent-duplicate skip must apply to the **concatenated**
list, and it must be **sorted by (start, end)** before the adjacent-duplicate skip, which today relies on
`intron_set` already being ordered. ⛔ Do not re-key or merge `intron_set` upstream — the comment at
`bam_scanner.cpp` records that merging demotes `SPLICED_ANNOT` to `SPLICED_UNANNOT` and widens `t_inds`.

Then the three `FragmentPath` fields that assume exclusivity:

| field | today | after |
|---|---|---|
| `path.sj_implicit` | `st == SPLICE_IMPLICIT` | ⭐ **`!cr.implicit_introns.empty()`** — "this `L` depends on an intron that was never sequenced", which is what the flag is *for* |
| `path.path_ambiguous` | `implicit && cr.implicit_ambiguous` | ⭐ **`cr.implicit_ambiguous`** — if any gap's candidates disagree, `L` is undetermined regardless of whether the fragment also has an observed splice |
| `path.sj_strand` | replaced by the OR of implied strands **iff** implicit | ⚠ **decision D2 below** |

---

## §3 ⛔ THREE DECISIONS THE IMPLEMENTER MUST NOT MAKE SILENTLY

### D1 — do mixed fragments leave the `RNA_SPLICED` pool? ⭐ **Recommend: yes**

`_accumulator_reference._pool()`: `if spliced: return None if sj_implicit else FragmentPool.RNA_SPLICED`.
Setting `sj_implicit` per §2.3 therefore **removes mixed fragments from the pure RNA length pool**.

That is the recommendation, and the existing rationale carries: the pool is a **length** histogram used
to *fit* the fragment-length model, and a length partly inferred from the annotation is a product of the
model it would be used to fit. ⚠ The counter-argument is real and should be recorded rather than
dismissed — RNA *certification* comes from the observed splice and is not in doubt; only the length is.

⚠ **Measure the mass this moves before and after.** Both effects push the same way (the intron is cut, so
`L` shrinks; and the fragment leaves the pool), so the pool's mean will drop for two reasons at once and
they must be reported separately.

### D2 — what is `sj_strand` for a mixed fragment?

Today: an implicit fragment's `sj_strand` is the OR of the implied introns' transcript strands, because
no motif was sequenced. A mixed fragment **has** a sequenced motif.

⭐ **Recommend: keep the observed `cr.sj_strand` and do not OR in the implied strands.** An observed
motif is evidence; an implied strand is an inference from the annotation, and `sj_strand`'s stated job
(`_sj_edge_id`) is to resolve an observed intron against the annotation. ⚠ Mixing an inference into an
observation is how `primary` went wrong.

### D3 — a gap containing MORE THAN ONE intron is still not handled

`transcript_has_implicit_intron_in_gap` returns the **first** matching intron and stops
(`resolve_context.h:790-800`). A mate gap spanning two annotated introns of the same transcript therefore
has only one cut, and `L` remains too long by the second.

⭐ **Recommend: out of scope here, but MEASURE it and record the residual.** Fixing the gate is the large
effect; this is a second-order one and conflating them makes neither measurable. ⚠ It must be measured,
not assumed small — human introns are usually far longer than a mate gap, but short-intron species and
short-intron genes exist.

---

## §4 THE GATES — written first, verified failing

⚠ House rule: write each gate, watch it fail, then implement, then **perturb the code and watch the gate
fail again**. `falsification_needs_perturbation`.

| | gate | where |
|---|---|---|
| **U1** | ⭐ **A spliced fragment with an annotated intron in its mate gap has that intron cut.** Construct in the Python reference: blocks `[100,150) [400,450)` with an observed CIGAR-N `[150,400)` and a candidate transcript whose next intron is `[500,900)` inside a mate gap `[450,1000)`. `L` must exclude **both** | `tests/test_implicit_splice.py` |
| **U2** | **An observed CIGAR-N gap is NOT re-derived.** Same fragment, and the candidate's intron list also contains the observed intron: `cr.implicit_introns` must not contain it, and `L` is unchanged from U1 | `tests/test_implicit_splice.py` |
| **U3** | ⛔ **The near-match trap of §1.** A candidate whose intron is inside the observed gap ± K but *not equal* to it must **not** be emitted for that gap, and `L` must not shrink | `tests/test_implicit_splice.py` |
| **U4** | **`splice_type` does not move.** A `SPLICED_ANNOT` fragment that gains a gap intron is still `SPLICED_ANNOT` — the census, the strand training and the report do not shift | `tests/test_resolution.py` |
| **U5** | **`path_ambiguous` fires on a mixed fragment** whose gap candidates disagree, and the deposit is rejected and counted | `tests/native/test_implicit_splice_deposit.py` |
| **S3** | **Byte-identity**: the C++ matches `_accumulator_reference.py` exactly, and is bit-identical at 1/2/4/8 workers | `tests/native/test_accumulator_native_parity.py` (automatic) |
| ⭐ **G-tail** | **THE HEADLINE.** On `gdna_none ss0.99 capture_off` the anchor's support ceiling must fall to the library's true maximum — **713 bp**, read from `truth_fragment_lengths.tsv`, not chosen — and `qc.dropped_too_long` must collapse from **280,558**. ⚠ Score the mass ≥700 bp: truth is **0** | new, `scripts/design/fl_anchor_gap.py` extended |
| ⭐ **G-sd** | The anchor's sd against truth must fall from **+27.0 %**. This is `docs/accumulator/FRAGMENT_LENGTH_AUDIT.md` G3b, and it is **this work's** gate, not C3's | same |
| **G-gdna** | ⛔ **The control must not move.** `DNA_INTERGENIC` vs truth gdna stays exact (0.00023 vs 0.00024 at ≥600) — a change there means the fix reached fragments with no introns, which is impossible and therefore a bug | same |

### The perturbations, chosen to break each distinct claim

| | perturbation | must fail |
|---|---|---|
| **X1** | restore the `splice_type == SPLICE_UNSPLICED` gate | G-tail, U1 |
| **X2** | drop the §2.1 observed-gap filter | U2/U3 — and ⚠ if nothing fails, the filter is untested and the near-match trap is live |
| **X3** | union the intron lists but leave `path.sj_implicit = (st == SPLICE_IMPLICIT)` | D1's pool-purity gate |
| **X4** | leave `path_ambiguous = implicit && cr.implicit_ambiguous` | U5 |
| **X5** | concatenate the intron lists **without sorting** before the adjacent-duplicate skip | S3 byte-identity |

---

## §5 BLAST RADIUS

| | |
|---|---|
| **`payload_schema_digest`** | ⭐ **does not move** — no new field. The 8 pilot caches stay *loadable*, ⛔ **but their contents are now wrong** and must be rebuilt (~56 s) because the tally itself changes |
| **S3 byte-identity** | **reopens.** `_accumulator_reference.py` is the specification and moves first; the C++ follows |
| **goldens** | move a **third** time (P1 units → C2 FL models → this). ⛔ Still regenerate **once**, after C3, twice and diff |
| **`rigel report`** | splice census **unchanged** by design (§2.2). ⚠ The FL histograms move — that is the point |
| **every FL consumer** | the scorer, the calibration divisors, and by **D7** every transcript's effective length in the EM. This is upstream of all of them |
| **`docs/accumulator/JUNCTION_OPPORTUNITY.md` §3** | ⭐ **re-run the θ control afterwards.** Its +59.8 % sd figure is measured against a contaminated anchor; C3's true target is only knowable once this lands |

⚠ **`docs/SESSION_HANDOFF.md` §1 and the LEDGER's C1 entry both need a correction note**: C1 recorded the
support ceiling moving 713 → 1000 as a win. **713 was the library's true maximum.** The scanner's deleted
definition **B** took its length from the transcript, so it could not produce an uncut intron — it was
right about the ceiling and `L` was wrong. That does not undo C2, whose reasons all still hold.

---

## §6 HOW TO KNOW IT WORKED

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
pip install --no-build-isolation -e ".[dev]"
python -m pytest tests/ -q                        # ⛔ python -m, not bare pytest
OMP_NUM_THREADS=1 python scripts/design/build_scan_cache.py --force
OMP_NUM_THREADS=1 python scripts/design/fl_anchor_gap.py
```

**Done when:** the anchor's mass ≥ 700 bp is **0** on the zero-gDNA conditions (truth: 0),
`qc.dropped_too_long` has collapsed from 6.1 % of deposits, the gDNA control is unmoved, and the suite is
green apart from the goldens.

⛔ **Do not tune anything to reach those numbers.** Every target above is read from
`truth_fragment_lengths.tsv` or from a control that must not move. If the residual will not close, the
remaining mass is D3's multi-intron gaps or a third mechanism — **measure it and report it**, do not
close the gap with a constant (`docs/SESSION_HANDOFF.md` §3 trap 12).
