# The per-junction strand table — one small extension that fixes the RNA strand model

**Owner-directed design, 2026-07-28.** A deferred TODO, now load-bearing: the RNA strand overdispersion is
not currently a measurement (`strand_overdispersion_design.md` §3). This is the fix, and it is a strict
addition to the strand model — **it does not touch the accumulator, the mass channels, or the solver.**

---

## 1. THE DEFECT IN ONE TABLE

The two halves of the same Beta-Binomial are trained on **different populations**:

| | source | population | result |
|---|---|---|---|
| **κ** (the mean) | `StrandModels.exonic_spliced`, via `get_is_strand_qualified()` | **`SPLICED_ANNOT` only**, unique-mapper, unambiguous strand, non-chimeric | ✅ clean — reproduces the deep-junction rate (LBX0190 κ = 0.0023 vs 0.0029–0.0038 measured at depth ≥ 100) |
| **`od_r`** (the dispersion) | the accumulator's boundary spliced channels | `SPLICED_ANNOT` **+ `SPLICED_UNANNOT` + `SPLICED_IMPLICIT`** | ⛔ raw `od_mom` = **10.7–79.9**, an impossibility for an intraclass correlation, clipped to the 0.2 ceiling on **4/4** real libraries |

**The strand model's own docstring already states the correct rule** — *"Annotated splice junctions **prove
RNA origin**, making this an uncontaminated measure of library strand specificity"* — and
`get_is_strand_qualified()` already enforces it. We simply **discard the junction identity**, pushing two
parallel `int8` vectors (`exonic_spliced_obs`, `exonic_spliced_truth`), which leaves no way to compute a
dispersion *across junctions*. That is why the dispersion had to be scavenged from the accumulator, where
the population is wrong.

### Why the two contaminants are not equally harmful (owner, and it is the key asymmetry)

* **`SPLICED_UNANNOT` — near-zero effect.** Its endpoints are, by definition, not annotated splice sites, so
  its spliced flux credits boundaries **elsewhere**; it dilutes other places rather than the
  annotated-junction seeds.
* ⛔ **`SPLICED_IMPLICIT` — the damaging one.** It is a paired-end fragment whose **mate gap spans an
  annotated intron**, so it credits **exactly the annotated boundaries where the genuine junctions live**.
  And **the splice motif was never sequenced**: the scanner orients it by the annotated intron's transcript
  strand (`motif_strand |= isj.strand`). The motif is *the* independent ground truth — with no motif the
  aligner-vs-motif comparison is not noisy, it is **undefined**. For a gDNA fragment whose mate gap merely
  straddles an intron, the orientation is arbitrary ⇒ ~50/50, which is exactly what was measured
  (MO_3021 at depth ≥ 10 reads **0.4876**).

✅ **`SPLICE_ARTIFACT` is already fully excluded** (`bam_scanner.cpp:1443` — no deposit, no FL pool), so the
aligner's *known* systematic errors are not in the training set. What is in it is the classes never filtered.

## 2. THE DESIGN

> A splice junction is uniquely specified by **(reference, start, end, genomic splice-motif strand)**.
> Hash strand-qualified spliced fragments on that key and track **sense** (aligner orientation agrees with
> the motif) and **antisense** (disagrees) counts per junction.

From one table we get both quantities, **from the same population**:

```
    κ      = Σ_j sense_j / Σ_j (sense_j + antisense_j)          the global strand specificity
    od_r   = Beta-Binomial dispersion of (sense_j | n_j) across junctions, mean κ
```

### 2.1 What changes

**C++ (`bam_scanner.cpp`) — a strict addition:**

```cpp
struct SJStrandKey { int32_t ref_id; int32_t start; int32_t end; int8_t motif_strand; };

struct StrandObservations {
    std::vector<int8_t> exonic_spliced_obs;      // unchanged
    std::vector<int8_t> exonic_spliced_truth;    // unchanged
    std::vector<int8_t> exonic_obs;              // unchanged
    std::vector<int8_t> exonic_truth;            // unchanged
    // NEW — per-junction sense/antisense, keyed as above
    std::unordered_map<SJStrandKey, std::pair<uint32_t, uint32_t>, SJStrandKeyHash> sj_strand;
};
```

* Populated **inside the existing `if (result.get_is_strand_qualified())` branch** — so the gate, and
  therefore the population, is *unchanged and already correct*. Nothing new is admitted.
* `merge_strand_obs` gains a map merge (sum both counters per key) — same per-worker-then-merge pattern
  the vectors already use.
* Exported as six parallel arrays (`sj_ref_id`, `sj_start`, `sj_end`, `sj_strand`, `sj_n_sense`,
  `sj_n_antisense`) via the existing `vec_to_ndarray` idiom into `result["strand_observations"]`.

**Python (`strand_model.py`, `gdna_strand.py`):**

* `StrandModels` gains the per-junction table.
* `fit_strand_balance` is **unchanged** (κ still from the aggregate) — see the gate in §4.
* `fit_rna_strand_overdispersion` is fed the **per-junction** `(sense_j, n_j)` instead of boundary sides,
  at mean κ, through the estimator landed in `c99c399f` — which already weights by **information**
  (`I = 1/Var(od_mom)|₀`, not the pair count, since μ = κ ≠ ½) and shrinks to `od₀`.

**`fit_rna_strand_from_substrate` is deleted.** Its boundary-side seeds are the contaminated population.

### 2.2 What does NOT change

* **The accumulator** — no new channel, no spec change, no `tests/native/` change, no C++ ABI change to
  `deposit`. The `spliced` predicate stays as it is, because the spliced **mass** channel legitimately
  wants implicit and novel RNA for the peel, the graft and the mature-RNA floor. **Narrowing that predicate
  would silently delete real RNA from the solver's mass accounting — a far larger blast radius, and it is
  why the fix belongs in the strand model instead.**
* The solver, the deconvolution, the region/boundary map.
* κ, and therefore every κ-dependent quantity.

**The only behavioural change is `od_r`**, which enters `strand_loglik` on every node. That is the point.

## 3. DESIGN DECISIONS THAT ARE STILL OPEN

**D1 — multi-junction fragments: credit one junction or all?** A fragment crossing K junctions is one
observation today. Per-junction semantics naturally credit all K.
* ⚠ Crediting all K introduces **within-fragment correlation** — the same molecule, hence the same
  orientation, counted K times — which *inflates* the apparent dispersion. This is the same hazard the
  current RNA fit already documents as its "v1 approximation" (a fragment spanning K junctions credits ~K
  boundary sides).
* ⚠ It also changes κ if κ is recomputed from the table (multi-junction fragments gain weight).
* **Recommendation: credit all K in the table, but keep κ from the existing aggregate vectors**, so κ is
  bit-identical and the change is provably isolated to `od_r`. Then *measure* how far the table-derived κ
  differs; if it is negligible, collapse to one source in a follow-up.

**D2 — is the motif strand redundant in the key?** `(ref, start, end)` determines the motif, so the strand
is derivable. Keeping it in the key is harmless, self-documenting, and catches the pathological case of a
position serving junctions of both strands. **Recommendation: keep it, and assert consistency** — a key
collision with a differing strand is a real signal worth surfacing rather than silently merging.

**D3 — depth restriction.** At κ ≈ 0.99 (or ≈ 0.01) a junction needs many reads before a disagreeing read
is even possible. ⭐ **Probably needs no rule at all**: the information weighting landed in `c99c399f`
already gives a 1-read junction ~zero weight automatically. **Recommendation: no threshold; verify the
weighting suffices, and only revisit if it does not.**

**D4 — memory.** ~200–300 k annotated junctions in a human transcriptome × (2 × uint32 + key) is a few MB
per worker before merge. Not a concern; state it so nobody re-litigates it.

## 4. VALIDATION PLAN

**Gate 1 — κ is BIT-IDENTICAL.** The population and the aggregate are untouched; if κ moves at all, the
change is not isolated and something is wrong. This is the falsification test.

**Gate 2 — `od_r` on real data.** Pre-registered: it must **fall off the ceiling** on all four libraries and
land near the honest measurement, **0.001–0.016** (`strand_overdispersion_design.md` §3), i.e. `Beta(a,a)`
with a ≈ 31–464. ⭐ **This is the strongest available check, because that number was measured independently
— from deep junctions, before this design existed.** If the table reproduces it, two independent routes
agree.

**Gate 3 — the synthetic suite.** `od_r` is currently 0.0008–0.0017 there and the suite has true `od = 0` by
construction, so it must stay small. ⚠ It cannot validate a *non-zero* `od_r`; only real data can.

**Gate 4 — the 32-condition A/B, both refit settings.** `od_r` enters `strand_loglik` on every node, so this
is not a local change. **Pre-register stranded conditions as the falsifier**: they carry real strand
information, so a correct de-contamination should help or be neutral there. A large stranded regression
means we removed signal rather than contamination.

**Gate 5 — goldens LAST**, and re-record every A/B baseline from the current tree in the same session.

## 5. WHY THIS IS THE RIGHT SHAPE

* It **adds** a table rather than changing a rule — the population gate is already correct and is reused
  verbatim, so there is no new admission policy to get wrong.
* It puts κ and `od_r` on **one population**, which is the actual defect.
* It gives the **per-junction depth structure** the owner identified as necessary, for free.
* It leaves the accumulator and the solver alone, so the blast radius is one scalar.
* It removes the last consumer of the contaminated boundary-side spliced seeds.

**Not implemented.** Requires a C++ rebuild (`pip install --no-build-isolation -e .`).
