# The per-junction SJ strand table — design and implementation plan

**Owner-directed, 2026-07-28.** A deferred TODO, now load-bearing: the RNA strand overdispersion is not
currently a measurement. This is the fix. It is a **strict extension of the strand model** and it touches
neither the accumulator, the mass channels, nor the solver.

---

## 1. THE DEFECT

**The two halves of the same Beta-Binomial are trained on different populations.**

| | source | population | result |
|---|---|---|---|
| **κ** (the mean) | `StrandModels.exonic_spliced`, gated by `get_is_strand_qualified()` | `SPLICED_ANNOT` only, unique-mapper, unambiguous strand, non-chimeric | ✅ **clean** — reproduces the deep-junction rate (LBX0190 κ = 0.0023 vs 0.0029–0.0038 measured at depth ≥ 100) |
| **`od_r`** (the dispersion) | the accumulator's boundary spliced channels | `SPLICED_ANNOT` **+ `UNANNOT` + `IMPLICIT`** | ⛔ raw `od_mom` = **10.7–79.9** — impossible for an intraclass correlation — clipped to the 0.2 ceiling on **4/4** real libraries |

The strand model's own docstring already states the correct rule — *"Annotated splice junctions **prove RNA
origin**, making this an uncontaminated measure of library strand specificity"* — and
`resolve_context.h:97` already enforces it:

```cpp
bool get_is_strand_qualified() const {
    return splice_type == SPLICE_SPLICED_ANNOT      // annotated ONLY
        && !ambig_strand
        && (align_strand == STRAND_POS || align_strand == STRAND_NEG)
        && (sj_strand == STRAND_POS || sj_strand == STRAND_NEG)
        && chimera_type == CHIMERA_NONE;
}
```

**We simply discard the junction identity.** `bam_scanner.cpp:1687` pushes two parallel `int8` vectors of
strand labels, so there is no way to compute a dispersion *across* junctions — which is why it was
scavenged from the accumulator, where the population is wrong.

### 1.1 Why the two contaminants are not equally harmful (owner)

* **`SPLICED_UNANNOT` — ~zero effect.** Its endpoints are by definition not annotated splice sites, so its
  spliced flux credits boundaries **elsewhere**. It dilutes other places, not the annotated-junction seeds.
* ⛔ **`SPLICED_IMPLICIT` — the damaging one.** Its mate gap spans an **annotated intron**, so it credits
  **exactly the annotated boundaries where the genuine junctions live**, and **the motif was never
  sequenced** (`motif_strand |= isj.strand` uses the annotation instead). The motif is *the* independent
  ground truth; with no motif the aligner-vs-motif comparison is not noisy, it is **undefined**. For a gDNA
  fragment whose gap merely straddles an intron the orientation is arbitrary ⇒ ~50/50, exactly as measured
  (MO_3021 at depth ≥ 10 reads **0.4876**).
* ✅ **`SPLICE_ARTIFACT` is already fully excluded** (`bam_scanner.cpp:1443` — no deposit, no FL pool).

## 2. THE DESIGN

> A splice junction is uniquely specified by **(reference, start, end, genomic splice-motif strand)**.
> Hash strand-qualified spliced fragments on that key; track **sense** (aligner orientation agrees with the
> motif) and **antisense** counts per junction.

```
    κ     = Σⱼ senseⱼ / Σⱼ (senseⱼ + antisenseⱼ)                the global strand specificity
    od_r  = Beta-Binomial dispersion of (senseⱼ | nⱼ) across junctions, at mean κ
```

### 2.1 ⭐ The existing 2×2 is EXACTLY the marginal of this table

`StrandModel` is a 2×2 contingency table of `align_strand × sj_strand`. Writing `sense ≡ align == sj`:

| recovered from the table | equals |
|---|---|
| `Σ` over motif-**POS** junctions of **sense** | `pos_pos` |
| `Σ` over motif-**POS** junctions of **antisense** | `neg_pos` |
| `Σ` over motif-**NEG** junctions of **sense** | `neg_neg` |
| `Σ` over motif-**NEG** junctions of **antisense** | `pos_neg` |

**The per-junction table is a strict REFINEMENT: the current 2×2 is its marginal, exactly** — provided each
fragment credits exactly one junction (§2.3). So the 2×2 becomes a **derived view** and there is one source
of truth, which is the cleanup the owner asked for.

### 2.2 ⚠ Integer widths — the owner's challenge, answered

> *"You said int8 which doesn't seem right.. we can easily have more than 256 counts?!?"*

**The existing `int8` is correct, because it does not store counts.** `std::vector<int8_t>
exonic_spliced_obs` holds one **strand label** per fragment (`STRAND_POS = 1`, `STRAND_NEG = 2`); the
*count* is the vector's `size()`. So there is no overflow today. But the challenge is right in substance —
**the new per-junction fields ARE counts and must be sized as counts**:

| field | type | why |
|---|---|---|
| `ref_id` | `int32` | matches the existing reference id |
| `start`, `end` | `int64` | genomic coordinates; `int32` would suffice for human (chr1 = 248 Mb) but not for every genome, and the cost is nil |
| `motif_strand` | `int8` | a **label** (POS/NEG), like the existing vectors |
| **`n_sense`, `n_antisense`** | **`uint64`** | these are **counts**. `uint32` would be safe in practice (the deepest plausible single junction is ~10⁷ reads, 400× headroom) but at ~200–300 k human junctions `uint64` costs ~2.4 MB extra and removes the question permanently. **Take `uint64`.** |
| global totals in Python | `int64` | numpy default; a 10⁹-read library sums well past `uint32` |

### 2.3 ⭐ The multi-junction question — RESOLVED by how `sj_strand` is obtained

A fragment can cross K junctions. The decisive fact: **`sj_strand` is read from the BAM `XS`/`ts` TAG**
(`bam_scanner.cpp:600`, `read_sj_strand`) — it is **one value per FRAGMENT**, not per junction. So all K
junctions of a fragment necessarily share the **same** sense bit.

**Therefore crediting all K would duplicate a single observation K times: zero added information, maximal
within-fragment correlation, and a directly inflated dispersion — the very quantity we are trying to
measure honestly.**

**⇒ CREDIT EXACTLY ONE JUNCTION PER FRAGMENT — the leftmost (first) CIGAR-N junction.** Deterministic; one
fragment = one observation; the 2×2 stays exactly recoverable (§2.1); no correlation. A fragment's
orientation is a single draw regardless of how many junctions it spans, so nothing is lost.

⚠ Record in the code that this is a *choice*: attributing a multi-junction fragment to its first junction
is arbitrary but unbiased. If it ever matters, the alternative is a 1/K split — **do not** credit all K.

### 2.4 What changes, and what does not

**Changes:**
* `bam_scanner.cpp` — a per-junction map in `StrandObservations`, populated inside the **existing**
  `if (result.get_is_strand_qualified())` branch, merged per worker, exported as six parallel arrays.
* `strand_model.py` — `StrandModel` gains the table; the 2×2 becomes derived (§2.1); `StrandModels` exposes it.
* `gdna_strand.py` — `fit_rna_strand_overdispersion` is fed per-junction `(senseⱼ, nⱼ)` at mean κ, through
  the estimator landed in `c99c399f` (which already weights by **information**, not the pair count, since
  μ = κ ≠ ½).
* **`fit_rna_strand_from_substrate` is DELETED** — its boundary-side seeds are the contaminated population.

**Does NOT change:**
* ⭐ **The accumulator.** No new channel, no spec change, no `tests/native/` change, no `deposit` ABI
  change. The `spliced` predicate stays as-is because the spliced **mass** channel legitimately wants
  implicit and novel RNA for the peel, the graft and the mature-RNA floor. **Narrowing it would silently
  delete real RNA from the solver's mass accounting — a far larger blast radius, and the reason this fix
  belongs in the strand model.**
* The solver, the deconvolution, the region/boundary map, κ, and every κ-dependent quantity.

**The only behavioural change is `od_r`**, which enters `strand_loglik` on every node. That is the point.

## 3. THE STRAND-MODEL CLEANUP (owner-requested)

> *"I want the implementation to also clean up how we build and use the strand model. This can live inside
> the strand model, which should already be an elegantly designed data structure."*

1. ⭐ **One source of truth.** The 2×2 becomes a derived property of the junction table (§2.1), not four
   independently-maintained counters.
2. **Retire the observe/finalize lifecycle if it is now vestigial.** `observe()` / `observe_batch()` /
   `_cached_p_sense` / `_finalized` exist because the model was built incrementally. If it is now built
   once from arrays, a frozen constructor + computed properties is simpler and removes a mutable-state
   footgun. ⚠ Check `observe()`'s remaining callers first.
3. ⚠ **Rename `StrandBalance.rna_strand_overdispersion`.** It is `1/(n_obs+3)` — a **posterior width**, a
   QC-only thin-count power diagnostic — and it *collides by name* with the decode's genuine Beta-Binomial
   `rna_strand_overdispersion` from `gdna_strand.py`. CLAUDE.md already has to warn that these are
   "distinct". Two different quantities must not share a name. Suggest `posterior_width` or `kappa_se`.
4. **Keep the `exonic` diagnostic sub-model — it is LIVE**, not dead: `cli.py:493-497` reports
   `exonic_all_specificity` and the `exonic_spliced − exonic` gap as the **gDNA-contamination QC signal**.
5. **Surface the junction table in QC**: number of junctions, depth distribution, and the fitted `od_r`
   with its information — because "how many junctions are deep enough to measure this" is now a
   first-class question about a library.

## 4. IMPLEMENTATION PLAN (ordered, each step independently verifiable)

**S1 — C++ collection.** Add `SJStrandKey` + hash, the map in `StrandObservations`, population inside the
existing qualified branch (leftmost junction only), the merge, and the six-array export.
*Gate: rebuild (`pip install --no-build-isolation -e .`); κ **bit-identical**; the table's marginal 2×2
equals the existing 2×2 **exactly** on a real BAM and on the synthetic suite.* **This is the whole
correctness argument — if the marginal does not match, stop.**

**S2 — Python model.** Build the table into `StrandModel`; make the 2×2 derived; keep the public API
(`p_r1_sense`, `strand_specificity`, `n_observations`) unchanged.
*Gate: all existing strand tests pass untouched; κ bit-identical.*

**S3 — cleanup.** §3 items 2–5. *Gate: behaviour-neutral; goldens untouched.*

**S4 — the `od_r` refit.** Feed per-junction `(senseⱼ, nⱼ)` to `fit_rna_strand_overdispersion`; delete
`fit_rna_strand_from_substrate` and its now-dead boundary-side seed path.
*Gate: §5 — this is the only step that changes behaviour.*

**S5 — the 32-condition A/B, then goldens LAST.**

## 5. VALIDATION

**Gate 1 — κ is BIT-IDENTICAL** (S1, S2, S3). The population and the aggregate are untouched. **If κ moves
at all, the change is not isolated and something is wrong.** This is the falsification test.

**Gate 2 — ⭐ `od_r` on real data must land at 0.001–0.016.** That range was measured **independently, from
deep junctions, before this design existed** (`strand_overdispersion_design.md` §3): 0.0035 / 0.0020 /
0.0023 / 0.0011 on LBX0190 / LBX0588 / MO_3021 / vcap, i.e. `Beta(a,a)` with a ≈ 31–464. **If the table
reproduces those numbers, two unrelated routes agree — the strongest available check.** It must also fall
off the 0.2 ceiling on all four.

**Gate 3 — synthetic.** `od_r` is 0.0008–0.0017 today and the suite has true `od = 0` by construction, so
it must stay small. ⚠ It **cannot** validate a non-zero `od_r`; only real data can.

**Gate 4 — 32-condition A/B, both refit settings.** `od_r` enters `strand_loglik` on every node.
**Pre-register stranded conditions as the falsifier**: they carry real strand information, so a correct
de-contamination should help or be neutral there. A large stranded regression means signal was removed, not
contamination.

**Gate 5 — goldens LAST**, and re-record every A/B baseline from the current tree in the same session.

## 6. RISKS

* **The XS/ts tag may be absent or unreliable on some BAMs.** `sj_strand_tag = "auto"` already handles
  discovery, and `get_is_strand_qualified()` already requires a defined `sj_strand`, so such fragments are
  already excluded — but the junction table will then be small, and Gate 2 needs enough junctions to be
  meaningful. Report the junction count.
* **An unstranded library** (κ ≈ ½) makes the dispersion nearly unidentifiable, which is correct — the
  information weighting and the prior should handle it. Verify on `ss_0.50` synthetic conditions.
* **`od_r` falling from 0.2 to ~0.003 sharpens the RNA strand likelihood ~60×.** This is the intended
  correction, but it is a large change in one term. Gate 4 is not optional.
