# Counts, densities, and the conversions between them — the audit

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Written 2026-08-10 against `a56e139c`, at the owner's request: audit every place the tool uses a
    count, every place it uses a density, and every place it converts between them — and settle the
    paradigm. ✅ **The plan in §9 is COMPLETE (2026-08-11)** and its rulings have moved to `DESIGN.md`,
    its lessons to `TRAPS.md` and its numbers to `ROADMAP.md` §0. ⚠ Sections 1–8 are the ORIGINAL audit
    and describe the tree as it was; where they disagree with the permanent docs, the permanent docs win.

---

## 1. ⭐⭐⭐ THERE ARE THREE CURRENCIES, NOT TWO — AND THAT IS THE WHOLE SOURCE OF THE CONFUSION

The question was posed as counts *versus* densities. The tree actually carries **three** distinct
quantities, and every defect found in this campaign is a place where two of them were treated as one.

| currency | what it is | where it lives | conserved? |
|---|---|---|---|
| **INCIDENCE** | `+1` on every object the fragment touched | `node_contained_count`, `edge_unspliced_count`, `edge_spliced_count`, `sj_count` | ⛔ **NO** — one fragment books `max(K,1)` of them |
| **FRAGMENT** | a coverage-weighted share summing to **exactly 1** per fragment | `edge_unspliced_mass`, `edge_spliced_mass`, ⭐ `sj_mass` (added 2026-08-11 — without it the sum was 0.747 on RNA) | ✅ **YES**, by construction |
| **DENSITY** | fragments per base of opportunity, `ρ` | `edge_unspliced_inv_length_sum`, `sj_inv_length_sum` | n/a |

⭐ **The mass rule is coverage-weighted, not `1/K`**, and the executable specification says why: *"Both
conserve; only this one says WHERE the fragment sat, and only this one is expressible per base — which is
how the two are told apart at all."* So the mass is a fragment count **attributed by coverage**, which is
what makes it the right thing to split by a composition fraction.

⛔ **`count` is NOT a fragment count**, except at a NODE — where containment is exclusive, so a contained
fragment deposits on exactly one node and incidence and fragment coincide. That single exception is why
the confusion survives: half the tree's counts really are fragment counts, and half are not.

## 2. ⭐⭐ THE CONVERSION MATRIX — and which conversions are MODEL-FREE

| from | to | operator | model-free? | where it happens today |
|---|---|---|---|---|
| incidence | fragment | `× q`, `q = mass/count` | ✅ **measured** | `priors.py` — and **nowhere else** |
| incidence | density | `÷ eff_c` | ⛔ a fitted-pmf functional | `node_geometry`, `density_model`, `node_init`, `npmle`, `background_reference`, `gdna_landscape`, `messages/head` |
| density | incidence | `× eff_c` | ⛔ same | `capture_eff_length`, `priors.gdna_eff_len` |
| **incidence** | **density** | **`= inv_length_sum`** | ✅ **EXACT at an EDGE** | ⛔ `second_pass.py` only — calibration never reads it |
| fragment | density | `÷ (mass / inv_length_sum)` | ✅ measured at an EDGE | nowhere |

⭐⭐ **The fourth row is the finding.** At a contiguous edge the opportunity is `w−1` and the deposit is
`1/(w−1)`, so they cancel identically and `E[inv_length_sum] = ρ` **exactly, for any length distribution**
(`accumulator.h`). The tool therefore already measures the edge density with no model — and calibration
computes the same quantity as `count / eff_gdna`, a functional of the fitted gDNA pmf, instead.

## 3. ⭐⭐⭐ THE SYMMETRY IS COMPLEMENTARY, AND IT IS THE PARADIGM

|  | incidence → fragment | incidence → density |
|---|---|---|
| **NODE** | ✅ **FREE** — containment is exclusive, they are the same number | ⚠ needs a model *as deposited* — `1/w` does not cancel `(ell−w+1)₊`. ⭐ A `1/A` deposit would make it free too (§8) |
| **EDGE** | ⚠ needs `q` — but `q` is **MEASURED**, not modelled | ✅ **FREE** — `inv_length_sum` **is** `ρ`, exactly |

> **On each axis two of the three conversions are free, and they are a different two. There is no axis on
> which a model is actually REQUIRED — the node's density needs one only because of the deposit rule we
> chose, and §8 prices changing it.**

⛔ **And the tree currently takes the model-dependent path on BOTH axes.** That is the paradigm defect: not
that a choice was made badly, but that the free path exists on each axis and is unused.

⚠ **`inv_length_sum` AS DEPOSITED TODAY is not a density at a NODE and must never be used as one.**
`scan_payload.py` says so explicitly, which is why the bank is not called `density`. ⭐ **But that is a
property of the DEPOSIT RULE, not of the node**: the node deposits `1/L`, and `1/L` does not cancel
`(ell−w+1)₊`. A deposit of `1/(ell−w+1)₊` would cancel it exactly — see §8, where that formula is
re-derived and priced. So the NODE row of the table above reads "needs a model" **for the bank we
currently store**, and is a choice rather than a law.

## 4. THE SITE AUDIT

**Counts used as counts (correct, no conversion):** the strand Beta-Binomial (`gdna_strand`,
`strand_deconv`) reads per-strand incidence columns and forms a ratio, so inflation cancels; Poisson
precision `1/n` in `node_init` uses the raw count as a count of observations, which is what it is.

**Counts converted to densities (`÷ eff_c`, model-dependent):** `node_geometry._rate` is the single
choke point and it is used correctly — 0 where the divisor is 0, never floored. Its consumers are
`density_model`, `density_deconv`, `npmle`, `gdna_landscape`, `background_reference`, `node_init`'s
`rho_c = f_c·M/E_c`, and `messages/head`'s transfer variance. ✅ **All of these genuinely want a density**
(the enrichment landscape, the NegBinom peel, the relay's message currency), so the model dependence is
paid for something.

**Fractions applied to counts:** `sweep.py` `gdna_mass = f_g·count` / `rna_mass = (1−f_g)·count + spliced`.
✅ Dimensionally right — `f_g` is a COUNT share (established: `Σ_c rho_c·E_c = M` forces it) — ⚠ but the
result is an **incidence**, and the field is named `mass_*`, which is the naming defect that let the next
one through.

**Incidence → fragment (`× q`):** `priors.py` only.

**Density → length:** `priors.gdna_eff_len` contracts each object as `min(m/ρ_ref, S)`. ⚠ Uses a single
global `ρ_ref`; under capture the four gDNA pools disagree by 44 bp on their own mean, so `ρ_ref` inherits
a placement residual. Stated, not fixed.

## 5. THE DEFECTS THIS AUDIT FOUND

1. ⛔⛔ **`pipeline.py:993-998` sums incidences and calls them fragments.** Already recorded at
   `a56e139c`: applying `q` to the edge terms moves the reported library `f_gdna` **0.3781 → 0.4214**
   against a truth of 0.5085, +4.3 pp **toward** truth off capture, on all four conditions tested.
2. ⛔ **`mass_gdna_edge` / `mass_rna_edge` are incidences under a `mass_` name**, while `mass_gdna_node`
   really is a fragment count. Same prefix, different currency, and defect 1 is the direct consequence.
3. ⚠ **`q` is measured pooled and applied per component**, so splitting a pooled mass by a count share
   carries the `q_g ≠ q_r` error — `Δphi` ≤ 0.6 pp on the total prior, placement-driven, not repairable
   in production. Recorded, bounded, accepted.
4. ✅ **CLOSED 2026-08-10 — `sj_mass` was built and the library count is now computable.** The defect was
   larger than this entry knew: the junction axis was not merely unconvertible, it was the ONLY axis a
   spliced fragment with no line-crossing block reached at all, so **1,222,375 of 4,830,713 RNA fragments
   (25.3 %)** deposited on no conserved bank. gDNA read `1.000x` deposited and RNA `0.747x`. With the bank
   in place both are `1.000x` with **0** unaccounted. Lesson: `TRAPS: an-identity-with-a-qualifier`.

## 6. ⭐⭐⭐ THE RULING THIS NEEDS — one invariant, and a name for each field

> **A FRACTION IS A COUNT SHARE AND MULTIPLIES A CONSERVED FRAGMENT COUNT. A DENSITY IS DERIVED FROM
> THAT, NEVER THE REVERSE. AND EVERY FIELD NAMES ITS CURRENCY.**

* **FRAGMENT** wherever the answer is a count — the prior's pseudo-counts, the library summary, the QC
  report. It is model-free end to end: the mass is measured, needs no divisor, inherits no pmf error.
* **DENSITY** only where a density is genuinely the quantity — the enrichment landscape, the NegBinom
  peel, the `gdna_eff_len` contraction, the relay's currency. There the model dependence buys something.
* ⛔ **INCIDENCE is an internal representation and must never leave the calibration layer un-converted.**

**The three edits that implement it**, in order of value:

1. ✅ **DONE 2026-08-10, and edit 3 had to land first.** ⛔ Applying `q` alone was measured to FLIP the
   error's sign rather than remove it — at ladder g50 capture_off `f_gdna` goes −0.1234 → **+0.0722**,
   because converting the axes without `sj_mass` drops the 25.3 % of RNA that deposits no conserved mass.
   The shipped fix is `CalibrationResult.library_{gdna,rna}_fragments`, each axis converted by its own
   population's `q`, read by `pipeline.py` and `calibration_truth_ab.py`.
2. ⏸ **STILL OPEN** — rename `mass_*_edge` → `count_*_edge`. ⚠ THREE of the five `mass_*` fields are
   incidences (`mass_gdna_edge`, `mass_rna_edge`, `mass_rna_junction`) and two are fragment counts; this
   entry said two. Deliberately NOT bundled with edit 1 (`one-thing-varied`).
3. ✅ **DONE** — `sj_mass` landed, `JunctionEdge` 16 → 24 B exactly as costed.

## 7. ⭐⭐⭐ THE FREE GATE — the tool can measure its own length-model error, with no oracle

**The construction.** At a line whose two flanks are `{intron, exon}` or `{intergenic, exon}`, mature RNA
cannot cross — leaving an exon means using a junction, which puts the fragment in `edge_spliced`
(`TRAPS: mature-rna-never-crosses-a-seam`). So a seam line's unspliced crossings are gDNA (plus nascent
RNA, where any exists). Two independent estimates of the same `ρ` are then available at that slot:

    inv_length_sum        = rho                          EXACT, no model
    count / eff_gdna      = rho * E_true[w-1] / E_fit[w-1]   the fitted gDNA pmf

    =>   count / (eff_gdna * inv_length_sum)  =  E_true[w-1] / E_fit[w-1]

⭐ **A per-slot, oracle-free, simulator-free measurement of the fitted gDNA fragment-length model's bias
— computable on real data.**

**MEASURED** (count-weighted over seam lines, cached `_main` payloads):

| condition | seam lines | fragments | ratio | implied bias |
|---|---|---|---|---|
| ladder g50 capture_off | 22,227 | 241,419 | 0.9937 | **+1.4 bp** |
| ladder g50 capture_on | 14,661 | 1,387,267 | 0.9891 | **+2.8 bp** |
| ladder g98 capture_off | 22,230 | 472,394 | 0.9937 | **+1.4 bp** |
| ladder g98 capture_on | 17,375 | 2,713,225 | 0.9891 | **+2.8 bp** |

⭐⭐ **It is internally consistent in exactly the way it should be**: identical at g50 and g98 (it is a
property of the library preparation, not of the gDNA level), the fitted pmf reads **long** in both arms,
and the bias **doubles under capture** — the direction and the capture-dependence both match `ROADMAP.md`
§0's recorded capture defect (+13.6 bp at a 330 bp mean, +3.5 bp at 120 bp; these libraries sit at ~235 bp).

⛔ **Its precondition is NO NASCENT RNA**, which holds on every panel here (`nrna_none`) and on cfRNA
(cell-free means no active transcription) but **not on tissue**. State it, do not assume it.

⭐ **Why this matters beyond the audit.** `ROADMAP.md` §0 calls the gDNA pmf under capture "the last
unfixed length defect" and rules that it must not be attacked by editing `gdna_opportunity` because it is
a PLACEMENT problem. This measures that residual directly, at runtime, on the user's own library — which
is the instrument that problem has been missing. It is also the natural home for a shipped QC number.


---

## 8. ⭐⭐⭐ THE MODEL-FREE NODE DENSITY — it exists, it was derived, and the answer is "yes but not as an addition"

**Owner, 2026-08-10:** *"Previously we derived the model-free density for NODES. It's a different derived
formula."* ⭐ It is, and it is already gated: `scripts/design/node_density_derivation.py` (T0–T6, each
perturbed) and the `reciprocal-opportunity-deposit` memory. Re-derived here independently and it matches.

**The formula.** The deposit weight is `1/OPPORTUNITY`, not `1/length`. At an object whose opportunity for
a length-`w` fragment is `A(w)`, depositing `1/A(w)` gives

    E[ SUM 1/A ]  =  SUM_w  rho * A(w) * f(w) * (1/A(w))  =  rho     exactly, for ANY pmf

At an EDGE `A(w) = w−1`, a function of `w` alone — which is what ships. **At a NODE
`A(w) = (ell − w + 1)₊`, which depends on the NODE LENGTH as well** — and the accumulator knows `ell` at
deposit time, so it is depositable. That is the different formula: an edge is the `ell → 0` limit of a
node. ⭐ Measured in the derivation: at a 1000 bp node with 200 bp fragments `count/ell` reads **0.80×**
truth while `Σ1/A` reads **1.000000×**.

### 8.1 ⛔ TWO RECORDED VERDICTS DECIDE THIS, AND NEITHER IS ABOUT THE DENSITY

**T3 — model-free and informative are MUTUALLY EXCLUSIVE, and it is a theorem.** Model-free means both
components share the coefficient, so the channel's row is `(K, K)` and its determinant against any other
row is zero. **A channel cannot both be model-free and tell gDNA from RNA.** Measured: `Σ1/A` alone scores
0.113 discrimination efficiency at a 151 bp node and **0.000** at 1000 bp.

⭐ **So it buys a LEVEL, never a SPLIT** — and the split is what calibration is for.

**And the second one is the honest limit of the whole model-free idea**: `Σ1/A = rho` is model-free with
respect to the **LENGTH distribution**, not with respect to **PLACEMENT**. The cancellation assumes `rho`
is uniform across the admissible start positions. At an EDGE that window is about one fragment length; at
a NODE it is the whole node, which is kilobases. ⛔ **So the node version does NOT solve the capture
placement problem** — it removes the length-model layer sitting on top of it, and leaves the placement
residual exactly where `EQUATIONS.md` §4.4 says it is.

### 8.2 WHAT IT WOULD ACTUALLY BUY, END TO END

| need | does the model-free node density help? |
|---|---|
| **the EM prior's pseudo-counts** | ⛔ **No.** At a NODE incidence = fragment already; the node term needs no conversion and no density. |
| **the composition split `f_g`** | ⛔ **No, provably** — T3. |
| **the density surface** (`density_model`'s anchor, `background_reference`, `npmle`, `gdna_landscape`, `node_init` source 1, `priors`' `rho_ref`) | ⚠ **Yes, but by ~1 %.** The shipped path is `count / eff_gdna`, which ALREADY applies the effective-length correction — it just uses a FITTED pmf to do it. So the gain is exactly that pmf's error: **0.9937 off capture, 0.9891 under it** (§7). |
| **diagnostics** | ⭐⭐ **Yes, and this is the real prize.** It extends §7's oracle-free gate from seam LINES to every structurally pure-gDNA NODE — the intergenic and intronic anchors — which is a far larger population and is where the prior-free pass's anchor actually lives. |

⭐ **The one number that frames the decision**: the model-free node density is worth about **1.1 %** on the
density surface, because the shipped divisor is already doing the right correction with a slightly wrong
pmf. It is not a correctness fix; it is the removal of a dependency, plus a much better diagnostic.

### 8.3 ⛔⛔ THE SEQUENCING POINT THAT CHANGES THE CLEANUP PLAN

`node_contained_inv_length_sum` deposits `1/L`. That is **not** `1/A` at a node, which is exactly why it
is model-free at an edge and not at a node, and exactly why it ended up with no consumer.

> **So the choice is not "keep or delete". It is "delete it, or fix its deposit rule" — the same 8 bytes
> either way.** Changing `1/L` to `1/(ell − w + 1)₊` turns a dead bank into the model-free node density
> at zero net memory.

⛔ **And both are cache-invalidating, so doing them separately costs TWO full re-scans of every panel.**
Deleting the bank now and re-adding a differently-deposited one later is the expensive path.

**Revised bank ledger** if the deposit rule is fixed rather than the bank deleted:

| | delete all three | fix the deposit instead |
|---|---|---|
| `Node` | 24 → **8 B** | 24 → **16 B** |
| `ContiguousEdge` | 48 → 40 B | 48 → 40 B |
| model-free node density | ⛔ gone, and re-adding costs a second invalidation | ✅ available |

⛔ **ONE BLOCKER IS RECORDED AND UNRESOLVED**: fixed-point headroom. `A` can be **1** (a fragment exactly
filling its node), so the quantum can be `2³²` — far above the `L ∈ [20, 2000]` range the `uint64`
accumulator was sized for. ⚠ `DESIGN.md` §3.1's ~800× headroom figure is computed for `1/L` and does not
carry. **Re-price it before any of this ships.** And `Σ1/A` is ~3× noisier at the crossover `ell ≈ E[L]`,
which is where exon nodes sit (~98 bp median against ~235 bp fragments) — ⭐ but the anchors are
intergenic and intronic nodes, which are long, so the crossover lands away from the population that
matters.

### 8.4 RECOMMENDATION

⭐ **Fix the deposit rule rather than delete the bank — but do it as ONE event with the other two
removals, and clear the fixed-point blocker first.** The correctness gain is small (~1 %) and honest;
the diagnostic gain is large; the memory cost is zero; and the alternative sequencing pays the
cache-invalidation cost twice.

⛔ **If the fixed-point headroom cannot be cleared cheaply, delete all three and accept that the
model-free node density costs a second invalidation later.** It is worth ~1 % and a diagnostic, not a
re-scan of every panel on its own.


---

## 9. ✅ THE IMPLEMENTATION PLAN — COMPLETE (2026-08-11)

All six phases landed, plus two the plan did not anticipate. The rulings moved to `DESIGN.md`, the lessons
to `TRAPS.md` (`an-identity-with-a-qualifier`, and a second occurrence of `could-the-arm-have-fired`), the
numbers to `ROADMAP.md` §0. ⛔ Nothing here is authoritative; read those.

* the node deposit is `1/A`, so the node channel is a density
* `sj_mass` exists, so a spliced fragment that crosses no line is countable — and the rule was then
  extended so a junction edge is a BOUNDARY like any other, claiming at both its positions
* the library figure is a conserved FRAGMENT count, derived on read, with one home
* ⭐ the deposit-behaviour digest closes `a-hash-that-misses-its-artifact` for deposit-rule changes
* ⭐ ONE NUMERIC CONVENTION — a COUNT is an integer, a FRACTION is float64; the fixed point is gone
* the panels are re-scanned and gated

⚠ **What the plan got wrong, worth keeping:** it ranked "apply `q` to the edge terms" first at +4.3 pp.
Measured against per-fragment truth that fix alone FLIPS the error's sign (−0.1234 → +0.0722), because
25.3 % of RNA deposits no conserved mass for `q` to convert. The cheap fix was not a partial version of
the right one.
