# Counts, densities, and the conversions between them — the audit

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Written 2026-08-10 against `a56e139c`, at the owner's request: audit every place the tool uses a
    count, every place it uses a density, and every place it converts between them — and settle the
    paradigm.

---

## 1. ⭐⭐⭐ THERE ARE THREE CURRENCIES, NOT TWO — AND THAT IS THE WHOLE SOURCE OF THE CONFUSION

The question was posed as counts *versus* densities. The tree actually carries **three** distinct
quantities, and every defect found in this campaign is a place where two of them were treated as one.

| currency | what it is | where it lives | conserved? |
|---|---|---|---|
| **INCIDENCE** | `+1` on every object the fragment touched | `node_contained_count`, `edge_unspliced_count`, `edge_spliced_count`, `sj_count` | ⛔ **NO** — one fragment books `max(K,1)` of them |
| **FRAGMENT** | a coverage-weighted share summing to **exactly 1** per fragment | `edge_unspliced_mass`, `edge_spliced_mass` | ✅ **YES**, by construction |
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
| **NODE** | ✅ **FREE** — containment is exclusive, they are the same number | ⛔ needs a model: the opportunity `(ell−w+1)₊` does not cancel `1/w` |
| **EDGE** | ⚠ needs `q` — but `q` is **MEASURED**, not modelled | ✅ **FREE** — `inv_length_sum` **is** `ρ`, exactly |

> **On each axis two of the three conversions are free, and they are different two. There is no axis on
> which all three are free — and no axis on which a model is actually required.**

⛔ **And the tree currently takes the model-dependent path on BOTH axes.** That is the paradigm defect: not
that a choice was made badly, but that the free path exists on each axis and is unused.

⚠ **`inv_length_sum` at a NODE is not a density and must never be used as one.** `scan_payload.py` says
so explicitly, which is why the bank is not called `density`. Do not "extend the model-free density to
nodes" — it does not exist there.

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
4. ⛔ **There is no `sj_mass`**, so the junction axis cannot be converted to fragments at all and a fully
   conserved library fragment count is **not computable today**.

## 6. ⭐⭐⭐ THE RULING THIS NEEDS — one invariant, and a name for each field

> **A FRACTION IS A COUNT SHARE AND MULTIPLIES A CONSERVED FRAGMENT COUNT. A DENSITY IS DERIVED FROM
> THAT, NEVER THE REVERSE. AND EVERY FIELD NAMES ITS CURRENCY.**

* **FRAGMENT** wherever the answer is a count — the prior's pseudo-counts, the library summary, the QC
  report. It is model-free end to end: the mass is measured, needs no divisor, inherits no pmf error.
* **DENSITY** only where a density is genuinely the quantity — the enrichment landscape, the NegBinom
  peel, the `gdna_eff_len` contraction, the relay's currency. There the model dependence buys something.
* ⛔ **INCIDENCE is an internal representation and must never leave the calibration layer un-converted.**

**The three edits that implement it**, in order of value:

1. `pipeline.py` applies `q` to its edge terms, and `q_spliced = edge_spliced_mass / edge_spliced_count`
   to the spliced half. ⭐ **Free, arithmetic, worth 4.3 pp of the deliverable.**
2. Rename `mass_*_edge` → `count_*_edge` (or convert at the source in `sweep.py` so the name becomes
   true). Defect 2 is what made defect 1 invisible.
3. Decide on `sj_mass` — an accumulator ADDITION, `JunctionEdge` 16 → 24 B. Without it the library
   fragment count has an un-convertible third term.

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
