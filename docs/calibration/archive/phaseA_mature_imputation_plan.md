# Phase A — mature-density imputation: dry-run implementation plan

**Status:** dry-run plan, 2026-06-12. First phase of
[`deconvolution_roadmap.md`](deconvolution_roadmap.md). Goal: build + validate (in isolation, against the
oracle) the **per-strand mature-density imputation** — the RNA mirror of `density_model`. **No** subtraction
and **no** strand-resolution yet (Phases B/C). This doc walks the implementation end-to-end, documents
exactly what must be done, and flags every issue that could arise during execution.

## 1. Goal & success criteria

Produce, per region `r` and strand `s ∈ {+,−}`, a predicted **contained-unspliced mature count** `M_{r,s}`,
imputed from the spliced crossings at eligible splice-junction boundaries. Validate it *standalone* against
oracle truth before it touches `calibrate`.

**Success:**
1. On the flagship failing regions (locus 21: 226/231/236/240), predicted `M_s` matches the oracle
   contained-unspliced mature within ~10–15%, so the implied residual fraction `(U − ΣM)/U` recovers the
   oracle gDNA fraction (0.56–0.72, vs the current 0.37–0.40).
2. On the gDNA-dominated regions (sig=6: 224/234/242), `M_s ≈ 0` (no spurious mature → no over-subtraction).
3. The toy AMBIG scenario (`tests/calibration/test_ambig_scenario.py` genome) behaves consistently.

## 2. The key realization — this *is* `density_model`, on spliced crossings

The mature imputation is **structurally identical** to the existing gDNA density imputation
(`density_model.node_gdna_density`); only the inputs change. This is exactly the symmetry the design calls
for ("the same imputation method"):

| element | gDNA density (existing) | **mature density (Phase A)** |
|---|---|---|
| anchor boundary | count-observable (no shared exon) | **eligible splice junction** (per strand) |
| crossing count | unspliced (`pos+neg`) | **spliced** (the junction-strand's mature) |
| crossing eff-length | `fl_mean` (gDNA FL) | `fl_mean_RNA` (RNA FL) |
| region eff-length | `region_eff_length(L, gdna_pmf)` | `E_mu = region_eff_length(L, rna_pmf)` |
| interior fill | `runfill_bidirectional` | **same** `runfill_bidirectional` |
| output | gDNA density → `count_gdna_frac` | per-strand mature density → `M_{r,s}` |

Crucially, **averaging the two bounding boundaries' densities** (what `density_model` already does:
`density[i] = mean(left, right)`) reproduces the validated `÷(n_junctions·fl_mean)` geometry — a region
with two bounding junctions gets `(S_left + S_right)/(2·fl_mean)`, the back-solved factor. So there is **no
new geometry to invent**; the divisor falls out of the existing per-boundary-density averaging.

## 3. Inputs (all already available)

- `substrate.{contained,left,right}`: `n_unspliced_pos/neg`, **`n_spliced_sense/antisense`**,
  `mass_unspliced`, `mass_spliced` (per region; left/right are the boundary sides inside region `r`).
- `region_arrays`: `signature`, `strand_class`, `ref_id`, `region_size_bp`.
- `κ = rna_sense_frac` (from `fit_strand_balance`).
- RNA FL pmf → `fl_mean_RNA = boundary_eff_length(rna_pmf)`; `E_mu = region_eff_length(region_size_bp,
  rna_pmf)`.
- `splice_junction_eligibility(sig_l, sig_r) → SpliceJunction(exon_side, strands)` — already identifies
  which boundaries are clean junctions **and the strand(s) that splice there**.

## 4. Dry-run walkthrough

### Step 4.1 — per-strand junction anchors

For each internal boundary `(r, r+1)` (same ref), call `splice_junction_eligibility(sig[r], sig[r+1])`.
If eligible, it names the **exon side** (`"L"` → region `r`; `"R"` → region `r+1`) and the **splicing
strand(s)** (`SpliceJunction.strands`, a tuple of `TS_POS`/`TS_NEG`). Record, per strand `s`, a boundary
anchor for the exon-side region carrying the spliced crossing count on that region's side:

- exon side `"L"` (region `r`): crossing count = `substrate.right[r].spliced_total` (region `r`'s right side).
- exon side `"R"` (region `r+1`): crossing count = `substrate.left[r+1].spliced_total`.

where `spliced_total = n_spliced_sense + n_spliced_antisense`. **Anchor density** = `spliced_total /
fl_mean_RNA` (the mature rate at that junction, on strand `s`).

→ **ISSUE A (spliced strand semantics).** Confirm in execution: at an eligible junction `strands=(s,)`, is
`spliced_total` entirely strand-`s` mature? The eligibility's matched-set rule yields a *single* strand for
a clean junction, and the splice motif (GT-AG) is directional, so the spliced reads there should all be
that junction's mature regardless of read mapping orientation. Verify `n_spliced_sense`/`n_spliced_antisense`
semantics (sense relative to *what*) and that their sum is the junction-strand mature crossing. If a junction
can legitimately carry both strands (overlapping antiparallel junctions), handle per-strand separately.

### Step 4.2 — per-strand mature density field (mirror `node_gdna_density`)

For each strand `s` independently, build a density array `m_s[R]`:
- a region that is the **exon side of an eligible strand-`s` junction** takes the **mean** of its available
  strand-`s` anchor densities (left and/or right) — exactly `density_model`'s `mean(left,right)`.
- regions with no strand-`s` junction anchor: `runfill_bidirectional(m_s, ref_id)` carries the nearest
  anchored neighbour inward (interior exon sub-regions — see ISSUE C).
- regions still unset (no strand-`s` anchor anywhere in the ref): `m_s = 0` (no mature evidence → predict no
  mature; this is the safe direction — never invent mature to subtract).

### Step 4.3 — predicted contained-unspliced mature count

`M_{r,s} = m_s[r] · E_mu[r]`. This is the count to (eventually, Phase B) subtract — distributed across the
observed strand counts by `κ` (Phase B, not here): `κ·M_{r,s}` to strand-`s` sense, `(1−κ)·M_{r,s}` to its
antisense. **Phase A stops here** — it produces `M_{+}`, `M_{−}` per region and validates them.

## 5. Module shape (code-health)

New module `calibration/mature_density.py` — the RNA mirror, single-purpose, no gDNA entanglement:

```
def mature_density(substrate, region_arrays, rna_region_eff_len, fl_mean_rna)
    -> MatureDensity        # dataclass: m_pos[R], m_neg[R], mature_pos[R], mature_neg[R],
                            #            junction_anchor_pos[bool R], junction_anchor_neg[bool R]
```

Refactor opportunity (do if clean, not forced): the anchor→mean→runfill→global-fallback skeleton is shared
with `density_model`. If a small private helper (e.g. `_impute_field(anchor_density, anchor_mask, ref_id)`)
can serve both without contortion, extract it; otherwise keep them parallel and readable. Decide during
execution — **readability over premature DRY.**

`splice_junction.py`: Phase A only *reads* `splice_junction_eligibility` (keep). Its
`region_splice_gdna_frac` (the fraction-at-boundary) is **untouched in Phase A**, retired in Phase B.

## 6. Validation harness (standalone, before wiring)

Reuse the locus-21 forensic scaffold (`scripts/debug/phase3_flagship_debug.py` helpers).

1. **Oracle per-strand contained-unspliced mature.** Extend `oracle_contained_unspliced` to split the
   mature (`mrna` origin) count **by transcript strand** → `oracle_mature_pos[R]`, `oracle_mature_neg[R]`
   (read origin → transcript → strand via `idx.t_df`).
2. **Compare** `M_{r,s}` (predicted) vs `oracle_mature_{s}[r]`, per region, on the flagship condition +
   the toy AMBIG genome. Report: the four failing regions, the three good regions, and aggregate MARD.
3. **Residual check (the payoff):** `(U_r − M_{+,r} − M_{−,r}) / U_r` vs the oracle gDNA fraction — this is
   the number that must move 0.37→0.56 on region 226. (Phase A computes it as a *diagnostic*; the actual
   subtraction + resolution is Phases B/C.)
4. New debug script `scripts/debug/phaseA_mature_imputation.py` (committed); a focused unit test on the toy
   genome locking M≈0 on gDNA-only regions and M≈oracle on mature-rich regions.

## 7. Issues & risks (to resolve during execution)

- **A. Spliced strand semantics** (§4.1) — confirm `n_spliced_sense+antisense` at an eligible junction = the
  junction-strand mature crossing; confirm orientation convention. *Highest-priority verification.*
- **B. Per-junction eff-length = `fl_mean_RNA`** — first cut; the back-solve showed ~10% scatter (510–584
  vs 514). Likely second-order (junction position within FL of a short exon, FL truncation). Defer
  refinement unless §6 validation exceeds the 10–15% bar; if so, derive the junction-crossing eff-length
  carefully (pure geometry, no constant).
- **C. Interior exon sub-regions** — a long exon split into pieces (by overlap with the opposite-strand
  gene) has interior regions with **no bounding splice junction** (their boundaries are exon↔exon, not
  exon↔intron → `eligibility` returns None). These must inherit `m_s` via `runfill_bidirectional` from the
  exon's junction-anchored edge regions. Verify the run-fill reaches them and stays within the exon (it is
  local + same-ref, capture-agnostic — no gate, per roadmap §6).
- **D. Substrate fractional crossing vs true crossing** — the accumulator deposits *fractional* mass; the
  substrate's spliced crossing (3,146) ≠ the oracle's strict crossings (1,944) at one junction. The density
  interpretation needs `spliced/fl_mean` to be an **unbiased mature-rate estimator** (the same property
  `density_model` relies on for gDNA: "the accumulator deposits the molecule's true span, so the one-side
  crossing flux is an unbiased density estimator"). Verify this holds for the **spliced** channel; if the
  spliced crossing is binned differently (e.g. by exonic length vs genomic span), correct the eff-length to
  match (cf. the prior gDNA FL span fix, `accumulator_fragment_span_redesign.md`). *This is the most likely
  real bug source* — the 2× was clean, but the per-region scatter may live here.
- **E. Per-strand separation in AMBIG** — a region near both a `+` junction and a `−` junction must anchor
  `m_+` and `m_−` independently from the respective junctions. The `strands` tuple drives this; keep the two
  fields fully separate.
- **F. `n_junctions_s = 0`** — region with no strand-`s` junction and no run-fill reach → `m_s = 0`. Safe
  (predict no mature), but note it: a genuinely mature-rich region with no nearby eligible junction (rare:
  single-exon transcript, or all junctions blacklisted) gets no mature estimate and will fall to Phase C's
  strand resolution.
- **G. `E_mu` region vs exon length** — `E_mu = region_eff_length(region_size_bp, rna_pmf)` uses the
  **region** length, consistent with how `mass_unspliced` is accumulated per region. Confirm the predicted
  `M` (region-scoped) is compared to region-scoped oracle counts (it is, in §6). No exon-length leakage.
- **H. True AMBIG-*exon* regions** — `splice_junction_eligibility` rejects a junction whose exon side has
  *both* exon strands (`exon_strands != matched`). Locus 21's failing regions are single-exon-strand
  (sig=9 = `E−|I+`: a `−` exon over a `+` intron), so eligibility **works** there (confirmed in the
  dissection trace: region 226 left boundary → `SpliceJunction(exon_side='R', strands=(−1,))`). A region with
  genuinely overlapping `+` *and* `−` exons (`E+|E−`) is a harder case deferred to Phase C; flag if it
  appears in the suite.

## 7.5. Execution outcome (2026-06-12 — COMPLETE)

**Built:** `calibration/mature_density.py` (`mature_density() → MatureDensity`, per-strand). Run-fill stays
within a contiguous same-strand-exon run (`_exon_run_segments`), never bridging an intron.

**Issue D — RESOLVED, no fix needed.** `scripts/debug/phaseA_issueD_mature_density.py` (clean multi-exon
gene, no gDNA/nascent): `M̂ = mean(bounding-junction flux)/fl_mean_RNA · E_mu` recovers the
contained-unspliced mature within ~5% across exons 400–1500 bp (ratio 0.96–1.05), and predicts ~0 for the
sub-FL 150 bp exon (`E_mu→0`). The accumulator spliced flux `÷ fl_mean_RNA` **is** an unbiased mature-rate
estimator. The earlier 3,146-vs-1,944 gap was fragment-level (accumulator) vs read1-only (crude oracle)
counting — not a bias. Issue A confirmed by the user (one boundary = all spliced crossings into the region).

**Payoff validated** on the flagship locus 21 (`scripts/debug/phaseA_mature_imputation.py`): predicted
per-strand mature matches the oracle, and the residual fraction `(U−M₊−M₋)/U` moves the failing regions
from the broken **0.37–0.40 → 0.54–0.57**, matching the oracle gDNA fraction (0.57–0.60); the
gDNA-dominated sig=6 regions stay at 0.97–0.98 (no over-subtraction). Locked by
`tests/calibration/test_mature_density.py` (2 tests). 209 calibration tests pass.

**Carry-forward to Phase B/C:**
- A small ~5–10% systematic **under**-prediction of mature remains (within the 10–15% bar) — the residual
  strand resolution (#3) + the gDNA crossing (#1) should tighten it; revisit only if it limits the leak fix.
- **Zero-mass slivers** (e.g. locus-21 regions 227/232, `U≈0`) get a tiny run-filled `M>0` → a nonsense
  residual when divided by `U`. **Phase B must clamp `M₊+M₋ ≤ U`** (and treat `U≈0` regions as no-op).

## 8. Explicitly out of scope for Phase A

- The strand-aware **subtraction** producing `U'_pos/U'_neg` → **Phase B**.
- Retiring `region_splice_gdna_frac` → **Phase B**.
- **Strand resolution** of the residual / AMBIG informativeness / nascent prior → **Phase C**.
- Wiring into `calibrate`, golden, suite → **Phase D**.
