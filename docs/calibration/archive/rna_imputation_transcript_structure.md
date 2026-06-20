# RNA imputation must respect transcript structure (the corrected RNA design)

**Status:** authoritative design correction (2026-06-17). Supersedes the RNA-imputation framing in
`CALIBRATION_PLAN_v6.md` §8 (which imputed *total* RNA across *all* edges and is now known wrong). Grounded in
the diagnostic + code-trace + data-audit workflow (`scripts/debug/rna_imputation_diagnostic.py`). This is the
reference for *what the RNA model should be* and *what the next implementation steps are*.

---

## 0. TL;DR

The cross-node RNA imputation in v6 Phase B was built as the structural twin of the gDNA imputation: carry a
neighbour's **total** same-strand RNA density (unspliced + spliced) across **any** flanking edge. That premise
is wrong, because **RNA is not a single smooth field** — it decomposes into two species that flow differently:

- **Mature RNA** is *contained* in exons (dense in the exon body) and crosses a splice junction **only as
  spliced** (one-sided — it jumps to the next exon). It is a **fixed, directly-observed anchor**, never
  imputed, never deconvolved.
- **Nascent RNA** (unspliced pre-mRNA) is the genuinely *contiguous* species — it flows exon→boundary→intron
  as **unspliced**. This — and only this — is what may impute across an exon↔intron edge.

Transcript structure gates the imputation in three ways:
- **(a) TSS/TES edges (intergenic↔exon) are zero-RNA black holes** — no RNA imputation crosses a transcript end.
- **(b) Exon↔intron edges carry the UNSPLICED (nascent) component only** — never the one-sided spliced.
- **(c) An exon's contained mature is anchored by the spliced floor (a direct observation), not imputed from
  the thin junction crossing.**

The v6 Phase-B imputation violates all three. The measured consequences (zero-gDNA library, true f_g=0
everywhere): **introns fabricate ~75% gDNA** (the RNA prior is gated OFF there because the spliced skips the
intron, so the nascent is invisible), and in dense data **single-strand exons over-call gDNA** (the imputation
under-counts the contained mature and over-rides the reliable strand).

---

## 1. The RNA model (what the pie's RNA actually is)

Calibration solves the **unspliced** pie `{f₊, f₋, f_g}` per node (RNA+ / RNA− / gDNA fractions of the node's
unspliced mass). Spliced mass is separate — a fixed RNA floor. The unspliced RNA (`f₊+f₋`) at a node is the sum
of two species:

| species | where it lives | how it crosses a splice junction | role in calibration |
|---|---|---|---|
| **mature** (exon body) | *contained* in exons (dense) | **spliced only** (one-sided, exon flank) | a **fixed floor**/anchor — directly observed, never deconvolved |
| **nascent** (pre-mRNA) | exon **and** intron, contiguous | **unspliced** (contiguous, both sides) | **imputable** across exon↔intron (the only RNA that flows) |
| gDNA | genomic, smooth | unspliced (two-sided) | the smooth field (existing gDNA count prior) |

Per node type:
- **Exon region, contained unspliced** = mature-body + nascent (+ gDNA). Dense mature-body.
- **Intron region, contained unspliced** = nascent (+ gDNA). No mature (spliced out).
- **Exon↔intron boundary, UNSPLICED crossing** = nascent (+ gDNA). *Mature cannot cross unspliced here.*
- **Exon↔intron boundary, SPLICED crossing** = mature, **one-sided** (exon flank).

**Consequence:** the boundary's **unspliced crossing density IS the clean nascent signal** (mature is excluded
by construction — it would be spliced). This is the quantity to impute onto the intron.

---

## 2. The three transcript-structure rules

Worked from the user's examples (single-exon, two-exon transcripts):

### (a) TSS/TES (intergenic↔exon) = a zero-RNA black hole
A transcript is a linear node span bounded by a TSS and a TES; **the anchoring TSS/TES boundaries are not part
of the transcript.** Fragments crossing them belong to other things (other exons, overlapping introns, gDNA).
So **no RNA imputation flows across a TSS/TES edge** — the intergenic side is a *precise zero-RNA* prior, not
merely "no positive coupling." (Today the sweep happens to set `Q=∞` across such edges, which blocks the
log-odds *propagation* but is NOT an enforced zero-RNA prior on the *imputation* — weaker than rule (a) wants.)

### (b) Exon↔intron = impute the UNSPLICED (nascent) only
At a splice-junction boundary:
- **Spliced** (one-sided, exon flank) = **fixed mature anchor**. Never imputed, never deconvolved. It *anchors*
  the partition of the boundary's unspliced crossings (the spliced floor) and the adjacent exon's contained
  mature (via the existing `region_splice_gdna_frac` count-fraction).
- The **unspliced** crossing (nascent + gDNA) is the thing solved — and the **nascent** component imputes
  **contiguously** across the edge (exon↔intron, both directions). **The spliced does NOT cross into the
  intron** (it jumps to the next exon), so imputing *total* RNA across exon↔intron is wrong.

### (c) An exon's contained mature is anchored, not imputed
An exon's mature-body is *dense and contained*; only a thin slice crosses the junction as spliced. So the
junction crossing **under-represents** the exon's RNA (measured ~4×). The exon's mature magnitude comes from
(i) the **strand likelihood** on its contained reads (reliable at good SS) and (ii) the **spliced floor** (the
junction spliced says "≥ this much mature"), **not** from imputing the thin junction crossing onto the exon.

---

## 3. The evidence (diagnostic + code trace + data audit)

### 3.1 Diagnostic (zero-gDNA, true f_g=0 everywhere) — `scripts/debug/rna_imputation_diagnostic.py`
Genome: a single-exon transcript, a two-exon transcript, and multi-exon +/− genes; gDNA=0, nascent~30, SS=0.99.
Result: total RNA 5135, **but 1032 mass fabricated as gDNA — ALL on the 5 intron nodes** (f_g 0.70–0.75); every
exon and intergenic node reads f_g=0.

| node | sig | f_g | why |
|---|---|---|---|
| exon (single-exon, two-exon, multi-exon) | EXON± | **0.000** ✓ | strand likelihood at SS=0.99 drives RNA; RNA prior weak here |
| intergenic | NONE | **0.000** ✓ | no mass |
| **intron** (×5) | INTRON± | **0.70–0.75** ✗ | RNA prior **OFF** (no spliced on intron sides) → count clue (gDNA-by-signature) + global foundation fabricate gDNA from real nascent |

**The mechanism (traced):** every intron has `spliced_sense=0` on **both** intron-facing sides — the one-sided
spliced intron-skip lands on the **exon** flanks and skips *over* the intron. `node_rna_density`'s `*_ok` gate
requires same-strand `spliced>0` on the facing side, so the intron gets `rho_hat=NaN, τ_rna=0` → **RNA prior
entirely off**. The intron carries genuine nascent (own RNA density ~0.23–0.26, contained unspliced 173–688)
that the imputation can't see. **The intron-facing unspliced-RNA crossing densities (~0.18–0.25) ARE computed
but DISCARDED by the spliced-gate** — they are exactly the nascent signal rule (b) wants.

*(Separately, in dense real data — quick_3to1 — single-strand exons over-call gDNA (f_g~0.76): there the RNA
prior's τ is high enough to override the strand with the under-counted μ_s — rule (c) + the τ-query bug.)*

### 3.2 Code violations (current `rna_density_model` / `calibrate` / `simplex_sweep`)
1. **Rule b** — `rna_density_model.py:122-123`: `rho_hat_s` **includes** the spliced count (`+ l_spl`/`+ r_spl`).
   The imputed density carries the one-sided spliced across the junction.
2. **Rule b (the OFF-on-introns face)** — `rna_density_model.py:124-125`: `left_ok/right_ok` gate on
   `spliced_sense>0`, so the imputation only fires where a spliced anchor faces the dest. Introns (no spliced on
   their sides) get the prior turned off — the nascent is never imputed onto them.
3. **Rule a** — no structural TSS/TES exclusion; intergenic↔exon is `count_observable`, only the accidental
   spliced-gate blocks it. No precise zero-RNA prior.
4. **Rule c** — the imputed magnitude is the thin junction-crossing slice, used as the exon's RNA target →
   ~4× under-count → exon false gDNA in dense data.
5. **BUG (τ-query axis)** — `calibrate._rna_prior_fraction` queries the var~mean at the **side** density
   `rho_hat_s`, but the curve was **fit** at the **region** density (the gDNA path fits and queries at region
   density). Off-support → mis-calibrated, over-confident τ_rna → over-rides the reliable strand.
6. **BUG (strand collapse)** — `rho_hat` is built from `c_total` (pos+neg) + `n_spliced_sense` and shared across
   both strands; combined with the spliced inclusion, the own-strand target is anchored to the junction slice.

### 3.3 Data capabilities — all derivable, NO C++/payload change
- **TSS/TES identifiable**: `is_intergenic = (signature==0)` on one flank + an exon bit on the other; equivalently
  `splice_junction_eligibility(sig_l,sig_r) is None` with the non-exon side intergenic.
- **Exon↔intron identifiable**: `splice_junction.splice_junction_eligibility` returns `SpliceJunction(exon_side,
  strands)` — the existing, parameter-free predicate (it already names the exon side).
- **Unspliced/spliced separable**: channels 0/1 (unspliced) vs 2/3 (spliced); `SubstrateView.mass_unspliced`
  vs `mass_spliced`, per boundary side.
- **Spliced one-sided**: confirmed in the accumulator deposit + measured (100% boundary-owned, 97.8% one-sided).
- **Edge classification** ({TSS/TES black-hole | exon↔intron unspliced-only | within-transcript}) is fully
  derivable from the two flanking signatures + `boundary_region_indices` + `splice_junction_eligibility`. The
  classifier + the directed operators are the build; the data is sufficient.

---

## 4. The corrected imputation mechanism

**One principle:** the cross-node RNA imputation operates on the **nascent (unspliced-RNA) density only**, gated
by a per-edge transcript-structure classifier; mature is a fixed anchor (spliced floor + strand), never imputed.

1. **Per-edge classifier** (new, pure signatures + adjacency): each region↔boundary edge →
   `{TSS_TES_BLACKHOLE, EXON_INTRON_UNSPLICED, WITHIN_TRANSCRIPT}`.
2. **The imputed RNA density carries the UNSPLICED-RNA crossing density only** — drop the `+ spliced` term from
   `rho_hat`. Spliced stays the fixed floor (already correct in `simplex_sweep` as the one-sided lower bound).
3. **Un-gate the nascent across exon↔intron**: the imputation onto an intron uses the flanking boundary's
   **unspliced-RNA** crossing density (which is present), NOT gated on a spliced anchor (which is absent there
   by construction). This is the fix for the intron false-gDNA.
4. **TSS/TES black hole**: no RNA imputation across an intergenic↔exon edge; the intergenic side is a precise
   zero-RNA prior.
5. **Exon mature anchor (rule c)**: the exon's RNA magnitude is the strand likelihood + the spliced floor
   (the existing `region_splice_gdna_frac` count-fraction already converts the junction spliced to the exon's
   contained mature density) — the imputation must not under-write it with the thin junction crossing.
6. **Fix the τ-query** to the region density (consistent with the gDNA path + the fit axis), and ensure the
   imputation **defers to the strand** where the strand is informative (honest precision; the imputation is
   unreliable and must not over-ride a sharp strand likelihood — the user's core principle).

**Symmetry / traversal:** the imputation remains node-pair, bidirectional, traversing boundary↔region↔boundary;
the edge classifier just gates *what* (nascent only) and *whether* (not across TSS/TES) per edge.

---

## 5. Next steps (the implementation plan)

In order; each independently testable; the diagnostic scenario is the regression lock (true f_g=0 everywhere).

- **Step 1 — fix the two plain bugs** (cheap, isolated): (i) query τ_rna at the region density in
  `_rna_prior_fraction`; (ii) the strand-collapse / spliced-inclusion in `rho_hat` (folds into Step 2).
- **Step 2 — the per-edge transcript-structure classifier** (new helper, pure signatures + adjacency):
  `{TSS_TES_BLACKHOLE | EXON_INTRON_UNSPLICED | WITHIN_TRANSCRIPT}` from `boundary_region_indices` +
  `splice_junction_eligibility` + the `is_intergenic` predicate. Unit-test on the diagnostic genome.
- **Step 3 — nascent-only imputation across exon↔intron** (rule b): drop spliced from `rho_hat`; un-gate the
  nascent onto introns (use the unspliced-RNA crossing density, not a spliced anchor). *Gate: the diagnostic
  introns drop from f_g~0.75 toward ~0.*
- **Step 4 — TSS/TES black hole** (rule a): no RNA imputation across intergenic↔exon; precise zero-RNA on the
  intergenic side.
- **Step 5 — exon mature anchor honesty** (rule c): the imputation must not under-write the exon's contained
  mature; lean on the strand + the spliced-floor count-fraction; ensure the prior defers to a sharp strand.
- **Validation gate** (every step): the diagnostic (true f_g=0 ⇒ introns AND exons ≈0) + the net-flow
  before/after on the quick_3to1 capture-on conditions (gdna300 win preserved, zero-DNA regression resolved) +
  the suite + goldens. The complex AMBIG battery is Phase-C (AMBIG-RNA needs boundary nodes).

**Relation to v6 phases:** this corrects the RNA-imputation *content*; the bipartite **boundary nodes** (Phase
C) + **AMBIG-RNA** (which needs the per-junction motif strand) remain as planned. The single-strand nascent
imputation here is the correct Phase-B content.

---

## 6. Open questions
- **Q1 (nascent direction onto the exon):** the exon's own contained unspliced mixes mature-body + nascent; do
  we impute a nascent *prior* onto the exon at all, or is the exon fully handled by strand + spliced-floor and
  the nascent imputation is only *intron-ward*? (Leaning: the exon's RNA is strand+floor-anchored; the nascent
  imputation's main job is rescuing the intron.)
- **Q2 (precise zero-RNA at TSS/TES):** is "no imputation across" sufficient, or do we add an explicit zero-RNA
  prior on the intergenic side (it already reads f_g≈0 via the count clue — intergenic gDNA — so maybe moot)?
- **Q3 (Phase-B prior status meanwhile):** keep the committed (broken) Phase-B prior on `main` while building
  the correction, or neutralize it first for a clean baseline? (It helps gdna-present, hurts zero-DNA.)
