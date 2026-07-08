# Three-channel calibration: mature RNA, nascent RNA, and gDNA

**Status:** design v3, for review (2026-07-07). Full rewrite incorporating two rounds of maintainer review + a
third-party critique. Supersedes `disagreement_shrinkage_prior_design_v2.md` (the two-component variance that
regressed) and folds in `mature_rna_channel_design.md`. Grounded on `~/Downloads/rigel_runs/ambig_dense_10mb`.
Memories: `three_component_mature_nascent_design`, `mature_channel_proof_and_design`,
`rna_imputation_transcript_structure`, `region236_ambig_tau_nullspace_intron_ceiling`.

> **TL;DR — one framework, one new field.**
> RNA travels in **two channels**, gDNA in a third:
> * **Mature** is the *default* RNA channel. It **is today's RNA message** — pass-through plus a spliced
>   source/sink at junctions — merely **re-gated to flow only between mature-eligible (exon) nodes.** It never
>   enters an intron.
> * **Nascent** is the one genuinely new field. It **originates only inside introns**, flows contiguously in
>   its own channel until a TSS/TES, and **defaults to zero**. It can never become mature (the one-way rule).
> * **gDNA** is unchanged (genomic).
>
> The node belief stays the 3-term simplex `(f_pos, f_neg, f_g)`; the grid solver, the calibration output, the
> EM, and the priors are all **unchanged**. The split lives entirely in the message layer. This fixes the v2
> regression at the source: in a `nrna_none` condition the introns hold no nascent → the nascent channel is 0
> → an exon's mature stays in the mature channel instead of bleeding into the intron as fake nascent.

---

## 1. The framework

Every RNA molecule is either **mature** (spliced product, exon-only) or **nascent** (unspliced pre-mRNA, whole
gene body). These are physically distinct populations, so they get distinct message channels with distinct flow
rules. The design is the set of rules for the two channels:

| | **mature channel** | **nascent channel** |
|---|---|---|
| what it carries | spliced product (default RNA) | pre-mRNA |
| **originates** | everywhere (the default) — an exon's own RNA is mature unless nascent claims it | **only inside introns** (mature-ineligible RNA-bearing regions) |
| **flows** between | mature-eligible (exon) nodes | any strand-continuous nodes (exon↔intron↔exon) |
| **splice junctions** | **source/sink** (the observed spliced mass), strand-locked by motif | pass through (a pre-mRNA read spans the junction unspliced) |
| **exon↔exon boundaries** | pass through — ordinary BP propagation, *not* a source/sink | pass through |
| **introns** | **cannot enter** (an intron has no mature) | present (this is where it lives) |
| **TSS / TES** | sink | sink |
| default value | the residual: `total RNA − nascent` | **0** |
| ↔ opposite strand | **never** (mature is motif-stranded, strand-locked) | never |

The **one-way rule** (maintainer): once spliced, RNA never becomes unspliced. Nascent can turn into mature
(splicing) but never the reverse, and in the calibration model neither channel's mass ever crosses to the
opposite strand.

**Why this is the whole fix.** The v2 regression was a *lumped* RNA message carrying an exon's mature into
introns as fake nascent. Here mature simply cannot enter an intron (it is not mature-eligible), and nascent is
zero unless an intron sourced it — so in `nrna_none`, nothing flows into the introns. The gDNA/RNA solve and
the calibration output are untouched; only *where each RNA species is allowed to travel* changes.

## 2. Definitions — crisp at boundaries, ambiguous only in exons

* **Spliced ⇒ mature, always, and always stranded by the splice motif** (independent of library
  strand-specificity, `scoring.py:46`). Not all mature is spliced — a mature molecule's exon-body reads are
  unspliced and look identical to a nascent read.
* **The two species are perfectly distinguished at a boundary:** the *spliced* crossing mass is mature; the
  *unspliced* crossing mass is nascent. This is observed, and it is where the mature source/sink lives.
* **Ambiguous only inside an exon:** an exon-body unspliced read could be mature or nascent — the strand tilt
  cannot tell them apart. So the split is a message-layer quantity, never a per-node likelihood term (§3).
* **By region type, per strand:** intergenic = gDNA; intron = nascent | gDNA (no mature); exon = mature |
  nascent | gDNA.

## 3. Node belief — the 3-term simplex is unchanged (agreed)

The belief stays `(f_pos, f_neg, f_g)` over the unspliced mass, solved by the unchanged per-node grid
(`simplex_logodds`). The strand tilt — the only intrinsic signal — identifies only the totals `(R₊, R₋, g)` and
is exactly flat along the mature/nascent split, so a joint 4-simplex would grid a locally-flat axis and blow up
AMBIG nodes to a 4-D cube (the retired-lattice OOM). The five-term latent `{m₊,m₋,n₊,n₋,g}` **factors** into the
strand-identified simplex × a per-strand mature/nascent split that lives entirely in the **messages**. The
solver receives one combined RNA density per strand (mature + nascent), exactly as it does one RNA density
today.

## 4. Node structural properties — derived at index build, read as flags (Phase 1)

The routing must not be signature bit-twiddling inside the `_scan` inner loop — it is unreadable, untestable,
and *wrong* at overlap seams (the signature is only 4 bits and cannot express transcript multiplicity). **All
topology is derived from the annotation at index build, stored in the index, wired through to calibration, and
read by `_scan` as flat boolean flags.** This is a standalone, unit-tested phase with no behavior change.
Decision (was open issue O1): **annotation at index build** — count-zero-info-correct (structure never depends
on counts) and robust at the retained-intron/AMBIG seams that defeat both the signature bit-pattern and the
read-observed `junction_strand`.

**Per region — stored in `region_df` (new columns), computed from the transcripts:**
* `mature_eligible[s]` — **the region overlaps an exon of a *multi-exon* (spliced) transcript on strand `s`.**
  This is the elegant statement of "spliced ⇒ mature": every hard case falls out — a single-exon transcript's
  exon is not multi-exon ⇒ nascent (`mature ≡ nascent`, the maintainer's rule); a retained-intron region
  (intron of a spliced tx + a single-exon tx) is not multi-exon-covered ⇒ nascent; the three-transcript
  `(5000,9000)` is covered by a multi-exon exon ⇒ mature-eligible. **The 4-bit signature cannot express this**
  (it conflates single- and multi-exon exons — the same blind spot that hides the retained-intron junction).
  * **Overlap, not full coverage** — a region has constant *signature* but its multi-exon-vs-single-exon
    coverage can vary within it (a multi-exon terminal exon merges with an adjacent single-exon exon when no
    intron separates them). "Overlaps a multi-exon exon" over-approximates safely: it never *under*-calls
    mature-eligibility, so it can never re-introduce a bleed (a mislabeled single-exon fragment travels the
    mature channel, which cannot enter an intron — benign; the aggregate EM disentangles it downstream).
* `nascent_eligible[s]` — the region carries strand-`s` pre-mRNA = it has *any* strand-`s` signature bit. This
  is derived from the (already-stored) signature at calibration setup — **not** a new stored column.

**Per boundary — stored in a new `boundaries.feather`** (one row per boundary, ref-major, `k+1` per `k`-region
reference, aligning by construction with the accumulator's boundary indexing; validated by count):
* `genomic_sj_strand` — the splice-junction motif strand (`0`/`POS`/`NEG`), set where a boundary position is a
  transcript intron endpoint. **Used** to route the mature (spliced) mass to the correct strand's exon flank
  from the annotation (replacing the read-observed `junction_strand`, a follow-up C++ cleanup).
* `is_splice_junction` — `genomic_sj_strand != 0`. **Used** as the mature source/sink topology marker.
* `is_tss` / `is_tes` — a transcript's 5′/3′ terminus at this position (orientation by strand). **Used** in the
  calibration diagnostic output (transcript-boundary QC); *not* consumed by the flow gates (continuity is
  `free_s`, which correctly handles a nested TSS inside another transcript's body).

**The emission gates then read as flags** (per strand `s`, edge `src → dst`):
```python
emit_gdna       = facing_mass > 0                                             # genomic, unchanged
emit_nascent[s] = free_s(edge) and src.nascent_eligible[s] and dst.nascent_eligible[s]
emit_mature[s]  = free_s(edge) and src.mature_eligible[s]  and dst.mature_eligible[s]
```
`emit_nascent` reduces to the existing `free_s` topology (continuity ⇒ both flanks strand-present ⇒ both
nascent-eligible). `emit_mature` is the new restriction: mature flows only mature-eligible ↔ mature-eligible.

> **Correction to the third-party review.** Its gate `emit_mature = … and (not is_splice_junction)` is **wrong**
> — verified against the three-transcript example (§9). A splice junction does **not** block the mature channel;
> mature must still flow across it into a mature-eligible destination (TB's mature continuing through TA's
> donor). The junction only source/sinks the spliced mass. The gate keys on **`mature_eligible` (multi-exon
> coverage)**, never on junction-presence.

## 5. The two RNA channels in `_scan`

Alongside the unchanged gDNA message, the single per-strand RNA message becomes two. Everything is per strand
`s` and stays on strand `s`.

**Mature channel `ρ_mat,s` — today's RNA message, re-gated.** Source density = `fb_mat·M_src/E_rna` (the
running mature belief) `+ SPs/E_spl` (junction source) `− SPd/E_spl` (junction sink), i.e. the *current*
`_scan` RNA formula (`bp_solver.py:1216-1220`) with `fb_mat` in place of the lumped `fbp`. Emitted only where
`emit_mature[s]`. The spliced source/sink terms are non-zero only at junctions (the geometry routes spliced to
the exon flank), so they fire automatically and correctly — including the retained-intron seam.

**Nascent channel `ρ_nasc,s` — new.** Source density = `fb_nasc·M_src/E_rna` (no spliced terms — nascent has no
junction source/sink). Emitted where `emit_nascent[s]`. Sourced only inside introns (§6); defaults to 0.

**Into the solver:** the RNA⁺ message is the combined density `ρ_mat,s + ρ_nasc,s`, one Gaussian on `log f_pos`
— the solver is untouched. An intron destination receives `ρ_mat = 0` (mature can't enter), so its RNA is
nascent-only; an exon receives both. Mature therefore still lowers `f_g` (RNA presence) exactly where it is
allowed to exist.

**Removed:** the current `rho_mat_dst` *subtraction hack* (`bp_solver.py:1218`) is gone — mature no longer needs
to be subtracted out of a lumped message going into an intron, because the mature channel simply does not emit
into an intron.

## 6. The node split — mature is the default, nascent is the residual's complement

Each node carries two running beliefs per strand, `fb_mat` and `fb_nasc`, with the invariant
`fb_mat + fb_nasc = f_pos` (the solved total RNA fraction). The split is set by structure + the nascent channel:

* **Nascent source — mature-ineligible RNA-bearing regions (pure introns, single-exon exons):** local
  `fb_nasc = f_pos` (all nascent), `fb_mat = 0`. This is where nascent is *conceived*.
* **Mature default — mature-eligible regions:** local `fb_nasc = 0`, `fb_mat = f_pos` (all mature by default).
  Nascent reaches them **only** through the nascent channel from adjacent introns; their RNA *excess* over that
  arriving nascent level stays mature (`fb_mat = f_pos − fb_nasc`).
* **TSS/TES:** source neither channel (both sink). A terminal exon's RNA is mature by default and gets nascent
  only if an interior intron drives it up.

This is the maintainer's conservative rule made mechanical: **`nascent = 0` unless an intron lifts it.** No node
stores a mature/nascent *fraction* to solve — mature is the residual of the total after the nascent channel.

## 7. The bleed, resolved — the worked example

`TA⁺` exons `(1000,5000),(10000,14000)`; first exon density 10; junction at 5000, spliced density 2.

**Buggy behavior (removed):** emit `total − spliced = 10 − 2 = 8` as nascent into the intron. The spliced count
is capture-biased and *under-counts* mature, so the 8 is fake nascent — the regression.

**This design:**
* `nrna_none`: the intron `(5000,10000)` is a pure intron with RNA ≈ 0 → the nascent channel is **0** there →
  nothing to emit into the intron. The exon is mature-eligible → its density-10 RNA is mature (`fb_mat = 10`,
  `fb_nasc = 0`) and travels the mature channel, which **cannot enter the intron**. **Emit 0 nascent → no
  bleed.** The intron stays `f_g ≈ 1`.
* **nascent present:** the intron carries real nascent → the nascent channel lifts it into the exon; the exon's
  excess is mature. Correct, same machinery.

So the answer to "do we emit 8 as nascent?" is **no** — we emit the nascent channel value (0 here); mature never
enters the intron.

## 8. Initialization / direction-dependence — resolved (former open Q3)

The knot was: a TSS sources total RNA with no mature/nascent split, and a forward pass cannot partition it until
it reaches a junction. It dissolves: **the TSS/TES source neither channel — there is no unsplit-total message to
partition.** Nascent is sourced *only inside introns* and defaults to 0; mature is the default everywhere else.
A terminal exon receives mature from its junction and nascent from its intron, both arriving from the interior
via whichever forward/backward pass carries them; the forward pass *from the TSS* carries zero of both,
correctly. The forward/backward combine (`α ⊗ β`, `bp_solver.py:1259-1270`) integrates the two directions
per channel, as it already does per component. Because nascent defaults to 0 across the *entire* sweep unless an
intron with real counts lifts it, the forward and backward passes cannot disagree on a phantom baseline.

## 9. Overlapping transcripts & intron retention — handled by the flags, no special cases

The stress example: `TA⁺(1000-5000 ^ 20000-24000)` and `TB⁺(1000-9000 ^ 16000-24000)`, junctions at 5000, 9000,
16000, 20000. Region eligibility (verified 2026-07-07):

| region | overlap | `nascent_elig⁺` | `mature_elig⁺` | role |
|---|---|---|---|---|
| (1000,5000) | EX(TA)+EX(TB) | ✓ | ✓ (exon-flank of 5000) | mature default |
| (5000,9000) | **IN(TA)+EX(TB)** | ✓ | ✓ (exon-flank of **9000**, TB donor) | mature; nascent via channel |
| (9000,16000) | IN(TA)+IN(TB) | ✓ | ✗ (intron flank only) | **nascent source** |
| (16000,20000) | **IN(TA)+EX(TB)** | ✓ | ✓ (exon-flank of 16000) | mature; nascent via channel |
| (20000,24000) | EX(TA)+EX(TB) | ✓ | ✓ (exon-flank of 20000) | mature default |

* **Boundary 5000** `EX⁺ ↔ (EX⁺|IN⁺)`: destination mature-eligible ⇒ mature **flows** (TB continues) while the
  spliced mass source/sinks TA's donor. Both behaviors, one gate. "Mature can splice *or* continue."
* **Boundary 9000** `(EX⁺|IN⁺) ↔ IN⁺`: destination is a pure intron ⇒ mature-**ineligible** ⇒ mature **sinks**
  entirely (TB's donor), only nascent survives into `(9000,16000)`.
* `(5000,9000)` and `(16000,20000)` are mature-eligible (TB exons) **and** intron for TA. They are *not* nascent
  sources (they are mature-eligible); TA's nascent reaches them through the nascent channel from the pure intron
  `(9000,16000)`, and their excess over that nascent level is TB's mature. The mixed `EX|IN` region falls out of
  the same flags with **no branching** in the solver loop.

Only the one pure-intron region `(9000,16000)` conceives nascent for this locus; everything else is mature by
default, lifted by the channel where TA's pre-mRNA is real.

## 10. Sourcing nascent — the integrated strand + density likelihood

Nascent is *conceived* only at **intron regions** and **exon-intron boundaries** (never exons — that was the
bleed). Its source must be the node's **message-free local solve** (like the strand solve); messages only
*relay* it. Sourcing it from the post-message `f_pos` was the Phase-2 over-source bug (an unstranded gDNA
intron's ambiguous ~0.25 `f_pos` leaked as phantom nascent). The fix: two **likelihoods** that resolve the
node's own gDNA-vs-nascent split before any message.

### 10.1 The regime map

Strand is the one **universal, reliable** signal; density is a **regime-dependent** secondary whose *control*
changes:

| regime | intron region | exon-intron boundary |
|---|---|---|
| **stranded** (± capture) | strand tilt → nascent ✓ | strand tilt → nascent ✓ |
| **unstranded, no capture** | density vs **intergenic** background (clean, depleted control) ✓ | (density confounded) |
| **unstranded + capture** | introns depleted → density weak ⚠ | **only** source, but confounded ⚠ (deferred) |

### 10.2 Two likelihoods that integrate themselves (no explicit weight)

Both are log-likelihoods on `log f_g`, summed in `ψ` (the solver already does this). The **strand** Beta-Binomial
carries Fisher information `∝ N·(2κ−1)²` — it **vanishes at κ=½** (gDNA-mean = RNA-mean = ½), so it contributes
exactly nothing to an unstranded node. The **density-vs-background** term is always present. So the `(2κ−1)²`
weighting is *built into the strand likelihood's Fisher information* — we never write it explicitly. Stranded ⇒
strand dominates; unstranded ⇒ density dominates. Automatic.

### 10.3 The intron density likelihood (intergenic control)

Introns and intergenic are **both off-target/depleted** under capture, so the intergenic background is a
*clean, un-confounded* control for intron density (the enrichment confound lives at boundaries, near probes —
§10.5). Model the intergenic background as a log-density Gaussian `log ρ_bg ~ N(μ_bg, σ²_bg)` (a distribution,
absorbing capture non-uniformity). For a mature-ineligible **region** the density term on `log f_g` is:

```
target    = log( ρ_bg · E / M )          # pin the gDNA component to the background FRACTION ρ_bg/ρ
precision = 1 / ( σ²_bg + 1/N )          # data-driven: background spread ⊕ the node's Poisson count info
```

* intron at background density → `target ≈ 0` → `f_g≈1`, **no nascent**;
* intron at 10× background → `target ≈ log 0.1` → excess is nascent;
* between → smooth (no threshold). "How likely by chance" **is** the precision — the inverse variance of that
  Gaussian. Deep-in-bulk pins `f_g≈1` tightly; far-tail pins `f_g` low tightly; ambiguous → low precision.

This **replaces the capped floor** for mature-ineligible regions (below); it is present in **pass-1**, so the
nascent source is calibrated before any message.

### 10.4 Precision is bounded — "uncapped" ≠ unbounded (resolves the cap concern)

The 1-pseudo-observation cap (`_GLOBAL_STAB_PREC`) exists so a **prior** cannot override a node's own data. The
density term is a **likelihood** (the node's count vs the population), so it is *not* subject to that cap — but
it is **hard-bounded by the data**: `1/(σ²_bg + 1/N) ≤ N` (never exceeds the node's own Poisson information) and
`≤ 1/σ²_bg` (a finite empirical spread). It **cannot collapse** like the retired `vbg→0` runaway (both `σ²_bg`
and `N` are fixed — nothing sharpens). And it only sets the node's **own** belief: the **outgoing message**
precision remains ceilinged by the disagreement variance `pr = n/(n·σ²_imp+1) ≤ 1/σ²_imp` (§11), so no
overconfident message can ever leave a node regardless of its local certainty.

### 10.5 Boundaries — strand-only first (density nudge deferred)

At an exon-intron boundary the unspliced crossing is nascent+gDNA, and under capture an enriched boundary pulls
down captured-gDNA **and** captured-nascent — indistinguishable **by density** (the intergenic control is wrong
here: the boundary is enriched, intergenic depleted). So **strand is the only reliable boundary signal for
v1**, and the nascent source at a mature-ineligible **boundary** is **gated by the strand discriminability**:

```
precision(nascent source | boundary) = (2κ−1)² · p_loc
```

Unstranded (κ=½) ⇒ precision 0 ⇒ the boundary **self-sources nothing** (it becomes a pure relay of nascent
arriving from adjacent introns) — no phantom nascent from an enriched boundary. Stranded ⇒ the boundary sources
via strand. The `(2κ−1)²` here is the derived strand discriminability (not a tuned constant); it is the interim
stand-in for the deferred boundary density nudge.

*Deferred boundary density nudge (O5):* a self-referential **exon-intron-boundary density distribution** (built
in pass-2, after the KDE, since boundaries are capture-enriched). Its broad spread — depleted-off-probe ⊕
enriched-on-probe — **automatically gives the density term low precision** via `1/(σ²_ctrl + 1/N)` (no threshold,
no control-region set needed: the confound expresses itself as a wide σ²_ctrl → weak nudge → strand dominates).
Under the sparse-nascent assumption the distribution ≈ the gDNA-boundary control. Build only if the A/B shows
real nascent left on the table under unstranded capture; the risk is false-positive nascent on strongly-captured
gDNA boundaries, so strand-first is the conservative default.

### 10.6 count-zero-information — reconciled

A count *in isolation* carries no gDNA/RNA information (the principle). But a count *relative to the population
background* is the **global-prior source** — the third architectural signal. The density term is that source,
made data-driven; the only departure is that for unstranded introns it is allowed to **dominate** (there is no
intrinsic strand signal to override). This is exactly the maintainer's "differential density becomes the
dominant signal when strand is weak."

### 10.7 Strand assignment of intronic nascent (unchanged)

Single-strand intron → assign to that strand. AMBIG unstranded intron → the strand is genuinely
under-determined; the existing τ-marginal splits it, resolved downstream by the EM. We do not invent a strand
model for that case — the conservative default (nascent 0) is why it stays safe.

## 11. Per-channel disagreement variance + diagnostics

The three-channel extension of the shipped v1/v2 estimator (`_poisson_moment_var`,
`adjacent_component_disagreement_variance`) yields `σ²_imp` for gDNA, nascent, and mature over adjacent edges,
staged as v2 (pass-1 scalar → refit per-channel on the solved belief, AMBIG excluded → pass-2 → refit incl.
AMBIG). **We do not assume `σ²_mature ≫ σ²_nascent`** — the maintainer is right that the large-variance signal
is a capture **depleted↔enriched crossing** phenomenon afflicting *every* channel across a capture boundary, so
we **measure** all three and let the A/B say whether mature separates. Whatever the values, the algebraic
`pr = n/(n·σ² + 1)` caps every message precision, killing the region-503 override structurally without a bespoke
strand-gate.

**Diagnostics (maintainer's #7):** surface the three `σ²_imp`, the message-precision distribution, and the
nascent/mature routing-fraction distributions in the standard calibration diagnostic output — a genuinely
useful calibration-health readout. (A summary-diagnostic bug fix just landed on `main`; wire the new lines into
that same path.)

## 12. Scope — message layer only, output unchanged (confirmed)

Nothing downstream reads a mature/nascent split (subagent-confirmed): `CalibrationResult` stores only
`mass_gdna_*`/`mass_rna_*`; `assemble_priors` collapses RNA to one `rna_prior_count`; the C++ EM prior is two
aggregate scalars and reconstructs mature-vs-nascent itself from reads. The region/boundary projections still
emit `gdna_mass`/`rna_mass` (RNA re-summed over the two channels), so there are **zero** changes to `result.py`,
`priors.py`, `estimator.py`, `em_solver.cpp`, or the grid solver `simplex_logodds.py`. Blast radius: the
structural-flags pass (new) + `bp_solver._scan` + the estimator + `calibrate.py` refit — behind the existing
`sweep_disagreement_shrinkage` flag.

## 13. Open issues / decisions

* **O1 — junction-flag source. RESOLVED → annotation at index build** (count-zero-info-correct, robust at the
  retained-intron/AMBIG seams). Shipped in Phase 1.
* **O2 — boundary density nudge. DEFERRED (maintainer): strand-only at boundaries for v1.** The self-referential
  boundary-density distribution (§10.5) is added only if the A/B shows real nascent left on the table under
  unstranded capture.
* **O3 — local-exon control for boundaries. DEFERRED (maintainer).** Using the adjacent exon's gDNA as the
  boundary's gDNA baseline is really an *imputation* (nascent bleeding from a gDNA neighbour) — we already have
  imputation infrastructure, but it removes the concept of a node as a *self-defining* nascent source, so we do
  not use it.
* **O4 — exon-intergenic (TSS/TES) boundaries as control. HELD OUT (agreed).** Too rare + degradation-biased
  under capture.
* **O5 — double-counting at the solve.** The combined RNA message sums mature + nascent densities from disjoint
  neighbour sets; verify on a node flanked by *both* an intron and an exon that `ρ_mat + ρ_nasc` does not
  over-impute the total.
* **O6 — AMBIG-unstranded nascent strand (§10.7)** is under-determined; equal-split-then-EM is a fallback, not a
  model. Watch.
* **O7 — AMBIG τ (secondary, watch).** Motif-stranded mature may help break the region-236 τ over-hedge; measure,
  don't engineer.
* **O8 — boundary `f_g` under capture.** The genome-wide density prior over-attributes RNA to *enriched*
  boundaries (wrong control) — a pre-existing effect the deferred boundary density nudge (O2) addresses. The
  Phase-2c `(2κ−1)²` source-gate prevents this from leaking into the *nascent source*; the boundary's own `f_g`
  output is out of scope until O2.

## 14. Implementation plan (phased; behind `sweep_disagreement_shrinkage`)

**Phase 1 — structural properties at index build → calibration. ✅ SHIPPED (green, 1082 tests).**
Per-region `mature_eligible_{pos,neg}` (multi-exon overlap) in `region_df`; new `boundaries.feather`
(`is_tss`/`is_tes`/`is_splice_junction`/`genomic_sj_strand`); load + validate + schema guard;
`RegionArrays.mature_eligible` + new `BoundaryArrays`; `INDEX_FORMAT_VERSION` 5→6; `test_structural_masks.py`.

**Phase 2 — the two RNA channels in `_scan`. ✅ IMPLEMENTED (legacy flag-off byte-identical; routing unit-tested).**
`_chain_mature_eligible`; `nasc_p`/`nasc_n` running beliefs; mature portion (`fbp−nasc_p` + spliced source `SPs`)
gated by **destination** `mature_eligible`; `rho_mat_dst` removed. A/B: total error improves, beats v2 on the
intron bleed (0.758→0.806) — but a **residual bleed remains** (0.993→0.806) because the nascent *source* is the
naive post-message `f_pos`. Phase 2c fixes that.

**Phase 2c — the calibrated nascent source (THIS PLAN; the fix for the residual bleed).** The nascent source
becomes the node's message-free strand+density solve (§10), not the naive `f_pos`.

1. **Intergenic background density model** (`regions`/`bp_solver`): compute `(μ_bg = ⟨log ρ_g⟩, σ²_bg =
   Var(log ρ_g))` from **intergenic** regions' gDNA density (clean control). Refine `_floor_estimate` /
   `_gdna_seed_estimate` to expose `μ_bg, σ²_bg` (intergenic-only for the shape; the pooled rate may keep intron
   gDNA). No new tuned constant.
2. **Data-driven intron density likelihood** (`_global_logprior`): for mature-ineligible **region** nodes,
   set the density term `target = log(ρ_bg·E/M)`, `precision = 1/(σ²_bg + 1/N)` (N = node unspliced count),
   **not clamped to `_GLOBAL_STAB_PREC`** (it is a likelihood, bounded by `N` — §10.4). Threads a
   mature-ineligible-region mask + `N` + `(μ_bg,σ²_bg)` into `_global_logprior`; present in **pass-1**. This
   calibrates introns' `f_g` (= their nascent) so a background intron pins to `f_g≈1`.
3. **Nascent-source precision by node type** (`node_sweep`). The intergenic control is valid only for
   **off-target (depleted)** nodes, so:
   * **pure-intron regions** (coarse-type intron — off-target, intergenic control valid): `pnasc = p_loc`
     (density-calibrated by step 2);
   * **all other mature-ineligible sources — exon-bearing single-exon/retained-intron regions *and*
     exon-intron boundaries** (on-target/capture-enriched → intergenic control confounded): `pnasc =
     (2κ−1)²·p_loc` (strand-only; density deferred, O2). Unstranded ⇒ precision 0 ⇒ self-source nothing (pure
     relay), so an enriched node cannot inject phantom nascent;
   * **mature-eligible nodes**: `pnasc = 0` (relay).
   `(2κ−1)²` is the derived strand discriminability (`rna_sense_frac`/`gdna_strand`), not a tuned constant.
4. **Validate** (`ab_bleed.py` + gene-level): `nrna_none/*/on` introns back to `f_g≈0.99`; `nrna_present`
   nascent still flows; stranded introns unchanged; the message ceiling `1/σ²_imp` is unchanged (local-only
   change). Legacy flag-off stays byte-identical.

**Phase 3 — per-channel variance, diagnostics, dead-code.** Extend the estimator to
`(σ²_gDNA, σ²_nascent, σ²_mature)`; the RNA message precision uses the channel(s) crossing each edge; refit
staged as v2. Wire the three `σ²_imp` + routing fractions into the diagnostic output. **Strip the dead spliced
lower-bound floor** (`simplex_logodds.py:191-194` + the `build_node_statics` zeroing, maintainer's #6).

**No changes** to `result.py`, `priors.py`, `estimator.py`, `em_solver.cpp`, or the grid solver.

## 15. Validation

In-process A/B (`total_threads=1`) on `ambig_dense_10mb`. Targets/guards:
* **Regression gone:** truly-gDNA introns in `nrna_none/*/on` stay `f_g ≈ 0.99`; gene-level mature
  `spearman_abund` on `gdna100/ss0.50/on` ≥ baseline 0.873 (beating v1's 0.892).
* **Region 503:** `f_g ≈ 0.63`, mature message precision bounded (`≈ 1/σ²_mature`, measured).
* **Overlap/retention (§9):** the `EX|IN` regions stay RNA (not under-called to gDNA); the pure intron correctly
  carries nascent (present) or gDNA (`nrna_none`).
* **nascent present:** real nascent still flows into introns (no over-correction).
* **Capture-OFF / no-gDNA / stranded:** flat; the three `σ²_imp` reported and sane; EX+IN− and mature-AMBIG
  classes recover (v1 regressed them).
* Full 24-condition net-flow + mature accuracy: flat-or-better. Ship one phase at a time on a clean result.
```
