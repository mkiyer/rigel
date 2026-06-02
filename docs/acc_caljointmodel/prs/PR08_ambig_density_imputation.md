# PR 8 — AMBIG-region density imputation: decodable-node ρ₀, no own-node "false gDNA"

Surfaced by the post-PR07 deep dive
([calibration_post_pr07_diagnosis.md](../calibration_post_pr07_diagnosis.md)).
**Dominant** cause of the remaining failures; first of the 3-PR robustness split
(PR8 AMBIG → PR9 strand knife-edge → PR10 spliced double-count).

## Definitions (carefully)

A calibration **region** is a maximal span of constant signature; an internal
**boundary** lies between two adjacent regions of the same reference. Mass is
deposited per region (contained) and per boundary side (crossing flux), so for the
E-step each region carries three views: *contained*, *left* (its left boundary's
near side), *right* (its right boundary's near side).

- **AMBIG region** ⟺ `ts_class == TS_AMBIG`: annotated transcripts on **both**
  strands over the span ([signature.py](../../../src/rigel/calibration/signature.py)).
  Every read is sense for one gene and antisense for the other ⇒ **no valid sense
  split** ⇒ the strand channel is invalid. (Note `TS_NONE`/intergenic is **not**
  AMBIG — gDNA is unstranded there, so an arbitrary sense is neutral and the strand
  channel stays valid.)
- **Decodable region** ⟺ `ts_class != TS_AMBIG` (i.e. `NONE`, `POS`, or `NEG`). Its
  gDNA/RNA split can be recovered (strand channel valid, or intergenic).
- **Decodable boundary** ⟺ **at least one** of its two flanking regions is
  decodable. A boundary-crossing fragment can be deconvolved (FL/count channel +
  the decodable neighbour's strand context) whenever one side is decodable — this
  is the premise of using boundary-crossing reads, and for **hybrid-capture data
  the calibration signal lives in boundaries**, so we must not discard them.
- **Undecodable** ⟺ an AMBIG region's *contained* mass (no strand, no crossing), or
  a boundary whose *both* neighbours are AMBIG. These carry no recoverable gDNA/RNA
  information.

## Root cause (audited; full detail in the diagnosis)

An AMBIG region's contained mass is undecodable, so on iteration 1 it falls back to
the uninformative `π_g = 0.5` ([calibrate.py:118](../../../src/rigel/calibration/calibrate.py#L118); strand masked, count silent) ⇒ **half its reads become
"false gDNA" with zero evidence.** That false gDNA then (a) drives the global ρ₀
(`update_rho_0` sums over **all** regions), (b) raises `pi_g_prior = ω·ρ₀·L/n_u`
which feeds back as more AMBIG gDNA, while (c) the correctly-RNA decodable regions'
ω → 0 collapses the denominator. The result is a circular `ρ₀ → prior → AMBIG-gDNA
→ ρ₀` that runs away geometrically (ρ₀ `0.2 → 20390`, never converges) and poisons
the per-locus prior so the EM dumps nascent → gDNA. The sweep compounds it by
imputing an AMBIG region's exposure partly from its *own* false gDNA
([sweep.py:163](../../../src/rigel/calibration/sweep.py#L163), `a_reg`).

## The fix

**Principle (per the design discussion):** ρ₀ **and** the exposure dispersion are
global library parameters and must be estimated from **decodable evidence only** —
decodable regions and decodable boundaries. AMBIG-undecodable contributions are
*withheld*. AMBIG regions are then *recovered* (imputed), never used to *determine*
the global parameters, and never rely on their own-node signal. "False gDNA"
disappears: an AMBIG region's contained mass no longer manufactures gDNA from a
50/50 prior — its gDNA comes only from the imputed density `μ_g = ω_imputed·ρ₀·L`.

### F1 — decodable-node `ρ₀` (essential)

Replace the all-region `Σ M_g / Σ(ω·L)` with a per-node-masked sum. Per region `r`,
with `Dr = (ts_class[r] != AMBIG)`, and the boundary side attributed to `r`
decodable iff `r` is decodable **or** its same-reference far neighbour is
(`DlB[r] = Dr or (internal-left and Dr[r-1])`, similarly `DrB`). A decodable
region's own boundary sides — including reference terminals — are therefore always
kept, so a no-AMBIG library is scored identically to before:

```
num = Σ_r [ Dr·contained.m_g[r] + DlB[r]·left.m_g[r] + DrB[r]·right.m_g[r] ]
den = Σ_r [ (Dr or DlB[r] or DrB[r]) · ω[r]·L_phys[r] ]
ρ₀  = num / den           (existing den≤0 fallback retained)
```

- The **contained** mass of an AMBIG region is dropped (undecodable).
- The **boundary** mass at a decodable boundary is kept — including the side
  attributed to an AMBIG region, because the boundary itself is decodable. This is
  the "use all decodable evidence / keep boundaries for capture" requirement.
- The denominator keeps a region's exposure capacity iff **any** of its nodes is
  decodable (so we never put boundary mass in the numerator without its capacity in
  the denominator — that would re-inflate ρ₀). A boundary with *both* neighbours
  AMBIG is fully withheld.

**Accounting note (deliberately conservative).** This keeps the existing physical
length `L_phys` denominator (so **non-AMBIG scenarios are numerically unchanged** —
every node is decodable, the sums reduce to the current formula, goldens stable).
The only approximation is that an AMBIG region with a decodable boundary keeps its
*full* `ω·L_phys` capacity in the denominator while contributing only its boundary
mass to the numerator — a small, conservative downward bias on ρ₀ (less gDNA, the
safe direction). A fully node-based denominator (contained `region_eff_len` +
boundary `μ_FL`, per [effective_length.py](../../../src/rigel/calibration/effective_length.py)) would be cleaner but changes ρ₀'s scale for *all* scenarios and
couples in the orthogonal `L_phys`-vs-FL-corrected question — deferred out of PR8.

### F1b — decodable-region exposure dispersion (essential)

`update_exposure_dispersion` is fit from the **pre-sweep** per-region exposure
posteriors. AMBIG regions' pre-sweep ω is the garbage that the sweep overwrites, so
fit the dispersion over decodable regions only (`ts_class != AMBIG`). (The
dispersion is a per-region exposure variance; only regions carry ω, so this is
region-level. It also removes the AMBIG outliers that drove the dispersion to its
ceiling in the runaway.)

### F2 — impute AMBIG exposure from neighbours only (recommended)

In `sweep_ambig_exposure`, drop the AMBIG region's **own** `(a_reg, b_reg)` from its
imputed posterior (use the propagated `from_left + from_right` neighbour evidence
only), and mask the propagation's `(a_reg, b_reg)` to **decodable** regions so one
AMBIG region's false gDNA cannot flow into another's imputation. An isolated AMBIG
region with no decodable neighbour evidence then imputes `ω → 1` (neutral), the
correct no-information fallback.

All three are **structural — no new constants** (Q6-clean): masks and dropped terms.

## Why this preserves boundary evidence (capture)

For hybrid-capture libraries the gDNA enrichment is read out at captured-exon
boundaries. F1 keeps **every decodable boundary's** crossing gDNA in ρ₀ — including
boundaries that flank an AMBIG region, as long as the *other* side is decodable. We
withhold only the genuinely uninformative contained-AMBIG mass and both-AMBIG
boundaries. Validated on `anti_intron_nrna_30`: 11 regions, 1 AMBIG (region 5);
region 5's contained mass is withheld but both its boundaries are decodable (POS
neighbours), so its boundary evidence is retained.

## Test plan — empirically prove correctness across AMBIG topologies

Define and sweep a matrix of AMBIG scenarios (gDNA ∈ {0,20,100}, ss ∈ {0.65,0.9,1.0}):

1. **Contained AMBIG** (antisense single-exon gene inside another's intron) — the
   `anti_intron` failures. Expect: zero-gDNA ⇒ ρ₀≈0, converged, nRNA detected,
   gDNA≈0; real-gDNA ⇒ AMBIG region still receives correct gDNA.
2. **Overlapping AMBIG** (two genes overlapping on opposite strands) —
   `overlap_antisense`. Boundaries span the overlap; confirm RNA recovered.
3. **Both-neighbours-AMBIG region** (undecodable boundaries) — a constructed case;
   confirm graceful imputation (ω→1) and no runaway.
4. **AMBIG + real gDNA capture-like** (enriched boundary, depleted intergenic) —
   confirm boundary gDNA anchors ρ₀ and AMBIG gDNA is recovered, not zeroed.
5. **Regression:** single_exon / non_overlapping / paralogs (no AMBIG) numerically
   unchanged; zero-gDNA still recovers; goldens for AMBIG scenarios regenerate.

Each must show `converged=True` and ρ₀ of the right order. The matrix is the proof.

## Rollback

F1, F1b, F2 independently revertable.

## Open questions resolved (this round)

- *Dispersion exclusion?* **Yes** (F1b) — decodable regions/boundaries feed both
  density and dispersion; undecodable nodes are withheld from both.
- *Eliminate boundaries for simplicity?* **No** — per-node decodability keeps
  decodable boundaries (≥1 decodable neighbour); capture signal lives there.
- *"False gDNA"?* **Removed** — AMBIG regions rely on sweep + global ρ₀, not
  own-node 50/50.

## Resolution + validation

Denominator: **conservative `L_phys`** (decision confirmed). Validated on the
scenario suite + goldens:

- The three runaway cases converge with the correct ρ₀: `anti_intron nrna30/70`
  (ρ₀ `20386 → ~2e-7`, RNA 2000/2000, gDNA 0, nRNA detected) and `overlap_anti fc4`
  (ρ₀ `6.2 → ~3e-8`, RNA fully recovered). `test_antisense_intronic` and
  `test_overlapping_antisense` now pass.
- **Non-AMBIG scenarios numerically unchanged** — confirmed: only the two
  AMBIG-bearing goldens moved (`antisense_overlap_ss90` gDNA `434.9 → 5.9e-6`;
  `antisense_contained_ss90`), every other golden byte-identical, and `nrna_dc`
  ρ₀ matches pre-PR8 exactly.
- The remaining 3 `nrna_dc g20` failures are Mechanism 2 (real-gDNA strand
  weakness), the target of PR9 — out of scope here.

The fully node-based `region_eff_len`/`μ_FL` denominator (which would also require
reworking the exposure posterior's `β_post`) is deferred to its own effective-length
PR.
