# Enrichment-prior redesign — gDNA-specific predictor, single-strand teacher, depleted floor

**Status:** design, 2026-06-29. Targets flagship defect (c) — the dominant ~77% of the AMBIG-exon gDNA
under-call — established by the flagship dive: on a balanced-tilt AMBIG node the strand likelihood is flat
in `f_g`, so the per-node prior MEAN is the answer, and today that mean comes from `ê(z)` which under-predicts
enriched gDNA (truth f_g 0.757, anchor 0.533). Companion to `dispersion_aware_message_precision.md` (which
shipped defect (a)/(b)'s message layer). Defects (a) one-fused-prior and (b) message re-anchoring remain
separate, after this.

## 0. The diagnosis, in one paragraph
`ê(z) = E[ρ_g | z]` is the enrichment-aware transfer that lifts a node's gDNA prior above the depleted
baseline under capture. Its quality is entirely the quality of the predictor `z`. Today `z` = boundary-crossing
flux density, which **saturates** at high enrichment (measured: z pins at ~5 while true ρ_g runs to 30–37, so
`ê` caps at ~16). The naive replacement — the node's contained density `M/E` — is **ratio-inflated**: it equals
`gDNA + RNA`, so it tracks total coverage, not gDNA. Measured: `corr(contained, gDNA)` decays 0.998→0.900 as
gDNA:RNA falls 6.8:1→0.07:1, and the deeper problem is the *level* bias (contained over-predicts gDNA by the
RNA amount, worst where RNA dominates). The fix is a **gDNA-specific predictor** (subtract the RNA), a
**non-circular teacher** (single-strand exons, where gDNA is strand-known), and a **depleted floor** (the
truly off-target minimum, not the pooled mean).

## 1. The four changes (all inside `bp_solver.fit_enrichment_transfer` + `_gdna_seed_estimate`, anchored pre-sweep)

### Change 1 — PREDICTOR: boundary-flux `z` → the depleted-residual `ρ_resid` (gDNA-specific, strand-free)
Per exon region:
```
ρ_mature = M_spliced(motif-stranded, flanking clean boundary) / E_spliced     # mature RNA density
ρ_resid  = clip( M_unspliced − ρ_mature · E_rna_contained , 0 ) / E_gdna_contained
```
- Strand-FREE (the splice motif gives the strand even on AMBIG / unstranded) ⇒ defined on every exon node.
- gDNA-SPECIFIC (the mature RNA is subtracted) ⇒ NOT ratio-inflated ⇒ the SS→AMBIG transfer is clean.
- Already exists as `rho_spliced` in `fit_enrichment_transfer` (currently used in the *response* blend, not as
  the predictor). The change is to promote it to the predictor `z`.
- Assumes nascent RNA is sparse (residual = gDNA + nascent ≈ gDNA). See open issue O1.

### Change 2 — TEACHER/RESPONSE: the clean strand-solved gDNA on single-strand exons
```
ê = MonotoneMean.fit( ρ_resid[ss_exon] , gdna_density_strand[ss_exon] , weight = (2κ−1)² )
    where gdna_density_strand = f_g_strand · M_unspliced / E_gdna  (the data-pinned strand solve)
```
- The response is now the CLEAN teacher label (strand-solved gDNA), not the old strand/spliced blend — the
  spliced-derived quantity moved to the predictor (Change 1). `ê` thus learns to **de-bias `ρ_resid` against
  the strand truth** (correcting the mature-eff-len / nascent / FL biases systematically).
- Single-strand exons are the only non-circular teacher (AMBIG gDNA is the unknown). Disjoint teacher (SS) vs
  target (AMBIG) ⇒ cross-fitting is structural (no node informs its own prior).
- Spans the spectrum: the fit ranges over all `ρ_resid` from depleted→enriched single-strand exons (Q3). The
  enriched end is pinned only by the available enriched single-strand exons (open issue O3).
- `ê`'s residual variance (`sigma2_level`) is the EBCF reliability → the per-node prior precision.

### Change 3 — FLOOR: pooled MEAN `ρ_global` → the truly-depleted minimum `ρ_floor`
```
ρ_floor = the intergenic (truly off-target) gDNA density   # the depleted baseline / "minimum"
```
- Replaces `ρ_global` (intergenic+intron+strand-seed exposure-pooled MEAN, which averages depleted+enriched and
  sits above the floor) as the SHRINKAGE TARGET / backstop. Intergenic is the purest depleted class; introns may
  be pooled in for coverage but carry nascent risk (open issue O5).
- Used as: the low anchor of `ê`, and the backstop mean where enrichment evidence is absent.

### Change 4 — APPLY + SELF-COLLAPSE (the two gates)
Per exon node, prior gDNA density:
```
ρ_prior = w_enrich · [ w_strand · ê(ρ_resid) + (1 − w_strand) · ρ_resid ] + (1 − w_enrich) · ρ_floor
```
with precision from `ê`'s residual variance (EBCF). Two derivable gates, no detector, no new constant:
- `w_strand = (2κ−1)²` — strand discriminability. κ→1: use the teacher-calibrated `ê(ρ_resid)`. κ→½
  (unstranded): the teacher degenerates, fall back to the RAW strand-free `ρ_resid` (still a gDNA estimate,
  ~0.9 corr). This is the graceful unstranded handoff (open issue O4).
- `w_enrich` = the existing slope-evidence weight. Off-capture (uniform gDNA): `ρ_resid ≈ ρ_floor` everywhere ⇒
  `ê` flat, `w_enrich`→0 ⇒ collapses to `ρ_floor` ⇒ same code, correct behaviour.
The per-node grid solve is UNCHANGED: `ρ_prior` → implied f_g prior, times the tilt likelihood. On a balanced
AMBIG node the tilt is flat ⇒ f_g = the (now-correct) prior; on a single-strand node the tilt dominates ⇒
unchanged (stays perfect). The `(2κ−1)²` switch makes the prior fill EXACTLY the tilt's null space.

## 2. Why this lifts the flagship anchor 0.533 → ~0.757
The 0.533 was `ê(saturated z)` capping at density ~16 vs truth ~30. With `ρ_resid` (monotone, gDNA-specific) as
the predictor and the strand-truth teacher, `ê` maps each AMBIG node's residual to its calibrated gDNA level
with no saturation ceiling and no RNA-ratio inflation. The transfer bias that remained under the contained
predictor (SS 0.59 vs AMBIG 0.73 at matched contained density) is removed because the residual is gDNA-specific.

## 3. Boundaries (do NOT forget — ~half the flagship leak)
The dive showed AMBIG-ADJACENT boundary nodes share the under-call (mean f_g 0.530). Boundaries carry crossing
mass, not contained mass, so they need the boundary analog of `ρ_resid` (crossing gDNA density after subtracting
the motif-spliced crossing). The plan must extend the predictor + prior to boundary nodes, not only regions.

## 4. Open issues
- **O1 (nascent contamination):** `ρ_resid = gDNA + nascent`. Sparse-nascent assumption; fails in nrna_rnd
  conditions (where the shipped message change already regressed). EBCF down-weights a noisy predictor, and the
  teacher carries nascent too (so `ê` partly learns it), but a real nascent model may be needed.
- **O2 (mature-estimate quality) — PREREQUISITE:** `ρ_resid` subtracts `ρ_mature`; a wrong mature density
  (the known spliced half-triangle eff-len ~2× issue) corrupts the predictor. Must verify/fix the mature
  density before trusting `ρ_resid`.
- **O3 (teacher spectrum coverage):** `ê`'s enriched end is pinned only by enriched single-strand exons; AMBIG
  targets above the teacher's `ρ_resid` range are extrapolations. Measure the two histograms.
- **O4 (unstranded handoff):** validate the `w_strand` raw-residual fallback improves (not regresses) unstranded;
  this is the regime the shipped change destabilized.
- **O5 (floor estimation):** intergenic-only may be low-count/noisy; introns add coverage but nascent risk.
- **O6 (validation infra):** the 0.995 was 3:1-inflated ⇒ the fix MUST be validated across gDNA levels
  (gdna25/100/400/1000). The `gdna_benchmark_5mb` BAMs were cleaned — regenerate the suite.
- **O7 (defects a/b remain):** the one-fused-prior (a) and the message re-anchoring (b) are the other ~23%,
  separate from this. Once `ê` is right, (b) — unblocking the correctly-solved single-strand neighbour message —
  compounds with it.
- **O8 (circularity/anchoring):** `ρ_resid`, `ê`, `ρ_floor` are all fit ONCE pre-sweep on belief-independent
  inputs ⇒ non-circular, consistent with the shipped anchored-global architecture. Confirm.
