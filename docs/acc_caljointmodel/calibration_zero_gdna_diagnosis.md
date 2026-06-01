# Zero-gDNA edge case — why exposure weights don't shrink to 1.0 (root-cause diagnosis)

**Context.** Exposure weights ω exist for **hybrid capture**: probes enrich some
regions ~1000× and deplete others, so gDNA distributes very unevenly. ω = 1.0 is
neutral; ω > 1 enriched, ω < 1 depleted. `gdna_eff_len = Σ φ·ω·L_phys`
(exposure-weighted) is therefore **correct and necessary** — without it, a 2 kb
captured region inside a 100 kb locus would be drowned by 98 kb of depleted
space. (An earlier "make `gdna_eff_len` geometric" idea was wrong: it fixes the
symptom below but destroys capture compensation and regresses gDNA-bearing
scenarios 12→18.)

**Symptom.** With **zero real gDNA** (`single_exon` baseline: ss=1.0, no gDNA, a
spliced helper gene), the single-exon gene t1 loses 298/339 reads to gDNA
(87.9% error). The exposure weights should shrink to ~1.0 (no signal ⇒ uniform),
but t1's region gets ω=0.039.

## The chain (traced — `scripts/debug/trace_zero_gdna.py`)

The exposure posterior is a **Gamma shrinkage** toward the prior mean 1.0:

```
ω = (1/φ + M_g) / (1/φ + ρ₀·L)        # shrinks to 1.0 ⟺ prior pseudocount 1/φ ≫ ρ₀·L
```

It *would* give ω→1 under true zero gDNA (ρ₀→0). It doesn't, because of a cascade
seeded by **phantom gDNA**:

1. **`eps_s = 1.0`** (the gDNA splice-artifact rate). `update_eps_s` =
   `(1 + Σπ_g·n_s) / (1 + Σn_s)` is evaluated on the **contained** view, but in
   single_exon **all 39 spliced reads are junction-crossing → boundary**, so
   contained `Σn_s = 0` ⇒ `(1+0)/(1+0) = 1`. (Even with data this prior is
   mis-centered: the doc calls it "Beta(1,1)" but the formula's no-data value is
   1, not 0.5 — and biologically a spliced read is RNA, so the rate should be
   ≈0, not 1.)
2. **eps_s=1 ⇒ all 39 junction reads → gDNA** (`m_g += eps_s·m_spliced`).
   `total_M_g` jumps from ~0 to 39 (= exactly `n_spliced_obs`). The *contained*
   reads are still correctly all-RNA (`pi_g=0`); the phantom gDNA is entirely the
   boundary/spliced mass.
3. **Phantom gDNA ⇒ ρ₀ ≠ 0** (0.00077) and **spurious exposure heterogeneity**:
   t_helper's exon regions (which carry the junction reads) get ω≈52; t1's region
   gets ω≈0.04.
4. **φ M-step fits φ=39.77** from that spurious ω spread.
5. **Large φ ⇒ weak shrinkage prior** (1/φ=0.025 ≪ ρ₀·L=0.77) ⇒ t1's ω stays at
   0.039 instead of shrinking to 1.
6. **Low ω ⇒ small `gdna_eff_len`** (ω·L≈39) ⇒ in the per-locus EM a tiny eff-len
   = high per-fragment density ⇒ the gDNA component becomes a fragment
   **attractor** ⇒ t1's reads → gDNA (the 87.9% error).

## Confirmation

Forcing `eps_s = 0.01` (spliced reads are RNA): `total_M_g` 39→0.39, ρ₀→5e-5,
φ→2.9, and **t1's ω 0.039→0.854** (shrinks toward 1; `gdna_eff_len≈850 ≈ t1's
~800`, no attractor). The exposure shrinkage works once the phantom gDNA is gone.

## Root cause & fix direction (to discuss — Q6)

**Taproot: the `eps_s` mechanism is degenerate.** It (a) defaults to 1.0 with no
data, (b) is mis-centered (a spliced read carries a junction ⇒ it is RNA; the
gDNA splice-artifact rate is biologically ≈0), and (c) is estimated from
*contained* spliced reads only, which are ~absent whenever splicing implies
boundary crossing. Candidate fixes (need your call on the parameterization):
- center the `eps_s` prior **low** (spliced ⇒ RNA), e.g. a Beta with small mean,
  so the no-data / sparse-data default is ≈0 rather than 1;
- estimate `eps_s` from **all** spliced reads (incl. boundary), not contained only.

**Secondary (related, optional):** even with `eps_s` fixed, the φ M-step fits the
exposure variance from whatever ω spread exists; under genuinely sparse gDNA it
could be regularized toward small (strong shrinkage) so "sparse evidence ⇒ trust
ω≈1" holds without relying on ρ₀→0 exactly.

## Conceptual issues to fix (spliced-fragment handling)

1. **"Contained spliced" is a near-empty, wrong category.** A spliced fragment
   jumps exons across an intron, so by definition it spans **multiple regions**
   → it lands in the *boundary* views, never *contained*. Contained spliced can
   only arise from **unannotated** transcripts (novel intronic/intergenic
   splicing absent from the GTF) — a tiny fraction. Yet `update_eps_s` trains
   `eps_s` from the **contained** spliced counts (≈0) while the E-step applies
   `eps_s` to the **boundary** spliced mass (the real junctions). Trained from
   the empty place, applied to the populated place, with a degenerate value.

2. **`eps_s` is essentially untrainable from observation.** A splice artifact is
   a gDNA fragment that *misaligns* as spliced. The reliable detector is the
   `alignable` splice-junction blacklist (index `splice_blacklist.feather`,
   applied at scan by `filter_blacklisted_sjs`) — it removes artifact-junction
   reads upstream so they are never counted as spliced. Without a blacklist,
   there is no robust signal: modeling would need a spliced↔gDNA association per
   region, and at **zero gDNA there is no data**. The honest default is a small
   conservative constant (a "splicing error rate", ~1e-3), not a trained value.

## Strand finding (`fit_strand_balance`)

`p_r1_sense = 0` **exactly** (all 39 spliced reads measured antisense) →
`kappa_rna` clamped to `_BB_FLOOR = 1e-6` → the RNA strand Beta-Binomial is
razor-sharp (a≈1, b≈1e6), so *any* sense read scores as gDNA. Two concerns:
(a) **perfect-stranding degeneracy** — a perfectly stranded library shouldn't
collapse the strand channel to a knife-edge (and `rho_r_bb` also floors at 1e-6,
removing the softening); (b) **boundary double-count** — `fit_strand_balance`
pooled `n_reads=124` from only 39 spliced fragments, because a junction read is
counted in *both* flanking boundary views. (This does not cause the single_exon
over-call — contained reads are still correctly RNA — but it is a real bug.)
Also verify the `p_r1_sense=0` *direction* is the simulator's true R1-antisense
convention and not an orientation flip.
