# Post-PR07 calibration deep-dive: the 6 remaining scenario failures

After PR07 (eps_s removal, spliced→RNA) the scenario suite dropped from 12 → 6
failures but PR07 *unmasked* deeper issues. This is the investigation record; the
fixes are split across PR8 (AMBIG ρ₀ runaway), PR9 (strand knife-edge), PR10
(spliced double-count).

## The 6 failures and their two mechanisms

Consolidating run (full pipeline; ρ₀ and `conv` separate the two classes):

| case | ρ₀ | disp | conv | κ_rna | ρ_r | RNAobs/exp | gDNA | nRNA |
|---|---|---|---|---|---|---|---|---|
| anti_intron nrna30 ss100 | **2.0e4** | 100 | ✗ | 1e-6 | 1e-6 | 384/2000 | 1616 | 0 |
| anti_intron nrna70 ss100 | **2796** | 100 | ✗→✓ | 1e-6 | 1e-6 | 209/2000 | 1791 | 0 |
| overlap_anti fc4 ss100 | **6.2** | 27.6 | ✗ | 1e-6 | 1e-6 | 602/1000 | 398 | 0 |
| nrna_dc g20 n70 ss65 | 0.039 | 5.6 | ✓ | 0.36 | 0.014 | 283/1404 | 1717 | 0 |
| nrna_dc g20 n0 ss90 | 0.041 | 0.27 | ✓ | 0.13 | 1e-6 | 1110/983 | 890 | 118 |
| nrna_dc g20 n70 ss90 | 0.036 | 2.1 | ✗ | 0.11 | 1e-6 | 1386/1404 | 614 | 771 |

**Mechanism 1 — ρ₀ runaway (dominant; PR8).** The zero-real-gDNA scenarios with an
AMBIG region (top 3) show ρ₀ exploding (true ρ₀ ≈ 0), dispersion at ceiling,
non-convergence, nRNA=0. See the full audit below.

**Mechanism 2 — weak strand separation at low ss (PR9).** `nrna_dc g20_n70_s65` has
a healthy converged ρ₀ but κ_rna=0.36 sits close to gDNA's 0.5 → the strand channel
barely separates and RNA leaks to gDNA. (`g20_n70_s90` is essentially a pass failing
only a borderline neg-control check, 28 vs limit 25.)

## Audit: how AMBIG regions seed the ρ₀ runaway (PR8)

Per-iteration ρ₀ on `anti_intron_nrna_30`: `0.20, 0.31, 0.48, 0.74, 1.14, 1.77,
2.75, … , 502 (disp hits ceiling 100), … , 20390` — geometric, never converges
(delta stuck ≈ 0.019).

- An AMBIG region (genes on both strands, e.g. g2's exon inside g1's intron) is
  **doubly undeconvolvable**: strand channel correctly masked (estep `_llr_strand`
  excludes `TS_AMBIG`); no boundary-crossing fragments (contained).
- Iteration 1: no evidence → `π_g = pi_g_prior = 0.5` (init), count channel silent
  (`m_d_prev=0`). **Half its reads become false gDNA.** In a zero-gDNA library this
  is the only "gDNA" present.
- `update_rho_0 = Σ M_g/Σ(ω·L)` sums over **all** regions → the false gDNA drives ρ₀.
- `pi_g_prior = ω·ρ₀·L/n_u` feeds ρ₀ back into the AMBIG region's prior → more gDNA.
- Decodable regions correctly resolve to RNA → their ω → 0 → denominator collapses.
- Numerator pinned by the AMBIG region, denominator collapsing ⇒ ρ₀ underdetermined
  ⇒ geometric runaway. Circular: **ρ₀ → prior → AMBIG gDNA → ρ₀.**
- The sweep compounds it: `alpha_swept` includes the region's own `a_reg`
  (self-referential imputation).

The converged-state per-region table confirms every decodable g1 region is RNA
(`gfrac=0`); only the AMBIG region carries gDNA (`gfrac=0.58`). PR07 unmasked this —
the old `eps_s=1` phantom exon gDNA had incidentally anchored ρ₀.

Fix → PR8: estimate ρ₀ from strand-decodable regions only (AMBIG regions are
*imputed*, not voters); impute AMBIG exposure from neighbours (drop self `a_reg`).

## Reframe: the strand "inversion" is NOT a bug

The old PR8 design called `p_r1_sense = 0` (κ_rna=1e-6) "the clear bug." The audit
shows it is a **self-consistent convention**. The R2-flip (bam_scanner.cpp:647)
puts `align_strand` in the R1-orientation frame; both κ_rna (spliced) and the
unspliced `k_sense` (confirmed `n+=0, n−=n_u` for POS regions) live in that *same*
frame, so the BB channel still separates RNA (align-match ≈ κ_rna, an extreme) from
gDNA (= 0.5), and single_exon resolves correctly. The actionable item is only the
misleading "both report the correct transcript strand" docstring in strand_model.py
(fold into PR9).

The real strand bug is the **knife-edge** (PR9): at perfect stranding κ_rna and
ρ_r_bb both floor at 1e-6 → razor-sharp BB → any off-strand read scores gDNA. Plus
the **double-count** (PR10): junction reads pooled on both boundary sides (124 reads
from 39 fragments).

## Diagnostic scripts (scripts/debug/, untracked)

`summarize_failing_scenarios.py` (the table above), `trace_antisense_nrna.py`
(per-region channels), `trace_zero_gdna.py`, `trace_single_exon.py`.
