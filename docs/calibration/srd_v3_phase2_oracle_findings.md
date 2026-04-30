# SRD v3 — Oracle benchmark findings & next-phase plan

**Status:** post-categorization-fix; pre Phase-2 (joint FL × strand intronic mixture)
**Inputs:** VCaP simulated oracle BAMs (40M PE fragments per condition, R1-antisense, SS=0.99)
**Index:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/rigel_index` (rebuilt — no synthetic-nRNA mega-exons)

## 1. Summary of fixes landed in this iteration

| # | Bug | Fix |
|---|---|---|
| 1 | Stale index treated synthetic nRNAs as `EXON`-type intervals → `INTRONIC` always 0 | Index rebuilt at new path |
| 2 | R1-antisense library produced inverted SENSE/ANTISENSE labels (99:1 reversed) | Threaded `read1_sense` from `StrandModel` through `_compute_fragment_strand` |
| 3 | Per-block aggregate `exon_bp_pos / _neg` summed across cgranges hits and across paired-end blocks; used vs `genomic_footprint` (incl. mate gap) → 1.85M false `EXON_INCOMPATIBLE` on perfect oracle data | New per-T `max_T(exon_bp[T])` check via `np.maximum.at` scatter; default `exon_fit_tolerance_bp = 0` |
| 4 | `gdna_fraction` reporting halved by intergenic split | Redefined as `(gdna_em + intergenic) / total`; added `gdna_em_fraction` and `intergenic_fraction` |

## 2. Validated outcomes (oracle benchmarks)

### Pristine — `gdna_none_ss_0.99_nrna_none` (truth: 100% mRNA)

| Metric | Pre-fix | Post-fix | Change |
|---|---:|---:|---|
| `EXON_INCOMPAT` | 1,848,420 | **0** | ✅ −100% |
| Calibration `n_pool` | 1,953,964 | **0** | ✅ correctly empty |
| Quality | "good" (false signal) | **"fallback"** | ✅ |
| `mrna_fraction` | 0.99763 | 0.99734 | ✓ |
| Strand on `EXON_CONTAINED` | — | 99:1 SENSE:ANTI | ✓ matches SS=0.99 |

### Dna20m — `gdna_dna20m_ss_0.99_nrna_none` (truth: 50% mRNA, 50% gDNA)

| Metric | Pre-fix | Post-fix | Truth |
|---|---:|---:|---:|
| `EXON_INCOMPAT` | 2,929,930 | **103,830** | ~0 (residue is real exon-edge gDNA) |
| Mixture-EM `pi_pool` | 0.395 | 0.4548 | (gDNA share of pool) |
| `mrna_fraction` | 0.4987 | **0.4986** | 0.50 ✓ |
| `gdna_fraction` (incl. intergenic) | 0.4222 | 0.4182 | 0.50 (16% under) |
| `nrna_fraction` (siphon!) | 0.0791 | **0.0832** | 0.0 |

The intronic pool in dna20m is **exactly 50:50 strand-symmetric** (4,686,883 SENSE : 4,687,486 ANTISENSE) — the unmistakable signature of pure-gDNA origin — yet ~3.3M reads still leak to nRNA at the locus EM stage.

## 3. Open issues & recommendations

### P1 — nRNA siphon (highest-impact remaining bug, ~8% mass-loss)

The locus-level EM is misclassifying gDNA as nascent RNA even when calibration knows the truth. Suspected contributing factors:

- **(a) Strand symmetry not enforced at the locus level.** Calibration sees the 50:50 intronic strand split as "gDNA-like" but per-locus EM has no penalty for breaking that symmetry. With stranded data, EM can siphon SENSE fragments to nRNA, leaving an antisense-skewed gDNA residue that the model accepts because no symmetry constraint is in place.
- **(b) gDNA effective-length is wrong.** Current code uses a **harmonic mean of locus transcript lengths** as a proxy for the gDNA effective length. This is fundamentally incorrect: gDNA effective length is determined by the genomic span (and the gDNA fragment-length model), not by transcript length. Wrong effective lengths → biased likelihoods → systematic misallocation. For *unstranded* data this is the only way nRNA siphon can occur.
- **(c) nRNA shadow priors per locus are independent of strand-symmetry evidence.** `compute_locus_priors_from_partitions` derives per-locus Dirichlet priors from per-fragment SRD posteriors. If those posteriors don't encode strand symmetry of the intronic sub-pool, the locus prior is uninformative and EM can drift.

### P2 — Categorization tolerance for real (STAR/minimap-aligned) data

Default is now `exon_fit_tolerance_bp = 0`. Correct for oracle, possibly too strict for soft-clipped real data. Re-test on STAR-aligned VCaP runs and expose as user-tunable.

### P3 — Pristine fallback path

`quality="fallback"`, `n_pool=0`, `pi_pool=0`, `mixture_iterations=0` is the correct response. Add a unit test asserting downstream priors don't NaN in this regime.

### P4 — Cosmetic
- `strand_model.label` missing from summary JSON.
- Promote `gdna_fraction / gdna_em_fraction / intergenic_fraction` triple to top of summary.

## 4. Phase 2 plan — defeat the nRNA siphon

The investigation will proceed bottom-up from a single-transcript mini-genome to confirm each contributing factor in isolation, then up to multi-locus. Detailed plan and progress to be tracked in `srd_v3_phase2_nrna_siphon.md` (next document).
