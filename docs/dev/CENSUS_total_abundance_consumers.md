# CENSUS — every place an abundance / density / level is formed, and what a measured-total swap does to it

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** Rulings live in `DESIGN.md` §3.1a-ii,
    derivations in `EQUATIONS.md` §2.3b, lessons in `TRAPS.md`, current numbers in `ROADMAP.md` §1
    rank 3. This file exists because the owner asked for the ENUMERATION (2026-08-21) and an
    enumeration is a working artefact: it goes stale the moment a site moves. **Re-derive it, never
    carry it.**

## 0. Why this exists, and the error it corrects

A prior session declared Part C rung ② ("replace existing abundance readings per consumer") BLOCKED,
on the grounds that `ROADMAP.md` §4.3 had refused the swap. ⛔ **That was wrong twice over**, and the
owner caught it:

1. **§4.3 is scoped to ONE consumer.** It refused `region_start_count / ℓ` substituted for
   `RegionGeometry.inv_abundance` *inside `messages/currency.py`'s `enrichment_ratio`* — unmasked, with
   no side selection and no model-free flag — judged on that policy's own deliverable. It is a verdict
   about a policy, not about the quantity.
2. **"Readers of `inv_abundance`" is not "places total abundance is used."** The first set has one live
   member; the second has **240**.

## 1. The census — 240 sites, and how to regenerate it

Eight parallel readers over module clusters plus three adversarial completeness critics; every site
classified by axis / form / component / population / consumer, with a swap verdict. Machine-readable
dump: rebuild it rather than trusting a copy.

| | count |
|---|---|
| sites total | **240** (194 in `src/`, the rest in `scripts/design/`) |
| `NO_SWAP` | 166 in `src/` |
| `NEEDS_MEASUREMENT` | 20 in `src/` |
| `SWAP` | 8 in `src/` |

**The four FORMS, and only two of them are defects:**

| form | expectation | verdict |
|---|---|---|
| `contained_inv_bank` — `Σ 1/(ℓ−w+1)` | `ρ·P(w ≤ ℓ)` — **TRUNCATED** | model-free but biased; the defect `EQUATIONS.md` §2 names |
| `count_over_eff` — `count / E_contained` | `ρ` — **unbiased** | correct, but the fl pmf enters the DIVISOR |
| `boundary_inv_bank` / `sj_inv_bank` — `Σ 1/(w−1)` | `ρ` exactly | ⭐ already model-free AND unbiased; leave alone |
| `mass_over_eff` where numerator and divisor are COMPONENT-MATCHED | `ρ_c` | ⛔ **NOT a total-abundance site.** A deconvolved gDNA mass over the gDNA opportunity is legitimate; swapping a TOTAL in would add RNA |

⭐⭐⭐ **WHERE TOTAL ABUNDANCE IS ACTUALLY NEEDED — the owner's list, 2026-08-21, and it is short:**
① the DENSITY MODEL (the gDNA background and the intron deconvolution); ② the TOTAL-ABUNDANCE LANDSCAPE
(planned, rung ③); ③ ENRICHMENT RATIOS, wherever they are formed (a ratio of two totals). ⛔ Everything
POST-calibration is out of scope by construction: once beliefs exist, each component has its own fl
distribution and its own opportunity, so a per-component estimate is the right instrument and a total
is the wrong one.

⭐⭐ **The distinction that governs every verdict**: a site is a swap candidate only where the quantity
wanted is a TOTAL (or a gDNA level at a structurally pure-gDNA object, where total ≡ gDNA). 61 sites
are `gDNA_only_structural_selection` — those qualify; 50 are `gDNA_only_deconvolved` and 46
`component_split_by_belief` — those do not.

## 2. What was SWAPPED, and what it measured

Landed behind `CalibrationConfig.background_abundance` (`"contained"` default ⇒ bit-identical;
`"measured_total"` ⇒ the START/END pair). Both consumers keep their own estimator — only the
`(counts, exposure)` pair changes, so the Gamma conjugacy, the pool predicate and the MAD trim are
untouched, and a double-walled region excludes itself through the existing zero-support predicate.

* `density_deconv.fit_intron_background` → the intron λ-factor → **ψ**
* `background_reference.measure_background` → the refit floor → the npmle floor

⭐⭐ **Priced at the ESTIMATOR level on all 32 conditions** (`total_abundance_audit.py` arm ⓕ), against
the `gdna` origin partition's own start rate — a truth neither estimator's machinery targets:

| pool | capture-OFF shipped / new | capture-ON shipped / new |
|---|---|---|
| intergenic (reaches ψ) | 1.0000–1.0001 / 1.0000 | **0.5627–0.5644 / 1.0818–1.0830** |
| +introns | 1.1821–1.1823 / 1.2055–1.2058 | **0.2344–0.2348 / 1.0649–1.0654** |

**The verdict: a TIE off capture and a 1.8–4.3× repair under it.** The mechanism is that the gDNA fl
pmf is itself capture-distorted, so the shipped divisor is mis-estimated; a pmf-free exposure is
immune. ⛔ The ONE cost is visible on the +introns rows and it is not the estimator's: those pools
carry nascent RNA, both forms over-read, and the START form takes ~2 pp more of it because a fragment
starting in an intron and reaching into an exon books a START there. Where the pool is clean the term
is absent.

⚠ **The exon population is where the truncation actually lives**: the shipped reciprocal bank reads
**0.3197–0.3439** of the true total on every stratum. Exons are in no pure-gDNA pool, so no consumer
above touches them — the one consumer that does is the currency channel, which is §4.3's territory.

## 3. ⛔⛔⛔ RETRACTED — `capture_eff_length.py:194` / `priors.py:353` ARE NOT TOTAL-ABUNDANCE SITES

This section claimed they were "the highest-value target" for the measured total. **The owner rejected
that the same day and was right; the section is kept as a retraction rather than deleted, because the
mistake is the instructive part.**

The readers flagged those two lines because the NUMBER is dramatic (the EM ruler contracts by a mean
factor **0.345**, and **0.0951** at `g00`, where the correct factor is exactly **1.000**). But:

* `transcript_capture_eff_lengths` takes a **solved** `CalibrationResult`. It runs POST-calibration.
* So `_global_reference_density(contained_m, contained_S)` is `mass_gdna_region / gdna_region_eff_len`
  — a DECONVOLVED gDNA mass over the gDNA opportunity. **Component-matched.**
* ⭐ **Once beliefs exist there is no total-abundance question left**: each component has its own fl
  distribution and its own opportunity, and using them is correct — which is exactly what those two
  files do.

That is §1's own form table (row 4: a component-matched `mass_over_eff` is NOT a total-abundance site)
applied to the site the readers were loudest about. ⭐ **The 0.0951 is a COMPOSITION defect**, and its
repair is the one `ROADMAP.md` §1 rank 3 already records — substitute the composition arrays and the
SHIPPED shrinkage returns 1.000 — so it belongs to whatever fixes `f_g`, not to this channel.

⚠ **The lesson for the next census**: a subagent can rank a site by the MAGNITUDE of its error, and
magnitude does not tell you which repair applies. The form rule does, and it has to be applied to the
loudest finding rather than around it.

## 4. The remaining `NEEDS_MEASUREMENT` set in `src/`, and what each needs

Each is a site where the total-abundance estimator is admissible on a SUBSET of the population and
forbidden on the rest, so the swap needs a mask rather than a substitution:

* `landscape.py:205/349` + `calibrate.py:265` (`_fit_gdna_hyperprior`) — the substrate is a MIXTURE:
  single-strand slots are component-matched (leave), `~free_pos & ~free_neg` slots are structurally
  pure gDNA (admissible). The Poisson kernel's shape survives the swap unchanged —
  `P(count | ρ·E)` takes `(region_start_count, ℓ)` as readily as `(f_g·M, eff_gdna)`.
* `region_init.py:367` (`rho_g`) — deconvolved for `f_g < 1`, but `f_g ≡ 1` on the pure-gDNA subset.
* `messages/relay.py:312` (`region_total_density`, the reframe frame pair) — flagged as *"the
  highest-value site §4.3 does NOT cover"*: it is a TOTAL formed as `mass_over_eff`, so its divisor is
  a function of the composition being solved for. ⚠ Policy work is frozen pending this thread.
* `calibrate.py:194/196` (`_build_intron_prior`) — both sides live in CONTAINED-count units, so
  substituting a START count on one side alone would compare two different populations. Needs the
  pair on BOTH sides or neither.
* `calibrate.py:603` (the npmle enrichment fit) — `region_gdna_geometry` hands it total mass over the
  gDNA opportunity, which IS the forbidden form; its consumers are the QC report and the
  toy-injection substrate, so the blast radius is a verdict rather than the tool. ⚠ Toy injection
  means a wrong landscape there can mislead a toy experiment.

## 4b. ⭐⭐⭐ WHAT THE ADVERSARIAL PASS FOUND THAT THE READERS MISSED — three critics, and they earned it

The eight readers covered `src/rigel/calibration/` well and missed whole TREES. The critics' finds, in
value order:

1. ⛔⛔⛔ **`measure_background`'s output is DEAD for the solve** — assigned in `calibrate`, then read
   only by the `_debug["calibration_priors"]` bundle (`calibrate.py:740`). `background=` occurs exactly
   once in all of `src/rigel/calibration/` outside its own module; the sole `DensityNPMLE.fit` call site
   omits it; `fit_landscape` has no such parameter. Its "refit floor" role was retired when
   `landscape.py` replaced `npmle`'s hyperprior role. **So one of the two landed swaps has no live
   consumer, and the ~1 % end-to-end move came entirely from `fit_intron_background`.** An owner call:
   wire it, or converge-and-delete.
2. **`region_geometry.py:470/479` `region_total_density` — verdict CORRECTED to SWAP.** Its own
   docstring calls it "the LAZY, composition-aware total density": a TOTAL formed as `mass_over_eff`,
   so the divisor is a function of the composition being solved for. ⭐ And `_flux` at :479 has its
   exact model-free replacement already on the same dataclass at the same shape —
   `geometry.inv_sj_lo` / `inv_sj_hi`.
3. **`relay.py:669` (scalar) and `relay.py:812` (vector) — the reframe RATIO, missed entirely.** The
   readers found where `rho_lo`/`rho_hi` are BUILT and not where the ratio is FORMED, and the ratio is
   the quantity. ⛔ The two are gated on bit-for-bit equality, so a one-sided edit is the failure mode.
   ⚠ Policy work is frozen pending this thread.
4. **Whole trees absent from the census**: `src/rigel/estimator.py` (the transcript AND locus axes — a
   module literally named for abundance), `src/rigel/report/`, `src/rigel/cli.py`, and every C++ site
   (`em_solver.cpp`). All classified NO_SWAP on inspection — component-matched by construction — which
   is the reassuring outcome, but the census could not have said so before they were read.
5. Two incidental defects worth their own tickets: `estimator.py:976`'s `gdna_eff_len_per_bp` is DEAD
   (`pipeline.py` writes only `gdna_eff_len_em`, so the golden `_loci_df.tsv` carries a column of
   zeros), and `report/model.py:450` divides a deconvolved gDNA mass by RAW bp where
   `calibration/track.py:37` divides the SAME mass by an effective length — two densities under one
   name.
6. ⚠ A doc self-contradiction: `EQUATIONS.md` §1.5 gives the crossing opportunity as
   `E_f[min(w−1, R_lo, R_hi, R_lo+R_hi−w+1)₊]` and warns of a 4× error at a first exon, while §2.1
   substitutes `(w−1)` unconditionally. Both are in scope of the same claim and they disagree.

## 5. `scripts/design/` — 46 sites, and a broken form here corrupts a VERDICT

`object_composition.py` (6 sites), `anchor_opportunity_census.py`, `worst_objects.py`. These are
measurement code: they do not ship, but they are what the project's rulings are read off. None was
found to be using a form its own question forbids; the census records them so a future reader does not
have to re-derive that.
