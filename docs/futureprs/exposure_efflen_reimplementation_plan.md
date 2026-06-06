# Re-implementation plan — IPR exposure-weighted gDNA effective length

**Status:** authoritative spec for re-establishing the **inverse-participation-ratio (IPR)
exposure-weighted gDNA effective-LENGTH** — the capability that contracts the gDNA component's
effective support onto capture-enriched exons. Some of it already survives in `priors.py`
(encouraging). Written to survive a context compaction.

## 0. Scope correction (read first — supersedes earlier framing)
- **The IPR was the real solution** (user is certain), and its skeleton already lives in
  `priors.py` (the `(Σg)²/Σ(g²/geom)` block). We are not starting from zero.
- **Phase 1 = the IPR exposure-weighted effective length WITHOUT shrinkage / regularization.**
  Get the un-shrunk contraction provably correct first.
- **Shrinkage is deferred and its form is undecided.** There were (at least) two shrinkage
  approaches historically — an early **Gamma–Poisson EB posterior** (`exposure.py @ c0af5b3^`) and
  a more-recent **concentration-trust shrink** (the `gdna_locus/(gdna_locus+shrink_count)` blend
  currently in `priors.py`). The Gamma posterior was **not** replaced-then-lost as a regression; it
  was an earlier form. **Both are parked** (see §7). We will decide *whether and what* shrinkage is
  needed empirically, from the toy sweep (§5) — sparse/noisy nodes may over-contract; that is the
  symptom shrinkage would treat.

## 1. The problem this fixes (measured)
Benchmark `scripts/sim/bench_calibration.py` on the priority suite
`/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_20mb_ss99_gdna4x`
(config `scripts/sim/configs/hyb_capture_20mb_ss99_gdna4x.yaml`; 20 Mb, ss 0.99, capture on,
gDNA 4:1):
- **gDNA→RNA leak 36.2%** (434,660 frags, ALL in-locus/on-exon); **pool gDNA 52% vs 80% truth**.
- A/B proof (temp scale on `gdna_eff_len`): contracting eff-len 10× drops the leak 36→13% and
  lifts pool 52→76%. So the dominant lever is the **gDNA component's per-position density =
  gdna_prior_count / gdna_eff_len**; the eff-len is ~the exonic footprint, too diffuse for the
  capture-enriched probed exons — the gDNA component can't win on-exon fragments per-fragment.
- 50%-gDNA / non-capture cases are already good (memory `calibration_benchmark_findings`).

## 2. The model — exposure field + its IPR (un-shrunk)
Per calibration **node** n (region AND each boundary side — substrate 3-view contained/left/right,
kept separate), define the **exposure** (raw local enrichment, NO shrinkage):
```
ω_n = (g_n / L_n) / gdna_density_global        # local gDNA density ÷ global gDNA density
```
`g_n` = deconvolved gDNA mass on the node (post boundary-flux transport); `L_n` = node length;
`gdna_density_global = Σg / ΣL`. `derive.py` already computes exactly this as the (currently
QC-only) `region_exposure` / `left_exposure` / `right_exposure` — un-shrunk.

**Per-locus gDNA effective support length = IPR of the exposure field:**
```
gdna_eff_len = (Σ_n∈L  ω_n · L_n)²  /  Σ_n∈L (ω_n² · L_n)
```
This is the "effective number of participating bases" of the exposure-weighted field. Limits
(all verified by hand):
- **Uniform exposure** (ω_n ≡ c): `= (cΣL)²/(c²ΣL) = ΣL` → full gDNA support (no contraction).
- **Perfect capture** (all mass on one probed node of length L\*, others ω≈0): `= (ωL*)²/(ω²L*) =
  L*` → contracts to the probed footprint.
- Monotone in between: more enrichment → shorter eff-len → higher gDNA per-position density →
  gDNA competes for on-exon fragments at its true local density.

**Algebraic identity (why the skeleton is already in `priors.py`):** substitute ω_n=(g_n/L_n)/ρ₀:
```
Σ ω L  = Σ g/ρ₀ ;   Σ ω²L = Σ g²/(ρ₀²L)  ⟹  (Σ ω L)²/Σ(ω²L) = (Σ g)² / Σ(g²/L)
```
ρ₀ cancels. So the exposure-field IPR **equals the mass-IPR `(Σg)²/Σ(g²/L)` on bare node lengths
L**. `priors.py` already computes `(Σg)²/Σ(g²/geom)` — the SAME form, but (a) it divides by
`geom = gdna_geom_len` (FL-aware support ≈ R+L̄) instead of bare `L`, and (b) it then shrinks toward
`span` via `concentration_trust`. **Phase 1 = drop (b), and decide (a) empirically.**

## 3. The fix (surgical changes to `priors.assemble_priors`)
1. **Remove the concentration-trust shrink** (defer to §7): use `gdna_support_len` directly as
   `gdna_eff_len` (still floored by `_GDNA_EFF_LEN_FLOOR=1.0`, a numerical guard, not a shrink).
2. **Length basis — `geom` vs bare `L`:** the theory/IPR uses bare physical node length `L`; the
   current code uses FL-aware `geom`. **Resolve on the toy (§5):** bare L gives sharper contraction
   for short exons; geom is FL-consistent with the EM's per-position rate. Pick the one that makes
   the toy's perfect-capture eff-len land on the probed footprint. (Decision recorded in §6.)
3. **Mass is NOT scaled by exposure** — `gdna_prior_count` stays the raw deconvolved gDNA count.
   Exposure modulates LENGTH only.
4. **Keep** the boundary-flux transport upstream (`_transport_boundary_flux`, present) and the
   overlap-share projection (`_project_regions_to_loci`).
5. **Anti-collapse:** the numerical floor only. Over-contraction on *sparse* nodes is the deferred
   shrinkage's job — first measure whether it actually happens on the toy.

No new constants in Phase 1 (honors `calibration_no_magic_numbers`). `exposure_dispersion` (φ) is
NOT restored — it belongs to the deferred Gamma-posterior shrinkage. `gdna_eff_len_shrink` becomes
unused in Phase 1 (set to 0 / bypass the blend; do not delete the config field yet — discuss).

## 4. What survives / what's parked (git forensic, corrected)
- **Survives (Phase-1 reuse):** `priors._transport_boundary_flux` (boundary transport iterations,
  present); `_project_regions_to_loci`; the `(Σg)²/Σ(g²/geom)` IPR skeleton; `derive.py`'s un-shrunk
  `region/left/right_exposure`. The rename refactor was pure.
- **Parked, kept handy (NOT Phase 1):** the Gamma–Poisson EB posterior at `c0af5b3^`
  (`src/rigel/calibration/exposure.py`, `tests/calibration/test_exposure.py`) and
  `CalibrationConfig.exposure_dispersion`; the current `concentration_trust` shrink. See §7.

## 5. Toy scenario — derive-from-theory + empirical + contractual tests (the instrument)
Build a tiny, transparent scenario where the correct eff-len is known by construction, then make
the behavior contractual.

**Geometry:** 1 Mb genome, **10 genes / transcripts**. Each transcript = **two 1 kb exons** with a
long intron between them (~20 kb → genomic span ≈ 22 kb; "spanning 2 kb" = 2 kb *exonic*). Genes
spaced across 1 Mb with intergenic gaps (~22 kb span × 10 ≈ 220 kb used; ample spacing). Single
isoform per gene. (Geometry is a fixture detail — confirm/trivially adjust at build; it does not
change the method. The point is a stark exon/intron/intergenic contrast so contraction is dramatic
and assertable.)

**Capture-energy sweep** (one suite, `capture.configs` list — knobs confirmed in
`src/rigel/sim/capture.py`: mass `= off_target_weight·eff_len + binding_per_base·probe_overlap`):
- `off`     — `enabled: false` (uniform; no capture).
- `weak`    — `off_target_weight: 1, binding_per_base: 1`.
- `strong`  — `off_target_weight: 1, binding_per_base: 10`.
- `extreme` — `off_target_weight: 0, binding_per_base: 100` → **zero off-target baseline ⇒ ALL gDNA
  on probed exons, zero intron/intergenic reads** = the user's "perfect capture."
`probe_density` sets the probe-covered fraction of each exon (the perfect-capture footprint).

**Contractual expectations (the unit tests):**
- `off`: `gdna_eff_len ≈ full gDNA support` (≈ per-locus span; uniform exposure → IPR = ΣL).
- `weak → strong → extreme`: `gdna_eff_len` **decreases monotonically**.
- `extreme`: `gdna_eff_len ≈ probed-exon footprint` (≈ probe-covered bases of the 2 kb exonic),
  **NOT the ~22 kb span** — i.e. ~10× contraction. gDNA→RNA leak ≈ 0 at extreme capture.
- Off-target (intron/intergenic) gDNA is always classified gDNA (the easy direction).

**Workflow:** write a sim config (mirror `hyb_capture_20mb_ss99_gdna4x.yaml`, shrunk to the toy +
the 4-point capture sweep) → `simulate_suite.py` → `bench_calibration.py --run` → inspect
`gdna_eff_len` per locus across the sweep (add a debug dump if needed). Iterate the §3 length-basis
choice until the contraction matches. The toy should show **excellent** performance.

## 6. Implementation steps
1. **Exposure per node (un-shrunk):** ensure `ω_n` (region + left + right) is available to the prior
   — either consume `derive.py`'s existing `*_exposure` or compute the mass-IPR form directly in
   `priors.py` (identical; simpler). No Gamma posterior.
2. **`priors.assemble_priors`:** `gdna_eff_len = (Σ ω L)²/Σ(ω²L)` per locus (= `(Σg)²/Σ(g²/L)`);
   drop the `concentration_trust` blend; resolve bare-L vs `geom` on the toy (record the decision in
   the docstring). Keep transport + projection + numerical floor.
3. **Config:** bypass `gdna_eff_len_shrink` (no shrink in Phase 1); do NOT restore
   `exposure_dispersion`. Discuss before adding/removing any constant.
4. **Mass unchanged.**

## 7. Deferred — shrinkage candidates (documented, kept handy; NOT Phase 1)
Decide IF/WHAT after the toy reveals whether sparse nodes over-contract. Two documented options:
- **Gamma–Poisson EB posterior** (early form): `M_g ~ Poisson(ρ₀·ω·L)`, `ω ~ Gamma(s,s)`, `s=1/φ` ⇒
  `ω̂ = (s+M_g)/(s+ρ₀L)` (convex pull of raw enrichment toward 1; sparse→1, dense→m),
  `Var(log ω̂)=1/(s+M_g)`. Code: `git show c0af5b3^:src/rigel/calibration/exposure.py` and
  `tests/calibration/test_exposure.py`; theory: `docs/futureprs/exposure_prior_robustness_theory.md`,
  `exposure_prior_noise_floor_theory.md`. Would feed shrunk `ω̂` into the §2 IPR in place of raw ω.
- **Concentration-trust shrink** (more-recent form, currently in `priors.py`): blend
  `gdna_eff_len = (1-trust)·span + trust·gdna_support_len`, `trust = gdna_locus/(gdna_locus +
  shrink_count)`. A per-locus count-based shrink toward the uniform span.

## 8. Tests + validation
1. **Toy-sweep contractual unit tests** (§5) — the primary new tests; make the contraction +
   monotonicity + perfect-capture footprint contractual.
2. **Capture-leak regression** on `hyb_capture_20mb_ss99_gdna4x`: assert leak drops well below 36%
   (target ≤10%, iterate) and pool gDNA approaches 80%.
3. **No-regression:** 500 kb suite (pool ±2%, FP ~0) + unit + golden suites stay green; regen
   goldens only after the calibration is finalized.
4. **Validate** via `scripts/sim/bench_calibration.py --sim-base <suite> --run --force` on BOTH the
   toy and the 20 Mb 4:1 capture suite (+ 500 kb), reporting pool + per-fragment confusion.

## 9. Git artifacts to re-read in the fresh research phase
- current: `src/rigel/calibration/priors.py` (the IPR skeleton + transport + projection),
  `derive.py` (un-shrunk `*_exposure`), `config.py`; `src/rigel/sim/capture.py` (sweep knobs);
  `scripts/sim/bench_calibration.py`, `simulate_suite.py`, `scripts/sim/configs/hyb_capture_20mb_ss99_gdna4x.yaml`.
- parked shrinkage (only if §7 is revisited): `git show c0af5b3^:src/rigel/calibration/exposure.py`,
  `git show c0af5b3^:tests/calibration/test_exposure.py`, `git show 266ba70 -- src/rigel/calibration/`.
- theory: `docs/futureprs/phase6_exposure_efflen_theory.md`, `exposure_prior_robustness_theory.md`,
  `exposure_prior_noise_floor_theory.md`, `exposure_prior_findings_report.md`.
