# Re-implementation plan — EB-shrunk exposure as a gDNA effective-LENGTH factor

**Status:** authoritative spec for re-implementing a capability that was implemented + tested in
the cyclic calibration, then **lost in the cyclic→acyclic teardown** and never restored. Written
to survive a context compaction. Recovered entirely from git + the theory docs (refs below).

## 0. The problem this fixes (measured)
Benchmark `scripts/sim/bench_calibration.py` on the priority suite
`/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_20mb_ss99_gdna4x`
(config `scripts/sim/configs/hyb_capture_20mb_ss99_gdna4x.yaml`; 20 Mb, ss 0.99, capture on,
gDNA 4:1):
- **gDNA→RNA leak 36.2%** (434,660 frags, ALL in-locus/on-exon); **pool gDNA 52% vs 80% truth**.
- A/B proof (temp scale on `gdna_eff_len`): contracting eff-len 10× drops the leak 36→13% and
  lifts pool 52→76%. So the dominant lever is the **gDNA component's per-position density =
  abundance / gdna_eff_len**; the eff-len is ~the exonic footprint, far too diffuse for the
  capture-enriched probed exons.
- 50%-gDNA / non-capture cases are good (see memory `calibration_benchmark_findings`).

## 1. What was lost (git forensic)
- `src/rigel/calibration/exposure.py` @ **`c0af5b3^`** — the EB-shrunk exposure (Gamma posterior).
  Deleted in `c0af5b3` (Phase-6 teardown).
- `tests/calibration/test_exposure.py` @ **`c0af5b3^`** — the EB-shrinkage unit tests.
- `CalibrationConfig.exposure_dispersion` (φ) — purged in the acyclic rebuild.
- Exposure-as-LENGTH was then turned off by **"Option A"** (dangling WIP `266ba70`
  "option-A-geometric-efflen", 2026-06-04): `gdna_exposure_len = Σ ω·L` → `gdna_geom_len = Σ L`
  (geometric, ω-free; "exposure lives only in the mass, never the length").
- The acyclic `derive.py` kept a per-node exposure ω but as an **un-shrunk ratio**
  (`region_exposure = (gdna_mass/eff_len)/gdna_density_global`), **consumed by nothing**
  (QC-only). The current eff-len (`priors.py`) is the **mass**-IPR on **geometric** lengths.
- NOT lost: boundary-flux transport iterations (`priors._transport_boundary_flux`, present);
  the rename refactor (A–I) was pure (eff-len logic identical pre/post). No un-merged calibration
  work on any branch; dangling commits are only WIP stashes (all recoverable).

## 2. The recovered model (Gamma–Poisson EB exposure)
From `exposure.py @ c0af5b3^` + `docs/futureprs/exposure_prior_robustness_theory.md`:
```
M_g  ~  Poisson(ρ₀ · ω · L)          # per-node gDNA mass given local exposure ω
ω    ~  Gamma(s, s),  s = 1/φ         # prior mean E[ω]=1 (uniform), φ = exposure_dispersion
ω̂    = E[ω|data] = (s + M_g) / (s + ρ₀·L)        # EB-shrunk posterior mean
Var(log ω̂) = 1 / (s + M_g)
```
`ω̂` is a convex pull of the raw enrichment `m = M_g/(ρ₀·L)` toward 1:
`ω̂ = [s·1 + (ρ₀L)·m] / [s + ρ₀L]`. **Sparse support (small ρ₀L) → ω̂→1 (uniform); dense → ω̂→m.**
(The "1/2, 10/20, 100/200 same ratio but shrink the sparse ones" behaviour.) Computed **per node
— region AND each boundary side, kept separate** (substrate 3-view: contained / left / right).

## 3. The fix (exposure modulates LENGTH, not mass)
**Mass is NOT multiplied by exposure.** `gdna_prior_count` stays the raw deconvolved gDNA count.
Exposure modulates the **gDNA component effective length**. The gDNA per-position density the EM
sees is `gdna_prior_count / gdna_eff_len`, so a shorter eff-len at enriched loci → higher density →
competitive on the captured exons.

**Per-locus gDNA effective length = IPR of the EB-shrunk exposure field** (user + 
`phase6_exposure_efflen_theory.md`):
```
gdna_eff_len[L] = (Σ_nodes∈L  ω̂·L_node)²  /  Σ_nodes∈L (ω̂² · L_node)
```
- Contracts toward the enriched (high-ω̂) nodes (ω² weighting); uniform ω̂=1 → ΣL (full support).
- Robust because ω̂ is EB-shrunk: a sparse/noisy node's enrichment shrinks to 1 and does NOT
  spuriously contract (or inflate) the support.
- Keep an **anti-collapse guard** (Option A's valid concern): floor the eff-len (and/or ω̂ within
  `[ε, 1/ε]`) so a blind/zero-gDNA node can't drive eff-len→0 or ∞.

**Do NOT** use the old additive `Σ ω·L` (pre-Option-A): with mass∝ω·L it cancels to global density
ρ₀ (no enrichment benefit). Use the IPR form above.

Contrast w/ current `priors.py`: it computes `(Σg)²/Σ(g²/geom)` from the **un-shrunk mass** `g` on
**geometric** lengths `geom=gdna_geom_len` (R+L̄). Algebraically that's the exposure-field IPR with
`geom` instead of bare `L` AND no EB shrinkage. The fix swaps in EB-shrunk `ω̂` and bare `L`.

## 4. Old → new variable map (lost in the rename)
- `ω` / `omega` (EB-shrunk) ← was `exposure.Exposure.omega`; now `derive.region/left/right_exposure`
  (un-shrunk, QC-only) → make these the EB-shrunk `ω̂`.
- `exposure_dispersion` (φ) ← purged `CalibrationConfig` field → restore.
- `gdna_exposure_len` (Σω·L, pre-Option-A) → `gdna_geom_len` (ΣL, Option A, current) → **target:
  exposure-corrected IPR** above.
- `gdna_support_len` (current `priors.py` local) = the mass-IPR = the value passed to the EM as
  `LocusPriors.gdna_eff_len`. This is the line to replace.

## 5. Implementation steps
1. **EB-shrunk exposure per node** (`derive.py`, replacing the un-shrunk ratio): `ω̂ = (s+M_g)/(s+
   ρ₀·L)`, `s=1/φ`, per region + left + right, using `regions/left/right.gdna_mass`,
   `region_eff_len`/`boundary_eff_len`, `gdna_density_global`. Carry `ω̂` (and optionally
   `Var(log ω̂)`) on `CalibrationResult` (replace the QC-only `exposure_*`).
2. **Exposure-corrected eff-len** (`priors.assemble_priors`): per locus, project the nodes'
   `(ω̂·L)` and `(ω̂²·L)` to loci (overlap-share, like the current `_project_regions_to_loci`),
   then `gdna_eff_len = (Σ ω̂L)² / Σ(ω̂²L)`, floored. Replace the current mass-IPR block. Keep the
   boundary-flux transport upstream (unchanged).
3. **Config:** restore `CalibrationConfig.exposure_dispersion` (φ). Recover the old default
   (check `exposure.py` callers @ `c0af5b3^` / old config); start there, treat as a tunable (honor
   the `calibration_no_magic_numbers` memory — discuss before any new constant).
4. **Mass unchanged:** `gdna_prior_count` = raw deconvolved gDNA (no ω scaling).

## 6. Tests + validation
1. **Restore** `tests/calibration/test_exposure.py` from `c0af5b3^` (EB shrinkage: hand-calc ω̂,
   sparse→1, vectorized), adapted to the acyclic per-node API.
2. **Capture-leak regression test** (new): assert the 4:1 capture-on/ss99 gDNA→RNA leak is below a
   threshold (target: well under the current 36%; aim ≤10% then iterate) and pool gDNA within X%
   of 80%. Drive via `bench_calibration.py` or a small fixture.
3. **No-regression:** 500 kb suite (pool ±2%, FP ~0) + the unit + golden suites stay green; regen
   goldens only after the calibration is finalized.
4. **Validate** with `scripts/sim/bench_calibration.py --sim-base <suite> --run --force` on BOTH
   suites (20 Mb 4:1 capture AND 500 kb), reporting pool + per-fragment confusion.

## 7. Git artifacts to re-read in the fresh research phase
- `git show c0af5b3^:src/rigel/calibration/exposure.py`  (EB Gamma posterior)
- `git show c0af5b3^:tests/calibration/test_exposure.py`  (EB tests)
- `git show 266ba70 -- src/rigel/calibration/`            (the Option-A exposure→geometric diff)
- `git show 266ba70^1:src/rigel/calibration/priors.py`    (pre-Option-A exposure-as-length assembly)
- docs/futureprs/: `exposure_prior_robustness_theory.md`, `exposure_prior_noise_floor_theory.md`,
  `exposure_prior_findings_report.md`, `phase6_exposure_efflen_theory.md`
- current: `src/rigel/calibration/derive.py`, `priors.py`, `config.py`; `scripts/sim/bench_calibration.py`

## 8. Open decisions (resolve during implementation)
- `exposure_dispersion` (φ) default + whether per-node `L` should be the FL-aware support or bare
  region length in the IPR (doc says bare L; current uses `gdna_geom_len`).
- The anti-collapse guard form (floor eff-len, clamp ω̂, or blend toward span like the current
  trust-shrink).
- Whether to down-weight nodes by `Var(log ω̂)` (the noise-floor doc) beyond the EB mean shrinkage.
- Re-validate that this does not regress the FP-averse direction at non-capture / 50% gDNA.
