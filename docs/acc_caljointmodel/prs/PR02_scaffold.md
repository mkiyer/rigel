# PR 2 — Calibrator scaffold (types, substrate, placeholder)

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §7 (PR 2).
**Type:** Python-only (no index rebuild, no C++).
**Build required:** no.
**Status:** Draft for critique. Introduces **0 new calibrator tunables**;
net **−1** config knob (removes the obsolete `boundary_split_factor`). Open
magic-number reconciliations are collected in §II.7 for your decision before
coding (per the pause-and-discuss rule).

PR 2 turns the skeletons from PR 1 into the real **data contract** of the
calibrator: the output schema (`CalibrationResult`), the adapter that turns
the raw accumulator payload into the model's sufficient statistics
(`CalibrationSubstrate`), and a *placeholder* `calibrate()` that returns a
valid, trivial result (no contamination inferred). **No inference math lands
here** — the E-step/exposure (PR 4) and M-step/outer loop (PR 5) consume what
this PR builds. The point of a scaffold PR is that everything downstream can
be wired and tested against a real, validated result object before the
deconvolution is written.

---

# Part I — Theoretical design overview

This part is the "why": the generative model the calibrator implements, and
where PR 2's two objects (substrate, result) sit inside it. Full derivations
live in [`01_generative_model.md`](../../caljointmodel/01_generative_model.md)
and [`03_inference.md`](../../caljointmodel/03_inference.md); this is the
operational summary.

## I.1 What calibration is for (G2 / G3 / G4)

Rigel must separate three things that all deposit reads on the genome:
mature **mRNA**, **nascent RNA**, and **genomic-DNA contamination (gDNA)**.
gDNA is unstranded, roughly uniform per base, and produces no spliced reads;
RNA is stranded, concentrated on exons, and produces spliced reads. The
calibrator's job is to recover, library-wide and per-region, *how much of the
observed signal is gDNA* so the downstream locus EM gets an informed prior.

Three goals, each an output of this module:

- **G2 — contained-fragment deconvolution.** For each region, split the
  contained fragment mass into gDNA vs RNA.
- **G3 — boundary-flux deconvolution.** Do the same for fragments crossing
  region boundaries (where hybrid-capture enrichment spills gDNA out of
  captured exons — boundaries are where that contamination is *measured*).
- **G4 — per-region exposure.** Recover a per-region gDNA exposure factor
  `ω_r` (local accessibility/copy-number relative to the library mean), with
  an uncertainty, so noisy regions are down-weighted downstream.

## I.2 Per-region sufficient statistics, not per-fragment

The calibrator never iterates fragments. The downstream locus EM already does
per-fragment work at higher precision; re-doing it here would be wasteful and
would block the natural overdispersion treatment that only exists at the
region scale. Instead the model is stated over **per-region (and
per-boundary) sufficient statistics** — three integer counts per region:

| Symbol | Meaning | From the 4-channel payload |
|---|---|---|
| `n_u` | unspliced fragment count | `ch0 + ch1` |
| `n_s` | spliced fragment count | `ch2 + ch3` |
| `k_+` | sense-strand count among unspliced | `ch0` (+ region) or `ch1` (− region) |

where the accumulator channel encoding is `ch = (spliced?2:0) + (strand_pos?0:1)`:

```
ch0 = unspliced & strand+      ch2 = spliced & strand+
ch1 = unspliced & strand−      ch3 = spliced & strand−
```

These three numbers per region are *exactly* what an NB-count + BB-strand +
deterministic-spliced model needs. Everything else is wasted resolution.

## I.3 Two latents, five library hyperparameters

Per region the model has two latents:
- `π_g[r]` — the gDNA mixing proportion (P(fragment is gDNA | region `r`));
- `ω_r` — the gDNA exposure factor (multiplicative on the library density).

Five **library-wide** hyperparameters are fitted (PR 5):

| Param | Role |
|---|---|
| `ρ_0` | library mean gDNA density (fragments per bp) |
| `φ` | overdispersion — unifies NB count dispersion **and** the Gamma exposure-prior spread |
| `ρ_d_bb` | gDNA strand Beta-Binomial dispersion (mean **fixed** `κ_d = 0.5` — gDNA is unstranded by construction) |
| `ρ_r_bb` | RNA strand Beta-Binomial dispersion (mean `κ_rna` from the strand-balance fit, PR 3) |
| `ε_s` | gDNA splice-artifact rate (failsafe; ~`1e-4`–`1e-3`) |

`κ_d = 0.5` is not a free parameter (biology), and the RNA strand mean
`κ_rna` is not fitted in the deconvolution — it is the strand-balance model's
job (PR 3). That leaves only the *dispersions* and the density to fit, which
is what keeps the parameter count tiny.

## I.4 Three channels combine into one per-region split

The E-step (PR 4) combines three independent likelihood channels by adding
log-Bayes-factors:

1. **Count (Negative-Binomial).** Is `n_u` consistent with the gDNA-only mean
   `μ_g = ω_r ρ_0 L_eff`, or does it exceed what gDNA can explain? Via the
   Gamma-Poisson identity, the same `φ` that sets NB dispersion sets the prior
   spread of `ω_r`, giving a **closed-form Gamma posterior** for exposure (G4)
   — no Newton, no quadrature.
2. **Strand (Beta-Binomial vs Beta-Binomial).** gDNA strand ~ `BB(n_u, 0.5,
   ρ_d_bb)`; RNA strand ~ `BB(n_u, κ_rna, ρ_r_bb)`. Overdispersion is the
   robustness mechanism: a moderate strand asymmetry on a contaminated region
   no longer reads as decisive RNA (this is the paralog-phantom fix).
3. **Spliced (deterministic).** A spliced read is RNA except for the tiny
   artifact rate `ε_s`, so `n_s` routes almost entirely into RNA mass without
   a Bayes step.

The count + strand log-BFs update `π_g[r]`; the unspliced count is split
`M_g = n_u · π_g`, `M_d = n_u − M_g`; spliced mass is added to RNA. That
yields the G2 masses; the same code on boundary statistics yields G3.

## I.5 Boundaries are first-class — the D1 side-attribution (locked)

A boundary sits between two regions and stores, per channel: a **fractional
`mass_left` / `mass_right`** (the block-side mass contributed by the left /
right region — *disjoint, asymmetric* observations) and **one shared integer
`flux`** (fragment-events crossing). The earlier "½ half-split" of boundary
mass into regions is **abandoned**; the accumulator already partitions mass by
side, so a region receives the natural side-attributed contributions:

```
region r  ←  (its LEFT boundary's right-side mass)   [r is that boundary's right region]
          +  (its RIGHT boundary's left-side mass)   [r is that boundary's left region]
```

The integer `flux` (shared by both neighbours) drives the count/strand
likelihood → statistical power; the fractional `mass_*` is the magnitude that
power deconvolves → density. **Both are used.** This is implemented on top of
the boundary↔region index mapping built in PR 1
(`region_boundary_indices`).

## I.6 Where PR 2 sits: the substrate and the result

The inference (PR 4/5) is pure scalar algebra over arrays. PR 2 builds the two
arrays-objects it reads and writes:

- **`CalibrationSubstrate` — the bridge from raw payload to sufficient
  statistics.** It reduces the 4 channels into `(n_u, n_s, k_+)` and projects
  the boundary arrays onto regions via the D1 side-attribution, exposing
  **three parallel per-region views** — *contained*, *left*, *right* — so the
  E-step runs the identical code three times (doc 03 §2/§6). Each view also
  carries the float **mass magnitude** to be deconvolved (for *contained* this
  equals the integer count; for the boundary views it is the fractional
  `mass_*`). This is the only object that knows the channel encoding and the
  payload topology — the inference stays encoding-agnostic.
- **`CalibrationResult` — the container for the model's outputs.** Per-region
  deconvolved mass (G2/G3), the exposure posterior `(ω, Var[log ω])` (G4), the
  five fitted hyperparameters, and convergence diagnostics. Its
  `__post_init__` enforces the *intrinsic* invariants (non-negativity, valid
  parameter ranges, monotone mass-change) so a malformed result can never
  reach the downstream prior bridge.

The **placeholder `calibrate()`** returns the trivial, valid result —
**no gDNA inferred**: all contained/boundary mass attributed to RNA, unit
exposure `ω = 1`, hyperparameters at their documented initial values. This is
the "identity" baseline: it lets PR 3's strand model and (eventually) the
`quant_from_buffer` rewrite wire against a real `CalibrationResult` before the
deconvolution exists. PR 4 replaces the body; the schema and substrate do not
change.

---

# Part II — Implementation plan

## II.1 Target surface (end of PR 2)

```
src/rigel/calibration/
  result.py      # CalibrationResult: real schema (doc 04 §5) + __post_init__
  substrate.py   # CalibrationSubstrate + SubstrateView + from_payload (real adapter)
  calibrate.py   # calibrate(): placeholder returning the no-gDNA result
  config.py?     # (in rigel/config.py) CalibrationConfig reconciled (see II.4)
  __init__.py    # add SubstrateView to re-exports
```

No new modules (uses PR 1's `region_arrays.py` + `signature.py` + `errors.py`).
The package stays at 10 calibration modules.

## II.2 — T1: `CalibrationResult` (real schema)

Replace the empty placeholder with doc 04 §5, frozen `slots=True`:

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    # G2 — contained mass (float64[R])
    mass_g_contained, mass_d_contained
    # G3 — boundary mass attributed to each region (float64[R])
    mass_g_left, mass_d_left, mass_g_right, mass_d_right
    # G4 — exposure posterior (float64[R])
    omega, log_omega_var
    # library hyperparameters
    rho_0, phi, rho_d_bb, rho_r_bb, eps_s            # floats
    # convergence diagnostics
    n_iterations: int; converged: bool
    mass_change_history: np.ndarray                  # float64[n_iterations]
    # provenance
    n_regions: int; config: CalibrationConfig
```

`__post_init__` enforces only the **intrinsic** invariants — those checkable
without the substrate:

- every `mass_*` / `omega` / `log_omega_var` is `float64`, length `n_regions`,
  finite, and `≥ 0` (mass) / `> 0` (`omega`, `log_omega_var`);
- `rho_0 > 0`, `phi > 0`, and `0 < rho_d_bb, rho_r_bb, eps_s < 1`;
- `mass_change_history` has length `n_iterations` and is non-increasing within
  a `1e-9` scale tolerance.

> **Deliberate scope note.** Doc 04 §5.1 also lists *mass conservation*
> (`mass_g + mass_d == fragment count`). The result does **not** carry the
> substrate, so it cannot self-check that here; conservation is verified in
> `calibrate()` and in the substrate↔result conservation test (II.6). Keeping
> `__post_init__` to intrinsic checks avoids smuggling the substrate into the
> result type.

## II.3 — T2: `CalibrationSubstrate` (the adapter)

A uniform per-view container plus the top-level substrate:

```python
@dataclass(frozen=True, slots=True)
class SubstrateView:
    n_unspliced: np.ndarray   # int64[R]  — count channel (drives NB + BB)
    n_spliced:   np.ndarray   # int64[R]
    k_plus:      np.ndarray   # int64[R]  — sense among unspliced
    mass_unspliced: np.ndarray  # float64[R] — magnitude to deconvolve
    mass_spliced:   np.ndarray  # float64[R]

@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    n_regions: int
    L_eff:    np.ndarray      # float64[R] — v1: physical region length (region_size_bp)
    ts_class: np.ndarray      # int8[R]    — region transcript-strand class
    contained: SubstrateView
    left:      SubstrateView
    right:     SubstrateView

    @classmethod
    def from_payload(cls, payload: AccumulatorPayload,
                     region_arrays: RegionArrays) -> "CalibrationSubstrate": ...
```

`from_payload` steps:

1. **Validate alignment** (doc 04 §4.1): `region_arrays.n_regions ==
   payload.r_total` and `ref_region_offsets` agree, else
   `CalibrationSubstrateError`. (Reuses the PR 1 guard logic.)
2. **Contained view** — reduce `payload.region_contained[R,4]`:
   `n_unspliced = ch0+ch1`, `n_spliced = ch2+ch3`,
   `k_plus = sense(ch0 vs ch1 by ts_class)`. `mass_* = n_*` (cast float).
3. **Boundary projection** — `lb, rb = region_boundary_indices(
   payload.ref_region_offsets, payload.ref_boundary_offsets)` (PR 1).
   - **left view** of region `r` ← boundary `lb[r]` (r is its *right* region):
     counts from `boundary_flux[lb]`, mass from `boundary_mass_right[lb]`.
   - **right view** of region `r` ← boundary `rb[r]` (r is its *left* region):
     counts from `boundary_flux[rb]`, mass from `boundary_mass_left[rb]`.
   - reduce each the same way (`n_u/n_s/k_+` from flux; `mass_*` from the
     fractional side mass), orienting `k_plus` by **region `r`'s** `ts_class`.
4. `L_eff = region_arrays.region_size_bp`; `ts_class = region_arrays.ts_class`.

**Strand orientation (Q2).** `sense` uses `ts_class[r]`: `+` → `ch0`/`ch2`,
`−` → `ch1`/`ch3`. For `TS_NONE` / `TS_AMBIG` regions there is no defined
sense; per Q2 we assign a **fixed arbitrary** convention (sense = `strand+`).
This is harmless: such regions carry `κ_rna ≈ 0.5` / neutral strand evidence
downstream, so the choice cannot manufacture RNA. (Not a tunable.)

## II.4 — T3: `CalibrationConfig` reconciliation

`CalibrationConfig` already exists (`rigel/config.py`) but its values drifted
from the spec, and one field is now obsolete. **These are existing-constant
reconciliations, not new tunables — but they are magic-number changes, so they
are listed for your decision in §II.7 and applied only once you confirm.**

| Field | Current code | Spec (doc 03 §8 / 04 §3) | Proposed |
|---|---|---|---|
| `max_outer_iterations` | 50 | 25 | **25** |
| `mass_rel_tol` | 1e-4 | 1e-4 | unchanged |
| `phi_floor` | 1e-9 | 1e-6 | **1e-6** |
| `boundary_split_factor` | 1.0 | 0.5 | **remove** (D1 kills the half-split) |

Removing `boundary_split_factor` requires updating its `__post_init__` guard
and any reference (none outside config after D1).

## II.5 — T4: placeholder `calibrate()`

Real signature (already wired in PR 1), real return — trivial result:

```python
def calibrate(*, payload, region_arrays, strand_model, config) -> CalibrationResult:
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    R = substrate.n_regions
    zeros, ones = np.zeros(R), np.ones(R)
    return CalibrationResult(
        mass_g_contained=zeros, mass_d_contained=substrate.contained.mass_unspliced
                                                + substrate.contained.mass_spliced,
        mass_g_left=zeros,  mass_d_left=substrate.left.mass_unspliced
                                        + substrate.left.mass_spliced,
        mass_g_right=zeros, mass_d_right=substrate.right.mass_unspliced
                                         + substrate.right.mass_spliced,
        omega=ones, log_omega_var=np.full(R, config.phi),   # prior variance
        rho_0=<init>, phi=<init>, rho_d_bb=0.01, rho_r_bb=0.01, eps_s=1e-3,
        n_iterations=0, converged=True, mass_change_history=np.empty(0),
        n_regions=R, config=config,
    )
```

i.e. **no gDNA inferred**: all mass → RNA, `ω = 1`, hyperparameters at the
doc 03 §7 initial values. This satisfies conservation trivially (`M_g = 0`,
`M_d = total`). `rho_0`/`phi` use the doc 03 §7 data-driven initializers (a
closed-form over the substrate, not constants).

> The placeholder consuming `strand_model` is intentional even though it does
> not use it yet — it keeps the live call site (PR 1) and signature exercised.

**`run_pipeline` follow-through (required).** In PR 1, `run_pipeline` did
`calibrate(...)` then `raise AssertionError("unreachable")` — fine while
`calibrate` *raised*. Now that it *returns*, that `AssertionError` would fire
as a **new** failure mode. PR 2 captures the result and replaces the assertion
with the PR 6 boundary:

```python
    try:
        calibration = calibrate(payload=…, region_arrays=…,
                                strand_model=strand_models, config=config.calibration)
    finally:
        buffer.cleanup()
    # payload → calibrate → priors → EM (quant_from_buffer) is rewritten in PR 6.
    raise NotImplementedError(
        "run_pipeline post-calibration wiring (quant_from_buffer) lands in PR 6 …")
```

So the pipeline still stops with `NotImplementedError` (same mode as PR 1),
but only **after** `calibrate()` has run and produced a valid result — no
`AssertionError`, no new mode.

## II.6 — T5: tests (`tests/calibration/`)

- `test_result_schema.py` — `__post_init__` accepts a valid result; rejects
  negative mass, `omega ≤ 0`, out-of-range hyperparameters, increasing
  `mass_change_history`, length mismatches.
- `test_substrate.py` — on a hand-built tiny payload + `RegionArrays`:
  channel reductions (`n_u/n_s/k_+`) correct for `+`, `−`, and `NONE`/`AMBIG`
  regions; boundary projection picks the correct side (`mass_right` of the
  left boundary, `mass_left` of the right boundary) and the correct shared
  flux; terminals contribute zero; misaligned payload → `CalibrationSubstrateError`.
- `test_substrate_conservation.py` — `from_payload` over a real sim index
  (mini scenario): per region, contained `mass_u + mass_s` equals the raw
  contained fragment total; boundary mass attributed to regions equals the
  payload boundary mass (no double-count, no loss across ref seams).
- `test_calibrate_placeholder.py` — `calibrate(...)` returns a valid
  `CalibrationResult` with `mass_g ≡ 0`, `mass_d == substrate totals`,
  `omega ≡ 1`, `converged`, `n_regions == payload.r_total`.
- `test_package_surface.py` — extend for `SubstrateView`; `calibrate` no
  longer raises `NotImplementedError`.

## II.7 Open decisions (magic numbers — pause & discuss)

Per the hard rule, nothing below is coded until you confirm:

1. **`phi_floor` 1e-9 → 1e-6?** The spec (doc 03 §8 `_PHI_FLOOR`, doc 04 §3)
   says `1e-6`; the live code has an unexplained `1e-9`. Proposed: adopt the
   spec `1e-6`.
2. **`max_outer_iterations` 50 → 25?** Spec says 25 (converges in 3–10).
   Proposed: 25.
3. **Remove `boundary_split_factor`?** D1 abandoned the half-split, so this
   knob is dead. Proposed: remove it from `CalibrationConfig`. (Confirm you do
   not want to keep a side-attribution weight knob instead — the D1 model has
   none.)
4. **Placeholder semantics = "no gDNA" (all-RNA, ω=1).** Confirm this is the
   intended reading of the plan's "zero-mass/unit-exposure" (zero *gDNA* mass,
   not zero total mass — the latter would violate conservation).
5. **`CalibrationResult.__post_init__` intrinsic-only** (conservation checked
   in `calibrate`/tests, not the dataclass). Confirm.

## II.8 Acceptance gate

- [ ] `import rigel`; `pytest --collect-only tests/` clean.
- [ ] New `tests/calibration/` pass; PR 1 tests still green.
- [ ] `calibrate(...)` returns a schema-valid `CalibrationResult` on a real
      sim index; substrate conservation holds across ref seams.
- [ ] Full suite: failure **mode** stays `NotImplementedError` — now from the
      post-calibrate `quant_from_buffer`/PR 6 boundary, reached only *after*
      `calibrate()` returns a valid result. Count ≈ unchanged (end-to-end
      `run_pipeline` tests still stop pre-EM); **no new modes** (no
      `AssertionError` from the old "unreachable" guard).
- [ ] `ruff check src/ tests/` clean.
- [ ] Magic-number budget: 0 new tunables (net −1 config knob).

## II.9 Rollback

Revert the PR. `calibrate()` returns to the PR 1 `NotImplementedError` stub;
no schema persisted to disk, so nothing else changes.
