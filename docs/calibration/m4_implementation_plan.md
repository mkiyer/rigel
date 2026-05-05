# M4 — Global gDNA Densities + κ Implementation Plan

**Predecessor:** [calibration_v6_plan.md §M4](calibration_v6_plan.md) and
[m3_implementation_plan.md](m3_implementation_plan.md) (M3 cleanup must
land first so the accumulator surface is stable).

**Goal:** convert the M3 `CalibrationScanPayload` plus the M2 region
table into the three locked global gDNA densities
(`ρ̂_INTERGENIC`, `ρ̂_INTRON`, `ρ̂_EXON-INTRON`) and one negative-binomial
overdispersion estimate `κ_t` per density.  All three densities are
fragments-per-bp scalars; `κ` is a unitless shape parameter consumed by
the locoregional EB shrinkage in M6.

**Non-goal:** any per-`Locus` or per-`MultiLocus` quantity.  M4 produces
exactly one `GlobalDensityTable`, which is the input to M6.

---

## 1. Locked formulas (carried verbatim from §2.6)

For each density type `t ∈ {INTERGENIC, INTRON, EXON-INTRON}`:

$$\hat\rho_t \;=\; \frac{N_t}{D_t}, \qquad \hat\kappa_t \;=\; \frac{\sum_R \mu_R^2}{\sum_R (N_R - \mu_R)^2 - \sum_R \mu_R}$$

where $\mu_R = \hat\rho_t \cdot L_{\mathrm{eff}}(R)$ and the sums in $\kappa$
are over the same regions used for $\rho$.

**Numerator / denominator (the only definition that matters):**

| `t`            | $N_t$ (numerator) | $D_t$ (denominator) | Region filter |
|----------------|-------------------|---------------------|---------------|
| `INTERGENIC`   | $\sum_R$ `per_region_counts[R, mask=0b100]` | $\sum_R L_{\mathrm{eff}}(R)$ | `region_df.type == INTERGENIC` |
| `INTRON`       | $\sum_R$ `per_region_counts[R, mask=0b010]` | $\sum_R L_{\mathrm{eff}}(R)$ | `region_df.type == INTRON` |
| `EXON-INTRON`  | $\sum_R \big(u_L(R)\cdot 1_L(R) + u_R(R)\cdot 1_R(R)\big)$ | $\sum_R \big(1_L(R) + 1_R(R)\big) \cdot \bar L_{\mathrm{gDNA}}$ | `region_df.type == EXON` AND $1_L \lor 1_R$ |

with $1_L(R) = $ `region_df.boundary_flux_left[R]` (bool), likewise $1_R$.

**$L_{\mathrm{eff}}$ — overlap geometry, locked:**

$$L_{\mathrm{eff}}(R) = (\mathrm{end} - \mathrm{start}) + \bar L_{\mathrm{gDNA}} - 1$$

This is the **only** geometry permitted in the calibration codebase.  The
M4 module enforces it at the source via a single helper
`_l_eff_overlap(span_bp, fl_mean)` so M6 cannot drift.

**Fallback:** any density with $D_t = 0$ → $\hat\rho_t = 0$ AND
$\hat\kappa_t = $ `KAPPA_DEFAULT`.

---

## 2. Module layout

Two new modules, both pure Python (numpy + pandas, no native code):

```
src/rigel/calibration/
├── density_global.py    # GlobalDensityTable + compute_global_densities
├── _kappa.py            # KappaEstimate + estimate_kappa
└── ...
```

Both are pure functions over the M2 `region_df` and the M3
`CalibrationScanPayload`.  No I/O.  No side effects.

---

## 3. `density_global.py`

### 3.1 Public surface

```python
@dataclass(frozen=True, slots=True)
class GlobalGdnaDensity:
    """One density estimate for one mask category."""
    type: Literal["INTERGENIC", "INTRON", "EXON-INTRON"]
    rho: float                  # fragments / bp
    n_fragments: int            # numerator
    eff_length_bp: float        # denominator
    n_regions_used: int         # rows that contributed > 0 to denominator
    kappa: KappaEstimate        # NB overdispersion + fallback diagnostic


@dataclass(frozen=True, slots=True)
class GlobalDensityTable:
    """The complete global-density block of CalibrationResult."""
    intergenic:    GlobalGdnaDensity
    intron:        GlobalGdnaDensity
    exon_intron:   GlobalGdnaDensity
    gdna_fl_mean:  float        # echo of the input; carried for downstream

    def to_summary_dict(self) -> dict[str, dict[str, float]]: ...
```

### 3.2 Orchestrator

```python
def compute_global_densities(
    region_df: pd.DataFrame,        # the M2 region table
    payload:   CalibrationScanPayload,
    *,
    gdna_fl_mean: float,            # from M7 PoolFLModels.gdna_fl_model.mean
) -> GlobalDensityTable:
    ...
```

**Invariants checked at the top:**

* `len(region_df) == payload.per_region_counts.shape[0]` (else: hard
  raise with the canonical "Rebuild the index" suffix).
* `gdna_fl_mean > 0` (raise with actionable message — caller bug).
* `payload.u_left.shape == payload.u_right.shape == (len(region_df),)`.

### 3.3 Vectorized implementation

The implementation is three numpy reductions, no loops:

```python
EXON       = int(RegionType.EXON)
INTRON     = int(RegionType.INTRON)
INTERGENIC = int(RegionType.INTERGENIC)

types = region_df["type"].to_numpy()
spans = (region_df["end"].to_numpy() - region_df["start"].to_numpy())
leff  = spans + gdna_fl_mean - 1.0   # overlap geometry, locked

prc = payload.per_region_counts        # (R, 8), int64

# --- INTERGENIC density ---
m_ig = types == INTERGENIC
n_ig = int(prc[m_ig, mask.INTERGENIC].sum())
d_ig = float(leff[m_ig].sum())
rho_ig = n_ig / d_ig if d_ig > 0 else 0.0

# --- INTRON density (analogous, mask.INTRON) ---

# --- EXON-INTRON density: capture-aware boundary flux ---
m_ex = types == EXON
fl   = region_df["boundary_flux_left"].to_numpy()
fr   = region_df["boundary_flux_right"].to_numpy()
sides = fl.astype(np.int8) + fr.astype(np.int8)        # 0, 1, or 2
eligible = m_ex & (sides > 0)
num_ex = int(
    (payload.u_left  * fl)[eligible].sum() +
    (payload.u_right * fr)[eligible].sum()
)
den_ex = float((sides[eligible] * gdna_fl_mean).sum())
rho_ex = num_ex / den_ex if den_ex > 0 else 0.0
```

Each density's $\kappa$ comes from `estimate_kappa(per_region_counts,
per_region_eff_lengths, rho_hat)` over the same filtered region set.

For EXON-INTRON the per-region count is $u_L \cdot 1_L + u_R \cdot 1_R$
and the per-region effective length is $(1_L + 1_R) \cdot \bar L_{\mathrm{gDNA}}$.

### 3.4 Algorithmic complexity

Three numpy passes over `(R,)` arrays.  R is on the order of `n_genes ×
n_exons_per_gene ≈ 200_000` for human; total work < 1 ms.  No
materialization of (R, 8) sub-arrays beyond `prc[:, single_column]`
which is a strided view.

---

## 4. `_kappa.py`

### 4.1 Public surface

```python
KAPPA_DEFAULT       = 100.0
KAPPA_MIN           = 1.0
KAPPA_MAX           = 1.0e6
MIN_REGIONS_FOR_MOM = 5


@dataclass(frozen=True, slots=True)
class KappaEstimate:
    value: float            # always finite, always in [KAPPA_MIN, KAPPA_MAX]
    n_regions: int          # how many regions contributed to MoM
    fallback_used: bool     # True ⇒ value == KAPPA_DEFAULT
    fallback_reason: str    # "" if fallback_used is False


def estimate_kappa(
    counts:      np.ndarray,    # (R,) int64 per-region observed counts
    eff_lengths: np.ndarray,    # (R,) float64 per-region L_eff
    rho_hat:     float,         # the global density estimate
) -> KappaEstimate:
    ...
```

### 4.2 Implementation

```python
def estimate_kappa(counts, eff_lengths, rho_hat):
    n = int(counts.size)
    if n < MIN_REGIONS_FOR_MOM:
        return _fallback(n, "n_regions < MIN_REGIONS_FOR_MOM")
    if rho_hat <= 0.0:
        return _fallback(n, "rho_hat <= 0")
    if not np.any(eff_lengths > 0):
        return _fallback(n, "all eff_lengths == 0")

    mu = rho_hat * eff_lengths
    excess = float(np.sum((counts - mu) ** 2) - np.sum(mu))
    if excess <= 0.0:
        return _fallback(n, "excess variance <= 0 (Poisson or under-dispersed)")

    kappa = float(np.sum(mu ** 2) / excess)
    kappa = float(np.clip(kappa, KAPPA_MIN, KAPPA_MAX))
    return KappaEstimate(value=kappa, n_regions=n,
                         fallback_used=False, fallback_reason="")


def _fallback(n_regions: int, reason: str) -> KappaEstimate:
    return KappaEstimate(value=KAPPA_DEFAULT, n_regions=n_regions,
                         fallback_used=True, fallback_reason=reason)
```

### 4.3 Why per-region MoM and not a single-pooled estimator

The naïve pooled MoM `κ̂ = μ²/(σ² − μ)` over a single mean μ̄ conflates
true overdispersion with **exposure heterogeneity** (regions of very
different lengths).  The per-region formula above is the correct MoM
under a NB($\mu_R, \kappa$) model with shared $\kappa$ — it weights each
region by its expected mean and recovers κ unbiased even when exposures
span 3–4 orders of magnitude (which they do: short exons vs whole
introns differ by ~10³×).

A future plan can replace MoM with a single-parameter Newton iteration
on the NB log-likelihood if MoM proves unstable on real data.  Tracked
in §8 of the parent plan.

---

## 5. Test plan

Two new test files, both pure unit (no scanner, no fixtures from disk).

### 5.1 `tests/test_density_global.py` (≥ 11 tests)

* **Hand-counted 3-region scenario:** one INTERGENIC, one INTRON, one
  EXON region with hand-set `per_region_counts`, `u_left`, `u_right`,
  flags.  Assert all three densities to within 1e-12.
* **Pure mRNA library:** `payload.fl_hist[exon_only_mask, :]` populated;
  all gDNA-pool masks empty; all three densities < 1e-9.
* **EXON-INTRON with all flags False:** `boundary_flux_left == False AND
  boundary_flux_right == False` everywhere → `rho_ex == 0.0`,
  `n_regions_used == 0`, no division by zero.
* **INTRON-only library:** all counts in mask `0b010`; INTERGENIC = 0,
  INTRON > 0, EXON-INTRON = 0.
* **50/50 RNA + gDNA library:** densities consistent across types
  within a 30% tolerance band (random-seed integration test).
* **Single-side eligible exon:** `bf_left=True, bf_right=False`; assert
  $u_R$ does NOT contribute to the numerator.
* **Empty `region_df` per type:** if no INTRON regions exist at all
  (single-exon-genome test), `intron.rho == 0`, `intron.n_regions_used
  == 0`, no exception.
* **Schema validation:** mismatched `len(region_df)` vs
  `per_region_counts.shape[0]` → `ValueError("Rebuild the index ...")`.
* **`gdna_fl_mean <= 0` rejection.**
* **`L_eff` overlap geometry pinned:** synthesize one INTERGENIC region
  of span `L=100`, `gdna_fl_mean=350`; assert `eff_length_bp == 100 + 350
  - 1 == 449`.  (This is the M6-protection test — if anyone changes the
  formula, this test catches it.)
* **`GlobalDensityTable.to_summary_dict()` round-trip:** all six
  scalars per type emitted with correct keys.

### 5.2 `tests/test_estimate_kappa.py` (≥ 10 tests)

* **Degenerate inputs:**
  * `n_regions < 5` → fallback with reason "n_regions < MIN_REGIONS_FOR_MOM".
  * `rho_hat == 0` → fallback with reason "rho_hat <= 0".
  * All `eff_lengths == 0` → fallback with reason "all eff_lengths == 0".
  * Underdispersed counts (constant) → fallback with reason "excess
    variance <= 0".
* **Heterogeneous-exposure recovery:** synthesize 50 regions with
  $L_{\mathrm{eff}} \sim \mathrm{LogNormal}$, draw $N_R \sim
  \mathrm{NB}(\mu_R, \kappa^*)$, assert `estimate_kappa(...)` recovers
  $\kappa^*$ within ±30% for $\kappa^* \in \{2, 20, 200\}$ at
  `n_regions=500`.
* **Clipping at `KAPPA_MIN`:** synthesize NB with $\kappa^* = 0.5$;
  assert returned value `>= KAPPA_MIN` (clipped).
* **Clipping at `KAPPA_MAX`:** synthesize Poisson-like data (κ→∞);
  assert returned value `<= KAPPA_MAX` (clipped).
* **Determinism:** same input twice → bit-identical output.
* **`fallback_used=False` ⇒ `fallback_reason == ""`** (schema
  invariant; trivial assert).

---

## 6. Cleanup opportunities (collateral wins)

While reading the code for this plan I noticed two pre-existing items
that are worth piggybacking on the M4 commit (they are zero-risk and
share the same review surface):

1. **`_C_BASE_DEFAULT = 10.0` lives in `locus.py`** as a private module
   constant.  M6 needs it as a public, importable name.  Move it to
   `src/rigel/calibration/locus_prior.py` (which doesn't exist yet — M6
   will create it; M4 doesn't touch this).
2. **`kFlBins = 1024` is duplicated** in
   [src/rigel/native/calibration/accumulator.h](src/rigel/native/calibration/accumulator.h) and as
   `(8, 1024)` in [src/rigel/calibration/scan_payload.py](src/rigel/calibration/scan_payload.py).  A future
   widening (e.g., 2048 bins for ultra-long-fragment libraries) needs
   to touch both.  Add a Python module-level `FL_HIST_N_BINS = 1024`
   that mirrors the C++ side and have `scan_payload.py` use it; M4
   ships this in the same commit because `density_global.py` does not
   need to hard-code the shape.

Both are 1–5 LoC.  Defer 1; do 2 in this milestone.

---

## 7. Exit gate

* New tests green: ≥ 11 in `tests/test_density_global.py`, ≥ 10 in
  `tests/test_estimate_kappa.py` (≥ 21 total, beating the v6-plan
  budget of 11).
* `pytest tests/` full suite green (969+21 = 990 passing, the 21
  pre-existing golden-output failures unchanged).
* `density_global.py` ≤ 200 LoC; `_kappa.py` ≤ 100 LoC.
* Public surface exported via `src/rigel/calibration/__init__.py`:
  `compute_global_densities`, `GlobalDensityTable`, `GlobalGdnaDensity`,
  `estimate_kappa`, `KappaEstimate`, `KAPPA_DEFAULT`, `KAPPA_MIN`,
  `KAPPA_MAX`.
* `grep -rn "max(0, .* - .* + 1)" src/rigel/calibration/` returns no
  hits (the containment formula is forbidden in the calibration
  package).

---

## 8. What this milestone deliberately does NOT do

* Does not call into the EM solver.
* Does not consume `region_df["tx_pos_bp"]` / `tx_neg_bp` /
  `exon_pos_bp` / `exon_neg_bp` — those are M6 inputs.
* Does not write to `summary.json` — that wiring is M8b.
* Does not replace any SRD-v1 module.  `density_global.py` and
  `_kappa.py` are additive.
* Does not introduce a per-reference (per-chromosome) EB level —
  explicitly forbidden by the parent plan §2.6.1.
