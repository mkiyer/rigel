# SRD v1 Implementation Plan (Consolidated, Final)

**Status:** APPROVED — phases below are the contract for execution.
**Supersedes:** `simplification_overhaul_v1.md`,
`simplification_plan_other_claude.md`, `srd_v1_plan.md`, and the prior
cluster-channel-bearing version of this file (archived as
`srd_v1_implementation_v0_clusters.md.bak`).
Phase 0 audit findings live in [srd_v1_audit.md](srd_v1_audit.md).

This document is the single reference point. It folds in:

- Substantive content from both the user's original plan and the agent's
  revised plan.
- The user-locked decision to **drop the cluster channel** (multimapper
  pollution + exon-boundary effective-length mismatch made it medium-weight
  and brittle).
- The corrected antisense-pool inclusion rule (gated by SS, since at SS=0.5
  antisense exonic is 1000:1 RNA-contaminated for highly expressed
  transcripts and provides zero gDNA enrichment).
- Phase 0 audit findings (locked).
- Module names actually present in the repo.

---

## 1 — Context and goals

Rigel's calibration step (`src/rigel/calibration/`) currently runs a v5
three-channel mixture EM (density + strand + FL) over annotation regions.
It fails on hybrid-capture data because the density channel assumes a
single global gDNA rate `λ_G` while real capture libraries exhibit 100–140×
spatial variation. The fix is not more parameters; it is **changing what we
classify**: stop classifying *regions* with density, classify individual
*fragments* with cheap geometric and (when available) strand signals.

Calibration's only load-bearing output to downstream is the gDNA
fragment-length model, `gdna_fl_model`, consumed by the C++ scorer. SRD v1
preserves that contract and discards everything else.

### 1.1 Design constraints (user-locked)

- **Incredibly simple, blazing fast, library-agnostic.**
- **No cross-region density.**
- **No SS dependence in the pool definition.** Pool is a pure geometric
  construct (EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC). SS is tracked
  as a diagnostic but does not gate or weight pool membership.
- **First-principles fragment pools** derived from genomic geometry.
- **High-confidence fragments only.** We are building a fragment-length
  model; one bad fragment costs more than ten missing ones. No
  marginally-informative fragment classes (e.g. antisense-exonic) enter
  the pool.
- **Graceful degradation** for QC-fail libraries (gDNA≈0, RNA≈0) via a
  shared global-FL empirical-Bayes prior.
- **Delete, don't deprecate.** Pre-1.0 single-repo project; no zombie
  fields on `CalibrationResult`.
- **Rip out region partitioning entirely** (region tables / `_annotate.py` /
  `region_counts` / `fl_table` plumbing). NOTE: `partition.py` is *not*
  region partitioning — it is the per-locus CSR scatter and is renamed to
  `locus_partition.py` for clarity (Phase 3.5).
- **Per-locus prior is rebuilt, not deleted.** v5's region-derived
  `compute_locus_priors` dies. In its place, a Pass-4 step computes a
  per-locus Dirichlet prior from per-fragment SRD posteriors
  (`γ_f = P(gDNA | features)`) aggregated to per-locus `γ_ℓ`. This
  preserves the load-bearing functional role of the v5 prior — a
  persistent M-step bias that keeps gDNA mass concentrated in
  gDNA-rich loci across every EM iteration — without re-introducing
  cross-locus density comparisons. See §2.6.

### 1.2 Headline Phase 0 audit findings (locked)

1. **C++ requires zero source changes.** The scorer reads only
   `gdna_fl_model`. The C++ EM solver reads only per-locus `α_gdna`/`α_rna`
   arrays — agnostic to how they are computed. SRD computes them in
   Python from per-fragment posteriors before invoking the EM (§2.6).
2. **`partition.py` stays but is renamed** to `locus_partition.py`. It is
   the per-locus CSR scatter for the per-locus EM, *not* region
   partitioning.
3. **`compute_locus_priors` is replaced, not removed wholesale.** The v5
   region-aggregating body dies; a new ~30-line function in
   `locus.py` computes the same `(α_gdna, α_rna)` arrays from already-
   scored per-locus partition data. Deleting the v5 body cleanly severs
   the entire region data flow; the new body has no region dependency.
4. **`FragmentLengthModels.global_model` already exists** and
   `FragmentLengthModel.finalize(prior_counts=..., prior_ess=...)` already
   implements EB shrinkage. `_fl_empirical_bayes.py` is a thin wrapper.
5. **`region_id` is NOT in `FragmentBuffer`.** Region IDs only live in the
   optional `fl_table` returned by the scanner. Removing region machinery
   has zero buffer-schema impact.

### 1.3 What changed from the previous (cluster-channel) plan

| Removed | Reason |
|---|---|
| Pass 1 overlap-cluster builder (`_cluster.py`) | Multimapper pollution + exon-boundary effective-length mismatch |
| Per-cluster γ closed form | Same |
| Fit A (strand-based weighted FL) (`_fit_a_strand.py`) | Same |
| Fit A vs Fit B fusion + KL divergence | Single pool needs no fusion |
| Cluster-channel CalibrationResult fields | No cluster channel |
| Knobs `min_cluster_N`, `gamma_min`, `kl_tolerance` | No cluster channel |

| Added / changed vs. user's "no cluster" plan | Reason |
|---|---|
| `UNSPLICED_ANTISENSE_EXONIC` excluded from the pool entirely | Building a fragment-length model demands high-confidence fragments. At any imperfect SS, antisense-exonic is heavily contaminated by RNA from highly-expressed loci (worked example: at SS=0.5 a 1000×-RNA transcript dumps ~500 RNA fragments per ~0.5 gDNA fragment into the pool). Tracked as diagnostic only. |

---

## 2 — SRD v1 design

### 2.1 Signals

| Signal | Strength | Library | Used in |
|---|---|---|---|
| Spliced N-junction CIGAR | Hard — 100% RNA | All | SPLICED category (trains `RNA_FL`) |
| Exon-fit geometry (with bp tolerance) | Strong — excludes mature mRNA | All | Pool (EXON_INCOMPATIBLE) |
| No-exon-overlap geometry | Strong — gDNA + nRNA only | All | Pool (INTRONIC, INTERGENIC) |
| Antisense vs tx_strand | Diagnostic only — confirms stranding, surfaces unannotated antisense transcription | Stranded | `n_antisense_exonic` count (NOT in pool) |
| Fragment length | Soft — posterior weight in mixture | All | Pool (mixture decomposition) |

**Explicitly dropped:** cross-region density, Poisson/NegBin rate fits,
LogNormal expression priors, Beta-Binomial dispersion, regional Dirichlet
warm-starts, overlap clusters, per-cluster strand statistics, capture-class
composite Poisson.

### 2.2 Pass 0 — Fragment categorization

Walk `FragmentBuffer` once. For each fragment:

1. **Apply unique-map filter** (`num_hits == 1`). Multimappers excluded
   from category statistics — their true origin is unresolved. Counted
   in `n_multimap_excluded` for diagnostics.
2. **Categorize** into one of seven mutually exclusive `FragmentCategory`
   values (uint8 enum):
   - `SPLICED` — any N-CIGAR (read directly from `splice_type`).
   - `UNSPLICED_SENSE_EXONIC` — fits inside one annotated exon within
     `exon_fit_tolerance_bp`; strand matches a unique tx_strand.
   - `UNSPLICED_ANTISENSE_EXONIC` — fits inside one annotated exon;
     strand opposite to a unique tx_strand.
   - `UNSPLICED_EXONIC_AMBIG` — fits inside an exon but tx_strand is
     undefined: unstranded library, mixed-strand exon overlap (`ambig_strand`),
     or unknown read strand.
   - `EXON_INCOMPATIBLE` — unspliced; does not fit any annotated exon
     within tolerance (overhangs an exon-intron boundary, or is longer than
     any containing exon). Cannot be mature mRNA.
   - `INTRONIC` — unspliced; ≥1 transcript candidate, no exon overlap.
   - `INTERGENIC` — unspliced; no transcript candidates at all.

Implementation notes:
- Reuses C++ scanner outputs (`splice_type`, `t_indices`, `intron_bp[i]`,
  `ambig_strand`, `exon_strand`) — no new geometry code in C++.
- Exon-fit test reduces to `min(intron_bp[i]) <= exon_fit_tolerance_bp`
  via `np.minimum.reduceat` over the CSR `intron_bp` array.
- `exon_fit_tolerance_bp` config knob, default 5 bp.
- All strand assignments inherit the existing R2-flip convention from the
  BAM scanner. No new strand semantics.

### 2.3 Pass 1 — Pool assembly

A single unified pool of gDNA-candidate fragments, defined by genomic
geometry alone:

```
Pool := EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC
```

**Library-agnostic. No SS dependence.** Identical pool definition for
stranded and unstranded libraries.

**Why `UNSPLICED_ANTISENSE_EXONIC` is excluded** (worked example with a
highly-expressed transcript at RNA:gDNA = 1000:1):

| SS | Antisense-exonic pool composition | Signal-to-noise |
|---|---|---|
| 0.50 | ~500 RNA + 0.5 gDNA | 1000:1 RNA-dominated |
| 0.70 | ~150 RNA + 0.5 gDNA | 300:1 RNA-dominated |
| 0.90 | ~50 RNA + 0.5 gDNA | 100:1 RNA-dominated |
| 0.95 | ~25 RNA + 0.5 gDNA | 50:1 RNA-dominated |
| 1.00 | 0 RNA + 0.5 gDNA | clean (but rare in practice) |

At every realistic SS the antisense-exonic pool is RNA-dominated. The
1-D mixture would still recover π correctly (it's well-posed), but
`gDNA_FL` would be dominated by soft-residual noise from mis-assigned
RNA fragments. Since calibration's only product is a fragment-length
*shape*, any RNA contamination directly degrades the output.

**`UNSPLICED_SENSE_EXONIC` and `UNSPLICED_EXONIC_AMBIG`** are also
excluded (sense is dominated by mature mRNA; ambig is unresolvable).

**Diagnostic tracking.** Per-category counts are still emitted in
`category_counts`; downstream tooling can use `n_antisense_exonic` to
verify stranding, surface unannotated antisense transcription, etc.
But none of these counts feed the pool.

### 2.4 Pass 2 — 1-D mixture EM (the only iterative step)

Under the assumption `nRNA_FL ≈ RNA_FL` (both reflect post-library-prep
fragmentation of long precursor molecules), fit:

$$
\text{pool\_FL}(l) = \pi \cdot \text{gDNA\_FL}(l) + (1-\pi) \cdot \text{RNA\_FL}(l)
$$

- `RNA_FL(l)` is **fixed** (from SPLICED fragments, smoothed by global-FL
  prior — see Pass 3).
- `gDNA_FL(l)` is free as a Dirichlet-smoothed histogram.
- `π` is a scalar in [0, 1]; initialized to 0.1.
- Soft-EM, one E-step + one M-step per iteration; <30 iterations to
  convergence; tolerance on `|Δπ|` (default 1e-4).
- Convex given fixed `RNA_FL`.

**Skip conditions:** none. The mixture EM is convex given fixed
`RNA_FL` and runs unconditionally. Small or empty pools converge
immediately to π → 0 (since the smoothed-pool initial gDNA component
matches RNA), and EB shrinkage in Pass 3 collapses `gDNA_FL` to global
FL shape. No `min_pool_size` gate is required; robustness is
guaranteed by the prior, not by a threshold.

### 2.5 Pass 3 — Empirical-Bayes smoothing and graceful fallback

**Global FL.** Histogram over *all uniquely-mapped* fragments (unfiltered
by category). This is the prior anchor for all downstream FL models.
Reuses `FragmentLengthModels.global_model`.

**`RNA_FL`.** Histogram of SPLICED fragment lengths, Dirichlet-smoothed
toward global FL with ESS `fl_prior_ess` (default 500). Near-zero SPLICED
count (QC-fail) → collapses to global FL shape.

**`gDNA_FL`.** Histogram returned by Pass 2, scaled to `π · n_pool`,
Dirichlet-smoothed toward global FL with the same ESS. If Pass 2 was
skipped or failed, set to global FL shape.

**`gdna_fl_model_quality` enum** (purely behavioral; no pool-size gate):
- `"good"` — mixture converged AND `π > 0.02` (meaningful gDNA signal).
- `"weak"` — converged but `π ≤ 0.02` (no detectable gDNA, or pool so
  small that EB shrunk π to ~0).
- `"fallback"` — mixture failed to converge (degenerate input). `gDNA_FL`
  falls back to global FL shape; per-locus FL LLR ≈ 0 everywhere;
  downstream EM relies on splice + strand alone.

**QC-fail paths:**
- **gDNA ≈ 0**: `π → 0`, quality = "weak" or "fallback". Correct: no gDNA
  signal to find.
- **RNA ≈ 0** (QC-fail library): SPLICED ≈ 0 → `RNA_FL` collapses to
  global FL → mixture is degenerate → `gDNA_FL` falls back. Loud `WARN`
  log emitted (threshold: SPLICED < 100). Tool completes; user sees the
  signal.

### 2.6 Pass 4 — Per-locus Dirichlet prior from per-fragment posteriors

**Purpose.** Restore the load-bearing functional role of v5's
`compute_locus_priors` — a persistent M-step Dirichlet pseudocount that
biases each per-locus EM toward its data-supported gDNA share — without
any cross-locus density comparison or region table.

**Why this is needed.** The C++ EM solver consumes per-locus `α_gdna` /
`α_rna` arrays in two distinct ways: as a warm-start for `θ` at iter 0
(visible) and as a Dirichlet pseudocount in every M-step (persistent
bias). When the calibrated FL channel is weak — e.g. small libraries
where `RNA_FL` and `gDNA_FL` both collapse to global FL — zero priors
leave the EM free to absorb gDNA-looking unspliced fragments into any
RNA component in the locus, including silent negative-control
transcripts. The persistent bias is what prevented this in v5.

**Algorithm (single O(n_units) reduction over per-locus partitions, run
immediately before the per-locus EM):**

For every fragment `f` in every per-locus partition `ℓ`, compute the
per-fragment gDNA posterior under a 2-class model whose prior is the
library-wide `π_pool` from Pass 2 and whose likelihoods are the same
calibrated `gDNA_FL` and `RNA_FL` already used for scoring:

$$
\gamma_f \;=\; \sigma\!\Big(
  \log\tfrac{\pi_{\text{pool}}}{1-\pi_{\text{pool}}}
  + \log p_{\text{gDNA}}(f)
  - \log\max_{t} p_{\text{RNA}}(f \mid t)
\Big)
$$

The gDNA log-likelihood `ℓog p_gDNA(f)` already encodes splice
incompatibility (gDNA splice penalty applied during scoring), so SPLICED
fragments collapse to `γ_f ≈ 0` automatically. Aggregate to per-locus:

$$
\gamma_\ell \;=\; \frac{1}{|U_\ell|}\sum_{f \in U_\ell} \gamma_f, \quad
\alpha_{\text{gDNA}}[\ell] = c_{\text{base}}\,\gamma_\ell, \quad
\alpha_{\text{RNA}}[\ell] = c_{\text{base}}\,(1 - \gamma_\ell)
$$

where `c_base = 5.0` matches the v5 default. The data **already lives
in the per-locus partitions** (`gdna_log_liks` per-unit, `log_liks`
per-candidate); no extra buffer pass is required.

**Inputs.** Per-locus partitions (post-`partition_and_free`), `pi_pool`
from the calibration result.

**Outputs.** `α_gdna[ℓ]`, `α_rna[ℓ]` arrays consumed by the existing C++
EM entry points (no C++ change).

**Properties.**

- *Library-agnostic.* `γ_f` uses only fragment-local features (FL,
  splice, optionally strand). No cross-locus density comparison.
- *Hybrid-capture-immune.* Per-locus aggregation is independent of
  any genome-wide rate.
- *Graceful degradation.* When `gdna_fl_quality == "fallback"`,
  `gDNA_FL` falls back to global FL shape; `gDNA_FL[len]` and
  `RNA_FL[len]` are nearly identical, so `γ_f ≈ π_pool` everywhere and
  the per-locus prior collapses to a uniform `π_pool` baseline. This is
  exactly the "capitalism" behavior the QC-fail case wants.
- *Single knob.* `c_base` controls the prior strength. Kept as an
  internal constant (5.0) per the simplicity constraint; promoted to
  `CalibrationConfig` only if Phase 5 benchmarks demand it.

**Implementation site.** `src/rigel/locus.py` exposes a new function
`compute_locus_priors_from_partitions(partitions, loci, pi_pool, *,
c_base=5.0)`. The pipeline calls it between `partition_and_free` and
`_run_locus_em_partitioned`, replacing the temporary `_zero_locus_priors`
stub introduced during the Phase 3 wiring.

### 2.7 New `CalibrationResult` schema

**Kept:**

```python
@dataclass(frozen=True)
class CalibrationResult:
    # Models
    gdna_fl_model: FragmentLengthModel
    rna_fl_model: FragmentLengthModel
    global_fl_model: FragmentLengthModel
    gdna_fl_quality: Literal["good", "weak", "fallback"]

    # Library-level
    strand_specificity: float

    # Per-category counts (length 7; matches FragmentCategory enum order)
    category_counts: np.ndarray[int64]
    n_multimap_excluded: int

    # Pool diagnostics
    n_pool: int
    pi_pool: float
    mixture_converged: bool
    mixture_iterations: int

    # Config echo
    exon_fit_tolerance_bp: int
    fl_prior_ess: float

    # Free-form (warnings, debug strings)
    extra: dict[str, Any]
```

**Deleted outright** (not deprecated): `region_e_gdna`, `region_n_total`,
`lambda_gdna`, `region_gamma`, `region_gamma_strand`, `mu_R`, `sigma_R`,
`mixing_pi`, `mixing_pi_soft`, `kappa_G`, `kappa_R`, `em_n_iter`,
`em_converged`, `n_eligible`, `n_soft`, `n_spliced_hard`, `lam_G_on`,
`lam_G_off`, `capture_class_mode`. Plus all cluster-related fields from
the prior plan: `n_pool_a`, `pool_a_weight`, `pi_pool_b`,
`fit_ab_kl_divergence`, `n_clusters_total`, `n_clusters_gdna_dominant`,
`fused`, `used_fit_a`, `used_fit_b`.

### 2.8 Why this works where v5 failed

| v5 failure | SRD v1 behavior |
|---|---|
| Single global `λ_G` misspecified under hybrid capture | No density channel. Spatial variation irrelevant. |
| Poisson / LogNormal / Beta-Binomial / κ dispersion | None. Histograms + a single 1-D mixture. |
| 50-iter fused-LLR EM with quadrature | Single 1-D EM, <30 iters. No quadrature. |
| Brittle hyperparameters | Two robust knobs: `exon_fit_tolerance_bp` (5 bp) and `fl_prior_ess` (500). Robustness from EB shrinkage, not from gates. |
| No viable unstranded signal | EXON_INCOMPATIBLE + INTRONIC + INTERGENIC pool: library-agnostic, first-principles. |
| Regional prior load-bearing but misfit | No regional prior. v5's `compute_locus_priors` body deleted. The Pass-4 replacement (§2.6) computes the same `(α_gdna, α_rna)` arrays from per-fragment SRD posteriors aggregated to the per-locus EM's own loci — zero region dependency. |
| Unstable on QC-fail libraries | EB global-FL fallback; graceful, loud, visible. |

### 2.9 Trade-off accepted

Stranded libraries get no antisense-derived gDNA enrichment in the pool.
**Acceptable** because:

- The exon-incompatibility + nonexonic pool is library-agnostic and
  carries the bulk of the gDNA signal.
- Mixing in antisense-exonic at any realistic SS injects RNA noise that
  directly degrades the recovered FL shape — the only thing the
  calibration step produces.
- High-SS libraries do still benefit indirectly: cleaner SPLICED →
  sharper `RNA_FL` → better mixture separation in Pass 2.
- Cluster-channel complexity (per-cluster effective-length correction,
  splice-block accounting, multimap handling) contradicts the "simple and
  fast" constraint.
- v2 may revisit if Phase 5 validation on real data shows the pool is
  starved on stranded libraries.

---

## 3 — Phased execution plan

### Phase 0 — Audit (✅ COMPLETE)

Deliverable: [srd_v1_audit.md](srd_v1_audit.md). Findings locked in §1.2.

### Phase 1 — Park v5 in legacy branch (skipped per user direction)

User authorized free deletion; rollback path preserved by git history alone.

### Phase 2 — Build SRD core (✅ COMPLETE; revised, no cluster channel)

**Files in `src/rigel/calibration/`:**

| File | Purpose |
|---|---|
| `_categorize.py` | Pass 0 categorizer. `FragmentCategory` enum; `categorize_chunk(...)`; `build_t_to_ref_id(...)`. |
| `_fl_mixture.py` | Pass 2 1-D mixture EM. `fit_fl_mixture(...)` returning `FLMixtureResult`. |
| `_fl_empirical_bayes.py` | Pass 3 EB FL builders: `build_global_fl`, `build_rna_fl`, `build_gdna_fl`. Thin wrappers over `FragmentLengthModel.finalize(prior_counts=..., prior_ess=...)`. |
| `_result_srd.py` | New minimal `CalibrationResult` (Phase 3 will replace `_result.py` with this content). |
| `_simple.py` | Orchestrator `calibrate_gdna(...)`. Single buffer walk → unique-map filter → categorize → SS-gated pool → mixture → EB FL models. |

**Tests in `tests/`:**

| File | Coverage |
|---|---|
| `test_categorize.py` | 11 tests. All 7 categories covered, including within-tolerance overhang, beyond-tolerance overhang, multi-candidate fit, ambig via mixed strand, ambig via unknown strand. |
| `test_fl_mixture.py` | 7 tests. Known-mixture recovery, π=0, π=1, identical FLs, empty pool, mass conservation, length mismatch raises. |
| `test_fl_empirical_bayes.py` | 6 tests. Global finalize correctness, RNA collapse to global shape on empty spliced, gDNA fallback shape, distinct gDNA when signal strong, ESS effect on shrinkage, max_size propagation. |

**Decisions baked in:**
- No SS-threshold branch in the orchestrator's control flow.
- Antisense-exonic inclusion gated by continuous `min_ss_for_antisense`.
- `gdna_fl_model_quality` enum emitted.
- `exon_fit_tolerance_bp` default 5; `fl_prior_ess` default 500;
  `min_pool_size` default 500; `min_ss_for_antisense` default 0.85.
- Unique-map filter mandatory.

**Status:** All 24 tests green. v5 still wired into pipeline.

**Outstanding from this phase:** add the SS-gating to `_simple.py` (the
current Phase 2 implementation includes antisense unconditionally; needs
the new gate per §2.3) and add a regression test for the gate.

### Phase 3 — Wire SRD; replace `CalibrationResult`; gut deprecated knobs

**Phase 3.5 — Rename `partition.py` → `locus_partition.py`** (one focused
commit, done first):

- `git mv src/rigel/partition.py src/rigel/locus_partition.py`
- `git mv tests/test_partition.py tests/test_locus_partition.py`
- Update import in `pipeline.py`:
  `from .partition import partition_and_free` → `from .locus_partition import partition_and_free`
- Update docstring to clarify scope ("Per-locus CSR scatter for the
  per-locus EM. Unrelated to region/annotation partitioning.").

**Modified files:**

- `src/rigel/calibration/_result.py` — replace contents per §2.6. Drop
  every v5 field. Keep class name `CalibrationResult`. Delete `_result_srd.py`.
- `src/rigel/calibration/__init__.py` — export `calibrate_gdna` from
  `_simple`. Drop v5 exports.
- `src/rigel/config.py` — remove dead v5 knobs (κ bracket, eligibility
  floor, capture-class flags, μ_R/σ_R seeds). Add `CalibrationConfig`
  dataclass:
  ```python
  @dataclass(frozen=True)
  class CalibrationConfig:
      exon_fit_tolerance_bp: int = 5
      fl_prior_ess: float = 500.0
      max_iter: int = 50
      tol: float = 1e-4
  ```
- `src/rigel/pipeline.py` — replace `calibrate_gdna(...)` call site with
  the new SRD signature. Remove `region_counts` / `fl_table` plumbing in
  `scan_and_buffer`. Remove `_compute_priors` helper. In
  `quant_from_buffer`, pass zero `α_gdna` / `α_rna` arrays to the C++ EM.
- `src/rigel/estimator.py` — remove any reads of `region_e_gdna` /
  `lambda_gdna`.
- `src/rigel/locus.py` — **delete v5 `compute_locus_priors`** body and add
  the SRD `compute_locus_priors_from_partitions(partitions, loci, pi_pool, *, c_base=5.0)`
  per §2.6.
- `src/rigel/native/*` — **no source changes** per Phase 0 audit.

**Tests:**

- `tests/test_calibration_simple.py` — seven scenarios:
  1. **Stranded pure-RNA** (SS=0.95, zero gDNA): `π_pool ≈ 0`, quality
     `"weak"` or `"fallback"`.
  2. **Stranded 20% gDNA** (SS=0.95, distinct gDNA FL): `π_pool ≈ 0.20 ± 0.02`;
     `gDNA_FL` matches injection.
  3. **Unstranded 20% gDNA** (SS=0.5, distinct gDNA FL): `π_pool ≈ 0.20 ± 0.03`;
     `gDNA_FL` matches injection. Critically: same pool definition as
     stranded; verifies library-agnosticism.
  4. **Hybrid-capture-like** (stranded, 100× density variation, uniform
     10% gDNA): SRD unaffected by density; consistent `gDNA_FL`. v5 fails.
  5. **QC-fail near-empty RNA**: graceful fallback, loud warning,
     `quality = "fallback"`.
  6. **Antisense-exclusion regression** (SS=0.5, highly-expressed RNA loci
     synthesizing many antisense-exonic fragments, low real gDNA): `n_pool`
     does NOT include antisense; `π_pool` stays low and accurate; verify
     `category_counts[UNSPLICED_ANTISENSE_EXONIC] > 0` (diagnostic
     tracking still works).
  7. **Schema contract**: `CalibrationResult` fields match new schema; no
     v5 fields present.

**Goldens:** regenerate `tests/golden/` with `pytest --update-golden` after
eyeballing diffs. Splice-anchored loci must not regress; hybrid-capture
loci should improve.

**Exit criteria:** `pytest tests/ -v` green; smoke `rigel quant` runs on a
stranded BAM and an unstranded BAM; eyeball `summary.json`.

### Phase 4 — Bulk deletion

**Move to `src/rigel/_deprecated/calibration_v5/`:**

- `src/rigel/calibration/_em.py`
- `src/rigel/calibration/_stats.py`
- `src/rigel/calibration/_fl_model.py` (old; superseded by `_fl_empirical_bayes.py`)
- `src/rigel/calibration/_calibrate.py`
- `src/rigel/calibration/_annotate.py`

**Delete from `src/rigel/index.py`:**

- `build_region_table` function (lines 525–641).
- `region_df` property (line 713).
- `region_cr` property (lines 725–738).
- `uniform_region_exposures` if unique-consumer (review line 891–916).
- The region build invocation in `TranscriptIndex` constructor.

**Delete from `src/rigel/native/bam_scanner.cpp`:**

- The `fl_table` / `region_id` accumulation paths.
- **KEEP** all per-fragment exon resolution, strand classification, and
  splice-type detection — these are SRD primitives.
- **Recompile** with `pip install --no-build-isolation -e .` and run
  `pytest tests/ -v`.

**Archive to `docs/deprecated/calibration_v5/`:**

- `docs/calibration/strand_channel_theory_and_redesign.md`
- `docs/calibration/strand_channel_implementation_plan.md`
- `docs/calibration/strand_first_hybrid_capture_plan.md`
- `docs/calibration/density_nb_model_plan.md`
- `docs/calibration/gdna_locus_prior_deep_diagnostic.md`
- `docs/calibration/capture_class_density_plan.md`
- `docs/calibration/channel_fusion_hybrid_capture.md`
- `docs/calibration/effective_length_audit.md` (review for mappability content)
- `docs/calibration/gdna_implementation_log.md`
- `docs/calibration/gdna_length_model_options.md`
- `docs/calibration/root_cause_locus_em_nrna_siphon.md`
- `docs/calibration/SESSION_HANDOFF_2026-04-23.md`
- `docs/calibration/simplification_overhaul_v1.md`
- `docs/calibration/simplification_plan_other_claude.md`
- `docs/calibration/srd_v1_implementation_v0_clusters.md.bak`
- `docs/calibration/srd_v1_plan.md`

Keep on `docs/calibration/`: this file (`srd_v1_implementation.md`),
`srd_v1_audit.md`, plus the new `srd_v1_results.md` from Phase 5.

**Archive scripts:**
- `scripts/calibration/stress_mini/` (uses v5 API).
- `scripts/benchmark/compare_vcap_optionb.py` — either update to read new
  SRD diagnostics or archive.

**Verify:**
- `pytest tests/ -v` green.
- `ruff check src/ tests/` clean.
- `grep -r "from rigel._deprecated" src/` returns nothing.
- `src/rigel/calibration/` LOC dropped from ~2100 to ~300–500.

**Exit criteria:** all of the above; commit pushed.

### Phase 5 — Integration validation

**Tasks:**

1. Run all 13 Armis2 conditions:
   `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml`
2. Analyze vs v5 baseline:
   `python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml -o results/srd_v1_report`
3. **Pass criteria:**
   - Stranded pure-RNA: relative-error change ≤ +5% vs v5.
   - Stranded high-contamination: improvement vs v5.
   - Unstranded with gDNA: SRD produces non-trivial gDNA separation
     (v5 was essentially untested here).
   - Hybrid-capture-like synthetic: SRD recovers near-uniform gDNA;
     `quant.gdna_total / strand_inferred_gdna < 1.1×`.
   - VCaP titration (if present): same metric.
4. **Profile:**
   - Pass 0 (categorize) < 1 s / M fragments.
   - Mixture fit < 100 ms.
   - Total calibration < 25% of v5 wallclock.

**Deliverable:** `docs/calibration/srd_v1_results.md`.

**Decision gate for v2 (cluster channel):** if Phase 5 shows
`gdna_fl_model_quality = "good"` on most benchmarks AND moderate-SS
libraries are within target, **v2 does not add the cluster channel**. Add
it only if specific benchmarks reveal a quality gap SRD v1 cannot close
via simpler tuning.

### Phase 6 — Cleanup & docs

- Update `CLAUDE.md`, `README.md`, `docs/MANUAL.md`, `docs/METHODS.md` to
  describe SRD v1; remove v5 references.
- Bump version. `CHANGELOG.md` entry for breaking changes:
  - `CalibrationResult` schema replaced.
  - `compute_locus_priors` removed.
  - `partition.py` renamed to `locus_partition.py`.
  - v5 config knobs removed; new `CalibrationConfig` knobs added.
- Archive v5 docs in `docs/deprecated/` (per Phase 4 list).
- **Phase 4b carryover (deferred from Phase 4):** delete now-dead region machinery:
  - `src/rigel/index.py`: drop `build_region_table`, `region_df` property,
    `region_cr` property, and the constructor invocation that produces
    `regions.feather`. Index format bump (or backwards-compat shim that
    treats missing `regions.feather` as `region_df = None`).
  - `src/rigel/native/bam_scanner.cpp`: remove `fl_table` and `region_id`
    accumulation paths (Python already ignores `region_evidence` post-SRD).
    Recompile via `pip install --no-build-isolation -e .`.
  - `src/rigel/pipeline.py`: remove the dead `build_region_index` C++ call
    (currently around lines 243-251).
  - `src/rigel/mappability.py`: delete `uniform_region_exposures` and any
    region-exposure helpers left without consumers after the above.
  - `tests/test_regions.py`: delete (pure region-structure tests).
  - `tests/test_mappability.py`: prune region-exposure tests; keep any
    transcript-level mappability coverage that survives.
  - Confirmed safe: `src/rigel/sim/reads.py` builds its own
    `_region_boundaries` internally and does **not** depend on
    `index.region_df`.

### Phase 7 (deferred / separate effort) — Mappability-corrected effective length

Not part of this overhaul. Tracked as a follow-up: incorporate mappability
into the per-transcript effective length used by `_scoring_impl`'s
log-likelihood. Independent of calibration.

---

## 4 — Risk register

| Risk | Mitigation |
|---|---|
| Pass 0 misclassifies overlapping-gene fragments | Exhaustive `test_categorize.py`; reuse C++-side mixed-strand detection |
| Overhang tolerance too tight or too loose | Unit-tested at multiple values; default 5 bp; config knob |
| Pool starvation on highly-spliced libraries (most reads SPLICED, few EXON_INCOMPATIBLE/INTRONIC/INTERGENIC) | EB shrinkage to global FL handles small pools naturally; `pi_pool ≈ 0` and `gdna_fl_model_quality = "weak"` flagged; downstream EM falls back to splice + strand signals; Phase 5 measures incidence |
| Unstranded `gDNA_FL` fit unstable on real data | `gdna_fl_model_quality` flag + global-FL fallback; Phase 5 validates |
| nRNA contamination in pool biases `gDNA_FL` toward `RNA_FL` | Accepted trade-off; documented; the 1-D mixture still recovers `gDNA_FL` when it differs from `RNA_FL` |
| Multimap filter removes too much pool | Phase 0 audit quantifies multimap fraction; threshold tunable (presently mandatory) |
| Strand-channel gap on stranded libs (no antisense, no cluster) | Phase 5 decision gate; v2 reconsiders only with ground-truth evidence |
| Golden test churn hides real regressions | Eyeball every golden diff in Phase 3 before regenerating |
| Region-partition removal breaks something unexpected | Phase 4 done as a separate commit; easy revert |
| C++ smoke recompile reveals a hidden region dependency | Recompile and run full suite at end of Phase 4; revert if broken |

---

## 5 — Out of scope (deliberate)

- Overlap-cluster channel (deferred to v2 decision gate).
- Cross-region density signals.
- C++ per-locus EM math changes (unless Phase 0 forces warm-start removal).
- Mappability-corrected effective length (Phase 7 follow-up).
- Strand-model / RNA-FL training changes.
- Backward-compat shims for v5 fields. No zombies.
- Protocol-detection (capture vs non-capture) heuristic — SRD is immune
  by construction.

---

## 6 — Phase status snapshot

| Phase | Status | Notes |
|---|---|---|
| 0 | ✅ Complete | `srd_v1_audit.md` |
| 1 | Skipped | User authorized free deletion |
| 2 | ✅ Complete (24 tests green) | Modules: `_categorize.py`, `_fl_mixture.py`, `_fl_empirical_bayes.py`, `_result_srd.py`, `_simple.py`. Pool = EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC; library-agnostic. |
| 3.5 | Next | `partition.py` → `locus_partition.py` rename |
| 3 | ✅ Done | Wired SRD into pipeline; `_result.py` replaced; v5 `compute_locus_priors` body deleted; SRD `compute_locus_priors_from_partitions` added (§2.6); 10 new tests in `test_calibration_simple.py` (1 conditional skip); v5 tests deleted; goldens regenerated. 957 pass, 1 skip. |
| 4a | ✅ Done | v5 calibration modules deleted (`_calibrate.py`, `_em.py`, `_fl_model.py`, `_stats.py`, `_annotate.py`); v5 debug + benchmark scripts deleted; `scripts/calibration/stress_mini/` removed; `profiler.py` updated to SRD signature; v5 docs archived to `docs/deprecated/calibration_v5/`. |
| 4b | Deferred to Phase 6 | Index/region machinery removal carried over: `index.build_region_table` + `region_df`/`region_cr`; `tests/test_regions.py` + `test_mappability.py` region portions; C++ `bam_scanner.cpp` `fl_table`/`region_id` accumulation; `mappability.uniform_region_exposures`; `pipeline.py` region-index build call. Skipped now to avoid index-format/C++ churn before SRD v1 benchmark validation. |
| 5 | In progress | Armis2 benchmark validation; results doc |
| 6 | Pending | CLAUDE.md / README / CHANGELOG **+ Phase 4b cleanup** |
