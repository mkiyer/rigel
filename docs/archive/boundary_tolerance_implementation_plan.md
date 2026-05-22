# Boundary Overlap Tolerance (K) — Implementation Plan

**Date:** 2026-05-08
**Status:** Reviewed and implementation-ready
**Scope:** Add a CLI- and config-configurable boundary-overlap tolerance
`K` (in base pairs, default 3) that gates **both** the scanner-side
boundary-flux numerator and the per-region mask classification, with a
matched correction in the `B_cross` denominator. Here `K` means the
minimum bp clearance required on each side of a boundary; `K=3` accepts
a fragment only if at least 3 bp lie on each side. Wire end-to-end so
the calibration prior is K-consistent at every step.

---

## 1. Motivation and design rule

A 1-bp overhang of an alignment past an exon boundary is overwhelmingly
alignment noise (soft-clip drift, indel slippage, end-position rounding
near GT-AG splice motifs), not a true fragment straddling an exon-intron
junction. The scanner currently uses K=0 — every 1-bp overhang counts as
a boundary-crossing event and contaminates `EXON-INTRON` density
estimation (and the mask state used to stratify the FL histogram and
per-region counts).

A small tolerance (K = 2–3 bp) removes the artefact-dominated tail at the
cost of a small, *known and documented* loss of true short-overhang
fragments — exactly the trade you described.

### Invariant — the only thing we cannot get wrong

> **The same K must be applied to the numerator and the denominator.**

If we change one and not the other we introduce a predictable exposure
bias of approximately `B_cross(K) / B_cross(0)` (or worse). For the
recommended `K=3`, that is roughly `1 - 4 / (E[L] - 1)`. Concretely:

- **Numerator (scanner) side** — three places, all gated on K:
  1. `obs_mask` build: a region whose intersection with the fragment is
    `< q(K)` bp does **not** contribute its type bit to `obs_mask`,
    where `q(K)=max(K,1)`. The intersection must be the exact sum of
    aligned exon-block overlap for that region, not the outer genomic
    span of the fragment.
  2. Boundary-flux counters `u_left[rid]` / `u_right[rid]`: require both
    fragment endpoints to clear the boundary by ≥ `q(K)` bp on the
    inside and ≥ `q(K)` bp on the outside.
  3. `n_unannotated_ref` accounting must respect that `obs_mask == 0`
    is now reachable for fragments that overlap regions only below
    tolerance; treat those as `n_below_tolerance` rather than as
    missing annotation. `n_unobserved` remains reserved for qname groups
    that never produce an observation candidate before `observe()`.
- **Denominator (`_exposure.py`) side**:
  `q(K) = max(K, 1)` and
  `B_cross(K) = Σ pmf[ℓ] · max(ℓ - 2q(K) + 1, 0)`.
  The `max(K, 1)` term preserves the current strict-crossing behavior
  at `K=0`: a fragment must still have at least 1 bp on each side of
  the boundary, yielding `B_cross(0)=E[L]-1`.
- **`contained_exposure_clipped`**: unchanged. K only affects the
  fraction of fragments that *cross* a boundary, not the fraction
  contained inside a region.

### Step-1 disambiguation result (read before reading the rest)

Read the region-index code at
[src/rigel/native/calibration/region_index.h](src/rigel/native/calibration/region_index.h)
and [src/rigel/calibration/regions.py](src/rigel/calibration/regions.py):

- Each region carries a **single** type bit (EXON | INTRON | INTERGENIC).
- The fragment-level `EXON-INTRON` (mask `0b011`) state is **emergent** —
  it appears only when a fragment's exon-blocks span an EXON region and
  an adjacent INTRON region simultaneously, OR'd into `obs_mask`.
- There is **no** "explicit boundary region" type whose width would
  depend on K.

Therefore K is purely a runtime calibration parameter. **The index does
not need to be rebuilt for different K values.** This is the key
simplification driving the design.

---

## 2. Default value

`boundary_tolerance: int = 3`

Justification:

- Mapping noise is empirically ≤ 2 bp in the vast majority of aligner
  outputs. K=3 catches a comfortable margin.
- Loss of true boundary-crossing exposure relative to current `K=0`
  behavior is ≈ `2(K-1) / E[L_gdna] = 4/350 ≈ 1.1%` for `K=3` in a
  typical short-fragment library. Documented and accepted. If we later
  redefine K as "ignore up to and including K bp" rather than
  "require at least K bp", this formula must change to
  `2K / E[L_gdna]` and the CLI wording must change with it.
- Setting K=0 reproduces the current behavior bit-for-bit (regression
  guard, see §6).

Validation domain: `K ≥ 0`. We do not cap K above; a user who runs
K = 50 with a 350 bp FL is intentionally suppressing a large share of
boundary signal, and that is their choice.

---

## 3. Touch-point matrix

| # | Layer | File | Change |
|---|---|---|---|
| 1 | Config | [src/rigel/config.py](src/rigel/config.py) | Add `BamScanConfig.boundary_tolerance: int = 3` (frozen, validated). |
| 2 | CLI | [src/rigel/cli.py](src/rigel/cli.py) | New flag `--boundary-tolerance INT` (default 3); plumbed into `BamScanConfig`. |
| 3 | Pipeline | [src/rigel/pipeline.py](src/rigel/pipeline.py) | Pass `K` from `BamScanConfig` into `scanner.set_regions(...)` and into all density / locus-prior helpers. |
| 4 | Native ABI — accumulator | [src/rigel/native/calibration/accumulator.h](src/rigel/native/calibration/accumulator.h) | Add `int32_t boundary_tolerance_` member; constructor `(int64_t n_regions, int32_t boundary_tolerance = 0)`; add `n_below_tolerance` to `CalibrationPayload`. |
| 5 | Native ABI — accumulator | [src/rigel/native/calibration/accumulator.cpp](src/rigel/native/calibration/accumulator.cpp) | In `observe(...)`: gate `obs_mask |= type_mask(b)` on exact per-region aligned-block overlap ≥ `q(K)`; gate `u_left/u_right` on full endpoint clearance ≥ `q(K)`. |
| 6 | Native ABI + bindings — scanner | [src/rigel/native/bam_scanner.cpp](src/rigel/native/bam_scanner.cpp) | Extend `_bam_impl.BamScanner.set_regions(..., boundary_tolerance=0)`; thread into per-worker `cal_acc` construction and payload export. |
| 7 | Python native shim | [src/rigel/native.py](src/rigel/native.py) | Usually no code change; update module docs only if we document the new `BamScanner.set_regions` argument there. |
| 8 | Density denominator | [src/rigel/calibration/_exposure.py](src/rigel/calibration/_exposure.py) | `boundary_crossing_exposure(fl, *, boundary_tolerance: int = 0)` returning `Σ pmf[ℓ] · max(ℓ - 2q(K) + 1, 0)`. |
| 9 | Global density estimator | [src/rigel/calibration/density_global.py](src/rigel/calibration/density_global.py) | Pass `K` through to `boundary_crossing_exposure(...)`. |
| 10 | Calibration orchestrator | [src/rigel/calibration/_orchestrator.py](src/rigel/calibration/_orchestrator.py) | Plumb `K` through the orchestrator entry point so density and prior helpers see the same value the scanner used. |
| 11 | Locus prior helper | [src/rigel/calibration/locus_prior.py](src/rigel/calibration/locus_prior.py) | Pass `K` through `assemble_priors`, `estimate_locus_gdna`, and `expected_gdna_count_global`. |
| 12 | Scan payload + diagnostics | [src/rigel/calibration/scan_payload.py](src/rigel/calibration/scan_payload.py), [src/rigel/calibration/_diagnostics.py](src/rigel/calibration/_diagnostics.py) | Accept/export `n_below_tolerance`; keep `global_counts.sum() == n_observed` and `Diagnostics.total() == n_observed` by treating below-tolerance as an annotation/QC subcounter, not a ninth partition bucket. |
| 13 | Calibration result | [src/rigel/calibration/_result.py](src/rigel/calibration/_result.py) | Persist `boundary_tolerance` in `CalibrationResult.to_summary_dict` and in the per-pipeline `summary.json` so analyses can verify what K was used. |
| 14 | Tests | new + existing | See §6. |
| 15 | Docs | [docs/MANUAL.md](docs/MANUAL.md), [docs/parameters.md](docs/parameters.md), [CHANGELOG.md](CHANGELOG.md) | Document the parameter, the default, the trade-off, and the K=0 reproducibility note. |

### Minimum complete diff (illustrative, not literal)

**`accumulator.h`**

```cpp
class CalibrationAccumulator {
public:
    explicit CalibrationAccumulator(int64_t n_regions,
                                    int32_t boundary_tolerance = 0)
        : boundary_tolerance_(boundary_tolerance) {
        if (boundary_tolerance_ < 0)
            throw std::invalid_argument("boundary_tolerance must be >= 0");
        // ... existing member init ...
    }
    int32_t boundary_tolerance() const { return boundary_tolerance_; }
private:
    int32_t boundary_tolerance_ = 0;
    // ... rest unchanged ...
};
```

**`accumulator.cpp` — the only place K touches the hot path**

```cpp
const int32_t K = boundary_tolerance_;
const int32_t q = std::max<int32_t>(K, 1);

// Build obs_mask only from regions whose exact aligned-block overlap is >= q.
// Implementation detail: replace the current `hits_`-only merge with a
// parallel `hit_overlap_bp_` vector, or with a tiny Hit{rid, bp} vector.
// For duplicate region IDs across exon blocks, sum bp during the sorted merge.
// Do not use [frag_start, frag_end) as a proxy for spliced fragments; that
// overestimates region overlap and defeats the purpose of K.
auto bit_for_hit = [&](int32_t rid, int64_t overlap_bp) -> uint8_t {
    return (overlap_bp >= q) ? regions.type_mask(rid) : uint8_t{0};
};

// Per-region fan-out must use the same qualified-hit list that built obs_mask.
// Raw hits are still tracked so obs_mask==0 can be split into true
// n_unannotated_ref versus n_below_tolerance.

// Boundary flux: require q-clearance on BOTH sides of the boundary.
// For K=0, q=1, so this is equivalent to the existing strict crossing
// predicates `frag_start < boundary && frag_end > boundary` on integer
// half-open coordinates.
if (flux_eligible && (regions.type_mask(rid) & mask::EXON) != 0) {
    const int64_t rs = regions.start(rid);
    const int64_t re = regions.end(rid);
    if (frag_start + q <= rs && frag_end >= rs + q) payload_.u_left[rid]++;
    if (frag_start + q <= re && frag_end >= re + q) payload_.u_right[rid]++;
}
```

> Endpoint convention: the existing K=0 code uses strict `<` / `>`. The
> correct rewrite is **not** `frag_start + K <= b && frag_end >= b + K`,
> because that includes zero-overhang boundary-touching fragments when
> `K=0`. Use `q=max(K,1)` and shifted integer inequalities as above.
> Confirm with the regression test in §6.1.

**`_exposure.py`**

```python
def boundary_crossing_exposure(
    fl: FragmentLengthModel, *, boundary_tolerance: int = 0
) -> float:
    """B_cross(K) = Σ pmf[ℓ] · max(ℓ - 2q(K) + 1, 0)."""
    if boundary_tolerance < 0:
        raise ValueError("boundary_tolerance must be >= 0")
    pmf = fl.pmf
    ell = np.arange(pmf.size, dtype=np.float64)
    q = max(float(boundary_tolerance), 1.0)
    val = float((pmf * np.maximum(ell - 2.0 * q + 1.0, 0.0)).sum())
    return val if val > 0.0 else 0.0
```

At K=0, `q=1` and this is bit-identical to the current implementation:
`max(ℓ - 1, 0)`. At K=3, the accepted start-position count is
`max(ℓ - 5, 0)`, corresponding to at least 3 bp of fragment span on
both sides of the boundary.

---

## 4. Corner cases and explicit non-fixes

1. **Spliced multi-block fragments.** Do not use the union span
  `[frag_start, frag_end)` for K-gated mask construction. It can
  overestimate region overlap for spliced fragments and falsely keep
  low-overlap type bits. The accumulator should sum exact per-block
  overlap by region while it does the existing sorted merge. Boundary
  flux remains restricted to UNSPLICED fragments, but `obs_mask` feeds
  global counts and FL pools for all fragments, so this precision is
  worth the small hot-path bookkeeping cost.
2. **K beyond region width.** A region of width < `q(K)` never contributes
   its type bit and never receives boundary-flux counts. Correct
   behavior: micro-regions are below tolerance by definition.
3. **K vs `n_unannotated_ref`.** A fragment whose raw region hits are
  all below tolerance becomes `obs_mask = 0`. This must increment
  `n_below_tolerance`, not `n_unannotated_ref`. The latter remains a QC
  signal for BAM/index reference mismatches where no region was hit at
  all. `global_counts[0]`, `fl_hist[0]`, and `n_observed` still include
  below-tolerance fragments so existing balance checks remain valid.
4. **`B_cross(K) = 0`.** When the FL PMF is concentrated below `2q(K)-1`,
   `B_cross` underflows to 0. The downstream callers already handle
   this ("no boundary information available; fall back accordingly").
   Add an explicit summary-stat warning: "boundary_tolerance suppresses
   100% of expected boundary exposure under the gDNA FL — calibration
   will not use boundary evidence."
5. **Index-time vs scan-time semantics.** No index format change. Old
   indices remain compatible. The K parameter lives only in
   `BamScanConfig` and `summary.json`.

### Explicit non-fixes

- **Index-builder K.** The region partition is K-independent (single
  type bit per region), so the index needs no K. Do not add one.
- **Scoring-side overhang.** `scoring.cpp` already uses an "overhang"
  concept (`scoring.cpp:384`); it is **unrelated** — that is a per-fragment
  alignment-quality penalty applied during fragment-vs-transcript
  scoring, not a boundary-event tolerance. Leave it untouched.
- **Per-region K.** Do not add per-region K overrides. K must be
  globally uniform or the numerator/denominator invariant breaks.
- **Asymmetric K.** Do not allow inside-K ≠ outside-K. The same reason.

---

## 5. Migration sequence

The change is independently bisectable in three commits:

### Commit 1 — Native plumbing (no behavior change)

- Add `boundary_tolerance` member + constructor kwarg (default 0) to
  `CalibrationAccumulator`.
- Add `n_below_tolerance` to `CalibrationPayload` and merge/export it,
  but it should remain zero while all callers use K=0.
- Add the same kwarg to `BamScanner::set_regions(...)`; do not put it on
  the `BamScanner` constructor because K only matters when calibration
  regions are installed. Thread it into all per-worker `cal_acc`
  constructions.
- Update bindings.
- **No call site changes yet.** All defaults stay at K=0. Test suite
  passes unchanged.

Acceptance: full test suite green; numerical quantification outputs are
bit-identical. If any golden pins the exact `summary.json` schema, delay
public summary exposure of the new zero-valued counter until Commit 3.

### Commit 2 — Apply K in numerator and denominator together

- Implement the `bit_for_hit` gate and the K-clearance test in
  `observe(...)`.
- Add `boundary_tolerance` kwarg to `boundary_crossing_exposure(...)`.
- Plumb K through `density_global.py`, `_orchestrator.py`, `locus_prior.py`,
  `_result.py`.
- Default is still K=0. **Bit-identical** behavior for existing test
  goldens at K=0.

Acceptance: full test suite green; existing numerical goldens unchanged
at K=0.

### Commit 3 — Surface as CLI parameter; flip default to K=3

- Add `BamScanConfig.boundary_tolerance: int = 3`.
- Add `--boundary-tolerance` CLI argument (default 3).
- Plumb from CLI → `PipelineConfig` → `BamScanConfig` → scanner.
- **Regenerate goldens** in the same commit. Document the regeneration
  in `CHANGELOG.md`.

Acceptance: full benchmark sweep at K=3; calibration-quality metrics
documented in §6.4.

---

## 6. Tests

### 6.1 Bit-identical regression at K=0

- Run the entire current test suite with `boundary_tolerance=0`. All
  goldens must be unchanged.

### 6.2 Accumulator unit tests (new)

Add `tests/test_boundary_tolerance.py`:

- **mask gating**: a single unspliced fragment overlapping EXON by
  100 bp and INTRON by 1 bp → `obs_mask == 0b001` (EXON only) at K=2,
  `obs_mask == 0b011` (EXON|INTRON) at K=0.
- **boundary flux gating**: 1-bp overhang at K=3 does **not** fire
  `u_left/u_right`; 3-bp overhang does.
- **endpoint convention**: K=0 produces bit-identical counters to the
  pre-change implementation on a small synthetic fixture.
- **below-tolerance accounting**: a fragment with raw region hits but no
  qualified region hits increments `n_below_tolerance` and `global_counts[0]`,
  not `n_unannotated_ref`.
- **micro-region**: a region of width 2 bp under K=3 contributes no
  mask bit and no boundary count, regardless of fragment overlap.
- **B_cross(K) table**: assert `B_cross(K=0) - B_cross(K=3) ≈ 4` for
  FL = N(350, 100) (within numerical tolerance).
- **invariant**: at K=0, `B_cross(fl, boundary_tolerance=0) == E[L] - 1`
  (existing invariant).

### 6.3 Pipeline integration test (new)

- Run two synthetic mini-pipelines that differ only in K (0 vs 3).
- Assert `summary.json["boundary_tolerance"]` matches the input.
- On a pure uniform gDNA fixture with no simulated boundary noise, assert
  the boundary-density estimate is stable between K=0 and K=3 within
  Poisson noise. The numerator and denominator should scale together;
  mismatched movement is the sentinel for a numerator/denominator bug.

### 6.4 Synthetic-sweep acceptance gate (the optimality lock)

Run the gDNA-only uniform-density acceptance fixture from the previous
plan iteration:

1. Synthetic sim: pure uniform gDNA, no RNA, no nRNA, FL = δ(350) (or
   known PMF). Single chromosome with varied locus geometries.
2. Run pipeline at K=3.
3. Assert:
   - $\rho_{ig} = \rho_{in} = \rho_b$ within $3\sqrt{N}/L$ Poisson noise
     on each branch.
   - For each locus, $(\eta_g - \text{realized\_gdna}) / \sqrt{\text{realized\_gdna}} < 3$.
   - $\Sigma\alpha / \Sigma\text{gdna} \in [0.99, 1.01]$.
4. Record this as a permanent regression gate under
   `tests/golden/calibration/uniform_gdna_K3/`.

### 6.5 Real-data benchmark sweep

Re-run the full synthetic sim suite (`/Users/mkiyer/Downloads/rigel_runs/sim_synthetic`)
at K=3 and compare to the existing K=0 results:

| Metric | K=0 (current) | K=3 (target) |
|---|---|---|
| ρ_ex/ρ_ig | 0.92 | should approach 1.0 if K artefact dominated; document residual |
| Σα/Σgdna | 0.957 | depends on outcome of ρ ratio |
| Pearson(α, gdna) | 0.987–0.996 | should not regress |
| RNA→gDNA misclass | 0.13–3.03% | should not regress materially |
| Spurious nRNA | 0–0.36% | unchanged (Issue #2 is K-independent) |

If `ρ_ex/ρ_ig` does not move toward 1.0 at K=3, then the K hypothesis
is partial, and the residual must be explained by region geometry or
FL asymmetry (next investigation step). Either outcome is publishable
diagnostic data.

---

## 7. Configuration and surface

### CLI

```bash
rigel quant --bam X.bam --index idx/ -o out/ \
    --boundary-tolerance 3   # default
```

Help text (paste into `cli.py`):

> Minimum number of base pairs a fragment must clear an exon-intron
> boundary on each side to be counted as a boundary-crossing event
> (default 3 bp). Removes alignment-noise artefacts at the cost of a
> small loss of true short-overhang fragments. Set to 0 to reproduce
> pre-2026.05 behavior.

### Config

```python
@dataclass(frozen=True)
class BamScanConfig:
    ...
    boundary_tolerance: int = 3
    """Min bp a fragment must clear an exon/intron boundary on each
    side to be counted as a boundary-crossing event (and to contribute
    its EXON|INTRON bit to the per-fragment region mask). Applied
    consistently in the scanner and in the matched ``B_cross`` exposure
    denominator."""
```

### Persisted in summary

```json
{
  "calibration": {
    "boundary_tolerance": 3,
    ...
  }
}
```

---

## 8. Documentation deltas

- [docs/MANUAL.md](docs/MANUAL.md): new "Boundary tolerance" subsection
  under the calibration section. Explain the trade-off and the K=0
  reproducibility note.
- [docs/parameters.md](docs/parameters.md): add the parameter to the
  reference table.
- [CHANGELOG.md](CHANGELOG.md): document the new flag, the default, and
  the golden regeneration in commit 3.
- [docs/METHODS.md](docs/METHODS.md): add the `B_cross(K)` formula and
  the numerator-side K-clearance rule. State the invariant explicitly.

---

## 9. Resolved accounting decision

Do **not** let K change the meaning of `n_unannotated_ref`. That counter
continues to mean "observed fragment, but no raw region hit at all,"
which preserves its value as a BAM/index reference mismatch QC signal.

Add `n_below_tolerance` for the new K-specific case: raw region hits
exist, but no hit has enough aligned-block overlap to contribute a mask
bit. These fragments still increment `global_counts[0]`, `fl_hist[0]`,
and `n_observed`, so the existing balance invariants remain exact.
Because they are already included in `global_counts[0]`, do not add
`n_below_tolerance` into `Diagnostics.total()` as an extra bucket; expose
it as an auxiliary QC field in the calibration summary.

---

## 10. Acceptance order

1. **Commit 1** lands and is bisectable (no behavior change).
2. **Commit 2** lands; full suite + goldens unchanged at K=0.
3. **Commit 3** lands; goldens regenerated at K=3; benchmark sweep
   recorded; uniform-gDNA fixture from §6.4 added as a permanent gate.
4. **Optimality declaration:** if §6.4 passes, the Bayesian prior is
   mathematically locked under the global-only model. Phase 4
   (independent-flank prior) becomes the next focus, not a fix for
   Phase 2 issues.

---

## 11. Out of scope (deferred)

- **gDNA → synthetic-nRNA leakage** (Issue #2 from the previous report).
  K does not affect read-level component competition. Sentinel-only.
- **Sparse-locus over-prior** (Issue #3). Architecturally addressed by
  Phase 4 independent-flank.
- **LOO mega-locus guard** (v3 plan §2.1). Independent work; landing K
  does not block it.
