# Accumulator + Calibrator — Consolidated Implementation Plan

**Date:** 2026-05-29
**Status:** Ready for implementation. Supersedes
[`../accumulator/01_implementation.md`](../accumulator/01_implementation.md)
and [`../caljointmodel/05_implementation_plan.md`](../caljointmodel/05_implementation_plan.md).

This document merges the fractional-accumulator rewrite and the
calibration-v6 (jointmodel) rebuild into one phased plan. The two
efforts share a substrate (`AccumulatorPayload`) and were drafted
separately; running them sequentially under a single plan removes
double-work, dropped overlaps, and ambiguous handoffs.

## 0. Document scope and cross-references

This file replaces the *phase plans* of the two predecessor docs. It
does **not** restate the design specs. Authoritative specs live in:

| Topic | Authoritative doc |
|---|---|
| Fractional accumulator data structures, deposit algorithm, mass conservation invariants | [`../accumulator/00_design.md`](../accumulator/00_design.md) |
| Legacy accumulator audit (no bugs; locked decisions) | [`../accumulator/audit_phase1.md`](../accumulator/audit_phase1.md) |
| Calibration goals (G2/G3/G4), substrate, principles | [`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) |
| Generative model (NB count, BB strand, deterministic spliced RNA, Gamma exposure) | [`../caljointmodel/01_generative_model.md`](../caljointmodel/01_generative_model.md) |
| Why the legacy calibration must be burned | [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md) |
| Inference algorithm (E-step, M-step, closed forms, Newton) | [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) |
| Public API (`CalibrationConfig`, `CalibrationResult`, locus prior pseudocount formulas) | [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) |
| Validation plan (unit, synthetic, scenario, armis2) | [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) |

Any conflict between this plan and a spec doc is resolved by
updating the spec doc, not by deviating in implementation.

## 1. State at start of plan

- **Phase 1 (audit + spec) of the accumulator: DONE.** Memo at
  [`../accumulator/audit_phase1.md`](../accumulator/audit_phase1.md);
  Python reference at `tests/native/_accumulator_reference.py`; 10
  spec tests under `tests/native/test_accumulator_spec.py`.
- **Phase 2 (native rewrite) of the accumulator: DONE.** New
  `src/rigel/native/calibration/accumulator.{h,cpp}` per spec §3, §7.1.
- **Phase 3 (nanobind binding + Python wrapper) of the accumulator:
  DONE.** Exposed as `rigel.native.Accumulator` with zero-copy
  ndarray views; thin Python wrapper at `src/rigel/_accumulator.py`
  provides dataclass-style indexed access (`acc.regions[i].contained[ch]`,
  `acc.boundaries[i].mass_left[ch]`, etc.).
- **20/20 spec tests pass** (`pytest tests/native/`).
- **Calibration v5 partially severed.** `bam_scanner.cpp` no longer
  references `CalibrationAccumulator`/`CalibrationPayload`. The
  scanner still emits `result["calibration"] = None`. Python
  calibration-v5 modules (`_arrays.py`, `_diagnostics.py`,
  `_orchestrator.py`, `scan_payload.py`, `fractional_evidence.py`,
  `region_count_ledger.py`, `signature.py`, `_simple.py`,
  `_categorize.py`, etc.) and their tests are still in the tree but
  expected to be import-broken — to be removed in Phase A below.
- **`docs/accumulator/01_implementation.md` and
  `docs/caljointmodel/05_implementation_plan.md` are superseded by
  this document** (they may stay in the tree for history).

## 2. Cross-cutting conventions

- **Build between every phase.** `pip install --no-build-isolation -e .`
  must succeed before any test command runs. Phases B and D
  recompile the native module; Phases A, C, E–G are Python-only.
- **One PR per phase.** Each phase ends at a checkpoint with its own
  acceptance gate. The next phase begins only after the previous
  gate is green.
- **No premature deletions.** Do not delete a file outside the
  current phase's file-touch matrix.
- **No heredoc in terminal.** Per workspace convention. Multi-line
  diagnostics go in `scripts/debug/`.
- **Calibration is goal-directed.** Every line of new calibration
  code must trace to G2, G3, or G4
  ([`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) §1).
- **Sufficient statistics, not fragments.** Calibration is $O(R)$
  scalar algebra over `AccumulatorPayload` views. Per-fragment
  iteration in calibration is a bug.
- **FL is not a calibration channel.** FL stays in the downstream
  EM scorer
  ([`scoring.py`](../../src/rigel/scoring.py),
  [`frag_length_model.py`](../../src/rigel/frag_length_model.py))
  and is consumed per-fragment. See
  [`../accumulator/audit_phase1.md`](../accumulator/audit_phase1.md)
  decision #6 and
  [`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) §1.

---

## 3. Phase overview

```
Phase A — Burn legacy calibration       (Python-only deletion, ~1 PR)
Phase B — Wire accumulator into scanner (C++ + Python substrate, ~1 PR)
Phase C — Calibrator scaffold           (Public types, validation, ~1 PR)
Phase D — Calibrator implementation     (E-step, M-step, outer loop, 3 sub-PRs)
Phase E — Integrate                     (Locus prior consumer rewrite)
Phase F — Validate                      (Synthetic, scenario, armis2)
Phase G — Cleanup                       (Goldens, ruff, docs, postmortem)
```

| Phase | Touches | Build needed | Output |
|---|---|---|---|
| A | Python only | no | Clean slate; `rigel quant` runs through a stub |
| B | C++ + Python | **yes** | Native scanner emits `AccumulatorPayload` |
| C | Python only | no | Public `CalibrationConfig` / `CalibrationResult`; stub `calibrate()` |
| D | Python only | no | Working calibrator; aborts at locus prior consumer |
| E | Python only | no | End-to-end `rigel quant` runs on a scenario BAM |
| F | Python + benchmarks | no | Synthetic + scenario + armis2 validated |
| G | Mechanical | no | Goldens regenerated; docs updated |

---

## 4. Phase A — Burn legacy calibration

**Goal.** Excise every line of calibration-v5 from the live tree so
Phases B–E build cleanly on top of a known-empty slate. The
fractional accumulator (live in `rigel.native.Accumulator`) is *not*
yet wired into the scanner — that is Phase B. During Phase A the
scanner continues to emit `result["calibration"] = None`.

### A.1 Archive salvageable algorithms

Move to `archive/calibration_legacy_2026_05/` (preserving relative
paths under `src/rigel/calibration/...`). These algorithms may be
useful as references for downstream EM work later:

| File | Why archive |
|---|---|
| `src/rigel/calibration/_fl_empirical_bayes.py` | EB Dirichlet smoothing pattern |
| `src/rigel/calibration/_fl_mixture.py` | Closed-form FL mixture EM |
| `src/rigel/calibration/fl.py` | Public FL surface |
| `src/rigel/calibration/_simple.py` | SRD gDNA FL recovery |

Use `git mv` (not plain `mv`) so history is preserved. Write
`archive/calibration_legacy_2026_05/README.md` with the SHA, the
ARCHIVE-vs-DELETE distinction, and a no-resurrection warning citing
this plan.

### A.2 Delete everything else under `src/rigel/calibration/`

Per [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md)
§4. Concrete delete list (final):

```
src/rigel/calibration/_arrays.py
src/rigel/calibration/_categorize.py
src/rigel/calibration/_diagnostics.py
src/rigel/calibration/_orchestrator.py
src/rigel/calibration/_result.py
src/rigel/calibration/accumulator_view.py
src/rigel/calibration/adaptive_prior.py
src/rigel/calibration/background_model.py
src/rigel/calibration/bootstrap.py
src/rigel/calibration/boundaries.py
src/rigel/calibration/boundary_model.py
src/rigel/calibration/boundary_sweep.py
src/rigel/calibration/compartment_strand_deconv.py
src/rigel/calibration/coverage_weight.py
src/rigel/calibration/density_global.py
src/rigel/calibration/density_model.py
src/rigel/calibration/density_observation.py
src/rigel/calibration/exposure.py
src/rigel/calibration/fractional_evidence.py
src/rigel/calibration/fusion.py
src/rigel/calibration/latent_states.py
src/rigel/calibration/locus_partition.py    # only if not consumed by locus.py — verify
src/rigel/calibration/prior.py
src/rigel/calibration/region_count_ledger.py
src/rigel/calibration/region_partition.py   # only if not consumed by locus.py — verify
src/rigel/calibration/regions.py
src/rigel/calibration/scan_payload.py
src/rigel/calibration/signature.py
src/rigel/calibration/strand_deconv.py
```

Also delete `src/rigel/native/calibration/region_signature.h` and
`src/rigel/native/calibration/region_index.h` if no remaining caller
needs them. The new accumulator does not. Verify via grep before
deleting.

### A.3 Delete dependent tests

```
tests/test_adaptive_prior.py
tests/test_background_model.py
tests/test_bayesian_prior_acceptance.py
tests/test_boundaries.py
tests/test_boundary_model.py
tests/test_boundary_sweep.py
tests/test_calibrate.py
tests/test_calibration_accumulator.py
tests/test_calibration_integration.py
tests/test_calibration_iteration.py
tests/test_calibration_prior.py
tests/test_calibration_result.py
tests/test_compartment_strand_deconv.py
tests/test_coverage_weight.py
tests/test_density_global.py
tests/test_density_model.py
tests/test_density_observation.py
tests/test_exposure.py
tests/test_fl_eff_len_cache.py
tests/test_fl_models.py
tests/test_latent_states.py
tests/test_locus_partition.py
tests/test_region_count_ledger.py
tests/test_region_index_native.py
tests/test_region_partition.py
tests/test_region_persist.py
tests/test_region_unspliced_mass.py
```

Plus any `tests/test_strand_*` that reference deleted modules
(case-by-case grep).

**Keep.** `tests/native/test_accumulator_spec.py` (the new
substrate's spec), `tests/test_pipeline_smoke.py` (Phase E will fix
it up), all non-calibration tests (`test_em_impl.py`, etc.).

### A.4 Add the stub module

Create `src/rigel/calibration/_stub.py`:

```python
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

@dataclass(frozen=True, slots=True)
class CalibrationConfig:
    max_outer_iterations: int = 25
    mass_rel_tol: float = 1e-4
    phi_floor: float = 1e-6
    boundary_split_factor: float = 0.5

@dataclass(frozen=True, slots=True)
class CalibrationResult:
    # All zero / empty placeholders matching the Phase C schema.
    n_regions: int = 0
    converged: bool = False
    n_iterations: int = 0
    # ... (all fields from doc 04 §5, defaulted to zeros / np.array([]))

def calibrate(*_args, **_kwargs) -> CalibrationResult:
    raise NotImplementedError(
        "Calibration burn-down; see docs/acc_caljointmodel/ and "
        "docs/caljointmodel/. Implementation lands in Phase C/D."
    )
```

Rewrite `src/rigel/calibration/__init__.py` to re-export exactly
those three symbols and nothing else.

### A.5 Excise call sites

Files to touch:

- **`src/rigel/pipeline.py`** — strip the entire legacy calibration
  block (~150 lines) from `quant_from_buffer`. Replace with a
  single call to `calibrate(...)` that will raise until Phase D.
  (The pipeline does not reach the calibration call in current main
  because the scanner emits `result["calibration"] = None`; this
  edit makes the intent explicit.)
- **`src/rigel/estimator.py`** — drop calibration imports. The EM
  no longer reads anything from a `CalibrationResult` directly; it
  receives priors from the locus-prior consumer (whose Phase A stub
  raises `NotImplementedError`).
- **`src/rigel/locus.py::compute_locus_priors_from_partitions`** —
  replace body with `raise NotImplementedError("Locus prior consumer
  lands in Phase E.")`. Preserve signature.
- **`src/rigel/config.py`** — delete all legacy calibration
  sub-configs. Replace with the single stub `calibration:
  CalibrationConfig` field on `PipelineConfig`
  ([`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §8).
- **`src/rigel/cli.py`** — strip calibration-related CLI flags
  (everything that referenced the v5 sub-configs).
- **`src/rigel/native/bam_scanner.cpp`** — remove the now-dead
  `obs_set` / `obs_splice` / `obs_ref` / `obs_fs` / `obs_fe` /
  `obs_fragment_strand` / `obs_exons` capture locals from the
  fragment-finalize block (they were only used by the deleted
  `cal_acc.observe(...)` call). Keep `set_regions(...)` and
  `RegionIndex` for now — Phase B replaces them with the
  accumulator's `region_edges` API.

### A.6 CHANGELOG

```
2026-05-29: Calibration v5 burn-down. Replacement is the joint
fractional-accumulator + calibration-v6 rewrite per
docs/acc_caljointmodel/.
```

### A.7 Acceptance gate

- [ ] `python -c "import rigel"` succeeds.
- [ ] `pytest --collect-only tests/` exits clean (no import errors).
- [ ] `pytest tests/native/` → 20 passed (accumulator spec still green).
- [ ] `git grep -nE 'background_model|fusion|latent_states|boundary_sweep|strand_deconv|adaptive_prior|p_unexpressed|fit_status|fused_soft_label|fractional_evidence|RegionCountLedger' src/` returns zero hits.
- [ ] `wc -l src/rigel/calibration/*.py` is ≤ 200 lines (stub + `__init__` only).
- [ ] `rigel quant ...` on a scenario BAM scans successfully and aborts at the calibration stub with `NotImplementedError`.

### A.8 Rollback

Revert the Phase A PR.

---

## 5. Phase B — Wire accumulator into the scanner

**Goal.** Make `BamScanner` build per-reference `Accumulator`
instances during the BAM scan and emit the resulting
`AccumulatorPayload` as `result["calibration"]`. After Phase B the
substrate is real (no longer `None`) and Phase C can build on it.

This phase combines the original accumulator-plan Phase 4 (cutover)
with the original caljointmodel Phase 3.1 (native substrate). The
substrate emitted matches the schema in
[`../accumulator/00_design.md`](../accumulator/00_design.md) §7.3,
not the older "extension" framing in
[`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §2.

### B.1 Native: per-reference accumulator collection

In `src/rigel/native/calibration/accumulator.{h,cpp}` add (or extend
an existing helper) a small wrapper that owns one `Accumulator` per
reference and supports thread-safe merging:

```cpp
class AccumulatorSet {
 public:
  // Construct one Accumulator per reference from a flat edge layout:
  //   edges      : int64[E]   sorted edges across all refs, contiguous per ref
  //   ref_offsets: int64[F+1] first edge of each ref (ref f uses
  //                            edges[ref_offsets[f] .. ref_offsets[f+1]))
  AccumulatorSet(const int64_t* edges, size_t n_edges,
                 const int64_t* ref_offsets, size_t n_refs);

  Accumulator& at(int32_t ref_id);
  const Accumulator& at(int32_t ref_id) const;
  size_t n_refs() const noexcept;

  // Element-wise sum of another AccumulatorSet into this one
  // (per-reference, per-Region / per-Boundary).
  // Used to merge per-worker AccumulatorSets after the scan.
  void merge_from(const AccumulatorSet& other);
};
```

`merge_from` requires identical edges (asserts at start). It does
element-wise `+=` on the `Region.contained[c]` / `Boundary.mass_*[c]`
/ `Boundary.flux[c]` channels.

### B.2 Native: scanner integration

In `src/rigel/native/bam_scanner.cpp`:

1. Replace `BamScanner::set_regions(...)`'s 9-array signature with
   a 3-array `set_region_edges(edges, ref_offsets, n_refs)`. Keep
   the old `RegionIndex` if some other path still needs it; if not,
   delete it (and `region_index.h`, `region_signature.h`) here.
2. `WorkerState` holds `std::unique_ptr<AccumulatorSet> acc_set` (one
   per worker, sized from the partition at scan start).
3. At the fragment observation site (currently the dead `obs_set`
   block — removed in Phase A.5), call:

   ```cpp
   if (acc_set && !is_multimap &&
       any_non_chimeric_resolved && obs_blocks_valid) {
       Accumulator& acc = acc_set->at(obs_ref);
       acc.deposit(block_starts, block_ends, n_blocks,
                   spliced /* per-fragment */,
                   strand_pos /* per-fragment */);
   }
   ```

   `block_starts`/`block_ends` come from the resolved fragment's
   `exons` member; `spliced` is the per-fragment splice flag
   (per audit memo §5); `strand_pos` from `obs_fragment_strand`
   (`STRAND_POS` → true, `STRAND_NEG` → false; ambiguous → return
   early per audit memo §6).
4. After scan, merge per-worker `acc_set` into a `merged_acc_set` on
   the `BamScanner`.
5. In `BamScanner::build_result()`, populate
   `result["calibration"]` with the schema from B.3 (zero-copy
   ndarray views over the merged accumulator's buffers, owned by
   capsules tied to a Python-side owner that keeps the merged
   accumulator alive).

### B.3 Native: nanobind binding for the payload

Expose an `AccumulatorPayload` view via the existing `_bam_impl`
module. Per-reference layout:

```python
result["calibration"] = {
    # Topology
    "ref_region_offsets":   np.ndarray,  # int64[F+1]
    "ref_boundary_offsets": np.ndarray,  # int64[F+1]
    "region_edges":         np.ndarray,  # int64[sum(n_edges_per_ref)]

    # Region-keyed: shape (N_total, 4)  uint32
    # Channels: 0=unspl_pos, 1=unspl_neg, 2=spl_pos, 3=spl_neg
    "region_contained":     np.ndarray,

    # Boundary-keyed: shape (B_total, 4)  channels same as above
    "boundary_mass_left":   np.ndarray,  # float32
    "boundary_mass_right":  np.ndarray,  # float32
    "boundary_flux":        np.ndarray,  # uint32
}
```

This is the **only** schema the Python side must consume. The
column-per-channel `int8` style of the old `RegionCountLedger` is
gone.

### B.4 Python adapter: `AccumulatorPayload`

Create `src/rigel/scan_payload.py` (resurrected fresh; the legacy
version was deleted in Phase A):

```python
@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    """View over the native accumulator's per-reference buffers.

    Channel encoding (axis -1, size 4):
      0 = unspliced, sense
      1 = unspliced, antisense
      2 = spliced,   sense
      3 = spliced,   antisense
    """
    ref_region_offsets:   np.ndarray  # int64[F+1]
    ref_boundary_offsets: np.ndarray  # int64[F+1]
    region_edges:         np.ndarray  # int64[E]

    region_contained:     np.ndarray  # uint32[N_total, 4]
    boundary_mass_left:   np.ndarray  # float32[B_total, 4]
    boundary_mass_right:  np.ndarray  # float32[B_total, 4]
    boundary_flux:        np.ndarray  # uint32[B_total, 4]

    @classmethod
    def from_scan_result(cls, scan_result: dict) -> "AccumulatorPayload":
        cal = scan_result["calibration"]
        if cal is None:
            raise ValueError(
                "scan_result['calibration'] is None; "
                "BamScanner.set_region_edges(...) was not called."
            )
        return cls(**{k: cal[k] for k in cls.__annotations__})
```

### B.5 Pipeline integration

In `src/rigel/pipeline.py`:

1. Add `build_region_partition_from_index(index)` helper returning
   `(edges, ref_offsets)` arrays from the existing reference index.
   The partition design lives elsewhere (mirrors today's
   `RegionArrays`/`region_partition.py` logic — those modules were
   deleted in Phase A, so the partition-builder needs to live in
   `index.py` going forward). Use whatever per-reference region
   layout the calibrator needs (exonic vs intronic vs intergenic
   regions per the existing partition scheme).
2. Before scan: call `scanner.set_region_edges(edges, ref_offsets, n_refs)`.
3. After scan: `payload = AccumulatorPayload.from_scan_result(result)`.
4. Pass `payload` to the calibration stub (which still raises in
   Phase B; Phase D fills it in).

### B.6 Tests

- `tests/test_scanner_accumulator_integration.py` — build a tiny
  synthetic name-sorted BAM (or reuse an existing scenario BAM); run
  the full scan with a hand-crafted region partition; assert the
  resulting `AccumulatorPayload` matches a reference computed by
  running `tests/native/_accumulator_reference.Accumulator.deposit(...)`
  on the same fragment stream. ~10 fragments, 2 regions, easy to
  verify by hand.
- `tests/test_accumulator_payload.py` — payload-level invariants:
  shape, dtype, contiguity, length consistency
  (`region_contained.shape[0] == ref_region_offsets[-1]` etc.); GC
  test (numpy views keep underlying C++ buffers alive).
- `tests/native/test_accumulator_spec.py` — still 20 passing.

### B.7 Acceptance gate

- [ ] Native build succeeds (`pip install --no-build-isolation -e .`).
- [ ] `pytest tests/native/` → 20 passed.
- [ ] `pytest tests/test_scanner_accumulator_integration.py -v` passes.
- [ ] `pytest tests/test_accumulator_payload.py -v` passes.
- [ ] `pytest --collect-only tests/` clean.
- [ ] `rigel quant ...` reaches the calibration stub with a
      non-None payload; aborts on `NotImplementedError`.
- [ ] No remaining grep hits for legacy region/signature symbols in
      `src/`.

### B.8 Rollback

Revert the Phase B PR. The accumulator spec tests still pass (the
new module is standalone).

---

## 6. Phase C — Calibrator scaffold

**Goal.** Land the real `CalibrationConfig` and `CalibrationResult`
types, the substrate-validation layer, and a placeholder
`calibrate(...)` that returns zero-mass / unit-exposure. No
inference yet. After Phase C the pipeline runs through calibration
and aborts at the Phase E locus-prior consumer.

### C.1 Module surface

Per [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §1:

```
src/rigel/calibration/calibrate.py       # new (single file)
src/rigel/calibration/__init__.py        # re-export CalibrationConfig, CalibrationResult, calibrate
```

Delete `src/rigel/calibration/_stub.py` (Phase A's placeholder).

### C.2 Calibration substrate adapter

The calibrator consumes a **per-substrate-set view** of the
`AccumulatorPayload`. Per
[`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §2,
the calibrator runs three times — once on contained regions, once on
left boundaries, once on right boundaries — sharing the global
hyperparameters.

Add to `src/rigel/scan_payload.py`:

```python
@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Per-substrate-set view (contained, left, or right).

    All arrays length |S|. Channel reductions of the
    AccumulatorPayload follow:
        n_unspl[s] = region_contained[s, 0] + region_contained[s, 1]
                   (or boundary_flux[s, 0] + boundary_flux[s, 1] for boundary sets)
        n_spliced[s] = region_contained[s, 2] + region_contained[s, 3]
        k_plus[s]    = region_contained[s, 0]   (or boundary_flux[s, 0])
    """
    n_unspl:   np.ndarray  # int64[|S|]
    n_spliced: np.ndarray
    k_plus:    np.ndarray
    L_eff:     np.ndarray  # float64
    kappa_rna: np.ndarray  # float64; from StrandModel.p_r1_sense
    # ... plus topology fields to project boundary masses back to regions
```

Note the boundary substrates expose **`flux`** as the count signal
(not mass) — per accumulator design §6 the boundary `mass_*` fields
are per-block-side weights, not per-fragment counts. The integer
`flux` field is the correct per-fragment-event tally for the
calibrator's NB / BB likelihoods. Boundary mass channels feed only
into the final per-region exposure aggregation
(§4 of [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md)).

> **Open: confirm at Phase C kickoff** that `flux` is the right
> count signal for the boundary substrate's NB likelihood. If
> per-block-side mass is preferred for some reason, the calibrator's
> count likelihood needs to switch to a continuous-mass analog.
> Decision goes in [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §2.

### C.3 Public types

Per [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §3, §5:

- `CalibrationConfig` — 4 fields, all numeric, no decision thresholds.
- `CalibrationResult` — per-region mass arrays (G2), boundary mass
  arrays (G3), exposure posterior (G4), 5 library hyperparameters,
  convergence diagnostics, provenance.
- `__post_init__` invariants per §5.1.

### C.4 Exceptions

```python
class CalibrationSubstrateError(ValueError): ...
class CalibrationConvergenceError(RuntimeError): ...
```

### C.5 Placeholder `calibrate()`

```python
def calibrate(substrate, strand_model, config=CalibrationConfig()):
    _validate_substrate(substrate, strand_model)  # raises CalibrationSubstrateError
    R = substrate.n_regions
    return CalibrationResult(
        mass_g_contained=np.zeros(R), mass_d_contained=np.zeros(R),
        mass_g_left=np.zeros(R),  mass_d_left=np.zeros(R),
        mass_g_right=np.zeros(R), mass_d_right=np.zeros(R),
        omega=np.ones(R),
        log_omega_var=np.full(R, config.phi_floor),
        rho_0=0.001, phi=0.1,
        rho_d_bb=0.01, rho_r_bb=0.01, eps_s=1e-3,
        n_iterations=0, converged=False,
        mass_change_history=np.zeros(0),
        n_regions=R, config=config,
    )
```

### C.6 Tests

`tests/calibration/`:
- `test_config_defaults.py` — defaults match doc 04 §3.
- `test_result_invariants.py` — each `__post_init__` check fires.
- `test_substrate_invariants.py` — each invariant violation in
  doc 04 §4.1 raises `CalibrationSubstrateError`.

### C.7 Acceptance gate

- [ ] `pytest tests/calibration/ -v` passes.
- [ ] `pytest tests/native/` still 20 passed.
- [ ] `rigel quant ...` runs through calibration returning the
      placeholder; aborts at the locus prior `NotImplementedError`
      from Phase A.5.

### C.8 Rollback

Revert Phase C PR. Phase A's `_stub.py` returns.

---

## 7. Phase D — Calibrator implementation

**Goal.** Replace the placeholder `calibrate()` with the real
inference per [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md).

Three sub-phases, each independently committable.

### D1 — Per-region E-step (G2/G3 unified deconvolution)

Implements [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §3.

- `_per_region_estep(substrate_set, pi_g_prior, rho_d_bb, rho_r_bb,
  kappa_rna, eps_s, omega, rho_0, L_eff, phi, M_d_unspl_prev) ->
  (M_g, M_d, k_plus_g_hat, k_plus_d_hat)`. Vectorized over `|S|`.
  `kappa_d = 0.5` is a module-level constant.
- `_boundary_half_split(M_g_L, M_g_R, ref_offsets) ->
  M_g_boundary_contribution`. Symmetric split per
  [`../caljointmodel/01_generative_model.md`](../caljointmodel/01_generative_model.md) §7.

Tests in `tests/calibration/test_g2_g3_deconvolution.py`:
- Recover known $\pi_r^{(g)}$ within 5% for regions with ≥ 20 fragments.
- **Hybrid-capture sanity** (captured exon + depleted flank).
- **Paralog sanity** (126/26 strand split → $M_r^{(d, \text{cont})} < 10$;
  risk-flagged; see audit memo §2.5).
- Mass conservation: $M_r^{(g)} + M_r^{(d)} =$ total count to $10^{-9}$.
- Vectorized output matches scalar reference on 1000 random regions.

### D2 — G4 closed-form exposure posterior

Implements [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §4.

- `_update_exposure(M_g_tot, rho_0, L_eff, phi) -> (omega, log_omega_var)`.

Tests in `tests/calibration/test_g4_exposure.py`:
- Recover known $\omega_r$ within 10% for regions with ≥ 10 expected gDNA fragments.
- Empty region: $\omega = 1$, $\log\sigma^2 = \phi$.
- Variance scales as $1/(1/\phi + M)$.

### D3 — Global M-step + outer loop

Implements [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §5, §1.

- `_m_step_rho_0`, `_m_step_eps_s` — closed forms.
- `_m_step_phi`, `_m_step_rho_d_bb`, `_m_step_rho_r_bb` —
  `scipy.optimize.minimize_scalar` with moment-estimator warm start.
- `_update_pi_g_prior`.
- `calibrate(...)` outer loop per [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §1.

Tests:
- `tests/calibration/test_m_step.py` — each M-step recovers its
  truth within tolerance (per
  [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §1.3).
- `tests/calibration/test_outer_loop.py` — mass-change monotone
  decreasing; convergence in ≤ 25 iterations.
- `tests/calibration/test_e2e_synthetic.py` —
  $R = 1000$ synthetic substrate; recover hyperparameters and per-region
  exposure within tolerances from
  [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §2.

### Phase D acceptance gate

- [ ] All `tests/calibration/test_g2_g3_*`, `test_g4_*`,
      `test_m_step*`, `test_outer_loop*`, `test_e2e_synthetic*` pass.
- [ ] Magic-number audit on `calibrate.py`: ≤ 8 numeric literals
      (per [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) §8).
- [ ] Mass-change diagnostic monotone-decreasing on every synthetic run.
- [ ] `rigel quant ...` runs through calibration successfully;
      still aborts at the locus prior consumer.

### Phase D rollback

Revert per-sub-phase. D1 must land before D3 (D3 calls D1).

---

## 8. Phase E — Integrate

**Goal.** Wire `CalibrationResult` into the locus-prior consumer so
the pipeline runs end-to-end.

### E.1 Rewrite `compute_locus_priors_from_partitions`

Per [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) §6:

```python
def compute_locus_priors_from_partitions(
    locus_partitions, region_arrays, calibration, *, kappa
) -> LocusPriors: ...
```

Implement the §6.2 pseudocount formulas verbatim:

$$
\alpha_t^{(d)} = \kappa \sum_{r \in r(t)} \phi_{t,r} w_r \bigl[M_r^{(d, \text{cont})} + \tfrac{1}{2}(M_r^{(d,L)} + M_r^{(d,R)})\bigr]
$$
$$
\alpha_t^{(g)} = \kappa \sum_{r \in r(t)} \phi_{t,r} w_r \hat{\omega}_r \rho_0 L_{t,r}^{\text{eff}}
$$

Delete the Phase A `NotImplementedError` stub.

### E.2 Tests

`tests/calibration/test_locus_prior_consumer.py`:
- `test_pseudocount_formula_exact` — matches §6.2 formulas to $10^{-10}$.
- `test_mass_conservation_invariant` — §6.3 invariant holds.
- `test_symmetric_paralog_locus_symmetric_pseudocount`.
- `test_empty_region_only_gdna_prior_contribution`.

### E.3 Smoke

- `pytest tests/test_pipeline_smoke.py` — `rigel quant` on a minimal
  scenario BAM produces a non-empty `quant.feather`.

### E.4 Acceptance gate

- [ ] All `tests/calibration/test_locus_prior_consumer.py` pass.
- [ ] `pytest tests/test_pipeline_smoke.py` passes.
- [ ] `rigel quant` runs end-to-end on at least one scenario BAM.
- [ ] No `NotImplementedError` raises from any calibration / locus
      prior path on the smoke BAM.

---

## 9. Phase F — Validate

**Goal.** Per
[`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md):
prove the rebuild fixes the paralog and hybrid-capture cases without
regressing on standard sweeps.

### F.1 Layer 3 — Scenario tests

- **Paralog** — `tests/scenarios_aligned/test_paralogs.py`. Expected
  $t_1 \approx t_2 \approx 500$
  ([`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §3.1).
- **Hybrid-capture** — new scenario per
  [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §3.2.

### F.2 Layer 4 — Synthetic benchmark sweeps

Per [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §4.

Run `scripts/benchmark/configs/locus_simple_*.yaml` on the
post-rebuild binary, compare against a clean snapshot of the
pre-burn baseline (capture this snapshot at the **start** of Phase A
into `scratch/preburn/` so Phase F has a comparison point).

Acceptance:
- mRNA recovery: median relative error < 0.15.
- nRNA siphon reduced from baseline ~37% to < 20% worst case.
- gDNA recovery at low contamination: relative error < 1.3×.

### F.3 Layer 4 — Armis2 real-data smoke

Per [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §5.

Conditions:
- `gdna_none_ss_1.00_nrna_none`
- `gdna_none_ss_0.50_nrna_rand`
- `gdna_high_ss_1.00_nrna_none`
- `gdna_high_ss_0.50_nrna_rand`

Acceptance per §5.3.

### F.4 Numerical-stability stress

`tests/calibration/test_numerical_stability.py` per
[`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §6.

### F.5 Acceptance gate

- [ ] Paralog test passes naturally (no tolerance loosening).
- [ ] Hybrid-capture scenario test passes.
- [ ] Synthetic sweeps acceptance thresholds met (or any regression
      root-caused and documented).
- [ ] Armis2 corner conditions acceptance per §5.3.
- [ ] All `CalibrationResult` fields finite across all runs.
- [ ] Mass-change history monotone in every run.

### F.6 Validation report

Write `docs/acc_caljointmodel/validation_report.md` per
[`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §7.

---

## 10. Phase G — Cleanup

**Goal.** Mechanical: regenerate goldens, ruff, audit, doc updates.

### G.1 Tasks

1. Regenerate `tests/golden/*` with `pytest tests/ --update-golden`.
2. `ruff check src/ tests/` and `ruff format src/ tests/`.
3. Magic-number audit across the whole calibration module: ≤ 8
   literals, each annotated with a comment citing the relevant spec
   doc.
4. Update `CLAUDE.md` and `.github/copilot-instructions.md`:
   - "Python Module Roles" table: remove deleted v5 modules, add
     `scan_payload.py` (AccumulatorPayload), `calibrate.py`.
   - Architecture section: reflect the substrate (per-ref
     `AccumulatorPayload`) and the calibrator's three-substrate-set
     loop.
5. Update `docs/accumulator/00_design.md` §7.3: replace the stale
   `region_fl_hist` / `boundary_fl_hist` fields in the
   `AccumulatorPayload` example with the schema actually shipped in
   Phase B (no FL).
6. Mark `docs/accumulator/01_implementation.md` and
   `docs/caljointmodel/05_implementation_plan.md` as superseded by
   this document (front-matter banner; do not delete — they retain
   historical context).
7. Write `docs/acc_caljointmodel/postmortem.md` with what worked,
   what didn't, and what to do differently next time.
8. Archive retention decision for `archive/calibration_legacy_2026_05/`.
   Recommendation: keep until the next minor release.

### G.2 Acceptance gate

- [ ] `pytest tests/` all green.
- [ ] `ruff check src/ tests/` clean.
- [ ] Magic-number audit ≤ 8, each annotated.
- [ ] Docs updated.

---

## 11. File-touch matrix

Rows = files; columns = phases. `**X**` = create; `X` = edit;
`del` = delete; `mv` = move/rename.

| File | A | B | C | D | E | F | G |
|---|---|---|---|---|---|---|---|
| `src/rigel/native/calibration/accumulator.{h,cpp}` | — | edit (`AccumulatorSet` + thread-safe merge) | — | — | — | — | — |
| `src/rigel/native/bam_scanner.cpp` | strip dead `obs_*` locals | replace `set_regions` → `set_region_edges`; deposit at fragment site; emit payload | — | — | — | — | — |
| `src/rigel/native.py` | — | re-export `AccumulatorPayload` symbols (if any new) | — | — | — | — | — |
| `src/rigel/_accumulator.py` | — | possibly extend with payload helpers | — | — | — | — | — |
| `src/rigel/scan_payload.py` | del (legacy) | **add** (`AccumulatorPayload`, `CalibrationSubstrate`) | extend | — | — | — | — |
| `src/rigel/calibration/_stub.py` | **add** | — | del | — | — | — | — |
| `src/rigel/calibration/calibrate.py` | — | — | **add** (placeholder) | edit (real impl, 3 sub-PRs) | — | — | — |
| `src/rigel/calibration/__init__.py` | rewrite to re-export 3 stub symbols | — | rewrite to re-export real types | — | — | — | — |
| `src/rigel/calibration/*.py` (everything else under v5) | **del** | — | — | — | — | — | — |
| `src/rigel/native/calibration/region_index.h`, `region_signature.h` | conditionally del (grep first) | conditionally del | — | — | — | — | — |
| `src/rigel/pipeline.py` | strip v5 block; call stub | wire `AccumulatorPayload` build into scan | wire substrate into `calibrate(...)` | — | wire `CalibrationResult` into locus prior | — | — |
| `src/rigel/estimator.py` | drop calibration imports | — | — | — | — | — | — |
| `src/rigel/locus.py::compute_locus_priors_from_partitions` | replace with `NotImplementedError` stub | — | — | — | rewrite per doc 04 §6 | — | — |
| `src/rigel/config.py` | strip v5 sub-configs; add stub `CalibrationConfig` | — | replace stub with real `CalibrationConfig` | — | — | — | — |
| `src/rigel/cli.py` | strip v5 flags | — | — | — | — | — | — |
| `src/rigel/index.py` | — | add `build_region_partition` helper (edges, ref_offsets) | — | — | — | — | — |
| `archive/calibration_legacy_2026_05/` | **mv** salvageable files here | — | — | — | — | — | — |
| `tests/test_*.py` (v5 tests) | **del** ~28 files | — | — | — | — | — | — |
| `tests/native/test_accumulator_spec.py` | — | — | — | — | — | — | — |
| `tests/test_scanner_accumulator_integration.py` | — | **add** | — | — | — | — | — |
| `tests/test_accumulator_payload.py` | — | **add** | — | — | — | — | — |
| `tests/calibration/test_config_defaults.py` | — | — | **add** | — | — | — | — |
| `tests/calibration/test_result_invariants.py` | — | — | **add** | — | — | — | — |
| `tests/calibration/test_substrate_invariants.py` | — | — | **add** | — | — | — | — |
| `tests/calibration/test_g2_g3_deconvolution.py` | — | — | — | **add** (D1) | — | — | — |
| `tests/calibration/test_g4_exposure.py` | — | — | — | **add** (D2) | — | — | — |
| `tests/calibration/test_m_step.py`, `test_outer_loop.py`, `test_e2e_synthetic.py` | — | — | — | **add** (D3) | — | — | — |
| `tests/calibration/test_locus_prior_consumer.py` | — | — | — | — | **add** | — | — |
| `tests/calibration/test_numerical_stability.py` | — | — | — | — | — | **add** | — |
| `tests/scenarios_aligned/test_hybrid_capture.py` | — | — | — | — | — | **add** | — |
| `tests/scenarios_aligned/test_paralogs.py` | — | — | — | — | — | edit (expected outputs) | — |
| `tests/test_pipeline_smoke.py` | — | — | — | — | edit | — | — |
| `tests/golden/*` | — | — | — | — | — | — | regenerate |
| `CHANGELOG.md` | edit | — | — | — | — | — | edit |
| `CLAUDE.md` | — | — | — | — | — | — | edit |
| `.github/copilot-instructions.md` | — | — | — | — | — | — | edit |
| `docs/accumulator/00_design.md` | — | — | — | — | — | — | edit (§7.3 schema) |
| `docs/accumulator/01_implementation.md` | — | — | — | — | — | — | mark superseded |
| `docs/caljointmodel/05_implementation_plan.md` | — | — | — | — | — | — | mark superseded |
| `docs/caljointmodel/04_interface_contract.md` | — | — | — | — | — | — | drop "Status (2026-05-29)" pending banner |
| `docs/acc_caljointmodel/validation_report.md` | — | — | — | — | — | **add** | — |
| `docs/acc_caljointmodel/postmortem.md` | — | — | — | — | — | — | **add** |

---

## 12. Risk register

| Phase | Top risk | Mitigation |
|---|---|---|
| A | Hidden cross-module import couples a "deleted" module to something live | Run `pytest --collect-only` after every deletion batch; fix imports as they surface |
| A | The pre-burn synthetic baseline snapshot is forgotten until Phase F needs it | Capture it as the **first** action of Phase A, before any deletion |
| B | Per-worker `AccumulatorSet` merge has data race on the shared buffer | Worker-local sets are merged only after all workers join; same pattern as `stats_.merge_from` |
| B | The region partition builder doesn't exist after Phase A deleted `region_partition.py` | Re-implement minimally in `index.py` per Phase B.5; spec lives in this doc |
| C | Boundary substrate's count signal: `flux` vs `mass`? | Open question flagged in C.2; resolve at Phase C kickoff and lock the decision in doc 04 §2 |
| D1 | Three-leg paralog rescue empirically weaker than the prior FL-bearing design | Phase F paralog regression must verify. If it fails: (a) reintroduce FL only for paralog-flagged regions, (b) tighten BB priors, (c) document partial regression. See [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md) §2.5 |
| D3 | Per-region BB strand likelihood unstable at very small $\rho_d^{\text{BB}}$ | `scipy.betabinom.logpmf` handles the limit; `_BB_FLOOR` Newton bound prevents drift |
| D3 | Mass-change diagnostic increases between iterations (EM violation) | Indicates bug. Raise `CalibrationConvergenceError`; do not silently continue |
| E | Locus prior pseudocount underflows or overflows on extreme regions | Consumer unit tests cover extreme-mass cases; clip $\alpha_t \in [\alpha_{\min}, \alpha_{\max}]$ |
| F | Regression in armis2 sweeps with no clear root cause | `scripts/debug/dump_calibration_state.py` (write in Phase D) is the debugger |

---

## 13. Sign-off

Approval to begin **Phase A** = approval of this consolidation. Each
phase gate is its own checkpoint and is mergeable as a standalone
PR. No phase begins until the previous gate is green.
