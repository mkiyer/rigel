# M5 — EM `prior_weight_rna` ABI + `Locus` / `MultiLocus` Rename

**Predecessor:** [calibration_v6_plan.md §M5](calibration_v6_plan.md).
M4 should land first so that this milestone's EM ABI extension can be
benchmarked against an unchanged density layer.

**Goal:** two structural changes that together unlock M6 (per-`Locus`
gDNA mass) without re-naming churn during the numerical work:

1. **EM ABI:** native solver accepts continuous nRNA suppression via a
   per-component `prior_weight_rna` weight (default = all-ones, identical
   to current behaviour).
2. **Locus nomenclature:** the existing `Locus` class is renamed to
   `MultiLocus` (because it is a connected component of transcripts —
   the EM-execution unit, not a single contiguous interval), and a new
   tiny `Locus` dataclass represents one contiguous genomic interval
   (the calibration-estimation unit).

This is **invasive** but **mechanical**: zero numerical change at
default settings (`prior_weight_rna = None` ≡ all ones), plus a
codebase-wide find-and-replace for the rename.  No legacy alias.  No
deprecation warning.  Hard cut.

---

## 1. Rationale

### 1.1 Why split `Locus` into `MultiLocus` + `Locus`

The current [src/rigel/scored_fragments.py](src/rigel/scored_fragments.py) `Locus` carries
`merged_intervals: list[tuple[str, int, int]]` — the maximal merged
genomic spans of its constituent transcripts.  Most components have
exactly one merged interval; paralog clusters across references have 2+.

M6 needs **per-interval** gDNA estimates because the gDNA density is a
property of contiguous genomic territory, not of a transcript-graph
component.  Continuing to use `merged_intervals: list[tuple]` for this
forces every M6 call site to pay the unpacking tax and ad-hoc validate
whether the "tuple" is `(str, int, int)`.  A first-class type fixes this
in one stroke.

The new contract:

* **`Locus`** (new, frozen dataclass): one `(ref, ref_id, start, end)`
  contiguous interval.  The unit M6 estimates gDNA mass on.
* **`MultiLocus`** (renamed from `Locus`): one connected component of
  transcripts linked by shared fragments.  The unit the EM is run on.
  Carries `loci: tuple[Locus, ...]` instead of `merged_intervals`.

### 1.2 Why the `prior_weight_rna` ABI

Under the M6 + M8 design, every locus EM has at least two RNA components
(mRNA + nRNA per transcript) and one gDNA component.  The default policy
suppresses nRNA by giving it `prior_weight_rna = 0` while mRNA keeps
weight 1.  But hard-coding a binary "is this component nRNA?" decision
in C++ leaks knowledge of the EM component layout into the solver.

A continuous `float[n_components]` weight vector keeps the solver
agnostic and lets future extensions (e.g., per-isoform priors, weighted
mRNA grouping) pass arbitrary fractional weights without an ABI break.
Default `nullptr` ≡ all ones preserves bit-identical EM output for
everything that doesn't opt in.

---

## 2. Part A — `prior_weight_rna` ABI

### 2.1 C++ surface

In [src/rigel/native/em_solver.cpp](src/rigel/native/em_solver.cpp), the
`compute_ovr_prior_and_warm_start` static helper currently takes:

```cpp
static void compute_ovr_prior_and_warm_start(
    const std::vector<EmEquivClass>& ec_data,
    const double* unambig_totals,
    const double* eligible,
    double        alpha_gdna,
    double        alpha_rna,
    int           gdna_idx,
    double*       prior_out,
    double*       theta_init_out,
    int           n_components,
    bool          use_vbem)
```

Add **one** parameter immediately after `gdna_idx`:

```cpp
const float*  prior_weight_rna,  // nullptr ⇒ all ones (default)
```

Implementation change (the only places that touch `coverage_totals`):

```cpp
// total_rna_coverage now weights each RNA component
double total_rna_coverage = 0.0;
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] > 0.0 && i != gdna_idx) {
        const float w = prior_weight_rna ? prior_weight_rna[i] : 1.0f;
        total_rna_coverage += w * coverage_totals[i];
    }
}

// per-component prior fan-out
for (int i = 0; i < n_components; ++i) {
    if (eligible[i] <= 0.0) {
        prior_out[i] = 0.0;
    } else if (i == gdna_idx) {
        prior_out[i] = baseline + std::max(alpha_gdna, EM_LOG_EPSILON);
    } else if (total_rna_coverage > 0.0) {
        const float w = prior_weight_rna ? prior_weight_rna[i] : 1.0f;
        prior_out[i] = baseline + std::max(
            alpha_rna * w * coverage_totals[i] / total_rna_coverage,
            EM_LOG_EPSILON);
    } else if (n_rna_eligible > 0) {
        prior_out[i] = baseline + std::max(alpha_rna / n_rna_eligible,
                                           EM_LOG_EPSILON);
    } else {
        prior_out[i] = baseline;
    }
}
```

The warm-start block (`alpha_gdna / alpha_rna * others`) is **not**
modified.  Warm-start is mass-based, not per-component.

**Why `float` and not `double`:** the weights are user-policy values in
`[0, 1]` with no need for double precision.  `float` halves the memory
footprint of the per-locus weight vector (matters at batch scale: 200K
loci × 1KB → 200MB).

### 2.2 nanobind exposures

Both `run_locus_em_native` and `batch_locus_em_partitioned` gain one
new argument.  The pattern matches existing optional ndarray args:

```cpp
m.def("run_locus_em_native", &run_locus_em_native,
      ...,
      nb::arg("prior_weight_rna").none() = nb::none(),
      "...");

m.def("batch_locus_em_partitioned", &batch_locus_em_partitioned,
      ...,
      nb::arg("locus_prior_weight_rna").none() = nb::none(),
      "...");
```

`run_locus_em_native` accepts `np.ndarray[float32, (n_components,)] |
None`.  `batch_locus_em_partitioned` accepts `list[np.ndarray[float32]]
| None`; if `None`, every locus is treated as all-ones.  If a
non-`None` list is shorter than `n_loci`, hard error (caller bug).

### 2.3 Python call sites

* `run_locus_em_native` is called from `src/rigel/native.py` thunks
  used by `tests/test_em_impl.py::_make_locus`.  Add the new optional
  kwarg with default `None`; tests pass through unchanged.
* `batch_locus_em_partitioned` is called from `src/rigel/scan.py` (the
  per-batch runner used by `pipeline.quant_from_buffer`) and possibly
  from `src/rigel/estimator.py`.  In both, plumb a new optional kwarg
  `locus_prior_weight_rna=None` from the caller down to the C++ entry.
  No call site sets it in this milestone — that is M8b's job after
  M6 produces `PriorTable.prior_weight_rna`.

### 2.4 The `build_prior_weight_rna` helper

```python
def build_prior_weight_rna(
    multi_locus: MultiLocus,
    em_data: ScoredFragments,
    *,
    nrna_weight: float = 0.0,
) -> np.ndarray:
    """Construct the per-component nRNA-suppression weight vector.

    Components are laid out by the EM as `[mRNA_t0, nRNA_t0, mRNA_t1,
    nRNA_t1, ..., gDNA]`.  This helper returns a `float32` array of
    length `2 * n_t + 1` with mRNA entries set to 1.0, nRNA entries
    set to ``nrna_weight``, and the gDNA entry set to 1.0 (the gDNA
    component is not affected by ``prior_weight_rna``; the solver
    routes it through ``alpha_gdna`` / ``gdna_idx``).
    """
```

Lives in `src/rigel/calibration/locus_prior.py` (created in M6); for M5
we ship a stub that returns `np.ones(2 * n_t + 1, dtype=np.float32)`
plus a TODO comment pointing at M6 / `nrna_weight`.

### 2.5 Tests — `tests/test_em_prior_weight.py` (≥ 5 cases)

* `test_none_argument_bit_identical_to_baseline` — pass `None` vs omit;
  EM output bit-identical (master regression guard).
* `test_all_ones_bit_identical_to_none` — pass an explicit
  `np.ones(n_components, dtype=np.float32)`; EM output bit-identical.
* `test_full_nrna_suppression` — set every nRNA component's weight to
  `0.0`, mRNA to `1.0`; assert EM output's nRNA components have
  `theta < EM_LOG_EPSILON_DEPENDENT_BOUND` after convergence.
* `test_half_share` — mRNA weight 1.0, nRNA weight 0.5; assert
  posterior mass ratio between mRNA and nRNA shifts by the expected
  factor in a controlled-coverage scenario.
* `test_vbem_baseline_preserved` — same input through MAP and VBEM
  modes with and without weights; both modes respect the weights
  identically.
* `test_degenerate_all_zero_weights` — every RNA weight = 0.0;
  `total_rna_coverage == 0` branch fires; EM falls back to uniform RNA
  prior and converges without NaN.

---

## 3. Part B — `Locus` / `MultiLocus` rename

### 3.1 New `Locus` dataclass

In [src/rigel/scored_fragments.py](src/rigel/scored_fragments.py):

```python
@dataclass(frozen=True, slots=True)
class Locus:
    """One contiguous genomic interval — the calibration-estimation unit.

    A ``MultiLocus`` is composed of one or more ``Locus`` intervals.
    Most ``MultiLocus``es have exactly one ``Locus``; paralog clusters
    across references have several.
    """
    ref:    str
    ref_id: int
    start:  int
    end:    int

    @property
    def span(self) -> int:
        return self.end - self.start
```

`ref_id` is carried alongside `ref` to avoid every `Locus` consumer
re-translating via `index.resolver.get_ref_to_id()`.

### 3.2 Renamed `MultiLocus`

The existing `Locus` class is renamed `MultiLocus` and gains
`loci: tuple[Locus, ...]` in place of `merged_intervals`:

```python
@dataclass(slots=True)
class MultiLocus:
    """Connected component of transcripts linked by shared fragments.

    The unit the EM is run on.  Composed of one or more ``Locus``
    intervals (the calibration-estimation unit).
    """
    multi_locus_id:     int           # was: locus_id
    transcript_indices: np.ndarray    # int32, global transcript IDs
    unit_indices:       np.ndarray    # int32, EM unit row indices
    gdna_span:          int           # Σ (l.end - l.start) for l in loci
    loci:               tuple[Locus, ...]
```

* `merged_intervals` is **removed** entirely.  No alias.
* `locus_id` → `multi_locus_id`.  No alias.
* `gdna_span` is now redundant with `sum(l.span for l in loci)` but is
  retained as a precomputed cache because EM hot paths read it; it is
  asserted-equal to the sum at construction time.

### 3.3 `build_loci` becomes `build_multi_loci`

[src/rigel/locus.py](src/rigel/locus.py) `build_loci()` is renamed to
`build_multi_loci()` and emits `MultiLocus` objects whose `loci` field
is constructed from the merge loop output:

```python
loci_list: list[Locus] = []
for ref_name, s, e in merged:
    loci_list.append(Locus(
        ref=ref_name,
        ref_id=index.resolver.get_ref_to_id()[ref_name],
        start=s, end=e,
    ))

multi_loci.append(MultiLocus(
    multi_locus_id=lid,
    transcript_indices=t_idx,
    unit_indices=...,
    gdna_span=span,
    loci=tuple(loci_list),
))
```

The old `merged_intervals=[(ref, s, e), ...]` construction is deleted
in the same edit.

### 3.4 Call-site migration

The following files reference `Locus` or `merged_intervals` and must be
updated in lockstep:

| File | Change |
|---|---|
| `src/rigel/scored_fragments.py` | Define new `Locus`, rename old to `MultiLocus`, drop `merged_intervals` |
| `src/rigel/locus.py` | `build_loci` → `build_multi_loci`; emit `loci=tuple(...)` |
| `src/rigel/locus_partition.py` | Update import + `LocusPartition` doc references |
| `src/rigel/pipeline.py` | `from .locus import build_multi_loci, ...`; rename local `loci` → `multi_loci` |
| `src/rigel/scan.py` | Update type hints / iteration variables |
| `src/rigel/estimator.py` | Update type hints / iteration variables |
| `tests/conftest.py` | `Locus(...) → MultiLocus(...)`; `merged_intervals=[("chr1", 0, 10000)]` → `loci=(Locus(ref="chr1", ref_id=0, start=0, end=10000),)` |
| `tests/test_locus_partition.py` | same |
| `tests/test_em_impl.py` | uses `_make_locus` helper; verify it doesn't reference `merged_intervals` |
| Any other test that constructs `Locus(...)` or reads `.merged_intervals` | grep + fix |

The migration is mechanical: one PR, one commit, no compatibility
shim.  Reviewer-friendly because the diff is concentrated and
type-checker (mypy / pyright) flags every missed site.

### 3.5 Why no deprecation alias

Per the user's mandate ("no legacy code, no backwards compatibility"):
shipping a deprecation alias `Locus = MultiLocus` for one release would
mask exactly the bugs M6 needs to catch (per-`Locus` vs per-`MultiLocus`
confusion).  Hard cut now; everyone reading `Locus` going forward gets
the contiguous-interval semantics.

### 3.6 Tests — `tests/test_locus_rename.py` (≥ 3 cases)

* `test_single_locus_multi_locus` — two-exon, single-ref transcript
  produces a `MultiLocus` with `len(loci) == 1`.
* `test_multi_ref_multi_locus` — paralog cluster across `chr1` + `chr2`
  produces a `MultiLocus` with `len(loci) == 2`, `loci` sorted by
  `(ref_id, start)`.
* `test_gdna_span_matches_sum_of_loci` — assert `multi_locus.gdna_span
  == sum(l.span for l in multi_locus.loci)` for every `MultiLocus`
  produced by `build_multi_loci` on the `mini_index` fixture.

---

## 4. Cleanup opportunities (collateral wins)

While reading these surfaces I noticed three items worth doing in this
milestone (each touches the same files, so cost is incremental):

1. **`gdna_span: int` is computed during `build_loci` and never
   updated.**  Make `MultiLocus` frozen too — pass `gdna_span=sum(l.span
   for l in loci)` at construction.  Frozen + slots is the consistent
   pattern across the new calibration result types (M4 + M6 + M7).
   Mutability of `MultiLocus` is not used anywhere (verified by grep).
2. **`scored_fragments.py` mixes EM-data CSR types and the locus
   container.**  The `Locus` / `MultiLocus` types belong with the locus
   builder, not in the CSR module.  Move `MultiLocus` and the new
   `Locus` to `src/rigel/locus.py` (next to `build_multi_loci`); leave
   `LocusPartition` (which is a CSR view) where it is.  All call sites
   already say `from .locus import build_loci`; the new
   `from .locus import MultiLocus, Locus` is the natural follow-up.
3. **Drop `Locus` from `scored_fragments.py`'s `__all__`** (it's not in
   one today; pin it isn't there post-rename either).

Items 1 and 2 ship in the rename commit; item 3 is implicit.

---

## 5. Commit ordering

Two commits, in this order:

### 5.1 Commit `m5-em-abi` — Part A only

* C++ `compute_ovr_prior_and_warm_start` gains `const float*
  prior_weight_rna` parameter.
* Both nanobind entry points (`run_locus_em_native`,
  `batch_locus_em_partitioned`) gain the optional kwarg.
* Python thunks in `src/rigel/native.py` and call sites in `scan.py` /
  `estimator.py` get the optional pass-through.
* Tests: `tests/test_em_prior_weight.py` (5 cases).
* **Exit gate:** all existing EM tests bit-identical with `None`
  argument; new tests green.

### 5.2 Commit `m5-locus-rename` — Part B only

* New `Locus` + renamed `MultiLocus` in `src/rigel/locus.py` (moved out
  of `scored_fragments.py`).
* `build_multi_loci` emits the new types.
* All call sites updated (single mechanical pass).
* Tests: `tests/test_locus_rename.py` (3 cases) + every existing test
  that constructed `Locus(...)` with `merged_intervals=` is updated to
  the new constructor.
* **Exit gate:** full test suite green; `grep -rn "merged_intervals"
  src/ tests/` returns zero hits; `grep -rn "from .scored_fragments
  import.*Locus" src/` returns zero hits (Locus moved to `.locus`).

Splitting reduces blast radius — if Part A regresses an EM test, we
revert one commit, not the whole milestone.

---

## 6. Combined exit gate

* `tests/test_em_prior_weight.py` ≥ 5 green.
* `tests/test_locus_rename.py` ≥ 3 green.
* Full `pytest tests/` green (969+8 = 977+ passing, the 21
  pre-existing golden-output failures unchanged).
* EM `tests/test_em_impl.py` byte-identical against pre-M5 outputs
  when caller passes `None` for `prior_weight_rna`.
* `grep -rn "merged_intervals\|build_loci\b" src/ tests/` returns zero
  hits.
* `grep -rn "Locus" src/rigel/scored_fragments.py` returns zero hits.

---

## 7. What this milestone deliberately does NOT do

* **Does not** populate `prior_weight_rna` from any calibration source.
  The wiring `(MultiLocusPrior → np.ndarray → batch_locus_em_partitioned)`
  is M6 + M8b.  M5 ships the ABI only.
* **Does not** split the EM gDNA component per-`Locus`.  One gDNA
  component per `MultiLocus` remains.  Per-`Locus` decomposition is
  noted as a future-plan deferral in the parent §8.
* **Does not** introduce `MultiLocusPrior` or `LocusGdnaEstimate` —
  those are M6 types.
* **Does not** migrate any SRD-v1 calibration path.  M5 is independent
  of the calibration data flow; it touches the EM solver and the locus
  type only.
