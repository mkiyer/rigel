# gDNA Effective Length Normalization Plan v2

Status: v2 implementation plan, code-checked, not yet applied

This plan consolidates the deferred gDNA work from
`effective_length_consistency_2026_05.md`, Phase 6 of
`effective_length_implementation_plan_2026_05.md`, and the current audit in
`gdna_effective_length_audit_2026-05-15.md`.

The target change is narrow but important: move gDNA effective-length
normalization into the same EM contract used by RNA components. The first
implementation should be geometry-only. It does not need mappability,
capture-panel, or alignment-access weights, but it should be shaped so those
weights can be inserted later without another redesign.

Naming convention: use `gdna_eff_len` everywhere new code needs this quantity.
That matches the existing `log_eff_len` vocabulary in `em_solver.cpp` and keeps
grep results clean. In formulas below, `Ltilde_gdna(M)` is the mathematical
object and `gdna_eff_len[M]` is its implementation field.

## V2 Code-Checked Changes

The review feedback was checked against the current codebase and incorporated
with the following decisions:

- Accepted: `ref_lengths_arr` must be built unconditionally in
   `assemble_priors`; the current code builds it only under
   `intergenic_flank_bp > 0`, but `gdna_eff_len_for_loci` needs it for every
   assembly call.
- Accepted: the native EM ABI migration should be a single-step required C++
   positional argument plus a Python-wrapper `None -> ones` default. Tests call
   `AbundanceEstimator.run_batch_locus_em_partitioned`, not the native binding
   directly, so the C++ deprecation/default path is unnecessary.
- Accepted: removing the scorer-side `e_h` denominator makes both `gdna_flank`
   and `t_span_` dead scoring state. Remove them in the same implementation
   series.
- Accepted with caveat: the scorer currently emits one replicated
   `gdna_log_lik(u)` scalar per EM unit. Do not refactor the scorer to
   per-`MultiLocus` `H_{u,M}` subsets in the first patch; keep the existing
   scalar replication semantics and test them.
- Accepted: add `gdna_eff_len` diagnostics to `loci.feather` and summary
   aggregate statistics under `summary.json`, not to `quant.feather` or
   `gene_quant.feather`.
- Accepted cleanup: retire the deprecated `n_gdna_exon_intron` compatibility
   aggregate after the effective-length behavior is validated.

## Decision

Use a per-`MultiLocus` FL-marginal gDNA effective length based on overlap
support:

```text
P(f | gDNA_M) ∝ h_G(ell_f) q_G(f) / Ltilde_gdna(M)
Ltilde_gdna(M) = sum_ell h_G(ell) N_valid_starts(M, ell)
```

where `M` is the `MultiLocus` EM component, `h_G` is the gDNA fragment-length
PMF, `q_G(f)` is the non-length gDNA score term, and
`N_valid_starts(M, ell)` is the number of genome-valid fragment start positions
whose length-`ell` fragment overlaps the `MultiLocus` transcript-footprint
union.

This FL-marginal rewrite is exact, not an approximation. Under the standard
conditional generative model,

```text
P(ell, s | M, gDNA) = h_G(ell) I[s is valid for (M, ell)] / Ltilde_gdna(M)
P(ell | M, gDNA)   = h_G(ell) N_valid_starts(M, ell) / Ltilde_gdna(M)
P(s | ell, M)      = 1 / N_valid_starts(M, ell)
```

Multiplying `P(ell | M)` and `P(s | ell, M)` cancels
`N_valid_starts(M, ell)`, leaving the per-fragment likelihood contribution
`log h_G(ell_f) - log Ltilde_gdna(M)` plus non-length score terms. Moving
`-log Ltilde_gdna(M)` from the scorer into EM component metadata therefore
changes the parameterization, not the intended probability model.

This replaces the current per-hit scorer correction:

```text
-log(max(t_span(anchor_t) + gdna_flank - genomic_footprint + 1, 1))
```

The current correction is real length normalization, but it is hidden in the
scorer, anchored to one candidate transcript, and parameterized differently
from RNA. It also uses a mean-flank containment-style denominator rather than
the actual overlap support of the gDNA component that enters EM.

## Why Overlap Support

gDNA fragments are genomic intervals. A fragment can support the gDNA component
for a transcript locus even when its start or end lies outside the transcript
footprint, as long as the aligned fragment overlaps the footprint and is routed
into the locus candidate set. Therefore the geometric opportunity space should
count overlap starts, not contained starts.

> RNA containment and gDNA overlap are intentionally asymmetric. RNA fragments
> are generated from transcript sequence, so their observed insert must lie
> within the transcript's spliced molecule; bases outside that molecule are not
> part of the RNA template and would be clipped or absent. gDNA fragments are
> genomic intervals, so they can extend past either transcript-footprint
> boundary by up to `ell - 1` bases and still be valid evidence for the local
> gDNA component if the aligned fragment overlaps the locus and is routed into
> its candidate set.

Coordinates are half-open, and integer starts are 0-based. For one locus
interval `[a, b)` and a fragment length `ell`, a fragment starting at `s`
overlaps the locus when:

```text
s < b and s + ell > a
```

Equivalently, valid starts are:

```text
s in [a - ell + 1, b)       # half-open start-position window
```

or `a - ell + 1` through `b - 1` inclusive. Away from contig ends, the count is:

```text
(b - a) + ell - 1
```

This agrees with the intuition that possible overlapping fragments can extend
outside both locus boundaries. The coordinate nuance is that
`locus_end + ell - 1` describes the farthest right base touched by a possible
overlapping fragment, not the farthest valid start position. The denominator
counts start positions.

For a single interval away from contig ends:

```text
Ltilde_gdna([a,b)) = sum_ell h_G(ell) ((b - a) + ell - 1)
                   = (b - a) + E_G[ell] - 1
```

There is no separate `gdna_flank` term in this model. The boundary opportunity
is created directly by the fragment-length marginalization.

## MultiLocus Geometry

The native EM has one gDNA component per `MultiLocus`, not per transcript and
not per contiguous `Locus` interval. Therefore the effective length must be
computed for the `MultiLocus` genomic footprint union.

`MultiLocus.loci` already stores merged non-overlapping intervals, sorted by
`(ref_id, start)`, and `MultiLocus.gdna_span` is the total merged footprint
span. For a multi-interval locus, the exact overlap denominator is not always
`gdna_span + E[ell] - 1`. Instead, for each fragment length `ell`:

1. Convert every interval `[a_i, b_i)` on the same reference into a valid-start
   window `[a_i - ell + 1, b_i)`.
2. Clip each window to valid genomic starts `[0, ref_len - ell + 1)`.
3. Merge overlapping start windows per reference.
4. Sum the merged lengths across references.

In formula form:

```text
N_valid_starts(M, ell) = sum_ref | union_{i in M, ref_i=ref}
    ([a_i - ell + 1, b_i) ∩ [0, ref_len(ref) - ell + 1)) |
```

The union step matters. If two loci on the same contig are separated by less
than `ell`, a single start position can produce a fragment overlapping both;
that start should be counted once for the shared gDNA component.

## Relationship To RNA Effective Length

RNA currently follows the desired contract:

```text
score_RNA(f,t) = log h_R(ell_f) + log q_R(f,t)
EM uses       = log theta_t - log Ltilde_t + score_RNA(f,t)
```

gDNA should follow the same shape:

```text
score_gDNA(f,M) = log h_G(ell_f) + log q_G(f,M)
EM uses         = log theta_g - log Ltilde_gdna(M) + score_gDNA(f,M)
```

The practical rule is simple: fragment-length PMF terms stay in the scorer;
effective-length denominators live in EM component metadata.

## Expected EM Impact

This should be framed as a material model fix, not only as infrastructure. The
old scorer-side denominator is:

```text
e_old(h) = max(t_span(anchor_t) + gdna_flank - genomic_footprint_h + 1, 1)
```

For the unstranded RCA's locus-15-style cases, `e_old` can be `1` to a few
hundred because the anchor transcript span is one transcript's genomic span and
the genomic footprint of an unspliced/nRNA-like fragment can approach that span.
The new denominator is approximately:

```text
Ltilde_gdna(M) ≈ locus_span(M) + E_G[ell] - 1
```

For a 40 kb multi-isoform locus, that is roughly `40,000`. The per-fragment
shift in the gDNA assignment objective is:

```text
delta log P_gDNA ≈ log(e_old / Ltilde_gdna(M))
```

So `e_old = 1, 100, 500` against `Ltilde_gdna = 40,000` gives shifts of about
`-10.6`, `-6.0`, and `-4.4` nats per fragment, respectively. The sign is
unambiguous: removing the old tiny anchor-transcript denominator should reduce
gDNA posterior odds in exactly the short-anchor/large-locus cases where nRNA
was being siphoned. Calibration prior strength may still matter, but this fix
has a real chance to close a substantial fraction of the unstranded-nRNA gap on
its own.

Before changing calibration priors, run the synthetic 24 benchmark with this
fix alone and measure how much of the unstranded nRNA loss is recovered. This
first measurement should determine how aggressive later prior retuning needs to
be.

## Relationship To Calibration Priors

No calibration-density denominator change is required in the first patch, but
that deferral must be measured rather than assumed safe.

The new `Ltilde_gdna(M)` is the likelihood normalizer for the local gDNA EM
component. The existing `gdna_prior_count[M]` is a calibration-derived prior
expected count, assembled from locoregional gDNA estimates and eligibility
rules. Keeping that prior unchanged initially isolates the likelihood
normalization fix and makes the benchmark delta interpretable, but the old
prior's practical leverage was tuned against the old gDNA likelihood scale. A
multi-nat per-fragment shift can make the same numeric prior behave weaker or
stronger relative to the likelihood.

The acceptance check for the first patch is therefore:

```text
sum_M E_EM[N_gdna(M)]  vs.  pi_pool * n_em_fragments
```

computed in fractional assignment mode for no-nRNA scenarios where the pool
calibration is already known to be reliable. If the EM gDNA total drifts
materially from the calibrated pool expectation, do not hide that under the
"no calibration change" label. Either re-derive the per-locus prior under the
new likelihood scale or introduce an explicit, benchmarked prior scale factor
as a second-stage change.

Longer term, the same overlap-support primitive should become the geometry
kernel for a fully unified exposure model:

```text
Ltilde_gdna(M) = sum_ell h_G(ell) sum_s W_M(s, ell)
```

where `W_M(s, ell)` can later encode mappability, capture-panel accessibility,
alignment survival, duplicate/read-start constraints, or category-specific
routing. In the initial implementation, `W_M(s, ell)` is just the 0/1 indicator
that `[s, s + ell)` is genome-valid and overlaps the `MultiLocus` footprint.

This keeps intergenic calibration concerns separate: capture data may make
intergenic density a poor prior source, but that is a prior-calibration problem,
not a reason to keep a per-hit anchor-transcript denominator in the scorer.

## Implementation Plan

### 1. Add a gDNA overlap exposure primitive

Add a public helper in `src/rigel/calibration/_exposure.py`:

```python
def gdna_eff_len_for_loci(
    loci: tuple[Locus, ...],
    ref_lengths: Mapping[str | int, int],
    fl: FragmentLengthModel,
    *,
    min_value: float = 1.0,
) -> float:
    ...
```

Implementation notes:

- Bundle the single-interval closed-form fast path in the first patch. If `loci` contains
  exactly one interval and the expanded start window for every nonzero `ell`
  stays within contig bounds, return `span + fl.mean - 1` floored by
  `min_value`.
- Fall back to the exact union path for contig-edge intervals,
  multi-interval loci, and cross-reference `MultiLocus` records.
- In the exact path, iterate only over fragment lengths with nonzero `fl.pmf`.
- For each `ell`, group start windows by `ref_id`, clip to
  `[0, ref_len - ell + 1)`, sort/merge, and accumulate union length.
- Weight the union length by `fl.pmf[ell]`.
- Return `max(total, min_value)` for EM use.
- Support `min_value=0.0` for future density/exposure callers.

In `assemble_priors`, hoist this construction out of the
`intergenic_flank_bp > 0` conditional:

```python
ref_lengths_arr = np.asarray(list(index.ref_lengths.values()), dtype=np.int64)
```

The current conditional was sufficient when reference lengths were only used to
clamp the diagnostic intergenic flank query, but `gdna_eff_len_for_loci` needs
reference lengths unconditionally. Once the array is unconditional, pass it to
both the existing intergenic-flank path and the new gDNA exposure helper.

The fast path matters because `assemble_priors` can touch tens of thousands of
loci, and most `MultiLocus` objects are single genomic intervals. The exact
path should still be the reference implementation and should be used whenever
edge clipping or interval union could change the answer.

Tests should cover:

- Single interval away from contig ends gives `span + mean - 1`.
- Intervals at the left and right contig edges are clipped correctly.
- Two far-apart intervals add independently.
- Two nearby intervals merge expanded start windows and avoid double-counting.
- Multi-reference loci are summed per reference.

### 2. Store per-locus gDNA effective lengths with priors

Extend `PriorTable` in `src/rigel/calibration/locus_prior.py`:

```python
gdna_eff_len: np.ndarray  # float64, shape (n_loci,)
```

Populate it in `assemble_priors` while iterating over `multi_loci`:

```python
gdna_eff_len_arr[idx] = gdna_eff_len_for_loci(
    ml.loci,
    index.ref_lengths,
    gdna_fl,
    min_value=1.0,
)
```

This field belongs next to `gdna_prior_count` and `enable_gdna` because all
three are per-`MultiLocus` component metadata consumed by the EM.

Keep `gdna_eff_len` flat-array-only on `PriorTable`; do not duplicate it onto
`MultiLocusPrior`. `MultiLocusPrior.gdna_prior_count` is already duplicated for
historical diagnostic reasons, but copying that pattern would add another sync
hazard. The EM needs only the array aligned by `multi_locus_id`.

`ref_lengths` plumbing is already available: `TranscriptIndex.load()` populates
`index.ref_lengths`, and `assemble_priors` already derives reference-length
arrays from it for per-locus estimates. Reuse that source of truth, indexed by
`Locus.ref_id`/`Locus.ref`, rather than adding a second reference-length path.

### 3. Thread the new array through Python EM dispatch

Add `gdna_eff_len` to the partitioned EM path:

- `pipeline._run_locus_em_partitioned(...)`
- `AbundanceEstimator.run_batch_locus_em_partitioned(...)`
- the call to native `_batch_locus_em_partitioned(...)`

Use a contiguous `float64` array and validate length equals `len(partition_tuples)`.
Production calls should pass `prior_table.gdna_eff_len` aligned by
`multi_locus_id`.

### 4. Add native EM ABI support with a Python-wrapper default

Extend `batch_locus_em_partitioned` in `src/rigel/native/em_solver.cpp` with a
new required positional argument immediately after `locus_enable_gdna`:

```cpp
f64_1d locus_gdna_eff_lens
```

Then:

- Capture `const double* gel_ptr = locus_gdna_eff_lens.data();`.
- Validate `locus_gdna_eff_lens.shape(0) == n_loci` before releasing the GIL.
- Add `double gdna_eff_len` to `extract_locus_sub_problem_from_partition`.
- Set `sub.log_eff_len[sub.gdna_idx] = log(max(gdna_eff_len, 1.0));`.
- Update comments that currently say gDNA `log Ltilde = 0` because scorer
  likelihoods are pre-normalized.
- Update nanobind argument registration.

This changes both the EM E-step and `assign_posteriors`, because both already
consume `sub.log_eff_len`.

To preserve existing Python-level tests, make the wrapper default live in
`AbundanceEstimator.run_batch_locus_em_partitioned`, not in C++:

```python
def run_batch_locus_em_partitioned(..., gdna_eff_len: np.ndarray | None = None, ...):
   if gdna_eff_len is None:
      gdna_eff_len = np.ones(len(partition_tuples), dtype=np.float64)
```

Production callers in `pipeline._run_locus_em_partitioned` should pass real
`prior_table.gdna_eff_len`; tests that exercise the wrapper can keep omitting
the new argument until they are updated. The codebase currently has no direct
test calls to `_batch_locus_em_partitioned`, so a C++ `None` fallback and
deprecation cycle would add complexity without buying compatibility.

### 5. Remove the scorer-side gDNA denominator

In `src/rigel/native/scoring.cpp`, remove the per-hit `-log(e_h)` term from
both gDNA paths:

- multimapper path using `anchor_t = cp.t_ind[start]`
- non-MM path using `best_t`

The scorer must still perform MM hit aggregation and NH normalization. The v2
implementation should preserve the current partitioner contract: the scorer
emits one per-unit scalar over the full multimapping group, and the partitioner
replicates that same scalar into every `MultiLocus` touched by that unit.

For each gDNA-eligible hit `h` in the full unit-level set `H_u`:

```text
B_h = log h_G(ell_h) + gdna_log_sp + LOG_HALF + log_nm_h + penalties_h
```

The emitted scalar should be:

```text
gdna_log_lik(u) = logsumexp_{h in H_u} B_h - log |H_u|
```

This is exactly the current online `logsumexp(...) - log(nh_gdna)` structure
with only the per-hit `-log(e_h)` removed. If every eligible hit has the same
`B_h`, then:

```text
gdna_log_lik(u) = B
```

That means multiple MM hits of the same fragment do not multiply the fragment's
mass. If a unit is split into two `MultiLocus` partitions, both partitions see
the same `gdna_log_lik(u)` scalar and `nh_gdna` reflects the full MM group, not
the per-component subset. Mass conservation across `MultiLocus` components is
therefore a partitioning/EC responsibility in the current architecture, not a
scorer-side `H_{u,M}` responsibility. Do not refactor the scorer toward
per-`MultiLocus` hit subsets in the first patch.

After the change, the non-MM single-hit gDNA log-likelihood should include:

```text
gdna_frag_len_log_lik(genomic_footprint)
+ gdna_log_sp
+ LOG_HALF
+ log_nm
+ penalties already intended for gDNA
```

It should not depend on anchor transcript span. Remove `gdna_flank` fully in
the same implementation series rather than leaving a no-op. The bounded churn
is in `scoring.cpp`, the `NativeFragmentScorer` nanobind constructor,
`FragmentScorer.from_models`, and the small set of scoring call sites/tests.
Keeping a no-op `gdna_flank` would invite accidental reintroduction of the old
model.

Dead code to remove in the same patch series:

- `src/rigel/native/scoring.cpp`: `t_span_`, `t_span_arr` constructor
   parameter, its assignment/copy block, and `nb::arg("t_span_arr")`.
- `src/rigel/native/scoring.cpp`: `gdna_flank_`, `gdna_flank` constructor
   parameter, nanobind argument, and all Option-B harmonic-mean comments.
- `src/rigel/scoring.py`: `FragmentScorer.t_span_arr`, the local
   `t_span_arr` derivation, and pass-throughs into `FragmentScorer` and
   `NativeFragmentScorer`.
- `src/rigel/scoring.py`: `gdna_flank` derivation and pass-through.

This removal is safe because RNA coordinate mapping uses `t_length_`,
`t_start_`, and exon CSR arrays; no non-gDNA scoring path reads `t_span_`.

After any C++ edit, rebuild with:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

### 6. Emit lightweight diagnostics

Assignment semantics stay the same, except that posterior odds can change
because the gDNA component now uses the correct denominator. Add two diagnostic
fields to `loci.feather` via `AbundanceEstimator.get_loci_df` and the
`pipeline._build_locus_meta` payload:

- `gdna_eff_len`
- `gdna_eff_len_per_bp = gdna_eff_len / gdna_span`

The per-bp ratio is a cheap runtime sanity check. For large single-interval
loci, it should be close to `1 + (E_G[ell] - 1) / gdna_span`; for nearby
multi-interval loci it also verifies that the union-merge logic is not
double-counting expanded start windows.

Do not add these fields to `quant.feather`, `nrna_quant.feather`, or
`gene_quant.feather`; they are locus diagnostics, not transcript/gene abundance
outputs. In `summary.json`, add a new `gdna_eff_len` block with aggregate
distribution statistics over loci, for example `min`, `median`, `p95`, and
`max` for both `gdna_eff_len` and `gdna_eff_len_per_bp`. That lets benchmark
reports flag pathological loci without loading the full locus table.

Optional profiling hook: if `emit_locus_stats=True`, add `gdna_eff_len` and
`gdna_log_eff_len` to the native `LocusProfile` dictionary. This is not needed
for correctness, but it is cheap and helpful for drill-down profiling.

## Validation Plan

Unit tests:

1. Geometry tests for `gdna_eff_len_for_loci` using degenerate
   fragment-length PMFs where exact integer counts are obvious.
2. Closed-form fast-path tests showing a single interval away from contig ends
   gives `span + E_G[ell] - 1` without per-length merging.
3. Edge and multi-interval tests showing clipped and unioned start windows
   match hand-counted answers.
4. `PriorTable.empty()` and `assemble_priors()` include a finite positive
   `gdna_eff_len` array aligned by `multi_locus_id`.
5. Native EM smoke test showing that increasing `locus_gdna_eff_lens[i]`
   decreases gDNA posterior odds when log-likelihoods and priors are fixed.
6. Python wrapper compatibility test showing omitted `gdna_eff_len` behaves like
   `np.ones(n_loci)` so existing wrapper-level tests keep passing, while the
   native binding itself requires `locus_gdna_eff_lens`.
7. Scoring regression showing gDNA log-likelihood no longer changes when the
   same unspliced fragment is anchored to transcripts with different genomic
   spans.
8. MM aggregation regression showing that removing `-log(e_h)` preserves the
   `logsumexp(hit_scores) - log(nh_gdna)` normalization. Include a fragment
   with hits in two `MultiLocus` components and assert that each component sees
   the same replicated `gdna_log_lik(u)` scalar, with `nh_gdna` reflecting the
   full MM group rather than a per-component subset.
9. Synthetic locus-15 regression, reconstructed from the
   `gdna_zero_ss_0.50_nrna_high` RCA geometry with deterministic FL PMF:
   `gdna_eff_len` is within a few percent of `locus_span + E_G[ell] - 1`, the
   gDNA posterior contribution per nRNA-like fragment drops by approximately
   `log(e_old / gdna_eff_len)`, and fractional-mode nRNA mass for the locus
   increases relative to the old denominator.

Integration tests:

1. Re-run the existing pipeline smoke and golden-output tests after intentional
   golden updates.
2. Re-run the synthetic 24 benchmark with this fix alone before any calibration
   retuning, and compare against the existing report.
3. Specifically inspect `gdna_zero_ss_0.50_nrna_high` in fractional mode and
   sample mode. This change is expected to reduce gDNA dominance in short-anchor
   large-locus cases and may materially improve unstranded nRNA recovery.
4. In no-nRNA scenarios, compare `sum_M E_EM[N_gdna(M)]` to
   `pi_pool * n_em_fragments`. If the EM gDNA total drifts materially from the
   calibrated pool expectation, treat that as evidence that the prior needs
   re-derivation or an explicit scale factor.
5. Stratify benchmark deltas by `gdna_span`, transcript span, and observed
   fragment length to verify that changes are largest where the old
   anchor-transcript denominator was most discrepant.
6. Inspect `gdna_eff_len_per_bp` distributions in `summary.json` and locus
   output to catch clipping or union bugs at runtime.

Acceptance criteria:

- gDNA effective length is visible as EM component metadata, not hidden in
  scorer log-likelihoods.
- gDNA scorer output is independent of anchor transcript span.
- MM gDNA scorer aggregation is still normalized by `nh_gdna` and does not
  double-count multiple alignments of the same fragment.
- For a single interval away from contig edges, the denominator equals
  `span + E_G[ell] - 1`.
- Multi-interval `MultiLocus` denominators use exact union-of-start support and
  do not double-count nearby intervals.
- `gdna_eff_len` and `gdna_eff_len_per_bp` diagnostics are emitted.
- The native EM binding requires `locus_gdna_eff_lens`; the Python wrapper is
   the only compatibility layer that supplies `ones` when omitted.
- Synthetic 24 fix-only results quantify both nRNA recovery and no-nRNA gDNA
  calibration drift before any prior tuning is attempted.
- Existing tests pass after recompilation and intentional golden updates.

## Cleanup And Commit Sequence

Recommended implementation sequence:

1. Add `gdna_eff_len_for_loci` plus unit tests, including the single-interval
   closed-form fast path.
2. Extend `PriorTable` with flat-array-only `gdna_eff_len`; update
   `assemble_priors`, including unconditional `ref_lengths_arr` construction.
3. Add the native EM required positional argument, add the Python wrapper
   `None -> ones` default, and thread real `prior_table.gdna_eff_len` through
   `pipeline` to `estimator` to native EM.
4. Remove scorer-side `e_h`, `gdna_flank`, and `t_span_`; add the replicated-MM
   scalar regression and the locus-15-style behavioral regression.
5. Add `loci.feather`, `summary.json`, and optional `locus_stats.feather`
   diagnostics.
6. Retire the deprecated `n_gdna_exon_intron` compatibility aggregate after the
   core effective-length tests are green. The canonical fields are
   `n_gdna_boundary_observed` and `n_gdna_exon_only`.

## Expected Effects And Risks

Expected effects:

- Clearer and more auditable gDNA/RNA symmetry in EM.
- Less anchor-transcript dependence for gDNA multimappers.
- A much larger gDNA denominator for short-anchor/large-locus cases compared
   with the current mean-flank per-hit correction, reducing gDNA likelihood by
   several to roughly ten nats per affected fragment.
- A plausible material recovery of unstranded nRNA mass in loci where the RCA
   showed nRNA-like fragments being siphoned into gDNA.
- Cleaner future path for mappability/capture-aware gDNA exposure.

Risks:

- Posterior shifts may expose existing calibration-prior strength issues,
  especially in unstranded nRNA-high cases.
- Golden outputs will likely change.
- If capture data has strong on-target enrichment, geometry-only overlap
  support remains incomplete; that should be addressed in a later weighted
  exposure phase rather than folded into this first correction.

The main point: this is the most likely fix for the per-fragment objective
asymmetry identified in the unstranded RCA. It may not fully solve nRNA/gDNA
identifiability, and calibration-prior tuning remains a separate lever, but the
fix should be evaluated as a potentially substantive deconvolution improvement
rather than only as code cleanup.