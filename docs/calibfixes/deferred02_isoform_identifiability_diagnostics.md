# Deferred PR: Isoform Identifiability Diagnostics

## Status

Deferred from the v5 calibration-prior repair. Do not fix isoform collapse by adding RNA prior floors. A transcript dropping to zero under stronger gDNA/RNA priors may reflect an equivalence-class identifiability problem, not a calibration bug.

## Goal

Emit locus/transcript diagnostics that explain whether a zeroed or unstable transcript has unique fragment evidence. The diagnostic should identify rank-deficient or near-duplicate isoform columns in the EM compatibility matrix.

## Non-goals

- Do not add transcript-level pseudocount floors.
- Do not change EM assignments in this diagnostic PR.
- Do not collapse transcripts into genes internally.
- Do not use gene-level summaries to influence inference.

## Files likely to edit

| File | Purpose |
|---|---|
| `src/rigel/estimator.py` | Collect or pass EM unit compatibility data around locus EM. |
| `src/rigel/scored_fragments.py` or related EM data structure | Expose per-unit transcript compatibility needed for diagnostics. |
| `src/rigel/locus.py` | Attach diagnostic summaries to `MultiLocus` outputs. |
| `src/rigel/output` or CLI result writer | Emit optional diagnostic table when requested. |
| `tests/test_estimator.py` | Unit/integration tests for diagnostic calculation. |
| `tests/test_golden_output.py` | Update only if default output changes; prefer opt-in output. |

## Diagnostic matrix

For each MultiLocus, build a transcript compatibility matrix:

```text
M[u, t] = compatibility or likelihood contribution of EM unit u to transcript t
W[u, u] = unit weight, for example fragment count or total posterior/evidence weight
S = M.T @ W @ M
```

Use the same transcript columns that enter RNA EM components. gDNA should be excluded from transcript-rank diagnostics, but report whether gDNA is a strong competing component separately.

## Metrics to emit

Per locus:

```text
multi_locus_id
n_transcripts
n_units
matrix_rank
effective_rank
condition_number
smallest_singular_value
n_near_duplicate_pairs
has_gdna_candidate
prior_alpha_gdna_add
prior_alpha_rna_add
```

Per transcript:

```text
multi_locus_id
transcript_id
gene_id
estimated_count
unique_unit_count
unique_unit_mass
max_column_correlation
nearest_transcript_id
near_duplicate_group_id
isoform_ambiguous
zeroed_with_no_unique_evidence
```

## Implementation recipe

1. Add a pure helper module, for example `src/rigel/isoform_identifiability.py`.
2. Input should be plain arrays so tests can construct small matrices without a full pipeline run.
3. Compute singular values with `numpy.linalg.svd` on a dense matrix for small loci.
4. For large loci, use a sparse or randomized fallback only if needed; first implementation can skip huge loci with a flag and count them.
5. Effective rank can use entropy of normalized singular values:

```text
effective_rank = exp(-sum p_i log p_i)
```

6. Near-duplicate transcript pairs can use column cosine similarity above a strict diagnostic threshold such as `0.999`. This threshold is for reporting only, not inference.
7. Add an opt-in CLI flag, for example `--emit-isoform-diagnostics`, so default outputs are not churned.

## Tests

1. Full-rank toy matrix reports rank equal to transcript count and no ambiguity.
2. Duplicate columns report rank deficiency and a near-duplicate pair.
3. A transcript with a unique unit has `unique_unit_count > 0` and is not marked `zeroed_with_no_unique_evidence`.
4. A zero-estimated transcript with no unique evidence is marked ambiguous, not treated as a calibration-prior failure.
5. Large-locus skip/fallback path emits a clear flag and does not crash.

## Benchmark use

Use this diagnostic after PR 3 ESS sweeps. For every transcript that newly drops to zero under a stronger prior policy, report whether it had unique evidence. A stronger ESS policy should not be merged as default if well-conditioned transcripts collapse without a corresponding truth improvement.
