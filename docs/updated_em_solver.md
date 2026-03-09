

This is an excellent insight. Let me study the current implementation to validate the proposal and plan the rewrite. 

Starting: *Study current EM solver code* (1/3)

Study current EM solver architecture

Searched for files matching `**/em_solver.cpp`, 1 match

Searched for files matching `**/estimator.py`, 1 match

Searched for files matching `**/locus.py`, 1 match

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1-1), lines 1 to 100

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/estimator.py#1-1), lines 1 to 100

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/locus.py#1-1), lines 1 to 100

Searched for regex `^(void|bool|double|int|py::object).*map_em_step|run_squarem|linked_map_em_step|linked_run_squarem|extract_locus_sub_problem|assign_posteriors|batch_locus_em|run_locus_em_native` (`**/em_solver.cpp`), 20 results

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#400-400), lines 400 to 814

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#615-615), lines 615 to 700

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#815-815), lines 815 to 920

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1081-1081), lines 1081 to 1200

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1325-1325), lines 1325 to 1510

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1608-1608), lines 1608 to 1750

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1750-1750), lines 1750 to 1930

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1930-1930), lines 1930 to 2075

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#2075-2075), lines 2075 to 2250

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#2250-2250), lines 2250 to 2550

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/estimator.py#100-100), lines 100 to 1519

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/estimator.py#1-1), lines 1 to 450

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/locus.py#100-100), lines 100 to 600

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/locus.py#600-600), lines 600 to 859

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#2700-2700), lines 2700 to 2800

Created memory file [](file:///Users/mkiyer/Library/Application%20Support/Code/User/workspaceStorage/9ccf7fd595e4847da09af2f90d434095/GitHub.copilot-chat/memory-tool/memories/N2QxOWJiOWMtZGZmYi00Y2M0LThiMjgtNzdmNmRiZWZmZTc3/rigel_em_solver_architecture.md)

Created memory file [](file:///Users/mkiyer/Library/Application%20Support/Code/User/workspaceStorage/9ccf7fd595e4847da09af2f90d434095/GitHub.copilot-chat/memory-tool/memories/N2QxOWJiOWMtZGZmYi00Y2M0LThiMjgtNzdmNmRiZWZmZTc3/rigel_em_solver_detailed_reference.md)

Starting: *Analyze mathematical correctness* (2/3)

Let me read the critical sections directly to verify the details before forming my analysis. 

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#520-520), lines 520 to 580

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#615-615), lines 615 to 820

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#820-820), lines 820 to 1080

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1081-1081), lines 1081 to 1340

Now let me read `run_locus_em_native` and the Python-side `LocusEMInput` and the prior computation: 

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1340-1340), lines 1340 to 1520

Now let me read the Python-side LocusEMInput and the nrna_frac prior computation: 

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/estimator.py#195-195), lines 195 to 330

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/estimator.py#530-530), lines 530 to 700

Updated todo list

Now let me read a few more critical sections — the assign_posteriors, extract_locus_sub_problem, and how results flow back: 

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1930-1930), lines 1930 to 2080

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1608-1608), lines 1608 to 1930

Now let me check the LocusSubProblem struct and how `batch_locus_em` calls the solver:  

Searched for text `struct LocusSubProblem` (`**/src/rigel/native/em_solver.cpp`), 1 result

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#1567-1567), lines 1567 to 1610

Now let me read how batch_locus_em calls the solver to see the process_locus lambda:   

Read [](file:///Users/mkiyer/proj/rigel/src/rigel/native/em_solver.cpp#2300-2300), lines 2300 to 2550

Starting: *Write implementation plan* (3/3)

I've now thoroughly studied the codebase. Your insight is **correct and I fully agree**. Here's my analysis and implementation plan.

---

## Analysis: Why the Linked Solver Is Now Unnecessary

You're right on the core insight. The linked solver (`linked_map_em_step` + `linked_run_squarem`) exists to solve a problem that no longer exists:

**Old world (2T+1):** 5 transcripts → 5 identical nRNA shadows with identical effective lengths, identical scores, identical priors. The EM's posterior assignment is perfectly symmetric — it can never learn which shadow "owns" the intronic reads. The linked solver breaks this by maintaining a separate parameter space (θ_t, η_t) and using the per-transcript Beta prior to force apportionment.

**New world (T+N+1):** 5 transcripts → 1 shared nRNA. The intronic reads unambiguously map to that single nRNA component. The EM assigns counts directly. No symmetry to break.

The current linked solver carries significant complexity:
- **State space mismatch**: SQUAREM operates on θ_t (T+1 dims) but the E-step operates on T+N+1 components. The decomposition `θ[n_t+n] = Σ(θ_t[t] × η_t)` and the inverse (collapse warm start → θ_t space) is error-prone.
- **η_t not extrapolated**: SQUAREM only accelerates θ_t, while η_t is updated in each M-step but not extrapolated. This limits acceleration because the nascent component can only change as fast as the M-step moves η.
- **Complex pruning**: Must check if ALL transcripts sharing an nRNA are pruned before zeroing the nRNA prior.
- **Apportionment arithmetic**: The M-step's mass-share computation `W_t = θ_t × η_t / α_n` requires careful epsilon handling.

With the flat hierarchical solver, SQUAREM operates directly on the full T+N+1 state vector. Every component — mRNA, nRNA, and gDNA — gets extrapolated. The hierarchical M-step enforces the η constraint naturally within the EM fixed-point iteration.

**Mathematical verification that post-hoc η is exact:**

After normalization: `θ[n_t+n] = η_n × C_group / Z` and `Σ_{t∈n} θ[t] = (1-η_n) × C_group / Z`, so `θ[n_t+n] / (θ[n_t+n] + Σ θ[t∈n]) = η_n`. The global normalization cancels out — we can recover η losslessly from the flat θ.

**One important refinement:** The Beta prior (α, β) for η should be **per-nRNA**, not per-transcript. This is biologically correct — the nascent fraction is a property of the transcription unit (defined by the nRNA span), not individual splice isoforms. It also simplifies the prior computation hierarchy.

---

## Implementation Plan: Hierarchical EM Solver

### Phase A: C++ em_solver.cpp

#### A1. New `hierarchical_map_em_step()` replacing both `map_em_step()` and `linked_map_em_step()`

```cpp
static void hierarchical_map_em_step(
    const double* theta,          // [n_components]: flat T+N+1 weights
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,    // [n_components]
    const double* unambig_totals, // [n_components]
    const double* prior,          // [n_components]
    double*       em_totals,      // [n_components]: zeroed then accumulated
    double*       theta_new,      // [n_components]: output (normalized)
    int           n_components,
    // Hierarchical nRNA parameters (all NULL/0 for pure unlinked mode)
    int           n_transcripts,  // 0 if no nRNA
    int           n_nrna,         // 0 if no nRNA
    const int32_t* t_to_nrna,    // [n_t]: local transcript → local nRNA
    const int32_t* nrna_to_t_off, // [n_nrna + 1]: CSR offsets
    const int32_t* nrna_to_t_idx, // CSR flat
    const double* nrna_frac_alpha, // [n_nrna]: per-nRNA Beta alpha
    const double* nrna_frac_beta,  // [n_nrna]: per-nRNA Beta beta
    int           estep_threads,
    rigel::EStepThreadPool* pool);
```

**E-step:** Identical to current `map_em_step()` — compute log_weights from θ, accumulate em_totals.

**M-step branches:**
- If `n_nrna == 0`: standard MAP-EM (`θ_new[i] = em_totals[i] + unambig[i] + prior[i]`, normalize).
- If `n_nrna > 0`: hierarchical update per the user's algorithm:
  1. For each nRNA group: compute C_nrna, C_mrna, C_group (including prior and unambig)
  2. Update η_n via Beta MAP
  3. Set θ_new[n_t+n] = η_n × C_group
  4. Set θ_new[t] = (1-η_n) × φ_t × C_group
  5. Handle gDNA independently
  6. Normalize

#### A2. Modify `run_squarem()`

Add optional hierarchical parameters. The SQUAREM loop stays structurally identical — two EM steps, extrapolation, stabilization step — but calls `hierarchical_map_em_step()` instead of `map_em_step()`. State vector size = T+N+1 (same as current unlinked). Remove VBEM path (it was never used with linked mode; if still needed for unlinked, we keep it but it's separate).

Pruning: Operates per-component on the flat θ. When an nRNA component is pruned, its prior is zeroed and mass redistributes through the group in the post-prune EM step. No special "all-transcripts-pruned" check needed — the standard per-component evidence check handles it.

Remove `nrna_frac` from `EMResult`.

#### A3. Delete `linked_map_em_step()` and `linked_run_squarem()` entirely

~470 lines of complex code deleted.

#### A4. Modify `LocusSubProblem`

```cpp
// Change:
std::vector<double> nrna_frac_alpha;  // [n_t] → [n_nrna]
std::vector<double> nrna_frac_beta;   // [n_t] → [n_nrna]
```

In `extract_locus_sub_problem()`: gather per-nRNA α/β from a global per-nRNA array (instead of per-transcript).

#### A5. Simplify `run_locus_em_native()`

- Remove linked/unlinked branching — always call `run_squarem()` with hierarchical params
- Remove "collapse warm start to theta_t space" logic — warm start is already T+N+1
- Change signature: per-nRNA α/β instead of per-transcript
- Return 3-tuple (θ, α, em_totals) instead of 4-tuple (drop nrna_frac)

#### A6. Simplify `batch_locus_em()`

- Remove linked/unlinked branching in `process_locus`
- Remove "collapse warm start" and "nrna_frac_init" setup
- LocusSubProblem already carries per-nRNA α/β
- Pass nRNA arrays from sub-problem to `run_squarem()`

#### A7. Update nanobind bindings

- `run_locus_em_native`: per-nRNA α/β, return 3-tuple
- `batch_locus_em`: per-nRNA α/β arrays (global)

### Phase B: Python estimator.py

#### B1. `LocusEMInput`

```python
# Change:
nrna_frac_alpha: np.ndarray  # float64[n_transcripts] → float64[n_nrna]
nrna_frac_beta: np.ndarray   # float64[n_transcripts] → float64[n_nrna]
```

#### B2. `AbundanceEstimator`

- `nrna_frac_alpha`, `nrna_frac_beta`: size N (per-nRNA) instead of T (per-transcript)
- Remove `_t_to_nrna` as a constructor parameter for nrna_frac purposes (still needed for output fan-out)
- η computed post-hoc from θ in output methods

#### B3. `compute_hybrid_nrna_frac_priors()`

Rewrite to produce per-nRNA priors. The hierarchy simplifies:
1. **Per-nRNA** ← shrink toward locus-strand
2. **Locus-strand** ← shrink toward global empirical

The TSS-group level can be kept by grouping nRNAs that share the same start position (nRNA start = transcription start). The intronic evidence is already per-nRNA, so the base-level computation is natural.

#### B4. `run_batch_locus_em()`

- Pass per-nRNA α/β arrays to C++ `batch_locus_em()`
- Handle 3-value return (no more nrna_frac)

#### B5. `get_counts_df()` / `get_gene_counts_df()`

- Compute η from θ: `η_n = nrna_em_counts[n] / (nrna_em_counts[n] + Σ mrna_counts[t ∈ n])`
- Fan-out logic stays similar (nRNA counts to per-transcript via `t_to_nrna`)

### Phase C: Python `locus.py`

#### C1. `build_locus_em_data()`

- `nrna_frac_alpha/beta`: slice from per-nRNA arrays using `local_to_global_nrna`

### Phase D: Pipeline integration

- `pipeline.py`, `profiler.py`: Update `compute_hybrid_nrna_frac_priors()` call to produce per-nRNA priors
- Remove per-transcript nrna_frac plumbing

### Phase E: Tests

- `test_em_impl.py`: Update `_run_em` and `_run_linked_em` helpers for unified solver
- `test_bias.py`: Update `run_locus_em_native` calls (3-tuple return, per-nRNA α/β)
- `test_estimator.py`: Update mock for per-nRNA priors
- `test_gdna.py`: Minimal changes
- Regenerate golden output

---

### Summary of deletions vs additions

| Metric | Before | After |
|--------|--------|-------|
| `linked_map_em_step()` | ~200 lines | **deleted** |
| `linked_run_squarem()` | ~270 lines | **deleted** |
| `hierarchical_map_em_step()` | (new) | ~80 lines |
| SQUAREM state dimensions | T+1 (linked) or T+N+1 (unlinked) | T+N+1 (unified) |
| EM functions calling conventions | 2 branches (linked/unlinked) in run_locus_em_native, batch_locus_em | 1 unified path |
| nrna_frac tracking | per-transcript array in EMResult | not returned; computed post-hoc from θ |
| nrna_frac prior space | per-transcript (size T) | per-nRNA (size N) |

Net: **~400 lines deleted**, ~80 added in C++. Significant simplification of control flow in `run_locus_em_native()`, `batch_locus_em()`, and the Python calling code.

---

Shall I proceed with implementation?