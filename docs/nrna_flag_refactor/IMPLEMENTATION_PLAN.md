# Refactor nRNA Flags: `is_nrna` + `is_synthetic`

## Problem Statement

The current codebase uses two boolean flags on transcripts to model nascent RNA:

- **`is_synthetic_nrna`**: True for RIGEL-generated synthetic single-exon transcripts that shadow multi-exon loci.
- **`is_nascent_equiv`**: True for annotated single-exon transcripts that happen to cover a merged nRNA span.

These names conflate two orthogonal concepts ("Is this transcript nRNA?" vs. "Did rigel synthesize it?") and create confusion:

1. **`is_synthetic_nrna`** is the *only* flag consulted for pool separation (mRNA vs nRNA) in scoring, routing, FL training exclusion, EM count attribution, and annotation. But annotated single-exon transcripts (`is_nascent_equiv=True`) are also indistinguishable from nascent RNA — they should participate in the same nRNA pool logic.

2. **`is_nascent_equiv`** is poorly named and documented. The MANUAL.md column description is wrong (says "synthetic nRNA span" when it means the opposite). It appears in only 4 code sites, 3 of which are output columns and 1 is the awkward `is_syn | is_equiv` union in `get_nrna_counts_df()` — a code smell that exists because the naming doesn't capture the actual concept.

3. The current approach creates an artificial asymmetry: annotated single-exon transcripts that cover nRNA spans get mRNA counts, while synthetic nRNA transcripts covering the same span get nRNA counts. The flag refactor makes this distinction explicit and principled.

## New Flag Semantics

Replace two flags with two orthogonal flags:

| New Flag | Meaning |
|----------|---------|
| `is_nrna` | This transcript is indistinguishable from nascent RNA (single-exon, overlaps multi-exon loci). Both annotated and synthetic transcripts can be nRNA. |
| `is_synthetic` | This transcript was created by rigel (not present in the annotation GTF). |

### Mapping from old to new

| Old state | New state |
|-----------|-----------|
| `is_synthetic_nrna=True` | `is_nrna=True, is_synthetic=True` |
| `is_nascent_equiv=True` | `is_nrna=True, is_synthetic=False` |
| Both False (regular transcript) | `is_nrna=False, is_synthetic=False` |

### Design rules

1. **All annotated single-exon transcripts** that cover a merged nRNA span are flagged `is_nrna=True`. They are "nascent RNA compatible."
2. **Synthetic nRNA transcripts** are created only when no annotated single-exon transcript covers the merged span (existing behavior, unchanged).
3. **`is_nrna`** is the single flag for: pool separation in scoring/C++, FL training exclusion, annotation pool assignment, and filtering nRNA entities for output.
4. **`is_synthetic`** determines: whether a transcript appears in `quant.feather` (synthetics are excluded), and whether its counts are attributed as mRNA or nRNA.

### Output behavior

| Output file | What appears | New flag columns |
|-------------|-------------|-----------------|
| `quant.feather` | All annotated transcripts (not synthetics). Includes annotated nRNA-compatible transcripts. | `is_nrna` column added |
| `nrna_quant.feather` | All nRNA entities (both annotated and synthetic). | `is_synthetic` column (already exists) |
| `gene_quant.feather` | All genes. nRNA counts are sum of `is_nrna` transcript counts per gene. | No new columns |
| `loci.feather` | All loci. nRNA counts are sum of `is_nrna` transcript counts per locus. | No new columns |

### Critical semantic change

With the new flags, `is_nrna` replaces the current `is_synthetic_nrna` *everywhere* that the code asks "is this transcript nascent RNA?". This means:

- **C++ pool separation** (`scoring.cpp`): Uses `is_nrna` array. Both annotated and synthetic nRNA transcripts compete in the nRNA pool, separately from mRNA. *This is the correct behavior* — an annotated single-exon transcript that covers a multi-exon locus should not compete with mRNA splice isoforms in the same likelihood pool.
- **FL training exclusion**: nRNA transcripts (both annotated and synthetic) are excluded from fragment-length model training.
- **Annotation pool code**: Fragments assigned to nRNA transcripts get `ZP=nRNA`.

---

## Scope

### Behavioral changes

1. **C++ pool separation now includes annotated nRNA transcripts.** Previously, only synthetic nRNA transcripts were separated into the nRNA pool. Now annotated single-exon transcripts flagged `is_nrna=True` also get nRNA pool separation. This should improve likelihood competition in loci with both spliced and unspliced isoforms.

2. **FL training excludes annotated nRNA transcripts.** Previously, unique fragments mapping to annotated nRNA-compatible transcripts were included in FL training. Now they are excluded (consistent with synthetic nRNA exclusion).

3. **Annotation pool assignment includes annotated nRNA.** Fragments assigned to annotated nRNA transcripts now get `ZP=nRNA` instead of `ZP=mRNA`.

4. **`quant.feather` gains `is_nrna` column, loses `is_nascent_equiv` column.**

5. **`nrna_quant.feather`**: annotated nRNA transcripts now report actual counts (unchanged for synthetics). Previously, annotated nascent-equiv transcripts appeared in `nrna_quant.feather` with count=0. Now they report actual nRNA counts since they participate in the nRNA pool.

### NOT changed

- The EM solver itself (C++ `em_solver.cpp`) does not reference these flags.
- Locus construction (`locus.py`) does not reference these flags.
- The synthetic nRNA creation algorithm in `build_nrna_entities()` is unchanged — synthetics are still only created when no annotated transcript covers the span.
- The `nrna_t_index` cross-reference from multi-exon transcripts to their nRNA entity is unchanged.

---

## File-by-File Changes

### Phase 1: Transcript dataclass and index

#### 1A. `src/rigel/transcript.py` (L42-43, L132-133)

Rename the two dataclass fields:

```python
# OLD:
is_synthetic_nrna: bool = False   # L42
is_nascent_equiv: bool = False    # L43

# NEW:
is_nrna: bool = False
is_synthetic: bool = False
```

Update `to_dict()` (L132-133):

```python
# OLD:
'is_synthetic_nrna': self.is_synthetic_nrna,
'is_nascent_equiv': self.is_nascent_equiv,

# NEW:
'is_nrna': self.is_nrna,
'is_synthetic': self.is_synthetic,
```

#### 1B. `src/rigel/index.py` — `build_nrna_entities()` (L145-300)

**Nascent-equiv detection (L269):**

```python
# OLD:
s_tx.is_nascent_equiv = True

# NEW:
s_tx.is_nrna = True
```

**Synthetic creation (L291):**

```python
# OLD:
is_synthetic_nrna=True,

# NEW:
is_nrna=True,
is_synthetic=True,
```

**Docstring (L153-155):** Update references from `is_nascent_equiv` / `is_synthetic_nrna` to `is_nrna` / `is_synthetic`.

#### 1C. `src/rigel/index.py` — Synthetic filtering (L774)

```python
# OLD:
t for t in transcripts if not t.is_synthetic_nrna

# NEW:
t for t in transcripts if not t.is_synthetic
```

#### 1D. `src/rigel/index.py` — Backward compatibility on load (L860-863)

```python
# OLD:
if "is_synthetic_nrna" not in self.t_df.columns:
    self.t_df["is_synthetic_nrna"] = False
if "is_nascent_equiv" not in self.t_df.columns:
    self.t_df["is_nascent_equiv"] = False

# NEW:
# Backward compat: rename old columns to new names
if "is_synthetic_nrna" in self.t_df.columns and "is_nrna" not in self.t_df.columns:
    # Old index: is_synthetic_nrna → is_nrna=True, is_synthetic=True
    #            is_nascent_equiv  → is_nrna=True, is_synthetic=False
    old_syn = self.t_df["is_synthetic_nrna"].values.astype(bool)
    old_equiv = self.t_df.get("is_nascent_equiv", pd.Series(False, index=self.t_df.index)).values.astype(bool)
    self.t_df["is_nrna"] = old_syn | old_equiv
    self.t_df["is_synthetic"] = old_syn
    self.t_df.drop(columns=["is_synthetic_nrna"], inplace=True)
    if "is_nascent_equiv" in self.t_df.columns:
        self.t_df.drop(columns=["is_nascent_equiv"], inplace=True)
if "is_nrna" not in self.t_df.columns:
    self.t_df["is_nrna"] = False
if "is_synthetic" not in self.t_df.columns:
    self.t_df["is_synthetic"] = False
```

---

### Phase 2: Pipeline and scoring (nRNA mask propagation)

#### 2A. `src/rigel/pipeline.py` — BAM scan setup (L242-243)

The nRNA mask passed to C++ resolver for FL training exclusion:

```python
# OLD:
nrna_arr = index.t_df["is_synthetic_nrna"].values.astype("uint8")

# NEW:
nrna_arr = index.t_df["is_nrna"].values.astype("uint8")
```

#### 2B. `src/rigel/pipeline.py` — Estimator constructor (L365)

```python
# OLD:
is_synthetic_nrna=index.t_df["is_synthetic_nrna"].values,

# NEW:
is_nrna=index.t_df["is_nrna"].values,
```

#### 2C. `src/rigel/pipeline.py` — `_populate_em_annotations()` (L442)

```python
# OLD:
is_syn = index.t_df["is_synthetic_nrna"].values

# NEW:
is_nrna = index.t_df["is_nrna"].values
```

And the pool assignment logic (around L462):

```python
# OLD:
is_nrna = is_syn[best_tid[valid_t]].astype(bool)
pool[valid_t] = np.where(is_nrna, POOL_CODE_NRNA, POOL_CODE_MRNA).astype(np.uint8)

# NEW (variable rename to avoid shadow):
nrna_mask = is_nrna[best_tid[valid_t]].astype(bool)
pool[valid_t] = np.where(nrna_mask, POOL_CODE_NRNA, POOL_CODE_MRNA).astype(np.uint8)
```

#### 2D. `src/rigel/scoring.py` — NativeFragmentScorer (L209-212)

```python
# OLD:
# Build is_synthetic_nrna array for pool separation
t_is_nrna_arr = np.zeros(index.num_transcripts, dtype=np.uint8)
if "is_synthetic_nrna" in index.t_df.columns:
    t_is_nrna_arr = index.t_df["is_synthetic_nrna"].values.astype(np.uint8)

# NEW:
# Build is_nrna array for pool separation
t_is_nrna_arr = np.zeros(index.num_transcripts, dtype=np.uint8)
if "is_nrna" in index.t_df.columns:
    t_is_nrna_arr = index.t_df["is_nrna"].values.astype(np.uint8)
```

#### 2E. `src/rigel/scan.py` — Deterministic-unambig pool routing (L196)

```python
# OLD:
is_syn = index.t_df["is_synthetic_nrna"].values
...
pool = POOL_CODE_NRNA if is_syn[tid] else POOL_CODE_MRNA

# NEW:
is_nrna_arr = index.t_df["is_nrna"].values
...
pool = POOL_CODE_NRNA if is_nrna_arr[tid] else POOL_CODE_MRNA
```

---

### Phase 3: Estimator output methods

#### 3A. `src/rigel/estimator.py` — Constructor (L114, L122-124)

```python
# OLD:
is_synthetic_nrna: np.ndarray | None = None,
...
self._synthetic_mask = (
    np.asarray(is_synthetic_nrna, dtype=bool)
    if is_synthetic_nrna is not None
    else np.zeros(num_transcripts, dtype=bool)
)

# NEW:
is_nrna: np.ndarray | None = None,
...
self._nrna_mask = (
    np.asarray(is_nrna, dtype=bool)
    if is_nrna is not None
    else np.zeros(num_transcripts, dtype=bool)
)
```

#### 3B. `src/rigel/estimator.py` — `get_transcript_counts_df()` (L355, L386, L393)

```python
# OLD (L355):
is_synthetic = index.t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in index.t_df.columns else ...
mrna = np.where(is_synthetic, 0.0, t_total)

# NEW:
is_synthetic = index.t_df["is_synthetic"].values if "is_synthetic" in index.t_df.columns else ...
is_nrna = index.t_df["is_nrna"].values if "is_nrna" in index.t_df.columns else ...
mrna = np.where(is_nrna, 0.0, t_total)
nrna = np.where(is_nrna, t_total, 0.0)
```

Output DataFrame (L386):

```python
# OLD:
"is_nascent_equiv": index.t_df["is_nascent_equiv"].values,

# NEW:
"is_nrna": is_nrna,
"nrna": nrna,
```

Filter (L393):

```python
# OLD:
df = df[~is_synthetic].reset_index(drop=True)

# NEW:
df = df[~is_synthetic].reset_index(drop=True)
# (unchanged — synthetics are still excluded from quant.feather)
```

**Key output change:** `quant.feather` gains `is_nrna` and `nrna` columns, loses `is_nascent_equiv`. Annotated nRNA transcripts now show `mrna=0, nrna=<count>` instead of `mrna=<count>`.

#### 3C. `src/rigel/estimator.py` — `get_gene_counts_df()` (L431)

```python
# OLD:
is_synthetic = index.t_df["is_synthetic_nrna"].values if ...
t_nrna = np.where(is_synthetic, t_counts_all, 0.0)

# NEW:
is_nrna = index.t_df["is_nrna"].values if "is_nrna" in index.t_df.columns else ...
is_synthetic = index.t_df["is_synthetic"].values if "is_synthetic" in index.t_df.columns else ...
t_nrna = np.where(is_nrna, t_counts_all, 0.0)
```

Also update `n_annotated` computation and `nrna` column in the gene output:

```python
# OLD:
annotated_mask = (~is_synthetic).astype(np.float64)

# NEW (unchanged logic but uses new column name):
annotated_mask = (~is_synthetic).astype(np.float64)
```

Gene-level nRNA now includes counts from both annotated and synthetic nRNA transcripts.

#### 3D. `src/rigel/estimator.py` — `get_nrna_counts_df()` (L502-503)

```python
# OLD:
is_syn = t_df["is_synthetic_nrna"].values if ...
is_equiv = t_df["is_nascent_equiv"].values if ...
mask = is_syn | is_equiv
...
counts = np.where(is_syn[mask], t_total[idx], 0.0)

# NEW:
is_nrna = t_df["is_nrna"].values if "is_nrna" in t_df.columns else ...
is_synthetic = t_df["is_synthetic"].values if "is_synthetic" in t_df.columns else ...
mask = is_nrna  # single flag — no union needed
...
counts = t_total[idx]  # all nRNA entities get their actual counts
```

The `is_synthetic` column in the output is unchanged:

```python
# OLD:
"is_synthetic": is_syn[mask],

# NEW:
"is_synthetic": is_synthetic[mask],
```

**Key semantic change:** Annotated nRNA transcripts now report actual counts in `nrna_quant.feather` (previously reported 0). This is correct because with the new pool separation, their fragments will be routed to the nRNA pool.

#### 3E. `src/rigel/estimator.py` — `get_loci_df()` (L590)

```python
# OLD:
is_syn = index.t_df["is_synthetic_nrna"].values if ...
syn_in_locus = in_locus & is_syn
nrna = float(t_total[syn_in_locus].sum())
n_syn = int(syn_in_locus.sum())
n_annot = int(in_locus.sum()) - n_syn

# NEW:
is_nrna = index.t_df["is_nrna"].values if "is_nrna" in index.t_df.columns else ...
is_synthetic = index.t_df["is_synthetic"].values if "is_synthetic" in index.t_df.columns else ...
nrna_in_locus = in_locus & is_nrna
nrna = float(t_total[nrna_in_locus].sum())
n_syn = int((in_locus & is_synthetic).sum())
n_annot = int(in_locus.sum()) - n_syn
```

Output column `n_nrna_entities` should count all nRNA entities (not just synthetics):

```python
# OLD:
"n_nrna_entities": n_syn,

# NEW:
"n_nrna_entities": int(nrna_in_locus.sum()),
```

---

### Phase 4: Test fixtures and golden outputs

#### 4A. `tests/conftest.py` (L292)

```python
# OLD:
"is_synthetic_nrna": np.zeros(num_transcripts, dtype=bool),

# NEW:
"is_nrna": np.zeros(num_transcripts, dtype=bool),
"is_synthetic": np.zeros(num_transcripts, dtype=bool),
```

#### 4B. `tests/test_estimator.py` (L68-69, L584)

```python
# OLD:
"is_nascent_equiv": [False] * num_transcripts,
"is_synthetic_nrna": [False] * num_transcripts,

# NEW:
"is_nrna": [False] * num_transcripts,
"is_synthetic": [False] * num_transcripts,
```

And the column check (L584):

```python
# OLD:
"is_nascent_equiv",

# NEW:
"is_nrna",
```

#### 4C. `tests/test_pipeline_routing.py` (L133)

```python
# OLD:
"is_synthetic_nrna": np.zeros(self.num_transcripts, dtype=bool),

# NEW:
"is_nrna": np.zeros(self.num_transcripts, dtype=bool),
"is_synthetic": np.zeros(self.num_transcripts, dtype=bool),
```

#### 4D. `tests/test_transcript_space_fl.py` (L631-633)

```python
# OLD:
nrna_mask = t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in t_df.columns else None

# NEW:
nrna_mask = t_df["is_nrna"].values if "is_nrna" in t_df.columns else None
```

#### 4E. Golden output regeneration

```bash
pytest tests/ --update-golden
```

Golden files in `tests/golden/` will be regenerated with the new column names.

---

### Phase 5: Scripts and documentation

#### 5A. `scripts/profiling/profiler.py` (L577)

```python
# OLD:
is_synthetic_nrna=index.t_df["is_synthetic_nrna"].values,

# NEW:
is_nrna=index.t_df["is_nrna"].values,
```

#### 5B. `scripts/profiling/profile_partition_capan1.py` (L134)

```python
# OLD:
is_synthetic_nrna=index.t_df["is_synthetic_nrna"].values,

# NEW:
is_nrna=index.t_df["is_nrna"].values,
```

#### 5C. `docs/MANUAL.md` (L242)

Update the `quant.feather` column table:

```markdown
# OLD:
| `is_nascent_equiv` | `1` if transcript is a synthetic nRNA span |

# NEW:
| `is_nrna` | `1` if transcript is compatible with nascent RNA |
```

Add `nrna` column to the table if not already present.

---

### Phase 6: C++ code

**No C++ source code changes required.** The C++ modules (`scoring.cpp`, `resolve_context.h`, `bam_scanner.cpp`) store the nRNA mask as `std::vector<uint8_t>` with generic names like `t_is_nrna_`. The Python side already passes the array by value — changing which column is used (Phase 2D, 2A) is sufficient.

---

## Verification Plan

### Step 1: Compile and unit tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/ -v
```

All 1029+ tests must pass.

### Step 2: Regenerate golden outputs

```bash
pytest tests/ --update-golden
```

Inspect diffs to verify column renames appear correctly:
- `is_nascent_equiv` → `is_nrna` in transcript output
- `is_synthetic_nrna` should not appear in any output

### Step 3: Lint

```bash
ruff check src/ tests/
ruff format src/ tests/
```

### Step 4: Index backward compatibility

Load an existing index (built before the refactor) and verify:
- Old `is_synthetic_nrna` and `is_nascent_equiv` columns are auto-mapped to `is_nrna` and `is_synthetic`
- No errors on load

### Step 5: Full-scale benchmark

```bash
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --force
python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml -o results/benchmark_nrna_refactor
```

Compare to pre-refactor baseline. Expected changes:
- Pool-level: nRNA counts should increase (annotated nRNA transcripts now contribute)
- mRNA counts should decrease correspondingly
- Transcript-level correlation should be similar or better (cleaner pool separation)
- The nRNA siphon effect in `nrna_rand` conditions may change magnitude

### Step 6: Annotated BAM spot-check

```bash
rigel quant --bam <test.bam> --index <index> -o /tmp/test --annotated-bam /tmp/test/annotated.bam
# Verify annotated nRNA transcripts get ZP=nRNA
```

---

## Risk Assessment

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| Golden test diffs from column renames | Certain | `--update-golden`; review diffs |
| Benchmark regression from pool separation change | Moderate | The change is directionally correct; compare before/after |
| Old indexes fail to load | Low | Backward compat code handles remapping |
| Scripts break on renamed columns | Low | grep search found all sites |
| nRNA siphon increases | Possible | More transcripts in nRNA pool = more absorption. Monitor in benchmarks. |

---

## Summary of All File Changes

| File | Lines affected | Nature of change |
|------|---------------|-----------------|
| `src/rigel/transcript.py` | L42-43, L132-133 | Rename fields |
| `src/rigel/index.py` | L153-155, L269, L291, L774, L860-863 | Flag assignment + backward compat |
| `src/rigel/pipeline.py` | L242, L365, L442, L462 | Column references |
| `src/rigel/scoring.py` | L209-212 | Column reference |
| `src/rigel/scan.py` | L196 | Column reference |
| `src/rigel/estimator.py` | L114, L122-124, L355, L386, L393, L431, L502-503, L590 | Constructor + 4 output methods |
| `tests/conftest.py` | L292 | Fixture |
| `tests/test_estimator.py` | L68-69, L584 | Mock index + column check |
| `tests/test_pipeline_routing.py` | L133 | Mock index |
| `tests/test_transcript_space_fl.py` | L631-633 | FL test |
| `scripts/profiling/profiler.py` | L577 | Column reference |
| `scripts/profiling/profile_partition_capan1.py` | L134 | Column reference |
| `docs/MANUAL.md` | L242 | Documentation |
| `tests/golden/*` | All golden files | Regenerated |

**Total: 13 source files + golden outputs. Zero C++ changes.**
