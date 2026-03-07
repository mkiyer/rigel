# rigel Second-Order Cleanup Plan

Date: 2026-03-06

Follows the completed [cleanup_plan.md](cleanup_plan.md) (Phases 1–4 done,
748 tests passing, 21 golden scenarios bit-exact).

This plan addresses findings from a second external code review plus
internal analysis.  The goal is to close the remaining maintainability
gaps so that adding or renaming a parameter is a one-file change rather
than a synchronized three-file edit.

---

## Principles

1. Every fix must pass the full 748-test suite and 21 golden scenarios.
2. Prefer eliminating code over adding abstractions.
3. A new parameter should require **one** addition (the dataclass field)
   plus its argparse line — everything else should be derived.
4. No behavioral or numerical changes.

---

## Phase 5 — Collapse CLI Parameter Projection (HIGH)

**Problem:** config.py advertises itself as the "single source of truth"
(line 3), but adding a parameter still requires synchronized edits in
three manually-maintained projection layers:

| Layer | Location | What it does |
|-------|----------|-------------|
| 1. `_build_quant_defaults()` | cli.py:415–527 | Instantiates default dataclasses, manually copies ~30 fields into a flat dict (with boolean inversions and type coercions) |
| 2. `_resolve_quant_args()` | cli.py:530–580 | Merges CLI → YAML → defaults dict onto `argparse.Namespace`, field by field |
| 3. `_build_pipeline_config()` | cli.py:125–190 | Manually copies ~30 fields back from Namespace into frozen dataclasses (with inverse boolean coercions) |

The keep\_duplicates ↔ skip\_duplicates double inversion is the canonical
example: `_build_quant_defaults` does `not scan.skip_duplicates` (line 445),
then `_build_pipeline_config` does `not args.keep_duplicates` (line 166).
The mapping is correct today, but the indirection invites future drift.

Additionally, `_write_run_config()` (cli.py:192–235) manually lists
every parameter a **fourth** time for config.json serialization.

**Fix — declarative parameter registry:**

Create a module-level list of `_ParamSpec` records that is the single
description of **how CLI flags map to config fields**.  Each record carries:

```python
@dataclass(frozen=True)
class _ParamSpec:
    cli_dest: str              # argparse dest, e.g. "em_prior_alpha"
    config_path: str           # dotted path, e.g. "em.prior_alpha"
    transform: str = "direct"  # "direct" | "invert_bool" | "log_penalty" | ...
```

Then:

1. **`_build_quant_defaults()`** becomes a loop: instantiate the default
   `PipelineConfig()`, walk `_PARAM_SPECS`, extract each field via
   `config_path`, apply the inverse transform, build the dict.  ~15 lines.

2. **`_build_pipeline_config()`** becomes a loop: walk `_PARAM_SPECS`,
   read `getattr(args, spec.cli_dest)`, apply the forward transform,
   `setattr` onto a mutable builder, freeze.  ~15 lines.

3. **`_write_run_config()`** becomes: serialize `args.__dict__` filtered
   to the known `cli_dest` keys.  ~5 lines.

4. **`_resolve_quant_args()`** stays the same (it already loops over the
   defaults dict and is transform-agnostic).

Special cases that need `transform` variants:
- `keep_duplicates` ↔ `scan.skip_duplicates` → `"invert_bool"`
- `overhang_alpha` → `scoring.overhang_log_penalty` → `"log_penalty"`
- `mismatch_alpha` → `scoring.mismatch_log_penalty` → `"log_penalty"`
- `gdna_splice_penalty_unannot` → `scoring.gdna_splice_penalties` → custom
- `sj_strand_tag` → list ↔ str normalization
- `threads` → fans out to `em.n_threads` + `scan.n_scan_threads`

The argparse `add_argument` calls stay as-is — they define help text and
types, which aren't derivable from the dataclass.  But the mapping between
arg dest and config field is no longer hand-maintained in four places.

**Validation:** Add a startup `assert` (or test) that every `_ParamSpec.cli_dest`
exists in argparse and every `config_path` resolves on `PipelineConfig`.

**Files:** cli.py (rewrite helpers), new `_param_specs` list (in cli.py or config.py)
**Risk:** Medium — pure plumbing, but touches the entire CLI→config pipeline.
Must verify YAML override, CLI precedence, and config.json output.
**Tests:** Existing test_cli.py covers argument parsing and precedence.
Add one test asserting round-trip: `_build_pipeline_config(_resolve(...))` yields
a config equal to `PipelineConfig()` when no CLI or YAML overrides are given.

---

## Phase 6 — Version Metadata Deduplication (MEDIUM)

**Problem:** The version string "0.1.0" is defined in three places:

| Location | Purpose |
|----------|---------|
| `pyproject.toml:7` | Build system / PyPI |
| `src/rigel/__init__.py:1` | `__version__` at runtime |
| `cli.py:292` (inside `_write_quant_outputs`) | Hardcoded in `summary.json` |

A version bump requires updating all three.  The summary.json copy is
particularly dangerous because it produces incorrect provenance silently.

**Fix:**

1. **`__init__.py` already reads correctly** — `get_version()` in cli.py
   imports `__version__` from `__init__.py`.  No change needed there.

2. **`pyproject.toml` → `__init__.py`:** Use `importlib.metadata` to
   derive `__version__` from the installed package metadata at import time:
   ```python
   from importlib.metadata import version, PackageNotFoundError
   try:
       __version__ = version("rigel")
   except PackageNotFoundError:
       __version__ = "0.0.0+unknown"
   ```
   This eliminates the `__init__.py` hardcoded string.  `pyproject.toml`
   becomes the sole authority.

3. **`cli.py:292`:** Replace the hardcoded `"0.1.0"` with:
   ```python
   from . import __version__
   summary = {
       "rigel_version": __version__,
       ...
   }
   ```

**Files:** `__init__.py`, `cli.py`
**Risk:** Low.  `importlib.metadata` is stdlib since Python 3.8.
**Tests:** Add a test that `summary.json["rigel_version"]` matches
`rigel.__version__`.

---

## Phase 7 — Python/C++ Constant Parity (MEDIUM)

**Problem:** `scoring.cpp` re-defines three constants that already exist
on the Python side (and in `constants.h` for `STRAND_NEG`):

| Constant | Python source | C++ source | Risk |
|----------|--------------|-----------|------|
| `LOG_HALF` | `math.log(0.5)` in scoring.py:34 | `-0.6931471805599453` in scoring.cpp:32 | Low (Python computes dynamically) |
| `TAIL_DECAY_LP` | `math.log(0.99)` in frag\_length\_model.py:33 | `-0.01005033585350145` in scoring.cpp:33 | Low (Python computes dynamically) |
| `STRAND_NEG` | `int(Strand.NEG)=2` in types.py:82 | `2` in scoring.cpp:34 | Medium (already in constants.h:56 but scoring.cpp doesn't use it) |

The resolve module already exports constants cleanly from `constants.h`
via `resolve.cpp:120–131`.  The scoring module doesn't follow this pattern.

**Fix:**

1. **scoring.cpp:** Replace the three local `constexpr` definitions with
   includes from `constants.h` (for `STRAND_NEG`) and a new `scoring_constants.h`
   (for `LOG_HALF` and `TAIL_DECAY_LP`).  Or just add `LOG_HALF` and
   `TAIL_DECAY_LP` to `constants.h` since they're shared.

2. **Export from C++:** Add `m.attr("LOG_HALF") = ...` and
   `m.attr("TAIL_DECAY_LP") = ...` to the scoring module's nanobind init.

3. **Python side:** `scoring.py` and `frag_length_model.py` can optionally
   import from the C++ module with a Python-computed fallback, or keep
   computing dynamically (since `math.log` is exact at IEEE 754 level).

4. **Parity test:** Add a test that asserts `_scoring_impl.LOG_HALF == math.log(0.5)`
   and `_scoring_impl.TAIL_DECAY_LP == math.log(0.99)` to catch any future
   divergence.

The simplest effective fix: keep Python computing dynamically, make
scoring.cpp `#include "constants.h"` for `STRAND_NEG`, and add a parity
test for the floating-point constants.

**Files:** `scoring.cpp`, `constants.h`, new test in `test_types.py`
**Risk:** Low — no behavioral change, just source reorganization + assertion.

---

## Phase 8 — FragmentBuffer Encapsulation (MEDIUM)

**Problem:** `pipeline.py` directly mutates `FragmentBuffer` private
fields in two places:

| Location | Mutation |
|----------|----------|
| pipeline.py:371–378 | `buffer._chunks.append(chunk)`, `buffer._total_size += size`, `buffer._memory_bytes += chunk.memory_bytes`, `buffer._spill_oldest()` |
| pipeline.py:502–503 | `buffer._chunks.clear()`, `buffer._memory_bytes = 0` |

The first block (scan\_and\_buffer) constructs chunks externally and injects
them.  This exists because the native scanner produces raw dicts that
pipeline.py converts to `_FinalizedChunk` objects — the buffer's own
`_finalize_native()` does the same conversion but from its internal
`FragmentAccumulator`.

The second block (quant\_from\_buffer) is a post-scan cleanup to release
memory before the EM phase.

**Fix:**

1. **Add `FragmentBuffer.inject_chunk(chunk)`:** Public method that
   appends a `_FinalizedChunk`, updates `_total_size` and `_memory_bytes`,
   and triggers spill-if-needed.  This is exactly the logic currently in
   `_finalize_native()` lines 457–468, just without the `_native_acc`
   finalization step.  Extract to a shared `_accept_chunk()` internal,
   expose as `inject_chunk()`.

2. **Add `FragmentBuffer.release()`:** Public method that clears
   `_chunks`, resets `_memory_bytes` to 0, and calls `cleanup()` for
   spilled files.  Replaces the two-line hack in `quant_from_buffer`.

3. **pipeline.py:** Replace the direct mutations with:
   ```python
   buffer.inject_chunk(chunk)      # replaces lines 371-378
   ...
   buffer.release()                # replaces lines 502-503
   ```

**Files:** `buffer.py`, `pipeline.py`
**Risk:** Low — extracting existing logic into methods with identical semantics.
**Tests:** Existing buffer tests + golden scenarios.

---

## Phase 9 — Unify Fragment-Length Defaults (MEDIUM)

**Problem:** Two competing defaults for maximum fragment length:

| Location | Default | Used by |
|----------|---------|---------|
| `frag_length_model.py:42` `DEFAULT_MAX_FRAG_SIZE` | 2000 | `FragmentLengthModels()` constructor |
| `config.py:194` `BamScanConfig.max_frag_length` | 1000 | Production pipeline via `scan_and_buffer()` |

The production pipeline always passes `max_size=scan.max_frag_length` (1000),
but standalone construction `FragmentLengthModels()` silently uses 2000.
Tests like `test_gdna.py:34` instantiate `FragmentLengthModels()` without
arguments, getting the 2000 default.  This means tests run with a different
statistical envelope than production — fragment lengths 1001–2000 get
individual bins in tests but hit the overflow tail in production.

**Fix:**

1. **Change `DEFAULT_MAX_FRAG_SIZE` to 1000** in `frag_length_model.py`
   to match the production config default.

2. **Or delete `DEFAULT_MAX_FRAG_SIZE`** and change the
   `FragmentLengthModels.__init__` signature to require `max_size`
   (no default), forcing callers to be explicit.

Option 1 is simpler and less disruptive.  The 2000 value was a
conservative over-estimate; production already proved 1000 is sufficient
(RNA-seq fragments are typically 100–500 bp).

3. **Audit test helpers** that instantiate `FragmentLengthModels()` to
   verify they don't depend on the 1001–2000 range.

**Files:** `frag_length_model.py`, possibly test helpers
**Risk:** Low — only affects overflow bin boundary, no numerical change
for fragments ≤1000 (which is all of production and nearly all tests).

---

## Phase 10 — Documentation Truthfulness (LOW)

**Problem:** Several documentation mismatches:

1. **BamScanConfig docstring** (config.py:172) says `include_multimap`
   defaults to `False`, but the actual field default is `True` (line 193).
   `parameters.md` correctly says "yes".

2. **README.md** is effectively empty (one line: "rna-seq transcript
   quantification tool").  Users have no quick-start guide, forcing them
   into scattered docs.

**Fix:**

1. **Fix the docstring:** Change "default False" → "default True" in the
   `BamScanConfig` docstring.

2. **Write a minimal README** with: one-paragraph description, install
   instructions, hello-world quant command, pointer to `docs/parameters.md`.
   This is a documentation task, not a code task.

**Files:** `config.py`, `README.md`
**Risk:** None.

---

## Execution Order

| Phase | Severity | Deps | Summary |
|-------|----------|------|---------|
| **5** | HIGH | — | Declarative CLI↔config mapping (eliminates 4× manual sync) |
| **6** | MEDIUM | — | Version deduplication (importlib.metadata) |
| **7** | MEDIUM | — | C++ constant parity (include constants.h, add parity test) |
| **8** | MEDIUM | — | Buffer encapsulation (inject\_chunk + release methods) |
| **9** | MEDIUM | — | Fragment-length default unification (2000→1000) |
| **10** | LOW | 6 | Docstring fix + README |

Phases 6–10 are independent and can proceed in any order.
Phase 5 is the highest-value change and should go first.

---

## What We're NOT Doing (and Why)

| Suggestion | Why skip |
|---|---|
| Auto-generate argparse from dataclass fields | Help text, type coercions, and grouping are CLI-specific; they don't map cleanly from dataclass metadata. The `_ParamSpec` registry is the right middle ground. |
| Move `_ParamSpec` to a separate module | It's ~50 lines of data; keeping it in cli.py avoids another import and another file to maintain. |
| Export ALL C++ constants to Python | Only the 3 scoring constants are duplicated. The resolve module already exports strand/splice constants correctly. |
| Fragment-length: require max\_size everywhere | Too disruptive for test ergonomics; changing the default to 1000 achieves parity with zero test rewrites. |
| Full README rewrite | Out of scope for a code cleanup plan; a minimal quick-start section is enough. |
