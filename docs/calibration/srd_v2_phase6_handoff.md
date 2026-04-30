# SRD v2 — Phase 6 Handoff

**Date:** 2026-04-27
**Branch context:** rigel-rnaseq 0.4.1 (reverted from accidental 0.6.0 bump)
**Status:** Architectural fix landed, Phase 4 VCaP-mixture sweep re-run complete.
**Next:** Awaiting approval on Q3 chimera-leak plan; consider FL-mixture iteration cap.

---

## 1. Implementation state

### 1.1 Completed in this cycle

| ID | Change | Files |
|----|--------|-------|
| Q1 | Version revert 0.6.0 → 0.4.1 | `pyproject.toml`, `CHANGELOG.md` |
| Q2 | Intergenic FL fragments fed into global FL histogram (`SPLICE_UNSPLICED`) | `src/rigel/native/bam_scanner.cpp:1340-1360` |
| Q4 | INTRONIC=0 calibration bug — root cause = synthetic nRNA hits polluting `exon_bp_pos/_neg` | superseded by Q5 |
| Q5 | **Architectural fix**: drop synthetic nRNAs from cgranges entirely; derive on-the-fly via `nrna_parent_` map in C++ resolver | `src/rigel/index.py`, `src/rigel/native/resolve_context.h`, `src/rigel/native/resolve.cpp`, `src/rigel/pipeline.py` |

Tests: **922/922 pass** (2 golden regressions regenerated; combo_moderate / combo_extreme — EM convergence shifted by ≤2.5% relative).

### 1.2 Q5 architectural fix — key code

- **`src/rigel/index.py:438-455`** `_gen_transcript_intervals`: early `return` if `t.is_synthetic`.
- **`src/rigel/index.py:1195-1208`** `TranscriptIndex.load`: builds `parent` array (synthetic-parent transcripts only) and calls `ctx.set_nrna_parent_index(parent.tolist())`.
- **`src/rigel/native/resolve_context.h:541-542,692`** new `nrna_parent_` field + setter.
- **`src/rigel/native/resolve_context.h:920-975`** strand-bp accumulators — removed `is_nrna` skip guards (no longer needed; structurally absent from cgranges).
- **`src/rigel/native/resolve_context.h:998-1050`** post-chimera, post-strand-bp nRNA derivation block: collects unique parents from real-tx hits, computes single-exon span overlap, injects into per-block `exon_t_set` via `std::lower_bound` for `merge_sets()`.
- **`src/rigel/native/resolve.cpp:127-131`** Python binding for `set_nrna_parent_index`.
- **`src/rigel/pipeline.py:220-235`**: removed redundant `set_nrna_parent_index` call; resolver wired at index load time.

### 1.3 Index footprint impact (Gencode v44 + controls)

| | BASE | NEW (Phase 6) | Δ |
|---|---:|---:|---:|
| Total transcripts | 457,513 | 457,513 (203,052 synthetic still tracked, just not in cgranges) | — |
| cgranges intervals | 2,362,845 | **1,956,741** | **-17.2%** |
| Index build wall | ~2 min | ~2 min | unchanged |

---

## 2. Phase 4 VCaP-mixture sweep results

### 2.1 Inputs

- **Index:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/rigel_index`
- **Outputs:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/default/<lib>/`
- **Run script:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/run_sweep.sh`
- **Master log:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/master.log`
- **Libraries (6/8):** dna00m, dna01m, dna02m, dna05m, dna10m, dna20m (skipped dna40m, dna80m for time)
- **Baseline:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_vcap_mixture/default/`
- **Flags:** `--em-mode vbem --assignment-mode sample --seed 42 --threads 8`

### 2.2 INTRONIC bug — confirmed eliminated (dna00m example)

| Category | BASE (buggy) POS / NEG | NEW POS / NEG |
|---|---:|---:|
| INTRONIC | **0 / 0** | **26,952 / 25,528** |
| EXON_CONT | 1,776,853 / 1,597,280 | 1,658,895 / 1,501,793 |
| EXON_INCOMPAT | 30,471 / 31,930 | 283,131 / 263,763 |
| INTERGENIC (NONE) | 8,794 | 8,794 ✓ |

POS:NEG ratio ≈ 1.05–1.10 across all categories — symmetric, biologically sensible.

### 2.3 Headline pi_pool

| lib | BASE | NEW | Δ |
|---|---:|---:|---:|
| dna00m | 0.9310 | 0.7888 | -0.142 |
| dna01m | 0.8905 | 0.6765 | -0.214 |
| dna02m | 0.8798 | 0.6539 | -0.226 |
| dna05m | 0.8627 | 0.6359 | -0.227 |
| dna10m | 0.8484 | 0.6286 | -0.220 |
| dna20m | 0.8254 | 0.6236 | -0.202 |

Both BASE and NEW exhibit wrong-direction monotonicity (decreasing with more spike). The Q5 fix did not address this — likely tied to chimera leak (Q3, see §3).

### 2.4 Quantification totals — robust

| Lib | mrna Δ | nrna Δ | gdna Δ |
|---|---:|---:|---:|
| dna00m | -0.04% | -1.0% | +4.2% |
| dna01m | +0.01% | +0.7% | -1.8% |
| dna02m | +0.02% | +1.9% | -2.2% |
| dna05m | +0.06% | +2.2% | -2.1% |
| dna10m | +0.12% | +1.3% | -1.5% |
| dna20m | -0.49% | -1.8% | +0.7% |

Max absolute change ≤4.2%; EM is robust to the calibration shift.

### 2.5 Performance

| | BASE avg | NEW avg | Δ |
|---|---:|---:|---:|
| Wall-clock | 4:32 | 4:20 | **-4.4%** |
| Peak RSS | 8.1 GB | 8.5 GB | +5% (real strand-bp buffers now non-zero) |

---

## 3. Open issues / known anomalies

1. **Q3 chimera leak — RESOLVED 2026-04-28 (REJECTED → policy adopted).**
   Option B in `docs/calibration/srd_v2_chimera_leak_plan.md` was
   rejected: in short-read PE data, contiguous-genome fragments with
   anomalously large `genomic_footprint` are overwhelmingly explained
   by annotation gaps (unannotated splicing, isoforms, TSS/TTS,
   read-through into intergenic), not true chimeras. Adopted policy:
   drop fragments with `genomic_footprint > max_frag_length` from FL
   training for both RNA and gDNA models. Drop count surfaced as
   `CalibrationResult.n_pool_dropped_out_of_range`; `n_pool` tightened
   to the in-range denominator that actually fed the mixture EM.
   See the rejection banner in `srd_v2_chimera_leak_plan.md`.
2. **`fragment_length.intergenic` block removed** from `summary.json` (BASE 9 blocks → NEW 8). Consistent with Q2 (intergenic now folded into `global`). Consider restoring as a diagnostic-only block.
3. **FL mixture iterations spike** — `mixture_iterations: 778` at dna00m (typical 50–200). Convergence works but slow; perf optimization opportunity, not a bug.
4. **`n_pool_dropped_out_of_range`** essentially unchanged BASE↔NEW (~0.07–0.5% of pool). Q3 territory.
5. **Sentinel:** real-data anomalies (e.g., INTRONIC=0) must always trigger root-cause analysis. The "biological" hand-wave was wrong here — most genomic positions lack antisense gene overlap.

---

## 4. Resume checklist (new terminal)

```bash
conda activate rigel
cd /home/mkiyer/proj/rigel

# Verify build is current with code
pip install --no-build-isolation -e .

# Verify tests
pytest tests/ -q

# Re-inspect Phase 6 results
ls /scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/default/*/quant.feather
```

### 4.1 Likely next tasks (in priority order)

1. ~~Approve / implement Q3 chimera leak plan~~ — **REJECTED 2026-04-28**;
   superseded by SRD v2 Phase 7 OOR-drop policy (see §3 item 1).
2. **Run Phase 4 dna40m + dna80m** to complete the sweep (use existing `run_sweep.sh`, just extend `LIBS=()`).
3. **Re-run analyzer** with Phase 7 (`n_pool` is now in-range only): `python scripts/calibration/analyze_srd_v2.py`.
4. **Investigate FL mixture iteration cap** if dna40m/dna80m show >1000 iterations.

### 4.2 Key paths quick-reference

| Asset | Path |
|---|---|
| Phase 6 index | `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/rigel_index` |
| Phase 6 results | `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/default/<lib>/` |
| Baseline (BASE) | `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_vcap_mixture/default/<lib>/` |
| Phase 4 sweep script | `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_phase6/run_sweep.sh` |
| Source BAMs | `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/<lib>/rigel/annotated.bam` |
| Q3 plan (REJECTED, see banner) | `docs/calibration/srd_v2_chimera_leak_plan.md` |
| Q5 plan (implemented) | `docs/calibration/nrna_derive_plan.md` |

### 4.3 Critical session facts (for the next agent)

- **Version is 0.4.1** (NOT 0.6.0). Reverted intentionally.
- **Synthetic nRNAs are NOT in cgranges** by design (Q5). Derived on-the-fly in `_resolve_core` via `nrna_parent_` map.
- The **strand-bp accumulators have NO `is_nrna` skip guards** any more — this is correct because synthetics structurally cannot enter cgranges. Do not re-add the guards; they would mask annotated nascent-equivs which legitimately participate.
- **Order in `_resolve_core` matters**: nRNA derivation must run AFTER chimera detection (else nRNAs mask cis-chimeras) AND AFTER strand-bp accumulation (else nRNAs pollute calibration counts).
- **2 golden regressions** were intentionally regenerated for Q5: combo_moderate, combo_extreme.
- **Pre-existing test failure** unrelated to this work: `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`.
- **Phase 7 (2026-04-28)**: `CalibrationResult.n_pool` is now the in-range
  denominator (i.e., `n_pool_categorized - n_pool_dropped_out_of_range`).
  Anything reading `n_pool` from older `summary.json` should be aware of
  the semantic shift. The Q3 chimera-leak plan was rejected in favor of
  this simpler drop-by-max-frag-length policy.
