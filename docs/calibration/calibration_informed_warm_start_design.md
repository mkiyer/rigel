# Calibration-Informed Per-Transcript Warm Start — Design & Implementation

> ## ⏸️ STATUS: DEFERRED (2026-07-10) — not in production
>
> Built end-to-end and validated against ground truth on branch **`calib-warmstart-v1`** (off v0.7.0),
> behind the **default-off** `RIGEL_EM_WARMSTART=calib` flag. **Deferred, not merged.** It added
> substantial calibration→EM plumbing and, on the ground-truth A/B, **helped the flagship but net-
> regressed across the matrix** — not worth the added complexity vs. current priorities (the flagship is
> already ~3% error; higher-priority streamlining: code review, perf, output, docs, real-data
> benchmarking). **Production is v0.7.0 (`origin/main`), which contains none of this.** The branch and
> this doc preserve everything for a clean future resume — see *§0 Outcome & Resumption*.

---

## 0. Outcome, diagnosis & how to resume  *(read this first)*

### What was built (branch `calib-warmstart-v1`, off v0.7.0 `ba6f1f10`)
| commit | phase |
|---|---|
| `8131a558` | design doc |
| `a3e76c7e` | **A** — per-strand node densities (`RnaWarmStart`) on `CalibrationResult` |
| `1f3f2328` + `8d2b5a4d` | **B** — `build_transcript_warm_start` (capture-corrected density bottleneck) + per-strand-spliced AMBIG fix |
| `47d67cd7` + `4d02411d` | **C** — `WS_MODE==4` (`RIGEL_EM_WARMSTART=calib`) in `em_solver.cpp` + length-guard |

Files: `result.py` (`RnaWarmStart`), `calibrate.py` (`_build_rna_warm_start`), `capture_eff_length.py`
(`build_transcript_warm_start`, `_capture_correction`, `_global_depleted_density`, `_pooled_seam_mass`),
`config.py` (`TranscriptGeometry.warm_start_counts`), `pipeline.py` (wiring), `em_solver.cpp`
(`t_warm_start` arg + `WS_MODE==4`), `estimator.py` (`_t_warm_start`). Tests:
`tests/calibration/test_rna_warm_start.py`, `test_transcript_warm_start.py`. **Default path byte-identical
(goldens pass); the flag is default-off.** Two adversarial code reviews (Phase B, Phase C) found no
confirmed bugs.

### Results — ground-truth 3-pool A/B (`quick_3to1_5mb`, `baseline` vs `calib`)
Scored against simulation truth via `scripts/debug/benchmark_ab_report.py` (per-pool surplus =
assigned − observed-true, for gDNA / mRNA / **synthetic** nRNA; + mature-transcript accuracy).

- **Flagship** (`gdna300 ss0.99 nnone capON`): nRNA siphon **103k → 83k (−19%)**, mature `abs_err` better. **Win.**
- **gDNA-free `nrna_none`**: siphon driven to **~0**. **Win.**
- **Net across all 16 conditions: Σ|surplus| +4.1% WORSE.** `calib` is a **nascent-suppressor** — it helps
  where nascent is *over*-called (the siphon) but *deepens* the under-call where real nascent exists
  (`nrna_rnd`), and the freed mass inflates the already-over-called **mRNA** pool. Capture-OFF mRNA is
  slightly worse (the seed perturbs it). The unstranded `ss0.50 capON` 120k siphon stays **flat**.

### Root-cause diagnosis (the key insights to resume with)
1. **Hard `min` is a ceiling estimator** → it *under-seeds* any transcript with one scarce node (a
   heterogeneous real-nascent span; a weak-exon mature under capture-off).
2. **Nascent is pure-ambiguous by construction** — every nascent fragment is co-compatible with gDNA
   (unspliced genomic sequence), so `unambig_totals[nascent] ≡ 0`. Nascent is the **only floor-less pool**
   (mRNA is floored by unambiguous splice-junction reads, gDNA by intergenic reads). Therefore the
   **seed/prior is the sole control** of the gDNA↔nascent split, and **zeroing nascent = permanent death**
   (no data revives it). `calib` *replaces* the coverage seed with the bottleneck, so a `0` bottleneck
   **kills** a real nascent the production default would have kept alive — the `nrna_rnd` regression.

### How to resume — the fix direction
1. **Soft (quantile) bottleneck** in place of `min` (`q` swept, `q=0 ≡` today's hard min): a real
   multi-node nascent gets a nonzero *survival* seed while a fake nascent still → ~0. This is *both* the
   heterogeneity fix *and* the survival floor. **Do NOT "max with the coverage seed"** — that would forbid
   lowering the fake siphon (the coverage seed *is* the over-seed).
2. **Never zero a live component** — keep a small survival seed so the strand/FL likelihood adjudicates
   real-vs-fake nascent rather than the seed pre-deciding death.
3. **Regime-gate** to capture-ON / gDNA-present, so the seed is truly inert on capture-OFF.
4. **Calibration-confidence gating** (`var_pos/var_neg` on `NodeBelief`) for the unstranded-AMBIG regime.
5. **Ceiling:** the seed is only as accurate as calibration's intronic gDNA-vs-nascent deconvolution —
   the unstranded-AMBIG siphon is a *calibration-accuracy* problem, not a warm-start one.

### Resume mechanics
```bash
git checkout calib-warmstart-v1 && pip install --no-build-isolation -e .   # recompile — Phase C touched C++
# ground-truth 3-pool A/B:
python scripts/debug/benchmark_ab_report.py <suite> --out ab.json \
    --arms baseline:0 calib:0::RIGEL_EM_WARMSTART=calib          # render: benchmark_ab_render.py
```
`RIGEL_EM_WARMSTART=calib` enables it; unset = production behavior (byte-identical to v0.7.0).

---

*The full design and phase plan below is unchanged — it documents the approach as built.*

**Status:** design agreed; external critique incorporated; simplified to a minimal "does-it-help"
v1 (2026-07-10); **implemented behind the flag, then DEFERRED (see §0).**

**v1 scope (minimal, but the full node taxonomy — all four node types are load-bearing):**
- Nodes: **region (exon + intron) contained + exon↔intron boundary (nascent crossing) + splice
  junction (mature spliced)**. Both boundary densities are essential:
  - the **exon↔intron boundary crossing** is the captured window into the depleted intron (nascent);
  - the **splice-junction spliced mass** is what recovers mature and separates exon-sharing isoforms.
  They are kept as **separate terms** (crossing vs spliced), routed by node role, so a neighbor's
  mature spliced mass never inflates a nascent shadow's bottleneck (which would reinstate the siphon).
- Operator: **hard `min`** over observable nodes.
- Correction: capture-corrected via gDNA enrichment (`ε`), floored by the depleted gDNA mode.
- Seed: **replace** the RNA-internal coverage seed behind `RIGEL_EM_WARMSTART=calib`; simple
  coverage fallback when a transcript has no observable node.

**Deferred to v2 (add only if the A/B earns it):** transcript-level evidence shrinkage,
calibration-confidence gating (`var_pos/var_neg`), windowed / copy-number-aware `ρ`, explicit
`unambig` floor, soft/quantile bottleneck.

---

## 1. Problem overview

### 1.1 The zombie / siphon problem

Rigel's per-locus EM has `n_t + 1` components — one per transcript row (annotated mRNA **and**
synthetic nascent-shadow spans alike) plus one gDNA component. Under hybrid capture with genomic
contamination, mass that is truly gDNA (or truly mature) leaks onto the **synthetic nascent shadow**
components — the "siphon." Two mechanisms feed it:

1. **The warm start is structure-blind.** The production warm start
   (`compute_grouped_warm_start`, `em_solver.cpp:913`) seeds each component from its unambiguous
   fragment count plus coverage-weighted shares of the ambiguous fragments, then projects through
   the pool-level gDNA/RNA prior. A nascent shadow spans every exon of its gene, so it grabs
   coverage shares of all exonic fragments and is seeded with real mass **even when its true
   abundance is zero**. Effective length never enters the init.
2. **The prior is rich-get-richer.** `apply_grouped_prior_update` (`em_solver.cpp:714`) re-applies
   every EM step and distributes the RNA pseudocount proportional to the current counts
   (`out_counts[i] = R · raw_counts[i] / n_rna`). It **cannot revive a zero-count transcript** (good)
   but **amplifies whatever the warm start seeded** (bad). A wrongly-seeded shadow snowballs.

The calibration prior today is only two per-locus scalars — the aggregate gDNA-vs-RNA split
(`assemble_priors`, `priors.py:196`). The per-transcript distribution of the RNA budget is left
entirely to fragment coverage. **Calibration produces no per-transcript abundance today.**

### 1.2 Goal

Replace the RNA-internal seed with a **calibration-informed per-transcript abundance** that
distributes the RNA prior mass non-uniformly (keeping the RNA total), built from the per-node
deconvolution calibration already computes — so structurally unsupported components start near zero
and the rich-get-richer prior amplifies the *right* transcripts.

---

## 2. The core principle — the density bottleneck

A node's density is the sum of **all** sources at that node. A transcript `T` contributes its own
density to every node it covers, so for every node `n ∈ T`:

```
node_density[n] ≥ density(T)   ⇒   density(T) ≤ min over n∈T of node_density[n]
```

The **minimum node density is a true ceiling** on the transcript's density (a max-flow / bottleneck
bound), and — because when `T` dominates its nodes the min ≈ `density(T)` — it is also a good point
estimate. This is the warm-start value.

### 2.1 Why the bottleneck is the nascent gate — and why boundaries are essential

The node taxonomy carries the mature/nascent signal, and `min` selects each transcript's binding
constraint:

| Node type (v1) | Density term | Bottlenecks |
|---|---|---|
| **Exon, contained** | `f_s·M_u / eff_rna` (unspliced) | mature **and** nascent (shared ceiling) |
| **Intron, contained** | `f_s·M_u / eff_rna` (unspliced) | nascent only |
| **Exon↔intron boundary** | **nascent crossing** `f_s·M_side / eff_rna_side` | nascent only |
| **Splice junction** | **mature spliced** `spliced_s / eff_spl` (half-triangle) | mature only |

- **Nascent N** → `min(exon contained, intron contained, exon↔intron boundary crossing)`. The exon
  term is high (that density is really M's), so `min` picks an intron/boundary — ~0 without real
  nascent.
- **Mature M** → `min(exon contained, splice-junction spliced)`. The **splice-junction spliced mass**
  is what recovers mature and — because exon-sharing isoforms differ only in their junctions — what
  separates isoforms that a contained-exon bottleneck alone could not distinguish.

**The two boundary densities are kept as separate terms, never summed.** A splice donor/acceptor site
is *also* an exon↔intron boundary, so a single boundary node carries both a nascent crossing term
(`f_s·M_side/eff_rna_side`, ≈0 without real nascent) and a mature spliced term (`spliced_s/eff_spl`,
high). If a nascent shadow's bottleneck saw the *sum*, a neighbor's mature spliced mass would inflate
that node → and with a depleted intron interior the shadow would be over-seeded → the catastrophic
siphon. So nascent bottlenecks on the crossing term only; mature bottlenecks on the spliced term. Both
come from the same `belief + geometry` (`node_densities` computes exactly these two pieces,
[node_geometry.py:237](../../src/rigel/calibration/node_geometry.py) — we expose them separately
rather than as its sum). The `eff_rna`-vs-`eff_spl` basis difference is handled by the
factor-1-under-uniform eff-length construction (each density estimates the true ρ), so they are
commensurable in the `min`.

**Boundaries are load-bearing under capture.** The intron *interior* is off-probe → depleted, often
to zero total mass → **excluded by the observability mask**. If the bottleneck saw only contained
regions, a nascent shadow with a fully-depleted intron would lose its only distinguishing node and
fall back to the shared exons → seeded high → **catastrophic siphon**. The exon↔intron **boundary**
survives because its exon flank is captured, so a genuine unspliced (nascent) read crossing it is
pulled down and *observable* even when the intron interior is empty. That captured crossing is the
nascent signal under capture.

The transcript→node incidence (`_transcript_node_incidence`, `capture_eff_length.py:43`) already
produces these node sets: a multi-exon mRNA → its exon regions; a synthetic single-exon nRNA span →
all its regions (incl introns) **and** the interior boundaries between them.

---

## 3. The complication — hybrid capture

Under capture, off-probe regions are **depleted**: their observed density is low regardless of the
underlying RNA. A naive `min` would bottleneck every transcript at its most-depleted node. We
correct for depletion using an *endogenous* control.

### 3.1 gDNA enrichment *is* the local capture efficiency

gDNA is genomic substrate — present genome-wide at ~uniform concentration, no exon/strand structure.
So the observed gDNA density at a node is a direct readout of how efficiently capture pulls down that
locus:

```
ε_node = gdna_density[node] / rho_ref      ∈ (0, 1]     (1 = fully on-target)
corrected_rna[node] = rna_density[node] / ε_node = θ   (true, capture-removed)
```

`rho_ref` is the **enriched** gDNA density mode, reused from `_global_reference_density`
(`capture_eff_length.py:155`). Because `ε` is measured from gDNA **at the same node** as the RNA, the
correction is local — and copy-number-robust (§4.4).

---

## 4. Robustness

### 4.1 Observability (the two zero-types)

- **gDNA present, RNA ≈ 0** → observable, genuinely RNA-empty → a **real bottleneck** (zeros a fake
  nascent).
- **gDNA ≈ 0 and RNA ≈ 0** → we observed nothing → **unobservable** → **excluded** (micro-region /
  low-mappability / sparse-count zero). Threshold-free: a node enters iff **total observed mass**
  (`gDNA + RNA`) `> 0`.

### 4.2 Two correction guards (the dangerous direction is *over*-correction)

The failure that keeps the siphon alive is **over-correction** in low-background regions
(`ε → 0` → `1/ε` huge → RNA noise blown up → fake nascent). Under-correction (high background → little
lift) is self-safe. Both guards target over-correction; neither adds a new constant:

| Failure mode | Guard | Source |
|---|---|---|
| Depleted node, **few** gDNA fragments → noisy density | `w = C/(C+1)`, C = node total observed mass | reuses existing shrinkage (`capture_eff_length.py:319`) |
| **Long** depleted region, moderate count but tiny density → `1/ε` blows up | `ε_floor = rho_depleted/rho_ref` | `rho_depleted` = depleted gDNA KDE mode (new `_global_depleted_density`, symmetric to `rho_ref`) |

### 4.3 Zero / sparse-gDNA — graceful by construction

`gDNA → 0` happens only in clean libraries — where we neither need the correction nor have a siphon.
The algorithm degrades to the no-capture baseline (`min` of raw stranded RNA density) automatically:
`rho_ref` is `None` when `< 5` gDNA-observed nodes (`capture_eff_length.py:169`) → no correction;
a unimodal gDNA field ⇒ `rho_depleted → rho_ref` ⇒ `ε_floor → 1` ⇒ no-op. **No synthetic gDNA** — it
would corrupt the real gDNA masses that *are* `gdna_prior_count`, and its magnitude is an
unprincipled knob.

### 4.4 Copy-number variation (tumor) — cancels via per-locus renormalization

`gdna_density = copy_number × capture_efficiency × const`; a `k×` amplification inflates gDNA density
→ `ε` overestimated → RNA under-corrected. **This largely cancels**: each locus is an independently
renormalized EM subproblem, and CNV segments are far larger than a gene, so within a locus copy
number is ~constant:

```
under k×:  gdna_density_n → k·gdna_density_n,  observed rna_density_n → m·rna_density_n
  ⇒  ceiling(Tᵢ) → (m/k)·ceiling(Tᵢ)  — a common factor across ALL transcripts in the locus
  ⇒  cancels in the per-locus RNA renormalization.
```

It works *because* `ε` is measured at the same nodes as the RNA: capture varies *within* a locus
(exon vs intron — the signal, doesn't cancel), CNV is constant *across* a locus (a nuisance, cancels).
Residuals: `ε`-clipping (amplification → introns under-lifted → nascent under-seeded, safe;
**deletions** → over-correction, bounded by `ε_floor`/`C/(C+1)`); intra-locus CNV breakpoints (rare);
the *pool* gDNA-vs-RNA split is CNV-sensitive but that is existing calibration, not this warm start.

---

## 5. The algorithm (v1, precise)

For each transcript `T`, strand `s` from `t_to_strand_arr`, over `T`'s region + interior-boundary
incidence nodes:

```
# --- per-node capture-corrected density on strand s (term chosen by node ROLE) ---
rho_rna    = region contained     : f_s·M_u    / eff_rna        (mature + nascent)
             exon↔intron boundary : f_s·M_side / eff_rna_side   (nascent crossing)
             splice junction      : spliced_s  / eff_spl        (mature spliced)
             # crossing and spliced are SEPARATE terms on the same boundary node — never summed
rho_gdna   = mass_gdna_* / gdna_*_eff_len                     (contained or per-side, on the result)
ε          = clip(rho_gdna / rho_ref, ε_floor, 1),   ε_floor = rho_depleted / rho_ref
w          = C/(C+1),   C = total observed mass (gDNA+RNA) at the node
corr       = rho_rna · (1 + w·(1/ε − 1))
observable = (total observed mass at node) > 0

# --- transcript ceiling and warm-start count ---
ceiling(T)   = min over observable nodes of corr              # hard min
warm_raw[T]  = ceiling(T) · E_T^em                            # density → observed count
```

**Count conversion is exact.** `corr ≈ θ_T` (true abundance, capture removed) and `E_T^em`
(`effective_lengths_em`) is the *capture-contracted* observable effective length, so
`θ_T · E_T^em` = the expected **observed** fragment count the VBEM posterior accumulates — the right
scale to seed the Dirichlet α. Only relative `warm_raw` matters (the pool prior renormalizes RNA to `R`).

**Fallbacks.** `rho_ref is None` → `ε = 1` → bottleneck of raw stranded RNA density (no-capture
baseline). A transcript with **no observable node** → keep the existing coverage seed (never zeroed).

**No new magic numbers:** `_KDE_BW=0.4`, `_KDE_PROM=0.05`, the `<5`-node guard, and `C/(C+1)` are all
reused; `ε_floor` is data-derived from the depleted mode.

---

## 6. Implementation plan

All calibration work is Python; the only C++ change passes and consumes one per-transcript vector,
mirroring `t_eff_lens`.

### Phase A — expose the signals (Python calibration; default output byte-identical)

- **A1 · Per-strand node densities on the result — three terms kept separate.** At result-construction
  in `calibrate.py` (`:232-248`), from the final `belief.f_pos/f_neg` (raw for **all** nodes) + the
  static `geometry`, compute and persist region-indexed per-strand fields on `CalibrationResult`
  (`result.py:38-84`, validate at `:91-102`):
  1. **contained** `f_s·M_u/eff_rna` (region node);
  2. **boundary crossing** `f_s·M_side/eff_rna_side` (nascent — the unspliced first term of
     `node_densities`);
  3. **junction spliced** `spliced_s/eff_spl` (mature — the second term of `node_densities`,
     motif-stranded, half-triangle basis).
  `node_densities` (`node_geometry.py:237`) computes (2)+(3) summed; we expose them **separately** so a
  neighbor's mature spliced mass cannot inflate a nascent shadow's crossing bottleneck. **No
  `bp_solver` change** — we read the raw belief, not the (boundary-zeroed) `NodeDeconv` projection.
- **A2 · Depleted mode.** Add `_global_depleted_density` (leftmost significant KDE peak, mirrors
  `_global_reference_density`; same `<5` → `None` guard). `rho_ref` and `rho_depleted` are recomputed
  in the builder from result fields — no persistence needed.
- **A3 · Regression gate.** New fields are unread by the default path; the `node_densities` call is
  additive and does not alter existing arrays. Confirm `pytest tests/ -q` (1086) green, golden
  feather/TSV/JSON **byte-identical**.

### Phase B — the warm-start builder (Python)

- `build_transcript_warm_start(calibration, region_arrays, index, effective_lengths_em) →
  float64[n_transcripts]`, beside `transcript_capture_eff_lengths`, sharing all its
  `_transcript_node_incidence` arrays (region, interior-boundary, **and splice-junction**; avoid a
  second `intervals.feather` read). Implements §5: strand-select; per-node density chosen by **node
  role** — contained (region) + crossing (exon↔intron boundary, nascent) + spliced (splice junction,
  mature); capture correction (`rho_ref` reused, `ε_floor` from `rho_depleted`, node `C/(C+1)`); hard
  `min`; `× E_T^em`. A transcript with no observable node → emit a NaN/sentinel meaning "use coverage
  seed."
- Called in `pipeline._setup_geometry_and_estimator` (`pipeline.py:425`, before `assemble_priors`);
  the vector stored on the estimator like `_t_eff_len_em`.

### Phase C — plumb to C++ and add the mode

- `estimator.run_batch_locus_em_partitioned` (`estimator.py:282`): add a `t_warm_start` keyword, cast
  beside `t_eff_lens` (`:341`).
- `em_solver.cpp batch_locus_em_partitioned` (`:1935`, binding `:2678`): add one `f64_1d
  t_warm_start_arr` + `nb::arg`; gather per component in `extract_locus_sub_problem_from_partition`
  (`:1760`, mirror the eff-len gather at `:1897`) into `sub.warm_seed[nc]`. The array is **read-only**
  during the parallel EM (no false-sharing on the read path); keep it contiguous, indexed identically
  to `t_eff_lens`.
- `compute_grouped_warm_start` (`:913`): add `WS_MODE == 4` (`"calib"`) — for each RNA component with a
  finite `warm_seed`, use it as `warm_raw` (replacing the coverage seed); otherwise keep the coverage
  seed; the gDNA component (index `n_t`) keeps its coverage seed. Then the existing
  `apply_grouped_prior_update` projection. All other modes and the default stay byte-identical.
  Recompile: `pip install --no-build-isolation -e .`.

### Phase D — A/B (see §8)

---

## 7. Testing

**Unit (`tests/calibration/`):** `build_transcript_warm_start` — (a) mature-uniform-exon /
nascent-zero-intron toy → shadow ceiling ≈ 0, mature ≈ exon density; (b) **capture case: depleted
intron *interior* excluded, but the exon↔intron boundary carries the nascent signal** → shadow
ceiling tracks the boundary; (c) capture correction lifts a real depleted boundary, `ε_floor` caps a
deep off-target region; (d) `rho_ref=None` → raw-`min`; (e) strand selection (a `+` transcript ignores
`−` density in an AMBIG overlap); (f) no-observable-node → coverage fallback (never 0); (g) CNV
invariance — scaling a locus's gDNA by `k` leaves the *relative* per-transcript seed unchanged (up to
`ε`-clipping). `_global_depleted_density` — unimodal → `≈ rho_ref`; bimodal → distinct mode; `<5` →
`None`.

**Regression & toy:** `pytest tests/ -q` green; goldens **byte-identical** with the mode default-off.
`scripts/debug/toy_prod.py` production-faithful toys (incl. ample single-stranded nodes).

---

## 8. Validation (A/B protocol)

- **Baseline:** `RIGEL_EM_WARMSTART` unset. **Arm:** `RIGEL_EM_WARMSTART=calib`. Separate processes
  (config file-flips don't reach a subprocess); delete `rigel_out` to force re-quant.
- **Metric:** soft 3-pool — synthetic nascent EM count (`est.nrna_em_count`, the siphon) + mature
  accuracy (abundant-tx Spearman, n_FP). Never `net_flow`.
- **Suites (easy → hard):** `quick_3to1_5mb` (incl. the flagship
  `gdna_gdna300_ss_0.99_nrna_none_capture_on`) → `ambig_dense_10mb`. Use the oracle
  (`scripts/debug/oracle.py`) to separate a calibration error from a warm-start effect.
- **Success:** the flagship siphon drops materially, mature accuracy flat-or-better, **no regression**
  on capture-off / gDNA-free / stranded-clean conditions.

---

## 9. Risks & mitigations

| Risk | Mitigation |
|---|---|
| **Mixed-basis min** — contained (`÷eff_rna`), crossing (`÷eff_rna_side`), spliced (`÷eff_spl` half-triangle) mixed in one `min` | Handled by the factor-1-under-uniform eff-length construction (each density estimates the true ρ). Crossing (nascent) and spliced (mature) also never mix *within one transcript's* `min` — routed by node role. Validate on toys. |
| **Spliced mass contaminating nascent** — a splice site is also an exon↔intron boundary | Crossing and spliced kept as separate terms; nascent uses crossing only, mature uses spliced only. |
| **`rho_ref` mode-snapping** — discrete jumps at low N | Keep both modes on the same mass-weighted contained basis; the `ε_floor` *ratio* is steadier than either endpoint; validate `ε_floor→1` on capture-off toys. |
| **Noise amplification in deep depletion** | Bounded by `ε_floor` + node `C/(C+1)`. |
| **Copy-number variation (tumor)** | Locus-constant CNV cancels in per-locus renormalization (§4.4); amplification → under-seed (safe); deletions bounded. |
| **Inherited calibration error** (AMBIG-unstranded, −260k) — the seed trusts the deconvolution | v2: gate node contributions by `var_pos/var_neg` (on `NodeBelief`, free). |
| **Observability cliff** for thin transcripts | Accepted for v1; v2 transcript-level shrinkage smooths it. |

---

## 10. Deferred to v2 (add only if the A/B earns it)

- **Transcript-level evidence shrinkage** — blend the ceiling toward the coverage seed by
  `w_T = C_T/(C_T+1)`; removes the observability cliff, protects thin transcripts (either direction).
- **Calibration-confidence gating** — down-weight a node by its deconvolution variance
  (`var_pos/var_neg` on `NodeBelief`), limiting garbage-in from AMBIG-unstranded loci.
- **Windowed / copy-number-aware `ρ`** — chromosome / multi-locus-scale `rho_ref`/`rho_depleted` for
  macroscopic background & CN variation.
- **Explicit `unambig` floor** — seed in `[unambig-implied density, ceiling]` so it never falls below
  hard evidence.
- **Soft bottleneck** — a swept mass-weighted low quantile (`q=0 ≡ hard min`) if a real transcript is
  crushed by one node.

---

## 11. Open questions

1. **Mixed-basis commensurability** (§9 row 1) — confirmed contained + boundary share the `eff_rna`
   basis; revisit before adding splice junctions.
2. **Calibration-confidence gating** — `var_pos/var_neg` are on `NodeBelief`; is variance-weighting the
   right form, and does it need to reach the result?
3. **gDNA-component seed consistency** — RNA transcripts seed from the bottleneck, gDNA from coverage;
   verify the magnitudes compose sanely through `apply_grouped_prior_update` (likely fine — the pool
   prior sets gDNA independently).

---

## Appendix — key code references

| Concern | Location |
|---|---|
| Current warm start | `em_solver.cpp:913` `compute_grouped_warm_start` |
| Pool prior projection | `em_solver.cpp:714` `apply_grouped_prior_update` |
| C++ EM entry / extract / gather | `em_solver.cpp:1935` (binding `:2678`) / `:1760` / `:1897` |
| Warm-start mode dispatch | `em_solver.cpp:924` (`RIGEL_EM_WARMSTART`) |
| Per-strand beliefs / densities / variances | `node_geometry.py:203` `NodeBelief` (incl. `var_pos/var_neg`), `:218` `node_densities` |
| Enriched reference density | `capture_eff_length.py:155` `_global_reference_density` |
| Transcript→node incidence | `capture_eff_length.py:43` `_transcript_node_incidence` |
| Result schema / construction | `result.py:38-84` (validate `:91-102`) / `calibrate.py:232-248` |
| Per-transcript EM eff-length | `pipeline.py:425` `_setup_geometry_and_estimator` → `effective_lengths_em` |
| Estimator → C++ | `estimator.py:282` / `:341` / `:394` |
| Transcript strand | `index.t_to_strand_arr` (POS=1/NEG=2, defined for nRNA) |
