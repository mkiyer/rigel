# Calibration code audit

## Executive summary

The calibration package is **algorithmically sound and functionally healthy** post-collapse. Every auditor confirms the same thing: the core data flow (payload → substrate → chain → geometry/statics → forward-backward sweep → projections → priors) is correct, the single belief-free Poisson disagreement-variance precision path is cleanly implemented, and the FB message combination does no double-counting. **No correctness bugs surfaced during this (structural) pass** — a clean bill on structure, not a formal correctness proof. The debt is almost entirely *residue from the collapse*: dead code, stale docstrings that describe removed machinery, and precomputed index structures whose only consumers were deleted. The user's "dire need of a cleanup" is justified — the cruft actively misleads (docstrings assert the opposite of what the code does) — but the fixes are overwhelmingly low-risk deletions and doc rewrites.

**The 3–5 biggest wins:**

1. **Kill the dead spliced-floor path.** `build_node_statics` hard-zeroes `statics.spliced_pos/neg` (bp_solver.py:454) with a self-admitted "dead-code cleanup is a follow-up" comment; the entire node-local spliced-lower-bound term threads zeros through `NodeStatics`, `_region/_boundary_strand_stats`, and all four `simplex_logodds` solver functions including the AMBIG cube. Complete the promised follow-up. (Raised by 3 clusters.)
2. **Delete the dead index precompute.** `BoundaryArrays` has **zero callers in src/**; the per-boundary structural flags (`is_tss/is_tes/is_splice_junction/genomic_sj_strand`) and per-region `mature_eligible_{pos,neg}` are built, persisted, validated, and threaded through `RegionArrays` but read only by tests — they fed the removed mature/nascent overlay. (Raised by 3 clusters.)
3. **Rewrite the stale algorithm docstrings.** `calibrate.py`'s module header and `CalibrationConfig` still describe "Gauss-Seidel passes", "per-pass frozen-snapshot var~mean reliability", "ρ_global re-fit each pass", and a "2-simplex lattice" — none of which the single-FB-pass code does. Multiple other docstrings (bp_solver ê(z), gdna_density_prior "not wired in yet", effective_length/fl.py "RNA cannot be FL-corrected") assert removed or inverted facts.
4. **Split bp_solver.py (1328 lines).** A cohesive ~350-line gDNA global-prior subsystem is a natural separate module; extracting it also dissolves the awkward `gdna_density_prior.py → bp_solver.py` back-import.
5. **One source of truth for coarse type.** `_node_region_type` (bp_solver.py:525-534) reimplements `signature.coarse_type_array` inline; its own docstring asserts the equivalence.

**INDEX-LEVERAGE VERDICT: PARTIAL, and mostly the inverse of the hypothesis — REMOVE, don't wire in.** Exactly one runtime derivation should switch to the canonical column: `_node_region_type` → `signature.coarse_type_array` (behavior-preserving, equivalence already asserted). Everything else runs opposite: the v6 boundary columns (`is_tss/is_tes/is_splice_junction/genomic_sj_strand`) and `mature_eligible_{pos,neg}` are **dead precompute the solver already correctly bypasses** and should be deleted, not read. The two derivations the prompt flags as switch candidates — the `_spliced_faces` junction routing and the TSS/TES continuity gate — must **NOT** switch: the accumulator's motif-observed `junction_strand` (domain {0,1,2}) is deliberately more correct at AMBIG seams than the annotation-derived `genomic_sj_strand` ({0,1,2,3}), and the continuity gate `free_s = has_s(left)&has_s(right)` is strictly richer than an `is_tss/is_tes` lookup.

---

## Priority 1 — high-value, low-risk (do first)

| Area | Issue | Recommendation | Files | Effort |
|---|---|---|---|---|
| Docs — orchestrator | Module docstring + "THE SOLVE" comment describe retired per-pass var~mean / Gauss-Seidel loop and omit the actual two-pass (pass-1 no-prior / pass-2 KDE-prior) structure | Rewrite flow to: strand balance → overdispersions → build chain/geometry/init → compute scalar σ²_imp once → pass-1 (prior=None) → train `GdnaDensityPrior` → pass-2 (mixture prior) → deconv → `gdna_density_global`. Drop all "var~mean reliability"/"ρ_global re-fit each pass" language | calibrate.py:1-28, 149-154 | low |
| Docs — config | `CalibrationConfig` class docstring says "2-simplex lattice"/"Gauss-Seidel passes"; `sweep_max_passes`/`sweep_convergence_delta` document a nested refit loop | Rewrite to single FB pass w/ belief-free σ²_imp; verify `sweep_n_tilt` stays (live: AMBIG cube) | config.py (CalibrationConfig) | low |
| Docs — bp_solver | Docstrings reference removed ê(z), mature/nascent overlay, "Gauss-Seidel", and claim the spliced floor "is a node-local term" (contradicts the zeroing) | Delete ê(z) mentions (21-24, 970), fix the "node-local floor" clause at 979 to "message-carried only" | bp_solver.py:21-24, 968-970, 979 | low |
| Docs — prior | `gdna_density_prior.py` says "Nothing here is wired into the solve yet (P2.4)" — it IS the production Phase-2 prior consumed by `_kde_logprior` | Rewrite preamble to state it is the production Phase-2 prior; drop P2.2/P2.4 roadmap framing | gdna_density_prior.py:15-17 | low |
| Docs — efflen/fl | Both module docstrings assert only gDNA consumes a modelled FL and "RNA cannot be FL-corrected / PR-5 future" — false: `build_node_geometry` calls the eff-length primitives with `rna_fl_pmf` | State BOTH gDNA and RNA FL pmfs drive per-node eff-lengths in the sweep | effective_length.py:6-8, fl.py:9-11 | low |
| Docs — node solver | `simplex.py`'s 34-line docstring describes a 4-term per-node solver that no longer lives there and cites an archived spec; `_solve_nodes_logodds_all` says messages are "still linear-fraction" when the code evaluates log-fraction Gaussians | Rewrite `simplex.py` docstring to the surviving strand-mixture primitive only; correct the "still linear-fraction" clause | simplex.py:1-34, simplex_logodds.py:404-412 | low |
| Docs — priors invariant | "bedrock invariant" paragraph states `S_s = E[ℓ]` / deposits `ρ·E[ℓ]`, contradicting the code (`seam_support = 0.5·(side_density_len...)`, line 251) and the same file's own text | Replace with `S_s = ½·(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})])` — pure doc fix, code is right | priors.py:303-316 | low |
| Index-leverage / consolidation | `_node_region_type` reimplements `signature.coarse_type_array` inline; local `_EXON_BITS/_INTRON_BITS` exist only for this duplicate | Call `signature.coarse_type_array(region_arrays.signature)`; keep only the node-scatter local. Behavior-preserving | bp_solver.py:525-534 (vs signature.py:152-165) | low |
| Dead-code — params | `node_sweep` `max_passes`/`convergence_delta` are documented "unused (FB is single-pass)" but threaded live from config | Drop both params + call-site args; retire/deprecate the two config fields + validators; fix config docstring | bp_solver.py:937,949-950,985; calibrate.py:181-182; config.py | low–med |
| Dead-code — BoundaryArrays | `BoundaryArrays` + `from_boundary_df` have zero callers in src/ or tests | Delete the dataclass + factory (pure dead code, zero production risk) | region_arrays.py:132-193 | low |
| Dead-code — orphan field | `NodeDensity.global_density` FIELD is dead (decl :56-58, assignment :174) | ⚠️CRITIC-CORRECTED: delete the FIELD only. **KEEP** the local `global_density` computation at :160 — it IS consumed at :161 (`np.where(anchored, density, global_density)`) to fill un-anchored regions. Deleting the computation breaks the density fill. | density_model.py:56-58, 174 | low |
| Dead-code — unreachable branch | `_global_logprior`: `node_sweep` always passes `s2_bg` non-None (a float from `_floor_estimate`), so the `is not None` guards are always True and the "legacy capped floor (s2_bg is None ⇒ byte-identical)" branch is UNREACHABLE dead code with a stale comment asserting the dead equivalence — in the hot prior-assembly path | Delete the unreachable `s2_bg is None` branch, collapse the always-true guards, drop the stale "byte-identical" comment | bp_solver.py:782, 799, 802 | low |
| Dead-code — QC field | `StrandBalance.rna_strand_overdispersion` (1/(n+3) QC diagnostic) never read in production; costs ~20 docstring lines managing a naming collision with the live overdispersion | Delete field + computation (66-71) + the "do not confuse" docstring block; update tests | strand_balance.py:47, 66-71, 12-30 | low |
| Dead-code — struct fields | `_SideQuantities` populates `antisense`, `mass`, `mass_spliced`, `count_gdna_frac_var` (a whole Poisson-floor block); only `sense/n_side/strand_observable/count_gdna_frac/count_observable` are read | Drop the four unread fields + the variance block + its docstring; `_side_strand_orientation` returns `(sense, n_side, strand_observable)` | strand_deconv.py:46-58, 124-127, 141-145 | low |
| Dead-code — field | `CalibrationSubstrate.region_len` populated, never read (solver reads `region_size_bp` directly) | Remove field + assignment; adjust the one test | substrate.py:88,118 | low |
| Dead-code — export | `_floor_log_density` defined + in `__all__`, never called anywhere | Delete it + its `__all__` entry; move the 1/(n+1) rationale to the inline comment at simplex_logodds.py:335 if canonical | simplex_logodds.py:73,30 | low |
| Algorithm — hotspot | `node_gdna_density` uses a Python `for i in np.where(~own)[0]` loop over the *majority* of nodes (all exonic/AMBIG) at genome scale | Vectorize with mask arithmetic (la/ra anchors → 0.5·(l+r), l, or r; else NaN). Behavior-identical | density_model.py:141-148 | low |
| Consolidation | `(2κ−1)²` strand-discriminability weight written inline twice + referenced in a comment | Add one helper `strand_discriminability(rna_sense_frac)`; call from both sites | bp_solver.py:586,711; calibrate.py:107 | low |
| Dead-code — constant | `_MSG_PSEUDOCOUNT = 1.0` defined, never referenced; the real pseudocount is a hardcoded `+ 1.0` | Use the named constant in the `pr = n_src/(n_src·σ²_imp + 1.0)` denominators (per no-magic-number ethos) | bp_solver.py:76-77, 1169/1194/1210 | low |
| Docs — dead xrefs | `_fit_seed_varmean` "twin of `_edge_varmean`" (removed); `_poisson_moment_var` "shared by v1/v2 per-component estimators" (v2 deleted); "v1" qualifiers with no v2 | Delete the dead cross-references and orphaned "v1"/"v2" labels | bp_solver.py:652, 890-895, 922 | low |

---

## Priority 2 — high-value, higher-effort-or-risk (worth doing, plan carefully)

- **Split bp_solver.py (1328 lines).** Two auditors converge on extracting the cohesive gDNA global-prior subsystem (`_gdna_seed_estimate`, `_LogLinearVarMean`, `_fit_seed_varmean`, `_floor_estimate`, `_global_logprior`, `_kde_logprior`, and the shared `_node_region_type`) into a new module (e.g. `gdna_prior_terms.py` / `gdna_global_prior.py`). This drops bp_solver to ~950 lines **and dissolves the `gdna_density_prior.py → bp_solver` back-import** (it currently reaches in for `_node_region_type` and `node_densities`). A second auditor additionally proposes hiving off the geometry layer (`NodeGeometry`/`build_node_geometry`/`node_densities`/`NodeStatics`/strand_stats/init_beliefs). **Recommendation:** do the prior-subsystem extraction first (highest leverage, kills the cycle); treat the geometry split as an optional follow-up. Pure code move; risk is import churn + a golden re-check. See Module decomposition below. *(Effort med, Risk med.)*

- **Remove the dead spliced-floor path (the #1 residue item).** `statics.spliced_pos/neg` are zeroed (bp_solver.py:454); complete the follow-up: drop the fields from `NodeStatics`, remove the derivations that feed them (`_boundary_strand_stats._floor`/ls/rs, `_region_strand_stats` `spl`/spliced_*), remove `_zero_spl`, and remove the `spliced_pos/spliced_neg` params + `short_p/short_n` floor terms from the four `simplex_logodds` functions (1-D and 2-D AMBIG cube). **Scope guard:** the *geometry* `spliced_*` (via `node_densities` + `_scan`'s SPs/SNs message terms) is LIVE and separate — touch only the STATICS/node-local floor. Expected byte-identical (zeros in ⇒ zero contribution) — but PROVE it with a golden re-check, don't assert it up front. *(Effort med, Risk med — threads several signatures.)*

- **Delete `mature_eligible_{pos,neg}`** (`_compute_mature_eligible` ~60 lines, `REGION_MATURE_COLUMNS`, the `RegionArrays` fields + `_mature_col`, feather schema). Dead post-collapse (fed the removed mature/nascent split); read only by tests. **Index-schema/golden caveat:** bump the region feather schema and regenerate `regions.feather` goldens; runtime behavior unchanged. *(Effort med, Risk med.)*

- **Retire the boundary structural-flag columns** (`is_tss/is_tes/is_splice_junction/genomic_sj_strand`) and decide the `boundaries.feather` lifecycle. `validate_boundaries_against_regions` already asserts boundary positions equal region interfaces (so positions carry no new info), and the flags have no runtime consumer. Drop the flag columns from `build_boundary_partition`/`BOUNDARY_COLUMNS`; then decide whether the whole `build_boundary_partition`/`load_boundaries`/validate lifecycle is worth keeping. Larger index-schema change; sequence after the BoundaryArrays delete. *(Effort med, Risk med.)*

- **Precompute the transcript→node incidence at index build.** `capture_eff_length._transcript_node_incidence` reads `intervals.feather` at runtime and recomputes each transcript's region/boundary incidence via per-exon searchsorted; its own docstring concedes it is "annotation-only — could be precomputed at index build." Move it to a companion feather keyed by `t_index`. Behavior-preserving (a divergence surfaces as a golden change); at minimum hoist the feather read out of any repeated call path. *(Effort med, Risk med.)*

- **Consolidate the triple-documented region+seam gDNA node model** (`_gdna_region_node_arrays` docstring, `assemble_priors` docstring, the inline comment block — ~90 lines of near-duplicate prose that already drifted into the `E[ℓ]` contradiction). Keep the authoritative derivation in `_gdna_region_node_arrays` only; reduce the other two to a one-liner + "see `_gdna_region_node_arrays`". *(Effort low, Risk low — but grouped here as it pairs with the invariant fix.)* priors.py:190-232, 285-323, 330-336.

- **Decide the fate of the parked experiments** (`RIGEL_RHOSTAR=kmeans` — a 64-line `_per_locus_kmeans_rhostar` + branch that also re-implements the log-space shrinkage blend; `RIGEL_ORACLE_CALIB` override in the hot return path). Both are never-default env-gated study code. **Confirm with the user** each study is closed, then delete (git preserves them); if kept, move behind a clearly-labeled experimental section / debug helper. priors.py:51-114, 414-429; calibrate.py:249-264. *(Effort low, Risk low, but needs a user decision.)*

---

## Priority 3 — nice-to-have / low-value (candidates to skip)

- **`_posterior_median_fg` study branches** (mean/parabola/cdf gated on `RIGEL_MEDIAN_MODE`, never set; ~25 lines + two bare constants `1e-300`, `-1e-12`). If the Fix-3 grid study is closed, collapse to `snap`. simplex_logodds.py:90-131.
- **`simplex.py` module merge** — 64-line file, one function (`_mixture_strand_loglik`), one caller. Move into `simplex_logodds.py` and delete, or rename. Only `from .simplex import ...` changes.
- **`strand_loglik` vs `_mixture_strand_loglik`** — two Beta-Binomial-moment implementations of the same math; `strand_loglik` has no production caller (tests only). Consolidate only if touching the area. strand_likelihood.py:27 vs simplex.py:45.
- **`beta_concentration()` duplicated** byte-identically on both strand models, test-only callers → single free function beside `overdispersion_for_beta`. gdna_strand.py:89-96, 110-113.
- **`strand_summary.py`** unused mirror properties (`strand_specificity`, `read1_sense`), an over-general `strand_contrast_identifiable` reached only via `is_identifiable`, and a heavy 25-line `__post_init__` on a QC-only object built from an already-trusted model. Trim.
- **Dead scalar signature helpers** — `is_ambiguous_signature` + `coarse_strand_from_signature` (a fully dead pair), scalar `coarse_type_from_signature` (docstring-only reference), possibly `validate_signature`. Inflate `__all__`. signature.py:110-149.
- **`_capture` diagnostic hook** in `node_sweep` (~25 lines, runs a second full strand-only solve; default None). Confirm no shipped path sets it, then remove. Medium risk (external debug tooling). bp_solver.py:1240-1263.
- **`node_chain.build_node_chain` Python double loop** — O(n_nodes) pure Python; `region_arrays` already encodes the B-R interleave vectorized. Defer unless genome-scale build time bites. node_chain.py:80-103.
- **`_component_message` consolidation** — `emit_p`/`emit_n` are near-identical 15-line copies; the precision formula + density-message projection appear three times. Factor one `_component_message` helper (best folded into the module split). bp_solver.py:1158-1214.
- **`_clean_exon_boundary` duplicate** of the clean-boundary derivation in `_gdna_seed_estimate` (docstring admits "consolidate when Phase 2 lands" — it has). Extract one shared helper. gdna_density_prior.py:69-88.
- **`_build_global_logprior` extraction** from `node_sweep` (lines 1042-1092) — subsumed by the module split.
- **Extract `_local_loglik_logodds` micro-redundancy** — `_log_fg`/`_log1m_fg` evaluated twice per call; reuse the line-210/211 results. Bit-identical, one fewer log_expit pair. simplex_logodds.py:210-227.
- **Magic-number discussions** (per ethos, flag-for-discussion, not urgent): `_MIN_KDE_TEACHERS = 10` (move to config alongside other KDE knobs), `_kde_logprior` `0.5`/`99.5` percentile bridge tails, `capture_eff_length` shrinkage `c_ev+1` + `1e-9` floors (the "no new constant" docstring claim is inaccurate). calibrate.py:78; bp_solver.py:882; capture_eff_length.py:18,202-212.
- **Stale doc-path references** — repoint or drop: nonexistent `docs/acc_caljointmodel/prs/PR06_integrate.md`, `docs/caljointmodel/04_interface_contract.md`, `effective_length_redesign_plan.md`; bare PR/section tags ("PR 4b §I.1", "doc 01 §9", "PR 5", "Fix-3", "D2/D3/D5"). priors.py:9-10,190,288,330; efflen/fl docstrings; regions.py/region_arrays.py/signature.py/substrate.py mature-nascent citations.
- **`regions.py` cohesion** (674 lines owns both region + boundary partitions; docstring only lists the region half). Largely self-resolves once the boundary code is deleted.

---

## Index-leverage — the detailed verdict

**Overall: PARTIAL — and the actionable move is REMOVAL, not adoption.** The solver is *not* sloppily re-deriving structure it should read from v6. It correctly reads the canonical per-region `signature` (via `region_arrays.strand_class = transcript_strand_class`), and it deliberately bypasses two annotation-derived columns in favor of stronger observed sources. The defect is the inverse: v6 precomputes a whole tier of columns whose consumers were deleted.

Per-derivation:

| Derivation | Verdict | Why | Behavior-preserving? / Divergence risk |
|---|---|---|---|
| `_node_region_type` coarse exon>intron>intergenic (bp_solver.py:525-534) | **SWITCH** to `signature.coarse_type_array` | Third hand-rolled copy of the same bit test; equivalence asserted in its own docstring; derives from the same `signature` column | Yes, byte-identical (no golden change). No divergence risk — same source column |
| Spliced/mature junction routing `_spliced_faces` / `junction_strand` (bp_solver.py:159-173) | **KEEP** observed accumulator strand; **do NOT** read `genomic_sj_strand` | Motif-observed strand (domain {0,1,2}, one strand/boundary) is correct at AMBIG / exon↔exon seams where flanking signatures cannot orient the junction; the index column is annotation-derived and can be 3=both — not substitutable and strictly weaker | Switching WOULD change behavior. Add a one-line note recording the deliberate divergence |
| TSS/TES "sink" via continuity gate `free_s = has_s(left)&has_s(right)` (bp_solver.py:334-378) | **KEEP** the continuity gate; **do NOT** read `is_tss`/`is_tes` | The gate is strictly richer — it also blocks the non-shared strand at exon↔AMBIG seams; `is_tss/is_tes` would under-block | Switching WOULD change behavior. Add a one-line note |
| `genomic_sj_strand` per-boundary column (regions.py:361-377) | **DELETE** (fold into boundary-flags cleanup) | Redundant, weaker duplicate of `junction_strand`; no runtime consumer | Removal only. Index-schema/golden change (detectable) |
| `is_tss`/`is_tes`/`is_splice_junction` per-boundary columns | **DELETE** | Built, persisted, validated, but read by nothing in production (fed the removed overlay) | Removal only. Index-schema change |
| `mature_eligible_{pos,neg}` per-region columns | **DELETE** | Structural input to the removed mature/nascent message split; no runtime reader | Removal only. Index-schema/golden regen |
| `BoundaryArrays` / `from_boundary_df` | **DELETE** | Zero callers in src/ or tests | Pure dead code, zero production risk |

**Bottom line:** one clean adoption (`coarse_type_array`), two deliberate non-switches to document, and a tier of dead precompute to strip. Do not wire any v6 boundary column into the solver to "justify" it — that re-introduces removed structure.

---

## Module decomposition — is bp_solver.py (1328 lines) too big?

**Yes.** It spans five distinct concerns: (a) per-face geometry + belief/density types + `node_densities` (~160L), (b) statics + strand_stats + init (~200L), (c) the gDNA global-prior/floor/KDE stack (~340L), (d) the disagreement-variance moment (~50L), (e) the FB sweep + the two `chain_*_deconv` projections (~390L). Two independent auditors flagged the same overload, and there is hard evidence it wants to split: `gdna_density_prior.py` reaches back INTO bp_solver for the privates `_node_region_type` and `node_densities`, creating an import cycle.

**Proposed split (do the first; the second is optional):**

1. **Extract the global-prior subsystem → new `gdna_prior_terms.py`** (highest leverage): move `_gdna_seed_estimate`, `_LogLinearVarMean`, `_fit_seed_varmean`, `_floor_estimate`, `_global_logprior`, `_kde_logprior`, `_poisson_moment_var`/`adjacent_disagreement_variance`, and the shared `_node_region_type`. `gdna_density_prior.py` then imports from this module instead of from bp_solver — **the back-import cycle dissolves** and the prior machinery sits next to `GdnaDensityPrior`. Drops bp_solver to ~950 lines.
2. **Optional follow-up — extract the geometry layer → `node_geometry.py`**: `NodeGeometry`/`build_node_geometry`/`node_densities`/`NodeStatics`/`build_node_statics`/`_*_strand_stats`/`init_beliefs`/`_type_belief`. Leaves bp_solver.py as `node_sweep` + the two projections (~450L).

Within `node_sweep` (~330 lines doing five jobs), also **extract the global-prior assembly block (bp_solver.py:1042-1092) into `_build_global_logprior(...)`** so the function reads as the BP skeleton (local solve → scan → combine → final solve → writeback). This move folds naturally into split #1.

Both are behavior-preserving code moves; the only risk is import churn + a golden re-check. Recommend split #1 + the `node_sweep` extraction now; defer #2.

---

## Dead code & stale docs — concrete removal list

**Dead code (post-collapse residue):**
- `BoundaryArrays` + `from_boundary_df` — region_arrays.py:132-193 (zero callers).
- `_compute_mature_eligible` + `REGION_MATURE_COLUMNS` + `RegionArrays.mature_eligible_*` + `_mature_col` — regions.py:245-306,50-52; region_arrays.py:56-59,98-124.
- Boundary structural-flag columns `is_tss/is_tes/is_splice_junction/genomic_sj_strand` (+ possibly the whole `build_boundary_partition`/`load_boundaries`/validate lifecycle) — regions.py:309-383,458-516.
- `statics.spliced_pos/neg` + `_zero_spl` + the `_boundary_strand_stats._floor` closure + region `spl`/spliced_* + the `spliced_pos/spliced_neg` params and `short_p/short_n` terms in the four `simplex_logodds` functions — bp_solver.py:454,1019,1026; simplex_logodds.py:168-194,322-329.
- `node_sweep` params `max_passes`/`convergence_delta` + config fields `sweep_max_passes`/`sweep_convergence_delta` + validators — bp_solver.py:937,949-950,985; config.py; calibrate.py:181-182.
- `NodeDensity.global_density` FIELD ONLY (decl :56-58, assignment :174) — **KEEP** the computation at :160 (consumed at :161 to fill un-anchored regions; deleting it breaks the fill — critic-corrected).
- `StrandBalance.rna_strand_overdispersion` field + computation — strand_balance.py:47,66-71.
- `_SideQuantities` fields `antisense`/`mass`/`mass_spliced`/`count_gdna_frac_var` + the variance block — strand_deconv.py:46-58,124-127,141-145.
- `CalibrationSubstrate.region_len` — substrate.py:88,118.
- `_floor_log_density` + `__all__` entry — simplex_logodds.py:73,30.
- `_MSG_PSEUDOCOUNT` (or wire it in) — bp_solver.py:76-77.
- `is_ambiguous_signature` + `coarse_strand_from_signature` (dead pair), scalar `coarse_type_from_signature` — signature.py:123-149.
- Study branches: `_posterior_median_fg` non-`snap` branches (`RIGEL_MEDIAN_MODE`); `_per_locus_kmeans_rhostar` + `RIGEL_RHOSTAR=kmeans` branch; `RIGEL_ORACLE_CALIB` override; `GdnaDensityPrior.include_boundaries`/`_lscv_bandwidth`/`_find_modes` if the studies are closed — simplex_logodds.py:90-131; priors.py:51-114,414-429; calibrate.py:249-264; gdna_density_prior.py:157-163,270-303.
- `_capture` hook (if no shipped path sets it) — bp_solver.py:1240-1263.

**Stale docs to rewrite/delete:** calibrate.py:1-28,149-154; config.py CalibrationConfig; bp_solver.py:21-24,968-970,979,652,890-895,922; gdna_density_prior.py:15-17; effective_length.py:6-8; fl.py:9-11; simplex.py:1-34; simplex_logodds.py:404-412 + phase labels; priors.py:303-316 (invariant),190-232/285-323/330-336 (triple-doc); dead doc-path refs (priors.py:9-10,190,288,330; regions/region_arrays/signature/substrate mature-nascent + old-PR citations).

---

## Consolidation — duplicated logic and its single home

- **Coarse exon>intron>intergenic type** — 3 copies (`_node_region_type` bit test, `signature.coarse_type_from_signature`, `signature.coarse_type_array`). **Home:** `signature.coarse_type_array`; callers use it.
- **`(2κ−1)²` strand discriminability** — inline at bp_solver.py:586 and :711 + a calibrate.py comment. **Home:** a `strand_discriminability(rna_sense_frac)` helper (natural in the strand cluster, e.g. strand_likelihood.py; or module-local during the prior-subsystem extraction).
- **The per-component message math** (precision `n_src/(n_src·σ²_imp+1)`, density-message projection, precision-weighted log-fraction combine) — 3 near-identical copies (gDNA, RNA+, RNA−). **Home:** one `_component_message(...)` helper in bp_solver's scan.
- **Clean intron/intergenic↔exon boundary derivation** — `gdna_density_prior._clean_exon_boundary` duplicates `_gdna_seed_estimate`'s. **Home:** one shared `clean_exon_boundary(...)` in the geometry/prior module.
- **`beta_concentration`** — byte-identical on both strand models. **Home:** one free function beside `overdispersion_for_beta`.
- **`strand_loglik` vs `_mixture_strand_loglik`** — two Beta-Binomial-moment implementations. **Home:** `_mixture_strand_loglik` (single-strand specialization re-exported for tests) — only if touching the area.
- **Region+pooled-seam gDNA node-model prose** — documented 3×. **Home:** `_gdna_region_node_arrays` docstring; others get a one-liner + pointer.

---

## What we are NOT changing — scope guard

- **Do not re-introduce removed complexity.** The collapse to one belief-free Poisson disagreement-variance precision path was deliberate and validated. Do not resurrect the mature/nascent 5-term overlay, per-component A/B variance functions, polymorphic dispatch, the per-pass var~mean-refit outer loop, or the spliced-sink/Q2-floor experiment. All doc rewrites must describe the *single* FB path, not restore the old one.
- **House style stays.** Dense semicolon one-liners are the codebase idiom — do not flag `E702` or "simplify" them. Conventional short names (`kappa` in comments/where conventional) are fine; only genuinely obscure Greek-letter *variable names* warrant a rename.
- **Keep the deliberate index-leverage non-switches.** The observed `junction_strand` (not `genomic_sj_strand`) and the `free_s` continuity gate (not `is_tss/is_tes`) are load-bearing correctness choices — document them, do not "consolidate" them against v6 columns.
- **Don't over-refactor the healthy corners.** `gdna_strand.py`'s shared `_fit_overdispersion` core, `strand_likelihood.py`, `effective_length.py` as single FL-marginal owner, and the structure-inputs live paths are exemplary — leave them.
- **Magic-number items are discussion flags, not silent edits.** Per project ethos, surface `_MIN_KDE_TEACHERS`, the KDE bridge percentiles, and the capture shrinkage constant for the user; do not just hardcode-around them.
- **Index-schema removals are gated changes.** Deleting feather columns (`mature_eligible_*`, boundary flags) requires a schema bump + golden index regen; sequence and land these deliberately, not as drive-by deletes.