# Unified solver — build handoff & bug-fix continuation (2026-07-23)

**Branch `calib-ambig-init-wip`. Nothing committed — the working tree holds it all. DO NOT REVERT: we fix
forward.** The one-mode "unified" solver is BUILT, flag-gated (byte-identical off), and its **message
arithmetic is now correct** (§4 — the composition invariant is restored). What remains is a *precision*
problem, not an arithmetic one, and it is the task the design already names R2/R3.

Read `docs/calibration/unified_solver_design.md` first (the full design; §2 the theorem, §6 the mature
routing, §10 the exact implementation). This doc = state + what was fixed + where to go next.

---

## 1. What is built, and where

| piece | location | status |
|---|---|---|
| **`_scan_unified`** — the one-mode relay + two-iteration combine | `bp_solver.py`, nested in `node_sweep` (branched at `if _UNIFIED:`) | runs; **byte-identical off** |
| `node_total_density` — lazy, composition-aware `ρ_tot` (per-side spliced) | `bp_solver.py` after `node_global_geometry` | tested |
| `reframe_density`, `density_mode_logfrac`, `total_density`, `k_*` | `calibration/enrichment_frame.py` | tested (52) |
| the §2 theorem + primitives | `tests/calibration/test_enrichment_frame.py` | 52 pass |
| **`unified_message_audit.py`** — fused message vs TRUTH per node class | `scripts/debug/` | the primary probe |
| **`unified_relay_trace.py`** — per-hop drift localization | `scripts/debug/` | how §4 was found |

**Flags** (all default to the current path): `RIGEL_UNIFIED` (on = the new solver), `RIGEL_UNIFIED_ROUTE`
(mature peel/graft, default 1), `RIGEL_UNIFIED_RHO_ITERS` (lazy-`ρ_tot` iterations, default 2),
**`RIGEL_UNIFIED_S2T`** (new — the R3 σ²_transfer ablation: `1` keep / `flat` drop the cliff-height term /
`0` retire entirely).

**Standing green gates:** flag-off `gdna_none` guard = **3,704,635** (byte-identical); `pytest
tests/calibration/ tests/native/` = **371 pass, 1 xfail, 2 xpass**; `ruff check src/ tests/` clean.
*(Pre-existing and NOT from this work: 5 `tests/test_golden_output.py` failures and 3 `ruff` findings in
`scripts/debug/{confident_fp_trace,enrichment_relay_validation}.py` — both reproduce at `HEAD`. Goldens are
regenerated last, per the build plan.)*

---

## 2. The scoreboard

`pass0_oracle_bench.py`, 32 scenarios, mass-weighted:

| arm | mwae ↓ | corr ↑ | `gdna_none` guard ↓ |
|---|---|---|---|
| base (the shipping `_scan`) | **0.1361** | **0.688** | **3,704,635** |
| unified — at the previous handoff | 0.1865 | 0.524 | 4,352,156 |
| unified — now (default `S2T=1`) | 0.1864 | **0.617** | 5,456,663 |
| unified — now, `RIGEL_UNIFIED_S2T=flat` (⚠ breaks R4, §5) | **0.1573** | 0.597 | 4,075,777 |

Per node class on `gdna_gdna300_ss_0.50_nrna_present_capture_on` (`unified_node_type_diag.py`), mean `f_g`:

| | intron | exon | boundary |
|---|---|---|---|
| oracle | 0.848 | 0.677 | 0.646 |
| unified, previous handoff | 0.309 | 0.379 | 0.303 |
| unified, now (`S2T=1`) | 0.595 | 0.431 | 0.567 |
| unified, now (`S2T=flat`) | 0.630 | 0.497 | 0.572 |

The "washes to ~0.3 regardless of truth" symptom is **gone**. Correlation is up 0.524 → 0.617 and the
mwae gap to base has closed ~58 % (0.1865 → 0.1573). It is not yet a win, so **do not flip the default.**

---

## 3. ⭐ THE ROOT CAUSE — the relay was not carrying a composition

**The message factor per component is `f_c = ρ_c·E_c/M_dst`, and `Σ_c f_c` MUST be ≈1 — it *is* a
composition, which is exactly what design §2 proves the reframe buys.** It was not. Measured:

```
node class    Σ_c f_c  (must be 1)
intron            75.6
boundary          42.6
exon               8.5
```

The own (message-free) seed was a clean 1.00–1.01 at every class, so the defect was entirely in the relay.

**Why.** `_unified_solve` relayed three *independent* per-component densities with nothing tying them to
the node's own observed mass. The reframe `r = ρ_tot(dst)/ρ_tot(src)` cannot supply that tie: `ρ_tot` is
built from the node's *belief* `f_g` while the message carries its own composition, so every hop multiplies
in an uncancelled residual. Geometrically the residual telescopes (`exon→bnd` 0.431 × `bnd→exon` 2.298 ≈ 1),
but its per-hop spread is enormous — p90 = **209** on `bnd→exon`, because `ρ_tot` itself swings ~500×
between an intron and its neighbouring exon. So the relay is a **multiplicative random walk**, and over the
measured mean adopt-run of 9 hops (max 31) its arithmetic mean ratchets to 6–20×.

**This is a known, previously-retired defect class.** `_scan`'s docstring says it outright: *"The old
three-independent-log-fraction relay violated all three (a boundary relayed `fbp = 51.9`)"* — which is why
the shipping path stores the running belief as `(μ_λ, σ²_λ)`, a composition by construction. The unified
solver reintroduced the old shape and reproduced the identical failure at the identical magnitude.

**The fix (`_pin` / `_pin_v`)** is the design's own ÷`M_dst`, applied at **every** node rather than only at
the final combine: a belief about node `i` is a claim about `i`'s own observed mass, so scale it to that
mass. An **uninformed** component (precision 0) is supplied by `i`'s own density in the normalizer — that
is what keeps a partial message partial (a seam sending gDNA only still gives `f_g < 1`, §2) instead of
renormalizing it into the shift's "the missing component is absent". A structurally-dead strand has own
density 0 and contributes nothing, correctly. Result: **Σ_c f_c 75.6 → 1.02.**

---

## 4. The other five fixes (all measured, all in the tree)

1. **`_relay` reframe denominator.** `src_face` was a **dead parameter**; the code divided by
   `rho_node0[s]` (the mature-free NODE `ρ_tot`) instead. Since the previous hop scaled the density *into*
   `dst_face[s]`, every acceptor boundary injected an uncancelled `(1+ρ_spl/ρ_unspl)` that compounded.
   Relay gDNA inflation 15× → 7×.
2. **`_transport` frameless-source fallback.** `rho_src = 1.0` made `r` the destination's **absolute**
   density (10³ on a short node) — a raw scale masquerading as a ratio. Now `r = 1` pass-through, matching
   the relay's own guard and §5.
3. **The mature routing is per-FACE and exon-only** (§6). The graft fired into `is_reg_a[i]` — *any*
   region, so a junction's whole mature flux was grafted into every flanking **intron**, which carries no
   mature (intron `Σf_r` 70.9 → 15.8 with the route off). It now fires only into `ex_a[i]`, and both peel
   and graft use the face that faces the other endpoint (`spl_*_f[df]` / `[sf]`) rather than the
   node-pooled sum, which double-counted at exon↔exon junctions — the A3 per-face fix `_scan` already
   carries in `mature_dilution[df]`. The graft is also applied **before** the reframe (it is a density
   measured at the source); the peel stays after (measured at the destination).
4. **Coherent own precision on both arms.** `_own_pg` is built from the REFERENCE-FREE evidence
   `tau0_lam`; `_own_pp`/`_own_pn` were built from `vp_loc`/`vn_loc` — the local solve's **posterior**
   variance, which pools the shared reference (the phantom `message_precision_derivation.md` §2 forbids as
   a message precision). On unstranded data that asymmetry is total: τ=0 ⇒ the gDNA arm is silent at every
   non-locked node (measured `p_own_g` = 0.00) while the RNA arm still asserts "f_r = ½" at precision 0.35.
   Because `_pin` normalizes the pie, that one-sided claim pushed `f_g` down everywhere. Rebuilt from the
   same τ with `_scan`'s Jacobian (`Var(log f_r) = f_g²/τ_λ`). **The single largest modelling gain**
   (intron 0.415 → 0.540, exon 0.400 → 0.532).
5. **The graft's MEASUREMENT precision is not τ-gated.** Fix 4 alone drove every RNA message precision to
   0 — i.e. it re-created the **τ-gag bug** (a junction COUNT is not an imputation; gating it on τ silenced
   52 % of spliced junctions unstranded). The prediction channel stays τ-gated; the graft now carries
   `S/(1+S·σ²_T)` on top, exactly as `_scan` does.

---

## 5. R3 — measured, and it is real

Design §8 R3 predicted that once the mode is enrichment-correct, σ²_transfer's cliff-height damping
`(μ_dst−μ_src)²` becomes the double-count of the reframe and should be retired *with* this change. That is
now measured (`RIGEL_UNIFIED_S2T`):

| | mwae | corr |
|---|---|---|
| `1` (keep — default) | 0.1864 | 0.617 |
| **`flat`** (drop the cliff term, keep the `var_proj` floor) | **0.1573** | 0.597 |
| `0` (retire entirely) | worse per-class; not benched |

**`flat` is a −0.029 mwae move, 19 scenarios better / 11 worse**, and it pulls the phantom guard back from
5,456,663 to 4,075,777. It helps regions (intron |err| 0.327→0.306, exon 0.288→0.261) and hurts boundaries
(0.319→0.346).

**⚠ But `flat` VIOLATES R4 and is NOT shippable.** Stranded capture-ON regresses badly:

| stranded scenario | base | `S2T=1` | `S2T=flat` |
|---|---|---|---|
| `gdna300_ss_0.99_nrna_present_capture_on` mwae | 0.0283 | 0.0593 | **0.1191** |
| — corr | 0.986 | 0.952 | **0.907** |
| `gdna100_ss_0.99_nrna_present_capture_on` mwae | 0.0717 | 0.0793 | **0.1402** |
| — corr | 0.946 | 0.945 | **0.842** |

The default arm holds stranded corr ≥ 0.945 across the board; `flat` drops it to 0.842–0.907. This is
exactly the risk the previous handoff flagged — *"the σ²_transfer damping is load-bearing until R2/R3; do
not remove it before the transfer-variance rework."* So the aggregate mwae win is bought from the stranded
arm, and it is **left off by default**: `flat` is the *starting point* for R2/R3, not a fix. `Var(log r)`
(R2) must land **with** it, since retiring the damping without supplying the reframe's own variance is only
half the trade — the messages then arrive undamped AND unqualified, and a strong strand likelihood gets
diluted (R4 by construction).

---

## 6. NEXT STEP

The arithmetic is sound and the theorem holds; the residual is a **systematic compression toward the
middle** — `f_g` tracks the truth (corr 0.617) but is pulled to ~0.5 from both directions. That shape is a
precision problem, so the next task is the one the design already scopes:

1. **R2 — `Var(log r)`.** The reframe is applied as exact. `enrichment_frame.composition_logvar` already
   derives it; feed it into the fused precision. This is the missing half of §5's trade.
2. **R3 — retire σ²_transfer** together with (1). `RIGEL_UNIFIED_S2T=flat` is the measured starting point;
   understand the boundary regression it causes before flipping.
3. **The exon RNA over-call is the NASCENT channel, not the graft.** ✅ *Measured and closed —*
   `scripts/debug/unified_graft_check.py` compares the grafted mature mass against the oracle's contained
   `mat_uns_*` at every exon: **median 1.070, geomean 0.898, mass-weighted mean 1.195**; as a fraction of
   the exon's unspliced mass, grafted `f_mature` = 0.180 vs oracle 0.222. So the graft is **unbiased and
   if anything slightly low** — `E_spl` and `E_r` ARE commensurate, and the standing spliced-frame warning
   at `bp_solver.py` ~L1069 does **not** bite here. (Per-node spread is wide, p10 0.24 / p90 2.59, so it is
   noisy, not biased.) Since exon `msg_f_r` is 0.550 against a truth of 0.323 while the mature part of that
   is right, **the over-call is in the relayed NASCENT (unspliced-imputed) RNA** — look at what the peel
   leaves behind on the `exon → boundary` direction and at how much RNA arrives from non-junction
   neighbours, not at the spliced frame.
4. **The anchor structure.** On unstranded data the *only* gDNA composition anchor is the intergenic
   `struct_lock` (precision 71071 vs 0 everywhere else), so one node's claim is smeared genome-wide with no
   restoring force. `_scan` bounds this with the `_gdna_wall` (intergenic SINKS incoming gDNA, §12.9);
   `_unified_solve` has **no wall**. Worth an A/B.

---

## 7. Commands

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
cd /Users/mkiyer/proj/rigel

# byte-identity + phantom guard (base 3,704,635)
OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py
RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py

# THE PRIMARY PROBE — fused message vs truth; watch the SUMf column (must be ~1)
RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_message_audit.py
# per-hop drift localization (which edge type injects an uncancelled factor)
RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_relay_trace.py

# the battery (base arm already in /tmp/pass0_oracle_bench.tsv; delete to reset)
RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --arm uni
OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --report --vs base --new uni

# per-node-type localization; ablations
RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py uni
RIGEL_UNIFIED=1 RIGEL_UNIFIED_S2T=flat OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py flat
RIGEL_UNIFIED=1 RIGEL_UNIFIED_ROUTE=0 OMP_NUM_THREADS=1 python scripts/debug/unified_node_type_diag.py nr

OMP_NUM_THREADS=1 python -m pytest tests/calibration/ tests/native/ -q   # 371 pass
```

**Gates on every fix:** flag-off byte-identical (3,704,635); **Σ_c f_c ≈ 1** in the message audit (the new
cheap invariant — check it *first*, it catches an arithmetic regression in one second); stranded must not
regress (R4); battery per condition, stranded AND unstranded; goldens last.

---

## 8. Where the code lives

`_unified_solve()` is nested in `node_sweep` (`bp_solver.py`). Inside it:
* `_rho_faces(fgc)` — lazy `ρ_tot` (node / left-face / right-face, WITH spliced at the acceptor).
* `_pin` / `_pin_v` — **the composition pin** (§3). The invariant lives here.
* `_relay(seq, nbr, dst_face, src_face, df, sf)` — forward/backward context accumulation; `_damp`
  (σ²_transfer) + the per-face mature routing.
* `_transport(src, valid, df, sf, arrs, dst_face_v, src_face_v)` — the vectorized combine reframe + graft +
  peel + filter + pin.
* `_fuse_v` — precision-weighted density fuse (α ⊗ β).
* the two-iteration loop → `density_mode_logfrac` (inline as `mo_g/mo_p/mo_n`) → `_local_solve` (the ψ solve).

An inert `_capture` hook (`_capture["_uni"]`, `_capture["_uni_static"]`) dumps every intermediate the two
probes read; it is `None` in production.

**Design record:** `unified_solver_design.md` §2 (theorem, ÷M_dst), §6 (mature routing), §10 (the precise
implementation), §8 (R1–R4 open risks — R3 now measured, see §5 above).
