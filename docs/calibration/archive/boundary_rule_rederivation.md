# The boundary rule, re-derived — the "inversion" was a Simpson metric artifact

**Status (2026-07-22): DERIVE step complete — the premise is REFUTED.** The task was to "re-derive the inverted
boundary solvability rule" (`solve_gate_design.md` §4: DOF calls coin-flip boundaries solvable and meaningful
boundaries unsolvable). Careful dissection shows there is **no inversion to fix**: the apparent inversion was a
**between-structure (Simpson) confound** in the pooled correlation metric. Under an honest per-position metric,
**boundaries resolve iff STRANDED (corr 0.85) and are a pure identifiability floor when UNSTRANDED (corr 0.04)** —
which is exactly what the strand evidence `tau0_lam` (already in the DOF criterion) reflects. The owner's
enrichment-cliff mental model is a real *mechanism* but is **not** what discriminates resolving from non-resolving
boundaries — strand is.

Tools: `scripts/debug/boundary_dissect.py` (the boundary cuts), `scripts/debug/pos_corr_by_type.py` (the honest
per-position metric across all node types + the stranded control). Both pool the cached `ambig_dense_10mb`
substrate (`_selfsolve_cache`), refit=0, `OMP_NUM_THREADS=1`.

---

## 1. Reproduce the "inversion" (pooled metric — the previous session's basis)

Pooled corr(f_g, oracle) over the 20 **unstranded** (ss0.50) scenarios, split by the DOF verdict:

| node type | DOF | pooled corr | previous reading |
|---|---|---|---|
| single-strand region | SOLVABLE | 0.633 | resolves ✓ |
| AMBIG region | SOLVABLE | 0.689 | resolves ✓ |
| boundary | SOLVABLE | **0.129** | "coin-flip — wrongly solved" |
| boundary | UNSOLVABLE | **0.675** | "meaningful — wrongly withheld" |

This is what `solve_gate_design.md` §3–4 called the boundary inversion.

## 2. The cliff hypothesis FAILS to separate (owner's model)

Dissecting the unstranded boundaries by the owner's proposed discriminators — none separates resolving from
coin-flip (`boundary_dissect.py`, pooled corr):

- **CUT 3 — flank enrichment gap** (the cliff magnitude `|mu_proj[L]−mu_proj[R]|`): quartiles 0.32 / 0.45 / 0.38 /
  0.45. A **bigger cliff is not worse.** The cliff does not discriminate.
- **CUT 5 — flank gDNA-message disagreement** `|log f_g^L − log f_g^R|`: agree 0.463 vs disagree 0.385 — barely
  differs. "Discordant messages muddle it" is not supported.
- **CUT 6 — one-sided spliced present**: has 0.054 vs no 0.514 — the *cleanest* pooled cut, but see §3: it is the
  same Simpson artifact (spliced boundaries sit at different structural positions than no-spliced ones).

## 3. The confound — the pooled corr is Simpson-contaminated

**CUT 7** — split the DOF verdict *within* each flank-type pair:

| flank pair | DOF-solvable corr | DOF-unsolvable corr |
|---|---|---|
| intron\|exon | 0.116 | −0.031 |
| exon\|exon | 0.060 | −0.059 |
| intgc\|exon | — | (oracle f_g≡1, corr undefined) |

*Within a fixed flank-type, NEITHER DOF group is meaningful (all ≤ 0.12).* So the pooled "boundary UNSOLVABLE
0.675" is **entirely between-structure**: DOF-unsolvable boundaries happen to concentrate at
intergenic-adjacent seams where the coarse structural signal is f_g≈1, and that coarse signal — not any
per-boundary deconvolution — drives the pooled correlation up. Classic Simpson's paradox.

## 4. The honest metric — per-position correlation across scenarios

Hold a node **fixed** (same chain id, stable across scenarios) and vary the gDNA level across the 20 unstranded
scenarios; corr_i = corr(f_g, oracle) over that node's scenarios. This removes the between-structure confound.
Aggregate corr_i by type × DOF (`pos_corr_by_type.py`):

| node type | DOF | **per-position corr_i** | pooled (for contrast) |
|---|---|---|---|
| single-strand region | SOLVABLE | **0.434 RESOLVES** | 0.633 |
| single-strand region | UNSOLVABLE | 0.153 coin-flip | −0.037 |
| AMBIG region | SOLVABLE | **0.503 RESOLVES** | 0.689 |
| AMBIG region | UNSOLVABLE | 0.064 coin-flip | −0.066 |
| **boundary** | **SOLVABLE** | **0.039 coin-flip** | 0.129 |
| **boundary** | **UNSOLVABLE** | **0.007 coin-flip** | 0.675 ← the artifact |

**The metric is validated by its own contrasts:** it gives high corr for the nodes that *should* resolve
(regions the DOF criterion solves) and ~0 for those that should not (the withheld regions), and it exposes the
0.675 boundary artifact as ~0. Mass-weighting the boundary corr makes it *worse* (−0.05), so this is not a
low-count-noise artifact — high-mass unstranded boundaries do not resolve either.

## 5. The stranded control — boundaries DO resolve when stranded

Same metric, same code, on the 12 **stranded** (ss0.99) scenarios (`pos_corr_by_type.py 0.99`):

| node type | DOF | per-position corr_i |
|---|---|---|
| single-strand region | SOLVABLE | 0.966 |
| AMBIG region | SOLVABLE | 0.844 |
| **boundary** | **SOLVABLE** | **0.845 RESOLVES** (92 % of positions corr>0.5) |

So a boundary's identifiability is governed by **strand** (`tau0_lam>0` ⇒ resolves; `≈0` ⇒ floor), exactly like a
region. The DOF criterion already keys on `tau0_lam`. **There is no boundary-specific inversion.**

---

## 6. The reframe (part 1) — no inversion; the DOF "inversion" was a Simpson artifact

1. **No inversion to fix.** The previous session's "boundary rule inverted (a real bug)" was the third
   metric-artifact of this arc (after mwae-counts-defaults and the emission-golden): a pooled correlation
   confounded by coarse structure. Corrected here and in `solve_gate_design.md`.
2. The honest per-position metric shows unstranded boundaries do not currently resolve (corr 0.04) — but that is
   *not* an identifiability floor. §7 shows the information IS there and the arithmetic destroys it.

## 7. The real finding — the information IS received; the message arithmetic SATURATES it (owner's thesis, confirmed)

Owner's thesis (2026-07-22): "boundaries receive the information they need to solve, but the solution is bad …
unreliable, but not necessarily unsolvable." An information-chain decomposition (`boundary_info_chain.py`,
per-position across the 20 unstranded scenarios, all corr vs the boundary's ORACLE f_g) confirms it exactly:

| stage | corr | reading |
|---|---|---|
| **CEILING** corr(bnd_oracle, flank_oracle) | **0.886** | the answer IS in the neighbour ✓ |
| SRC-QUAL corr(flank_oracle, flank_solved) | 0.352 | the neighbour only partly solves itself (unstranded region) |
| SRC→BND corr(bnd_oracle, flank_solved) | 0.224 | the transferable signal is positive |
| **TRANSFER1** corr(flank_solved, flank_message) | **0.212** | the message barely carries the source's own belief |
| DELIVERED corr(bnd_oracle, combined_message) | −0.071 | ≈0 by the time it arrives |
| SOLVED corr(bnd_oracle, bnd_solved) | 0.040 | |

**TRANSFER1 = 0.21 is the bug.** The message is a *deterministic function* of the source belief (the mode
formula), so `corr(source_solved_fg, message_implied_fg)` should be ≈ 1.0. It is 0.21 — the mode arithmetic
scrambles the source's gDNA fraction.

**The mechanism (trace, `boundary_mode_trace.py` node 3236):** the boundary's gDNA message **saturates to
f_g ≈ 1.0 (all-gDNA over-call)** in nearly every scenario, decoupled from the source's belief (~0.4). Two ways,
both routing to the same clamp:
- **Density mode** `implied = clip(max(ρ_g, 1/E_g^dst)·E_g^dst / md)`: divides the source-imputed gDNA by the
  boundary's **sparse, noisy crossing mass `md`** (often 0–3 fragments). When `ρ_g` is small it floors at the
  **one-fragment floor** (`imputed = max(ρ_g·E_g, 1) = 1`), and `1/md` with `md<1` → clamps to 1. When `ρ_g` is
  larger it **over-imputes** (`ρ_g·E_g^dst` up to 906 vs `md`=297 → ratio 3 → clamp 1).
- **Shift mode** (`RIGEL_GATE_SHIFT=all`): *worse* (message corr −0.206). Here `f_g = imputed_gDNA / (imputed_gDNA
  + imputed_RNA)`, and the imputed RNA **collapses** because the mature subtraction (`rho_pos/neg − absorb`)
  drives `ρ_r → 0` → `f_g = gDNA/(gDNA+0) → 1`. Same saturation.

**Common root:** at a boundary the RNA is systematically **under-imputed** / the gDNA **over-imputed**, saturating
`f_g → 1`. The density mode does it via the sparse-`md` division + one-fragment floor; the shift does it via the
mature-subtraction RNA collapse. Both feed the boundary a biased, near-constant all-gDNA message that carries
none of the source's real gDNA belief. The fold then discounts it (solved corr 0.04, not 1) — which is why the
boundary is *unreliable* rather than *confidently-wrong-at-1*.

## 8. Fix direction (to design next — NOT yet implemented)

The theory (`message_propagation_arithmetic.md` §3) says the recipient rule is the SHIFT — `f_g^dst =
ρ_g E_g^dst / (ρ_g E_g^dst + (ρ_+ + ρ_−) E_r^dst)` — normalize by the IMPUTED total, not the observed `md`. The
shift is the right frame; it fails today because the **imputed RNA at a boundary is wrong** (the mature
reconciliation collapses it) and the **sparse-crossing statistics** are unstable. So the boundary fix is really
the **RNA-imputation / mature-reconciliation at boundaries** (`pass0_roadmap.md` §4-B — the mature capture-scale
correction), plus a sparse-`md` guard. Concretely to derive next:
1. Why does the imputed RNA collapse at a boundary (the `−absorb` mature subtraction; the `rho_pos/neg` sign)?
2. Can the shift be made to carry the source's f_g faithfully (TRANSFER1 → ~1) once the RNA imputation is fixed?
3. The one-fragment floor `1/E_g^dst` with a large boundary `E_g^dst` is a *large absolute* gDNA floor — is that
   the correct resolution floor for a sparse crossing, or should the floor scale with the crossing opportunity?

This is a genuine pass-0 **arithmetic** task (owner's priority), not a solve-gate and not the hyperprior.

## 9. The RNA-imputation bug, oracle-confirmed (`boundary_rna_trace.py`)

Owner's model (2026-07-22): the exon→boundary hop is the reliable data; the boundary "absorbs the spliced density
and must remove it … the remaining unspliced fragments → gDNA + nascent." The trace of an intron|exon boundary
(node 3236, exon-flank source) confirms the model and pinpoints the break, using the oracle's nascent/mature
split (`nas_uns` / `mat_uns`):

- **`or_mat` (oracle mature *unspliced* at the boundary crossing) = 0 in EVERY scenario.** Mature never straddles
  the junction unspliced — it is entirely in the spliced channel. So the boundary crossing is *physically* just
  gDNA + nascent, exactly the owner's model.
- **The current split under-removes mature.** The code does `rho_r = rna_src − mat_abs` = (exon RNA density) −
  (junction-spliced density). At an enriched exon: `rna_src`=3.28 but `mat_abs`=0.21 → residual "nascent" = 3.07,
  while the boundary's true nascent ≈ 0 (oracle f_g=0.99). Shift-implied f_g = **0.53 vs oracle 0.99**.
- **Root cause:** the exon's `rna_src` is dominated by **within-exon mature** (mature mRNA inside a single exon has
  no visible junction → counted *unspliced* → it sits in `rna_src`, not in the junction `mat_abs`). The junction
  spliced count captures only the mature that crosses the junction — a small fraction — so the subtraction leaves
  ~all the mature mislabeled as nascent. And calibration *cannot* separate within-exon mature from nascent at the
  exon (that is the per-locus EM's job) — so importing the exon's "RNA" density inherently drags its mature along.

**⟹ the reframe (owner's model, made precise):** do NOT import the exon's mature-contaminated RNA to impute the
boundary's nascent. The boundary's OWN unspliced crossing (`md`) is already mature-free (oracle-confirmed), so
split *that* into gDNA + nascent: take the gDNA density from the reliable exon flank, and the **nascent as the
residual** `md/E_r − ρ_g`. No nascent/mature separation at the exon is needed.

## 10. The gDNA transfer + σ²_transfer — the cliff piece (`boundary_transfer_var.py`)

The reframe still needs the exon's gDNA density transferred to the boundary's crossing frame — and the exon
gDNA is *enriched* while the boundary crossing gDNA is only *partially* enriched (it straddles). That is why
`ρ_g·E_g^dst` over-imputes (906 vs md=297 → saturates). σ²_transfer is supposed to damp cliff-crossing messages,
but at intron|exon boundaries it does the opposite of what's wanted:

| capture | exon→bnd σ²_T (median) | intron→bnd σ²_T | intron flank mass |
|---|---|---|---|
| off | 3.93 | **0.13** | 763 |
| on | 3.47 | **0.39** | 47 |
| verystrong | 6.37 | 8.76 | **1.0** |

- **σ²_transfer damps the reliable exon edge, trusts the sparse intron edge.** The boundary's own projected
  enrichment `μ_proj[bnd]` usually sits near the depleted intron, so the mode-gap `(μ_bnd − μ_exon)²` is large →
  the dense, reliable exon message is gagged while the sparse intron is trusted. Backwards for the boundary.
- **The intron degrades under capture, exactly as the owner predicted:** intron flank mass 763 → 47 → 1 as
  capture strengthens; under verystrong capture the intron is empty and only *then* its σ²_T blows up.

## 11. The proposed architecture (to design → plan → execute, owner-gated)

Two separable pieces:
1. **RNA/nascent = residual of the boundary's own mature-free crossing.** Stop importing the exon RNA density;
   `ρ_nascent = md/E_r − ρ_g`. Removes the mature-contamination bug (§9).
2. **gDNA density transfer across the partial-capture cliff.** Transfer the exon's (enriched) gDNA density to the
   boundary's crossing frame with a capture-scale correction, and fix σ²_transfer so it trusts the dense exon
   rather than gagging it (§10). This is the hard, owner-intuition piece (the cliff).

## 14. Fix #2 LANDED — the GEO-MEAN crossing mode (derive-boundary-mode workflow)

The transfer mode (§11.1 piece 2) was derived by a 5-hypothesis multi-agent workflow over the cached
`ambig_dense_10mb` **panel** (`scripts/debug/boundary_panel.py` — 23 high-error splice-junction exon→boundary→
intron triplets from the worst scenarios, tracked across all 32; evaluator `boundary_mode_lib.py`).

**The answer — the geometric-mean cliff-interpolated, mature-free crossing decomposition.** Capture enrichment is
multiplicative, so a partial-capture crossing sits at the **log-midpoint** of its two flank gDNA densities:
`ρ_g^cross = √(ρ_g^exon · ρ_g^intron)` (single-flank fallback). The crossing is mature-free, so `f_g =
ρ_g^cross·E_g/md` and **nascent = residual**. **ZERO free constants; count-legal** (per-node source densities;
`md` only as the mature-free total being split). The workflow REJECTED the higher-scoring candidates
(own-crossing-anchor, exon-gdna-guarded) as **count-zero-info violations** — they let the crossing *mass* `md`
vote composition (`corr(md,f_g)=−0.25`, verified). The geo-mean sits **on the identifiability frontier**
(`under_fix+over_guard ≈ 1.0`, a genuine wall); the residual (per-boundary composition + the low-gDNA tail) is
the **global gDNA prior's** job, not the mode's.

**Implementation** (`bp_solver`, three coupled changes on boundary-dst edges, now default):
1. **Suppress** the mature-contaminated exon unspliced-RNA import (`fp_s·sm/_er`) → nascent = residual. Keep the
   spliced MEASUREMENT + the mature absorption. The intron flank (mature-free) keeps its clean nascent.
2. **Geo-mean** crossing gDNA `rho_g_cross` (precomputed, vectorized) as the gDNA density.
3. **Density frame** (observed-`md` mature-free anchor).

**Measured (Fix#1 + geo-mode):** the antisense-intronic leak (§13.1) FIXED (t2 93 → **0**, phantom locus 61 → 6);
boundary per-position corr 0.105 → **0.233**; regions +; off = byte-identical. Offline panel score 0.19 → 0.67
under-fix, panel_err 0.81 → 0.33.

**Interrogation 1 (change 1 is load-bearing).** The geo-mean gDNA (2+3) WITHOUT change 1 regresses (boundary corr
−0.04, phantom up, leak persists): the geo-mean gDNA fights the still-mature-contaminated exon RNA message.
Change 1 supersedes the mature-gate-dismantle's exon→intron *direct* nascent emission — but the net nascent is
preserved (it now flows from the **mature-free intron** flank via the boundary's residual relay; the
intron→exon relay + mature-no-hallucination tests still pass; only the one *mechanism* test was updated →
`test_geomode_suppresses_exon_unspliced_rna_into_boundary`).

**Interrogation 2 (the quintuple grid — DEFERRED to a fresh thread).** The owner's `intron ↔ IE ↔ exon ↔ EI ↔
intron` grid (`scripts/debug/quintuple_grid.py`, the exon is a RELAY, not a source — its gDNA is imputed from the
introns, the 3-hop telephone) exposes a **separate** failure: when **nascent ≫ mature**, the solver **over-calls
gDNA** (middle-exon f_g 0.81 vs oracle 0.31; 0.69 vs 0.05 at near-zero gDNA). The geo-mode gives *identical*
results there — this is the **unstranded nascent-vs-gDNA identifiability floor** (nascent is unspliced ⇒ locally
indistinguishable from gDNA; the only RNA proof is the scarce spliced/mature), NOT the boundary mode. The lever is
the **intron→exon nascent relay** (the introns are transcribed and well-defined) or the global gDNA prior. This
is the next investigation; the grid is its regression harness.

## 12. σ²_transfer bugs (owner-flagged, 2026-07-22) — both confirmed in the code

`node_global_geometry` ([bp_solver.py:174](../../src/rigel/calibration/bp_solver.py#L174)) builds the total
density the enrichment NPMLE / σ²_transfer is fit + projected on. `mass = msl + msr` (UNSPLICED crossing only),
`eff = egl + egr` (gDNA eff-length). Two owner concerns, both confirmed:

1. **Spliced density is EXCLUDED.** The boundary's density counts only its unspliced crossing (gDNA + nascent),
   not its spliced (mature). The exon flank's density, by contrast, counts its contained unspliced mass — which
   *includes* within-exon mature (unspliced, no junction). So the two densities are on **asymmetric bases**: the
   exon's includes mature, the boundary's does not. This inflates the mode-gap `(μ_exon − μ_bnd)²` → **σ²_transfer
   gags the exon→boundary edge** (the backwards behavior in §10). Owner's fix: include the boundary spliced
   density so both totals are on the same (mature-inclusive) basis. **Directly connects the §10 backwards-σ²_T to
   the spliced-exclusion.**
2. **The total density is composition-CONDITIONAL (the deeper issue).** `eff = eff_gdna` makes the density a
   gDNA-frame quantity `mass/E_gdna`, but the mass is a gDNA+RNA mixture and gDNA/RNA have **different fragment
   lengths → different eff-length denominators**. So the "total density" the enrichment landscape is built on
   *depends on the very composition we are solving for*. Modeling a transfer variance on a composition-conditional
   density is ill-posed at a boundary (where composition is the unknown). This is a foundational concern, not a
   tweak.

**Owner's further point (no single capture scheme):** probe placement varies — sparse probes (many exons
uncaptured), or a probe *over the splice junction* (the boundary can enrich HIGHER than the exon). So the
exon→boundary capture-scale is not monotone and cannot be a single formula; every junction may differ by
proximity to the nearest probe. ⟹ the derivation needs a **diverse boundary substrate** (probe-placement regimes
× gDNA × expression, with oracle) and must be proven **one step at a time**, possibly via a multi-agent
derivation workflow over many worked examples.

## 13. Fix #1 LANDED (measured) — the DIRECTIONAL spliced-density σ²_transfer (`RIGEL_BND_SPLICED_DENSITY`)

The first σ²_transfer bug (§12.1) is fixed and measured. Include the boundary's spliced (mature) **density**
`ρ_spl = Σ_side spl/E_spl` (per-frame, not raw mass — avoids the §12.2 composition-conditional trap) in the
boundary's total **only on the exon↔boundary edge** (owner's correction: the intron has no mature, so
intron↔boundary stays mature-free — hence per-EDGE, not per-node). Code: `_boundary_spliced_mass_increment` +
a second `mu_proj_mat`/`var_proj_mat` projection selected per-edge in `_scan` (env-gated, off = byte-identical).

**Effect (per-position corr across the 20 unstranded scenarios; `pos_corr_by_type.py`, `gdna_none_guard.py`):**

| metric | baseline | directional fix |
|---|---|---|
| exon→bnd σ²_T (median, cap_on) | 3.47 (gagged) | **1.35** (trusted) |
| intron→bnd σ²_T (median, cap_on) | 0.39 | 4.82 (correctly gagged — sparse) |
| boundary solve corr (mass-wtd) | −0.05 | **+0.12** |
| single-strand region corr | 0.434 | 0.511 |
| AMBIG region corr | 0.503 | 0.609 |
| **gdna_none phantom gDNA (total mass)** | 5.05 M | **3.51 M (−30%)** |

Three-front win: the exon edge is now trusted, the boundary solve ~3× better (still below the 0.35 "resolves"
bar), regions +20% via spillover (better boundary messages → better flank regions), and the zero-gDNA phantom
DROPS 30% (the hard guard). This is the **precision** piece. The **mode** piece — the boundary gDNA message still
saturating to f_g≈1 (the nascent=residual reframe, §8/§11.1) — is still open and is the next step.

### 13.1 The precision fix UNMASKS the mode bug (owner's diagnosis, 2026-07-22) — coupled bugs

Shipping Fix #1 crossed one scenario tripwire: `test_antisense_intronic` ss_0.65 (a NEG single-exon `g2` inside
`g1`'s intron, abundance 0) — t2 mRNA leak **0 → 93** (limit 20). Repro (`anti_intron_repro.py`, gdna=0):

| | baseline | Fix #1 |
|---|---|---|
| t2 leak | 0 | 93 |
| phantom gDNA | 22 | **61 ↑** |
| nRNA | 1424 | 1321 |

The phantom gDNA going **UP** (not down) is the diagnostic. **Fix #1 (precision) is correct** — excluding the
spliced density from σ²_transfer was the bug; including it (directionally) is the clean fix. But it gave the
exon→boundary message the precision it deserves while its **MODE still saturates to f_g≈1** (§8/§11 — the
all-gDNA over-call). Before Fix #1 the wrong mode was gagged by the buggy low precision, so it didn't propagate;
now it's trusted, so the **saturated all-gDNA composition bleeds from the exon boundary into the flanking intron**
(phantom gDNA ↑), and the displaced RNA rebalances onto the antisense single-exon (t2 leak). **The precision was
masking the mode bug.** ⟹ the antisense regression is not caused by Fix #1 — it is *exposed* by it. The real fix
is the **transfer MODE** (the nascent-as-residual / correct composition transfer, §8/§11.1), which must resolve
BOTH the boundary saturation AND the antisense leak. Fix #1 is held (uncommitted, in place) until the mode fix
lands, so the coherent whole commits with a passing suite.
