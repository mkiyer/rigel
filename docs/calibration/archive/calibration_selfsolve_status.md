# Calibration node self-solve — state, open issues, and next steps

**Status:** living status doc (branch `calib-ambig-init-wip`). Companion to `calibration_design.md` (the design +
equations). This doc captures **where the boundary/node self-solve initiative stands**, the **open issues and
unknowns**, and **how we intend to proceed** — a single place to re-orient from. No production code has changed;
everything below is theory + sandbox-validated (`scripts/debug/*`).

---

## 0. The goal

Fix the **unstranded + hybrid-capture gDNA→RNA collapse** (`ambig_dense_10mb`, `gdna_*_ss_0.50_*_capture_on`):
enriched exons that are ~75–99 % gDNA in truth calibrate to `f_g ≈ 0`, leaking ~13 M gDNA fragments to RNA. The
lever is a correct, honestly-weighted **boundary self-solve** that makes junctions confident gDNA emitters, so the
message pass can carry gDNA into the blank enriched exons.

---

## 1. Settled design (validated)

### 1.1 The model (locked to the biology)
- **Splicing happens at boundaries.** A **region** (exon/intron interior) has **no spliced** fragments.
- **Fragment compatibility:**
  - exon **unspliced-contained** → mature RNA **or** nascent RNA **or** gDNA (3-way);
  - intron **unspliced-contained** → nascent RNA **or** gDNA (2-way);
  - boundary **unspliced-crossing** → nascent RNA **or** gDNA only (**no mature** — mature crossing a junction is
    *spliced*);
  - **spliced** → mature RNA, strand fixed by the genomic motif.
- **`nascent ≈ 0` is the initialization assumption** (a boundary's unspliced-crossing → gDNA), **revised during the
  sweep**. It is *exact* for `nrna_none`; it is the accepted, revisable initial estimate otherwise.
- **Consequence:** only **boundaries** can self-solve gDNA-vs-RNA from spliced/unspliced (their crossing has no
  mature). A **region cannot** — its unspliced carries contained-mature; it relies on strand balance + the gDNA
  prior + boundary messages. *(This corrects an earlier mis-scoped region-projected A/B; see §5.)*

### 1.2 Boundary self-solve — what works
Validated in `scripts/debug/boundary_ab.py` against oracle boundary-crossing truth
(`gdna = gdna_uns / (gdna_uns + nascent_uns)`):

| scenario | single-strand boundary `f_g` MAE (prod → §6) | reading |
|---|---|---|
| flagship (unstranded `nrna_none`) | **0.455 → 0.006** (0.545 → 0.994, true 1.0) | decisive win |
| stranded `nrna_none` | 0.132 → 0.138 (≈ tie) | both correct |
| stranded `nrna_present` | 0.043 → 0.036 | strand peels nascent |

**Confident-emitter mechanism (the point):** production boundaries emit gDNA at precision **0.41** — effectively
silent, *why* the enriched exons never hear "gDNA here." The honest count precision makes the same boundary emit
at **163** (example `N_uns=517`, cap `1/φ=238`). Where production is *over*confident (stranded grid var 2281), the
honest precision caps at `1/φ`. AMBIG boundaries are perfect (`f_g=1`) with the peel gated to single-strand.

### 1.3 Honest precision — derived (§6.2 of `calibration_design.md`)
The message currency is per-component **density**, so each component carries its own **count** precision — not a
composition of two:
```
precision(log ρ_g)  = N_u / (1 + N_u·φ)      (gDNA, from the unspliced-crossing count)  → cap 1/φ
precision(log ρ_mat) = N_s / (1 + N_s·φ)      (mature, from the spliced count; separate channel)
```
Unified across node types: `Var(log ρ_g) = Var(log f_g) + (1/N_u + φ)` — the composition variance plus the
**count floor** the current code omits (it uses the bare fraction variance → would emit at infinite precision when
`Var(log f_g)→0`). Boundaries: `Var(log f_g)=0` under nascent≈0 → the count precision. Regions unstranded:
`Var(log f_g)=∞` → placeholder, deferred to the sweep.

### 1.4 Theory foundations (validated)
- **Count conservation** (`count_conservation.py`): the node solve conserves *counts* (`ΣN_c=N`), not densities
  (density = count/E, E per-component). Self-defense KKT preserved (`δ_c ∝ N_c/π_c`), reduces to the old solve at
  equal FL.
- **Uniqueness under count conservation** (`uniqueness_count.py`): unimodal in all realistic regimes; the one
  multimodal case (all-targets-tiny "push-up") cannot occur at init, is eliminated by an honest prior, and has a
  unique global min a grid resolves deterministically.
- **Self-defense / message currency / recovery** (`node_solver.py`, `bp_reconcile.py`, `bp_dependency.py`): the
  earlier BP theory (S1–S8).

### 1.5 The unified-solver plan (agreed in principle)
There are two redundant solve blocks today — init (`init_beliefs` → grid + `_type_belief`) and the sweep
(`_local_solve` → grid + `solvable`-gate). Collapse to **one solver** used by init (no prior/messages) and the
sweep (prior + messages), passing the **real spliced** (not `_zero_spl`), dissolving `_type_belief`, the
`solvable`-gate, and the mature channel into it.

---

## 2. The value default — RESOLVED (it was tilt-vs-tie-break, not the frame)

Prototyped in `scripts/debug/density_frame.py`. The finding:

"Initial unspliced composition = 100 % gDNA" is a hard default (`f_g=1`) and needs **no tuning constant**. The
magic number I briefly introduced was the wrong *device*, not the wrong frame:

- A **tilt** `+c·log(f_g)` drives `f_g→1` with a strength that scales with `c` (`f_g` slides 0.83 → 0.95 → 1 as
  c = 1 → 5 → 10) — a magic number.
- A **tie-break** `−ε·(1−f_g)` merely breaks the flat-direction degeneracy toward the `f_g=1` vertex, and the
  result is **INSENSITIVE to ε across 7 orders of magnitude** (`1e-2 … 1e-9`): absent strand → `f_g=1`
  (nascent≈0), strand imbalance → peeled, all ε-invariant. A tie-break, not a strength.

**It works in the EXISTING log-odds solver** — the tie-break gives the identical ε-invariant result there as in a
per-component count/density frame. So **"log-odds vs density frame" was a red herring**: the real issue was the
DEFAULT (uniform-`f_g` reference → 0.5 vs. a **nascent≈0 vertex tie-break** → 1) and the *device* (tilt vs
tie-break). **No re-architecture to a density solver is needed.**

**Second finding — a real count-zero-info violation to fix alongside.** The strand **overdispersion** term
(`od_g`/`od_r`) prefers `f_g=0.5` at κ=½ via its `−½·log(var)` normalization: a 50/50 split has lower
Beta-Binomial variance (od ∝ component-count²), so the *count magnitude* tilts the gDNA-vs-RNA VALUE — forbidden
by count-zero-info. It scales with `od` (→ 0.5 at od=0.05; harmless at the flagship's fitted `od_g≈0`), and it is
present in production (`simplex.py`). **`od` must set the precision, not the value.**

**So the value default is settled** — two targeted, magic-free changes in the existing per-node solver:
(1) replace the uniform reference with the **nascent≈0 vertex tie-break**; (2) fix the **od count-zero-info
violation** (od → precision only). No frame re-architecture; §3/#A closes.

---

## 3. Open issues & unknowns

| # | issue | status | approach |
|---|---|---|---|
| **A** | **Value default** (nascent≈0 without a magic number) | **CLOSED (§2)** | It was tilt-vs-tie-break, not the frame. Use a `−ε·(1−f_g)` **tie-break** (ε-invariant) in the existing log-odds solver + fix the **od count-zero-info violation** (od → precision, not value). No re-architecture. |
| **B** | **`φ` (density overdispersion)** — the honest precision cap | proxy only | A *clean* `φ` needs the *deconvolved* gDNA densities (a solve) — the same circularity as #3. The no-solve `estimate_phi` uses **total** density → conflated (0.004 `nrna_none` vs 0.023 `nrna_present`, RNA leaking in). **Ship `φ=0` at init** (honest structure, uncapped); derive the real `φ` post-solve alongside #3. |
| **C** | **`σ²_imp` (message transfer variance)** — over-inflated | deferred (#3) | Total-density, **pooled over regimes** → dominated by enrichment jumps → messages over-weakened. Needs the **enriched-vs-depleted regime-stratified** model (within-regime → the reliable message variance; across-regime → the jump). **This is the end-to-end flagship bottleneck** — the boundary self-solve makes junctions emit correctly, but the message only *lands* once `σ²_imp` is de-conflated. |
| **D** | **Unstranded + nascent (`nrna_present`)** | fundamental limit | gDNA and nascent are unspliced-identical without strand → nascent≈0 overshoots and *nothing* at init (or unstranded sweep) can peel it. Accepted, revisable; the honest edge of the design. |
| **E** | **AMBIG magnitude self-solve** | open sub-problem | `p` depends only on the tilt, so strand can't peel an AMBIG node's `f_g` (holds `f_g=1`, needs the sweep). **With strand-specific data an AMBIG node may still admit a defined precision** — flagged by the user as a distinct, interesting problem to design separately. |
| **F** | **Region contained-mature** | sweep-resolved | An exon's unspliced includes mature; only removable by the **boundary spliced-mature message** (sweep) or strand. Region init `f_g=1` is a placeholder (`var=∞`), not a claim. |
| **G** | **Sweep architecture** — re-solve-from-scratch vs anchored-at-init | resolved by #A | The sweep re-solves each node from the uniform reference (re-hedging to 0.5, discarding the init). With #A's **tie-break default** it re-solves *toward nascent≈0* instead, so boundaries keep emitting their gDNA — no separate anchoring machinery needed. |

---

## 4. How we proceed (next steps)

1. ~~**Settle the frame (#A).**~~ **DONE (§2):** tie-break, not tilt; existing solver; no re-architecture.
2. **Implement the unified solver** (the immediate next step):
   - **value default** = the `−ε·(1−f_g)` nascent≈0 **tie-break** replacing the uniform reference (Jeffreys +
     symmetric Jacobian);
   - **fix the od count-zero-info violation** (od → precision, not the value's `−½·log(var)` preference);
   - **precision** = the count floor `1/(Var(log f_g) + 1/N_u + φ)`, `φ=0`;
   - init = the solve with no prior/messages (dissolve `_type_belief`); sweep = the same solve + prior + messages +
     real spliced + the mature channel (dissolve `_zero_spl` + the `solvable`-gate);
   - **Validate**: boundary A/B (from production), goldens, and no catastrophic end-to-end regression.
3. **De-conflate `σ²_imp` (#C / #3).** The regime-stratified transfer variance — the step that makes the confident
   boundary emissions actually *land* in the enriched exons (the end-to-end flagship fix). Derive the clean `φ`
   here too (same post-solve, regime-stratified machinery).
4. **Then**: the AMBIG self-solve (#E, with the user's design), and the `nrna_present`-unstranded question (#D).

**Sequencing note:** step 2 (boundary self-solve + honest precision) sets up correct, confident *emissions*; step
3 (`σ²_imp`) is what lets them *land*. Both are needed for the end-to-end flagship recovery — step 2 alone
improves the init/boundary state but may not move the end-to-end until step 3.

---

## 5. Evidence base & corrections log

**Sandboxes (`scripts/debug/`):**
- `boundary_ab.py` — boundary self-solve A/B vs oracle boundary-crossing truth (§1.2).
- `node_solver.py` — the general node solver S1–S8 (self-defense, composition, revision).
- `count_conservation.py` — count-conservation self-defense + FL bias.
- `uniqueness_count.py` — unimodality under count conservation.
- `mom_vs_grid.py` — the **mis-scoped region A/B** (kept as the record of the error in §5 below).
- `scratchpad/pin_value.py` — the `+c·log f_g` experiment that exposed the magic-number/frame problem (§2).

**Docs:** `calibration_design.md` (design + equations); `CALIBRATION_ARCHITECTURE.md` (count-zero-info,
authoritative); `calibration_prior_production_reference.md` (what ships on `main`).

**Corrections made this initiative (so they don't recur):**
1. **Region vs boundary conflation.** A "general node solver / no branch" framing led to judging *boundary*
   self-solve by a *region*-projected metric. Regions have no spliced; only boundaries self-solve. (§1.1)
2. **Composition harmonic precision → count precision.** The honest precision is the per-component *count*
   precision, not `1/(1/N_u+1/N_s)`; the latter conflated the gDNA-vs-mature *split* precision with the
   gDNA-density *message* precision. (§1.3)
3. **MoM value solver retired.** A closed-form method-of-moments *value* regressed (κ≈½ ill-conditioning); the grid
   value is robust. The value question is now the frame question (§2).
4. **Magic-number `φ`/`c` → frame artifacts.** `σ²_imp` is a sweep quantity, not an init overdispersion; the
   `nascent≈0` tuning constant is an artifact of the log-odds frame (§2), not a real parameter.
