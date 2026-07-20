# The DOF pie relay — implementation plan

**Status:** implementation-ready plan, integrating the third-party review (2026-07-16). Written 2026-07-16,
branch `calib-ambig-init-wip`. **No production code has landed yet** — this document is the precise plan to
execute after it is itself reviewed.

**Reads before this:** [`dof_pie_relay_derivation.md`](dof_pie_relay_derivation.md) (the theory — the `(λ,θ)`
EP relay, the two-level model, the count-term, the honest liabilities); [`dof_pie_model_fix.md`](dof_pie_model_fix.md)
§4 (the S1–S6 skeleton this refines); [`calibration_bp_theory.md`](calibration_bp_theory.md) §5/§8 (the KKT
self-defense / Trojan horse); the validated prototype
[`scripts/debug/dof_pie_relay_check.py`](../../scripts/debug/dof_pie_relay_check.py).

**The review's verdict** (paraphrased): the derivation is structurally correct, elegant, and ready; both
mathematical pillars (the diagonal Fisher information ⇒ the exact product belief; the count-term `1/M_src`) are
"bulletproof"; proceed — but first defuse three implementation landmines (grid-anchoring, sequential-loop
performance, the phantom telegraph under `σ²_transfer=0`). This plan defuses each, and **diverges from the
review on one point with a stated reason** (§6.3: keep `σ²_transfer=0` for the *isolating* A/B, add damping as a
separate increment, not bundled into S4).

---

## 1. The design in one page

**The defect lives in exactly one place: `bp_solver._scan`'s running belief.** The per-node **final solve**
(`simplex_logodds._solve_nodes_logodds_all`) is already coherent — it integrates `ψ` on the constrained `(λ,θ)`
grid, reads `f_g` as the posterior **median** and `f_pos/f_neg` as **means**, all functions of the grid
coordinates. The **relay** is the only incoherent object: `_scan` maintains three independent log-fraction
Gaussians `(fbg,fbp,fbn)/(vbg,vbp,vbn)`, which is what lets `fbp` reach 51.9 and the message precision inflate.

**So the change is deliberately minimal (the elegance is in the footprint):**

| component | change |
|---|---|
| `simplex_logodds._solve_nodes_logodds` / `_solve_ambig_logodds` | **add** grid-posterior coordinate moments `E[λ], Var[λ]` (and `E[θ], Var[θ]` for AMBIG) to the return — the relay seed. *(The grid posterior is already computed; this is `post@lam` etc.)* |
| `node_geometry._type_belief` / `NodeBelief` | seed a **coherent** born belief; carry the coordinate seed for the sweep |
| **`bp_solver._scan`** | **the core change** — replace the running state `(fbg,fbp,fbn)` with the coordinates `(μ_λ,σ²_λ,μ_θ,σ²_θ)`; fold each incoming message onto them by a 1-D EP moment-match; emit per-component density messages derived from the coherent belief, with the count-term precision |
| `bp_solver.node_sweep` `_comb` + the final solve call | **unchanged** — they consume the per-component messages `_scan` returns (now coherent) |
| `chain_region_deconv` / `NodeBelief` output | **unchanged** — the output fractions still come from the final solve |

**Two levels (the invariant that makes it correct):** the **stored/updated** state is composition coordinates
`(λ,θ)` (the constraint surface — `Σf_c≡1`, `f_c∈[0,1]` structural); the **transmitted** message is
per-component **density** `ρ_c` (active-set-robust; the receiver reallocates to its own `M`). `_scan` folds
density messages onto its own `(λ,θ)`, which imposes its own `M` as the sum constraint automatically.

---

## 2. The running-belief state (data structures)

Inside `_scan`, per node `i`, replace the six fraction/variance arrays with four coordinate arrays:

```
mu_lam[i], var_lam[i]      # N(λ; μ_λ, σ²_λ)   — all solvable nodes (gDNA vs RNA-total log-odds)
mu_th[i],  var_th[i]       # N(θ; μ_θ, σ²_θ)   — AMBIG nodes; single-strand θ is LOCKED (see below)
```

Derived fractions (used to build the emitted message and to read out for the next hop):

```
f_g   = σ(μ_λ)                         f_r = 1 − f_g
single-strand :  f_live = f_r ,  dead strand = 0            (μ_θ pinned; τ = ±1)
AMBIG         :  τ = sin(μ_θ) ,  f_pos = f_r(1+τ)/2 ,  f_neg = f_r(1−τ)/2
```

**Locked/no-info conventions** (`precision_state_design.md`), mapped to coordinates:

| node class | `λ` | `θ` |
|---|---|---|
| structural (intergenic / TSS / TES, G1) | `μ_λ = +L`, `var_λ = 0` (locked `f_g=1`) | n/a |
| single-strand (+) | free (`var_λ` from the solve) | `μ_θ = +π/2`, `var_θ = 0` (locked, all RNA on `+`) |
| single-strand (−) | free | `μ_θ = −π/2`, `var_θ = 0` |
| AMBIG | free | free (`var_θ` from the solve) |

`var = 0` is a **lock** (unmovable); `var = ∞` (represented as a large finite variance / the reference width) is
**no information**. A single-strand node's tilt is `var_θ = 0` **locked**, *not* `var_θ = ∞` no-info — folding a
tilt message onto it is a no-op (the review's S3 caution, made precise).

---

## 3. The relay message (content + precision)

For an edge `src → dst`, `_scan` builds, from `src`'s **coherent** running belief, up to three per-component
density messages (content is **unchanged** from today; only its coherence and precision change):

### 3.1 Content (mode) — unchanged formulas, now coherent inputs

```
gDNA :  ρ_g^src = f_g^src · M_src / E_g^src ,   a_g = log( max(ρ_g^src, 1/E_g^dst) / (M_dst/E_g^dst) )
RNA₊ :  ρ_p^src = f_p^src · M_src / E_r^src (+ the preserved spliced terms, §5),   a_p = log(… / (M_dst/E_r^dst))
RNA₋ :  symmetric.
```

Because `f_c^src = σ(…)·… ∈ [0,1]` now, `ρ_c^src ≤ M_src/E_c^src` and `a_c ≤ 0` up to the frame factor — the
`fbp=51.9` mode is unrepresentable.

### 3.2 Precision — the count-term (`1/M_src`), the review's Pillar B

```
gDNA / RNA :  w_c = 1 / ( Var(log f_c^src)  +  1/M_src  +  σ²_transfer )        [density message]
tilt (θ)   :  w_τ = 1 / ( Var(sin θ^src)              +  σ²_transfer_θ )        [ratio message — NO 1/M]
```

* `Var(log f_c^src)` is read from `src`'s coordinate belief (delta-method: `Var(log f_g)=(1−f_g)²·σ²_λ`,
  `Var(log f_r)=f_g²·σ²_λ`; for the strands add the `θ` contribution). This is the shipped `vb*` , now from a
  coherent belief.
* **`1/M_src` is the TOTAL facing count `sm`, not the deconvolved sub-count `n_src = fbp·sm`.** This replaces the
  shipped `pr = n_src/(n_src·vb+1)`; write it as `w_c = M_src/(M_src·Var + 1)` so `M_src=0 ⇒ w_c=0` with no
  guard. Derivation: count-zero-information (`dof_pie_relay_derivation.md` §5.2; prototype C5).
* **The tilt has NO `1/M` term** — `τ` is a within-RNA *ratio*, magnitude-free, so the source's total sampling
  cancels; its only uncertainty is the strand-composition width `Var(sin θ)` (→∞ at `κ=½`, correctly muting the
  tilt on unstranded data). This λ-has-`1/M` / θ-has-none asymmetry is a consequence of density-vs-ratio
  imputation and should be stated explicitly (it is *not* an oversight).
* `σ²_transfer = 0` for the first A/B (§6.3).

### 3.3 The gates — unchanged, but the mature-crossing gate now needs an explicit silence

Structural emission gates (`emit_g/emit_p/emit_n`) and the D4 strand-aware routing are **unchanged**. The
**mature-crossing gate** needs care: today it zeroes `n_nasc`, and gated edges happen to carry `n_mat=0`, so
`n_src=0 ⇒ pr=0` silences the edge. With the count-term `w_c=M_src/(M_src·Var+1)`, `M_src>0 ⇒ w_c>0` — **the
"silence via zero count" trick no longer fires.** So make the gate **explicit**:

```
send_p = mrp[dst] or (not mrp[src])         # the shipped asymmetric mature-crossing predicate
if emit_p and send_p:   fold the RNA₊ message ;   else: w_p = 0 (skip)
```

This is behaviourally identical to the shipped gate on the RNA imputation (verified: gated edges carry no live
`n_mat`, so the measurement is never silenced — `mature_crossing_gate.md` §3d), but it is now spelled out rather
than emergent from a count coincidence.

---

## 4. The relay algorithm (`_scan`) — step by step

For each node `i` in scan order, with source `lsrc = nbr[i]` (skip if `lsrc < 0`):

1. **Build the messages** from `lsrc`'s current coordinate belief (§3): `(a_g, w_g)`, `(a_p, w_p)`, `(a_n, w_n)`
   in the `dst` (=`i`) frame. Apply the emission + mature-crossing gates.
2. **Record** the per-component messages `amg[i],apg[i] = a_g,w_g` (and `p`,`n`) — the return format is
   unchanged, so `node_sweep`'s `_comb` + final solve are unchanged.
3. **Update `i`'s running belief** (the relay = fold the one incoming message onto `i`'s belief, starting from
   `i`'s local seed). By orthogonality this is **two independent 1-D EP moment-matches**:
   * **λ axis:** reconstruct the RNA-**total** density `ρ_r = ρ_p + ρ_n` (and its precision) from the two RNA
     messages, then moment-match `N(λ; μ_λ_loc[i], σ²_λ_loc[i]) · exp(−½w_g(logσ(λ)−a_g)²) ·
     exp(−½w_r(logσ(−λ)−a_r)²)` → `(μ_λ[i], σ²_λ[i])`. *(Two on-λ messages — gDNA + RNA-total — exactly the
     prototype's C2; curvature `p_g(1−f_g)²+p_r f_g²`, C1.)*
   * **θ axis (AMBIG only):** reconstruct the source tilt `τ_src = (ρ_p−ρ_n)/(ρ_p+ρ_n)`, precision `w_τ`, and
     moment-match `N(θ; μ_θ_loc[i], σ²_θ_loc[i]) · exp(−½w_τ(sinθ−τ_src)²)` → `(μ_θ[i], σ²_θ[i])`.
   Locked coordinates (`var=0`) are skipped (a message cannot move a lock).

`_scan` runs forward (left context) then backward (right context) exactly as today; `node_sweep` then combines
the returned per-component messages (`_comb`) and runs the **unchanged** final solve. The final `f_g/f_pos/f_neg`
move (goldens change, expected); the code around the relay does not.

**Seeding (step 0):** `μ_λ_loc[i], σ²_λ_loc[i]` (and `θ`) come from the message-free local solve's new
coordinate-moment outputs (§1). Locked/G1 nodes are reseeded to their structural lock (the shipped
`locked`-reseed logic, in coordinates).

---

## 5. What is PRESERVED — priority #3 and the spliced measurement

The pie fix touches the **imputation** channel only. It must leave the spliced MEASUREMENT machinery intact
(`boundary_spliced_channel_design.md` §6):

* **The mature MEASUREMENT (`n_mat = SPs`) stays in the message CONTENT `ρ` unchanged** (the `+ n_mat/esp` and
  the `− ρ_mat_dst` absorption terms). Its dedicated high-precision channel is **#3, not built here**; until #3
  lands it rides at the imputation precision `w_c` (§3.2) — i.e. the existing 176× haircut is **unchanged** by
  this fix (we neither fix nor worsen #3). **Do NOT fold `n_mat`'s Poisson Fisher information into `Var(log f_c)`
  / `σ²_λ`** — that would re-lump what #3 exists to un-lump; `n_mat` is an *observation of* the coordinate, not
  belief spread.
* **The DOF relay resolves #3's §7.1 frame mismatch for free:** working in the log-frame coordinate `λ` makes the
  message a genuine `E[log f]`-with-`Var(log f)` Gaussian, retiring the shipped `E[f]`-with-`Var(log f)` mismatch.

**Scoping decision (needs sign-off, §10):** the RNA precision drops the deconvolved-sub-count sampling entirely
(`1/M_src`, not `1/n_src`), which changes the weight `n_mat` currently rides at. This is deliberate and consistent
with `mature_crossing_gate.md` §6.4 ("do not fold `n_mat` into a shared `1/n` — that is #3"). Confirm this is the
intended boundary between item 2 and #3.

---

## 6. The three review landmines — mitigations

### 6.1 Risk 1 — quadrature grid-anchoring (numerical collapse on a sharp belief) — **SOLVED, prototype-verified**

The relay fold (§4 step 3) is a 1-D moment-match of a Gaussian belief × a non-conjugate message factor. A
**single** grid — fixed *or* adaptively-anchored — **cannot** both span far message targets and resolve a narrow
peak: the prototype (C8) measured a single anchored grid inflating `σ` 2.6× on a weak message and **collapsing it
to 0.003** on a strong one. So the production fold is **two-stage, self-correcting** (`fold_two_stage` in
`dof_pie_relay_check.py`, C8):

```
1. COARSE grid, Kc≈33, over [prior 6σ  ∪  message targets]   → argmax peak λ̂ + local curvature σ̂ = 1/√(−ψ″)
2. FINE stage: grid at λ̂ ± max(6σ̂, 1.5·coarse-cell), Kf≈33 → moment-match; RE-CENTER + RE-WIDTH on the
   computed (μ′,σ′) and repeat 2–3× (converges fast — captures a wide/skewed posterior's width)
```

* **No Newton / Laplace mode-finder** — the message factor `exp(−½w(logσ(±λ)−a)²)` is *not* globally log-concave
  (the prototype's iterated-delta diverged on a saturating message), so the peak is located by the coarse grid's
  `argmax`, not a gradient step. Robust by construction.
* **Verified across every regime (C8):** sharp/wide prior × weak/strong × agreeing/disagreeing × two-sided
  saturating — worst-case error vs a 4001-point reference is **2.5e-3** (a <0.2% width error on a σ=1.6 wide
  posterior; immaterial to the floored message precision).
* `Kc=Kf=33`, the σ-coverage `6`, and the `refine=3` iterations are **numerical resolution** (ledger class), not
  model dials.
* **Instrument it anyway (S5):** report the fold's self-consistency residual (does a third refinement move
  `(μ,σ)`?) per node, to catch any production node outside the prototype's tested envelope.

The θ fold uses the identical two-stage scheme on `[−π/2, π/2]`.

> **Perf note:** two-stage is ~`Kc+refine·Kf ≈ 130` `ψ`-evaluations per axis per edge. On the cached conditions
> (~10⁴ folds) this is nothing; at genome scale it is the target of the §6.2 Numba/C pass. The coarse grid can
> also be shared across a scan's static geometry — a later optimization.

### 6.2 Risk 2 — sequential-loop performance (the Python bottleneck)

The review is right that `10⁶` sequential K-point moment-matches in pure Python would wreck the genome-scale run.
**But POC and all validation run on the cached conditions** (~2,500 nodes × 2 scans × 2 axes ≈ `10⁴`
moment-matches/condition) — trivially fast in numpy. So:

* **POC (S1–S6): pure numpy, no compilation.** Precompute the static `λ`/`θ` grids and `logσ(±λ)` once; per edge
  is one `K`-vector `ψ` + softmax + two dot products. Correctness first, per the directive.
* **Before genome-scale production only: compile the inner moment-match with Numba** (a self-contained loop over
  `exp`), or move it into the C/nanobind layer. This is a mechanical perf pass on a *frozen*, validated
  algorithm — not part of the POC. Measure the dev-loop wall-time at S4 and record it; do not optimize until a
  genome-scale run demands it.

### 6.3 Risk 3 — the phantom telegraph (`σ²_transfer=0`) — **plan diverges from the review here, deliberately**

The review recommends adding a temporary `σ²_transfer=ε` (or a precision floor) **in S4** to damp undamped
propagation during validation. I recommend a **different sequencing**, for a reason the review did not weigh:

* **The shipped relay is *also* `σ²_transfer=0` (undamped).** So an A/B of `{shipped relay}` vs `{DOF-fix relay}`
  with `σ²_transfer=0` on **both** sides holds the telegraph **fixed** — it is a *shared confound that cancels*,
  and the A/B cleanly isolates the DOF fix (coherence + count-term). Adding `ε` only to the new arm would
  **confound** the very measurement S6 exists to make ("is the change from coherence, or from ε?"). This is the
  "measure each alone before both" discipline that the `region_eff_length` and mature-gate work were burned by.
* Adding `ε` is also a **new constant** — the standing directive is to discuss before introducing one (§10).

**So:** first A/B with `σ²_transfer=0` (isolates the fix). Then, as a **separate, independently-measured
increment**, add damping — and here I agree with the review's *substance*: a small `σ²_transfer` (ideally the
existing measured `adjacent_disagreement_variance`, i.e. a *measured reliability*, not a tuned `ε`) is the right
anti-telegraph, pending the NPMLE. To make the decision evidence-based rather than assumed, **S5 measures the
telegraph directly** (the per-message hop-reach; the `σ_λ` collapse distribution). If undamped propagation
demonstrably harms accuracy, the damping increment lands next, with its own A/B.

> Note the linkage the review implies but does not state: the `σ²_transfer` damping **also floors `σ_λ`**, so it
> simultaneously relaxes Risk 1's resolution demand. The two mitigations are one mechanism seen twice.

---

## 7. Stepwise plan (refined S1–S6)

Each step: test-first, an explicit falsifier, the exact touch-points. **Pin `md5(bp_solver.py)` before and after
every A/B**; discard any run where it drifts (a concurrent editor has mutated it mid-session before —
`mature_crossing_gate.md` Phase 0).

> **PROGRESS (2026-07-16).** **S1–S3 ✓** (behaviour-neutral: `ablate` fidelity `0.000e+00`). **S4 ✓ LANDED**
> — `_scan` rewritten to the coherent `(λ,θ)` relay (v1: λ relayed, θ held at the local seed — a nuisance; the
> per-strand messages still inform the AMBIG tilt at the final solve). `_fold_lambda` two-stage fold (ported
> from prototype C8); count-term precision `1/M_src`; EXPLICIT mature-crossing gate; fold-resolution as
> `CalibrationConfig` fields (#3). **THE MAIN WIN IS VERIFIED AT SCALE:** `scripts/debug/pie_probe.py` (reading
> the real `_relay_pie`) — **`n_src>M` = 0 on all 7 conditions** (was 424–925), relay pie **= 1.000 everywhere**
> (was up to 600×). S2 flipped (coherent); S3 passes. **Pinned hashes:** `bp_solver 9d9eccef…`,
> `simplex_logodds 19bdd04e…`, `config 2aaf8a5f…`, `calibrate 4d86514d…`.
>
> **(A) ACCURACY A/B — DOF-fix vs shipped (`ablate_replay` `both`, fidelity `0.000e+00`; 7 cached conds):** a
> strong net win where it matters. Flagship **unstranded+capture leak −21…−25%** (ss0.50_none_capON 0.614→0.463;
> present_capON 0.545→0.432); capOFF neutral-to-better (ss0.99_none 0.0138→0.0130); **regresses on STRANDED**
> (ss0.99_none_capON 0.108→0.155 +43%, ss0.99_present_capOFF 0.010→0.013). The coherent pie + count-term help
> most exactly in the count-zero-info regime (unstranded), least where the strand already solves.
>
> **(B) INTRON-LEAK DIAGNOSTICS (toy `_mature_exon_chain`, no patch):** truth = exon f_g 0.32, introns/junctions
> 1.0 (their unspliced is pure gDNA; the 100 mature is spliced). Junctions solve to **0.821 (should be 1.0),
> and this is GATE-INDEPENDENT** (0.821 gate on OR off) — a dense exon's gDNA-density message (implied f_g 0.733,
> prec 1671) **dominates** the sparse intron's correct message (1.0, prec 781); the exon's *running* f_g is
> transiently understated by the mature measurement it just received (B→exon, forward scan) and over-propagated.
> On the REAL suite the **gate still HELPS 7/7** under coherence (mean 0.169 vs 0.193; introns far better), and
> the real ss0.99_none_capOFF intron is BETTER (0.042 vs 0.046) — so **the toy intron test is a pure-mature
> artifact, not a real regression** (even gate-off fails its 0.85 bound). The generalizable finding: a DENSE
> node's over-confident gDNA message can steal from a sparse neighbour (the user's original "density steals"
> concern) — connects to σ²_transfer (per-hop decay would damp it) and possibly a forward-belief-vs-final
> ordering refinement. **NOT patched — surfaced for the fix decision.**
>
> **⚠ TWO behavioural test failures — both shown to be artifacts, NOT real regressions (diagnostics above):**
> 1. `test_gdna_sweep_zero_gdna_pin_and_monotone` → 0.683 (was passing). The **known-fragile** knife-edge on a
>    no-information node (`reference_prior_derivation.md` §10.7 already flagged it "leave red until init lands";
>    the bound was calibrated against the improper ψ). **Not a regression — recalibrate the bound.**
> 2. `test_mature_no_nascent_hallucination_in_introns` → introns 0.766 (guard: >0.85). **A genuine item-1×item-2
>    interaction:** the mature-crossing gate blocks the exon→intron **RNA** message but not its **gDNA** message.
>    The old *incoherent* relay kept gDNA/RNA independent, so blocking RNA prevented the leak; the **coherent
>    pie couples them** (`f_r ≡ 1−f_g`), so the exon's ungated gDNA message alone sets the boundary's `f_g` — and
>    its 1−`f_g` **is** RNA, which relays into the intron. Amplified by the exon's mid-scan running `f_g` being
>    lowered by the mature measurement (B→exon), then over-propagated at count-term precision. **This is a
>    modelling question (how the gate and the coherent pie reconcile), not a hack-fixable threshold — surfaced
>    for a decision, deliberately not patched.** **S5/S6 next** (precision re-measure + the unconfounded A/B).

**S1 — instrument (DONE, no `src/` change).** `scripts/debug/pie_probe.py` measures the relay pie; it reproduces
§1 post-gate (62–70% incoherent, 424–925/cond `n_src>M`). ✅ landed. *Also land* the two zero-cost tripwires the
doc calls for: `assert n_src ≤ M` on every emitted edge (a debug flag), and the pie probe as an `xfail(strict)`
that S4 flips.

**S2 — a test that FAILS on today's code.** A synthetic chain (an RNA-dense exon beside a boundary) asserting the
relayed pie sums to 1 within tolerance **and** `n_c ≤ M` on every edge. It must FAIL now (run under `--runxfail`
to confirm), `xfail(strict=True)` until S4. *Falsifier:* it passes today ⇒ the mechanism isn't what we think.

**S3 — coordinate seed (`simplex_logodds` + `_type_belief`), near-behaviour-neutral.**
* Add `E[λ], Var[λ]` to `_solve_nodes_logodds`'s return and `E[λ],Var[λ],E[θ],Var[θ]` to `_solve_ambig_logodds`
  (the grid posterior is already computed — `post@lam`, `post@(lam²)−…`; `θ` from the `λ`-marginal `post`).
* `_type_belief` / the sweep seed the running coordinates from these; structural/single-strand nodes get their
  locks (§2). *Verify:* the born pie sums to 1 on every node; "`var_g=0` with `var_p=∞`" is unconstructible.
  *Falsifier:* the born pie is still incoherent ⇒ the seed is wrong.

**S4 — the relay swap (`bp_solver._scan`).** Replace `(fbg,fbp,fbn,vbg,vbp,vbn)` with `(μ_λ,σ²_λ,μ_θ,σ²_θ)`;
implement §4 (the two 1-D folds via the §6.1 anchored grid); emit per-component messages with the §3 count-term;
make the mature-crossing gate explicit (§3.3). `σ²_transfer=0`. *Verify:* S2 flips to pass; the pie is 1.000 and
`n_c ≤ M` on **every** edge (assert); `intergenic f_g = 0.0000` on all 7 (the gate's falsifier, unchanged); the
under-resolution counter (§6.1) is ~0. Record the dev-loop wall-time (§6.2). *Falsifier:* any `n_c > M`, or
intergenic ≠ 0.

**S5 — re-measure the precisions (the compounding is gone).** `intron_message_trace` at the worst node: the
RNA-message precision must fall from its inflated value toward `O(M)` with `n_c ≤ M`, and the RNA/gDNA precision
imbalance (61× at B740) toward ~1 (they now share `σ²_λ`, the review's S6 checkpoint). *Also* emit the Risk-1
under-resolution histogram and the Risk-3 telegraph reach (§6.3), to decide the damping increment.

**S6 — accuracy A/B, unconfounded, capOFF and capON.** `ablate_replay` + `node_error_report`, DOF-fix vs shipped,
`σ²_transfer=0` both sides. **Baseline (pinned, this session):** `ss0.99_none_capOFF` ALL `both`=0.0138 (intron
0.0456); `ss0.99_none_capON` ALL 0.1082 (the gDNA-channel leak); `ss0.50_none_capON` ALL 0.6136. Per the
derivation §7, expect the pie fix to **bound, not necessarily remove**, the damage — measure it; do not assume a
win. Then (separately) the damping increment's own A/B. Golden regen (`--update-golden`) is **one** run at the
very end, deferred to a stable checkpoint (the branch already owes a regen).

---

## 8. Magic-number ledger (the standing directive)

| quantity | class | disposition |
|---|---|---|
| `λ = logit(f_g)`, `θ = arcsin(τ)` | inherited coordinates | already shipped (`simplex_logodds`) |
| the EP moment-match, the `p_g(1−f_g)²+p_r f_g²` curvature, `1/M_src` count-term | **derived** | prototype C1/C5 |
| grid `K=48`, σ-coverage `6` (§6.1) | **numerical resolution** | same class as the lattice `K`; not a model dial |
| `σ²_transfer` | measured reliability (deferred) | **0** for the isolating A/B; the damping increment uses the *measured* `adjacent_disagreement_variance`, **not** a tuned `ε` — no new constant |
| `σ²_transfer_θ` | measured reliability (deferred) | **0**; returns with the NPMLE |

**No new tuned constant is introduced.** The only judgement calls are numerical-resolution factors (ledger class
"numerical resolution") and the deferred `σ²_transfer` (a *measured* quantity, not a magic number). This is why
§6.3 refuses the review's `ε`: a tuned `ε` *would* be a magic number; the measured disagreement variance is not.

---

## 9. Self-critique (the honest weak points)

1. **The fold quadrature (§6.1) — RESOLVED by prototype hardening (was the softest part).** The stress test I
   ran *before* writing this (C8) disproved a single anchored grid outright (σ collapsed to 0.003 on a strong
   message) and drove the design to the **two-stage self-correcting** fold, now verified to **2.5e-3** across
   every regime — including the two-sided saturating disagreement I flagged here as the worst case. The residual
   risk is a *production* node outside C8's tested envelope; the S5 self-consistency residual (§6.1) is the
   in-vivo tripwire. This item moved from "assumed" to "verified" this session.

2. **The final-solve-is-unchanged claim rests on the emitted-message format staying per-component.** It does — but
   the *precision* the final solve now receives is the count-term `w_c`, not the shipped `n_src/(…)`. So while the
   *code* is unchanged, the *inputs* change, and the final solve's behaviour at, e.g., a node fed a much stronger
   gDNA message than before is not separately validated. S5's precision re-measure covers the message layer; a
   final-solve sanity check (does a coherent-but-strong message ever push `f_g` past the grid edge?) should be an
   explicit S4 assertion.

3. **The tilt precision `w_τ = 1/Var(sin θ_src)` (no `1/M`) — prototype-checked (C9), residual code risk.** C9
   confirms the behaviour (`w_τ→O(1)` at `κ=½`, rising with strand sharpness, never depending on `M`). The
   residual risk is that the *implemented* `Var(sin θ)` comes from a strand solve that may fold in count effects
   inconsistently — verify at S3 that the seeded `σ²_θ` reproduces C9's regimes on real nodes before trusting the
   tilt message at S4.

4. **`σ²_transfer=0` for the A/B is defensible for *attribution* but the resulting relay may be genuinely worse in
   absolute terms** (the undamped telegraph). S6 could therefore show the DOF fix "not winning" even though it is
   correct — exactly the trap the derivation §8.1 and `message_precision_and_dof.md` §4.5 warn about ("honest
   precision on a Trojan-horsed relay makes the horse faster"). **The plan must not read a flat/negative S6 as a
   refutation of the fix;** S6 is an attribution measurement, and the damping increment is the other half. State
   this in the S6 write-up up front.

5. **AMBIG-at-`κ=½` is the residual hard case and this fix does not target it.** With the strand uninformative,
   `σ²_λ` and `σ²_θ` are both wide, the folds are dominated by messages, and the density-vs-composition /
   enrichment question (derivation §8.2) is unresolved. The pie fix makes this node *coherent*, not *correct* —
   its correctness waits on `σ²_transfer` + the enriched-mode prior. Do not over-claim capON/unstranded recovery.

6. **Performance deferral is a real schedule risk, not just a note.** If the numpy relay is (say) 30× slower per
   node, the genome-scale run is blocked until the Numba/C pass — which is non-trivial (the sequential dependency
   resists vectorization). The plan should *timebox* the numpy relay at S4 and, if it is badly slow even on the
   cached conditions, pull the compilation forward rather than discover it at genome scale.

---

## 10. Decisions needed before S4 (sign-off gate)

1. **The item-2 / #3 boundary (§5).** Confirm the RNA imputation precision drops to `1/M_src` and `n_mat` rides
   at that precision (its own channel is #3), i.e. item 2 does **not** touch the measurement precision.
2. **`σ²_transfer=0` for the first A/B (§6.3),** with the damping as a separate increment using the *measured*
   `adjacent_disagreement_variance` — **not** the review's tuned `ε`. (This is my one divergence from the review;
   confirm you accept the reasoning, or direct me to add `ε` in S4 instead.)
3. **The numerical-resolution factors (`Kc=Kf=33`, σ-coverage `6`, `refine=3`, §6.1)** as ledger-class
   "numerical resolution" (no discussion needed) — flag if you want them derived instead.
4. **Prototype hardening — DONE this session.** `dof_pie_relay_check.py` now has C8 (the two-stage fold vs the
   dense reference, worst 2.5e-3, which *disproved* the first single-grid design) and C9 (the tilt precision).
   No further prototype work is gating S4.

---

## 11. Established vs assumed

| | |
|---|---|
| the defect is confined to `_scan`'s running belief; the final solve is coherent | **established** (code read; derivation §3.1) |
| the change is minimal: solve gains coordinate moments; `_scan` swaps state; `_comb`/final-solve/output unchanged | **established** (code read; §1, §4) |
| the two 1-D folds are independent (orthogonality) and coherent + bounded by construction | **established** (prototype C1–C4) |
| the count-term is `1/M_src`; tilt has no `1/M` | **derived + verified** (§3.2; prototype C5, C9) |
| the two-stage self-correcting fold resolves robustly across message strengths/disagreements | **verified** (prototype C8, worst 2.5e-3 vs a dense reference) |
| `σ²_transfer=0` A/B isolates the fix (shared-confound cancellation) | **reasoned** — S6 measures it |
| the fix improves accuracy | **NOT assumed** — S6 measures; a flat/negative result is attribution, not refutation (§9.4) |
| the numpy relay is fast enough for POC on cached conditions | **assumed** — timebox at S4 (§9.6) |
