# The DOF pie relay — an elegant derivation of coherent belief propagation on the composition manifold

**Status:** theory + design derivation. Written 2026-07-16, branch `calib-ambig-init-wip`.
**Audience:** a reviewer who does not know the codebase. Everything needed to check the argument is stated
here; the numbers are reproduced by a standalone script (§12). **No production code changes with this
document** — it establishes the theoretical ground and a validated prototype, per the agreed sequencing
("solid on theoretical grounds and prototypes first").

**Companion documents.**
* [`dof_pie_model_fix.md`](dof_pie_model_fix.md) — the **measurement + prescription** (item 2). Its §1 measures
  the defect; this document derives the fix it prescribes. *(That doc's stated justification — "it bounds item
  1's residual" — is stale; the case is re-derived here on the **post-gate** residual, §3.3.)*
* [`calibration_bp_theory.md`](calibration_bp_theory.md) — the consolidated BP theory. Its §5 (the joint
  constrained solve / KKT self-defense), §6 (per-component density messages), and §8/§10.4/§11 (**"the
  per-component running belief must be sum-reconciled before it is relayed"**) predicted exactly this fix. This
  document is its coordinate-space realization.
* [`reference_prior_derivation.md`](reference_prior_derivation.md) §10 — the resolved `c=½` reference and the
  `(λ,θ)` parameterisation with the verified orthogonality `I_{f_g,τ}=0`. This document **inherits** those
  coordinates; it does not re-open them.
* [`boundary_spliced_channel_design.md`](boundary_spliced_channel_design.md) — priority #3 (the spliced
  MEASUREMENT vs nascent IMPUTATION split). §6 states what this fix must **preserve**.
* **Prototype:** [`scripts/debug/dof_pie_relay_check.py`](../../scripts/debug/dof_pie_relay_check.py) — numpy
  only, no biology, runs standalone. Every quantitative claim below is a labelled check in it.

---

## 1. The problem in one paragraph

At each of ~10⁶ genomic nodes the calibrator infers a **composition** `(f_pos, f_neg, f_g)` — sense-RNA /
antisense-RNA / gDNA — of the node's observed fragment mass, and passes belief between adjacent nodes by a
forward-backward belief-propagation sweep. The **per-node solve is coherent** (it integrates the log-posterior
on a constrained grid where `f_pos+f_neg+f_g≡1`). The **relay is not**: the running belief that a node emits and
carries is maintained as **three independent Gaussians** on `(log f_g, log f_pos, log f_neg)`, updated one
component at a time. Nothing ties them to a composition. The consequence, measured in production (§3.3): on
**62–70 %** of solvable nodes the relayed "pie" does not sum to 1, and on **424–925 nodes per condition** a
single component fraction exceeds 1 — a boundary relaying that it is "5191 % RNA". This document derives the
correct relay: **maintain the running belief in the free coordinates of the composition manifold** (`λ` for a
single-strand node, `(λ,θ)` for an AMBIG node) and **fold the messages onto those coordinates**. Coherence,
boundedness (`f_c ∈ [0,1]`), and `n_c ≤ M` then hold **by construction**.

---

## 2. Setup: the composition manifold and its coordinates

### 2.1 Components, densities, fractions

A node observes total mass `M` and, per component `c`, has an effective length `E_c` (a precomputed geometric
constant). The latent quantity is the per-component **density** `ρ_c` (fragments per effective base). The
**composition** is `f_c = ρ_c E_c / Σ_{c'} ρ_{c'} E_{c'}`, so `f_c > 0`, `Σ_c f_c = 1`. Two node classes:

| class | live components | composition d.o.f. |
|---|---|---|
| **single-strand** | gDNA `g`, RNA on the one live strand `r` | **1** |
| **AMBIG** | gDNA `g`, RNA₊ `p`, RNA₋ `m` | **2** |
| **structural** (intergenic / TSS / TES) | gDNA only | **0** (locked `{0,0,1}`) |

### 2.2 The coordinates — inherited, not proposed

The composition manifold is a simplex; its free coordinates are already chosen, validated, and **in production**
(`simplex_logodds.py`; `reference_prior_derivation.md` §10.7):

```
    λ = logit(f_g) = log(f_g / f_r) ∈ ℝ           f_g = σ(λ),   f_r = 1 − f_g = σ(−λ)
    θ = arcsin(τ)  ∈ [−π/2, π/2]                   τ = sin θ ,   f_p = f_r(1+τ)/2,  f_m = f_r(1−τ)/2
```

Two facts make `(λ,θ)` the **natural** coordinates — not merely a convenient chart — and they are the load
bearing reason the belief factorizes (§5.1):

1. **Information-orthogonality.** For the multinomial observation the Fisher information is **diagonal** in
   these coordinates: `I_{f_g,τ} = 0` exactly (verified to ~1e-11, `reference_prior_bb_check.py` C1;
   `reference_prior_derivation.md` §10.1). `λ` (gDNA-vs-RNA-total) and `θ` (the strand tilt) carry *disjoint*
   information.
2. **The reference factorizes in them.** The Berger–Bernardo reference prior is
   `Beta(½,½)_λ ⊗ Beta(½,½)_θ`, and under `θ = arcsin(τ)` its tilt term **vanishes identically**
   (`reference_prior_derivation.md` §10.3). The relay and the reference live in the same frame.

`λ` is the axis calibration actually delivers (the gDNA-vs-RNA split); `θ` is the strand-tilt **nuisance**.

---

## 3. The defect, precisely

### 3.1 The solve is constrained; the relay is not

The per-node solve (`_solve_nodes_logodds_all`) integrates `ψ` on the `(λ[,θ])` grid, reading `f_g = σ(λ)` and
`f_active = 1−f_g` — one free parameter, both fractions functions of it, `f_g + f_active ≡ 1` by construction.
**The coherent, 1-d.o.f. concept is already implemented — in the solve.**

The relay (`bp_solver._scan`) is three independent updates with nothing linking them (schematically):

```
fbg ← geometric-mean(fbg_local, gDNA message)        # a Gaussian on log f_g
fbp ← geometric-mean(fbp_local, RNA₊ message)        # a Gaussian on log f_p   (INDEPENDENT)
fbn ← geometric-mean(fbn_local, RNA₋ message)        # a Gaussian on log f_m   (INDEPENDENT)
```

`log f_g`, `log f_p`, `log f_m` are treated as **three free variables**. They are not: on a single-strand node
they are two deterministic functions of the one coordinate `λ`; on an AMBIG node, three functions of `(λ,θ)`.
Treating them as free is the whole defect. This is `calibration_bp_theory.md` §8's **Trojan horse** — "the
per-component running belief is updated one component at a time and relayed without being sum-reconciled" — and
its §10.4 requirement: *"the per-component running belief must be **sum-reconciled** (solved jointly) before it
is relayed."*

### 3.2 Two distinct failures, one cause

**(a) Incoherence.** `median(f_g) + mean(f_p) + mean(f_m) ≠ 1` even at birth (median-vs-mean), and the
independent relay compounds it. A relayed triple is not a composition, so `f_p` can exceed 1 and
`n_p = f_p·M` can exceed `M`.

**(b) Curvature in the wrong coordinate.** A gDNA message constrains `log f_g`; an RNA-total message constrains
`log f_r = log(1−f_g)`. Both are functions of the **single** `λ`. Folding them independently ADDS their
precisions, `p_g + p_r`. Folding them correctly onto `λ` uses the chain rule:

```
    d log f_g / dλ = 1 − f_g          d log f_r / dλ = −f_g
    λ-curvature   =  p_g·(1−f_g)²  +  p_r·f_g²     ≤  p_g + p_r ,   =  (p_g+p_r)/2  at  f_g = ½.
```

The independent relay over-states the `λ`-precision **by up to 2×** (prototype **C1**). *This is the one
correct, non-retracted identity salvaged from the superseded `message_precision_and_dof.md`* (its retraction
banner: "two messages contribute `p·[(1−f_g)²+f_g²]` of λ-curvature, which is ≤ p and exactly p/2 at f_g=½ — the
pair **under**-counts, never doubles"). **We frame the fix as "curvature applied in the wrong coordinate," NOT
as the retracted "double-counted measure / extra Jacobian" story** (which `calibration_roadmap.md` §1.2/§4
explicitly reject: "No Jacobian… Do not restore the Jacobian").

### 3.3 The case, re-derived on the POST-GATE residual (task caveat)

`dof_pie_model_fix.md`'s stated justification was conditional on a now-retired predecessor. Re-measuring on the
**shipped** tree (the mature-crossing gate already in), reconstructing the running-belief pie exactly from the
sweep's own combine (`scratchpad/pie_probe.py`, no production patch — the reconstruction is the solver's own
arithmetic, validated by `intron_message_trace.py`):

```
relay pie  f_g+f_p+f_m  over ~2500 solvable nodes, per condition (worst of forward/backward scan):
  ss0.50_none_capOFF : p50 1.42  MAX 38.5   incoherent(>1.01) 76%   components>1 (n_src>M): 717
  ss0.99_none_capOFF : p50 1.05  MAX 39.5   incoherent        70%                            614
  ss0.50_none_capON  : p50 1.31  MAX 6273   incoherent        58%                            765
  ss0.99_none_capON  : p50 1.43  MAX 54588  incoherent        75%                            925   (B616: fbg=603 on 1 fragment)
  ss0.99_present_capOFF: p50 1.05 MAX 18.9  incoherent        63%                            424
```

The **born** (local) pie is nearly coherent (p50 ≈ 1.005, MAX ≈ 1.04–1.17); essentially all the damage is in the
**relay**. The defect is real, large, and **independent of the gate**. Under capture it reaches **600×** (node
B616 relays `fbg = 603` on a single fragment) — the flagship 13 M-fragment leak regime.

---

## 4. Two levels: what is STORED vs what is TRANSMITTED

This is the reconciliation a reviewer must see first, because the two governing theory documents appear to
disagree and do not:

* `calibration_bp_theory.md` §6: **"composition is not transmissible; the currency is per-component DENSITY."**
  (An intron at 100 % gDNA abuts a boundary where RNA turns on — the *active set changes*, so copying a
  composition is wrong; the gDNA **density field** is what is continuous, and the receiver reallocates to its
  **own** observed total.)
* This fix: **maintain the running belief in composition coordinates `(λ,θ)`.**

These operate at **different levels** and are jointly required:

| level | object | why |
|---|---|---|
| **STORED / updated** (per node) | the free coordinates `(λ[,θ])` | the composition manifold's constraint surface; makes `f_c ∈ [0,1]`, `Σf_c=1` structural |
| **TRANSMITTED** (on an edge) | per-component **density** `ρ_c` | active-set-robust; the receiver reallocates to its own `M` |

The receiver **folds the incoming density message onto its own `(λ,θ)`** — and because `f_c = σ(λ)…` are
fractions of the receiver's **own** total, the fold **automatically imposes the receiver's `M` as the sum
constraint** (`Σf_c = 1`). The `(λ,θ)` fold *is* `calibration_bp_theory.md` §5's joint constrained solve,
expressed in the free coordinates. Storing coordinates ≠ transmitting composition; `calibration_roadmap.md`
Phase 5 endorses storing the free-axis number and does **not** license transmitting composition. (Prototype
**C7** shows the density message + the receiver's own `M` reconciles a badly-undersampling source correctly.)

---

## 5. The framework

### 5.1 The belief family — a product of independent Gaussians, legitimately

The running belief is a Gaussian per free coordinate:

```
    single-strand :  N(λ; μ_λ, σ²_λ)
    AMBIG         :  N(λ; μ_λ, σ²_λ) · N(θ; μ_θ, σ²_θ)
```

The **product form is exact, not an approximation**: `λ ⊥ θ` (`I_{f_g,τ}=0`), so the second-order belief on the
manifold is a product of independent marginals — the Fisher metric is diagonal, so the natural
exponential-family belief factorizes. "The RNA precision" and "the gDNA precision" cease to be two facts: they
are two views of the single `σ²_λ` (dof §2.1). Asking "does the solver respect the gDNA precision when it steals
RNA?" dissolves — there is one number `σ²_λ`, and moving it moves both fractions. That **is** the pie.

**Locked/no-info conventions** (`precision_state_design.md`): `var = 0 ⇒ locked`, `var = ∞ ⇒ no information`.
Mapped onto coordinates: a structural node locks `λ` (`f_g=1`); a **single-strand** node has `θ` **structurally
locked** (`τ = ±1`, all RNA on the live strand — a `var_θ = 0` lock, not a `var_θ = ∞` no-info); an AMBIG node
has both `λ` and `θ` free.

### 5.2 The message and its precision — the count-term, derived

The imputation model is unchanged (`CALIBRATION_ARCHITECTURE.md` §1.2): impute the **density**,
`ρ_c^{dst} ≈ ρ_c^{src}` (identity mean; no learned slope — a slope ≠ 1 would launder an eff-length bug,
`imputation_variance_model.md` §5). The source emits `ρ_c^{src} = f_c^{src}·M_{src}/E_c^{src}`, which the
receiver reads as a Gaussian observation of its own `log f_c^{dst}`:

```
    mode  a_c = log ρ_c^{src} + log(E_c^{dst}/M_{dst})          (the receiver-frame implied log f_c)
    prec  w_c = 1 / ( Var(log f_c^{src})  +  1/M_{src}  +  σ²_transfer )
                       ╰── composition ──╯   ╰ sampling ╯   ╰ adjacent-pair overdispersion (a PRIOR; =0 now) ╯
```

**The count-term is `1/M_{src}` — the TOTAL count — not `1/n_c` (the sub-count).** This is forced by
count-zero-information and is the point most likely to be challenged (the pre-`simplex_logodds` archive used
`1/n_c`). The derivation: `log ρ_c = log f_c + log M − log E_c`, and composition ⊥ total (count-zero-info), so

```
    Var(log ρ_c^{src}) = Var(log f_c^{src})  +  Var(log M_{src})  =  Var(log f_c^{src}) + 1/M_{src}.
```

The count `M` enters only the **total** sampling. It carries **zero** composition information: the split
`f_c` is *deconvolved*, and its uncertainty `Var(log f_c)` is a separate term set by the strand likelihood and
prior. Using `1/n_c` would model the sub-count as an **independent Poisson observation** — i.e. treat the count
as directly voting on the composition, the exact violation `CALIBRATION_ARCHITECTURE.md` §0 forbids. Prototype
**C5** demonstrates both: when `n_c = f_c·Poisson(M)` (deconvolved) the sampling variance is `1/M = 0.001`,
**50× tighter** than `1/n_c = 0.05`; and at `κ = ½` the composition variance `Var(log f_c) ≈ 2.8` dominates so
`M` buys nothing (`pr ≈ 0.357` for `M=100` and `M=10⁴` alike) — count-zero-information, exactly.

This precision **lowers** with `+1/M` (the correct direction; it avoids the retracted "×0.5 doubles precision"
sign error). Written `pr = M/(M·Var + 1)` it gives `pr = 0` at `M = 0` with no guard and saturates at
`1/Var(log f_c)` — a source can never send more confidence than its own composition belief supports. **This is
algebraically the shipped `bp_solver._scan` message precision with two corrections:** (i) the sampling count is
the **total** `M_{src}`, not the deconvolved `n_src = fbp·sm` (dof §6.4's count-term); and (ii) `Var(log f_c)` is
now read from a **coherent** belief where `f_c ≤ 1`.

> **Provisional, per roadmap ordering.** `calibration_roadmap.md` Phase 6 sequences message precision **last**
> ("every prior σ²_imp verdict is void on the old substrate"). The formula above is derived on the **current**
> `c=½`-Jeffreys substrate (Phase 1, landed) and stands or falls with it. `σ²_transfer` and any overdispersion
> `φ` remain **0** pending the NPMLE fit — see the honest liability in §8.1.

### 5.3 The fold — Expectation Propagation = the coordinate chain rule

The message is a Gaussian in `log f_c(λ)`, a nonlinear function of `λ`; the fold onto the Gaussian belief is
therefore **non-conjugate**. The principled projection is **Expectation Propagation**: multiply the belief by
the message factor and moment-match the result back to a Gaussian.

```
    p̃(λ)  ∝  N(λ; μ_λ, σ²_λ) · ∏_c exp( −½ w_c (log f_c(λ) − a_c)² )
    (μ_λ', σ²_λ')  =  ( E_{p̃}[λ],  Var_{p̃}[λ] )        ← the exact moment-match
```

Because each message informs exactly one coordinate (§5.4), this is a **1-D** integral per axis — `O(K)` with
`K ≈ 15–30` points, cheap at genome scale, and the two axes fold independently by orthogonality.

**The delta-method is the first-order EP** (linearize `log f_c` at `μ_λ`): its curvature contribution is exactly
`w_c·(d log f_c/dλ)² = p_g(1−f_g)² + p_r f_g²` — the §3.2 identity — and its mean update is the message's implied
`λ`, precision-weighted. **But the prototype shows a single delta-step is inadequate**: it under-shoots a strong
message (C2: `f_g` update 0.79 vs exact 2.47), and even iterating the linearization mishandles a **saturating**
message (`log f_r = +2`, implied `f_r = 7.4` — the `fbp=51.9` pathology — has no `λ` root; iterated delta gives
`−0.25` vs exact `−3.37`). **Design verdict: the production fold is the EXACT 1-D quadrature moment-match**; the
delta-method is a small-message diagnostic and the derivation of the curvature identity, not the shipping code.

**Boundedness is structural.** For any `(μ_λ, θ)`, `f_g = σ(μ_λ) ∈ (0,1)`, `Σf_c = 1`, `n_c = f_c·M ≤ M`. A
message that implies `f_c > 1` becomes a bounded pull toward `λ = ±∞` — never a fraction `> 1` (prototype **C2**:
every case, including the saturating one, keeps the pie at exactly 1.000).

### 5.4 The `(RNA-total, tilt)` reframe — one message per coordinate

The two per-strand RNA messages are transformed **at the source** into the two orthogonal coordinates before
sending:

```
    RNA-total density  ρ_r^{src} = ρ_p^{src} + ρ_m^{src}     →  a message on  log f_r  →  folds onto λ
    tilt               τ^{src}   = (ρ_p − ρ_m)/(ρ_p + ρ_m)   →  a message on  sin θ     →  folds onto θ
```

`τ` is a scale-free ratio: it transports with **no eff-length offset and no magnitude dependence**
(`τ^{dst} = τ^{src}` under identity-density imputation, both strands sharing `E_r`). So on the `λ` axis a node
folds two density messages (gDNA + RNA-total), and on the `θ` axis one (the tilt) — each a 1-D EP update,
independent by `λ ⊥ θ`. This aligns the messages with `reference_prior_derivation.md` §10.4's **two-group
(gDNA-vs-RNA-total) axis + nuisance tilt** decomposition exactly.

### 5.5 Self-defense emerges (the mynotes §2 / bp_theory §8 principle)

A weak RNA message (small `w_r`) contributes `w_r·f_g²` of `λ`-curvature — small — so it barely moves a confident
`λ` (large existing precision). Deviation `∝ w_msg/π_total`: the KKT self-defense law `δ_c ∝ ρ_c/π_c`
(`calibration_bp_theory.md` §5). Prototype **C3**: a confident `f_g = 0.900` moves `+0.008` under an honest weak
wrong message (`w=2`) but `+0.231` under an overconfident one (`w=90`) — *"overconfidence is the one thing that
breaks self-defense"*, which is why the honest precision (§5.2) is not optional. On the tilt axis, a strong tilt
message leaves `λ` **untouched** (orthogonality; C3). This is `mynotes.md` §2 realized: *"weak RNA precision must
not move a high DNA precision"* — the 1-DOF and 2-DOF self-defense both fall out of the coordinate fold.

### 5.6 1-DOF reduces to AMBIG-restricted-to-`λ`

A single-strand node is an AMBIG node with one strand structurally dead: `τ` pinned to `±1`, no live tilt d.o.f.,
`RNA-total = the one live strand`. The `λ` fold is then **bit-identical** to the single-strand path (prototype
**C4**, `max|Δ| = 0`). One codepath; the class is a `θ`-lock.

---

## 6. What this PRESERVES — the priority-#3 boundary

The pie fix touches only the **imputation** channel's coherence and precision. It must **not** disturb the
spliced MEASUREMENT machinery (`boundary_spliced_channel_design.md` §6):

* **The mature MEASUREMENT stays an OBSERVATION of the coordinate, with its own count precision.** The
  motif-stranded spliced count `SPs` is a direct, strand-free observation of the destination exon's mature RNA;
  its precision is a **count** (`1/n_eff_mat`), *not* the running belief's `Var(log f_c)`. The fix must **not**
  fold the mature Fisher information into `σ²_λ`, or it re-lumps exactly what #3 un-lumps. In coordinate terms:
  the mature message is one more Gaussian factor on `log f_p` (or `log f_m` by motif strand) folded onto `(λ,θ)`
  — carrying tilt information — but its `w` is a count, contributed by the measurement, not the belief spread.
* **Absorption** (`−ρ_mat_dst` at exon→boundary) and the **forward-backward two-channel combine** are untouched.
* **The DOF relay RESOLVES #3's §7.1 frame mismatch as a by-product.** Today `fbp` is a mean `E[f]` paired with
  `Var(log f)` — two frames. Working in the log-frame coordinate `λ = logit(f_g)` makes the message a genuine
  log-frame Gaussian (`E[log f]` with `Var(log f)`), which §7.1 recommends as the cure.

Priority #3 ships **default-OFF** until the NPMLE `σ²_transfer` lands, so it is not live during this work; the
requirement is only that the pie fix leaves its channels intact.

---

## 7. Properties and falsification tests

Binary, checkable claims (the prototype checks the theory; production tests check the implementation later):

| claim | status |
|---|---|
| the relay pie sums to 1 (`|Σf_c − 1| < 1e-9`) on every solvable node | **by construction** (C2) |
| `f_g, f_p, f_m ∈ [0,1]` for every `(λ,θ) ∈ ℝ × [−π/2, π/2]` | **by construction** (C2) |
| `n_c = f_c·M ≤ M` on every edge (today violated on 424–925 nodes/cond) | **by construction** (C1/C2) |
| two on-`λ` messages contribute `p_g(1−f_g)²+p_r f_g² ≤ p_g+p_r`, `=½` at `f_g=½` | **verified** (C1) |
| the fold's exact quadrature moment-match; delta under-shoots strong/saturating msgs | **verified** (C2) |
| self-defense `δ ∝ w/π`; `λ` untouched by a tilt message | **verified** (C3) |
| 1-DOF ≡ AMBIG-restricted-to-`λ`, bit-identical | **verified** (C4) |
| count-term `1/M_src` (count-zero-info), not `1/n_c`; vanishes at `κ=½` | **verified** (C5) |
| a Gaussian message decays precision without moving its mean | **verified** (C6) |
| density message + receiver's own `M` reconciles an undersampling source | **verified** (C7) |
| orthogonality `I_{f_g,τ}=0` on the implemented grids | **verified** (`reference_prior_bb_check.py` C1) |

---

## 8. Honest open liabilities (what this does NOT settle)

A reviewer should weigh these; none is introduced by the fix, but the fix interacts with them.

### 8.1 `σ²_transfer = 0` removes per-hop decay — the largest liability

Every prior precision document calls the adjacent-pair overdispersion `σ²_transfer` (a.k.a. `σ²_bio`, the
"communication variance") a **prerequisite**, not optional (`honest_precision_unified_design.md` §3.4;
`imputation_variance_model.md` §2; `precision_state_design.md` §4). It is currently **0** (deferred to the NPMLE
projection). Consequences (prototype **C6**): with `σ²_transfer = 0` the relay precision **never decays per hop**
— a message propagates **undamped** across arbitrarily many hops, which `count_space_relay_implementation_plan.md`
§9 flags as the phantom-gDNA-laundering risk, and which over-precises RNA relay (RNA `σ²_bio` is genuinely
`> 0`, "spiky"). **The DOF fix neither creates nor closes this** — it is inherited. What it *does* do is kill the
**amplitude** inflation (C1/C2: `n_c ≤ M`, curvature `≤ p`) that made an undamped message *loud*. The honest
sequencing: the coherent relay and the honest precision must land **together** —
`message_precision_and_dof.md` §4.5 measured that "honest precision on a Trojan-horsed relay makes the horse
faster" (0.118 → 0.154). The pie fix is the reconciliation that makes the precision safe; `σ²_transfer` is the
separate NPMLE work that adds the decay.

### 8.2 Density vs composition under capture ENRICHMENT — scoped out

`bp_sandbox.py` shows density transfer can *undersample* (a sparse boundary imputing a low `ρ_g` into a dense
enriched exon → `f_g` collapse), while composition transfer over-copies at active-set changes. The `(λ,θ)` fold
+ the receiver's own `M` reconciles the **magnitude** (C7), but **what direction to send under capture
enrichment** — where `ρ_g` genuinely jumps at a probe edge — is the `σ²_transfer` / enrichment-stratification
problem (`calibration_bp_theory.md` §6.1, §11), **not** the pie fix. The pie fix makes the relay coherent; it
does not resolve enrichment.

### 8.3 Forward-backward is not exact here

`count_space_relay_implementation_plan.md` §9 demoted "FB is exact on a path": the strand Beta-Binomial is
non-Gaussian and `f_g`-degenerate at `κ=½`, and the per-strand `free_s` gate splits the chain. The FB sweep with
Gaussian `(λ,θ)` beliefs is a principled **EP smoother**, not an exact inference — claim it as such.

### 8.4 The tilt-axis reference sign — already resolved, but watch the test

`calibration_roadmap.md` Phase 1 flagged a counterintuitive `τ`-term sign that once regressed
`test_gdna_sweep_zero_gdna_pin_and_monotone`. `reference_prior_derivation.md` §10 **resolved** this: under
`θ = arcsin(τ)` the tilt reference term vanishes identically (no term written). The relay's tilt **message** fold
(§5.4) is separate from the reference and introduces no such term. No new sign trap — but the failing test is a
known landmine to re-derive when production init lands, not to trip over.

---

## 9. Relationship to the shipped code (formalization, not novelty)

A reviewer who knows the codebase must see this stated plainly: **the coordinate and the precision formula are
essentially the shipped ones**; the novelty is **coherence** and the **count currency**, not a new latent.

| shipped today (`bp_solver._scan`, `simplex_logodds`, `_type_belief`) | this document |
|---|---|
| latent `λ = logit(f_g)`, `f_g = σ(λ)` (solve); `θ = arcsin(τ)` grid | **same** — inherited (§2.2) |
| running belief = 3 independent `(fb_c, vb_c)` on `log f_c` | **replaced** by `(μ_λ,σ²_λ[,μ_θ,σ²_θ])` — the defect fix (§5.1) |
| `pr = n_src/(n_src·vb+1)`, `n_src = fbp·sm` (deconvolved sub-count) | count-term → **`1/M_src`** total (§5.2, dof §6.4) |
| combine per component in `log f_c` space (independent) | **fold onto `λ[,θ]`** by EP moment-match (§5.3) |
| `_type_belief` seeds `var_g/var_p/var_n` independently | seed one coordinate `(μ_λ,σ²_λ[,μ_θ,σ²_θ])`; derive the three (dof S3) |
| final grid solve already coherent (median readout) | **unchanged** — the message operator (moment-match) ≠ the readout (median) |

The **readout stays the posterior median** from the (unchanged) final grid solve; only the **message/relay**
operator becomes the coordinate moment-match (`honest_precision_unified_design.md` §1.5: matched moments for the
message, median for the readout — two different operators). This closes the standing `honest_precision` OQ-C
("renormalize the pie") **by construction**.

---

## 10. Implementation sketch (deferred to post-review — no code lands now)

For orientation only; the stepwise plan is `dof_pie_model_fix.md` §4 (S1–S6), to be executed **after** this
derivation is reviewed. The two touch-points:

* **`node_geometry._type_belief` (dof S3):** produce the seed `(μ_λ, σ²_λ)` (+ `(μ_θ, σ²_θ)` for AMBIG) from the
  same local grid solve it already runs (the grid posterior's `λ`/`θ` moments), and derive `var_g/var_p/var_n`
  from them. The born pie sums to 1; "`var_g=0` with `var_p=∞`" becomes unconstructible.
* **`bp_solver._scan` (dof S4):** carry `(μ_λ, σ²_λ[, μ_θ, σ²_θ])` as the running belief; on each edge, form the
  per-component density message (unchanged content, incl. the preserved spliced terms §6), transform the two RNA
  messages to `(RNA-total, tilt)`, and fold via the **exact 1-D quadrature** onto `λ` (gDNA + RNA-total) and `θ`
  (tilt). Recover `fb_g/fb_p/fb_n = σ(μ_λ)…` wherever the loop needs a fraction. `n_src = f_c·M ≤ M` — assert it.

Magic-number ledger: **none added.** `λ`, `θ`, the EP moment-match, the chain-rule curvature, and the count-term
`1/M` are all derived; `σ²_transfer` stays the existing `0`; the quadrature `K` is a numerical resolution.

---

## 11. Established vs assumed

| | |
|---|---|
| the relay pie is incoherent on 62–70 % of solvable nodes; `n_src>M` on 424–925/cond; MAX 600× (capture) | **measured** (§3.3, post-gate) |
| the solve is constrained; the relay is not (the Trojan horse) | **read from code + measured** (`calibration_bp_theory.md` §8) |
| `(λ,θ)` are information-orthogonal (`I_{f_g,τ}=0`) ⇒ the product belief is exact | **verified** (`reference_prior_bb_check.py`) |
| two on-`λ` messages contribute `p_g(1−f_g)²+p_r f_g²` (≤ p, ½ at ½) | **verified** (C1); the kept identity from the retracted DOF doc |
| the fold = EP moment-match; exact 1-D quadrature required (delta under-shoots) | **verified** (C2) |
| coherence + boundedness + `n_c ≤ M` hold by construction | **structural** (C2) |
| self-defense + orthogonality + 1-DOF reduction | **verified** (C3, C4) |
| count-term `1/M_src` (count-zero-info) not `1/n_c` | **derived + verified** (C5) |
| `σ²_transfer=0` ⇒ undamped propagation (RNA over-precision, phantom-gDNA risk) | **inherited liability, documented** (§8.1, C6) |
| the fix improves production accuracy | **NOT established** — measure at dof S6, unconfounded, after implementation |

---

## 12. Numerical validation

`scripts/debug/dof_pie_relay_check.py` (numpy/scipy, standalone, no rigel imports) — reproduces every check
above:

```
C1 wrong-coordinate curvature (p_g(1-f_g)^2 + p_r f_g^2, <= p, = p/2 at 1/2) ......... CONFIRMED
C2 EP fold: exact quadrature required; pie coherent + bounded by construction ........ CONFIRMED
C3 self-defense (delta ~ w/pi) + tilt/lambda orthogonality ........................... CONFIRMED
C4 1-DOF reduces to AMBIG-restricted-to-lambda, bit-identical ........................ CONFIRMED
C5 count-term = 1/M_src (count-zero-info), not 1/n_c; vanishes at kappa=1/2 .......... CONFIRMED
C6 precision-decays-without-moving-the-mean; sigma^2_transfer=0 = undamped ........... DOCUMENTED
C7 density message + receiver's own D: robust undersampling reconciliation ........... CONFIRMED
```

The measurement of the production defect (§3.3) is `scratchpad/pie_probe.py` (reconstructs the running-belief
pie from the sweep's `_capture` hook, no production patch).
