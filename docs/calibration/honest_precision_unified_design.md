<!-- title: Honest precision — unified design (AUTHORITATIVE synthesis) -->
# Honest precision — the unified design (synthesis)

**Status:** authoritative synthesis, code-audited 2026-06-20 against `main`. Supersedes the *conflicting*
precision claims in `honest_precision_message_design.md`, `global_prior_count_space_derivation.md`,
`count_space_execution_plan.md` §2.2, `count_space_solve_impl_plan.md`, and `calibration_next_steps.md` §2–3 —
this doc reconciles them and states what is true **in the shipped code** plus the one correct way forward.
Every formula below was read verbatim from `src/rigel/calibration/` (file:line cited); the design docs zigzag
and several are stale, so trust this audit over their prose.

---

## 0. What "honest precision" means (the standard)

A node's evidence about its pie `(f₊, f₋, f_g)` must respect three things, and nothing must fake confidence it
does not have:

1. **COUNT** — a component is a discrete count `n_c` of fragments. Its precision is count-space: a density
   `ρ_c = n_c / E_c` has sampling variance `Var(ρ_c) = (n_c + a)/E_c²` (Poisson, `a` = 1 pseudo-obs). A bare
   density erases the count that governs this.
2. **LENGTH** — `E_c` is the bp interval (effective length) the count is measured over. A short `E_c` means a
   density is known imprecisely **even at `n=0`** (the `1/E²` blows up). A bare density erases the length too.
3. **BIOLOGICAL variance** — imputing a density from a **source** node to a **destination** node adds a
   between-node biological spread `σ²_bio(μ)`, learned by the `var~mean` fit. Total imputed-density variance =
   `σ²_bio(μ) + (n_src + a)/E_src²`.

**The one binding rule that ties them together (the "one-currency" principle).** Every term that constrains the
pie must be a **count-space log-likelihood on the simplex, competing one-for-one**. The strand likelihood is `N`
*real* fragments. A message must enter as `N_src` *pseudo*-fragments — **governed by the SOURCE's evidence**,
never amplified by the DESTINATION's mass. The named anti-patterns:
- a delivered precision that scales with **destination mass²** (a `(M_dst/E_dst)²` Jacobian); and
- a variance floor that **vanishes at the simplex walls** (`f → 0/1`) — i.e. lets a "clean 0/1" message pin a
  node with unbounded confidence exactly where pinning is decided.

---

## 1. Current state — audited against code

| Layer | Honest? | Where |
|---|---|---|
| **Strand likelihood** | ✅ **Honest** | `simplex._mixture_strand_loglik` / `strand_likelihood.strand_loglik` |
| **Imputation message** | ❌ **Dishonest delivery** | `bp_solver._message` (575, 581) → `simplex_sweep._local_loglik` (84, 95) |
| **Global gDNA prior** | ❌ **Dishonest delivery** | `bp_solver.node_sweep` (636, 710) → `simplex_sweep._local_loglik` (76) |
| **σ²_bio input** | ⚠️ **Circular** | `bp_solver.fit_gdna_varmean`/`fit_rna_varmean` fed the belief snapshot |
| **Node state** | ⚠️ **Bare fractions** | `bp_solver.NodeBelief` — no `(n_c, τ_c)` |

### 1.1 Strand likelihood — honest ✅ (the exemplar to copy)
`simplex.py:61-79` (3-component) / `strand_likelihood.py:66-76` (2-component collapse):
```
p    = ½·f_g + κ·f₊ + (1−κ)·f₋          # depends ONLY on the tilt t = f₊−f₋ (an identity)
mean = N·p ,   N = u_pos + u_neg
var  = N·p(1−p) + (N·f_g)²·¼·od_g + (N·f₊)²·κ(1−κ)·od_r + (N·f₋)²·κ(1−κ)·od_r
loglik = −½(u_pos − mean)²/var − ½·log(var)
```
Re-derived Fisher information in `f_g`: `I(f_g) = N·(κ−½)²/[p(1−p)] → N·(2κ−1)²` at `p=½`. **The count `N`
enters strictly as a linear multiplier on the strand discriminability `(2κ−1)²`** — the count-zero-information
principle, exactly. Length never enters (the strand split is a count-ratio, not a density). No `(M/E)²` factor;
the variance floor is a constant `1e-9` (not a `μ(1−μ)` floor), and the **overdispersion *inflates* variance at
the walls** (`(N·f_g)²·od_g` grows as `f_g→1`) — correctly *deflating* precision for a confident-looking vertex.
The two overdispersions are fit by pooled method-of-moments (`gdna_strand.py:116-337`) with a **count-independent
gDNA weight** (breaks circularity) and applied **symmetrically** so `κ=½ ⇒ uninformative`. This is the standard.
*(Minor: it is a normal-moment approximation to the exact Beta-Binomial, "small and ~constant in N" per the
docstring — acceptable.)*

### 1.2 Imputation message — honest density-precision, **dishonest delivery** ❌
`bp_solver._message` (`bp_solver.py:559-583`), verbatim:
```
mu_f    = clip((rho_src·E_dst − spliced_dst)/M_dst, 0, 1)          # :566  identity-mean — SOURCE-governed ✅
mu_q    = max(½(rho_src + rho_dst_cur), _EPS)                       # :569  query point for σ²_bio
sigma2  = varmean.predict(mu_q)                                    # :570
pois    = (rho_src + 1/E_src)/E_src                                # :573  SOURCE Poisson sampling noise ✅
tau_rho = 1/(sigma2 + pois)                                        # :574  honest DENSITY precision ✅
tau_f   = tau_rho · (M_dst/E_dst)²                                 # :575  ❌ DESTINATION-MASS² JACOBIAN
var_floor = mu_f·(1−mu_f)·(1/M_src + 1/M_dst)                      # :581  ❌ floor VANISHES at the walls
tau_f   = 1/(1/tau_f + var_floor)                                  # :582
```
applied as `psi −= ½·tau_f·(f − mu_f)²` (`simplex_sweep.py:84` gDNA, `:95` per-strand RNA).

- **The mean and the density precision are honest.** `mu_f` is the source's facing density matched at the dest
  (with the dest's known spliced subtracted); `tau_rho = 1/(σ²_bio + ρ_src/E_src)` is the learned biological
  spread plus the **source's** Poisson noise. This is the `(density, precision)` carrier the docs call shipped —
  and it *is* honest **in density space**.
- **VIOLATION A — the `(M_dst/E_dst)²` Jacobian (`:575`).** Converting a density-precision to a *fraction*-
  precision is algebraically a correct change-of-variables (`f = ρ·E/M ⇒ Var(f) = Var(ρ)·(E/M)²`). But it makes
  the **delivered** precision scale with `M_dst²`, while the destination's own strand Fisher info scales only
  **linearly** in `M_dst` (§1.1). So the message's *relative* weight is `∝ M_dst²/M_dst = M_dst` — **unbounded
  in destination mass**. A few-fragment source pins a fragment-rich destination. This is the exact amplification
  the standard forbids (causally confirmed at node463: a `μ=0, τ=1e7` message overrode `u_pos=16336` of the
  node's own data → phantom `f_g=1`).
- **VIOLATION B — the wall-vanishing floor (`:581`).** `μ_f(1−μ_f) → 0` at `f → 0/1`, so the only count-based
  cap **disappears at the simplex walls** — exactly where pinning is decided. The docstring defends this
  ("a clean 0/1 message still PINS"); it is the named anti-pattern.

### 1.3 Global gDNA prior — honest mean, **same dishonest delivery** ❌ (the user's instinct, confirmed)
`bp_solver.node_sweep` (`:632-710`), recomputed each pass over the `self_solvable` nodes:
```
rho_obs_n = M/E                                                    # :698  raw, f_g-INDEPENDENT axis ✅
rho_mean  = mass-weighted Σ(f_g·M/E) over self_solvable           # :702  exposure-pooled gDNA density ✅ (mean)
gprior_vm = fit_offset(rho_obs_n, (rho_g_n − rho_mean)², …)       # :703  width fit on the raw axis ✅
mu_global  = clip(rho_mean·E/M, 0, 1)                             # :708
geom2_global = (E/M)²                                             # :636  ❌ same Jacobian
tau_global = max(1/(sigma2_n · geom2_global), 1.0)                # :710  = (M/E)²/σ²_n  ❌  + magic 1.0 floor
```
applied as `psi += −½·tau_global·(f_g − mu_global)²` on `~strand_obs` nodes (`simplex_sweep.py:76`).

- The **mean** is an honest count/length object (exposure-pooled rate), and the width's **μ-axis is raw `M/E`**
  (non-circular). So the global got *closer* to honest than the messages — but its **precision is the identical
  `(M/E)²` Jacobian** (`:710` via `:636`), plus a bare `max(…, 1.0)` magic floor.
- **The count-space prototype was reverted.** `grep -rn RIGEL_GLOBAL_ALLNODES src/` → **nothing** (the flag and
  the Poisson-Gamma application live only in `docs/`). The live default is the Gaussian-`(M/E)²` path. So
  the "recent attempt at honest precision" for the global **did not ship** — the user is right.
- **Stale prose:** `CLAUDE.md` calls this a "MAD-spread precision" — false; the live width is a `MonotoneVarMean`
  density-spread × `(E/M)²`.

### 1.4 σ²_bio — fit on the **circular** belief snapshot ⚠️
`bp_solver.py:683-690`: `snap = node_densities(NodeBelief(f_pos,f_neg,f_g,…))` then `fit_gdna_varmean(…, snap)`.
`node_densities` (`:225`) builds `ρ_g = f_g·M/E` — **belief-derived**. The `var~mean` axis *and* response are
both functions of the current `f`. At pass 0 the init is the all-gDNA seed (`f_g=1`, `:266`), and `_edge_varmean`
drops `ρ_g=0` pairs (`:498`), so the surviving gDNA edges are predominantly all-gDNA↔all-gDNA — an **artificially
tight spread** in the early passes that set the basin → `τ_ρ` is over-confident exactly when it most shapes the
solution. The `var~mean` *machine* (`variance_model.MonotoneVarMean.fit_offset`) is correct; only its **input**
is circular. (The global already uses the raw `M/E` axis — the message fits simply don't.)

### 1.5 State + readout (secondary)
- **State** (`NodeBelief`, `bp_solver.py:189-199`) is bare fractions `+ var_*` scalars; the discrete count and
  length are re-derived from geometry each pass, not carried. The "bare density erases count+length" shape,
  applied to the persisted state.
- **Readout** (`simplex_sweep.py:151-154`): `f_g` = posterior **median** (`_fg_median`), `f±` = posterior
  **means** (`_axis_mean`) — they need not sum to 1, and re-emitted messages mix median-`f_g` with mean-`f±`.
  The median is deliberate (avoids a documented +8.7-pt grid-MAP skew) but the projection is an *incidental*
  median/mean mix, not a principled message-vs-readout operator split.
- **EM boundary** (`priors.assemble_priors`): consumes deconvolved **masses only**; `gdna_frac_var` is dropped.
  So all this precision is **internal to the sweep** — it sharpens the calibration *masses*, then the EM gets
  point masses. (Fine for now; flagged so no one expects per-node precision downstream.)

---

## 2. The root cause is shared, and the masking is load-bearing

Both the message (`:575`) and the global (`:636/:710`) carry the **identical `(M/E)²` factor**. At any node the
two pulls are *both* multiplied by `(M_i/E_i)²`, so **their ratio cancels it** — global-vs-message *looks*
balanced only because both are inflated by the same destination-mass² factor. Against the strand likelihood
(which has **no** `(M/E)²`), both over-steamroll the node's own fragments.

**This is why every piecemeal fix regressed** (measured, in the docs):
- count-space **global alone** (messages left at `(M/E)²`): UNSTRAND-0 167K→16K ✅ but STRANDED-0 24→27,933 ❌
  and CAPTURE −4% — the now-honest weak global is steamrolled by the still-amplified messages.
- capping the **wall floor alone** (Jacobian left in): CAPTURE AMBIG 221,742→208,022 ❌ — strips legitimate
  high-`N_src` capture pins while leaving the low-`N_src`/high-`M_dst` phantom.

⇒ **The message delivery and the global delivery must become count-space TOGETHER, in one pass.** Do not trust
the current global's apparent calibration — it is a cancellation of two errors.

---

## 3. The honest form (the target)

One currency: the strand likelihood is `N` real fragments; **every prior is `N_pseudo` pseudo-fragments**,
distributed across the pie, bounded by its own evidence and **independent of `M_dst`**.

### 3.1 Message = `(μ_c, N_src)`
```
μ_c    = clip((ρ_src·E_dst − spliced_dst)/M_dst, 0, 1)         # unchanged — the imputed composition
N_src  = ρ_src² / (σ²_bio(μ_q) + ρ_src/E_src)  = ρ_src²·τ_ρ    # the SOURCE effective count (= 1/CV²_total)
```
`N_src` is the source's count-currency capacity: a sparse or biologically-noisy source carries few pseudo-
fragments; a deep, homogeneous source carries many — but **never more than its own evidence**, and it is
**`M_dst`-independent** (re-derived: `μ_f`- and `M_dst`-invariant). **Delete the `(M_dst/E_dst)²` Jacobian
(`:575`) and the wall-vanishing binomial floor (`:581-582`)** — `N_src` is the honest, non-vanishing
source-evidence bound that replaces both.

### 3.2 Application = the **FULL-simplex** Dirichlet pseudo-count
This is the crux that the prior docs got wrong, **resolved**:
```
log-prior  =  Σ_c α_c · log f_c ,     α_c = N_src · μ_c    for ALL THREE c ∈ {+, −, g}
```
- The **broken** form (rejected in `count_space_solve_impl_plan.md`) was **one-sided**: it added `α_g·log f_g`
  *only* and left the RNA channels at `α=0`. `∂/∂f_g (α_g log f_g) = α_g/f_g > 0` everywhere ⇒ a **monotone
  push to `f_g=1`**, regardless of `μ_g` (measured: a `μ_g=0.2` message → median `f_g≈0.77`).
- The **correct** form spreads `N_src` across **all three** components (`μ₊+μ₋+μ_g = 1`). This is `log Dir(α+1)`
  whose **mode is `f_c = α_c/Σα = μ_c`** — two-sided, pulls *toward* the message content. A zero channel
  (`μ_c=0`, or `free_s` false ⇒ `N_src=0`) drops out harmlessly (`0·log f_c = 0`). At node463: `N_src ≪ M_dst`
  ⇒ the strand likelihood dominates ⇒ no phantom. At a real capture crossing: the source's large genuine
  `N_src` pins — preserving the capture recovery the wall-floor was protecting, **without** a Jacobian and
  **without** a bound that vanishes at the walls.

> Equivalent Gaussian view (if a Gaussian-on-fraction representation is kept instead): `τ_f = N_src/(μ(1−μ))`
> — i.e. the binomial information of `N_src` pseudo-observations, **`M_dst`-independent**. The Dirichlet is
> cleaner (no linearization point, directly in the count currency) and is the recommendation; the *binding
> property either way* is "total pseudo-count = `N_src`, not `M_dst²`."

### 3.3 Global = exposure-pooled rate + count-space confidence, applied as the same Dirichlet
Keep the live exposure-pooled mean (`rho_mean`, `:702` — it already converges to 0 in a zero-gDNA library by
intergenic-length dilution, and to the genome-average under capture). Replace the `(M/E)²`-Jacobian precision
with a **count-space confidence `N_global`** and apply it as the **same** Dirichlet pseudo-count
`α_g += N_global·μ_global` (with the RNA complement `N_global·(1−μ_global)` split across the admissible strands,
so it is two-sided like §3.2). `N_global` scales honestly: large in a confident zero/low-spread library, small
under capture (large between-node spread) ⇒ **emergent deference**, no hand-tuned defer. Drop the magic `1.0`
floor (`:710`). **Land this with §3.1–3.2 in one change** (the masking, §2).

**Refinement (Step-0, the one-currency principle dissolves the `strand_obs` hard fork).** The live code applies
the node-class prior as a **hard XOR** (`simplex_sweep.py:66-77`): `strand_obs ⇒` the U-shaped Beta(½,½)
**Jeffreys**; `~strand_obs ⇒` the **global** Gaussian. That fork is the artifact behind **Bug C** (the 133K
single-strand-exon over-call in unstranded-0, §5): a single-strand exon at κ≈½ has a *flat* likelihood, and the
Jeffreys reference then **vertex-pushes** it to phantom gDNA — and it never sees the global (which would pull it
to 0). In one count-currency there is **no fork**: every node gets the strand likelihood (`N` real fragments)
**+** the global pseudo-count (`N_global·μ_global`) **+** a Jeffreys ½-pseudo-count, all added, all competing.
Deference is then emergent — the strand wins where decisive, the global (→0 in a zero-gDNA library) wins where
the strand is flat — so the κ≈½ single-strand exon settles at `f_g≈0` instead of pinning a vertex. **This is the
mechanism that brings Bug C (the larger half of unstranded-0) into the honest-precision program** rather than
leaving it orphaned. It is a natural part of Step 2 (the global is already recomputed for all nodes; the change
is to *apply* it everywhere as a pseudo-count and demote the Jeffreys from a hard alternative to a ½-count
baseline). Validate it does not regress the *informative* single-strand nodes (κ far from ½), where the strand
likelihood must still dominate the now-also-present global.

### 3.4 σ²_bio on the **raw** observable axis (the prerequisite) — sharpened by Step-0
Step-0 read `_edge_varmean` (`bp_solver.py:466-505`) closely: the message `var~mean` has **three** circularity
vectors, not one, and only two are cleanly fixable. Each fit point is
`(index μ = ½(ρ_g,dst + ρ_g,src), response = (ρ_g,dst − ρ_g,src)², offset = ρ_g/E pair)` over edges kept by
`ok = (ρ_g,dst>0) & (ρ_g,src>0)`, and `_message` then queries the curve at `μ_q = ½(ρ_src + ρ_dst_cur)`.
- **(a) the fit INDEX** `μ` is the belief `ρ_g` pair → **fix:** index by the **raw `M/E`** pair (the gDNA channel's
  `f_g`-independent observable — exactly what the global already uses at `:698`). ✅ clean.
- **(b) the QUERY point** `μ_q` in `_message` is the belief `ρ` → **fix:** query at raw `M/E`. ✅ clean.
- **(c) the fit SET** `ok = (ρ_g>0)` drops nodes the *current belief* has driven to `f_g≈0` (at init nothing
  drops; as RNA nodes converge they fall out) → **fix:** filter on the raw `M>0`, not `ρ_g>0`. ✅ clean.
- **the RESPONSE** `(ρ_g,dst − ρ_g,src)²` is a between-node **gDNA-density** difference — inherently belief-
  derived (there is *no* `f_g`-independent gDNA-density difference). ⚠️ **stays** as the frozen-snapshot value;
  the global's width fit has the identical residual-side circularity (`:705`), so this is consistent, not worse.

So **Step 1 = de-circularize the index, query, and fit-set (a/b/c) to the raw `M/E` axis** (aligning the message
σ²_bio with the global's already-raw axis), accepting the response stays snapshot-derived. The **RNA** channel is
the harder half: its "raw" per-strand observable `(u_s + spliced_s)/E_rna` is contaminated by gDNA on strand `s`
(half of gДНК reads), so there is no clean `f`-independent RNA density — flag the exact RNA construction as an
open item (it may have to keep the snapshot index, or use `u_s` with a gDNA-rate subtraction).

`N_src` is only as honest as `σ²_bio`, so this is a **prerequisite** for §3.1. It is **shippable on its own**
(correctness/robustness; the index/query/set fix is small and principled; **honest expectation = small headline
move** — measure against `honest_precision_baseline.py`, judge with Step 2).

### 3.5 State `(n_c, τ_c)` per component (the door)
Migrate `NodeBelief` to per-component `(n_c, τ_c)` (carry `n₊/n₋/n_g` and their precisions). This (a) stops the
"bare density erases count+length" loss between passes, (b) makes `N_src` a stored quantity rather than re-
derived, and (c) keeps the **forward-backward / relay door open** (a future FB swap forwards `(density,
precision)` = `(n_c/E, τ_c)`; the Dirichlet *mode* `μ_c` is relay-invariant under `N`-discounting, so the GS
Dirichlet and the FB carrier are the same object — no second substrate rewrite). Keep `_fg_median` for the
**readout** (distinct operator from the matched-moment message), and renormalize the readout pie.

---

## 4. The Dirichlet ⇄ (density,precision) zigzag — reconciled

The three docs are not in genuine conflict once two distinctions are made:
1. **one-sided vs full-simplex** — §3.2. The rejection in `count_space_solve_impl_plan.md` is correct *only* for
   the one-sided add-only form; the full-simplex `α_c = N_src·μ_c` is two-sided and is the honest delivery its
   own §1–2 advocate.
2. **GS-now vs FB-later** — `count_space_relay_implementation_plan.md` invariant #3 ("never regress to
   Dirichlet; keep `(density, precision)`") is a **forward-backward relay** requirement, not a GS one. The
   per-hop GS sweep does **not** relay, so the Dirichlet is correct now; and because the readout is the median
   (≈ mode = `μ_c`), the Dirichlet *moves the mean toward uniform under discounting but not the mode* — so it
   does **not** actually break the FB door. Storing `(n_c, τ_c)` (§3.5) satisfies both invariants at once.

`count_space_execution_plan.md` §2.2 ("message precision is already the count-space form, nothing left to
build") is **wrong against code** — the `(M/E)²` Jacobian (`:575`) and the wall-vanishing floor (`:581`) are
live. It is true only *asymptotically* (as `M_dst→∞` the floor's `1/M_src` term dominates), which is not the
operative regime and collapses at the walls anyway.

---

## 5. Scope — what honest precision can and cannot do

- **Reachable:** the capture path (AMBIG-exon `221,421` / total `1,694,088`, oracle `f_g≈0.84`) and the stranded
  message paths — these run through `_message`/global and are governed by the Jacobian + floor. Honest delivery
  is expected to remove the low-`N_src`/high-`M_dst` phantom while preserving the high-`N_src` capture pins.
- **The `κ≈½` identifiability floor — two-faced, and only PARTLY reachable.** At `κ=½`, `p ≡ ½` for every point
  on the simplex (an algebraic identity) — the strand mean carries **zero** information about the gDNA/RNA split;
  only a weak second-order `od`-ratio variance term remains. Step-0 measured both faces:
  - **over-call** in zero-gDNA unstranded (167,791 total), and
  - **under-call** in real-gDNA-under-capture unstranded (`capture-unstranded`: **−367,064 siphon**, 22%).

  The residue that **no prior/precision change can cross** needs *new information* — the **fragment-length
  signal** (gDNA fragments run longer than RNA). BUT a large part of the *over-call* is **not** pure
  identifiability — it is **Bug C**: the **single-strand-exon Jeffreys vertex-push** (133K of the 167K). A
  POS/NEG-class exon at κ≈½ has a flat likelihood (`p = ½f_g + ½f_pos ≡ ½` with `f_neg` locked), so the
  U-shaped Beta(½,½) Jeffreys reference (applied because `strand_obs=True`) pushes it toward a vertex and
  manufactures phantom gDNA. This is reachable: see the §3.3 refinement (dissolve the `strand_obs` hard fork —
  let the count-space global, which → 0 in a zero-gDNA library, govern these flat-likelihood nodes instead of the
  vertex-pushing Jeffreys). So the honest claim: Step 2 + the §3.3 refinement should move **both** the Bug-C
  over-call **and** the capture-unstranded siphon; the *irreducible* κ≈½ residue (the truly balanced AMBIG node
  with no neighbour and no FL) is the FL track's job. **Set expectations per-condition against the recorded
  baseline — do not promise zero unstranded-0.**
- **Out of scope (bigger downstream prize, separate track):** the EM-side eff-len / FL deficit under capture —
  off-capture the EM recovers a 2×-wrong prior to near-perfect, but under capture it goes blind. That, and the
  FL signal, are the higher-ROI work once calibration precision is locked.

---

## 6. Implementation plan (sequenced; each step gated + measured)

**Step 0 — validation harness FIRST. ✅ BUILT (`scripts/debug/honest_precision_baseline.py`).** *(The plan
originally named `gmean_study.py`, which does not exist; the real by-origin oracle harness is
`dissect_regions.py`, and the new harness reuses its cache.)* It runs the production `calibrate()` on cached
payloads (fast, no quant) for the 4 study conditions and reports, **signed**, per condition + per
(strand_class × region_type): `phantom` (Σ max(prod−oracle,0)), `siphon` (Σ max(oracle−prod,0)), `net`,
totals, and the two scalars — vs a **by-origin oracle** (BAM split by the simulator read name). It records a
baseline JSON and, on re-run, **diffs** against it (deterministic: a re-run shows 0 metrics moved) and flags
regressions. This is a **LOCAL dev gate** (the ~500 MB suite lives under `~/Downloads`, not the repo) — *not*
pytest; the portable invariants (factor-1-uniform, the Step-2 node-level unit tests) stay in pytest. Run it with
`OMP_NUM_THREADS=1`. The post-EM `net_gdna_to_rna` gate is the **secondary, slower** layer (regenerate via the
`calibration-benchmark` skill / `evaluate_suite.py`); the harness does not re-quant.

**Step-0 RECORDED BASELINES (current `main`, the CI targets — these REPLACE the stale doc numbers):**

| condition | total prod gDNA | oracle | **phantom** | **siphon** | net | AMBIG-exon |
|---|---|---|---|---|---|---|
| stranded-0 (none, ss.99) | 24 | 0 | 24 | 0 | +24 | **0** |
| unstranded-0 (none, ss.50) | 167,791 | 0 | 167,791 | 0 | +167,791 | 34,455 |
| capture (gdna300, ss.99) | 1,694,088 | 1,662,770 | 37,247 | 5,929 | +31,318 | 221,421 |
| capture-unstranded (gdna300, ss.50) | 1,295,864 | 1,662,928 | 164,573 | 531,637 | **−367,064** | 190,716 |

**Step-0 findings (the stale-number corrections — these refine the plan):**
1. **`capture ≈ 359,766` is NOT reproducible** — it is from an older code state. The real current AMBIG-exon
   capture gDNA is **221,421** (= `count_space_execution_plan` §2.2's 221,742); total contained = 1,694,088,
   net over-call only +31,318 (+1.9%). The CI target is the recorded baseline, not the doc number.
2. **`unstranded-0` is 167,791, not 34,886.** The 34,455 is only the **AMBIG-exon** subset. The dominant
   component is **133,336 on *single-strand* (POS/NEG) exons** — the Jeffreys vertex-push under a flat (κ≈½)
   likelihood ("Bug C"), a *distinct* mechanism (§5).
3. **`stranded-0` AMBIG-exon = 0** ✅ (the intrinsic-solve invariant holds), but total contained = **24** (tiny
   POS/NEG-exon residual) — assert AMBIG-exon = 0 hard, track total = 24 loosely.
4. **capture-unstranded was NOT missing** (it exists in the suite) — and it is a **major, opposite-direction
   error**: a **−367,064 gDNA siphon** (22% *under*-call). The κ≈½ floor cuts both ways (§5).
5. `dissect_regions.py`'s `--soften-global`/`--global-all-regions` flags are **dead** (they monkeypatch
   `bp.global_gdna_prior`, removed) — clean up.

**Step 1 — σ²_bio raw-density fit (§3.4). Shippable alone.** Change only the *input* to
`fit_gdna_varmean`/`fit_rna_varmean` (gDNA: `M/E`; RNA: strand-split `u_s/E_rna`). No new constant; keep
`_MSG_PSEUDOCOUNT=1`. **Gate:** factor-1-uniform holds; `honest_precision_baseline.py` runs and its 4-condition
signed phantom/siphon is recorded (movement noted, **not judged** — Step 1 + Step 2 are judged together per the
user); full pytest suite green. *(De-risks Step 2; honest expectation = small headline move on its own.)*

**Step 2 — count-space delivery, messages + global TOGETHER (§3.1–3.3, §3.5).** Migrate `NodeBelief → (n_c,
τ_c)`; `_message → (μ_c, N_src)` (delete the `(M/E)²` Jacobian + wall-floor); `_local_loglik` applies the gDNA +
per-strand-RNA messages **and** the global as the **full-simplex Dirichlet** `Σ_c (N·μ_c)·log f_c` for **every**
node (dissolving the `strand_obs` Jeffreys/global fork, §3.3); global mean unchanged, global precision →
`N_global`. Keep `_fg_median` readout (renormalize the pie). **Gates:** portable pytest — node463 unit test
(`N_src ≪ M_dst ⇒ f_pos stays ≈1`, no phantom), **two-sided Dirichlet** unit test (`μ_g=0.2 ⇒ f_g pulled toward
0.2`, NOT 0.77), factor-1-uniform, full suite green; local `honest_precision_baseline.py` gate vs the recorded
baseline — stranded-0 AMBIG-exon stays 0; **capture net stays ≈ +31,318** (no new siphon — the masking fix must
hold both directions); **unstranded-0 and capture-unstranded should *improve*** (the Bug-C over-call and the
−367K siphon are the §3.3/§5 targets), but **judged with Step 1, not alone**; goldens regen only with a
reviewed, attributable diff.

**Step 3 (deferred / gated) — forward-backward relay.** Only if a measured non-uniform-gDNA *step* scenario
shows relay beating the global + a windowed-pooled control by a margin that matters in quant error. The `(n_c,
τ_c)` state + the `refit_pass(snapshot)` hook keep this door open; do not build it now.

---

## 7. New-constant ledger (no-magic-numbers review)
- `a = 1` — sampling pseudocount (`_MSG_PSEUDOCOUNT`, already decided).
- **Removed:** the global's bare `1.0` precision floor (`:710`) → replaced by the count-derived `N_global`.
- **No others.** `N_src`, `N_global`, `μ_c` are all derived; `σ²_bio` is the fitted curve.

## 8. Open questions
- **OQ-A (global zero-pin).** In the pure full-simplex Dirichlet, a zero-gDNA library gives `μ_g→0 ⇒ α_g=0`, so
  the global drops out and an *isolated balanced-AMBIG* node has nothing pinning `f_g→0`. Measure whether this
  regresses vs the current Gaussian global; it is likely the **honest `κ≈½` floor** (§5), but if a rate-style
  zero-pin is needed it is a *global-only* principled exception, not full unification.
- **OQ-B (RNA raw axis).** The honest RNA `σ²_bio` axis needs strand-split observed counts, not `M/E` — confirm
  the exact construction (`u_s/E_rna`, spliced handling).
- **OQ-C (readout consistency).** Keep median for `f_g` but renormalize `(f₊,f₋,f_g)` so the pie sums to 1 and
  the re-emitted message currency is consistent (currently median-`f_g` + mean-`f±`).

## 9. Stale-prose reconciliation (fix while here)
- `CLAUDE.md`: global is **not** "MAD-spread precision" (live = `MonotoneVarMean` density-spread × `(E/M)²`);
  `sweep_max_passes` default is **6**, not 4 (`config.py:269`).
- `simplex.py` docstring (`:22-30`): describes a `count` Gaussian-pull term and a `gdna_prior_count` Dirichlet
  that are **not wired** into `_local_loglik` — stale/aspirational; the live terms are strand + sided-spliced +
  Jeffreys/global + imputation only.
- `variance_model.py` docstring: the "power-law fallback" does not exist; the real fallback is the flat
  `_constant_offset`.
