# Calibration initialization & the prior-free first pass

**Status:** design + audit (branch `calib-ambig-init-wip`). Motivated by the unstranded + hybrid-capture
gDNA→RNA collapse (`ambig_dense_10mb`, `gdna_*_ss_0.50_*_capture_on`), where enriched exons that are 70–99 %
gDNA in truth are calibrated to f_g ≈ 0 (13 M-fragment leak; see `scripts/debug/calib_pool_benchmark.py`).
This document specifies **how calibration should initialize every node — value *and* precision — before any
message propagation**, derives the honest density precision, states the prior-free first-pass goal, and
audits the current code against the spec. A standalone math/stats prototype (`scripts/debug/bp_sandbox.py`)
develops and validates the theory with no tool/biology baggage.

---

## 0. The shift in one paragraph

Calibration deconvolves each node's fragment mass into the composition `(f_RNA+, f_RNA−, f_gDNA)` on the
2-simplex. The **value** of a node's composition is worth nothing without its **precision** — a node at
`{0,0,1}` with zero precision is a blank that *will move*; a node at `{0,0,1}` with infinite precision is a
locked structural fact. Everything in calibration is a precision-weighted combination, so the whole method
stands or falls on assigning **honest precisions at initialization**. The failure mode we are fixing is that
the enriched exons currently get pinned to f_g ≈ 0 by an overconfident depleted-background prior *before* the
messages that carry the real gDNA signal can act — and those messages are themselves impotent or silent. The
cure starts at initialization: give every node an honest self-solved precision, let boundaries self-solve
their composition from the spliced-vs-unspliced structure (a strand-free gDNA/RNA readout), and run the first
message pass **prior-free**, propagating **composition** (scale-free) rather than absolute density.

---

## 1. Precision is the currency

Work per component in **log-density** (`log ρ_c`, `ρ_c = f_c · D`, `D` = observed total density) and per node
in **composition log-odds** (`λ = log(ρ_g / ρ_c)`; for a single-strand `+` node `f_g = σ(λ)`, `f_+ = 1−σ(λ)`).
A belief is a Gaussian `(mean, precision)` in these coordinates; beliefs combine by precision-weighted sums.

Two poles of precision, both used at init:

| precision | meaning | node examples |
|---|---|---|
| `0` (variance ∞) | **no information** — the value is a placeholder and *will move* | default nodes; unstranded strand tilt |
| `∞` (variance `0`) | **locked structural fact** — unmovable, blocks change | intergenic gDNA; an off / unexpressed RNA strand |

> **Code convention (load-bearing).** The implementation stores the precision state as a **variance**:
> `var = 0` ⇒ **infinite precision** (locked), `var = ∞` ⇒ **zero precision** (blank). Read every "precision"
> claim here through that inversion (`node_geometry.py` `_type_belief`).

---

## 2. The honest density precision (derivation)

A component density is estimated from a fragment count over an effective length:

```
ρ̂ = n / L_eff
```

- **`L_eff` is FL-based.** For a **boundary crossing**, a fragment of length ℓ crosses a genomic point iff its
  start lies in a window of width ℓ, so the expected crossings at density ρ are `ρ·μ_FL` (μ = mean FL): the
  crossing effective length is `L_eff = μ_FL`. For a **contained region** of length L, a fragment fits iff it
  is short enough, so `L_eff = E[max(0, L−ℓ)]` (the full FL distribution). Spliced (mature RNA) uses the RNA
  FL; unspliced (gDNA) uses the gDNA FL — **different distributions, different L_eff**.

- **Leading order (Poisson).** `n | ρ ~ Poisson(ρ·L_eff)`, so
  ```
  Var(log ρ̂) = Var(n)/E[n]² = 1/n     ⇒     precision(log ρ) = n     (the COUNT)
  ```

- **Overdispersion (Poisson is too generous).** Real counts over-disperse relative to Poisson (GC/mappability,
  local heterogeneity, FL-model error). Add a relative-variance floor φ (a negative-binomial-style dispersion):
  ```
  Var(log ρ̂) = 1/n + φ     ⇒     precision = n / (1 + n·φ)   →   saturates at 1/φ
  ```
  A confident (large-n) measurement **cannot exceed the overdispersion ceiling** `1/φ`. This is the honest cap
  that a bare Poisson `precision = n` would violate.

- **Where fragment length enters precision — the right reading of "length matters."** The FL **mean** sets
  `L_eff` and therefore the density *value*. Because a boundary crossing has a **short** `L_eff` (≈ μ_FL, a few
  hundred bp), it collects **few counts**, so its **count precision is inherently modest** — this is the real
  sense in which "boundary length is short → the estimate is variable." The FL **spread** enters through
  (a) the region `L_eff = E[max(0,L−ℓ)]` shape (it matters most for short regions, where long fragments can't
  fit) and (b) the overdispersion φ (FL-model error). It does **not** independently inflate the count variance:
  given a global FL model, the per-fragment lengths are *ancillary* to n, so the naive "1 count from a 100 bp
  vs a 300 bp fragment ⇒ 3× density" is an artifact of a per-fragment `1/ℓ` estimator; the correct `n/μ`
  estimator has variance `1/n` regardless. Net: **precision is count-driven, the (short, FL-scale) crossing
  length limits the count, and φ caps it.**

- **Composition precision (the key quantity).** The composition log-odds `λ = log(ρ_g/ρ_c)` combines two
  independent log-densities, so its variance is the **sum**:
  ```
  Var(λ) = Var(log ρ_g) + Var(log ρ_c) = (1/n_g + φ) + (1/n_c + φ)
  ```
  **This precision is nonzero even for unstranded data** whenever both counts exist — it comes from the *counts
  of the two components*, not from strand. That is the whole reason a boundary can self-solve its gDNA/RNA
  split with no strand information (§3.4).

---

## 3. Initialization

Initialization runs in three layers over the two poles of precision, before any message.

### 3.1 Default — the blank node
Every node starts at **`{f_+, f_−, f_g} = {0, 0, 1}` (100 % gDNA) at zero precision.** The density value is
irrelevant: precision 0 means the node carries no information and will be moved entirely by its self-solve and
messages. (Matches current code.)

### 3.2 Structure init — the locks (from the region signature, never the counts)
- **Intergenic regions & intergenic↔exon boundaries (TSS / TES):** pure gDNA `{0,0,1}` at **infinite
  precision** — unmovable. They are **RNA sinks** (an RNA message dies at them — RNA cannot exist intergenic)
  and **gDNA sources** for the *population background* (their density anchors the depleted mode of the gDNA
  prior, §6). In message propagation they are **barriers**: they neither move nor relay composition across
  themselves, so one gene's composition cannot leak into the next and the background cannot be injected into an
  expressed exon. (The sandbox intergenic scenario demonstrates the barrier.)
- **Single-strand nodes:** any node where one RNA strand is structurally absent → that strand's density is set
  to `0` at **infinite precision** — a permanent sink for RNA on the off strand. The remaining axis (the ON
  strand vs gDNA) is the 1-D quantity to solve.

### 3.3 Strand-tilt init — the self-solve from per-strand counts
Each node self-solves an initial gDNA/RNA **tilt** from the two-component (gDNA-vs-RNA) strand Beta-Binomial
likelihood of its per-strand counts, with a precision that **scales with strand-specificity**:
- unstranded (κ ≈ ½): the strand likelihood is flat ⇒ **precision → 0** (the tilt carries no information);
- strand-specific (κ ≫ ½): the imbalance is informative ⇒ **high precision** (scaling with the discriminability
  `w(κ) = (2κ−1)²` and the count).

For the flagship unstranded case this contributes *nothing* — which is correct, and is exactly why the boundary
self-solve (§3.4), which is strand-free, is the load-bearing source.

### 3.4 Boundary init — self-solve from spliced vs unspliced (**NEW**)
A boundary is special: it directly observes two physically distinct fragment classes at one genomic position.

- **Spliced fragments = mature RNA, with known strand.** A splice junction carries a single-stranded genomic
  motif (canonical GT/AG), so the spliced fragments' strand is *fixed by the motif*, not inferred. They are
  placed directly into the correct RNA-strand density: `ρ_+ = n_spl+ / μ_RNA`, `ρ_− = n_spl− / μ_RNA`.
- **Unspliced fragments = gDNA (+ sparse nascent).** Nascent RNA crossing an unspliced junction is rare in real
  biology, so at init we assign the unspliced density to gDNA: `ρ_g = n_uns / μ_gDNA`. (With strand-specific
  data the unspliced count can be strand-deconvolved into nascent vs gDNA via the Beta-Binomial, which carries
  its own precision; unstranded, it is all gDNA.)
- **The boundary therefore self-solves a full composition `{ρ_+, ρ_−, ρ_g}` with honest density precision**
  (§2), *without any prior and without strand*. Example (`+` junction, μ_RNA = 200, μ_gDNA = 100, φ = 0.02):

  ```
  100 spliced(+)  → ρ_+ = 100/200 = 0.5
  300 unspliced   → ρ_g = 300/100 = 3.0
  λ = log(ρ_g/ρ_+) = log 6 = 1.79   ⇒   f_g = σ(1.79) = 0.857
  Var(λ) = (1/300 + φ) + (1/100 + φ) = 0.053   ⇒   precision(λ) ≈ 18.8   (NONZERO, no strand used)
  ```

This is the crux: **the spliced/unspliced structure is a strand-free readout of the gDNA/RNA split**, and its
precision comes honestly from the two counts. Boundaries are the only nodes that can localize the gDNA/RNA
composition without a prior in unstranded data.

---

## 4. The mandatory message-free self-solve

Before any message, **every** node must run its local self-solve (§3.3–3.4) to set its internal
`(composition, precision)`. This step establishes the precisions that the first message pass then combines.
It must run on **all** nodes (see the audit §7 item 6 — the current local solve skips G1-locked and empty
nodes, which is a bug for the emitter story).

---

## 5. The first message pass — propagate composition, not density (no prior)

Consider a `+` splice-junction boundary **B** adjacent to its exon region **R**:

```
B: ρ_+ = 0.5, ρ_g = 3.0  → composition f_g = 0.857, precision(λ) ≈ 18.8   (self-solved, §3.4)
R: 10 000 unspliced over L_eff ≈ 500 − μ_FL  → total density D_R ≈ 33, composition UNKNOWN (precision 0)
```

R's total density (33) is a precise direct count; its **composition** is a blank. B carries the composition at
real precision. The message must impute R's **composition**, and R keeps its **own** total density:

- **Composition transfer (correct).** Send B's `λ` (scale-free fraction). R adopts `f_g = 0.857` and computes
  its gDNA density as `0.857 · 33 ≈ 28` — the right answer. Prior-free.
- **Density transfer (the current tool — wrong).** Send B's *absolute* gDNA density `ρ_g = 3.0` and impute it
  as R's gDNA density ⇒ `f_g = 3/33 = 0.09`, gDNA density ≈ 3 — the **undersampling collapse**. B's crossing
  density (3.0) is far below R's interior density (33) because enrichment and crossing geometry differ, so B's
  *absolute density* does not transfer; only its *composition* does.

The sandbox reproduces both exactly (`f_g` → 0.857 vs 0.091). It also shows the two safety properties that make
this work prior-free:

- **Self-gating on gDNA-free libraries.** With `ρ_g = 0` at the boundary (no unspliced gDNA), `f_g → 0`, so R is
  imputed as RNA — no false gDNA, no prior needed.
- **Structural barriers hold.** An intergenic/TSS/TES lock stays f_g = 1 (unmovable) and does not inject its
  gDNA into an adjacent expressed exon (it is a barrier, not an emitter into genes).

> **Caveat — composition transfer is valid between physically related nodes.** B → its own adjacent exon shares
> the same transcript coverage and the same local capture enrichment, so the composition transfers. It does
> **not** transfer across a gene boundary (an intergenic lock and an expressed exon have different compositions
> — hence the barrier), nor is a *single* junction's composition automatically the whole exon's when the exon
> has multiple, disagreeing junctions. The message-precision / edge-gating model that governs *when* a
> composition transfers is the next phase; this document fixes initialization and the composition currency.

---

## 6. Pass-0 KDE from the self-solved boundaries (**NEW**)

Because boundaries self-solve an honest, prior-free composition (§3.4), we can build the gDNA-density prior's
**enriched mode directly from them, before the first message pass** — a "pass-0 KDE". Select the boundaries
that are structurally clean gDNA-enrichment readouts — **intron↔exon, single-strand, enriched (high unspliced
density), RNA-free (low spliced density)** — and use their adjacent exon densities as the enriched training
sample. Validated on the flagship (`boundary_landscape.py`, `stratified_kde_test.py`): 62–94 such boundaries,
whose adjacent exons are 99.7 % gDNA in truth at density ≈ 26 (true plateau ≈ 30); training the KDE on them
moves its enriched mode from 1.58 → 30.1 and recovers the flagship pool gDNA from 10 % → 53 % of truth, while
`gdna_none`/capture-off select **zero** such boundaries (perfectly self-gating). This replaces the current
depleted-only KDE (trained on the already-crushed pass-1 densities), which has no enriched sample and so walls
off the true enriched density by ≈ 354 nats.

The elegant target is to run the **first message pass with no gDNA hyperprior at all** (pure Bayesian:
structure locks + strand tilt + boundary self-solve + composition messages). The depleted background mode
(intergenic + intron) and the enriched mode (boundary self-solves) are both **observed**, not assumed, so a
learned prior — if used — should be a *consequence* of the self-solves, not a hand-set floor imposed before
them.

---

## 7. Audit — current code vs this spec

Audited files under `src/rigel/calibration/` (`node_geometry.py`, `simplex_logodds.py`, `bp_solver.py`,
`strand_likelihood.py`, `gdna_strand.py`, `signature.py`). Verdicts:

| # | spec item | verdict | note |
|---|---|---|---|
| 1 | default `{0,0,1}` at zero precision | **match** | `_type_belief`: `f_g=1`, `var_g=∞` |
| 2 | intergenic/TSS/TES: locked gDNA, unmovable, barrier, background source | **partial** | value+lock+unmovable match; but the *block* is not the init `var=0` — the message-free solve **skips** G1 → `fg_loc=0`, effective precision `1/ε≈1e9`; the node is a **silent wall pinned at f_g≈0**, not the confident gDNA emitter the docstring intends. The emitter behavior exists only behind the default-off `emit_locked` flag. |
| 3 | single-strand off-strand = 0 at infinite precision | **match** | `var_p/var_n = 0` on the off strand |
| 4 | strand-tilt precision vanishes at κ=½ | **partial** | precision rises with κ but does not literally vanish (Jeffreys/Jacobian reference floor gives finite precision); κ=½ flatness holds only when the independently-fit `od_g == od_r` |
| 5 | **boundary self-solve from spliced/unspliced** | **missing** | boundaries init identically to regions off *unspliced* counts; spliced mass never enters ψ or any variance (only the output `rna_mass` and the sweep messages). No density-precision self-solve. **This is the biggest gap.** |
| 6 | mandatory message-free local solve on all nodes | **partial** | exists (`fg_loc/vg_loc`), but **skips G1-locked and empty nodes**, and carries the global prior (not a pure likelihood-only solve) |
| 7 | first pass prior-free (no gDNA hyperprior) | **diverge** | only the Phase-2 **KDE** is pass-2-gated; the global stability floor `ρ_global` + the depleted `ρ_floor` override + the uncapped intron density likelihood are baked into **every** pass. The floor is documented as "the sole bleed-stopper" — removing it naively lets unstranded/AMBIG introns hallucinate nascent RNA. |

**Five biggest divergences to fix (in order):**
1. **No boundary self-solve (item 5).** Implement §3.4: boundaries produce `{ρ_+, ρ_−, ρ_g}` + density precision
   from spliced (motif-stranded RNA) and unspliced (gDNA) counts. This is the load-bearing strand-free signal.
2. **First pass is not prior-free (item 7).** The depleted floor is what kills unstranded enriched exons; the
   goal is to remove the pre-imposed floor and let the depleted+enriched modes emerge from the self-solves
   (§6). The floor's bleed-stopping role must be replaced by the honest composition messages, not a hyperprior.
3. **G1 blocking is a silent wall pinned at f_g≈0 (item 2), not a confident gDNA emitter.** Two coincidental
   code paths (`init var=0` vs the skip-path `1/ε`) both produce "huge precision"; the emitter intent needs the
   `emit_locked` fix promoted from default-off, and the local solve must not skip G1.
4. **The message currency is absolute density, not composition (§5).** The B→R message imputes the source's
   gDNA density onto the destination; it should impute the composition (fraction) scaled to the destination's
   own total density.
5. **Strand-tilt precision floor at κ=½ (item 4)** — small but should be made to vanish cleanly (or documented
   as intended), since the whole unstranded story rests on the strand carrying *zero* information there.

---

## 8. The sandbox — `scripts/debug/bp_sandbox.py`

A self-contained (numpy-only, no biology, no rigel imports) prototype implementing §2–§5: the density
precision, the three-layer init, the boundary self-solve, and a composition-vs-density first message pass on
the toy example, plus the `gdna_none`, stranded, and intergenic-barrier scenarios. Use it to develop and
adversarially test the message theory before touching the tool. Current output confirms: boundary self-solve
precision 18.8 with no strand; composition transfer → R f_g 0.857 (correct) vs density transfer → 0.091
(collapse); `gdna_none` self-gates to 0; strand precision 0 at κ=½ and 47.8 at κ=0.99; intergenic stays locked
and does not corrupt its neighbour.

---

## 9. Open problems / next

- **The message-precision / edge-gating model** (when does a composition transfer between two nodes?) — the
  next phase after this init spec. B→own-exon transfers; across a gene boundary it must not.
- **Nascent in `nrna_present`** breaks "unspliced crossing = gDNA"; the unspliced leg then needs a
  nascent-vs-gDNA split (strand-deconvolved when stranded; a harder problem unstranded).
- **Coverage:** ≈ 18 % of enriched exons have no adjacent clean-gDNA boundary; they rely on the pass-0 KDE mode
  rather than a direct message.
- **Removing the floor without RNA bleed:** item 7 is the highest-risk change; the honest composition messages
  must demonstrably replace the floor's bleed-stopping role across all 24 conditions before it is removed.
