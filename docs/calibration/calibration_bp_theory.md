# Calibration belief propagation — the theory

**Status:** theory reference (branch `calib-ambig-init-wip`). This is the current, consolidated theory for the
calibration deconvolution and its message passing. It **supersedes the message-currency sections of**
`calibration_initialization.md` (kept for its detailed code audit) and reflects everything learned tracing the
unstranded + hybrid-capture gDNA→RNA collapse. **No code should change until this theory is complete and proven
in the sandbox.** Every claim below is validated numerically by three standalone prototypes (numpy only, no
tool/biology):

- `scripts/debug/bp_theory.py` — the consolidated solver (node solve + chain BP + the hard tests)
- `scripts/debug/bp_reconcile.py` — total-density reconciliation + flagship recovery
- `scripts/debug/bp_dependency.py` — the interdependency / Trojan-horse

---

## 1. The problem

A single-pass scan deposits fragment mass into **nodes** — genomic **regions** (exon / intron / intergenic
interiors) and the **boundaries** between them — arranged as a linear **chain** per reference
(…region — boundary — region…). Each node observes a **total fragment density** `D` (mass ÷ an FL-based
effective length). Calibration must **deconvolve** each node's density into up to three **components**:

- **gDNA** (`g`) — genomic DNA contamination; present *everywhere*, its density modulated by capture enrichment.
- **RNA+ / RNA−** (`p` / `n`) — sense / antisense RNA (mature + nascent), present *only where a transcript of
  that strand is active*.

The components **sum to the observed total**: `ρ_g + ρ_p + ρ_n = D`. The output is the split (equivalently the
composition `f_c = ρ_c/D`), which feeds the per-locus EM prior.

**The hard case this theory targets.** Under **unstranded** data (κ≈½, no strand clue) with **hybrid capture**,
gDNA is enriched on-target (exons) exactly like RNA, so an enriched exon is ~70–99 % gDNA in truth yet is
currently calibrated to `f_g ≈ 0` — a 13 M-fragment gДNA→RNA leak (`scripts/debug/calib_pool_benchmark.py`,
`ambig_dense_10mb`). The root causes traced this session: a depleted-background prior that pins enriched exons,
silent locked-seam emission, and — the subject of this doc — a **message model that transmits the wrong thing
at the wrong precision**.

---

## 2. State and the precision currency

Each node carries:
- the **observed total density `D`** (a direct count over an effective length — *high precision, a hard fact*);
- a **belief about the composition**, represented per active component `c` as a Gaussian on **log-density**
  `log ρ_c` with mean `m_c` and **precision `π_c`** (`= 1/Var(log ρ_c)`).

Precision is the currency; the value is worthless without it. Two poles, both used:

| precision | meaning | examples |
|---|---|---|
| `π = 0` (Var = ∞) | **no information** — placeholder, *will move* | a blank exon; the unstranded strand tilt |
| `π = ∞` (Var = 0) | **locked structural fact** — unmovable, absolute defense | intergenic gDNA; an off / inactive RNA strand |

> **Code convention:** the implementation stores the state as a **variance** (`var=0` ⇒ π=∞ ; `var=∞` ⇒ π=0).

The components are **not** free — they live on the scaled simplex `Σρ_c = D` (§5). The number of **degrees of
freedom** is therefore (#active components − 1): intergenic **0** (locked), single-strand **1**, AMBIG **2**.

---

## 3. Honest density precision (derivation)

A component density is `ρ̂ = n / L_eff` (count over an FL-based effective length; boundary crossing
`L_eff = μ_FL`, contained region `L_eff = E[max(0, L−ℓ)]`; gDNA and RNA use their own FL distributions).

- **Poisson leading order.** A Poisson molecular field of rate ρ gives `n ~ Poisson(ρ·L_eff)`, so
  `Var(log ρ̂) = Var(n)/E[n]² = 1/n` ⇒ **precision(log ρ) = n**, the count. (Fragment length is *ancillary*
  to `n` given a global FL model, so the naive per-fragment `1/ℓ` estimate is a red herring; the FL **mean**
  is the normalizer, the FL **spread** shapes `L_eff` and feeds the overdispersion below, and a **short** FL
  gives a boundary its inherently modest count — the true sense in which "length limits precision".)
- **Overdispersion cap.** Real counts over-disperse (GC/mappability, clustering, FL-model error):
  `Var(log ρ̂) = 1/n + φ` ⇒ `precision = n/(1 + nφ)`, saturating at `1/φ`. Poisson alone is too generous.
- **Composition precision.** A ratio of two components has `Var(log(ρ_a/ρ_b)) = Var(log ρ_a) + Var(log ρ_b)`
  — so it is **nonzero for unstranded data** whenever both counts exist. This is why a boundary can self-solve
  its gDNA/RNA split from spliced-vs-unspliced counts with **no strand information** (§4).

---

## 4. Initialization / node self-solve (condensed; full spec + audit in `calibration_initialization.md`)

Before any message, every node self-solves its local belief:

- **Default:** `{g,p,n} = {0,0,1}` at π=0 (blank; the value is irrelevant).
- **Structure locks (π=∞):** intergenic regions and intergenic↔exon (TSS/TES) boundaries → pure gDNA, RNA
  sinks, **barriers** (they neither move nor inject their composition into genic neighbours; their gDNA anchors
  the *population* background, not a per-node message). Single-strand nodes → the off strand locked at 0.
- **Strand tilt:** a gDNA-vs-RNA tilt from the per-strand counts, precision scaling with strand-specificity
  (→0 at κ≈½, high when stranded).
- **Boundary self-solve (the load-bearing prior-free signal):** spliced fragments are **mature RNA with known
  strand** (the splice motif is single-stranded) → placed into the correct RNA strand; unspliced fragments are
  **gDNA** (nascent is sparse) → placed in gDNA; each with its **honest density precision** (§3). A boundary
  thus produces a full `{ρ_p, ρ_n, ρ_g}` with real precision **without strand and without a prior**.

---

## 5. The node solve — a joint constrained MAP (and the self-defense derivation)

A node combines, **per component**, its local belief with all incoming messages (precision-weighted Gaussian
product on `log ρ_c`) to get a per-component target `(m_c, π_c)`. It then solves the composition **jointly under
the sum constraint** — never component-by-component:

```
ρ* = argmin_{ρ_c ≥ 0}  Σ_c π_c (log ρ_c − m_c)²      subject to      Σ_c ρ_c = D.
```

**Derivation (KKT).** With multiplier `μ` for the constraint, `L = Σ_c π_c(log ρ_c − m_c)² + μ(Σ_c ρ_c − D)`:

```
∂L/∂ρ_c = 2π_c (log ρ_c − m_c)/ρ_c + μ = 0
  ⇒   δ_c ≡ (log ρ_c − m_c) = −(μ/2) · ρ_c / π_c .
```

`μ` is shared by all components (it enforces the single sum). So the **deviation of each component from its
target is `δ_c ∝ ρ_c / π_c`** — inversely proportional to its precision. This is the mathematical statement of
**self-defense**:

- a **high-π (confident)** component barely deviates from its target — it defends itself;
- a **low-π (weak)** component deviates freely and **absorbs the constraint** (the excess `D − Σ`targets`).

The `ρ_c` factor is the **log-lever**: a *large* confident component can still supply/absorb some density at a
*small* log-deviation, so finite-precision protection is **proportional, not absolute**. Absolute protection is
exactly what a **π=∞ structural lock** gives (`δ_c → 0`).

**How components "communicate" inside a node.** Only through this constraint. A gDNA message and an RNA message
do not update separate scalars; they jointly move the 1–2 DoF, and the sum reallocates. This is why you **cannot
update one component and residual the rest** — doing so is the *Trojan horse* (§8).

**How differing total densities reconcile.** The messages set per-component *targets*; the receiver's own `D` is
the hard constraint; the KKT redistributes the excess to the lowest-π components. A boundary with total 10 does
not overwrite a region's total 100 — it informs the *split*, and the region allocates its *own* 100. Validated
across three precision regimes in `bp_reconcile.py`.

---

## 6. Message currency — per-component density fields, not composition

**Composition is not transmissible.** The active component set changes along the chain (intergenic 1 →
single-strand 2 → AMBIG 3), so composition jumps discontinuously wherever a strand switches on/off (an intron at
gDNA 100 % abuts an intron→exon boundary at {RNA+ 80 %, gDNA 20 %}). And absolute *totals* are not transmissible
either (the receiver must reallocate to its own `D`). The currency is **per-component DENSITY**, each a field
with its own continuity:

- **gDNA** — active everywhere; a smooth density field in log-space **modulo a multiplicative capture
  enrichment factor**. Its message is `ρ_g`, transferable across any edge, carrying an **enrichment-transfer
  variance** (the density can jump at a probe edge).
- **RNA+/RNA−** — coverage fields active **only along their transcript**, exactly 0 elsewhere. The message is
  `ρ_s`, transferable only between neighbours where strand `s` is active on **both** ends (a deactivation is a
  sink), carrying an **along-transcript continuity variance** (small along a transcript).

**Composition is an emergent OUTPUT** of the node solve (§5), never a message. This is what makes the
active-set change graceful: at an intron→exon seam, gDNA transfers as a *density* from the intron field (not as
the intron's 100 % composition), and RNA+ is a newly-active component informed by the boundary's own spliced
coverage; the composition (≈ {gDNA 45 %, RNA+ 55 %}) falls out of the fields + the sum-to-`D` solve
(`bp_reconcile.py`).

**Message precision is capped and honest:**

```
π_msg(c, i→j) = 1 / ( 1/π_c^belief(i)  +  σ²_transfer(c, i→j) )
```

— the harmonic combination of (a) the source's own belief precision on `c` and (b) the per-edge transfer
variance for `c` (enrichment-transfer for gDNA; along-transcript continuity for RNA). A message can therefore
**never be more precise than the source knows, nor than the edge reliably carries.** This honesty is not
optional: §8 shows self-defense fails against an *overconfident* message.

---

## 7. Message propagation — forward-backward; messages change and decay

The chain is a forest of linear paths, so BP is exact in **one forward + one backward pass**. In the forward
pass a node's forward-belief is `local ⊗ (forward message from the left)`; it emits to the right the message
derived from that forward-belief (so it never feeds a neighbour its own message back — true tree BP). The
backward pass is symmetric. The final per-node belief is `local ⊗ left-message ⊗ right-message`, followed by
the joint constrained solve (§5).

**Messages change as they propagate.** A relayed message accumulates the intermediate nodes' beliefs (its mean
updates) and its **precision decays** — because the per-hop transfer variance `σ²_transfer` **adds** each hop
(`π` after `k` hops `= 1/(1/π_src + k·σ²_hop)`). So a far node is rightly *less* sure than a near one. Validated
in `bp_theory.py` TEST 4: a gDNA source relays rightward while its precision decays `≈ 6.3 → 3.5 → 2.4 → 1.9`,
and against a fixed weak RNA prior the gDNA fraction correctly falls with distance (0.875 → 0.785).

---

## 8. The self-defense principle (the crux)

> **A node's confident components must not be dominated by weaker, less confident messages.** If this holds, a
> node can receive *incorrect* messages at *honest* (therefore low) precision and survive; if it fails,
> incorrect low-precision messages override high-precision belief and the whole chain can be corrupted.

Self-defense is delivered by **three** mechanisms, all necessary:

1. **The joint constrained solve (§5).** By the KKT result `δ_c ∝ ρ_c/π_c`, a confident component barely
   deviates; the weak components absorb the constraint. The **Trojan-horse failure** is *not* doing this —
   updating the addressed (weak) component independently and making a confident component the residual. Demo
   (`bp_dependency.py`): a confident gDNA=30 is dragged to **25.6** by the residual under the naive update, but
   only to **27.3** under the joint solve; in the AMBIG case an RNA+ message is absorbed by the *other* weak
   component (RNA− 4.0→2.5), sparing gDNA.
2. **Honest precision (§3, §6).** Protection is *proportional to the precision ratio*, so a wrong message must
   be *low precision* — which it is, if precision is honest (a wrong gDNA density arises from a low count, an
   unreliable edge, or an uncertain source, all of which cap `π_msg` low). Demo (`bp_theory.py` TEST 1): a
   confident gDNA (π=50) survives an incorrect RNA+ message at honest π=2 (f_g 0.80), but is **dominated** by
   the same wrong message at overconfident π=80 (f_g 0.61). *Overconfidence is the one thing that breaks
   self-defense* — which is why the precision model must never claim more than it earns.
3. **Structural locks (π=∞).** For genuine certainties (intergenic gDNA, an inactive RNA strand), `δ_c → 0`:
   absolute, unconditional defense. Demo: a locked gDNA is unmoved (30.00) even by an overconfident wrong
   message.

**Practical corollary — the current tool's exposure.** The production *per-node* solve is joint
(`_solve_nodes_logodds_all` on the log-odds simplex — good), **but** the sweep's *per-component running belief*
(`fbg/fbp/fbn` in `bp_solver._scan`) is updated one component at a time and **relayed without being
sum-reconciled** — a live Trojan-horse pathway. Fixing that relay is a requirement of the message-passing
rebuild.

---

## 9. What the sandbox proves (the hard problems)

| test (`bp_theory.py`) | question | result |
|---|---|---|
| **Self-defense** | can a confident component survive a wrong weak message? | yes at honest π (f_g 0.80); dominated only if overconfident (0.61); lock = unmoved |
| **Enriched-exon recovery** | does a blank exon recover from gDNA-dominant junctions with **no prior**? | f_g 0.001 → **0.986** |
| **gdna_none self-gating** | do junctions with no gDNA leave the exon as RNA? | f_g **0.001** (no false positive) |
| **Relay + decay** | do messages change and confidence attenuate with distance? | gDNA fraction 0.875 → 0.785, precision decays per hop |
| **Total-density reconciliation** (`bp_reconcile.py`) | how does a boundary (D=10) inform a region (D=100)? | region allocates its own 100; split set by precisions; flagship recovers f_g 0.98, robust to undersampling |
| **Trojan horse** (`bp_dependency.py`) | does a weak message on one component flip a confident one? | naive: yes (30→25.6); joint solve: no (30→27.3) |

**The theory is internally consistent and reproduces both the desired recoveries and the known failure modes.**

---

## 10. Implications for the code (do not implement yet)

When we do touch code, the theory prescribes:

1. **Node solve:** keep the joint constrained solve; ensure it is the *only* place composition is formed.
2. **Message currency:** per-component **density** with the enrichment (gDNA) / continuity (RNA) transforms;
   composition is never transmitted.
3. **Message precision:** `1/(1/π_src + σ²_transfer)` — honest, capped, decaying per hop. Never overconfident.
4. **Fix the relay:** the per-component running belief must be **sum-reconciled** (solved jointly) before it is
   relayed, or the Trojan horse re-enters through the message path.
5. **Prior-free first pass:** the depleted-background floor should not pre-pin enriched exons; the depleted mode
   (intergenic/intron) and the enriched mode (boundary self-solves / pass-0 KDE) are *observed*, and the
   composition messages + self-defense should carry the deconvolution, not a hand-set prior.

## 11. Open theory items (before code)

- **The per-edge transfer variances** `σ²_transfer(g)` (enrichment) and `σ²_transfer(s)` (RNA continuity) — how
  are they estimated from the data (not hand-set)? This is the remaining quantitative piece; it decides, at
  every edge, which component is the confidently-imputed one.
- **AMBIG (2 DoF)** stress cases: opposite-strand overlapping transcripts where both RNA strands *and* gDNA are
  active — confirm the joint solve + self-defense behave under two simultaneous weak/strong messages.
- **Nascent in `nrna_present`** breaks "unspliced = gDNA"; the unspliced leg needs a nascent-vs-gDNA split
  (strand-deconvolved when stranded; open when unstranded).
- **Convergence / order-independence** of the forward-backward on realistic multi-junction exons (the flagship
  503 had disagreeing junctions) — verify the joint solve is stable to message order.
