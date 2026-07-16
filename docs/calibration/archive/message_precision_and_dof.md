# Message precision and the degrees-of-freedom problem

> ## ⚠ PARTIALLY RETRACTED — 2026-07-15
>
> **§5 (the DOF "one bug in three costumes" thesis) is WITHDRAWN. It contains a sign error.**
>
> §5.1 argued that two messages on one axis make messages *~2× too strong*, and that this predicts the
> `σ²_imp` dial should want ×0.5. **That is inverted.** `σ²_imp` is a **variance**: halving it *doubles*
> message precision. A term that is 2× too strong wants ×**2**. §4.4's own table header says so
> (*"message strength ↑ as scale ↓"*). The two dials are not commensurable — ×0.5 *halves* the Jeffreys term
> and *doubles* the message term — so **"both want ½" was numerology**, not evidence.
>
> The premise is also false independently: two messages contribute `p·[(1−f_g)² + f_g²]` of λ-curvature,
> which is **≤ p** and exactly **p/2** at f_g=½ — the pair *under*-counts, never doubles.
>
> **What actually causes the over-call:** the 1-DOF prior is exactly `kde + 0.5·λ` — an improper linear ramp
> whose strength is set by the grid half-width `_DEFAULT_L`. The prior alone returns **f_g = 0.9994**.
> §0.1 and §4's measurements stand and were the clue; §5's explanation of them was wrong.
>
> **Read [`prior_ramp_and_bp_roadmap.md`](prior_ramp_and_bp_roadmap.md) instead** for the corrected
> derivation, the verdicts, and the roadmap.

**Status:** working document, circulated for review and input (branch `calib-ambig-init-wip`, 2026-07-14).
**Purpose:** consolidate what we now know — theoretically and empirically — about **belief-propagation message
precision** in Rigel calibration, and the **degrees-of-freedom (DOF)** problem we believe underlies it.
**Ask of the reader:** §8 (open questions). We think we have found *one* bug wearing three costumes; we want
that challenged.

**Companion documents** (this doc integrates, does not replace them):
`CALIBRATION_ARCHITECTURE.md` (the authoritative information/precision model) ·
`calibration_bp_theory.md` (the BP theory + sandbox proofs — **most relevant**) ·
`calibration_initialization.md` (the node self-solve) · `npmle_struggles.md` + `npmle_roadmap.md` (the pass-0
gDNA-rate prior that replaced the KDE/floor) · `oracle_and_benchmarking.md` (how we measure truth).

---

## 0. Executive summary

1. **Messages are not the villain — the local prior is.** Tracing per-node error against oracle truth: on the
   hardest condition (zero-gDNA, unstranded, nascent, capture) the strand-only error is **0.51**, adding the
   local **prior makes it 0.97** (near-maximally wrong), and the **messages rescue it to 0.28**. Messages have
   been doing the peeling all along.
2. **The messages are gagged.** The single global `σ²_imp` measured on **total** density is 0.96–3.54, so the
   message precision `pr = n_src/(n_src·σ²_imp+1)` **saturates at 1/σ²_imp ≈ 0.28–1.05 pseudo-observations** —
   no message can outweigh ~1 read, however confident its source. `calibration_bp_theory.md` §6.1 predicted
   this and put the inflation at **~6.2×** in the sandbox.
3. **The Trojan horse is live and measurable.** `calibration_bp_theory.md` §8 warned that the sweep's
   per-component running belief is updated **one component at a time and relayed without sum-reconciliation**.
   We measured it: **100 % of corrupted single-strand nodes receive BOTH a gDNA message and an RNA message on
   their ONE free axis, and the two precisions simply add** (1.16 + 1.22 = 2.39).
4. **Two independent dials both want exactly ½** — halving `σ²_imp` and halving the RNA-parsimony Jeffreys each
   improve the benchmark. We read this as the signature of a **double-count on one-DOF nodes**, not as two
   coincidences.
5. **The theory's own message-precision formula fails when adopted piecemeal.** `π_msg = 1/(1/π_src +
   σ²_transfer)` (theory §6) is *worse* than today's count-based rule in production — because it presupposes
   the **joint constrained solve + sum-reconciled relay** that the code does not implement.

**Working thesis:** *everything is ~2× too strong on one-DOF nodes because we treat one degree of freedom as
if it were two independent ones — in the prior (two density priors multiplied) and in the messages (two
precisions added).* Fix the DOF treatment first; only then are honest, stronger messages safe.

---

## 1. Background: the model in one page

From `CALIBRATION_ARCHITECTURE.md`: a fragment **count carries no intrinsic gDNA/RNA information**; a node's
composition is set by exactly three sources — the **strand likelihood** (the only intrinsic signal; the count
enters only as Fisher information), **cross-node imputation** (messages), and the **population gDNA prior**.

Each node carries an observed **total density** `D = M/E` (a hard fact) and a belief about its **composition**
`(f₊, f₋, f_g)` — sense-RNA / antisense-RNA / gDNA. Per `calibration_bp_theory.md` §2, beliefs live on
**log-density** with **precision `π = 1/Var(log ρ)`** as the currency, with two poles: `π=0` (no information,
will move) and `π=∞` (a structural lock, immovable).

### 1.1 Degrees of freedom — the organising fact

The components are **not free**: they live on the scaled simplex `Σ_c ρ_c = D`. So the number of free
parameters is `(#active components − 1)`:

| node class | active components | **DOF** | consequence |
|---|---|---|---|
| **intergenic / TSS-TES seam (G1)** | gDNA only | **0** | locked; `π=∞`; a barrier, not a source of composition |
| **single-strand (G2)** | gDNA + one RNA strand | **1** | `f_active = 1 − f_g`. **One number describes the node.** |
| **AMBIG (G3)** | gDNA + RNA₊ + RNA₋ | **2** | the 2-simplex; strand-degenerate when κ≈½ |

This is the crux. **A single-strand node has ONE free number.** If gDNA moves, RNA *must* move — they are not
two facts, they are one fact viewed twice. Everything below follows from taking that seriously (or failing to).

---

## 2. What the theory already prescribes (`calibration_bp_theory.md`)

The theory is coherent, sandbox-validated, and **largely not implemented**. Its four load-bearing claims:

**(a) The node solve is a joint constrained MAP (§5).** Combine per component, then solve *jointly* under the
sum constraint — never component-by-component:

```
ρ* = argmin_{ρ_c ≥ 0}  Σ_c π_c (log ρ_c − m_c)²      s.t.      Σ_c ρ_c = D
```

KKT gives the **self-defense law**:

```
δ_c ≡ (log ρ_c − m_c) = −(μ/2) · ρ_c / π_c        (μ shared — it enforces the ONE sum)
```

*A component's deviation from its target is inversely proportional to its precision.* A confident component
barely moves; weak components absorb the constraint. Components communicate **only** through this constraint.

**(b) The message currency is per-component DENSITY, not composition (§6).** Composition is not transmissible
(the active set changes along the chain); totals are not transmissible (the receiver must reallocate to its own
`D`). Composition is an **emergent output** of the node solve, never a message.

**(c) Message precision is honest and capped (§6):**

```
π_msg(c, i→j) = 1 / ( 1/π_c^belief(i)  +  σ²_transfer(c, i→j) )
```

— harmonic in (the source's own belief precision) and (the per-edge transfer variance). *A message can never be
more precise than the source knows, nor than the edge reliably carries.*

**(d) Self-defense has three mechanisms, all necessary (§8):** the joint constrained solve; honest precision;
and structural locks. And the explicit warning:

> **The Trojan-horse failure** is updating the addressed (weak) component independently and making a confident
> component the residual.

with the **practical corollary** naming our code:

> The production per-node solve is joint (`_solve_nodes_logodds_all` — good), **but** the sweep's per-component
> running belief (`fbg/fbp/fbn` in `bp_solver._scan`) is updated one component at a time and **relayed without
> being sum-reconciled** — a live Trojan-horse pathway.

**And a prediction we can now confirm (§6.1):** the global `σ²_imp` is estimated on **total** density, which
conflates (1) real gDNA enrichment transfer, (2) *structural* boundary-vs-region composition differences, and
(3) Poisson sampling. Sandbox: `Var(Δlog total)=0.56` vs `Var(Δlog gDNA)=0.09` — **6.2× inflation**, i.e. we
over-weaken gDNA messages ~6×.

---

## 3. What production actually implements (the gap)

| theory prescribes | production does |
|---|---|
| joint constrained MAP, composition emergent | per-node ψ solve **is** joint (log-odds simplex) ✅ … |
| relay a **sum-reconciled** belief | …but `_scan` updates `fbg`/`fbp`/`fbn` **independently** (three separate 1-D Gaussian products) and relays them **unreconciled** ❌ |
| `π_msg = 1/(1/π_src + σ²_transfer)`, per-component | `pr = n_src/(n_src·σ²_imp + 1)` = `1/(σ²_imp + 1/n_src)` — uses the **count** as the source's uncertainty, ignoring `belief.var_gdna` ❌ |
| per-component transfer variance | **ONE global scalar** `σ²_imp` from **total**-density disagreement ❌ |
| prior-free first pass | a pass-0 prior **including an unbounded RNA-parsimony Jeffreys** ❌ |

Concretely, in `bp_solver._scan` each channel does its own 1-D update, e.g. for gDNA:

```python
pr      = n_src / (n_src * sig_imp + 1.0)      # count-based, global scalar
pt      = pg_loc[i] + pr                        # precisions ADD
fbg[i]  = exp((pg_loc[i]*lfg_loc[i] + pr*mo) / pt)
vbg[i]  = 1.0 / pt
```

and identically for `fbp` (RNA₊) and `fbn` (RNA₋). **Nothing enforces `Σ_c ρ_c = D` on the running belief**, and
the relayed message is built from these unreconciled components.

---

## 4. New empirical evidence

All numbers: the 24-condition `ambig_dense_10mb` oracle suite, calibration-basis (pre-EM), metric
`mwae_fg` = mass-weighted |Δf_g| vs oracle (0 = perfect). Tooling in §9; the full loop now runs in **~9 s**.

### 4.1 The impotence, quantified

`adjacent_disagreement_variance` measured per condition (`scripts/debug/msg_precision_diag.py`):

| condition | σ²_imp | precision cap `1/σ²_imp` |
|---|---|---|
| gdna100 ss0.50 nrna_none cap **off** | 1.599 | 0.63 |
| gdna100 ss0.50 nrna_none cap **on** | **3.541** | **0.28** |
| gdna100 ss0.50 nrna_present cap off | 1.245 | 0.80 |
| none ss0.50 nrna_none cap on | 0.956 | 1.05 |

**A message is worth at most ~0.3–1.0 reads** — and *capture makes it worse* (cap_on 3.54 ⇒ 0.28), precisely
where peeling matters most. This is theory §6.1's prediction, measured in production.

### 4.2 The message trace — messages rescue, the prior corrupts

`scripts/debug/msg_trace.py` compares three states per region node against oracle truth: **strand only**
(no prior, no messages) → **local** (strand + prior, message-free) → **FINAL** (local ⊗ messages).

| condition | strand only | **local (strand+prior)** | FINAL (+messages) | verdict |
|---|---|---|---|---|
| none, ss0.50, nrna_present, cap_on | 0.5072 | **0.9749** | **0.2780** | 0 corrupted, **30 rescued** (7.9e5 mass) |
| gdna100, ss0.50, nrna_none, cap_on | 0.2928 | 0.1626 | 0.1197 | 24 corrupted, 48 rescued |
| gdna100, ss0.99, nrna_none, cap_off | 0.0164 | **0.0025** | 0.0066 | 3 corrupted, 1 rescued |

Read this carefully:

* **Unstranded, zero gDNA:** the prior turns a 0.51 into a **0.97** — near-maximally wrong — and the messages
  *rescue* it to 0.28. **The prior manufactures the over-call; the messages fight it.**
* **Stranded:** local is nearly perfect (0.0025) and messages *mildly hurt* (0.0066) — where the strand already
  knows, messages only add noise.

This directly corroborates theory §8: *"The flagship collapse was in fact a blank exon dominated by a
confident-but-wrong local prior … which is why the first pass must be prior-free."*

**The corrupting term is the RNA-parsimony Jeffreys** `−log(1−f_g)` inside the prior projection: it pushes
`f_g→1`, and at κ≈½ **nothing opposes it**.

| Jeffreys strength | all | ss0.50 | ss0.99 | zero-gDNA | worst corner | over-call | under-call |
|---|---|---|---|---|---|---|---|
| **full (today)** | 0.1181 | 0.2178 | 0.0184 | 0.2238 | 0.2780 | 7.8 M | 1.8 M |
| **half** | **0.1020** | **0.1903** | **0.0137** | 0.1459 | 0.0786 | 3.2 M | 6.6 M |
| **none** | 0.1303 | 0.2416 | 0.0190 | **0.0917** | **0.0217** | **1.5 M** | 14.4 M |

Full ⇒ over-call; none ⇒ under-call (the *old* collapse); **half is better on every aggregate metric**
(all −14 %, stranded −26 %, worst corner −72 %).

### 4.3 The Trojan horse is live (the DOF double-count, measured)

At **corrupted** nodes (local was right, final is wrong), `msg_trace.py` reports:

```
corrupted single-strand nodes: 3   (ss0.99)   /   13   (ss0.50 cap_on)
of these, receiving BOTH a gDNA and an RNA message on that ONE axis:  100 %
their precision split:  gDNA 1.16  +  RNA 1.22  =  2.39     vs local 6.92
```

**Every** corrupted one-DOF node was hit by two messages on one axis, and the two precisions **add**. That is
theory §8's Trojan horse, in production, quantified.

Secondary mechanism: corruption is delivered by **extreme message modes**. Node 1449 (true 1.000, local 0.999,
local precision 76.3) was dragged to **0.695** by a message of precision 2.4 — only **3 %** relative weight —
because the message asserted `f_g ≈ 7e-6` and the combination is a precision-weighted average in **log** space
(an unbounded pull). Conversely a node with local precision 4990 barely moved (0.999 → 0.910): **self-defense
does work where precision is genuinely high** (theory §5's `δ_c ∝ ρ_c/π_c`).

### 4.4 Two independent dials, both wanting ½

σ²_imp scaling (message strength ↑ as scale ↓):

| σ²_imp × | all | ss0.50 | ss0.99 | worst | over-call | under-call |
|---|---|---|---|---|---|---|
| 1.0 (today) | 0.1181 | 0.2178 | **0.0184** | 0.2780 | 7.8 M | 1.8 M |
| **0.5** | **0.1075** | **0.1953** | 0.0196 | 0.1545 | 6.0 M | 3.2 M |
| 0.25 | 0.1169 | 0.2097 | 0.0242 | 0.1023 | 4.9 M | 6.4 M |
| 0.1 | 0.1356 | 0.2420 | 0.0292 | 0.0587 | 4.4 M | 10.2 M |
| 0.03 | 0.1529 | 0.2689 | 0.0369 | **0.0510** | 4.1 M | 13.5 M |

Two facts: (i) **×0.5 is a free win** (all −9 %, worst −44 %) — messages *are* throttled; (ii) beyond that a
**trade-off** appears — the worst corner keeps improving (0.278 → 0.051, −82 %: unstranded RNA **can** be
peeled by messages) while the **under-call explodes** (1.8 → 13.5 M) and **stranded degrades monotonically**
(0.018 → 0.043). *Unstranded wants loud messages; stranded wants quiet ones.* No global constant wins.

So **two separate dials — `σ²_imp` and the Jeffreys — each independently want ×0.5.** We do not think that is a
coincidence (§5).

### 4.5 What failed, and why it matters

| experiment | result | reading |
|---|---|---|
| **Refit `P(ρ)` on the peeled belief** (shipped) | 0.151 → **0.118** (−22 %); worst 0.52 → 0.28 | the pass-0 prior should bootstrap only |
| **Honest σ²_imp(ρ)** — fixed-kernel regression of adjacent disagreement on the source's *solved* rate | 0.118 → **0.189** ❌ (worst 0.52 → **0.95**) | trained on a **not-yet-peeled** belief, neighbours "agree" on their *wrong* rates ⇒ σ² small ⇒ confident messages **propagate the error**. *Honesty measured against a wrong belief is not honesty.* |
| **V1 = the theory's own formula** `1/(π_src^belief + σ²_imp)` | 0.118 → **0.154** ❌ | gives **loud** messages to confident (stranded) sources that don't need them and **silences** uncertain sources exactly where peeling is needed |
| V2/V4/V5 (belief+count, conservative min, max) | 0.123–0.124 ❌ | `V4 = min(baseline, belief)` can only *weaken* — and it is worse ⇒ **messages are already too weak** |

The V1 result is the important one: **the theory's prescribed message precision is *worse* in production.** We
do not read this as the theory being wrong — we read it as **the theory not being adoptable piecemeal**. Its
`π_msg` presupposes the joint constrained solve *and* the sum-reconciled relay (§2a, §2d). Bolting honest
precision onto a Trojan-horsed relay makes the horse faster.

---

## 5. The DOF analysis — one bug in three costumes

### 5.1 One degree of freedom, two "independent" constraints

Take a **single-strand** node: `ρ_g + ρ_r = D`, i.e. `f_r = 1 − f_g`. **One free number.** Now count what we
apply to that one number:

**In the prior.** We multiply a **gDNA-density prior** `P(ρ_g)` by an **RNA-density Jeffreys** `p(ρ_r) ∝ 1/ρ_r`
as though `ρ_g` and `ρ_r` were independent draws. They are not — they are `f_g·D` and `(1−f_g)·D`, two views of
one number. Multiplying two independent-looking densities on one axis **double-counts the measure**. The
Jacobian `|dρ_r/df_g| = D` is what produces the `−log(1−f_g)` push; applying it *on top of* a prior that already
constrains the same axis over-counts.

**In the messages.** The node receives a gDNA message on `log f_g` **and** an RNA message on `log(1−f_g)`. Both
are added to the same ψ; both derive from the **same source belief** (the source's own `f_g`, since the source
is also one-DOF). Their precisions **add**: `pt = π_loc + π_g + π_rna`. So the messages get **~2× the weight**
of a single honest message about the single axis.

**Prediction.** If one number is being constrained twice, then everything that constrains it is ~2× too strong
⇒ **both** the message precision **and** the Jeffreys should want ×0.5. **That is exactly what we measure**
(§4.4, §4.2). We consider this the strongest evidence for the DOF thesis.

### 5.2 Two degrees of freedom (AMBIG)

An AMBIG node has three components and **2 DOF** on the 2-simplex. Here the picture differs:

* the three messages (gDNA, RNA₊, RNA₋) address a 2-dimensional space, so they are **not** trivially redundant;
* theory §8 (`bp_theory.py` TEST 6) shows the healthy behaviour: a wrong RNA₊ message is **absorbed by the other
  weak component** (RNA₋), *sparing* a confident gDNA — the joint solve routes the trade;
* but at κ≈½ RNA₊/RNA₋ are **strand-degenerate**, so the effective DOF collapses toward 1 for the gDNA-vs-RNA
  question — which is why unstranded AMBIG is where everything hurts.

So the double-count is **worst at 1 DOF** and **re-emerges at 2 DOF under strand degeneracy** — consistent with
the benchmark (ss0.99 healthy, ss0.50 pathological).

### 5.3 Why the additive-ψ relay is the mechanism

The per-node ψ solve *is* joint — but it is fed messages whose **precisions were already summed** across
components, and it emits a **running belief updated per component without sum reconciliation** (§3). So:

* on the way **in**, one axis is constrained twice;
* on the way **out**, the relayed components need not satisfy `Σρ_c = D`, so the next node receives a belief
  that is not a valid composition — the Trojan horse propagates.

This also explains the **extreme modes** of §4.3: an unreconciled component can drift to `f_g ≈ 7e-6`, and a
log-space precision-weighted average then transmits that absurdity with only 3 % of the weight.

---

## 6. Synthesis — the working thesis

> We treat **one degree of freedom as if it were two independent ones**, in both the prior (two density priors
> multiplied) and the messages (two precisions added). Everything that touches that axis is therefore ~2× too
> strong. The two independent ×0.5 optima are the symptom; the Trojan-horsed relay is the mechanism; the
> catastrophic unstranded over-call is the consequence (the Jeffreys wins the axis when the strand cannot
> object).

This is a **single** bug, not three. It predicts, and is consistent with, every measurement above:

| observation | explained by |
|---|---|
| σ²_imp wants ×0.5 | one axis, two message precisions added |
| Jeffreys wants ×0.5 | one axis, two density priors multiplied |
| 100 % of corrupted 1-DOF nodes get both messages | the double-count is the corruption path |
| stranded fine, unstranded catastrophic | the strand is the only thing that can out-vote a double-counted prior |
| the theory's honest π_msg makes things *worse* | it presupposes the joint/reconciled machinery we lack |
| stronger messages fix the worst corner but over-peel globally | messages are honest-ish; the *solve* mis-allocates |

---

## 7. What we propose (for critique)

1. **Fix the DOF treatment first** — implement theory §5/§10: one prior on the free axis (not two multiplied),
   and a **sum-reconciled** relay so a node emits a valid composition. Expect both magic ½'s to disappear.
2. **Then** de-conflate `σ²_transfer` per component (theory §6.1) — this is where the honest ~6× precision is,
   and it is only safe once (1) holds.
3. **Then** revisit `π_msg = 1/(1/π_src + σ²_transfer)`; we expect it to *win* once the relay is reconciled.
4. **Keep the prior out of the first pass** (theory §10.5) — our trace shows it is the corrupting term, not the
   messages.

We deliberately have **not** shipped either ×0.5; they are diagnostics, not fixes (`CALIBRATION_ARCHITECTURE.md`
§7 magic-number ledger).

---

## 8. Open questions — where we want input

1. **Is the double-count real, or are we over-reading two coincidental ½'s?** The cleanest falsification: if
   the DOF thesis is right, a correct one-DOF treatment should make *both* dials want **1.0**. Is there a
   sharper test?
2. **What is the correct one-DOF prior?** On a 1-DOF node, is the right object a single prior on `f_g` derived
   from `P(ρ_g)` with the *correct* change of variables to `D`-scaled composition — and does the RNA Jeffreys
   then vanish, or reappear as a Jacobian we must keep exactly once?
3. **Should message precision depend on the recipient?** Our data says *unstranded wants loud, stranded wants
   quiet*, which is a **recipient** property — yet BP says the recipient should already discount via
   precision-weighting, and the sender cannot know the recipient's state. Is "stranded degrades under loud
   messages" simply proof that the strand's precision is **understated**, or is a recipient-aware gate
   legitimate?
4. **How do we estimate `σ²_transfer` without circularity?** Every belief-free per-node gDNA rate is either
   RNA-contaminated (raw density) or prior-circular (projection); the honest one needs a solved belief, which
   needs messages. Theory §6.1 proposes gDNA-identifiable edges (RNA-free enriched exons). Is that enough
   support at genome scale, and does it generalise off-capture?
5. **Extreme modes:** should a message mode be bounded (theory says "verify-don't-clip"), or does
   sum-reconciliation alone remove the pathology?
6. **AMBIG at κ≈½:** when RNA₊/RNA₋ are strand-degenerate, is the honest statement "1.x DOF"? Should the solve
   explicitly collapse the degenerate pair rather than carry two near-unidentifiable components?

---

## 9. Reproducing this

```bash
conda activate rigel
SC=<scratch>;  SUITE=~/Downloads/rigel_runs/ambig_dense_10mb

# one-time: cache the per-condition scan + oracle (calibration-independent)   [~13 min]
python scripts/debug/calib_pool_benchmark.py --suite $SUITE --cache-dir $SC/bench_cache --out $SC/base.tsv
# thereafter: the full 24-condition oracle benchmark                          [~9 s]
python scripts/debug/calib_pool_benchmark.py --suite $SUITE --cache-dir $SC/bench_cache --out $SC/x.tsv

python scripts/debug/msg_precision_diag.py --cache-dir $SC/bench_cache          # σ²_imp + its precision cap
python scripts/debug/msg_trace.py --cache-dir $SC/bench_cache --condition <c>   # per-message corruption trace
```

| tool | what it answers |
|---|---|
| `scripts/debug/calib_pool_benchmark.py` | calibration-basis error vs oracle, per condition (`--cache-dir` ⇒ 9 s loop) |
| `scripts/debug/msg_trace.py` | strand-only vs local vs final; corrupted/rescued nodes; local-vs-message precision; the DOF double-count check |
| `scripts/debug/msg_precision_diag.py` | σ²_imp as production computes it, and its precision cap |
| `src/rigel/calibration/message_precision.py` | the honest σ²_imp(ρ) estimator — **built, ablated, currently UNUSED** (§4.5) |

Relevant code: `bp_solver.node_sweep` / `_scan` (the relay), `simplex_logodds._solve_nodes_logodds_all` (the
per-node ψ solve), `gdna_rate_prior.GdnaRatePrior.logprior` (the prior projection + the Jeffreys),
`bp_solver.adjacent_disagreement_variance` (σ²_imp).

---

## 10. Status ledger

| item | state |
|---|---|
| pass-0 gDNA-rate NPMLE prior (replaces KDE + floor) | **shipped** (`npmle_struggles.md`) |
| refit loop (`calib_refit_iters=1`) | **shipped**, −22 % mean error |
| 9-second oracle dev loop | **shipped** |
| honest per-node σ²_imp(ρ) | built, **ablated, not shipped** (backfires) |
| Jeffreys ×0.5 / σ²_imp ×0.5 | **diagnostics only — not shipped** (magic numbers) |
| joint constrained solve + sum-reconciled relay (theory §5/§10) | **not implemented — the proposed next step** |
| per-component `σ²_transfer` (theory §6.1) | **not implemented** |
