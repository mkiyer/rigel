<!-- title: AMBIG zero-gDNA phantom — the no-neighbour intrinsic-solve fix -->
# Fix (B): the zero-gDNA AMBIG phantom is a sweep-structure bug, not a precision bug

**Status:** derived + prototyped + measured; awaiting review before shipping.
**One line:** strand-decisive AMBIG region nodes whose *both* chain neighbours are empty are **never solved** —
the directional sweep only calls `_solve` when a neighbour message exists — so they keep their signature
init `{0,0,1}` (f_g = 1.0) and manufacture phantom gDNA. The fix solves every solvable node on its own
intrinsic evidence (strand likelihood + global prior) even when it has no neighbour.

---

## 1. The problem

In a **stranded, zero-gDNA** library (the truth is *no gDNA anywhere*) the calibration still reports
**~35,400 fragments of gDNA on AMBIG exon nodes** — pure phantom. This is the long-standing zero-gDNA
AMBIG over-call. It also drives ~half of the analogous unstranded-zero phantom.

I had attributed this to the **global gDNA prior being too weak at AMBIG nodes** (a precision problem —
that the variance curve, queried at the RNA-inflated observed density `ρ_obs`, went humble exactly where we
need to pin gDNA to zero). **That diagnosis was wrong.** The data below shows the global prior at these
nodes is in fact *strong and correct*; the phantom has nothing to do with precision.

## 2. The investigation that overturned the precision story

I tracked, per sweep pass, the estimated global gDNA-density mean against an **oracle** (built by scanning a
by-origin split of the same BAM through the production accumulator), and the per-node global precision.

**The global mean is right.** Estimated vs oracle, mass-weighted over the self-solvable nodes:

| scenario | oracle ρ̄ | estimated ρ̄ | verdict |
|---|---|---|---|
| capture (gdna300, ss0.99) | 13.42 | 13.48 (stable across passes) | ✅ right, converges |
| stranded-0 (none, ss0.99) | 0.00 | 0.0001 | ✅ right, converges |
| unstranded-0 (none, ss0.50) | 0.00 | 1.96 → 3.20 | ❌ diverges (a *separate* bug — §6) |

**The global precision at the phantom nodes is huge, and points the right way.** Dumping the per-node global
term at the high-mass AMBIG nodes in stranded-0:

```
node463:  ρ_obs=12.2  σ²_n=5.5e-7  τ_global=2.7e8  μ_global≈0   → f_g=0.000   ✅
node453:  ρ_obs=12.5  σ²_n=5.7e-7  τ_global=2.7e8  μ_global≈0   → f_g=0.000   ✅
node473:  ρ_obs=12.3  σ²_n=5.6e-7  τ_global=2.7e8  μ_global≈0   → f_g=1.000   ✗  phantom
node2167: ρ_obs=13.1  σ²_n=6.0e-7  τ_global=2.9e8  μ_global≈0   → f_g=1.000   ✗  phantom
```

All four nodes carry an **identical** strong, correct global pull toward f_g = 0 (τ ≈ 2.7e8). Two obey it;
two ignore it and sit pinned at the init f_g = 1.0. A correct, strong prior that some nodes simply *ignore*
is not a precision problem — something is preventing those two nodes from ever being solved.

(My earlier "AMBIG τ is weak" reading was a **median artifact**: the median over *all* AMBIG nodes was 1.0
because it was dominated by zero-*mass* AMBIG nodes floored to τ=1; the AMBIG nodes that actually *hold the
phantom mass* have τ ≈ 1e8.)

## 3. Root cause

The two stuck nodes are not balanced-count nodes — their strand counts are **decisive**:

```
node473:  upos=8200,  uneg=71     → overwhelmingly one strand ⇒ strand likelihood wants f_g≈0
node2167: upos=256,   uneg=22945  → overwhelmingly the other strand ⇒ strand likelihood wants f_g≈0
```

So the strand likelihood, *if consulted*, would pin them to ≈0 — exactly as it does for the structurally
identical node463 (upos=16336, uneg=149 → f_g=0.000). The difference is purely their neighbours:

```
node : left  right_face_mass(L)   right  left_face_mass(R)   has message neighbour?
 473    472        0.00            474         0.00            NO   ← never solved
2167   2166        0.00           2168         0.00            NO   ← never solved
 463    462        0.00            464       233.49            YES  ← solved (→ f_g=0)
 453    452        0.00           454        58.22            YES  ← solved (→ f_g=0)
```

The sweep is a **directional Gauss-Seidel chain walk**: in the L→R pass each node is pulled by its *left*
neighbour's message, in R→L by its *right* neighbour's. The inner loop is

```python
for i in ...:
    if not solvable[i]:
        continue
    src = nbr_arr[i]
    if src < 0 or MS[sf][src] <= _EPS:
        continue          # ← reference terminal / empty neighbour: SKIP
    ... build message ...
    _solve(i, message..., global...)   # ← _solve lives INSIDE the message branch
```

`_solve` — which integrates the strand likelihood, the global prior, the Jeffreys prior **and** the message
— is only ever reached *after* a valid neighbour message is constructed. A node whose both neighbours are
empty hits the `continue` in **both** directions, so `_solve` is **never called for it**. It is left at its
`init_beliefs` value, and the signature-binary init for an AMBIG node is `{f₊,f₋,f_g} = {0,0,1}` — i.e.
**f_g = 1.0**. The node's own strand counts (8200:71) and the correct global prior (μ≈0, τ≈1e8) are never
looked at.

This is the bug: **`_solve` is gated on the existence of a neighbour message, but the strand likelihood and
the global prior are node-*intrinsic* and need no neighbour.** An isolated node (both neighbours empty) is
the one case where the gate silently drops a node on the floor at its gDNA-favouring init.

Why "isolated" is common in zero-gDNA — and the reframing it forces: a boundary node only carries mass when
fragments *cross* it (unspliced fragments spanning the seam). In a clean RNA library the only fragments are
RNA, which is **spliced** inside exons and does not cross exon↔intron seams as unspliced mass; gDNA, which is
the thing that *would* cross seams genomically, is absent. So nearly every boundary is empty. Measured on the
benchmark chains (1221 region + 1222 boundary nodes):

| scenario | empty boundaries | message-isolated solvable nodes | their share of solvable mass |
|---|---|---|---|
| capture (gdna300) | 13 / 1222 (1%) | 22 / 1924 (1.1%) | **0.0%** |
| stranded-0 (none) | 1195 / 1222 (**98%**) | 275 / 325 (**85%**) | **80%** |

This is the important point: in the zero-gDNA regime the message graph is **almost entirely disconnected** —
80% of the solvable mass sits in nodes with *no* neighbour to hear from. So an intrinsic self-solve is not an
edge case there; it is the **dominant** path. Calibration of a clean library is mostly a per-node
strand+global solve, with propagation only kicking in where unspliced mass actually crosses seams (i.e. where
gDNA is present). The same AMBIG exon under capture has crossing gDNA on its boundaries → it gets solved by a
message → the bug is invisible, which is exactly why this surfaced only on the zero-gDNA conditions.

## 4. The fix

Solve every solvable node on its **intrinsic evidence** at least once per pass, regardless of neighbour
messages. Concretely: after the two directional passes, sweep the nodes that received no message this pass
(both neighbours empty) and `_solve` them with **zero messages** — strand likelihood + global prior +
Jeffreys only:

```python
# precomputed once (masses are frozen across passes):
lv = (left  >= 0) & (MS[1][left]  > _EPS)   # left neighbour's right face has mass
rv = (right >= 0) & (MS[0][right] > _EPS)   # right neighbour's left face has mass
has_msg_nbr = lv | rv

# after the L→R and R→L directional passes, inside the per-pass loop:
for i in order:
    if solvable[i] and not has_msg_nbr[i]:
        _solve(i, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, mu_global[i], tau_global[i])
```

**Why this formulation (and not the alternatives):**

- **It respects the directional / "never blend two boundaries onto one region" ruling.** The intrinsic solve
  fires *only* for nodes with no neighbour at all, and passes **zero** messages — there is nothing to blend.
  Nodes with ≥1 neighbour are untouched here; the directional loop already solves them exactly as before.
- **No double-solve, no overwrite.** `has_msg_nbr` is the exact complement of "solved by the directional
  loop": a node with even one non-empty neighbour is solved in that direction (the empty direction's
  `continue` then correctly preserves that result). So the intrinsic pass adds solves *only* for the nodes
  that would otherwise get none.
- **It fixes the init, not the prior.** The architecture deliberately makes `init_beliefs` signature-binary
  and count-blind (AMBIG → `{0,0,1}`) and delegates the count/strand refinement to the *solve*. The bug is
  that the solve never runs; the fix restores the solve. We do **not** make the init count-aware (that would
  duplicate the strand machinery in the wrong place).
- **Convergence is safe.** An isolated node is no one's source (its empty-mass neighbours are themselves
  unsolvable) and receives no message, so its intrinsic solve depends only on frozen strand counts and the
  per-pass global prior. It cannot feed back; it simply tracks the (converging) global each pass.

## 5. Measured impact

AMBIG-node estimated gDNA mass `Σ f_g·M` (oracle truth: 0 for both zero-gDNA cases; ~real for capture):

| scenario | baseline | with fix | note |
|---|---|---|---|
| **stranded-0** (none, ss0.99) | 35,400 | **0** | phantom eliminated — strand pins 8200:71 → f_g=0 |
| **unstranded-0** (none, ss0.50) | 64,038 | 34,886 (**−46%**) | residual is the separate Bug C (§6) |
| **capture** (gdna300, ss0.99) | 360,298 | 359,766 (−0.15%) | unchanged — this AMBIG mass is *real* (oracle f_g=0.84) |

Confirmed on the production-facing `CalibrationResult` (clean, de-gated code): AMBIG-region contained gDNA
mass in stranded-0 is **0** (was 35,400; oracle 0), total contained gDNA 24 ≈ 0; capture and unstranded
match the table above. Test suite: **full suite 1050/1050 pass** with the fix (incl. the factor-1-uniform
invariant and all golden regressions — no golden shifts: the existing golden scenarios contain no isolated
no-neighbour AMBIG node, so the fix is a strict no-op wherever every node already has a neighbour).

## 6. What it does NOT fix (honest scope)

- **The unstranded-zero divergence (Bug C).** At κ≈0.5 the strand likelihood is flat, so the intrinsic solve
  for an isolated node falls back on the global prior alone — and in the unstranded case the global *mean*
  itself diverges upward (1.96 → 3.20 vs truth 0), because the "self-solvable" single-strand nodes that feed
  the mean are corrupted by the Jeffreys Beta(½,½) vertex-push when the strand is flat. The fix still halves
  this phantom (the isolated nodes are now pulled toward the wrong-but-bounded mean instead of stuck at 1.0),
  but the root divergence is a distinct bug (the Jeffreys readout under a flat likelihood) and is tracked
  separately.
- **The capture AMBIG residual.** Under capture the AMBIG mass is governed by neighbour messages with a
  (correctly) weak global; the small remaining error there is message-side and is the target of the
  count-space pseudo-count pivot, not this fix.

## 7. Relation to the belief-propagation relay model (the deeper question)

A clean way to state the intended BP semantics: a message has a **density** (signal) and a **precision**
(reliability); a node integrates incoming messages weighted by precision and forwards its result.

- A node with **zero internal precision** on a channel (no data — e.g. the f_g of a zero-count node, or the
  f₊ of an un-sampled exon) should **trust and relay** an incoming message on that channel ~unmodified.
- A node with **infinite internal precision** on a channel (a structural lock — e.g. f₋ = 0 on a +-strand
  region) should **ignore** the incoming message and forward **zero** on that channel.

Where do we stand against this?

- **The "infinite precision ⇒ block" half is implemented**, via the `free_s` strand gate. A +-strand region
  has `free_neg = False`, so the antisense message is never formed in or out of it — it kills the f₋ channel
  and forwards zero, exactly the locked-`f₋=0` case. ✓
- **The "zero precision ⇒ trust + relay" half is NOT implemented**, for two verified reasons:
  1. The solve gate `solvable = (fp|fn) & mass>0` skips zero-count nodes as a *destination* — so they never
     adopt an incoming message into their belief in the first place.
  2. The message precision is **re-derived per hop from the immediate sender's mass**, and the binomial floor
     `var_floor = f(1−f)·(1/M_src + 1/M_dst)` sends `τ→0` as `M_src→0`. So even ungated, a zero-count node
     emits a zero-precision message — it cannot **forward** a confident upstream signal. Our messages are
     *re-derived*, not *forwarded*; a low-count node is a wall, not a transparent conduit.

  This is a genuine architectural inconsistency, but its **measured footprint is ~zero** (§3 table): under
  capture (gDNA present, where relay *would* matter) the chain is 99% connected, so the per-hop mechanism
  relays fine; under zero-gDNA (chain 98% fragmented) there is no gDNA to relay. The regime that fragments
  the chain is the regime with nothing to forward. True forward-relay would matter only for **sparse-but-
  present** gDNA (gDNA regions separated by empty gaps) — not stressed by the current suite. Note this is
  **orthogonal to the count-space pivot**: a count-space message `α = N_pseudo·f[src]` is *also* re-derived
  from the sender's own `N_pseudo`, so the pivot does not add forward-relay either. Forward-relay is a third,
  independent axis, deferred until a scenario demonstrates it matters.

## 8. Why this is the right amount of gating (not too little, not too much)

The fix is the **minimal complete** realization of "every node solves itself": a message-solve already
includes the intrinsic terms (strand + global + Jeffreys), so an extra unconditional intrinsic solve would be
*redundant* for any node that has a message — it only ever changes the outcome for a node with **no** message,
which is exactly `not has_msg_nbr`. We should **not** drop the remaining gates:

- `mass > 0`: a zero-count node's own composition is moot (no fragments to classify in the EM), so solving it
  is wasted work — and it can't usefully relay anyway (§7). Leave it unsolved at its init.
- `(fp|fn)` (excludes intergenic / NONE): an intergenic node's signature init `{0,0,1}` (≈ all gDNA) *is* the
  correct gDNA-prior answer; "solving" it against the global *fraction* (< 1 under capture) would wrongly pull
  its f_g down. The gate keeps it pinned as the G1 sink it is. ✓

So the answer to "do we need extra gating if nodes solve elegantly?" is **no extra gating, and no fewer
gates** — the existing gates are correct, and the one missing piece was the intrinsic solve for the
no-message case, which this fix adds.

## 9. The diff

Two hunks in `calibration/bp_solver.py::node_sweep`, ~7 lines: (1) precompute `has_msg_nbr` next to the
existing per-node global arrays; (2) the intrinsic-solve loop after the two directional passes. No new
constants, no signature changes, no C++.
