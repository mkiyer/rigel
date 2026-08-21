# The spike-and-slab reference — a primer with a worked example

    ⚠ **A DEV DOC, written 2026-08-20 as teaching material for the owner.** The authoritative
    derivation is `EQUATIONS.md` §9d (the digits) and the ruling is `DESIGN.md` §0c.3 (the shape);
    where this primer and those disagree, they win. Delete this file once it has served its purpose.

## 1. The problem, in one sentence

ψ solves one object's composition from its own counts **plus a prior guess** — the "reference", a
Beta(a, b) whose mean `m` says *"before looking at this object's evidence, what share of its
unspliced mass do I expect to be gDNA?"* — and today that guess is a **constant** (½ at any slot
mature RNA can reach; 0.75 where it cannot), which is wrong at almost every object in a real
library.

Why this is the thing to fix (measured this session, 16-condition ladder, muted rung, certified
truth): ~90 % of the in-scope calibration frontier is a diffuse under-call at exons and their
boundaries, and feeding the reference the condition's true mean melts the near-vertex population
0.029× at g98. The reference is the root cause of the frontier — but no single number can be the
fix, because the same truth-fed mean reads **2.65× WORSE** at `g05 stranded × capture-ON`. One
mean cannot sit where every object's truth sits. **The reference must be per-object.**

## 2. One decomposition carries the whole question

The per-object mean is `m_i = ρ_g,i · E_g,i / M_i` (§9c): the gDNA the object's own density
predicts, as a share of the unspliced mass it actually holds. `E_g,i` (gDNA opportunity) and
`M_i` (unspliced mass) are exactly known per object. **The only unknown is `ρ_g,i` — the gDNA
density at THIS object.** Decompose it:

    ρ_g,i = ρ_0 · ε_i        ρ_0 = the un-enriched (background) gDNA density
                             ε_i = the capture enrichment at object i, ε ≥ 1

Off capture, `ε ≡ 1` everywhere and `ρ_0` is directly measurable from the objects mature RNA
cannot occupy (intergenic, introns, and the boundaries against them) — that case is solved: the
per-object form beats the shipped constant 0.40–0.50× on contaminated capture-OFF rows and takes
both zero controls to exactly 0.

Under capture, everything hinges on `ε_i`, and the measured field says how hard that is
(§9d.2, g98 ss0.99, fragments/bp):

| | intergenic | intron | exon\|intergenic | exon\|intron | span |
|---|---|---|---|---|---|
| capture OFF | 0.1005 | 0.1005 | 0.0983 | 0.0984 | 1.0× |
| capture ON | 0.0042 | 0.0042 | 0.510 | 0.478 | **122×** |

## 3. Why the distribution of ε is a SPIKE plus a SLAB

A probe panel is a **finite list of targets**. An exon is either in the panel or it is not:

* **Unprobed ⇒ ε = 1 exactly.** Not "small" — *exactly* one: the same molecules, no pull-down.
  Across many exons that is a **point mass** (the SPIKE).
* **Probed ⇒ ε is a large factor** set by probe count, overlap, and chemistry — things the object
  itself says almost nothing about. Across many exons that is a **broad distribution** (the SLAB).

So the prior on `log ρ_g,i` is a two-component mixture — and this is the *physical structure of a
capture panel*, not a modelling convenience:

    p(log ρ_g,i) = π · N(log ρ_0, σ_0²)  +  (1 − π) · Unif[ log ρ_anchor,i , log ρ_max,i ]

Every term is measured or structural — **no new constant** (§9d.4's table):

| term | source |
|---|---|
| `ρ_0`, `σ_0²` | the off-target anchors (intergenic + intron REGIONs) — measured exactly, before any solve |
| slab lower `ρ_anchor,i` | the adjacent `exon\|intron` flank's density — enrichment is **monotone in probe proximity**, so the flank (at the footprint's edge) is a floor for the exon inside it |
| slab upper `ρ_max,i` | the object's **own total density less the RNA the adjacent junction proves is there** (a spliced fragment cannot be gDNA); degenerates to "an object cannot hold more gDNA than the mass it holds" where there is no flux |
| slab shape | uniform in `log ρ` — the scale-invariant ignorance statement for a positive rate on a known interval |
| `π` | the unprobed fraction. At pass-0: the Jeffreys ½ (the same convention that already sets the strength). And it is **measurable**: the off-target/in-gene anchor ratio is 0.98 without probes and 113–114 with, so capture *detection* is free; with a supplied panel list, π is **counted** |

Two properties give the mixture its power:

* **A Beta cannot put its median at a vertex; a mixture with an atom can.** The near-vertex
  population (true `f_g ≥ 0.999`) carries 60–65 % of the frontier at g98, and only an atom
  reaches it honestly.
* **A location near a vertex is a STRONG claim regardless of the strength budget.** With
  `a + b = 1` but `m → 0`, the term `b·log(1 − f)` still diverges against the opposite vertex —
  logit(m) *is* the claim's size in nats (DESIGN §6b). This is why every single-mode per-object
  candidate failed under capture (§6 below): pricing an *enriched* slot with the *depleted*
  background puts a near-vertex location on the wrong side, and that is a loud wrong claim, not
  a gentle one.

## 4. The two bounds are a SUPPORT, not an estimate

    lower:  f_g ≥ ρ_anchor,i · E_g,i / M_i                    (monotone flank)
    upper:  f_g ≤ 1 − (ρ_r·E_r,i − S) / M_i                   (junction flux; S = certified spliced)

Both come from **neighbours** — which is §0c's structural fact in arithmetic: an exon's gDNA is
not measurable at the exon, at any depth. Measured violation mass for the upper bound is 0–19.5 %
across conditions (`TRAPS: an-upper-bound-is-not-an-estimate` — reading an endpoint as the answer
is the refused move; a bound with a measured violation rate is *a support with a tail, which is
exactly what a slab is*). The bound is loose where RNA is abundant and tight where it is scarce —
the right way round: at `g98 ss0.99 ON` the upper endpoint alone reads 0.9918 against a truth of
0.9817, already excluding the neutral ½.

## 5. The worked example

Library: `g98 ss0.99 capture-ON` (the worst in-scope condition; anchor table in §2). Background
`ρ_0 = 0.0042`; in-gene flank anchor `0.478`.

**Object A — a probed exon.** `E_g = 2,000 bp`, unspliced mass `M = 3,000` fragments, true
`f_g = 0.98` (so its true `ρ_g = 0.98·3000/2000 = 1.47` — the flank's 0.478 under-reads it 3.1×,
inside the documented 2.6–3.6× band). The adjacent junction's flux proves ~60 unspliced-RNA
fragments.

| candidate reference | m_i | in nats (logit) | what happens |
|---|---|---|---|
| shipped constant | 0.500 | 0.0 | holds the slot 2–3 % off the vertex; tiny per slot × every exon = **the P1 mechanism, 60–65 % of the frontier** |
| single-mode, background density (`pooled`) | `0.0042·2000/3000` = **0.0028** | **−5.9** | a near-vertex location on the WRONG side — measured 5–8× regressions, `B exon\|exon` crushed 11× |
| single-mode, flank density (`local`) | `0.478·2000/3000` = 0.319 | −0.76 | right order for THIS slot — but a corridor slot with no `exon\|intron` flank falls back to the row above, and that fallback is where the damage was (identical to 0.04 % between the two variants) |
| **the mixture** | spike at f ≈ 0.0028 (π = ½) + slab on [0.319, `(3000−60)/3000` = **0.980**] | median at π = ½ → **0.319**, honestly weak | the support *contains the truth* (0.98 ≤ upper endpoint); the location is mild (−0.76 nats, not −5.9); the slot's own 3,000 fragments dominate, as they should |

**Object B — an unprobed exon, same library.** Its flank anchor also reads the depleted level, so
slab-lower ≈ spike ≈ 0.0042·E/M: the mixture *collapses toward the spike* and correctly says
"almost no gDNA here" — the SAME formula, adapted per object by its own measured endpoints.

**The three limits (each checked in §9d.4):**

1. **capture-OFF** — anchor ratio 0.98 ⇒ no enriched population exists ⇒ the slab collapses onto
   the spike ⇒ exactly the validated capture-OFF per-object form. A strict generalisation, by
   construction rather than tuning.
2. **g00 (zero-gDNA control)** — the anchors measure `ρ_0 = 0`, the flank bound is 0, spike and
   slab-lower coincide at 0, and the median reaches the vertex **because a measurement put the
   atom there**. (One care: `log ρ` breaks at 0 — take the slab in `ρ` or clamp, §9c.1's clamp.)
3. **g98 capture-ON** — the slab's upper endpoint alone (0.9918 vs truth 0.9817) is informative
   before any mixing.

What ψ finally receives is the mixture's **median** as its location, at the **unchanged strength
of one pseudo-observation** (§9c.1 — overturned by ~1.5 stranded fragments; nothing about
strength is re-opened).

## 6. This session's refutations, as design constraints (all at the muted rung, current tree)

1. **Library-wide mean**: melts g98's vertex population 0.029× but reads 2.65× at g05
   stranded × ON, and at `g98 ss0.50 off` it is a *cancelling pair inside a 1.001× total*
   (P1 0.029× vs mid-truth 2.574×). ⇒ per-object is REQUIRED.
2. **Single-mode per-object (`pooled`, and the previously unpriced `local`)**: 5–8× on all three
   stranded × capture-ON conditions, with `B exon|exon` damage identical to 0.04 % between them —
   deep corridor slots have no `exon|intron` flank, fall back to the depleted background, and get
   a loud wrong near-vertex location. ⇒ a point estimate of `ρ_g` cannot serve both enrichment
   populations; **the mixture is not optional under capture**.
3. The corridor population (no flank, no pure-gDNA object in reach — S3c's shape, the hop
   census's "neither currency serves `R exon ← B gene edge[term]` under capture") is exactly the
   population whose **enrichment state** the mixture must infer. That inference is the open
   design question below.

## 7. Where you can contribute — the genuinely open questions

* **What updates π per object at pass-0?** Candidates, cheapest first: (a) leave it at ½ and let
  the slot's own evidence decide (the mixture's medians already adapt through the endpoints);
  (b) the slot's **total** abundance projected onto the bimodal total-abundance landscape — your
  own robustification idea: capture enriches gDNA and RNA alike, so the *total* is informative
  about ε without touching the composition being solved; (c) where a panel's target list is
  supplied, π is counted, not inferred. (b) is attractive and needs a derivation for its
  variance; it is the same object as your "project local disagreement onto the global landscape"
  note.
* **The second solve**: once the gDNA landscape is trained on pass-0's solved exons, does it
  *replace* the slab or *multiply* the mixture? (The circularity discipline from §0c.0d applies:
  too strong a reference and the landscape learns the prior back.)
* **The strength on structural slots** — the one degree of freedom every measurement so far has
  held fixed (`a + b = 1`); the hard-vs-soft-floor record in `vertex_ceiling.py` says choosing it
  by panel number alone would be tuning.
* **σ_0² and the slab's tails**: the upper bound has 4.6–19.5 % measured violation mass, so the
  slab needs a tail policy (clamp vs widen), and the g00 coordinate care.

## 8. The gates it must pass (already written, rank 3)

Reduce EXACTLY to the shipped form at capture-OFF (limit ① is the test); both zero controls; a
SHUFFLE control; and — from this session — the three refutations above re-run as arms that must
now come out the other way, plus the evidence split on every row
(`TRAPS: an-imputation-must-cost-something-every-hop`'s discipline applied to a prior).
