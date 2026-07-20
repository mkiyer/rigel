# The Spliced Message Precision (`I_motif`) — a belief-propagation synopsis for external review

**Purpose.** A self-contained statement of *what the spliced-junction message is* in Rigel's calibration
belief-propagation (BP) pass, written for external mathematics / statistics reviewers. The narrow question we
need settled: **is the spliced-fragment message a legitimate BP message, and what precision does it carry?**
An internal adversarial review rejected an earlier formulation on "short-exon" grounds; we believe that
rejection mis-modelled a BP *message* as a hard *constraint*, and we want an authoritative outside read before
we implement.

This document assumes no Rigel-specific background beyond §1.

---

## 1. The setting (minimal)

Rigel deconvolves a hybrid RNA-seq + genomic-DNA (gDNA) library. The genome is partitioned into a **linear
chain of nodes** — alternating **regions** (segments between transcript features) and **boundaries** (splice
junctions and gene edges): `… region ↔ boundary ↔ region ↔ …`. Each node owns a fragment count and we solve,
per node, the composition of its **unspliced** fragment mass on the axis

```
    λ = logit(f_g),   f_g = P(a node's unspliced fragment is gDNA),   f_r = 1 − f_g = P(RNA).
```

The solve is a single forward-backward **belief-propagation pass** over the chain: each node holds a Gaussian
belief `N(μ_λ, σ²_λ)`, receives one message from each neighbour, folds them, and relays. **Everything below is
about one kind of message.** There is no other machinery in play.

**The identification problem (why RNA evidence is scarce).** The only *intrinsic* per-node signal for `f_g` is
the **genome-strand tilt**: RNA is transcribed from one strand (library sense fraction κ), gDNA is
strand-symmetric (½). The strand Fisher information for `λ` is `∝ N·(2κ−1)²`, which is **identically 0 when the
library is unstranded (κ = ½)** — the regime we must handle. In an unstranded library, then:

- **gDNA evidence exists**: intergenic regions are structurally pure gDNA (a composition certainty).
- **RNA evidence is almost absent** — *except* at splice junctions.

A **spliced** fragment (one that aligns across an exon-exon junction with the intron removed) is **guaranteed
pure RNA**: only a processed mRNA can be spliced; gDNA cannot. So **spliced fragments at a boundary are the
principal — often the only — source of RNA evidence in an unstranded library.** Without this evidence, the pass
has gDNA sources and no RNA counterweight, and composition propagates toward gDNA unchecked. Supplying that RNA
evidence as a BP message is the subject of this note.

---

## 2. What the boundary observes (all directly measured, no inference)

A boundary node `b` that is a splice junction directly measures, on its own effective length:

| quantity | meaning | how obtained |
|---|---|---|
| `S_b` | spliced fragment count (or fractional mass) crossing the junction | the accumulator, at deposit |
| `E_spl,b` | the boundary's **spliced effective length** (the genomic measure of start-positions that produce a spliced crossing) | a pure geometric function of the fragment-length pmf; **always present and computable** |
| `strand_b` | the junction's motif strand (+/−) | the splice-site dinucleotide motif, at deposit |

From these, the boundary computes a **spliced-RNA density** and its **count-based precision**:

```
    ρ_spl,b  = S_b / E_spl,b                       (a measured pure-RNA density)
    Var(log ρ_spl,b) = 1 / S_eff,b                 (Poisson / count precision;  S_eff = Kish effective count)
```

`S_eff` is the fragment count's statistical power — the honest amount of information in "there are `S` spliced
reads here." Nothing here is inferred; it is a direct measurement with a known variance.

**Owner's principle (the crux of our position).** *The boundary sends the message that it observes spliced
fragments. It does not need to know the recipient's exon length, coverage, or configuration — it reports what
it measured and lets the recipient handle the message.* This is ordinary BP: a factor node emits its local
evidence; the variable node it speaks to does the reconciliation.

---

## 3. The message (the proposal)

`I_motif` is the message a spliced boundary `b` emits **to its transcript-continuous exon neighbour** `e`
(single-stranded: routed only to the flank carrying `strand_b`'s exon; one-directional along the junction).
Following §2, it is a Gaussian factor on `e`'s RNA coordinate:

```
    m_{b→e}(λ_e)  ∝  exp[ −½ · pr_motif · ( c_e(λ_e) − mode )² ],
      mode        = log ρ_spl,b            (the measured spliced-RNA log-density),
      pr_motif    = 1 / ( 1/S_eff,b  +  σ²_transfer(b→e) ),
```

where `c_e(·)` is `e`'s RNA-density coordinate (`log f_r,e + const`) and `σ²_transfer` is the belief-free
enrichment-crossing transfer variance already used by every message in the pass (it damps a message that
crosses a capture-enrichment discontinuity; ≈ 0 within a uniform regime). **The precision is the spliced
count's power `S_eff`, attenuated only by the honest transfer noise** — exactly the two things the owner named:
Poisson counting over the boundary's (short) support, and the boundary→exon transfer variance.

`I_motif` enters **only** the RNA channel (a spliced read is evidence of RNA presence, never of gDNA), and it
is a **standing additive factor** seeded from the *raw* measured `S_b` — so no downstream subtraction (the
mature-vs-nascent absorption at the junction) can cancel it. `S_b = 0 ⇒ pr_motif = 0` (no spliced, no message —
count-zero-information preserved).

This is a **message, not a constraint.** It asserts "RNA is present here at density `ρ_spl`, with `S_eff`
counts of support." The recipient folds it against its own belief and its own data, and *combines* — it is not
pinned.

---

## 4. How the recipient handles it — including the tiny-exon case (the rebuttal to the rejection)

BP is explicit about what a node does with an incoming message: it multiplies (in log-space, adds) the message
into its belief, folds, and relays the product (minus the return edge) onward.

Two regimes for the recipient exon `e`:

**(a) `e` carries its own unspliced fragments (a normal or long exon).** `e` has its own likelihood over
`λ_e`; it *combines* the spliced message with that likelihood by precision-weighting. The spliced message
raises `e`'s RNA evidence; `e`'s own data tempers it. Standard BP reconciliation.

**(b) `e` is a tiny exon with ~no unspliced fragments (the case the rejection worried about).** A short internal
exon (length ≪ fragment length) contributes **~no unspliced fragments of its own** — its mature reads all cross
the junction and are counted as *spliced*, not as its unspliced mass. In BP, **a node with no local likelihood
passes the incoming message through unchanged**: its cavity toward the next node *is* the incoming message. So
the spliced RNA message **relays onward** through the tiny exon to the next boundary/region — to wherever there
*are* unspliced RNA fragments for it to inform. The tiny exon is a relay hop, not a sink. **This relay is the
mechanism already built and in production** (the forward-backward chain relay); it is not new machinery.

**Why we believe the rejection is mis-framed.** The internal adversarial review modelled `I_motif` as a
one-sided *hard bound* `f_g,e ≤ 1 − ρ_spl/ρ_tot,e` and showed that, applied as a hard cap to a short exon's
unspliced belief, it drives a nearly-fragment-free exon toward RNA. But `I_motif` is a **soft Gaussian message
with finite precision**, folded — not a cap. On a node with negligible own mass, a finite-precision message
does not "pin" a wrong answer; it either (a) relays onward (regime b), or (b) is out-weighed once the node
*does* have data (regime a). The concrete failing examples in the review were short exons with `mat_uns = 0`
(no unspliced mass) — precisely the nodes that, in BP, have no likelihood and therefore **relay** rather than
absorb. We do not think the message needs a geometry-dependent down-weight to be honest; the relay already
supplies the length-awareness for free.

---

## 5. The precise question(s) we need reviewers to settle

We want the mechanism of §3–§4 confirmed or corrected on these specific points. Each is stated so an outside
reviewer can adjudicate without Rigel internals.

**Q1 — Pool legitimacy.** The spliced fragments (`S`, mature RNA that *left* via the junction) are a *disjoint
pool* from the exon's *unspliced* crossing fragments (nascent RNA + gDNA), whose composition `f_g,e` is what we
solve. Is it legitimate for a message derived from the spliced pool to inform `f_g,e`? Our position: yes,
because the spliced pool is **direct proof that this locus is transcribed**, and the message only asserts *RNA
presence* (a soft push on the RNA channel), relayed to where unspliced RNA exists — it does **not** claim to
quantify the exon's nascent fraction (which is separately un-inferable: there is no fixed ratio of nascent to
mature). Is "RNA presence, soft, relayed" a legitimate BP message, or does the pool disjointness forbid it
even in soft form?

**Q2 — The relay through a fragment-free node.** In our BP implementation the emitted message precision carries
a `1/N_src` sampling term (the emitting node's own count). For a *relay* through a fragment-free tiny exon, is
the correct behaviour that the message passes **undamped** by that node's (absent) count — i.e. the `1/N`
term must not re-apply at a pure relay — and if so, what is the exact cavity expression that guarantees a
tiny exon relays the upstream spliced evidence without attenuation? (This is the one place our current message
formula may under-propagate the RNA evidence, and we want it derived rather than patched.)

**Q3 — The precision value.** Is `pr_motif = 1/(1/S_eff + σ²_transfer)` the correct precision for a
count-`S_eff` density measurement relayed one hop, or should the boundary's short spliced effective length
enter the precision (not only the density value)? We currently put the eff-length only in the density `ρ_spl =
S/E_spl` (the mode) and the count `S_eff` in the precision. Is that the right split?

**Q4 — Coordinate / frame.** The message is a factor on the RNA log-density; the running belief is on
`λ = logit f_g`. We recently found (and fixed) that accumulating a message's precision requires a Jacobian
conversion between the `log f_c` frame and the `λ` frame — without it the chain relay self-amplifies by
`1/(1−f_g)²` per hop (a geometric bootstrap that saturates at `1/σ²_transfer`). We want confirmation that the
spliced message's precision, once folded and relayed, is frame-consistent under this correction.

---

## 6. Why this must resolve to a *solution* (not a rejection)

The unstranded, hybrid-capture regime is the one Rigel exists to handle. In it, the strand channel is silent
and the structural channel supplies only gDNA. **Spliced fragments are the sole direct evidence of RNA.** A
calibration pass that cannot use them will propagate gDNA unchecked and mis-classify expressed loci as genomic
background — the failure mode we observe. So "the spliced message is refuted" cannot be the end state: some
correct form of this message must carry the RNA evidence. We are asking reviewers to help us find the correct
form of `I_motif`, on the four precise points of §5, rather than to accept or reject a single formulation.

---

## Appendix — notation

```
  λ = logit(f_g)            gDNA-vs-RNA log-odds of a node's UNSPLICED mass
  f_g / f_r = 1 − f_g       gDNA / RNA fraction of the unspliced mass
  κ                         library sense fraction (½ = unstranded ⇒ strand Fisher info = 0)
  S_b, S_eff,b              spliced count / its Kish effective count at boundary b (pure RNA)
  E_spl,b                   boundary b's spliced effective length (FL-geometric; always computable)
  ρ_spl,b = S_b / E_spl,b   measured spliced-RNA density (the message MODE)
  pr_motif                  the spliced message PRECISION (this note's subject)
  σ²_transfer(b→e)          belief-free enrichment-crossing transfer variance (existing machinery)
  N_src                     an emitting node's own fragment count (Poisson sampling term)
```

*Companion (full architecture, for a reviewer who wants the surrounding design):*
`docs/calibration/message_precision_derivation.md`.
