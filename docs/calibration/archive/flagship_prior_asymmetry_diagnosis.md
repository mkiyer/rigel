# Flagship dissection — the error is the gDNA prior's CLAMPED LEFT TAIL (an asymmetric ψ), 2026-07-17

**Status:** debugging result (meticulous localization before any premise change). Tools:
`scripts/debug/flagship_dissect.py` (three-stage decomposition via the `_capture` hook),
`scripts/debug/prior_crush_probe.py` (per-node ψ balance). Condition: `ss_0.50_nrna_present_capture_on` (the
unstranded+capture DIAGNOSTIC — not a practical library, but it exposes what strand-specific data hides by
self-solving).

## What we set out to answer
The unstranded-capture error (mwae ≈ 0.57) cannot be uniform — which node subset breaks, and why? Is it a
strand-deposit bug, a message bug, or a modelling error?

## The three-stage decomposition (mass-weighted |f_g − oracle|)
The `_capture` hook exposes the belief at three stages: strand-ALONE (`fg_strand`), strand+prior (`fg_loc`,
no messages), and final (`f_g`, +messages).

| class | mass frac | strand-only | +prior | +msgs (final) |
|---|---|---|---|---|
| ALL | 1.000 | **0.329** | **0.622** | 0.570 |
| AMBIG | 0.189 | 0.320 | 0.701 | 0.671 |
| single | 0.691 | 0.284 | 0.646 | 0.587 |
| gonly (locked gDNA) | 0.041 | 1.000 | 0.000 | 0.000 |

**The prior nearly DOUBLES the error** (0.329 → 0.622); the messages then partially RECOVER it (0.622 → 0.570).
So strand-only (the symmetric Beta(½,½) Jeffreys, no NPMLE) is far better than strand+NPMLE. The error localises
in the strand-blind classes (AMBIG, single); the locked-gDNA class (`gonly`) is served correctly by the prior.

## Ruled out (not a bug there)
* **Strand deposit** — oracle gDNA pos-frac = 0.500, RNA pos-frac = 0.500. Unstranded ⇒ 50/50 for both
  components, exactly as it must be. No wrong-strand deposit.
* **Messages** — on the 3132 nodes that receive a gDNA message, mass·help = 469 308 vs mass·hurt = 13 737
  (**34:1 in favour of HELP**). The messages are the RESCUE, not the disease; σ²_transfer gagging them is a
  secondary loss on top of the prior's primary damage.
* **Direction** — under-calling gDNA (boundaries signed −0.437; truly-98%-gDNA single-strand nodes solve to
  f_g ≈ 0.00). Enriched gDNA is being called RNA.

## Root cause — the asymmetric ψ (a CLAMPED left tail)
Per-node ψ balance for node 1909 (truly 98% gDNA, M=66 547, log10 ρ_tot=1.53, sits on the enriched mode 1.48):

```
gDNA NPMLE prior term : f_g 0.01:-7.42  0.1:-5.64  0.5:-5.49  0.9:-5.04  0.99:-5.06   (favours high f_g by +0.60)
RNA ref ½log(1−f_g)   : f_g 0.1:-0.05   0.5:-0.39  0.9:-1.15  0.98:-1.96             (penalises high f_g by -1.10)
NET ψ (prior + RNAref) argmax f_g = 0.002        ← crushed to zero
+restore ½log(f_g)    argmax f_g = 0.700          ← moves to the truth
```

The strand-only solve uses the **symmetric** Beta(½,½) Jeffreys `½log f_g + ½log(1−f_g)`: the `½log f_g` arm →
−∞ at f_g→0 and PROTECTS gDNA nodes from the RNA arm's `½log(1−f_g)` barrier at f_g→1. The NPMLE was designed to
REPLACE the `½log f_g` arm — but it **clamps its left tail flat** (`gdna_rate_prior.logprior` uses
`np.interp(..., left=logP[0])`, a constant). A proper prior `P(log ρ_g)` must decay to −∞ as `ρ_g→0` (f_g→0);
the clamp makes it a finite floor (−7.42 at f_g=0.01) with essentially no gradient. So:

* the gDNA side lost its f_g→0 barrier (clamped flat, +0.60 nats of weak preference);
* the RNA side kept its full f_g→1 barrier (−1.10 nats);
* net → strand-blind nodes slide to f_g ≈ 0 (called RNA). Restoring the `½log f_g` barrier moves node 1909's
  argmax 0.002 → 0.700, node 2999 (truly 21% gDNA) 0.000 → 0.300 — both toward truth.

## This is a KNOWN deferred item, now proven load-bearing
`CLEANUP_LOG.md` already lists: *"the logP_g left tail must land with the reference replacement (Phase 3) —
replace the clamped `np.interp(..., left=logP[0])` with a genuinely decaying left tail
(reference_prior_derivation.md §10.8)."* And `calibration_roadmap.md` flags `logP_r` as UNWRITTEN (the RNA arm is
still the Jeffreys reference, un-refit). The dissection shows this pair — clamped gDNA left tail + un-refit RNA
arm — is not a tail nicety; it is the DOMINANT flagship error. Fixing the reference asymmetry, not band-aiding
with strand-gating, is the path.

## Candidate fix (to validate next, carefully — it touches the ψ reference)
Give the gDNA prior a genuinely decaying left tail so `logP_g → −∞` as `ρ_g → 0` (equivalently, retain a
`½log f_g` vertex barrier alongside the NPMLE population shape). Symmetrically, the principled endpoint is to
WRITE `logP_r` (a fitted RNA prior) so both arms are proper and balanced, retiring the lopsided Jeffreys RNA
reference. Validate on the flagship at the FULL-solve level (not just the local ψ argmax) and confirm it does NOT
manufacture gDNA in the zero-gDNA / stranded controls before adopting.

## Full 24-scenario A/B — V4 is REJECTED (2026-07-17, `scripts/debug/test_v4_symmetric.py`)
Ran the candidate fix (V4: NPMLE + restored symmetric Jeffreys barrier) across ALL 24 ambig scenarios
(gdna_none / gdna100 / gdna300 × ss × nrna × capture) — the `gdna_none` set is the zero-gDNA false-positive
control. **V4 is rejected:** it fixes the gDNA-present flagship (gdna300 unstranded capON 0.63→0.20, 0.57→0.17;
gdna100 0.45→0.15) but CATASTROPHICALLY manufactures false gDNA in the **zero-gDNA unstranded** rows
(gdna_none ss0.50 none: 0.004→0.40 / 0.01→0.43 — a 40% false-positive gDNA). Net over 24: 0.105→0.097 (flat).
No-prior (symmetric Jeffreys + strong messages, no NPMLE) ALSO fails gdna_none unstranded (0.49). Stranded is
fine under every variant.

**The fundamental structure (the real finding).** In the UNSTRANDED conditions the truth is OPPOSITE by gDNA
level: gDNA-present wants f_g→1, gDNA-absent wants f_g→0. A fixed prior bias is right for one and wrong for the
other — V0 is RNA-biased (good gdna_none, crushes gdna300); V4/no-prior are neutral-or-up (good gdna300, false
gDNA in gdna_none). **The prior cannot decide f_g in the unstranded case; the gDNA LEVEL is a strand-independent
DATA quantity** — the pure-gDNA intergenic/intron density (dense in gdna300, empty in gdna_none). The current
NPMLE is fit on TOTAL density (gDNA+RNA conflated), so it is blind to this discriminator. **Correct direction
(the owner's long-standing directive `empirical_gdna_prior_from_intron_intergenic_regions`): anchor the gDNA
abundance on the pure-gDNA intergenic+intron regions** so the prior defaults to RNA when intergenic is empty and
permits gDNA (messages lift f_g) when intergenic is dense. The asymmetric-ψ / clamped-tail fix (V4) is necessary
but NOT sufficient — it must be paired with a data-anchored gDNA level, or it manufactures gDNA where there is
none.

## Reproduce
```
OMP_NUM_THREADS=1 python scripts/debug/flagship_dissect.py --top 25
OMP_NUM_THREADS=1 python scripts/debug/prior_crush_probe.py --nodes 1909,1055,2999
OMP_NUM_THREADS=1 python scripts/debug/test_v4_symmetric.py          # all 24, V4 vs current
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --no-prior # all 24, symmetric Jeffreys baseline
```
