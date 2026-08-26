# THE ANCHOR IS A MESSAGE — the integration design (owner ruling, 2026-08-25)

**The ruling.** The spliced-fragment anchor is message propagation: a one-hop imputation from the
flank boundary into the exon. It may not exist as an appendage beside the message framework — it
must live inside it. The sender publishes its spliced-fragment observation unchanged; the
RECIPIENT decides how to use it, and the recipient's arithmetic is the route-sum + NB marginal.
Whether it lands in an existing policy or a new one is left to this design.

**The decision: no new policy.** The anchor is a new *stream* — the CERTIFIED-FLUX stream — not a
new composition scheme. Policies differ by how they compose neighbour claims; the anchor is
evidence content. Forking a fourth policy to carry it would repeat the pattern that produced the
relay's operator pile. It lands in `RelayPolicy` behind one switch, `SilentPolicy` stays silent
(the control is restored — pre-2026-08-25 the anchor rode the local solve of every policy,
control included), and `CurrencyPolicy` stays frozen untouched.

**What the survey found, and why this is a homecoming rather than a port.** The framework already
carried every piece the appendage rebuilt:

* `StepContext.inv_sj_lo/hi` is the per-face, per-strand, model-free certified-RNA measurement —
  and it is a per-fragment reciprocal-opportunity SUM, i.e. route-sum-correct by construction.
  The appendage rebuilt this observation from a RouteTable.
* The relay's splice-in stream already delivers a certified-RNA claim to exons — built from
  `sj_count/eff_sj` as a per-face density, which is the RATIO-OF-SUMS pooling the round-2 review
  refuted (a k-route face reads ~k× low). ⭐ So the route-sum fix is a repair to the framework's
  observation surface, not an anchor-private trick, and the anchor stream is the UPGRADE of the
  splice-in claim: same certified observation, correctly pooled, delivered as a count likelihood
  instead of a Gaussian.
* `rna_one_sided` is plumbed through ψ but never populated by the relay — the one-sided
  certified-RNA bound is this stream's ancestor and becomes an obsolescence candidate once the
  two-sided NB claim ships through the same channel.

## The contract changes, smallest set

1. **`PsiMessage.lam_rows`** — one new field: an `(n_slots, K)` λ-factor row array, or `None`.
   This is the message layer adopting ψ's existing general evidence currency (`lam_logprior`
   rows: θ-independent, finite, zero row = inert, row-max 0). The backbone asserts finiteness and
   shape beside the existing message assertions. `is_silent` lists it.
2. **The backbone sums it into the FINAL solve only**: `sweep.solve_chain`'s final ψ call passes
   `lam_logprior = intron_prior ⊕ msg.lam_rows`. Phase-A (`build_region_init`) and the
   own-evidence precision (`density_factor_precision`) never see it — that is the citizenship
   move. The intron factory keeps own-evidence citizenship (it is the slot's own density against
   a library constant, no hop); the anchor loses it (it is a neighbour's observation, a hop, an
   imputation).
3. **DEFERRED (recorded, not built): `RegionGeometry.route_rate_lo/hi`** — folding the
   route-summed certified rate per boundary face into the geometry (and fixing the relay
   splice-in's ratio-of-sums pooling at the source) is the follow-up, priced with the
   weak-message re-price. Round 1 keeps the `RouteTable` inside the prepared evidence bundle so
   the priced arithmetic moves byte-identically.
4. **`RelayPolicy(switches, flux=...)`** — constructor injection of the prepared evidence
   bundle (`rna_anchor.prepare_flux_evidence`: the selections, the route-summed rates, the
   nascent frame arrays, the intergenic background). The policy's grid-keyed memo builds the
   rows (`rna_anchor.flux_rows` — the recipient arithmetic: estimator fits, NB quadrature at
   claimed exons, guarded Gaussian at eligible boundaries) and `deliver` ships them as
   `lam_rows` behind the `certified_flux` switch. The arithmetic stays in `rna_anchor`
   (layer 5); `messages/` (layer 6) imports DOWN into it. `build_rna_anchor_factor` survives as
   the composed arithmetic REFERENCE the parity gate holds the policy to.
5. **`calibrate.py` loses the appendage**: `_anchor_at` / the factor summing are deleted; the
   policy constructor receives the bundle. `config.rna_anchor` keeps its name and gates the
   stream: it is live iff `rna_anchor ∧ message_propagation ∧ message_policy == "relay"`.

## What changes behaviourally, and the acceptance bar

* The anchor leaves phase-A and the own-evidence precision: anchored slots revert to
  no-own-composition-evidence status, which changes the relay's mismatch-deflation inputs and the
  refit training populations. It still shapes refits through the delivered beliefs (the fits
  train on post-sweep beliefs), which is the honest route.
* SILENT REVERTS TO THE PRE-ANCHOR CONTROL. Its ladder rows will regress to pre-anchor levels by
  construction; that is the control being restored, not a loss.
* The shipped candidate is relay-with-stream. Acceptance: per in-scope stratum, within noise of
  or better than the graduated relay+appendage readings (`CLEAN_LIBRARY_RESIDUAL.md` §5d /
  the 16-row all-policies table, 2026-08-24), both zero controls held, pure-gDNA/RNA-bearing split
  never pooled. The arithmetic-parity gate (relay rows ≡ appendage rows on the anchor toy)
  separates a citizenship effect from an arithmetic bug by construction.

## What this kills or marks

* DEAD NOW: `calibrate._anchor_at`, the factor-sum contract comment, the appendage's RouteTable
  memo (the route table's consumers move into geometry/the bundle builder).
* OBSOLESCENCE CANDIDATES (measured later, in the weak-message re-price): the relay's Gaussian
  splice-in certified-RNA claim (superseded by the stream's NB claim — same observation, better
  law), `ONE_SIDED_RNA` and the `rna_one_sided` plumbing, `certified_q` machinery.
* THE RECORDED DEBTS CARRY OVER UNCHANGED: single-flank centre bias, sj strand-column matching,
  shared-fragment pairs, the capture-ON +12 % mature offset, the transport-dispersion
  decomposition (the "why do two flanks disagree at infinite depth" derivation) — now all
  first-class message-layer questions, which is where they always belonged.
