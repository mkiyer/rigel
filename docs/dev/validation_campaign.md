# Validation campaign — handoff

    ⚠ **A DEV DOC. Provisional, not authoritative, and nothing may cite it.** When a finding here settles,
    MOVE it to its permanent home (`ROADMAP.md` for a number, `TRAPS.md` for a lesson, `DESIGN.md` for a
    ruling, `EQUATIONS.md` for a derivation) and delete it here in the same edit. See `docs/dev/README.md`.

    Opened 2026-08-07. ⭐ The settled version of everything below is `ROADMAP.md` §1.

## The goal, in one sentence

**Characterise the whole tool end to end before adding anything to it** — with message propagation OFF and
`length_likelihood` OFF, so every number is attributable to the tool as it stands.

⛔ Both switches stay down for the entire campaign. `length_likelihood`'s panels and oracle caches are built
and ready and it is still not time: turning it on mid-campaign makes every other number unattributable.

## Where things stand

The tool **runs end to end** — 51 s on a 10 M-fragment condition, all outputs written, and **counts conserve
exactly** (`mrna + nrna + gdna == n_unambig + n_em` to the fragment; the only excluded class is `chimeric`).
So there is no plumbing to fix. What is unknown is how *accurate* it is, and where the error lives.

Calibration is accurate on three of four strata (median library `f_gdna` error **0.005–0.012**) and
**blind on the fourth**: unstranded × capture-ON reports **0.033–0.058** while the truth spans **0.00 → 0.98**.
That is not a noisy estimate, it is a flat line, and it is the whole open problem.

## The four items

### 1. calibration-prior-vs-oracle — ⛔ no instrument exists

Calibration's endpoint is the EM's prior: `LocusPriors` in `calibration/priors.py`, three float64 arrays
indexed by `multi_locus_id`:

```
gdna_prior_count   the gDNA component's Dirichlet pseudocount
rna_prior_count    the RNA group's pseudocount (the EM splits it by evidence)
gdna_eff_len       the IPR-contracted effective length of the gDNA component
```

⭐ **These are FRAGMENT COUNTS — the EM adds them straight to its own counts (`G = n_gdna + a_g`)** — which
is what makes an oracle version computable at all: the origin-split truth gives every locus its true gDNA
and RNA fragment counts, so "the correct prior" is a well-defined number rather than a modelling opinion.

⚠ **Nothing compares `LocusPriors` to anything.** Grep confirms: every consumer is a unit test of its
arithmetic, never its accuracy. That instrument is the first build of the campaign.

Suggested shape — `scripts/design/prior_vs_oracle.py`, reusing `pass0_vs_oracle.py`'s oracle plumbing and
the existing oracle caches:
- build `LocusPriors` from the shipped calibration;
- build the oracle `LocusPriors` by projecting origin-split truth onto the same loci;
- report per locus and per stratum, **mass-weighted and unweighted**, and separately for the gDNA and RNA
  arms — they are different questions and pooling them would hide a sign flip.

### 2. tool-absolute-accuracy — ⚠ partial

The panels carry `truth_abundances_nrna_none.tsv` (per transcript: `mrna_abundance`, `nrna_abundance`,
`total_rna`, `spliced_length`), and `rigel.sim.analysis` exists. What does not exist is a per-condition
transcript-level accuracy table across a whole panel. ⛔ Score against TRUTH, not against a previous run.

### 3. error-downstream-of-calibration — ⛔ needs item 1 first, and it is the one that matters

Inject the ORACLE prior from item 1 and re-quantify. **This decides the project's direction.** Every
campaign for months has assumed the error is calibration's:

- if a perfect prior removes most of the error → calibration is confirmed as the bottleneck and
  `length_likelihood` becomes the next build;
- if most of the error survives → that assumption is **wrong**, and the work belongs in the EM.

⚠ It is a re-solve ceiling, so it needs a `noop` arm that is byte-identical (`TRAPS: byte-identity-gate`),
and `vertex_ceiling.py` already carries reusable override plumbing of exactly this shape.

### 4. performance — ⚠ measured on the wrong substrate

Measured on one 10 M-fragment ladder condition (35,135-node chr22 index):

```
per-locus EM        15.9 s   47 %
native scan          6.5 s   19 %
calibration          6.5 s   19 %
second-pass drain    3.5 s   10 %
                    ------
total               33.5 s
```

⭐ **Message propagation costs nothing measurable** — 33.7 s ON vs 33.5 s OFF. A hypothesis that the relay
was the expensive part is refuted, so do not start there.

⛔⛔ **This profile is upside down from the real one and both are correct.** Calibration cost is
depth-INDEPENDENT — every node in the index is solved regardless of depth — so it scales with the INDEX
while the EM scales with the DATA. On real cfRNA at genome scale (~1.5 M nodes) calibration has measured
~66 s against the EM's ~24 s, the reverse of the table above. **Optimising against this panel would tune
the wrong stage** (`TRAPS: toys-rank-hotspots-backwards`).

Both substrates are on disk:
- real cfRNA — `~/Downloads/rigel_runs/cfrna/` (4 samples, e.g. `mctp_vcap_rna20m_dna05m`, 11 G)
- the genome-scale index — `~/Downloads/rigel_runs/refs/rigel_index`

## Substrate, all cached and verified

| panel | conditions | oracle cache | what it resolves |
|---|---|---|---|
| `suite/ladder` | 36 | 341 M ✅ | gDNA 0→98 % × stranded/unstranded × capture on/off |
| `suite/flgap_short` | 4 | 32 M ✅ | gDNA ≪ RNA (realised gap −41.0 %) |
| `suite/flgap_long` | 4 | 32 M ✅ | gDNA ≫ RNA (+40.4 %) |
| `suite/pilot` | 8 | scan cache only | chr22 + ERCC, 0 %/100 % gDNA |

⭐ The flgap caches were verified by a cached re-run that was **byte-identical on all 208 scored fields**
(2m39s cold → 21 s warm), so they are genuinely hit and not silently stale.

## Open questions I could not answer

- **Is `nrna_parent_count` (44,030) meant to differ from `nrna_total` (30,125)?** Not obviously a defect,
  but nobody has stated the relation.
- **`tpm_total_rna` + the nRNA table's `tpm` sums to 1,001,068, not 1,000,000** — 0.107 %, too large for
  rounding, so the two tables normalise over slightly different denominators
  (`estimator.py:504-510`). Small, real, and unexplained.
