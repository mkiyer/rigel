# ROADMAP — the short ranked view

**What this file is.** The one-line-per-claim state of the tool and the ordered next steps — nothing
else. Three rules keep it short: the substance of every item lives in `ISSUES.md` (the open entries plus
the append-only CLOSED / REFUSED record); the changelog is git, so this file records no history; and no
figure lives here — a claim names the instrument that re-derives it (owner, 2026-08-22). How performance
is judged is `SUCCESS.md`; rulings are `DESIGN.md`; lessons are `TRAPS.md`, cited by name.

## The frame (owner, 2026-09-22 and 2026-09-23)

The version on disk is `pyproject.toml`'s; the target is 0.8.0, A RELEASE OF THE TOOL (`DESIGN.md` §0b).
Two numbers are primary and answer different questions: the transcript table against per-transcript truth is
what the release ships on (`quant_accuracy.py`, per stratum, above `--arm base_reseed`, under fractional
assignment), and the calibration result against an oracle calibration is what ranks a calibration mechanism
(`calibration_vs_oracle.py`, `solvability_audit.py`, `prior_vs_oracle.py`). Three strata are in scope;
unstranded × capture-ON is deferred and never ranked on a pooled total (`TRAPS: never-pool-the-strata`). The
fragment-length composition channel is retired until after 0.8.0.

The cleanup is closed: the tree is the shipped infrastructure with nothing dead, duplicated or narrated, and what
it deliberately left is `ISSUES: hygiene-ledger`. The capture-contracted length of the locus gDNA component,
the synthetic spans and the annotated transcripts is kept as the one shared rule (owner, 2026-09-23,
`DESIGN.md` §7.2); robustness over synthetic accuracy is `DESIGN.md` §0b. The work, in order, each step derived
on one page and A/B'd against what ships: **the junction price's precision within a gene**, then **the
pseudocount**, now uncancelled on capture.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `calibration_oracle.py`.
- **Library gDNA fraction**: calibration's is accurate on the three in-scope strata and structurally blind on
  the deferred one — `solvability_audit.py`, `policy_benchmark.py --by-class`; the transcript table biases it
  toward gDNA off capture through the pseudocount, and on capture that bias no longer has the gDNA component's
  length to cancel it — `quant_accuracy.py` (the pools per row),
  `ISSUES: the-pseudocount-prior-is-biased-toward-gdna`.
- **The deliverable, end to end**: measured per stratum on the rebuilt ladder under the shared length; what
  remains in scope is the isoform split under capture, where the junction price is noisy within a gene, and the
  gDNA split the pseudocount biases — `quant_accuracy.py --arm base` above `--arm base_reseed`,
  `--arm oracle_ruler`, `ISSUES: the-junction-price-is-noisy-within-a-gene`.
- **Fragment lengths**: closed, both halves — `calibration/fl.py`, `gdna_density.py`; watch
  `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: robust to the annotation (`EQUATIONS.md` §6a–§6c); on real data read
  `clamped_at_ceiling` and `effective_seeds`, never the bare value.
- **The message layer**: `transfer` ships on the two-phase backbone (`DESIGN.md` §6b.12–§6b.14); `silent` is the
  floor; the bar — win on unstranded, minimal harm on stranded, never pooled — is `policy_benchmark.py --panel
  ladder`.
- **The gDNA landscape prior**: done for 0.8.0 (`DESIGN.md` §7.1) — `calibration_vs_oracle.py`.
- **ψ**: the composition closes structurally on every published object (`test_vertex_reference.py`); the tilt is
  integrated on derived nodes (`EQUATIONS.md` §9e–§9f); the λ bracket follows the landscape prior's demand.
- **The prior assembler**: with perfect masses its own error is negligible — `prior_vs_oracle.py`; a perfect
  `LocusPriors` is worth little in scope — `quant_accuracy.py --arm oracle`; the per-transcript lane the EM never
  receives is worth far more — `quant_accuracy.py --arm oracle_alloc_seed`, `ISSUES: per-transcript-prior-lane`.
- **The capture-contracted length**: one shared rule for every EM component — each object's conserved share at
  that object's own capture efficiency, read against the landscape's located enriched mode, a junction priced
  from its neighbours by conservation of bases (`DESIGN.md` §7.2, `EQUATIONS.md` §11); the classes sit near one
  scale and the junction price is noisy within a gene — `ruler_vs_truth.py --scale` (the class means and the
  within-gene spread), `ISSUES: the-junction-price-is-noisy-within-a-gene`; what the gDNA witness cannot see of a
  transcript-designed panel is declared — `ISSUES: ruler-witness-geometry-on-transcript-panels`.
- **Performance**: the sweep is one native call, bit-identical at every thread count; the thread is PARKED —
  `ISSUES: performance-memory-bounded-solve`, `profiling/profiler.py`, `profiling/sweep_replay.py`.
- **Panels**: the 16-condition ladder is rebuilt under the corrected capture physics, cached and certified; the
  test chromosome is cached and certified; the junction-probed twin and the fl-gap side panels are stale —
  `panel.py status`, `ISSUES: flgap-panels-stale-nascent-model`.
- **Attribution floor**: the shipped assignment is a sampled draw; every `quant_accuracy` arm runs fractional
  and no delta below `--arm base_reseed` is attributable (`TRAPS: the-deliverable-is-not-reproducible-by-default`).
- **Reading rules**: rank per stratum; quote `mwae_all` / Σ|err| and the shipped column, never pass-0
  (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the order

1. **The junction price's precision within a gene** — the capture efficiency of a splice junction inside a
   multi-exonic transcript, imputed only from the regions and boundaries within one fragment of it and never
   pooled across junctions, because some junctions are probed and some are not (owner, 2026-09-23;
   `ISSUES: pooling-junctions`); judged with no EM by the within-gene spread beside the class scale
   (`ruler_vs_truth.py --scale`, `TRAPS: judge-a-ruler-by-its-within-gene-spread`), then on the transcript table
   per stratum (`ISSUES: the-junction-price-is-noisy-within-a-gene`).
2. **The pseudocount** — now uncancelled on capture; its repair is taken with the per-transcript prior, through
   the lane that exists (`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`,
   `ISSUES: per-transcript-prior-lane`), judged on the pools per row first and the transcript table second, `g98`
   and the zero controls beside.
3. **The release** — `docs/PUBLISHING.md` is the procedure; what gates it is the state: the deliverable
   measured per stratum, the zero rows clean, the suite at its standing count, `preflight.py --full` green, the
   standing risks re-read, and the manual true of what ships.

**Parked, each with its entry**: the AMBIG slots' remaining defects (`ISSUES: gdna-landscape-trains-on-false-positives`,
`ISSUES: the-tilt-census-as-an-instrument`, `ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`) · the intron's own
solve on unstranded capture-OFF · the message policy's open cases (`ISSUES: two-sided-exon-row`,
`ISSUES: flux-floor-dispersion`, `ISSUES: message-layer-open-cases`, `ISSUES: refit-vs-message-arbitration`) ·
`expand-the-gdna-spectrum` · `transfer-variance-premise` · `nascent-stress-sensitivity` · `hygiene-ledger` ·
`splice-out-premise-bias-uncorrected` · `drain-contaminates-certified-rna` · `crossing-pool-contrast` ·
`pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk` · `performance-memory-bounded-solve`.

## Deliberately not next

The length composition channel (retired until after 0.8.0) · a capture efficiency the EM re-reads as it runs
(deferred by the owner, `DESIGN.md` §7.2) · anything whose only target is the deferred stratum · every mechanism
in `ISSUES.md`'s CLOSED / REFUSED section — read it before proposing anything, because each entry is a build that
was measured and turned down, with the number that killed it.
