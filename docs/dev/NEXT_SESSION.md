# NEXT SESSION — the capture efficiency of a splice junction (owner, 2026-09-23)

Problem 2 is ported and kept: one shared rule prices every component's contracted length (owner, 2026-09-23:
"a step forward in the right direction"). The next session's prompt:

> Read `CLAUDE.md`, then `ISSUES: the-junction-price-is-noisy-within-a-gene` and `EQUATIONS.md` §11, the conserved
> frame. The task: **improve the capture efficiency calculation for splice junctions within multi-exonic
> transcripts.** Use the owner's words — capture efficiency, contracted length, region, boundary, junction,
> fragment, conserved share — and never opportunity, tilt, knife edge, family or pin in prose; if the derivation
> needs more than a page, it is off track.
>
> What ships: every EM component — the locus gDNA component, every synthetic nascent span, every annotated
> transcript — has a contracted length that is the sum, over every object its fragments deposit on, of its
> conserved share of that object times the object's capture efficiency. The shares are the deposit rule's and add
> up to the fl-marginal length exactly, per object; each region and boundary reads its own gDNA count. A junction
> is priced by conservation of bases from the objects beside it, `c_lo + c_hi − ½ (c_intron,lo + c_intron,hi)` —
> the efficiencies at the junction's low and high boundaries less the intron pieces just inside them — never below
> 0 and not clipped at 1. That price is what put the classes on one scale (18.7 % → 2.6 % apart at `g50 ss.99 ON`,
> 3.9 % at `g05`) and it is close on average (0.959 of the simulator's own junction capture over 395 junctions),
> but it is a difference of noisy boundary posteriors and each junction's error is its own, so isoforms that
> differ by a junction differ by independent noise. The within-gene sd of log(L / Y) is 0.072 / 0.057 at `g05` /
> `g50 ss.99 ON` against the old per-base ruler's 0.052 / 0.040, and the transcript table regressed there
> (stranded × capture ON 3.46 → 5.59 % and 5.50 → 6.18 %; `g98` 31.05 → 26.24 %): at `g05` the transcript error
> rose 2.13 points while the gene error rose 0.20, so it is isoform allocation. The junction at its adjacent
> pieces is as precise (0.057 / 0.041) but prices junction capture at 1/1.73 of the truth; even the true gDNA on
> every object leaves the sum at 0.064 / 0.046, so part of the noise is the rule's and part calibration's. The
> cuts beside an exon edge are not the blocker.
>
> Start also from what `ISSUES: the-junction-price-is-noisy-within-a-gene` records as unexplained or unfixed, each
> with its number there. The unprobed multi-exon transcripts at `g05 ss.99 ON` read L / Y far above the synthetic
> spans and the gDNA footprints of their own class, and `g50` does not; the stated, unmeasured hypothesis is that an
> object with almost no gDNA reads a posterior mean above its near-zero truth and the junction sum adds two such
> means. And the review's four findings: the short-intron fallback keys on a contained support of exactly 0, so an
> intron piece holding less than one expected start reads its own posterior, about the population's mean; the clip
> at 0 leaves a few transcripts carrying most of their share on junctions priced 0, contracted far below the
> adjacent-piece price; a transcript's contracted length may exceed its fl-marginal length, and its parent span's;
> and gDNA's conserved share takes an unbounded reach at every boundary, even within a fragment of a reference end,
> where the accumulator clips (consequential only on short contigs). No threshold is proposed for the first — a
> threshold is a constant.
>
> The constraints, all the owner's. Capture is LOCAL: a junction's efficiency comes only from the regions and
> boundaries within one fragment length of it (exon regions and exon|exon boundaries hold unspliced fragments
> compatible with the transcript), and a transcript-level anchor does not help predict one junction. NEVER POOL
> JUNCTIONS: "Some junctions are probed (probe spanning the junction itself) and some junctions are not probed.
> For this reason, you cannot pool junctions. It would be completely theoretically invalid and if it improves the
> results it would be purely by chance. A different type of probe panel could look much worse." ROBUSTNESS over
> synthetic accuracy: "Robustness is valued over accuracy on a synthetic simulated dataset, because real data has
> far more complexity." Conservation holds PER OBJECT, gated against the reference accumulator (`docs/TESTING.md`
> §5). ONE shared rule for every component: the junction's price may change, the rule it sits in may not. The
> length stays static: capture efficiency inside the EM is a separate branch, deferred. Refused already, each with
> its number — read them before proposing: `ISSUES: the-spliced-read-junction-price`,
> `ISSUES: a-junction-price-clipped-at-one`, `ISSUES: the-junction-price-at-its-neighbours-mean`,
> `ISSUES: pooling-junctions`.
>
> Do, in order: (1) derive on one page where the junction price's variance comes from — the boundary posteriors it
> adds, the intron pieces it subtracts, the realized counts under them — and what a local, per-junction price with
> the same mean and less variance would be. (2) Prototype it outside `src/`, A/B'd against the shipped rule on the
> same conditions: a `--module` ruler in `ruler_vs_truth.py` for the per-transcript read, a worktree for `--scale`
> and the EM, which read the shipped lengths only. (3) Judge with no EM first: `ruler_vs_truth.py --scale` on
> `g50 ss.99 ON`, `g05 ss.99 ON` and one capture-OFF row, the class means AND the within-gene spread, the
> junction-probed transcripts and the unprobed class apart (`TRAPS: judge-a-ruler-by-its-within-gene-spread`) — a
> price that restores the spread and loses the scale is the adjacent-piece price again. (4) Only then the EM:
> `quant_accuracy.py` per stratum under fractional assignment, above `base_reseed`, the pools per row before the
> transcript table, the genes beside the transcripts. The pseudocount's gDNA bias is uncancelled on capture now
> (`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`: at `g50 ss.99 ON` gDNA +172.9k, the synthetic pool 0.28×),
> so the transcript table on capture carries both errors; the gene error beside it separates the isoform allocation
> from the pools. DERIVE → PROTOTYPE → A/B; stop and discuss before anything enters `src/`; the owner drives
> commits.
