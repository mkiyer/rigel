# ROADMAP — the short ranked view

**What this file is.** The one-line-per-claim state of the tool and the ordered next steps — nothing
else. Three rules keep it short: the substance of every item lives in `ISSUES.md` (the open entries plus
the append-only CLOSED / REFUSED record); the changelog is git, so this file records no history; and no
figure lives here — a claim names the instrument that re-derives it (owner, 2026-08-22). How performance
is judged is `SUCCESS.md`; rulings are `DESIGN.md`; lessons are `TRAPS.md`, cited by name.

## The 0.8.0 frame

The version on disk is `pyproject.toml`'s; the target is 0.8.0, A RELEASE OF THE TOOL (owner, 2026-09-19,
`DESIGN.md` §0b's amendment). Two numbers are primary and they answer different questions: the transcript
table against per-transcript truth is what the release ships on (`quant_accuracy.py`, read only above its
reseed floor), and the calibration result against an oracle calibration is what ranks a calibration
mechanism (`calibration_vs_oracle.py`, `solvability_audit.py`, `prior_vs_oracle.py`). Neither stands in for
the other: the transcript table cannot say whether calibration or the EM moved, and a calibration figure
cannot say whether the user's number improved. Three strata are in scope (unstranded × capture-OFF, stranded × capture-OFF, stranded ×
capture-ON); unstranded × capture-ON is deferred — reported on every benchmark, never a development
target, and never ranked on a pooled total (`TRAPS: never-pool-the-strata`). The fragment-length
composition channel is retired until after 0.8.0. The full ruling, including why the ladder gives gDNA
and RNA equal fragment lengths, is `DESIGN.md` §0b.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Library gDNA fraction**: calibration's is accurate on the three in-scope strata and structurally blind on the
  deferred one (at κ = ½ no channel reaches an AMBIG slot; the θ-independent-channel search is closed) —
  `solvability_audit.py`, `policy_benchmark.py --by-class`. The transcript table the user reads does not keep it,
  in either direction: at capture-OFF the EM assigns to gDNA unspliced RNA the simulator drew as nascent
  (measured only at the panel's nascent stress share), and under capture it gives probed gDNA to isoform-rich
  genes, both against a calibration that had the split right (`ISSUES: em-overturns-the-calibrated-gdna-split`).

- **The deliverable, end to end** (2026-09-19; `quant_accuracy.py --set em.assignment_mode=fractional`,
  transcript-level Σ|Δ| as a share of the true RNA at `g00` / `g05` / `g50` / `g98`): unstranded OFF 3.1 / 3.4 /
  4.2 / 35.0 %, stranded OFF 3.3 / 2.9 / 3.9 / 24.8 %, stranded ON 7.4 / 10.5 / 13.9 / 98.8 %, deferred 9.0 / 16.7 /
  28.3 / 804 %; gene level in scope 1.3–5.2 % below `g98`, where RNA is 2 % of the library and the error is gDNA
  on transcripts. A perfect prior recovers at most 5 % in scope and nothing measurable at `g05`: the residual is
  the EM's and the assignment's (`ISSUES: end-to-end-error-unattributed`, CLOSED). Stranded × capture-ON, the
  stratum above the owner's 5 % bar, is the capture ruler under junction-spanning probes — the simulator's own
  capture-aware length takes it to 2.9 / 3.6 / 8.4 % (`ISSUES: ruler-witness-geometry-on-transcript-panels`) —
  and what remains at `g50` is the EM's gDNA under-call (`ISSUES: em-overturns-the-calibrated-gdna-split`). Truth
  as the per-transcript weights removes about half of the error everywhere (`ISSUES: per-transcript-prior-lane`).
- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `calibration_oracle.py`.
- **Fragment lengths**: closed, both halves — gDNA by the two-pool contrast (`calibration/fl.py`,
  `gdna_density.py`; gates `test_fl.py`, `test_gdna_density.py`), RNA sound as shipped
  (`ISSUES: the-rna-length-law-fix`, CLOSED). Watch: `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: robust to the annotation (`EQUATIONS.md` §6a–§6c, `DESIGN.md` §3.3a);
  on real data read `clamped_at_ceiling` and `effective_seeds`, never the bare value.
- **The message layer**: `transfer` ships on the two-phase backbone with the level lanes
  (`DESIGN.md` §6b.12–§6b.14); `silent` is the measured floor. The bar — win on unstranded, minimal harm
  on stranded, never pooled — is `policy_benchmark.py --panel ladder`; the zero rows are solved by the
  prior under both policies, so the "beats silence" count is read on the contaminated rows, where every
  row favours `transfer`; `calibration_walk.py` says the messages still carry the stranded capture-ON
  rows and are essential on the deferred stratum.
- **The gDNA landscape prior**: done for 0.8.0 (`DESIGN.md` §7.1); the zero controls are solved on the
  metric (`calibration_vs_oracle.py`) and the
  in-scope per-object composition error did not move.
- **ψ**: the composition closes structurally on every published object (`test_vertex_reference.py`);
  the reference location is deleted (`DESIGN.md` §6b.1); the tilt is integrated on nodes that follow each
  slot's strand term, a derived count and no lattice (`DESIGN.md` §6b.15.11, `EQUATIONS.md` §9e), exact at any
  depth, with the tilt's hypothesis space {pure +, pure −, mixed} (`EQUATIONS.md` §9f); the λ bracket follows
  the landscape prior's derived demand (`landscape.required_logodds_window`).
- **The prior assembler**: with perfect masses its own error is negligible — `prior_vs_oracle.py` — and a
  perfect `LocusPriors` is worth 3–5 % of the in-scope transcript error — `quant_accuracy.py --arm oracle`. The
  lane the EM never receives is worth far more: `rna_prior_weight` is built end to end and `pipeline.py` omits
  it, so the shipped EM carries no per-transcript information at all (`ISSUES: per-transcript-prior-lane`).
- **The ruler is the transcript's bases at their pieces' capture efficiencies, against the landscape's
  located enriched mode** (`DESIGN.md` §7.2, `EQUATIONS.md` §11): each efficiency a posterior mean from the
  piece's own count and its edge crossings, no floor and no junction object, so the unprobed class reads
  within ±0.2 nat of the simulator's truth where the floor read +3.4 (`ruler_vs_truth.py`); with no enriched
  gDNA mode nothing contracts, so the zero controls and both capture-OFF strata read a factor of exactly
  1.000 with nothing moved; a mode's members are kernels with a location, so a sparse library is told it
  has no reference rather than handed one read off anchors' walls (`gdna_reference_members` is the
  regime); what the gDNA witness cannot see of a transcript-designed panel is declared
  (`ISSUES: ruler-witness-geometry-on-transcript-panels`) — and on the ladder's junction-spanning panel it costs
  more than it corrects: switched off, stranded × capture-ON at `g05` falls from 10.5 to 7.8 %
  (`quant_accuracy.py --arm oracle_ruler` prices it end to end); the never-passed per-transcript prior lane
  (`ISSUES: per-transcript-prior-lane`) is the other pre-EM item.
- **Performance**: the port and the work outside it are done, and a deep run is 0.87 of what it was
  (`DESIGN.md` §6b.15). The sweep is ONE native call over a pool of threads, bit-identical at every thread
  count; the scan's split, the second pass's lookups, the fragment-length fits and the sweep's arena followed.
  What remains is ranked and resumable in `ISSUES: performance-memory-bounded-solve`, whose next item — the
  short-template taper table — is already derived and measured two ways, one bit-identical and one priced. The
  thread is PARKED by the owner, not finished — `profiling/profiler.py`, `profiling/sweep_replay.py --threads`.
- **Panels**: the sparse-nascent 16-condition ladder and the 30-condition test chromosome, both cached
  and certified — `panel.py status`; the fl-gap side panels carry a different nascent model —
  `ISSUES: flgap-panels-stale-nascent-model`. The ladder's nascent level is a development stress
  (`DESIGN.md` §0b).
- **Oracle FIELD certification**: every ladder row is stamped, but the uniformity gate is vacuous on
  capture-ON and zero-gDNA rows — read the stamp with its vacuity flag — `calibration_oracle.py`.
- **Attribution floor**: the shipped assignment is a sampled draw, so every `quant_accuracy` arm runs with
  `--set em.assignment_mode=fractional` (owner, 2026-09-19) and the report refuses to mix modes; no delta below
  the reseed floor is attributable — re-derive `--arm base_reseed` in the same session
  (`TRAPS: the-deliverable-is-not-reproducible-by-default`).
- **Reading rules**: rank per stratum; quote `mwae_all` / Σ|err| and the shipped column, never `solv%`
  or pass-0 (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the recommended order

The owner's standing order (2026-09-19): THE METHOD, back to the front. Performance ran until the tool was
fast enough to iterate on — a deep run is 0.87 of what it was and every step of it was bit-identical
(`DESIGN.md` §6b.15) — and it is now a PARKED THREAD of its own kind: machine work, judged on seconds and
bytes, with its own ranked list and its own instruments in `ISSUES: performance-memory-bounded-solve`.
Nothing it did moved a number, so every accuracy measurement recorded before it still stands.

The ranked list below is the method: what the tool ANSWERS, judged on 0.8.0's metric per stratum. It is the
substance of each item.

The method is the dissection loop: run the panel → worst in-scope scenario → rank its objects by error
mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated fix → add the offending
transcripts to the test chromosome → re-run → repeat. The facts this ranking leans on, each named with
its instrument: the ruler reads 1.000 at `g00` and off capture, so the metric page is the composition's
(`calibration_vs_oracle.py`); a perfect prior is worth 3–5 % of the in-scope transcript error end to end and
nothing measurable at `g05` (`quant_accuracy.py`); by class the in-scope residual sits on the
intron's own solve (unstranded OFF) and on exon|exon boundaries and walled exons (stranded ON)
(`policy_benchmark.py --by-class`).

1. **Stranded × capture-ON: the capture ruler under junction-spanning probes** (owner, 2026-09-19: the focus,
   with a target below 5 %) — `ISSUES: ruler-witness-geometry-on-transcript-panels`. Root-caused: a probe across a
   junction makes capture isoform-specific, the EM's isoform split reads the ratio of capture-aware lengths, and
   the gDNA witness cannot see a junction (at zero gDNA there is no witness at all); the simulator's own length
   takes the stratum to 2.9 / 3.6 % at `g00` / `g05`. The repair needs an owner decision, because the one
   observable that sees isoform-specific capture is the probe design and Rigel reads no panel — and the first
   question is whether the panels Rigel will meet span junctions. Judged on `quant_accuracy.py` with
   `oracle_ruler` as the ceiling and the test chromosome's benign and junction-probed panels as the control.

2. **The EM's residual: the per-transcript allocation** — `ISSUES: per-transcript-prior-lane`. `rna_prior_weight`
   is built end to end and never passed; truth as the allocation weights removes about half the transcript
   error in every stratum, the gDNA-free `g00` rows included, and still halves what `oracle_ruler` leaves on
   stranded × capture-ON. A wiring gap plus a support decision: that arm hands over the true support, two
   weightings are refused with their numbers, and the next candidate is a sparsity mechanism.
   `ISSUES: nested-antisense-leak-under-the-sane-ruler` sits in it.

3. **The EM does not hold calibration's gDNA split** — `ISSUES: em-overturns-the-calibrated-gdna-split`, in both
   directions: at capture-OFF it takes unspliced RNA into gDNA (`g05 ss.50 OFF` 0.104 against 0.05, measured only
   at the nascent stress share, so SIZED first at the realistic share, `ISSUES: nascent-stress-sensitivity`), and
   under capture it gives probed gDNA to isoform-rich genes (`g50 ss.99 ON` −3 % against a calibration at +0.4 %),
   under a perfect prior, a true ruler and a true gDNA length alike.

4. **The pre-EM prior chain** — what the oracle arms price it at on the deliverable. `rna_prior_count` carries
   almost all of a perfect prior's value, at `g98`; `gdna_eff_len` at truth moves nothing above the floor in
   scope, so `ISSUES: capture-blind-gdna-divisor` and `ISSUES: eb-shrinkage-magic-ess` are near zero on the
   ladder's table (the ESS still dominates the fl-gap arm); `ISSUES: antisense-prior-assembly-casualty` owns an
   xfail. Each judged on the deliverable AND on `prior_vs_oracle.py`, so a repair that moves the prior and not
   the user's number is visible as that.

5. **Calibration accuracy where the strand tilt matters** — the AMBIG slots with RNA on both strands
   (`DESIGN.md` §6b.15.12–§6b.15.13). The tilt atom and the strand channel's protocol decision landed 2026-09-14 (the
   strand-pure under-call and the gDNA-free deadband CLOSED); the θ measure is settled (both flattenings
   REFUSED, `ISSUES: strand-marginal-volume-factor`); the lanes' own defects are fixed and gated by
   `test_encompassing_locus.py`. What is left is the landscape estimator's vertex bias
   (`ISSUES: gdna-landscape-trains-on-false-positives` (d) — the lever on the stranded zero controls and
   on every unwitnessed both-strand slot), the census as an instrument
   (`ISSUES: the-tilt-census-as-an-instrument`), and a known limit to watch rather than build against
   (`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`). Each judged on the metric per stratum, both
   zero controls and the shared-exon stress at depth, never on the ladder alone.
6. **The intron's own solve on unstranded capture-OFF** — the intron class carries the largest share of
   the in-scope error there (`policy_benchmark.py --by-class`): the factory profile's resolution against
   the intergenic background (`density_deconv`); dissect with `worst_objects.py`.
7. **The vertex atom** — on silent genes and nascent-free introns; a
   mechanism for it is the prior's reference (`ISSUES: reference-prior-refuted-at-concept-level`
   constrains the form) or the intron's own solve, not a message.
8. **The message policy, only where a row is above the bar**: one prototype mechanism at a time, in C++ in
   a worktree, the two trees scored with `policy_benchmark.py --by-class`, halves apart, pass zero beside the pipeline:
   `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`,
   `ISSUES: message-layer-open-cases`.

Then, in standing order: `ISSUES: refit-vs-message-arbitration` (re-read under the E-step: the walk now says the prior does the
unstranded rows and the messages the stranded capture-ON ones).

9. **The release itself** — `docs/PUBLISHING.md` is the procedure and it is two commands plus a wait. What
   gates it is not the procedure but the state: the deliverable measured and not regressed per stratum, the
   zero controls at 0.000 and 1.000, the suite at its standing count, `preflight.py --full` green, the
   standing risks re-read (`ISSUES: capture-degeneracy-standing-risk`,
   `ISSUES: flgap-panels-stale-nascent-model`), and the manual true of what ships.

**The other kind of work, parked and resumable**: `ISSUES: performance-memory-bounded-solve` carries the
machine thread — what a deep run costs now, what is ranked next with its measured price, and the two
candidates already researched and not taken. It resumes without re-deriving anything, and it is judged on
seconds, bytes and bit-identity rather than on the metric.

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `transfer-variance-premise` ·
`nascent-stress-sensitivity` · `hygiene-ledger` · `flgap-panels-stale-nascent-model` ·
`background-abundance-pair-unruled` · `splice-out-premise-bias-uncorrected` ·
`drain-contaminates-certified-rna` (the ceiling refused the in-solve correction; two recorded
follow-ups) · `crossing-pool-contrast` (blocked) · `pure-rna-mirror-asymmetry` ·
`capture-degeneracy-standing-risk`.

## Deliberately not next

The length composition channel (retired until after 0.8.0) · anything whose only target is the deferred
stratum · every mechanism in `ISSUES.md`'s CLOSED / REFUSED section — read it before proposing anything,
because each entry is a build that was measured and turned down, with the number that killed it.
