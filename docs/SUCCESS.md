# SUCCESS — how performance is judged

How Rigel's performance is measured: how every number is read, the two primary numbers and the instruments behind
them, the bars 0.8.0 is held to, and the instruments in the order you run them. Not here: the scope and its rulings
(`DESIGN.md` §0b), the panels and gates (`TESTING.md`), where the tool stands now (`ROADMAP.md`; run the instrument
for a current number), and lessons (`TRAPS.md`).

---

## How every number is read

- **Per stratum, never pooled.** Three strata are in scope — unstranded × capture-OFF, stranded × capture-OFF and
  stranded × capture-ON. Unstranded × capture-ON is deferred: it is reported on every table and is never a target,
  and because it carries most of the error, a pooled total mostly reports it (`TRAPS: never-pool-the-strata`).
- **In three pools.** gDNA, annotated RNA (the reference transcripts) and nascent RNA (Rigel's synthetic spans),
  never RNA against DNA. The annotated pool is what a user counts, so DNA called annotated RNA is the severe error and
  the annotated pool is read first; DNA called nascent RNA is reported beside it and is the milder one.
- **Against truth.** Every scorer reads the simulator's per-fragment truth — the oracle BAM's read names, the oracle
  caches split from them, the truth tables beside each condition — never the previous run.
- **In a pinned, fractional pair.** An end-to-end arm runs with `--set scan.total_threads=1` (`OMP_NUM_THREADS=1` does
  not pin the scan) and `--set em.assignment_mode=fractional` (the shipped assignment is a sampled draw), so an A/B
  pair is exactly reproducible; `panel.py score` sets both. Each row records its assignment mode, and the report
  refuses to set two modes side by side. An effect is judged by its size, with genes and pools beside transcripts.
- **Against a fresh baseline.** Both arms of a comparison come from this session's tree, varying one thing
  (`TRAPS: re-record-the-baseline`).
- **Zero controls are controls.** The `g00` rung has no gDNA, so every gDNA fragment it reports is a false positive.
  It is read on every stratum and never tuned on: at truth zero, anything that lowers the estimate "improves" it
  (`TRAPS: zero-target-guards-are-one-sided`).

## Two primary numbers

| the question | the number | the instruments |
|---|---|---|
| is CALIBRATION right? — what ranks a calibration mechanism | the calibration result against an oracle calibration | `calibration_vs_oracle.py`, `prior_vs_oracle.py` |
| is THE TOOL right? — what the release ships on | the transcript table against per-transcript truth | `quant_accuracy.py`, through `panel.py score` |

Neither stands in for the other. The transcript table is downstream of calibration, the prior assembler, the
effective lengths and the EM, so it cannot say which of them moved, and ranking a calibration mechanism on it is
refused; a calibration figure cannot say whether the user's number improved.

## The calibration metric

| | what is scored | against | instrument |
|---|---|---|---|
| the result | `P`, the `CalibrationResult`'s six deconvolved arrays | `O`, the same result with only those arrays replaced by the origin-split truth | `calibration_vs_oracle.py` |
| the prior | `LocusPriors`: each locus's `gdna_count`, which the EM reads against its own fragment count | `O`, the same assembler fed the true per-object masses | `prior_vs_oracle.py` |
| one number | the library gDNA fraction | the per-fragment truth | `calibration_vs_oracle.py`, each row's `pools` block |
| the length | the capture-contracted length the EM divides by | the simulator's own capture-aware yield | `ruler_vs_truth.py --scale` |

`P − O` is calibration's own error and nothing else, because `O` is calibration done perfectly with the shipped
assembler. `prior_vs_oracle.py` also reports `O − S`, where `S` converts the crossing term by gDNA's true share
instead of the pooled one: the assembler's part, a different repair in a different file. `--set SECTION.FIELD=VALUE`
applies to both arms, so a config value is priced on the metric itself.

⛔ `O` keeps `P`'s capture efficiencies and the gDNA component's length, so it cannot see the capture-contracted length
or the efficiencies' own errors: read `ruler_vs_truth.py` beside every capture-ON arm. Its `--scale` read-out prints,
with no EM, L / Y for every hypothesis class (the locus gDNA component, the synthetic spans, the annotated
transcripts by probed class) and the within-gene spread of log(L / Y). Read both: the class means must sit on one
scale, because the gDNA-versus-RNA split reads their ratio, and the spread must not grow, because the isoform split
reads the ratios inside a gene (`TRAPS: judge-a-ruler-by-its-within-gene-spread`).

**The message policy** is scored on its own by `policy_benchmark.py`: each policy's whole-library gDNA error against
the certified per-slot truth (`slot_truth.npz`), per condition, with `--by-class` naming the node class where the
error sits. Unstranded rows must beat `silent`, the measured floor; stranded rows must do as little harm as possible
(`TESTING.md` §0a).

## The end-to-end number

`quant_accuracy.py` scores the transcript, gene and pool tables against each condition's realised fragment counts.
Every arm but `base` runs the identical pipeline with one thing substituted, so the difference from `base` is what
that one thing is worth:

| arm | substitutes | read |
|---|---|---|
| `base` | nothing: the shipped pipeline | the deliverable |
| `no_capture` | every capture efficiency read as 1; calibration's bulk kept | what reading capture is worth: a row where `base` is worse is a row where the capture reader cost more than it gave. It substitutes no truth, so it runs under `--bam` too |
| `oracle` | the `LocusPriors` built from the origin-split truth | what a perfect prior is worth; what is left belongs to the EM, the lengths and the annotation |
| `oracle_ruler` | the simulator's capture-aware effective length for every transcript and span, anchored on the fully probed transcripts | what the lengths do to the isoform split; the gDNA component keeps the shipped length |
| `oracle_alloc_seed` | the true relative abundances as the EM's per-transcript allocation weights, with the shipped coverage seed | a capability proof, never a headroom claim |
| `noop`, `oracle_ruler_noop` | nothing, through the substitution path | must be byte-identical to `base` |

Every prior arm wraps `assemble_priors`, which runs after the effective lengths are built, so only `oracle_ruler`
reaches the lengths; say which call an arm patches and check it sits downstream of everything you mean to price
(`TRAPS: a-symptom-is-not-a-second-defect`). `--bam` scores an aligner's BAM of the same reads against the same truth:
the ladder's STAR rung (`TESTING.md` §0d).

**The screen a change is read on** is every row there is, per stratum: the ladder's 16 conditions on the simulator's
alignments and on the STAR rung; the low-gDNA side panel and the fl-gap side panels (`TESTING.md` §0); and the test
chromosome's three probe designs at stranded capture-ON (`test_reference.yaml` and its two probe panels,
`TESTING.md` §0a). Every capture-ON row is read beside `no_capture`.

## 0.8.0 is done when

Each bar names its quantity, its strata and its truth; the numeric threshold on each is the owner's call.

1. **`P − O` is small on the three in-scope strata**, and the residual is attributed: to the assembler
   (`prior_vs_oracle.py`'s `O − S`), to the composition, or to a node class (`policy_benchmark.py --by-class`).
2. **The zero controls read zero**: the `g00` rung, on every instrument that reads it. A stratum can read healthy on
   every contaminated row and still claim gDNA in a library holding none.
3. **The capture-contracted length is right because the composition is**, never patched: at `g00`, where no gDNA
   measures capture, it reads 1 (`em_effective_length` equals `effective_length`), and elsewhere it tracks the
   simulator's yield (`ruler_vs_truth.py --scale`). A separate correction of the length is a defect, not a fix.
4. **The three in-scope strata improve, or at worst hold, on the transcript table**, in a pinned A/B pair decomposed
   by the arms, with the deferred stratum and the rest of the screen (the fl-gap panels, the low-gDNA panel, the
   probe designs) reported beside them.

The accumulator's own correctness (fidelity to the specification, worker-order independence, conservation) is a
closed regression check, gated in the suite (`TESTING.md` §5).

## The instruments, in the order you run them

Steps 0–1 build a panel once; after that a calibration change is re-scored in minutes, because neither cache depends
on `calibration/`. `python scripts/design/preflight.py --full` first, in a new session: it checks every instrument
imports and runs its `--self-test`.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=$RIGEL_RUNS/suite LADDER=$RIGEL_RUNS/suite/ladder INDEX=$RIGEL_RUNS/suite/rigel_index
CFG=scripts/sim/configs/gdna_ladder.yaml

# 0. WHERE AM I?  Names the next stage; every stage is resumable.
python scripts/sim/panel.py status --config $CFG

# 1. BUILD — once. `cache` builds both caches (scan and oracle), the step that makes the rest cheap.
python scripts/sim/panel.py build    --config $CFG
python scripts/sim/panel.py simulate --config $CFG --jobs 8
python scripts/sim/panel.py cache    --config $CFG --jobs 8

# 2. CALIBRATION AGAINST ORACLE CALIBRATION, per stratum — no EM.
python scripts/design/calibration_vs_oracle.py --suite $LADDER --index $INDEX --oracle-cache $LADDER/oracle_cache
python scripts/design/prior_vs_oracle.py --suite $LADDER --index $INDEX --oracle-cache $LADDER/oracle_cache --jobs 6
python scripts/design/ruler_vs_truth.py --panel ladder --condition gdna_g05_ss_0.99_nrna_mid_capture_on --scale
python scripts/design/policy_benchmark.py --panel ladder --policies silent transfer --by-class

# 3. THE TOOL END TO END — fractional, scan pinned; --jobs 2 (a 10 M-fragment condition holds ~8.5 GB).
python scripts/sim/panel.py score  --config $CFG --arms base no_capture oracle oracle_ruler --jobs 2
python scripts/sim/panel.py report --config $CFG --arms base no_capture oracle oracle_ruler
#    ...and on STAR alignments of the same reads (aligned once with `panel.py align`).
python scripts/sim/panel.py score  --config $CFG --arms base no_capture --bam rustar_recommended --jobs 2
#    ...and the rest of the screen: the same two arms on each side panel and probe design.
for C in gdna_ladder_lowg flgap_rna_long flgap_rna_short test_reference test_reference_probes_junction test_reference_probes_sparse; do
  python scripts/sim/panel.py score --config scripts/sim/configs/$C.yaml --arms base no_capture --jobs 2
done

# 4. A RENAME, REFACTOR OR SPEED-UP must be a numeric no-op against a frozen reference (TESTING.md §6).
python scripts/design/rename_identity.py --check --suite $LADDER --index $INDEX \
       --condition gdna_g05_ss_0.99_nrna_mid_capture_off --reference <frozen reference>
```

A scan change is decided on the STAR rung (`TESTING.md` §0d), and once the scan is final every panel is re-cached
(`panel.py cache --force`, which reaches `build_scan_cache.py`) before calibration is read. When dissecting rather than
scoring: run the panel → take the worst **in-scope** scenario → dissect it to the highest-error node class
(`policy_benchmark.py --by-class`) → find the mechanism → fix → re-run the panel. The worst scenario overall is in the
deferred stratum, and picking it is how the ranking gets quietly re-inverted.
