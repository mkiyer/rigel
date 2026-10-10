# SUCCESS — how performance is judged, and what "done" means for 0.8.0

**What this file is.** How Rigel's performance is measured: the metric and its instruments, Stage A (the
accumulator — is the tally faithful, unbiased and sufficient?), Stage B (calibration scored against an
oracle calibration — the 0.8.0 work), the conditions under which 0.8.0 is done, which numbers are targets
and which are thermometers, and the instruments in the order you run them. Not here: the scope and
its rulings (`DESIGN.md` §0b), how the panels are built (`TESTING.md`), the ranked next steps
(`ROADMAP.md`), current measurements (run the named instrument — a figure that appears here is a
worked example, never current behaviour), and lessons (`TRAPS.md`).

---

## The scope

The version on disk is `pyproject.toml`'s; the target is 0.8.0, a calibration release. Three strata are
the target (unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON) and unstranded ×
capture-ON is deferred — reported on every table, never a development target — and fragment length is in
scope (owner, 2026-10-02). The rulings, their reasons and the equal-fragment-length
forcing function are `DESIGN.md` §0b; every table in this file is read per stratum, never pooled
(`TRAPS: never-pool-the-strata`), because the deferred stratum carries most of the pooled error.

---

## Two primary numbers, and what each one answers (owner, 2026-09-19: `DESIGN.md` §0b's amendment)

0.8.0 is a release of the TOOL, so the transcript table is a first-class number beside the calibration
metric rather than a thermometer under it. They answer different questions and neither stands in for the
other: the transcript table is downstream of both calibration and the EM, so a single end-to-end figure
cannot say which of the two moved — and a calibration figure cannot say whether the user's number improved.

| the question | the number | the instrument | how it is read |
|---|---|---|---|
| **is CALIBRATION right?** — the number that ranks a calibration mechanism | the `CalibrationResult` against an oracle calibration | `calibration_vs_oracle.py` | per stratum, never pooled |
| **is THE TOOL right?** — the number the release ships on | the transcript table against per-transcript truth | `quant_accuracy.py --arm base --set em.assignment_mode=fractional` | per stratum, by size, in a pinned A/B pair |

⛔ THE END-TO-END NUMBER IS READ UNDER FRACTIONAL ASSIGNMENT, IN A PINNED PAIR, WITH A DECOMPOSITION (owner,
2026-09-19 and 2026-09-28). The shipped EM assigns by a sampled draw, so every arm runs with `--set
em.assignment_mode=fractional`, which removes the draw from the comparison; each row records its mode and the
report refuses to set two modes side by side. An A/B pair runs with the scan pinned (`--set scan.total_threads=1`;
`OMP_NUM_THREADS=1` does not pin it), so it is exactly reproducible, and an effect is judged by its size, with genes
and pools read beside the transcript table. The arms decompose it: `oracle` is what a perfect prior is worth end to
end, so what remains under it belongs to the EM and the assignment rather than to calibration; every prior arm
wraps `assemble_priors`, so none reaches a transcript's length. `oracle_ruler` does: it hands the EM the
simulator's own capture-aware length for every transcript and synthetic span, anchored on the fully probed
transcripts, in place of the shipped one, and leaves the locus gDNA component on the shipped rule — so under
capture it prices end to end what the transcripts' lengths do to the isoform split, and not their scale beside
the gDNA component, which its anchor sets (`ISSUES: ruler-witness-geometry-on-transcript-panels`).

---

## The calibration metric — the result against oracle calibration

The number that ranks a calibration mechanism is calibration scored against an oracle calibration. Ranking a
calibration mechanism on the transcript table stays REFUSED for the reason above.

| | what is scored | against | instrument |
|---|---|---|---|
| **primary** | the `CalibrationResult` itself: the six deconvolved arrays | `O`, the same result with only the deconvolved arrays replaced by the origin-split truth (at capture-OFF it carries no enriched mode and is the no-enrichment null) | `calibration_vs_oracle.py`. `O` keeps `P`'s efficiencies, so it cannot see the capture-contracted length: read `ruler_vs_truth.py` beside every capture-ON arm |
| **primary, the prior** | what calibration hands the EM's prior — `gdna_count`, `gdna_eff_len` per multi-locus (the EM forms the pseudocounts from `gdna_count` and its own fragment count) | `O`, the same assembler fed the origin-split truth masses | `prior_vs_oracle.py` (`P − O`, the count; `gdna_eff_len`'s truth is `ruler_vs_truth.py`) |
| **primary, one number** | the library `f_gdna` | the simulator's per-fragment truth | `calibration_vs_oracle.py` — each row's `pools` block, `P_gdna` against `true_gdna` |
| **controls** | zero gDNA, where truth is a constant | 0.000 exactly | the `g00` rung of the ladder, read per condition by `calibration_vs_oracle.py` and `quant_accuracy.py` |
| **the deliverable** | the transcript table a user reads | `truth_abundances.tsv` | `quant_accuracy.py --arm base` |

Why `P − O` and not the transcript number: attribution. `O` is calibration done perfectly with the
shipped assembler, so `P − O` is calibration's own error and nothing else; the transcript number adds
the assembler, the effective-length model, the EM's ambiguity and the annotation. `prior_vs_oracle.py`
reports `O − S` (the assembler's pooled crossing share) beside it, because they are different repairs in different files.
Steer a CALIBRATION change by `P − O`, not by the transcript table: most of the stranded × capture-OFF
misassignment is ordinary isoform ambiguity, and a calibration change that improves `P − O` and leaves the
transcript table flat has done its job. The transcript table is what the RELEASE ships on, and it steers the
work downstream of calibration — the assembler, the ruler, the EM.

### The ruler — the capture-contracted length, and why a transcript's sits outside every prior arm's patch point

`effective_lengths_em` is built inside `_setup_geometry_and_estimator` before `pipeline.py` calls
`assemble_priors`, and an arm that patches `assemble_priors` leaves the shipped shrinkage installed —
so a ceiling measured that way has never priced the ruler. The length's truth is `ruler_vs_truth.py`: the
simulator's own capture-aware yield, the one instrument that scores a length against what generated the reads
(`docs/TESTING.md` §0c). Its `--scale` read-out is how the length is judged, with no EM: L / Y on the shipped
lengths for every component — the locus gDNA component, the synthetic spans, and the annotated transcripts by
probed class with the junction-probed ones apart — and the within-gene spread of the annotated transcripts'
log(L / Y). Read both: the class means must sit on one scale, because the gDNA-versus-RNA split reads the
ratio of the gDNA component's length to the RNA's, and the within-gene spread must not grow, because the
isoform split reads the ratios inside a gene and a repair of the scale can cost it
(`TRAPS: judge-a-ruler-by-its-within-gene-spread`); then price the change on the transcript table, per stratum.
⛔ Say which call your arm patches, and check it sits downstream of everything you mean to price.

Every EM component's length is one shared rule: the component's conserved share of each region and boundary
its fragments deposit on (a junction, in a spliced transcript's own coordinates) times that object's capture
weight (`DESIGN.md` §7.2, `EQUATIONS.md` §11). An object's weight is the posterior mode of its own gDNA density
under the fitted gDNA landscape, read inside the block solve from its strand columns with the RNA amount
integrated out and the factors delivered to it, relative to the typical read object's; a junction, where gDNA never
deposits, is priced from the objects within one fragment of it and capped at the most captured of them. There
is no detector and no reference: a capture-OFF library is read the same way, and the Poisson noise its
weights carry is priced per stratum (`ISSUES: the-capture-weights-are-noisy-on-a-capture-off-library`). A
second lane is built and not wired: the per-transcript RNA prior (`rna_prior_weight`) is never passed in
production (`ISSUES: per-transcript-prior-lane`).

---

## STAGE A — the accumulator: done, and that is a measurement

Kept as a record and a regression check, not as work. Stage A asks whether the information is *there*;
the 0.8.0 work asks whether calibration *finds* it. The three criteria: **fidelity** — the tally
reproduces the specification exactly; **bias** — no channel is systematically off against per-fragment
truth; **sufficiency** — perfecting the stored length information changes nothing downstream. Do not chase the gDNA length model's residual under
capture: its cause is known and is not a divisor (both opportunity functions assume uniform placement
and capture does not).

The five identities that keep it done, all gated:

| | | |
|---|---|---|
| C++ vs the executable specification | byte-identical | `tests/native/test_accumulator_spec.py` |
| the same BAM at 1/2/4/8 workers | bit-identical, every bank | `tests/test_scan_order_independence.py` |
| `Σ node_start_count == deposited` | exact | same |
| `deposited + deferred + dropped_* == offered` | exact, deferred non-empty | same |
| the origin partitions sum to the full payload | exact, every channel | `tests/calibration/_oracle.py` |

`ReadSimConfig.r1_sense` selects the protocol direction independently of `strand_specificity`, gated
both ways in `test_strand_sense_convention.py`; the panel does not expose the knob (`TESTING.md` §1).

---

## STAGE B — calibration: the 0.8.0 work

Calibration is iterative; **pass-0** is the prior-free solve, the first thing after initialisation,
before any fitted prior exists. It is the right place to start because it has no feedback loop in it —
every later iteration is conditional on pass-0 being sane, and `TRAPS: variance-fitted-on-the-belief` is
what happens when a prior is fitted on an unsound belief.

### The oracle, and what it buys

`tests/calibration/_oracle.py` produces the production accumulator run on the BAM split by true origin,
so every region and boundary has its true gDNA and RNA counts on each strand. Three quantities follow,
and the differences between them are the whole diagnostic:

| | quantity | what a gap to the next one means |
|---|---|---|
| **T** | the truth: each object's real gDNA/RNA split | — |
| **C** | a ceiling: the best answer reachable under stated conditions | `T − C` is information the accumulator destroyed → Stage A work |
| **P** | what calibration actually produces (`calibration_vs_oracle.py`'s `P`) | `C − P` is the solver gap → the 0.8.0 work |

The decomposition is a measurement, not a build — the per-object form of the ceiling discipline
(`TRAPS: measure-the-ceiling-first`) — and it is subject to the ruler caveat above.

### The oracle arms — `scripts/design/_oracle_arms.py`

A helper, not an instrument: the oracle cache builder `calibration_oracle.py --build` uses (which is why
`panel.py cache` runs that) and the per-object scorer `calibration_vs_oracle.py` reads. T is in the DRAINED
frame: draining three origin partitions separately is not the same operation as draining the whole, so the
partitions are lifted by replaying the whole's choices (`lift_drain_parts`) and the sum-to-full identity is
asserted on the drained frame (gates: `tests/calibration/test_oracle_arms.py`).

---

## 0.8.0 is done when

Each row names the quantity that is the bar, on which strata, against what truth. The numeric threshold
on each is an owner call and is not invented here.

1. **`P − O` is small on all three in-scope strata**, and the residual is attributed — to the assembler
   (`prior_vs_oracle.py`'s `O − S`), to the composition, or to a node class (`policy_benchmark.py --by-class`).
2. **The zero controls read zero**: the `g00` rung of the ladder, on every instrument that reads it.
   An in-scope stratum can read healthy on every contaminated row and still claim gDNA in a library
   containing none; only a zero control finds that.
3. **The effective-length shrinkage is correct because the composition is**, not because it was patched:
   at `g00` the factor reads 1.000, and elsewhere it tracks the simulator's truth, read with `ruler_vs_truth.py` —
   the oracle calibration keeps `P`'s efficiencies and cannot see it. A separate shrinkage correction is a defect,
   not a fix.
4. **The three in-scope strata improve, or at worst hold, on the DELIVERABLE** — the transcript table in a
   pinned A/B pair, decomposed by the arms so that what moved is attributable — and the deferred stratum is
   reported on every table.

The gate that held this work back, kept because it will apply again: a solve tuned against a wrong
input on exactly the conditions it is meant to rescue is tuned against a manufactured discriminant
(`TRAPS: prove-the-substrate`). Its live form is the ruler above.

---

## Reporting numbers vs target numbers

The target is the calibration result against the oracle calibration, per stratum, on the three in-scope
cells. The library gDNA fraction and the transcript table are the products, reported every session so
their trajectory is visible, and neither is the steering wheel: the library figure mixes accumulator,
solver and unidentifiability into one number, and the transcript figure adds the EM and the annotation.
Score the contaminated conditions — zero-gDNA rows are saturated at truth = 0, so anything that lowers
the estimate "improves" them (`TRAPS: zero-target-guards-are-one-sided`); they are controls, never
targets. Quote the shipped result, and never a pooled total
(`TRAPS: the-intermediate-is-not-the-deliverable`).

---

## The instruments, in the order you would run them

Steps 0–1 build the scenarios once; after that a calibration change is re-scored in minutes, because
both caches survive it (neither depends on `calibration/`).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite
LADDER=$SUITE/ladder
INDEX=$SUITE/rigel_index
CFG=scripts/sim/configs/gdna_ladder.yaml

# 0. WHERE AM I?  Every stage is expensive and resumable, and this names the next one.
python scripts/sim/panel.py status --config $CFG

# 1. BUILD THE SCENARIOS — once.  `cache` builds BOTH caches and is the step that makes the rest cheap.
#    The reference carve is manual: it needs the SOURCE genome/GTF, which a panel config does not name.
python scripts/sim/panel.py build    --config $CFG
python scripts/sim/panel.py simulate --config $CFG --jobs 8
python scripts/sim/panel.py cache    --config $CFG --jobs 8

# 2. THE PRIMARY METRIC — CALIBRATION AGAINST ORACLE CALIBRATION.
#    (a) the calibration result, P vs O, per stratum        ~5-12 s/condition, no EM
#        It cannot see the capture-contracted length: read (c) beside every capture-ON arm.
python scripts/design/calibration_vs_oracle.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache
#    (b) the PRIOR the EM actually reads, P vs O, per stratum   ~40 s/condition with the cache warm
python scripts/design/prior_vs_oracle.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache --jobs 6
#    (c) the CAPTURE-CONTRACTED LENGTH the EM divides by, against the simulator's own yield, no EM: every
#        component's class mean on one scale AND the within-gene spread — read both.
for C in gdna_g05_ss_0.99_nrna_mid_capture_on gdna_g50_ss_0.99_nrna_mid_capture_on; do
  python scripts/design/ruler_vs_truth.py --panel ladder --condition $C --scale
done

# 3. THE NUMBER THE RELEASE SHIPS ON — the tool end to end, with the ceiling arms above it; the `g00` rows are the
#    zero controls, read per condition.
#    `panel.py score` reads every arm under fractional assignment (the protocol above) with the scan pinned.
#    --jobs 2, not more: run_pipeline holds 7-8.5 GB per 10 M-fragment condition.
python scripts/sim/panel.py score  --config $CFG --arms base oracle oracle_ruler --jobs 2
python scripts/sim/panel.py report --config $CFG --arms base oracle oracle_ruler

# 4. STAGE A is CLOSED — this block is a REGRESSION check, run it after an accumulator or native change.
python -m pytest tests/native tests/calibration -q     # FIDELITY
```

Steps 0 and 2(a)–(c) take about 15 minutes on a built panel. Run the set together and record it
together (`TRAPS: re-record-the-baseline`). When dissecting rather than scoring: run the panel → take
the worst **in-scope** scenario → dissect it to the highest-error node class (`policy_benchmark.py --by-class`) → find the
cause → fix → repeat. The worst scenario overall is the deferred stratum, and picking it is how the
ranking gets quietly re-inverted.
