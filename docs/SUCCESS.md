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
capture-ON is deferred — reported on every table, never a development target — with the fragment-length
composition channel retired until after 0.8.0. The ruling, its reasons and the equal-fragment-length
forcing function are `DESIGN.md` §0b; every table in this file is read per stratum, never pooled
(`TRAPS: never-pool-the-strata`), because the deferred stratum carries most of the pooled error.

---

## The metric — the calibration result against oracle calibration

The primary number is calibration scored against an oracle calibration; the end-to-end transcript
number is a thermometer. The transcript table is downstream of calibration and of the EM, and a single
end-to-end figure cannot say which of the two moved.

| | what is scored | against | instrument |
|---|---|---|---|
| **primary** | the `CalibrationResult` itself: the six deconvolved arrays and the effective-length ruler derived from them | `O`, the same result with only the deconvolved arrays replaced by the origin-split truth; `U`, the no-enrichment null | `calibration_vs_oracle.py` — read `ruler_n_moved`, never the aggregate |
| **primary, the prior** | the prior calibration ships — `gdna_prior_count`, `rna_prior_count`, `gdna_eff_len` per multi-locus | `O`, the same assembler fed the origin-split truth masses | `prior_vs_oracle.py` (`P − O`) |
| **primary, per object** | each region's and boundary's own `f_g`, and whether it is confidently wrong | the oracle payload: the production accumulator run on the BAM split by true origin | `solvability_audit.py` |
| **primary, one number** | the library `f_gdna` | the simulator's per-fragment truth | `calibration_vs_oracle.py` — each row's `pools` block, `P_gdna` against `true_gdna` |
| **controls** | zero-gDNA and zero-RNA, where truth is a constant | 0.000 and 1.000 exactly | `zero_controls.py`, and the `g00` rung |
| **thermometer** | the transcript table a user reads | `truth_abundances.tsv` | `quant_accuracy.py --arm base` |

Why `P − O` and not the transcript number: attribution. `O` is calibration done perfectly with the
shipped assembler, so `P − O` is calibration's own error and nothing else; the transcript number adds
the assembler, the effective-length model, the EM's ambiguity and the annotation. `prior_vs_oracle.py`
reports `O − Fo` (the assembler's error) beside it, because they are different repairs in different files.
Read the thermometer, do not steer by it: most of the stranded × capture-OFF misassignment is ordinary
isoform ambiguity, and a calibration change that improves `P − O` and leaves the transcript table flat
has done its job. Its noise floor is measured, not assumed — `quant_accuracy.py --arm base_reseed`
prints it beside the effect and must be re-run in the same session (`TRAPS: re-record-the-baseline`).

### The ruler — the effective-length shrinkage sits outside every ceiling's patch point

`effective_lengths_em` is built inside `_setup_geometry_and_estimator` before `pipeline.py` calls
`assemble_priors`, and an arm that patches `assemble_priors` leaves the shipped shrinkage installed —
so a ceiling measured that way has never priced the ruler. `calibration_vs_oracle.py` substitutes at the
`calibrate` boundary and is the one instrument that reaches it. ⛔ Say which call your arm patches, and
check it sits downstream of everything you mean to price.

The ruler's defect is a symptom, not an independent bug: at zero gDNA the shrinkage contracts
transcripts when the correct factor is exactly 1 (even off capture, where its own contract says so),
because `rho_ref` is fabricated from false-positive gDNA — and substituting only the composition arrays
into the shipped shrinkage gives the truth factor. ⛔ Fix the composition and check that the factor
follows; a separate shrinkage repair would be half of a cancelling pair
(`TRAPS: a-cancelling-defect-pair`; `ISSUES: g00-shrinkage-upstream-repair`). A second lane is built and
not wired: the per-transcript RNA prior (`rna_prior_weight`) is never passed in production
(`ISSUES: per-transcript-prior-lane`).

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
| **P** | what pass-0 actually produces | `C − P` is the solver gap → the 0.8.0 work |

The decomposition is a measurement, not a build — the per-object form of the ceiling discipline
(`TRAPS: measure-the-ceiling-first`) — and it is subject to the ruler caveat above.

### How pass-0 is scored

Pass-0's job is not accuracy but to produce a substrate the gDNA hyperprior can be fitted against. An
object with no own evidence reporting `f_g ≈ ½` at zero precision is correct, so scoring every object
that carries mass counts honest ignorance as error. The measurement is a partition
(`scripts/design/solvability_audit.py`):

1. **Undetermined** — no own evidence; excluded from the error denominator. Its only failure mode is
   claiming a precision it has not earned, and that check must exist or the exclusion hides the largest
   error in the library (`TRAPS: excluding-a-population-hides-it`): `undetermined_overreach_rows`
   buckets the class by `|f_pred − ½|`; the correct answer for the class is ½ at `sd = ∞`.
2. **Solvable and right.**
3. **Solvable and wrong**, split by confidence. Confidently wrong is the defect — a wrong value with a
   tight variance outvotes correct neighbours and anchors the prior, so it propagates. The comparison is
   in log space (`var_gdna` is `Var(log f_g)`, `TRAPS: log-variance-is-not-linear`), and the headline is
   a calibration curve, which needs no threshold.

"No own evidence" is not a binary: the strand arm's information `I(f_g) ∝ (2κ−1)²` is exactly zero only
at κ = ½ while κ is fitted, so a threshold on `tau_lam` promotes objects to "solvable" whose own
statement has an `sd(λ)` of thousands of nats (`TRAPS: a-threshold-on-a-fitted-residue`). Strength is
therefore reported as a curve over `sd(λ) = 1/√τ` decades, and the panel table carries `weak%`, the
share of scored error above 10 nats; a better threshold was refuted, τ being continuous across the
region. Read `weak%` before `mwae`: a row near 100 is reporting the messages and the reference, not a
solve. `locked` is the structurally pure-gDNA class on both axes (`region_geometry.g1_locked`), never
`~solvable & is_region` (`TRAPS: two-masks-one-name`).

### The instrument — `scripts/design/pass0_vs_oracle.py`

It scans, builds T, runs its arms, and scores every one per object and per class; it also populates
the oracle cache every other scorer reads, which is why `panel.py cache` runs it (gates:
`tests/calibration/test_pass0_vs_oracle.py`). C is two ceilings, each defined by a lever that already
exists, never an estimator (`TRAPS: no-magic-numbers`):

| | what it is | what its gap to P means |
|---|---|---|
| **C_input** | `calibrate` handed the simulator's own post-capture length pmfs at both solve depths | how much of the error is wrong inputs rather than wrong solving |
| **C_info** | a classification, per object: is the 2×2 of `EQUATIONS.md` §3.1 identified from this object's own stored channels at all? | not a gap — C_info ignores neighbours and the sweep does not, so it can be "worse" than P |

C_input is a length-input ceiling and under the 0.8.0 scope a diagnostic, not a route (the other
library-level inputs are injectable via `InjectedCalibrationPriors` but the simulator writes no truth
for them; κ is not free either, `TRAPS: specificity-and-sense-are-complements`). ⛔ Do not read it as an
argument for the length channel, which is deferred by ruling.

The cross-tab is the point: objects undetermined by C_info *and* carried entirely by the messages have
no answer of their own, and that cell is reported with its mass share and its error share, needing no
threshold because it is a cell of a partition. The two classifications are each exhaustive (gated): the
solver's own — `own_evidence` / `message_only` / `struct_lock`, from
`region_init.has_own_composition_evidence` and `region_geometry.g1_locked` — and C_info's —
`identified` / `undet_no_separation` / `undet_out_of_range` / `absent`. Every arm, T included, is in
the DRAINED frame: draining three origin partitions separately is not the same operation as draining the
whole, so the partitions are lifted by replaying the whole's choices (`lift_drain_parts`) and the
sum-to-full identity is asserted on the drained frame.

---

## 0.8.0 is done when

Each row names the quantity that is the bar, on which strata, against what truth. The numeric threshold
on each is an owner call and is not invented here.

1. **`P − O` is small on all three in-scope strata**, and the residual is attributed — to the assembler
   (`O − Fo`), to the composition, or to a class that is provably undetermined. It is not done while the
   residual sits on objects `C_info` calls identified (`solvability_audit.py` re-derives the share).
2. **The zero controls read zero**: `zero_controls.py` on both arms, and the `g00` rung of the ladder.
   An in-scope stratum can read healthy on every contaminated row and still claim gDNA in a library
   containing none; only a zero control finds that.
3. **The effective-length shrinkage is correct because the composition is**, not because it was patched:
   at `g00` the factor reads 1.000, and elsewhere it tracks the truth factor when only the composition
   arrays are substituted. A separate shrinkage correction is a defect, not a fix.
4. **Pass-0 is monotone**: adding real evidence to an object never moves its answer away from truth.
5. **Pass-0 depends on no quantity a later iteration produces** — no feedback in the first solve.
6. **The three in-scope strata do not regress on the thermometer**, and the deferred stratum is reported
   on every table.

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
targets. Quote the shipped column, not pass-0, and never a pooled total
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

# 2. IS THE SUBSTRATE SOUND?  (TRAPS: prove-the-substrate — prove the simulator before the code)
python scripts/design/simulator_gates.py --suite $LADDER --reference $SUITE/reference
python scripts/design/suite_resolves.py $INDEX --suite $LADDER

# 3. THE PRIMARY METRIC — CALIBRATION AGAINST ORACLE CALIBRATION.
#    (a) the calibration result and the ruler, P vs O vs U, per stratum   ~5-12 s/condition, no EM
python scripts/design/calibration_vs_oracle.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache
#    (b) the PRIOR the EM actually reads, P vs O, per stratum   ~50 s/condition with the cache warm
python scripts/design/prior_vs_oracle.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache --jobs 6
#    (c) per OBJECT: solvable, solved wrong, and CONFIDENTLY wrong.  Read `weak%` before `mwae`.
python scripts/design/solvability_audit.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache

# 4. THE TWO ZERO CONTROLS — owner-required on EVERY experiment, both arms.
python scripts/design/zero_controls.py

# 5. THE THERMOMETER — the tool end to end, and the prior-injection ceiling above it.
#    --jobs 2, not more: run_pipeline holds 7-8.5 GB per 10 M-fragment condition.
#    `base_reseed` is the noise floor; any arm delta inside it is a sampling draw.
python scripts/sim/panel.py score  --config $CFG --arms base base_reseed oracle --jobs 2
python scripts/sim/panel.py report --config $CFG --arms base base_reseed oracle

# 6. STAGE A is CLOSED — this block is a REGRESSION check, run it after an accumulator or native change.
python -m pytest tests/native tests/calibration -q     # FIDELITY
python scripts/design/fl_pool_purity.py                #       are the pools pure gDNA (only where the lengths differ)
```

Steps 0 and 2–4 take about 15 minutes on a built panel; `suite_resolves.py`'s requirement (c),
replicate pairs, is the one known failure, deferred by owner ruling. Run the set together and record it
together (`TRAPS: re-record-the-baseline`). When dissecting rather than scoring: run the panel → take
the worst **in-scope** scenario → dissect it to the highest-error object (`worst_objects.py`) → find the
cause → fix → repeat. The worst scenario overall is the deferred stratum, and picking it is how the
ranking gets quietly re-inverted.
