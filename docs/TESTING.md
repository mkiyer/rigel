# TESTING — the panels, the gates, and what the suite can judge

**What this file is.** The manual for Rigel's test substrates: the simulated panels and what each can
judge (§0), how the test chromosome is built, swept and read (§0a), the toy harness (§0b), the
simulator's backbone and its gaps (§1), the build commands (§2), the simulator's gates (§3), how results
are evaluated (§4), the test suite (§5) and the profiling substrate (§6). Not here: the 0.8.0 scope, the equal-fragment-length ruling and the nascent scope ruling
(`DESIGN.md` §0b); how performance is judged and the instruments' run order (`SUCCESS.md`); the suite's
standing pass count (`CLAUDE.md`); lessons (`TRAPS.md`, cited by name); open problems and refusals
(`ISSUES.md`). The panels' history is git.

Two rules govern the whole file (`TRAPS: prove-the-substrate`, `TRAPS: can-the-benchmark-resolve-it`):

> ⛔ Prove the substrate before you prove the code: when a simulated axis is the axis you are judging,
> gate the simulator on it. And prove the suite can resolve the axis you are changing before quoting a
> number from it.

---

## 0. The panels — one ladder, two fl-gap side panels, and a toy harness

Everything simulated lives under `~/Downloads/rigel_runs/suite/` (the ladder, the two side panels, the
carved reference and the index) and `~/Downloads/rigel_runs/test_reference/` (the test chromosome, §0a).
The panel says how much error there is and where; a toy (§0b) says why.

**The ladder** (`scripts/sim/configs/gdna_ladder.yaml`, `ladder/`) is the only panel the tool is ranked on.

| | |
|---|---|
| judges | Stage B, the calibration solver (`SUCCESS.md`) |
| conditions | 16: gDNA `g00 / g05 / g50 / g98` × strand specificity `0.50 / 0.99` × capture off / on |
| depth | 10 M fragments total per condition, fixed; the gDNA rate decides only the split |
| fragment lengths | identical for both origins (206 ± 98, truncated to [50, 500], read 100) — a forcing function, `DESIGN.md` §0b |
| nascent RNA | sparse on every row: `mode: sparse`, `on_fraction 0.50` of gene spans on, level logU(1, 100) independent of the mature level. 0.50 is a development stress level, not real data (`DESIGN.md` §0b) |

The 0.8.0 scope — three strata in scope, unstranded × capture-ON deferred but reported on every table,
the length composition channel retired, and why the ladder gives gDNA and RNA equal fragment lengths —
is ruled in `DESIGN.md` §0b. Report per stratum, never pooled (`TRAPS: never-pool-the-strata`); the
`g00` rung is the false-positive check on every stratum, and a control is never tuned on.

### The gDNA rungs

`f_gdna = rate/(1+rate)`, and the RNA share thins as gDNA rises because the total is held:

| rung | g00 | g05 | g50 | g98 |
|---|---|---|---|---|
| rate | 0.0 | 0.052632 | 1.0 | 49.0 |
| f_gdna | 0 | 0.05 | 0.50 | 0.98 |
| n_rna | 10.0 M | 9.50 M | 5.00 M | 0.20 M |

`g00` is the required zero-gDNA control, `g98` the top of the range, and `g05` exists because
real libraries live at 1–10 % gDNA and that corner needs a rung with capture on (`0 < rate ≤ 0.10`). Three
levels is a floor (`TRAPS: a-single-level-panel-cannot-see-a-constant`). Real libraries run from almost
zero gDNA to over 98 %, and the RNA-side accuracy that thins at the top rungs is a property of such
libraries, not an artefact.

### The fl-gap side panels — two arms of opposite sign

A side panel, never a ladder rung: the ladder equalises fragment lengths, so it is a null for anything
whose mechanism is a gDNA-vs-RNA length difference, and these two arms exert exactly that. RNA is not
reliably longer than gDNA (true for cfRNA, false elsewhere), and `E_r/E_g − 1` flips sign between the
arms, so a real repair must move them in opposite directions.

| | `flgap_rna_long/` | `flgap_rna_short/` |
|---|---|---|
| config | `scripts/sim/configs/flgap_rna_long.yaml` | `scripts/sim/configs/flgap_rna_short.yaml` |
| gDNA / RNA fl configured | 75 ± 20 / 250 ± 60 | 250 ± 60 / 75 ± 20 |
| gDNA / RNA fl measured | 78.58 / 247.62 (gap +169 bp) | 249.59 / 78.43 (gap −171 bp) |
| conditions | `g50` × ss {0.50, 0.99} × capture {off, on} | the same 4 |
| nascent RNA | `mode: fragment_share`, `shares: [0.20]` | the same |

The configured lengths are a configuration and the gap is a measurement (the sampler is a rejection draw
truncated to `[frag_min, frag_max]`): read it off the payload's deposit histograms. What may be read off these arms is everything that stops
before the EM (`calibration_vs_oracle.py`) and the library gDNA fraction; the
transcript-level number is not a calibration result here, because a length gap hands the EM the answer.

⛔ **Both side panels carry a different nascent model from the ladder** — `fragment_share` at a flat
0.20 on every expressed multi-exon span against the ladder's `sparse`. Each panel's data matches its own
config, so every measurement taken on either arm stands, but a ladder-vs-side-panel comparison varies
two things and no claim may be carried across them (`ISSUES: flgap-panels-stale-nascent-model`).

### The third fl panel — on the test chromosome, on the current nascent model

| | `scenarios_fl_rna_long/` | `scenarios_fl_gdna_long/` | `scenarios_fl_equal200/` |
|---|---|---|---|
| config (`scripts/sim/configs/`) | `test_reference_fl_rna_long.yaml` | `test_reference_fl_gdna_long.yaml` | `test_reference_fl_equal200.yaml` |
| gDNA / RNA fl configured | 100 / 250 | 250 / 100 | 200 / 200 — the control |
| conditions | the full 30 of §0a | the same 30 | the same 30 |

These carry the current nascent model structurally: the test reference is `abundance.mode: file`, so
each config's `nrna:` block is dead by design and nascent comes from `test_abundances.tsv`. The suite
arms' stale nascent model therefore does not block work on the fl gap. The equal-200 arm is the control:
an fl-gap result that also moves on the equal-length arm is an artefact.

### Every scenario must be cached — a requirement of the 0.8.0 loop

The loop is *change calibration → re-run calibration → score it against oracle calibration*, and nothing
in it may re-scan a BAM. Two caches carry it and `panel.py cache` builds both:

| cache | holds | invalidated by |
|---|---|---|
| scan (`build_scan_cache.py`) | the accumulator payload — scan once, calibrate many times | any accumulator change (`payload_schema_digest`), the index (`graph_hash`, `reach_digest`) |
| oracle (`<panel>/oracle_cache/`) | the origin-split truth every scorer reads: five partitions per condition (`gdna` / `mrna` / `nrna` and the per-strand `rna_pos` / `rna_neg`), the undrained `_main` payload, and the certified `slot_truth.npz` | the accumulator or the index, and nothing else |

Neither is invalidated by a calibration change. A scenario without both caches is not usable for
development: `panel.py status` names what is missing, and an instrument fed a stale cache refuses it
(the oracle cache is keyed by the scan cache's own key). The toy harness's donor bundle is the deliberate
exception (§0b): it is a function of the calibration code that fit it, so caching it would serve a stale
answer. How `status` counts a `g00` row is a §2 gotcha.

---

## 0a. The method-development test chromosome — where the message policy is developed

A synthetic chromosome the owner designs one structure at a time, described in one YAML file, swept over
30 conditions, cached, and scored in seconds — so a policy change is scored on every condition between
edits, and the ladder is kept for the shipping judgement.

### The sweep — 30 conditions

`scripts/sim/configs/test_reference.yaml` sweeps three axes; two more requirements are carried by the
substrate rather than the config, so that no condition is uniform in either:

| axis | values | where it lives |
|---|---|---|
| gDNA fraction | `g00` `g05` `g25` `g50` `g98` | the config (`rate = f/(1−f)`) |
| strand specificity | `0.50` unstranded · `0.70` · `0.99` strand-specific | the config |
| capture | off · on | the config |
| which transcripts are captured | about half probed, half with no probe | `test_probes.bed` |
| nascent RNA | sparse — up to about half the transcripts carry it, the rest exactly none, at levels below, near and above mature | the `nrna_abundance` column |

`g00` is the zero control (every reported gDNA fragment is a false positive) and `ss 0.70` is the
transition rung between the two regimes the bars below are written for. The nascent pattern is
hand-authored, not drawn: on a chromosome this small a draw is mostly sampling noise, and the point of a
debug substrate is knowing which structure carries nascent RNA before reading a number. Depth is
`n_total_fragments: 1170000` per condition, raised with each block so per-transcript depth holds. The
config's `nrna:` block is dead by design (`abundance.mode: file` makes the rendered TSV the one source).

### The commands — from the YAML to a scored benchmark

Everything derived is rebuilt from the one YAML. Run in order; each stage is resumable and
`panel.py status` names the next one. The `cache` stage builds both caches, pre-warms the `g00` rows
(held out of the oracle sweep) and copies their `_main` from the scan cache, then certifies
(`calibration_oracle.py` writes `slot_truth.npz`; a failed FIELD gate is reported, not fatal — the row is
COMPOSITION-certified). One row by hand is a §2 gotcha.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
T=~/Downloads/rigel_runs/test_reference

# 1. renders (GTFs, abundances, probe panels — into the repo AND the runs dir) + the two-contig FASTA;
#    `--check` alone says whether the checked-in renders still match the YAML
python scripts/sim/build_test_reference.py
rigel index --fasta $T/test_chr.fa --gtf $T/test_chr.gtf --collapse-duplicate-transcripts --no-mappability --no-tsv -o $T/idx

# 2. per panel (seven configs: the benign panel, the two adversarial probe panels, od05, three fl arms)
#    ⛔ ONE panel's `cache` at a time, or give each its own RIGEL_SCRATCH: the origin split writes
#    $RIGEL_SCRATCH/rigel_pass0_oracle/<condition>.<origin>.bam, and six of the seven configs share their
#    condition names, so two caching side by side corrupt each other's partition (sum-to-full fails).
CFG=scripts/sim/configs/test_reference.yaml
python scripts/sim/panel.py simulate --config $CFG
python scripts/sim/panel.py cache    --config $CFG

# 3. score every policy on every condition — seconds
python scripts/design/policy_benchmark.py --panel test
```

⛔ **Rebuild after every edit to the YAML.** The FASTA, index, reads, caches and `slot_truth.npz` all
describe the annotation they were built from; a benchmark against stale caches answers a different
question (`TRAPS: a-green-suite-hid-five-dead-instruments` in another dress). Move a superseded derived
set out of the way the same day.

### How to read the benchmark — two halves, two bars

The bar is not "beat silent". `SilentPolicy` is the measured floor, and on strand-specific data a
sighted exon's own strand solve is excellent, so a message can mostly only disturb it. Message
propagation exists for unstranded data, where the strand channel is dead and the local answer is a
default rather than a measurement.

| rows | the bar |
|---|---|
| unstranded (`ss 0.50`) | the policy must win — this is what the message layer is for |
| stranded (`ss 0.99`) | the policy must do as little harm as possible — near 1.00× of silence is a pass |

⛔ Never pool the two halves, and never pool conditions within them: a total hides a sign flip between
strata, and the halves are judged against different bars. A toy and the panel can disagree in rank
(`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`): develop on the test chromosome, confirm on the ladder.

### The substrate — one hand-edited file, twelve blocks

`scripts/sim/test_reference/test_chr.yaml` is the one hand-edited file: the `rigel sim` scenario schema
(`genes → {gene_id, strand, transcripts: [{t_id, exons, abundance, nrna_abundance}]}`, exons 0-based
half-open) plus `probed` per gene and a `shadow_genes` list. Its header is the design record of every
block; this table says only what each one stresses.

| block | what it stresses | gene types |
|---|---|---|
| anchored twin | the message layer's five controls on one shape (exon 1 kb · intron 7 kb · exon 1 kb · intron 7 kb · exon 1 kb), each type varying one bit — probed, nascent, silent — × five abundance blocks; plus 8 shadow transcripts the index never sees | `clean` `nasc` `cap` `capnasc` `silent` |
| mono | single-exon: no sj, the gene edge is the only message source | `mono` `capmono` + two silent controls |
| isoform | host + one second isoform, grown one structure at a time, every structure replicated across A ≫ B, A ≪ B, A ≈ B | `altstart` `altss` (`nest` queued, `ISSUES: message-layer-open-cases`) |
| walled | exon pieces with no licensed face — the scan's stress test, the four classes the ladder's walled-exon census found | `chain` `tssalt` `tandem` `altlast` |
| terminus-cluster | ten transcript ends 126–147 bp into a shared last exon (mirrored from MIR99AHG) — the empty exon pieces the level lane crosses | `cluster` `capcluster` |
| both-stranded | two genes per locus on opposite strands, the host and its antisense (mirrored from TTC28-AS1, PPM1F-AS1, and a convergent pair) — the AMBIG nodes' substrate | `asin` `asinrev` `span` `conv` |
| sj+terminus | one boundary carrying a junction and a terminus of the same strand (mirrored from RUNX1 and LARGE1) | `sjterm` `capsjterm` |
| encompassing | a 20 kb single-exon antisense over the whole twin-shape host (mirrored from ENSG00000280007 over TUBA8): every host slot admits both strands, every boundary carries only the host's bits, the antisense's level lives in its two 1.5 kb single-strand flanks — the owner's encompassing locus on the panel | `enc` `capenc` |
| in-exon | a 500 bp single-exon antisense wholly inside the host's 3 kb last exon (mirrored from ENSG00000273300 in UFD1's 3' UTR): no junction and no single-strand piece of its own, the tilt atom's accepted limit | `inexon` `capinexon` |
| shared-exon | two spliced genes sharing one exact 10 kb last exon (the W12 deep stress; the tail-to-tail class of SMARCB1 × DERL3): a level on each strand into one walled AMBIG exon whose faces each pair a junction with a terminus | `shared` `capshared` |
| in-intron | a 600 bp single-exon gene centred in the opposite strand's first intron (mirrored from the ladder's 203 intronic pseudogenes and lncRNAs): the majority AMBIG class, walled, its only message source its own termini | `inintron` `capinintron` |
| head-to-head | `conv` with the strands swapped and nothing else: two 5' ends overlap in a 2 kb exon∩exon piece (mirrored from TRMT2A ⟷ RANBP1), a TSS and a donor of different strands on one boundary | `div` `capdiv` |

265 genes on a 7.637 Mb chromosome (`genome_length` in the YAML). Every gene carries an explicit strand
and the chromosome keeps equal + / − representation — a sign error is invisible on one strand, and the
builder refuses an imbalance. Abundances are molar ladders in half-decade steps with mature up the
blocks and nascent down (10/30/100/300/1000 against 100/30/10/3/1), independent levels, as the ladder
draws them.

`shadow_genes` sit on the blank contig `test_blank`: unannotated transcription the simulator draws from
and the index never sees — the control for anything that trusts the annotation's notion of "pure gDNA".
The GTFs, `test_abundances.tsv` and the three probe BEDs are rendered by `build_test_reference.py` and
versioned beside the YAML (`--check` and `tests/test_test_reference_renders.py` refuse a drifted render);
the benign panel tiles each probed gene's exon union 8 × 125 bp, never spanning an sj, the sparse panel
centres one probe per exon, the junction panel places one two-block probe per sj. The FASTA
(`test_chr.fa`, both contigs), the index, the reads and the caches are derived, never hand-edited: a
spliced transcript needs a GT..AG at every intron or the aligner and the simulator disagree with the
annotation, and generating the chromosome from a fixed seed with a motif injected at every declared
intron makes that impossible to get wrong. The builder refuses a reference it cannot simulate and reports
every problem at once (an exon off the end, overlapping exons, an intron under 4 bp, a duplicate id, an
abundance row without a transcript or the reverse, the wrong reference name); `--self-test` perturbs each.

⚠ κ is fitted from spliced reads, so a transcriptome of only single-exon transcripts cannot calibrate
(`CalibrationStrandError`): at least one multi-exon transcript with real depth is needed before a number
means anything. Nascent RNA needs no declaration: `rigel index` creates a single-exon nascent entity
over each multi-exon transcript, fed from the contributor's `nrna_abundance` column.

## 0b. The toy harness — a mini chromosome you define, calibrated in under a second

`scripts/design/toy_harness.py` is the third substrate, for isolating one mechanism on geometry you fully
control: small enough to read every object, cheap enough to sweep one variable seven times. Prefer §0a
for developing a policy — a toy is too small to fit the library-level quantities calibration needs, so
it harvests them from a panel condition (`DonorGlobals`, ~30 s per session), and a toy result is
conditional on its donor; the test chromosome is a real simulated library of its own.

### The one idea that makes it work

A real cached condition acts as donor: it is calibrated once and its fitted bundle is injected
(`InjectedCalibrationPriors`). The toy supplies only the controlled per-region geometry.

| the donor supplies | so the toy never invents |
|---|---|
| κ and the spliced sample behind it, both strand overdispersions | the strand channel's protocol decision (`region_init.strand_discriminability`) reads as it does on real data — a gDNA-free donor keeps its channel |
| the intron background and the pre-solve total-density landscape (`abundance_landscape`) | a handful of regions cannot fit these |
| both fragment-length pmfs | passed as `calibrate` kwargs, not part of the priors bundle |
| capture on/off and its numeric knobs | reproduced in the toy's own simulation, with probes written from the spec |
| frag mean/sd/min/max, read length, strand specificity | read from the donor's post-capture truth |
| gDNA density per base | derived, not chosen — below |

⛔ The gDNA level is derived, not chosen: the injected landscape is an absolute log-density model, so a
toy at the wrong depth is a different library, not a small one. The harness measures the donor's gDNA
counts-per-base on its own structurally pure-gDNA regions (`Σcount / ΣE`, `EQUATIONS.md` §7.1) and
simulates the toy to match. `ToySpec` has no gDNA field, and a gate asserts it never grows one.

Vocabulary is `DESIGN.md` §0 (counts; density = abundance, counts per base). A region's stored counts
are contained counts, so `counts = density × effective_length`, never `density × bp`.

### Running it

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite

python scripts/design/toy_harness.py --list            # the spec ladder, simplest first

# one spec against one donor — prints EVERY object beside per-object truth   (~0.1-5 s)
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g50_ss_0.50_nrna_mid_capture_off

# sweep the transcript's RNA density; the gDNA background stays pinned, so one variable moves
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g50_ss_0.50_nrna_mid_capture_off --sweep-density

python scripts/design/toy_harness.py --spec all --donor <cond>     # the whole ladder
```

Which ladder condition is the donor is an experimental variable: it sets capture, the strand regime and
the gDNA level, so it selects a stratum (`ss_0.50 × capture_on` is the deferred one). Harvesting is
deliberately not cached; harvest once per session and run many toys against it (`harvest()` then
`run_toy()` in a loop).

### `spliced_exons` — the rung to understand first

One gene, one transcript, `TA+ (1,000, 2,000) (9,000, 10,000)` on 12 kb — `nested_exons`' twin, with an
intron and an sj where the nesting was, so the two rungs differ by exactly one structure.

```
REGION intergenic [0, 1000)         BOUNDARY @1,000   intergenic|exon, pure gDNA (TSS+)   G1 object
REGION exon  [1000, 2000)   TA e1   BOUNDARY @2,000   intron|exon, the DONOR+ side        NOT a G1 object
REGION intron [2000, 9000)  TA i1   BOUNDARY @9,000   intron|exon, the ACCEPTOR+ side     NOT a G1 object
REGION exon  [9000, 10000)  TA e2   BOUNDARY @10,000  intergenic|exon, pure gDNA (TES+)   G1 object
REGION intergenic [10000, 12000)
SJ BOUNDARY 2,000 → 9,000 (+), pure mature RNA, NOT a chain slot
```

The hard part is the two exon|intron boundaries. Mature RNA cannot cross one contiguously
(`TRAPS: mature-rna-never-crosses-a-boundary`), but the solver's continuity gate says a strand is
admissible there (RNA that has not spliced *there* could cross), so they are not G1 objects
(`DESIGN.md` §0) and the solver has to derive what the structure implies. The ladder's sparse nascent
model holds both faces at once: introns with nascent RNA in them and introns whose truth is exactly pure
gDNA. Unlike `nested_exons` there is own evidence inside the gene — the 7 kb intron is where the intron
factory lives — while the two exons have essentially no own evidence on an unstranded library and are
carried by messages. The derivation this rung exists to land is `EQUATIONS.md` §3.6.

### `splice_both_strands` — the both-stranded rung

Four transcripts, both strands, overlapping exons and overlapping introns, two sj pointing opposite ways,
on the same 12 kb chromosome:

```
TA+ (2,000, 3,000) (9,000, 10,000)     2 exons, + strand, intron 3,000–8,999
TB+ (2,000, 10,000)                    1 exon,  + strand, spans TA's intron
TC− (1,000, 11,000)                    1 exon,  − strand, spans everything
TD− (1,000, 2,500) (8,500, 11,000)     2 exons, − strand, intron 2,500–8,499
```

Every earlier rung let a boundary ask "is my neighbour an exon?" and get a yes or a no. Here the 4-bit
signature the index stores is `{intron₊, intron₋, exon₊, exon₋}` and three regions carry both kinds —
[2,500, 3,000) `0111`, [3,000, 8,500) `1111` (all four), [8,500, 9,000) `1011`; the rest are `0001` or
`0011`. `coarse_type_array` reports every one of them as `exon` (exon wins, the strand is collapsed), so
on this rung the string `intron|exon` never appears and any rule phrased on the coarse type is silent
exactly where it is needed. The sj axis (sj 0: 2,500 → 8,500 on −; sj 1: 3,000 → 9,000 on +) and the
boundary flags carry the rest:

| boundary | 1,000 | 2,000 | 2,500 | 3,000 | 8,500 | 9,000 | 10,000 | 11,000 |
|---|---|---|---|---|---|---|---|---|
| flag | `TES_NEG` | `TSS_POS` | `DONOR_NEG` | `DONOR_POS` | `ACCEPTOR_NEG` | `ACCEPTOR_POS` | `TES_POS` | `TSS_NEG` |

⚠ Check the `_NEG` convention before relying on it. Boundary 2,500 is flagged `DONOR_NEG`, but on a −
transcript the molecule runs right-to-left, so TD−'s biological donor is at 8,500. Whether the bits mean
"genomic-low end of a − intron" or "the transcript's actual donor" decides the sign of a derivation;
`EQUATIONS.md` §3.5b rules that this family of predicates is written in genomic terms.

### Writing a new spec

Add a `ToySpec` to `SPECS` in the harness. Each rung adds exactly one structure to the one before it, so
when a row goes wrong the thing that changed is the thing to look at.

```python
"my_case": ToySpec(
    name="my_case",
    what_it_probes="one sentence: which mechanism this isolates, and why this structure isolates it",
    genome_length=5_000,
    genes=[{"gene_id": "TA", "strand": "+",
            "transcripts": [{"t_id": "TA", "exons": [(1_000, 3_000)], "abundance": 100.0}]}],
    n_rna_fragments=1_000,      # the RNA knob; gDNA is pinned by the donor
    nrna_abundance=0.0,         # nascent RNA; the ONLY way an intron carries RNA
    captured=None,              # transcript ids to probe when the donor is capture-ON; None = all
    seed=7,
),
```

There is no gene-free rung (`TranscriptIndex` requires a transcript). Use a silent gene
(`abundance=0.0`), a better first rung anyway: every object is then structurally pure gDNA and any
deviation from `f_g = 1` is a false positive with nothing to cancel against it.

### Capture-ON needs care, and the harness says when it is starved

Under capture the donor's off-target density is far below its capture-OFF twin and capture depletes
intergenic space; the sweep prints a STARVED banner naming which object has bitten:

| starved object | does `--genome-length` help? | why |
|---|---|---|
| intergenic region | barely | you would need megabases, which is why the donor injects the background; under capture the toy's own intergenic regions are decoration |
| a boundary, capture off | no | a 0-bp boundary's counts are `density × mean_FL`, independent of chromosome length; the only lever is library depth |
| a boundary, capture on | yes — lengthen it | the gDNA budget is `rate × genome_length` while the probe footprint is fixed, so a longer chromosome hands capture a bigger budget to concentrate onto the same probes |
| what else lifts a boundary? | probe geometry, then binding strength | probes must tile per exon (below); raising `binding_per_base` un-matches the toy from the donor's chemistry, lengthening keeps every harvested global intact |

The mechanism is not panel-specific: under capture the gDNA signal leaves the intergenic and intronic
regions and arrives at the boundaries abutting an exon, so the capture-ON rung is the one the
`intron|exon` boundary actually matters in, and the intron region is dead at every chromosome length —
the well-counted side inverts with capture (`TRAPS: capture-inverts-the-counted-side`). Lengthening
earns a solve on the higher-gDNA capture-ON donors and never on `g00`; say which rows carried a
capture-ON aggregate, and re-measure counts before quoting one. ⚠ Two length slips in
`_donor_sim_params`: `frag_mean` is read from `truth_summary.json`'s `all` row (the mixture), and that
post-truncation realised mean is fed back as the pre-truncation generating mean, so the toy runs a few
bp longer than its donor.

**Probes tile per exon.** Probes are written in transcript space, so a probe spanning an internal sj has
a two-block genomic footprint, and `sim/capture/sampler._split_scale` multiplies every gDNA fragment
overlapping it by `gdna_split_penalty`; tiling across the whole transcript suppressed exactly the
population that spans an `intron|exon` boundary. `_toy_probes` tiles within each exon, so every probe
is unsplit and ends on the boundary. The split-probe case is real in a real panel and deserves its own rung.

### What a toy cannot judge

Magnitudes do not transfer between donors (direction is preserved, size is not — gate on direction and
ordering, quote magnitudes with their donor); it cannot rank defects (five objects cannot say what
fraction of a library's error a mechanism owns — localise on the panel, isolate on a toy, re-measure on
the panel); and it is a correctness instrument, never a profiling one
(`TRAPS: toys-rank-hotspots-backwards`). Gates: `tests/calibration/test_toy_harness.py`, each with its
own perturbation; the donor is a scenario the gates build, so none silently skips without the panel.

---

## 1. The simulated panel (shared backbone)

A real human backbone, not a generated mini-genome: `chr21` + `chr22` + the 92 ERCC spike-in references,
carved from the same GRCh38 + GENCODE v46 sources the production index is built from. Neither is filler:
a single-reference synthetic index once hid a reference-id-space mismatch that silently dropped nearly
every fragment inside `deposit()` while every golden test passed
(`TRAPS: one-reference-hides-refid-bugs`); the spike-ins make the id space non-trivial for RNA, two
genomic chromosomes make it non-trivial on the gDNA path too.

Every fragment-length parameter is measured — the count-weighted mean and sd of a real cfRNA library's
own RNA pool, 206.1 ± 98.3, handed to both origins (`DESIGN.md` §0b). Those are pre-capture parameters
and the truth files are post-capture: hybrid capture selects for length, so the simulator draws the
length marginal as `f_post(w) ∝ f_pre(w) · total_eff(w)`. ⛔ Score against `truth_fragment_lengths.tsv`,
never against `frag_mean` (`TRAPS: capture-selects-for-length`).

### The simulator's transcriptome is a Rigel index — including the nascent entities

`rigel sim` builds (or is given, config key `index:`) a rigel index and simulates its transcript list, so
what is simulated is exactly what `rigel quant` reads: annotated transcripts plus the synthetic nascent
entities `index.create_nrna_transcripts` makes (one single-exon transcript over each multi-exon span,
TSS/TES clustered within `NRNA_MERGE_TOLERANCE`, an annotated single-exon transcript adopted where one
covers the span). `transcript_filter` is refused — filter the GTF before building the index.

| | |
|---|---|
| nascent RNA is a transcript, not a parallel space | its molecules are `entity.nrna_abundance`, sampled on its own template, its reads keep the `nrna_` origin tag. Under `sparse` the level is drawn per entity — off with probability `1 − on_fraction`, else log-uniform over `abundance_ranges`, independent of the mature level; under `additive_ratio` and `fragment_share` it is `Σ abundance × nrna_ratio` over contributors |
| one multinomial over every RNA row | mature and entity rows together, `prob ∝ abundance × capture-aware effective length`; the mature/nascent split follows from molecules and lengths, so a nascent-on condition and its nascent-off twin do not share a bit-identical mature stream |
| capture binds by genomic overlap | every probe → genomic blocks → gDNA, projected onto every transcript whose exons it touches, either strand (ds-cDNA at capture); intron-split pieces take `gdna_split_penalty`. gDNA and nascent are enriched at the same rate under one probe, which is the physics |

Two consequences ruled correct (owner, 2026-08-19): capture depletes nascent RNA about 8× (roughly 20 %
of RNA fragments off capture, under 3 % on — read a condition's `truth_summary.json`), because probes
tile exons and a pre-mRNA is mostly intron, and adding nascent after the capture simulation is wrong —
molecules exist first; and the post-capture fl distribution and abundances are the ground truth
(`truth_kind: post_capture_empirical`), which is why gate G-S5 (§3) is stated per RNA population.

### What the simulator still does not do

| | consequence |
|---|---|
| counts are Poisson by construction | nothing dispersion-dependent validates here |
| the panel is all R1-antisense | the engine can emit either (`ReadSimConfig.r1_sense`, gated in `test_strand_sense_convention.py`) but `orchestrator.run_condition_grid` does not expose it. Real cfRNA is dUTP, so this is not urgent |
| the tool's gDNA reach assumption (`taper_g = 1`) is untested | latent: gDNA is not simulated on the spike-ins |
| each population is written as one contiguous block of read names | a per-fragment truth join checked by "does an impossible label appear?" is nearly blind. ⛔ Gate such a join on a count identity against the scanner's own `stats.total` / `stats.n_read_names` (`_oracle.check_walk_alignment`); `tests/calibration/test_prior_vs_oracle.py` pins both halves |

---

## 2. Building it — `panel.py`, one command per stage

`scripts/sim/panel.py` is the whole loop: every path derived from one panel YAML, every stage resumable,
each stage shelling out to the instrument that owns it. It adds no measurement code.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
CFG=scripts/sim/configs/gdna_ladder.yaml       # the panel the tool is ranked on

python scripts/sim/panel.py status   --config $CFG              # always first
python scripts/sim/panel.py build    --config $CFG              # index + capture probes
python scripts/sim/panel.py simulate --config $CFG --jobs 8     # ~21 min, ~16 GB
python scripts/sim/panel.py cache    --config $CFG --jobs 8     # BOTH caches — mandatory, §0
python scripts/sim/panel.py score    --config $CFG --jobs 8 --arms base oracle noop
python scripts/sim/panel.py report   --config $CFG --arms base oracle
```

`status` prints what exists, what is missing, and which stage to run next. `build` cannot carve the
reference — that needs the source genome/GTF, which the config does not name — so it prints the exact
command and stops:

```bash
python scripts/sim/build_suite_reference.py \
    --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf \
    --refs chr21 chr22 --ercc -o $SUITE/reference
```

The gates are not stages of `panel.py`; run both before quoting anything:

```bash
python scripts/design/simulator_gates.py --suite $SUITE/ladder --reference $SUITE/reference
```

### Two workflow gotchas

* ⛔ **`--force` does not reach the `simulate` stage.** `cmd_cache` passes it to `build_scan_cache.py`
  and `cmd_build` honours it for the index and probes, but `cmd_simulate` shells out to
  `simulate_reads.py` with no such flag, and the simulator skips a condition whose oracle BAM already
  exists (`skip_existing`). So editing a config and re-running `panel.py simulate --force` reports
  success and reproduces the old reads. Delete the condition directories first.
* ⛔ **The zero-gDNA rows are held out of the oracle sweep** (`pass0_vs_oracle.py` scores no row whose
  truth is exactly zero), and `status` counts an oracle condition complete only when `gdna`, `mrna`,
  `nrna` and `_main` are all present. `panel.py cache` fills the `g00` rows with the per-condition
  prewarm and copies `_main` from the scan cache; a panel cached only by running `pass0_vs_oracle.py`
  directly reads ✘ on those rows while being complete for every scorer. One row by hand:

  ```bash
  python scripts/design/pass0_vs_oracle.py --suite $SUITE/ladder --index $SUITE/rigel_index \
      --oracle-cache $SUITE/ladder/oracle_cache --_prewarm gdna_g00_ss_0.50_nrna_mid_capture_off
  ```

  `--_prewarm` is the worker half of `--jobs`, hidden from `--help`, and reuses the shipped loader so a
  stale cache is still refused; do not also pass `--conditions <that row>`, or the hold-out exits first
  with "no contaminated conditions found". A `_main` beside a `g00` row proves nothing about what has
  been measured there (`TRAPS: shard-an-arm-sweep-by-condition`), and `status`'s four-part count is not
  a scorer's requirement — `calibration_oracle.py` needs the five partitions and refuses without them.

* Every panel config states `gdna.genomic_refs: [chr21, chr22]` explicitly; the engine does not infer
  which references carry genomic DNA (`TRAPS: annotated-is-not-genomic`).
* Export `RIGEL_SCRATCH` before a sweep, or the instruments write their per-condition work under `/tmp`
  (tens of GB after a ladder rebuild).

---

## 3. The simulator's own gates — G-S1…G-S6

`scripts/design/simulator_gates.py`, scored on the panel's per-fragment truth (the oracle BAM's read
names). Every gate is directional or an absolute count; none carries a threshold — a pass mark on "how
much longer is a captured fragment" would be inventing the capture efficiency curve.

| | gate | form |
|---|---|---|
| G-S1 | gDNA fragments on an RNA-only reference | absolute count, must be 0 |
| G-S2 | genomic references carrying gDNA | ≥ 2, each non-zero, on every gDNA condition |
| G-S3 | gDNA mean length, capture off → on | strictly greater under capture |
| G-S4 | on-target vs off-target gDNA mean length | on-target strictly longer. A regression guard, not a falsification: it passed with the capture defect present, because the conditional was right and only the marginal was discarded (`TRAPS: a-gate-that-already-passed`) |
| G-S5 | mean length of each RNA population (mature, nascent) separately, capture off → on | strictly greater under capture, per condition × pool — per population, so it is invariant to how much nascent there is |
| G-S6 | gDNA fragments longer than their own reference | 0 |

"On-target" means overlaps a probe, not "its start lands in an exon": the start-territory version is
geometry-confounded and stays inverted under any correct capture model
(`TRAPS: on-target-by-start-is-geometry`); the script prints it underneath as the diagnostic it is.

---

## 4. How results are evaluated

The 0.8.0 metric is the calibration result scored against oracle calibration, not the end-to-end
transcript number (`DESIGN.md` §0b; the judging rules are `SUCCESS.md`). Three questions, three
instruments:

| question | instrument |
|---|---|
| how wrong is calibration, against oracle calibration? — the 0.8.0 metric | `calibration_vs_oracle.py` · `prior_vs_oracle.py` (the `LocusPriors` the EM reads) · `solvability_audit.py` (pass-0) · `pass0_vs_oracle.py` for the T/C/P decomposition |
| how wrong is the end-to-end answer? — per transcript, per pool, against per-fragment truth | `panel.py score` / `report`, i.e. `quant_accuracy.py` — the only end-to-end scorer, a thermometer rather than the target |

⛔ A ceiling only prices what its arm can reach: the effective-length shrinkage is built before
`assemble_priors`, which every measurement arm patches (`SUCCESS.md`, the ruler). Say which call your
arm patches and check it sits downstream of everything you mean to price.

Hard per-fragment label recovery is the wrong target: an unspliced RNA fragment and a gDNA fragment from
the same locus can be sequence-identical. A net-flow reduction, `net(a→b) = flow[a][b] − flow[b][a]`, cancels
the unrecoverable part and keeps only systematic bias; the module that computed it per locus was retired
2026-09-13 for want of an entry point (in git), and the absolute per-transcript error stays the report. Hard-label
metrics are nearly blind to a
calibration-prior change (`TRAPS: hard-labels-miss-soft-change`); the soft 3-pool surplus is not built.

---

## 5. The test suite

```bash
python -m pytest tests/ -q
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
```

`CLAUDE.md` is the home of the standing pass count and accounts every delta: re-derive it, never adjust
it, because several gates are parametrised over the files on disk. Any failure is a regression.
`tests/scenarios_aligned/test_multimap_counting.py::TestParalogMultimapping::test_gdna_sweep[gdna_100]`
is a real EM unidentifiability (`TRAPS: identical-paralogs-are-bimodal`); if it fails again, do not fix
it by moving a seed.

The goldens run under the default sampling mode (`EMConfig.seed = None`, `assignment_mode = "sample"`),
and two runs of the identical pipeline on the identical BAM return different transcript counts
(`TRAPS: the-deliverable-is-not-reproducible-by-default`). Regenerate the goldens twice and diff, and pin
`EMConfig.seed` in any instrument that compares two end-to-end runs (`quant_accuracy.py` does, and prints
a `base_reseed` noise floor beside the effect).

---

## 6. Development discipline for test substrates

Develop on controlled toys, validate on real data — both, in that order: a big suite has confounds that
hide mechanisms, and a toy ranks hotspots backwards (`TRAPS: toys-rank-hotspots-backwards`).

**Profiling is judged on a deep, real, high-complexity RNA-seq library** (owner, 2026-08-17) — never a
panel condition and never a toy. The cfRNA libraries on disk are sparse and small, so they under-fill
exactly the structures whose cost the performance work is about (the grid solve, the per-slot arrays, the
fragment buffer); they are smoke tests that find a defect cheaply, and the deep library decides. This
leaves `TRAPS: real-data-is-a-test-input` and the accuracy panel (the ladder) untouched.

**Read a timing only from back-to-back A/B pairs.** Stages nobody touched drift 25–40 % between runs
taken at different times — page cache, concurrent work — so a before-and-after taken hours apart credits
the machine's mood to the change. Alternate the two arms in one sitting (stash `src/` only) and check
that the untouched stages read 1.00. **Every speed-up is proven a numeric no-op** before it is believed:
`design/rename_identity.py --bam` end to end on a real library and on both ladder references, and
`profiling/sweep_replay.py` for anything inside the sweep.

The instruments are `scripts/profiling/profiler.py` (the whole pipeline as a tree of named stages, with
per-stage peak and held RSS; `--set` for any config field, `--scan-only` for the scan alone, `--compare`
for two reports) and `scripts/profiling/sweep_replay.py` (one calibration sweep, replayed and compared
bit for bit; `--block-slots N|none` replays it at another locus-block size, which must move nothing —
the chunk-exactness of the whole sweep on real data); set `OMP_NUM_THREADS` deliberately. The frozen
references are `~/Downloads/rigel_runs/arms/review_identity_*.json` (two ladder conditions and the LBX0190
library, frozen on the tree of 2026-09-14 that carries the lanes worklist and the landscape's location floor,
and bit-identical through the cleanup that followed) and the captured sweeps
`~/Downloads/rigel_runs/perf/sweeps_MO_3021_step6` (the lanes worklist's tree; `DESIGN.md` §6b.15), with the
deep library's timing baseline beside them in `perf/baseline_2026-09-14/` (two back-to-back pairs at
8 threads). The
refit sweeps' captures carry their message cache, so their replay exercises the cache path and ψ; sweep 0
exercises the whole message layer. `ISSUES: performance-memory-bounded-solve` carries the work.

A both-strand stress test needs ample single-stranded regions (the population prior trains on them).
How to A/B honestly: in-process, opposite extremes, never on a saturated condition, one thing varied,
both arms sharing their random input (`TRAPS: perturb-every-gate`) — a byte-identical hard-label result
is no evidence.
