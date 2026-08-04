# TESTING — the substrate, the gates, and what they can judge

**Two rules govern this whole file**, and both came from real failures (`TRAPS.md` A7, A8):

> ⛔ **Prove the substrate before you prove the code.** When a simulated axis is the axis you are judging,
> gate the simulator on it.
>
> ⛔ **Prove the suite can resolve the axis you are changing, before quoting a number from it.**

---

## 0. ⭐⭐ THERE ARE TWO PANELS AND A TOY HARNESS, AND THEY ANSWER DIFFERENT QUESTIONS

Both panels live under `~/Downloads/rigel_runs/suite/` and share one reference and one index. ⛔ **Neither
replaces the other — deleting either loses a question.** ⭐⭐ **And §0b's TOY HARNESS is the third
substrate: the panels say how much error there is and where, a toy says WHY.** Localise on a panel,
isolate on a toy, re-measure on the panel.

| | `pilot/` — 8 conditions | ⭐ `ladder/` — 36 conditions |
|---|---|---|
| **judges** | **Stage A**, the accumulator | **Stage B**, the calibration solver |
| fragment lengths | REALISTIC and DIFFERENT per origin (RNA 206 ± 98, gDNA 157 ± 125) | ⭐ **IDENTICAL for both origins** (206 ± 98, [50, 500], read 100) |
| why that way | a length model can only be judged where the components' lengths actually differ | with the length channel neutralised, composition can only come from **density and strand**, so residual error is attributable to those |
| gDNA axis | {0, 100 %} rate ⇒ f_gdna ∈ {0, 0.5} | ⭐ **9 rungs, f_gdna 0 → 0.98** |
| depth | 10 M RNA, gDNA **added on top** (so 10–20 M total) | ⭐ **10 M TOTAL, fixed**; the rate decides only the SPLIT |
| config | `scripts/sim/configs/pilot.yaml` | `scripts/sim/configs/gdna_ladder.yaml` |

### The gDNA ladder, in detail

`f_gdna = rate/(1+rate)`, and the RNA share thins as gDNA rises because the **total** is held:

| rung | g00 | g01 | g05 | g10 | g25 | g50 | g75 | g90 | g98 |
|---|---|---|---|---|---|---|---|---|---|
| f_gdna | 0 | 0.01 | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.98 |
| n_rna | 10.0 M | 9.90 M | 9.50 M | 9.00 M | 7.50 M | 5.00 M | 2.50 M | 1.00 M | 0.20 M |

× strand {0.50, 0.99} × capture {off, on} = **36**. Owner ruling: real libraries run from almost zero
gDNA to **> 98 %** of all fragments and Rigel must be robust across all of it; 10 M total per condition
is the ceiling, and the RNA-side abundance accuracy that thins at the top rungs is an accepted
trade-off — it is a property of such libraries, not an artifact.

⚠ **The focus condition is `ss_0.50 + capture_on`** (owner). Everything else is a CONTROL that must not
regress, and the `g00` rung is the false-positive check. ⛔ Do not tune on the controls.

⛔ **"EQUAL" IS CONFIGURED, NOT ACHIEVED — AND THE RESIDUAL IS MEASURED, NOT ASSUMED.** Identical
parameters still leave a realised gap of **+4.68 bp off capture / +3.57 on** (`TRAPS.md` F11: a mature
fragment must fit inside its transcript, gDNA need not, so transcript truncation pulls RNA down). What
that residual is WORTH was measured before it was accepted, one thing varied: replacing both pmfs by
their pooled average moves the per-object error by **−0.0002** off capture and **−0.0054** under it —
under 2.5 % of the error, with a *fitted* gap 3–5× larger than the residual. Owner ruling: carry it.

⭐ **An ORACLE CACHE ships beside the panel** at `ladder/oracle_cache/` (36 conditions, 157 MB). It holds
the per-origin split payloads, so `--oracle-cache` turns a 4-minute-per-condition oracle build into
seconds. ⭐⭐ **It stays valid across every CALIBRATION change** — the oracle depends only on the
accumulator and the index — so one cache serves an entire solver-debugging campaign. It is keyed by the
scan cache's own key, so a stale one is refused rather than silently used. ⛔ Rebuild it after any
**accumulator** change (delete the directory; it repopulates).

---

## 0b. ⭐⭐ THE TOY HARNESS — a mini chromosome you define, calibrated in under a second

`scripts/design/toy_harness.py`. ⭐ **This is the third substrate and the one to reach for FIRST when a
mechanism needs isolating.** The two panels above tell you *how much* error there is and *where*; a toy
tells you *why*, because it is small enough to read every object and cheap enough to sweep one variable
seven times.

⛔ **What it replaces.** Every defect found before it was found by sifting a 36-condition, 10-million
fragment panel for badly-behaved objects and reasoning backwards from a node id. That is slow, the
example is never quite the one you wanted, and a fix can only be measured in aggregate where two errors
cancel. The harness root-caused C1 in **five objects and 0.1 s** after weeks of panel-level work.

### The one idea that makes it work

A toy cannot fit the library-level quantities calibration needs — one transcript has no population to
estimate a strand balance, an enrichment landscape or an intergenic background from. So **a real cached
condition acts as DONOR**: it is calibrated once, and its fitted bundle is injected
(`InjectedCalibrationPriors`, whose docstring already specified this use). The toy supplies only the
controlled per-node **geometry**, which is the thing under study.

| the donor supplies | so the toy never invents |
|---|---|
| κ, both strand overdispersions, the Fisher noise-floor sample sizes | the strand deadband behaves exactly as it does on real data |
| the enrichment NPMLE, the intron background, ρ_bg | a handful of nodes cannot fit these |
| both fragment-length pmfs | passed as `calibrate` kwargs, not part of the priors bundle |
| capture on/off + its numeric knobs | reproduced in the toy's own **simulation**, with probes written from the spec |
| frag mean/sd/min/max, read length, strand specificity | read from the donor's post-capture truth |
| ⭐ **gDNA density per base** | see below — this one is not a knob |

⭐⭐ **THE gDNA LEVEL IS DERIVED, NOT CHOSEN, AND THAT IS NOT NEGOTIABLE.** The injected enrichment
landscape is an **absolute** log-density model, so a toy at the wrong depth is a *different* library, not
a small one. The harness measures the donor's gDNA counts-per-base on the donor's own
structurally-pure-gDNA nodes (`Σcount / ΣE`, `EQUATIONS.md` §7.1) and simulates the toy to match it.
`ToySpec` therefore has **no gDNA field**, and a gate asserts it never grows one.

### ⭐ Terminology — one word per concept (owner, 2026-08-04)

| | |
|---|---|
| **counts** | discrete **integer** fragment counts. What the accumulator stores, and the solver's Poisson `n` |
| **density** = **abundance** | **counts per base.** The two words mean the same thing and are used interchangeably |

⛔ Not the simulator's molar `abundance=` field, which is a per-transcript weight, not a density.
⚠ And a node's stored counts are **contained** counts, so `counts = density × effective_length`, never
`density × bp`. The sweep reports the density it *asked* for and the counts each object *received*, side
by side, and never converts one into the other.

### Running it

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite

python scripts/design/toy_harness.py --list            # the spec ladder, simplest first

# one spec against one donor — prints EVERY object beside per-object truth   (~0.1-5 s)
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g25_ss_0.50_nrna_none_capture_off

# ⭐ sweep the transcript's RNA density; the gDNA background stays PINNED, so one variable moves
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g25_ss_0.50_nrna_none_capture_off --sweep-density

python scripts/design/toy_harness.py --spec all --donor <cond>     # the whole ladder
```

⚠ Any of the 36 ladder conditions can be the donor, and **which one you pick is an experimental
variable** — it sets capture on/off, the strand regime and the gDNA level. Harvesting costs one scan +
one calibrate (~30 s) and is deliberately **not cached**: the bundle is a function of the calibration
code that fit it, so a stored copy would go stale on exactly the changes the harness exists to test.
⭐ Harvest once per session and run many toys against it (`harvest()` then `run_toy()` in a loop).

### Writing a new spec

Add a `ToySpec` to `SPECS` in the harness. The ladder is ordered simplest-first and **each rung adds
exactly ONE structure to the one before it** — that is the point: when a row goes wrong, the thing that
changed is the thing to look at.

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

⚠ **There is no gene-free rung**: `TranscriptIndex` requires at least one transcript, so an unannotated
chromosome cannot be indexed. Use a **silent** gene (`abundance=0.0`) instead — which is a better first
rung anyway, since every object is then structurally pure gDNA and any deviation from `f_g = 1` is a
false positive with nothing to cancel against it.

### ⚠⚠ CAPTURE-ON needs care, and the harness will tell you when it is starved

Under capture the donor's **off-target** density is ~24× lower than its capture-OFF twin
(0.00106 vs 0.02566 counts/bp on `g25`), and capture actively *depletes* intergenic space. Three
consequences, all measured, and the sweep prints a ⛔ STARVED banner naming which one has bitten:

| starved object | does `--genome-length` help? | why |
|---|---|---|
| intergenic **node** | ⚠ barely | counts scale with bp at fixed density, but capture depletes off-target so hard that 57 kb yields **2 counts**. You would need megabases — which is exactly why the real panel is 93 Mb of chr21+chr22, and exactly why the donor **injects** `background` / `intron_background` / ρ_bg. ⭐ Under capture the toy's own intergenic nodes are decoration; the background comes from the donor by design |
| intergenic\|exon **edge** | ⛔ **no, not at all** | a 0-bp line's counts are `density × mean_FL`, **independent of the chromosome length**. At off-target density that is ~0.26 counts however long you make it |
| … so what lifts an edge? | ⭐ **a PROBE reaching the transcript end** | which is why toy probes are **TILED** across each transcript rather than centred. A single central probe leaves a 2 kb transcript's ends at off-target density and the whole chain starves. This is the owner's C6 geometry: a probe on a first or last exon partially overlaps the boundary and enriches it |

⭐ With tiled probes the capture-ON arm works and shows the **same** defect as capture-OFF (exon truth
0.049, `fg_loc` 0.510, pred 0.771) — which is what makes "one rule for both arms" testable at all.

### ⚠ What a toy CANNOT judge

* **Magnitudes do not transfer between donors.** The same `nrna` contrast is a factor of **31** against
  the `g25` ladder donor and **1.25** against a six-gene synthetic donor: direction preserved, size not.
  A small donor cannot determine the enrichment landscape or the intron background, so its injected
  globals are a different regime. ⛔ Gate on direction and ordering; quote magnitudes with their donor.
* **It cannot rank defects.** Five objects cannot tell you what fraction of a real library's error a
  mechanism owns. Localise on the panel, isolate on a toy, then re-measure on the panel.
* **It shares `TRAPS.md` A9's warning** — a toy ranks performance hotspots backwards. It is a
  *correctness* instrument, never a profiling one.

Gates: `tests/calibration/test_toy_harness.py` — 7, each carrying its own perturbation, and the donor is
a scenario the gates **build** so none of them silently skips when the panel is absent.

---

## 1. The simulated panel (shared backbone)

**A real human backbone, not a generated mini-genome:** `chr21` + `chr22` + the **92 ERCC** spike-in
references, carved out of the same GRCh38 + GENCODE v46 sources the production index is built from.

⚠ **Neither the ERCC references nor chr21 is filler.** A single-reference synthetic index once hid a
reference-id-space mismatch that silently dropped 476,719 of 476,732 fragments inside `deposit()` while
every golden test passed (`TRAPS.md` E1). 92 spike-ins make the reference-id space non-trivial **for
RNA** — but they are RNA-only, and gDNA takes a *different* branch through the scanner. **Two genomic
chromosomes is what makes the space non-trivial on the gDNA path too.**

### The grid

The owner's 8-condition minimum, a clean 2³: **gDNA {none, 100 %} × strand {0.50, 0.99} × capture
{off, on}**, nascent held at `none` so every `gdna_none` cell is a pure false-positive test.
Config: `scripts/sim/configs/pilot.yaml`.

⚠ **`10,000,000` is the RNA depth**, and gDNA is added **on top** at `rate × n_rna` — so a `gdna100`
condition is 20 M fragments and a `gdna_none` condition is 10 M. The depth is set so that fragments per
annotated base — and therefore per **node** — matches the chr22-only panel it replaced; calibration
accuracy depends on per-node counts, so holding that fixed is what keeps a re-baseline attributable to a
code change.

⚠ **Every fragment-length parameter is measured, not chosen** — the count-weighted mean and sd of a real
cfRNA library's own pools: **RNA 206.1 ± 98.3**, **gDNA 156.5 ± 124.6**.

⭐⭐ **Those are PRE-capture parameters and the truth files are POST-capture.** Hybrid capture selects for
length, so the simulator draws the length marginal proportional to the capture-weighted opportunity,
`f_post(w) ∝ f_pre(w) · total_eff(w)`. ⛔ **Score against `truth_fragment_lengths.tsv`, never against
`frag_mean`** — the configured parameters describe a library that was never sequenced (`TRAPS.md` F5).

### ⚠ What the simulator still does NOT do

| | consequence |
|---|---|
| counts are **Poisson by construction** — a multinomial at fixed abundance, measured ω < 5e-5 | nothing dispersion-dependent validates here. Real junction overdispersion is ≤ 0.02–0.03 |
| only **R1-antisense** libraries | `strand_specificity` is a swap probability about a *fixed* orientation, never a choice of orientation, so no condition exercises the R1-sense (KAPA-style) branch. Real cfRNA is dUTP, so this is not urgent |
| the tool's **gDNA reach** assumption is untested | `node_geometry` says gDNA's template is the chromosome, so `taper_g = 1`. True for 50 Mb, false for a 273 bp contig. Latent: it goes live only when a short reference has two nodes, and gDNA is no longer simulated on the spike-ins at all |

---

## 2. Building it

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite
REFS=~/Downloads/rigel_runs/refs

# 1. carve the backbone out of the real genome + annotation   (~2 s)
python scripts/sim/build_suite_reference.py \
    --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf \
    --refs chr21 chr22 --ercc -o $SUITE/reference

# 2. index it   (~2 s, deterministic)
rigel index --fasta $SUITE/reference/genome.fa --gtf $SUITE/reference/genes.gtf \
    --collapse-duplicate-transcripts --no-mappability --no-tsv -o $SUITE/rigel_index

# 3. design the capture panel — this IS the density step   (~17 s)
python scripts/sim/design_suite_probes.py --gtf $SUITE/reference/genes.gtf \
    -o $SUITE/reference/capture_panel.tsv --capture-fraction 0.5

# 4. simulate the 8 conditions   (~21 min, 8 workers, ~16 GB)
python scripts/sim/simulate_reads.py --config scripts/sim/configs/pilot.yaml -j 8

# 5. cache the scans so calibration can be iterated without rescanning   (~2 min, 198x on reload)
python scripts/design/build_scan_cache.py --index $SUITE/rigel_index --suite $SUITE/pilot

# 6. ⛔ THE GATES — run BOTH before quoting anything
python scripts/design/simulator_gates.py --suite $SUITE/pilot --reference $SUITE/reference
python scripts/design/suite_resolves.py $SUITE/rigel_index --suite $SUITE/pilot
```

⚠ **`pilot.yaml` states `gdna.genomic_refs: [chr21, chr22]` explicitly.** The engine does **not** infer
which references carry genomic DNA, and a config that asks for gDNA without stating it is rejected — "has
an annotation" is not "is genomic" (`TRAPS.md` E2).

⚠ **The scan cache is refused, not silently accepted, when it does not describe the index it is loaded
against** — keyed on `graph_hash`, `reach_digest` and `payload_schema_digest`. Any accumulator change
invalidates every cache by design, so re-run step 5.

---

## 3. The simulator's own gates — G-S1…G-S6

`scripts/design/simulator_gates.py`, scored on the panel's per-fragment truth (the oracle BAM's read
names). ⭐ **Every one is directional or an absolute count. Not one carries a threshold** — a pass mark on
"how much longer is a captured fragment" would be inventing the capture efficiency curve.

| | gate | form |
|---|---|---|
| **G-S1** | gDNA fragments on an RNA-only reference | absolute count, must be **0** |
| **G-S2** | genomic references carrying gDNA | **≥ 2**, each non-zero, on every gDNA condition |
| **G-S3** | gDNA mean length, capture off → on | **strictly greater** under capture |
| **G-S4** | on-target vs off-target gDNA mean length | on-target **strictly longer**. ⚠ A regression guard, not a falsification — it passed *with* the capture defect present, because the conditional was right and only the marginal was discarded (`TRAPS.md` A3) |
| **G-S5** | \|μ_g − μ_r\|, capture off → on | **strictly narrower** under capture |
| **G-S6** | gDNA fragments longer than their own reference | **0** |

⛔ **"On-target" means OVERLAPS A PROBE, not "its start lands in an exon."** The start-territory version is
geometry-confounded and stays inverted under any correct capture model — `TRAPS.md` F6. The script prints
the start-territory table underneath as the diagnostic it is.

---

## 4. Can the suite resolve the axis? — `suite_resolves.py`

⭐ **There are no tuned thresholds in it, deliberately.** Every requirement is scored against its
**degenerate value** — the number a structurally blind suite scores — and passes iff it lands strictly on
the non-degenerate side. Those are boundaries, not choices: an unresolvable partition is *exactly* 1.000×,
a suite with no length variation has variance *exactly* 0, a Poisson simulator has ω *exactly* 0, and a
capture arm that is length-neutral narrows the length gap by *exactly* 0.00.

The requirements: (a) a capture **density step**, (b) fragment-length **variance**, (c) non-Poisson
**counts**, (d) termini strictly **inside an exon**, (e) ample **single-stranded** nodes, (f) a low-gDNA ×
strong-capture **corner**, (g) **partition resolution**, (h) a **narrowed length-gap** regime.

⚠ **Two are known-failing and both are named work**, not suite defects:
* **(c)** needs the overdispersion mechanism built in *and* replicate conditions.
* **(f)** needs one gDNA rate in the 1–10 % band real libraries live in — one line of `pilot.yaml` plus a
  re-run of the affected conditions. It opens exactly the regime where `TRAPS.md` F4 says the hardest
  failure mode lives.

⭐ **The gate's teeth are proven on three degenerate inputs**, each failing for its own reason: a reference
shaped like the deleted 10 Mb suite (121 nodes == 121 merged regions), a "starved toy" with both-strand
nodes and no single-stranded ones, and truth files written to sd 0 / capture off / no replicates.

---

## 5. How results are evaluated

⭐ **Net fragment flow is the primary metric** — `rigel.sim.analysis` (`analyze_net_flow`, `FlowData`), via
`scripts/sim/evaluate_suite.py`.

**Hard per-fragment label recovery is the wrong target.** An unspliced RNA fragment and a gDNA fragment
from the same locus can be sequence-identical and genuinely unrecoverable. Build the per-locus
`flow[true][assigned]` matrix and reduce to `net(a→b) = flow[a][b] − flow[b][a]`: symmetric,
unrecoverable misassignment cancels, and **only systematic bias survives**.

⚠ **Report absolute per-transcript error alongside the net**, because net cancels.
⚠ And hard-label metrics are nearly blind to a calibration-prior change — `TRAPS.md` B4.

⛔ **Missing:** the **soft** 3-pool surplus. The surviving code computes only the hard-label version, and
the soft one is the metric that actually sees a prior change.

---

## 6. The test suite

```bash
python -m pytest tests/ -q                     # ⛔ never bare `pytest` — the repo root must be on sys.path
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
```

⚠ **The standing baseline is 21 `test_golden_output` failures plus one.** The 21 are stale expectations,
expected, and they have moved seven times; ⛔ **regenerate them once, at the end of the accumulator work,
twice, and diff** — they run under the default sampling mode, so a flaky expectation baked in now is
permanent, and regenerating is **not** validating.

The 22nd is `tests/scenarios_aligned/test_paralogs.py::test_gdna_sweep[gdna_100]` — a real EM
unidentifiability, not flakiness (`TRAPS.md` D9). ⛔ Do not fix it by moving a seed.

⭐ **"22 failures, 21 goldens and the paralog row" is the baseline. A 23rd failure, or any other
non-golden name in the list, is a regression.**

---

## 7. Development discipline for test substrates

**Develop on controlled toys, validate on real data.** A big suite has confounds that hide mechanisms; a
toy ranks hotspots backwards (`TRAPS.md` A9). Both, in that order.

**Any both-strand stress test needs ample single-stranded nodes** — the population prior trains on them,
and a "starved toy" is one of the three degenerate inputs `suite_resolves.py` is proven against.

**How to A/B honestly:** in-process, opposite extremes, never on a saturated condition, one thing varied,
and both arms sharing their random input (`TRAPS.md` A2) — a byte-identical hard-label result is **no
evidence**.
