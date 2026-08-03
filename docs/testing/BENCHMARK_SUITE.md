# The benchmark suite

What it is, how to build it, how to run it, and — the part that matters — **what it can and cannot
judge**. Built 2026-07-30 to replace the 32-condition `ambig_dense_10mb` suite, which was deleted
because it could not judge what it was being used to judge.

> ⛔ **Before quoting any number from this suite, run `scripts/design/suite_resolves.py`.** That is not
> a formality. The suite it replaces was used for months to evaluate a partition change it was
> *structurally incapable of seeing* (`docs/SESSION_HANDOFF.md` §3 trap 15).

---

## 1. Why the old one was replaced

| | `ambig_dense_10mb` (deleted) | why it mattered |
|---|---|---|
| fine nodes vs merged regions | **1,698 == 1,698** | could not see a partition change **at all** |
| fragment-length variance | **exactly 0** — every fragment 200 bp | nothing length-dependent was exercised, and length is the core of the design |
| count overdispersion | **0 by construction** | nothing dispersion-dependent validated |
| terminus+splice-site seams | 11.8 % vs the real 1.0 % | **12x** over-represented |
| performance hotspots | ranked **backwards** vs real cfRNA | a whole analysis was spent on the toy's #1 hotspot |

---

## 2. What it is

**A real human backbone, not a generated mini-genome.** `chr22` plus the **92 ERCC** spike-in
references, carved out of the same GRCh38 + GENCODE v46 sources the production index is built from.

⚠ **The ERCC references are not filler.** `docs/SESSION_HANDOFF.md` §3 trap 20: a single-reference synthetic
index once hid a reference-id-space mismatch that silently dropped **476,719 of 476,732** real fragments
inside `deposit()` while every golden test passed. 92 extra references cost ~83 kb and make that space
non-trivial.

⭐ **It is structurally human at 2 % of the size** — which is the whole argument for a real chromosome:

| | human index | **chr22 + ERCC** | ⛔ deleted suite |
|---|---|---|---|
| nodes | 1,043,881 | **22,231** | 1,698 |
| median node | 151 bp | **147 bp** | — |
| **nodes / merged regions** | 1.387 x | **1.386 x** | **1.000 x** |
| **termini the merge hid** | 53.42 % | **52.72 %** | ~0 |
| seams: terminus / splice / BOTH | 40.70 / 58.31 / 0.99 % | **40.82 / 57.96 / 1.22 %** | 11.8 % BOTH |
| 1 bp nodes | 15,687 | **332** | — |

### The pilot grid

The owner's 8-condition minimum, a clean 2³: **gDNA {none, 100 %} × strand {0.50, 0.99} × capture
{off, on}**, with nascent held at `none`.

⚠ **`5,000,000` is the RNA depth.** gDNA is *added on top* at `rate × n_rna`, so a `gdna100` condition
is 5 M RNA **+ 5 M gDNA = 10 M fragments**; a `gdna_none` condition is 5 M.

⚠ **Every fragment-length parameter is measured, not chosen** — the count-weighted mean and sd of
MO_3021's own pools (`rigel/fragment_lengths.feather`, the library that is also the S3 byte-identity
gate): **RNA 206.1 ± 98.3**, **gDNA 156.5 ± 124.6**. The deleted suite used `frag_std: 0`.

| | |
|---|---|
| build time | **741 s**, 8 workers |
| peak RSS | **4.09 GB** |
| on disk | 8.3 GB |
| scan cache | **45.7 s** once, **0.16 s** to reload — 291 x |

---

## 3. How to build it

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite
REFS=~/Downloads/rigel_runs/refs

# 1. carve the backbone out of the real genome + annotation
python scripts/sim/build_suite_reference.py \
    --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf \
    --refs chr22 --ercc -o $SUITE/reference

# 2. index it  (1.4 s)
rigel index --fasta $SUITE/reference/genome.fa --gtf $SUITE/reference/genes.gtf \
    --collapse-duplicate-transcripts --no-mappability --no-tsv -o $SUITE/rigel_index

# 3. design the capture panel — this IS requirement (a), the density step
python scripts/sim/design_suite_probes.py --gtf $SUITE/reference/genes.gtf \
    -o $SUITE/reference/capture_panel.tsv --capture-fraction 0.5

# 4. simulate the 8 conditions  (741 s)
python scripts/sim/simulate_reads.py --config scripts/sim/configs/chr22_pilot.yaml -j 8

# 5. cache the scans so calibration can be iterated without rescanning  (45.7 s)
python scripts/design/build_scan_cache.py --index $SUITE/rigel_index --suite $SUITE/pilot
```

⭐ The **capture panel is not scaffolding** — it covers 720 of 1,440 chr22 gene groups, so the
captured/uncaptured boundary is a sharp density step with real transcripts on both sides. Hybrid capture
is ~1000× on exons (`docs/SESSION_HANDOFF.md` §1 fact 15), so that boundary is a 3-decade cliff in gDNA density.
A uniform background cannot distinguish "the relay works" from "the global prior reached it".

---

## 4. ⛔ The gate: prove it resolves the axis BEFORE quoting it

```bash
python scripts/design/suite_resolves.py $SUITE/rigel_index \
    --compare-index $REFS/rigel_index --suite $SUITE/pilot
```

⭐ **There are no tuned thresholds in it, deliberately** — a pass mark would be a magic number. Every
requirement is scored against its **degenerate value**: the number a structurally blind suite scores.
Those are boundaries, not choices — an unresolvable partition is *exactly* 1.000×, a suite with no
length variation has variance *exactly* 0, a Poisson simulator has ω *exactly* 0. A requirement passes
iff it lands strictly on the non-degenerate side. The human index prints alongside for scale, never as
a threshold.

### Current verdict — 6 of 8 resolve

| | requirement | measured | blind |
|---|---|---|---|
| (g) | nodes / merged regions | **1.386 x** | 1.000 |
| (g) | termini the merge hid | **4,906** (52.72 %) | 0 |
| (d) | termini strictly inside an exon | **4,958** (53.28 %) | 0 |
| (e) | single-stranded nodes | **18,127** (81.5 %) | 0 |
| (e) | both-stranded nodes | **3,313** | 0 |
| (e) | single-stranded per both-stranded | **5.47** | 0 |
| (a) | capture density step, p90/p10 | **2,123 x** | 1.0 |
| (b) | fragment-length sd, worst pool | **87.4 bp** | 0 |
| **(c)** | replicate pairs (needed for ω) | ⛔ **0** | 0 |
| **(f)** | conditions at 0 < gDNA ≤ 10 % with capture | ⛔ **0** | 0 |

**Neither failure is a defect in the suite; both are named work** — see `TODO.md`. (c) needs the
overdispersion mechanism *and* replicate conditions. (f) needs one more gDNA rate in the band real
libraries live in; the pilot grid is `none`/`100 %` by the owner's 8-condition spec.

⭐ **The gate's teeth are proven on three degenerate inputs**, each failing for its own reason: a
reference shaped like the deleted suite (121 nodes == 121 merged regions), a "starved toy" with
both-strand nodes and zero single-stranded ones, and truth files written to the deleted suite's exact
shape (sd 0, capture off, no replicates).

---

## 5. How results are evaluated

⭐ **Net fragment flow is the primary metric, and it is already implemented** — `rigel.sim.analysis`
(`analyze_net_flow`, `FlowData`), reachable via `scripts/sim/evaluate_suite.py`. The methodology and the
reasoning behind it are in `docs/testing/BENCHMARKING.md`; the short version:

* **Hard per-fragment label recovery is the wrong target.** An unspliced RNA fragment and a gDNA
  fragment from the same locus can be sequence-identical and genuinely unrecoverable. Build the per-locus
  `flow[true][assigned]` matrix and reduce to `net(a→b) = flow[a][b] − flow[b][a]`: symmetric,
  unrecoverable misassignment cancels and **only systematic bias survives**.
* ⚠ **The hard-label metric is nearly blind to a calibration-prior change** — a real change can move the
  soft 3-pool counts by tens of thousands of fragments while the hard-label net is byte-identical. Treat
  a byte-identical hard-label result as **no evidence**, not as no change.
* ⚠ **Report absolute per-transcript error alongside the net**, because net cancels.
* ⚠ **`gdna_none` conditions are pure false-positive tests** — any net gDNA→RNA there is false by
  construction, and unstranded (ss = 0.50) is the hard case because there is no strand clue at all.

⛔ **The soft 3-pool surplus does not exist yet.** `docs/testing/BENCHMARKING.md` names it as the *primary* pool
metric and the surviving code computes only the hard-label version. Logged in `TODO.md`.

---

## 6. ⛔ What this suite cannot do today

**It cannot produce a calibration number**, because `calibrate()` does not run. ⚠ The reason has
changed and the old one is stale: `fl.py` and `substrate.py` were rewired in S5.b and S5.d. What remains
is **S5.e** (`build_node_geometry` and `bp_solver`'s per-face consumers) and **S5.f** (`calibrate` itself,
`CalibrationResult`, `priors`, `pipeline`). See `docs/calibration/S5_DESIGN_LOG.md` §2.

⚠ **Every scan cache built before 2026-07-30 is now correctly REFUSED** — S5.a added `length_sum` to the
payload and S5.a-3 added a `payload_schema_digest` to the cache key, so a stale cache fails at the door
instead of deep in the loader. Re-scan: 45.7 s for all 8 conditions.

Everything else above — the backbone, the index, the pilot, the resolution census, the scan cache
machinery — is built and verified without calibration.
