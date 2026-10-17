# Rigel 0.8.0 — how we got here, what it cost, and the way to a release (2026-10-12)

Sandbox, not authoritative. Written for the owner, who asked why we left the tree that was nearly ready to ship, why
the intron factory was removed, why capture was rebuilt, and whether the work since is progress.

All numbers below were measured today under one protocol: scan pinned to one thread, fractional assignment, EM seed 0.
Every tree ran on the same simulated panels and real libraries. Figures are transcript error, Σ|est − truth| / truth
in %; gene error is in parentheses. Receipts and drivers are in `$RIGEL_RUNS/prototypes/2026-10-12_reboot/`.

| Name | Commit | What it is |
|---|---|---|
| **tag** | v0.7.1 `ff26d8e5` | the released product |
| **10-05** | `2c482115` (Oct 5, 22:28) | the last commit before the rewrite, junction clip included |
| **foundation** | `3e7e0fbc` (Oct 10, 07:13) | 10-05 plus the count foundation (factory removed, observed strand claims); same capture reader as 10-05 |
| **main** | `624dfade` (Oct 11) | the foundation plus the per-object capture reader (`7240e929`, Oct 10, 12:50) and the multimapper second pass |
| **candidate** | main + 19 lines | the per-object reader on own evidence (split·DNA, this morning's prototype) |

## 1. The short answer

**The problem that started all this was solved on October 10 at 07:13**, by the count foundation, with the OLD capture
reader:
- RNA-short, stranded, capture-off goes from 42.30 % on the 10-05 tree to **3.54 %** on the foundation.
- The foundation is the best tree we have measured: best or tied on every in-scope stratum of the ladder and both
  fl-gap panels, and on most robustness panels.

**Five hours later the per-object capture reader landed on top of it.** That reader was right in principle: your
spectrum ruling said capture must be local, and it is measurably needed for very small panels. But it degraded almost
everything else:
- ladder stranded capture-on, 3.87 → 6.05 %;
- ladder deferred stratum, 7.63 → 62.96 %;
- capture-off zero-gDNA rows, 1.93 → 13.61 %;
- the strand-noise and weak-strand panels.

**The next day's report credited the RNA-short fix to the new reader**, and blamed the ladder's deferred collapse on
the factory removal. Both attributions were wrong, and they made the rewrite look necessary.

**Since then every session has used main as its baseline**, and tried to repair the per-object reader's noise with
more machinery: averages, local priors, own-evidence splits, a shared field, today's pooling. None of those was ever
compared against the foundation on the robustness panels.

## 2. The machine, in four steps

1. **Scan.** Every fragment is tallied onto the genome's regions and boundaries, by strand, spliced or not.
2. **Calibration.** For each object it estimates how many unspliced fragments are gDNA, sense RNA or antisense RNA.
   It combines strand deconvolution, messages between neighbouring objects, and a landscape prior: the population
   distribution of gDNA density over the objects.
3. **Capture weights.** gDNA lands uniformly before capture, so an object's gDNA density relative to a reference is its
   capture efficiency. The efficiencies contract every EM component's effective length.
4. **EM.** Per locus, transcripts, nascent spans and the gDNA component compete for fragments.

**The OLD reader (10-05 and foundation).**
- Its reference is the landscape's located enriched mode: the density of the captured population.
- An object's efficiency is the posterior mean of its gDNA density, clipped at that reference.
- With no located mode it applies no capture at all.
- It is robust because efficiencies are bounded, smooth, and exactly uniform when the library shows no captured
  population. It is library-level because it needs a captured population large enough to form a mode.

**The NEW reader (main).**
- An object's weight is the posterior mode of its own gDNA density under the landscape, relative to the typical object.
- It is local and unbounded, so the noise in every object's evidence reaches the EM.

## 3. What happened, and why

| Date | Event | Reason |
|---|---|---|
| 10-01 | RNA-short stranded capture-off fails, 42.3 % | 393 exons shorter than a gDNA fragment hold strand-noise "gDNA"; the old reader's census counted them as a captured population |
| 10-02–04 | Narrow repairs tried; each fixed the row and broke another stratum | — |
| 10-05 eve | **Your ruling: capture is a spectrum; no detector** | The old reader's "no mode → no capture" is a library-level switch; it read the two sparse plasma panels as uncaptured |
| 10-05–07 | Per-object reader prototypes; the factory audit | The record twice concluded that the reader was fine given correct per-object counts, and that the counts were the problem |
| 10-09 | Count foundation accepted | Factory removed; observed strand claims; conservative allocation under uncertainty |
| **10-10 07:13** | **Count foundation committed (`3e7e0fbc`)** | **This fixed the RNA-short row.** Its counts no longer give the tiny exons false gDNA, so no false mode forms. |
| 10-10 a.m. | Capture-reader benchmark: "ship 0.8.0 with the shipped reader" | That reader was the foundation's |
| **10-10 12:50** | **Per-object reader committed (`7240e929`)** | Your instruction after rejecting the on/off switch. Its capture-off cost was accepted as "the price of no detector". |
| 10-10/11 | Multimapper second pass (`d87313ab`); reader speed | Independent of the reader; VCaP gDNA fraction 0.240 → 0.253 (truth ≈ 0.252) |
| 10-11 | Ladder re-scored against the 09-30 arms | Misattributed: RNA-short fix credited to the reader; the ladder's deferred collapse blamed on the factory |
| 10-11–12 | Repair attempts on the per-object reader; reassessment; candidate; pooling | Each judged against main, never against the foundation |

## 4. The measurements

### Release panels

| Stratum | Tag | 10-05 | Foundation | Main | Candidate |
|---|---:|---:|---:|---:|---:|
| Ladder stranded OFF, contaminated / g00 | 3.95 / 3.70 | 2.02 / 1.93 | **2.01 / 1.93** | 3.37 / 13.61 | 2.62 / 3.81 |
| Ladder stranded ON, contaminated / g00 | 9.91 / 20.22 | 3.85 / 6.46 | **3.87 / 6.47** | 6.05 / 15.10 | 6.24 / 6.61 |
| Ladder unstranded OFF, contaminated | 10.29 | 2.29 | 2.30 | 2.40 | 2.34 |
| Ladder unstranded ON (deferred) | 33.21 | 10.42 | **7.63** | 62.96 | 52.73 |
| RNA-short stranded OFF / ON | 6.31 / 20.33 | 42.30 / 19.97 | **3.54 / 13.32** | 4.99 / 14.24 | 4.61 / 20.09 |
| RNA-short unstranded OFF / ON (deferred) | 6.00 / 29.76 | 3.78 / 19.41 | 3.80 / **18.47** | 3.95 / 68.06 | 3.94 / 63.65 |
| RNA-long stranded OFF / ON | 5.23 / 7.38 | 2.26 / 4.48 | **2.26 / 4.25** | 3.52 / 6.09 | 3.46 / 6.33 |
| RNA-long unstranded OFF / ON (deferred) | 5.29 / 28.99 | 2.13 / **6.86** | 2.16 / 43.22 | 2.55 / 34.11 | 2.21 / 14.77 |

### Robustness panels (test chromosome, contaminated rows)

| Row | Tag | 10-05 | Foundation | Main | Candidate |
|---|---:|---:|---:|---:|---:|
| test stranded OFF / ON | 8.20 / 8.58 | 7.39 / 7.19 | 7.61 / 7.32 | 9.35 / 8.17 | 8.03 / 8.27 |
| test strand 0.70 ON | 10.28 | **6.91** | 7.76 | 7.88 | 9.87 |
| strand noise (odg05) stranded ON / strand 0.70 ON | 14.52 / 16.61 | 13.44 / 12.90 | **13.15 / 12.87** | 16.99 / 21.62 | 17.36 / 27.59 |
| weak binding (bind01) strand 0.70 ON | 22.79 | 11.65 | **11.37** | 12.20 | 25.12 |
| scarce gDNA (lowg) ON | 8.72 | 22.50 | **12.92** | 23.66 | 25.20 |
| equal lengths (fle200) stranded ON | 9.24 | 6.82 | **6.78** | 8.41 | 8.04 |
| gDNA long (flgl) stranded ON | 12.83 | 12.39 | **9.07** | 9.76 | 12.08 |
| sparse probes stranded ON / strand 0.70 ON | 74.76 / 55.26 | 24.37 / 27.47 | 24.43 / **26.08** | **21.20** / 30.88 | 22.36 / 33.57 |
| junction-spanning probes stranded ON / 0.70 ON | 40.01 / 40.88 | 32.16 / 32.41 | 32.05 / 32.37 | **23.39 / 21.22** | 23.09 / 24.11 |
| single probe, stranded ON (whole library) | 9.88 | 8.75 | 8.82 | 11.46 | 9.22 |

### The probed gene itself on the single-probe panel

This is the case your spectrum ruling was about. Net mRNA error on the one targeted gene (gB1_cap) at g50:

| | 10-05 | Foundation | Main | Candidate |
|---|---:|---:|---:|---:|
| stranded | +64.7 % | +64.7 % | −11.0 % | −8.7 % |
| strand 0.70 | +116.2 % | +116.7 % | −2.9 % | −0.1 % |

- **Old reader:** it finds no captured population with one probe, so it applies no capture. Captured gDNA in the
  targeted gene is called RNA.
- **New reader:** it reads the probe locally.
- **At g05:** the old reader's error on that gene is 1–7 %, the same range as the new reader's.
- **Larger panels:** on the sparse panel's 112 probed genes the old reader is within about a point of the new one.

### Real libraries

The VCaP mix carries per-read truth. DNA recovery is the mix's gDNA minus the transcriptome-only run's, over the DNA
reads' fragments.
- 10-05 recovers 91.5 % (494 k DNA fragments leak to nascent RNA).
- The candidate recovers 98.5 %.
- Most of that gain is the multimapper second pass, which is independent of the reader: the reader switch alone
  moved VCaP's gDNA fraction 0.2395 → 0.2401, and the second pass moved it to 0.2528.

## 5. What was progress and what was not

**Progress, keep:**
- **The count foundation.** It fixed the original fragment-length problem and is neutral or better almost
  everywhere. Its one cost is the deferred stratum on RNA-long (6.9 → 43.2 %) and on the sparse and junction-probe
  layouts.
- **The multimapper second pass**, the real-data truth gain.
- **The junction clip.**
- **The strand and constants work.**
- **The instrument fixes.**

**Right principle, wrong trade so far:**
- **The per-object reader.** Small panels do need local capture, and it delivers that. But on the matrix it costs far
  more than it gains: capture-off noise, strand noise, weak strand, the ladder's stranded ON and deferred rows.

**Not progress:**
- **The repair variants since 10-11**, including today's pooling. Pooling broke the probe-placement physics outright:
  probes mark specific exons, and a gene's exons need not share a capture level.

**The intron factory, specifically.**
- **What it was:** it gave every intron a composition claim equal to the library's intergenic gDNA background.
- **Why it was removed:** that is a library-level prior presented as local evidence, counted on top of the landscape
  prior. Under capture an intron within a fragment of a probe is not at background.
- **What removing it cost:** today's numbers say removal was mostly right. The cost is confined to unstranded ×
  capture-ON, the deferred stratum.

## 6. Why it became a snowball

1. **Two changes landed on the same day without an A/B of each alone on the full matrix**: the count foundation and
   the reader.
2. **The next day's attribution was wrong in both directions.** So the reader looked necessary for the RNA-short fix,
   and the factory looked guilty for the deferred collapse.
3. **The baseline moved to main.** Every repair was judged against a degraded tree. The last good tree was never
   carried as a standing column, and the robustness panels were built after it.
4. **The aggregate scores hid the case that motivated the rewrite.** The single probe's target gene doesn't show up in
   a whole-library sum.

Lessons:
- Carry the last good tree in every A/B.
- Land one mechanism at a time against the full matrix.
- Score targeted genes separately from the library.

## 7. The way to a release

Only one component stands between the foundation and a release: **capture for very small panels.** Two options.

**A. Release 0.8.0 from the foundation.**
- **Base:** `3e7e0fbc`, plus the changes since that are independent of the reader: the multimapper second pass and
  the ruler instrument fix. The per-object reader and its speed and repair work leave the release branch.
- **What it gives:** the strongest tree we have measured on every in-scope stratum. It solves the problem we set out
  to solve.
- **Its documented limits:**
  - A library whose captured population is too small to locate (a handful of probes) gets no capture correction;
    measured, the targeted gene is +65–117 % at g50 and +1–7 % at g05.
  - Junction-spanning probe layouts (32 % against main's 23 %).
  - The deferred stratum on RNA-long.
- **The catch:** it keeps the library-level switch your 10-05 ruling rejected. That is the decision.

**B. Hold 0.8.0 for one new reader that is both local and robust.**
- **Acceptance test:**
  - no worse than the foundation on any stratum of today's matrix, within seed noise;
  - the single probe's target gene correct;
  - RNA-short stranded off below 5 %.
- **What it takes:** a better idea than the variants already refused; it starts from what each existing reader does
  right. Time-box it, and fall back to A if it does not pass.

**Recommendation: A for 0.8.0, then B as 0.8.1's single project, with the foundation as its bar.**
- A ships the robust tree we already have, and states its one principled gap with measured numbers.
- B then adds locality without being allowed to give robustness away again.

If the small-panel case must ship in 0.8.0, then B with a short time box.

Either way, nothing else gets built until this choice is made.
