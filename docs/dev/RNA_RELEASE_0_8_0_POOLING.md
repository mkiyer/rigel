# Rigel 0.8.0 — capture shrinkage by exonic-footprint pooling (2026-10-12)

Sandbox, not authoritative. This responds to the owner's comments on the stage 1–4 decision in
`RNA_RELEASE_0_8_0_RESULTS.md`:
- shrink continuously to uniform capture weights when there is too little gDNA to set them;
- accept short-RNA if it degrades gracefully;
- weak-strand and unstranded captured rows must degrade gracefully;
- run a competitive head-to-head.

Prototype, receipts and outputs are in `$RIGEL_RUNS/prototypes/2026-10-12_capture_shrinkage/` (README).
Nothing changed in `src/` on main, and nothing was committed.

Protocol for every Rigel arm: scan pinned, fractional assignment, the candidate reader (`split·DNA`, policy 3)
underneath. Each arm changes only the capture weights after `capture_weights`. The tag is v0.7.1, main is `624dfade`,
and the candidate is stage 3's `split·DNA`.

## 1. The failure the owner pointed at is a consistency failure, not a variance one

Case: lowg g001, stranded, capture ON, 1,170 gDNA fragments in 1.17 M. Score: candidate 39.71 %, tag 7.68 %.

**What happens in a probed gene.**
- About 16 true gDNA fragments sit across the gene's exonic objects.
- The RNA-free boundaries hold 1–2 reads each. The reader resolves them to the captured level, about 10,000× the
  typical object.
- The RNA-rich exons hold about 3 gDNA reads each among about 14,000 RNA reads. The strand channel cannot see them,
  so they read the landscape's spike. That spike is a library-level fallback.
- Junction prices inherit the boundary level. Each transcript's length becomes junction-dominated.

**Why it is catastrophic.** Isoforms with nested ends need their length ratios right to about 1 % (main collapses
an isoform at a 1 % error). The ratios become 1 : 0.99999 : 0.99995 against a truth of 1 : 0.978 : 0.926, and the
isoforms collapse.

**Non-monotone in gDNA.** With less gDNA (d10 g001) nothing reads, the field stays uniform, and the row is fine
(13.7 %). With more gDNA the exons resolve and the field is consistent again.

**How the obvious forms of "shrink to uniform" fare** (transcripts %, genes % in brackets; full tables in the
prototype):

| Arm | Rule | lowg g001 ON | d10 g01 ON | d10 g50 ON | test 0.99 ON contaminated |
|---|---|---:|---:|---:|---:|
| candidate | — | 39.71 | 38.29 | 21.45 (2.88) | 8.27 (0.93) |
| uniform weights | every weight 1 | 7.24 | 16.30 | 41.60 (29.32) | 24.99 (19.85) |
| Dirichlet | per locus, λ = K/(K+O), one pseudo-fragment per object | 36.50 | 37.36 | 21.68 | — |

Uniform weights fix scarce gDNA and break every captured row with gDNA. A count-based rate barely moves the failure:
16 reads on the boundaries are real information, and the defect is that the exon pieces beside them are blind.

## 2. The mechanism: pool a gene's exonic footprint, with a spread its own evidence sets

The owner's 2026-10-07 footprint ruling says enrichment is shared by a region, its boundaries and the connected
junctions, while an intron's interior lies beyond a probe's reach. That fixes the pooling unit: per locus (the run
between message terminals), the exon regions and every boundary beside one. Introns and intergenic regions keep the
reader's mode.

**Per member o.**
- m_o is the reader's posterior mode under the landscape, unchanged; it carries the landscape's veto of
  strand-noise peaks.
- I_o is the likelihood's curvature at m_o: the data's own information at the readout. It is 0 where the likelihood
  is flat there, as for an RNA-rich exon or a bound lying far above the mode. It is positive where reads or a binding
  bound constrain it.

**Per group.**
- ĝ = Σ I·m / Σ I.
- τ̂² by Paule–Mandel on (m, 1/I) over the informative members. It is 0 when they agree within their own noise or
  fewer than two are informative.
- Each member reads the mode of log L_o(x) + log π(x) − (x − ĝ)²/(2τ̂²), its own likelihood under the landscape
  localised to its group. At τ̂ = 0 every member is at ĝ.

**Limits.**
- No gDNA in a gene, or no disagreement beyond noise: the gene is uniform. This is the owner's principle, applied
  continuously.
- Abundant gDNA with real exon-specific capture: the resolved members keep their own levels.
- A blind member of a dispersed group falls back toward the landscape, in proportion to the spread.

**Constraints.** No detector, no library reference, no new constant. The spread is a standard random-effects moment
estimator; the τ and grid lattices are numerical only.

## 3. Results

### Graceful degradation toward unstranded

Seven test panels (test, odg05, fle200, bind01, bind1, bind100, single), contaminated conditions, ratio of sums;
transcripts % (genes %).

| Stratum | Tag | Main | Candidate | **Pooling (PM)** | salmon | kallisto |
|---|---:|---:|---:|---:|---:|---:|
| strand 0.99 ON | 10.64 (3.04) | 10.90 (1.54) | 10.94 (1.60) | **9.34 (1.55)** | 54.04 | 47.34 |
| strand 0.70 ON | 13.69 (4.45) | 12.85 (2.68) | 16.66 (2.97) | **10.44 (2.78)** | 64.96 | 87.46 |
| unstranded ON | 19.48 (11.29) | 15.60 (4.46) | 23.38 (10.13) | **10.12 (3.91)** | 73.93 | 87.23 |
| strand 0.99 OFF | 8.60 | 9.84 | 8.32 | **7.75** | 12.58 | 10.97 |
| strand 0.70 OFF | 8.78 | 9.16 | 8.21 | **8.13** | 25.00 | 15.43 |
| unstranded OFF | 9.04 | 8.55 | 8.24 | **8.12** | 14.18 | 15.06 |

### Scarce gDNA under capture

| Row | Tag | Main | Candidate | **Pooling (PM)** |
|---|---:|---:|---:|---:|
| lowg g001 ON | 7.68 | 38.71 | 39.71 | **5.97** |
| lowg g01 ON | 9.76 | 8.47 | 10.56 | **5.72** |
| d10 g01 ON | 19.23 | 38.06 | 38.29 | **16.26** |
| d10 g05 ON | 23.12 | 19.60 | 23.84 | **13.22** |

### Exon-specific probe designs: the remaining cost

| Panel (capture ON) | Tag | Main | Candidate | **Pooling (PM)** |
|---|---:|---:|---:|---:|
| sparse 0.99 / 0.70 | 74.76 / 55.26 | 21.20 / 30.88 | 22.36 / 33.57 | **22.16 / 29.12** |
| junction 0.99 / 0.70 | 40.01 / 40.88 | 23.39 / 21.22 | 23.09 / 24.11 | **28.42 / 27.50** |

The junction-probe genes carry real within-transcript capture structure: probes span junctions, so an isoform's
end zone is uncaptured. Pooling still flattens part of it. The cost is +3 to +5 against the candidate and main;
pooling stays below the tag.

### The ladder and the fl-gap panels

The mode-based arm with the DerSimonian–Laird spread was measured first. PM is the same mechanism with the better
spread estimator.

| Ladder stratum | Tag | Main | Candidate | Pooling (mode-based, DL) |
|---|---:|---:|---:|---:|
| stranded OFF, contaminated / g00 | 3.95 / 3.70 | 3.37 / 13.61 | 2.62 / 3.81 | **2.20 / 2.08** |
| stranded ON, contaminated / g00 | 9.91 / 20.22 | 6.05 / 15.10 | 6.24 / 6.61 | 6.29 (genes 2.05 vs 1.39) / 6.58 |
| unstranded OFF, contaminated / g00 | 10.29 / 20.09 | 2.40 / 1.59 | 2.34 / 1.59 | 2.35 / 1.60 |
| unstranded ON (deferred), contaminated | 33.21 (14.25) | 62.96 (23.69) | 52.73 (19.53) | **11.34 (7.62)** |

(PM ladder, fl-gap and real-library rows: §7.)

## 4. Refused variants, each with the number that refused it

| Variant | What it did | Refusing measurement |
|---|---|---|
| per-object ι shrinkage | each object toward its locus's evidence-weighted level | breaks the intermediate rows: lowg g01 ON 10.56 → 23.61, d10 g05 ON 23.84 → 41.59 |
| harmonic fill along the chain | blind objects take their neighbours' level | ramps between anchors; the landscape's spike counts 0-read objects as resolved |
| local-curvature fill | fill from the likelihood's own curvature | strand noise makes false peaks (an exon with 0 gDNA peaks at −3.0) |
| pooling with DL spread on pooled likelihoods | spread by DerSimonian–Laird | ladder stranded ON contaminated 7.97 (g98 57.8) |
| pooling with type-II ML spread | spread by type-II maximum likelihood on pooled likelihoods | ladder unstranded OFF 2.34 → 5.03; summing members' RNA-integration slopes makes false group levels of ±2 nats |
| joint ML of level and spread | the landscape cancels at τ = 0 | noise peaks set the level: test g00 ON 6.81 → 19.29 |
| footprint centring in dispersed groups | blind members centred on neighbours | sparse 0.99 ON 25.67 → 33.13 |
| two-number (Gaussian) readout | (I·m + ĝ/τ²)/(I + 1/τ²), no curve | junction 0.99 ON 36.61, sparse 0.99 ON 26.34; blind members of dispersed groups must fall back toward the landscape |

## 5. Competitive head-to-head

Tools: salmon 1.10.3 (`-l A`) and kallisto 0.51.1 (`--rf-stranded` only at strand 0.99). They ran on each test
condition's simulated FASTQs and were scored by `quant_accuracy.py`'s own transcript and gene metric.

**gDNA present.** Rigel is far ahead: 7 % against 47–57 % at test 0.99 ON contaminated, and 10 % against 74–87 %
unstranded ON.

**No gDNA.** salmon and kallisto lead on transcripts by 1–4 points; Rigel leads on genes. The entire transcript gap
is the `cluster` families: ten isoforms whose ends differ by 1–2 bp, with equal shares in the truth. Rigel's EM
concentrates mass there, while salmon's prior and kallisto's EM stay near equal shares. Excluding those families,
Rigel is best:

| Zero gDNA, excluding cluster families | Rigel | salmon | kallisto |
|---|---:|---:|---:|
| test 0.99 OFF | **3.19** | 4.07 | 5.04 |
| test 0.99 ON | **2.48** | 3.76 | 4.77 |
| test unstranded OFF | **3.03** | 3.58 | 5.15 |

This is `ISSUES: per-transcript-prior-lane`, an EM decision outside capture.

**Weak strand.** salmon's automatic library type picks a stranded protocol at strand 0.70 and drops to 29 % at g00.
kallisto run unstranded stays at 7.6 %. Rigel degrades smoothly.

## 6. What remains for the owner

1. **Pooling is not ready to port as built.** It fixes the scarce-gDNA, weak-strand, unstranded and short-exon
   classes on every simulated panel. On real captured libraries it moves gDNA into RNA (§7). Pooling exon regions
   only removes the part into nascent RNA but not the part into mRNA, and gives back the simulated gains (§8). The
   next step is a modelling ruling: a boundary's capture derived from its exon's and its intron's, plus a per-locus
   dissection of the leak into mRNA.
2. **Port, once a variant passes.** Pooling needs each exonic slot's likelihood at a second centre. Groups never
   cross a native block, so it fits inside the reader's block loop: fit each group after the loop, then re-read only
   the members of dispersed groups. The two-number shortcut is refused (§4).
3. **The EM on near-identical isoforms** (§5). This is the only place salmon and kallisto beat Rigel, and an owner
   ruling on the per-transcript prior. It is not part of the capture work.

## 7. Release-panel and real-library results for the PM arm

### Ladder and fl-gap panels

| Ladder stratum | Tag | Main | Candidate | **Pooling (PM)** |
|---|---:|---:|---:|---:|
| stranded OFF contaminated / g00 | 3.95 / 3.70 | 3.37 / 13.61 | 2.62 / 3.81 | **2.20 / 2.08** |
| stranded ON contaminated | 9.91 (2.52) | 6.05 (1.37) | 6.24 (1.39) | **5.95** (1.88) |
| stranded ON g00 | 20.22 | 15.10 | 6.61 | 6.58 |
| unstranded OFF contaminated / g00 | 10.29 / 20.09 | 2.40 / 1.59 | 2.34 / 1.59 | 2.35 / 1.60 |
| unstranded ON (deferred), contaminated | 33.21 (14.25) | 62.96 (23.69) | 52.73 (19.53) | **10.97 (7.28)** |
| RNA-short 0.99 ON / unstranded ON | 20.33 / 29.76 | 14.24 / 68.06 | 20.09 / 63.65 | **14.79 / 16.30** |

The stranded ON gene cost sits at g05: genes 0.84 → 1.59. At g50 and g98 PM matches the candidate.

### Real libraries: a systematic shift of gDNA into RNA

VCaP DNA recovery is measured as the mix minus the transcriptome-only run, over the DNA reads' 4.66 M fragments.

| Library | Candidate | Pooling (PM) |
|---|---:|---:|
| VCaP mix: DNA fragments recovered as gDNA | 98.5 % | **94.6 %** (to nascent RNA +206 k vs +63 k; to mRNA +45 k vs +8 k) |
| VCaP mix: gDNA fraction (origin ≈ 0.2518) | 0.2533 | 0.2438 |
| cfRNA-B (94 % gDNA): mRNA | 41,003 | **88,267** |
| cfRNA-A / cfRNA-C: gDNA fraction | 0.0990 / 0.1651 | 0.0945 / 0.1566 |

**Not artifacts.** The transcriptome-only library's pooled lifts lie 99.7 % inside the mix's own pooled footprint.
Pooling reads the same exome capture in both libraries; the candidate only ever saw the objects it could resolve
individually.

**The cost is a modelling error in what pooling equates.** In the EM a locus's gDNA density inside its exons is read
off its boundaries' and introns' gDNA, scaled by their capture-weighted opportunities. An exon|intron boundary is
less captured than the exon interior, because a crossing fragment half-overlaps the probe.
- Pooling sets the boundary to the interior's level. The boundary reads then imply too little gDNA inside the exon,
  and the surplus goes to RNA. FLG and SOX9 in cfRNA-B go from 0 to about 200 mRNA fragments each.
- The candidate errs the other way. It reads those boundaries at the typical object's level, so it over-infers
  exonic gDNA.

The ladder's g05 stranded ON gene cost is the same mechanism in simulation.

## 8. Pooling exon regions only (boundaries inform, keep their own readings)

| Measure | Candidate | Pooling (PM) | Regions only |
|---|---:|---:|---:|
| lowg g001 ON / d10 g01 ON | 39.71 / 38.29 | 5.97 / 16.26 | 25.05 / 25.20 |
| captured strand 0.70 / unstranded (6 panels) | 16.66 / 23.38 | 10.44 / 10.12 | 13.47 / 24.90 |
| ladder stranded OFF g00 | 3.81 | 2.08 | 4.48 |
| VCaP DNA recovered | 98.5 % | 94.6 % | 97.3 % |
| VCaP leak to nascent / to mRNA | 63 k / 8 k | 207 k / 45 k | 77 k / 48 k |
| cfRNA-B gDNA fraction / mRNA | 0.9373 / 41 k | 0.8814 / 88 k | 0.8889 / 62 k |

Regions-only pooling gives back most of the simulated gains. The consistency they need ties the exon|intron boundaries
to their exons, because the junction prices are built from those boundaries. It also halves VCaP's leak, removing the
part into nascent RNA, which comes from the pooled boundaries. The part into mRNA remains: lifting blind exon pieces
to their gene's level leaks gDNA to mRNA on a real captured library.

So neither "pool the boundaries with their exons" nor "leave them out" is right. A boundary's capture lies between its
exon's and its intron's. The coherent next step is to derive the boundary's level from those two, the relation the
junction price already assumes, and to find why lifted exon pieces send gDNA to mRNA. That needs a per-locus dissection
on a real captured library, or on a simulated analogue with cfDNA-like fragment lengths that reproduces the shift.
It is a modelling ruling for the owner, not a variant to tune.
