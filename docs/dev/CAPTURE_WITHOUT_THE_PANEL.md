# Capture the data cannot see — junction-spanning probes, the isoform split, and what to do without the panel

A review for the owner, written 2026-09-19 after the stranded × capture-ON dissection. It is a sandbox note:
provisional, cited by nothing, and meant to be argued with. The permanent record of what was measured is
`ISSUES: ruler-witness-geometry-on-transcript-panels` and `ISSUES: em-overturns-the-calibrated-gdna-split`; the
data and scratch tools are in `~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/`.

**The constraints this review works under** (owner, 2026-09-19): probes often span splice junctions; the probe
panel is usually NOT available; stranded × capture-ON should get below 5 %.

**STATUS (2026-09-19, later): the physics question of §4 is SETTLED.** The owner ruled that gDNA and cDNA
bind alike — a split probe's geometry is gDNA's whole disadvantage, since only a transcript holding the junction
can bind the probe whole — and the simulator now binds every fragment through one contiguous part of a probe,
with no penalty (`docs/TESTING.md`, the capture row). The ladder was deleted and is being re-simulated under it;
§7's first step is done, and what remains below is the question of junction-targeted designs and the ideas of §6.

---

## The short version

1. **What the EM needs is one number per transcript** — its capture-aware effective length — and only the
   RATIOS of those numbers inside a locus matter. A 10–20 % error in one isoform's number relative to its
   sibling's decides whether a silent isoform receives tens of thousands of fragments.
2. **A probe across a splice junction makes capture isoform-specific**, and the gDNA witness the shipped ruler
   reads cannot see a junction (gDNA has none), while at zero gDNA there is no witness at all. With the
   simulator's own lengths the stranded × capture-ON error falls from 10.5 % to 3.6 % at `g05` and from 7.4 %
   to 2.9 % at `g00`. Averaging those true lengths within each gene throws the gain away (8.4 %), so wherever a
   correction is needed it has to be isoform-resolved: no gene-level correction can reach 5 %.
3. **Most of the ladder's error is a physics assumption, not a property of junction probes as such.** The
   simulator binds a half-matched probe at full strength when the fragment is cDNA and at a fifth when it is
   gDNA — an asymmetry I know no physical reason for. Re-simulated with the two binding alike, the SHIPPED ruler
   reads 4.7 % at `g05` (from 10.5 %) and 8.2 % at `g50` (from 13.9 %): below the 5 % target at `g05` with no
   change to Rigel. On a panel that targets junctions deliberately (the test chromosome's twin) the problem
   survives symmetric binding. So the real-world risk is set by two facts about real panels — the half-match
   physics, and how often a probe sits on an isoform-discriminating junction — and both should be settled before
   anything is built.
4. **The information needed is in the data, without the panel.** RNA coverage has a characteristic shape around
   every probe; spliced reads are certified RNA and carry the junction-spanning capture directly; gDNA marks
   where probes sit even where it misreads their strength; and a kit captures identically in every library
   that uses it, which a cohort can exploit. If the asymmetry is real, or real panels target junctions, the build
   is a **capture field learned from the data** — one per-base capture efficiency, from which RNA, gDNA and
   nascent lengths are all derived consistently — with a cohort-level "panel profile" as the product path. In
   every case a sparsity prior on isoforms is the safety net, and the `g50` residual is the EM's gDNA split.

---

## 1. What the EM needs, and why a small error costs so much

A fragment `f` compatible with transcripts `t₁…tₖ` is scored against each with likelihood `P(f | t) = c(f) / E_t`,
where `c(f)` is the fragment's capture weight — a property of its own sequence — and `E_t = Σ_{f' ⊂ t} c(f')` is
the transcript's capture-aware effective length. The posterior over the compatible transcripts is
`θ_t c(f)/E_t ÷ Σ_{t'} θ_{t'} c(f)/E_{t'}`, and `c(f)` cancels. So **the EM needs `E_t` and nothing else**, and only
its ratios within a locus (and against the locus's gDNA component and nascent entities). This is why the
`oracle_ruler` arm — the simulator's `E_t` handed to the shipped EM — brings the stranded capture-ON rows to the
capture-OFF level: the rest of the machinery is already right.

The cancellation has one condition that matters later: `c(f)` must be the same whichever component produced
`f`. That holds among RNA isoforms (a fragment's sequence is the same in each), and it fails between an isoform
and any other component that can hold the same unspliced fragment — gDNA, a nascent entity — wherever the two
are captured differently at the same position.

**The existence test.** Take a silent isoform `A` whose fragments are all compatible with an expressed isoform
`B`, and let `n_s` be the fragments in the shared region and `n_u` those in `B`'s unique region. The EM's fixed
point gives `A` mass exactly when

    n_s / E_A  >  n_u / (E_B − E_A),

that is, when the shared region's fragment density per effective base exceeds the unique region's. With plain
lengths CACNA1I reads 5,725 effective bases shared against 3,940 unique, and `A` correctly gets nothing. The
shipped ruler contracts the shared region to 3,605 while leaving the unique one at 3,935; the shared density now
looks higher, and the silent isoform takes 38,493 fragments. The true factors of the two isoforms differ by 3 %
(0.872 / 0.900); the shipped ones by 19 % (0.630 / 0.780).

So the yardstick for any capture correction is **the within-gene spread of its error**, not its per-transcript
accuracy (`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The shipped ruler is closer to the truth transcript by
transcript than no ruler at all, and worse end to end.

## 2. What capture does to `E_t`

A captured fragment's weight rises with how much of a probe its sequence matches. Around a probe that spans the
junction `J` between exon `a` and exon `b`, three kinds of fragment compete:

| fragment | matches | carried by |
|---|---|---|
| spans `J` | the whole probe | only the isoforms that hold `J` |
| ends in `a` (or starts in `b`) without spanning `J` | one half of the probe | every isoform holding that exon end — and the nascent span |
| gDNA over the exon end | one half of the probe | gDNA |

Isoform specificity lives in the first row. The gDNA witness sees only the third. An ordinary probe inside an
exon is isoform-shared: every isoform holding the exon gets the same boost, and gDNA sees it — so the shipped
design, which reads capture off gDNA piece by piece, is right on a panel whose probes never cross a junction.
That is what the test chromosome's benign panel shows (§3).

## 3. What was measured

All numbers: transcript-level Σ|Δ| as a share of the true RNA, fractional assignment, stranded × capture-ON.
(Σ|Δ| counts a misplaced fragment twice — once where it went, once where it belonged — so 10 % is about 5 % of
fragments on the wrong transcript.)

**The ladder** (its panel: 13,704 probes placed along transcripts, 24 % spanning a junction):

| lengths the EM divides by | `g00` | `g05` | `g50` |
|---|---|---|---|
| shipped ruler | 7.4 | 10.5 | 13.9 |
| no ruler (plain lengths) | 7.4 | 7.8 | 10.1 |
| shipped formula, certified true gDNA counts | — | 13.0 | 14.8 |
| shipped ruler averaged within each gene | 7.4 | 8.7 | 14.1 |
| the simulator's lengths averaged within each gene | 7.7 | 8.4 | 12.8 |
| **the simulator's own lengths (`oracle_ruler`)** | **2.9** | **3.6** | **8.4** |
| re-simulated with gDNA binding a half-match as cDNA does (§4), shipped | — | **4.7** | **8.2** |
| the same, `oracle_ruler` | — | 3.4 | 8.6 |
| capture-OFF, for reference (stranded) | 3.3 | 2.9 | 3.9 |

Read it four ways. A perfect gDNA witness makes things WORSE — the witness measures the wrong quantity. The
truth averaged within genes is no better than no ruler — the value is entirely in isoform-resolved lengths. Let
gDNA bind a half-match as cDNA does and the shipped ruler recovers almost all of it at `g05` — the "wrong
quantity" is mostly the simulator's asymmetry (§4). And the error concentrates: at `g05` 14 genes carry half of it, genes with four or more isoforms 93 %, their totals
right to a fraction of a percent. What the true lengths leave at `g50` is mostly the EM giving captured gDNA to
isoform-rich genes against a calibration that had the split right (`ISSUES: em-overturns-the-calibrated-gdna-split`).

**The controlled case — the test chromosome, the same libraries under different panels:**

| stranded × capture-ON | `g00` | `g05` | `g25` | `g50` |
|---|---|---|---|---|
| capture-OFF | 7.5 | 7.2 | 7.3 | 9.5 |
| benign panel (4 of 2,469 probes cross a junction), shipped | 7.3 | 5.3 | 7.1 | 11.9 |
| junction-probed twin, shipped | 40.8 | 34.9 | 34.9 | 35.7 |
| junction-probed twin, `oracle_ruler` | 14.0 | 12.5 | 16.6 | 21.3 |
| junction-probed twin with symmetric physics (§4), shipped | 38.8 | 33.8 | 34.3 | 35.1 |
| junction-probed twin with symmetric physics, `oracle_ruler` | 11.3 | 11.0 | 11.7 | 16.6 |

(The test chromosome's isoform structures are hand-built stress cases, which is why its capture-OFF floor is
higher than the ladder's.) With no junction-spanning probes, capture costs nothing. With a panel that targets
junctions, the shipped ruler loses ~30 points at every gDNA level, including `g00`, under either physics.

## 4. The physics question the ladder cannot answer

The simulator's capture model (`rigel.sim.capture.sampler`): a fragment weighs
`off_target_weight + binding_per_base × overlap`, the overlap being the fragment's match to its best probe; a
probe's pieces that sit contiguously in the fragment's template bind at full scale, and pieces separated by an
intron bind at `gdna_split_penalty` = 0.2. Walk that through the table of §2:

* cDNA spanning `J` — contiguous, full scale, overlap up to 120 bp.
* cDNA holding one exon end — a single contiguous piece, **full scale**, overlap up to 60 bp.
* gDNA over the same exon end — a piece of a split probe, **0.2**, overlap up to 60 bp.

So the same 60-base half-match binds five times harder when it is cDNA than when it is gDNA (or nascent RNA). I
know of no physical reason for that asymmetry: at hybridisation both are denatured library fragments offering the
same 60 complementary bases. Three regimes are possible, and they call for different fixes:

| regime | half-match binds | gDNA witness at a junction probe's flanks | junction-specific boost |
|---|---|---|---|
| **as simulated** | cDNA full, gDNA ⅕ | reads a fifth of the RNA's capture — wrong for every isoform | modest (≤ 2×, linear in overlap) |
| **symmetric, strong half-match** | both nearly full (binding saturates with overlap) | right | small |
| **symmetric, weak half-match** | both weak | right | large — and invisible to gDNA |

Two measurements already bear on it:

* **In the simulator the junction-specific boost is modest.** In pure-RNA capture-ON single-isoform genes a
  junction under a junction-spanning probe sits +0.04 (log, median) against its gene's junctions, one near an
  ordinary exon probe −0.02, one with no probe nearby −0.11; separating junction-probed from exon-probed junctions
  by their spliced-read count reaches an AUC of 0.60 (`junction_witness_feasibility.txt`).
* **On a junction-targeted panel the problem survives symmetric physics.** Setting `gdna_split_penalty` to 1.0
  on the test chromosome's junction twin moves the shipped ruler by 1–2 points (§3): there the isoform-specific
  boost dominates, and no gDNA witness can see it under any physics. The true ruler improves by 2–5 points under
  symmetric physics — at `g00` too, where there is no gDNA — because the simulator's rule also applies to the
  nascent entities' split pieces: the asymmetry breaks the cancellation of §1 between an isoform and every
  component that holds the same unspliced fragment, gDNA and the nascent span alike. What the true ruler still
  leaves on the junction twin (11–17 % against the benign panel's 5–12 %) is information a junction-targeted
  panel spends: it concentrates the library on some features and starves the rest.

* **On the ladder the problem is mostly the asymmetry.** Two ladder rungs re-simulated with `gdna_split_penalty`
  1.0 (the same seed and abundances; `~/Downloads/rigel_runs/suite/ladder_split1/`): the shipped ruler reads
  4.7 % at `g05` and 8.2 % at `g50` against 10.5 % and 13.9 % as simulated, while the true ruler barely moves
  (3.4 % and 8.6 %). With the witness able to read the half-matched flanks, the gDNA-based ruler does its job on a
  panel whose junction probes fall where a random tiling puts them; what `g50` keeps is the EM's gDNA under-call
  (library gDNA 0.479 against 0.50), not the ruler.

So the three regimes of the table are not equally likely to hurt. Under symmetric binding the gDNA witness is
right at every flank and the only thing it cannot see is the junction-specific boost, which matters in
proportion to how often a probe sits on a junction that discriminates isoforms — on the ladder's random tiling
it costs 1.3 points at `g05` (4.7 % against the true ruler's 3.4 %), on the test chromosome's twin, built to put
probes there, it is most of the error, and on a real panel it is unknown.

**What would settle the physics.** Your knowledge of the capture chemistry first: 120-mer baits, hybridisation
temperature, whether 50–60 bp partial matches are captured near-fully (the "near-target" reads in exome data
suggest partial overlap captures efficiently, but I have not verified the literature). Then data: if the probe
design of even ONE panel is known — your own cfRNA design, say — a captured library with gDNA gives capture as a
function of overlap for cDNA and gDNA directly. Without any design: in a real captured library with gDNA, compare
the gDNA enrichment at exon ends flanking a strongly enriched junction with the gDNA enrichment inside the same
gene's probed exon bodies. And in the simulator, the half-match rule should be ONE parameter applied to cDNA and
gDNA alike, so the ladder can be run across the plausible range rather than at one unvalidated corner of it.

## 5. What the data contains about capture, with no panel

| witness | what it sees | what it cannot see |
|---|---|---|
| **gDNA coverage** | the genomic footprint of every probe, junction probes included (as paired half-footprints at two exon ends across one intron) | a junction probe's strength for cDNA; anything at `g00`; much at `g05` in one library (sparse) |
| **RNA exon coverage — its level** | relative capture of pieces shared by all isoforms of a gene (the isoform mix is the same across them, so it cancels) | capture of isoform-unique pieces, which is confounded with the isoform's abundance |
| **RNA exon coverage — its shape** | where probes sit, from shape alone: a probe leaves a bump of a known form (the probe interval convolved with the fragment-length distribution), while an isoform change leaves a step at an exon boundary | a junction probe's bump straddles the boundary and must be told from a step by shape |
| **spliced reads** | certified RNA (no gDNA confound); a junction's own capture; the offsets of junction-spanning fragments around the junction | abundance and capture multiply at an alternative junction |
| **redundancy** | an isoform with several distinguishing features (unique exons, unique junctions) has several estimates of its abundance; a probed feature is an outlier among them | an isoform with one distinguishing feature |
| **the cohort** | a kit captures identically in every library that uses it, while isoform usage varies between samples — so capture is the stable part | a single library |
| **the annotation** | where every junction and exon end is, which features are isoform-specific, which isoforms nest | where probes are |

The strongest single fact: **probe locations are recoverable from shape**, and with a physical model with a
handful of library-level parameters, locations are all `E_t` needs.

## 6. Ways to overcome it

Each idea: why it could work, what it needs, how it fails, how it would be judged. `oracle_ruler` is the ceiling
for every one of them; the test chromosome's benign and junction-probed twins (and their symmetric-physics
variants) are the controlled substrate; `quant_accuracy.py --set em.assignment_mode=fractional` and
`ruler_vs_truth.py` (read for its within-gene spread) are the instruments.

### 6.1 Settle the physics first — the prerequisite, not a fix

Make the simulator's half-match rule physical and parameterised (§4), re-simulate the test chromosome's twins and
two or three ladder rungs across its range, and let every candidate below be scored across it. Cost: a
simulator change and a few hours of simulation. Without it, the next build may be tuned to an artefact.

### 6.2 Stop the harm where the witness is blind

Under the simulator's asymmetric physics the shipped ruler loses to no ruler on a junction-spanning panel; under
symmetric binding it does not (§4), so this matters only if the asymmetry is real or a panel targets junctions.
Cheapest: do not let the gDNA witness create
within-gene differences it cannot support — shrink each isoform's factor toward its gene's in proportion to the
witness's reliability, and detect junction-spanning panels from their gDNA signature (paired half-footprints
flanking introns). Measured bound: the shipped ruler averaged within genes reaches 8.7 % at `g05` (from 10.5 %);
no ruler 7.8 %. **It cannot reach 5 %** — the truth averaged within genes is 8.4 % — so this is a guard, not a
solution. It also shows the contraction shifting the RNA/gDNA balance at `g50` (14.1 % gene-averaged against
10.1 % plain), which 6.3 has to get right as well.

### 6.3 A capture field learned from the data — the recommended build

Estimate ONE per-base capture efficiency field `c(x)` over the transcriptome's sequences — exon bodies shared by
every isoform that holds them, junction windows belonging to the isoforms that hold the junction — and derive
every component's length from it: each isoform's `E_t`, the nascent entities', and the gDNA component's (so RNA
and gDNA cancel as §1 requires). Parameterise `c` physically: an off-target level, a binding strength and a
half-match response shared by the whole library, plus a sparse set of probe locations. Fit it jointly with the
abundances, in the way Cufflinks' fragment-bias correction, alpine and Salmon's bias models learn a bias from the
data and fold it into effective lengths — capture is a fragment-level bias whose covariate (probe overlap) is
latent rather than computable, which is the new part.

Why it could work: probe locations are identifiable from coverage SHAPE whatever the abundance (§5); exon-body
capture is identifiable from constitutive pieces; junction-window capture from spliced reads and their offsets;
and the library-level parameters are fitted on thousands of probes at once. Where an isoform has no reads, its
unobserved junction windows take the library's prior probability of a junction probe — the same empirical-Bayes
move the gDNA landscape makes.

What it needs: positional coverage (a per-region profile at modest resolution, built during the scan) and the
spliced-fragment offsets around each junction; neither is in the payload today. How it fails: shallow
libraries (bumps undetectable), probe-to-probe efficiency variation (a second latent per probe, shrunk to the
library's distribution), dense tiling where bumps merge into plateaus (then levels, not shapes, carry it), and
the physics regime (the half-match response must be fitted, which 6.1 makes testable). Judge: within-gene spread
against the simulator's truth, then the transcript table against `oracle_ruler`, on both test-chromosome twins
and the ladder, halves apart.

A cheaper first prototype: two passes. Quantify with plain lengths, fit the capture field to the residual pattern
of each gene's pieces and junctions (a probed feature is a positive residual shared by every isoform using it;
an abundance error is not), recompute `E_t`, and quantify again.

### 6.4 Infer the panel itself — "the data is its own panel"

A special case of 6.3 worth separating because it is how capture-based CNV callers already work without a
design: CNVkit's `guess_baits` infers the baited regions from coverage. Call probe footprints from pooled
coverage (RNA and gDNA), recognise junction probes as paired half-footprints across an intron with enriched
spliced reads, and compute `E_t` from the inferred map with the fitted physics. A further twist: vendor panels
are tilings — probes every `s` bases along each target transcript's spliced sequence — so a target's footprints
fix its tiling phase, and the phase predicts the probes that land on junctions of isoforms too rare to show a
footprint of their own. Speculative until a real panel is looked at.

### 6.5 Learn the kit across a cohort — the product path

Capture is a property of the kit, identical in every library captured with it; isoform usage is a property of
the sample. Across tens to hundreds of libraries the stable, per-feature component of the coverage is capture —
the "panel of normals" that exome and panel CNV callers (GATK gCNV, CNVkit) build for exactly this problem.
Product shape: `rigel panel learn --bams … → a panel profile`, then `rigel quant --panel-profile`. It fits the
MCTP cfRNA use (one kit, many libraries), it is exact where one library is not (a junction's capture is
separable from its usage once usage varies), and pooling gDNA across the cohort turns a sparse witness into a
deep one. How it fails: few libraries per kit; batch effects between capture lots; a kit revision. A profile is
learned data, not a design file, so it needs its own key (the kit, the annotation) and its own staleness rule.

### 6.6 Make the EM robust to what capture does — the safety net

The damage is phantom isoforms: a silent isoform that absorbs coverage a wrong length makes look like excess.
A sparsity prior on isoform support — already the next candidate for `ISSUES: per-transcript-prior-lane` —
resists that under any length error. Headroom, measured with true weights (a capability proof that hands over
the true support): the error at `g05` falls from 10.5 % to 3.3 % with the shipped lengths, and from 3.6 % to 1.6 %
with true ones. This does not fix capture; it makes the answer less sensitive to it, and it serves every stratum.

### 6.7 Say what cannot be known

Some splits will stay uncertain whatever is built. Quantify each locus under the range of plausible capture
fields (plain, the witness's, the learned one) and publish the spread per isoform: an isoform whose count moves
several-fold across them is capture-sensitive, and a user should know. Gene-level counts are robust under capture
(1.9 % at `g05`) and are the honest primary output where the spread is wide. `count_unambig` beside `count`
already tells a user which isoform calls rest on shared fragments.

### 6.8 When a design IS available

An optional probe BED makes 6.4's map exact; the physics parameters are still fitted. Cheap to support once 6.3
exists, and the one panel whose design you do hold is the calibration set for 6.1.

### Two ideas considered and set aside

* **A fragment-length witness of capture.** Overlap-weighted binding favours longer fragments near a probe, so the
  local fragment-length distribution carries capture strength independently of abundance. It would feed the
  opportunity model, not the retired composition channel — but it rests on the binding being linear in overlap,
  which the physics question leaves open, and it is weak where binding saturates. (In the simulator, which binds
  linearly, capture moves the ladder's RNA fragment-length mean from 216.75 to 229.09 bp.) Parked behind 6.1.
* **Spike-ins.** Probed spike-ins of known abundance calibrate the library-level capture parameters directly,
  if a kit carries them; they say nothing about where a sample's probes sit.

## 7. A recommended order

1. **Physics and panel design** (6.1) — your read of the chemistry; a parameterised half-match rule applied to
   cDNA and gDNA alike; and, on a real panel, how often probes sit on isoform-discriminating junctions. This is
   now the deciding step: under symmetric binding the shipped ruler already meets the target at `g05` on the
   ladder, and the ladder's 10.5 % is the simulator's asymmetry.
2. **If binding is symmetric** (my expectation): make it the simulator's default, re-measure the ladder, and move
   stranded × capture-ON's remaining error where it now sits — the EM's gDNA split at `g50`
   (`ISSUES: em-overturns-the-calibrated-gdna-split`). Junction-targeted designs stay a declared risk, with the
   capture field (6.3) as the insurance, prototyped on the test chromosome's twin.
3. **If the asymmetry is real**: the capture field (6.3, prototyped first as the two-pass form) becomes the main
   build, judged by within-gene spread and against `oracle_ruler`, with RNA, gDNA and nascent lengths from the
   one field; the guard of 6.2 ships meanwhile.
4. **In either case**: the sparsity prior (6.6), already the per-transcript lane's next build and a protection for
   every stratum; the cohort panel profile (6.5) as the product path for kits without a design file; capture
   sensitivity in the output (6.7).

## 8. Questions only the owner can answer

1. The chemistry: how strongly does a 50–60 base partial match to a 120-mer bait capture, relative to a full
   match — and is there any reason gDNA and cDNA would differ?
2. Is there one panel whose design you hold, and a captured library on it with gDNA? It calibrates everything.
3. Is a cohort (one kit, many libraries) the typical deployment? How many libraries per kit?
4. Would you accept a learned panel profile as an optional input, given the ruling that Rigel reads no panel?
5. Under capture, is gene-level the acceptable primary output where the isoform split is capture-sensitive?

## 9. Where the evidence is

`~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/`: `all_scenarios.txt` (the fractional panel),
`stranded_on_arms.txt` and `tables/` (every ruler and prior arm, per transcript), `testchr/` (benign against
junction), `physics/` (the symmetric-physics twin), `junction_witness_feasibility.txt`, `confusion_*`
(per-fragment truth against assignment at `g50`), and the scratch runner (`dissect_run.py`, `dissect_analyze.py`,
`confusion.py`). The ladder variant under symmetric physics: `~/Downloads/rigel_runs/suite/ladder_split1/`.
