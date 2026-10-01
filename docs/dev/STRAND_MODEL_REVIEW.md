# Strand model robustness: overdispersion, strand specificity and splice artifacts

*A brief for a review panel, 2026-09-30. Rigel `main` at 35e254ef. Nothing here is implemented; the measurements
come from read-only prototypes run on per-object tables extracted from real and simulated libraries (the appendix
lists the paths).*

## 0. What we are asking for

Rigel ships its 0.8.0 release soon, and one piece of the strand model is not robust on real data. We want your
ideas on two problems, the first primary:

1. **Estimate the strand overdispersion robustly**, so that calibration (the stage that separates RNA from
   genomic-DNA contamination) stays accurate on real libraries. Two kinds of contamination inflate today's
   estimates. Splice-alignment artifacts inflate the RNA side. Antisense RNA inflates the gDNA side.
2. **Begin designing a system that identifies and prunes splice artifacts.** Spliced reads whose strand
   contradicts the library protocol are the most obvious artifacts. They are also what corrupts the strand
   model's training set.

We think the problem is very solvable. **A simple, elegant design wins over a clever one.** Section 6 lists the
candidates we have considered, with the evidence for each. Section 10 asks specific questions, and section 11
gives the response format we would like.

---

## 1. Background

### 1.1 The tool

Rigel quantifies transcripts from RNA-seq while separating RNA from genomic-DNA (gDNA) contamination.
1. A single-pass C++ BAM scan tallies fragments into genomic objects: exon pieces, introns, and the boundaries
   between them.
2. **Calibration** deconvolves each object into gDNA and RNA.
3. A per-locus EM assigns the RNA to transcripts, using calibration's result as its prior.

There are exactly three populations of fragments:
- **gDNA**: double-stranded, so its fragments land on either strand with probability ½;
- **RNA from the + strand**;
- **RNA from the − strand**.

A spliced fragment is certified RNA. An unspliced fragment could be any of the three, and separating them is the
whole problem. The populations allowed at an object come from the annotation: gDNA always, plus RNA on each
strand that has a transcript there.

### 1.2 The strand channel

In a stranded library, RNA reads land on a fixed orientation relative to their transcript, with high fidelity.
We write **κ** for the fraction of RNA reads in the "sense" column under the protocol's convention:
- dUTP / R1-antisense protocols give κ ≈ 0.002 (99.8 % fidelity);
- unstranded libraries give κ = ½.

gDNA is ½ on every library. **On a stranded library, the strand split at an object is the strongest evidence
of how much gDNA it holds.** Calibration's per-object strand log-likelihood is (`docs/EQUATIONS.md` §5.1):

```
p    = ½·f_g + κ·(1 − f_g)                       f_g = the object's gDNA fraction, N = its fragments
var  = N·p(1−p) + (N·f_g)²·¼·od_g + (N·(1−f_g))²·κ(1−κ)·od_r
loglik = −½·(sense − N·p)² / var − ½·log var
```

**od** (overdispersion) is the intra-class correlation ρ of a Beta-Binomial: each object's true strand rate
varies around its population mean, so the split is wider than Binomial. `od_g` belongs to gDNA and `od_r` to
RNA. Two consequences matter:
- **od caps the evidence.** An object of N fragments carries the strand evidence of only `N_eff = N/(1 + (N−1)·od)`
  coin flips. At od = 0.2, a 1,523-fragment object is worth five. The information about `f_g` saturates at
  `(½−κ)²/(p(1−p)·od)`: **dispersion caps it, not depth.**
- **od only matters on stranded libraries.** At κ = ½ the strand term is flat whatever the od.

### 1.3 Why a wrong od is expensive

Measured on our simulated panels, by forcing ("injecting") a value into calibration. The error is calibration's
summed error against an oracle calibration:

| injected od | true od | where | calibration error |
|---|---|---|---|
| 0.2 | 0 | 16-condition benchmark ladder, stranded × capture ON | **+200 %** (436k → 1,308k fragments) |
| 0.2 | 0 | ladder, stranded × capture OFF | +18 % |
| 0.2 | 0 | small test panel, stranded × capture ON | +255 % |
| 0 | 0.05 (gDNA only) | odg05 panel, stranded × capture OFF | +20 % (transcripts +1 %) |

- **Too high** mutes the strand channel. Calibration then falls back on read density, which hybrid capture
  distorts.
- **Too low** makes the channel overconfident.

The costs are lopsided: an od that is too high is far worse than one that is too low.

---

## 2. What ships today, and what is wrong with it

`src/rigel/calibration/calibrate.py:_fit_strand`, `strand_balance.py` and `gdna_strand.py`.

**κ.** The pooled sense fraction over uniquely mapped spliced fragments at annotated junctions, one observation per
fragment (at its leftmost annotated junction), as a Beta(1,1) posterior mean. Nothing removes artifacts except the
index-time junction blacklist (§8.1).

**The RNA od.** A pooled method-of-moments fit over the per-junction strand table
(`SJStrandTable`: per junction, the count of reads agreeing and disagreeing with the motif strand), at mean κ.

**The gDNA od.** The "away-half" moment over gDNA seeds.
- **Seeds:** intron regions, exon|intron boundaries and gene edges of single-strand genes whose counts cannot hold
  unspliced mature RNA.
- **How it works:** each seed is oriented so that RNA of its own gene pulls the split toward κ, and only seeds
  that fall on the far side are kept.
- **The argument:** this is unbiased under any amount of *sense* RNA. But antisense RNA pushes seeds onto the far
  side, which is exactly what inflates it.

**The reconcile.** The two values are blended by "information". The inputs are mismatched:
- gDNA's value comes from one estimator and its precision from another, which is computed and then discarded;
- RNA's precision is computed as if its od were exactly 0. That credits it with up to **~700,000× too much
  evidence** on a real library (7.5×10⁷ against an honest ~10²).

So RNA's number always wins:
- on a simulated panel with gDNA od 0.05, the tool uses ≤ 0.0007 for gDNA;
- within a single 1 % subsample of a real library, the gDNA and RNA values sit up to 0.19 apart.

**Two constants remain:**
- every value is clipped at 0.2, the Beta(2,2) "ceiling", labelled a clamp until it can be retired;
- the no-evidence fallback is 0.2 (it has never fired). The owner has ruled it becomes 0.

### Measured instability on real libraries

The shipped values, and the joint fit we prototyped (§6, option 3):

| library | depth | shipped (RNA / gDNA) | joint fit | note |
|---|---|---|---|---|
| VCaP, RNA-only half | full | 0.047 / 0.047 | **0.2 (clip)** | gDNA seeds are antisense RNA: raw moment 0.42–0.86 |
| VCaP, RNA-only half | 1 % (3 draws) | 0.020–0.065 / **0.121–0.200** | 0.2 | |
| MO_3021 | full | 0.008 (gDNA) | 0.121 | gDNA side alone 0.133 |
| MO_3021 | 1 % (3 draws) | gDNA 0.12–0.20 | 0.2 | |
| LBX0588 (~90 % gDNA) | full | 0.055 / 0.017 | 0.017 | clean gDNA seeds |
| LBX0588 | 10 % (3 draws) | gDNA 0–0.019 | 0.011–0.019 | gDNA total varies 3.4 % across draws (shipped) against a 0.7–1.5 % floor |

---

## 3. What real data shows

All the tables below come from per-object censuses extracted from real scans.

### 3.1 Junctions (the RNA side): splice artifacts inflate κ

The wrong-strand fraction of junction reads, by junction depth:

| reads at the junction | LBX0588 (~90 % gDNA) | MO_3021 | LBX0190 | VCaP RNA, 10 % |
|---|---|---|---|---|
| 1 | 12.0 % | 0.21 % | 0.28 % | 0.06 % |
| 2–3 | 13.8 % | 0.17 % | 0.34 % | 0.01 % |
| 4–7 | 6.2 % | 0.22 % | 0.18 % | 0.01 % |
| 8–31 | 1.15 % | 0.27 % | 0.24 % | 0.00 % |
| 32–127 | 0.52 % | 0.28 % | 0.15 % | 0.00 % |
| ≥ 128 | 0.24 % | 0.23 % | 0.29 % | 0.00 % |

- **On the gDNA-heavy library, shallow junctions are largely artifacts.** gDNA misaligned as spliced has a
  strand split near ½. As depth rises, the rate falls to the true κ.
- **RNA-rich libraries are flat with depth.**
- So the pooled κ (0.064 on LBX0588) is about 20× the deep-junction value.

An older case shows the same thing:
- LBX0077 (96 % gDNA, 2026-08) had a strand specificity of 0.786, against 0.991–0.999 for its 20 same-flowcell
  siblings.
- Two independent estimates put **42–43 % of its spliced "RNA"** in the artifact category.
- κ and the RNA fragment-length law are both trained on that population.

### 3.2 Seeds (the gDNA side): antisense RNA inflates the gDNA od

- **LBX0588:** away fractions centre on 0.49 and are symmetric at every depth. These seeds are gDNA.
- **MO_3021, VCaP RNA and LBX0190:** away fractions pile up at 0 (host-strand RNA) and at 1 (antisense RNA). The
  away-half estimator keeps exactly the pile at 1.

### 3.3 A mixture prototype: every object belongs to one of the three populations

**The model.**
- Each junction or seed takes its strand count from one population, with a mean fixed by physics: gDNA ½,
  host-strand RNA κ, opposite-strand RNA 1−κ.
- A junction in the gDNA class **is** a splice artifact. A junction in the opposite-strand class is
  *reverse-stranded*.
- Contamination enters as a different **mean**, never as extra variance.
- It is fit by maximum likelihood with no tuning constant.

**Junctions alone** (artifact and reversed classes fit, κ re-estimated). "Evidence" is the log-likelihood gain of
the best od over od = 0, with κ and the class weights re-maximised at each od:

| library | shipped κ | κ, artifacts out | artifact share | reversed share | best od | evidence |
|---|---|---|---|---|---|---|
| LBX0588, full depth | 0.064 | **0.0029** | 22 % | 0.2 % | 0.03 | 16.3 |
| LBX0588, 10 % (3 draws) | 0.060–0.070 | 0.0024–0.0064 | 1–14 % | 3–8 % | 0 / 0.05 / 0.2 | 0.0 / 0.4 / 0.3 |
| MO_3021, full / 10 % / 1 % | 0.0023 | 0.0023 / 0.0021 / 0.0012 | 0 | 0 | 0.02 / 0.02 / 0 | 215 / 5.2 / 0.0 |
| VCaP RNA, 10 % / 1 % | 0.0001 | ≈ 0 | 0 | ≈ 0 | 0.03 / 0.999 | 1.1 / 0.0 |
| LBX0190 | 0.0025 | 0.0024 | 0.1 % | ≈ 0 | 0.003 | 27.3 |
| every simulated row (true od 0) | ≈ 0.010 | exact to 4 digits | 0 | 0 | 0 | 0.0 |

**κ hardly depends on the od.** On LBX0588 it stays at 0.0024–0.0030 for any od from 0 to 0.03.

**Junctions and seeds together**, with one od for each population:
- **LBX0588:** gDNA 0.015, RNA 0.024. One shared od fits both (likelihood-ratio statistic 0.6).
- **MO_3021:** gDNA 0.068, RNA 0.014. The data reject one shared value (statistic 86), but the gDNA seeds there
  are partly antisense RNA (§4).

### 3.4 Where the prototype fails

1. **Walking from od = 0 finds the wrong answer.** Iterating "label, then estimate" from binomial on LBX0588's
   junctions settled on κ 0.11 and od 0.52 with no artifacts. That is 27 log-likelihood units worse than the
   global fit (κ 0.0029, 22 % artifacts). The artifact class has to be found by a global fit.
2. **Relabelling explains away every bit of spread.** On LBX0588's seeds, the same walk labelled 26 % of seeds
   antisense RNA and 28 % host RNA, and returned od = 0 exactly, on a library that is about 90 % gDNA.
3. **Partly mixed objects read as od.**
   - The ladder's g05 capture-ON row has a true od of 0, yet its seeds read 0.021, and 0.134 when the RNA inside
     the seeds is forced to be clean.
   - odg05's g05 row (gDNA od 0.05) reads 0.073.
4. **The od is undetermined at κ ≈ 0.** On VCaP at 1 %, the junction likelihood is flat to within 0.005 units
   anywhere from od = 0 to od = 0.999. The maximum is arbitrary (0.999), driven by 2 reversed junctions out of
   44,241.
5. **Low depth.** On LBX0588 at 10 % (mostly 2–3 fragments per seed), two fits of the same seed model from
   different starting weights land far apart:
   - one labels 62–81 % of seeds gDNA, with od 0.003–0.015;
   - the other labels 4–6 % gDNA, with od ≤ 0.003;
   - at full depth, 94 % are gDNA.

   At this depth the labellings are barely distinguishable.

---

## 4. The identifiability problem at the centre

**Composition heterogeneity is indistinguishable from overdispersion.** Suppose a seed's gDNA fraction f varies
across seeds. Its strand rate is then a mixture, `½·f + κ·(1−f)` with sense RNA, or `½·f + (1−κ)·(1−f)` with
antisense RNA.
- If seeds carry sense RNA on some and antisense RNA on others, the rates spread on both sides of ½, exactly as a
  Beta around ½ does.
- From the strand count alone, "gDNA with a noisy strand rate" and "gDNA plus a little RNA of either strand" have
  the same likelihood.
- **Whole-object contamination is identifiable,** because it sits at a different mean (a seed that is all
  antisense RNA, a junction that is all misaligned DNA). **Partial mixing is not.**

**Junctions are close to whole-object.** An unexpressed junction hit by misaligned DNA is pure artifact. An
expressed junction's artifact reads are a small share of a deep junction. Intron seeds are often partly mixed:
gDNA plus some unspliced host RNA.

**The ICC parameterisation is degenerate at extreme κ.** With mean κ ≈ 10⁻⁴, a Beta-Binomial with ρ → 1 is a
two-point mixture: a junction reads either all sense or all antisense. That is the same likelihood as "a few fully
reversed junctions", so a handful of reads can move the RNA od anywhere in [0, 1). Yet the solve's RNA width
`κ(1−κ)·od_r·N²` does depend on it. At N = 1,000 and κ = 10⁻⁴, od_r = 1 lets RNA explain about ±10 wrong-strand
reads, against ±0.3 at od_r = 0.

**Calibration's own labels cannot break this for antisense.** At an object where the annotation admits RNA on only
one strand, the solve reads antisense RNA as gDNA. It would certify exactly the contaminated seeds as pure.

---

## 5. Rulings and constraints (binding)

- **One shared od** for gDNA and RNA (owner ruling). The physics need not force equality, but a difference between
  two noisy estimates biases the strand model, and robustness wins.
  - On the only simulated panel with a planted difference (gDNA 0.05, RNA 0), a shared value was within noise or
    better on both headline numbers.
  - A two-value model cannot be validated: the simulator cannot plant an RNA od, and real data has no truth.
- **The fallback is 0 (binomial).** With no evidence, the od is 0. "Start at binomial and accumulate evidence"
  (owner).
- **No estimate may be able to sabotage the tool,** but no arbitrary hard cap either. The 0.2 clip stays,
  labelled a clamp, until a design retires it.
- **No magic numbers.** Every constant, threshold or tunable is derived, or explicitly ruled by the owner with
  its reason. Numeric tolerances that move no result are exempt.
- **Robustness first, then accuracy.** A design must survive these before any accuracy reading counts:
  - zero gDNA;
  - zero RNA;
  - low counts (1 % and 10 % subsamples).
- **Real data is a test input, never a design input.** Design on derivations and simulations; use real libraries
  to check stability and to falsify.
- **Fragment-length composition is retired until after 0.8.0** and may not be proposed as a channel.
- **Scope.** Three strata are in scope: unstranded × capture OFF, stranded × capture OFF, stranded × capture ON.
  Unstranded × capture ON is reported but not a target. Every score is read per stratum, never pooled.
- **Every change** goes through a falsification test first, verified failing, and an A/B against what ships under
  pinned conditions.

---

## 6. Candidate designs for the od

**Option 1: od from genuine junctions only, starting at 0.** This is the current lead.
- **How:**
  - Fit the three-class junction mixture globally. That removes artifacts and reversed junctions.
  - Fix κ from the genuine class.
  - Take the od from the genuine class: 0 unless the likelihood gain clears an evidence bar.
  - Apply that one value to gDNA and RNA. gDNA seeds become a consistency check only.
- **For it:**
  - Junctions are certified RNA once artifacts are out, so they are the cleanest pure population.
  - It never touches the antisense-contaminated seeds.
  - It is exactly 0 on every simulated row.
  - The real-data instabilities fall away: the three LBX0588 10 % draws and VCaP 1 % all carry under 0.4 units
    of evidence, so they stay at 0.
  - Where evidence exists, the value holds across depth (MO_3021 reads 0.02 at full depth and at 10 %).
  - On unstranded libraries no artifact detection is needed: artifacts and genuine junctions both sit at ½, so
    neither κ nor the od is biased.
- **Against:**
  - If gDNA's true od exceeds RNA's, the tool uses RNA's value. On odg05 that costs about 20 % of stranded
    capture-OFF calibration error.
  - The one real library with clean gDNA seeds (LBX0588) shows no such gap: seeds 0.008–0.015, junctions 0.03.
  - It needs one constant, the evidence bar (question 3).

**Option 2: one joint three-class mixture over junctions and seeds, with a shared od.**
- **For it:** it uses all the evidence and handles whole-seed antisense contamination.
- **Against:** the failures in §3.4: relabelling, partial mixing (+0.02 to +0.13 where the truth is 0), and swings
  at low depth.

**Option 3: the joint influence-weighted moment** (prototyped and measured 2026-09-27). One estimating equation
sums gDNA seeds and junctions, with each object weighted by its inverse variance.
- **For it:**
  - On odg05: calibration −10 % / −19 % and transcripts −3 % / −7 % (stranded OFF / ON).
  - Steadiest on LBX0588 at 10 %.
- **Against:**
  - On antisense-contaminated real libraries it sits at the 0.2 clip (VCaP RNA, MO_3021 subsamples).
  - At low gDNA it has two roots, so it needs a root rule.
  - Artifacts inflate its junction side.

**Option 4: robust M-estimation of a single Beta-Binomial.** For example, minimum density-power-divergence, or
Huber-type weights on standardised residuals.
- **For it:** bounded influence; contaminated objects are down-weighted automatically.
- **Against:**
  - It needs a tuning constant (the divergence α or the Huber cutoff).
  - It must stay Fisher-consistent at the Beta-Binomial, where heavy tails are the very thing being measured.
  - It does nothing for partial mixing.

**Option 5: two passes.** Calibrate at od = 0, estimate the od from objects the solve calls pure, then calibrate
once more.
- **For it:** composition comes from the full deconvolution, density included.
- **Against:**
  - It is circular for antisense: the solve calls antisense-contaminated objects gDNA (§4).
  - Capture makes density-based composition unreliable.
  - It doubles calibration's cost.

**Option 6: binomial everywhere (od ≡ 0).**
- **For it:** the simplest possible; it cannot be sabotaged; it is what Option 1 returns whenever evidence is weak.
- **Against:** where a real od exists, the channel is overconfident. That costs +20 % on odg05, and on MO_3021,
  where junctions show an od of 0.02 with overwhelming evidence (215 units).

**Option 7: change what od means.**
- **Variants:**
  - Estimate the quantity the solve consumes (the strand-count variance at an object of size N) directly, not
    the ICC.
  - Or use a parameterisation that stays identifiable at extreme κ, such as a logit-normal rate or a dispersion
    on the log-odds scale.
- **For it:** it removes the κ → 0 degeneracy (§4).
- **Against:** it means changing the native strand term and its derivation.

We have not found a candidate that measures a gDNA-only od robustly under antisense contamination without some
external composition information. Section 10 asks whether one exists.

---

## 7. Strand specificity κ

**What we found.** Splice artifacts bias κ on gDNA-heavy stranded libraries (§3.1), and κ enters every
object's strand likelihood. The three-class junction fit recovers it, exact on every simulated row, and on real
libraries:
- stable across depth where artifacts are absent;
- 0.002–0.006 where the pooled value is 0.06–0.07.

**Alternatives we can see:**
- weight each junction by its posterior probability of being genuine;
- estimate κ from deep junctions only (a depth threshold, which is a constant);
- extrapolate the depth trend.

**Not yet measured:** how much the corrected κ moves calibration and the transcript table on a real library. Our
simulator writes no splice artifacts, so the simulated panels cannot see this defect at all.

---

## 8. Splice artifacts: the secondary goal

### 8.1 What happens today

- **The blacklist is built at index time.** `alignable` simulates genomic reads, aligns them, and lists junctions
  the aligner wrote at least twice, with the longest anchor seen.
- **The reject rule runs per alignment record during the scan.** It matches `(ref, start, end)` ignoring strand
  and rejects when either anchor is at or below the stored maximum. This happens before mates are joined, so a
  junction rejected on one mate and kept on the other survives.
- **A fragment whose every junction is rejected** is labelled ARTIFACT:
  - calibration holds it out;
  - the EM lets gDNA explain it, but at a footprint that still spans the rejected intron.
- **A fragment with one junction rejected and one kept** stays spliced.
- **The strand model and the RNA fragment-length law** train only on uniquely mapped fragments at annotated
  junctions. So the contamination in §3.1 is uniquely mapped artifacts at annotated junctions.

### 8.2 Known failures (ISSUES: splicing-artifacts)

- **It fails both ways.**
  - On the VCaP RNA-only half, the blacklist rejected 598,010 records where a pure-gDNA library's rate predicts
    about 20,000.
  - On gDNA-heavy libraries, artifacts escape (LBX0077: 42–43 %).
- **The largest artifact class is decided by the reference sequence.** At annotated junctions with a 1–10 bp short
  side, 96 % of DNA-library artifacts have zero reference mismatches between their spliced and unspliced
  placements.
  - The aligner's junction-database bonus breaks the tie toward splicing.
  - Rigel reads NM only, never MD, the sequence or the reference.
- **Long-anchor artifacts are probably of remote origin** (retrocopies, paralogs). That hypothesis is not yet
  measured.
- **The blacklist builder applies its minimum count per read length before aggregating,** and it stores no count.

### 8.3 The strand signal: two obvious artifact types

On a stranded library:
- **Unstranded junction:** reads split near ½. This is misaligned gDNA: the gDNA class in §3.3.
- **Reverse-stranded junction:** reads consistently opposite the junction's motif strand. A genuine read would
  need the transcript to splice at the reverse-complement motif, so this is almost certainly an alignment
  artifact. It is the opposite-strand class in §3.3.

The same per-junction fit that gives κ and the od gives each junction a posterior probability for each class.
**Limits:**
- **Shallow artifacts mostly look genuine.** An n-read artifact junction has every read on the right strand with
  probability 2⁻ⁿ, so it shows no strand signal at all. With LBX0588's fit, P(artifact) is:
  - 0.97 for a 1-read junction whose read is on the wrong strand;
  - 0.12 for a 1-read junction whose read is on the right strand;
  - 0.07 for a 2-read junction with both reads on the right strand.
- **κ and the od need only the library-level share and per-junction weights. Pruning needs per-junction
  decisions,** and for shallow junctions strand alone cannot make them.
- **On unstranded libraries strand sees nothing.** There, artifacts bias neither κ nor the od, but they still
  corrupt quantification, because a spliced fragment is certified RNA.

### 8.4 Other junction-level signals, available once junctions are aggregated

- the anchor-length distribution across the junction's reads;
- mismatches near the junction (needs MD or the reference, which Rigel does not read);
- annotated versus novel;
- the multimapping share;
- motif class;
- **depth relative to the local gDNA density.** The expected artifact count at a junction is roughly the local
  gDNA coverage times that junction's misalignment propensity. A junction whose depth that alone explains, in an
  otherwise unexpressed region, is suspect;
- the aligner's own junction table (STAR's `SJ.out.tab`).

### 8.5 The architectural problem

A junction can be judged only after the scan has aggregated it. But in one pass the scan already:
1. deposits calibration's tallies, spliced fragments included;
2. fills the per-junction strand table;
3. buffers the EM's fragments, which carry **no junction identity** (`src/rigel/buffer.py`: splice type,
   candidates, footprint, NM).

So pruning a junction after aggregation is not possible today without one of:
- **(a)** carrying a junction id on each buffered spliced fragment, and deferring the spliced fragments'
  calibration deposits until junction decisions are made;
- **(b)** a junction pre-pass before the main scan: read only the spliced records, or the aligner's junction table;
- **(c)** keeping decisions at index time only (the status quo, which fails both ways);
- **(d)** scan, decide, then re-scan.

**What a pruned fragment becomes.** The owner's leading candidate is to edit the alignment: delete the rejected
junction and keep the largest aligned block. The fragment is then unspliced at its true length in both calibration
and the EM.

### 8.6 How the two problems couple

**The owner's proposed procedure:**
1. identify and remove splicing artifacts;
2. estimate strand specificity and fix it;
3. estimate each object's gDNA/RNA mixture;
4. estimate the od from pure objects.

Iterated to self-consistency, steps 1, 2 and 4 on junctions are the mixture fit of §3.3. Step 3 on seeds is where
§4's identifiability problem lives.

**One aggregated junction table feeds four outputs:**
- the library's artifact share;
- per-junction artifact posteriors, for pruning;
- κ;
- the od (Option 1).

---

## 9. What we can test cheaply

- **Per-object censuses.**
  - Six real libraries at several depths, and simulated panels with known truth.
  - Junction strand counts, seed strand counts and κ.
  - Any estimator can be scored in seconds to minutes, with no pipeline run (appendix).
- **Simulated panels.**
  - True od 0 everywhere, except odg05 (gDNA 0.05, RNA 0).
  - The simulator writes no splice artifacts, no RNA od and no antisense RNA on gene seeds. It is therefore blind
    to the κ defect and to the gDNA-side contamination.
  - Building new simulated panels is expensive; designs testable on existing data are strongly preferred.
- **The VCaP mix.** Real alignments of an exome-DNA library mixed with a transcriptome library, with per-fragment
  truth from read names. It is our one real library with truth: real artifacts, and real gDNA against real RNA.
  It costs one whole-genome job per arm.
- **The release instruments.**
  - calibration against an oracle calibration, per stratum;
  - the transcript table against per-transcript truth;
  - gene-level and pool-level error;
  - all in pinned A/B pairs.

---

## 10. Questions for the panel

1. **Option 1.** Is "od from genuine junctions only, applied to both populations" sound? What breaks it? In
   particular:
   - a library with real gDNA strand noise but clean RNA;
   - a junction-poor library;
   - capture.
2. **The gDNA side.** Is there a simple estimator of a gDNA od that survives antisense contamination without
   external composition information? Or should gDNA seeds only check the junction value? If a check, what should
   happen when they disagree?
3. **The evidence bar.** How should "enough evidence for od > 0" be decided?
   - Candidates: a boundary likelihood-ratio test at a stated level (a 5 % bar is a gain of 1.35); an
     information criterion; a Bayes factor.
   - The codebase already decides whether the strand channel is live by a Bayes factor at equal prior odds, with
     no constant (`docs/EQUATIONS.md` §5.2b). Is there an equally constant-free rule for the od? What prior over
     the od would it need, and is that prior a hidden constant?
4. **Parameterisation.** Is the ICC the right od for RNA at κ ≈ 10⁻³ to 10⁻⁴, given the degeneracy in §4? Would
   a log-odds-scale dispersion, or estimating the solve's variance directly, be simpler and safer?
5. **κ.** Mixture-based κ, posterior-weighted κ, or something simpler? Does κ need its own uncertainty in the
   strand likelihood?
6. **Shallow artifacts.** For pruning, how should junctions too shallow for strand evidence be treated? Pool them
   by depth? Use the local gDNA density? Down-weight rather than prune?
7. **Pruning architecture.** Which of §8.5's options (a)–(d), or another? What should a pruned fragment become?
8. **Unstranded libraries.** Which non-strand signals (§8.4) are worth building first?
9. **One shared od.** Does the data (§3.3) argue against it?
10. **What have we missed?** Failure modes, simpler designs, or a reason this framing is wrong.

## 11. Response format

Please return, in this order:
1. **A verdict on each option in §6**, one or two sentences each: keep, fix (how), or drop (why).
2. **Your proposed design**, at most one page:
   - the estimator for κ and for the od;
   - the order of operations;
   - what happens with no evidence;
   - how the artifact posterior feeds pruning.
3. **Every constant your design introduces,** and how it is derived or why it is unavoidable.
4. **Failure modes,** each with a cheap test, preferably on the existing censuses (§9).
5. **Answers to the questions in §10,** as short as possible, and only where you have a view.

---

## Appendix: data and code

**Prototypes (read-only, numpy/scipy):** `~/Downloads/rigel_runs/prototypes/2026-09-30_robust_od/`

| file | contents |
|---|---|
| `look.py` | depth tables of junction wrong-strand fractions and seed away fractions |
| `mix.py` | the three-class Beta-Binomial mixture; a stable log-pmf gated against scipy |
| `em.py` | the owner's four steps as label/estimate iterations from od = 0 |
| `junction_evidence.py` | the global junction fit and the profile-likelihood evidence for od > 0 |
| `real3.txt` / `synth3.txt` | the mixture results |

**Censuses:** `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/scratch/`
- `census_<library>.npz` (real): per-junction `sj_sense`/`sj_anti`; per-seed `r_*` (introns) and `b_*`
  (boundaries) sense/total; `kappa`.
- `seeds/*.npz` (simulated, with known truth): `g_sense`/`g_total` (seeds), `r_sense`/`r_total` (junctions).

**The 2026-09-27 measurements of the shipped value, the joint fit, binomial and the ceiling:**
`~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/REPORT.md`. Earlier derivations:
`~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/` (`00_synthesis.md` first).

**Code:**

| file | what it holds |
|---|---|
| `src/rigel/calibration/calibrate.py` (`_fit_strand`) | the strand fit's order |
| `src/rigel/calibration/strand_balance.py` | κ |
| `src/rigel/calibration/gdna_strand.py` | both od fits and the reconcile |
| `src/rigel/strand_model.py` (`SJStrandTable`) | the per-junction table |
| `src/rigel/native/transfer_rows.h` | the native strand term |
| `src/rigel/native/bam_scanner.cpp` | the junction reject rule |
| `src/rigel/splice_blacklist.py` | the blacklist build |
| `src/rigel/buffer.py` | the EM's fragment record |

**Docs:**
- `docs/EQUATIONS.md` §5–§6: the strand likelihood and od. §6a–§6c describe a weighted fit that never shipped.
- `docs/ISSUES.md`: `strand-overdispersion-one-shared-value`, `the-overdispersion-design-refusals` (designs
  already refused, each with the number that killed it), `splicing-artifacts`.
- `docs/dev/GDNA_SPLICE_ARTIFACTS_*.md`: the artifact investigation.
