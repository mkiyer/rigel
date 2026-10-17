# Rigel 0.8.0 — reader results, plan stages 1–4 (2026-10-12)

Sandbox, not authoritative. Executes stages 1–4 of `RNA_RELEASE_0_8_0_PLAN.md`. No `src/` change, no
installation change, no commit. Everything ran from an isolated wheel built from a worktree at `624dfade` plus a
19-line reader-policy patch. Prototype, receipts, logs and outputs:
`$RIGEL_RUNS/prototypes/2026-10-12_reader_split/` (README there). Protocol for every arm: scan pinned to
one thread, fractional assignment, EM seed 0, other stages on their normal threads. Accuracy runs only, so there is
no speed claim.

## 1. The candidate

**Reader rule.** Each object's capture weight is read from its own strand evidence. At an object admitting both RNA
strands (`statics.free_pos & statics.free_neg`) the held DNA level also enters. The composition row and the held RNA
levels no longer enter the reader anywhere. In the prototype this is policy 3 ("split·DNA").

**What stays the same.** The count solve, the landscape, the median unit and the EM are unchanged. Calibration's
six count arrays are identical to main's in all **554** compared conditions. The calibration metric `P − O` is
therefore unchanged by construction.

**Why this rule.**
- At a single-strand object, two strand counts identify the gDNA amount for any κ ≠ ½.
- At a both-strand object they cannot: gDNA is confounded with equal sense and antisense RNA at every κ.
- Stage 1 found that the held DNA level is the one neighbour factor that resolves that confound.

## 2. Stage 1 — attribution and identity (calibration only, 11 conditions, 9 policies)

- **Identities, exact on every condition:**
  - the six count arrays are identical across policies;
  - policy 0 reproduces the archived main weights to 0.0;
  - policy 1 reproduces the archived own-evidence weights to 0.0;
  - the split's raw modes equal own evidence at single-strand slots and shipped at both-strand slots.
- **Zero-gDNA false weights (ladder g00 stranded OFF):** 696 under shipped, 647 of them single-strand.
  - Dropping the composition row everywhere leaves 132.
  - Dropping the DNA level leaves 582, and removes all 49 both-strand ones.
  - The RNA levels create none.
- **Captured both-strand exons (test g98/g50 stranded ON):**
  - A both-strand slot keeping only the DNA level keeps 100% of captured both-strand slots at g98. The
    composition row alone or the RNA levels alone keep 91%, the same as own evidence.
  - Against the simulator's own yield (host 197.3, antisense 90.9 for all three pairs), the DNA-level-only rule
    reads the g50 antisense factors 88.3 / 90.1 / 134.4. The full shipped factor set reads 85.2 / 22.5 / 24.2.

## 3. Release table — in-scope strata against the tag

Transcript % (gene %), ratio of sums; g00 kept apart. Tag = v0.7.1. Candidate = split·DNA. A × marks a stratum
where the candidate exceeds the tag.

| Panel · stratum | Tag | Main | Candidate |
|---|---:|---:|---:|
| Ladder str OFF contaminated / g00 | 3.95 / 3.70 | 3.37 / 13.61 | 2.62 / **3.81 ×** (genes 0.34 vs tag 0.47) |
| Ladder str ON contaminated / g00 | 9.91 / 20.22 | 6.05 / 15.10 | 6.24 / 6.61 |
| Ladder unstr OFF contaminated / g00 | 10.29 / 20.09 | 2.40 / 1.59 | 2.34 / 1.59 |
| RNA-short str OFF / ON / unstr OFF | 6.31 / 20.33 / 6.00 | 4.99 / 14.24 / 3.95 | 4.61 / 20.09 / 3.94 |
| RNA-long str OFF / ON / unstr OFF | 5.23 / 7.38 / 5.29 | 3.52 / 6.09 / 2.55 | 3.46 / 6.33 / 2.21 |
| Test str OFF / ON / unstr OFF, contaminated | 8.20 / 8.58 / 9.25 | 9.35 / 8.17 / 9.18 | 8.03 / 8.27 / 8.64 |
| Test g00 str OFF / str ON / unstr OFF | 11.35 / 14.22 / 9.38 | 7.66 / 8.22 / 5.62 | 7.70 / 6.81 / 5.61 |

The ladder g00 stranded-OFF excess is 0.11 point (about 2,300 fragments of 2.07 M). The ladder has no replicate
seeds. On the test chromosome, three g00 seeds vary by 0.5–1.0 point at a tenth of the depth, which suggests a
replicate spread of about 0.2–0.3 point at ladder depth; that is an estimate, not a measurement.

**Per rung against the tag:** every in-scope ladder and gap rung is at or below the tag, except ladder g00
stranded OFF (+0.11).

**Rungs worse than main:**
- RNA-short stranded ON: +5.85.
- Ladder g05 stranded ON: +0.33.
- RNA-long stranded ON: +0.24.
- Ladder g98 stranded OFF: +0.08.

**Test g98 stranded ON:** the candidate reads 53.16 % (genes 22.12), against tag 50.71, main 52.69 and own
evidence 72.21. Its antisense overcalls are gone: the ten largest gene overcalls sum to 2,353, against main's
2,262 and own evidence's 6,685.

## 4. Robustness panels (test chromosome, 259 conditions, tag run fresh where the October 5 sweep lacked it)

| Class | Candidate against tag | Candidate against main |
|---|---|---|
| Capture-OFF rows, any κ, fragment-length law, depth, gDNA overdispersion | at or below the tag everywhere except odg05 κ = 0.70 OFF (+0.3) | better on nearly every row |
| Stranded ON, equal or longer RNA, probes (junction, sparse, single) | below the tag | from 1.2 worse to 2.2 better (sparse +1.2, bind01 +0.6; better on junction, single, d100, fle200) |
| Stranded ON, gDNA longer than RNA (short exons) | gdna-long 12.08 vs 12.83: below | worse: +2.3 (and RNA-short +5.9) |
| **Stranded ON, scarce gDNA under capture** | **lowg g001 39.7 vs 7.7; d10 g01 38.3 vs 19.2** | equal to main (38.7, 38.1) |
| Stranded ON, overdispersed gDNA or intermediate capture | odg05 17.4 vs 14.5; bind1 10.3 vs 9.3 | about equal to main (17.0, 10.8) |
| Weak strand κ = 0.70 ON | mixed: bind01 25.1 vs 22.8, odg05 27.6 vs 16.6, fle200 12.8 vs 11.9; below elsewhere | worse: +1.5 to +13 |
| Unstranded × capture-ON (deferred) | mixed | worse on test-chromosome panels, better on the ladder (52.7 vs 63.0) and RNA-long (14.8 vs 34.1) |

## 5. Real libraries

| Library | Arm | gDNA fraction | Weights above e | Note |
|---|---|---:|---:|---|
| VCaP mix (origin ≈ 0.2518) | main / candidate | 0.2528 / 0.2533 | — / 650,663 | exome-captured DNA |
| VCaP transcriptome reads only (new real near-gDNA-free control) | main / candidate | 0.0078 / 0.0077 | 108,398 / 49,283 | in the mix's captured footprint: 53 % / 81 % |
| cfRNA-A · cfRNA-C · cfRNA-B (no truth) | main / candidate | 0.0997 / 0.0990 · 0.1653 / 0.1651 · 0.9373 / 0.9373 | 10,349 / 5,115 · 45,269 / 28,723 · 67,224 / 48,214 | pools unchanged |

The transcriptome-only library is itself captured: its extreme weights concentrate in the exome footprint, which
covers 25 % of objects. The weights main adds over the candidate (58,783) lie 70 % outside that footprint, so they
are mostly false capture. All runs exited 0, at 6–11 GB.

## 6. Failure classes, each with its mechanism

1. **Ladder g00 stranded OFF, +0.11 against the tag.** At both-strand objects the held DNA level imports a
   neighbour's strand-noise "gDNA" as a lower bound. These are 49 weights, none among the 93 most extreme. Removing
   them restores the antisense failure.
2. **Short exons under long gDNA (stranded ON, RNA-short and gdna-long).** gDNA cannot be contained in a piece shorter
   than its fragments, so own evidence is empty there. The composition row and RNA levels that used to substitute
   are gone. The candidate stays at or below the tag.
3. **Scarce gDNA under capture.** A few hundred to a few thousand gDNA fragments in the library: own evidence reads
   capture from noise, and isoform splits suffer while gene error stays under 1 %. Main fails identically; the
   tag's detector fell back to flat weights. This gap has existed since 2026-10-10 and the candidate does not
   change it.
4. **Weak strand (κ = 0.70) and unstranded captured rows.** Own evidence weakens as κ → ½, and the neighbour
   information that carried capture there is removed at single-strand objects.

**Common root of classes 1–4.** At zero or scarce gDNA, own evidence is flat at low density and cannot contradict
any neighbour's level. A neighbour's level is false when its "gDNA" is strand noise, made firm by the binomial
strand model (od = 0 by ruling). So no recipient-side admission rule separates true from false neighbour levels.
The fix, if one is wanted, is upstream: in how a level is formed, or in shrinking a contrast toward flat in
proportion to the gDNA behind it.

## 7. The decision this leaves to the owner

Accept the candidate as 0.8.0's reader, with classes 1–4 listed in the manual as limitations and as ISSUES entries,
then proceed to the port (stage 5) and release (stage 6)? Or require one of classes 2–4 to be fixed first, choosing
the mechanism from its cause above?

In favour of accepting:
- It removes the release-blocking flat-field failure: ladder g00 stranded OFF 13.61 → 3.81 %, and on a real library
  it removes about 59,000 mostly false extreme weights.
- It improves almost every capture-OFF row across ten panels.
- It keeps the antisense safety that own evidence alone loses.
- It leaves calibration untouched.
- It meets the tag on every in-scope stratum within 0.11 point.

Against: class 3 is a real-data-relevant regime, captured RNA with little gDNA, where both main and the candidate
are far worse than the tag on isoforms.
