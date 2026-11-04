# The capture reader: working memory

This file was kept at the 2026-11-03 reset as the one record of the capture-reader work from 2026-10-29 to 2026-11-03.
It is not a design. It is where the owner and Claude start developing the reader together.

**Conventions**
- Numbers are from the ladder's STAR rung on the source at 81fc67b6 (identical to ec996bba), unless marked.
- **c1** means every capture efficiency is set to 1.
- **×a / ×b** is a ratio of transcript / gene error, against the shipped reader unless marked "own c1".
- The prototypes are in `$RIGEL_RUNS/prototypes/<campaign>/`. The retired design documents and reviews are in
  `$RIGEL_RUNS/archive/cleanup_2026-11-03/docs/dev/`, mainly `CAPTURE_READER_DESIGN.md`.

## 1. What ships

- **The reader** (`calibration/capture_efficiency.py`) reads each object's calibrated gDNA booking, `k = f_g·M`, as a
  Poisson count on its gDNA support `S`.
  - Its prior is the gDNA landscape split into two groups. The largest basin collapses onto its peak, the *bulk*; the
    rest is the population's own density.
  - Each enriched basin's peak caps the reading.
  - The efficiency is `exp(E[level] − bulk)`.
- **The landscape** (`calibration/landscape.py`, fitted in `calibrate._fit_gdna_hyperprior`) is trained on the booked
  gDNA of locked regions, located single-strand regions and zero-count regions. It is refitted three times; each refit
  re-solves every object under the previous landscape.
- **The efficiencies** scale the EM's effective lengths for gDNA and for every transcript
  (`calibration/capture_eff_length.py`, `calibration/priors.py`).

## 2. Why it fails

**1. A posterior is read back as data, three times. Measured.**
- The three places:
  - the landscape is trained on `f_g·M`;
  - the E-step starts from the previous landscape;
  - the reader takes the booking `k` as a Poisson count.
- At one g00 boundary, ψ prefers 20 gDNA fragments over none by e^4.6, and the reader turns that into e^250.

**2. The bulk is a resolution wall, not a measured density, when gDNA is scarce. Measured.**
- At g00 the bulk is 5.0e-7 per bp, against the grid's deepest wall of 4.3e-7. At that bulk, 3.6 booked fragments read
  as 10⁴×.
- At sparse capture the bulk is 7.5e-7 against a true 2.9e-6 to 1.6e-5, so captured contrasts read 1.3–3.0 nat too high.

**3. The second pass's draw manufactures the false gDNA witnesses. Measured.**
- 97 % of g00's false gDNA is drawn. Calibrating on pass one alone gives c1: 202,229 / 26,423 against 201,046 / 27,710.
- The draw's gDNA term ignores source and strand: it sums a copy's traffic over both strands. As a result:
  - RNA multimappers split evenly among copies that have no unique traffic, and become intergenic gDNA "seeds";
  - antisense multimappers at inverted paralogs become gDNA too.
- The seeds are 92 of the 107 located training regions.

**4. Counting basins acts as a library-level switch. Measured.**
- About 100 false fragments create a second basin, and then 74 objects whose bookings never changed read about 5.5×10⁴.
- The refit loop's E-step grows that basin even when the training data are held fixed (it appears at refit 3).

**5. Where an object's own evidence is flat, only the population can read capture. Measured.**
- In the deferred stratum, a captured single-strand object's own strand evidence for its true level is −0.2 nat
  (regions) and −1.1 nat (boundaries). The population-shaped booking carries +167 and +1,314 nat for it.
- Every reader built on an object's own evidence lost there, by ×7.5 to ×9.7.
- So the population is both the cause of the g00 tail and the source of capture's benefit.

**6. The global landscape spreads a local change. Hypothesis.**
- The exact-copy control changed 980 sites. At stranded capture ON, gDNA called RNA rose ×2.7–2.8, mostly more than
  500 bp from any changed site.
- 0.7.1 was flat on the simulator's alignments and rose on STAR (1.19 M → 1.61 M).
- The landscape is the leading suspect, but sampling, exposures and possibly the length law changed too.
- It conflicts with the owner's ruling that capture inference is local.

**7. At duplicated families, capture's location is mostly the annotation's partition. Measured on one family.**
- Strand agrees across the copies in 99.9 % of shared records.
- The simulator captures by coordinates, so how this transfers to real hybridisation is a hypothesis.

**8. The reader amplifies small placement changes. Measured.**
- At LINC00319, about 25 changed placements tripled the booked gDNA (34 → 101), and the EM moved about 5,500 RNA
  fragments onto the gene (4,976 → 10,500; truth 91).

**9. Unannotated transcription is booked as gDNA. Measured.**
- This is Axiom 0's lock: where the annotation admits no RNA, every fragment is gDNA.
- With lncRNA withheld, 56 % of its RNA is called gDNA.
- On VCaP, 222 regions and 2,960 boundaries read above e×.

## 3. Facts any new reader must respect

**Capture is worth a lot where there is gDNA:**
- c1 costs stranded capture ON 4.3–5.0 M transcript errors, against the shipped reader's 0.53–0.75 M.
- In the deferred stratum, shipped reads 0.15 / 0.09 of c1.
- The benefit is −0.41 M to −4.09 M transcripts per capture-ON condition at g05 and above.

**It is harmful where there is none:**
- At g00 κ0.99 OFF, shipped makes 1.90 M transcript errors against c1's 0.20 M and 0.7.1's 1.38 M.
- At g00 κ0.99 ON: 1.19 M against 0.65 M.

**The g00 damage is a small tail:**
- 229 regions and 581 boundaries read above 10× carry 99 % / 100 % of the gap; the 57 k small lifts carry nothing.
- On the simulator's alignments HEAD equals c1 on every g00 row: the tail needs real multimappers.

**On a correct tally** (truth-placed), the reader sits near its own c1 where there is no captured gDNA (×1.02–1.08 /
×1.08–1.14), and beats c1 by ×2.3–31 where there is.

**Sparse capture:**
- c1 beats shipped on transcripts at 4 of 5 sparse rows, and loses genes ×3–5.
- On the single-probe screen, every reader built reads the target gene within 2.2 % of truth; c1 reads +5.5 % to +115 %.

**Strand:**
- Spliced strand counts are binomial (z-score sd 1.00).
- At a busy single-strand object, strand cannot see a gDNA share below κ's flip rate (about 1 %).
- At κ = ½, RNA and gDNA have identical emissions at an object.
- Pure gDNA at identical copies is not identifiable.

**Noise:**
- One draw seed moves g05 κ0.99 ON genes between 88,898 and 107,022 (sd 8.6 %).
- c1 moves with the tally, so compare every arm with its own c1.

**Multimapped gDNA under capture is real:** pass one alone costs g05 κ0.99 ON 380,977 / 101,845 → 445,971 / 142,391.

**Plentiful gDNA:** where gDNA is plentiful the bulk is right (1.90e-4 against 1.87e-4 per bp at g05 κ0.99 ON).

**Panel limits:** the simulator captures by coordinate overlap, and the ladder has equal gDNA and RNA lengths and no
unannotated transcription.

## 4. The ideas with the most merit

Ranked by potential to succeed. Each would replace an operation, not add a layer.

**1. One population, fitted once from observation likelihoods** (the owner's "fewer statistical handoffs").
- **The idea.** Fit the gDNA population once, from each object's own count likelihood, and let it enter once. This
  deletes the three re-reads of item 1, the refit loop, the located-variance admission and the reliability weights.
- **Evidence.** The first candidate cut the switch row ×0.44 / ×0.36, and RNA called gDNA fell 13,140 → 1,331.
- **Why it failed.** It lost the deferred stratum ×7–9.5, for two named reasons:
  - its unsmoothed atoms became 16–21 reader basins;
  - its likelihood was strand-only, and at κ ≈ ½ that cannot hold an enriched mode.

  The owner ruled that this refutes that one fit, not the direction.
- **Next.** A continuous population fitted on each object's full evidence, including the neighbour rows that carry
  capture in the deferred stratum, with a reader that does not count basins. It must first survive the cheap contrast
  between g00 and unstranded capture.

**2. Treat the bulk as an estimated nuisance, and target the yield contrast** (an outside reviewer's idea).
- **The idea.** The EM uses only ratios of `Σ_o A_to c_o`. Conditioned on the total gDNA count, the bulk cancels and an
  object's share is `S_o c_o / Σ S_j c_j`. Zero counts then give an upper limit on the bulk, not a floor at the wall.
  One bulk estimate would serve the draw, the reader and the intergenic rate.
- **Evidence.** The wall is the amplifier at g00 (item 2).
- **Why it has not landed.** It was never tested.
- **Next.** The cheapest refutation: change only the grid's extent, with the origin fixed, and see whether the g00 tail
  moves. The risk is that the sparse over-read is part of shipped's sparse gene advantage.

**3. Shrink each object toward c = 1 by its own uncertainty** (the owner's sparse-capture ruling).
- **The idea.** Per object and continuous: `E[z | ẑ] = τ²/(τ² + v)·ẑ`.
- **Evidence.** Using 0.7.1's moderate factors, HEAD lands on 0.7.1's switch-row error. The g00 tail's bookings sit at
  14–17 nat² against the deferred captured exons' 4.5–5.2 nat², which is different information; nothing has shown it
  separates the two.
- **The constraint.** It must shrink an honest centre on the contrasts the EM uses. Shrinking the booking itself is the
  refused widened reader.

**4. Remove the false witnesses at their source: the draw's gDNA term.**
- **The idea.** The likelihood draw, `ld_both`, is one stated likelihood over placement, reading and source.
- **Evidence.** Off capture it lands the switch row on its own c1 (×0.106 / ×0.064). Unstranded OFF goes ×0.823 / ×0.624
  and stranded OFF ×0.977 / ×0.936, nearly all of the truth-placed gain.
- **Why it has not landed.** It loses g05 κ0.99 ON ×1.041 / ×1.194, entirely through the reader at duplicated families.
  Fixing those families' gDNA density at the truth recovers 93 % of that loss.
- **Next.** It waits for a reader that does not amplify placement changes (idea 1, 2 or 3).

**5. Keep the draw's deposits out of the population's training.**
- **The idea.** Train the landscape on pass one's unique mass. It is one mechanism with no constants, and an exact
  no-op on the simulator's alignments.
- **Evidence.**
  - g00 against own c1: κ0.99 OFF 9.44 → 4.63, κ0.99 ON 1.82 → 1.16, and κ0.50 → 1.00.
  - The contaminated strata stay at about ×1.00 / ×0.94–1.00.
- **Costs.** g05 κ0.50 OFF 1.155 → 1.190 and g50 κ0.99 OFF 0.999 → 1.057; −51 k gDNA on the unstranded rows; VCaP
  unchanged.
- **Next.** Re-run it on HEAD. It is a correction layer: delete it if idea 1 or 4 makes it redundant.

**6. Local capture: a footprint within an exonic unit.**
- **Evidence.**
  - 97–99 % of truly captured exons have both boundaries at a tenth of their level or more, against 15–31 % at g00.
  - As a prior, it rescues ambiguous captured objects from −2.23 to −0.61 nat.
- **Why it has not landed.** Every footprint reader built so far amplified the false witnesses.
- **Next.** Try it only after ideas 4 or 5 remove those witnesses. It must account for one molecule incrementing several
  boundaries.

**7. An explicit explanation for unannotated transcription,** in place of "anything unexplained is gDNA".
- **Evidence.**
  - Withheld lncRNA is called 56 % gDNA, and retained genes are damaged ×1.4–3.0.
  - VCaP's false population is unannotated.
- **Next.** Score it on the annotation-omission panel, by where the hidden RNA goes.

**8. Stable aggregate outputs with explicit uncertainty** (the owner's direction).
- **Evidence.** Within-gene reallocation is 49–92 % of HEAD's transcript error.
- **The caveat.** Genes fail too: the switch row has 414 k gene errors against c1's 27.7 k.
- **Next.** Define which aggregate is identified, then test it across draw seeds and omission.

**Parked**
- **The share frame with a count-noise term in ψ's arm.** It took g98 κ0.99 OFF genes 37,304 → 7,595, but raised gDNA
  called RNA in all 12 contaminated STAR rows. Deferred past 0.8.0.
- **The posterior-median readout.** It removes prior-only lifts (8,014 → 1,147 empty regions above 2×), but it misses
  the switch and costs genes ×1.17 on a correct tally. Owner: never bundle it.
- **The strand frame by an unbiased moment estimate.** Paused.

**Suggested first questions**, each a day or less with a stop rule:
1. On the exact-copy pair (g50 κ0.99 ON on the simulator's alignments), find the first quantity that goes wrong:
   landscape basins, bulk, saturation, per-object factors, ψ shares, or the length law.
2. Run idea 2's grid-extent test at g00.
3. Re-run idea 5 on HEAD.

## 5. The owner's rulings that bind this work

- **Capture is a spectrum.** It is per object and continuous, and local: a single-probe panel must not move a global
  landscape. No detector, status, library-level switch or gate.
- **Sparse capture: robustness before accuracy.** Shrink toward c = 1 per object, from its own evidence.
- **A repair replaces an existing operation and names what it deletes.** No correction layers, no filters for
  troublesome objects, no count threshold, no maximum-factor clamp, no invented acceptance margin.
- **Proven exact ties** get canonical arithmetic and the exposure-weighted equal-density convention. Near-ties stay in
  the likelihood. Identical emissions get an explicit convention.
- **Judging a change.**
  - gDNA called RNA ranks first, but a rise in it is not by itself a veto.
  - Read Brier and a calibration check beside the linear score, fractionally by read name, on identical read IDs.
  - A truth-fed arm is an oracle intervention, never a ceiling.
- **An exclusion argument with unproved steps** leaves an open question, not a refusal.
- **Still refused:**
  - fewer refits;
  - training frozen on the initial solve;
  - the split reader;
  - direction 1 (each object read once from its own evidence);
  - the widened reader with a new weight.

## 6. Dead ends

| tried | what killed it |
|---|---|
| Widened booking | closed only 24 % / 26 % of the g00 tail |
| Each object read once from its own evidence | stranded ON ×2.74 / ×2.00; deferred ×8.06 / ×7.23 |
| Split reader (own evidence at single-strand objects) | deferred ×7.5 / ×5.7 |
| Footprint readers (pooled, robust consensus) | 1,205–1,572 objects above 10× at g00 |
| Markov-coupled co-enrichment | 224 / 305 objects above 10×; invents capture on unstranded OFF (×2.03 / ×2.21) |
| Population NPMLE on strand-only likelihoods | deferred ×7–9.5; 15,677 regions read as "depleted" |
| One refit instead of three | g05 κ0.50 OFF +11 % / +16 % |
| Training frozen on the initial solve | the basin stays |
| Yield read from the unique count | capture-ON genes ×2 |
| Source-aware draw prototypes | g05 κ0.99 ON +9.6 % / +18.7 % |
| gDNA free per object in the draw (`ld_free`) | g05 κ0.99 ON ×1.054 / ×1.275; uncertified |
| A fixed landscape prior inside the draw | multimodal; stopped 1,086 nats below the best basin |
| Changing only the draw's floor value | 1,897,003 → 1,901,007 |
| Strand overdispersion fitted on spliced RNA | z-score sd 1.00: nothing to fit |
| Oracle ruler used as a ceiling | the gDNA pool collapses 500 k → 19 k at g05 κ0.99 ON |

## 7. Reusable tools

All are under `$RIGEL_RUNS/prototypes/`.

**References and dissection**
- `2026-10-30_plan_review/truth_draw.py CONDITION OUT.jsonl` is the truth-placed reference. Its arms are `truth_place`,
  `truth_full` and `truth_matched`. `--trace` prints each refit's basins, and `--dump` writes per-object capture and
  booking. `seed_spread.py` gives a draw A/B's noise floor.
- `2026-10-29_reader_census/`:
  - `reader_census.py`: shipped against c1;
  - `reader_dump.py`: per object, f_g, Var(log f_g), k, S and the efficiency;
  - `reader_dissect.py`: objects far from 1, and the training set by basin;
  - `reader_ablate.py`: one population set to 1 per arm.
- `2026-10-29_reader_once/`:
  - `reader_bundle.py` / `bundle_lib.py`: ψ's evidence with the landscape divided out, refused unless it reproduces ψ;
  - `trace_false_gdna.py`: names every fragment in a flagged object by read name;
  - `ab_once.py --train-unique`: idea 5.
- `2026-10-30_switch_reader/deferred_check.py`: per object, the posterior mass on the bulk and the evidence at truly
  captured and truly empty objects.

**Families and the draw**
- `2026-10-31_likelihood_draw/`:
  - `families.py`, `capture_family.py`, `dissect_family.py`: one duplicated family as an explicit problem;
  - `ldraw.py`, `vsolve.py`, `vcert.py`: the likelihood draw and its certified solver, with gates 1–15;
  - `cache/`: STAR pass-one caches for four conditions.

**Scoring and decomposition**
- `2026-11-02_count_deficit/decomp.py`: the per-locus split of the EM prior share's error by class, refit and object
  size, about 25 s per condition.
- `2026-11-02_release_comparison/`: one common scorer for HEAD and 0.7.1, including 0.7.1's capture factors on HEAD.
- `2026-11-02_preview_validation/`: the paired scorer (`paired.py`), the omission and copy-control panels (`omission.py`,
  `copies.py`) and their tables.
