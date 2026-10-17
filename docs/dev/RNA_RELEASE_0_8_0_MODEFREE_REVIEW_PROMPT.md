# Reviewer prompt: the mode-free capture reader (copy everything below the line)

---

You are an independent technical reviewer. Your background is empirical-Bayes statistics, mixture models and RNA-seq
quantification. You are reviewing one component of **Rigel**, a Bayesian RNA-seq quantifier that separates RNA from
genomic-DNA (gDNA) contamination: its **capture reader**.

Hybrid capture enriches fragments that overlap probes by orders of magnitude. gDNA is a uniform template before
capture, so calibration reads each genomic object's **capture efficiency** from its deconvolved gDNA density, against
a population prior: the gDNA-density "landscape". The EM then contracts every component's effective length by those
efficiencies.

**The decision already taken.** The owner will ship Rigel 0.8.0 from the "count foundation" (commit `3e7e0fbc`), with
a new **mode-free reader** that replaces the old reader's dependence on a located enriched mode. Integration waits on
your review.

Your job:

1. Find the weaknesses and the improvement opportunities in the mode-free reader before it is integrated.
2. Assess a proposed follow-up, the **probe-footprint mechanism**.

## Materials

- **The review packet:** `docs/dev/RNA_RELEASE_0_8_0_MODEFREE_REVIEW.md`. It covers the theory, the exact
  implementation, the full measured matrix, the losses with their mechanisms, the refused variants, the footprint
  proposal, the integration plan and eight questions. **Read it first and completely.**
- **The evidence bundle:** `rigel_modefree_review.zip`. `HANDOFF.md` describes it and `HANDOFF_SHA256.json` verifies
  it. It contains:
  - the foundation's source;
  - the shipped arm's exact diff and a standalone listing of the reader, verified equal to what ran;
  - every arm's per-condition JSONL results and receipts;
  - the real-library results;
  - the scripts that generated every table;
  - the project's permanent docs and `CLAUDE.md`, which define the vocabulary and the rules.
- **If you have access to the machine the runs used**, the packet's §11 lists the live paths. The prototype worktree
  is `<worktree>` and the run directory is `$RIGEL_RUNS/prototypes/2026-10-13_modefree/`.

## Rules your proposals must respect

These are the project owner's standing rulings. You may argue that a ruling costs accuracy, but only with a measured
or derivable cost. Do not quietly design around one.

1. **Axiom 0.** Every genomic object holds at most three populations: gDNA, RNA on the + strand and RNA on the −
   strand. "Mature" and "nascent" are not populations.
2. **Capture is a spectrum.** No capture detector and no library-level on/off gate. Capture must be read per object,
   continuously and locally. Probes can be sparse: one probed exon in a ten-exon gene leaves the other nine uncaptured.
   Gene-level pooling was refused for exactly that reason.
3. **No new empirical constant.** Every number must be derived or already defined in the code. If a proposal needs a
   constant, name it explicitly as a discussion item; do not hide it in a functional form.
4. **No tuning to the benchmark panels.** Real data are test inputs, never design inputs.
5. **Rulings need benchmarks.** Every proposed change must state how it would be A/B'd:
   - pinned scan (`--set scan.total_threads=1`), fractional assignment, `OMP_NUM_THREADS=1`;
   - one mechanism at a time;
   - read per stratum, never pooled;
   - on the deciding panels (ladder, both fragment-length gap panels), the robustness panels and the real libraries.
6. **In-scope strata:** unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON. Unstranded ×
   capture-ON is deferred: reported, not a release requirement.

Do not recommend a variant the packet already refused (§7) unless you explain why its measured failure would not
recur. Do not propose redesigning calibration's solve for 0.8.0; the reader consumes the solve's outputs as they are.

## What to deliver

A markdown report of about 3–6 pages, with these sections.

1. **Verdict on the mode-free reader, option (b)**, in at most 10 lines. Is it sound to ship as described? List every
   **blocker**: a defect that must be fixed before 0.8.0 ships, with the evidence that it is one.
2. **Findings**, numbered. Each one has:
   - a severity: **blocker**, **fix before integration**, or **consider later**;
   - the specific packet section, code line or table row it concerns;
   - the reasoning, labelled as *derived*, *measured* (cite the number) or *conjecture*;
   - a concrete, falsifiable test: what to run, and what result would confirm or refute it.
3. **Improvements to option (b)**, ranked by expected benefit against risk and effort. Address at least:
   - the null that collapses the bulk basin, and the flank flattening it causes (packet §6.1);
   - the log readout against the EM's linear use of efficiencies (§3.2, §6.1);
   - per-level saturation's dependence on the census basins (§6.3);
   - whether the solve's uncertainty can be used without the failure in §3.3, and to remove §6.2's false gDNA.
4. **Option (c), the probe-footprint mechanism** (packet §8):
   - Is the hypothesis right, that the in-scope loss is flank flattening and the missing information is local?
   - If you would build it: how `q_o` and the shoulder mass should be **derived**, with no new constant, and how it
     must be gated and tested.
   - The expected gain against its development cost and risk, given that previous local-coupling readers were refused.
   - A go/no-go recommendation, and whether it belongs before or after 0.8.0.
5. **Integration risks** in the packet's §9: unit semantics, consumers that may assume efficiencies ≤ 1, edge cases
   (no landscape, one basin, empty or below-bulk basins, `S = 0`), and test strategy.
6. **Evidence check.** Do the tables support the claims in §1? Could anything in the protocol mislead: ratio-of-sums
   strata, realized-count truth in §6.1, the simulator's capture model, or comparison arms on different code paths?
7. **Answers to the packet's §10 questions**, one short paragraph each.

Be specific and terse. A finding without a test is a conjecture: say so. Where you are uncertain, give the experiment
that would settle it rather than a judgement.
