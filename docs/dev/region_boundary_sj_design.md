# The per-transcript RNA prior — design, findings, and what is left

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When a section settles, MOVE it — the vocabulary
    to `DESIGN.md` §0, a derivation to `EQUATIONS.md`, a lesson to `TRAPS.md`, a number to `ROADMAP.md`
    — and delete it here in the same edit.

    ⭐ Rewritten end to end 2026-08-13. Everything below is either MEASURED on the shipped tree or is an
    explicit open question. Nothing here is aspiration.

---

## 0. ⭐⭐⭐ THE ONE RESULT — a correct per-transcript prior WORKS

`g00 ss0.99 capture_off`, true relative abundances as the per-transcript allocation weights:

| | base | **truth, WARM START ZEROED** | truth, seed kept | FLIPPED, seed kept |
|---|---|---|---|---|
| transcript Σ\|err\| | 4,275,988 | **696,900 (0.16×)** | 858,395 (0.20×) | 5,398,909 (1.26×) |
| transcript fp_mass | 987,034 | **153 (0.00×)** | 36,775 (0.04×) | 226,666 |
| transcript fn_mass | 507,099 | **1,726 (0.00×)** | 2,109 (0.00×) | 101,596 |
| gene Σ\|err\| | 445,944 | **57,154 (0.13×)** | 260,603 (0.58×) | 888,783 (1.99×) |
| gene fp_mass | 23,529 | **0** | 2,209 | 1,951 |

**84 % of transcript misassignment removed, 87 % at gene level, with correct weights; 26 % worse with
deliberately wrong ones.** The machinery converts good weights into good answers, and a weighting
function is therefore worth designing. That was the whole question stage 5 existed to answer.

⭐ **The owner's original specification — zero the warm start and let the prior carry it — is the BEST
arm**, and by a clear margin at gene level (0.13× against 0.20×). The two arms differ in one thing: with
the seed kept, the coverage-weighted guess is added to the truth and dilutes it.

⚠ `warm_start="prior"` is only meaningful WITH a weight vector: under the shipped evidence-proportional
rule an all-zero seed leaves the RNA pool at zero, because `out[i]` is proportional to `raw[i]`.

⛔ **THREE LIMITS, AND THEY MUST TRAVEL WITH THE NUMBER.**

1. **It is a CAPABILITY proof, not a ceiling.** It does not price what a real weighting function could
   earn. Doing that needs controls this arm deliberately omits (§6).
2. **True weights also hand over the true SUPPORT.** A zero weight is EXACTLY absorbing (measured
   `0.0000` in all four mode × warm-start combinations), so every silent transcript is switched off for
   free — which is most of the fp_mass collapse to 0.04×. A real weighting function does not get that
   without solving a sparsity problem. 4,579 of 8,750 annotated transcripts are silent at `g00` and
   3,465 of them still compete for fragments.
3. **One condition, zero gDNA, stranded.** The blind stratum (unstranded × capture-ON) is untested, and
   it is the one carrying most of the tool's error.

---

## 1. ⛔⛔⛔ THE BUG THAT ALMOST INVERTED THE RESULT — read this before trusting any arm

The first two runs of stage 5 said the allocation was **inert**: truth weights matched `base` to within
one fragment. The conclusion drafted from that was *"the prior is too weak to matter."*

It was a plumbing bug. `pipeline._run_locus_em_partitioned` **accepted `rna_prior_weight` and never
passed it on** — the line was deleted while removing an unrelated caller, leaving a parameter that is
accepted and silently ignored. Every allocation, however extreme, produced byte-identical output.

⛔ **A FIRE COUNTER CANNOT SEE THIS.** The counter counted nonzero entries in the weight array — a
property of the array the harness was handed, true whatever the pipeline then did with it. What caught
it was **`--arm oracle_alloc_flip`**: inject a maximally WRONG allocation, which *must* move the answer,
and observe that it did not. ⭐ A correct-looking result would have been invisible; only the arm that had
to fail could reveal it.

⚠ Two defences now exist and both must stay: the flip arm, and a fire counter that watches the
**estimator** receive the array rather than the harness send it.

---

## 2. WHAT IS BUILT, AND WHAT EACH PIECE IS FOR

| stage | | state |
|---|---|---|
| **1** | transcript → ordered `(kind, id)` PATH over REGION / BOUNDARY / SPLICE JUNCTION, in TRANSCRIPTION order | ✅ `splice_graph.build_transcript_path`, 8 gates |
| **2** | per-object EVIDENCE: mass, count, `(f_g, f_pos, f_neg)` | ✅ published on `CalibrationResult` — ⛔ but see §3, it does not close |
| **3** | the per-transcript weight LANE into the EM | ✅ 24 gates, 12 perturbations, 0 holes |
| **4** | `EMConfig.warm_start` = `coverage` / `prior` / `uniform` | ✅ 6 gates |
| **5** | the PROOF — oracle weights | ✅ **§0** |
| **6** | the weighting function | ⛔ **not started — this is the next real work** |

### The transcript path (`build_transcript_path`)

A CSR of `offsets` / `kind` / `obj_id`, following the `(kind, obj_idx)` idiom `NodeChain` already uses,
with a third kind because a splice junction is neither a position nor an interval. Verified on the
shipped index: **0 of 15,669** transcripts differ from the independently-checked
`_transcript_node_incidence` on the region axis, **0** on the boundary axis, **45,609** splice steps —
exactly the annotated (transcript, intron) count — and transcription order holds on **7,312/7,312**
minus-strand and **8,357/8,357** plus-strand transcripts.

⭐ The splice-junction join goes through INTRON COORDINATES, never the flanking region pair. An intron
that resolves to no slot RAISES rather than silently shortening the path, because a shorter path still
reads as a well-formed walk.

### The lane

`float64[n_transcripts]` → pipeline (**flat**, never `[ids]`) → estimator → subproblem remap →
`out[i] = raw[i] + rna_prior·w_i/Σw`. Absent ⇒ the shipped path, unchanged.

⭐⭐ **The shipped rule is ALREADY this update at `w_i = raw[i]`** — `out[i] = raw[i]·(1 + P/A)` is
identically `raw[i] + P·raw[i]/A`. A prior that echoes the posterior carries no information, which is
exactly why it was neutral. `EQUATIONS.md` §9c has the derivation. **The two designs differ in the
weights alone.**

---

## 3. ⛔⛔ THE ψ DEFECT — a separate work path, deferred by owner ruling

The three-way composition does not close. Measured on `g00 ss0.99 capture_off`, every object addressed
by a chain slot:

| axis | sums to 1 | sums into (0,1) | sums > 1 |
|---|---|---|---|
| REGION | 74.72 % | **25.25 %** — median 0.978, p5 **0.869** | 12 |
| BOUNDARY | 77.24 % | **22.71 %** — median 0.979, p5 **0.850** | 16 |

`NodeDeconv` asserts `f_pos + f_neg + gdna_frac = 1`. It is false for a quarter of every axis. The
mechanism is visible at `sweep.py`'s write-back: the three posterior means are `np.clip(·, 0, 1)`-ed
**independently**, and an unsolvable slot keeps an init instead. ⚠ By linearity of expectation three
posterior means over ONE lattice should close, so the clip is a symptom rather than the whole cause.

⛔ It was invisible because nothing consumed the strand split: `chain_edge_deconv` published `0` for both
RNA strands, so the crossing axis's composition summed to `f_g` alone. That projection is fixed; the
non-closure is not.

⛔ **Do NOT repair it by renormalising at publication** — that makes a 15 %-short object
indistinguishable from a solved one. `ROADMAP.md` §2 item 0b is the tracked path.

---

## 4. IDENTIFIABILITY — the owner's objection, measured

*"A fragment sampled from TA may be exactly compatible with TA and TB. Your truth cannot be TA."*
Correct in kind. Measured in magnitude:

| | |
|---|---|
| fragments compatible with 2+ annotated transcripts | **70.43 %** |
| true mass in **likelihood-flat** directions (no estimator can see the split) | **0.29 %** |
| independent measure — null directions of the class-incidence matrix | **0.96 %** |
| exactly-exchangeable transcript groups, census over 952 loci / 132 M candidates | **zero** |
| true mass on transcripts with **no fragment of their own** | **23.99 %** (4,277 transcripts) |
| true mass with direct private evidence | **75.71 %** |

⭐ So per-fragment ambiguity is near-universal and per-transcript non-identifiability is ~0.3–1 %. An
oracle per-transcript prior is **not** meaningfully fake. The real difficulty is the middle tier: a
quarter of the mass rides on transcripts whose abundance is recoverable only through model shape.

⚠ **Gene collapse is the wrong identifiability axis in both directions** — too fine for 7.68 % of units
(compatibility classes spanning more than one gene), too coarse for 62.3 % of mass (multi-isoform genes
whose split IS recoverable). And a strict compatibility-class axis would silently BECOME the gene axis,
since the strict class is empty here.

---

## 5. THE GRAPH MODEL — vocabulary and structures (owner rulings, 2026-08-12)

**REGION** (was node) · **BOUNDARY** (was edge / line / cut / seam) · **SPLICE JUNCTION** (`sj`).
⭐ A convergence, not an invention: all three axes already carry two live names each (`n_regions` 9 files
vs `n_nodes` 20; `n_edges` 11 vs `n_cuts` 5; `n_junctions` 11 vs `n_sj` 5), and `RegionArrays`' own
docstring calls its rows "nodes". ⚠ `DESIGN.md` §0 is not updated until the bulk rename lands.

### A splice junction is anchored at TWO BOUNDARIES

    region r  ->  boundaries (r−1, r)          ARITHMETIC
    boundary b ->  regions   (b, b+1)          ARITHMETIC
    sj j      ->  boundaries (b_left, b_right) ARITHMETIC
    boundary  ->  splice junctions             CSR, BOTH WAYS

⭐⭐ The two CSRs (`sj_at_left`, `sj_at_right`) give BOTH things needed at once: the per-boundary total
is a CSR range sum, and the full `(src,dst)` pair table *is* the slot. ⚠ Only the left direction exists
today; the right-hand CSR is index-derived, no re-scan.

⛔⛔ **`donor` / `acceptor` are BANNED — they are 5'/3' and therefore strand-dependent.** Measured: the
code is genomically ordered and CORRECT (`src < dst` on all 13,482 junctions including all 6,527
minus-strand ones), so for **48.4 %** of junctions the identifier `donor_cut` holds the ACCEPTOR. The
names lie about correct code. Use `boundary_left` / `boundary_right`.

### Boundary EVENTS are a bitmask; a splice junction is not an event

`TSS | TES | STRAND_CHANGE | EXON_INTRON_CHANGE` are reasons a boundary exists and can co-occur ⇒ a
bitmask. A splice junction is a different object that REFERENCES two boundaries. Once said, the
`edges_df.kind == JUNCTION` discriminator disappears rather than being renamed.

### Five sj structures → two tables, one CSR, one flag

`junctions_df` (identity, with `ANNOTATED | BLACKLISTED` as flags — the blacklist stops being a separate
keyed structure) · `transcript_sj` CSR · the boundary event bitmask · the two CSRs above.

⛔⛔ **BLOCKING for the SJ OUTPUT (not for the prior):** `SJStrandTable` and the accumulator disagree
**1.576×** on the same condition — 8,582 junctions / 3,597,564 obs against 8,846 / 5,668,526 — because
`SJStrandTable` counts the leftmost junction only, unique mappers only. Publishing either first ships
two answers to one question.

### The mass ruling, reversed — with its premise

`accumulator.h` ruled a mass bank is ONE value because *"nothing reads a mass per strand"*. ⭐ That
premise is now false: artifactual splice junctions accumulate SYMMETRICALLY on both strands like gDNA
and are detectable by the existing strand model, and a per-strand mass is what makes artifact filtering
single-pass. **Record the reversal with the premise that changed**, or it gets re-litigated in both
directions. Price: the payload schema digest moves, every scan cache rebuilds — **8.3 s/condition, ~6
min for the ladder**. ⭐ Bundle it with the two dead banks `ROADMAP.md` §2.7 already wants removed.

### Two hazards any new design must not inherit

⛔ The `sj_id` is a **DENSE RANK** — dropping one annotated intron renumbered 13,476 of 13,481 surviving
slots. A published id is a within-run join key, never a cross-version identifier.
⛔ **Two reference orderings exist**: the load path groups `sj.feather` ALPHABETICALLY, the graph axis
uses FASTA order. Identical on this 2-chromosome index, divergent on any genome with chr1/chr2/chr10.
Every join must go through the coordinate key.

### The sharing rule — SETTLED, keep the shipped one

Every crossed boundary consumes mass. Each base votes exactly once, split between the at-most-two
objects bounding its slice, so mass is a **partition of the fragment's bases** (normaliser `L`). The
proposed alternative (weight an object by its adjacent bases) agrees at 1 and 2 objects and diverges
from 3 up — at `n=3` equal segments, `[0.375, 0.25, 0.375]` against `[⅓, ⅓, ⅓]`. Owner ruled for the
shipped rule. ⭐ Region mass stays implicit: containment is exclusive, so the contained count IS the mass.

---

## 6. WHAT IS LEFT, IN ORDER

1. ⭐⭐⭐ **STAGE 6 — the weighting function.** This is the work. The foundation is proven (§0), the
   evidence is published (§2), the path exists. ⛔ Before designing: re-read §0's three limits and
   `TRAPS: a-trap-names-the-defect-not-the-repair`.
2. **Re-run stage 5 on the blind stratum** (unstranded × capture-ON) and one mid-gDNA condition before
   generalising the 80 %.
3. **Price what a REAL weighting function could earn** — the ceiling, with the controls §0 omits: a
   support-only ablation (inject the oracle only where truth is nonzero) to separate deconvolution from
   sparsity, and a gene-aggregated oracle to price the within-gene split.
4. **ψ's non-closure** (§3) — owner ruled this a separate branch, after the prior lands.
5. **The graph refactor** (§5) — the rename, the two sj CSRs, the boundary-event bitmask. `sj_mass[2]`
   is `ROADMAP.md` §2 item 0c: it moves the payload schema digest, so bundle it with the two dead banks
   §2.7 wants removed and pay the ~6 min re-scan once. ⛔ One schema change at a time, or
   `rescan_panels.py`'s byte-identity gate cannot attribute the delta.
6. **The SJ output**, once `SJStrandTable` is reconciled.
7. **The index duplicate map** — `ROADMAP.md` §2 item 0d. ⚠ It needs an index rebuild but **NOT** a
   panel re-scan: the cache is keyed on `graph_hash` / `reach_digest` / `payload_schema_digest` /
   `scan_config_digest`, and an alias map on the same partition changes none of them. So it does NOT
   have to wait for the `sj_mass[2]` re-scan — the two were conflated once.
