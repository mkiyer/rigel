# Plan: nRNA Architecture Redesign — Single-Exon Transcripts as Nascent RNA Equivalents

Redesign nRNA modeling so that annotated single-exon transcripts naturally serve as nascent RNA equivalents for overlapping multi-exon transcripts. Only create synthetic nRNA shadows when no annotated equivalent exists. Single-exon transcripts never emit separate nRNA candidates — they ARE both mature and nascent RNA. Recommended approach: diagnostic-first (Phase A) to quantify scope before implementation.

---

### Advantages

**1. Eliminates redundancy at annotated loci.** Currently, if single-exon S spans multi-exon T, the EM sees THREE competing entities: mRNA(T), mRNA(S), nRNA(shadow). Under the proposed model: only mRNA(T) and mRNA(S). The nRNA shadow was scoring identically to mRNA(S) for unspliced fragments, so the redundancy was a siphon with zero information gain.

**2. Scoring equivalence verified.** I traced the exact C++ code paths — for unspliced fragments, mRNA scoring against single-exon S produces the *same* log-likelihood as nRNA scoring against the shadow:
- **Fragment length**: `frag_len_log_lik(flen)` where `flen` for S = genomic footprint (exon = span) — identical to nRNA's `frag_len_log_lik(genomic_footprint)`
- **Overhang**: `oh = rl - ebp` for S vs `oh = rl - (ebp + ibp)` for nRNA. For fragments fully within S, `ebp(S) = ebp_T + ibp_T`, so overhangs match
- **Coverage weight**: S's exonic length = S's span = nRNA's span when they share coordinates
- **Strand and NM**: identical in both paths

**3. Intronic fragment routing is naturally correct.** Fragments in T's introns get `ebp=0` against T (rejected as mRNA for T) but `ebp>0` against S (accepted as mRNA for S). This works because `_gen_transcript_intervals()` in index.py generates both EXON and TRANSCRIPT intervals, and S's EXON interval covers the entire span including T's introns. Verified in resolve_context.h (UNION merge).

**4. Eliminates all zeroed-out nRNA entries.** Currently, single-exon-only genes create nRNA entries that are then zeroed in locus.py. That entire class of waste disappears.

**5. Reduces EM component count.** Each annotated-equivalent locus loses one nRNA component. Single-exon-only genes lose their dummy nRNA entirely.

**6. More principled locus connectivity.** Currently, nRNA sharing via union-find fanout ([em_solver.cpp lines 2160-2185](src/rigel/native/em_solver.cpp#L2160)) artificially links transcripts by nRNA index. Under the proposed model, connectivity comes from actual fragment co-occurrence in cgranges overlaps — biologically more meaningful, avoids mega-loci from coincidentally similar coordinates.

**7. Clean mental model.** "A single-exon transcript IS both mature and nascent RNA" is a simple, biologically correct statement. No more create-skip-zero dance.

---

### Disadvantages / Challenges

**1. `t_to_nrna_arr` dense-array invariant breaks.** Currently every transcript has a valid `nrna_idx` in `[0, num_nrna)`. Under this model, some transcripts have no nRNA entity. Best handled with a -1 sentinel + ~11 guard-clause insertions across Python and C++ (`if idx < 0: continue`). Feasible but not trivial.

**2. Index-time containment detection (~20-40 new lines).** A new step in `compute_nrna_table()` must check whether each multi-exon T is fully contained by a single-exon S on the same strand/ref. Coordinate comparison after sorting — straightforward but new logic.

**3. Output interpretation changes.** For annotated-equivalent loci, the "nascent RNA signal" is embedded in the single-exon transcript's mRNA abundance. Needs an `is_nascent_equivalent` flag in output and documentation, otherwise users may misinterpret single-exon abundance.

**4. Locus topology changes (regression risk).** Removing nRNA fan-out changes locus groupings. Changes should be *more* principled, but constitute a regression risk requiring benchmark validation.

**5. Partial coverage of the problem.** The MAJORITY of multi-exon genes likely do NOT have annotated single-exon equivalents. For those genes (the common case), synthetic shadows are still needed — behavior unchanged. The benefit is concentrated on: (a) genes with single-exon variants (MALAT1-like), and (b) pure single-exon genes (eliminate dummy nRNA).

**6. Two-type nRNA system.** The codebase must consistently handle both annotated equivalents (regular transcripts, not in nRNA table) and synthetic shadows (in nRNA table). While conceptually clean, both paths must be tested.

---

### Steps

**Phase A: Diagnostic** (no code changes)
1. Write `scripts/debug/nrna_diagnostic.py` — load GENCODE index, compute: count of single-exon transcripts fully containing multi-exon transcripts, count of nRNA entries eliminated vs retained, locus size distribution changes
2. Run on a real reference index
3. Report findings → go/no-go decision

**Phase B: Index changes** (*depends on Phase A go*)
1. Modify `compute_nrna_table()` in index.py — containment detection
2. Add `is_nascent_equivalent` field to `Transcript` in transcript.py
3. Use -1 sentinel in `t_to_nrna_arr` for transcripts without synthetic nRNA
4. Update nRNA feather serialization

**Phase C: Scoring and EM changes** (*depends on Phase B*)
1. Add -1 guards in scoring.cpp (~2 locations: singlemapper + multimapper paths)
2. Add -1 guard in em_solver.cpp (union-find fanout, 1 location)
3. Add -1 filters in locus.py (2 locations), estimator.py (2 locations)
4. Simplify prior zeroing in locus.py (remove single-exon check)

**Phase D: Validation** (*parallel with late Phase C*)
1. Recompile, `pytest tests/ -v`, `pytest tests/ --update-golden`
2. Run benchmark suite, compare nRNA siphon / mRNA / gDNA errors against baseline
3. Spot-check: gene with known single-exon variant → verify locus construction and EM output

**Relevant files**: index.py (`compute_nrna_table` lines 114-159), transcript.py (line 28), scoring.cpp (lines 1092-1212), em_solver.cpp (lines 2160-2185), locus.py (lines 71-145, 420-433), estimator.py (~line 386)

---

### Verdict

**Favorable.** The architecture is biologically sound, scoring-equivalent (verified at the C++ level), and produces a cleaner codebase. Implementation is tractable (~50-70 net modified lines + ~11 guard insertions). The primary uncertainty is **scope of impact** — how many loci are actually affected in a real annotation — which Phase A addresses before committing to implementation.

**Decisions captured**:
- Full containment only (S.start ≤ T.start AND T.end ≤ S.end), same strand/ref
- Multiple single-exon candidates compete naturally in EM
- -1 sentinel in `t_to_nrna_arr` (simplest semantics)
- EM solver itself unchanged — just fewer components