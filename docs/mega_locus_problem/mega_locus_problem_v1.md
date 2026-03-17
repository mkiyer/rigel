Created memory file gdna-siphon-root-cause-2026-03.md

Starting: *Propose fix approaches* (5/5)

---

## Complete Root Cause Analysis: gDNA Siphon with minimap2

### The Headline

**ONE mega-locus** containing 217,246 transcripts from 43,049 genes absorbs 10.6M EM fragments. The single gDNA component in this locus receives **4,719,192 fragments (44.5%)** — accounting for **99.2% of ALL minimap2 EM gDNA**. Without this mega-locus, minimap2's small loci produce only 38K gDNA (comparable to oracle).

### The Mechanism: Universal Candidate θ Accumulation

The EM has `2×n_t + 1` components per locus (mRNA + nRNA per transcript, plus ONE gDNA). In the mega-locus this is **434,493 components**. The structural asymmetry:

| Component | Candidate for | θ accumulation |
|-----------|---------------|----------------|
| nRNA_j | ~100-1000 fragments overlapping transcript j | small |
| mRNA_j | ~100-1000 fragments overlapping transcript j's exons | small |
| **gDNA** | **ALL ~8.4M unspliced fragments in the entire locus** | **massive** |

The gDNA component is the **only universal candidate** — it competes for **every** unspliced fragment. This creates a positive feedback loop:

1. gDNA absorbs antisense reads easily (gDNA: `log(0.5)` = −0.693 vs nRNA: `log(0.1)` = −2.303 → gDNA wins by 1.6 nats)
2. θ_gDNA grows from these accumulations
3. For sense reads, nRNA only has a 1.8× likelihood advantage (`log(0.9)` vs `log(0.5)`)
4. But θ_gDNA >> θ_nRNA_j, so the product θ×L favors gDNA even for sense reads
5. More sense reads flow to gDNA → θ_gDNA grows further → cycle amplifies

### What's NOT the Cause

| Hypothesis | Test | Result |
|-----------|------|--------|
| nRNA sparsity prior (α=0.9) | Changed α to 1.0 | −12K shift (negligible) |
| VBEM vs MAP mode | Tested both | −35K shift (negligible) |
| **gDNA EB prior (gdna_init=868K)** | **Reduced scale 100×** | **−388K shift (only 8%)** |
| Alignment-level scoring | Checked non-mega loci | **Normal gDNA rates (7-17%)** |

The siphon is **LIKELIHOOD-driven, not prior-driven**. Even with near-zero gDNA prior, the EM converges to 4.33M gDNA in the mega-locus.

### Key Numbers

```
TRUTH:   mRNA=4,280,764   nRNA=5,719,236   gDNA=2,000,000

MEGA-LOCUS (minimap2):
  Transcripts: 217,246 (of which 180,938 are UNEXPRESSED)
  Genes: 43,049
  EM fragments: 10,596,856
  EM components: 434,493
  mRNA predicted: 3,539,404 (truth for these tx: 4,019,742)
  nRNA predicted: 2,338,260 (truth for these tx: 5,677,485) ← 3.3M DEFICIT
  gDNA predicted: 4,719,192 ← 3.85M EXCESS over gdna_init

  nRNA deficit (3.34M) ≈ gDNA excess (3.85M)

NON-MEGA LOCI (minimap2):
  8,755 loci, all < 100 tx
  gDNA total: 38,631 (normal)
  nRNA total: 44,445 (vs truth 41,751 — essentially perfect)
```

### Prior Mass Distribution (mega-locus)

```
gDNA:    868,328   (95.5% of total prior!)
nRNA:    ~38,700   (0.9 × ~43K nRNA components)
mRNA:     ~2,200   (0.01 × 217K + tiny OVR)
```

The EB formula: `gdna_init = shrunk_density × Σ(exonic_bp)`. For 217K transcripts, the sum of exonic lengths is enormous, giving gDNA a **95.5% share of the initial prior mass**.

### Proposed Fix Approaches

The fix must address the **θ imbalance**, not just the prior (which only accounts for 8% of the problem).

**Option A: Per-Gene gDNA Components**
Change the EM from `2×n_t + 1` to `n_t + n_nrna + n_genes` components. Each gene gets its own gDNA component that only serves fragments overlapping that gene's region. This eliminates the universal-candidate advantage and is the most principled fix.

**Option B: Mega-Locus EM Decomposition**
For loci exceeding a gene count threshold (e.g., >100 genes), decompose the EM into sub-problems (e.g., per-chromosome-arm or per-gene-cluster). Each sub-problem gets its own gDNA component. Shared multimapper fragments are allocated in a second pass.

**Option C: θ_gDNA Normalization**
In the E-step, replace θ_gDNA with `θ_gDNA / n_genes` when computing responsibilities. This is mathematically equivalent to per-gene gDNA components with equal θ, but requires no structural change to the CSR matrix.

**Option D: Per-Fragment gDNA Rate Constraint**
Instead of a single θ_gDNA in the EM, fix the per-fragment gDNA probability at the EB-estimated rate (e.g., ~8%) and run a modified EM where only mRNA/nRNA proportions are free. gDNA absorbs a fixed fraction.

My recommendation is **Option A** (per-gene gDNA) as the principled long-term fix, or **Option C** (θ normalization) as a quick interim fix to validate the hypothesis before committing to a larger refactor.