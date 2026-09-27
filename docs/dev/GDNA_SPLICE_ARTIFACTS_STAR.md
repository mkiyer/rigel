# STAR Parameters and Spurious Spliced Alignment of gDNA

**Owner:** `hulkrna` pipeline (`resources/star_arriba_parameters.txt`)
**Status:** analysis and recommendations; **no parameter changes made**
**Parent:** [GDNA_SPLICE_ARTIFACTS_OVERVIEW.md](GDNA_SPLICE_ARTIFACTS_OVERVIEW.md)
**STAR version:** 2.7.11b · index built with `--sjdbGTFfile genes_controls.gtf --sjdbOverhang 149`

---

## 0. Position statement

**STAR tuning is a mitigation, not the fix.** The design goal is that
`rigel` + `alignable` make the system robust to aligner settings; every
parameter we tighten here is a hole that reopens on the next STAR version, the
next chemistry change, or the next parameter edit. This BAM is also consumed by
arriba, which mandates several of the current settings.

Accordingly this document recommends a **small number of changes defensible on
their own merits**, and explicitly identifies changes that are *not* worth making
even though they would reduce false positives.

## 1. Effective parameter set

STAR's `Log.out` prints only *re-defined* parameters, so effective values were
obtained by combining `Log.out` with `STAR --help` defaults for 2.7.11b.

### Explicitly set (`resources/star_arriba_parameters.txt`)

```
genomeLoad NoSharedMemory
outSAMtype BAM Unsorted
outSAMstrandField intronMotif
outSAMattributes NH HI AS NM MD MC
outSAMunmapped Within KeepPairs
outFilterMultimapNmax 50                    # default 10
peOverlapNbasesMin 10                       # default 0
alignSplicedMateMapLminOverLmate 0.5
alignSJstitchMismatchNmax 5 -1 5 5          # default 0 -1 0 0   <-- LOOSENED
chimSegmentMin 10
chimOutType WithinBAM HardClip
chimJunctionOverhangMin 10
chimScoreDropMax 30
chimScoreJunctionNonGTAG 0
chimScoreSeparation 1
chimSegmentReadGapMax 3
chimMultimapNmax 50
quantMode GeneCounts
```

### Defaults in force (not overridden) — the ones that matter here

| parameter | default | relevance |
|---|---|---|
| **`alignSJDBoverhangMin`** | **3** | **annotated-junction anchor floor — primary culprit** |
| `alignSJoverhangMin` | 5 | novel-junction anchor floor |
| `outFilterType` | `Normal` | `BySJout` would re-filter reads by surviving junctions |
| `scoreGenomicLengthLog2scale` | -0.25 | bias toward short genomic span |
| `scoreGapNoncan` / `scoreGapGCAG` / `scoreGapATAC` | -8 / -4 / -8 | motif penalties |
| `scoreGap` | 0 | no flat gap penalty |
| `outFilterMismatchNmax` | 10 | absolute mismatch cap |
| `outFilterMismatchNoverLmax` | 0.3 | mismatch rate cap |
| `outFilterScoreMinOverLread` | 0.66 | |
| `outFilterMatchNminOverLread` | 0.66 | |
| `alignIntronMin` | 21 | |
| `alignIntronMax` | 0 → ~589,824 | `(2^winBinNbits)*winAnchorDistNbins` |
| `winAnchorMultimapNmax` | 50 | |
| `outSJfilterOverhangMin` | 30 12 12 12 | **`SJ.out.tab` only — does NOT filter the BAM** |
| `outSJfilterCountUniqueMin` | 3 1 1 1 | `SJ.out.tab` only |
| `outSJfilterDistToOtherSJmin` | 10 0 5 10 | `SJ.out.tab` only |

> **Important trap:** the `outSJfilter*` family filters the *`SJ.out.tab` report*,
> not the alignments in the BAM — **unless** `outFilterType BySJout` is set, which
> makes STAR re-filter reads to those whose junctions survived. Tightening
> `outSJfilter*` alone changes nothing that `rigel` sees.

## 2. The primary mechanism: `alignSJDBoverhangMin = 3`

An unspliced gDNA fragment that happens to match **3 bases** on the far side of an
annotated junction can be stitched across it. Three bases match by chance
frequently, and the resulting alignment scores better than the unspliced
alternative because `scoreGap = 0` and `scoreGenomicLengthLog2scale = -0.25`
rewards the shorter genomic span implied by removing an intron.

### Evidence that this is happening, not theorised

Read-level measurement on the index case (see the parent document for the full
table). Surviving junction rows at **annotated** junctions with no blacklist entry:

```
anchor_min bin      n      mean NM    % mismatch within 5bp    median intron
1-3              5,472       0.14              0.24%               2,655
4-5                847       0.39             14.05%                 848
6-10             1,387       0.29              6.13%                 759
11-20            2,570       0.20              0.86%                 759
21-40            4,471       0.23              1.07%                 823
41+              6,861       0.22              0.93%                 811
```

The `1-3` bin is a large, sharply-delimited population with **pristine sequence**
(mean NM 0.14, essentially no mismatches near the junction). That is the
signature of a 3-base *perfect chance match*, not of a real junction and not of a
mis-aligned read.

**Consequence for detection strategy:** mismatch-based heuristics cannot
discriminate this class. Anchor length is the only signal. This is why the fix
belongs in `alignable`/`rigel` (catalogue coverage and anchor thresholds) or in
this parameter — not in mismatch filtering.

The `4-5` bin's jump to **14.05%** near-junction mismatches is the second
mechanism at work: `alignSJstitchMismatchNmax` permitting mismatched stitches once
the anchor is long enough to need them.

### Cohort-scale confirmation

35.7% of surviving annotated-junction rows have `anchor_min ≤ 10`. In the
`SJ.out.tab` view, 920 of 2,661 non-blacklisted annotated junctions (34.6%) have
a maximum overhang ≤ 10 across *all* crossing reads, versus **0.52–1.53%** in
clean controls — a 20–60× enrichment.

## 3. Recommended changes

### R1 — `alignSJDBoverhangMin`: 3 → 8 (easy win)

```
alignSJDBoverhangMin 8
```

**Mechanism:** requires 8 rather than 3 bases past an annotated junction.
Chance-match probability falls by roughly `4^-5` ≈ 1/1000 per candidate site.

**Expected benefit:** directly targets the 5,472 pristine `anchor_min ≤ 3` rows
and much of the 847 + 1,387 in the 4–10 bins. Together these are 7,706 junction
rows / 1,398 junctions / **5,163 fragments** in the index case
(`evidence/bug4_sjdb_stitch_no_entry.tsv`), plus a share of the blacklisted-but-
escaping population.

**Expected cost:** genuine junctions supported *only* by reads with < 8 bp
overhang are lost. With 138 bp mates and 278 bp merged reads this is a small
population, and such junctions are poorly-supported evidence in any case. The
`11-20` and `41+` bins show real junctions are not concentrated at tiny anchors.
**Quantify before adopting** — this is the one number this investigation did not
measure directly.

**Value choice:** 8 rather than 10 or 12 keeps the change conservative. STAR's own
`--twopassMode`-oriented recommendations and ENCODE long-RNA settings use values
in the 1–8 range; 8 is at the strict end of conventional practice rather than
beyond it.

**Consistency requirement:** if adopted, the `alignable` simulation must be
regenerated with the same value, or the blacklist will contain artifact classes
that production can no longer produce (harmless) while its anchor statistics
become miscalibrated relative to production (not harmless). See
[GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md](GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md).

### R2 — `alignSJstitchMismatchNmax`: `5 -1 5 5` → `3 -1 3 3` (easy win)

**Current setting is looser than STAR's default** (`0 -1 0 0`) on three of four
motif classes. The four positions are
`(non-canonical, GT/AG, GC/AG, AT/AC)`; `-1` means unlimited.

**Mechanism:** we currently permit up to 5 mismatches when stitching across
non-canonical, GC/AG and AT/AC junctions. This directly enables the mismatched
stitches visible as the 14.05% spike in the `4-5` anchor bin.

**Supporting evidence — the GC/AG blind spot.** Among the index case's novel
junctions, GC/AG-family motifs (`intron_motif ∈ {3,4,5,6}`) are **0.03%** of
blacklisted (1 of 3,000) but **8.62%** of non-blacklisted (452 of 5,247) — a
**287×** asymmetry. Non-canonical junctions are near-absent everywhere because
`outSAMstrandField intronMotif` suppresses alignments with undefined strand.

**Caveat — arriba dependency.** This value came from arriba's recommended STAR
invocation. Arriba's documentation specifies loose stitch settings to retain
fusion-spanning reads. **Do not change this without checking arriba's current
recommendation and re-running fusion calling on a known-positive library.** Going
to `0 -1 0 0` (STAR default) is likely too aggressive for arriba; `3 -1 3 3` is
the compromise proposed here. Treat as medium-confidence pending that check.

### R3 — `outFilterType BySJout` (evaluate; do not adopt blind)

**Mechanism:** STAR performs a second pass keeping only reads whose junctions
survived the `outSJfilter*` family. This is the *only* way to make
`outSJfilterOverhangMin 30 12 12 12` and `outSJfilterCountUniqueMin 3 1 1 1`
affect the BAM.

**Why it is attractive here:** `outSJfilterCountUniqueMin` implements
"a junction needs ≥N independently-supporting unique reads." Artifact junctions in
this library are overwhelmingly thin: 78.7% of surviving novel junctions have
**zero** uniquely-mapped crossing reads, and 45.1% of non-blacklisted annotated
junctions have exactly one.

**Why to be cautious:** it changes the read set globally, interacts with the
chimeric detection arriba depends on, and makes the BAM's junction content
dependent on library depth — a property that undermines cross-library
comparability. It also violates the robustness principle in §0 more than R1/R2 do.

**Recommendation:** measure it on the retained FASTQs, but prefer the equivalent
logic implemented in `rigel` where it can be applied with full knowledge of the
gDNA model. See [GDNA_SPLICE_ARTIFACTS_RIGEL.md](GDNA_SPLICE_ARTIFACTS_RIGEL.md) §5.

### R4 — `scoreGenomicLengthLog2scale`: leave alone

Default `-0.25` rewards shorter genomic span, which mildly favours introducing an
intron. Making it less negative (e.g. `-0.5` is *more* negative; `0` removes the
term) would reduce that bias but affects **all** alignment scoring including
mate-gap handling and chimeric detection. The blast radius is large and the
expected specific benefit is small. **Not recommended.**

### R5 — `alignSJoverhangMin` (novel junctions): 5 → 12, low priority

The novel-junction artifact class does **not** have a short-anchor signature —
only 5.8% have `anchor_min ≤ 10`, and its median anchor is 40. Raising this floor
would cost true novel junctions while barely touching the artifacts. **Not
recommended** as an artifact control, though 8–12 is defensible on general
grounds if novel-junction precision is independently desired.

## 4. What will NOT help

| candidate | why not |
|---|---|
| tightening `outSJfilter*` alone | filters `SJ.out.tab`, not the BAM — `rigel` never sees it |
| mismatch filters for the class-A artifacts | the `anchor_min ≤ 3` population has mean NM **0.14** — there is nothing to filter on |
| `alignIntronMax` reduction | escaping novel junctions have *longer* introns (median 5,433) than blacklisted ones (2,037), but real junctions span the same range; no separation |
| `outFilterMismatchNoverLmax` reduction | would hit class B but also genuine variant-rich reads; `rigel` can weight this more precisely |
| lowering `outFilterMultimapNmax` from 50 | `rigel` runs with `include_multimap = true` and *uses* multi-placements; suppressing them at the aligner loses information rigel wants. Handle in rigel instead |

## 5. Testing protocol

The exact STAR input FASTQs are retained, so parameter variants can be tested
without re-fetching or re-trimming:

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/runs/human/\
mctp_LBX0077_SI_44263_HKWNYDRX7/fastq/k2_prefilter_{1,2}.dedup.fastq.gz
```

### Sweep design

Run STAR with the production parameter file plus one override per arm. Use the
production index (`refs/human/star_index`, `sjdbOverhang 149`) — no rebuild needed
for any of R1–R3, R5.

```bash
export PATH=/home/mkiyer/sw/miniforge3/envs/star/bin:$PATH
IDX=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/star_index
PF=/home/mkiyer/proj/hulkrna/resources/star_arriba_parameters.txt
FQ=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/runs/human/mctp_LBX0077_SI_44263_HKWNYDRX7/fastq

STAR --runThreadN 8 --genomeDir $IDX --parametersFiles $PF \
     --readFilesIn $FQ/k2_prefilter_1.dedup.fastq.gz $FQ/k2_prefilter_2.dedup.fastq.gz \
     --readFilesCommand "gunzip -c" \
     --outSAMattributes NH HI AS NM MD MC jM jI \
     --outFileNamePrefix out_<ARM>/ --outTmpDir <TMP>/<ARM> \
     <ARM OVERRIDE>
```

Adding `jM jI` to `outSAMattributes` makes STAR emit per-junction motif and
annotation status directly on each record, which simplifies per-read scoring.
(Production does not set these; they are for analysis only.)

**Arms:** `prod` (baseline) · `sjdbOH8` · `sjdbOH10` · `stitch333`
(`3 -1 3 3`) · `stitch000` (`0 -1 0 0`) · `bySJout` · `combo`
(`sjdbOH8` + `stitch333`).

Use SLURM account `mkiyer99` (the pipeline's account) for consistency.

### Metrics per arm

Compute both a sensitivity and a specificity axis. **A single-library sweep can
only measure the specificity axis** — the index case has almost no genuine
splicing, so true-positive loss must be measured on a clean control. Run the same
arms on `mctp_MI_1089_SI_40085_HKF7JDRX5` (gDNA 0.001) — its FASTQs are not
retained, so this requires a second targeted rerun (see the parent document §3.1
for the command pattern).

| axis | library | metric |
|---|---|---|
| false positives removed | index case | `Log.final.out` "Number of splices: Total"; count of spliced records with `anchor_min ≤ 10` at annotated junctions; `rigel` `strand_specificity` and `splice_artifact/all_spliced` after re-quant |
| true positives retained | clean control | "Number of splices: Annotated (sjdb)"; junctions with `maxoh > 20 & uniq ≥ 5`; `mrna_fraction` stability |
| fusion calling intact | a known-positive fusion library | arriba call concordance — required for R2 |

**Primary decision metric:** index-case `strand_model.strand_specificity` moving
from 0.7856 toward the 0.991–0.999 range seen in clean siblings, **while** the
clean control's annotated-splice count and `mrna_fraction` stay within ~1%.

## 6. Recommendation summary

| ID | change | confidence | cost | adopt? |
|---|---|---|---|---|
| R1 | `alignSJDBoverhangMin 3 → 8` | high mechanism, unmeasured TP cost | low | **yes, after measuring TP cost** |
| R2 | `alignSJstitchMismatchNmax 5 -1 5 5 → 3 -1 3 3` | high mechanism, arriba risk | low | **yes, after arriba check** |
| R3 | `outFilterType BySJout` | plausible, broad side effects | medium | measure; prefer rigel-side equivalent |
| R4 | `scoreGenomicLengthLog2scale` | low benefit, wide blast radius | — | no |
| R5 | `alignSJoverhangMin 5 → 12` | targets wrong population | low | no |

**Both R1 and R2 should be treated as consistency-coupled to `alignable`:** if
production stitching rules change, the artifact catalogue must be regenerated
under the same rules or its anchor statistics no longer describe production.
