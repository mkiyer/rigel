"""
Synthetic multimapping test framework for hulkrna.

Creates a genome with duplicated segments that produce secondary alignments
when reads are aligned with minimap2. Validates that hulkrna correctly
groups, pairs, and resolves multimapping fragments.

Genome layout
-------------
One "true" locus on seg_true with a 3-exon gene (t1).  Additional
duplicate segments (seg_dup1, seg_dup2, ...) that copy portions of
seg_true to create alignment ambiguity.  Only seg_true has GTF
annotation — the duplicates are "unannotated paralogs".

Exon coordinates on seg_true (0-based half-open):
  exon1: 500–700   (200 bp)
  exon2: 1000–1200 (200 bp)
  exon3: 1500–1700 (200 bp)
  intron1: 700–1000 (300 bp)
  intron2: 1200–1500 (300 bp)

Each seg_dup copies a contiguous region of seg_true that includes one
or more exons, causing reads from t1 to produce secondary alignments.
"""

import csv
import logging
import subprocess
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.bam import parse_bam_file, _group_records_by_hit
from hulkrna.fragment import Fragment
from hulkrna.index import HulkIndex, read_transcripts, write_bed12
from hulkrna.pipeline import run_pipeline
from hulkrna.resolution import resolve_fragment
from hulkrna.sim.genome import MutableGenome, reverse_complement
from hulkrna.sim.reads import SimConfig, ReadSimulator
from hulkrna.transcript import Transcript
from hulkrna.types import Interval, Strand

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TRUE_SEG = "seg_true"
TRUE_SEG_LEN = 3000
EXON1 = (500, 700)
EXON2 = (1000, 1200)
EXON3 = (1500, 1700)
INTRON1 = (700, 1000)
INTRON2 = (1200, 1500)
GENE_ID = "g1"
TRANSCRIPT_ID = "t1"
ABUNDANCE = 1000.0


def build_genome_and_gtf(work_dir: Path, n_duplicates: int = 1,
                         dup_region: tuple[int, int] = (400, 1800),
                         asymmetric: bool = False):
    """Build a FASTA with seg_true + duplicated segments, and a GTF.

    Parameters
    ----------
    work_dir : Path
        Output directory.
    n_duplicates : int
        Number of duplicate segments to add (1 = 2x multimapping, etc.)
    dup_region : tuple[int, int]
        Region of seg_true to copy into duplicates (default: covers
        all 3 exons with flanking sequence).
    asymmetric : bool
        If True, duplicate only the R2-side (exon1 region 400-800)
        so that R1 maps uniquely and R2 maps to multiple locations.

    Returns
    -------
    fasta_path, gtf_path : Path, Path
    """
    rng = np.random.default_rng(42)
    dna_bases = np.array(list("ACGT"), dtype="<U1")

    # Generate seg_true sequence
    true_seq_arr = dna_bases[rng.integers(0, 4, size=TRUE_SEG_LEN)]
    true_seq = list(true_seq_arr)

    # Inject splice motifs (GT-AG for + strand)
    # intron1: 700-1000 → GT at 700, AG at 998
    true_seq[INTRON1[0]] = "G"
    true_seq[INTRON1[0] + 1] = "T"
    true_seq[INTRON1[1] - 2] = "A"
    true_seq[INTRON1[1] - 1] = "G"
    # intron2: 1200-1500 → GT at 1200, AG at 1498
    true_seq[INTRON2[0]] = "G"
    true_seq[INTRON2[0] + 1] = "T"
    true_seq[INTRON2[1] - 2] = "A"
    true_seq[INTRON2[1] - 1] = "G"

    true_seq_str = "".join(true_seq)

    # Extract the region to duplicate
    dup_start, dup_end = dup_region
    if asymmetric:
        # Only duplicate the exon1 region (R2 side for FR library)
        # R2 maps to 5' of transcript → exon1. R1 maps to 3' → exon3.
        # Duplicating only exon1 region makes R2 multimap, R1 unique.
        dup_start, dup_end = 400, 800
    dup_fragment = true_seq_str[dup_start:dup_end]

    # Write multi-segment FASTA
    fasta_path = work_dir / "genome.fa"
    with open(fasta_path, "w") as f:
        # seg_true
        f.write(f">{TRUE_SEG}\n")
        for i in range(0, len(true_seq_str), 80):
            f.write(true_seq_str[i:i + 80] + "\n")
        # Duplicate segments
        for d in range(n_duplicates):
            seg_name = f"seg_dup{d + 1}"
            # Add small padding on each side to avoid edge effects
            pad_left = "".join(dna_bases[rng.integers(0, 4, size=200)])
            pad_right = "".join(dna_bases[rng.integers(0, 4, size=200)])
            dup_seq = pad_left + dup_fragment + pad_right
            f.write(f">{seg_name}\n")
            for i in range(0, len(dup_seq), 80):
                f.write(dup_seq[i:i + 80] + "\n")

    # Index the FASTA
    pysam.faidx(str(fasta_path))

    # Write GTF — only for seg_true
    gtf_path = work_dir / "genes.gtf"
    with open(gtf_path, "w") as f:
        attrs = (f'gene_id "{GENE_ID}"; transcript_id "{TRANSCRIPT_ID}"; '
                 f'gene_name "TestGene"; gene_type "protein_coding"; '
                 f'tag "basic";')
        for start, end in [EXON1, EXON2, EXON3]:
            # GTF is 1-based inclusive
            f.write(f"{TRUE_SEG}\ttest\texon\t{start + 1}\t{end}\t.\t+\t.\t{attrs}\n")

    return fasta_path, gtf_path


def simulate_reads(work_dir: Path, fasta_path: Path, gtf_path: Path,
                   n_fragments: int = 2000, seed: int = 42):
    """Simulate paired-end reads from the true transcript.

    Uses the hulkrna ReadSimulator to generate reads from t1.

    Returns
    -------
    r1_path, r2_path : Path, Path
    """
    # Build a MutableGenome-like object for the true segment
    genome = MutableGenome(TRUE_SEG_LEN, seed=seed, name=TRUE_SEG)

    # We need to overwrite with our actual sequence that has splice motifs
    fai = pysam.FastaFile(str(fasta_path))
    true_seq = fai.fetch(TRUE_SEG)
    fai.close()
    for i, base in enumerate(true_seq):
        genome._seq[i] = base

    # Build transcript object
    t = Transcript(
        t_id=TRANSCRIPT_ID,
        g_id=GENE_ID,
        ref=TRUE_SEG,
        strand=Strand.POS,
        exons=[Interval(*EXON1), Interval(*EXON2), Interval(*EXON3)],
        t_index=0,
        g_index=0,
        g_name="TestGene",
        g_type="protein_coding",
        is_basic=True,
        abundance=ABUNDANCE,
    )

    config = SimConfig(
        frag_mean=300.0,
        frag_std=60.0,
        frag_min=100,
        frag_max=600,
        read_length=150,
        error_rate=0.0,
        strand_specificity=1.0,
        seed=seed,
    )
    sim = ReadSimulator(genome, [t], config=config)
    r1_path, r2_path = sim.write_fastq(work_dir, n_fragments)
    return r1_path, r2_path


def align_reads(work_dir: Path, fasta_path: Path, r1_path: Path,
                r2_path: Path, gtf_path: Path) -> Path:
    """Align reads with minimap2 (--secondary=yes) → name-sorted BAM."""
    bam_path = work_dir / "reads.bam"

    # Build BED12 annotation for minimap2
    bed_path = work_dir / "annotation.bed"
    transcripts = read_transcripts(gtf_path)
    write_bed12(transcripts, bed_path)

    minimap2_cmd = [
        "minimap2", "-ax", "splice:sr",
        "--secondary=yes",
        "-j", str(bed_path),
        str(fasta_path),
        str(r1_path),
        str(r2_path),
    ]
    sort_cmd = ["samtools", "sort", "-n", "-o", str(bam_path), "-"]

    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    _, p2_err = p2.communicate()
    p1_err = p1.stderr.read()
    p1.wait()

    if p2.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {p2_err.decode()}")

    logger.info(f"Aligned → {bam_path}")
    return bam_path


def analyze_bam(bam_path: Path):
    """Analyze the BAM for multimap characteristics."""
    print("\n" + "=" * 70)
    print("BAM ALIGNMENT ANALYSIS")
    print("=" * 70)

    bamfh = pysam.AlignmentFile(str(bam_path), "rb")

    # Count reads by type
    total = 0
    primary = 0
    secondary = 0
    supplementary = 0
    unmapped = 0
    n_with_nh = 0
    nh_values = Counter()

    # Track per read-name statistics
    qname_hits = defaultdict(lambda: {"r1_primary": 0, "r2_primary": 0,
                                       "r1_secondary": 0, "r2_secondary": 0,
                                       "r1_supplementary": 0, "r2_supplementary": 0})

    for read in bamfh.fetch(until_eof=True):
        total += 1
        qname = read.query_name

        if read.is_unmapped:
            unmapped += 1
            continue

        if read.is_supplementary:
            supplementary += 1
            key = "r1_supplementary" if read.is_read1 else "r2_supplementary"
            qname_hits[qname][key] += 1
        elif read.is_secondary:
            secondary += 1
            key = "r1_secondary" if read.is_read1 else "r2_secondary"
            qname_hits[qname][key] += 1
        else:
            primary += 1
            key = "r1_primary" if read.is_read1 else "r2_primary"
            qname_hits[qname][key] += 1

        if read.has_tag("NH"):
            n_with_nh += 1
            nh_values[read.get_tag("NH")] += 1

    bamfh.close()

    print(f"\nTotal records: {total}")
    print(f"  Primary:       {primary}")
    print(f"  Secondary:     {secondary}")
    print(f"  Supplementary: {supplementary}")
    print(f"  Unmapped:      {unmapped}")
    print(f"\nNH tag presence: {n_with_nh}/{total} ({100*n_with_nh/max(total,1):.1f}%)")
    print(f"NH value distribution: {dict(nh_values.most_common())}")

    # Per read-name grouping analysis
    n_qnames = len(qname_hits)
    multimap_qnames = sum(1 for q in qname_hits.values()
                          if q["r1_secondary"] > 0 or q["r2_secondary"] > 0)
    asym_qnames = 0

    for qname, counts in qname_hits.items():
        r1_locs = counts["r1_primary"] + counts["r1_secondary"]
        r2_locs = counts["r2_primary"] + counts["r2_secondary"]
        if r1_locs != r2_locs:
            asym_qnames += 1

    print(f"\nRead names: {n_qnames}")
    print(f"  With secondaries: {multimap_qnames} ({100*multimap_qnames/max(n_qnames,1):.1f}%)")
    print(f"  Asymmetric (R1/R2 different #locs): {asym_qnames} ({100*asym_qnames/max(n_qnames,1):.1f}%)")

    # Show first few asymmetric examples
    if asym_qnames > 0:
        print("\n  Examples of asymmetric multimapping:")
        shown = 0
        for qname, counts in sorted(qname_hits.items()):
            r1_locs = counts["r1_primary"] + counts["r1_secondary"]
            r2_locs = counts["r2_primary"] + counts["r2_secondary"]
            if r1_locs != r2_locs and shown < 5:
                print(f"    {qname}: R1={r1_locs} locs, R2={r2_locs} locs")
                shown += 1


def analyze_hit_grouping(bam_path: Path):
    """Analyze how bam.py groups reads into hits."""
    print("\n" + "=" * 70)
    print("HIT GROUPING ANALYSIS (bam.py)")
    print("=" * 70)

    bamfh = pysam.AlignmentFile(str(bam_path), "rb")
    stats = {}

    n_groups = 0
    n_unique = 0
    n_multi = 0
    hits_per_group = Counter()
    nh_vs_hits = []  # (nh, len(hits)) pairs

    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        bamfh.fetch(until_eof=True), stats,
        include_multimap=True, skip_duplicates=True,
    ):
        n_groups += 1
        n_hits = len(hits) + len(sec_r1) + len(sec_r2)
        hits_per_group[n_hits] += 1
        nh_vs_hits.append((nh, n_hits))

        if n_hits == 1:
            n_unique += 1
        else:
            n_multi += 1

    bamfh.close()

    print(f"\nTotal read-name groups: {n_groups}")
    print(f"  Unique (1 hit):  {n_unique} ({100*n_unique/max(n_groups,1):.1f}%)")
    print(f"  Multi (>1 hits): {n_multi} ({100*n_multi/max(n_groups,1):.1f}%)")
    print(f"\nHits per group distribution:")
    for n_hits in sorted(hits_per_group):
        count = hits_per_group[n_hits]
        print(f"  {n_hits} hits: {count} groups ({100*count/max(n_groups,1):.1f}%)")

    # Check NH tag vs actual hits agreement
    nh_correct = sum(1 for nh, h in nh_vs_hits if nh == h)
    nh_under = sum(1 for nh, h in nh_vs_hits if nh < h)
    nh_over = sum(1 for nh, h in nh_vs_hits if nh > h)
    print(f"\nNH tag vs actual hits agreement:")
    print(f"  Correct (NH == len(hits)): {nh_correct}")
    print(f"  Under   (NH < len(hits)):  {nh_under} (minimap2 NH bug)")
    print(f"  Over    (NH > len(hits)):  {nh_over}")

    # Check for unpaired hits (singleton r1 or r2)
    n_unpaired_hits = 0
    n_total_hits = 0
    n_sec_r1_total = 0
    n_sec_r2_total = 0
    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        pysam.AlignmentFile(str(bam_path), "rb").fetch(until_eof=True),
        {}, include_multimap=True,
    ):
        for r1_reads, r2_reads in hits:
            n_total_hits += 1
            if not r1_reads or not r2_reads:
                n_unpaired_hits += 1
        n_sec_r1_total += len(sec_r1)
        n_sec_r2_total += len(sec_r2)

    print(f"\nTotal hits (primary): {n_total_hits}")
    print(f"  Unpaired (missing R1 or R2): {n_unpaired_hits} ({100*n_unpaired_hits/max(n_total_hits,1):.1f}%)")
    print(f"  Secondary R1 locations: {n_sec_r1_total}")
    print(f"  Secondary R2 locations: {n_sec_r2_total}")

    print(f"\nBAM stats from parse_bam_file:")
    for key in sorted(stats):
        print(f"  {key}: {stats[key]}")


def analyze_fragment_resolution(bam_path: Path, index: HulkIndex):
    """Analyze how fragments resolve against the index."""
    print("\n" + "=" * 70)
    print("FRAGMENT RESOLUTION ANALYSIS")
    print("=" * 70)

    bamfh = pysam.AlignmentFile(str(bam_path), "rb")
    stats = {}

    n_frags = 0
    n_resolved = 0
    n_intergenic = 0
    n_correct = 0  # truth transcript in candidates
    n_wrong = 0
    n_hits_dist = Counter()
    truth_tid = TRANSCRIPT_ID

    # Find t_index for the truth transcript
    truth_t_idx = None
    for _, row in index.t_df.iterrows():
        if row["t_id"] == truth_tid:
            truth_t_idx = int(row["t_index"])
            break

    if truth_t_idx is None:
        print(f"WARNING: Could not find {truth_tid} in index!")
        return

    print(f"Truth transcript: {truth_tid} (t_index={truth_t_idx})")

    for nh, hits, sec_r1, sec_r2 in parse_bam_file(
        bamfh.fetch(until_eof=True), stats,
        include_multimap=True,
    ):
        all_hits = list(hits)
        # For analysis, add sec locations as singletons
        for r1_reads in sec_r1:
            all_hits.append((r1_reads, []))
        for r2_reads in sec_r2:
            all_hits.append(([], r2_reads))
        num_hits = max(nh, len(all_hits))
        n_hits_dist[num_hits] += 1

        for r1_reads, r2_reads in all_hits:
            frag = Fragment.from_reads(r1_reads, r2_reads, sj_strand_tag="ts")
            n_frags += 1

            if not frag.exons:
                continue

            result = resolve_fragment(frag, index)
            if result is None:
                n_intergenic += 1
                continue

            n_resolved += 1
            if truth_t_idx in result.t_inds:
                n_correct += 1
            else:
                n_wrong += 1

    bamfh.close()

    print(f"\nTotal fragments: {n_frags}")
    print(f"  Resolved:    {n_resolved}")
    print(f"  Intergenic:  {n_intergenic}")
    if n_resolved > 0:
        print(f"  Correct (truth in candidates): {n_correct} ({100*n_correct/n_resolved:.1f}%)")
        print(f"  Wrong (truth NOT in candidates): {n_wrong} ({100*n_wrong/n_resolved:.1f}%)")

    print(f"\nnum_hits distribution:")
    for nh in sorted(n_hits_dist):
        count = n_hits_dist[nh]
        print(f"  {nh} hits: {count} read-name groups")


def run_pipeline_and_compare(bam_path: Path, index: HulkIndex,
                             n_simulated: int):
    """Run the full pipeline and compare to ground truth."""
    print("\n" + "=" * 70)
    print("PIPELINE COUNT COMPARISON")
    print("=" * 70)

    pipe = run_pipeline(
        bam_path, index,
        seed=42,
        sj_strand_tag="ts",
        gdna_threshold=0.0,
        em_pseudocount=0.5,
        include_multimap=True,
    )

    counts_df = pipe.counter.get_counts_df(index)
    pred = {r.transcript_id: float(r.count)
            for r in counts_df.itertuples(index=False)}

    truth = {TRANSCRIPT_ID: n_simulated}

    print(f"\nGround truth: {TRANSCRIPT_ID} = {n_simulated} fragments")
    print(f"\nPipeline predictions:")
    total_pred = 0
    for tid in sorted(pred, key=lambda t: pred[t], reverse=True):
        count = pred[tid]
        if count > 0.1:
            err = count - truth.get(tid, 0)
            total_pred += count
            print(f"  {tid:<20} {count:>8.1f}  (error: {err:>+8.1f})")

    print(f"\nTotal predicted: {total_pred:.1f}")
    print(f"Total truth:     {n_simulated}")
    print(f"Difference:      {total_pred - n_simulated:+.1f}")

    t1_pred = pred.get(TRANSCRIPT_ID, 0.0)
    print(f"\n{TRANSCRIPT_ID} accuracy: {t1_pred:.1f} / {n_simulated} "
          f"({100 * t1_pred / max(n_simulated, 1):.1f}%)")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Synthetic multimapping test for hulkrna",
    )
    parser.add_argument("--n-dups", type=int, default=2,
                        help="Number of duplicate segments (default: 2)")
    parser.add_argument("--n-frags", type=int, default=2000,
                        help="Number of fragments to simulate (default: 2000)")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="Working directory (default: temp)")
    parser.add_argument("--keep", action="store_true",
                        help="Keep working directory after completion")
    parser.add_argument("--asymmetric", action="store_true",
                        help="Duplicate only exon1 region (R2 multimaps, R1 unique)")
    args = parser.parse_args()

    if args.work_dir:
        work_dir = args.work_dir
        work_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_multimap_"))
        cleanup = not args.keep

    print(f"Working directory: {work_dir}")
    print(f"Duplicates: {args.n_dups}, Fragments: {args.n_frags}")

    try:
        # 1. Build genome with duplicated segments
        logger.info("Building genome and GTF...")
        fasta_path, gtf_path = build_genome_and_gtf(
            work_dir, n_duplicates=args.n_dups,
            asymmetric=args.asymmetric,
        )

        # 2. Simulate reads from true transcript
        logger.info("Simulating reads...")
        r1_path, r2_path = simulate_reads(
            work_dir, fasta_path, gtf_path,
            n_fragments=args.n_frags,
        )

        # 3. Align with minimap2
        logger.info("Aligning reads...")
        bam_path = align_reads(work_dir, fasta_path, r1_path, r2_path, gtf_path)

        # 4. Build index
        logger.info("Building HulkIndex...")
        index_dir = work_dir / "index"
        HulkIndex.build(fasta_path, gtf_path, index_dir, write_tsv=False)
        index = HulkIndex.load(index_dir)

        print(f"\nIndex: {index.num_transcripts} transcripts, {index.num_genes} genes")
        for _, row in index.t_df.iterrows():
            print(f"  {row['t_id']}: {row['ref']}:{row['start']}-{row['end']}")

        # 5. Analyze BAM
        analyze_bam(bam_path)

        # 6. Analyze hit grouping
        analyze_hit_grouping(bam_path)

        # 7. Analyze fragment resolution
        analyze_fragment_resolution(bam_path, index)

        # 8. Run pipeline and compare
        run_pipeline_and_compare(bam_path, index, args.n_frags)

    finally:
        if cleanup:
            import shutil
            shutil.rmtree(work_dir, ignore_errors=True)
            print(f"\nCleaned up {work_dir}")
        else:
            print(f"\nArtifacts preserved at {work_dir}")


if __name__ == "__main__":
    main()
