"""Diagnose why INTRONIC count == 0 in real-data calibration.

Synthetic test: build a tiny index with a 2-exon transcript that has
a long intron, simulate a paired-end fragment landing entirely in the
intron, run the C++ scanner, and check what tx_bp_pos / exon_bp_pos /
splice_type / category come out as.
"""
import numpy as np
import pysam
from pathlib import Path
import tempfile
import subprocess

from rigel.index import TranscriptIndex
from rigel.calibration._categorize import (
    categorize_chunk, FragmentCategory, StrandLabel, N_STRAND_LABELS,
)
from rigel.splice import SPLICE_UNSPLICED


def main():
    tmp = Path(tempfile.mkdtemp(prefix="rigel_intronic_bug_"))
    print("tmp:", tmp)

    # Tiny FASTA: 1 chromosome, 100kb of As
    fasta = tmp / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 100000 + "\n")
    pysam.faidx(str(fasta))

    # Tiny GTF: 1 gene, 1 transcript, 2 exons with a 5kb intron
    gtf = tmp / "ann.gtf"
    gtf.write_text(
        "\n".join([
            'chr1\trigel\tgene\t1000\t8000\t.\t+\t.\tgene_id "g1"; gene_name "G1";',
            'chr1\trigel\ttranscript\t1000\t8000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
            'chr1\trigel\texon\t1000\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";',
            'chr1\trigel\texon\t7000\t8000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";',
        ]) + "\n"
    )

    idx_dir = tmp / "index"
    TranscriptIndex.build(str(fasta), str(gtf), str(idx_dir), write_tsv=False)

    # Simulate a tiny BAM with paired-end reads landing in the INTRON
    # (positions 3000-3100 R1, 3200-3300 R2; both inside intron 2000-7000)
    sam = tmp / "reads.sam"
    sam_lines = [
        "@HD\tVN:1.6\tSO:queryname",
        "@SQ\tSN:chr1\tLN:100000",
    ]
    # Pair 1: both reads INTRONIC (no exon overlap, but inside transcript span)
    sam_lines += [
        "r_intron\t99\tchr1\t3001\t60\t100M\t=\t3201\t300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
        "r_intron\t147\tchr1\t3201\t60\t100M\t=\t3001\t-300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
    ]
    # Pair 2: both reads INTERGENIC (after transcript end at 8000)
    sam_lines += [
        "r_inter\t99\tchr1\t9001\t60\t100M\t=\t9201\t300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
        "r_inter\t147\tchr1\t9201\t60\t100M\t=\t9001\t-300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
    ]
    # Pair 3: both reads EXONIC inside exon 1
    sam_lines += [
        "r_exon\t99\tchr1\t1100\t60\t100M\t=\t1300\t300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
        "r_exon\t147\tchr1\t1300\t60\t100M\t=\t1100\t-300\t" + "A"*100 + "\t" + "I"*100 + "\tNM:i:0\tNH:i:1",
    ]
    sam.write_text("\n".join(sam_lines) + "\n")

    bam = tmp / "reads.bam"
    pysam.view("-bS", "-o", str(bam), str(sam), catch_stdout=False)

    # Run the C++ scanner
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer

    index = TranscriptIndex.load(str(idx_dir), retain_test_structures=True)
    # Dump the overlap index intervals
    print("Index intervals:")
    import pyarrow.feather as pf
    iv_df = pf.read_feather(str(idx_dir / "intervals.feather"))
    print(iv_df.to_string())
    cfg = BamScanConfig()
    stats, strand_model, fl_models, buffer = scan_and_buffer(
        bam_path=str(bam), index=index, scan=cfg,
    )

    print(f"\nbuffer.total_fragments = {buffer.total_fragments}")
    for chunk in buffer.iter_chunks():
        print(f"\nchunk.size = {chunk.size}")
        for i in range(chunk.size):
            print(
                f"  [{i}] splice_type={chunk.splice_type[i]} "
                f"num_hits={chunk.num_hits[i]} "
                f"genomic_start={chunk.genomic_start[i]} "
                f"genomic_footprint={chunk.genomic_footprint[i]} "
                f"exon_bp_pos={chunk.exon_bp_pos[i]} "
                f"exon_bp_neg={chunk.exon_bp_neg[i]} "
                f"tx_bp_pos={chunk.tx_bp_pos[i]} "
                f"tx_bp_neg={chunk.tx_bp_neg[i]} "
                f"chimera={chunk.chimera_type[i]}"
            )
        cc = categorize_chunk(chunk)
        for i in range(chunk.size):
            cat = cc.category[i]
            cat_name = (
                FragmentCategory(cat).name if cat != 255 else "FILTERED"
            )
            print(f"  [{i}] category={cat_name} keep={cc.keep[i]}")


if __name__ == "__main__":
    main()
