"""A two-gene fixture for the BAM passes and the output directory: genes A and B, one two-exon transcript each, with
identical layouts 1,000 bases apart, and three spliced read pairs placed wholly at either one."""

from __future__ import annotations


def two_gene_index(root):
    """Genes A and B, one two-exon transcript each, identical layouts 1,000 bases apart."""
    import pysam

    from rigel.index import TranscriptIndex

    fasta, gtf = root / "genome.fa", root / "genes.gtf"
    fasta.write_text(">chr1\n" + "ACGT" * 625 + "\n")
    pysam.faidx(str(fasta))
    rows = [
        f'chr1\tfixture\texon\t{s + off}\t{e + off}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}";\n'
        for g, off in (("A", 0), ("B", 1000))
        for s, e in ((101, 300), (501, 700))
    ]
    gtf.write_text("".join(rows))
    TranscriptIndex.build(fasta, gtf, root / "index", write_tsv=False)
    return TranscriptIndex.load(root / "index")


def three_pairs(path, offset):
    """Three spliced pairs placed wholly at A (offset 0) or wholly at B (offset 1000): the same names, flags, record
    counts and placement outcomes either way."""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 2500}]}
    mates = [(250, 99, "50M200N50M", 100, 555), (555, 147, "80M", 80, 250)]
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for pair in range(3):
            for start, flag, cigar, length, mate in mates:
                read = pysam.AlignedSegment(out.header)
                read.query_name = f"pair{pair}"
                read.reference_id = read.next_reference_id = 0
                read.reference_start = start + offset
                read.next_reference_start = mate + offset
                read.mapping_quality, read.flag, read.cigarstring = 60, flag, cigar
                read.query_sequence = "A" * length
                read.set_tags([("NH", 1), ("NM", 0), ("XS", "+", "A")])
                out.write(read)
