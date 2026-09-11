"""The synthetic index builder every structural test starts from.

`build_test_index` turns a GTF string into a loaded `TranscriptIndex`, writing an all-N FASTA of the
requested references beside it. It lives in its own module rather than in `conftest.py` because it is a
plain builder, not a fixture: a test that needs a bespoke index imports and calls it, and importing a
conftest by module name is ambiguous as soon as a sub-directory has one of its own.
"""

from __future__ import annotations

import textwrap

from rigel.index import TranscriptIndex


# ---------------------------------------------------------------------------
# Two genes on chr1, one on + and one on - strand.
# Gene A (g1): two transcripts (t1 has 3 exons, t2 has 2 exons)
# Gene B (g2): one transcript (t3 has 2 exons), negative strand
#
# Coordinates (1-based inclusive as in GTF):
#   t1: exons at 100-200, 300-400, 500-600
#   t2: exons at 100-200, 500-600
#   t3: exons at 1000-1100, 1200-1300 (neg strand)

MINI_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t1000\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
    chr1\ttest\texon\t1200\t1300\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
""")

# Expected 0-based half-open coordinates after GTF parsing:
#   t1: exons (99,200), (299,400), (499,600)
#   t2: exons (99,200), (499,600)
#   t3: exons (999,1100), (1199,1300)


def build_test_index(tmp_path_factory, gtf_text, genome_size=2000, name="idx", refs=None):
    """Build a TranscriptIndex from a GTF string (session/module-scoped helper).

    Parameters
    ----------
    tmp_path_factory : pytest.TempPathFactory
        Pytest factory for session/module-scoped temp dirs.
    gtf_text : str
        GTF content string (GENCODE-style, 1-based inclusive).
    genome_size : int
        Length of chr1 in the synthetic FASTA (default: 2000).
    name : str
        Sub-directory name for this index (keeps separate indexes apart).
    refs : dict[str, int] | None
        Optional ``{reference_name: length}`` for a MULTI-reference index. ``None`` (the default) builds
        the single ``chr1`` of ``genome_size`` bp that every existing caller expects, so this parameter
        is purely additive. The splice-graph matrix needs it: region and boundary ids must stay
        contiguous per reference, and no boundary may cross references.

    Returns
    -------
    TranscriptIndex
        Loaded index with C++ FragmentResolver ready.
    """
    import pysam

    base = tmp_path_factory.mktemp(name)
    gtf_path = base / "test.gtf"
    gtf_path.write_text(gtf_text)

    fasta_path = base / "genome.fa"
    with open(fasta_path, "w") as f:
        for ref_name, ref_len in (refs or {"chr1": genome_size}).items():
            f.write(f">{ref_name}\n")
            seq = "N" * ref_len
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(fasta_path))

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir, retain_test_structures=True)
