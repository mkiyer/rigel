"""Every junction the aligner reports is read as it reported it, whatever the alignable store holds.

An alignable store also holds the junctions its aligner spliced the store's own genomic fragments at. Rigel reads none
of them. At the recommended aligner settings a genomic fragment's contiguous reading is reported beside its spliced
one, so both readings reach calibration and the EM, which weigh them on the whole likelihood (sequence, fragment
length, strand, placement and abundance), and no reading is rejected before either stage sees it.

Each index here is built from a stand-in store whose junction table names every junction these reads cross, on the
genome's reference and on a decoy the BAM carries and the FASTA lacks, with anchors at least the reads' own, so a scan
that applied the table would reject every one. Then:
- a pair through such a junction is read spliced at it, by calibration and by the EM;
- a two-hit template whose one hit crosses such a junction is two placements to calibration, as it is to the EM;
- a pair on the decoy is off the index, alone or beside a hit on the genome.
"""

from __future__ import annotations

import sys
import types

import pyarrow as pa
import pytest
from _bam_builder import R1_FWD, R2_REV, SECONDARY, record, write_bam
from _index_builder import StandInStore, write_fasta

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType

GENOME = 3000
#: The BAM's references: the genome, and a decoy the FASTA lacks (the GRCh38 analysis-set layout).
BAM_REFS = {"chr1": GENOME, "chrDecoy": GENOME}
#: One two-exon transcript, its intron [300, 400).
GTF = (
    'chr1\tt\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tt\texon\t401\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)
INTRON = (300, 400)
#: The anchor every junction row below carries, and every spliced read's shorter anchor: at or under it, a junction
#: the table names would be rejected.
ANCHOR = 50


class StoreWithJunctions(StandInStore):
    """A stand-in store that also carries the junction table a real store does: one row per (reference, intron,
    read length), with the most each anchor reached among the store's spliced genomic fragments."""

    def __init__(self, chromosomes, junctions):
        super().__init__(chromosomes)
        self._junctions = junctions

    def splice_blacklist_table(self):
        rows = [
            {
                "chrom": ref,
                "intron_start": start,
                "intron_end": end,
                "strand": ".",
                "read_length": 100,
                "count": 5,
                "max_anchor_left": ANCHOR,
                "max_anchor_right": ANCHOR,
            }
            for ref, start, end in self._junctions
        ]
        return pa.Table.from_pylist(rows)


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    """The index of `GTF`, built from a store whose junction table names `INTRON` on chr1 and on the decoy. The
    store is reached through the companion package's one entry point, ``alignable.open``, so the build reads
    whatever of it the code reads."""
    base = tmp_path_factory.mktemp("aligned_junctions")
    (base / "test.gtf").write_text(GTF)
    fasta = write_fasta(base / "genome.fa", {"chr1": GENOME})
    store_dir = base / "store"
    store_dir.mkdir()
    store = StoreWithJunctions({"chr1": GENOME}, [("chr1", *INTRON), ("chrDecoy", *INTRON)])
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(sys.modules, "alignable", types.SimpleNamespace(open=lambda path, **_: store))
        TranscriptIndex.build(
            fasta, base / "test.gtf", base / "index", write_tsv=False, alignable_zarr_path=store_dir
        )
    return TranscriptIndex.load(base / "index")


def _pair(qname, f_pos, f_cigar, r_pos, *, ref=0, nh=1, hi=None, secondary=False):
    extra = SECONDARY if secondary else 0
    return [
        record(qname, R1_FWD | extra, f_pos, f_cigar, r_pos, ref=ref, nh=nh, hi=hi, nm=0),
        record(qname, R2_REV | extra, r_pos, "100M", f_pos, ref=ref, nh=nh, hi=hi, nm=0),
    ]


#: Read 1 crosses `INTRON` with 50 bases on either side, read 2 lies in the second exon: fragment [250, 600).
SPLICED = ("50M100N50M", 250, 500)


def _scan(tmp_path, index, records):
    """The scan's stats, its accumulator payload and the splice type of every fragment it buffered for the EM."""
    bam = write_bam(tmp_path / "reads.bam", records, BAM_REFS)
    stats, _strand, buffer, payload = scan_and_buffer(bam, index, BamScanConfig(total_threads=1))
    splice_types = [int(x) for chunk in buffer.iter_chunks_consuming() for x in chunk.splice_type]
    return stats, payload, splice_types


def test_a_junction_the_store_names_is_read_spliced_by_both_stages(index, tmp_path):
    """The pair is read as the aligner reported it: spliced at the annotated intron it crossed. Calibration deposits it
    as certified RNA on that junction, and the EM scores it spliced at its sequenced junction, not as a gap its blocks
    leave open for the annotation to fill."""
    cigar, f_pos, r_pos = SPLICED
    stats, payload, splice_types = _scan(tmp_path, index, _pair("spliced", f_pos, cigar, r_pos))
    assert stats.n_sj_observed == 1
    assert stats.n_census_spliced_annot == 1 and stats.n_census_unspliced == 0
    assert stats.n_census_spliced_implicit == 0
    assert payload.qc.deposited == 1
    assert int(payload.sj_count.sum()) == 1
    assert splice_types == [SpliceType.SPLICED_ANNOT]


def test_a_template_with_a_hit_through_such_a_junction_is_two_placements_to_calibration(
    index, tmp_path
):
    """One hit crosses the junction, the other lies contiguous in intergenic sequence. Both hits are usable, so the
    EM weighs both, and calibration holds the template on both placements rather than depositing it as a certain
    fragment at the contiguous one."""
    cigar, f_pos, r_pos = SPLICED
    records = _pair("two_places", f_pos, cigar, r_pos, nh=2, hi=1)
    records += _pair("two_places", 2000, "100M", 2150, nh=2, hi=2, secondary=True)
    stats, payload, _ = _scan(tmp_path, index, records)
    assert stats.n_placed_multimapping == 1
    assert payload.qc.deposited == 0
    assert payload.qc.deferred_multiple_placements == 1


def test_a_pair_on_a_decoy_the_store_names_is_off_the_index(index, tmp_path):
    """The decoy is in the BAM and the store and not in the FASTA, so it has no place in the index: a pair there is
    withheld as off the index, and beside a hit on the genome it is a rejected hit, leaving the template unique."""
    stats, payload, _ = _scan(tmp_path, index, _pair("decoy", 1000, "100M", 1150, ref=1))
    assert stats.n_withheld_off_index == 1
    assert payload.qc.deposited == 0

    records = _pair("genome_and_decoy", 150, "100M", 300, nh=2, hi=1)
    records += _pair("genome_and_decoy", 1000, "100M", 1150, ref=1, nh=2, hi=2, secondary=True)
    stats, payload, _ = _scan(tmp_path, index, records)
    assert stats.n_placed_unique == 1 and stats.n_hits_rejected_off_index == 1
    assert payload.qc.deposited == 1
