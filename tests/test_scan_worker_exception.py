"""A C++ exception thrown inside a scan worker thread reaches Python as an ordinary exception.

The scan runs in a child process: an exception that escapes a ``std::thread`` calls
``std::terminate``, which would abort the pytest process itself.

The throw site is the fragment buffer's uint16 ``read_length`` column. A ``D`` operation advances the
reference inside one aligned block, so a read with a 70,000 bp deletion makes an unspliced intergenic
fragment whose aligned length does not fit in 16 bits. It is the first read name in the BAM, and the
ordinary pairs after it, one read-name batch each on one worker, outnumber the input queue's capacity,
so the reader finishes only if the failing worker releases it.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

GENOME_LENGTH = 100_000
N_TRAILING_PAIRS = 64


def _build_inputs(tmp_path):
    import pysam

    from rigel.index import TranscriptIndex

    fasta, gtf = tmp_path / "g.fa", tmp_path / "a.gtf"
    fasta.write_text(">chr1\n" + "A" * GENOME_LENGTH + "\n")
    pysam.faidx(str(fasta))
    gtf.write_text(
        'chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    index_dir = tmp_path / "idx"
    TranscriptIndex.build(str(fasta), str(gtf), str(index_dir), write_tsv=False)

    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": GENOME_LENGTH}],
    }

    def _read(qname, pos, mate_pos, cigar, is_r1):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = 60
        a.flag = 0x1 | (0x40 if is_r1 else 0x80) | (0x20 if is_r1 else 0x10)
        a.cigar = cigar
        a.query_sequence = "A" * 100
        a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        a.next_reference_id = 0
        a.next_reference_start = mate_pos
        a.set_tags([("NH", 1, "i")])
        return a

    # R1 covers [1000, 71100) as one block (50M 70000D 50M); R2 lies inside it.
    reads = [
        _read("long", 1000, 71000, [(0, 50), (2, 70_000), (0, 50)], True),
        _read("long", 71000, 1000, [(0, 100)], False),
    ]
    for k in range(N_TRAILING_PAIRS):
        pos = 80_000 + 200 * k
        reads.append(_read(f"pair{k:03d}", pos, pos + 100, [(0, 100)], True))
        reads.append(_read(f"pair{k:03d}", pos + 100, pos, [(0, 100)], False))
    bam_path = tmp_path / "long.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    return bam_path, index_dir


def test_worker_exception_propagates_to_python(tmp_path):
    bam_path, index_dir = _build_inputs(tmp_path)
    script = textwrap.dedent(f"""
        from rigel.config import BamScanConfig
        from rigel.index import TranscriptIndex
        from rigel.pipeline import scan_and_buffer

        index = TranscriptIndex.load({str(index_dir)!r})
        scan = BamScanConfig(total_threads=1, read_name_batch_size=1)
        assert scan.resolved_scan_threads() == (1, 0)
        scan_and_buffer({str(bam_path)!r}, index, scan)
    """)
    proc = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, timeout=60
    )
    assert proc.returncode == 1, (
        f"the child exited with {proc.returncode} (a negative code is a signal); stderr:\n"
        f"{proc.stderr}"
    )
    assert (
        "RuntimeError: Fragment buffer column 'read_length' cannot store value 70100 as uint16"
        in proc.stderr
    ), proc.stderr
