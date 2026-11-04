"""A C++ exception thrown inside a scan worker thread reaches Python as an ordinary exception.

The scan runs in a child process: an exception that escapes a ``std::thread`` calls
``std::terminate``, which would abort the pytest process itself.

The throw site is the fragment buffer's uint16 ``read_length`` column. A ``D`` operation advances the
reference inside one aligned block, so a read with a 70,000 bp deletion makes an unspliced intergenic
fragment whose aligned length does not fit in 16 bits. It is the first read name in the BAM, and the
ordinary pairs after it, one read-name batch each, outnumber the input queue's capacity. On one worker the
reader finishes only if the failing worker releases it; on several, the others drain the queue while one
holds the error, and the scan must still stop and raise it. Each run records the worker count the native
scan was handed, so a scan quietly run on one worker cannot stand in for several.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

import pytest
from _bam_builder import R1_FWD, R2_REV, record, write_bam
from _index_builder import build_index

GENOME_LENGTH = 100_000
N_TRAILING_PAIRS = 64


def _build_inputs(tmp_path):
    gtf = (
        'chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\ttest\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    index_dir = build_index(tmp_path, gtf, {"chr1": GENOME_LENGTH}).index_dir

    # R1 covers [1000, 71100) as one block (50M 70000D 50M); R2 lies inside it.
    reads = [
        record("long", R1_FWD, 1000, "50M70000D50M", 71000),
        record("long", R2_REV, 71000, "100M", 1000),
    ]
    for k in range(N_TRAILING_PAIRS):
        pos = 80_000 + 200 * k
        reads += [
            record(f"pair{k:03d}", R1_FWD, pos, "100M", pos + 100),
            record(f"pair{k:03d}", R2_REV, pos + 100, "100M", pos),
        ]
    bam_path = write_bam(tmp_path / "long.bam", reads, {"chr1": GENOME_LENGTH})
    return bam_path, index_dir


@pytest.mark.parametrize("n_workers", [1, 4])
def test_worker_exception_propagates_to_python(tmp_path, n_workers):
    bam_path, index_dir = _build_inputs(tmp_path)
    script = textwrap.dedent(f"""
        import sys

        import rigel.pipeline as pipeline
        from rigel.config import BamScanConfig
        from rigel.index import TranscriptIndex

        native_scanner = pipeline._native_scanner

        class Recording:
            def __init__(self, scanner):
                self.scanner = scanner

            def scan(self, *args, **kwargs):
                print(f"native scan n_workers={{kwargs['n_workers']}}", file=sys.stderr, flush=True)
                return self.scanner.scan(*args, **kwargs)

        pipeline._native_scanner = lambda index, scan: Recording(native_scanner(index, scan))
        index = TranscriptIndex.load({str(index_dir)!r})
        # a limit above the 70,100-bp molecule, so the pair is a hit and the buffer, not the length rule, refuses it
        scan = BamScanConfig(
            total_threads={n_workers}, bgzf_threads=0, read_name_batch_size=1, max_frag_length=100_000
        )
        assert scan.resolved_scan_threads() == ({n_workers}, 0)
        pipeline.scan_and_buffer({str(bam_path)!r}, index, scan)
    """)
    proc = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, timeout=60
    )
    assert proc.returncode == 1, (
        f"the child exited with {proc.returncode} (a negative code is a signal); stderr:\n"
        f"{proc.stderr}"
    )
    assert f"native scan n_workers={n_workers}" in proc.stderr, proc.stderr
    assert (
        "RuntimeError: Fragment buffer column 'read_length' cannot store value 70100 as uint16"
        in proc.stderr
    ), proc.stderr
