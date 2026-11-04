"""Synthetic alignments for a test that scans a BAM it writes: the read-pair flags, one record builder, and a
name-sorted BAM writer.

A record's query is ``A`` repeated to the length its CIGAR consumes, at full quality, so a test states only what the
scanner reads: the position, the CIGAR, the mate and the tags.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import pysam

#: Paired, properly paired, read 1 forward with its mate reverse; and read 2 reverse with its mate forward.
R1_FWD = 0x1 | 0x2 | 0x40 | 0x20
R2_REV = 0x1 | 0x2 | 0x80 | 0x10
#: The other orientation: read 1 reverse, read 2 forward.
R1_REV = 0x1 | 0x2 | 0x40 | 0x10
R2_FWD = 0x1 | 0x2 | 0x80 | 0x20
SECONDARY = 0x100


def record(
    qname: str,
    flag: int,
    pos: int,
    cigar,
    mate_pos: int,
    *,
    ref: int = 0,
    mapq: int = 60,
    nh: int | None = 1,
    hi: int | None = None,
    nm: int | None = None,
    xs: str | None = None,
    ts: str | None = None,
    sa: str | None = None,
) -> pysam.AlignedSegment:
    """One aligned record on reference *ref* with its mate on the same reference. *cigar* is a string or pysam's
    ``(op, length)`` tuples. A tag left ``None`` is not written."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = flag
    a.reference_id = ref
    a.reference_start = pos
    a.mapping_quality = mapq
    if isinstance(cigar, str):
        a.cigarstring = cigar
    else:
        a.cigar = cigar
    n = a.infer_query_length()
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = ref
    a.next_reference_start = mate_pos
    tags = [(tag, v, "i") for tag, v in (("NH", nh), ("HI", hi), ("NM", nm)) if v is not None]
    tags += [(tag, v, "A") for tag, v in (("XS", xs), ("ts", ts)) if v is not None]
    if sa is not None:
        tags.append(("SA", sa, "Z"))
    a.set_tags(tags)
    return a


def write_bam(path, records: Iterable[pysam.AlignedSegment], refs: Mapping[str, int]) -> str:
    """*records* as a name-sorted BAM over *refs* (``{name: length}``, in reference-id order); returns the path."""
    path = str(path)
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": name, "LN": int(length)} for name, length in refs.items()],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for r in records:
            out.write(r)
    pysam.sort("-n", "-o", path, path)
    return path
