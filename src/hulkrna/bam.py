"""
hulkrna.bam — BAM read-pair iterator and CIGAR parsing.

Requires a **name-sorted** (``samtools sort -n``) or **collated**
(``samtools collate``) BAM file so that all alignments for a read
name arrive consecutively.  This eliminates the unbounded-memory
pair buffer required for coordinate-sorted input.

Read-name groups are processed atomically: all alignments for a
query name are collected, categorised, grouped into *hits*, and
yielded as ``(nh, hits)``.

A **hit** represents one mapping location for the fragment and
comprises all BAM records (primary + supplementary) for that location.
Each hit is a ``(r1_reads, r2_reads)`` tuple of lists.

Hit grouping strategy
---------------------
* **HI tag present** (STAR, some minimap2 configs): records are
  grouped by their HI (Hit Index) value.
* **HI tag absent** (minimap2 default): the primary pair (primary +
  supplementary records) is returned as a single hit.  Secondary
  R1 and R2 locations are returned separately so the pipeline can
  perform transcript-aware pairing (resolve-then-pair) rather than
  blind cross-product pairing.  This avoids creating false chimeric
  pairs between paralogs on different chromosomes.

For unique mappers (NH = 1) there is exactly one hit.  For
multimappers (NH > 1, when ``include_multimap=True``) there may
be several.
"""

import logging
from itertools import groupby
from pathlib import Path

import pysam

from .types import Strand

logger = logging.getLogger(__name__)

# N: intron (reference skip)
# M, D, =, X: advance reference within exon
CIGAR_OPS_ADVANCE_REF = {
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
    pysam.CDEL,
}


def detect_sj_strand_tag(
    bam_path: str | Path,
    candidates: tuple[str, ...] = ("XS", "ts"),
    max_spliced_reads: int = 1000,
) -> tuple[str, ...]:
    """Scan a BAM file to detect which splice-junction strand tags are present.

    Examines reads containing CIGAR reference-skip (``N``) operations
    and checks whether they carry any of the *candidates* tags.

    Parameters
    ----------
    bam_path : str or Path
        Path to the BAM file.
    candidates : tuple of str
        Tag names to look for, in priority order (default
        ``("XS", "ts")``).
    max_spliced_reads : int
        Stop scanning after this many spliced reads (default 1000).

    Returns
    -------
    tuple[str, ...]
        Detected tag names, preserving the order of *candidates*.
        Empty tuple if none of the candidates were found.
    """
    detected: set[str] = set()
    n_spliced = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.cigartuples is None:
                continue

            # Only inspect spliced reads (those with CIGAR N ops)
            has_splice = any(
                op == pysam.CREF_SKIP for op, _ in read.cigartuples
            )
            if not has_splice:
                continue

            n_spliced += 1
            for tag in candidates:
                if read.has_tag(tag):
                    detected.add(tag)

            if n_spliced >= max_spliced_reads:
                break

    # Preserve candidate priority order
    result = tuple(t for t in candidates if t in detected)

    if result:
        logger.info(
            "Auto-detected SJ strand tag(s): %s (scanned %d spliced reads)",
            ", ".join(result), n_spliced,
        )
    else:
        # TODO: Implement strand inference from reference annotation when
        # no splice-junction strand tag is present in the BAM file.
        # For now all splice junctions will receive Strand.NONE, which
        # disables strand-aware counting.
        logger.warning(
            "No SJ strand tags (%s) found after scanning %d spliced reads. "
            "All splice junctions will have strand=NONE.",
            ", ".join(candidates), n_spliced,
        )

    return result


def parse_read(read, sj_strand_tag: str | tuple[str, ...] = "XS"):
    """Parse a single alignment's CIGAR into exon blocks and splice junctions.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A single aligned read.
    sj_strand_tag : str or tuple of str
        BAM tag(s) to read splice-junction strand from.  When a
        tuple is given the tags are checked in order and the first
        one present on the read is used.  Examples: ``"XS"`` (STAR),
        ``"ts"`` (minimap2), ``("XS", "ts")`` (check both).

    Returns
    -------
    tuple
        ``(ref, ref_strand, exons, sjs)`` where *exons* is a list of
        ``(start, end)`` tuples and *sjs* is a list of
        ``(start, end, sj_strand)`` tuples.
    """
    exons = []
    sjs = []
    ref = read.reference_name
    start = read.reference_start
    pos = start
    ref_strand = Strand.from_is_reverse(read.is_reverse)

    # Determine SJ strand from the configured tag(s)
    _sj_strand = Strand.NONE
    tags = (sj_strand_tag,) if isinstance(sj_strand_tag, str) else sj_strand_tag
    for tag in tags:
        if read.has_tag(tag):
            _sj_strand = Strand.from_str(read.get_tag(tag))
            # minimap2's 'ts' tag is alignment-relative: the sign
            # indicates whether the donor-acceptor motif is on the
            # same (+) or opposite (-) strand as the alignment.
            # Convert to reference-relative by flipping for
            # reverse-mapped reads.
            if tag == "ts" and read.is_reverse:
                _sj_strand = _sj_strand.opposite()
            break

    for (cig_op, cig_len) in read.cigartuples:
        # handle reference skip (intron)
        if cig_op == pysam.CREF_SKIP:
            if pos > start:
                # add the exon from start to pos
                exons.append((start, pos))
            # add the splice junction from pos to (pos + cig_len)
            sjs.append((pos, pos + cig_len, _sj_strand))
            start = pos + cig_len
            pos = start
        elif cig_op in CIGAR_OPS_ADVANCE_REF:
            # advance reference position
            pos += cig_len

    # add the last exon if it exists
    if pos > start:
        exons.append((start, pos))
    return ref, ref_strand, exons, sjs


def _pair_reads(r1_reads, r2_reads):
    """Pair read-1 and read-2 alignments by mate position.

    Matching uses ``(reference_id, reference_start)`` of each read
    against the mate position encoded in ``(next_reference_id,
    next_reference_start)`` of its partner.

    Reads whose mate is unmapped (``mate_is_unmapped``) are yielded
    as singletons: ``(read, None)`` or ``(None, read)``.

    Parameters
    ----------
    r1_reads : list[pysam.AlignedSegment]
        Read-1 alignments.
    r2_reads : list[pysam.AlignedSegment]
        Read-2 alignments.

    Returns
    -------
    list[tuple]
        Each element is ``(r1, r2)``, ``(r1, None)``, or ``(None, r2)``.
    """
    pairs = []

    # Separate mate-unmapped singletons from pairable reads
    paired_r1 = []
    for r1 in r1_reads:
        if r1.mate_is_unmapped:
            pairs.append((r1, None))
        else:
            paired_r1.append(r1)

    paired_r2 = []
    for r2 in r2_reads:
        if r2.mate_is_unmapped:
            pairs.append((None, r2))
        else:
            paired_r2.append(r2)

    # Index r2 reads by their genomic position for O(n) matching
    r2_by_pos: dict[tuple[int, int], list[int]] = {}
    for i, r2 in enumerate(paired_r2):
        key = (r2.reference_id, r2.reference_start)
        r2_by_pos.setdefault(key, []).append(i)

    used_r2: set[int] = set()
    for r1 in paired_r1:
        mate_key = (r1.next_reference_id, r1.next_reference_start)
        candidates = r2_by_pos.get(mate_key, [])
        matched = False
        for idx in candidates:
            if idx in used_r2:
                continue
            r2 = paired_r2[idx]
            # Verify reciprocal mate pointers
            if (r2.next_reference_id == r1.reference_id
                    and r2.next_reference_start == r1.reference_start):
                pairs.append((r1, r2))
                used_r2.add(idx)
                matched = True
                break
        if not matched:
            # r1 with no matching r2 → singleton
            pairs.append((r1, None))

    # Remaining unmatched r2 → singletons
    for i, r2 in enumerate(paired_r2):
        if i not in used_r2:
            pairs.append((None, r2))

    return pairs


# ---------------------------------------------------------------------------
# Hit grouping
# ---------------------------------------------------------------------------


def _group_records_by_hit(usable_records):
    """Group BAM records for one query name into alignment hits.

    Each hit is a ``(r1_reads, r2_reads)`` tuple where both elements
    are lists of ``pysam.AlignedSegment`` (primary + supplementary
    records for that mapping location).

    Strategy
    --------
    * **HI tag present** — group records by their HI value.  All
      hits are returned in ``hits``; secondary lists are empty.
    * **HI tag absent** — the primary pair (primary + supplementary)
      is the sole entry in ``hits``.  Mate-unmapped singletons are
      also appended.  Pairable secondary R1/R2 locations are
      returned separately for transcript-aware pairing by the
      pipeline.

    Parameters
    ----------
    usable_records : list[pysam.AlignedSegment]
        Filtered records (no QC-fail, no unmapped, no skipped
        duplicates) for one query name.

    Returns
    -------
    tuple[list[tuple[list, list]], list[list], list[list]]
        ``(hits, sec_r1_locs, sec_r2_locs)`` where *hits* contains
        the primary pair (and mate-unmapped singletons for the
        no-HI path, or all HI-grouped hits for the HI path),
        and *sec_r1_locs* / *sec_r2_locs* are lists of secondary
        R1/R2 locations (each a single-element list of BAM records)
        for transcript-aware pairing.
    """
    # Detect whether HI tags are available
    has_hi = any(r.has_tag('HI') for r in usable_records)

    if has_hi:
        groups: dict[int, tuple[list, list]] = {}
        for r in usable_records:
            hi = r.get_tag('HI') if r.has_tag('HI') else 0
            if hi not in groups:
                groups[hi] = ([], [])
            if r.is_read1:
                groups[hi][0].append(r)
            elif r.is_read2:
                groups[hi][1].append(r)
        return [groups[hi] for hi in sorted(groups)], [], []

    # Fallback: separate primary from secondary locations
    #
    # minimap2 secondaries have RNEXT/PNEXT that point back to the
    # PRIMARY mate, not to a reciprocal secondary.  Instead of blind
    # cross-product pairing (which creates false chimeric pairs between
    # paralogs on different chromosomes), we return secondary R1/R2
    # locations separately so the pipeline can perform transcript-aware
    # pairing (resolve each half independently, then pair by transcript
    # set intersection).
    #
    # Primary + supplementary records form the primary hit (always
    # included).  Secondaries whose mate is explicitly unmapped are
    # emitted as singleton hits.  All other secondaries go into the
    # sec_r1_locs / sec_r2_locs lists for resolve-then-pair.

    primary_r1: list = []     # primary + supplementary R1
    primary_r2: list = []     # primary + supplementary R2
    sec_r1_locs: list[list] = []  # secondary R1 locations for pairing
    sec_r2_locs: list[list] = []  # secondary R2 locations for pairing
    singleton_hits: list[tuple[list, list]] = []

    for r in usable_records:
        if r.is_supplementary:
            # Supplementary stays with primary (hit 0)
            if r.is_read1:
                primary_r1.append(r)
            else:
                primary_r2.append(r)
        elif r.is_secondary:
            # Secondaries whose mate is explicitly unmapped cannot
            # participate in pairing — emit as singleton hits.
            if r.mate_is_unmapped:
                if r.is_read1:
                    singleton_hits.append(([r], []))
                else:
                    singleton_hits.append(([], [r]))
            elif r.is_read1:
                sec_r1_locs.append([r])
            else:
                sec_r2_locs.append([r])
        else:
            if r.is_read1:
                primary_r1.append(r)
            else:
                primary_r2.append(r)

    # Primary hit: always included
    hits: list[tuple[list, list]] = []
    if primary_r1 or primary_r2:
        hits.append((primary_r1, primary_r2))

    # Mate-unmapped singletons
    hits.extend(singleton_hits)

    return hits, sec_r1_locs, sec_r2_locs


def parse_bam_file(
    bam_iter,
    stats,
    *,
    skip_duplicates=True,
    include_multimap=False,
):
    """Parse a name-sorted / collated BAM, yielding hit groups.

    All alignments for a read name (query name) must arrive
    consecutively — guaranteed by ``samtools sort -n`` or
    ``samtools collate``.

    Records are grouped into **hits** (mapping locations).  Each hit
    contains all BAM records — primary, secondary, and supplementary —
    for that location.  Secondary records are retained (not filtered)
    so that multimappers are correctly represented.  Supplementary
    records are retained so that chimeric / split alignments can be
    merged into a single ``Fragment``.

    Parameters
    ----------
    bam_iter : iterator
        Iterator over ``pysam.AlignedSegment`` objects.
    stats : dict
        Mutable dict populated with filtering statistics.
    skip_duplicates : bool
        If True (default), discard reads marked as PCR/optical
        duplicates.
    include_multimap : bool
        If True, yield multimapping read-name groups (NH > 1).
        If False (default), skip them (but still count in stats).

    Yields
    ------
    tuple[int, list[tuple[list, list]], list[list], list[list]]
        ``(nh, hits, sec_r1_locs, sec_r2_locs)`` where *nh* is the NH
        tag value, *hits* is a list of ``(r1_reads, r2_reads)`` tuples
        for the primary pair (and mate-unmapped singletons / HI-grouped
        hits), and *sec_r1_locs* / *sec_r2_locs* are lists of secondary
        R1/R2 locations for transcript-aware pairing by the pipeline.
    """
    stat_keys = [
        'total', 'qc_fail', 'unmapped', 'secondary',
        'supplementary', 'duplicate',
        'n_read_names', 'unique', 'multimapping',
        'proper_pair', 'improper_pair', 'mate_unmapped',
    ]
    for key in stat_keys:
        stats.setdefault(key, 0)

    for _qname, group_iter in groupby(
        bam_iter, key=lambda r: r.query_name
    ):
        usable: list = []
        nh = 1

        for read in group_iter:
            stats['total'] += 1

            # Skip unusable alignments
            if read.is_qcfail:
                stats['qc_fail'] += 1
                continue
            if read.is_unmapped:
                stats['unmapped'] += 1
                continue

            # Duplicate handling
            if read.is_duplicate:
                stats['duplicate'] += 1
                if skip_duplicates:
                    continue

            # Enforce paired-end
            if not read.is_paired:
                raise ValueError(
                    f"Input BAM must be paired-end, but read "
                    f"{read.query_name!r} is unpaired"
                )

            # Count secondary / supplementary for stats
            if read.is_secondary:
                stats['secondary'] += 1
            if read.is_supplementary:
                stats['supplementary'] += 1

            # Read NH from first usable alignment in the group
            if read.has_tag('NH'):
                nh = read.get_tag('NH')

            usable.append(read)

        # Skip empty groups (everything filtered)
        if not usable:
            continue

        stats['n_read_names'] += 1

        # Detect multimapper: NH tag or presence of secondary records
        has_secondary = any(r.is_secondary for r in usable)
        is_multimap = nh > 1 or has_secondary

        if is_multimap:
            stats['multimapping'] += 1
            if not include_multimap:
                continue
        else:
            stats['unique'] += 1

        # Group records into hits + separate secondary locations
        hits, sec_r1, sec_r2 = _group_records_by_hit(usable)

        # Track hit-level stats (primary hits only;
        # secondary locations are paired later by the pipeline)
        for r1_reads, r2_reads in hits:
            if not r1_reads or not r2_reads:
                stats['mate_unmapped'] += 1
            else:
                # Check proper-pair flag on the primary R1
                r1_primary = [
                    r for r in r1_reads
                    if not r.is_supplementary and not r.is_secondary
                ]
                if r1_primary and r1_primary[0].is_proper_pair:
                    stats['proper_pair'] += 1
                else:
                    stats['improper_pair'] += 1

        yield nh, hits, sec_r1, sec_r2
