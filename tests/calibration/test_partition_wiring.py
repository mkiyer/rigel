"""What the scanner sees: the partition wiring.

``build_region_partition_arrays`` decides which genomic partition the C++ accumulator deposits into, and
:meth:`RegionArrays.from_index` reads the same frame, so the calibration geometry and the scanner cannot address
different partitions. The case exercised throughout is an alternative TSS strictly inside another transcript's exon:
both flanks carry the same ``exon_pos`` signature, so a signature-merged partition would delete that region bound and
hide the terminus from the scanner, and such termini are a large share of human transcript termini. A benchmark cannot
replace these tests, since an annotation with no mergeable adjacency gives the same partition either way; so they assert
directly that the bound is emitted, the scanner is handed it, the two builders agree, and the cache key moves with it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_region_partition_arrays
from rigel.calibration.signature import BIT_EXON_POS, coarse_type_array

from _index_builder import build_test_index


GENOME = 2000

#: t1's 5' end at 0-based 250 lies strictly inside t0's exon [200, 400). Both sides are exon_pos, so the signature does
#: not change there, which is why a signature merge would delete this bound.
GTF_ALT_TSS = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t251\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
"""

#: The same annotation with t1's alternative start removed (it now shares t0's): one region fewer, the only difference
#: being the region bound a signature merge cannot see. Under such a merge the two annotations would produce an
#: identical partition and hash the same.
GTF_NO_ALT_TSS = GTF_ALT_TSS.replace(
    'chr1\ttest\texon\t251\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
    'chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
)

ALT_TSS_POS = 250


@pytest.fixture(scope="module")
def alt_tss_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF_ALT_TSS, genome_size=GENOME, name="w1b_alt")


@pytest.fixture(scope="module")
def no_alt_tss_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF_NO_ALT_TSS, genome_size=GENOME, name="w1b_noalt")


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# The scanner sees the region bound a signature merge would delete.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_the_region_bound_is_signature_invisible(alt_tss_index):
    """This region bound carries no signature change, which is what makes it invisible to a merged partition; if this
    fails the fixture no longer exercises the case and the other tests here are vacuous.
    """
    regions = alt_tss_index.regions_df
    i = int(np.flatnonzero(regions["start"].to_numpy() == ALT_TSS_POS)[0])
    assert regions["signature"].iloc[i] == regions["signature"].iloc[i - 1] == BIT_EXON_POS
    assert int(regions["end"].iloc[i - 1]) == ALT_TSS_POS


def test_scanner_partition_carries_the_alternative_tss(alt_tss_index):
    """The positions handed to ``BamScanner.set_regions`` must include the alternative TSS, which a signature-merged
    partition would delete."""
    positions, _offsets, _types = build_region_partition_arrays(alt_tss_index)
    assert ALT_TSS_POS in set(positions.tolist())


def test_scanner_partition_is_exactly_the_region_region_bound_array(alt_tss_index):
    """The emitted triple is the partition frame, per reference, in ``ref_names`` order, not a re-derivation that could
    drift."""
    positions, offsets, types = build_region_partition_arrays(alt_tss_index)
    regions = alt_tss_index.regions_df

    n_refs = len(alt_tss_index.ref_names)
    assert offsets.shape == (n_refs + 1,)
    assert int(offsets[-1]) == positions.shape[0]
    # k regions per reference ⇒ k+1 region-bound positions; region_types is 1:1 with regions.
    assert positions.shape[0] == len(regions) + n_refs
    assert types.shape[0] == len(regions)

    for i, ref in enumerate(alt_tss_index.ref_names):
        sub = regions[regions["ref_name"] == ref]
        want = np.append(sub["start"].to_numpy(np.int64), sub["end"].to_numpy(np.int64)[-1])
        np.testing.assert_array_equal(positions[offsets[i] : offsets[i + 1]], want)
    np.testing.assert_array_equal(types, coarse_type_array(regions["signature"].to_numpy()))


def test_region_geometry_follows_the_scanner_partition(alt_tss_index):
    """``RegionArrays.from_index`` and ``build_region_partition_arrays`` must read the same frame: the calibration
    geometry addresses the payload the scanner produced, and a disagreement would surface only as a shape error far from
    the cause.
    """
    ra = RegionArrays.from_index(alt_tss_index)
    regions = alt_tss_index.regions_df
    assert ra.n_regions == len(regions)
    np.testing.assert_array_equal(ra.start, regions["start"].to_numpy(np.int64))
    np.testing.assert_array_equal(ra.end, regions["end"].to_numpy(np.int64))
    np.testing.assert_array_equal(ra.signature, regions["signature"].to_numpy(np.uint8))

    positions, offsets, _types = build_region_partition_arrays(alt_tss_index)
    # k+1 region-bound positions per reference against k regions: the payload's own relation.
    np.testing.assert_array_equal(
        np.diff(offsets.astype(np.int64)) - 1, np.diff(ra.ref_offsets.astype(np.int64))
    )


def test_adjacent_regions_may_share_a_signature(alt_tss_index):
    """Adjacent regions may share a signature: a large share of human transcript termini fall at a region bound where
    the signature does not change, so merging equal-signature neighbours would delete the bounds that make them visible.
    """
    ra = RegionArrays.from_index(alt_tss_index)
    i = int(np.flatnonzero(ra.start == ALT_TSS_POS)[0])
    assert ra.signature[i] == ra.signature[i - 1] == BIT_EXON_POS


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# The cache key must track the partition the scanner actually sees.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_partition_hash_distinguishes_a_signature_invisible_region_bound(
    alt_tss_index, no_alt_tss_index
):
    """Two annotations differing by exactly one signature-invisible region bound must not share a cache key; under a
    merged partition they would hash the same, and a payload scanned against one would load silently against the other.
    """
    assert len(alt_tss_index.regions_df) == len(no_alt_tss_index.regions_df) + 1
    assert alt_tss_index.partition_hash != no_alt_tss_index.partition_hash
