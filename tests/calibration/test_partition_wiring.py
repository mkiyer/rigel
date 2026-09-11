"""What the scanner actually sees: the partition wiring.

``build_region_partition_arrays`` is the one function that decides which genomic partition the C++
accumulator deposits into, and :meth:`RegionArrays.from_index` reads the same frame, so the
calibration geometry and the scanner cannot address different partitions. The case exercised
throughout is an alternative TSS lying strictly INSIDE another transcript's exon: both flanks carry
the same ``exon_pos`` signature, so a signature-merged partition deletes that region bound and the
scanner never learns the terminus is there, and it is a large share of human transcript termini
rather than a corner case. A benchmark cannot substitute for these gates, because an annotation with
no mergeable adjacency produces the same partition either way and shows no effect at all — so the
assertions here are direct: the bound is emitted, the scanner is handed it, the two builders agree,
and the cache key moves with it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_region_partition_arrays
from rigel.calibration.signature import BIT_EXON_POS, coarse_type_array

from _index_builder import build_test_index


GENOME = 2000

#: t1's 5' end at 0-based 250 lies strictly INSIDE t0's exon [200, 400). Both sides are exon_pos —
#: the signature does not change there, which is precisely why a signature merge deletes this bound.
GTF_ALT_TSS = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t251\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
"""

#: The SAME annotation with t1's alternative start removed (it now shares t0's) — one region fewer,
#: and the ONLY difference is the region bound a signature merge cannot see. Under such a merge the
#: two annotations produce an identical partition and hash the same.
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
# The gate: the scanner sees the region bound a signature merge would delete.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_the_region_bound_is_signature_invisible(alt_tss_index):
    """The premise, asserted rather than assumed: this region bound carries NO signature change.

    That is what makes it invisible to a merged partition, and it is the whole reason the fixture is
    built this way. If this ever fails the fixture has stopped exercising the case and every other
    test here is vacuous.
    """
    regions = alt_tss_index.regions_df
    i = int(np.flatnonzero(regions["start"].to_numpy() == ALT_TSS_POS)[0])
    assert regions["signature"].iloc[i] == regions["signature"].iloc[i - 1] == BIT_EXON_POS
    assert int(regions["end"].iloc[i - 1]) == ALT_TSS_POS


def test_scanner_partition_carries_the_alternative_tss(alt_tss_index):
    """The position handed to ``BamScanner.set_regions`` must include the alternative TSS. Against a
    signature-merged partition this fails, because the merge deleted it."""
    positions, _offsets, _types = build_region_partition_arrays(alt_tss_index)
    assert ALT_TSS_POS in set(positions.tolist())


def test_scanner_partition_is_exactly_the_region_region_bound_array(alt_tss_index):
    """The emitted triple IS the partition frame, per reference, in ``ref_names`` order — not a
    re-derivation of it that could drift."""
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
    """``RegionArrays.from_index`` and ``build_region_partition_arrays`` must read the SAME frame.

    They are the two halves of one contract — the calibration geometry addresses the payload the
    scanner produced — and nothing downstream can detect them disagreeing except as a shape error
    far from the cause.
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
    """The invariant that replaces "no two adjacent regions share a signature".

    That was a signature-merged partition's defining property, and it is exactly what has to go: a
    large share of human transcript termini fall at a region bound where the signature does not
    change, so the merge deletes the bound that makes them visible.
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
    """Two annotations differing by exactly one signature-invisible region bound must not share a
    cache key.

    Under a merged partition they produce identical arrays and hash the same, so a payload scanned
    against one would load silently against the other. That is the stale-cache failure
    ``partition_hash`` exists to make impossible.
    """
    assert len(alt_tss_index.regions_df) == len(no_alt_tss_index.regions_df) + 1
    assert alt_tss_index.partition_hash != no_alt_tss_index.partition_hash
