"""WHAT THE SCANNER SEES — the v8 partition wiring (accumulator plan W1b).

     (W1b) · Graph:

``build_region_partition_arrays`` is the one function that decides which genomic partition the C++
accumulator deposits into: the v8 region bound array. :meth:`RegionArrays.from_index` reads the same
frame, so the calibration geometry and the scanner cannot address different partitions.

⚠ **This module exists because the standard 32-condition bench CANNOT see that change.** Measured
2026-07-29: the ``ambig_dense_10mb`` annotation has ZERO mergeable adjacencies, so its v8 region
partition is byte-identical to its v7 region partition — 1,698 == 1,698, equal row-for-row on
``(ref_name, start, end, length, signature)``. The partition effect is **+25.0 %** on
``quick_3to1_5mb`` and **+38.7 %** on the human annotation, and **0.0 %** on the bench. So the arm
needs a direct assertion that the scanner now sees the v8 region_bound set, and this is it.

The case throughout is **TRAPS: the-source-does-not-cite-docs — an alternative TSS strictly interior to another transcript's exon**.
Both flanks are ``exon_pos``, so the v7 merge deletes the region_bound and the scanner never learns the
terminus is there. That is the defect v8 fixes, and **53.4 %** of human transcript termini sit in it
(232,451 of 435,291). ⚠ The 59.5 % this file used to quote was computed under the buggy
``~is_synthetic & ~is_nrna`` transcript filter, i.e. it was reading the annotation the way the bug did.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_region_partition_arrays
from rigel.calibration.signature import BIT_EXON_POS, coarse_type_array

from conftest import build_test_index


GENOME = 2000

#: t1's 5' end at 0-based 250 lies strictly INSIDE t0's exon [200, 400). Both sides are exon_pos —
#: the signature does not change there, which is precisely why the retired v7 merge deleted this region_bound.
GTF_ALT_TSS = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t251\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
"""

#: The SAME annotation with t1's alternative start removed (it now shares t0's) — one region fewer,
#: and the ONLY difference is the region_bound a signature merge cannot see. Measured under v7: the two
#: annotations produced an identical partition and hashed the same (``645a0dd8aa560236``).
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
# The arm's own gate: the scanner sees the region_bound the merge deletes.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_the_region_bound_is_signature_invisible(alt_tss_index):
    """The premise, asserted rather than assumed: this region_bound carries NO signature change.

    That is what made it invisible to a merged partition, and it is what makes the fixture a TRAPS: the-source-does-not-cite-docs case.
    If this ever fails the fixture has stopped exercising TRAPS: the-source-does-not-cite-docs and every other test here is vacuous.
    """
    regions = alt_tss_index.regions_df
    i = int(np.flatnonzero(regions["start"].to_numpy() == ALT_TSS_POS)[0])
    assert regions["signature"].iloc[i] == regions["signature"].iloc[i - 1] == BIT_EXON_POS
    assert int(regions["end"].iloc[i - 1]) == ALT_TSS_POS


def test_scanner_partition_carries_the_alternative_tss(alt_tss_index):
    """⭐ THE W1b ASSERTION. The position handed to ``BamScanner.set_regions`` must include the
    alternative TSS. Against the pre-W1b tree this FAILS: the merge deleted it."""
    positions, _offsets, _types = build_region_partition_arrays(alt_tss_index)
    assert ALT_TSS_POS in set(positions.tolist())


def test_scanner_partition_is_exactly_the_region_region_bound_array(alt_tss_index):
    """The emitted triple IS the v8 graph, per reference, in ``ref_names`` order — not a
    re-derivation of it that could drift."""
    positions, offsets, types = build_region_partition_arrays(alt_tss_index)
    regions = alt_tss_index.regions_df

    n_refs = len(alt_tss_index.ref_names)
    assert offsets.shape == (n_refs + 1,)
    assert int(offsets[-1]) == positions.shape[0]
    # k regions per reference ⇒ k+1 region_bound positions; region_types is 1:1 with regions.
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
    # k+1 region_bound positions per reference against k regions: the payload's own boundary/region relation.
    np.testing.assert_array_equal(
        np.diff(offsets.astype(np.int64)) - 1, np.diff(ra.ref_offsets.astype(np.int64))
    )


def test_adjacent_regions_may_share_a_signature(alt_tss_index):
    """⭐ The invariant that REPLACES neighbour-differs.

    A signature-merged partition's defining property was that no two adjacent regions share a
    signature. v8 drops exactly that — 59.5 % of human transcript termini fall at a region_bound where the
    signature does not change, so the merge was deleting the region_bound that made them visible.
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
    """⭐ Two annotations differing by exactly one signature-invisible region_bound must not share a cache key.

    Under the retired merged partition they produced identical arrays and hashed the SAME
    (``645a0dd8aa560236``, measured) — so a payload scanned against one would have loaded silently
    against the other. That is the stale-cache failure ``partition_hash`` exists to make impossible,
    and it was live until the partition became the region bound array.
    """
    assert len(alt_tss_index.regions_df) == len(no_alt_tss_index.regions_df) + 1
    assert alt_tss_index.partition_hash != no_alt_tss_index.partition_hash
