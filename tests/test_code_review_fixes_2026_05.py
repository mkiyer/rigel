"""Regression tests for code-review findings (May 2026).

1. ``Locus.ref_id`` must hold the canonical resolver/BAM ref id
   (defined by ``index.ref_lengths`` insertion order via
   ``index.ref_name_to_id``), NOT the pandas categorical code of
   ``index.t_df["ref"]``.  The two id spaces only coincide when
   reference names are inserted in lexical order.

"""

import numpy as np
import pysam
import pytest

from rigel.index import TranscriptIndex
from rigel.locus import build_multi_loci
from rigel.scored_fragments import ScoredFragments


# ---------------------------------------------------------------------------
# Fixture: 2-ref index where pandas categorical ordering DIFFERS from
# the canonical (FASTA insertion) ordering.
#
# FASTA order: chrZ first (canonical id 0), chrA second (canonical id 1).
# Pandas ``astype("category")`` alphabetizes → category codes:
#   chrA → 0, chrZ → 1   (the OPPOSITE of the canonical mapping)
#
# The pre-fix code stored the categorical code in ``Locus.ref_id``,
# which silently routed downstream per-Locus gDNA estimation to the
# wrong contig.
# ---------------------------------------------------------------------------


_TWO_REF_GTF = """\
chrZ\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g_z"; transcript_id "t_z";
chrA\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g_a"; transcript_id "t_a";
"""


@pytest.fixture(scope="module")
def two_ref_index(tmp_path_factory):
    """Build a 2-ref index with chrZ before chrA in the FASTA."""
    base = tmp_path_factory.mktemp("two_ref")
    gtf_path = base / "two_ref.gtf"
    gtf_path.write_text(_TWO_REF_GTF)

    fasta_path = base / "genome.fa"
    seq = "N" * 2000
    with open(fasta_path, "w") as f:
        f.write(">chrZ\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
        f.write(">chrA\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(fasta_path))

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir, retain_test_structures=True)


def test_two_ref_index_has_id_space_mismatch(two_ref_index):
    """Sanity check: confirm the fixture exposes the ID-space mismatch
    the bug used to silently traverse.
    """
    idx = two_ref_index
    # Canonical ids: FASTA-insertion order.
    assert idx.ref_name_to_id["chrZ"] == 0
    assert idx.ref_name_to_id["chrA"] == 1
    # Pandas categorical codes: alphabetic.
    cats = list(idx.t_df["ref"].cat.categories)
    assert cats == ["chrA", "chrZ"], (
        f"Expected pandas to alphabetize categories; got {cats}."
    )


def test_locus_ref_id_is_canonical_not_categorical(two_ref_index):
    """``Locus.ref_id`` must equal ``index.ref_name_to_id[Locus.ref]``
    for every Locus produced by ``build_multi_loci``.

    Pre-fix, ``Locus.ref_id`` carried the pandas categorical code,
    which differed from the canonical id whenever the FASTA was not
    in lexical reference-name order (as in this fixture).
    """
    idx = two_ref_index

    # Two unambiguous fragments — one per ref — to force two MultiLoci.
    em_data = ScoredFragments(
        offsets=np.array([0, 1, 2], dtype=np.int64),
        t_indices=np.array([0, 1], dtype=np.int32),
        log_liks=np.zeros(2, dtype=np.float64),
        count_cols=np.zeros(2, dtype=np.uint8),
        coverage_weights=np.ones(2, dtype=np.float64),
        locus_t_indices=np.array([0, 1], dtype=np.int32),
        locus_count_cols=np.array([0, 0], dtype=np.uint8),
        is_spliced=np.zeros(2, dtype=bool),
        gdna_log_liks=np.zeros(2, dtype=np.float64),
        frag_ids=np.array([0, 1], dtype=np.int64),
        frag_class=np.zeros(2, dtype=np.int8),
        splice_type=np.zeros(2, dtype=np.uint8),
        n_units=2,
        n_candidates=2,
    )

    multi_loci = build_multi_loci(em_data, idx)
    assert len(multi_loci) >= 2, "Expected at least one MultiLocus per ref."

    saw_chrA = saw_chrZ = False
    for ml in multi_loci:
        for loc in ml.loci:
            assert loc.ref_id == idx.ref_name_to_id[loc.ref], (
                f"Locus(ref={loc.ref!r}) has ref_id={loc.ref_id}, "
                f"but canonical id is {idx.ref_name_to_id[loc.ref]}."
            )
            if loc.ref == "chrA":
                saw_chrA = True
            elif loc.ref == "chrZ":
                saw_chrZ = True
    assert saw_chrA and saw_chrZ, "Expected loci on both chrA and chrZ."
