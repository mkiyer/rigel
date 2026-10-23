"""`rigel.locus` — the `Locus` / `MultiLocus` surface, and the reference id a `Locus` carries. The
first block pins the dataclass invariants: `Locus` is a frozen contiguous interval with a `span`, a
`MultiLocus` carries one or more of them (a paralog cluster spans two references), and the
`span_bp` `build_multi_loci` computes is the sum of its loci's spans. The second pins that
`Locus.ref_id` is the canonical resolver/BAM ref id
given by `index.ref_lengths` insertion order through `index.ref_name_to_id`, and not the pandas
categorical code of `index.t_df["ref"]` — the two id spaces coincide only when reference names are
inserted in lexical order, so the fixture builds an index where they disagree and the wrong one would
route per-locus gDNA estimation to the wrong contig.
"""

import dataclasses

import numpy as np
import pytest
from _index_builder import write_fasta

from rigel.index import TranscriptIndex
from rigel.locus import Locus, MultiLocus, build_multi_loci
from rigel.scored_fragments import ScoredFragments


# ---------------------------------------------------------------------------
# Locus (contiguous interval) — frozen dataclass invariants
# ---------------------------------------------------------------------------


def test_locus_is_frozen():
    loc = Locus(ref="chr1", ref_id=0, start=10, end=20)
    with pytest.raises(dataclasses.FrozenInstanceError):
        loc.start = 5  # type: ignore[misc]


def test_locus_span_property():
    loc = Locus(ref="chr2", ref_id=1, start=100, end=350)
    assert loc.span == 250


# ---------------------------------------------------------------------------
# MultiLocus invariants
# ---------------------------------------------------------------------------


def _make_multi_locus(loci_seq, *, multi_locus_id=0):
    return MultiLocus(
        multi_locus_id=multi_locus_id,
        transcript_indices=np.array([0], dtype=np.int32),
        unit_indices=np.array([0], dtype=np.int32),
        span_bp=sum(loc.span for loc in loci_seq),
        loci=tuple(loci_seq),
    )


def test_multi_locus_single_interval():
    ml = _make_multi_locus([Locus("chr1", 0, 0, 10000)])
    assert len(ml.loci) == 1
    assert ml.span_bp == ml.loci[0].span == 10000


def test_multi_locus_multi_ref():
    """Paralog cluster spanning two refs produces multiple Locus intervals."""
    loci = (
        Locus("chr1", 0, 100, 600),
        Locus("chr2", 1, 50, 250),
    )
    ml = _make_multi_locus(loci)
    assert len(ml.loci) == 2
    assert ml.span_bp == sum(loc.span for loc in ml.loci) == 700
    # ref_id ordering assertion: the constructor accepts them sorted by
    # (ref_id, start); confirm we can iterate in declared order.
    assert [loc.ref for loc in ml.loci] == ["chr1", "chr2"]


# ---------------------------------------------------------------------------
# build_multi_loci integration
# ---------------------------------------------------------------------------


def test_build_multi_loci_span_bp_matches_sum_of_loci(mini_index):
    """``span_bp`` should equal Σ(end-start) over ``loci``."""
    from rigel.scored_fragments import ScoredFragments

    # Synthetic em_data: one ambiguous unit linking t0 and t1 (both on chr1).
    em_data = ScoredFragments(
        offsets=np.array([0, 2], dtype=np.int64),
        t_indices=np.array([0, 1], dtype=np.int32),
        log_liks=np.zeros(2, dtype=np.float32),
        count_cols=np.zeros(2, dtype=np.uint8),
        coverage_weights=np.ones(2, dtype=np.float32),
        is_spliced=np.zeros(1, dtype=bool),
        gdna_log_liks=np.zeros(1, dtype=np.float32),
        frag_ids=np.array([0], dtype=np.int64),
        frag_class=np.zeros(1, dtype=np.int8),
        splice_type=np.zeros(1, dtype=np.uint8),
        n_units=1,
        n_candidates=2,
    )
    multi_loci = build_multi_loci(em_data, mini_index)
    assert len(multi_loci) > 0
    for ml in multi_loci:
        loci_sum = sum(loc.span for loc in ml.loci)
        assert ml.span_bp == loci_sum
        for loc in ml.loci:
            assert isinstance(loc, Locus)
            assert loc.start <= loc.end


# ── ``Locus.ref_id`` is the canonical FASTA/BAM ref id, not a pandas categorical code ────────────


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

    fasta_path = write_fasta(base / "genome.fa", {"chrZ": 2000, "chrA": 2000})

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir, retain_test_structures=True)


def test_two_ref_index_has_id_space_mismatch(two_ref_index):
    """Non-vacuity: the fixture must actually expose the id-space mismatch, or the gate below holds
    over two id spaces that agree and proves nothing.
    """
    idx = two_ref_index
    # Canonical ids: FASTA-insertion order.
    assert idx.ref_name_to_id["chrZ"] == 0
    assert idx.ref_name_to_id["chrA"] == 1
    # Pandas categorical codes: alphabetic.
    cats = list(idx.t_df["ref"].cat.categories)
    assert cats == ["chrA", "chrZ"], f"Expected pandas to alphabetize categories; got {cats}."


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
        log_liks=np.zeros(2, dtype=np.float32),
        count_cols=np.zeros(2, dtype=np.uint8),
        coverage_weights=np.ones(2, dtype=np.float32),
        is_spliced=np.zeros(2, dtype=bool),
        gdna_log_liks=np.zeros(2, dtype=np.float32),
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


def test_the_components_of_nothing_have_the_shape_of_any_components():
    """The native union-find returned (labels, 0) when there was no unit or no transcript, and five CSR arrays
    otherwise; its one caller unpacks five, so an empty input would fail as an unpacking error far from the cause."""
    from rigel.native import connected_components

    out = connected_components(
        np.array([0], dtype=np.int64), np.array([], dtype=np.int32), np.int32(3)
    )
    assert len(out) == 5
    n_comp, t_off, t_flat, u_off, u_flat = out
    assert (
        n_comp == 0
        and list(t_off) == [0]
        and len(t_flat) == 0
        and list(u_off) == [0]
        and len(u_flat) == 0
    )
