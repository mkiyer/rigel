"""A fragment observed at several places: every place it could have come from enters the EM, gDNA sums over its
places, each place weighs its candidates by its capture, and an intergenic place's gDNA comes at a fixed rate no
locus carries (`ISSUES: multimapper-intergenic-alignments`, `ISSUES: the-gdna-term-is-a-mean-over-places`,
`ISSUES: mismatch-gate-nan-in-multimapper-gdna`).

The model. A fragment is observed as its usable hits ``S``; the molecule came from one component at one place, and
the aligner reported the set, so ``P(fragment | c) = Σ_{p ∈ S(c)} w(p) · f_c(ℓ_p) · s_c(p) · e(p) / Z_c``, ``S(c)`` the
places ``c`` can be at and ``Z_c`` the EM's capture-contracted length of ``c``. Three consequences:

* Every component's term is the SUM over its places, never the mean: a place is the hits whose extents overlap on one
  reference, alternatives read at the most likely of them, and distinct places add, for a transcript and for gDNA
  alike. A term at ``−∞`` adds nothing, so the hits' order cannot matter, in a library with no mismatch too.
* ``w(x) = Σ_o m_o(x) · c_o``, the deposit rule's mass of the placement at each object it deposits on times the
  object's capture efficiency, is the weight whose sum over a component's placements, on the law the efficiencies are
  read on, is that component's normaliser. Every candidate at a place carries it; on a unique fragment it is common
  to every candidate and is not applied. An intergenic place weighs 1: intergenic sequence is uncaptured (owner,
  2026-10-27).
* An intergenic place is in no locus. Its gDNA is a per-unit background ``λ · b_u`` beside the locus's components,
  whose absolute rates are ``M θ_c / L_c``, with ``M`` the locus's own data count re-estimated at every M-step: the
  maximiser of ``Σ_u log(Σ_c M θ_c q_uc / L_c + λ b_u) + Σ_d log(M θ_d) − M`` under the grouped prior on ``θ``. Its
  mass leaves the locus and joins the intergenic gDNA. It competes at its rate, so no difference of scores prunes it;
  and a multimapper whose only places the gates leave are intergenic is intergenic gDNA outright.

The owner's five cases (2026-10-19): intergenic only, unique, unspliced — certified gDNA; intergenic only, unique,
spliced — withheld; a multimapper all of whose alignments are intergenic and unspliced — certified gDNA; one with
intergenic spliced and unspliced alignments — the spliced ones rejected, gDNA by the unspliced; one with genic and
intergenic alignments, spliced or not — the transcript alignments and the intergenic unspliced ones compete in the EM.
"""

from __future__ import annotations

import dataclasses
import itertools
import math
from dataclasses import dataclass

import numpy as np
import pysam
import pytest
from _index_builder import write_fasta

from rigel.buffer import FRAG_MULTIMAPPER, FRAG_UNAMBIG
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig
from rigel.estimator import AbundanceEstimator
from rigel.frag_length_model import FragmentLengthModel
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.routing import FragmentRouter
from rigel.scoring import build_scorer, sequence_log_penalties
from rigel.splice import SpliceStrandCol, SpliceType
from rigel.stats import PipelineStats
from rigel.strand_model import StrandModel, StrandModels
from rigel.types import Strand

# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  1. THE SCAN — the owner's five cases on hand-built multimappers
# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#
# The toy reference: one gene, tA (+), exons [1000, 1300) and [2000, 2300) with the intron [1300, 2000) between; every
# other base intergenic. Its synthetic span covers [1000, 2300).

GENOME = 20_000
R1_FWD, R2_REV = 0x1 | 0x2 | 0x40 | 0x20, 0x1 | 0x2 | 0x80 | 0x10
SECONDARY = 0x100

GENIC = (1050, "100M", 1150, "100M")  # [1050, 1250): inside tA's first exon
GENIC_SPLICED = (1250, "50M700N50M", 2100, "100M")  # across tA's annotated intron
INTERGENIC = (8050, "100M", 8150, "100M")  # [8050, 8250)
INTERGENIC_2 = (12050, "100M", 12150, "100M")  # [12050, 12250)
INTERGENIC_SPLICED = (14050, "50M500N50M", 14700, "100M")  # a junction no transcript holds


def _record(qname, flag, ref, pos, cigar, mate_pos, nh, hi, xs=None, nm=0):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.flag = flag
    a.reference_id = ref
    a.reference_start = pos
    a.mapping_quality = 60 if nh == 1 else 3
    a.cigarstring = cigar
    n = sum(int(x) for x in __import__("re").findall(r"(\d+)[MIS=X]", cigar))
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = ref
    a.next_reference_start = mate_pos
    tags = [("NH", nh, "i"), ("HI", hi, "i"), ("NM", nm, "i")]
    if xs is not None:
        tags.append(("XS", xs, "A"))
    a.set_tags(tags)
    return a


def _template(qname, *places, nm=None):
    """One read pair's records at each of its places, the first primary and the rest secondary. A place is
    ``(f_pos, f_cigar, r_pos, r_cigar)`` on the first reference, or ``(ref, f_pos, f_cigar, r_pos, r_cigar)``;
    ``nm`` each place's mismatches (none by default)."""
    nh, records = len(places), []
    nm = nm or (0,) * nh
    for hi, place in enumerate(places, start=1):
        ref, (f_pos, f_cigar, r_pos, r_cigar) = (
            (place[0], place[1:]) if len(place) == 5 else (0, place)
        )
        extra = 0 if hi == 1 else SECONDARY
        xs = "+" if "N" in f_cigar or "N" in r_cigar else None
        records += [
            _record(qname, R1_FWD | extra, ref, f_pos, f_cigar, r_pos, nh, hi, xs, nm[hi - 1]),
            _record(qname, R2_REV | extra, ref, r_pos, r_cigar, f_pos, nh, hi, xs, nm[hi - 1]),
        ]
    return records


def _build_index(base, refs):
    """``refs``: every reference's name, each holding tA's exon layout as its own gene."""
    fasta, gtf = base / "g.fa", base / "a.gtf"
    write_fasta(fasta, {ref: GENOME for ref in refs})
    gtf.write_text(
        "".join(
            f'{ref}\tt\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g{ref}"; transcript_id "t{ref}";\n'
            for ref in refs
            for s, e in [(1000, 1300), (2000, 2300)]
        )
    )
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    return TranscriptIndex.load(str(base / "idx"))


@pytest.fixture(scope="module")
def toy_index(tmp_path_factory):
    return _build_index(tmp_path_factory.mktemp("places"), ["chr1"])


@pytest.fixture(scope="module")
def two_ref_index(tmp_path_factory):
    return _build_index(tmp_path_factory.mktemp("places2"), ["chr1", "chr2"])


def _write_bam(tmp_path, index, records):
    """``records`` as a name-sorted BAM over the index's references."""
    bam = str(tmp_path / "t.bam")
    sq = [{"SN": ref, "LN": GENOME} for ref in index.ref_names]
    with pysam.AlignmentFile(
        bam, "wb", header={"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": sq}
    ) as out:
        for r in records:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)
    return bam


def _scan_chunks(tmp_path, index, records, **scan):
    """The scan's stats and its buffered chunks, consumed."""
    bam = _write_bam(tmp_path, index, records)
    stats, _strand, buffer, _payload = scan_and_buffer(
        bam, index, BamScanConfig(total_threads=1, **scan)
    )
    return stats, list(buffer.iter_chunks_consuming())


def _scan(tmp_path, index, records):
    """The scan's stats, and per buffered row of the template its count of transcript candidates."""
    stats, chunks = _scan_chunks(tmp_path, index, records)
    return stats, [int(n) for chunk in chunks for n in np.diff(chunk.t_offsets)]


def test_case_1_a_unique_intergenic_unspliced_fragment_is_certified_gdna(toy_index, tmp_path):
    stats, rows = _scan(tmp_path, toy_index, _template("c1", INTERGENIC))
    assert stats.n_fragments == 1 and stats.n_intergenic == 1
    assert rows == [0]  # buffered for calibration's geometry; the scorer gives it no unit


def test_case_2_a_unique_intergenic_spliced_fragment_is_withheld(toy_index, tmp_path):
    stats, rows = _scan(tmp_path, toy_index, _template("c2", INTERGENIC_SPLICED))
    assert stats.n_fragments == 0 and stats.n_withheld_unannotated_splice == 1
    assert rows == []


def test_case_3_an_all_intergenic_multimapper_is_certified_gdna_and_reaches_no_em_unit(
    toy_index, tmp_path
):
    stats, rows = _scan(tmp_path, toy_index, _template("c3", INTERGENIC, INTERGENIC_2))
    assert stats.n_fragments == 1 and stats.n_intergenic == 1 and stats.n_placed_multimapping == 1
    assert rows == []


def test_case_4_a_multimappers_intergenic_splice_is_rejected_and_its_unspliced_place_is_gdna(
    toy_index, tmp_path
):
    stats, rows = _scan(tmp_path, toy_index, _template("c4", INTERGENIC_SPLICED, INTERGENIC))
    assert stats.n_fragments == 1 and stats.n_intergenic == 1
    assert rows == [0]  # one usable hit left: a unique intergenic fragment, as case 1


def test_case_5_a_genic_and_intergenic_multimapper_brings_both_places_to_the_em(
    toy_index, tmp_path
):
    """The intergenic hit is buffered beside the genic one, the template counted once and not intergenic: the
    EM sees the place gDNA most likely came from. Today the scan drops it, and the fragment reaches the EM at
    its genic place alone."""
    stats, rows = _scan(tmp_path, toy_index, _template("c5", GENIC, INTERGENIC))
    assert stats.n_fragments == 1 and stats.n_intergenic == 0 and stats.n_placed_multimapping == 1
    assert sorted(rows) == [
        0,
        2,
    ]  # the intergenic place, and the genic one with tA and its synthetic span


def test_case_5_with_a_spliced_genic_place_brings_its_intergenic_place_to_the_em(
    toy_index, tmp_path
):
    stats, rows = _scan(tmp_path, toy_index, _template("c5s", GENIC_SPLICED, INTERGENIC))
    assert stats.n_fragments == 1 and stats.n_intergenic == 0
    assert sorted(rows) == [0, 1]  # the intergenic place, and the spliced one with tA alone


def _place_masses(chunks):
    """``{row: [(kind, flat object id, mass)]}`` over every buffered row with a transcript candidate."""
    out, row = {}, 0
    for chunk in chunks:
        off = chunk.place_offsets
        for i, n in enumerate(np.diff(chunk.t_offsets)):
            if n > 0:
                sl = slice(int(off[i]), int(off[i + 1]))
                out[row] = list(
                    zip(
                        chunk.place_kind[sl].tolist(),
                        chunk.place_index[sl].tolist(),
                        chunk.place_mass[sl].tolist(),
                    )
                )
            row += 1
    return out


def test_a_genic_place_on_a_later_reference_names_its_objects_on_the_flat_axes(
    two_ref_index, tmp_path
):
    """The masses are written at the objects' FLAT ids, the axes the efficiencies are indexed on: a place inside the
    second reference's first exon names that exon's region, past every region of the first reference."""
    import rigel.native as N
    from rigel.calibration.region_arrays import RegionArrays

    stats, chunks = _scan_chunks(tmp_path, two_ref_index, _template("c5r", (1, *GENIC), INTERGENIC))
    assert stats.n_fragments == 1 and stats.n_intergenic == 0
    ra = RegionArrays.from_index(two_ref_index)
    chr2 = two_ref_index.ref_name_to_id["chr2"]
    start, ref_id = np.asarray(ra.start), np.asarray(ra.ref_id)
    exon = int(np.flatnonzero((ref_id == chr2) & (start == 1000))[0])
    assert exon >= int(np.count_nonzero(ref_id != chr2))
    assert list(_place_masses(chunks).values()) == [[(int(N.PLACE_REGION), exon, 1.0)]]


def test_the_place_masses_survive_a_spill_to_disk(toy_index, tmp_path):
    """A buffer forced to disk hands back every row's masses as the in-memory buffer holds them."""
    records = _template("c5", GENIC_SPLICED, INTERGENIC) + _template("c5b", GENIC, INTERGENIC_2)
    _stats, in_memory = _scan_chunks(tmp_path, toy_index, records)
    _stats, spilled = _scan_chunks(
        tmp_path, toy_index, records, buffer_size_bytes=1, spill_dir=tmp_path
    )
    masses = _place_masses(in_memory)
    assert masses and all(masses.values())
    assert _place_masses(spilled) == masses


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  2. THE SCORE — gDNA sums over its places; the intergenic term; each place's weight
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class Hit:
    """One buffered row: its transcript candidates (none for an intergenic place) and what the scorer reads."""

    t: tuple[int, ...]
    splice_type: SpliceType = SpliceType.UNSPLICED
    footprint: int = 250
    nm: int = 0
    clip: int = 0
    log_weight: float = 0.0
    start: int | None = (
        None  # its genomic start; none given: a place of its own, apart from the group's others
    )
    exon_bp: int = 100  # each candidate's exonic bases, of 100: fewer is overhang


class _Chunk:
    """The ``size`` and ``to_scoring_arrays`` the router reads, from rows ``(frag_id, class, Hit)``, each row's place
    weight among them."""

    def __init__(self, rows):
        hits = [h for _f, _c, h in rows]
        self.size = len(hits)
        offsets = np.concatenate([[0], np.cumsum([len(h.t) for h in hits])])
        n_cand = int(offsets[-1])
        frag_ids = [f for f, _c, _h in rows]
        fclass = [c for _f, c, _h in rows]
        self._arrays = {
            "t_offsets": np.asarray(offsets, dtype=np.int32),
            "t_indices": np.asarray([t for h in hits for t in h.t], dtype=np.int32),
            "frag_lengths": np.full(n_cand, 200, dtype=np.int32),
            "exon_bp": np.asarray([h.exon_bp for h in hits for _ in h.t], dtype=np.uint16),
            "splice_type": np.asarray([int(h.splice_type) for h in hits], dtype=np.uint8),
            "align_strand": np.full(len(hits), int(Strand.POS), dtype=np.uint8),
            "fragment_classes": np.asarray(fclass, dtype=np.uint8),
            "frag_id": np.asarray(frag_ids, dtype=np.int64),
            "read_length": np.full(len(hits), 100, dtype=np.uint16),
            "genomic_footprint": np.asarray([h.footprint for h in hits], dtype=np.int32),
            "genomic_start": np.asarray([h.start for h in hits], dtype=np.int32),
            "nm": np.asarray([h.nm for h in hits], dtype=np.uint16),
            "clip": np.asarray([h.clip for h in hits], dtype=np.uint16),
            "log_weight": np.asarray([h.log_weight for h in hits], dtype=np.float32),
        }

    def to_scoring_arrays(self):
        return self._arrays


class _Buffer:
    def __init__(self, chunks):
        self._chunks = list(chunks)
        self.total_fragments = sum(c.size for c in chunks)

    def iter_chunks_consuming(self):
        while self._chunks:
            yield self._chunks.pop(0)


class _Stamps(dict):
    """The annotation table the router writes, read back as ``{frag_id: (ZF, best_tid)}``."""

    def add(self, *, frag_id, best_tid, tx_flags, **_):
        self[int(frag_id)] = (int(tx_flags), int(best_tid))

    def add_batch(self, *, frag_ids, best_tids, tx_flags, **_):
        for f, t, z in zip(frag_ids, best_tids, tx_flags):
            self.add(frag_id=f, best_tid=t, tx_flags=z)


class _Index:
    def __init__(self, n_t):
        import pandas as pd

        self.t_to_g_arr = np.zeros(n_t, dtype=np.int64)
        self.t_to_strand_arr = np.full(n_t, int(Strand.POS), dtype=np.int8)
        self.num_transcripts = n_t
        self.num_genes = 1
        self.t_df = pd.DataFrame(
            {
                "t_id": [f"t{i}" for i in range(n_t)],
                "ref": ["chr1"] * (n_t // 2)
                + ["chr2"] * (n_t - n_t // 2),  # transcripts 0, 1 and 2, 3
                "start": np.zeros(n_t, dtype=np.int32),
                "end": np.full(n_t, 10000, dtype=np.int32),
                "length": np.full(n_t, 1000, dtype=np.int32),
                "is_single_exon": np.zeros(n_t, dtype=bool),
                "is_synthetic": np.zeros(n_t, dtype=bool),
            }
        )

    def build_exon_csr(self):
        empty = np.empty(0, dtype=np.int32)
        return np.zeros(self.num_transcripts + 1, dtype=np.int32), empty, empty, empty


N_T = 4

#: A multimapper's hit's sequence prices at a typical library's mismatch rate, per mismatch and per clipped base.
MISMATCH_LP, CLIP_LP = sequence_log_penalties(0.002)


def _route_all(groups, classes=None, *, fl=None, mismatch_log_penalty=None, split=None):
    """Score fragments, each a group of hits, through the shipped scorer and router; ``split`` cuts the rows into two
    chunks before that row. Returns the scored units, the router's stats and the rows it stamped."""
    classes = classes or [FRAG_MULTIMAPPER] * len(groups)
    index = _Index(N_T)
    fl = fl if fl is not None else _live_fl()
    scorer = build_scorer(
        StrandModels(exonic_spliced=StrandModel(pos_pos=20)),
        fl,
        fl,
        index,
        overhang_log_penalty=FragmentScoringConfig().overhang_log_penalty,
        mismatch_log_penalty=MISMATCH_LP if mismatch_log_penalty is None else mismatch_log_penalty,
        clip_log_penalty=CLIP_LP,
    )
    estimator = AbundanceEstimator(N_T, em_config=EMConfig(seed=1))
    stats, stamps = PipelineStats(), _Stamps()
    router = FragmentRouter(scorer, estimator, stats, index, annotations=stamps)
    rows = [
        (f, c, h if h.start is not None else dataclasses.replace(h, start=1000 + 10_000 * k))
        for f, (g, c) in enumerate(zip(groups, classes))
        for k, h in enumerate(g)
    ]
    chunks = [rows] if split is None else [rows[:split], rows[split:]]
    em = router.route(_Buffer([_Chunk(c) for c in chunks]), 1_000_000)
    return em, stats, stamps


def _route(groups, classes=None, **kw):
    """The scored units alone (:func:`_route_all`)."""
    return _route_all(groups, classes, **kw)[0]


def _live_fl():
    """A length law that varies with length on [50, 500] and has no mass beyond."""
    pmf = np.zeros(1001)
    pmf[50:501] = np.linspace(1.0, 2.0, 451)
    return FragmentLengthModel.from_pmf(pmf, 1000)


def _term(footprint, nm=0, *, fl=None, mismatch_log_penalty=None):
    """One place's gDNA term: the law at the footprint and ½, read off a unique unspliced fragment (one hit, whose
    sequence is common to its candidates and not read), and the hit's mismatches at a multimapper's price."""
    em = _route([[Hit((0,), footprint=footprint)]], [FRAG_UNAMBIG], fl=fl)
    price = MISMATCH_LP if mismatch_log_penalty is None else mismatch_log_penalty
    return float(em.gdna_log_liks[0]) + (nm * price if nm else 0.0)


def _intergenic(em):
    """The per-unit intergenic term the scorer emits beside the locus's gDNA term."""
    out = getattr(em, "intergenic_log_liks", None)
    assert out is not None, "the scorer emits no intergenic term"
    return out


def test_gdna_sums_over_two_places_never_their_mean():
    em = _route([[Hit((0,), footprint=250), Hit((1,), footprint=300, nm=1)]])
    expected = np.logaddexp(_term(250), _term(300, nm=1))
    assert float(em.gdna_log_liks[0]) == pytest.approx(expected, abs=1e-5)


def test_gdna_sums_over_three_places():
    hits = [
        Hit((0,), footprint=250),
        Hit((1,), footprint=300, nm=1),
        Hit((2,), footprint=350, nm=2),
    ]
    em = _route([hits])
    expected = np.logaddexp.reduce([_term(250), _term(300, nm=1), _term(350, nm=2)])
    assert float(em.gdna_log_liks[0]) == pytest.approx(expected, abs=1e-5)


def test_two_identical_places_give_gdna_twice_one_places_term():
    """Two places the molecule could have come from, each as likely: gDNA's term doubles, while each RNA candidate
    keeps its own place's term. The mean would hold gDNA at one place's term and bias the fragment toward RNA."""
    em = _route([[Hit((0,), footprint=250), Hit((1,), footprint=250)]])
    assert float(em.gdna_log_liks[0]) == pytest.approx(_term(250) + math.log(2.0), abs=1e-5)


@pytest.mark.parametrize("order", list(itertools.permutations(range(3))))
def test_in_a_library_with_no_mismatch_the_gdna_term_is_the_same_in_every_hit_order(order):
    """In a library with no mismatch a mismatched place's term is ``−∞`` and adds nothing to the sum, so the
    fragment's gDNA term is its one perfect place's, whichever order the hits arrive in. Two ``−∞`` terms ahead of a
    finite one once made the running log-sum-exp NaN and gDNA was dropped, and in the other orders the mean divided by
    three."""
    hits = [Hit((0,), nm=2), Hit((1,), nm=3), Hit((2,), nm=0)]
    gate = -math.inf
    em = _route([[hits[i] for i in order]], mismatch_log_penalty=gate)
    assert not em.is_spliced[0]
    assert float(em.gdna_log_liks[0]) == pytest.approx(
        _term(250, nm=0, mismatch_log_penalty=gate), abs=1e-5
    )


def test_a_case5_units_intergenic_place_is_its_intergenic_term_and_not_the_loci():
    """The genic place's gDNA is the locus component's; the intergenic place's is the intergenic term's."""
    em = _route([[Hit((0,), footprint=250), Hit((), footprint=300)]])
    assert float(em.gdna_log_liks[0]) == pytest.approx(_term(250), abs=1e-5)
    assert float(_intergenic(em)[0]) == pytest.approx(_term(300), abs=1e-5)


def test_a_case5_unit_with_a_spliced_genic_place_competes_rather_than_being_certified():
    """RNA's only reading is the spliced place, gDNA's the intergenic one: no locus gDNA candidate, an intergenic
    term, and the unit is not certified RNA (owner, 2026-10-19: case 5, spliced or not, competes in the EM)."""
    em = _route([[Hit((0,), SpliceType.SPLICED_ANNOT, footprint=400), Hit((), footprint=300)]])
    assert em.gdna_log_liks[0] == -np.inf
    assert float(_intergenic(em)[0]) == pytest.approx(_term(300), abs=1e-5)
    assert not em.is_spliced[0]


def test_a_unit_with_no_intergenic_place_has_no_intergenic_term():
    em = _route([[Hit((0,), footprint=250), Hit((1,), footprint=300)]])
    assert _intergenic(em)[0] == -np.inf


def _geometric_fl():
    """A law whose mass falls geometrically with length, so a footprint moves a gDNA term by whole nats."""
    pmf = np.zeros(1001)
    lengths = np.arange(50, 1001)
    pmf[50:] = np.exp(-lengths / 20.0)
    return FragmentLengthModel.from_pmf(pmf, 1000)


def test_the_intergenic_term_is_kept_wherever_its_place_survives():
    """It competes at the intergenic rate against the locus's rates and lengths, none of which are in the scores, so
    no difference of scores bounds its posterior (the EM gate scaling a weight and a length together shows it): it is
    kept at every footprint, those the locus's gDNA bound would prune included, and the unit is never certified."""
    fl = _geometric_fl()
    delta = -math.log(FragmentScoringConfig().pruning_min_posterior)
    beyond = 0
    for footprint in range(200, 700, 7):
        em = _route(
            [[Hit((0,), SpliceType.SPLICED_ANNOT, footprint=200), Hit((), footprint=footprint)]],
            fl=fl,
        )
        term = fl._log_prob[footprint] + math.log(
            0.5
        )  # the place's term itself: a unique one's is pruned beyond
        assert float(_intergenic(em)[0]) == pytest.approx(term, abs=1e-5), footprint
        assert not em.is_spliced[0]
        beyond += float(np.max(em.log_liks)) - term > delta
    assert beyond


#: A genic place that fails a gate beside an intergenic one that does not, with the gate's penalty.
GATED_GENIC_PLACE = {
    "zero capture weight": ([Hit((0,), log_weight=-math.inf), Hit(())], None),
    "a mismatch in a library with none": ([Hit((0,), nm=1), Hit((), nm=0)], -math.inf),
}


@pytest.mark.parametrize("case", list(GATED_GENIC_PLACE))
@pytest.mark.parametrize("reverse", [False, True], ids=["genic first", "intergenic first"])
@pytest.mark.parametrize("split", [None, 1], ids=["one chunk", "across chunks"])
def test_a_multimapper_whose_only_surviving_place_is_intergenic_is_intergenic_gdna(
    case, reverse, split
):
    """Every RNA candidate fails and the intergenic place does not: the fragment is gDNA there, in no locus, counted and
    stamped as a multimapper none of whose hits resolves is and never gated, whatever the hits' order and wherever a
    chunk boundary falls. It reaches the summary's intergenic gDNA."""
    from rigel.annotate import AF_GDNA_INTERGENIC
    from rigel.cli import _quantification_pools

    hits, penalty = GATED_GENIC_PLACE[case]
    em, stats, stamps = _route_all(
        [hits[::-1] if reverse else hits], mismatch_log_penalty=penalty, split=split
    )
    assert em.n_units == 0
    assert stats.n_intergenic == 1 and stats.n_gated_out == 0
    assert stamps == {0: (AF_GDNA_INTERGENIC, -1)}
    pools = _quantification_pools(stats, count_annotated=0.0, count_nascent=0.0, count_gdna_em=0.0)
    assert pools["count_gdna"] == 1.0 and pools["fraction_intergenic"] == 1.0


def test_a_multimapper_no_place_of_which_survives_is_gated():
    from rigel.annotate import AF_UNRESOLVED

    em, stats, stamps = _route_all(
        [[Hit((0,), nm=1), Hit((), nm=1)]], mismatch_log_penalty=-math.inf
    )
    assert em.n_units == 0
    assert stats.n_gated_out == 1 and stats.n_intergenic == 0
    assert stamps == {0: (AF_UNRESOLVED, -1)}


def test_every_candidate_at_a_place_carries_that_places_weight():
    """RNA's candidates and gDNA's term at a place carry its log-weight, each place its own; the intergenic place
    carries the weight it is given (1 for every intergenic place, the weight function's rule)."""
    unweighted = [Hit((0,), footprint=250), Hit((1,), footprint=250), Hit((), footprint=300)]
    weighted = [
        Hit((0,), footprint=250, log_weight=math.log(3.0)),
        Hit((1,), footprint=250, log_weight=math.log(0.5)),
        Hit((), footprint=300),
    ]
    base, em = _route([unweighted]), _route([weighted])
    by_t = lambda e: dict(zip(e.t_indices.tolist(), e.log_liks.tolist()))  # noqa: E731
    assert by_t(em)[0] - by_t(base)[0] == pytest.approx(math.log(3.0), abs=1e-5)
    assert by_t(em)[1] - by_t(base)[1] == pytest.approx(math.log(0.5), abs=1e-5)
    expected = np.logaddexp(_term(250) + math.log(3.0), _term(250) + math.log(0.5))
    assert float(em.gdna_log_liks[0]) == pytest.approx(expected, abs=1e-5)
    assert float(_intergenic(em)[0]) == pytest.approx(_term(300), abs=1e-5)


def _rna(em, t=0):
    """Transcript ``t``'s term in the unit."""
    return dict(zip(em.t_indices.tolist(), em.log_liks.tolist()))[t]


def test_a_transcript_sums_over_its_distinct_places():
    """Two hits holding one transcript at disjoint extents are two places it could have come from: its term is their
    sum, as gDNA's is, never one hit's."""
    one_place = _route([[Hit((0,), start=1000), Hit((1,), start=5000)]])
    two_places = _route([[Hit((0,), start=1000), Hit((0,), start=5000)]])
    assert _rna(two_places) == pytest.approx(_rna(one_place) + math.log(2.0), abs=1e-5)


def test_a_transcripts_overlapping_hits_are_one_place_read_at_the_more_likely():
    """Hits whose extents overlap are alternatives at one place: they count once, at the more likely, whose overhang
    and mismatches are already in its score. The least-overhang rule took the hit with no overhang and NM 3 over the
    one with overhang 1 and NM 0, 4.6 nats more likely (`ISSUES: merge-keeps-least-overhang`)."""
    em = _route([[Hit((0,), start=1000, nm=3), Hit((0,), start=1020, exon_bp=99)]])
    alone = _route([[Hit((0,), start=1020, exon_bp=99), Hit((1,), start=9000)]])
    assert _rna(em) == pytest.approx(_rna(alone), abs=1e-5)


def test_in_a_library_with_no_mismatch_a_transcript_keeps_its_surviving_alternative():
    """A library with no mismatch takes the hit with one, and the transcript is still read at its other alternative;
    the least-overhang rule kept the excluded one and lost the transcript."""
    em = _route(
        [[Hit((0,), start=1000, nm=1), Hit((0,), start=1020, exon_bp=99)]],
        mismatch_log_penalty=-math.inf,
    )
    assert 0 in em.t_indices.tolist() and math.isfinite(_rna(em))


def test_gdnas_overlapping_hits_are_one_place_and_the_same_coordinates_on_two_references_two():
    """gDNA follows the transcripts' place definition: two of its hits overlapping on one reference are alternatives at
    one place and count once, at the more likely; at the same coordinates on two references they are two places."""
    one_place = _route([[Hit((0,), start=1000), Hit((1,), start=1020, nm=1)]])
    two_places = _route([[Hit((0,), start=1000), Hit((2,), start=1000)]])
    assert float(one_place.gdna_log_liks[0]) == pytest.approx(_term(250), abs=1e-5)
    assert float(two_places.gdna_log_liks[0]) == pytest.approx(_term(250) + math.log(2.0), abs=1e-5)


def test_a_unique_fragment_is_untouched_by_a_place_weight():
    """One place: its weight is common to every candidate and cancels, so it is not applied and the fragment
    scores bit for bit as it does without one (a guard for the change; it holds today)."""
    for hit_kw, fclass in ((dict(t=(0, 1)), FRAG_UNAMBIG), (dict(t=(0,)), FRAG_UNAMBIG)):
        plain = _route([[Hit(**hit_kw)]], [fclass])
        weighted = _route([[Hit(**hit_kw, log_weight=math.log(5.0))]], [fclass])
        for name in ("log_liks", "gdna_log_liks", "is_spliced", "t_indices", "offsets"):
            assert np.array_equal(getattr(plain, name), getattr(weighted, name)), name


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  3. THE PLACE WEIGHT — the deposit rule's mass at each object, times its efficiency
# ═══════════════════════════════════════════════════════════════════════════════════════════════════

#: Two outer regions longer than every fragment, so no placement meets a wall the shares do not see, and a cluster
#: of short regions between them, so one placement crosses several boundaries. The junction [1000, 1100) jumps
#: two of them.
BOUNDS = [0, 1000, 1060, 1100, 1190, 2400]
SJ = [(0, 1000, 1100, Strand.POS)]
LENGTHS = (60, 150, 300)
MAX_LENGTH = 1000


def _accumulators():
    from native._accumulator_reference import Accumulator as Reference, Partition

    from rigel.native import Accumulator as Native

    partition = Partition.from_region_bounds(
        [BOUNDS], region_types=[[0] * (len(BOUNDS) - 1)], sj=SJ
    )
    native = Native(
        region_bounds=np.asarray(BOUNDS, dtype=np.int64),
        region_types=np.zeros(len(BOUNDS) - 1, dtype=np.uint8),
        max_length=MAX_LENGTH,
        ref=0,
    )
    native.set_sj(
        np.ascontiguousarray(partition.sj_offsets, dtype=np.int32),
        np.ascontiguousarray(partition.sj_boundary_right, dtype=np.int32),
        np.ascontiguousarray(partition.sj_strand, dtype=np.int8),
    )
    return partition, Reference, native


def _native_masses(native, start, end, introns=()):
    """``{(kind, object): mass}`` of one placement on its observed path, from the native dry run."""
    import rigel.native as N

    kinds, ids, masses = native.placement_masses(
        start=start, end=end, observed_introns=tuple(introns), sj_strand=int(Strand.POS)
    )
    names = {
        int(N.PLACE_REGION): "region",
        int(N.PLACE_BOUNDARY): "boundary",
        int(N.PLACE_SJ): "sj",
    }
    return {(names[int(k)], int(i)): float(m) for k, i, m in zip(kinds, ids, masses)}


def _reference_masses(partition, Reference, start, end, introns=()):
    """The same placement deposited into a fresh specification, read back off every mass bank it touched."""
    acc = Reference(partition, max_fragment_length=MAX_LENGTH)
    acc.deposit(0, start, end, observed_introns=tuple(introns), sj_strand=int(Strand.POS))
    t, out = acc.tally, {}
    for r in np.flatnonzero(t.region_contained_count.sum(axis=1)):
        out[("region", int(r))] = 1.0
    for bank in (t.boundary_unspliced_mass, t.boundary_spliced_mass):
        for e in np.flatnonzero(bank):
            out[("boundary", int(e))] = out.get(("boundary", int(e)), 0.0) + float(bank[e])
    for j in np.flatnonzero(t.sj_mass.sum(axis=1)):
        out[("sj", int(j))] = float(t.sj_mass[j].sum())
    return out


@pytest.mark.parametrize(
    "start,end,introns",
    [
        (100, 400, ()),  # contained in the first region
        (950, 1050, ()),  # across one boundary
        (900, 1300, ()),  # across four boundaries
        (950, 1150, [(1000, 1100)]),  # spliced, the junction its only object
        (900, 1300, [(1000, 1100)]),  # spliced, beside crossed boundaries
    ],
)
def test_a_placements_masses_are_what_the_deposit_rule_gives_each_object(start, end, introns):
    partition, Reference, native = _accumulators()
    got = _native_masses(native, start, end, introns)
    want = _reference_masses(partition, Reference, start, end, introns)
    assert got.keys() == want.keys()
    for key in want:
        assert got[key] == pytest.approx(want[key], abs=1e-12), key
    assert sum(got.values()) == pytest.approx(1.0, abs=1e-12)


def test_the_place_weight_sums_to_gdnas_normaliser_over_every_placement():
    """``Σ_x f(ℓ_x) w(x) = Σ_r S_r c_r + Σ_e M_e c_e`` over every placement on the reference: the weight and the
    capture-contracted length the EM divides gDNA by are one rule, at any efficiencies."""
    from rigel.calibration.capture_eff_length import place_log_weights
    from rigel.calibration.effective_length import conserved_cut_shares, contained_eff_length

    _partition, _reference, native = _accumulators()
    pmf = np.zeros(MAX_LENGTH + 1)
    pmf[list(LENGTHS)] = [0.5, 0.3, 0.2]
    rng = np.random.default_rng(7)
    n_regions, n_boundaries = len(BOUNDS) - 1, len(BOUNDS) - 2
    c_region = rng.uniform(0.2, 9.0, n_regions)
    c_boundary = rng.uniform(0.2, 9.0, n_boundaries)

    offsets, kinds, ids, masses, f = [0], [], [], [], []
    for length, p in zip(LENGTHS, pmf[list(LENGTHS)]):
        for start in range(0, BOUNDS[-1] - length + 1):
            k, i, m = native.placement_masses(
                start=start, end=start + length, sj_strand=int(Strand.POS)
            )
            kinds += list(k)
            ids += list(i)
            masses += list(m)
            offsets.append(len(kinds))
            f.append(p)
    log_w = place_log_weights(
        np.asarray(offsets, dtype=np.int64),
        np.asarray(kinds, dtype=np.int8),
        np.asarray(ids, dtype=np.int64),
        np.asarray(masses, dtype=np.float64),
        intergenic=np.zeros(len(f), dtype=bool),
        c_region=c_region,
        c_boundary=c_boundary,
        c_sj=np.ones(len(SJ)),
    )
    log_w = np.asarray(log_w, dtype=np.float64)
    lhs = float(np.sum(np.asarray(f) * np.exp(log_w)))

    bounds = np.asarray(BOUNDS, dtype=np.float64)
    region_len = np.diff(bounds)
    S = contained_eff_length(region_len, pmf)
    left, right = conserved_cut_shares(
        pmf, region_len[:-1], region_len[1:], bounds[1:-1], bounds[-1] - bounds[1:-1]
    )
    rhs = float(np.sum(S * c_region) + np.sum((left + right) * c_boundary))
    # The scorer reads ``log w`` in float32: rounding it moves each weight by at most half an ulp of ``log w``, so
    # the sum by at most that share of itself. Derived from the representation, not fitted.
    tolerance = 0.5 * float(np.finfo(np.float32).eps) * float(np.max(np.abs(log_w)))
    assert lhs == pytest.approx(rhs, rel=tolerance)


def test_a_place_weighs_its_objects_by_mass_and_an_intergenic_place_weighs_one():
    from rigel.calibration.capture_eff_length import place_log_weights
    import rigel.native as N

    R, B, J = int(N.PLACE_REGION), int(N.PLACE_BOUNDARY), int(N.PLACE_SJ)
    # rows: a unique fragment (no objects), one contained place, one crossing two boundaries, one spliced, and an
    # intergenic place whose region reads captured
    offsets = np.asarray([0, 0, 1, 3, 5, 6], dtype=np.int64)
    kinds = np.asarray([R, B, B, B, J, R], dtype=np.int8)
    ids = np.asarray([2, 0, 1, 1, 0, 3], dtype=np.int64)
    masses = np.asarray([1.0, 0.4, 0.6, 0.5, 0.5, 1.0])
    log_w = place_log_weights(
        offsets,
        kinds,
        ids,
        masses,
        intergenic=np.asarray([False, False, False, False, True]),
        c_region=np.asarray([1.0, 1.0, 3.0, 17.0]),
        c_boundary=np.asarray([2.0, 5.0]),
        c_sj=np.asarray([11.0]),
    )
    assert np.exp(np.asarray(log_w, dtype=np.float64)) == pytest.approx(
        [1.0, 3.0, 0.4 * 2.0 + 0.6 * 5.0, 0.5 * 5.0 + 0.5 * 11.0, 1.0], rel=1e-6
    )


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  4. THE EM — the intergenic term is a fixed rate in absolute units
# ═══════════════════════════════════════════════════════════════════════════════════════════════════

_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)


def _locus(n_units, *, gdna_ll, intergenic_ll, deterministic, mode="map", assignment="fractional"):
    """One locus: transcript 0 with ``deterministic`` certain fragments, and ``n_units`` units each with a
    negligible RNA reading, a locus gDNA reading and an intergenic one."""
    em = AbundanceEstimator(1, em_config=EMConfig(mode=mode, assignment_mode=assignment, seed=1))
    em.unambig_counts[0, int(SpliceStrandCol.SPLICED_ANNOT_SENSE)] = deterministic
    part = (
        np.arange(n_units + 1, dtype=np.int64),
        np.zeros(n_units, dtype=np.int32),
        np.full(n_units, -1000.0, dtype=np.float32),
        np.ones(n_units, dtype=np.float32),
        np.full(n_units, _SENSE, dtype=np.uint8),
        np.zeros(n_units, dtype=np.uint8),
        np.full(n_units, gdna_ll, dtype=np.float32),
        np.full(n_units, intergenic_ll, dtype=np.float32),
    )
    return em, part


def _run(em, part, *, gdna_eff_len, rate):
    out = em.run_batch_locus_em_partitioned(
        [part],
        [np.zeros(1, dtype=np.int32)],
        np.zeros(1),
        rna_prior_count=np.zeros(1),
        gdna_eff_len=np.asarray([gdna_eff_len]),
        intergenic_rate=rate,
        em_iterations=100_000,
        em_convergence_delta=1e-13,
    )
    return float(out[1][0]), float(out[2][0]), float(em.gdna_intergenic_em_count)


@pytest.mark.parametrize("n_units", [50, 500])
def test_the_intergenic_term_takes_rate_times_b_over_a_whatever_the_loci_size(n_units):
    """K units each read by the locus's gDNA at ``a = q / L`` per unit of its absolute rate and by the intergenic
    background at ``λ b``: the maximiser gives the background ``λ b / a`` fragments, the locus the rest, whatever
    ``K`` and whatever the locus's other fragments (``M θ_G a + λ b = K a``). A component whose share the M-step
    did not update would give the background a share of ``K`` instead, and the locus's total read with its
    pseudocounts or with the background's mass would move the split."""
    q = b = math.exp(-5.0)
    length, rate = 1000.0, 0.02
    em, part = _locus(n_units, gdna_ll=math.log(q), intergenic_ll=math.log(b), deterministic=100)
    rna, gdna, intergenic = _run(em, part, gdna_eff_len=length, rate=rate)
    expected = rate * b / (q / length)
    assert intergenic == pytest.approx(expected, rel=1e-6)
    assert gdna == pytest.approx(n_units - expected, rel=1e-6)
    assert rna == pytest.approx(0.0, abs=1e-9)


def test_under_the_grouped_prior_the_em_maximises_the_derived_objective():
    """The grouped prior sits on the locus's composition alone: the EM's answer is the maximiser of
    ``K log(M θ_G a + λ b) + D log(M (1 − θ_G)) − M + P_g log θ_G + P_R log(1 − θ_G)``, found here directly. A prior
    pseudocount reaching the intergenic term, or ``M`` read with the pseudocounts in it, moves the split."""
    from scipy.optimize import minimize

    n_units, deterministic, length, rate = 100, 50, 1000.0, 0.02
    q = b = math.exp(-5.0)
    a = q / length
    gdna_prior, rna_prior = 30.0, 70.0

    def negative(x):
        m, theta = math.exp(x[0]), 1.0 / (1.0 + math.exp(-x[1]))
        return -(
            n_units * math.log(m * theta * a + rate * b)
            + deterministic * math.log(m * (1.0 - theta))
            - m
            + gdna_prior * math.log(theta)
            + rna_prior * math.log(1.0 - theta)
        )

    best = minimize(
        negative,
        x0=[math.log(150.0), 0.0],
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20_000},
    )
    m, theta = math.exp(best.x[0]), 1.0 / (1.0 + math.exp(-best.x[1]))
    expected = n_units * rate * b / (m * theta * a + rate * b)

    em, part = _locus(
        n_units, gdna_ll=math.log(q), intergenic_ll=math.log(b), deterministic=deterministic
    )
    out = em.run_batch_locus_em_partitioned(
        [part],
        [np.zeros(1, dtype=np.int32)],
        np.asarray([gdna_prior]),
        rna_prior_count=np.asarray([rna_prior]),
        gdna_eff_len=np.asarray([length]),
        intergenic_rate=rate,
        em_iterations=100_000,
        em_convergence_delta=1e-13,
    )
    assert float(em.gdna_intergenic_em_count) == pytest.approx(expected, rel=1e-6)
    assert float(out[2][0]) == pytest.approx(n_units - expected, rel=1e-6)


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_sampled_assignment_rounds_the_intergenic_term_within_the_locus(mode):
    em, part = _locus(
        50, gdna_ll=-5.0, intergenic_ll=-5.0, deterministic=100, mode=mode, assignment="sample"
    )
    rna, gdna, intergenic = _run(em, part, gdna_eff_len=1000.0, rate=0.02)
    assert intergenic == round(intergenic) and gdna == round(gdna)
    assert rna + gdna + intergenic == 50


def test_with_no_intergenic_place_the_rate_is_inert():
    """A locus none of whose units has an intergenic place is solved exactly as it is without the term."""
    results = []
    for rate in (0.02, 5.0):
        em, part = _locus(50, gdna_ll=-5.0, intergenic_ll=-math.inf, deterministic=100)
        results.append(_run(em, part, gdna_eff_len=1000.0, rate=rate))
    assert results[0] == results[1]
    assert results[0][2] == 0.0


def test_scaling_a_places_weight_and_its_transcripts_length_together_changes_nothing():
    """``(q_t, L_t) → (a q_t, a L_t)`` leaves the transcript's rate ``M θ_t q_t / L_t`` where it was, so the answer must
    not move. A prune of the intergenic term by its score against the transcript's put 14.26 of these 50 units on the
    background at ``a = 1`` and none at ``a = 10⁶`` (the independent review's replay, 2026-10-27)."""
    from rigel.estimator import TranscriptGeometry

    n_units, rate = 50, 0.1
    results, scores = [], []
    for a in (1.0, 1e6):
        scored = _route([[Hit((0,), SpliceType.SPLICED_ANNOT, log_weight=math.log(a)), Hit(())]])
        scores.append(float(scored.log_liks[0]))
        em = AbundanceEstimator(
            1,
            geometry=TranscriptGeometry(
                effective_lengths=np.array([1000.0]), effective_lengths_em=np.array([1000.0 * a])
            ),
            em_config=EMConfig(mode="map", assignment_mode="fractional", seed=1, n_threads=1),
        )
        em.unambig_counts[0, int(SpliceStrandCol.SPLICED_ANNOT_SENSE)] = 100
        part = (
            np.arange(n_units + 1, dtype=np.int64),
            np.zeros(n_units, dtype=np.int32),
            np.full(n_units, scored.log_liks[0], dtype=np.float32),
            np.ones(n_units, dtype=np.float32),
            np.full(n_units, int(SpliceStrandCol.SPLICED_ANNOT_SENSE), dtype=np.uint8),
            np.full(n_units, scored.is_spliced[0], dtype=np.uint8),
            np.full(n_units, scored.gdna_log_liks[0], dtype=np.float32),
            np.full(n_units, scored.intergenic_log_liks[0], dtype=np.float32),
        )
        results.append(_run(em, part, gdna_eff_len=1.0, rate=rate))
    assert 1.0 < results[0][2] < n_units - 1.0
    # Each score is float32: half an ulp of each moves the transcript's rate by that share, and a unit's split by no
    # more. Derived from the representation, not fitted.
    tolerance = 0.5 * float(np.finfo(np.float32).eps) * sum(abs(s) for s in scores)
    assert results[1] == pytest.approx(results[0], rel=tolerance)


# --- an independent reference: the joint fixed point of (θ, M), plain iteration from the equations ---------------


@dataclass(frozen=True)
class _Problem:
    """One locus: per unit each transcript's log-likelihood (``−∞``: no candidate), the locus gDNA's and the
    background's, float32 as the scorer hands them over; each transcript's certain fragments and length; gDNA's
    length; the rate ``λ``."""

    rna: np.ndarray
    gdna: np.ndarray
    background: np.ndarray
    deterministic: np.ndarray
    lengths: np.ndarray
    gdna_length: float
    rate: float


def _problem(name):
    rng = np.random.default_rng(11)
    none = np.float32(-np.inf)

    def ll(n, k=None):
        return rng.uniform(-8.0, -2.0, n if k is None else (n, k)).astype(np.float32)

    lengths, certain = np.array([900.0, 1400.0]), np.array([30.0, 5.0])
    if name == "mixed":
        return _Problem(ll(40, 2), ll(40), ll(40), certain, lengths, 4e3, 0.02)
    if name == "no locus gDNA":  # spliced genic places beside intergenic ones
        return _Problem(ll(30, 2), np.full(30, none), ll(30), certain, lengths, 4e3, 0.02)
    if name == "strong certain RNA":
        return _Problem(ll(20, 2), ll(20), ll(20), np.array([5e3, 2e3]), lengths, 4e3, 0.02)
    if name == "sparse":
        return _Problem(ll(3, 2), ll(3), ll(3), np.array([1.0, 0.0]), lengths, 4e3, 0.02)
    raise KeyError(name)


def _reference(p, *, mode, gdna_prior, rna_prior):
    """``Σ_u`` of the final assignment's RNA, locus gDNA and background posteriors at the fixed point, by plain
    iteration of the model's equations: the E-step at weights ``θ_c / L_c`` (VBEM: ``exp(ψ(α_c) − ψ(Σα)) / L_c``) and
    ``λ / M``, the grouped update on the counts (gDNA its count plus ``P_g``; RNA's total plus ``P_R`` shared by each
    component's count), ``M`` the locus's own data count; the assignment read at ``θ̄ = α / Σα`` under VBEM."""
    from scipy.special import digamma

    q = np.exp(np.column_stack([p.rna, p.gdna]).astype(np.float64))
    b = np.exp(p.background.astype(np.float64))
    L = np.r_[p.lengths, p.gdna_length]
    det = np.r_[p.deterministic, 0.0]
    has_gdna = bool(np.isfinite(p.gdna).any())
    n_t = p.rna.shape[1]

    def grouped(n):
        out = np.zeros_like(n)
        out[n_t] = n[n_t] + gdna_prior if has_gdna else 0.0
        out[:n_t] = (n[:n_t].sum() + rna_prior) * n[:n_t] / n[:n_t].sum()
        return out

    def posteriors(weights, m):
        num = q * weights
        bg = b * p.rate / m
        tot = num.sum(1) + bg
        return num / tot[:, None], bg / tot

    state = np.full(n_t + 1, 1.0 / (n_t + 1)) if mode == "map" else np.ones(n_t + 1)
    m = len(q) + det.sum()
    for _ in range(100_000):
        floored = np.maximum(
            state, np.finfo(np.float64).tiny
        )  # the solver's floor: ψ has a pole at 0
        weights = (
            state / L if mode == "map" else np.exp(digamma(floored) - digamma(state.sum())) / L
        )
        r, _ = posteriors(weights, m)
        counts = det + r.sum(0)
        new = grouped(counts)
        new = new / new.sum() if mode == "map" else new
        # a step of an ulp or so of the largest component: at these rates the error is far below the comparison's
        done = np.abs(new - state).max() <= 1e-15 * new.max() and abs(counts.sum() - m) <= 1e-13 * m
        state, m = new, counts.sum()
        if done:
            break
    else:
        raise AssertionError("the reference did not converge")
    r, r_bg = posteriors(state / state.sum() / L, m)
    return float(r[:, :n_t].sum()), float(r[:, n_t].sum()), float(r_bg.sum())


def _native(p, *, mode, gdna_prior, rna_prior):
    from rigel.estimator import TranscriptGeometry

    n_u, n_t = p.rna.shape
    em = AbundanceEstimator(
        n_t,
        geometry=TranscriptGeometry(effective_lengths=p.lengths, effective_lengths_em=p.lengths),
        em_config=EMConfig(mode=mode, assignment_mode="fractional", seed=1, n_threads=1),
    )
    em.unambig_counts[:, int(SpliceStrandCol.SPLICED_ANNOT_SENSE)] = p.deterministic
    cand = [np.flatnonzero(np.isfinite(p.rna[u])) for u in range(n_u)]
    n_cand = sum(c.size for c in cand)
    part = (
        np.concatenate([[0], np.cumsum([c.size for c in cand])]).astype(np.int64),
        np.concatenate(cand).astype(np.int32),
        np.concatenate([p.rna[u, c] for u, c in enumerate(cand)]),
        np.ones(n_cand, dtype=np.float32),
        np.full(n_cand, _SENSE, dtype=np.uint8),
        (~(np.isfinite(p.gdna) | np.isfinite(p.background))).astype(np.uint8),
        p.gdna,
        p.background,
    )
    out = em.run_batch_locus_em_partitioned(
        [part],
        [np.arange(n_t, dtype=np.int32)],
        np.asarray([gdna_prior]),
        rna_prior_count=np.asarray([rna_prior]),
        gdna_eff_len=np.asarray([p.gdna_length]),
        intergenic_rate=p.rate,
        em_iterations=1_000_000,
        em_convergence_delta=1e-14,
    )
    return float(out[1][0]), float(out[2][0]), float(em.gdna_intergenic_em_count)


@pytest.mark.parametrize("prior", [(0.0, 0.0), (6.0, 20.0)], ids=["no prior", "grouped prior"])
@pytest.mark.parametrize("mode", ["map", "vbem"])
@pytest.mark.parametrize("name", ["mixed", "no locus gDNA", "strong certain RNA", "sparse"])
def test_the_em_reaches_the_reference_fixed_point(name, mode, prior):
    """The shipped solver, SQUAREM over ``(θ, M)`` in one state, lands where plain iteration of the model's equations
    does, under MAP and under the shipped VBEM, with the grouped prior or none: the RNA, locus-gDNA and background
    totals of the final assignment, for a mixed locus, one with no locus gDNA candidate, one dominated by certain RNA
    and a sparse one."""
    p = _problem(name)
    gdna_prior, rna_prior = prior
    got = _native(p, mode=mode, gdna_prior=gdna_prior, rna_prior=rna_prior)
    want = _reference(p, mode=mode, gdna_prior=gdna_prior, rna_prior=rna_prior)
    assert sum(got) == pytest.approx(len(p.rna), rel=1e-12)
    assert got == pytest.approx(want, rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_rate_that_leaves_the_locus_nothing_takes_every_unit(mode):
    """A background that outweighs every composition of a locus with no certain fragment drives ``M`` toward 0, where
    the column's weight ``λ / M`` grows without bound: every unit goes to the background, with no NaN on the way."""
    p = _problem("mixed")
    p = dataclasses.replace(
        p, rna=p.rna - 7.0, gdna=p.gdna - 7.0, deterministic=np.zeros(2), rate=1e3
    )
    got = _native(p, mode=mode, gdna_prior=0.0, rna_prior=0.0)
    assert all(math.isfinite(x) for x in got)
    assert got[2] == pytest.approx(len(p.rna), rel=1e-9)


def test_at_rate_zero_the_intergenic_term_is_inert():
    """A rate of 0 puts nothing on a background the units carry: the locus is solved as it is with no term at all."""
    p = _problem("mixed")
    with_term = _native(
        dataclasses.replace(p, rate=0.0), mode="vbem", gdna_prior=6.0, rna_prior=20.0
    )
    no_place = np.full(len(p.rna), -np.inf, dtype=np.float32)
    without = _native(
        dataclasses.replace(p, background=no_place), mode="vbem", gdna_prior=6.0, rna_prior=20.0
    )
    assert with_term[2] == 0.0
    assert with_term == pytest.approx(without, rel=1e-9)


# ═══════════════════════════════════════════════════════════════════════════════════════════════════
#  5. THE ACCOUNTING — the intergenic term's mass is intergenic gDNA, outside the loci
# ═══════════════════════════════════════════════════════════════════════════════════════════════════


def test_the_intergenic_terms_mass_joins_the_summarys_intergenic_gdna():
    from rigel.cli import _quantification_pools

    stats = PipelineStats()
    stats.n_intergenic = 300
    pools = _quantification_pools(
        stats,
        count_annotated=500.0,
        count_nascent=100.0,
        count_gdna_em=80.0,
        count_gdna_intergenic_em=20.0,
    )
    assert pools["count_gdna_em"] == 80.0
    assert pools["count_gdna"] == 80.0 + 300 + 20.0
    assert pools["fraction_intergenic"] == pytest.approx(320.0 / 1000.0)
    assert pools["fraction_annotated"] + pools["fraction_nascent"] + pools[
        "fraction_gdna"
    ] == pytest.approx(1.0)


def test_a_unit_the_intergenic_term_wins_is_stamped_intergenic_gdna():
    """The EM's winner code for the intergenic term reaches the annotated BAM as ``AF_GDNA_INTERGENIC``, the stamp a
    template no hit of which resolves already carries; the locus's own gDNA keeps ``AF_GDNA_EM``."""
    import pandas as pd

    from rigel.annotate import AF_ANNOTATED, AF_GDNA_EM, AF_GDNA_INTERGENIC
    from rigel.native import EM_WINNER_INTERGENIC
    from rigel.pipeline import _populate_em_annotations

    class _Part:
        frag_ids = np.asarray([7, 8, 9], dtype=np.int64)
        frag_class = np.full(3, FRAG_MULTIMAPPER, dtype=np.int8)
        locus_id = 0

    class _Ann:
        def add_batch(self, **kw):
            self.kw = kw

    class _Idx:
        t_to_g_arr = np.zeros(1, dtype=np.int64)
        t_df = pd.DataFrame({"is_synthetic": [False]})

    ann = _Ann()
    _populate_em_annotations(
        [_Part()],
        np.asarray([0, -2, int(EM_WINNER_INTERGENIC)], dtype=np.int32),
        np.ones(3, dtype=np.float32),
        np.full(3, 2, dtype=np.int16),
        ann,
        _Idx(),
    )
    assert ann.kw["tx_flags"].tolist() == [AF_ANNOTATED, AF_GDNA_EM, AF_GDNA_INTERGENIC]
    assert ann.kw["best_tids"].tolist() == [0, -1, -1]


def test_the_annotated_bam_writes_a_multimapper_called_intergenic_at_no_one_place(
    toy_index, tmp_path
):
    """Through the whole pipeline: in a library whose unique unspliced hits show no mismatch, a multimapper whose genic
    place has one and whose intergenic one has none is counted intergenic, never gated, and written ``ZF 0x25`` on
    every record with no hit marked the winning placement (``ZH 1``): the call is the fragment's, at no one place, as a
    locus's gDNA winner's is."""
    from rigel.annotate import AF_GDNA_INTERGENIC
    from rigel.config import PipelineConfig
    from rigel.pipeline import run_pipeline

    records = [
        r
        for i in range(60)
        for r in _template(f"u{i:03d}", (1010 + 3 * i, "100M", 1110 + 3 * i, "100M"))
    ]
    records += _template("zz_case5", GENIC, INTERGENIC, nm=(1, 0))
    annotated = tmp_path / "annotated.bam"
    result = run_pipeline(
        _write_bam(tmp_path, toy_index, records),
        toy_index,
        PipelineConfig(scan=BamScanConfig(total_threads=1), annotated_bam_path=annotated),
    )
    assert result.stats.mismatch_rate == 0.0
    assert result.stats.n_intergenic == 1 and result.stats.n_gated_out == 0
    with pysam.AlignmentFile(str(annotated), "rb", check_sq=False) as f:
        out = [r for r in f.fetch(until_eof=True) if r.query_name == "zz_case5"]
    assert len(out) == 4
    assert {r.get_tag("ZF") for r in out} == {AF_GDNA_INTERGENIC}
    assert {r.get_tag("ZH") for r in out} == {0}
