"""A hit's sequence: what tells a fragment's alignments apart beyond where they sit.

Each read base of a hit is aligned or soft-clipped. At the library's mismatch rate ε an aligned base reads its
reference base with probability 1 − ε and each of the three others with ε/3; a clipped base is placed nowhere and reads
at the composition of the four bases, ¼. Against an aligned base that matches, a mismatch therefore costs
log[ε / (3(1 − ε))] and a clipped base log[1 / (4(1 − ε))] (``rigel.scoring.sequence_log_penalties``).

* ε is the library's own: NM over the read bases (M, =, X and I) of its unique unspliced hits. A junction's alignment
  is the aligner's choice of where to cut, and a multimapper's hits are the question the prices answer.
* A hit carries its records' NM and soft clips, both mates'.
* The prices are read only on a multimapper's hits, which they tell apart: one hit's sequence is common to all its
  candidates and cancels (``tests/test_pipeline_routing.py::test_a_one_hit_fragment_reads_no_sequence``).
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest
from test_multimapper_places import (
    CLIP_LP,
    GENIC,
    GENIC_SPLICED,
    INTERGENIC,
    MISMATCH_LP,
    Hit,
    _build_index,
    _route,
    _scan_chunks,
    _template,
    _term,
    _write_bam,
)

from rigel.scoring import sequence_log_penalties

BASES = "ACGT"


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    """One gene, tA's exons at [1000, 1300) and [2000, 2300) on chr1 (``test_multimapper_places``)."""
    return _build_index(tmp_path_factory.mktemp("sequence"), ["chr1"])


# ── the prices ──────────────────────────────────────────────────────────────────────────────────────────────────────


def _emissions(eps):
    """P(read base | reference base), by enumeration: an aligned base, and a clipped one."""
    aligned = {(r, b): (1.0 - eps) if b == r else eps / 3.0 for r in BASES for b in BASES}
    clipped = {(r, b): 0.25 for r in BASES for b in BASES}
    return aligned, clipped


@pytest.mark.parametrize("eps", [1e-4, 0.002, 0.01, 0.1, 0.5, 0.7])
def test_the_prices_are_the_emission_ratios(eps):
    aligned, clipped = _emissions(eps)
    for r in BASES:  # each a distribution over the read base
        assert sum(aligned[r, b] for b in BASES) == pytest.approx(1.0)
        assert sum(clipped[r, b] for b in BASES) == pytest.approx(1.0)
    mismatch, clip = sequence_log_penalties(eps)
    for r, b in itertools.product(BASES, BASES):
        if b != r:
            assert mismatch == pytest.approx(math.log(aligned[r, b] / aligned[r, r]), rel=1e-12)
        assert clip == pytest.approx(math.log(clipped[r, b] / aligned[r, r]), rel=1e-12)


def test_a_library_with_no_mismatch_makes_one_impossible_and_one_with_no_unique_unspliced_hit_prices_nothing():
    assert sequence_log_penalties(0.0) == (-math.inf, -math.log(4.0))
    assert sequence_log_penalties(None) == (0.0, 0.0)


@pytest.mark.parametrize("eps", [-1e-9, 0.75, 0.9, math.nan])
def test_a_rate_no_library_has_is_refused(eps):
    with pytest.raises(ValueError, match="no library's"):
        sequence_log_penalties(eps)


# ── the scan: a hit's clips, and the library's rate ────────────────────────────────────────────────────────────────


def _clips(chunks):
    return [int(c) for chunk in chunks for c in chunk.clip]


def test_a_hit_carries_the_soft_clips_of_all_its_records(index, tmp_path):
    """Both ends of both mates: a unique hit, and a multimapper's two hits, each its own."""
    records = _template("a_unique", (1050, "10S90M", 1150, "92M8S"))
    records += _template(
        "b_multimapper", (1050, "3S95M2S", 1150, "100M"), (8050, "100M", 8150, "4S96M")
    )
    _stats, chunks = _scan_chunks(tmp_path, index, records)
    assert _clips(chunks) == [10 + 8, 3 + 2, 4]


def test_the_mismatch_rate_reads_the_unique_unspliced_hits_alone(index, tmp_path):
    """NM over the read bases a hit aligns — M and I, never a deletion's or a clip's — summed over the unique
    unspliced hits: a spliced hit and a multimapper's hits stay out, an intergenic unique hit comes in."""
    records = _template("a_unique", (1050, "5S45M2I48M", 1150, "50M3D50M"), nm=(2,))
    records += _template("b_intergenic", INTERGENIC, nm=(1,))
    records += _template("c_spliced", GENIC_SPLICED, nm=(3,))
    records += _template("d_multimapper", GENIC, INTERGENIC, nm=(4, 4))
    stats, _chunks = _scan_chunks(tmp_path, index, records)
    assert (
        stats.n_unique_unspliced_mismatches == 2 * 2 + 2 * 1
    )  # each record of a place carries the place's NM
    assert stats.n_unique_unspliced_bases == (45 + 2 + 48) + (50 + 50) + 2 * 100
    assert (
        stats.mismatch_rate == stats.n_unique_unspliced_mismatches / stats.n_unique_unspliced_bases
    )


def test_a_library_with_no_unique_unspliced_hit_has_no_rate(index, tmp_path):
    stats, _chunks = _scan_chunks(tmp_path, index, _template("s", GENIC_SPLICED, nm=(1,)))
    assert stats.n_unique_unspliced_bases == 0 and stats.mismatch_rate is None


# ── the scorer: a multimapper's hits ────────────────────────────────────────────────────────────────────────────────


def _rna(em):
    return dict(zip(em.t_indices.tolist(), em.log_liks.tolist()))


@pytest.mark.parametrize("swap", [False, True])
def test_a_multimappers_hits_differ_by_their_sequence(swap):
    """Two hits alike but for their sequence, each its own place: every candidate at a hit carries its sequence price,
    RNA's and gDNA's alike, so the transcripts differ by it and gDNA sums the two prices' terms. One mismatch and one
    clipped base, 8.7 nats, stay within the pruning bound."""
    plain, dirty = Hit((0,), nm=0, clip=0), Hit((1,), nm=1, clip=1)
    em = _route([[dirty, plain] if swap else [plain, dirty]])
    price = MISMATCH_LP + CLIP_LP
    assert _rna(em)[1] - _rna(em)[0] == pytest.approx(price, abs=1e-5)
    assert float(em.gdna_log_liks[0]) == pytest.approx(
        np.logaddexp(_term(250), _term(250) + price), abs=1e-5
    )


def test_the_hit_that_explains_every_base_is_the_more_likely_place():
    """The splice artifact's case: a spliced hit that aligns every base beside a contiguous one that clips the bases
    it cannot place. The contiguous reading's candidates pay ¼ per clipped base, and only they do."""
    em = _route([[Hit((0,), start=1000, clip=0), Hit((1,), start=1000, clip=3)]])
    assert _rna(em)[0] - _rna(em)[1] == pytest.approx(-3 * CLIP_LP, abs=1e-5)


# ── the pipeline: the scan's rate prices the scorer ─────────────────────────────────────────────────────────────────


def test_the_pipeline_prices_a_hit_at_the_librarys_own_rate(index, tmp_path, monkeypatch):
    import rigel.pipeline as PL
    from rigel.config import BamScanConfig, PipelineConfig

    records = [
        r
        for i in range(40)
        for r in _template(f"u{i:03d}", (1010 + 3 * i, "100M", 1110 + 3 * i, "100M"), nm=(i % 3,))
    ]
    records += _template("zz_multimapper", GENIC, INTERGENIC, nm=(1, 0))
    seen = {}
    shipped = PL.build_scorer

    def spy(*args, **kw):
        seen.update(kw)
        return shipped(*args, **kw)

    monkeypatch.setattr(PL, "build_scorer", spy)
    result = PL.run_pipeline(
        _write_bam(tmp_path, index, records),
        index,
        PipelineConfig(scan=BamScanConfig(total_threads=1)),
    )
    eps = sum(2 * (i % 3) for i in range(40)) / (40 * 200)
    assert result.stats.mismatch_rate == pytest.approx(eps, rel=1e-12)
    mismatch, clip = sequence_log_penalties(eps)
    assert seen["mismatch_log_penalty"] == pytest.approx(mismatch, rel=1e-12)
    assert seen["clip_log_penalty"] == pytest.approx(clip, rel=1e-12)
