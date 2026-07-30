"""The same BAM must give the same tally at any worker count — bit for bit.

    Plan: ``docs/IMPLEMENTATION_PLAN.md`` §0 step 8   ·   Design: ``ACCUMULATOR_DESIGN.md`` §10.1

⭐ THIS IS NEWLY ACHIEVABLE, and it is the reason every channel is an integer. The float accumulator this
replaced differed by only ~3.7e-7 per cell between worker counts — and that propagated to a **~2.6 %
difference in the calibration output**, so the same BAM gave different answers on different machines.
Integer addition is associative, so the per-worker merge is now exact.

⚠ WHY IT IS NOT ENOUGH TO TEST ``Accumulator.merge_from`` DIRECTLY (the parity module does that). This
exercises the *scanner's* worker path: each worker builds its own ``AccumulatorSet`` from the scanner's
members, and the chunk→worker split is data-dependent. A worker whose set was constructed differently —
junctions installed on the template but not on the copies, say — is invisible to a merge test on one
accumulator and shows up here.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel import scan_payload
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario


SEED = 20260730


@pytest.fixture
def oracle(tmp_path):
    """A scenario built so that **every** bank of the tally receives something.

    ⚠ The obvious scenario does not. Two single-isoform genes deposit into the nodes and the junction
    edges and leave **both contiguous-edge banks identically zero**: every cut is an exon boundary, so a
    mature fragment either fits inside an exon (contained) or splices across the gap (junction edge), and
    it never has bases on both sides of a line. A bit-identity gate over an all-zero array passes for the
    wrong reason — this project has already had one report "32/32 IDENTICAL" on an arm with zero rows. So:

    * ``t2`` starts at 500, **inside** ``t1``'s first exon. That makes 500 a cut, and a ``t1`` fragment
      spanning it has bases on both sides — a contiguous crossing. If that fragment also uses the
      junction, it is a crossing in the SPLICED bank, which is the channel the old design merged away.
    * ``t4`` ends at 650, which cuts a **50 bp** node out of ``[500, 700)``. Nothing else here can be
      SPANNED: spanning needs one segment covering a node whole, so at 200 bp nodes and 220 bp fragments
      it essentially never happens, and mature RNA can never span the node before a junction at all —
      it has no base past the exon end, it splices there. A 50 bp node is spanned by ordinary gDNA.
    * ``gdna_fraction`` puts genomic fragments in the intronic and intergenic nodes, which is the
      unspliced bank and the two pure gDNA length pools.

    ⚠ Every one of those three was added because the assertion below found the array empty. That is the
    fixture doing its job: each is a population the tally has, and a determinism gate that never sees one
    is not testing it.
    """
    scenario = Scenario(
        "worker_determinism",
        genome_length=6000,
        seed=SEED,
        work_dir=tmp_path / "worker_determinism",
    )
    scenario.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 700), (1200, 1700)], "abundance": 60},
            {"t_id": "t2", "exons": [(500, 700), (1200, 1700)], "abundance": 25},
            {"t_id": "t4", "exons": [(200, 650), (1200, 1700)], "abundance": 15},
        ],
    )
    scenario.add_gene(
        "g2",
        "-",
        [{"t_id": "t3", "exons": [(3000, 3400), (4000, 4400)], "abundance": 40}],
    )
    result = scenario.build_oracle(
        n_fragments=4000,
        gdna_fraction=0.3,
        sim_config=ReadSimConfig(
            frag_mean=220,
            frag_std=60,
            frag_min=80,
            frag_max=600,
            read_length=100,
            strand_specificity=1.0,
            seed=SEED,
        ),
    )
    yield result
    scenario.cleanup()


def _raw_tally(result, n_workers: int, monkeypatch) -> dict:
    """Scan at ``n_workers`` and return the accumulator's raw payload dict.

    ⚠ The monkeypatch is a TEST-LOCAL adapter around a boundary that is mid-move: S3 changed the payload
    keys and S4 rewrites ``AccumulatorPayload`` to read them, so ``from_scan_result`` cannot parse this
    scan yet. Deleting these two lines at S4 and comparing payload attributes instead is the whole of the
    follow-up. It is deliberately here rather than in the pipeline: production code does not get a shim to
    make a test run.
    """
    monkeypatch.setattr(
        scan_payload.AccumulatorPayload,
        "from_scan_result",
        classmethod(lambda cls, scan_result: scan_result["calibration"]),
    )
    _, _, _, _, payload = scan_and_buffer(
        str(result.bam_path),
        result.index,
        BamScanConfig(sj_strand_tag="auto", total_threads=n_workers),
    )
    assert payload is not None, "the scan produced no accumulator payload"
    return payload


def test_the_tally_is_bit_identical_at_1_2_4_and_8_workers(oracle, monkeypatch):
    baseline = _raw_tally(oracle, 1, monkeypatch)

    # A tally of zeros would satisfy bit-identity trivially, so the baseline has to be shown to contain
    # something first. ⚠ A bit-identity gate in this project has already lied in exactly this way: an arm
    # with ZERO rows scored "32/32 IDENTICAL" because the comparison looped over the empty arm's rows.
    for key, why in [
        ("node_start_count", "nothing was deposited at all"),
        ("node_contained_count", "no fragment fitted inside a node"),
        ("node_spanning_count", "no fragment covered a node whole"),
        (
            "edge_unspliced_count",
            "no unspliced contiguous crossing — the mixture being deconvolved",
        ),
        ("edge_spliced_count", "no SPLICED contiguous crossing — the certified-RNA channel"),
        ("sj_count", "no annotated junction was used"),
        ("pool_lengths", "no fragment entered a length pool"),
    ]:
        assert int(np.asarray(baseline[key]).sum()) > 0, f"{key}: {why}"

    array_keys = [key for key, value in baseline.items() if isinstance(value, np.ndarray)]
    assert len(array_keys) >= 12, (
        f"only {len(array_keys)} arrays in the payload; the gate is too narrow"
    )

    for n_workers in (2, 4, 8):
        other = _raw_tally(oracle, n_workers, monkeypatch)
        assert set(other) == set(baseline), f"{n_workers} workers: different payload keys"
        for key in array_keys:
            expected = np.asarray(baseline[key])
            actual = np.asarray(other[key])
            assert actual.dtype == expected.dtype, f"{n_workers} workers: {key} dtype"
            assert np.array_equal(actual, expected), (
                f"{n_workers} workers: {key} is not bit-identical to the single-worker run — "
                f"{int(np.count_nonzero(actual != expected))} of {expected.size} cells differ"
            )
        assert dict(other["qc"]) == dict(baseline["qc"]), f"{n_workers} workers: qc"


def test_the_start_count_invariant_holds_on_a_real_scan(oracle, monkeypatch):
    """``sum(node_start_count) == deposited`` — the accumulator's one non-tautological invariant.

    The three "conservation identities" it replaced could only be evaluated by re-running the deposit, so a
    deliberately broken replay satisfied all three while 91 % of the crossings were junk. This one is
    checkable against a number the deposit counts independently.
    """
    tally = _raw_tally(oracle, 1, monkeypatch)
    assert int(np.asarray(tally["node_start_count"]).sum()) == int(tally["qc"]["deposited"])
