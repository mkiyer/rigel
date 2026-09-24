"""The EM's two pseudocounts: calibration's gDNA count, restated at the locus's own scale.

The locus EM's theta is each component's share of EVERY fragment in the locus: the E-step units and the
deterministic spliced fragments, which skip the E-step but are added to every M-step
(``em_solver.cpp:map_em_step`` / ``vbem_step``, ``raw = unambig_totals + em_totals``). The grouped prior update
is exact EM for ``P_g log theta_g + P_R log(1 - theta_g)``, so the prior leaves the EM's own answer in place
iff ``P_g : P_R = G : (N - G)`` over that whole pool. Calibration's unspliced split ``G : U_R`` states the
unspliced fragments' split as the whole pool's and leans toward gDNA in proportion to the spliced share.

The falsification test drives the SHIPPED solver through the pipeline's own locus runner on a locus whose
counts are exact expectations, so the prior-free fixed point is the truth, and hands it calibration's exact gDNA
count: gDNA and every transcript must land on the truth whatever the spliced share and however the spliced
fragments split between the E-step and the deterministic path.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.locus import Locus, MultiLocus
from rigel.locus_partition import partition_and_free
from rigel.pipeline import _run_locus_em_partitioned, em_pseudocounts
from rigel.scored_fragments import ScoredFragments
from rigel.splice import SpliceStrandCol

# ── the rule ─────────────────────────────────────────────────────────────────────────────────────


def test_the_odds_are_calibrations_gdna_share_of_the_fragments_it_counted():
    pg, pr = em_pseudocounts(np.array([300.0]), np.array([1000.0]), np.array([600.0]))
    assert pg[0] / (pg[0] + pr[0]) == pytest.approx(0.3)


def test_a_fragment_calibration_never_saw_carries_the_locus_share_not_rna():
    """Calibration deposits no multimapper, so its count covers the fragments it saw. Reading that count
    against the whole locus instead would call every unseen fragment RNA: here gDNA is 3/10 of what
    calibration saw, and the prior must say 3/10 whether or not the locus also holds unseen fragments."""
    seen_only = em_pseudocounts(np.array([300.0]), np.array([1000.0]), np.array([1000.0]))
    with_unseen = em_pseudocounts(np.array([300.0]), np.array([1000.0]), np.array([1500.0]))
    for pg, pr in (seen_only, with_unseen):
        assert pg[0] / (pg[0] + pr[0]) == pytest.approx(0.3)


def test_the_strength_is_one_pseudocount_per_gdna_eligible_unit():
    pg, pr = em_pseudocounts(
        np.array([300.0, 5.0]), np.array([1000.0, 50.0]), np.array([600.0, 40.0])
    )
    np.testing.assert_allclose(pg + pr, [600.0, 40.0])


def test_a_gdna_count_above_the_locus_total_is_a_share_of_one_never_a_negative_rna_pseudocount():
    pg, pr = em_pseudocounts(np.array([1200.0]), np.array([1000.0]), np.array([600.0]))
    assert (pg[0], pr[0]) == (600.0, 0.0)


def test_a_locus_with_no_gdna_eligible_unit_gets_no_prior():
    pg, pr = em_pseudocounts(np.array([80.0]), np.array([100.0]), np.array([0.0]))
    assert (pg[0], pr[0]) == (0.0, 0.0)


def test_a_locus_calibration_counted_none_of_gets_no_prior():
    pg, pr = em_pseudocounts(np.array([0.0]), np.array([0.0]), np.array([40.0]))
    assert (pg[0], pr[0]) == (0.0, 0.0)


# ── the falsification test: neutral on the shipped solver ────────────────────────────────────────

# Unit-length bins. exon1 intron1 exon2 intron2 exon3 are genomic; j12 j23 j13 are junction start positions
# (spliced only). t0 = exon1-j12-exon2-j23-exon3, t1 = exon1-j13-exon3, t2 = exon1-j12-exon2, and t3 is the
# synthetic span exon1..exon3 (unspliced RNA).
_SUPPORT = {
    0: ("exon1", "j12", "exon2", "j23", "exon3"),
    1: ("exon1", "j13", "exon3"),
    2: ("exon1", "j12", "exon2"),
    3: ("exon1", "intron1", "exon2", "intron2", "exon3"),
}
_GDNA = ("exon1", "intron1", "exon2", "intron2", "exon3")
_JUNCTIONS = ("j12", "j23", "j13")
# Each component's fragments per bin: uniform, so the model is exact.
_PER_BIN = {0: 20, 1: 10, 2: 10, 3: 10}
_GDNA_PER_BIN = 10


def _bins(junction_len: int) -> dict[str, int]:
    return {
        "exon1": 10,
        "intron1": 30,
        "exon2": 10,
        "intron2": 30,
        "exon3": 10,
        "j12": junction_len,
        "j23": junction_len,
        "j13": 2 * junction_len,
    }


def _toy_locus(ss: float, junction_len: int):
    """Exact expected counts: returns (ScoredFragments, deterministic per transcript, lengths, truth)."""
    bins = _bins(junction_len)
    n_t = len(_SUPPORT)
    offsets, cand_t, cand_ll, is_spliced, gdna_ll = [0], [], [], [], []
    deterministic = np.zeros(n_t)
    truth_gdna = 0.0

    def unit(cands, ll, spliced, g):
        cand_t.extend(cands)
        cand_ll.extend([ll] * len(cands))
        offsets.append(len(cand_t))
        is_spliced.append(spliced)
        gdna_ll.append(g)

    for b, nb in bins.items():
        cands = [t for t in range(n_t) if b in _SUPPORT[t]]
        spliced = b in _JUNCTIONS
        g_here = _GDNA_PER_BIN if (b in _GDNA and not spliced) else 0
        for _ in range(nb):
            n_rna = sum(_PER_BIN[t] for t in cands)
            if spliced:
                if len(cands) == 1:
                    # one candidate and an annotated splice: never in the E-step
                    deterministic[cands[0]] += n_rna
                else:
                    for _k in range(n_rna):
                        unit(cands, 0.0, True, -np.inf)
                continue
            truth_gdna += g_here
            for sense, n in (
                (True, n_rna * ss + g_here * 0.5),
                (False, n_rna * (1 - ss) + g_here * 0.5),
            ):
                assert abs(n - round(n)) < 1e-9
                for _k in range(int(round(n))):
                    unit(cands, np.log(ss if sense else 1 - ss), False, np.log(0.5))

    n_units, n_cand = len(is_spliced), len(cand_t)
    em = ScoredFragments(
        offsets=np.array(offsets, np.int64),
        t_indices=np.array(cand_t, np.int32),
        log_liks=np.array(cand_ll, np.float64),
        count_cols=np.zeros(n_cand, np.uint8),
        coverage_weights=np.ones(n_cand, np.float64),
        locus_t_indices=np.array([cand_t[offsets[u]] for u in range(n_units)], np.int32),
        locus_count_cols=np.zeros(n_units, np.uint8),
        is_spliced=np.array(is_spliced, bool),
        gdna_log_liks=np.array(gdna_ll, np.float64),
        frag_ids=np.arange(n_units, dtype=np.int64),
        frag_class=np.zeros(n_units, np.int8),
        splice_type=np.zeros(n_units, np.uint8),
        n_units=n_units,
        n_candidates=n_cand,
    )
    t_len = np.array([sum(bins[b] for b in _SUPPORT[t]) for t in range(n_t)], np.float64)
    truth_t = np.array([_PER_BIN[t] * t_len[t] for t in range(n_t)])
    gdna_len = float(sum(bins[b] for b in _GDNA))
    return em, deterministic, t_len, gdna_len, truth_gdna, truth_t


def _run(ss: float, junction_len: int, mode: str):
    em, deterministic, t_len, gdna_len, truth_gdna, truth_t = _toy_locus(ss, junction_len)
    n_t = len(t_len)
    est = AbundanceEstimator(
        n_t, em_config=EMConfig(seed=1, mode=mode, assignment_mode="fractional", n_threads=1)
    )
    est._t_eff_len_em = t_len
    est.unambig_counts[:, int(SpliceStrandCol.SPLICED_ANNOT_SENSE)] = deterministic
    locus = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.arange(n_t, dtype=np.int32),
        unit_indices=np.arange(em.n_units, dtype=np.int32),
        gdna_span=int(gdna_len),
        loci=(Locus(ref="chr1", ref_id=0, start=0, end=int(gdna_len)),),
    )
    index = SimpleNamespace(
        t_to_g_arr=np.zeros(n_t, np.int64), g_df=pd.DataFrame({"is_synthetic": [False]})
    )
    _run_locus_em_partitioned(
        est,
        partition_and_free(em, [locus]),
        [locus],
        index,
        gdna_count=np.array([truth_gdna]),  # calibration exact
        gdna_eff_len=np.array([gdna_len]),
        em_config=EMConfig(mode=mode, iterations=20000, convergence_delta=1e-12),
    )
    t = est.unambig_counts.sum(axis=1) + est.em_counts.sum(axis=1)
    return est.gdna_em_count, truth_gdna, t, truth_t


@pytest.mark.parametrize("mode", ["map", "vbem"])
@pytest.mark.parametrize("ss", [0.9, 0.5])
@pytest.mark.parametrize("junction_len", [5, 20, 60])
def test_an_exact_calibration_leaves_gdna_and_every_transcript_at_the_truth(mode, ss, junction_len):
    """The spliced share runs from 11 % to 60 % of the locus. Calibration's unspliced split ``G : U_R`` over-calls
    gDNA here by +70 to +944 fragments against a truth of 900, growing with the spliced share, and moves the
    transcripts by as much (measured on the pre-change solver path, 2026-09-24)."""
    gdna, truth_gdna, t, truth_t = _run(ss, junction_len, mode)
    assert gdna == pytest.approx(truth_gdna, abs=0.5), f"gDNA {gdna:.2f} against {truth_gdna:.0f}"
    if mode == "map":
        np.testing.assert_allclose(t, truth_t, atol=0.5)
    else:
        # VBEM moves the split between two isoforms that share positions by O(1) fragment at these counts
        # whatever the prior (``exp ψ(α) ≈ α − ½``): 1.1–1.8 of 250–800 here. The shipped odds moved them by
        # 67–900.
        np.testing.assert_allclose(t, truth_t, rtol=0.01)
