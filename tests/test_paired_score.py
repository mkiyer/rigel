"""Gate for ``scripts/design/paired_score.py``: the paired proper score, recorded through the real pipeline.

The instrument's ``--self-test`` plants each defect in hand-built arrays. This gate runs the record step end to end on
a small simulated library, so the hooks, the annotated-BAM join and the closure are exercised on a real run: the class
sums close on the pools ``quant_accuracy.score_library`` reports, the hooks change no count, every read pair is
accounted for, a pinned arm re-recorded pairs with itself at a difference of exactly zero, and a class swapped or a
read name mis-joined on real records is refused.
"""

from __future__ import annotations

import numpy as np
import pytest

import paired_score as PS
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

#: Named as a panel condition is, so the record and the pair read its stratum from it.
CONDITION = "toy_ss_0.99_capture_off"


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """Three genes with staggered isoforms, gDNA at 1x the RNA fragment count and nascent RNA in the introns."""
    sc = Scenario("ps", genome_length=9000, seed=17, work_dir=tmp_path_factory.mktemp("ps_sim"))
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 30},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 40}])
    sc.add_gene(
        "g3",
        "+",
        [
            {"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 25},
            {"t_id": "t3b", "exons": [(6600, 7100), (7400, 7440), (7800, 8300)], "abundance": 12},
        ],
    )
    return sc.build_oracle(
        n_rna_fragments=3000,
        gdna_fraction=1.0,
        nascent_abundance=15.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            read_length=90,
            strand_specificity=0.99,
            seed=17,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


def _record(toy, tmp_path, settings=()):
    cfg = PS.protocol_config(list(settings))
    return PS.record_bam("arm", toy.bam_path, toy.index, cfg, tmp_path / "work", {}, CONDITION)


@pytest.fixture(scope="module")
def shipped(toy, tmp_path_factory):
    return _record(toy, tmp_path_factory.mktemp("ps_shipped"))


def _as_record(arrays, row, arm):
    return {"meta": {**row, "arm": arm}, **arrays}


def test_a_real_run_is_accepted_and_its_class_sums_close_on_the_pools(shipped):
    arrays, row = shipped
    assert row["accepted"], row["checks"]
    assert set(row["closure"]) == {"annotated", "nascent", "gdna", "intergenic"}
    for name, c in row["closure"].items():
        assert c["diff"] <= c["bound"], (name, c)
        assert c["pool"] > 0, f"the toy holds no {name} mass, so its closure checks nothing"
    assert row["checks"]["ledger_mismatches"] == 0 and row["checks"]["zf_not_the_tables_flag"] == 0
    assert row["counts"]["fast_path_intergenic"] > 0 and row["counts"]["deterministic"] > 0
    assert row["counts"]["em_units"] > 0
    assert row["settings"] == {"scan.total_threads": 1, "em.assignment_mode": "fractional"}


def test_every_read_pair_has_one_status_and_only_the_forecast_ones_are_scored(shipped):
    arrays, row = shipped
    n = arrays["key"].size
    assert sum(sum(v.values()) for v in row["status"].values()) == n
    assert sum(row["in_bam"].values()) == n
    forecast = arrays["status"] == PS.FORECAST
    assert row["scores"]["n_forecast"] == int(forecast.sum())
    assert np.all(arrays["mass"][forecast].sum(1) > 0) and np.all(arrays["mass"][~forecast] == 0)


def test_the_hooks_change_no_count(toy, tmp_path):
    """The record's run writes the annotated BAM and asks the EM for its class masses; neither may move a number."""
    cfg = PS.protocol_config([])
    hooked, units, mass, _ = PS.capture_run(toy.bam_path, toy.index, cfg, tmp_path / "a.bam")
    plain = run_pipeline(str(toy.bam_path), toy.index, cfg)
    a, b = hooked.estimator, plain.estimator
    np.testing.assert_array_equal(
        a.get_counts_df(toy.index)["count"], b.get_counts_df(toy.index)["count"]
    )
    assert (a.gdna_em_count, a.synthetic_count, a.gdna_intergenic_em_count) == (
        b.gdna_em_count,
        b.synthetic_count,
        b.gdna_intergenic_em_count,
    )
    assert mass.shape == (units["fid"].size, 4)


def test_a_pinned_arm_re_recorded_pairs_with_itself_at_exactly_zero(toy, shipped, tmp_path):
    again = _record(toy, tmp_path)
    row = PS.pair_records(_as_record(*shipped, "a"), _as_record(*again, "a2"))
    assert row["delta"]["brier"] == 0.0 and row["delta"]["mcb"] == 0.0
    assert sum(row["moved"].values()) == 0
    assert sum(row["sets"]["both"].values()) == row["common"]["a"]["n_forecast"]


def test_two_settings_pair_on_every_read_name(toy, shipped, tmp_path):
    silent = _record(toy, tmp_path, ["calibration.message_policy=silent"])
    row = PS.pair_records(_as_record(*shipped, "shipped"), _as_record(*silent, "silent"))
    assert row["varies"] == {"calibration.message_policy": [None, "silent"]}
    total = sum(sum(row["sets"][s].values()) for s in ("both", "a_only", "b_only", "neither"))
    assert total == shipped[0]["key"].size
    assert row["stamps"]["a"]["revision"] == row["stamps"]["b"]["revision"]


def test_the_record_files_round_trip_and_a_changed_file_is_refused(shipped, tmp_path):
    arrays, row = shipped
    PS.write_record(tmp_path, CONDITION, arrays, row)
    back = PS.load_record(tmp_path, CONDITION)
    for k in ("key", "kind", "status", "mass"):
        np.testing.assert_array_equal(back[k], arrays[k])
    # another run's arrays under this record's row
    PS.write_record(tmp_path / "other", CONDITION, {**arrays, "mass": arrays["mass"][::-1]}, row)
    (tmp_path / "other" / f"{CONDITION}.npz").replace(tmp_path / f"{CONDITION}.npz")
    with pytest.raises(SystemExit, match="digest"):
        PS.load_record(tmp_path, CONDITION)


def test_a_swapped_class_on_a_real_run_is_refused(toy, tmp_path, monkeypatch):
    """Break the fixed code: book the annotated column as gDNA and the gDNA column as annotated. The closure and the
    winner check must fire."""
    ann, gdna = PS.ANN, PS.GDNA
    monkeypatch.setattr(PS, "ANN", gdna)
    monkeypatch.setattr(PS, "GDNA", ann)
    _, row = _record(toy, tmp_path)
    assert not row["accepted"]
    assert row["checks"]["closure_over_bound"] >= 2
    assert row["checks"]["winner_above_its_class"] > 0


def _plant_capture(monkeypatch):
    """Plant a capture reading on the capture-OFF toy, whose reader reads none: the reader's own answer with half its
    regions at efficiency 3, saturated at 3, at weight 1. It moves the toy's pools, so a record that never reached
    the ``no_capture`` arm reads it."""
    import importlib

    calibrate = importlib.import_module("rigel.calibration.calibrate")
    inner = calibrate._capture_fields

    def planted(*args, **kwargs):
        fields = inner(*args, **kwargs)
        eff = np.array(fields["gdna_capture_efficiency_region"], dtype=np.float64)
        eff[: eff.size // 2] = 3.0
        return dict(
            fields,
            gdna_capture_efficiency_region=eff,
            gdna_capture_saturation=3.0,
            gdna_capture_weight=1.0,
        )

    monkeypatch.setattr(calibrate, "_capture_fields", planted)


def test_a_no_capture_record_reads_no_capture_and_closes_on_quant_accuracys_no_capture_pools(
    toy, tmp_path, monkeypatch
):
    """Every capture-ON row is read beside ``no_capture``, so a record can run under it: every capture efficiency
    read as 1, through ``quant_accuracy.install_no_capture`` and refused by its ``check_no_capture``. Under a planted
    capture reading the shipped record reads the plant and the ``no_capture`` record does not: its pools are
    ``quant_accuracy.py --arm no_capture``'s on the same run, bit for bit, and the pair names the arm as what varies."""
    qa = PS.sibling("quant_accuracy.py")
    cfg = PS.protocol_config([])
    restore, _ = qa.install_no_capture()
    try:
        plain = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    want = qa.score_library(plain, plain.estimator.get_counts_df(toy.index), {})

    _plant_capture(monkeypatch)
    base = PS.record_bam("base", toy.bam_path, toy.index, cfg, tmp_path / "work", {}, CONDITION)
    nc = PS.record_bam(
        "nc", toy.bam_path, toy.index, cfg, tmp_path / "work", {}, CONDITION, qa_arm="no_capture"
    )
    for arrays, row in (base, nc):
        assert row["accepted"], row["checks"]
    pools = ("annotated_est", "nascent_est", "gdna_est", "n_intergenic_gdna")
    assert [nc[1]["library"][k] for k in pools] == [want[k] for k in pools]
    assert base[1]["library"]["annotated_est"] != want["annotated_est"], (
        "the plant moved nothing, so this test checks nothing"
    )
    assert (base[1]["qa_arm"], nc[1]["qa_arm"]) == ("base", "no_capture")
    row = PS.pair_records(_as_record(*base, "base"), _as_record(*nc, "nc"))
    assert row["varies"] == {"qa_arm": ["base", "no_capture"]}
    assert sum(row["moved"].values()) > 0


def test_a_no_capture_record_whose_run_read_capture_is_refused(toy, tmp_path, monkeypatch):
    """Break the arm: its wrapper installed, but the reading it replaces handed back untouched. The record must refuse
    the run rather than book the shipped reading as ``no_capture``."""
    qa = PS.sibling("quant_accuracy.py")
    _plant_capture(monkeypatch)
    real = qa.install_no_capture

    def inert():
        restore, fired = real()
        restore()
        return (lambda: None), {"n": 1}

    monkeypatch.setattr(qa, "install_no_capture", inert)
    with pytest.raises(RuntimeError, match="no_capture"):
        PS.record_bam(
            "nc",
            toy.bam_path,
            toy.index,
            PS.protocol_config([]),
            tmp_path / "work",
            {},
            CONDITION,
            qa_arm="no_capture",
        )


def test_a_mis_joined_read_name_on_real_records_is_refused(shipped):
    arrays, row = shipped
    shifted = {**arrays, "key": np.roll(arrays["key"], 1)}
    with pytest.raises(SystemExit, match="same_source"):
        PS.pair_records(_as_record(arrays, row, "a"), _as_record(shifted, row, "b"))
