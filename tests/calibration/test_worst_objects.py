"""Falsification gates for ``scripts/design/worst_objects.py`` — step 3 of the debug loop.

The instrument decides **which objects a human then spends hours staring at**. Every way of getting
that ranking or its context subtly wrong is a way of sending the reader after the wrong mechanism,
and none of them would look like a failure — the table would still be full of plausible rows. These
gates are those ways.

⛔ Each gate carries its own perturbation, in the same test.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from rigel.config import CalibrationConfig, PipelineConfig
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        path = Path(__file__).resolve().parents[2] / "scripts" / "design" / name
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


WO = _sibling("worst_objects.py")
P0 = _sibling("pass0_vs_oracle.py")


@pytest.fixture(scope="module")
def measured(tmp_path_factory):
    wd = tmp_path_factory.mktemp("wo_sim")
    sc = Scenario("wo", genome_length=9000, seed=23, work_dir=wd / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 80},
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 20},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 30}])
    sc.add_gene(
        "g3",
        "+",
        [
            {"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 25},
            {"t_id": "t3b", "exons": [(6600, 7100), (7400, 7440), (7800, 8300)], "abundance": 10},
        ],
    )
    res = sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            read_length=90,
            strand_specificity=0.99,
            seed=23,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )
    m = P0.measure_condition(
        bam=str(res.bam_path),
        index=res.index,
        pipeline_config=PipelineConfig(),
        calibration_config=CalibrationConfig(),
        work_dir=tmp_path_factory.mktemp("wo_split"),
        tag="wo",
    )
    return m, res.index


def _dissect(measured, axis="region", top=25):
    m, index = measured
    return m, WO.dissect(m, axis=axis, arm="pass0", top=top, index=index)


# ── GATE 1: the ranking is by error MASS ──────────────────────────────────────────────────────────


def test_the_ranking_is_by_error_MASS_and_a_RATE_ranking_would_differ(measured):
    """⭐ The whole instrument turns on this. A 1 bp region with two fragments can carry ``|Δf_g| = 1``
    and be worth two fragments of error; ranking by rate puts it above an exon carrying thousands.

    PERTURBATION: rank the same objects by ``|Δf_g|`` instead and require a DIFFERENT leader. If the
    two orderings agreed on this scenario the gate would be vacuous and would pass with the defect.
    """
    m, d = _dissect(measured)
    errs = [abs(r["err"]) for r in d["rows"]]
    assert errs == sorted(errs, reverse=True), "rows are not ordered by |error mass|"

    live, err = d["live"], d["err"]
    total = d["total"]
    rate = np.zeros_like(err)
    np.divide(np.abs(err), total, out=rate, where=live & (total > 0))
    by_rate = int(np.argmax(np.where(live, rate, -1.0)))
    by_mass = d["rows"][0]["obj"]
    assert by_rate != by_mass, (
        "the mass-ranked and rate-ranked leaders coincide on this scenario, so this gate cannot "
        "detect a rate ranking — the fixture needs an object that is small but badly wrong"
    )
    assert abs(err[by_mass]) > abs(err[by_rate]), "the mass ranking did not pick the larger error"


# ── GATE 2: objects with no mass are never ranked ─────────────────────────────────────────────────


def test_an_object_with_NO_MASS_never_appears_in_the_table(measured):
    """⛔ Most of any real index carries no fragments. A ranking that admitted them would fill with
    ``0/0`` rows whose ``pred_fg`` is NaN and whose error is exactly zero — harmless-looking padding
    that pushes real objects off the end of the list.

    PERTURBATION: the scenario must actually CONTAIN empty objects, or the gate proves nothing.
    """
    m, d = _dissect(measured)
    assert (~d["live"]).sum() > 0, "no empty objects in the fixture; the gate would be vacuous"
    for r in d["rows"]:
        assert r["true_g"] + r["true_r"] > 0
        assert np.isfinite(r["pred_fg"]) and np.isfinite(r["true_fg"])
    assert d["n_live"] == int(d["live"].sum())


# ── GATE 3: the neighbours are adjacent CHAIN SLOTS ───────────────────────────────────────────────


def test_the_neighbour_errors_are_ADJACENT_CHAIN_SLOTS_not_adjacent_objects(measured):
    """⭐ On an object with no own evidence the neighbour columns are the whole explanation, so they
    must be the objects that actually message it. The chain is ``N E N E … N``, so a REGION's
    neighbours are two BOUNDARIES — **not** the regions at ``obj ± 1``, which are two slots away and do not
    message it directly.

    PERTURBATION: compute the same-axis ``obj ± 1`` errors and require them to DIFFER, so the gate
    would fail if the instrument used that (much more natural-looking) indexing.
    """
    m, d = _dissect(measured)
    chain = m.debug_pass0["chain"]
    n_regions, n_boundaries = int(m.payload.n_regions), int(m.payload.n_boundaries)
    region_slot, boundary_slot = WO._slot_lookup(chain, n_regions, n_boundaries)
    kind = np.asarray(chain.kind)
    obj_idx = np.asarray(chain.obj_idx, dtype=np.int64)

    boundary_err = np.asarray(m.arms["pass0"].mass_gdna_boundary, np.float64) - np.asarray(
        m.truth.mass_gdna_boundary, np.float64
    )
    region_err = np.asarray(m.arms["pass0"].mass_gdna_region, np.float64) - np.asarray(
        m.truth.mass_gdna_region, np.float64
    )

    checked = differed = 0
    for r in d["rows"]:
        s = int(region_slot[r["obj"]])
        for side, nb in zip((s - 1, s + 1), r["nb_err"]):
            if np.isnan(nb) or not 0 <= side < kind.shape[0]:
                continue  # a reference boundary has no neighbour on that side
            assert kind[side] == WO.BOUNDARY, "a REGION's chain neighbour must be an BOUNDARY slot"
            np.testing.assert_allclose(nb, boundary_err[obj_idx[side]], rtol=1e-9, atol=1e-6)
            checked += 1
        # the same-axis alternative must not coincide, or the gate cannot see the difference
        for other in (r["obj"] - 1, r["obj"] + 1):
            if 0 <= other < n_regions and not np.isclose(region_err[other], r["nb_err"][0]):
                differed += 1
    assert checked > 0, "no neighbours were checked; the gate would be vacuous"
    assert differed > 0, "same-axis indexing gives the same answer here; the gate cannot detect it"


# ── GATE 4: the concentration curve is monotone and complete ──────────────────────────────────────


def test_the_concentration_curve_is_MONOTONE_and_reaches_the_whole_error(measured):
    """The curve is the first thing a reader acts on — concentrated sends them hunting a mechanism,
    diffuse sends them after a systematic bias. It must therefore be a real cumulative curve.

    PERTURBATION: a curve computed on UNSORTED error would not be monotone; assert monotonicity and
    that the top-k share of a deliberately concentrated input is ~1.
    """
    m, d = _dissect(measured)
    shares = [s for _, s, _ in d["concentration"]]
    assert shares == sorted(shares), "the cumulative share decreased"
    assert all(0.0 <= s <= 1.0 + 1e-12 for s in shares)

    spike = np.zeros(5000)
    spike[7] = 1000.0
    spike[11] = -900.0
    curve = WO.concentration(spike)
    assert curve[0][1] > 0.999, "a two-object spike must put ~all the error in the top 10"
    flat = np.full(5000, 3.0)
    flat_curve = WO.concentration(flat)
    assert flat_curve[0][1] == pytest.approx(10 / 5000, rel=1e-6), "a flat input must be diffuse"


# ── GATE 5: the profile reports the BACKGROUND rate, not just the top's ───────────────────────────


def test_the_profile_reports_BOTH_the_top_share_and_the_background_share(measured):
    """⭐ "80 % of the worst objects are exons" means nothing if exons are 80 % of everything. The
    profile must carry both numbers so the ratio is visible.

    PERTURBATION: each set of shares must sum to 1 over the class partition — a profile that dropped
    a class, or double-counted one, breaks that identity in a way a single share could hide.
    """
    m, d = _dissect(measured)
    for masks, names in ((d["solver"], P0.SOLVER_CLASSES), (d["info"], P0.INFO_CLASSES)):
        rows = WO._profile(d, masks, names)
        assert {r[0] for r in rows} == set(names)
        assert sum(r[1] for r in rows) == pytest.approx(1.0, abs=1e-9)
        assert sum(r[2] for r in rows) == pytest.approx(1.0, abs=1e-9)
        assert all(len(r) == 3 for r in rows), "the background share is not reported"


# ── GATE 6: the classes are the ones pass0_vs_oracle scored ───────────────────────────────────────


def test_the_classes_are_the_SCORED_ones_not_a_second_computation(measured):
    """⛔ Two definitions of one class is how they drift. The dissection must consume the masks the
    measurement already scored, not recompute its own from the capture.

    PERTURBATION: recompute them here from the raw capture and require byte-identity — if the
    instrument ever forks its own copy, this fails.
    """
    m, d = _dissect(measured)
    fresh = P0.solver_class_masks(
        m.debug_pass0["capture"],
        m.debug_pass0["chain"],
        int(m.payload.n_regions),
        int(m.payload.n_boundaries),
    )["region"]
    for name in P0.SOLVER_CLASSES:
        np.testing.assert_array_equal(d["solver"][name], fresh[name])
        np.testing.assert_array_equal(d["solver"][name], m.solver_masks["region"][name])
