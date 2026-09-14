"""The encompassing locus — a single-exon TB− over a two-exon TA+ — at four abundance regimes (owner,
2026-09-13). The RNA level lanes must pass TB−'s level, measured in its single-strand exons, through every
one of TA+'s boundaries into the both-stranded exons, and TA+'s junction flux must be a + source there;
the both-stranded exons must then solve to their truth, and with both genes silent everything reads gDNA.

The lanes' rule this gates: a strand's level crosses a face unless the boundary carries that strand's own
bits; here every boundary carries only TA+'s bits and every region admits −. Found broken three ways by
the audit that made this file (the lanes worklist, 2026-09-13): the layer was off without an intron factory,
the RNA lanes died with the gDNA lane, and a strand with no single-strand exon had no coordinate and hence
no flux source.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

_HERE = Path(__file__).resolve()
_DESIGN = _HERE.parents[2] / "scripts" / "design"


def _load(name: str):
    key = f"_encompass_{name}"
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, _DESIGN / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


TH = _load("toy_harness")

#: TA+ against TB−, in the simulator's abundance units: TA+ ≫ TB−, TA+ ≪ TB−, TA+ ≈ TB−, both silent
REGIMES = {
    "ta_high_tb_low": (1000.0, 30.0),
    "ta_low_tb_high": (30.0, 1000.0),
    "ta_equal_tb": (500.0, 500.0),
    "both_silent": (0.0, 0.0),
}


@pytest.fixture(scope="module")
def donor(tmp_path_factory):
    """A stranded (``ss_0.99``) donor with gDNA, so the harness has a strand model and a gDNA density to
    match the toy to."""
    wd = tmp_path_factory.mktemp("encompass_donor")
    sc = Scenario("donor_ss_0.99_capture_off", genome_length=120_000, seed=11, work_dir=wd / "sim")
    for i in range(6):
        base = 10_000 + i * 18_000
        sc.add_gene(
            f"g{i}",
            "+" if i % 2 == 0 else "-",
            [
                {
                    "t_id": f"g{i}_t1",
                    "exons": [(base, base + 2_500), (base + 6_000, base + 8_500)],
                    "abundance": 40.0 * (i + 1),
                }
            ],
        )
    res = sc.build_oracle(
        n_rna_fragments=60_000,
        gdna_fraction=0.6,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=60,
            frag_min=80,
            frag_max=400,
            read_length=100,
            strand_specificity=0.99,
            seed=11,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=200, frag_std=60),
    )
    return TH.harvest(
        wd, res.index, bam=str(res.bam_path), name="donor_ss_0.99_nrna_none_capture_off"
    )


def _spec(regime: str) -> TH.ToySpec:
    ta, tb = REGIMES[regime]
    silent = ta == 0.0 and tb == 0.0
    return TH.ToySpec(
        name=f"encompass_{regime}",
        what_it_probes="a single-exon TB− encompassing a two-exon TA+",
        genome_length=40_000,
        genes=[
            TH._gene("gA", "+", [(11_000, 12_000), (19_000, 20_000)], ta if not silent else 1.0),
            TH._gene("gB", "-", [(10_000, 30_000)], tb if not silent else 1.0),
        ],
        # both silent: a handful of RNA fragments against the donor's gDNA depth, so the locus is gDNA
        n_rna_fragments=200 if silent else 20_000,
        seed=7,
    )


@pytest.fixture(scope="module")
def runs(donor, tmp_path_factory):
    wd = tmp_path_factory.mktemp("encompass")
    return {regime: TH.run_toy(_spec(regime), donor, wd / regime) for regime in REGIMES}


def _slots(r):
    """The toy's slots by role: TA+'s two exons (exon on both strands), the region between them (TA+'s
    intron, TB−'s exon), TB−'s single-strand flanks, and the AMBIG boundaries."""
    rows = TH.object_rows(r)
    by = {}
    for row in rows:
        if row["axis"] == "region":
            start = int(row["where"].split("–")[0].replace(",", ""))
            by.setdefault(start, row)
    exons = [by[11_000], by[19_000]]
    between = by[12_000]
    flanks = [by[10_000], by[20_000]]
    ambig_boundaries = [
        row for row in rows if row["axis"] == "boundary" and row["where"] in ("@12,000", "@19,000")
    ]
    return exons, between, flanks, ambig_boundaries, rows


def test_tb_minus_level_reaches_every_both_stranded_slot(runs):
    """With TB− expressed its level, measured in its single-strand flanks, is delivered into the cube of
    every both-stranded slot with counts: TA+'s exons, the region between them and the AMBIG boundaries."""
    for regime in ("ta_high_tb_low", "ta_low_tb_high", "ta_equal_tb"):
        r = runs[regime]
        exons, between, _flanks, bnds, _rows = _slots(r)
        cube = r.capture.cube_rows or {}
        for row in [*exons, between, *bnds]:
            assert row["n"] > 0, (regime, row["where"])
            got = cube.get(int(row["slot"]))
            assert got is not None and got.profile_neg is not None, (
                f"{regime}: no − level delivered at {row['where']} (slot {row['slot']})"
            )


def test_ta_plus_junction_flux_is_a_source_at_its_exons(runs):
    """TA+ has no single-strand exon anywhere, so on a per-strand coordinate its junction's certified flux
    was silently unused; the + level must be delivered at both of TA+'s exons whenever TA+ is expressed."""
    for regime in ("ta_high_tb_low", "ta_low_tb_high", "ta_equal_tb"):
        r = runs[regime]
        exons, _between, _flanks, _bnds, _rows = _slots(r)
        cube = r.capture.cube_rows or {}
        for row in exons:
            got = cube.get(int(row["slot"]))
            assert got is not None and got.profile_pos is not None, (
                f"{regime}: no + level delivered at TA+'s exon {row['where']}"
            )


@pytest.mark.xfail(
    strict=True,
    reason="ISSUES: capture-on-strand-pure-ambig-undercall — the region between TA+'s exons (TB−'s exon, TA+'s "
    "intron) is strand-pure and mostly gDNA when TB− is low, and the tilt continuum's median sits below the "
    "strand cap; closes with the witnessed atom (the lanes worklist's L5)",
)
def test_the_both_stranded_exons_solve_to_their_truth(runs):
    """The exon-on-exon slots — mostly +, mostly −, or balanced — read their gDNA fraction within 0.05 of
    the truth once the levels reach them (0.29–0.44 against 0.00 when they did not)."""
    for regime in ("ta_high_tb_low", "ta_low_tb_high", "ta_equal_tb"):
        r = runs[regime]
        exons, between, flanks, _bnds, _rows = _slots(r)
        for row in [*exons, between, *flanks]:
            assert abs(row["pred_fg"] - row["true_fg"]) < 0.05, (
                f"{regime}: {row['where']} reads f_g {row['pred_fg']:.3f} against {row['true_fg']:.3f}"
            )


def test_with_both_genes_silent_the_locus_reads_gdna(runs):
    """No RNA to speak of: every slot with counts is gDNA and must read so."""
    r = runs["both_silent"]
    _exons, _between, _flanks, _bnds, rows = _slots(r)
    counted = [row for row in rows if row["n"] >= 20 and np.isfinite(row["true_fg"])]
    assert counted, "the silent regime has no counted slot: the gate is vacuous"
    assert min(row["true_fg"] for row in counted) > 0.8, "the premise: the locus is gDNA"
    for row in counted:
        # 0.10: a gene-edge boundary here holds a few dozen fragments, and its solve is that shallow
        assert abs(row["pred_fg"] - row["true_fg"]) < 0.10, (
            f"{row['where']} reads f_g {row['pred_fg']:.3f} against {row['true_fg']:.3f} on a gDNA locus"
        )
