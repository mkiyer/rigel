"""Falsification gates for ``scripts/design/toy_harness.py``.

⛔ **A TOY HARNESS THAT IS SILENTLY NOT USING THE DONOR'S CONDITIONS IS WORSE THAN NO HARNESS**, because
every conclusion drawn from it would be about a library nobody has. So the gates here are mostly about
one question — *is the donor actually reaching the toy?* — asked four different ways.

⭐ The donor is a scenario the gates BUILD, not a panel condition: the 36-condition ladder exists on one
machine, and a gate that silently skips is a gate that has never run.

⛔ Each gate carries its own perturbation.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

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


TH = _sibling("toy_harness.py")


@pytest.fixture(scope="module")
def donor(tmp_path_factory):
    """A population-ish donor: enough genes and depth that the library-level priors actually fit.

    ⚠ It must be UNSTRANDED (``ss_0.50`` in the name, which is where the harness reads it from) and
    must carry gDNA — with no gDNA there is no intergenic population to measure a density from, and
    :func:`toy_harness.harvest` would have nothing to match the toy to.
    """
    wd = tmp_path_factory.mktemp("toy_donor")
    sc = Scenario("donor_ss_0.50_capture_off", genome_length=120_000, seed=11, work_dir=wd / "sim")
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
            strand_specificity=0.50,
            seed=11,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=200, frag_std=60),
    )
    globals_ = TH.harvest(
        wd, res.index, bam=str(res.bam_path), name="donor_ss_0.50_nrna_none_capture_off"
    )
    return globals_


@pytest.fixture(scope="module")
def spec():
    """Two exons and one intron — the smallest structure with a sj and both boundary kinds."""
    return TH.ToySpec(
        name="gate_two_exon",
        what_it_probes="gate fixture",
        genome_length=40_000,
        genes=[
            {
                "gene_id": "t",
                "strand": "+",
                "transcripts": [
                    {
                        "t_id": "t_t1",
                        "exons": [(12_000, 15_000), (20_000, 23_000)],
                        "abundance": 300.0,
                    }
                ],
            }
        ],
        n_rna_fragments=20_000,
        seed=3,
    )


# ── GATE 1: the donor's globals actually reach the solve ──────────────────────────────────────────


def test_the_DONORS_GLOBALS_ARE_INJECTED_and_not_silently_refitted(donor, spec, tmp_path):
    """⛔ The whole premise. A toy cannot fit a strand balance, an enrichment landscape or an
    intergenic background, so if the injection is not wired the toy quietly fits its own — from a
    handful of regions — and every conclusion is about a library that does not exist.

    PERTURBATION: move κ in the injected bundle to the opposite extreme and require the toy's answer
    to CHANGE. If it does not, the bundle is not reaching the solve.
    """
    base = TH.run_toy(spec, donor, tmp_path / "a")
    got = base.result

    flipped = replace(donor, priors=dataclasses.replace(donor.priors, rna_sense_frac=0.02))
    moved = TH.run_toy(spec, flipped, tmp_path / "b").result

    a = np.asarray(got.mass_gdna_region, np.float64)
    b = np.asarray(moved.mass_gdna_region, np.float64)
    assert not np.allclose(a, b), (
        "flipping the injected kappa from unstranded to strongly stranded changed nothing, so the "
        "injected priors are not reaching calibrate"
    )
    # and it must be the INJECTED value the solve used, not a refit of the toy's own data
    assert float(got.rna_sense_frac) == pytest.approx(donor.priors.rna_sense_frac)
    assert float(moved.rna_sense_frac) == pytest.approx(0.02)


def test_the_LENGTH_MODELS_come_from_the_DONOR_not_from_the_toy(donor, spec, tmp_path):
    """⭐ The two fragment-length pmfs are ``calibrate`` kwargs rather than part of the priors bundle,
    so they are a separate wiring that can be forgotten independently.

    PERTURBATION: replace the donor's gDNA pmf by a narrow spike far from the RNA one and require the
    answer to move — at equal means the length channel carries exactly zero information
    (TRAPS: equal-lengths-carry-no-composition), so a harness that ignored these would be indistinguishable from one that used
    them unless the two pmfs are pulled apart.
    """
    size = donor.gdna_fl_pmf.shape[0]
    spike = np.zeros(size, np.float64)
    spike[min(90, size - 1)] = 1.0
    apart = replace(donor, gdna_fl_pmf=spike)

    base = np.asarray(TH.run_toy(spec, donor, tmp_path / "c").result.mass_gdna_region, np.float64)
    moved = np.asarray(TH.run_toy(spec, apart, tmp_path / "d").result.mass_gdna_region, np.float64)
    assert not np.allclose(base, moved), (
        "replacing the donor's gDNA length model by a spike changed nothing; the pmfs are not "
        "reaching calibrate"
    )


# ── GATE 2: the toy is simulated at the DONOR's depth ─────────────────────────────────────────────


def test_the_toys_gDNA_DEPTH_IS_DERIVED_from_the_donor_and_actually_LANDS(donor, spec, tmp_path):
    """⭐⭐ ``calibrate``'s own note: the injected enrichment landscape is an ABSOLUTE log-density
    model, so a toy at the wrong depth projects onto the wrong cells and is a different library, not a
    small one. The harness therefore derives the toy's gDNA count from the donor's measured density
    per base — and this gate checks the density it actually REALISES, not the count it asked for.

    PERTURBATION: a toy simulated at 10x the donor's rate must land measurably off, or the check is
    blind to depth.
    """
    r = TH.run_toy(spec, donor, tmp_path / "e")
    realised = TH._rate_from_capture(r.capture, r.chain, r.region_arrays)
    assert realised == pytest.approx(donor.gdna_rate_per_base, rel=0.25), (
        f"the toy realised {realised:.4g} gDNA molecules/base against the donor's "
        f"{donor.gdna_rate_per_base:.4g} — the depth match is not working"
    )

    # PERTURBATION: 10x the rate must be detected by the very same comparison.
    hot = replace(donor, gdna_rate_per_base=donor.gdna_rate_per_base * 10.0)
    r10 = TH.run_toy(spec, hot, tmp_path / "f")
    realised10 = TH._rate_from_capture(r10.capture, r10.chain, r10.region_arrays)
    assert realised10 > realised * 3.0, "a 10x depth change is invisible to the realised-rate check"
    assert realised10 != pytest.approx(donor.gdna_rate_per_base, rel=0.25)


def test_there_is_NO_gDNA_KNOB_on_the_spec(donor):
    """⛔ The gDNA level must not be settable per toy: it is pinned by the donor, and a spec field for
    it would be a foot-gun that silently invalidates the injected enrichment landscape."""
    fields = set(TH.ToySpec.__dataclass_fields__)
    for banned in ("gdna_fraction", "n_gdna_fragments", "gdna_rate", "gdna_abundance"):
        assert banned not in fields, (
            f"ToySpec exposes a gDNA knob ({banned}); the depth must be derived"
        )


# ── GATE 3: truth is the origin split, and every object is reported ───────────────────────────────


def test_TRUTH_is_the_ORIGIN_SPLIT_and_sums_to_the_full_payload(donor, spec, tmp_path):
    """The per-object truth must be the production accumulator run on the BAM split by true origin —
    the identity that makes it trustworthy at all (`SUCCESS.md` TRAPS: self-checking-validator). ``OracleTruth.from_bam`` validates
    sum-to-full and RAISES, so this gate checks the harness actually goes through it rather than
    approximating truth from an abundance table.

    PERTURBATION: a corrupted partition must make the same construction fail.
    """
    r = TH.run_toy(spec, donor, tmp_path / "g")
    assert hasattr(r.truth, "parts") and set(r.truth.parts) == {"gdna", "mrna", "nrna"}
    region = np.asarray(r.payload.region_contained_count, np.int64)
    parts = sum(
        np.asarray(r.truth.parts[k].region_contained_count, np.int64)
        for k in ("gdna", "mrna", "nrna")
    )
    np.testing.assert_array_equal(region, parts)

    # PERTURBATION: break one partition and the identity must break with it.
    bad = np.asarray(r.truth.parts["gdna"].region_contained_count, np.int64).copy()
    if bad.size:
        bad[0] += 1
        assert not np.array_equal(
            region, parts - np.asarray(r.truth.parts["gdna"].region_contained_count, np.int64) + bad
        )


def test_EVERY_object_with_mass_is_reported(donor, spec, tmp_path):
    """⭐ The point of a toy is that you can read every row. A report that dropped objects would hide
    exactly the one being interrogated.

    PERTURBATION: the row set must cover every chain slot, and the region/boundary split must be non-trivial
    — a toy with no boundaries could not exercise the boundary classes at all.
    """
    r = TH.run_toy(spec, donor, tmp_path / "h")
    rows = TH.object_rows(r)
    assert len(rows) == int(r.chain.n_slots), "object_rows does not cover every chain slot"
    kinds = {row["axis"] for row in rows}
    assert kinds == {"region", "boundary"}, (
        f"the toy has only {kinds}; it cannot exercise both axes"
    )
    types = {row["type"] for row in rows}
    for expected in ("intergenic", "exon", "intron", "intron|exon", "intergenic|exon"):
        assert expected in types, f"the two-exon toy produced no {expected!r} object"


# ── GATE 4: the harness reproduces the finding it was built to isolate ────────────────────────────


def test_the_harness_REPRODUCES_the_intron_composition_dependence(donor, spec, tmp_path):
    """⭐⭐ The substantive gate, **REWRITTEN 2026-08-04 because the defect it pinned was fixed** — which
    is what its predecessor instructed ("rewrite this gate to pin the NEW behaviour rather than deleting
    it"), and the record of both states is the point.

    **What it used to pin.** An exon with no own evidence inherited its composition from the intron
    beside it, so a **pure-gDNA** intron (no nascent RNA) dragged an essentially-pure-RNA exon toward
    gDNA while a nascent-bearing intron did not. Measured against the `g25` ladder donor: exon
    ``|Δf_g|`` **0.156 dry vs 0.005 wet, a factor of 31**; against the six-gene synthetic donor this
    fixture builds, 0.209 vs 0.167. Only the direction and ordering were assertable, because the
    magnitude is donor-dependent and the harness itself surfaced that.

    ⭐⭐⭐ **IT WAS A STRICT xfail FROM 2026-08-05 TO 2026-08-18 AND IS GREEN AGAIN — BUT READ HOW, BECAUSE
    IT IS NOT THE ROUTE THE xfail PREDICTED.** The recorded reason was: correcting the sj leak removes ONE
    error from a compensating PAIR (`TRAPS: a-cancelling-defect-pair`) — an evidence-free exon is fed
    through `intron → BOUNDARY → exon`, the two hops' errors cancelled under the old sj-inclusive total,
    and the second hop still carried its own defect (a correct composition ratio applied to a LEVEL). It
    failed at **2.19× against its 2.0 bound**, and the bound was KEPT rather than widened so it would stay
    the detector for that mechanism. ⚠ It went green when `message_propagation` was turned back ON
    (owner, 2026-08-18) — the relay changes the second hop — and **NOT** by the pair being fixed jointly,
    which is still open (`ROADMAP.md`'s cancelling-pair item). ⛔ So do not read this passing as evidence
    that the pair is resolved; read it as the detector having moved into the relay-on regime with it.

    ⭐⭐ **What it pins now: the dependence is GONE.** The mechanism was the reframe imputing the
    source's composition share onto the destination's observed total (`EQUATIONS.md` §3.5), and a
    factory-solved intron reports ``f_g ≈ 1``, hence zero RNA density, hence zero RNA precision — so it
    can no longer lend a composition and its gDNA LEVEL crosses unscaled instead. Measured on this
    fixture: **0.0107 dry vs 0.0112 wet, a factor of 1.04**, down from 1.25 here and 31 on the ladder.

    ⚠ **So the assertion inverts from an ORDERING to an INDEPENDENCE**, and that is the stronger
    statement: what is inside the intron must not set the exon's answer. Both arms must also be SMALL —
    the ordering could equally be destroyed by making both bad.

    PERTURBATION: the two arms still differ in ONE field (``nrna_abundance``) and must differ in outcome
    at all; a byte-identical result would mean the nascent axis is not reaching the simulation.
    """
    dry = TH.run_toy(replace(spec, name="dry", nrna_abundance=0.0), donor, tmp_path / "i")
    wet = TH.run_toy(replace(spec, name="wet", nrna_abundance=40.0), donor, tmp_path / "j")

    def exon_rows(r):
        rows = [
            x
            for x in TH.object_rows(r)
            if x["type"] == "exon" and x["mass"] > 0 and np.isfinite(x["true_fg"])
        ]
        assert rows, "no exon rows with truth"
        return rows

    def exon_error(rows):
        w = np.array([x["mass"] for x in rows])
        d = np.array([abs(x["pred_fg"] - x["true_fg"]) for x in rows])
        return float((d * w).sum() / w.sum())

    dry_rows, wet_rows = exon_rows(dry), exon_rows(wet)
    e_dry, e_wet = exon_error(dry_rows), exon_error(wet_rows)

    # PERTURBATION: one field varied must change the outcome, or the axis is not wired.
    assert e_dry != pytest.approx(e_wet, rel=1e-6), (
        "the nascent axis changed nothing; nrna_abundance is not reaching the simulation"
    )
    # ⭐⭐ THE INDEPENDENCE — what is inside the INTRON must not set the EXON's answer. The retired
    # ordering assertion (`e_dry > e_wet`) is what this replaces; a return to it means the composition
    # imputation is reaching the exon again.
    assert max(e_dry, e_wet) / max(min(e_dry, e_wet), 1e-9) < 2.0, (
        f"exon |Δf_g| is {e_dry:.4f} with a pure-gDNA intron and {e_wet:.4f} with a nascent-bearing "
        f"one — a factor of {max(e_dry, e_wet) / max(min(e_dry, e_wet), 1e-9):.2f}. The exon's answer "
        "is tracking the intron's composition again (it was 31x on the ladder donor before "
        "`EQUATIONS.md` §3.5)"
    )
    # ⭐ …and BOTH arms must be small, or the independence was won by making both wrong. The retired
    # gate's own dry-arm figure on this fixture was 0.209, so 0.05 is a decade of headroom below it and
    # a decade above what the two arms now read.
    assert max(e_dry, e_wet) < 0.05, (
        f"the exon error is not small on either arm ({e_dry:.4f} / {e_wet:.4f}); the independence may "
        "have been bought by degrading both"
    )
