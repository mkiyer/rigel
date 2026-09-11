"""An AMBIG region must not read its own RNA as gDNA, and must still see gDNA when there is gDNA.

An AMBIG exon region is one where a `+` and a `-` transcript overlap, so its own strand is undefined
and the strand channel tells it nothing: on the count alone, RNA it cannot attribute to either strand
looks exactly like genomic background. The answer has to come from the region's strand-clean flanking
boundaries, and this file is the end-to-end lock on that (sim → scan → calibrate on a controlled toy
genome: two overlapping opposite-strand transcripts forming one AMBIG region, single-strand flanks,
and separate multi-exon genes so the strand model trains). Two behaviours are held together, and
neither means much alone — at gDNA=0 with nascent RNA present the AMBIG region's gDNA fraction is
near zero, and with real gDNA and no nascent it reads substantial gDNA, which is what stops the first
arm from passing trivially by always answering zero.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.fl import build_fl_models
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _ambig_gdna_fraction(work_dir, *, gdna_abundance: int, nrna_abundance: float) -> float:
    """Build the toy AMBIG scenario, calibrate, and return the AMBIG region's contained gDNA fraction."""
    sc = Scenario("ambig_reg", genome_length=30000, seed=7, work_dir=work_dir)
    # Overlapping opposite-strand pair → an AMBIG region at 5000-6000 with single-strand flanks.
    sc.add_gene(
        "gA", "+", [{"t_id": "TA", "exons": [(1000, 1500), (4000, 6000)], "abundance": 100}]
    )
    sc.add_gene(
        "gB", "-", [{"t_id": "TB", "exons": [(5000, 7000), (10000, 10500)], "abundance": 100}]
    )
    # Standard multi-exon genes (both strands) so the strand model trains.
    sc.add_gene(
        "s1",
        "+",
        [
            {
                "t_id": "S1",
                "exons": [(12000, 12500), (13500, 14000), (15000, 15500)],
                "abundance": 120,
            }
        ],
    )
    sc.add_gene(
        "s2",
        "-",
        [
            {
                "t_id": "S2",
                "exons": [(17000, 17500), (18500, 19000), (20000, 20500)],
                "abundance": 120,
            }
        ],
    )
    gd = (
        GDNAConfig(
            abundance=gdna_abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000
        )
        if gdna_abundance > 0
        else None
    )
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(
            frag_mean=250,
            frag_std=50,
            frag_min=80,
            frag_max=600,
            read_length=100,
            strand_specificity=0.99,
            seed=7,
        ),
        gdna_config=gd,
        nrna_abundance=float(nrna_abundance),
    )
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_frame(idx.regions_df, idx.ref_name_to_id)
    # The same call production makes; see tests/calibration/_oracle.py.
    fl = build_fl_models(pl)
    # The sj axis and the boundary flags are BOTH required against the same index the payload was
    # scanned on: an axis addressing a different graph would place every splice on the wrong boundary, and
    # calibrate refuses rather than proceeding.
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )

    result = calibrate(
        pl,
        ra,
        sm,
        fl.gdna_pmf,
        fl.rna_pmf,
        cfg.calibration,
        sj=build_sj_geometry_arrays(idx),
        boundary_flags=build_boundary_flags_array(idx),
    )
    sc.cleanup()

    ambig = np.flatnonzero(np.asarray(ra.strand_class) == TS_AMBIG)
    assert ambig.size >= 1, "the overlapping pair did not form an AMBIG region"
    g = np.asarray(result.mass_gdna_region)[ambig]
    r = np.asarray(result.mass_rna_region)[ambig]
    return float(g.sum() / max(g.sum() + r.sum(), 1e-9))


def test_ambig_no_false_gdna_from_nascent(tmp_path):
    """This is the standing detector that message propagation keeps supplying an AMBIG slot's answer.

    An AMBIG slot has no own composition evidence: κ = ½ makes the strand λ-term identically 0, and
    the Schur complement on a both-strand region is exactly 0. With the messages muted the slot falls
    back to ψ's uninformative reference and reads roughly half gDNA against a truth of zero, so the
    bound below is a direct measurement of whether the neighbours are still being heard. Giving the
    slot its own θ-independent composition channel instead was priced and refused
    (`TRAPS: a-linear-likelihood-emits-a-sign`).
    """
    # gDNA=0 + nascent: the AMBIG region must NOT read the nascent/mature RNA as gDNA. Imputation
    # from the strand-cleaned flanking boundaries is what drives it to ~0.
    frac = _ambig_gdna_fraction(tmp_path / "none", gdna_abundance=0, nrna_abundance=30.0)
    assert frac < 0.08, (
        f"AMBIG gDNA fraction {frac:.3f} too high at gDNA=0+nascent (the fix should be ~0)"
    )


def test_ambig_reads_gdna_when_present(tmp_path):
    # Sanity (so the test above is not trivially satisfied by always reading 0): with real gDNA and no
    # nascent, the AMBIG region — imputed from its cleaned flanks — reads substantial gDNA.
    frac = _ambig_gdna_fraction(tmp_path / "gdna", gdna_abundance=200, nrna_abundance=0.0)
    assert frac > 0.2, f"AMBIG gDNA fraction {frac:.3f} too low with real gDNA present"
