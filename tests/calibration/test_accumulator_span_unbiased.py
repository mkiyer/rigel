"""Behavioral guard: the accumulator deposits the molecule's span, not read blocks.

Phase A (test-first) of the accumulator span redesign
(``docs/calibration/accumulator_fragment_span_redesign.md``). On a uniform-density
(no-capture) genome the boundary-crossing density estimator (one-side gDNA flux /
fl_mean) must agree with the exact contained estimator (contained count /
region_eff_len). Today it over-counts by ~1.2–1.56× (paired-end mate-gap slice
crediting), so this is ``xfail(strict=True)`` until the span fix lands — at which
point it flips to a hard failure to prompt removing the marker.

Kept small + deterministic (fixed seed, ~60k fragments) so it runs fast.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import _check_region_payload_alignment, scan_and_buffer
from rigel.splice import SpliceType


def _crossing_vs_contained_ratio(bam_path, index) -> float:
    """Return crossing-ρ / contained-ρ over count-observable nodes (≈1.0 if unbiased)."""
    _s, sm, fla, _b, pl = scan_and_buffer(str(bam_path), index, BamScanConfig(sj_strand_tag="auto"))
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    _check_region_payload_alignment(ra, pl)
    flm = build_fl_models(
        global_counts=fla.global_model.counts,
        rna_counts=fla.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl),
        max_size=fla.max_size,
    )
    gpmf = flm.gdna_pmf
    reg_eff = region_eff_length(ra.region_size_bp, gpmf)
    fl_mean = boundary_eff_length(gpmf)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    rids = np.asarray(ra.ref_id)
    reg_obs, bnd_obs = count_observable_masks(np.asarray(ra.signature), rids)

    c = sub.contained
    cont_cnt = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.float64)
    obs_reg = reg_obs & (reg_eff > 1.0)
    rho_contained = float(cont_cnt[obs_reg].sum() / reg_eff[obs_reg].sum())

    flux = []
    for r in np.where(bnd_obs)[0]:
        rs = r + 1
        if rs < ra.n_regions and rids[r] == rids[rs]:
            flux.append((float(sub.right.n_unspliced[r]) + float(sub.left.n_unspliced[rs])) / 2.0)
    rho_crossing = float(np.mean(flux)) / fl_mean
    return rho_crossing / rho_contained


def test_crossing_density_unbiased(tmp_path):
    from rigel.sim import Scenario
    from rigel.sim.reads import GDNAConfig, ReadSimConfig

    sc = Scenario("span_unbiased", genome_length=40000, seed=29, work_dir=str(tmp_path / "sim"))
    sc.add_gene(
        "G", "+",
        [{"t_id": "G.1", "exons": [(12000, 13000), (16000, 17000), (20000, 21000)], "abundance": 40}],
    )
    res = sc.build_oracle(
        n_fragments=60000,
        sim_config=ReadSimConfig(strand_specificity=0.99, frag_mean=250, frag_std=50, seed=29),
        gdna_config=GDNAConfig(abundance=4000.0, frag_mean=350, frag_std=100),
    )
    ratio = _crossing_vs_contained_ratio(res.bam_path, res.index)
    # Unbiased target: the crossing estimator recovers the same uniform density.
    assert 0.9 <= ratio <= 1.1, f"crossing/contained density ratio = {ratio:.3f} (expected ≈1.0)"


def test_implicit_splice_routes_to_spliced_channel(tmp_path):
    """An implicitly-spliced fragment (annotated intron inside the mate gap, no
    CIGAR-N) must deposit on the SPLICED channels (ch2/3) and CUT the intron —
    not fill across it on the unspliced channels. Pre-Phase-C it was mis-channeled
    as unspliced and filled [min,max], depositing unspliced mass through the intron.
    """
    import numpy as np

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.signature import coarse_type_array
    from rigel.sim import Scenario
    from rigel.sim.reads import ReadSimConfig

    sc = Scenario("implicit_chan", genome_length=6000, seed=7, work_dir=str(tmp_path / "sim"))
    # Short intron (50 bp ≪ insert) → fragments spanning both exons place the
    # junction in the unsequenced mate gap = implicit splice (no CIGAR-N).
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(1000, 1150), (1200, 1500)], "abundance": 80}])
    res = sc.build_oracle(
        n_fragments=40000,
        sim_config=ReadSimConfig(
            frag_mean=250, frag_std=40, frag_min=180, frag_max=400,
            read_length=100, strand_specificity=1.0, seed=7,
        ),
    )
    _s, sm, _fl, _b, pl = scan_and_buffer(
        str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto")
    )
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    ctype = coarse_type_array(np.asarray(ra.signature))
    rid = np.asarray(ra.ref_id)
    gene = rid == res.index.ref_name_to_id["implicit_chan"]

    # (1) Implicit splices route to the SPLICED channel: the exon-flanking
    #     boundaries carry spliced crossing flux.
    spliced_flux = int((sub.left.n_spliced[gene] + sub.right.n_spliced[gene]).sum())
    assert spliced_flux > 1000, f"expected substantial spliced crossing flux, got {spliced_flux}"

    # (2) The intron is CUT, not filled: the intron region carries NO unspliced
    #     mass (contained or crossing) — the implicit molecules skip it.
    intron = gene & (ctype == 1)
    assert intron.any()
    intron_unspliced = float(
        sub.contained.n_unspliced[intron].sum()
        + sub.left.n_unspliced[intron].sum()
        + sub.right.n_unspliced[intron].sum()
    )
    assert intron_unspliced == 0.0, f"intron carries unspliced mass {intron_unspliced} (not cut)"


@pytest.mark.skip(reason="Phase D: artifact-splice hold-out not yet implemented")
def test_artifact_splice_held_out_and_mass_conserved():
    """``SPLICE_ARTIFACT`` fragments are not deposited (no contained/boundary mass,
    no FL pool), and total deposited mass == n_deposited (artifacts excluded)."""
    raise NotImplementedError
