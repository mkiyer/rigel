"""Phase A: mature_density — per-strand contained-unspliced mature imputation from spliced crossings.

End-to-end on a controlled toy genome (sim → scan → mature_density): a clean multi-exon ``+`` gene (no
gDNA, no nascent) so the contained-unspliced count IS the mature. Locks:

  1. a substantial internal exon's predicted ``M_pos`` ≈ its contained-unspliced count (the imputation
     recovers the mature from the bounding splice junctions);
  2. intron / intergenic regions get ``M = 0`` (mature is predicted only where exons carry it — never
     invented where there is no junction-anchored exon);
  3. the ``−`` strand field is ~0 for a ``+``-only gene (per-strand separation).

See docs/calibration/phaseA_mature_imputation_plan.md (Issue D validated in
scripts/debug/phaseA_issueD_mature_density.py).
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.mature_density import mature_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType

_EXON = BIT_EXON_POS | BIT_EXON_NEG


def _run(work_dir):
    sc = Scenario("mature_density", genome_length=40000, seed=11, work_dir=work_dir)
    # One + gene with a LONG internal exon (the imputation target) bounded by two splice junctions.
    sc.add_gene("gP", "+", [{"t_id": "TP", "exons": [
        (2000, 2500), (4000, 4150), (6000, 7500), (9000, 9500), (11000, 11400), (13000, 13500),
    ], "abundance": 300}])
    # A − gene so the strand model + balance train (needs spliced reads on both strands).
    sc.add_gene("gN", "-", [{"t_id": "TN", "exons": [(20000, 20500), (22000, 22500), (24000, 24500)],
                             "abundance": 120}])
    res = sc.build_oracle(
        n_fragments=20000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=None, nrna_abundance=0.0,
    )
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, _sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    md = mature_density(sub, ra, e_mu, boundary_eff_length(fl.rna_pmf))
    sc.cleanup()
    return md, sub, ra, idx


def test_internal_exon_mature_recovered(tmp_path):
    md, sub, ra, idx = _run(tmp_path / "clean")
    df = idx.region_df
    sig = np.asarray(ra.signature)
    U = (sub.contained.n_unspliced_pos + sub.contained.n_unspliced_neg).astype(float)
    # the long internal exon 6000-7500 (a single region in this layout); find it by coordinate +
    # the + exon bit (the toy genome is one reference, so the coordinate is unambiguous).
    start = df["start"].to_numpy()
    end = df["end"].to_numpy()
    cand = np.flatnonzero((start <= 6500) & (end > 6500) & ((sig & BIT_EXON_POS) != 0))
    assert cand.size == 1, f"expected one + exon region over 6500, found {cand.size}"
    i = int(cand[0])
    pred = md.mature_pos[i]
    actual = U[i]  # no gDNA / nascent ⇒ contained-unspliced is all mature
    assert actual > 100, f"sanity: the long exon should carry substantial mature (U={actual})"
    assert 0.80 <= pred / actual <= 1.20, f"M_pos/U = {pred / actual:.2f} (predicted {pred:.0f} vs {actual:.0f})"


def test_no_mature_off_exon_and_wrong_strand(tmp_path):
    md, sub, ra, idx = _run(tmp_path / "clean2")
    sig = np.asarray(ra.signature)
    non_exon = (sig & _EXON) == 0
    # intron / intergenic regions must carry NO predicted mature (no junction-anchored exon there).
    assert np.allclose(md.mature_pos[non_exon], 0.0), "mature predicted off-exon"
    assert np.allclose(md.mature_neg[non_exon], 0.0), "mature predicted off-exon (−)"
    # the + gene's exon regions must have ~0 on the − field (per-strand separation); allow tiny leakage.
    pos_exon = (sig & BIT_EXON_POS) != 0
    neg_exon = (sig & BIT_EXON_NEG) != 0
    pos_only = pos_exon & ~neg_exon
    assert np.allclose(md.mature_neg[pos_only], 0.0), "− mature predicted on a +-only exon"
