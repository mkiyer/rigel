"""Behavioral guards for the accumulator span redesign.

(``docs/CARRY_FORWARD.md``.) The accumulator
deposits the MOLECULE's contiguous genomic span(s), not the sequenced read blocks:
  - unspliced  → one [min,max] span (mate gap filled) ⇒ the boundary-crossing
    density estimator agrees with the exact contained estimator (no over-count);
  - implicit splice → spliced channel, intron cut (not filled);
  - artifact splice → held out entirely (unrecoverable true span).

Kept small + deterministic (fixed seeds) so they run fast.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType


def _crossing_vs_contained_ratio(bam_path, index) -> float:
    """Return crossing-ρ / contained-ρ over count-observable objects (≈1.0 if unbiased).

    ⭐ **Both are RATIOS OF SUMS**, pooled over their own axis (`CARRY_FORWARD.md` §2,
    ``ρ_bg = Σg/ΣE``). The predecessor took the MEAN of the per-boundary fluxes and divided by
    ``fl_mean``, which is a mean of ratios — a different number whenever the supports differ.

    ⭐ **And a line is ONE number.** The predecessor averaged a boundary's two faces,
    ``(right[r] + left[r+1]) / 2``, because the old accumulator split one crossing across them. A
    contiguous edge is a 0-bp line with one count, so the ½ and the loop over flanks both go.
    """
    _s, sm, fla, _b, pl = scan_and_buffer(str(bam_path), index, BamScanConfig(sj_strand_tag="auto"))
    ra = RegionArrays.from_frame(index.nodes_df, index.ref_name_to_id)
    CalibrationSubstrate._check_alignment(pl, ra)
    flm = build_fl_models(
        global_counts=fla.global_model.counts,
        rna_counts=fla.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl),
        max_size=fla.max_size,
    )
    gpmf = flm.gdna_pmf
    node_eff = contained_eff_length(ra.region_size_bp, gpmf)
    # gDNA's template is the chromosome, so a line's divisor is the unbounded-reach limit mu_g − 1 —
    # the SAME number at every line, which is why the crossing pool needs only a count of lines.
    edge_eff = float(crossing_eff_length(gpmf, [UNBOUNDED_REACH], [UNBOUNDED_REACH])[0])
    sub = CalibrationSubstrate.from_payload(pl, ra)
    node_obs, edge_obs = count_observable_masks(np.asarray(ra.signature), np.asarray(ra.ref_id))

    live_nodes = node_obs & (node_eff > 1.0)
    contained = np.asarray(sub.node_contained.count, np.float64).sum(1)
    rho_contained = float(contained[live_nodes].sum() / node_eff[live_nodes].sum())

    crossing = np.asarray(sub.edge_unspliced.count, np.float64).sum(1)
    rho_crossing = float(crossing[edge_obs].sum() / (edge_eff * int(edge_obs.sum())))
    return rho_crossing / rho_contained


def test_crossing_density_unbiased(tmp_path):
    from rigel.sim import Scenario
    from rigel.sim.reads import GDNAConfig, ReadSimConfig

    sc = Scenario("span_unbiased", genome_length=40000, seed=29, work_dir=str(tmp_path / "sim"))
    sc.add_gene(
        "G",
        "+",
        [
            {
                "t_id": "G.1",
                "exons": [(12000, 13000), (16000, 17000), (20000, 21000)],
                "abundance": 40,
            }
        ],
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
            frag_mean=250,
            frag_std=40,
            frag_min=180,
            frag_max=400,
            read_length=100,
            strand_specificity=1.0,
            seed=7,
        ),
    )
    _s, sm, _fl, _b, pl = scan_and_buffer(
        str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto")
    )
    ra = RegionArrays.from_frame(res.index.nodes_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    ctype = coarse_type_array(np.asarray(ra.signature))
    rid = np.asarray(ra.ref_id)
    gene = rid == res.index.ref_name_to_id["implicit_chan"]

    # (1) Implicit splices route to the JUNCTION axis — the molecule JUMPED, it did not cross. ⭐ In
    #     the new model a splice deposits on its junction edge ONLY, never on the contiguous lines it
    #     splices over, so the evidence is `sj_count` rather than a spliced channel on a boundary.
    junction_flux = int(np.asarray(sub.junction.count, np.int64).sum())
    assert junction_flux > 1000, f"expected substantial junction flux, got {junction_flux}"

    # (2) The intron is CUT, not filled: the intron NODE carries no contained mass, and neither of the
    #     lines bounding it carries an unspliced crossing — the implicit molecules skip both.
    from rigel.calibration.region_arrays import node_right_edge

    intron = gene & (ctype == 1)
    assert intron.any()
    contained = np.asarray(sub.node_contained.count, np.float64).sum(1)
    crossing = np.asarray(sub.edge_unspliced.count, np.float64).sum(1)
    right_edge = node_right_edge(np.asarray(ra.ref_id))
    bounding = np.zeros(crossing.shape[0], dtype=bool)
    for r in np.flatnonzero(intron):
        for e in (right_edge[r], right_edge[r - 1] if r > 0 else -1):
            if e >= 0:
                bounding[e] = True
    intron_unspliced = float(contained[intron].sum() + crossing[bounding].sum())
    assert intron_unspliced == 0.0, f"intron carries unspliced mass {intron_unspliced} (not cut)"


def _write_pair(out, qn, *, r1_pos, r1_cigar, r2_pos, r2_cigar, xs):
    import re

    import pysam

    def seg(is_r1, pos, cigar, mpos):
        a = pysam.AlignedSegment()
        a.query_name = qn
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = 60
        a.cigarstring = cigar
        a.next_reference_id = 0
        a.next_reference_start = mpos
        a.query_sequence = "A" * sum(int(x) for x in re.findall(r"(\d+)[MIS]", cigar))
        a.flag = (0x1 | 0x2) | (0x40 if is_r1 else 0x80) | (0x20 if is_r1 else 0x10)
        tags = [("NH", 1)] + ([("XS", "+", "A")] if (is_r1 and xs) else [])
        a.set_tags(tags)
        return a

    out.write(seg(True, r1_pos, r1_cigar, r2_pos))
    out.write(seg(False, r2_pos, r2_cigar, r1_pos))


def test_artifact_splice_held_out_and_mass_conserved(tmp_path):
    """``SPLICE_ARTIFACT`` fragments (a blacklisted CIGAR-N junction) are held out
    of the accumulator entirely — no contained/boundary mass — while non-artifact
    fragments deposit normally, so total mass == n_deposited (artifacts excluded).

    Construction: 20 plain unspliced pairs (positive control) + 50 pairs spanning
    an ANNOTATED junction. Without a blacklist all 70 deposit (the 50 as spliced);
    blacklisting that junction turns the 50 into artifacts (junction removed →
    n_sj_blacklisted > 0 → SPLICE_ARTIFACT) which must then deposit nothing.
    """
    import numpy as np
    import pandas as pd
    import pysam

    from rigel.index import SJ_BLACKLIST_FEATHER, TranscriptIndex

    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 1000 + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "a.gtf"
    # 2-exon transcript → annotated junction at 0-based [200, 300).
    gtf.write_text(
        'chr1\tsrc\texon\t1\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t301\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    idx_dir = tmp_path / "idx"
    TranscriptIndex.build(fasta_file=str(fasta), gtf_file=str(gtf), output_dir=str(idx_dir))

    bam = tmp_path / "art.bam"
    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 4000}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(20):  # unspliced positive control (exon 1)
            _write_pair(
                out, f"u{i:03d}", r1_pos=20, r1_cigar="100M", r2_pos=120, r2_cigar="80M", xs=False
            )
        for i in range(50):  # span the annotated junction [200, 300)
            _write_pair(
                out,
                f"r{i:03d}",
                r1_pos=150,
                r1_cigar="50M100N50M",
                r2_pos=350,
                r2_cigar="80M",
                xs=True,
            )

    def total_mass(payload) -> float:
        """Every deposit on every axis. ⭐ ``node_start_count`` is the fragment tally — one per accepted
        fragment — and is what "mass == n_deposited" means now that mass IS a count."""
        return float(np.asarray(payload.node_start_count, np.int64).sum())

    cfg = BamScanConfig(sj_strand_tag="XS")

    # Control: no blacklist → all 70 deposit (the 50 as annotated splices).
    idx = TranscriptIndex.load(str(idx_dir))
    s0, _, _, _, pl0 = scan_and_buffer(str(bam), idx, cfg)
    assert s0.n_sj_blacklisted == 0
    assert s0.n_with_annotated_sj == 50
    assert total_mass(pl0) > 65.0  # 20 contained + ~50 spliced crossing mass

    # Blacklist the annotated junction → the 50 become SPLICE_ARTIFACT → held out.
    pd.DataFrame(
        {
            "ref": ["chr1"],
            "start": np.array([200], np.int32),
            "end": np.array([300], np.int32),
            "max_anchor_left": np.array([10000], np.int32),
            "max_anchor_right": np.array([10000], np.int32),
        }
    ).to_feather(idx_dir / SJ_BLACKLIST_FEATHER)
    idx2 = TranscriptIndex.load(str(idx_dir))
    s1, _, _, _, pl1 = scan_and_buffer(str(bam), idx2, cfg)
    assert s1.n_sj_blacklisted == 50, "blacklist did not flag the junction"
    # Only the 20 non-artifact unspliced controls remain → mass == n_deposited.
    assert total_mass(pl1) == 20.0, f"artifacts not held out (mass={total_mass(pl1)})"
