"""The accumulator deposits the molecule's contiguous genomic spans, not the sequenced read blocks.

That distinction makes the two density estimators agree, and each test holds one case of it: an unspliced pair becomes
one [min, max] span with the mate gap filled, so the boundary-crossing estimator recovers the same uniform density as
the exact contained estimator; an implicitly spliced pair (the sj in the unsequenced gap, so no CIGAR N) routes to the
sj axis with the intron skipped rather than filled; and an artifact splice is held out of the tally, because its true
span is unrecoverable. Each scenario is small and seeded so the suite can afford a real sim-scan per test.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.fl import build_fl_models
from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG


def count_observable_masks(signature, ref_id) -> tuple[np.ndarray, np.ndarray]:
    """Which objects hold a count that only gDNA can make: a region with no exon bit, and a contiguous
    boundary whose two flanks share none (a shared exon strand continues across it, so unspliced RNA
    crosses it inside that exon). ``(region_mask[N], boundary_mask[E])``, each on its own axis."""
    sig = np.asarray(signature).astype(np.int64)
    lo, hi = boundary_region_indices(ref_id)
    return (sig & _EXON_BITS) == 0, (sig[lo] & sig[hi] & _EXON_BITS) == 0


def _crossing_vs_contained_ratio(bam_path, index) -> float:
    """Return crossing-ρ / contained-ρ over count-observable objects (≈1.0 if unbiased). Both densities are ratios of
    sums pooled over their own axis (``rho = Σ count / Σ opportunity``); a mean of per-boundary ratios would not be
    comparable with the contained density. A contiguous boundary is a 0-bp boundary carrying one crossing count, so it
    is one number.
    """
    _s, sm, _b, pl = scan_and_buffer(str(bam_path), index, BamScanConfig(sj_strand_tag="auto"))
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    CalibrationSubstrate._check_alignment(pl, ra)
    # The length models from the payload alone, with no index inputs, so the gDNA pmf is fitted from
    # the contained pair (production's `rigel.pipeline.library_fl_models` also passes the index's
    # opportunities and partition).
    flm = build_fl_models(pl)
    gpmf = flm.gdna_pmf
    region_eff = contained_eff_length(ra.region_size_bp, gpmf)
    # gDNA's template is the chromosome, so a boundary's divisor is the unbounded-reach limit mu_g − 1, the same at
    # every boundary, so the crossing pool needs only a count of boundaries.
    boundary_eff = float(crossing_eff_length(gpmf, [UNBOUNDED_REACH], [UNBOUNDED_REACH])[0])
    sub = CalibrationSubstrate.from_payload(pl, ra)
    region_obs, boundary_obs = count_observable_masks(
        np.asarray(ra.signature), np.asarray(ra.ref_id)
    )

    live_regions = region_obs & (region_eff > 1.0)
    contained = np.asarray(sub.region_contained.count, np.float64).sum(1)
    rho_contained = float(contained[live_regions].sum() / region_eff[live_regions].sum())

    crossing = np.asarray(sub.boundary_unspliced.count, np.float64).sum(1)
    rho_crossing = float(crossing[boundary_obs].sum() / (boundary_eff * int(boundary_obs.sum())))
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
    """An implicitly spliced fragment (an annotated intron inside the mate gap, so no CIGAR N) must deposit on the
    spliced axis and skip the intron rather than fill [min, max] across it on the unspliced channels, which would read
    downstream as intronic gDNA.
    """
    import numpy as np

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.signature import coarse_type_array
    from rigel.sim import Scenario
    from rigel.sim.reads import ReadSimConfig

    sc = Scenario("implicit_chan", genome_length=6000, seed=7, work_dir=str(tmp_path / "sim"))
    # Short intron (50 bp ≪ insert) → fragments spanning both exons place the
    # sj in the unsequenced mate gap = implicit splice (no CIGAR-N).
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
    _s, sm, _b, pl = scan_and_buffer(
        str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto")
    )
    ra = RegionArrays.from_frame(res.index.regions_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    ctype = coarse_type_array(np.asarray(ra.signature))
    rid = np.asarray(ra.ref_id)
    gene = rid == res.index.ref_name_to_id["implicit_chan"]

    # (1) Implicit splices route to the sj axis: the molecule jumped, it did not cross. A splice deposits on its sj
    #     boundary only, never on the contiguous boundaries it splices over, so the evidence is `sj_count` rather than a
    #     spliced channel on a boundary.
    sj_flux = int(np.asarray(sub.sj.count, np.int64).sum())
    assert sj_flux > 1000, f"expected substantial sj flux, got {sj_flux}"

    # (2) The intron is skipped, not filled: the intron region carries no contained mass, and neither boundary bounding
    #     it carries an unspliced crossing.
    from rigel.calibration.region_arrays import region_right_boundary

    intron = gene & (ctype == 1)
    assert intron.any()
    contained = np.asarray(sub.region_contained.count, np.float64).sum(1)
    crossing = np.asarray(sub.boundary_unspliced.count, np.float64).sum(1)
    right_boundary = region_right_boundary(np.asarray(ra.ref_id))
    bounding = np.zeros(crossing.shape[0], dtype=bool)
    for r in np.flatnonzero(intron):
        for e in (right_boundary[r], right_boundary[r - 1] if r > 0 else -1):
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
    """``SPLICE_ARTIFACT`` fragments (a blacklisted CIGAR-N sj) are held out of the accumulator entirely (no contained
    or boundary mass) while other fragments deposit normally, so total mass == n_deposited, artifacts excluded.

    Construction: 20 plain unspliced pairs (positive control) plus 50 pairs spanning an annotated sj. Without a
    blacklist all 70 deposit (the 50 as spliced); blacklisting that sj turns the 50 into artifacts (sj removed →
    n_sj_blacklisted > 0 → SPLICE_ARTIFACT), which must then deposit nothing.
    """
    import numpy as np
    import pandas as pd
    import pysam

    from _index_builder import rebuild_with_splice_blacklist
    from rigel.index import TranscriptIndex

    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 1000 + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "a.gtf"
    # 2-exon transcript → annotated sj at 0-based [200, 300).
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
        for i in range(50):  # span the annotated sj [200, 300)
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
        """Every deposit on every axis. ``region_start_count`` is the fragment tally, one per accepted fragment, and is
        what "mass == n_deposited" means, mass being a count."""
        return float(np.asarray(payload.region_start_count, np.int64).sum())

    cfg = BamScanConfig(sj_strand_tag="XS")

    # Control: no blacklist → all 70 deposit (the 50 as annotated splices).
    idx = TranscriptIndex.load(str(idx_dir))
    s0, _, _, pl0 = scan_and_buffer(str(bam), idx, cfg)
    assert s0.n_sj_blacklisted == 0
    assert s0.n_with_annotated_sj == 50
    assert total_mass(pl0) > 65.0  # 20 contained + ~50 spliced crossing mass

    # Blacklist the annotated sj → the 50 become SPLICE_ARTIFACT → held out.
    rebuild_with_splice_blacklist(
        idx_dir,
        pd.DataFrame(
            {
                "ref": ["chr1"],
                "start": np.array([200], np.int32),
                "end": np.array([300], np.int32),
                "max_anchor_left": np.array([10000], np.int32),
                "max_anchor_right": np.array([10000], np.int32),
            }
        ),
    )
    idx2 = TranscriptIndex.load(str(idx_dir))
    s1, _, _, pl1 = scan_and_buffer(str(bam), idx2, cfg)
    assert s1.n_sj_blacklisted == 50, "blacklist did not flag the sj"
    # Only the 20 non-artifact unspliced controls remain → mass == n_deposited.
    assert total_mass(pl1) == 20.0, f"artifacts not held out (mass={total_mass(pl1)})"
