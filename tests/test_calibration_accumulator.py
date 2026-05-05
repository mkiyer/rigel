"""Tests for the M3 calibration accumulator path (C++ + Python).

Covers:

* ``BamScanner.set_regions`` binding contract
* End-to-end pipeline payload attachment
* ``CalibrationScanPayload.from_scan_dict`` validation (shape, dtype,
  balance assertion)
* Worker-merge equality (1 vs 4 workers)

See ``docs/calibration/m3_implementation_plan.md`` §5.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.scan_payload import CalibrationScanPayload
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.native import BamScanner
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import Scenario, SimConfig

SEED = 42


# ---------------------------------------------------------------------------
# Fixture: small oracle scenario for end-to-end + worker-merge tests
# ---------------------------------------------------------------------------


@pytest.fixture
def calib_scenario(tmp_path):
    """Mini oracle scenario; enough genes/reads to populate calibration."""
    sc = Scenario(
        "calib_acc_test",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "calib_acc",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    sim_config = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=SEED,
    )
    result = sc.build_oracle(n_fragments=500, sim_config=sim_config)
    yield sc, result
    sc.cleanup()


# ---------------------------------------------------------------------------
# 1. set_regions binding contract
# ---------------------------------------------------------------------------


def _make_scanner(index):
    return BamScanner(index.resolver, "XS", True, False)


def _basic_region_arrays(index):
    """Build a minimal valid (ref_ids, starts, ends, type_masks) input."""
    resolver_map = index.resolver.get_ref_to_id()
    # One INTERGENIC region [0, 100) on every known ref.
    ref_ids = np.array(sorted(resolver_map.values()), dtype=np.int32)
    n = ref_ids.size
    starts = np.zeros(n, dtype=np.int64)
    ends = np.full(n, 100, dtype=np.int64)
    type_masks = np.full(n, 0b100, dtype=np.uint8)  # INTERGENIC bit
    return ref_ids, starts, ends, type_masks


class TestSetRegions:
    def test_basic_install(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        n_refs = int(ri.max()) + 1
        scanner.set_regions(ri, s, e, tm, n_refs)  # should not raise

    def test_length_mismatch_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        with pytest.raises(Exception):
            # ends array truncated by one
            scanner.set_regions(ri, s, e[:-1], tm, int(ri.max()) + 1)

    def test_double_set_rejected(self, calib_scenario):
        _, result = calib_scenario
        scanner = _make_scanner(result.index)
        ri, s, e, tm = _basic_region_arrays(result.index)
        n_refs = int(ri.max()) + 1
        scanner.set_regions(ri, s, e, tm, n_refs)
        with pytest.raises(Exception):
            scanner.set_regions(ri, s, e, tm, n_refs)


# ---------------------------------------------------------------------------
# 2. Pipeline end-to-end payload attachment + balance
# ---------------------------------------------------------------------------


class TestPipelinePayload:
    def test_scan_returns_payload(self, calib_scenario):
        _, result = calib_scenario
        scan_cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
        _stats, _sm, _flm, _buf, payload = scan_and_buffer(
            str(result.bam_path), result.index, scan_cfg,
        )
        assert payload is not None
        # global_counts sums to n_observed by validator construction
        assert int(payload.global_counts.sum()) == payload.n_observed
        # at least some EXON-only observations on a clean strand-1 sim
        EXON_ONLY = 0b001
        assert payload.global_counts[EXON_ONLY] > 0
        # FL histogram has shape (8, 1024) and sums to <= n_observed
        assert payload.fl_hist.shape == (8, 1024)
        assert int(payload.fl_hist.sum()) == payload.n_observed
        # per-region rows == n_regions
        n_regions = len(result.index.region_df)
        assert payload.per_region_counts.shape == (n_regions, 8)
        assert payload.u_left.shape == (n_regions,)
        assert payload.u_right.shape == (n_regions,)

    def test_run_pipeline_attaches_payload(self, calib_scenario, tmp_path):
        _, result = calib_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(str(result.bam_path), result.index, config)
        assert pr.calibration_payload is not None
        assert pr.calibration_payload.n_observed > 0


# ---------------------------------------------------------------------------
# 3. Payload validation (shape, dtype, balance)
# ---------------------------------------------------------------------------


def _good_payload_dict(n_regions: int = 3, n_observed: int = 10) -> dict:
    """A schema-correct, internally consistent calibration dict."""
    global_counts = np.zeros(8, dtype=np.int64)
    global_counts[1] = n_observed  # all EXON_ONLY
    fl_hist = np.zeros((8, 1024), dtype=np.int64)
    fl_hist[1, 200] = n_observed
    return {
        "global_counts": global_counts,
        "per_region_counts": np.zeros((n_regions, 8), dtype=np.int64),
        "fl_hist": fl_hist,
        "u_left": np.zeros(n_regions, dtype=np.int64),
        "u_right": np.zeros(n_regions, dtype=np.int64),
        "n_observed": n_observed,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_unobserved": 0,
        "n_unannotated_ref": 0,
    }


class TestPayloadValidation:
    def test_good_dict_round_trips(self):
        d = _good_payload_dict()
        p = CalibrationScanPayload.from_scan_dict(d, n_total=10)
        assert p.n_observed == 10

    def test_balance_violation_raises(self):
        d = _good_payload_dict(n_observed=10)
        # Pretend total was 15 (5 unaccounted)
        with pytest.raises(ValueError, match="balance assertion failed"):
            CalibrationScanPayload.from_scan_dict(d, n_total=15)

    def test_dtype_mismatch_raises(self):
        d = _good_payload_dict()
        d["global_counts"] = d["global_counts"].astype(np.int32)
        with pytest.raises(ValueError, match="dtype"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_shape_mismatch_raises(self):
        d = _good_payload_dict()
        d["fl_hist"] = np.zeros((8, 512), dtype=np.int64)
        with pytest.raises(ValueError, match="shape"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_global_counts_must_sum_to_n_observed(self):
        d = _good_payload_dict(n_observed=10)
        d["n_observed"] = 5  # tampered
        with pytest.raises(ValueError, match="global_counts"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_n_unannotated_ref_le_n_observed(self):
        d = _good_payload_dict(n_observed=10)
        d["n_unannotated_ref"] = 99
        with pytest.raises(ValueError, match="n_unannotated_ref"):
            CalibrationScanPayload.from_scan_dict(d)

    def test_none_dict_raises(self):
        with pytest.raises(ValueError, match="set_regions"):
            CalibrationScanPayload.from_scan_dict(None)


# ---------------------------------------------------------------------------
# 4. Worker-merge equality (1 vs N workers must be byte-identical)
# ---------------------------------------------------------------------------


class TestWorkerMergeEquality:
    def test_one_vs_four_workers_identical(self, calib_scenario):
        _, result = calib_scenario

        def _run(n: int) -> CalibrationScanPayload:
            cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=n)
            _, _, _, _, p = scan_and_buffer(
                str(result.bam_path), result.index, cfg,
            )
            return p

        a = _run(1)
        b = _run(4)
        assert a is not None and b is not None
        np.testing.assert_array_equal(a.global_counts, b.global_counts)
        np.testing.assert_array_equal(a.per_region_counts, b.per_region_counts)
        np.testing.assert_array_equal(a.fl_hist, b.fl_hist)
        np.testing.assert_array_equal(a.u_left, b.u_left)
        np.testing.assert_array_equal(a.u_right, b.u_right)
        assert a.n_observed == b.n_observed
        assert a.n_excluded_multimap == b.n_excluded_multimap
        assert a.n_excluded_chimera == b.n_excluded_chimera
        assert a.n_excluded_artifact == b.n_excluded_artifact
        assert a.n_unobserved == b.n_unobserved
        assert a.n_unannotated_ref == b.n_unannotated_ref
# ---------------------------------------------------------------------------
# 5. Hand-built BAM contract tests against `mini_index`
#
# These exercise the mask-correctness, boundary-flux, and observation-policy
# contracts that the M3 plan exit gate calls out but that the macro-level
# fixtures above cannot pin in isolation.
#
# `mini_index` (see conftest) lays chr1 (length 2000) out as the per-ref
# tile below.  Region indices 0..10 follow the canonical sort order:
#
#     0  INTERGENIC  [   0,   99)   mask 0b100
#     1  EXON        [  99,  200)   mask 0b001    flux: left=False, right=True
#     2  INTRON      [ 200,  299)   mask 0b010
#     3  EXON        [ 299,  400)   mask 0b001    flux: left=True,  right=True
#     4  INTRON      [ 400,  499)   mask 0b010
#     5  EXON        [ 499,  600)   mask 0b001    flux: left=True,  right=False
#     6  INTERGENIC  [ 600,  999)   mask 0b100
#     7  EXON        [ 999, 1100)   mask 0b001    flux: left=False, right=True
#     8  INTRON      [1100, 1199)   mask 0b010
#     9  EXON        [1199, 1300)   mask 0b001    flux: left=True,  right=False
#    10  INTERGENIC  [1300, 2000)   mask 0b100
# ---------------------------------------------------------------------------

import pysam  # noqa: E402

# 8-state mask bit layout (matches accumulator.cpp):
EXON_BIT = 0b001
INTRON_BIT = 0b010
INTERGENIC_BIT = 0b100


def _build_bam(tmp_path, ref_lengths, reads):
    """Write a query-name-sorted BAM containing `reads` and return its path."""
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": ref, "LN": L} for ref, L in ref_lengths],
    }
    bam_path = str(tmp_path / "hand.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam_path, bam_path)
    return bam_path


def _make_pair(
    qname, ref_id, r1_pos, r2_pos, *,
    r1_cigar=None, r2_cigar=None,
    r1_is_reverse=False, r2_is_reverse=True,
    r2_ref_id=None, nh=1, extra_tags=None,
):
    """Construct an FR-oriented paired-end read.

    Defaults: 50 bp R1 forward, 50 bp R2 reverse, both unspliced.  Pass
    ``r2_ref_id`` to make a chimeric pair.
    """
    if r1_cigar is None:
        r1_cigar = [(0, 50)]
    if r2_cigar is None:
        r2_cigar = [(0, 50)]
    if r2_ref_id is None:
        r2_ref_id = ref_id

    seq_len_r1 = sum(n for op, n in r1_cigar if op in (0, 1, 4))
    seq_len_r2 = sum(n for op, n in r2_cigar if op in (0, 1, 4))

    def _ref_span(cigar):
        return sum(n for op, n in cigar if op in (0, 2, 3, 7, 8))

    r1 = pysam.AlignedSegment()
    r1.query_name = qname
    r1.reference_id = ref_id
    r1.reference_start = r1_pos
    r1.mapping_quality = 60
    flag1 = 0x1 | 0x2 | 0x40
    if r1_is_reverse:
        flag1 |= 0x10
    if r2_is_reverse:
        flag1 |= 0x20
    r1.flag = flag1
    r1.cigar = r1_cigar
    r1.query_sequence = "A" * seq_len_r1
    r1.query_qualities = pysam.qualitystring_to_array("I" * seq_len_r1)
    r1.next_reference_id = r2_ref_id
    r1.next_reference_start = r2_pos
    tlen = (r2_pos + _ref_span(r2_cigar)) - r1_pos if r2_ref_id == ref_id else 0
    r1.template_length = tlen

    r2 = pysam.AlignedSegment()
    r2.query_name = qname
    r2.reference_id = r2_ref_id
    r2.reference_start = r2_pos
    r2.mapping_quality = 60
    flag2 = 0x1 | 0x2 | 0x80
    if r2_is_reverse:
        flag2 |= 0x10
    if r1_is_reverse:
        flag2 |= 0x20
    r2.flag = flag2
    r2.cigar = r2_cigar
    r2.query_sequence = "A" * seq_len_r2
    r2.query_qualities = pysam.qualitystring_to_array("I" * seq_len_r2)
    r2.next_reference_id = ref_id
    r2.next_reference_start = r1_pos
    r2.template_length = -tlen

    tags1 = [("NH", nh)]
    tags2 = [("NH", nh)]
    if extra_tags:
        tags1.extend(extra_tags)
        tags2.extend(extra_tags)
    r1.set_tags(tags1)
    r2.set_tags(tags2)
    return [r1, r2]


def _scan(bam_path, index, n_threads=1):
    cfg = BamScanConfig(sj_strand_tag="auto", n_scan_threads=n_threads)
    _, _, _, _, payload = scan_and_buffer(bam_path, index, cfg)
    return payload


class TestMaskCorrectness:
    """One unspliced fragment per intended mask → exactly one global_counts bin."""

    def _run_single(self, mini_index, tmp_path, r1_pos, r2_pos):
        ref_lens = [(name, int(L)) for name, L in mini_index.ref_lengths.items()]
        chr1_id = mini_index.resolver.get_ref_to_id()["chr1"]
        reads = _make_pair("frag", chr1_id, r1_pos, r2_pos)
        bam = _build_bam(tmp_path, ref_lens, reads)
        return _scan(bam, mini_index)

    def test_exon_only(self, mini_index, tmp_path):
        # R1 [110,160), R2 [145,195) — both inside EXON [99,200).
        p = self._run_single(mini_index, tmp_path, r1_pos=110, r2_pos=145)
        assert p.global_counts[EXON_BIT] == 1
        assert int(p.global_counts.sum()) == 1

    def test_intron_only(self, mini_index, tmp_path):
        # Both mates inside INTRON [200,299).  R1 [210,260), R2 [240,290).
        p = self._run_single(mini_index, tmp_path, r1_pos=210, r2_pos=240)
        assert p.global_counts[INTRON_BIT] == 1
        assert int(p.global_counts.sum()) == 1

    def test_intergenic_only(self, mini_index, tmp_path):
        # Both mates inside INTERGENIC [600,999).  R1 [650,700), R2 [700,750).
        p = self._run_single(mini_index, tmp_path, r1_pos=650, r2_pos=700)
        assert p.global_counts[INTERGENIC_BIT] == 1
        assert int(p.global_counts.sum()) == 1

    def test_exon_intron_straddle(self, mini_index, tmp_path):
        # Fragment overlaps EXON [99,200) and INTRON [200,299).
        # R1 [180,230) crosses 200; R2 [240,290) inside INTRON.
        p = self._run_single(mini_index, tmp_path, r1_pos=180, r2_pos=240)
        assert p.global_counts[EXON_BIT | INTRON_BIT] == 1
        assert int(p.global_counts.sum()) == 1

    def test_intergenic_exon_straddle(self, mini_index, tmp_path):
        # Fragment overlaps INTERGENIC [0,99) and EXON [99,200).
        # R1 [70,120) crosses 99; R2 [130,180) inside EXON.
        p = self._run_single(mini_index, tmp_path, r1_pos=70, r2_pos=130)
        assert p.global_counts[INTERGENIC_BIT | EXON_BIT] == 1
        assert int(p.global_counts.sum()) == 1


class TestBoundaryFlux:
    """Pass-D contract: an unspliced exon-touching fragment that crosses
    region edge `s` (resp. `e`) must increment u_left[r] (resp. u_right[r])
    exactly once for that exon region `r`."""

    def _run(self, mini_index, tmp_path, r1_pos, r2_pos):
        ref_lens = [(name, int(L)) for name, L in mini_index.ref_lengths.items()]
        chr1_id = mini_index.resolver.get_ref_to_id()["chr1"]
        bam = _build_bam(
            tmp_path, ref_lens,
            _make_pair("frag", chr1_id, r1_pos, r2_pos),
        )
        return _scan(bam, mini_index)

    def test_left_straddle_only(self, mini_index, tmp_path):
        # Fragment [70,180) straddles left edge (99) of EXON region 1
        # but does not reach its right edge (200).
        p = self._run(mini_index, tmp_path, r1_pos=70, r2_pos=130)
        assert p.u_left[1] == 1
        assert p.u_right[1] == 0

    def test_right_straddle_only(self, mini_index, tmp_path):
        # Fragment [180,250) straddles right edge (200) of EXON region 1
        # but does not reach its left edge (99).
        p = self._run(mini_index, tmp_path, r1_pos=180, r2_pos=200)
        assert p.u_left[1] == 0
        assert p.u_right[1] == 1

    def test_full_span_both_flags(self, mini_index, tmp_path):
        # Fragment [70,250) straddles BOTH edges of EXON region 1.
        p = self._run(mini_index, tmp_path, r1_pos=70, r2_pos=200)
        assert p.u_left[1] == 1
        assert p.u_right[1] == 1

    def test_interior_fragment_no_flux(self, mini_index, tmp_path):
        # Fragment fully inside EXON region 1 → no flux.
        p = self._run(mini_index, tmp_path, r1_pos=110, r2_pos=145)
        assert p.u_left[1] == 0
        assert p.u_right[1] == 0


class TestObservationPolicy:
    """Pin the per-class scan-time exclusion counters."""

    def test_excluded_multimap(self, mini_index, tmp_path):
        # NH=2 → multimapper.  Counted at early dispatch, never reaches observe().
        ref_lens = [(name, int(L)) for name, L in mini_index.ref_lengths.items()]
        chr1_id = mini_index.resolver.get_ref_to_id()["chr1"]
        reads = _make_pair("mm", chr1_id, 110, 145, nh=2)
        bam = _build_bam(tmp_path, ref_lens, reads)
        p = _scan(bam, mini_index)
        assert p.n_excluded_multimap == 1
        assert p.n_observed == 0

    def test_excluded_chimera(self, tmp_path_factory, tmp_path):
        # Build a 2-ref index (chr1 + chr2 each with a transcript) so the
        # resolver knows both refs.  Then place R1 on chr1 and R2 on chr2
        # — a trans-chromosomal pair that resolves only as chimeric.
        import pysam
        from rigel.index import TranscriptIndex
        gtf = (
            'chr1\tt\texon\t100\t200\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "t1"; gene_name "GA"; '
            'gene_type "protein_coding"; tag "basic";\n'
            'chr2\tt\texon\t100\t200\t.\t+\t.\t'
            'gene_id "g2"; transcript_id "t2"; gene_name "GB"; '
            'gene_type "protein_coding"; tag "basic";\n'
        )
        base = tmp_path_factory.mktemp("chim_idx")
        gtf_path = base / "test.gtf"
        gtf_path.write_text(gtf)
        fasta_path = base / "genome.fa"
        with open(fasta_path, "w") as f:
            for ref, L in [("chr1", 1000), ("chr2", 1000)]:
                f.write(f">{ref}\n")
                seq = "N" * L
                for i in range(0, L, 80):
                    f.write(seq[i:i + 80] + "\n")
        pysam.faidx(str(fasta_path))
        idx_dir = base / "index"
        TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
        idx = TranscriptIndex.load(idx_dir, retain_test_structures=True)

        ref_lens = [(name, int(L)) for name, L in idx.ref_lengths.items()]
        ref_to_id = idx.resolver.get_ref_to_id()
        chr1_id, chr2_id = ref_to_id["chr1"], ref_to_id["chr2"]
        reads = _make_pair(
            "chim", chr1_id, 120, 130,
            r2_ref_id=chr2_id,
        )
        bam = _build_bam(tmp_path, ref_lens, reads)
        p = _scan(bam, idx)
        assert p.n_excluded_chimera == 1
        assert p.n_observed == 0

    @pytest.mark.skip(
        reason="SPLICE_ARTIFACT requires splice_blacklist.feather plumbing; "
               "covered indirectly by tests/test_splice_blacklist.py and will "
               "be pinned end-to-end in M9 once the blacklist is wired into "
               "mini_index."
    )
    def test_excluded_artifact(self, mini_index, tmp_path):
        ...
