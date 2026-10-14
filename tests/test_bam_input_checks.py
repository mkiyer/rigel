"""A BAM the scan would read wrong is refused before its records are trusted.

Each case once produced a full output set and exit 0: a coordinate-sorted BAM, whose mates are not adjacent, so every
mate was its own fragment; references the index does not hold, so every fragment was intergenic; a reference the index
holds at another length, another assembly; and a truncated BAM, whose read error ended the scan like the end of file.
"""

import logging

import pysam
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import BamInputError, run_pipeline
from rigel.sim import ReadSimConfig, Scenario

SEED = 11


@pytest.fixture
def built(tmp_path):
    sc = Scenario("bam_input", genome_length=5000, seed=SEED, work_dir=tmp_path / "sc")
    sc.add_gene(
        "g1",
        "+",
        [{"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80}],
    )
    result = sc.build(
        n_fragments=400,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=SEED,
        ),
    )
    yield result
    sc.cleanup()


def _run(bam, index):
    return run_pipeline(
        bam,
        index,
        config=PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="ts")),
    )


def _reheader(src, dst, edit):
    """Copy a BAM with ``edit(header_dict)`` applied to its header; records keep their reference ids."""
    with pysam.AlignmentFile(str(src), "rb") as fin:
        header = fin.header.to_dict()
        edit(header)
        with pysam.AlignmentFile(str(dst), "wb", header=header) as fout:
            for r in fin.fetch(until_eof=True):
                fout.write(r)


def test_a_coordinate_sorted_bam_is_refused(built, tmp_path):
    coordinate = tmp_path / "coordinate.bam"
    pysam.sort("-o", str(coordinate), str(built.bam_path))
    with pytest.raises(BamInputError, match="sorted by coordinate"):
        _run(coordinate, built.index)


def test_references_the_index_does_not_hold_are_refused(built, tmp_path):
    def prefix(header):
        for sq in header["SQ"]:
            sq["SN"] = "chr" + sq["SN"]

    renamed = tmp_path / "renamed.bam"
    _reheader(built.bam_path, renamed, prefix)
    with pytest.raises(BamInputError, match="None of"):
        _run(renamed, built.index)


def test_a_reference_at_another_length_is_refused(built, tmp_path):
    def lengthen(header):
        header["SQ"][0]["LN"] += 1

    relengthed = tmp_path / "relengthed.bam"
    _reheader(built.bam_path, relengthed, lengthen)
    with pytest.raises(BamInputError, match="another length"):
        _run(relengthed, built.index)


def test_references_beyond_the_index_are_named_in_a_warning(built, tmp_path, caplog):
    def add_decoy(header):
        header["SQ"].append({"SN": "chrUn_decoy", "LN": 1000})

    extended = tmp_path / "extended.bam"
    _reheader(built.bam_path, extended, add_decoy)
    with caplog.at_level(logging.WARNING):
        _run(extended, built.index)
    assert any("chrUn_decoy" in r.getMessage() for r in caplog.records)


def test_a_truncated_bam_is_refused(built, tmp_path):
    truncated = tmp_path / "truncated.bam"
    data = built.bam_path.read_bytes()
    truncated.write_bytes(data[: len(data) // 2])
    with pytest.raises(BamInputError, match="truncated"):
        _run(truncated, built.index)


def test_a_read_error_inside_the_bam_ends_the_run_not_the_scan(built, tmp_path):
    """A BAM corrupt inside, its end-of-file marker intact, opens; its read error must fail the run rather than end
    the scan on the records read before it. Decompression runs on the reader's own thread here: threaded read-ahead
    can meet the corrupt block while the header is read, which fails the run too, but not through the record loop
    this gate holds."""
    corrupt = tmp_path / "corrupt.bam"
    data = bytearray(built.bam_path.read_bytes())
    middle = len(data) // 2
    data[middle : middle + 4096] = bytes(4096)
    corrupt.write_bytes(bytes(data))
    config = PipelineConfig(
        em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="ts", bgzf_threads=0)
    )
    with pytest.raises(RuntimeError, match="truncated or corrupt"):
        run_pipeline(corrupt, built.index, config=config)


def test_the_cli_reports_a_refused_bam_as_an_error(built, tmp_path):
    from rigel.cli import build_parser

    coordinate = tmp_path / "coordinate.bam"
    pysam.sort("-o", str(coordinate), str(built.bam_path))
    args = build_parser().parse_args(
        [
            "quant",
            "--bam",
            str(coordinate),
            "--index",
            str(built.index_dir),
            "-o",
            str(tmp_path / "out"),
        ]
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert str(exited.value.code).startswith("Error: ") and "sorted by coordinate" in str(
        exited.value.code
    )
