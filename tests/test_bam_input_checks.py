"""A BAM the scan would read wrong is refused before its records are trusted.

Each case would otherwise produce a full output set and exit 0: a coordinate-sorted BAM, whose mates are not adjacent,
so every mate would be its own fragment; references the index does not hold, so every fragment would be intergenic; a
reference the index holds at another length, another assembly; and a truncated BAM, whose read error would end the scan
like the end of file.
"""

import logging
import re

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


@pytest.mark.parametrize("hard_link", [False, True], ids=["same-path", "hard-link"])
def test_an_annotated_bam_that_is_the_input_is_refused_before_the_scan(built, tmp_path, hard_link):
    """`run_pipeline` refuses it too, for a caller that is not the CLI: the annotation pass would truncate the input
    it reads, and its check against the scan would then delete it. The input keeps its bytes."""
    import os
    import shutil

    bam = tmp_path / "in.bam"
    shutil.copyfile(built.bam_path, bam)
    before = bam.read_bytes()
    annotated = bam
    if hard_link:
        annotated = tmp_path / "linked.bam"
        os.link(bam, annotated)
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="ts"),
        annotated_bam_path=annotated,
    )
    with pytest.raises(BamInputError, match="is the input BAM"):
        run_pipeline(bam, built.index, config=config)
    assert bam.read_bytes() == before


def _rewrite(src, dst, edit):
    """Copy a BAM with each record passed through ``edit(record)``, which returns it, changed or not, or ``None`` to
    drop it."""
    with pysam.AlignmentFile(str(src), "rb") as fin:
        with pysam.AlignmentFile(str(dst), "wb", header=fin.header) as fout:
            for r in fin.fetch(until_eof=True):
                r = edit(r)
                if r is not None:
                    fout.write(r)


def _first_scanned(r) -> bool:
    """Whether the scan reads this record past its filter: mapped, not QC-failed, not a duplicate."""
    return not (r.is_unmapped or r.is_qcfail or r.is_duplicate)


def test_an_unpaired_record_is_refused_as_an_error_naming_it(built, tmp_path):
    """Rigel reads paired-end libraries only. A record the flag calls unpaired stops the scan, and the CLI reports it
    as an error naming the record, not as a traceback."""
    from rigel.cli import build_parser

    named: list[str] = []

    def unpair(r):
        if not named and _first_scanned(r):
            r.flag &= ~0x1
            named.append(r.query_name)
        return r

    unpaired = tmp_path / "unpaired.bam"
    _rewrite(built.bam_path, unpaired, unpair)
    with pytest.raises(BamInputError, match=re.escape(named[0])):
        _run(unpaired, built.index)
    args = build_parser().parse_args(
        [
            "quant",
            "--bam",
            str(unpaired),
            "--index",
            str(built.index_dir),
            "-o",
            str(tmp_path / "out"),
        ]
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    message = str(exited.value.code)
    assert message.startswith("Error: ") and "paired" in message and named[0] in message


def _pairs_lacking_a_mate(bam) -> int:
    """The read pairs (runs of records sharing a name) a record of which calls its mate mapped, while the run holds no
    record of that mate: counted from the BAM itself, apart from the scan."""
    n, name, seen, expected = 0, None, 0, 0
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if r.query_name != name:
                n += bool(expected & ~seen)
                name, seen, expected = r.query_name, 0, 0
            segment = 1 if r.is_read1 else 2 if r.is_read2 else 0
            seen |= segment
            if r.is_paired and not r.mate_is_unmapped:
                expected |= {1: 2, 2: 1}.get(segment, 0)
    return n + bool(expected & ~seen)


def _scan_counts(bam, index):
    from rigel.pipeline import scan_and_buffer

    stats, _strand_models, buffer, _payload = scan_and_buffer(
        str(bam), index, BamScanConfig(sj_strand_tag="ts")
    )
    buffer.cleanup()
    return stats


def test_a_name_grouped_bam_lacks_no_mate(built, caplog):
    with caplog.at_level(logging.WARNING):
        stats = _scan_counts(built.bam_path, built.index)
    assert stats.n_read_pairs_lacking_mate == 0 == _pairs_lacking_a_mate(built.bam_path)
    assert "lack" not in caplog.text


def _drop_some_second_mates():
    """An edit for `_rewrite`: every seventh read name, in the BAM's order, loses its second mate's records."""
    order: dict[str, int] = {}

    def edit(r):
        k = order.setdefault(r.query_name, len(order))
        return None if r.is_read2 and k % 7 == 0 else r

    return edit


def _unsorted_header_coordinate_order(src, dst, tmp_path):
    """The BAM sorted by coordinate, its header no longer saying so: the case the header check cannot see."""
    coordinate = tmp_path / "coordinate.bam"
    pysam.sort("-o", str(coordinate), str(src))
    _reheader(coordinate, dst, lambda header: header.get("HD", {}).pop("SO", None))


@pytest.mark.parametrize("case", ["mates removed", "sorted by coordinate, unsaid"])
def test_read_pairs_lacking_a_mate_their_flags_call_mapped_are_counted_and_warned(
    built, tmp_path, caplog, case
):
    """A record whose flag calls its mate mapped, in a run of records holding no record of that mate: the mates are
    elsewhere in the BAM (sorted by coordinate, with no @HD SO to say so) or were filtered out record by record. The
    scan counts such read pairs and warns with the count; no threshold decides anything."""
    bam = tmp_path / "input.bam"
    if case == "mates removed":
        _rewrite(built.bam_path, bam, _drop_some_second_mates())
    else:
        _unsorted_header_coordinate_order(built.bam_path, bam, tmp_path)
    expected = _pairs_lacking_a_mate(bam)
    assert expected > 0
    with caplog.at_level(logging.WARNING):
        stats = _scan_counts(bam, built.index)
    assert stats.n_read_pairs_lacking_mate == expected
    assert f"{expected:,} of {stats.n_read_pairs:,} read pairs lack a mate" in caplog.text
