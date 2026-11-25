"""The output directory is one run's (``rigel.outputs``): a run retires the previous run's completion marker and
products before writing, writes each file whole and its summary last; export converts only a complete run's own
products; the annotated BAM may name no file the run reads or writes, and is published only once it agrees with the
scan."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from _two_genes import three_pairs, two_gene_index

from rigel import cli, quant_writer
from rigel.report.substrate import SubstrateError, load_substrate


@pytest.fixture
def run(tmp_path):
    """A two-gene index, a BAM of three pairs, and ``quant(offset, *options)``: the pairs placed at A (offset 0) or B
    (offset 1000), quantified into one output directory through the CLI."""
    index = two_gene_index(tmp_path)
    reads, out = tmp_path / "reads.bam", tmp_path / "results"
    parser = cli.build_parser()

    def quant(offset: int, *options: str) -> int:
        three_pairs(reads, offset)
        argv = [
            "quant",
            "--bam",
            str(reads),
            "--index",
            str(index.index_dir),
            "-o",
            str(out),
            "--threads",
            "1",
        ]
        return cli.quant_command(parser.parse_args([*argv, *options]))

    return quant, out, index, parser


def _counts(path) -> dict:
    return pd.read_feather(path) if str(path).endswith(".feather") else pd.read_parquet(path)


@pytest.mark.parametrize("target", ["summary.json", "quant.feather", "index"])
def test_an_annotated_bam_naming_a_file_the_run_reads_or_writes_is_refused(run, tmp_path, target):
    """Before anything is written: a summary, a table or an index file as the annotated BAM. A previous valid
    summary keeps its bytes."""
    quant, out, index, _ = run
    quant(0)
    before = (out / "summary.json").read_bytes()
    path = sorted(Path(index.index_dir).iterdir())[0] if target == "index" else out / target
    with pytest.raises(SystemExit, match="a file this run reads or writes"):
        quant(0, "--annotated-bam", str(path))
    assert (out / "summary.json").read_bytes() == before


def test_an_annotated_bam_through_a_link_to_a_product_is_refused(run, tmp_path):
    quant, out, _index, _ = run
    quant(0)
    link = tmp_path / "link.bam"
    link.symlink_to(out / "gene_quant.feather")
    with pytest.raises(SystemExit, match="a file this run reads or writes"):
        quant(0, "--annotated-bam", str(link))


def test_an_interrupted_rerun_leaves_no_run_a_reader_accepts(run):
    """A rerun stopped while writing its track: the previous run's summary was retired first, so the report refuses
    the directory as incomplete instead of reading this run's counts beside the previous run's track."""
    quant, out, _index, _ = run
    quant(1000)
    with patch.object(quant_writer, "_write_calibration_track", side_effect=OSError("interrupted")):
        with pytest.raises(OSError, match="interrupted"):
            quant(0)
    with pytest.raises(SubstrateError, match="no complete run"):
        load_substrate(out)


def test_an_export_does_not_outlive_its_run(run):
    """Export, then a rerun on other reads: the rerun retires the export with its other products, so no Parquet of
    the previous run survives beside the new Feather."""
    quant, out, _index, parser = run
    quant(0)
    assert cli.export_command(parser.parse_args(["export", str(out), "-f", "parquet"])) == 0
    assert (out / "quant.parquet").is_file()
    quant(1000)
    assert not (out / "quant.parquet").exists()


def test_export_converts_only_the_runs_own_products(run):
    """A user's own Feather in the directory is neither converted nor removed."""
    quant, out, _index, parser = run
    quant(0)
    pd.DataFrame({"x": [1]}).to_feather(out / "mine.feather")
    assert cli.export_command(parser.parse_args(["export", str(out), "-f", "tsv"])) == 0
    assert not (out / "mine.tsv").exists()
    quant(1000)
    assert (out / "mine.feather").is_file()


def test_export_refuses_a_directory_without_a_complete_run(tmp_path):
    parser = cli.build_parser()
    (tmp_path / "quant.feather").write_bytes(b"")
    assert cli.export_command(parser.parse_args(["export", str(tmp_path), "-f", "tsv"])) == 1


def test_a_failed_second_pass_leaves_a_good_annotated_bam_as_it_was(run, tmp_path):
    """The annotated BAM is written beside its path and moved into place only once it agrees with the scan: a pass
    that fails leaves the good one already there byte for byte."""
    from rigel.config import BamScanConfig, CalibrationConfig, EMConfig, PipelineConfig
    from rigel.pipeline import run_pipeline

    quant, out, index, _ = run
    tagged = tmp_path / "tagged.bam"
    quant(0, "--annotated-bam", str(tagged))
    good = tagged.read_bytes()
    reads, moved = tmp_path / "reads.bam", tmp_path / "moved.bam"
    three_pairs(moved, 1000)
    config = PipelineConfig(
        scan=BamScanConfig(total_threads=1),
        em=EMConfig(n_threads=1),
        calibration=CalibrationConfig(n_threads=1),
        annotated_bam_path=tagged,
    )
    with pytest.raises(RuntimeError, match="did not see what the scan saw"):
        run_pipeline(
            reads, index, config=config, before_annotated_bam=lambda _result: moved.replace(reads)
        )
    assert tagged.read_bytes() == good
