"""Common bad inputs end in a short error naming the file, through the real CLI: a manifest that is not JSON or not a
mapping, a corrupt table, a text file given as a BAM. Each escaped as a raw parser traceback."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from _two_genes import three_pairs, two_gene_index

from rigel import cli


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("named")
    index = two_gene_index(root)
    reads = root / "reads.bam"
    three_pairs(reads, 0)
    return Path(index.index_dir), reads


def _quant(bam, index_dir, out):
    parser = cli.build_parser()
    return cli.quant_command(
        parser.parse_args(
            [
                "quant",
                "--bam",
                str(bam),
                "--index",
                str(index_dir),
                "-o",
                str(out),
                "--threads",
                "1",
            ]
        )
    )


@pytest.mark.parametrize(
    "damage, said",
    [
        (
            lambda d: (d / "manifest.json").write_text("{"),
            "manifest.json is not a readable JSON manifest",
        ),
        (lambda d: (d / "manifest.json").write_text("3"), "manifest.json holds no mapping"),
        (
            lambda d: (d / "transcripts.feather").write_bytes(b"not a feather file"),
            "could not be read",
        ),
    ],
    ids=["truncated-json", "a-json-scalar", "a-corrupt-table"],
)
def test_a_damaged_index_is_a_named_error(built, tmp_path, damage, said):
    index_dir, reads = built
    copy = tmp_path / "index"
    shutil.copytree(index_dir, copy)
    damage(copy)
    with pytest.raises(SystemExit, match=said):
        _quant(reads, copy, tmp_path / "out")


def test_a_text_file_given_as_a_bam_is_a_named_error(built, tmp_path):
    index_dir, _ = built
    text = tmp_path / "reads.bam"
    text.write_text("this is not a BAM\n")
    with pytest.raises(SystemExit, match="cannot be read as a BAM"):
        _quant(text, index_dir, tmp_path / "out")
