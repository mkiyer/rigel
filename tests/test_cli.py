"""CLI parser tests."""

from hulkrna.cli import build_parser


def test_index_gtf_parse_mode_default_strict():
    parser = build_parser()
    args = parser.parse_args([
        "index",
        "--fasta", "a.fa",
        "--gtf", "a.gtf",
        "--output-dir", "out",
    ])
    assert args.gtf_parse_mode == "strict"


def test_index_gtf_parse_mode_warn_skip():
    parser = build_parser()
    args = parser.parse_args([
        "index",
        "--fasta", "a.fa",
        "--gtf", "a.gtf",
        "--output-dir", "out",
        "--gtf-parse-mode", "warn-skip",
    ])
    assert args.gtf_parse_mode == "warn-skip"
