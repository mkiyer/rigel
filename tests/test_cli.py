"""The command line is the only interface most runs have, so every flag must reach the config it
names.

`rigel index`'s GTF parse mode defaults to strict and can be switched to warn-and-skip; `rigel
quant`'s defaults, its boolean flags and the resolution of its arguments into a `PipelineConfig` are
gated flag by flag, and the config survives a write/read round trip unchanged; an unknown YAML key
is refused by `rigel quant` and `rigel sim`. A flag parsed into a field nothing reads is invisible
at runtime, which is why these are checked against the resolved config rather than against the
parser's namespace.
"""

import textwrap

import pytest

from rigel.cli import build_parser, _resolve_quant_args, _build_quant_defaults
from rigel.config import BamScanConfig


# ---------------------------------------------------------------------------
# Helper: parse quant subcommand with minimal required args
# ---------------------------------------------------------------------------

_QUANT_REQ = ["quant", "--bam", "x.bam", "--index", "idx", "-o", "out"]


def _parse_quant(*extra_args):
    """Parse quant subcommand with required I/O + optional extras."""
    parser = build_parser()
    return parser.parse_args([*_QUANT_REQ, *extra_args])


# ---------------------------------------------------------------------------
# Index tests
# ---------------------------------------------------------------------------


def test_index_gtf_parse_mode_default_strict():
    parser = build_parser()
    args = parser.parse_args(
        [
            "index",
            "--fasta",
            "a.fa",
            "--gtf",
            "a.gtf",
            "--output-dir",
            "out",
        ]
    )
    assert args.gtf_parse_mode == "strict"


def test_index_gtf_parse_mode_warn_skip():
    parser = build_parser()
    args = parser.parse_args(
        [
            "index",
            "--fasta",
            "a.fa",
            "--gtf",
            "a.gtf",
            "--output-dir",
            "out",
            "--gtf-parse-mode",
            "warn-skip",
        ]
    )
    assert args.gtf_parse_mode == "warn-skip"


def test_sim_unknown_scenario_key_is_refused(tmp_path):
    """``n_fragment`` for ``n_fragments`` must not simulate the default fragment count."""
    cfg = tmp_path / "scenario.yaml"
    cfg.write_text("genome_length: 5000\nn_fragment: 10\n")
    args = build_parser().parse_args(["sim", "--config", str(cfg), "-o", str(tmp_path / "out")])
    with pytest.raises(ValueError, match=r"\['n_fragment'\]"):
        args.func(args)


def test_sim_accepts_every_documented_scenario_key(tmp_path):
    """Every top-level key docs/MANUAL.md lists for a ``rigel sim`` scenario is accepted."""
    cfg = tmp_path / "scenario.yaml"
    cfg.write_text(
        textwrap.dedent(
            """\
            name: documented
            ref_name: chr1
            genome_length: 5000
            seed: 7
            n_fragments: 50
            frag_mean: 250
            frag_std: 50
            frag_min: 50
            frag_max: 1000
            read_length: 150
            error_rate: 0.0
            genes:
              - gene_id: g1
                strand: "+"
                transcripts:
                  - {t_id: t1, exons: [[100, 300], [500, 700]], abundance: 100}
            """
        )
    )
    args = build_parser().parse_args(["sim", "--config", str(cfg), "-o", str(tmp_path / "out")])
    assert args.func(args) == 0


# ---------------------------------------------------------------------------
# Quant defaults (all overridable args should be None before resolution)
# ---------------------------------------------------------------------------


class TestQuantDefaults:
    """Before _resolve_quant_args, overridable args are None."""

    def test_include_multimap_default_none(self):
        args = _parse_quant()
        assert args.include_multimap is None

    def test_keep_duplicates_default_none(self):
        args = _parse_quant()
        assert args.keep_duplicates is None

    def test_config_default_none(self):
        args = _parse_quant()
        assert args.config is None

    def test_scan_read_name_batch_size_default_none(self):
        args = _parse_quant()
        assert args.scan_read_name_batch_size is None

    def test_scan_bgzf_threads_default_none(self):
        args = _parse_quant()
        assert args.scan_bgzf_threads is None

    def test_scan_buffer_size_default_none(self):
        args = _parse_quant()
        assert args.scan_buffer_size is None

    def test_tsv_default_none(self):
        args = _parse_quant()
        assert args.tsv is None

    def test_emit_locus_stats_default_none(self):
        args = _parse_quant()
        assert args.emit_locus_stats is None


# ---------------------------------------------------------------------------
# Boolean optional action
# ---------------------------------------------------------------------------


class TestBooleanFlags:
    """--flag / --no-flag correctly set True / False."""

    def test_include_multimap_explicit_true(self):
        args = _parse_quant("--include-multimap")
        assert args.include_multimap is True

    def test_include_multimap_explicit_false(self):
        args = _parse_quant("--no-include-multimap")
        assert args.include_multimap is False

    def test_keep_duplicates_explicit_true(self):
        args = _parse_quant("--keep-duplicates")
        assert args.keep_duplicates is True

    def test_keep_duplicates_explicit_false(self):
        args = _parse_quant("--no-keep-duplicates")
        assert args.keep_duplicates is False


# ---------------------------------------------------------------------------
# _resolve_quant_args: CLI > YAML > defaults
# ---------------------------------------------------------------------------


class TestResolveQuant:
    """_resolve_quant_args merges CLI > YAML > defaults."""

    def test_hardcoded_defaults_applied(self):
        args = _parse_quant()
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.include_multimap is True
        assert args.sj_strand_tag == ["auto"]

    def test_cli_overrides_default(self):
        args = _parse_quant("--no-include-multimap")
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.include_multimap is False

    def test_yaml_overrides_default(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            textwrap.dedent("""\
            em_iterations: 500
            include_multimap: false
        """)
        )
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.em_iterations == 500
        assert args.include_multimap is False

    def test_quant_yaml_block_overrides_default(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            textwrap.dedent("""\
            quant:
              em_iterations: 500
        """)
        )
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.em_iterations == 500

    def test_cli_overrides_yaml(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            textwrap.dedent("""\
            em_iterations: 500
        """)
        )
        args = _parse_quant(
            "--config",
            str(cfg),
            "--em-iterations",
            "200",
        )
        _resolve_quant_args(args, _build_quant_defaults())
        # CLI wins
        assert args.em_iterations == 200

    def test_yaml_hyphens_normalised(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("em-iterations: 500\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.em_iterations == 500

    def test_yaml_unknown_key_is_refused(self, tmp_path):
        """A misspelt key must not run on the default of the key it meant, at the top level or in
        the ``quant`` block."""
        for body in ("em_iteration: 500\n", "quant:\n  em_iteration: 500\n"):
            cfg = tmp_path / "cfg.yaml"
            cfg.write_text(body)
            args = _parse_quant("--config", str(cfg))
            with pytest.raises(ValueError, match=r"\['em_iteration'\]"):
                _resolve_quant_args(args, _build_quant_defaults())

    def test_yaml_sj_strand_tag_string(self, tmp_path):
        """YAML with scalar sj_strand_tag works."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("sj_strand_tag: XS\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.sj_strand_tag == "XS"

    def test_yaml_sj_strand_tag_list(self, tmp_path):
        """YAML with list sj_strand_tag works."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("sj_strand_tag: [XS, ts]\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.sj_strand_tag == ["XS", "ts"]

    def test_yaml_scan_read_name_batch_size(self, tmp_path):
        """YAML can set the advanced scanner queue batch size."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("scan_read_name_batch_size: 128\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.scan_read_name_batch_size == 128

    def test_cli_scan_read_name_batch_size_overrides_yaml(self, tmp_path):
        """CLI qname batch size wins over YAML like other parameters."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("scan_read_name_batch_size: 128\n")
        args = _parse_quant("--config", str(cfg), "--scan-read-name-batch-size", "256")
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.scan_read_name_batch_size == 256

    def test_yaml_tsv_and_emit_locus_stats_apply_when_the_flags_are_absent(self, tmp_path):
        """A boolean set in the YAML reaches the run unless the command line names the flag."""
        from rigel.cli import _build_pipeline_config

        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("tsv: true\nemit_locus_stats: true\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.tsv is True
        assert args.emit_locus_stats is True
        assert _build_pipeline_config(args).emit_locus_stats is True

        args = _parse_quant("--config", str(cfg), "--no-tsv", "--no-emit-locus-stats")
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.tsv is False
        assert args.emit_locus_stats is False

    def test_config_yaml_rerun_keeps_tsv_and_emit_locus_stats(self, tmp_path):
        """The config.yaml a run writes reproduces its --tsv and --emit-locus-stats."""
        from rigel.cli import _write_config_yaml

        args = _parse_quant("--tsv", "--emit-locus-stats")
        _resolve_quant_args(args, _build_quant_defaults())
        written = tmp_path / "config.yaml"
        _write_config_yaml(written, args)

        rerun = build_parser().parse_args(["quant", "--config", str(written)])
        _resolve_quant_args(rerun, _build_quant_defaults())
        assert rerun.tsv is True
        assert rerun.emit_locus_stats is True


# ---------------------------------------------------------------------------
# Config round-trip: defaults → resolve → build should match PipelineConfig()
# ---------------------------------------------------------------------------


class TestConfigRoundTrip:
    """Registry-driven CLI ↔ config round-trip."""

    def test_default_config_round_trip(self):
        """No CLI/YAML overrides → config matches PipelineConfig() defaults."""
        import dataclasses
        from rigel.config import PipelineConfig
        from rigel.cli import _build_pipeline_config

        args = _parse_quant()
        _resolve_quant_args(args, _build_quant_defaults())

        result = _build_pipeline_config(args)
        ref = PipelineConfig()

        # EM fields, the seed included: an unset --seed is the config's fixed seed, never a timestamp
        for f in dataclasses.fields(ref.em):
            assert getattr(result.em, f.name) == getattr(ref.em, f.name), f.name

        # Scan fields, sj_strand_tag included: it round-trips through its transform
        for f in dataclasses.fields(ref.scan):
            assert getattr(result.scan, f.name) == getattr(ref.scan, f.name), f.name

        # Scoring: log penalties match exactly
        assert result.scoring.overhang_log_penalty == ref.scoring.overhang_log_penalty
        assert result.scoring.mismatch_log_penalty == ref.scoring.mismatch_log_penalty

    def test_param_specs_cover_all_defaults(self):
        """Every key in _build_quant_defaults matches a _ParamSpec or is CLI-only."""
        from rigel.cli import _PARAM_SPECS

        spec_dests = {s.cli_dest for s in _PARAM_SPECS}
        cli_only: set[str] = set()
        defaults = _build_quant_defaults()
        for key in defaults:
            assert key in spec_dests or key in cli_only, (
                f"Default key {key!r} not in _PARAM_SPECS or cli_only"
            )

    def test_scan_read_name_batch_size_flows_to_config(self):
        """``--scan-read-name-batch-size`` reaches ``BamScanConfig``."""
        from rigel.cli import _build_pipeline_config

        args = _parse_quant("--scan-read-name-batch-size", "256")
        _resolve_quant_args(args, _build_quant_defaults())
        cfg = _build_pipeline_config(args)
        assert cfg.scan.read_name_batch_size == 256

    def test_scan_performance_flags_flow_to_config(self):
        """The renamed scan performance flags reach ``BamScanConfig``."""
        from rigel.cli import _build_pipeline_config

        args = _parse_quant(
            "--threads",
            "8",
            "--scan-bgzf-threads",
            "2",
            "--scan-buffer-size",
            "1.5",
            "--scan-fragments-per-chunk",
            "1234",
        )
        _resolve_quant_args(args, _build_quant_defaults())
        cfg = _build_pipeline_config(args)
        assert cfg.em.n_threads == 8
        assert cfg.scan.total_threads == 8
        assert cfg.calibration.n_threads == 8
        assert cfg.scan.bgzf_threads == 2
        assert cfg.scan.buffer_size_bytes == int(1.5 * 1024**3)
        assert cfg.scan.fragments_per_chunk == 1234

    def test_scan_buffer_default_is_two_gib(self):
        """The default scan buffer cap is 2 GiB, as the manual and --help state."""
        assert BamScanConfig().buffer_size_bytes == 2 * 1024**3
