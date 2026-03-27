"""CLI parser tests."""

import textwrap


from rigel.cli import build_parser, _resolve_quant_args, _build_quant_defaults


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

    def test_prior_pseudocount_default_none(self):
        args = _parse_quant()
        assert args.prior_pseudocount is None

    def test_config_default_none(self):
        args = _parse_quant()
        assert args.config is None


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
        assert args.prior_pseudocount == 1.0
        assert args.sj_strand_tag == ["auto"]

    def test_cli_overrides_default(self):
        args = _parse_quant("--no-include-multimap")
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.include_multimap is False

    def test_yaml_overrides_default(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(textwrap.dedent("""\
            prior_pseudocount: 2.0
            include_multimap: false
        """))
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.prior_pseudocount == 2.0
        assert args.include_multimap is False

    def test_cli_overrides_yaml(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(textwrap.dedent("""\
            prior_pseudocount: 2.0
        """))
        args = _parse_quant(
            "--config", str(cfg),
            "--prior-pseudocount", "3.0",
        )
        _resolve_quant_args(args, _build_quant_defaults())
        # CLI wins
        assert args.prior_pseudocount == 3.0

    def test_yaml_hyphens_normalised(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("prior-pseudocount: 2.0\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _build_quant_defaults())
        assert args.prior_pseudocount == 2.0

    def test_yaml_unknown_keys_logged(self, tmp_path, caplog):
        """Unknown YAML keys are logged as warnings, not raised."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("bogus_key: 99\n")
        args = _parse_quant("--config", str(cfg))
        import logging
        with caplog.at_level(logging.WARNING):
            _resolve_quant_args(args, _build_quant_defaults())
        assert "bogus_key" in caplog.text

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

        result = _build_pipeline_config(args, seed=42, sj_strand_tag="auto")
        ref = PipelineConfig()

        # EM fields (except overridden seed)
        for f in dataclasses.fields(ref.em):
            if f.name == "seed":
                assert result.em.seed == 42
                continue
            assert getattr(result.em, f.name) == getattr(ref.em, f.name), f.name

        # Scan fields (except overridden sj_strand_tag)
        for f in dataclasses.fields(ref.scan):
            if f.name == "sj_strand_tag":
                assert result.scan.sj_strand_tag == "auto"
                continue
            assert getattr(result.scan, f.name) == getattr(ref.scan, f.name), f.name

        # Scoring: log penalties match exactly
        assert result.scoring.overhang_log_penalty == ref.scoring.overhang_log_penalty
        assert result.scoring.mismatch_log_penalty == ref.scoring.mismatch_log_penalty
        # gdna_splice_penalties: result has explicit dict, ref has None
        assert result.scoring.gdna_splice_penalties is not None

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
