"""CLI parser tests."""

import textwrap
from pathlib import Path

import pytest

from hulkrna.cli import build_parser, _resolve_quant_args


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

    def test_em_prior_alpha_default_none(self):
        args = _parse_quant()
        assert args.em_prior_alpha is None

    def test_config_default_none(self):
        args = _parse_quant()
        assert args.config is None

    def test_nrna_frac_kappa_min_default_none(self):
        args = _parse_quant()
        assert args.nrna_frac_kappa_min is None


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
# Advanced parameters
# ---------------------------------------------------------------------------


class TestAdvancedParams:
    """Advanced parameters are parsed correctly."""

    def test_tss_window(self):
        args = _parse_quant("--tss-window", "500")
        assert args.tss_window == 500

    def test_nrna_frac_mom_min_evidence_global(self):
        args = _parse_quant("--nrna-frac-mom-min-evidence-global", "100")
        assert args.nrna_frac_mom_min_evidence_global == 100.0

    def test_nrna_frac_kappa_fallback(self):
        args = _parse_quant("--nrna-frac-kappa-fallback", "10")
        assert args.nrna_frac_kappa_fallback == 10.0

    def test_nrna_frac_kappa_min_obs(self):
        args = _parse_quant("--nrna-frac-kappa-min-obs", "50")
        assert args.nrna_frac_kappa_min_obs == 50


# ---------------------------------------------------------------------------
# _resolve_quant_args: CLI > YAML > defaults
# ---------------------------------------------------------------------------


_SMALL_DEFAULTS = {
    "seed": None,
    "em_prior_alpha": 0.01,
    "em_prior_gamma": 1.0,
    "include_multimap": True,
    "keep_duplicates": False,
    "no_tsv": False,
    "tss_window": 200,
    "em_iterations": 1000,
    "em_convergence_delta": 1e-6,
    "confidence_threshold": 0.95,
    "sj_strand_tag": ["auto"],
    "annotated_bam": None,
    "overhang_alpha": 0.01,
    "mismatch_alpha": 0.1,
    "gdna_splice_penalty_unannot": 0.01,
    "nrna_frac_kappa_global": None,
    "nrna_frac_kappa_locus": None,
    "nrna_frac_kappa_tss": None,
    "nrna_frac_mom_min_evidence_global": 50.0,
    "nrna_frac_mom_min_evidence_locus": 30.0,
    "nrna_frac_mom_min_evidence_tss": 20.0,
    "nrna_frac_kappa_min": 2.0,
    "nrna_frac_kappa_max": 200.0,
    "nrna_frac_kappa_fallback": 5.0,
    "nrna_frac_kappa_min_obs": 20,
}


class TestResolveQuant:
    """_resolve_quant_args merges CLI > YAML > defaults."""

    def test_hardcoded_defaults_applied(self):
        args = _parse_quant()
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.include_multimap is True
        assert args.em_prior_alpha == 0.01
        assert args.tss_window == 200
        assert args.nrna_frac_kappa_min == 2.0
        assert args.sj_strand_tag == ["auto"]

    def test_cli_overrides_default(self):
        args = _parse_quant("--no-include-multimap", "--tss-window", "500")
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.include_multimap is False
        assert args.tss_window == 500

    def test_yaml_overrides_default(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(textwrap.dedent("""\
            em_prior_alpha: 0.05
            tss_window: 300
            include_multimap: false
        """))
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.em_prior_alpha == 0.05
        assert args.tss_window == 300
        assert args.include_multimap is False

    def test_cli_overrides_yaml(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(textwrap.dedent("""\
            em_prior_alpha: 0.05
            tss_window: 300
        """))
        args = _parse_quant(
            "--config", str(cfg),
            "--em-prior-alpha", "0.1",
        )
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        # CLI wins
        assert args.em_prior_alpha == 0.1
        # YAML wins over default
        assert args.tss_window == 300

    def test_yaml_hyphens_normalised(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("nrna-frac-kappa-min: 4.0\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.nrna_frac_kappa_min == 4.0

    def test_yaml_unknown_keys_logged(self, tmp_path, caplog):
        """Unknown YAML keys are logged as warnings, not raised."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("bogus_key: 99\n")
        args = _parse_quant("--config", str(cfg))
        import logging
        with caplog.at_level(logging.WARNING):
            _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert "bogus_key" in caplog.text

    def test_yaml_sj_strand_tag_string(self, tmp_path):
        """YAML with scalar sj_strand_tag works."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("sj_strand_tag: XS\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.sj_strand_tag == "XS"

    def test_yaml_sj_strand_tag_list(self, tmp_path):
        """YAML with list sj_strand_tag works."""
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("sj_strand_tag: [XS, ts]\n")
        args = _parse_quant("--config", str(cfg))
        _resolve_quant_args(args, _SMALL_DEFAULTS)
        assert args.sj_strand_tag == ["XS", "ts"]
