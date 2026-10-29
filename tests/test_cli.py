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

import json

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


def test_index_refuses_a_negative_nascent_tolerance_as_an_error_before_writing(tmp_path):
    """A negative tolerance clustered as 0 does, and the manifest recorded it as given. `build` validates its settings
    as an `IndexConfig`, so the refusal is an error before the output directory is made."""
    from _index_builder import MINI_GTF, write_fasta

    fasta = write_fasta(tmp_path / "genome.fa", {"chr1": 2000})
    gtf = tmp_path / "test.gtf"
    gtf.write_text(MINI_GTF)
    out = tmp_path / "idx"
    args = build_parser().parse_args(
        ["index", "--fasta", str(fasta), "--gtf", str(gtf), "-o", str(out), "--no-mappability"]
        + ["--nascent-tolerance", "-5"]
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    message = str(exited.value.code)
    assert message.startswith("Error: ") and "nascent_tolerance must be >= 0; got -5" in message
    assert not out.exists()


def test_sim_unknown_scenario_key_is_refused(tmp_path):
    """``n_fragment`` for ``n_fragments`` must not simulate the default fragment count."""
    cfg = tmp_path / "scenario.yaml"
    cfg.write_text("genome_length: 5000\nn_fragment: 10\n")
    args = build_parser().parse_args(["sim", "--config", str(cfg), "-o", str(tmp_path / "out")])
    with pytest.raises(SystemExit, match=r"Error: .*\['n_fragment'\]"):
        args.func(args)
    assert not (tmp_path / "out").exists()  # checked before anything is written


_SIM_GENE = """genes:
  - gene_id: g1
    strand: "+"
    transcripts:
      - t_id: t1
        exons: [[100, 300], [500, 700]]
"""


@pytest.mark.parametrize(
    "scenario, refusal",
    [
        (_SIM_GENE + "        abundnace: 50\n", r"abundnace"),
        (_SIM_GENE.replace("    strand:", "    gene_nmae: G1\n    strand:"), r"gene_nmae"),
        ("genome_length: 5000\ngenes:\n", r"at least one gene"),
    ],
    ids=["misspelled-transcript-key", "misspelled-gene-key", "no-gene"],
)
def test_sim_refuses_a_scenario_it_would_read_wrong(tmp_path, scenario, refusal):
    """A misspelled gene or transcript key ran on its default (an abundance of 100), and an empty genes section failed
    as a TypeError; each is refused by name before the output directory is made."""
    cfg = tmp_path / "scenario.yaml"
    cfg.write_text(
        "genome_length: 5000\n" + scenario if not scenario.startswith("genome_length") else scenario
    )
    args = build_parser().parse_args(["sim", "--config", str(cfg), "-o", str(tmp_path / "out")])
    with pytest.raises(SystemExit, match=refusal):
        args.func(args)
    assert not (tmp_path / "out").exists()


def test_quant_refuses_an_unknown_config_key_with_an_error_not_a_traceback(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("em_iteration: 500\n")
    args = _parse_quant("--config", str(cfg))
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert str(exited.value.code).startswith("Error: unknown config key(s)")


@pytest.mark.parametrize(
    "flag, value",
    [("--overhang-alpha", "nan"), ("--overhang-alpha", "-0.5"), ("--overhang-alpha", "1.5")],
)
def test_quant_refuses_an_alpha_outside_zero_to_one_before_loading_the_index(tmp_path, flag, value):
    """An alpha is the probability kept per base of overhang: NaN gated almost every fragment and wrote NaN into
    summary.json, and one above 1 rewarded overhang. It is refused as an error, before the index loads."""
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"")
    index = tmp_path / "index"
    index.mkdir()
    args = build_parser().parse_args(
        [
            "quant",
            "--bam",
            str(bam),
            "--index",
            str(index),
            "-o",
            str(tmp_path / "out"),
            flag,
            value,
        ]
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert (
        str(exited.value.code) == f"Error: {flag} must be a number in [0, 1], got {float(value)!r}"
    )


def _quant_with_empty_inputs(tmp_path, *extra):
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"")
    index = tmp_path / "index"
    index.mkdir(exist_ok=True)
    return build_parser().parse_args(
        ["quant", "--bam", str(bam), "--index", str(index), "-o", str(tmp_path / "out"), *extra]
    )


@pytest.mark.parametrize("value", ["nan", "-0.5", "1", "1.5"])
def test_quant_refuses_a_pruning_threshold_outside_zero_to_one_before_loading_the_index(
    tmp_path, value
):
    """A candidate is pruned when its likelihood falls below this fraction of its fragment's best. At 1 or above, or
    NaN, the margin became 0 and a fragment kept only its tied-best candidates; below 0 nothing was pruned. 0 keeps
    every candidate, so the range is [0, 1), refused as an error before the index loads."""
    args = _quant_with_empty_inputs(tmp_path, "--pruning-min-posterior", value)
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert str(exited.value.code) == (
        f"Error: --pruning-min-posterior must be a number in [0, 1), got {float(value)!r}"
    )


@pytest.mark.parametrize("value", [False, "1e-4"])
def test_a_pruning_threshold_that_is_not_a_number_is_refused(value):
    """YAML reads ``1e-4`` as a string, which reached the scorer and failed there as a TypeError, after the scan; a
    bool is an int to Python but never a threshold (False, as 0, would pass the range check alone)."""
    from rigel.config import FragmentScoringConfig

    with pytest.raises(ValueError, match=r"--pruning-min-posterior must be a number in \[0, 1\)"):
        FragmentScoringConfig(pruning_min_posterior=value)


@pytest.mark.parametrize("value", [0, 0.0, 1e-4, 0.999])
def test_a_pruning_threshold_in_zero_to_one_is_accepted(value):
    from rigel.config import FragmentScoringConfig

    assert FragmentScoringConfig(pruning_min_posterior=value).pruning_min_posterior == value


def test_quant_refuses_an_annotated_bam_it_could_not_write_before_the_run(tmp_path):
    """The annotated BAM is written by a second pass after the quantification; a directory that does not exist used to
    fail only then."""
    args = _quant_with_empty_inputs(
        tmp_path, "--annotated-bam", str(tmp_path / "missing" / "annotated.bam")
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    message = str(exited.value.code)
    assert message.startswith("Error: --annotated-bam") and "missing" in message


def test_quant_reports_an_index_it_cannot_load_as_an_error(tmp_path):
    args = _quant_with_empty_inputs(tmp_path)
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert str(exited.value.code).startswith("Error: Index at")


def test_quant_reports_an_index_missing_a_file_as_an_error_not_a_traceback(
    tmp_path, tmp_path_factory
):
    """`rigel quant` catches the index's load refusal; a missing file used to escape it as another exception."""
    from pathlib import Path

    from _index_builder import MINI_GTF, build_test_index

    index_dir = Path(build_test_index(tmp_path_factory, MINI_GTF, name="cli_missing").index_dir)
    (index_dir / "regions.feather").unlink()
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"")
    args = build_parser().parse_args(
        ["quant", "--bam", str(bam), "--index", str(index_dir), "-o", str(tmp_path / "out")]
    )
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    message = str(exited.value.code)
    assert message.startswith(f"Error: Index at {index_dir} is missing regions.feather")
    assert "Rebuild" in message


def test_quant_refuses_a_config_that_is_not_a_mapping(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("- em_iterations\n- 500\n")
    args = _parse_quant("--config", str(cfg))
    with pytest.raises(SystemExit) as exited:
        args.func(args)
    assert str(exited.value.code).startswith("Error: ") and "mapping" in str(exited.value.code)


def test_config_yaml_records_every_path_absolute(tmp_path, monkeypatch):
    """A config.yaml reruns from any directory: the BAM, index and output were resolved, but the annotated BAM and the
    spill directory were kept as typed, so a rerun from elsewhere wrote them elsewhere."""
    import yaml

    from rigel.cli import _write_config_yaml

    monkeypatch.chdir(tmp_path)
    args = _parse_quant("--annotated-bam", "ann/annotated.bam", "--tmpdir", "spill")
    _resolve_quant_args(args, _build_quant_defaults())
    _write_config_yaml(tmp_path / "config.yaml", args)
    quant = yaml.safe_load((tmp_path / "config.yaml").read_text())["quant"]
    assert quant["annotated_bam"] == str((tmp_path / "ann/annotated.bam").resolve())
    assert quant["tmpdir"] == str((tmp_path / "spill").resolve())


def test_a_hard_scoring_gate_dumps_as_strict_json():
    """An alpha of 0 is the hard gate, whose log-penalty is −inf. The configuration summary.json records holds the
    alpha, which is finite, so the file stays JSON (``-Infinity`` is not)."""
    import math

    from rigel.config import FragmentScoringConfig, PipelineConfig

    cfg = PipelineConfig(scoring=FragmentScoringConfig(overhang_alpha=0.0))
    json.dumps(cfg.to_dict(), allow_nan=False)
    assert cfg.scoring.overhang_log_penalty == -math.inf


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

    def test_an_unknown_key_is_refused_with_the_keys_this_version_knows(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text("em_iteration: 500\n")
        args = _parse_quant("--config", str(cfg))
        with pytest.raises(ValueError, match="knows:") as refused:
            _resolve_quant_args(args, _build_quant_defaults())
        known = str(refused.value).split("knows:")[1]
        for key in [*_build_quant_defaults(), "bam_file", "index_dir", "output_dir", "tsv"]:
            assert key in known
        assert "written by" not in str(refused.value)

    def test_a_config_another_version_wrote_is_refused_naming_that_version(
        self, tmp_path, monkeypatch
    ):
        """Every 0.7.1 ``config.yaml`` carries four keys removed since and one renamed since, so rerunning one is
        refused: say why.

        The release reads a 0.7.1 file (the tree on disk says 0.7.1 until the bump). The file is the release's
        own ``config.yaml`` with 0.7.1 in its header, the four keys added and the renamed one under its 0.7.1
        name, so a change to the header the writer emits breaks this test, not the hint."""
        import rigel
        from rigel.cli import _write_config_yaml

        monkeypatch.setattr(rigel, "__version__", "0.8.0")
        written = _parse_quant()
        _resolve_quant_args(written, _build_quant_defaults())
        cfg = tmp_path / "config.yaml"
        _write_config_yaml(cfg, written)
        ours = "Generated by 0.8.0 on "
        assert ours in cfg.read_text()
        removed = (
            "  assignment_min_posterior: 0.01\n  gdna_em_llr_bias: 0.0\n"
            "  gdna_prior_mixture_bridge: true\n  sweep_n_grid_single_strand: 9\n"
        )
        cfg.write_text(
            cfg.read_text()
            .replace(ours, "Generated by 0.7.1 on ")
            .replace("\nquant:\n", f"\nquant:\n{removed}")
            .replace("\n  intron_overhang_tolerance:", "\n  splicing_anchor_tolerance:")
        )
        args = _parse_quant("--config", str(cfg))
        with pytest.raises(ValueError, match="written by rigel 0.7.1") as refused:
            _resolve_quant_args(args, _build_quant_defaults())
        for key in (
            "assignment_min_posterior",
            "gdna_em_llr_bias",
            "gdna_prior_mixture_bridge",
            "sweep_n_grid_single_strand",
            "splicing_anchor_tolerance",
        ):
            assert key in str(refused.value).split("knows:")[0]

    def test_a_config_this_version_wrote_is_refused_with_no_version_hint(self, tmp_path):
        """A key misspelt into this version's own ``config.yaml`` is a typo, not a removed key."""
        from rigel.cli import _write_config_yaml

        written = _parse_quant()
        _resolve_quant_args(written, _build_quant_defaults())
        cfg = tmp_path / "config.yaml"
        _write_config_yaml(cfg, written)
        cfg.write_text(cfg.read_text().replace("\nquant:\n", "\nquant:\n  em_iteration: 500\n"))
        args = _parse_quant("--config", str(cfg))
        with pytest.raises(ValueError, match=r"\['em_iteration'\]") as refused:
            _resolve_quant_args(args, _build_quant_defaults())
        assert "written by" not in str(refused.value)

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

        # Scoring: the log penalty matches exactly
        assert result.scoring.overhang_log_penalty == ref.scoring.overhang_log_penalty

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
