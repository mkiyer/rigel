"""The splice-junction strand tag a library is read with: what the user names is checked, what ``auto`` finds is said
and recorded, and a library whose spliced reads carry no tag is not quantified as unstranded in silence.

Junction strands come only from the tag (XS, or minimap2's ts). ``auto`` reads the first spliced reads and falls to
``none`` when none carries a tag — then no fragment trains the strand model and a stranded library is quantified as
unstranded. An unknown spec used to become ``XS,ts`` in the scanner without a word, and the resolved mode was written
nowhere, so a run could not be reproduced from its outputs.
"""

import logging

import pytest
import yaml

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario


@pytest.mark.parametrize(
    "spec", ["auto", "none", "XS", "ts", "XS,ts", "ts,XS", ("XS",), ("XS", "ts"), ()]
)
def test_every_documented_tag_spec_is_accepted(spec):
    assert BamScanConfig(sj_strand_tag=spec).sj_strand_tag == spec


@pytest.mark.parametrize("spec", ["XZ", "xs", "XS;ts", ("XS", "XZ"), ("auto", "XS")])
def test_an_unknown_tag_spec_is_refused(spec):
    with pytest.raises(ValueError, match="sj_strand_tag"):
        BamScanConfig(sj_strand_tag=spec)


@pytest.mark.parametrize("spec", ["auto", "XZ", "xs", "XS;ts", "ts,ts"])
def test_the_native_readers_refuse_a_spec_they_cannot_read(spec):
    """One layer down, the same rule: the scanner and the annotation writer refuse any spec but the ones
    they read, ``auto`` included (the caller resolves it first), rather than reading it as ``XS,ts``."""
    from rigel.native import BamAnnotationWriter, BamScanner, FragmentResolver

    for reader in (BamScanner, BamAnnotationWriter):
        with pytest.raises(ValueError, match="sj_tag_spec must be"):
            reader(FragmentResolver(), spec, True, True)


@pytest.mark.parametrize("spec", ["XS", "ts", "XS,ts", "ts,XS", "none", ""])
def test_the_native_readers_accept_every_spec_they_read(spec):
    from rigel.native import BamAnnotationWriter, BamScanner, FragmentResolver

    for reader in (BamScanner, BamAnnotationWriter):
        reader(FragmentResolver(), spec, True, True)


@pytest.fixture(scope="module")
def spliced_library(tmp_path_factory):
    sc = Scenario("sjtag", genome_length=4000, seed=7, work_dir=tmp_path_factory.mktemp("sjtag"))
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sim = ReadSimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450, read_length=100, seed=7
    )
    result = sc.build_oracle(n_fragments=300, sim_config=sim)
    yield result
    sc.cleanup()


def _run(library, spec):
    config = PipelineConfig(em=EMConfig(seed=7), scan=BamScanConfig(sj_strand_tag=spec))
    return run_pipeline(library.bam_path, library.index, config=config)


def test_auto_is_resolved_said_and_recorded(spliced_library, caplog):
    """The result carries the mode the scan read with, not ``auto``, and the log says which it was."""
    with caplog.at_level(logging.INFO, logger="rigel.pipeline"):
        pr = _run(spliced_library, "auto")
    resolved = pr.pipeline_config.scan.sj_strand_tag
    assert resolved in ("XS", "ts", "XS,ts", "none")
    assert any(
        "splice-junction strand tag" in r.getMessage() and resolved in r.getMessage()
        for r in caplog.records
    )


def test_auto_is_resolved_in_one_place_and_only_once(spliced_library):
    """``resolve_scan_config`` is where ``auto`` becomes a tag, for the pipeline, ``scan_and_buffer`` and every
    instrument that keys a scan cache by its config: it names what the detector reads, however ``auto`` is spelled,
    and hands any other config back as it is, so a resolved config resolves to itself and a named tag stands."""
    from rigel.native import detect_sj_strand_tag
    from rigel.pipeline import resolve_scan_config

    bam = spliced_library.bam_path
    detected = detect_sj_strand_tag(str(bam))
    for spelled in ("auto", ("auto",)):
        resolved = resolve_scan_config(BamScanConfig(sj_strand_tag=spelled, total_threads=1), bam)
        assert resolved == BamScanConfig(sj_strand_tag=detected, total_threads=1)
        assert resolve_scan_config(resolved, bam) is resolved
    named = BamScanConfig(sj_strand_tag="none")
    assert resolve_scan_config(named, bam) is named


def test_spliced_reads_with_no_tag_are_not_quantified_unstranded_in_silence(
    spliced_library, caplog
):
    """With no tag the strand model cannot train on the junctions, and the library reads as unstranded: say so."""
    with caplog.at_level(logging.WARNING, logger="rigel.pipeline"):
        _run(spliced_library, "none")
    warnings = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("strand tag" in m and "unstranded" in m for m in warnings), warnings


def test_config_yaml_records_the_resolved_tag(spliced_library, tmp_path):
    """The config.yaml a run writes names the tag the scan read with, so a rerun from it reads the same."""
    from rigel.cli import (
        _build_quant_defaults,
        _resolve_quant_args,
        _write_config_yaml,
        build_parser,
    )

    pr = _run(spliced_library, "auto")
    args = build_parser().parse_args(["quant", "--bam", "x.bam", "--index", "idx", "-o", "out"])
    _resolve_quant_args(args, _build_quant_defaults())
    written = tmp_path / "config.yaml"
    _write_config_yaml(written, args, pr.pipeline_config)
    recorded = yaml.safe_load(written.read_text())["quant"]["sj_strand_tag"]
    assert ",".join(recorded) == pr.pipeline_config.scan.sj_strand_tag
