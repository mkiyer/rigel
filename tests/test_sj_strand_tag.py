"""The splice-junction strand tag a library is read with: what the user names is checked, what the tags gave is
counted and said, and a library whose spliced reads carry no tag is not quantified as unstranded in silence.

Junction strands come only from the tag (XS, or minimap2's ts). ``auto`` reads both on every spliced record, so a
library whose first spliced reads lack the tag still has its strand: choosing one mode from the first spliced reads
read none of a late-tagged library's tags. An unknown spec used to become ``XS,ts`` in the scanner without a word.
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


@pytest.mark.parametrize("spec", ["XZ", "xs", "XS;ts", "ts,ts"])
def test_the_native_readers_refuse_a_spec_they_cannot_read(spec):
    """One layer down, the same rule: the scanner and the annotation writer refuse any spec but the ones
    they read, rather than reading it as ``XS,ts``."""
    from rigel.native import BamAnnotationWriter, BamScanner, FragmentResolver

    for reader in (BamScanner, BamAnnotationWriter):
        with pytest.raises(ValueError, match="sj_tag_spec must be"):
            reader(FragmentResolver(), spec, True, True)


@pytest.mark.parametrize("spec", ["auto", "XS", "ts", "XS,ts", "ts,XS", "none", ""])
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


def test_auto_is_kept_and_what_the_tags_gave_is_said(spliced_library, caplog):
    """``auto`` is a mode the scan reads with, so the result keeps it; the log says what the tags gave, record by
    record, and every spliced record is counted once."""
    with caplog.at_level(logging.INFO, logger="rigel.pipeline"):
        pr = _run(spliced_library, "auto")
    assert pr.pipeline_config.scan.sj_strand_tag == "auto"
    st = pr.stats
    assert st.n_sj_tag_tagged > 0 and st.n_sj_tag_untagged == 0 and st.n_sj_tag_conflicting == 0
    said = [
        r.getMessage()
        for r in caplog.records
        if "splice-junction strand tag auto" in r.getMessage()
    ]
    assert said and f"{st.n_sj_tag_tagged} with a strand" in said[0], said


def test_spliced_reads_with_no_tag_are_not_quantified_unstranded_in_silence(
    spliced_library, caplog
):
    """With no tag the strand model cannot train on the junctions, and the library reads as unstranded: say so."""
    with caplog.at_level(logging.WARNING, logger="rigel.pipeline"):
        _run(spliced_library, "none")
    warnings = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("strand tag" in m and "unstranded" in m for m in warnings), warnings


def test_config_yaml_records_the_tag_spec(spliced_library, tmp_path):
    """The config.yaml a run writes names the tag spec the scan read with, so a rerun from it reads the same."""
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


def _strip_tags_from_first_spliced(src, dst, n):
    """``src`` rewritten in its own record order with XS and ts removed from its first ``n`` spliced records."""
    import pysam

    seen = 0
    with (
        pysam.AlignmentFile(str(src), "rb") as inp,
        pysam.AlignmentFile(str(dst), "wb", template=inp) as out,
    ):
        for r in inp:
            if seen < n and any(op == 3 for op, _ in (r.cigartuples or ())):
                seen += 1
                r.set_tag("XS", None)
                r.set_tag("ts", None)
            out.write(r)
    return seen


@pytest.fixture(scope="module")
def late_tagged_library(tmp_path_factory):
    """A stranded library whose first 1,200 spliced records carry no strand tag and whose later ones carry XS."""
    work = tmp_path_factory.mktemp("late_tags")
    sc = Scenario("late_tags", genome_length=20_000, seed=11, work_dir=work)
    sc.add_gene(
        "g1",
        "+",
        [
            {
                "t_id": "t1",
                "exons": [(1000, 1300), (2000, 2300), (3000, 3300), (4000, 4300)],
                "abundance": 100,
            }
        ],
    )
    sim = ReadSimConfig(
        frag_mean=250,
        frag_std=40,
        frag_min=120,
        frag_max=450,
        read_length=100,
        strand_specificity=0.95,
        seed=11,
    )
    result = sc.build_oracle(n_fragments=6000, sim_config=sim)
    late = work / "late_tagged.bam"
    stripped = _strip_tags_from_first_spliced(result.bam_path, late, 1200)
    assert stripped == 1200, stripped
    yield late, result.index
    sc.cleanup()


def test_auto_reads_the_tag_on_every_record_not_from_the_first_spliced_reads(late_tagged_library):
    """``auto`` reads each record's own tag. Choosing one mode from the first spliced reads finds no tag in a library
    whose first 1,200 spliced records lack one, and then reads none of the later, tagged records: a stranded library
    trains no strand observation and is quantified as unstranded."""
    bam, index = late_tagged_library
    pr = run_pipeline(
        str(bam),
        index,
        config=PipelineConfig(em=EMConfig(seed=7), scan=BamScanConfig(sj_strand_tag="auto")),
    )
    sm = pr.strand_models.exonic_spliced
    assert sm.n_observations > 0, (
        "the tagged spliced records after the first 1,200 trained no strand observation"
    )
    assert sm.strand_specificity > 0.8
