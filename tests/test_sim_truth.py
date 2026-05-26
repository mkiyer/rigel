"""Tests for shared simulator truth and manifest helpers."""

import json
from dataclasses import dataclass

from rigel.sim.manifest import (
    condition_dir_name,
    condition_manifest_map,
    gdna_label_for_rate,
    load_manifest,
    write_manifest,
)
from rigel.sim.truth import count_origins_from_fastq, parse_origin


def test_parse_mrna_fastq_origin():
    origin = parse_origin("@tx1:10-250:f:7/1")

    assert origin.kind == "mrna"
    assert origin.transcript_id == "tx1"
    assert origin.start == 10
    assert origin.end == 250
    assert origin.strand == "f"
    assert origin.index == 7


def test_parse_nrna_bam_origin():
    origin = parse_origin("nrna_tx2:4-300:r:11")

    assert origin.kind == "nrna"
    assert origin.transcript_id == "tx2"
    assert origin.start == 4
    assert origin.end == 300
    assert origin.strand == "r"
    assert origin.index == 11


def test_parse_gdna_origin_without_ref():
    origin = parse_origin("gdna:100-390:f:3/2")

    assert origin.kind == "gdna"
    assert origin.transcript_id is None
    assert origin.ref is None
    assert origin.start == 100
    assert origin.end == 390
    assert origin.strand == "f"
    assert origin.index == 3


def test_parse_gdna_origin_with_ref():
    origin = parse_origin("gdna:chrSynthetic:100-390:r:3")

    assert origin.kind == "gdna"
    assert origin.ref == "chrSynthetic"
    assert origin.start == 100
    assert origin.end == 390
    assert origin.strand == "r"
    assert origin.index == 3


def test_count_origins_from_fastq(tmp_path):
    fastq = tmp_path / "reads.fq"
    fastq.write_text(
        "@tx1:0-100:f:0/1\nACGT\n+\nIIII\n"
        "@nrna_tx1:0-100:f:1/1\nACGT\n+\nIIII\n"
        "@gdna:chr1:0-100:r:2/1\nACGT\n+\nIIII\n"
    )

    counts = count_origins_from_fastq(fastq)

    assert counts["mrna"] == 1
    assert counts["nrna"] == 1
    assert counts["gdna"] == 1


def test_condition_dir_name_and_gdna_labels():
    assert condition_dir_name("low", 0.9, "rand") == "gdna_low_ss_0.90_nrna_rand"
    assert (
        condition_dir_name("low", 0.9, "rand", "on")
        == "gdna_low_ss_0.90_nrna_rand_capture_on"
    )
    assert gdna_label_for_rate(0, None, 0) == "r0"
    assert gdna_label_for_rate(0.3, None, 1) == "r0.3"
    assert gdna_label_for_rate(0.3, ["none", "low"], 1) == "low"


def test_manifest_round_trip(tmp_path):
    @dataclass
    class SimSection:
        frag_mean: int = 250

    class Config:
        simulation = SimSection()
        gdna = {"frag_mean": 350}
        nrna = {"modes": ["none"]}
        abundance = {"kind": "lognormal"}

    conditions = [{"name": "gdna_none_ss_1.00_nrna_none", "n_rna": 10}]

    path = write_manifest(tmp_path, Config(), conditions)
    loaded = load_manifest(tmp_path)

    assert path == tmp_path / "manifest.json"
    assert json.loads(path.read_text()) == loaded
    assert loaded["simulation"]["frag_mean"] == 250
    assert condition_manifest_map(loaded) == {conditions[0]["name"]: conditions[0]}