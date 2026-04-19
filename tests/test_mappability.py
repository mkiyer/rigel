"""
Tests for ``rigel.mappability`` and the index ``mappable_effective_length``
column / manifest plumbing.

The tests deliberately avoid the alignable Zarr opener — that path is
covered by manual integration runs against the real GRCh38/STAR store
in ``scripts/debug/alignable_landscape_probe.py``.  Here we exercise:

* :func:`rigel.mappability.uniform_region_exposures` — the
  ``--no-mappability`` fallback (assume ``f(p) = 1``).
* The index build with no Zarr writes a sensible regions.feather and
  manifest, with ``mappable_effective_length == length`` per region.
* The manifest JSON records the format version, the rigel version,
  and ``"mappability": null`` in the no-Zarr case.
* :func:`rigel.index.load_manifest` round-trips and returns ``None``
  when the file is absent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from rigel.index import (
    INDEX_FORMAT_VERSION,
    MANIFEST_JSON,
    REGIONS_FEATHER,
    TranscriptIndex,
    load_manifest,
)
from rigel.mappability import uniform_region_exposures


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _mini_inputs(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "ann.gtf"

    with open(fasta, "w") as fh:
        fh.write(">chr1\n")
        fh.write("A" * 2000 + "\n")
    pysam.faidx(str(fasta))

    with open(gtf, "w") as fh:
        fh.write(
            'chr1\tsrc\texon\t1\t300\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tsrc\texon\t401\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            'chr1\tsrc\texon\t1\t1000\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
        )
    return fasta, gtf


# ----------------------------------------------------------------------
# uniform_region_exposures
# ----------------------------------------------------------------------


class TestUniformExposures:
    def test_matches_length(self) -> None:
        df = pd.DataFrame({
            "ref": ["chr1", "chr1", "chr2"],
            "start": [0, 100, 50],
            "end": [50, 250, 200],
        })
        out = uniform_region_exposures(df)
        assert out.dtype == np.float32
        np.testing.assert_array_equal(out, np.array([50, 150, 150], dtype=np.float32))

    def test_empty(self) -> None:
        df = pd.DataFrame({"ref": [], "start": [], "end": []})
        out = uniform_region_exposures(df)
        assert out.shape == (0,)
        assert out.dtype == np.float32

    def test_negative_or_zero_clamped(self) -> None:
        df = pd.DataFrame({
            "ref": ["chr1", "chr1"],
            "start": [10, 100],
            "end": [10, 50],   # zero-length and inverted
        })
        out = uniform_region_exposures(df)
        np.testing.assert_array_equal(out, np.array([0, 0], dtype=np.float32))


# ----------------------------------------------------------------------
# Index build with --no-mappability path
# ----------------------------------------------------------------------


class TestIndexBuildNoMappability:
    def test_regions_have_mappable_effective_length_column(self, tmp_path: Path) -> None:
        fasta, gtf = _mini_inputs(tmp_path)
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(
            fasta_file=str(fasta), gtf_file=str(gtf), output_dir=str(out_dir),
        )
        regions = pd.read_feather(out_dir / REGIONS_FEATHER)
        assert "mappable_effective_length" in regions.columns
        assert regions["mappable_effective_length"].dtype == np.float32
        # No Zarr -> exposure equals length exactly.
        np.testing.assert_array_equal(
            regions["mappable_effective_length"].to_numpy(),
            regions["length"].astype(np.float32).to_numpy(),
        )

    def test_manifest_written(self, tmp_path: Path) -> None:
        fasta, gtf = _mini_inputs(tmp_path)
        out_dir = tmp_path / "idx"
        TranscriptIndex.build(
            fasta_file=str(fasta), gtf_file=str(gtf), output_dir=str(out_dir),
        )
        path = out_dir / MANIFEST_JSON
        assert path.exists()
        manifest = load_manifest(out_dir)
        assert manifest is not None
        assert manifest["format_version"] == INDEX_FORMAT_VERSION
        assert "rigel_version" in manifest
        assert manifest["mappability"] is None

    def test_load_manifest_returns_none_when_absent(self, tmp_path: Path) -> None:
        assert load_manifest(tmp_path) is None


# ----------------------------------------------------------------------
# Zarr-driven exposures (skipped unless alignable + a real store are
# available).  Smoke-tests the documented API path end-to-end.
# ----------------------------------------------------------------------


_ZARR_PATH = "/scratch/mkiyer_root/mkiyer0/shared_data/alignable/alignable_grch38_star.zarr.zip"


@pytest.mark.skipif(
    not Path(_ZARR_PATH).exists(),
    reason="alignable Zarr store not available on this host",
)
def test_compute_region_exposures_smoke(tmp_path: Path) -> None:
    pytest.importorskip("alignable")

    from rigel.mappability import compute_region_exposures

    df = pd.DataFrame({
        "ref": ["chr22", "chr22"],
        "start": [16_000_000, 20_000_000],
        "end":   [16_001_000, 20_001_000],
        "length": [1000, 1000],
    })
    exposures, prov = compute_region_exposures(_ZARR_PATH, df, read_length=100)
    assert exposures.shape == (2,)
    assert exposures.dtype == np.float32
    assert (exposures >= 0).all()
    assert (exposures <= 1000).all()
    assert prov.read_length == 100
    assert prov.aligner == "star"
