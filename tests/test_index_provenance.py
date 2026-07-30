"""An index must be rebuildable from its own manifest — sources, their hashes, and the build flags.

    TODO item 1   ·   Ledger: "An index cannot be rebuilt from its own manifest" (2026-07-30)

⛔ WHY THIS EXISTS. `manifest.json` recorded `format_version` and `rigel_version` and nothing else — not
the FASTA, not the GTF, not the flags. Rebuilding the human index meant *inferring* the source by matching
node counts, and `--collapse-duplicate-transcripts` was discovered from a build failure. An artifact that
cannot be reproduced from its own provenance is one nobody can safely re-derive later.

⚠ **This is not `CARRY_FORWARD.md` §3 trap 25.** That trap forbids storing a hash of an artifact *beside*
that artifact, because the two can drift apart and the stale hash then verifies clean. These hashes are of
**external inputs the index cannot recompute from itself** — the genome and the annotation are not in the
index. Provenance, not a cache key. `partition_hash` and `graph_hash` remain computed on demand.

The teeth are in three places, and each is a different failure this catches:

* the expected flag set is read off `inspect.signature(TranscriptIndex.build)`, never written out here, so
  a new build parameter that does not reach the manifest fails this file rather than silently escaping it;
* the digests are re-derived by a **different algorithm** (whole-file read here, streamed chunks in the
  implementation) — `CARRY_FORWARD.md` §3 trap 1;
* a one-byte edit to the GTF must move its recorded digest, which is what separates hashing the *content*
  from recording the *path*.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import pysam
import pytest

from rigel import index as index_module
from rigel.index import MANIFEST_JSON, INDEX_FORMAT_VERSION, TranscriptIndex, source_record

#: `build()` parameters that name a source file rather than tune the build. They are recorded under
#: `sources` with a content digest, so they must NOT also appear under `build_flags`.
SOURCE_PARAMETERS = frozenset({"fasta_file", "gtf_file", "output_dir", "alignable_zarr_path"})


def sha256_of(path: Path) -> str:
    """Digest by reading the whole file at once — deliberately NOT the implementation's chunked loop."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expected_flag_names() -> set[str]:
    """The build parameters that are neither a source nor an output, read off the signature itself."""
    parameters = inspect.signature(TranscriptIndex.build).parameters
    return {name for name in parameters if name not in SOURCE_PARAMETERS}


@pytest.fixture
def sources(tmp_path: Path) -> tuple[Path, Path]:
    """A minimal FASTA + GTF pair that `TranscriptIndex.build` accepts."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "A" * 2000 + "\n")
    pysam.faidx(str(fasta))

    gtf = tmp_path / "ann.gtf"
    gtf.write_text(
        'chr1\tsrc\texon\t1\t300\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t401\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\texon\t1\t1000\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
    )
    return fasta, gtf


def build_and_read_manifest(fasta: Path, gtf: Path, out_dir: Path, **kwargs) -> dict:
    TranscriptIndex.build(str(fasta), str(gtf), str(out_dir), write_tsv=False, **kwargs)
    return json.loads((out_dir / MANIFEST_JSON).read_text())


class TestSources:
    def test_the_fasta_and_the_gtf_are_recorded_with_their_content_digests(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        fasta, gtf = sources
        manifest = build_and_read_manifest(fasta, gtf, tmp_path / "idx")

        assert manifest["sources"]["fasta"]["path"] == str(fasta.resolve())
        assert manifest["sources"]["fasta"]["bytes"] == fasta.stat().st_size
        assert manifest["sources"]["fasta"]["sha256"] == sha256_of(fasta)

        assert manifest["sources"]["gtf"]["path"] == str(gtf.resolve())
        assert manifest["sources"]["gtf"]["bytes"] == gtf.stat().st_size
        assert manifest["sources"]["gtf"]["sha256"] == sha256_of(gtf)

    def test_an_absent_alignable_store_is_recorded_as_null_not_omitted(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        """`null` says "no store was supplied"; a missing key says "this manifest predates the field"."""
        fasta, gtf = sources
        manifest = build_and_read_manifest(fasta, gtf, tmp_path / "idx")
        assert "alignable_zarr" in manifest["sources"]
        assert manifest["sources"]["alignable_zarr"] is None

    def test_a_one_byte_edit_to_the_gtf_moves_its_recorded_digest(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        """The digest is of the CONTENT. Recording the path alone would pass every other test here."""
        fasta, gtf = sources
        first = build_and_read_manifest(fasta, gtf, tmp_path / "a")

        gtf.write_text(gtf.read_text().replace("\t401\t700\t", "\t402\t700\t"))
        second = build_and_read_manifest(fasta, gtf, tmp_path / "b")

        assert first["sources"]["gtf"]["path"] == second["sources"]["gtf"]["path"]
        assert first["sources"]["gtf"]["sha256"] != second["sources"]["gtf"]["sha256"]
        assert second["sources"]["gtf"]["sha256"] == sha256_of(gtf)

    def test_the_same_sources_record_the_same_digests_into_two_output_dirs(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        fasta, gtf = sources
        first = build_and_read_manifest(fasta, gtf, tmp_path / "a")
        second = build_and_read_manifest(fasta, gtf, tmp_path / "b")
        assert first["sources"] == second["sources"]


class TestTheDigestItself:
    """⛔ Added after perturbation: the build-level tests above could NOT see either of these.

    Every fixture here is ~2 KB against a 1 MB read chunk, and every fixture path is already absolute —
    so a digest loop that read only the first chunk, and a record that stored the raw input string,
    both passed all eight tests. On the real 1.6 GB GTF the first would hash one megabyte of a 1.6 GB
    file, and the second would record `../refs/genes.gtf`.
    """

    def test_the_digest_covers_a_file_larger_than_one_read_chunk(self, tmp_path: Path) -> None:
        """Two files with an IDENTICAL first chunk must still get different digests."""
        head = b"A" * (index_module._DIGEST_CHUNK_BYTES + 100)
        first = tmp_path / "first.bin"
        second = tmp_path / "second.bin"
        first.write_bytes(head)
        second.write_bytes(head[:-1] + b"B")

        assert (
            first.read_bytes()[: index_module._DIGEST_CHUNK_BYTES]
            == (second.read_bytes()[: index_module._DIGEST_CHUNK_BYTES])
        )
        assert source_record(first)["sha256"] != source_record(second)["sha256"]
        assert source_record(first)["sha256"] == sha256_of(first)

    def test_the_read_chunk_size_cannot_change_the_digest(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """It is an I/O buffer. If the answer moves with it, it is not one."""
        payload = tmp_path / "payload.bin"
        payload.write_bytes(bytes(range(256)) * 5000)
        at_default = source_record(payload)["sha256"]

        monkeypatch.setattr(index_module, "_DIGEST_CHUNK_BYTES", 7)
        assert source_record(payload)["sha256"] == at_default == sha256_of(payload)

    def test_a_relative_source_path_is_recorded_as_an_absolute_one(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A path recorded relative to whatever directory the build ran in points nowhere later."""
        target = tmp_path / "nested" / "source.txt"
        target.parent.mkdir()
        target.write_text("provenance")

        monkeypatch.chdir(tmp_path)
        recorded = source_record(Path("nested") / "source.txt")["path"]
        assert Path(recorded).is_absolute()
        assert Path(recorded) == target.resolve()

    def test_a_directory_source_is_digested_over_its_whole_tree(self, tmp_path: Path) -> None:
        """The alignable store may be an unpackaged `.zarr` directory; nothing else reaches that branch."""
        store = tmp_path / "store.zarr"
        (store / "deep").mkdir(parents=True)
        (store / "top.bin").write_bytes(b"top")
        (store / "deep" / "leaf.bin").write_bytes(b"leaf")

        before = source_record(store)
        assert before["bytes"] == len(b"top") + len(b"leaf")

        (store / "deep" / "leaf.bin").write_bytes(b"leaX")
        after = source_record(store)
        assert after["sha256"] != before["sha256"]

    def test_a_directory_digest_depends_on_the_names_not_only_the_contents(
        self, tmp_path: Path
    ) -> None:
        """Two stores holding the same bytes under different names are different stores."""
        left = tmp_path / "left"
        right = tmp_path / "right"
        left.mkdir()
        right.mkdir()
        (left / "a.bin").write_bytes(b"payload")
        (right / "b.bin").write_bytes(b"payload")
        assert source_record(left)["sha256"] != source_record(right)["sha256"]


class TestBuildFlags:
    def test_every_build_parameter_that_is_not_a_source_reaches_the_manifest(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        """⭐ The expected set is read off the signature, so a NEW parameter fails here by construction.

        A hand-written list in this file would drift the moment someone adds a build knob — which is the
        exact way the manifest got thin in the first place.
        """
        fasta, gtf = sources
        manifest = build_and_read_manifest(fasta, gtf, tmp_path / "idx")
        assert set(manifest["build_flags"]) == expected_flag_names()

    def test_the_recorded_flags_are_the_values_the_build_actually_used(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        fasta, gtf = sources
        manifest = build_and_read_manifest(
            fasta,
            gtf,
            tmp_path / "idx",
            collapse_duplicate_transcripts=True,
            gtf_parse_mode="warn-skip",
            nrna_tolerance=37,
            feather_compression="zstd",
        )
        flags = manifest["build_flags"]
        assert flags["collapse_duplicate_transcripts"] is True
        assert flags["gtf_parse_mode"] == "warn-skip"
        assert flags["nrna_tolerance"] == 37
        assert flags["feather_compression"] == "zstd"
        assert flags["write_tsv"] is False

    def test_defaults_are_recorded_explicitly_rather_than_left_out(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        """A rebuilder must not have to know this rigel version's defaults to reproduce the artifact."""
        fasta, gtf = sources
        manifest = build_and_read_manifest(fasta, gtf, tmp_path / "idx")
        flags = manifest["build_flags"]
        assert flags["collapse_duplicate_transcripts"] is False
        assert flags["gtf_parse_mode"] == "strict"
        assert flags["splice_blacklist_min_count"] == 2


class TestExistingKeysSurvive:
    def test_format_version_and_rigel_version_are_still_recorded(
        self, sources: tuple[Path, Path], tmp_path: Path
    ) -> None:
        """The version gate at load time reads these; provenance is additive to it, not a replacement."""
        fasta, gtf = sources
        manifest = build_and_read_manifest(fasta, gtf, tmp_path / "idx")
        assert manifest["format_version"] == INDEX_FORMAT_VERSION
        assert isinstance(manifest["rigel_version"], str)
        assert manifest["rigel_version"] != ""
