"""The index records the unobservable genome (`rigel.observable`): the runs of the reference no read pair is observed
at, written as ``unobservable.feather`` and checked at load.

Without an alignable store the mask is the FASTA's unknown runs (``N`` or ``n``); with one it is the store's zero set
at its longest read-length bin — every base no covering pair was placed back at — so a short N gap the store's pairs
span stays observable and the read length the store measured best is the one the mask reads. Every gate below holds the
index against an enumeration of the same rule, one base at a time.

The store is `_index_builder.StandInStore`, the package's public surface over zero runs given per bin, so these gates
hold the index's windowing, stitching and reference policy, not the store's measure.
"""

from __future__ import annotations

import json
import logging

import numpy as np
import pandas as pd
import pytest
from _index_builder import MINI_GTF, StandInStore, stand_in_store

import rigel.observable as OBS
import rigel.splice_blacklist as SB
from rigel.index import MANIFEST_JSON, UNOBSERVABLE_FEATHER, IndexLoadError, TranscriptIndex
from rigel.observable import unknown_bases


def _runs_of(flags) -> list[tuple[int, int]]:
    """The ``[start, end)`` runs of True, one base at a time."""
    flags = list(flags)
    runs, start = [], None
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        if not flag and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(flags)))
    return runs


def _sequence(length: int, edits: dict[int, str]) -> str:
    """``ACGT`` repeated to *length*, with each edit written over it at its position."""
    seq = list(("ACGT" * (length // 4 + 1))[:length])
    for pos, text in edits.items():
        seq[pos : pos + len(text)] = list(text)
    return "".join(seq)


#: A genome that exercises the unknown-base rule: a run at a reference's start, a single N, a lower-case run, a mixed
#: case run, an IUPAC code that is NOT unknown, a run ending at the reference's end, a reference with none, and one
#: that is unknown end to end.
GENOME = {
    "chr1": _sequence(
        2000, {0: "NNN", 100: "N", 500: "nnnn", 700: "NnN", 800: "R", 1990: "N" * 10}
    ),
    "chr2": _sequence(1000, {}),
    "chr3": "N" * 300,
}


def _write_genome(base, genome=GENOME):
    import pysam

    fasta = base / "genome.fa"
    with open(fasta, "w") as f:
        for ref, seq in genome.items():
            f.write(f">{ref}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(fasta))
    gtf = base / "test.gtf"
    gtf.write_text(MINI_GTF)
    return fasta, gtf


def _expected(runs_by_ref: dict[str, list[tuple[int, int]]]) -> list[tuple[str, int, int]]:
    return [(ref, s, e) for ref, runs in runs_by_ref.items() for s, e in runs]


def _rows(index) -> list[tuple[str, int, int]]:
    df = index.unobservable_df
    return list(
        zip(df["ref"].astype(str), df["start"].astype(int), df["end"].astype(int), strict=True)
    )


def _manifest(index_dir) -> dict:
    return json.loads((index_dir / MANIFEST_JSON).read_text())


# ── the rule for a base ───────────────────────────────────────────────────────────────────────────


def test_a_base_is_unknown_iff_it_is_N_or_n():
    seq = "ACGTNnRrYyKkMmSsWwBbDdHhVv-.Xx"
    got = unknown_bases(seq)
    assert got.tolist() == [c in "Nn" for c in seq]


# ── without a store ───────────────────────────────────────────────────────────────────────────────


def test_without_a_store_the_index_records_the_fastas_unknown_runs(tmp_path):
    fasta, gtf = _write_genome(tmp_path)
    out = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, out, write_tsv=False)
    expected = _expected({ref: _runs_of(c in "Nn" for c in seq) for ref, seq in GENOME.items()})
    assert expected == [
        ("chr1", 0, 3),
        ("chr1", 100, 101),
        ("chr1", 500, 504),
        ("chr1", 700, 703),
        ("chr1", 1990, 2000),
        ("chr3", 0, 300),
    ]
    first = TranscriptIndex.load(out)
    assert _rows(first) == expected
    # and again through a second load: the file is the record, not the build's memory
    assert _rows(TranscriptIndex.load(out)) == expected

    observable = _manifest(out)["observable"]
    masked = sum(e - s for _r, s, e in expected)
    genome = sum(len(seq) for seq in GENOME.values())
    assert observable["read_length"] is None
    assert observable["unobservable_bases"] == masked
    assert observable["observable_bases"] == genome - masked


# ── with a store ──────────────────────────────────────────────────────────────────────────────────


#: The store's measure of chr1 and chr2. At its longest bin (100) chr1's zero set holds its N runs but NOT the single N
#: at 100, a short gap the store's pairs span with both mates on known sequence, and a repeat at [400, 650) the shorter
#: bin (50) also reaches past; chr2 has a zero run across [930, 1000). chr3 is not in the store.
STORE_ZERO = {
    "chr1": {
        100: [(0, 3), (400, 650), (500, 504), (700, 703), (1990, 2000)],
        50: [(0, 3), (100, 101), (380, 700), (1990, 2000)],
    },
    "chr2": {100: [(930, 1000)], 50: [(900, 1000)]},
}


def _store(*, lengths=None, status="complete"):
    return StandInStore(
        lengths or {ref: len(GENOME[ref]) for ref in STORE_ZERO},
        bins=(50, 100),
        zero=STORE_ZERO,
        status=status,
    )


def _merged(runs):
    """Runs with overlaps and touches joined, ascending."""
    out = []
    for s, e in sorted(runs):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


@pytest.fixture
def store_build(tmp_path, monkeypatch):
    """Build with a stand-in store; the blacklist comes back empty, as a store without artifacts gives it."""
    store_dir = tmp_path / "store"
    store_dir.mkdir()
    (store_dir / "zarr.json").write_text(
        "{}"
    )  # the manifest digests the store's path as a real build's does
    blacklist = pd.DataFrame(
        {
            c: pd.Series([], dtype=t)
            for c, t in (
                ("ref", object),
                ("start", np.int32),
                ("end", np.int32),
                ("max_anchor_left", np.int32),
                ("max_anchor_right", np.int32),
                ("count", np.int64),
            )
        }
    )
    monkeypatch.setattr(SB, "load_splice_blacklist_from_zarr", lambda path, min_count=2: blacklist)

    def build(store, out):
        stand_in_store(monkeypatch, store)
        fasta, gtf = _write_genome(tmp_path)
        TranscriptIndex.build(fasta, gtf, out, write_tsv=False, alignable_zarr_path=store_dir)
        return TranscriptIndex.load(out)

    return build


def test_with_a_store_the_mask_is_its_zero_set_at_its_longest_bin(tmp_path, store_build, caplog):
    store = _store()
    with caplog.at_level(logging.WARNING, logger="rigel.observable"):
        index = store_build(store, tmp_path / "idx")
    expected = _expected(
        {
            "chr1": _merged(STORE_ZERO["chr1"][100]),
            "chr2": STORE_ZERO["chr2"][100],
            # not in the store: its N runs
            "chr3": [(0, 300)],
        }
    )
    assert _rows(index) == expected
    assert ("chr1", 100, 101) not in _rows(index)  # an N gap pairs span is observable to the store
    assert {rl for *_x, rl in store.calls} == {100}  # the longest bin, and no other
    assert "chr3" in caplog.text  # the reference the store lacks is named

    observable = _manifest(tmp_path / "idx")["observable"]
    assert observable["read_length"] == 100
    assert observable["read_length_bins"] == [50, 100]
    assert observable["references_from_fasta"] == ["chr3"]
    masked = sum(e - s for _r, s, e in expected)
    assert observable["unobservable_bases"] == masked


@pytest.mark.parametrize("window", [1, 7, 64, 333, 10_000])
def test_the_window_moves_no_run(window, tmp_path, monkeypatch):
    """Runs cross window edges at every one of these windows (chr1's [400, 650) and [1990, 2000), chr2's [930, 1000));
    the mask read a window at a time is the mask read whole."""
    stand_in_store(monkeypatch, _store())
    fasta, _gtf = _write_genome(tmp_path)
    lengths = {ref: len(seq) for ref, seq in GENOME.items()}
    whole, _ = OBS.mask_from_store("stand-in", lengths, fasta, window=max(lengths.values()))
    got, _ = OBS.mask_from_store("stand-in", lengths, fasta, window=window)
    assert got.equals(whole)


@pytest.mark.parametrize(
    "store_kw, refusal",
    [
        ({"lengths": {"chr1": 1999, "chr2": 1000}}, "another length"),
        ({"status": "partial"}, "not 'complete'"),
    ],
)
def test_a_store_for_another_assembly_or_unfinished_is_refused_before_anything_is_written(
    store_kw, refusal, tmp_path, store_build
):
    out = tmp_path / "idx"
    # an index already there, built from the same sources with no store
    fasta, gtf = _write_genome(tmp_path)
    TranscriptIndex.build(fasta, gtf, out, write_tsv=False)
    before = _rows(TranscriptIndex.load(out))
    with pytest.raises(ValueError, match=refusal):
        store_build(_store(**store_kw), out)
    # refused before the first write: the index it would have replaced still loads, unchanged
    assert _rows(TranscriptIndex.load(out)) == before


# ── at load ───────────────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "corrupt, message",
    [
        (lambda df: df.iloc[[1, 0, 2, 3, 4, 5]], "ascending and merged"),
        (lambda df: df.assign(end=np.where(df.index == 0, 101, df["end"])), "ascending and merged"),
        (lambda df: df.assign(end=np.where(df.index == 0, 100, df["end"])), "ascending and merged"),
        (
            lambda df: df.assign(end=np.where(df.index == 4, 2001, df["end"])),
            "outside its reference",
        ),
        (lambda df: df.assign(end=np.where(df.index == 0, 0, df["end"])), "empty or outside"),
        (
            lambda df: df.assign(ref=np.where(df.index == 5, "chrX", df["ref"])),
            "references the index lacks",
        ),
        (lambda df: df.rename(columns={"end": "stop"}), "columns"),
    ],
)
def test_the_load_refuses_a_mask_that_is_not_merged_runs_inside_the_references(
    corrupt, message, tmp_path
):
    fasta, gtf = _write_genome(tmp_path)
    out = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, out, write_tsv=False)
    path = out / UNOBSERVABLE_FEATHER
    corrupt(pd.read_feather(path)).reset_index(drop=True).to_feather(path)
    with pytest.raises(IndexLoadError, match=message):
        TranscriptIndex.load(out)


def test_an_index_without_the_mask_is_refused(tmp_path):
    fasta, gtf = _write_genome(tmp_path)
    out = tmp_path / "idx"
    TranscriptIndex.build(fasta, gtf, out, write_tsv=False)
    (out / UNOBSERVABLE_FEATHER).unlink()
    with pytest.raises(IndexLoadError, match=UNOBSERVABLE_FEATHER):
        TranscriptIndex.load(out)
    # and one of the format before the mask, whatever files it holds
    TranscriptIndex.build(fasta, gtf, out, write_tsv=False)
    manifest = _manifest(out)
    manifest["format_version"] = 10
    (out / MANIFEST_JSON).write_text(json.dumps(manifest))
    with pytest.raises(IndexLoadError, match="format_version=10"):
        TranscriptIndex.load(out)


# ── what a run reports of it ──────────────────────────────────────────────────────────────────────


def _per_reference(index) -> dict[str, dict[str, int]]:
    """The enumeration: each reference's masked bases, one run at a time, against its length."""
    df = index.unobservable_df
    masked = {ref: 0 for ref in index.ref_lengths}
    for ref, s, e in zip(df["ref"].astype(str), df["start"], df["end"], strict=True):
        masked[ref] += int(e - s)
    return {
        ref: {"observable_bases": int(n) - masked[ref], "unobservable_bases": masked[ref]}
        for ref, n in index.ref_lengths.items()
    }


def test_a_run_reports_the_mask_without_a_store(tmp_path):
    """`summary.json`'s observable genome: the index's own record of what the mask is, its bases per reference as an
    enumeration counts them, and no store."""
    from rigel import cli

    fasta, gtf = _write_genome(tmp_path)
    TranscriptIndex.build(fasta, gtf, tmp_path / "idx", write_tsv=False)
    index = TranscriptIndex.load(tmp_path / "idx")
    block = cli._observable_genome_summary(index)
    record = _manifest(tmp_path / "idx")["observable"]
    assert block["rule"] == record["rule"] and block["read_length"] is None
    assert block["store"] is None
    assert block["by_reference"] == _per_reference(index)
    assert (
        block["observable_bases"]
        == record["observable_bases"]
        == sum(v["observable_bases"] for v in block["by_reference"].values())
    )
    assert block["unobservable_bases"] == record["unobservable_bases"]


def test_a_run_reports_the_mask_and_its_store(tmp_path, store_build):
    """With a store, the block names the store the mask was read from: its aligner, the bins and the bin read, the
    fragment model it tiled, and the references it lacks, masked at their unknown bases."""
    from rigel import cli

    index = store_build(_store(), tmp_path / "idx")
    block = cli._observable_genome_summary(index)
    record = _manifest(tmp_path / "idx")["observable"]
    assert block["rule"] == record["rule"] and block["read_length"] == 100
    assert block["store"] == {
        "aligner": "star",
        "aligner_version": "2.7.11b",
        "read_length_bins": [50, 100],
        "frag_len_mean": 250,
        "frag_len_sd": 50,
        "error_rate": 0.01,
        "references_from_fasta": ["chr3"],
    }
    assert block["by_reference"] == _per_reference(index)
    assert block["unobservable_bases"] == record["unobservable_bases"]
