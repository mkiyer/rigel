"""The synthetic index builder every structural test starts from.

`build_test_index` turns a GTF string into a loaded `TranscriptIndex`, writing a FASTA of known sequence for
the requested references beside it (`write_fasta`, which every test that builds an index writes its genome with);
`rebuild_with_splice_blacklist` gives an existing index a splice blacklist the way a real build does. It lives in its own module rather than in `conftest.py` because it is a
plain builder, not a fixture: a test that needs a bespoke index imports and calls it, and importing a
conftest by module name is ambiguous as soon as a sub-directory has one of its own.
"""

from __future__ import annotations

import textwrap
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from rigel.index import REF_LENGTHS_FEATHER, TranscriptIndex, load_manifest, load_ref_lengths


# ---------------------------------------------------------------------------
# Two genes on chr1, one on + and one on - strand.
# Gene A (g1): two transcripts (t1 has 3 exons, t2 has 2 exons)
# Gene B (g2): one transcript (t3 has 2 exons), negative strand
#
# Coordinates (1-based inclusive as in GTF):
#   t1: exons at 100-200, 300-400, 500-600
#   t2: exons at 100-200, 500-600
#   t3: exons at 1000-1100, 1200-1300 (neg strand)

MINI_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t1000\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
    chr1\ttest\texon\t1200\t1300\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
""")

# Expected 0-based half-open coordinates after GTF parsing:
#   t1: exons (99,200), (299,400), (499,600)
#   t2: exons (99,200), (499,600)
#   t3: exons (999,1100), (1199,1300)


def write_fasta(path: str | Path, refs: Mapping[str, int]) -> Path:
    """Write an indexed FASTA of KNOWN sequence: each reference is ``ACGT`` repeated to its length.

    A structural test's genome holds no unknown base, so an index built from it observes every base and its
    unobservable set is empty, as a real genome's is outside its gaps. An all-``N`` genome is unobservable end to end.
    Returns *path*.
    """
    import pysam

    path = Path(path)
    with open(path, "w") as f:
        for ref_name, ref_len in refs.items():
            f.write(f">{ref_name}\n")
            seq = ("ACGT" * (ref_len // 4 + 1))[:ref_len]
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(path))
    return path


def build_test_index(tmp_path_factory, gtf_text, genome_size=2000, name="idx", refs=None):
    """Build a TranscriptIndex from a GTF string (session/module-scoped helper).

    Parameters
    ----------
    tmp_path_factory : pytest.TempPathFactory
        Pytest factory for session/module-scoped temp dirs.
    gtf_text : str
        GTF content string (GENCODE-style, 1-based inclusive).
    genome_size : int
        Length of chr1 in the synthetic FASTA (default: 2000).
    name : str
        Sub-directory name for this index (keeps separate indexes apart).
    refs : dict[str, int] | None
        Optional ``{reference_name: length}`` for a MULTI-reference index. ``None`` (the default) builds
        the single ``chr1`` of ``genome_size`` bp that every existing caller expects, so this parameter
        is purely additive. The splice-graph matrix needs it: region and boundary ids must stay
        contiguous per reference, and no boundary may cross references.

    Returns
    -------
    TranscriptIndex
        Loaded index with C++ FragmentResolver ready.
    """
    base = tmp_path_factory.mktemp(name)
    gtf_path = base / "test.gtf"
    gtf_path.write_text(gtf_text)

    fasta_path = write_fasta(base / "genome.fa", refs or {"chr1": genome_size})

    idx_dir = base / "index"
    TranscriptIndex.build(fasta_path, gtf_path, idx_dir, write_tsv=False)
    return TranscriptIndex.load(idx_dir, retain_test_structures=True)


class StandInStore:
    """The public surface of an ``alignable`` store the index reads (`rigel.observable`): ``status``,
    ``chromosomes``, ``read_length_bins``, ``unmappable_regions`` and the provenance attributes, so the suite needs
    neither the package nor a store.

    *zero* gives, per reference and read-length bin, the runs no covering read pair was placed back at — what
    ``AlignableStore.unmappable_regions`` returns at threshold 0 — and every other base is placed back. ``calls``
    records each read, ``(chrom, start, end, read_length)``.
    """

    def __init__(self, chromosomes, *, bins=(100,), zero=None, status="complete"):
        self.chromosomes = {str(k): int(v) for k, v in chromosomes.items()}
        self.read_length_bins = sorted(int(b) for b in bins)
        self.status = status
        self.aligner, self.aligner_version = "star", "2.7.11b"
        self.frag_len_mean, self.frag_len_sd, self.error_rate = 250, 50, 0.01
        self._zero = zero or {}
        self.calls = []

    def unmappable_regions(self, chrom, start, end, read_length, *, threshold=0.0):
        import numpy as np

        if threshold != 0.0:
            raise AssertionError("the index reads the zero set, threshold 0")
        self.calls.append((chrom, start, end, read_length))
        flags = np.zeros(end - start, dtype=bool)
        for s, e in self._zero.get(chrom, {}).get(read_length, ()):
            lo, hi = max(s, start), min(e, end)
            if lo < hi:
                flags[lo - start : hi - start] = True
        edges = np.diff(np.concatenate(([0], flags.astype(np.int8), [0])))
        return [
            (chrom, start + int(a), start + int(b))
            for a, b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1), strict=True)
        ]


def stand_in_store(monkeypatch, store: StandInStore) -> StandInStore:
    """Make *store* the one the index opens for its unobservable mask, whatever path it is given."""
    import rigel.observable

    monkeypatch.setattr(rigel.observable, "open_alignable_store", lambda path: store)
    return store


def rebuild_with_splice_blacklist(index_dir: str | Path, blacklist: pd.DataFrame) -> None:
    """Rebuild the index in *index_dir* from the sources and flags its manifest records, now with
    *blacklist* as its splice blacklist.

    The blacklist goes through `build`'s own alignable-store path, so the manifest records its source
    as a real build's does. The store loader is patched to return *blacklist* (the persisted columns:
    ``ref``, ``start``, ``end``, ``max_anchor_left``, ``max_anchor_right``), the store's mask is a
    :class:`StandInStore` that places every base back, and the store is an empty stand-in directory beside the
    index, so the suite needs neither the ``alignable`` package nor a store.
    """
    import pytest

    import rigel.splice_blacklist

    index_dir = Path(index_dir)
    manifest = load_manifest(index_dir)
    store = index_dir.parent / f"{index_dir.name}_alignable_stand_in"
    store.mkdir(exist_ok=True)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(
            rigel.splice_blacklist,
            "load_splice_blacklist_from_zarr",
            lambda store_path, *, min_count: blacklist,
        )
        stand_in_store(mp, StandInStore(load_ref_lengths(index_dir / REF_LENGTHS_FEATHER)))
        TranscriptIndex.build(
            manifest["sources"]["fasta"]["path"],
            manifest["sources"]["gtf"]["path"],
            index_dir,
            alignable_zarr_path=store,
            **manifest["build_flags"],
        )
