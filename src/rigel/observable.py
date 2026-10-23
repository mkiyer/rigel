"""The observable genome — the bases a read pair can be observed at.

A fragment is observed only where a read pair drawn there is placed back there. A base the reference does not know
(``N``) is never observed, and neither is one an aligner cannot place; both are mappability. The index records the
complement of the observable genome as a MASK: half-open runs on each reference, merged, written beside the regions
(``unobservable.feather``). Its rule:

* **With an alignable store** — the companion tool that tiles a genome with read pairs, aligns them and counts, per
  base and read-length bin, the pairs covering the base and those placed back where they came from — the mask is the
  store's ZERO set at its longest read-length bin: every base no covering pair was placed back at
  (`AlignableStore.unmappable_regions` at threshold 0). N is in it by its own measure, since no read maps to N, while
  a short gap whose pairs span it on known sequence stays observable, as it is to an aligner. The longest bin is the
  universal set: a longer read maps more, so what it cannot place is unobservable to the shorter reads too, but for a
  sliver. A reference the store does not hold takes its N runs, named in a warning; one it holds at another length is
  another assembly, and refused.
* **Without a store** — synthetic genomes, and panels whose simulator places every fragment by truth — the mask is the
  reference's unknown runs, where a base is unknown iff it is ``N`` or ``n``: the one unobservability known with no
  aligner's word, and the rule the simulator draws by (`rigel.sim.capture.sampler`,
  `tests/test_sim_unknown_sequence.py`).

Gate: ``tests/test_index_observable.py``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from .config import CONSTANTS

__all__ = [
    "UNOBSERVABLE_COLUMNS",
    "check_mask",
    "mask_from_fasta",
    "mask_from_store",
    "open_alignable_store",
    "true_runs",
    "unknown_bases",
    "unknown_runs",
]

logger = logging.getLogger(__name__)

#: The mask's on-disk columns: one half-open ``[start, end)`` run per row.
UNOBSERVABLE_COLUMNS = ("ref", "start", "end")

#: ASCII separates a letter's two cases by one bit, so ORing it in lower-cases a letter, and only ``N`` and ``n`` become
#: ``n``.
_CASE_BIT = ord("a") - ord("A")
_N_LOWER = ord("n")


def unknown_bases(seq: str | bytes | np.ndarray) -> np.ndarray:
    """``bool[len(seq)]`` — True where the base is unknown (``N`` or ``n``)."""
    if isinstance(seq, str):
        seq = seq.encode()
    raw = (
        np.frombuffer(seq, dtype=np.uint8)
        if isinstance(seq, bytes)
        else np.asarray(seq, dtype=np.uint8)
    )
    return (raw | _CASE_BIT) == _N_LOWER


def true_runs(mask: np.ndarray) -> np.ndarray:
    """``int64[n, 2]`` — the ``[start, end)`` runs where a boolean array is true, ascending."""
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.empty((0, 2), dtype=np.int64)
    edges = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    return np.stack([np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)], axis=1).astype(
        np.int64
    )


def unknown_runs(seq: str | bytes | np.ndarray) -> np.ndarray:
    """``int64[n, 2]`` — the ``[start, end)`` runs of unknown bases in a sequence, ascending."""
    return true_runs(unknown_bases(seq))


def _frame(runs_by_ref: Mapping[str, np.ndarray]) -> pd.DataFrame:
    """The mask as a frame, in the references' order and ascending within each."""
    refs, starts, ends = [], [], []
    for ref, runs in runs_by_ref.items():
        runs = np.asarray(runs, dtype=np.int64).reshape(-1, 2)
        refs.extend([ref] * runs.shape[0])
        starts.append(runs[:, 0])
        ends.append(runs[:, 1])
    return pd.DataFrame(
        {
            "ref": pd.Series(refs, dtype=object),
            "start": np.concatenate(starts) if starts else np.zeros(0, np.int64),
            "end": np.concatenate(ends) if ends else np.zeros(0, np.int64),
        }
    )


def _fasta_runs(fasta_file: str | Path, refs) -> dict[str, np.ndarray]:
    """Each named reference's unknown runs, from the FASTA's own sequence."""
    import pysam

    with pysam.FastaFile(str(fasta_file)) as fa:
        return {ref: unknown_runs(fa.fetch(ref)) for ref in refs}


def mask_from_fasta(fasta_file: str | Path, ref_lengths: Mapping[str, int]) -> pd.DataFrame:
    """The mask with no store: every reference's unknown runs."""
    return _frame(_fasta_runs(fasta_file, list(ref_lengths)))


def open_alignable_store(path: str | Path):
    """The alignable store at *path*, through the companion package's public API (``alignable.open``)."""
    try:
        import alignable
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError(
            "Reading an alignable store requires the 'alignable' package. Install it into the rigel "
            "environment, or rebuild the index with --no-mappability."
        ) from exc
    # no chunk cache: the windows below read each chunk once, in order
    return alignable.open(str(path), cache_size=0)


def _store_runs(store, ref: str, length: int, read_length: int, window: int) -> np.ndarray:
    """The store's zero set on one reference, read a window at a time and stitched where a run reaches a window's
    edge, so the window sets memory and never a run."""
    out: list[list[int]] = []
    for lo in range(0, length, window):
        hi = min(lo + window, length)
        for _ref, start, end in store.unmappable_regions(ref, lo, hi, read_length, threshold=0.0):
            if out and out[-1][1] == start:
                out[-1][1] = end
            else:
                out.append([int(start), int(end)])
    return np.asarray(out, dtype=np.int64).reshape(-1, 2)


def mask_from_store(
    store_path: str | Path,
    ref_lengths: Mapping[str, int],
    fasta_file: str | Path,
    *,
    window: int = CONSTANTS.resources.alignable_window_bp,
) -> tuple[pd.DataFrame, dict]:
    """The mask from an alignable store: its zero set at its longest read-length bin on every reference it holds at
    the FASTA's length, and the FASTA's unknown runs on any it does not hold. Returns ``(mask, record)``, the record
    being what the manifest keeps of the store's measure.

    Refused, before anything is written: a store that is not complete (a chunk never computed reads as unmappable
    everywhere it lies) and one that holds a reference at another length (another assembly).
    """
    store = open_alignable_store(store_path)
    status = str(store.status)
    if status != "complete":
        raise ValueError(
            f"The alignable store {store_path} is {status!r}, not 'complete': a chunk it never computed reads as "
            "unmappable, so its mask would hide sequence an aligner can place. Finish the store "
            "(alignable finalize) or rebuild the index with --no-mappability."
        )
    chroms = {str(k): int(v) for k, v in store.chromosomes.items()}
    differ = [r for r, n in ref_lengths.items() if r in chroms and chroms[r] != int(n)]
    if differ:
        shown = ", ".join(f"{r} ({ref_lengths[r]:,} against {chroms[r]:,})" for r in differ[:5])
        raise ValueError(
            f"The alignable store {store_path} holds {len(differ):,} of the FASTA's references at another length "
            f"({shown}): it was built for another assembly."
        )
    absent = [r for r in ref_lengths if r not in chroms]
    if absent:
        logger.warning(
            "%s of the FASTA's %s references are not in the alignable store (%s): their unobservable runs are "
            "their N runs.",
            f"{len(absent):,}",
            f"{len(ref_lengths):,}",
            ", ".join(absent[:5]) + (f" and {len(absent) - 5:,} more" if len(absent) > 5 else ""),
        )
    read_length = int(max(store.read_length_bins))
    from_fasta = _fasta_runs(fasta_file, absent) if absent else {}
    runs = {
        ref: (
            _store_runs(store, ref, int(n), read_length, window)
            if ref in chroms
            else from_fasta[ref]
        )
        for ref, n in ref_lengths.items()
    }
    record = {
        "rule": "alignable store: no covering read pair placed back, at the store's longest read-length bin",
        "read_length": read_length,
        "read_length_bins": [int(b) for b in store.read_length_bins],
        "aligner": str(store.aligner),
        "aligner_version": None if store.aligner_version is None else str(store.aligner_version),
        "frag_len_mean": int(store.frag_len_mean),
        "frag_len_sd": int(store.frag_len_sd),
        "error_rate": float(store.error_rate),
        "references_from_fasta": absent,
    }
    return _frame(runs), record


def check_mask(mask: pd.DataFrame, ref_lengths: Mapping[str, int]) -> pd.DataFrame:
    """The mask as :meth:`TranscriptIndex.load` takes it, or ``ValueError``: exactly the three columns, every reference
    one of the index's, every run inside its reference, and the runs of each reference ascending and merged (no two
    touch). Returns the frame with ``start``/``end`` as int64."""
    if tuple(mask.columns) != UNOBSERVABLE_COLUMNS:
        raise ValueError(
            f"the unobservable mask's columns are {list(mask.columns)}, not {list(UNOBSERVABLE_COLUMNS)}"
        )
    ref = mask["ref"].astype(str).to_numpy()
    start = mask["start"].to_numpy(np.int64)
    end = mask["end"].to_numpy(np.int64)
    order = {name: i for i, name in enumerate(ref_lengths)}
    unknown = sorted(set(ref) - set(order))
    if unknown:
        raise ValueError(f"the unobservable mask names references the index lacks: {unknown[:5]}")
    length = np.array([ref_lengths[r] for r in ref], dtype=np.int64)
    if np.any(start < 0) or np.any(end <= start) or np.any(end > length):
        raise ValueError("the unobservable mask holds a run that is empty or outside its reference")
    rid = np.array([order[r] for r in ref], dtype=np.int64)
    same = rid[1:] == rid[:-1]
    if np.any(rid[1:] < rid[:-1]) or np.any(same & (start[1:] <= end[:-1])):
        raise ValueError(
            "the unobservable mask is not in the references' order with each reference's runs ascending and merged"
        )
    return pd.DataFrame({"ref": pd.Series(ref, dtype=object), "start": start, "end": end})
