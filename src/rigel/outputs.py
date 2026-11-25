"""``rigel quant``'s output directory as one contract, read by the writer (``rigel quant``), the exporter (``rigel
export``) and the reader (``rigel report``): every product's file name, the derivatives a product may have, the summary's
schema version, and the rules that keep a directory one run's.

- ``summary.json`` is the completion marker. A run removes it before anything else, then every product and derivative
  an earlier run may have left, and writes it last; a directory without it holds no complete run, and the report and
  the exporter refuse it.
- Every file is written to a temporary sibling and moved into place, so a reader never sees a partial write.
- A product's derivatives (its ``--tsv`` mirror, what ``rigel export`` writes) are this module's to remove; a file it
  does not name is never touched.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path

#: summary.json's layout version, which a consumer reads before parsing; ``rigel report`` reads this one only.
SUMMARY_SCHEMA_VERSION = 11

#: The completion marker and the run's own configuration.
SUMMARY = "summary.json"
CONFIG = "config.yaml"

#: The quantification tables every run writes, in the order written.
QUANT_TABLES = ("quant", "gene_quant", "nascent_quant", "loci")
#: The raw fragment-length histograms, the genome-wide gDNA track and the per-locus EM profile (written when the run
#: collects it).
FRAGMENT_LENGTHS = "fragment_lengths"
CALIBRATION_TRACK = "calibration_track"
LOCUS_STATS = "locus_stats"
#: Every product written as Feather.
TABLES = (*QUANT_TABLES, FRAGMENT_LENGTHS, CALIBRATION_TRACK, LOCUS_STATS)
#: The track's genome-browser view.
BEDGRAPH = f"{CALIBRATION_TRACK}.bedgraph"
#: A Feather product's derivatives: its ``--tsv`` mirror, and what ``rigel export`` writes.
DERIVED_SUFFIXES = (".tsv", ".parquet")


def table_path(output_dir, name: str) -> Path:
    """A product's Feather file."""
    return Path(output_dir) / f"{name}.feather"


def managed_paths(output_dir) -> list[Path]:
    """Every file ``rigel quant`` writes or ``rigel export`` derives in *output_dir*, the completion marker first."""
    d = Path(output_dir)
    paths = [d / SUMMARY, d / CONFIG, d / BEDGRAPH]
    for name in TABLES:
        paths.append(table_path(d, name))
        paths += [d / f"{name}{suffix}" for suffix in DERIVED_SUFFIXES]
    return paths


def begin_run(output_dir) -> None:
    """Retire the previous run before anything of this one is written: its completion marker first, so a reader never
    takes a half-rewritten directory for a run, then every product and derivative it may have left, so nothing of
    another run survives beside this one's (a mirror, an export, a profile this run does not collect)."""
    for path in managed_paths(output_dir):
        path.unlink(missing_ok=True)


def publish(path, write: Callable[[Path], None]) -> None:
    """Write *path* by ``write(temporary sibling)``, then move the sibling into place: a reader sees no file or the
    whole one, never a partial write, and a failed write leaves no sibling behind."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        write(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def is_complete(output_dir) -> bool:
    """Whether *output_dir* holds a complete run: its completion marker is there."""
    return (Path(output_dir) / SUMMARY).is_file()
