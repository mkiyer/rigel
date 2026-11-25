"""Load the report substrate from a ``rigel quant`` output directory written by this version of rigel.

The substrate is ``summary.json``, read and checked eagerly, and the feather tables the report draws
(:data:`REPORT_TABLES`), read lazily on first access; each is required. The report reads the summary only through
:class:`SummaryReader`, and only the keys :data:`SUMMARY_KEYS` lists, every one of which must be there: a summary of
another schema or one missing a key is refused, never read with a default.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import pandas as pd

from .. import outputs

#: Every ``summary.json`` key the report reads, by section: ``""`` is the top level and ``fragment_length.*`` each
#: fragment-length category. Each must be in the summary; a missing one refuses the report.
SUMMARY_KEYS: dict[str, tuple[str, ...]] = {
    "": ("schema_version", "rigel_version", "timestamp", "configuration"),
    "command.arguments": ("bam_file", "index_dir"),
    "alignment_stats": (
        "total_reads",
        "mapped_reads",
        "secondary_reads",
        "supplementary_reads",
        "proper_pairs",
        "duplicate_reads",
        "qc_fail_reads",
    ),
    "read_pairs": ("total", "excluded", "excluded_by_option"),
    "read_pairs.withheld": ("unpaired", "unannotated_splice", "off_index", "chimeric"),
    "read_pairs.placed": ("unique", "multimapping"),
    "fragment_stats": ("total", "genic", "intergenic"),
    "fragment_stats.splice": ("unspliced", "spliced_annotated", "spliced_implicit"),
    "strand_model": ("strand_specificity", "p_r1_sense", "read1_sense", "ci_95", "strand_live"),
    "strand_model.diagnostics": ("exonic_all_specificity", "contamination_gap"),
    "calibration": (
        "gdna_density_global",
        "gdna_bulk_density",
        "gdna_capture_saturation",
        "fraction_rna_sense",
        "n_regions",
    ),
    "fragment_length.*": (
        "n_observations",
        "mean",
        "std",
        "median",
        "mode",
        "overflow_fraction",
    ),
    "quantification": (
        "n_transcripts",
        "n_genes",
        "n_loci",
        "count_annotated",
        "count_nascent",
        "count_gdna_em",
        "fraction_annotated",
        "fraction_nascent",
        "fraction_gdna",
        "fraction_gdna_em",
        "count_intergenic",
        "fraction_intergenic",
    ),
}

#: The tables ``rigel quant`` writes that the report draws.
REPORT_TABLES = tuple(
    f"{name}.feather"
    for name in (outputs.FRAGMENT_LENGTHS, "gene_quant", outputs.CALIBRATION_TRACK)
)


class SummaryReader:
    """``summary.json`` as the report reads it: a section reads as exactly the keys :data:`SUMMARY_KEYS` lists, each
    the summary's value. A section or key the summary lacks raises :class:`SubstrateError`; one the list does not
    name raises ``KeyError``."""

    def __init__(self, summary: dict):
        self._summary = summary

    def __getitem__(self, name: str) -> dict:
        keys = SUMMARY_KEYS[name]
        if name.endswith(".*"):
            return {
                cat: _listed(held, keys, f"{name[:-2]}.{cat}")
                for cat, held in self._section(name[:-2]).items()
            }
        return _listed(self._section(name), keys, name)

    def _section(self, name: str) -> dict:
        node = self._summary
        for step in name.split(".") if name else ():
            if not (isinstance(node, dict) and isinstance(node.get(step), dict)):
                raise SubstrateError(f"summary.json has no `{name}` section")
            node = node[step]
        return node


def _listed(held: dict, keys: tuple[str, ...], where: str) -> dict:
    """The listed keys of one section, each the summary's value."""
    missing = [key for key in keys if key not in held]
    if missing:
        raise SubstrateError(f"summary.json's `{where or 'top level'}` lacks {missing}")
    return {key: held[key] for key in keys}


class SubstrateError(Exception):
    """Raised when the report substrate is missing or unreadable."""


@dataclass
class ReportSubstrate:
    """Everything the report builder needs from a quant output directory: ``summary`` loaded up front, and each table
    :data:`REPORT_TABLES` names read on first use."""

    output_dir: Path
    summary: dict

    def _feather(self, name: str) -> pd.DataFrame:
        path = self.output_dir / name
        try:
            return pd.read_feather(path)
        except Exception as exc:  # pragma: no cover - a corrupt file is unusual
            raise SubstrateError(f"Could not read {path.name}: {exc}") from exc

    # -- tables the report draws --
    @cached_property
    def fragment_lengths(self) -> pd.DataFrame:
        return self._feather("fragment_lengths.feather")

    @cached_property
    def gene_quant(self) -> pd.DataFrame:
        return self._feather("gene_quant.feather")

    @cached_property
    def calibration_track(self) -> pd.DataFrame:
        return self._feather("calibration_track.feather")

    @property
    def sample_name(self) -> str:
        """A human label for the sample: the input BAM's stem."""
        return Path(SummaryReader(self.summary)["command.arguments"]["bam_file"]).stem


def load_substrate(output_dir: str | Path) -> ReportSubstrate:
    """Load the report substrate from a ``rigel quant`` output directory written by this version of rigel: its
    ``summary.json``, at this version's schema with every key the report reads, and every table it draws. The tables
    load lazily.

    Raises
    ------
    SubstrateError
        If the directory, ``summary.json`` or a table is missing or unreadable, the summary is of another schema, or it
        lacks a key the report reads.
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise SubstrateError(f"Not a directory: {output_dir}")

    summary_path = output_dir / outputs.SUMMARY
    if not outputs.is_complete(output_dir):
        raise SubstrateError(
            f"No {outputs.SUMMARY} in {output_dir}: it holds no complete run (a `rigel quant` writes its summary "
            "last). Point `rigel report` at a completed `rigel quant` output directory."
        )
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SubstrateError(f"Could not parse {summary_path}: {exc}") from exc
    if not isinstance(summary, dict):
        raise SubstrateError(f"{summary_path} is not a summary: it holds no mapping")

    version = summary.get("schema_version")
    if version != outputs.SUMMARY_SCHEMA_VERSION:
        raise SubstrateError(
            f"{summary_path} is summary schema {version}, and this rigel reads schema "
            f"{outputs.SUMMARY_SCHEMA_VERSION}: run this version's `rigel quant` to report on it."
        )
    # every listed key read once here, so a missing one refuses the report before it is built
    reader = SummaryReader(summary)
    for name in SUMMARY_KEYS:
        reader[name]
    missing = [name for name in REPORT_TABLES if not (output_dir / name).exists()]
    if missing:
        raise SubstrateError(f"{output_dir} lacks {missing}, which every `rigel quant` run writes")
    return ReportSubstrate(output_dir=output_dir, summary=summary)
