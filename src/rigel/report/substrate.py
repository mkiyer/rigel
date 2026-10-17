"""Load the ``rigel quant`` report substrate into a structured bundle.

The substrate is what the report reads from a ``rigel quant`` output directory:

* ``summary.json`` — the lean run manifest, read eagerly.
* the companion feather tables the report draws (``fragment_lengths``,
  ``gene_quant``, ``calibration_track``) — read lazily on first access, so a
  report only pays I/O and memory for the tables it actually uses.

Only ``summary.json`` is required; a missing companion table empties the
charts and tables drawn from it. The report reads the summary only through
:class:`SummaryReader`, and only the keys :data:`SUMMARY_KEYS` lists.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import pandas as pd

# The summary.json schema this builder targets is the one its writer writes,
# ``cli.SUMMARY_SCHEMA_VERSION``, read at call time and never restated here. Reports still build on
# other versions, but a mismatch is surfaced as a warning banner.
from .. import cli

#: Every ``summary.json`` key the report reads, by section, with the value a summary without it reads: the one
#: list. The reads are tolerant, so a summary of another schema version still builds a report behind the warning
#: banner; that is also why a key renamed or misspelled here would read its default without a sound, and
#: ``tests/test_report.py`` holds every one to what the schema-4 writer writes. ``""`` is the top level and
#: ``fragment_length.*`` each fragment-length category.
SUMMARY_KEYS: dict[str, dict] = {
    "": {
        "schema_version": None,
        "rigel_version": None,
        "timestamp": None,
        "configuration": None,
        "calibration": None,
    },
    "input": {"bam_file": None, "index_dir": None},
    "alignment_stats": {
        "total_reads": 0,
        "mapped_reads": 0,
        "secondary_reads": 0,
        "supplementary_reads": 0,
        "proper_pairs": 0,
        "duplicate_reads": 0,
        "qc_fail_reads": 0,
        "read_groups": 0,
        "unique_reads": 0,
        "multimapping_reads": 0,
    },
    "fragment_stats": {
        "total": 0,
        "genic": 0,
        "intergenic": 0,
        "chimeric": 0,
        "chimeric_trans": 0,
        "chimeric_cis_same": 0,
        "chimeric_cis_diff": 0,
    },
    "fragment_stats.splice": {
        "unspliced": 0,
        "spliced_annotated": 0,
        "spliced_unannotated": 0,
        "spliced_implicit": 0,
        "splice_artifact": 0,
        "sj_blacklisted": 0,
        "sj_blacklist_size": None,
        "sj_blacklist_loaded": None,
    },
    "strand_model": {
        "protocol": "?",
        "strand_specificity": 0.5,
        "p_r1_sense": 0.5,
        "read1_sense": True,
        "ci_95": (0.0, 1.0),
        "strand_live": False,
    },
    "strand_model.diagnostics": {"exonic_all_specificity": None, "contamination_gap": None},
    "calibration": {
        "gdna_density_global": 0,
        "gdna_bulk_density": None,
        "gdna_capture_saturation": None,
        "fraction_rna_sense": 0,
        "n_regions": 0,
    },
    "fragment_length.*": {
        "n_observations": 0,
        "mean": 0,
        "std": 0,
        "median": 0,
        "mode": 0,
        "overflow_fraction": 0.0,
    },
    "quantification": {
        "n_transcripts": 0,
        "n_genes": 0,
        "n_loci": 0,
        "count_annotated": 0.0,
        "count_nascent": 0.0,
        "count_gdna_em": 0.0,
        "count_intergenic_rna": 0,
        "fraction_annotated": 0.0,
        "fraction_nascent": 0.0,
        "fraction_gdna": 0.0,
        "fraction_gdna_em": 0.0,
        "fraction_intergenic": 0.0,
        "fraction_intergenic_rna": 0.0,
    },
}


class SummaryReader:
    """``summary.json`` as the report reads it. A section reads as exactly the keys :data:`SUMMARY_KEYS` lists for
    it, each the summary's own value or, where the summary holds none, the listed default; a missing or ``null``
    section reads as all defaults. A section or a key the list does not name raises ``KeyError``, so the report
    cannot read one the list leaves out."""

    def __init__(self, summary: dict):
        self._summary = summary

    def __getitem__(self, name: str) -> dict:
        keys = SUMMARY_KEYS[name]
        if name.endswith(".*"):
            return {cat: _listed(held, keys) for cat, held in self._section(name[:-2]).items()}
        return _listed(self._section(name), keys)

    def holds(self, path: str) -> bool:
        """Whether the summary holds the listed key ``path`` names (``section.key``), whatever its value."""
        name, _, key = path.rpartition(".")
        if key not in SUMMARY_KEYS[name]:
            raise KeyError(path)
        return key in self._section(name)

    def _section(self, name: str) -> dict:
        node = self._summary
        for step in name.split(".") if name else ():
            node = node.get(step) if isinstance(node, dict) else None
        return node if isinstance(node, dict) else {}


def _listed(held, keys: dict) -> dict:
    """The listed keys of one section, each the summary's value or its default."""
    held = held if isinstance(held, dict) else {}
    return {key: held.get(key, default) for key, default in keys.items()}


class SubstrateError(Exception):
    """Raised when the report substrate is missing or unreadable."""


@dataclass
class ReportSubstrate:
    """Everything the report builder needs from a quant output directory.

    ``summary`` is loaded up front; the feather tables are lazy ``cached_property``
    accessors that read from ``output_dir`` on first use and return ``None`` when
    the file is absent.
    """

    output_dir: Path
    summary: dict
    #: Non-fatal issues to surface in the report (schema drift, missing tables).
    warnings: list[str] = field(default_factory=list)

    def _feather(self, name: str) -> pd.DataFrame | None:
        path = self.output_dir / name
        if not path.exists():
            return None
        try:
            return pd.read_feather(path)
        except Exception as exc:  # pragma: no cover - corrupt file is unusual
            raise SubstrateError(f"Could not read {path.name}: {exc}") from exc

    # -- tables the report currently consumes --
    @cached_property
    def fragment_lengths(self) -> pd.DataFrame | None:
        return self._feather("fragment_lengths.feather")

    @cached_property
    def gene_quant(self) -> pd.DataFrame | None:
        return self._feather("gene_quant.feather")

    @cached_property
    def calibration_track(self) -> pd.DataFrame | None:
        return self._feather("calibration_track.feather")

    @property
    def schema_version(self) -> int | None:
        return SummaryReader(self.summary)[""]["schema_version"]

    @property
    def sample_name(self) -> str:
        """A human label for the sample, derived from the input BAM stem."""
        bam = SummaryReader(self.summary)["input"]["bam_file"]
        return Path(bam).stem if bam else self.output_dir.name


def load_substrate(output_dir: str | Path) -> ReportSubstrate:
    """Load the report substrate from a ``rigel quant`` output directory.

    Reads ``summary.json`` and checks which companion tables are present (for
    warnings); the tables themselves load lazily on first access.

    Raises
    ------
    SubstrateError
        If the directory or ``summary.json`` is missing / unreadable.
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise SubstrateError(f"Not a directory: {output_dir}")

    summary_path = output_dir / "summary.json"
    if not summary_path.exists():
        raise SubstrateError(
            f"No summary.json in {output_dir}. Point `rigel report` at a "
            f"`rigel quant` output directory."
        )
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SubstrateError(f"Could not parse {summary_path}: {exc}") from exc

    warnings: list[str] = []
    sv = SummaryReader(summary)[""]["schema_version"]
    if sv != cli.SUMMARY_SCHEMA_VERSION:
        warnings.append(
            f"summary.json schema_version is {sv}; this builder targets "
            f"v{cli.SUMMARY_SCHEMA_VERSION}. Some panels may be incomplete."
        )
    if not (output_dir / "fragment_lengths.feather").exists():
        warnings.append(
            "fragment_lengths.feather not found; the fragment-length "
            "distribution charts will be omitted."
        )

    return ReportSubstrate(output_dir=output_dir, summary=summary, warnings=warnings)
