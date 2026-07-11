"""Load the ``rigel quant`` report substrate into a structured bundle.

The substrate is everything ``rigel quant`` writes to its output directory that
the report needs:

* ``summary.json`` — the lean run manifest (schema v2).
* ``fragment_lengths.feather`` — tidy ``(category, length, count)`` histograms.
* ``gene_quant.feather`` / ``quant.feather`` / ``nrna_quant.feather`` /
  ``loci.feather`` — the expression + pool tables.

Only ``summary.json`` is required; every companion table is optional so the
report degrades gracefully on partial or older outputs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

#: summary.json schema this builder targets. Reports still build on other
#: versions, but a mismatch is surfaced as a warning banner.
EXPECTED_SCHEMA_VERSION = 2


class SubstrateError(Exception):
    """Raised when the report substrate is missing or unreadable."""


def _read_feather(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_feather(path)
    except Exception as exc:  # pragma: no cover - corrupt file is unusual
        raise SubstrateError(f"Could not read {path.name}: {exc}") from exc


@dataclass
class ReportSubstrate:
    """Everything the report builder needs, loaded from a quant output directory."""

    output_dir: Path
    summary: dict
    fragment_lengths: pd.DataFrame | None = None
    quant: pd.DataFrame | None = None
    gene_quant: pd.DataFrame | None = None
    nrna_quant: pd.DataFrame | None = None
    loci: pd.DataFrame | None = None
    #: Non-fatal issues to surface in the report (schema drift, missing tables).
    warnings: list[str] = field(default_factory=list)

    @property
    def schema_version(self) -> int | None:
        return self.summary.get("schema_version")

    @property
    def sample_name(self) -> str:
        """A human label for the sample, derived from the input BAM stem."""
        bam = self.summary.get("input", {}).get("bam_file")
        if bam:
            return Path(bam).stem
        return self.output_dir.name


def load_substrate(output_dir: str | Path) -> ReportSubstrate:
    """Load the report substrate from a ``rigel quant`` output directory.

    Parameters
    ----------
    output_dir
        Directory containing ``summary.json`` and its companion feather tables
        (i.e. the ``--output-dir`` of a ``rigel quant`` run).

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
    sv = summary.get("schema_version")
    if sv is None:
        warnings.append(
            "summary.json has no schema_version (pre-v2 run); fragment-length "
            "distributions may be unavailable."
        )
    elif sv != EXPECTED_SCHEMA_VERSION:
        warnings.append(
            f"summary.json schema_version is {sv}; this builder targets "
            f"v{EXPECTED_SCHEMA_VERSION}. Some panels may be incomplete."
        )

    fl = _read_feather(output_dir / "fragment_lengths.feather")
    if fl is None:
        warnings.append(
            "fragment_lengths.feather not found; the fragment-length "
            "distribution charts will be omitted."
        )

    return ReportSubstrate(
        output_dir=output_dir,
        summary=summary,
        fragment_lengths=fl,
        quant=_read_feather(output_dir / "quant.feather"),
        gene_quant=_read_feather(output_dir / "gene_quant.feather"),
        nrna_quant=_read_feather(output_dir / "nrna_quant.feather"),
        loci=_read_feather(output_dir / "loci.feather"),
        warnings=warnings,
    )
