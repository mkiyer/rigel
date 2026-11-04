"""The simulated molecules: each transcript's base abundance, and the nascent RNA each condition adds. A helper
module, not a script.

The base abundances are set once per run, by a log-uniform draw (:func:`assign_random_abundances`) or from a salmon,
kallisto or sim file (:func:`assign_file_abundances`). The nascent mode is applied per condition, onto the index's
nascent entities: at a molecular ratio to the mature level (:func:`apply_nrna_ratio`) or sparse and independent of it
(:func:`apply_sparse_nrna`).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.sim.transcriptome import assign_nrna_to_entities
from rigel.transcript import Transcript

from _sim_config import AbundanceConfig

logger = logging.getLogger(__name__)


# ── the base abundances ──────────────────────────────────────────────────────────────────────────


def assign_base_abundances(transcripts: list[Transcript], config: AbundanceConfig) -> bool:
    """Set every transcript's base abundance, once for every condition, as ``config.mode`` says. Returns whether
    the abundance file supplied the nascent levels, which replace the nascent sweep."""
    if config.mode == "random":
        assign_random_abundances(transcripts, config)
        return False
    if config.mode != "file":
        raise ValueError(f"Unknown abundance mode: {config.mode}")
    if not config.file:
        raise ValueError("abundance.file must be specified for mode='file'")
    has_file_nrna = assign_file_abundances(transcripts, config.file)
    if has_file_nrna:
        logger.info("Abundance file provided nascent abundances — skipping the nascent sweep")
    return has_file_nrna


def assign_random_abundances(transcripts: list[Transcript], config: AbundanceConfig) -> None:
    """Draw each annotated transcript's ``abundance``: expressed with probability ``frac_expressed``, then
    LogUniform(min, max). ``nrna_abundance`` is set to 0; the nascent mode is applied per condition.
    """
    if config.min <= 0:
        raise ValueError(f"abundance min must be > 0 for log-uniform, got {config.min}")

    rng = np.random.default_rng(config.seed)
    n = len(transcripts)
    # The draw runs over the annotated rows only, in index order: a nascent entity has no mature abundance, and
    # the synthetic rows must not shift the stream.
    annotated = np.array([not t.is_synthetic for t in transcripts], dtype=bool)

    expressed = np.zeros(n, dtype=bool)
    expressed[annotated] = rng.random(int(annotated.sum())) < config.frac_expressed

    total_rna = np.zeros(n)
    n_expr = int(expressed.sum())
    total_rna[expressed] = np.exp(rng.uniform(np.log(config.min), np.log(config.max), size=n_expr))

    for i, t in enumerate(transcripts):
        t.abundance = float(total_rna[i])
        t.nrna_abundance = 0.0

    total_mrna = sum(t.abundance for t in transcripts)
    logger.info("Log-uniform abundances: %d/%d expressed, total mRNA=%.1f", n_expr, n, total_mrna)


def _detect_abundance_format(path: str | Path) -> tuple[str, pd.DataFrame]:
    """Detect an abundance file's format, plain or gzipped, and return ``(format_name, dataframe)``:

    salmon   — ``quant.sf``: ``Name``, ``TPM``.
    kallisto — ``abundance.tsv``: ``target_id``, ``tpm``.
    sim      — the molar truth table (``_sim_truth.write_molar_truth``): ``transcript_id``, ``mrna_abundance`` and
               optionally ``nrna_abundance``.

    A condition's truth table is refused: it holds realised fragments, not molecules.
    """
    # gzip is detected by its magic bytes, whatever the extension
    compression: str = "infer"
    try:
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                compression = "gzip"
    except OSError:
        pass

    df = pd.read_csv(path, sep="\t", compression=compression)
    cols = set(df.columns)

    # A condition's table is known by its realised column, whatever else it carries.
    if "observed_mrna_fragments" in cols:
        raise ValueError(
            f"{path} is a condition's truth table, the fragments one condition realised, not an "
            "abundance input; the molar table is the run's truth_abundances_nrna_<label>.tsv"
        )
    if {"Name", "TPM"}.issubset(cols):
        return "salmon", df
    if {"target_id", "tpm"}.issubset(cols):
        return "kallisto", df
    if {"transcript_id", "mrna_abundance"}.issubset(cols):
        return "sim", df
    # kallisto, case-insensitively
    if {"target_id", "tpm"}.issubset({c.lower() for c in df.columns}):
        df.columns = [c.lower() for c in df.columns]
        return "kallisto", df

    raise ValueError(
        f"Unrecognized abundance file format.  Columns: {sorted(cols)}.\n"
        "Expected one of:\n"
        "  salmon quant.sf  — Name, TPM\n"
        "  kallisto abundance.tsv — target_id, tpm\n"
        "  sim TSV — transcript_id, mrna_abundance [, nrna_abundance]"
    )


def _load_abundance_map(path: str | Path) -> tuple[dict[str, tuple[float, float | None]], str]:
    """Load ``{transcript_id: (abundance, nrna or None)}`` and the detected format; the nRNA value is ``None`` when
    the file has none."""
    fmt, df = _detect_abundance_format(path)
    abund_map: dict[str, tuple[float, float | None]] = {}
    if fmt == "salmon":
        for _, row in df.iterrows():
            abund_map[str(row["Name"])] = (float(row["TPM"]), None)
    elif fmt == "kallisto":
        for _, row in df.iterrows():
            abund_map[str(row["target_id"])] = (float(row["tpm"]), None)
    elif fmt == "sim":
        has_nrna = "nrna_abundance" in df.columns
        for _, row in df.iterrows():
            nrna: float | None = float(row["nrna_abundance"]) if has_nrna else None
            abund_map[str(row["transcript_id"])] = (float(row["mrna_abundance"]), nrna)
    return abund_map, fmt


def assign_file_abundances(transcripts: list[Transcript], tsv_path: str | Path) -> bool:
    """Load abundances from a salmon, kallisto or sim file (:func:`_detect_abundance_format`); a transcript the file
    does not name gets 0. Returns whether the file supplied nRNA values; when it did not, the caller applies the
    configured nascent sweep.
    """
    abund_map, fmt = _load_abundance_map(tsv_path)
    logger.info("Detected abundance format: %s (%d entries)", fmt, len(abund_map))

    matched = 0
    has_nrna_data = False
    # A file's nRNA value for a multi-exon transcript is pooled onto its entity, as in the ratio mode; a value on an
    # entity or single-exon row is kept directly.
    per_contributor = np.zeros(len(transcripts))

    for i, t in enumerate(transcripts):
        if t.t_id in abund_map:
            total_or_mrna, nrna = abund_map[t.t_id]
            t.abundance = total_or_mrna
            if nrna is not None:
                has_nrna_data = True
                if len(t.exons) > 1:
                    per_contributor[i] = nrna
                else:
                    t.nrna_abundance = nrna
            matched += 1
        else:
            t.abundance = 0.0
    if per_contributor.any():
        direct = {
            i: transcripts[i].nrna_abundance
            for i in range(len(transcripts))
            if transcripts[i].nrna_abundance > 0
        }
        assign_nrna_to_entities(transcripts, per_contributor)
        for i, amount in direct.items():
            transcripts[i].nrna_abundance += amount

    logger.info(
        "File abundances (%s): matched=%d/%d, has_nrna=%s",
        fmt,
        matched,
        len(transcripts),
        has_nrna_data,
    )
    return has_nrna_data


# ── the nascent modes ────────────────────────────────────────────────────────────────────────────


def apply_nrna_ratio(transcripts: list[Transcript], ratio: float) -> None:
    """Nascent RNA at a molecular ratio: ``ratio`` nascent molecules per molecule of every expressed multi-exon
    transcript, pooled onto its entity (``assign_nrna_to_entities``)."""
    per = np.array(
        [
            (t.abundance or 0.0) * ratio if (t.abundance or 0.0) > 0 and len(t.exons) > 1 else 0.0
            for t in transcripts
        ]
    )
    n_entities = assign_nrna_to_entities(transcripts, per)
    logger.info(
        "Set nRNA molecular ratio %.3g: %d contributors onto %d entities (mRNA=%.1f, nRNA=%.1f)",
        ratio,
        int(np.sum(per > 0)),
        n_entities,
        sum(t.abundance or 0.0 for t in transcripts),
        sum(t.nrna_abundance for t in transcripts),
    )


def _eligible_entities(transcripts: list[Transcript]) -> list[Transcript]:
    """The nascent entities an expressed multi-exon transcript names, each once, in transcript order: a silent
    gene's span must stay exactly 0."""
    by_index = {t.t_index: t for t in transcripts}
    eligible: list[Transcript] = []
    seen: set[int] = set()
    for t in transcripts:
        if (t.abundance or 0.0) <= 0.0 or len(t.exons) <= 1:
            continue
        entity = by_index.get(t.span_t_index)
        if entity is None:
            raise ValueError(
                f"transcript {t.t_id} has no nascent entity (span_t_index={t.span_t_index}); a "
                f"transcript list without entities is not a rigel index's"
            )
        if entity.t_index not in seen:
            seen.add(entity.t_index)
            eligible.append(entity)
    return eligible


def apply_sparse_nrna(
    transcripts: list[Transcript],
    abundance_range: tuple[float, float],
    *,
    on_fraction: float,
    seed: int,
) -> float:
    """Sparse nascent RNA: each eligible entity is transcribed with probability ``on_fraction`` and, if so, draws an
    absolute molecular abundance log-uniformly over ``abundance_range``. Returns the realised nascent:mature
    molecular ratio. Requires ``0 < lo <= hi`` and ``0 <= on_fraction <= 1``.

    The level is independent of the mature level, because the nascent:mature ratio reflects stability, not
    expression: a stable abundant transcript shows little nascent signal and an unstable rare one can show more
    nascent than mature. The draw is log-uniform because present levels span decades. The unit is the entity, not
    the transcript, so ``on_fraction`` is the fraction of gene spans: drawn per isoform, a gene with ``n`` isoforms
    would carry nascent with probability ``1 - (1 - p)^n``.
    """
    lo, hi = float(abundance_range[0]), float(abundance_range[1])
    if lo <= 0.0 or hi < lo:
        raise ValueError(
            f"abundance_range must satisfy 0 < lo <= hi for a LOG-uniform draw; got ({lo}, {hi})"
        )
    if not 0.0 <= on_fraction <= 1.0:
        raise ValueError(f"on_fraction must be between 0 and 1; got {on_fraction}")

    rng = np.random.default_rng(seed)
    for t in transcripts:
        t.nrna_abundance = 0.0
    eligible = _eligible_entities(transcripts)

    log_lo, log_hi = np.log(lo), np.log(hi)
    n_on = 0
    for entity in eligible:
        if rng.random() >= on_fraction:
            continue
        entity.nrna_abundance = float(np.exp(rng.uniform(log_lo, log_hi)) if hi > lo else lo)
        n_on += 1

    total_mrna = sum(t.abundance or 0.0 for t in transcripts)
    total_nrna = sum(t.nrna_abundance for t in transcripts)
    realized_ratio = total_nrna / total_mrna if total_mrna > 0 else 0.0
    logger.info(
        "Sparse nRNA: %d/%d gene spans ON (on_fraction=%.3g), level ~ logU(%.3g, %.3g), "
        "mRNA=%.1f, nRNA=%.1f, realised molar ratio=%.4g",
        n_on,
        len(eligible),
        on_fraction,
        lo,
        hi,
        total_mrna,
        total_nrna,
        realized_ratio,
    )
    return realized_ratio
