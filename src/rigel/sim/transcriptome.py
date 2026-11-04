"""The transcriptome the simulator draws from: the rigel index's rows as transcripts, and nascent molecules pooled
onto the index's nascent entities.

Nascent abundance lives on the entities: every multi-exon transcript links through ``span_t_index`` to the entity
over its TSS/TES cluster and keeps ``nrna_abundance = 0`` itself. Sampling is one multinomial over all RNA rows,
transcript and entity alike, in proportion to abundance x effective length, so the nascent fragment share follows
from molecules and lengths and is never imposed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ..index import TRANSCRIPT_ALIASES_FEATHER, TranscriptIndex
from ..transcript import Transcript
from ..types import Interval, Strand

logger = logging.getLogger(__name__)

__all__ = ["assign_nrna_to_entities", "load_transcripts_from_index"]


def _aliases(index_dir: Path) -> dict[str, list[str]]:
    """``{kept t_id: [dropped t_ids]}``: the duplicates the index folded into each transcript it kept."""
    aliases: dict[str, list[str]] = {}
    alias_path = index_dir / TRANSCRIPT_ALIASES_FEATHER
    if alias_path.exists():
        columns = pd.read_feather(alias_path)[["dropped_t_id", "kept_t_id"]]
        for dropped, kept in columns.itertuples(index=False):
            aliases.setdefault(str(kept), []).append(str(dropped))
    return aliases


def load_transcripts_from_index(index_dir: str | Path) -> list[Transcript]:
    """One :class:`Transcript` per row of the rigel index, annotated transcripts and nascent entities, with the
    index's structure flags and ``t_index`` equal to the row, so the simulator and `rigel quant` see the same set.

    The index, not the GTF, is read because only it knows the duplicate collapse and the nascent entities. Each
    transcript carries the ids of the duplicates folded into it (``aliases``), so a probe panel designed on the raw
    annotation still binds. Raises if the rows are not contiguous ``t_index`` 0..n-1 or a row has no exons.
    """
    index = TranscriptIndex.load(index_dir, retain_lookup_maps=True)
    aliases = _aliases(Path(index_dir))
    transcripts: list[Transcript] = []
    for row in index.t_df.itertuples(index=False):
        exons = index.get_exon_intervals(int(row.t_index))
        if exons is None or len(exons) == 0:
            raise ValueError(f"index row {row.t_index} ({row.t_id}) has no exon intervals")
        t = Transcript(
            ref=str(row.ref),
            strand=Strand(int(row.strand)),
            exons=[Interval(int(s), int(e)) for s, e in np.asarray(exons).tolist()],
            t_id=str(row.t_id),
            g_id=str(row.g_id),
            g_name=str(row.g_name) if row.g_name is not None else None,
            g_type=str(row.g_type) if row.g_type is not None else None,
            t_index=int(row.t_index),
            g_index=int(row.g_index),
            is_basic=bool(row.is_basic),
            is_mane=bool(row.is_mane),
            is_ccds=bool(row.is_ccds),
            is_single_exon=bool(row.is_single_exon),
            is_synthetic=bool(row.is_synthetic),
            span_t_index=int(row.span_t_index),
            span_n_contributors=int(row.span_n_contributors),
            aliases=tuple(sorted(aliases.get(str(row.t_id), ()))),
        )
        t.length = t.compute_length()
        transcripts.append(t)
    if [t.t_index for t in transcripts] != list(range(len(transcripts))):
        raise ValueError("index rows are not contiguous t_index 0..n-1")
    n_syn = sum(1 for t in transcripts if t.is_synthetic)
    logger.info(
        "Loaded %d transcripts from index %s (%d synthetic nascent entities)",
        len(transcripts),
        index_dir,
        n_syn,
    )
    return transcripts


def assign_nrna_to_entities(transcripts: list[Transcript], per_contributor: np.ndarray) -> int:
    """Pool each contributor's nascent molecules onto its nascent entity; return how many entities received any.

    ``per_contributor[i]`` is transcript ``i``'s nascent molecular abundance. Each multi-exon transcript links
    through ``span_t_index`` to one single-exon entity (synthetic, or an annotated transcript covering the span),
    whose ``nrna_abundance`` becomes the sum over its contributors; every other row ends at 0. Raises when a
    contributor is single-exon or has no entity: the list did not come from a rigel index.
    """
    by_index = {t.t_index: t for t in transcripts}
    for t in transcripts:
        t.nrna_abundance = 0.0
    n_entities = 0
    for t, amount in zip(transcripts, per_contributor):
        if amount <= 0.0:
            continue
        if len(t.exons) <= 1:
            raise ValueError(f"{t.t_id} is single-exon and cannot contribute nascent RNA")
        if t.span_t_index < 0 or t.span_t_index not in by_index:
            raise ValueError(
                f"{t.t_id} is multi-exon but links to no nascent entity (span_t_index="
                f"{t.span_t_index}); the transcript list did not come from a rigel index"
            )
        entity = by_index[t.span_t_index]
        if entity.nrna_abundance == 0.0:
            n_entities += 1
        entity.nrna_abundance += float(amount)
    return n_entities
