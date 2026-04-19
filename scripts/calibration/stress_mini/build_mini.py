"""Build a mini-genome with controlled repetitive-block structure.

Design
------
The genome is a stitching of 500-bp "blocks":

  - **Unique** blocks: filled with random DNA (independent per block).
    Mappability ≈ 1.  Two flavours:

        - ``U-tx``   — hosts a transcript (exon-bearing)
        - ``U-int``  — intergenic unique filler

  - **Repetitive** blocks: a single 500-bp sequence copied N times
    throughout the genome, with N ∈ {2, 5, 10, 50}.  Reads falling
    on any copy multimap with NH ≈ N.  Mappability per-copy ≈ 1/N.

Transcripts are placed inside ``U-tx`` blocks.  Each transcript has
1–10 exons (100–1000 bp), with introns 500–10000 bp.  Exons + introns
are all placed inside the same ``U-tx`` block, so every transcript
lives in a fully-unique region.  Inter-transcript spacing is filled
by a random stitching of the remaining block types.

Truth knobs we keep track of:

  * ``repeat_copy_id``     — identifies the rep-motif group (0..3)
  * ``repeat_copies``      — 2/5/10/50
  * ``block_type``         — ``U-tx`` / ``U-int`` / ``R2`` / ``R5`` / ``R10`` / ``R50``
  * per-block genomic interval

The output is:

  * ``genome.fa`` — mini genome
  * ``genes.gtf`` — transcripts
  * ``blocks.tsv`` — truth block table (for analysis)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.sim.genome import MutableGenome
from rigel.sim.annotation import GeneBuilder

logger = logging.getLogger(__name__)

BLOCK_LEN = 500
REPEAT_COUNTS = (2, 5, 10, 50)


@dataclass
class Block:
    """One 500-bp block in the genome."""
    start: int
    end: int                     # exclusive
    kind: str                    # 'U-tx' | 'U-int' | 'R2' | 'R5' | 'R10' | 'R50'
    rep_id: int | None = None    # group id for repetitive blocks
    rep_copy: int | None = None  # 0..N-1 within the group
    t_id: str | None = None      # transcript carried (only U-tx)


@dataclass
class TranscriptSpec:
    """A transcript placed in a U-tx block."""
    t_id: str
    g_id: str
    strand: str
    exons: list[tuple[int, int]]      # absolute genome coords
    abundance: float
    n_exons: int
    length: int


def _rand_dna(rng: np.random.Generator, n: int) -> str:
    return "".join("ACGT"[i] for i in rng.integers(0, 4, size=n))


def _place_transcript(
    rng: np.random.Generator,
    block_start: int,
    block_len: int,
    t_id: str,
    g_id: str,
    strand: str,
    abundance: float,
) -> TranscriptSpec | None:
    """Pack a random transcript inside a single block of size `block_len`.

    Returns None if we can't fit even a 1-exon tx in the block.
    """
    n_exons = int(rng.integers(1, 11))              # 1..10
    exons: list[tuple[int, int]] = []
    pos = int(rng.integers(0, 200))                 # small padding at start
    remaining = block_len - pos - 200               # leave end-padding
    for _ in range(n_exons):
        if remaining < 120:
            break
        ex_len = int(rng.integers(100, min(1001, remaining) + 1))
        exons.append((block_start + pos, block_start + pos + ex_len))
        pos += ex_len
        remaining -= ex_len
        # intron (skip if last exon)
        intr_len = int(rng.integers(500, 10_001))
        if remaining < intr_len + 120:
            break
        pos += intr_len
        remaining -= intr_len
    if not exons:
        return None
    length = sum(e - s for s, e in exons)
    return TranscriptSpec(
        t_id=t_id, g_id=g_id, strand=strand, exons=exons,
        abundance=abundance, n_exons=len(exons), length=length,
    )


def build_mini_genome(
    *,
    n_transcripts: int = 60,
    n_unique_intergenic: int = 20,
    tx_block_len: int = 60_000,   # big enough for 10 exons × 1000 bp + introns
    seed: int = 42,
) -> tuple[MutableGenome, GeneBuilder, list[Block], list[TranscriptSpec], pd.DataFrame]:
    """Stitch the mini genome.

    Layout (randomly shuffled):

        [U-tx]  [U-tx]  [R2]  [R2]                  <- each R-copy is 500 bp
        [U-tx]  [R5]×5  [U-int]  [R10]×10  [U-tx]   <- so R5×5 adds 2500 bp etc.
        …
    """
    rng = np.random.default_rng(seed)

    # --- decide the block playlist ------------------------------------------
    playlist: list[tuple[str, int | None]] = []  # (kind, rep_id)
    for k in range(n_transcripts):
        playlist.append(("U-tx", None))
    for _ in range(n_unique_intergenic):
        playlist.append(("U-int", None))
    # One repeat group per N; the group contributes N copies of the same seq
    rep_groups: list[tuple[int, int]] = []   # (rep_id, copies)
    rep_motifs: dict[int, str] = {}
    next_rep_id = 0
    # Create two groups for each repetition level so we have more samples
    for n_copies in REPEAT_COUNTS:
        for _ in range(2):
            motif = _rand_dna(rng, BLOCK_LEN)
            rep_motifs[next_rep_id] = motif
            for _copy in range(n_copies):
                playlist.append((f"R{n_copies}", next_rep_id))
            rep_groups.append((next_rep_id, n_copies))
            next_rep_id += 1

    rng.shuffle(playlist)

    # --- compute positions and the overall length ---------------------------
    blocks: list[Block] = []
    cur = 0
    rep_copy_counter: dict[int, int] = {}
    for kind, rep_id in playlist:
        if kind == "U-tx":
            blen = tx_block_len
        else:
            blen = BLOCK_LEN
        b = Block(start=cur, end=cur + blen, kind=kind)
        if rep_id is not None:
            b.rep_id = rep_id
            b.rep_copy = rep_copy_counter.get(rep_id, 0)
            rep_copy_counter[rep_id] = b.rep_copy + 1
        blocks.append(b)
        cur += blen

    genome_length = cur
    logger.info("Mini genome length = %d bp across %d blocks",
                genome_length, len(blocks))

    # --- build the MutableGenome (random DNA) then overwrite R-blocks -------
    genome = MutableGenome(genome_length, seed=seed, name="minichr")
    # Overwrite repetitive blocks with the motif so all copies match
    for b in blocks:
        if b.rep_id is not None:
            genome.edit(b.start, rep_motifs[b.rep_id])

    # --- place transcripts in U-tx blocks -----------------------------------
    builder = GeneBuilder(genome, ref_name="minichr")
    tx_specs: list[TranscriptSpec] = []
    tx_idx = 0
    for b in blocks:
        if b.kind != "U-tx":
            continue
        strand = "+" if rng.random() < 0.5 else "-"
        abund = float(np.exp(rng.uniform(np.log(20.0), np.log(1000.0))))
        t_id = f"t{tx_idx:03d}"
        g_id = f"g{tx_idx:03d}"
        spec = _place_transcript(
            rng, b.start, b.end - b.start, t_id, g_id, strand, abund,
        )
        if spec is None:
            continue
        b.t_id = t_id
        tx_specs.append(spec)
        builder.add_gene(g_id, strand, [{
            "t_id": t_id,
            "exons": spec.exons,
            "abundance": abund,
        }])
        tx_idx += 1

    logger.info("Placed %d transcripts on %d U-tx blocks",
                len(tx_specs), sum(1 for b in blocks if b.kind == "U-tx"))

    # --- truth table --------------------------------------------------------
    block_df = pd.DataFrame([
        {
            "start": b.start, "end": b.end, "length": b.end - b.start,
            "kind": b.kind, "rep_id": b.rep_id, "rep_copy": b.rep_copy,
            "t_id": b.t_id,
        }
        for b in blocks
    ])

    return genome, builder, blocks, tx_specs, block_df


def expected_mappability(kind: str) -> float:
    """Truth per-block mappability: 1/N for Rn, 1 otherwise."""
    if kind.startswith("R"):
        return 1.0 / int(kind[1:])
    return 1.0


__all__ = [
    "BLOCK_LEN", "REPEAT_COUNTS",
    "Block", "TranscriptSpec",
    "build_mini_genome", "expected_mappability",
]
