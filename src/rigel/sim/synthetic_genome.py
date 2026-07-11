#!/usr/bin/env python3
"""Generate a synthetic mini-genome with realistic gene/isoform structure.

Creates:
  - Random FASTA genome (default 10 Mb, single chromosome)
  - GTF annotation with ~50 genes, ~250 transcripts
    - 1-10 isoforms per gene (avg ~5)
    - Sense/antisense overlapping gene pairs
    - Multi-exon and single-exon transcripts
    - GT-AG splice sites injected into FASTA

Output structure:
  <outdir>/
    genome.fa
    genome.fa.fai
    genes.gtf

Reference generation is normally driven by ``scripts/sim/simulate_suite.py``.
Use ``simulate_suite.py --reference-only`` for a reference-only debugging run.

This is the **function path** used by the whole-genome suite + the reference-gen scripts. It shares
its low-level mechanics — DNA generation (``genome.random_dna_array``), FASTA writing
(``genome.write_fasta_file``), and splice-motif placement (``splice_motif.place_intron_motif``) —
with the **class path** (``MutableGenome`` + ``GeneBuilder``, used by ``Scenario``). The two paths
keep DISTINCT output contracts on purpose: this path writes ``genome.fa`` + a transcript+exon
``genes.gtf`` (the interface ``suite.py`` / ``analysis.py`` / ``build_toy_2exon_reference.py`` and the
pre-generated benchmarks depend on), whereas ``GeneBuilder`` writes ``{name}.fa`` + an exon-only
``annotations.gtf``. Only the GTF writers below (:func:`write_gtf`) differ by design.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .genome import random_dna_array, write_fasta_file
from .splice_motif import place_intron_motif

logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────

GENOME_LENGTH = 10_000_000   # 10 Mb
REF_NAME = "chr_syn"
SEED = 42
N_GENES = 50
MIN_ISOFORMS = 1
MAX_ISOFORMS = 10
TARGET_TRANSCRIPTS = 250

# Gene spacing and size
MIN_GENE_SPAN = 5_000
MAX_GENE_SPAN = 80_000
MIN_INTERGENIC_GAP = 10_000

# Exon structure
MIN_EXONS = 1
MAX_EXONS = 10
MIN_EXON_LEN = 80
MAX_EXON_LEN = 2000
MIN_INTRON_LEN = 500
MAX_INTRON_LEN = 10_000

# Fraction of genes that have an antisense overlapping partner
ANTISENSE_OVERLAP_FRAC = 0.3


@dataclass
class ExonDef:
    start: int
    end: int


@dataclass
class TranscriptDef:
    t_id: str
    gene_id: str
    gene_name: str
    strand: str
    exons: list[ExonDef]

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end


@dataclass
class GeneDef:
    gene_id: str
    gene_name: str
    strand: str
    start: int
    end: int
    transcripts: list[TranscriptDef] = field(default_factory=list)


# ── Genome generation ─────────────────────────────────────────────────────

def generate_random_genome(length: int, seed: int) -> str:
    """Generate a random DNA sequence (shared generator — see ``genome.random_dna_array``)."""
    return "".join(random_dna_array(length, np.random.default_rng(seed)))


def inject_splice_sites(seq_array: np.ndarray, transcripts: list[TranscriptDef]):
    """Inject GT-AG splice sites at intron boundaries.

    Modifies ``seq_array`` (a uint8 byte array) in-place. For + strand: GT at intron start, AG at
    intron end. For - strand: CT / AC (reverse complement of GT-AG). Introns shorter than 4 bp are
    skipped (the function-path policy). Placement is the shared :func:`place_intron_motif`.
    """
    def _edit(pos: int, bases: str) -> None:
        for i, b in enumerate(bases):
            seq_array[pos + i] = ord(b)

    for tx in transcripts:
        exons = tx.exons
        for i in range(len(exons) - 1):
            place_intron_motif(
                _edit, exons[i].end, exons[i + 1].start, tx.strand, on_short="skip"
            )


# ── Gene/transcript structure generation ──────────────────────────────────

def _generate_exon_structure(
    rng: np.random.Generator,
    gene_start: int,
    gene_end: int,
    n_exons: int,
    strand: str,
) -> list[ExonDef]:
    """Generate a set of non-overlapping exons within [gene_start, gene_end]."""
    available_span = gene_end - gene_start

    # Generate exon lengths
    exon_lens = rng.integers(MIN_EXON_LEN, MAX_EXON_LEN + 1, size=n_exons)

    # Generate intron lengths (n_exons - 1 introns)
    if n_exons > 1:
        intron_lens = rng.integers(MIN_INTRON_LEN, MAX_INTRON_LEN + 1, size=n_exons - 1)
    else:
        intron_lens = np.array([], dtype=int)

    # Total required length
    total_needed = int(exon_lens.sum() + intron_lens.sum())

    # Scale down if needed
    if total_needed > available_span:
        scale = available_span / total_needed * 0.9
        exon_lens = np.maximum(MIN_EXON_LEN, (exon_lens * scale).astype(int))
        if n_exons > 1:
            intron_lens = np.maximum(MIN_INTRON_LEN, (intron_lens * scale).astype(int))
        total_needed = int(exon_lens.sum() + intron_lens.sum())

    # If still too large, reduce exon count
    if total_needed > available_span:
        n_exons = max(1, n_exons // 2)
        exon_lens = exon_lens[:n_exons]
        intron_lens = intron_lens[:max(0, n_exons - 1)]
        total_needed = int(exon_lens.sum() + intron_lens.sum())

    # Place within gene region
    slack = available_span - total_needed
    offset = gene_start + (rng.integers(0, max(1, slack)) if slack > 0 else 0)

    exons = []
    pos = offset
    for i in range(n_exons):
        exon_start = pos
        exon_end = pos + int(exon_lens[i])
        exons.append(ExonDef(exon_start, exon_end))
        pos = exon_end
        if i < n_exons - 1:
            pos += int(intron_lens[i])

    return exons


def _generate_isoform(
    rng: np.random.Generator,
    base_exons: list[ExonDef],
    isoform_idx: int,
) -> list[ExonDef]:
    """Generate an isoform by modifying a base exon structure.

    Modifications:
    - Skip random internal exons (exon skipping)
    - Extend/shorten first/last exon (alternative TSS/TES)
    - Add alternative 5'/3' splice sites (shift boundaries)
    """
    n_base = len(base_exons)
    if n_base <= 1 or isoform_idx == 0:
        return list(base_exons)

    exons = list(base_exons)
    modification = rng.integers(0, 4)

    if modification == 0 and n_base > 2:
        # Exon skipping: remove 1-2 internal exons
        n_skip = min(rng.integers(1, 3), n_base - 2)
        skip_indices = sorted(rng.choice(range(1, n_base - 1), size=n_skip, replace=False))
        exons = [e for i, e in enumerate(exons) if i not in skip_indices]

    elif modification == 1:
        # Alternative first exon: shift start of first exon
        shift = rng.integers(-200, 200)
        new_start = max(0, exons[0].start + shift)
        if new_start < exons[0].end - MIN_EXON_LEN:
            exons[0] = ExonDef(new_start, exons[0].end)

    elif modification == 2:
        # Alternative last exon: shift end of last exon
        shift = rng.integers(-200, 500)
        new_end = exons[-1].end + shift
        if new_end > exons[-1].start + MIN_EXON_LEN:
            exons[-1] = ExonDef(exons[-1].start, new_end)

    elif modification == 3 and n_base > 1:
        # Alternative splice site: shift an internal boundary
        idx = rng.integers(0, n_base - 1)
        shift = rng.integers(-50, 50)
        new_end = exons[idx].end + shift
        if new_end > exons[idx].start + MIN_EXON_LEN:
            exons[idx] = ExonDef(exons[idx].start, new_end)

    return exons


def _exon_signature(exons: list[ExonDef]) -> tuple[tuple[int, int], ...]:
    """Canonical exon-coordinate tuple used to detect structural duplicates."""
    return tuple((e.start, e.end) for e in sorted(exons, key=lambda e: (e.start, e.end)))


def _dedup_isoforms(
    rng: np.random.Generator,
    isoforms: list[list[ExonDef]],
    genome_length: int,
    max_attempts: int = 25,
) -> list[list[ExonDef]]:
    """Ensure no two isoforms share an identical exon coordinate tuple.

    Isoforms with the same intron chain but differing UTR boundaries (i.e.
    differing only in the start of exon 0 and/or the end of the last exon)
    are *allowed* — that is a realistic, biologically-relevant ambiguity
    case that downstream logic must cope with.

    On collision the offending isoform's first/last exon boundary is
    randomly shifted (by at least ±1 bp, up to ±300 bp) until either the
    exon tuple becomes unique or ``max_attempts`` is exhausted, in which
    case a deterministic ``isoform_idx``-based shift is applied as a
    fallback so the function always returns a duplicate-free list.
    """
    seen: dict[tuple, int] = {}
    out: list[list[ExonDef]] = []
    for idx, exons in enumerate(isoforms):
        sig = _exon_signature(exons)
        if sig not in seen:
            seen[sig] = idx
            out.append(exons)
            continue

        # Collision — perturb UTR boundaries until unique
        cur = list(exons)
        for attempt in range(max_attempts):
            cur = _perturb_utr(rng, cur, genome_length)
            new_sig = _exon_signature(cur)
            if new_sig not in seen:
                break
        else:
            # Deterministic last-resort: shift exon-0 start by idx+1 bp
            shift = idx + 1
            first = cur[0]
            new_start = max(0, first.start + shift if first.start + shift < first.end - MIN_EXON_LEN
                            else first.start - shift)
            cur[0] = ExonDef(new_start, first.end)
            new_sig = _exon_signature(cur)
            # If still colliding (extremely unlikely), nudge the last exon too.
            if new_sig in seen:
                last = cur[-1]
                cur[-1] = ExonDef(last.start, min(genome_length, last.end + shift))
                new_sig = _exon_signature(cur)

        seen[new_sig] = idx
        out.append(cur)
    return out


def _perturb_utr(
    rng: np.random.Generator,
    exons: list[ExonDef],
    genome_length: int,
) -> list[ExonDef]:
    """Shift the 5' or 3' UTR boundary by a non-zero offset.

    Preserves the intron chain and total exon count; only the start of the
    first exon and/or end of the last exon change.
    """
    new_exons = list(exons)
    # Choose which end(s) to perturb
    end = rng.integers(0, 3)  # 0=first, 1=last, 2=both
    if end in (0, 2):
        first = new_exons[0]
        delta = int(rng.integers(-300, 301))
        if delta == 0:
            delta = 1
        new_start = max(0, first.start + delta)
        if new_start >= first.end - MIN_EXON_LEN:
            new_start = max(0, first.end - MIN_EXON_LEN - 1)
        new_exons[0] = ExonDef(new_start, first.end)
    if end in (1, 2):
        last = new_exons[-1]
        delta = int(rng.integers(-300, 301))
        if delta == 0:
            delta = 1
        new_end = min(genome_length, last.end + delta)
        if new_end <= last.start + MIN_EXON_LEN:
            new_end = last.start + MIN_EXON_LEN + 1
        new_exons[-1] = ExonDef(last.start, new_end)
    return new_exons


def generate_genes(
    genome_length: int,
    n_genes: int,
    seed: int,
    *,
    min_isoforms: int = MIN_ISOFORMS,
    max_isoforms: int = MAX_ISOFORMS,
    target_transcripts: int | None = TARGET_TRANSCRIPTS,
    antisense_overlap_frac: float = ANTISENSE_OVERLAP_FRAC,
) -> list[GeneDef]:
    """Generate gene definitions with isoforms.

    Produces approximately ``n_genes`` genes. When ``target_transcripts`` is
    provided, isoform counts use the historical exponential-like distribution
    scaled to that target. When it is ``None``, each gene gets a uniform random
    isoform count in ``[min_isoforms, max_isoforms]``.
    Some genes get antisense overlapping partners.
    """
    if n_genes <= 0:
        raise ValueError("n_genes must be > 0")
    if min_isoforms <= 0:
        raise ValueError("min_isoforms must be > 0")
    if max_isoforms < min_isoforms:
        raise ValueError("max_isoforms must be >= min_isoforms")
    if not 0.0 <= antisense_overlap_frac <= 1.0:
        raise ValueError("antisense_overlap_frac must be between 0 and 1")
    rng = np.random.default_rng(seed)
    genes: list[GeneDef] = []

    # Decide isoform counts per gene.
    if target_transcripts is None:
        isoform_counts = rng.integers(
            min_isoforms, max_isoforms + 1, size=n_genes,
        ).astype(int).tolist()
    else:
        # Draw from the original shifted exponential-like distribution, but
        # scale it to the requested target and always emit n_genes primary genes.
        scale = max(1.0, float(target_transcripts) / float(n_genes))
        isoform_counts = [
            min(max_isoforms, max(min_isoforms, int(rng.exponential(scale)) + 1))
            for _ in range(n_genes)
        ]

    actual_n_genes = len(isoform_counts)

    # Decide which genes get antisense partners
    n_antisense = int(actual_n_genes * antisense_overlap_frac)
    antisense_indices = set(rng.choice(actual_n_genes, size=n_antisense, replace=False))

    # Place genes along the genome
    # Reserve space for genes + intergenic gaps
    total_gene_slots = actual_n_genes + n_antisense  # antisense genes add slots
    avg_slot_size = genome_length // (total_gene_slots + 1)

    pos = MIN_INTERGENIC_GAP
    gene_idx = 0

    for i in range(actual_n_genes):
        n_isoforms = isoform_counts[i]
        strand = "+" if rng.random() < 0.5 else "-"
        gene_span = rng.integers(MIN_GENE_SPAN, min(MAX_GENE_SPAN, avg_slot_size))

        if pos + gene_span >= genome_length - MIN_INTERGENIC_GAP:
            break

        gene_id = f"GENE{gene_idx + 1:04d}"
        gene_name = f"Gene{gene_idx + 1}"

        # Generate base exon structure for the gene
        n_base_exons = rng.integers(MIN_EXONS, MAX_EXONS + 1)
        base_exons = _generate_exon_structure(rng, pos, pos + gene_span, n_base_exons, strand)

        gene_def = GeneDef(
            gene_id=gene_id,
            gene_name=gene_name,
            strand=strand,
            start=pos,
            end=pos + gene_span,
        )

        # Generate isoforms (collect exon lists first so we can dedup before
        # constructing TranscriptDefs).
        iso_exon_lists: list[list[ExonDef]] = []
        for iso_idx in range(n_isoforms):
            if iso_idx == 0:
                exons = list(base_exons)
            else:
                exons = _generate_isoform(rng, base_exons, iso_idx)

            # Validate exon bounds
            exons = [ExonDef(max(0, e.start), min(genome_length, e.end)) for e in exons]
            exons = [e for e in exons if e.end > e.start + 10]
            if not exons:
                exons = list(base_exons)
            iso_exon_lists.append(exons)

        iso_exon_lists = _dedup_isoforms(rng, iso_exon_lists, genome_length)
        for iso_idx, exons in enumerate(iso_exon_lists):
            t_id = f"{gene_id}.{iso_idx + 1}"
            tx = TranscriptDef(
                t_id=t_id,
                gene_id=gene_id,
                gene_name=gene_name,
                strand=strand,
                exons=exons,
            )
            gene_def.transcripts.append(tx)

        genes.append(gene_def)
        gene_idx += 1

        # Generate antisense partner if selected
        if i in antisense_indices:
            # Antisense gene overlaps part of the sense gene
            as_strand = "-" if strand == "+" else "+"
            as_start = pos + rng.integers(0, gene_span // 3)
            as_end = min(pos + gene_span, as_start + rng.integers(MIN_GENE_SPAN, gene_span))

            if as_end > as_start + MIN_GENE_SPAN:
                as_gene_id = f"GENE{gene_idx + 1:04d}"
                as_gene_name = f"Gene{gene_idx + 1}as"
                n_as_exons = rng.integers(2, min(8, MAX_EXONS + 1))
                n_as_isoforms = rng.integers(1, min(4, MAX_ISOFORMS + 1))

                as_base_exons = _generate_exon_structure(
                    rng, as_start, as_end, n_as_exons, as_strand
                )

                as_gene_def = GeneDef(
                    gene_id=as_gene_id,
                    gene_name=as_gene_name,
                    strand=as_strand,
                    start=as_start,
                    end=as_end,
                )

                as_iso_lists: list[list[ExonDef]] = []
                for iso_idx in range(n_as_isoforms):
                    if iso_idx == 0:
                        exons = list(as_base_exons)
                    else:
                        exons = _generate_isoform(rng, as_base_exons, iso_idx)
                    exons = [ExonDef(max(0, e.start), min(genome_length, e.end)) for e in exons]
                    exons = [e for e in exons if e.end > e.start + 10]
                    if not exons:
                        exons = list(as_base_exons)
                    as_iso_lists.append(exons)

                as_iso_lists = _dedup_isoforms(rng, as_iso_lists, genome_length)
                for iso_idx, exons in enumerate(as_iso_lists):
                    t_id = f"{as_gene_id}.{iso_idx + 1}"
                    tx = TranscriptDef(
                        t_id=t_id,
                        gene_id=as_gene_id,
                        gene_name=as_gene_name,
                        strand=as_strand,
                        exons=exons,
                    )
                    as_gene_def.transcripts.append(tx)

                genes.append(as_gene_def)
                gene_idx += 1

        # Advance position
        pos += gene_span + rng.integers(MIN_INTERGENIC_GAP, MIN_INTERGENIC_GAP * 3)

    return genes


# ── File writers ──────────────────────────────────────────────────────────

def write_fasta(genome_seq: str, ref_name: str, outdir: Path) -> Path:
    """Write the genome FASTA (as ``genome.fa``) with a samtools index.

    Uses the shared :func:`genome.write_fasta_file`; the only difference from
    ``MutableGenome.write_fasta`` is the chosen filename (``genome.fa`` — the suite interface).
    """
    return write_fasta_file(genome_seq, ref_name, Path(outdir) / "genome.fa")


def write_gtf(genes: list[GeneDef], ref_name: str, outdir: Path) -> Path:
    """Write a transcript+exon GTF annotation as ``genes.gtf``.

    Intentionally distinct from ``GeneBuilder.write_gtf`` (which writes an exon-only
    ``annotations.gtf`` via ``GTFRecord``): this format + filename are the interface the suite,
    ``analysis.py``'s ``rigel index`` call, ``build_toy_2exon_reference.py``, and the pre-generated
    benchmarks depend on, so the two GTF writers are kept separate by design (see the module
    docstring). The dinucleotide convention + splice-motif placement ARE shared (``splice_motif``).
    """
    gtf_path = outdir / "genes.gtf"
    with open(gtf_path, "w") as f:
        for gene in genes:
            for tx in gene.transcripts:
                # Transcript line
                tx_start = tx.start + 1  # GTF is 1-based
                tx_end = tx.end
                f.write(
                    f"{ref_name}\trigel_sim\ttranscript\t{tx_start}\t{tx_end}\t.\t"
                    f"{tx.strand}\t.\t"
                    f'gene_id "{tx.gene_id}"; transcript_id "{tx.t_id}"; '
                    f'gene_name "{tx.gene_name}";\n'
                )
                # Exon lines
                for i, exon in enumerate(tx.exons, 1):
                    e_start = exon.start + 1  # GTF is 1-based
                    e_end = exon.end
                    f.write(
                        f"{ref_name}\trigel_sim\texon\t{e_start}\t{e_end}\t.\t"
                        f"{tx.strand}\t.\t"
                        f'gene_id "{tx.gene_id}"; transcript_id "{tx.t_id}"; '
                        f'gene_name "{tx.gene_name}"; exon_number "{i}";\n'
                    )
    return gtf_path
