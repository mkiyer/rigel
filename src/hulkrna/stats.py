"""
hulkrna.stats — Structured pipeline statistics tracking.

Replaces ad-hoc counter variables with a structured dataclass,
reducing boilerplate and enabling consistent serialization.
"""

from dataclasses import dataclass, field, asdict


@dataclass
class PipelineStats:
    """Structured counters for the hulkrna pipeline.

    BAM-level stats are populated by ``parse_bam_file()``.
    Resolution and model-training stats are populated by the pipeline
    phases.

    Use ``to_dict()`` for JSON serialization.
    """

    # --- BAM-level (populated by parse_bam_file) ---
    total: int = 0
    qc_fail: int = 0
    unmapped: int = 0
    secondary: int = 0
    supplementary: int = 0
    duplicate: int = 0
    n_read_names: int = 0
    unique: int = 0
    multimapping: int = 0
    proper_pair: int = 0
    improper_pair: int = 0
    mate_unmapped: int = 0

    # --- Resolution-level ---
    n_fragments: int = 0
    n_chimeric: int = 0
    n_chimeric_interchrom: int = 0
    n_chimeric_strand_same: int = 0
    n_chimeric_strand_diff: int = 0
    n_intergenic_unspliced: int = 0
    n_intergenic_spliced: int = 0
    n_with_exon: int = 0
    n_with_intron_fallback: int = 0
    n_with_annotated_sj: int = 0
    n_with_unannotated_sj: int = 0
    n_unique_gene: int = 0
    n_multi_gene: int = 0

    # --- Strand model training ---
    n_strand_trained: int = 0
    n_strand_skipped_no_sj: int = 0
    n_strand_skipped_multi_gene: int = 0
    n_strand_skipped_ambiguous: int = 0

    # --- Insert size model training ---
    n_insert_unambiguous: int = 0
    n_insert_ambiguous: int = 0
    n_insert_intergenic: int = 0

    # --- Unique counting (1 gene, 1 transcript, NH=1) ---
    n_counted_truly_unique: int = 0

    # --- Isoform-ambiguous (1 gene, N transcripts, NH=1) ---
    n_counted_isoform_ambig: int = 0

    # --- Gene-ambiguous (N genes, NH=1) ---
    n_counted_gene_ambig: int = 0

    # --- Multimapper (NH > 1) ---
    n_multimapper_groups: int = 0
    n_multimapper_alignments: int = 0
    n_counted_multimapper: int = 0

    # --- gDNA contamination ---
    n_gdna_em: int = 0

    # --- Per-gene gDNA evidence (for hierarchical shadow model) ---
    n_intronic_antisense: int = 0
    n_exonic_antisense: int = 0

    @property
    def n_intergenic(self) -> int:
        """Total intergenic fragments (backward-compatible)."""
        return self.n_intergenic_unspliced + self.n_intergenic_spliced

    @property
    def n_gdna_unique(self) -> int:
        """Intergenic fragments assigned deterministically to gDNA."""
        return self.n_intergenic

    @property
    def n_gdna_total(self) -> int:
        """Total fragments assigned to gDNA (intergenic + EM)."""
        return self.n_gdna_unique + self.n_gdna_em

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary.

        Includes computed properties (``n_intergenic``, ``n_gdna_*``)
        for convenience alongside the raw dataclass fields.
        """
        d = asdict(self)
        # Add computed properties
        d["n_intergenic"] = self.n_intergenic
        d["n_gdna_unique"] = self.n_gdna_unique
        d["n_gdna_total"] = self.n_gdna_total
        return d

    def as_bam_stats_dict(self) -> dict:
        """Return a mutable dict compatible with parse_bam_file's stats arg.

        ``parse_bam_file()`` expects a plain dict with string keys.
        This method returns one backed by this dataclass's BAM-level fields.
        """
        return _BamStatsProxy(self)


class _BamStatsProxy(dict):
    """Dict-like proxy that reads/writes BAM-level fields on a PipelineStats.

    Used so that ``parse_bam_file()`` can populate PipelineStats directly
    without knowing about the dataclass.
    """

    _BAM_KEYS = {
        'total', 'qc_fail', 'unmapped', 'secondary',
        'supplementary', 'duplicate',
        'n_read_names', 'unique', 'multimapping',
        'proper_pair', 'improper_pair', 'mate_unmapped',
    }

    def __init__(self, stats: PipelineStats):
        super().__init__()
        self._stats = stats
        # Initialize from current dataclass values
        for key in self._BAM_KEYS:
            super().__setitem__(key, getattr(stats, key))

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if key in self._BAM_KEYS:
            setattr(self._stats, key, value)

    def setdefault(self, key, default=None):
        if key not in self:
            self[key] = default
        return self[key]

    def __getitem__(self, key):
        if key in self._BAM_KEYS:
            return getattr(self._stats, key)
        return super().__getitem__(key)
