"""
rigel.stats — Structured pipeline statistics tracking.

Replaces ad-hoc tally variables with a structured dataclass,
reducing boilerplate and enabling consistent serialization.
"""

from dataclasses import dataclass, asdict


@dataclass
class PipelineStats:
    """Structured counters for the rigel pipeline.

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
    n_chimeric_trans: int = 0
    n_chimeric_cis_strand_same: int = 0
    n_chimeric_cis_strand_diff: int = 0
    n_intergenic_unspliced: int = 0
    n_intergenic_spliced: int = 0
    n_with_exon: int = 0
    n_with_annotated_sj: int = 0
    n_with_unannotated_sj: int = 0
    n_same_strand: int = 0
    n_ambig_strand: int = 0

    # --- Strand model training ---
    n_strand_trained: int = 0
    n_strand_skipped_no_sj: int = 0
    n_strand_skipped_ambig_strand: int = 0
    n_strand_skipped_ambiguous: int = 0

    # --- Fragment length model training ---
    n_frag_length_unambiguous: int = 0
    n_frag_length_ambiguous: int = 0
    n_frag_length_intergenic: int = 0

    # --- Routing counters (authoritative) ---
    deterministic_unambig_units: int = 0
    em_routed_unambig_units: int = 0
    em_routed_ambig_same_strand_units: int = 0
    em_routed_ambig_opp_strand_units: int = 0
    em_routed_multimapper_units: int = 0

    # --- Gating ---
    n_gated_out: int = 0

    # --- Multimapper scan stats (NH > 1) ---
    n_multimapper_groups: int = 0
    n_multimapper_alignments: int = 0

    # --- gDNA contamination ---
    n_gdna_em: int = 0

    @property
    def n_intergenic(self) -> int:
        """Total intergenic fragments (unspliced + spliced)."""
        return self.n_intergenic_unspliced + self.n_intergenic_spliced

    @property
    def n_gdna_unambig(self) -> int:
        """Intergenic fragments assigned deterministically to gDNA."""
        return self.n_intergenic

    @property
    def n_gdna_total(self) -> int:
        """Total fragments assigned to gDNA (intergenic + EM)."""
        return self.n_gdna_unambig + self.n_gdna_em

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary.

        Includes computed properties (``n_intergenic``, ``n_gdna_*``)
        for convenience alongside the raw dataclass fields.
        """
        d = asdict(self)
        # Add computed properties
        d["n_intergenic"] = self.n_intergenic
        d["n_gdna_unambig"] = self.n_gdna_unambig
        d["n_gdna_total"] = self.n_gdna_total
        return d
