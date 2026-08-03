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

    # --- Splice-artifact blacklist (per-read, across all alignments) ---
    n_sj_observed: int = 0
    n_sj_blacklisted: int = 0

    # --- The per-fragment splice census (scanner QC; ``rigel report``'s splice breakdown) ---
    #
    # ⭐ ONE observation per fragment the scanner offers the accumulator — unique mapper, resolved,
    # non-chimeric. That population is STATED, which the predecessor's was not: these counts used to
    # be read off the fragment-length category models, so they silently counted only the fragments
    # that also yielded a length observation.
    #
    # ⚠ The field names are derived, not chosen — ``rigel.splice.census_field`` builds each one from
    # its :class:`~rigel.splice.SpliceType` member name, and the C++ builds the same key from
    # ``splice_type_label``. Renaming one of these by hand breaks that correspondence silently,
    # because the copy below uses ``dict.get(key, 0)``.
    n_census_unspliced: int = 0
    n_census_spliced_unannot: int = 0
    n_census_spliced_annot: int = 0
    n_census_spliced_implicit: int = 0
    n_census_splice_artifact: int = 0

    #: Censused but never offered to the accumulator's ``deposit()``: the fragment is not one
    #: molecule on one cut axis (chiefly blocks on more than one reference). It exists so that
    #:
    #:   Σ census − n_census_splice_artifact
    #:       == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered
    #:
    #: closes exactly. Those returns were a silent loss before this counter.
    n_deposit_not_offered: int = 0

    @property
    def n_intergenic(self) -> int:
        """Total intergenic fragments (unspliced + spliced)."""
        return self.n_intergenic_unspliced + self.n_intergenic_spliced

    @property
    def n_gdna_unambig(self) -> int:
        """Intergenic fragments assigned deterministically to gDNA."""
        return self.n_intergenic

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary.

        Includes computed properties (``n_intergenic``, ``n_gdna_unambig``)
        for convenience alongside the raw dataclass fields.
        """
        d = asdict(self)
        # Add computed properties
        d["n_intergenic"] = self.n_intergenic
        d["n_gdna_unambig"] = self.n_gdna_unambig
        return d
