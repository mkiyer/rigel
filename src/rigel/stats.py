"""
rigel.stats — Structured pipeline statistics tracking.

Replaces ad-hoc tally variables with a structured dataclass,
reducing boilerplate and enabling consistent serialization.
"""

from dataclasses import dataclass, asdict


@dataclass
class PipelineStats:
    """Structured counters for the rigel pipeline.

    The scan counters are copied from the native scanner's stats by
    ``pipeline._apply_scan_stats``; the routing and gating counters are
    added by ``FragmentRouter.route``.

    Use ``to_dict()`` for JSON serialization.
    """

    # --- BAM-level ---
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

    # --- The junction-strand tag on each spliced record (``sj_strand_tag``; read per record) ---
    n_sj_tag_tagged: int = 0  # spliced records given a strand
    n_sj_tag_untagged: int = 0  # spliced records no tag gave a strand
    n_sj_tag_conflicting: int = 0  # spliced records whose XS and ts disagree (ambiguous)
    n_sj_tag_invalid: int = 0  # XS or ts tags holding neither '+' nor '-' (read as absent)

    # --- Hits whose two mates are not one molecule, by reason, and the templates left with none ---
    n_pairs_rejected_missing_mate: int = 0  # one mate only: no fragment length at all
    n_pairs_rejected_reference: int = 0  # the mates on two references
    n_pairs_rejected_strand: int = 0  # on one strand
    n_pairs_rejected_outward: int = (
        0  # facing away from each other, or a mate past the other's 5' end
    )
    n_pairs_rejected_length: int = 0  # no path gives a fragment length in [1, max_frag_length]
    n_withheld_no_valid_pair: int = (
        0  # templates withheld: no hit left (no fragment, not intergenic)
    )

    # --- The per-fragment splice census (scanner QC; ``rigel report``'s splice breakdown) ---
    #
    # ONE observation per placement the scanner hands the accumulator's deposit adapter — every
    # non-chimeric hit, intergenic ones included: once for a unique mapper, once per placement for a
    # multimapping fragment. That population is stated; reading these off a fragment-length category
    # model instead would silently count only the fragments that also yielded a length observation.
    #
    # The field names are derived, not chosen — ``rigel.splice.census_field`` builds each one from
    # its :class:`~rigel.splice.SpliceType` member name, and the C++ builds the same key from
    # ``splice_type_label``. Renaming one of these by hand breaks that correspondence silently,
    # because ``pipeline._apply_scan_stats`` copies them with ``dict.get(key, 0)``.
    n_census_unspliced: int = 0
    n_census_spliced_unannot: int = 0
    n_census_spliced_annot: int = 0
    n_census_spliced_implicit: int = 0
    n_census_splice_artifact: int = 0

    #: Censused but never offered to the accumulator's ``deposit()``: the fragment is not one
    #: molecule on one region_bound axis (chiefly blocks on more than one reference). It exists so that
    #:
    #:   Σ census − n_census_splice_artifact
    #:       == qc.deposited + Σ qc.dropped_* + qc.deferred_undetermined_gap
    #:          + qc.deferred_multiple_placements + n_deposit_not_offered
    #:
    #: closes exactly when no fragment offers more than one placement. The census counts placements
    #: and the ``qc`` counters count fragments, so a multimapping fragment that offers several
    #: placements is censused once per placement and counted once. Those returns were a silent loss
    #: before this counter.
    n_deposit_not_offered: int = 0

    @property
    def n_intergenic(self) -> int:
        """Total intergenic fragments (unspliced + spliced)."""
        return self.n_intergenic_unspliced + self.n_intergenic_spliced

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary.

        Includes the computed ``n_intergenic`` alongside the raw dataclass
        fields.
        """
        d = asdict(self)
        d["n_intergenic"] = self.n_intergenic
        return d
