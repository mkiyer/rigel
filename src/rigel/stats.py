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

    # --- Records: each flag counter counts every record carrying the flag ---
    total: int = 0
    qc_fail: int = 0
    unmapped: int = 0
    secondary: int = 0
    supplementary: int = 0
    duplicate: int = 0
    proper_pair: int = 0  # primary records flagged as a proper pair
    filtered: int = 0  # records the scan drops: unmapped, QC-failed, or a skipped duplicate

    # --- One outcome per read pair (a run of records sharing a name), one status per eligible record ---
    #   n_read_pairs = n_excluded + the n_withheld_* + n_excluded_by_option + n_placed_unique + n_placed_multimapping
    #   total = filtered + n_records_used + n_records_unused + n_records_unplaced
    n_read_pairs: int = 0
    n_excluded: int = 0  # every record filtered
    n_excluded_by_option: int = 0  # dropped by --no-include-multimap
    n_placed_unique: int = 0  # one usable hit
    n_placed_multimapping: int = 0  # more than one usable hit
    n_records_used: int = 0  # the alignment of a usable hit's location
    n_records_unused: int = 0  # in a placed read pair, in no usable hit
    n_records_unplaced: int = 0  # in a read pair withheld or excluded by the option
    template_digest: int = (
        0  # every read pair with an eligible record folded into 64 bits, for the annotation pass
    )

    # --- Resolution-level ---
    n_fragments: int = 0
    n_intergenic: int = 0  # templates that reach no locus: unspliced, gDNA by construction
    n_with_exon: int = 0
    n_with_annotated_sj: int = 0
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

    # --- Multimapping fragments with a hit in a transcript, and those hits ---
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

    # --- NM and NH on each record: absent, or not an integer of 0..2^31-1 (read as absent: NM 0, NH 1) ---
    n_nm_tag_missing: int = 0
    n_nm_tag_invalid: int = 0
    n_nm_conflicting: int = 0  # alignments whose records disagree on NM (each keeps the larger)
    n_nh_tag_missing: int = 0
    n_nh_tag_invalid: int = 0

    # --- Hits that are not usable, by reason, and the templates left with none (no fragment, not intergenic) ---
    n_hits_rejected_missing_mate: int = 0  # one mate only: no fragment length at all
    n_hits_rejected_reference: int = 0  # the mates on two references
    n_hits_rejected_strand: int = 0  # on one strand
    n_hits_rejected_outward: int = 0  # their 5' ends facing away from each other
    n_hits_rejected_length: int = 0  # no path gives a fragment length in [1, max_frag_length]
    n_hits_rejected_contradictory: int = (
        0  # a molecule whose junctions overlap unequal or name both strands
    )
    n_hits_rejected_novel_junction: int = (
        0  # a molecule with a junction no annotated transcript holds
    )
    n_hits_rejected_unsupported: int = 0  # junctions no one transcript over both mates holds
    n_hits_rejected_off_index: int = 0  # both mates on one reference the index does not hold
    n_hits_rejected_chimeric: int = 0  # a mate aligned in pieces (SA, or a 5' hard clip)
    n_withheld_unpaired: int = (
        0  # templates none of whose hits is a molecule, a template with no hit included
    )
    # ... by the furthest check any of its hits reached
    n_withheld_unpaired_missing_mate: int = 0
    n_withheld_unpaired_reference: int = 0
    n_withheld_unpaired_strand: int = 0
    n_withheld_unpaired_outward: int = 0
    n_withheld_unpaired_length: int = 0
    n_withheld_unannotated_splice: int = (
        0  # templates whose molecules no annotated transcript holds
    )
    n_withheld_off_index: int = 0  # templates with a hit on a reference outside the index
    n_withheld_chimeric: int = 0  # templates with a chimeric hit, none usable
    n_nh_above_reported: int = 0  # NH above a mate's reported records (a truncated output)

    # --- The per-fragment splice census (scanner QC; ``rigel report``'s splice breakdown) ---
    #
    # ONE observation per placement the scanner hands the accumulator's deposit adapter — every
    # resolved hit, intergenic ones included: once for a unique mapper, once per placement for a
    # multimapping fragment. That population is stated; reading these off a fragment-length category
    # model instead would silently count only the fragments that also yielded a length observation.
    #
    # The field names are derived, not chosen — ``rigel.splice.census_field`` builds each one from
    # its :class:`~rigel.splice.SpliceType` member name, and the C++ builds the same key from
    # ``splice_type_label``. ``pipeline._apply_scan_stats`` refuses a scan key with no field, so renaming
    # one of these by hand fails the scan rather than dropping the count.
    n_census_unspliced: int = 0
    n_census_spliced_annot: int = 0
    n_census_spliced_implicit: int = 0
    n_census_splice_artifact: int = 0

    #: Censused but never offered to the accumulator's ``deposit()``: no accumulator was installed.
    #: It exists so that
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

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary."""
        return asdict(self)
