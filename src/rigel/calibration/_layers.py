"""The layering of the calibration package — the declaration of where a change goes.

The one rule: an import may point DOWN a layer or SIDEWAYS within one, never UP. A module that
needs something from a higher layer is telling you the thing belongs lower, or that the module
itself belongs higher. A type almost always belongs lower.

``LAYERS`` below is that declaration, lowest layer first, and it is authoritative: every module in
the package must appear in it exactly once. ``tests/calibration/test_layering.py`` enforces both
halves — the direction rule against the real imports, and the requirement that nothing is unplaced,
since an unplaced module is a file whose home nobody decided. The graph is re-derived from the imports,
so no count belongs here.

Where to put a change
---------------------
=========================================  ==========================================================
if the change is about…                    it goes in layer
=========================================  ==========================================================
what a fragment tally MEANS                1 · the payload view
how many places a fragment COULD have sat  2 · opportunity
one slot's own numbers, and psi            3 · geometry and the per-slot solve
which strand a fragment came from          4 · strand
how dense a component is, and the priors   5 · density and prior
what one neighbour tells another           6 · the solve (and ``messages/`` inside it)
turning the solve into a result            7 · assemble
=========================================  ==========================================================

A layer says where, not how big. It is not a promise that its modules are the right size or the
right count; whether a layer should be fewer files is a separate question this file does not answer.
"""

from __future__ import annotations

__all__ = ["LAYERS", "layer_of"]

#: ``(number, title, modules)`` — the layer a module belongs to, lowest first. Every module in the
#: package must appear exactly once; the gate fails on an unplaced module.
LAYERS: tuple[tuple[int, str, tuple[str, ...]], ...] = (
    (
        0,
        "vocabulary — no calibration deps",
        # The words everything else is written in: `signature` is the region bitmask, `region_chain` is
        # the N E N E … N sequence, `errors` the exception types. Nothing here knows what a solve is.
        ("errors", "signature", "region_chain", "_layers"),
    ),
    (
        1,
        "the payload view",
        # What a fragment tally MEANS. `splice_graph` is the index's; `substrate` and `region_arrays` are
        # the accumulator's banks presented as per-object arrays.
        ("splice_graph", "substrate", "region_arrays"),
    ),
    (
        2,
        "opportunity — how many places a fragment COULD have sat",
        # The deposit weight is 1/opportunity, so every divisor in the tool is derived here and nowhere
        # else. `fl` is the entry point the scanner and the second pass call.
        # `gdna_density` is the gDNA background RATE — a count divided by an opportunity, which is this
        # layer's job. It owns BOTH estimators of that one quantity (the naive pooled rate and the
        # contamination-robust one-sided rate), so layer 5's `density_deconv` calls DOWN to it rather
        # than carrying a second implementation.
        (
            "effective_length",
            "capture_eff_length",
            "sj_opportunity",
            "gdna_opportunity",
            "gdna_density",
            "fl",
        ),
    ),
    (
        3,
        "geometry and the per-slot solve",
        # One slot's own numbers, and psi — the log-density log-odds posterior over (f_pos, f_neg, f_g),
        # which `simplex_logodds` owns and which is the single densest thing in the package.
        # `total_abundance` is the composition-FREE region count and exposure (the START/END banks
        # side-selected by the wall rule) — geometry work, and it reads the geometry.
        ("region_geometry", "simplex_logodds", "total_abundance"),
    ),
    (
        4,
        "strand — which strand a fragment came from",
        # `strand_summary` is the dependency-light QC view the pipeline reads without importing
        # calibration; the other two are production.
        ("strand_balance", "strand_summary", "gdna_strand"),
    ),
    (
        5,
        "density and prior",
        # How dense a component is, and every fitted population prior.
        (
            "density_model",
            "density_deconv",
            "landscape",
            # `capture_efficiency` is the per-piece capture efficiency: a posterior under the landscape
            # prior from the deconvolved gDNA, read by the ruler and the locus prior through the result.
            "capture_efficiency",
            # `abundance_landscape` is the pre-pass-0 TOTAL-density field + mode census — it reuses
            # `landscape`'s estimator sideways and reads `total_abundance` (layer 3) down.
            "abundance_landscape",
        ),
    ),
    (
        6,
        "the solve — what one neighbour tells another",
        # The backbone and the message policy. `sweep` owns the shape of the solve and its assertions and
        # runs it in one native call; `messages/` owns every argument about what a message should say.
        (
            # `region_init` is the strand protocol decision and the own-evidence predicate
            "region_init",
            # `sweep` is the backbone; `blocks` is the diagnostics capture and the chain's view;
            "blocks",
            "sweep",
            "messages",
            # `messages/__init__` is the policy's interface (a name, a strand model, a library) and the
            # chain view; `messages/silent` is the measured floor every policy is judged against;
            # `messages/transfer` is the shipped composition-transfer policy's name, strand model and
            # library — its builders, passes and solve are the kernel's (`native/transfer_kernel.h`).
            "messages/silent",
            "messages/transfer",
        ),
    ),
    (
        7,
        "assemble — turning the solve into a result",
        (
            "derive",
            "priors",
            "result",
            "diagnostics",
            "track",
            "calibrate",
            "<pkg>",
        ),
    ),
)

_OF: dict[str, int] = {m: num for num, _t, members in LAYERS for m in members}


def layer_of(module: str) -> int | None:
    """The layer of a module, by its package-relative name (``"sweep"``, ``"messages/transfer"``).

    ``None`` means unplaced, which the gate treats as a failure rather than as a default: a module with
    no declared home is exactly the state this file exists to end.
    """
    return _OF.get(module)
