"""The layering of the calibration package: the declaration of where a change goes.

An import may point down a layer or sideways within one, never up. A module that needs something from a higher layer
is saying that the thing belongs lower, or that the module itself belongs higher; a type almost always belongs lower.
``LAYERS`` is that declaration, lowest layer first, and every module in the package appears in it exactly once.

=========================================  ==========================================================
if the change is about…                    it goes in layer
=========================================  ==========================================================
what a fragment tally means                1 · the payload view
how many places a fragment could have sat  2 · opportunity
one slot's own numbers, and psi            3 · geometry and the per-slot solve
which strand a fragment came from          4 · strand
how dense a component is, and the priors   5 · density and prior
what one neighbour tells another           6 · the solve (and ``messages/`` inside it)
turning the solve into a result            7 · assemble
=========================================  ==========================================================

Gate: ``tests/calibration/test_layering.py`` checks the direction rule against the real imports and refuses an
unplaced module.
"""

from __future__ import annotations

__all__ = ["LAYERS", "layer_of"]

#: ``(number, title, modules)``: the layer each module belongs to, lowest first.
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
        # What a fragment tally means. `splice_graph` is the index's; `substrate` and `region_arrays` present the
        # accumulator's banks as per-object arrays.
        ("splice_graph", "substrate", "region_arrays"),
    ),
    (
        2,
        "opportunity — how many places a fragment COULD have sat",
        # The deposit weight is 1/opportunity, so every divisor in the tool is derived here. `gdna_density` is the gDNA
        # background rate, a count over an opportunity: the contamination-robust one-sided rate, reported beside the
        # naive pooled one. `shared_length_law` decides the EM's two length tables from the two laws `fl` measures.
        (
            "effective_length",
            "capture_eff_length",
            "sj_opportunity",
            "gdna_opportunity",
            "gdna_density",
            "fl",
            "shared_length_law",
        ),
    ),
    (
        3,
        "geometry and the per-slot solve",
        # One slot's own numbers, and psi: the log-density log-odds posterior over (f_pos, f_neg, f_g), which
        # `simplex_logodds` owns.
        ("region_geometry", "simplex_logodds"),
    ),
    (
        4,
        "strand — which strand a fragment came from",
        # The library's sense fraction and the strand protocol decision read from it.
        ("strand_balance",),
    ),
    (
        5,
        "density and prior",
        # How dense a component is, and every fitted population prior.
        (
            "landscape",
            # `capture_efficiency` is every region's and boundary's capture efficiency: a posterior under the
            # landscape prior from the object's own deconvolved gDNA count, read by the ruler and the locus prior.
            "capture_efficiency",
            # `intron_departure` is the introns' departure prior ε for the library, from single-strand introns'
            # antisense counts under the background law (2) and the library's sense fraction (4).
            "intron_departure",
        ),
    ),
    (
        6,
        "the solve — what one neighbour tells another",
        # The backbone and the message policy. `sweep` owns the shape of the solve and its assertions and runs it
        # in one native call; `messages/` owns what a message says.
        (
            # `sweep` is the backbone;`diagnostics` is the instruments' record of a solve, which production never passes
            "diagnostics",
            "sweep",
            "messages",
            # `messages/__init__` is the policy interface (a name, a strand model, a library), the chain view and the
            # per-slot arrays it is built from;
            # `messages/silent` is the floor every policy is judged against; `messages/transfer` is the shipped
            # composition-transfer policy, whose builders, passes and solve are the kernel's
            # (`native/transfer_kernel.h`).
            "messages/silent",
            "messages/transfer",
        ),
    ),
    (
        7,
        "assemble — turning the solve into a result",
        (
            "priors",
            "result",
            "track",
            "calibrate",
            "<pkg>",
        ),
    ),
)

_OF: dict[str, int] = {m: num for num, _t, members in LAYERS for m in members}


def layer_of(module: str) -> int | None:
    """The layer of a module by its package-relative name (``"sweep"``, ``"messages/transfer"``), or ``None`` if
    unplaced, which the gate treats as a failure.
    """
    return _OF.get(module)
