"""⭐⭐⭐ THE LAYERING — **where does a change go?** This file is the answer, and a test enforces it.

       Gate: ``tests/calibration/test_layering.py``
       Census: ``python scripts/design/module_census.py``

⛔⛔ **THE PROBLEM THIS FILE EXISTS FOR.** The package had no stated shape, and the shape is not what a
flat module list suggests: measured from the AST when this file landed (2026-08-07, 35 modules), there were
**no import cycles** and **18 of the 35 modules had exactly one importer**. It was never a knot. It was a
**FLAT PILE of peers** — and a flat pile is the one structure that cannot tell you where to add anything,
because every file is equally plausible. ⚠ The counts move whenever a module is added or retired; the
census re-derives them, and no number here is maintained by hand.

⭐ The layers below were not invented. They were **read off the existing import boundaries**: every one of them
was already there, and the only thing missing was a name and a gate. That is why declaring them costs no
behaviour change — ``test_layering.py`` passed on the tree as it stood.

**THE ONE RULE: an import may point DOWN a layer or SIDEWAYS within one. Never UP.** A module that needs
something from a higher layer is telling you the thing belongs lower, or that your module belongs higher.

Where to put a change
---------------------
=========================================  ==========================================================
if the change is about…                    it goes in layer
=========================================  ==========================================================
what a fragment tally MEANS                1 · the payload view
how many places a fragment COULD have sat  2 · opportunity
one slot's own numbers, and ψ              3 · geometry and the per-slot solve
which strand a fragment came from          4 · strand
how dense a component is, and the priors   5 · density and prior
what one neighbour tells another           6 · the solve  (⭐ and `messages/` inside it)
turning the solve into a result            7 · assemble
=========================================  ==========================================================

⚠ **A layer is not a promise that its modules are the right SIZE.** Layer 4 is five modules for one concept
and layer 5 is six. Whether those should be fewer files is a separate question with its own risk, and this
file deliberately does not answer it — it answers *where*, which is what was blocking a reader. ⛔ The
per-layer line counts that used to be quoted here (939 for layer 4, 1,858 for layer 5, measured 2026-08-07)
had drifted and are dropped rather than re-stated: the census prints them live, which is the only place a
count of a moving thing belongs.
"""

from __future__ import annotations

__all__ = ["LAYERS", "layer_of"]

#: ``(number, title, modules)`` — the layer a module belongs to, lowest first. ⛔ Every module in the
#: package must appear exactly once; the gate fails on an unplaced module, because an unplaced module is a
#: new file nobody decided the home of.
LAYERS: tuple[tuple[int, str, tuple[str, ...]], ...] = (
    (
        0,
        "vocabulary — no calibration deps",
        # The words everything else is written in. ⭐ `signature` is the region bitmask, `region_chain` is the
        # N E N E … N sequence, `errors` is the one exception type. Nothing here knows what a solve is.
        ("errors", "signature", "region_chain", "_layers"),
    ),
    (
        1,
        "the payload view",
        # What a fragment tally MEANS. `splice_graph` is the v8 index; `substrate` and `region_arrays` are
        # the accumulator's banks presented as per-object arrays.
        ("splice_graph", "substrate", "region_arrays"),
    ),
    (
        2,
        "opportunity — how many places a fragment COULD have sat",
        # ⭐ The deposit weight is 1/OPPORTUNITY, so every divisor in the tool is derived here and nowhere
        # else. `fl` is the entry point the scanner and the second pass call.
        (
            "effective_length",
            "capture_eff_length",
            "sj_opportunity",
            "gdna_opportunity",
            "fl",
        ),
    ),
    (
        3,
        "geometry and the per-slot solve",
        # One slot's own numbers, and ψ — the log-density log-odds posterior over (f_pos, f_neg, f_g).
        # ⚠ `simplex_logodds` is ψ and was 784 boundaries when this file landed (2026-08-07); it is the
        # single densest thing in the package. ⛔ Like every other count once quoted here that number has
        # DRIFTED and is not maintained by hand — the census re-derives it.
        ("region_geometry", "simplex_logodds"),
    ),
    (
        4,
        "strand — which strand a fragment came from",
        # ⚠ FIVE modules for one concept. `strand_likelihood` is the TWO-component executable REFERENCE
        # that ψ's three-component form is gated against; the other four are production.
        (
            "strand_likelihood",
            "strand_deconv",
            "strand_balance",
            "strand_summary",
            "gdna_strand",
        ),
    ),
    (
        5,
        "density and prior",
        # How dense a component is, and every fitted population prior. ⚠ SIX modules; 1,858 boundaries when
        # this file landed (2026-08-07) — the same drifted count the module docstring above records, kept
        # here with its date rather than restated as current. The census prints it live.
        (
            "run_fill",
            "density_model",
            "density_deconv",
            "npmle",
            "landscape",
            "background_reference",
        ),
    ),
    (
        6,
        "the solve — what one neighbour tells another",
        # ⭐ The backbone and the message policy. `sweep` owns the shape of the solve and five assertions;
        # `messages/` owns every argument about what a message should say. `DESIGN.md` §6.1.
        (
            "region_init",
            "sweep",
            "messages",
            "messages/variance",
            "messages/silent",
            "messages/head",
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
    """The layer of a module, by its package-relative name (``"sweep"``, ``"messages/head"``).

    ``None`` means UNPLACED, which the gate treats as a failure rather than as a default — a module with no
    declared home is exactly the flat-pile state this file exists to end.
    """
    return _OF.get(module)
