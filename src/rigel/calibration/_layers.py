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
        # `total_abundance` is the composition-FREE per-slot total (the START/END banks side-selected
        # by the wall rule, plus the exact boundary banks) — geometry work, and it reads the geometry.
        # `structural_claims` is the first pass's substrate — per-slot structural classes derived from
        # the statics alone (sideways of `region_geometry`, whose statics and `g1_locked` it reads).
        ("region_geometry", "simplex_logodds", "total_abundance", "structural_claims"),
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
        # How dense a component is, and every fitted population prior. ⚠ It was SIX modules until
        # 2026-08-21, when `background_reference` was converge-and-deleted: it computed a second pooled
        # intergenic background on the SAME pool as `density_deconv.fit_intron_background` (measured
        # identical, n = 1,298) and no caller ever consumed it. The census prints the live count.
        (
            "run_fill",
            "density_model",
            "density_deconv",
            # `rna_anchor` is the certified-flux stream's arithmetic (owner design 2026-08-24;
            # ruled A MESSAGE 2026-08-25): the sender-side evidence bundle and the recipient-side
            # rows, anchored on certified splice flux + adjacent-intron nascent rates. It reuses
            # `density_deconv`'s NegBinomial SIDEWAYS and reads layer-3 geometry/claims DOWN; the
            # RELAY (layer 6) imports it DOWN to deliver the claim as `PsiMessage.lam_rows`.
            "rna_anchor",
            "landscape",
            # `abundance_landscape` is the pre-pass-0 TOTAL-density field + mode census — it reuses
            # `landscape`'s estimator SIDEWAYS and reads `total_abundance` (layer 3) DOWN.
            "abundance_landscape",
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
            # `messages/foundation` is THE FOUNDATION SPEC (owner architecture, 2026-08-26):
            # the Message type (unspliced + spliced lanes, provenance structural), the
            # propagation-time PropagationModel and the solve-time SolveModel — the contracts
            # every variance model implements by overriding their narrow extension points.
            "messages/foundation",
            "messages/policy",
            # `messages/policy` is THE MESSAGE POLICY skeleton (owner tear-down, 2026-08-27):
            # the foundation spec's (PropagationModel, SolveModel) pair running on the backbone's
            # Policy protocol. The other policies are DONORS scheduled for extract-then-delete.
            "messages/variance",
            "messages/silent",
            "messages/relay",
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
    """The layer of a module, by its package-relative name (``"sweep"``, ``"messages/relay"``).

    ``None`` means UNPLACED, which the gate treats as a failure rather than as a default — a module with no
    declared home is exactly the flat-pile state this file exists to end.
    """
    return _OF.get(module)
