"""THE BACKBONE — the self-solve, two directional passes, one ψ solve, one write-back, four assertions.

       Gate: ``tests/calibration/test_sweep_backbone.py``

Each slot's unspliced fragment mass is deconvolved into a pie ``(f_pos, f_neg, f_g)`` — sense-RNA /
antisense-RNA / gDNA — over the ``N E N E … N`` chain (`region_chain`), on the TWO-PHASE shape: every
node's own claim (`prepare`), a forward pass then a backward pass in which each RECIPIENT receives what
its neighbour sends (`propagate`), and ONE solve per node from its own evidence, the two held messages
and the prior (`solve`). The chain is a forest of linear paths, so that is exact belief propagation, not
an iteration.

This file knows nothing about capture, splice in, levels, lanes or enrichment — those words do not appear
in it, and that is the design rather than tidiness. Everything about *what a message says* is a
:mod:`~.messages` policy. What is left here is the shape of the solve and the four invariants no policy may
break:

===================================================  ====================================================
the backbone asserts                                 it would have caught
===================================================  ====================================================
the kernel sees only the two NEIGHBOUR states        TRAPS: a-message-from-the-destinations-belief
every delivered row is one row per slot, finite      a row array off the solve grid, or a NaN reaching ψ
``|T| <= 3``                                         AXIOM 0, made executable
the write-back touches only ``solvable`` slots       a replay that read the untouched mask as a difference
===================================================  ====================================================

The assertions live HERE, not in the policy, and that is the entire point: a future policy can be as
wrong as it likes and still cannot commit any of them, and each has shipped at least once.

And ONE structural rule, which is what lets the chain be solved a locus at a time: a TERMINAL — a REGION
that admits no RNA strand, so it is structurally pure gDNA, solved and fixed before any message exists —
RECEIVES NOTHING. The pass never asks the kernel for a hop into one; the terminal holds silence from a side
it has a neighbour on. Nothing can therefore cross a terminal, and the chain breaks into independent loci
at every one (`region_chain.locus_blocks`). Measured on the human chain before it was made structural: of
1,206,202 composition faces and 4,621,302 lane faces the shipped policy built, none delivered into a
terminal, so the rule moved no number; it is written so that no future policy can move the boundary
condition a block solve depends on.

Two gates on this slot in the pipeline
--------------------------------------
Both come from the region SIGNATURE and never from the counts:

* the SOLVE gate (``solvable``) — a slot deconvolves its own split iff it admits >= 1 RNA strand and has
  unspliced mass. A slot with no admissible RNA strand (an intergenic region, or a gene boundary) is a
  LOCKED all-gDNA object; it is not solved and keeps its signature-binary init, because RNA cannot cross a
  gene boundary so its unspliced mass is purely gDNA.
* the EMISSION gate — which MESSAGES a slot sends. That is a policy question and lives there.
"""

from __future__ import annotations

import numpy as np

from .blocks import block_slice, gather, view_fields
from .message_cache import MessageCache
from .messages import BlockContext, ChainView, PsiMessage, Received
from .messages.silent import SilentPolicy
from .region_geometry import (
    RegionBelief,
    RegionGeometry,
    RegionStatics,
    g1_locked,
    region_gdna_geometry,
    region_rna_geometry,
)
from .region_init import (
    build_region_init,
    has_own_composition_evidence,
    strand_discriminability,
)
from .signature import BIT_EXON_NEG, BIT_EXON_POS, coarse_type_array
from .simplex_logodds import CompositionPriors, _logodds_grid, _solve_regions_logodds_all
from .region_chain import BOUNDARY, REGION, RegionChain, RegionDeconv, locus_blocks

__all__ = [
    "AssertionCounts",
    "chain_boundary_deconv",
    "chain_region_deconv",
    "solve_chain",
]


# ⛔ ASSERTIONS A SHIPPED POLICY IS KNOWN TO VIOLATE, each entered with the measurement that proved it.
# An entry is COUNTED and PUBLISHED rather than raised, because widening an assertion to fit a defect is
# how a gate becomes vacuous (TRAPS: perturb-every-gate / TRAPS: a-gate-that-reconstructs); an entry also
# carries a STRICT xfail in the gate file, this project's convention for a PROVEN defect whose fix is
# panel-negative on its own. The dict is EMPTY, so anything violated raises.
#: ``name -> why it is not fatal yet``.
_KNOWN_VIOLATIONS: dict[str, str] = {}


class AssertionCounts(dict):
    """How many slots violated each backbone assertion, published into the diagnostics capture.

    A count rather than a bool, because TRAPS: could-the-arm-have-fired is the rule: before believing "the
    arm changed nothing", check it COULD have changed something. An assertion reporting 0 violations on a substrate
    where the predicate can never fire is not evidence, so the report also carries how many slots were
    ELIGIBLE for each check.
    """

    def note(self, name: str, violated, eligible) -> None:
        n_v = int(np.count_nonzero(violated))
        n_e = int(np.count_nonzero(eligible)) if eligible is not None else int(np.size(violated))
        self._add(name, n_v, n_e)
        if n_v and name not in _KNOWN_VIOLATIONS:
            raise AssertionError(
                f"backbone assertion {name!r} violated at {n_v:,} of {n_e:,} eligible slots. "
                f"The policy may not do this — see rigel.calibration.sweep's module docstring for what "
                f"each assertion catches. If this is a NEW and PROVEN defect whose fix is panel-negative "
                f"alone, add it to _KNOWN_VIOLATIONS with the measurement, never widen the predicate."
            )

    def _add(self, name: str, n_v: int, n_e: int) -> None:
        cur = self.get(name, {"violations": 0, "eligible": 0})
        self[name] = {"violations": cur["violations"] + n_v, "eligible": cur["eligible"] + n_e}

    def absorb(self, other: "AssertionCounts") -> None:
        """One block's counts into the chain's: the sweep checks every block as it solves it and
        publishes the sums, so the report reads as it did when the chain was one block."""
        for name, c in other.items():
            self._add(name, int(c["violations"]), int(c["eligible"]))


def _check_message(
    msg: PsiMessage, ctx: BlockContext, counts: AssertionCounts, n_owned: int | None = None
) -> None:
    """The assertions on what the policy actually delivered: the population axiom, and every row
    channel one row per slot on the solve grid and finite.

    TRAPS: a-message-from-the-destinations-belief is not checked here because it is enforced BY
    CONSTRUCTION: the propagate kernel is called with two INDICES and builds the message into the
    destination from the source's claim and what the source holds; the backbone writes ``held`` and the
    policy never reaches past its hop. A structural impossibility beats a check. The write-back
    assertion is checked at the write-back, where its basis lives.

    ``n_owned`` restricts the COUNTS to the block's own slots (its first ``n_owned``); the shape checks
    still hold over every slot of the context, since the policy delivers one row per slot it saw.
    """
    n = int(ctx.n_slots) if n_owned is None else int(n_owned)
    # ── AXIOM 0, made executable: |T(slot)| = 1 + free_pos + free_neg, and it is <= 3 ALWAYS ───────────
    # There are THREE populations and there is no fourth. This is a function of TWO BITS, which is what
    # makes it structural rather than something to remember — and the message packet carries exactly three
    # component channels (gDNA, RNA+, RNA-) for the same reason.
    pop = ctx.population_size()[:n]
    counts.note("population_at_most_three", pop > 3, np.ones_like(pop, bool))
    # ── the λ rows: one row per slot in psi's general evidence currency, or absent ───────────────────
    if msg.lam_rows is not None:
        rows = np.asarray(msg.lam_rows)
        if rows.shape[0] != ctx.n_slots or rows.ndim != 2:
            raise ValueError(
                f"lam_rows has shape {rows.shape}; expected ({ctx.n_slots}, K) — a policy must "
                "deliver one row per slot on the solve grid"
            )
        counts.note("lam_rows_finite", ~np.isfinite(rows[:n]).all(axis=1), np.ones(n, bool))
    # ── the cube channel: a (K, K_t) row per AMBIG slot, or absent ──────────────────────────────────
    if msg.cube_rows is not None:
        amb = np.asarray(ctx.free_pos, bool) & np.asarray(ctx.free_neg, bool)
        for slot, row in msg.cube_rows.items():
            row = np.asarray(row)
            if not (0 <= int(slot) < ctx.n_slots) or not amb[int(slot)]:
                raise ValueError(
                    f"cube_rows carries slot {slot}, which is not an AMBIG slot — a cube exists only "
                    "where both strands are live"
                )
            if row.ndim != 2 or row.shape[0] != int(ctx.n_grid):
                raise ValueError(
                    f"cube_rows[{slot}] has shape {row.shape}; expected ({ctx.n_grid}, K_t) — one row "
                    "over the (λ, θ) cube on the solve grid"
                )
        bad = np.array(
            [not np.isfinite(np.asarray(r)).all() for k, r in msg.cube_rows.items() if int(k) < n],
            bool,
        )
        counts.note("cube_rows_finite", bad, np.ones(bad.shape[0], bool))
    # ⛔ TRAPS: could-the-arm-have-fired's anti-degeneracy clause, the half that makes the gate mean
    # anything: on a chain
    # where NO slot admits both RNA strands, ``|T| <= 3`` is satisfied by a substrate that never had a
    # three-population slot to test. That is not the axiom holding, it is the check never running — so the
    # eligible set is the slots that actually reach 3, and a substrate with none of them says so.
    counts.note("population_reaches_three", np.array([]), pop >= 3)


def solve_chain(
    chain: RegionChain,
    statics: RegionStatics,
    geometry: RegionGeometry,
    belief: RegionBelief,
    region_arrays,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_gdna_obs: float = 0.0,
    n_rna_obs: float = 0.0,
    n_grid: int,
    logodds_window: float = 10.0,
    n_tilt: int | None = None,
    n_grid_ss: int | None = None,
    gdna_prior=None,
    rna_prior=None,
    intron_prior=None,
    policy=None,
    block_slots: int | None = None,
    message_cache: "MessageCache | None" = None,
    _capture: dict | None = None,
) -> RegionBelief:
    """One forward-backward sweep over the chain. Returns the resolved :class:`RegionBelief`.

    ``message_cache`` shares the message layer's output between sweeps whose inputs to it are identical
    (:class:`MessageCache`): a block whose digest is held is served its messages and pays only its two ψ
    solves. ``None`` runs the layer for every block, as does a diagnostic capture.

    THE SWEEP IS SOLVED A LOCUS BLOCK AT A TIME. The chain breaks at every TERMINAL — a region that
    admits no RNA strand, structurally pure gDNA, which receives nothing (the module docstring) — into
    loci, and `region_chain.locus_blocks` merges consecutive loci into blocks of up to ``block_slots``
    slots (``None``: the whole chain as one block). Each block is solved on its own slice of every
    input — its own claims, the policy's claims and rules, the two passes, the solve and the
    write-back — reading one slot beyond its own range where that slot is the terminal its last node
    receives from. The only information that crosses a block boundary is the policy's LIBRARY
    (:meth:`~.messages.Policy.library`), reduced once over the whole chain from observations and
    geometry, and the answer is the same for every ``block_slots`` because ψ's read-out is
    chunk-exact (`simplex_logodds`); the block size sets only the working set — a block's
    ``(slots, K)`` arrays instead of the chain's — against the per-block overhead.

    ``policy`` is the message-composition policy (:mod:`~.messages`). It defaults to
    :class:`~.messages.silent.SilentPolicy`, which sends nothing — so a reader of this file plus five
    boundaries holds the whole working system. The shipped answer is
    :class:`~.messages.transfer.TransferPolicy`, which ``calibrate`` passes explicitly.

    ``gdna_prior=None`` is a first-class PRIOR-FREE solve: ψ then carries the Jeffreys reference
    measure alone on both arms. Prior-free is not reference-free. That pass's only job is to be a
    training substrate for the population gDNA hyperprior — it is not the deliverable, and it does not
    have to answer objects it cannot solve.

    ⛔ ψ carries NO reference location: the reference is the symmetric Jeffreys measure and asserts
    nothing; background information enters as the ``intron_prior`` λ-factor, a likelihood whose
    precision scales with counts.
    """
    policy = policy if policy is not None else SilentPolicy()
    n = int(chain.n_slots)
    fp, fn = np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool)
    is_region = np.asarray(chain.kind) == REGION
    # the structural classes, from the signature and once for the chain: the EXON region (the SPLICE
    # IN's destination class, a policy input rather than a gate) and the signature's two EXON bits per
    # slot — a strand's opportunity geometry, which the RNA level lanes read to tell a strand's intron
    # from its exon inside an overlapping locus
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    ri = np.clip(np.asarray(chain.obj_idx, dtype=np.int64), 0, rtype.shape[0] - 1)
    is_exon_region = is_region & (rtype[ri] == 2)
    sig = np.asarray(region_arrays.signature).astype(np.int64)[ri]
    exon_pos = is_region & ((sig & BIT_EXON_POS) > 0)
    exon_neg = is_region & ((sig & BIT_EXON_NEG) > 0)
    # a TERMINAL receives nothing (the structural rule in the module docstring): a region with no
    # admissible RNA strand, the same predicate the SOLVE gate locks — and where the chain breaks
    terminal = is_region & g1_locked(fp, fn)

    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    scalars = dict(
        n_grid=int(n_grid),
        logodds_window=float(logodds_window),
        n_tilt=None if n_tilt is None else int(n_tilt),
        strand_live=strand_discriminability(kappa, od_g, od_r, n_gdna_obs, n_rna_obs) > 0.0,
    )
    # THE LIBRARY — the policy's reductions over the WHOLE chain, once, from observations and geometry
    # alone: the only information a message may carry across a locus boundary
    library = policy.library(
        ChainView(
            **view_fields(chain, statics, geometry, is_exon_region, exon_pos, exon_neg),
            factory_rows=intron_prior,
            **scalars,
        )
    )

    out = {
        k: np.asarray(getattr(belief, k), dtype=np.float64).copy()
        for k in ("f_pos", "f_neg", "f_g", "var_pos", "var_neg", "var_gdna")
    }
    has_composition = np.zeros(n, dtype=bool)
    counts = AssertionCounts()
    diagnostics: list = []
    rowless = 0  # slots owned by blocks whose policy delivered no row array
    for b in locus_blocks(chain, terminal, block_slots):
        sl = slice(b.start, b.end)
        res = _solve_block(
            block_slice(chain, sl),
            block_slice(statics, sl),
            block_slice(geometry, sl),
            block_slice(belief, sl),
            is_exon_region[sl],
            exon_pos[sl],
            exon_neg[sl],
            terminal[sl],
            n_owned=b.stop - b.start,
            kappa=kappa,
            od_g=od_g,
            od_r=od_r,
            n_gdna_obs=n_gdna_obs,
            n_rna_obs=n_rna_obs,
            n_grid_ss=n_grid_ss,
            gdna_prior=gdna_prior,
            rna_prior=rna_prior,
            factory_rows=None if intron_prior is None else intron_prior[sl],
            policy=policy,
            library=library,
            scalars=scalars,
            cache=None if _capture is not None else message_cache,
            _capture={} if _capture is not None else None,
        )
        own = slice(b.start, b.stop)
        for k, arr in res["belief"].items():
            out[k][own] = arr[: b.stop - b.start]
        has_composition[own] = res["has_composition"][: b.stop - b.start]
        counts.absorb(res["counts"])
        if "lam_rows_finite" not in res["counts"]:
            rowless += b.stop - b.start
        if _capture is not None:
            diagnostics.append((b, res["diagnostics"]))
    # the row check's ELIGIBLE set reads as the chain's: where any block delivered rows, the chain's
    # row array covers every slot (a silent block's rows are zero rows, `blocks.gather` fills them so), and
    # a zero row is a finite row that was checked — so the published count does not depend on how the
    # chain was cut. Where no block delivered, the check never ran and the key stays absent.
    if rowless and "lam_rows_finite" in counts:
        counts._add("lam_rows_finite", 0, rowless)

    if _capture is not None:  # inert diagnostic hook
        _capture.update(gather(diagnostics, n))
        _capture.update(
            backbone_assertions=counts,
            order=np.arange(n, dtype=np.int64),
            left=np.asarray(chain.left, np.int64),
            right=np.asarray(chain.right, np.int64),
            # which policy ran, read off the artifact — the witness an instrument's "the arm ran"
            # assertion needs (TRAPS: an-ablation-that-never-ran), never a config flag it did not thread
            policy_name=str(getattr(policy, "name", type(policy).__name__)),
            solve_grid=_logodds_grid(int(n_grid), float(logodds_window))[1],
            intron_prior=intron_prior,
        )

    return RegionBelief(**out, has_composition=has_composition)


def _solve_block(
    chain: RegionChain,
    statics: RegionStatics,
    geometry: RegionGeometry,
    belief: RegionBelief,
    is_exon_region,
    exon_pos,
    exon_neg,
    terminal,
    *,
    n_owned: int,
    kappa: float,
    od_g: float,
    od_r: float,
    n_gdna_obs: float,
    n_rna_obs: float,
    n_grid_ss: int | None,
    gdna_prior,
    rna_prior,
    factory_rows,
    policy,
    library,
    scalars: dict,
    cache: "MessageCache | None",
    _capture: dict | None,
) -> dict:
    """One block of the chain, solved end to end on its own slice of every input: the self-solve, the
    policy's claims and rules, the two passes, the solve, the write-back — today's whole sweep, on
    ``chain.n_slots`` slots of which the first ``n_owned`` are the block's own (the rest is the terminal
    it reads). Returns the block's belief arrays, its ``has_composition`` predicate, its assertion counts and
    its diagnostic capture (``None`` unless asked for), each over every slot of the block; the caller
    keeps the owned prefix."""
    n_grid, logodds_window, n_tilt = scalars["n_grid"], scalars["logodds_window"], scalars["n_tilt"]
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fp, fn = statics.free_pos, statics.free_neg
    f_pos = np.asarray(belief.f_pos, dtype=np.float64).copy()
    f_neg = np.asarray(belief.f_neg, dtype=np.float64).copy()
    f_g = np.asarray(belief.f_g, dtype=np.float64).copy()
    var_pos = np.asarray(belief.var_pos, dtype=np.float64).copy()
    var_neg = np.asarray(belief.var_neg, dtype=np.float64).copy()
    var_g = np.asarray(belief.var_gdna, dtype=np.float64).copy()
    # the INCOMING belief, kept for the diagnostic capture: it is the ``fg_ref`` the final solve freezes
    # its variance at, so a channel-ablation replay must pass the SAME reference to be faithful.
    _fg_init, _fp_init, _fn_init = f_g.copy(), f_pos.copy(), f_neg.copy()

    fields = view_fields(chain, statics, geometry, is_exon_region, exon_pos, exon_neg)
    # per-slot "global" gDNA support — the basis the rate prior is fit and projected on.
    mass_global, eff_global = region_gdna_geometry(geometry)

    _, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))

    def _psi(
        g_arr, msg: PsiMessage, *, fg_ref, fpos_ref, fneg_ref, extra_lam_rows=None, cube_rows=None
    ):
        """The per-slot solve (the log-density log-odds backend). Phase A calls it with a silent message;
        the final call passes the policy's rows (the λ rows and the cube rows).

        ``fg_ref`` is the count-zero-information variance freeze: the reference is the incoming belief,
        so the variance — hence the message precision — is evaluated near the truth and not at a flat 1/2.
        It is passed EXPLICITLY rather than closed over, because the write-back rebinds the belief and one
        diagnostic solve below deliberately runs after that.
        """
        return _solve_regions_logodds_all(
            u_pos,
            u_neg,
            fp,
            fn,
            n_slot,
            spliced_slot,
            kappa=kappa,
            od_g=od_g,
            od_r=od_r,
            n_grid=int(n_grid),
            L=float(logodds_window),
            n_tilt=n_tilt,
            n_grid_ss=n_grid_ss,
            priors=g_arr,
            # the gDNA intron-factory λ-factor (anchored, per-intron, 0 elsewhere): deconvolves confident gDNA
            # from introns against the intergenic background BEFORE the sweep resolves the pie. Added to ψ,
            # distinct from the gDNA arm; participates in the local solve AND the message layer.
            lam_logprior=(
                factory_rows
                if extra_lam_rows is None
                else (
                    np.asarray(extra_lam_rows, np.float64)
                    if factory_rows is None
                    else factory_rows + np.asarray(extra_lam_rows, np.float64)
                )
            ),
            fg_ref=fg_ref,
            fpos_ref=fpos_ref,
            fneg_ref=fneg_ref,
            # the CUBE channel: the RNA level lanes delivered at AMBIG slots, final solve only
            cube_rows=cube_rows,
        )

    # THE gDNA ARM of ψ — the COMPOSITION prior, and ONLY that. A total-density model is an ENRICHMENT
    # model, not a DNA composition prior: letting it vote a slot's f_g is the count-votes-composition
    # regression.
    # ONE construction site for ψ's composition arms. The RNA member stays ``None`` until something
    # fits an RNA landscape; ``None`` there means "that arm takes its derived reference", which is the
    # shipped behaviour and a first-class configuration rather than a gap.
    # The RNA arm asks the SAME landscape about the OTHER component: the complementary fraction
    # `1 - f_g` against RNA's own opportunity. `mass_global` is shared because both components split one
    # unspliced population (`region_rna_geometry`).
    _rna_mass, _eff_rna = region_rna_geometry(geometry)
    global_lp = CompositionPriors(
        gdna=gdna_prior.logprior(solve_grid, mass_global, eff_global)
        if gdna_prior is not None
        else None,
        rna=rna_prior.logprior(1.0 - solve_grid, _rna_mass, _eff_rna)
        if rna_prior is not None
        else None,
    )

    # Slot ids ARE the genomic visiting order, so the order is ``arange`` and the chain does not store
    # it. The scans are sequential, so iterate as a Python list of ints.
    order_list = list(range(int(chain.n_slots)))
    # ── the per-slot message-free SELF-SOLVE ──────────────────────────────────────────────────────────
    own = build_region_init(
        statics,
        geometry,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        n_gdna_obs=n_gdna_obs,
        n_rna_obs=n_rna_obs,
        n_grid=int(n_grid),
        logodds_window=float(logodds_window),
        n_tilt=n_tilt,
        n_grid_ss=n_grid_ss,
        belief=belief,
        priors=global_lp,
        intron_prior=factory_rows,
    )

    has_own_composition = np.asarray(own.tau_lam, np.float64) > 0.0
    ctx = BlockContext(
        **fields,
        **scalars,
        # the intron factory's rows are an observation on the context: the one array that is both
        # ψ's λ-factor and the intron's own claim
        factory_rows=factory_rows,
        # beliefs — SOURCE-SIDE ONLY (TRAPS: a-message-from-the-destinations-belief)
        has_own_composition=has_own_composition,
        belief_fg=f_g,
    )
    CNT, n_slot, spliced_slot = ctx.unspliced_count, ctx.n_slot, ctx.spliced_slot
    u_pos, u_neg = CNT[:, 0], CNT[:, 1]
    # ── the SOLVE gate. Structural, from the signature, never from the counts ──────────────────────────
    solvable = (fp | fn) & (n_slot > 0.0)

    # ── THE MESSAGE LAYER: served from the cache where its every input is unchanged, else run ─────────
    key = None if cache is None else cache.key(ctx, library, policy)
    entry = None if key is None else cache.get(key)
    if entry is not None:
        msg = cache.message(entry)
        from_left = from_right = None
        held_composition = entry[2].copy()
        counts = AssertionCounts(entry[3])
    else:
        prepared = policy.prepare(ctx, library)

        # ── PHASE 1, PROPAGATE: the FORWARD pass L→R and the BACKWARD pass R→L ───────────────────────
        # ⛔ ONE pass each, in chain order, which on a chain IS forward-backward. It is not an iterative
        # scheme (TRAPS: a-comment-quoted-as-a-finding). When both passes end every node holds one
        # message from each neighbour it has: the recipient's kernel wrote its row, or silence stands.
        term = np.asarray(terminal, bool).tolist()
        from_left = _pass(
            order_list, left.tolist(), prepared, n_grid, backward=False, terminal=term
        )
        from_right = _pass(
            order_list[::-1], right.tolist(), prepared, n_grid, backward=True, terminal=term
        )

        # ── PHASE 2, SOLVE: the policy's half — the two held messages into ψ's channels ──────────────
        msg = prepared.solve(from_left, from_right)
        counts = AssertionCounts()
        _check_message(msg, ctx, counts, n_owned)
        # a COMPOSITION row received from a neighbour, read off the two tables (see the has_composition
        # predicate below for why not `msg.lam_rows`)
        held_composition = from_left.has_composition | from_right.has_composition
        if key is not None:
            cache.put(key, msg, held_composition, counts)

    # THE CITIZENSHIP SEAM: a delivered claim joins the λ-factor rows in the FINAL solve only. Phase-A
    # (`build_region_init`) and the own-evidence precision never see it — an imputation may inform the
    # fused answer, never masquerade as the slot's own evidence.
    dc_fin = _psi(
        global_lp,
        msg,
        fg_ref=f_g,
        fpos_ref=f_pos,
        fneg_ref=f_neg,
        extra_lam_rows=msg.lam_rows,
        cube_rows=msg.cube_rows,
    )

    # ── THE WRITE-BACK — only SOLVABLE slots ──────────────────────────────────────────────────────────
    # A locked slot (no admissible RNA strand) or an empty one keeps its signature-binary init. ⛔ Do
    # not extend the skip to UNIDENTIFIED slots and defer them to the prior: that arm is refuted, because
    # the prior resolves an imperfectly-solved slot better than a deferred ``f_g = 1``.
    mg_, mp_, mn_ = dc_fin.gdna_frac, dc_fin.rna_pos_frac, dc_fin.rna_neg_frac
    vg_, vp_, vn_ = dc_fin.gdna_frac_var, dc_fin.rna_pos_frac_var, dc_fin.rna_neg_frac_var
    out_fg = np.where(solvable, np.clip(mg_, 0.0, 1.0), f_g)
    out_fpos = np.where(solvable, np.clip(mp_, 0.0, 1.0), f_pos)
    out_fneg = np.where(solvable, np.clip(mn_, 0.0, 1.0), f_neg)
    out_vg = np.where(solvable, vg_, var_g)
    out_vpos = np.where(solvable, vp_, var_pos)
    out_vneg = np.where(solvable, vn_, var_neg)
    # ── the write-back touched ONLY solvable slots ───────────────────────────────────────────────────
    # ⛔ Without this the mask is invisible to a replay, which then compares the solve's raw output
    # against the shipped belief and reads the mask as a difference (TRAPS: byte-identity-gate).
    # Reproducing a pipeline stage means reproducing its WRITE-BACK.
    o = slice(0, n_owned)
    untouched = ~np.asarray(solvable, bool)[o]
    counts.note(
        "writeback_only_solvable",
        untouched
        & (
            (out_fg[o] != _fg_init[o])
            | (out_fpos[o] != _fp_init[o])
            | (out_fneg[o] != _fn_init[o])
            | (out_vg[o] != np.asarray(belief.var_gdna, np.float64)[o])
        ),
        untouched,
    )
    f_g, f_pos, f_neg = out_fg, out_fpos, out_fneg
    var_g, var_pos, var_neg = out_vg, out_vpos, out_vneg

    # ── ``has_composition`` — does this slot hold a COMPOSITION, or only a bound? ─────────────────
    # An own composition channel (the solver's own precision), structural certainty, or a
    # COMPOSITION row received from a neighbour (`Received.has_composition`). A level lane, a ceiling
    # (the RNA lanes' or the node's own flux's) and a cube row are BOUNDS: one-sided, so the value the
    # solve settles on within the admitted half-line is the prior's, and a slot with nothing at all
    # believes the prior outright. A bound-only node does NOT train the landscape
    # (`calibrate._fit_gdna_hyperprior`). ⛔ Read this off the HELD MESSAGES, not `msg.lam_rows`, which
    # fuses compositions and bounds into one row: the bound-only slots with a non-flat row are exactly
    # the ones the training rule excludes.
    has_composition = (
        has_own_composition_evidence(own.tau_lam) | g1_locked(fp, fn) | held_composition
    )

    if _capture is not None:  # inert diagnostic hook
        # strand-ONLY local belief (no global prior, no messages) — to split the local error into the
        # strand likelihood vs the global gDNA prior contribution. Same solver, global=None.
        fg_strand = _solve_regions_logodds_all(
            u_pos,
            u_neg,
            fp,
            fn,
            n_slot,
            spliced_slot,
            kappa=kappa,
            od_g=od_g,
            od_r=od_r,
            n_grid=int(n_grid),
            L=float(logodds_window),
            n_tilt=n_tilt,
            n_grid_ss=n_grid_ss,
            priors=None,
        ).gdna_frac
        # the message-free self-solve variances, for the local-error attribution — a debug-only solve, so
        # the production path carries none of it (the self-solve fractions come from ``own``).
        # Deliberately AFTER the write-back, so its reference is the OUTGOING belief — an instrument
        # comparing against the shipped solve depends on the same reference.
        _dc_loc = _psi(global_lp, PsiMessage.silent(), fg_ref=f_g, fpos_ref=f_pos, fneg_ref=f_neg)
        _capture.update(
            n_slot=n_slot.copy(),
            fg_loc=own.f_g,
            fg_strand=fg_strand,
            fp_loc=own.f_pos,
            fn_loc=own.f_neg,
            vg_loc=_dc_loc.gdna_frac_var,
            vp_loc=_dc_loc.rna_pos_frac_var,
            f_g=f_g.copy(),
            f_pos=f_pos.copy(),
            f_neg=f_neg.copy(),
            var_g=var_g.copy(),
            solvable=solvable,
            count=CNT,
            spliced=spliced_slot,
            mature=fields["sj_count"].sum(axis=1),
            free_pos=np.asarray(fp, bool),
            free_neg=np.asarray(fn, bool),
            eff_global=eff_global,
            mass_global=mass_global,
            eff_gdna=fields["eff_gdna"],
            eff_rna=fields["eff_rna"],
            # the full per-slot global prior term on the solve grid, so a diagnostic can replay the solve
            # with message channels ablated (the message help/hurt attribution).
            global_lp=global_lp,
            _tau0_lam=own.tau_lam,
            # the incoming belief (the final solve's ``fg_ref``) + the intron-factory λ arm, so an ablation
            # replay reproduces the shipped f_g exactly BEFORE ablating. ⛔ The final solve's row
            # factor is ``intron_prior`` PLUS the delivered rows (`msg.lam_rows`) — a replay passing the
            # bare ``intron_prior`` is unfaithful whenever the policy delivers anything, so both are
            # published and a faithful replay sums them.
            fg_init=_fg_init,
            fpos_init=_fp_init,
            fneg_init=_fn_init,
            lam_rows=msg.lam_rows,
            cube_rows=msg.cube_rows,
            # the two tables, the instruments' view of what each node heard
            from_left=from_left,
            from_right=from_right,
            held_composition=held_composition,
            solvable_mask=solvable,
        )

    return dict(
        belief=dict(
            f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_pos, var_neg=var_neg, var_gdna=var_g
        ),
        has_composition=has_composition,
        counts=counts,
        diagnostics=_capture,
    )


def _pass(seq, nbr, prepared, n_grid: int, *, backward: bool, terminal=None) -> Received:
    """ONE directional pass — phase 1 of the two-phase solve: for each node in ``seq``, in chain order,
    the node RECEIVES from its neighbour of the other kind, into its row of the pass's table.

    The whole direction dependence is which neighbour array is read. The backbone owns the table
    (:class:`~.messages.Received`) and writes ``has_neighbour``; the policy's kernel writes the lanes.
    ``-1`` is a reference terminal: the node has NO NEIGHBOUR there — the chain's two end nodes hear
    from one side, every other node from two. A node with a neighbour whose row stays empty holds
    SILENCE: delivered, and nothing to say. The two states are the table's, so a kernel cannot leave a
    node "never spoken to" — it either writes the row or it does not. A policy that sends nothing at
    all returns no kernel, and every node then holds silence from this side.

    ``terminal`` (a bool per node, or ``None`` for none) marks the nodes that RECEIVE NOTHING: the hop
    into one is never asked of the kernel and the node holds silence — delivered, and empty. What a
    terminal SENDS is the policy's business as for any node; what it hears is not, and that is the
    boundary condition the locus solve stands on.
    """
    received = Received.empty(len(seq), int(n_grid))
    receive = prepared.propagate(received, backward=backward)
    has_neighbour = received.has_neighbour
    for i in seq:
        s = nbr[i]
        if s < 0:
            continue
        has_neighbour[i] = True
        if terminal is not None and terminal[i]:
            continue
        if receive is not None:
            receive(s, i)
    return received


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# THE OUTPUT CONTRACT — projecting the chain belief back onto the two payload axes.
# Not part of the solve: this is what ``CalibrationResult`` / ``priors`` / ``derive`` consume.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def chain_region_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's REGION slots back onto the REGION axis as a :class:`RegionDeconv` — what
    ``CalibrationResult`` / ``priors`` / ``derive`` consume.

    A region's contained population carries no spliced term, and that is structural: the accumulator
    credits ``region_contained`` only when the fragment used no sj, so a contained fragment is unspliced
    by construction. ⛔ Do not add a ``+ mass_spliced`` term here — the quantity is identically zero on
    the region axis, so adding it would be adding a channel that cannot exist.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    reg = kind == REGION
    count = np.asarray(substrate.region_contained.count, dtype=np.float64).sum(axis=1)
    n = count.shape[0]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    ri = idx[reg]
    f_g[ri] = belief.f_g[reg]
    f_pos[ri] = belief.f_pos[reg]
    f_neg[ri] = belief.f_neg[reg]
    return RegionDeconv(
        gdna_mass=f_g * count,
        rna_mass=(1.0 - f_g) * count,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )


def chain_boundary_deconv(chain: RegionChain, belief: RegionBelief, substrate) -> RegionDeconv:
    """Project the chain belief's BOUNDARY slots onto the CONTIGUOUS-BOUNDARY axis — the crossing flux that
    ``priors`` and ``derive`` consume.

    ONE per-boundary result, not a ``(left, right)`` pair of per-region ones: splitting each boundary's
    flux onto its two flanking regions only for ``priors`` to pool the halves back together is a no-op,
    and that sum-then-halve pattern is what hides a factor of 2. ``CalibrationResult`` carries
    per-boundary arrays for this reason.

    The RNA mass is spliced-inclusive: a boundary's certified-RNA crossings (``boundary_spliced``) are RNA
    whatever the unspliced mixture resolves to, since gDNA cannot be spliced.
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    boundary = kind == BOUNDARY
    unspliced = np.asarray(substrate.boundary_unspliced.count, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(substrate.boundary_spliced.count, dtype=np.float64).sum(axis=1)
    n = unspliced.shape[0]
    ei = idx[boundary]
    f_g = np.zeros(n)
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g[ei] = np.asarray(belief.f_g, dtype=np.float64)[boundary]
    # THE PER-STRAND RNA SPLIT, PROJECTED ON THIS AXIS TOO. ψ solves the simplex
    # ``(f_g, f_pos, f_neg)`` at EVERY slot — AXIOM 0's `T(slot)`, which is a function of the two
    # `free_*` bits and never of the slot's kind — so a BOUNDARY slot has the same three-way composition
    # a REGION slot does. ⛔ Emitting zeros for both RNA strands here would publish a composition on the
    # crossing axis that sums to ``f_g`` alone, which a per-transcript prior reading composition per
    # object then believes.
    f_pos[ei] = np.asarray(belief.f_pos, dtype=np.float64)[boundary]
    f_neg[ei] = np.asarray(belief.f_neg, dtype=np.float64)[boundary]
    return RegionDeconv(
        gdna_mass=f_g * unspliced,
        rna_mass=(1.0 - f_g) * unspliced + spliced,
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
    )
