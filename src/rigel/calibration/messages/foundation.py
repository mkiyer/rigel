"""THE FOUNDATION SPEC — the message-propagation architecture, embedded as code
(owner design, ratified 2026-08-26).

**The problem, in one sentence.** Solve one NODE (a REGION or a BOUNDARY) given the node's own
belief (its density, strand and structural evidence), TWO incoming messages (one from each
direction along the chain), and the global inputs (fragment-length pmfs, the intergenic
background, the annotation's structure).

**The two timepoints.** Everything the message layer does happens at exactly one of two moments,
and no piece of arithmetic may straddle them:

1. **PROPAGATION TIME** — the forward scan, then the backward scan. At each hop a node sees ONE
   incoming message and its own belief, and decides what to relay forward. The rule is fixed by
   this spec: *a node with no own belief relays the message unchanged; a node with its own belief
   integrates the two and relays the blend.* What a hop may COST the message (the propagation
   variance model) is the open problem `PropagationModel.attenuate` reserves.
2. **SOLVE TIME** — after both scans, each node holds two travelled messages plus its own belief,
   and decides how much CREDIT each deserves. The recipient decides — never the sender. What the
   solve does with each lane (the solve variance model, including the certified-flux treatment)
   is the open problem `SolveModel` reserves.

**The message.** One type, and the lanes are separated BY PROVENANCE, because the reconcile and
the solve are entitled to treat a measurement differently from a belief:

* the UNSPLICED lanes — one claim per population (AXIOM 0: gDNA, RNA+, RNA−; a fourth population
  is inexpressible). These are deconvolved BELIEFS: products of a solve, imputations, the things
  a premise must be charged on.
* the SPLICED lanes — certified RNA per strand (gDNA cannot splice, so a spliced gDNA lane does
  not exist as a field). These are MEASUREMENTS: certified observations whose weight at solve
  time may be a count-likelihood rather than a Gaussian — the certified-flux treatment is a
  `solve_spliced` implementation, not a side channel.

**The rules the base classes ENFORCE (no implementation may break them):**

* A sender publishes its claim unchanged; deciding how much of an arriving claim to believe is
  the recipient's job (the reconciling node at propagation time; the solved node at solve time).
* A single-source transform may only LOWER a precision — `reconcile` refuses an amplifying
  attenuation at runtime. Precision rises only by the additive fusion of independent witnesses
  (a node's own belief is independent of what its neighbours relayed).
* Lanes never mix in transit: spliced stays spliced, unspliced stays unspliced. Only a SOLVE
  model — the recipient, at solve time — may convert a lane into solve evidence.
* A node with nothing of its own passes a message through byte-identically.

**What this spec deliberately does NOT decide** — the two open problems, each an override point:

* `PropagationModel.attenuate(claim, lane)` — the propagation variance model: what one hop costs.
* `SolveModel.solve_unspliced` / `SolveModel.solve_spliced` — the solve variance model: how
  the recipient converts two travelled messages plus its own belief into evidence.

The backbone (`sweep.solve_chain`) already runs the two timepoints — the scans are propagation
time and `deliver` + the final ψ are solve time — and its `Policy` protocol is where a
spec-conforming implementation plugs in: a policy's scan step is a `PropagationModel`, its deliver
is a `SolveModel` returning the `PsiMessage` ψ consumes. The shipped policies predate this spec
and are grandfathered; the unified policy is built AGAINST it.

Gate: ``tests/calibration/test_message_foundation.py`` — every rule above is asserted there and
was watched firing against a deliberately broken build.
"""

from __future__ import annotations

import dataclasses
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from . import PsiMessage

__all__ = ["Claim", "Hop", "Message", "PassThroughPropagation", "PropagationModel", "SolveModel"]


@dataclass(frozen=True, slots=True)
class Hop:
    """One hop's identity: source slot, destination slot, scan direction. Handed to
    `propagate`/`attenuate` so a model can read its own per-hop tables (frames, licences) —
    the model builds those in its `prepare` hook, from the context."""

    src: int
    dst: int
    backward: bool


#: Tolerance for the no-amplification check — floating-point headroom, not a policy knob.
_AMPLIFY_TOL = 1.0 + 1e-12


@dataclass(frozen=True, slots=True)
class Claim:
    """One population's claim on one lane: an abundance (counts per unit opportunity, in the
    carrier's frame), a precision (inverse log-variance; exactly 0 means NO claim), and the
    MEASUREMENT stream ``measured`` — the share of that weight independent witnesses actually
    COUNTED (a struct-locked gDNA anchor's fragments, a certified sj flux), as opposed to what
    a solve merely believes.

    ⭐ MEASUREMENT CITIZENSHIP (the shipped relay's law, ported 2026-08-26): ``precision``
    weights VALUE fusion in transit, but what a solve may DELIVER as a psi channel precision is
    only ever ``measured`` — an imputation may inform the blend, never masquerade as counted
    evidence. ``measured <= precision`` holds by construction (every seed feeds both streams or
    only ``precision``; every transit cost damps both) and is a property, not an enforced
    invariant."""

    abundance: float
    precision: float
    measured: float = 0.0

    @classmethod
    def silent(cls) -> "Claim":
        return cls(abundance=0.0, precision=0.0, measured=0.0)

    @property
    def is_silent(self) -> bool:
        return float(self.precision) <= 0.0


@dataclass(frozen=True, slots=True)
class Message:
    """THE message — the only thing that travels, and the only shape a node's own belief needs.

    Three unspliced lanes (the populations of AXIOM 0) and two spliced lanes (certified RNA per
    strand; a spliced gDNA lane is structurally inexpressible). Provenance is the field name."""

    unspliced_gdna: Claim
    unspliced_rna_pos: Claim
    unspliced_rna_neg: Claim
    spliced_rna_pos: Claim
    spliced_rna_neg: Claim

    #: the lane names, in one place — iteration order for every lane-wise rule
    LANES = (
        "unspliced_gdna",
        "unspliced_rna_pos",
        "unspliced_rna_neg",
        "spliced_rna_pos",
        "spliced_rna_neg",
    )

    @classmethod
    def silent(cls) -> "Message":
        return cls(*[Claim.silent() for _ in cls.LANES])

    @property
    def is_silent(self) -> bool:
        return all(self.lane(name).is_silent for name in self.LANES)

    def lane(self, name: str) -> Claim:
        return getattr(self, name)

    def with_lane(self, name: str, claim: Claim) -> "Message":
        return dataclasses.replace(self, **{name: claim})


def _fuse(own: Claim, incoming: Claim) -> Claim:
    """Additive fusion of two INDEPENDENT witnesses: precision-weighted mean, precisions sum —
    the one licensed way a precision may rise. Both streams add: the belief stream weights the
    value; the measurement stream simply accumulates its independent counted witnesses."""
    p = own.precision + incoming.precision
    a = (own.precision * own.abundance + incoming.precision * incoming.abundance) / p
    return Claim(abundance=a, precision=p, measured=own.measured + incoming.measured)


class PropagationModel(ABC):
    """PROPAGATION TIME — how a node relays one incoming message.

    The skeleton is the spec and is FINAL; a variance model implements only `attenuate` (how much
    one hop weakens the arriving claim, per lane). The base enforces: pass-through when the node has
    nothing, lane isolation, and the no-amplification rule."""

    def propagate(self, own: Message, incoming: Message, hop: "Hop | None" = None) -> Message:
        out = {}
        for lane in Message.LANES:
            arriving = incoming.lane(lane)
            weakened = self.attenuate(arriving, lane, hop) if not arriving.is_silent else arriving
            if float(weakened.precision) > float(arriving.precision) * _AMPLIFY_TOL or float(
                weakened.measured
            ) > float(arriving.measured) * _AMPLIFY_TOL + (
                0.0 if float(arriving.measured) > 0.0 else 1e-30
            ):
                raise ValueError(
                    f"the attenuation amplified lane {lane!r} "
                    f"({arriving.precision}/{arriving.measured} -> "
                    f"{weakened.precision}/{weakened.measured}): a single-source "
                    "transform may only LOWER a precision (the foundation spec's rule, "
                    "both streams)"
                )
            mine = own.lane(lane)
            out[lane] = weakened if mine.is_silent else _fuse(mine, weakened)
        return Message(**out)

    def prepare(self, ctx) -> None:
        """Optional hook: build per-hop tables (frames, licences, the fitted premise) from the
        context, once per sweep. The default needs nothing."""

    @abstractmethod
    def attenuate(self, claim: Claim, lane: str, hop: "Hop | None") -> Claim:
        """OPEN PROBLEM (propagation variance): how much one hop weakens an arriving claim —
        a signal loses strength in transit, and this says how much."""


class PassThroughPropagation(PropagationModel):
    """The reference reconcile: a free hop. The spec's own control implementation — it exists so
    the rules are testable, not as a production model (a free hop lets an imputation arrive at
    full strength beside a real measurement, the measured failure attenuation exists to price)."""

    def attenuate(self, claim: Claim, lane: str, hop: "Hop | None" = None) -> Claim:
        return claim


class SolveModel(ABC):
    """SOLVE TIME — how the solved node converts its two travelled messages and its own belief
    into the evidence ψ consumes.

    The recipient decides; the template only fixes the ASSEMBLY: the unspliced solve supplies
    ψ's Gaussian channels, the spliced solve supplies the row-factor likelihood (the
    certified-flux treatment lives there), and the base joins them into one `PsiMessage`."""

    def prepare(self, ctx) -> None:
        """Optional hook: build per-sweep tables (the node totals, licences, the grid domain)
        from the context, once per sweep. The default needs nothing."""

    def solve(self, own: Message, forward: Message, backward: Message) -> PsiMessage:
        gaussian = self.solve_unspliced(own, forward, backward)
        rows = self.solve_spliced(own, forward, backward)
        if rows is None:
            return gaussian
        return dataclasses.replace(gaussian, lam_rows=np.asarray(rows, np.float64))

    @abstractmethod
    def solve_unspliced(self, own: Message, forward: Message, backward: Message) -> PsiMessage:
        """OPEN PROBLEM (solve variance, belief lanes): the Gaussian channels ψ receives."""

    @abstractmethod
    def solve_spliced(self, own: Message, forward: Message, backward: Message):
        """OPEN PROBLEM (solve variance, measurement lanes): the row-factor likelihood, or
        ``None`` for no claim — the certified-flux arithmetic is an implementation of THIS."""
