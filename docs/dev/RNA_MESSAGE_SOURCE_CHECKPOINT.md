# Message-source certification checkpoint

*2026-10-08. Prototype implementation and review record. Production source and installed
native library preserved; no commit or push.*

The numerical density reader is now precise enough to audit its inputs. Precision alone
does not make those inputs honest likelihoods. This checkpoint addresses one arithmetic
loss, isolates one source-coordinate mismatch, and identifies a failed assumption in
the existing edge transfer. Neither a new capture reader nor a production change is ready.

**Correction after the placement check:** the proposed absolute-density forwarding below
has been withdrawn. The edge's own DNA measurement remains valid, but its use as a lower
bound on the whole adjacent exon is now authorized for isolated withdrawal. Settled counterexample results
have one home: [the-edge-density-floor-under-capture](../ISSUES.md#the-edge-density-floor-under-capture).

## Implemented arithmetic correction

The isolated worktree's `transfer_rows.h::blur_row` retains the existing finite Gaussian
convolution, edge padding, variance and normalization. It now uses local log-sum-exp
where an ordinary probability sum falls below the floating-point type's normal range.
It retains log kernel weights as well. The ordinary path stays available; the old
probability floor no longer flattens distinct tails.

The code and receipts are in `.cache/rigel_runs/2026-10-08_blur_tails/`. The previous
native package in `2026-10-08_inner_budget/site` remains frozen; the new package is
`blur_tails/site`. `reference.py` computes the same discrete convolution independently.
`test_blur.py` includes deep source tails, narrow kernels, coordinate shifts and zero
width. `mutate.py` compiles changed copies of the actual header. `count_screen.py`
records fresh, separately scored cached-condition pairs. `verify.py` checks receipts,
source identities and preservation of the main tree's pre-existing patch.

Settled test and count results have one home: **BLURRED-LIKELIHOOD TAILS** in
[ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice). The identity is
**A convolution preserves relative likelihood in its tails** in [EQUATIONS.md](../EQUATIONS.md).
This was not an end-to-end transcript benchmark or a new reader integration.

## The source-coordinate problem, in plain language

A gene-edge boundary can supply DNA evidence from its own observed reads and its DNA
opportunity. The count solver expresses that evidence as a possible DNA share of the
recipient's observed reads. That is its current count-facing representation.

The new density reader asks a different question: what expected DNA density could have
generated the reads, allowing the expected RNA amount and observed total to fluctuate?
Reading the old projected row as though it described an intensity changes its meaning.
Even with DNA density held fixed, changing the hypothetical RNA amount changes what
this DNA-only neighbour appears to say. Keeping the source in DNA-density units avoids
that substitution **only under the retained transfer assumption**. It does not establish
that the boundary's density applies to a neighbouring object or can travel across later
faces. The diagnostic uses an adjacent boundary's evidence; it adds no distant RNA
expression or probe information to inference.

The exact algebra and retained one-sided assumption are **An edge's count projection
is not its intensity coordinate** in EQUATIONS. **EDGE-PRODUCER COORDINATE AUDIT** in
ISSUES owns the counterexample, real-table census and attribution measurements.

## Implemented diagnostic adapter

All code is in `.cache/rigel_runs/2026-10-08_edge_source/`:

| File | Responsibility |
|---|---|
| `audit.py` | Check the native EDGE constructor against its formula, reconstruct normalized delivery, and count active direct edges on saved contexts. |
| `edge_sources.py` | For a directly delivered EDGE only, retain source count/opportunity, remove its projected composition term and evaluate the same lower-bound Poisson factor at candidate DNA density. |
| `control.py`, `test_edge.py` | Unchanged-interface falsification, independent Poisson reference, old count projection, source units, source multiplicity and unchanged other channels. |
| `mutate.py` | Break the adapter itself, check each new gate is exercised, and restore it. |
| `replay.py` | Compare two preselected recipients on their existing density grids without fitting a population or changing count-facing delivery. |
| `trace.py` | Remove one direct EDGE source with every other rule fixed, identifying copies already sent through the native passes. |

The other composition terms remain; already extracted DNA and RNA objects are retained
unchanged. RNA exclusions still use the original delivery, so relabelling an edge does
not accidentally enable a second use of neighbouring observations. The adapter accepts
the existing three- or four-factor tuples without creating a second nuisance model.
It is deliberately outside production and is not attached to a population or capture
reader. The census finds no active edge source in the saved zero control; this adapter
cannot repair that control's false capture.

## Why the proposed general repair was withdrawn

The source can travel beyond its first recipient. The native trace shows copies under
forwarding and splice-transfer rules. Relabelling only the final direct EDGE message
leaves those copies in the old coordinate. Nor would adding its DNA factor beside the
old composition copy be acceptable: that would use the source twice.

The earlier proposed next step was to preserve the DNA density through the existing
routes. That was too broad. DESIGN's **The paradigm** explicitly refuses absolute DNA
levels as the basis for transport across capture locales: their levels can differ.
Reusing an existing channel does not make that a numerical or coordinate-only change.

There is a more immediate problem. A probe inside an exon end can overlap most fragments
crossing that boundary while reaching only a small part of the fragments contained
throughout the exon. Both objects can be enriched, just as the owner's locality argument
predicts, while the boundary is more enriched than the exon average. Neighbouring capture
does not imply the density ordering used by the current lower bound.

This concerns the **claim about the neighbouring exon**, not the identity of reads at the
edge. In the model the edge remains pure DNA. Its own count/opportunity evidence remains
useful at the edge and in the existing population learning. No distant RNA information
is needed to make or test this distinction.

## A smaller counting correction was checked and rejected as a complete fix

Let `M` be the boundary count, `N` the recipient's total count, `f` its candidate expected
DNA fraction, and `k = f * Eb/Ee`. The old factor uses the fixed mean `k*N`. Boundary
crossing and adjacent-region containment are disjoint observation banks under the
scanner's deposit rule. Under the same one-sided density ordering, their two Poisson
means can instead be conditioned on the observed sum. For `M>0` and `k*N<M`, the relative
log factor is

```
T = (M + N) / (1 + k)
Q = N * log(T/N) + M * log(k*T/M)
```

Use the zero-count limits; the factor is zero when `M=0` or `k*N>=M`, and minus infinity
when `k=0<M`. This is the corresponding one-sided binomial factor and removes the
fixed-total substitution without a new prior or constant. Its independent binomial and
Poisson-optimization checks are in `.cache/rigel_runs/2026-10-08_edge_contract/calculate.py`.

It is not a selected implementation: it still assumes the false density ordering.
It also cannot be multiplied separately for two edge sources without accounting for
reuse of the recipient total; distinct boundaries can themselves count the same fragment.
The first source-conversion audit therefore did not establish a general likelihood
repair. The successful arithmetic checks and failed generality check are both retained.

`verify_geometry.py` independently uses the actual simulator probe projection and fragment
weights, checks mirrored footprints, and tests overlap survival probabilities. Its checks
verify the counterexample; they do not certify either transfer factor. No calibration,
population fit, benchmark or native modification was made for this follow-up.

## Authorized next experiment

I recommend first pricing the unsupported claim by withdrawing just the EDGE-derived
cross-object floor. This is the smallest experiment that can tell us how much useful
information would be lost. It is not yet a recommendation to ship the ablation.

1. Freeze the existing isolated candidate as control. In the worktree, make the published
   EDGE contribution uninformative at its producer. Keep the EDGE face role and existing
   lane licences so removing it does not silently activate a replacement LEVEL path.
   Preserve the boundary's own observations, opportunities and pure-DNA solve. Preserve
   every other source, face map, hop price and pass. Add no upper side or substitute prior.
2. Write the failing specifications before that change: the physical end-probe example
   must not create this false floor; the pure-DNA source must still constrain its own
   density; other messages must remain available. Check both ends, both strand mirrors,
   zero counts, both length-gap directions and the existing unstranded forwarding trace.
   A source removed at publication must leave no displaced copies downstream.
3. Deliberately restore the floor, leave an onward copy, erase the edge's own measurement,
   enable the fallback path and suppress unrelated messages. Every new gate must reject
   a relevant actual defect. This work does not change the production face references.
4. On cached conditions, compare pass zero and the full calibration separately, including
   per-object region/boundary census, the weak-strand OFF stress case, the zero control,
   and the stranded ON and RNA-long deferred cases. Price this mechanism without changing
   population weights, admission, prior multiplicity, capture reader or numerical support.
5. If the information loss is tolerable, progress to the existing per-stratum test/panel
   gates. If it is not, do not compensate with a fitted cutoff or an upper bound. Return
   with the measured loss and a derivation of any proposed weaker relation before adding
   it. Retain the prior-once and reader-integration checkpoints in either case.

The direct-edge adapter remains an attribution instrument, not the candidate for this
experiment. The source is absent in the saved zero control, so withdrawing it is not a
solution to that control's false capture. The frozen count candidate's previous panel
and real-library results have not changed.

## Owner authorization

The owner authorized this isolated experiment on 2026-10-08. The interpretation and
scope have moved to **The edge floor is an approximation** in [DESIGN.md](../DESIGN.md).
Keep the useful local-imputation role in view; neither the counterexample nor this
permission establishes a large accuracy benefit or makes removal the shipping choice.
The independent message-support proposal remains pending. The separately approved
population-weight comparison is complete in the admission checkpoint; neither experiment
selects a capture prior or authorizes reader integration. Both-strand
consumer cost remains unsuitable for production.

## Completed isolated withdrawal

The native prototype omits the EDGE row at its producer while retaining the EDGE face
role. It changes no pass, map or lane. The package is frozen in
`.cache/rigel_runs/2026-10-08_edge_ablation/site`; its Python assembly still selects the
previous factory removal, so it does not bundle the separate direct-assembly cleanup.

`inputs.py` constructs end-probe examples with equal lengths and both gap directions,
both ends and both strand orientations. `freeze.py` separately zeros only the old
stored EDGE rows to define the expected downstream pass result. `test_edge.py` checks
the absent floor, preserved pure-DNA source and unchanged other sources/licences;
`mutate.py` rebuilds actual changed native code in the worktree and restores it.
`count_screen.py` reads cached observations and read-name truth at pass zero and after
refits. `stress.py` replays the weak-strand OFF fixture against its validated origin
partition. No genome process, new probe input or panel simulation is involved.

The settled census, falsification results and interpretation have moved to
[the-edge-density-floor-under-capture](../ISSUES.md#the-edge-density-floor-under-capture).
The experiment exposes a real tradeoff: the approximation can bias a weak-strand
object, but blanket removal loses valuable unstranded imputation. The ablation is
retained for attribution and is not selected for shipping. Its package stays frozen;
the working native source returns to the pre-ablation implementation. Keep the
zero-gDNA evidence problem and existing release priorities ahead of a broader boundary
model redesign. Any weaker relation would need a derivation and its own experiment.
