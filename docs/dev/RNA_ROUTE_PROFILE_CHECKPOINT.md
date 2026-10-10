# Local route profiling: implementation and robustness checkpoint

*2026-10-07. Development record for owner and external review. No production source,
native binary, existing test or golden was changed. No commit or push.*

**Owner clarification after the first experiment:** capture creates a local footprint,
including partial-probe binding by fragments on either side of a junction. The comparison
below rejects exact equality at one isolated face; it does not reject useful approximate
co-enrichment across a region/boundary neighborhood. The earlier blanket stop was too broad
for that larger local model. A complete-footprint overlap calculation now confirms related
but unequal enrichment of both flanking regions, ordinary boundaries and the junction.
The ruling and measurements have moved to DESIGN and ISSUES. No numerical coupling or new
production rule has been implemented.

## Decision reached

The one-route profiling calculation is small and correct under its stated model. The
model is **not ready to supply general density evidence**: its shared-capture assumption
fails a noiseless change in relative junction capture. The experiment stopped before RNA
integration and calibration, as specified before any measurements. No simulated-panel
accuracy was used to choose this decision, the model or a numerical setting.

This narrows the next problem to opportunities and identifiability. It does not reject
evidence curves, removal of posterior-variance admission, the factory-removal count
candidate, or the goal of applying the explicit DNA landscape once. It does not change
the deferred status of unstranded capture-ON.

## What the method does, in plain language

At a splice boundary we observe unspliced reads on the two strands and reads that cross
the junction. We do not know how much of the local RNA used that splice route. The
prototype tries every permitted route share for each proposed DNA/RNA mixture and keeps
the best explanation of these local observations. It uses no distant RNA abundance.

The calculation can avoid a large search: after normalization, the allowed source
probabilities lie on one line. A single concave maximization on that line gives the
answer. There is no routing prior, hard admission cutoff, fitted correction or detector.

The weakness is what normalization cancels. It removes one multiplier shared by all
source observations. A probe can preferentially enrich the spliced fragments, however,
so their number can grow while continuous fragments stay unchanged. The model reads
that change as a different mixture unless the opportunities describe this relative
capture. Profiling the RNA route share cannot distinguish those explanations.

The noiseless counterexample and placement measurements have one authoritative home in
[ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice). The exact equivalence, including
multiple routes, is in [EQUATIONS.md](../EQUATIONS.md), **Capture-opportunity nonidentification
at a splice face**. These concern the local factor, not a measured release regression.

## Exactly what was implemented

All runnable files are in `.cache/rigel_runs/2026-10-07_route_profile/`:

| File | Implementation |
|---|---|
| `PROTOCOL.md` | Predeclared stages and stop before integration if relative capture confounds density |
| `route_profile.py` | One face, one junction: construct normalized endpoint means; preserve impossible support; maximize the conditional multinomial likelihood with endpoints and a scalar derivative root |
| `test_route_profile.py` | Independent maximization in the original RNA route coordinate; analytic strand/unstranded/zero limits; opportunity, unit and source-scaling contrasts |
| `capture_mismatch.py` | Noiseless relative-capture contrast, oracle opportunity control and independent exhaustive fragment-placement calculation; its generality condition intentionally exits with failure |
| `test_capture_mismatch.py` | Check the counterexample algebra, depth scaling and placement enumeration; passing these characterizes the failure, not model suitability |
| `mutate.py` | Deliberately break units, routing, opportunities, support, empty observations, depth scaling and splice geometry; require every arithmetic/characterization gate to fire |

The route-share placeholder failed 22 of 33 initial gates before the implementation.
The final arithmetic and failure-characterization suite has 40 passing cases. Nine
deliberate defects exercise all 40. One characterization assertion initially required
bit equality of two differently rounded ratios; it now allows floating-point roundoff.
The generality condition remains failed and is not converted into a green gate.

The placement assay enumerates fragments on two exons separated by an intron, moving
one probe between the intron and the second exon. Length gaps run in both directions,
and binding strength varies separately. These are symbolic/finite observation checks,
not panel simulations or probe inputs to the quantifier. They establish a possible
failure mode; they do not estimate its prevalence in real libraries.

No RNA integral, multiple-route optimizer, neighboring-face product, fitted landscape,
capture readout, native code or production setting was added. The existing production
suite and preflight are not claimed as release validation for a disconnected reference.

## Specific correction to the implementation plan

**Reuse the existing messages.** The previous description of a new local-sharing model
was unnecessary. `transfer_kernel.h` already sends composition profiles and RNA/DNA
levels along licensed faces, incorporates junction flux, applies discrepancy widths and
keeps directional propagation. `transfer_rows.h::hop_price` already accounts for counting
noise and excess disagreement for the relevant level hops; alternative-splice messages
have their own per-pair discrepancy calculation. These mechanisms are not blanket proof
that every current approximation is correct, but they must be inspected and reused before
introducing another neighborhood model.

The concrete repair is at the observation/prior boundary and the consumer interface:
retain uncertainty about DNA/RNA origin instead of reducing it to a solved DNA count and
then constructing a Poisson likelihood as though that count were observed. Some primitive
messages also need rebuilding because their widths or source rows contain prior information;
the factory-removal and exact-strand prototypes address demonstrated instances. Existing
messages already carry curves, so “replace messages with curves” is not an accurate task.
The isolated raw splice factor omitted the existing discrepancy treatment; its exactness
failure does not by itself justify replacing the existing propagation policy.

Do not add free capture factors to this face and fit them from the same counts. That
would add unidentifiable parameters. Do not broaden the likelihood until a particular
panel score improves. Integrating rather than profiling the same incorrect opportunity
assumption also leaves the defect in place.

The proposed next checkpoint is an **observation and opportunity audit**, before
restarting integration:

1. Map each currently licensed local comparison to the physical fragment observations
   it uses. Record opportunities, repeated observations and the exact shared-selection
   assumption required for cancellation. Distinguish a common coordinate from a common
   set of capturable fragments.
2. Determine whether the existing recorded data identify that comparison's relative
   opportunities. Verify candidate cancellations by enumerating fragment placements,
   changing probe position and both length laws independently. Expected-yield inputs may
   be used as an oracle control, never silently supplied to the estimator.
3. For a comparison whose selection really cancels, derive the smallest observation
   likelihood from those existing inputs and test it independently. If the data cannot
   identify the transported density, preserve that ambiguity instead of manufacturing
   a point estimate or adding a tuning parameter. The retained DNA landscape is explicit
   prior information, not a substitute observation.
4. Bring that concrete evidence contract back for review before changing the licensed
   assembly. Only then resume the RNA integral, frozen-input replacement, admission
   contrast and prior-once comparison, one mechanism at a time.

The audit's deliverable is one table of observation/selection dependencies plus a small
independent reference for any justified transport. This is bounded work, not a new
general message framework. Own-strand evidence alone remains insufficient as a release
capture reference; omitting all useful local information is not a finished solution.

Following the owner's clarification, the candidate should be evaluated on the **whole
local footprint**. Move one physical probe and let all affected observation banks change
consistently. Keep the earlier independently varied bank multipliers as a test of the
exact-cancellation claim, not the sole acceptance criterion for a neighborhood approximation.
Use the existing region/boundary graph, including splice connectivity; genomic distance
alone does not describe RNA adjacency. No probe coordinates enter the estimator.

The intended locality is already served by the message system. Preserve its justified
neighbor information while retaining distinct opportunities and uncertainty. A long
region can be enriched at one end without enriching its far boundary; a region with zero
contained opportunity supplies no contradictory observation. Neighbor corroboration is
never a capture gate. First establish whether the existing per-face discrepancy treatment
is sufficient for the proposed evidence output; only a demonstrated failure would justify
changing it. Do not start by fitting a new numerical correlation or discrepancy variance
from panel accuracy, or disguise a local regularization prior as independent observation
evidence. Preserve the no-echo accounting when adjacent objects reuse fragments.

The owner clarification supplies the direction for further exploration. A concrete
coupling that changes the licensed assembly still needs its specified review before
integration. Existing measured count improvements remain available; large panel runs
are not the next instrument for choosing this local model.
