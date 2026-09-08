# THE BOTH-STRANDED LOCUS — the design and the implementation plan (2026-09-06)

The last message-propagation case: REGIONS and BOUNDARIES that admit RNA on BOTH strands — overlapping
transcripts on opposite strands (AXIOM 0's two bits both set, "AMBIG"). Sandbox status: this is the
working design; what settles moves to `DESIGN.md` under the MOVE rule.

## 0. The size of the case, measured (2026-09-06, the landed policy, the ladder)

`ambig_census.py` / `ambig_geometry.py`, full pipeline, error in fragments:

| ladder row | AMBIG nodes | AMBIG share of the row's remaining error |
|---|---|---|
| g50 ss.99 ON (stranded, in scope) | 9,912 of 70,176 (14 %) | **38 %** (83,910 of 220,404) |
| g98 ss.99 ON (stranded, in scope) | 9,912 | **50 %** (100,291 of 201,578) |
| g50 ss.99 OFF (stranded, in scope) | 9,912 | 16.5 % |
| g50 ss.50 OFF (unstranded, in scope) | 9,912 | 13.5 % |

Where it sits: the overlapping loci's exon|exon boundaries (terminus and junction) and their walled
exons. What made the nodes both-stranded, on `g50 ss.99 ON`, in nearly equal thirds:

| overlap geometry (+ strand / − strand) | nodes | AMBIG error | exemplar locus on the ladder |
|---|---|---|---|
| **an exon of one gene inside an INTRON of the other** (+exon/−intron 31 %, +intron/−exon 27 %) | 6,089 | 58 % | TTC28-AS1 (+) inside TTC28 (−), chr22:27,978,013–28,008,581; MYO18B-AS1 (−) inside MYO18B (+) |
| **an exon of one gene over an EXON of the other, no terminus in the overlap** ("spanning") | 1,316 | 21 % | PPM1F-AS1 (+), whose 36 kb second exon covers TOP3B (−) exon by exon; PPIL2/YPEL1; SERHL2/RRP7BP |
| convergent 3' ends (both TES in the overlap) | 321 | 5 % | — |
| divergent 5' ends (both TSS in the overlap) | 270 | 3 % | LINC00649 (+) / ATP5PO (−), chr21:33,915,533–33,915,780 |
| an intron of each inside the other (nested introns) | 1,168 | 2 % | — |

A block that mirrors the first two rows covers four fifths of the case.

## 1. What an AMBIG node is, and why it is a message problem

A node's composition is a point on the 2-simplex ``(f_g, f_+, f_−)`` with TWO degrees of freedom — ψ
parametrises them as ``λ`` (gDNA against RNA) and the tilt ``θ`` (RNA+ against RNA−), and
``_compose`` makes closure structural. A node's own evidence is its unspliced total ``n`` and the split
of those fragments by genome strand, ``u_+`` and ``u_−``. With the library's strand specificity ``κ``:

    E[u_+ / n]  =  ½·f_g  +  κ·f_+  +  (1 − κ)·f_−

ONE equation. At a single-strand node the annotation sets the other strand's share to zero, so the
equation pins ``f_g``: that is why a sighted exon's own solve is excellent on stranded data. At an
AMBIG node the equation pins only a LINE in ``(f_g, f_+)``; ψ marginalises the tilt and the strand
term is Schur-cancelled on the λ axis (`simplex_logodds._solve_ambig_logodds`). ⭐ **So at an AMBIG node
the gDNA share comes from the messages and the prior alone** — today the level lane's gDNA lower bound
(2026-09-05), rung 1's factory row at AMBIG introns, and nothing about the second degree of freedom.

⭐⭐⭐ **THE KEY FACT: a message about the RNA split IS a message about the gDNA share.** Because the one
equation links them, supplying either unknown from outside pins the other through the node's own
strand counts. And the supply need not be two-sided:

    the lane's gDNA level        f_g ≥ b_g                     (a lower bound)
    an RNA+ level from the + gene  f_+ ≥ b_+                    (a lower bound)
    an RNA− level from the − gene  f_− ≥ b_−                    (a lower bound)
    the node's own strand split    ½·f_g + κ·f_+ + (1−κ)·f_− = p_+  (an equation, stranded data)

At ``κ ≈ 1`` the equation reads ``f_g = 2(p_+ − f_+)``, so ``f_+ ≥ b_+`` gives ``f_g ≤ 2(p_+ − b_+)``:
an UPPER bound on the gDNA share from a LOWER bound on RNA — the side the gDNA lane cannot give and
every upper side refused so far could not give honestly. Three one-sided claims plus the node's own
equation bracket ``f_g`` from both sides, and every one of the three is honest under capture ("at least
this much"). On unstranded data (``κ = ½``) the equation says nothing, but ``f_+ ≥ b_+`` and ``f_− ≥ b_−``
still give ``f_g ≤ 1 − b_+ − b_−``: the RNA lanes carry the upper side even where no strand channel
exists. This is the design's centre: **the RNA levels are how the message layer gets a second side without
assuming anything about enrichment.**

## 2. The sources — where an RNA level is MEASURED

The measurement-only law holds (a level is made from a node's own measurement, never from what it
holds as a composition). Two measured sources of a strand's RNA density exist:

1. **A single-strand node's own strand profile.** Its live strand ``s`` carries ``1 − f_g``, so its RNA
   density is ``ρ_s = (1 − σ(λ)) · n / a_r`` — the composition profile read in the RNA coordinate through
   the node's own total and its RNA opportunity ``a_r`` (`eff_rna`, one number per slot: the same geometry
   serves both strands). The same `level_of_profile` map with ``1 − σ`` in place of ``σ``. An intron's
   factory profile is such a source for its gene's unspliced RNA (the nascent level), which is exactly
   what an antisense exon nested in that intron needs to hear.
2. **The certified spliced flux** at an exon's junctions: spliced fragments are RNA of a KNOWN strand
   (the transcript's), measured, never solved (the owner's spliced law). The flux's route rate gives a
   lower bound on that strand's RNA density at the exon — the certified-flux stream's arithmetic
   (`rna_anchor`: route sum, the NB marginal), which the relay delivers as a one-sided RNA claim today
   (`PsiMessage.rna_one_sided`, "the destination holds AT LEAST the RNA the bound accounts for"). This
   source works on UNSTRANDED data and is the deferred stratum's only RNA witness. ⭐ Bringing it into the
   lane retires the relay's last feature the transfer policy lacks.

What is NOT a source: an AMBIG node's own strand split (it fixes a line, not a level: it cannot be
emitted without knowing ``f_g``); a held composition (the law).

## 3. The design — THE RNA LEVEL LANES

⭐ **Owner rulings, 2026-09-08.** (1) THE THREE LEVELS TRAVEL TOGETHER: a message carries gDNA, RNA+ and
RNA− as abundance levels in one object, each component present only where a measurement exists and
empty otherwise, and an empty component is forwarded like the rest. (2) THE CERTIFIED FLUX JOINS PHASE 1
as an RNA source: a strand's spliced fragments at an exon's junctions are that strand's RNA level at
the exon (one hop for the spliced claim itself; the lane's laws after). (3) ONE REPRESENTATION
EVERYWHERE: every message is a profile on the solve grid, a level is a profile over its log density,
and the solve evaluates a held level at the density each grid cell implies — no parametric summary of
any message anywhere in the transfer policy. (4) THE BAR: a step is refused only for harm above about
one percent of a row; differences inside that band are noise unless they are a bug; the calibration
target is a few percent of error, not less.

The gDNA lane's design, per strand, with one change in the witness.

* **Representation.** `Message.level_rna_pos` / `level_rna_neg` are `Level` profiles over
  ``u_s = log(ρ_s / ρ_ref,s)`` on the solve grid — absolute, no map, no knowledge of the recipient, a
  profile because the claims are one-sided. ``ρ_ref,s`` is a coordinate choice per strand (the library's
  strand-``s`` unspliced density over its single-strand exons: Σ counts / Σ RNA opportunity), so the
  window covers expressed and silent genes alike; no constant.
* **Faces.** Strand ``s``'s level crosses a face iff the boundary carries NONE of strand ``s``'s four bits
  (`FLAG_TSS_s`, `FLAG_TES_s`, `FLAG_DONOR_s`, `FLAG_ACCEPTOR_s`) and the recipient admits ``s``
  (`free_s`) — "the strand's population is unchanged at this face", derived from the index, regardless of
  what the OTHER strand does there. In the overlap loci nearly every boundary inside the AMBIG stretch is
  the other gene's feature, so the host's RNA level crosses them all and the antisense's own termini and
  junctions stop its own. The gDNA lane keeps its faces (every face without a composition rule).
* **Emission.** An empty node forwards; a full node emits the INTERSECTION of its own RNA_s level's lower
  side and the RNA_s level it holds (bounds intersect, they do not multiply). ⭐ Refined by phase 0's
  first finding (§4a): a source emits its level on EVERY face, composition or not — levels always travel,
  a composition rides alongside where a map exists — and the recipient's solve reads the composition
  from a side that sent both, the level otherwise. This holds for the gDNA lane too and is the first
  measured step of phase 1.
* **The first hop from an intron into its own boundary is TWO-SIDED** (phase 0's second finding, §4a):
  rung 1's law — one shared unspliced population — makes the intron's RNA_s level exact at its boundary
  whatever the junction bit says; the exon beyond is where the population changes. Lower-only everywhere
  else.
* **The hop price — the one change.** The gDNA lane is priced by the totals' discrepancy because totals
  are all it can see. An RNA_s level is priced by ITS OWN witness: both nodes' strand-``s`` counts'
  counting (``trigamma(u_s + ½)`` each) plus the discrepancy of strand-``s`` count densities beyond
  counting, ``max(0, log(r_s)² − (1/u_s,src + 1/u_s,x))`` with ``r_s`` the ratio of ``u_s / a_r``. The totals
  would be the wrong witness: entering an overlap the total jumps because the OTHER strand's RNA joins,
  while this strand's density is unchanged — the owner's rule applied to the population the lane carries.
  Under capture the same caveat as the gDNA lane (a probe edge inside the locus), the same law: LOWER-ONLY,
  priced, nothing assumed.
* **The recipient — ONE rule for every lane, no Gaussian anywhere (owner's consistency requirement,
  2026-09-08).** Every message in this system is a max-normalised log-profile on the solve grid; a
  level's coordinate is a log density. The solve consumes a held level by EVALUATING it at the density
  each grid cell implies: at a cell ``(λ, θ)`` the three shares are ``f_g = σ(λ)``,
  ``f_+ = (1 − σ)(1 + τ)/2``, ``f_− = (1 − σ)(1 − τ)/2``; the implied density of population ``s`` is
  ``ρ_s = f_s · n / a_s``; the held profile is read at ``log(ρ_s / ρ_ref,s)`` and added to ψ. On the λ
  axis this is exactly what `profile_of_level` does for the gDNA lane today; on the AMBIG cube it is the
  same map with the tilt inside. So the delivery is a per-strand row of the cube, the analogue of
  `lam_rows` — not the relay's Gaussian `rna_mode / rna_prec / rna_one_sided` channels, which are the
  relay's and retire with it. A one-sided profile stays one-sided through the map (the map is monotone),
  so "at least this much RNA+" arrives as a wall in the cube without any parametric summary.
* **The tilt needs no lane of its own.** With both strands' RNA bounds inside the cube, ``θ`` is
  constrained through ``f_+`` and ``f_−``; the `theta` channel stays unused. The tracker's open question
  ("rule whether the tilt has any message") is answered: the tilt's information travels as the two RNA
  levels, and a separate tilt profile would count the same witnesses twice.
* **Single-strand recipients** (phase 3, measured apart): an RNA_s level arriving at a single-strand node
  of strand ``s`` says ``1 − f_g ≥ b_s``, an upper bound on the gDNA share from the same transcript's
  adjacent piece. Honest under capture-OFF; a probe edge inside a transcript makes it over-claim exactly
  as the gDNA lane's cliff did, so it is judged on the sparse and junction panels before it is kept.

**The laws, restated for the lanes.** A level is absolute. A level that crosses a face says "at least this
much" and nothing more. Bounds intersect. A level is made from a measurement, never from a held
composition. Each lane is priced per hop by the discrepancy of its own witness. Nothing pooled.

## 4. The substrate — both-stranded blocks on the test chromosome (the initial phase)

Mirrored from the ladder's loci (the owner's rule), grown ONE structure at a time, each measured before
the next. Every locus is TWO genes on opposite strands (ids ``gB2_{cap}{type}_{regime}_H`` for the host and ``_A`` for
the antisense; the render gate reads the type without the role suffix), the host ``U`` on strand ``h`` in the twin shape
(exon 1 kb · intron 7 kb · exon 1 kb · intron 7 kb · exon 1 kb, span 17 kb) and an antisense ``A`` on the
other strand; the pair shares one type name (``asin_ab`` for both genes), so the render gate's rule that
every type sits on both strands is met inside the locus, and the locus is replicated with the host on +
and on − so termini and junctions are exercised in both orientations.

1. **`asin` — an antisense's exons inside the host's intron** (TTC28-AS1 in TTC28; MYO18B-AS1 in MYO18B:
   58 % of the AMBIG error). ``A`` has two 300 bp exons and a 1.5 kb intron, all inside the host's first
   intron: ``A = [s+3000, s+3300) ∪ [s+4800, s+5100)`` for a + host (mirrored within the span for a −
   host). Nodes it makes: the host intron's pieces outside ``A`` (single-strand, the RNA_h source: the
   host's unspliced level), ``A``'s exons (h-intron / a-exon), ``A``'s intron (h-intron / a-intron), and
   boundaries that are ``A``'s termini and junctions only — faces the host's level crosses and ``A``'s own
   level stops at. FIRST. Built 2026-09-07 as `asin` (host on +) and `asinrev` (host on −, the antisense
   mirrored within the span).
2. **`span` — an antisense exon over a host exon** (PPM1F-AS1 over TOP3B; 21 %). ``A`` has a 300 bp first
   exon in the host's first intron at ``[s+5000, s+5300)`` and a second exon ``[s+7500, s+9500)`` that
   covers the host's middle exon ``[s+8000, s+9000)`` with 500 bp overhangs into both host introns. Nodes:
   the host exon (h-exon / a-exon, walled by ``A``), the overhangs (h-intron / a-exon), and the host's
   own acceptor and donor now INSIDE ``A``'s exon — faces ``A``'s level crosses and the host's mature level
   does not. SECOND.
3. **`conv` — a simple overlap at the ends of two opposite-stranded transcripts** (the owner's
   structure, 2026-09-07; convergent 3' ends, 5 % of the AMBIG error on the ladder). Two two-exon
   transcripts: ``TF`` on + with exons ``[s+1000, s+2000) ∪ [s+10000, s+16000)`` and ``TR`` on − with exons
   ``[s+14000, s+16000) ∪ [s+20000, s+21000)``. Their 3' ends overlap in ``[s+14000, s+16000)``
   (h-exon / a-exon); the boundary at ``s+16000`` carries TF's TES+ and TR's acceptor on −, a junction and
   a terminus on ONE boundary from DIFFERENT strands, so it is also the sj+terminus case (D) in its
   both-stranded form; the boundary at ``s+14000`` is TR's TES− inside TF's exon. Span 21 kb; one
   orientation as specified. IN PHASE 0 with the first two (owner, 2026-09-07).
4. **`div`** — divergent 5' overlaps (3 %; LINC00649/ATP5PO). QUEUED.

Regimes and twins, under the replication rule: host : antisense abundance ``ab`` 90 : 10 (the realistic
lncRNA case), ``ba`` 10 : 90 (a near-silent host under an expressed antisense — the case where the host's
strand carries almost no RNA and a false RNA_h claim would show), ``eq`` 50 : 50; pair total 100 (block 2);
unprobed and HOST-ONLY probed (``cap`` prefix: the probes tile the host's exons and nothing of the
antisense — an enrichment cliff between the two strands inside one node, the adversarial case a real
capture design makes). Per structure: 3 regimes × 2 probing × 2 orientations = 12 loci, 24 genes, 444 kb
at the 37 kb pitch for `asin`; `span` and `conv` in one orientation each (6 loci, 12 genes) in phase 0 —
24 loci, 48 genes, 888 kb in all (owner, 2026-09-07: phase 0 is the infrastructure, all three at once). The chromosome grows and the budget with it as the cluster block did (RNA depth per
existing transcript held within a few percent). No nascent on the hosts in the first cut (the walled
block's convention); the host's nascent level as the RNA_h source is the variant to add when `asin` is
measured, because then the host's intron carries RNA the antisense's exons must not claim.

What the block must show before phase 1 starts (the baseline): the AMBIG census on the new chromosome
(reach and error at AMBIG nodes by class), silent against the landed policy on all 30 conditions and the
three panels, halves apart, pass zero beside the pipeline.

### 4a. Phase 0's baseline (2026-09-07, the benign panel; the junction and sparse panels follow the rebuild)

Every number on the test chromosome moved with the new reads; compare policies within this substrate.
`policy_benchmark.py --panel test`: unstranded 14/20 rows below silence (worst 1.32× at `g50 ss.50 OFF`,
the recorded one-sided first-pass weakness), stranded 7/10 (worst 1.01×), the six zero controls within
1.01×. The AMBIG census (`ambig_census.py test`): 120 both-stranded nodes; on `g50 ss.99 ON` they carry
**15.2 %** of the row's error (1,666 of 10,994) — 1,242 at twelve `B exon|exon [sj]` boundaries and 403
at eighteen walled exons, the `span` loci's host exon walled by the antisense and the host's junctions
inside the antisense's exon: the ladder's signature, reproduced in miniature. On `g98 ss.99 ON` 2.3 %,
on `g50 ss.50 OFF` 1.9 %. By type on `g50 ss.99 ON` (silence → landed): `capasin_ab_H` 130 → 70,
`capasinrev_ab_H` 59 → 25, `capconv_eq_H` 42 → 18 (the gDNA lane's lower bounds already help at the
probed hosts), but `capspan_eq_H` 818 → 957 and `capspan_ba_H` 710 → 704: the host-only-probed spanning
locus is where the landed lane already over-claims — the cliff between the probed host exon and the
unprobed antisense pieces around it, the risk §6 names. On the unstranded in-scope row most
both-stranded types read worse than silence (`capconv_ab_H` 77 → 191, `asinrev_ab_H` 36 → 92): the
plateau at RNA-rich nodes, as everywhere. The junction-probed panel: unstranded 8/10 rows at or below
silence (the capture-OFF rows are the benign panel's), stranded 14/20, worst 1.175× at `g50 ss.70 OFF`,
every capture-ON row below silence. The sparse-probed panel: unstranded 8/10, stranded 12/20, worst
1.175× at `g50 ss.70 OFF` and 1.097× at `g50 ss.99 ON` (11,751 → 12,896) — the latter is the
terminus-cluster block's recorded cliff (`capcluster_*` +1,161 of the +1,145; the both-stranded types
move by less than 20 fragments there).

**What the block showed on its first day (`entrance_probe.py`, `span_ab` at pass zero; the dissect of
`capspan_eq_H`), two findings the phase 1 plan must absorb:**

1. **The stretch's entrance starves the lane.** At the first boundary of the AMBIG stretch the host
   intron's information arrives as a COMPOSITION (rung 1's FORWARD), because the lane serves only faces
   without a composition rule; the next face inward is a lane face, a held composition is never re-issued
   as a level (the law, measured), so the boundary sends SILENCE and every node of the stretch beyond it
   holds nothing — the census's 15.2 % unreached on `g50 ss.99 ON` is exactly the AMBIG nodes. The same
   at the far end: an EMPTY boundary that holds only a composition cannot forward. The fix is not a
   conversion but a change of the lane's faces: **a source emits its level on EVERY face** (levels always
   travel; a composition rides alongside where a map exists), the recipient's own solve reads the
   composition from a side that sent both and the level otherwise, and forwarding is by levels. This
   touches the gDNA lane's reach on the ladder and is measured on its own before anything else in
   phase 1.
2. **The AMBIG boundary needs a two-sided RNA level from its own intron.** The `span` host's acceptor
   and donor sit inside the antisense's exon: AMBIG boundaries with 373 and 395 crossings, true gDNA
   share 0.78, solved at 0.002 by silence AND the landed policy (290 and 307 fragments each — two thirds
   of the locus's error). Their strand split is ``p_+ ≈ ½·f_g`` (the host's unspliced RNA+ at its own
   junction is nascent, ≈ 0; the antisense's RNA− reads on −), a line the local solve resolves toward the
   prior. What pins it is ``f_+ ≈ 0``: an UPPER bound on one strand's RNA, which a lower-only lane cannot
   carry. It is honest here by rung 1's own law — an intron and its boundary share ONE unspliced
   population, so the intron's RNA_s level enters its own boundary TWO-SIDED whatever the junction bit
   says (the unspliced crossing at an acceptor IS the intron's unspliced population; what changes is the
   exon beyond) — and lower-only from there on. §3's face rule is refined accordingly: a strand's level
   crosses into a BOUNDARY from its own intron two-sided, crosses a face into a REGION only where the
   strand carries no feature, and is lower-only everywhere but that first hop.

## 5. The implementation plan

| phase | what lands | gates (falsification first, each watched firing) | judged by |
|---|---|---|---|
| 0 | the both-stranded block — `asin`, `span`, `conv` (owner, 2026-09-07): YAML, renders, index, the seven panels rebuilt, the baseline | the render gate; the strand-balance gate | the AMBIG census on the block; the 30-condition table |
| 1 | the RNA level lanes: sources (own strand profiles; the certified flux at an exon's junctions), AMBIG recipients, delivery as per-strand rows evaluated on the cube — profiles everywhere, no Gaussian | per-strand continuity faces from the flag bits (a junction of + stops RNA+ and not RNA−); the RNA coordinate round trip; the hop price by strand counts; an AMBIG recipient's delivered claim is one-sided in ψ; THE BRACKET THEOREM on a hand-built node (three lower bounds + the strand equation ⇒ a two-sided ``f_g``, and removing any one opens a side); no echo; every lane survives the passes | the block's AMBIG nodes at their destinations; the ladder's AMBIG classes (the census by geometry); pass zero and full; the three panels |
| 2 | source 2: the certified flux as an RNA_s level at exons (the anchor's arithmetic inside the lane) — PROPOSED FOR PHASE 1 (2026-09-08): it is the one RNA witness on unstranded data, and the design is not complete for both data kinds without it | the flux level is a lower bound at the route rate's Poisson width; strictly one hop for the flux itself (the spliced law) | the same, plus the unstranded rows — the deferred stratum reported |
| 3 | RNA_s levels at single-strand recipients (the upper side of ``f_g``) | the sparse and junction panels' cliff rows | refused if a probe panel's stranded row moves more than the noise |
| 4 | `span`, then `conv` / `div` | as phase 0 per structure | as phase 1 |
| 5 | the tilt ruling: the `theta` channel retired from the message contract, or kept for an AMBIG boundary's own split as a local term only | — | a byte-identity check where nothing changes |

Order inside phase 1, one mechanism at a time: (1) LEVELS ALWAYS TRAVEL — the gDNA lane emits on every
face, the solve prefers a composition where both arrive (a change to today's lane, measured alone on the
ladder and the panels: it changes reach); (2) the RNA lanes' faces and plumbing (silent by construction:
no sources yet, byte-identical); (3) the sources — own strand profiles, two-sided into the own boundary,
lower-only beyond; (4) the recipient delivery at AMBIG nodes. Each step A/B'd against the landed policy
before the next.

## 6. Risks, and the refusal criteria

* **The cliff again.** Host-only probing puts a probe edge at every host exon inside the overlap; an RNA_a
  level crossing it over-claims. The lower-only law and the strand-count price are the defence; the
  measurement decides, as for the gDNA lane.
* **The noise ratchet.** A near-silent host (`ba`) has RNA_h modes that are noise around zero; their
  lower sides must not become false floors at the antisense's exons. The intersection form limits this to
  the tightest single witness; the `ba` regime is the control that shows it.
* **The two witnesses of one degree of freedom.** The tilt must not be constrained twice (an RNA level
  and a tilt profile from the same source); phase 5's ruling closes that door.
* **Nothing pooled, no constant.** ``ρ_ref,s`` is a coordinate, the grid is the solve grid, every width is
  counting or a measured discrepancy.

The bar is the completion contract's: improves or stays stable with minimal harm on the test
chromosome, the ladder and every panel, halves apart, pass zero beside the pipeline, zero controls
reported on every experiment.

## 7. Decisions for the owner

1. Host-only probing as the `cap` twin (recommended: it is what a capture design does).
2. Two orientations per structure (recommended for `asin`, the first structure; the render gate would
   accept one).
3. No nascent on the hosts in the first cut, the nascent variant added when `asin` is measured
   (recommended).
4. The half-Gaussian one-sided delivery first, profile rows only if it proves lossy (recommended).
5. The certified flux joins the lane as an RNA source in phase 2 (recommended: it retires the relay's
   last unique feature and serves the unstranded rows).
