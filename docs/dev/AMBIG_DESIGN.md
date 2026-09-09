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

### 4b. Phase 1 as built and measured (2026-09-08, prototype `lanes_proto.py`, not in `src/`)

Four arms, each stacked on the last, each gated with its perturbation watched firing, each A/B'd against
the landed policy on the three test panels, pass zero beside the full pipeline, halves apart:

1. **`lat` — levels always travel.** The gDNA lane serves every directed face not into intergenic; the
   solve reads the composition from a side that sent both, the level otherwise. Gate: the `span_ab`
   stretch's nine AMBIG nodes hold a level under the arm (1 of 9 under the landed policy); a node holding
   both from one side reads the composition (the level-first twin fires). ⭐ **Measured alone it does NOT
   meet the bar**, and the mechanism is the one §6 named: a one-sided floor reaching a node with no
   channel of its own for the gDNA share moves it UP the line. Three faces of the same thing: (i) the
   NOISE RATCHET — on `g00 ss.70 OFF` one exon's strand profile at κ = 0.7 reads noise as gDNA
   (slot 512, f_g ≈ 0.08) and, through the FORWARD faces the landed lane could not cross (an empty
   boundary holding only a composition sends SILENCE), becomes a floor at every exon of its gene
   (+712 on a 9,340 row, all but 148 of which is the shadow transcripts' constant); (ii) on UNSTRANDED
   rows the same floor at channel-free single-strand exons, diffuse (366 slots, none above 35 fragments,
   `g50 ss.50 OFF` +3.6 %); (iii) on the LADDER's `g05 ss.99 OFF` (+7.1 %) and `g05 ss.50 OFF` (+7.8 %) the
   floors land on AMBIG walled exons and exon|exon boundaries with near-zero true gDNA (slot 37345:
   10,188 fragments, truth 0, 3,673 → 3,934). Where it wins it wins big: every capture-ON ladder row of
   both halves (unstranded × ON −15 to −17 %, the deferred stratum; stranded × ON −1 to −3 %), and the
   `capspan` types on `g50 ss.99 ON` (957 → 641, 704 → 567). Faces (iii) are the RNA lanes' destinations.
2. **`rna_plumb` — the RNA lanes' faces, no sources.** Per strand from the flag bits: a junction of `s`
   stops `s` into the exon and is two-sided into the intron; a terminus of `s` stops `s` both ways; no
   face admits a node that does not admit `s`. ⭐ **The intron test must be PER STRAND**: the h-intron ∩
   a-exon piece carries the antisense's EXON signature, so the coarse `is_exon_region` calls it an exon
   for both strands and the host's level stops one hop short of the host's own acceptor — exactly the
   node finding 2 is about. `StepContext` carries no per-strand exon bits; the prototype derives
   `(exon_pos, exon_neg)` from the region signature in the harness (168 two-sided faces on the test
   chromosome exist only under the per-strand test). A landing adds the two bits to `StepContext` as
   geometry. Byte-identical to `lat` (gated).
3. **`rna_src` — the sources.** A single-strand node's own claim as its live strand's RNA level
   (`rna_level_of_profile`: the composition profile read at ``1 − sigma``, the total's Poisson tail above;
   the round trip holds where the coordinate resolves — near ``f_r → 1`` the level saturates at the total
   and ten λ cells share one u cell, the same limit the gDNA level has at ``f_g → 1``); the certified flux
   at each of an exon's junctions as that strand's level at the exon (the count's Poisson likelihood at
   the route rate's own opportunity, lower side; a zero count claims nothing); the pair of an exon's two
   junctions charged their own disagreement beyond counting (item 7's rule). Gates: the round trip; the
   price by strand counts (a hop with equal totals and strand counts 30 vs 3 is priced 5.3 against
   0.0007 by the totals); no echo (a sharpened own level leaves what the node holds unchanged and reaches
   its neighbours); one hop (a junction boundary with flux and no own claim has no RNA level of its own).
   Byte-identical to `rna_plumb` (nothing reads the lanes yet).
4. **`rna_amb` — delivery at AMBIG nodes on the cube.** The held RNA levels (both sides intersected,
   plus the exon's OWN flux level — the spliced claim's one hop, boundary → exon, read at the exon) as a
   ``(K, K_t)`` row: at each cell ``f_s = (1 − sigma)(1 ± tau)/2``, the density ``f_s n / a_r``, the
   profile read at ``log(rho_s / rho_ref,s)``. Gates: THE BRACKET THEOREM on a hand-built node (truth
   0.5/0.3/0.2, n = 400: three lower bounds give a 90 % interval [0.458, 0.542] at κ = 0.99 and at κ = 0.5;
   removing the gDNA bound opens the lower side to 0.003, removing either RNA bound opens the upper
   side); a lower-only profile is monotone on the cube (a two-sided one is not); with the channel disarmed
   the arm is byte-identical to `rna_src`. ⭐ The solver has no cube channel: the prototype wraps the
   AMBIG solve and adds the rows to ψ's strand mixture; a landing adds an inert `PsiMessage` field.

**Two prices were wrong on first contact, and each correction is a ruling already made, applied.**
(a) Across the TWO-SIDED face (an intron into its own boundary) the strand-count discrepancy is huge
(3 + fragments at the intron against 229 at the gDNA-rich boundary) and blurred the sharp upper side to a
slope; but rung 1's identity says those two nodes share one unspliced population, so the difference in
their strand counts is gDNA's half and the other strand, never this population's change — the identity
hop charges COUNTING ALONE. With it the span host's acceptor and donor read 0.775 and 0.758 at pass zero
against 0.780 (0.406 / 0.462 under `lat`, 0.002 under silence). (b) The flux floor at an exon's junction
over-reads the exon's true RNA body density by 0–40 % (mean ≈ 15 %, nine readings; stage 0 measured the
route rate unbiased against mature RNA at −3 % with a 5–9 % scatter beyond counting, so this is the
scatter's tail); at κ = 0.7 a 19 % over-read on one junction collapsed the `capspan_ab` host exon to
f_g = 0.000 (truth 0.342, +458 alone). Two junctions of one exon are two witnesses of one density: their
disagreement beyond counting, charged to both floors, removes it (11,981 against `lat`'s 11,947). A
single-junction exon still pays counting only — the transport dispersion beyond counting has no local
witness there (the owner's owed decomposition).

**Standings on the test chromosome, the final stack against the LANDED policy, full pipeline:**

| panel | stranded (minimal harm) | unstranded (must win) | zero controls |
|---|---|---|---|
| benign | wins 13/16, worst 1.011× (`g05 ss.70 ON`); `g50 ss.99 ON` 10,994 → 9,722 (0.884×) | wins 4/8, worst 1.032× (`g50 ss.50 OFF`, diffuse) | `g00 ss.70 OFF` 1.061× (+570 of which 9,192 is the shadow constant), `g00 ss.99 *` 1.00 |
| junction-probed | wins 16/16, worst 1.000× | wins 6/8, worst 1.032× | `g00 ss.70 OFF` 1.061× |
| sparse-probed | wins 15/16, worst 1.003× | wins 4/8, worst 1.079× (`g05 ss.50 ON`, the deferred stratum) | `g00 ss.70 OFF` 1.061× |

At the design's target: `capspan_eq_H` on `g50 ss.99 ON`, the two junction boundaries 290 / 307 → 7 / 2
fragments and the walled host exon 271 → 126 (the landscape prior gives the lower side there; the flux
gives the upper, ~6 points under truth from the flux floor's over-read). What remains is step 1's residue
on the unstranded capture-OFF rows: a diffuse +1–3 points at channel-free single-strand exons, which is
exactly what phase 3's single-strand recipients (``1 − f_g ≥ b_s``, the upper side) are for.

**The ladder (16 rows, full pipeline, against the landed policy; pass zero in brackets):** the full
stack `rna_amb` — unstranded 4/6 non-zero rows won, `g05 ss.50 OFF` 1.124× (1.076×), `g50 ss.50 OFF`
1.014× (1.141×), the capture-ON rows 0.80–0.83×; stranded 5/6, `g05 ss.99 OFF` 1.031× (0.906×), the rest
0.86–0.99×. The CONTROL `rna_amb_gkeep` (the RNA lanes everywhere, the gDNA lane on the LANDED faces) —
unstranded 6/6, worst 1.000× (0.98×); stranded 6/6, worst 0.994× (0.97×); `g98 ss.99 ON` 0.877×
(0.875×), `g05 ss.99 ON` 0.990× (0.818×), `g50 ss.99 ON` 0.992× (0.908×); the g00 rows identical. ⭐ **So
step 1 is REFUSED by the bar and steps 2–4 on the landed gDNA lane are LANDED (2026-09-08,
`DESIGN.md` §6b.13; `ISSUES: levels-always-travel-for-the-gdna-lane`)**; the landed `transfer` reproduces
the control to the fragment on all 60 (condition, frame) pairs of the benign panel.


### 4c. The owner's two corrections, measured (2026-09-08, `p1b_proto.py`, arms on the landed policy)

**The ruling, in the owner's words.** The splice-junction rate is an ESTIMATE of the adjacent exon's RNA
abundance, not a lower bound (a probe across the junction can put it far above an unprobed exon); the
price of the junction → exon message is the disagreement between the two nodes of that pair, no special
case for an exon with two junctions. And: gDNA measurements as absolute abundance levels can be
propagated — a both-stranded node whose gDNA share is determined may pass its gDNA level on.

**What the measurement said, form by form (three panels + the ladder, halves apart, both frames, against
the landed policy):**

* **Two-sided estimate** (any witness): REFUSED. On the sparse-probe panel's zero control `g00 ss.99 ON`
  it adds 354 false fragments (94 → 448) and `g05 ss.99 ON` loses 4–6 %; the upper side over-claims at a
  probe cliff exactly as every upper side before it (`DESIGN.md` §6b.12). The floor stays lower-sided.
* **The pair's TOTAL abundance as the witness**: worse than the strand's abundance on every row that
  moves (pass zero `g05 ss.99 OFF` +7 %). At an overlap the exon's total holds the other transcript's RNA
  while the junction sees one transcript; the rule charges that as uncertainty and softens a correct floor
  25-fold at the equal-abundance overlap exon (u = 0: −146 → −5.9 nats). ⚠ Two errors were made and
  corrected on the way: the boundary's "total" must include the junction flux (a junction boundary's
  contiguous crossing is ~1 fragment; `StepContext`'s own rule), and the strand count a junction is priced
  against must be read on the genome-strand column that strand's RNA READS on (`read_column`: under an
  antisense protocol, κ ≈ 0.01 here, a 90 % − transcript's reads sit on the + column; the wrong column read
  it as 10 % and priced its floor away).
* ⭐ **LANDED: the lower-sided flux level priced by the node pair's disagreement in the STRAND's
  abundance** — the junction's spliced count at its route rate against the exon's count of that strand
  (through `read_column`) per RNA opportunity, `count_price`'s form, one price for one or two junctions
  (`flux_level(..., v)`). Test chromosome, full pipeline: worst 1.005× / 1.002× / 1.018× (sparse
  `g05 ss.99 ON`, +70 fragments), zero controls unmoved; pass zero worst 1.016×. Ladder, full pipeline against
  the junction-pair form it replaces: stranded worst 1.000× (`g50 ss.99 ON` 0.985×, `g98 ss.99 ON` 0.985×),
  unstranded worst 1.011× (`g05 ss.50 OFF`, +584 on 50,923), the g00 rows identical; pass zero worst 1.018×
  (`g05 ss.99 ON`). The two-sided form on the ladder: `g05 ss.99 ON` 1.016× full and 1.035× at pass zero —
  worse wherever it differs. The both-stranded gDNA source on the ladder: `g05 ss.50 OFF` 1.745×,
  `g05 ss.99 OFF` 1.619×, `g05 ss.99 ON` 1.465× through the pipeline — refused beyond doubt.
* **A both-stranded node emitting its gDNA level** (its own strand counts on the cube plus the RNA levels
  it holds and its own flux estimate, marginalised over the tilt): REFUSED AS BUILT. It costs 2.6 % on
  `g05 ss.70 OFF` on all three panels (11,247 → 11,534), up to 13.7 % on the sparse panel's `g05 ss.99 ON`,
  and +4…+38 % on the weak-κ zero controls — at κ = 0.7 and low gDNA the bracket has little leverage,
  the emitted level's mode is noise, and a noisy level travels as a floor (the ratchet of
  `ISSUES: the-lower-bound-noise-ratchet`, from a new source). Gated on a hand-built node it does what it
  says (mode at the bracketed density); on the chain the nodes that emit are mostly not the ones the rule
  was for. Recorded in `ISSUES: ambig-node-as-a-gdna-source`; the walled overlap exon keeps its lower side
  from the landscape prior for now.

### 4d. Phase 2 — the ceiling at single-strand nodes, built and measured (2026-09-08, `p2_proto.py`)

**The check before building.** The closest refusal on record (`ISSUES: the-certified-flux-row-as-a-level`:
a flux level into an exon, 1.5–4× worse on the junction-probed panel's stranded capture-ON rows) is the
mechanism phase 2 extends to every single-strand exon, so it was measured first where it failed. Two
findings. (1) At every probed junction the flux over-reads the exon body 20–40× (`capspan_eq_H`'s host
junctions: ×21–42 on the junction panel), and at a gDNA-rich exon neither the strand count nor the total
can see it — gDNA fills both. Phase 1's landed flux level at both-stranded exons already pays this there
(`g98 ss.99 ON` 0.988 read as 0.640, ~175 fragments; `g98 ss.50 ON` 0.978 → 0.098); the rows still
improve overall, so it is recorded, not blocking (`ISSUES: flux-floor-dispersion`). At single-strand
exons the strand channel defends on stranded rows and capture-OFF has no over-read, so the exposure
is the deferred stratum. (2) ⭐ At a LICENSED face the composition already carries the flux: rung 2's
splice-in map has the route rate as the spliced density joining at the exon, which CAPS the gDNA share
(`face_map_lambda` saturates at the flux ceiling). Reading the flux again as a ceiling there counts it
twice. Hence the rule: **the ceiling is read only from a face that sent no composition** — the gDNA
lane's own solve rule.

**Census (test chromosome, pass zero).** Under the rule the ceiling reaches single-strand exons where
rung 2 is silent: `g00 ss.99 OFF` 530 exons from their own flux (every junction crossing is empty at
g00), `g50 ss.50 OFF` 48 from held levels and 0 from flux, `g50 ss.99 ON` 58 + 24. The naive form reads
400–650 more.

**Gates** (four, each perturbation watched): the ceiling row is non-increasing in λ and round-trips
with `rna_level_of_profile` (the gDNA map does not); on the zero-gDNA row the own flux ceiling is read
at 472 single-strand exons and every delivered row is a ceiling; NO DOUBLE COUNT — at 281 exons whose
face sent a composition with a flux the `ceil` arm leaves the row identical to the policy before it
while the naive arm changes them; a node's own strand profile never enters its own ceiling. Breaking
the rule in the source fires three suite gates (the new hand-built one, the per-face toy gate, and the
existing recursive-reference gate).

**Standings (against the policy before it, full pipeline):**

| panel | stranded (minimal harm) | unstranded (must win) | zero controls |
|---|---|---|---|
| benign | 14/16, worst 1.001× | 6/8, worst 1.000× | `g00 ss.50 OFF` 0.854×, `g00 ss.70 ON` 0.976× |
| junction-probed | 14/16, worst 1.000× | 7/8, worst 1.032× (`g05 ss.50 ON`, deferred) | `g00 ss.50 OFF` 0.854×, `g00 ss.70 ON` 0.764× |
| sparse-probed | 14/16, worst 1.005× | 6/8, worst 1.012× (`g05 ss.50 ON`, deferred) | `g00 ss.50 OFF` 0.854× |

Pass zero: unstranded 8/8, 8/8, 7/8 and stranded 16/16, 14/16, 13/16, the unstranded zero controls
0.43–0.46×. The naive form (`ceil_all`): `g00 ss.70 ON` 36.6× through the pipeline, `g05 ss.99 ON`
+8.9 % — refused. **Ladder,** full pipeline against the policy before it: stranded 6/6 non-zero rows won, worst 1.000× (`g98 ss.99 ON` 0.994×); unstranded 4/6, every in-scope capture-OFF row won (`g05` 0.997×, `g50` 0.995×, `g98` 0.999×), the two losses in the DEFERRED stratum (`g05 ss.50 ON` 1.012×, `g50 ss.50 ON` 1.011×); the g00 rows identical; at pass zero every non-zero row of both halves won (`g05 ss.50 OFF` 0.958×, `g05 ss.50 ON` 0.939×). The naive form on the ladder: stranded `g05 ss.99 ON` 1.026× — refused.

## 5. The implementation plan

| phase | what lands | gates (falsification first, each watched firing) | judged by |
|---|---|---|---|
| 0 | the both-stranded block — `asin`, `span`, `conv` (owner, 2026-09-07): YAML, renders, index, the seven panels rebuilt, the baseline | the render gate; the strand-balance gate | the AMBIG census on the block; the 30-condition table |
| 1 ✅ LANDED 2026-09-08 (steps 2–4; step 1 refused by the bar, §4b; the two corrections §4c) | the RNA level lanes: sources (own strand profiles; the certified flux at an exon's junctions), AMBIG recipients, delivery as per-strand rows evaluated on the cube — profiles everywhere, no Gaussian | per-strand continuity faces from the flag bits (a junction of + stops RNA+ and not RNA−); the RNA coordinate round trip; the hop price by strand counts; an AMBIG recipient's delivered claim is one-sided in ψ; THE BRACKET THEOREM on a hand-built node (three lower bounds + the strand equation ⇒ a two-sided ``f_g``, and removing any one opens a side); no echo; every lane survives the passes | the block's AMBIG nodes at their destinations; the ladder's AMBIG classes (the census by geometry); pass zero and full; the three panels |
| 2 | source 2: the certified flux as an RNA_s level at exons (the anchor's arithmetic inside the lane) — PROPOSED FOR PHASE 1 (2026-09-08): it is the one RNA witness on unstranded data, and the design is not complete for both data kinds without it | the flux level is a lower bound at the route rate's Poisson width; strictly one hop for the flux itself (the spliced law) | the same, plus the unstranded rows — the deferred stratum reported |
| 3 ✅ LANDED 2026-09-08 as phase 2 (§4d; the ceiling read only from a face without a composition) | RNA_s levels at single-strand recipients (the upper side of ``f_g``) | the sparse and junction panels' cliff rows | refused if a probe panel's stranded row moves more than the noise |
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
