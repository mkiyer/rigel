# NEXT SESSION — the reference location is DELETED and priced; build the learned-enrichment channel, then the weak-message re-contrast

    ⚠ **A DEV DOC, and it is a HANDOFF.** Read it, do the work, MOVE what settles into the
    permanent docs, and DELETE this file in the same session. (Precedent: eight previous
    `NEXT_SESSION.md` files were each deleted per this instruction.)

## ① WHAT THE 2026-08-24 SESSION SHIPPED (working tree; the owner drives the commit)

1. ⛔⛔⛔ **THE REFERENCE LOCATION IS DELETED FROM `src/`** (owner ruling, in-session):
   `_location_term`, `structural_reference_location`, `measured_reference_location`, both config
   flags, the `location` member of `CompositionPriors`, and the threading through
   `calibrate`/`sweep`. Three gate files deleted with it; one strict xfail died with them — its
   content (the relay transports a structural claim across a population change) stays recorded as
   a relay-rebuild constraint (`DESIGN.md` §0c.2/§6b.1). Ruling + price recorded in **`DESIGN.md`
   §6b.1**; the stale ROADMAP/EQUATIONS mentions repaired; goldens updated (magnitudes: gdna
   scenarios < 0.03 fragments; nascent scenarios up to 4.5 fragments toward MORE RNA — the
   deleted `m = 0.75` was suppressing intron RNA). **Suite: 0 failed / 3,664 / 8 xfail**
   (CLAUDE.md's baseline line carries the full accounting). `preflight` ✔; the three touched
   instruments' `--self-test` all pass; `ruff` clean.
2. **The measured price, 0.8.0's own metric** (`calibration_vs_oracle.py`, all 16, same-session
   before/after — the table is in `DESIGN.md` §6b.1): stranded × capture-ON +4–8 %, smaller rows
   up to +42 % (`g98` ss.99 OFF region), `g00` capture-OFF IMPROVES, `g00` capture-ON +30–45 %.
   Accepted by the owner's frame: the score was partly measuring the term.
3. **The "0.001 pseudocounts" ask, resolved by measurement, not argument**: sub-Jeffreys exponents
   make the WINDOW the prior — at `a = b = 0.001`, one read's answer moves 0.097 → 0.010 → 0.0001
   as L goes 10 → 20 → 40 (L-dependence ~0.10 in f), and one read snaps to near-vertex certainty;
   Jeffreys is L-invariant (0.0017) and reads one read as its true 2:1 odds. The ½,½ measure is
   the weakest stable point and it asserts NOTHING (it is the arcsine grid measure); with the
   location term gone, "any nonzero data overcomes the prior" already holds — one fragment moves
   the posterior to the likelihood's own odds (`THETA_COORDINATE_PROTOTYPE.md` §4b).
4. **The owner's point-2 measurement frame, settled**: `pass0_claimed_ab.py` IS the
   fanout-substrate-only instrument (it scores the claimed populations and nothing else,
   DELIVER/REFUTE split). ⛔ A finer "solved subset" partition via
   `has_own_composition_evidence` ∪ message-touched was tried and COLLAPSED — the predicate is a
   divide-by-zero guard, not a resolving-power test (its own docstring), so it admits every slot
   with counts and the "unsolved" remainder is just empty slots. Do not rebuild that partition;
   the honest-ignorance frame is `solvability_audit.py`'s fixed denominator.

## ② THE STANDING OWNER RULINGS (2026-08-24) — BINDING

1. ⛔⛔ **Capture gating is refuted as a concept** — capture is a SPECTRUM; enrichment must be
   LEARNED, universally. No mechanism may gate on a capture boolean.
2. ⭐⭐⭐ **The enrichment estimate comes from EXONS.** Stranded data solves exons independently
   (free, no messages needed). Unstranded exons solve only through boundaries, and boundary gDNA
   UNDERESTIMATES exonic gDNA — the core limitation of unstranded × capture-ON (= the deferred
   stratum; scope and analysis agree).
3. ⭐⭐⭐ **Messages extremely WEAK, always propagated** — never able to overpower a live strand
   channel. ⛔ Historical "dampening measured worse" evidence is CONFOUNDED (measured with the
   now-deleted location term live) and may not be cited; re-price on the current tree.
4. **The θ = arcsin(√f_g) coordinate is the architectural direction** (prototype measured
   equivalent; sparse behavior measured desirable — `THETA_COORDINATE_PROTOTYPE.md`).

## ③ THE ORDER TOWARD 0.8.0

1. ⭐⭐⭐ **Build the learned-enrichment density channel** — the win-back for the deletion's price
   and the universal replacement for both the deleted location and the refuted gate. The channel
   form ships already (`density_lambda_factor`); what is missing is the RATE for enriched slots.
   Measured anchors: the intergenic rate is EXACT at unprobed boundaries (μ-check median ratio
   1.005 capture-OFF) and 38–50× low under capture; the boundary-frame factor at the intergenic
   rate took the `g00` zero controls to ~0 and beat the old location 1.9× on the unstranded
   contaminated stratum. Ruling 2 names the estimand: exon-derived enrichment where exons solve
   (stranded), transported with the boundary-underestimates caveat where they do not. The
   acceptance instrument is the μ-check: predicted rate over certified gDNA ≈ 1 on BOTH capture
   states, per stratum, never pooled. Zero controls ~0 stay the free falsification.
2. ⭐⭐ **The weak-message re-contrast** (ruling 3) on the repaired tree: sweep message dampening;
   the stranded strata must be UNHARMED by propagation (the falsification), unstranded helped.
   This is `ROADMAP.md` §1's reference-first ordering with the reference now actually fixed.
3. **The θ coordinate as its own campaign** (prototype doc §8 lists the build items; the
   reference-term test scaffolding it needed replaced is already gone — this session deleted it).

## ④ OPEN ITEMS CARRIED

* The fan-out's high-gDNA residual needs an INFORMATION CEILING (`PLAN_first_pass_redesign.md`
  §H.5 ⑤); blunt form measured-DISQUALIFIED.
* The fl-gap side panels still carry the RETIRED uniform nascent model (owner decision);
  `object_composition.PURE_GDNA_STRATA` still includes `R intron` (recorded, deliberate).
* `test_flank_to_exon_variance_matches_monte_carlo` passes at `tau_b = 30` vs production ~0.17
  (§H.4 ④c) — right assertion, unrepresentative fixture.
* `ROADMAP.md` §0's numbers predate the deletion — re-derive before quoting any of them.

## ⑤ SESSION HANDOFF PROMPT

```
Rigel — the reference location is deleted and priced; build the learned-enrichment density
channel (owner rulings in docs/dev/NEXT_SESSION.md ②), then the weak-message re-contrast.
Branch fragment-length-gold-standard.

⛔ FIRST TWO COMMANDS:
     python scripts/design/preflight.py          (~2 s; --full is ~1 hour, don't)
     python -m pytest tests/ -q
   Baseline: the count in CLAUDE.md's standing-baseline line. ANY failure is a regression.

READ, IN THIS ORDER:
  1. docs/dev/NEXT_SESSION.md                    ← this handoff; delete it at the end
  2. docs/DESIGN.md §6b.1                        ← the deletion ruling and its measured price
  3. docs/dev/THETA_COORDINATE_PROTOTYPE.md      ← the coordinate + the boundary-factor evidence
  4. CLAUDE.md's AXIOM 0 and the 0.8.0 SCOPE

⛔ No capture gating (refuted — enrichment is a spectrum, learn it). Messages weak, always on.
The owner drives commits.
```
