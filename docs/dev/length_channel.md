# The fragment-length channel — the current state

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-09, rewritten 2026-08-10 after the A/B. The derivations have moved to
    `EQUATIONS.md` §3d–§3e, the rulings to `DESIGN.md` §3.1c, the lessons to `TRAPS.md`, and the
    numbers to `ROADMAP.md` §0/§1.4. What is left is state and open questions.

---

## 0. ⭐⭐⭐ THE VERDICT

**`length_likelihood` is INADMISSIBLE as it stands** and the switch must stay down. On the `g00`
zero-gDNA control it reports **54–57 %** gDNA in a library containing none. `ROADMAP.md` §1.4 has the
table and the mechanism; `TRAPS: a-single-level-panel-cannot-see-a-constant` has the lesson.

⭐ **The implementation is faithful, and that is worth knowing** — the failure is not a bug hunt:

* the model matches the accumulator's deposit rule exactly (`u = 1/w` at a node,
  `1/(w−1)` at a line, versus `inv_length_quantum(length)` and `inv_length_quantum(length−1)`)
* the 2³² fixed-point descale happens in `substrate.py`; the strand-column collapse is handled
* all three wiring parts exist — bank, moments, precision
* ⭐ its null is EXACT: one pmf on both arms gives `max |ptp| = 0.000e+00` on every condition
* ⭐ the low-`N` Gaussian worry is real but small — 93–98 % of fragments sit at slots with `N >= 50`
  ⚠ substrate-specific: at genome scale (~1.5 M nodes) per-slot depth falls ~40× and this does not carry

---

## 1. WHAT IS OPEN, RANKED

1. ⛔⛔ **The near-flat row.** With no length gap the channel still speaks: at `g00` the two fitted pmfs
   are **1.2 bp** apart (there is no gDNA to fit one on, so it falls back to the anchor), yet the row's
   argmax carries bias **+0.66** with median `|Δ| = 1.0000`. Enabling the channel lets that speak at
   **100 %** of slots. ⭐ `EQUATIONS.md` §3d fixes the PRECISION side of this and is derived, measured and
   ready; it is **not sufficient**, because this failure is in the MODE. ⚠ Whatever fixes the mode must
   shrink the ROW as the pmf gap closes, and no derivation for that exists yet.
2. ⚠ **The RNA moments are mis-specified against their own banks** by 2.3 % (nodes) and 7.5 % (lines) —
   `EQUATIONS.md` §3e. ⛔ The obvious repair is parked: it assumes all RNA at a line is mature mRNA
   confined to an exon, which nascent RNA violates and which is circular besides (owner, 2026-08-10).
3. ⚠ **The gDNA pmf under capture** reads long (node moment ratio 0.976 against 1.000 off capture) — the
   known capture-placement defect, separate and unfixed.
4. ⚠ **Unexplained**: a +2.3 % node residual that survives a perfect pmf; and why the derived line
   opportunity helps off capture and hurts under it while improving the moments on both
   (`TRAPS: a-moment-match-is-not-sufficient`).

---

## 2. ⛔ THE SUBSTRATE GAPS THAT MADE THIS HARD

* **Every condition on both flgap panels is `g50`.** One gDNA level cannot validate a composition
  estimator — it nearly landed this feature on an 87 % improvement.
* **Every condition on both panels is `nrna_none`.** So nothing nascent-related is testable there,
  including item 1's row question and the whole nascent channel in `ROADMAP.md` §0.
* ⚠ Cached payloads are `drain: null` while production drains, and on this channel the drain is not
  small — it adds tail mass to the very pool `rna_pmf` is fitted on (217.93 → 226.10 bp).
  `length_channel_census.py` drains by default; `--undrained` is the diagnostic arm.

---

## 3. THE INSTRUMENTS

| | |
|---|---|
| `length_channel_census.py` | the channel isolated from the solver — does it fire, is the covariance conditioned, is its claim TRUE by count stratum, the smooth-shrinkage arm, and the per-population moment audit |
| `short_exon_fl_probe.py` | offline from read names + index; its PREMISE GATE is what caught the truth-table aggregate bug |
| `length_likelihood_ab.py` | the on/off A/B against real truth, both cache layouts. ⛔ Run the `g00` arm |

---

## 4. MY NOTES

Analysis of the zero gDNA case.

The tool does not know that we have zero gDNA. When we build our gDNA and RNA FL distribution, we do get an initial hypothesis/sense. We may see very low intergenic and intronic fragment counts. This can happen with zero gDNA OR with very strong hybrid capture. Strong hybrid capture of exons will deplete intergenic/intronic regions. So gDNA can appear to be very low in intergenic/intronic regions while it is enriched in exons.

We build our gDNA FL distribution using
- intergenic regions
- intergenic-exon edges (these may be partially enriched by hybrid capture)
- intronic regions (accept the nascent RNA contamination)
- intron-exon edges (may be partially enriched by hybrid capture)

When we have truly ZERO gDNA, all of these fragment sources will either be sparse (close to zero) or contaminated with RNA. 

So in the zero gDNA case, our gDNA FL distribution will be the same as the RNA FL distribution.

Also remember that we SHRINK to the global FL distribution intentionally, so that sparse FL distributions don't mislead us.

So zero gDNA ~ global FL distribution == RNA FL distribution.

So I would argue that in the zero gDNA case, the length likelihood signal should completely disappear! There is no discriminating signal.

So the error that is occuring and causing a catastrophic failure of the zero gDNA scenarios must be a BUG. I argue that the length likelihood model is currently buggy. There is a bug somewhere. Theoretically, the length likelihood channel should elegantly handle the zero gDNA case (the FL distributions shoudl be equivalent) and should elegantly handle the zero RNA case as well for the same reason.

This issue should not become a bug hunt. Why doesn't the length likelihood model handle this elegantly as it should?
