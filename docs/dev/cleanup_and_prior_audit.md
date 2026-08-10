# The purge, the simplification, and the density→count audit

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Written 2026-08-10 against `f470a570`. Two things: (1) the audit the owner asked for — is the EM
    prior formed correctly from calibration's fractional output — and (2) the cleanup plan.

---

## 0. ⭐⭐⭐ THE AUDIT — AND THE ANSWER IS YES, WITH ONE EXCEPTION WORTH THE WHOLE DOC

**The question.** Calibration operates on DENSITY (counts/bp) and outputs fractional composition
(`f_g`, `f_pos`, `f_neg`). Those fractions apply to densities, not counts — so converting them to a
pseudo-COUNT prior by multiplying a count by `f_g` would be wrong whenever the two components have
different fragment-length distributions, because then `E_g ≠ E_r` and a density share is not a count
share.

**The answer: `f_g` is a COUNT share everywhere in this tree, never a density share.** The worry is
exactly right in general and does not bite here, because the tree computes the count share *first* and
derives the density *from* it — which is the safe direction. Traced end to end:

| # | where | quantity | units |
|---|---|---|---|
| 1 | `accumulator.h` | `count` (+1 per touched object), `mass` (conserved, 1 per fragment) | fragments |
| 2 | `node_init.py:335-343` | `rho_g = f_g·M/E_g`, `rho_rna = f_rna·M/E_r` | density |
| 3 | `density_model.py:183` | `count_gdna_frac = rho_g·E_g / N` | **count share** |
| 4 | `sweep.py:644,674` | `gdna_mass = f_g·count` | count |
| 5 | `priors.py:329-341` | `a_g = Σ_nodes f_g·count + Σ_edges f_g·unspliced·q` | fragments |
| 6 | `em_solver.cpp:699` | `G = n_gdna + a_g` | fragments |

**Step 2 is the load-bearing one.** `node_gdna_geometry` (renamed from `node_global_geometry` in this
pass — see §1.3) returns `geometry.eff_gdna`, so `rho_g = f_g·M/E_g` and `rho_rna = f_rna·M/E_r`, giving

    Σ_c rho_c · E_c  =  (f_g + f_pos + f_neg) · M  =  M

⭐ **That identity IS the statement that the fractions apportion the COUNT.** The fragment-length
difference is carried entirely by `E_g` vs `E_r` and never by the fraction — which is precisely why
multiplying a count by `f_g` is legitimate.

**Step 3 is the only place a density becomes a fraction, and it does the conversion properly.**
`count_gdna_frac = density · region_eff_len / contained_gdna`, with `region_eff_len = geometry.eff_gdna`
(`density_model.py:119`) — i.e. `rho_g·E_g/N`, the gDNA's own effective length. ⭐ Its own docstring says
so: "the gDNA fraction of the contained unspliced population implied by the local gDNA density." And it
is not even a vote in the solve — it seeds the gDNA strand fit (`gdna_strand.py:379`) and nothing else.

**So the owner's worked example resolves as follows.** "100 counts, 80 % gDNA — is it 80 gDNA counts?"
If the 0.80 were `rho_g/(rho_g+rho_r)`, then **no**: the count share would be
`rho_g·E_g / (rho_g·E_g + rho_r·E_r)`, which differs from 0.80 whenever `E_g ≠ E_r`. But the solver's
0.80 is already the count share, so **yes, 80 counts is right** — and the density, if you want it, is
`0.80·100/E_g`.

### 0.1 ⛔⛔ THE ONE EXCEPTION, AND IT IS A REAL RESIDUAL FL DEPENDENCE IN THE PRIOR

`assemble_priors` converts the EDGE term from a crossing count to a fragment count with a single scalar
per line (`priors.py:329-338`):

    q  =  calibration.edge_mass_per_crossing        # the accumulator's own mass / count at that line
    gdna_edge = mass_gdna_edge · q
    rna_edge  = (mass_rna_edge − mass_rna_spliced_edge) · q

⭐ The logic is right — a fragment crossing `K` lines is `+1` on each, so the count must be divided by
the mean number of lines crossed, and `mass/count` is exactly that, measured rather than modelled.

⛔ **But `q` is measured on the POOLED population and applied to the gDNA and RNA parts separately.**
`q = [min(w−1,a) + min(w−1,b)] / 2(w−1)` under a uniform field — an explicit function of `w`. A longer
fragment crosses more lines, so it has a *lower* mass per crossing. **So gDNA and RNA have different true
`q` whenever their length distributions differ, and the assembler gives them the same one.**

This is the owner's worry, one level deeper than where it was being looked for: the *fractions* are
dimensionally safe, and the *count-per-crossing conversion* is not.

**Direction:** if gDNA is shorter than RNA, gDNA's true `q` is higher than the pooled `q`, so the gDNA
edge term is **under-counted** and the RNA edge term over-counted — and vice versa. It vanishes exactly
when the two length distributions coincide, which is why the ladder (equal configured lengths) cannot
see it and why nothing has ever measured it.

### ⭐⭐⭐ MEASURED 2026-08-10 — `edge_q_population.py`, DRAINED, with a perfect `f_g`

**The defect is REAL, it is NOT driven by fragment length, and it is SMALL on the deliverable.**

| condition | `mu_g−mu_r` | `q_g` | `q_r` | edge share | EDGE err | **TOTAL err** | **`Δphi` TOTAL** |
|---|---|---|---|---|---|---|---|
| flgap_short capture_off | −83.94 | 0.7361 | 0.5243 | 3.3 % | +0.629 % | +0.021 % | **+0.00013** |
| flgap_short capture_on | −96.23 | 0.6771 | 0.5890 | 26.0 % | −1.017 % | −0.264 % | **−0.00158** |
| flgap_long capture_off | +102.21 | 0.5690 | 0.5234 | 6.2 % | +4.794 % | +0.295 % | **+0.00183** |
| flgap_long capture_on | +105.96 | 0.5262 | 0.5889 | 51.2 % | +1.950 % | +0.998 % | **+0.00596** |
| ⭐ ladder capture_off (**the null**) | +4.98 | 0.6330 | 0.5233 | 4.8 % | +3.182 % | +0.153 % | **+0.00095** |
| ladder capture_on | +15.33 | 0.5783 | 0.5890 | 40.3 % | +1.011 % | +0.407 % | **+0.00244** |

⛔ **THE HYPOTHESIS THAT MOTIVATED THIS IS REFUTED.** `q_g ≠ q_r` on *every* condition including the
equal-length null, where a 4.98 bp gap still gives `q_g` 0.6330 against `q_r` 0.5233 — a 21 % relative
difference, and an EDGE error of **+3.18 %**, *five times larger* than flgap_short's +0.63 % at a
17× larger gap. The sign does not track the gap either: flgap_short (−84 bp) and flgap_long (+102 bp)
both give `q_g > q_r` off capture, and flgap_long *reverses* under capture.

⭐⭐ **THE REAL MECHANISM IS PLACEMENT, NOT LENGTH — and that is a better explanation.** `q` falls toward
the node spacing where the flanking nodes are shorter than the fragment. gDNA is genomically uniform, so
it crosses lines in intergenic regions where nodes are long and `q → 1`; RNA is confined to transcripts,
where exon nodes are short, so an RNA fragment crosses more lines and carries less mass per crossing.
**The two populations occupy different parts of the partition, so they have different `q` even at
identical fragment lengths.**

⚠ **This kills the repair as well as the hypothesis.** A production `q_c` cannot be modelled from the
fitted pmfs, because the driver is *where* each population sits, not *how long* it is — and the pooled
bank is the only per-line evidence production has.

✅ **AND THE MAGNITUDE SAYS DO NOTHING.** The EM consumes the composition of the TOTAL prior, and the
node term (which needs no conversion at all) dilutes the edge error heavily: `Δphi` on the total prior is
**+0.00013 … +0.00596**, i.e. at most **0.6 percentage points**, and ≤ 0.24 pp on every capture-off arm.
Against calibration's own median library `f_gdna` error of 0.005–0.012 on the healthy strata, this is at
or below the noise floor of the thing it would be correcting.

⭐ **The one structural statement worth keeping**: the error scales with the EDGE SHARE of the prior
(3.3 % → 51.2 % across these conditions), and hybrid capture inflates that share, so the worst case is
capture-ON at a large length gap. It is bounded, it is measured, and it is not worth a mechanism.

**VERDICT: record the bound, build nothing.** `TRAPS.md` gets one line — *a conversion measured on a
pooled population and applied per component is population-blind, and here the blindness is placement,
not length* — and `DESIGN.md` §3.1c gets the ≤ 0.6 pp bound beside its existing `Δphi ≤ 5e-4` claim, so
the two are not confused (they are different quantities: that one is the shipped `f_g` against the
unspliced pool, this one is a perfect `f_g` isolating the `q` conversion alone).

### 0.2 Two smaller things the same audit turned up

* **`rho_ref` is one global reference** (`priors.py:353`, `_global_reference_density`) and the gDNA
  effective-length contraction is measured against it per object. Under capture the four gDNA pools
  disagree by **44 bp** on their own mean (`length_salvage_ab.py`, ladder g98 capture_on: 251.67 /
  251.80 / 279.92 / 295.71) — a truth-free measure of non-uniform placement. A single `rho_ref`
  inherits that. Not a dimensional error; a stated assumption that is measurably false under capture.
* **The prior's strength is 1.000 pseudo-fragment per real unspliced fragment by construction**
  (`DESIGN.md` §3.1c) and there is no knob. That is settled and this audit found nothing against it.

---

## 1. ✅ THE PURGE — DONE for everything that does not touch the accumulator

**Landed 2026-08-10, byte-identical on every production path.** The flag defaulted `False`, so removal had
to change nothing — verified by hashing `mass_gdna_node/edge`, `mass_rna_node/edge`,
`mass_rna_spliced_edge`, `edge_mass_per_crossing` and both `gdna_*_eff_len` on two cached ladder
conditions before and after: identical through the wiring removal, through the moments move, AND through
the rename.

Removed: `length_likelihood.py`, its test, the `CalibrationConfig` flag, `_build_length_loglik`, the
parameter from `sweep` and from **five** `simplex_logodds` signatures, the `tau_len` block in `node_init`,
the two dead `NodeGeometry` length channels, the `_layers` entry, the `_synthetic` fixture helper, and the
seven instruments that existed only to price the feature. Suite 3,288 → 3,219, and the −69 reconciles
exactly: 43 (the channel's test file) + 21 (7 scripts × 3 file-parametrised suites) + 2 (the src module)
+ 2 (its test file) + 1 (`test_no_import_points_UP_a_layer`).

⭐ `node_init`'s docstring now says "four information sources" and means it.

### 1.1 ⭐ The moments were MIS-FILED, not dead — and moving them is the simplification

`contained_moments` / `crossing_moments` / `build_slot_moments` compute the OPPORTUNITY-TILTED length
moments, and `pass0_vs_oracle.py` reads them for its `pi_hat` identification classifier. They are a
geometry — functionals of the same opportunities `contained_eff_length` and `crossing_eff_length` return —
and they sat in a layer-5 composition module only because that was their first caller. They now live in
`effective_length.py` (layer 2), beside the opportunities they are functionals of, and the layering test
passes. ⚠ **They arrived without their tests**, which went with the deleted file: restoring moment-only
coverage is OWED.

### 1.2 ⛔ WHAT REMAINS, AND IT HAS A COST NOBODY HAS PRICED

The accumulator BANKS are still stored and now read by nobody: `node_contained_length_sum`,
`edge_unspliced_length_sum`, `node_contained_inv_length_sum`. Removing them touches `bam_scanner.cpp`,
`accumulator.h/.cpp`, `scan_payload.py` and `tests/native/_accumulator_reference.py` (which leads), and
**invalidates every scan cache by schema digest, by design**.

⛔ **That is the decision to take before doing it**: ~40 ladder + 8 flgap conditions × 4 payloads each
would need re-scanning — hours of compute, on the substrate every measurement in this campaign runs on.
The saving is real (three uint64 banks; `Node` and `ContiguousEdge` both shrink) but it is a MEMORY win,
not a correctness one.

⭐ **`edge_unspliced_inv_length_sum` and `sj_inv_length_sum` STAY** — both are live in `second_pass`. And
`PopulationView.inv_length_sum` stays: `pass0_vs_oracle.py` reads it. Only `mean_length` and
`total_inv_length_sum` are genuinely dead there.

### 1.3 ✅ The rename that mattered

`node_global_geometry` → **`node_gdna_geometry`**, `eff_global` → **`eff_gdna_slot`**, 25 sites in 9
files. It returns `eff_gdna` and nothing else, which is what makes `Σ_c rho_c·E_c = M` hold and therefore
what makes `f_g` a COUNT share — the exact question §0 was commissioned to answer, and a reader had to
open two files to learn that "global" meant "gDNA". The docstring now says so.

⚠ **Still owed** (§2 items 2-3): `mass_gdna_edge` is `f_g × unspliced COUNT` while `mass_gdna_node` is a
per-fragment count — same prefix, different units — and the `second_pass.py:443-445` comment is measured
false in premise (harmless in effect). Both have a wider blast radius than this pass took on.

## 2. ⭐⭐ THE SIMPLIFICATION — the naming defects this audit found, ranked by how much harm they can do

⛔ These are exactly the shape of `TRAPS: two-masks-one-name`, and the first two are in the prior chain.

1. **`node_global_geometry` returns `eff_gdna`, and `node_init` binds it as `eff_global`.**
   (`node_geometry.py:329-341`, `node_init.py:315`.) A reader auditing "does the prior mix a density
   share with a count share?" must open two files to learn that "global" means "gDNA". ⭐ **Rename to
   `node_gdna_geometry` / `eff_gdna`.** This is the single highest-value rename in the package: it is
   the one place where a wrong reading produces the exact bug this audit was commissioned to look for.
2. **`mass_gdna_edge` is `f_g × unspliced COUNT`, not a conserved mass** — it becomes a mass only after
   `priors.py` multiplies by `q`. `mass_gdna_node` *is* a count-per-fragment, so the two fields share a
   prefix and differ in units. ⭐ Either rename the edge fields (`count_gdna_edge`) or convert at the
   source. `DESIGN.md` §0 already bans synonyms; this is the same rule applied to units.
3. **`second_pass.py:443-445`** asserts `sj_inv_length_sum` and `edge_unspliced_inv_length_sum` are "the
   same quantity on the same scale". True of the deposit (`1/(L−1)` both), false of the opportunity
   (`w−1` vs the tapered exonic crossing count, so `E/rho` = 0.0319 at a 50 bp reach against exactly
   1.000). ⭐ Measured this session: it does **not** flip choices (0.3-0.5 % of truly-spliced records
   lost, flat-to-inverted in reach), so this is a comment correction, not a code change — but the
   comment as written would license a future reader to build on a false premise.
4. **`inv_length_sum` is deliberately not called `density`** (`scan_payload.py:31`) and that discipline
   is why the model-free/not-model-free distinction survived. Items 1-3 are the same discipline applied
   where it was skipped.

### 2.1 Structural cleanups that fall out of the purge

* `substrate.PopulationView` loses `length_sum`, `mean_length`, `total_inv_length_sum` → the class
  becomes count + inv_length + mass, which is what its docstring already describes.
* `node_init.build_node_init` loses a parameter and a source, so "three live sources" is literal.
* `simplex_logodds` loses a parameter from **five** signatures — it is threaded through the 1-DOF path,
  the 2-DOF path, the block solver and two regrids. That is the largest single reduction in argument
  count available in the package.
* Re-run `python scripts/design/module_census.py` afterwards: it re-derives the layering, flags upward
  imports, and lists docstrings naming a sibling with no import edge (6 were genuinely stale at the last
  measurement). ⭐ Do this **after** the purge, not before — the purge will create new stale references.

---

## 2b. ⭐⭐⭐ THE BANK LEDGER — what is live, what is dead, and what the mass is for

**Two different things are called "mass" and only one of them is a bank.**

* **`CalibrationResult.mass_*`** — `mass_gdna_node/edge`, `mass_rna_node/edge`, `mass_rna_spliced_edge`,
  `edge_mass_per_crossing`. COMPUTED by calibration (`f_g × count`), consumed by `assemble_priors`.
  ⭐ Untouched by the purge and verified byte-identical. These are the fraction→count conversion.
* **`edge_unspliced_mass`** — the accumulator BANK, the conserved fragment mass, fixed point. ⭐⭐ **This
  is the one that makes the conversion possible**: `q = mass/count` is what turns an object-INCIDENCE
  total into a FRAGMENT count. ⛔ It must never be removed.

| bank | bytes | consumer | verdict |
|---|---|---|---|
| `node_contained_count[2]` | 8 | geometry, strand model | ✅ LIVE |
| `node_contained_inv_length_sum` | 8 | — | ⛔ **DEAD** |
| `node_contained_length_sum` | 8 | — | ⛔ **DEAD** |
| `edge_unspliced_count[2]` | 8 | geometry, strand model | ✅ LIVE |
| `edge_spliced_count[2]` | 8 | `calibrate.py` → `mass_rna_spliced_edge` | ✅ LIVE |
| `edge_unspliced_inv_length_sum` | 8 | `second_pass.py` genomic hypothesis | ✅ LIVE |
| `edge_unspliced_length_sum` | 8 | — | ⛔ **DEAD** |
| ⭐ `edge_unspliced_mass` | 8 | `mass_per_crossing` → `q` → `assemble_priors` | ✅ **LOAD-BEARING** |
| ⚠ `edge_spliced_mass` | 8 | — | ⛔ **DEAD — but confirm intent, see below** |
| `sj_count[2]` | 8 | geometry + aligner-artifact detection | ✅ LIVE |
| `sj_inv_length_sum` | 8 | `second_pass.py` spliced hypothesis | ✅ LIVE |

**Struct impact if the four dead banks go:** `Node` **24 → 8 B** (a 3× shrink — it becomes the two count
columns and nothing else), `ContiguousEdge` **48 → 32 B**, `JunctionEdge` unchanged at 16 B. At human
scale that is roughly a **40 % cut** in the accumulator's ~85 MB. ⚠ Re-measure rather than trust that
arithmetic.

### ⚠ `edge_spliced_mass` IS DEAD, AND THAT IS SUSPICIOUS RATHER THAN OBVIOUS — A DECISION FOR THE OWNER

It was added deliberately, with a careful docstring calling it "a per-LINE certified-RNA term,
commensurate with the unspliced mass at the same line". But `calibrate.py:619` builds
`mass_rna_spliced_edge` from `substrate.edge_spliced.COUNT`, not from the mass. So the intent and the
wiring disagree, and the question is which is right.

⭐ **The current wiring is self-consistent and I believe correct.** `mass_rna_edge` is an INCIDENCE
(`(1−f_g)·unspliced_count + spliced_count`); subtracting the spliced *count* leaves
`(1−f_g)·unspliced_count`, which is then multiplied by `q` — and `q` is the UNSPLICED population's own
`mass/count`. Right population, right units, at every step. The spliced mass is not needed because
spliced crossings are removed *before* the conversion rather than converted.

⛔ **So the decision is: was `edge_spliced_mass` built for a job that has since been done a better way
(delete it), or is `calibrate.py:619` using the count where the mass was intended (wire it)?** The units
argument above says the former, but the bank was not added by accident and the owner should confirm before
8 bytes per edge are removed on my reading alone.

## 2c. ⭐ THE REST OF THE LEDGER — everything still owed, ranked by cost

| # | item | cost | risk | reversible |
|---|---|---|---|---|
| 1 | **Dead substrate surface**: `PopulationView.length_sum`, `.mean_length`, `.total_inv_length_sum` | minutes | none — no consumer in `src/` | trivially |
| 2 | **Moment tests**: `contained_moments`/`crossing_moments`/`build_slot_moments` moved to `effective_length.py` WITHOUT their tests, which went with the deleted file | ~1 h | ⛔ untested geometry in a live layer-2 module — this is the highest-value item on the list | n/a |
| 3 | **The four dead banks** (native + schema + reference spec) | ~2 h edit, **hours of re-scanning** | schema digest invalidates every cache | only by re-scanning |
| 4 | **`mass_gdna_edge` units rename** — it is `f_g × unspliced COUNT` while `mass_gdna_node` is a per-fragment count; same prefix, different units | ~1 h | touches `result`, `priors`, `derive`, `track`, goldens | yes |
| 5 | **`second_pass.py:443-445` comment** — measured false in premise, harmless in effect | minutes | none | yes |
| 6 | **`module_census.py`** re-run; fix stale sibling references the purge created | minutes | none | yes |
| 7 | **`UNDOCUMENTED_DEBT`** — 8 scripts on disk and unlisted, pre-dating this campaign | ~1 h | none | yes |

⭐ **Do 1, 2, 5, 6 first — they are cheap, safe and reversible, and item 2 is a real gap.** Item 3 is the
only one with a cost that needs a decision. Item 4 is worth doing but not urgent.

## 3. THE ORDER TO DO IT IN, AND THE GATE AT EACH STEP

| step | gate |
|---|---|
| ~~1. Measure `q_gdna` vs `q_rna`~~ | ✅ **DONE 2026-08-10** — real, placement-driven, ≤ 0.6 pp on the total prior. Record the bound, build nothing (§0.1) |
| ~~2. Purge the channel: python side~~ | ✅ **DONE** — suite green, and byte-identical on every scored calibration field (§1) |
| 3. Purge the banks: native + reference spec + schema (§1.2) | `_accumulator_reference.py` byte-identity; re-measure struct sizes and memory |
| 4. The renames (§2) | ⭐ item 1 **DONE** and byte-identical (§1.3); items 2-3 still owed |
| 5. `module_census.py`, then update `DESIGN.md` §3.1's table and `CLAUDE.md`'s instrument rows | no upward imports, no stale sibling references |

⭐ **Step 2 has a free correctness gate that is worth stating explicitly**: `length_likelihood` defaults
`False`, so removing it must be *byte-identical* on every production path. If any arm moves, the flag
was not the only thing the code was doing — and that is a finding, not a merge conflict.

⛔ **Step 1 comes first and is not optional.** It is the only open question in the prior chain, it is
measurable today from the oracle caches, and it is the one thing in this document that could change
what the tool computes rather than how it reads.
