# Node-state representation: fractions, counts, densities, and the vestigial variance

**Status:** design note (2026-06-23). Pre-implementation. Prompted by the question *"is the per-node
solver operating in fraction space, and should the stored node solution be counts or densities instead
of fractions?"* This note records the careful code audit that answers it and scopes the cleanup.

**TL;DR.** The solver is sound: cross-node messages already travel as **density**, precision as **counts**,
output as **mass**. The fraction is only the per-node solve latent. A careful audit then surfaced two things:
(1) the per-node variance machinery (`NodeBelief.var_*`, `NodeDeconv.gdna_frac_var`) is **provably dead** —
it never affects any point estimate or any output — and should be retired; (2) the per-node *composition*
is the one quantity that is **face-invariant**, and the fraction is its natural representation — so replacing
the stored fraction with a single count or density would be a **regression on boundary nodes** (which have two
faces). The recommended cleanup therefore **retires the dead variance** and **reframes** the fraction's role in
the docs/types, but does **not** re-represent the composition. The user's underlying instinct — *density is the
node↔node currency, counts carry the power* — is already true in the code; this note shows where, and why the
composition stays a fraction.

---

## 1. The three currencies, by layer (what the audit found)

| Layer | Currency | Where | Status |
|---|---|---|---|
| Per-node solve latent | **fraction** `(f₊,f₋,f_g)` on the 2-simplex | `simplex._simplex_lattice`; the strand mixture `p=½f_g+κf₊+(1−κ)f₋` | correct — the strand likelihood is intrinsically a proportion |
| Node↔node messages | **density** `ρ = mass/eff_len` | `_message` (from `ρ_src`), `node_densities`, `_edge_varmean`, `_gdna_seed_estimate`, `fit_enrichment_transfer`, `_global_logprior` | correct — this is the "true currency" the question asks for; already implemented |
| Statistical power / precision | **count** | strand `N·(2κ−1)²`; messages `N_src=ρ²·τ_ρ`; global `n_glob=ρ²/σ²_prior` | correct — the `(M/E)²` Jacobian is gone (count-space refactor shipped) |
| Calibration output | **mass** | `CalibrationResult.mass_*`; `assemble_priors` → per-locus counts + eff-len | correct — no per-node fraction reaches quant/EM |

No site differences, averages, or compares a fraction across two nodes without first converting to density
(verified by an exhaustive sweep). **The fear that "fractions are the message currency" is not what the code
does.**

---

## 2. The dead-variance audit (confirmed, with the full reader/writer map)

The per-node variance state — `NodeBelief.var_pos`, `NodeBelief.var_neg`, `NodeBelief.var_gdna`, and
`NodeDeconv.gdna_frac_var` (produced by `simplex_sweep._fg_var`) — is a leftover from an earlier design where
`1/var` was meant to be the message precision. The count-space message refactor replaced that with the
`var~mean` curves (`gdna_vm`/`rna_vm`) + the effective counts `N_src`/`n_glob`. The variance is now computed and
threaded but **never read as precision, and never reaches any output**.

**Complete reader/writer map (repo-wide grep, src + tests + scripts):**

- `NodeBelief.var_pos`, `var_neg` — written in `_type_belief` (init) and passed through `node_sweep`
  (`bp_solver.py:997`, `:1060–1061`). **No reader anywhere** except `tests/calibration/test_bp_solver.py`.
  Pure pass-through. `node_densities` ignores all variances.
- `NodeBelief.var_gdna` — written in init and updated each pass (`var_g[i] = _fg_var(...)`,
  `bp_solver.py:992`). Its **only** escape is `chain_region_deconv` (`bp_solver.py:1087`) → `NodeDeconv.gdna_frac_var`.
  It is **not** read by the sweep's message/precision/`has_msg_nbr` logic.
- `NodeDeconv.gdna_frac_var` — written in `_solve_nodes` (`simplex_sweep.py:182`), `chain_region_deconv`
  (`:1092`), `chain_boundary_side_deconv` (`:1124`, zeros). Read in **exactly one place**: `_type_belief`
  (`bp_solver.py:278`), where the *strand-only init* deconv's `gdna_frac_var` seeds the G2 init `var_g`.
  That is circular within init. It does **not** appear in `result.py`, `priors.py`, `derive.py`,
  `estimator.py`, `pipeline.py`, or `calibrate.py` — confirmed: `calibrate.py:188–193` reads only
  `.gdna_mass`/`.rna_mass` off the region/boundary deconvs.
- `_fg_var` — called only in `_solve_nodes` and `node_sweep._solve` for the above. Also imported by several
  `scripts/debug/*` for their own diagnostics (keep the helper; it is a legitimate inspection primitive).
- Unrelated namesakes (NOT this state, do not touch): `_SideQuantities.count_gdna_frac_var` (count-module
  Poisson floor, used by `diag_count_posterior_calibration.py`) and `node_pair_fit_over_iters.py`'s
  `median_raw_var_gdna` (a var~mean residual diagnostic).

**Why retiring it is output-preserving (bit-identical).** The point estimates `f_g`/`f₊`/`f₋` come from
`_fg_median`/`_axis_mean`, which do **not** depend on `_fg_var`. The variance never feeds back into the
fractions, the messages, the global, or `has_msg_nbr`. Therefore removing the variance fields changes **no**
mass in `CalibrationResult` and **no** golden. The variance's only "consumers" are the assertions in
`test_bp_solver.py` (the G1=0 / G2-shares-`f_g`-var / G3=∞ convention) — those tests are testing a convention
that nothing relies on, and they get removed/rewritten as part of the cleanup. **Note** the stale comment block
`bp_solver.py:236–245` ("var=0 locked / var=∞ no-info ... the sweep's precision `1/var`") is aspirational and no
longer describes the code — it should be deleted with the fields.

---

## 3. Should the *composition* be stored as counts or densities instead of a fraction?

This is the substantive design question, and the careful audit changes the earlier (naive) recommendation.

**The decisive fact: a node's composition is face-invariant; its mass and density are not.**
- A **region** node presents one geometry both ways: one mass `M`, one eff-len `E`. Here a count triple
  `(n_g,n₊,n₋)=f·M` or a density triple `ρ=f·M/E` is unambiguous.
- A **boundary** node is solved from its *summed* two-side mass, but it has **two faces** with different
  crossing mass `M_left≠M_right` and different eff-len `E_left≠E_right`. Its **composition is one number**
  `f_g` (the same on both faces — that is the whole point of a single per-node solve), but its **mass** is two
  numbers (`f_g·M_left`, `f_g·M_right`) and its **density** is two numbers (`f_g·M_left/E_left`,
  `f_g·M_right/E_right`).

The boundary projection makes this concrete: `chain_boundary_side_deconv` (`bp_solver.py:1111–1126`) takes the
boundary's `f_g` and **re-bases** it onto each flanking region's side mass (`side_fg · side_mass`). That is an
inherently *fractional* operation — same composition, different mass base. `node_densities` likewise multiplies
the one `f_g` by the per-face `(M,E)` to emit the two per-face densities the messages need.

**Consequences for each candidate primary representation:**
- **Fraction (current).** One number per node, face-invariant. Mass = `f·M_face`, density = `f·M_face/E_face`,
  both one multiply away, per face. Minimal and exact. The boundary projection and `node_densities` are direct.
- **Counts.** To keep one value per node you'd store `n_g=f_g·M_total`. Every downstream use then divides back
  out: `node_densities` needs `n_g/M_total·M_face/E_face` (must now thread `M_total` for boundaries), and the
  side projection needs `n_g/M_total·M_side`. The `n_g/M_total` round-trip is **not** bit-identical to the
  original `f_g` (floating-point), so it forces a golden regen for ULP-level drift while adding an `M_total`
  field and a division. Net: more state, more ops, a re-derived fraction anyway.
- **Density.** Same problem, worse: density is *also* face-dependent on a boundary, so a single stored density
  cannot represent the node; you'd store two per-face densities (double the state) — which is exactly what
  `node_densities` already derives on demand.

**Therefore: keep the composition as a fraction.** It is the unique face-invariant per-node state; counts and
density are correctly kept as *derived, per-face views* (mass for output, density for the wire). This does **not**
conflict with the "operate in count space" goal — that goal is about *precision*, and precision is already in
count units (`N_src`, `n_glob`, the strand Fisher info), independent of how the composition is stored.

---

## 4. Recommended cleanup (scoped, low-risk)

**A. Retire the dead variance (the unambiguous win — bit-identical).**
1. `NodeBelief`: drop `var_pos`, `var_neg`, `var_gdna`; delete the stale `236–245` variance-convention comment.
2. `NodeDeconv`: drop `gdna_frac_var`.
3. `_type_belief`: stop reading/setting `gdna_frac_var`/`var_*`; it keeps setting the fractions from the
   strand-only deconv exactly as now.
4. `node_sweep`: drop the `var_g` thread (the `_solve` call no longer writes a variance; `_fg_var` is no longer
   called in the production path — keep the helper in `simplex_sweep` for the debug scripts).
5. `chain_region_deconv`/`chain_boundary_side_deconv`/`_solve_nodes`: drop the `gdna_frac_var=` argument.
6. Tests: remove/rewrite the variance-convention assertions in `tests/calibration/test_bp_solver.py`
   (`~lines 162–302`); keep the `node_densities` and `_fg_median` tests (update the `NodeBelief(...)`
   constructors to the slimmer signature).
7. Verify: full suite + goldens green with **no** `--update-golden` (the proof that it is output-preserving).

**B. Reframe the fraction's role (clarity, no behavior change).**
- At the `NodeBelief` definition, document the layering explicitly: *the belief stores the node's
  face-invariant gDNA/RNA composition (a fraction); the node↔node message currency is density
  (`node_densities`); the calibration output is mass (`NodeDeconv`/`CalibrationResult`); precision is carried in
  count units (the `var~mean` curves + `N_src`/`n_glob`), not in this struct.* This kills the recurring
  "are fractions the currency?" confusion at the source.
- Consider renaming the `NodeDeconv` fraction fields to read as *reported composition* (debug/inspection) vs
  the consumed `*_mass` — optional, only if it doesn't churn the debug scripts.

**C. Do NOT migrate the composition to a count/density triple.** Justified in §3 (the two-faced boundary). If a
future need arises (e.g., a per-face solve), revisit — but it would *add* per-face state, not replace the
fraction.

**Out of scope / unchanged:** the lattice (`_simplex_lattice`), the strand likelihood, the message/global
count-space precision, `node_densities`, `assemble_priors`, the EM. None of these touch the variance and none
depend on the composition being a fraction vs a derived view.

---

## 5. Open decision for the user

The dead-variance retirement (§4.A) is a clean, output-preserving win and is recommended unconditionally.
The §3 finding **overturns** the earlier "migrate to a count triple" line item: the careful audit shows the
fraction is the correct face-invariant node state, so the recommendation is to **reframe, not re-represent**.
If you still want a count/density primary despite the boundary cost (extra `M_total` state + a golden regen for
ULP drift + a re-derived fraction at the boundary projection), that is a deliberate choice — say so and §3 lists
exactly what it entails. Otherwise the plan is §4.A + §4.B.
