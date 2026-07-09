# Two-component {RNA, gDNA} calibration rewrite — implementation plan (v2)

**Status:** plan, ready for review (2026-07-08). Supersedes `two_component_message_passing_redesign.md`
(the design contract) and the first-pass `two_component_implementation_plan.md`. This is a fresh, complete
document; it folds in a v0.6.4-vs-branch code audit and an adversarial design review that corrected two
material errors in the first pass (see §0 and §5).

**Baselines.** v0.6.4 release = `28961cb2` (read via `git show 28961cb2:<path>`). Design-in-progress HEAD =
`7bb58c5b` on branch `calib-disagreement-shrinkage-v1`. HEAD is **3 commits past** the v0.6.4 tag; the module
_set_ is identical. All line numbers are HEAD unless prefixed `v064:`.

**Target modules.** `calibration/bp_solver.py` (the sweep — the bulk of the work), `calibration/calibrate.py`
(orchestration + the new background loop), `config.py` + `cli.py` (flag removal), plus an index-format audit
(`calibration/regions.py`, `calibration/region_arrays.py`, `index.py`) and the dev-script/test fallout. The
per-node solver `calibration/simplex_logodds.py` is **leveraged verbatim** (see Q3). **No C++ changes.**

**Governing principle.** `CALIBRATION_ARCHITECTURE.md` — the count-zero-information axiom: a fragment count
carries no intrinsic gDNA/RNA signal; it enters only as the overdispersed Fisher information of the strand
likelihood. A node's composition is set by three sources: the strand likelihood, the cross-node messages, and
the global/background prior.

**Code-organization principle (user directive — shapes the whole rewrite).** *Leverage the index; do not
re-derive structure at solve time.* The v6 index (`INDEX_FORMAT_VERSION=6`, `boundaries.feather`) already
precomputes the structural facts the solver needs — `is_tss / is_tes / is_splice_junction / genomic_sj_strand`
per boundary (`regions.py:70-80`) and `signature` / `mature_eligible_{pos,neg}` per region. Yet the solver
re-derives equivalents at runtime from signature-bit arithmetic and `if`-statements
(`build_node_geometry`'s `any_exon_l = (sig & (BIT_EXON_POS|BIT_EXON_NEG))` and `_spliced_faces`,
`bp_solver.py:164-196`; the TSS/TES zero-RNA sinks in `init_beliefs`; the intergenic/intron classification in
`_floor_estimate`). The rewrite should **read those precomputed columns** instead — fewer branches, one source of
truth, and the index becomes the clean structural substrate it was built to be. This is why the v6 index is
**kept and leveraged, not reverted** (§7); the format is a genuine forward improvement in organization, and this
principle is how we avoid the "complex flag logic and if-statements to re-derive" the format was meant to retire.

---

## 0. Verdict: are we going in circles?

**No — the rewrite changes production in four deliberate ways that v0.6.4 never had. But the first-pass plan
mis-stated *which* changes are net-zero, and that error would have made validation inconclusive. The corrected
picture:**

Production today (`sweep_disagreement_shrinkage=False`, the default at `config.py:290`) is **byte-identical to
v0.6.4**: it runs the lumped RNA message with the **legacy `σ²_edge` precision** (`bp_solver.py:1444-1448` /
`1490-1494`) and the capped intron floor. `git show 28961cb2:config.py` has no shrinkage flag at all — v0.6.4
**is** the σ²_edge path. The entire +446-line `bp_solver` branch delta (the 5-term mature/nascent overlay + the
A/B scaffolding) lives behind that default-off flag and **has never run in a shipped build**.

So the work splits cleanly into one net-zero part and four genuine forward changes:

**Net-zero for production (pure dead-code cleanup — removing what production never ran):**
- The 5-term mature/nascent overlay: `nasc_p/nasc_n` running beliefs (`1310-1314`), the mature sub-messages
  (`1398-1428` / `1454-1477`), the `mat_elig` gate + `_chain_mature_eligible` (`903-936`), `pnasc_*`/`w_str`
  (`1288-1296`). The `NodeBelief` was already 3-term `{f_pos,f_neg,f_g}` at v0.6.4 and never changed —
  "collapse to 3 terms" is a *no-op on the belief*; the work is entirely in `_scan`'s message arithmetic.
- The shelved A/B variants: `adjacent_disagreement_local` (`966-1010`), `adjacent_component_disagreement_variance`
  (`1013-1062`), and the `sweep_disagreement_pass2` knob. Never promoted to default; pure deletion.

**Four forward changes that make the rewrite ≠ v0.6.4 (each needs its own validation gate):**
1. **Precision basis: `σ²_edge` → belief-free `σ²_imp + 1/n_src`.** This is the change the first-pass plan wrongly
   called "no behavior change." Production runs `σ²_edge` (`1447`), a squared message-vs-belief residual measured
   against the **global-prior-inclusive** local belief `_local_solve(global_lp)` (`1237`) — an echo chamber that
   suppresses legitimate nascent messages into floor-pinned introns. Replacing it with `σ²_imp + 1/n_src`
   (`adjacent_disagreement_variance:951-963`; the `pr = n_src/(n_src·σ²_imp + 1)` form) removes that circularity.
   **This moves the shipped default output.** It is the motivation for the whole disagreement-shrinkage line of
   work, but it has never been the production default — so §4/§6 gate it against the σ²_edge baseline *by itself*,
   before anything is layered on top.
2. **The crux spliced-sink attenuation (Q1).** A per-boundary precision down-weight that attacks the confirmed
   exon→intron through-relay bleed. New mechanism; a wager, not a defect-fix (§1 Q1).
3. **The intron-density background likelihood: uncapped, intergenic-first, softly iterative (Q2).** v0.6.4 caps
   the floor at `_GLOBAL_STAB_PREC` (`783`, too weak to hold the bleed) and lumps intergenic+intron into one
   self-contaminating background (`s2_bg = s2_floor`, `736`). New mechanism (§1 Q2).
4. **A4 exon-nascent deferral** — nascent over an exon has no home in either version; deferring it to messages
   (recovered by the forward+backward `_comb`, `1505-1511`) is new.

**The one way this becomes a real circle:** ship the lumped RNA message *without* landing the crux (2) and the
intron-floor rewrite (3). That is lumped-RNA + `σ²_imp`, which reintroduces the through-relay bleed on top of a
precision basis that trusts messages *more* than σ²_edge did — arguably worse than v0.6.4 in the bleed regime.
The sequencing (§4) forbids regenerating goldens on any such intermediate tree.

---

## 1. The three resolved open questions

### Q1 — capture on the spliced sink → a per-boundary precision likelihood (a wager, well-posed)

**Resolution — resolves *against* the redesign doc §3 as written.** The doc's `M_mat = S·(E_exon-body/E_spliced)`
is a *global* eff-length ratio hard-subtracted from the exon→intron transfer. That is wrong on the user's own
argument: spliced-vs-unspliced capture differs **per capture probe** — a probe spanning the junction captures
spliced preferentially; a probe inside the exon captures both equally — and the probe geometry is unknown, so
**no global ratio exists and each boundary is an independent observation** (user directive). It also carries a
mode/precision defect: inflating the subtraction ~66× would clamp the message mode to the min-observable floor
`max(rho, 1/erd)` (`1440`) while the precision `n_src` stays at full count — a fully-absorbed exon would still
pin the intron at full confidence.

**What we do instead.** Keep the *measured* subtraction `rho_mat_dst = SPd[i]/ESPd[i]` (`1438`) as the hard,
capture-agnostic, per-boundary evidence. Express the *unmeasured* exon-body-mature excess as a **per-boundary
precision down-weight** on the forward RNA message — a variance addition, mode/precision-consistent by
construction (it can never clamp a mode at full precision). On any edge whose destination face carries spliced
mass (`SPd[i] > _EPS` — structurally an exon→boundary sink edge; regions carry `spliced_*=0`, so the gate is
automatic via `build_node_geometry`), replace

```
pr = n_src / (n_src·sig_rna_e + 1)                       # 1443 / 1489
```

with

```
pr = n_src / (n_src·(sig_rna_e + sink_var_i) + 1)
```

where `sink_var_i` is a boundary-local log-density variance built **only from that junction's own spliced
observation** — the squared log gap between the crossing RNA density and the measured mature density `S/E_spl`.
Plus the `−log1p(−f)` Jeffreys parsimony already at `_global_logprior:805`, so the presumption defaults to
mature-sink and only the excess the intron's own background likelihood supports crosses as nascent (recovered by
the existing forward+backward `_comb`, `1505-1511` — the A4 path).

**This is a wager, correctly gated — not a "fix to a v0.6.4 defect" (correction from the first pass).** The
down-weight lowers the message *precision* but not its *mode* (`mo` at `1440` still carries the full exon-body
mature density). It works through the **2-hop chain geometry** (`node_chain.py:89-97`: regions never neighbor
regions — the interleave is B-R-B-R, so exon→intron is always `exon-region → boundary → intron-region`):

- Q1 fires on **hop 1** (exon→boundary): `sink_var ≈ (log 66)² ≈ 17`, so `pr ≈ 1/17 ≈ 0.06` — the exon→boundary
  RNA message is effectively silenced against a boundary whose local precision is floored by the stability prior
  (`_GLOBAL_STAB_PREC=1.0`): a ~5% log-space pull leaves the boundary belief ≈ its local solve.
- The **hop-2 edge** (boundary→intron) has `SPd[intron]=0`, so it runs at **full precision** — but this is
  *correct*: the boundary's own belief is a legitimate nascent signal, because a contiguous unspliced read
  straddling a splice site is pre-mRNA/nascent (mature would be spliced and booked to `mass_spliced`;
  `substrate.py:154`). So Q1 attenuates the mature contamination *entering* the boundary; the clean boundary
  belief relays onward.
- **The last-line bleed-stopper is unchanged**: whatever survives the ~5% residual pull enters the intron at
  full precision and is countered only by the intron-density floor (Q2) — which is exactly the mechanism MEMORY
  records as +12–13% worse calib / ~2× phantom bleed on the 16-suite. **Q1 reduces the intron's input; it does
  not strengthen the defense at the intron.** That defense is the intron-floor rewrite (Q2), and together they
  are the plan's central bet (§5.2).

**Well-posedness (bug fixed vs the first pass — both operands, not just one).** `sink_var_i` takes a log of a
ratio of two densities; **both** can blow up:
- `rho` at `1439` is *post-subtraction* (`n_nasc/er + n_mat/esp − rho_mat_dst`) and can be ≤ 0 or below the
  observable floor → clamp its density input to `max(rho, 1/erd)` (the same floor `mo` uses at `1440`).
- `S/E_spl` is the *more likely* blow-up: as `E_spl → 0`, `log(S/E_spl) → +∞`. Floor `E_spl` with the same
  `min_eff_length` guard used elsewhere, and floor `S` at the observable minimum.
- Explicit short-circuit: `_sink_logdensity_gap` returns `0.0` immediately when `SPd[i] ≤ _EPS` (mature density
  `S/E_spl → −∞` otherwise), *before* any log arithmetic — belt-and-suspenders with the upstream `SPd[i] > _EPS`
  edge gate, so the helper is safe to call unconditionally.

**Honest framing (no over-claim).** `sink_var_i = (log ρ_cross − log(S/E_spl))²` is a squared belief-free
log-residual. It is *not* a sampling variance like `_poisson_moment_var` (`883`), and it is structurally the same
"surprise" shape as the `σ²_edge` we are retiring — the crucial difference is that it is **belief-free** (built
from the boundary's own spliced measurement, not from a prior-inclusive local belief), so it does not have the
echo-chamber circularity. We are **not** claiming it is derived from first principles; the squared-log-gap form is
a principled heuristic and its exact estimator is an **OPEN DECISION** (§7), to be confirmed on `toy_msg.py`
before it is locked. No hidden constant is introduced.

**Code basis (all present, reused):** the spliced **source** half-Gaussian floor
`psi -= 0.5·spliced_pos·max(spliced_pos/n − f_pos, 0)²` (`simplex_logodds.py:191-194`); the **sink** subtraction
(`1438`/`1484`); the precision hook (`1443`/`1489`); the per-boundary spliced density and its half-triangle
eff-length via `node_densities` and `NodeGeometry.spliced_*`.

### Q2 — soft iterative background: intergenic-first, grown by deconvolved gDNA mass

**Resolution (user directive: admission is soft, not binary).** Three changes plus an honest convergence story.

1. **Soft admission by deconvolved gDNA mass.** In `_floor_estimate` (`676`) replace the strand-gated background
   density
   ```
   gdna_frac = (1−w) + w·f_g_init;  dens_g = (M/E)·gdna_frac     # w = (2κ−1)²
   ```
   with the **deconvolved** gDNA density `dens_g = f_g_init·(M/E)`. The strand gate is the toy root cause: for
   unstranded data (`w=0`) `gdna_frac=1`, so a nascent-heavy intron's *full* raw `M/E` re-enters the background
   every pass (`toy_msg.py`: 5× excess, true `f_g=0.20`, ~80% nascent, called ≈gDNA). The deconvolved density
   caps each intron's contribution at its own solved `f_g·M/E`, which is exactly the user's "use the intron's
   deconvolved gDNA mass as the background" — soft and continuous, not include/exclude. Apply the same fix to the
   strand-gated seed weight in `_gdna_seed_estimate` (`585-587`). This generalizes what the KDE prior already does
   (`gdna_density_prior.build_training_substrate` trains on each node's `ρ_g = f_g·M/E`); we extend that rule to
   the floor/`s2_bg` estimators, which today do not use it.

2. **Intergenic-first seed.** Round 0 fits the background from **intergenic regions alone** (`node_rtype==0`) —
   gDNA-clean, nascent-free. Both anchors already exist (`_gdna_seed_estimate` pins intergenic `f_g=1` at `583`;
   `_floor_estimate` locks intergenic `f_g_init=1` at `689`). Retire the intergenic+intron lump `s2_bg = s2_floor`
   (`736`) — that lumped background *is* the toy failure.

3. **Explicit outer background loop in `calibrate.py`.** Today the background is rebuilt *inside* `node_sweep`
   (`1174-1182`) from the belief passed in, and "iteration" is a fixed 2-pass (`194`/`226`) with no convergence
   check. Lift the background to an outer loop of 2–3 rounds: round `k+1` grows it by each intron's *deconvolved*
   gDNA mass, stopping on relative `|Δρ_bg|/ρ_bg < config.sweep_convergence_delta` (existing; no new constant).
   **This requires a real refactor, not just a wrapper** (correction from the first pass): `ρ_bg`/`s2_bg` must
   become a value **passed into** `node_sweep` (and `_global_logprior`) with the internal per-call recompute
   suppressed — a new `node_sweep` parameter that the touch list (§2) now records.

**Convergence — the property to study (favorable but NOT a proven contraction; correction from the first pass).**
The idealized map `T(ρ_bg) = [G_intergenic + Σ_introns min(ρ_bg, d_i)·E_i] / E_total` with the intron densities
`d_i` held **fixed** is a Banach contraction (monotone, 1-Lipschitz `min`, intergenic exposure pinned in the
denominator ⇒ `dT/dρ_bg < 1`), approaching its unique fixed point *from below* (intergenic-first ⇒ the most
conservative, nascent-preserving fixed point). **But `d_i` is not fixed:** each round re-solves `f_g` through the
full nonlinear sweep (strand × messages × KDE) and, if the loop wraps the full pass, **retrains the KDE**
(`calibrate.py:222`) — positive feedback that MEMORY documents as the pass-2 message crush
(`node1085_kde_valley_boundary_crush`). So the clean 1-Lipschitz bound does not transfer. Design accordingly:
   - **Keep a damping / relaxation step and an oscillation guard** — do NOT assume "no damping, contraction
     proven." Under-relax `ρ_bg^{k+1} = (1−a)·ρ_bg^k + a·T(ρ_bg^k)` if the raw trace overshoots.
   - Keep intergenic pinned as the lower anchor every round; keep the Jeffreys brake (`805`/`880`).
   - **Loop scope is an EMPIRICAL STUDY, not a pre-commit (user directive, §7).** Compare **pass-1-only** (wrap
     just the background update, retrain the KDE once at the end — cheaper, avoids the KDE-retrain feedback)
     against the **full iterative** pass-1/refit/KDE/pass-2, on the toy + the benchmark ρ_bg trace. This is where
     the convergence risk lives, so measure it rather than assume; pass-1-only is the prior, not the verdict.
   - **Fundamental limit (A2):** density alone cannot separate a truly gDNA-dense intron (CNV / low mappability)
     from a nascent-heavy one; iteration only grows sample size where nascent is sparse. Correctness on
     adversarial unstranded high-CNV introns still rides on strand (`κ≠½`) or the nascent-sparse prior. Document,
     do not pretend to solve.

### Q3 — leverage the existing strand/message unification (do NOT redesign)

**It already exists, is clean, and is unchanged from v0.6.4.** `git diff 28961cb2..HEAD -- simplex_logodds.py
strand_likelihood.py simplex.py` is **empty** (verified). `bp_solver._local_solve` (`1149-1160`) calls
`simplex_logodds._solve_nodes_logodds_all` (`398`), which builds **one** log-posterior ψ over `λ = logit(f_g)`
(`_local_loglik_logodds:168-228`) summing, on the same grid:
- the strand Beta-Binomial `_mixture_strand_loglik` (`189`; count enters only as overdispersed Fisher info —
  count-zero-info),
- the sided spliced floor (`191-194`),
- the Jeffreys strand reference (`196-202`),
- the count-space global gDNA prior (`204-205`),
- the cross-node messages as log-fraction Gaussians (gDNA `212-215`, per-strand RNA `216-222`) + the
  change-of-variable Jacobian (`227`).

Strand precision `(2κ−1)² → 0` for unstranded data, so a `κ=½` node **defers to messages + prior with no
branch** — the unified strand-and-message solve the redesign §7 asked for **already holds**. There is no
strand-vs-message dichotomy to build.

**The only problem is upstream of this solver, and we are already fixing it:** what the *legacy `σ²_edge`*
measures message surprise against — `resid = mo − lfg_loc[i]` compares to the **global-prior-inclusive** local
belief `_local_solve(global_lp)` (`1237`), coupling message trust to the floor and suppressing legitimate nascent
in floor regions. Retiring the `σ²_edge` branches (forward change #1, §0) structurally removes this. **Leave
`_solve_nodes_logodds_all` untouched.** (Non-blocking note: `f_g` is a posterior *median* while `f_pos/f_neg` are
*means*, so the three need not sum to 1 — a deliberate marginal-estimator choice, not a defect.)

---

## 2. Touch list

### 2a. Exhaustive touch table

| File | Function(s) / symbol | Action | What / why |
|---|---|---|---|
| `calibration/bp_solver.py` | lumped RNA⁺ branch `1429-1452`; RNA⁻ `1478-1498` | modify | Make lumped the **sole** RNA path (drop the `_mat_split if`). Add the Q1 per-boundary `sink_var` precision term on `SPd[i]>_EPS` edges, with the `max(rho,1/erd)` clamp + Jeffreys. |
| `calibration/bp_solver.py` | gDNA msg `1358-1370`; RNA `1442-1448`/`1488-1494` | modify | Delete the `else:` legacy `σ²_edge` branches; keep only `pr = n_src/(n_src·σ²_imp + 1)`. (Forward change #1 — a live production change; gate it, §4/§6.) |
| `calibration/bp_solver.py` | `node_sweep` sig dispatch `1245-1263` | modify | Collapse `float\|tuple\|ndarray` polymorphism to one `disagreement_sigma2: float`; remove `_sig_g/_sig_rna/_sig_edge/_shrink_on/_per_edge`. |
| `calibration/bp_solver.py` | `node_sweep` signature | modify | **Add a `background` parameter** (`ρ_bg`/`s2_bg`) and suppress the internal `_floor_estimate`/`_gdna_seed_estimate` recompute (`1174-1182`) — required by the Q2 outer loop. |
| `calibration/bp_solver.py` | `_mat_split` `1269`; 3-channel state `1271-1296`; `nasc_p/nasc_n` `1310-1314`; per-edge sig `1334-1337`; RNA⁺ 3ch `1398-1428`; RNA⁻ 3ch `1454-1477` | remove | Delete the 5-term overlay entirely (restores the lumped `_scan`). |
| `calibration/bp_solver.py` | `mature_nascent_split` kwarg `1084` | remove | Overlay gone. |
| `calibration/bp_solver.py` | `_chain_mature_eligible` `903-936` | remove | Remove the 3-channel **consumer** of `mature_eligible_*` (the mature message gate). The index **columns stay** (kept for the leverage principle / future use); only this dead consumer goes. |
| `calibration/bp_solver.py` | `adjacent_disagreement_local` `966-1010`; `adjacent_component_disagreement_variance` `1013-1062` | remove | A/B variants; never promoted. |
| `calibration/bp_solver.py` | `_floor_estimate` `676-737`; `_gdna_seed_estimate` `585-587` | modify | Deconvolved `dens_g = f_g·M/E`; intergenic-first `s2_bg` (retire `736`); unify the `min_eff_length` guard with `build_training_substrate`. |
| `calibration/bp_solver.py` | `_global_logprior` `740-806`; `_kde_logprior` floor exclusion `1212-1220` | modify | Consume the passed-in `ρ_bg`/`s2_bg`; make the floor-vs-KDE partition **structural** (one density prior + one Jeffreys per node) instead of the post-hoc `kde_lp[floor_mask]=0` at `1220`. |
| `calibration/bp_solver.py` | `_poisson_moment_var` `883-900`; `_adjacent_edges` `939-948`; `adjacent_disagreement_variance` `951-963` | keep | The surviving σ²_imp estimator (Q3 KEEP). |
| `calibration/bp_solver.py` | new `_sink_logdensity_gap` helper | add | The only new helper (§3.1): boundary-local squared-log-gap, clamped; returns 0 when no spliced mass. |
| `calibration/calibrate.py` | module imports `42-43` | remove | `adjacent_component_disagreement_variance`, `adjacent_disagreement_local` imports (else `ImportError`). |
| `calibration/calibrate.py` | pass dispatch `174-212`; sweeps `194`/`226` | modify | Remove `shrink`/`pass2` branching; always compute `σ²_imp`; add the Q2 `_iterate_background` outer loop (§3.2). |
| `config.py` | `sweep_disagreement_shrinkage` `290`; `sweep_disagreement_pass2` `298` + validator `374-377`; `sweep_mature_nascent_split` `306`; comments `292`/`295` | remove | One production path; no flag ladder. Clean the now-stale comments too. |
| `cli.py` | `_ParamSpec`s `625-626`; `--sweep-disagreement-shrinkage` `1214-1223`; `--sweep-disagreement-pass2` `1225-1232` | remove | Match config. (FL-report additions elsewhere in `cli.py` are orthogonal — keep.) |
| `calibration/simplex_logodds.py` | `_solve_nodes_logodds_all` / `_local_loglik_logodds` | keep | The unification (Q3) — verbatim. |
| `index.py`, `calibration/regions.py`, `calibration/region_arrays.py` | `INDEX_FORMAT_VERSION=6`; `boundaries.feather` (`is_tss/is_tes/is_splice_junction/genomic_sj_strand`); `signature`; `mature_eligible_{pos,neg}`; `BoundaryArrays`/`from_boundary_df` | **keep** | **v6 is kept and leveraged, not reverted** (§0 principle, §7). No index rebuild, no fixture regen. The columns are the structural substrate the solver should read. |
| `calibration/bp_solver.py` | `build_node_geometry` `138-196` (`any_exon_l/r`, `_spliced_faces` `164-196`) | modify | **Leverage the index:** route the spliced source/sink from the boundary's precomputed `is_splice_junction`/`genomic_sj_strand` instead of re-deriving exon/junction structure from signature bits. Replaces branches with a column read. |
| `calibration/bp_solver.py` | `init_beliefs` (TSS/TES zero-RNA sinks); `_floor_estimate` node classification `676-737` | modify | **Leverage the index:** read `is_tss`/`is_tes` for the zero-RNA sinks and the region type for intergenic/intron classification from the index, not from runtime signature-bit `if`s. Bounded to swapping re-derivation for reads — must net-remove branches, not add. |
| `pipeline.py` (+13); `sim/locus_sweep.py` (+10); `frag_length_model.py` (~200); `calibration/fl.py` (+72) | — | keep | Orthogonal FL-report work (commit d2a22571); not overlay wiring. Confirm, keep. |
| `scripts/debug/benchmark_ab_report.py` `37`,`39`,`115-138` | modify | Remove the `--sweep-disagreement-shrinkage`/`--sweep-disagreement-pass2` arg construction/parsing (else the harness breaks at runtime). |
| `scripts/debug/pass_trace.py` `79`,`118`,`142` | modify | Drop references to `sweep_disagreement_shrinkage` + `adjacent_component_disagreement_variance`. |
| `scripts/debug/toy_msg.py`, `toy_sweep.py` (scratchpad) | modify | Drop the removed config kwargs from the toy configs (dev instruments). |
| `tests/calibration/test_bp_solver.py` `687-718` incl. `test_three_channel_legacy_path_unchanged_when_flag_off` `714` | modify | **Repurpose, don't blind-delete** `714` (it is a flag-off byte-identity test referencing `mature_eligible_pos`); delete the 3-channel routing cases; update remaining calls to the single-float `disagreement_sigma2`; add a through-relay + iterative-background assertion. |
| `tests/calibration/test_structural_masks.py` | keep | v6 retained → the structural-mask index asserts stay valid. |
| `tests/calibration/test_priors.py` `110-111` | keep | `mature_eligible_*` fields retained; kwargs unchanged. |
| `tests/golden/*` | regen | `pytest --update-golden` — only after the §4 gates pass. |

### 2b. REMOVE ledger (deletion-only; all default-off scaffolding)

`bp_solver.py`: `_chain_mature_eligible` (`903-936`), `adjacent_disagreement_local` (`966-1010`),
`adjacent_component_disagreement_variance` (`1013-1062`), 3-channel state (`1271-1296`), `nasc_p/nasc_n`
(`1310-1314`), RNA⁺/RNA⁻ 3-channel branches (`1398-1428`, `1454-1477`), legacy `σ²_edge` `else:` branches
(`1363-1370`, `1444-1448`, `1490-1494`), polymorphic dispatch (`1245-1263`), `mature_nascent_split` kwarg
(`1084`). `calibrate.py`: imports `42-43`. `config.py`: three flags + validator + stale comments. `cli.py`: two
`_ParamSpec`s + two flag definitions. Dev scripts: the flag references above. Tests: 3-channel routing cases.

**NOT removed (v6 kept + leveraged, §0 principle):** the index format, `boundaries.feather`, `mature_eligible_*`
columns, `BoundaryArrays`, `test_structural_masks.py`. Only the dead 3-channel *consumer* (`_chain_mature_eligible`)
is removed; the structural columns become solver inputs (§2a `build_node_geometry`/`init_beliefs` rows).

---

## 3. New code to write

Two pieces, both by-purpose-named, both grounded in existing machinery.

### 3.1 The spliced-sink message precision (`_scan`, lumped RNA branches)

A pure local computation inside `_scan`; no new call signature beyond one small helper. On the `SPd[i] > _EPS`
edge, after the existing `rho` (`1439`), mode `mo` (`1440`), and count `n_src` (`1441`):

```
# Presumption: RNA crossing a spliced-junction boundary is mature (measured by S), not nascent.
# Down-weight the forward RNA message by THIS boundary's own spliced evidence — a per-boundary
# likelihood (no global E_exon/E_spl ratio, no on/off gate). Jeffreys defaults to mature-sink;
# only the intron's own background-likelihood excess crosses as nascent (recovered via _comb).
sink_var = _sink_logdensity_gap(rho, erd, SPd[i], ESPd[i])   # boundary-local; 0 when no spliced mass
pr = n_src / (n_src * (sig_rna_e + sink_var) + 1.0)
```

`_sink_logdensity_gap` is the only new helper — a few lines: short-circuit `return 0.0` when `SPd[i] ≤ _EPS`;
else clamp the crossing density to `max(rho, 1/erd)` and floor `E_spl`/`S` (the well-posedness fixes above),
compare to the measured mature density `S/E_spl`, return the squared log gap floored at 0. The Jeffreys parsimony
reuses `−np.log1p(−f)` already at `_global_logprior:805`.
**Nothing hard-subtracts an inflated `M_mat`** — the mode keeps the measured-`S` subtraction (`1438`), so no
clamp/precision mismatch. The exact gap estimator is flagged OPEN (§7) for confirmation on `toy_msg.py`.

### 3.2 The iterative background loop (`calibrate.py`)

Lift the background out of `node_sweep` into an explicit, inspectable outer loop, and pass it back in:

```
def _iterate_background(chain, geometry, region_arrays, belief0, kappa, config):
    """Grow the gDNA-background rate from INTERGENIC alone, then add each intron's DECONVOLVED
    gDNA mass (f_g·M/E) softly (min-contribution). Under-relax + stop on relative change
    < config.sweep_convergence_delta (cap config.sweep_max_passes). Returns (rho_bg, s2_bg)."""
```

Round 0 = `_gdna_seed_estimate` on the intergenic seed set only (already `struct_seed`, `583`). Each round runs
the sweep with the current `background`, reads back `belief.f_g`, recomputes `dens_g = f_g·M/E`
(§1 Q2 fix), re-pools with intergenic pinned, under-relaxes, checks the relative-Δ stop. No new constant: the
stop reuses `sweep_convergence_delta`; the cap reuses `sweep_max_passes`. This makes `ρ_bg`/`s2_bg` a first-class
value threaded into `node_sweep` → `_global_logprior`, rather than rebuilt per call. **Loop scope (pass-1-only vs
full iterative) is decided by the empirical study in §1 Q2 / §6**, not pre-committed; the interface above works
for either — pass-1-only simply omits the KDE retrain inside the loop.

---

## 4. Implementation sequencing

Each step leaves a runnable tree and is independently testable. **No `pip install` recompile at any step** —
nothing under `src/rigel/native/` changes; the accumulator already emits `boundary_junction_strand`
(`accumulator.cpp:211` → `substrate.py:188` → the spliced source/sink placement `bp_solver.py:163`),
independent of `boundaries.feather`.

1. **Precision-basis flip — as its OWN gated step (the correction).** Delete `adjacent_disagreement_local`,
   `adjacent_component_disagreement_variance`, the polymorphic dispatch, the `pass2` knob; make
   `disagreement_sigma2: float`; retire the `σ²_edge` `else:` branches so `σ²_imp + 1/n_src` is the only message
   precision. Clean the `calibrate.py:42-43` imports. **Gate (blocking):** calib-vs-oracle on 16 + 24-AMBIG for
   `σ²_imp`-only **vs the current σ²_edge production default** — because this moves the shipped output (§0 #1).
   Record the delta; it must be neutral-or-better before proceeding, and it must be *separable* from the crux.
2. **Remove the 5-term overlay.** Delete the 3-channel state and RNA branches, `_chain_mature_eligible`, the
   `mature_nascent_split` kwarg; the lumped `elif emit_p/emit_n` become the sole `if`. *Test:* solver runs; this
   is now lumped-RNA on σ²_imp (expect the through-relay bleed — that is the point of step 3).
3. **Land the crux (Q1).** Add `_sink_logdensity_gap` + the precision term (§3.1) with the clamp. **Gate
   (blocking):** `toy_msg.py` — (a) a mature-only intron must NOT bleed (recover the confirmed gated behavior),
   *and* (b) a **hop-2-isolating** case (a thin, low-nascent boundary between a high-mature exon and a
   genuinely-gDNA intron) must not bleed at the full-precision entry edge, *and* (c) a nascent-dominant intron
   must still recover nascent.
4. **Iterative background (Q2).** The `_floor_estimate`/`_gdna_seed_estimate` deconvolved-mass fix; intergenic-
   first `s2_bg`; the `node_sweep` `background` parameter + suppressed internal recompute; the `_iterate_background`
   loop (§3.2) with under-relaxation; unify `min_eff_length`; make floor-vs-KDE structural. **Gate:** `toy_msg.py`
   5× excess → nascent recovered; log the `ρ_bg` trace to confirm 2–3-round convergence with no oscillation.
5. **Config/CLI/scripts cleanup.** Remove the three flags + CLI + validator + stale comments; fix
   `benchmark_ab_report.py`, `pass_trace.py`, and the toy scripts. *Test:* `CalibrationConfig` construction/
   validation; the dev scripts import and run.
6. **Leverage the index (v6 kept).** No revert. Wire the boundary columns
   (`is_tss/is_tes/is_splice_junction/genomic_sj_strand`) into `build_node_geometry`'s spliced source/sink routing
   and `init_beliefs`' TSS/TES sinks, replacing the signature-bit `if`-derivation (§0 principle, §2a). *Gate:* the
   refactor is *intended* to be behavior-preserving — assert bit-identical node geometry vs the signature-derived
   path on the toy + a fixture. **If numbers DO move, this is a discovery, not a failure — and the response is
   neither "revert" nor "trust the index blindly."** A divergence means the indexer (`_compute_mature_eligible`/
   `build_boundary_partition`, from the GTF) and the solver's signature-bit derivation disagree on an edge case
   (candidates: coincident opposite junctions `genomic_sj_strand=3`, TSS==TES points). Isolate the diverging
   nodes, verify against the actual transcript coordinates which representation is correct — it may be *either*
   side — fix the wrong one, and document it as a structural-drift **fix**, not a regression. Do not fold such a
   fix into this refactor's golden diff silently; land it as its own reviewed change.
7. **Goldens + full validation.** `pytest tests/ --update-golden` (only after steps 1 & 3 gates pass — never on a
   bleed-reintroducing intermediate); then §6. Mark the superseded design docs
   (`three_component_mature_nascent_design.md`, `mature_rna_channel_design.md`,
   `disagreement_shrinkage_prior_design_v2.md`) as superseded.

---

## 5. Prospective problems & mitigations (prioritized)

1. **Validation-attribution conflation (highest — the first-pass error).** The σ²_imp precision flip is itself a
   live production change (production is σ²_edge = v0.6.4, §0). If it is not gated *separately* from the crux, the
   §6 numbers cannot separate "the crux worked" from "the precision flip moved everything," and the falsifiable
   16-suite bet (below) becomes unattributable. *Mitigation:* §4 step 1 is its own gate vs the σ²_edge baseline.
2. **Losing the stranded-exon 3-channel win (the central wager).** MEMORY: 2-component was +12–13% worse calib,
   ~2× phantom bleed on stranded/clean loci; the 3-channel's per-strand mature gate earned its keep there. The
   rewrite bets crux (Q1) + intron-floor (Q2) recover that **without** the channel — and Q1 only reduces the
   intron's *input*, leaving the intron floor as the sole defense (§1 Q1). *Mitigation:* §6 validates
   calib-vs-oracle on the 16-suite **before the channel is deleted for good**; if stranded-exon accuracy is not
   recovered, this is a **stop-and-rethink gate**, not a merge. Genuinely unproven — the plan's one bet.
3. **The circle.** Overlay removed + RNA reverted to lumped + crux not landed ⇒ the v0.6.4 through-relay bleed on
   a message-trusting precision basis. *Mitigation:* §4 forbids goldens before step 3's toy gate.
4. **`sink_var` well-posedness.** `rho` (`1439`) is post-subtraction and can be ≤ 0 / sub-floor ⇒ `log` undefined.
   *Mitigation:* clamp to `max(rho, 1/erd)` inside `_sink_logdensity_gap` (§3.1).
5. **Iterative-background non-convergence / KDE feedback.** `d_i = f_g·M/E` is re-solved nonlinearly each round;
   a full-pass loop retrains the KDE (the node-1085 crush). The clean Banach bound does not transfer. *Mitigation:*
   pass-1-only loop scope, under-relaxation, intergenic lower anchor, Jeffreys brake, relative-Δ stop, oscillation
   guard; log the trace (§1 Q2, §3.2). *Fundamental limit (A2):* adversarial unstranded CNV/mappability introns
   are unresolvable by density alone — documented, not solved.
6. **`s2_bg = s2_floor` silently surviving (`736`).** If the region set is not actually changed, the fix is
   illusory. *Mitigation:* step 4 asserts intergenic-first in code + the 5×-excess toy gate.
7. **Double-Jeffreys.** Floor nodes must not carry `−log(1−f)` from both `_global_logprior:805` and
   `_kde_logprior:880`. *Mitigation:* step 4 makes the floor-vs-KDE partition structural (replace the post-hoc
   `kde_lp[floor_mask]=0` at `1220`).
8. **`min_eff_length` applied only to the KDE.** A sub-fragment intron (`E→0`) blows up `1/E` into the floor
   (`_floor_estimate` has no guard; `build_training_substrate:150` does). *Mitigation:* unify (step 4).
9. **Index-leverage refactor drift (medium).** Reading boundary columns (`is_splice_junction`/`genomic_sj_strand`/
   `is_tss`/`is_tes`) instead of re-deriving from signature bits could subtly change node geometry if the index
   annotation and the bit derivation disagree at edge cases (coincident opposite junctions `genomic_sj_strand=3`,
   TSS==TES points). *Mitigation:* §4 step 6 asserts bit-identical node geometry vs the signature-derived path on
   the toy + a fixture before the derivation code is deleted — it is a simplification, so numbers must not move.
10. **Dev-script breakage (low, but a footgun).** `benchmark_ab_report.py`/`pass_trace.py` construct/parse the
    removed flags and `import` removed symbols. *Mitigation:* step 5 fixes them alongside the CLI.
11. **Validation confound (methodology).** Better calibration can look worse post-EM (the separate eff-length
    collapse bug). *Mitigation:* judge calib-vs-oracle only (§6), never post-EM gene counts.

---

## 6. Validation plan

**Toy gates (`scratchpad/toy_msg.py`; alignment-based oracle; must include ample single-strand nodes — never
AMBIG-in-isolation, per the memory directive):**
- (a) **Nascent recovery:** an intron with a genuine density excess over the intergenic-first background → nascent
  (the current failure: 5× excess, true `f_g=0.20`, must not read ≈gDNA).
- (b) **No through-relay (hop-1):** a mature-only intron must not inherit exon-body mature (recover the gated
  behavior).
- (c) **No through-relay (hop-2):** a thin, low-nascent boundary between a high-mature exon and a genuinely-gDNA
  intron — the full-precision entry edge must not bleed (the case Q1 does *not* directly gate).
- (d) **No exon-nascent orphaning:** nascent identified at the intron, propagated back to the unstranded exon via
  `_comb` (A4).
- (e) **Convergence:** log the `ρ_bg` trace; assert relative-Δ < `sweep_convergence_delta` within
  `sweep_max_passes` (expect 2–3 rounds), no oscillation.

These gates are also where the two §7 empirical studies are run: gates (a)/(b)/(c) instrument the **`sink_var`
behavior** study (attenuation strength + gap-form comparison); gate (e) plus the benchmark ρ_bg trace run the
**pass-1 vs iterative loop** study (and the damping question). Neither form is locked until its gate decides.

**Benchmarks — calib-vs-oracle, NOT post-EM** (`scratchpad/calib_oracle_accuracy.py` corrected; the
`calibration-benchmark` skill's net fragment-flow metric):
- **Step-1 attribution gate:** `σ²_imp`-only vs the current `σ²_edge` production default (§4 step 1) — establish
  the precision-flip delta in isolation.
- **16-suite** (`~/Downloads/rigel_runs/quick_3to1_5mb`): the decisive stranded/clean-locus test for the wager
  (§5.2). Must **beat or match** the 3-channel on Σ|err|, intron-bleed, AMBIG-exon.
- **24-AMBIG** (`ambig_dense_10mb`): the 2-D-solver stress suite; 2-component historically ties/wins.
- Report **per-class** (stranded exon / AMBIG exon / intron nascent) so the wager is falsifiable, not
  aggregate-masked.
- **Directional confusion cell — the decisive metric for the wager:** track the **stranded-exon → intron bleed
  rate** specifically (mass leaking from a clean, highly-expressed stranded exon into its adjacent introns), not
  just symmetric per-class Σ|err|. This is the exact localized regression the 3-channel's per-strand mature gate
  used to prevent (memory: the 16-suite +12–13% / ~2× phantom bleed). Aggregate calibration can look perfect while
  this cell spikes. If it does, Q1's down-weight is not aggressive enough — re-examine the `sink_var` scale (the
  §7-1 study) **before** the 3-channel code is deleted for good. This cell, not the aggregate, is the §5.2
  stop-and-rethink trigger.

**Determinism:** `total_threads=1`, confirm bit-identical calibration across processes (the C++ parallel-FP
nondeterminism is orthogonal, independently tracked). Collapsing to one path should reduce, not add, the
nondeterminism surface.

**Goldens:** regenerate only after the step-1 and step-3 gates pass; review the diff for the expected direction
(less intron gDNA over-call; no new no-gDNA false positives).

---

## 7. Scope guard & open decisions

**Explicitly NOT doing:**
- **No mature/nascent channel in calibration.** Calibration emits per-region gDNA-vs-RNA mass only; the per-locus
  EM reconstructs mature/nascent (`n_t+1` components). Output schema (`CalibrationResult`, `assemble_priors`, the
  two per-locus Dirichlet scalars) and the EM are **unchanged**.
- **No hard `M_mat = S·(E_exon/E_spl)` subtraction** (redesign doc §3) — replaced by the per-boundary precision
  likelihood (§1 Q1).
- **No `σ²_edge` surprise-anchor, no per-component/per-edge variants, no A/B flags left behind** — one belief-free
  `σ²_imp + 1/n_src` path.
- **No redesign of the per-node ψ solve** (`simplex_logodds`) — leveraged verbatim (Q3).
- **No C++ changes.**
- **No eff-length-collapse fix** — a separate, confirmed bug that inverts post-EM benchmarks (hence calib-vs-oracle
  validation). Out of scope here.
- **No post-EM tuning; no new magic constants** — every introduced quantity is either a node's own observation
  (`sink_var`) or reuses an existing knob (`sweep_convergence_delta`, `sweep_max_passes`).

**EMPIRICAL STUDIES (user directive — decided by measurement during implementation, not pre-committed):**
1. **`sink_var` behavior** (§1 Q1, §3.1). Study the empirical behavior of the per-boundary down-weight before
   locking its form — raw squared log-density gap vs a Poisson-corrected gap (à la `_poisson_moment_var`), and how
   strongly it attenuates across real boundaries. **Also study whether an unbounded `sink_var` (which can reach
   ≈20 at a high-mature-exon / silent-intron seam, driving `pr → 0`) over-silences the hop-1 message and leaves
   the boundary hyper-dependent on the stability prior — i.e. whether saturation is needed.** If it is, the
   saturation must be *derived* (e.g. from the observable-count Fisher information, capping precision at what the
   count can support), **not a hardcoded ceiling** — a bare cap would be exactly the kind of magic number this
   effort forbids. Instrument on `toy_msg.py` + the benchmark; the §4-step-3 gate is the arbiter.
2. **Pass-1 vs iterative background loop** (§1 Q2, §3.2). Study pass-1-only vs the full iterative loop empirically
   (the ρ_bg convergence trace + calib-vs-oracle), including whether under-relaxation damping is needed. Do not
   assume; measure. Pass-1-only is the prior, not the decision. **Study via a throwaway harness, NOT a shipped
   switch:** implement pass-1-only as the mainline default (its interface, §3.2, already subsumes both — the
   full-iterative variant is just the KDE retrain moved inside the loop) and test full-iterative via a local scratch
   modification. Commit only the winner as the single path. A permanent `config` switch here would recreate exactly
   the A/B flag-cruft this rewrite exists to delete (the `sweep_disagreement_*` ladder is the cautionary tale).

**RESOLVED:**
3. **Index format: KEEP v6 and LEVERAGE it (§0 principle).** Do not revert. The v6 index is the cleaner,
   future-facing structural substrate; the rewrite reads its precomputed boundary/region columns instead of
   re-deriving structure via signature-bit `if`-statements (§2a `build_node_geometry`/`init_beliefs`; §4 step 6).
   `mature_eligible_*` and `boundaries.feather` are retained (only the dead 3-channel consumer
   `_chain_mature_eligible` is removed). This directly serves the goal of leaner, index-driven code over flag
   logic — and is itself part of why the rewrite improves on v0.6.4 rather than circling back to it.
