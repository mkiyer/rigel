# Getting the accumulator's information into the solver — and back out correctly

    Status:   PROPOSED, 2026-07-31. Nothing here has landed.
    Supersedes: `docs/calibration/S5_DESIGN_LOG.md` §2 Phase 2 rows 5 and 7 (S5.a2). Row 4 (S5.g / A7) is already done.
    Companion: `docs/accumulator/DESIGN.md` (what is deposited) · `docs/accumulator/NODE_DENSITY_DERIVATION.md` (why)

⭐ **Read §1 before touching anything.** It is the whole theory in one page and every step below is a
direct consequence of it. §2 states the two defects and the evidence for them. §§4–7 are the executable
steps.

⛔ **This plan does NOT change the deposit rule, the strand likelihood, or the message passing.** If you
find yourself editing `bam_scanner.cpp`, `strand_likelihood.py` or `bp_solver.node_sweep`, stop — you
have left the plan.

---

## §0 STATUS AND ORDER

| step | what | size | status |
|---|---|---|---|
| **P0** | Measure the hole. No code changes | ¼ day | ✅ **DONE 2026-07-31** — premise confirmed; see §4 and `docs/WIP.md` |
| **P1** | The units fix in `assemble_priors` — a **bug** | 1–2 days | ✅ **DONE 2026-07-31** — see `docs/WIP.md`. ⛔ T3 caught the old prior EXCEEDING the whole library |
| **P2** | Wire the length likelihood into the per-node solve — the **information fix** | 3–5 days | ⚠ **BUILT + GATED OFF 2026-07-31.** Mechanism proven (blind mass 100 % → 0 %); ⛔ blocked on the FL pools — see `docs/WIP.md` |
| **P3** | Change the node deposit weight `L → A` | 2–3 days | ⚠ **decide after P2, not before** |

⚠ **Do P1 before P2.** P2's end-to-end A/B is read through P1's code; leaving the units error in means
every P2 number is measured through a biased lens.

⚠ **Regenerate `tests/golden/` ONCE, after P2** — not after P1, and not twice.

---

## §1 THE THEORY, FROM THE GROUND UP

### 1.1 What the solver is trying to do

At every object — a **node** (a genomic interval) or an **edge** (a 0-bp line between two adjacent
nodes) — answer one question:

> Of the fragments here, how many are gDNA and how many are RNA?

Two unknowns: `rho_g` and `rho_r`, the number of fragments of each component **starting per genomic
base**. Two unknowns need **two independent equations**.

### 1.2 Where equations come from

A fragment of length `w` lands on a given object from a specific number of start positions. Call it the
**opportunity**, `A(w)`. There are only two opportunity formulas in the entire design, both verified by
exact enumeration (`scripts/design/node_density_derivation.py`, T0):

```
    contained in a node of length ell :   A(w) = max(0, ell - w + 1)
    crossing a 0-bp line              :   A(w) = max(0, w - 1)
```

*(A third, `A(w) = max(0, w - ell - 1)`, is the node-**spanning** population. It is stored and unused;
see §8.)*

The accumulator deposits some weight `h(w)` each time a fragment qualifies. Because fragment starts are
a Poisson process, the expectation of any such accumulator is:

```
    E[ Sum h ]  =  rho_g · SUM_w f_g(w)·A(w)·h(w)   +   rho_r · SUM_w f_r(w)·A(w)·h(w)
                =  rho_g ·          a_g              +   rho_r ·          a_r
```

where `f_c(w)` is component `c`'s fitted fragment-length pmf.

**Every stored channel is linear in `rho`, always.** Each choice of `h` gives one row `(a_g, a_r)` of a
linear system. That is the whole theory.

| deposit weight `h(w)` | stored as | the row it contributes | separates the components on |
|---|---|---|---|
| `1` | `count` | `(E_g[A], E_r[A])` | mean opportunity |
| `1/L` | `inv_length_sum` | `(E_g[A/w], E_r[A/w])` | **mean length** |
| `1/A` | ⛔ *not stored* | `(P_g(A>0), P_r(A>0))` ≈ `(1,1)` | **nothing — deliberately flat** |
| `L` | `length_sum` | `(E_g[Aw], E_r[Aw])` | **spread (2nd moment)** |

Pick two non-parallel rows, solve the 2×2. Done.

### 1.3 Worked example — a 500 bp node

gDNA fragments all 50 bp, RNA all 200 bp, both at `rho = 10` fragments per base.

```
opportunity     A_g = 500-50+1 = 451        A_r = 500-200+1 = 301
true counts     N_g = 4510                  N_r = 3010          N = 7520

observables     Sum 1/L = 4510/50 + 3010/200 = 105.25
                Sum 1/A = 4510/451 + 3010/301 = 20.000  == rho_g + rho_r, EXACTLY, with no model

invert (count, Sum 1/L):  N_g + N_r = 7520 ;  N_g/50 + N_r/200 = 105.25   ->  N_g = 4510  ✓
invert (count, Sum 1/A):  451·rho_g + 301·rho_r = 7520 ; rho_g + rho_r = 20  ->  rho_g = 10 ✓
```

Both work. ⚠ Note `Sum 1/L / 500 = 0.21`, **not** 10 — `Sum 1/L` has units of *fragments per base of
fragment*, not per base of node. The quantity that is a genuine node density is `Sum 1/A`
(`docs/accumulator/NODE_DENSITY_DERIVATION.md` T1: `1/A` is the **unique** weight with that property).

### 1.4 Why one channel is not enough — the singular case

Take gDNA `150 ± 20` and RNA `150 ± 80`, means matched to machine precision:

| frame | `det(count, Sum1/L)` | `cond(count, Sum1/L)` | `cond(+ SumL)` |
|---|---|---|---|
| **contiguous edge** | **2.0e-14** | **3.3e+18** ⛔ singular | **1.5e+03** |
| node ell=151 | 3.3 | 371 | 377 |
| node ell=1000 | 2.7e+03 | 546 | 6052 |

At an edge the `count` row is `(mu_g - 1, mu_r - 1)` and the `Sum1/L` row is `(1, 1)`, so the
determinant is exactly `mu_g - mu_r`. **Equal mean lengths ⇒ the two channels are the same equation, at
any depth.** `SumL` reads the second moment instead (row ratio 1.263 where `Sum1/L`'s is 1.000000) and
rescues it.

⚠ Equal means is not a corner. The four real cfRNA libraries: gDNA `146 / 102 / 157 / 309`, RNA
`227 / 270 / 206 / 246` — vcap already has gDNA **longer** than RNA, so the difference passes through
zero somewhere inside the space Rigel ships into. Under capture, on-target gDNA runs ~42 bp shorter than
off-target, so the separation is compartment-specific *within* one library.

### 1.5 Three jobs, three numbers

| channel | job | droppable? |
|---|---|---|
| `count` | statistical power. `Var(log rho_c) = 1/(f_c·n)`; a Beta-Binomial needs an integer | **no** |
| one **level/mean** channel (`Sum1/L` today, `Sum1/A` under P3) | the split's first discriminant | **no** |
| one **spread** channel (`SumL`) | what survives equal means | only if that failure is accepted |

---

## §2 THE TWO DEFECTS

### 2.1 Defect A — the solver never reads the length equations

**Verify it yourself before believing it:**

```bash
grep -rn "inv_length_sum\|length_sum" src/ | grep -v "scan_payload\|native/\|substrate.py"
# -> no hits. Both channels are decoded in calibration/substrate.py and read by NOBODY.
```

`build_node_geometry` ([node_geometry.py:186-193](../src/rigel/calibration/node_geometry.py#L186-L193))
pulls **only** `.count` from the substrate. The fitted length models enter solely as *divisors*
(`eff_gdna` / `eff_rna`) converting a count into a density.

So the solver's second equation comes from [node_init.py](../src/rigel/calibration/node_init.py)'s four
sources: (1) structural lock, (2) a population gDNA-density prior, (3) the strand Beta-Binomial,
(4) nothing — "100 % gDNA at zero precision". **Fragment length is not on that list.**

⭐ **Why that is the hole.** The strand channel carries information `I(f_g) ∝ (2κ-1)²`, which is
**exactly zero at κ = ½** — the gDNA fraction cancels identically out of the Beta-Binomial mean
(`docs/SESSION_HANDOFF.md` §2, verified to 5.6e-17). That is why
[node_init.py:250-251](../src/rigel/calibration/node_init.py#L250-L251) credits both-strand (AMBIG)
nodes **zero** strand evidence. On unstranded data, or at any both-strand locus, sources 1–3 are all
silent and the solver falls through to source 4.

**The evidence is in the project's own baseline** (`docs/WIP.md`, S5.f, truth `f_gdna = 0` exactly on the
`none` arm):

| condition | `f_gdna` | |
|---|---|---|
| none **ss0.99** capture_off | **0.0030** | strand alive |
| none **ss0.99** capture_on | **0.0016** | strand alive |
| none **ss0.50** capture_off | **0.0793** | ⛔ **strand dead — 26× worse** |
| gdna100 **ss0.50** capture_on | **0.3754** | ⛔ 25 % low; the ledger records it as "unexplained" |

**The two worst rows in the entire baseline are the two where the strand channel is dead.** 443,277
phantom gDNA fragments in a library that contains none.

### 2.2 Defect B — the output sums an intensive quantity

[priors.py:252-292](../src/rigel/calibration/priors.py#L252-L292) turns per-object answers into the two
per-locus EM pseudocounts by **summing counts across objects**:

```python
gdna_prior_count = SUM_objects  f_g · count          # raw sum
rna_prior_count  = SUM_objects  (1-f_g) · count      # raw sum, and NO divisor at all
```

Those go straight into the EM as **additive pseudocounts in fragment units**, competing with the locus's
real fragment counts (`G = n_gdna + a_g`, `R = n_rna + a_r`,
[em_solver.cpp:735-741](../src/rigel/native/em_solver.cpp#L735-L741)).

But a fragment deposits on **`max(K, 1)` objects**, where `K` is the number of lines it crosses
(`_accumulator_reference.deposit`: contained iff `K = 0`, otherwise one deposit per crossed line). For a
partition of spacing `s`:

```
    incidences(w) = max( 1 , (w-1)/s )
```

⭐ **Counts are conserved exactly where every node is longer than every fragment, and become a
length-weighted count where they are not.** 56.7 % of human nodes are shorter than one 200 bp fragment.

The bias in the ratio is exactly `SUM A_g / SUM A_r`. Measured on the chr22 pilot scan caches:

| | gDNA | RNA |
|---|---|---|
| incidences per fragment (mass-weighted, measured) | **1.031** | **≈1.17–1.24** |
| share of incidences landing on edges | **8.1 %** | **38.4 %** |

→ **the prior's g : r ratio under-calls gDNA by ~13–19 %**, and the error grows as the partition gets
finer. This is `docs/SESSION_HANDOFF.md` §3 trap 4 in the mirror: the old design's fractional mass depended on
region sizes; the new design's incidence *multiplicity* does.

⚠ **The baseline's own headline is that sum.** `5,370,056 + 5,373,548 = 10,743,604` equals
`node_contained_count + edge_unspliced_count + edge_spliced_count` in the payload exactly. The reported
`f_gdna` is incidence-weighted, not fragment-weighted.

**Reproduce both numbers:**

```bash
python - <<'EOF'
import numpy as np, os
base = os.path.expanduser("~/Downloads/rigel_runs/suite/pilot/scan_cache")
def tot(d):
    z = np.load(os.path.join(base, d, "payload.npz"))
    return (int(z['node_start_count'].sum()),
            int(z['node_contained_count'].sum()) + int(z['edge_unspliced_count'].sum())
            + int(z['edge_spliced_count'].sum()))
fA, iA = tot("gdna_gdna100_ss_0.50_nrna_none_capture_off")   # mixture
fB, iB = tot("gdna_none_ss_0.50_nrna_none_capture_off")      # zero gDNA -> pure RNA
print("RNA  incidences/fragment =", iB/fB)
print("gDNA incidences/fragment =", (iA-iB)/(fA-fB))   # differencing isolates the gDNA component
EOF
```

---

## §3 THE INVARIANT YOU ALREADY HAVE

`node_start_count` (`uint32[n_nodes]`) is incremented **once per accepted fragment**, at the node
containing its first covered base. `SUM node_start_count == qc.deposited`, exactly, externally
checkable. It is the one axis on which every fragment appears exactly once.

⚠ It is currently used only as an assertion in two unit tests. **P1's T3 makes it a real gate.**

---

## §4 P0 — MEASURE THE HOLE ✅ **DONE, 2026-07-31 — see `docs/WIP.md`'s P0 entry**

    Tool: scripts/design/composition_evidence_census.py   (no production code was changed)

⭐ **The premise survived falsification, and the finding is larger than the prediction.** Headline, mass-
weighted share of unspliced mass whose own composition precision is zero (structural locks excluded):

| condition | κ | no-evidence | `no-ev \| AMBIG` | `no-ev \| single-strand` |
|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 0.5000 | **49.4 %** | 93.3 % | **63.6 %** |
| gdna100 ss0.50 capture_on | 0.4990 | **98.2 %** | 99.9 % | **98.9 %** |
| gdna100 ss0.99 capture_off | 0.0101 | **12.7 %** | 93.3 % | **0.0 %** |
| gdna100 ss0.99 capture_on | 0.0101 | **30.0 %** | 99.9 % | **0.0 %** |
| **none** × 4 | both | **100.0 %** | 100 % | 100 % |

⭐ **Three findings that were not predicted and that change how P2 is scored:**

1. **AMBIG mass is blind in every condition, stranded or not** (93.3–100 %). The Schur gate silences
   strand on both-strand nodes by design, so **13.3–40.1 % of library mass has never had own composition
   evidence in any library Rigel has run.** This is not an unstranded-data problem. ⭐ It is the strongest
   argument for P2: the length likelihood has no dependence on the tilt `θ`, so the Schur argument that
   kills strand there does not apply to it.
2. **A zero-gDNA library is 100 % blind at any κ** — `strand_evidence`'s noise floor divides by
   `N_gdna`, so `N_gdna = 0 ⇒ disc = 0`. Yet `none ss0.50 capture_off` reports `f_gdna = 0.0793`.
3. **Capture roughly doubles blindness by removing the anchor**: structurally-locked mass collapses
   **28.7 % → 1.0 %**.

⛔ **Consequence for P2's gate (§6.4): do not score on the `none` arm alone.** All four zero-gDNA
conditions are saturated at 100 % blind, so anything scores "better" there — trap 19. The clean
comparison is the **gdna100** arm, and `gdna100 ss0.50 capture_on` (98.2 % blind, `f_gdna` 0.3754 against
~0.50, the row S5.f recorded as "unexplained") is the sharpest single target in the suite.

**The falsification.** `--inject-kappa 0.5` on `gdna100 ss0.99 capture_off` — one field replaced on the
condition's own fitted priors, nothing else — moves it 12.7 % → **49.4 %** and 0.0 % → **63.6 %**,
reproducing the independently-simulated `ss0.50` condition **to the last digit**. The census tracks the
strand channel and nothing else.

```bash
python scripts/design/composition_evidence_census.py \
    --index ~/Downloads/rigel_runs/suite/rigel_index \
    --cache-root ~/Downloads/rigel_runs/suite/pilot/scan_cache [--inject-kappa 0.5]
```

<details><summary>The original P0 specification, kept for provenance</summary>

**Purpose:** falsify the *premise* of P2 before building it. If the mass at `tau_lam == 0` is small, the
diagnosis in §2.1 is wrong and P2 must not be built.

**How.** Write `scripts/design/composition_evidence_census.py`. For each of the 8 pilot scan caches, run
`calibrate()` with the `_debug` hook (`calibrate.py` already populates it with `chain`, `geometry`,
`statics`, `belief`), recompute `build_node_init`'s `tau_lam`, and report — **weighted by unspliced
mass, never by node count**:

| output | why |
|---|---|
| share of unspliced mass on slots with `tau_lam == 0` | mass solved by prior + messages alone |
| the same, split by `kappa = 0.5` vs `0.99` | isolates the strand channel's contribution |
| the same, split by single-strand vs AMBIG (`free_pos ^ free_neg`) | isolates the Schur gate |
| the same, split by NODE vs EDGE slot | tells you which axis P2 must serve first |

⚠ **Mass-weight it.** `docs/SESSION_HANDOFF.md` §1 fact 6 records a claim that survived for months because it
was bp-weighted (0.8738) when the estimator is mass-weighted (0.9596). Do not repeat that.

**Recorded prediction, to be written down before the run:** the `ss0.50` conditions show a materially
larger `tau_lam == 0` mass share than the `ss0.99` conditions, concentrated on AMBIG slots.

**Gate.** None — this is a measurement. But record it in `docs/WIP.md`; it is P2's before-picture.

</details>

---

## §5 P1 — THE UNITS FIX

### 5.1 What changes, mathematically

Replace "sum the counts" with "pool the density, then integrate it over the span":

```
    rho_g(locus) = SUM_locus share·m_g  /  SUM_locus share·A_g      <- ratio of SUMS
    rho_r(locus) = SUM_locus share·m_r  /  SUM_locus share·A_r
    gdna_prior_count = rho_g · span_bp
    rna_prior_count  = rho_r · span_bp
```

where

* `m_g`, `m_r` are the per-object deconvolved masses already in `CalibrationResult`;
* `A_g` is `gdna_node_eff_len` / `gdna_edge_eff_len` — **already carried**;
* `A_r` is the RNA twin — **computed in `build_node_geometry` as `eff_rna` and currently discarded**;
* `span_bp` is the locus's **genomic** span, `SUM_locus share · region_size_bp`. ⚠ The *same* span for
  both components — `rho_c` is defined as fragments of component `c` starting per genomic base, so the
  span cancels out of the ratio and only sets the scale.

⚠ **Ratio of sums, never mean of ratios** (`docs/SESSION_HANDOFF.md` §2, `rho_bg = Sum g / Sum E`).

⚠ Edges contribute mass and support but **zero** `span_bp` — an edge is a 0-bp line. That is correct and
not an omission.

### 5.2 Files, in dependency order

**(1) [`src/rigel/calibration/calibrate.py`](../src/rigel/calibration/calibrate.py), line 135.**

Rename `_project_eff_gdna(chain, geometry, payload)` → `_project_eff(chain, eff_slots, payload)`, taking
the per-slot array as an argument instead of reaching into `geometry.eff_gdna`. Call it twice:

```python
node_eff_gdna, edge_eff_gdna = _project_eff(chain, geometry.eff_gdna, payload)
node_eff_rna,  edge_eff_rna  = _project_eff(chain, geometry.eff_rna,  payload)
```

⚠ Keep the existing docstring's point: **this is a projection, not a recomputation.** Do not call
`contained_eff_length` / `crossing_eff_length` again here — two implementations of one quantity is
`docs/SESSION_HANDOFF.md` §3 traps 2 and 27.

Pass both new arrays into the `CalibrationResult(...)` constructor at the end of `calibrate`.

**(2) [`src/rigel/calibration/result.py`](../src/rigel/calibration/result.py).**

Add two fields beside their gDNA twins:

```python
rna_node_eff_len: np.ndarray   # float64[n_nodes] -- contained_eff_length on the RNA pmf
rna_edge_eff_len: np.ndarray   # float64[n_edges] -- crossing_eff_length  on the RNA pmf
```

Add `"rna_node_eff_len"` to the `n_nodes` loop and `"rna_edge_eff_len"` to the `n_edges` loop in
`__post_init__`. ⚠ The shape check is the load-bearing one: `E = N - n_refs` differs from `N` by only a
few hundred genome-wide, so a mis-keyed array is a *plausible* length and nothing downstream would fault.

**(3) [`src/rigel/calibration/priors.py`](../src/rigel/calibration/priors.py).**

Generalize `_gdna_node_arrays` to take the component's mass and support arrays:

```python
def _component_node_arrays(mass_node, mass_edge, eff_node, eff_edge, calibration, region_arrays):
    """Per-node (mass_total, support_total, edge_mass, edge_support) for ONE component.
    Identical to the former _gdna_node_arrays, with the arrays passed in rather than read off
    the gDNA fields. `edge_owner_nodes` is unchanged."""
```

Call it twice in `assemble_priors`. Then project **five** arrays instead of the current set:

| projected name | value |
|---|---|
| `gdna_mass` | `mass_g_node + mass_g_edge` (edges via `edge_owner_nodes`) |
| `gdna_support` | `A_g_node + A_g_edge`, same attribution |
| `rna_mass` | `mass_r_node + (mass_rna_edge - mass_rna_spliced_edge)`, same attribution |
| `rna_support` | `A_r_node + A_r_edge`, same attribution |
| `span_bp` | `region_arrays.region_size_bp` (nodes only) |

and finish:

```python
rho_g = np.divide(proj["gdna_mass"], proj["gdna_support"],
                  out=np.zeros_like(proj["gdna_mass"]), where=proj["gdna_support"] > 0.0)
rho_r = np.divide(proj["rna_mass"], proj["rna_support"],
                  out=np.zeros_like(proj["rna_mass"]),  where=proj["rna_support"]  > 0.0)
gdna_prior_count = rho_g * proj["span_bp"]
rna_prior_count  = rho_r * proj["span_bp"]
```

⛔ **Never floor a zero support to epsilon.** `docs/SESSION_HANDOFF.md` §3 trap 23: an object with no
opportunity for a component must emit **nothing at zero precision**, never a floored division. That
default once seeded false gDNA into neighbouring exons. Use the `where=` form above.

⚠ `gdna_eff_len` (the third `LocusPriors` field, the IPR contraction) is **unchanged by P1**. It is
already a density-based quantity and its `span` is the *effective* support, deliberately. Do not
"consistency-fix" it into the genomic span.

⚠ `mass_rna_spliced_edge` stays withheld from `rna_prior_count`, and `mass_rna_junction` stays unread.
Both are settled owner rulings (`result.py`'s docstring); P1 does not reopen them.

**(4) Fixtures that hand-build a `CalibrationResult`** — all will fail on the new kwargs and must be
updated: `tests/calibration/test_result_schema.py`, `test_priors.py`, `test_capture_eff_length.py`,
`tests/calibration/_oracle.py`.

### 5.3 Gates — write these FIRST and verify them failing

⛔ **`docs/SESSION_HANDOFF.md` §3 trap 1 and the owner's standing rule: a test written after the code proves
nothing, and a test that is never seen to fail proves less.** Run each one against unmodified `HEAD` and
record the failure before writing any implementation.

| id | test | asserts | current behaviour |
|---|---|---|---|
| **T1** | **partition invariance**. Build two indexes from the same GTF, the second with extra cuts inserted *inside* exons (physically nothing changes). Scan the same BAM against both, calibrate, compare `gdna_prior_count` per locus | agreement within Poisson noise | ⛔ **moves** |
| **T2** | **length sweep**. Fixed true g:r; sweep `mu_g/mu_r` over `{0.5, 0.75, 1.0, 1.5, 2.0, 3.0}` — **both directions**, per the owner's ruling that there is no rule RNA is longer | the prior's `g:r` is flat in `mu_g/mu_r` | ⛔ drifts by `SUM A_g / SUM A_r` |
| **T3** | **conservation**. `SUM_loci (gdna_prior_count + rna_prior_count)` against `SUM node_start_count` over those loci's nodes, minus the spliced and junction populations | agreement to the spliced/junction accounting | never checked |
| **T4** | **zero support**. A locus whose nodes are all shorter than one fragment (`A_r == 0`) | `rna_prior_count == 0`, not `inf`, not `nan` | — |

**Perturbations that must break them** (`falsification_needs_perturbation`: breaking the code and
watching the test fail is the other half of the discipline):

| | perturbation | must fail |
|---|---|---|
| P1a | pass `eff_gdna` where `eff_rna` belongs | T2 |
| P1b | drop the `* span_bp` | T3 |
| P1c | revert to the raw sums | T1, T2 |
| P1d | use `mean of ratios` instead of `ratio of sums` | T1 |
| P1e | floor the support to `1e-9` instead of `where=` | T4 |

### 5.4 ⚠ THE SCORING TRAP — read before running the A/B

On `gdna100 ss0.50 capture_off` the two defects currently **partly cancel**:

```
fragment truth                        f_gdna ~ 0.519
incidence-weighted truth              f_gdna ~ 0.480      (the 1.031 vs ~1.17 tilt)
S5.f reports                          f_gdna   0.4998     <- the solver over-calls by ~4% on top
```

⚠ The 0.519 is **derived, not read off the manifest**: `node_start_count` totals 9,634,502 on this
condition against 4,636,413 on its zero-gDNA twin, so 4,998,089 accepted fragments are gDNA. The
manifest's `gdna_rate = 1.0` implies "≈ 0.5" *before* the `max_fragment_length` filter, which drops
280,100 fragments non-uniformly. ⛔ **Re-derive it from the payloads for the arm you actually score** —
these are approximate and the sign of the P1 move is the load-bearing part, not the third digit.

**Fix the aggregation alone and the sign of the error flips.** Expect the headline to move to roughly
0.53 against a truth of 0.519 — i.e. *further from 0.50* while *closer to truth*.

⛔ **Score P1 against the simulator's own fragment counts, never against the S5.f table.** This is
trap 2's family (cancelling errors hiding a scale error) and trap 19's (a one-sided metric). Write the
prediction into `docs/WIP.md` **before** the run so it cannot be rationalised afterwards.

⚠ Run the end-to-end arm under `EMConfig.assignment_mode = "map"` or `"fractional"`, held fixed across
both arms — calibration is bit-identical run to run but the EM samples from the posterior by design
(`docs/calibration/S5_DESIGN_LOG.md` §0).

---

## §6 P2 — THE LENGTH LIKELIHOOD (this is S5.a2)

### 6.1 The maths, stated completely

Fix one object `o` and one population `p`. Define, for fragment length `w`:

```
    A(w)   the opportunity        (ell - w + 1)+  contained at a node
                                  (w - 1)+        crossing at an edge
    u(w)   the stored inv weight  1/w             at a NODE     (matches inv_length_quantum(length))
                                  1/(w-1)         at an EDGE    (matches inv_length_quantum(length-1))
```

⚠ Both `inv_length_sum` weights come straight from `_accumulator_reference.deposit`; do not re-derive
them from prose. `length_sum` deposits the full `w` in **both** frames.

A component-`c` fragment that **landed here** has the opportunity-tilted pmf

```
    g_c(w) = f_c(w)·A(w) / E_c[A]           and   E_c[A]  IS  eff_gdna / eff_rna at this slot
```

Precompute five scalars per object per component (all pure functions of the fitted pmf and the object's
geometry):

```
    m1_c  = SUM_w g_c(w)·u(w)          matches inv_length_sum
    m2_c  = SUM_w g_c(w)·w             matches length_sum
    q1_c  = SUM_w g_c(w)·u(w)^2
    q2_c  = SUM_w g_c(w)·w^2
    q12_c = SUM_w g_c(w)·u(w)·w
```

Let `pi` be the gDNA share of this object's landed count. ⭐ **`pi` IS `f_g`** — `chain_node_deconv`
computes `gdna_mass = f_g · count`, so the λ grid already parameterizes exactly this quantity. Then,
conditional on the observed integer `N` (the count summed over both strand columns):

```
    mean:   mu(pi)   = N · ( pi·m1_g + (1-pi)·m1_r ,  pi·m2_g + (1-pi)·m2_r )

    2nd:    E_pi[u^2]  = pi·q1_g  + (1-pi)·q1_r
            E_pi[w^2]  = pi·q2_g  + (1-pi)·q2_r
            E_pi[uw]   = pi·q12_g + (1-pi)·q12_r

    cov:    Sigma(pi) = N · [[ E_pi[u^2] - E_pi[u]^2 ,  E_pi[uw] - E_pi[u]E_pi[w] ],
                             [      (sym)           ,  E_pi[w^2] - E_pi[w]^2      ]]

    loglik: l(pi) = -0.5·(x - mu)^T Sigma^-1 (x - mu)  -  0.5·log det Sigma
                    with x = ( inv_length_sum , length_sum ), both summed over the strand axis
```

⚠ **Why the Gaussian is legitimate here.** `x` is a **sum of `N` i.i.d. draws**, not a ratio. The heavy
tail recorded in `docs/calibration/S5_DESIGN_LOG.md` §3.4 (realised sd 187.5 against a predicted 0.375) is the *ratio
estimator* `phi-hat`'s tail, and the robust scale there matched prediction to ~10 %. This is the same
quasi-likelihood `scripts/design/observable_efficiency.py` uses to produce **every** efficiency number
in the design log — it is being lifted into the solver, not invented.

⚠ **Strand-agnostic.** The gDNA/RNA question does not depend on which genome strand the read aligned to,
so sum both columns. The strand Beta-Binomial keeps its own separate columns and is untouched.

### 6.2 ⭐ The integration point already exists

`simplex_logodds` accepts `lam_logprior`, an `(m, K)` array on the σ(λ) grid, and adds it straight into
`psi`:

* 1-D solve: [simplex_logodds.py:295-296](../src/rigel/calibration/simplex_logodds.py#L295-L296)
* 2-D (AMBIG) solve: [simplex_logodds.py:506-507](../src/rigel/calibration/simplex_logodds.py#L506-L507)

**The length likelihood has exactly that shape.** Add a sibling parameter — do **not** fold it into
`lam_logprior`, because a separate argument is what makes the A/B vary one thing (the same pattern as
`edge_rna_reach`, the A7 switch: `None` ⇒ byte-identical to today).

⚠ Thread it through the single-strand regrid path too — see `_regrid_global(_s(lam_logprior, bidx), ...)`
at [simplex_logodds.py:684-689](../src/rigel/calibration/simplex_logodds.py#L684-L689). Missing this
leaves the term silently absent on exactly the single-strand nodes, which is the half of the genome
where it is easiest to validate.

⭐ **And its precision needs no new code either.** `density_factor_precision(lam_logprior, lam_grid)`
([density_deconv.py:186](../src/rigel/calibration/density_deconv.py#L186)) reads the curvature of *any*
`(m, K)` λ-factor and returns `tau`. It is named for its first caller but the maths is generic. So
[node_init.py:251](../src/rigel/calibration/node_init.py#L251) becomes:

```python
tau_lam = np.where(single_strand, i_strand, 0.0)
tau_fac = density_factor_precision(intron_prior, lam_grid)
if tau_fac is not None:
    tau_lam = tau_lam + tau_fac
tau_len = density_factor_precision(length_loglik, lam_grid)      # <- source 4, finally live
if tau_len is not None:
    tau_lam = tau_lam + tau_len
```

⚠ Unlike `i_strand`, `tau_len` is **not** gated to single-strand nodes. The Schur argument that kills the
strand term on AMBIG nodes is that the strand likelihood is rank-1 in the tilt `theta`; the length
likelihood does not depend on `theta` at all, so it informs `lambda` directly on every node. **That is
the entire point of the step** — it is the only source that speaks on an AMBIG node.

### 6.3 Files

| file | change |
|---|---|
| **new** `src/rigel/calibration/length_likelihood.py` | the tilted moments `m1/m2/q1/q2/q12` per slot per component, and the `(n_slots, K)` log-likelihood on the λ grid. Pure; imports only `numpy` and `effective_length` |
| [`node_geometry.py`](../src/rigel/calibration/node_geometry.py) §`build_node_geometry` | gather `inv_length_sum` and `length_sum` onto the chain beside `unspliced_count`, same three-line pattern (`node_contained` at NODE slots, `edge_unspliced` at EDGE slots). Add two `NodeGeometry` fields |
| [`simplex_logodds.py`](../src/rigel/calibration/simplex_logodds.py) | new `length_loglik=None` parameter, added to `psi` at the two sites above and threaded through `_solve_nodes_logodds_all` + `_regrid_global` |
| [`node_init.py`](../src/rigel/calibration/node_init.py) | build `length_loglik`, pass it to `_solve_nodes_logodds_all`, add `tau_len` to `tau_lam`; update the module docstring's "four sources" list — source 4 is no longer "nothing" |
| [`bp_solver.py`](../src/rigel/calibration/bp_solver.py) §`node_sweep` | pass the array through to the per-node solve. **No other change**; the message model is untouched |
| [`config.py`](../src/rigel/config.py) | one bool on `CalibrationConfig`, default **off** for the first landing so the A/B has an arm |

### 6.4 Gates — write these FIRST and verify them failing

**Unit, against brute force, no tolerance:**

| id | test |
|---|---|
| **U1** | `m1_c`, `m2_c`, `q1_c`, `q2_c`, `q12_c` against **exact enumeration** over integer start positions on a small node (`ell` 5, 10, 25, 151) and a line, for a point-mass pmf and a two-point pmf. ⛔ Enumerate; do not call the module's own helper (`docs/SESSION_HANDOFF.md` §3 trap 1) |
| **U2** | `E_c[A]` computed inside `length_likelihood` equals `eff_gdna`/`eff_rna` from `build_node_geometry`, to the last bit. If it does not, there are two implementations of the effective length |
| **U3** | `l(pi)` is maximised at the true `pi` on synthetic data drawn from a known mixture (1000 draws, `mu_g != mu_r`) |
| **U4** | with `f_g == f_r` (identical pmfs) the log-likelihood is **flat in `pi`** and `tau_len == 0`. A flat factor must carry no information, never the grid's own width |
| **U5** | `N == 0` ⇒ the row is flat and `tau_len == 0`, no `nan`, no division by zero |

**⭐ The end-to-end gate, and its data is already on disk:**

```
gdna_none_ss_0.50_nrna_none_capture_off :  f_gdna  0.0793  ->  must fall toward 0
gdna_gdna100_ss_0.50_capture_off        :  f_gdna  0.4998  ->  must not degrade
gdna_gdna100_ss_0.50_capture_on         :  f_gdna  0.3754  ->  the second target
```

Those conditions are unstranded, so the strand channel is **provably** zero and the length channel is the
only thing that can move them.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
OMP_NUM_THREADS=1 python scripts/design/build_scan_cache.py     # already cached; --force to rebuild
# then calibrate() per condition, both arms, 2.6-6.4 s each
```

⛔ **Score all 8 conditions, both arms.** `docs/SESSION_HANDOFF.md` §3 trap 19: on a zero-gDNA library *any*
change that lowers the gDNA fraction scores better, and that one-sidedness has reversed a published
verdict in this project once (reported −13.1 % win; the full battery made it worse). Score per node and
per condition on **soft** quantities, never on a pooled hard label.

**Perturbations:**

| | perturbation | must fail |
|---|---|---|
| P2a | set `m1_g = m1_r` and `m2_g = m2_r` | the entire end-to-end gain must vanish |
| P2b | use the **untilted** pmf `f_c` instead of `g_c` | a measurable bias on U1/U3 |
| P2c | swap `length_sum` into the `inv_length_sum` slot | loudly, at U1 |
| P2d | drop the `-0.5 log det Sigma` term | U3 (the composition-dependent covariance is what makes it a likelihood) |
| P2e | gate `tau_len` to single-strand nodes (i.e. copy `i_strand`'s gate) | the end-to-end gate — AMBIG nodes are the whole point |
| P2f | omit the `_regrid_global` threading | a single-strand-only stratum must show the term absent |

⚠ Include a **low-count stratum** in the scoring (objects with `N <= 3`): 94 % of confident-gDNA nodes on
LBX0190 carry zero counts and 80.5 % of partition nodes carry none. The Gaussian is asymptotic; confirm
it does not turn confident on thin data. If it does, the honest fix is that `own_precision(n, v, live)`
already damps by the count — not a threshold constant (⛔ **no magic numbers**).

---

## §6.5 ⛔ P2 IS BLOCKED — THE ROAD OUT

    Measured 2026-07-31; `docs/WIP.md`'s P2 entry has the numbers and its own correction note.

⭐ **The likelihood is not the defect.** Two owner points, both verified:

* **EB shrinkage to the global histogram is correct**, and at ``n_gdna = 0`` it shrinks all the way —
  ``gdna_pmf == global_pmf`` byte-identically.
* **Identical pmfs for the two components make this channel exactly inert** — proven byte-identically
  end to end. The channel is a *comparator*; give it one ruler and it says nothing.

⛔ **So the defect is that the two rulers differ.** On a zero-gDNA library every fragment is RNA, so
``global`` and ``rna_fl_pmf`` describe **one population** — and they disagree:

| | global (the EB anchor) | RNA_SPLICED pool | gap |
|---|---|---|---|
| mean | 210.1 | 234.5 | **+11.6 %** |
| sd | 85.5 | 146.2 | **+71.1 %** |
| support | [50, **713**] | [50, **1000**] | ⛔ |
| mass > 500 bp | 0.10 % | **5.53 %** | **53×** |

Two causes, and they are separable:

1. **The junction-opportunity tilt** (`docs/accumulator/DESIGN.md` §8.1(b)). A fragment enters `RNA_SPLICED`
   by crossing ≥1 annotated junction, which longer fragments do more often. The measured +11.6 % / +71 %
   sits against `docs/calibration/S5_DESIGN_LOG.md` §3.6's independently-predicted **+14 % / +50 %**, and the
   ``rna/global`` ratio rises with length (log-ratio vs length, corr **+0.70**) — the opportunity
   signature. ⚠ §8.1(b) has been **"not yet decided"** since S5.b. **It is now the blocker.**
2. ⭐ **A FRAME mismatch, which is new and is not a tilt.** ``global`` is the **scanner's** histogram;
   the five pools are the **accumulator's** ``pool_lengths``. Different code paths, and the support
   ceilings differ (713 vs 1000). A 53× excess of >500 bp molecules is not a smooth opportunity effect.
   ⚠ `scan_cache.calibration_inputs`'s own docstring already half-names this — "Only the unconditional
   ``global`` histogram still comes from the scan, because no pool is unconditional" — and **EB-shrinking
   accumulator-frame pools toward a scanner-frame anchor is `docs/SESSION_HANDOFF.md` §3 trap 27**: two
   implementations of one quantity, disagreeing.

⭐ **The full audit of this area is `docs/accumulator/FRAGMENT_LENGTH_AUDIT.md`** — F1/F2 below are C1/C3
there, and the audit found four more defects (a dead cache field, the scanner's silent ambiguous-length
drop, and the shrunk pmf reaching the EM's transcript effective lengths).

### The steps, in dependency order

| | step | why it is where it is |
|---|---|---|
| **F1** | **Give the EB anchor the accumulator's own frame.** Either an unconditional length bin in the accumulator (⚠ reopens the S3 byte-identity gate) or an anchor built from the pools themselves. **Owner call** — the first is correct, the second is cheap | until the anchor and the pools measure the same thing, no downstream length comparison means anything, and F2 cannot be measured |
| **F2** | **`docs/accumulator/DESIGN.md` §8.1(b): divide each pooled histogram by its own opportunity before normalising** — ``placements(w)`` for a crossing pool, ``(ell−w+1)+`` for a contained pool, the transcript-level count for the junction pool | the named blocker. ⚠ The junction pool's opportunity is the hard one and §8.1(b) says so |
| **F3** | **Gate the length channel on BOTH pools having data**, the `strand_evidence` analogue (its noise floor divides by ``n_gdna_obs``, so ``N_gdna = 0 ⇒ disc = 0``). Even with F1+F2 correct, a zero-gDNA library has no gDNA length model and the discriminant is **undefined**, not small | independent of F1/F2 and worth doing regardless; it is the difference between "no evidence" and "confident nonsense" |
| **F4** | **Re-run `scripts/design/length_likelihood_ab.py`.** Before/after pictures are both recorded | — |

### ⭐ The falsification test F2 never had

> **On a zero-gDNA library every fragment is RNA. So after F1+F2, the de-tilted ``rna_fl_pmf`` must equal
> the unconditional histogram — and the length likelihood must go EXACTLY inert.**

Pass/fail, no tuning, on four conditions that already exist on disk. §8.1(b) has been an open item for
lack of a way to judge it; this is one.

⛔ **Do not damp the channel to make the numbers look better.** It is reporting a real disagreement
between two length models. Damping hides an upstream defect behind a tuned constant —
`docs/SESSION_HANDOFF.md` §3 trap 12, recorded three times over.

---

## §7 P3 — THE DEPOSIT WEIGHT `L -> A` (decide after P2)

⛔ **Do not start this until P2 is measured.** It reopens the S3 byte-identity gate against
`tests/native/_accumulator_reference.py` and it may turn out to be unnecessary.

**What.** At a node, deposit `round(2^32 / A)` where `A = (ell - L + 1)` for contained and
`(L - ell - 1)` for spanning, instead of `round(2^32 / L)`. The edge is untouched — `1/(L-1)` already
**is** the reciprocal opportunity. `density_quantum` already takes a placement count and the node length
is in the cut array the deposit already binary-searches. Overflow is priced and closed
(`docs/accumulator/NODE_DENSITY_DERIVATION.md` §7.3: even 10^8 fragments at `A = 1` reaches 2.3 % of uint64).

**Why it might be worth it, and this argument is NOT in the existing ruling.** A4 was ruled R-b (do not
store `Sum1/A`) on the grounds that it buys almost no **composition** information — median efficiency
0.953 → 0.960. That pricing scored the *split*. Two consumers of the **level** have since appeared:

1. **P1 needs `rho_c` per locus**, and today `rho_c = count / E_f[A]` depends on the fitted length models
   in exactly the frame where they are weakest (`docs/accumulator/DESIGN.md` §8.1(a): at a contained node the
   naive ratio converges to the **harmonic** mean, biased low — and the only pure gDNA pool *is*
   intergenic-contained).
2. ⭐ **`node_total_density` is the message-reframe currency for the entire belief propagation**
   ([bp_solver.py:534](../src/rigel/calibration/bp_solver.py#L534): `r = rho_tot(dst)/rho_tot(src)`).
   Today it is `node_total_density(geometry, f_g)` — a function of **the belief**, which is why the frame
   had to be frozen at the init belief with a 20-line comment explaining that a second iteration was a BP
   violation. `Sum1/A` makes `rho_tot` a **measured constant**, so `docs/SESSION_HANDOFF.md` §3 trap 11 ("a
   message may use the destination's CONSTANTS, never its BELIEFS") holds by construction rather than by
   a freeze.

**Before ruling.** Re-run `scripts/design/observable_efficiency.py` with the **level** scored as well as
the split. The A4 ruling was made against half the objective.

**The open sub-question.** If `Sum1/A` is stored, is `Sum1/L` still worth its bytes? `(count, Sum1/A)`
already identifies the split on the mean (determinant `A_g - A_r`), so `Sum1/L` is a better-conditioned
version of the *same* discriminant at short nodes, not an independent one. Measured cost of dropping it:
`0.960 -> 0.918` at 151 bp but `0.782 -> 0.749` at 1000 bp and `0.832 -> 0.804` at 3000 bp
(`docs/accumulator/NODE_DENSITY_DERIVATION.md` §6). ⛔ `SumL` is **not** droppable — it is the only channel that survives
equal means.

---

## §8 WHAT THIS PLAN DELIBERATELY DOES NOT DO

| | why |
|---|---|
| touch the deposit rule's **partition** (no `1/K` split) | measured: `1/K` reads up to **3.6×** low, and agrees with the correct rule only where nodes are coarser than fragments (`docs/accumulator/DESIGN.md` §4.2). The per-object estimator is correct and partition-free; only the *aggregation* was wrong |
| touch the strand Beta-Binomial | it reads integers per strand and is intact under the new accumulator |
| touch the message-passing model | P2 adds a term to the per-node ψ solve; the relay, the reframe and the variance laws are untouched |
| wire `node_spanning` | it is stored, unused, and `docs/calibration/S5_DESIGN_LOG.md` §1 A3 measures it as the largest single win available (0.000 → 0.758 at a 25 bp node). ⚠ But A6 must come first: spanning is a **subset** of edge-crossing, so `observable_efficiency.var_set`'s zero cross-population covariance is exactly wrong for that pair. Do it after P2 |
| revisit the κ mirror | `docs/SESSION_HANDOFF.md` §0 C4: the mirror is **consistent** and cancels, so the inference is right and only the exported scalar is mis-labelled. `TODO.md` §6 |
| adopt `assignment_mode="map"` to go green on the one failing test | ⛔ a negative control is one-sided (trap 19) and MAP is the mode that most suppresses assignment. `TODO.md` §7 |

---

## §9 THE EVIDENCE THIS PLAN RESTS ON — re-derive, do not quote

| claim | how to re-derive |
|---|---|
| the two channels have no consumers | `grep -rn "inv_length_sum\|length_sum" src/ \| grep -v "scan_payload\|native/\|substrate.py"` |
| incidences per fragment: gDNA 1.031, RNA ≈1.17–1.24 | the snippet in §2.2, over the pilot scan caches |
| the baseline's `f_gdna` column is an incidence sum | `5,370,056 + 5,373,548` vs `contained + edge_unspliced + edge_spliced` in `payload.npz` |
| `(count, Sum1/L)` is singular at equal means | build two pmfs with **realised** means matched (truncation shifts them — match after truncating), form the 2×2 from §1.2, read `cond` |
| the incidence multiplier `max(1, (w-1)/s)` | `SUM_nodes (ell-w+1)+ + n_lines·(w-1)`, all over `SUM ell`, from `nodes.feather` |
| the baseline `f_gdna` table | `docs/WIP.md`, the S5.f entry |
| 56.7 % of human nodes shorter than one fragment | `scripts/design/index_census.py INDEX --gtf GTF` |
| the efficiency tables | `scripts/design/observable_efficiency.py` |
| the opportunity formulas | `scripts/design/node_density_derivation.py`, T0 |

⭐ **Every number in this file is a pointer to a script. Run the script.**
