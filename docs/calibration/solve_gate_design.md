# The §6B Solve-Gate — derived, implemented, and EMPIRICALLY REFUTED

**Status (2026-07-21): REFUTED — do not re-attempt.** The DOF solve-gate ([`message_system_derivation.md`](message_system_derivation.md)
§6B; roadmap item A) was taken through the full mantra — design & derive, implement, measure — and the
measurement refuted its central premise. This note is the durable record so the idea is not re-tried.

---

## 1. The idea (§6B)

Pass-0 today solves every **structurally** solvable node (`solvable = (fp|fn) & mass>0`) — a node with a strand
and mass deconvolves its gDNA-vs-RNA split. §6B proposed a stricter **DOF** gate: a node is solvable iff every
**free axis** is *identified* (has ≥1 nonzero-precision source among {strand tilt, messages, prior}) —
single-strand needs `λ`; AMBIG needs `λ` and `θ`. An **un**identified node would **skip** its solve and **keep
its signature-binary init** (`f_g=1`, max variance), "deferring to the Phase-2 gDNA hyperprior" without a pass-0
phantom to fight.

## 2. Design analysis (what "skip" actually produces)

Grounding it in the code surfaced the load-bearing facts:
- The **AMBIG (G3) init is `f_g=1` at max variance** (`node_geometry.init_beliefs`) — so "skip" reverts an
  unidentified AMBIG node to **all-gDNA**.
- `solvable` also drives `locked` → `struct_lock` + the sweep, so the DOF gate had to be a **separate write-back
  gate** (`write_solvable`), never fed into `solvable` (else unidentified nodes would be wrongly marked
  composition-*certain*).
- `tau_th` is **not relayed** (θ evidence = strand tilt only); the identification signals available at the
  write-back are `tau0_lam/tau0_th` (the I_strand seed) + the combined message precisions `prec_g/prec_p/prec_n`.
- On **zero-gDNA** data `I_strand=0` even when stranded (`N_gdna=0 ⇒ σ²_d=∞ ⇒ disc=0`), so those nodes are
  genuinely unidentifiable in pass-0 → they skip → `f_g=1`. **On a zero-gDNA library, `f_g=1` is the phantom.**

This predicted a `gdna_none` regression standalone — the design's honest cost, "paid back by the prior."

## 3. Implementation (flag-gated, `RIGEL_SOLVE_GATE`)

A separate write-back gate: `write_solvable = solvable & dof_id`, with `λ` identified ⟺ `tau0_lam>0 OR prec_g>0
OR prec_p+prec_n>0`, `θ` identified ⟺ `tau0_th>0 OR prec_p>0 OR prec_n>0`, AMBIG needing both. Byte-identical
with the flag off.

## 4. Measurement (`pass0_bench` on `ambig_dense_10mb`, 32 scenarios) — the refutation

| arm | mean mwae | gdna_none | verdict |
|---|---|---|---|
| refit=0, gate **OFF** (baseline) | **0.2030** | ~0.01–0.47 | — |
| refit=0, gate **ON** | 0.2135 (**+0.010**) | **rises** (0.047→0.096) | regresses standalone |
| **refit=1, gate OFF** (hyperprior) | **0.0998** | ~0.00–0.01 | **the hyperprior HALVES the error** |
| refit=1, gate **ON** | 0.1244 (**+0.025**) | **rises hard** (0.009→0.119) | regresses *with* the prior too |

15/32 scenarios worse, **0 better**, at refit=1. **The §6B premise is empirically false:** deferring a node to
`f_g=1` does **not** let the prior resolve it better — the prior resolves an imperfectly-**solved** node (which
carries the observed-density information) *far* better than a deferred all-gDNA init. Even the confident-wrong
slice §6B targeted is better handled by *solving then prior-refining* than by *deferring*.

## 5. Implications

1. **A is dead.** Reverted (no lingering knob). `solvable` stays structural; pass-0 solves every strand+mass node.
2. **§6B's skip rule is retracted** — flag it refuted in `message_system_derivation.md`.
3. **The hyperprior is the lever, and it is already wired** (`calib_refit_iters`): refit=1 alone takes mwae
   0.2030 → **0.0998**. It works **best with pass-0 solving everything**, not with a solve-gate deferral. So the
   real next step is improving the **hyperprior itself** (the DNA-prior track, `dna_prior_session_resume.md`) —
   not gating pass-0.
4. **"Destination decides" is retired as a pass-0 mechanism.** Emission is honest (always emit, `pr→0` when
   vacuous); a vacuous node already sends weak messages and its own solve is reference-anchored — that is
   sufficient. The over-commitment §6B feared does not, in practice, beat solving + the prior.
