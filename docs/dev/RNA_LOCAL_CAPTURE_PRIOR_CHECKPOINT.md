# Local capture-prior implementation checkpoint

2026-10-09. Isolated prototype approved by the owner. No change to production source,
tests, goldens or the installed extension; no commit or push. This is a reference
implementation and a bounded screen, not a release candidate.

## What was implemented

The prototype is in `.cache/rigel_runs/2026-10-09_capture_prior/`.

| File | Responsibility |
|---|---|
| `local_prior.py` | Integrate local density evidence against a fixed Gamma background and the proper enrichment prior; return the continuous correction score, posterior component probability and numerical error estimate. |
| `evidence.py` | Adapt analytic curves and the existing exact latent-count observation reference. A small independent RNA-assay factor exercises neighbour information on both-strand and unstranded controls. |
| `test_local_prior.py` | Mathematical and locality falsifications, including an independent physical-density integral and a direct raw-observation integral. |
| `mutate.py` | Execute faulty copies of the reader and record the tests that catch them. |
| `screen.py` | Separate the slab change from the odds change on frozen own evidence; calculate exact Poisson sampling risks without random simulation or parameter fitting. |
| `freeze_messages.py` | Reconstruct current typed factors from the fresh saved count contexts, checking existing single-strand assembly identity. No calibration refit or alternative message policy. |
| `complete.py` | Feed those factors into the existing isolated native density oracle and price a complete readout, with a per-object wall-time budget. |
| `verify.py` | Check receipts, source identity and the retained unresolved numerical gates. |

The reader accepts **a likelihood function**, not an inferred DNA count. Its production
candidate entry point receives only that function, background mean and background shape.
Annotation count and remote RNA maxima are absent. Diagnostic controls retain the old
bounded slab and old odds only to isolate the two changes; they are not proposed
production options.

Final checks: 28 prototype tests pass; all 13 deliberately faulty implementations
are caught; the documentation boundary gate passes all 11 tests. The production
suite was not rerun for this source-preserving prototype; its last verified count
remains 2,814. `verification.json` checks production source and native identity and
retains the three unresolved source-interpolation failures outside production.

The slab uses the background **mean**, with density equal to mean times capture factor.
It does not sample a separate background rate and multiply that by capture. Background
estimation is held fixed and is not implemented here. Gamma backgrounds, including their
point-mass limit, are the numerical reference interface, not a new background-fitting
decision. An exact zero or unidentified background is rejected explicitly; no floor is
substituted and no library capture status is created.

The reference uses scalar SciPy quadrature on log density. It splits at supplied curve
features but integrates both tails; a split is not a statistical ceiling. Known flat
evidence uses its analytic limit. An evidence curve's arbitrary multiplicative scale
cancels from the score. The raw observation reference integrates RNA out rather than
turning a posterior DNA estimate into a Poisson observation.

The independent small RNA-assay control means `J ~ Poisson(a*r)`, where `r` is the target's
RNA amount and `a` is the measured opportunity ratio. Its contribution is integrated
with the same nuisance reference. This is an exactly checkable local-information fixture,
not a proposal to replace the existing messages with a new measurement model.

## What the measurements establish

The single home for numerical outcomes is **PROPER LOCAL CAPTURE PRIOR** in
[ISSUES](../ISSUES.md); the derivation and limits are **Proper local capture reference
and continuous correction** in [EQUATIONS](../EQUATIONS.md).

The three failed-first tests reproduce substantive defects in the archived bounded
reader on actual pure-DNA observations, where its Poisson likelihood is legitimate.
They do not use a missing function or an import failure as evidence of a model defect.
The candidate is then checked on no evidence, zero counts, missing own opportunity with
informative neighbours, concentrated capture, uncertain background, continuity, units,
likelihood scaling and annotation/disconnected-RNA invariance. Mutation receipts contain
actual assertion failures, not predicted test sensitivity.

The archive screen selects **every** object meeting its declared small-count diagnostic
criteria, plus the previously implicated zero-DNA boundary. Capture targets come from
archived simulator expected yields. Realized DNA counts appear only as diagnostic labels.
The separate pure-DNA risk calculation averages over the Poisson distribution instead of
selecting a fortunate random sample. No setting is optimized to those measurements.

The complete-consumer screen reuses current observations and message factors. The old
native intensity integrator is loaded in an ordinary isolated process using `python -S`
and explicit module paths; the installed count extension is never replaced. Both-strand
factors retain their witness exclusions and separate physical coordinates. The deferred
unstranded cost case uses the same fixed background as the stranded control; it cannot
validate an unstranded background estimator.

## What it does not establish

The scalar arithmetic passes, but the complete both-strand consumer exceeds the declared
cost budget. The results do not authorize copying the nested reference integrator into
production. The small completed cases are not a speed comparison with the current reader.

Current message tables still carry the known finite-support and narrow-source interpolation
defects. Outer quadrature refinement does not repair them. The prior itself has no remote
RNA or annotation input; this is conditional locality, not a proof that every upstream
message builder is local. Unit covariance is also not a fragment-length-gap benchmark.

Sparse genuine capture is conservatively estimated, and small false corrections remain.
No claim of zero false capture, copy-number immunity, or full recovery from one DNA read
is made. The mixture odds are a reference assumption, not learned probe prevalence.

No panel calibration, EM, transcript/gene comparison or real-library quantification is
run for this prototype. The accepted count foundation and its previous receipts remain
unchanged. Background fitting, capture transfer between components, junction pricing and
the production integration checkpoint remain open.

## Reproduction

Use the `rigel` environment and `OMP_NUM_THREADS=1`. Check `uptime` before each experiment.
The reference gates and archive screen take seconds. Each complete-consumer invocation
has a 60-second diagnostic limit; do not extend that into a whole-panel run.

```bash
python -m pytest .cache/rigel_runs/2026-10-09_capture_prior/test_local_prior.py -q
READER_ARM=old python -m pytest .cache/rigel_runs/2026-10-09_capture_prior/test_local_prior.py -q -k 'disconnected or annotation or ceiling'
python .cache/rigel_runs/2026-10-09_capture_prior/mutate.py
python .cache/rigel_runs/2026-10-09_capture_prior/screen.py
python .cache/rigel_runs/2026-10-09_capture_prior/freeze_messages.py
```

The `READER_ARM=old` command deliberately fails its three tests. For a complete readout,
use the saved evidence extension while avoiding the environment's editable-import finder:

```bash
PYTHONPATH="$PWD/.cache/rigel_runs/2026-10-08_inner_budget/site:$CONDA_PREFIX/lib/python3.12/site-packages" \
  python -S .cache/rigel_runs/2026-10-09_capture_prior/complete.py \
  --case test_str_on --slot 1000
```

Use `--label refined --rtol 1e-8` for the completed-case outer-tolerance check. The native
inner tolerance remains unchanged; this is not a certificate of total error. The recorded
source hashes identify the exact extension and typed input arrays.

## Next implementation boundary

Keep this reader as an oracle. The authorized shared outer-quadrature experiment is now
complete and rejected for cost; see **SHARED CAPTURE INTEGRATION** in ISSUES. The owner's
[external review packet](RNA_CAPTURE_RELEASE_REVIEW.md) asks for the smallest complete
consumer before more numerical machinery is developed. No more informative prior,
strand-only shortcut, confidence cutoff or simpler-but-different both-strand model is
selected here. Retain source precision, background uncertainty and expected-yield checks
before a release A/B. The active work order is [the fix plan](RNA_SHORT_FIX_PLAN.md).
