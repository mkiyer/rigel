export const meta = {
  name: 'cliff-precision-design',
  description: 'Derive the honest cross-enrichment-cliff precision law + produce the implementation plan',
  phases: [
    { title: 'Derive', detail: '3 independent derivations of the honest cross-cliff precision, each MC-validated' },
    { title: 'Critique', detail: '2 adversaries stress-test the derivations against the constraints' },
    { title: 'Plan', detail: 'synthesize ONE settled law + a concrete implementation plan' },
  ],
}

const ENV = 'source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel && OMP_NUM_THREADS=1'
const REPO = '/Users/mkiyer/proj/rigel'

const SUBSTRATE = `
Rigel calibration pass-0 message-passing solver. We have ROOT-CAUSED a precision bug and need to DERIVE the fix.

THE BUG (root-caused, reproduced to machine precision by a prior audit workflow). When a message is transported
across a large enrichment cliff, its PRECISION is preserved when it should decay. Concretely (bp_solver.py,
_relay/_transport): a message is reframed by r = ρ_tot(dst)/ρ_tot(src) — this scales the MODE (density) so
compositions are comparable. σ²_transfer = Var(log r) is meant to damp the precision across the cliff, but:
  (a) it is EXEMPTED (=0) on matched/graft reframes (correct for the mode's enrichment cancellation, WRONG for
      the precision), and
  (b) even where applied it is ~1/n (negligible on well-counted nodes),
so a message crossing a 407× cliff arrives at FULL precision. Anchor case (gdna_gdna300_ss_0.99_capture_on, node
1909): an exon, oracle f_g=0.985 (mostly gDNA), mass 66,544, its own message-free solve fg_loc=0.953 CORRECT but
WEAK (τ=1.6). Its flanking boundary 1910 (mass 14, more RNA) has a real 47-count spliced junction; the reframe
r≈407 (ρ_tot(1909)=34 vs ρ_tot(1910)=0.08, dominated by 1909's GDNA) stretches that boundary's RNA density UP and
delivers it as a confident RNA+ message (precision ~17-26, mode f_pos=0.718) → f_g collapses to 0.51. A channel
ablation shows BOTH the RNA-measurement stream (cm_p) AND the composition λ-stream (c_tau) are corrupted the same
way; only dropping BOTH recovers node 1909 to 0.91.

THE OWNER'S INVARIANT (the design target). Precision is fundamentally constrained by (1) discrete COUNTS over (2)
a discrete LENGTH in base pairs. Enrichment ratios (up to 1000×) scale the MODES so compositions are comparable
— but enrichment must NEVER preserve/scale the PRECISION as if the measurement were made at the destination's
scale. "Going across a big enrichment cliff must change the scale of precisions": a small-count / small-density
source, stretched by a large r onto a large node, must arrive WEAKENED in proportion to the stretch — because a
node and its 1000×-denser neighbour are almost certainly DIFFERENT kinds of region (enriched/expressed vs
depleted), so the imputation "we share a composition" is far more suspect the bigger the cliff. The reframe r is
dominated by the MAJORITY component (gDNA at 1909), so applying it to a MINORITY component (RNA) inflates the
minority — the composition-mismatch is the heart of it.

HARD CONSTRAINTS on the fix:
- Respect the counts/length the measurement ORIGINATED from (honest precision). NO magic numbers.
- Must NOT over-correct: on UNSTRANDED AMBIG nodes the composition τ_own=0 and cross-node messages are the ONLY
  information (this is where the σ²_transfer win — unstr-capON 0.387→0.163 — comes from). A blanket own-count cap
  kills that. The fix must be composition/mass-mismatch-aware, damping the CLIFF-crossing confidence without
  killing legitimate propagation on genuinely-similar neighbours.
- The MODE reframe (enrichment cancellation) is CORRECT and stays; this is about the PRECISION.

SUBSTRATE (do NOT read docs/calibration/archive/):
- ${REPO}/src/rigel/calibration/bp_solver.py — _relay (~line 395), _transport (~495), the reframe r, s2t /
  transfer_logvar, logvar_tot, the THREE streams (full pg/pp/pn, measurement mg/mp/mn, composition τ), the
  additive fuses (:478, _fuse_add :545/:566). The graft s2t=0 exemption (:421 relay, :509 transport).
- ${REPO}/src/rigel/calibration/enrichment_frame.py — composition_logvar (= 1/n + [(1/E_g−1/E_r)/B]²·Var(f_g)),
  reframe_density, message_precision, transfer_logvar, the message-variance laws.
- ${REPO}/docs/calibration/message_variance_derivation.md — the settled M1-M5 variance laws (Var(log r) = M5).
- ${REPO}/docs/calibration/SESSION_2026_07_24_HANDOFF_4.md §6 — the finding.
- MC template: ${REPO}/scripts/debug/message_variance_mc.py (pure numpy, the arbiter for variance laws).

METHOD: DERIVE the honest cross-cliff precision law, then MC-VALIDATE it (write a pure-numpy MC that constructs a
source measurement at a small scale, transports it across a large enrichment cliff to a large node, and shows the
delivered precision matches the honest value — near-zero for a big mismatched cliff, full for a matched same-scale
transport). Run it: ${ENV} python <your_script>. Ground every claim in numbers.
`

phase('Derive')
const framings = [
  'Derive it as MODEL/IMPUTATION uncertainty: the message assumes "src and dst share composition", whose validity degrades with the enrichment/expression difference. Quantify Var(log ρ_c^dst | src) as a function of the cliff — what honest term makes precision decay with the composition-mismatch (NOT merely Var(log r) of well-counted totals)? Derive from the reframe r being dominated by the majority component so the minority is stretched by an unrelated ratio.',
  'Derive it from COUNTS-OVER-LENGTH first principles: the source measured n_src counts over its length; when reframed to a dst of very different density, what is the honest effective count of that measurement AT the destination (the discrete resolving power)? Show that a measurement stretched by a large r loses effective count in proportion to the stretch, and that a matched same-density transport keeps it. Connect to the intron-factory / spliced count semantics.',
  'Derive it as the correct σ²_transfer: re-examine whether σ²_transfer=Var(log r) is the RIGHT object or whether the graft exemption + the 1/n-smallness are the two bugs. Should the graft exemption apply to the mode ONLY (keep composition cancellation) while the PRECISION always pays Var(log r) — AND should Var(log r) itself carry a composition-mismatch term (not just the counting 1/n)? Derive the composition-mismatch component of Var(log r).',
]
const derivations = await parallel(
  framings.map((f, i) => () =>
    agent(
      `${SUBSTRATE}\n\nYou are independent derivation agent #${i + 1}. ${f}\n\n`
      + `Produce: (1) the DERIVED honest cross-cliff precision law (a formula, with each term's meaning and its `
      + `count/length origin); (2) an MC validation you WROTE and RAN (${REPO}/scratchpad/cliff_mc_${i + 1}.py) — `
      + `report the measured numbers, including the anchor-like case (small mismatched source across a ~400× cliff `
      + `→ delivered precision should be ~0) AND a matched same-scale transport (precision preserved); (3) how it `
      + `avoids over-damping the unstranded/AMBIG propagation where messages are the only info. Return a thorough `
      + `PROSE report (no rigid schema) with the formula, the MC numbers, and the exact bp_solver change points.`,
      { label: `derive#${i + 1}`, phase: 'Derive', effort: 'high' }
    )
  )
)

phase('Critique')
const derJson = derivations.filter(Boolean).map((d, i) => `\n===== DERIVATION #${i + 1} =====\n${d}`).join('\n')
const critics = await parallel([
  'Attack each derived law on the HARD CONSTRAINTS: does it introduce a magic number? does it over-damp the unstranded/AMBIG arm (τ_own=0, messages-only — the M5 win)? does it break the count-honesty at the origin? MC-test the failure modes you find (write + run your own script).',
  'Attack the COMPLETENESS: the ablation showed BOTH cm_p and c_tau are corrupted, and node 1909 also has a stretched-wrong MODE (f_pos=0.718). Does the proposed precision fix ALONE recover the anchor, or is a MODE fix also needed? Which of the 3 derivations, if implemented, actually recovers node 1909 to ~0.95 without a mode change — verify by a numeric estimate. Flag if a mode fix is unavoidable.',
].map((angle, i) => () =>
  agent(
    `${SUBSTRATE}\n\nYou are CRITIC #${i + 1}. ${angle}\n\nThe derivations:\n${derJson}\n\n`
    + `Instrument/MC where needed (${REPO}/scratchpad/cliff_crit_${i + 1}.py, RUN it). Return a PROSE critique with `
    + `numbers: which derivation survives, what breaks, and whether a mode fix is also required.`,
    { label: `critic#${i + 1}`, phase: 'Critique', effort: 'high' }
  )
))

phase('Plan')
const PLAN_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['law', 'rationale', 'mc_validated', 'implementation_steps', 'mode_fix_needed', 'over_correction_guard', 'risks'],
  properties: {
    law: { type: 'string', description: 'the settled honest cross-cliff precision law (formula + each term)' },
    rationale: { type: 'string', description: 'why it is honest (counts/length origin) and matches the owner invariant' },
    mc_validated: { type: 'string', description: 'the MC evidence it reproduces the anchor (big mismatched cliff → ~0 precision) and preserves a matched transport' },
    implementation_steps: { type: 'array', items: { type: 'string' }, description: 'ordered concrete bp_solver / enrichment_frame changes (file:line-level)' },
    mode_fix_needed: { type: 'string', description: 'is a MODE fix also required to recover the anchor, or does the precision fix alone suffice? evidence.' },
    over_correction_guard: { type: 'string', description: 'how it avoids regressing the unstranded/AMBIG arm (the M5 win)' },
    risks: { type: 'array', items: { type: 'string' } },
  },
}
const critJson = critics.filter(Boolean).map((c, i) => `\n===== CRITIQUE #${i + 1} =====\n${c}`).join('\n')
const plan = await agent(
  `${SUBSTRATE}\n\nYou are the ADJUDICATOR. Synthesize the 3 derivations and 2 critiques into ONE settled honest `
  + `cross-cliff precision LAW + a concrete implementation plan. Resolve disagreements by re-running MC yourself `
  + `(${ENV} python ...). The law must be honest (counts/length), magic-number-free, must recover the anchor (node `
  + `1909 → ~0.95), and must NOT regress the unstranded/AMBIG arm. State clearly whether a MODE fix is also needed.\n\n`
  + `DERIVATIONS:\n${derJson}\n\nCRITIQUES:\n${critJson}`,
  { label: 'plan', phase: 'Plan', schema: PLAN_SCHEMA, effort: 'high' }
)

return { derivations: derivations.filter(Boolean), critiques: critics.filter(Boolean), plan }
