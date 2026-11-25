// em_prior.h — the locus EM's grouped gDNA/RNA calibration prior.
//
// Calibration hands each locus two pseudo-fragment counts, one of gDNA and one of RNA. Every M-step adds them to the
// locus's gDNA and RNA pools and shares the RNA pseudocount over the RNA components.
#pragma once

#include <algorithm>
#include <cmath>

#include "constants.h"

namespace rigel {

// The calibration prior for one locus: pseudo-fragments of gDNA and of RNA, and how the RNA share is allocated.
struct AggregatePrior {
    double gdna_prior_fragments = 0.0;
    double rna_prior_fragments  = 0.0;
    int    gdna_index = -1;
    bool   has_gdna_candidate = false;

    // Per-component allocation weight for the RNA prior; nullptr allocates in proportion to the evidence each
    // component already carries. Both are the additive pseudocount `a_i = rna_prior · w_i / Σw`; the evidence rule
    //     out[i] = raw[i] · (1 + rna_prior/rna_count) = raw[i] + rna_prior · raw[i]/rna_count
    // is the case `w_i = raw[i]`. They differ at `raw[i] == 0`: evidence weights make `out[i] = 0` an absorbing
    // state, so a component with no evidence is never revived, and a strictly positive weight has no such state.
    const double* component_rna_prior_weight = nullptr;
};

inline double nonnegative_finite(double x) {
    return (std::isfinite(x) && x > 0.0) ? x : 0.0;
}

// The sum of nonnegative_finite(v[i]) over the RNA components (every i but the gDNA index).
inline double rna_sum(const double* v, int gdna_index, int n_components) {
    double s = 0.0;
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_index) continue;
        s += nonnegative_finite(v[i]);
    }
    return s;
}

// Split the locus into gDNA and RNA, add the calibration prior to each, and share the RNA prior over the RNA
// components in proportion to the evidence each already carries (or to the caller's allocation weight).
//
// Every RNA component, annotated transcript or synthetic span, receives its share: under evidence weights the prior
// is the common factor `(1 + rna_prior/rna_count)` and moves only the gDNA:RNA split, while a withheld share would
// move RNA between components the data cannot tell apart. Within one call the RNA components sum to
// `rna_count + rna_prior` and gDNA to `gdna_count + gdna_prior`, exactly.
inline void apply_grouped_prior_update(
    const double* raw_counts,
    const double* carried_state,
    const AggregatePrior& aggregate_prior,
    double* out_counts,
    int n_components)
{
    const int gdna_index = aggregate_prior.gdna_index;
    const bool has_gdna = aggregate_prior.has_gdna_candidate && gdna_index >= 0 && gdna_index < n_components;

    // The RNA pool's evidence and carried totals, and its allocation weight when the caller supplies one.
    const double rna_count = rna_sum(raw_counts, gdna_index, n_components);
    const double rna_carried = carried_state != nullptr ? rna_sum(carried_state, gdna_index, n_components) : 0.0;
    const double* rna_prior_weight = aggregate_prior.component_rna_prior_weight;
    const double rna_weight = rna_prior_weight != nullptr ? rna_sum(rna_prior_weight, gdna_index, n_components) : 0.0;
    const bool weighted = rna_weight > EM_LOG_EPSILON;

    // The gDNA pseudocount needs a gDNA component to land on.
    const double gdna_prior = has_gdna ? aggregate_prior.gdna_prior_fragments : 0.0;
    // The RNA pseudocount lands on the RNA components, so it is not gated on a gDNA component (an all-spliced locus
    // has none). It is dropped only when nothing can receive it, tested on the denominator the chosen branch divides
    // by (the weight total, else the count, else the carried state): any other total could keep a live prior at a
    // locus with zero count but carried alpha (VBEM), multiply it by zero and move the gDNA:RNA split.
    double rna_prior = aggregate_prior.rna_prior_fragments;
    const double prior_recipients = weighted ? rna_weight : ((rna_count > EM_LOG_EPSILON) ? rna_count : rna_carried);
    if (prior_recipients <= EM_LOG_EPSILON) rna_prior = 0.0;

    const double gdna_count = has_gdna ? nonnegative_finite(raw_counts[gdna_index]) : 0.0;
    const double rna_total = rna_count + rna_prior;

    std::fill(out_counts, out_counts + n_components, 0.0);
    if (has_gdna) out_counts[gdna_index] = gdna_count + gdna_prior;

    // Weighted allocation: out[i] = raw[i] + rna_prior · w_i / Σw, so the pool sums to rna_count + rna_prior. A zero
    // weight total is not `weighted` and takes the evidence rule below, which keeps the prior and the split.
    if (weighted) {
        const double scale = rna_prior / rna_weight;
        for (int i = 0; i < n_components; ++i) {
            if (i == gdna_index) continue;
            out_counts[i] = nonnegative_finite(raw_counts[i]) + nonnegative_finite(rna_prior_weight[i]) * scale;
        }
        return;
    }

    // Evidence weights: each RNA component's evidence scaled by rna_total / rna_count. With zero evidence, reachable
    // under VBEM (which passes `alpha` as the carried state), the prior-only pool is shared in proportion to the
    // carried alpha; the warm start passes nullptr and never reaches that.
    const bool by_evidence = rna_count > EM_LOG_EPSILON;
    if (!by_evidence && !(rna_carried > EM_LOG_EPSILON && carried_state != nullptr)) return;
    const double* share = by_evidence ? raw_counts : carried_state;
    const double inv = 1.0 / (by_evidence ? rna_count : rna_carried);
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_index) continue;
        out_counts[i] = rna_total * nonnegative_finite(share[i]) * inv;
    }
}

}  // namespace rigel
