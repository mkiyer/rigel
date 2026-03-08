// fast_exp.h — Vectorized exp() for the EM E-step hot path
//
// Exploits three properties of the call site in em_step_kernel_range():
//   1. Input x is always <= 0  (log-sum-exp normalization subtracts max_val)
//   2. Most inputs are deeply negative (x < -708 underflows to 0.0)
//   3. Moderate accuracy suffices (posteriors are normalized; ~25 ULP is fine)
//
// Three-layer optimization:
//   Layer 1 — Early-zero skip:  (Handled by caller) skip full blocks if x < -708.0
//   Layer 2 — Inline function:  eliminates DYLD stub + enables interleaving
//   Layer 3 — NEON 2-wide:      processes two doubles per vector iteration
//
// Algorithm: Cody-Waite range reduction + degree-11 Horner polynomial (FMA)
//   x = n*ln2 + r,  |r| <= ln2/2
//   exp(r) evaluated via Horner chain with Taylor coefficients (1/k!)
//   Result scaled by 2^n via IEEE-754 exponent bit manipulation
//
// Coefficients: Taylor series 1/k! for k=2..11 on |r| < ln2/2 ≈ 0.347.
// Max truncation error: r^12/12! ≈ 5.7e-15 (~25 ULPs).
//
// References:
//   - Cody & Waite, "Software Manual for the Elementary Functions", 1980
//   - Cephes Mathematical Library (Stephen Moshier), public domain

#pragma once

#include "simd_detect.h"
#include <cmath>
#include <cstdint>
#include <cstring>  // memcpy for type-punning

#if RIGEL_HAS_NEON
#include <arm_neon.h>
#endif

namespace rigel {

// ================================================================
// Constants
// ================================================================

namespace detail {

// Range reduction: x = n*ln2 + r
// Split ln2 into high + low parts for precision (Cody-Waite)
inline constexpr double EXP_LOG2E  =  1.4426950408889634074;   // 1/ln2
inline constexpr double EXP_LN2_HI =  6.93147180369123816490e-01;
inline constexpr double EXP_LN2_LO =  1.90821492927058500170e-10;

// Early-zero cutoff: exp(x) underflows to IEEE 754 +0.0 for x < -708.39
// Use -708.0 as a conservative cutoff (bit-compatible with std::exp)
inline constexpr double EXP_CUTOFF = -708.0;

// Horner polynomial coefficients: c_k = 1/k! for k = 2..11
inline constexpr double C2  = 0.5;                        // 1/2!
inline constexpr double C3  = 1.6666666666666666e-01;     // 1/3!
inline constexpr double C4  = 4.1666666666666664e-02;     // 1/4!
inline constexpr double C5  = 8.3333333333333332e-03;     // 1/5!
inline constexpr double C6  = 1.3888888888888889e-03;     // 1/6!
inline constexpr double C7  = 1.9841269841269841e-04;     // 1/7!
inline constexpr double C8  = 2.4801587301587302e-05;     // 1/8!
inline constexpr double C9  = 2.7557319223985893e-06;     // 1/9!
inline constexpr double C10 = 2.7557319223985888e-07;     // 1/10!
inline constexpr double C11 = 2.5052108385441720e-08;     // 1/11!

} // namespace detail

// ================================================================
// Scalar fast_exp — works on all architectures
// ================================================================

inline double fast_exp_scalar(double x) {
    if (x < detail::EXP_CUTOFF) return 0.0;

    double n = std::round(x * detail::EXP_LOG2E);
    double r = x - n * detail::EXP_LN2_HI - n * detail::EXP_LN2_LO;

    // Horner evaluation of exp(r) = 1 + r*(1 + r*(C2 + ...))
    double p = detail::C11;
    p = p * r + detail::C10;
    p = p * r + detail::C9;
    p = p * r + detail::C8;
    p = p * r + detail::C7;
    p = p * r + detail::C6;
    p = p * r + detail::C5;
    p = p * r + detail::C4;
    p = p * r + detail::C3;
    p = p * r + detail::C2;
    p = p * r + 1.0;
    p = p * r + 1.0;

    // Scale by 2^n via IEEE-754 exponent bit manipulation
    uint64_t bits;
    std::memcpy(&bits, &p, sizeof(bits));
    bits += static_cast<uint64_t>(static_cast<int64_t>(n)) << 52;
    std::memcpy(&p, &bits, sizeof(p));
    return p;
}

// ================================================================
// NEON 2-wide fast_exp (ARM64)
// ================================================================

#if RIGEL_HAS_NEON

inline float64x2_t fast_exp_neon(float64x2_t x) {
    const float64x2_t log2e  = vdupq_n_f64(detail::EXP_LOG2E);
    const float64x2_t ln2_hi = vdupq_n_f64(detail::EXP_LN2_HI);
    const float64x2_t ln2_lo = vdupq_n_f64(detail::EXP_LN2_LO);
    const float64x2_t cutoff = vdupq_n_f64(detail::EXP_CUTOFF);

    // Constant hoisting for Horner polynomial
    const float64x2_t c10 = vdupq_n_f64(detail::C10);
    const float64x2_t c9  = vdupq_n_f64(detail::C9);
    const float64x2_t c8  = vdupq_n_f64(detail::C8);
    const float64x2_t c7  = vdupq_n_f64(detail::C7);
    const float64x2_t c6  = vdupq_n_f64(detail::C6);
    const float64x2_t c5  = vdupq_n_f64(detail::C5);
    const float64x2_t c4  = vdupq_n_f64(detail::C4);
    const float64x2_t c3  = vdupq_n_f64(detail::C3);
    const float64x2_t c2  = vdupq_n_f64(detail::C2);
    const float64x2_t c1  = vdupq_n_f64(1.0);

    // Underflow mask
    uint64x2_t uf_mask = vcltq_f64(x, cutoff);

    // Range reduction
    float64x2_t n = vrndnq_f64(vmulq_f64(x, log2e));
    float64x2_t r = vfmsq_f64(x, n, ln2_hi);
    r = vfmsq_f64(r, n, ln2_lo);

    // Horner polynomial evaluation (FMA chain)
    float64x2_t p = vdupq_n_f64(detail::C11);
    p = vfmaq_f64(c10, p, r);
    p = vfmaq_f64(c9,  p, r);
    p = vfmaq_f64(c8,  p, r);
    p = vfmaq_f64(c7,  p, r);
    p = vfmaq_f64(c6,  p, r);
    p = vfmaq_f64(c5,  p, r);
    p = vfmaq_f64(c4,  p, r);
    p = vfmaq_f64(c3,  p, r);
    p = vfmaq_f64(c2,  p, r);
    p = vfmaq_f64(c1,  p, r);
    p = vfmaq_f64(c1,  p, r);

    // Scale by 2^n (unsigned shift avoids signed-left-shift UB)
    int64x2_t ni = vcvtnq_s64_f64(n);
    uint64x2_t nu = vreinterpretq_u64_s64(ni);
    nu = vshlq_n_u64(nu, 52);
    int64x2_t p_bits = vreinterpretq_s64_f64(p);
    p_bits = vaddq_s64(p_bits, vreinterpretq_s64_u64(nu));
    p = vreinterpretq_f64_s64(p_bits);

    // Apply underflow mask (zeros out lanes where x < -708.0)
    p = vreinterpretq_f64_u64(
        vbicq_u64(vreinterpretq_u64_f64(p), uf_mask));

    return p;
}

#endif // RIGEL_HAS_NEON

} // namespace rigel
