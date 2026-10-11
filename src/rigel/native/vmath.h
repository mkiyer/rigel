/**
 * vmath.h — `log` and `exp` over an array, by their series, written as plain loops the compiler vectorises.
 *
 * The honest reader's node loops spend their time in two logs and one exp per lattice node (the Poisson pair
 * and the trapezoid sum's exponentials); libm's scalar calls neither vectorise nor pipeline across nodes.
 * These are the textbook series with their natural truncation, nothing fitted:
 *
 *   log x = 2·atanh(s) + e·ln 2,  x = 2^e·m, m ∈ [√½, √2), s = (m − 1)/(m + 1), |s| ≤ 3 − 2√2 ≈ 0.1716,
 *           2·atanh(s) = 2s·Σ_{k≥0} s^{2k}/(2k + 1), kept to k = 11 (the next term is below 2^-54 of the sum);
 *   exp x = 2^k·e^r,  k = round(x/ln 2), |r| ≤ ½ ln 2 ≈ 0.347,  e^r = Σ_{j≤14} r^j/j! (the next term is below 2^-54);
 *           x below log(DBL_MIN) is 0, and so is −∞. The reduction r = x − k·ln 2 uses ln 2 in two parts, the
 *           first with 21 significant bits so that k·hi is exact, the second their exact difference; what remains
 *           is the double M_LN2's own distance from ln 2, δ ≈ 2.3e-17, multiplied by k.
 *
 * `vlog` is within 2 ulp of libm; `vexp`'s relative error is within 3·2^-52 + (|x|/ln 2 + ½)·δ (the gate:
 * tests/native/test_vmath.py, against numpy on a million points, δ derived there from ln 2 to forty digits). The
 * reader's sums take exponentials within T = 36 nats of their maximum, where that is under 2e-15 of a term; the
 * reader's gates hold the kernel to its specification at 1e-7, nine orders above.
 * Arguments to `vlog` are positive and normal (the Poisson rates are); `vexp` takes any real or −∞.
 * Pure C++17, no Python.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace rigel::vmath {

constexpr int kLogTerms = 12;  // s^0 … s^22 of 2·atanh(s)/(2s)
constexpr int kExpTerms = 15;  // r^0 … r^14 of e^r

struct Series {
    double log_c[kLogTerms];  // 1/(2k + 1)
    double exp_c[kExpTerms];  // 1/j!
    constexpr Series() : log_c{}, exp_c{} {
        for (int k = 0; k < kLogTerms; ++k) log_c[k] = 1.0 / double(2 * k + 1);
        double f = 1.0;
        for (int j = 0; j < kExpTerms; ++j) { exp_c[j] = 1.0 / f; f *= double(j + 1); }
    }
};
inline constexpr Series kSeries{};

inline void vlog(const double* x, double* out, int n) {
    for (int i = 0; i < n; ++i) {
        uint64_t bits;
        std::memcpy(&bits, &x[i], sizeof bits);
        double e = double(int64_t((bits >> 52) & 0x7ffULL) - 1023);
        uint64_t mb = (bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL;  // the mantissa in [1, 2)
        double m;
        std::memcpy(&m, &mb, sizeof m);
        const bool high = m > M_SQRT2;
        m = high ? 0.5 * m : m;  // now in [√½, √2)
        e = high ? e + 1.0 : e;
        const double s = (m - 1.0) / (m + 1.0), s2 = s * s;
        double p = kSeries.log_c[kLogTerms - 1];
        for (int k = kLogTerms - 2; k >= 0; --k) p = p * s2 + kSeries.log_c[k];
        out[i] = 2.0 * s * p + e * M_LN2;
    }
}

// ln 2 in two parts: `hi` is M_LN2 with its low 32 mantissa bits cleared (21 significant bits), `lo` the exact remainder
struct Ln2Split {
    double hi, lo;
    Ln2Split() {
        uint64_t bits;
        double v = M_LN2;
        std::memcpy(&bits, &v, sizeof bits);
        bits &= ~0xffffffffULL;
        std::memcpy(&hi, &bits, sizeof hi);
        lo = v - hi;
    }
};
inline const Ln2Split kLn2{};

inline void vexp(const double* x, double* out, int n) {
    const double floor_x = std::log(std::numeric_limits<double>::min());
    for (int i = 0; i < n; ++i) {
        const double xi = x[i] < floor_x ? floor_x : x[i];  // below the smallest normal: the result rounds to 0 below
        const double k = std::nearbyint(xi * (1.0 / M_LN2));
        const double r = (xi - k * kLn2.hi) - k * kLn2.lo;
        double p = kSeries.exp_c[kExpTerms - 1];
        for (int j = kExpTerms - 2; j >= 0; --j) p = p * r + kSeries.exp_c[j];
        const int64_t ki = int64_t(k);
        uint64_t bits = uint64_t(ki + 1023) << 52;  // 2^k, k ≥ −1022 by the floor
        double scale;
        std::memcpy(&scale, &bits, sizeof scale);
        out[i] = x[i] < floor_x ? 0.0 : p * scale;
    }
}

}  // namespace rigel::vmath
