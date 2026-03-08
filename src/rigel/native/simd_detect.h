// simd_detect.h — Portable SIMD architecture detection and dispatch helpers
//
// Provides compile-time macros and runtime CPU detection for selecting
// the optimal SIMD code path across ARM64 (NEON) and x86_64 (SSE2/AVX2/AVX-512).
//
// Compile-time macros (set by -march=native or explicit -m flags):
//   RIGEL_ARCH_ARM64      — AArch64 (Apple Silicon, AWS Graviton, etc.)
//   RIGEL_ARCH_X86_64     — x86-64 (Intel / AMD)
//   RIGEL_HAS_NEON        — ARM NEON (always 1 on AArch64)
//   RIGEL_HAS_SSE2        — x86 SSE2 (always 1 on x86-64)
//   RIGEL_HAS_AVX2        — x86 AVX2 (set by -march=native on capable CPUs)
//   RIGEL_HAS_FMA         — x86 FMA  (set by -march=native on capable CPUs)
//   RIGEL_HAS_AVX512F     — x86 AVX-512 Foundation
//   RIGEL_SIMD_DOUBLES    — doubles per SIMD register (1/2/4/8)
//
// Runtime detection (for portable builds without -march=native):
//   rigel::cpu_has_avx2()    — true if runtime CPU supports AVX2+FMA
//   rigel::cpu_has_avx512f() — true if runtime CPU supports AVX-512F
//
// Usage for Phase 3 vectorized exp:
//
//   #include "simd_detect.h"
//
//   // Native build (-march=native): compile-time selection
//   #if RIGEL_HAS_AVX2
//       void exp_vec(double* d, int n) { /* AVX2 path, 4-wide */ }
//   #elif RIGEL_HAS_NEON
//       void exp_vec(double* d, int n) { /* NEON path, 2-wide */ }
//   #else
//       void exp_vec(double* d, int n) { /* scalar fallback */ }
//   #endif
//
//   // Portable build: runtime dispatch via CPUID
//   // Use __attribute__((target("avx2,fma"))) on x86 to compile AVX2
//   // code even without global -mavx2 flag.

#pragma once

// ================================================================
// 1. Architecture detection
// ================================================================

#if defined(__aarch64__) || defined(_M_ARM64)
    #define RIGEL_ARCH_ARM64  1
    #define RIGEL_ARCH_X86_64 0
#elif defined(__x86_64__) || defined(_M_X64)
    #define RIGEL_ARCH_ARM64  0
    #define RIGEL_ARCH_X86_64 1
#else
    #define RIGEL_ARCH_ARM64  0
    #define RIGEL_ARCH_X86_64 0
#endif

// ================================================================
// 2. ISA feature detection (compile-time)
// ================================================================
//
// These reflect what the compiler is currently targeting.
// With -march=native they match the host CPU.
// Without it (portable), they reflect the baseline (NEON or SSE2).

// ARM NEON — mandatory on AArch64, always available
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define RIGEL_HAS_NEON 1
#else
    #define RIGEL_HAS_NEON 0
#endif

// x86 SSE2 — mandatory on x86-64, always available
#if defined(__SSE2__) || (RIGEL_ARCH_X86_64 && defined(_MSC_VER))
    #define RIGEL_HAS_SSE2 1
#else
    #define RIGEL_HAS_SSE2 0
#endif

// x86 AVX2 — Haswell (2013), Zen1 (2017), and later
#ifdef __AVX2__
    #define RIGEL_HAS_AVX2 1
#else
    #define RIGEL_HAS_AVX2 0
#endif

// x86 FMA — typically available alongside AVX2
#ifdef __FMA__
    #define RIGEL_HAS_FMA 1
#else
    #define RIGEL_HAS_FMA 0
#endif

// x86 AVX-512 Foundation — Skylake-X (2017), Zen4 (2022), and later
#ifdef __AVX512F__
    #define RIGEL_HAS_AVX512F 1
#else
    #define RIGEL_HAS_AVX512F 0
#endif

// ================================================================
// 3. SIMD vector width
// ================================================================
//
// Best available vector width as number of double-precision floats.
// Used for loop tiling and buffer alignment decisions.

#if RIGEL_HAS_AVX512F
    inline constexpr int RIGEL_SIMD_DOUBLES = 8;   // 512-bit
#elif RIGEL_HAS_AVX2
    inline constexpr int RIGEL_SIMD_DOUBLES = 4;   // 256-bit
#elif RIGEL_HAS_NEON || RIGEL_HAS_SSE2
    inline constexpr int RIGEL_SIMD_DOUBLES = 2;   // 128-bit
#else
    inline constexpr int RIGEL_SIMD_DOUBLES = 1;   // scalar
#endif

// ================================================================
// 4. Runtime CPU feature detection (x86_64 only)
// ================================================================
//
// For portable builds (RIGEL_PORTABLE=ON, no -march=native), the
// compile-time macros above reflect the baseline ISA (SSE2 only).
// Use these runtime checks to dispatch to AVX2/AVX-512 code paths
// compiled via __attribute__((target(...))) or separate TUs.
//
// On ARM64, runtime detection is unnecessary — NEON is always on.

#if RIGEL_ARCH_X86_64
    #if defined(__GNUC__) || defined(__clang__)
        #include <cpuid.h>
    #elif defined(_MSC_VER)
        #include <intrin.h>
    #endif
#endif

namespace rigel {

#if RIGEL_ARCH_X86_64

inline bool cpu_has_avx2() {
#if defined(__GNUC__) || defined(__clang__)
    unsigned int eax, ebx, ecx, edx;
    // Leaf 7, subleaf 0: extended feature flags
    if (!__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx))
        return false;
    if (!((ebx >> 5) & 1))  // bit 5 = AVX2
        return false;
    // Leaf 1: FMA is ECX bit 12
    if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx))
        return false;
    return (ecx >> 12) & 1;
#elif defined(_MSC_VER)
    int info[4];
    __cpuidex(info, 7, 0);
    if (!((info[1] >> 5) & 1))
        return false;
    __cpuid(info, 1);
    return (info[2] >> 12) & 1;
#else
    return false;
#endif
}

inline bool cpu_has_avx512f() {
#if defined(__GNUC__) || defined(__clang__)
    unsigned int eax, ebx, ecx, edx;
    if (!__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx))
        return false;
    return (ebx >> 16) & 1;  // bit 16 = AVX-512F
#elif defined(_MSC_VER)
    int info[4];
    __cpuidex(info, 7, 0);
    return (info[1] >> 16) & 1;
#else
    return false;
#endif
}

#else  // ARM64 or unknown — no runtime detection needed

inline bool cpu_has_avx2()    { return false; }
inline bool cpu_has_avx512f() { return false; }

#endif

// Human-readable string for the active SIMD path (compile-time)
inline const char* simd_description() {
#if RIGEL_HAS_AVX512F
    return "x86_64 AVX-512F (8-wide double)";
#elif RIGEL_HAS_AVX2 && RIGEL_HAS_FMA
    return "x86_64 AVX2+FMA (4-wide double)";
#elif RIGEL_HAS_AVX2
    return "x86_64 AVX2 (4-wide double)";
#elif RIGEL_HAS_SSE2
    return "x86_64 SSE2 (2-wide double)";
#elif RIGEL_HAS_NEON
    return "ARM64 NEON (2-wide double)";
#else
    return "scalar (no SIMD)";
#endif
}

} // namespace rigel
