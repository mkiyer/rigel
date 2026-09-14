// simd_detect.h — compile-time SIMD architecture detection
//
// Macros for selecting the SIMD code path across ARM64 (NEON) and x86_64 (AVX2/AVX-512), set by
// -march=native or explicit -m flags:
//   RIGEL_ARCH_ARM64      — AArch64 (Apple Silicon, AWS Graviton, etc.)
//   RIGEL_ARCH_X86_64     — x86-64 (Intel / AMD)
//   RIGEL_HAS_NEON        — ARM NEON (always 1 on AArch64)
//   RIGEL_HAS_AVX2        — x86 AVX2 (set by -march=native on capable CPUs)
//   RIGEL_HAS_FMA         — x86 FMA  (set by -march=native on capable CPUs)
//   RIGEL_HAS_AVX512F     — x86 AVX-512 Foundation

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
