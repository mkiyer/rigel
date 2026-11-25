// solve_bindings.h — the numpy array types and the ψ argument unpacking that `_solve_impl`'s two binding TUs share
// (solve_kernel.cpp and testing/solve.cpp).
#pragma once

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "psi_kernel.h"

namespace nb = nanobind;

namespace solve_bindings {

namespace P = psi_kernel;

using Vec = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Mat = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using Cube = nb::ndarray<double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using BoolVec = nb::ndarray<bool, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using BoolMat = nb::ndarray<bool, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using IdxVec = nb::ndarray<int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using IdxMat = nb::ndarray<int32_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using KindMat = nb::ndarray<int8_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using FlagVec = nb::ndarray<uint16_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// An array the kernel writes into, or whose data it keeps past the cast: the caller's own array, never a converted
// copy.
template <class A, class H>
A exact(H&& h) {
    return nb::cast<A>(std::forward<H>(h), false);
}

// The fitted gDNA arm's curve, (log_rho, logP): two (G,) arrays, G >= 1, the caller's own, since the kernel keeps
// their data. `who` names the call in a refusal.
inline void check_arm_curve(const Vec& log_rho, const Vec& logP, const char* who) {
    if (log_rho.shape(0) < 1 || logP.shape(0) != log_rho.shape(0))
        throw std::invalid_argument(std::string(who) + ": the gDNA arm's curve is not two (G,) arrays");
}

// A gDNA arm as a call receives it, None or a tuple of `n` arrays, the curve first; its arrays, or none.
inline std::vector<Vec> arm_tuple(nb::object gdna, size_t n, const char* who) {
    if (gdna.is_none()) return {};
    nb::tuple t = nb::cast<nb::tuple>(gdna);
    if (nb::len(t) != n) throw std::invalid_argument(std::string(who) + ": the gDNA arm is not a " +
                                                     std::to_string(n) + "-tuple of arrays");
    std::vector<Vec> arrays;
    for (size_t i = 0; i < n; ++i) arrays.push_back(exact<Vec>(t[i]));
    check_arm_curve(arrays[0], arrays[1], who);
    return arrays;
}

// The fitted gDNA arm's arrays as a call receives them: None, or (log_rho, logP, mass, eff).
struct ArmArrays {
    Vec log_rho, logP, mass, eff;
    P::Arm arm;
    ArmArrays(nb::object gdna, int m) {
        const std::vector<Vec> a = arm_tuple(gdna, 4, "psi");
        if (a.empty()) return;
        log_rho = a[0]; logP = a[1]; mass = a[2]; eff = a[3];
        if (static_cast<int>(mass.shape(0)) != m || static_cast<int>(eff.shape(0)) != m)
            throw std::invalid_argument("psi: the gDNA arm's support is not one (mass, eff) per slot");
        arm = P::Arm{log_rho.data(), logP.data(), static_cast<int>(log_rho.shape(0)), mass.data(), eff.data()};
    }
    ArmArrays(const ArmArrays&) = delete;
};

// The delivered composition rows as a call receives them: None or (m, K).
struct RowPriors {
    bool present;
    Mat rows;
    RowPriors(nb::object row_logprior, int K) : present(!row_logprior.is_none()) {
        if (present) rows = nb::cast<Mat>(row_logprior);
        if (present && static_cast<int>(rows.shape(1)) != K)
            throw std::invalid_argument("psi: a prior row is not on the solve grid");
    }
    const double* row(int i) const { return present ? rows.data() + static_cast<size_t>(i) * rows.shape(1) : nullptr; }
};

// The delivered RNA level rows as a call receives them (`simplex_logodds.CubeRows`), unpacked over m slots.
inline P::Rows unpack_rows(int m, int K, IdxVec& cube_slot, Mat& cube_pos, BoolVec& has_pos, Mat& cube_neg,
                           BoolVec& has_neg, Vec& cube_u, Vec& cube_total, Vec& cube_opportunity, Vec& cube_rho) {
    const int d = static_cast<int>(cube_slot.shape(0));
    if (d && (static_cast<int>(cube_pos.shape(1)) != K || static_cast<int>(cube_u.shape(0)) != K))
        throw std::invalid_argument("psi: a delivered row is not on the solve grid");
    return P::unpack_rows(m, K, d, cube_slot.data(), cube_pos.data(), has_pos.data(), cube_neg.data(), has_neg.data(),
                          cube_u.data(), cube_total.data(), cube_opportunity.data(), cube_rho.data());
}

// ψ's inputs as `psi_solve` and the gates' `psi_cube` receive them: the per-slot counts, admissibility and reference
// composition; the fitted gDNA arm; the delivered composition rows and RNA level rows; and the grid they are read on.
struct PsiInputs {
    Vec u_pos, u_neg, fg_ref, fpos_ref, fneg_ref, lam;
    BoolVec allow_pos, allow_neg;
    // the delivered RNA level rows' arrays, which `delivered` points into
    IdxVec cube_slot; Mat cube_pos, cube_neg; BoolVec cube_has_pos, cube_has_neg; Vec cube_u, cube_total, cube_opp, cube_rho;
    ArmArrays arm;
    RowPriors priors;
    P::Rows delivered;
    P::Grid grid;

    PsiInputs(Vec u_pos_, Vec u_neg_, BoolVec allow_pos_, BoolVec allow_neg_, Vec fg_ref_, Vec fpos_ref_,
              Vec fneg_ref_, double kappa, Vec lam_, nb::object gdna, nb::object row_logprior, IdxVec cube_slot_,
              Mat cube_pos_, BoolVec cube_has_pos_, Mat cube_neg_, BoolVec cube_has_neg_, Vec cube_u_,
              Vec cube_total_, Vec cube_opportunity_, Vec cube_rho_, int n_tilt)
        : u_pos(u_pos_), u_neg(u_neg_), fg_ref(fg_ref_), fpos_ref(fpos_ref_), fneg_ref(fneg_ref_), lam(lam_),
          allow_pos(allow_pos_), allow_neg(allow_neg_), cube_slot(cube_slot_), cube_pos(cube_pos_),
          cube_neg(cube_neg_), cube_has_pos(cube_has_pos_), cube_has_neg(cube_has_neg_), cube_u(cube_u_),
          cube_total(cube_total_), cube_opp(cube_opportunity_), cube_rho(cube_rho_),
          arm(gdna, slots()), priors(row_logprior, K()),
          delivered(unpack_rows(slots(), K(), cube_slot, cube_pos, cube_has_pos, cube_neg, cube_has_neg, cube_u,
                                cube_total, cube_opp, cube_rho)),
          grid(lam.data(), K(), kappa, n_tilt, &arm.arm) {}
    PsiInputs(const PsiInputs&) = delete;

    int slots() const { return static_cast<int>(u_pos.shape(0)); }
    int K() const { return static_cast<int>(lam.shape(0)); }

    // Slot i's inputs to ψ.
    P::SlotInputs slot(int i) const {
        const bool built = arm.arm.built();
        return P::SlotInputs{u_pos.data()[i], u_neg.data()[i], fg_ref.data()[i], fpos_ref.data()[i], fneg_ref.data()[i],
                             allow_pos.data()[i], allow_neg.data()[i], built, built ? arm.arm.shift(i) : 0.0,
                             priors.row(i), delivered.index[i] >= 0 ? &delivered.rows[delivered.index[i]] : nullptr};
    }
};

}  // namespace solve_bindings
