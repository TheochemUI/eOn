/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "eon/ForceNorm.h"

#include <algorithm>
#include <cstddef>

#ifdef WITH_HIGHWAY

// foreach_target.h re-includes this file once per SIMD target. The path is
// relative to the project include root.
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "client/ForceNorm.cpp"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace eonc {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

double MaxFreeAtomForceNorm(const double *HWY_RESTRICT forces,
                            const double *HWY_RESTRICT fixed, long nAtoms) {
  if (forces == nullptr || nAtoms <= 0) {
    return 0.0;
  }
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d);
  const size_t n = static_cast<size_t>(nAtoms);
  const auto half = hn::Set(d, 0.5);
  auto vmax = hn::Zero(d);

  size_t i = 0;
  for (; i + lanes <= n; i += lanes) {
    auto x = hn::Zero(d);
    auto y = hn::Zero(d);
    auto z = hn::Zero(d);
    hn::LoadInterleaved3(d, forces + 3 * i, x, y, z);
    // Mul then Add, not MulAdd. A contracted FMA rounds once; the scalar
    // tail rounds each product and the sum, and both paths must agree.
    const auto nrm =
        hn::Sqrt(hn::Add(hn::Add(hn::Mul(x, x), hn::Mul(y, y)), hn::Mul(z, z)));
    if (fixed != nullptr) {
      auto fx = hn::Zero(d);
      auto fy = hn::Zero(d);
      auto fz = hn::Zero(d);
      hn::LoadInterleaved3(d, fixed + 3 * i, fx, fy, fz);
      const auto allFixed = hn::And(
          hn::Gt(fx, half), hn::And(hn::Gt(fy, half), hn::Gt(fz, half)));
      vmax = hn::Max(vmax, hn::IfThenElseZero(hn::Not(allFixed), nrm));
    } else {
      vmax = hn::Max(vmax, nrm);
    }
  }

  return std::max(hn::ReduceMax(d, vmax),
                  detail::maxFreeAtomForceNormScalar(
                      forces, fixed, static_cast<long>(i), nAtoms));
}

} // namespace HWY_NAMESPACE
} // namespace eonc
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace eonc {

HWY_EXPORT(MaxFreeAtomForceNorm);

double maxFreeAtomForceNorm(const double *forces, const double *fixed,
                            long nAtoms) {
  return HWY_DYNAMIC_DISPATCH(MaxFreeAtomForceNorm)(forces, fixed, nAtoms);
}

} // namespace eonc
#endif // HWY_ONCE

#else // !WITH_HIGHWAY

namespace eonc {

double maxFreeAtomForceNorm(const double *forces, const double *fixed,
                            long nAtoms) {
  return detail::maxFreeAtomForceNormScalar(forces, fixed, 0, nAtoms);
}

} // namespace eonc

#endif // WITH_HIGHWAY
