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
#include "eon/PbcSimd.h"

#include <cmath>
#include <cstddef>

#ifndef EON_PBC_SIMD_SCALAR_DEFINED
#define EON_PBC_SIMD_SCALAR_DEFINED
namespace {

void wrapMinimumImageScalar(double *data, size_t begin, size_t n) {
  for (size_t i = begin; i < n; ++i) {
    const double x = data[i];
    data[i] = x - std::floor(x + 0.5);
  }
}

void wrapLegacyUnitScalar(double *data, size_t begin, size_t n) {
  for (size_t i = begin; i < n; ++i) {
    const double y = data[i] + 1.0;
    data[i] = y - std::trunc(y);
  }
}

} // namespace
#endif // EON_PBC_SIMD_SCALAR_DEFINED

#ifdef WITH_HIGHWAY

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "client/PbcSimd.cpp"
#include "hwy/foreach_target.h" // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace eonc::pbc {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

void wrapMinimumImage(double *HWY_RESTRICT data, size_t n) {
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d);
  const auto half = hn::Set(d, 0.5);
  size_t i = 0;
  if (lanes > 0) {
    for (; i + lanes <= n; i += lanes) {
      const auto v = hn::LoadU(d, data + i);
      hn::StoreU(hn::Sub(v, hn::Floor(hn::Add(v, half))), d, data + i);
    }
  }
  wrapMinimumImageScalar(data, i, n);
}

void wrapLegacyUnit(double *HWY_RESTRICT data, size_t n) {
  const hn::ScalableTag<double> d;
  const size_t lanes = hn::Lanes(d);
  const auto one = hn::Set(d, 1.0);
  size_t i = 0;
  if (lanes > 0) {
    for (; i + lanes <= n; i += lanes) {
      const auto y = hn::Add(hn::LoadU(d, data + i), one);
      hn::StoreU(hn::Sub(y, hn::Trunc(y)), d, data + i);
    }
  }
  wrapLegacyUnitScalar(data, i, n);
}

} // namespace HWY_NAMESPACE
} // namespace eonc::pbc
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace eonc::pbc {

HWY_EXPORT(wrapMinimumImage);
HWY_EXPORT(wrapLegacyUnit);

void wrapMinimumImage(double *data, size_t n) {
  if (n == 0) {
    return;
  }
  HWY_DYNAMIC_DISPATCH(wrapMinimumImage)(data, n);
}

void wrapLegacyUnit(double *data, size_t n) {
  if (n == 0) {
    return;
  }
  HWY_DYNAMIC_DISPATCH(wrapLegacyUnit)(data, n);
}

} // namespace eonc::pbc
#endif // HWY_ONCE

#else // !WITH_HIGHWAY

namespace eonc::pbc {

void wrapMinimumImage(double *data, size_t n) {
  wrapMinimumImageScalar(data, 0, n);
}

void wrapLegacyUnit(double *data, size_t n) {
  wrapLegacyUnitScalar(data, 0, n);
}

} // namespace eonc::pbc

#endif // WITH_HIGHWAY
