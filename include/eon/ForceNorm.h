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
#pragma once

#include <algorithm>
#include <cmath>

namespace eonc {

namespace detail {

/// Scalar max of per-atom Euclidean norms on `[begin, nAtoms)`.
/// Fully fixed atoms (`fixed` entries all > 0.5) contribute nothing.
inline double maxFreeAtomForceNormScalar(const double *forces,
                                         const double *fixed, long begin,
                                         long nAtoms) {
  double best = 0.0;
  if (forces == nullptr || nAtoms <= 0 || begin >= nAtoms) {
    return best;
  }
  if (begin < 0) {
    begin = 0;
  }
  for (long i = begin; i < nAtoms; ++i) {
    if (fixed != nullptr) {
      const double *mask = fixed + 3 * i;
      if (mask[0] > 0.5 && mask[1] > 0.5 && mask[2] > 0.5) {
        continue;
      }
    }
    const double *row = forces + 3 * i;
    const double xx = row[0] * row[0];
    const double yy = row[1] * row[1];
    const double zz = row[2] * row[2];
    best = std::max(best, std::sqrt((xx + yy) + zz));
  }
  return best;
}

} // namespace detail

/// Max Euclidean norm over N x 3 row-major force rows.
/// An atom is ignored when `fixed` is non-null and all three of its mask
/// entries are greater than 0.5 (fully fixed). `fixed == nullptr` keeps
/// every atom. Returns 0 when `nAtoms` is not positive.
[[nodiscard]] double maxFreeAtomForceNorm(const double *forces,
                                          const double *fixed, long nAtoms);

} // namespace eonc
