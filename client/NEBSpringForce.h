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

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"
#include "SafeMath.h"
#include <variant>
#include <vector>

namespace eonc::neb {

/// Result of spring force computation for a single image.
struct SpringResult {
  AtomMatrix forceSpringPar; // Parallel (tangent) spring force
  AtomMatrix forceSpring;    // Full spring force (needed for DNEB)
};

/// Uniform spring constant for all images.
struct UniformSpring {
  double ksp;

  SpringResult compute(long i, const AtomMatrix &tangent, double distNext,
                       double distPrev, const AtomMatrix &posDiffNext,
                       const AtomMatrix &posDiffPrev,
                       const std::shared_ptr<Matter> &image) const;
};

/// Energy-weighted spring constants (variable per segment).
struct WeightedSpring {
  std::vector<double> springConstants;

  SpringResult compute(long i, const AtomMatrix &tangent, double distNext,
                       double distPrev) const;
};

/// Onsager-Machlup action-based springs (Mandelli & Parrinello 2021).
struct OnsagerMachlupSpring {
  double base_k;
  std::vector<AtomMatrix> L_vecs;

  SpringResult compute(long i, const AtomMatrix &tangent,
                       const AtomMatrix &posNext, const AtomMatrix &posPrev,
                       const AtomMatrix &pos,
                       const std::shared_ptr<Matter> &image) const;
};

using SpringStrategy =
    std::variant<UniformSpring, WeightedSpring, OnsagerMachlupSpring>;

/// Build the appropriate spring strategy from parameters and current path
/// state.
SpringStrategy
buildSpringStrategy(const Parameters &params,
                    const std::vector<std::shared_ptr<Matter>> &path,
                    long numImages, int atoms, double maxEnergy, double E_ref);

} // namespace eonc::neb
