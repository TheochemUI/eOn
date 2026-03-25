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

namespace eonc {

/// Base for eigenmode solvers. Holds shared state (pot, params, stats).
/// Concrete classes (Dimer, ImprovedDimer, Lanczos, AtomicGPDimer) inherit
/// this and provide compute(), getEigenvalue(), getEigenvector() as regular
/// (non-virtual) member functions. Dispatch is via std::variant
/// (EigenmodeStrategy) rather than virtual dispatch.
class LowestEigenmode {
protected:
  std::shared_ptr<Potential> pot;
  const Parameters &params;

public:
  long totalForceCalls{0};
  double statsTorque{0.0};
  double statsCurvature{0.0};
  double statsAngle{0.0};
  long statsRotations{0};
  long totalIterations{0};
  static const char MINMODE_DIMER[];
  static const char MINMODE_GPRDIMER[];
  static const char MINMODE_LANCZOS[];

  LowestEigenmode(std::shared_ptr<Potential> potPassed,
                  const Parameters &parameters)
      : pot{potPassed},
        params{parameters} {}
  ~LowestEigenmode() = default;
};

} // namespace eonc

using eonc::LowestEigenmode;
