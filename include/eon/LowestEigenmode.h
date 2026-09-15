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
/// Dispatch is virtual so EigenmodeStrategy does not change layout when
/// WITH_GPRD adds AtomicGPDimer.
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
  static const char MINMODE_DAVIDSON[];

  LowestEigenmode(std::shared_ptr<Potential> potPassed,
                  const Parameters &parameters)
      : pot{potPassed}, params{parameters} {}
  virtual ~LowestEigenmode() = default;

  virtual void compute(std::shared_ptr<Matter> matter,
                       AtomMatrix initialDirection) = 0;
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;
};

} // namespace eonc
