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
#include "EigenmodeStrategy.h"
#include "Matter.h"
#include "Optimizer.h"
#include "SaddleSearchMethod.h"
#include "StatusTypes.h"

#include <string>

namespace eonc {

class MinModeSaddleSearch : public SaddleSearchMethod {

public:
  /// Status codes / messages live in StatusTypes.h
  /// (eonc::SaddleStatus). Use SaddleStatus::Good etc. and
  /// statusMessage(SaddleStatus) at call sites.

  MinModeSaddleSearch(std::shared_ptr<Matter> matterPassed,
                      const AtomMatrix &modePassed, double reactantEnergyPassed,
                      const Parameters &parametersPassed,
                      std::shared_ptr<Potential> potPassed);
  ~MinModeSaddleSearch() = default;
  AtomMatrix getEigenvector() override; // lowest eigenmode
  double getEigenvalue() override;      // estimate for the lowest eigenvalue
  SaddleStatus getStatus() const override { return status; }
  int getIterationCount() const override { return iteration; }
  int getForceCalls() const override { return forcecalls; }

  SaddleStatus run() override;
  SaddleStatus run(long max_iterations_override);

  int forcecalls{0};
  int iteration{0};
  SaddleStatus status{SaddleStatus::Good};

private:
  AtomMatrix mode;
  AtomMatrix initialTangent_;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<EigenmodeStrategy>
      minModeMethod; // shared with the objective func
  double reactantEnergy;
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::MinModeSaddleSearch;
