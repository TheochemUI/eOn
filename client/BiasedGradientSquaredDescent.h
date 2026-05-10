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
#include "EonLogger.h"

#include "Eigen.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "SaddleSearchMethod.h"
#include <vector>

namespace eonc {

class BiasedGradientSquaredDescent : public SaddleSearchMethod {
public:
  BiasedGradientSquaredDescent(const std::shared_ptr<Matter> &matterPassed,
                               double reactantEnergyPassed,
                               const Parameters &parametersPassed)
      : SaddleSearchMethod(matterPassed->getPotential(), parametersPassed),
        eigenvalue{0.0},
        saddle{matterPassed},
        reactantEnergy{reactantEnergyPassed} {
    eigenvector.resize(saddle->numberOfAtoms(), 3);
    eigenvector.setZero();
  }
  ~BiasedGradientSquaredDescent() = default;

  SaddleStatus run() override;
  double getEigenvalue() override;
  AtomMatrix getEigenvector() override;
  SaddleStatus getStatus() const override { return status; }

  double eigenvalue;
  AtomMatrix eigenvector;

  std::shared_ptr<Matter> saddle;

  SaddleStatus status{SaddleStatus::Good};

private:
  double reactantEnergy;
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::BiasedGradientSquaredDescent;
