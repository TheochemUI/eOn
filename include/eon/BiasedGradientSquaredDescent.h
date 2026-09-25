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
#include <utility>
#include <vector>

namespace eonc {

class BiasedGradientSquaredDescent : public SaddleSearchMethod {
public:
  BiasedGradientSquaredDescent(std::shared_ptr<Matter> matterPassed,
                               double reactantEnergyPassed,
                               const Parameters &parametersPassed)
      : SaddleSearchMethod(matterPassed->getPotential(), parametersPassed),
        saddle{std::move(matterPassed)}, reactantEnergy{reactantEnergyPassed} {
    eigenvector.resize(saddle->numberOfAtoms(), 3);
    eigenvector.setZero();
  }
  ~BiasedGradientSquaredDescent() = default;

  int run();
  double getEigenvalue();
  AtomMatrix getEigenvector();
  std::string_view describeStatus(int status) const override {
    return MinModeSaddleSearch::statusMessage(status);
  }
  int getStatus() const override { return status; }

  double eigenvalue{0.0};
  AtomMatrix eigenvector;

  std::shared_ptr<Matter> saddle;

  int status{0};

private:
  double reactantEnergy;
  eonc::log::Scoped log;
};

} // namespace eonc
