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
#include "EonLogger.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "SaddleSearchMethod.h"
#include <vector>

namespace eonc {

class DynamicsSaddleSearch : public SaddleSearchMethod {
public:
  DynamicsSaddleSearch(std::shared_ptr<Matter> matterPassed,
                       const Parameters &parametersPassed)
      : SaddleSearchMethod(nullptr, parametersPassed),
        product{std::make_shared<Matter>(*matterPassed)},
        reactant{std::make_shared<Matter>(*matterPassed)},
        saddle{matterPassed} {
    this->pot = matterPassed->getPotential();
    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
  };
  ~DynamicsSaddleSearch() = default;

  int run();
  double getEigenvalue();
  AtomMatrix getEigenvector();
  std::string_view describeStatus(int status) const override {
    return MinModeSaddleSearch::statusMessage(status);
  }

  int refineTransition(const std::vector<std::shared_ptr<Matter>> &snapshots,
                       std::shared_ptr<Matter> product);

  double eigenvalue{0.0};
  AtomMatrix eigenvector;

  double time{0.0};

  std::shared_ptr<Matter> product;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;

  int status{0};

private:
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::DynamicsSaddleSearch;
