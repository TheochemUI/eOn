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
#include "MinModeSaddleSearch.h"
#include "SaddleSearchMethod.h"
#include <vector>

class DynamicsSaddleSearch : public SaddleSearchMethod {
public:
  DynamicsSaddleSearch(std::shared_ptr<Matter> matterPassed,
                       std::shared_ptr<Parameters> parametersPassed)
      : SaddleSearchMethod(nullptr, parametersPassed),
        product{std::make_shared<Matter>(*matterPassed)},
        reactant{std::make_shared<Matter>(*matterPassed)},
        saddle{matterPassed} {
    this->pot = matterPassed->getPotential();
    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
    log = spdlog::get("combi");
  };
  ~DynamicsSaddleSearch() = default;

  int run(void);
  double getEigenvalue();
  AtomMatrix getEigenvector();

  int refineTransition(std::vector<std::shared_ptr<Matter>> MDSnapshots,
                       std::shared_ptr<Matter> product);

  double eigenvalue;
  AtomMatrix eigenvector;

  double time;

  std::shared_ptr<Matter> product;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;

  int status;

private:
  std::shared_ptr<spdlog::logger> log;
};
