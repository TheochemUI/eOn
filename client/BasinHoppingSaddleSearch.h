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

class BasinHoppingSaddleSearch : public SaddleSearchMethod {
public:
  BasinHoppingSaddleSearch(std::shared_ptr<Matter> reactant,
                           std::shared_ptr<Matter> displacement,
                           std::shared_ptr<Potential> potPassed,
                           std::shared_ptr<Parameters> parametersPassed)
      : SaddleSearchMethod(potPassed, parametersPassed),
        reactant{std::make_shared<Matter>(*reactant)},
        saddle{displacement} {
    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
    log = spdlog::get("combi");
  }
  ~BasinHoppingSaddleSearch() = default;

  int run(void);
  double getEigenvalue();
  AtomMatrix getEigenvector();

  double eigenvalue;
  AtomMatrix eigenvector;

  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;
  std::shared_ptr<Matter> product;

  int status;

private:
  std::shared_ptr<spdlog::logger> log;
};
