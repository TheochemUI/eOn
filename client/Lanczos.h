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
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"

// Lanczos method to find the lowest curvature mode
class Lanczos : public LowestEigenmode {

public:
  Lanczos(std::shared_ptr<Matter> matter, std::shared_ptr<Parameters> params,
          std::shared_ptr<Potential> pot);
  ~Lanczos() = default;
  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  AtomMatrix lowestEv;
  double lowestEw;
  shared_ptr<spdlog::logger> log;
};
