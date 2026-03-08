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

#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"

// dimer method to find the lowest curvature mode
class MonteCarlo {

public:
  MonteCarlo(std::shared_ptr<Matter> const matterIn, const Parameters &paramsIn)
      : matter{matterIn},
        params{paramsIn} {}
  ~MonteCarlo() = default;

  void run(int numSteps, double temperature, double stepSize);

private:
  std::shared_ptr<Matter> matter;
  const Parameters &params;
  eonc::log::Scoped log;
};
