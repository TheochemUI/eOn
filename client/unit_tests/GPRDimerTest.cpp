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

#include "GPRDimerTest.h"
#include "eon/AtomicGPDimer.h"
#include "eon/GPRHelpers.h"
#include "eon/HelperFunctions.h"
#include "eon/Job.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"
#include "eonc_test_aliases.hpp"

namespace tests {

TEST_F(GPRDimerTest, TestMatter) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");
  AtomMatrix mode;
  Parameters parameters;
  parameters.load("config.ini");
  auto pot = eonc::helpers::makePotential(parameters);
  auto initial = std::make_shared<Matter>(pot, parameters);
  auto saddle = std::make_shared<Matter>(pot, parameters);
  initial->con2matter(reactantFilename);
  saddle->con2matter(displacementFilename);
  mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
  auto saddleSearch = std::make_unique<MinModeSaddleSearch>(
      saddle, mode, initial->getPotentialEnergy(), parameters, pot);
  auto minModeMethod = std::make_unique<AtomicGPDimer>(saddle, parameters, pot);
  minModeMethod->compute(saddle, mode);
  cout << minModeMethod->getEigenvalue();
}

} /* namespace tests */
