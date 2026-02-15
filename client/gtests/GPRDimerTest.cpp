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
/*
 * GPRDimerTest.cpp
 *
 *  Created on: 07 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../AtomicGPDimer.h"
#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Job.h"
#include "../Matter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "GPRDimerTest.h"

namespace tests {

GPRDimerTest::GPRDimerTest() {
  // TODO Auto-generated constructor stub
}

GPRDimerTest::~GPRDimerTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRDimerTest, TestMatter) {
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");
  AtomMatrix mode;
  Parameters parameters;
  parameters.load("config.ini");
  auto pot = helper_functions::makePotential(parameters);
  auto initial = std::make_shared<Matter>(pot, parameters);
  auto saddle = std::make_shared<Matter>(pot, parameters);
  initial->con2matter(reactantFilename);
  saddle->con2matter(displacementFilename);
  mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
  auto saddleSearch = std::make_unique<MinModeSaddleSearch>(
      saddle, mode, initial->getPotentialEnergy(), parameters, pot);
  auto minModeMethod = std::make_unique<AtomicGPDimer>(saddle, parameters, pot);
  minModeMethod->compute(saddle, mode);
  cout << minModeMethod->getEigenvalue();
}

} /* namespace tests */
