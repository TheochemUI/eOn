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
 * ImpDimerTest.cpp
 *
 *  Created on: 07 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../HelperFunctions.h"
#include "../ImprovedDimer.h"
#include "../Job.h"
#include "../Matter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "ImpDimerTest.h"

namespace tests {

ImpDimerTest::ImpDimerTest() {
  // TODO Auto-generated constructor stub
}

ImpDimerTest::~ImpDimerTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(ImpDimerTest, TestMatter) {
  double optConvergedForce = 0.001;
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");
  AtomMatrix mode;
  Parameters parameters;
  parameters.potential_options.potential = PotType::LJ;
  parameters.optimizer_options.method = OptType::CG;
  parameters.optimizer_options.converged_force = optConvergedForce;
  parameters.dimer_options.converged_angle = 0.001;
  parameters.saddle_search_options.minmode_method =
      LowestEigenmode::MINMODE_DIMER;
  auto pot = eonc::helpers::makePotential(parameters);
  auto initial = std::make_shared<Matter>(pot, parameters);
  auto displacement = std::make_shared<Matter>(pot, parameters);
  auto saddle = std::make_shared<Matter>(pot, parameters);
  initial->con2matter(reactantFilename);
  saddle->con2matter(displacementFilename);
  mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
  auto saddleSearch = std::make_unique<MinModeSaddleSearch>(
      saddle, mode, initial->getPotentialEnergy(), parameters, pot);
  auto minModeMethod = std::make_unique<ImprovedDimer>(saddle, parameters,
                                                       saddle->getPotential());
  minModeMethod->compute(saddle, mode);
  cout << minModeMethod->getEigenvalue();
}

} /* namespace tests */
