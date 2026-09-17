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
 * ObsTest.cpp
 *
 *  Created on: 05 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "ObsTest.h"
#include "eon/GPRHelpers.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eonc_test_aliases.hpp"

namespace tests {

TEST_F(ObsTest, TestMatter) {
  string confile("pos.con");
  Parameters parameters;
  ParametersLoadAccess::potential_options(parameters).potential =
      PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(parameters);
  auto matter = std::make_shared<Matter>(pot, parameters);
  matter->con2matter(confile);
  gpr::Observation o = eonc::helpers::eon_matter_to_init_obs(matter.get());
  EXPECT_EQ(o.R.extractEigenMatrix(), matter->getPositions())
      << "Positions do not match";
  EXPECT_EQ(o.G.extractEigenMatrix() * -1, matter->getForces())
      << "Forces do not match";
  EXPECT_EQ(o.E[0], matter->getPotentialEnergy())
      << "Potential energy does not match";
}

} /* namespace tests */
