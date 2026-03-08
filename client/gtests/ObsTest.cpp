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

#include <algorithm>

#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "ObsTest.h"

namespace tests {

ObsTest::ObsTest() {
  // TODO Auto-generated constructor stub
}

ObsTest::~ObsTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(ObsTest, TestMatter) {
  string confile("pos.con");
  Parameters parameters;
  parameters.potential_options.potential = PotType::MORSE_PT;
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
