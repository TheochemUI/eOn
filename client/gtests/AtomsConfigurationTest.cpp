/*
 * AtomsConfigurationTest.cpp
 *
 *  Created on: 04 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "AtomsConfigurationTest.h"

namespace tests {

AtomsConfigurationTest::AtomsConfigurationTest() {
  // TODO Auto-generated constructor stub
}

AtomsConfigurationTest::~AtomsConfigurationTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(AtomsConfigurationTest, TestMatter) {
  Parameters *parameters;
  Matter *matter = new Matter(parameters);
  matter->con2matter("pos.con");
  gpr::AtomsConfiguration a = helper_functions::eon_matter_to_atmconf(matter);
  int getFixed = count_if(a.is_frozen.getInternalVector().begin(),
                          a.is_frozen.getInternalVector().end(),
                          [](int i) { return i == 1; });
  int getFree = count_if(a.is_frozen.getInternalVector().begin(),
                         a.is_frozen.getInternalVector().end(),
                         [](int i) { return i == 0; });

  EXPECT_EQ(getFixed, matter->numberOfFixedAtoms())
      << "AtomsConf has the wrong number of fixed/frozen atoms";
  EXPECT_EQ(getFree, matter->numberOfFreeAtoms())
      << "AtomsConf has the wrong number of free atoms";
  int posindex = 0;
  for (auto i = 0; i < matter->numberOfAtoms(); i++) {
    for (auto dim = 0; dim < 3; dim++) {
      EXPECT_EQ(a.positions[posindex + dim], matter->getPosition(i, dim))
          << "AtomsConf does not match positions from Matter";
    }
    posindex = posindex + 3;
  }
  delete matter;
}

} /* namespace tests */
