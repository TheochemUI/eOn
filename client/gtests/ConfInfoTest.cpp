/*
 * ConfInfoTest.cpp
 *
 *  Created on: 04 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../HelperFunctions.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "ConfInfoTest.h"

namespace tests {

ConfInfoTest::ConfInfoTest() {
  // TODO Auto-generated constructor stub
}

ConfInfoTest::~ConfInfoTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(ConfInfoTest, TestMatter) {
  Parameters *parameters;
  Matter *matter = new Matter(parameters);
  matter->con2matter("./data/pos.con");
  ConfInfo a = helper_functions::eon_matter_to_init_confinfo(matter);
  EXPECT_EQ(a.conf_fro.getNumPoints(), matter->numberOfFixedAtoms())
      << "Frozen config has the wrong number of fixed/frozen atoms";
  delete matter;
}

} /* namespace tests */
