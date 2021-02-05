/*
 * ObsTest.cpp
 *
 *  Created on: 05 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

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
  Parameters *parameters;
  Matter *matter = new Matter(parameters);
  matter->con2matter("./data/pos.con");
  Observation o = helper_functions::eon_matter_to_init_obs(matter);
  delete matter;
}

} /* namespace tests */
