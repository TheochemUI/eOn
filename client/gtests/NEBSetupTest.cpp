/*
 * NEBSetupTest.cpp
 *
 *  Created on: 1 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Log.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "NEBSetupTest.h"

namespace tests {

NEBSetupTest::NEBSetupTest() {
    reactantFilename = helper_functions::getRelevantFile("reactant.con");
    productFilename = helper_functions::getRelevantFile("product.con");

    parameters = std::make_unique<Parameters>();
    parameters->potential = "morse_pt";
    parameters->nebImages = 7;
    parameters->LogPotential = false;
    log_init(parameters.get(), (char *)"test.log");

    initmatter = std::make_unique<Matter>(parameters.get());
    finalmatter = std::make_unique<Matter>(parameters.get());

    initmatter->con2matter(reactantFilename);
    finalmatter->con2matter(productFilename);
}

NEBSetupTest::~NEBSetupTest() {
}

TEST_F(NEBSetupTest, TestMatter) {
  // Setup the run
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  // Setup counters
  EXPECT_EQ(imgArray.back().getForces(), this->finalmatter->getForces())
      << "Forces do not match";
}

} /* namespace tests */
