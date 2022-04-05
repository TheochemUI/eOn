/*
 * SCGMatterTest.cpp
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
#include "NEBTest.h"

namespace tests {

NEBTest::NEBTest() {
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

NEBTest::~NEBTest() {
}

TEST_F(NEBTest, TestMatter) {
  // Setup the run
    auto *params = new Parameters;
    params->potential = "morse_pt";
    params->nebImages = 7;
    params->LogPotential = false;
  auto initPath = helper_functions::prepInitialPath(params);
  auto imgArray = std::get<std::vector<Matter> >(initPath);
  auto tangentArray = std::get<std::vector<AtomMatrix> >(initPath);
  auto projForceArray = tangentArray; // Initially the same
  // Setup counters
  const int nimages = this->parameters->nebImages;
  auto extremumPosition = std::vector<double>(2*nimages+1, 0);
  auto extremumEnergy = std::vector<double>(2*nimages+1, 0);
  auto extremumCurvature = std::vector<double>(2*nimages+1, 0);
  int numExtrema = 0;
  std::cout<<imgArray[2].getPositions();
}

} /* namespace tests */
