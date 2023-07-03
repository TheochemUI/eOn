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
  double optConvergedForce = 0.001;
  string potential("morse_pt");
  string reactantFilename("pos.con");
  string displacementFilename("displacement.con");
  string modeFilename("direction.dat");
  string optimizer("lbfgs");
  SaddleSearchMethod *saddleSearch;
  AtomMatrix mode;
  LowestEigenmode *minModeMethod;
  Parameters *parameters = new Parameters;
  parameters->load("config.ini");
  // parameters->potential = potential;
  // parameters->optMethod = optimizer;
  // parameters->optConvergedForce = optConvergedForce;
  // parameters->dimerConvergedAngle = 0.0873;
  // parameters->saddleMinmodeMethod = LowestEigenmode::MINMODE_GPRDIMER;
  Matter *initial = new Matter(parameters);
  Matter *saddle = new Matter(parameters);
  initial->con2matter(reactantFilename);
  saddle->con2matter(displacementFilename);
  mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
  saddleSearch = new MinModeSaddleSearch(
      saddle, mode, initial->getPotentialEnergy(), parameters);
  minModeMethod = new AtomicGPDimer(saddle, parameters);
  minModeMethod->compute(saddle, mode);
  cout << minModeMethod->getEigenvalue();
  delete minModeMethod;
  delete initial;
  delete saddle;
  delete parameters;
}

} /* namespace tests */
