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
  string optimizer("cg");
  SaddleSearchMethod *saddleSearch;
  AtomMatrix mode;
  LowestEigenmode *minModeMethod;
  Parameters *parameters = new Parameters;
  // parameters->load("matter.ini");
  parameters->potential = PotType::LJ;
  parameters->optMethod = optimizer;
  parameters->optConvergedForce = optConvergedForce;
  parameters->dimerConvergedAngle = 0.001;
  parameters->saddleMinmodeMethod = LowestEigenmode::MINMODE_DIMER;
  Matter *initial = new Matter(parameters);
  Matter *displacement = new Matter(parameters);
  Matter *saddle = new Matter(parameters);
  initial->con2matter(reactantFilename);
  saddle->con2matter(displacementFilename);
  mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
  saddleSearch = new MinModeSaddleSearch(
      saddle, mode, initial->getPotentialEnergy(), parameters);
  minModeMethod = new ImprovedDimer(saddle, parameters, saddle->getPotential());
  minModeMethod->compute(saddle, mode);
  cout << minModeMethod->getEigenvalue();
  delete minModeMethod;
  delete initial;
  delete displacement;
  delete saddle;
  delete parameters;
}

} /* namespace tests */
