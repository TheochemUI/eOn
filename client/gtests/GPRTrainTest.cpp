/*
 * GPRTrainTest.cpp
 *
 *  Created on: 27 Mar 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/GPRPotential/GPRPotential.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "GPRTrainTest.h"
#include "../subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"

namespace tests {

GPRTrainTest::GPRTrainTest() {
  // TODO Auto-generated constructor stub
}

GPRTrainTest::~GPRTrainTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRTrainTest, TestMatter) {
  aux::ProblemSetUp problem_setup;
  string confile{"pos.con"};
  Parameters *parameters = new Parameters;
  parameters->potential = "morse_pt";
  parameters->potential = "morse_pt";
  log_init(parameters, (char *)"test.log");
  Matter *matter = new Matter(parameters);
  matter->con2matter("pos.con");
  auto config_data = helper_functions::eon_matter_to_frozen_conf_info(matter,  1);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  auto R_init = std::get<gpr::Coord>(config_data);
  // Setup GPR
  // GPRPotential pot{parameters};
  gpr::GPRSetup gpr_parameters;
  aux::AuxiliaryFunctionality aux_func;
  gpr::GaussianProcessRegression *gprfunc = new gpr::GaussianProcessRegression();
  auto all_obs = helper_functions::eon_matter_to_init_obs(matter);

  gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gpr_parameters.jitter_sigma2 = 0.;
  gprfunc->setParameters(gpr_parameters);

  gprfunc->getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
  gprfunc->getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
  gprfunc->getSexpAtCovarianceFunction()->setConfInfo(atoms_config);

  gprfunc->getConstantCovarianceFunction()->setConstSigma2(1.);

  auto  p = helper_functions::eon_parameters_to_gpr(parameters);
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = matter->getCell()(i);
  }

  gprfunc->initialize(p, atoms_config);
  gprfunc->setHyperparameters(all_obs, atoms_config);
  gprfunc->calculatePotential(all_obs);
  gprfunc->optimize(all_obs);
  // gprfunc->evaluateEnergyAndGradient();
  EXPECT_EQ(gprfunc->energy_and_gradient[0], 3)
      << "Energy does not match";

  // double energ = gprfunc->evaluateEnergy(gprfunc->R_matrix, gprfunc->R_indices, matter->getPositionsV());
  // gprfunc->calculateMeanPrediction();
  // gprfunc->calculatePosteriorMeanPrediction();
  // EXPECT_EQ(energ, 3)
  //     << "Covariance does not match";
  // EXPECT_EQ(gprfunc->energy_and_gradient[0], 3)
  //     << "Potential energy does not match";
  delete matter;
  delete parameters;
  delete gprfunc;
}

} /* namespace tests */
