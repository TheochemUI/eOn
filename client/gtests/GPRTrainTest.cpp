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

GPRTrainTest::~GPRTrainTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRTrainTest, TestMatter) {
  aux::ProblemSetUp problem_setup;
  auto config_data = helper_functions::eon_matter_to_frozen_conf_info(this->initmatter.get(),  0.5);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  auto R_init = std::get<gpr::Coord>(config_data);
  // Setup the observations
  auto initPath = helper_functions::prepInitialPath(this->parameters.get());
  auto imgArray = std::get<std::vector<Matter> >(initPath);
  auto tangentArray = std::get<std::vector<AtomMatrix> >(initPath);
  auto projForceArray = tangentArray; // Initially the same
  auto obspath = helper_functions::prepInitialObs(imgArray);
  // Setup GPR
  // GPRPotential pot{parameters};
  gpr::GPRSetup gpr_parameters;
  aux::AuxiliaryFunctionality aux_func;
  gpr::GaussianProcessRegression *gprfunc = new gpr::GaussianProcessRegression();

  gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gpr_parameters.jitter_sigma2 = 0.;
  gprfunc->setParameters(gpr_parameters);

  gprfunc->getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
  gprfunc->getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
  gprfunc->getSexpAtCovarianceFunction()->setConfInfo(atoms_config);

  gprfunc->getConstantCovarianceFunction()->setConstSigma2(1.);

  auto  p = helper_functions::eon_parameters_to_gpr(this->parameters.get());
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = this->initmatter->getCell()(i);
  }

  gprfunc->initialize(p, atoms_config);
  gprfunc->setHyperparameters(obspath, atoms_config);
  gprfunc->calculatePotential(obspath);
  gprfunc->optimize(obspath);
  // gprfunc->calculatePosteriorMeanPrediction();
  GPRPotential pot{this->parameters.get()};
  pot.registerGPRObject(gprfunc);
  auto egf_gpro = helper_functions::energy_and_forces(this->initmatter.get(), &pot);
  auto energy_gpro = std::get<double>(egf_gpro);
  auto forces_gpro = std::get<AtomMatrix>(egf_gpro);
  // // gprfunc->evaluateEnergyAndGradient();
  // EXPECT_EQ(energy_gpro, 3)
  //     << "Energy does not match";

  // double energ = gprfunc->evaluateEnergy(gprfunc->R_matrix, gprfunc->R_indices, matter->getPositionsV());
  // gprfunc->calculateMeanPrediction();
  // gprfunc->calculatePosteriorMeanPrediction();
  // EXPECT_EQ(energ, 3)
  //     << "Covariance does not match";
  // EXPECT_EQ(gprfunc->energy_and_gradient[0], 3)
  //     << "Potential energy does not match";
  delete gprfunc;
}

} /* namespace tests */
