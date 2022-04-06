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

    gprfunc = std::make_unique<gpr::GaussianProcessRegression>();
}

GPRTrainTest::~GPRTrainTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRTrainTest, TestMatter) {
  // Constants
  const auto init_eref = this->initmatter.get()->getPotentialEnergy();
  const auto init_frcsref = this->initmatter.get()->getForcesFree();
  aux::ProblemSetUp problem_setup;
  auto config_data = helper_functions::eon_matter_to_frozen_conf_info(this->initmatter.get(),  5);
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
  this->gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  this->gprfunc->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gpr_parameters.jitter_sigma2 = 0.;
  this->gprfunc->setParameters(gpr_parameters);

  this->gprfunc->getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
  this->gprfunc->getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
  this->gprfunc->getSexpAtCovarianceFunction()->setConfInfo(atoms_config);

  this->gprfunc->getConstantCovarianceFunction()->setConstSigma2(1.);

  auto  p = helper_functions::eon_parameters_to_gprpot(this->parameters.get());
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = this->initmatter->getCell()(i);
  }

  this->gprfunc->initialize(p, atoms_config);
  this->gprfunc->setHyperparameters(obspath, atoms_config);
  this->gprfunc->optimize(obspath);
  // Multiple observations
  // auto oo = obspath;
  // this->gprfunc->calculatePotential(oo);
  // std::cout<<"Energy: "<<oo.E.extractEigenMatrix()<<std::endl;
  // std::cout<<"Correct Energy: "<<obspath.E.extractEigenMatrix()<<std::endl;

  // Single observation
  gpr::Observation o;
  o.clear();
  o.R.resize(init_frcsref.rows(), init_frcsref.cols()); // Also works with getPositionsFree()
  o.G.resize(init_frcsref.rows(), init_frcsref.cols());
  o.R.assignFromEigenMatrix(this->initmatter.get()->getPositionsFree());
  o.E.resize(1);
  this->gprfunc->calculatePotential(o);
  ASSERT_NEAR(o.E.extractEigenMatrix()(0), init_eref, this->threshold*1e3)
      << "Energy does not match";
  EXPECT_TRUE((o.G.extractEigenMatrix() * -1).isApprox(init_frcsref, this->threshold))
      << "Forces do not match";
  // std::cout<<"Energy: "<<o.E.extractEigenMatrix()<<std::endl;
  // std::cout<<"Correct Energy: "<<this->initmatter.get()->getPotentialEnergy()<<std::endl;
  // std::cout<<"Forces : \n"<<o.G.extractEigenMatrix()*-1<<std::endl;
  // std::cout<<"Correct Forces: \n"<<init_frcsref<<std::endl;

  // Function calls
  GPRPotential pot{this->parameters.get()};
  pot.registerGPRObject(this->gprfunc.get());
  auto egf_gpro = helper_functions::gpr_energy_and_forces(this->initmatter.get(), &pot);
  auto energy_gpro = std::get<double>(egf_gpro);
  auto forces_gpro = std::get<AtomMatrix>(egf_gpro);
  ASSERT_NEAR(energy_gpro, init_eref, this->threshold*1e3)
      << "Energy does not match";
  EXPECT_TRUE(forces_gpro.isApprox(init_frcsref, this->threshold))
      << "Forces do not match";
}

} /* namespace tests */
