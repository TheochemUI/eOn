/*
 * GPRPotTest.cpp
 *
 *  Created on: 27 Mar 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "GPRPotTest.h"

namespace tests {

GPRPotTest::GPRPotTest() {
  this->reactantFilename = helper_functions::getRelevantFile("reactant.con");
  this->productFilename = helper_functions::getRelevantFile("product.con");

  this->parameters = std::make_unique<Parameters>();
  this->parameters->potential = "morse_pt";
  this->parameters->nebImages = 7;
  this->parameters->LogPotential = false;
  log_init(parameters.get(), (char *)"test.log");

  this->initmatter = std::make_unique<Matter>(parameters.get());
  this->finalmatter = std::make_unique<Matter>(parameters.get());

  this->initmatter->con2matter(reactantFilename);
  this->finalmatter->con2matter(productFilename);

  this->gprfunc = std::make_unique<gpr::GaussianProcessRegression>();

  // Constants
  this->init_eref = this->initmatter.get()->getPotentialEnergy();
  this->init_frcsref = this->initmatter.get()->getForcesFree();
  this->config_data = helper_functions::eon_matter_to_frozen_conf_info(
      this->initmatter.get(), 5);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  auto R_init = std::get<gpr::Coord>(config_data);
  // Setup the observations
  this->imgArray = helper_functions::prepInitialPath(this->parameters.get());
  this->obspath = helper_functions::prepInitialObs(imgArray);
  // Setup GPR
  auto eondat = std::make_pair(*this->parameters, *this->initmatter);
  *this->gprfunc = helper_functions::initializeGPR(*this->gprfunc, atoms_config,
                                                   obspath, eondat);
  this->gprfunc->setHyperparameters(obspath, atoms_config);
  this->gprfunc->optimize(obspath);
  // Comparer
  this->comparer = [&](const gpr::EigenMatrix &lhs,
                       const gpr::EigenMatrix &rhs) -> bool {
    return lhs.isApprox(rhs, 1e-2); // FIXME: lowered because of the gradients
  };
}

GPRPotTest::~GPRPotTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRPotTest, TestMatter) {
  // Matter calls
  GPRPotential gprpot{this->parameters.get()};
  gprpot.registerGPRObject(this->gprfunc.get());
  auto matterClone = std::make_unique<Matter>(this->parameters.get());
  matterClone->con2matter(this->reactantFilename);
  matterClone->setPotential(&gprpot);
  ASSERT_NEAR(matterClone->getPotentialEnergy(), init_eref, 1e-3)
      << "Energy does not match";
  ASSERT_PRED2(comparer, matterClone->getForcesFree(),
               init_frcsref)
      << "Forces don't match";
}

} /* namespace tests */
