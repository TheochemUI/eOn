/*
 * NEBJobGPRTest.cpp
 *
 *  Created on: 6 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "NEBJobGPRTest.h"

namespace tests {

NEBJobGPRTest::NEBJobGPRTest() {
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

NEBJobGPRTest::~NEBJobGPRTest() {}

TEST_F(NEBJobGPRTest, TestMatter) {
  // Prepare GPR potential
  auto gprpot = new GPRPotential {this->parameters.get()};
  gprpot->registerGPRObject(this->gprfunc.get());
  // // Setup matter objects
  auto matterOne = std::make_unique<Matter>(this->parameters.get());
  matterOne->con2matter(this->reactantFilename);
  matterOne->setPotential(gprpot);
  auto matterFin = std::make_unique<Matter>(this->parameters.get());
  matterFin->con2matter(this->productFilename);
  matterFin->setPotential(gprpot);
  auto matterTest = std::make_unique<Matter>(this->parameters.get());
  matterTest->con2matter(this->productFilename);
  matterTest->setPositions(
      (matterOne->getPositions() * 1.01 + matterFin->getPositions() * 0.03) /
      2);
  matterTest->setPotential(gprpot);
  double blah = matterTest->getPotentialEnergy();
  // RUN NEB!!!
  NudgedElasticBand *neb = new NudgedElasticBand(
      matterOne.get(), matterFin.get(), this->parameters.get());
  neb->compute();
  for (long idx {0}; idx <= neb->images+1; idx++){// excludes final, initial
    neb->image[idx]->setPotential(gprpot);
    // IC(idx, neb->image[idx]->getPotential()->getName());
  }
  obspath.printSizes();
  bool mustUpdate = helper_functions::maybeUpdateObs(*neb, obspath, *this->parameters);
  // obspath.printSizes();
  // this->gprfunc->setHyperparameters(obspath, atoms_conf, false);
  this->gprfunc->optimize(obspath);
  gprpot->registerGPRObject(this->gprfunc.get());
  matterTest->setPotential(gprpot);
  IC(mustUpdate, blah, matterTest->getPotentialEnergy());
  ASSERT_NE(blah, matterTest->getPotentialEnergy())
      << "Energy not changed after updating gpr";
  while(mustUpdate){
    // this->gprfunc->setHyperparameters(obspath, atoms_conf, false);
    this->gprfunc->optimize(obspath);
    auto nebTwo = helper_functions::prepGPRNEBround(*gprpot,
                                                *matterOne, *matterFin,
                                                *this->parameters);
    nebTwo->compute();
    mustUpdate = helper_functions::maybeUpdateObs(*nebTwo, obspath, *this->parameters);
    gprpot->registerGPRObject(this->gprfunc.get());
    matterTest->setPotential(gprpot);
    IC(mustUpdate, blah, matterTest->getPotentialEnergy());
  };
  // // Final round
  //   auto nebFin = helper_functions::prepGPRNEBround(*this->gprfunc,
  //                                               *matterOne, *matterFin,
  //                                               *this->parameters);
  //   nebFin->compute();
  delete neb;
}

} /* namespace tests */
