/*
 * NEBJobGPRTest.cpp
 *
 *  Created on: 6 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../NudgedElasticBandJob.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Log.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "NEBJobGPRTest.h"

namespace tests {

NEBJobGPRTest::NEBJobGPRTest() {
    reactantFilename = helper_functions::getRelevantFile("reactant.con");
    productFilename = helper_functions::getRelevantFile("product.con");

    parameters = std::make_unique<Parameters>();
    parameters->potential = "morse_pt";
    parameters->nebImages = 7;
    parameters->LogPotential = false;
    log_init(parameters.get(), (char *)"test.log");

    gprparameon = std::make_unique<Parameters>();
    gprparameon->potential = "gpr_pot";
    gprparameon->nebImages = 7;
    gprparameon->LogPotential = false;

    initmatter = std::make_unique<Matter>(parameters.get());
    finalmatter = std::make_unique<Matter>(parameters.get());

    initmatter->con2matter(reactantFilename);
    finalmatter->con2matter(productFilename);

    gprfunc = std::make_unique<gpr::GaussianProcessRegression>();
}

NEBJobGPRTest::~NEBJobGPRTest() {
}

TEST_F(NEBJobGPRTest, TestMatter) {
  const auto init_eref = this->initmatter.get()->getPotentialEnergy(); // Surprisingly still needed, buggy?
  const auto init_frcsref = this->initmatter.get()->getForcesFree();
  auto config_data = helper_functions::eon_matter_to_frozen_conf_info(this->initmatter.get(),  2);
  auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
  // Setup the run
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  auto obspath = helper_functions::prepInitialObs(imgArray);
  // Setup GPR
  gpr::GPRSetup gpr_parameters;
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
  // Prepare GPR potential
  GPRPotential gprpot{this->gprparameon.get()};
  gprpot.registerGPRObject(this->gprfunc.get());
  // Setup matter objects
  auto matterOne = std::make_unique<Matter>(this->gprparameon.get());
  matterOne->con2matter(this->reactantFilename);
  matterOne->setPotential(&gprpot);
  auto matterFin = std::make_unique<Matter>(this->gprparameon.get());
  matterFin->con2matter(this->productFilename);
  matterFin->setPotential(&gprpot);
  auto matterTest = std::make_unique<Matter>(this->gprparameon.get());
  matterTest->con2matter(this->productFilename);
  matterTest->setPositions((matterOne->getPositions()*1.01+matterFin->getPositions()*0.03) / 2);
  matterTest->setPotential(&gprpot);
  double blah = matterTest->getPotentialEnergy();
  std::cout<<matterTest->getPotentialEnergy()<<" Matter at Test point\n";
  // RUN NEB!!!
  NudgedElasticBand *neb = new NudgedElasticBand(matterOne.get(), matterFin.get(), this->gprparameon.get());
  neb->compute();
  bool mustUpdate = helper_functions::maybeUpdateObs(*neb, obspath);
  this->gprfunc->setHyperparameters(obspath, atoms_config, false);
  this->gprfunc->optimize(obspath);
  gprpot.registerGPRObject(this->gprfunc.get());
  matterTest->setPotential(&gprpot);
  ASSERT_NE(blah, matterTest->getPotentialEnergy())<<"Energy not changed after updating gpr";
  while(mustUpdate){
    this->gprfunc->setHyperparameters(obspath, atoms_config, false);
    this->gprfunc->optimize(obspath);
    gprpot.registerGPRObject(this->gprfunc.get());
    auto matterThree = std::make_unique<Matter>(this->gprparameon.get());
    matterThree->con2matter(this->reactantFilename);
    matterThree->setPotential(&gprpot);
    matterTest->setPotential(&gprpot);
    std::cout<<matterTest->getPotentialEnergy()<<" Matter 3\n";
    ASSERT_NE(blah, matterTest->getPotentialEnergy())<<"Energy not changed after updating gpr";
    auto matterFour = std::make_unique<Matter>(this->gprparameon.get());
    matterFour->con2matter(this->productFilename);
    matterFour->setPotential(&gprpot);
    NudgedElasticBand *nebTwo = new NudgedElasticBand(matterThree.get(), matterFour.get(), this->gprparameon.get());
    nebTwo->compute();
    mustUpdate = helper_functions::maybeUpdateObs(*nebTwo, obspath);
    delete nebTwo;
  };
    gprpot.registerGPRObject(this->gprfunc.get());
    auto matterThree = std::make_unique<Matter>(this->gprparameon.get());
    matterThree->con2matter(this->reactantFilename);
    matterThree->setPotential(&gprpot);
    matterTest->setPotential(&gprpot);
    std::cout<<matterTest->getPotentialEnergy()<<" Matter 3\n";
    ASSERT_NE(blah, matterTest->getPotentialEnergy())<<"Energy not changed after updating gpr";
    auto matterFour = std::make_unique<Matter>(this->gprparameon.get());
    matterFour->con2matter(this->productFilename);
    matterFour->setPotential(&gprpot);
    NudgedElasticBand *nebTwo = new NudgedElasticBand(matterThree.get(), matterFour.get(), this->gprparameon.get());
    nebTwo->compute();
    delete nebTwo;
  delete neb;
}

} /* namespace tests */
