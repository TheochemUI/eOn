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

// TEST_F(NEBJobGPRTest, TestMatter) {
//   const auto init_eref = this->initmatter.get()->getPotentialEnergy(); // Surprisingly still needed, buggy?
//   const auto init_frcsref = this->initmatter.get()->getForcesFree();
//   auto config_data = helper_functions::eon_matter_to_frozen_conf_info(this->initmatter.get(),  2);
//   auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
//   Parameters eonp = *this->parameters;
//   // Setup the run
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto obspath = helper_functions::prepInitialObs(imgArray);
//   // Setup GPR
//   auto eondat = std::make_pair(eonp, *this->initmatter);
//   *this->gprfunc = helper_functions::initializeGPR(*this->gprfunc, atoms_config, obspath, eondat);
//   this->gprfunc->setHyperparameters(obspath, atoms_config);
//   this->gprfunc->optimize(obspath);
//   // Prepare GPR potential
//   GPRPotential gprpot{this->gprparameon.get()};
//   gprpot.registerGPRObject(this->gprfunc.get());
//   // Setup matter objects
//   auto matterOne = std::make_unique<Matter>(this->gprparameon.get());
//   matterOne->con2matter(this->reactantFilename);
//   matterOne->setPotential(&gprpot);
//   auto matterFin = std::make_unique<Matter>(this->gprparameon.get());
//   matterFin->con2matter(this->productFilename);
//   matterFin->setPotential(&gprpot);
//   auto matterTest = std::make_unique<Matter>(this->gprparameon.get());
//   matterTest->con2matter(this->productFilename);
//   matterTest->setPositions((matterOne->getPositions()*1.01+matterFin->getPositions()*0.03) / 2);
//   matterTest->setPotential(&gprpot);
//   double blah = matterTest->getPotentialEnergy();
//   std::cout<<matterTest->getPotentialEnergy()<<" Matter at Test point\n";
//   // RUN NEB!!!
//   NudgedElasticBand *neb = new NudgedElasticBand(matterOne.get(), matterFin.get(), this->gprparameon.get());
//   neb->compute();
//   bool mustUpdate = helper_functions::maybeUpdateObs(*neb, obspath, eonp);
//   this->gprfunc->setHyperparameters(obspath, atoms_config, false);
//   this->gprfunc->optimize(obspath);
//   gprpot.registerGPRObject(this->gprfunc.get());
//   matterTest->setPotential(&gprpot);
//   ASSERT_NE(blah, matterTest->getPotentialEnergy())<<"Energy not changed after updating gpr";
//   while(mustUpdate){
//     this->gprfunc->setHyperparameters(obspath, atoms_config, false);
//     this->gprfunc->optimize(obspath);
//     auto nebTwo = helper_functions::prepGPRNEBround(*this->gprfunc,
//                                                 *matterOne, *matterFin,
//                                                 *this->gprparameon);
//     nebTwo->compute();
//     mustUpdate = helper_functions::maybeUpdateObs(*nebTwo, obspath, eonp);
//   };
//   // Final round
//     auto nebFin = helper_functions::prepGPRNEBround(*this->gprfunc,
//                                                 *matterOne, *matterFin,
//                                                 *this->gprparameon);
//     nebFin->compute();
//     delete neb;
// }

} /* namespace tests */
