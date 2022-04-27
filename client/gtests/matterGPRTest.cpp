
/*
 * matterGPRTest.cpp
 *
 *  Created on: 15 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>
#include <memory>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/LJ/LJ.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Job.h"
#include "../Log.h"
#include "../GPRMatter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "matterGPRTest.h"

namespace tests {

matterGPRTest::matterGPRTest() {
    reactantFilename = helper_functions::getRelevantFile("reactant.con");
    productFilename = helper_functions::getRelevantFile("product.con");

    params.potential = "morse_pt";
    params.nebImages = 7;
    params.LogPotential = false;
    log_init(&params, (char *)"test.log");
}

matterGPRTest::~matterGPRTest() {
  // TODO Auto-generated destructor stub
}

    TEST_F(matterGPRTest, trainGPRobj){
    // Setup the Matter objects
    Matter reactant{&params}, product{&params}, testp{&params};
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    testp.setPositions((reactant.getPositions()*10.01) / 2);
    // Setup the observations
    auto imgArray = helper_functions::prepInitialPath(&params);
    // Setup GPR
    GPRobj gprob(reactant, params);
    gprob.trainGPR(imgArray);
    // gprob.curpath.printSizes();
    EXPECT_EQ(gprob.curpath.R.getSize(), 27) << "Position vector has wrong size.";
    EXPECT_EQ(gprob.curpath.G.getSize(), 27) << "Gradient vector has wrong size.";
    EXPECT_EQ(gprob.curpath.E.getSize(), 9) << "Energy vector has wrong size.";
    EXPECT_EQ(gprob.curpath.R.getNumRows(), 9) << "Position vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.E.getNumRows(), 9) << "Energy has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.G.getNumRows(), 9) << "Gradient vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.R.getNumCols(), 3) << "Position vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.E.getNumCols(), 1) << "Energy has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.G.getNumCols(), 3) << "Gradient vector has wrong number of elements.";
    }

    TEST_F(matterGPRTest, retrainGPRobj){
    // Setup the Matter objects
    Matter reactant{&params}, product{&params}, testp{&params};
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    testp.con2matter(productFilename);
    testp.setPositions((reactant.getPositions()*10.01) / 2);
    // Setup the observations
    auto imgArray = helper_functions::prepInitialPath(&params);
    // Setup GPR
    GPRobj gprob(reactant, params);
    gprob.trainGPR(imgArray);
    std::vector<Matter> ia2{testp};
    gprob.retrainGPR(ia2);
    EXPECT_EQ(gprob.curpath.R.getSize(), 30) << "Position vector has wrong size.";
    EXPECT_EQ(gprob.curpath.G.getSize(), 30) << "Gradient vector has wrong size.";
    EXPECT_EQ(gprob.curpath.E.getSize(), 10) << "Energy vector has wrong size.";
    EXPECT_EQ(gprob.curpath.R.getNumRows(), 10) << "Position vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.E.getNumRows(), 10) << "Energy has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.G.getNumRows(), 10) << "Gradient vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.R.getNumCols(), 3) << "Position vector has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.E.getNumCols(), 1) << "Energy has wrong number of elements.";
    EXPECT_EQ(gprob.curpath.G.getNumCols(), 3) << "Gradient vector has wrong number of elements.";
    }

    TEST_F(matterGPRTest, makeGPRMatter){
    // Setup the Matter objects
    Matter reactant{&params}, product{&params}, testp{&params};
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    testp.con2matter(productFilename);
    testp.setPositions((reactant.getPositions()*10.01) / 2);
    // Setup the observations
    auto imgArray = helper_functions::prepInitialPath(&params);
    // Setup GPR
    auto gpf = std::make_shared<GPRobj>(reactant, params);
    gpf->trainGPR(imgArray);
    // Setup GPRMatter
    GPRMatter gp_testp{testp, gpf}, gp_reactant{reactant, gpf}, gp_prod{product, gpf};
    auto [energy, forces] = gp_reactant.true_free_energy_forces();
    AtomMatrix trueForces = reactant.getForcesFree();
    size_t nfree(reactant.numberOfFreeAtoms());
    // Test true energies
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"True reactant energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True reactant forces don't match";
    trueForces = product.getForcesFree();
    std::tie(energy, forces) = gp_prod.true_free_energy_forces();
    EXPECT_NEAR(energy, product.getPotentialEnergy(), this->threshold)<<"True product energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True product forces don't match";
    std::tie(energy, forces) = gp_testp.true_free_energy_forces();
    trueForces = testp.getForcesFree();
    ASSERT_NEAR(energy, testp.getPotentialEnergy(), this->threshold)<<"True testp energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True testp forces don't match";
    std::tie(energy, forces) = gp_reactant.gpr_energy_forces();
    trueForces = reactant.getForcesFree();
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"GPR reactant energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR reactant forces don't match";
    trueForces = product.getForcesFree();
    std::tie(energy, forces) = gp_prod.gpr_energy_forces();
    EXPECT_NEAR(energy, product.getPotentialEnergy(), this->threshold)<<"GPR product energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR product forces don't match";
    trueForces = testp.getForcesFree();
    std::tie(energy, forces) = gp_testp.gpr_energy_forces();
    EXPECT_NE(energy, testp.getPotentialEnergy())<<"GPR testp energies should not match";
    EXPECT_FALSE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"GPR testp forces should not match";
    gp_testp.updateMatter(reactant);
    trueForces = reactant.getForcesFree();
    std::tie(energy, forces) = gp_testp.gpr_energy_forces();
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"GPR reactant energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR reactant forces don't match";
    std::tie(energy, forces) = gp_testp.true_free_energy_forces();
    trueForces = reactant.getForcesFree();
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"GPR reactant energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR reactant forces don't match";
    }

    TEST_F(matterGPRTest, testTruePotMatter){
    // Setup the Matter objects
    Matter reactant{&params}, product{&params}, testp{&params};
    reactant.con2matter(reactantFilename);
    product.con2matter(productFilename);
    testp.con2matter(productFilename);
    testp.setPositions((reactant.getPositions()*10.01) / 2);
    // Setup the observations
    auto imgArray = helper_functions::prepInitialPath(&params);
    // Setup GPR
    auto gpf = std::make_shared<GPRobj>(reactant, params);
    gpf->trainGPR(imgArray);
    // Setup GPRMatter
    GPRMatter gp_testp{testp, gpf}, gp_reactant{reactant, gpf}, gp_prod{product, gpf};
    AtomMatrix trueForces = reactant.getForcesFree();
    size_t nfree(reactant.numberOfFreeAtoms());
    EXPECT_NEAR(gp_reactant.truePotMatter.getPotentialEnergy(), reactant.getPotentialEnergy(), this->threshold)<<"GPRMatter internal object energies don't match";
    EXPECT_NEAR((gp_reactant.truePotMatter.getForcesFree()-reactant.getForcesFree()).norm(), 0, this->threshold)<<"GPRMatter internal object forces don't match true forces";
    EXPECT_NEAR(gp_prod.truePotMatter.getPotentialEnergy(), product.getPotentialEnergy(), this->threshold)<<"GPRMatter internal object energies don't match";
    EXPECT_NEAR((gp_prod.truePotMatter.getForcesFree()-product.getForcesFree()).norm(), 0, this->threshold)<<"GPRMatter internal object forces don't match true forces";
    EXPECT_NEAR(gp_testp.truePotMatter.getPotentialEnergy(), testp.getPotentialEnergy(), this->threshold)<<"GPRMatter internal object energies don't match";
    EXPECT_NEAR((gp_testp.truePotMatter.getForcesFree()-testp.getForcesFree()).norm(), 0, this->threshold)<<"GPRMatter internal object forces don't match true forces";
    }

    // TEST_F(matterGPRTest, insertIntoMatter){
    // // Setup the Matter objects
    // Matter reactant{&params}, product{&params}, testp{&params};
    // reactant.con2matter(reactantFilename);
    // product.con2matter(productFilename);
    // testp.con2matter(productFilename);
    // testp.setPositions((reactant.getPositions()*10.01) / 2);
    // // Setup the observations
    // auto imgArray = helper_functions::prepInitialPath(&params);
    // // Setup GPR
    // auto gpf = std::make_shared<GPRobj>(reactant, params);
    // gpf->trainGPR(imgArray);
    // // Setup GPRMatter
    // GPRMatter gp_testp{testp, gpf}, gp_reactant{reactant, gpf}, gp_prod{product, gpf};
    // auto [energy, forces] = gp_reactant.true_free_energy_forces();
    // auto gppot = gpf->yieldGPRPot();
    // reactant.setPotential(&gppot);
    // AtomMatrix trueForces = reactant.getForcesFree();
    // size_t nfree(reactant.numberOfFreeAtoms());
    // EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"GPR reactant energies don't match";
    // EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR reactant forces don't match";
    // }

} /* namespace tests */
