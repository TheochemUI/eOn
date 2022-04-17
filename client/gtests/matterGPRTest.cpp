
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
    auto obspath = helper_functions::prepInitialObs(imgArray);
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
    auto obspath = helper_functions::prepInitialObs(imgArray);
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
    auto obspath = helper_functions::prepInitialObs(imgArray);
    // Setup GPR
    auto gpf = std::make_shared<GPRobj>(reactant, params);
    gpf->trainGPR(imgArray);
    // Setup GPRMatter
    GPRMatter gp_testp{testp, gpf}, gp_reactant{reactant, gpf}, gp_prod{product, gpf};
    auto reactant_pe_forces = gp_reactant.true_free_energy_forces();
    double energy = std::get<double>(reactant_pe_forces);
    AtomMatrix forces = std::get<AtomMatrix>(reactant_pe_forces);
    AtomMatrix trueForces = reactant.getForcesFree();
    size_t nfree(reactant.numberOfFreeAtoms());
    // Test true energies
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"True reactant energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True reactant forces don't match";
    trueForces = product.getForcesFree();
    auto product_pe_forces = gp_prod.true_free_energy_forces();
    forces = std::get<AtomMatrix>(product_pe_forces);
    energy = std::get<double>(product_pe_forces);
    EXPECT_NEAR(energy, product.getPotentialEnergy(), this->threshold)<<"True product energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True product forces don't match";
    auto testp_pe_forces = gp_testp.true_free_energy_forces();
    trueForces = testp.getForcesFree();
    forces = std::get<AtomMatrix>(testp_pe_forces);
    energy = std::get<double>(testp_pe_forces);
    ASSERT_NEAR(energy, testp.getPotentialEnergy(), this->threshold)<<"True testp energies don't match";
    EXPECT_TRUE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"True testp forces don't match";
    reactant_pe_forces = gp_reactant.gpr_energy_forces();
    trueForces = reactant.getForcesFree();
    forces = std::get<AtomMatrix>(reactant_pe_forces);
    energy = std::get<double>(reactant_pe_forces);
    EXPECT_NEAR(energy, reactant.getPotentialEnergy(), this->threshold)<<"GPR reactant energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR reactant forces don't match";
    trueForces = product.getForcesFree();
    product_pe_forces = gp_prod.gpr_energy_forces();
    forces = std::get<AtomMatrix>(product_pe_forces);
    energy = std::get<double>(product_pe_forces);
    EXPECT_NEAR(energy, product.getPotentialEnergy(), this->threshold)<<"GPR product energies don't match";
    EXPECT_NEAR((trueForces-forces).norm(), 0, this->threshold)<<"GPR product forces don't match";
    trueForces = testp.getForcesFree();
    testp_pe_forces = gp_testp.gpr_energy_forces();
    forces = std::get<AtomMatrix>(testp_pe_forces);
    energy = std::get<double>(testp_pe_forces);
    EXPECT_NE(energy, testp.getPotentialEnergy())<<"GPR testp energies should not match";
    EXPECT_FALSE((trueForces-forces).isMuchSmallerThan(this->threshold))<<"GPR testp forces should not match";
    // EXPECT_EQ(energy, reactant.getPotentialEnergy())<<"True energies don't match";
    // std::cout<<energy<<"\n";
    // std::cout<<testp.getPotentialEnergy()<<"\n";
    // std::cout<<"\nAre energies close? "<<gpm.areEnergiesCloseToTrue();
    // std::cout<<"\nAre forces close? "<<gpm.areForcesCloseToTrue();
    }

// TEST_F(matterGPRTest, TestGPRTrain) {

//     Matter initmatter{&params};
//     initmatter.con2matter(reactantFilename);
//     auto config_data = helper_functions::eon_matter_to_frozen_conf_info(&initmatter,  5);
//     auto atoms_config = std::get<gpr::AtomsConfiguration>(config_data);
//     // Setup the observations
//     auto imgArray = helper_functions::prepInitialPath(&params);
//     auto obspath = helper_functions::prepInitialObs(imgArray);
//     // Setup GPR
//     auto eondat =  std::make_pair(params,initmatter);
//     Morse pot_morse;
//     GPRobj gprob(initmatter, params);
//     gprob.trainGPR(imgArray, &pot_morse);
//     // gprob.curpath.printSizes();
//     Matter matterTest{&params};
//     matterTest.con2matter(productFilename);
//     matterTest.setPositions((initmatter.getPositions()*10.01) / 2);
//     std::vector<Matter> ia2{initmatter, matterTest};
//     gprob.retrainGPR(ia2, &pot_morse);
//     // imgArray.push_back(initmatter);
//     // imgArray.push_back(matterTest);
//     // gprob.trainGPR(imgArray);
//     // gprob.curpath.printSizes();
//     // Matter matterFin{&params};
//     // matterFin.con2matter(productFilename);
//     // std::vector<Matter> ia3{matterFin};
//     // // gprob.trainGPR(ia3);
//     // gprob.curpath.printSizes();
//     // ptrGPRM->trainGPR(obspath);
//     // auto pe_forces = ptrGPRM->gpr_energy_forces(initmatter);
//     // double energy = std::get<double>(pe_forces);
//     // AtomMatrix gradEig = std::get<AtomMatrix>(pe_forces);
//     // std::cout<<"\nGPR Energy:\n";
//     // std::cout<<energy;
//     // std::cout<<"\nGPR Gradients:\n";
//     // std::cout<<gradEig;
//     // auto pe_true_forces = ptrGPRM->true_free_energy_forces(initmatter);
//     // double true_energy = std::get<double>(pe_true_forces);
//     // AtomMatrix true_gradEig = std::get<AtomMatrix>(pe_true_forces);
//     // std::cout<<"\nTrue Energy:\n";
//     // std::cout<<true_energy;
//     // std::cout<<"\nTrue Gradients:\n";
//     // std::cout<<true_gradEig;
//     // auto matterTest = std::make_unique<Matter>(&params);
//     // matterTest->con2matter(productFilename);
//     // matterTest->setPositions((initmatter.getPositions()*10.01) / 2);
//     // pe_forces = ptrGPRM->gpr_energy_forces(*matterTest);
//     // energy = std::get<double>(pe_forces);
//     // gradEig = std::get<AtomMatrix>(pe_forces);
//     // std::cout<<"\nGPR Energy:\n";
//     // std::cout<<energy;
//     // std::cout<<"\nGPR Gradients:\n";
//     // std::cout<<gradEig;
//     // pe_true_forces = ptrGPRM->true_free_energy_forces(*matterTest);
//     // true_energy = std::get<double>(pe_true_forces);
//     // true_gradEig = std::get<AtomMatrix>(pe_true_forces);
//     // std::cout<<"\nTrue Energy:\n";
//     // std::cout<<true_energy;
//     // std::cout<<"\nTrue Gradients:\n";
//     // std::cout<<true_gradEig;
//     // std::cout<<"\nTrue Energy:\n";
//     // std::cout<<matterTest->getPotentialEnergy(); //true_energy;
//     // std::cout<<"\nTrue Gradients:\n";
//     // std::cout<<matterTest->getForcesFree();
//     // std::cout<<ptrGPRM->isCloseTo(*matterTest);
//     // std::cout<<ptrGPRM->isCloseTo(initmatter);
//     // obspath.append(helper_functions::eon_matter_to_init_obs(*matterTest));
//     // ptrGPRM->trainGPR(obspath);
//     // pe_forces = ptrGPRM->gpr_energy_forces(*matterTest);
//     // energy = std::get<double>(pe_forces);
//     // gradEig = std::get<AtomMatrix>(pe_forces);
//     // std::cout<<"\nGPR Energy:\n";
//     // std::cout<<energy;
//     // std::cout<<"\nGPR Gradients:\n";
//     // std::cout<<gradEig;
// }

} /* namespace tests */
