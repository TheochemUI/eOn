/*
 * GPRHelpersTest.cpp
 *
 *  Created on: 07 Feb 2021
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>
#include <numeric>
#include <memory>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/LJ/LJ.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Job.h"
#include "../Log.h"
#include "../Matter.h"
#include "../MinModeSaddleSearch.h"
#include "../Parameters.h"
#include "GPRHelpersTest.h"

namespace tests {

GPRHelpersTest::GPRHelpersTest() {
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

GPRHelpersTest::~GPRHelpersTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRHelpersTest, TestMatter) {
  // TODO: Clean up
  string confile("pos.con");
  double energy_morse_ref{-1775.79}; // matter->getPotentialEnergy()
  Parameters *pmorse = new Parameters;
  pmorse->potential = "morse_pt";
  pmorse->LogPotential = false;
  Matter *pm = new Matter(pmorse);
  pm->con2matter(confile);
  Morse pot_morse;
  pot_morse.setParams(pmorse); // *mostly* pointless
  EXPECT_EQ(pot_morse.getName(), "morse_pt"s)
    << "Potential name is incorrect [Morse]";
  double energy_morse{0};
  AtomMatrix forces_morse = AtomMatrix::Constant(pm->numberOfAtoms(), 3, 0);
  auto egf_morse = helper_functions::energy_and_forces(pm, &pot_morse);
  energy_morse = std::get<double>(egf_morse);
  forces_morse = std::get<AtomMatrix>(egf_morse);
  EXPECT_EQ(forces_morse, pm->getForces())
      << "Forces do not match [Morse]";
  EXPECT_EQ(energy_morse, pm->getPotentialEnergy())
      << "Potential energy does not match [Morse]";
  Parameters *plj = new Parameters;
  plj->potential = "lj";
  plj->LogPotential = false;
  Matter *pljm = new Matter(plj);
  pljm->con2matter(confile);
  LJ ljpot;
  ljpot.setParams(plj);
  EXPECT_EQ(ljpot.getName(), "lj"s)
    << "Potential name is incorrect [LJ]";
  double energy_lj{0};
  AtomMatrix forces_lj = AtomMatrix::Constant(pljm->numberOfAtoms(), 3, 0);
  auto egf_lj = helper_functions::energy_and_forces(pm, &ljpot);
  energy_lj = std::get<double>(egf_lj);
  forces_lj = std::get<AtomMatrix>(egf_lj);
  pljm->setPotential(&ljpot);
  EXPECT_EQ(pljm->getPotential()->getName(), "lj"s)
    << "Potential name is incorrect in matter [LJ]";
  EXPECT_EQ(forces_lj, pljm->getForces())
      << "Forces do not match [LJ]";
  EXPECT_EQ(energy_lj, pljm->getPotentialEnergy())
      << "Potential energy does not match [LJ]";
  delete pljm;
  delete plj;
  delete pm;
  delete pmorse;
}

TEST_F(GPRHelpersTest, TestNEBInitPath) {
  // Setup the run
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  EXPECT_EQ(imgArray.back().getForces(), this->finalmatter->getForces())
      << "Forces do not match";
}

TEST_F(GPRHelpersTest, TestNEBInitObs) {
  // Setup the observations
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  auto obspath = helper_functions::prepInitialObs(imgArray);
  // Assertions
  const int nimages = this->parameters->nebImages;
  const int totImages = nimages + 2; // Final and end
  EXPECT_EQ(obspath.E.getSize(), totImages)
      << "Number of elements of energies do not match";
  EXPECT_EQ(obspath.R.getNumPoints(), totImages)
      << "Number of elements of R matrices do not match";
  EXPECT_EQ(obspath.G.getNumPoints(), totImages)
      << "Number of elements of G matrices do not match";
  // obspath.printSizes();
}

TEST_F(GPRHelpersTest, TestPathLength) {
  // Setup the observations
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  double plen = helper_functions::get_path_length(imgArray);
  double true_plen = imgArray.front().distanceTo(imgArray.back());
  EXPECT_NEAR(plen, true_plen, this->threshold)
      << "Path length is incorrect";
}

TEST_F(GPRHelpersTest, TestUnravelFreeV) {
  // Setup the observations
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  auto ufreecoords = helper_functions::unravel_free_coordsV(imgArray);
  const size_t nfree = imgArray.front().numberOfFreeAtoms();
  const size_t nimgs = imgArray.size();
  EXPECT_EQ(ufreecoords.size(), nfree*nimgs*3)
    << "Unraveled coordinates have the wrong number of elements!\n";
  for (size_t idx{0}; idx < nimgs; idx++){
    EXPECT_EQ(ufreecoords.segment(idx * nfree * 3, nfree * 3), imgArray[idx].getPositionsFreeV())
      << "Point isn't the same";
  }
}

TEST_F(GPRHelpersTest, TestUnravelV) {
  // Setup the observations
  auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
  auto ucoords = helper_functions::unravel_coordsV(imgArray);
  const size_t npoints = imgArray.front().numberOfAtoms();
  const size_t nimgs = imgArray.size();
  EXPECT_EQ(ucoords.size(), npoints*nimgs*3)
    << "Unraveled coordinates have the wrong number of elements!\n";
  for (size_t idx{0}; idx < nimgs; idx++){
    EXPECT_EQ(ucoords.segment(idx * npoints * 3, npoints * 3), imgArray[idx].getPositionsV())
      << "Point isn't the same";
  }
}

// TEST_F(GPRHelpersTest, TestUnravelFree) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ufreecoords = helper_functions::unravel_free_coords(imgArray);
//   const size_t nfree = imgArray.front().numberOfFreeAtoms();
//   const size_t nimgs = imgArray.size();
//   size_t nimgrows{0};
//   for (auto img:imgArray){
//     nimgrows += img.numberOfFreeAtoms();
//   }
//   // IC(ufreecoords.rows(), ufreecoords.cols(), ufreecoords.size());
//   EXPECT_EQ(ufreecoords.rows(), nimgrows)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   EXPECT_EQ(ufreecoords.cols(), nfree*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ufreecoords.row(idx), imgArray[idx].getPositionsFree())
//       << "Point isn't the same";
//     // IC(ufreecoords.row(idx), imgArray[idx].getPositionsFree());
//   }
// }

// TEST_F(GPRHelpersTest, TestUnravel) {
//   // Setup the observations
//   auto imgArray = helper_functions::prepInitialPath(this->parameters.get());
//   auto ucoords = helper_functions::unravel_coords(imgArray);
//   const size_t npositions = imgArray.front().numberOfAtoms();
//   const size_t nimgs = imgArray.size();
//   size_t nimgrows{0};
//   // for (auto img:imgArray){
//   //   nimgrows += img.numberOfAtoms();
//   // }
//   IC(ucoords.rows(), ucoords.cols(), ucoords.size());
//   // EXPECT_EQ(ucoords.rows(), nimgrows)
//   //   << "Unraveled coordinates have the wrong number of elements!\n";
//   EXPECT_EQ(ucoords.cols(), npositions*3)
//     << "Unraveled coordinates have the wrong number of elements!\n";
//   for (size_t idx{0}; idx < nimgs; idx++){
//     EXPECT_EQ(ucoords.row(idx), imgArray[idx].getPositions())
//       << "Point isn't the same";
//     // IC(ucoords.row(idx), imgArray[idx].getPositions());
//   }
// }

} /* namespace tests */
