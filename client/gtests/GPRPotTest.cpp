/*
 * GPRPotTest.cpp
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
#include "GPRPotTest.h"

namespace tests {

GPRPotTest::GPRPotTest() {
  // TODO Auto-generated constructor stub
}

GPRPotTest::~GPRPotTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(GPRPotTest, TestMatter) {
  string confile{"pos.con"};
  Parameters *parameters = new Parameters;
  parameters->potential = "morse_pt";
  Matter *matter = new Matter(parameters);
  // Setup GPR
  GPRPotential pot{parameters};
  gpr::GPRSetup gpr_parameters;
  gpr::GaussianProcessRegression gprfunc;
  gprfunc.getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gprfunc.getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
  gpr_parameters.jitter_sigma2 = 0.;
  gprfunc.setParameters(gpr_parameters);
  
  gprfunc.getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
  gprfunc.getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
  gprfunc.getSexpAtCovarianceFunction()->setConfInfo(conf_info);

  gprfunc.getConstantCovarianceFunction()->setConstSigma2(1.);
  // gprfunc.evaluateTrainingCovarianceMatrix(x1, x1_ind, K);

  // Prepare GPR inputs
  Morse trupot;
  int nAtoms = matter->numberOfAtoms();
  auto posdata = matter->getPositions();
  auto celldat = matter->getCell();
  AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  trupot.force(nAtoms, pos, nullptr, frcs, &energy, bx, 1);
  // TODO: Find a less hacky way
  AtomMatrix finForces{forces};
  for (int i = 0; i <nAtoms; i++){
    if(matter->getFixed(i)){
      finForces.row(i).setZero();
    }
  }
  // Indices
  aux::AuxiliaryFunctionality aux_func;
  gpr::Observation observation_all;
  auto R_indices = posdata;
  VectorXd vec_joined(sizeof(double) + finForces.size());
  vec_joined << energy, finForces;
  observation_all.R.assignFromEigenMatrix(posdata);
  // aux_func.assembleMatrixOfRepetitiveCoordinates(
  //   observation_all.R, posdata, R_indices);
  // Setup potential
  // gprfunc.decomposeCovarianceMatrix(posdata, R_indices); // - takes covariance matrix and vector of repetitive indices
  gprfunc.calculateMeanPrediction(vec_joined); // - takes a vector of combined energy and force
  gprfunc.calculatePosteriorMeanPrediction(); // - no arguments
  gpr::Observation o = helper_functions::eon_matter_to_init_obs(matter);
  gprfunc.calculatePotential(o);
  for (int i = 0; i < o.E.getNumRows(); i++){
    std::cout<<std::endl;
    for (int j=0; j < o.E.getNumCols(); j++){
      std::cout<<o.E(i,j);
    }
  }
  pot.registerGPRObject(&gprfunc);
  // Call potential
  double energy_gpro{10};
  AtomMatrix forces_gpro = AtomMatrix::Constant(matter->numberOfAtoms(), 3, 0);
  auto egf_gpro = helper_functions::energy_and_forces(matter, &pot);
  energy_gpro = std::get<double>(egf_gpro);
  forces_gpro = std::get<AtomMatrix>(egf_gpro);
  EXPECT_EQ(forces_gpro, matter->getForces())
      << "Forces do not match";
  EXPECT_EQ(energy_gpro, 3)
      << "Potential energy does not match";
  delete matter;
  delete parameters;
}

} /* namespace tests */
