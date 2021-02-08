// An interface to the GPDimer library

#include "AtomicGPDimer.h"
#include "HelperFunctions.h"
#include "potentials/Morse/Morse.h"
#include "Log.h"
#include <cassert>
#include <cmath>

#include "gprdimer/gpr/AtomicDimer.h"
#include "gprdimer/gpr/Enums.h"
#include "gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "gprdimer/gpr/ml/GaussianProcessRegression.h"

using namespace helper_functions;

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(Matter *matter, Parameters *params) {
  parameters = params;
  matterCenter = new Matter(parameters);
  matterDimer = new Matter(parameters);
  *matterCenter = *matter;
  *matterDimer = *matter;
  InputParameters p = eon_parameters_to_gpr(params);
  atmd::AtomicDimer atomic_dimer;
  aux::ProblemSetUp problem_setup;
}

AtomicGPDimer::~AtomicGPDimer() {
  delete matterCenter;
  delete matterDimer;
}

void AtomicGPDimer::compute(Matter *matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  AtomsConfiguration atoms_config;
  Observation init_observations, init_middle_point;
  gpr::Coord orient_init, R_init;
  aux::ProblemSetUp problem_setup;
  atoms_config = eon_matter_to_atmconf(matter);
  VectorXd initialDirection = VectorXd::Map(initialDirectionAtomMatrix.data(),
                                            3 * matter->numberOfAtoms());
  VectorXd tau = initialDirection.array() * matter->getFreeV().array();
  tau = initialDirection;
  tau.normalize();
  *matterCenter = *matter;
  *matterDimer = *matter;
  VectorXd x0_r = matterCenter->getPositionsV();
  matterDimer->setPositionsV(x0_r + parameters->finiteDifference * tau);
  R_init.resize(matterDimer->getPositionsV().rows(),
                matterDimer->getPositionsV().cols());
  R_init.assignFromEigenMatrix(matterDimer->getPositionsV());
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  problem_setup.activateFrozenAtoms(R_init, parameters->gprActiveRadius,
                                    atoms_config);
  // FIXME: does this work?
  orient_init.clear();
  orient_init.resize(1, 3 * matter->numberOfFreeAtoms());
  vector<double> vec(initialDirection.data(),
                     initialDirection.data() + initialDirection.size());
  orient_init.assignToSlice(0, vec);
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);
  Potential *potential = Potential::getPotential(parameters);
  Morse *mot;
  // FIXME: Broken call
  atomic_dimer.execute(*mot);
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
