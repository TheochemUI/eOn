// An interface to the GPDimer library

#include "AtomicGPDimer.h"
#include "HelperFunctions.h"
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
  *matterCenter = *matter;
  p = eon_parameters_to_gpr(params);
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = matter->getCell()[i];
  }
}

AtomicGPDimer::~AtomicGPDimer() { delete matterCenter; }

void AtomicGPDimer::compute(Matter *matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  atoms_config = eon_matter_to_atmconf(matter);
  *matterCenter = *matter;
  R_init.resize(matterCenter->getPositionsFree().rows(),
                matterCenter->getPositionsFree().cols());
  R_init.assignFromEigenMatrix(matterCenter->getPositionsFree());
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  problem_setup.activateFrozenAtoms(R_init, parameters->gprActiveRadius,
                                    atoms_config);
  // FIXME: does this work?
  orient_init.clear();
  orient_init.resize(matterCenter->getPositionsFree().rows(),
                     matterCenter->getPositionsFree().cols());
    AtomMatrix freeOrient(matterCenter->numberOfFreeAtoms(),3);
    int i,j = 0;
    for(i=0; i<matterCenter->numberOfAtoms(); i++)
    {
        if(!matterCenter->getFixed(i))
        {
            initialDirectionAtomMatrix.row(j) = freeOrient.row(i);
            j++;
        }
    }
  orient_init.assignFromEigenMatrix(freeOrient);
  // orient_init.resize(1, 3 * matter->numberOfFreeAtoms());
  // for (auto i = 0; i < matter->numberOfAtoms(); i++) {
  //   if (matter->getFixed(i)==0) {
  //     orient_init.set(0,)=initialDirectionAtomMatrix[i];
  //   }

  // }
  // vector<double> vec(initialDirection.data(),
  //                    initialDirection.data() + initialDirection.size());
  // orient_init.assignToSlice(0, vec);
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);

  Potential *potential= Potential::getPotential(parameters);
  atomic_dimer.execute(*potential);
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
