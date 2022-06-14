// an interface to the gpdimer library

#include "AtomicGPDimer.h"
#include "HelperFunctions.h"
#include "GPRHelpers.h"
#include "Log.h"
#include <cassert>
#include <cmath>

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(Matter *matter, Parameters *params) {
  parameters = params;
  matterCenter = new Matter(parameters);
  *matterCenter = *matter;
  p = helper_functions::eon_parameters_to_gprd(params);
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = matter->getCell()(i);
  }
  this->freeOrient.resize(matterCenter->numberOfFreeAtoms(), 3);
}

AtomicGPDimer::~AtomicGPDimer() { delete matterCenter; }

void AtomicGPDimer::compute(Matter *matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  *matterCenter = *matter;
  for (size_t idx{0}, jdx{0}; idx < static_cast<size_t>(matterCenter->numberOfAtoms()); idx++)
    {
        if(!matterCenter->getFixed(idx))
        {
           freeOrient.row(jdx) = initialDirectionAtomMatrix.row(idx);
            jdx++;
            if (jdx==static_cast<size_t>(matterCenter->numberOfFixedAtoms())){
              break;
            }
        }
    }
  auto [atoms_config, R_init] = helper_functions::eon_matter_to_frozen_conf_info(matterCenter,
                                                             parameters->gprActiveRadius);
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  orient_init.clear();
  orient_init.resize(1, matterCenter->getPositionsFree().size());
  orient_init.resize(1, this->freeOrient.size());
  for (size_t idx{0}; idx < freeOrient.size(); ++idx) {
    orient_init(0, idx) = freeOrient.reshaped<Eigen::RowMajor>()[idx];
  }
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);

  Potential *potential= Potential::getPotential(parameters);
  atomic_dimer.execute(*potential);
  // Forcefully set the right positions
  matter->setPositionsFreeV(atomic_dimer.getFinalCoordOfMidPoint());
  this->totalIterations=atomic_dimer.getIterations();
  this->totalForceCalls=atomic_dimer.getTotalForceCalls();
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
