/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
// An interface to the GPDimer library

#include "AtomicGPDimer.h"
#include "GPRHelpers.h"
#include "HelperFunctions.h"
#include "fpe_handler.h"
#include <cassert>
#include <cmath>

#include "subprojects/gpr_optim/gpr/AtomicDimer.h"
#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/structures/Structures.h"

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(std::shared_ptr<Matter> matter,
                             std::shared_ptr<Parameters> params,
                             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  matterCenter = std::make_shared<Matter>(pot, params);
  *matterCenter = *matter;
  p = helper_functions::eon_parameters_to_gpr(params.get());
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = matter->getCell()(i);
  }
}

void AtomicGPDimer::compute(std::shared_ptr<Matter> matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  atoms_config = helper_functions::eon_matter_to_atmconf(matter.get());
  // R_init.resize(1, matterCenter->getPositionsFree().size());
  // R_init.assignFromEigenMatrix(matterCenter->getPositionsFreeV());
  R_init.resize(1, matterCenter->getPositionsFree().rows() *
                       matterCenter->getPositionsFree().cols());
  int counter = 0;
  for (int i = 0; i < matterCenter->getPositionsFree().rows(); ++i) {
    for (int j = 0; j < matterCenter->getPositionsFree().cols(); ++j) {
      R_init[counter++] = matterCenter->getPositionsFree()(i, j);
    }
  }
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  problem_setup.activateFrozenAtoms(R_init, params->gprActiveRadius,
                                    atoms_config);
  orient_init.clear();
  orient_init.resize(matterCenter->getPositionsFree().rows(),
                     matterCenter->getPositionsFree().cols());
  AtomMatrix freeOrient(matterCenter->numberOfFreeAtoms(), 3);
  int i, j = 0;
  for (i = 0; i < matterCenter->numberOfAtoms(); i++) {
    if (!matterCenter->getFixed(i)) {
      freeOrient.row(j) = initialDirectionAtomMatrix.row(i);
      j++;
      if (j == matterCenter->numberOfFixedAtoms()) {
        break;
      }
    }
  }
  // orient_init.assignFromEigenMatrix(freeOrient);
  orient_init.resize(1, freeOrient.rows() * freeOrient.cols());
  counter = 0;
  for (int i = 0; i < freeOrient.rows(); ++i) {
    for (int j = 0; j < freeOrient.cols(); ++j) {
      orient_init[counter++] = freeOrient(i, j);
    }
  }
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);

  // Potential *potential = Potential::getPotential(parameters);
  auto potential = helper_functions::makePotential(params);
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();
  atomic_dimer.execute(*potential.get());
  fpeh.restore_fpe();
  // Forcefully set the right positions
  matter->setPositionsFreeV(atomic_dimer.getFinalCoordOfMidPoint());
  this->totalIterations = atomic_dimer.getIterations();
  this->totalForceCalls = atomic_dimer.getTotalForceCalls();
  pot->forceCallCounter = atomic_dimer.getTotalForceCalls();
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
