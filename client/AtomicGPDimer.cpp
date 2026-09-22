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

#include "eon/AtomicGPDimer.h"
#include "eon/GPRHelpers.h"
#include "eon/HelperFunctions.h"
#include "eon/fpe_handler.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "subprojects/gpr_optim/gpr/AtomicDimer.h"
#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gpr_optim/structures/Structures.h"

namespace eonc {

namespace {

// AtomMatrix is row-major N×3; gpr::Coord is row-major 1×(3N) with the same
// flat packing [x0,y0,z0,x1,...].
void copyAtomMatrixToCoord(const AtomMatrix &src, gpr::Coord &dst) {
  dst.resize(1, static_cast<Eigen::Index>(src.size()));
  if (src.size() > 0) {
    std::memcpy(dst.data(), src.data(),
                static_cast<size_t>(src.size()) * sizeof(double));
  }
}

} // namespace

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(std::shared_ptr<Matter> matter,
                             const Parameters &params,
                             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  if (!matter) {
    throw std::invalid_argument("AtomicGPDimer: null Matter");
  }
  matterCenter = std::make_shared<Matter>(pot, params);
  *matterCenter = *matter;
  p = eonc::helpers::eon_parameters_to_gpr(params);
  const Matrix3d cell = matter->getCell();
  for (int i = 0; i < 9; i++) {
    p.cell_dimensions.value[i] = cell.data()[i];
  }
}

void AtomicGPDimer::compute(std::shared_ptr<Matter> matter,
                            AtomMatrix initialDirectionAtomMatrix) {
  atoms_config = eonc::helpers::eon_matter_to_atmconf(matter.get());
  copyAtomMatrixToCoord(matterCenter->getPositionsFree(), R_init);
  init_middle_point.clear();
  init_middle_point.R = R_init;
  init_observations.clear();
  problem_setup.activateFrozenAtoms(
      R_init, params.gpr_dimer_options().active_radius, atoms_config);
  AtomMatrix freeOrient(matterCenter->numberOfFreeAtoms(), 3);
  int j = 0;
  for (int i = 0; i < matterCenter->numberOfAtoms(); i++) {
    if (!matterCenter->getFixed(i)) {
      freeOrient.row(j) = initialDirectionAtomMatrix.row(i);
      j++;
      if (j == matterCenter->numberOfFreeAtoms()) {
        break;
      }
    }
  }
  copyAtomMatrixToCoord(freeOrient, orient_init);
  atomic_dimer.initialize(p, init_observations, init_middle_point, orient_init,
                          atoms_config);

  auto potential = eonc::helpers::makePotential(params);
  pot::PotentialWrapper wrapper(
      [&potential](long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
        potential->force(N, R, atomicNrs, F, U, variance, box);
      });
  // Restore traps if execute throws. The saddle-search catch must not
  // leave later force calls running with traps still masked.
  {
    eonc::FPEGuard fpe;
    atomic_dimer.execute(wrapper);
  }
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
  const gpr::Coord &orient = atomic_dimer.getFinalOrientation();
  const long nFree = matterCenter->numberOfFreeAtoms();
  const long nAtoms = matterCenter->numberOfAtoms();
  if (nFree <= 0 || orient.size() != 3 * nFree) {
    return AtomMatrix::Zero(nAtoms, 3);
  }
  AtomMatrix freeMode = Eigen::Map<const AtomMatrix>(orient.data(), nFree, 3);
  if (nFree == nAtoms) {
    return freeMode;
  }
  AtomMatrix full = AtomMatrix::Zero(nAtoms, 3);
  long k = 0;
  for (long i = 0; i < nAtoms; ++i) {
    if (!matterCenter->getFixed(i)) {
      full.row(i) = freeMode.row(k++);
    }
  }
  return full;
}

} // namespace eonc
