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
#include "eon/BasinHoppingSaddleSearch.h"
#include "eon/Dimer.h"
#include "eon/ImprovedDimer.h"
#include "eon/Lanczos.h"
#include "eon/LowestEigenmode.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/NudgedElasticBand.h"
#include <cmath>
#include <cstddef>
#include <cstdio>

namespace eonc {

int BasinHoppingSaddleSearch::highestEnergyInteriorImage(
    const std::vector<std::shared_ptr<Matter>> &path, long numImages) {
  double emax = -1e100;
  int highest = 0;
  for (long i = 1; i <= numImages; i++) {
    double etest = path[static_cast<size_t>(i)]->getPotentialEnergy();
    QUILL_LOG_DEBUG(eonc::log::get(), "i: {} Etest: {:.1f}", i, etest);
    if (etest > emax) {
      emax = etest;
      highest = static_cast<int>(i);
    }
  }
  return highest;
}

int BasinHoppingSaddleSearch::run() {
  // minimize "saddle"
  saddle->relax(false, true, false, "displacementmin");
  product = std::make_shared<Matter>(pot, params);
  *product = *saddle;
  // accept or reject based on boltzman
  // exp(-de/(kB*params.main_options().temperature))
  double eproduct, ereactant, de;
  eproduct = product->getPotentialEnergy();
  ereactant = reactant->getPotentialEnergy();
  de = eproduct - ereactant;
  double kB = params.constants().kB;
  double Temperature = params.main_options().temperature;
  double arg = -de / (kB * Temperature);
  double p = std::exp(arg);
  double r = eonc::rng::random();
  if (ereactant < eproduct) {
    if (r > p) { // reject
      status = 1;
      return status;
    }
  }
  // NEB reactant to minimized "saddle"
  NudgedElasticBand neb(reactant, product, params, pot);
  if (!eonc::io::io_ok(
          neb.path[0]->matter2con("neb_initial_band.con", false))) {
    QUILL_LOG_WARNING(log, "Failed to write neb_initial_band.con");
  }
  // Reactant is frame 0. Interiors are 1..numImages, including the last.
  for (int j = 1; j <= neb.numImages; j++) {
    if (!eonc::io::io_ok(
            neb.path[j]->matter2con("neb_initial_band.con", true))) {
      QUILL_LOG_WARNING(log, "Failed to append neb_initial_band frame");
    }
  }
  neb.compute();
  // pick the maximum energy image along the band
  int HighestImage = highestEnergyInteriorImage(neb.path, neb.numImages);
  if (HighestImage < 1) {
    QUILL_LOG_WARNING(log, "No interior NEB image for basin hopping");
    status = MinModeSaddleSearch::STATUS_BAD_NO_BARRIER;
    return status;
  }
  // do dimer
  // Calculate initial direction
  AtomMatrix r_1 = neb.path[HighestImage - 1]->getPositions();
  AtomMatrix r_3 = neb.path[HighestImage + 1]->getPositions();
  AtomMatrix direction =
      initialDimerDirection(*neb.path[HighestImage], r_1, r_3);
  MinModeSaddleSearch dim(neb.path[HighestImage], direction.normalized(),
                          ereactant, params, pot);
  // ProcessSearchJob treats STATUS_GOOD as a saddle. The climb's own
  // status is that decision; discarding it records a failed climb as found.
  status = dim.run();
  *saddle = *neb.path[HighestImage];
  eigenvalue = dim.getEigenvalue();
  eigenvector = dim.getEigenvector();
  return status;
}

double BasinHoppingSaddleSearch::getEigenvalue() { return eigenvalue; }

AtomMatrix BasinHoppingSaddleSearch::getEigenvector() { return eigenvector; }

} // namespace eonc
