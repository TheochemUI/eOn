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
#include "FiniteDifferenceJob.h"
#include "EonLogger.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Matter.h"

#include <format>
#include <fstream>
#include <stdexcept>

using namespace eonc::helpers;

std::vector<std::string> FiniteDifferenceJob::run(void) {
  auto reactant = std::make_unique<Matter>(pot, params);
  if (!eonc::io::io_ok(reactant->con2matter("pos.con"))) {
    EONC_LOG_CRITICAL("Failed to load pos.con");
    throw std::runtime_error("failed to load pos.con");
  }
  AtomMatrix posA = reactant->getPositions();

  double dRs[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, -1};

  AtomMatrix forceA = reactant->getForces();

  long epicenter = eonc::EpiCenters::minCoordinatedEpiCenter(
      reactant.get(), params.structure_comparison_options.neighbor_cutoff);
  AtomMatrix displacement;
  displacement.resize(reactant->numberOfAtoms(), 3);
  displacement.setZero();
  printf("displacing atoms:");
  for (int i = 0; i < reactant->numberOfAtoms(); i++) {
    if (reactant->distance(epicenter, i) <= 3.3) {
      printf(" %i", i);
      for (int j = 0; j < 3; j++) {
        if (!reactant->getFixed(i)) {
          displacement(i, j) = randomDouble(1.0);
        }
      }
    }
  }
  printf("\n");
  displacement.normalize();

  std::ofstream results("results.dat");
  results << std::format("{:>14s}    {:>14s}\n", "dR", "curvature");
  printf("%14s    %14s\n", "dR", "curvature");
  AtomMatrix posB;
  AtomMatrix forceB;
  double curvature = 0.0;
  for (int dRi = 0; dRs[dRi] != -1; dRi++) {
    posB = posA + displacement * dRs[dRi];
    reactant->setPositions(posB);
    forceB = reactant->getForces();
    curvature = matDot(forceB - forceA, displacement) / dRs[dRi];
    results << std::format("{:14.8f}    {:14.8f}\n", dRs[dRi], curvature);
    printf("%14.8f    %14.8f\n", dRs[dRi], curvature);
    results.flush();
  }

  std::vector<std::string> empty;
  return empty;
}
