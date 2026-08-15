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
#include "eon/FiniteDifferenceJob.h"
#include "eon/EonLogger.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"

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
  results << "0 termination_reason\n";
  results << "GOOD termination_reason_text\n";
  results << "finite_difference job_type\n";
  results << std::format("{} total_force_calls\n",
                         PotRegistry::get().total_force_calls());

  std::ofstream table("curvature.dat");
  table << std::format("{:>14s}    {:>14s}\n", "dR", "curvature");
  printf("%14s    %14s\n", "dR", "curvature");
  AtomMatrix posB;
  AtomMatrix forceB;
  double curvature = 0.0;
  for (int dRi = 0; dRs[dRi] != -1; dRi++) {
    posB = posA + displacement * dRs[dRi];
    reactant->setPositions(posB);
    forceB = reactant->getForces();
    curvature = matDot(forceB - forceA, displacement) / dRs[dRi];
    table << std::format("{:14.8f}    {:14.8f}\n", dRs[dRi], curvature);
    results << std::format("{:.12e} dR_{}\n", dRs[dRi], dRi);
    results << std::format("{:.12e} curvature_{}\n", curvature, dRi);
    printf("%14.8f    %14.8f\n", dRs[dRi], curvature);
    table.flush();
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back("results.dat");
  returnFiles.push_back("curvature.dat");
  return returnFiles;
}
