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
#include "eon/BaseStructures.h"
#include "eon/EonLogger.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"

#include <format>
#include <fstream>
#include <stdexcept>

std::vector<std::string> FiniteDifferenceJob::run(void) {
  auto reactant = std::make_unique<Matter>(pot, params);
  const std::string posFile = eonc::helpers::getRelevantFile("pos.con");
  if (!eonc::io::io_ok(reactant->con2matter(posFile))) {
    EONC_LOG_CRITICAL("Failed to load {}", posFile);
    throw std::runtime_error("failed to load " + posFile);
  }
  AtomMatrix posA = reactant->getPositions();

  double dRs[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, -1};

  AtomMatrix forceA = reactant->getForces();

  const double cutoff = params.structure_comparison_options.neighbor_cutoff;
  long epicenter =
      eonc::EpiCenters::minCoordinatedEpiCenter(reactant.get(), cutoff);
  AtomMatrix displacement;
  displacement.resize(reactant->numberOfAtoms(), 3);
  displacement.setZero();
  printf("displacing atoms:");
  for (int i = 0; i < reactant->numberOfAtoms(); i++) {
    if (reactant->distance(epicenter, i) <= cutoff) {
      printf(" %i", i);
      for (int j = 0; j < 3; j++) {
        if (!reactant->getFixed(i)) {
          displacement(i, j) = eonc::rng::randomDouble(1.0);
        }
      }
    }
  }
  printf("\n");
  const double dispNorm = displacement.norm();
  if (!(dispNorm > 0.0)) {
    throw std::runtime_error(
        "FiniteDifferenceJob: no free atoms in the epicenter neighborhood");
  }
  displacement /= dispNorm;

  auto env = JobResultEnvelope::fromMinimization(
      RunStatus::GOOD, params.potential_options.potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  env.job_type = "finite_difference";

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
    env.extras.emplace_back(std::format("dR_{}", dRi), dRs[dRi]);
    env.extras.emplace_back(std::format("curvature_{}", dRi), curvature);
    printf("%14.8f    %14.8f\n", dRs[dRi], curvature);
    table.flush();
  }

  env.force_calls = PotRegistry::get().total_force_calls();
  env.writeResultsDat("results.dat");
  std::vector<std::string> returnFiles;
  returnFiles.push_back("results.dat");
  returnFiles.push_back("curvature.dat");
  return returnFiles;
}
