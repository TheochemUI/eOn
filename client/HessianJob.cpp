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
#include "HessianJob.h"
#include "EonLogger.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"

#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

std::vector<std::string> HessianJob::run(void) {
  std::string matter_in("pos.con");

  std::vector<std::string> returnFiles;

  auto matter = std::make_unique<Matter>(pot, params);

  if (!eonc::io::io_ok(matter->con2matter(matter_in))) {
    EONC_LOG_CRITICAL("Failed to load {}", matter_in);
    throw std::runtime_error("failed to load " + matter_in);
  }

  Hessian hessian(params, matter.get());
  long nAtoms = matter->numberOfAtoms();

  VectorXi moved(nAtoms);
  moved.setConstant(-1);

  int nMoved = 0;
  for (int i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      moved[nMoved] = i;
      nMoved++;
    }
  }
  moved = moved.head(nMoved);
  hessian.getFreqs(matter.get(), moved);

  std::string results_file("results.dat");
  returnFiles.push_back(results_file);

  std::ofstream out(results_file, std::ios::binary);
  if (out) {
    out << std::format("{} force_calls\n",
                       PotRegistry::get().total_force_calls());
  }

  return returnFiles;
}
