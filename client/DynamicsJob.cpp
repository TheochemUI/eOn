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
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

#include "eon/BaseStructures.h"
#include "eon/Dynamics.h"
#include "eon/DynamicsJob.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "magic_enum/magic_enum.hpp"

std::vector<std::string> DynamicsJob::run(void) {
  auto R = std::make_shared<Matter>(pot, params);
  const std::string posFile = eonc::helpers::getRelevantFile("pos.con");
  if (!eonc::io::io_ok(R->con2matter(posFile))) {
    EONC_LOG_CRITICAL("Failed to load {}", posFile);
    throw std::runtime_error("failed to load " + posFile);
  }

  auto d = std::make_unique<Dynamics>(R.get(), params);
  d->run();

  std::string productFilename("final.con");
  if (!eonc::io::io_ok(R->matter2con(productFilename))) {
    EONC_LOG_ERROR("Failed to write {}", productFilename);
  }

  const std::string resultsFilename("results.dat");
  std::ofstream out(resultsFilename, std::ios::binary);
  if (!out) {
    throw std::runtime_error("failed to open " + resultsFilename);
  }
  out << std::format("{} termination_reason\n",
                     static_cast<int>(RunStatus::GOOD));
  out << std::format("{} termination_reason_text\n",
                     magic_enum::enum_name<RunStatus>(RunStatus::GOOD));
  out << "dynamics job_type\n";
  out << std::format("{} total_force_calls\n",
                     PotRegistry::get().total_force_calls());
  out << std::format("{:.12f} potential_energy\n", R->getPotentialEnergy());
  out.close();
  if (!out) {
    throw std::runtime_error("failed to write " + resultsFilename);
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(productFilename);
  returnFiles.push_back(resultsFilename);
  return returnFiles;
}
