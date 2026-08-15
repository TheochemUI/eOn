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
#include "eon/HessianJob.h"
#include "eon/BaseStructures.h"
#include "eon/EonLogger.h"
#include "eon/Hessian.h"
#include "eon/Matter.h"
#include "eon/MobileAtoms.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "magic_enum/magic_enum.hpp"

#include <filesystem>
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

  // [Hessian] phva_atoms = PHVA mobile/active set (displaced in FD). free/fixed
  // is the optimizer mask; resolveMobileAtoms intersects the list with free.
  const VectorXi mobile =
      eonc::resolveMobileAtoms(matter.get(), params.hessian_options.phva_atoms);
  const bool no_mobile = mobile.size() == 0;
  if (!no_mobile) {
    hessian.getFreqs(matter.get(), mobile);
  }

  std::string results_file("results.dat");
  returnFiles.push_back(results_file);

  std::ofstream out(results_file, std::ios::binary);
  if (out) {
    const auto status =
        no_mobile ? RunStatus::FAIL_MAX_ITERATIONS : RunStatus::GOOD;
    out << std::format("{} termination_reason\n", static_cast<int>(status));
    out << std::format("{} termination_reason_text\n",
                       magic_enum::enum_name<RunStatus>(status));
    out << "hessian job_type\n";
    out << std::format("{} force_calls\n",
                       PotRegistry::get().total_force_calls());
    out << std::format("{} total_force_calls\n",
                       PotRegistry::get().total_force_calls());
  }
  if (std::filesystem::exists("hessian.dat")) {
    returnFiles.push_back("hessian.dat");
  }

  return returnFiles;
}
