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
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/MobileAtoms.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"

#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

namespace eonc {

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
  bool freqs_ok = false;
  if (!no_mobile) {
    const VectorXd freqs = hessian.getFreqs(matter.get(), mobile);
    freqs_ok = freqs.size() > 0;
  }

  std::string results_file("results.dat");
  returnFiles.push_back(results_file);

  const auto status =
      freqs_ok ? RunStatus::GOOD : RunStatus::FAIL_POTENTIAL_FAILED;
  auto env = JobResultEnvelope::fromMinimization(
      status, params.potential_options.potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  env.job_type = "hessian";
  env.extras.emplace_back("force_calls", static_cast<double>(env.force_calls));
  env.writeResultsDat(results_file);
  if (std::filesystem::exists("hessian.dat")) {
    returnFiles.push_back("hessian.dat");
  }

  return returnFiles;
}

} // namespace eonc
