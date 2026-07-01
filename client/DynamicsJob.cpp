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
#include <stdexcept>
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "EonLogger.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "Potential.h"

std::vector<std::string> DynamicsJob::run(void) {
  auto R = std::make_shared<Matter>(pot, params);
  auto F = std::make_shared<Matter>(pot, params);
  if (!eonc::io::io_ok(R->con2matter("pos.con"))) {
    EONC_LOG_CRITICAL("Failed to load pos.con");
    throw std::runtime_error("failed to load pos.con");
  }
  *F = *R;

  auto d = std::make_unique<Dynamics>(R.get(), params);
  d->run();

  *F = *R;
  std::string productFilename("final.con");
  if (!eonc::io::io_ok(F->matter2con(productFilename))) {
    EONC_LOG_ERROR("Failed to write {}", productFilename);
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(productFilename);
  return returnFiles;
}
