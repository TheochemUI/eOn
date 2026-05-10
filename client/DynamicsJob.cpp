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
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "Potential.h"

std::vector<std::string> DynamicsJob::run() {
  auto R = std::make_shared<Matter>(pot, params);
  auto F = std::make_shared<Matter>(pot, params);
  R->con2matter("pos.con");
  *F = *R;

  auto d = std::make_unique<Dynamics>(R.get(), params);
  d->run();

  *F = *R;
  std::string productFilename("final.con");
  F->matter2con(productFilename);

  std::vector<std::string> returnFiles;
  returnFiles.push_back(productFilename);
  return returnFiles;
}
