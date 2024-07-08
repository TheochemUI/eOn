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
#include <stdio.h>
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "Potential.h"

using namespace helper_functions;

std::vector<std::string> DynamicsJob::run(void) {
  auto R = std::make_shared<Matter>(pot, params);
  auto F = std::make_shared<Matter>(pot, params);
  R->con2matter("pos.con");
  *F = *R;

  Dynamics *d = new Dynamics(R.get(), params.get());
  d->run();

  *F = *R;
  FILE *fileProduct;
  std::string productFilename("final.con");
  returnFiles.push_back(productFilename);

  fileProduct = fopen(productFilename.c_str(), "wb");
  F->matter2con(fileProduct);
  fclose(fileProduct);

  delete d;

  std::vector<std::string> returnFiles;
  return returnFiles;
}
