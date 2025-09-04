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
#pragma once

#include "Job.h"
#include "Parameters.h"

class PrefactorJob : public Job {
public:
  PrefactorJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~PrefactorJob(void) = default;
  std::vector<std::string> run(void);
  // Ugly but OK for now I guess
  static const char PREFACTOR_REACTANT[];
  static const char PREFACTOR_SADDLE[];
  static const char PREFACTOR_PRODUCT[];
};
