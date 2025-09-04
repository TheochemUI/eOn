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

class FiniteDifferenceJob : public Job {
public:
  FiniteDifferenceJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~FiniteDifferenceJob(void) = default;
  std::vector<std::string> run(void);
};
