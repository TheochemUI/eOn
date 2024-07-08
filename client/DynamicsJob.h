/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
#pragma once
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class DynamicsJob : public Job {

public:
  DynamicsJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~DynamicsJob(void) = default;
  std::vector<std::string> run(void);
  std::vector<std::string> returnFiles;
};
