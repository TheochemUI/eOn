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

class MinimizationJob : public Job {
public:
  MinimizationJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        fcalls{0} {
    log = spdlog::get("combi");
  }
  ~MinimizationJob(void) = default;
  std::vector<std::string> run(void);

private:
  size_t fcalls;
  RunStatus status;
  shared_ptr<spdlog::logger> log;
};
