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
#include "Parameters.h"

class PointJob : public Job {
public:
  PointJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~PointJob(void) = default;
  std::vector<std::string> run(void);

private:
  std::shared_ptr<spdlog::logger> log;
};
