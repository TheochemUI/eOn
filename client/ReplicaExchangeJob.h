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

#include "Dynamics.h"
#include "Job.h"
#include "Parameters.h"

class ReplicaExchangeJob : public Job {
public:
  ReplicaExchangeJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        forceCalls{0} {
    log = spdlog::get("combi");
  }
  ~ReplicaExchangeJob(void) = default;
  std::vector<std::string> run(void);

private:
  void saveData();

  size_t forceCalls;
  //        std::shared_ptr<Matter>*replica;
  std::shared_ptr<Matter> pos;
  //        Dynamics **replicaDynamics;
  //        double *replicaTemperature;
  std::vector<std::string> returnFiles;
  shared_ptr<spdlog::logger> log;
};
