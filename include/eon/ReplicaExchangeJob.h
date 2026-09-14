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
#include "EonLogger.h"

#include "Dynamics.h"
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

namespace eonc {

class ReplicaExchangeJob : public Job {
public:
  using Job::Job;
  ~ReplicaExchangeJob() = default;
  std::vector<std::string> run(void);
  /// Matter-first; returns replica-0 Matter after sampling.
  std::shared_ptr<Matter> runFromMatter(std::shared_ptr<Matter> initial);

private:
  void saveData();

  size_t forceCalls{0};
  //        std::shared_ptr<Matter>*replica;
  std::shared_ptr<Matter> pos;
  //        Dynamics **replicaDynamics;
  //        double *replicaTemperature;
  std::vector<std::string> returnFiles;
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::ReplicaExchangeJob;
