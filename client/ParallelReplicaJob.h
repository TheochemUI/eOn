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

#include "Job.h"
#include "Matter.h"

#include <memory>
#include <vector>

namespace eonc {

class ParallelReplicaJob : public Job {
public:
  ParallelReplicaJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~ParallelReplicaJob() = default;
  std::vector<std::string> run() override;

private:
  std::vector<std::string> returnFiles;
  std::shared_ptr<Matter> reactant;

  void dephase(Matter &trajectory);
  int refineTransition(const std::vector<std::shared_ptr<Matter>> &snapshots,
                       bool fake = false);
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::ParallelReplicaJob;
