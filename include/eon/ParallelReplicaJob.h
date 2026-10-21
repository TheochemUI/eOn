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

#include <cstddef>
#include <memory>
#include <vector>

namespace eonc {

class ParallelReplicaJob : public Job {
public:
  using Job::Job;
  ~ParallelReplicaJob() = default;
  std::vector<std::string> run() override;
  /// Matter-first entry; returns final trajectory Matter.
  std::shared_ptr<Matter> runFromMatter(std::shared_ptr<Matter> initial);

private:
  std::vector<std::string> returnFiles;
  std::shared_ptr<Matter> reactant;

  void dephase(Matter &trajectory, BondBoost *bias);
  // CG basin check. The secant step and the early stop apply only when the
  // configured minimizer is conjugate gradients. Other methods take the
  // ordinary relax path.
  void basinRelax(Matter &matter, bool quiet);
  int refineTransition(const std::vector<std::shared_ptr<Matter>> &snapshots,
                       bool fake = false);

  struct CallProfile {
    size_t initial{0};
    size_t dephaseMd{0};
    size_t dephaseQuench{0};
    size_t md{0};
    size_t stateCheck{0};
    size_t decorrelate{0};
    size_t product{0};
    long stateChecks{0};
    long hotInside{0};
    long earlyAccepts{0};
  };
  CallProfile profile_{};

  eonc::log::Scoped log;
};

} // namespace eonc
