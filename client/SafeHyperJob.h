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
#include "ReplicaDynamicsJob.h"

#include <vector>

namespace eonc {

class SafeHyperJob : public ReplicaDynamicsJob {
public:
  SafeHyperJob(std::unique_ptr<Parameters> parameters)
      : ReplicaDynamicsJob(std::move(parameters)) {}
  ~SafeHyperJob(void) = default;

private:
  int dynamics() override;
  void reportResults() override;

  double transitionPot{0.0};
  std::vector<double> timeBuffer;
  std::vector<double> biasBuffer;
};

} // namespace eonc

using eonc::SafeHyperJob;
