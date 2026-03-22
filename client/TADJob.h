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

namespace eonc {

class TADJob : public ReplicaDynamicsJob {
public:
  TADJob(std::unique_ptr<Parameters> parameters)
      : ReplicaDynamicsJob(std::move(parameters)) {}
  ~TADJob() = default;

private:
  int dynamics() override;
  void initExtra() override;
  void reportResults() override;

  std::shared_ptr<Matter> crossing;
  double barrier{0.0};
  std::vector<double> timeBuffer;
};

} // namespace eonc

using eonc::TADJob;
