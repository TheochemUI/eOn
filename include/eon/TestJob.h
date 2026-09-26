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

#include <memory>
#include <utility>

namespace eonc {

class TestJob : public Job {
public:
  TestJob(std::unique_ptr<Parameters> params)
      : Job(std::move(params)), tolerance{0.01} {}
  ~TestJob() = default;
  std::vector<std::string> run();

private:
  double tolerance{0.01};
  void checkFullSearch();
  void checkPotentials();
  double getEnergyDiff(const std::string &potTag, double refEnergy);
  double getForceDiff(const std::string &potTag, double refForce);
};

} // namespace eonc
