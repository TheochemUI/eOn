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
namespace eonc {
class TestJob : public Job {
public:
  TestJob(std::unique_ptr<Parameters> params)
      : Job(std::move(params)),
        tolerance{0.01} {}
  ~TestJob(void) = default;
  std::vector<std::string> run(void);

private:
  double tolerance;
  void checkFullSearch(void);
  void checkMode(void);
  void checkPotentials(void);
  double getEnergyDiff(std::string potTag, double refEnergy);
  double getForceDiff(std::string potTag, double refForce);
};

} // namespace eonc
