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

/// Ring-polymer instanton ([Instanton]): the tunnelling splitting between
/// two minima beyond one-dimensional WKB (mode splitting), or the thermal
/// rate through a saddle (mode rate): the ring below the crossover
/// temperature, and the parabolic barrier factor times harmonic
/// transition-state theory above it, written to instanton.con and
/// results.dat.
class InstantonJob : public Job {
public:
  InstantonJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt) {}
  ~InstantonJob(void) = default;
  std::vector<std::string> run(void) override;
};

} // namespace eonc
