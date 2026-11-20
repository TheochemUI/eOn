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

/// Energy-continuation line-integral NEB paths and harmonic quantum
/// transition-state rates ([Line Integral Tunneling]).
class LineIntegralTunnelingJob final : public Job {
public:
  LineIntegralTunnelingJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt) {}
  ~LineIntegralTunnelingJob() = default;
  std::vector<std::string> run() override;
};

} // namespace eonc
