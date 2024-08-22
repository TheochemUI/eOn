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

#include "Optimizer.h"
#include "client/matter/Matter.h"

namespace eonc {
class RelaxJob {
public:
  RelaxJob(const OptimBase::Params &_p)
      : m_p{_p} {
    m_log = spdlog::get("combi");
  }
  ~RelaxJob(void) = default;
  bool runImpl(Matter &);

private:
  OptimBase::Params m_p;
  std::shared_ptr<spdlog::logger> m_log;
};

} // namespace eonc
