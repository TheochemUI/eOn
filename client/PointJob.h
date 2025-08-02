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

#include "client/matter/Matter.h"
namespace eonc {
class PointJob {

public:
  PointJob() { m_log = spdlog::get("combi"); }
  ~PointJob(void) = default;
  bool runImpl(Matter &);

private:
  std::shared_ptr<spdlog::logger> m_log;
};

} // namespace eonc
