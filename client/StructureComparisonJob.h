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
#include "client/matter/StructComparer.hpp"
namespace eonc {
class StructureComparisonJob {
public:
  StructureComparisonJob(eonc::mat::StructComparer &sc)
      : m_sc{sc} {
    m_log = spdlog::get("combi");
  }
  ~StructureComparisonJob(void) = default;
  bool runImpl(Matter &, Matter &);

private:
  eonc::mat::StructComparer m_sc;
  std::shared_ptr<spdlog::logger> m_log;
};

} // namespace eonc
