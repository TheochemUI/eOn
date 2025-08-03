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
#include "client/StructureComparisonJob.h"

namespace eonc {
bool StructureComparisonJob::runImpl(Matter &m1, Matter &m2) {
  bool res = m_sc.compare(m1, m2, /*indistinguishable=*/true);
  if (res) {
    SPDLOG_LOGGER_INFO(m_log, "Structures match");
  } else {
    SPDLOG_LOGGER_INFO(m_log, "Structures do not match");
  }
  return res;
}

} // namespace eonc
