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
#include "Matter.h"
namespace eonc {
class PointJob : public Job<PointJob> {
  Matter &_mat;

public:
  PointJob(Matter &mat_a)
      : _mat{mat_a} {
    log = spdlog::get("combi");
  }
  ~PointJob(void) = default;
  std::vector<std::string> run(void) override;

private:
  std::shared_ptr<spdlog::logger> log;
};

} // namespace eonc
