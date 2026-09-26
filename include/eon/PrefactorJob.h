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

#include <string_view>

namespace eonc {

class PrefactorJob : public Job {
public:
  PrefactorJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~PrefactorJob() = default;
  std::vector<std::string> run() override;
  static constexpr std::string_view PREFACTOR_REACTANT{"reactant"};
  static constexpr std::string_view PREFACTOR_SADDLE{"saddle"};
  static constexpr std::string_view PREFACTOR_PRODUCT{"product"};
};

} // namespace eonc
