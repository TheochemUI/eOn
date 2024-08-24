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
#include "client/parsers/ParseOptim.hpp"
#include "client/BaseStructures.h"
#include "magic_enum/magic_enum.hpp"

namespace eonc::opt {
void from_toml(OptimBase::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  const auto &config = tbl.at_path("Optimizer");
  params.optM =
      magic_enum::enum_cast<OptType>(config["method"].value_or("none"s),
                                     magic_enum::case_insensitive)
          .value_or(params.optM);
  params.optCM = magic_enum::enum_cast<ConvergenceMeasure>(
                     config["convergence_metric"].value_or("none"s),
                     magic_enum::case_insensitive)
                     .value_or(params.optCM);
  params.optConvergedForce =
      config["converged_force"].value_or(params.optConvergedForce);
  params.optMaxIter = config["max_iterations"].value_or(params.optMaxIter);
  params.optMaxMove = config["max_move"].value_or(params.optMaxMove);
}
} // namespace eonc::opt
