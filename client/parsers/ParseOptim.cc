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
#include "client/Parser.hpp"
#include "magic_enum/magic_enum.hpp"
#include <iostream>

namespace eonc::opt {
void from_toml(OptimBase::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  params.optM = get_enum_toml<OptType>(tbl["method"]).value_or(params.optM);
  params.optCM = magic_enum::enum_cast<ConvergenceMeasure>(
                     tbl["convergence_metric"].value_or("none"s),
                     magic_enum::case_insensitive)
                     .value_or(params.optCM);
  params.optConvergedForce =
      tbl["converged_force"].value_or(params.optConvergedForce);
  params.optMaxIter = tbl["max_iterations"].value_or(params.optMaxIter);
  params.optMaxMove = tbl["max_move"].value_or(params.optMaxMove);
}

void from_toml(ConjugateGradients::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  from_toml(*static_cast<OptimBase::Params *>(&params), tbl); // OptimBase
  const auto &config = tbl.at_path("CG");
  params.max_iter_before_reset =
      config["max_iter_before_reset"].value_or(params.max_iter_before_reset);
}
} // namespace eonc::opt
