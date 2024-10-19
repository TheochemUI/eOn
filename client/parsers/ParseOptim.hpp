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

#include "client/ConjugateGradients.h"
#include "client/OptParams.hpp"
#include "client/Parser.hpp"

namespace eonc::opt {
// General template for extracting parameters from TOML
template <typename ParamsType>
void extract_common_params(ParamsType &params,
                           const toml::node_view<const toml::node> &tbl) {
  params.optM = get_enum_toml<OptType>(tbl["method"]).value_or(params.optM);
  params.optCM = get_enum_toml<ConvergenceMeasure>(tbl["convergence_metric"])
                     .value_or(params.optCM);
  params.optConvergedForce =
      tbl["converged_force"].value_or(params.optConvergedForce);
  params.optMaxIter = tbl["max_iterations"].value_or(params.optMaxIter);
  params.optMaxMove = tbl["max_move"].value_or(params.optMaxMove);
}

void from_toml(OptParams &, const toml::node_view<const toml::node> &);
void from_toml(ConjugateGradients::Params &,
               const toml::node_view<const toml::node> &);

struct OptBaseVisitor {
  BaseOptParams operator()(std::monostate);
  BaseOptParams operator()(const ConjugateGradients::Params &);
};

} // namespace eonc::opt
