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
#include "client/OptParams.hpp"
#include <variant>

namespace eonc::opt {
void from_toml(ConjugateGradients::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  extract_common_params(params, tbl); // Extract common params
  const auto &config = tbl.at_path("CG");
  params.max_iter_before_reset =
      config["max_iter_before_reset"].value_or(params.max_iter_before_reset);
}

void from_toml(std::monostate &params,
               const toml::node_view<const toml::node> &tbl) {
  throw std::logic_error("Incorrect OptParams variant for from_toml");
}

void from_toml(OptParams &params,
               const toml::node_view<const toml::node> &tbl) {

  if (std::holds_alternative<std::monostate>(params)) {
    // Make parameters
    auto otype = get_enum_toml<OptType>(tbl["method"]).value();
    switch (otype) {
    case OptType::CG: {
      params = ConjugateGradients::Params();
      break;
    }
    case OptType::Unknown: {
      throw std::runtime_error("[Optimizer] Cannot create Unknown");
    }
    default: {
      throw std::runtime_error("Unsupported optimization type");
    }
    }
  }
  std::visit([&](auto &optp) -> void { return from_toml(optp, tbl); }, params);
}

BaseOptParams OptBaseVisitor::operator()(std::monostate) {
  throw std::logic_error("Incorrect OptParams variant for OptBaseVisitor");
}

BaseOptParams
OptBaseVisitor::operator()(const ConjugateGradients::Params &params) {
  // intentional slicing to reduce typing
  return BaseOptParams(static_cast<BaseOptParams>(params));
}

} // namespace eonc::opt
