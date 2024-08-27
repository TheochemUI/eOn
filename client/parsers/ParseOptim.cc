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

namespace eonc::opt {
void from_toml(ConjugateGradients::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  extract_common_params(params, tbl); // Extract common params
  const auto &config = tbl.at_path("CG");
  params.max_iter_before_reset =
      config["max_iter_before_reset"].value_or(params.max_iter_before_reset);
}
} // namespace eonc::opt
