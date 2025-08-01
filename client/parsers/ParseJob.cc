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
#include "client/parsers/ParseJob.hpp"
#include "client/parsers/ParseOptim.hpp"

namespace eonc::job {
void from_toml(mat::StructComparer::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  const auto &config = tbl.at_path("Structure_Comparison");
  params.distanceDifference =
      config["distance_difference"].value_or(params.distanceDifference);
  params.neighborCutoff =
      config["neighbor_cutoff"].value_or(params.neighborCutoff);
  params.checkRotation =
      config["check_rotation"].value_or(params.checkRotation);
  params.indistinguishableAtoms =
      config["indistinguishable_atoms"].value_or(params.indistinguishableAtoms);
  params.energyDifference =
      config["energy_difference"].value_or(params.energyDifference);
  params.removeTranslation =
      config["remove_translation"].value_or(params.removeTranslation);
}

void from_toml(MinimizationJob::Params &params,
               const toml::node_view<const toml::node> &tbl) {
  std::cout << tbl;
  opt::from_toml(params.oparams, tbl);
}
} // namespace eonc::job
