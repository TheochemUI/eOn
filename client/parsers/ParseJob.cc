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
#include "client/BaseStructures.h"
#include "magic_enum/magic_enum.hpp"

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
} // namespace eonc::job
