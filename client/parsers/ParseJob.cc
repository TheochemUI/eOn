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

namespace eonc::job {
void from_toml(mat::StructComparer::Params &scparams,
               const toml::node_view<const toml::node> &tbl) {
  const auto &config = tbl.at_path("Structure_Comparison");
  scparams.distanceDifference =
      config["distance_difference"].value_or(scparams.distanceDifference);
  scparams.neighborCutoff =
      config["neighbor_cutoff"].value_or(scparams.neighborCutoff);
  scparams.checkRotation =
      config["check_rotation"].value_or(scparams.checkRotation);
  scparams.indistinguishableAtoms = config["indistinguishable_atoms"].value_or(
      scparams.indistinguishableAtoms);
  scparams.energyDifference =
      config["energy_difference"].value_or(scparams.energyDifference);
  scparams.removeTranslation =
      config["remove_translation"].value_or(scparams.removeTranslation);
}
} // namespace eonc::job
