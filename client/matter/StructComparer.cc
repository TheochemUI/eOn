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
#include "client/matter/StructComparer.hpp"
#include "client/Parser.hpp"
#include "client/matter/MatterHelpers.hpp"

namespace eonc::mat {

void StructComparer::fromTOML(const toml::table &tbl) {
  config_section(tbl, "Structure_Comparison");
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

void StructComparer::setupCompareFunc() {
  compareFunc = [this](const Matter &m1Ref, const Matter &m2Ref,
                       double distanceDifference) -> bool {
    if (scparams.checkRotation && scparams.indistinguishableAtoms) {
      return sortedR(m1Ref, m2Ref, distanceDifference);
    } else if (scparams.indistinguishableAtoms) {
      return identical(m1Ref, m2Ref, distanceDifference);
    } else if (scparams.checkRotation) {
      return rotationMatch(m1Ref, m2Ref, distanceDifference);
    } else {
      return distanceDifference > m1Ref.perAtomNorm(m2Ref);
    }
  };
}

bool StructComparer::compare(const Matter &m1, const Matter &m2,
                             const bool indistinguishable) {
  if (m1.numberOfAtoms() != m2.numberOfAtoms()) {
    return false;
  }

  std::optional<std::pair<Matter, Matter>> maybeCopies;

  if (m1.numberOfFixedAtoms() == 0 && scparams.removeTranslation) {
    maybeCopies.emplace(m1, m2);
    translationRemove(maybeCopies->first, maybeCopies->second);
  }

  const auto &m1Ref = maybeCopies ? maybeCopies->first : m1;
  const auto &m2Ref = maybeCopies ? maybeCopies->second : m2;

  return compareFunc(m1Ref, m2Ref, scparams.distanceDifference);
}

} // namespace eonc::mat
