#include "MatterHelpers.hpp"
#include "HelperFunctions.h"

namespace eonc {

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

bool StructComparer::compare(const Matter &m1, const Matter &m2,
                             const bool indistinguishable) {

  if (m1.numberOfAtoms() != m2.numberOfAtoms())
    return false;
  if (scparams.checkRotation && indistinguishable) {
    return helper_functions::sortedR(m1, m2, scparams.distanceDifference);
  } else if (indistinguishable) {
    // XXX: This needs to work!
    // if (m1.numberOfFixedAtoms() == 0 and scparams.removeTranslation) {
    //   helper_functions::translationRemove(m1, m2);
    // }
    return helper_functions::identical(m1, m2, scparams.distanceDifference);
  } else if (scparams.checkRotation) {
    return helper_functions::rotationMatch(m1, m2, scparams.distanceDifference);
  } else {
    // XXX: This needs to work!
    // if (m1.numberOfFixedAtoms() == 0 and scparams.removeTranslation) {
    //   helper_functions::translationRemove(m1, m2);
    // }
    return (scparams.distanceDifference) > m1.perAtomNorm(m2);
  }
  return false;
}

} // namespace eonc
