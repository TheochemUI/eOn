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
#include "client/matter/MatterHelpers.hpp"

namespace eonc::mat {

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
