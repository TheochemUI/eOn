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
#include "NEBTangent.h"

#include <algorithm>
#include <cmath>

namespace eonc::neb {

namespace {
/// Normalize tangent with fallback to posDiffNext if near-zero.
AtomMatrix normalizeTangent(AtomMatrix tang, const AtomMatrix &posDiffNext) {
  double norm = tang.norm();
  if (norm > 1e-10) {
    tang /= norm;
  } else {
    tang = posDiffNext;
    norm = tang.norm();
    if (norm > 1e-10) {
      tang /= norm;
    }
  }
  return tang;
}
} // namespace

AtomMatrix SimpleTangent::compute(const AtomMatrix &posDiffNext,
                                  const AtomMatrix & /*posDiffPrev*/,
                                  double /*energy*/, double /*energyPrev*/,
                                  double /*energyNext*/) const {
  return normalizeTangent(posDiffNext, posDiffNext);
}

AtomMatrix ImprovedTangent::compute(const AtomMatrix &posDiffNext,
                                    const AtomMatrix &posDiffPrev,
                                    double energy, double energyPrev,
                                    double energyNext) const {
  AtomMatrix tang;

  if (energyNext > energy && energy > energyPrev) {
    tang = posDiffNext;
  } else if (energy > energyNext && energyPrev > energy) {
    tang = posDiffPrev;
  } else {
    // Extremum: energy-weighted combination
    double energyDiffPrev = energyPrev - energy;
    double energyDiffNext = energyNext - energy;
    double minDiffEnergy =
        std::min(std::abs(energyDiffPrev), std::abs(energyDiffNext));
    double maxDiffEnergy =
        std::max(std::abs(energyDiffPrev), std::abs(energyDiffNext));

    if (energyDiffPrev > energyDiffNext) {
      tang = posDiffNext * minDiffEnergy + posDiffPrev * maxDiffEnergy;
    } else {
      tang = posDiffNext * maxDiffEnergy + posDiffPrev * minDiffEnergy;
    }
  }

  return normalizeTangent(tang, posDiffNext);
}

TangentStrategy buildTangentStrategy(const Parameters &params) {
  if (params.neb_options.climbing_image.use_old_tangent) {
    return SimpleTangent{};
  }
  return ImprovedTangent{};
}

} // namespace eonc::neb
