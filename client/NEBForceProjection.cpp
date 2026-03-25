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
#include "NEBForceProjection.h"
#include "HelperFunctions.h"

#include <algorithm>
#include <cmath>

namespace eonc::neb {

// matDot is now in Eigen.h (shared across all files)

AtomMatrix computeTangent(const AtomMatrix &posDiffNext,
                          const AtomMatrix &posDiffPrev, double energy,
                          double energyPrev, double energyNext,
                          bool use_old_tangent) {
  AtomMatrix tang;

  if (use_old_tangent) {
    tang = posDiffNext;
  } else {
    // Improved tangent scheme
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
  }

  // Normalize with safety check
  double norm = tang.norm();
  if (norm > 1e-10) {
    tang /= norm;
  } else {
    // Fallback: use direction to next image
    tang = posDiffNext;
    norm = tang.norm();
    if (norm > 1e-10) {
      tang /= norm;
    }
  }

  return tang;
}

AtomMatrix forcePerp(const AtomMatrix &force, const AtomMatrix &tangent) {
  return force - matDot(force, tangent) * tangent;
}

AtomMatrix climbingImageForce(const AtomMatrix &force,
                              const AtomMatrix &tangent,
                              const AtomMatrix &forceDNEB) {
  return force - 2.0 * matDot(force, tangent) * tangent + forceDNEB;
}

AtomMatrix computeDNEB(const AtomMatrix &forceSpring, const AtomMatrix &tangent,
                       const AtomMatrix &fPerp) {
  AtomMatrix forceSpringPerp =
      forceSpring - matDot(forceSpring, tangent) * tangent;

  const double forceSpringPerpNorm = forceSpringPerp.norm();
  const double forcePerpNorm = fPerp.norm();

  if (forceSpringPerpNorm > 1e-10 && forcePerpNorm > 1e-10) {
    AtomMatrix forcePerpNormalized = fPerp / forcePerpNorm;
    AtomMatrix dneb =
        forceSpringPerp -
        matDot(forceSpringPerp, forcePerpNormalized) * forcePerpNormalized;

    double switching = 2.0 / eonc::helpers::pi *
                       std::atan(forcePerpNorm * forcePerpNorm /
                                 (forceSpringPerpNorm * forceSpringPerpNorm));
    dneb *= switching;
    return dneb;
  }

  return AtomMatrix::Zero(tangent.rows(), tangent.cols());
}

void zeroTranslation(AtomMatrix &projectedForce, int nFreeAtoms, int nAtoms) {
  if (nFreeAtoms == nAtoms) {
    for (int j = 0; j <= 2; j++) {
      double translationMag = projectedForce.col(j).sum();
      int natoms = projectedForce.col(j).size();
      projectedForce.col(j).array() -=
          translationMag / static_cast<double>(natoms);
    }
  }
}

} // namespace eonc::neb
