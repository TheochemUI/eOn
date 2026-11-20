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
#include "eon/NEBForceProjection.h"
#include "eon/HelperFunctions.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace eonc::neb {

// matDot is now in Eigen.h (shared across all files)

namespace {

template <typename Coordinate>
Coordinate tangentImpl(const Coordinate &posDiffNext,
                       const Coordinate &posDiffPrev, double energy,
                       double energyPrev, double energyNext,
                       bool useOldTangent) {
  Coordinate tangent;
  if (useOldTangent) {
    tangent = posDiffNext;
  } else if (energyNext > energy && energy > energyPrev) {
    tangent = posDiffNext;
  } else if (energy > energyNext && energyPrev > energy) {
    tangent = posDiffPrev;
  } else {
    const double energyDiffPrev = energyPrev - energy;
    const double energyDiffNext = energyNext - energy;
    const double minDiffEnergy =
        std::min(std::abs(energyDiffPrev), std::abs(energyDiffNext));
    const double maxDiffEnergy =
        std::max(std::abs(energyDiffPrev), std::abs(energyDiffNext));
    if (energyDiffPrev > energyDiffNext) {
      tangent = posDiffNext * minDiffEnergy + posDiffPrev * maxDiffEnergy;
    } else {
      tangent = posDiffNext * maxDiffEnergy + posDiffPrev * minDiffEnergy;
    }
  }
  return tangent;
}

} // namespace

AtomMatrix computeTangent(const AtomMatrix &posDiffNext,
                          const AtomMatrix &posDiffPrev, double energy,
                          double energyPrev, double energyNext,
                          bool use_old_tangent) {
  AtomMatrix tangent = tangentImpl(posDiffNext, posDiffPrev, energy, energyPrev,
                                   energyNext, use_old_tangent);
  double norm = tangent.norm();
  if (norm > 1e-10) {
    tangent /= norm;
  } else {
    tangent = posDiffNext;
    norm = tangent.norm();
    if (norm > 1e-10) {
      tangent /= norm;
    }
  }
  return tangent;
}

VectorXd computeTangent(const VectorXd &posDiffNext,
                        const VectorXd &posDiffPrev, double energy,
                        double energyPrev, double energyNext,
                        bool use_old_tangent) {
  VectorXd tangent = tangentImpl(posDiffNext, posDiffPrev, energy, energyPrev,
                                 energyNext, use_old_tangent);
  const double norm = tangent.norm();
  if (!(norm > 0.0)) {
    throw std::invalid_argument("NEB tangent has zero length");
  }
  return tangent / norm;
}

AtomMatrix forcePerp(const AtomMatrix &force, const AtomMatrix &tangent) {
  return force - matDot(force, tangent) * tangent;
}

VectorXd forcePerp(const VectorXd &force, const VectorXd &tangent) {
  return force - force.dot(tangent) * tangent;
}

AtomMatrix climbingImageForce(const AtomMatrix &force,
                              const AtomMatrix &tangent,
                              const AtomMatrix &forceDNEB) {
  return force - 2.0 * matDot(force, tangent) * tangent + forceDNEB;
}

AtomMatrix computeDNEB(const AtomMatrix &forceSpring, const AtomMatrix &tangent,
                       const AtomMatrix &fPerp, bool useSwitching) {
  AtomMatrix forceSpringPerp =
      forceSpring - matDot(forceSpring, tangent) * tangent;

  const double forceSpringPerpNorm = forceSpringPerp.norm();
  const double forcePerpNorm = fPerp.norm();

  if (forceSpringPerpNorm > 1e-10 && forcePerpNorm > 1e-10) {
    AtomMatrix forcePerpNormalized = fPerp / forcePerpNorm;
    AtomMatrix dneb =
        forceSpringPerp -
        matDot(forceSpringPerp, forcePerpNormalized) * forcePerpNormalized;

    if (useSwitching) {
      double switching = 2.0 / eonc::helpers::pi *
                         std::atan(forcePerpNorm * forcePerpNorm /
                                   (forceSpringPerpNorm * forceSpringPerpNorm));
      dneb *= switching;
    }
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
