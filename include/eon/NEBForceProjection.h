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
#pragma once

#include "Eigen.h"
#include "Matter.h"

namespace eonc::neb {

/// Compute the tangent vector at image i using the improved tangent scheme.
/// Returns a normalized tangent vector.
AtomMatrix computeTangent(const AtomMatrix &posDiffNext,
                          const AtomMatrix &posDiffPrev, double energy,
                          double energyPrev, double energyNext,
                          bool use_old_tangent);

/// Compute the perpendicular component of force relative to the tangent.
AtomMatrix forcePerp(const AtomMatrix &force, const AtomMatrix &tangent);

/// Compute the climbing image projected force.
/// F_CI = F - 2*(F.t)*t + forceDNEB
AtomMatrix climbingImageForce(const AtomMatrix &force,
                              const AtomMatrix &tangent,
                              const AtomMatrix &forceDNEB);

/// Compute the doubly-nudged elastic band perpendicular spring force.
AtomMatrix computeDNEB(const AtomMatrix &forceSpring, const AtomMatrix &tangent,
                       const AtomMatrix &forcePerp);

/// Zero net translational force for fully free systems.
void zeroTranslation(AtomMatrix &projectedForce, int nFreeAtoms, int nAtoms);

} // namespace eonc::neb
